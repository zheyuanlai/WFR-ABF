"""Batched PyTorch ABF(+Fisher--Rao) engine (GPU backend).

This is the fast backend for the 2D ``xi(x, y) = x`` study.  It simulates a
*batch* of independent runs at once: particle tensors have shape ``(B, N)`` with
``B`` = runs/seeds/configs in the batch and ``N`` = particles per run.  The
science is identical to the CPU reference (:mod:`abffr.simulation`):

* same potential / reaction coordinate (:func:`abffr.potentials.*_torch`),
* same biased overdamped dynamics,
* same ABF target ``F'(x) = E[dV/dx | X = x]``,
* same Fisher--Rao target densities and score,
* same metrics (each run produces a ``diag`` dict with the *identical* structure
  to :func:`abffr.simulation.run_simulation`, so :mod:`abffr.metrics` and
  :mod:`abffr.diagnostics` are reused verbatim).

Two estimator modes are provided:

``binned_smooth`` (production)
    Bin particles onto the x-grid (``scatter_add``) and Gaussian-smooth the
    count/force histograms -- an ``O(N + G)`` per-step estimator that is a
    Riemann-sum approximation of the CPU kernel estimator (see
    :mod:`abffr.torch_utils`).
``kernel_reference`` (validation only, slow)
    The exact ``O(G*N)`` Nadaraya--Watson kernel estimator and reflected KDE,
    a faithful GPU port of the CPU engine, used to bound the binning error.

A batch must share ``target_type``, ``eta``, ``fr_every`` and ``burnin_fraction``
(so the FR firing schedule and smoothing bandwidth are common); ``gamma`` and
``seed`` may vary per row.  Grouping is handled by :mod:`abffr.parallel`.
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import torch

from . import potentials, torch_utils as tu
from .io_utils import RunSpec, make_rng_streams
from .simulation import _init_positions  # reuse CPU init for matched seeds

EPS = 1e-12


@dataclass
class BatchResult:
    diags: List[Dict]          # one CPU-style diag dict per row
    runtime_seconds: float     # wall-clock for the whole batch
    device: str


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _init_batch_positions(specs, n_particles, domain, x_mode, y_mode,
                          device, dtype):
    """Per-row initial conditions, matched to the CPU init stream by seed.

    Using :func:`abffr.io_utils.make_rng_streams` (the same SeedSequence split
    as the CPU engine) means a given ``seed`` has *identical* initial conditions
    across methods and across the CPU/torch backends -- the basis for clean
    matched-seed comparisons.
    """
    xmin, xmax = domain["x_min"], domain["x_max"]
    ymin, ymax = domain["y_min"], domain["y_max"]
    X = np.empty((len(specs), n_particles), dtype=np.float64)
    Y = np.empty((len(specs), n_particles), dtype=np.float64)
    for b, spec in enumerate(specs):
        rng_init, _, _ = make_rng_streams(spec.seed)
        X[b] = _init_positions(rng_init, n_particles, xmin, xmax, x_mode)
        Y[b] = _init_positions(rng_init, n_particles, ymin, ymax, y_mode)
    return (torch.as_tensor(X, device=device, dtype=dtype),
            torch.as_tensor(Y, device=device, dtype=dtype))


def _empty_diag(target_type):
    return dict(
        target_type=target_type,
        steps=[], times=[],
        Fprime_hat=[], F_hat=[], p_hat_grid=[], q_target_grid=[],
        X_snap=[], Y_snap=[],
        barrier_crossings=[], n_unique_ancestors=[], gamma_eff=[],
        fr_applied=[], fr_event_fraction=[], fr_event_fraction_max=[],
        fr_events_total=[], score_mean=[], score_std=[], score_min=[],
        score_max=[],
        # GPU-only extra ancestor diagnostics (optional columns downstream).
        ancestor_ess=[], max_clone_multiplicity=[], target_l2=[],
    )


def _kernel_estimator(x_grid_t, X, Y, h):
    """Exact O(G*N) Nadaraya--Watson contribution for one step (per row).

    Returns ``(num_contrib, den_contrib)`` each ``(B, G)`` matching the CPU
    ``weights`` accumulation, for the ``kernel_reference`` validation mode.
    """
    # (B, G, N): kernel of every grid node against every particle.
    diff = x_grid_t.view(1, -1, 1) - X.unsqueeze(1)
    w = torch.exp(-0.5 * (diff / h) ** 2) / (h * np.sqrt(2.0 * np.pi))
    dvdx = potentials.dVdx_xy_torch(X, Y)
    return (w * dvdx.unsqueeze(1)).sum(-1), w.sum(-1)


def _kde_reflected(x_grid_t, X, eta, xmin, xmax):
    """Reflected-boundary KDE marginal on the grid (``kernel_reference`` mode).

    Mirror images about ``xmin``/``xmax`` remove KDE edge bias, matching
    :func:`abffr.simulation.kde_marginal`.  Returns ``(B, G)``.
    """
    N = X.shape[1]
    X_all = torch.cat([2.0 * xmin - X, X, 2.0 * xmax - X], dim=1)
    diff = x_grid_t.view(1, -1, 1) - X_all.unsqueeze(1)
    w = torch.exp(-0.5 * (diff / eta) ** 2) / (eta * np.sqrt(2.0 * np.pi))
    return w.sum(-1) / N


# --------------------------------------------------------------------------- #
# Fisher--Rao target densities and score (batched, on grid)
# --------------------------------------------------------------------------- #
def _build_target(target_type, Fhat_target, B_grid, F_ref_grid, p_hat,
                  beta, dx, width):
    """Batched FR target ``q`` (``(B, G)``, normalised on the grid)."""
    if target_type == "uniform":
        B, G = p_hat.shape
        return torch.full((B, G), 1.0 / max(width, EPS),
                          device=p_hat.device, dtype=p_hat.dtype)
    if target_type == "self":
        return tu.normalize_density(p_hat, dx)
    # estimated / oracle: q ~ exp(-beta (F_target - B)), normalised.
    F_target = F_ref_grid if target_type == "oracle" else Fhat_target
    exponent = -beta * (F_target - B_grid)
    exponent = exponent - exponent.max(dim=1, keepdim=True).values
    q = torch.exp(exponent)
    q = tu.normalize_density(q, dx)
    q = q.clamp_min(EPS)
    return tu.normalize_density(q, dx)


def _fr_score(X, p_hat, q_grid, x0, dx, beta, score_clip):
    """Batched Fisher--Rao score (mean-zero), matching ``simulation.fr_score``."""
    Zp = tu.trapezoid(p_hat, dx).clamp_min(EPS).unsqueeze(1)
    Zq = tu.trapezoid(q_grid, dx).clamp_min(EPS).unsqueeze(1)
    p_g = p_hat / Zp
    q_g = q_grid / Zq
    log_ratio_grid = torch.log(p_g.clamp_min(EPS)) - torch.log(q_g.clamp_min(EPS))
    baseline = tu.trapezoid(p_g * log_ratio_grid, dx).unsqueeze(1)  # KL(p||q)

    p_part = tu.interp1d(p_g, X, x0, dx).clamp_min(EPS)
    q_part = tu.interp1d(q_g, X, x0, dx).clamp_min(EPS)
    S = torch.log(p_part) - torch.log(q_part) - baseline
    if score_clip is not None:
        S = S.clamp(-float(score_clip), float(score_clip))
    return S


def _resample(X, Y, ancestors, S, gamma_eff, dt, max_event_fraction, gen,
              jitter, noise_scale):
    """Vectorised fixed-N Fisher--Rao birth--death.

    Per row: a particle with ``S_i > 0`` dies with prob ``1 - exp(-gamma S_i dt)``
    (same as the CPU engine); each dead slot is overwritten by a *clone* sampled
    from the under-represented pool ``{i : S_i < 0}`` with weights ``|S_i|``.
    This keeps ``N`` fixed, is fully GPU-vectorised, and matches the study spec's
    resampling design.  ``max_event_fraction`` caps the per-row replaced fraction
    by randomly sub-sampling deaths.  Returns ``(X, Y, ancestors, num_deaths)``.

    NB: this differs from the CPU ``resample_fixed_N`` (which allows
    ``n_die != n_clone`` with a random top-up); both are gentle, capped
    birth--death schemes and are validated to give same-order metrics.
    """
    B, N = X.shape
    g_dt = (gamma_eff.unsqueeze(1) * dt)
    pos = S > 0
    neg = S < 0
    p_death = torch.where(pos, (1.0 - torch.exp(-g_dt * S)).clamp(0.0, 1.0),
                          torch.zeros_like(S))
    u = torch.rand(S.shape, generator=gen, device=S.device, dtype=S.dtype)
    die = pos & (u < p_death)

    # Cap deaths per row at floor(max_event_fraction * N) via random sub-sample.
    if max_event_fraction is not None:
        cap = int(np.floor(float(max_event_fraction) * N))
        if cap <= 0:
            die = torch.zeros_like(die)
        else:
            key = torch.where(die, torch.rand(S.shape, generator=gen,
                                              device=S.device, dtype=S.dtype),
                              torch.full_like(S, float("inf")))
            sorted_key, _ = torch.sort(key, dim=1)
            thresh = sorted_key[:, min(cap - 1, N - 1)].unsqueeze(1)
            die = die & (key <= thresh)

    # Birth pool: under-represented particles weighted by |S| (= -S for neg).
    w_birth = torch.where(neg, -S, torch.zeros_like(S))
    has_birth = w_birth.sum(dim=1) > 0
    # Rows with no birth candidates cannot resample: cancel their deaths.
    die = die & has_birth.unsqueeze(1)
    # Safe weights for multinomial (rows without births fall back to uniform but
    # have no deaths to fill, so the draw is unused).
    w_safe = torch.where(has_birth.unsqueeze(1), w_birth,
                         torch.ones_like(w_birth))
    sources = torch.multinomial(w_safe, N, replacement=True, generator=gen)

    identity = torch.arange(N, device=X.device).expand(B, N)
    gather_idx = torch.where(die, sources, identity)

    X_new = torch.gather(X, 1, gather_idx)
    Y_new = torch.gather(Y, 1, gather_idx)
    anc_new = torch.gather(ancestors, 1, gather_idx)
    if jitter and jitter > 0.0:
        born = die
        X_new = X_new + born * jitter * noise_scale * torch.randn(
            X.shape, generator=gen, device=X.device, dtype=X.dtype)
        Y_new = Y_new + born * jitter * noise_scale * torch.randn(
            Y.shape, generator=gen, device=Y.device, dtype=Y.dtype)
    return X_new, Y_new, anc_new, die.sum(dim=1)


# --------------------------------------------------------------------------- #
# Core batched engine
# --------------------------------------------------------------------------- #
def run_batch(
    specs: List[RunSpec],
    *,
    cfg: Dict,
    x_grid: np.ndarray,
    F_ref: np.ndarray,
    Fprime_ref: np.ndarray,
    ev,
    device: torch.device,
    dtype: torch.dtype,
    estimator: str = "binned_smooth",
    base_seed: int = 0,
) -> BatchResult:
    """Simulate a batch of runs sharing (target_type, eta, fr_every, burnin)."""
    if not specs:
        return BatchResult([], 0.0, str(device))
    target_type = specs[0].target_type
    eta = float(specs[0].eta)
    fr_every = int(specs[0].fr_every)
    burnin_fraction = float(specs[0].burnin_fraction)
    for s in specs:  # invariant guarded by the batching layer
        if (s.target_type, float(s.eta), int(s.fr_every),
                float(s.burnin_fraction)) != (
                target_type, eta, fr_every, burnin_fraction):
            raise ValueError("run_batch requires a (target_type, eta, fr_every, "
                             "burnin_fraction)-homogeneous batch.")

    sim = cfg["simulation"]
    abf = cfg["abf"]
    fr = cfg.get("fr", {})
    beta = float(sim["beta"]); dt = float(sim["dt"])
    n_steps = int(sim["n_steps"]); n_particles = int(sim["n_particles"])
    eval_every = int(sim["eval_every"])
    domain = cfg["domain"]
    xmin, xmax = float(domain["x_min"]), float(domain["x_max"])
    ymin, ymax = float(domain["y_min"]), float(domain["y_max"])
    width = xmax - xmin

    h = float(abf["h"])
    update_every = max(1, int(abf.get("update_every", 1)))
    min_count = float(abf.get("min_count", 1.0))
    # Estimated-target free-energy EMA: the study spec names this
    # ``fr.target_ema_alpha``; fall back to the CPU engine's ``abf.ema_alpha``
    # (both default to 0.05, so the CPU/torch estimated target stays consistent).
    ema_alpha = float(fr.get("target_ema_alpha", abf.get("ema_alpha", 0.05)))
    # EMA cadence correction: applying the per-step alpha only every
    # ``update_every`` steps would slow the EMA; use the matched decay.
    ema_alpha_eff = 1.0 - (1.0 - ema_alpha) ** update_every

    ramp_fraction = float(fr.get("ramp_fraction", 0.1))
    score_clip = fr.get("score_clip", 5.0)
    max_event_fraction = fr.get("max_event_fraction", 0.10)
    jitter = float(fr.get("jitter", 0.0))
    x_init_mode = sim.get("x_init_mode", "mixed")
    y_init_mode = sim.get("y_init_mode", "mixed")
    x_barrier = float(getattr(ev, "x_barrier", 0.0))

    B = len(specs)
    G = len(x_grid)
    x_grid_t = torch.as_tensor(np.asarray(x_grid), device=device, dtype=dtype)
    x0 = float(x_grid[0]); dx = tu.grid_spacing(x_grid_t)
    idx0 = int(np.argmin(np.abs(np.asarray(x_grid))))
    F_ref_t = torch.as_tensor(np.asarray(F_ref), device=device,
                              dtype=dtype).view(1, G)

    fr_enabled = (target_type != "none")
    fr_burnin = int(round(burnin_fraction * n_steps))
    ramp_steps = int(round(ramp_fraction * n_steps))
    gamma_vec = torch.as_tensor([float(s.gamma) for s in specs], device=device,
                                dtype=dtype)

    # Smoothing kernels for the binned-smooth estimator.
    k_h, r_h = tu.gaussian_kernel1d(h, dx, device, dtype)
    k_eta, r_eta = tu.gaussian_kernel1d(eta, dx, device, dtype)
    use_kernel_ref = (estimator == "kernel_reference")

    # RNG: per-row matched initial conditions; one reproducible noise stream for
    # the whole batch (seeded from the batch's run_ids, so results are
    # deterministic given the config). See module/README notes on RNG.
    gen = tu.make_generator(
        tu.stable_seed(base_seed, *[s.run_id for s in specs]), device)
    X, Y = _init_batch_positions(specs, n_particles, domain, x_init_mode,
                                 y_init_mode, device, dtype)
    ancestors = torch.arange(n_particles, device=device).expand(B, n_particles).contiguous()
    noise_scale = float(np.sqrt(2.0 * dt / beta))

    # ABF accumulators and current grid estimates.
    C_acc = torch.zeros((B, G), device=device, dtype=dtype)
    S_acc = torch.zeros((B, G), device=device, dtype=dtype)
    Fprime_hat = torch.zeros((B, G), device=device, dtype=dtype)
    F_hat = torch.zeros((B, G), device=device, dtype=dtype)
    Fhat_target = torch.zeros((B, G), device=device, dtype=dtype)

    barrier_crossings = torch.zeros(B, device=device, dtype=torch.long)
    prev_sign = torch.sign(X - x_barrier)

    # Windowed FR accumulators (reset every snapshot), all (B,) device tensors.
    win_n = torch.zeros(B, device=device, dtype=torch.long)
    win_events = torch.zeros(B, device=device, dtype=torch.long)
    win_frac_sum = torch.zeros(B, device=device, dtype=dtype)
    win_frac_max = torch.zeros(B, device=device, dtype=dtype)
    win_smean_sum = torch.zeros(B, device=device, dtype=dtype)
    win_sstd_sum = torch.zeros(B, device=device, dtype=dtype)
    win_smin = torch.full((B,), float("inf"), device=device, dtype=dtype)
    win_smax = torch.full((B,), float("-inf"), device=device, dtype=dtype)
    win_target_l2 = torch.zeros(B, device=device, dtype=dtype)

    diags = [_empty_diag(target_type) for _ in range(B)]

    def recompute_grid():
        nonlocal Fprime_hat, F_hat
        if use_kernel_ref:
            # S_acc / C_acc hold the exact kernel-accumulated numerator
            # (force-weighted) and denominator (weights); no smoothing needed.
            Fprime_hat = S_acc / (C_acc + min_count + EPS)
        else:
            num_s = tu.smooth_grid(S_acc, k_h, r_h, dx)
            den_s = tu.smooth_grid(C_acc, k_h, r_h, dx)
            Fprime_hat = num_s / (den_s + min_count + EPS)
        F_hat = tu.center_at_index(tu.cumulative_trapezoid(Fprime_hat, dx), idx0)

    def current_p_hat(Xc):
        if use_kernel_ref:
            p = _kde_reflected(x_grid_t, Xc, eta, xmin, xmax)
        else:
            hist = tu.scatter_grid(tu.nearest_index(Xc, x0, dx, G), G)
            p = tu.smooth_grid(hist, k_eta, r_eta, dx) / n_particles
        return tu.normalize_density(p, dx)

    t0 = time.time()
    for step in range(n_steps):
        dvdx = potentials.dVdx_xy_torch(X, Y)
        dvdy = potentials.dVdy_xy_torch(X, Y)

        # --- ABF accumulation (every step, from current X, Y) --------------- #
        if use_kernel_ref:
            num_c, den_c = _kernel_estimator(x_grid_t, X, Y, h)
            S_acc += num_c
            C_acc += den_c
        else:
            idx = tu.nearest_index(X, x0, dx, G)
            C_acc += tu.scatter_grid(idx, G)
            S_acc += tu.scatter_grid(idx, G, dvdx)

        if step % update_every == 0:
            recompute_grid()
            Fhat_target = (1.0 - ema_alpha_eff) * Fhat_target + ema_alpha_eff * F_hat

        # --- Langevin + ABF proposal --------------------------------------- #
        abf_at_X = tu.interp1d(Fprime_hat, X, x0, dx)
        noise_x = torch.randn(X.shape, generator=gen, device=device, dtype=dtype)
        noise_y = torch.randn(Y.shape, generator=gen, device=device, dtype=dtype)
        X_prop = tu.reflect_into(X + (-dvdx + abf_at_X) * dt + noise_scale * noise_x,
                                 xmin, xmax)
        Y_prop = tu.reflect_into(Y + (-dvdy) * dt + noise_scale * noise_y,
                                 ymin, ymax)
        tu.assert_finite("X_prop", X_prop)
        tu.assert_finite("Y_prop", Y_prop)

        # --- Barrier crossings on the Langevin move ------------------------ #
        new_sign = torch.sign(X_prop - x_barrier)
        crossed = (new_sign != prev_sign) & (new_sign != 0) & (prev_sign != 0)
        barrier_crossings += crossed.sum(dim=1)

        # --- Fisher--Rao birth--death -------------------------------------- #
        gamma_pos = (step > fr_burnin) if ramp_steps > 0 else (step >= fr_burnin)
        do_fr = (fr_enabled and step >= fr_burnin
                 and ((step - fr_burnin) % fr_every == 0) and gamma_pos)
        if fr_enabled:
            s_ramp = max((step - fr_burnin) / ramp_steps, 0.0) if ramp_steps > 0 else None
            if step < fr_burnin:
                gamma_eff = torch.zeros_like(gamma_vec)
            elif ramp_steps <= 0:
                gamma_eff = gamma_vec
            else:
                gamma_eff = gamma_vec * (1.0 - np.exp(-s_ramp))
        else:
            gamma_eff = torch.zeros_like(gamma_vec)

        if do_fr:
            p_hat = current_p_hat(X_prop)
            q_grid = _build_target(target_type, Fhat_target, F_hat, F_ref_t,
                                   p_hat, beta, dx, width)
            S = _fr_score(X_prop, p_hat, q_grid, x0, dx, beta, score_clip)
            X, Y, ancestors, ndeath = _resample(
                X_prop, Y_prop, ancestors, S, gamma_eff, dt,
                max_event_fraction, gen, jitter, noise_scale)
            # Windowed diagnostics.
            frac = ndeath.to(dtype) / n_particles
            win_n += 1
            win_events += ndeath
            win_frac_sum += frac
            win_frac_max = torch.maximum(win_frac_max, frac)
            win_smean_sum += S.mean(dim=1)
            win_sstd_sum += S.std(dim=1)
            win_smin = torch.minimum(win_smin, S.amin(dim=1))
            win_smax = torch.maximum(win_smax, S.amax(dim=1))
            win_target_l2 += _marginal_l2(p_hat, q_grid, dx)
        else:
            X, Y = X_prop, Y_prop

        prev_sign = torch.sign(X - x_barrier)

        # --- Snapshot ------------------------------------------------------- #
        if step % eval_every == 0 or step == n_steps - 1:
            recompute_grid()
            p_snap = current_p_hat(X)
            q_snap = _build_target(target_type, Fhat_target, F_hat, F_ref_t,
                                   p_snap, beta, dx, width)
            _record_snapshot(
                diags, step, step * dt, Fprime_hat, F_hat, p_snap, q_snap,
                X, Y, ancestors, barrier_crossings, gamma_eff,
                win_n, win_events, win_frac_sum, win_frac_max, win_smean_sum,
                win_sstd_sum, win_smin, win_smax, win_target_l2, n_particles)
            # reset window
            win_n = torch.zeros(B, device=device, dtype=torch.long)
            win_events = torch.zeros(B, device=device, dtype=torch.long)
            win_frac_sum = torch.zeros(B, device=device, dtype=dtype)
            win_frac_max = torch.zeros(B, device=device, dtype=dtype)
            win_smean_sum = torch.zeros(B, device=device, dtype=dtype)
            win_sstd_sum = torch.zeros(B, device=device, dtype=dtype)
            win_smin = torch.full((B,), float("inf"), device=device, dtype=dtype)
            win_smax = torch.full((B,), float("-inf"), device=device, dtype=dtype)
            win_target_l2 = torch.zeros(B, device=device, dtype=dtype)

    if device.type == "cuda":
        torch.cuda.synchronize()
    runtime = time.time() - t0

    for d in diags:
        d["fr_burnin"] = fr_burnin
        d["fr_every"] = int(fr_every)
        d["n_steps"] = int(n_steps)
        d["dt"] = float(dt)
    return BatchResult(diags, runtime, str(device))


def _marginal_l2(p_hat, q_grid, dx):
    """L2 distance between (renormalised) ``p_hat`` and ``q`` per row -> (B,)."""
    p = p_hat / tu.trapezoid(p_hat, dx).clamp_min(EPS).unsqueeze(1)
    q = q_grid / tu.trapezoid(q_grid, dx).clamp_min(EPS).unsqueeze(1)
    width = (p.shape[1] - 1) * dx
    return torch.sqrt((tu.trapezoid((p - q) ** 2, dx) / max(width, EPS)).clamp_min(0.0))


def _record_snapshot(diags, step, t, Fprime_hat, F_hat, p_snap, q_snap, X, Y,
                     ancestors, barrier_crossings, gamma_eff, win_n, win_events,
                     win_frac_sum, win_frac_max, win_smean_sum, win_sstd_sum,
                     win_smin, win_smax, win_target_l2, n_particles):
    """Move per-row snapshot data to numpy and append to each run's diag.

    This is the *only* GPU->CPU transfer inside the time loop and it happens at
    ``eval_every`` cadence on grid-sized arrays plus the particle snapshot.
    """
    Fp = Fprime_hat.detach().cpu().numpy()
    F = F_hat.detach().cpu().numpy()
    P = p_snap.detach().cpu().numpy()
    Q = q_snap.detach().cpu().numpy()
    Xn = X.detach().cpu().numpy()
    Yn = Y.detach().cpu().numpy()
    anc = ancestors.detach().cpu().numpy()
    bc = barrier_crossings.detach().cpu().numpy()
    ge = gamma_eff.detach().cpu().numpy()
    wn = win_n.detach().cpu().numpy()
    we = win_events.detach().cpu().numpy()
    wfs = win_frac_sum.detach().cpu().numpy()
    wfm = win_frac_max.detach().cpu().numpy()
    wsm = win_smean_sum.detach().cpu().numpy()
    wss = win_sstd_sum.detach().cpu().numpy()
    wmin = win_smin.detach().cpu().numpy()
    wmax = win_smax.detach().cpu().numpy()
    wtl = win_target_l2.detach().cpu().numpy()

    for b, d in enumerate(diags):
        nfr = int(wn[b])
        d["steps"].append(int(step)); d["times"].append(float(t))
        d["Fprime_hat"].append(Fp[b].copy()); d["F_hat"].append(F[b].copy())
        d["p_hat_grid"].append(P[b].copy()); d["q_target_grid"].append(Q[b].copy())
        d["X_snap"].append(Xn[b].copy()); d["Y_snap"].append(Yn[b].copy())
        d["barrier_crossings"].append(int(bc[b]))
        anc_b = anc[b]
        counts = np.bincount(anc_b, minlength=n_particles).astype(np.float64)
        nz = counts[counts > 0]
        d["n_unique_ancestors"].append(int((counts > 0).sum()))
        d["ancestor_ess"].append(float(nz.sum() ** 2 / np.maximum((nz ** 2).sum(), 1.0)))
        d["max_clone_multiplicity"].append(int(counts.max()))
        d["gamma_eff"].append(float(ge[b]))
        d["fr_applied"].append(bool(nfr > 0))
        d["fr_event_fraction"].append(float(wfs[b] / nfr) if nfr else 0.0)
        d["fr_event_fraction_max"].append(float(wfm[b]) if nfr else 0.0)
        d["fr_events_total"].append(int(we[b]))
        d["score_mean"].append(float(wsm[b] / nfr) if nfr else float("nan"))
        d["score_std"].append(float(wss[b] / nfr) if nfr else float("nan"))
        d["score_min"].append(float(wmin[b]) if nfr else float("nan"))
        d["score_max"].append(float(wmax[b]) if nfr else float("nan"))
        d["target_l2"].append(float(wtl[b] / nfr) if nfr else float("nan"))
