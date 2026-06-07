"""Production core for the WCA-dimer ABF vs ABF+Fisher-Rao study.

Refactored and extended from ``scratch/abffr_core.py`` (the validated demo core).
Adds the diagnostic FR target variants (``fr_uniform``, ``fr_oracle``), Stage-C
mechanism diagnostics (RC marginal, FR target, per-bin effective counts, region
fractions, birth/death locations, ancestor ESS) and a no-oracle-leakage audit.

Method names (``method`` argument to :func:`run_sampler_gpu`):
    abf           ABF only baseline.
    fr_estimated  ABF + FR with an ONLINE EMA estimated target. The practical method.
    fr_uniform    ABF + FR with a uniform-in-z target. Naive ablation (diagnostic).
    fr_oracle     ABF + FR with the TI reference target. Positive control (diagnostic).

CRITICAL ALGORITHMIC RULE
-------------------------
Only ``fr_oracle`` may receive ``oracle_free_energy``. ``fr_estimated`` (and ``abf``
and ``fr_uniform``) MUST NOT receive it; the sampler asserts this. The TI reference
is used ONLY for post-hoc L2 evaluation, and for the ``fr_oracle`` diagnostic.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time
from dataclasses import dataclass, replace, asdict

import numpy as np
import torch
import torch.nn.functional as F

EPS = 1.0e-12

FR_METHODS = ("fr_estimated", "fr_uniform", "fr_oracle")
ALL_METHODS = ("abf",) + FR_METHODS


def choose_device():
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


DEVICE = choose_device()
DTYPE = torch.float32
if DEVICE.type == "cuda":
    torch.backends.cuda.matmul.allow_tf32 = True
    torch.backends.cudnn.allow_tf32 = True


@dataclass(frozen=True)
class DimerWCAParams:
    n_dim: int = 10
    a: float = 1.5
    sigma: float = 1.0
    epsilon: float = 1.0
    h: float = 2.0
    w: float = 2.0
    beta: float = 1.0
    min_r: float = 0.65
    force_clip: float = 250.0

    @property
    def n_particles(self):
        return self.n_dim * self.n_dim

    @property
    def box_length(self):
        return self.n_dim * self.a

    @property
    def r0(self):
        return 2.0 ** (1.0 / 6.0) * self.sigma

    @property
    def wca_cutoff(self):
        return self.r0


@dataclass(frozen=True)
class SimConfig:
    n_replicas: int = 1024
    dt: float = 2.0e-3
    n_steps: int = 250_000
    save_every: int = 2500
    seed: int = 42

    z_min: float = -0.2
    z_max: float = 1.2
    n_grid: int = 160

    abf_bandwidth: float = 0.025
    kde_bandwidth: float = 0.070
    abf_smooth_sigma: float = 0.50

    mean_force_sample_clip: float = 500.0
    use_clipped_force_for_mean_force: bool = True
    abf_bias_scale: float = 1.0
    abf_edge_extrapolate: bool = True
    boundary_wall_strength: float = 80.0
    abf_force_clip: float = 40.0
    abf_warmup_steps: int = 10_000
    estimator_burn_in_steps: int = 10_000

    # Fisher-Rao (online estimated target + diagnostic variants).
    fr_rate: float = 0.10
    score_clip: float = 2.0
    fr_start_steps: int = 20_000
    fr_every: int = 5
    target_ema_rate: float = 0.005
    max_event_fraction: float = 0.02

    # Region boundaries (compact / transition / stretched) for diagnostics.
    transition_lo: float = 0.25
    transition_hi: float = 0.75
    # Interior evaluation window (exclude wall-affected edges) as a fraction.
    eval_z_lo: float = 0.0
    eval_z_hi: float = 1.0

    def config_hash(self) -> str:
        return hashlib.md5(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()[:12]


@dataclass(frozen=True)
class TIConfig:
    z_min: float = -0.2
    z_max: float = 1.2
    n_z: int = 51
    n_replicas: int = 128
    dt: float = 2.0e-3
    n_thermalization: int = 5_000
    n_steps: int = 50_000
    sample_every: int = 20
    seed: int = 314159
    smooth_sigma: float = 1.0
    z_chunk_size: int = 8
# === END PART 1 ===


def as_tensor(x, device=DEVICE, dtype=DTYPE):
    return torch.as_tensor(x, device=device, dtype=dtype)


def wrap_positions(q, L):
    return torch.remainder(q, L)


def minimum_image(delta, L):
    return delta - L * torch.round(delta / L)


def clip_forces(forces, force_clip):
    norm = torch.linalg.norm(forces, dim=-1, keepdim=True)
    scale = torch.clamp(force_clip / torch.clamp(norm, min=EPS), max=1.0)
    return forces * scale


class WCADimerEngine:
    """GPU pair-list force engine for the dimer in WCA solvent."""

    def __init__(self, params, device=DEVICE, dtype=DTYPE):
        self.params = params
        self.device = device
        self.dtype = dtype
        self.N = params.n_particles
        self.L = float(params.box_length)
        self.sigma = float(params.sigma)
        self.epsilon = float(params.epsilon)
        self.h = float(params.h)
        self.w = float(params.w)
        self.r0 = float(params.r0)
        self.cutoff = float(params.wca_cutoff)
        self.min_r = float(params.min_r)

        pair_i, pair_j = torch.triu_indices(self.N, self.N, offset=1, device=device)
        keep = ~((pair_i == 0) & (pair_j == 1))
        self.pair_i = pair_i[keep].long()
        self.pair_j = pair_j[keep].long()
        self.n_pairs = int(self.pair_i.numel())

    def force(self, q, compute_energy=False):
        B = q.shape[0]
        qi = q.index_select(1, self.pair_i)
        qj = q.index_select(1, self.pair_j)
        delta = minimum_image(qi - qj, self.L)
        r = torch.linalg.norm(delta, dim=-1)
        r_safe = torch.clamp(r, min=self.min_r * self.sigma)

        active = r <= self.cutoff
        inv = self.sigma / r_safe
        inv6 = inv ** 6
        inv12 = inv6 ** 2
        dVdr = 4.0 * self.epsilon * (-12.0 * inv12 / r_safe + 6.0 * inv6 / r_safe)
        dVdr = torch.where(active, dVdr, torch.zeros_like(dVdr))
        f_pair = (-dVdr / torch.clamp(r_safe, min=EPS)).unsqueeze(-1) * delta

        forces = torch.zeros_like(q)
        idx_i = self.pair_i.view(1, -1, 1).expand(B, -1, 2)
        idx_j = self.pair_j.view(1, -1, 1).expand(B, -1, 2)
        forces.scatter_add_(1, idx_i, f_pair)
        forces.scatter_add_(1, idx_j, -f_pair)

        d01 = minimum_image(q[:, 0, :] - q[:, 1, :], self.L)
        r01 = torch.linalg.norm(d01, dim=1).clamp_min(EPS)
        u = (r01 - self.r0 - self.w) / self.w
        dVdr_dim = -4.0 * self.h * u * (1.0 - u ** 2) / self.w
        f01 = (-dVdr_dim / r01).unsqueeze(-1) * d01
        forces[:, 0, :] += f01
        forces[:, 1, :] -= f01

        if not compute_energy:
            return forces

        V_wca = 4.0 * self.epsilon * (inv12 - inv6) + self.epsilon
        V_wca = torch.where(active, V_wca, torch.zeros_like(V_wca)).sum(dim=1)
        V_dim = self.h * (1.0 - u ** 2) ** 2
        return V_wca + V_dim, forces


# ---------------------------------------------------------------------------
# Geometry, local mean force, initial conditions
# ---------------------------------------------------------------------------
def dimer_displacement_and_length(q, params):
    d01 = minimum_image(q[:, 0, :] - q[:, 1, :], params.box_length)
    r01 = torch.linalg.norm(d01, dim=1).clamp_min(EPS)
    return d01, r01


def reaction_coordinate(q, params):
    _, r01 = dimer_displacement_and_length(q, params)
    return (r01 - params.r0) / (2.0 * params.w)


def grad_xi_dimer(q, params):
    d01, r01 = dimer_displacement_and_length(q, params)
    grad0 = d01 / (2.0 * params.w * r01[:, None])
    return grad0, -grad0


def local_mean_force(q, physical_forces, params):
    d01, r01 = dimer_displacement_and_length(q, params)
    grad_difference = physical_forces[:, 1, :] - physical_forces[:, 0, :]
    energetic = (params.w / r01) * torch.sum(d01 * grad_difference, dim=1)
    entropic = (2.0 * params.w) / (params.beta * r01)
    return energetic - entropic


def add_abf_force(q, forces, mean_force_at_z, params):
    grad0, grad1 = grad_xi_dimer(q, params)
    out = forces.clone()
    out[:, 0, :] += mean_force_at_z[:, None] * grad0
    out[:, 1, :] += mean_force_at_z[:, None] * grad1
    return out


def add_reaction_coordinate_wall_force(q, forces, z, sim, params):
    if sim.boundary_wall_strength <= 0.0:
        return forces
    upper_excess = torch.clamp(z - sim.z_max, min=0.0)
    lower_excess = torch.clamp(z - sim.z_min, max=0.0)
    dU_wall_dz = sim.boundary_wall_strength * (upper_excess + lower_excess)
    return add_abf_force(q, forces, -dU_wall_dz, params)


def lattice_initial_conditions(params, n_replicas, device=DEVICE, dtype=DTYPE, seed=0, jitter=0.015):
    g = torch.Generator(device=device)
    g.manual_seed(seed)
    coords = [((0.5 + i) * params.a, (0.5 + j) * params.a)
              for i in range(params.n_dim) for j in range(params.n_dim)]
    base = torch.tensor(coords, device=device, dtype=dtype)
    q = base.unsqueeze(0).repeat(n_replicas, 1, 1)
    shifts = torch.rand((n_replicas, 1, 2), device=device, dtype=dtype, generator=g) * params.box_length
    q = wrap_positions(q + shifts, params.box_length)
    if jitter > 0.0:
        q[:, 2:, :] = wrap_positions(
            q[:, 2:, :] + jitter * torch.randn(q[:, 2:, :].shape, device=device, dtype=dtype, generator=g),
            params.box_length,
        )
    q[:, 1, :] = wrap_positions(q[:, 0, :] + torch.tensor([0.0, params.r0], device=device, dtype=dtype), params.box_length)
    return q


def project_dimer_to_z(q, z, params):
    z = torch.as_tensor(z, device=q.device, dtype=q.dtype)
    if z.ndim == 0:
        z = z.expand(q.shape[0])
    target_r = params.r0 + 2.0 * params.w * z
    d01, r01 = dimer_displacement_and_length(q, params)
    direction = d01 / r01[:, None]
    midpoint = q[:, 1, :] + 0.5 * d01
    target_d = target_r[:, None] * direction
    out = q.clone()
    out[:, 0, :] = midpoint + 0.5 * target_d
    out[:, 1, :] = midpoint - 0.5 * target_d
    return wrap_positions(out, params.box_length)


# ---------------------------------------------------------------------------
# KDE / ABF helpers
# ---------------------------------------------------------------------------
def gaussian_kernel_torch(diff, bandwidth):
    bw = float(max(bandwidth, EPS))
    return torch.exp(-0.5 * (diff / bw) ** 2) / (bw * math.sqrt(2.0 * math.pi))


def smooth_profile_torch(y, sigma_grid_points):
    sigma = float(sigma_grid_points)
    if sigma <= 0:
        return y.clone()
    radius = max(1, int(math.ceil(4.0 * sigma)))
    x = torch.arange(-radius, radius + 1, device=y.device, dtype=y.dtype)
    kernel = torch.exp(-0.5 * (x / sigma) ** 2)
    kernel = kernel / kernel.sum()
    yp = F.pad(y.view(1, 1, -1), (radius, radius), mode="replicate")
    return F.conv1d(yp, kernel.view(1, 1, -1)).view(-1)


def cumulative_trapezoid_torch(y, x):
    out = torch.zeros_like(y)
    out[1:] = torch.cumsum(0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1]), dim=0)
    return out


def normalize_profile_zero_at_midpoint_torch(profile, grid, midpoint=0.5):
    idx = torch.argmin(torch.abs(grid - midpoint))
    return profile - profile[idx]


def trapz_torch(y, x):
    return torch.sum(0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1]))


def normalize_density_on_grid_torch(p_grid, grid, eps=EPS):
    p_grid = torch.clamp(p_grid, min=eps)
    mass = trapz_torch(p_grid, grid)
    return p_grid / torch.clamp(mass, min=eps)


def kde_1d_torch(eval_points, samples, bandwidth, z_min, z_max):
    reflected = torch.cat([samples, 2.0 * z_min - samples, 2.0 * z_max - samples])
    diff = eval_points[:, None] - reflected[None, :]
    p = gaussian_kernel_torch(diff, bandwidth).sum(dim=1) / max(samples.numel(), 1)
    return torch.clamp(p, min=EPS)


def free_energy_from_density_torch(p_grid, grid, beta):
    Fz = -(1.0 / beta) * torch.log(torch.clamp(p_grid, min=EPS))
    return normalize_profile_zero_at_midpoint_torch(Fz, grid)


def mean_force_from_free_energy_torch(free_energy, grid):
    out = torch.empty_like(free_energy)
    out[1:-1] = (free_energy[2:] - free_energy[:-2]) / (grid[2:] - grid[:-2])
    out[0] = (free_energy[1] - free_energy[0]) / (grid[1] - grid[0])
    out[-1] = (free_energy[-1] - free_energy[-2]) / (grid[-1] - grid[-2])
    return out


def interp_uniform_grid(profile, grid, z, outside_value=0.0):
    dz = grid[1] - grid[0]
    x = (z - grid[0]) / dz
    idx0 = torch.floor(x).long()
    inside = (idx0 >= 0) & (idx0 < grid.numel() - 1)
    idx0c = idx0.clamp(0, grid.numel() - 2)
    frac = (x - idx0c.to(z.dtype)).clamp(0.0, 1.0)
    val = (1.0 - frac) * profile[idx0c] + frac * profile[idx0c + 1]
    return torch.where(inside, val, torch.full_like(val, outside_value))


def interp_uniform_grid_edge(profile, grid, z):
    dz = grid[1] - grid[0]
    x = (z - grid[0]) / dz
    idx0 = torch.floor(x).long()
    idx0c = idx0.clamp(0, grid.numel() - 2)
    frac = (x - idx0c.to(z.dtype)).clamp(0.0, 1.0)
    val = (1.0 - frac) * profile[idx0c] + frac * profile[idx0c + 1]
    val = torch.where(z < grid[0], profile[0].expand_as(val), val)
    val = torch.where(z > grid[-1], profile[-1].expand_as(val), val)
    return val


class TorchKernelABFEstimator:
    def __init__(self, z_grid, bandwidth, smooth_sigma=0.0, edge_extrapolate=False):
        self.z_grid = z_grid
        self.bandwidth = float(bandwidth)
        self.smooth_sigma = float(smooth_sigma)
        self.edge_extrapolate = bool(edge_extrapolate)
        self.num = torch.zeros_like(z_grid)
        self.den = torch.zeros_like(z_grid)
        self.n_updates = 0

    def update(self, z_samples, force_samples):
        weights = gaussian_kernel_torch(self.z_grid[:, None] - z_samples[None, :], self.bandwidth)
        self.num += torch.sum(weights * force_samples[None, :], dim=1)
        self.den += torch.sum(weights, dim=1)
        self.n_updates += int(z_samples.numel())

    def mean_force_profile(self):
        raw = self.num / torch.clamp(self.den, min=EPS)
        raw = torch.where(self.den > EPS, raw, torch.zeros_like(raw))
        return smooth_profile_torch(raw, self.smooth_sigma)

    def evaluate(self, z_samples):
        profile = self.mean_force_profile()
        if self.edge_extrapolate:
            return interp_uniform_grid_edge(profile, self.z_grid, z_samples)
        return interp_uniform_grid(profile, self.z_grid, z_samples, outside_value=0.0)

    def pmf_profile(self):
        pmf = cumulative_trapezoid_torch(self.mean_force_profile(), self.z_grid)
        return normalize_profile_zero_at_midpoint_torch(pmf, self.z_grid)

    def effective_counts(self):
        """Kernel-accumulated weight per grid bin (proxy for N_eff(z_j))."""
        return self.den.clone()


# ---------------------------------------------------------------------------
# Fisher-Rao target densities (estimated / uniform / oracle) and score
# ---------------------------------------------------------------------------
def recentered_clipped_score_torch(raw_score, score_clip):
    """Mean-zero, clipped score. Idempotent: applying twice == applying once.

    Ends on a recenter (not a bare clamp) so a second application is a guaranteed
    no-op; the score is normalized once in fr_score_torch and the birth-death step
    is defensive without shifting the death/birth partition near S_i=0.
    """
    score = raw_score - raw_score.mean()
    for _ in range(3):
        score = torch.clamp(score, -score_clip, score_clip)
        score = score - score.mean()
    return score


def fr_target_estimated_torch(grid, F_target_ema, current_bias_profile, beta, eps=EPS):
    """ESTIMATED-target marginal q_n(z) ~ exp[-beta (F_target_ema(z) - B_n(z))].

    F_target_ema is the ONLINE EMA of the unscaled ABF free energy A_hat_n, and
    B_n = current_bias_profile = abf_scale * A_hat_n. NO TI / oracle input.
    """
    if F_target_ema is None:
        raise ValueError("fr_estimated requires F_target_ema (init at fr_start_steps).")
    if current_bias_profile is None:
        raise ValueError("fr_estimated requires the current ABF bias B_n.")
    log_q = -beta * (F_target_ema - current_bias_profile)
    log_q = log_q - torch.max(log_q)
    return normalize_density_on_grid_torch(torch.exp(log_q), grid, eps=eps)


def fr_target_uniform_torch(grid, eps=EPS):
    """Naive ablation: uniform target on the z-grid."""
    return normalize_density_on_grid_torch(torch.ones_like(grid), grid, eps=eps)


def fr_target_oracle_torch(grid, oracle_free_energy, current_bias_profile, beta, eps=EPS):
    """DIAGNOSTIC oracle target q_n(z) ~ exp[-beta (F_ref(z) - B_n(z))].

    Uses the (unknown in practice) TI reference free energy. Positive control only.
    """
    if oracle_free_energy is None:
        raise ValueError("fr_oracle requires oracle_free_energy (the TI reference).")
    if current_bias_profile is None:
        raise ValueError("fr_oracle requires the current ABF bias B_n.")
    log_q = -beta * (oracle_free_energy - current_bias_profile)
    log_q = log_q - torch.max(log_q)
    return normalize_density_on_grid_torch(torch.exp(log_q), grid, eps=eps)


def fr_score_torch(z_samples, grid, sim, q_grid, eps=EPS):
    """Log-density-ratio FR score against a prebuilt target q_n.

    S_i = log p_hat(Z_i) - log q_n(Z_i) - KL(p_hat || q_n), then clip + recenter.
    Returns (score, p_grid, q_grid, kl_pq).
    """
    p_grid = normalize_density_on_grid_torch(
        kde_1d_torch(grid, z_samples, sim.kde_bandwidth, sim.z_min, sim.z_max), grid, eps=eps
    )
    p_at = interp_uniform_grid_edge(p_grid, grid, z_samples)
    q_at = interp_uniform_grid_edge(q_grid, grid, z_samples)
    log_ratio_grid = torch.log(torch.clamp(p_grid, min=eps)) - torch.log(torch.clamp(q_grid, min=eps))
    kl_pq = trapz_torch(p_grid * log_ratio_grid, grid)
    raw_score = (
        torch.log(torch.clamp(p_at, min=eps))
        - torch.log(torch.clamp(q_at, min=eps))
        - kl_pq
    )
    score = recentered_clipped_score_torch(raw_score, sim.score_clip)
    return score, p_grid, q_grid, kl_pq


def fixed_population_birth_death_torch(q, score, sim, ancestors=None, fr_interval=None):
    """Fixed-population birth-death with a max_event_fraction safeguard.

    Over-represented replicas (positive score) die; each is replaced by a clone of
    an under-represented replica (birth weights ~ max(-S_i, 0)). Deaths are capped
    at max_event_fraction * n_replicas (random subselection). The whole replica
    configuration is copied. If ``ancestors`` (long tensor, one label per replica)
    is given it is updated in place-style and returned so lineage can be tracked.
    Returns (q_new, ancestors_new, stats) where stats holds replacement count and
    the z-free birth/death replica indices.
    """
    R = q.shape[0]
    # score is already mean-zero+clipped by fr_score_torch; recenter is now
    # idempotent so this is a defensive no-op (safe for any ad-hoc caller).
    score = recentered_clipped_score_torch(score, sim.score_clip)
    interval = sim.fr_every if fr_interval is None else fr_interval
    dt_eff = sim.dt * max(int(interval), 1)

    empty = torch.empty(0, dtype=torch.long, device=q.device)
    # max_event_fraction <= 0 (or so small that floor(frac*R) == 0) disables all
    # events, matching the convention "fraction -> 0 ⇒ zero replacements".
    max_events = int(sim.max_event_fraction * R)
    no_event = {"replacement": 0, "death_idx": empty, "birth_src": empty}
    if max_events < 1:
        return q.clone(), (None if ancestors is None else ancestors.clone()), no_event

    death_weights = torch.clamp(score, min=0.0)
    birth_weights = torch.clamp(-score, min=0.0)
    death_mass = torch.sum(death_weights)
    birth_mass = torch.sum(birth_weights)
    if bool(((death_mass <= EPS) | (birth_mass <= EPS)).detach().cpu().item()):
        return q.clone(), (None if ancestors is None else ancestors.clone()), no_event

    death_prob = torch.where(
        death_weights > 0.0,
        1.0 - torch.exp(-sim.fr_rate * death_weights * dt_eff),
        torch.zeros_like(death_weights),
    )
    death_indices = torch.nonzero(
        torch.rand(R, device=q.device, dtype=q.dtype) < death_prob, as_tuple=False
    ).flatten()
    n_events = int(death_indices.numel())
    if n_events == 0:
        return q.clone(), (None if ancestors is None else ancestors.clone()), no_event

    if n_events > max_events:
        perm = torch.randperm(n_events, device=q.device)[:max_events]
        death_indices = death_indices[perm]
        n_events = max_events

    clone_sources = torch.multinomial(birth_weights, n_events, replacement=True)
    q_new = q.clone()
    q_new[death_indices] = q.index_select(0, clone_sources)
    ancestors_new = None
    if ancestors is not None:
        ancestors_new = ancestors.clone()
        ancestors_new[death_indices] = ancestors.index_select(0, clone_sources)
    return q_new, ancestors_new, {
        "replacement": int(n_events), "death_idx": death_indices, "birth_src": clone_sources}


# ---------------------------------------------------------------------------
# Diagnostics helpers
# ---------------------------------------------------------------------------
def to_numpy(x):
    if isinstance(x, torch.Tensor):
        return x.detach().cpu().numpy()
    return np.asarray(x)


def profile_l2_error_np(profile, reference, grid, mask=None):
    p = np.asarray(profile, dtype=float)
    r = np.asarray(reference, dtype=float)
    g = np.asarray(grid, dtype=float)
    if mask is not None:
        p, r, g = p[mask], r[mask], g[mask]
    return math.sqrt(np.trapezoid((p - r) ** 2, g) / (g[-1] - g[0]))


def align_additive_constant_np(profile, reference, grid, mask=None):
    profile = np.asarray(profile, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if mask is None:
        return profile - np.mean(profile - reference)
    return profile - np.mean(profile[mask] - reference[mask])


def region_masks_np(grid, sim):
    g = np.asarray(grid, dtype=float)
    compact = g < sim.transition_lo
    transition = (g >= sim.transition_lo) & (g <= sim.transition_hi)
    stretched = g > sim.transition_hi
    return {"compact": compact, "transition": transition, "stretched": stretched}


def eval_window_mask_np(grid, sim):
    g = np.asarray(grid, dtype=float)
    return (g >= sim.eval_z_lo) & (g <= sim.eval_z_hi)


def region_fractions_torch(z, sim):
    """Fraction of replicas in compact/transition/stretched regions."""
    total = float(z.numel())
    compact = float((z < sim.transition_lo).sum().item()) / total
    transition = float(((z >= sim.transition_lo) & (z <= sim.transition_hi)).sum().item()) / total
    stretched = float((z > sim.transition_hi).sum().item()) / total
    return compact, transition, stretched


def ancestor_ess(ancestors, n_replicas):
    """ESS_ancestor = 1 / sum_a w_a^2 where w_a is the fraction descended from a."""
    if ancestors is None:
        return float("nan"), 0
    counts = torch.bincount(ancestors, minlength=n_replicas).to(torch.float64)
    w = counts / counts.sum().clamp_min(1.0)
    ess = 1.0 / torch.clamp((w * w).sum(), min=EPS)
    n_unique = int((counts > 0).sum().item())
    return float(ess.item()), n_unique


def region_l2_errors(profile, reference, grid, sim):
    """L2 of (profile-reference) over compact/transition/stretched regions.

    Each region mask is intersected with the interior evaluation window so all
    three regions share the SAME interior domain as the headline l2_f / l2_fp
    (otherwise compact/stretched would include wall-affected edge bins that the
    headline metric deliberately excludes). profile/reference are arrays on `grid`.
    """
    masks = region_masks_np(grid, sim)
    eval_mask = eval_window_mask_np(grid, sim)
    out = {}
    for name, m in masks.items():
        m = m & eval_mask
        if m.sum() < 2:
            out[name] = float("nan")
            continue
        out[name] = profile_l2_error_np(profile, reference, grid, mask=m)
    return out


# ---------------------------------------------------------------------------
# No-oracle-leakage audit
# ---------------------------------------------------------------------------
def assert_no_oracle_leakage(method, oracle_free_energy):
    """Hard guard: only fr_oracle may receive an oracle/TI reference target.

    fr_estimated (the practical method), abf, and fr_uniform MUST NOT receive it.
    """
    if method == "fr_oracle":
        if oracle_free_energy is None:
            raise ValueError("fr_oracle requires oracle_free_energy (the TI reference).")
        return
    if oracle_free_energy is not None:
        raise AssertionError(
            f"NO-ORACLE-LEAKAGE VIOLATION: method='{method}' received oracle_free_energy. "
            "Only 'fr_oracle' may use the TI reference as a target.")


# ---------------------------------------------------------------------------
# Main sampler. abf / fr_estimated / fr_uniform / fr_oracle.
# ---------------------------------------------------------------------------
@torch.inference_mode()
def run_sampler_gpu(method, params, sim, engine, initial_q=None,
                    oracle_free_energy=None, collect_diagnostics=True, verbose=True):
    """Run one sampler on the GPU.

    method in {abf, fr_estimated, fr_uniform, fr_oracle}. The FR target is built
    per-variant; only fr_oracle reads oracle_free_energy (the TI reference). The
    EMA estimated target is maintained ONLY for fr_estimated and never sees TI.
    """
    if method not in ALL_METHODS:
        raise ValueError(f"method must be one of {ALL_METHODS}, got {method!r}")
    assert_no_oracle_leakage(method, oracle_free_energy)
    is_fr = method in FR_METHODS

    torch.manual_seed(sim.seed)
    q = (lattice_initial_conditions(params, sim.n_replicas, engine.device, engine.dtype, seed=sim.seed)
         if initial_q is None else initial_q.clone())
    grid = torch.linspace(sim.z_min, sim.z_max, sim.n_grid, device=engine.device, dtype=engine.dtype)
    oracle_t = None
    if method == "fr_oracle":
        oracle_t = torch.as_tensor(oracle_free_energy, device=engine.device, dtype=engine.dtype)
        oracle_t = normalize_profile_zero_at_midpoint_torch(oracle_t, grid)

    bias_estimator = TorchKernelABFEstimator(grid, sim.abf_bandwidth, sim.abf_smooth_sigma, edge_extrapolate=sim.abf_edge_extrapolate)
    production_estimator = TorchKernelABFEstimator(grid, sim.abf_bandwidth, sim.abf_smooth_sigma, edge_extrapolate=sim.abf_edge_extrapolate)
    noise_scale = math.sqrt(2.0 * sim.dt / params.beta)
    total_replacement_events = 0
    F_target_ema = None
    ancestors = (torch.arange(sim.n_replicas, device=engine.device, dtype=torch.long) if is_fr else None)

    # birth/death z-location histograms accumulated over the run
    birth_hist = np.zeros(sim.n_grid, dtype=np.float64)
    death_hist = np.zeros(sim.n_grid, dtype=np.float64)

    diag = {k: [] for k in ["steps", "times", "mean_force", "pmf",
                             "p_hat", "q_target", "pq_l2", "kl_pq", "eff_counts",
                             "frac_compact", "frac_transition", "frac_stretched",
                             "ancestor_ess", "n_unique_ancestor", "repl_cumulative"]}

    t0 = time.perf_counter()
    for step in range(sim.n_steps + 1):
        forces_raw = engine.force(q, compute_energy=False)
        forces_physical = clip_forces(forces_raw, params.force_clip)
        z = reaction_coordinate(q, params)

        mean_force_input = forces_physical if sim.use_clipped_force_for_mean_force else forces_raw
        f_local = local_mean_force(q, mean_force_input, params)
        f_local = torch.clamp(f_local, -sim.mean_force_sample_clip, sim.mean_force_sample_clip)
        bias_estimator.update(z, f_local)
        if step >= sim.estimator_burn_in_steps:
            production_estimator.update(z, f_local)
        ramp = min(1.0, step / max(sim.abf_warmup_steps, 1))
        abf_scale = sim.abf_bias_scale * ramp
        abf_at_z = abf_scale * torch.clamp(bias_estimator.evaluate(z), -sim.abf_force_clip, sim.abf_force_clip)
        A_hat = bias_estimator.pmf_profile()                 # unscaled A_hat_n(z)
        current_bias_profile = abf_scale * A_hat             # B_n(z)
        # Online EMA target maintained ONLY for fr_estimated, starting at fr_start.
        if method == "fr_estimated" and (step + 1) >= sim.fr_start_steps:
            if F_target_ema is None:
                F_target_ema = A_hat.clone()
            else:
                r = sim.target_ema_rate
                F_target_ema = (1.0 - r) * F_target_ema + r * A_hat
            F_target_ema = normalize_profile_zero_at_midpoint_torch(F_target_ema, grid)

        transport = clip_forces(add_abf_force(q, forces_physical, abf_at_z, params), params.force_clip)
        transport = clip_forces(add_reaction_coordinate_wall_force(q, transport, z, sim, params), params.force_clip)

        if step % sim.save_every == 0 or step == sim.n_steps:
            report_estimator = production_estimator if production_estimator.n_updates > 0 else bias_estimator
            diag["steps"].append(step)
            diag["times"].append(step * sim.dt)
            diag["mean_force"].append(to_numpy(report_estimator.mean_force_profile()))
            diag["pmf"].append(to_numpy(report_estimator.pmf_profile()))
            diag["repl_cumulative"].append(total_replacement_events)
            if collect_diagnostics:
                fc, ft, fs = region_fractions_torch(z, sim)
                diag["frac_compact"].append(fc)
                diag["frac_transition"].append(ft)
                diag["frac_stretched"].append(fs)
                diag["eff_counts"].append(to_numpy(report_estimator.effective_counts()))
                p_grid = normalize_density_on_grid_torch(
                    kde_1d_torch(grid, z, sim.kde_bandwidth, sim.z_min, sim.z_max), grid)
                diag["p_hat"].append(to_numpy(p_grid))
                # build the *current* FR target for logging (no dynamics effect here)
                q_grid = _build_fr_target(method, grid, F_target_ema, current_bias_profile,
                                          oracle_t, params.beta)
                if q_grid is not None:
                    log_ratio = torch.log(torch.clamp(p_grid, min=EPS)) - torch.log(torch.clamp(q_grid, min=EPS))
                    diag["q_target"].append(to_numpy(q_grid))
                    diag["pq_l2"].append(float(math.sqrt(trapz_torch((p_grid - q_grid) ** 2, grid).item())))
                    diag["kl_pq"].append(float(trapz_torch(p_grid * log_ratio, grid).item()))
                else:
                    diag["q_target"].append(np.full(sim.n_grid, np.nan))
                    diag["pq_l2"].append(float("nan"))
                    diag["kl_pq"].append(float("nan"))
                if is_fr:
                    ess, nuq = ancestor_ess(ancestors, sim.n_replicas)
                    diag["ancestor_ess"].append(ess)
                    diag["n_unique_ancestor"].append(nuq)
                else:
                    diag["ancestor_ess"].append(float("nan"))
                    diag["n_unique_ancestor"].append(sim.n_replicas)

        if step == sim.n_steps:
            break

        q = wrap_positions(q + sim.dt * transport + noise_scale * torch.randn_like(q), params.box_length)
        if is_fr:
            next_step = step + 1
            do_fr = (next_step >= sim.fr_start_steps
                     and (next_step - sim.fr_start_steps) % max(int(sim.fr_every), 1) == 0)
            if do_fr:
                z_new = reaction_coordinate(q, params)
                q_grid = _build_fr_target(method, grid, F_target_ema, current_bias_profile,
                                          oracle_t, params.beta)
                if q_grid is not None:
                    score, p_fr, q_fr, kl_pq = fr_score_torch(z_new, grid, sim, q_grid)
                    q, ancestors, stats = fixed_population_birth_death_torch(
                        q, score, sim, ancestors=ancestors, fr_interval=sim.fr_every)
                    total_replacement_events += stats["replacement"]
                    if collect_diagnostics and stats["replacement"] > 0:
                        # record births/deaths at their PRE-replacement z (z_new).
                        zb = z_new.index_select(0, stats["birth_src"])
                        zd = z_new.index_select(0, stats["death_idx"])
                        birth_hist += np.histogram(to_numpy(zb), bins=sim.n_grid,
                                                   range=(sim.z_min, sim.z_max))[0]
                        death_hist += np.histogram(to_numpy(zd), bins=sim.n_grid,
                                                   range=(sim.z_min, sim.z_max))[0]

    diag["runtime_seconds"] = time.perf_counter() - t0
    diag["method"] = method
    diag["grid"] = to_numpy(grid)
    diag["total_replacement_events"] = total_replacement_events
    diag["F_target_ema"] = to_numpy(F_target_ema) if F_target_ema is not None else None
    diag["birth_hist"] = birth_hist
    diag["death_hist"] = death_hist
    diag["hist_edges"] = np.linspace(sim.z_min, sim.z_max, sim.n_grid + 1)
    for key in ["steps", "times", "mean_force", "pmf", "p_hat", "q_target", "eff_counts",
                "pq_l2", "kl_pq", "frac_compact", "frac_transition", "frac_stretched",
                "ancestor_ess", "n_unique_ancestor", "repl_cumulative"]:
        diag[key] = np.asarray(diag[key])
    if verbose:
        extra = f", replacements {total_replacement_events}" if is_fr else ""
        print(f"{method:13s}: {diag['runtime_seconds']:.1f}s{extra}")
    return diag


def _build_fr_target(method, grid, F_target_ema, current_bias_profile, oracle_t, beta):
    """Construct q_n(z) for the given FR method, or None for abf / not-yet-started."""
    if method == "abf":
        return None
    if method == "fr_uniform":
        return fr_target_uniform_torch(grid)
    if method == "fr_estimated":
        if F_target_ema is None:
            return None
        return fr_target_estimated_torch(grid, F_target_ema, current_bias_profile, beta)
    if method == "fr_oracle":
        return fr_target_oracle_torch(grid, oracle_t, current_bias_profile, beta)
    raise ValueError(method)


# ---------------------------------------------------------------------------
# TI reference (EVALUATION ONLY) + cache
# ---------------------------------------------------------------------------
@torch.inference_mode()
def constrained_ti_reference_gpu(params, sim, ti, engine, eval_grid_np=None, verbose=True):
    torch.manual_seed(ti.seed)
    z_all = torch.linspace(ti.z_min, ti.z_max, ti.n_z, device=engine.device, dtype=engine.dtype)
    mf = torch.zeros(ti.n_z, device=engine.device, dtype=engine.dtype)
    count = torch.zeros(ti.n_z, device=engine.device, dtype=engine.dtype)
    noise_scale = math.sqrt(2.0 * ti.dt / params.beta)

    t0 = time.perf_counter()
    for start in range(0, ti.n_z, ti.z_chunk_size):
        end = min(start + ti.z_chunk_size, ti.n_z)
        z_chunk = z_all[start:end]
        C = z_chunk.numel()
        q0 = lattice_initial_conditions(params, C * ti.n_replicas, engine.device, engine.dtype, seed=ti.seed + start)
        z_batch = z_chunk.repeat_interleave(ti.n_replicas)
        q = project_dimer_to_z(q0, z_batch, params)
        total_steps = ti.n_thermalization + ti.n_steps
        for step in range(total_steps):
            forces_raw = engine.force(q, compute_energy=False)
            transport = clip_forces(forces_raw, params.force_clip)
            q = wrap_positions(q + ti.dt * transport + noise_scale * torch.randn_like(q), params.box_length)
            q = project_dimer_to_z(q, z_batch, params)
            if step >= ti.n_thermalization and (step - ti.n_thermalization) % ti.sample_every == 0:
                sample_forces = engine.force(q, compute_energy=False)
                if sim.use_clipped_force_for_mean_force:
                    sample_forces = clip_forces(sample_forces, params.force_clip)
                f_local = local_mean_force(q, sample_forces, params).view(C, ti.n_replicas)
                mf[start:end] += f_local.sum(dim=1)
                count[start:end] += ti.n_replicas
        if verbose:
            print(f"TI chunk {start:02d}:{end:02d} done, elapsed {(time.perf_counter()-t0)/60:.1f} min")

    mean_force = mf / torch.clamp(count, min=1.0)
    mean_force_smooth = smooth_profile_torch(mean_force, ti.smooth_sigma)
    free_energy = normalize_profile_zero_at_midpoint_torch(cumulative_trapezoid_torch(mean_force_smooth, z_all), z_all)

    if eval_grid_np is None:
        eval_grid = torch.linspace(sim.z_min, sim.z_max, sim.n_grid, device=engine.device, dtype=engine.dtype)
    else:
        eval_grid = torch.as_tensor(eval_grid_np, device=engine.device, dtype=engine.dtype)
    mf_eval = interp_uniform_grid(mean_force_smooth, z_all, eval_grid, outside_value=0.0)
    fe_eval = interp_uniform_grid(free_energy, z_all, eval_grid, outside_value=float(free_energy[0].item()))
    fe_eval = normalize_profile_zero_at_midpoint_torch(fe_eval, eval_grid)
    return {
        "label": "thermodynamic integration reference (evaluation only)",
        "grid": to_numpy(eval_grid),
        "mean_force": to_numpy(mf_eval),
        "free_energy": to_numpy(fe_eval),
        "z_ti": to_numpy(z_all),
        "runtime_seconds": time.perf_counter() - t0,
    }


def load_or_compute_ti_reference(path, params, sim, ti, engine, verbose=True):
    """Load the cached TI reference if the grid matches, else compute and cache it.

    EVALUATION ONLY. Never feeds the Fisher-Rao target (except fr_oracle diagnostic).
    """
    eval_grid_np = np.linspace(sim.z_min, sim.z_max, sim.n_grid)
    if os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        grid_ok = (d["grid"].shape[0] == sim.n_grid
                   and abs(float(d["grid"][0]) - sim.z_min) < 1e-5
                   and abs(float(d["grid"][-1]) - sim.z_max) < 1e-5)
        if grid_ok:
            if verbose:
                print(f"Loaded cached TI reference from {path}")
            return {"label": str(d["label"]), "grid": d["grid"],
                    "mean_force": d["mean_force"], "free_energy": d["free_energy"], "z_ti": d["z_ti"]}
        if verbose:
            print(f"Cache {path} grid mismatch; recomputing.")
    ref = constrained_ti_reference_gpu(params, sim, ti, engine, eval_grid_np=eval_grid_np, verbose=verbose)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, label=ref["label"], grid=ref["grid"], mean_force=ref["mean_force"],
             free_energy=ref["free_energy"], z_ti=ref["z_ti"])
    if verbose:
        print(f"Saved TI reference to {path} (runtime {ref['runtime_seconds']/60:.1f} min)")
    return ref


# ---------------------------------------------------------------------------
# Evaluation metrics vs reference (additive-constant aligned on interior window)
# ---------------------------------------------------------------------------
def final_l2_errors(diag, reference, sim):
    grid = reference["grid"]
    mask = eval_window_mask_np(grid, sim)
    mf, fe = diag["mean_force"][-1], diag["pmf"][-1]
    fe_al = align_additive_constant_np(fe, reference["free_energy"], grid, mask=mask)
    l2_fp = profile_l2_error_np(mf, reference["mean_force"], grid, mask=mask)
    l2_f = profile_l2_error_np(fe_al, reference["free_energy"], grid, mask=mask)
    reg_f = region_l2_errors(fe_al, reference["free_energy"], grid, sim)
    reg_fp = region_l2_errors(mf, reference["mean_force"], grid, sim)
    return {"l2_f": l2_f, "l2_fp": l2_fp,
            "l2_f_compact": reg_f["compact"], "l2_f_transition": reg_f["transition"],
            "l2_f_stretched": reg_f["stretched"],
            "l2_fp_compact": reg_fp["compact"], "l2_fp_transition": reg_fp["transition"],
            "l2_fp_stretched": reg_fp["stretched"]}


def timeseries_l2(diag, reference, sim):
    """L2(F) and L2(F') at every saved time; plus integrated L2(F) over time."""
    grid = reference["grid"]
    mask = eval_window_mask_np(grid, sim)
    times = np.asarray(diag["times"], dtype=float)
    l2_f_t, l2_fp_t = [], []
    for k in range(len(times)):
        fe_al = align_additive_constant_np(diag["pmf"][k], reference["free_energy"], grid, mask=mask)
        l2_f_t.append(profile_l2_error_np(fe_al, reference["free_energy"], grid, mask=mask))
        l2_fp_t.append(profile_l2_error_np(diag["mean_force"][k], reference["mean_force"], grid, mask=mask))
    l2_f_t = np.asarray(l2_f_t)
    l2_fp_t = np.asarray(l2_fp_t)
    integrated_f = float(np.trapezoid(l2_f_t, times)) if len(times) > 1 else float("nan")
    return {"times": times, "l2_f_t": l2_f_t, "l2_fp_t": l2_fp_t, "integrated_l2_f": integrated_f}
