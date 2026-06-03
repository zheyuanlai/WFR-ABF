#!/usr/bin/env python3
"""Validate the torch backend before spending GPU hours.

Checks (all on tiny, fast sims):

  1. Analytic torch gradient vs finite differences (+ torch-vs-numpy potential).
  2. numpy-CPU engine vs torch-CPU engine (metrics same order of magnitude).
  3. torch-CPU vs torch-CUDA (only if CUDA is available locally).
  4. Batch-size invariance: batch_size_configs = 1 vs 4 (same order).
  5. binned_smooth estimator vs kernel_reference estimator (binning error small).

Strict gates: no NaN/inf anywhere; FR event fraction <= max_event_fraction.
Stochastic gates use generous "same order of magnitude" tolerances because RNG
and GPU math differ between backends (exact trajectory equality is NOT required).

Outputs:
  results/two_dim_xi_x/validation/validation_report.json
  results/two_dim_xi_x/validation/validation_summary.csv

Prints PASS/FAIL and exits non-zero if any check fails.

Example:
  python scripts/validate_torch_backend.py --config configs/two_dim_xi_x_smoke_gpu.yaml
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from abffr import (io_utils, metrics, potentials, reference,  # noqa: E402
                   simulation, simulation_torch, torch_utils as tu)
from abffr.io_utils import RunSpec, make_rng_streams  # noqa: E402

# Tiny sim overrides used for all stochastic checks (fast, deterministic shape).
TINY = dict(n_particles=200, n_steps=1500, eval_every=300)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True)
    p.add_argument("--device", default=None, help="Override device (cpu/cuda).")
    p.add_argument("--output-root", default=None)
    p.add_argument("--n-points", type=int, default=20000,
                   help="Random points for the gradient finite-difference check.")
    return p.parse_args(argv)


def _tiny_cfg(cfg):
    c = {k: (dict(v) if isinstance(v, dict) else v) for k, v in cfg.items()}
    c["simulation"] = dict(cfg["simulation"])
    c["simulation"].update(TINY)
    return c


def _l2_profile(a, b, x_grid, mask, center=False):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if center:
        a = a - a[mask].mean(); b = b - b[mask].mean()
    xa = x_grid[mask]
    return float(np.sqrt(np.trapezoid((a[mask] - b[mask]) ** 2, xa) / (xa[-1] - xa[0])))


# --------------------------------------------------------------------------- #
# Check 1: gradient
# --------------------------------------------------------------------------- #
def check_gradient(cfg, device, dtype, n_points, checks):
    d = cfg["domain"]
    g = tu.make_generator(12345, device)
    x = (torch.rand(n_points, generator=g, device=device, dtype=torch.float64)
         * (d["x_max"] - d["x_min"]) + d["x_min"])
    y = (torch.rand(n_points, generator=g, device=device, dtype=torch.float64)
         * (d["y_max"] - d["y_min"]) + d["y_min"])
    dVdx_a, dVdy_a = potentials.grad_potential_xy_torch(x, y)
    eps = 1e-5
    dVdx_fd = (potentials.potential_xy_torch(x + eps, y)
               - potentials.potential_xy_torch(x - eps, y)) / (2 * eps)
    dVdy_fd = (potentials.potential_xy_torch(x, y + eps)
               - potentials.potential_xy_torch(x, y - eps)) / (2 * eps)
    ex = (dVdx_a - dVdx_fd).abs(); ey = (dVdy_a - dVdy_fd).abs()
    # torch vs numpy potential/gradient at the same points.
    xn, yn = x.cpu().numpy(), y.cpu().numpy()
    dV_np = potentials.potential_xy(xn, yn)
    dnp = float(np.max(np.abs(dV_np - potentials.potential_xy_torch(x, y).cpu().numpy())))
    gx_np, gy_np = potentials.grad_potential_xy(xn, yn)
    dgx = float(np.max(np.abs(gx_np - dVdx_a.cpu().numpy())))
    dgy = float(np.max(np.abs(gy_np - dVdy_a.cpu().numpy())))

    max_fd = float(max(ex.max(), ey.max()))
    passed = (max_fd < 1e-3 and dnp < 1e-8 and max(dgx, dgy) < 1e-8
              and torch.isfinite(dVdx_a).all() and torch.isfinite(dVdy_a).all())
    checks.append(dict(
        name="gradient_finite_difference", passed=bool(passed),
        max_fd_error=max_fd, mean_fd_error_dx=float(ex.mean()),
        mean_fd_error_dy=float(ey.mean()),
        max_torch_vs_numpy_V=dnp, max_torch_vs_numpy_gradx=dgx,
        max_torch_vs_numpy_grady=dgy, tol="max_fd<1e-3, torch_vs_numpy<1e-8"))


# --------------------------------------------------------------------------- #
# Helpers to run one numpy / torch run and return final profiles + summary
# --------------------------------------------------------------------------- #
def _run_numpy(spec, cfg, x_grid, ref, ev):
    sim = cfg["simulation"]; abf = cfg["abf"]; fr = cfg["fr"]
    r_init, r_noise, r_fr = make_rng_streams(spec.seed)
    diag = simulation.run_simulation(
        target_type=spec.target_type, beta=float(sim["beta"]), dt=float(sim["dt"]),
        n_steps=int(sim["n_steps"]), n_particles=int(sim["n_particles"]),
        eval_every=int(sim["eval_every"]), x_grid=x_grid, F_ref=ref["F_ref"],
        Fprime_ref=ref["Fprime_ref"], domain=cfg["domain"], h=float(abf["h"]),
        eta=float(spec.eta), min_count=float(abf.get("min_count", 1.0)),
        ema_alpha=float(abf.get("ema_alpha", 0.05)), gamma=float(spec.gamma),
        burnin_fraction=float(spec.burnin_fraction),
        ramp_fraction=float(fr.get("ramp_fraction", 0.1)), fr_every=int(spec.fr_every),
        score_clip=fr.get("score_clip", 5.0),
        max_event_fraction=fr.get("max_event_fraction", 0.10),
        x_init_mode=sim.get("x_init_mode", "mixed"),
        y_init_mode=sim.get("y_init_mode", "mixed"), x_barrier=ev.x_barrier,
        rng_init=r_init, rng_noise=r_noise, rng_fr=r_fr)
    return diag, metrics.final_summary(diag, x_grid, ref["F_ref"], ref["Fprime_ref"], ev)


def _run_torch(specs, cfg, x_grid, ref, ev, device, dtype, estimator="binned_smooth"):
    res = simulation_torch.run_batch(
        specs, cfg=cfg, x_grid=x_grid, F_ref=ref["F_ref"],
        Fprime_ref=ref["Fprime_ref"], ev=ev, device=device, dtype=dtype,
        estimator=estimator, base_seed=0)
    summaries = [metrics.final_summary(d, x_grid, ref["F_ref"], ref["Fprime_ref"], ev)
                 for d in res.diags]
    return res, summaries


def _same_order(a, b, lo=0.2, hi=5.0):
    if not (np.isfinite(a) and np.isfinite(b)):
        return False
    if max(abs(a), abs(b)) < 1e-6:
        return True
    if min(abs(a), abs(b)) < 1e-9:
        return False
    r = abs(a) / abs(b)
    return lo <= r <= hi


# --------------------------------------------------------------------------- #
# Checks 2-5
# --------------------------------------------------------------------------- #
def check_numpy_vs_torch(cfg, x_grid, ref, ev, device, dtype, mask, checks):
    spec = RunSpec("abf_fr_estimated", "estimated", seed=0, gamma=0.05, eta=0.10,
                   burnin_fraction=0.0, fr_every=5)
    _, s_np = _run_numpy(spec, cfg, x_grid, ref, ev)
    res, [s_t] = _run_torch([spec], cfg, x_grid, ref, ev, device, dtype)
    finite = not (s_np["any_nan"] or s_t["any_nan"])
    same_F = _same_order(s_np["final_l2_F"], s_t["final_l2_F"])
    same_Fp = _same_order(s_np["final_l2_Fprime"], s_t["final_l2_Fprime"])
    checks.append(dict(
        name="numpy_vs_torch_cpu",
        passed=bool(finite and same_F and same_Fp),
        numpy_final_l2_F=s_np["final_l2_F"], torch_final_l2_F=s_t["final_l2_F"],
        numpy_final_l2_Fprime=s_np["final_l2_Fprime"],
        torch_final_l2_Fprime=s_t["final_l2_Fprime"],
        tol="same order of magnitude (0.2x-5x), no NaN"))


def check_cuda_consistency(cfg, x_grid, ref, ev, dtype, mask, checks):
    if not torch.cuda.is_available():
        checks.append(dict(name="torch_cpu_vs_cuda", passed=True, skipped=True,
                           note="CUDA unavailable locally; run this check on the "
                                "remote GPU."))
        return
    spec = RunSpec("abf_fr_estimated", "estimated", seed=0, gamma=0.05, eta=0.10,
                   burnin_fraction=0.0, fr_every=5)
    _, [s_cpu] = _run_torch([spec], cfg, x_grid, ref, ev, torch.device("cpu"), dtype)
    _, [s_gpu] = _run_torch([spec], cfg, x_grid, ref, ev, torch.device("cuda"), dtype)
    ok = (_same_order(s_cpu["final_l2_F"], s_gpu["final_l2_F"])
          and not s_cpu["any_nan"] and not s_gpu["any_nan"])
    checks.append(dict(name="torch_cpu_vs_cuda", passed=bool(ok),
                       cpu_final_l2_F=s_cpu["final_l2_F"],
                       gpu_final_l2_F=s_gpu["final_l2_F"],
                       tol="same order of magnitude, no NaN"))


def check_batch_invariance(cfg, x_grid, ref, ev, device, dtype, checks):
    specs = [RunSpec("abf_fr_estimated", "estimated", seed=s, gamma=0.05, eta=0.10,
                     burnin_fraction=0.0, fr_every=5) for s in range(4)]
    # Batch of 4.
    _, s4 = _run_torch(specs, cfg, x_grid, ref, ev, device, dtype)
    # Each run alone (batch of 1).
    s1 = [_run_torch([sp], cfg, x_grid, ref, ev, device, dtype)[1][0] for sp in specs]
    per = []
    ok = True
    for sp, a, b in zip(specs, s1, s4):
        same = _same_order(a["final_l2_F"], b["final_l2_F"], lo=0.1, hi=10.0)
        ok = ok and same and not a["any_nan"] and not b["any_nan"]
        per.append(dict(seed=sp.seed, bs1=a["final_l2_F"], bs4=b["final_l2_F"],
                        same_order=bool(same)))
    checks.append(dict(name="batch_size_invariance", passed=bool(ok), per_run=per,
                       tol="same order of magnitude per matched run (0.1x-10x)"))


def check_binned_vs_kernel(cfg, x_grid, ref, ev, device, dtype, mask, checks):
    """Static + dynamic comparison of binned_smooth vs the exact kernel estimator.

    Static: on one *fixed* particle cloud both estimators see identical data, so
    the profile difference is the pure discretisation (binning) error.  Dynamic:
    an ``abf_only`` run with each estimator should reach the same-order final
    L2(F).
    """
    import torch as _t

    from abffr import potentials as _pot

    G = len(x_grid)
    x_grid_t = _t.as_tensor(np.asarray(x_grid), device=device, dtype=dtype)
    dx = tu.grid_spacing(x_grid_t); x0 = float(x_grid[0])
    h = float(cfg["abf"]["h"]); min_count = float(cfg["abf"].get("min_count", 1.0))
    rng = np.random.default_rng(0)
    Xn = rng.uniform(ev.eval_x_min, ev.eval_x_max, (1, 400))
    Yn = rng.uniform(cfg["domain"]["y_min"], cfg["domain"]["y_max"], (1, 400))
    X = _t.as_tensor(Xn, device=device, dtype=dtype)
    Y = _t.as_tensor(Yn, device=device, dtype=dtype)
    num_k, den_k = simulation_torch._kernel_estimator(x_grid_t, X, Y, h)
    Fp_k = (num_k / (den_k + min_count)).cpu().numpy()[0]
    k_h, r_h = tu.gaussian_kernel1d(h, dx, device, dtype)
    idx = tu.nearest_index(X, x0, dx, G)
    C = tu.scatter_grid(idx, G)
    S = tu.scatter_grid(idx, G, _pot.dVdx_xy_torch(X, Y))
    Fp_b = (tu.smooth_grid(S, k_h, r_h, dx)
            / (tu.smooth_grid(C, k_h, r_h, dx) + min_count)).cpu().numpy()[0]
    static_diff = _l2_profile(Fp_b, Fp_k, x_grid, mask)
    static_scale = _l2_profile(Fp_k, np.zeros_like(Fp_k), x_grid, mask) + 1e-9
    static_rel = static_diff / static_scale

    # Dynamic same-order check.
    spec = RunSpec("abf_only", "none", seed=0, gamma=0.0, eta=0.10,
                   burnin_fraction=0.0, fr_every=1)
    _, [s_b] = _run_torch([spec], cfg, x_grid, ref, ev, device, dtype, "binned_smooth")
    _, [s_k] = _run_torch([spec], cfg, x_grid, ref, ev, device, dtype, "kernel_reference")
    final_close = _same_order(s_b["final_l2_F"], s_k["final_l2_F"], lo=0.2, hi=5.0)
    ok = (static_rel < 0.10 and final_close
          and np.isfinite(Fp_b).all() and not s_b["any_nan"] and not s_k["any_nan"])
    checks.append(dict(
        name="binned_smooth_vs_kernel_reference", passed=bool(ok),
        static_relative_Fprime_l2=static_rel, static_Fprime_l2_diff=static_diff,
        binned_final_l2_F=s_b["final_l2_F"], kernel_final_l2_F=s_k["final_l2_F"],
        tol="static relative Fprime L2 < 0.10; final L2(F) same order"))


def check_event_fraction(cfg, x_grid, ref, ev, device, dtype, checks):
    cap = float(cfg["fr"].get("max_event_fraction", 0.10))
    specs = [RunSpec("abf_fr_uniform", "uniform", seed=0, gamma=g, eta=0.10,
                     burnin_fraction=0.0, fr_every=5) for g in (0.1, 0.5)]
    res, summaries = _run_torch(specs, cfg, x_grid, ref, ev, device, dtype)
    worst = max(s["max_fr_event_fraction"] for s in summaries)
    any_nan = any(s["any_nan"] for s in summaries)
    checks.append(dict(name="fr_event_fraction_cap", passed=bool(worst <= cap + 1e-6
                                                                 and not any_nan),
                       max_event_fraction=worst, cap=cap,
                       tol="max_fr_event_fraction <= cap + 1e-6"))


def main(argv=None):
    args = parse_args(argv)
    cfg = io_utils.load_config(args.config)
    if args.output_root:
        cfg["output_root"] = args.output_root
    if args.device:
        cfg["device"] = args.device
    cfg = _tiny_cfg(cfg)

    device = tu.resolve_device(cfg.get("device"))
    dtype = tu.resolve_dtype(cfg.get("dtype"))
    cuda_local = torch.cuda.is_available()
    print(f"[validate] device={device} dtype={dtype} cuda_available={cuda_local}")

    # Tiny on-grid reference (no CSV gate: validation must run standalone).
    x_grid, ref, _ = reference.load_reference_for_run(cfg, require_csv=False)
    ev = metrics.EvalConfig.from_domain(cfg["domain"])
    mask = ev.eval_mask(x_grid)

    checks = []
    check_gradient(cfg, device, dtype, args.n_points, checks)
    check_numpy_vs_torch(cfg, x_grid, ref, ev, device, dtype, mask, checks)
    check_cuda_consistency(cfg, x_grid, ref, ev, dtype, mask, checks)
    check_batch_invariance(cfg, x_grid, ref, ev, device, dtype, checks)
    check_binned_vs_kernel(cfg, x_grid, ref, ev, device, dtype, mask, checks)
    check_event_fraction(cfg, x_grid, ref, ev, device, dtype, checks)

    all_pass = all(c["passed"] for c in checks)
    out_dir = io_utils.ensure_dir(os.path.join(io_utils.output_root(cfg), "validation"))
    report = dict(device=str(device), dtype=str(dtype), cuda_available=cuda_local,
                  config=args.config, all_pass=all_pass, checks=checks)
    io_utils.save_json(os.path.join(out_dir, "validation_report.json"), report)
    summ = pd.DataFrame([{"check": c["name"], "passed": c["passed"],
                          "skipped": c.get("skipped", False),
                          "detail": c.get("tol", "")} for c in checks])
    summ.to_csv(os.path.join(out_dir, "validation_summary.csv"), index=False)

    print("\n==================== VALIDATION ====================")
    for c in checks:
        tag = "SKIP" if c.get("skipped") else ("PASS" if c["passed"] else "FAIL")
        print(f"  [{tag}] {c['name']}")
        for k, v in c.items():
            if k in ("name", "passed", "skipped", "per_run"):
                continue
            print(f"           {k}: {v}")
        if "per_run" in c:
            for r in c["per_run"]:
                print(f"           {r}")
    print("====================================================")
    print(f"OVERALL: {'PASS' if all_pass else 'FAIL'}")
    if not cuda_local:
        print("NOTE: CUDA was unavailable locally -- the CPU-vs-CUDA check was "
              "skipped. Re-run this script on the remote GPU before production.")
    print(f"Report: {os.path.relpath(os.path.join(out_dir, 'validation_report.json'))}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
