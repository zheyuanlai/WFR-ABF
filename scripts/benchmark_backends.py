#!/usr/bin/env python3
"""Benchmark the CPU/numpy and torch backends; recommend a batch size.

Measures throughput (particle-steps/second) for:
  1. the existing CPU/numpy engine (one config),
  2. torch on the active device, swept over batch_size_configs 1, 4, 8, 16, 32,
     and 64 (64 only if it fits in memory),
  3. torch CUDA if available.

On a machine without CUDA this still benchmarks numpy + torch-CPU and prints a
clear note that the CUDA benchmark must be run remotely.

Output:
  results/two_dim_xi_x/benchmark/benchmark_backends.csv
with columns:
  backend,device,batch_size,n_particles,n_steps,n_configs,
  runtime_seconds,steps_per_second,particle_steps_per_second,
  final_l2_F,final_l2_Fprime,max_gpu_memory_mb

Example (LOCAL smoke; keep n-steps small on CPU):
  python scripts/benchmark_backends.py --config configs/two_dim_xi_x_smoke_gpu.yaml \
      --n-steps 1000
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from abffr import (io_utils, metrics, reference, simulation,  # noqa: E402
                   simulation_torch, torch_utils as tu)
from abffr.io_utils import RunSpec, make_rng_streams  # noqa: E402


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True)
    p.add_argument("--device", default=None)
    p.add_argument("--dtype", default=None, choices=["float32", "float64"])
    p.add_argument("--n-steps", type=int, default=None,
                   help="Override n_steps (recommended on CPU to keep it quick).")
    p.add_argument("--n-particles", type=int, default=None)
    p.add_argument("--batch-sizes", default="1,4,8,16,32,64",
                   help="Comma-separated batch sizes to sweep.")
    p.add_argument("--skip-numpy", action="store_true",
                   help="Skip the (slow) numpy backend timing.")
    p.add_argument("--output-root", default=None)
    return p.parse_args(argv)


def _final(diag_or_summary):
    return diag_or_summary["final_l2_F"], diag_or_summary["final_l2_Fprime"]


def _reset_gpu_mem(device):
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)


def _gpu_mem_mb(device):
    if device.type == "cuda":
        return float(torch.cuda.max_memory_allocated(device) / 1e6)
    return float("nan")


def bench_numpy(cfg, x_grid, ref, ev):
    import time
    sim = cfg["simulation"]; abf = cfg["abf"]; fr = cfg["fr"]
    spec = RunSpec("abf_fr_estimated", "estimated", seed=0, gamma=0.05, eta=0.10,
                   burnin_fraction=0.0, fr_every=5)
    r_init, r_noise, r_fr = make_rng_streams(spec.seed)
    t0 = time.time()
    diag = simulation.run_simulation(
        target_type="estimated", beta=float(sim["beta"]), dt=float(sim["dt"]),
        n_steps=int(sim["n_steps"]), n_particles=int(sim["n_particles"]),
        eval_every=int(sim["eval_every"]), x_grid=x_grid, F_ref=ref["F_ref"],
        Fprime_ref=ref["Fprime_ref"], domain=cfg["domain"], h=float(abf["h"]),
        eta=0.10, min_count=float(abf.get("min_count", 1.0)),
        ema_alpha=float(abf.get("ema_alpha", 0.05)), gamma=0.05,
        ramp_fraction=float(fr.get("ramp_fraction", 0.1)), fr_every=5,
        score_clip=fr.get("score_clip", 5.0),
        max_event_fraction=fr.get("max_event_fraction", 0.10),
        x_init_mode=sim.get("x_init_mode", "mixed"),
        y_init_mode=sim.get("y_init_mode", "mixed"), x_barrier=ev.x_barrier,
        rng_init=r_init, rng_noise=r_noise, rng_fr=r_fr)
    runtime = time.time() - t0
    summ = metrics.final_summary(diag, x_grid, ref["F_ref"], ref["Fprime_ref"], ev)
    n_steps = int(sim["n_steps"]); n_p = int(sim["n_particles"])
    return dict(backend="numpy", device="cpu", batch_size=1, n_particles=n_p,
                n_steps=n_steps, n_configs=1, runtime_seconds=runtime,
                steps_per_second=n_steps / runtime,
                particle_steps_per_second=n_p * n_steps / runtime,
                final_l2_F=summ["final_l2_F"], final_l2_Fprime=summ["final_l2_Fprime"],
                max_gpu_memory_mb=float("nan"))


def bench_torch(cfg, x_grid, ref, ev, device, dtype, batch_size):
    specs = [RunSpec("abf_fr_estimated", "estimated", seed=s, gamma=0.05, eta=0.10,
                     burnin_fraction=0.0, fr_every=5) for s in range(batch_size)]
    _reset_gpu_mem(device)
    res = simulation_torch.run_batch(
        specs, cfg=cfg, x_grid=x_grid, F_ref=ref["F_ref"],
        Fprime_ref=ref["Fprime_ref"], ev=ev, device=device, dtype=dtype,
        estimator=cfg.get("abf", {}).get("estimator", "binned_smooth"), base_seed=0)
    summ = [metrics.final_summary(d, x_grid, ref["F_ref"], ref["Fprime_ref"], ev)
            for d in res.diags]
    n_steps = int(cfg["simulation"]["n_steps"]); n_p = int(cfg["simulation"]["n_particles"])
    rt = res.runtime_seconds
    return dict(backend="torch", device=str(device), batch_size=batch_size,
                n_particles=n_p, n_steps=n_steps, n_configs=batch_size,
                runtime_seconds=rt, steps_per_second=n_steps / rt,
                particle_steps_per_second=batch_size * n_p * n_steps / rt,
                final_l2_F=float(np.mean([s["final_l2_F"] for s in summ])),
                final_l2_Fprime=float(np.mean([s["final_l2_Fprime"] for s in summ])),
                max_gpu_memory_mb=_gpu_mem_mb(device))


def main(argv=None):
    args = parse_args(argv)
    cfg = io_utils.load_config(args.config)
    io_utils.apply_cli_overrides(cfg, device=args.device, dtype=args.dtype,
                                 n_steps=args.n_steps, n_particles=args.n_particles,
                                 output_root=args.output_root)
    device = tu.resolve_device(cfg.get("device"))
    dtype = tu.resolve_dtype(cfg.get("dtype"))
    print(f"[benchmark] device={device} dtype={dtype} "
          f"n_steps={cfg['simulation']['n_steps']} "
          f"n_particles={cfg['simulation']['n_particles']}")

    x_grid, ref, _ = reference.load_reference_for_run(cfg, require_csv=False)
    ev = metrics.EvalConfig.from_domain(cfg["domain"])
    rows = []

    if not args.skip_numpy:
        print("[benchmark] numpy backend (1 config)...")
        try:
            rows.append(bench_numpy(cfg, x_grid, ref, ev))
            print(f"   numpy: {rows[-1]['runtime_seconds']:.2f}s, "
                  f"{rows[-1]['particle_steps_per_second']:.3e} particle-steps/s")
        except Exception as exc:
            print(f"   numpy backend failed: {exc}")

    batch_sizes = [int(b) for b in args.batch_sizes.split(",")]
    for bs in batch_sizes:
        try:
            r = bench_torch(cfg, x_grid, ref, ev, device, dtype, bs)
            rows.append(r)
            mem = "" if np.isnan(r["max_gpu_memory_mb"]) else f", {r['max_gpu_memory_mb']:.0f} MB"
            print(f"   torch bs={bs:>3}: {r['runtime_seconds']:.2f}s, "
                  f"{r['particle_steps_per_second']:.3e} particle-steps/s{mem}")
        except RuntimeError as exc:
            if "out of memory" in str(exc).lower():
                print(f"   torch bs={bs}: OOM (skipping larger sizes)")
                if device.type == "cuda":
                    torch.cuda.empty_cache()
                break
            raise

    out_dir = io_utils.ensure_dir(os.path.join(io_utils.output_root(cfg), "benchmark"))
    df = pd.DataFrame(rows)
    cols = ["backend", "device", "batch_size", "n_particles", "n_steps", "n_configs",
            "runtime_seconds", "steps_per_second", "particle_steps_per_second",
            "final_l2_F", "final_l2_Fprime", "max_gpu_memory_mb"]
    df = df[cols]
    dest = os.path.join(out_dir, "benchmark_backends.csv")
    df.to_csv(dest, index=False)

    torch_rows = df[df["backend"] == "torch"]
    print("\n==================== BENCHMARK ====================")
    print(df.to_string(index=False))
    if not torch_rows.empty:
        best = torch_rows.loc[torch_rows["particle_steps_per_second"].idxmax()]
        print(f"\nRecommended batch_size_configs for this machine: "
              f"{int(best['batch_size'])} "
              f"({best['particle_steps_per_second']:.3e} particle-steps/s on "
              f"{best['device']})")
    if device.type != "cuda":
        print("NOTE: CUDA unavailable locally -- this is the CPU benchmark only. "
              "Run on the remote GPU to choose the production batch_size_configs.")
    print(f"Saved: {os.path.relpath(dest)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
