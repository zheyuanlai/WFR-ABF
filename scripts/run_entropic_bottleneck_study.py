"""Run the entropic-bottleneck ABF+FR mechanism study.

Single-GPU batched execution: each stage builds a list of (config, seed) rows
and runs them — across all methods — in wide `simulate_batch` calls, chunked
only to bound GPU memory.  One H200 is plenty; no multi-GPU machinery.

Stages
------
  0  reproduce  : notebook setting, 5 seeds, [abf, fr_estimated, fr_uniform, fr_oracle]
  1  seeds      : notebook setting, 20 seeds, [abf, fr_estimated]
  2  omega      : omega_in in {1,5,10,25,50}, 10 seeds, [abf, fr_estimated]
  3  beta       : beta in {2,4,8,12}, 10 seeds, [abf, fr_estimated]
  4  gamma      : gamma in {1,3,5,10,15,25,50}, 10 seeds, [abf, fr_estimated]

Raw output: results/entropic_bottleneck/raw/<stage>/<runkey>.npz  (idempotent;
skipped if a valid file exists unless --overwrite).

Usage
-----
  CUDA_VISIBLE_DEVICES=2 python -u scripts/run_entropic_bottleneck_study.py --stage all
  CUDA_VISIBLE_DEVICES=2 python -u scripts/run_entropic_bottleneck_study.py --stage 2
  # optional 2-GPU sharding of (config,seed) rows:
  CUDA_VISIBLE_DEVICES=2 python -u ... --stage 4 --shard 0 --num-shards 2
  CUDA_VISIBLE_DEVICES=3 python -u ... --stage 4 --shard 1 --num-shards 2
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import asdict, replace

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import eb_abffr_core as eb  # noqa: E402

RAW_ROOT = os.path.join(os.path.dirname(__file__), "..", "results", "entropic_bottleneck", "raw")

BASE = eb.PhysConfig()  # notebook Stage 0 defaults

# method bundles
MAIN_METHODS = [eb.ABF, eb.FR_ESTIMATED]
ALL_METHODS = [eb.ABF, eb.FR_ESTIMATED, eb.FR_UNIFORM, eb.FR_ORACLE]


# -----------------------------------------------------------------------------
# stage definitions: each returns (stage_name, methods, list of (runkey, config, seed))
# -----------------------------------------------------------------------------
def stage_rows(stage):
    rows = []
    if stage == "0":
        methods = ALL_METHODS
        for sd in range(5):
            rows.append((f"omega_in25_beta8_gamma15_seed{sd}", BASE, sd))
        return "stage0_reproduce", methods, rows
    if stage == "1":
        methods = MAIN_METHODS
        for sd in range(20):
            rows.append((f"omega_in25_beta8_gamma15_seed{sd}", BASE, sd))
        return "stage1_seeds", methods, rows
    if stage == "2":
        methods = MAIN_METHODS
        # Lower dt (and raise n_steps to keep total time T=40 fixed) so EVERY
        # omega_in is integrated identically and stably.  The y-channel is an OU
        # process y <- y(1 - omega^2 dt) + noise; explicit Euler needs
        # |1 - omega^2 dt| < 1, i.e. omega < sqrt(2/dt).  At dt=1e-3 omega_in=50
        # diverges (mult = -1.5); dt=4e-4 gives omega^2 dt = 1.0 (stable) and a
        # uniform timestep across the whole sweep avoids a dt confound in the
        # headline mechanism figure.
        dt2, ns2 = 4e-4, 100000
        for oin in (1.0, 5.0, 10.0, 25.0, 50.0):
            cfg = replace(BASE, omega_in=oin, dt=dt2, n_steps=ns2, save_every=1000)
            for sd in range(10):
                rows.append((f"omega_in{oin:g}_beta8_gamma15_seed{sd}", cfg, sd))
        return "stage2_omega", methods, rows
    if stage == "3":
        methods = MAIN_METHODS
        for beta in (2.0, 4.0, 8.0, 12.0):
            cfg = replace(BASE, beta=beta)
            for sd in range(10):
                rows.append((f"omega_in25_beta{beta:g}_gamma15_seed{sd}", cfg, sd))
        return "stage3_beta", methods, rows
    if stage == "4":
        methods = MAIN_METHODS
        for gamma in (1.0, 3.0, 5.0, 10.0, 15.0, 25.0, 50.0):
            cfg = replace(BASE, gamma=gamma)
            for sd in range(10):
                rows.append((f"omega_in25_beta8_gamma{gamma:g}_seed{sd}", cfg, sd))
        return "stage4_gamma", methods, rows
    raise ValueError(f"unknown stage {stage}")


def run_path(stage_name, method_name, runkey):
    return os.path.join(RAW_ROOT, stage_name, f"{method_name}__{runkey}.npz")


def save_record(path, rec):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    flat = {}
    for k, v in rec.items():
        if k == "config":
            for ck, cv in v.items():
                flat[f"cfg__{ck}"] = cv
        elif v is None:
            flat[f"none__{k}"] = True
        else:
            flat[k] = v
    np.savez_compressed(path, **flat)


def valid_existing(path):
    if not os.path.exists(path):
        return False
    try:
        with np.load(path, allow_pickle=True) as d:
            return "final_l2_f" in d and np.isfinite(float(d["final_l2_f"]))
    except Exception:
        return False


# structural params that must be uniform within a batched call
STRUCT_KEYS = ("N", "dt", "n_steps", "save_every", "fr_every", "fr_burnin",
               "ramp_fraction", "h", "eta", "min_count", "ess_window_steps")


def struct_sig(cfg):
    return tuple(getattr(cfg, k) for k in STRUCT_KEYS)


def run_stage(stage, methods, rows, overwrite, max_rows_per_call, batch_seed_base,
              shard, num_shards):
    """Execute one stage. rows = list of (runkey, config, seed).

    All `methods` for a given (config, seed) row share init+noise (matched seed).
    Rows are grouped by structural signature, sharded, then chunked into wide
    `simulate_batch` calls.  Results saved per (method, runkey); existing valid
    files are skipped unless overwrite.
    """
    stage_name = stage
    # optional sharding across GPUs: take every num_shards-th row
    if num_shards > 1:
        rows = [r for i, r in enumerate(rows) if i % num_shards == shard]
        print(f"[shard {shard}/{num_shards}] {len(rows)} rows", flush=True)

    # which rows still need ANY method computed?
    pending = []
    for runkey, cfg, sd in rows:
        need = [m for m in methods if overwrite or not valid_existing(run_path(stage_name, m.name, runkey))]
        if need:
            pending.append((runkey, cfg, sd, need))
    n_done = len(rows) - len(pending)
    print(f"stage {stage_name}: {len(rows)} rows, {n_done} already complete, "
          f"{len(pending)} to run", flush=True)
    if not pending:
        return

    # group by (struct_sig, frozenset of needed-method names) so a batched call
    # is uniform in structure and method set
    groups = {}
    for runkey, cfg, sd, need in pending:
        key = (struct_sig(cfg), tuple(m.name for m in need))
        groups.setdefault(key, []).append((runkey, cfg, sd, need))

    call_idx = 0
    t_stage = time.time()
    for (sig, mnames), grp in groups.items():
        method_objs = [eb.METHOD_REGISTRY[n] for n in mnames]
        # chunk rows to bound memory: R = chunk * len(methods)
        for c0 in range(0, len(grp), max_rows_per_call):
            chunk = grp[c0:c0 + max_rows_per_call]
            configs = [cfg for (_, cfg, _, _) in chunk]
            seeds = [sd for (_, _, sd, _) in chunk]
            spec = eb.BatchSpec(configs=configs, seeds=seeds, methods=method_objs,
                                batch_seed=batch_seed_base + call_idx)
            t0 = time.time()
            recs = eb.simulate_batch(spec)
            if torch.cuda.is_available():
                torch.cuda.synchronize()
            # recs are ordered (b major, m minor) == (chunk row, method)
            M = len(method_objs)
            for bi, (runkey, _, _, _) in enumerate(chunk):
                for mi, m in enumerate(method_objs):
                    rec = recs[bi * M + mi]
                    save_record(run_path(stage_name, m.name, runkey), rec)
            dt_call = time.time() - t0
            print(f"  call {call_idx}: R={len(chunk)*M} "
                  f"({len(chunk)} rows x {M} methods) in {dt_call:.0f}s", flush=True)
            call_idx += 1
    print(f"stage {stage_name} done in {time.time()-t_stage:.0f}s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    help="0,1,2,3,4 or 'all' or comma list e.g. 1,2")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--max-rows-per-call", type=int, default=70,
                    help="(config,seed) rows per batched call; R = this * n_methods")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--num-shards", type=int, default=1)
    args = ap.parse_args()

    if args.stage == "all":
        stages = ["0", "1", "2", "3", "4"]
    else:
        stages = args.stage.split(",")

    print(f"device={eb.DEVICE} dtype={eb.DTYPE} "
          f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES','?')}", flush=True)
    for st in stages:
        stage_name, methods, rows = stage_rows(st.strip())
        run_stage(stage_name, methods, rows, args.overwrite,
                  args.max_rows_per_call, batch_seed_base=10000 + int(st.strip()) * 1000,
                  shard=args.shard, num_shards=args.num_shards)


if __name__ == "__main__":
    main()
