#!/usr/bin/env python3
"""Run WCA-dimer production jobs from a YAML config (or a single ad-hoc job).

One .npz per run is written to results/wca_production/raw/, so an interrupted
sweep never loses completed work. Runs are idempotent: a valid existing result
is skipped unless --overwrite. Sharding splits a stage across GPUs/processes.

Examples
--------
# whole stage on the visible GPU
CUDA_VISIBLE_DEVICES=1 python scripts/run_wca_production.py \
    --config configs/wca_production.yaml --stage main

# stage split across 7 GPUs (run each line on its own GPU)
CUDA_VISIBLE_DEVICES=1 python scripts/run_wca_production.py \
    --config configs/wca_production.yaml --stage main --shard 0 --num-shards 7
...
CUDA_VISIBLE_DEVICES=7 python scripts/run_wca_production.py \
    --config configs/wca_production.yaml --stage main --shard 6 --num-shards 7

# one ad-hoc job
CUDA_VISIBLE_DEVICES=1 python scripts/run_wca_production.py \
    --method fr_estimated --seed 42 --n-steps 250000 --n-replicas 1024 \
    --fr-rate 0.1 --target-ema-rate 0.005 --max-event-fraction 0.02
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import numpy as np  # noqa: E402
import wca_abffr_core as core  # noqa: E402
import wca_jobs  # noqa: E402
from wca_jobs import RunSpec  # noqa: E402


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="configs/wca_production.yaml")
    p.add_argument("--stage", default=None, help="stage name in the YAML stages block")
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--num-shards", type=int, default=1)
    p.add_argument("--overwrite", action="store_true", help="recompute even valid results")
    p.add_argument("--dry-run", action="store_true", help="list jobs, do not run")
    p.add_argument("--limit", type=int, default=None, help="cap number of jobs (debug)")
    p.add_argument("--verbose", action="store_true")
    # single ad-hoc job
    p.add_argument("--method", default=None, choices=list(core.ALL_METHODS))
    p.add_argument("--name", default=None)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n-steps", type=int, default=250000)
    p.add_argument("--n-replicas", type=int, default=1024)
    p.add_argument("--a", type=float, default=1.5)
    p.add_argument("--fr-rate", type=float, default=0.10)
    p.add_argument("--target-ema-rate", type=float, default=0.005)
    p.add_argument("--max-event-fraction", type=float, default=0.02)
    p.add_argument("--fr-every", type=int, default=5)
    p.add_argument("--fr-start-steps", type=int, default=20000)
    p.add_argument("--score-clip", type=float, default=2.0)
    p.add_argument("--save-every", type=int, default=2500)
    return p.parse_args(argv)


def single_spec(args) -> RunSpec:
    name = args.name or args.method
    return RunSpec(
        stage="single", name=name, method=args.method, seed=args.seed,
        n_steps=args.n_steps, n_replicas=args.n_replicas, a=args.a,
        fr_rate=args.fr_rate, target_ema_rate=args.target_ema_rate,
        max_event_fraction=args.max_event_fraction, fr_every=args.fr_every,
        fr_start_steps=args.fr_start_steps, score_clip=args.score_clip,
        save_every=args.save_every)


def main(argv=None):
    args = parse_args(argv)
    cfg = wca_jobs.load_yaml(args.config)
    raw_dir = os.path.join(cfg["output_root"], "raw")
    cache_dir = cfg.get("cache_dir", "cache")
    os.makedirs(raw_dir, exist_ok=True)

    if args.stage:
        specs = wca_jobs.expand_stage(cfg, args.stage)
    elif args.method:
        specs = [single_spec(args)]
    else:
        raise SystemExit("provide either --stage or --method")

    # deterministic order, then shard by stride
    specs = sorted(specs, key=lambda s: s.run_id())
    if args.num_shards > 1:
        specs = [s for i, s in enumerate(specs) if i % args.num_shards == args.shard]
    if args.limit:
        specs = specs[:args.limit]

    cvd = os.environ.get("CUDA_VISIBLE_DEVICES", "unset")
    print(f"[run] stage={args.stage} shard={args.shard}/{args.num_shards} "
          f"jobs={len(specs)} device={core.DEVICE} CUDA_VISIBLE_DEVICES={cvd}")

    if args.dry_run:
        for s in specs:
            path = wca_jobs.run_npz_path(raw_dir, s)
            status = "VALID" if wca_jobs.run_is_valid(path) else "todo"
            print(f"  [{status:5s}] {s.run_id()}")
        return 0

    # build one engine per distinct system `a` lazily inside the loop
    engines: dict = {}
    base = wca_jobs.effective_base(cfg, args.stage) if args.stage else cfg.get("base", {})
    n_done = n_skip = n_fail = n_nan = 0
    for k, spec in enumerate(specs):
        path = wca_jobs.run_npz_path(raw_dir, spec)
        if not args.overwrite and wca_jobs.run_is_valid(path):
            n_skip += 1
            continue
        params = wca_jobs.build_params(spec)
        if spec.a not in engines:
            engines[spec.a] = core.WCADimerEngine(params, core.DEVICE, core.DTYPE)
        engine = engines[spec.a]
        try:
            out = wca_jobs.execute_run(spec, base, engine, cache_dir=cache_dir, verbose=args.verbose)
            wca_jobs.save_run(path, out)
            n_done += 1
            if bool(out["had_nan"]):
                n_nan += 1
            print(f"  [{k+1}/{len(specs)}] {spec.name} seed{spec.seed} "
                  f"N{spec.n_replicas} T{spec.n_steps} a{spec.a:g}: "
                  f"L2(F)={out['l2_f']:.4f} L2(Fp)={out['l2_fp']:.4f} "
                  f"intF={out['integrated_l2_f']:.3f} repl={int(out['total_replacement_events'])} "
                  f"essA={out['final_ancestor_ess']:.0f} "
                  f"{'NAN!' if out['had_nan'] else ''} ({out['wall_seconds']:.0f}s)")
        except Exception as exc:  # keep the sweep alive
            n_fail += 1
            print(f"  [{k+1}/{len(specs)}] {spec.run_id()} FAILED: {exc!r}")
    print(f"[run] DONE done={n_done} skipped={n_skip} failed={n_fail} nan={n_nan}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
