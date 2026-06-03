#!/usr/bin/env python3
"""Split a GPU study's run list into shard files for multi-GPU execution.

Each shard is a JSON file listing the runs to execute; one shard runs on one
GPU via ``scripts/run_gpu_shard.py``.  Runs that *batch together* (same
target_type, eta, fr_every, burnin_fraction) are kept in the same shard via
greedy bin-packing so each shard's GPU batches stay full.

CLI:
  python scripts/make_gpu_shards.py \
      --config configs/two_dim_xi_x_production_gpu.yaml \
      --stage production_gpu --num-shards 8

Writes:
  results/two_dim_xi_x/production_gpu/shards/shard_000.json
  ...
and prints the estimated number of runs.
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from abffr import io_utils, parallel  # noqa: E402

STAGES = ["smoke_gpu", "tuning_gpu", "production_gpu"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True)
    p.add_argument("--stage", required=True, choices=STAGES)
    p.add_argument("--num-shards", type=int, required=True)
    p.add_argument("--seeds", default=None, help="Comma-separated seed override.")
    p.add_argument("--output-root", default=None)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    cfg = io_utils.load_config(args.config)
    if args.output_root:
        cfg["output_root"] = args.output_root
    seeds = ([int(s) for s in args.seeds.split(",")] if args.seeds
             else [int(s) for s in cfg["simulation"]["seeds"]])

    specs = io_utils.build_run_specs(cfg, seeds)
    if not specs:
        print("[make_gpu_shards] ERROR: no runs built from config.")
        return 1

    # Group by batch key, then greedy bin-pack groups into shards (largest
    # first) to balance run counts while keeping batchable runs together.
    groups = defaultdict(list)
    for s in specs:
        groups[parallel.batch_key(s)].append(s)
    group_list = sorted(groups.values(), key=len, reverse=True)

    n = max(1, int(args.num_shards))
    shards = [[] for _ in range(n)]
    load = [0] * n
    for grp in group_list:
        j = min(range(n), key=lambda i: load[i])
        shards[j].extend(grp)
        load[j] += len(grp)

    stage_root = io_utils.stage_dir(cfg, args.stage)
    shard_dir = io_utils.ensure_dir(os.path.join(stage_root, "shards"))
    # Clear any stale shards so a re-shard is clean.
    for old in os.listdir(shard_dir):
        if old.startswith("shard_") and old.endswith(".json"):
            os.remove(os.path.join(shard_dir, old))

    written = 0
    for i, shard in enumerate(shards):
        path = os.path.join(shard_dir, f"shard_{i:03d}.json")
        io_utils.save_json(path, dict(
            stage=args.stage, config=os.path.abspath(args.config),
            shard_index=i, num_shards=n, num_runs=len(shard),
            runs=[s.to_row() for s in shard]))
        written += 1
        print(f"   shard_{i:03d}.json: {len(shard)} runs")

    print(f"[make_gpu_shards] stage={args.stage} total_runs={len(specs)} "
          f"shards={written} out={os.path.relpath(shard_dir)}")
    print(f"[make_gpu_shards] ESTIMATED NUMBER OF RUNS: {len(specs)} "
          f"(methods x gamma x eta x burnin x fr_every x seeds)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
