#!/usr/bin/env python3
"""Run one shard of a GPU study on one GPU (checkpoint/resume).

A shard is a JSON file produced by ``scripts/make_gpu_shards.py``.  This script
runs only that shard's runs, writing tag-suffixed CSVs
(``<prefix>_<kind>__shard_NNN.csv``) and per-run completed/failed markers.  It
does NOT merge or build config summaries -- that happens once, after all shards
finish, via ``run_abf_fr_grid_torch.py --merge-only`` (see
``scripts/run_full_gpu_study.sh``).

The device is whatever torch sees; set ``CUDA_VISIBLE_DEVICES`` to pin a GPU:

  CUDA_VISIBLE_DEVICES=0 python scripts/run_gpu_shard.py \
      --config configs/two_dim_xi_x_production_gpu.yaml \
      --stage production_gpu \
      --shard results/two_dim_xi_x/production_gpu/shards/shard_000.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from abffr import io_utils, parallel  # noqa: E402
from abffr.io_utils import RunSpec  # noqa: E402

STAGES = ["smoke_gpu", "tuning_gpu", "production_gpu"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", required=True)
    p.add_argument("--stage", required=True, choices=STAGES)
    p.add_argument("--shard", required=True, help="Path to a shard JSON file.")
    p.add_argument("--device", default=None)
    p.add_argument("--dtype", default=None, choices=["float32", "float64"])
    p.add_argument("--batch-size-configs", type=int, default=None)
    p.add_argument("--n-steps", type=int, default=None)
    p.add_argument("--n-particles", type=int, default=None)
    p.add_argument("--base-seed", type=int, default=0)
    p.add_argument("--force", action="store_true")
    p.add_argument("--no-resume", dest="resume", action="store_false", default=True)
    p.add_argument("--conditional-snapshots", default="final",
                   choices=["final", "all"])
    p.add_argument("--output-root", default=None)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    cfg = io_utils.load_config(args.config)
    io_utils.apply_cli_overrides(
        cfg, device=args.device, dtype=args.dtype,
        batch_size_configs=args.batch_size_configs, n_steps=args.n_steps,
        n_particles=args.n_particles, output_root=args.output_root)

    with open(args.shard) as fh:
        data = json.load(fh)
    runs = data.get("runs", data)
    specs = [RunSpec(method=r["method"], target_type=r["target_type"],
                     seed=int(r["seed"]), gamma=float(r["gamma"]),
                     eta=float(r["eta"]), burnin_fraction=float(r["burnin_fraction"]),
                     fr_every=int(r["fr_every"])) for r in runs]
    tag = os.path.splitext(os.path.basename(args.shard))[0]

    setup = parallel.prepare_stage(cfg, args.stage, require_csv=True)
    stage_root = io_utils.stage_dir(cfg, args.stage)
    prefix = io_utils.stage_prefix(args.stage)
    print(f"[run_gpu_shard] shard={tag} device={setup['device']} "
          f"runs={len(specs)} (CUDA_VISIBLE_DEVICES="
          f"{os.environ.get('CUDA_VISIBLE_DEVICES', 'unset')})")

    summary = parallel.run_specs(
        specs, cfg=cfg, stage_root=stage_root, prefix=prefix,
        x_grid=setup["x_grid"], ref=setup["ref"], ev=setup["ev"],
        device=setup["device"], dtype=setup["dtype"],
        estimator=cfg.get("abf", {}).get("estimator", "binned_smooth"),
        batch_size=int(cfg.get("batch_size_configs", 16)), base_seed=args.base_seed,
        tag=tag, resume=args.resume, force=args.force,
        conditional=args.conditional_snapshots)

    print(f"[run_gpu_shard] shard={tag} DONE done={summary['n_done']} "
          f"failed={summary['n_failed']} skipped={summary['n_skipped']} "
          f"nan={summary['n_nan']} in {summary['wall_seconds']:.1f}s")
    return 1 if summary["n_failed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
