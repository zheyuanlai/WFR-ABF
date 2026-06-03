#!/usr/bin/env python3
"""Run the ABF / ABF+FR grid on the PyTorch (GPU) backend.

This is the GPU counterpart of ``scripts/run_abf_fr_grid.py``.  It batches
independent runs onto one device, checkpoints/resumes per run, and writes the
same CSV schema (a superset) under the stage directory so the existing plotting
and report-table scripts work with ``--stage tuning_gpu`` / ``production_gpu``.

The science is identical to the CPU reference (same potential, reaction
coordinate, ABF/FR targets and metrics); the GPU backend uses a fast
binned-smoothed 1D estimator validated against the CPU/kernel estimator by
``scripts/validate_torch_backend.py``.

Examples
--------
  # Local CPU-torch smoke (writes results/two_dim_xi_x/tuning_gpu/)
  python scripts/run_abf_fr_grid_torch.py \
      --config configs/two_dim_xi_x_smoke_gpu.yaml --stage tuning_gpu --device cpu

  # One shard of a sharded production run (usually via run_gpu_shard.py)
  python scripts/run_abf_fr_grid_torch.py \
      --config configs/two_dim_xi_x_production_gpu.yaml --stage production_gpu \
      --shard results/two_dim_xi_x/production_gpu/shards/shard_000.json

  # Merge tagged shard CSVs + write config summaries after all shards finish
  python scripts/run_abf_fr_grid_torch.py \
      --config configs/two_dim_xi_x_production_gpu.yaml --stage production_gpu \
      --merge-only

This script never launches the full production grid implicitly -- you must pass
the production config and stage explicitly.  Use --dry-run to preview.
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
    p.add_argument("--shard", default=None,
                   help="Optional shard JSON (list of run specs) to run instead "
                        "of the full grid.")
    p.add_argument("--tag", default=None,
                   help="CSV/output tag (default: derived from --shard or 'main').")
    p.add_argument("--merge-only", action="store_true",
                   help="Skip simulation: merge tagged CSVs + write config "
                        "summaries (run after all shards finish).")
    # Resume / safety.
    p.add_argument("--resume", dest="resume", action="store_true", default=True)
    p.add_argument("--no-resume", dest="resume", action="store_false")
    p.add_argument("--force", action="store_true",
                   help="Re-run even completed run_ids.")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--max-runs", type=int, default=None)
    p.add_argument("--conditional-snapshots", default="final",
                   choices=["final", "all"])
    p.add_argument("--base-seed", type=int, default=0)
    # Config overrides.
    p.add_argument("--device", default=None)
    p.add_argument("--dtype", default=None, choices=["float32", "float64"])
    p.add_argument("--batch-size-configs", type=int, default=None)
    p.add_argument("--n-steps", type=int, default=None)
    p.add_argument("--n-particles", type=int, default=None)
    p.add_argument("--eval-every", type=int, default=None)
    p.add_argument("--estimator", default=None,
                   choices=["binned_smooth", "kernel_reference"])
    p.add_argument("--seeds", default=None, help="Comma-separated seed override.")
    p.add_argument("--output-root", default=None)
    return p.parse_args(argv)


def _load_shard_specs(path):
    with open(path) as fh:
        data = json.load(fh)
    runs = data["runs"] if isinstance(data, dict) and "runs" in data else data
    return [RunSpec(method=r["method"], target_type=r["target_type"],
                    seed=int(r["seed"]), gamma=float(r["gamma"]),
                    eta=float(r["eta"]), burnin_fraction=float(r["burnin_fraction"]),
                    fr_every=int(r["fr_every"])) for r in runs]


def main(argv=None):
    args = parse_args(argv)
    cfg = io_utils.load_config(args.config)
    seeds = ([int(s) for s in args.seeds.split(",")] if args.seeds else None)
    io_utils.apply_cli_overrides(
        cfg, device=args.device, dtype=args.dtype,
        batch_size_configs=args.batch_size_configs, n_steps=args.n_steps,
        n_particles=args.n_particles, eval_every=args.eval_every, seeds=seeds,
        output_root=args.output_root, estimator=args.estimator)

    stage_root = io_utils.stage_dir(cfg, args.stage)
    prefix = io_utils.stage_prefix(args.stage)

    if args.merge_only:
        parallel.merge_stage_csvs(stage_root, prefix)
        parallel.write_config_summaries(stage_root, prefix, cfg)
        print(f"[run_abf_fr_grid_torch] merge-only done; outputs under "
              f"{os.path.relpath(stage_root)}/")
        return 0

    # Reference gate + device/dtype/eval setup.
    setup = parallel.prepare_stage(cfg, args.stage, require_csv=True)

    # Build the run list (full grid or a shard subset).
    if args.shard:
        specs = _load_shard_specs(args.shard)
        tag = args.tag or os.path.splitext(os.path.basename(args.shard))[0]
    else:
        run_seeds = seeds if seeds is not None else [int(s) for s in cfg["simulation"]["seeds"]]
        specs = io_utils.build_run_specs(cfg, run_seeds)
        tag = args.tag or "main"
    if args.max_runs is not None:
        specs = specs[:args.max_runs]

    estimator = cfg.get("abf", {}).get("estimator", "binned_smooth")
    print(f"[run_abf_fr_grid_torch] stage={args.stage} device={setup['device']} "
          f"dtype={setup['dtype']} estimator={estimator} "
          f"batch_size_configs={cfg.get('batch_size_configs', 16)} "
          f"n_runs={len(specs)} tag={tag} out={os.path.relpath(stage_root)}")
    if args.dry_run:
        for s in specs:
            print("   ", s.run_id)
        print(f"[run_abf_fr_grid_torch] dry-run: {len(specs)} runs planned "
              f"(not executed).")
        return 0

    summary = parallel.run_specs(
        specs, cfg=cfg, stage_root=stage_root, prefix=prefix,
        x_grid=setup["x_grid"], ref=setup["ref"], ev=setup["ev"],
        device=setup["device"], dtype=setup["dtype"], estimator=estimator,
        batch_size=int(cfg.get("batch_size_configs", 16)), base_seed=args.base_seed,
        tag=tag, resume=args.resume, force=args.force,
        conditional=args.conditional_snapshots)

    # Single-process full grid: merge + config summaries here.  Sharded runs do
    # this once afterwards via --merge-only (see run_full_gpu_study.sh).
    if not args.shard:
        parallel.merge_stage_csvs(stage_root, prefix)
        parallel.write_config_summaries(stage_root, prefix, cfg)

    print(f"[run_abf_fr_grid_torch] DONE  done={summary['n_done']} "
          f"failed={summary['n_failed']} skipped={summary['n_skipped']} "
          f"nan={summary['n_nan']} in {summary['wall_seconds']:.1f}s")
    if summary["n_failed"]:
        print(f"[run_abf_fr_grid_torch] WARNING: {summary['n_failed']} run(s) "
              f"failed; see {os.path.relpath(os.path.join(stage_root, 'failed'))}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
