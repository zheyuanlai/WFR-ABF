#!/usr/bin/env python3
"""Tiny smoke test for the WCA production pipeline.

Runs the `smoke` stage (few replicas, few steps) for all four methods, verifies:
  * every method executes with no NaNs in the final profiles;
  * fr_estimated runs WITHOUT any oracle/TI input (no-leakage guard holds);
  * the no-leakage guard actually fires when an oracle is wrongly passed;
  * per-run .npz output files are created and reload with required keys.

This is a sanity check only -- NOT a scientific result (2k steps << 10k warmup,
so the gains here are meaningless). Run on any free GPU:

  CUDA_VISIBLE_DEVICES=7 python scripts/smoke_wca_production.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import numpy as np  # noqa: E402
import wca_abffr_core as core  # noqa: E402
import wca_jobs  # noqa: E402

CONFIG = "configs/wca_production.yaml"
REQUIRED_KEYS = ["l2_f", "l2_fp", "times", "l2_f_t", "final_pmf", "final_mean_force",
                 "final_p_hat", "birth_hist", "death_hist", "config_hash", "spec_json"]


def main():
    cfg = wca_jobs.load_yaml(CONFIG)
    raw_dir = os.path.join(cfg["output_root"], "raw")
    cache_dir = cfg.get("cache_dir", "cache")
    base = wca_jobs.effective_base(cfg, "smoke")
    specs = wca_jobs.expand_stage(cfg, "smoke")
    print(f"[smoke] device={core.DEVICE} jobs={len(specs)}")

    engine = core.WCADimerEngine(wca_jobs.build_params(specs[0]), core.DEVICE, core.DTYPE)
    failures = []

    # 1. the no-leakage guard must FIRE for fr_estimated + oracle.
    try:
        core.assert_no_oracle_leakage("fr_estimated", np.zeros(4))
        failures.append("no-leakage guard did NOT fire for fr_estimated+oracle")
    except AssertionError:
        print("[smoke] PASS: no-leakage guard fires for fr_estimated+oracle")

    # 2. run every method; check files + finiteness.
    for spec in specs:
        out = wca_jobs.execute_run(spec, base, engine, cache_dir=cache_dir, verbose=True)
        path = wca_jobs.run_npz_path(raw_dir, spec)
        wca_jobs.save_run(path, out)
        if not os.path.exists(path):
            failures.append(f"{spec.name}: output file missing")
            continue
        if not wca_jobs.run_is_valid(path):
            failures.append(f"{spec.name}: run_is_valid False")
        reloaded = wca_jobs.load_run(path)
        missing = [k for k in REQUIRED_KEYS if k not in reloaded]
        if missing:
            failures.append(f"{spec.name}: missing keys {missing}")
        if bool(out["had_nan"]):
            failures.append(f"{spec.name}: had_nan True")
        if not np.isfinite(out["l2_f"]):
            failures.append(f"{spec.name}: l2_f not finite")
        print(f"  {spec.method:13s} L2(F)={out['l2_f']:.4f} L2(Fp)={out['l2_fp']:.4f} "
              f"repl={int(out['total_replacement_events'])} "
              f"essA={out['final_ancestor_ess']:.0f} had_nan={bool(out['had_nan'])} "
              f"file_ok={os.path.exists(path)}")

    # 3. fr_estimated must have produced a real estimated target (no NaNs in q).
    est = [s for s in specs if s.method == "fr_estimated"]
    if est:
        r = wca_jobs.load_run(wca_jobs.run_npz_path(raw_dir, est[0]))
        if not np.isfinite(np.asarray(r["final_q_target"], dtype=float)).any():
            failures.append("fr_estimated: q_target all-NaN (target never built)")
        else:
            print("[smoke] PASS: fr_estimated built a finite estimated target q_n(z)")

    print("=" * 60)
    if failures:
        print("[smoke] FAILED:")
        for f in failures:
            print("  -", f)
        return 1
    print("[smoke] ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
