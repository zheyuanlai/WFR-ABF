#!/usr/bin/env python3
"""Aggregate ABF--FR results into report tables.

Usage:
    python scripts/make_report_tables.py --stage tuning
    python scripts/make_report_tables.py --stage eval

Outputs:
    tuning:  results/two_dim_xi_x/tuning/table_tuning_top_configs.csv
    eval:    results/two_dim_xi_x/eval/table_main_results.csv

`table_main_results.csv` reports, per method/target, the median and IQR (over
seeds) of the final/integrated L2 errors, the matched-seed probability of
beating ABF-only, median barrier crossings and median FR event fraction.

Missing inputs produce a clear warning instead of a crash.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from abffr import io_utils  # noqa: E402

DEFAULT_ROOT = "results/two_dim_xi_x"


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", required=True,
                   choices=["tuning", "eval", "smoke_gpu", "tuning_gpu",
                            "production_gpu"])
    p.add_argument("--output-root", default=DEFAULT_ROOT)
    p.add_argument("--config", default=None,
                   help="Optional config (only used to read output_root).")
    p.add_argument("--top-k", type=int, default=20,
                   help="Rows in the tuning top-configs table.")
    return p.parse_args(argv)


def _warn(msg):
    print(f"[make_report_tables] WARNING: {msg}")


def _iqr(s):
    s = np.asarray(s, dtype=float)
    s = s[np.isfinite(s)]
    return float(np.percentile(s, 75) - np.percentile(s, 25)) if s.size else float("nan")


def make_tuning_table(tuning_dir, top_k, prefix="tuning"):
    path = os.path.join(tuning_dir, f"{prefix}_config_summary.csv")
    if not os.path.exists(path):
        _warn(f"{os.path.relpath(path)} missing; run the tuning/smoke stage first.")
        return
    cfg = pd.read_csv(path)
    cfg = cfg.sort_values("median_integrated_l2_F").reset_index(drop=True)
    cols = ["method", "target_type", "gamma", "eta", "burnin_fraction", "fr_every",
            "n_seeds", "median_integrated_l2_F", "iqr_integrated_l2_F",
            "median_final_l2_F", "iqr_final_l2_F",
            "median_final_l2_Fprime", "iqr_final_l2_Fprime",
            "median_mean_fr_event_fraction", "median_max_fr_event_fraction",
            "median_barrier_crossings", "frac_nan"]
    cols = [c for c in cols if c in cfg.columns]
    out = cfg[cols].head(top_k).copy()
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    dest = os.path.join(tuning_dir, f"table_{prefix}_top_configs.csv")
    out.to_csv(dest, index=False)
    print(f"[make_report_tables] wrote {os.path.relpath(dest)} ({len(out)} rows)")
    print(out.to_string(index=False))


def _best_config_per_method(final_df):
    """For each method pick the config with the lowest median final L2(F)."""
    chosen = {}
    for method, g in final_df.groupby("method"):
        med = g.groupby("config_id")["final_l2_F"].median().sort_values()
        chosen[method] = med.index[0]
    return chosen


def make_eval_table(eval_dir, prefix="eval"):
    path = os.path.join(eval_dir, f"{prefix}_final_summary.csv")
    if not os.path.exists(path):
        _warn(f"{os.path.relpath(path)} missing; run the eval/GPU stage first.")
        return
    final = pd.read_csv(path)

    # ABF-only per-seed final L2 for matched-seed comparison.
    abf = final[final["method"] == "abf_only"]
    abf_by_seed = abf.set_index("seed")["final_l2_F"].to_dict() if not abf.empty else {}

    chosen = _best_config_per_method(final)
    rows = []
    for method, cid in chosen.items():
        g = final[final["config_id"] == cid]
        target = g["target_type"].iloc[0]
        # Matched-seed probability of beating ABF-only on final L2(F).
        if method == "abf_only" or not abf_by_seed:
            prob = float("nan")
        else:
            wins, n = 0, 0
            for _, r in g.iterrows():
                a = abf_by_seed.get(r["seed"])
                if a is None or not np.isfinite(r["final_l2_F"]):
                    continue
                n += 1
                wins += int(r["final_l2_F"] < a)
            prob = wins / n if n else float("nan")
        rows.append(dict(
            method=method, target_type=target,
            median_final_l2_Fprime=float(np.nanmedian(g["final_l2_Fprime"])),
            iqr_final_l2_Fprime=_iqr(g["final_l2_Fprime"]),
            median_final_l2_F=float(np.nanmedian(g["final_l2_F"])),
            iqr_final_l2_F=_iqr(g["final_l2_F"]),
            median_integrated_l2_F=float(np.nanmedian(g["integrated_l2_F"])),
            iqr_integrated_l2_F=_iqr(g["integrated_l2_F"]),
            prob_beats_abf=prob,
            median_barrier_crossings=float(np.nanmedian(g["barrier_crossings"])),
            median_fr_event_fraction=float(np.nanmedian(g["mean_fr_event_fraction"])),
        ))
    out = pd.DataFrame(rows)
    # Order: ABF-only first, then FR methods.
    order = {"abf_only": 0, "abf_fr_estimated": 1, "abf_fr_uniform": 2,
             "abf_fr_oracle": 3, "abf_fr_self": 4}
    out["_o"] = out["method"].map(order).fillna(9)
    out = out.sort_values("_o").drop(columns="_o").reset_index(drop=True)
    dest = os.path.join(eval_dir, "table_main_results.csv")
    out.to_csv(dest, index=False)
    print(f"[make_report_tables] wrote {os.path.relpath(dest)} ({len(out)} rows)")
    print(out.to_string(index=False))


def main(argv=None):
    args = parse_args(argv)
    root = args.output_root
    if args.config:
        cfg = io_utils.load_config(args.config)
        root = cfg.get("output_root", root)
    sub = io_utils.STAGE_TO_DIR.get(args.stage, args.stage)
    prefix = io_utils.stage_prefix(args.stage)
    data_dir = os.path.join(root, sub)
    if args.stage == "tuning":
        make_tuning_table(data_dir, args.top_k, prefix=prefix)
    elif args.stage == "eval":
        make_eval_table(data_dir, prefix=prefix)
    else:  # GPU stages: full grid -> both the top-configs and main-results tables
        make_tuning_table(data_dir, args.top_k, prefix=prefix)
        make_eval_table(data_dir, prefix=prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
