#!/usr/bin/env python3
"""Aggregate entropic-bottleneck raw runs into summary CSVs + NPZ.

Reads results/entropic_bottleneck/raw/<stage>/*.npz and writes, under
results/entropic_bottleneck/summaries/:
  summary.csv          one row per run (scalar metrics + config metadata)
  config_summary.csv   per (stage, config, method) median/IQR over seeds + win rate
  arrays.npz           seed-stacked time series + final profiles + conditional
                       diagnostics, keyed by (stage|config|method), for plotting

Usage:
  python scripts/analyze_entropic_bottleneck.py
  python scripts/analyze_entropic_bottleneck.py --stages stage2_omega stage4_gamma
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

import numpy as np

RAW_ROOT = os.path.join(os.path.dirname(__file__), "..", "results", "entropic_bottleneck", "raw")
SUM_ROOT = os.path.join(os.path.dirname(__file__), "..", "results", "entropic_bottleneck", "summaries")

CFG_COLS = ["beta", "H", "omega_out", "omega_in", "s", "gamma", "N", "n_steps",
            "fr_every", "target_ema_rate", "score_clip", "max_event_fraction",
            "ess_window_steps"]
SCALAR_COLS = (["stage", "method", "target_mode", "seed"] + CFG_COLS +
               ["final_l2_f", "final_l2_fp", "int_l2_f", "final_ess",
                "n_die", "n_clone", "repl_fraction", "n_fr_apply",
                "cond_abserr_mean", "cond_abserr_max"])


def _v(d, k):
    x = d[k]
    if isinstance(x, np.ndarray) and x.ndim == 0:
        x = x.item()
    return x


def _uses_fr(method):
    return str(method).startswith("fr_")


def load_run(path, stage):
    with np.load(path, allow_pickle=True) as d:
        rec = {"stage": stage,
               "method": str(_v(d, "method")),
               "target_mode": str(_v(d, "target_mode")),
               "seed": int(_v(d, "seed"))}
        for c in CFG_COLS:
            rec[c] = _v(d, f"cfg__{c}")
        for k in ("final_l2_f", "final_l2_fp", "int_l2_f", "final_ess",
                  "n_die", "n_clone", "repl_fraction", "n_fr_apply"):
            rec[k] = float(_v(d, k))
        if not _uses_fr(rec["method"]):
            rec["n_die"] = 0.0
            rec["n_clone"] = 0.0
            rec["repl_fraction"] = 0.0
            rec["n_fr_apply"] = 0.0
        ce = d["cond_abs_err"]
        rec["cond_abserr_mean"] = float(np.nanmean(ce))
        rec["cond_abserr_max"] = float(np.nanmax(ce))
        # arrays kept for plotting
        rec["_t"] = d["t"]
        rec["_l2_f_t"] = d["l2_f_t"]
        rec["_l2_fp_t"] = d["l2_fp_t"]
        rec["_ess_t"] = d["ess_t"]
        rec["_x_grid"] = d["x_grid"]
        rec["_F_hat"] = d["F_hat"]
        rec["_Fp_hat"] = d["Fp_hat"]
        rec["_F_ref"] = d["F_ref"]
        rec["_Fp_ref"] = d["Fp_ref"]
        rec["_p_hat"] = d["p_hat"]
        rec["_cond_centers"] = d["cond_centers"]
        rec["_cond_emp_var"] = d["cond_emp_var"]
        rec["_cond_ref_var"] = d["cond_ref_var"]
        rec["_cond_abs_err"] = d["cond_abs_err"]
        rec["_cond_count"] = d["cond_count"]
    return rec


def config_key(rec):
    """Identify a config cell (everything but seed/method)."""
    return (rec["stage"], rec["beta"], rec["H"], rec["omega_out"], rec["omega_in"],
            rec["s"], rec["gamma"], rec["N"], rec["n_steps"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stages", nargs="*", default=None,
                    help="stage dir names; default = all under raw/")
    args = ap.parse_args()
    os.makedirs(SUM_ROOT, exist_ok=True)

    stage_dirs = args.stages
    if stage_dirs is None:
        stage_dirs = sorted(os.path.basename(p) for p in glob.glob(os.path.join(RAW_ROOT, "*"))
                            if os.path.isdir(p))
    print("stages:", stage_dirs)

    runs = []
    for stage in stage_dirs:
        for path in sorted(glob.glob(os.path.join(RAW_ROOT, stage, "*.npz"))):
            runs.append(load_run(path, stage))
    print(f"loaded {len(runs)} runs")

    # ---- summary.csv (per run scalars) ----
    with open(os.path.join(SUM_ROOT, "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=SCALAR_COLS)
        w.writeheader()
        for r in runs:
            w.writerow({k: r[k] for k in SCALAR_COLS})
    print(f"wrote summary.csv ({len(runs)} rows)")

    write_config_summary(runs)
    write_arrays(runs)


def write_config_summary(runs):
    """Per (stage, config, method): median/IQR over seeds + matched-seed win rate
    of each FR method vs abf (computed on shared seeds within the same config)."""
    # group by config cell, then method
    cells = {}
    for r in runs:
        cells.setdefault(config_key(r), {}).setdefault(r["method"], {})[r["seed"]] = r

    cols = ["stage", "beta", "H", "omega_out", "omega_in", "s", "gamma", "N",
            "n_steps", "method", "target_mode", "n_seeds",
            "med_l2_f", "iqr_l2_f", "med_l2_fp", "iqr_l2_fp",
            "med_int_l2_f", "med_final_ess", "med_repl_fraction",
            "med_n_die", "med_n_clone", "med_cond_abserr_mean",
            "gain_l2_f_vs_abf", "winrate_vs_abf"]

    def med_iqr(vals):
        a = np.asarray(vals, float)
        q1, q2, q3 = np.percentile(a, [25, 50, 75])
        return q2, q3 - q1

    rows = []
    for ck, by_method in sorted(cells.items()):
        abf_by_seed = by_method.get("abf", {})
        abf_med_l2f = (np.median([rr["final_l2_f"] for rr in abf_by_seed.values()])
                       if abf_by_seed else float("nan"))
        for method, by_seed in sorted(by_method.items()):
            recs = list(by_seed.values())
            ex = recs[0]
            ml2f, il2f = med_iqr([r["final_l2_f"] for r in recs])
            ml2fp, il2fp = med_iqr([r["final_l2_fp"] for r in recs])
            # matched-seed gain + win rate vs abf
            shared = sorted(set(by_seed) & set(abf_by_seed))
            if method != "abf" and shared:
                wins = np.mean([by_seed[s]["final_l2_f"] < abf_by_seed[s]["final_l2_f"]
                                for s in shared])
                gain = (abf_med_l2f - ml2f) / abf_med_l2f if abf_med_l2f == abf_med_l2f else float("nan")
            else:
                wins, gain = float("nan"), (0.0 if method == "abf" else float("nan"))
            rows.append({
                "stage": ex["stage"], "beta": ex["beta"], "H": ex["H"],
                "omega_out": ex["omega_out"], "omega_in": ex["omega_in"], "s": ex["s"],
                "gamma": ex["gamma"], "N": ex["N"], "n_steps": ex["n_steps"],
                "method": method, "target_mode": ex["target_mode"], "n_seeds": len(recs),
                "med_l2_f": ml2f, "iqr_l2_f": il2f, "med_l2_fp": ml2fp, "iqr_l2_fp": il2fp,
                "med_int_l2_f": np.median([r["int_l2_f"] for r in recs]),
                "med_final_ess": np.median([r["final_ess"] for r in recs]),
                "med_repl_fraction": np.median([r["repl_fraction"] for r in recs]),
                "med_n_die": np.median([r["n_die"] for r in recs]),
                "med_n_clone": np.median([r["n_clone"] for r in recs]),
                "med_cond_abserr_mean": np.median([r["cond_abserr_mean"] for r in recs]),
                "gain_l2_f_vs_abf": gain, "winrate_vs_abf": wins,
            })
    with open(os.path.join(SUM_ROOT, "config_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote config_summary.csv ({len(rows)} rows)")


def write_arrays(runs):
    """Seed-stacked time series + final profiles + conditional diagnostics, keyed
    'stage|method|beta|omega_in|gamma'.  Used by the plotting script."""
    groups = {}
    for r in runs:
        key = (f"{r['stage']}|{r['method']}|beta{r['beta']:g}"
               f"|oin{r['omega_in']:g}|gamma{r['gamma']:g}")
        groups.setdefault(key, []).append(r)

    out = {}
    for key, recs in groups.items():
        recs = sorted(recs, key=lambda r: r["seed"])
        out[f"{key}::t"] = recs[0]["_t"]
        out[f"{key}::x_grid"] = recs[0]["_x_grid"]
        out[f"{key}::F_ref"] = recs[0]["_F_ref"]
        out[f"{key}::Fp_ref"] = recs[0]["_Fp_ref"]
        out[f"{key}::cond_centers"] = recs[0]["_cond_centers"]
        out[f"{key}::cond_ref_var"] = recs[0]["_cond_ref_var"]
        for arr in ("l2_f_t", "l2_fp_t", "ess_t", "F_hat", "Fp_hat", "p_hat",
                    "cond_emp_var", "cond_abs_err", "cond_count"):
            out[f"{key}::{arr}"] = np.stack([r[f"_{arr}"] for r in recs], axis=0)
        out[f"{key}::seeds"] = np.array([r["seed"] for r in recs])
    np.savez_compressed(os.path.join(SUM_ROOT, "arrays.npz"), **out)
    print(f"wrote arrays.npz ({len(groups)} groups)")


if __name__ == "__main__":
    main()
