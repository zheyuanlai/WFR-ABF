#!/usr/bin/env python3
"""Aggregate WCA production raw runs into summary CSVs / NPZ.

Reads results/wca_production/raw/*.npz and writes, under
results/wca_production/summaries/ (override with --out):
  summary.csv               one row per run (scalar metrics + metadata)
  config_summary.csv        per (stage,name) median/IQR over seeds
  timeseries_summary.csv    median/IQR of L2(F),L2(Fp) over seeds, per time
  profiles_summary.npz      seed-mean final F(z),F'(z),p(z),q(z),N_eff,birth/death
  winrates.csv              matched-seed win rate of each FR config vs abf

Usage:
  python scripts/analyze_wca_production.py
  python scripts/analyze_wca_production.py --stages main pilot
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import numpy as np  # noqa: E402
import wca_jobs  # noqa: E402

SCALAR_COLS = [
    "run_id", "stage", "name", "method", "seed", "n_steps", "n_replicas", "a",
    "fr_rate", "target_ema_rate", "max_event_fraction", "fr_every", "fr_start_steps",
    "score_clip", "config_hash", "core_version", "had_nan", "total_replacement_events",
    "l2_f", "l2_fp", "integrated_l2_f",
    "l2_f_compact", "l2_f_transition", "l2_f_stretched",
    "l2_fp_compact", "l2_fp_transition", "l2_fp_stretched",
    "final_ancestor_ess", "final_n_unique_ancestor", "runtime_seconds", "wall_seconds",
]


def _val(d, k):
    v = d[k]
    if isinstance(v, np.ndarray) and v.ndim == 0:
        v = v.item()
    return v


def load_summary_rows(raw_dir, stages=None):
    rows = []
    for path in sorted(glob.glob(os.path.join(raw_dir, "*.npz"))):
        try:
            d = wca_jobs.load_run(path)
        except Exception as exc:
            print(f"  skip unreadable {path}: {exc!r}")
            continue
        if "l2_f" not in d:
            continue
        if stages and str(_val(d, "stage")) not in stages:
            continue
        rows.append({k: _val(d, k) for k in SCALAR_COLS if k in d})
    return rows


def write_summary_csv(rows, out_path):
    if not rows:
        print("  no rows -> no summary.csv")
        return
    cols = [c for c in SCALAR_COLS if c in rows[0]]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  wrote {out_path} ({len(rows)} runs)")


def _iqr(a):
    a = np.asarray(a, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan"), float("nan"), float("nan")
    return float(np.median(a)), float(np.percentile(a, 25)), float(np.percentile(a, 75))


def write_config_summary(rows, out_path):
    """Per (stage,name,n_steps,n_replicas,a) median/IQR over seeds."""
    groups = {}
    for r in rows:
        key = (r["stage"], r["name"], r["method"], int(r["n_steps"]),
               int(r["n_replicas"]), float(r["a"]))
        groups.setdefault(key, []).append(r)
    out = []
    for key, rs in sorted(groups.items()):
        rec = dict(stage=key[0], name=key[1], method=key[2], n_steps=key[3],
                   n_replicas=key[4], a=key[5], n_seeds=len(rs))
        for metric in ["l2_f", "l2_fp", "integrated_l2_f", "l2_f_transition",
                       "l2_fp_transition", "final_ancestor_ess", "total_replacement_events"]:
            med, lo, hi = _iqr([r.get(metric, float("nan")) for r in rs])
            rec[f"{metric}_median"] = med
            rec[f"{metric}_q25"] = lo
            rec[f"{metric}_q75"] = hi
        out.append(rec)
    if out:
        cols = list(out[0].keys())
        with open(out_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(out)
        print(f"  wrote {out_path} ({len(out)} configs)")
    return out


def write_winrates(rows, out_path):
    """Matched-seed win rate of each FR config vs the abf baseline.

    Matched within (stage, n_steps, n_replicas, a, seed). Win = lower final L2(F).
    Also reports median paired % gain = 100*(abf - fr)/abf.
    """
    # index abf baselines by matched key+seed
    def mkey(r):
        return (r["stage"], int(r["n_steps"]), int(r["n_replicas"]), float(r["a"]), int(r["seed"]))
    abf = {}
    for r in rows:
        if r["method"] == "abf":
            abf[mkey(r)] = r
    groups = {}
    for r in rows:
        if r["method"] == "abf":
            continue
        gk = (r["stage"], r["name"], r["method"], int(r["n_steps"]), int(r["n_replicas"]), float(r["a"]))
        base = abf.get(mkey(r))
        if base is None:
            continue
        groups.setdefault(gk, []).append((r, base))
    out = []
    for gk, pairs in sorted(groups.items()):
        gains_f, gains_int, wins = [], [], 0
        for fr, base in pairs:
            if base["l2_f"] > 0:
                gains_f.append(100.0 * (base["l2_f"] - fr["l2_f"]) / base["l2_f"])
            if base.get("integrated_l2_f", 0) and base["integrated_l2_f"] > 0:
                gains_int.append(100.0 * (base["integrated_l2_f"] - fr["integrated_l2_f"]) / base["integrated_l2_f"])
            wins += int(fr["l2_f"] < base["l2_f"])
        n = len(pairs)
        rec = dict(stage=gk[0], name=gk[1], method=gk[2], n_steps=gk[3], n_replicas=gk[4], a=gk[5],
                   n_pairs=n, n_wins=wins, win_rate=(wins / n if n else float("nan")),
                   median_gain_pct_F=float(np.median(gains_f)) if gains_f else float("nan"),
                   mean_gain_pct_F=float(np.mean(gains_f)) if gains_f else float("nan"),
                   median_gain_pct_intF=float(np.median(gains_int)) if gains_int else float("nan"))
        out.append(rec)
    if out:
        cols = list(out[0].keys())
        with open(out_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(out)
        print(f"  wrote {out_path} ({len(out)} FR configs)")
    return out


def write_timeseries_summary(raw_dir, out_path, stages=None):
    """Median/IQR of L2(F),L2(Fp) over seeds, per (stage,name,...) and saved time."""
    series = {}
    for path in sorted(glob.glob(os.path.join(raw_dir, "*.npz"))):
        d = wca_jobs.load_run(path)
        if "times" not in d:
            continue
        if stages and str(_val(d, "stage")) not in stages:
            continue
        key = (str(_val(d, "stage")), str(_val(d, "name")), str(_val(d, "method")),
               int(_val(d, "n_steps")), int(_val(d, "n_replicas")), float(_val(d, "a")))
        series.setdefault(key, []).append(
            (np.asarray(d["times"], float), np.asarray(d["l2_f_t"], float),
             np.asarray(d["l2_fp_t"], float), np.asarray(d["repl_cumulative"], float)))
    rows = []
    for key, runs in sorted(series.items()):
        L = min(len(t) for t, _, _, _ in runs)
        times = runs[0][0][:L]
        F = np.stack([r[1][:L] for r in runs])
        Fp = np.stack([r[2][:L] for r in runs])
        repl = np.stack([r[3][:L] for r in runs])
        for j in range(L):
            rows.append(dict(
                stage=key[0], name=key[1], method=key[2], n_steps=key[3], n_replicas=key[4], a=key[5],
                t=float(times[j]),
                l2_f_median=float(np.median(F[:, j])), l2_f_q25=float(np.percentile(F[:, j], 25)),
                l2_f_q75=float(np.percentile(F[:, j], 75)),
                l2_fp_median=float(np.median(Fp[:, j])), l2_fp_q25=float(np.percentile(Fp[:, j], 25)),
                l2_fp_q75=float(np.percentile(Fp[:, j], 75)),
                repl_cumulative_median=float(np.median(repl[:, j])), n_seeds=len(runs)))
    if rows:
        cols = list(rows[0].keys())
        with open(out_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        print(f"  wrote {out_path} ({len(series)} configs x times)")


def write_profiles_summary(raw_dir, out_path, stages=None):
    """Seed-mean final profiles per (stage,name,...): F,F',p,q,N_eff,birth,death."""
    groups = {}
    grid = ref_F = ref_Fp = edges = None
    for path in sorted(glob.glob(os.path.join(raw_dir, "*.npz"))):
        d = wca_jobs.load_run(path)
        if "final_pmf" not in d:
            continue
        if stages and str(_val(d, "stage")) not in stages:
            continue
        key = "|".join([str(_val(d, "stage")), str(_val(d, "name")), str(_val(d, "method")),
                        str(int(_val(d, "n_steps"))), str(int(_val(d, "n_replicas"))), f"{float(_val(d,'a')):g}"])
        groups.setdefault(key, []).append(d)
        if grid is None:
            grid, ref_F, ref_Fp = np.asarray(d["grid"]), np.asarray(d["ref_free_energy"]), np.asarray(d["ref_mean_force"])
            edges = np.asarray(d["hist_edges"])
    out = {"grid": grid, "ref_free_energy": ref_F, "ref_mean_force": ref_Fp, "hist_edges": edges}
    for key, ds in groups.items():
        def stk(field):
            return np.stack([np.asarray(x[field], float) for x in ds])
        import warnings
        with warnings.catch_warnings():  # abf has all-NaN q_target; mean is NaN by design
            warnings.simplefilter("ignore", category=RuntimeWarning)
            out[f"{key}@F"] = np.nanmean(stk("final_pmf"), axis=0)
            out[f"{key}@Fp"] = np.nanmean(stk("final_mean_force"), axis=0)
            out[f"{key}@p"] = np.nanmean(stk("final_p_hat"), axis=0)
            out[f"{key}@q"] = np.nanmean(stk("final_q_target"), axis=0)
            out[f"{key}@neff"] = np.nanmean(stk("final_eff_counts"), axis=0)
            out[f"{key}@birth"] = np.nansum(stk("birth_hist"), axis=0)
            out[f"{key}@death"] = np.nansum(stk("death_hist"), axis=0)
    np.savez_compressed(out_path, **out)
    print(f"  wrote {out_path} ({len(groups)} configs)")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/wca_production.yaml")
    ap.add_argument("--raw", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--stages", nargs="*", default=None)
    args = ap.parse_args(argv)
    cfg = wca_jobs.load_yaml(args.config)
    raw_dir = args.raw or os.path.join(cfg["output_root"], "raw")
    out_dir = args.out or os.path.join(cfg["output_root"], "summaries")
    os.makedirs(out_dir, exist_ok=True)
    print(f"[analyze] raw={raw_dir} out={out_dir} stages={args.stages or 'ALL'}")

    rows = load_summary_rows(raw_dir, args.stages)
    write_summary_csv(rows, os.path.join(out_dir, "summary.csv"))
    write_config_summary(rows, os.path.join(out_dir, "config_summary.csv"))
    write_winrates(rows, os.path.join(out_dir, "winrates.csv"))
    write_timeseries_summary(raw_dir, os.path.join(out_dir, "timeseries_summary.csv"), args.stages)
    write_profiles_summary(raw_dir, os.path.join(out_dir, "profiles_summary.npz"), args.stages)
    print("[analyze] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
