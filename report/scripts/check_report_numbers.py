#!/usr/bin/env python3
r"""Verify that the numbers used in the report match the CSV-derived values.

This is an independent cross-check of ``build_report_assets.py``: it recomputes
the headline quantities straight from the merged production CSVs and asserts
that

  1. ``report/tables/report_numbers.json`` agrees with the recomputation;
  2. the LaTeX macros in ``report/tables/numbers.tex`` round to those values;
  3. the rounded medians appear verbatim in the main result tables
     (``best_configs_by_integrated_F.tex`` / ``best_configs_by_final_F.tex``).

Exit status is non-zero if any check fails, so it can gate ``make``.
"""
from __future__ import annotations

import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPORT_DIR = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(REPORT_DIR)
TAB_DIR = os.path.join(REPORT_DIR, "tables")
PROD = os.path.join(REPO_ROOT, "results", "two_dim_xi_x", "production_gpu")
EB_DIR = os.path.join(REPO_ROOT, "results", "entropic_bottleneck", "summaries")
WCA_DIR = os.path.join(REPO_ROOT, "results", "wca_production", "summaries")

FINAL = "final_l2_F"
INTEGRATED = "integrated_l2_F"

_fails = []
_checks = 0


def check(cond, msg):
    global _checks
    _checks += 1
    if cond:
        print(f"  [ok]   {msg}")
    else:
        print(f"  [FAIL] {msg}")
        _fails.append(msg)


def approx(a, b, rtol=1e-4, atol=1e-9):
    return bool(np.isfinite(a) and np.isfinite(b)
                and abs(a - b) <= atol + rtol * abs(b))


def load_merged(base):
    merged = os.path.join(PROD, f"{base}.csv")
    if os.path.exists(merged):
        return pd.read_csv(merged)
    shard = os.path.join(PROD, f"{base}__shard_000.csv")
    return pd.read_csv(shard)


def med(s):
    return float(np.nanmedian(np.asarray(s, dtype=float)))


def best_config(g, metric):
    return g.groupby("config_id")[metric].median().sort_values().index[0]


def prob_beats(g, abf_by_seed):
    wins = n = 0
    for _, r in g.iterrows():
        a = abf_by_seed.get(r["seed"])
        if a is None or not np.isfinite(r[FINAL]):
            continue
        n += 1
        wins += int(r[FINAL] < a)
    return wins / n if n else float("nan")


def _eb_row(df, **kw):
    m = df
    for k, v in kw.items():
        m = m[m[k] == v]
    return m.iloc[0]


def check_eb(nums, macros):
    """Recompute EB headline numbers from config_summary medians.

    Note: EB/WCA summaries are pre-aggregated medians-over-seeds (not per-run
    like the metastability study), so we cross-check against the published
    config-summary rows rather than recomputing from raw runs.
    """
    print("\n[entropic bottleneck]")
    path = os.path.join(EB_DIR, "config_summary.csv")
    if not os.path.exists(path):
        check(False, "EB config_summary.csv exists")
        return
    cs = pd.read_csv(path)
    s1_fr = _eb_row(cs, stage="stage1_seeds", method="fr_estimated")
    gain = 100 * float(s1_fr["gain_l2_f_vs_abf"])
    win = int(round(s1_fr["winrate_vs_abf"] * s1_fr["n_seeds"]))
    ebn = nums.get("entropic_bottleneck", {})
    check(approx(ebn.get("gain_strong_pct", float("nan")), gain, rtol=1e-3),
          f"json EB cold gain% matches CSV ({gain:.1f})")
    check(ebn.get("cold_win") == win, f"json EB cold wins matches CSV ({win})")
    check(macros.get("EBgainStrong") == f"{gain:.1f}",
          "macro EBgainStrong rounds to CSV")
    check(macros.get("EBcoldWin") == str(win), "macro EBcoldWin matches CSV")
    # sign-flip sanity: beta=4 harmful, beta=8 helpful
    b4 = 100 * float(_eb_row(cs, stage="stage3_beta", method="fr_estimated",
                             beta=4.0)["gain_l2_f_vs_abf"])
    check(b4 < 0 and gain > 0, f"EB gain flips sign (beta4={b4:.1f}, cold={gain:.1f})")


def check_wca(nums, macros):
    print("\n[WCA dimer]")
    cpath = os.path.join(WCA_DIR, "config_summary.csv")
    wpath = os.path.join(WCA_DIR, "winrates.csv")
    if not (os.path.exists(cpath) and os.path.exists(wpath)):
        check(False, "WCA config_summary.csv + winrates.csv exist")
        return
    wr = pd.read_csv(wpath)
    w = wr[(wr["stage"] == "main") & (wr["name"] == "fr_est_tuned")].iloc[0]
    gain = float(w["median_gain_pct_F"])
    win = int(w["n_wins"])
    wcan = nums.get("wca", {})
    check(approx(wcan.get("tuned_gain_pct", float("nan")), gain, rtol=1e-3),
          f"json WCA tuned gain% matches CSV ({gain:.1f})")
    check(wcan.get("tuned_win") == win, f"json WCA tuned wins matches CSV ({win})")
    check(macros.get("WCAtunedGainPct") == f"{gain:.1f}",
          "macro WCAtunedGainPct rounds to CSV")
    check(macros.get("WCAtunedWin") == str(win), "macro WCAtunedWin matches CSV")
    # failure boundary: fr_rate 0.5 reverses (gain < tuned, and ESS collapses)
    aggr = float(wr[(wr["stage"] == "main")
                    & (wr["name"] == "fr_est_aggressive")].iloc[0]["median_gain_pct_F"])
    check(aggr < 0 < gain, f"WCA aggressive reverses (aggr={aggr:.1f}, tuned={gain:.1f})")


def main():
    print("[check_report_numbers] recomputing from CSVs in",
          os.path.relpath(PROD, REPO_ROOT))
    fs = load_merged("production_gpu_final_summary")

    # ---- structural checks ------------------------------------------------ #
    n_runs = len(fs)
    n_configs = fs["config_id"].nunique()
    n_seeds = fs["seed"].nunique()
    n_methods = fs["method"].nunique()
    n_nan = int(fs[[FINAL, "final_l2_Fprime", INTEGRATED]].isna().to_numpy().sum())
    print("\n[structure]")
    check(n_runs == 545, f"n_runs == 545 (got {n_runs})")
    check(n_configs == 109, f"n_configs == 109 (got {n_configs})")
    check(n_seeds == 5, f"n_seeds == 5 (got {n_seeds})")
    check(n_methods == 4, f"n_methods == 4 (got {n_methods})")
    check(n_nan == 0, f"no NaNs in key metrics (got {n_nan})")

    # ---- recompute headline numbers --------------------------------------- #
    abf = fs[fs["method"] == "abf_only"]
    abf_by_seed = abf.set_index("seed")[FINAL].to_dict()
    recomputed = {}
    recomputed["abf"] = dict(
        integ=med(abf[INTEGRATED]), final=med(abf[FINAL]),
        finalp=med(abf["final_l2_Fprime"]))
    for m in ["abf_fr_estimated", "abf_fr_uniform", "abf_fr_oracle"]:
        g = fs[fs["method"] == m]
        ci, cf = best_config(g, INTEGRATED), best_config(g, FINAL)
        gi, gf = g[g["config_id"] == ci], g[g["config_id"] == cf]
        recomputed[m] = dict(
            int_integ=med(gi[INTEGRATED]), int_final=med(gi[FINAL]),
            fin_final=med(gf[FINAL]), fin_finalp=med(gf["final_l2_Fprime"]),
            fin_integ=med(gf[INTEGRATED]),
            int_prob=prob_beats(gi, abf_by_seed),
            fin_prob=prob_beats(gf, abf_by_seed))

    # ---- compare against report_numbers.json ------------------------------ #
    jpath = os.path.join(TAB_DIR, "report_numbers.json")
    print("\n[report_numbers.json]")
    if not os.path.exists(jpath):
        check(False, "report_numbers.json exists (run build_report_assets.py)")
        return _finish()
    nums = json.load(open(jpath))
    check(nums["meta"]["n_runs"] == n_runs, "json n_runs matches CSV")
    check(nums["meta"]["n_configs"] == n_configs, "json n_configs matches CSV")
    check(approx(nums["abf_only"]["median_integrated_l2_F"], recomputed["abf"]["integ"]),
          "json ABF integrated L2(F) matches CSV")
    check(approx(nums["abf_only"]["median_final_l2_F"], recomputed["abf"]["final"]),
          "json ABF final L2(F) matches CSV")
    check(approx(nums["estimated_by_integrated"]["median_integrated_l2_F"],
                 recomputed["abf_fr_estimated"]["int_integ"]),
          "json estimated(by-integrated) integrated L2(F) matches CSV")
    check(approx(nums["estimated_by_final"]["median_final_l2_F"],
                 recomputed["abf_fr_estimated"]["fin_final"]),
          "json estimated(by-final) final L2(F) matches CSV")
    check(approx(nums["estimated_by_final"]["prob_beats_abf"],
                 recomputed["abf_fr_estimated"]["fin_prob"]),
          "json estimated(by-final) prob_beats matches CSV")
    check(approx(nums["uniform_by_integrated"]["median_integrated_l2_F"],
                 recomputed["abf_fr_uniform"]["int_integ"]),
          "json uniform(by-integrated) integrated L2(F) matches CSV")
    check(approx(nums["oracle_best"]["median_integrated_l2_F"],
                 recomputed["abf_fr_oracle"]["int_integ"]),
          "json oracle integrated L2(F) matches CSV")

    # improvement percentages
    imp = nums["improvements"]
    exp_est_int = 100 * (recomputed["abf"]["integ"]
                         - recomputed["abf_fr_estimated"]["int_integ"]) / recomputed["abf"]["integ"]
    check(approx(imp["est_int_integratedF_pct"], exp_est_int),
          "json estimated integrated improvement % matches CSV")
    exp_orac_int = 100 * (recomputed["abf"]["integ"]
                          - recomputed["abf_fr_oracle"]["int_integ"]) / recomputed["abf"]["integ"]
    check(approx(imp["oracle_integratedF_pct"], exp_orac_int),
          "json oracle integrated improvement % matches CSV")

    # ---- numbers.tex macros round to the CSV values ----------------------- #
    print("\n[numbers.tex macros]")
    mpath = os.path.join(TAB_DIR, "numbers.tex")
    macros = {}
    if os.path.exists(mpath):
        for line in open(mpath):
            m = re.match(r"\\newcommand\{\\(\w+)\}\{(.*)\}", line.strip())
            if m:
                macros[m.group(1)] = m.group(2)
    check(macros.get("NRuns") == str(n_runs), "macro NRuns matches CSV")
    check(macros.get("NConfigs") == str(n_configs), "macro NConfigs matches CSV")
    check(macros.get("AbfIntegF") == f"{recomputed['abf']['integ']:.2f}",
          "macro AbfIntegF rounds to CSV")
    check(macros.get("EstIntegF") == f"{recomputed['abf_fr_estimated']['int_integ']:.2f}",
          "macro EstIntegF rounds to CSV")
    check(macros.get("EstFinFinalF") == f"{recomputed['abf_fr_estimated']['fin_final']:.4f}",
          "macro EstFinFinalF rounds to CSV")
    check(macros.get("OracIntegF") == f"{recomputed['abf_fr_oracle']['int_integ']:.2f}",
          "macro OracIntegF rounds to CSV")

    # ---- rounded medians appear in the main result tables ----------------- #
    print("\n[main tables]")
    for fname, key, valfmt in [
        ("best_configs_by_integrated_F.tex",
         f"{recomputed['abf']['integ']:.2f}", "ABF integrated"),
        ("best_configs_by_integrated_F.tex",
         f"{recomputed['abf_fr_estimated']['int_integ']:.2f}", "estimated integrated"),
        ("best_configs_by_final_F.tex",
         f"{recomputed['abf_fr_estimated']['fin_final']:.4f}", "estimated final"),
    ]:
        path = os.path.join(TAB_DIR, fname)
        txt = open(path).read() if os.path.exists(path) else ""
        check(key in txt, f"table {fname} contains {valfmt} value {key}")

    # ---- EB + WCA case numbers -------------------------------------------- #
    check_eb(nums, macros)
    check_wca(nums, macros)

    return _finish()


def _finish():
    print(f"\n[check_report_numbers] {_checks - len(_fails)}/{_checks} checks passed.")
    if _fails:
        print("FAILED:")
        for f in _fails:
            print("  -", f)
        return 1
    print("All report numbers are consistent with the CSV-derived values.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
