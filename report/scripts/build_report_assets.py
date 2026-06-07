#!/usr/bin/env python3
r"""Build all report assets from the merged production-GPU CSVs.

This script is **lightweight**: it never re-runs a simulation.  It loads the
already-computed result CSVs (and the reference grid / profile), verifies the
study metadata (number of runs / seeds / methods / configs and absence of
NaNs), and emits

  * figures   -> report/figures/fig0?_*.png
  * tables    -> report/tables/*.csv  and  report/tables/*.tex   (booktabs)
  * numbers   -> report/tables/report_numbers.json  (machine-checkable)
                 report/tables/numbers.tex          (LaTeX \newcommand macros)

The two selection rules requested in the brief are kept strictly separate:

  best_configs_by_integrated_F  : per method, the config minimising the
                                  *median-over-seeds* integrated L2(F)
                                  (rewards fast / anytime convergence);
  best_configs_by_final_F       : per method, the config minimising the
                                  *median-over-seeds* final L2(F)
                                  (rewards final-budget accuracy).

The main-text table uses the integrated rule; the final-error rule is the
supplement.  ``check_report_numbers.py`` re-derives every number in
report_numbers.json straight from the CSVs and fails on any mismatch.

Usage
-----
    python report/scripts/build_report_assets.py
    python report/scripts/build_report_assets.py --results-root results/two_dim_xi_x
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import report_cases  # noqa: E402  (EB + WCA case assets)

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
HERE = os.path.dirname(os.path.abspath(__file__))
REPORT_DIR = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(REPORT_DIR)
FIG_DIR = os.path.join(REPORT_DIR, "figures")
TAB_DIR = os.path.join(REPORT_DIR, "tables")

# --------------------------------------------------------------------------- #
# Method styling (kept identical to src/abffr/plotting.py so the report figures
# match the in-repo figures; inlined so this script has no torch dependency).
# --------------------------------------------------------------------------- #
METHOD_ORDER = ["abf_only", "abf_fr_estimated", "abf_fr_uniform", "abf_fr_oracle"]
METHOD_COLOR = {
    "abf_only": "#222222",
    "abf_fr_estimated": "#1f77b4",
    "abf_fr_uniform": "#ff7f0e",
    "abf_fr_oracle": "#2ca02c",
}
METHOD_LABEL = {
    "abf_only": "ABF only",
    "abf_fr_estimated": "ABF+FR (estimated)",
    "abf_fr_uniform": "ABF+FR (uniform)",
    "abf_fr_oracle": "ABF+FR (oracle, diagnostic)",
}
METHOD_LABEL_SHORT = {
    "abf_only": "ABF only",
    "abf_fr_estimated": "FR estimated",
    "abf_fr_uniform": "FR uniform",
    "abf_fr_oracle": "FR oracle*",
}
INTEGRATED = "integrated_l2_F"
FINAL = "final_l2_F"
EVAL_LO, EVAL_HI = -2.5, 2.5  # interior evaluation window (metrics.py)


def _set_style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 160, "font.size": 11,
        "axes.titlesize": 12, "axes.labelsize": 11, "legend.fontsize": 9,
        "axes.grid": True, "grid.alpha": 0.25, "lines.linewidth": 1.8,
    })


def _save(fig, name):
    fig.tight_layout()
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"   [fig] wrote figures/{name}")
    return path


# --------------------------------------------------------------------------- #
# IO helpers
# --------------------------------------------------------------------------- #
def _warn(msg):
    print(f"[build_report_assets] WARNING: {msg}")


def load_merged(prod_dir, base):
    """Prefer the merged CSV; fall back to the single shard if it is missing."""
    merged = os.path.join(prod_dir, f"{base}.csv")
    if os.path.exists(merged):
        return pd.read_csv(merged), os.path.relpath(merged, REPO_ROOT)
    shard = os.path.join(prod_dir, f"{base}__shard_000.csv")
    if os.path.exists(shard):
        _warn(f"{base}.csv missing; using shard {os.path.basename(shard)}.")
        return pd.read_csv(shard), os.path.relpath(shard, REPO_ROOT)
    raise FileNotFoundError(f"Neither {base}.csv nor a shard exists in {prod_dir}")


# --------------------------------------------------------------------------- #
# Numeric helpers
# --------------------------------------------------------------------------- #
def med(s):
    return float(np.nanmedian(np.asarray(s, dtype=float)))


def iqr(s):
    s = np.asarray(s, dtype=float)
    s = s[np.isfinite(s)]
    return float(np.percentile(s, 75) - np.percentile(s, 25)) if s.size else float("nan")


def best_config(method_df, metric):
    """config_id minimising the median-over-seeds ``metric`` for one method."""
    g = method_df.groupby("config_id")[metric].median().sort_values()
    return g.index[0]


def matched_seed_prob_beats(cfg_rows, abf_by_seed, col=FINAL):
    """Matched-seed P(method beats ABF-only) on ``col`` (lower is better)."""
    wins = n = 0
    for _, r in cfg_rows.iterrows():
        a = abf_by_seed.get(r["seed"])
        if a is None or not np.isfinite(r[col]):
            continue
        n += 1
        wins += int(r[col] < a)
    return (wins / n) if n else float("nan")


def config_record(fs, cid, abf_by_seed):
    """Median/IQR summary for one config_id, plus matched-seed win rate."""
    g = fs[fs["config_id"] == cid]
    r0 = g.iloc[0]
    # ABF-only is the baseline; a self vs self win-rate is meaningless ("--").
    prob = (float("nan") if r0["method"] == "abf_only"
            else matched_seed_prob_beats(g, abf_by_seed, FINAL))
    return dict(
        config_id=cid,
        method=r0["method"], target_type=r0["target_type"],
        gamma=float(r0["gamma"]), eta=float(r0["eta"]),
        burnin_fraction=float(r0["burnin_fraction"]), fr_every=int(r0["fr_every"]),
        n_seeds=int(g["seed"].nunique()),
        median_integrated_l2_F=med(g["integrated_l2_F"]),
        iqr_integrated_l2_F=iqr(g["integrated_l2_F"]),
        median_final_l2_F=med(g["final_l2_F"]),
        iqr_final_l2_F=iqr(g["final_l2_F"]),
        median_final_l2_Fprime=med(g["final_l2_Fprime"]),
        iqr_final_l2_Fprime=iqr(g["final_l2_Fprime"]),
        median_barrier_crossings=med(g["barrier_crossings"]),
        median_mean_fr_event_fraction=med(g["mean_fr_event_fraction"]),
        prob_beats_abf=prob,
    )


# --------------------------------------------------------------------------- #
# LaTeX formatting
# --------------------------------------------------------------------------- #
def f_F(v):       # free-energy L2 (~1e-2)
    return "--" if not np.isfinite(v) else f"{v:.4f}"


def f_int(v):     # integrated L2(F) (~15)
    return "--" if not np.isfinite(v) else f"{v:.2f}"


def f_Fp(v):      # mean-force L2 (~7e-2)
    return "--" if not np.isfinite(v) else f"{v:.4f}"


def f_prob(v):
    return "--" if not np.isfinite(v) else f"{v:.2f}"


def f_evt(v):     # FR event fraction (~1e-5); LaTeX math
    if not np.isfinite(v) or v == 0:
        return "$0$"
    exp = int(np.floor(np.log10(abs(v))))
    mant = v / 10 ** exp
    return f"${mant:.1f}\\times10^{{{exp}}}$"


def f_barrier(v):  # barrier crossings (~5e5); LaTeX math, x10^5
    return "--" if not np.isfinite(v) else f"${v/1e5:.2f}\\times10^{{5}}$"


def f_gamma(v):
    return f"{v:g}"


def latex_table(df, columns, headers, fmts, caption, label, col_align=None,
                note=None):
    """Render a booktabs LaTeX table from ``df`` rows."""
    ncol = len(columns)
    align = col_align or ("l" + "r" * (ncol - 1))
    lines = [r"\begin{table}[t]", r"\centering",
             r"\small",
             rf"\caption{{{caption}}}",
             rf"\label{{{label}}}",
             rf"\begin{{tabular}}{{{align}}}",
             r"\toprule",
             " & ".join(headers) + r" \\",
             r"\midrule"]
    for _, row in df.iterrows():
        cells = [fmts[c](row[c]) if c in fmts else str(row[c]) for c in columns]
        lines.append(" & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    if note:
        lines.append(rf"\par\smallskip\footnotesize {note}")
    lines.append(r"\end{table}")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def fig01_geometry(ref_dir):
    prof_path = os.path.join(ref_dir, "reference_profile.csv")
    npz_path = os.path.join(ref_dir, "reference_grid.npz")
    if not (os.path.exists(prof_path) and os.path.exists(npz_path)):
        _warn("reference files missing; skipping fig01.")
        return
    prof = pd.read_csv(prof_path)
    g = np.load(npz_path)
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.4))
    cf = ax[0, 0].contourf(g["x_grid"], g["y_grid"], g["V_grid"], levels=40,
                           cmap="RdYlBu_r")
    ax[0, 0].contour(g["x_grid"], g["y_grid"], g["V_grid"], levels=12,
                     colors="k", linewidths=0.4, alpha=0.4)
    ax[0, 0].set_title(r"(a) Potential $V(x,y)$")
    fig.colorbar(cf, ax=ax[0, 0])
    cf = ax[0, 1].contourf(g["x_grid"], g["y_grid"], g["rho_grid"], levels=40,
                           cmap="viridis")
    ax[0, 1].set_title(r"(b) Canonical density $\pi\propto e^{-\beta V}$")
    fig.colorbar(cf, ax=ax[0, 1])
    for a in (ax[0, 0], ax[0, 1]):
        a.set_xlabel("x"); a.set_ylabel("y"); a.grid(False)
    ax[1, 0].plot(prof["x"], prof["F_ref"], "k", lw=2)
    ax[1, 0].axvline(0, color="gray", ls=":", lw=0.8)
    ax[1, 0].set_title(r"(c) Reference free energy $F_{\rm ref}(x)$")
    ax[1, 0].set_xlabel("x"); ax[1, 0].set_ylabel(r"$F_{\rm ref}$")
    ax[1, 1].plot(prof["x"], prof["Fprime_ref"], "k", lw=2)
    ax[1, 1].axhline(0, color="gray", ls="--", lw=0.7)
    ax[1, 1].set_title(r"(d) Reference mean force $F'_{\rm ref}(x)$")
    ax[1, 1].set_xlabel("x"); ax[1, 1].set_ylabel(r"$F'_{\rm ref}$")
    for a in (ax[1, 0], ax[1, 1]):
        a.set_xlim(-3, 3)
    _save(fig, "fig01_reference_geometry.png")


def _median_iqr_band(ax, sub, col, color, label, logy=False):
    """Median (line) and IQR (band) over seeds vs time for one config."""
    grp = sub.groupby("t")[col]
    t = np.array(sorted(sub["t"].unique()))
    q50 = grp.median().reindex(t).values
    q25 = grp.quantile(0.25).reindex(t).values
    q75 = grp.quantile(0.75).reindex(t).values
    ax.plot(t, q50, color=color, label=label)
    ax.fill_between(t, q25, q75, color=color, alpha=0.18, linewidth=0)
    if logy:
        ax.set_yscale("log")


def fig02_convergence(long, selected):
    """Median+IQR over seeds for the selected (best-by-integrated) config/method."""
    panels = [("l2_F", r"$\|\hat F_t-F_{\rm ref}\|_{L^2}$", True),
              ("l2_Fprime", r"$\|\hat F'_t-F'_{\rm ref}\|_{L^2}$", True),
              ("barrier_crossings", "cumulative barrier crossings", False),
              ("fr_event_fraction", "FR event fraction", False)]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.6))
    for ax, (col, ylab, logy) in zip(axes.ravel(), panels):
        if col not in long.columns:
            ax.text(0.5, 0.5, f"'{col}' missing", ha="center"); ax.set_axis_off()
            continue
        for m in METHOD_ORDER:
            cid = selected.get(m)
            if cid is None:
                continue
            sub = long[long["config_id"] == cid]
            if sub.empty:
                continue
            _median_iqr_band(ax, sub, col, METHOD_COLOR[m], METHOD_LABEL[m], logy)
        ax.set_xlabel("time $t$"); ax.set_ylabel(ylab); ax.set_title(ylab)
        ax.legend(fontsize=8)
    fig.suptitle("Convergence of the selected (best integrated $L^2(F)$) "
                 "configuration: median and IQR over 5 seeds", y=1.005)
    _save(fig, "fig02_convergence_curves.png")


def _center_window(x, y, lo=EVAL_LO, hi=EVAL_HI):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = (x >= lo) & (x <= hi)
    return y - (np.mean(y[m]) if np.any(m) else np.mean(y))


def fig03_profiles(prof, ref, selected):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    for m in METHOD_ORDER:
        cid = selected.get(m)
        sub = prof[prof["config_id"] == cid] if cid else prof[prof["method"] == m]
        if sub.empty:
            continue
        x = np.array(sorted(sub["x"].unique()))
        gp = sub.groupby("x")
        Fp = gp["Fprime_hat"].mean().reindex(x).values
        F = _center_window(x, gp["F_hat"].mean().reindex(x).values)
        p = gp["p_hat"].mean().reindex(x).values
        c = METHOD_COLOR[m]
        axes[0].plot(x, Fp, color=c, label=METHOD_LABEL[m])
        axes[1].plot(x, F, color=c, label=METHOD_LABEL[m])
        axes[2].plot(x, p, color=c, label=METHOD_LABEL[m])
    axes[0].plot(ref["x"], ref["Fprime_ref"], "k--", lw=1.4, label="reference")
    axes[1].plot(ref["x"], _center_window(ref["x"], ref["F_ref"]), "k--",
                 lw=1.4, label="reference")
    axes[2].plot(ref["x"], ref["p_ref"], "k:", lw=1.3,
                 label=r"$p_{\rm ref}$ (unbiased)")
    for ax, title, yl in zip(
            axes, [r"(a) Mean force $\hat F'(x)$", r"(b) Free energy $\hat F(x)$",
                   r"(c) $x$-marginal $\hat p(x)$"],
            [r"$F'$", r"$F$", "density"]):
        ax.set_xlabel("x"); ax.set_ylabel(yl); ax.set_title(title)
        ax.set_xlim(-2.7, 2.7); ax.legend(fontsize=8)
    _save(fig, "fig03_final_profiles.png")


def fig04_target_ablation(fs):
    """Per-target distribution over configs (seed-medians), with ABF baseline."""
    metrics = [("integrated_l2_F", r"$\int_0^T\|\hat F_t-F_{\rm ref}\|\,dt$"),
               ("final_l2_F", r"final $\|\hat F_T-F_{\rm ref}\|$"),
               ("final_l2_Fprime", r"final $\|\hat F'_T-F'_{\rm ref}\|$")]
    fr_methods = ["abf_fr_estimated", "abf_fr_uniform", "abf_fr_oracle"]
    abf = fs[fs["method"] == "abf_only"]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    for ax, (col, title) in zip(axes, metrics):
        data, labels, colors = [], [], []
        for m in fr_methods:
            g = fs[fs["method"] == m]
            per_cfg = g.groupby("config_id")[col].median().values
            data.append(per_cfg[np.isfinite(per_cfg)])
            labels.append(METHOD_LABEL_SHORT[m]); colors.append(METHOD_COLOR[m])
        bp = ax.boxplot(data, tick_labels=labels, patch_artist=True,
                        showmeans=True)
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.55)
        if not abf.empty:
            base = med(abf[col])
            ax.axhline(base, color=METHOD_COLOR["abf_only"], ls="--", lw=1.4,
                       label="ABF only (median)")
            ax.legend(fontsize=8)
        ax.set_title(title); ax.set_ylabel(title)
        for tk in ax.get_xticklabels():
            tk.set_rotation(15); tk.set_ha("right")
    fig.suptitle("Target-type ablation: distribution over the 36 configurations "
                 "per target (each point = seed-median)", y=1.02)
    _save(fig, "fig04_target_ablation.png")


def fig05_mechanism(long, fr_events, selected):
    """Early-time / redistribution mechanism for the estimated FR vs ABF-only."""
    est = selected.get("abf_fr_estimated")
    abf = selected.get("abf_only")
    sel = [("ABF only", abf, METHOD_COLOR["abf_only"]),
           ("ABF+FR (estimated)", est, METHOD_COLOR["abf_fr_estimated"])]
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    # (a) early-time L2(F)
    ax = axes[0, 0]
    for lab, cid, c in sel:
        if cid is None:
            continue
        sub = long[(long["config_id"] == cid) & (long["t"] <= 40)]
        grp = sub.groupby("t")["l2_F"]
        t = np.array(sorted(sub["t"].unique()))
        ax.plot(t, grp.median().reindex(t).values, color=c, label=lab, marker="o",
                ms=3)
    ax.set_yscale("log"); ax.set_xlabel("time $t$ (early window)")
    ax.set_ylabel(r"$\|\hat F_t-F_{\rm ref}\|_{L^2}$")
    ax.set_title("(a) Early-time free-energy error"); ax.legend(fontsize=8)
    # (b) region fractions over time (estimated FR config)
    ax = axes[0, 1]
    regcols = [("left_frac", "left basin", "#1f77b4"),
               ("barrier_frac", "barrier", "#7f7f7f"),
               ("right_frac", "right basin", "#d62728")]
    if est is not None and all(c in long.columns for c, _, _ in regcols):
        sub = long[long["config_id"] == est]
        t = np.array(sorted(sub["t"].unique()))
        for col, lab, c in regcols:
            ax.plot(t, sub.groupby("t")[col].median().reindex(t).values,
                    color=c, label=lab)
        ax.axhline(1/3, color="k", ls=":", lw=0.8)
        ax.set_ylabel("particle fraction"); ax.set_xlabel("time $t$")
        ax.set_title("(b) Region occupancy (estimated FR)"); ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "region-fraction columns missing", ha="center")
        ax.set_axis_off()
    # (c) cumulative barrier crossings
    ax = axes[1, 0]
    for lab, cid, c in sel:
        if cid is None:
            continue
        sub = long[long["config_id"] == cid]
        t = np.array(sorted(sub["t"].unique()))
        ax.plot(t, sub.groupby("t")["barrier_crossings"].median().reindex(t).values,
                color=c, label=lab)
    ax.set_xlabel("time $t$"); ax.set_ylabel("cumulative barrier crossings")
    ax.set_title("(c) Barrier crossings"); ax.legend(fontsize=8)
    # (d) FR score std over time (estimated FR)
    ax = axes[1, 1]
    scol = "score_std" if "score_std" in long.columns else (
        "fr_score_std" if "fr_score_std" in long.columns else None)
    if est is not None and scol is not None:
        sub = long[long["config_id"] == est]
        t = np.array(sorted(sub["t"].unique()))
        ax.plot(t, sub.groupby("t")[scol].median().reindex(t).values,
                color=METHOD_COLOR["abf_fr_estimated"], label="score std")
        ax.set_xlabel("time $t$"); ax.set_ylabel(r"$\mathrm{std}_i\, S_i$")
        ax.set_title("(d) FR score dispersion (estimated FR)"); ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "score-std column missing", ha="center")
        ax.set_axis_off()
    _save(fig, "fig05_mechanism.png")


def fig06_conditional(cond):
    if cond is None or cond.empty:
        _warn("conditional diagnostics missing; skipping fig06.")
        return
    t_final = cond["t"].max()
    c = cond[cond["t"] == t_final]
    bins = sorted(c["x_bin_center"].unique())
    methods = [m for m in METHOD_ORDER if m in set(c["method"])]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.7))
    width = 0.8 / max(len(methods), 1)
    for i, m in enumerate(methods):
        sub = c[c["method"] == m]
        l2 = [np.nanmean(sub[sub["x_bin_center"] == b]["conditional_l2_y"]) for b in bins]
        kl = [np.nanmean(sub[sub["x_bin_center"] == b]["conditional_kl_y"]) for b in bins]
        xs = np.arange(len(bins)) + i * width
        axes[0].bar(xs, l2, width=width, color=METHOD_COLOR[m], label=METHOD_LABEL[m])
        axes[1].bar(xs, kl, width=width, color=METHOD_COLOR[m], label=METHOD_LABEL[m])
    for ax, ylab, title in zip(
            axes,
            [r"$\|\hat p(y|x)-p_{\rm ref}(y|x)\|_{L^2}$",
             r"$D_{\rm KL}(\hat p(y|x)\,\|\,p_{\rm ref}(y|x))$"],
            [r"(a) Conditional $Y\,|\,X$ $L^2$ error",
             r"(b) Conditional $Y\,|\,X$ KL divergence"]):
        ax.set_xticks(np.arange(len(bins)) + 0.4 - width / 2)
        ax.set_xticklabels([f"{b:g}" for b in bins])
        ax.set_xlabel("$x$-bin centre"); ax.set_ylabel(ylab); ax.set_title(title)
        ax.legend(fontsize=8)
    fig.suptitle(r"Conditional sampling fidelity at $T$ "
                 r"(aggregated over configs and seeds; not config-resolved)",
                 y=1.02)
    _save(fig, "fig06_conditional_y_given_x.png")


def fig07_hyperparameters(cs):
    """Estimated-target sensitivity: gamma x eta and gamma x burn-in heatmaps."""
    sub = cs[cs["target_type"] == "estimated"]
    if sub.empty:
        _warn("no estimated configs; skipping fig07.")
        return
    val = "median_integrated_l2_F"
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    for ax, (idx, lab) in zip(axes, [("eta", r"$\eta$"),
                                     ("burnin_fraction", "burn-in fraction")]):
        piv = (sub.groupby([idx, "gamma"])[val].min().reset_index()
               .pivot(index=idx, columns="gamma", values=val))
        im = ax.imshow(piv.values, origin="lower", aspect="auto", cmap="viridis_r")
        ax.set_xticks(range(len(piv.columns)))
        ax.set_xticklabels([f"{c:g}" for c in piv.columns])
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels([f"{r:g}" for r in piv.index])
        ax.set_xlabel(r"$\gamma$"); ax.set_ylabel(lab); ax.grid(False)
        vmid = np.nanmedian(piv.values)
        for ii in range(piv.shape[0]):
            for jj in range(piv.shape[1]):
                v = piv.values[ii, jj]
                if np.isfinite(v):
                    ax.text(jj, ii, f"{v:.2f}", ha="center", va="center",
                            fontsize=8, color="white" if v < vmid else "black")
        fig.colorbar(im, ax=ax, label=r"median $\int\|\hat F-F_{\rm ref}\|\,dt$")
    fig.suptitle("Estimated-target FR hyperparameter sensitivity "
                 "(lower = faster convergence)", y=1.02)
    _save(fig, "fig07_hyperparameter_heatmaps.png")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-root",
                    default=os.path.join(REPO_ROOT, "results", "two_dim_xi_x"))
    args = ap.parse_args(argv)

    _set_style()
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(TAB_DIR, exist_ok=True)
    prod = os.path.join(args.results_root, "production_gpu")
    ref_dir = os.path.join(args.results_root, "reference")

    print("[build_report_assets] loading merged CSVs ...")
    fs, fs_src = load_merged(prod, "production_gpu_final_summary")
    cs, _ = load_merged(prod, "production_gpu_config_summary")
    long, _ = load_merged(prod, "production_gpu_runs_long")
    prof, _ = load_merged(prod, "production_gpu_profiles")
    try:
        fr_events, _ = load_merged(prod, "production_gpu_fr_events")
    except FileNotFoundError:
        fr_events = None
    try:
        cond, _ = load_merged(prod, "production_gpu_conditional_diagnostics")
    except FileNotFoundError:
        cond = None
    ref = pd.read_csv(os.path.join(ref_dir, "reference_profile.csv"))

    # ----------------------------------------------------------------- checks
    methods = sorted(fs["method"].unique())
    seeds = sorted(int(s) for s in fs["seed"].unique())
    n_runs = len(fs)
    n_configs = fs["config_id"].nunique()
    key_cols = ["final_l2_F", "final_l2_Fprime", "integrated_l2_F"]
    n_nan = int(fs[key_cols].isna().to_numpy().sum())
    n_anynan_flag = int(fs["any_nan"].sum()) if "any_nan" in fs.columns else 0
    print(f"  runs={n_runs} configs={n_configs} seeds={seeds} methods={methods}")
    print(f"  NaNs in key metrics={n_nan} ; any_nan flags set={n_anynan_flag}")
    assert n_nan == 0, "NaNs present in key metrics!"
    assert n_anynan_flag == 0, "any_nan flag set on some run!"

    abf = fs[fs["method"] == "abf_only"]
    abf_by_seed = abf.set_index("seed")["final_l2_F"].to_dict()
    abf_cid = abf["config_id"].iloc[0]

    fr_methods = ["abf_fr_estimated", "abf_fr_uniform", "abf_fr_oracle"]
    sel_int = {"abf_only": abf_cid}
    sel_fin = {"abf_only": abf_cid}
    rec_int, rec_fin = {}, {}
    rec_int["abf_only"] = config_record(fs, abf_cid, abf_by_seed)
    rec_fin["abf_only"] = rec_int["abf_only"]
    for m in fr_methods:
        g = fs[fs["method"] == m]
        ci = best_config(g, INTEGRATED)
        cf = best_config(g, FINAL)
        sel_int[m] = ci
        sel_fin[m] = cf
        rec_int[m] = config_record(fs, ci, abf_by_seed)
        rec_fin[m] = config_record(fs, cf, abf_by_seed)

    # ----------------------------------------------------------------- tables
    def records_to_df(rec):
        return pd.DataFrame([rec[m] for m in METHOD_ORDER])

    df_int = records_to_df(rec_int)
    df_fin = records_to_df(rec_fin)

    csv_cols = ["method", "target_type", "gamma", "eta", "burnin_fraction",
                "fr_every", "n_seeds", "median_integrated_l2_F",
                "median_final_l2_F", "median_final_l2_Fprime",
                "median_barrier_crossings", "median_mean_fr_event_fraction",
                "prob_beats_abf"]
    df_int[csv_cols].to_csv(os.path.join(TAB_DIR, "best_configs_by_integrated_F.csv"),
                            index=False)
    df_fin[csv_cols].to_csv(os.path.join(TAB_DIR, "best_configs_by_final_F.csv"),
                            index=False)

    # LaTeX: main results table (selection by integrated L2(F)).
    df_int = df_int.copy()
    df_int["label"] = df_int["method"].map(METHOD_LABEL_SHORT)
    df_fin = df_fin.copy()
    df_fin["label"] = df_fin["method"].map(METHOD_LABEL_SHORT)

    tbl_cols = ["label", "median_integrated_l2_F", "median_final_l2_F",
                "median_final_l2_Fprime", "prob_beats_abf",
                "median_mean_fr_event_fraction"]
    tbl_head = [r"Method", r"$\int_0^T\!\|\hat F_t-F_{\rm ref}\|\,dt$",
                r"$\|\hat F_T-F_{\rm ref}\|$", r"$\|\hat F'_T-F'_{\rm ref}\|$",
                r"$P(\text{beats ABF})$", r"FR event frac."]
    fmts = {"median_integrated_l2_F": f_int, "median_final_l2_F": f_F,
            "median_final_l2_Fprime": f_Fp, "prob_beats_abf": f_prob,
            "median_mean_fr_event_fraction": f_evt}
    note_oracle = (r"$^{*}$\,The oracle target uses the reference free energy "
                   r"$F_{\rm ref}$ (the unknown being computed); it is a "
                   r"diagnostic positive control, \emph{not} a usable algorithm. "
                   r"Per-config medians over $5$ seeds; one configuration per "
                   r"method, selected by the smallest median integrated $L^2(F)$.")
    with open(os.path.join(TAB_DIR, "best_configs_by_integrated_F.tex"), "w") as fh:
        fh.write(latex_table(
            df_int, tbl_cols, tbl_head, fmts,
            caption=(r"Main comparison, configurations selected by smallest "
                     r"median \emph{integrated} $L^2(F)$ (convergence speed). "
                     r"Lower is better except $P(\text{beats ABF})$."),
            label="tab:main_integrated",
            col_align="lrrrrr", note=note_oracle))
    note_final = (r"$^{*}$\,Oracle is a diagnostic control (uses $F_{\rm ref}$), "
                  r"not a deployable method. Configurations here are selected by "
                  r"smallest median \emph{final} $L^2(F)$; cf.\ "
                  r"Table~\ref{tab:main_integrated}, which need not pick the same "
                  r"hyperparameters.")
    with open(os.path.join(TAB_DIR, "best_configs_by_final_F.tex"), "w") as fh:
        fh.write(latex_table(
            df_fin, tbl_cols, tbl_head, fmts,
            caption=(r"Supplement: configurations selected by smallest median "
                     r"\emph{final} $L^2(F)$ (final-budget accuracy)."),
            label="tab:main_final",
            col_align="lrrrrr", note=note_final))

    # Selected-configuration hyperparameters table (integrated rule).
    hp_cols = ["label", "target_type", "gamma", "eta", "burnin_fraction",
               "fr_every", "n_seeds"]
    hp_head = [r"Method", r"target", r"$\gamma$", r"$\eta$", r"burn-in",
               r"$n_{\rm FR}$", r"seeds"]
    hp_fmts = {"gamma": f_gamma, "eta": f_gamma, "burnin_fraction": f_gamma,
               "fr_every": lambda v: f"{int(v)}", "n_seeds": lambda v: f"{int(v)}"}
    with open(os.path.join(TAB_DIR, "selected_configs.tex"), "w") as fh:
        fh.write(latex_table(
            df_int, hp_cols, hp_head, hp_fmts,
            caption=(r"Hyperparameters of the selected configurations "
                     r"(integrated-$L^2(F)$ rule). $n_{\rm FR}$ is the "
                     r"Fisher--Rao stride \texttt{fr\_every}."),
            label="tab:selected_configs", col_align="llrrrrr"))

    # Experimental-design summary table.
    fr_only = fs[fs["method"] != "abf_only"]
    design_rows = [
        ("Inverse temperature $\\beta$", "4.0"),
        ("Time step $\\Delta t$", "0.002"),
        ("Steps per run", "100{,}000"),
        ("Particles per run $N$", "1000"),
        ("Seeds", f"{len(seeds)} ({', '.join(map(str, seeds))})"),
        ("Methods", "4 (ABF; FR estimated / uniform / oracle)"),
        ("FR rate $\\gamma$", ", ".join(f"{g:g}" for g in
                                        sorted(fr_only['gamma'].unique()))),
        ("FR bandwidth $\\eta$", ", ".join(f"{e:g}" for e in
                                           sorted(fr_only['eta'].unique()))),
        ("FR burn-in fraction", ", ".join(f"{b:g}" for b in
                                          sorted(fr_only['burnin_fraction'].unique()))),
        ("FR stride \\texttt{fr\\_every}", ", ".join(f"{int(v)}" for v in
                                          sorted(fr_only['fr_every'].unique()))),
        ("Total configurations", str(n_configs)),
        ("Total runs", str(n_runs)),
    ]
    lines = [r"\begin{table}[t]", r"\centering", r"\small",
             r"\caption{Production study design (read from the result CSVs and "
             r"\texttt{configs/two\_dim\_xi\_x\_production\_gpu.yaml}).}",
             r"\label{tab:design}", r"\begin{tabular}{ll}", r"\toprule",
             r"Quantity & Value \\", r"\midrule"]
    for k, v in design_rows:
        lines.append(f"{k} & {v} " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    with open(os.path.join(TAB_DIR, "design_summary.tex"), "w") as fh:
        fh.write("\n".join(lines))

    # ----------------------------------------------------------------- numbers
    def pct(a, b):  # percent reduction from a (baseline) to b
        return 100.0 * (a - b) / a if a else float("nan")

    abf_r = rec_int["abf_only"]
    est_i, est_f = rec_int["abf_fr_estimated"], rec_fin["abf_fr_estimated"]
    uni_i, uni_f = rec_int["abf_fr_uniform"], rec_fin["abf_fr_uniform"]
    ora = rec_int["abf_fr_oracle"]

    numbers = dict(
        meta=dict(n_runs=n_runs, n_configs=n_configs, n_seeds=len(seeds),
                  n_methods=len(methods), seeds=seeds, methods=methods,
                  n_nan=n_nan, source=fs_src,
                  beta=4.0, dt=0.002, n_steps=100000, n_particles=1000,
                  n_fr_configs_per_target=int(
                      fr_only[fr_only["method"] == "abf_fr_estimated"]
                      ["config_id"].nunique()),
                  gamma_values=sorted(float(g) for g in fr_only["gamma"].unique()),
                  eta_values=sorted(float(e) for e in fr_only["eta"].unique()),
                  burnin_values=sorted(float(b) for b in
                                       fr_only["burnin_fraction"].unique())),
        abf_only=abf_r,
        estimated_by_integrated=est_i, estimated_by_final=est_f,
        uniform_by_integrated=uni_i, uniform_by_final=uni_f,
        oracle_best=ora,
        improvements=dict(
            est_int_integratedF_pct=pct(abf_r["median_integrated_l2_F"],
                                        est_i["median_integrated_l2_F"]),
            est_fin_finalF_pct=pct(abf_r["median_final_l2_F"],
                                   est_f["median_final_l2_F"]),
            est_fin_finalFp_pct=pct(abf_r["median_final_l2_Fprime"],
                                    est_f["median_final_l2_Fprime"]),
            uni_int_integratedF_pct=pct(abf_r["median_integrated_l2_F"],
                                        uni_i["median_integrated_l2_F"]),
            uni_fin_finalF_pct=pct(abf_r["median_final_l2_F"],
                                   uni_f["median_final_l2_F"]),
            oracle_integratedF_pct=pct(abf_r["median_integrated_l2_F"],
                                       ora["median_integrated_l2_F"]),
            oracle_finalF_pct=pct(abf_r["median_final_l2_F"],
                                  ora["median_final_l2_F"]),
            est_vs_oracle_integrated_gap_pct=pct(ora["median_integrated_l2_F"],
                                                 est_i["median_integrated_l2_F"]),
        ),
    )
    # EB + WCA case assets (figures, tables) and their CSV-derived numbers.
    print("[build_report_assets] building EB + WCA case assets ...")
    meta_int_pct = numbers["improvements"]["est_int_integratedF_pct"]
    case_macros, case_json = report_cases.build_cases(FIG_DIR, TAB_DIR, meta_int_pct)
    numbers.update(case_json)

    with open(os.path.join(TAB_DIR, "report_numbers.json"), "w") as fh:
        json.dump(numbers, fh, indent=2)
    print("   [num] wrote tables/report_numbers.json")
    def m_F(v): return f"{v:.4f}"
    def m_int(v): return f"{v:.2f}"
    def m_Fp(v): return f"{v:.4f}"
    def m_pct(v): return f"{v:.1f}"
    def m_evt(v):
        if v == 0:
            return r"0"
        exp = int(np.floor(np.log10(abs(v)))); mant = v / 10 ** exp
        return f"{mant:.1f}\\times10^{{{exp}}}"
    macros = {
        "NRuns": str(n_runs), "NConfigs": str(n_configs),
        "NSeeds": str(len(seeds)), "NMethods": str(len(methods)),
        "NSteps": "100{,}000", "NParticles": "1000", "BetaVal": "4",
        "NFRconfigs": str(numbers["meta"]["n_fr_configs_per_target"]),
        "AbfIntegF": m_int(abf_r["median_integrated_l2_F"]),
        "AbfFinalF": m_F(abf_r["median_final_l2_F"]),
        "AbfFinalFp": m_Fp(abf_r["median_final_l2_Fprime"]),
        "AbfBarrier": f_barrier(abf_r["median_barrier_crossings"]).strip("$"),
        # estimated (by integrated)
        "EstIntegF": m_int(est_i["median_integrated_l2_F"]),
        "EstIntFinalF": m_F(est_i["median_final_l2_F"]),
        "EstIntGamma": f"{est_i['gamma']:g}",
        "EstIntBurnin": f"{est_i['burnin_fraction']:g}",
        "EstIntProb": m_pct(100 * est_i["prob_beats_abf"]),
        "EstIntEvt": m_evt(est_i["median_mean_fr_event_fraction"]),
        "EstIntegImpPct": m_pct(numbers["improvements"]["est_int_integratedF_pct"]),
        # estimated (by final)
        "EstFinFinalF": m_F(est_f["median_final_l2_F"]),
        "EstFinIntegF": m_int(est_f["median_integrated_l2_F"]),
        "EstFinFinalFp": m_Fp(est_f["median_final_l2_Fprime"]),
        "EstFinGamma": f"{est_f['gamma']:g}",
        "EstFinBurnin": f"{est_f['burnin_fraction']:g}",
        "EstFinProb": m_pct(100 * est_f["prob_beats_abf"]),
        "EstFinEvt": m_evt(est_f["median_mean_fr_event_fraction"]),
        "EstFinFinalImpPct": m_pct(numbers["improvements"]["est_fin_finalF_pct"]),
        "EstFinFinalFpImpPct": m_pct(numbers["improvements"]["est_fin_finalFp_pct"]),
        # uniform
        "UniIntegF": m_int(uni_i["median_integrated_l2_F"]),
        "UniIntegImpPct": m_pct(numbers["improvements"]["uni_int_integratedF_pct"]),
        "UniFinalF": m_F(uni_f["median_final_l2_F"]),
        "UniFinFinalImpPct": m_pct(numbers["improvements"]["uni_fin_finalF_pct"]),
        # oracle (diagnostic)
        "OracIntegF": m_int(ora["median_integrated_l2_F"]),
        "OracFinalF": m_F(ora["median_final_l2_F"]),
        "OracFinalFp": m_Fp(ora["median_final_l2_Fprime"]),
        "OracProb": m_pct(100 * ora["prob_beats_abf"]),
        "OracIntegImpPct": m_pct(numbers["improvements"]["oracle_integratedF_pct"]),
        "OracFinalImpPct": m_pct(numbers["improvements"]["oracle_finalF_pct"]),
        "EstOracGapPct": m_pct(numbers["improvements"]["est_vs_oracle_integrated_gap_pct"]),
    }
    macros.update(case_macros)   # EB + WCA case macros (namespaced)
    with open(os.path.join(TAB_DIR, "numbers.tex"), "w") as fh:
        for k, v in macros.items():
            fh.write(rf"\newcommand{{\{k}}}{{{v}}}" + "\n")
    print("   [num] wrote tables/numbers.tex")

    # ----------------------------------------------------------------- figures
    print("[build_report_assets] rendering figures ...")
    fig01_geometry(ref_dir)
    fig02_convergence(long, sel_int)
    fig03_profiles(prof, ref, sel_int)
    fig04_target_ablation(fs)
    fig05_mechanism(long, fr_events, sel_int)
    fig06_conditional(cond)
    fig07_hyperparameters(cs)

    print("[build_report_assets] done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
