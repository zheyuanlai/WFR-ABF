#!/usr/bin/env python3
"""Report-quality plots for the WCA dimer production study (Matplotlib only).

Reads results/wca_production/summaries/* and writes PNGs to
results/wca_production/plots/. Figures:
  fig1_convergence.png        L2(F),L2(Fp),repl-fraction vs time (main stage)
  fig2_final_profiles.png     F(z),F'(z),p(z) vs reference (seed-mean)
  fig3_seed_performance.png   per-seed final L2(F), integrated L2(F), final L2(Fp)
  fig4_mechanism.png          N_eff(z), per-bin F'(z) error, birth/death, fractions
  fig5_failure_boundary.png   fr_rate vs final L2(F) / ancestor ESS / repl fraction
  fig6_difficulty.png         % gain vs budget / replicas / crowding

Each figure is skipped (with a message) if its inputs are absent, so partial
results still plot. Usage:
  python scripts/plot_wca_production.py
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import wca_jobs  # noqa: E402

# stable per-method colours/labels
STYLE = {
    "abf": ("#444444", "ABF only"),
    "fr_est_gentle": ("#74add1", "FR est. gentle (rate 0.05)"),
    "fr_est_tuned": ("#1f78b4", "FR est. tuned (rate 0.10)"),
    "fr_est_strong": ("#6a3d9a", "FR est. strong (rate 0.20)"),
    "fr_est_aggressive": ("#e31a1c", "FR est. aggressive (rate 0.50)"),
    "fr_uniform": ("#33a02c", "FR uniform target"),
    "fr_oracle": ("#ff7f00", "FR oracle target (diagnostic)"),
}


def style_for(name):
    return STYLE.get(name, ("#888888", name))


def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh))


def f(row, key, default=np.nan):
    try:
        return float(row[key])
    except (KeyError, ValueError, TypeError):
        return default


def load_timeseries(summ_dir, stage):
    """Return {name: dict(t, l2_f_median/q25/q75, l2_fp_*, repl_med, method)}."""
    rows = read_csv(os.path.join(summ_dir, "timeseries_summary.csv"))
    series = defaultdict(lambda: defaultdict(list))
    meta = {}
    for r in rows:
        if r["stage"] != stage:
            continue
        nm = r["name"]
        meta[nm] = r["method"]
        for col in ["t", "l2_f_median", "l2_f_q25", "l2_f_q75",
                    "l2_fp_median", "l2_fp_q25", "l2_fp_q75", "repl_cumulative_median"]:
            series[nm][col].append(f(r, col))
    out = {}
    for nm, d in series.items():
        arr = {k: np.asarray(v, float) for k, v in d.items()}
        order = np.argsort(arr["t"])
        out[nm] = {k: v[order] for k, v in arr.items()}
        out[nm]["method"] = meta[nm]
    return out


def fig1_convergence(summ_dir, plot_dir, stage="main"):
    ts = load_timeseries(summ_dir, stage)
    if not ts:
        print(f"  [fig1] no timeseries for stage={stage}; skip")
        return
    show = ["abf", "fr_est_tuned", "fr_uniform", "fr_oracle", "fr_est_aggressive"]
    names = [n for n in show if n in ts] or list(ts.keys())
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for nm in names:
        c, lab = style_for(nm)
        d = ts[nm]
        axes[0].plot(d["t"], d["l2_f_median"], color=c, label=lab, lw=1.8)
        axes[0].fill_between(d["t"], d["l2_f_q25"], d["l2_f_q75"], color=c, alpha=0.15)
        axes[1].plot(d["t"], d["l2_fp_median"], color=c, lw=1.8)
        axes[1].fill_between(d["t"], d["l2_fp_q25"], d["l2_fp_q75"], color=c, alpha=0.15)
        if str(d["method"]) != "abf" and np.nanmax(d["repl_cumulative_median"]) > 0:
            axes[2].plot(d["t"], d["repl_cumulative_median"], color=c, lw=1.8, label=lab)
    axes[0].set(xlabel="time", ylabel=r"$L^2(F)$", title="Free-energy error vs time (log y)")
    axes[1].set(xlabel="time", ylabel=r"$L^2(F')$", title="Mean-force error vs time (log y)")
    axes[2].set(xlabel="time", ylabel="cumulative replacements", title="Birth-death activity")
    axes[0].set_yscale("log")   # warm-up transient ~1.0 dwarfs the converged ~0.04 on a linear axis
    axes[1].set_yscale("log")
    axes[0].legend(fontsize=7.5, loc="upper right")
    if axes[2].get_legend_handles_labels()[0]:
        axes[2].legend(fontsize=7.5, loc="upper left")
    for ax in axes:
        ax.grid(alpha=0.3, which="both")
    fig.suptitle(f"Figure 1: convergence (stage={stage}, median + IQR over seeds)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(plot_dir, "fig1_convergence.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  wrote {out}")


def load_profiles(summ_dir):
    path = os.path.join(summ_dir, "profiles_summary.npz")
    if not os.path.exists(path):
        return None
    return dict(np.load(path, allow_pickle=True))


def _profile_key(prof, stage, name, field, n_steps=None, n_replicas=None, a=None):
    """Find the first key matching stage|name|method|...|@field."""
    for k in prof:
        if "@" not in k:
            continue
        head, fld = k.rsplit("@", 1)
        parts = head.split("|")
        if fld != field or parts[0] != stage or parts[1] != name:
            continue
        if n_steps is not None and int(parts[3]) != int(n_steps):
            continue
        if n_replicas is not None and int(parts[4]) != int(n_replicas):
            continue
        return k
    return None


def fig2_final_profiles(summ_dir, plot_dir, stage="main"):
    prof = load_profiles(summ_dir)
    if prof is None:
        print("  [fig2] no profiles_summary.npz; skip")
        return
    grid = prof["grid"]
    refF, refFp = prof["ref_free_energy"], prof["ref_mean_force"]
    show = ["abf", "fr_est_tuned", "fr_uniform", "fr_oracle"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    axes[0].plot(grid, refF, "k--", lw=2, label="TI reference")
    axes[1].plot(grid, refFp, "k--", lw=2, label="TI reference")
    plotted = 0
    for nm in show:
        kF = _profile_key(prof, stage, nm, "F")
        if kF is None:
            continue
        plotted += 1
        c, lab = style_for(nm)
        axes[0].plot(grid, prof[kF], color=c, lw=1.6, label=lab)
        kFp = _profile_key(prof, stage, nm, "Fp")
        if kFp is not None:
            axes[1].plot(grid, prof[kFp], color=c, lw=1.6)
        kp = _profile_key(prof, stage, nm, "p")
        if kp is not None:
            axes[2].plot(grid, prof[kp], color=c, lw=1.6, label=lab)
    if plotted == 0:
        print(f"  [fig2] no profiles for stage={stage}; skip")
        plt.close(fig)
        return
    for ax in (axes[0], axes[1], axes[2]):
        ax.axvspan(0.25, 0.75, color="gray", alpha=0.08)
        ax.grid(alpha=0.3)
    axes[0].set(xlabel="z", ylabel="F(z)", title="Free energy vs reference")
    axes[1].set(xlabel="z", ylabel="F'(z)", title="Mean force vs reference")
    axes[2].set(xlabel="z", ylabel=r"$\hat p(z)$", title="RC marginal (shaded = transition)")
    axes[0].legend(fontsize=7.5)
    fig.suptitle(f"Figure 2: final profiles (stage={stage}, seed-mean)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(plot_dir, "fig2_final_profiles.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  wrote {out}")


def fig3_seed_performance(summ_dir, plot_dir, stage="main"):
    rows = [r for r in read_csv(os.path.join(summ_dir, "summary.csv")) if r["stage"] == stage]
    if not rows:
        print(f"  [fig3] no summary rows for stage={stage}; skip")
        return
    order = ["abf", "fr_est_gentle", "fr_est_tuned", "fr_est_strong",
             "fr_est_aggressive", "fr_uniform", "fr_oracle"]
    names = [n for n in order if any(r["name"] == n for r in rows)]
    metrics = [("l2_f", "final $L^2(F)$"), ("integrated_l2_f", "integrated $L^2(F)$"),
               ("l2_fp", "final $L^2(F')$")]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (key, title) in zip(axes, metrics):
        data, labs, cols = [], [], []
        for nm in names:
            vals = [f(r, key) for r in rows if r["name"] == nm]
            vals = [v for v in vals if np.isfinite(v)]
            if not vals:
                continue
            data.append(vals)
            labs.append(style_for(nm)[1].split(" (")[0])
            cols.append(style_for(nm)[0])
        if not data:
            continue
        bp = ax.boxplot(data, vert=True, patch_artist=True, showmeans=True,
                        medianprops=dict(color="black"))
        for patch, col in zip(bp["boxes"], cols):
            patch.set_facecolor(col)
            patch.set_alpha(0.5)
        for i, vals in enumerate(data, start=1):
            ax.scatter(np.full(len(vals), i) + np.random.uniform(-0.08, 0.08, len(vals)),
                       vals, s=12, color="black", alpha=0.6, zorder=3)
        ax.set_xticks(range(1, len(labs) + 1))
        ax.set_xticklabels(labs, rotation=35, ha="right", fontsize=7.5)
        ax.set(ylabel=title, title=title)
        ax.grid(alpha=0.3, axis="y")
    fig.suptitle(f"Figure 3: per-seed performance (stage={stage})", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(plot_dir, "fig3_seed_performance.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  wrote {out}")


def fig4_mechanism(summ_dir, plot_dir, stage="main"):
    prof = load_profiles(summ_dir)
    ts = load_timeseries(summ_dir, stage)
    if prof is None:
        print("  [fig4] no profiles; skip")
        return
    grid = prof["grid"]
    refFp = prof["ref_mean_force"]
    edges = prof["hist_edges"]
    centers = 0.5 * (edges[:-1] + edges[1:])
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    show = ["abf", "fr_est_tuned", "fr_oracle"]
    # (a) effective counts N_eff(z)
    for nm in show:
        k = _profile_key(prof, stage, nm, "neff")
        if k is None:
            continue
        c, lab = style_for(nm)
        axes[0, 0].plot(grid, prof[k], color=c, lw=1.6, label=lab)
    axes[0, 0].axvspan(0.25, 0.75, color="gray", alpha=0.08)
    axes[0, 0].set(xlabel="z", ylabel=r"$N_{\rm eff}(z)$ (kernel weight)",
                   title="(a) effective sample coverage")
    axes[0, 0].legend(fontsize=8)
    # (b) per-bin |F'(z) - F'_ref|
    for nm in show:
        k = _profile_key(prof, stage, nm, "Fp")
        if k is None:
            continue
        c, lab = style_for(nm)
        axes[0, 1].plot(grid, np.abs(prof[k] - refFp), color=c, lw=1.6, label=lab)
    axes[0, 1].axvspan(0.25, 0.75, color="gray", alpha=0.08)
    axes[0, 1].set(xlabel="z", ylabel=r"$|\hat F'(z)-F'_{\rm ref}(z)|$",
                   title="(b) per-bin mean-force error")
    # (c) birth/death locations for the tuned FR
    kb = _profile_key(prof, stage, "fr_est_tuned", "birth")
    kd = _profile_key(prof, stage, "fr_est_tuned", "death")
    if kb is not None and kd is not None:
        width = centers[1] - centers[0]
        axes[1, 0].bar(centers, prof[kb], width=width, color="#1f78b4", alpha=0.6, label="births (clone source)")
        axes[1, 0].bar(centers, -prof[kd], width=width, color="#e31a1c", alpha=0.6, label="deaths")
        axes[1, 0].axvspan(0.25, 0.75, color="gray", alpha=0.08)
        axes[1, 0].axhline(0, color="black", lw=0.8)
        axes[1, 0].legend(fontsize=8)
    axes[1, 0].set(xlabel="z", ylabel="count", title="(c) FR birth/death locations (tuned)")
    # (d) state fractions over time (tuned vs abf)
    rows = read_csv(os.path.join(summ_dir, "summary.csv"))
    have_frac = _plot_state_fractions(axes[1, 1], summ_dir, stage)
    if not have_frac:
        axes[1, 1].text(0.5, 0.5, "state fractions:\nsee raw frac_* arrays", ha="center")
    axes[1, 1].set(xlabel="time", ylabel="fraction in transition region",
                   title="(d) transition-region occupancy")
    fig.suptitle(f"Figure 4: mechanism diagnostics (stage={stage}, seed-mean)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(plot_dir, "fig4_mechanism.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  wrote {out}")


def _plot_state_fractions(ax, summ_dir, stage):
    """Plot transition-region fraction over time from raw npz (tuned vs abf)."""
    import glob
    raw_dir = os.path.join(os.path.dirname(summ_dir), "raw")
    series = defaultdict(list)
    for path in glob.glob(os.path.join(raw_dir, f"{stage}__*.npz")):
        d = wca_jobs.load_run(path)
        nm = str(d["name"]) if d["name"].ndim == 0 else str(d["name"].item())
        if nm not in ("abf", "fr_est_tuned"):
            continue
        series[nm].append((np.asarray(d["times"], float), np.asarray(d["frac_transition"], float)))
    if not series:
        return False
    for nm, runs in series.items():
        L = min(len(t) for t, _ in runs)
        t = runs[0][0][:L]
        fr = np.stack([r[1][:L] for r in runs])
        c, lab = style_for(nm)
        ax.plot(t, np.median(fr, axis=0), color=c, lw=1.8, label=lab)
        ax.fill_between(t, np.percentile(fr, 25, axis=0), np.percentile(fr, 75, axis=0),
                        color=c, alpha=0.15)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    return True


def fig5_failure_boundary(summ_dir, plot_dir, stage="failure"):
    """fr_rate ladder: final L2(F), ancestor ESS, replacement count vs fr_rate."""
    rows = [r for r in read_csv(os.path.join(summ_dir, "summary.csv"))
            if r["stage"] == stage and r["method"] == "fr_estimated"]
    abf_rows = [r for r in read_csv(os.path.join(summ_dir, "summary.csv"))
                if r["stage"] == stage and r["method"] == "abf"]
    if not rows:
        print(f"  [fig5] no failure-stage FR rows; skip")
        return
    by_rate = defaultdict(lambda: defaultdict(list))
    for r in rows:
        rate = f(r, "fr_rate")
        by_rate[rate]["l2_f"].append(f(r, "l2_f"))
        by_rate[rate]["ess"].append(f(r, "final_ancestor_ess"))
        by_rate[rate]["repl"].append(f(r, "total_replacement_events"))
    rates = sorted(by_rate.keys())
    med = lambda rate, key: np.nanmedian(by_rate[rate][key])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    l2 = [med(r, "l2_f") for r in rates]
    ess = [med(r, "ess") for r in rates]
    repl = [med(r, "repl") for r in rates]
    axes[0].plot(rates, l2, "o-", color="#1f78b4")
    if abf_rows:
        abf_med = np.nanmedian([f(r, "l2_f") for r in abf_rows])
        axes[0].axhline(abf_med, color="#444444", ls="--", label="ABF only")
        axes[0].legend(fontsize=8)
    axes[0].set(xlabel="fr_rate", ylabel=r"final $L^2(F)$", title="(a) accuracy vs fr_rate")
    axes[1].plot(rates, ess, "o-", color="#33a02c")
    axes[1].set(xlabel="fr_rate", ylabel="ancestor ESS", title="(b) diversity vs fr_rate")
    axes[2].plot(rates, repl, "o-", color="#e31a1c")
    axes[2].set(xlabel="fr_rate", ylabel="total replacements", title="(c) birth-death activity")
    for ax in axes:
        ax.set_xscale("log")
        ax.grid(alpha=0.3)
    fig.suptitle(f"Figure 5: failure boundary (stage={stage}, seed-median)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(plot_dir, "fig5_failure_boundary.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  wrote {out}")


def _gain_vs_axis(summ_dir, stage, axis_key):
    """Matched-seed % gain of fr_est_tuned vs abf, grouped by axis_key value."""
    rows = read_csv(os.path.join(summ_dir, "summary.csv"))
    rows = [r for r in rows if r["stage"] == stage]
    abf = {}
    for r in rows:
        if r["method"] == "abf":
            abf[(f(r, axis_key), int(float(r["seed"])))] = f(r, "l2_f")
    by_ax = defaultdict(list)
    for r in rows:
        if r["name"] != "fr_est_tuned":
            continue
        key = (f(r, axis_key), int(float(r["seed"])))
        if key in abf and abf[key] > 0:
            by_ax[f(r, axis_key)].append(100.0 * (abf[key] - f(r, "l2_f")) / abf[key])
    xs = sorted(by_ax.keys())
    med = [np.median(by_ax[x]) for x in xs]
    q25 = [np.percentile(by_ax[x], 25) for x in xs]
    q75 = [np.percentile(by_ax[x], 75) for x in xs]
    return xs, med, q25, q75


def fig6_difficulty(summ_dir, plot_dir):
    panels = [("difficulty_budget", "n_steps", "budget (steps)"),
              ("difficulty_replicas", "n_replicas", "number of replicas"),
              ("difficulty_crowding", "a", "lattice spacing a (smaller = denser)")]
    available = []
    for stage, axkey, _ in panels:
        xs, med, _, _ = _gain_vs_axis(summ_dir, stage, axkey)
        if xs:
            available.append((stage, axkey))
    if not available:
        print("  [fig6] no difficulty-stage data; skip")
        return
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for ax, (stage, axkey, xlabel) in zip(axes, panels):
        xs, med, q25, q75 = _gain_vs_axis(summ_dir, stage, axkey)
        if not xs:
            ax.text(0.5, 0.5, f"no data\n({stage})", ha="center", va="center")
            ax.set(title=xlabel)
            continue
        ax.plot(xs, med, "o-", color="#1f78b4")
        ax.fill_between(xs, q25, q75, color="#1f78b4", alpha=0.18)
        ax.axhline(0, color="black", lw=0.8, ls=":")
        ax.set(xlabel=xlabel, ylabel="% gain in final $L^2(F)$ vs ABF",
               title=f"FR tuned advantage vs {xlabel.split(' (')[0]}")
        ax.grid(alpha=0.3)
        if axkey == "a":
            ax.invert_xaxis()  # denser (smaller a) on the right
    fig.suptitle("Figure 6: difficulty dependence (fr_est_tuned vs ABF, matched-seed)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(plot_dir, "fig6_difficulty.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  wrote {out}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/wca_production.yaml")
    ap.add_argument("--summaries", default=None)
    ap.add_argument("--plots", default=None)
    ap.add_argument("--main-stage", default="main")
    args = ap.parse_args(argv)
    cfg = wca_jobs.load_yaml(args.config)
    summ_dir = args.summaries or os.path.join(cfg["output_root"], "summaries")
    plot_dir = args.plots or os.path.join(cfg["output_root"], "plots")
    os.makedirs(plot_dir, exist_ok=True)
    print(f"[plot] summaries={summ_dir} plots={plot_dir}")
    np.random.seed(0)  # jitter in fig3 only; deterministic
    fig1_convergence(summ_dir, plot_dir, stage=args.main_stage)
    fig2_final_profiles(summ_dir, plot_dir, stage=args.main_stage)
    fig3_seed_performance(summ_dir, plot_dir, stage=args.main_stage)
    fig4_mechanism(summ_dir, plot_dir, stage=args.main_stage)
    fig5_failure_boundary(summ_dir, plot_dir, stage="failure")
    fig6_difficulty(summ_dir, plot_dir)
    print("[plot] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
