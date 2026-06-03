#!/usr/bin/env python3
"""Plot the ABF--FR study results (tuning or eval stage).

Usage:
    python scripts/plot_abf_fr_study.py --stage tuning
    python scripts/plot_abf_fr_study.py --stage eval

Tuning figures (results/two_dim_xi_x/figures_tuning/):
    fig_tuning_gamma_eta_heatmap.png
    fig_tuning_gamma_burnin_heatmap.png
    fig_tuning_top_configs_boxplot.png

Eval figures (results/two_dim_xi_x/figures_eval/):
    fig01_reference_geometry.png
    fig02_eval_time_curves.png
    fig03_eval_final_profiles.png
    fig04_target_ablation.png
    fig05_good_vs_bad_fr.png
    fig06_conditional_y_given_x.png

Missing inputs produce a clear warning and that figure is skipped (no crash).
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from abffr import io_utils, plotting  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

DEFAULT_ROOT = "results/two_dim_xi_x"


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", required=True,
                   choices=["tuning", "eval", "smoke_gpu", "tuning_gpu",
                            "production_gpu"])
    p.add_argument("--output-root", default=DEFAULT_ROOT,
                   help="Experiment output root (default: %(default)s).")
    p.add_argument("--config", default=None,
                   help="Optional config (only used to read output_root).")
    return p.parse_args(argv)


def _warn(msg):
    print(f"[plot_abf_fr_study] WARNING: {msg}")


def _load(path):
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path)
    except Exception as exc:  # pragma: no cover
        _warn(f"could not read {os.path.relpath(path)}: {exc}")
        return None


# --------------------------------------------------------------------------- #
# Tuning figures
# --------------------------------------------------------------------------- #
def _pivot_heatmap(cfg_df, index, columns, value, target="estimated"):
    sub = cfg_df[cfg_df["target_type"] == target]
    if sub.empty:
        return None
    # Collapse any other free hyperparameters by taking the best (min) value.
    grp = sub.groupby([index, columns])[value].min().reset_index()
    return grp.pivot(index=index, columns=columns, values=value)


def fig_tuning_heatmaps(tuning_dir, fig_dir, prefix="tuning"):
    cfg_df = _load(os.path.join(tuning_dir, f"{prefix}_config_summary.csv"))
    if cfg_df is None:
        _warn("tuning_config_summary.csv missing; skipping heatmaps.")
        return
    value = "median_integrated_l2_F"

    # gamma vs eta
    piv = _pivot_heatmap(cfg_df, "eta", "gamma", value)
    path = os.path.join(fig_dir, "fig_tuning_gamma_eta_heatmap.png")
    if piv is None or piv.size == 0:
        _warn("no estimated-target configs for gamma-eta heatmap.")
    else:
        fig, ax = plt.subplots(figsize=(6.0, 4.6))
        if piv.shape[0] < 2 or piv.shape[1] < 2:
            _warn(f"gamma-eta grid is degenerate {piv.shape}; drawing a small map.")
        plotting.heatmap(ax, piv, xlabel=r"$\gamma$", ylabel=r"$\eta$",
                         title=r"Estimated-target FR: median $\int\|\hat F-F_{\rm ref}\|\,dt$",
                         cbar_label=value)
        plotting.save_fig(fig, path)
        print("   wrote", os.path.relpath(path))

    # gamma vs burn-in
    piv = _pivot_heatmap(cfg_df, "burnin_fraction", "gamma", value)
    path = os.path.join(fig_dir, "fig_tuning_gamma_burnin_heatmap.png")
    if piv is None or piv.size == 0:
        _warn("no estimated-target configs for gamma-burnin heatmap.")
    else:
        fig, ax = plt.subplots(figsize=(6.0, 4.6))
        if piv.shape[0] < 2 or piv.shape[1] < 2:
            _warn(f"gamma-burnin grid is degenerate {piv.shape}; drawing a small map.")
        plotting.heatmap(ax, piv, xlabel=r"$\gamma$", ylabel="burn-in fraction",
                         title=r"Estimated-target FR: median $\int\|\hat F-F_{\rm ref}\|\,dt$",
                         cbar_label=value)
        plotting.save_fig(fig, path)
        print("   wrote", os.path.relpath(path))


def fig_tuning_boxplot(tuning_dir, fig_dir, prefix="tuning", top_k=6):
    cfg_df = _load(os.path.join(tuning_dir, f"{prefix}_config_summary.csv"))
    final_df = _load(os.path.join(tuning_dir, f"{prefix}_final_summary.csv"))
    if cfg_df is None or final_df is None:
        _warn("config/final summary missing; skipping top-config boxplot.")
        return
    ranked = cfg_df.sort_values("median_integrated_l2_F").head(top_k)
    groups, colors = {}, {}
    for _, r in ranked.iterrows():
        cid = r["config_id"]
        vals = final_df[final_df["config_id"] == cid]["final_l2_F"].values
        label = _short_label(r)
        groups[label] = vals
        colors[label] = plotting.method_color(r["method"])
    if not groups:
        _warn("no configs to plot in boxplot.")
        return
    fig, ax = plt.subplots(figsize=(max(7.0, 1.3 * len(groups)), 4.8))
    plotting.boxplot_by_group(ax, groups, ylabel=r"final $\|\hat F_T-F_{\rm ref}\|_{L^2}$",
                              title=f"Top {len(groups)} configs by median integrated $L^2(F)$"
                                    " (over seeds)", colors=colors)
    path = os.path.join(fig_dir, "fig_tuning_top_configs_boxplot.png")
    plotting.save_fig(fig, path)
    print("   wrote", os.path.relpath(path))


def _short_label(r):
    if r["method"] == "abf_only":
        return "ABF only"
    tt = {"estimated": "est", "uniform": "unif", "oracle": "orac",
          "self": "self"}.get(r["target_type"], r["target_type"])
    return f"{tt} g={r['gamma']:g}\n eta={r['eta']:g} bi={r['burnin_fraction']:g}"


# --------------------------------------------------------------------------- #
# Eval figures
# --------------------------------------------------------------------------- #
def fig01_reference_geometry(ref_dir, fig_dir):
    prof = _load(os.path.join(ref_dir, "reference_profile.csv"))
    grid_path = os.path.join(ref_dir, "reference_grid.npz")
    if prof is None or not os.path.exists(grid_path):
        _warn("reference files missing; skipping fig01 (run run_reference_2d.py).")
        return
    g = np.load(grid_path)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    cf = plotting.potential_contour(axes[0, 0], g["x_grid"], g["y_grid"], g["V_grid"])
    fig.colorbar(cf, ax=axes[0, 0])
    cf = plotting.density_contour(axes[0, 1], g["x_grid"], g["y_grid"], g["rho_grid"])
    fig.colorbar(cf, ax=axes[0, 1])
    plotting.reference_profile(axes[1, 0], prof["x"], prof["F_ref"],
                               r"$F_{\rm ref}(x)$", r"Reference free energy")
    plotting.reference_profile(axes[1, 1], prof["x"], prof["Fprime_ref"],
                               r"$F'_{\rm ref}(x)$", r"Reference mean force")
    path = os.path.join(fig_dir, "fig01_reference_geometry.png")
    plotting.save_fig(fig, path)
    print("   wrote", os.path.relpath(path))


def _method_order(df):
    order = ["abf_only", "abf_fr_estimated", "abf_fr_uniform", "abf_fr_oracle",
             "abf_fr_self"]
    return [m for m in order if m in set(df["method"])]


def _center_on_window(x, y, lo=-2.5, hi=2.5):
    """Subtract the mean over the evaluation window [lo, hi] (display centring).

    The L2 metric centres free-energy profiles by their mean over this same
    window, so centring the plotted curves the same way makes the visual
    comparison consistent with the reported errors.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = (x >= lo) & (x <= hi)
    return y - (np.mean(y[m]) if np.any(m) else np.mean(y))


def fig02_time_curves(eval_dir, fig_dir, prefix="eval"):
    long = _load(os.path.join(eval_dir, f"{prefix}_runs_long.csv"))
    if long is None:
        _warn("eval_runs_long.csv missing; skipping fig02.")
        return
    panels = [("l2_Fprime", r"$\|\hat F'-F'_{\rm ref}\|_{L^2}$", True),
              ("l2_F", r"$\|\hat F-F_{\rm ref}\|_{L^2}$", True),
              ("barrier_crossings", "cumulative barrier crossings", False),
              ("fr_event_fraction", "FR event fraction", False)]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    for ax, (col, ylab, logy) in zip(axes.ravel(), panels):
        for method in _method_order(long):
            sub = long[long["method"] == method]
            agg = sub.groupby("t")[col].mean()
            plotting.time_curve(ax, agg.index.values, agg.values,
                                plotting.method_label(method),
                                plotting.method_color(method), logy=logy)
        ax.set_xlabel("time"); ax.set_ylabel(ylab); ax.set_title(ylab)
        ax.legend(fontsize=8)
    path = os.path.join(fig_dir, "fig02_eval_time_curves.png")
    plotting.save_fig(fig, path)
    print("   wrote", os.path.relpath(path))


def fig03_final_profiles(eval_dir, ref_dir, fig_dir, prefix="eval"):
    prof = _load(os.path.join(eval_dir, f"{prefix}_profiles.csv"))
    ref = _load(os.path.join(ref_dir, "reference_profile.csv"))
    if prof is None or ref is None:
        _warn("eval_profiles.csv or reference_profile.csv missing; skipping fig03.")
        return
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    methods = _method_order(prof)
    # Average profiles over seeds for each method.
    for method in methods:
        sub = prof[prof["method"] == method]
        gp = sub.groupby("x")
        x = np.array(sorted(sub["x"].unique()))
        Fp = gp["Fprime_hat"].mean().reindex(x).values
        F = _center_on_window(x, gp["F_hat"].mean().reindex(x).values)
        p = gp["p_hat"].mean().reindex(x).values
        col = plotting.method_color(method)
        axes[0].plot(x, Fp, color=col, label=plotting.method_label(method))
        axes[1].plot(x, F, color=col, label=plotting.method_label(method))
        axes[2].plot(x, p, color=col, label=plotting.method_label(method))
    axes[0].plot(ref["x"], ref["Fprime_ref"], "k--", lw=1.5, label="reference")
    # Centre F_ref the same way as the estimates (mean over the eval window).
    Fref_c = _center_on_window(ref["x"].values, ref["F_ref"].values)
    axes[1].plot(ref["x"], Fref_c, "k--", lw=1.5, label="reference")
    axes[2].plot(ref["x"], ref["p_ref"], "k:", lw=1.3, label=r"$p_{\rm ref}$ (unbiased)")
    for ax, title, yl in zip(
            axes, [r"Mean force $\hat F'(x)$", r"Free energy $\hat F(x)$",
                   r"x-marginal $\hat p(x)$"],
            [r"$F'$", r"$F$", "density"]):
        ax.set_xlabel("x"); ax.set_ylabel(yl); ax.set_title(title)
        ax.set_xlim(-2.7, 2.7); ax.legend(fontsize=8)
    path = os.path.join(fig_dir, "fig03_eval_final_profiles.png")
    plotting.save_fig(fig, path)
    print("   wrote", os.path.relpath(path))


def fig04_target_ablation(eval_dir, fig_dir, prefix="eval"):
    final = _load(os.path.join(eval_dir, f"{prefix}_final_summary.csv"))
    if final is None:
        _warn("eval_final_summary.csv missing; skipping fig04.")
        return
    metrics_cols = [("final_l2_F", r"final $\|\hat F-F_{\rm ref}\|$"),
                    ("integrated_l2_F", r"$\int\|\hat F-F_{\rm ref}\|\,dt$"),
                    ("final_marginal_l2_uniform", r"marginal $\|\hat p-{\rm Unif}\|$"),
                    ("barrier_crossings", "barrier crossings")]
    methods = _method_order(final)
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.4))
    for ax, (col, title) in zip(axes, metrics_cols):
        groups = {plotting.method_label(m): final[final["method"] == m][col].values
                  for m in methods}
        colors = {plotting.method_label(m): plotting.method_color(m) for m in methods}
        plotting.boxplot_by_group(ax, groups, ylabel=title, title=title, colors=colors)
    path = os.path.join(fig_dir, "fig04_target_ablation.png")
    plotting.save_fig(fig, path)
    print("   wrote", os.path.relpath(path))


def fig05_good_vs_bad(eval_dir, ref_dir, fig_dir, prefix="eval"):
    final = _load(os.path.join(eval_dir, f"{prefix}_final_summary.csv"))
    long = _load(os.path.join(eval_dir, f"{prefix}_runs_long.csv"))
    prof = _load(os.path.join(eval_dir, f"{prefix}_profiles.csv"))
    cond = _load(os.path.join(eval_dir, f"{prefix}_conditional_diagnostics.csv"))
    ref = _load(os.path.join(ref_dir, "reference_profile.csv"))
    if final is None or long is None or prof is None:
        _warn("eval summary/long/profiles missing; skipping fig05.")
        return
    # Identify good vs bad FR configs by median final L2(F) among FR methods.
    fr = final[final["method"] != "abf_only"]
    if fr.empty:
        _warn("no FR runs in eval; skipping fig05.")
        return
    med = fr.groupby("config_id")["final_l2_F"].median().sort_values()
    good_cid = med.index[0]
    bad_cid = med.index[-1]
    abf_cid = final[final["method"] == "abf_only"]["config_id"].iloc[0] \
        if (final["method"] == "abf_only").any() else None
    sel = {"good FR": good_cid, "bad FR": bad_cid}
    if abf_cid is not None:
        sel["ABF only"] = abf_cid
    colmap = {"good FR": "#2ca02c", "bad FR": "#d62728", "ABF only": "#222222"}

    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))
    # (a) L2(F) over time
    ax = axes[0, 0]
    for lab, cid in sel.items():
        sub = long[long["config_id"] == cid]
        agg = sub.groupby("t")["l2_F"].mean()
        ax.plot(agg.index, agg.values, color=colmap[lab], label=lab)
    ax.set_yscale("log"); ax.set_xlabel("time")
    ax.set_ylabel(r"$\|\hat F-F_{\rm ref}\|_{L^2}$")
    ax.set_title("Free-energy error over time"); ax.legend(fontsize=8)
    # (b) final F profile
    ax = axes[0, 1]
    for lab, cid in sel.items():
        sub = prof[prof["config_id"] == cid]
        x = np.array(sorted(sub["x"].unique()))
        F = _center_on_window(x, sub.groupby("x")["F_hat"].mean().reindex(x).values)
        ax.plot(x, F, color=colmap[lab], label=lab)
    if ref is not None:
        ax.plot(ref["x"], _center_on_window(ref["x"].values, ref["F_ref"].values),
                "k--", label="reference")
    ax.set_xlim(-2.7, 2.7); ax.set_xlabel("x"); ax.set_ylabel(r"$\hat F(x)$")
    ax.set_title("Final free-energy profile"); ax.legend(fontsize=8)
    # (c) FR event fraction over time
    ax = axes[1, 0]
    for lab, cid in sel.items():
        sub = long[long["config_id"] == cid]
        agg = sub.groupby("t")["fr_event_fraction"].mean()
        ax.plot(agg.index, agg.values, color=colmap[lab], label=lab)
    ax.set_xlabel("time"); ax.set_ylabel("FR event fraction")
    ax.set_title("Birth--death event fraction"); ax.legend(fontsize=8)
    # (d) conditional Y|X L2 by bin
    ax = axes[1, 1]
    if cond is not None:
        _conditional_bars(ax, cond, sel)
    else:
        ax.text(0.5, 0.5, "conditional diagnostics missing",
                ha="center", va="center"); ax.set_axis_off()
    path = os.path.join(fig_dir, "fig05_good_vs_bad_fr.png")
    plotting.save_fig(fig, path)
    print("   wrote", os.path.relpath(path))


def _conditional_bars(ax, cond, sel):
    """Grouped bars of conditional Y|X L2 by x-bin for selected configs.

    The conditional CSV identifies runs by (method, target_type, seed) rather
    than config_id, so selected configs are matched back via their method.
    """
    t_final = cond["t"].max()
    cond = cond[cond["t"] == t_final]
    bins = sorted(cond["x_bin_center"].unique())
    width = 0.8 / max(len(sel), 1)
    for i, (lab, _cid) in enumerate(sel.items()):
        method = "abf_only" if lab == "ABF only" else _label_to_method(lab, cond)
        sub = cond[cond["method"] == method]
        vals = [np.nanmean(sub[sub["x_bin_center"] == b]["conditional_l2_y"])
                for b in bins]
        xs = np.arange(len(bins)) + i * width
        ax.bar(xs, vals, width=width, label=lab,
               color={"good FR": "#2ca02c", "bad FR": "#d62728",
                      "ABF only": "#222222"}.get(lab, "#888"))
    ax.set_xticks(np.arange(len(bins)) + 0.4 - width / 2)
    ax.set_xticklabels([f"{b:g}" for b in bins])
    ax.set_xlabel("x-bin centre")
    ax.set_ylabel(r"$\|\hat p(y|x)-p_{\rm ref}(y|x)\|_{L^2}$")
    ax.set_title(r"Conditional $Y\,|\,X$ error"); ax.legend(fontsize=8)


def _label_to_method(lab, cond):
    # Fall back to any FR method present (good/bad are both FR).
    fr_methods = [m for m in cond["method"].unique() if m != "abf_only"]
    return fr_methods[0] if fr_methods else "abf_only"


def fig06_conditional(eval_dir, fig_dir, prefix="eval"):
    cond = _load(os.path.join(eval_dir, f"{prefix}_conditional_diagnostics.csv"))
    if cond is None:
        _warn("eval_conditional_diagnostics.csv missing; skipping fig06.")
        return
    t_final = cond["t"].max()
    cond = cond[cond["t"] == t_final]
    bins = sorted(cond["x_bin_center"].unique())
    methods = [m for m in ["abf_only", "abf_fr_estimated", "abf_fr_uniform",
                           "abf_fr_oracle"] if m in set(cond["method"])]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    width = 0.8 / max(len(methods), 1)
    for i, m in enumerate(methods):
        sub = cond[cond["method"] == m]
        l2 = [np.nanmean(sub[sub["x_bin_center"] == b]["conditional_l2_y"]) for b in bins]
        kl = [np.nanmean(sub[sub["x_bin_center"] == b]["conditional_kl_y"]) for b in bins]
        xs = np.arange(len(bins)) + i * width
        axes[0].bar(xs, l2, width=width, label=plotting.method_label(m),
                    color=plotting.method_color(m))
        axes[1].bar(xs, kl, width=width, label=plotting.method_label(m),
                    color=plotting.method_color(m))
    for ax, ylab, title in zip(
            axes, [r"$\|\hat p(y|x)-p_{\rm ref}\|_{L^2}$", r"$D_{\rm KL}(\hat p\,\|\,p_{\rm ref})$"],
            ["Conditional Y|X  L2 error", "Conditional Y|X  KL divergence"]):
        ax.set_xticks(np.arange(len(bins)) + 0.4 - width / 2)
        ax.set_xticklabels([f"{b:g}" for b in bins])
        ax.set_xlabel("x-bin centre"); ax.set_ylabel(ylab); ax.set_title(title)
        ax.legend(fontsize=8)
    path = os.path.join(fig_dir, "fig06_conditional_y_given_x.png")
    plotting.save_fig(fig, path)
    print("   wrote", os.path.relpath(path))


def main(argv=None):
    args = parse_args(argv)
    root = args.output_root
    if args.config:
        cfg = io_utils.load_config(args.config)
        root = cfg.get("output_root", root)
    plotting.set_style()
    ref_dir = os.path.join(root, "reference")

    # Resolve the stage's data dir, figure dir and CSV prefix.  GPU stages get
    # BOTH the tuning-style (heatmaps/boxplot) and eval-style (time-curves,
    # profiles, ablation, conditional) figures since they have all the data.
    sub = io_utils.STAGE_TO_DIR.get(args.stage, args.stage)
    prefix = io_utils.stage_prefix(args.stage)
    data_dir = os.path.join(root, sub)
    fig_dir = io_utils.ensure_dir(os.path.join(
        root, io_utils.STAGE_TO_FIGDIR.get(args.stage, f"figures_{sub}")))
    is_gpu = args.stage in ("smoke_gpu", "tuning_gpu", "production_gpu")
    print(f"[plot_abf_fr_study] stage={args.stage}  in={data_dir}  out={fig_dir} "
          f"prefix={prefix}")

    if args.stage in ("tuning",) or is_gpu:
        fig_tuning_heatmaps(data_dir, fig_dir, prefix=prefix)
        fig_tuning_boxplot(data_dir, fig_dir, prefix=prefix)
    if args.stage in ("eval",) or is_gpu:
        fig01_reference_geometry(ref_dir, fig_dir)
        fig02_time_curves(data_dir, fig_dir, prefix=prefix)
        fig03_final_profiles(data_dir, ref_dir, fig_dir, prefix=prefix)
        fig04_target_ablation(data_dir, fig_dir, prefix=prefix)
        fig05_good_vs_bad(data_dir, ref_dir, fig_dir, prefix=prefix)
        fig06_conditional(data_dir, fig_dir, prefix=prefix)
    print("[plot_abf_fr_study] done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
