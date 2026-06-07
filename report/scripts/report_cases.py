#!/usr/bin/env python3
r"""Case II (entropic bottleneck) and Case III (WCA dimer) report assets.

Self-contained so that both ``build_report_assets.py`` and
``check_report_numbers.py`` can import it without triggering the metastability
build. It renders the EB/WCA figures in the report's matplotlib style and
returns CSV/npz-derived LaTeX macros + a JSON block, so no prose number can
drift from the data. Unlike the metastability pipeline (which recomputes from
per-run CSVs), the EB/WCA numbers are read from the already-aggregated
``config_summary.csv`` / ``winrates.csv`` medians.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
EB_DIR = os.path.join(REPO_ROOT, "results", "entropic_bottleneck", "summaries")
WCA_DIR = os.path.join(REPO_ROOT, "results", "wca_production", "summaries")

CASE_COLOR = {"abf": "#222222", "fr_estimated": "#1f77b4",
              "fr_uniform": "#ff7f0e", "fr_oracle": "#2ca02c"}
CASE_LABEL = {"abf": "ABF only", "fr_estimated": "ABF+FR (estimated)",
              "fr_uniform": "ABF+FR (uniform)",
              "fr_oracle": "ABF+FR (oracle, diagnostic)"}


def _warn(msg):
    print(f"[report_cases] WARNING: {msg}")


def _set_style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 160, "font.size": 11,
        "axes.titlesize": 12, "axes.labelsize": 11, "legend.fontsize": 9,
        "axes.grid": True, "grid.alpha": 0.25, "lines.linewidth": 1.8,
    })


def _save(fig, fig_dir, name):
    fig.tight_layout()
    path = os.path.join(fig_dir, name)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"   [fig] wrote figures/{name}")
    return path


# --------------------------------------------------------------------------- #
# Entropic-bottleneck loaders
# --------------------------------------------------------------------------- #
def load_eb():
    """Return (config_summary df, arrays dict-of-dicts keyed by 'stage|cfg')."""
    cs = pd.read_csv(os.path.join(EB_DIR, "config_summary.csv"))
    z = np.load(os.path.join(EB_DIR, "arrays.npz"), allow_pickle=True)
    arr = {}
    for key in z.files:
        head, field = key.split("::")
        arr.setdefault(head, {})[field] = z[key]
    return cs, arr


def _eb_pick(arr, prefix):
    """First arrays key whose head starts with ``prefix`` (stage|method|...)."""
    for head in arr:
        if head.startswith(prefix):
            return arr[head]
    return None


# --------------------------------------------------------------------------- #
# Entropic-bottleneck figures
# --------------------------------------------------------------------------- #
def _eb_stage_cfg(arr, stage):
    """ABF and FR config dicts for a stage (prefers stage1_seeds, falls back)."""
    abf = _eb_pick(arr, f"{stage}|abf|")
    fr = _eb_pick(arr, f"{stage}|fr_estimated|")
    return abf, fr


def fig_eb_01_convergence(arr, fig_dir):
    abf, fr = _eb_stage_cfg(arr, "stage1_seeds")
    if abf is None:
        abf, fr = _eb_stage_cfg(arr, "stage0_reproduce")
    if abf is None or fr is None:
        _warn("EB convergence arrays missing; skipping fig_eb_01.")
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7))
    for ax, field, ylab in zip(
            axes, ["l2_f_t", "l2_fp_t"],
            [r"$\|\hat F_t-F_{\rm ref}\|_{L^2}$",
             r"$\|\hat F'_t-F'_{\rm ref}\|_{L^2}$"]):
        for d, m in [(abf, "abf"), (fr, "fr_estimated")]:
            t = d["t"]
            y = d[field]                       # (seeds, T)
            q50 = np.nanmedian(y, axis=0)
            q25 = np.nanpercentile(y, 25, axis=0)
            q75 = np.nanpercentile(y, 75, axis=0)
            ax.plot(t, q50, color=CASE_COLOR[m], label=CASE_LABEL[m])
            ax.fill_between(t, q25, q75, color=CASE_COLOR[m], alpha=0.18, lw=0)
        ax.set_yscale("log"); ax.set_xlabel("time $t$"); ax.set_ylabel(ylab)
        ax.set_title(ylab); ax.legend(fontsize=8)
    fig.suptitle("Entropic bottleneck: convergence (median and IQR over seeds)",
                 y=1.02)
    _save(fig, fig_dir, "fig_eb_01_convergence.png")


def fig_eb_02_omega(cs, fig_dir):
    s = cs[(cs["stage"] == "stage2_omega") & (cs["method"] == "fr_estimated")]
    if s.empty:
        _warn("EB omega sweep missing; skipping fig_eb_02.")
        return
    s = s.sort_values("omega_in")
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.plot(s["omega_in"], 100 * s["gain_l2_f_vs_abf"], "o-",
            color=CASE_COLOR["fr_estimated"])
    ax.axhline(0, color="k", ls="--", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel(r"bottleneck stiffness $\omega_{\rm in}$ (log scale)")
    ax.set_ylabel(r"FR gain over ABF, $\%$ of $\|F-F_{\rm ref}\|$")
    ax.set_title(r"Gain is large but flat in $\omega_{\rm in}$ "
                 r"(energetic barrier dominates)")
    _save(fig, fig_dir, "fig_eb_02_omega.png")


def fig_eb_03_gamma(cs, fig_dir):
    s = cs[(cs["stage"] == "stage4_gamma") & (cs["method"] == "fr_estimated")]
    if s.empty:
        _warn("EB gamma sweep missing; skipping fig_eb_03.")
        return
    s = s.sort_values("gamma")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    axes[0].plot(s["gamma"], 100 * s["gain_l2_f_vs_abf"], "o-",
                 color=CASE_COLOR["fr_estimated"])
    axes[0].axhline(0, color="k", ls="--", lw=0.8)
    axes[0].set_xlabel(r"FR rate $\gamma$"); axes[0].set_ylabel(r"gain over ABF, $\%$")
    axes[0].set_title("(a) Gain rises monotonically (no failure boundary)")
    axes[1].plot(s["gamma"], s["med_final_ess"], "s-", color="#9467bd",
                 label="ancestor ESS")
    axes[1].set_xlabel(r"FR rate $\gamma$"); axes[1].set_ylabel("ancestor ESS")
    ax2 = axes[1].twinx()
    ax2.plot(s["gamma"], s["med_repl_fraction"], "^-", color="#8c564b",
             label="repl. fraction")
    ax2.set_ylabel("replacement fraction"); ax2.grid(False)
    axes[1].set_title("(b) Diversity falls but accuracy keeps improving")
    lines = axes[1].get_lines() + ax2.get_lines()
    axes[1].legend(lines, [l.get_label() for l in lines], fontsize=8)
    _save(fig, fig_dir, "fig_eb_03_gamma.png")


def fig_eb_04_conditional(arr, fig_dir):
    abf, fr = _eb_stage_cfg(arr, "stage1_seeds")
    if abf is None:
        abf, fr = _eb_stage_cfg(arr, "stage0_reproduce")
    if abf is None or "cond_centers" not in abf:
        _warn("EB conditional arrays missing; skipping fig_eb_04.")
        return
    centers = abf["cond_centers"]
    ref_var = abf["cond_ref_var"]
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.plot(centers, ref_var, "k--", lw=1.6,
            label=r"analytic $1/(\beta\omega(x)^2)$")
    for d, m in [(abf, "abf"), (fr, "fr_estimated")]:
        emp = np.nanmedian(d["cond_emp_var"], axis=0)  # (seeds, bins) -> bins
        ax.plot(centers, emp, "o-", color=CASE_COLOR[m], label=CASE_LABEL[m])
    ax.set_yscale("log"); ax.set_xlabel("$x$-bin centre")
    ax.set_ylabel(r"$\mathrm{Var}(Y\mid X=x)$")
    ax.set_title("Conditional fidelity vs the analytic law")
    ax.legend(fontsize=8)
    _save(fig, fig_dir, "fig_eb_04_conditional.png")


# --------------------------------------------------------------------------- #
# WCA loaders + figures
# --------------------------------------------------------------------------- #
def load_wca():
    cs = pd.read_csv(os.path.join(WCA_DIR, "config_summary.csv"))
    wr = pd.read_csv(os.path.join(WCA_DIR, "winrates.csv"))
    ts = pd.read_csv(os.path.join(WCA_DIR, "timeseries_summary.csv"))
    prof = np.load(os.path.join(WCA_DIR, "profiles_summary.npz"), allow_pickle=True)
    pf = {k: prof[k] for k in prof.files}
    return cs, wr, ts, pf


# main-stage methods, in display order
WCA_MAIN = [("abf", "abf"), ("fr_est_tuned", "fr_estimated"),
            ("fr_uniform", "fr_uniform"), ("fr_oracle", "fr_oracle")]


def _wca_method_color(method):
    return CASE_COLOR.get(method, "#888888")


def fig_wca_01_convergence(ts, fig_dir):
    s = ts[ts["stage"] == "main"]
    if s.empty:
        _warn("WCA main timeseries missing; skipping fig_wca_01.")
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7))
    for ax, col, lo, hi, ylab in [
            ("l2_f_median", "l2_f_q25", "l2_f_q75",
             None, r"$\|\hat F_t-F_{\rm ref}\|_{L^2}$"),
            ("l2_fp_median", "l2_fp_q25", "l2_fp_q75",
             None, r"$\|\hat F'_t-F'_{\rm ref}\|_{L^2}$")]:
        # col/lo/hi packed oddly; fix mapping below
        pass
    # clean re-do with explicit triples
    panels = [("l2_f", r"$\|\hat F_t-F_{\rm ref}\|_{L^2}$"),
              ("l2_fp", r"$\|\hat F'_t-F'_{\rm ref}\|_{L^2}$")]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7))
    for ax, (base, ylab) in zip(axes, panels):
        for name, method in WCA_MAIN:
            sub = s[s["name"] == name].sort_values("t")
            if sub.empty:
                continue
            ax.plot(sub["t"], sub[f"{base}_median"], color=_wca_method_color(method),
                    label=CASE_LABEL.get(method, name))
            ax.fill_between(sub["t"], sub[f"{base}_q25"], sub[f"{base}_q75"],
                            color=_wca_method_color(method), alpha=0.15, lw=0)
        ax.set_yscale("log"); ax.set_xlabel("step"); ax.set_ylabel(ylab)
        ax.set_title(ylab); ax.legend(fontsize=8)
    fig.suptitle("WCA dimer: convergence (median and IQR over seeds)", y=1.02)
    _save(fig, fig_dir, "fig_wca_01_convergence.png")


def _wca_prof(pf, stage, name, method, field, nsteps=250000, nrep=1024, a=1.5):
    key = f"{stage}|{name}|{method}|{nsteps}|{nrep}|{a}@{field}"
    return pf.get(key)


def fig_wca_02_profiles(pf, fig_dir):
    grid = pf.get("grid"); Fref = pf.get("ref_free_energy"); Fpref = pf.get("ref_mean_force")
    if grid is None:
        _warn("WCA profiles npz missing grid; skipping fig_wca_02.")
        return
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    for name, method in WCA_MAIN:
        F = _wca_prof(pf, "main", name, method, "F")
        Fp = _wca_prof(pf, "main", name, method, "Fp")
        p = _wca_prof(pf, "main", name, method, "p")
        if F is None:
            continue
        c = _wca_method_color(method)
        axes[0].plot(grid, F, color=c, label=CASE_LABEL.get(method, name))
        axes[1].plot(grid, Fp, color=c)
        if p is not None:
            axes[2].plot(grid, p, color=c)
    axes[0].plot(grid, Fref, "k--", lw=1.4, label="TI reference")
    axes[1].plot(grid, Fpref, "k--", lw=1.4)
    for ax, title, yl in zip(
            axes, [r"(a) Free energy $\hat F(z)$", r"(b) Mean force $\hat F'(z)$",
                   r"(c) $z$-marginal"], [r"$F$", r"$F'$", "density"]):
        ax.set_xlabel("$z$"); ax.set_ylabel(yl); ax.set_title(title)
    axes[0].legend(fontsize=8)
    fig.suptitle("WCA final profiles (median over seeds) vs TI reference", y=1.02)
    _save(fig, fig_dir, "fig_wca_02_profiles.png")


def fig_wca_03_seed(cs, fig_dir):
    s = cs[cs["stage"] == "main"].copy()
    if s.empty:
        _warn("WCA main config_summary missing; skipping fig_wca_03.")
        return
    order = ["abf", "fr_est_gentle", "fr_est_tuned", "fr_uniform", "fr_oracle",
             "fr_est_strong", "fr_est_aggressive"]
    s["ord"] = s["name"].map({n: i for i, n in enumerate(order)})
    s = s.sort_values("ord")
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    for ax, (col, lo, hi, title) in zip(axes, [
            ("l2_f_median", "l2_f_q25", "l2_f_q75", r"final $\|\hat F-F_{\rm ref}\|$"),
            ("integrated_l2_f_median", "integrated_l2_f_q25", "integrated_l2_f_q75",
             r"integrated $\|\hat F-F_{\rm ref}\|$"),
            ("l2_fp_median", "l2_fp_q25", "l2_fp_q75",
             r"final $\|\hat F'-F'_{\rm ref}\|$")]):
        xs = np.arange(len(s))
        yerr = np.vstack([(s[col] - s[lo]).values, (s[hi] - s[col]).values])
        colors = ["#222222" if n == "abf" else "#1f77b4" for n in s["name"]]
        ax.bar(xs, s[col], yerr=yerr, color=colors, alpha=0.8, capsize=3)
        ax.set_xticks(xs)
        ax.set_xticklabels(s["name"], rotation=30, ha="right", fontsize=8)
        ax.set_title(title)
    fig.suptitle("WCA per-configuration error (median and IQR over seeds)", y=1.04)
    _save(fig, fig_dir, "fig_wca_03_seed.png")


def fig_wca_04_mechanism(pf, fig_dir):
    grid = pf.get("grid")
    if grid is None:
        _warn("WCA profiles npz missing; skipping fig_wca_04.")
        return
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    for name, method in [("abf", "abf"), ("fr_est_tuned", "fr_estimated")]:
        neff = _wca_prof(pf, "main", name, method, "neff")
        Fp = _wca_prof(pf, "main", name, method, "Fp")
        if neff is not None:
            axes[0].plot(grid, neff, color=_wca_method_color(method),
                         label=CASE_LABEL.get(method, name))
        if Fp is not None:
            axes[1].plot(grid, np.abs(Fp - pf["ref_mean_force"]),
                         color=_wca_method_color(method))
    birth = _wca_prof(pf, "main", "fr_est_tuned", "fr_estimated", "birth")
    death = _wca_prof(pf, "main", "fr_est_tuned", "fr_estimated", "death")
    if birth is not None:
        axes[2].plot(grid, birth, color="#2ca02c", label="clone (birth)")
        axes[2].plot(grid, death, color="#d62728", label="death")
        axes[2].legend(fontsize=8)
    axes[0].set_title(r"(a) Effective ABF samples $N_{\rm eff}(z)$"); axes[0].legend(fontsize=8)
    axes[1].set_title(r"(b) Mean-force error $|\hat F'-F'_{\rm ref}|$")
    axes[2].set_title("(c) Birth / death locations")
    for ax in axes:
        ax.set_xlabel("$z$")
    fig.suptitle("WCA mechanism: FR raises sampling quality, not the marginal", y=1.02)
    _save(fig, fig_dir, "fig_wca_04_mechanism.png")


def fig_wca_05_failure(cs, fig_dir):
    s = cs[cs["stage"] == "failure"].copy()
    if s.empty:
        _warn("WCA failure stage missing; skipping fig_wca_05.")
        return
    abf = s[s["name"] == "abf"]
    fr = s[s["name"].str.startswith("sweep_fr_rate_")].copy()
    fr["rate"] = fr["name"].str.replace("sweep_fr_rate_", "").astype(float)
    fr = fr.sort_values("rate")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    axes[0].plot(fr["rate"], fr["l2_f_median"], "o-", color="#1f77b4")
    if not abf.empty:
        axes[0].axhline(abf["l2_f_median"].iloc[0], color="#222222", ls="--",
                        label="ABF only")
        axes[0].legend(fontsize=8)
    axes[0].set_xscale("log"); axes[0].set_xlabel(r"\texttt{fr\_rate}")
    axes[0].set_ylabel(r"final $\|\hat F-F_{\rm ref}\|$")
    axes[0].set_title("(a) Error reverses above the sweet spot")
    axes[1].plot(fr["rate"], fr["final_ancestor_ess_median"], "s-", color="#9467bd")
    axes[1].set_xscale("log"); axes[1].set_xlabel(r"\texttt{fr\_rate}")
    axes[1].set_ylabel("ancestor ESS")
    axes[1].set_title("(b) Diversity collapses as the rate grows")
    _save(fig, fig_dir, "fig_wca_05_failure.png")


def fig_wca_06_difficulty(wr, fig_dir):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    b = wr[wr["stage"] == "difficulty_budget"].sort_values("n_steps")
    axes[0].plot(b["n_steps"], b["median_gain_pct_F"], "o-", color="#1f77b4")
    axes[0].set_xlabel("steps"); axes[0].set_ylabel(r"gain over ABF, $\%$")
    axes[0].set_title("(a) Gain grows with budget")
    r = wr[wr["stage"] == "difficulty_replicas"].sort_values("n_replicas")
    axes[1].plot(r["n_replicas"], r["median_gain_pct_F"], "s-", color="#2ca02c")
    axes[1].set_xscale("log", base=2); axes[1].set_xlabel("replicas $N$")
    axes[1].set_title("(b) Flat in replica count")
    c = wr[wr["stage"] == "difficulty_crowding"].sort_values("a")
    axes[2].plot(c["a"], c["median_gain_pct_F"], "^-", color="#d62728")
    axes[2].set_xlabel("lattice spacing $a$ (smaller = denser)")
    axes[2].set_title("(c) Largest gain in dense solvent")
    for ax in axes:
        ax.axhline(0, color="k", ls="--", lw=0.7)
    fig.suptitle("WCA difficulty dependence (matched-seed gain over ABF)", y=1.02)
    _save(fig, fig_dir, "fig_wca_06_difficulty.png")


# --------------------------------------------------------------------------- #
# Numbers (macros + JSON), recomputed from config_summary / winrates medians
# --------------------------------------------------------------------------- #
def _row(df, **kw):
    """Single matching row as a Series; raises if not exactly one."""
    m = df
    for k, v in kw.items():
        m = m[m[k] == v]
    if len(m) != 1:
        raise ValueError(f"expected 1 row for {kw}, got {len(m)}")
    return m.iloc[0]


def eb_numbers(cs):
    """Headline EB numbers from config_summary medians."""
    s1_abf = _row(cs, stage="stage1_seeds", method="abf")
    s1_fr = _row(cs, stage="stage1_seeds", method="fr_estimated")
    s0_uni = _row(cs, stage="stage0_reproduce", method="fr_uniform")
    s0_ora = _row(cs, stage="stage0_reproduce", method="fr_oracle")
    b4 = _row(cs, stage="stage3_beta", method="fr_estimated", beta=4.0)
    g = cs[(cs["stage"] == "stage4_gamma") & (cs["method"] == "fr_estimated")]
    g_max = g.sort_values("gamma").iloc[-1]
    g_lo = g.sort_values("gamma").iloc[0]
    return dict(
        abf_l2f=float(s1_abf["med_l2_f"]),
        fr_l2f=float(s1_fr["med_l2_f"]),
        gain_strong_pct=100 * float(s1_fr["gain_l2_f_vs_abf"]),
        cold_win=int(round(s1_fr["winrate_vs_abf"] * s1_fr["n_seeds"])),
        cold_nseeds=int(s1_fr["n_seeds"]),
        uniform_gain_pct=100 * float(s0_uni["gain_l2_f_vs_abf"]),
        oracle_gain_pct=100 * float(s0_ora["gain_l2_f_vs_abf"]),
        warm_gain_beta4_pct=100 * float(b4["gain_l2_f_vs_abf"]),
        gamma_max=float(g_max["gamma"]),
        gamma_max_gain_pct=100 * float(g_max["gain_l2_f_vs_abf"]),
        ess_hi=float(g_lo["med_final_ess"]),
        ess_lo=float(g_max["med_final_ess"]),
    )


def wca_numbers(cs, wr):
    """Headline WCA numbers from config_summary + winrates medians."""
    abf = _row(cs, stage="main", name="abf", method="abf")
    tuned = _row(cs, stage="main", name="fr_est_tuned", method="fr_estimated")
    wr_tuned = _row(wr, stage="main", name="fr_est_tuned")
    wr_uni = _row(wr, stage="main", name="fr_uniform")
    wr_ora = _row(wr, stage="main", name="fr_oracle")
    wr_aggr = _row(wr, stage="main", name="fr_est_aggressive")
    bud = wr[wr["stage"] == "difficulty_budget"].sort_values("n_steps")
    crow_abf = _row(cs, stage="difficulty_crowding", name="abf", method="abf", a=1.35)
    crow_fr = _row(cs, stage="difficulty_crowding", name="fr_est_tuned",
                   method="fr_estimated", a=1.35)
    return dict(
        abf_l2f=float(abf["l2_f_median"]),
        tuned_l2f=float(tuned["l2_f_median"]),
        tuned_gain_pct=float(wr_tuned["median_gain_pct_F"]),
        tuned_win=int(wr_tuned["n_wins"]),
        tuned_nseeds=int(wr_tuned["n_pairs"]),
        uniform_gain_pct=float(wr_uni["median_gain_pct_F"]),
        oracle_gain_pct=float(wr_ora["median_gain_pct_F"]),
        aggressive_gain_pct=float(wr_aggr["median_gain_pct_F"]),
        budget_lo_gain_pct=float(bud.iloc[0]["median_gain_pct_F"]),
        budget_hi_gain_pct=float(bud.iloc[-1]["median_gain_pct_F"]),
        dense_abf_l2f=float(crow_abf["l2_f_median"]),
        dense_fr_l2f=float(crow_fr["l2_f_median"]),
    )


# --------------------------------------------------------------------------- #
# LaTeX tables
# --------------------------------------------------------------------------- #
def _simple_table(rows, header, caption, label, align, note=None):
    lines = [r"\begin{table}[t]", r"\centering", r"\small",
             rf"\caption{{{caption}}}", rf"\label{{{label}}}",
             rf"\begin{{tabular}}{{{align}}}", r"\toprule",
             " & ".join(header) + r" \\", r"\midrule"]
    for r in rows:
        lines.append(" & ".join(r) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    if note:
        lines.append(rf"\par\smallskip\footnotesize {note}")
    lines += [r"\end{table}", ""]
    return "\n".join(lines)


def write_eb_tables(cs, tab_dir):
    # Design table
    s1 = _row(cs, stage="stage1_seeds", method="fr_estimated")
    design = [
        (r"Inverse temperature $\beta$", "8"),
        (r"Barrier height $H$", "2.5"),
        (r"Channel stiffness $\omega_{\rm out},\omega_{\rm in}$", "1, 25"),
        (r"Channel width $s$", "0.25"),
        (r"Walkers $N$", "256"),
        (r"Time step $\Delta t$", r"$10^{-3}$"),
        (r"Steps per run", "40{,}000"),
        (r"FR rate $\gamma$", f"{s1['gamma']:g}"),
        (r"FR stride \texttt{fr\_every}", "10"),
        (r"Score clip $S_{\max}$", "3"),
        (r"Event cap $f_{\max}$", "0.08"),
        (r"Target EMA rate", "0.005"),
    ]
    with open(os.path.join(tab_dir, "eb_design.tex"), "w") as fh:
        fh.write(_simple_table(
            design, [r"Quantity", r"Value"],
            r"Entropic-bottleneck base regime and Fisher--Rao hyperparameters.",
            "tab:eb_design", "ll"))
    # Beta sweep table
    b = cs[(cs["stage"] == "stage3_beta")].copy()
    rows = []
    for beta in sorted(b["beta"].unique()):
        abf = _row(b, beta=beta, method="abf")
        fr = _row(b, beta=beta, method="fr_estimated")
        rows.append([f"{beta:g}", f"{abf['med_l2_f']:.3f}", f"{fr['med_l2_f']:.3f}",
                     f"{100*fr['gain_l2_f_vs_abf']:+.1f}",
                     f"{fr['winrate_vs_abf']:.2f}"])
    with open(os.path.join(tab_dir, "eb_beta.tex"), "w") as fh:
        fh.write(_simple_table(
            rows, [r"$\beta$", r"ABF $\|F\!-\!F_{\rm ref}\|$",
                   r"FR $\|F\!-\!F_{\rm ref}\|$", r"gain (\%)", r"win rate"],
            r"Temperature sweep at $\omega_{\rm in}=25$: the FR gain flips sign "
            r"as ABF crosses from well-sampled (hot) to starved (cold).",
            "tab:eb_beta", "lrrrr"))


def write_wca_tables(cs, wr, tab_dir):
    design = [
        (r"Inverse temperature $\beta$", "1"),
        (r"Time step $\Delta t$", "0.002"),
        (r"Steps per run", "250{,}000"),
        (r"Replicas $N$", "1024"),
        (r"Lattice spacing $a$", "1.5"),
        (r"Dimer barrier $h$, width $w$", "2, 2"),
        (r"FR rate \texttt{fr\_rate}", "0.10"),
        (r"FR onset \texttt{fr\_start\_steps}", "20{,}000"),
        (r"FR stride \texttt{fr\_every}", "5"),
        (r"Score clip $S_{\max}$", "2.0"),
        (r"Event cap $f_{\max}$", "0.02"),
        (r"Target EMA rate", "0.005"),
    ]
    with open(os.path.join(tab_dir, "wca_design.tex"), "w") as fh:
        fh.write(_simple_table(
            design, [r"Quantity", r"Value"],
            r"WCA dimer base system and tuned Fisher--Rao hyperparameters.",
            "tab:wca_design", "ll"))
    # Main comparison table
    order = [("abf", "ABF only", "abf"),
             ("fr_est_tuned", "FR estimated (tuned)", "fr_estimated"),
             ("fr_uniform", "FR uniform", "fr_uniform"),
             ("fr_oracle", r"FR oracle$^{*}$", "fr_oracle"),
             ("fr_est_strong", "FR estimated (strong)", "fr_estimated"),
             ("fr_est_aggressive", "FR estimated (aggressive)", "fr_estimated")]
    rows = []
    for name, label, method in order:
        c = _row(cs, stage="main", name=name, method=method)
        w = wr[(wr["stage"] == "main") & (wr["name"] == name)]
        if name == "abf":
            gain, win = "--", "--"
        else:
            gain = f"{w['median_gain_pct_F'].iloc[0]:+.1f}"
            win = f"{int(w['n_wins'].iloc[0])}/{int(w['n_pairs'].iloc[0])}"
        rows.append([label, f"{c['l2_f_median']:.4f}", f"{c['l2_fp_median']:.4f}",
                     gain, win, f"{c['final_ancestor_ess_median']:.0f}"
                     if np.isfinite(c['final_ancestor_ess_median']) else "--"])
    note = (r"$^{*}$\,The oracle target uses the thermodynamic-integration "
            r"reference $F_{\rm ref}$; it is a diagnostic control, \emph{not} a "
            r"deployable method. Medians over $10$ seeds.")
    with open(os.path.join(tab_dir, "wca_main.tex"), "w") as fh:
        fh.write(_simple_table(
            rows, [r"Method", r"$\|F\!-\!F_{\rm ref}\|$",
                   r"$\|F'\!-\!F'_{\rm ref}\|$", r"gain (\%)", r"wins",
                   r"anc. ESS"],
            r"WCA main comparison (250k steps, $N=1024$, $a=1.5$). The headline "
            r"gain is in the free energy; the mean-force gain is modest.",
            "tab:wca_main", "lrrrrr", note=note))


def write_synthesis_table(ebn, wcan, meta_int_pct, tab_dir):
    rows = [
        ["Metastability", "ABF strong", f"{meta_int_pct:+.1f} (int.)", "--",
         "oracle best"],
        ["Entropic bottleneck (hot, $\\beta{=}4$)", "ABF easy",
         f"{ebn['warm_gain_beta4_pct']:+.1f}", "0/10", "n/a"],
        ["Entropic bottleneck (cold, $\\beta{=}8$)", "ABF starved",
         f"{ebn['gain_strong_pct']:+.1f}",
         f"{ebn['cold_win']}/{ebn['cold_nseeds']}", "oracle worse"],
        ["WCA dimer", "ABF starved", f"{wcan['tuned_gain_pct']:+.1f}",
         f"{wcan['tuned_win']}/{wcan['tuned_nseeds']}", r"oracle $\approx$ est."],
    ]
    with open(os.path.join(tab_dir, "synthesis.tex"), "w") as fh:
        fh.write(_simple_table(
            rows, [r"Case", r"Regime", r"FR gain (\%)", r"wins",
                   r"target effect"],
            r"Cross-case synthesis: the estimated-target FR gain grows with how "
            r"sample-starved ABF is, while the value of the target \emph{shape} "
            r"(oracle vs.\ estimated) shrinks. Metastability gain is on the "
            r"integrated error; the others on final $\|F-F_{\rm ref}\|$.",
            "tab:synthesis", "lllll"))


def case_macros(ebn, wcan):
    def p(v): return f"{v:.1f}"
    def f4(v): return f"{v:.4f}"
    def f3(v): return f"{v:.3f}"
    return {
        "EBabfLtwoF": f3(ebn["abf_l2f"]), "EBfrLtwoF": f3(ebn["fr_l2f"]),
        "EBgainStrong": p(ebn["gain_strong_pct"]),
        "EBcoldWin": str(ebn["cold_win"]),
        "EBuniGain": p(ebn["uniform_gain_pct"]),
        "EBoracleGain": p(ebn["oracle_gain_pct"]),
        "EBwarmGainBetaFour": p(ebn["warm_gain_beta4_pct"]),
        "EBgammaMaxGain": p(ebn["gamma_max_gain_pct"]),
        "EBessHi": f"{ebn['ess_hi']:.0f}", "EBessLo": f"{ebn['ess_lo']:.0f}",
        "WCAabfLtwoF": f4(wcan["abf_l2f"]), "WCAtunedLtwoF": f4(wcan["tuned_l2f"]),
        "WCAtunedGainPct": p(wcan["tuned_gain_pct"]),
        "WCAtunedWin": str(wcan["tuned_win"]),
        "WCAuniGainPct": p(wcan["uniform_gain_pct"]),
        "WCAoracleGainPct": p(wcan["oracle_gain_pct"]),
        "WCAaggrGainPct": p(wcan["aggressive_gain_pct"]),
        "WCAbudgetLoGain": p(wcan["budget_lo_gain_pct"]),
        "WCAbudgetHiGain": p(wcan["budget_hi_gain_pct"]),
        "WCAdenseAbfErr": f3(wcan["dense_abf_l2f"]),
        "WCAdenseFrErr": f3(wcan["dense_fr_l2f"]),
    }


def build_cases(fig_dir, tab_dir, meta_int_pct):
    """Render all EB/WCA figures + tables; return (macros dict, json block)."""
    _set_style()
    eb_cs, eb_arr = load_eb()
    wca_cs, wca_wr, wca_ts, wca_pf = load_wca()

    fig_eb_01_convergence(eb_arr, fig_dir)
    fig_eb_02_omega(eb_cs, fig_dir)
    fig_eb_03_gamma(eb_cs, fig_dir)
    fig_eb_04_conditional(eb_arr, fig_dir)
    fig_wca_01_convergence(wca_ts, fig_dir)
    fig_wca_02_profiles(wca_pf, fig_dir)
    fig_wca_03_seed(wca_cs, fig_dir)
    fig_wca_04_mechanism(wca_pf, fig_dir)
    fig_wca_05_failure(wca_cs, fig_dir)
    fig_wca_06_difficulty(wca_wr, fig_dir)

    ebn = eb_numbers(eb_cs)
    wcan = wca_numbers(wca_cs, wca_wr)
    write_eb_tables(eb_cs, tab_dir)
    write_wca_tables(wca_cs, wca_wr, tab_dir)
    write_synthesis_table(ebn, wcan, meta_int_pct, tab_dir)
    return case_macros(ebn, wcan), {"entropic_bottleneck": ebn, "wca": wcan}





