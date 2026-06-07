#!/usr/bin/env python3
"""Generate the four entropic-bottleneck figures from summaries/arrays.npz.

Fig 1  convergence: L2(F_t), L2(F'_t) vs time, median+IQR over seeds, ABF vs FR
       (stage1_seeds; falls back to stage0_reproduce).
Fig 2  entropic strength sweep: % gain in final L2(F) vs omega_in, + win rate
       (stage2_omega).
Fig 3  FR-rate failure boundary: final L2(F) and ancestor ESS / repl vs gamma
       (stage4_gamma).
Fig 4  conditional fidelity: Var(Y|X in bin) vs analytic 1/(beta omega^2) at
       x in {-1,-0.5,0,0.5,1}, ABF vs FR side by side (stage1_seeds).

Usage: python scripts/plot_entropic_bottleneck.py
"""
from __future__ import annotations

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SUM = os.path.join(os.path.dirname(__file__), "..", "results", "entropic_bottleneck", "summaries")
PLOTS = os.path.join(os.path.dirname(__file__), "..", "results", "entropic_bottleneck", "plots")

C_ABF, C_FR = "C3", "C0"


def load():
    return dict(np.load(os.path.join(SUM, "arrays.npz"), allow_pickle=True))


def g(A, key):
    return A[key] if key in A else None


def med_iqr(stack):
    """stack:(n_seeds, T) -> median, q1, q3 over seeds."""
    q1, q2, q3 = np.percentile(stack, [25, 50, 75], axis=0)
    return q2, q1, q3


def figure1(A):
    """Convergence: prefer stage1_seeds, fall back to stage0_reproduce."""
    stage = "stage1_seeds"
    if f"{stage}|abf|beta8|oin25|gamma15::t" not in A:
        stage = "stage0_reproduce"
    pref = f"{stage}|{{m}}|beta8|oin25|gamma15"
    t = g(A, pref.format(m="abf") + "::t")
    if t is None:
        print("fig1: no data"); return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, key, ttl in [
        (axes[0], "l2_f_t", r"Free-energy error  $\|\hat F_t - F_{\rm ref}\|_{L^2}$"),
        (axes[1], "l2_fp_t", r"Mean-force error  $\|\hat F'_t - F'_{\rm ref}\|_{L^2}$")]:
        for method, c, lbl in [("abf", C_ABF, "ABF only"), ("fr_estimated", C_FR, "ABF + Fisher--Rao")]:
            stk = g(A, pref.format(m=method) + f"::{key}")
            if stk is None:
                continue
            m, q1, q3 = med_iqr(stk)
            ax.fill_between(t, q1, q3, color=c, alpha=0.20)
            ax.plot(t, m, color=c, lw=2.4, label=lbl)
        ax.set(xlabel="time", ylabel="RMS error", yscale="log", title=ttl)
        ax.legend()
    n_seeds = g(A, pref.format(m="abf") + "::seeds")
    fig.suptitle(f"Convergence (median + IQR over {len(n_seeds)} matched seeds, {stage})", y=1.02)
    plt.tight_layout()
    out = os.path.join(PLOTS, "fig1_convergence.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


def figure2(A):
    """Entropic strength sweep: % gain in final L2(F) and win rate vs omega_in."""
    oins = [1.0, 5.0, 10.0, 25.0, 50.0]
    gains, winrates, oks = [], [], []
    for oin in oins:
        pa = f"stage2_omega|abf|beta8|oin{oin:g}|gamma15::l2_f_t"
        pf = f"stage2_omega|fr_estimated|beta8|oin{oin:g}|gamma15::l2_f_t"
        if pa not in A or pf not in A:
            continue
        fa = np.median(A[pa][:, -1]); ff = np.median(A[pf][:, -1])
        gains.append(100.0 * (fa - ff) / fa)
        # matched-seed win rate
        sa = A[f"stage2_omega|abf|beta8|oin{oin:g}|gamma15::seeds"]
        sf = A[f"stage2_omega|fr_estimated|beta8|oin{oin:g}|gamma15::seeds"]
        af = {int(s): A[pa][i, -1] for i, s in enumerate(sa)}
        frd = {int(s): A[pf][i, -1] for i, s in enumerate(sf)}
        shared = sorted(set(af) & set(frd))
        winrates.append(np.mean([frd[s] < af[s] for s in shared]))
        oks.append(oin)
    if not oks:
        print("fig2: no data"); return
    fig, ax1 = plt.subplots(figsize=(7.2, 5.0))
    ax1.plot(oks, gains, "o-", color=C_FR, lw=2.2, ms=8, label="% gain in final $L^2(F)$")
    ax1.axhline(0, color="k", lw=0.8, ls=":")
    ax1.set_xscale("log")
    ax1.set(xlabel=r"$\omega_{\rm in}$ (bottleneck stiffness)",
            ylabel="percentage gain in final $L^2(F)$  (FR vs ABF)")
    ax1.set_xticks(oks); ax1.set_xticklabels([f"{o:g}" for o in oks])
    ax2 = ax1.twinx()
    ax2.plot(oks, winrates, "s--", color="C2", lw=1.6, ms=6, alpha=0.8, label="win rate")
    ax2.set_ylabel("matched-seed win rate", color="C2"); ax2.set_ylim(0, 1.05)
    ax1.set_title("Entropic strength sweep: FR gain grows with the bottleneck")
    l1, lb1 = ax1.get_legend_handles_labels(); l2, lb2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, lb1 + lb2, loc="best")
    plt.tight_layout()
    out = os.path.join(PLOTS, "fig2_omega_sweep.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


def figure3(A):
    """FR-rate failure boundary: final L2(F) and ancestor ESS vs gamma."""
    gammas = [1.0, 3.0, 5.0, 10.0, 15.0, 25.0, 50.0]
    gv, med_l2f, q1_l2f, q3_l2f, med_ess, med_repl = [], [], [], [], [], []
    abf_ref = None
    pa = "stage4_gamma|abf|beta8|oin25|gamma15::l2_f_t"  # abf is gamma-independent
    if pa in A:
        abf_ref = np.median(A[pa][:, -1])
    for ga in gammas:
        pf = f"stage4_gamma|fr_estimated|beta8|oin25|gamma{ga:g}::l2_f_t"
        pe = f"stage4_gamma|fr_estimated|beta8|oin25|gamma{ga:g}::ess_t"
        if pf not in A:
            continue
        fl = A[pf][:, -1]
        q1, q2, q3 = np.percentile(fl, [25, 50, 75])
        gv.append(ga); med_l2f.append(q2); q1_l2f.append(q1); q3_l2f.append(q3)
        med_ess.append(np.median(A[pe][:, -1]) if pe in A else np.nan)
    if not gv:
        print("fig3: no data"); return
    fig, ax1 = plt.subplots(figsize=(7.6, 5.0))
    ax1.fill_between(gv, q1_l2f, q3_l2f, color=C_FR, alpha=0.18)
    ax1.plot(gv, med_l2f, "o-", color=C_FR, lw=2.2, ms=7, label="FR final $L^2(F)$")
    if abf_ref is not None:
        ax1.axhline(abf_ref, color=C_ABF, lw=1.8, ls="--", label="ABF baseline")
    ax1.set_xscale("log")
    ax1.set(xlabel=r"$\gamma$ (Fisher--Rao birth--death rate)",
            ylabel="final $L^2(F)$")
    ax1.set_xticks(gv); ax1.set_xticklabels([f"{x:g}" for x in gv])
    ax2 = ax1.twinx()
    ax2.plot(gv, med_ess, "s--", color="C2", lw=1.6, ms=6, alpha=0.85,
             label="ancestor ESS (windowed)")
    ax2.set_ylabel("median final ancestor ESS", color="C2")
    ax1.set_title("FR-rate safe regime and failure boundary")
    l1, lb1 = ax1.get_legend_handles_labels(); l2, lb2 = ax2.get_legend_handles_labels()
    ax1.legend(l1 + l2, lb1 + lb2, loc="best")
    plt.tight_layout()
    out = os.path.join(PLOTS, "fig3_gamma_failure.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


def figure4(A):
    """Conditional fidelity: empirical Var(Y|X in bin) vs analytic, ABF vs FR."""
    stage = "stage1_seeds"
    pref = f"{stage}|{{m}}|beta8|oin25|gamma15"
    if pref.format(m="abf") + "::cond_emp_var" not in A:
        stage = "stage0_reproduce"
        pref = f"{stage}|{{m}}|beta8|oin25|gamma15"
    centers = g(A, pref.format(m="abf") + "::cond_centers")
    ref = g(A, pref.format(m="abf") + "::cond_ref_var")
    if centers is None:
        print("fig4: no data"); return
    K = len(centers)
    xpos = np.arange(K)
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    width = 0.36
    for off, method, c, lbl in [(-width/2, "abf", C_ABF, "ABF"),
                                (+width/2, "fr_estimated", C_FR, "ABF + FR")]:
        stk = g(A, pref.format(m=method) + "::cond_emp_var")  # (seeds, K)
        if stk is None:
            continue
        m = np.nanmedian(stk, axis=0)
        lo = np.nanpercentile(stk, 25, axis=0); hi = np.nanpercentile(stk, 75, axis=0)
        ax.bar(xpos + off, m, width, color=c, alpha=0.85, label=lbl,
               yerr=[m - lo, hi - m], capsize=3, ecolor="k", error_kw=dict(lw=1))
    ax.plot(xpos, ref, "k_", ms=22, mew=2.6, label=r"analytic $1/(\beta\omega(x)^2)$")
    ax.set_xticks(xpos); ax.set_xticklabels([f"x={c:g}" for c in centers])
    ax.set(ylabel=r"$\widehat{\rm Var}(Y\mid X\in I_j)$",
           title="Conditional fidelity: FR preserves the analytic $Y\\mid X$")
    ax.legend()
    plt.tight_layout()
    out = os.path.join(PLOTS, "fig4_conditional_fidelity.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


def main():
    os.makedirs(PLOTS, exist_ok=True)
    A = load()
    figure1(A)
    figure2(A)
    figure3(A)
    figure4(A)


if __name__ == "__main__":
    main()
