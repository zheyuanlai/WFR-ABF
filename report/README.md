# ABF–Fisher-Rao report

A manuscript-style report answering, on a controlled 2D example with a
quadrature-exact reference free energy:

1. **Does adding Fisher–Rao (FR) to ABF help?**
2. **How should FR be added so that it helps?**
3. **What are the limitations, failure modes, and future directions?**

The compiled PDF is `main.pdf` (built from `main.tex`).

---

## 1. How the report was generated

Nothing here re-runs a simulation. Every number and figure is derived from the
already-computed production-study CSVs and the reference grid under
`../results/two_dim_xi_x/`. The pipeline is:

```
results/two_dim_xi_x/*.csv,*.npz
        │
        ▼
scripts/build_report_assets.py   →  figures/*.png, tables/*.{tex,csv,json}
        │
        ▼
scripts/check_report_numbers.py  →  verifies tables/numbers match the CSVs
        │
        ▼
main.tex (+ sections/, tables/, figures/, refs.bib)  →  main.pdf
```

* `scripts/build_report_assets.py` loads the merged CSVs (falling back to the
  `__shard_000` files only if a merged file is missing), checks the study
  metadata (number of runs/configs/seeds/methods; absence of NaNs), computes the
  summary tables under **two explicit selection rules**, renders the seven
  figures, and writes:
  * `tables/numbers.tex` — LaTeX `\newcommand` macros for every in-text number,
  * `tables/report_numbers.json` — the same numbers in machine-readable form.
  The report's prose uses the macros, so no figure can drift from the data.
* `scripts/check_report_numbers.py` independently recomputes the headline
  quantities from the CSVs and fails (non-zero exit) if `report_numbers.json`,
  the `numbers.tex` macros, or the main tables disagree with the recomputation.

## 2. Which result files are used

From `../results/two_dim_xi_x/`:

| File | Used for |
| --- | --- |
| `reference/reference_profile.csv` | reference `F_ref`, `F'_ref`, `p_ref` profiles (Fig. 1, Fig. 3) |
| `reference/reference_grid.npz` | potential `V(x,y)` and Boltzmann density (Fig. 1) |
| `production_gpu/production_gpu_final_summary.csv` | per-run final + integrated metrics → all tables, Fig. 4 |
| `production_gpu/production_gpu_config_summary.csv` | per-config medians/IQR → Fig. 7 heatmaps |
| `production_gpu/production_gpu_runs_long.csv` | per-snapshot time series → Fig. 2, Fig. 5 |
| `production_gpu/production_gpu_profiles.csv` | final F/F'/p profiles per run → Fig. 3 |
| `production_gpu/production_gpu_fr_events.csv` | FR birth–death event stats (loaded; diagnostics) |
| `production_gpu/production_gpu_conditional_diagnostics.csv` | conditional `Y\|X` errors → Fig. 6 |

The metastability case (§1) uses `../results/two_dim_xi_x/`. The two additional
cases re-render their figures in the report style (via `scripts/report_cases.py`)
and read their headline numbers from already-aggregated summaries:

| File | Used for |
| --- | --- |
| `../results/entropic_bottleneck/summaries/config_summary.csv` | EB β/ω/γ-sweep medians → `tab:eb_*`, EB macros |
| `../results/entropic_bottleneck/summaries/arrays.npz` | EB convergence + analytic-conditional arrays → `fig_eb_01`, `fig_eb_04` |
| `../results/wca_production/summaries/config_summary.csv` | WCA per-config medians → `tab:wca_*`, `fig_wca_03/05` |
| `../results/wca_production/summaries/winrates.csv` | WCA matched-seed gains/winrates → WCA macros, `fig_wca_06` |
| `../results/wca_production/summaries/timeseries_summary.csv` | WCA convergence bands → `fig_wca_01` |
| `../results/wca_production/summaries/profiles_summary.npz` | WCA final profiles + mechanism → `fig_wca_02/04` |

Unlike the metastability pipeline (which recomputes from per-run CSVs), the
EB/WCA numbers are read from the published config-summary/winrate medians;
`check_report_numbers.py` cross-checks them against those same summaries.

Merged files are preferred; the `*__shard_000.csv` shards are a fallback only.
`best_configs.csv`, `table_main_results.csv` and `table_production_gpu_top_configs.csv`
already exist in the repo but are **not** trusted blindly: the report recomputes
its tables so that the two selection rules (below) are never silently mixed.

## 3. How to rebuild figures and tables

```bash
cd report
python scripts/build_report_assets.py     # regenerate figures/ and tables/
python scripts/check_report_numbers.py     # verify numbers vs the CSVs
```

Both scripts only need `numpy`, `pandas`, and `matplotlib` (no PyTorch). They
resolve paths from their own location, so they can be run from any directory.

## 4. How to compile the PDF

```bash
cd report
make            # = make assets check pdf  (regenerate, verify, compile)
# or just compile, assuming assets already exist:
make pdf
```

`make` uses `latexmk` if available, otherwise falls back to
`pdflatex` + `bibtex` + two more `pdflatex` passes. Only standard LaTeX packages
are required (`amsmath`, `graphicx`, `booktabs`, `subcaption`, `hyperref`,
`cleveref`, `siunitx`, `natbib`, `geometry`, `xcolor`). `make clean` removes
build artefacts; `make distclean` also removes generated figures/tables/PDF.

## 5. Which numbers are *diagnostic* vs *practical* (read this)

This distinction is essential to interpreting the report correctly.

* **Practical methods** — the only claims about a usable algorithm compare
  **ABF only** vs **ABF + FR with the *estimated* target**
  (`q_t^est ∝ exp(-β(F̂_t^target − B_t))`). The **uniform** target is treated as
  a *naive target ablation* (does flattening the marginal suffice?), not as the
  main method.
* **Diagnostic control — the oracle target.** `q_t^oracle ∝ exp(-β(F_ref − B_t))`
  uses the **reference free energy `F_ref`**, i.e. *exactly the unknown the
  method is trying to compute*. The oracle is therefore **not a usable
  algorithm**. It is reported only as a **positive / ideal-information control**
  that shows whether the FR mechanism has potential when the target marginal
  direction is accurate. The gap between the estimated and oracle targets
  measures the headroom available from a better online target estimate; it is
  **never** a recommendation to "use the oracle".

The same care applies to the **two selection rules**, kept in separate tables:

* `tables/best_configs_by_integrated_F.{csv,tex}` — selected by smallest median
  **integrated** `L²(F)` → rewards **fast convergence** (anytime performance).
  This is the **main-text** table (`\ref{tab:main_integrated}`).
* `tables/best_configs_by_final_F.{csv,tex}` — selected by smallest median
  **final** `L²(F)` → rewards **final-budget accuracy**. This is the
  **supplement** table (`\ref{tab:main_final}`).

These rules need not select the same hyperparameters (e.g. the estimated target
is best served by an aggressive rate for integrated error but a gentle, delayed
rate for final error), and the report is explicit that final-budget improvements
are modest.

A further honest caveat: `conditional_diagnostics.csv` is keyed only by
`(method, target_type, seed)` — it is **not configuration-resolved**, so Fig. 6
supports only the *aggregate* statement that successful FR does not corrupt
`Y|X`, not a per-hyperparameter good-vs-bad attribution.
