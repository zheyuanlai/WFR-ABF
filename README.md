# Adaptive Biasing Force (ABF) with Fisher-Rao for Free Energy Computation

This repository contains exploratory notebooks (`ABF-FR.ipynb`, `eABF.ipynb`,
`WCA_dimer.ipynb`) and a reproducible, script-based study of whether a
Fisher--Rao (FR) birth--death/resampling correction helps or hurts ABF on a 2D
model potential. The notebooks are kept as reference only; the study below is
driven entirely from the command line.

## 2D ABF-FR Fisher--Rao Ablation Study

### Scientific goal

We study, on a controlled 2D example with a known reference free energy, when a
Fisher--Rao birth--death correction **helps** ABF, when it **hurts**, and *why*.
Concretely we ask:

* Does FR improve the reaction-coordinate (`x`) marginal exploration while
  **preserving** the conditional distribution `Y | X = x` that ABF relies on?
* When FR hurts, is the cause a **noisy/inaccurate estimated target**, or is the
  birth--death **mechanism itself** harmful (e.g. it distorts `Y | X = x` or
  collapses particle diversity)?

The oracle and uniform targets are diagnostic controls that bracket the
estimated (practical) method between an ideal target and a naive one.

### Potential and reaction coordinate

Reaction coordinate `xi(x, y) = x`. The local mean force along the reaction
coordinate is therefore `f(x, y) = dV/dx`. The potential is

```
V(x, y) =  3 e^{-x^2} ( e^{-(y - 1/3)^2} - e^{-(y - 5/3)^2} )
         - 5 e^{-y^2} ( e^{-(x - 1)^2} + e^{-(x + 1)^2} )
         + 0.2 x^4 + 0.2 (y - 1/3)^4
```

It has two metastable basins near `x = -1` and `x = +1` separated by a barrier
near `x = 0` (plus an upper sub-basin near `(0, 5/3)`), so plain dynamics is
metastable in `x`. Default inverse temperature `beta = 4`.

### Reference quantities (by `y`-quadrature)

Free energy is defined up to a constant; all profiles are centred
(`F <- F - mean(F)`, over the evaluation window) before computing L2 errors.

```
F_ref(x)   = -(1/beta) * log integral_y e^{-beta V(x,y)} dy   + C
F'_ref(x)  = [ integral_y dV/dx e^{-beta V} dy ] / [ integral_y e^{-beta V} dy ]
p_ref(x)   proportional to integral_y e^{-beta V(x,y)} dy      (unbiased x-marginal)
```

### Biased dynamics and ABF

`n_particles` independent overdamped Langevin replicas share one kernel ABF
mean-force estimate `F'_hat_t(x) ~= E[dV/dx | X = x]`. With bias potential
`B_t(x) = integral F'_hat_t`:

```
dX_t = [ -dV/dx(X_t,Y_t) + F'_hat_t(X_t) ] dt + sqrt(2/beta) dW^x_t
dY_t = [ -dV/dy(X_t,Y_t)                 ] dt + sqrt(2/beta) dW^y_t
```

### Methods compared

| method              | FR target `q_t(x)`                                             |
| ------------------- | -------------------------------------------------------------- |
| `abf_only`          | none (baseline)                                                |
| `abf_fr_estimated`  | `proportional to exp(-beta (F_hat_target(x) - B_t(x)))` (EMA)   |
| `abf_fr_uniform`    | `Uniform` on the `x`-grid                                       |
| `abf_fr_oracle`     | `proportional to exp(-beta (F_ref(x) - B_t(x)))`  *(diagnostic)*|

FR nudges the particle `x`-marginal `p_hat_t(x)` toward `q_t(x)` by gentle
fixed-`N` birth--death: `S_i = log(p_hat(X_i)/q(X_i)) - KL(p_hat||q)`; particles
with `S_i > 0` die and with `S_i < 0` clone, with rate `gamma` (soft-ramped after
a burn-in), score clipping and a `max_event_fraction` safety cap.

**Why oracle is only a diagnostic.** The oracle target is allowed to use the
reference free energy `F_ref` (which is exactly what ABF is trying to estimate),
so it is *not* a usable method. It is an upper bound on how much an ideal,
perfectly-informed FR target could help, used to separate "FR target is noisy"
from "FR mechanism is harmful". The practical method is `abf_fr_estimated`.

### Layout

```
configs/   two_dim_xi_x_smoke.yaml      tiny, fast (pipeline check)
           two_dim_xi_x_tuning.yaml     moderate-but-large production grid
src/abffr/ potentials, reference, simulation, metrics, diagnostics, plotting, io_utils
scripts/   run_reference_2d.py, run_abf_fr_grid.py, plot_abf_fr_study.py, make_report_tables.py
results/two_dim_xi_x/{reference,tuning,eval,figures_tuning,figures_eval}/
```

### How to run

All commands are deterministic given the seeds in the config. Outputs go under
`results/two_dim_xi_x/`.

1. Reference (fast):

   ```bash
   python scripts/run_reference_2d.py --config configs/two_dim_xi_x_smoke.yaml
   ```

   Writes `reference/reference_profile.csv`, `reference/reference_grid.npz` and
   four reference figures.

2. Smoke test of the full pipeline (a few minutes):

   ```bash
   python scripts/run_abf_fr_grid.py  --config configs/two_dim_xi_x_smoke.yaml --stage smoke
   python scripts/plot_abf_fr_study.py --stage tuning
   python scripts/make_report_tables.py --stage tuning
   ```

   `--stage smoke` writes into `results/two_dim_xi_x/tuning/` so it can be
   inspected with the same `--stage tuning` plotting/table commands.

3. Tuning grid (LONG, run explicitly only when you want it):

   ```bash
   python scripts/run_abf_fr_grid.py  --config configs/two_dim_xi_x_tuning.yaml --stage tuning
   python scripts/plot_abf_fr_study.py --stage tuning
   python scripts/make_report_tables.py --stage tuning
   ```

   The tuning grid is `4 methods x gamma x eta x burn-in x 5 seeds` at
   `100k` steps with `1000` particles. The kernel ABF estimator is `O(n_grid x
   n_particles)` per step, so this is a multi-CPU-hour production run; preview it
   with `--dry-run` and cap it with `--max-runs`. It writes `best_configs.csv`.

4. Eval of the selected configs (reads `tuning/best_configs.csv`):

   ```bash
   python scripts/run_abf_fr_grid.py  --config configs/two_dim_xi_x_tuning.yaml --stage eval
   python scripts/plot_abf_fr_study.py --stage eval
   python scripts/make_report_tables.py --stage eval
   ```

### Output files

```
reference/reference_profile.csv        x, F_ref, Fprime_ref, p_ref
reference/reference_grid.npz           x_grid, y_grid, V_grid, rho_grid
tuning/tuning_runs_long.csv            per (run, snapshot) time series
tuning/tuning_final_summary.csv        per run final + time-integrated metrics
tuning/tuning_config_summary.csv       per config, median/IQR over seeds
tuning/tuning_profiles.csv             per run final F'(x), F(x), p(x), q(x)
tuning/tuning_fr_events.csv            per (run, snapshot) FR birth-death stats
tuning/tuning_conditional_diagnostics.csv   conditional Y|X errors per x-bin
tuning/best_configs.csv                best config per target type (+ ABF-only)
tuning/table_tuning_top_configs.csv    ranked top configs
eval/eval_*.csv                        same schema, for the selected configs
eval/table_main_results.csv           main method comparison table
figures_tuning/, figures_eval/         manuscript-style figures
```

### Tuning criterion

The best estimated-target config is chosen by the **primary** metric
`median_seeds integral_0^T ||F_hat_t - F_ref||_{L2} dt`, tie-broken by the
**secondary** metric `median_seeds ||F_hat_T - F_ref||_{L2}`, subject to safety
filters (no NaNs; typical and maximum FR event fraction within
`max_event_fraction`). See `select_best_configs` in `scripts/run_abf_fr_grid.py`.

### Interpretation of key metrics

* `l2_F`, `l2_Fprime` -- centred L2 error of the free energy / mean-force profile
  against the reference, over the interior window `[-2.5, 2.5]` (the reflecting
  walls and the `x^4` confinement make the edges uninformative). Lower is better.
* `integrated_l2_F` -- the primary tuning metric: time-integral of `l2_F`,
  rewarding *fast* convergence, not just the final value.
* `marginal_l2_uniform` / `marginal_l2_target` -- how flat the `x`-marginal is /
  how close it is to the FR target.
* `barrier_crossings` -- cumulative `x = 0` crossings (a direct exploration
  measure). FR should *increase* this if it helps exploration.
* `fr_event_fraction` -- fraction of particles undergoing birth/death per FR
  application; must stay small (`<= max_event_fraction`) for FR to be "gentle".
* `conditional_l2_y` / `conditional_kl_y` -- the crux: distance between the
  empirical `p_hat(y | x in bin)` and the reference `p_ref(y | x)`. If FR flattens
  the `x`-marginal but **inflates** this, the birth--death mechanism is corrupting
  the conditional sampling that ABF depends on -- FR is hurting for a structural
  reason, not just a noisy target.
* `prob_beats_abf` -- matched-seed probability that a method's final `l2_F` beats
  ABF-only.

### Notes / caveats

* The **smoke** settings (small `gamma`, `5000` steps) intentionally produce only
  tiny FR event rates, so FR is nearly indistinguishable from ABF there -- the
  smoke run validates the pipeline, not the science. Meaningful FR effects appear
  only at the tuning scale (`gamma` up to `0.1`, `100k` steps).
* `tqdm` is an optional progress-bar dependency; the scripts run without it.
* Results files (`*.npz`, `*.png`, `*.csv` under `results/`) are git-ignored by
  default; force-add a specific artifact with `git add -f` if you want it tracked.
