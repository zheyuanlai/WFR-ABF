# WCA dimer: does Fisher--Rao birth--death help ABF compute the free energy?

A production study of **adaptive biasing force (ABF)** versus **ABF + Fisher--Rao
(FR) birth--death with an online estimated target**, on the Weeks--Chandler--Andersen
(WCA) dimer in a dense solvent.

> **Scope and honesty.** The practical claim is always **ABF only vs. ABF + FR
> estimated target**. `fr_uniform` (uniform-in-z target) and `fr_oracle` (target
> built from the TI reference) are *diagnostics*: `fr_oracle` uses the unknown
> answer and is **not a deployable algorithm**; it bounds the headroom only.
> Negative and failure results are reported, not hidden.

---

## 1. Problem setup: WCA dimer and reaction coordinate

A dimer of two tagged particles sits in a periodic 2D bath of `n_dim^2 - 2` WCA
solvent particles (`n_dim=10`, lattice spacing `a`, box `L = n_dim * a`). All
particles interact through the purely repulsive WCA potential (Lennard--Jones cut
and shifted at its minimum `r0 = 2^{1/6} sigma`). The two dimer particles
additionally feel a symmetric double-well bond

  V_dim(r) = h * (1 - ((r - r0 - w)/w)^2)^2,

with barrier height `h=2` and width `w=2`. The reaction coordinate is the scaled
bond length

  z(q) = (|q_1 - q_2| - r0) / (2 w).

`z ≈ 0` is the **compact** state, `z ≈ 1` the **stretched** state, and
`z ∈ [0.25, 0.75]` the **transition** region. The free-energy profile `F(z)` has a
barrier near `z = 0.5`; the dimer must push solvent aside to cross it, so the
barrier is partly **entropic / crowding**-driven. Dynamics are overdamped Langevin
at `beta = 1`, `dt = 0.002`, with periodic boundaries and the minimum-image
convention. Default system: `a = 1.5`, `sigma = eps = 1`.

The local mean force along `z` is estimated from the constrained projection of the
physical forces (energetic term) plus the geometric/entropic `2w/(beta r)` term;
ABF accumulates a kernel estimate `Â(z)` of `F(z)` and biases the dynamics by
`B_n(z) = s_n Â_n(z)` with a warm-up ramp `s_n`.

## 2. Methods

All methods share the same ABF machinery, initial condition, and seed. They differ
only in the FR birth--death correction layered on top.

- **`abf`** — ABF only. Baseline.
- **`fr_estimated`** (the practical method) — ABF + FR birth--death. The FR target
  marginal is built from an **online EMA estimate** of the ABF free energy,
  maintained separately from the bias:

    F^target_{n+1}(z) = (1-r) F^target_n(z) + r Â_n(z),   (r = `target_ema_rate`)
    q_n(z) ∝ exp[ -beta ( F^target_n(z) - B_n(z) ) ],   normalized on the z-grid.

  No TI reference enters this target. As `Â_n → F` and the bias ramp saturates,
  `q_n` flattens, so FR pressure fades as ABF converges.
- **`fr_uniform`** (diagnostic ablation) — same machinery, target `q_n(z) = const`.
- **`fr_oracle`** (diagnostic positive control) — `q_n(z) ∝ exp[-beta(F_ref(z) -
  B_n(z))]` using the TI reference `F_ref`. **Uses the unknown answer; not
  deployable.** Bounds how much an ideal target could help.

**FR score and birth--death.** At each FR step the per-replica score is

  S_i = log p̂_n(Z_i) - log q_n(Z_i) - KL(p̂_n || q_n),

(p̂ = reflected-boundary KDE marginal of z), clipped to `±score_clip` and recentered
to mean zero. Birth--death is **fixed population**: replicas with `S_i > 0` (over-
represented) die with probability `1 - exp(-fr_rate · S_i · dt_eff)`; each death is
replaced by cloning a replica drawn with weight `∝ max(-S_i, 0)`; the number of
simultaneous events is capped at `max_event_fraction · N`. The **entire** WCA
configuration is copied on a clone, not just `z`. Replica lineage is tracked so we
can report the **ancestor effective sample size** `ESS = 1 / Σ_a w_a^2`.

**No-leakage guarantee.** `run_sampler_gpu` calls `assert_no_oracle_leakage`: only
`fr_oracle` may receive the reference; `abf`/`fr_estimated`/`fr_uniform` raise if a
reference target is passed. The TI reference is used only to (a) build the
`fr_oracle` diagnostic and (b) score final accuracy. An 8-agent adversarial review
confirmed no path by which TI reaches the `fr_estimated` target.

## 3. Experimental design

- **TI reference (evaluation only).** Constrained-MD thermodynamic integration on a
  51-point z-grid, smoothed and interpolated to the 160-point evaluation grid,
  computed once per system `a` and cached. Free energies aligned by subtracting the
  mean offset on the **interior window** `z ∈ [-0.1, 1.1]` (wall-affected edge bins
  excluded) before computing any `L^2` error.
- **Metrics.** Primary: final `L^2(F)`. Secondary: integrated `∫ L^2(F_t) dt`,
  final `L^2(F')`. Region-wise `L^2` over compact/transition/stretched (all on the
  interior window). Matched-seed win rate `P_seed(L^2_F(FR) < L^2_F(ABF))`.
- **Stages.** (0) reproduce the demo; (A) staged pilot tuning around the demo
  winner; (B) main comparison, 10 seeds, all methods; (failure) `fr_rate` ladder up
  to 2.0; (D) difficulty axes — budget `n_steps`, replica count, solvent crowding
  `a`.

Default tuned FR config (`fr_est_tuned`): `fr_rate=0.10`, `target_ema_rate=0.005`,
`max_event_fraction=0.02`, `fr_every=5`, `fr_start_steps=20000`, `score_clip=2.0`.

## 4. Main results

Stage B: 10 seeds (`0,1,2,3,4,5,6,7,42,123`), `N=1024`, `n_steps=250000`,
`dt=0.002`, matched initial condition per seed. Errors on the interior window
`z ∈ [-0.1, 1.1]`, additive constant aligned. Median over seeds; gain = matched-seed
median of `100·(L²_F(ABF) − L²_F(method))/L²_F(ABF)`.

| method | final L²(F) | gain vs ABF | win rate | int. L²(F) gain | final L²(F') | ancestor ESS | replacements |
|---|---|---|---|---|---|---|---|
| `abf` (baseline) | 0.0835 | — | — | — | 0.322 | — | 0 |
| `fr_est_gentle` (rate 0.05) | 0.0562 | **+32.8%** | 10/10 | +17.1% | 0.276 | 157 | 2.7k |
| **`fr_est_tuned` (rate 0.10)** | **0.0361** | **+55.0%** | **10/10** | +30.6% | 0.251 | 86 | 5.0k |
| `fr_est_strong` (rate 0.20) | 0.0329 | **+58.9%** | 10/10 | +40.5% | 0.259 | 45 | 8.9k |
| `fr_est_aggressive` (rate 0.50) | 0.0935 | **−13.1%** | 1/10 | +11.9% | 0.388 | 22 | 17.6k |
| `fr_uniform` (diagnostic) | 0.0372 | +54.6% | 10/10 | +31.3% | 0.253 | 90 | 5.0k |
| `fr_oracle` (diagnostic) | 0.0387 | +52.7% | 10/10 | +25.8% | 0.248 | 100 | 4.1k |

**Headline.** The practical method **ABF + FR estimated target reliably beats ABF**:
`fr_est_tuned` cuts the median final `L²(F)` from **0.0835 to 0.0361 (+55.0%)** and
wins on **all 10/10 matched seeds**. `fr_est_strong` (rate 0.20) is marginally better
still (+58.9%) but with lower diversity (ESS 45 vs 86). We report `fr_est_tuned`
(rate 0.10) as the recommended operating point: near-maximal gain with a comfortable
diversity margin. (Figures 1, 3.)

**Where the gain lives.** The improvement is largest in **final `L²(F)`** (+55%) and
**integrated `L²(F)`** (+31%); the **mean-force `L²(F')`** improves only modestly
(0.322 → 0.251, ≈ 22%). We therefore do **not** claim a large mean-force improvement;
the headline is the free-energy error.

**Almost no headroom to the oracle.** The estimated target (0.0361) is **as accurate
as the oracle target** (0.0387) — the online EMA target already extracts essentially
all the benefit an ideal target could provide on this system. The oracle is a
diagnostic that uses the unknown answer and is **not deployable**; it does not beat
the practical method here.

**The uniform target is nearly as good.** `fr_uniform` (0.0372, +54.6%) ≈
`fr_est_tuned`. On this WCA dimer the dominant effect is **balanced birth--death
resampling itself** — pulling replicas out of the over-sampled compact basin — rather
than the precise shape of the target. This is an honest limitation of the "smart
target" story for this particular system; see §5.

## 5. Mechanism analysis

The core hypothesis was that FR helps by reallocating replicas to under-sampled RC
regions, improving the ABF mean-force estimate in the transition window. The
diagnostics support a **refined** version of this.

- **Coverage is not visibly reshaped at the endpoint.** Final RC-marginal region
  fractions (compact / transition / stretched) and the fraction of effective ABF
  samples falling in the transition window `[0.25,0.75]` are nearly identical across
  `abf`, `fr_est_tuned`, `fr_uniform`, `fr_oracle` (transition occupancy ≈ 0.31–0.33,
  `N_eff` in transition ≈ 35% for all). FR does not converge to a different stationary
  RC distribution.
- **The gain appears immediately after FR activates and is sustained.** ABF and FR
  trajectories are identical until FR turns on at `t = 40` (step 20k); from `t = 80`
  onward FR is consistently below ABF (e.g. `t = 80`: ABF 0.114 vs FR 0.079; `t = 500`:
  ABF 0.083 vs FR 0.036). (Figure 1.)
- **Interpretation.** Birth--death is acting as a **variance-reduction resampler for
  the ABF estimator**: replicas stuck in over-represented configurations are replaced
  by whole-configuration clones of under-represented ones, raising the number of
  effectively-independent force samples per `z`-bin even when the marginal looks the
  same. Because the *whole WCA configuration* is copied (not just `z`), clones inject
  decorrelated solvent arrangements. This accelerates convergence of `F̂` (final and
  integrated `L²(F)`) without requiring a different stationary coverage — consistent
  with `fr_uniform ≈ fr_estimated ≈ fr_oracle`. (Figure 4: `N_eff(z)`, per-bin
  `F'(z)` error, birth/death locations, transition occupancy over time.)

So the mechanism is **resampling-driven variance reduction**, of which "better
transition coverage" is a contributing but not dominant component on this system.

## 6. Failure modes

The `fr_rate` ladder (pilot + failure stages) maps the safe operating regime cleanly
(Figure 5):

- **Too gentle** (`fr_rate ≤ 0.05`): few birth--death events, small gain (+24–33%).
- **Sweet spot** (`fr_rate ≈ 0.1–0.2`): maximal gain (+55–59%), ancestor ESS stays
  healthy (45–90 of 1024 distinct ancestors).
- **Too aggressive** (`fr_rate ≥ 0.5`): **gain reverses**. At rate 0.50, median final
  `L²(F)` is 0.0935 — **worse than ABF (−13%)**, winning only 1/10 seeds — with ~17.6k
  replacements and ancestor ESS collapsing to ~22. Over-aggressive birth--death
  destroys the replica diversity the ABF mean-force estimator depends on, so the
  estimator degrades faster than coverage improves. This is the predicted failure mode
  and is included as evidence, not hidden.

## 7. Difficulty dependence

Stage D varies one difficulty axis at a time (5 seeds each), comparing `abf` vs
`fr_est_tuned` on matched seeds. Gain = matched-seed median of
`100·(L²_F(ABF) − L²_F(FR))/L²_F(ABF)`; win rate over the 5 seeds. (Figure 6.)

**Budget — `n_steps` (N=1024, a=1.5).** FR wins at every budget and the *relative*
advantage **grows** with budget rather than washing out:

| n_steps | ABF L²(F) | FR L²(F) | gain | win | int. L²(F) gain | ESS |
|---|---|---|---|---|---|---|
| 50k | 0.1077 | 0.0636 | +35.5% | 5/5 | +5.5% | 397 |
| 100k | 0.0915 | 0.0476 | +47.4% | 5/5 | +22.9% | 211 |
| 150k | 0.0867 | 0.0399 | +53.9% | 5/5 | +25.6% | 131 |
| 250k | 0.0792 | 0.0358 | +54.4% | 5/5 | +33.5% | 89 |
| 400k | 0.0815 | 0.0342 | +56.9% | 5/5 | +42.8% | 60 |

ABF's final error **plateaus** around 0.08 after ≈150k, whereas FR keeps descending
(0.064 → 0.034). So FR is **not merely an early-time accelerator** that ABF later
catches: the gap widens with budget. (Ancestor ESS falls with budget because more
total birth–death events accumulate over a longer run, but stays healthy.)

**Replica count — `n_replicas` (250k, a=1.5).** FR wins 5/5 at every population size,
including small ensembles:

| n_replicas | ABF L²(F) | FR L²(F) | gain | win | ESS |
|---|---|---|---|---|---|
| 256 | 0.0851 | 0.0395 | +51.2% | 5/5 | 21 |
| 512 | 0.0838 | 0.0377 | +54.9% | 5/5 | 45 |
| 1024 | 0.0833 | 0.0343 | +58.8% | 5/5 | 80 |
| 2048 | 0.0850 | 0.0370 | +56.2% | 5/5 | 172 |

The gain is **roughly flat in replica count** (≈ +51–59%). FR helps even at N=256,
where the ancestor ESS is only ~21 — the resampling benefit does not require a large
population. ABF's own accuracy is nearly independent of N here (it is limited by
per-replica decorrelation, not ensemble size), so FR's relative edge persists.

**Solvent crowding — `a` (250k, N=1024; smaller `a` = denser bath, harder barrier).**
This axis changes the *absolute* difficulty by more than an order of magnitude:

| a | regime | ABF L²(F) | FR L²(F) | gain | win | int. L²(F) gain |
|---|---|---|---|---|---|---|
| 1.35 | dense / hard | 0.3649 | 0.2246 | +38.3% | 5/5 | +28.9% |
| 1.5 | default | 0.0817 | 0.0374 | +54.7% | 5/5 | +34.8% |
| 1.7 | loose / easy | 0.0228 | 0.0115 | +49.6% | 5/5 | +25.2% |

At `a=1.35` the crowded solvent makes the entropic barrier severe and ABF alone is
**~4× worse in absolute error** (0.365 vs 0.082 at default). FR still wins 5/5 and
delivers by far the **largest absolute error reduction** (0.365 → 0.225, Δ ≈ 0.14 —
versus Δ ≈ 0.044 at default and Δ ≈ 0.011 in the easy case). The *percentage* gain is
slightly lower at `a=1.35` (+38%) because the hard problem leaves a larger residual
error for both methods, but in the units that matter for a practitioner — absolute
free-energy error removed — the crowded regime is where FR helps most.

**Summary.** Across all three difficulty axes and **every one of the 60 Stage-D runs**,
`fr_est_tuned` beats matched-seed ABF (**12/12 cells at 5/5 wins**). The relative gain
grows with budget, is insensitive to replica count, and the largest absolute gains
appear in the hardest (densest) regime — consistent with the variance-reduction
mechanism of §5 mattering most when the ABF estimator is otherwise starved of
decorrelated transition-region samples.

## 8. Conclusion

On the WCA dimer, **ABF + Fisher--Rao birth--death with an online estimated target
reliably and substantially improves ABF free-energy accuracy** — a median **+55%**
reduction in final `L²(F)` at the tuned operating point, winning on **10/10** seeds,
with the gain concentrated in the free-energy error (final and integrated) rather than
the mean force. The benefit is driven by **balanced birth--death resampling acting as
variance reduction for the ABF estimator** after warm-up; on this system the precise
target shape matters little (uniform ≈ estimated ≈ oracle), and the estimated target
already reaches the oracle's accuracy, so there is little headroom for a smarter
target. The method has a clear safe regime: `fr_rate ≈ 0.1–0.2`. **Too-aggressive
birth--death (`fr_rate ≥ 0.5`) reverses the benefit**, dropping below ABF as replica
diversity collapses. Across all three Stage-D difficulty axes the advantage is robust:
it **grows with simulation budget** (+36% at 50k → +57% at 400k, with ABF plateauing
while FR keeps improving), is **insensitive to replica count** (+51–59% from N=256 to
2048), and produces its **largest absolute error reduction in the densest, hardest
solvent** (`a=1.35`: 0.365 → 0.225), winning **5/5 seeds in every one of the 12 cells**.
The estimated-target method is the deployable recommendation; the oracle is diagnostic
only and is never used outside evaluation.

---

## Reproduction

```bash
source /home/zheyuanlai/miniconda3/etc/profile.d/conda.sh && conda activate abffr
# smoke
CUDA_VISIBLE_DEVICES=7 python scripts/smoke_wca_production.py
# a stage across GPUs 1-7
bash scripts/launch_stage_multigpu.sh main
# analyze + plot
python scripts/analyze_wca_production.py
python scripts/plot_wca_production.py
```

Raw per-run `.npz`: `results/wca_production/raw/`. Summaries:
`results/wca_production/summaries/`. Figures: `results/wca_production/plots/`.
