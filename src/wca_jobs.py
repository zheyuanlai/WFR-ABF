"""Job expansion, run IO, and the single-run pipeline for the WCA production study.

Shared by run_wca_production.py / analyze_wca_production.py / plot_wca_production.py.
A *job* is a fully-specified single run (one method, one seed, one budget/replica/
crowding setting). Jobs are deterministic in their RunSpec, hashed for idempotency,
and saved one-.npz-per-run so interrupted sweeps never lose completed work.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from dataclasses import dataclass, asdict, replace as dc_replace

import numpy as np

import wca_abffr_core as core


# ---------------------------------------------------------------------------
# RunSpec: a single fully-specified run.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RunSpec:
    stage: str
    name: str            # human label of the method config (e.g. "fr_est_tuned")
    method: str          # sampler type: abf / fr_estimated / fr_uniform / fr_oracle
    seed: int
    n_steps: int
    n_replicas: int
    a: float
    fr_rate: float
    target_ema_rate: float
    max_event_fraction: float
    fr_every: int
    fr_start_steps: int
    score_clip: float
    save_every: int

    def spec_hash(self) -> str:
        return hashlib.md5(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()[:12]

    def run_id(self) -> str:
        return (f"{self.stage}__{self.name}__seed{self.seed}"
                f"__N{self.n_replicas}__T{self.n_steps}__a{self.a:g}__{self.spec_hash()}")


def build_params(spec: RunSpec) -> "core.DimerWCAParams":
    return core.DimerWCAParams(a=float(spec.a))


def build_sim(spec: RunSpec, base: dict) -> "core.SimConfig":
    """SimConfig from the YAML base block, overridden by this spec."""
    kw = dict(base)
    kw.update(
        n_replicas=int(spec.n_replicas), n_steps=int(spec.n_steps), seed=int(spec.seed),
        save_every=int(spec.save_every), fr_rate=float(spec.fr_rate),
        target_ema_rate=float(spec.target_ema_rate),
        max_event_fraction=float(spec.max_event_fraction), fr_every=int(spec.fr_every),
        fr_start_steps=int(spec.fr_start_steps), score_clip=float(spec.score_clip),
    )
    valid = core.SimConfig.__dataclass_fields__.keys()
    kw = {k: v for k, v in kw.items() if k in valid}
    return core.SimConfig(**kw)


# ---------------------------------------------------------------------------
# TI reference cache (keyed by crowding `a`; eval-only, never enters fr_estimated).
# ---------------------------------------------------------------------------
def ti_cache_path(cache_dir: str, a: float, n_grid: int) -> str:
    if abs(a - 1.5) < 1e-9:
        # default system keeps the original cache name (shared with the demo)
        return os.path.join(cache_dir, "wca_ti_reference.npz")
    return os.path.join(cache_dir, f"wca_ti_reference_a{a:g}_g{n_grid}.npz")


_REF_CACHE: dict = {}


def get_reference(spec: RunSpec, base: dict, engine, cache_dir="cache", verbose=False):
    """Load/compute the TI reference for this spec's system (memoised per process)."""
    sim = build_sim(spec, base)
    key = (round(float(spec.a), 6), int(sim.n_grid), float(sim.z_min), float(sim.z_max))
    if key in _REF_CACHE:
        return _REF_CACHE[key]
    params = build_params(spec)
    ti = core.TIConfig(z_min=sim.z_min, z_max=sim.z_max, dt=sim.dt)
    path = ti_cache_path(cache_dir, spec.a, sim.n_grid)
    ref = core.load_or_compute_ti_reference(path, params, sim, ti, engine, verbose=verbose)
    _REF_CACHE[key] = ref
    return ref


# ---------------------------------------------------------------------------
# Per-run output files (one .npz per run -> interrupt-safe, idempotent).
# ---------------------------------------------------------------------------
def run_npz_path(raw_dir: str, spec: RunSpec) -> str:
    return os.path.join(raw_dir, spec.run_id() + ".npz")


def run_is_valid(path: str) -> bool:
    if not os.path.exists(path):
        return False
    try:
        d = np.load(path, allow_pickle=True)
        ok = ("l2_f" in d.files and "times" in d.files
              and np.isfinite(float(d["l2_f"])) and not bool(d.get("had_nan", np.array(False))))
        d.close()
        return bool(ok)
    except Exception:
        return False


def execute_run(spec: RunSpec, base: dict, engine, cache_dir="cache", verbose=False):
    """Run one job and return a flat dict of scalars + arrays ready to save."""
    params = build_params(spec)
    sim = build_sim(spec, base)
    ref = get_reference(spec, base, engine, cache_dir=cache_dir, verbose=verbose)
    # Matched IC across methods: lattice IC depends only on (seed, n_replicas, params).
    ic = core.lattice_initial_conditions(params, sim.n_replicas, engine.device, engine.dtype, seed=sim.seed)
    oracle_fe = ref["free_energy"] if spec.method == "fr_oracle" else None

    t0 = time.perf_counter()
    diag = core.run_sampler_gpu(spec.method, params, sim, engine, initial_q=ic,
                                oracle_free_energy=oracle_fe, collect_diagnostics=True, verbose=verbose)
    fin = core.final_l2_errors(diag, ref, sim)
    ts = core.timeseries_l2(diag, ref, sim)

    had_nan = bool(np.isnan(diag["mean_force"][-1]).any() or np.isnan(diag["pmf"][-1]).any())
    out = {
        # identity / metadata
        "run_id": spec.run_id(), "spec_hash": spec.spec_hash(),
        "spec_json": json.dumps(asdict(spec), sort_keys=True),
        "stage": spec.stage, "name": spec.name, "method": spec.method,
        "seed": spec.seed, "n_steps": spec.n_steps, "n_replicas": spec.n_replicas,
        "a": spec.a, "fr_rate": spec.fr_rate, "target_ema_rate": spec.target_ema_rate,
        "max_event_fraction": spec.max_event_fraction, "fr_every": spec.fr_every,
        "fr_start_steps": spec.fr_start_steps, "score_clip": spec.score_clip,
        "config_hash": sim.config_hash(), "core_version": "wca_prod_v1",
        "runtime_seconds": diag["runtime_seconds"], "wall_seconds": time.perf_counter() - t0,
        "had_nan": had_nan, "total_replacement_events": diag["total_replacement_events"],
        # final scalar metrics
        "l2_f": fin["l2_f"], "l2_fp": fin["l2_fp"], "integrated_l2_f": ts["integrated_l2_f"],
        "l2_f_compact": fin["l2_f_compact"], "l2_f_transition": fin["l2_f_transition"],
        "l2_f_stretched": fin["l2_f_stretched"],
        "l2_fp_compact": fin["l2_fp_compact"], "l2_fp_transition": fin["l2_fp_transition"],
        "l2_fp_stretched": fin["l2_fp_stretched"],
        "final_ancestor_ess": float(diag["ancestor_ess"][-1]),
        "final_n_unique_ancestor": int(diag["n_unique_ancestor"][-1]),
        # grid + reference (small; duplicated per run for self-containment)
        "grid": ref["grid"], "ref_free_energy": ref["free_energy"], "ref_mean_force": ref["mean_force"],
        # timeseries
        "times": ts["times"], "l2_f_t": ts["l2_f_t"], "l2_fp_t": ts["l2_fp_t"],
        "repl_cumulative": diag["repl_cumulative"],
        "frac_compact": diag["frac_compact"], "frac_transition": diag["frac_transition"],
        "frac_stretched": diag["frac_stretched"],
        "ancestor_ess_t": diag["ancestor_ess"], "pq_l2_t": diag["pq_l2"], "kl_pq_t": diag["kl_pq"],
        # final profiles + mechanism
        "final_mean_force": diag["mean_force"][-1], "final_pmf": diag["pmf"][-1],
        "final_p_hat": (diag["p_hat"][-1] if len(diag["p_hat"]) else np.full_like(ref["grid"], np.nan)),
        "final_q_target": (diag["q_target"][-1] if len(diag["q_target"]) else np.full_like(ref["grid"], np.nan)),
        "final_eff_counts": (diag["eff_counts"][-1] if len(diag["eff_counts"]) else np.full_like(ref["grid"], np.nan)),
        "F_target_ema": (diag["F_target_ema"] if diag["F_target_ema"] is not None else np.full_like(ref["grid"], np.nan)),
        "birth_hist": diag["birth_hist"], "death_hist": diag["death_hist"], "hist_edges": diag["hist_edges"],
    }
    return out


def save_run(path: str, out: dict):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, **out)
    os.replace(tmp, path)  # atomic -> a half-written file is never seen as valid


def load_run(path: str) -> dict:
    d = np.load(path, allow_pickle=True)
    return {k: d[k] for k in d.files}


# ---------------------------------------------------------------------------
# Stage -> RunSpec expansion.
# ---------------------------------------------------------------------------
def effective_base(cfg: dict, stage: str) -> dict:
    """Base SimConfig block with this stage's `base_overrides` applied."""
    base = dict(cfg.get("base", {}))
    st = cfg.get("stages", {}).get(stage, {})
    base.update(st.get("base_overrides", {}))
    return base


def _method_knobs(cfg: dict, name: str) -> dict:
    """Return (type, knobs) for a named method, falling back to base defaults.

    knobs holds fr_rate, target_ema_rate, max_event_fraction, fr_every,
    fr_start_steps, score_clip.
    """
    base = cfg.get("base", {})
    m = cfg["methods"][name]
    defaults = dict(
        fr_rate=m.get("fr_rate", 0.0),
        target_ema_rate=m.get("target_ema_rate", base.get("target_ema_rate", 0.005)),
        max_event_fraction=m.get("max_event_fraction", base.get("max_event_fraction", 0.02)),
        fr_every=m.get("fr_every", base.get("fr_every", 5)),
        fr_start_steps=m.get("fr_start_steps", base.get("fr_start_steps", 20000)),
        score_clip=m.get("score_clip", base.get("score_clip", 2.0)),
    )
    return m["type"], defaults


def _mk_spec(stage, name, mtype, knobs, seed, n_steps, n_replicas, a, save_every):
    return RunSpec(
        stage=stage, name=name, method=mtype, seed=int(seed), n_steps=int(n_steps),
        n_replicas=int(n_replicas), a=float(a), fr_rate=float(knobs["fr_rate"]),
        target_ema_rate=float(knobs["target_ema_rate"]),
        max_event_fraction=float(knobs["max_event_fraction"]),
        fr_every=int(knobs["fr_every"]), fr_start_steps=int(knobs["fr_start_steps"]),
        score_clip=float(knobs["score_clip"]), save_every=int(save_every))


def expand_stage(cfg: dict, stage: str) -> list:
    """Expand one stage's YAML definition into a deduplicated list of RunSpecs.

    Stage-level `knob_overrides` override per-method FR knobs (e.g. fr_start_steps
    for tiny smoke runs); `base_overrides` (handled by effective_base) override the
    SimConfig base block. Both keep tiny stages self-consistent.
    """
    st = cfg["stages"][stage]
    base = effective_base(cfg, stage)
    save_every = int(base.get("save_every", 2500))
    knob_ov = dict(st.get("knob_overrides", {}))
    a_default = float(cfg.get("system", {}).get("a", 1.5))
    seeds = list(st.get("seeds", [42]))
    specs = []

    # axes that can be vectorised in a stage
    n_steps_list = st.get("n_steps_values", [st.get("n_steps", 250000)])
    n_rep_list = st.get("n_replicas_values", [st.get("n_replicas", 1024)])
    a_list = st.get("a_values", [a_default])

    base_methods = st.get("methods", [])
    for seed in seeds:
        for n_steps in n_steps_list:
            for n_rep in n_rep_list:
                for a in a_list:
                    for mname in base_methods:
                        mtype, knobs = _method_knobs(cfg, mname)
                        knobs.update(knob_ov)
                        specs.append(_mk_spec(stage, mname, mtype, knobs, seed,
                                              n_steps, n_rep, a, save_every))
                    # one-axis sweeps off fr_est_tuned (pilot / failure stages)
                    for sweep in st.get("sweeps", []):
                        axis = sweep["axis"]
                        mtype, knobs0 = _method_knobs(cfg, "fr_est_tuned")
                        knobs0.update(knob_ov)
                        for val in sweep["values"]:
                            knobs = dict(knobs0)
                            knobs[axis] = val
                            label = f"sweep_{axis}_{val:g}" if isinstance(val, float) else f"sweep_{axis}_{val}"
                            specs.append(_mk_spec(stage, label, mtype, knobs, seed,
                                                  n_steps, n_rep, a, save_every))
    # dedupe by run_id (a sweep may reproduce the tuned point)
    seen, uniq = set(), []
    for s in specs:
        rid = s.run_id()
        if rid not in seen:
            seen.add(rid)
            uniq.append(s)
    return uniq


def load_yaml(path: str) -> dict:
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh)
