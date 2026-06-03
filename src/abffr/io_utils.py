"""I/O, configuration and reproducibility helpers for the ABF--FR study.

All experiment outputs live under ``<output_root>`` (default
``results/two_dim_xi_x``) in fixed sub-directories so that scripts can find each
other's outputs without extra wiring.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List

import numpy as np

try:
    import yaml
except Exception:  # pragma: no cover - yaml is a hard dependency in practice
    yaml = None


# Mapping from a CLI ``--stage`` to the results sub-directory it reads/writes.
# ``smoke`` deliberately writes into the ``tuning`` directory so that a tiny
# smoke run can be inspected with the same ``--stage tuning`` plotting/table
# commands as a real tuning run.
STAGE_TO_DIR = {
    "smoke": "tuning",
    "tuning": "tuning",
    "eval": "eval",
    # GPU/torch backend stages (kept separate from the CPU outputs so the CPU
    # reference results are never overwritten).  ``smoke_gpu`` writes into the
    # ``tuning_gpu`` directory, mirroring how ``smoke`` writes into ``tuning``.
    "smoke_gpu": "tuning_gpu",
    "tuning_gpu": "tuning_gpu",
    "production_gpu": "production_gpu",
}
STAGE_TO_FIGDIR = {
    "smoke": "figures_tuning",
    "tuning": "figures_tuning",
    "eval": "figures_eval",
    "smoke_gpu": "figures_tuning_gpu",
    "tuning_gpu": "figures_tuning_gpu",
    "production_gpu": "figures_production_gpu",
}


def load_config(path: str) -> Dict[str, Any]:
    """Load a YAML config file into a plain dict."""
    if yaml is None:
        raise RuntimeError(
            "pyyaml is required to load configs. Install with `pip install pyyaml`."
        )
    with open(path, "r") as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise ValueError(f"Config {path!r} did not parse to a mapping.")
    return cfg


def ensure_dir(path: str) -> str:
    """Create ``path`` (and parents) if needed; return it."""
    os.makedirs(path, exist_ok=True)
    return path


def output_root(cfg: Dict[str, Any]) -> str:
    return cfg.get("output_root", "results/two_dim_xi_x")


def reference_dir(cfg: Dict[str, Any]) -> str:
    return ensure_dir(os.path.join(output_root(cfg), "reference"))


def stage_dir(cfg: Dict[str, Any], stage: str) -> str:
    """Results directory for a run stage (``smoke``/``tuning``/``eval``)."""
    sub = STAGE_TO_DIR.get(stage)
    if sub is None:
        raise ValueError(f"Unknown stage {stage!r}; expected one of {list(STAGE_TO_DIR)}")
    return ensure_dir(os.path.join(output_root(cfg), sub))


def figure_dir(cfg: Dict[str, Any], stage: str) -> str:
    sub = STAGE_TO_FIGDIR.get(stage)
    if sub is None:
        raise ValueError(f"Unknown stage {stage!r}; expected one of {list(STAGE_TO_FIGDIR)}")
    return ensure_dir(os.path.join(output_root(cfg), sub))


def stage_prefix(stage: str) -> str:
    """File-name prefix used by a stage's CSV outputs (``tuning``/``eval``)."""
    return STAGE_TO_DIR.get(stage, stage)


def apply_cli_overrides(cfg: Dict[str, Any], *, device=None, dtype=None,
                        batch_size_configs=None, n_steps=None, n_particles=None,
                        eval_every=None, seeds=None, output_root=None,
                        estimator=None) -> Dict[str, Any]:
    """Apply optional CLI overrides onto a loaded config (in place).

    Supports the GPU-script overrides the study spec asks for
    (``--device``/``--batch-size-configs``/``--n-steps``/``--n-particles`` and a
    few more) without editing YAML files.
    """
    if output_root is not None:
        cfg["output_root"] = output_root
    if device is not None:
        cfg["device"] = device
    if dtype is not None:
        cfg["dtype"] = dtype
    if batch_size_configs is not None:
        cfg["batch_size_configs"] = int(batch_size_configs)
    if estimator is not None:
        cfg.setdefault("abf", {})["estimator"] = estimator
    sim = cfg.setdefault("simulation", {})
    if n_steps is not None:
        sim["n_steps"] = int(n_steps)
    if n_particles is not None:
        sim["n_particles"] = int(n_particles)
    if eval_every is not None:
        sim["eval_every"] = int(eval_every)
    if seeds is not None:
        sim["seeds"] = [int(s) for s in seeds]
    return cfg


def make_rng_streams(seed: int):
    """Return three independent, reproducible RNG streams for one run.

    Splitting the seed into independent streams for initialisation, Langevin
    noise and Fisher--Rao resampling means that, for a fixed ``seed``, the
    initial conditions *and* the Langevin noise realisation are identical across
    methods regardless of whether Fisher--Rao is active.  Only the birth--death
    resampling consumes the third (separate) stream, giving clean matched-seed
    comparisons between ABF-only and ABF+FR.
    """
    ss = np.random.SeedSequence(int(seed))
    child_init, child_noise, child_fr = ss.spawn(3)
    return (
        np.random.default_rng(child_init),
        np.random.default_rng(child_noise),
        np.random.default_rng(child_fr),
    )


def save_json(path: str, obj: Dict[str, Any]) -> None:
    ensure_dir(os.path.dirname(path) or ".")
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=_json_default)


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


@dataclass
class RunSpec:
    """A single simulation run (one method/target/hyperparameter/seed point)."""

    method: str
    target_type: str          # "none", "estimated", "uniform", "oracle", "self"
    seed: int
    gamma: float = 0.0
    eta: float = 0.10
    burnin_fraction: float = 0.0
    fr_every: int = 5
    extra: Dict[str, Any] = field(default_factory=dict)

    @property
    def config_id(self) -> str:
        """Identifier shared by all seeds of the same hyperparameter point."""
        return (f"{self.method}|tt={self.target_type}|g={self.gamma:g}"
                f"|eta={self.eta:g}|bi={self.burnin_fraction:g}|fe={self.fr_every}")

    @property
    def run_id(self) -> str:
        return f"{self.config_id}|seed={self.seed}"

    def to_row(self) -> Dict[str, Any]:
        return dict(
            run_id=self.run_id,
            config_id=self.config_id,
            method=self.method,
            target_type=self.target_type,
            seed=int(self.seed),
            gamma=float(self.gamma),
            eta=float(self.eta),
            burnin_fraction=float(self.burnin_fraction),
            fr_every=int(self.fr_every),
        )


def build_run_specs(cfg: Dict[str, Any], seeds: List[int]) -> List[RunSpec]:
    """Expand a config into the list of :class:`RunSpec` to simulate.

    ``abf_only`` is run once per seed (no FR hyperparameters).  Each
    ``abf_fr_*`` method is crossed with the FR hyperparameter grid.  The target
    type is inferred from the method name and validated against
    ``fr.target_types``.
    """
    method_to_target = {
        "abf_only": "none",
        "abf_fr_estimated": "estimated",
        "abf_fr_uniform": "uniform",
        "abf_fr_oracle": "oracle",
        "abf_fr_self": "self",
    }
    fr = cfg.get("fr", {})
    allowed_targets = set(fr.get("target_types", ["estimated", "uniform", "oracle"]))
    gamma_values = list(fr.get("gamma_values", [0.02]))
    eta_values = list(fr.get("eta_values", [cfg.get("abf", {}).get("eta", 0.10)]))
    burnin_fractions = list(fr.get("burnin_fractions", [0.0]))
    fr_every_values = list(fr.get("fr_every_values", [5]))

    # ABF-only uses a single eta for its diagnostic KDE; pick the first listed.
    default_eta = eta_values[0] if eta_values else 0.10

    specs: List[RunSpec] = []
    for method in cfg.get("methods", ["abf_only"]):
        target = method_to_target.get(method)
        if target is None:
            raise ValueError(f"Unknown method {method!r} in config.")
        if method == "abf_only":
            for seed in seeds:
                specs.append(RunSpec(method=method, target_type="none", seed=seed,
                                     gamma=0.0, eta=default_eta,
                                     burnin_fraction=0.0, fr_every=1))
            continue
        if target not in allowed_targets:
            # Method requested but its target type was disabled in fr.target_types.
            continue
        for gamma in gamma_values:
            for eta in eta_values:
                for burnin in burnin_fractions:
                    for fr_every in fr_every_values:
                        for seed in seeds:
                            specs.append(RunSpec(
                                method=method, target_type=target, seed=seed,
                                gamma=float(gamma), eta=float(eta),
                                burnin_fraction=float(burnin), fr_every=int(fr_every),
                            ))
    return specs
