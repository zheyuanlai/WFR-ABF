"""Batching, checkpoint/resume and CSV assembly for the GPU backend.

This module turns a list of :class:`abffr.io_utils.RunSpec` into batched calls
to :func:`abffr.simulation_torch.run_batch` and writes outputs whose schema is a
*superset* of the CPU runner's (so the existing plotting/table scripts and the
study-spec columns are both satisfied).  Metric and conditional-diagnostic rows
are produced with the *same* :mod:`abffr.metrics` / :mod:`abffr.diagnostics`
code the CPU runner uses.

Resumption
----------
Each run has a unique ``run_id``.  A finished run drops a marker
``<stage>/completed/<safe_run_id>.done``; a crashed run drops
``<stage>/failed/<safe_run_id>.json`` with the error and config.  On restart,
completed run_ids (markers *and* any rows already in the final-summary CSV) are
skipped unless ``force=True``.  Each process writes tag-suffixed CSVs
(``<prefix>_<kind>__<tag>.csv``); :func:`merge_stage_csvs` concatenates the tags
into the canonical ``<prefix>_<kind>.csv`` (deduped by ``run_id``) for the
plotting/table scripts.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import re
import time
import traceback
from collections import defaultdict
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import torch

from . import (diagnostics, io_utils, metrics, reference, simulation_torch,
               torch_utils as tu)
from .io_utils import RunSpec

CSV_KINDS = ("runs_long", "final_summary", "profiles", "fr_events",
             "conditional_diagnostics")

# Columns that uniquely identify a row of each CSV kind, used to de-duplicate
# when merging tag-suffixed shard CSVs (a profile has one row per (run, x); a
# long row one per (run, snapshot); etc.).  ``keep="last"`` lets a forced re-run
# override an earlier row.
DEDUP_SUBSET = {
    "runs_long": ["run_id", "step"],
    "final_summary": ["run_id"],
    "profiles": ["run_id", "x"],
    "fr_events": ["run_id", "step"],
    # conditional_diagnostics rows carry only (method, target_type, seed) as
    # identity (not the full config), so several configs share a key -- like the
    # CPU runner we keep them all (plotting averages over them); no dedup.
}


# --------------------------------------------------------------------------- #
# Stage setup (device/dtype/reference/eval), shared by all GPU scripts
# --------------------------------------------------------------------------- #
def prepare_stage(cfg: Dict, stage: str, require_csv: bool = True, logger=print
                  ) -> Dict:
    """Resolve device/dtype, the stage output dir, the reference (gated) and the
    evaluation window for ``stage``."""
    device = tu.resolve_device(cfg.get("device"))
    dtype = tu.resolve_dtype(cfg.get("dtype"))
    stage_root = io_utils.stage_dir(cfg, stage)
    prefix = io_utils.stage_prefix(stage)
    x_grid, ref, csv_path = reference.load_reference_for_run(
        cfg, require_csv=require_csv, logger=logger)
    ev = metrics.EvalConfig.from_domain(cfg["domain"])
    return dict(device=device, dtype=dtype, stage_root=stage_root, prefix=prefix,
                x_grid=x_grid, ref=ref, ev=ev, csv_path=csv_path)


# --------------------------------------------------------------------------- #
# Batching
# --------------------------------------------------------------------------- #
def batch_key(spec: RunSpec):
    """Runs are batched within a common (target_type, eta, fr_every, burnin).

    These four must be constant inside a batch so the FR firing schedule and the
    smoothing bandwidth are shared; ``gamma`` and ``seed`` vary per row.
    """
    return (spec.target_type, float(spec.eta), int(spec.fr_every),
            float(spec.burnin_fraction))


def build_batches(specs: List[RunSpec], batch_size: int) -> List[List[RunSpec]]:
    groups: Dict = defaultdict(list)
    for s in specs:
        groups[batch_key(s)].append(s)
    batches: List[List[RunSpec]] = []
    for key in sorted(groups):
        gs = sorted(groups[key], key=lambda s: (float(s.gamma), int(s.seed)))
        for i in range(0, len(gs), batch_size):
            batches.append(gs[i:i + batch_size])
    return batches


# --------------------------------------------------------------------------- #
# Resume bookkeeping
# --------------------------------------------------------------------------- #
def _safe_run_id(run_id: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]", "_", run_id)
    return f"{s[:120]}__{hashlib.sha1(run_id.encode()).hexdigest()[:8]}"


def completed_dir(stage_root: str) -> str:
    return io_utils.ensure_dir(os.path.join(stage_root, "completed"))


def failed_dir(stage_root: str) -> str:
    return io_utils.ensure_dir(os.path.join(stage_root, "failed"))


def load_completed(stage_root: str, prefix: str) -> set:
    """Set of run_ids already finished (markers + rows in any final CSV)."""
    done = set()
    cdir = os.path.join(stage_root, "completed")
    if os.path.isdir(cdir):
        for p in glob.glob(os.path.join(cdir, "*.done")):
            try:
                with open(p) as fh:
                    done.add(json.load(fh)["run_id"])
            except Exception:
                pass
    for path in glob.glob(os.path.join(stage_root, f"{prefix}_final_summary*.csv")):
        try:
            df = pd.read_csv(path, usecols=["run_id"])
            done.update(df["run_id"].astype(str).tolist())
        except Exception:
            pass
    return done


def _write_marker(stage_root: str, spec: RunSpec, summary: Dict) -> None:
    path = os.path.join(completed_dir(stage_root), _safe_run_id(spec.run_id) + ".done")
    io_utils.save_json(path, dict(run_id=spec.run_id, config_id=spec.config_id,
                                  final_l2_F=summary.get("final_l2_F"),
                                  final_l2_Fprime=summary.get("final_l2_Fprime")))


def _write_failure(stage_root: str, spec: RunSpec, err: str) -> None:
    path = os.path.join(failed_dir(stage_root), _safe_run_id(spec.run_id) + ".json")
    io_utils.save_json(path, dict(run_id=spec.run_id, **spec.to_row(), error=err))


# --------------------------------------------------------------------------- #
# Row assembly (reuses metrics.py / diagnostics.py, like the CPU runner)
# --------------------------------------------------------------------------- #
def _rows_for_run(spec: RunSpec, diag: Dict, cfg: Dict, x_grid, ref, ev,
                  runtime_seconds: float, conditional: str):
    meta = spec.to_row()
    h = float(cfg["abf"]["h"])
    ramp_fraction = float(cfg.get("fr", {}).get("ramp_fraction", 0.1))
    beta = float(cfg["simulation"]["beta"])

    ts = metrics.time_series_metrics(diag, x_grid, ref["F_ref"], ref["Fprime_ref"], ev)
    long_rows, fr_rows = [], []
    for k, r in enumerate(ts):
        rf = metrics.region_fractions(diag["X_snap"][k], ev)
        long_rows.append({
            **meta, **r, "h": h, "ramp_fraction": ramp_fraction,
            # Study-spec aliases for the same quantities.
            "x_l2_to_target": r["marginal_l2_target"],
            "x_l2_to_uniform": r["marginal_l2_uniform"],
            "fr_score_std": r["score_std"], "fr_score_max": r["score_max"],
            "left_frac": rf["frac_left"], "barrier_frac": rf["frac_barrier"],
            "right_frac": rf["frac_right"],
        })
        fr_rows.append({
            **meta,
            "step": r["step"], "t": r["t"], "gamma_eff": r["gamma_eff"],
            "fr_applied": r["fr_applied"],
            "fr_event_fraction": r["fr_event_fraction"],
            "fr_event_fraction_max": r["fr_event_fraction_max"],
            "fr_events_total": r["fr_events_total"],
            "num_deaths": int(diag["fr_events_total"][k]),
            "num_births": int(diag["fr_events_total"][k]),
            "event_fraction": r["fr_event_fraction"],
            "score_mean": r["score_mean"], "score_std": r["score_std"],
            "score_min": r["score_min"], "score_max": r["score_max"],
            "target_l2": float(diag["target_l2"][k]),
            "n_unique_ancestors": r["n_unique_ancestors"],
            "ancestor_ess": float(diag["ancestor_ess"][k]),
            "max_clone_multiplicity": int(diag["max_clone_multiplicity"][k]),
        })

    summary = metrics.final_summary(diag, x_grid, ref["F_ref"], ref["Fprime_ref"], ev)
    final_row = {**meta, **summary, "runtime_seconds": float(runtime_seconds),
                 "total_barrier_crossings": int(summary["barrier_crossings"])}

    Fp = diag["Fprime_hat"][-1]; F = diag["F_hat"][-1]
    p = diag["p_hat_grid"][-1]; q = diag["q_target_grid"][-1]
    profile_rows = []
    for j in range(len(x_grid)):
        profile_rows.append({
            **meta, "x": float(x_grid[j]),
            "F_ref": float(ref["F_ref"][j]), "Fprime_ref": float(ref["Fprime_ref"][j]),
            "Fprime_hat": float(Fp[j]), "F_hat": float(F[j]),
            "p_hat": float(p[j]), "q_target": float(q[j]),
            "p_hat_x": float(p[j]), "target_x": float(q[j]),
        })

    cmeta = dict(method=spec.method, target_type=spec.target_type, seed=int(spec.seed))
    snap_idx = None if conditional == "final" else list(range(len(diag["steps"])))
    cond_rows = diagnostics.conditional_diagnostics(
        diag, cmeta, beta, cfg["domain"], snapshot_indices=snap_idx)

    return dict(runs_long=long_rows, final_summary=[final_row],
                profiles=profile_rows, fr_events=fr_rows,
                conditional_diagnostics=cond_rows), summary


# --------------------------------------------------------------------------- #
# Tag-suffixed CSV writer (process-local, append + flush)
# --------------------------------------------------------------------------- #
class _CsvBuffer:
    def __init__(self, stage_root: str, prefix: str, tag: str,
                 seed_from_disk: bool = True):
        self.stage_root, self.prefix, self.tag = stage_root, prefix, tag
        self.rows = {k: [] for k in CSV_KINDS}
        # Seed from this tag's prior partial output (resume within a process).
        # Skipped on --force so a forced re-run overwrites cleanly rather than
        # appending duplicate rows (conditional_diagnostics is not deduped).
        if not seed_from_disk:
            return
        for k in CSV_KINDS:
            path = self._path(k)
            if os.path.exists(path):
                try:
                    self.rows[k] = pd.read_csv(path).to_dict("records")
                except Exception:
                    self.rows[k] = []

    def _path(self, kind: str) -> str:
        return os.path.join(self.stage_root, f"{self.prefix}_{kind}__{self.tag}.csv")

    def extend(self, payload: Dict[str, list]) -> None:
        for k, rows in payload.items():
            self.rows[k].extend(rows)

    def flush(self) -> None:
        for k in CSV_KINDS:
            if self.rows[k]:
                pd.DataFrame(self.rows[k]).to_csv(self._path(k), index=False)


# --------------------------------------------------------------------------- #
# Main driver
# --------------------------------------------------------------------------- #
def run_specs(
    specs: List[RunSpec],
    *,
    cfg: Dict,
    stage_root: str,
    prefix: str,
    x_grid: np.ndarray,
    ref: Dict,
    ev,
    device: torch.device,
    dtype: torch.dtype,
    estimator: str = "binned_smooth",
    batch_size: int = 16,
    base_seed: int = 0,
    tag: str = "main",
    resume: bool = True,
    force: bool = False,
    conditional: str = "final",
    logger=print,
) -> Dict:
    """Run all ``specs`` on one ``device`` with checkpoint/resume.

    Returns a small summary dict (counts, runtimes, NaN count).
    """
    completed = load_completed(stage_root, prefix) if (resume and not force) else set()
    todo = [s for s in specs if s.run_id not in completed]
    skipped = len(specs) - len(todo)
    batches = build_batches(todo, batch_size)
    logger(f"[parallel] device={device} estimator={estimator} "
           f"runs={len(specs)} todo={len(todo)} skipped(resume)={skipped} "
           f"batches={len(batches)} batch_size<={batch_size} tag={tag}")

    buf = _CsvBuffer(stage_root, prefix, tag, seed_from_disk=not force)
    n_done, n_failed, n_nan = 0, 0, 0
    t_start = time.time()

    for bi, batch in enumerate(batches):
        try:
            res = simulation_torch.run_batch(
                batch, cfg=cfg, x_grid=x_grid, F_ref=ref["F_ref"],
                Fprime_ref=ref["Fprime_ref"], ev=ev, device=device, dtype=dtype,
                estimator=estimator, base_seed=base_seed)
        except Exception as exc:  # whole-batch failure
            err = f"{exc}\n{traceback.format_exc()}"
            for spec in batch:
                _write_failure(stage_root, spec, err)
            n_failed += len(batch)
            logger(f"[parallel] BATCH {bi+1}/{len(batches)} FAILED "
                   f"({len(batch)} runs): {exc}")
            continue

        per_run_runtime = res.runtime_seconds / max(len(batch), 1)
        for spec, diag in zip(batch, res.diags):
            try:
                payload, summary = _rows_for_run(
                    spec, diag, cfg, x_grid, ref, ev, per_run_runtime, conditional)
                buf.extend(payload)
                _write_marker(stage_root, spec, summary)
                n_done += 1
                n_nan += int(bool(summary.get("any_nan")))
            except Exception as exc:
                _write_failure(stage_root, spec, f"{exc}\n{traceback.format_exc()}")
                n_failed += 1
        buf.flush()
        logger(f"[parallel] batch {bi+1}/{len(batches)} done "
               f"({len(batch)} runs, {res.runtime_seconds:.1f}s, "
               f"{res.runtime_seconds/max(len(batch),1):.2f}s/run); "
               f"cumulative done={n_done} failed={n_failed} nan={n_nan}")

    buf.flush()
    return dict(n_runs=len(specs), n_done=n_done, n_failed=n_failed,
                n_skipped=skipped, n_nan=n_nan, wall_seconds=time.time() - t_start,
                tag=tag)


# --------------------------------------------------------------------------- #
# Merge + config aggregation (canonical CSVs for plotting / tables)
# --------------------------------------------------------------------------- #
def merge_stage_csvs(stage_root: str, prefix: str, logger=print) -> None:
    """Concatenate all tag-suffixed CSVs into canonical ``<prefix>_<kind>.csv``."""
    for kind in CSV_KINDS:
        parts = sorted(glob.glob(os.path.join(stage_root, f"{prefix}_{kind}__*.csv")))
        if not parts:
            continue
        frames = []
        for p in parts:
            try:
                frames.append(pd.read_csv(p))
            except Exception:
                pass
        if not frames:
            continue
        df = pd.concat(frames, ignore_index=True)
        subset = [c for c in DEDUP_SUBSET.get(kind, []) if c in df.columns]
        if subset:
            df = df.drop_duplicates(subset=subset, keep="last")
        out = os.path.join(stage_root, f"{prefix}_{kind}.csv")
        df.to_csv(out, index=False)
        logger(f"[parallel] merged {len(parts)} part(s) -> "
               f"{os.path.relpath(out)} ({len(df)} rows)")


def _iqr(s):
    s = np.asarray(s, dtype=float)
    s = s[np.isfinite(s)]
    return float(np.percentile(s, 75) - np.percentile(s, 25)) if s.size else float("nan")


def summarize_configs(final_df: pd.DataFrame) -> pd.DataFrame:
    """Per-config median/IQR over seeds (mirrors the CPU runner)."""
    rows = []
    keys = ["config_id", "method", "target_type", "gamma", "eta",
            "burnin_fraction", "fr_every"]
    for cid, g in final_df.groupby("config_id"):
        row = {k: g[k].iloc[0] for k in keys}
        row["n_seeds"] = len(g)
        for col in ["final_l2_F", "final_l2_Fprime", "integrated_l2_F",
                    "integrated_l2_Fprime", "final_marginal_l2_uniform",
                    "final_marginal_l2_target", "mean_fr_event_fraction",
                    "max_fr_event_fraction", "barrier_crossings",
                    "frac_left", "frac_barrier", "frac_right"]:
            if col in g:
                row[f"median_{col}"] = float(np.nanmedian(g[col]))
                row[f"iqr_{col}"] = _iqr(g[col])
        row["frac_nan"] = float(np.mean(g["any_nan"])) if "any_nan" in g else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def select_best_configs(config_df: pd.DataFrame, cfg: Dict) -> pd.DataFrame:
    """Best config per target type by median integrated L2(F) (mirror CPU)."""
    cap = float(cfg.get("fr", {}).get("max_event_fraction", 0.10))
    rows = []
    for target in ["none", "estimated", "uniform", "oracle", "self"]:
        sub = config_df[config_df["target_type"] == target].copy()
        if sub.empty:
            continue
        passed = (
            (sub.get("frac_nan", 0.0) == 0.0)
            & (sub.get("median_mean_fr_event_fraction", 0.0) <= cap + 1e-9)
            & (sub.get("median_max_fr_event_fraction", 0.0) <= 1.5 * cap + 1e-9)
        )
        sub["passed_safety"] = passed
        pool = sub[passed] if passed.any() else sub
        pool = pool.sort_values(
            ["median_integrated_l2_F", "median_final_l2_F"]).reset_index(drop=True)
        for rank, (_, r) in enumerate(pool.iterrows(), start=1):
            rows.append(dict(
                rank_within_target=rank, selected=(rank == 1),
                method=r["method"], target_type=r["target_type"],
                gamma=float(r["gamma"]), eta=float(r["eta"]),
                burnin_fraction=float(r["burnin_fraction"]), fr_every=int(r["fr_every"]),
                median_integrated_l2_F=float(r["median_integrated_l2_F"]),
                median_final_l2_F=float(r["median_final_l2_F"]),
                median_final_l2_Fprime=float(r["median_final_l2_Fprime"]),
                median_mean_fr_event_fraction=float(r["median_mean_fr_event_fraction"]),
                passed_safety=bool(r["passed_safety"]), n_seeds=int(r["n_seeds"]),
            ))
    return pd.DataFrame(rows)


def write_config_summaries(stage_root: str, prefix: str, cfg: Dict, logger=print):
    """Write ``<prefix>_config_summary.csv`` and ``best_configs.csv`` if possible."""
    path = os.path.join(stage_root, f"{prefix}_final_summary.csv")
    if not os.path.exists(path):
        logger(f"[parallel] {os.path.relpath(path)} missing; skipping config summary.")
        return
    final_df = pd.read_csv(path)
    if final_df.empty:
        return
    config_df = summarize_configs(final_df)
    config_df.to_csv(os.path.join(stage_root, f"{prefix}_config_summary.csv"), index=False)
    best_df = select_best_configs(config_df, cfg)
    best_df.to_csv(os.path.join(stage_root, "best_configs.csv"), index=False)
    logger(f"[parallel] wrote {prefix}_config_summary.csv ({len(config_df)} configs) "
           f"and best_configs.csv "
           f"({int(best_df['selected'].sum()) if len(best_df) else 0} selected)")
