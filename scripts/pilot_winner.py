#!/usr/bin/env python3
"""Summarize the Stage-A pilot: matched-seed median gain per swept config.

For each pilot config (sweep_<axis>_<value>), reports the matched-seed median and
IQR of % gain in final L2(F) vs the abf baseline (matched on seed), plus median
final L2(F), ancestor ESS, and replacement count. Helps lock the production FR
operating regime. Reads results/wca_production/raw/pilot__*.npz.

Usage: python scripts/pilot_winner.py
"""
from __future__ import annotations

import glob
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

import numpy as np  # noqa: E402
import wca_jobs  # noqa: E402


def val(d, k):
    v = d[k]
    return v.item() if hasattr(v, "ndim") and v.ndim == 0 else v


def main():
    raw = "results/wca_production/raw"
    rows = []
    for p in sorted(glob.glob(os.path.join(raw, "pilot__*.npz"))):
        d = wca_jobs.load_run(p)
        rows.append(dict(name=str(val(d, "name")), method=str(val(d, "method")),
                         seed=int(val(d, "seed")), l2_f=float(val(d, "l2_f")),
                         ess=float(val(d, "final_ancestor_ess")),
                         repl=int(val(d, "total_replacement_events")),
                         fr_rate=float(val(d, "fr_rate")),
                         ema=float(val(d, "target_ema_rate")),
                         mef=float(val(d, "max_event_fraction")),
                         fr_every=int(val(d, "fr_every")),
                         fr_start=int(val(d, "fr_start_steps"))))
    if not rows:
        print("no pilot runs found yet")
        return 1
    abf = {r["seed"]: r["l2_f"] for r in rows if r["method"] == "abf"}
    print(f"ABF baselines by seed: " + ", ".join(f"{s}:{v:.4f}" for s, v in sorted(abf.items())))
    groups = defaultdict(list)
    for r in rows:
        if r["method"] == "abf":
            continue
        groups[r["name"]].append(r)
    out = []
    for name, rs in groups.items():
        gains = [100.0 * (abf[r["seed"]] - r["l2_f"]) / abf[r["seed"]]
                 for r in rs if r["seed"] in abf and abf[r["seed"]] > 0]
        wins = sum(1 for r in rs if r["seed"] in abf and r["l2_f"] < abf[r["seed"]])
        ex = rs[0]
        out.append(dict(name=name, n=len(rs), wins=wins,
                        med_gain=float(np.median(gains)) if gains else float("nan"),
                        med_l2f=float(np.median([r["l2_f"] for r in rs])),
                        med_ess=float(np.median([r["ess"] for r in rs])),
                        med_repl=float(np.median([r["repl"] for r in rs])),
                        rate=ex["fr_rate"], ema=ex["ema"], mef=ex["mef"],
                        every=ex["fr_every"], start=ex["fr_start"]))
    out.sort(key=lambda r: -r["med_gain"])
    print(f"\n{'config':28s}{'n':>3}{'wins':>5}{'medGain%':>9}{'medL2F':>9}"
          f"{'medESS':>8}{'medRepl':>9}")
    for r in out:
        print(f"{r['name']:28s}{r['n']:>3}{r['wins']:>5}{r['med_gain']:>+9.1f}"
              f"{r['med_l2f']:>9.4f}{r['med_ess']:>8.0f}{r['med_repl']:>9.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
