"""Conditional-sampling diagnostics: does FR preserve ``Y | X = x``?

ABF accuracy hinges on sampling the correct conditional distribution
``p(y | x = x*)``, not merely on flattening the x-marginal.  For a set of x-bin
centres we compare the empirical conditional ``p_hat(y | x in I_j)`` (the y
histogram of particles whose x falls in a narrow window) to the reference
conditional ``p_ref(y | x = centre) proportional to exp(-beta V(centre, y))``.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np

from . import reference

DEFAULT_X_BIN_CENTERS = (-1.5, -0.5, 0.0, 0.5, 1.5)
EPS = 1e-30


def _empirical_y_density(Y_sel: np.ndarray, edges: np.ndarray):
    """Normalised histogram density of ``Y_sel`` on bin centres of ``edges``."""
    counts, _ = np.histogram(Y_sel, bins=edges, density=False)
    centers = 0.5 * (edges[:-1] + edges[1:])
    dy = edges[1] - edges[0]
    total = counts.sum()
    if total == 0:
        return centers, np.zeros_like(centers, dtype=float)
    dens = counts / (total * dy)
    return centers, dens


def conditional_diagnostics(
    diag: Dict,
    row_meta: Dict,
    beta: float,
    domain: Dict[str, float],
    x_bin_centers: Sequence[float] = DEFAULT_X_BIN_CENTERS,
    bin_half_width: float = 0.25,
    n_y_bins: int = 41,
    min_samples: int = 10,
    snapshot_indices: Optional[Sequence[int]] = None,
) -> List[Dict]:
    """Conditional Y|X diagnostics for one run (long format rows).

    ``row_meta`` supplies the identifying columns (method, target_type, seed).
    ``snapshot_indices`` selects which saved snapshots to evaluate (default: the
    final snapshot only, to keep output compact).
    """
    edges = np.linspace(domain["y_min"], domain["y_max"], int(n_y_bins) + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    dy = edges[1] - edges[0]

    # Reference conditionals are time-independent: precompute once per x-bin.
    ref_cond = {}
    for c in x_bin_centers:
        pref = reference.conditional_y_density(c, centers, beta)
        pref = pref / max(np.sum(pref) * dy, EPS)  # normalise on the bin centres
        ref_cond[c] = pref

    n_snap = len(diag["steps"])
    if snapshot_indices is None:
        snapshot_indices = [n_snap - 1]

    rows: List[Dict] = []
    for k in snapshot_indices:
        if k < 0:
            k = n_snap + k
        if k < 0 or k >= n_snap:
            continue
        X = diag["X_snap"][k]
        Y = diag["Y_snap"][k]
        t = float(diag["times"][k])
        for c in x_bin_centers:
            sel = np.abs(X - c) <= bin_half_width
            n_sel = int(np.count_nonzero(sel))
            if n_sel < min_samples:
                l2, kl = np.nan, np.nan
            else:
                _, p_emp = _empirical_y_density(Y[sel], edges)
                p_ref = ref_cond[c]
                l2 = float(np.sqrt(np.sum((p_emp - p_ref) ** 2) * dy))
                # KL(p_emp || p_ref) with clipping for empty empirical bins.
                mask = p_emp > 0
                kl = float(np.sum(p_emp[mask] * np.log(
                    (p_emp[mask] + EPS) / (p_ref[mask] + EPS))) * dy)
            # Column order matches the study spec.
            row = dict(row_meta)
            row.update(t=t, x_bin_center=float(c), conditional_l2_y=l2,
                       conditional_kl_y=kl, n_samples_in_bin=n_sel)
            rows.append(row)
    return rows
