"""Reference free energy / mean force for ``xi(x, y) = x`` by quadrature.

Definitions (free energy defined up to an additive constant ``C``):

    F_ref(x)      = -(1/beta) * log integral_y exp(-beta V(x, y)) dy + C
    F'_ref(x)     =  integral_y dV/dx(x, y) exp(-beta V(x, y)) dy
                     -----------------------------------------------------
                            integral_y          exp(-beta V(x, y)) dy
    p_ref(x)      proportional to integral_y exp(-beta V(x, y)) dy   (unbiased x-marginal)

The y-integral is evaluated by the trapezoidal rule on a fine grid with a
log-sum-exp shift for numerical stability.
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from . import potentials

EPS = 1e-300


def compute_reference(x_grid, y_grid, beta,
                      V_func=potentials.potential_xy,
                      dVdx_func=potentials.dVdx_xy):
    """Compute reference profiles on ``x_grid`` using y-quadrature on ``y_grid``.

    Returns a dict with keys ``log_Z``, ``F_ref``, ``Fprime_ref``, ``p_ref``,
    all 1-D arrays on ``x_grid``.  ``F_ref`` is centred to have zero mean over
    ``x_grid`` (the additive constant is fixed by the convention
    ``F <- F - mean(F)``); ``p_ref`` integrates to 1 over ``x_grid``.
    """
    x_grid = np.asarray(x_grid, dtype=float)
    y_grid = np.asarray(y_grid, dtype=float)

    # Shape (n_x, n_y): each row is a fixed x swept over the y quadrature nodes.
    xx = x_grid[:, None]
    yy = y_grid[None, :]

    phi = beta * V_func(xx, yy)               # beta V(x, y)
    dvdx = dVdx_func(xx, yy)                   # dV/dx(x, y)

    m = phi.min(axis=1, keepdims=True)        # per-x shift for stability
    w = np.exp(-(phi - m))                     # exp(-(beta V - m)) in [0, 1]

    Z_stab = np.trapezoid(w, y_grid, axis=1)   # integral of shifted weights
    log_Z = -m[:, 0] + np.log(np.maximum(Z_stab, EPS))

    Fprime_ref = (np.trapezoid(dvdx * w, y_grid, axis=1)
                  / np.maximum(Z_stab, EPS))

    F_ref = -(1.0 / beta) * log_Z
    F_ref = F_ref - np.mean(F_ref)             # F <- F - mean(F)

    # Unbiased x-marginal p_ref(x) proportional to Z(x) = exp(log_Z).
    lz = log_Z - log_Z.max()
    p_unnorm = np.exp(lz)
    p_ref = p_unnorm / np.maximum(np.trapezoid(p_unnorm, x_grid), EPS)

    return dict(log_Z=log_Z, F_ref=F_ref, Fprime_ref=Fprime_ref, p_ref=p_ref)


def conditional_y_density(x0, y_grid, beta, V_func=potentials.potential_xy):
    """Reference conditional density ``p_ref(y | x = x0)`` on ``y_grid``.

    Normalised to integrate to 1 over ``y_grid``.
    """
    y_grid = np.asarray(y_grid, dtype=float)
    phi = beta * V_func(np.full_like(y_grid, float(x0)), y_grid)
    phi = phi - phi.min()
    w = np.exp(-phi)
    Z = np.maximum(np.trapezoid(w, y_grid), EPS)
    return w / Z


def build_reference_grid(cfg: Dict, beta: float):
    """Build the 2-D reference grid (potential and Boltzmann density).

    Returns ``(x_grid, y_grid, V_grid, rho_grid)`` where ``V_grid`` and
    ``rho_grid`` have shape ``(ny, nx)`` (``indexing="xy"``) and ``rho_grid``
    integrates to 1 over the 2-D domain.
    """
    d = cfg["domain"]
    x_grid, y_grid, XX, YY = potentials.make_grid(
        d["x_min"], d["x_max"], d["y_min"], d["y_max"],
        d["nx_ref"], d["ny_ref"],
    )
    V_grid = potentials.potential_xy(XX, YY)
    phi = beta * (V_grid - V_grid.min())
    rho = np.exp(-phi)
    norm = np.trapezoid(np.trapezoid(rho, x_grid, axis=1), y_grid)
    rho_grid = rho / np.maximum(norm, EPS)
    return x_grid, y_grid, V_grid, rho_grid


def profile_grid(cfg: Dict):
    """1-D reaction-coordinate grid used for ABF profiles and FR targets."""
    d = cfg["domain"]
    return np.linspace(d["x_min"], d["x_max"], int(d["nx_profile"]))


def load_reference_for_run(cfg: Dict, require_csv: bool = True, logger=print):
    """Reference profiles on the simulation grid, with a fail-loud CSV gate.

    Like the CPU runner, the reference used for metrics and the oracle FR target
    is *recomputed on the simulation x-grid* by y-quadrature (grid-exact).  In
    addition -- per the GPU safety rules -- ``reference/reference_profile.csv``
    must exist (raise :class:`FileNotFoundError` otherwise); when present it is
    interpolated onto the grid and cross-checked, warning loudly on a mismatch
    (e.g. a stale reference computed for a different beta/domain).

    Returns ``(x_grid, ref_dict, csv_path)``.
    """
    import os

    from . import io_utils

    beta = float(cfg["simulation"]["beta"])
    x_grid = profile_grid(cfg)
    d = cfg["domain"]
    ny = int(d.get("ny_ref", 801))
    y_quad = np.linspace(d["y_min"], d["y_max"], ny)
    ref = compute_reference(x_grid, y_quad, beta)

    csv_path = os.path.join(io_utils.reference_dir(cfg), "reference_profile.csv")
    if not os.path.exists(csv_path):
        if require_csv:
            raise FileNotFoundError(
                f"Reference file {csv_path!r} not found. Run\n"
                f"    python scripts/run_reference_2d.py --config <this-config>\n"
                f"first (the GPU backend requires the reference to exist before "
                f"production).")
        logger(f"[reference] NOTE: {os.path.relpath(csv_path)} not found; using "
               f"the on-grid quadrature reference only.")
        return x_grid, ref, csv_path

    try:
        import pandas as pd
        df = pd.read_csv(csv_path)
        F_csv = np.interp(x_grid, df["x"].values, df["F_ref"].values)
        Fp_csv = np.interp(x_grid, df["x"].values, df["Fprime_ref"].values)
        dF = float(np.max(np.abs((F_csv - F_csv.mean())
                                 - (ref["F_ref"] - ref["F_ref"].mean()))))
        dFp = float(np.max(np.abs(Fp_csv - ref["Fprime_ref"])))
        if dF > 5e-2 or dFp > 5e-2:
            logger(f"[reference] WARNING: on-disk reference_profile.csv differs "
                   f"from the on-grid reference (max |dF|={dF:.3g}, "
                   f"max |dF'|={dFp:.3g}). This is expected for a different grid "
                   f"size but a large mismatch may mean a stale reference "
                   f"(different beta/domain). Re-run run_reference_2d.py.")
    except Exception as exc:  # pragma: no cover - cross-check is best-effort
        logger(f"[reference] WARNING: could not cross-check reference CSV: {exc}")
    return x_grid, ref, csv_path
