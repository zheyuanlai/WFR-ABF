"""2D model potential for the ABF--Fisher-Rao study.

The reaction coordinate is ``xi(x, y) = x`` throughout this package, so the
local mean force along the reaction coordinate is simply ``f(x, y) = dV/dx``.

Potential
---------
    V(x, y) =  3 * exp(-x^2) * ( exp(-(y - 1/3)^2) - exp(-(y - 5/3)^2) )
             - 5 * exp(-y^2) * ( exp(-(x - 1)^2) + exp(-(x + 1)^2) )
             + 0.2 * x^4 + 0.2 * (y - 1/3)^4

All functions are vectorised: ``x`` and ``y`` may be scalars or numpy arrays of
any broadcast-compatible shape, and the return value has the broadcast shape.
"""
from __future__ import annotations

import numpy as np

# Fixed centres used in the potential (kept as module constants so the gradient
# implementation cannot silently drift from the energy implementation).
_Y_A = 1.0 / 3.0   # centre of the first y Gaussian in the x-coupling term
_Y_B = 5.0 / 3.0   # centre of the second y Gaussian
_Y_Q = 1.0 / 3.0   # centre of the quartic confinement in y


def potential_xy(x, y):
    """Return the potential ``V(x, y)`` (vectorised)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    t1 = 3.0 * np.exp(-x**2) * (np.exp(-(y - _Y_A)**2) - np.exp(-(y - _Y_B)**2))
    t2 = -5.0 * np.exp(-y**2) * (np.exp(-(x - 1.0)**2) + np.exp(-(x + 1.0)**2))
    t3 = 0.2 * x**4 + 0.2 * (y - _Y_Q)**4
    return t1 + t2 + t3


def dVdx_xy(x, y):
    """Return ``dV/dx`` (the local mean force along ``xi(x, y) = x``)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    A = np.exp(-(y - _Y_A)**2) - np.exp(-(y - _Y_B)**2)
    d_t1 = 3.0 * (-2.0 * x) * np.exp(-x**2) * A
    B = np.exp(-y**2)
    d_t2 = -5.0 * B * ((-2.0 * (x - 1.0)) * np.exp(-(x - 1.0)**2)
                       + (-2.0 * (x + 1.0)) * np.exp(-(x + 1.0)**2))
    d_t3 = 0.8 * x**3
    return d_t1 + d_t2 + d_t3


def dVdy_xy(x, y):
    """Return ``dV/dy``."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    d_t1 = 3.0 * np.exp(-x**2) * (
        (-2.0 * (y - _Y_A)) * np.exp(-(y - _Y_A)**2)
        - (-2.0 * (y - _Y_B)) * np.exp(-(y - _Y_B)**2))
    C = np.exp(-(x - 1.0)**2) + np.exp(-(x + 1.0)**2)
    d_t2 = -5.0 * (-2.0 * y) * np.exp(-y**2) * C
    d_t3 = 0.8 * (y - _Y_Q)**3
    return d_t1 + d_t2 + d_t3


def grad_potential_xy(x, y):
    """Return the tuple ``(dV/dx, dV/dy)`` (vectorised)."""
    return dVdx_xy(x, y), dVdy_xy(x, y)


def make_grid(x_min, x_max, y_min, y_max, nx, ny):
    """Build evaluation grids for the potential.

    Returns
    -------
    x_grid : 1-D array of length ``nx``
    y_grid : 1-D array of length ``ny``
    XX, YY : 2-D arrays of shape ``(ny, nx)`` (``indexing="xy"``) so that
             ``XX[j, i] = x_grid[i]`` and ``YY[j, i] = y_grid[j]``.
    """
    x_grid = np.linspace(x_min, x_max, int(nx))
    y_grid = np.linspace(y_min, y_max, int(ny))
    XX, YY = np.meshgrid(x_grid, y_grid, indexing="xy")
    return x_grid, y_grid, XX, YY


def finite_difference_grad(x, y, eps=1e-6):
    """Central finite-difference gradient of :func:`potential_xy` (for tests)."""
    dVdx = (potential_xy(x + eps, y) - potential_xy(x - eps, y)) / (2.0 * eps)
    dVdy = (potential_xy(x, y + eps) - potential_xy(x, y - eps)) / (2.0 * eps)
    return dVdx, dVdy


# --------------------------------------------------------------------------- #
# PyTorch versions of the potential and its gradient
# --------------------------------------------------------------------------- #
# These mirror the numpy functions above term-for-term so that the GPU backend
# (:mod:`abffr.simulation_torch`) uses *exactly* the same energy/force as the
# reference CPU engine.  ``scripts/validate_torch_backend.py`` checks the torch
# gradient against finite differences of :func:`potential_xy_torch`.  The torch
# import is optional so that ``import abffr.potentials`` keeps working in a
# numpy-only environment (the CPU reference must never depend on torch).
try:  # pragma: no cover - exercised only when torch is installed
    import torch as _torch
except Exception:  # pragma: no cover
    _torch = None


def _require_torch():
    if _torch is None:
        raise ImportError(
            "PyTorch is required for the *_torch potential functions. "
            "Install a CPU build with e.g. "
            "`pip install torch --index-url https://download.pytorch.org/whl/cpu`."
        )
    return _torch


def potential_xy_torch(x, y):
    """Torch version of :func:`potential_xy` (matches the numpy implementation)."""
    torch = _require_torch()
    t1 = 3.0 * torch.exp(-x**2) * (
        torch.exp(-(y - _Y_A)**2) - torch.exp(-(y - _Y_B)**2))
    t2 = -5.0 * torch.exp(-y**2) * (
        torch.exp(-(x - 1.0)**2) + torch.exp(-(x + 1.0)**2))
    t3 = 0.2 * x**4 + 0.2 * (y - _Y_Q)**4
    return t1 + t2 + t3


def dVdx_xy_torch(x, y):
    """Torch version of :func:`dVdx_xy` (the local mean force ``f = dV/dx``)."""
    torch = _require_torch()
    A = torch.exp(-(y - _Y_A)**2) - torch.exp(-(y - _Y_B)**2)
    d_t1 = 3.0 * (-2.0 * x) * torch.exp(-x**2) * A
    B = torch.exp(-y**2)
    d_t2 = -5.0 * B * ((-2.0 * (x - 1.0)) * torch.exp(-(x - 1.0)**2)
                       + (-2.0 * (x + 1.0)) * torch.exp(-(x + 1.0)**2))
    d_t3 = 0.8 * x**3
    return d_t1 + d_t2 + d_t3


def dVdy_xy_torch(x, y):
    """Torch version of :func:`dVdy_xy`."""
    torch = _require_torch()
    d_t1 = 3.0 * torch.exp(-x**2) * (
        (-2.0 * (y - _Y_A)) * torch.exp(-(y - _Y_A)**2)
        - (-2.0 * (y - _Y_B)) * torch.exp(-(y - _Y_B)**2))
    C = torch.exp(-(x - 1.0)**2) + torch.exp(-(x + 1.0)**2)
    d_t2 = -5.0 * (-2.0 * y) * torch.exp(-y**2) * C
    d_t3 = 0.8 * (y - _Y_Q)**3
    return d_t1 + d_t2 + d_t3


def grad_potential_xy_torch(x, y):
    """Return ``(dV/dx, dV/dy)`` as torch tensors (vectorised)."""
    return dVdx_xy_torch(x, y), dVdy_xy_torch(x, y)
