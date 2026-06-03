"""PyTorch helpers for the GPU ABF--FR backend.

This module holds the device/dtype/seeding plumbing and the 1-D reaction-
coordinate grid primitives (binning, Gaussian smoothing, cumulative trapezoid,
linear interpolation) shared by :mod:`abffr.simulation_torch`.

Design note -- matching the CPU kernel estimator
------------------------------------------------
The CPU reference (:mod:`abffr.simulation`) estimates the conditional mean force
with a Nadaraya--Watson kernel regression evaluated *at every grid node*::

    Fprime_hat[j] = (sum_i K_h(x_j - X_i) dV/dx_i) / (sum_i K_h(x_j - X_i) + min_count)

with ``K_h`` a *density-normalised* Gaussian (``integral K_h = 1``).  The GPU
backend instead bins the particles onto the grid (counts ``C_j`` and force-sums
``S_j``) and convolves with a discretised, density-normalised Gaussian ``g``::

    Fprime_hat[j] = conv(S, g)[j] / (conv(C, g)[j] + min_count)

Because ``g`` is density-normalised (``sum_k g_k * dx = 1``), ``conv(C, g)[j]``
is a Riemann-sum approximation of ``sum_i K_h(x_j - X_i)`` -- the *same* units as
the CPU denominator, so ``min_count`` keeps its meaning -- and the estimator
converges to the CPU one as ``dx -> 0``.  This is the "binned-smooth" estimator;
:func:`abffr.simulation_torch` also offers an exact ``kernel_reference`` mode for
validation.
"""
from __future__ import annotations

from typing import Optional, Tuple

import torch
import torch.nn.functional as F


# --------------------------------------------------------------------------- #
# Device / dtype / RNG plumbing
# --------------------------------------------------------------------------- #
def resolve_device(requested: Optional[str] = None) -> torch.device:
    """Resolve a requested device string to an available :class:`torch.device`.

    ``"cuda"`` falls back to CPU with a printed warning when CUDA is
    unavailable (so a local CPU smoke run of a ``device: cuda`` config still
    works); ``"auto"`` picks CUDA when present.  Any explicit ``"cpu"`` is
    honoured as-is.
    """
    req = (requested or "auto").lower()
    has_cuda = torch.cuda.is_available()
    if req in ("auto", None):
        return torch.device("cuda" if has_cuda else "cpu")
    if req.startswith("cuda"):
        if has_cuda:
            return torch.device(req)
        print(f"[torch_utils] WARNING: device {req!r} requested but CUDA is "
              f"unavailable; falling back to CPU (run the CUDA benchmark "
              f"remotely).")
        return torch.device("cpu")
    return torch.device(req)


def resolve_dtype(name: Optional[str]) -> torch.dtype:
    """Map a dtype name (``"float32"``/``"float64"``) to a torch dtype."""
    return {"float32": torch.float32, "float": torch.float32,
            "float64": torch.float64, "double": torch.float64,
            None: torch.float32}.get(
                (name.lower() if isinstance(name, str) else name), torch.float32)


def make_generator(seed: int, device: torch.device) -> torch.Generator:
    """A seeded :class:`torch.Generator` on ``device`` (CUDA-safe)."""
    g = torch.Generator(device=device)
    g.manual_seed(int(seed) & ((1 << 63) - 1))
    return g


def stable_seed(*parts) -> int:
    """A deterministic non-negative 63-bit seed from arbitrary hashable parts.

    Used so a batch's RNG depends only on the runs it contains (and a base
    seed), making the whole backend reproducible across machines.
    """
    import hashlib
    h = hashlib.sha256("|".join(str(p) for p in parts).encode()).hexdigest()
    return int(h[:15], 16)  # 60 bits, comfortably < 2**63


# --------------------------------------------------------------------------- #
# Finite-value guards (fail loudly on NaN/inf, per the safety rules)
# --------------------------------------------------------------------------- #
class NonFiniteError(RuntimeError):
    """Raised when a tensor contains NaN or inf inside the integrator."""


def assert_finite(name: str, t: torch.Tensor) -> None:
    if not torch.isfinite(t).all():
        n_nan = int(torch.isnan(t).sum())
        n_inf = int(torch.isinf(t).sum())
        raise NonFiniteError(
            f"{name} contains non-finite values (nan={n_nan}, inf={n_inf}, "
            f"shape={tuple(t.shape)}).")


# --------------------------------------------------------------------------- #
# 1-D grid primitives
# --------------------------------------------------------------------------- #
def grid_spacing(x_grid: torch.Tensor) -> float:
    """Uniform spacing ``dx`` of a 1-D grid (assumed uniform, as built by
    :func:`abffr.reference.profile_grid`)."""
    return float((x_grid[-1] - x_grid[0]) / (x_grid.numel() - 1))


def gaussian_kernel1d(bandwidth: float, dx: float, device: torch.device,
                      dtype: torch.dtype, truncate: float = 4.0
                      ) -> Tuple[torch.Tensor, int]:
    """Density-normalised 1-D Gaussian kernel sampled on the grid.

    Returns ``(kernel, radius)`` where ``kernel`` has length ``2*radius + 1``
    and satisfies ``sum(kernel) * dx == 1`` (density normalisation), matching
    the CPU :func:`abffr.simulation.gaussian_kernel`.  ``bandwidth`` is in
    x-units (the config's ``abf.h`` or ``fr.eta``).
    """
    radius = max(1, int(round(truncate * bandwidth / dx)))
    offsets = torch.arange(-radius, radius + 1, device=device, dtype=dtype) * dx
    k = torch.exp(-0.5 * (offsets / bandwidth) ** 2)
    k = k / (k.sum() * dx)  # density normalisation: sum(k)*dx = 1
    return k, radius


def smooth_grid(values: torch.Tensor, kernel: torch.Tensor, radius: int,
                dx: float) -> torch.Tensor:
    """Convolve ``(B, G)`` grid values with a density-normalised ``kernel``.

    This is a *plain* discrete convolution (no extra ``dx`` factor): because the
    kernel samples the density-normalised Gaussian ``k[m] ~= K_h(m*dx)``,
    ``conv(C, k)[j] ~= sum_k C[k] K_h(x_j - x_k) ~= sum_i K_h(x_j - X_i)`` for a
    histogram-count grid ``C`` -- i.e. the *same* units as the CPU kernel
    denominator, so ``abf.min_count`` keeps its meaning and the binned numerator
    and denominator divide to the CPU estimator (the ``dx`` cancels in the
    ratio).  For a density (``p_hat``) the caller renormalises afterwards.

    Reflect padding at the domain edges mimics the CPU KDE's mirror-image
    boundary correction.  ``dx`` is accepted for signature symmetry but is not
    used here.
    """
    B, G = values.shape
    pad = min(radius, G - 1)  # torch reflect pad requires pad <= G-1
    xpad = F.pad(values.unsqueeze(1), (pad, pad), mode="reflect")
    w = kernel.to(values.dtype).view(1, 1, -1)
    if pad < radius:  # tiny grid: crop the kernel to the available pad
        crop = radius - pad
        w = w[..., crop:w.shape[-1] - crop]
    out = F.conv1d(xpad, w)
    return out.squeeze(1)


def cumulative_trapezoid(y: torch.Tensor, dx: float) -> torch.Tensor:
    """Cumulative trapezoid of ``(B, G)`` along the grid, with a leading 0.

    Matches ``scipy.integrate.cumulative_trapezoid(..., initial=0.0)`` used by
    the CPU engine to integrate the mean force into the bias potential.
    """
    trap = 0.5 * (y[:, 1:] + y[:, :-1]) * dx
    cum = torch.cumsum(trap, dim=1)
    zero = torch.zeros((y.shape[0], 1), device=y.device, dtype=y.dtype)
    return torch.cat([zero, cum], dim=1)


def trapezoid(y: torch.Tensor, dx: float) -> torch.Tensor:
    """Trapezoidal integral of ``(B, G)`` over the grid -> ``(B,)``."""
    return (y.sum(dim=1) - 0.5 * (y[:, 0] + y[:, -1])) * dx


def center_at_index(F_grid: torch.Tensor, idx0: int) -> torch.Tensor:
    """Subtract ``F[:, idx0]`` so the bias potential is pinned to 0 near x=0."""
    return F_grid - F_grid[:, idx0:idx0 + 1]


def nearest_index(X: torch.Tensor, x0: float, dx: float, G: int) -> torch.Tensor:
    """Nearest grid-node index for each particle position -> long ``(B, N)``."""
    idx = torch.round((X - x0) / dx).long()
    return idx.clamp_(0, G - 1)


def scatter_grid(idx: torch.Tensor, G: int, values: Optional[torch.Tensor] = None
                 ) -> torch.Tensor:
    """Scatter-add per-row particle contributions onto a ``(B, G)`` grid.

    ``values is None`` accumulates counts (histogram); otherwise accumulates the
    summed ``values`` (e.g. ``dV/dx``) per bin.  Uses ``scatter_add_`` as
    requested by the study spec.
    """
    B, N = idx.shape
    out = torch.zeros((B, G), device=idx.device,
                      dtype=(values.dtype if values is not None else torch.float32))
    src = torch.ones_like(idx, dtype=out.dtype) if values is None else values.to(out.dtype)
    out.scatter_add_(1, idx, src)
    return out


def interp1d(values: torch.Tensor, X: torch.Tensor, x0: float, dx: float
             ) -> torch.Tensor:
    """Linear interpolation of ``(B, G)`` grid values at particle positions.

    ``values[b]`` is sampled at ``X[b]`` (both row ``b``).  Mirrors
    ``numpy.interp`` with edge clamping; preferred over nearest-bin gather for
    the ABF bias force at particle positions.
    """
    G = values.shape[1]
    pos = (X - x0) / dx
    i0 = torch.floor(pos).long().clamp_(0, G - 2)
    frac = (pos - i0.to(values.dtype)).clamp_(0.0, 1.0)
    v0 = torch.gather(values, 1, i0)
    v1 = torch.gather(values, 1, i0 + 1)
    return v0 + frac * (v1 - v0)


def reflect_into(q: torch.Tensor, lo: float, hi: float) -> torch.Tensor:
    """Reflect coordinates into ``[lo, hi]`` (torch version of ``reflect_1d``)."""
    span = hi - lo
    z = torch.remainder(q - lo, 2.0 * span)
    z = torch.where(z > span, 2.0 * span - z, z)
    return z + lo


def normalize_density(p: torch.Tensor, dx: float, eps: float = 1e-12
                      ) -> torch.Tensor:
    """Normalise ``(B, G)`` densities to integrate to 1 over the grid."""
    Z = trapezoid(p, dx).clamp_min(eps).unsqueeze(1)
    return p / Z
