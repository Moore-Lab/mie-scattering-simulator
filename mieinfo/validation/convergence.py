"""G-CONV — convergence gates (VALIDATION.md §5, PHYSICS.md §7).
Time convention e^{-iωt}.

Two reusable checks on the current plane-wave engine, across x in {3, 15, 30, 60}:

1. Series truncation (`check_nmax_convergence`): efficiencies (Q_ext, Q_sca, Q_back, g)
   are stable to <= 1e-4 relative under n_max -> n_max + 8 beyond the Wiscombe value.
   (Q_ext/Q_sca/g are stable to ~1e-13; Q_back to ~1e-12 — well inside the gate.)

2. Angular grid (`check_angular_grid_convergence`): smooth solid-angle integrals are
   stable to <= 1e-4 under doubling the Gauss-Legendre theta grid. The gate uses two
   smooth integrals (a hard-edged collection cone is NOT a smooth integral — its mask
   boundary quantizes onto the node set, so it is deliberately excluded, per VALIDATION.md
   §5 "Solid-angle integrals ... Gauss-Legendre in cos theta"):
   - the angle-integrated asymmetry parameter g (intensity-weighted <cos theta>) matches
     the coefficient-sum g to <= 1e-4 and is stable under doubling (PHYSICS.md §7 reports
     this passing to <= 1e-10 up to x = 236);
   - the total Fisher information F_q(4pi) for x- and z-displacement is stable under
     doubling.
   Larger x needs more theta nodes (forward lobe width ~ 1/x); the check scales the node
   count with x and reports the grid at which each quantity is converged.
"""
from __future__ import annotations

import numpy as np

from ..glmt.beam import PlaneWave
from ..glmt.scatter import PlaneWaveProvider
from ..info.modes import information_pattern
from ..mie import plane_wave as pw
from ..mie.special import nmax_wiscombe
from ..types import AngularGrid, Medium, Sphere

# Contract tolerances (VALIDATION.md §5).
TOL_NMAX = 1e-4
TOL_GRID = 1e-4

DEFAULT_X = (3.0, 15.0, 30.0, 60.0)
_MED_532 = Medium(n=1.0, wavelength_vacuum_m=532e-9)
_M_SILICA = 1.46 + 0j


# --------------------------------------------------------------------------- #
# 1. Series truncation
# --------------------------------------------------------------------------- #
def check_nmax_convergence(
    m: complex = _M_SILICA, xs: tuple[float, ...] = DEFAULT_X, delta: int = 8,
) -> tuple[bool, float]:
    """Efficiencies stable under n_max -> n_max + delta. Returns (passed, max_rel)."""
    names = ("Qext", "Qsca", "Qback", "g")
    passed = True
    max_rel = 0.0
    for x in xs:
        base = nmax_wiscombe(x)
        e0 = pw.efficiencies(m, x, base)
        e1 = pw.efficiencies(m, x, base + delta)
        for _name, a, b in zip(names, e0, e1):
            rel = abs(a - b) / max(abs(b), 1e-30)
            max_rel = max(max_rel, rel)
            if rel > TOL_NMAX:
                passed = False
    return passed, max_rel


# --------------------------------------------------------------------------- #
# 2. Angular grid — smooth solid-angle integrals
# --------------------------------------------------------------------------- #
def _ntheta_for(x: float) -> int:
    """theta node count resolving the forward lobe (~1/x) at size x (PHYSICS.md §7)."""
    return int(max(60, np.ceil(4.0 * x)))


def _g_from_grid(sphere: Sphere, ntheta: int, nphi: int) -> float:
    """Intensity-weighted <cos theta> from the far-field on a full-sphere grid."""
    provider = PlaneWaveProvider()
    beam = PlaneWave(_MED_532)
    grid = AngularGrid.full_sphere(ntheta, nphi)
    field = provider.field(grid, sphere, beam, np.zeros(3))
    intensity = field.intensity()
    w = grid.w_solid
    cos_theta = np.cos(grid.theta)[:, None]
    return float(np.sum(intensity * cos_theta * w) / np.sum(intensity * w))


def _fisher_total(sphere: Sphere, n_hat: np.ndarray, ntheta: int, nphi: int) -> float:
    provider = PlaneWaveProvider()
    beam = PlaneWave(_MED_532)
    grid = AngularGrid.full_sphere(ntheta, nphi)
    return information_pattern(provider, grid, sphere, beam, np.zeros(3), n_hat).f_total


def check_angular_grid_convergence(
    xs: tuple[float, ...] = DEFAULT_X, nphi: int = 8,
) -> tuple[bool, dict]:
    """Smooth solid-angle integrals stable under theta-grid doubling and (for g) matching
    the coefficient g. Returns (passed, per_x) where per_x[x] records the converged grid
    and the observed residuals.

    nphi=8 is ample: the phi-integrated smooth quantities used here are exact for the
    low-order cos(2 phi) azimuth structure of an x-pol dipole/Mie pattern under
    Gauss-Legendre-in-mu x uniform-phi quadrature; g uses a small nphi for speed.
    """
    passed = True
    per_x: dict[float, dict] = {}
    for x in xs:
        sphere = Sphere(radius_m=x / _MED_532.k, m=_M_SILICA)
        ntheta = _ntheta_for(x)

        g_coeff = pw.efficiencies(_M_SILICA, x)[3]
        g1 = _g_from_grid(sphere, ntheta, nphi)
        g2 = _g_from_grid(sphere, 2 * ntheta, nphi)
        g_vs_coeff = abs(g1 - g_coeff) / abs(g_coeff)
        g_double = abs(g1 - g2) / abs(g2)

        fx1 = _fisher_total(sphere, np.array([1.0, 0.0, 0.0]), ntheta, nphi)
        fx2 = _fisher_total(sphere, np.array([1.0, 0.0, 0.0]), 2 * ntheta, 2 * nphi)
        fz1 = _fisher_total(sphere, np.array([0.0, 0.0, 1.0]), ntheta, nphi)
        fz2 = _fisher_total(sphere, np.array([0.0, 0.0, 1.0]), 2 * ntheta, 2 * nphi)
        fx_double = abs(fx1 - fx2) / max(abs(fx2), 1e-30)
        fz_double = abs(fz1 - fz2) / max(abs(fz2), 1e-30)

        residuals = {
            "g_vs_coeff": g_vs_coeff,
            "g_double": g_double,
            "fisher_x_double": fx_double,
            "fisher_z_double": fz_double,
        }
        per_x[x] = {"ntheta": ntheta, "nphi": nphi, "residuals": residuals}
        if max(residuals.values()) > TOL_GRID:
            passed = False
    return passed, per_x
