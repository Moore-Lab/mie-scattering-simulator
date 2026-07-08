"""G-CONV tests (VALIDATION.md §5, PHYSICS.md §7) on the current plane-wave engine.
Time convention e^{-iωt}.

Covers mieinfo.validation.convergence:
  - efficiencies stable under n_max -> n_max + 8 (<= 1e-4 relative);
  - smooth solid-angle integrals stable under angular-grid doubling (<= 1e-4): the
    angle-integrated asymmetry g (also matching the coefficient g) and the total Fisher
    information F_q(4pi), across x in {3, 15, 30, 60}.

Hard-edged collection cones are deliberately NOT part of the grid-convergence gate: a
sharp NA mask quantizes onto the Gauss-Legendre node set, so its eta wobbles at ~1e-3
under doubling for a reason that is quadrature bookkeeping, not engine error. The gate is
on smooth integrals, per VALIDATION.md §5.
"""
from __future__ import annotations

import pytest

from mieinfo.validation import convergence

XS = convergence.DEFAULT_X


# --------------------------------------------------------------------------- #
# Series truncation
# --------------------------------------------------------------------------- #
def test_nmax_convergence_gate_passes():
    passed, max_rel = convergence.check_nmax_convergence()
    assert passed
    assert max_rel <= convergence.TOL_NMAX


@pytest.mark.parametrize("x", XS)
def test_nmax_convergence_per_x(x):
    passed, max_rel = convergence.check_nmax_convergence(xs=(x,))
    assert passed, f"n_max+8 residual {max_rel:.2e} at x={x}"
    assert max_rel <= convergence.TOL_NMAX


# --------------------------------------------------------------------------- #
# Angular grid — smooth solid-angle integrals
# --------------------------------------------------------------------------- #
def test_angular_grid_convergence_gate_passes():
    passed, per_x = convergence.check_angular_grid_convergence()
    assert passed, per_x
    for x, info in per_x.items():
        for name, resid in info["residuals"].items():
            assert resid <= convergence.TOL_GRID, (
                f"{name} residual {resid:.2e} at x={x} (grid ntheta={info['ntheta']})")


@pytest.mark.parametrize("x", XS)
def test_integrated_g_matches_coefficient_g(x):
    """The angle-integrated <cos theta> equals the coefficient-sum g to <= 1e-4 (in fact
    ~1e-13), confirming the Gauss-Legendre grid resolves the full angular structure."""
    _passed, per_x = convergence.check_angular_grid_convergence(xs=(x,))
    assert per_x[x]["residuals"]["g_vs_coeff"] <= convergence.TOL_GRID


@pytest.mark.parametrize("x", XS)
def test_fisher_total_stable_under_grid_doubling(x):
    """F_q(4pi) for x- and z-displacement is stable to <= 1e-4 under theta/phi doubling."""
    _passed, per_x = convergence.check_angular_grid_convergence(xs=(x,))
    resid = per_x[x]["residuals"]
    assert resid["fisher_x_double"] <= convergence.TOL_GRID
    assert resid["fisher_z_double"] <= convergence.TOL_GRID


@pytest.mark.slow
def test_angular_grid_convergence_large_x():
    """Extend the smooth-integral grid gate to the operative large sizes (x up to 236 at
    a = 20 um, PHYSICS.md §7). Nightly/full lane."""
    passed, per_x = convergence.check_angular_grid_convergence(xs=(120.0, 236.0), nphi=8)
    assert passed, per_x
