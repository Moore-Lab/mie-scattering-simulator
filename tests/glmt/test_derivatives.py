"""Tests for glmt.derivatives — analytic displacement derivatives dE_s/dr_j, and the
glmt.translation axial translation-coefficient kernel. Time convention e^{-iωt}.

Central gate: G-DERIV — analytic vs finite-difference dE/dr_j <= 1e-6 over
|r_s| <= 0.3 lambda (VALIDATION.md §4). For a plane wave through GLMT the analytic
derivative is the exact phase gradient; for a focused Gaussian it is the product-rule
form (glmt.derivatives), validated against the finite-difference reference here.
"""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.glmt.beam import GaussianParaxial, PlaneWave
from mieinfo.glmt.scatter import GLMTProvider, field_derivative
from mieinfo.types import AngularGrid, Medium, Sphere

MED = Medium(n=1.0, wavelength_vacuum_m=532e-9)
K = MED.k
LAM = 2.0 * np.pi / K


def _grid():
    return AngularGrid.full_sphere(120, 24)


def _sphere():
    return Sphere(radius_m=1e-6, m=1.4607 + 0j)   # x ~ 11.8


# ---------------------------------------------------------------------------
# G-DERIV: plane wave through GLMT — exact analytic phase gradient
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("r_s", [
    np.zeros(3),
    np.array([0.1e-6, 0.0, 0.0]),
    np.array([0.0, 0.0, 0.12e-6]),
    np.array([0.05e-6, -0.07e-6, 0.09e-6]),
])
def test_glmt_plane_wave_analytic_vs_finite_difference(r_s):
    """analytic == finite_difference for GLMT(plane wave) to <= 1e-6 (G-DERIV)."""
    grid, sph, beam = _grid(), _sphere(), PlaneWave(MED)
    prov = GLMTProvider()
    da = field_derivative(prov, grid, sph, beam, r_s, method="analytic")
    dn = field_derivative(prov, grid, sph, beam, r_s, method="finite_difference")
    for A, N in ((da.dE_theta, dn.dE_theta), (da.dE_phi, dn.dE_phi)):
        rel = np.max(np.abs(A - N)) / (np.max(np.abs(N)) + 1e-30)
        assert rel <= 1e-6, f"analytic vs FD rel {rel:.2e}"


def test_glmt_plane_wave_analytic_matches_phase_gradient_formula():
    """The GLMT analytic derivative for a plane wave equals i k (k_hat - s_hat)_j E
    (PHYSICS.md §3.1)."""
    grid, sph, beam = _grid(), _sphere(), PlaneWave(MED)
    prov = GLMTProvider()
    f = prov.field(grid, sph, beam, np.zeros(3))
    d = field_derivative(prov, grid, sph, beam, np.zeros(3), method="analytic")
    TH = grid.theta[:, None]; PH = grid.phi[None, :]
    sinT = np.sin(TH)
    fac_z = 1j * K * (1.0 - np.cos(TH) * np.ones_like(PH))
    assert np.allclose(d.dE_theta[2], fac_z * f.E_theta, rtol=1e-10, atol=0)
    fac_x = 1j * K * (-sinT * np.cos(PH))
    assert np.allclose(d.dE_theta[0], fac_x * f.E_theta, rtol=1e-10, atol=0)


# ---------------------------------------------------------------------------
# G-DERIV: focused Gaussian — product-rule analytic vs finite-difference
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("r_s", [
    np.array([0.0, 0.0, 0.0]),
    np.array([0.08e-6, -0.05e-6, 0.06e-6]),
    np.array([0.0, 0.0, 0.15e-6]),      # |r_s| ~ 0.28 lambda
])
def test_glmt_gaussian_analytic_vs_finite_difference(r_s):
    """G-DERIV for a focused Gaussian: product-rule analytic derivative matches the
    finite-difference reference to <= 1e-6 (VALIDATION.md §4)."""
    grid, sph = _grid(), _sphere()
    beam = GaussianParaxial(MED, waist_m=3 * LAM)
    prov = GLMTProvider(bsc_method="quadrature")
    da = field_derivative(prov, grid, sph, beam, r_s, method="analytic")
    dn = field_derivative(prov, grid, sph, beam, r_s, method="finite_difference")
    for A, N in ((da.dE_theta, dn.dE_theta), (da.dE_phi, dn.dE_phi)):
        rel = np.max(np.abs(A - N)) / (np.max(np.abs(N)) + 1e-30)
        assert rel <= 1e-6, f"analytic vs FD rel {rel:.2e}"


def test_axial_info_backward_weighted_for_focused_beam():
    """dF_z/dOmega ∝ |dE/dr_z|^2 is backward-weighted (PHYSICS.md §4.1): the axial-
    displacement information density is larger in backscatter than forward, even for a
    focused beam."""
    grid, sph = _grid(), _sphere()
    beam = GaussianParaxial(MED, waist_m=3 * LAM)
    prov = GLMTProvider(bsc_method="quadrature")
    d = field_derivative(prov, grid, sph, beam, np.zeros(3), method="analytic")
    dens_z = np.abs(d.dE_theta[2]) ** 2 + np.abs(d.dE_phi[2]) ** 2
    fwd = grid.theta < np.pi / 6
    bwd = grid.theta > 5 * np.pi / 6
    assert dens_z[fwd].max() < dens_z[bwd].max()


def test_derivative_shape_and_axes():
    grid, sph, beam = _grid(), _sphere(), PlaneWave(MED)
    d = field_derivative(GLMTProvider(), grid, sph, beam, np.zeros(3),
                         method="analytic")
    assert d.dE_theta.shape == (3,) + grid.w_solid.shape
    assert d.dE_phi.shape == (3,) + grid.w_solid.shape


# ---------------------------------------------------------------------------
# Scalar axial translation-addition kernel — PENDING (raises NotImplementedError). Two
# projection drafts failed a genuine radius-independent, all-degree reconstruction, so
# glmt.translation raises rather than shipping wrong numbers (tests/glmt/test_translation.py).
# The validated GLMT displacement derivative (G-DERIV, above) uses the finite-difference /
# product-rule BSC gradient, never a translation matrix.
# ---------------------------------------------------------------------------


def test_glmt_derivative_does_not_regress_with_translation_kernel():
    """Wiring note (honesty): the scalar translation kernel does NOT feed the vector
    displacement derivative. field_derivative(method='analytic') for a focused Gaussian
    still matches finite_difference to <= 1e-6 (G-DERIV) via the product-rule BSC
    gradient — the kernel's existence must not change that validated path."""
    grid, sph = _grid(), _sphere()
    beam = GaussianParaxial(MED, waist_m=3 * LAM)
    prov = GLMTProvider(bsc_method="quadrature")
    r_s = np.array([0.0, 0.0, 0.12e-6])
    da = field_derivative(prov, grid, sph, beam, r_s, method="analytic")
    dn = field_derivative(prov, grid, sph, beam, r_s, method="finite_difference")
    for A, N in ((da.dE_theta, dn.dE_theta), (da.dE_phi, dn.dE_phi)):
        rel = np.max(np.abs(A - N)) / (np.max(np.abs(N)) + 1e-30)
        assert rel <= 1e-6
