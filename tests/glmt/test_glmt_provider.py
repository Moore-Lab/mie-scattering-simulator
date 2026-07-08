"""Tests for glmt.scatter.GLMTProvider — the GLMT FieldProvider (INTERFACES.md §4).

Gates covered:
- G-GOLD: a wide-waist Gaussian (and a plane-wave beam) through GLMT reproduces
  mie.plane_wave S1/S2 to <= 1e-6 (PHYSICS.md §2.1, VALIDATION.md §3).
- Dipole limit: at x = 0.05 the scattered field is the Rayleigh dipole field.
- Displacement: GLMTProvider.field(r_s) for a plane-wave beam equals the exact
  analytic phase shift (PHYSICS.md §3.1) to machine precision, and
  field_derivative(method='finite_difference') differentiates it.
Time convention e^{-iωt}.
"""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.glmt.beam import GaussianParaxial, PlaneWave, RichardsWolfFocus
from mieinfo.glmt.scatter import (GLMTProvider, PlaneWaveProvider,
                                  field_derivative)
from mieinfo.mie import plane_wave as pw
from mieinfo.types import AngularGrid, Medium, Sphere

MED = Medium(n=1.0, wavelength_vacuum_m=532e-9)
K = MED.k
LAM = 2.0 * np.pi / K


def _s1s2_field(sphere, grid):
    x = pw.size_parameter(sphere, MED)
    S1, S2 = pw.scattering_amplitudes(sphere.m, x, grid.theta)
    cph = np.cos(grid.phi)[None, :]
    sph = np.sin(grid.phi)[None, :]
    return S2[:, None] * cph, -S1[:, None] * sph, max(np.max(np.abs(S1)),
                                                       np.max(np.abs(S2)))


# ---------------------------------------------------------------------------
# G-GOLD: plane-wave limit reproduces S1/S2
# ---------------------------------------------------------------------------

def test_plane_wave_beam_through_glmt_reproduces_s1s2():
    sph = Sphere(radius_m=1e-6, m=1.4607 + 0j)   # x ~ 11.8
    grid = AngularGrid.full_sphere(120, 16)
    vf = GLMTProvider().field(grid, sph, PlaneWave(MED), np.zeros(3))
    et_ref, ep_ref, scale = _s1s2_field(sph, grid)
    assert np.max(np.abs(vf.E_theta - et_ref)) / scale <= 1e-6
    assert np.max(np.abs(vf.E_phi - ep_ref)) / scale <= 1e-6


def test_wide_waist_gaussian_reproduces_s1s2_ggold():
    """G-GOLD: wide-waist Gaussian GLMT field matches plane-wave S1/S2 (<= 1e-6)."""
    sph = Sphere(radius_m=1e-6, m=1.4607 + 0j)
    grid = AngularGrid.full_sphere(120, 16)
    g = GaussianParaxial(MED, waist_m=2000 * LAM)
    vf = GLMTProvider().field(grid, sph, g, np.zeros(3))
    et_ref, ep_ref, scale = _s1s2_field(sph, grid)
    assert np.max(np.abs(vf.E_theta - et_ref)) / scale <= 1e-6
    assert np.max(np.abs(vf.E_phi - ep_ref)) / scale <= 1e-6


@pytest.mark.parametrize("radius_um", [0.5, 2.0])
def test_glmt_qsca_matches_plane_wave(radius_um):
    sph = Sphere(radius_m=radius_um * 1e-6, m=1.4607 + 0j)
    x = pw.size_parameter(sph, MED)
    q_ref = pw.efficiencies(sph.m, x)[1]
    assert abs(GLMTProvider().q_sca(sph, PlaneWave(MED)) - q_ref) <= 1e-12


# ---------------------------------------------------------------------------
# Dipole (Rayleigh) limit at x = 0.05
# ---------------------------------------------------------------------------

def test_dipole_limit_field_shape():
    """At x = 0.05 the field is the Rayleigh dipole: E_theta ∝ cos(theta) cos(phi),
    E_phi ∝ -sin(phi) (PHYSICS.md §6, G-LIMIT dipole)."""
    a = 0.05 / K
    sph = Sphere(radius_m=a, m=1.5 + 0j)
    grid = AngularGrid.full_sphere(60, 24)
    vf = GLMTProvider().field(grid, sph, PlaneWave(MED), np.zeros(3))

    # phi = 0 column: E_theta ∝ cos(theta); normalize and compare shape. At x = 0.05
    # the leading dipole is contaminated by an O(x^2) quadrupole (a_2 term), so the
    # residual from a pure cos(theta) is ~1e-3, not zero — assert the dipole shape to
    # that physically honest level.
    i0 = int(np.argmin(np.abs(grid.phi - 0.0)))
    col = vf.E_theta[:, i0]
    shape = col / col[0]
    cos_shape = np.cos(grid.theta) / np.cos(grid.theta[0])
    assert np.max(np.abs(shape - cos_shape)) < 3e-3

    # E_phi ∝ -sin(phi) at fixed theta (dipole: S1 ~ const in theta)
    it = int(np.argmin(np.abs(grid.theta - np.pi / 3)))
    row = vf.E_phi[it, :]
    ref = -np.sin(grid.phi)
    # proportional up to a complex constant
    c = row[np.argmax(np.abs(ref))] / ref[np.argmax(np.abs(ref))]
    assert np.max(np.abs(row - c * ref)) / (np.max(np.abs(row)) + 1e-30) < 1e-3


# ---------------------------------------------------------------------------
# Displacement path (PHYSICS.md §3.1 / §3.2)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("r_s", [
    np.array([0.1e-6, -0.05e-6, 0.08e-6]),
    np.array([0.0, 0.0, 0.12e-6]),
    np.array([0.2e-6, 0.0, 0.0]),
])
def test_plane_wave_displacement_equals_exact_phase_shift(r_s):
    """GLMT(plane wave) displaced field == PlaneWaveProvider analytic phase shift
    (PHYSICS.md §3.1) to machine precision — validates BSC(center=r_s) + translation
    phase exp(-i k s.r_s)."""
    sph = Sphere(radius_m=1e-6, m=1.4607 + 0j)
    grid = AngularGrid.full_sphere(120, 24)
    beam = PlaneWave(MED)
    vf_g = GLMTProvider().field(grid, sph, beam, r_s)
    vf_p = PlaneWaveProvider().field(grid, sph, beam, r_s)
    scale = np.max(np.abs(vf_p.E_theta)) + np.max(np.abs(vf_p.E_phi))
    assert np.max(np.abs(vf_g.E_theta - vf_p.E_theta)) / scale <= 1e-10
    assert np.max(np.abs(vf_g.E_phi - vf_p.E_phi)) / scale <= 1e-10


def test_displacement_preserves_intensity_for_plane_wave():
    """A plane-wave displacement is a pure phase => |E| unchanged."""
    sph = Sphere(radius_m=1e-6, m=1.4607 + 0j)
    grid = AngularGrid.full_sphere(80, 16)
    prov = GLMTProvider()
    f0 = prov.field(grid, sph, PlaneWave(MED), np.zeros(3))
    fr = prov.field(grid, sph, PlaneWave(MED), np.array([0.1e-6, -0.03e-6, 0.05e-6]))
    assert np.allclose(fr.intensity(), f0.intensity(), rtol=1e-9)


def test_finite_difference_through_glmt_matches_plane_wave_analytic():
    """field_derivative(method='finite_difference') on GLMT(plane wave) matches the
    exact analytic phase-gradient derivative (G-DERIV, VALIDATION.md §4)."""
    sph = Sphere(radius_m=1e-6, m=1.4607 + 0j)
    grid = AngularGrid.full_sphere(120, 24)
    beam = PlaneWave(MED)
    r_s = np.array([0.05e-6, -0.03e-6, 0.06e-6])
    dfd = field_derivative(GLMTProvider(), grid, sph, beam, r_s,
                           method="finite_difference")
    dan = field_derivative(PlaneWaveProvider(), grid, sph, beam, r_s,
                           method="analytic")
    for A, N in ((dan.dE_theta, dfd.dE_theta), (dan.dE_phi, dfd.dE_phi)):
        rel = np.max(np.abs(A - N)) / (np.max(np.abs(N)) + 1e-30)
        assert rel <= 1e-6, f"FD(GLMT) vs analytic phase-grad rel {rel:.2e}"


def test_focused_beam_displacement_runs_and_changes_field():
    """A focused Gaussian gives a displacement-dependent far field (the info source)."""
    sph = Sphere(radius_m=1e-6, m=1.4607 + 0j)
    grid = AngularGrid.full_sphere(100, 24)
    prov = GLMTProvider(bsc_method="quadrature")
    g = GaussianParaxial(MED, waist_m=3 * LAM)
    f0 = prov.field(grid, sph, g, np.zeros(3))
    fr = prov.field(grid, sph, g, np.array([0.15e-6, 0.0, 0.0]))
    assert f0.E_theta.shape == (grid.theta.shape[0], grid.phi.shape[0])
    # transverse displacement changes the pattern (not a pure global phase)
    assert not np.allclose(fr.intensity(), f0.intensity(), rtol=1e-3)


def test_field_rejects_bad_rs_shape():
    sph = Sphere(radius_m=1e-6, m=1.4607 + 0j)
    grid = AngularGrid.full_sphere(40, 8)
    with pytest.raises(ValueError):
        GLMTProvider().field(grid, sph, PlaneWave(MED), np.zeros(2))


def test_bsc_cache_reused_across_calls():
    """BSCs depend only on (beam, center, n_max) and are cached (PHYSICS.md §7)."""
    sph = Sphere(radius_m=1e-6, m=1.4607 + 0j)
    grid = AngularGrid.full_sphere(60, 12)
    prov = GLMTProvider(bsc_method="quadrature")
    g = GaussianParaxial(MED, waist_m=3 * LAM)
    prov.field(grid, sph, g, np.zeros(3))
    n_before = len(prov._cache)
    prov.field(grid, sph, g, np.zeros(3))   # same center -> cache hit
    assert len(prov._cache) == n_before


def test_richards_wolf_through_glmt_runs():
    """A high-NA RW focus produces a finite GLMT far field (angular-spectrum BSCs)."""
    sph = Sphere(radius_m=1e-6, m=1.4607 + 0j)
    grid = AngularGrid.full_sphere(100, 24)
    rw = RichardsWolfFocus(MED, NA=0.6, n_quad=120)
    vf = GLMTProvider().field(grid, sph, rw, np.zeros(3))
    assert np.isfinite(vf.E_theta).all() and np.isfinite(vf.E_phi).all()
    assert np.max(np.abs(vf.E_theta)) > 0
