"""PlaneWaveProvider + phase-gradient field_derivative (INTERFACES.md §4).
Includes the G-DERIV check (analytic vs finite-difference) for the plane-wave provider,
which must agree to machine precision (VALIDATION.md §4). Time convention e^{-iωt}."""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider, field_derivative
from mieinfo.mie import plane_wave as pw
from mieinfo.types import AngularGrid, Medium, Sphere


def _setup():
    med = Medium(n=1.0, wavelength_vacuum_m=532e-9)
    sph = Sphere(radius_m=1e-6, m=1.46 + 0j)   # x ~ 11.8
    beam = PlaneWave(med)
    grid = AngularGrid.full_sphere(200, 64)
    return med, sph, beam, grid


def test_field_intensity_matches_amplitudes():
    med, sph, beam, grid = _setup()
    vf = PlaneWaveProvider().field(grid, sph, beam, np.zeros(3))
    x = pw.size_parameter(sph, med)
    S1, S2 = pw.scattering_amplitudes(sph.m, x, grid.theta)
    I_expect = ((np.abs(S2)[:, None] * np.cos(grid.phi)[None, :]) ** 2
                + (np.abs(S1)[:, None] * np.sin(grid.phi)[None, :]) ** 2)
    assert np.allclose(vf.intensity(), I_expect, rtol=1e-12, atol=0)


def test_qsca_matches_efficiencies():
    med, sph, beam, _ = _setup()
    x = pw.size_parameter(sph, med)
    assert abs(PlaneWaveProvider().q_sca(sph, beam) - pw.efficiencies(sph.m, x)[1]) <= 1e-12


def test_displacement_is_pure_phase():
    _, sph, beam, grid = _setup()
    prov = PlaneWaveProvider()
    f0 = prov.field(grid, sph, beam, np.zeros(3))
    fr = prov.field(grid, sph, beam, np.array([0.1e-6, -0.03e-6, 0.05e-6]))
    assert np.allclose(fr.intensity(), f0.intensity(), rtol=1e-10)  # |E| unchanged


@pytest.mark.parametrize("r_s", [
    np.zeros(3),
    np.array([0.1e-6, 0.0, 0.0]),
    np.array([0.0, 0.0, 0.12e-6]),
    np.array([0.05e-6, -0.07e-6, 0.09e-6]),
])
def test_analytic_vs_finite_difference(r_s):
    _, sph, beam, grid = _setup()
    prov = PlaneWaveProvider()
    da = field_derivative(prov, grid, sph, beam, r_s, method="analytic")
    dn = field_derivative(prov, grid, sph, beam, r_s, method="finite_difference")
    for A, N in ((da.dE_theta, dn.dE_theta), (da.dE_phi, dn.dE_phi)):
        rel = np.max(np.abs(A - N)) / (np.max(np.abs(N)) + 1e-30)
        assert rel <= 1e-6, f"analytic vs FD rel err {rel:.2e}"


def test_nondefault_polarization_raises():
    """PlaneWaveProvider must fail loud on a non-+x beam-frame polarization, not
    silently return the x-pol field (lab pol is the channel's job, PHYSICS §4.5)."""
    _, sph, _, grid = _setup()
    ybeam = PlaneWave(Medium(n=1.0, wavelength_vacuum_m=532e-9), polarization=(0.0, 1.0))
    with pytest.raises(NotImplementedError):
        PlaneWaveProvider().field(grid, sph, ybeam, np.zeros(3))


def test_info_density_reweights_intensity():
    """|dE/dr_z|² = k²(1−cosθ)²·|E|²: axial info vanishes forward, peaks backward."""
    _, sph, beam, grid = _setup()
    prov = PlaneWaveProvider()
    f = prov.field(grid, sph, beam, np.zeros(3))
    d = field_derivative(prov, grid, sph, beam, np.zeros(3), method="analytic")
    dens_z = np.abs(d.dE_theta[2]) ** 2 + np.abs(d.dE_phi[2]) ** 2
    k = beam.medium.k
    expect = (k * (1.0 - np.cos(grid.theta))[:, None]) ** 2 * f.intensity()
    assert np.allclose(dens_z, expect, rtol=1e-12, atol=0)
    # forward (small θ) is strongly suppressed vs backward
    fwd = grid.theta < np.pi / 6
    bwd = grid.theta > 5 * np.pi / 6
    assert dens_z[fwd].max() < dens_z[bwd].max()
