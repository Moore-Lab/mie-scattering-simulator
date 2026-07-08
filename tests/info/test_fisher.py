"""Fisher-information density/total and arbitrary-direction combination
(INTERFACES.md §5, PHYSICS.md §4.4). Time convention e^{-iωt}."""
from __future__ import annotations

import numpy as np

from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider, field_derivative
from mieinfo.info.fisher import fisher_total, info_density
from mieinfo.info.modes import combine_direction, information_pattern
from mieinfo.types import AngularGrid, Medium, Sphere


def _setup():
    med = Medium(n=1.0, wavelength_vacuum_m=532e-9)
    sph = Sphere(radius_m=1e-6, m=1.46 + 0j)
    return med, sph, PlaneWave(med), AngularGrid.full_sphere(120, 48)


def test_info_density_is_modulus_squared():
    _, _, _, grid = _setup()
    dEt = np.array([[1 + 1j, 2 + 0j]])
    dEp = np.array([[0 + 1j, 1 + 1j]])
    d = info_density((dEt, dEp), grid)
    assert np.allclose(d, np.abs(dEt) ** 2 + np.abs(dEp) ** 2)


def test_combine_direction_cartesian_axis_selects_component():
    _, sph, beam, grid = _setup()
    deriv = field_derivative(PlaneWaveProvider(), grid, sph, beam, np.zeros(3), method="analytic")
    dEt_x, dEp_x = combine_direction(deriv, np.array([1.0, 0.0, 0.0]))
    assert np.allclose(dEt_x, deriv.dE_theta[0]) and np.allclose(dEp_x, deriv.dE_phi[0])


def test_combine_direction_is_linear():
    _, sph, beam, grid = _setup()
    deriv = field_derivative(PlaneWaveProvider(), grid, sph, beam, np.zeros(3), method="analytic")
    n = np.array([1.0, 1.0, 0.0])
    dEt, dEp = combine_direction(deriv, n)
    expect_t = (deriv.dE_theta[0] + deriv.dE_theta[1]) / np.sqrt(2)
    expect_p = (deriv.dE_phi[0] + deriv.dE_phi[1]) / np.sqrt(2)
    assert np.allclose(dEt, expect_t) and np.allclose(dEp, expect_p)


def test_information_pattern_f_total_matches_integral():
    _, sph, beam, grid = _setup()
    pat = information_pattern(PlaneWaveProvider(), grid, sph, beam, np.zeros(3), np.array([1.0, 0, 0]))
    assert pat.f_total == fisher_total(pat.density, grid)
    assert pat.f_total > 0.0
    assert np.allclose(np.linalg.norm(pat.n_hat), 1.0)
