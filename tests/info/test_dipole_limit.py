"""M2 dipole-limit gate (VALIDATION.md §3, PHYSICS.md §4.1): at x≪1 the information
pattern reproduces the Tebbenjohanns-2019 dipole structure — transverse-motion info is
forward-heavy and peaks off-axis, axial-motion info is backward-weighted, and the
pattern is NOT the intensity pattern. Time convention e^{-iωt}."""
from __future__ import annotations

import numpy as np

from mieinfo.detect.optics import CollectionGeometry, collection_efficiency
from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider
from mieinfo.info.modes import information_pattern
from mieinfo.types import AngularGrid, Medium, Sphere


def _forward_fraction(pattern):
    fwd = pattern.grid.theta < np.pi / 2
    W = pattern.grid.w_solid
    return float(np.sum(pattern.density[fwd] * W[fwd]) / np.sum(pattern.density * W))


def _peak_theta(pattern):
    f_pol = np.sum(pattern.density * pattern.grid.w_solid, axis=1)
    return float(pattern.grid.theta[np.argmax(f_pol)])


def _dipole_patterns():
    med = Medium(n=1.0, wavelength_vacuum_m=532e-9)
    a = 0.05 / med.k                       # x = k a = 0.05 (deep Rayleigh)
    sph = Sphere(radius_m=a, m=1.46 + 0j)
    beam = PlaneWave(med)
    grid = AngularGrid.full_sphere(200, 96)
    prov = PlaneWaveProvider()
    px = information_pattern(prov, grid, sph, beam, np.zeros(3), np.array([1.0, 0, 0]))
    pz = information_pattern(prov, grid, sph, beam, np.zeros(3), np.array([0, 0, 1.0]))
    return prov, sph, beam, grid, px, pz


def test_axial_info_is_backward_weighted():
    _, _, _, _, px, pz = _dipole_patterns()
    assert _forward_fraction(pz) < 0.3                 # axial: backward-heavy
    assert _peak_theta(pz) > np.pi / 2                 # peaks in the backward hemisphere


def test_transverse_more_forward_than_axial():
    _, _, _, _, px, pz = _dipole_patterns()
    assert _forward_fraction(px) > _forward_fraction(pz)
    assert _forward_fraction(px) > 0.4


def test_transverse_peaks_off_axis():
    _, _, _, _, px, _ = _dipole_patterns()
    theta_pk = _peak_theta(px)
    assert 0.0 < theta_pk < np.pi / 2                  # forward but not on-axis (weight ∝ sin²θ)


def test_information_is_not_intensity():
    prov, sph, beam, grid, _, pz = _dipole_patterns()
    intensity = prov.field(grid, sph, beam, np.zeros(3)).intensity()
    i_peak = grid.theta[np.argmax(np.sum(intensity * grid.w_solid, axis=1))]
    assert i_peak < np.pi / 2 and _peak_theta(pz) > np.pi / 2   # intensity fwd, info_z bwd


def test_full_sphere_efficiency_is_one():
    _, _, _, _, px, pz = _dipole_patterns()
    full = CollectionGeometry(direction="both", NA=1.0)
    assert abs(collection_efficiency(px, full) - 1.0) < 1e-12
    assert abs(collection_efficiency(pz, full) - 1.0) < 1e-12
