"""G-LIMIT tests (VALIDATION.md §3, PHYSICS.md §4.1, §6) on the current plane-wave
engine. Time convention e^{-iωt}.

Covers the three analytic-limit families from mieinfo.validation.limits:
  - energy / optical theorem (real m: Q_ext == Q_sca; absorbing m: Q_ext > Q_sca);
  - dipole-limit information-pattern structure (transverse forward-heavy vs axial,
    axial backward-weighted) -- the Tebbenjohanns-2019 dipole result;
  - reciprocity/symmetry: phi-parity of the x-pol intensity and info patterns, and the
    exact x/y information cross identity (with the documented non-symmetry of the totals).
"""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider
from mieinfo.info.modes import information_pattern
from mieinfo.mie import plane_wave as pw
from mieinfo.types import AngularGrid, Medium, Sphere
from mieinfo.validation import limits

MED = Medium(n=1.0, wavelength_vacuum_m=532e-9)


# --------------------------------------------------------------------------- #
# Energy / optical theorem
# --------------------------------------------------------------------------- #
def test_optical_theorem_gate_passes():
    passed, max_rel = limits.check_energy_optical_theorem()
    assert passed
    assert max_rel <= limits.TOL_OPTICAL_THEOREM


@pytest.mark.parametrize("m,x", [(1.46 + 0j, 3.0), (1.46 + 0j, 15.0),
                                 (1.46 + 0j, 30.0), (1.33 + 0j, 3.0)])
def test_real_m_qext_equals_qsca(m, x):
    q_ext, q_sca, _, _ = pw.efficiencies(m, x)
    assert abs(q_ext - q_sca) / q_sca <= 1e-10


@pytest.mark.parametrize("m,x", [(1.5 + 0.01j, 10.0), (1.46 + 0.1j, 15.0),
                                 (1.6 + 0.5j, 5.0)])
def test_absorbing_m_qext_exceeds_qsca(m, x):
    q_ext, q_sca, _, _ = pw.efficiencies(m, x)
    assert q_ext > q_sca


# --------------------------------------------------------------------------- #
# Dipole-limit information-pattern structure
# --------------------------------------------------------------------------- #
def test_dipole_information_structure_gate_passes():
    passed, metrics = limits.check_dipole_information_structure(x=0.05)
    assert passed, metrics


def test_dipole_transverse_more_forward_than_axial():
    """Transverse (x) displacement info is far more forward-weighted than axial (z)."""
    _, metrics = limits.check_dipole_information_structure(x=0.05)
    assert metrics["x_fwd_frac"] > metrics["z_fwd_frac"] + 0.3


def test_dipole_axial_backward_weighted():
    """Axial (z) info peaks in backscatter and is dominantly behind the sphere
    (Tebbenjohanns 2019: backscatter-optimal axial detection)."""
    _, metrics = limits.check_dipole_information_structure(x=0.05)
    assert metrics["z_peak_deg"] > 90.0
    assert metrics["z_fwd_frac"] < 0.15   # -> axial backward fraction > 0.85


def test_dipole_transverse_peak_off_axis():
    """Transverse info vanishes forward (sin^2 theta weight) and peaks off-axis, not in
    backscatter."""
    _, metrics = limits.check_dipole_information_structure(x=0.05)
    assert 20.0 < metrics["x_peak_deg"] < 90.0


def test_axial_forward_weight_vanishes_on_axis():
    """The axial info weight (1 - cos theta)^2 is exactly zero in the exact forward
    direction (theta = 0): axial motion is invisible to a perfectly forward detector."""
    grid = AngularGrid.full_sphere(120, 4)
    sphere = Sphere(radius_m=0.05 / MED.k, m=1.46 + 0j)
    beam = PlaneWave(MED)
    provider = PlaneWaveProvider()
    pz = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([0.0, 0.0, 1.0]))
    # The most-forward node's axial info density is far below the pattern peak.
    fwd_node = int(np.argmin(grid.theta))
    peak = np.max(pz.density)
    assert pz.density[fwd_node].max() < 1e-2 * peak


# --------------------------------------------------------------------------- #
# Reciprocity / symmetry
# --------------------------------------------------------------------------- #
def test_phi_parity_gate_passes():
    passed, max_relerr = limits.check_symmetry_phi_parity(x=6.0)
    assert passed
    assert max_relerr <= limits.TOL_PARITY


@pytest.mark.parametrize("x", [0.5, 6.0, 20.0])
def test_phi_parity_holds_across_sizes(x):
    passed, max_relerr = limits.check_symmetry_phi_parity(x=x)
    assert passed, f"phi-parity residual {max_relerr:.2e} at x={x}"


def test_xy_cross_identity_gate_passes():
    passed, max_relerr = limits.check_xy_cross_identity(x=6.0)
    assert passed
    assert max_relerr <= limits.TOL_CROSS_IDENTITY


def test_xy_cross_identity_is_the_exact_symmetry():
    """density_y * cos^2(phi) == density_x * sin^2(phi) to machine precision, computed
    directly here (independent of the library helper) so the documented symmetry is
    pinned by the test, not just asserted by the module."""
    grid = AngularGrid.full_sphere(50, 72)
    sphere = Sphere(radius_m=6.0 / MED.k, m=1.46 + 0j)
    beam = PlaneWave(MED)
    provider = PlaneWaveProvider()
    px = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([1.0, 0.0, 0.0]))
    py = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([0.0, 1.0, 0.0]))
    cos2 = np.cos(grid.phi)[None, :] ** 2
    sin2 = np.sin(grid.phi)[None, :] ** 2
    scale = max(np.max(px.density), np.max(py.density))
    resid = np.max(np.abs(py.density * cos2 - px.density * sin2)) / scale
    assert resid < 1e-12


def test_naive_phi_relabel_does_not_hold():
    """Documented non-symmetry: the plain phi -> phi - pi/2 relabel does NOT map the
    x-pattern onto the y-pattern (the shared intensity I = |S2|^2 cos^2 + |S1|^2 sin^2 is
    not pi/2-periodic when S1 != S2). Guards against anyone re-asserting the naive claim."""
    nphi = 72
    grid = AngularGrid.full_sphere(50, nphi)
    sphere = Sphere(radius_m=6.0 / MED.k, m=1.46 + 0j)
    beam = PlaneWave(MED)
    provider = PlaneWaveProvider()
    px = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([1.0, 0.0, 0.0]))
    py = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([0.0, 1.0, 0.0]))
    shift = nphi // 4   # exactly pi/2 in phi index
    rolled = np.roll(px.density, shift, axis=1)
    reldiff = np.max(np.abs(rolled - py.density)) / np.max(py.density)
    assert reldiff > 0.1   # emphatically not a match


def test_xy_total_fisher_not_equal():
    """Documented non-symmetry of the TOTALS: F_x != F_y for a Mie sphere (they would be
    equal only if |S1| == |S2|). In the dipole limit F_y -> 2 F_x."""
    grid = AngularGrid.full_sphere(200, 64)
    sphere = Sphere(radius_m=0.05 / MED.k, m=1.46 + 0j)
    beam = PlaneWave(MED)
    provider = PlaneWaveProvider()
    fx = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([1.0, 0.0, 0.0])).f_total
    fy = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([0.0, 1.0, 0.0])).f_total
    # dipole limit: F_y == 2 F_x
    assert abs(fy / fx - 2.0) < 1e-2
