"""Tests for circular / elliptical polarization of RichardsWolfFocus (glmt.beam).

The Debye-Wolf focal field is LINEAR in the input transverse Jones vector, so an
arbitrary polarization (px, py) gives the exact superposition
    E(r) = px * E^{x-input}(r) + py * E^{y-input}(r),
with E^{y-input} the x-input field rotated 90 deg about +z. These tests pin that
relationship (previously the class silently returned the pure x-field for a circular
Jones vector — see docs/dissertation_comparison.md). Time convention e^{-iωt}.
"""
from __future__ import annotations

import numpy as np

from mieinfo.glmt.beam import RichardsWolfFocus
from mieinfo.types import Medium

MED = Medium(n=1.0, wavelength_vacuum_m=1064e-9)


def _pts():
    return np.array([
        [0.0, 0.0, 0.0],
        [0.2e-6, 0.1e-6, 0.3e-6],
        [0.5e-6, -0.3e-6, -0.2e-6],
        [-0.4e-6, 0.25e-6, 0.15e-6],
    ])


def test_xpol_unchanged_regression():
    """+x polarization is the unrotated Debye-Wolf field (bit-for-bit unchanged)."""
    rw = RichardsWolfFocus(MED, NA=0.63, polarization=(1.0, 0.0))
    pts = _pts()
    E = rw.focal_field(pts)
    E_ref = rw._focal_xyz(pts)  # the internal +x field
    assert np.allclose(E, E_ref, atol=0, rtol=0)


def test_circular_focal_field_relationship():
    """Right-circular (1, i)/sqrt2: Ey = i Ex at the focus, |E(0)| = 1."""
    rw = RichardsWolfFocus(MED, NA=0.63, polarization=(1.0, 1j))
    E0 = rw.focal_field(np.zeros((1, 3)))[0]
    # at the geometric focus only I0 survives -> transverse field is the pure Jones vec
    assert abs(np.sqrt(np.sum(np.abs(E0) ** 2)) - 1.0) < 1e-12
    assert abs(E0[1] - 1j * E0[0]) < 1e-12         # circular relationship Ey = i Ex


def test_unit_amplitude_all_polarizations():
    """|E| = 1 at the geometric focus for linear, circular, and elliptical Jones."""
    for pol in [(1, 0), (0, 1), (1, 1j), (1, -1j), (1, 1), (0.6, 0.8j)]:
        rw = RichardsWolfFocus(MED, NA=0.63, polarization=pol)
        E0 = rw.focal_field(np.zeros((1, 3)))[0]
        assert abs(np.sqrt(np.sum(np.abs(E0) ** 2)) - 1.0) < 1e-12


def test_linearity_superposition():
    """focal_field(px,py) == px*E_x + py*E_y for the (un-normalized) inputs."""
    rw = RichardsWolfFocus(MED, NA=0.63, polarization=(1.0, 1j))
    pts = _pts()
    Ex = rw._focal_xyz(pts)
    Ey = rw._focal_xyz_y(pts)
    # rw.polarization is normalized to (1, 1j)/sqrt2
    px, py = rw.polarization
    manual = px * Ex + py * Ey
    assert np.allclose(rw.focal_field(pts), manual)


def test_ypol_is_xpol_rotated_90deg():
    """+y input is the +x field rotated 90 deg about z (exact, no approximation)."""
    rw = RichardsWolfFocus(MED, NA=0.63, polarization=(1.0, 0.0))
    E0y = rw._focal_xyz_y(np.zeros((1, 3)))[0]
    # y-input focal field at origin is along +y with the same prefactor as x along +x
    E0x = rw._focal_xyz(np.zeros((1, 3)))[0]
    assert abs(E0y[1] - E0x[0]) < 1e-12
    assert abs(E0y[0]) < 1e-12 and abs(E0y[2]) < 1e-12


def test_circular_is_azimuthally_symmetric_intensity_on_ring():
    """Circular polarization gives an azimuthally-symmetric |E|^2 on a transverse ring
    at the focal plane (a signature that distinguishes it from linear, which has the
    cos-2phi lobes of I2)."""
    rw = RichardsWolfFocus(MED, NA=0.63, polarization=(1.0, 1j))
    r = 0.3e-6
    phi = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    pts = np.stack([r * np.cos(phi), r * np.sin(phi), np.zeros_like(phi)], axis=-1)
    I = np.sum(np.abs(rw.focal_field(pts)) ** 2, axis=-1)
    # linear would vary strongly in phi; circular is symmetric to a tight tolerance
    assert (I.max() - I.min()) / I.mean() < 1e-6
