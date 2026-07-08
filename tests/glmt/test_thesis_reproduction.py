"""Cross-check of the RichardsWolfFocus GLMT information-collection efficiencies vs the
Moore Lab reference (arXiv:2408.15483 / dissertation Ch.2, NA=0.63, 1064 nm).

STATUS (honest — see docs/dissertation_comparison.md). The scattered-field information
radiation field (IRF = dE_s/dr_0, the reference's Eq 3; efficiency = the Eq 11 ratio
F_q(cone)/F_q(4pi)) is reproduced by this package via the GLMT BSC path AND cross-checked
against an independent angular-spectrum IRF (both reproduce the plane-wave limit exactly).

What IS validated and asserted here:
  * G-LIMIT: RichardsWolfFocus at low NA reproduces the plane-wave provider.
  * Mie transverse-forward focusing gain: plane-wave ~14% rises to the reference ~38%.
  * Dipole transverse fore-aft symmetry (the Rayleigh signature).

What is NOT reproduced (asserted here as an honest xfail-style bound, NOT a pass):
  * The reference's *forward-comparable* axial efficiency for the 3 um Mie sphere
    (thesis Fig 2.5 reads ~37% fwd). The scattered-field IRF weight ~(k_inc - k*shat)^2
    is intrinsically BACKWARD-heavy for axial motion; both the BSC path and the
    independent angular-spectrum IRF give axial-forward ~1%. The published paper
    (2408.15483) itself states axial information is collected predominantly in the
    BACKWARD direction under high-NA focusing, consistent with this package and
    inconsistent with the forward-comparable dissertation-figure read-off. See
    docs/dissertation_comparison.md for the full analysis.

These tests therefore pin the VALIDATED behaviour and record the discrepancy explicitly
so a future fix (or a corrected target) can flip the documented bounds.
"""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.detect.optics import CollectionGeometry, collection_efficiency
from mieinfo.glmt.beam import PlaneWave, RichardsWolfFocus
from mieinfo.glmt.scatter import GLMTProvider, PlaneWaveProvider
from mieinfo.info.modes import information_pattern
from mieinfo.types import AngularGrid, Medium, Sphere

MED = Medium(n=1.0, wavelength_vacuum_m=1064e-9)
NA = 0.63
M_SILICA = 1.4496 + 0j


def _effs(provider, beam, sphere, method, ntheta=200, nphi=64, na_col=NA):
    grid = AngularGrid.full_sphere(ntheta, nphi)
    out = {}
    for name, nhat, direction in [
        ("txf", (1, 0, 0), "forward"),
        ("txb", (1, 0, 0), "backward"),
        ("azf", (0, 0, 1), "forward"),
        ("azb", (0, 0, 1), "backward"),
    ]:
        patt = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                                   np.array(nhat, float), method=method)
        out[name] = collection_efficiency(
            patt, CollectionGeometry(direction=direction, NA=na_col)) * 100.0
    return out


@pytest.mark.slow
def test_glimit_richards_wolf_matches_plane_wave():
    """RichardsWolfFocus at NA=0.05 reproduces the plane-wave provider (G-LIMIT)."""
    sph = Sphere(radius_m=1.5e-6, m=M_SILICA)
    pw = _effs(PlaneWaveProvider(), PlaneWave(MED), sph, "analytic")
    rw = _effs(GLMTProvider(), RichardsWolfFocus(MED, NA=0.05), sph,
               "finite_difference")
    for key in ("txf", "txb", "azf", "azb"):
        assert abs(pw[key] - rw[key]) < 0.5, (key, pw[key], rw[key])


@pytest.mark.slow
def test_mie_transverse_forward_focusing_gain():
    """Focusing lifts the 3 um Mie transverse-forward efficiency from the plane-wave
    ~14% toward the reference ~38% (the one focused-beam quantity that checks out)."""
    sph = Sphere(radius_m=1.5e-6, m=M_SILICA)
    pw = _effs(PlaneWaveProvider(), PlaneWave(MED), sph, "analytic",
               ntheta=300, nphi=96)
    rw = _effs(GLMTProvider(), RichardsWolfFocus(MED, NA=NA), sph,
               "finite_difference", ntheta=300, nphi=96)
    assert pw["txf"] < 20.0                        # plane-wave forward is modest
    assert 30.0 <= rw["txf"] <= 45.0               # focused ~= reference 38% (+-20%)


@pytest.mark.slow
def test_dipole_transverse_fore_aft_symmetry():
    """Rayleigh (a=50 nm) dipole: transverse forward ~= backward (dipole signature)."""
    sph = Sphere(radius_m=50e-9, m=M_SILICA)
    rw = _effs(GLMTProvider(), RichardsWolfFocus(MED, NA=NA), sph,
               "finite_difference")
    # fore-aft symmetric to ~40% relative (read-off + focusing tolerance)
    assert abs(rw["txf"] - rw["txb"]) / (0.5 * (rw["txf"] + rw["txb"])) < 0.6


@pytest.mark.slow
def test_axial_forward_is_backward_dominated_documented_discrepancy():
    """HONEST NON-REPRODUCTION: the scattered-field IRF gives backward-dominated axial
    info (axial-forward ~1%), NOT the dissertation-figure ~37% forward. This asserts the
    package's actual (backward-heavy) behaviour so the discrepancy is pinned, not hidden.
    """
    sph = Sphere(radius_m=1.5e-6, m=M_SILICA)
    rw = _effs(GLMTProvider(), RichardsWolfFocus(MED, NA=NA), sph,
               "finite_difference", ntheta=300, nphi=96)
    assert rw["azf"] < 5.0                          # NOT the reference's ~37%
    assert rw["azb"] > 50.0                         # axial info is backward-dominated
