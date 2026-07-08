"""G-GOLD (information pattern): the full pipeline (PlaneWaveProvider -> field_derivative
-> combine_direction -> info_density -> collection_efficiency) reproduces the operative
532 nm golden set (VALIDATION.md §2). @slow — nightly/full lane only (matching grids up
to 1400x420). Time convention e^{-iωt}.

Structure (forward fraction, peak angle) reproduces to machine precision. η reproduces to
~4e-3: the golden η were read off an 80-point NA interpolation (prototype/study_detection_532
eta_curve), whereas the package computes η exactly at each NA — so the residual is the
golden's interpolation artifact, not a package error, and the package value is the more
accurate one."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from mieinfo.detect.optics import CollectionGeometry, collection_efficiency
from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider
from mieinfo.info.modes import information_pattern
from mieinfo.types import AngularGrid, Medium, Sphere

GOLD = Path(__file__).resolve().parents[2] / "data" / "golden" / "information_pattern_results.json"
CASES = json.loads(GOLD.read_text())
MED = Medium(n=1.0, wavelength_vacuum_m=532e-9)


def _grid_for(a):
    return (1400 if a >= 12 else (1000 if a >= 8 else 700), 420)


def _forward_fraction(p):
    fwd = p.grid.theta < np.pi / 2
    W = p.grid.w_solid
    return float(np.sum(p.density[fwd] * W[fwd]) / np.sum(p.density * W))


def _peak_deg(p):
    fp = np.sum(p.density * p.grid.w_solid, axis=1)
    return float(np.degrees(p.grid.theta[np.argmax(fp)]))


@pytest.mark.slow
@pytest.mark.parametrize("r", CASES, ids=[f"a{r['a_um']:g}" for r in CASES])
def test_reproduces_golden_information_pattern(r):
    a = r["a_um"]
    sph = Sphere(radius_m=a * 1e-6, m=1.46 + 0j)
    nt, nph = _grid_for(a)
    grid = AngularGrid.full_sphere(nt, nph)
    beam = PlaneWave(MED)
    prov = PlaneWaveProvider()
    px = information_pattern(prov, grid, sph, beam, np.zeros(3), np.array([1.0, 0, 0]))
    pz = information_pattern(prov, grid, sph, beam, np.zeros(3), np.array([0, 0, 1.0]))

    # structure reproduces to machine precision
    assert abs(_forward_fraction(px) - r["info_x_fwd"]) < 1e-9
    assert abs(_forward_fraction(pz) - r["info_z_fwd"]) < 1e-9
    assert abs(_peak_deg(px) - r["x_peak_deg"]) < 1e-6
    assert abs(_peak_deg(pz) - r["z_peak_deg"]) < 1e-6

    # efficiencies reproduce within the golden's NA-interpolation artifact (see module doc)
    for axis, p in (("x", px), ("z", pz)):
        for direc, tag in (("forward", "for"), ("backward", "bac")):
            for na in (0.5, 0.8, 0.95):
                eta = collection_efficiency(p, CollectionGeometry(direction=direc, NA=na))
                assert abs(eta - r[f"eta_{axis}_{tag}_NA{na}"]) < 6e-3
