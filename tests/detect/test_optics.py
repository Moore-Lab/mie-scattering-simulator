"""Collection cones and efficiency (INTERFACES.md §6, PHYSICS.md §4.2)."""
from __future__ import annotations

import numpy as np

from mieinfo.detect.optics import CollectionGeometry, collection_efficiency, cone_mask
from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider
from mieinfo.info.modes import information_pattern
from mieinfo.types import AngularGrid, Medium, Sphere


def _pattern(axis):
    med = Medium(n=1.0, wavelength_vacuum_m=532e-9)
    sph = Sphere(radius_m=1e-6, m=1.46 + 0j)         # x ~ 11.8
    grid = AngularGrid.full_sphere(300, 96)
    return information_pattern(PlaneWaveProvider(), grid, sph, PlaneWave(med),
                               np.zeros(3), np.asarray(axis, dtype=float))


def test_cone_mask_forward_backward():
    grid = AngularGrid.full_sphere(180, 32)
    fwd = cone_mask(grid, CollectionGeometry(direction="forward", NA=np.sin(0.3)))
    bwd = cone_mask(grid, CollectionGeometry(direction="backward", NA=np.sin(0.3)))
    TH = np.broadcast_to(grid.theta[:, None], grid.w_solid.shape)
    assert np.all(TH[fwd] <= 0.3 + 1e-12)
    assert np.all(TH[bwd] >= np.pi - 0.3 - 1e-12)
    assert not np.any(fwd & bwd)


def test_efficiency_in_unit_interval_and_full_is_one():
    pat = _pattern([1, 0, 0])
    for na in (0.2, 0.5, 0.8, 0.95):
        for d in ("forward", "backward", "both"):
            eta = collection_efficiency(pat, CollectionGeometry(direction=d, NA=na))
            assert 0.0 <= eta <= 1.0
    assert abs(collection_efficiency(pat, CollectionGeometry(direction="both", NA=1.0)) - 1.0) < 1e-12


def test_efficiency_monotonic_in_na():
    pat = _pattern([1, 0, 0])
    etas = [collection_efficiency(pat, CollectionGeometry(direction="forward", NA=na))
            for na in (0.3, 0.5, 0.7, 0.9)]
    assert all(b >= a - 1e-12 for a, b in zip(etas, etas[1:]))


def test_transverse_forward_beats_backward_and_axial_backward_beats_forward():
    px, pz = _pattern([1, 0, 0]), _pattern([0, 0, 1])
    g_fwd = CollectionGeometry(direction="forward", NA=0.8)
    g_bwd = CollectionGeometry(direction="backward", NA=0.8)
    assert collection_efficiency(px, g_fwd) > collection_efficiency(px, g_bwd)   # transverse -> forward
    assert collection_efficiency(pz, g_bwd) > collection_efficiency(pz, g_fwd)   # axial -> backward
