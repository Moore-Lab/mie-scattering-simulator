"""Objective primitives for the detection optimizer: cached-field scoring and the
finite-difference sensitivity dict (INTERFACES.md §7; PHYSICS.md §4.2-4.5, §7).

Time convention e^{-iωt}. All fields are computed by the frozen PlaneWaveProvider.
"""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.detect.optics import CollectionGeometry
from mieinfo.detect.schemes import apply_scheme
from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider, field_derivative
from mieinfo.info.modes import combine_direction, information_pattern
from mieinfo.optimize.objective import (
    GeometryCache,
    build_geometry_cache,
    evaluate_geometry,
    score_geometry,
    sensitivity,
)
from mieinfo.types import AngularGrid, Medium, Sphere

MED = Medium(n=1.0, wavelength_vacuum_m=532e-9)
SPH = Sphere(radius_m=1e-6, m=1.46 + 0j)   # x ~ 11.8
PROV = PlaneWaveProvider()
GRID = AngularGrid.full_sphere(400, 96)

Z = np.array([0.0, 0.0, 1.0])
X = np.array([1.0, 0.0, 0.0])


def test_cache_is_geometry_independent_and_reused():
    # Building the cache solves the field ONCE; scoring many geometries reuses it.
    cache = build_geometry_cache(PROV, SPH, PlaneWave(MED), Z, GRID)
    assert isinstance(cache, GeometryCache)
    assert cache.pattern.f_total > 0.0
    assert cache.q_sca > 0.0
    # The field object handed to every geometry is literally the same array (no re-solve).
    f_id = id(cache.field)
    for na in (0.3, 0.6, 0.9):
        for direction in ("forward", "backward"):
            g = CollectionGeometry(direction=direction, NA=na, lo_mode="optimal")
            res = evaluate_geometry(cache, g)
            assert 0.0 <= res.eta_q <= 1.0 + 1e-12
    assert id(cache.field) == f_id


def test_score_matches_direct_apply_scheme():
    # score_geometry via the cache must equal a from-scratch apply_scheme call.
    n = X
    cache = build_geometry_cache(PROV, SPH, PlaneWave(MED), n, GRID)
    g = CollectionGeometry(direction="forward", NA=0.8, lo_mode="optimal")

    pat = information_pattern(PROV, GRID, SPH, PlaneWave(MED), np.zeros(3), n)
    deriv = field_derivative(PROV, GRID, SPH, PlaneWave(MED), np.zeros(3), method="analytic")
    dq = combine_direction(deriv, n)
    field = PROV.field(GRID, SPH, PlaneWave(MED), np.zeros(3))
    direct = apply_scheme(field, dq, pat, g, PROV.q_sca(SPH, PlaneWave(MED)))

    assert score_geometry(cache, g) == pytest.approx(direct.eta_q, rel=1e-12)


def test_eta_increases_with_NA_for_optimal_lo():
    # A larger cone collects more information with the optimal LO (monotone in NA).
    cache = build_geometry_cache(PROV, SPH, PlaneWave(MED), X, GRID)
    etas = [score_geometry(cache, CollectionGeometry(direction="forward", NA=na, lo_mode="optimal"))
            for na in (0.2, 0.4, 0.6, 0.8)]
    for a, b in zip(etas, etas[1:]):
        assert b >= a - 1e-9


def test_sensitivity_dict_is_populated():
    # INTERFACES §7: sensitivity carries d eta / d {radius, NA, waist, m}, all finite.
    g = CollectionGeometry(direction="backward", NA=0.8, lo_mode="optimal")
    sens = sensitivity(PROV, SPH, PlaneWave(MED), Z, GRID, g)
    assert set(sens.keys()) == {"radius", "NA", "waist", "m"}
    for k, v in sens.items():
        assert np.isfinite(v), (k, v)


def test_sensitivity_NA_positive_for_optimal_lo():
    # d eta / d NA > 0: widening the optimal-LO cone collects strictly more info. The
    # step is grid-aware so the hard-cone quantization (VALIDATION §5) does not zero it.
    g = CollectionGeometry(direction="backward", NA=0.6, lo_mode="optimal")
    sens = sensitivity(PROV, SPH, PlaneWave(MED), Z, GRID, g)
    assert sens["NA"] > 1e-3


def test_sensitivity_waist_zero_for_plane_wave():
    # A plane wave is waist-invariant (waist = inf); d eta / d waist is exactly 0.
    assert not np.isfinite(PlaneWave(MED).waist_m())
    g = CollectionGeometry(direction="forward", NA=0.8, lo_mode="optimal")
    sens = sensitivity(PROV, SPH, PlaneWave(MED), X, GRID, g)
    assert sens["waist"] == 0.0


def test_sensitivity_grid_common_mode_is_stable():
    # Because the grid is held common-mode across the ± evaluations, the NA sensitivity
    # is a smooth finite number rather than a spurious 0 or a huge grid-quantization spike.
    g = CollectionGeometry(direction="forward", NA=0.6, lo_mode="optimal")
    s1 = sensitivity(PROV, SPH, PlaneWave(MED), X, GRID, g)
    s2 = sensitivity(PROV, SPH, PlaneWave(MED), X, GRID, g)
    # Deterministic: same grid, same inputs -> identical result.
    assert s1["NA"] == pytest.approx(s2["NA"], rel=1e-12)
    assert s1["radius"] == pytest.approx(s2["radius"], rel=1e-12)
