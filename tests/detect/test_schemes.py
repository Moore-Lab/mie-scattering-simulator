"""Detection schemes: LO mode-overlap κ and apply_scheme (PHYSICS.md §4.2-4.3, §5)."""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.detect.optics import CollectionGeometry, DetectionResult, collection_efficiency
from mieinfo.detect.schemes import apply_scheme, lo_overlap
from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider, field_derivative
from mieinfo.info.modes import combine_direction, information_pattern
from mieinfo.types import AngularGrid, Medium, Sphere

MED = Medium(n=1.0, wavelength_vacuum_m=532e-9)
SPH = Sphere(radius_m=1e-6, m=1.46 + 0j)          # x ~ 11.8
PROV = PlaneWaveProvider()
BEAM = PlaneWave(MED)


def _setup(axis, ntheta=400, nphi=128):
    grid = AngularGrid.full_sphere(ntheta, nphi)
    n_hat = np.asarray(axis, dtype=float)
    pattern = information_pattern(PROV, grid, SPH, BEAM, np.zeros(3), n_hat)
    deriv = field_derivative(PROV, grid, SPH, BEAM, np.zeros(3), method="analytic")
    deriv_q = combine_direction(deriv, n_hat)
    field = PROV.field(grid, SPH, BEAM, np.zeros(3))
    q_sca = PROV.q_sca(SPH, BEAM)
    return field, deriv_q, pattern, q_sca


ALL_LO = ("optimal", "gaussian", "split", "quadrant", "self_homodyne")
ALL_DIR = ("forward", "backward", "both")


def test_result_type_and_fields():
    field, deriv_q, pattern, q_sca = _setup([1, 0, 0])
    g = CollectionGeometry(direction="forward", NA=0.8, lo_mode="optimal")
    r = apply_scheme(field, deriv_q, pattern, g, q_sca)
    assert isinstance(r, DetectionResult)
    assert r.geometry is g
    assert r.gamma_ba_rel == pytest.approx(q_sca)   # γ_ba ~ Q_sca


def test_eta_q_in_unit_interval_all_schemes_all_directions():
    for axis in ([1, 0, 0], [0, 0, 1]):
        field, deriv_q, pattern, q_sca = _setup(axis)
        for lo in ALL_LO:
            for d in ALL_DIR:
                for na in (0.3, 0.8, 0.95):
                    g = CollectionGeometry(direction=d, NA=na, lo_mode=lo)
                    r = apply_scheme(field, deriv_q, pattern, g, q_sca)
                    assert 0.0 <= r.eta_q <= 1.0 + 1e-12, (lo, d, na, axis, r.eta_q)


def test_kappa_in_unit_interval():
    field, deriv_q, _, _ = _setup([1, 0, 0])
    for lo in ALL_LO:
        for d in ("forward", "backward"):
            g = CollectionGeometry(direction=d, NA=0.8, lo_mode=lo)
            k = lo_overlap(field, deriv_q, g)
            assert 0.0 <= k <= 1.0 + 1e-12, (lo, d, k)


def test_optimal_lo_kappa_is_one():
    # Optimal LO ∝ ∂E/∂q saturates F_q(Ω): κ = 1 exactly.
    field, deriv_q, _, _ = _setup([1, 0, 0])
    for d in ("forward", "backward", "both"):
        g = CollectionGeometry(direction=d, NA=0.8, lo_mode="optimal")
        assert lo_overlap(field, deriv_q, g) == pytest.approx(1.0, abs=1e-10)


def test_optimal_eta_equals_collection_efficiency():
    # With the optimal LO (κ=1), apply_scheme's eta_q is exactly the cone efficiency.
    for axis in ([1, 0, 0], [0, 0, 1]):
        field, deriv_q, pattern, q_sca = _setup(axis)
        for d in ("forward", "backward", "both"):
            for na in (0.3, 0.8):
                g = CollectionGeometry(direction=d, NA=na, lo_mode="optimal")
                r = apply_scheme(field, deriv_q, pattern, g, q_sca)
                assert r.eta_q == pytest.approx(collection_efficiency(pattern, g), rel=1e-12)


def test_self_homodyne_shortfall_vs_optimal():
    # Self-homodyne (probe as its own LO) never beats the optimal LO for the same
    # geometry, and shows a genuine (nonzero) shortfall on the axes/directions it serves.
    saw_strict_shortfall = False
    for axis in ([1, 0, 0], [0, 0, 1]):
        field, deriv_q, pattern, q_sca = _setup(axis)
        for d in ("forward", "backward"):
            g_opt = CollectionGeometry(direction=d, NA=0.8, lo_mode="optimal")
            g_self = CollectionGeometry(direction=d, NA=0.8, lo_mode="self_homodyne")
            eta_opt = apply_scheme(field, deriv_q, pattern, g_opt, q_sca).eta_q
            eta_self = apply_scheme(field, deriv_q, pattern, g_self, q_sca).eta_q
            assert eta_self <= eta_opt + 1e-12                  # never exceeds optimal
            if eta_opt > 1e-3:                                  # a served configuration
                assert eta_self < eta_opt                       # strict shortfall
                assert eta_self > 0.0
                saw_strict_shortfall = True
    assert saw_strict_shortfall


def test_self_homodyne_direction_emphasis():
    # PHYSICS §5/§8: forward emphasizes transverse axes; backward recovers the collinear
    # axis. Carried by the cone efficiency, so it survives into the self-homodyne eta_q.
    fx, dqx, px, qs = _setup([1, 0, 0])   # transverse
    fz, dqz, pz, qsz = _setup([0, 0, 1])  # collinear (axial)
    gf = CollectionGeometry(direction="forward", NA=0.8, lo_mode="self_homodyne")
    gb = CollectionGeometry(direction="backward", NA=0.8, lo_mode="self_homodyne")
    eta_x_fwd = apply_scheme(fx, dqx, px, gf, qs).eta_q
    eta_x_bwd = apply_scheme(fx, dqx, px, gb, qs).eta_q
    eta_z_fwd = apply_scheme(fz, dqz, pz, gf, qsz).eta_q
    eta_z_bwd = apply_scheme(fz, dqz, pz, gb, qsz).eta_q
    assert eta_x_fwd > eta_x_bwd     # transverse -> forward
    assert eta_z_bwd > eta_z_fwd     # collinear  -> backward


def test_split_and_gaussian_are_azimuthal_parity_complementary():
    # A symmetric Gaussian bucket reads the axial (symmetric) signal; an antisymmetric
    # split reads the transverse (antisymmetric) signal. Each is parity-blind to the
    # other — a real single-mode limitation (documented in schemes.py).
    fx, dqx, px, qs = _setup([1, 0, 0])   # transverse -> antisymmetric far-field response
    fz, dqz, pz, qsz = _setup([0, 0, 1])  # axial -> symmetric response
    g_split = CollectionGeometry(direction="forward", NA=0.8, lo_mode="split")
    g_gauss = CollectionGeometry(direction="forward", NA=0.8, lo_mode="gaussian")
    kx_split = lo_overlap(fx, dqx, g_split)
    kx_gauss = lo_overlap(fx, dqx, g_gauss)
    kz_split = lo_overlap(fz, dqz, g_split)
    kz_gauss = lo_overlap(fz, dqz, g_gauss)
    assert kx_split > kx_gauss    # transverse favours split over Gaussian
    assert kz_gauss > kz_split    # axial favours Gaussian over split


def test_apodization_none_vs_aplanatic_both_valid():
    # Both apodization options run and stay in [0,1]; optimal-LO eta_q (κ=1) is
    # apodization-independent (the optimal LO reshapes to the delivered field).
    field, deriv_q, pattern, q_sca = _setup([1, 0, 0])
    g_apl = CollectionGeometry(direction="forward", NA=0.8, apodization="aplanatic",
                               lo_mode="optimal")
    g_none = CollectionGeometry(direction="forward", NA=0.8, apodization="none",
                                lo_mode="optimal")
    r_apl = apply_scheme(field, deriv_q, pattern, g_apl, q_sca)
    r_none = apply_scheme(field, deriv_q, pattern, g_none, q_sca)
    assert r_apl.eta_q == pytest.approx(r_none.eta_q, rel=1e-12)  # optimal LO: no apod dependence
    # A fixed-mode scheme's κ *does* respond to apodization (sanity: still in range).
    for apo in ("aplanatic", "none"):
        g = CollectionGeometry(direction="forward", NA=0.8, apodization=apo,
                               lo_mode="self_homodyne")
        assert 0.0 <= lo_overlap(field, deriv_q, g) <= 1.0


def test_unknown_lo_mode_raises():
    field, deriv_q, _, _ = _setup([1, 0, 0])
    g = CollectionGeometry(direction="forward", NA=0.8, lo_mode="bogus")
    with pytest.raises(ValueError):
        lo_overlap(field, deriv_q, g)
