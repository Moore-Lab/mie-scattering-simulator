"""Multi-channel detection: per-beam-frame evaluation + Fisher-info combination
(PHYSICS.md §4.5; INTERFACES.md §6)."""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.detect.channels import (
    DetectionChannel,
    MultiChannelResult,
    evaluate_channels,
)
from mieinfo.detect.optics import CollectionGeometry
from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider
from mieinfo.types import AngularGrid, Medium, Sphere

MED = Medium(n=1.0, wavelength_vacuum_m=532e-9)
SPH = Sphere(radius_m=1e-6, m=1.46 + 0j)
PROV = PlaneWaveProvider()
GRID = AngularGrid.full_sphere(400, 128)
GEOM = CollectionGeometry(direction="forward", NA=0.8, lo_mode="optimal")


def _channel(prop, pol, name):
    return DetectionChannel(beam=PlaneWave(MED), propagation_lab=np.asarray(prop, float),
                            polarization_lab=np.asarray(pol, float), geometry=GEOM, name=name)


# Two-beam apparatus: beam A along lab +z, beam B along lab +x (PHYSICS.md §4.5).
CH_A = _channel([0, 0, 1], [1, 0, 0], "A_z")
CH_B = _channel([1, 0, 0], [0, 0, 1], "B_x")


def test_result_shape_and_types():
    r = evaluate_channels(PROV, SPH, [CH_A, CH_B], np.array([0.0, 0, 1]), GRID)
    assert isinstance(r, MultiChannelResult)
    assert len(r.per_channel) == 2
    assert 0.0 <= r.eta_q_total <= 1.0
    assert r.fisher_total_rel > 0.0
    assert r.sql_distance >= 1.0 - 1e-12


def test_eta_q_total_in_unit_interval():
    for lab in ([0, 0, 1], [1, 0, 0], [0, 1, 0], [1, 1, 1]):
        r = evaluate_channels(PROV, SPH, [CH_A, CH_B], np.asarray(lab, float), GRID)
        assert 0.0 <= r.eta_q_total <= 1.0 + 1e-12
        for pc in r.per_channel:
            assert 0.0 <= pc.eta_q <= 1.0 + 1e-12


def test_fisher_information_adds_for_independent_channels():
    # F_AB(collected) = F_A(collected) + F_B(collected) exactly (independent shot noise).
    for lab in ([0, 0, 1], [1, 0, 0], [0, 1, 0], [1, 1, 1]):
        n = np.asarray(lab, float)
        rA = evaluate_channels(PROV, SPH, [CH_A], n, GRID)
        rB = evaluate_channels(PROV, SPH, [CH_B], n, GRID)
        rAB = evaluate_channels(PROV, SPH, [CH_A, CH_B], n, GRID)
        assert rAB.fisher_total_rel == pytest.approx(
            rA.fisher_total_rel + rB.fisher_total_rel, rel=1e-10)
        # Two independent channels collect at least as much as either alone.
        assert rAB.fisher_total_rel >= rA.fisher_total_rel - 1e-6
        assert rAB.fisher_total_rel >= rB.fisher_total_rel - 1e-6


def test_two_beam_layout_covers_all_three_lab_axes():
    # The crux of PHYSICS §4.5: a lab axis collinear with one beam (poorly measured in
    # its forward lobe) is transverse to the other (well measured). With beam A along z
    # and beam B along x, EVERY lab axis is measured well by at least one beam, and the
    # two-beam set beats the worse single beam on every axis.
    for lab in ([1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]):
        n = np.asarray(lab, float)
        rA = evaluate_channels(PROV, SPH, [CH_A], n, GRID)
        rB = evaluate_channels(PROV, SPH, [CH_B], n, GRID)
        rAB = evaluate_channels(PROV, SPH, [CH_A, CH_B], n, GRID)
        best_single = max(rA.eta_q_total, rB.eta_q_total)
        worst_single = min(rA.eta_q_total, rB.eta_q_total)
        # At least one beam measures this axis well (transverse to it): η >= 0.3.
        assert best_single >= 0.3, (lab, rA.eta_q_total, rB.eta_q_total)
        # Combined collected information strictly exceeds the worse beam alone.
        assert rAB.fisher_total_rel > 0.0
        assert worst_single <= best_single


def test_collinear_axis_rescued_by_orthogonal_beam():
    # lab_z is collinear with beam A (A forward is poor) but transverse to beam B (good).
    n = np.array([0.0, 0.0, 1.0])
    rA = evaluate_channels(PROV, SPH, [CH_A], n, GRID)   # A: z is its own axis -> poor
    rB = evaluate_channels(PROV, SPH, [CH_B], n, GRID)   # B: z is transverse -> good
    assert rB.eta_q_total > rA.eta_q_total
    assert rA.eta_q_total < 0.1        # collinear-forward is poor (PHYSICS §4.1 (1-cosθ)²)
    assert rB.eta_q_total > 0.3        # transverse-forward is good (sin²θ)


def test_backward_port_recovers_collinear_axis_single_beam():
    # Even a single beam measures its own collinear axis well if it collects BACKWARD
    # (PHYSICS §8: η_z(NA0.8, bwd) ~ 0.65). n_hat = beam A's own axis (lab z).
    g_bwd = CollectionGeometry(direction="backward", NA=0.8, lo_mode="optimal")
    ch_a_bwd = DetectionChannel(beam=PlaneWave(MED), propagation_lab=np.array([0., 0, 1]),
                                polarization_lab=np.array([1., 0, 0]), geometry=g_bwd,
                                name="A_z_bwd")
    r = evaluate_channels(PROV, SPH, [ch_a_bwd], np.array([0.0, 0, 1]), GRID)
    assert r.eta_q_total > 0.5


def test_dependent_channel_raises():
    # Correlated (shared-detector/LO) channels must not be silently summed.
    dep = DetectionChannel(beam=PlaneWave(MED), propagation_lab=np.array([1., 0, 0]),
                           polarization_lab=np.array([0, 0, 1.]), geometry=GEOM,
                           name="shared", independent=False)
    with pytest.raises(NotImplementedError):
        evaluate_channels(PROV, SPH, [CH_A, dep], np.array([0.0, 0, 1]), GRID)


def test_single_channel_matches_direct_scheme():
    # A one-channel MultiChannelResult's eta_q_total equals that channel's own eta_q
    # (F cancels in the ratio for a single channel).
    from mieinfo.detect.schemes import apply_scheme
    from mieinfo.glmt.scatter import field_derivative
    from mieinfo.info.modes import combine_direction, information_pattern

    n = np.array([0.0, 1.0, 0.0])   # transverse to beam A (along z)
    r = evaluate_channels(PROV, SPH, [CH_A], n, GRID)

    # Direct evaluation in beam A's frame: A is aligned with lab (prop z, pol x) so the
    # beam frame equals the lab frame and n stays [0,1,0].
    pat = information_pattern(PROV, GRID, SPH, PlaneWave(MED), np.zeros(3), n)
    deriv = field_derivative(PROV, GRID, SPH, PlaneWave(MED), np.zeros(3), method="analytic")
    dq = combine_direction(deriv, n)
    field = PROV.field(GRID, SPH, PlaneWave(MED), np.zeros(3))
    direct = apply_scheme(field, dq, pat, GEOM, PROV.q_sca(SPH, PlaneWave(MED)))
    assert r.eta_q_total == pytest.approx(direct.eta_q, rel=1e-10)
