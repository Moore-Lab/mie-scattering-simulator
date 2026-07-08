"""Detection-geometry optimizer: single-channel ranking, the axial/transverse
direction pick, the multi-channel two-beam win, and the sensitivity report
(INTERFACES.md §7; PHYSICS.md §4.1-4.5, §8; MASTER_PLAN F2/F3).

Time convention e^{-iωt}. Fields come from the frozen PlaneWaveProvider.
"""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.detect.channels import (
    DetectionChannel,
    MultiChannelResult,
    evaluate_channels,
)
from mieinfo.detect.optics import CollectionGeometry, DetectionResult
from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider
from mieinfo.optimize.search import (
    Constraints,
    OptResult,
    _best_channel_set,
    _channel_candidates,
    _na_grid,
    optimize_detection,
)
from mieinfo.types import AngularGrid, Medium, Sphere

MED = Medium(n=1.0, wavelength_vacuum_m=532e-9)
SPH = Sphere(radius_m=1e-6, m=1.46 + 0j)   # x ~ 11.8
PROV = PlaneWaveProvider()
GRID = AngularGrid.full_sphere(400, 96)

Z = np.array([0.0, 0.0, 1.0])
X = np.array([1.0, 0.0, 0.0])

# Forward/backward only, optimal LO: an unambiguous direction pick per axis.
C_FB = Constraints(na_max=0.8, directions_allowed=("forward", "backward"),
                   schemes_allowed=("optimal",))


# --------------------------------------------------------------------------- #
# Single-channel ranking and the axial/transverse direction pick
# --------------------------------------------------------------------------- #

def test_result_shape():
    r = optimize_detection(PROV, SPH, PlaneWave(MED), Z, C_FB, GRID)
    assert isinstance(r, OptResult)
    assert len(r.ranked) == len(_na_grid(C_FB.na_max)) * 2   # 2 directions x NA grid
    g, dr = r.best
    assert isinstance(g, CollectionGeometry)
    assert isinstance(dr, DetectionResult)
    assert r.best == r.ranked[0]
    # single-channel result: no channel set unless max_channels > 1
    assert r.best_channel_set is None
    assert r.best_multichannel is None


def test_ranked_sorted_by_eta_descending():
    for n in (Z, X, np.array([0.0, 1.0, 0.0])):
        r = optimize_detection(PROV, SPH, PlaneWave(MED), n, C_FB, GRID)
        etas = [dr.eta_q for _, dr in r.ranked]
        assert all(etas[i] >= etas[i + 1] - 1e-12 for i in range(len(etas) - 1)), etas
        assert r.best[1].eta_q == max(etas)


def test_axial_picks_backward():
    # PHYSICS §4.1/§8, MASTER_PLAN F2: axial (beam-collinear) information is
    # BACKWARD-weighted ((1-cosθ)² peaks in backscatter). For n_hat = z about a +z beam
    # the optimizer must pick 'backward' collection.
    r = optimize_detection(PROV, SPH, PlaneWave(MED), Z, C_FB, GRID)
    assert r.best[0].direction == "backward"
    # And backward genuinely beats the best forward candidate for this axis.
    best_fwd = max((dr.eta_q for g, dr in r.ranked if g.direction == "forward"))
    best_bwd = max((dr.eta_q for g, dr in r.ranked if g.direction == "backward"))
    assert best_bwd > best_fwd
    # Backward for the collinear axis is the ~0.65 regime of PHYSICS §8 at NA 0.8.
    assert r.best[1].eta_q > 0.5


def test_transverse_picks_forward():
    # PHYSICS §4.1/§8: transverse info has sin²θcos²φ weight — forward-weighted, near
    # zero backward. For n_hat = x about a +z beam the optimizer must pick 'forward'.
    r = optimize_detection(PROV, SPH, PlaneWave(MED), X, C_FB, GRID)
    assert r.best[0].direction == "forward"
    best_fwd = max((dr.eta_q for g, dr in r.ranked if g.direction == "forward"))
    best_bwd = max((dr.eta_q for g, dr in r.ranked if g.direction == "backward"))
    assert best_fwd > best_bwd
    # Backward is useless for transverse (PHYSICS §8: η_x(bwd) ~ 0.02).
    assert best_bwd < 0.1


def test_best_na_is_largest_allowed_for_optimal_lo():
    # With the optimal LO η is monotone in NA, so the recommended NA is the cap.
    r = optimize_detection(PROV, SPH, PlaneWave(MED), Z, C_FB, GRID)
    assert r.best[0].NA == pytest.approx(max(_na_grid(C_FB.na_max)))


def test_fixed_sphere_overrides_argument_sphere():
    # constraints.fixed_sphere, if given, is what gets optimized (INTERFACES §7).
    other = Sphere(radius_m=2e-6, m=1.50 + 0j)
    C = Constraints(na_max=0.8, directions_allowed=("forward", "backward"),
                    schemes_allowed=("optimal",), fixed_sphere=other)
    r_fixed = optimize_detection(PROV, SPH, PlaneWave(MED), Z, C, GRID)
    r_direct = optimize_detection(PROV, other, PlaneWave(MED), Z, C_FB, GRID)
    assert r_fixed.best[1].eta_q == pytest.approx(r_direct.best[1].eta_q, rel=1e-12)


# --------------------------------------------------------------------------- #
# Sensitivity report
# --------------------------------------------------------------------------- #

def test_sensitivity_populated_on_result():
    r = optimize_detection(PROV, SPH, PlaneWave(MED), Z, C_FB, GRID)
    assert set(r.sensitivity.keys()) == {"radius", "NA", "waist", "m"}
    for k, v in r.sensitivity.items():
        assert np.isfinite(v), (k, v)
    # NA sensitivity about the (large-NA) best geometry is a real, nonzero slope.
    assert abs(r.sensitivity["NA"]) > 1e-3
    # Plane-wave waist sensitivity is exactly zero (waist-invariant).
    assert r.sensitivity["waist"] == 0.0


# --------------------------------------------------------------------------- #
# Multi-channel search: the two-beam win
# --------------------------------------------------------------------------- #

AXES = (Z.copy(), X.copy())   # beam A along lab z, beam B along lab x
C_MULTI = Constraints(na_max=0.8, directions_allowed=("forward", "backward"),
                      schemes_allowed=("optimal",), beam_axes_lab=AXES, max_channels=2)


def test_multichannel_returns_channel_set_and_result():
    n = np.array([1.0, 1.0, 1.0])
    r = optimize_detection(PROV, SPH, PlaneWave(MED), n, C_MULTI, GRID)
    assert r.best_channel_set is not None
    assert r.best_multichannel is not None
    assert isinstance(r.best_multichannel, MultiChannelResult)
    # One collection port per distinct beam axis.
    axes_used = {tuple(np.round(ch.propagation_lab, 6)) for ch in r.best_channel_set}
    assert len(axes_used) == len(r.best_channel_set)


def _single_beam_collected(n, axis):
    """Best collected Fisher info for a one-beam apparatus on `axis` (relative units)."""
    C1 = Constraints(na_max=0.8, directions_allowed=("forward", "backward"),
                     schemes_allowed=("optimal",), beam_axes_lab=(np.asarray(axis, float),),
                     max_channels=1)
    cand = _channel_candidates(PROV, SPH, PlaneWave(MED), n, C1, GRID)
    _set, _eta, collected = _best_channel_set(cand, 1)
    return collected


def test_two_beam_set_covers_arbitrary_n_hat_better_than_single_beam():
    # THE POINT (PHYSICS §4.5, MASTER_PLAN F3): for an arbitrary n_hat that is neither
    # purely axial nor purely transverse to either beam, the best TWO-beam set collects
    # strictly more information than the best SINGLE beam on either axis. A lab axis
    # poorly measured in one beam's forward lobe is well measured by the complementary
    # beam, and independent-channel Fisher information adds.
    n = np.array([1.0, 1.0, 1.0])       # equal components on all three lab axes
    r = optimize_detection(PROV, SPH, PlaneWave(MED), n, C_MULTI, GRID)
    two_beam = r.best_multichannel.fisher_total_rel

    best_single = max(_single_beam_collected(n, Z), _single_beam_collected(n, X))
    assert two_beam > best_single * (1.0 + 1e-6)
    # The optimizer actually places a channel on BOTH complementary axes (that is how it
    # covers the arbitrary direction).
    assert len(r.best_channel_set) == 2


def test_two_beam_matches_evaluate_channels_exactly():
    # The returned MultiChannelResult must be the canonical evaluate_channels output for
    # the chosen set (self-consistency of the recommendation).
    n = np.array([0.0, 1.0, 1.0])
    r = optimize_detection(PROV, SPH, PlaneWave(MED), n, C_MULTI, GRID)
    ref = evaluate_channels(PROV, SPH, r.best_channel_set, n, GRID)
    assert r.best_multichannel.fisher_total_rel == pytest.approx(ref.fisher_total_rel, rel=1e-12)
    assert r.best_multichannel.eta_q_total == pytest.approx(ref.eta_q_total, rel=1e-12)


def test_two_beam_collected_is_sum_of_independent_channels():
    # Independent channels' collected Fisher information adds (PHYSICS §4.5): the
    # two-beam collected info equals the sum of the two single-beam collected infos when
    # each single beam is optimized independently and both land on their best port.
    n = np.array([1.0, 1.0, 1.0])
    r = optimize_detection(PROV, SPH, PlaneWave(MED), n, C_MULTI, GRID)
    total = r.best_multichannel.fisher_total_rel
    parts = _single_beam_collected(n, Z) + _single_beam_collected(n, X)
    assert total == pytest.approx(parts, rel=1e-9)


def test_dependent_channels_not_produced():
    # The optimizer only ever recommends independent channels (it never sets
    # independent=False), so the chosen set is always summable.
    n = np.array([1.0, 0.0, 1.0])
    r = optimize_detection(PROV, SPH, PlaneWave(MED), n, C_MULTI, GRID)
    assert all(ch.independent for ch in r.best_channel_set)


def test_empty_candidate_set_raises():
    C_bad = Constraints(na_max=0.0, directions_allowed=("forward",), schemes_allowed=("optimal",))
    with pytest.raises(ValueError):
        optimize_detection(PROV, SPH, PlaneWave(MED), Z, C_bad, GRID)
