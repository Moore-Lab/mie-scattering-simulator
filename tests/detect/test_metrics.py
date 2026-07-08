"""Detection metrics: imprecision, backaction, SQL distance (PHYSICS.md §4.2-4.3)."""
from __future__ import annotations

import math

import numpy as np
import pytest

from mieinfo.detect.metrics import backaction_rel, imprecision_rel, sql_distance


def test_imprecision_is_inverse_eta():
    for eta in (0.1, 0.25, 0.5, 0.9, 1.0):
        assert imprecision_rel(eta) == pytest.approx(1.0 / eta)


def test_imprecision_floor_is_one_at_full_collection():
    # η=1 (full-4π optimal LO) is the information floor: s_imp_rel == 1.
    assert imprecision_rel(1.0) == pytest.approx(1.0)
    # η < 1 always raises imprecision above the floor.
    assert imprecision_rel(0.5) > 1.0


def test_sql_distance_is_inverse_eta_and_at_least_one():
    for eta in np.linspace(0.05, 1.0, 20):
        d = sql_distance(float(eta))
        assert d == pytest.approx(1.0 / eta)
        assert d >= 1.0 - 1e-12
    # Exactly at the SQL when η -> 1.
    assert sql_distance(1.0) == pytest.approx(1.0)


def test_backaction_tracks_q_sca_and_is_collection_independent():
    # γ_ba depends only on Q_sca (total scattered power), not on any geometry/η.
    assert backaction_rel(2.9041) == pytest.approx(2.9041)
    assert backaction_rel(0.0) == 0.0
    # Monotonic in Q_sca.
    assert backaction_rel(3.0) > backaction_rel(1.0)


def test_zero_eta_saturates_metrics_without_nan():
    # Nothing collected -> unbounded imprecision and infinite SQL distance, not NaN.
    assert math.isinf(imprecision_rel(0.0))
    assert math.isinf(sql_distance(0.0))
    assert not math.isnan(imprecision_rel(0.0))
    assert not math.isnan(sql_distance(0.0))


def test_imprecision_backaction_product_minimized_at_sql():
    # S_imp * Γ_ba ∝ (1/η)*Q_sca is minimized (per fixed Q_sca) as η -> 1, i.e. at the
    # SQL, mirroring sql_distance = 1/η.
    q_sca = 2.5
    prod = [imprecision_rel(e) * backaction_rel(q_sca) for e in (0.2, 0.5, 1.0)]
    assert prod[0] > prod[1] > prod[2]
    assert prod[2] == pytest.approx(q_sca)  # η=1 -> product == Q_sca (the SQL value)
