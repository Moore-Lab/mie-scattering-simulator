"""Detection metrics: imprecision PSD, measurement backaction, and distance to the
standard quantum limit (PHYSICS.md §4.2-4.3, §5). Time convention e^{-iωt}.

The physics (PHYSICS.md §4.2-4.3):

  * Imprecision noise PSD  S_q^{imp}(ω) ∝ 1 / F_q(Ω).  Since a collection cone with a
    (possibly imperfect) LO collects F_q(Ω) = η_q · F_q(4π) = η_q · f_total, the
    imprecision *relative to* the total-information floor 1/f_total is

        s_imp_rel = S_q^{imp} / (1 / f_total) = f_total / F_q(Ω) = 1 / η_q .

    η_q here already folds in the LO mode-overlap κ (schemes.py); s_imp_rel ≥ 1, and
    equals 1 only for a full-4π optimal-LO measurement.

  * Backaction (photon-recoil heating) rate  Γ_ba ∝ P_scattered ∝ Q_sca, set by the
    *total* scattered power over all 4π and INDEPENDENT of what is collected. In the
    relative units the package reports, γ_ba_rel = Q_sca directly.

  * Standard quantum limit. Γ_ba is fixed by Q_sca; the imprecision–backaction product
    S_q^{imp}·Γ_ba is Heisenberg-bounded and the optimal scheme saturates it. Because
    Γ_ba does not depend on collection and S_q^{imp} ∝ 1/(η_q f_total), the product is
    minimised by η_q → 1. The distance to the SQL is thus the factor by which the
    product exceeds its minimum,

        sql_distance = 1 / η_q  (≥ 1; == 1 exactly at the SQL, η_q → 1).
"""
from __future__ import annotations

import numpy as np

# η below this is treated as "no information collected": metrics saturate rather than
# divide by ~0 (physically an infinitely imprecise, infinitely-far-from-SQL detector).
_ETA_FLOOR = 1e-300


def imprecision_rel(eta_q: float) -> float:
    """Relative imprecision PSD s_imp_rel = 1/η_q, in units of 1/f_total (≥ 1).

    η_q is the collected information fraction *including* the LO mode overlap κ.
    Returns +inf for η_q ≤ 0 (nothing collected → unbounded imprecision).
    """
    e = float(eta_q)
    if e <= _ETA_FLOOR:
        return float("inf")
    return 1.0 / e


def backaction_rel(q_sca: float) -> float:
    """Relative backaction rate γ_ba_rel ∝ Q_sca (PHYSICS.md §4.3).

    Backaction is set by the total 4π scattered power and is independent of the
    collection geometry, so it is simply Q_sca in the package's relative units.
    """
    return float(q_sca)


def sql_distance(eta_q: float) -> float:
    """Distance to the standard quantum limit, sql_distance = 1/η_q (≥ 1).

    1.0 means the configuration sits exactly at the SQL (η_q → 1, optimal LO over the
    full sphere); larger values are the factor by which the imprecision–backaction
    product exceeds the Heisenberg bound. Returns +inf for η_q ≤ 0.
    """
    e = float(eta_q)
    if e <= _ETA_FLOOR:
        return float("inf")
    return 1.0 / e
