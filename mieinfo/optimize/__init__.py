"""Detection-geometry optimizer — the project's recommendation (INTERFACES.md §7).

Given a scatterer, a probe beam, and a sensed displacement direction ``n_hat``, rank the
achievable collection geometries by the single figure of merit ``eta_q`` (collected /
total Fisher information, PHYSICS.md §4.2-4.3) and, for a multi-beam apparatus, pick the
SET of detection channels that jointly cover ``n_hat`` best (PHYSICS.md §4.5). The core
result reproduced here (PHYSICS.md §8, MASTER_PLAN F2): axial info is backward-weighted,
transverse info is forward-weighted, and a two-beam layout covers every lab axis.

The inner loop caches the field and only re-masks it per collection geometry (PHYSICS.md
§7) — no field is re-solved when only the cone changes. Time convention e^{-iωt}; SI
units unless a name says otherwise.
"""
from __future__ import annotations

from .objective import (
    GeometryCache,
    build_geometry_cache,
    evaluate_geometry,
    score_geometry,
    sensitivity,
)
from .search import Constraints, OptResult, optimize_detection

__all__ = [
    "Constraints",
    "OptResult",
    "optimize_detection",
    "GeometryCache",
    "build_geometry_cache",
    "evaluate_geometry",
    "score_geometry",
    "sensitivity",
]
