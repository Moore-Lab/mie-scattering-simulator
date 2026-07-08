"""Fisher-information density and total (INTERFACES.md §5, PHYSICS.md §4.1).

For an ideal shot-noise-limited coherent measurement that mode-matches a local
oscillator to ∂E_s/∂q, the classical Fisher information rate is
    F_q(Ω) ∝ ∫_Ω |∂E_s/∂q|² dΩ ,   density dF_q/dΩ ∝ |∂E_s/∂q|² .
Units are arbitrary but common across the field and its derivative (ratios — η, peak
angles, forward fractions — are what the package delivers). Time convention e^{-iωt}.
"""
from __future__ import annotations

import numpy as np

from ..types import AngularGrid


def info_density(deriv_q: tuple[np.ndarray, np.ndarray], grid: AngularGrid) -> np.ndarray:
    """(dE_theta_q, dE_phi_q) -> dF_q/dΩ = |dE/dq|² on the grid."""
    dE_theta_q, dE_phi_q = deriv_q
    return np.abs(dE_theta_q) ** 2 + np.abs(dE_phi_q) ** 2


def fisher_total(density: np.ndarray, grid: AngularGrid) -> float:
    """Integral of the information density over the full 4π solid angle."""
    return float(np.sum(density * grid.w_solid))
