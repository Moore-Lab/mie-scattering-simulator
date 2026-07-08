"""Sensed-parameter direction and the InformationPattern (INTERFACES.md §5).

The sensed parameter is displacement along an arbitrary unit vector n̂:
    ∂E_s/∂(r_s·n̂) = Σ_j n̂_j ∂E_s/∂r_{s,j}     (PHYSICS.md §4.4)
— a linear combination of the three Cartesian derivatives, so an arbitrary direction
costs no new field solves. Time convention e^{-iωt}.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..glmt.scatter import FieldProvider, field_derivative
from ..types import AngularGrid, FieldDerivative, Sphere
from .fisher import fisher_total, info_density


def combine_direction(deriv: FieldDerivative, n_hat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """dE/dq for q = r_s·n̂ as the linear combination of the Cartesian derivatives.
    Returns (dE_theta_q, dE_phi_q). n̂ is normalized to a unit vector."""
    n = np.asarray(n_hat, dtype=float)
    n = n / np.linalg.norm(n)
    dE_theta_q = np.tensordot(n, deriv.dE_theta, axes=(0, 0))  # (Ntheta, Nphi)
    dE_phi_q = np.tensordot(n, deriv.dE_phi, axes=(0, 0))
    return dE_theta_q, dE_phi_q


@dataclass(frozen=True)
class InformationPattern:
    grid: AngularGrid
    density: np.ndarray        # dF_q/dΩ
    n_hat: np.ndarray
    f_total: float


def information_pattern(provider: FieldProvider, grid: AngularGrid, sphere: Sphere,
                        beam, r_s_m: np.ndarray, n_hat: np.ndarray,
                        method: str = "analytic", n_max: int | None = None) -> InformationPattern:
    """Build the InformationPattern for displacement along n̂: field_derivative ->
    combine_direction -> info_density. Works for any FieldProvider (PHYSICS.md §4.4)."""
    deriv = field_derivative(provider, grid, sphere, beam, r_s_m, method=method, n_max=n_max)
    deriv_q = combine_direction(deriv, n_hat)
    density = info_density(deriv_q, grid)
    n = np.asarray(n_hat, dtype=float)
    return InformationPattern(grid=grid, density=density, n_hat=n / np.linalg.norm(n),
                              f_total=fisher_total(density, grid))
