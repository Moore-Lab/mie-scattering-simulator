"""Core value types (orchestrator-owned, frozen at M0 — INTERFACES.md §1).

Time convention e^{-iωt}. SI units unless a field name says otherwise. These are the
shared vocabulary every track codes against; changing a signature is a CONTRACT-CHANGE
(CONVENTIONS.md §5).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Sphere:
    radius_m: float
    m: complex  # relative refractive index n_particle / n_medium (Im(m) >= 0 absorbs)


@dataclass(frozen=True)
class Medium:
    n: float = 1.0                        # vacuum around the sphere
    wavelength_vacuum_m: float = 532e-9   # DETECTION default; the trap Medium uses 1064e-9

    @property
    def k(self) -> float:
        """Wavenumber in the medium, 2*pi*n / lambda_vacuum."""
        return 2.0 * np.pi * self.n / self.wavelength_vacuum_m


@dataclass(frozen=True)
class AngularGrid:
    """Gauss-Legendre in cos(theta), uniform in phi; carries solid-angle weights.

    theta = arccos(leggauss nodes) so theta runs high->low with the node order; every
    consumer meshes theta/phi as (Ntheta, Nphi) with indexing='ij'. w_solid sums to ~4*pi
    over the full sphere (dOmega = dphi * w_mu, w_mu already integrating dcos(theta)).
    """
    theta: np.ndarray     # (Ntheta,)
    phi: np.ndarray       # (Nphi,)
    w_solid: np.ndarray   # (Ntheta, Nphi)

    @classmethod
    def full_sphere(cls, ntheta: int, nphi: int) -> "AngularGrid":
        mu, w_mu = np.polynomial.legendre.leggauss(ntheta)
        theta = np.arccos(mu)
        phi = np.linspace(0.0, 2.0 * np.pi, nphi, endpoint=False)
        dphi = 2.0 * np.pi / nphi
        w_solid = np.outer(w_mu, np.full(nphi, dphi))
        return cls(theta=theta, phi=phi, w_solid=w_solid)


@dataclass(frozen=True)
class VectorField:
    """Far-field scattered E on an AngularGrid, spherical components (arbitrary units)."""
    grid: AngularGrid
    E_theta: np.ndarray   # (Ntheta, Nphi) complex
    E_phi: np.ndarray     # (Ntheta, Nphi) complex

    def intensity(self) -> np.ndarray:
        return np.abs(self.E_theta) ** 2 + np.abs(self.E_phi) ** 2


@dataclass(frozen=True)
class FieldDerivative:
    """dE_s/dr_j for j in (x, y, z) on the same grid (axis 0 = x,y,z)."""
    grid: AngularGrid
    dE_theta: np.ndarray  # (3, Ntheta, Nphi) complex
    dE_phi: np.ndarray    # (3, Ntheta, Nphi) complex
