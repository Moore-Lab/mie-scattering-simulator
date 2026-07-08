"""Plane-wave Mie coefficients, efficiencies, and far-field amplitudes
(Bohren & Huffman 1983, Ch. 4; PHYSICS.md §1). Time convention e^{-iωt}.

Numerics promoted verbatim from the validated prototype (prototype/mie_core.py) —
agrees with miepython to ~1e-12 across x = 0.3–236 (G-GOLD, VALIDATION.md §2). The
far-field amplitudes carry a global x2 constant vs miepython that cancels in every
delivered (ratio) quantity; the physical check is the normalized angular shape.
"""
from __future__ import annotations

import numpy as np

from ..types import Medium, Sphere
from .special import logarithmic_derivative, nmax_wiscombe, pi_tau, riccati_psi_chi


def size_parameter(sphere: Sphere, medium: Medium) -> float:
    """x = k a = 2*pi*n_medium*a / lambda_vacuum (detection light)."""
    return medium.k * sphere.radius_m


def mie_coefficients(m: complex, x: float, nmax: int | None = None
                     ) -> tuple[np.ndarray, np.ndarray]:
    """Mie coefficients a_n, b_n for n = 1..nmax (index 0 == n=1)."""
    if nmax is None:
        nmax = nmax_wiscombe(x)
    Dn = logarithmic_derivative(m, x, nmax)
    psi, chi = riccati_psi_chi(x, nmax)
    xi = psi - 1j * chi
    a = np.zeros(nmax, dtype=complex)
    b = np.zeros(nmax, dtype=complex)
    for n in range(1, nmax + 1):
        Dn_n = Dn[n - 1]
        fa = Dn_n / m + n / x
        fb = m * Dn_n + n / x
        a[n - 1] = (fa * psi[n] - psi[n - 1]) / (fa * xi[n] - xi[n - 1])
        b[n - 1] = (fb * psi[n] - psi[n - 1]) / (fb * xi[n] - xi[n - 1])
    return a, b


def efficiencies(m: complex, x: float, nmax: int | None = None
                 ) -> tuple[float, float, float, float]:
    """Q_ext, Q_sca, Q_back, and asymmetry g = <cos theta> (PHYSICS.md §1.2).

    Q_back is a coherent alternating sum: catastrophic cancellation floors its float64
    accuracy at ~1e-5–1e-6 for large x, while Q_ext/Q_sca/g stay ~1e-12 (PHYSICS.md §1.2).
    """
    a, b = mie_coefficients(m, x, nmax)
    n = np.arange(1, len(a) + 1)
    prefac = 2.0 / x ** 2
    q_ext = prefac * np.sum((2 * n + 1) * np.real(a + b))
    q_sca = prefac * np.sum((2 * n + 1) * (np.abs(a) ** 2 + np.abs(b) ** 2))
    nn = n[:-1]
    g_terms = (nn * (nn + 2) / (nn + 1)) * np.real(
        a[:-1] * np.conj(a[1:]) + b[:-1] * np.conj(b[1:])
    )
    g_terms2 = ((2 * n + 1) / (n * (n + 1))) * np.real(a * np.conj(b))
    g = (4.0 / (x ** 2 * q_sca)) * (np.sum(g_terms) + np.sum(g_terms2))
    q_back = (1.0 / x ** 2) * np.abs(np.sum((2 * n + 1) * (-1) ** n * (a - b))) ** 2
    return float(q_ext), float(q_sca), float(q_back), float(g)


def scattering_amplitudes(m: complex, x: float, theta: np.ndarray,
                          nmax: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Complex far-field amplitudes S1(theta), S2(theta) (B&H eq. 4.74).

    S1 governs perpendicular polarization, S2 parallel. `theta` may be an array; the
    return shape matches np.atleast_1d(theta).
    """
    if nmax is None:
        nmax = nmax_wiscombe(x)
    a, b = mie_coefficients(m, x, nmax)
    pi, tau = pi_tau(np.atleast_1d(theta), nmax)
    n = np.arange(1, nmax + 1)
    fac = (2 * n + 1) / (n * (n + 1))
    S1 = np.tensordot(fac * a, pi, axes=(0, 0)) + np.tensordot(fac * b, tau, axes=(0, 0))
    S2 = np.tensordot(fac * a, tau, axes=(0, 0)) + np.tensordot(fac * b, pi, axes=(0, 0))
    return S1, S2
