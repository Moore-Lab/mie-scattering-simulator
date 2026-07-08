"""Special functions for plane-wave Mie (PHYSICS.md §1.1, §1.3). Time convention
e^{-iωt}; Im(m) >= 0 for absorption (Bohren & Huffman).

Riccati-Bessel:  psi_n(z) = z j_n(z),  chi_n(z) = -z y_n(z),  xi_n = psi_n - i chi_n.
D_n(z) = psi_n'(z)/psi_n(z) by downward recurrence. Recurrence directions are load-
bearing (downward for D_n, upward for psi/chi); do not change them.
"""
from __future__ import annotations

import numpy as np


def nmax_wiscombe(x: float) -> int:
    """Wiscombe (1980) series truncation with the prototype's +3 safety margin."""
    return int(np.ceil(x + 4.05 * x ** (1.0 / 3.0) + 2)) + 3


def logarithmic_derivative(m: complex, x: float, nmax: int) -> np.ndarray:
    """D_n(m*x) for n = 1..nmax by downward recurrence (numerically stable).

    The burn-in margin above |mx| MUST scale with |mx| for the transient to decay:
    a fixed margin (e.g. +16) under-resolves D_n at large |mx| — at x=236, m=1.45
    (|mx|=342) it corrupts D_1 by ~7e-3 and every a_n by ~5e-4, propagating to
    Q_sca/Q_back and thus Gamma_ba. Textbook BHMIE rule: 15 + 15*sqrt(|mx|)
    (PHYSICS.md §1.1; regression-guarded by prototype/validate.py).
    """
    mx = m * x
    nstart = max(nmax, int(np.ceil(abs(mx)))) + 15 + int(np.ceil(15.0 * abs(mx) ** 0.5))
    D = np.zeros(nstart + 1, dtype=complex)
    for n in range(nstart, 0, -1):
        D[n - 1] = n / mx - 1.0 / (D[n] + n / mx)
    return D[1:nmax + 1]  # D_1 .. D_nmax


def riccati_psi_chi(x: float, nmax: int) -> tuple[np.ndarray, np.ndarray]:
    """psi_n(x), chi_n(x) for n = 0..nmax by upward recurrence (stable for real x)."""
    psi = np.zeros(nmax + 1)
    chi = np.zeros(nmax + 1)
    psi_prev = np.cos(x)    # psi_{-1}
    chi_prev = -np.sin(x)   # chi_{-1}
    psi[0] = np.sin(x)      # psi_0
    chi[0] = np.cos(x)      # chi_0
    for n in range(1, nmax + 1):
        psi[n] = (2 * n - 1) / x * psi[n - 1] - psi_prev
        chi[n] = (2 * n - 1) / x * chi[n - 1] - chi_prev
        psi_prev = psi[n - 1]
        chi_prev = chi[n - 1]
    return psi, chi


def pi_tau(theta: np.ndarray, nmax: int) -> tuple[np.ndarray, np.ndarray]:
    """Angle functions pi_n(cos theta), tau_n(cos theta) for n = 1..nmax.

    Returns arrays of shape (nmax,) + theta.shape (index 0 == n=1).
    """
    mu = np.cos(theta)
    pi = np.zeros((nmax + 1,) + np.shape(mu))
    tau = np.zeros((nmax + 1,) + np.shape(mu))
    pi[1] = np.ones_like(mu)
    tau[1] = mu * pi[1]
    for n in range(2, nmax + 1):
        pi[n] = ((2 * n - 1) / (n - 1)) * mu * pi[n - 1] - (n / (n - 1)) * pi[n - 2]
        tau[n] = n * mu * pi[n] - (n + 1) * pi[n - 1]
    return pi[1:], tau[1:]
