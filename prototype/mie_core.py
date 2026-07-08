"""
Reference plane-wave Mie scattering (Bohren & Huffman 1983 algorithm).

Purpose: a from-scratch, numerically-stable implementation used to (a) validate
against miepython and (b) generate golden reference values for the Claude Code
project. Notation follows Bohren & Huffman, "Absorption and Scattering of Light
by Small Particles" (1983), Ch. 4.

Conventions
-----------
x  = size parameter = 2*pi*a/lambda_medium = k*a
m  = relative refractive index = n_particle / n_medium (complex allowed)
psi_n(z) = z * j_n(z)              (Riccati-Bessel, 1st kind)
chi_n(z) = -z * y_n(z)             (Riccati-Bessel, 2nd kind)
xi_n(z)  = psi_n(z) - i*chi_n(z) = z * h_n^(1)(z)
D_n(z)   = psi_n'(z)/psi_n(z)      (logarithmic derivative, downward recurrence)
"""

import numpy as np


def nmax_wiscombe(x):
    """Wiscombe (1980) series-truncation criterion, with a small safety margin."""
    return int(np.ceil(x + 4.05 * x ** (1.0 / 3.0) + 2)) + 3


def logarithmic_derivative(m, x, nmax):
    """D_n(m*x) for n=1..nmax via downward recurrence (numerically stable)."""
    mx = m * x
    # Start the downward recurrence well above nmax for accuracy. The burn-in
    # margin above |mx| must scale with |mx| for the transient to decay: a fixed
    # +16 under-resolves D_n at large |mx| (e.g. x=236, m=1.45 -> |mx|=342 gives
    # ~7e-3 error in D_1 and ~5e-4 in every a_n). Textbook BHMIE rule 15+15*sqrt.
    nstart = max(nmax, int(np.ceil(abs(mx)))) + 15 + int(np.ceil(15.0 * abs(mx) ** 0.5))
    D = np.zeros(nstart + 1, dtype=complex)
    for n in range(nstart, 0, -1):
        D[n - 1] = n / mx - 1.0 / (D[n] + n / mx)
    return D[1:nmax + 1]  # D_1 .. D_nmax


def mie_coefficients(m, x, nmax=None):
    """Return Mie coefficients a_n, b_n for n=1..nmax."""
    if nmax is None:
        nmax = nmax_wiscombe(x)
    Dn = logarithmic_derivative(m, x, nmax)

    # Riccati-Bessel psi_n(x), chi_n(x) by upward recurrence (stable for real x).
    psi = np.zeros(nmax + 1)
    chi = np.zeros(nmax + 1)
    psi_prev = np.cos(x)   # psi_{-1}
    chi_prev = -np.sin(x)  # chi_{-1}
    psi[0] = np.sin(x)     # psi_0
    chi[0] = np.cos(x)     # chi_0
    for n in range(1, nmax + 1):
        psi[n] = (2 * n - 1) / x * psi[n - 1] - psi_prev
        chi[n] = (2 * n - 1) / x * chi[n - 1] - chi_prev
        psi_prev = psi[n - 1]
        chi_prev = chi[n - 1]

    xi = psi - 1j * chi  # xi_n(x)

    a = np.zeros(nmax, dtype=complex)
    b = np.zeros(nmax, dtype=complex)
    for n in range(1, nmax + 1):
        Dn_n = Dn[n - 1]
        fa = Dn_n / m + n / x
        fb = m * Dn_n + n / x
        a[n - 1] = (fa * psi[n] - psi[n - 1]) / (fa * xi[n] - xi[n - 1])
        b[n - 1] = (fb * psi[n] - psi[n - 1]) / (fb * xi[n] - xi[n - 1])
    return a, b


def efficiencies(m, x, nmax=None):
    """Q_ext, Q_sca, Q_back, and asymmetry parameter g = <cos theta>."""
    a, b = mie_coefficients(m, x, nmax)
    n = np.arange(1, len(a) + 1)
    prefac = 2.0 / x ** 2
    q_ext = prefac * np.sum((2 * n + 1) * np.real(a + b))
    q_sca = prefac * np.sum((2 * n + 1) * (np.abs(a) ** 2 + np.abs(b) ** 2))
    # Asymmetry parameter (B&H 4.62). First sum runs n=1..nmax-1.
    nn = n[:-1]
    g_terms = (nn * (nn + 2) / (nn + 1)) * np.real(
        a[:-1] * np.conj(a[1:]) + b[:-1] * np.conj(b[1:])
    )
    g_terms2 = ((2 * n + 1) / (n * (n + 1))) * np.real(a * np.conj(b))
    g = (4.0 / (x ** 2 * q_sca)) * (np.sum(g_terms) + np.sum(g_terms2))
    # Backscatter
    q_back = (1.0 / x ** 2) * np.abs(np.sum((2 * n + 1) * (-1) ** n * (a - b))) ** 2
    return q_ext, q_sca, q_back, g


def pi_tau(theta, nmax):
    """Angle functions pi_n(cos theta), tau_n(cos theta) for n=1..nmax."""
    mu = np.cos(theta)
    pi = np.zeros((nmax + 1,) + np.shape(mu))
    tau = np.zeros((nmax + 1,) + np.shape(mu))
    pi[1] = np.ones_like(mu)
    tau[1] = mu * pi[1]
    for n in range(2, nmax + 1):
        pi[n] = ((2 * n - 1) / (n - 1)) * mu * pi[n - 1] - (n / (n - 1)) * pi[n - 2]
        tau[n] = n * mu * pi[n] - (n + 1) * pi[n - 1]
    return pi[1:], tau[1:]


def scattering_amplitudes(m, x, theta, nmax=None):
    """
    Complex scattering amplitudes S1(theta), S2(theta).
    S1 governs perpendicular polarization, S2 parallel (B&H eq. 4.74).
    theta may be an array.
    """
    if nmax is None:
        nmax = nmax_wiscombe(x)
    a, b = mie_coefficients(m, x, nmax)
    pi, tau = pi_tau(np.atleast_1d(theta), nmax)
    n = np.arange(1, nmax + 1)
    fac = (2 * n + 1) / (n * (n + 1))
    # Broadcast: (nmax, Ntheta)
    S1 = np.tensordot(fac * a, pi, axes=(0, 0)) + np.tensordot(fac * b, tau, axes=(0, 0))
    S2 = np.tensordot(fac * a, tau, axes=(0, 0)) + np.tensordot(fac * b, pi, axes=(0, 0))
    return S1, S2


if __name__ == "__main__":
    # quick smoke test
    qext, qsca, qback, g = efficiencies(1.5 + 0j, 10.0)
    print(f"x=10, m=1.5:  Qext={qext:.6f} Qsca={qsca:.6f} Qback={qback:.6f} g={g:.6f}")
