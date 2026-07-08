"""Vector spherical wave functions (VSWFs) and vector angular harmonics for the
far-field scattered-field expansion (PHYSICS.md §1.3, §2.1). Time convention
e^{-iωt}; Im(m) >= 0 for absorption (Bohren & Huffman). SI units; angles in radians.

Scope of this module
--------------------
The scattered far field of a sphere illuminated by an arbitrary beam is a sum over
partial waves (n, m) of Mie coefficients times beam-shape coefficients (BSCs) times
the transverse *vector angular harmonics* of the outgoing VSWFs. For a plane wave
the BSCs collapse to the m = +-1 weights that reproduce the classic S1/S2 amplitudes
(PHYSICS.md §1.3). This module delivers:

- ``pi_tau_m``   : generalized angular functions pi_{nm}, tau_{nm} for a fixed |m|,
                   vectorized over the angular grid, stable to n_max ~ 267.
- ``far_field``  : assemble spherical far-field components E_theta, E_phi on an
                   AngularGrid from per-(n, m) scattered TM/TE amplitudes.
- ``plane_wave_bsc`` / ``plane_wave_far_field`` : the plane-wave BSC weights and the
                   far field they produce. This is the correctness anchor — it
                   reproduces ``mie.plane_wave.scattering_amplitudes`` (S1/S2) to
                   ~1e-15, and is the seam the future GLMTProvider (T1.7) plugs into.

Conventions (fixed here; documented so sibling BSC code packs identically)
--------------------------------------------------------------------------
Generalized angular functions, defined for degree n >= 1 and order m with |m| <= n::

    pi_{nm}(theta)  =  m * P_n^{|m|}(cos theta) / sin theta
    tau_{nm}(theta) =  d/dtheta  P_n^{|m|}(cos theta)

with P_n^{|m|} the associated Legendre function built WITHOUT the Condon-Shortley
phase and WITHOUT the (n-|m|)!/(n+|m|)! renormalization, so that at |m| = 1 they
equal Bohren & Huffman's mie angle functions exactly:
``pi_{n,+1} = pi_n``, ``tau_{n,+1} = tau_n`` (mie.special.pi_tau). Consequences of
the explicit ``m`` factor and |m|-order Legendre: pi is odd in the sign of m
(``pi_{n,-m} = -pi_{nm}``) and tau is even (``tau_{n,-m} = tau_{nm}``).

Far field (outgoing, e^{-iωt}); the transverse spherical components are::

    E_theta(theta, phi) = sum_{n,m} [ A_{nm} tau_{nm} + B_{nm} pi_{nm} ] e^{i m phi}
    E_phi(theta, phi)   = sum_{n,m} i [ A_{nm} pi_{nm} + B_{nm} tau_{nm} ] e^{i m phi}

where ``A_{nm}`` is the scattered TM partial-wave amplitude and ``B_{nm}`` the
scattered TE amplitude. In GLMT (PHYSICS.md §2.1) ``A_{nm} = a_n * g_{nm}^{TM}`` and
``B_{nm} = b_n * g_{nm}^{TE}``, i.e. the SAME Mie ``a_n, b_n`` (mie.plane_wave)
multiply the beam-specific BSCs. Field amplitudes carry the arbitrary common far-
field constant of ``mie.plane_wave.scattering_amplitudes`` (VectorField units are
"arbitrary common units", INTERFACES.md §1); only ratios/shapes are physical.

Packing of the coefficient arrays
---------------------------------
Per-(n, m) coefficient arrays are indexed ``[n_idx, m_idx]`` with ``n_idx = n - 1``
(n = 1..n_max) and ``m_idx = m + n_max`` (m = -n_max..n_max), i.e. a dense
(n_max, 2*n_max + 1) array with the unused |m| > n corners set to zero. Helpers
``m_index`` / ``m_values`` document the m axis. This is the packing the GLMT BSC
module (INTERFACES.md §3, glmt/bsc.py) mirrors.
"""
from __future__ import annotations

import numpy as np

from ..types import AngularGrid, VectorField
from .plane_wave import mie_coefficients
from .special import nmax_wiscombe

__all__ = [
    "pi_tau_m",
    "m_values",
    "m_index",
    "far_field",
    "plane_wave_bsc",
    "plane_wave_far_field",
]


def m_values(n_max: int) -> np.ndarray:
    """The azimuthal orders m = -n_max .. n_max carried by the packed [n, m] axis."""
    return np.arange(-n_max, n_max + 1)


def m_index(m: int, n_max: int) -> int:
    """Column index of order ``m`` in a packed (n_max, 2*n_max+1) coefficient array."""
    return m + n_max


def pi_tau_m(n_max: int, m: int, theta: np.ndarray
             ) -> tuple[np.ndarray, np.ndarray]:
    """Generalized angular functions pi_{nm}, tau_{nm} for a fixed order ``m``.

    Returns arrays of shape ``(n_max,) + theta.shape`` (index 0 == n = 1). Vectorized
    over ``theta``; stable to n_max ~ 267 (PHYSICS.md §7). Rows with n < |m| are zero
    (pi_{nm}, tau_{nm} are undefined/zero there).

    Definitions (module docstring): ``pi_{nm} = m P_n^{|m|}(cos t)/sin t``,
    ``tau_{nm} = d P_n^{|m|}(cos t)/dtheta``, no Condon-Shortley phase, so that
    ``pi_tau_m(n_max, 1, theta)`` reproduces ``mie.special.pi_tau`` exactly.

    The 1/sin(theta) factor is finite on a Gauss-Legendre cos-theta grid (nodes never
    land on theta = 0 or pi); callers evaluating exactly at the poles must pass nodes
    that avoid them (the far-field AngularGrid does — types.py).
    """
    am = abs(int(m))
    if am > n_max:
        shp = (n_max,) + np.shape(theta)
        return np.zeros(shp), np.zeros(shp)

    theta = np.asarray(theta, dtype=float)
    ct = np.cos(theta)
    st = np.sin(theta)

    shp = (n_max + 1,) + theta.shape
    # P[n] holds P_n^{|m|}(cos theta); rows below n = am stay zero.
    P = np.zeros(shp)

    # Seed P_{am}^{am} = (2*am - 1)!! * sin(theta)^am  (no Condon-Shortley phase).
    double_factorial = 1.0
    for k in range(1, 2 * am, 2):
        double_factorial *= k
    P[am] = double_factorial * st ** am
    if am + 1 <= n_max:
        P[am + 1] = ct * (2 * am + 1) * P[am]
    # Upward recurrence in n at fixed |m|:
    #   (n - am) P_n = (2n - 1) cos(t) P_{n-1} - (n - 1 + am) P_{n-2}
    for n in range(am + 2, n_max + 1):
        P[n] = ((2 * n - 1) * ct * P[n - 1] - (n - 1 + am) * P[n - 2]) / (n - am)

    # tau_{nm} = dP_n^{|m|}/dtheta = ( n cos(t) P_n - (n + am) P_{n-1} ) / sin(t).
    n_arr = np.arange(n_max + 1).reshape((-1,) + (1,) * theta.ndim)
    P_shift = np.zeros_like(P)
    P_shift[1:] = P[:-1]  # P_{n-1}
    dP = (n_arr * ct * P - (n_arr + am) * P_shift) / st

    pi_nm = int(m) * P[1:] / st  # explicit m factor -> odd in sign(m)
    tau_nm = dP[1:]              # even in sign(m)
    return pi_nm, tau_nm


def far_field(a_tm: np.ndarray, b_te: np.ndarray, grid: AngularGrid) -> VectorField:
    """Assemble the far-field VectorField from per-(n, m) scattered amplitudes.

    Parameters
    ----------
    a_tm, b_te : (n_max, 2*n_max + 1) complex
        Scattered TM amplitudes ``A_{nm}`` and TE amplitudes ``B_{nm}``, packed
        ``[n - 1, m + n_max]`` (see module docstring). In GLMT terms
        ``A_{nm} = a_n g_{nm}^{TM}``, ``B_{nm} = b_n g_{nm}^{TE}``.
    grid : AngularGrid
        Far-field evaluation grid (theta, phi meshed 'ij').

    Returns
    -------
    VectorField with E_theta, E_phi of shape (Ntheta, Nphi), in the same arbitrary
    common far-field units as mie.plane_wave.scattering_amplitudes.
    """
    a_tm = np.asarray(a_tm)
    b_te = np.asarray(b_te)
    n_max = a_tm.shape[0]
    if a_tm.shape != (n_max, 2 * n_max + 1) or b_te.shape != a_tm.shape:
        raise ValueError(
            f"a_tm/b_te must be (n_max, 2*n_max+1); got {a_tm.shape}, {b_te.shape}"
        )

    theta = np.asarray(grid.theta, dtype=float)   # (Ntheta,)
    phi = np.asarray(grid.phi, dtype=float)        # (Nphi,)
    n_theta = theta.shape[0]
    n_phi = phi.shape[0]

    e_theta = np.zeros((n_theta, n_phi), dtype=complex)
    e_phi = np.zeros((n_theta, n_phi), dtype=complex)

    ms = m_values(n_max)
    for m in ms:
        col = m_index(int(m), n_max)
        A_col = a_tm[:, col]  # (n_max,)  A_{nm}
        B_col = b_te[:, col]  # (n_max,)  B_{nm}
        if not (np.any(A_col) or np.any(B_col)):
            continue
        pi_nm, tau_nm = pi_tau_m(n_max, int(m), theta)  # (n_max, Ntheta)
        # Sum over n -> (Ntheta,)
        theta_sum = np.tensordot(A_col, tau_nm, axes=(0, 0)) \
            + np.tensordot(B_col, pi_nm, axes=(0, 0))
        phi_sum = 1j * (np.tensordot(A_col, pi_nm, axes=(0, 0))
                        + np.tensordot(B_col, tau_nm, axes=(0, 0)))
        emphi = np.exp(1j * int(m) * phi)  # (Nphi,)
        e_theta += theta_sum[:, None] * emphi[None, :]
        e_phi += phi_sum[:, None] * emphi[None, :]

    return VectorField(grid=grid, E_theta=e_theta, E_phi=e_phi)


def plane_wave_bsc(n_max: int) -> tuple[np.ndarray, np.ndarray]:
    """Plane-wave beam-shape weights g_{nm}^{TM}, g_{nm}^{TE} (x-polarized, +z).

    Packed (n_max, 2*n_max + 1) arrays. Only m = +-1 are nonzero (PHYSICS.md §1.3,
    §2.1): with ``C_n = (2n + 1) / (n(n + 1))``,

        g_{n,+1}^{TM} = g_{n,-1}^{TM} =  C_n / 2      (even in m)
        g_{n,+1}^{TE} = -g_{n,-1}^{TE} = C_n / 2      (odd in m: g^{TE}_{nm} = m C_n/2)

    Multiplying these by the Mie ``a_n, b_n`` and feeding to ``far_field`` reproduces
    ``E_theta = S2 cos(phi)``, ``E_phi = -S1 sin(phi)`` (PHYSICS.md §1.3). This is the
    plane-wave limit every BSC method must reproduce (VALIDATION.md §3, G-LIMIT).
    """
    n = np.arange(1, n_max + 1)
    C = (2 * n + 1) / (n * (n + 1))
    g_tm = np.zeros((n_max, 2 * n_max + 1), dtype=complex)
    g_te = np.zeros((n_max, 2 * n_max + 1), dtype=complex)
    ip = m_index(1, n_max)
    im = m_index(-1, n_max)
    g_tm[:, ip] = C / 2
    g_tm[:, im] = C / 2
    g_te[:, ip] = C / 2
    g_te[:, im] = -C / 2
    return g_tm, g_te


def plane_wave_far_field(m: complex, x: float, grid: AngularGrid,
                         n_max: int | None = None) -> VectorField:
    """Far field of an x-polarized plane wave via the VSWF reconstruction.

    Convenience wrapper: computes the Mie ``a_n, b_n`` (mie.plane_wave) at relative
    index ``m`` and size parameter ``x``, applies the plane-wave BSCs, and assembles
    the far field on ``grid``. Reproduces ``mie.plane_wave.scattering_amplitudes``
    S1/S2 to ~1e-15 (the correctness anchor and the seam GLMTProvider consumes).

    ``m`` is the relative refractive index (complex allowed, Im(m) >= 0); ``x = k a``
    the size parameter (dimensionless). ``n_max`` defaults to Wiscombe(x).
    """
    if n_max is None:
        n_max = nmax_wiscombe(x)
    a_n, b_n = mie_coefficients(m, x, n_max)  # (n_max,) each, index 0 == n = 1
    g_tm, g_te = plane_wave_bsc(n_max)
    a_tm = a_n[:, None] * g_tm  # A_{nm} = a_n g_{nm}^{TM}
    b_te = b_n[:, None] * g_te  # B_{nm} = b_n g_{nm}^{TE}
    return far_field(a_tm, b_te, grid)
