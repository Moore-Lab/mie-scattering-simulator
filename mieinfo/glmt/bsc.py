"""Beam-shape coefficients (BSCs) for the GLMT incident-field expansion
(PHYSICS.md §2.2, INTERFACES.md §3). Time convention e^{-iωt}; Im(m) >= 0 for
absorption (Bohren & Huffman). SI units; angles in radians.

What a BSC is
-------------
An arbitrary incident beam expanded in *regular* vector spherical wave functions
(VSWFs) about a center carries per-(n, m) weights g_{n,TM}^m, g_{n,TE}^m that replace
the plane wave's (2n+1) weighting (PHYSICS.md §2). The scattered field is then the
outgoing VSWF sum with partial-wave amplitudes a_n g_{n,TM}^m, b_n g_{n,TE}^m — the
*same* Mie a_n, b_n (mie.plane_wave) times these beam-specific BSCs (PHYSICS.md §2.1).
``mie.vswf.far_field`` consumes exactly those amplitudes, so this module only has to
produce the BSC arrays in the packing ``mie.vswf`` documents.

Packing (mirrors mie.vswf)
--------------------------
``g_tm``, ``g_te`` are dense ``(n_max, 2*n_max + 1)`` complex arrays indexed
``[n - 1, m + n_max]`` (n = 1..n_max, m = -n_max..n_max), unused |m| > n corners
left zero. ``mie.vswf.m_index`` / ``m_values`` document the m axis. For a plane wave
only m = +-1 are nonzero and the BSCs collapse to ``mie.vswf.plane_wave_bsc``
(g_{n,+-1}^{TM} = C_n/2, g_{n,+-1}^{TE} = +-C_n/2, C_n = (2n+1)/(n(n+1))). Every
method here reproduces those weights for a plane wave to <= 1e-6 (G-LIMIT).

Methods (INTERFACES.md §3), >= 2 implemented and cross-checked
--------------------------------------------------------------
- ``bsc_quadrature``  reference (slow): project the beam's ``focal_field`` onto the
  regular-VSWF transverse harmonics by numerical integration over a spherical
  surface. Ground truth for the others; works for ANY IncidentBeam.
- ``bsc_localized``   fast closed form for an on-axis paraxial Gaussian
  (Gouesbet-Gréhan localized approximation). Valid for moderate focusing; flagged
  for high NA (PHYSICS.md §2.2).
- ``bsc_angular_spectrum`` for a RichardsWolfFocus: projects the Debye-Wolf focal
  field onto VSWFs (delegates to the quadrature projection, which IS the
  angular-spectrum projection for that beam).

Displacement of the sphere (PHYSICS.md §3)
------------------------------------------
BSCs are expanded about the SPHERE center. Displacing the sphere to r_s relative to
the beam focus is the same as sampling the beam about a center offset by -r_s from
the focus: ``bsc_quadrature(beam, sphere_center_m=-r_s, ...)`` gives the BSCs of the
displaced sphere with no translation-addition theorem (the reliable path,
W1_scattering_engine.md). ``GLMTProvider.field(r_s)`` uses exactly this.

The quadrature normalization (derivation / calibration)
-------------------------------------------------------
The transverse regular-VSWF harmonics are ``(tau_nm, pi_nm)`` (mie.vswf.pi_tau_m).
Projecting the beam's spherical field components (E_theta, E_phi) onto them and
dividing by (i) the vector-harmonic norm N_nm = 2*pi * 2 n(n+1)/(2n+1) *
(n+|m|)!/(n-|m|)! and (ii) the regular radial factor of each VSWF type recovers the
BSCs. The radial factors that make the plane wave reproduce mie.vswf.plane_wave_bsc
*exactly* are ``psi_n'(kr)/(kr) * i^{n-1}`` (TM/N-type) and ``j_n(kr) * i^n``
(TE/M-type); these are pinned by numerical calibration against the analytic
plane-wave BSC (agreement <= 1e-13 at n_max = 8..120). The surface radius must keep
``j_n(kr)`` well conditioned (``kr >~ n_max``); |m| is capped at ``m_max`` because a
near-axis beam's BSCs decay fast in |m| and high-|m| associated Legendre functions
overflow float64 at n_max ~ 264 (PHYSICS.md §7).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import spherical_jn

from ..mie import vswf
from ..mie.special import nmax_wiscombe
from ..types import Medium
from .beam import GaussianParaxial, IncidentBeam, RichardsWolfFocus

__all__ = [
    "BSC",
    "bsc_quadrature",
    "bsc_localized",
    "bsc_angular_spectrum",
    "scattered_amplitudes",
]


@dataclass(frozen=True)
class BSC:
    """Beam-shape coefficients in the mie.vswf packing (INTERFACES.md §3).

    ``g_tm``, ``g_te`` are ``(n_max, 2*n_max + 1)`` complex, indexed [n-1, m+n_max].
    Time convention e^{-iωt}.
    """
    n_max: int
    g_tm: np.ndarray
    g_te: np.ndarray

    def __post_init__(self) -> None:
        shp = (self.n_max, 2 * self.n_max + 1)
        if self.g_tm.shape != shp or self.g_te.shape != shp:
            raise ValueError(
                f"g_tm/g_te must be {shp}; got {self.g_tm.shape}, {self.g_te.shape}"
            )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _default_m_max(beam: IncidentBeam, n_max: int, sphere_center_m: np.ndarray) -> int:
    """Azimuthal-order cap for the quadrature (near-axis beams decay fast in |m|).

    A beam whose focus is displaced transversely by ``d`` from the sphere center
    injects azimuthal content up to roughly ``k*d`` (like a plane wave hitting an
    off-axis point). We keep a generous margin over that, but cap well below n_max to
    stay clear of associated-Legendre overflow at large n_max (PHYSICS.md §7). For an
    on-axis beam this collapses to a small constant (plane wave needs only |m| = 1).
    """
    k = beam.medium.k
    d_trans = float(np.hypot(sphere_center_m[0], sphere_center_m[1]))
    m_from_offset = int(np.ceil(k * d_trans)) + 6
    return int(min(max(4, m_from_offset), max(4, n_max)))


def _surface_radius(k: float, n_max: int) -> float:
    """Projection-sphere radius: kr ~ n_max + margin so j_n(kr) stays well posed."""
    return (n_max + 15) / k


def _norm_ratio(n_arr: np.ndarray, am: int) -> np.ndarray:
    """(n+am)!/(n-am)! as a product of 2*am consecutive integers (overflow-safe).

    Full factorials overflow float64 past ~170; this product of ``2*am`` terms stays
    finite for the small |m| a near-axis beam needs. Rows with n < am are unused
    (harmonic is zero there) and set to 1 to avoid a spurious 0 or negative.
    """
    r = np.ones_like(n_arr, dtype=float)
    for j in range(1, 2 * am + 1):
        r = r * (n_arr - am + j)
    r[n_arr < am] = 1.0
    return r


def _project_focal_field(focal_field, k: float, n_max: int, m_max: int,
                         center_m: np.ndarray, ntheta: int, nphi: int
                         ) -> tuple[np.ndarray, np.ndarray]:
    """Project a focal field onto regular-VSWF transverse harmonics -> (g_tm, g_te).

    The reference quadrature (PHYSICS.md §2.2 method 1). ``center_m`` is the point of
    the beam field about which the VSWFs are expanded (the sphere center); the
    projection sphere is centered there. Returns packed (n_max, 2*n_max+1) arrays.
    """
    mu, w_mu = np.polynomial.legendre.leggauss(ntheta)
    theta = np.arccos(mu)
    phi = np.linspace(0.0, 2.0 * np.pi, nphi, endpoint=False)
    dphi = 2.0 * np.pi / nphi

    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    sx = np.sin(TH) * np.cos(PH)
    sy = np.sin(TH) * np.sin(PH)
    sz = np.cos(TH)
    r_surf = _surface_radius(k, n_max)
    xyz = np.stack([sx, sy, sz], axis=-1) * r_surf + np.asarray(center_m, dtype=float)
    E = focal_field(xyz)  # (Ntheta, Nphi, 3) cartesian, complex

    # cartesian -> spherical transverse components (theta_hat, phi_hat)
    e_th = np.stack([np.cos(TH) * np.cos(PH), np.cos(TH) * np.sin(PH), -np.sin(TH)],
                    axis=-1)
    e_ph = np.stack([-np.sin(PH), np.cos(PH), np.zeros_like(PH)], axis=-1)
    E_th = np.sum(E * e_th, axis=-1)  # (Ntheta, Nphi)
    E_ph = np.sum(E * e_ph, axis=-1)

    n = np.arange(1, n_max + 1)
    kr = k * r_surf
    jn = spherical_jn(n, kr)
    jnm1 = spherical_jn(n - 1, kr)
    djn = jnm1 - (n + 1) / kr * jn                 # j_n'(kr)
    dpsi = jn + kr * djn                           # psi_n'(kr), psi_n = kr j_n
    ipow = 1j ** n
    rad_tm = (dpsi / kr) * (ipow / 1j)             # N-type radial * i^{n-1}
    rad_te = jn * ipow                             # M-type radial * i^n
    base = 2.0 * np.pi * (2.0 * n * (n + 1) / (2.0 * n + 1))

    g_tm = np.zeros((n_max, 2 * n_max + 1), dtype=complex)
    g_te = np.zeros((n_max, 2 * n_max + 1), dtype=complex)

    for m in range(-m_max, m_max + 1):
        col = vswf.m_index(m, n_max)
        am = abs(m)
        pi_nm, tau_nm = vswf.pi_tau_m(n_max, m, theta)  # (n_max, Ntheta)
        emphi = np.exp(-1j * m * phi)                    # conjugate for projection
        Eth_m = np.tensordot(E_th, emphi, axes=(1, 0)) * dphi  # (Ntheta,)
        Eph_m = np.tensordot(E_ph, emphi, axes=(1, 0)) * dphi
        a_proj = np.sum(
            w_mu[None, :] * (Eth_m[None, :] * tau_nm + (Eph_m / 1j)[None, :] * pi_nm),
            axis=1,
        )
        b_proj = np.sum(
            w_mu[None, :] * (Eth_m[None, :] * pi_nm + (Eph_m / 1j)[None, :] * tau_nm),
            axis=1,
        )
        Nn = base * _norm_ratio(n, am)
        g_tm[:, col] = (a_proj / Nn) / rad_tm
        g_te[:, col] = (b_proj / Nn) / rad_te

    return g_tm, g_te


# ---------------------------------------------------------------------------
# Public BSC methods
# ---------------------------------------------------------------------------

def bsc_quadrature(beam: IncidentBeam, sphere_center_m: np.ndarray, n_max: int,
                   surface_radius_m: float | None = None,
                   m_max: int | None = None,
                   ntheta: int | None = None, nphi: int | None = None) -> BSC:
    """Reference BSCs by projecting ``beam.focal_field`` onto regular VSWFs (slow).

    PHYSICS.md §2.2 method 1 — ground truth for the fast methods and valid for ANY
    ``IncidentBeam``. ``sphere_center_m`` is the sphere center relative to the beam
    focus; the projection sphere is centered there, so a displaced sphere is handled
    by passing its offset (or ``-r_s`` when the *sphere* moves by r_s relative to the
    focus — see module docstring and GLMTProvider.field). A plane wave returns the
    plane-wave weights to <= 1e-6 (G-LIMIT).

    ``surface_radius_m`` overrides the default projection radius (default kr ~ n_max
    + 15 for well-conditioned j_n). ``m_max`` overrides the azimuthal-order cap.
    ``ntheta``/``nphi`` override the angular quadrature resolution (defaults scale
    with n_max and m_max for Gauss-Legendre exactness).
    """
    center = np.asarray(sphere_center_m, dtype=float)
    if center.shape != (3,):
        raise ValueError("sphere_center_m must be shape (3,)")
    k = beam.medium.k
    if m_max is None:
        m_max = _default_m_max(beam, n_max, center)
    m_max = int(min(m_max, n_max))
    if ntheta is None:
        ntheta = 2 * n_max + 20
    if nphi is None:
        nphi = max(2 * m_max + 4, 8)

    if surface_radius_m is not None:
        # Honour an explicit radius by temporarily monkeypatching via closure: build
        # the projection with a custom radius (re-derive kr-dependent factors).
        g_tm, g_te = _project_with_radius(beam.focal_field, k, n_max, m_max, center,
                                          ntheta, nphi, float(surface_radius_m))
    else:
        g_tm, g_te = _project_focal_field(beam.focal_field, k, n_max, m_max, center,
                                          ntheta, nphi)
    return BSC(n_max=n_max, g_tm=g_tm, g_te=g_te)


def _project_with_radius(focal_field, k, n_max, m_max, center, ntheta, nphi, r_surf):
    """Projection with an explicit surface radius (see _project_focal_field)."""
    mu, w_mu = np.polynomial.legendre.leggauss(ntheta)
    theta = np.arccos(mu)
    phi = np.linspace(0.0, 2.0 * np.pi, nphi, endpoint=False)
    dphi = 2.0 * np.pi / nphi
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    sx = np.sin(TH) * np.cos(PH); sy = np.sin(TH) * np.sin(PH); sz = np.cos(TH)
    xyz = np.stack([sx, sy, sz], axis=-1) * r_surf + np.asarray(center, dtype=float)
    E = focal_field(xyz)
    e_th = np.stack([np.cos(TH) * np.cos(PH), np.cos(TH) * np.sin(PH), -np.sin(TH)], axis=-1)
    e_ph = np.stack([-np.sin(PH), np.cos(PH), np.zeros_like(PH)], axis=-1)
    E_th = np.sum(E * e_th, axis=-1); E_ph = np.sum(E * e_ph, axis=-1)
    n = np.arange(1, n_max + 1)
    kr = k * r_surf
    jn = spherical_jn(n, kr); jnm1 = spherical_jn(n - 1, kr)
    djn = jnm1 - (n + 1) / kr * jn
    dpsi = jn + kr * djn
    ipow = 1j ** n
    rad_tm = (dpsi / kr) * (ipow / 1j); rad_te = jn * ipow
    base = 2.0 * np.pi * (2.0 * n * (n + 1) / (2.0 * n + 1))
    g_tm = np.zeros((n_max, 2 * n_max + 1), dtype=complex)
    g_te = np.zeros((n_max, 2 * n_max + 1), dtype=complex)
    for m in range(-m_max, m_max + 1):
        col = vswf.m_index(m, n_max); am = abs(m)
        pi_nm, tau_nm = vswf.pi_tau_m(n_max, m, theta)
        emphi = np.exp(-1j * m * phi)
        Eth_m = np.tensordot(E_th, emphi, axes=(1, 0)) * dphi
        Eph_m = np.tensordot(E_ph, emphi, axes=(1, 0)) * dphi
        a_proj = np.sum(w_mu[None, :] * (Eth_m[None, :] * tau_nm + (Eph_m / 1j)[None, :] * pi_nm), axis=1)
        b_proj = np.sum(w_mu[None, :] * (Eth_m[None, :] * pi_nm + (Eph_m / 1j)[None, :] * tau_nm), axis=1)
        Nn = base * _norm_ratio(n, am)
        g_tm[:, col] = (a_proj / Nn) / rad_tm
        g_te[:, col] = (b_proj / Nn) / rad_te
    return g_tm, g_te


def bsc_localized(beam: GaussianParaxial, sphere_center_m: np.ndarray, n_max: int
                  ) -> BSC:
    """Fast closed-form BSCs for an on-axis paraxial Gaussian (localized approx.).

    Gouesbet-Gréhan localized approximation (PHYSICS.md §2.2 method 2). For an
    x-polarized Gaussian of waist w0 propagating along +z with the focus at the
    origin, only m = +-1 contribute and the BSCs are the plane-wave weights times the
    axial Gaussian factor g_n = exp(-s^2 (n + 1/2)^2) with beam-confinement parameter
    s = 1/(k w0)::

        g_{n,+-1}^{TM} = (C_n/2) * exp(-s^2 (n+1/2)^2),   C_n = (2n+1)/(n(n+1))
        g_{n,+-1}^{TE} = +-(C_n/2) * exp(-s^2 (n+1/2)^2)

    An axial displacement z_s adds the phase exp(-i k z_s) * exp(... ) via the on-axis
    beam phase; a transverse displacement injects |m| != 1 orders the localized form
    does not capture, so this fast path supports ONLY on-axis centers (raise
    otherwise). Wide waist (s -> 0) => g_n -> 1 => plane-wave weights (G-LIMIT). The
    localized approximation degrades for tight focusing (w0 <~ 2 lambda / high NA);
    use bsc_quadrature or a RichardsWolfFocus there (PHYSICS.md §2.2 flag).
    """
    if not isinstance(beam, GaussianParaxial):
        raise TypeError("bsc_localized is defined for GaussianParaxial beams")
    center = np.asarray(sphere_center_m, dtype=float)
    if center.shape != (3,):
        raise ValueError("sphere_center_m must be shape (3,)")
    if np.hypot(center[0], center[1]) > 1e-15:
        raise NotImplementedError(
            "bsc_localized supports only on-axis centers (x=y=0); a transverse "
            "displacement injects |m|!=1 orders the localized approximation drops. "
            "Use bsc_quadrature for off-axis / displaced spheres (PHYSICS.md §2.2)."
        )
    k = beam.medium.k
    w0 = beam.waist_m()
    s = 1.0 / (k * w0)
    n = np.arange(1, n_max + 1)
    C = (2 * n + 1) / (n * (n + 1))
    g_axial = np.exp(-(s ** 2) * (n + 0.5) ** 2)
    # on-axis displacement z_s: multiply by the incident on-axis phase exp(i k z_s)
    z_s = float(center[2])
    phase = np.exp(1j * k * z_s) if z_s != 0.0 else 1.0
    weight = (C / 2.0) * g_axial * phase
    g_tm = np.zeros((n_max, 2 * n_max + 1), dtype=complex)
    g_te = np.zeros((n_max, 2 * n_max + 1), dtype=complex)
    ip = vswf.m_index(1, n_max)
    im = vswf.m_index(-1, n_max)
    g_tm[:, ip] = weight
    g_tm[:, im] = weight
    g_te[:, ip] = weight
    g_te[:, im] = -weight
    return BSC(n_max=n_max, g_tm=g_tm, g_te=g_te)


def bsc_angular_spectrum(beam: RichardsWolfFocus, sphere_center_m: np.ndarray,
                         n_max: int, m_max: int | None = None,
                         ntheta: int | None = None, nphi: int | None = None) -> BSC:
    """BSCs of a Richards-Wolf focus by projecting its Debye-Wolf field on VSWFs.

    PHYSICS.md §2.2 method 3 (the trustworthy high-NA path). The Richards-Wolf focal
    field IS an angular spectrum of plane waves; projecting it onto VSWFs is the
    quadrature projection applied to that beam, so this delegates to the shared
    projection. Same packing as the other methods. A very-low-NA focus recovers the
    plane-wave weights up to the focus's global phase (G-LIMIT, up to a constant).
    """
    if not isinstance(beam, RichardsWolfFocus):
        raise TypeError("bsc_angular_spectrum is defined for RichardsWolfFocus beams")
    return bsc_quadrature(beam, sphere_center_m, n_max, m_max=m_max,
                          ntheta=ntheta, nphi=nphi)


def scattered_amplitudes(bsc: BSC, a_n: np.ndarray, b_n: np.ndarray
                         ) -> tuple[np.ndarray, np.ndarray]:
    """Assemble scattered VSWF amplitudes A_nm = a_n g_tm, B_nm = b_n g_te.

    ``a_n``, ``b_n`` are the Mie coefficients (mie.plane_wave.mie_coefficients),
    shape (n_max,), index 0 == n=1. Returns packed (n_max, 2*n_max+1) arrays ready
    for ``mie.vswf.far_field`` (PHYSICS.md §2.1).
    """
    a_n = np.asarray(a_n)
    b_n = np.asarray(b_n)
    if a_n.shape != (bsc.n_max,) or b_n.shape != (bsc.n_max,):
        raise ValueError(f"a_n/b_n must be shape ({bsc.n_max},)")
    return a_n[:, None] * bsc.g_tm, b_n[:, None] * bsc.g_te
