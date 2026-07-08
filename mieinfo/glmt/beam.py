"""Incident beams (INTERFACES.md §2). Time convention e^{-iωt}.

M0 seeds only `PlaneWave` (the validation baseline that `PlaneWaveProvider` wraps).
W1b adds `GaussianParaxial`, `RichardsWolfFocus`, and `lab_from_beam_frame`. Each beam
carries its OWN Medium (hence its own wavelength): detection beams use a 532 nm Medium,
the 1064 nm trap Medium is separate. In a beam's own formulae k_hat = +z and the default
polarization is +x; a beam's lab-frame orientation is applied by the channel (§6).
"""
from __future__ import annotations

from typing import Protocol, runtime_checkable

import numpy as np

from ..types import Medium


@runtime_checkable
class IncidentBeam(Protocol):
    """Plane wave / paraxial Gaussian / Richards-Wolf focus behind one API."""
    medium: Medium
    polarization: np.ndarray  # Jones vector in the (x, y) transverse plane

    def focal_field(self, xyz: np.ndarray) -> np.ndarray:
        """E(r) at cartesian points xyz (..., 3) -> (..., 3) complex."""
        ...

    def waist_m(self) -> float:
        """Characteristic transverse scale (inf for a plane wave)."""
        ...


class PlaneWave:
    """x-polarized (frame convention) plane wave propagating along +z."""

    def __init__(self, medium: Medium, polarization=(1.0, 0.0)):
        self.medium = medium
        self.polarization = np.asarray(polarization, dtype=complex)

    def focal_field(self, xyz: np.ndarray) -> np.ndarray:
        xyz = np.asarray(xyz, dtype=float)
        phase = np.exp(1j * self.medium.k * xyz[..., 2])   # e^{ikz}, e^{-iωt}
        E = np.zeros(xyz.shape[:-1] + (3,), dtype=complex)
        E[..., 0] = self.polarization[0] * phase
        E[..., 1] = self.polarization[1] * phase
        return E

    def waist_m(self) -> float:
        return float("inf")


class GaussianParaxial:
    """Paraxial (fundamental TEM00) Gaussian beam, waist at the origin, along +z.

    Time convention e^{-iωt}. In the beam frame k_hat = +z and the transverse
    polarization is the Jones vector ``polarization`` (default +x). The scalar
    envelope uses the complex beam parameter q(z) = z + i z_R with z_R = k w0^2 / 2::

        E_t(r) = pol * (i z_R / q(z)) * exp(i k rho^2 / (2 q(z))) * exp(i k z)

    so |E| = 1 at the focus (rho = z = 0). This is the standard paraxial field; it is
    the reference the localized-approximation BSCs (bsc_localized) target, and its
    wide-waist limit (w0 >> lambda) reproduces the plane-wave BSCs (G-LIMIT). It is
    NOT valid at high NA (w0 ~ lambda): use RichardsWolfFocus there (PHYSICS.md §2.2).
    The paraxial field carries only transverse (x, y) components — no longitudinal
    E_z — which is the leading approximation; higher-order (Lax) corrections are
    omitted deliberately (flagged: bsc validity degrades for w0 <~ 2 lambda).
    """

    def __init__(self, medium: Medium, waist_m: float, polarization=(1.0, 0.0)):
        if waist_m <= 0:
            raise ValueError("waist_m must be positive")
        self.medium = medium
        self._w0 = float(waist_m)
        pol = np.asarray(polarization, dtype=complex)
        n = np.linalg.norm(pol)
        self.polarization = pol / n if n > 0 else pol

    def focal_field(self, xyz: np.ndarray) -> np.ndarray:
        xyz = np.asarray(xyz, dtype=float)
        k = self.medium.k
        z_R = 0.5 * k * self._w0 ** 2
        x = xyz[..., 0]
        y = xyz[..., 1]
        z = xyz[..., 2]
        rho2 = x * x + y * y
        q = z + 1j * z_R
        amp = (1j * z_R / q) * np.exp(1j * k * rho2 / (2.0 * q)) * np.exp(1j * k * z)
        E = np.zeros(xyz.shape[:-1] + (3,), dtype=complex)
        E[..., 0] = self.polarization[0] * amp
        E[..., 1] = self.polarization[1] * amp
        return E

    def waist_m(self) -> float:
        return self._w0


class RichardsWolfFocus:
    """High-NA aplanatic focus of a collimated, x-polarized beam (Richards & Wolf
    1959; PHYSICS.md §2.2, §2.3). Time convention e^{-iωt}; beam frame k_hat = +z.

    An aplanatic objective of numerical aperture ``NA`` (= n sin alpha_max, vacuum
    here so NA = sin alpha_max <= 1) focuses an incident collimated field. The focal
    field is the Debye-Wolf angular-spectrum integral over the convergent cone of
    plane waves. For an x-polarized input with an amplitude/apodization profile
    l(theta) (aplanatic sqrt(cos theta) factor and a Gaussian ``filling_factor``
    truncation), the Cartesian focal field at r = (rho, varphi, z) is

        E_x = -i A (I0 + I2 cos 2varphi)
        E_y = -i A  I2 sin 2varphi
        E_z = -2 A I1 cos varphi

    with (Novotny & Hecht, Principles of Nano-Optics, eq. 3.66)

        I0 = ∫ l(t) sqrt(cos t) sin t (1 + cos t) J0(k rho sin t) e^{i k z cos t} dt
        I1 = ∫ l(t) sqrt(cos t) sin^2 t         J1(k rho sin t) e^{i k z cos t} dt
        I2 = ∫ l(t) sqrt(cos t) sin t (1 - cos t) J2(k rho sin t) e^{i k z cos t} dt

    over t in [0, alpha_max]. l(t) = exp(-(sin t / (f0 sin alpha_max))^2) is the
    Gaussian filling of the pupil (``filling_factor`` f0; f0 -> inf is uniform
    illumination). The prefactor A is chosen so |E| = 1 at the geometric focus. This
    is the trustworthy high-NA field (the paraxial Gaussian is not valid at trap NA);
    its BSCs come from projecting this field onto VSWFs (bsc_angular_spectrum /
    bsc_quadrature). In the low-NA limit it collapses to the plane-wave / paraxial
    result (G-LIMIT). ``polarization`` rotates the transverse input Jones vector;
    only linear x is modelled exactly (y is the same field rotated 90 deg).

    VALIDATION CAVEAT (do not trust at high NA). This focused-beam path is gate-validated
    only in its plane-wave / low-NA limit. An independent cross-check against the Moore Lab
    dissertation (docs/dissertation_comparison.md) shows that at NA≈0.63 the GLMT engine does
    NOT reproduce the reference focused-beam information-density numbers — it misses the
    Gouy-phase axial-forward sensitivity and produces spurious transverse-backscatter. Use the
    plane-wave provider for weakly-focused beams; the high-NA focused-beam path needs fixing
    and validation before quantitative use.
    """

    def __init__(self, medium: Medium, NA: float, filling_factor: float = 1.0,
                 polarization=(1.0, 0.0), n_quad: int = 200):
        if not (0.0 < NA <= 1.0):
            raise ValueError("NA must be in (0, 1] (vacuum medium, NA = sin alpha)")
        self.medium = medium
        self.NA = float(NA)
        self.filling_factor = float(filling_factor)
        self._n_quad = int(n_quad)
        pol = np.asarray(polarization, dtype=complex)
        nrm = np.linalg.norm(pol)
        self.polarization = pol / nrm if nrm > 0 else pol
        self._alpha_max = float(np.arcsin(self.NA))
        # Precompute the aplanatic quadrature nodes in theta over [0, alpha_max].
        t, w = np.polynomial.legendre.leggauss(self._n_quad)
        self._t = 0.5 * self._alpha_max * (t + 1.0)          # nodes
        self._wt = 0.5 * self._alpha_max * w                  # weights
        # Normalize amplitude so |E(0)| = 1 at the focus.
        self._A = 1.0
        E0 = self._focal_xyz(np.zeros((1, 3)))
        amp0 = np.sqrt(np.sum(np.abs(E0[0]) ** 2))
        if amp0 > 0:
            self._A = 1.0 / amp0

    # -- internal: field for +x input polarization, before applying Jones rotation --
    def _debye_integrals(self, rho: np.ndarray, z: np.ndarray
                         ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        from scipy.special import jv
        k = self.medium.k
        t = self._t
        w = self._wt
        st = np.sin(t)
        ct = np.cos(t)
        sina = np.sin(self._alpha_max)
        f0 = self.filling_factor
        apod = np.sqrt(ct)                        # aplanatic sqrt(cos) factor
        fill = np.exp(-(st / (f0 * sina)) ** 2)   # Gaussian pupil filling
        l_t = fill * apod * w                     # combined weight per node
        # rho, z are (...,) arrays; broadcast against nodes on a trailing axis.
        rho_b = np.asarray(rho)[..., None]
        z_b = np.asarray(z)[..., None]
        krsin = k * rho_b * st                    # (..., n_quad)
        phase = np.exp(1j * k * z_b * ct)         # (..., n_quad)
        base = l_t * st * phase
        I0 = np.sum(base * (1.0 + ct) * jv(0, krsin), axis=-1)
        I1 = np.sum(base * st * jv(1, krsin), axis=-1)
        I2 = np.sum(base * (1.0 - ct) * jv(2, krsin), axis=-1)
        return I0, I1, I2

    def _focal_xyz(self, xyz: np.ndarray) -> np.ndarray:
        """+x-polarized aplanatic focal field at cartesian points (unrotated)."""
        xyz = np.asarray(xyz, dtype=float)
        x = xyz[..., 0]
        y = xyz[..., 1]
        z = xyz[..., 2]
        rho = np.hypot(x, y)
        varphi = np.arctan2(y, x)
        I0, I1, I2 = self._debye_integrals(rho, z)
        A = self._A
        E = np.zeros(xyz.shape[:-1] + (3,), dtype=complex)
        E[..., 0] = -1j * A * (I0 + I2 * np.cos(2 * varphi))
        E[..., 1] = -1j * A * (I2 * np.sin(2 * varphi))
        E[..., 2] = -2.0 * A * I1 * np.cos(varphi)
        return E

    def focal_field(self, xyz: np.ndarray) -> np.ndarray:
        """E(r) at cartesian points xyz (..., 3) -> (..., 3) complex.

        Computed for +x input polarization, then the transverse Jones vector rotates
        the field about +z (only real linear polarizations are exact)."""
        E = self._focal_xyz(xyz)
        px, py = self.polarization
        if abs(py) < 1e-15 and abs(px - 1.0) < 1e-15:
            return E
        # Rotate the field by the angle of the (real) Jones vector about z.
        ang = np.arctan2(np.real(py), np.real(px))
        c, s = np.cos(ang), np.sin(ang)
        # Evaluate the +x field in the rotated coordinate frame, then rotate back.
        xyz = np.asarray(xyz, dtype=float)
        R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        xyz_rot = xyz @ R  # rotate points by -ang (R^T applied on the right)
        E_rot = self._focal_xyz(xyz_rot)
        return E_rot @ R.T

    def waist_m(self) -> float:
        """Diffraction-limited focal spot scale ~ lambda / (pi NA)."""
        lam = 2.0 * np.pi / self.medium.k
        return lam / (np.pi * self.NA)


def lab_from_beam_frame(propagation_lab: np.ndarray,
                        polarization_lab: np.ndarray) -> np.ndarray:
    """3x3 rotation taking a beam-frame direction to the lab frame: maps the beam's
    +z to `propagation_lab` and +x to `polarization_lab` (Gram-Schmidt-orthogonalized
    against the propagation axis). Columns are the lab images of beam +x, +y, +z
    (INTERFACES.md §2, PHYSICS.md §4.5). Used by the channel to place a beam's pattern
    in the lab frame and to combine channels."""
    z = np.asarray(propagation_lab, dtype=float)
    z = z / np.linalg.norm(z)
    x = np.asarray(polarization_lab, dtype=float)
    x = x - np.dot(x, z) * z
    nx = np.linalg.norm(x)
    if nx < 1e-12:
        raise ValueError("polarization_lab must not be parallel to propagation_lab")
    x = x / nx
    y = np.cross(z, x)
    return np.column_stack([x, y, z])
