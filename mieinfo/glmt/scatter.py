"""FieldProvider — the cross-track seam (orchestrator-owned, INTERFACES.md §4).

`info/` and `detect/` depend on THIS protocol, never on a concrete provider (they take
the provider as an argument). Satisfied by `PlaneWaveProvider` now, `GLMTProvider` (W1)
later. Time convention e^{-iωt}.

Far field for x-polarized plane-wave incidence (beam frame, PHYSICS.md §1.3):
    E_theta ∝ S2(θ) cosφ,   E_phi ∝ −S1(θ) sinφ,
displaced by r_s via the exact plane-wave phase gradient (PHYSICS.md §3.1):
    E_s(ŝ; r_s) = E_s(ŝ; 0) · exp[i k (k̂ − ŝ)·r_s],   k̂ = +z.
Lab-frame polarization/orientation is applied by the channel (§4.5), not here.
"""
from __future__ import annotations

from typing import Protocol

import numpy as np

from ..mie import plane_wave as pw
from ..mie import vswf
from ..mie.special import nmax_wiscombe
from ..types import AngularGrid, FieldDerivative, Sphere, VectorField
from .beam import IncidentBeam, PlaneWave


class FieldProvider(Protocol):
    def field(self, grid: AngularGrid, sphere: Sphere, beam: IncidentBeam,
              r_s_m: np.ndarray, n_max: int | None = None) -> VectorField:
        ...

    def q_sca(self, sphere: Sphere, beam: IncidentBeam, n_max: int | None = None) -> float:
        ...


def _khat_minus_shat(grid: AngularGrid) -> np.ndarray:
    """(k̂ − ŝ)_j on the grid, j = x,y,z; shape (3, Ntheta, Nphi). k̂ = +z."""
    TH = grid.theta[:, None]
    PH = grid.phi[None, :]
    sinT = np.sin(TH)
    sx = sinT * np.cos(PH)
    sy = sinT * np.sin(PH)
    sz = np.cos(TH) * np.ones_like(PH)
    return np.stack([-sx, -sy, 1.0 - sz])


def _shat_dot(grid: AngularGrid, r: np.ndarray) -> np.ndarray:
    """s_hat . r on the grid; shape (Ntheta, Nphi)."""
    TH = grid.theta[:, None]
    PH = grid.phi[None, :]
    sinT = np.sin(TH)
    sx = sinT * np.cos(PH)
    sy = sinT * np.sin(PH)
    sz = np.cos(TH) * np.ones_like(PH)
    return sx * r[0] + sy * r[1] + sz * r[2]


def _apply_translation_phase(e_theta: np.ndarray, e_phi: np.ndarray,
                             grid: AngularGrid, k: float, r_s: np.ndarray
                             ) -> tuple[np.ndarray, np.ndarray]:
    """Multiply a sphere-centered far field by exp(-i k s_hat . r_s) (source at r_s)."""
    phase = np.exp(-1j * k * _shat_dot(grid, r_s))
    return e_theta * phase, e_phi * phase


class PlaneWaveProvider:
    """Plane-wave Mie FieldProvider; r_s via the analytic phase gradient (PHYSICS §3.1)."""

    def field(self, grid: AngularGrid, sphere: Sphere, beam: IncidentBeam,
              r_s_m: np.ndarray, n_max: int | None = None) -> VectorField:
        pol = np.asarray(beam.polarization, dtype=complex)
        if not (abs(abs(pol[0]) - 1.0) < 1e-9 and abs(pol[1]) < 1e-9):
            raise NotImplementedError(
                "PlaneWaveProvider computes in the beam frame with +x polarization "
                "(PHYSICS.md §0, §4.5); lab-frame polarization is applied by the detection "
                f"channel, not here. Got beam.polarization={pol}."
            )
        x = pw.size_parameter(sphere, beam.medium)
        S1, S2 = pw.scattering_amplitudes(sphere.m, x, grid.theta, n_max)
        cph = np.cos(grid.phi)[None, :]
        sph = np.sin(grid.phi)[None, :]
        E_theta = S2[:, None] * cph
        E_phi = -S1[:, None] * sph
        r_s = np.asarray(r_s_m, dtype=float)
        if np.any(r_s):
            dot = np.tensordot(r_s, _khat_minus_shat(grid), axes=(0, 0))  # (Ntheta,Nphi)
            phase = np.exp(1j * beam.medium.k * dot)
            E_theta = E_theta * phase
            E_phi = E_phi * phase
        return VectorField(grid=grid, E_theta=E_theta, E_phi=E_phi)

    def q_sca(self, sphere: Sphere, beam: IncidentBeam, n_max: int | None = None) -> float:
        x = pw.size_parameter(sphere, beam.medium)
        return pw.efficiencies(sphere.m, x, n_max)[1]


class GLMTProvider:
    """Generalized-Lorenz-Mie FieldProvider (INTERFACES.md §4, PHYSICS.md §2-3).

    Time convention e^{-iωt}. The scattered far field is the outgoing VSWF sum with
    partial-wave amplitudes A_nm = a_n g_{n,TM}^m, B_nm = b_n g_{n,TE}^m — the same
    Mie a_n, b_n (mie.plane_wave) times the beam's beam-shape coefficients (glmt.bsc),
    reconstructed by mie.vswf.far_field (PHYSICS.md §2.1). This is the physical trap
    model the PlaneWaveProvider approximates.

    Sphere displacement (PHYSICS.md §3.2, reliable path). Moving the sphere to r_s
    relative to the beam focus re-expands the incident beam about the new center. We
    realize this WITHOUT the translation-addition theorem in two exact pieces:
    (1) recompute the BSCs by quadrature for the beam sampled about the sphere center
        ``sphere_center_m = +r_s`` (this captures the incident phase/amplitude at the
        displaced sphere), and
    (2) multiply the reconstructed far field by the sphere->detector propagation phase
        ``exp(-i k s_hat . r_s)`` (the VSWFs are centered on the sphere, so a source at
        r_s radiates with that far-field phase).
    For a plane wave these two combine to the exact analytic phase shift
    ``exp(i k (k_hat - s_hat) . r_s)`` (PHYSICS.md §3.1) to machine precision, which is
    the correctness check for the displacement path. field(r_s) is therefore exact (to
    quadrature precision) for any beam, and the existing
    ``field_derivative(method='finite_difference')`` differentiates it directly.

    BSCs depend only on (beam, sphere_center, n_max) and are reused across the whole
    angular grid and across DOF, so they are cached (PHYSICS.md §7 cost lever). Auto-mode
    defaults to the EXACT quadrature BSC; the fast localized approximation is opt-in
    (``bsc_method='localized'``) and valid only for waist >> sphere, on-axis — it is NOT
    selected automatically because its error grows with sphere/waist (~74% at a≈0.8 w0).
    """

    def __init__(self, bsc_method: str = "auto", bsc_kwargs: dict | None = None):
        self.bsc_method = bsc_method
        self.bsc_kwargs = dict(bsc_kwargs or {})
        self._cache: dict = {}

    # -- BSC dispatch + cache ------------------------------------------------
    def _bsc(self, beam: IncidentBeam, center_m: np.ndarray, n_max: int):
        from . import bsc as bsc_mod
        from .beam import RichardsWolfFocus

        key = (id(beam), tuple(np.round(center_m, 18)), n_max, self.bsc_method,
               tuple(sorted(self.bsc_kwargs.items())))
        cached = self._cache.get(key)
        if cached is not None:
            return cached

        method = self.bsc_method
        if method == "auto":
            # Default to the EXACT quadrature BSC. The localized approximation is fast
            # but valid only when the beam varies slowly over the sphere (waist >>
            # sphere): its BSC error grows with sphere/waist (verified ~5% at a≈0.2 w0,
            # ~74% at a≈0.8 w0) and it drops the |m|!=1 orders for off-axis centers. So
            # auto-mode never selects it silently — request it explicitly
            # (bsc_method="localized") only in the wide-waist, on-axis regime.
            if isinstance(beam, RichardsWolfFocus):
                method = "angular_spectrum"
            else:
                method = "quadrature"

        if method == "localized":
            b = bsc_mod.bsc_localized(beam, center_m, n_max)
        elif method == "angular_spectrum":
            b = bsc_mod.bsc_angular_spectrum(beam, center_m, n_max, **self.bsc_kwargs)
        elif method == "quadrature":
            b = bsc_mod.bsc_quadrature(beam, center_m, n_max, **self.bsc_kwargs)
        else:
            raise ValueError(f"unknown bsc_method {method!r}")
        self._cache[key] = b
        return b

    def _resolve_nmax(self, sphere: Sphere, beam: IncidentBeam,
                      n_max: int | None) -> int:
        if n_max is not None:
            return int(n_max)
        x = pw.size_parameter(sphere, beam.medium)
        return nmax_wiscombe(x)

    def field(self, grid: AngularGrid, sphere: Sphere, beam: IncidentBeam,
              r_s_m: np.ndarray, n_max: int | None = None) -> VectorField:
        r_s = np.asarray(r_s_m, dtype=float)
        if r_s.shape != (3,):
            raise ValueError("r_s_m must be shape (3,)")
        n_max = self._resolve_nmax(sphere, beam, n_max)
        x = pw.size_parameter(sphere, beam.medium)
        a_n, b_n = pw.mie_coefficients(sphere.m, x, n_max)
        # (1) BSCs of the incident beam expanded about the displaced sphere center.
        bsc = self._bsc(beam, r_s, n_max)
        a_tm = a_n[:, None] * bsc.g_tm
        b_te = b_n[:, None] * bsc.g_te
        vf = vswf.far_field(a_tm, b_te, grid)
        # (2) sphere->detector propagation phase exp(-i k s_hat . r_s).
        if np.any(r_s):
            e_th, e_ph = _apply_translation_phase(vf.E_theta, vf.E_phi, grid,
                                                  beam.medium.k, r_s)
            return VectorField(grid=grid, E_theta=e_th, E_phi=e_ph)
        return vf

    def q_sca(self, sphere: Sphere, beam: IncidentBeam, n_max: int | None = None
              ) -> float:
        """Total scattered efficiency. Backaction is set by the TOTAL scattered power
        (PHYSICS.md §4.3); the plane-wave Q_sca (from a_n, b_n) is the leading value
        and is beam-independent for the same sphere, so it is used here (the focused
        beam re-weights the angular distribution, not the summed |a_n|^2+|b_n|^2 total
        used for Gamma_ba scaling)."""
        x = pw.size_parameter(sphere, beam.medium)
        return pw.efficiencies(sphere.m, x, n_max)[1]


def field_derivative(provider: FieldProvider, grid: AngularGrid, sphere: Sphere,
                     beam: IncidentBeam, r_s_m: np.ndarray,
                     method: str = "analytic", n_max: int | None = None) -> FieldDerivative:
    """dE_s/dr_j for j = x,y,z. 'analytic' (exact phase gradient for PlaneWaveProvider,
    PHYSICS §3.1) or 'finite_difference' (central, step ~1e-4·λ; the G-DERIV reference)."""
    r_s = np.asarray(r_s_m, dtype=float)

    if method == "analytic":
        if isinstance(provider, PlaneWaveProvider):
            f = provider.field(grid, sphere, beam, r_s, n_max)
            fac = 1j * beam.medium.k * _khat_minus_shat(grid)      # (3,Ntheta,Nphi)
            return FieldDerivative(grid=grid,
                                   dE_theta=fac * f.E_theta[None],
                                   dE_phi=fac * f.E_phi[None])
        if isinstance(provider, GLMTProvider):
            from . import derivatives as glmt_deriv
            return glmt_deriv.field_derivative_analytic(
                provider, grid, sphere, beam, r_s, n_max)
        raise NotImplementedError(
            "analytic field_derivative is implemented for PlaneWaveProvider (exact "
            "phase gradient, PHYSICS §3.1) and GLMTProvider (glmt.derivatives); got "
            f"{type(provider).__name__}."
        )

    if method == "finite_difference":
        lam_medium = 2.0 * np.pi / beam.medium.k
        h = 1e-4 * lam_medium
        shape = (3,) + grid.w_solid.shape
        dE_theta = np.empty(shape, dtype=complex)
        dE_phi = np.empty(shape, dtype=complex)
        for j in range(3):
            rp = r_s.copy(); rp[j] += h
            rm = r_s.copy(); rm[j] -= h
            fp = provider.field(grid, sphere, beam, rp, n_max)
            fm = provider.field(grid, sphere, beam, rm, n_max)
            dE_theta[j] = (fp.E_theta - fm.E_theta) / (2.0 * h)
            dE_phi[j] = (fp.E_phi - fm.E_phi) / (2.0 * h)
        return FieldDerivative(grid=grid, dE_theta=dE_theta, dE_phi=dE_phi)

    raise ValueError(f"unknown method {method!r} (use 'analytic' or 'finite_difference')")
