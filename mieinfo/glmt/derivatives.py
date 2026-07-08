"""Displacement derivatives dE_s/dr_{s,j} of the GLMT scattered far field
(PHYSICS.md §3, INTERFACES.md §4). Time convention e^{-iωt}; SI units.

The information source is the sphere's displacement r_s. GLMTProvider.field writes the
displaced far field as two exact pieces (glmt.scatter.GLMTProvider):

    E_s(s_hat; r_s) = FarField[ a_n g_{n,TM}^m(r_s),  b_n g_{n,TE}^m(r_s) ]
                      * exp(-i k s_hat . r_s)

where the BSCs g^{TM/TE}(r_s) are the incident beam expanded about the sphere center
r_s (glmt.bsc), and the scalar exp(-i k s_hat . r_s) is the sphere->detector
propagation phase. Differentiating by the product rule,

    dE_s/dr_{s,j} = [ FarField(dg/dr_j)  +  FarField(g) * (-i k s_hat_j) ]
                    * exp(-i k s_hat . r_s).                                (*)

Two regimes (PHYSICS.md §3):

- **Plane wave (exact, §3.1).** For a plane wave dg/dr_j reduces to a global phase and
  (*) collapses to the exact phase-gradient derivative
  ``dE_s/dr_j = i k (k_hat - s_hat)_j E_s`` — fully analytic, machine precise. This is
  the same result glmt.scatter uses for PlaneWaveProvider; here it is reached through
  the GLMT reconstruction so GLMTProvider has an analytic derivative too.

- **Focused beam (§3.2).** The BSC-gradient term dg/dr_j is the derivative of the
  incident-beam expansion under displacement. The fully-analytic route is the VSWF
  translation-addition theorem (Cruzan/Stein), which is HEAVY at n_max ~ 264 and is
  the stretch item (glmt.translation, partial — see module note there). This module
  instead computes dg/dr_j by a high-accuracy central difference of the (cheap,
  cached) BSC while keeping the physical translation-phase term analytic. That
  "semi-analytic" derivative matches the full finite-difference field derivative
  (glmt.scatter.field_derivative(method='finite_difference'), the G-DERIV reference)
  to <= 1e-6 over |r_s| <= 0.3 lambda. The purely-analytic translation-addition
  derivative is left pending (gaps_or_issues); the validated displacement path is the
  quadrature-BSC field + this derivative.
"""
from __future__ import annotations

import numpy as np

from ..mie import plane_wave as pw
from ..mie import vswf
from ..types import AngularGrid, FieldDerivative, Sphere
from .beam import IncidentBeam, PlaneWave

__all__ = ["field_derivative_analytic"]


def _shat_components(grid: AngularGrid) -> np.ndarray:
    """s_hat_j on the grid; shape (3, Ntheta, Nphi)."""
    TH = grid.theta[:, None]
    PH = grid.phi[None, :]
    sinT = np.sin(TH)
    sx = sinT * np.cos(PH)
    sy = sinT * np.sin(PH)
    sz = np.cos(TH) * np.ones_like(PH)
    return np.stack([sx, sy, sz])


def field_derivative_analytic(provider, grid: AngularGrid, sphere: Sphere,
                              beam: IncidentBeam, r_s_m: np.ndarray,
                              n_max: int | None = None) -> FieldDerivative:
    """Analytic dE_s/dr_j for a GLMTProvider (INTERFACES.md §4, PHYSICS.md §3).

    Plane-wave beams use the exact phase-gradient derivative (§3.1). Focused beams use
    the product-rule decomposition (*) with an analytic translation-phase term and a
    central-difference BSC gradient (semi-analytic; validated vs finite_difference to
    <= 1e-6, G-DERIV). Returns a FieldDerivative with axis 0 = (x, y, z).
    """
    r_s = np.asarray(r_s_m, dtype=float)
    if r_s.shape != (3,):
        raise ValueError("r_s_m must be shape (3,)")
    k = beam.medium.k

    # --- exact plane-wave path (fully analytic, machine precise) ---
    if isinstance(beam, PlaneWave):
        f = provider.field(grid, sphere, beam, r_s, n_max)
        TH = grid.theta[:, None]
        PH = grid.phi[None, :]
        sinT = np.sin(TH)
        khat_minus_shat = np.stack([
            -sinT * np.cos(PH),
            -sinT * np.sin(PH),
            1.0 - np.cos(TH) * np.ones_like(PH),
        ])
        fac = 1j * k * khat_minus_shat
        return FieldDerivative(grid=grid,
                               dE_theta=fac * f.E_theta[None],
                               dE_phi=fac * f.E_phi[None])

    # --- focused-beam path: product rule (*) ---
    n_max = provider._resolve_nmax(sphere, beam, n_max)
    x = pw.size_parameter(sphere, beam.medium)
    a_n, b_n = pw.mie_coefficients(sphere.m, x, n_max)

    bsc0 = provider._bsc(beam, r_s, n_max)
    vf0 = vswf.far_field(a_n[:, None] * bsc0.g_tm, b_n[:, None] * bsc0.g_te, grid)

    shat = _shat_components(grid)                      # (3, Ntheta, Nphi)
    sdotr = shat[0] * r_s[0] + shat[1] * r_s[1] + shat[2] * r_s[2]
    tphase = np.exp(-1j * k * sdotr)                   # (Ntheta, Nphi)

    lam = 2.0 * np.pi / k
    h = 1e-4 * lam
    shape = (3,) + grid.w_solid.shape
    dE_theta = np.empty(shape, dtype=complex)
    dE_phi = np.empty(shape, dtype=complex)

    for j in range(3):
        cp = r_s.copy(); cp[j] += h
        cm = r_s.copy(); cm[j] -= h
        bp = provider._bsc(beam, cp, n_max)
        bm = provider._bsc(beam, cm, n_max)
        dg_tm = (bp.g_tm - bm.g_tm) / (2.0 * h)
        dg_te = (bp.g_te - bm.g_te) / (2.0 * h)
        vf_dbsc = vswf.far_field(a_n[:, None] * dg_tm, b_n[:, None] * dg_te, grid)
        # product rule (*): BSC-gradient term + translation-phase-gradient term
        dE_theta[j] = (vf_dbsc.E_theta + vf0.E_theta * (-1j * k * shat[j])) * tphase
        dE_phi[j] = (vf_dbsc.E_phi + vf0.E_phi * (-1j * k * shat[j])) * tphase

    return FieldDerivative(grid=grid, dE_theta=dE_theta, dE_phi=dE_phi)
