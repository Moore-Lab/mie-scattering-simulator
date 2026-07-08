"""Tests for mieinfo.mie.vswf — vector angular harmonics and the far-field VSWF
reconstruction. Time convention e^{-iωt} (PHYSICS.md §0).

Correctness anchor (W1a Definition of Done): the plane-wave VSWF far field must
reproduce mie.plane_wave.scattering_amplitudes S1/S2 to <= 1e-8 at several x. This
is the seam the future GLMTProvider (T1.7) consumes, and the plane-wave limit every
BSC method must reproduce (VALIDATION.md §3, G-LIMIT).
"""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.mie import plane_wave as pw
from mieinfo.mie import vswf
from mieinfo.mie.special import nmax_wiscombe, pi_tau
from mieinfo.types import AngularGrid

M_SILICA = 1.4607  # silica relative index at 532 nm detection light (PHYSICS.md §0)


# ---------------------------------------------------------------------------
# Generalized angular functions
# ---------------------------------------------------------------------------

def test_pi_tau_m1_matches_mie():
    """pi_tau_m(n_max, 1, theta) reproduces mie.special.pi_tau exactly (|m|=1)."""
    n_max = 40
    theta = np.linspace(0.05, np.pi - 0.05, 37)
    pi_ref, tau_ref = pi_tau(theta, n_max)
    pi_m, tau_m = vswf.pi_tau_m(n_max, 1, theta)
    # The two use different recurrences (mie's direct pi/tau vs this module's
    # associated-Legendre P then d/dtheta), so they agree only to rounding.
    assert np.allclose(pi_m, pi_ref, rtol=1e-10, atol=1e-10)
    assert np.allclose(tau_m, tau_ref, rtol=1e-10, atol=1e-10)


def test_pi_tau_m_sign_symmetry():
    """pi is odd in sign(m) (pi_{n,-m} = -pi_{nm}); tau is even (tau_{n,-m}=tau_{nm})."""
    n_max = 25
    theta = np.linspace(0.1, np.pi - 0.1, 19)
    for m in (1, 2, 3):
        pi_p, tau_p = vswf.pi_tau_m(n_max, m, theta)
        pi_n, tau_n = vswf.pi_tau_m(n_max, -m, theta)
        assert np.allclose(pi_n, -pi_p, rtol=0, atol=1e-11)
        assert np.allclose(tau_n, tau_p, rtol=0, atol=1e-11)


def test_pi_tau_m_shape_and_below_order_zero():
    """Shape is (n_max,)+theta.shape and rows with n < |m| are zero."""
    n_max = 12
    theta = np.linspace(0.2, 2.9, 8)
    m = 5
    pi_m, tau_m = vswf.pi_tau_m(n_max, m, theta)
    assert pi_m.shape == (n_max,) + theta.shape
    assert tau_m.shape == (n_max,) + theta.shape
    # rows n = 1..m-1 (index 0..m-2) must vanish
    assert np.allclose(pi_m[: m - 1], 0.0)
    assert np.allclose(tau_m[: m - 1], 0.0)


def test_pi_tau_m_stable_high_order():
    """No overflow/NaN at n_max ~ 267 (x=236 detection regime, PHYSICS.md §7)."""
    n_max = 267
    theta = np.linspace(1e-3, np.pi - 1e-3, 500)
    for m in (1, 2):
        pi_m, tau_m = vswf.pi_tau_m(n_max, m, theta)
        assert np.isfinite(pi_m).all()
        assert np.isfinite(tau_m).all()


# ---------------------------------------------------------------------------
# Packing helpers
# ---------------------------------------------------------------------------

def test_m_index_and_values():
    n_max = 6
    ms = vswf.m_values(n_max)
    assert ms[0] == -n_max and ms[-1] == n_max
    assert len(ms) == 2 * n_max + 1
    assert vswf.m_index(0, n_max) == n_max
    assert vswf.m_index(1, n_max) == n_max + 1
    assert vswf.m_index(-1, n_max) == n_max - 1
    assert ms[vswf.m_index(3, n_max)] == 3


def test_plane_wave_bsc_only_pm1_nonzero():
    """Plane-wave BSCs live only on m = +-1 (PHYSICS.md §1.3)."""
    n_max = 10
    g_tm, g_te = vswf.plane_wave_bsc(n_max)
    assert g_tm.shape == (n_max, 2 * n_max + 1)
    ip, im = vswf.m_index(1, n_max), vswf.m_index(-1, n_max)
    mask = np.ones(2 * n_max + 1, dtype=bool)
    mask[[ip, im]] = False
    assert np.allclose(g_tm[:, mask], 0.0)
    assert np.allclose(g_te[:, mask], 0.0)
    # symmetry: g^TM even in m, g^TE odd in m
    assert np.allclose(g_tm[:, ip], g_tm[:, im])
    assert np.allclose(g_te[:, ip], -g_te[:, im])


# ---------------------------------------------------------------------------
# THE ANCHOR: plane-wave far field reproduces S1/S2 to <= 1e-8
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("x", [1.0, 15.0, 60.0])
def test_plane_wave_reconstructs_s1_s2(x):
    """E_theta = S2 cos(phi), E_phi = -S1 sin(phi) to <= 1e-8 (correctness anchor)."""
    n_max = nmax_wiscombe(x)
    grid = AngularGrid.full_sphere(ntheta=64, nphi=16)
    vf = vswf.plane_wave_far_field(M_SILICA, x, grid, n_max)

    S1, S2 = pw.scattering_amplitudes(M_SILICA, x, grid.theta, n_max)  # (Ntheta,)
    cphi = np.cos(grid.phi)[None, :]
    sphi = np.sin(grid.phi)[None, :]
    e_theta_ref = S2[:, None] * cphi
    e_phi_ref = -S1[:, None] * sphi

    err_theta = np.max(np.abs(vf.E_theta - e_theta_ref))
    err_phi = np.max(np.abs(vf.E_phi - e_phi_ref))
    scale = max(np.max(np.abs(S1)), np.max(np.abs(S2)), 1.0)
    assert err_theta / scale <= 1e-8, f"E_theta rel err {err_theta/scale:.2e} at x={x}"
    assert err_phi / scale <= 1e-8, f"E_phi rel err {err_phi/scale:.2e} at x={x}"


@pytest.mark.parametrize("x", [35.0, 95.0, 236.0])
def test_plane_wave_reconstructs_s1_s2_large_x(x):
    """Detection regime x in {35, 95, 236} => n_max up to ~267 (W1 brief)."""
    n_max = nmax_wiscombe(x)
    # forward lobe width ~1/x; use a finer theta grid at large x (PHYSICS.md §7).
    ntheta = int(max(200, 6 * n_max))
    grid = AngularGrid.full_sphere(ntheta=ntheta, nphi=12)
    vf = vswf.plane_wave_far_field(M_SILICA, x, grid, n_max)

    S1, S2 = pw.scattering_amplitudes(M_SILICA, x, grid.theta, n_max)
    cphi = np.cos(grid.phi)[None, :]
    sphi = np.sin(grid.phi)[None, :]
    err_theta = np.max(np.abs(vf.E_theta - S2[:, None] * cphi))
    err_phi = np.max(np.abs(vf.E_phi - (-S1[:, None] * sphi)))
    scale = max(np.max(np.abs(S1)), np.max(np.abs(S2)), 1.0)
    assert err_theta / scale <= 1e-8, f"E_theta rel err {err_theta/scale:.2e} at x={x}"
    assert err_phi / scale <= 1e-8, f"E_phi rel err {err_phi/scale:.2e} at x={x}"


def test_far_field_intensity_matches_physics_formula():
    """|E|^2 = |S2|^2 cos^2(phi) + |S1|^2 sin^2(phi) (PHYSICS.md §1.3)."""
    x = 12.0
    n_max = nmax_wiscombe(x)
    grid = AngularGrid.full_sphere(ntheta=48, nphi=20)
    vf = vswf.plane_wave_far_field(M_SILICA, x, grid, n_max)

    S1, S2 = pw.scattering_amplitudes(M_SILICA, x, grid.theta, n_max)
    cphi2 = np.cos(grid.phi)[None, :] ** 2
    sphi2 = np.sin(grid.phi)[None, :] ** 2
    I_ref = np.abs(S2)[:, None] ** 2 * cphi2 + np.abs(S1)[:, None] ** 2 * sphi2
    scale = np.max(I_ref)
    assert np.max(np.abs(vf.intensity() - I_ref)) / scale <= 1e-8


def test_far_field_absorbing_sphere():
    """Complex m (absorbing) also reconstructs S1/S2 (Im(m) >= 0, e^{-iωt})."""
    m = 1.4607 + 0.01j
    x = 8.0
    n_max = nmax_wiscombe(x)
    grid = AngularGrid.full_sphere(ntheta=40, nphi=12)
    vf = vswf.plane_wave_far_field(m, x, grid, n_max)
    S1, S2 = pw.scattering_amplitudes(m, x, grid.theta, n_max)
    cphi = np.cos(grid.phi)[None, :]
    sphi = np.sin(grid.phi)[None, :]
    scale = max(np.max(np.abs(S1)), np.max(np.abs(S2)), 1.0)
    assert np.max(np.abs(vf.E_theta - S2[:, None] * cphi)) / scale <= 1e-8
    assert np.max(np.abs(vf.E_phi + S1[:, None] * sphi)) / scale <= 1e-8


# ---------------------------------------------------------------------------
# far_field assembler: general (n, m) path and error handling
# ---------------------------------------------------------------------------

def test_far_field_returns_vectorfield_on_same_grid():
    n_max = 6
    grid = AngularGrid.full_sphere(ntheta=20, nphi=8)
    g_tm, g_te = vswf.plane_wave_bsc(n_max)
    a_n, b_n = pw.mie_coefficients(M_SILICA, 3.0, n_max)
    vf = vswf.far_field(a_n[:, None] * g_tm, b_n[:, None] * g_te, grid)
    assert vf.grid is grid
    assert vf.E_theta.shape == (grid.theta.shape[0], grid.phi.shape[0])
    assert vf.E_phi.shape == vf.E_theta.shape


def test_far_field_zero_coeffs_give_zero_field():
    n_max = 5
    grid = AngularGrid.full_sphere(ntheta=16, nphi=6)
    z = np.zeros((n_max, 2 * n_max + 1), dtype=complex)
    vf = vswf.far_field(z, z, grid)
    assert np.allclose(vf.E_theta, 0.0)
    assert np.allclose(vf.E_phi, 0.0)


def test_far_field_bad_shape_raises():
    n_max = 5
    grid = AngularGrid.full_sphere(ntheta=8, nphi=4)
    bad = np.zeros((n_max, n_max), dtype=complex)
    good = np.zeros((n_max, 2 * n_max + 1), dtype=complex)
    with pytest.raises(ValueError):
        vswf.far_field(bad, good, grid)


def test_far_field_convergence_under_nmax_bump():
    """Field stable to <=1e-4 under n_max -> n_max + 8 (G-CONV, VALIDATION.md)."""
    x = 15.0
    grid = AngularGrid.full_sphere(ntheta=64, nphi=12)
    n0 = nmax_wiscombe(x)
    vf0 = vswf.plane_wave_far_field(M_SILICA, x, grid, n0)
    vf1 = vswf.plane_wave_far_field(M_SILICA, x, grid, n0 + 8)
    scale = np.max(np.abs(vf0.E_theta)) + np.max(np.abs(vf0.E_phi))
    d = np.max(np.abs(vf1.E_theta - vf0.E_theta)) + np.max(np.abs(vf1.E_phi - vf0.E_phi))
    assert d / scale <= 1e-4
