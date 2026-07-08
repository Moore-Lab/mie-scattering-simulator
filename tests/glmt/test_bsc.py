"""Tests for glmt.bsc — beam-shape coefficients by quadrature / localized / angular
spectrum, and the beam classes they consume (glmt.beam GaussianParaxial /
RichardsWolfFocus). Time convention e^{-iωt} (PHYSICS.md §0).

Central gate here: G-LIMIT — a plane wave through any BSC method returns the
plane-wave weights (mie.vswf.plane_wave_bsc) to <= 1e-6, and a wide-waist Gaussian
approaches them (VALIDATION.md §3).
"""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.glmt import bsc as B
from mieinfo.glmt.beam import GaussianParaxial, PlaneWave, RichardsWolfFocus
from mieinfo.mie import vswf
from mieinfo.types import Medium

MED = Medium(n=1.0, wavelength_vacuum_m=532e-9)
K = MED.k
LAM = 2.0 * np.pi / K


# ---------------------------------------------------------------------------
# Beam classes
# ---------------------------------------------------------------------------

def test_gaussian_unit_amplitude_at_focus_and_no_ez():
    g = GaussianParaxial(MED, waist_m=3 * LAM)
    E0 = g.focal_field(np.zeros((1, 3)))
    assert abs(abs(E0[0, 0]) - 1.0) < 1e-12       # |E| = 1 at focus
    assert abs(E0[0, 1]) < 1e-12                   # y-component zero for x-pol
    assert abs(E0[0, 2]) < 1e-12                   # paraxial: no longitudinal E_z
    assert g.waist_m() == pytest.approx(3 * LAM)


def test_gaussian_rejects_bad_waist():
    with pytest.raises(ValueError):
        GaussianParaxial(MED, waist_m=0.0)


def test_richards_wolf_unit_amplitude_and_ez_grows_with_na():
    ez_frac = []
    for NA in (0.1, 0.5, 0.9):
        rw = RichardsWolfFocus(MED, NA=NA, n_quad=150)
        # off-axis point picks up the longitudinal component
        p = np.array([[0.3 * LAM, 0.0, 0.0]])
        E = rw.focal_field(p)[0]
        amp = np.sqrt(np.sum(np.abs(E) ** 2))
        ez_frac.append(abs(E[2]) / amp)
        E0 = rw.focal_field(np.zeros((1, 3)))
        assert abs(np.sqrt(np.sum(np.abs(E0[0]) ** 2)) - 1.0) < 1e-9   # |E(0)| = 1
    # longitudinal content increases with NA (high-NA focus is non-paraxial)
    assert ez_frac[0] < ez_frac[1] < ez_frac[2]


def test_richards_wolf_rejects_bad_na():
    with pytest.raises(ValueError):
        RichardsWolfFocus(MED, NA=1.5)
    with pytest.raises(ValueError):
        RichardsWolfFocus(MED, NA=0.0)


# ---------------------------------------------------------------------------
# G-LIMIT: plane wave -> plane-wave weights (all methods)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n_max", [8, 20, 40])
def test_plane_wave_quadrature_gives_plane_wave_weights(n_max):
    """bsc_quadrature(PlaneWave) reproduces mie.vswf.plane_wave_bsc (G-LIMIT)."""
    r_tm, r_te = vswf.plane_wave_bsc(n_max)
    scale = np.max(np.abs(r_tm))
    b = B.bsc_quadrature(PlaneWave(MED), np.zeros(3), n_max)
    assert b.g_tm.shape == (n_max, 2 * n_max + 1)
    assert np.max(np.abs(b.g_tm - r_tm)) / scale <= 1e-6
    assert np.max(np.abs(b.g_te - r_te)) / scale <= 1e-6


def test_plane_wave_quadrature_only_pm1_nonzero():
    n_max = 10
    b = B.bsc_quadrature(PlaneWave(MED), np.zeros(3), n_max)
    ip, im = vswf.m_index(1, n_max), vswf.m_index(-1, n_max)
    mask = np.ones(2 * n_max + 1, dtype=bool)
    mask[[ip, im]] = False
    assert np.max(np.abs(b.g_tm[:, mask])) <= 1e-9
    assert np.max(np.abs(b.g_te[:, mask])) <= 1e-9


def test_localized_wide_waist_gives_plane_wave_weights():
    """bsc_localized with waist -> inf reproduces the plane-wave weights (G-LIMIT)."""
    n_max = 12
    r_tm, r_te = vswf.plane_wave_bsc(n_max)
    scale = np.max(np.abs(r_tm))
    g = GaussianParaxial(MED, waist_m=1000 * LAM)
    b = B.bsc_localized(g, np.zeros(3), n_max)
    assert np.max(np.abs(b.g_tm - r_tm)) / scale <= 1e-6
    assert np.max(np.abs(b.g_te - r_te)) / scale <= 1e-6


def test_quadrature_wide_waist_gaussian_gives_plane_wave_weights():
    n_max = 12
    r_tm, r_te = vswf.plane_wave_bsc(n_max)
    scale = np.max(np.abs(r_tm))
    g = GaussianParaxial(MED, waist_m=1000 * LAM)
    b = B.bsc_quadrature(g, np.zeros(3), n_max)
    assert np.max(np.abs(b.g_tm - r_tm)) / scale <= 1e-6


def test_localized_matches_quadrature_at_moderate_waist():
    """The fast localized approximation agrees with the slow quadrature reference to
    a few percent at moderate focusing (localized validity, PHYSICS.md §2.2)."""
    n_max = 12
    g = GaussianParaxial(MED, waist_m=8 * LAM)
    bl = B.bsc_localized(g, np.zeros(3), n_max)
    bq = B.bsc_quadrature(g, np.zeros(3), n_max)
    scale = np.max(np.abs(bq.g_tm))
    rel = np.max(np.abs(bl.g_tm - bq.g_tm)) / scale
    assert rel < 5e-2, f"localized vs quadrature rel {rel:.2e}"


def test_localized_off_axis_raises():
    g = GaussianParaxial(MED, waist_m=5 * LAM)
    with pytest.raises(NotImplementedError):
        B.bsc_localized(g, np.array([0.1e-6, 0.0, 0.0]), 8)


# ---------------------------------------------------------------------------
# Convergence & high-n stability (detection regime)
# ---------------------------------------------------------------------------

def test_quadrature_stable_high_nmax():
    """No NaN/overflow at n_max ~ 264 (x=236 detection regime, PHYSICS.md §7)."""
    n_max = 264
    b = B.bsc_quadrature(PlaneWave(MED), np.zeros(3), n_max, m_max=2)
    r_tm, r_te = vswf.plane_wave_bsc(n_max)
    scale = np.max(np.abs(r_tm))
    assert np.isfinite(b.g_tm).all() and np.isfinite(b.g_te).all()
    assert np.max(np.abs(b.g_tm - r_tm)) / scale <= 1e-6


def test_quadrature_converges_under_resolution_bump():
    n_max = 10
    g = GaussianParaxial(MED, waist_m=3 * LAM)
    b0 = B.bsc_quadrature(g, np.zeros(3), n_max)
    b1 = B.bsc_quadrature(g, np.zeros(3), n_max, ntheta=4 * n_max + 40, nphi=16)
    scale = np.max(np.abs(b0.g_tm))
    assert np.max(np.abs(b1.g_tm - b0.g_tm)) / scale <= 1e-4


# ---------------------------------------------------------------------------
# Angular spectrum (Richards-Wolf) BSC
# ---------------------------------------------------------------------------

def test_angular_spectrum_low_na_is_plane_wave_up_to_phase():
    """A very-low-NA RW focus has BSCs proportional to the plane-wave weights on
    m = +-1 (same n-dependence) up to a global phase (PHYSICS.md §2.2)."""
    n_max = 8
    rw = RichardsWolfFocus(MED, NA=0.02, n_quad=250)
    b = B.bsc_angular_spectrum(rw, np.zeros(3), n_max)
    r_tm, _ = vswf.plane_wave_bsc(n_max)
    col = vswf.m_index(1, n_max)
    ratio = b.g_tm[:, col] / r_tm[:, col]
    # |ratio| ~ constant across n (only a global complex phase differs)
    good = np.abs(ratio) > 0
    spread = np.std(np.abs(ratio[good])) / np.mean(np.abs(ratio[good]))
    assert spread < 5e-3, f"|BSC/ref| not constant across n (spread {spread:.2e})"


def test_angular_spectrum_rejects_wrong_beam():
    with pytest.raises(TypeError):
        B.bsc_angular_spectrum(PlaneWave(MED), np.zeros(3), 8)


# ---------------------------------------------------------------------------
# BSC container + scattered_amplitudes helper
# ---------------------------------------------------------------------------

def test_bsc_shape_validation():
    n_max = 5
    good = np.zeros((n_max, 2 * n_max + 1), dtype=complex)
    bad = np.zeros((n_max, n_max), dtype=complex)
    B.BSC(n_max=n_max, g_tm=good, g_te=good.copy())
    with pytest.raises(ValueError):
        B.BSC(n_max=n_max, g_tm=bad, g_te=good)


def test_scattered_amplitudes_multiplies_mie_coeffs():
    from mieinfo.mie import plane_wave as pw
    n_max = 8
    b = B.bsc_quadrature(PlaneWave(MED), np.zeros(3), n_max)
    a_n, b_n = pw.mie_coefficients(1.4607, 12.0, n_max)
    a_tm, b_te = B.scattered_amplitudes(b, a_n, b_n)
    assert np.allclose(a_tm, a_n[:, None] * b.g_tm)
    assert np.allclose(b_te, b_n[:, None] * b.g_te)
    with pytest.raises(ValueError):
        B.scattered_amplitudes(b, a_n[:-1], b_n)
