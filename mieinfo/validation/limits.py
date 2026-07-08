"""G-LIMIT — analytic-limit gates (VALIDATION.md §3, PHYSICS.md §4.1, §6).
Time convention e^{-iωt}.

Three families of reusable checks on the current plane-wave engine:

1. Energy / optical theorem (`check_energy_optical_theorem`):
   for real m, Q_ext == Q_sca to <= 1e-10 (non-absorbing spheres scatter all they
   extinguish); for absorbing m (Im(m) > 0), Q_ext > Q_sca.

2. Dipole-limit information-pattern structure (`check_dipole_information_structure`):
   at x << 1 the scattered field -> a dipole field and the information radiation pattern
   -> the Tebbenjohanns-2019 dipole result (PHYSICS.md §4.1). Structural facts:
   - transverse (x) displacement info is forward-heavy relative to axial and peaks
     off-axis (the sin^2(theta) cos^2(phi) weight on a fore-aft-symmetric dipole gives a
     forward fraction >> the axial one and a peak near ~50 deg);
   - axial (z) displacement info is backward-weighted: the (1 - cos theta)^2 weight puts
     the peak in backscatter (theta > 90 deg) and >~ 85% of the info behind the sphere.
   This is the "backscatter-optimal axial detection" statement of Tebbenjohanns 2019.

3. Reciprocity / symmetry (`check_symmetry_phi_parity`, `check_xy_cross_identity`):
   for x-polarized incidence the intensity and information patterns are EVEN under both
   phi -> -phi and phi -> pi - phi (the x-polarization mirror planes). The y-displacement
   pattern is NOT a plain phi -> phi - pi/2 relabel of the x-displacement pattern, because
   the shared intensity factor I = |S2|^2 cos^2(phi) + |S1|^2 sin^2(phi) is itself not
   pi/2-periodic when S1 != S2. The EXACT relation that does hold (to machine precision) is

       density_y(theta, phi) * cos^2(phi) == density_x(theta, phi) * sin^2(phi),

   because both information densities share the common intensity factor and differ only by
   the transverse-displacement weight (cos^2(phi) for x, sin^2(phi) for y):
       density_x = sin^2(theta) cos^2(phi) * I,   density_y = sin^2(theta) sin^2(phi) * I.
   Cross-multiplying by the opposite weight cancels I and leaves the identity above. This
   is verified and documented per the directive.

Each check returns a plain result (bool, plus the observed margin where meaningful).
"""
from __future__ import annotations

import numpy as np

from ..glmt.beam import PlaneWave
from ..glmt.scatter import PlaneWaveProvider
from ..info.modes import InformationPattern, information_pattern
from ..mie import plane_wave as pw
from ..types import AngularGrid, Medium, Sphere

# Contract tolerances (VALIDATION.md §3).
TOL_OPTICAL_THEOREM = 1e-10  # real m: |Q_ext - Q_sca| / Q_sca
TOL_PARITY = 1e-9            # phi-parity relative residual (machine precision in practice)
TOL_CROSS_IDENTITY = 1e-9    # x/y cross identity relative residual

_MED_532 = Medium(n=1.0, wavelength_vacuum_m=532e-9)
_M_SILICA = 1.46 + 0j


# --------------------------------------------------------------------------- #
# 1. Energy / optical theorem
# --------------------------------------------------------------------------- #
def check_energy_optical_theorem(
    real_cases: tuple[tuple[complex, float], ...] = (
        (1.46 + 0j, 3.0), (1.46 + 0j, 15.0), (1.46 + 0j, 30.0),
        (1.33 + 0j, 3.0), (1.45 + 0j, 5.9),
    ),
    absorbing_cases: tuple[tuple[complex, float], ...] = (
        (1.5 + 0.01j, 10.0), (1.46 + 0.1j, 15.0), (1.6 + 0.5j, 5.0),
    ),
) -> tuple[bool, float]:
    """Real m => Q_ext == Q_sca (<= 1e-10); absorbing m => Q_ext > Q_sca.

    Returns (passed, max_rel_diff) where max_rel_diff is the worst |Q_ext - Q_sca|/Q_sca
    over the real-m cases (the quantity gated at 1e-10).
    """
    passed = True
    max_rel = 0.0
    for m, x in real_cases:
        q_ext, q_sca, _, _ = pw.efficiencies(m, x)
        rel = abs(q_ext - q_sca) / q_sca
        max_rel = max(max_rel, rel)
        if rel > TOL_OPTICAL_THEOREM:
            passed = False
    for m, x in absorbing_cases:
        q_ext, q_sca, _, _ = pw.efficiencies(m, x)
        if not (q_ext > q_sca):
            passed = False
    return passed, max_rel


# --------------------------------------------------------------------------- #
# 2. Dipole-limit information-pattern structure
# --------------------------------------------------------------------------- #
def _dipole_patterns(x: float = 0.05, ntheta: int = 200, nphi: int = 8):
    sphere = Sphere(radius_m=x / _MED_532.k, m=_M_SILICA)
    grid = AngularGrid.full_sphere(ntheta, nphi)
    beam = PlaneWave(_MED_532)
    provider = PlaneWaveProvider()
    px = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([1.0, 0.0, 0.0]))
    pz = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([0.0, 0.0, 1.0]))
    return px, pz


def _forward_fraction(pattern: InformationPattern) -> float:
    fwd = pattern.grid.theta < np.pi / 2
    w = pattern.grid.w_solid
    return float(np.sum(pattern.density[fwd] * w[fwd]) / np.sum(pattern.density * w))


def _peak_deg(pattern: InformationPattern) -> float:
    profile = np.sum(pattern.density * pattern.grid.w_solid, axis=1)
    return float(np.degrees(pattern.grid.theta[np.argmax(profile)]))


def check_dipole_information_structure(x: float = 0.05) -> tuple[bool, dict]:
    """Dipole-limit (x << 1) information-pattern structure vs Tebbenjohanns 2019.

    Returns (passed, metrics) with the transverse/axial forward fractions and peak
    angles used for the decision. Asserts:
    - transverse info is far more forward than axial: x_fwd_frac >> z_fwd_frac;
    - axial info is backward-weighted: z peak in backscatter (> 90 deg) and z forward
      fraction small (< ~0.15, i.e. axial backward fraction > ~0.85);
    - transverse peak is off-axis (not forward, not backward): 20 deg < x_peak < 90 deg.
    """
    px, pz = _dipole_patterns(x)
    metrics = {
        "x_fwd_frac": _forward_fraction(px),
        "z_fwd_frac": _forward_fraction(pz),
        "x_peak_deg": _peak_deg(px),
        "z_peak_deg": _peak_deg(pz),
    }
    passed = (
        metrics["x_fwd_frac"] > metrics["z_fwd_frac"] + 0.3   # transverse >> axial forward
        and metrics["z_fwd_frac"] < 0.15                      # axial backward-weighted
        and metrics["z_peak_deg"] > 90.0                      # axial peak in backscatter
        and 20.0 < metrics["x_peak_deg"] < 90.0               # transverse peak off-axis
    )
    return passed, metrics


# --------------------------------------------------------------------------- #
# 3. Reciprocity / symmetry
# --------------------------------------------------------------------------- #
def _symmetry_patterns(x: float = 6.0, ntheta: int = 40, nphi: int = 72):
    """Grid with an even Nphi divisible by 4 so both phi->-phi and phi->pi-phi index
    reflections (and the pi/2 quarter-turn used elsewhere) land exactly on nodes."""
    assert nphi % 4 == 0, "Nphi must be divisible by 4 for exact phi reflections"
    sphere = Sphere(radius_m=x / _MED_532.k, m=_M_SILICA)
    grid = AngularGrid.full_sphere(ntheta, nphi)
    beam = PlaneWave(_MED_532)
    provider = PlaneWaveProvider()
    px = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([1.0, 0.0, 0.0]))
    py = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([0.0, 1.0, 0.0]))
    pz = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([0.0, 0.0, 1.0]))
    field = provider.field(grid, sphere, beam, np.zeros(3))
    return grid, field, px, py, pz


def _phi_parity_relerr(density: np.ndarray, nphi: int) -> tuple[float, float]:
    """Max relative residual of density under phi -> -phi and phi -> pi - phi.

    Uniform phi with endpoint=False on [0, 2pi): phi_j = 2pi j / Nphi.
    phi -> -phi maps index j -> (-j) mod Nphi; phi -> pi - phi maps j -> (Nphi/2 - j) mod Nphi.
    """
    idx = np.arange(nphi)
    idx_neg = (-idx) % nphi
    idx_pi_minus = (nphi // 2 - idx) % nphi
    scale = np.max(density)
    e_neg = np.max(np.abs(density - density[:, idx_neg])) / scale
    e_pim = np.max(np.abs(density - density[:, idx_pi_minus])) / scale
    return e_neg, e_pim


def check_symmetry_phi_parity(x: float = 6.0) -> tuple[bool, float]:
    """phi-parity of the x-pol intensity and x/y/z information patterns.

    For x-polarized incidence every pattern is EVEN under phi -> -phi and phi -> pi - phi
    (the two mirror planes of x-polarization). Returns (passed, max_relerr).
    """
    grid, field, px, py, pz = _symmetry_patterns(x)
    nphi = len(grid.phi)
    max_relerr = 0.0
    for density in (field.intensity(), px.density, py.density, pz.density):
        e_neg, e_pim = _phi_parity_relerr(density, nphi)
        max_relerr = max(max_relerr, e_neg, e_pim)
    return (max_relerr <= TOL_PARITY), max_relerr


def check_xy_cross_identity(x: float = 6.0) -> tuple[bool, float]:
    """Exact x/y information-pattern symmetry (see module docstring):

        density_y(theta, phi) * cos^2(phi) == density_x(theta, phi) * sin^2(phi).

    Both equal sin^2(theta) cos^2(phi) sin^2(phi) * I, so cross-multiplying cancels the
    common intensity factor I and the identity holds to machine precision. This is the
    correct statement of the "y-displacement pattern relates to the x-displacement one by
    a phi -> phi - pi/2 azimuth relabeling" (the plain relabel fails because I is not
    pi/2-periodic; the relabel holds only jointly with an S1<->S2 swap). Returns
    (passed, max_relerr).
    """
    grid, _field, px, py, _pz = _symmetry_patterns(x)
    cos2 = np.cos(grid.phi)[None, :] ** 2
    sin2 = np.sin(grid.phi)[None, :] ** 2
    lhs = py.density * cos2
    rhs = px.density * sin2
    scale = max(np.max(px.density), np.max(py.density))
    max_relerr = float(np.max(np.abs(lhs - rhs)) / scale)
    return (max_relerr <= TOL_CROSS_IDENTITY), max_relerr


# NOTE (documented non-symmetry): the TOTAL Fisher information for x- and y-displacement
# are NOT equal for a Mie sphere. Integrating density_x = sin^2(theta) cos^2(phi) I and
# density_y = sin^2(theta) sin^2(phi) I over phi gives
#     F_x propto |S2|^2 * 3/4 + |S1|^2 * 1/4,   F_y propto |S2|^2 * 1/4 + |S1|^2 * 3/4
# per theta, which differ whenever |S1| != |S2| (i.e. always, except the isotropic
# forward/back points). In the dipole limit F_y -> 2 F_x. So there is NO x<->y total-info
# symmetry; only the pointwise cross identity above holds exactly. This is why the gate is
# the cross identity, not an equality of totals.
