"""G-GOLD — golden regression gates (VALIDATION.md §2). Time convention e^{-iωt}.

Two reusable checks:

- `check_efficiencies_golden` — `mie.plane_wave.efficiencies` reproduces
  `data/golden/golden_values.json` (`Q_ext, Q_sca, Q_back, g`) for every case.
- `check_information_pattern_golden` — the plane-wave information-pattern pipeline
  (`PlaneWaveProvider` → `field_derivative` → `combine_direction` → `info_density`
  → `collection_efficiency`) reproduces `data/golden/information_pattern_results.json`
  (the operative 532 nm set): per-axis forward-info fractions, peak polar angles, and
  η at NA 0.5/0.8/0.95 for x and z, forward and backward.

Both return `(passed: bool, max_rel_err: float)` so the pytest suite and
`mieinfo validate` reuse the same logic. `Q_back` is a coherent alternating sum with a
float64 cancellation floor (~1e-5–1e-6 at large x, PHYSICS.md §1.2), so it carries a
looser tolerance than `Q_ext/Q_sca/g` — that is a precision floor, not an error.

The information-pattern *structure* (forward fraction, peak angle) reproduces the golden
to machine precision; the golden η were read off an ~80-point NA interpolation, so η
carries a small absolute tolerance (the golden's interpolation artifact, not a package
error — the package value is the more accurate one; see tests/info/test_golden_pattern.py).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..detect.optics import CollectionGeometry, collection_efficiency
from ..glmt.beam import PlaneWave
from ..glmt.scatter import PlaneWaveProvider
from ..info.modes import InformationPattern, information_pattern
from ..mie import plane_wave as pw
from ..types import AngularGrid, Medium, Sphere

# Repo-root-relative data paths (this file is mieinfo/validation/golden.py).
_DATA = Path(__file__).resolve().parents[2] / "data" / "golden"
GOLDEN_VALUES = _DATA / "golden_values.json"
INFORMATION_PATTERN = _DATA / "information_pattern_results.json"

# Contract tolerances (VALIDATION.md §2).
TOL_EFF = 1e-9            # Q_ext, Q_sca, g relative
TOL_QBACK = 1e-5          # Q_back: alternating-sum cancellation floor at large x
TOL_STRUCT_FRAC = 1e-9    # forward-info fraction, absolute
TOL_STRUCT_PEAK_DEG = 1e-6  # peak polar angle, degrees
TOL_ETA = 6e-3           # η: golden NA-interpolation artifact (see module docstring)

# Detection golden set is at 532 nm in vacuum; silica m ~ 1.46 there.
_MED_532 = Medium(n=1.0, wavelength_vacuum_m=532e-9)
_M_SILICA_532 = 1.46 + 0j


def load_golden_values() -> list[dict]:
    """The `cases` list from golden_values.json (Q_ext/Q_sca/Q_back/g at fixed x)."""
    return json.loads(GOLDEN_VALUES.read_text())["cases"]


def load_information_pattern_golden() -> list[dict]:
    """The per-radius records from information_pattern_results.json (532 nm set)."""
    return json.loads(INFORMATION_PATTERN.read_text())


def check_efficiencies_golden() -> tuple[bool, float]:
    """`mie.plane_wave.efficiencies` vs golden_values.json.

    Returns (passed, max_rel_err) over all cases and the four quantities. Q_back uses
    the looser cancellation-floor tolerance; the others use TOL_EFF.
    """
    passed = True
    max_rel = 0.0
    for c in load_golden_values():
        m = complex(c["m_real"], c["m_imag"])
        q_ext, q_sca, q_back, g = pw.efficiencies(m, c["x"])
        got = {"Qext": q_ext, "Qsca": q_sca, "Qback": q_back, "g": g}
        for key, val in got.items():
            ref = c[key]
            rel = abs(val - ref) / max(abs(ref), 1e-30)
            max_rel = max(max_rel, rel)
            tol = TOL_QBACK if key == "Qback" else TOL_EFF
            if rel > tol:
                passed = False
    return passed, max_rel


def _grid_for(a_um: float) -> tuple[int, int]:
    """Matching angular grid for the 532 nm golden set (VALIDATION.md §2: 700–1400 θ
    nodes by size; larger x needs a finer forward lobe, PHYSICS.md §7)."""
    ntheta = 1400 if a_um >= 12 else (1000 if a_um >= 8 else 700)
    return ntheta, 420


def _forward_fraction(pattern: InformationPattern) -> float:
    fwd = pattern.grid.theta < np.pi / 2
    w = pattern.grid.w_solid
    return float(np.sum(pattern.density[fwd] * w[fwd]) / np.sum(pattern.density * w))


def _peak_deg(pattern: InformationPattern) -> float:
    profile = np.sum(pattern.density * pattern.grid.w_solid, axis=1)
    return float(np.degrees(pattern.grid.theta[np.argmax(profile)]))


def _patterns_for_case(record: dict) -> tuple[InformationPattern, InformationPattern]:
    a_um = record["a_um"]
    sphere = Sphere(radius_m=a_um * 1e-6, m=_M_SILICA_532)
    ntheta, nphi = _grid_for(a_um)
    grid = AngularGrid.full_sphere(ntheta, nphi)
    beam = PlaneWave(_MED_532)
    provider = PlaneWaveProvider()
    px = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([1.0, 0.0, 0.0]))
    pz = information_pattern(provider, grid, sphere, beam, np.zeros(3),
                             np.array([0.0, 0.0, 1.0]))
    return px, pz


def check_information_pattern_case(record: dict) -> tuple[bool, float]:
    """One radius of information_pattern_results.json. Returns (passed, max_eta_abs_err).

    Structure (forward fraction, peak angle) is asserted to machine precision inside;
    the returned scalar is the max absolute η error (the interpolation-limited quantity).
    """
    px, pz = _patterns_for_case(record)
    passed = True

    # Structure — reproduces to machine precision.
    if abs(_forward_fraction(px) - record["info_x_fwd"]) >= TOL_STRUCT_FRAC:
        passed = False
    if abs(_forward_fraction(pz) - record["info_z_fwd"]) >= TOL_STRUCT_FRAC:
        passed = False
    if abs(_peak_deg(px) - record["x_peak_deg"]) >= TOL_STRUCT_PEAK_DEG:
        passed = False
    if abs(_peak_deg(pz) - record["z_peak_deg"]) >= TOL_STRUCT_PEAK_DEG:
        passed = False

    # η at each NA / direction — within the golden NA-interpolation artifact.
    max_eta_err = 0.0
    for axis, pattern in (("x", px), ("z", pz)):
        for direction, tag in (("forward", "for"), ("backward", "bac")):
            for na in (0.5, 0.8, 0.95):
                eta = collection_efficiency(
                    pattern, CollectionGeometry(direction=direction, NA=na))
                ref = record[f"eta_{axis}_{tag}_NA{na}"]
                err = abs(eta - ref)
                max_eta_err = max(max_eta_err, err)
                if err >= TOL_ETA:
                    passed = False
    return passed, max_eta_err


def check_information_pattern_golden() -> tuple[bool, float]:
    """Full information_pattern_results.json set. Returns (passed, max_eta_abs_err).

    Heavy (matching grids up to 1400x420 over five radii) — the nightly/full lane runs
    it; the pytest wrapper marks it @slow.
    """
    passed = True
    max_eta_err = 0.0
    for record in load_information_pattern_golden():
        ok, err = check_information_pattern_case(record)
        passed = passed and ok
        max_eta_err = max(max_eta_err, err)
    return passed, max_eta_err
