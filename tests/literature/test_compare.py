"""G-LIT: simulation-vs-literature comparison (VALIDATION.md §6, INTERFACES.md §9).

Loads data/literature/benchmarks.json and runs ``compare_benchmark`` on every
benchmark. Hard success criteria (VALIDATION.md §6, M4):

  * the Tebbenjohanns-2019 dipole case (``tebbenjohanns2019_backward_z_NA08``)
    reproduces within tolerance — the analytic (1−cosθ)² backward-weighting of the
    axial-position information in the dipole limit;
  * at least one Mie-regime OR experimental efficiency benchmark reproduces
    numerically (here the Maurer-2023 small-particle backward-z case);
  * every benchmark RUNS and yields a documented ``discrepancy_note``, so the
    benchmarks that legitimately do NOT match numerically (definition gap /
    missing focused-beam physics) are explained, not silently failed.

"A pass by redefining the quantity is not a pass" (VALIDATION.md §6): the tests
assert the definition each benchmark uses, not just the number. Time convention
e^{-iωt}; SI units.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from mieinfo.detect.optics import CollectionGeometry, collection_efficiency
from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider
from mieinfo.info.modes import information_pattern
from mieinfo.literature.compare import ComparisonResult, compare_benchmark
from mieinfo.literature.schema import Benchmark, load_benchmarks
from mieinfo.types import AngularGrid, Medium, Sphere

REPO_ROOT = Path(__file__).resolve().parents[2]
BENCHMARKS_PATH = REPO_ROOT / "data" / "literature" / "benchmarks.json"


@pytest.fixture(scope="module")
def benchmarks() -> list[Benchmark]:
    return load_benchmarks(str(BENCHMARKS_PATH))


@pytest.fixture(scope="module")
def provider() -> PlaneWaveProvider:
    return PlaneWaveProvider()


@pytest.fixture(scope="module")
def results(benchmarks, provider) -> dict[str, ComparisonResult]:
    return {b.experiment_key: compare_benchmark(provider, b) for b in benchmarks}


# ---------------------------------------------------------------- runs & shape


def test_every_benchmark_runs_and_is_documented(results, benchmarks):
    """compare_benchmark returns a well-formed ComparisonResult with a non-trivial
    discrepancy_note for EVERY benchmark (matched or not)."""
    assert set(results) == {b.experiment_key for b in benchmarks}
    for b in benchmarks:
        r = results[b.experiment_key]
        assert isinstance(r, ComparisonResult)
        assert r.benchmark_key == b.experiment_key
        assert isinstance(r.predicted, float) and np.isfinite(r.predicted)
        assert isinstance(r.reported, float)
        assert r.reported == pytest.approx(b.target_value)
        assert isinstance(r.within_tolerance, bool)
        # A documented note is mandatory (VALIDATION.md §6).
        assert r.discrepancy_note and len(r.discrepancy_note.strip()) > 40
        # Prediction is an efficiency/fraction in [0, 1].
        assert -1e-9 <= r.predicted <= 1.0 + 1e-9


def test_notes_document_beam_substitution(results):
    """Every focused-beam run_config is run with a PlaneWave substitution, which
    must be disclosed in the note (not a silent modeling change)."""
    for key, r in results.items():
        assert "substitution" in r.discrepancy_note.lower(), key
        # And the note records the like-for-like efficiency definition claim.
        assert "η_q" in r.discrepancy_note or "eta" in r.discrepancy_note.lower(), key


# ------------------------------------------------- HARD criterion: dipole case


def test_tebbenjohanns_dipole_reproduces_within_tolerance(results):
    """VALIDATION.md §6 / M4: the Tebbenjohanns-2019 dipole benchmark
    (η_bw_z = 0.6 at NA=0.8) must reproduce within its stated tolerance."""
    r = results["tebbenjohanns2019_backward_z_NA08"]
    assert r.within_tolerance, r.discrepancy_note
    assert abs(r.predicted - r.reported) <= 0.1  # the benchmark's target_tolerance
    # Sanity: the backward-collected z-info dominates (backward-weighting).
    assert r.predicted > 0.4


def test_dipole_axial_info_has_analytic_one_minus_cos_squared_weighting():
    """The physics behind the pass: in the dipole limit the z-position information
    density equals |E_s|² · (1 − cosθ)² to machine precision — the analytic
    backward-weighting from the phase-gradient factor (k̂ − ŝ)_z = 1 − cosθ
    (PHYSICS.md §3.1, §4.1). This is WHY η_bw_z ≈ 0.6 and the backward
    hemisphere holds ~90% of the z-information."""
    med = Medium(n=1.0, wavelength_vacuum_m=1.55e-6)
    a = 0.01 / med.k                       # x = 0.01, deep Rayleigh/dipole limit
    sph = Sphere(radius_m=a, m=1.44 + 0j)
    beam = PlaneWave(med)
    grid = AngularGrid.full_sphere(300, 8)
    prov = PlaneWaveProvider()

    pz = information_pattern(prov, grid, sph, beam, np.zeros(3), np.array([0.0, 0, 1.0]))
    intensity = prov.field(grid, sph, beam, np.zeros(3)).intensity()
    weight = (1.0 - np.cos(grid.theta))[:, None] ** 2  # (Ntheta, 1) -> broadcast

    ratio = pz.density / (intensity * weight)
    # constant to machine precision -> density == intensity * (1 - cosθ)²
    assert np.std(ratio) / np.mean(ratio) < 1e-9

    # And that weighting makes the backward hemisphere hold >90% of the z-info.
    bwd = collection_efficiency(pz, CollectionGeometry(direction="backward", NA=1.0))
    assert bwd > 0.9


# ---------------------------- HARD criterion: a Mie/experimental efficiency pass


def test_at_least_two_benchmarks_pass_numerically(results):
    """M4: >= 2 independent results reproduced — the dipole case AND a Mie-regime
    or experimental efficiency benchmark (VALIDATION.md §6, §8)."""
    passed = {k for k, r in results.items() if r.within_tolerance}
    assert "tebbenjohanns2019_backward_z_NA08" in passed
    mie_or_experimental = {
        "maurer2023_smallparticle_z_backward",
        "maurer2023_gs_threshold",
        "magrini2021_forward_z_confocal",
    }
    assert passed & mie_or_experimental, (
        f"need a Mie/experimental pass beyond the dipole case; passed={sorted(passed)}"
    )
    assert len(passed) >= 2


def test_maurer_small_particle_backward_z_reproduces(results):
    """The Maurer-2023 small-particle (Rayleigh) backward-z detection efficiency
    η^d_z = 0.55 (Eq. 24 definition) reproduces within tolerance — an independent
    Mie-theory efficiency benchmark, same definition as mieinfo's η_q."""
    r = results["maurer2023_smallparticle_z_backward"]
    assert r.within_tolerance, r.discrepancy_note
    assert abs(r.predicted - r.reported) <= 0.08


def test_maurer_ground_state_threshold_is_cleared_one_sided(results):
    """The 1/9 ground-state-cooling threshold is a one-sided necessary condition
    (η^d > 1/9); the prediction must CLEAR it, and the note must say so."""
    r = results["maurer2023_gs_threshold"]
    assert r.within_tolerance, r.discrepancy_note
    assert r.predicted >= 1.0 / 9.0
    note = r.discrepancy_note.lower()
    assert "threshold" in note or "one-sided" in note or "necessary condition" in note


def test_tebbenjohanns_ideal_backward_fraction_one_sided(results):
    """The '>90% of the z-information is in the backward half-space' claim is a
    one-sided lower bound; predicted backward fraction must be >= 0.9 - tol, and
    the definition used is the backward HEMISPHERE (not the run_config NA)."""
    r = results["tebbenjohanns2019_ideal_z"]
    assert r.within_tolerance, r.discrepancy_note
    assert r.predicted >= 0.9 - 0.1
    # Cross-check against a direct backward-hemisphere computation of the fraction.
    med = Medium(n=1.0, wavelength_vacuum_m=1.55e-6)
    sph = Sphere(radius_m=5.0e-8, m=1.44 + 0j)
    grid = AngularGrid.full_sphere(160, 320)
    prov = PlaneWaveProvider()
    pz = information_pattern(prov, grid, sph, PlaneWave(med), np.zeros(3),
                             np.array([0.0, 0, 1.0]))
    direct = collection_efficiency(pz, CollectionGeometry(direction="backward", NA=1.0))
    assert r.predicted == pytest.approx(direct, abs=2e-3)


# ------------------------------- explained-not-matched: definition/model gaps


def test_magrini_is_explained_not_matched(results):
    """Magrini's reported η = Γ_meas/(Γ_ba+Γ_th) is a TOTAL efficiency (losses+QE);
    the plane-wave substitution also drops the focused beam's longitudinal E_z field
    that carries forward axial sensitivity. So this benchmark does NOT match
    numerically, and that must be DOCUMENTED (not silently passed)."""
    r = results["magrini2021_forward_z_confocal"]
    assert not r.within_tolerance
    note = r.discrepancy_note.lower()
    assert "explained-not-matched" in note
    # The note must name the two causes: total-vs-ideal definition gap and the
    # missing focused-beam longitudinal field.
    assert "total" in note and ("longitudinal" in note or "e_z" in note)


def test_no_pass_by_redefinition(results, benchmarks):
    """VALIDATION.md §6: a pass by redefining the quantity is not a pass. Every
    two-sided efficiency benchmark that passes must do so under the SAME collection
    geometry (direction + NA) the run_config specifies — not a relaxed cone."""
    by_key = {b.experiment_key: b for b in benchmarks}
    for key, r in results.items():
        if not r.within_tolerance:
            continue
        b = by_key[key]
        q = b.target_quantity.lower()
        # Threshold / lower-bound / hemisphere-fraction cases are legitimately
        # one-sided (documented); skip the strict-geometry check for those.
        if q == "backward_info_fraction":
            continue
        prov_txt = b.target_provenance.lower()
        if "threshold" in prov_txt or "necessary condition" in prov_txt:
            continue
        # Two-sided efficiency pass: recompute with the run_config's exact geometry
        # and confirm it is what compare_benchmark used.
        cfg = b.run_config
        med = Medium(n=1.0, wavelength_vacuum_m=float(cfg["medium"]["wavelength_vacuum_m"]))
        m = cfg["sphere"]["m"]
        m = complex(m[0], m[1]) if isinstance(m, (list, tuple)) else complex(m)
        sph = Sphere(radius_m=float(cfg["sphere"]["radius_m"]), m=m)
        g = cfg["grid"]
        grid = AngularGrid.full_sphere(int(g["ntheta"]), int(g["nphi"]))
        n_hat = np.asarray(cfg["n_hat"], float)
        pattern = information_pattern(PlaneWaveProvider(), grid, sph, PlaneWave(med),
                                      np.zeros(3), n_hat)
        geom = CollectionGeometry(direction=cfg["collection"]["direction"],
                                  NA=float(cfg["collection"]["NA"]))
        recomputed = collection_efficiency(pattern, geom)
        assert recomputed == pytest.approx(r.predicted, abs=2e-3), (
            f"{key}: passed prediction does not match the run_config geometry "
            f"(possible redefinition)"
        )
