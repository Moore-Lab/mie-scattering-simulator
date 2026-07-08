"""W3 data-file load tests (LITERATURE.md sec 6, INTERFACES.md sec 8).

Loads data/literature/experiments.json and benchmarks.json through the frozen
schema loaders and asserts the M3 acceptance counts and the required benchmark
cases. Time convention e^{-iomega t}; SI units. These tests validate the DATA,
not the (orchestrator-owned) schema.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from mieinfo.literature.schema import (
    Benchmark,
    Experiment,
    load_benchmarks,
    load_experiments,
)

# Repo root = .../mie_project ; this file is tests/literature/test_data_loads.py
REPO_ROOT = Path(__file__).resolve().parents[2]
EXPERIMENTS_PATH = REPO_ROOT / "data" / "literature" / "experiments.json"
BENCHMARKS_PATH = REPO_ROOT / "data" / "literature" / "benchmarks.json"

# Controlled vocabularies (LITERATURE.md sec 3). None is always allowed.
VALID_DIRECTIONS = {"forward", "backward", "split", "cavity", "both", "imaging"}
VALID_SCHEMES = {"homodyne", "heterodyne", "self-homodyne", "split", "imaging"}


@pytest.fixture(scope="module")
def experiments() -> list[Experiment]:
    return load_experiments(str(EXPERIMENTS_PATH))


@pytest.fixture(scope="module")
def benchmarks() -> list[Benchmark]:
    return load_benchmarks(str(BENCHMARKS_PATH))


# ---------------------------------------------------------------- experiments


def test_experiments_file_exists():
    assert EXPERIMENTS_PATH.is_file(), f"missing {EXPERIMENTS_PATH}"


def test_experiments_load_and_count(experiments):
    # M3 gate: >= 15 Experiment records.
    assert len(experiments) >= 15, f"need >=15 experiments, got {len(experiments)}"
    assert all(isinstance(e, Experiment) for e in experiments)


def test_experiment_keys_unique(experiments):
    keys = [e.key for e in experiments]
    dupes = {k for k in keys if keys.count(k) > 1}
    assert not dupes, f"duplicate experiment keys: {sorted(dupes)}"


def test_experiments_have_provenance(experiments):
    # Provenance is mandatory for every record (LITERATURE.md sec 3).
    for e in experiments:
        assert e.provenance and e.provenance.strip(), f"{e.key} missing provenance"
        assert e.reference and e.reference.strip(), f"{e.key} missing reference"


def test_refractive_index_convention(experiments):
    # Stored as [re, im]; parsed to complex; Im(m) >= 0 for absorption (PHYSICS.md sec 0).
    for e in experiments:
        if e.refractive_index is not None:
            assert isinstance(e.refractive_index, complex)
            assert e.refractive_index.imag >= 0.0, (
                f"{e.key}: Im(m) < 0 violates e^{{-iwt}} absorption convention"
            )
            assert e.refractive_index.real > 0.0, f"{e.key}: non-physical Re(m)"


def test_controlled_vocabularies(experiments):
    for e in experiments:
        if e.collection_direction is not None:
            assert e.collection_direction in VALID_DIRECTIONS, (
                f"{e.key}: bad collection_direction {e.collection_direction!r}"
            )
        if e.detection_scheme is not None:
            assert e.detection_scheme in VALID_SCHEMES, (
                f"{e.key}: bad detection_scheme {e.detection_scheme!r}"
            )


def test_reported_efficiencies_in_unit_range(experiments):
    # A detection efficiency is a fraction in [0, 1].
    for e in experiments:
        eff = e.reported_detection_efficiency
        if eff is not None:
            assert 0.0 <= eff <= 1.0, f"{e.key}: efficiency {eff} out of [0,1]"


def test_positive_physical_scalars(experiments):
    for e in experiments:
        if e.sphere_radius_m is not None:
            assert e.sphere_radius_m > 0.0, f"{e.key}: non-positive radius"
        if e.wavelength_m is not None:
            assert e.wavelength_m > 0.0, f"{e.key}: non-positive wavelength"
        if e.collection_NA is not None:
            assert 0.0 < e.collection_NA <= 1.0, f"{e.key}: NA out of (0,1]"


def test_tebbenjohanns_dipole_anchor_present(experiments):
    # The primary theory anchor must be in the DB (LITERATURE.md sec 2).
    keys = {e.key for e in experiments}
    assert "tebbenjohanns2019_backward_z_NA08" in keys
    assert "maurer2023_smallparticle_z_backward" in keys


# ----------------------------------------------------------------- benchmarks


def test_benchmarks_file_exists():
    assert BENCHMARKS_PATH.is_file(), f"missing {BENCHMARKS_PATH}"


def test_benchmarks_load_and_count(benchmarks):
    # M3 gate: >= 4 Benchmark records.
    assert len(benchmarks) >= 4, f"need >=4 benchmarks, got {len(benchmarks)}"
    assert all(isinstance(b, Benchmark) for b in benchmarks)


def test_benchmarks_reference_existing_experiments(benchmarks, experiments):
    exp_keys = {e.key for e in experiments}
    for b in benchmarks:
        assert b.experiment_key in exp_keys, (
            f"benchmark references unknown experiment_key {b.experiment_key!r}"
        )


def test_benchmarks_have_target_fields(benchmarks):
    for b in benchmarks:
        assert b.target_quantity and b.target_quantity.strip()
        assert isinstance(b.target_value, float)
        assert b.target_tolerance > 0.0, (
            f"{b.experiment_key}: tolerance must be > 0"
        )
        assert b.target_provenance and b.target_provenance.strip(), (
            f"{b.experiment_key}: missing target_provenance"
        )


def test_benchmark_run_configs_fully_specified(benchmarks):
    # run_config must carry enough to build a mieinfo RunConfig (INTERFACES.md sec 8,
    # ARCHITECTURE.md sec 5): sphere, medium, beam, collection, n_hat.
    required = {"sphere", "medium", "beam", "collection", "n_hat"}
    for b in benchmarks:
        missing = required - set(b.run_config)
        assert not missing, (
            f"{b.experiment_key}: run_config missing keys {sorted(missing)}"
        )
        sphere = b.run_config["sphere"]
        assert "radius_m" in sphere and "m" in sphere
        medium = b.run_config["medium"]
        assert "wavelength_vacuum_m" in medium


def test_required_tebbenjohanns_dipole_benchmark_present(benchmarks):
    # (a) The Tebbenjohanns-2019 dipole case is a hard requirement (LITERATURE.md sec 5,
    #     VALIDATION.md sec 3 G-LIMIT dipole gate).
    by_key = {b.experiment_key: b for b in benchmarks}
    assert "tebbenjohanns2019_backward_z_NA08" in by_key, (
        "the Tebbenjohanns-2019 dipole benchmark must be present"
    )
    dip = by_key["tebbenjohanns2019_backward_z_NA08"]
    assert dip.target_quantity == "detection_efficiency"
    assert dip.target_value == pytest.approx(0.6, abs=1e-9)
    assert "Tebbenjohanns" in dip.target_provenance


def test_required_mie_or_experimental_efficiency_benchmark_present(benchmarks):
    # (b) >= 1 Mie-regime OR experimental efficiency benchmark (LITERATURE.md sec 5).
    keys = {b.experiment_key for b in benchmarks}
    mie_or_experimental = {
        "maurer2023_smallparticle_z_backward",  # Mie-theory (Lorenz-Mie) efficiency
        "maurer2023_gs_threshold",              # Mie-theory ground-state threshold
        "wang2024_mie_backward_highNA",         # Mie-regime efficiency
        "magrini2021_forward_z_confocal",       # measured experimental efficiency
        "tebbenjohanns2021_backward_z_cryo",    # measured experimental efficiency
    }
    assert keys & mie_or_experimental, (
        "need >=1 Mie-regime or experimental-efficiency benchmark"
    )


def test_two_independent_results_covered(benchmarks):
    # M4 success criterion: the benchmark set spans the dipole case AND an
    # independent Mie/experimental efficiency (MASTER_PLAN.md sec 7).
    keys = {b.experiment_key for b in benchmarks}
    assert "tebbenjohanns2019_backward_z_NA08" in keys
    assert keys & {
        "maurer2023_smallparticle_z_backward",
        "wang2024_mie_backward_highNA",
        "magrini2021_forward_z_confocal",
    }
