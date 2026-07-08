"""Smoke tests for mieinfo.viz.plots (headless Agg; no display).

Each plotting function is called with a small grid and a tmp_path; we assert the
file is written and non-empty. Runtime is kept tiny (coarse grids, few NA points).
The comparison figure is fed a duck-typed ComparisonResult (INTERFACES.md §9) since
W4c's compare_benchmark is not part of M0.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider
from mieinfo.info.modes import information_pattern
from mieinfo.types import AngularGrid, Medium, Sphere
from mieinfo.viz.plots import (
    plot_comparison,
    plot_eta_vs_na,
    plot_information_pattern,
)


def _small_grid() -> AngularGrid:
    # Small but enough to exercise both θ and φ reductions; keeps runtime low.
    return AngularGrid.full_sphere(40, 16)


def _sphere() -> Sphere:
    return Sphere(radius_m=5e-7, m=1.46 + 0j)  # x ~ 5.9 at 532 nm — cheap Mie sum


def _beam() -> PlaneWave:
    return PlaneWave(Medium(n=1.0, wavelength_vacuum_m=532e-9))


def _nonempty(path) -> None:
    assert path.exists(), f"figure not written: {path}"
    assert path.stat().st_size > 0, f"figure is empty: {path}"


def test_plot_information_pattern_writes_file(tmp_path):
    grid = _small_grid()
    pattern = information_pattern(
        PlaneWaveProvider(), grid, _sphere(), _beam(), np.zeros(3), np.array([1.0, 0.0, 0.0])
    )
    out = tmp_path / "info_pattern.png"
    ret = plot_information_pattern(pattern, str(out))
    assert ret == str(out)
    _nonempty(out)


def test_plot_information_pattern_axial(tmp_path):
    # Axial n̂ exercises the (1−cosθ)² backscatter-weighted density path.
    grid = _small_grid()
    pattern = information_pattern(
        PlaneWaveProvider(), grid, _sphere(), _beam(), np.zeros(3), np.array([0.0, 0.0, 1.0])
    )
    out = tmp_path / "info_pattern_z.png"
    plot_information_pattern(pattern, str(out))
    _nonempty(out)


def test_plot_eta_vs_na_writes_file(tmp_path):
    grid = _small_grid()
    out = tmp_path / "eta_vs_na.png"
    ret = plot_eta_vs_na(
        PlaneWaveProvider(),
        _sphere(),
        _beam(),
        grid,
        axes=[[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        directions=("forward", "backward"),
        path=str(out),
        na_values=np.linspace(0.1, 0.95, 6),  # few points -> fast
    )
    assert ret == str(out)
    _nonempty(out)


def test_plot_eta_vs_na_single_axis_default_na(tmp_path):
    grid = _small_grid()
    out = tmp_path / "eta_single.png"
    # Default na_values sweep; single axis, single direction.
    plot_eta_vs_na(
        PlaneWaveProvider(),
        _sphere(),
        _beam(),
        grid,
        axes=[[1.0, 0.0, 0.0]],
        directions=("forward",),
        path=str(out),
    )
    _nonempty(out)


@dataclass
class _FakeComparison:
    """Duck-typed to INTERFACES.md §9 ComparisonResult for the smoke test."""

    benchmark_key: str
    predicted: float
    reported: float
    within_tolerance: bool
    discrepancy_note: str


def test_plot_comparison_writes_file(tmp_path):
    results = [
        _FakeComparison("tebbenjohanns2019_z_bwd", 0.67, 0.65, True, ""),
        _FakeComparison("magrini2021_x_fwd", 0.55, 0.42, False, "LO mode mismatch"),
        _FakeComparison("dania2022_split", 0.30, 0.31, True, ""),
    ]
    out = tmp_path / "comparison.png"
    ret = plot_comparison(results, str(out))
    assert ret == str(out)
    _nonempty(out)


def test_plot_comparison_single_result(tmp_path):
    results = [_FakeComparison("only_one", 0.5, 0.5, True, "")]
    out = tmp_path / "comparison_single.png"
    plot_comparison(results, str(out))
    _nonempty(out)


def test_plot_comparison_empty_raises(tmp_path):
    with pytest.raises(ValueError):
        plot_comparison([], str(tmp_path / "empty.png"))
