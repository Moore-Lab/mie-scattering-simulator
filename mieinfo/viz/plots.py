"""Matplotlib figures for mieinfo objects (ARCHITECTURE.md §5 provenance).

Every figure embeds its generating configuration in the title or a caption so the
artifact is reproducible from what it carries ("every figure stores its generating
config"). All plotting uses the non-interactive Agg backend — nothing is displayed;
figures are written to disk and closed. Time convention e^{-iωt}; SI units unless a
name says otherwise (``_um``, ``_deg``).

Public API
----------
- ``plot_information_pattern(pattern, path)`` — φ-averaged polar profile + a
  (θ, φ) heatmap of the information radiation pattern dF_q/dΩ (PHYSICS.md §4.1).
- ``plot_eta_vs_na(provider, sphere, beam, grid, axes, directions, path)`` —
  collected-information efficiency η_q(NA) curves per sensed axis and collection
  direction (PHYSICS.md §4.2, §8).
- ``plot_comparison(comparison_results, path)`` — predicted vs. reported quantities
  for a set of literature benchmarks (INTERFACES.md §9).
"""
from __future__ import annotations

from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")  # headless: never open a display

import matplotlib.pyplot as plt  # noqa: E402  (after backend selection)
import numpy as np  # noqa: E402

from ..detect.optics import CollectionGeometry, collection_efficiency  # noqa: E402
from ..info.modes import InformationPattern, information_pattern  # noqa: E402
from ..types import AngularGrid, Sphere  # noqa: E402

__all__ = [
    "plot_information_pattern",
    "plot_eta_vs_na",
    "plot_comparison",
]

# Caption styling shared across figures (small, unobtrusive provenance footer).
_CAPTION_KW = dict(ha="center", va="bottom", fontsize=7, color="0.35", wrap=True)


def _fmt_vec(v: np.ndarray) -> str:
    """Compact fixed-precision string for a small vector, e.g. '[1, 0, 0]'."""
    return "[" + ", ".join(f"{x:g}" for x in np.asarray(v, dtype=float).ravel()) + "]"


def _add_caption(fig, text: str) -> None:
    """Attach a provenance caption along the bottom of the figure."""
    fig.text(0.5, 0.005, text, **_CAPTION_KW)


def _phi_average(density: np.ndarray, grid: AngularGrid) -> np.ndarray:
    """φ-average of a (Ntheta, Nphi) density, weighted by the φ solid-angle spacing.

    Reduces to the mean over uniform φ nodes; robust to any per-column weighting the
    grid may carry. Returns a (Ntheta,) profile in θ order matching ``grid.theta``.
    """
    w = grid.w_solid
    wsum = np.sum(w, axis=1)
    wsum = np.where(wsum > 0.0, wsum, 1.0)
    return np.sum(density * w, axis=1) / wsum


def plot_information_pattern(pattern: InformationPattern, path: str) -> str:
    """Plot the information radiation pattern dF_q/dΩ (PHYSICS.md §4.1) and save to ``path``.

    Two panels:
      * a polar plot of the **φ-averaged** density vs. polar angle θ (0 = forward /
        +z, π = backward / −z), showing the forward-vs-backward weighting that makes
        transverse info forward-peaked and axial info backscatter-peaked;
      * a (θ, φ) heatmap of the full density.

    The sensed direction n̂, the total Fisher information F_q(4π), and the grid size
    are embedded in the title/caption for provenance. Returns ``path``.
    """
    grid = pattern.grid
    density = np.asarray(pattern.density, dtype=float)
    theta = np.asarray(grid.theta, dtype=float)
    phi = np.asarray(grid.phi, dtype=float)

    prof = _phi_average(density, grid)
    # Sort by θ so the polar/line trace is monotone regardless of node ordering.
    order = np.argsort(theta)
    th_s = theta[order]
    prof_s = prof[order]

    fig = plt.figure(figsize=(11, 4.6))
    ax_polar = fig.add_subplot(1, 2, 1, projection="polar")
    ax_polar.plot(th_s, prof_s, color="C0", lw=1.8)
    ax_polar.fill(th_s, prof_s, color="C0", alpha=0.15)
    # θ=0 (forward, +z) at top, increasing clockwise so backward (π) sits at bottom.
    ax_polar.set_theta_zero_location("N")
    ax_polar.set_theta_direction(-1)
    ax_polar.set_thetamin(0.0)
    ax_polar.set_thetamax(180.0)
    ax_polar.set_title("φ-averaged dF_q/dΩ  (θ: 0=fwd/+z, π=bwd/−z)", fontsize=9)

    ax_hm = fig.add_subplot(1, 2, 2)
    # imshow expects rows increasing downward; put θ on y, φ on x, sorted for a
    # readable axis. extent in degrees.
    phi_order = np.argsort(phi)
    Z = density[np.ix_(order, phi_order)]
    im = ax_hm.imshow(
        Z,
        aspect="auto",
        origin="lower",
        extent=[
            float(np.degrees(phi[phi_order][0])),
            float(np.degrees(phi[phi_order][-1])),
            float(np.degrees(th_s[0])),
            float(np.degrees(th_s[-1])),
        ],
        cmap="viridis",
    )
    ax_hm.set_xlabel("φ (deg)")
    ax_hm.set_ylabel("θ (deg)")
    ax_hm.set_title("dF_q/dΩ over the sphere", fontsize=9)
    fig.colorbar(im, ax=ax_hm, fraction=0.046, pad=0.04, label="dF_q/dΩ (arb.)")

    fig.suptitle(
        f"Information radiation pattern — n̂ = {_fmt_vec(pattern.n_hat)}", fontsize=11
    )
    caption = (
        f"config: n̂={_fmt_vec(pattern.n_hat)} | "
        f"grid={theta.size}×{phi.size} (θ×φ) | "
        f"F_q(4π)={pattern.f_total:.4g} (arb.)"
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    _add_caption(fig, caption)
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def plot_eta_vs_na(
    provider,
    sphere: Sphere,
    beam,
    grid: AngularGrid,
    axes: Sequence,
    directions: Iterable[str] = ("forward", "backward"),
    path: str = "eta_vs_na.png",
    na_values: Sequence[float] | None = None,
    r_s_m: np.ndarray | None = None,
    apodization: str = "aplanatic",
) -> str:
    """Plot collected-information efficiency η_q(NA) and save to ``path``.

    For each sensed axis in ``axes`` and each collection ``directions`` entry, builds
    the InformationPattern once (via ``information_pattern``) and sweeps NA through
    ``collection_efficiency`` (PHYSICS.md §4.2). One curve per (axis, direction);
    forward curves solid, backward dashed. Reproduces the §8 story: transverse info is
    forward-collected, axial info is backward-collected.

    ``axes`` are sensed directions n̂ (3-vectors); ``na_values`` defaults to a fine
    sweep in (0, 1). The sphere/beam/grid and the NA grid are embedded in the caption.
    Returns ``path``.
    """
    if na_values is None:
        na_values = np.linspace(0.05, 0.99, 24)
    na = np.asarray(na_values, dtype=float)
    r_s = np.zeros(3) if r_s_m is None else np.asarray(r_s_m, dtype=float)
    axes = [np.asarray(a, dtype=float) for a in axes]
    directions = list(directions)

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    linestyles = {"forward": "-", "backward": "--", "both": ":", "split": "-."}
    markers = {"forward": "o", "backward": "s", "both": "^", "split": "D"}

    for ai, axis in enumerate(axes):
        pattern = information_pattern(provider, grid, sphere, beam, r_s, axis)
        color = f"C{ai}"
        for d in directions:
            etas = [
                collection_efficiency(
                    pattern, CollectionGeometry(direction=d, NA=float(v), apodization=apodization)
                )
                for v in na
            ]
            ax.plot(
                na,
                etas,
                linestyle=linestyles.get(d, "-"),
                marker=markers.get(d, "o"),
                markersize=3.5,
                color=color,
                lw=1.6,
                label=f"n̂={_fmt_vec(axis)}, {d}",
            )

    ax.set_xlabel("collection NA  (vacuum: sin α)")
    ax.set_ylabel("η_q  = F_q(Ω) / F_q(4π)")
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlim(0.0, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="best")

    x = _size_parameter_safe(sphere, beam)
    lam = _wavelength_str(beam)
    ax.set_title("Collected-information efficiency vs. NA", fontsize=11)
    caption = (
        f"config: sphere a={sphere.radius_m:g} m, m={sphere.m} | {lam} | "
        f"x={x} | grid={np.asarray(grid.theta).size}×{np.asarray(grid.phi).size} (θ×φ) | "
        f"apodization={apodization} | r_s={_fmt_vec(r_s)} m | "
        f"provider={type(provider).__name__}"
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    _add_caption(fig, caption)
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def plot_comparison(comparison_results: Sequence, path: str) -> str:
    """Plot predicted vs. reported quantities for literature benchmarks (INTERFACES.md §9).

    ``comparison_results`` is a sequence of objects duck-typed to ``ComparisonResult``:
    each must expose ``benchmark_key`` (str), ``predicted`` (float), ``reported``
    (float), ``within_tolerance`` (bool), and ``discrepancy_note`` (str). Two panels:
      * a grouped bar chart of predicted vs. reported per benchmark, bars annotated
        pass/fail by ``within_tolerance``;
      * a predicted-vs-reported scatter with the y=x agreement line.

    The benchmark keys and pass/fail counts are embedded in the caption for
    provenance. Returns ``path``. Raises ValueError on an empty input.
    """
    results = list(comparison_results)
    if not results:
        raise ValueError("plot_comparison requires at least one ComparisonResult")

    keys = [str(getattr(r, "benchmark_key")) for r in results]
    predicted = np.array([float(getattr(r, "predicted")) for r in results], dtype=float)
    reported = np.array([float(getattr(r, "reported")) for r in results], dtype=float)
    ok = np.array([bool(getattr(r, "within_tolerance")) for r in results], dtype=bool)

    n = len(results)
    idx = np.arange(n)

    fig, (ax_bar, ax_sc) = plt.subplots(1, 2, figsize=(12, 5.0))

    width = 0.38
    ax_bar.bar(idx - width / 2, predicted, width, label="predicted", color="C0")
    ax_bar.bar(idx + width / 2, reported, width, label="reported", color="C1")
    ax_bar.set_xticks(idx)
    ax_bar.set_xticklabels(keys, rotation=30, ha="right", fontsize=8)
    ax_bar.set_ylabel("value")
    ax_bar.set_title("Predicted vs. reported", fontsize=10)
    ax_bar.legend(fontsize=8)
    top = max(float(np.max(np.abs(np.concatenate([predicted, reported])))), 1e-12)
    for i in range(n):
        mark = "✓" if ok[i] else "✗"
        col = "green" if ok[i] else "red"
        ax_bar.annotate(
            mark,
            (idx[i], max(predicted[i], reported[i])),
            textcoords="offset points",
            xytext=(0, 3),
            ha="center",
            fontsize=10,
            color=col,
        )
    ax_bar.set_ylim(top=top * 1.15 if top > 0 else 1.0)

    # Scatter: predicted (x) vs reported (y) with y=x line.
    lo = float(min(np.min(predicted), np.min(reported)))
    hi = float(max(np.max(predicted), np.max(reported)))
    if hi <= lo:
        hi = lo + 1.0
    pad = 0.05 * (hi - lo)
    line = np.array([lo - pad, hi + pad])
    ax_sc.plot(line, line, color="0.5", lw=1.0, ls="--", label="y = x (agreement)")
    ax_sc.scatter(
        predicted[ok], reported[ok], color="green", marker="o", s=45,
        label="within tolerance", zorder=3,
    )
    ax_sc.scatter(
        predicted[~ok], reported[~ok], color="red", marker="X", s=55,
        label="out of tolerance", zorder=3,
    )
    for i in range(n):
        ax_sc.annotate(
            keys[i], (predicted[i], reported[i]),
            textcoords="offset points", xytext=(4, 4), fontsize=7, color="0.3",
        )
    ax_sc.set_xlabel("predicted")
    ax_sc.set_ylabel("reported")
    ax_sc.set_xlim(line[0], line[1])
    ax_sc.set_ylim(line[0], line[1])
    ax_sc.set_aspect("equal", adjustable="box")
    ax_sc.grid(True, alpha=0.3)
    ax_sc.set_title("Agreement", fontsize=10)
    ax_sc.legend(fontsize=8, loc="best")

    n_pass = int(np.count_nonzero(ok))
    fig.suptitle("Literature comparison — predicted vs. reported", fontsize=12)
    caption = (
        f"config: {n} benchmark(s) [{', '.join(keys)}] | "
        f"within tolerance: {n_pass}/{n}"
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    _add_caption(fig, caption)
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


# --- provenance helpers (best-effort; never raise on unusual providers/beams) ---


def _size_parameter_safe(sphere: Sphere, beam) -> str:
    """x = k·a if the beam carries a Medium with a wavenumber; else 'n/a'."""
    try:
        k = beam.medium.k
        return f"{k * sphere.radius_m:.4g}"
    except Exception:  # pragma: no cover - provenance string only
        return "n/a"


def _wavelength_str(beam) -> str:
    try:
        lam = beam.medium.wavelength_vacuum_m
        n = beam.medium.n
        return f"λ_vac={lam * 1e9:g} nm, n={n:g}"
    except Exception:  # pragma: no cover - provenance string only
        return "λ=n/a"
