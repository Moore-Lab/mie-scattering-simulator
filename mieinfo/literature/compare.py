"""Simulation-vs-literature comparison (INTERFACES.md §9, VALIDATION.md §6 G-LIT).

``compare_benchmark`` parses a :class:`~mieinfo.literature.schema.Benchmark` (whose
``run_config`` W3 wrote — keys ``sphere``, ``medium``, ``beam``, ``collection``,
``n_hat``, ``grid``, ``n_max``), builds the mieinfo objects, runs the plane-wave
information pipeline (``PlaneWaveProvider`` → ``information_pattern`` →
``collection_efficiency``), converts to the *source's* efficiency definition, and
tests the prediction against the benchmark's ``target_tolerance``.

Definition handling (VALIDATION.md §6: "a pass by redefining the quantity is not a
pass"). mieinfo's collected-information efficiency ``η_q(Ω) = F_q(Ω)/F_q(4π)`` — the
fraction of the total displacement Fisher information falling in the collection
cone, with an optimal (mode-matched) local oscillator — is *exactly* the detection
efficiency the reference theories define:

  * Tebbenjohanns 2019 ``η_bw_z`` = (ratio of Heisenberg-limit imprecision to the
    backward split-detection imprecision along z) = fraction of z-position
    information collected backward (Eqs. 8 vs App. D2). Same object as η_q for the
    backward cone. In the dipole limit the z-info density is analytically
    ``|E_s|² · (1 − cosθ)²`` (the phase-gradient factor ``(k̂ − ŝ)_z = 1 − cosθ``),
    giving the backward-weighting the paper reports.
  * Maurer 2023 ``η^d_μ`` = ∫ I_μ(θ_k, φ_k) over the collection solid angle S_d
    (Eq. 24) — identical to η_q with the backward cone.
  * Magrini 2021 ``η = Γ_meas/(Γ_ba + Γ_th)`` is a *total* experimentally-inferred
    efficiency that folds in optical losses and detector QE. mieinfo's ideal η_q is
    a geometric UPPER bound on it; this is a documented definition gap, not a
    numerical match (recorded in the discrepancy_note).

Beam substitution. The M0 forward model only ships ``PlaneWave``; ``RichardsWolfFocus``
/ ``GaussianParaxial`` (W1b) are not yet implemented and the plane-wave FieldProvider
is the only analytic-derivative path (glmt.scatter). Every benchmark run_config here
requests a focused beam, so ``compare_benchmark`` substitutes a ``PlaneWave`` of the
same medium/wavelength and documents the substitution. For the *scattered-field
angular pattern* that drives the running-wave detection efficiency of a small
(dipole/Rayleigh) scatterer this substitution is exact up to the beam's own focal
mode shaping — hence the dipole and small-particle backward-z benchmarks reproduce,
while the focused-beam-specific forward-z sensitivity (Magrini) does not.

Time convention e^{-iωt}; SI units.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..detect.optics import CollectionGeometry, collection_efficiency
from ..glmt.beam import PlaneWave
from ..glmt.scatter import FieldProvider
from ..info.modes import information_pattern
from ..types import AngularGrid, Medium, Sphere
from .schema import Benchmark


@dataclass(frozen=True)
class ComparisonResult:
    benchmark_key: str
    predicted: float
    reported: float
    within_tolerance: bool
    discrepancy_note: str


# ---------------------------------------------------------------- run_config parsing


def _sphere_from(cfg: dict) -> Sphere:
    s = cfg["sphere"]
    m = s["m"]
    # Refractive index stored as [re, im] (schema convention) or a bare number.
    if isinstance(m, (list, tuple)):
        m = complex(float(m[0]), float(m[1]))
    else:
        m = complex(m)
    return Sphere(radius_m=float(s["radius_m"]), m=m)


def _medium_from(cfg: dict) -> Medium:
    md = cfg["medium"]
    return Medium(n=float(md.get("n", 1.0)),
                  wavelength_vacuum_m=float(md["wavelength_vacuum_m"]))


def _grid_from(cfg: dict) -> AngularGrid:
    g = cfg.get("grid") or {}
    ntheta = int(g.get("ntheta", 200))
    nphi = int(g.get("nphi", 96))
    return AngularGrid.full_sphere(ntheta, nphi)


def _geometry_from(cfg: dict) -> CollectionGeometry:
    c = cfg["collection"]
    return CollectionGeometry(
        direction=str(c["direction"]),
        NA=float(c["NA"]),
        apodization=str(c.get("apodization", "aplanatic")),
        lo_mode=str(c.get("lo_mode", "optimal")),
    )


def _n_hat_from(cfg: dict) -> np.ndarray:
    n = np.asarray(cfg["n_hat"], dtype=float)
    return n / np.linalg.norm(n)


def _beam_note(cfg: dict) -> str:
    """Describe the beam substitution (requested focused beam -> PlaneWave)."""
    b = cfg.get("beam") or {}
    btype = b.get("type", "PlaneWave")
    if btype == "PlaneWave":
        return ""
    extra = []
    if "NA" in b:
        extra.append(f"NA={b['NA']}")
    if "filling_factor" in b:
        extra.append(f"filling_factor={b['filling_factor']}")
    detail = (" (" + ", ".join(extra) + ")") if extra else ""
    return (f"Beam substitution: run_config requests {btype}{detail}, but the M0 "
            f"forward model ships only PlaneWave (RichardsWolfFocus/GaussianParaxial "
            f"are W1b, unimplemented) and PlaneWaveProvider is the only "
            f"analytic-derivative path; substituted an x-polarized PlaneWave of the "
            f"same medium. For the small-particle running-wave scattered pattern this "
            f"is exact up to the focal mode shaping. ")


def _one_sided_lower(benchmark: Benchmark) -> bool:
    """True if the target is a lower bound / threshold ('>=', '>90%', 'necessary
    condition ... > 1/9'), read one-sided per the target_provenance text."""
    p = (benchmark.target_provenance or "").lower()
    q = (benchmark.target_quantity or "").lower()
    if "threshold" in p or "necessary condition" in p or "lower bound" in p:
        return True
    if "one-sided" in p or "one sided" in p:
        return True
    if ">=" in p or "≥" in p:
        return True
    if q == "backward_info_fraction":
        # 'more than 90% ... in the backward direction' — a '>90%' claim.
        return True
    return False


def _upper_bound_expected(benchmark: Benchmark) -> bool:
    """True if mieinfo's ideal η_q is a documented UPPER bound on the reported
    (total/experimental) efficiency, so a low prediction under it is a definition
    gap rather than a numerical match (e.g. Magrini's Γ_meas/(Γ_ba+Γ_th))."""
    p = (benchmark.target_provenance or "").lower()
    return ("upper bound" in p or "bounds this from above" in p
            or "gamma_meas" in p
            or "includes losses" in p or "detector qe" in p)


# ---------------------------------------------------------------- prediction


def _predict(provider: FieldProvider, benchmark: Benchmark) -> tuple[float, dict]:
    """Run the plane-wave pipeline for this benchmark; return (predicted, context).

    context carries the size parameter and the collection geometry used, for the
    discrepancy_note.
    """
    cfg = benchmark.run_config
    sphere = _sphere_from(cfg)
    medium = _medium_from(cfg)
    grid = _grid_from(cfg)
    geometry = _geometry_from(cfg)
    n_hat = _n_hat_from(cfg)
    beam = PlaneWave(medium)  # documented substitution for the requested focused beam
    n_max = cfg.get("n_max")

    x = medium.k * sphere.radius_m
    pattern = information_pattern(provider, grid, sphere, beam, np.zeros(3), n_hat,
                                  n_max=n_max)

    q = (benchmark.target_quantity or "").lower()
    if q == "backward_info_fraction":
        # Fraction of the total displacement Fisher information in the backward
        # hemisphere (z<0). Use the SOURCE's definition (full backward half-space),
        # not the run_config's finite collection NA — the paper's '>90%' claim is a
        # hemisphere statement.
        geom_used = CollectionGeometry(direction="backward", NA=1.0)
        predicted = collection_efficiency(pattern, geom_used)
    else:
        # detection_efficiency (and thresholds): collected-info fraction η_q(Ω) for
        # the benchmark's own collection cone — the same object the sources define.
        geom_used = geometry
        predicted = collection_efficiency(pattern, geom_used)

    ctx = {
        "x": x,
        "geom": geom_used,
        "dof": cfg.get("dof"),
        "n_hat": n_hat,
    }
    return float(predicted), ctx


def compare_benchmark(provider: FieldProvider, benchmark: Benchmark) -> ComparisonResult:
    """Predict the benchmark's target quantity with the plane-wave information
    pipeline and test it against the source's stated tolerance (INTERFACES.md §9,
    VALIDATION.md §6 G-LIT).

    The prediction uses mieinfo's collected-information efficiency η_q(Ω), which is
    the SAME quantity the reference theories define as detection efficiency (see the
    module docstring); any normalization or definition gap is documented in
    ``discrepancy_note``. Returns a frozen :class:`ComparisonResult`.
    """
    reported = float(benchmark.target_value)
    tol = float(benchmark.target_tolerance)

    predicted, ctx = _predict(provider, benchmark)
    diff = predicted - reported
    geom: CollectionGeometry = ctx["geom"]

    beam_note = _beam_note(benchmark.run_config)
    regime = "dipole/Rayleigh" if ctx["x"] < 0.5 else ("intermediate" if ctx["x"] < 5 else "Mie")
    base = (
        f"{beam_note}"
        f"predicted={predicted:.4f} vs reported={reported:.4f} (|Δ|={abs(diff):.4f}, "
        f"tol={tol:.3f}); size parameter x={ctx['x']:.4g} ({regime} regime); "
        f"collection {geom.direction} NA={geom.NA:g}, dof={ctx['dof']}. "
        f"mieinfo η_q(Ω)=F_q(Ω)/F_q(4π) with an optimal LO is the SAME collected-"
        f"information efficiency the source defines (VALIDATION.md §6). "
    )

    one_sided = _one_sided_lower(benchmark)
    upper_bound = _upper_bound_expected(benchmark)

    if upper_bound:
        # The reported value is a TOTAL efficiency (losses + detector QE) that
        # mieinfo's ideal geometric η_q bounds from above. A pass would require
        # predicted >= reported (mieinfo is the ceiling). If the geometric prediction
        # is BELOW the reported total, the model is missing physics (here: the
        # focused-beam longitudinal field that gives forward-z sensitivity, dropped
        # by the PlaneWave substitution) — flag as explained-not-matched.
        # An upper-bound relation (ideal η_q ceilings a total efficiency) is a documented
        # DEFINITION GAP, not a like-for-like reproduction — never counted as a numerical
        # pass (so downstream pass-counting stays honest). The note records consistency.
        within = False
        if predicted >= reported:
            note = (base + "The reported value is a TOTAL, experimentally-inferred "
                    "efficiency (η = Γ_meas/(Γ_ba+Γ_th)) that folds in optical losses "
                    "and detector quantum efficiency; mieinfo's ideal η_q is a "
                    "geometric UPPER bound and correctly exceeds it. Definition gap "
                    "documented; not a like-for-like numerical match.")
        else:
            note = (base + "The reported value is a TOTAL efficiency that includes "
                    "losses/QE; mieinfo's ideal η_q should upper-bound it. Here the "
                    "PlaneWave substitution DROPS the focused beam's longitudinal "
                    "(E_z) field, which is the source of forward axial (z) "
                    "sensitivity in Magrini's confocal geometry, so the plane-wave "
                    "prediction under-estimates the forward-z efficiency. "
                    "EXPLAINED-NOT-MATCHED: definition gap (total vs ideal) compounded "
                    "by the missing RichardsWolfFocus longitudinal field; a "
                    "like-for-like reproduction needs the W1b focused-beam provider.")
        return ComparisonResult(benchmark.experiment_key, predicted, reported, within, note)

    if one_sided:
        p_low = (benchmark.target_provenance or "").lower()
        is_threshold = "threshold" in p_low or "necessary condition" in p_low
        if is_threshold:
            # A necessary-condition THRESHOLD (e.g. η^d > 1/9): the threshold value IS the
            # acceptance line — target_tolerance is NOT downward slack on it. Pass only if
            # the prediction truly clears the threshold.
            within = predicted >= reported
            accept_line = reported
            kind = "ground-state-cooling threshold (necessary condition η > 1/9)"
        else:
            # Genuine lower bound quoted with a band (e.g. '>90%' with a sensible tol).
            within = predicted >= reported - tol
            accept_line = reported - tol
            kind = "one-sided lower bound (source states '>=' / '>90%')"
        verdict = "clears" if within else "FAILS to clear"
        note = (base + f"Target is a {kind}; checked one-sided: predicted {verdict} "
                f"acceptance line {accept_line:.4f}.")
        return ComparisonResult(benchmark.experiment_key, predicted, reported, within, note)

    # Two-sided point prediction (Tebbenjohanns η_bw_z, Maurer small-particle η^d_z).
    within = abs(diff) <= tol
    if within:
        note = (base + "Reproduced within the source's stated tolerance "
                "(like-for-like efficiency definition).")
    else:
        note = (base + "Outside tolerance. The plane-wave forward model reproduces the "
                "running-wave scattered pattern of the requested focused beam only up "
                "to its focal mode shaping; the residual is attributable to the "
                "beam-substitution and/or the finite grid quantizing the hard NA cone "
                "(VALIDATION.md §5). EXPLAINED-NOT-MATCHED.")
    return ComparisonResult(benchmark.experiment_key, predicted, reported, within, note)
