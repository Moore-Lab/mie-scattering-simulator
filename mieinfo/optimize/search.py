"""Detection-geometry optimizer (INTERFACES.md §7; PHYSICS.md §4.2-4.5, §8).

This is the project's *point*: given a scatterer, a probe beam, and a sensed direction
``n_hat``, recommend the collection geometry (and, for a multi-beam apparatus, the SET
of channels) that collects the most information about ``q = r_s . n_hat``. The single
figure of merit is ``eta_q`` (PHYSICS §4.3). Time convention e^{-iωt}; SI units.

Single-channel:  rank every ``CollectionGeometry`` in
``directions_allowed x schemes_allowed x {NA grid up to na_max}`` by ``eta_q``, using
one CACHED field (only the cone mask changes per candidate — PHYSICS §7). Return the
ranked list (best first) and the best.

Multi-channel (``max_channels > 1``):  additionally search over SETS of
``DetectionChannel`` drawn from ``beam_axes_lab x directions x schemes``, maximising the
combined ``eta_q_total`` (PHYSICS §4.5). Because independent channels' Fisher
information ADDS, a set's collected info is the sum of its channels' collected info and
its total is the sum of their 4pi totals; the winning set is found by enumerating
combinations up to ``max_channels`` and is re-scored through ``evaluate_channels`` to
return the exact :class:`MultiChannelResult`.

The key finding this reproduces (PHYSICS §8, MASTER_PLAN F2): axial (beam-collinear)
information is BACKWARD-weighted, so for ``n_hat = z`` about a +z beam the optimizer
picks 'backward' collection; transverse information is forward-weighted, so for
``n_hat = x`` it picks 'forward'. A two-beam set (axes z and x) covers an arbitrary
``n_hat`` better than any single beam.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product

import numpy as np

from ..detect.channels import (
    DetectionChannel,
    MultiChannelResult,
    evaluate_channels,
)
from ..detect.optics import CollectionGeometry, DetectionResult
from ..glmt.beam import IncidentBeam, lab_from_beam_frame
from ..glmt.scatter import FieldProvider
from ..types import AngularGrid, Sphere
from .objective import build_geometry_cache, evaluate_geometry, sensitivity

# Default sampling of NA from just above 0 up to na_max. Coarse and shared across all
# candidates so the ~1e-3 hard-cone grid-quantization (VALIDATION.md §5) is common-mode.
_DEFAULT_NA_SAMPLES = 6


@dataclass(frozen=True)
class Constraints:
    """Achievable-configuration set the optimizer searches (INTERFACES.md §7)."""

    na_max: float
    directions_allowed: tuple[str, ...]        # subset of 'forward'/'backward'/'both'/'split'
    schemes_allowed: tuple[str, ...]           # lo_mode values, e.g. ('optimal','self_homodyne')
    radius_range_m: tuple[float, float] | None = None   # if sphere size is a design variable
    fixed_sphere: Sphere | None = None
    beam_axes_lab: tuple[np.ndarray, ...] | None = None  # candidate probe axes for channels
    max_channels: int = 1                       # >1 optimizes a SET of beams/ports (PHYSICS §4.5)


@dataclass(frozen=True)
class OptResult:
    """Optimizer output (INTERFACES.md §7). ``ranked`` is the single-channel view
    (best first); ``best_channel_set`` / ``best_multichannel`` are populated only when
    ``max_channels > 1``."""

    ranked: list[tuple[CollectionGeometry, DetectionResult]]   # best first
    sensitivity: dict[str, float]              # d eta / d param for radius, NA, waist, m
    best: tuple[CollectionGeometry, DetectionResult]
    best_channel_set: list[DetectionChannel] | None = None
    best_multichannel: MultiChannelResult | None = None


def _na_grid(na_max: float, n_samples: int = _DEFAULT_NA_SAMPLES) -> tuple[float, ...]:
    """NA candidates in (0, na_max], evenly spaced, capped just below 1. Excludes 0 (a
    zero-aperture cone collects nothing)."""
    cap = min(float(na_max), 1.0 - 1e-9)
    if cap <= 0.0:
        return ()
    return tuple(float(v) for v in np.linspace(cap / n_samples, cap, n_samples))


def _candidate_geometries(constraints: Constraints,
                          na_samples: int = _DEFAULT_NA_SAMPLES) -> list[CollectionGeometry]:
    """All ``CollectionGeometry`` in directions x schemes x NA-grid allowed by the
    constraints. The optimizer ranks ``eta_q`` over exactly this set."""
    nas = _na_grid(constraints.na_max, na_samples)
    geoms: list[CollectionGeometry] = []
    for direction, scheme, na in product(constraints.directions_allowed,
                                         constraints.schemes_allowed, nas):
        geoms.append(CollectionGeometry(direction=direction, NA=na, lo_mode=scheme))
    return geoms


def _rank_single_channel(provider: FieldProvider, sphere: Sphere, beam: IncidentBeam,
                         n_hat: np.ndarray, constraints: Constraints, grid: AngularGrid,
                         na_samples: int = _DEFAULT_NA_SAMPLES,
                         n_max: int | None = None
                         ) -> list[tuple[CollectionGeometry, DetectionResult]]:
    """Rank every candidate geometry for one beam by ``eta_q`` (descending), using a
    single cached field (PHYSICS §7). ``n_hat`` is in the beam frame here."""
    cache = build_geometry_cache(provider, sphere, beam, n_hat, grid, n_max)
    scored = [(g, evaluate_geometry(cache, g)) for g in _candidate_geometries(constraints, na_samples)]
    # Sort by eta descending; tie-break toward larger NA then direction name for stable order.
    scored.sort(key=lambda gr: (-gr[1].eta_q, -gr[0].NA, gr[0].direction, gr[0].lo_mode))
    return scored


# --------------------------------------------------------------------------- #
# Multi-channel search
# --------------------------------------------------------------------------- #

@dataclass
class _ChannelCandidate:
    """One candidate channel plus its cached collected/total Fisher info (relative
    units). Independent channels ADD, so a set's totals are sums of these — no field is
    re-solved when combining candidates into sets."""

    channel: DetectionChannel
    collected: float     # eta_{q,c} * F_c(4pi)   (relative)
    f_total: float       # F_c(4pi)               (relative)
    axis_key: bytes      # identifies the beam axis (one field solve per axis)


def _channel_candidates(provider: FieldProvider, sphere: Sphere, beam: IncidentBeam,
                        n_hat: np.ndarray, constraints: Constraints, grid: AngularGrid,
                        na_samples: int = _DEFAULT_NA_SAMPLES,
                        n_max: int | None = None) -> list[_ChannelCandidate]:
    """Build every candidate channel from ``beam_axes_lab x directions x schemes x NA``
    and cache its collected/total Fisher info.

    The field is solved ONCE per beam axis (PHYSICS §7): for a fixed axis, all
    direction/scheme/NA candidates re-mask the same cached beam-frame field. ``n_hat``
    is the LAB-frame sensed direction; each axis rotates it into its own beam frame with
    the transpose of ``lab_from_beam_frame`` (as ``evaluate_channels`` does).
    """
    axes = constraints.beam_axes_lab or (np.asarray(beam_prop_default(beam), float),)
    n_lab = np.asarray(n_hat, dtype=float)
    n_lab = n_lab / np.linalg.norm(n_lab)
    geoms = _candidate_geometries(constraints, na_samples)

    candidates: list[_ChannelCandidate] = []
    for axis in axes:
        prop = np.asarray(axis, dtype=float)
        prop = prop / np.linalg.norm(prop)
        pol = _polarization_for_axis(prop)
        R = lab_from_beam_frame(prop, pol)
        n_beam = R.T @ n_lab
        # ONE field solve per axis; every geometry below re-masks this cache.
        cache = build_geometry_cache(provider, sphere, beam, n_beam, grid, n_max)
        f_total = cache.pattern.f_total
        axis_key = np.round(prop, 9).tobytes()
        for g in geoms:
            res = evaluate_geometry(cache, g)
            ch = DetectionChannel(beam=beam, propagation_lab=prop, polarization_lab=pol,
                                  geometry=g,
                                  name=f"beam[{prop[0]:+.2g},{prop[1]:+.2g},{prop[2]:+.2g}]"
                                       f"_{g.direction}_{g.lo_mode}")
            candidates.append(_ChannelCandidate(channel=ch, collected=res.eta_q * f_total,
                                                f_total=f_total, axis_key=axis_key))
    return candidates


def _polarization_for_axis(prop: np.ndarray) -> np.ndarray:
    """A definite +x-image polarization for a beam along ``prop``: any unit vector
    transverse to ``prop``. Picks the lab axis least aligned with ``prop`` and
    orthogonalizes, so the choice is deterministic and never parallel to ``prop``."""
    prop = prop / np.linalg.norm(prop)
    trials = np.eye(3)
    idx = int(np.argmin(np.abs(prop @ trials.T)))
    v = trials[idx] - (trials[idx] @ prop) * prop
    return v / np.linalg.norm(v)


def beam_prop_default(beam: IncidentBeam) -> np.ndarray:
    """Default lab propagation axis for a single-beam multi-channel search when no
    ``beam_axes_lab`` is given: the beam's own +z."""
    return np.array([0.0, 0.0, 1.0])


def _best_channel_set(candidates: list[_ChannelCandidate], max_channels: int
                      ) -> tuple[list[DetectionChannel], float, float]:
    """Enumerate channel SETS (size 1..max_channels) and return the set that collects
    the most information about ``q`` (PHYSICS §4.5, MASTER_PLAN F3).

    The objective is the total COLLECTED Fisher information ``sum_c eta_{q,c} F_c(4pi)``
    (independent channels add), NOT the normalised ratio ``eta_q_total``. The ratio is
    degenerate for a symmetric multi-beam layout — two complementary beams each collect
    the same FRACTION of their own total, so ``eta_q_total`` is unchanged while the
    absolute collected information (hence the estimate variance, 1/F) improves. Maximising
    collected information is what makes a two-beam apparatus cover an arbitrary ``n_hat``
    better than any single beam. ``eta_q_total`` is reported alongside.

    At most one channel per beam axis per set (one collection port per beam axis), so
    the recommendation is a physical apparatus, not several cones stacked on one axis.
    Ties in collected information are broken toward FEWER channels (a simpler apparatus
    that collects the same information is preferred).
    """
    best_set: list[DetectionChannel] = []
    best_collected = -1.0
    best_eta = 0.0
    best_k = None
    n = len(candidates)
    max_k = min(max_channels, n)
    for k in range(1, max_k + 1):
        for combo in combinations(range(n), k):
            picked = [candidates[i] for i in combo]
            # one channel per beam axis
            axis_keys = [c.axis_key for c in picked]
            if len(set(axis_keys)) != len(axis_keys):
                continue
            f_sum = sum(c.f_total for c in picked)
            if f_sum <= 0.0:
                continue
            collected = sum(c.collected for c in picked)
            eta_total = collected / f_sum
            better = collected > best_collected * (1.0 + 1e-12)
            simpler_tie = (abs(collected - best_collected) <= abs(best_collected) * 1e-12
                           and best_k is not None and k < best_k)
            if better or simpler_tie:
                best_collected = collected
                best_eta = eta_total
                best_set = [c.channel for c in picked]
                best_k = k
    return best_set, best_eta, best_collected


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

def optimize_detection(provider: FieldProvider, sphere: Sphere, beam: IncidentBeam,
                       n_hat: np.ndarray, constraints: Constraints,
                       grid: AngularGrid, na_samples: int = _DEFAULT_NA_SAMPLES,
                       n_max: int | None = None) -> OptResult:
    """Recommend a detection geometry (single-channel) and, for ``max_channels > 1``, a
    channel SET (INTERFACES.md §7; PHYSICS.md §4.2-4.5).

    Single-channel: rank ``CollectionGeometry`` over ``directions x schemes x NA`` by
    ``eta_q`` about ``beam`` for the sensed direction ``n_hat`` (interpreted in the beam
    frame for the single-channel ranking, as the beam computes its pattern in its own
    frame). The field is solved once and only re-masked per candidate (PHYSICS §7).

    Multi-channel (``max_channels > 1``): also search over sets of
    :class:`DetectionChannel` drawn from ``constraints.beam_axes_lab x directions x
    schemes``, maximising ``eta_q_total``. Here ``n_hat`` is the LAB-frame sensed
    direction; each axis rotates it into its own beam frame. For an arbitrary ``n_hat``
    the best set covers ``n_hat``'s components across complementary beams.

    ``sensitivity`` reports d eta / d {radius, NA, waist, m} about the single-channel
    best geometry (finite differences, common grid).
    """
    sph = constraints.fixed_sphere if constraints.fixed_sphere is not None else sphere

    # --- single-channel ranking (always produced) ---
    ranked = _rank_single_channel(provider, sph, beam, n_hat, constraints, grid,
                                  na_samples, n_max)
    if not ranked:
        raise ValueError("no candidate geometries: check na_max/directions/schemes in constraints")
    best = ranked[0]

    # --- sensitivity about the single-channel best geometry ---
    sens = sensitivity(provider, sph, beam, n_hat, grid, best[0], n_max)

    best_channel_set: list[DetectionChannel] | None = None
    best_multichannel: MultiChannelResult | None = None
    if constraints.max_channels > 1:
        # Score candidates with the SAME truncation the returned MultiChannelResult uses:
        # evaluate_channels (below) resolves its own n_max (Wiscombe), so pass None here
        # too — otherwise an explicit caller n_max would rank on one truncation while the
        # re-scored winner reports another. (Single-channel ranking above honors n_max.)
        candidates = _channel_candidates(provider, sph, beam, n_hat, constraints, grid,
                                         na_samples, None)
        chosen, _eta, _collected = _best_channel_set(candidates, constraints.max_channels)
        if chosen:
            best_channel_set = chosen
            # Re-score the winning set through the canonical multi-channel evaluator so
            # the returned MultiChannelResult is exact and self-consistent.
            best_multichannel = evaluate_channels(provider, sph, chosen, n_hat, grid)

    return OptResult(ranked=ranked, sensitivity=sens, best=best,
                     best_channel_set=best_channel_set,
                     best_multichannel=best_multichannel)
