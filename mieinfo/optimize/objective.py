"""Objective evaluation for the detection-geometry optimizer (INTERFACES.md §7;
PHYSICS.md §4.2-4.5, §8). Time convention e^{-iωt}; SI units unless a name says
otherwise.

The single figure of merit is the collected-information fraction ``eta_q`` (PHYSICS
§4.3): maximising it minimises the imprecision-backaction product and drives the
configuration toward the SQL. This module provides the two evaluation primitives the
search layer ranks over:

  * :func:`score_geometry` — ``eta_q`` for one :class:`CollectionGeometry` about a
    single beam, using a CACHED field / derivative / information pattern. Only a cone
    MASK on the fixed field changes across candidates; the field is never re-solved per
    geometry (PHYSICS §7 "avoid recomputing BSCs when only the collection geometry
    changes"). :func:`build_geometry_cache` computes the reusable field once.

  * :func:`sensitivity` — d eta / d param (finite differences) for the design knobs
    ``radius``, ``NA``, ``waist``, and ``m`` about a chosen best geometry. Grid held
    common-mode across the ± evaluations (VALIDATION.md §5: hard-cone eta has
    grid-quantization ~1e-3, which cancels in the difference).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..detect.optics import CollectionGeometry, DetectionResult
from ..detect.schemes import apply_scheme
from ..glmt.beam import IncidentBeam
from ..glmt.scatter import FieldProvider, field_derivative
from ..info.fisher import fisher_total, info_density
from ..info.modes import InformationPattern, combine_direction
from ..types import AngularGrid, Sphere, VectorField


@dataclass(frozen=True)
class GeometryCache:
    """The beam-frame quantities that are INDEPENDENT of the collection geometry.

    Computed once per (provider, sphere, beam, n_hat, grid); every candidate
    :class:`CollectionGeometry` reuses these and only re-applies a cone mask through
    :func:`apply_scheme`. This is the cached-field pattern of PHYSICS §7: the field and
    its displacement derivative are solved a single time, not per NA / direction /
    scheme.
    """

    field: VectorField
    deriv_q: tuple[np.ndarray, np.ndarray]
    pattern: InformationPattern
    q_sca: float


def build_geometry_cache(provider: FieldProvider, sphere: Sphere, beam: IncidentBeam,
                         n_hat: np.ndarray, grid: AngularGrid,
                         n_max: int | None = None) -> GeometryCache:
    """Solve the fixed field / displacement-derivative / information pattern for one
    beam and sensed direction ``n_hat`` (in the beam frame). Reused across all
    collection geometries so the field is solved ONCE (PHYSICS §7)."""
    n = np.asarray(n_hat, dtype=float)
    n = n / np.linalg.norm(n)
    field = provider.field(grid, sphere, beam, np.zeros(3), n_max)
    deriv = field_derivative(provider, grid, sphere, beam, np.zeros(3),
                             method="analytic", n_max=n_max)
    deriv_q = combine_direction(deriv, n)
    density = info_density(deriv_q, grid)
    f_total = fisher_total(density, grid)
    pattern = InformationPattern(grid=grid, density=density, n_hat=n, f_total=f_total)
    q_sca = provider.q_sca(sphere, beam, n_max)
    return GeometryCache(field=field, deriv_q=deriv_q, pattern=pattern, q_sca=q_sca)


def evaluate_geometry(cache: GeometryCache, geometry: CollectionGeometry) -> DetectionResult:
    """Full :class:`DetectionResult` for one geometry from a cached field (only the cone
    mask changes). Thin wrapper over :func:`apply_scheme`."""
    return apply_scheme(cache.field, cache.deriv_q, cache.pattern, geometry, cache.q_sca)


def score_geometry(cache: GeometryCache, geometry: CollectionGeometry) -> float:
    """``eta_q`` for one geometry from the cached field — the ranking objective."""
    return evaluate_geometry(cache, geometry).eta_q


def _best_eta_for_sphere(provider: FieldProvider, sphere: Sphere, beam: IncidentBeam,
                         n_hat: np.ndarray, grid: AngularGrid,
                         geometry: CollectionGeometry,
                         n_max: int | None = None) -> float:
    """``eta_q`` of a fixed geometry for a (possibly perturbed) sphere/beam — the
    finite-difference kernel for the radius/m sensitivities. Rebuilds the cache because
    the field DOES change when the sphere changes, but keeps the grid common-mode."""
    cache = build_geometry_cache(provider, sphere, beam, n_hat, grid, n_max)
    return score_geometry(cache, geometry)


def sensitivity(provider: FieldProvider, sphere: Sphere, beam: IncidentBeam,
                n_hat: np.ndarray, grid: AngularGrid, geometry: CollectionGeometry,
                n_max: int | None = None) -> dict[str, float]:
    """d eta_q / d param about ``geometry`` for the design knobs, by central finite
    differences (PHYSICS §4.2, §8; INTERFACES §7 ``OptResult.sensitivity``).

    Keys ``'radius'`` (per metre), ``'NA'`` (per unit NA), ``'waist'`` (per metre), and
    ``'m'`` (per unit real refractive-index contrast). The grid is held fixed across
    every ± evaluation so the ~1e-3 hard-cone grid-quantization (VALIDATION.md §5)
    cancels common-mode in each difference.

      * radius: recompute the field for radius ± h_r (h_r = 1e-3 · radius).
      * NA:     re-mask the SAME cached field at NA ± h_NA (field is NA-independent).
      * waist:  d eta / d waist. A plane wave has waist = inf and its far-field pattern
                is waist-invariant, so the sensitivity is exactly 0.0 (reported as such
                rather than NaN). A finite-waist beam whose field genuinely depends on
                waist would finite-difference here once such a provider exists; at M0
                the field is waist-independent, hence 0.0.
      * m:      recompute for Re(m) ± h_m (h_m = 1e-3 · |Re(m)|), imaginary part fixed.
    """
    n = np.asarray(n_hat, dtype=float)
    out: dict[str, float] = {}

    # --- radius: field changes; rebuild cache at r ± h_r, common grid ---
    r0 = float(sphere.radius_m)
    h_r = 1e-3 * r0 if r0 > 0.0 else 1e-9
    sph_rp = Sphere(radius_m=r0 + h_r, m=sphere.m)
    sph_rm = Sphere(radius_m=r0 - h_r, m=sphere.m)
    eta_rp = _best_eta_for_sphere(provider, sph_rp, beam, n, grid, geometry, n_max)
    eta_rm = _best_eta_for_sphere(provider, sph_rm, beam, n, grid, geometry, n_max)
    out["radius"] = (eta_rp - eta_rm) / (2.0 * h_r)

    # --- NA: field is NA-independent; re-mask the ONE cached field (PHYSICS §7) ---
    cache = build_geometry_cache(provider, sphere, beam, n, grid, n_max)
    na0 = float(geometry.NA)
    # The cone edge is a hard mask on grid nodes, so d eta / d NA is a staircase in the
    # grid resolution (VALIDATION.md §5). A step smaller than the polar node spacing
    # near the cone edge would move no node and read a spurious 0. Size the step to span
    # several theta nodes: dNA = cos(alpha) * dalpha, dalpha ~ a few mean theta-spacings.
    alpha0 = np.arcsin(min(na0, 1.0 - 1e-9))
    dtheta_mean = float(np.pi / max(grid.theta.size, 1))
    h_na = max(1e-3, 3.0 * dtheta_mean * float(np.cos(alpha0)))
    na_p = min(na0 + h_na, 1.0 - 1e-9)
    na_m = max(na0 - h_na, 1e-9)
    g_p = CollectionGeometry(direction=geometry.direction, NA=na_p,
                             apodization=geometry.apodization, lo_mode=geometry.lo_mode)
    g_m = CollectionGeometry(direction=geometry.direction, NA=na_m,
                             apodization=geometry.apodization, lo_mode=geometry.lo_mode)
    eta_na_p = score_geometry(cache, g_p)
    eta_na_m = score_geometry(cache, g_m)
    out["NA"] = (eta_na_p - eta_na_m) / (na_p - na_m)

    # --- waist: plane wave is waist-invariant (inf) -> exactly 0.0; a finite-waist beam
    #     (e.g. GaussianParaxial via GLMTProvider) genuinely depends on waist, so FD it by
    #     rebuilding the beam at w0 ± h_w. If the beam type cannot be reconstructed with a
    #     waist argument (e.g. RichardsWolfFocus, parameterised by NA not waist), return NaN
    #     rather than a silently-wrong 0.0. ---
    w0 = float(beam.waist_m())
    if not np.isfinite(w0):
        out["waist"] = 0.0
    else:
        h_w = 1e-3 * w0
        try:
            beam_wp = type(beam)(beam.medium, w0 + h_w, beam.polarization)
            beam_wm = type(beam)(beam.medium, w0 - h_w, beam.polarization)
        except Exception:  # noqa: BLE001 — beam not reconstructable by (medium, waist, pol)
            out["waist"] = float("nan")
        else:
            eta_wp = _best_eta_for_sphere(provider, sphere, beam_wp, n, grid, geometry, n_max)
            eta_wm = _best_eta_for_sphere(provider, sphere, beam_wm, n, grid, geometry, n_max)
            out["waist"] = (eta_wp - eta_wm) / (2.0 * h_w)

    # --- m: field changes; rebuild at Re(m) ± h_m, imaginary part held ---
    m0 = complex(sphere.m)
    mre = m0.real
    h_m = 1e-3 * abs(mre) if abs(mre) > 0.0 else 1e-3
    sph_mp = Sphere(radius_m=r0, m=complex(mre + h_m, m0.imag))
    sph_mm = Sphere(radius_m=r0, m=complex(mre - h_m, m0.imag))
    eta_mp = _best_eta_for_sphere(provider, sph_mp, beam, n, grid, geometry, n_max)
    eta_mm = _best_eta_for_sphere(provider, sph_mm, beam, n, grid, geometry, n_max)
    out["m"] = (eta_mp - eta_mm) / (2.0 * h_m)

    return out
