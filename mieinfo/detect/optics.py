"""Collection optics: cones and collected-information efficiency (INTERFACES.md §6,
PHYSICS.md §4.2, §5). Collection is in vacuum here, so NA = sin(alpha) (capped < 1).
Time convention e^{-iωt}.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..info.modes import InformationPattern
from ..types import AngularGrid


@dataclass(frozen=True)
class CollectionGeometry:
    direction: str                    # 'forward' | 'backward' | 'both' | 'split'
    NA: float                         # numerical aperture n*sin(alpha); vacuum -> sin(alpha)
    apodization: str = "aplanatic"    # 'aplanatic' (sqrt(cos)) | 'none'
    lo_mode: str = "optimal"          # 'optimal' | 'gaussian' | 'self_homodyne' | 'quadrant'


@dataclass(frozen=True)
class DetectionResult:
    """Figure-of-merit for one collection geometry + scheme (INTERFACES.md §6,
    PHYSICS.md §4.2-4.3). Units are relative (ratios of Fisher information / Q_sca)."""
    geometry: CollectionGeometry
    eta_q: float          # collected info / total info, incl. LO overlap κ (∈ [0, 1])
    s_imp_rel: float      # imprecision PSD in units of 1/f_total (= 1/η_q, ≥ 1)
    gamma_ba_rel: float   # backaction rate ~ Q_sca (relative units)
    sql_distance: float   # ≥ 1; 1 == at the standard quantum limit


def cone_mask(grid: AngularGrid, geometry: CollectionGeometry) -> np.ndarray:
    """Boolean (Ntheta, Nphi) mask for the collection cone(s). Forward = cone about +z
    (small θ), backward = cone about −z (θ→π); 'both'/'split' collect the union."""
    alpha = np.arcsin(min(geometry.NA, 1.0))
    TH = np.broadcast_to(grid.theta[:, None], grid.w_solid.shape)
    d = geometry.direction
    if d == "forward":
        return TH <= alpha
    if d == "backward":
        return TH >= (np.pi - alpha)
    if d in ("both", "split"):
        return (TH <= alpha) | (TH >= (np.pi - alpha))
    raise ValueError(f"unknown collection direction {d!r}")


def apodization_weight(grid: AngularGrid, geometry: CollectionGeometry) -> np.ndarray:
    """Amplitude apodization applied by the collection lens across the aperture,
    combined with hard aperture truncation (PHYSICS.md §5, "Collection optics").

    Returns a real (Ntheta, Nphi) amplitude weight A(θ, φ) that multiplies the far
    field on its way to the detector:

      * ``apodization='aplanatic'`` — an energy-conserving aplanatic (Abbe-sine) lens
        maps the spherical far field onto a flat aperture with the classic √cos θ
        amplitude factor. Measured from each lobe's own axis: √cos θ forward, √|cos θ|
        (= √cos(π−θ)) backward, so both lobes are treated symmetrically.
      * ``apodization='none'`` — unit amplitude across the aperture (pure truncation).

    Aperture truncation is the hard cut of ``cone_mask``: outside the collected
    cone(s) the weight is zero. Because the aplanatic factor and truncation are real
    amplitude transforms applied to *both* the signal and any co-propagating LO, they
    cancel for the optimal LO (which reshapes to match) but reshape the achievable
    mode for a fixed-mode LO — that is where they change κ in :mod:`schemes`.
    """
    mask = cone_mask(grid, geometry)
    A = np.zeros(grid.w_solid.shape, dtype=float)
    apo = geometry.apodization
    if apo == "none":
        A[mask] = 1.0
        return A
    if apo == "aplanatic":
        TH = np.broadcast_to(grid.theta[:, None], grid.w_solid.shape)
        # √|cos θ|: forward lobe uses cos θ, backward lobe uses cos(π−θ)=−cos θ.
        A[mask] = np.sqrt(np.abs(np.cos(TH[mask])))
        return A
    raise ValueError(f"unknown apodization {apo!r} (use 'aplanatic' or 'none')")


def collection_efficiency(pattern: InformationPattern, geometry: CollectionGeometry,
                          use_apodization: bool = False) -> float:
    """η_q(Ω) = F_q(Ω) / F_q(4π) ∈ [0, 1] — the fraction of total information collected
    by the cone, with an optimal (∝ ∂E_s/∂q) local oscillator (PHYSICS.md §4.2, §4.3).

    The optimal LO saturates F_q(Ω) regardless of the lens amplitude profile (it
    reshapes to whatever the optics deliver), so ``use_apodization=False`` (the
    default, and the M0 behaviour) is the physically correct optimal-LO efficiency.
    ``use_apodization=True`` additionally applies the lens amplitude²-weighting to the
    collected information — the fraction a detector sees when the *aperture* (not the
    LO) is the mode-defining element, e.g. direct/intensity collection through the
    aplanatic lens. Both remain in [0, 1] and both give η=1 for a full 4π cone with no
    apodization loss.
    """
    mask = cone_mask(pattern.grid, geometry)
    W = pattern.grid.w_solid
    total = float(np.sum(pattern.density * W))
    if total <= 0.0:
        return 0.0
    if not use_apodization:
        collected = float(np.sum(pattern.density[mask] * W[mask]))
        return collected / total
    A2 = apodization_weight(pattern.grid, geometry) ** 2  # intensity apodization
    collected = float(np.sum(pattern.density * A2 * W))
    return collected / total
