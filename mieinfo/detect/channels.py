"""Multiple detection channels (PHYSICS.md §4.5; INTERFACES.md §6). Time convention
e^{-iωt}.

A real apparatus has several probe beams + collection ports. Each is a **channel** with
its own propagation axis, wavelength, polarization, focusing, and collection cone. Per
PHYSICS.md §4.5 the recipe is:

  1. Work in each channel's OWN beam frame (PHYSICS.md §0, k̂ = beam +z, pol = beam +x).
     The lab-frame sensed direction n̂ is rotated INTO the beam frame with the transpose
     of ``lab_from_beam_frame`` (whose columns are the lab images of beam x,y,z, so R
     maps beam→lab and Rᵀ maps lab→beam).
  2. Evaluate the information pattern for that beam-frame n̂ and apply the channel's cone.
  3. Combine. For INDEPENDENT channels (distinct beams/detectors → independent shot
     noise) the collected Fisher information ADDS:

         F_q^total = Σ_c  η_{q,c}(Ω_c) · F_{q,c}(4π_c)                       (PHYSICS 4.5)
         η_q^total = Σ_c η_{q,c}·F_{q,c}(4π_c) / Σ_c F_{q,c}(4π_c)

     This is the two-beam win: a lab axis collinear with beam A (poorly measured in A's
     forward lobe, (1−cosθ)² axial weight) is transverse to beam B (well measured, sin²θ
     weight). Channels flagged ``independent=False`` share a detector/LO — their noise is
     correlated and information does NOT simply add; that path raises NotImplementedError
     rather than over-adding.

The per-channel F_{q,c}(4π_c) are in the field's arbitrary-but-common units; because all
channels probe the same sphere with the same field normalization, they are directly
comparable and additive.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..glmt.beam import IncidentBeam, lab_from_beam_frame
from ..glmt.scatter import FieldProvider, field_derivative
from ..info.fisher import fisher_total, info_density
from ..info.modes import InformationPattern, combine_direction
from ..types import AngularGrid, Sphere
from .metrics import sql_distance
from .optics import CollectionGeometry, DetectionResult
from .schemes import apply_scheme


@dataclass(frozen=True)
class DetectionChannel:
    """One probe beam + one collection port, oriented in the lab frame (INTERFACES §6)."""
    beam: IncidentBeam                 # carries its own Medium (wavelength); +x-polarized
    propagation_lab: np.ndarray        # unit vector: beam's +z in the lab frame
    polarization_lab: np.ndarray       # unit vector: beam's +x in the lab frame
    geometry: CollectionGeometry
    name: str = ""
    independent: bool = True           # False if it shares a detector/LO with another channel


@dataclass(frozen=True)
class MultiChannelResult:
    per_channel: list[DetectionResult]   # aligned with the input channels
    eta_q_total: float                   # combined collected-info fraction of Σ_c F_c(4π)
    fisher_total_rel: float              # Σ_c η_c · F_c(4π), relative units
    sql_distance: float                  # ≥ 1; 1 == at the SQL


def evaluate_channels(provider: FieldProvider, sphere: Sphere,
                      channels: list[DetectionChannel], n_hat: np.ndarray,
                      grid: AngularGrid) -> MultiChannelResult:
    """Evaluate each channel in its own beam frame and combine (PHYSICS.md §4.5).

    ``n_hat`` is the sensed displacement direction in the LAB frame; each channel
    rotates it into its beam frame with the transpose of ``lab_from_beam_frame``.
    Independent channels' collected Fisher information adds; ``independent=False``
    raises (correlated noise is not modelled).
    """
    n_lab = np.asarray(n_hat, dtype=float)
    n_lab = n_lab / np.linalg.norm(n_lab)

    per_channel: list[DetectionResult] = []
    collected_sum = 0.0   # Σ_c η_c · F_c(4π)
    f_total_sum = 0.0     # Σ_c F_c(4π)

    for ch in channels:
        if not ch.independent:
            raise NotImplementedError(
                f"channel {ch.name!r} is flagged independent=False: correlated shot "
                "noise across a shared detector/LO is not modelled — Fisher information "
                "does not simply add (PHYSICS.md §4.5). Provide independent channels or "
                "implement the covariance combination."
            )

        # Lab → beam frame: Rᵀ, where R = lab_from_beam_frame maps beam → lab.
        R = lab_from_beam_frame(ch.propagation_lab, ch.polarization_lab)
        n_beam = R.T @ n_lab

        # Information pattern in the beam frame for this beam-frame direction.
        deriv = field_derivative(provider, grid, sphere, ch.beam, np.zeros(3),
                                 method="analytic")
        deriv_q = combine_direction(deriv, n_beam)
        density = info_density(deriv_q, grid)
        f_c = fisher_total(density, grid)                      # F_{q,c}(4π_c)
        pattern = InformationPattern(grid=grid, density=density, n_hat=n_beam, f_total=f_c)

        field = provider.field(grid, sphere, ch.beam, np.zeros(3))
        q_sca = provider.q_sca(sphere, ch.beam)
        result = apply_scheme(field, deriv_q, pattern, ch.geometry, q_sca)
        per_channel.append(result)

        collected_sum += result.eta_q * f_c
        f_total_sum += f_c

    eta_q_total = collected_sum / f_total_sum if f_total_sum > 0.0 else 0.0
    return MultiChannelResult(
        per_channel=per_channel,
        eta_q_total=eta_q_total,
        fisher_total_rel=collected_sum,
        sql_distance=sql_distance(eta_q_total),
    )
