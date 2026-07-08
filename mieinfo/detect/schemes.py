"""Detection schemes as LO mode-overlaps (PHYSICS.md §4.2-4.3, §5; INTERFACES.md §6).
Time convention e^{-iωt}.

The scheme → metrics map. For a coherent (quasi-classical) scattered field collected
over a solid angle Ω, an ideal shot-noise-limited homodyne measurement whose local
oscillator (LO) mode is matched to ∂E_s/∂q saturates the collected Fisher information
F_q(Ω); its efficiency is η = collection_efficiency (κ = 1). A *realistic* LO mode u
collects only the fraction

    κ = |⟨u, s⟩|² / (⟨u, u⟩ ⟨s, s⟩)     with   s ≡ A · ∂E_s/∂q  (signal through optics),

the normalised mode overlap between the LO and the (apodized) displacement-signal
field over the collection cone. The achieved information fraction is then

    η_q = κ · η   (η = F_q(Ω)/F_q(4π)).

Fields are 2-component (θ̂, φ̂) complex vectors on the AngularGrid; the inner product
integrates u_θ* v_θ + u_φ* v_φ over the cone with the solid-angle weights. The lens
apodization A (optics.apodization_weight) multiplies BOTH the signal and the LO on
their way to the detector, so it cancels for the optimal LO (which reshapes to match)
and only reshapes the achievable overlap for a fixed-mode LO.

Implemented LO modes (CollectionGeometry.lo_mode):
  * 'optimal'       — u ∝ A·∂E_s/∂q, κ = 1 (external homodyne with the ideal LO).
  * 'gaussian'      — a realistic single-detector Gaussian homodyne LO: a Gaussian
                      amplitude across the aperture carrying the probe polarization.
                      As a symmetric, single-detector coherent mode it reads the
                      symmetric (axial-type) signal and, by azimuthal parity, is blind
                      to the antisymmetric (transverse-type) signal — a real limitation
                      of a bucket homodyne with a single Gaussian LO.
  * 'quadrant'/'split' — balanced split/quadrant detection: the single effective mode is
                      the ANTISYMMETRIC spatial (sign) mode about the split axis (probe
                      polarization × sign(cos φ)). It reads the antisymmetric
                      (transverse-type) signal and is parity-blind to the symmetric
                      (axial-type) one — the complement of the Gaussian bucket.
  * 'self_homodyne' — the unscattered probe is the reference on a spatially-resolving
                      detector (split/quadrant/camera; the likely default hardware,
                      PHYSICS.md §5). Modelled as a camera that phase-references and
                      matched-filters spatially, so its only loss vs the optimal LO is
                      the CROSS-POLARIZED signal it cannot beat against: κ_self = the
                      fraction of collected |∂E_s/∂q|² that is co-polarized with the
                      probe. The forward/backward emphasis (transverse vs collinear axis)
                      is carried by the cone efficiency η, not by κ.

Modelling note (see gaps_or_issues): 'gaussian'/'split' are single-fixed-mode
overlaps, so their κ is honestly pessimistic (one spatial mode); 'self_homodyne' is a
camera (many modes, only cross-pol lost), so it is more optimistic. A single-element
self-homodyne would additionally be azimuthal-parity-selective like split/gaussian.
"""
from __future__ import annotations

import numpy as np

from ..info.modes import InformationPattern
from ..types import VectorField
from .metrics import backaction_rel, imprecision_rel, sql_distance
from .optics import CollectionGeometry, DetectionResult  # noqa: F401  (re-export DetectionResult)
from .optics import apodization_weight, collection_efficiency, cone_mask


def _cone_inner(u_theta: np.ndarray, u_phi: np.ndarray,
                v_theta: np.ndarray, v_phi: np.ndarray,
                mask: np.ndarray, w_solid: np.ndarray) -> complex:
    """⟨u, v⟩ = Σ_cone (u_θ* v_θ + u_φ* v_φ) dΩ over the collection cone."""
    integrand = (np.conj(u_theta) * v_theta + np.conj(u_phi) * v_phi) * w_solid
    return complex(np.sum(integrand[mask]))


def _probe_reference_mode(grid) -> tuple[np.ndarray, np.ndarray]:
    """Far-field polarization vector of the +x-polarized unscattered probe, projected
    onto (θ̂, φ̂):  u_θ = cosθ·cosφ,  u_φ = −sinφ  (PHYSICS.md §0 spherical basis).

    This is the angular structure of the imaging beam used as its own LO in
    self-homodyne, and the polarization carried by the Gaussian/split LOs. It is real and
    frame-defined (beam +x polarization); the channel rotates the geometry, not this.
    """
    TH = grid.theta[:, None]
    PH = grid.phi[None, :]
    u_theta = np.cos(TH) * np.cos(PH) * np.ones_like(PH)
    u_phi = -np.sin(PH) * np.ones_like(TH)
    return u_theta.astype(complex), u_phi.astype(complex)


def _probe_polarization_hat(grid) -> tuple[np.ndarray, np.ndarray]:
    """Unit (θ̂, φ̂) polarization direction of the probe reference at each point — the
    normalized :func:`_probe_reference_mode`. Points where the probe field vanishes
    (a measure-zero set, e.g. φ=π/2 on-axis) get a zero direction and contribute
    nothing to the co-polarized projection."""
    u_theta, u_phi = _probe_reference_mode(grid)
    norm = np.sqrt(np.abs(u_theta) ** 2 + np.abs(u_phi) ** 2)
    safe = norm > 0.0
    p_theta = np.zeros_like(u_theta)
    p_phi = np.zeros_like(u_phi)
    p_theta[safe] = u_theta[safe] / norm[safe]
    p_phi[safe] = u_phi[safe] / norm[safe]
    return p_theta, p_phi


def _aperture_radius(grid, geometry: CollectionGeometry) -> np.ndarray:
    """Aperture radial coordinate ρ = sinθ measured from each collected lobe's axis
    (forward: sinθ about +z; backward: sin(π−θ) = sinθ about −z). (Ntheta, Nphi)."""
    TH = np.broadcast_to(grid.theta[:, None], grid.w_solid.shape)
    return np.sin(TH)  # sin θ == sin(π−θ), same for forward and backward lobes


def _gaussian_lo_mode(grid, geometry: CollectionGeometry) -> tuple[np.ndarray, np.ndarray]:
    """Realistic Gaussian LO: a Gaussian amplitude across the aperture (in ρ = sinθ)
    carrying the probe polarization. The 1/e² aperture radius is sin(α) (the cone
    edge), so the Gaussian fills the collected cone with a soft, non-flat profile —
    the generic mismatch a real single-mode LO has against the structured signal."""
    rho = _aperture_radius(grid, geometry)
    rho_edge = min(geometry.NA, 1.0)                     # sin(α) at the cone edge
    w_ap = rho_edge if rho_edge > 0.0 else 1.0
    amp = np.exp(-(rho ** 2) / (w_ap ** 2))              # Gaussian aperture apodization
    u_theta, u_phi = _probe_reference_mode(grid)
    return amp * u_theta, amp * u_phi


def _split_lo_mode(grid, geometry: CollectionGeometry) -> tuple[np.ndarray, np.ndarray]:
    """Balanced split/quadrant LO: the antisymmetric spatial mode. Balanced difference
    of the two aperture halves = a sign(cos φ) sign mode (split axis along +x, the
    polarization/transverse axis) carrying the probe polarization. Its overlap with the
    displacement signal picks out the antisymmetric part — large for transverse motion,
    which produces an antisymmetric far-field response, near-zero for the symmetric
    forward axial response."""
    PH = grid.phi[None, :]
    sign_mode = np.sign(np.cos(PH)) * np.ones((grid.theta.size, grid.phi.size))
    u_theta, u_phi = _probe_reference_mode(grid)
    return sign_mode * u_theta, sign_mode * u_phi


def _single_mode_overlap(u_theta: np.ndarray, u_phi: np.ndarray,
                         sig_theta: np.ndarray, sig_phi: np.ndarray,
                         mask: np.ndarray, W: np.ndarray) -> float:
    """κ = |⟨u, s⟩|²/(⟨u,u⟩⟨s,s⟩) for a single fixed LO spatial mode u over the cone."""
    ss = _cone_inner(sig_theta, sig_phi, sig_theta, sig_phi, mask, W).real
    if ss <= 0.0:
        return 0.0
    uu = _cone_inner(u_theta, u_phi, u_theta, u_phi, mask, W).real
    if uu <= 0.0:
        return 0.0
    us = _cone_inner(u_theta, u_phi, sig_theta, sig_phi, mask, W)
    kappa = (abs(us) ** 2) / (uu * ss)
    return float(min(max(kappa, 0.0), 1.0))


def _self_homodyne_kappa(grid, geometry: CollectionGeometry,
                         sig_theta: np.ndarray, sig_phi: np.ndarray,
                         mask: np.ndarray, W: np.ndarray) -> float:
    """Camera self-homodyne (PHYSICS.md §5): the probe is the phase reference and the
    detector resolves space, so the only loss vs the optimal LO is the cross-polarized
    signal that cannot beat against the probe. κ = (co-polarized signal power)/(total
    signal power) over the cone, where "co-polarized" is the projection onto the probe
    polarization direction p̂:  κ = ⟨|p̂·s|²⟩ / ⟨|s|²⟩."""
    p_theta, p_phi = _probe_polarization_hat(grid)
    s_par = np.conj(p_theta) * sig_theta + np.conj(p_phi) * sig_phi  # complex scalar
    num = float(np.sum((np.abs(s_par) ** 2 * W)[mask]))
    den = float(np.sum(((np.abs(sig_theta) ** 2 + np.abs(sig_phi) ** 2) * W)[mask]))
    if den <= 0.0:
        return 0.0
    return float(min(max(num / den, 0.0), 1.0))


def lo_overlap(field: VectorField, deriv_q: tuple[np.ndarray, np.ndarray],
               geometry: CollectionGeometry) -> float:
    """Mode-overlap efficiency κ ∈ [0, 1] between the detection scheme and the (apodized)
    displacement-signal field over the collection cone (PHYSICS.md §5).

    For the single-fixed-mode schemes ('optimal', 'gaussian', 'split'/'quadrant'):
        κ = |⟨u, s⟩|² / (⟨u, u⟩ ⟨s, s⟩),   s = A · ∂E_s/∂q,   κ = 1 for 'optimal' (u = s).
    For 'self_homodyne' (a spatially-resolving camera): κ = co-polarized signal fraction.

    `field` (the unscattered/total field) is accepted for interface symmetry; the
    current models use the analytic probe reference, not `field`.
    """
    grid = field.grid
    mask = cone_mask(grid, geometry)
    W = grid.w_solid
    A = apodization_weight(grid, geometry)              # real amplitude, zero outside cone

    dE_theta_q, dE_phi_q = deriv_q
    sig_theta = A * dE_theta_q                          # signal through the optics
    sig_phi = A * dE_phi_q

    mode = geometry.lo_mode
    if mode == "optimal":
        return _single_mode_overlap(sig_theta, sig_phi, sig_theta, sig_phi, mask, W)
    if mode == "gaussian":
        u_theta, u_phi = _gaussian_lo_mode(grid, geometry)
        return _single_mode_overlap(u_theta, u_phi, sig_theta, sig_phi, mask, W)
    if mode in ("split", "quadrant"):
        u_theta, u_phi = _split_lo_mode(grid, geometry)
        return _single_mode_overlap(u_theta, u_phi, sig_theta, sig_phi, mask, W)
    if mode == "self_homodyne":
        return _self_homodyne_kappa(grid, geometry, sig_theta, sig_phi, mask, W)
    raise ValueError(
        f"unknown lo_mode {mode!r} (use 'optimal', 'gaussian', 'split', "
        "'quadrant', or 'self_homodyne')"
    )


def apply_scheme(field: VectorField, deriv_q: tuple[np.ndarray, np.ndarray],
                 pattern: InformationPattern, geometry: CollectionGeometry,
                 q_sca: float) -> DetectionResult:
    """Map (E_s, ∂E_s/∂q, information pattern, geometry, Q_sca) → DetectionResult
    (PHYSICS.md §4.2-4.3, §5). The achieved information fraction is

        η_q = κ · η,   η = collection_efficiency (optimal-LO cone fraction),
                       κ = lo_overlap (LO mode match; 1 for the optimal LO),

    and the metrics follow from η_q and Q_sca (metrics.py):
      s_imp_rel  = 1/η_q,   γ_ba_rel = Q_sca,   sql_distance = 1/η_q.
    """
    eta = collection_efficiency(pattern, geometry)      # optimal-LO cone fraction (κ=1)
    kappa = lo_overlap(field, deriv_q, geometry)
    eta_q = kappa * eta
    return DetectionResult(
        geometry=geometry,
        eta_q=eta_q,
        s_imp_rel=imprecision_rel(eta_q),
        gamma_ba_rel=backaction_rel(q_sca),
        sql_distance=sql_distance(eta_q),
    )
