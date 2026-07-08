"""VSWF translation-addition for GLMT sphere displacement along +z
(PHYSICS.md §3.2, W1_scattering_engine.md T1.8/T1.10). Time convention e^{-iωt}.

STATUS: NOT IMPLEMENTED — PENDING stretch item; raises NotImplementedError.

The reliable displacement path shipped by this engine does NOT use the translation-
addition theorem: ``GLMTProvider.field`` recomputes the beam-shape coefficients by
quadrature about the displaced sphere center (``glmt.bsc``, ``glmt.scatter``), and the
displacement derivative is the finite-difference / product-rule BSC gradient in
``glmt.derivatives``, validated against finite differences to ≤ 1e-6 (G-DERIV). That path
is exact (to quadrature precision) for any beam and is what M1 runs on.

Why this raises — two projection attempts, both caught by adversarial verification:
- The FIRST draft projected a shifted multipole onto the new-center basis using the
  ORIGIN sphere's Gauss–Legendre measure (under which the shifted-center harmonics are
  not orthogonal): reconstruction errors of tens to hundreds of percent.
- The SECOND draft used the correct O'-frame measure but extracted the coefficients from
  a SINGLE projection radius r'. A one-radius projection does NOT yield the true,
  radius-independent translation-addition coefficients: an independent field-point
  reconstruction (source degrees up to n_max, general radii) fails by ~3–14% off the
  narrow low-degree / small-radius regime its own self-test used. A degree-n multipole
  about O re-expands about O' with coefficients that must reproduce the field at EVERY
  radius and degree, which a single-radius fit cannot guarantee.

A correct scalar kernel needs the Cruzan/Stein Gaunt-coefficient recurrence. The full
VECTOR kernel that actually translates the GLMT beam-shape coefficients additionally
needs the M↔N mixing ("B") coefficient (a scalar kernel applied to g_tm/g_te alone
misses it, ~30% error) and is heavy at the operative n_max ~ 264 (PHYSICS.md §7). Both
are left for a future pass. Callers needing displaced BSCs must use
``glmt.bsc.bsc_quadrature`` about the displaced sphere center (the validated path).
"""
from __future__ import annotations

import numpy as np

__all__ = ["axial_translation_coefficients"]


def axial_translation_coefficients(k: float, d_m: float, n_max: int, m: int,
                                   r_prime_m: float | None = None,
                                   ntheta: int | None = None) -> np.ndarray:
    """Scalar axial (+z) regular→regular VSWF translation matrix — NOT IMPLEMENTED.

    Raises ``NotImplementedError``. A correct, radius-independent kernel needs the
    Cruzan/Stein Gaunt-coefficient recurrence (single-radius projection attempts did not
    converge to a true translation operator — see the module docstring); the validated
    displacement path is a quadrature-BSC recompute about the displaced center
    (``glmt.bsc.bsc_quadrature``), not this coefficient matrix.
    """
    raise NotImplementedError(
        "Analytic VSWF translation-addition is not implemented (PENDING; PHYSICS.md §7). "
        "A single-radius projection does not yield radius-independent translation "
        "coefficients (reconstruction fails ~3-14% off its narrow validated regime); a "
        "correct kernel needs the Cruzan/Stein Gaunt recurrence. Use "
        "glmt.bsc.bsc_quadrature about the displaced sphere center for displaced BSCs."
    )
