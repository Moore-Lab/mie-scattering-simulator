"""VSWF translation-addition (Cruzan/Stein) for GLMT sphere displacement
(PHYSICS.md §3.2, W1_scattering_engine.md T1.8). Time convention e^{-iωt}; SI units.

STATUS: NOT IMPLEMENTED — PENDING stretch item.

The *reliable* displacement path shipped by this engine does NOT use the translation-
addition theorem: ``GLMTProvider.field`` recomputes the beam-shape coefficients by
quadrature about the displaced sphere center (``glmt.bsc``, ``glmt.scatter``) and the
displacement derivative is the product-rule / finite-difference form in
``glmt.derivatives``, validated against finite differences to ≤ 1e-6 (G-DERIV). That
path is exact (to quadrature precision) for any beam and is what M1 runs on.

The fast analytic vector translation-addition at n_max ~ 264 (translation coefficients
scale steeply in order, PHYSICS.md §7) is left PENDING. An earlier draft here computed
"axial scalar translation coefficients" by projecting a shifted multipole onto the
new-center basis using the ORIGIN sphere's Gauss-Legendre weights — but that measure
does NOT make the shifted-center harmonics orthogonal, so the result was mathematically
wrong (reconstruction and group-property errors of tens to hundreds of percent). Rather
than ship silently-wrong numbers, the entry point raises ``NotImplementedError`` until a
correct Cruzan/Stein (Gaunt-recurrence) implementation lands with reconstruction- and
group-property tests. Callers needing displaced BSCs must use
``glmt.bsc.bsc_quadrature`` about the displaced center (the validated path).
"""
from __future__ import annotations

import numpy as np

__all__ = ["axial_translation_coefficients"]


def axial_translation_coefficients(k: float, d_m: float, n_max: int, m: int,
                                   ntheta: int | None = None) -> np.ndarray:
    """Axial (z-shift) scalar translation matrix ``alpha_{n,nu}^{(m)}(k d)``, fixed m.

    NOT IMPLEMENTED. The correct Cruzan/Stein axial translation kernel is a PENDING
    stretch item (see the module docstring); the validated displacement path is a
    quadrature-BSC recompute about the displaced center (``glmt.bsc.bsc_quadrature``),
    not this coefficient matrix. Raises rather than returning silently-wrong values.
    """
    raise NotImplementedError(
        "Analytic VSWF translation-addition is not implemented (PENDING stretch item, "
        "PHYSICS.md §7). Use glmt.bsc.bsc_quadrature about the displaced sphere center "
        "for displaced BSCs; see the mieinfo/glmt/translation.py module docstring."
    )
