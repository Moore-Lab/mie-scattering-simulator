"""glmt.translation is a PENDING stretch item (PHYSICS.md §3.2, §7). The analytic VSWF
translation-addition kernel raises NotImplementedError rather than returning numbers
that fail a genuine (radius-independent, all-degree) reconstruction — two projection
drafts were rejected by adversarial verification (the first used the wrong measure; the
second was a single-radius fit whose reconstruction self-test stayed in the narrow
low-degree / small-radius regime where it happened to work). The validated displacement
path is a quadrature-BSC recompute about the displaced center; its derivative is checked
in tests/glmt/test_derivatives.py. Time convention e^{-iωt}.
"""
from __future__ import annotations

import numpy as np
import pytest

from mieinfo.glmt.translation import axial_translation_coefficients
from mieinfo.types import Medium

K = Medium(n=1.0, wavelength_vacuum_m=532e-9).k


def test_translation_addition_not_implemented():
    """Must raise rather than return silently-wrong (radius-dependent) coefficients."""
    with pytest.raises(NotImplementedError):
        axial_translation_coefficients(K, 0.05 * (2 * np.pi / K), 10, m=1)
