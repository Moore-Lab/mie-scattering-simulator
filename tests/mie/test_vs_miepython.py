"""Cross-check the promoted engine against the miepython oracle across the operative
size range (x up to 236 at 532 nm), including the large-x burn-in regime. Mirrors
prototype/validate.py Gates 1-2. Skipped if miepython is unavailable. Time convention
e^{-iωt}."""
from __future__ import annotations

import numpy as np
import pytest

mp = pytest.importorskip("miepython")

from mieinfo.mie import plane_wave as pw

M = 1.45 + 0j
XS = [0.3, 1.0, 5.0, 14.7631, 29.5301, 59.05, 118.1, 236.2]


@pytest.mark.parametrize("x", XS, ids=[f"x{x:g}" for x in XS])
def test_efficiencies_vs_miepython(x):
    qe, qs, qb, g = pw.efficiencies(M, x)
    Qe, Qs, Qb, G = mp.efficiencies_mx(M, x)
    assert abs(qe - Qe) / abs(Qe) <= 1e-9
    assert abs(qs - Qs) / abs(Qs) <= 1e-9
    assert abs(g - G) / abs(G) <= 1e-9
    # Q_back is cancellation-limited at large x (PHYSICS.md §1.2) — honestly looser.
    assert abs(qb - Qb) / max(abs(Qb), 1e-30) <= 1e-5


@pytest.mark.parametrize("x", XS, ids=[f"x{x:g}" for x in XS])
def test_farfield_shape_vs_miepython(x):
    mu = np.cos(np.linspace(0.01, np.pi - 0.01, 200))
    S1, S2 = pw.scattering_amplitudes(M, x, np.arccos(mu))
    Ia = np.abs(S1) ** 2 + np.abs(S2) ** 2                  # our global x2 constant cancels
    iu = mp.i_unpolarized(M, x, mu, norm="bohren")
    assert np.max(np.abs(Ia / Ia.sum() - iu / iu.sum())) <= 1e-8
