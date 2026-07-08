"""G-GOLD: mieinfo.mie.plane_wave reproduces the seeded golden efficiencies
(VALIDATION.md §2). The golden set was validated against miepython to ~1e-12; the
1e-9 tolerance leaves margin for refactor noise. Time convention e^{-iωt}."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from mieinfo.mie import plane_wave as pw

GOLDEN = Path(__file__).resolve().parents[2] / "data" / "golden" / "golden_values.json"
CASES = json.loads(GOLDEN.read_text())["cases"]


@pytest.mark.parametrize("c", CASES, ids=[f"x{c['x']:.2f}" for c in CASES])
def test_efficiencies_match_golden(c):
    m = complex(c["m_real"], c["m_imag"])
    qext, qsca, qback, g = pw.efficiencies(m, c["x"])
    got = {"Qext": qext, "Qsca": qsca, "Qback": qback, "g": g}
    for key, val in got.items():
        rel = abs(val - c[key]) / max(abs(c[key]), 1e-30)
        assert rel <= 1e-9, f"{key}: rel err {rel:.2e} at x={c['x']}"


def test_energy_conservation_nonabsorbing():
    """Real m ⇒ Q_ext == Q_sca (optical theorem, VALIDATION.md §3)."""
    for c in CASES:
        if c["m_imag"] == 0.0:
            qext, qsca, _, _ = pw.efficiencies(complex(c["m_real"], 0.0), c["x"])
            assert abs(qext - qsca) / qsca <= 1e-10
