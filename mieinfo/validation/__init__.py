"""Validation gates (VALIDATION.md). Time convention e^{-iωt}; SI units.

Reusable gate functions grouped by taxonomy (VALIDATION.md §1):

- `golden`   — G-GOLD: regression vs the seeded validated numbers
  (`data/golden/*.json`).
- `limits`   — G-LIMIT: analytic limits (energy / optical theorem, the dipole-limit
  information-pattern structure, and reciprocity/symmetry φ-parity).
- `convergence` — G-CONV: stability of efficiencies under `n_max → n_max+8` and of
  smooth solid-angle integrals under angular-grid doubling.

Each function returns a plain result (a bool pass flag and, where meaningful, the
observed max relative error) so both the pytest suite and `make validate` /
`mieinfo validate` can reuse them without re-deriving tolerances. Tolerances are
contract values (VALIDATION.md preamble); loosening one is a CONTRACT-CHANGE.
"""
from __future__ import annotations

from . import convergence, golden, limits

__all__ = ["golden", "limits", "convergence"]
