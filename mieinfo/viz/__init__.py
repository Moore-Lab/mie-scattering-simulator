"""Visualization for mieinfo objects (ARCHITECTURE.md §5).

Matplotlib figures (Agg backend, headless) that take mieinfo objects and save
provenance-embedding figures to disk: the information radiation pattern, η_q(NA)
curves, and predicted-vs-reported literature comparisons.
"""
from __future__ import annotations

from .plots import (
    plot_comparison,
    plot_eta_vs_na,
    plot_information_pattern,
)

__all__ = [
    "plot_information_pattern",
    "plot_eta_vs_na",
    "plot_comparison",
]
