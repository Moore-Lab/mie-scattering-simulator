"""Interactive browser viewer for mieinfo information-radiation patterns.

Run with ``python -m mieinfo.webapp`` (or ``mieinfo serve``) and open the printed URL.
A small Flask app serves a Plotly 3D viewer with a parameter panel: selecting a
configuration loads a cached pattern from ``results/patterns/`` if it exists, otherwise
computes it with mieinfo and caches it. A sidebar lists every cached calculation.
"""
