"""Plane-wave Mie core (Bohren & Huffman 1983). Time convention e^{-iωt}.

Promoted from the validated prototype (prototype/mie_core.py, matches miepython to
~1e-12 across x = 0.3–236). This is the reference oracle every downstream track rests
on; preserve its numerics (PHYSICS.md §1).
"""
