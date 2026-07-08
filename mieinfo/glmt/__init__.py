"""Incident beams, beam-shape coefficients, and the FieldProvider seam.

At M0 the orchestrator seeds `beam.PlaneWave`, `scatter.FieldProvider` (protocol),
`scatter.PlaneWaveProvider`, and the phase-gradient `scatter.field_derivative` so the
info/detect tracks have a working provider from day one (ORCHESTRATOR.md §2.3). W1b/W1c
extend this package with GaussianParaxial/RichardsWolfFocus, the BSC methods, VSWF
translation-addition, and GLMTProvider. Time convention e^{-iωt}.
"""
