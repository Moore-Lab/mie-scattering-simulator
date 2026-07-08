# Prototype (validated seed / reference oracle)

Worked-through physics that grounds the project. Treat as the plane-wave oracle:
`mieinfo.mie.plane_wave` must reproduce these numbers (`VALIDATION.md §2`).

## Two wavelengths

Trap = 1064 nm (holds the sphere, not used for detection). **Detection = 532 nm**
imaging beams. All information/detection physics is at 532 nm; the size parameter
is `x = 2πa/λ_det`. Silica index: n≈1.4607 @532, n≈1.4496 @1064.

## Files

- `mie_core.py` — from-scratch Bohren-Huffman plane-wave Mie (coefficients,
  efficiencies, `S1/S2`, angle functions). Agrees with `miepython` 3.2 to ~1e-12
  across `x = 0.30–29.5`, including absorbing spheres. Promote into
  `mieinfo/mie/{plane_wave,special}.py` at M0, preserving numerics.
- `information_pattern.py` — plane-wave information radiation pattern and `η(NA)`
  helpers (`cone_mask`, `collection_efficiency`, `info_density`). Uses the exact
  plane-wave displacement phase gradient (`PHYSICS.md §3.1`).
- `study_detection_532.py` — the operative study: 532 nm, a=3–20 µm, forward AND
  backward collection, with a grid-convergence check (angle-integrated g vs
  coefficient g, passes to ≤1e-10 up to x=236).
- `golden_values.json` — validated `Q_ext, Q_sca, Q_back, g` (engine regression;
  wavelength-independent). G-GOLD source.
- `information_pattern_results.json` — **operative** info-pattern golden set at
  532 nm, a=3/5/8/12/20 µm, forward+backward η at NA 0.5/0.8/0.95 for x and z.
- `information_pattern_results_1064.json` — retained 1064 nm reference set.

## Reproduce / self-validate

```
pip install numpy scipy miepython        # (mpmath optional, for ground-truth checks)
python mie_core.py               # smoke test
python validate.py               # SELF-VALIDATION gate: mie_core vs miepython (~1e-12),
                                 #   golden round-trip, info-pattern regression. Exit 0 = OK.
python validate.py --full        # also regenerates + diffs the 532 nm set (slower)
python study_detection_532.py    # operative 532 nm study -> information_pattern_results_532.json
```

`validate.py` is the seed's proof-of-correctness — run it before building on the
seed. It earned its keep immediately: it caught a real large-`x` bug in `mie_core`
(insufficient `D_n` downward-recurrence burn-in margin above `|mx|`), which was
invisible while the golden set stopped at `x=29.5` but corrupted every `a_n` by
~5e-4 once the 532 nm detection wavelength pushed `x` to 236. Fixed in
`logarithmic_derivative` (margin now `15 + 15√|mx|`); see `PHYSICS.md §1.1`.

## Headline results (silica, detection at 532 nm, vacuum)

Per beam, in that beam's own frame ("axial" = collinear with the beam):

- **Transverse → forward.** Transverse-axis info is forward-heavy (fwd frac
  0.90–0.94) but peaks off-axis (40–55°). η_x(NA0.8, forward) ≈ 0.52 (a=3 µm) →
  0.59 (a=20 µm). Backward is useless for transverse (η ≈ 0.02).
- **Axial → backward.** Beam-collinear-axis info is backward-weighted (peak
  155–159°). Forward collection is poor (η_z(NA0.8, fwd) ≈ 0.055–0.069) but backward
  gives η_z(NA0.8, bwd) ≈ 0.65–0.69. The backward port is decisive for that axis.
- **Size is a weak knob across 3–20 µm.** η varies smoothly and slightly; the sharp
  Mie-resonance size dependence is a sub-micron effect (gone by x≳35).

Consequence for a multi-beam setup: every lab axis should be transverse to ≥1 beam
(measured forward there) or read out in a beam's backward port. Fisher information
adds across independent channels — this is what the two-beam imaging layout buys.

These are plane-wave numbers in each beam frame. The GLMT engine shifts them; its
plane-wave limit must reproduce them and its dipole limit (x≪1) must reproduce
Tebbenjohanns 2019.
