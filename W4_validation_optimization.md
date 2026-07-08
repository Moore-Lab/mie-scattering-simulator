# Session W4 — Validation, Optimization & Comparison

You build the layer that proves the model is right, searches for the best detection
geometry, and confronts predictions with the literature. Your validation harness
starts early against contracts; your optimizer and comparison run once W1/W2/W3
land.

Read first: `VALIDATION.md` (your primary directive), `PHYSICS.md §4–8`,
`INTERFACES.md §5–9`, `MASTER_PLAN.md §7–8`, `CONVENTIONS.md`.

## Scope

- `mieinfo/validation/golden.py`, `limits.py`, `convergence.py` — the gates. [W4a]
- `mieinfo/optimize/objective.py`, `search.py` — DOF-weighted objective, geometry/
  size/scheme search, sensitivity analysis. [W4b]
- `mieinfo/literature/compare.py` — simulate a benchmark config, compare to
  reported. [W4c]
- `mieinfo/viz/` — pattern plots, `η(NA)` curves, comparison figures (all plotting
  lives here).

Do **not** edit the physics implementations you're validating; if a gate fails, file
a `GATE-FAILURE` blocker with the diagnostic — do not weaken tolerances to pass.

## Tasks (DAG ids)

- T4.1 (early) `golden.py`: regression vs `data/golden/*` (G-GOLD). Build against
  contracts before the physics is final.
- T4.2 `limits.py`: plane-wave limit of GLMT, dipole limit, energy/optical theorem,
  reciprocity/symmetry (G-LIMIT). Needs W1 for the plane-wave limit and W2a for the
  dipole-pattern check.
- T4.3 `convergence.py`: `n_max` and grid convergence harness (G-CONV); report the
  coarsest converged grid for the optimizer to use.
- T4.4 `objective.py`: the efficiency objective under `Constraints`
  (`INTERFACES.md §7`); handles arbitrary direction `n_hat`, optional sphere-size-
  as-design-variable, and **channel-set** objectives (`eta_q_total` over a set of
  `DetectionChannel`, `PHYSICS.md §4.5`).
- T4.5 `search.py`: search over NA / direction / scheme, **over sets of beams/ports**
  when `max_channels > 1`, and radius if allowed, to maximize the objective;
  sensitivity `dη/dparam` for radius, NA, waist, `m`. Uses the cached-field pattern
  so geometry candidates don't resolve. Expect the optimizer to discover that the
  best set covers `n_hat`'s components across complementary beams (each lab axis is
  transverse to some beam) rather than piling NA onto one port.
- T4.6 `compare.py`: `compare_benchmark` — build the `RunConfig` from a benchmark,
  simulate, convert to the source's efficiency definition or flag the mismatch,
  report within-tolerance (G-LIT).

## Subdivision

W4a (validation) is independent and starts first — one session. W4b (optimizer) and
W4c (comparison) are separate sessions gated on W1/W2 (opt) and W2/W3 (compare).

## Contracts you satisfy

`optimize_detection`→`OptResult`, `compare_benchmark`→`ComparisonResult`
(`INTERFACES.md §7, §9`), plus the CLI targets `mieinfo validate|optimize|compare|
report` (`INTERFACES.md §10`). Consume providers/patterns/detection results through
their contracts; take the provider as an argument (plane-wave now, GLMT later).

## Gates (Definition of Done)

- All of §2–6 in `VALIDATION.md` implemented and wired into `make validate`.
- M4: `optimize_detection` returns a ranked config set for the silica setup with
  sensitivity; ≥2 independent literature/analytic results reproduced via G-LIT
  (Tebbenjohanns dipole + one Mie/experimental), discrepancies explained.
- Every result artifact carries its provenance (`ARCHITECTURE.md §5`); every figure
  stores its generating config.

## Watch-outs

- Reproduce two independent published results — this is a hard success criterion
  (`MASTER_PLAN.md §7`), not optional.
- This is a **general tool**; the Coriolis mode (transverse x/y sensing of the
  driven loop) is one input, not the only target. Support optimizing any `n_hat` and
  any channel set. Report the trade studies: NA↔efficiency, forward↔backward↔split
  (backward is essential for beam-collinear axes — `PHYSICS.md §8`), self-homodyne
  (probe LO)↔optimal LO, single-beam↔two-beam channel sets, and size↔efficiency
  (note it is nearly flat across 3–20 µm, `PHYSICS.md §8`).
- Tag every conclusion with the facts/assumptions (`F#`/`A#`) it depends on; the M5
  recommendation must state which findings are robust vs assumption-sensitive.
- Do not pass a G-LIT benchmark by redefining the quantity — alignment must be
  explicit and documented.
