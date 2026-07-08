# Mie Information-Pattern Optimization

Simulate how position information is encoded in light scattered by a levitated
microsphere, and evaluate/optimize detection geometry against that information
content. This is a **general-purpose instrument**; the silica microsphere setup
(Coriolis measurement) is the first validating instance and sets the defaults, but
nothing in the core is specialized to it.

**Core question:** for a given sphere (size, material), a set of probe beams (each
with its own wavelength, propagation axis, polarization, focusing), and a sensed
parameter (any displacement direction `n̂`), which collection geometry / detection
scheme / **combination of beams** captures the most Fisher information — i.e. gets
closest to the imprecision–backaction (Heisenberg) limit?

Two facts about this apparatus shape the defaults: detection uses a **532 nm
imaging probe, separate from the 1064 nm trap**, and there are **two probe beams**
along different axes, so the detector is **multi-channel** (Fisher information adds
across independent beams). A key consequence, already visible in the seeded physics:
a lab axis is measured well in the forward lobe of a beam *transverse* to it and
well in the *backward* lobe of a beam *collinear* with it — so the two-beam layout,
plus the available backward port, is what makes all three axes measurable. See
`PHYSICS.md §4.5, §8`.

This repository is a **build directive for Claude Code**, not finished software.
It is designed so that a small fleet of Claude Code sessions — one orchestrator
plus parallel workers — can develop the full simulation with minimal human
intervention. The physics has been worked through, a validated reference engine
is seeded, golden values are provided, and the work is decomposed into
parallelizable tracks with fixed interface contracts.

## Read order

1. `MASTER_PLAN.md` — scope, phases, the session map, milestones. Start here.
2. `PHYSICS.md` — the theory ground truth (equations, references, golden values).
   Every worker reads the sections relevant to its track. No re-derivation.
3. `ARCHITECTURE.md` — repo layout, package structure, tech stack, module DAG.
4. `INTERFACES.md` — the contracts between modules. **This is what makes parallel
   development possible.** Workers code against these, never against each other's
   live implementation.
5. `CONVENTIONS.md` — git workflow, Definition of Done, status protocol, how a
   session escalates a blocker instead of guessing.
6. `VALIDATION.md` — the correctness gates: analytic limits, golden values,
   convergence, literature benchmarks.
7. `LITERATURE.md` — seed bibliography and the extraction schema for the
   literature track.
8. `ORCHESTRATOR.md` — the manager session's playbook and the full task DAG.
9. `W*.md` — one brief per worker track. A session opens exactly one.

## Who runs what

| Session | Role | Opens |
|---|---|---|
| Orchestrator | Owns plan, contracts, integration, review, scheduling. Writes no feature code. | `ORCHESTRATOR.md` |
| W1 | Scattering engine (Mie + GLMT forward model). | `W1_scattering_engine.md` |
| W2 | Information & detection model. | `W2_information_detection.md` |
| W3 | Literature harvest, database, benchmarks. | `W3_literature_analysis.md` |
| W4 | Validation, optimization, simulation-vs-literature comparison. | `W4_validation_optimization.md` |

Each W-track is subdividable into concurrent sub-sessions (e.g. W1a/W1b/W1c) once
its internal interfaces are fixed; see the track brief for the split.

## Seed material

`prototype/` contains a **validated plane-wave Mie implementation**
(`mie_core.py`, agrees with `miepython` to ~1e-12 across x=0.3–236) and a working
information-pattern computation (`information_pattern.py`) that already reproduces
the qualitative result of Tebbenjohanns et al. (2019): axial-motion information is
backward-weighted. `golden_values.json` and `information_pattern_results.json` hold
validated reference numbers (the latter is the operative 532 nm, a=3–20 µm set).

**`prototype/validate.py` is the seed's self-validation gate** — run it first
(`python prototype/validate.py --full`, expect exit 0). It cross-checks mie_core
against `miepython`, round-trips the goldens, and regression-pins the info pattern.
It is not decoration: it already caught a real large-`x` bug in the seed (an
insufficient `D_n` recurrence burn-in margin that corrupted every `a_n` by ~5e-4
once the 532 nm wavelength pushed `x` past ~60; fixed, see `PHYSICS.md §1.1`).

The seed is **reconstructible and self-validating**: if `prototype/` is ever
missing, rebuild `mie_core` from `PHYSICS.md §1` and the info pattern from §3.1/§4.1,
then require `validate.py` to pass (`ORCHESTRATOR.md §2`). So its absence is not a
blocker — but when present, treat it as the reference oracle and the starting point
for `mie.plane_wave`; do not discard it.

## Non-negotiables

- Scientific code is worthless if wrong. Nothing is "done" until it passes its
  validation gate (`VALIDATION.md`). Golden values and analytic limits gate every
  physics module.
- Contracts in `INTERFACES.md` are frozen unless the orchestrator changes them.
  A worker that needs a contract change files a blocker; it does not fork the API.
- Determinism: pinned dependencies, seeded RNG, deterministic tests.
