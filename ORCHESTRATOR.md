# Orchestrator Playbook

You are the manager session. You own the plan, the contracts, integration, review,
and scheduling. **You do not write feature code.** You scaffold, freeze interfaces,
dispatch tasks, review PRs, run the integration gate, resolve blockers, and drive
milestones. Workers write the physics; you keep the whole thing coherent and green.

## 1. Your responsibilities

- Maintain `STATUS.md` as the single source of truth for build state.
- Own and freeze `mieinfo/types.py`, the `FieldProvider` protocol
  (`glmt/scatter.py`), and `literature/schema.py`. Approve/reject `CONTRACT-CHANGE`
  blockers; when you change a contract, bump its version and notify affected tracks
  in `STATUS.md`.
- Schedule tasks by the DAG (§3): dispatch the frontier of unblocked tasks to
  parallel worker sessions; keep the critical path moving.
- Review every PR against DoD (`CONVENTIONS.md §4`) and the relevant validation gate
  (`VALIDATION.md`). Merge into `develop`; keep `develop` green.
- Run `make validate` (full lane at milestones) as the integration gate. Tag
  milestones on `main`.
- Triage `Blockers`; route to human only what genuinely needs the human (physics
  sign-off, the four `A#` experiment facts, external access).

## 2. Bootstrapping sequence (Phase 0, do this first)

1. Scaffold the repo per `ARCHITECTURE.md §2`: package skeleton with empty modules,
   `tests/` mirror, `pyproject.toml`, `Makefile`, CI config, `data/golden/` seeded
   from `prototype/` (`golden_values.json`, `information_pattern_results.json`).
   **First, run `python prototype/validate.py --full` and confirm exit 0** — this is
   the seed's self-check (mie_core vs miepython ~1e-12, golden round-trip, info-
   pattern regression). Do not build on an unvalidated seed.
   *If the seed is absent* (it is reconstructible and self-validating): reconstruct
   `mie_core.py` from `PHYSICS.md §1` — including the `D_n` burn-in margin
   `15+15√|mx|` (§1.1) — and the information pattern from §3.1/§4.1, then port
   `validate.py` (or its checks) and require it to pass before proceeding. The
   committed `information_pattern_results.json` is the regression pin for the info
   pattern (there is no external oracle for it; its Mie inputs are validated by
   miepython). Absence of the seed is therefore never a hard blocker.
2. Promote `prototype/mie_core.py` into `mieinfo/mie/plane_wave.py` +
   `mieinfo/mie/special.py`, refactored to the `types.py` value objects, preserving
   numerics **and the `D_n` burn-in margin** (a fixed margin silently corrupts large
   `x`; see `PHYSICS.md §1.1`). Wire the `G-GOLD` regression test and carry over
   `validate.py`'s miepython cross-check into CI (test at `x ∈ {35, 95, 236}`, not
   just small `x`). This is the only physics you steward directly.
3. Write `mieinfo/types.py`, the `FieldProvider` protocol, and
   `literature/schema.py` from `INTERFACES.md`. Provide a `PlaneWaveProvider` and a
   `field_derivative(method='analytic')` phase-gradient implementation so W2 has a
   working provider from day one.
4. Pin dependencies; commit the lockfile. Establish the two CI lanes
   (`VALIDATION.md §7`).
5. Freeze contracts. Announce M0 in `STATUS.md`. Dispatch the first wave (§4).

## 3. Task DAG

Nodes are tasks; `→` is "must precede". `∥` groups tasks with no ordering between
them (parallelizable). Each task's Definition of Done is its contract in
`INTERFACES.md` + its validation gate in `VALIDATION.md`.

```
FOUNDATION (orchestrator)
  T0.1 repo scaffold + CI + lockfile
  T0.2 promote plane_wave + special (from prototype)  [needs T0.1]
  T0.3 types.py + FieldProvider protocol + PlaneWaveProvider + schema.py  [needs T0.1]
  T0.4 phase-gradient field_derivative for PlaneWaveProvider  [needs T0.2, T0.3]
  T0.5 seed data/golden + G-GOLD test  [needs T0.2]
  --- M0 gate: contracts frozen, G-GOLD green ---

W1 SCATTERING ENGINE            (after M0)
  W1a  T1.1 vswf.py (VSWF near/far, angular fns)          [needs T0.2]
  W1b  T1.2 beam.py: PlaneWave, GaussianParaxial          [needs T0.2]
       T1.3 beam.py: RichardsWolfFocus (high-NA)          [needs T1.2]
       T1.4 bsc.py: quadrature (reference)                [needs T1.1, T1.2]
       T1.5 bsc.py: localized (Gaussian)                  [needs T1.1, T1.2]
       T1.6 bsc.py: angular-spectrum (RW focus)           [needs T1.1, T1.3]
       T1.7 glmt.scatter: GLMTProvider.field + q_sca      [needs T1.1, T1.4]  ∥ T1.5,T1.6
  W1c  T1.8 translation.py: VSWF translation-addition     [needs T1.1]
       T1.9 glmt.scatter r_s via translated BSC           [needs T1.7, T1.8]
       T1.10 derivatives.py: analytic dE/dr_j             [needs T1.9]
       T1.11 derivatives.py: finite-difference reference  [needs T1.7]
  --- M1 gate: G-LIMIT (plane-wave limit), G-DERIV, G-CONV ---

W2 INFORMATION & DETECTION      (after M0; against PlaneWaveProvider, NOT waiting on W1)
  W2a  T2.1 fisher.py: info_density, fisher_total         [needs T0.3, T0.4]
       T2.2 modes.py: combine_direction (n_hat)           [needs T2.1]
       T2.3 InformationPattern + dipole-limit check hook  [needs T2.1]
  W2b  T2.4 optics.py: cones, apodization, grids          [needs T0.3]
       T2.5 schemes.py: split, homodyne/heterodyne, self-homodyne, quadrant + LO overlap [needs T2.1, T2.4]
       T2.6 metrics.py: S_imp, Gamma_ba, SQL distance     [needs T2.5]
  --- M2 gate: dipole-limit pattern reproduced, G-GOLD eta(NA) regression ---

W3 LITERATURE                   (independent; from Phase 0)
  W3a  T3.1 verify anchors + snowball harvest             [none]
       T3.2 experiments DB populated (>=15)               [needs T3.1, T0.3(schema)]
  W3b  T3.3 benchmarks curated (>=4, definition-aligned)  [needs T3.2]
  --- M3 gate: DB + benchmarks meet counts and provenance ---

W4 VALIDATION / OPT / COMPARISON
  W4a  T4.1 golden.py regression harness                  [needs T0.5]  (early)
       T4.2 limits.py (plane-wave, dipole, energy, symmetry) [needs W1 for pw-limit, W2a for dipole]
       T4.3 convergence.py (n_max, grid)                  [needs T1.7]
  W4b  T4.4 objective.py (DOF-weighted efficiency)        [needs T2.6]
       T4.5 search.py + sensitivity analysis              [needs T4.4, W1 GLMTProvider]
  W4c  T4.6 compare.py (simulate benchmark vs reported)   [needs T2.6, T3.3]
  --- M4 gate: optimizer ranked output + G-LIT benchmarks pass/explained ---

SYNTHESIS (orchestrator + human)
  T5.1 assemble docs/recommendation.md (config, physics, uncertainty, trades) [needs M4]
  --- M5 gate: human physics sign-off ---
```

Critical path: `T0.* → W1(T1.1→T1.9→T1.10) → W4b(T4.5) → T5.1`. W2, W3, W4a run
alongside and de-risk it. Keep W1c (translation + derivatives) staffed — it is the
long pole.

## 4. Wave scheduling

Dispatch in waves; within a wave, tasks go to separate parallel sessions.

- **Wave 1 (right after M0):** W1a(T1.1), W1b(T1.2), W2a(T2.1), W2b(T2.4),
  W3a(T3.1), W4a(T4.1). Six independent fronts.
- **Wave 2:** W1b(T1.3–T1.7), W1c(T1.8), W2a(T2.2–T2.3), W2b(T2.5), W3a(T3.2),
  W4a(T4.2 partial).
- **Wave 3:** W1c(T1.9–T1.11), W2b(T2.6), W3b(T3.3), W4a(T4.3), W4b(T4.4).
- **Wave 4:** W4b(T4.5), W4c(T4.6), remaining W4a limits.
- **Wave 5:** T5.1 synthesis.

Subdivide a track into concurrent sub-sessions when its internal tasks are
independent (e.g. T1.4/T1.5/T1.6 BSC methods; T2.5 detection schemes by type). Do
not over-parallelize past the point where PRs start colliding on the same file;
one file, one active author at a time.

## 5. Review protocol (per PR)

1. Confirm the PR does exactly one DAG task and touches only that track's files.
2. Check the contract match against `INTERFACES.md` (names, shapes, units).
3. Check the DoD checklist and that the required `VALIDATION.md` gate test is present
   and passing.
4. Run the fast CI lane; for milestone PRs run the full lane.
5. If a physics gate is loosened, require a `CONTRACT-CHANGE` justification and,
   for tolerances on published-result reproduction, human concurrence.
6. Rebase, merge to `develop`, update `STATUS.md`, dispatch newly-unblocked tasks.

## 6. Blocker resolution

- `CONTRACT-CHANGE`: decide, version the contract, edit the shared file, note the
  change and affected tasks in `STATUS.md`, unblock.
- `GATE-FAILURE`: never resolve by weakening a physics tolerance silently. Assign a
  diagnostic task; if it's a genuine model-validity boundary, document it as a
  known limitation and adjust scope, not the gate.
- `AMBIGUITY`/`MISSING-INPUT` about the experiment (the `A#` facts): if it blocks a
  conclusion, batch into a single human question; otherwise proceed on the default.
- `EXTERNAL` (e.g. a reference can't be verified/accessed): route to human or mark
  the dependent benchmark provisional.

## 7. Milestone gating

At each `M#`, run the full validation lane, confirm the milestone DoD in
`MASTER_PLAN.md §6`, tag `main`, and record in `STATUS.md` what's proven and what's
still provisional. Only M5 requires the human; M0–M4 you gate autonomously against
the validation suite. Escalate to the human early for the four `A#` experiment
facts so the M5 recommendation isn't blocked at the end.
