# Conventions

Git workflow, status protocol, Definition of Done, and the blocker-escalation
protocol that keeps parallel sessions from colliding or guessing. Every worker brief
says "read `CONVENTIONS.md` first"; the orchestrator enforces these at review
(`ORCHESTRATOR.md §5`). Section numbers here are referenced by other documents —
`§4` (Definition of Done, cited by `ORCHESTRATOR.md §5` and the `VALIDATION.md` intro)
and `§5` (blocker / `CONTRACT-CHANGE`, cited by the `INTERFACES.md` and `VALIDATION.md`
intros) are load-bearing; keep them at these numbers.

---

## 1. Git workflow

- **Two long-lived branches.** `develop` is the integration branch (kept green by the
  orchestrator); `main` carries milestone tags (`M0…M5`) only. No one commits directly
  to either.
- **One feature branch per DAG task,** named `<track>/<task-id>-<slug>`, e.g.
  `w1a/T1.1-vswf`, `w2b/T2.5-schemes`, `w4a/T4.1-golden`. A branch does exactly one
  DAG task (`ORCHESTRATOR.md §5.1`).
- **One PR per task, into `develop`.** The worker opens it; the orchestrator reviews
  against DoD (`§4`) and the relevant `VALIDATION.md` gate, rebases, and merges. Rebase
  (not merge-commit) to keep history linear.
- **One file, one active author at a time** (`ORCHESTRATOR.md §4`). A PR touches only
  its track's files (`§3`). If two ready tasks would edit the same file, they are
  serialized, not parallelized.
- Milestone tags are applied by the orchestrator on `main` after the full validation
  lane passes (`ORCHESTRATOR.md §7`).

## 2. Status protocol

`STATUS.md` is the single source of truth for build state; future sessions read it
first to resume. It uses the legend `☐ todo · ◐ in progress (owner) · ☑ done (PR) ·
⚠ blocked (blocker id)`.

- **Workers own their task rows.** On starting a task, set it `◐` with your session
  tag. On merge, set it `☑` with the PR reference. On hitting a blocker, set it `⚠`
  with the blocker id (`§5`) and add the entry under `Blockers`.
- **The orchestrator owns** the milestone rows, the `Blockers` triage, and the
  `Notes / decisions log`. It updates milestones only when the milestone DoD is met
  (`MASTER_PLAN.md §6`).
- Update `STATUS.md` **in the same PR** as the work it describes — a task is not done
  if its row still says `☐`/`◐`.

## 3. File ownership and track boundaries

- Each worker touches only the modules its brief lists. Cross-track edits are not
  allowed; you consume another track's output through its `INTERFACES.md` contract, not
  its live source.
- **Orchestrator-owned files** (a worker requests changes to these via `CONTRACT-CHANGE`,
  `§5`; it does not edit them): `mieinfo/types.py`, the `FieldProvider` protocol in
  `mieinfo/glmt/scatter.py`, `mieinfo/literature/schema.py`, `STATUS.md`, and the
  directive documents (`MASTER_PLAN.md`, `PHYSICS.md`, `ARCHITECTURE.md`,
  `INTERFACES.md`, `CONVENTIONS.md`, `VALIDATION.md`, `LITERATURE.md`,
  `ORCHESTRATOR.md`, the `W*` briefs).
- `mieinfo/glmt/scatter.py` is **co-owned**: the orchestrator owns the `FieldProvider`
  protocol and `PlaneWaveProvider` (T0.3); W1b writes `GLMTProvider` in the same file
  (T1.7), drawing its internals from the `glmt/` modules W1b owns (`bsc.py`,
  `translation.py`, `derivatives.py`). Coordinate edits to `scatter.py` through the
  orchestrator so the two owners don't collide.

## 4. Definition of Done

A task is **done** only when every box below is checked. This is the checklist the
orchestrator runs at review (`ORCHESTRATOR.md §5.3`); "scientific code that isn't
validated is a liability" (`VALIDATION.md`).

- [ ] Implements **exactly** its `INTERFACES.md` contract — names, argument order,
      return shapes, and units identical. No unilateral API changes (`§5`).
- [ ] The module docstring states the **`e^{-iωt}` time convention** (`PHYSICS.md §0`)
      and the units of every public quantity (SI unless the name carries `_um`/`_deg`).
- [ ] Any packing/index convention (e.g. BSC `[n, m]` packing) is documented in the
      module and identical across sibling implementations (`INTERFACES.md §3`).
- [ ] The task's **`VALIDATION.md` gate test is present, committed beside the module,
      and passing** — at least one gate per physics task (`VALIDATION.md §1`).
- [ ] Unit tests, `mypy`, and `ruff` pass on the **fast PR lane** (`VALIDATION.md §7`).
- [ ] Deterministic: seeded RNG, no wall-clock/order dependence.
- [ ] Result artifacts and figures carry their provenance / generating config
      (`ARCHITECTURE.md §5`).
- [ ] `STATUS.md` row updated to `☑` with the PR reference (`§2`), in the same PR.
- [ ] No re-derivation from memory: physics traces to a `PHYSICS.md` section, cited in
      the code where non-obvious.

## 5. Blockers and `CONTRACT-CHANGE`

When you are blocked, **escalate — do not guess, and do not fork the API.** File the
blocker; keep working the frontier of what is unblocked. A blocker is one entry
appended to `STATUS.md → Blockers` with an id `BLK-##`, its type, the owning task, and
exactly what is needed to clear it. Set the task row to `⚠ BLK-##`.

Blocker types (the orchestrator resolves per `ORCHESTRATOR.md §6`):

- **`CONTRACT-CHANGE`** — you need a frozen contract (`types.py`, `FieldProvider`,
  `schema.py`, or a signature in `INTERFACES.md`) to change. State the current
  signature, the needed one, and why. The orchestrator decides, versions the contract,
  edits the shared file, notes the change + affected tasks in `STATUS.md`, and unblocks.
  **Loosening any `VALIDATION.md` tolerance is a `CONTRACT-CHANGE`** (`VALIDATION.md`
  intro) — it is never done silently.
- **`GATE-FAILURE`** — a validation gate fails. Do **not** weaken the tolerance to pass.
  File the diagnostic; the orchestrator assigns a fix or, if it is a genuine
  model-validity boundary, documents it as a known limitation and adjusts scope, not the
  gate (`ORCHESTRATOR.md §6`).
- **`AMBIGUITY` / `MISSING-INPUT`** — an experiment fact (`A#`) or spec detail is
  undetermined. Proceed on the documented default (`MASTER_PLAN.md §8`) and file the
  blocker only if it gates a *conclusion*; the orchestrator batches genuine questions to
  the human.
- **`EXTERNAL`** — something outside the repo (a reference that can't be accessed/verified,
  external access). The orchestrator routes to the human or marks the dependent artifact
  provisional.

## 6. Coding standards

- **Time convention `e^{-iωt}`** stated in every module docstring; `Im(m) ≥ 0` for
  absorption (`PHYSICS.md §0`). A module that silently uses the other convention is a
  bug even if internally consistent.
- **Units:** SI unless a field/argument name says otherwise (`_m`, `_um`, `_deg`).
  Angles in radians internally; `_deg`-named quantities are degrees at the boundary only.
- **Numerics:** vectorized numpy over Python loops on the angular grid; the recurrence
  directions of `PHYSICS.md §1.1` (downward `D_n`, upward `ψ,χ`) are preserved; keep the
  promoted `mie.plane_wave` numerics bit-for-bit within the `G-GOLD` tolerance.
- **Types:** full type hints matching `INTERFACES.md`; `from __future__ import
  annotations` at the top of every module (`ARCHITECTURE.md §1`).
- **Providers are arguments, never imports:** `info/`, `detect/`, `optimize/`,
  `literature/compare.py` take a `FieldProvider` as a parameter; they never import a
  concrete provider (`W2` brief, `ARCHITECTURE.md §3`).

## 7. Commit and PR format

- **Commit subject:** `<task-id>: <imperative summary>` (e.g. `T1.1: VSWF near/far +
  angular functions`). Body explains non-obvious physics choices with a `PHYSICS.md`
  section reference.
- **PR body** states: (1) the single DAG task it closes; (2) the `INTERFACES.md`
  contract it satisfies; (3) the `VALIDATION.md` gate + the test that proves it; (4) the
  `§4` DoD checklist; (5) the `STATUS.md` rows updated. The orchestrator's review
  (`ORCHESTRATOR.md §5`) maps one-to-one onto these.
