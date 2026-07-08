# Session W3 — Literature Harvest, Database & Benchmarks

You build the empirical grounding for the project: a structured database of published
levitated-optomechanics detection geometries and results, and a curated benchmark set
that `W4c` confronts the simulation against under `G-LIT`. Your track is independent
of the physics implementation and **starts in Phase 0** — you do not wait on W1/W2.

Read first: `LITERATURE.md` (your primary directive), `PHYSICS.md §4, §9`,
`INTERFACES.md §8–9`, `CONVENTIONS.md`, `MASTER_PLAN.md §7–8`, `VALIDATION.md §6`.
Anchors: Tebbenjohanns 2019 (dipole), Maurer 2023 (Mie).

## Scope

- `data/literature/experiments.json` — the experiments database (`≥15`).
- `data/literature/benchmarks.json` — the curated benchmark set (`≥4`).

You **populate data**; you do **not** own the schema. `mieinfo/literature/schema.py`
(`Experiment`, `Benchmark`, `load_*`, `INTERFACES.md §8`) is orchestrator-owned —
consume it, and request any field change via `CONTRACT-CHANGE` (`CONVENTIONS.md §5`),
never by editing the schema or bending a record to fit. `literature/compare.py` is
W4c's, not yours. Do **not** touch `mie/`, `glmt/`, `info/`, `detect/`, `optimize/`,
`validation/`.

## Tasks (DAG ids from `ORCHESTRATOR.md §3`)

- **T3.1 verify anchors + snowball harvest.** Resolve every anchor DOI/arXiv
  (`LITERATURE.md §2`); confirm the specific claims this repo attributes to each source;
  snowball forward/backward from the anchors over the detection-efficiency /
  imprecision literature. Record provenance for every number. `[none]`
- **T3.2 experiments DB (≥15).** One record per distinct experiment or clearly
  parameterized theory result, across the schema fields, each with `provenance`.
  `[needs T3.1, T0.3 schema]`
- **T3.3 benchmarks (≥4, definition-aligned).** Promote the fully-parameterized cases to
  `Benchmark`s with `run_config`, `target_quantity/value/tolerance/provenance`. Must
  include the Tebbenjohanns dipole case and ≥1 Mie-regime/experimental efficiency.
  `[needs T3.2]`

## Subdivision

W3a = harvest + experiments DB (T3.1→T3.2); W3b = benchmark curation (T3.3), gated on
W3a. The snowball harvest can split across concurrent sub-sessions by sub-topic
(dipole/Rayleigh detection vs Mie-regime vs experimental-efficiency papers) once the
anchors are verified — merge into one `experiments.json` with no duplicate `key`s.

## Contracts you satisfy

Every record in `experiments.json` / `benchmarks.json` validates against
`schema.py` (`INTERFACES.md §8`) and loads via `load_experiments` /
`load_benchmarks`. `key` is unique. `run_config` in each `Benchmark` fully specifies a
`mieinfo` `RunConfig` (`ARCHITECTURE.md §5`) so `W4c.compare_benchmark`
(`INTERFACES.md §9`) can build and run it.

## Gates (Definition of Done)

- **M3:** ≥15 experiments with the required fields and provenance; ≥4 definition-aligned,
  reproducible benchmarks including (a) the Tebbenjohanns-2019 dipole case and (b) ≥1
  Mie-regime or experimental efficiency; all anchors verified or the unresolved ones
  filed as `EXTERNAL`. (`MASTER_PLAN.md §6`, `LITERATURE.md §6`.)
- Feeds **G-LIT** (`VALIDATION.md §6`): the benchmark set is what `W4c` must reproduce
  within tolerance or explain — two independent results is a hard M4 success criterion
  (`MASTER_PLAN.md §7`).

## Watch-outs

- **Provenance for every number** — figure/equation/table/page. No provenance, no entry
  (`LITERATURE.md §3`).
- **Do not invent parameters.** An unstated waist/NA/index is an explicit `null`, not a
  guess; a case that can't be fully pinned stays an `Experiment` (unpromoted) or is
  flagged `EXTERNAL`. Reproducibility of the DB depends on this.
- **Definition alignment is not yours to fake.** Record the target *as the source
  defines it*; the normalization/conversion to `mieinfo`'s definition is W4c's
  documented step (`VALIDATION.md §6`). Never pre-convert a number to make a benchmark
  "pass."
- **Conventions:** SI units (`_m`), complex refractive index with `Im(m) ≥ 0`
  (`PHYSICS.md §0`) in the `schema.py` complex encoding (`ARCHITECTURE.md §1`);
  controlled vocab for `collection_direction` / `detection_scheme` / `dof`
  (`LITERATURE.md §3`).
- **To-verify anchors** (arXiv:2409.00782, arXiv:2512.17894, Pluchar PRL 135 023601
  2025) must be confirmed before use; if a reference can't be accessed, file `EXTERNAL`
  rather than entering an unverified record.
