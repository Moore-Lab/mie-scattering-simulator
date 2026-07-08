# Literature

The seed bibliography and the extraction schema/method for the W3 track. W3 turns
this into a structured database of published detection geometries and results plus a
curated benchmark set that `W4c` (`compare.py`) confronts the simulation with under
the `G-LIT` gate. The frozen record shapes live in `INTERFACES.md §8`
(`Experiment`, `Benchmark`); this document is the reading list, the extraction rules,
and the acceptance counts.

---

## 1. Purpose and deliverables

W3 produces two data files, validated against the orchestrator-owned schema
(`mieinfo/literature/schema.py`, `INTERFACES.md §8`) and read by
`load_experiments` / `load_benchmarks`:

- `data/literature/experiments.json` — **≥ 15** `Experiment` records with the required
  fields and provenance (M3 count, `MASTER_PLAN.md §6`).
- `data/literature/benchmarks.json` — **≥ 4** `Benchmark` records, each with enough
  parameters to reproduce a `mieinfo` `RunConfig` numerically.

These feed the M4 success criterion of reproducing **≥ 2 independent published
results** (`MASTER_PLAN.md §7`, `VALIDATION.md §6`, §8): the Tebbenjohanns dipole case
plus at least one Mie-regime or experimental efficiency number.

## 2. Seed bibliography (anchors)

Tiered. **Verified** = the citation and the specific number/claim used elsewhere in
this repo have been confirmed against the source. **To-verify** = seeded from
`PHYSICS.md §9` but not yet checked; T3.1 resolves these (an unresolvable one becomes
an `EXTERNAL` blocker, `CONVENTIONS.md §5`).

**Primary theory anchors**
- Tebbenjohanns, Frimmer, Novotny, "Optimal position detection of a dipolar scatterer
  in a focused field," **PRA 100, 043821 (2019)** — arXiv:1907.12838. Dipole
  information pattern; backscatter-optimal axial detection; the imprecision–backaction
  (Heisenberg) limit. *The primary anchor and the dipole-limit benchmark
  (`VALIDATION.md §3`, §6).*
- Maurer, González-Ballestero, Romero-Isart, "Quantum theory of light interaction with
  a Lorenz–Mie particle: optical detection and 3D ground-state cooling," **PRA 108,
  033714 (2023)**; Erratum **PRA 109, 019901 (2024)**. The full-Mie generalization —
  the direct theoretical target and the Mie-regime benchmark source.

**Context / review**
- Gonzalez-Ballestero, Aspelmeyer, Novotny, Quidant, Romero-Isart, "Levitodynamics,"
  **Science 374, eabg3027 (2021)**. Field review; use for parameter ranges and to seed
  the snowball.

**Method references (not benchmarks; they ground the engine, `PHYSICS.md §1–2`)**
- Bohren & Huffman, *Absorption and Scattering of Light by Small Particles* (1983).
  Plane-wave Mie algorithm (`PHYSICS.md §1`).
- Gouesbet & Gréhan, *Generalized Lorenz–Mie Theories*, 2nd ed. (2017). GLMT/BSC
  (`PHYSICS.md §2.3`); confirm edition/chapter for the localized approximation and the
  `g_n` coefficients (also Doicu & Wriedt on the localized approximation).
- Richards & Wolf, **Proc. R. Soc. A 253, 358 (1959)**. High-NA aplanatic focal field
  for the Richards–Wolf beam (`PHYSICS.md §2.2`).

**Recent, on-topic — to-verify (T3.1)**
- arXiv:2409.00782, "Optimal displacement detection of arbitrarily-shaped levitated
  dielectric objects."
- arXiv:2512.17894, "Visualizing detection efficiency in optomechanical scattering."
- Pluchar et al., "Imaging-based quantum optomechanics," **PRL 135, 023601 (2025)**.

## 3. Extraction schema

Records validate against `INTERFACES.md §8` (frozen; changes go through
`CONTRACT-CHANGE`, `CONVENTIONS.md §5`). Extraction rules:

- **Provenance is mandatory.** Every `Experiment` carries `provenance` and every
  extracted number cites where in the paper it comes from (figure / equation / table /
  page). A number without provenance is not entered.
- **Every field is a value-with-provenance or an explicit `null`.** Do not infer or
  fill unstated parameters — a missing waist/NA is `null`, not a guess. `notes` records
  anything that constrains reproduction (e.g. "NA inferred from stated objective").
- **Units are SI**, field names ending `_m` are metres; convert reported values and
  note the original units in `notes` / `reported_imprecision` (which keeps the value +
  units *as reported*, per the schema).
- **Complex refractive index:** record per the `schema.py` encoding for
  `refractive_index: complex | None` (pydantic has no native complex — see
  `ARCHITECTURE.md §1`); convention `Im(m) ≥ 0` for absorption (`PHYSICS.md §0`).
- **Controlled vocabularies** align to the schema strings: `collection_direction ∈
  {forward, backward, split, cavity}`; `detection_scheme ∈ {homodyne, heterodyne,
  self-homodyne, split, imaging}`; `dof ∈ {x, y, z, rotation, …}`.

## 4. Harvest method

- **T3.1 — verify anchors + snowball.** Resolve every anchor's DOI/arXiv; confirm that
  each number this repo attributes to a source (e.g. the Tebbenjohanns backscatter
  claim, any Maurer efficiency) actually appears where claimed. Then snowball:
  forward/backward citations from the anchors, focused on levitated-optomechanics
  position-detection-efficiency and imprecision results. Record provenance as you go.
- **T3.2 — populate `experiments.json` (≥15).** One record per distinct experiment or
  clearly-parameterized theory result. Breadth over the schema fields (material, radius,
  wavelength, NA, direction, scheme, dof, reported efficiency/imprecision/backaction).
- Keep the harvest auditable: the snowball path (which anchor led to which record) is
  worth a line in `notes` so the DB is reproducible.

## 5. Benchmark curation

A `Benchmark` (`INTERFACES.md §8`) is an `Experiment` promoted to a *numerically
reproducible* case: a `run_config` that fully specifies a `mieinfo` `RunConfig`
(`ARCHITECTURE.md §5`), a `target_quantity` + `target_value` + `target_tolerance`, and
`target_provenance`.

- **Definition alignment is the crux.** The comparison must use the *same efficiency
  definition* as the source; any normalization conversion is documented and performed by
  `W4c`, not baked into the target (`VALIDATION.md §6`). "Passing" by redefining the
  quantity is not a pass.
- **≥ 4 benchmarks**, and the set must include (a) the **Tebbenjohanns-2019 dipole**
  case (also the `G-LIMIT` dipole gate, `VALIDATION.md §3`) and (b) **≥ 1 Mie-regime or
  experimental efficiency** (e.g. a Maurer-2023 number or a measured detection
  efficiency). These two are the independent-result requirement of `MASTER_PLAN.md §7`.
- A benchmark whose parameters cannot be fully pinned from the source is left as an
  `Experiment` (not promoted) or flagged `EXTERNAL` — do not invent parameters to make
  it runnable.

## 6. Gates and Definition of Done

- **M3 gate (orchestrator):** `experiments.json` has ≥ 15 records with the required
  fields and provenance; `benchmarks.json` has ≥ 4 definition-aligned, reproducible
  benchmarks including the two required independent cases; all anchors verified or the
  unresolved ones filed as `EXTERNAL`. (`MASTER_PLAN.md §6`, `STATUS.md`.)
- **Downstream (M4, `G-LIT`):** `W4c.compare_benchmark` predicts each benchmark's target
  within its `target_tolerance`, or the discrepancy is explained (regime, missing
  parameter, definition mismatch) in `discrepancy_note` (`VALIDATION.md §6`).
