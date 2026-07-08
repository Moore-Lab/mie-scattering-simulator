# Architecture

Repo layout, package structure, tech stack, and the module import DAG. This is what
the orchestrator scaffolds in Phase 0 (`ORCHESTRATOR.md §2.1`) before any worker
starts. Module paths and public signatures are fixed by `INTERFACES.md` — this
document places those modules in the tree, states the dependency graph between them,
and pins the numerical/performance and provenance conventions that the contracts
assume but do not themselves specify.

Nothing here overrides `INTERFACES.md`. Where a type is referenced but not given a
frozen signature there (notably `RunConfig`), this document says where it lives and
flags it as an M0 contract item to finalize during bootstrapping — it does not
freeze a competing signature.

---

## 1. Tech stack and dependencies

- **Python 3.11.** Union syntax (`int | None`, `complex | None`) and
  `from __future__ import annotations` are used throughout; every module begins with
  the future-annotations import so the contract signatures in `INTERFACES.md`
  evaluate lazily.
- **Runtime dependencies (the `mieinfo` package):**
  - `numpy` — arrays, vectorized VSWF/field evaluation, `numpy.polynomial.legendre.leggauss` for the Gauss–Legendre `cosθ` grid.
  - `scipy` — special functions where a stable library form is preferable to a hand recurrence (spherical Bessel/Hankel, Wigner-`d`/associated-Legendre helpers), and `scipy.optimize` for the geometry search (`optimize/search.py`).
  - `pydantic` (v2) — `literature/schema.py` (`Experiment`, `Benchmark` are `BaseModel`s, `INTERFACES.md §8`). Note: pydantic has no native `complex` type, so `refractive_index: complex | None` needs a custom type/validator (store as `"n+kj"` or `[re, im]`; document the encoding in `schema.py` and `LITERATURE.md §3`).
  - `matplotlib` — the only plotting dependency; all plotting is confined to `viz/`.
- **The core must not import `miepython` at runtime.** `miepython` is a *validation
  oracle only* — used in tests (`G-GOLD`, `G-LIMIT` cross-checks) and to regenerate
  `data/golden/golden_values.json`. It is a dev/test dependency, never a dependency of
  `mieinfo.mie` or anything downstream.
- **Dev / test dependencies:** `pytest`, `pytest-xdist` (parallel test lanes),
  `hypothesis` (property tests for the symmetry/reciprocity gates, `VALIDATION.md §3`),
  `miepython` (oracle), `mypy`, `ruff`.
- **Determinism:** all dependencies are pinned in a committed lockfile
  (`ORCHESTRATOR.md §2.4`); any stochastic step seeds its RNG explicitly (`§5`).

## 2. Repository layout (the Phase 0 scaffold)

The orchestrator creates this skeleton with empty/stubbed modules, a mirrored
`tests/` tree, and `data/golden/` seeded from `prototype/`. Each module is annotated
with its owning track and the `INTERFACES.md` section that freezes its API.

```
mie_project/
├── pyproject.toml            # package metadata, pinned deps, entry point `mieinfo`
├── requirements.lock         # pinned lockfile (determinism)
├── Makefile                  # install|test|validate|optimize|compare|report (INTERFACES §10)
├── .github/workflows/ci.yml  # two CI lanes (VALIDATION §7)
├── README.md  MASTER_PLAN.md  PHYSICS.md  ARCHITECTURE.md  INTERFACES.md
├── CONVENTIONS.md  VALIDATION.md  LITERATURE.md  ORCHESTRATOR.md  STATUS.md
├── W1_scattering_engine.md  W2_information_detection.md
├── W3_literature_analysis.md  W4_validation_optimization.md   # worker briefs (see note)
├── prototype/                # validated seed / reference oracle (notes.md describes it)
│   ├── mie_core.py           # plane-wave Mie oracle (miepython ~1e-12 across x=0.3–236; D_n burn-in 15+15√|mx|, PHYSICS §1.1)
│   ├── information_pattern.py # plane-wave info pattern + η(NA) helpers
│   ├── study_detection_532.py # operative 532 nm study (a=3–20 µm, fwd+bwd)
│   ├── validate.py           # SELF-VALIDATION gate: run `python prototype/validate.py --full` (exit 0); ORCHESTRATOR §2 runs it first
│   ├── notes.md              # prototype README / reproduce + self-validate
│   ├── golden_values.json                 # G-GOLD engine source
│   ├── information_pattern_results.json    # operative 532 nm info-pattern golden set
│   └── information_pattern_results_1064.json  # 1064 nm info-pattern regression pin
├── data/
│   ├── golden/               # seeded from prototype/ at scaffold time
│   │   ├── golden_values.json
│   │   └── information_pattern_results.json
│   └── literature/           # W3 output
│       ├── experiments.json  # ≥15 Experiment records (INTERFACES §8)
│       └── benchmarks.json   # ≥4 Benchmark records
├── docs/
│   └── recommendation.md     # M5 synthesis (T5.1)
├── mieinfo/
│   ├── __init__.py
│   ├── types.py              # [orchestrator] Sphere, Medium, AngularGrid, VectorField, FieldDerivative (INTERFACES §1); RunConfig (§5 below)
│   ├── cli.py                # validate|optimize|compare|report (INTERFACES §10)
│   ├── mie/
│   │   ├── special.py        # [orch T0.2, W1a extends] Riccati-Bessel, D_n, ψ/χ/ξ, π_n/τ_n (PHYSICS §1.1, §1.3)
│   │   ├── plane_wave.py     # [orchestrator, T0.2] promoted mie_core: a_n,b_n, efficiencies, S1/S2 (PHYSICS §1)
│   │   └── vswf.py           # [W1a] VSWF near/far, angular functions (T1.1)
│   ├── glmt/
│   │   ├── beam.py           # [W1b] IncidentBeam, PlaneWave, GaussianParaxial, RichardsWolfFocus, lab_from_beam_frame (INTERFACES §2)
│   │   ├── bsc.py            # [W1b] BSC + bsc_quadrature|localized|angular_spectrum (INTERFACES §3)
│   │   ├── scatter.py        # [orch: FieldProvider protocol + PlaneWaveProvider (T0.3) · W1b: GLMTProvider (T1.7)] + field_derivative (INTERFACES §4)
│   │   ├── translation.py    # [W1c] VSWF translation-addition (T1.8)
│   │   └── derivatives.py    # [W1c] analytic dE/dr_j + finite-difference reference (T1.10/1.11)
│   ├── info/
│   │   ├── fisher.py         # [W2a] info_density, fisher_total (INTERFACES §5)
│   │   └── modes.py          # [W2a] combine_direction, InformationPattern (INTERFACES §5)
│   ├── detect/
│   │   ├── optics.py         # [W2b] cone_mask, apodization, grids (INTERFACES §6)
│   │   ├── schemes.py        # [W2b] apply_scheme; split/homodyne/self-homodyne/quadrant (INTERFACES §6)
│   │   └── metrics.py        # [W2b] S_imp, Gamma_ba, SQL distance; DetectionChannel, evaluate_channels (INTERFACES §6)
│   ├── optimize/
│   │   ├── objective.py      # [W4b] Constraints objective, channel-set eta_q_total (INTERFACES §7)
│   │   └── search.py         # [W4b] optimize_detection + sensitivity (INTERFACES §7)
│   ├── literature/
│   │   ├── schema.py         # [orchestrator] Experiment, Benchmark, load_* (INTERFACES §8)
│   │   └── compare.py        # [W4c] compare_benchmark (INTERFACES §9)
│   ├── validation/
│   │   ├── golden.py         # [W4a] G-GOLD (VALIDATION §2)
│   │   ├── limits.py         # [W4a] G-LIMIT (VALIDATION §3)
│   │   └── convergence.py    # [W4a] G-CONV (VALIDATION §5)
│   └── viz/                  # [W4] all plotting: patterns, η(NA) curves, comparison figures
│       └── plots.py
└── tests/                    # mirrors mieinfo/ ; each gate test lives beside the module it gates
    ├── mie/ glmt/ info/ detect/ optimize/ literature/ validation/
    └── conftest.py           # shared fixtures (golden loaders, canonical grids)
```

**Worker-brief location (convention).** The four worker briefs live in the repo root
(`W1_scattering_engine.md`, `W2_information_detection.md`, `W3_literature_analysis.md`,
`W4_validation_optimization.md`); `README.md` references them there. (Earlier `README`
drafts used a `sessions/` path; that has been reconciled to the root layout — if a
future maintainer prefers a `sessions/` directory, move all four together and update
`README.md` in the same change.)

## 3. Module import DAG (must stay acyclic)

Imports flow one way. `types.py` sits at the root and imports nothing from the
package. The `FieldProvider` protocol (`glmt/scatter.py`) is the one seam the whole
downstream stack is written against — `info/` and `detect/` depend on the *protocol*,
never on a concrete provider (they take the provider as an argument, `W2` brief). This
is what lets W2 develop against `PlaneWaveProvider` while W1 is still building
`GLMTProvider`.

```
types.py  ────────────────────────────────────────────────┐  (leaf; imported by all)
   ▲                                                         │
mie/special ─► mie/plane_wave ─► mie/vswf                    │
                    │                 │                       │
                    │            glmt/beam ─► glmt/bsc ─► glmt/scatter (FieldProvider)
                    │                              │            │  ▲
                    │                       glmt/translation ─► glmt/derivatives
                    │                                            │
   PlaneWaveProvider (in scatter.py) wraps mie/plane_wave  ◄────┘
                    │
   ┌────────────────┴───────── depends only on FieldProvider protocol + types ─────────┐
   ▼                                                                                     ▼
info/fisher ─► info/modes                                            detect/optics ─► detect/schemes ─► detect/metrics
                    │                                                                     │
                    └──────────────────────────┬──────────────────────────────────────────┘
                                                ▼
                                       optimize/objective ─► optimize/search
literature/schema ─► literature/compare  (compare also uses FieldProvider + detect/)
validation/{golden,limits,convergence}   (import the modules they gate; no one imports validation/)
viz/                                     (imports types + result objects; nothing imports viz/)
cli.py                                   (top: wires validate/optimize/compare/report)
```

Rules the DAG enforces:
- `info/` and `detect/` **must not** import `glmt/` or `mie/` concretes — only
  `types.py` and the `FieldProvider`/`field_derivative` contract.
- `validation/` and `viz/` are sinks: they import inward, nothing imports them.
- `mie/plane_wave.py` is the base everything ultimately rests on and is the only
  physics the orchestrator stewards directly (it is already validated).

## 4. Numerical and performance architecture (see `PHYSICS.md §7`)

Detection is at `λ_det = 532 nm`, so the operative size parameters are large
(`x ≈ 35–236`, `n_max ≈ 51–264`) and the angular grid must resolve a `~1/x` forward
lobe (up to `1400 × 420` at `a = 20 µm`). Design consequences:

- **`AngularGrid` is Gauss–Legendre in `cosθ`, uniform in `φ`,** carrying solid-angle
  weights (`INTERFACES.md §1`). Angle-integrated quantities (`g`, `Q`, `η`) use these
  weights; the grid-convergence check (`∫`-derived `g` vs coefficient `g`) is the
  canary and passes to `≤1e-10` up to `x=236`.
- **BSC caching is the main cost lever.** BSCs depend only on `(beam, sphere_center,
  n_max)` and are reused across the whole angular grid, across DOF (x/y/z), across
  collection direction, and by the derivative. Cache them on that key
  (`W1` brief). The translation-addition step at `n_max ~ 264` is the dominant cost.
- **Vectorize VSWF evaluation over the grid; precompute angular functions once per
  `(n_max, grid)`.** Recurrences: downward for `D_n`, upward for `ψ, χ`
  (`PHYSICS.md §1.1`) — do not regress this when refactoring the promoted core.
- **Optimizer inner loop touches a cached field with a geometry mask; it never
  re-solves the scattering.** Changing NA / direction / scheme / collection cone is a
  mask or overlap on a fixed `VectorField` + `FieldDerivative`. Only a change of
  sphere, beam, or `r_s` invalidates the cache. (`W4b`, `PHYSICS.md §7`.)
- **Size sweep may sample `a` coarsely at the top end** — efficiencies are nearly flat
  across `3–20 µm` (`PHYSICS.md §8`), so resolution there buys little.

## 5. Artifacts, provenance, and reproducibility

Every result artifact (JSON, figure, comparison row) must be reproducible from what it
carries. This is the property `W4` relies on ("every result artifact carries its
provenance … every figure stores its generating config").

- **`RunConfig` — the serializable spec of one run.** It fully specifies a simulation:
  sphere, per-beam `Medium`/wavelength, beam(s), channel set, sensed direction
  `n_hat`, angular grid + `n_max`, and constraints. It is referenced by the literature
  schema (`run_config` field, `INTERFACES.md §8`), by `compare_benchmark`
  (`INTERFACES.md §9`), and by the CLI (`mieinfo optimize <config>`, `INTERFACES.md
  §10`) — but `INTERFACES.md` does not yet freeze its fields. **It lives in
  `mieinfo/types.py` and MUST be given a frozen signature at M0** (an orchestrator
  contract item, versioned like the others). Until then, treat `run_config: dict` as a
  placeholder whose keys mirror the constructor arguments of the value objects above.
- **Provenance envelope.** Each artifact embeds: the `RunConfig` (or its hash), the
  package version + git SHA, the grid + `n_max` actually used, the tolerances applied,
  and a UTC timestamp. Figures embed (or write alongside) the `RunConfig` that
  generated them.
- **Determinism.** Pinned deps + lockfile; any RNG is explicitly seeded; tests are
  deterministic (`VALIDATION.md`, `README.md` non-negotiables). Fresh clone →
  `make validate` → all applicable gates green is the reproducibility contract
  (`MASTER_PLAN.md §7.4`).

## 6. CLI and CI wiring

- **CLI (`mieinfo/cli.py`) and `Makefile`** expose `validate | optimize | compare |
  report`, exactly per `INTERFACES.md §10`. `make validate` is the integration gate
  the orchestrator runs before every merge to `develop`.
- **Two CI lanes** per `VALIDATION.md §7`: a fast PR lane (unit tests, `G-GOLD`,
  `G-DERIV` on small `x`, symmetry property tests, lint/type — must be green to merge)
  and a nightly/full lane (`@slow` BSC-quadrature cross-checks, full `G-LIMIT`,
  `G-CONV` sweep, all `G-LIT`, full-`x` golden regression). The orchestrator runs the
  full lane at each milestone gate before tagging `main`.
