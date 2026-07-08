# mie-scattering-simulator — `mieinfo`

Compute the **information radiation pattern** of a levitated dielectric microsphere and
optimize the optical detection geometry for position sensing.

Given a sphere, one or more probe beams, and a sensed displacement direction `n̂`, `mieinfo`
reports where the Fisher information about that motion lives across scattering angle, how
efficiently a given collection geometry captures it (`η`), and which detection configuration
gets closest to the imprecision–backaction (Heisenberg) limit. Built for the Moore Lab
levitated-optomechanics setup (fused-silica sphere, 1064 nm trap, dedicated **532 nm** imaging
readout, two probe beams, forward + backward ports), but the core is general.

## Status

Milestones **M0–M4 complete; 233 tests pass.** The plane-wave engine is validated against
`miepython` to 1e-12 (size parameter `x` up to 236); the information/detection pipeline
reproduces the Tebbenjohanns-2019 dipole result and a Maurer-2023 Mie-regime efficiency.
See **[`docs/recommendation.md`](docs/recommendation.md)** for the detection-geometry
recommendation and its honest validation ledger. Remaining: the full analytic VSWF
translation-addition (an optional speed path; the validated displacement path is
finite-difference), and human physics sign-off (M5).

## Install

```bash
pip install -e ".[dev]"     # runtime: numpy, scipy, pydantic; dev: miepython, matplotlib, pytest
```

Python ≥ 3.11.

## Command line

```bash
mieinfo validate            # run the validation gates (G-GOLD / G-LIMIT / G-CONV); nonzero exit on failure
mieinfo optimize            # rank detection geometries for the silica setup   -> results/optimize_result.json
mieinfo compare             # simulate the literature benchmarks (G-LIT)        -> results/comparison.json
mieinfo report              # optimize + compare + pointer to the recommendation
```

## Verifying the results

The engine is checked four independent ways — run any of these yourself:

```bash
python prototype/validate.py --full   # engine vs miepython (independent Mie impl): ~1e-12 across x=0.3-236
mieinfo validate                      # physics gates: optical theorem, dipole limit, phi-parity, convergence
mieinfo compare                       # reproduces published detection efficiencies (Tebbenjohanns 2019, Maurer 2023)
pytest -q                             # full test suite (234 tests)
```

Beyond the numerical checks, the physics is **hand-verifiable**: the information density is the
intensity reweighted by `(k̂ − ŝ)²` (PHYSICS.md §4.1) — so transverse-motion info ∝ `sin²θ`, axial ∝
`(1 − cosθ)²`. The code reproduces those analytic weights to machine precision, and the resulting
pattern (transverse peaks off-axis, axial in backscatter) is the Tebbenjohanns-2019 result. Make the
plots for any configuration with `mieinfo.viz.plot_information_pattern` / `plot_eta_vs_na`
(see `docs/figures/`).

## Python API

```python
import numpy as np
from mieinfo.types import Sphere, Medium, AngularGrid
from mieinfo.glmt.beam import PlaneWave
from mieinfo.glmt.scatter import PlaneWaveProvider
from mieinfo.info.modes import information_pattern
from mieinfo.detect.optics import CollectionGeometry, collection_efficiency

med = Medium(n=1.0, wavelength_vacuum_m=532e-9)     # 532 nm detection light, vacuum
sphere = Sphere(radius_m=5e-6, m=1.46 + 0j)         # fused silica, a = 5 µm
beam = PlaneWave(med)
grid = AngularGrid.full_sphere(500, 120)
provider = PlaneWaveProvider()

# Information pattern for axial (z) motion; backward-port efficiency at NA 0.8
pattern = information_pattern(provider, grid, sphere, beam, np.zeros(3), n_hat=[0, 0, 1])
eta = collection_efficiency(pattern, CollectionGeometry(direction="backward", NA=0.8))
print(eta)      # ~0.65 — axial-motion information is backward-weighted
```

Rank geometries / channel sets with `mieinfo.optimize.optimize_detection`; confront predictions
with published results via `mieinfo.literature.compare.compare_benchmark`.

## Key result

**Information ≠ intensity.** The scattered *intensity* is forward-peaked, but the *information*
about motion is reweighted by `(k̂ − ŝ)²`: transverse-motion info is forward/off-axis
(η ≈ 0.8 at NA 0.95), while axial-motion info is **backward-weighted** (η ≈ 0.72 backward vs
~6 % forward — the backward port is decisive for the beam-collinear axis). A two-beam layout
covers all three lab axes and Fisher information adds across beams. Trade studies in
[`docs/recommendation.md`](docs/recommendation.md).

## What's inside

| Package | What |
|---|---|
| `mieinfo.mie` | plane-wave Mie (Bohren–Huffman) + far-field VSWFs |
| `mieinfo.glmt` | focused beams, beam-shape coefficients, GLMT field provider, displacement derivatives |
| `mieinfo.info` | Fisher-information density, arbitrary-direction combination |
| `mieinfo.detect` | collection optics, detection schemes, imprecision/backaction/SQL metrics, multi-channel |
| `mieinfo.optimize` | detection-geometry optimizer + sensitivity analysis |
| `mieinfo.literature` | experiment/benchmark schema, database, simulation-vs-literature comparison |
| `mieinfo.validation` | G-GOLD / G-LIMIT / G-CONV gate harness |
| `mieinfo.viz` | pattern, η(NA), and comparison figures |

`prototype/` holds the validated reference-oracle seed — `python prototype/validate.py --full`
re-checks the engine against `miepython`, the golden values, and the information-pattern regression.

## Design docs (for contributors)

This package was built against a fixed set of design documents — the physics ground truth,
frozen module contracts, and validation gates. Read these before changing physics or interfaces:

- [`MASTER_PLAN.md`](MASTER_PLAN.md) — scope, phases, milestones
- [`PHYSICS.md`](PHYSICS.md) — theory ground truth (equations, references, golden values)
- [`ARCHITECTURE.md`](ARCHITECTURE.md) — repo layout, module DAG, numerics
- [`INTERFACES.md`](INTERFACES.md) — frozen module contracts
- [`VALIDATION.md`](VALIDATION.md) — the correctness gates
- [`CONVENTIONS.md`](CONVENTIONS.md), [`LITERATURE.md`](LITERATURE.md), [`ORCHESTRATOR.md`](ORCHESTRATOR.md), and `W1–W4_*.md` — workflow, bibliography, per-track build briefs
- [`STATUS.md`](STATUS.md) — current build state and decisions log

## License / contact

Moore Lab, Yale University. License: TBD.
