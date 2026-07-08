# Master Plan

## 1. Objective

Deliver a validated, documented Python package (`mieinfo`) that computes the
**information radiation pattern** of a levitated dielectric microsphere and uses
it to **evaluate and optimize detection geometry**. This is a **general-purpose
instrument**, not a one-off Coriolis calculation: given a sphere, one or more
probe beams (each with its own wavelength, propagation axis, polarization, and
focusing), and one or more collection channels (direction, NA, scheme), it reports
where the information about a chosen perturbation lives and how efficiently a given
optical configuration captures it. The silica Coriolis setup is the first
validating instance and drives the defaults, but nothing in the core is specialized
to it.

Concretely, the package must answer, for an arbitrary configuration:

- For a chosen sensed parameter — any displacement direction `n̂` (x, y, z, or any
  linear combination), with the derivative abstraction left open for other
  parameters — how is the Fisher information about it distributed over scattering
  direction, **per probe beam**?
- **Which direction is "axial" vs "transverse" is a property of a beam, not of the
  sphere.** A given lab axis is measured well in the forward lobe of a beam
  transverse to it and poorly in the forward lobe of a beam collinear with it. The
  tool computes this per beam and **combines information across channels** (Fisher
  information from independent channels adds).
- Given the collection optics (NA, forward / backward / split, homodyne /
  heterodyne / self-homodyne / imaging), what detection efficiency `η` (collected
  info / total available info) does each channel and each channel-set achieve, and
  how close does that put the measurement to the imprecision–backaction limit?
- What is the optimal configuration under real constraints, and how sensitive is it
  to sphere size, refractive index, and beam waist?
- Do the predictions agree with published detection-efficiency and
  imprecision-noise numbers where those exist?

## 2. Scope boundaries

**In scope:** classical/semiclassical Fisher-information treatment of coherent
optical position detection; plane-wave Mie and generalized Lorenz–Mie theory
(GLMT) for arbitrarily-focused probe beams; **detection wavelength decoupled from
trap wavelength** (the probe need not be the trap beam); **multiple simultaneous
probe/collection channels** along arbitrary axes, with Fisher information combined
across them; displacement derivatives of the scattered field along an arbitrary
direction; forward, backward, and split collection; collection-optics and
detection-scheme modelling (including imaging/self-homodyne of a dedicated probe);
optimization over geometry and channel sets; literature database and comparison.

**Out of scope (unless a later phase adds it):** full cavity-QED / master-equation
dynamics; feedback-loop and control design; thermal-force and gas-damping models;
manufacturing of optics; the xenon setup (this project is the silica setup —
"microsphere" without qualifier). Non-translational parameters (orientation of an
anisotropic/inhomogeneous particle) are a **stretch goal**: the parameter-derivative
interface is kept general so they can be added, but the homogeneous-sphere model
carries no orientation signal, so v1 implements translation only.

## 3. Physical regime (silica setup)

- Material: fused silica (medium = vacuum, so `m = n_particle`), `n ≈ 1.4496` at
  1064 nm and `n ≈ 1.4607` at 532 nm. Complex/dispersive index handled generally.
- **Two wavelengths.** Trap: `1064 nm` (holds the sphere, not used for detection).
  Detection: `532 nm` imaging beams. **All information/detection physics is at the
  detection wavelength.** The size parameter that sets the scattering is
  `x = 2πa n_med/λ_det`. Both wavelengths are parameters; defaults reflect this
  apparatus.
- Sphere radius: default sweep `a ∈ [3, 20] µm`. At `λ_det = 532 nm` this is
  `x ≈ 35–236`, `n_max ≈ 51–264` (Wiscombe `x + 4.05 x^{1/3} + 2`; see
  `PHYSICS.md §7`). The code must also cover Rayleigh→Mie for generality and for
  the dipole-limit validation. **Empirically, across 3–20 µm the collection
  efficiencies are nearly flat** — the sharp Mie-resonance size sensitivity is a
  sub-micron effect, so in this range size is a weak design knob (slightly favoring
  larger spheres for transverse readout).
- The probe beams are **focused beams**, not plane waves. Plane-wave Mie is the
  validation baseline and, for a weakly-focused/collimated probe over a small
  sphere, a good approximation via the phase-gradient derivative; GLMT (Richards–
  Wolf focus for high-NA probes) is the physical model. Because `x` reaches ~236,
  GLMT translation-addition at `n_max ~ 264` is the cost driver at the large end —
  see `PHYSICS.md §7` and `ARCHITECTURE.md`.

## 4. Phases

Phases are logical, not strictly sequential — the DAG in `ORCHESTRATOR.md` shows
what can run concurrently.

**Phase 0 — Foundation (orchestrator + W1a).**
Repo scaffold, CI, dependency lock, the special-function + plane-wave Mie core
promoted from `prototype/`, and the frozen interface contracts. Nothing else
starts until contracts are published.

**Phase 1 — Forward model (W1).**
GLMT beam-shape coefficients, focused-beam scattered field, translation-addition
machinery for arbitrary sphere position, and analytic displacement derivatives
`∂E_s/∂r_j`.

**Phase 2 — Information & detection (W2), concurrent with Phase 1.**
Fisher-information engine and information-pattern computation; detection-scheme
and collection-optics models; efficiency, imprecision PSD, backaction, SQL. Built
against the frozen `FieldProvider` contract and the plane-wave reference — does
not wait for GLMT.

**Phase 3 — Literature (W3), concurrent from Phase 0.**
Structured database of published detection geometries and results; curated
benchmark set with enough parameters to reproduce numerically.

**Phase 4 — Validation, optimization, comparison (W4).**
Validation harness (analytic limits, convergence, golden regression); geometry
optimizer; simulation-vs-literature comparison. Harness scaffolding starts early
against contracts; quantitative runs follow module completion.

**Phase 5 — Synthesis (orchestrator + human gate).**
The concrete recommendation for the silica setup, uncertainty-quantified, with
the trade studies behind it. Human physics sign-off here.

## 5. Session map

```
                          ┌───────────────────────┐
                          │     ORCHESTRATOR      │  plan · contracts · integrate · review
                          └───────────┬───────────┘
        ┌───────────────┬─────────────┼──────────────┬──────────────────┐
        ▼               ▼             ▼              ▼                  ▼
   ┌─────────┐    ┌───────────┐  ┌──────────┐  ┌──────────────┐   (human gate:
   │   W1    │    │    W2     │  │    W3    │  │      W4      │    Phase 5 sign-off)
   │scatter- │    │information│  │literature│  │validation +  │
   │ing eng. │    │+ detection│  │analysis  │  │optimization  │
   └────┬────┘    └─────┬─────┘  └────┬─────┘  └──────┬───────┘
   W1a plane-wave  W2a Fisher/    W3a harvest   W4a validation harness
   W1b GLMT/BSC    info-pattern   + database    W4b optimizer
   W1c displacement W2b detection  W3b benchmark W4c sim-vs-lit comparison
       derivatives     schemes         set
```

Dependency essentials (full DAG in `ORCHESTRATOR.md`):
- Everything depends on **W1a** (plane-wave core + special functions).
- **W2** depends on the `FieldProvider` *contract*, not the GLMT *implementation*
  — it develops against the plane-wave reference. This is the main concurrency win.
- **W3** is independent; starts immediately.
- **W4a** depends on contracts only (scaffold early). **W4b/W4c** need W1+W2
  implementations and W3's benchmark set.

## 6. Milestones (Definition of Done at each)

- **M0 Foundation:** `pip install -e .` works; CI green; `mieinfo.mie.plane_wave`
  reproduces `golden_values.json` to tolerance (`VALIDATION.md §2`); contracts
  published and frozen. Gate: orchestrator.
- **M1 Forward model:** GLMT scattered field reduces to plane-wave Mie in the
  wide-waist limit and to the dipole field in the small-`x` limit
  (`VALIDATION.md §3`); `∂E_s/∂r_j` matches finite-difference to `< 1e-6`
  relative. Gate: W4a + orchestrator.
- **M2 Information engine:** information pattern for a dipole reproduces
  Tebbenjohanns 2019 (transverse ⟂ axial structure; backscatter-optimal axial);
  `η(NA)` curves reproduce the seeded plane-wave numbers; energy/reciprocity
  checks pass. Gate: W4a + orchestrator.
- **M3 Literature base:** ≥ 15 experiments in the database with the required
  fields; ≥ 4 fully-parameterized benchmarks. Gate: orchestrator.
- **M4 Optimization + comparison:** optimizer returns a ranked configuration set
  for the silica setup with sensitivity analysis; predicted vs reported
  efficiency agrees within stated uncertainty for the benchmark set, discrepancies
  explained. Gate: W4 + orchestrator.
- **M5 Synthesis:** `docs/recommendation.md` with the configuration, the physics,
  the uncertainty, and the trade studies. Gate: human.

## 7. Success criteria

1. Every physics claim is backed by a passing validation gate.
2. The package reproduces at least two independent published results
   (Tebbenjohanns dipole limit; one Maurer-2023 Mie-regime number or one
   experimental efficiency from the benchmark set).
3. A concrete, defensible detection-geometry recommendation for the silica setup,
   with quantified sensitivity to the parameters that are actually uncertain.
4. Reproducible end to end: fresh clone → `make validate` → all gates green.

## 8. Apparatus facts (confirmed) and remaining open items

Confirmed apparatus facts — these are fixtures the defaults are built on, not open
questions (recorded as `FACT` entries in `STATUS.md`):

- **F1 Radius range:** `a ∈ [3, 20] µm`. At `λ_det = 532 nm`, `x ≈ 35–236`. Across
  this range detection efficiency is nearly flat (§3); resonance size-sensitivity is
  a sub-micron effect.
- **F2 Backward collection is available.** Both forward and backward ports are
  modelled and optimized. This is decisive: axial (beam-collinear) information is
  backward-weighted, so the backward port raises that axis's `η` from ~6% to
  ~65–70% (see `PHYSICS.md §8`).
- **F3 Detection is a dedicated 532 nm probe, not the trap.** Two imaging beams:
  (i) collinear with the 1064 nm trap, imaging the two directions transverse to it;
  (ii) a second beam propagating along a horizontal axis, sensitive to the remaining
  pair. Each lab axis is transverse to at least one beam. Self-homodyne, if used, is
  of the **probe**, not the trap. Model as a **multi-channel** detector; combine
  Fisher information across beams.
- **F4 Coriolis drive/readout:** the sphere's center of mass is driven around the
  horizontal loop N→E→S→W (CCW) and then CW; common mode is subtracted (drive-
  reversal isolates the velocity-odd Coriolis response from velocity-even
  backgrounds), and amplitudes are measured. This is **transverse (x, y) position
  sensing of a driven trajectory** — fully covered by the displacement machinery.
  It is one *application*; the tool is general (F-series defaults, not hardcoding).

Remaining genuinely-open items (defaults taken; confirm to tighten, non-blocking):

- A1: exact probe-beam waists / NA at the sphere (sets whether the plane-wave
  phase-gradient approximation suffices or GLMT-focus is required per beam).
  Default: provide both; treat the collinear beam as moderately focused and the
  horizontal beam similarly, and report the approximation error.
- A2: probe polarizations and exact propagation axes of the two beams in lab frame.
  Default: beam (i) along `z`, `x`-polarized; beam (ii) along `y`; parameterized so
  the real geometry drops in.
- A3: collection NA actually available at each port (forward/backward, each beam).
  Default: sweep NA up to 0.95 and report the curve; the recommendation reads off
  the achievable NA.
