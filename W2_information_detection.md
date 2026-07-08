# Session W2 — Information & Detection

You build the physics that turns a scattered field and its displacement derivative
into a detection figure of merit: the information radiation pattern, per-DOF Fisher
information, detection-scheme models, and the efficiency / imprecision / backaction /
SQL metrics the optimizer ranks on.

**Key point: you do not wait for W1.** You develop entirely against the
`FieldProvider` contract, satisfied today by `PlaneWaveProvider` +
`field_derivative(method='analytic')` (phase-gradient, `PHYSICS.md §3.1`), both
delivered at M0. When `GLMTProvider` lands, your code consumes it unchanged.

Read first: `PHYSICS.md §4–5, §8`, `INTERFACES.md §1, §4–6`, `CONVENTIONS.md`,
`VALIDATION.md §2–3`. Anchor: Tebbenjohanns 2019 (dipole); Maurer 2023 (Mie).

## Scope

- `mieinfo/info/fisher.py`, `info/modes.py` — information density, `F_q(Ω)`,
  arbitrary-direction `n_hat` combination. [W2a]
- `mieinfo/detect/optics.py`, `detect/schemes.py`, `detect/metrics.py` — collection
  cones/apodization, detection schemes with LO overlaps, and the metrics. [W2b]

Do **not** touch `glmt/`, `mie/`, `optimize/`, `literature/`.

Detection is at **λ_det = 532 nm** (a dedicated imaging probe, **not** the 1064 nm
trap), and the apparatus has **two probe beams** along different lab axes
(`MASTER_PLAN.md F2–F3`). So your detection layer is **multi-channel** from the
start, and "self-homodyne" means the *probe* beam as LO, not the trap.

## Tasks (DAG ids)

- T2.1 `info_density` (=|dE/dq|²) and `fisher_total` (4π integral).
- T2.2 `combine_direction`: `dE/dq` for `q = r_s·n_hat` as the linear combination of
  Cartesian derivatives (`PHYSICS.md §4.4`) — no new field solves. Keep the
  parameter behind the general abstraction (translation is the only implemented
  family for a homogeneous sphere).
- T2.3 `InformationPattern` + the dipole-limit reproduction hook used by G-LIMIT.
- T2.4 `optics.py`: forward/backward/split cones, aplanatic `√cosθ` apodization,
  aperture truncation, solid-angle grids (Gauss–Legendre in `cosθ`). Grids must go
  fine enough for large `x` (forward lobe `~1/x`; see `PHYSICS.md §7`).
- T2.5 `schemes.py`: split/quadrant, balanced homodyne/heterodyne (optimal LO ∝
  `dE/dq` saturates `F_q(Ω)`; also a realistic Gaussian LO), and **self-homodyne
  with the probe beam as LO** (per beam, per direction). Each returns the achieved
  fraction of `F_q(Ω)` via mode overlap.
- **T2.7 multi-channel** (`DetectionChannel`, `evaluate_channels`, `PHYSICS.md
  §4.5`): evaluate each beam in its own frame, rotate `s_hat`/`n_hat` into the lab
  frame via `lab_from_beam_frame`, apply that channel's cone, and **combine — Fisher
  info adds across independent channels**. Flag (don't silently sum) channels that
  share a detector/LO.
- T2.6 `metrics.py`: imprecision PSD `S_imp ∝ 1/F_q(Ω)`, backaction `Γ_ba ∝ Q_sca`
  (from the provider), and SQL distance.

## Subdivision

W2a and W2b are independent after T2.1; run them as two sessions. Within W2b, the
detection schemes (split / homodyne / self-homodyne) can be split further once
`optics.py` lands.

## Contracts you satisfy

`info_density`, `fisher_total`, `combine_direction`, `InformationPattern`,
`CollectionGeometry`, `cone_mask`, `apply_scheme`, `DetectionResult` — exactly per
`INTERFACES.md §5–6`. Consume `FieldProvider`/`field_derivative`/`VectorField`/
`FieldDerivative`; never import a concrete provider directly — take it as an argument.

## Gates (Definition of Done)

- Dipole-limit (G-LIMIT, M2): with `PlaneWaveProvider` at `x=0.05`, the information
  pattern matches Tebbenjohanns-2019 structure — transverse pattern peaks off-axis,
  axial pattern backward-weighted, ideal-scheme `η(4π)=1`. Assert normalized pattern
  shape and the axial-vs-transverse forward-fraction ordering.
- G-GOLD `η(NA)`: reproduce `data/golden/information_pattern_results.json` (the
  **532 nm**, a=3–20 µm set) `η` at NA 0.5/0.8/0.95 for x and z, **forward and
  backward**. Anchors: transverse forward `η_x(NA0.8,fwd)`≈0.52→0.59; axial backward
  `η_z(NA0.8,bwd)`≈0.65–0.69 (vs ≈0.06 forward). Match to grid tolerance.
- Property tests: `η∈[0,1]`, `η(4π)=1` with optimal LO, symmetry under `φ`-parity.

## Watch-outs

- The information pattern is **not** the intensity pattern — it's reweighted by the
  derivative (`PHYSICS.md §4.1`). Don't collapse the two.
- `η_q` is the single figure of merit (imprecision–backaction product); report
  `S_imp`, `Γ_ba`, and SQL distance alongside it, but optimize `η_q`.
- Self-homodyne uses the **probe** (532 nm) as LO, not the trap. Model its shortfall
  vs the optimal LO per beam and per direction. Forward collection is poor for the
  axis collinear with a beam; that axis is recovered in that beam's **backward** port
  *or* as a transverse axis of the other beam — the multi-channel combination is
  what makes all three axes measurable. Make sure that falls out.
- Fisher info adds across independent channels; do not add correlated channels.
- Keep everything a function of the passed provider so the plane-wave→GLMT swap is
  a one-line change at the call site.
