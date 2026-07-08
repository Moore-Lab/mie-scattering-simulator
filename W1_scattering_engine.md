# Session W1 — Scattering Engine

You build the optical forward model: from an incident trapping beam and a sphere to
the far-field scattered vector field and its displacement derivatives. You are the
long pole; W2 and W4b consume your output through the `FieldProvider` contract.

Read first: `PHYSICS.md §1–3, §6–8`, `INTERFACES.md §1–4`, `CONVENTIONS.md`,
`VALIDATION.md §2–5`. Seed oracle: `prototype/mie_core.py` (validated to 1e-12).

## Scope

- `mieinfo/mie/special.py`, `mie/vswf.py` — special functions, VSWFs. [W1a]
- `mieinfo/glmt/beam.py`, `glmt/bsc.py`, `glmt/scatter.py` — beams, beam-shape
  coefficients, the GLMT `FieldProvider`. [W1b]
- `mieinfo/glmt/translation.py`, `glmt/derivatives.py` — sphere displacement and
  analytic field derivatives. [W1c]

Do **not** touch `info/`, `detect/`, `optimize/`, `literature/`. The `FieldProvider`
protocol and `types.py` are orchestrator-owned — request changes via
`CONTRACT-CHANGE`.

## Tasks (DAG ids from `ORCHESTRATOR.md §3`)

- T1.1 VSWFs and angular functions (reuse prototype `pi_tau`; extend to VSWF
  near/far). Foundation for everything below.
- T1.2 `PlaneWave`, `GaussianParaxial` beams; T1.3 `RichardsWolfFocus` (the high-NA
  focus the real trap needs — paraxial Gaussian is not valid at trap NA).
- T1.4 BSC by quadrature (reference, `@slow`); T1.5 localized (fast); T1.6
  angular-spectrum (RW focus). All must return the plane-wave weights for a plane
  wave (G-LIMIT). T1.4/1.5/1.6 are independent → parallel sub-sessions.
- T1.7 `GLMTProvider.field` + `q_sca` (scattered amplitudes = `a_n g_TM`,
  `b_n g_TE`; `a_n,b_n` from `mie.plane_wave`).
- T1.8 VSWF translation-addition (Cruzan/Stein coefficients); T1.9 displaced-sphere
  field via translated BSCs; T1.10 analytic `dE/dr_j`; T1.11 finite-difference
  reference for G-DERIV.

## Subdivision

W1a (special/vswf) → gate for W1b/W1c. Within W1b, split BSC methods (T1.4/5/6)
across sessions once `vswf.py` and `beam.py` land. W1c (translation + derivatives)
is one focused session on the critical path — keep it staffed.

## Contracts you satisfy

`IncidentBeam`, `BSC`, `bsc_*`, `FieldProvider` (as `GLMTProvider`),
`field_derivative(method='analytic'|'finite_difference')` — exactly per
`INTERFACES.md §2–4`. Same BSC packing across all three methods.

## Gates (Definition of Done)

- G-GOLD: `GLMTProvider` in the wide-waist/plane-wave limit reproduces the
  prototype `S1/S2` to ≤ 1e-6.
- G-LIMIT: all BSC paths → plane-wave weights (≤1e-6); dipole-limit field at
  `x=0.05` matches the analytic dipole field.
- G-DERIV: analytic vs finite-difference `dE/dr_j` ≤ 1e-6 over `|r_s| ≤ 0.3λ`.
- G-CONV: fields/efficiencies stable under `n_max+8` and grid doubling (≤1e-4).

## Watch-outs

- `e^{-iωt}` convention throughout (`PHYSICS.md §0`); state it in every docstring.
- Cache BSCs by `(beam, sphere_center, n_max)` — they're reused across the whole
  angular grid, across DOF, and by the derivative. This is the main cost lever
  (`PHYSICS.md §7`). Vectorize VSWF evaluation over the grid.
- Localized-approximation validity degrades at high NA — flag it; the RW +
  angular-spectrum path is the trustworthy one at trap focusing.
- Downward recurrence for `D_n`, upward for `ψ,χ` (prototype pattern) — don't
  regress the numerics when refactoring.
- Each beam carries its **own** `Medium`/wavelength. Detection is at **532 nm**
  (not the 1064 nm trap), so the operative size parameter is `x = 2πa/λ_det` and
  reaches **~236** at a=20 µm ⇒ `n_max ~ 264`. Recurrences, VSWF evaluation, and
  translation-addition must stay stable and vectorized at that order; the
  translation step at large `n_max` is the cost driver (`PHYSICS.md §7`). Test the
  engine at `x ∈ {35, 95, 236}`, not just small `x`.
