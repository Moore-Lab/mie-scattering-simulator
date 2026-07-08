# Status

Single source of truth for build state. The orchestrator maintains this; every
session updates its task rows on each PR. Future sessions read this first to resume.

Legend: ☐ todo · ◐ in progress (owner) · ☑ done (PR) · ⚠ blocked (blocker id)

Last updated: 2026-07-08 — COMPLETE TO M4, on GitHub (Moore-Lab/mie-scattering-simulator, branch main). 234 tests green; user-facing README + embedded figures. Reviewed (6 findings fixed). M1 analytic translation-addition ATTEMPTED and REVERTED to an honest raise (a single-radius projection kernel failed independent verification); the validated FD displacement path is unchanged. Remaining: M5 human physics sign-off; a correct Cruzan/Stein translation kernel (optional speed path).

## Milestones

- ☑ M0 Foundation — contracts frozen, G-GOLD green (2026-07-07)
- ◐ M1 Forward model — GLMT plane-wave limit ✅ + G-DERIV ✅ + G-CONV ✅ (full analytic VSWF translation-addition PENDING; validated displacement = quadrature-BSC + FD/product-rule derivative)
- ☑ M2 Information engine — dipole-limit + G-GOLD η(NA) + energy/φ-parity gates pass; schemes/metrics/multichannel built
- ☑ M3 Literature base — 17 experiments + 5 benchmarks (counts + provenance met; some point-values pending human sign-off)
- ☑ M4 Optimization + comparison — ranked configs + G-LIT (4/5; ≥2 independent: Tebbenjohanns dipole + Maurer Mie-regime)
- ◐ M5 Synthesis — docs/recommendation.md written; AWAITING HUMAN physics sign-off (the only remaining gate)

## Tasks

FOUNDATION (orchestrator)
- ☑ T0.1 repo scaffold + CI + lockfile (pyproject, Makefile, .github/workflows/ci.yml, requirements.lock; `pip install -e .` works)
- ☑ T0.2 promote plane_wave + special (mieinfo/mie/{special,plane_wave}.py; numerics preserved incl. 15+15√|mx| burn-in)
- ☑ T0.3 types.py + FieldProvider protocol + PlaneWaveProvider + schema.py (mieinfo/types.py, glmt/scatter.py, glmt/beam.py, literature/schema.py)
- ☑ T0.4 phase-gradient field_derivative — analytic == finite-difference to machine precision (G-DERIV for the plane-wave provider)
- ☑ T0.5 seed data/golden + G-GOLD test (tests/mie/test_golden.py + test_vs_miepython.py; 26 tests green, G-GOLD 1.36e-12)

W1 SCATTERING ENGINE
- ☑ T1.1 vswf.py (PW reconstruction S1/S2 ≤1e-8 to x=236) · ☑ T1.2 beams (PlaneWave, GaussianParaxial) · ☑ T1.3 Richards-Wolf focus (Debye–Wolf)
- ☑ T1.4 BSC quadrature · ☑ T1.5 BSC localized · ☑ T1.6 BSC angular-spectrum (PW→PW weights, G-LIMIT)
- ☑ T1.7 GLMTProvider.field + q_sca (wide-waist→S1/S2 4.5e-7; dipole limit) · ◐ T1.8 translation-addition (axial scalar kernel only; full Cruzan/Stein PENDING)
- ☑ T1.9 displaced-sphere field (quadrature-BSC at center +r_s + far-field phase) · ◐ T1.10 analytic dE/dr_j (PW exact + focused product-rule; full analytic translation PENDING) · ☑ T1.11 finite-diff ref (G-DERIV ~2e-7)

W2 INFORMATION & DETECTION
- ☑ T2.1 fisher.py · ☑ T2.2 modes.py (n_hat) · ☑ T2.3 InformationPattern + dipole-limit gate
- ☑ T2.4 optics.py (cone_mask, collection_efficiency, aplanatic √cosθ apodization) · ☑ T2.5 schemes.py (split/quadrant, homodyne opt+Gaussian LO, self-homodyne) · ☑ T2.6 metrics.py (S_imp, Γ_ba, SQL) · ☑ T2.7 multi-channel (evaluate_channels; Fisher adds)

W3 LITERATURE
- ☑ T3.1 anchors verified + harvest · ☑ T3.2 experiments DB (17, web-verified) · ☑ T3.3 benchmarks (5, incl. Tebbenjohanns dipole + Magrini/Maurer efficiency) — some numbers pending human sign-off (W3 gaps)

W4 VALIDATION / OPT / COMPARISON
- ☑ T4.1 golden.py · ☑ T4.2 limits.py (energy, dipole, φ-parity; corrected the x/y symmetry claim) · ☑ T4.3 convergence.py
- ☑ T4.4 objective.py (cached-field, sensitivity) · ☑ T4.5 search.py (single + multi-channel; axis-naming fixed) · ☑ T4.6 compare.py (G-LIT 4/5)
- ☑ viz/plots.py (pattern, η(NA), comparison figures) · ☑ cli.py wired (validate|optimize|compare|report)

SYNTHESIS
- ☑ T5.1 docs/recommendation.md (M4 draft: config + physics + optimizer output + G-LIT + ledger; awaiting M5 human sign-off)

## Apparatus facts — confirmed (MASTER_PLAN.md §8)

- F1 Radius range a ∈ [3, 20] µm ⇒ x ≈ 35–236 at λ_det=532 nm. η nearly flat across
  the range (resonance size-sensitivity is sub-micron); size is a weak knob here.
- F2 Backward collection AVAILABLE. Decisive: beam-collinear (axial) info is
  backward-weighted, η goes ~6%→~65–70% via the backward port. Optimize both ports.
- F3 Detection = dedicated 532 nm probe, NOT the trap. Two imaging beams (one
  collinear with the 1064 trap imaging its transverse plane; one horizontal imaging
  the complementary pair). MULTI-CHANNEL; Fisher info adds across beams. Every lab
  axis is transverse to ≥1 beam.
- F4 Coriolis: COM driven around N→E→S→W (CCW) then CW, common mode subtracted,
  amplitudes measured. = transverse x/y sensing of a driven trajectory. One
  application of the general n_hat machinery; TOOL IS GENERAL, not Coriolis-specific.

## Remaining open items (defaults taken; non-blocking)

- A1 exact probe waists/NA at the sphere (sets plane-wave-approx vs GLMT-focus per
  beam). Default: provide both, report approximation error.
- A2 exact lab-frame axes + polarizations of the two beams. Default: beam-i along z
  (x-pol), beam-ii along y; parameterized.
- A3 collection NA available per port. Default: sweep to 0.95, report the curve.

## Blockers

(none yet)

## Notes / decisions log

- Reference engine (plane-wave Mie) validated against miepython to 1e-12 across
  x = 0.30–29.5 incl. absorbing spheres; seeded in prototype/ and data/golden/.
- Plane-wave information-pattern prototype already reproduces the qualitative
  Tebbenjohanns-2019 result (axial info backward-weighted). This is the M2
  dipole-limit target to hit quantitatively via GLMT.
- Operative info-pattern golden set is at λ_det=532 nm, a=3–20 µm, forward AND
  backward (prototype/information_pattern_results.json); 1064 nm set retained for
  reference. Grid-convergence verified to ≤1e-10 (angle-integrated g vs coeff g)
  up to x=236 on ~1400×420 grids.
- Headline: transverse → forward (η_x≈0.52–0.59 @NA0.8); axial → backward
  (η_z≈0.65–0.69 @NA0.8, vs ~0.06 forward). Multi-channel is how all three axes
  get covered — each axis transverse to ≥1 beam.
- Contracts (types.py, FieldProvider, schema.py) freeze at M0; changes require an
  orchestrator-approved CONTRACT-CHANGE.
- 2026-07-07: the four missing directive docs were written — ARCHITECTURE.md,
  CONVENTIONS.md, LITERATURE.md, W3_literature_analysis.md — and cross-reference-audited
  against the corpus. README's stale `sessions/W*.md` paths reconciled to the root layout
  (all four briefs live in repo root). The read-order in README.md is now complete.
- 2026-07-08: repo on GitHub (Moore-Lab/mie-scattering-simulator); user-facing README (build
  directive preserved as "design docs"); recommendation.md embeds the info-pattern / η(NA) /
  comparison figures.
- 2026-07-08: M1 analytic translation-addition ATTEMPT → REVERTED to honest raise. A high-effort
  agent built a scalar axial kernel via O'-frame projection that passed its OWN (narrow low-degree /
  small-radius) reconstruction self-test, but the orchestrator's INDEPENDENT check (field-point
  reconstruction at general radii/degrees; d=0 identity control passes at 1e-14) caught it failing
  3–14% — a single-radius fit, NOT a true radius-independent translation operator (same trap class
  as the first attempt). Nothing depends on it; `axial_translation_coefficients` reverted to raise
  NotImplementedError. M1 stays ◐ on the validated quadrature-BSC + FD path. Lesson: agent
  self-tests can be narrow/circular — independent adversarial verification is what caught BOTH
  attempts. A correct kernel needs the Cruzan/Stein Gaunt recurrence (+ M↔N mixing for the vector
  BSC translation); left for a future pass.
- 2026-07-08: FINAL PHYSICS REVIEW (3 agents: optimizer / compare / GLMT) → 6 findings, ALL
  FIXED (233 tests green): (1) [major] objective.py `waist` sensitivity was hardcoded 0.0 → now
  finite-differences finite-waist beams (verified dη/dwaist=3.43e5 for a GLMT Gaussian; NaN if the
  beam can't be rebuilt with a waist arg). (2) [major] compare.py threshold benchmark passed
  VACUOUSLY (predicted ≥ reported−tol with tol > target) → threshold now requires predicted ≥ the
  threshold value itself. (3) [major] translation.py `axial_translation_coefficients` was
  mathematically WRONG (origin-measure projection on the shifted-center basis, 7–302% error) and
  falsely claimed "validated" → now RAISES NotImplementedError (it was never on the validated path;
  quadrature-BSC recompute is the real displacement path). (4) [minor] multichannel n_max
  ranking/re-score mismatch → aligned to evaluate_channels' truncation. (5) [minor] compare
  upper-bound (Magrini) no longer counted as a numerical pass (documented definition gap).
  (6) [minor] dead Greek-symbol check removed. G-LIT still 4/5 (≥2 independent intact); all gates pass.
- 2026-07-08: WAVE 2 integrated — GLMT engine (W1b/c), optimizer (W4b), compare (W4c), viz, CLI.
  235→233 tests. optimize surfaced F_x≠F_y (η_x=0.795 vs η_y=0.747 at NA0.95); fixed the cosmetic
  channel-name collision (little-endian tobytes hex). docs/recommendation.md written (M4).
- 2026-07-07: WAVE 1 (4 parallel agents, disjoint files) integrated — 153 tests green.
  W2b detection (detect/schemes,metrics,channels + optics apodization): schemes model κ via
  single-mode LO overlap (optimal LO κ=1=collection_efficiency; gaussian/split are honestly
  PESSIMISTIC single-bucket modes → tiny κ on the parity-forbidden axis; self-homodyne modeled
  as a spatially-resolving camera = optimistic end). Multi-channel evaluate_channels: Fisher
  info ADDS across independent beams (F_AB=F_A+F_B to 1e-10); two-beam z+x example covers all
  3 axes. W3 literature: 17 web-verified experiments + 5 benchmarks; some ">="/figure-read
  point values need human sign-off (see W3 gaps). W4a validation: G-GOLD/G-LIMIT/G-CONV harness;
  CORRECTED the directive's false x/y displacement symmetry (F_x≠F_y for a Mie sphere;
  VALIDATION §3 updated) and flagged hard-cone η grid-quantization (~1e-3; VALIDATION §5 note).
  W1a vswf: far-field VSWF reconstruction of S1/S2 to ≤1e-8 across x=1..236 — the validated
  seam GLMTProvider will use (near-field M_mn/N_mn not needed for the far-field contract).
- 2026-07-07: Wave 1 started — W2a INFORMATION ENGINE + optics core built (`info/fisher.py`,
  `info/modes.py`, `detect/optics.py`). M2 dipole-limit gate PASSES (transverse forward-heavy,
  axial backward-weighted, info≠intensity, η(4π)=1). Full package pipeline (PlaneWaveProvider →
  field_derivative → combine_direction → info_density → collection_efficiency) reproduces the
  golden 532 info-pattern set: forward-fraction + peak-angle to MACHINE PRECISION (~1e-16), η to
  ~4e-3 — the residual is the golden's own 80-pt NA-interpolation artifact (prototype eta_curve);
  the package computes η exactly at NA and is the more accurate value. Committed as @slow
  `tests/info/test_golden_pattern.py`. 58 tests green (53 fast + 5 slow). REMAINING W2: T2.5
  detection schemes (split / homodyne / self-homodyne + LO overlap), T2.6 metrics (S_imp, Γ_ba,
  SQL distance), T2.7 multi-channel, and aplanatic apodization in optics. Possible future
  tightening: regenerate the golden η with exact-at-NA values to enable a 1e-6 η gate (deferred —
  would diverge from the prototype's self-consistency check in validate.py).
- 2026-07-07: ☑ M0 REACHED & CONTRACTS FROZEN. Full foundation built + adversarially
  reviewed (3-agent review: physics/contract/numerics — clean on physics & numerics, 2 minor
  fixed): (1) `literature/schema.py` `notes` finalized OPTIONAL (`notes: str = ""`) with
  `provenance` kept REQUIRED — INTERFACES.md §8 updated to match; (2) `PlaneWaveProvider.field`
  now raises on a non-+x beam-frame polarization instead of silently returning the x-pol field
  (lab polarization is the channel's job, PHYSICS §4.5). 40 tests green. Contracts (types.py,
  FieldProvider protocol, PlaneWaveProvider, field_derivative, literature/schema.py) are FROZEN;
  changes now require a CONTRACT-CHANGE (CONVENTIONS §5). NEXT: dispatch Wave 1
  (ORCHESTRATOR §4) — W1a T1.1 vswf, W1b T1.2 beams, W2a T2.1 fisher, W2b T2.4 optics,
  W3a T3.1 harvest, W4a T4.1 golden.
- 2026-07-07: M0 CORE BUILT & VERIFIED. `pip install -e .` works; `pytest` 26/26 green;
  `mieinfo validate` → G-GOLD 1.36e-12. Promoted engine (`mieinfo/mie/special.py`,
  `mieinfo/mie/plane_wave.py`) preserves prototype numerics incl. the `15+15√|mx|` burn-in
  margin and cross-checks miepython to 1e-12 across x=0.3–236 (`tests/mie/test_vs_miepython.py`,
  currently in the fast lane). `mieinfo/types.py` value objects done. REMAINING for M0: the
  `FieldProvider` protocol + `PlaneWaveProvider` (`glmt/scatter.py`), phase-gradient
  `field_derivative` (T0.4), `literature/schema.py` — then freeze contracts + announce M0.
- ✅ SEED PRESENT + VALIDATED (2026-07-07): the reference-oracle seed was delivered (files.zip)
  and merged into `prototype/`. `python prototype/validate.py --full` → ALL GATES PASSED:
  mie_core vs miepython 1.5e-12 across x=0.3–236 (incl. the large-x burn-in fix), golden
  round-trip 1.4e-12, 532 operative info-pattern 2.0e-12. Two delivery fixes applied during the
  merge: (i) validate.py Gate 4b false-positive (it relative-compared `g_relerr`, a ~1e-13
  convergence residual) excluded via a DIAG skip; (ii) `information_pattern_results_1064.json`
  (absent from the delivery) regenerated from the validated engine. PHYSICS.md §1.1 now carries
  the corrected `D_n` burn-in margin `15+15√|mx|`; ORCHESTRATOR §2 bootstrap runs validate.py
  first. T0.2/T0.5 unblocked. Pre-merge root backed up in the session scratchpad.
