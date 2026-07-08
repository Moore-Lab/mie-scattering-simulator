# Validation

Correctness gates. Scientific code that isn't validated is a liability. Every
physics task passes at least one gate here (`CONVENTIONS.md §4`). `make validate`
runs all gates whose modules exist and exits nonzero on any failure; it is the
integration gate before every merge to `develop`.

Tolerances are contract values — loosening one is a `CONTRACT-CHANGE`
(`CONVENTIONS.md §5`).

## 1. Gate taxonomy

- **G-GOLD** regression vs seeded validated numbers.
- **G-LIMIT** analytic limits (plane-wave limit of GLMT, dipole limit, energy,
  reciprocity/symmetry).
- **G-DERIV** analytic vs finite-difference field derivatives.
- **G-CONV** convergence in `n_max` and angular grid.
- **G-LIT** simulation vs published benchmark.

## 2. G-GOLD — golden regression (owner W4a, source: prototype)

- `mieinfo.mie.plane_wave` reproduces `data/golden/golden_values.json`
  (`Q_ext, Q_sca, Q_back, g`) for all 9+2 cases to **relative ≤ 1e-9** (the
  prototype matches `miepython` to 1.4e-12; 1e-9 leaves margin for refactor noise).
- `S1(θ), S2(θ)` at `x ≈ 14.76` reproduce the prototype far-field to **≤ 1e-9**
  over `θ ∈ [0,π]` (engine check; wavelength-independent — it fixes `x`).
- The plane-wave information pattern (`PlaneWaveProvider` + phase-gradient
  derivative) reproduces `data/golden/information_pattern_results.json` — the
  **532 nm**, a=3–20 µm set: per-axis forward-info fractions, peak polar angles, and
  `η` at NA 0.5/0.8/0.95 for x and z, **forward and backward** — to **relative
  ≤ 1e-6** (grid-dependent; use the matching grid, 700–1400 θ nodes by size).

Anchor numbers a session can eyeball (`PHYSICS.md §8`): transverse forward
`η_x(NA0.8,fwd)` ≈ 0.52 (a=3 µm) → 0.59 (a=20 µm); axial backward
`η_z(NA0.8,bwd)` ≈ 0.65–0.69 vs forward ≈ 0.055–0.069; z-info peak ≈ 155–159°.
Backward collection is decisive for the beam-collinear axis.

## 3. G-LIMIT — analytic limits (owner W4a)

- **Plane-wave limit of GLMT:** with `GaussianParaxial` waist → ∞ (or a plane-wave
  BSC path), `glmt.scatter` reproduces `mie.plane_wave` `S1/S2` to **≤ 1e-6**, and
  all three BSC methods return the plane-wave weights to **≤ 1e-6**. Gates M1.
- **Dipole (Rayleigh) limit:** for `x = 0.05` the scattered field → dipole field;
  the information pattern reproduces the Tebbenjohanns-2019 dipole result —
  transverse pattern ∝ dipole × `s_⟂²`, axial pattern backward-weighted, and the
  ideal-scheme `η → 1` at full `4π` with optimal LO. Assert the pattern shape
  (normalized) and the axial-vs-transverse forward-fraction ordering. Gates M2.
- **Energy / optical theorem:** `Q_ext` from the forward amplitude equals the
  coefficient-sum `Q_ext`; for real `m`, `Q_ext = Q_sca` to **≤ 1e-10**; for
  absorbing `m`, `Q_ext > Q_sca`.
- **Symmetry/reciprocity:** for `x`-polarized incidence the intensity and information
  patterns satisfy the expected `φ`-parity (even in `φ`, period `π`). **Correction
  (W4a, 2026-07-07):** a `y`-displacement pattern does **not** equal the `x`-pattern
  under a bare `φ→φ−π/2` relabel for a Mie sphere — the incident `x`-polarization breaks
  the transverse `x/y` symmetry, so `F_x ≠ F_y` (with `F_y → 2 F_x` in the dipole limit,
  since `∫sin²θ|S1|² ≠ ∫sin²θ|S2|²`). The true `x↔y` symmetry rotates the incident
  polarization together with the displacement. Property-tested (`φ`-parity holds;
  `F_x ≠ F_y` is asserted, not the false equality).

## 4. G-DERIV — displacement derivative (owner W1c + W4a)

- `field_derivative(method='analytic')` matches `method='finite_difference'`
  (central, step auto-scaled to `~1e-4·λ` and Richardson-checked) to **relative
  ≤ 1e-6** in `E_θ`, `E_φ`, for each Cartesian direction, over a grid of `r_s`
  spanning the trap region (e.g. `|r_s| ≤ 0.3 λ`). Gates M1.
- Consistency: for `PlaneWaveProvider`, the analytic phase-gradient derivative
  (`PHYSICS.md §3.1`) equals the finite-difference derivative of the analytic-phase
  field to machine precision.

## 5. G-CONV — convergence (owner W4a)

- Efficiencies, patterns, and `η(Ω)` stable to **relative ≤ 1e-4** under
  `n_max → n_max + 8` beyond the Wiscombe value, across `x ∈ {3, 15, 30, 60}`.
- Solid-angle integrals stable to **≤ 1e-4** under doubling `Ntheta, Nphi` (Gauss–
  Legendre in `cosθ`). Report the grid at which each quantity is converged; the
  optimizer uses the coarsest converged grid.
- **Caveat (W4a, 2026-07-07):** the `≤ 1e-4` grid-convergence holds for **smooth**
  solid-angle integrals (angle-integrated `g`, `F_q(4π)`). A **hard-edged NA collection
  cone** does *not* converge to `1e-4` under grid doubling — the sharp `cone_mask`
  boundary quantizes onto the Gauss–Legendre nodes, giving ~`1e-3`–`1e-2` wobble in
  `η(Ω)` at `x=60`. Mitigate by using a sufficiently fine grid for the reported `η`, or
  a soft-edged cone; the optimizer must fix its grid when ranking cone geometries so the
  quantization is common-mode across candidates.

## 6. G-LIT — literature benchmarks (owner W4c, data from W3)

- For each `Benchmark` in the DB, `compare_benchmark` predicts the target quantity
  within the benchmark's stated `target_tolerance`, **or** the discrepancy is
  explained (regime outside model validity, missing parameter, definition
  mismatch) in `discrepancy_note`. At least the Tebbenjohanns dipole case and one
  Mie-regime or experimental efficiency benchmark must pass numerically. Gates M4.
- Comparison uses the same efficiency *definition* as the source; W4c documents any
  normalization conversion. A "pass" by redefining the quantity is not a pass.

## 7. CI lanes

- **PR lane (fast):** unit tests, `G-GOLD`, `G-DERIV` on small `x`, symmetry
  property tests, lint/type. Must be green to merge.
- **Nightly/full lane:** `@slow` BSC-quadrature cross-checks, full `G-LIMIT`,
  `G-CONV` sweep, all `G-LIT` benchmarks, full-`x` golden regression.
- The orchestrator runs the full lane at each milestone gate before tagging.

## 8. What "validated" means for the final recommendation (M5)

The recommendation is only reportable if: (a) the forward model passed G-LIMIT +
G-DERIV + G-CONV, (b) the information engine passed the dipole-limit gate and the
golden `η(NA)` regression, (c) at least two independent literature/analytic results
are reproduced (G-LIT), and (d) the sensitivity analysis shows which conclusions
depend on unconfirmed assumptions (`A#`). Anything not meeting this is reported as
provisional with the gap named.
