# Physics Ground Truth

This is the authoritative theory reference. Implement against it; do not
re-derive from memory. Every equation here has an associated validation target in
`VALIDATION.md`. Notation is fixed in §0 and used identically across all modules.

Symbol conventions match Bohren & Huffman (B&H) for Mie and Gouesbet & Gréhan for
GLMT; where they differ, B&H wins and the difference is flagged.

---

## 0. Notation and conventions (frozen)

- `k = 2π/λ` wavenumber in the medium; medium is vacuum unless stated (`λ` = vacuum λ).
- **Two wavelengths.** `λ_trap = 1064 nm` holds the sphere; `λ_det = 532 nm` is the
  probe/imaging light. **All scattering, information, and detection physics uses
  `λ_det`** — the trap wavelength enters only if the trap beam is itself a probe
  (it is not, on this apparatus). Silica index is dispersive: `n(1064) ≈ 1.4496`,
  `n(532) ≈ 1.4607`. Use `λ_det` (and its `n`) throughout unless a quantity is
  explicitly about the trap.
- `a` sphere radius; `x = k a = 2πa n_med/λ_det` size parameter (detection light).
- **Per-beam frame.** All single-beam formulae below are written in that beam's own
  frame: `k_hat = +z` is *that beam's* propagation axis, `θ` measured from it.
  "Forward" = small `θ` (beam-collinear), "backward" = `θ→π`. A lab axis that is
  transverse to the beam behaves like `x`/`y` here; the axis collinear with the beam
  behaves like `z`. Multi-beam configurations rotate each beam's frame into the lab
  frame (`§4.5`).
- `m = n_particle / n_medium`, complex allowed; convention: `Im(m) ≥ 0` for
  absorption with time dependence `e^{-iωt}` (B&H). **State the `e^{-iωt}` vs
  `e^{+iωt}` convention in every module docstring; the code base uses `e^{-iωt}`.**
- Incident direction `k_hat = +z`. Default incident polarization `x_hat`.
- Observation direction `s_hat = (sinθ cosφ, sinθ sinφ, cosθ)`.
- Riccati–Bessel: `ψ_n(z) = z j_n(z)`, `χ_n(z) = −z y_n(z)`,
  `ξ_n(z) = ψ_n(z) − i χ_n(z) = z h_n^{(1)}(z)`.
- `D_n(z) = ψ_n'(z)/ψ_n(z)` logarithmic derivative (downward recurrence).
- Sphere displacement from beam focus / origin: `r_s = (x_s, y_s, z_s)`.
- Fisher information about parameter `q`: `F_q`. Detection efficiency `η ∈ [0,1]`.

---

## 1. Plane-wave Mie (baseline, VALIDATED in `prototype/`)

### 1.1 Coefficients (numerically stable B&H form)

Compute `D_n(mx)` by downward recurrence. **The burn-in margin above `|mx|` must
scale with `|mx|`** — use `n_start = max(n_max, ⌈|mx|⌉) + 15 + ⌈15√|mx|⌉` (textbook
BHMIE). A fixed margin (e.g. `+16`) silently under-resolves `D_n` at large `|mx|`:
at the 532 nm detection wavelength `x=236, m=1.45 ⇒ |mx|=342`, a `+16` margin gives
`D_1` wrong by ~7e-3 and *every* `a_n` wrong by ~5e-4 (uniform across order), which
propagates to `Q_sca`/`Q_back` at 1e-5–1e-4 and to `Γ_ba ∝ Q_sca` directly. This
was a real bug in the seed, caught only when the detection wavelength pushed `x`
past ~60; `prototype/validate.py` regression-guards it against miepython.

```
D_{n-1} = n/(mx) − 1/(D_n + n/(mx)),    D_{n_start} = 0
```

`ψ_n(x)`, `χ_n(x)` by upward recurrence with
`ψ_{-1}=cos x, ψ_0=sin x, χ_{-1}=−sin x, χ_0=cos x`:

```
ψ_n = (2n−1)/x · ψ_{n-1} − ψ_{n-2}    (same recurrence for χ)
ξ_n = ψ_n − i χ_n
```

`χ_n` grows with order, so in float64 the upward recurrence **overflows to NaN**
for very high `n` (≈ n=887 at x=236). Wiscombe `n_max` (≈267 there) is far below
that, so it is safe — but do not naively crank `n_max` past the overflow point;
if higher orders are ever needed, scale `χ_n` or use a stable library evaluation.

Then

```
a_n = [ (D_n/m + n/x) ψ_n − ψ_{n-1} ] / [ (D_n/m + n/x) ξ_n − ξ_{n-1} ]
b_n = [ (m D_n + n/x) ψ_n − ψ_{n-1} ] / [ (m D_n + n/x) ξ_n − ξ_{n-1} ]
```

### 1.2 Cross sections / efficiencies

```
Q_ext = (2/x²) Σ_{n≥1} (2n+1) Re(a_n + b_n)
Q_sca = (2/x²) Σ_{n≥1} (2n+1) (|a_n|² + |b_n|²)
Q_back = (1/x²) |Σ_{n≥1} (2n+1)(−1)^n (a_n − b_n)|²
g Q_sca = (4/x²)[ Σ_{n} n(n+2)/(n+1) Re(a_n a*_{n+1}+b_n b*_{n+1})
                  + Σ_n (2n+1)/(n(n+1)) Re(a_n b*_n) ]
```

`Q_back` is a coherent **alternating** sum; at large `x` catastrophic cancellation
limits it to ~1e-5–1e-6 relative in float64 (vs ~1e-12 for `Q_ext, Q_sca, g`).
`validate.py` gives it a correspondingly looser tolerance — that is a precision
floor, not an error.

For non-absorbing spheres `Q_ext = Q_sca` (checked to 1e-12 in the prototype).

### 1.3 Far-field amplitudes

```
S1(θ) = Σ_n (2n+1)/(n(n+1)) [ a_n π_n(cosθ) + b_n τ_n(cosθ) ]
S2(θ) = Σ_n (2n+1)/(n(n+1)) [ a_n τ_n(cosθ) + b_n π_n(cosθ) ]
```

with angle functions `π_1=1, π_0=0`,
`π_n = (2n−1)/(n−1) μ π_{n-1} − n/(n−1) π_{n-2}`, `τ_n = n μ π_n − (n+1) π_{n-1}`,
`μ = cosθ`. For `x_hat`-polarized incidence the far scattered field components are
`E_θ ∝ S2(θ) cosφ`, `E_φ ∝ −S1(θ) sinφ`, so the scattered intensity map is

```
I(θ,φ) ∝ |S2(θ)|² cos²φ + |S1(θ)|² sin²φ.
```

This exact code is in `prototype/mie_core.py` and matches `miepython` to 1e-12.
`mieinfo.mie.plane_wave` must reproduce it (regression target).

---

## 2. Generalized Lorenz–Mie theory (GLMT) — the physical trap model

The trap is a strongly focused beam. The scattered field is still expanded in
vector spherical wave functions (VSWFs), but the incident-field expansion carries
**beam-shape coefficients (BSCs)** `g_{n,TM}^m`, `g_{n,TE}^m` that replace the
plane-wave weights. For a plane wave the BSCs collapse to the `(2n+1)` weighting
of §1 (this is a validation limit).

### 2.1 Scattered-field expansion

With BSCs of the incident beam expanded about the sphere center, the scattered
partial-wave amplitudes are

```
scattered TM amplitude ∝ a_n g_{n,TM}^m,     scattered TE amplitude ∝ b_n g_{n,TE}^m,
```

i.e. the *same* Mie coefficients `a_n, b_n` from §1.1 multiply the beam-specific
BSCs. The far field is the VSWF sum; provide it as
`E_s(θ,φ; r_s, beam)` (complex vector) and the differential scattered power.

### 2.2 Beam-shape coefficients — implement ≥ 2 methods, cross-check

1. **Quadrature / projection (reference, slow):** project the focused-beam field
   onto VSWFs by numerical integration over a spherical surface. Ground truth for
   the others.
2. **Localized approximation (Gouesbet–Gréhan, fast):** closed-form BSCs for a
   Gaussian beam via the localized operator. Fast enough for optimization loops.
   Valid for moderate focusing; **flag validity limits** for high NA.
3. **Angular-spectrum / vector-diffraction (Richards–Wolf) for tight focus:**
   the correct high-NA focal field (an aplanatic lens focusing a plane wave) is the
   Richards–Wolf integral; project that onto VSWFs. Required because the silica
   trap is high-NA and the paraxial Gaussian localized approximation is not
   trustworthy there.

Represent the beam abstractly (see `INTERFACES.md IncidentBeam`) so plane wave,
paraxial Gaussian, and Richards–Wolf focus are interchangeable behind one API.

### 2.3 References for GLMT/BSC

Gouesbet & Gréhan, *Generalized Lorenz–Mie Theories*, 2nd ed. (2017); Gouesbet,
"T-matrix formulation and the g_n coefficients"; Doicu & Wriedt on the localized
approximation; Richards & Wolf, Proc. R. Soc. A 253, 358 (1959) for the focal
field. W3 confirms exact editions/volumes.

---

## 3. Displacement of the sphere — the information source

Two equivalent routes; implement the analytic one and check it against the
numerical one.

### 3.1 Plane-wave shortcut (exact, cheap — used in the prototype)

For a plane wave, translating the sphere by `r_s` multiplies the far-field
scattered wave by a direction-dependent phase:

```
E_s(s_hat; r_s) = E_s(s_hat; 0) · exp[ i k (k_hat − s_hat) · r_s ]
```

so

```
∂E_s/∂r_{s,j} = i k (k_hat − s_hat)_j · E_s(s_hat; 0).
```

This is the exact plane-wave displacement derivative and the seed for the
information-pattern demo. **It is only correct for a plane wave** — for a focused
beam the BSCs themselves depend on `r_s` and the derivative has extra terms.

### 3.2 GLMT translation-addition theorem (general)

Displacing the sphere relative to the beam focus re-expands the incident beam
about the new center: the BSCs transform via the **VSWF translation-addition
theorem** (Cruzan/Stein coefficients). Then

```
∂E_s/∂r_{s,j} = (∂/∂r_{s,j}) [ VSWF sum with translated BSCs ].
```

Two sub-strategies (implement both; the second gates M1):
- **Analytic:** differentiate the translated BSCs w.r.t. `r_s` using
  recurrence relations for the translation coefficients.
- **Numerical (validation):** central finite difference of `E_s(θ,φ; r_s)` in
  each Cartesian direction. The analytic derivative must match this to `< 1e-6`
  relative (`VALIDATION.md §3`), over a grid of `r_s` inside the trap.

Displacement scale of interest: motional amplitudes are ≪ λ (zero-point and
thermal), so the linear (first-derivative) response is what matters; second order
only if a nonlinearity study is requested.

---

## 4. Fisher information for coherent optical position detection

This is the heart of the project. Framework follows Tebbenjohanns–Frimmer–Novotny
(PRA 100, 043821, 2019) for the dipole and Maurer–González-Ballestero–Romero-Isart
(PRA 108, 033714, 2023) for the full Lorenz–Mie particle.

### 4.1 Information radiation pattern

The scattered field is a coherent (quasi-classical) field. Estimating a
displacement `q = r_{s,j}` from an ideal shot-noise-limited coherent measurement
that collects the field over solid angle `Ω` and mode-matches a local oscillator
to `∂E_s/∂q`, the classical Fisher information rate is

```
F_q(Ω) ∝ ∫_Ω |∂E_s/∂q (s_hat)|² dΩ          (per unit photon flux normalization)
```

Define the **information radiation pattern** as the integrand density

```
dF_q/dΩ ∝ |∂E_s/∂q (s_hat)|².
```

For the plane-wave case (§3.1) this is `k² (k_hat − s_hat)_j² |E_s(s_hat)|²`,
which is **not** the intensity pattern `|E_s|²`: it is reweighted by
`(k_hat − s_hat)_j²`. Consequences (reproduced numerically in `prototype/`):

- **Transverse (x, y):** weight `(s_hat_⟂)²` = `sin²θ cos²φ` (x). Vanishes forward,
  peaks off-axis (~40–62° in the prototype). Info is more side/forward-weighted.
- **Axial (z):** weight `(1 − cosθ)²`. Zero forward, maximal (=4) in backscatter.
  Axial-motion information is **backward-weighted**; forward collection captures
  only a few percent. This matches Tebbenjohanns 2019 ("backscattering detection
  provides sufficient information to cool axial motion below unity phonon").

### 4.2 Detection efficiency and imprecision

```
η_q(Ω) = F_q(Ω) / F_q(4π) ∈ [0,1]
S_q^{imp}(ω) ∝ 1 / F_q(Ω)            (imprecision noise PSD ∝ 1/collected info)
```

Real optics cap `Ω` (finite NA), and the LO mode-match may be imperfect
(`η → κ·η`, `κ ≤ 1`). The optimizer maximizes `η_q` (or a mode-weighted
combination) over the achievable configuration set.

### 4.3 Backaction and the standard quantum limit

Measurement backaction = photon-recoil heating, set by the **total** scattered
power (all `4π`), independent of what is collected:

```
Γ_ba ∝ P_scattered ∝ Q_sca.
```

The imprecision–backaction product is bounded (Heisenberg limit); the optimal
scheme saturates it. Since `Γ_ba` is fixed by `Q_sca` and imprecision `∝ 1/F_q(Ω)`,
the product is minimized by maximizing `η_q`. Hence **`η_q` is the single figure
of merit**, and reaching the SQL requires `η_q → 1` (with the optimal LO). The
package reports `η_q`, `S_q^{imp}`, `Γ_ba`, and the distance to the SQL for each
configuration.

### 4.4 The sensed parameter (general; arbitrary direction)

The sensed parameter `q` is general. The implemented family is **displacement along
an arbitrary unit vector `n_hat`**:
`∂E_s/∂(r_s·n_hat) = Σ_j n_hat_j ∂E_s/∂r_{s,j}` — a linear combination of the three
Cartesian derivatives, so an arbitrary direction costs no new field solves. Keep the
derivative behind a `parameter` abstraction (`INTERFACES.md §4`) so other parameters
(e.g. orientation of an anisotropic particle) can be added later; for a homogeneous
sphere only translation carries signal, so v1 implements translation.

*Coriolis instance.* The Coriolis observable is a velocity inferred from a position
time series; velocity sensitivity is set by position imprecision on the driven mode.
The drive is a horizontal COM loop (N→E→S→W, then reversed); the sensed quantity is
transverse `x`/`y` position (`MASTER_PLAN.md F4`). This is one call into the general
`n_hat` machinery, not a special code path.

### 4.5 Multiple detection channels (this apparatus has two beams)

Each probe beam + collection port is an independent **channel** `c`: its own
propagation axis `k_hat_c`, wavelength, polarization, focusing, and collection
geometry `Ω_c`. Compute the information pattern in each beam's own frame (`§0`),
then rotate `s_hat` into the lab frame to apply that beam's collection cone.

For **independent** channels (independent shot noise — distinct beams/detectors),
Fisher information **adds**:
```
F_q^total(config) = Σ_c  η_{q,c}(Ω_c) · F_{q,c}(4π_c)
```
i.e. total collected information about `q` is the sum of each channel's collected
information. This is why the two-beam layout works: a lab axis that is *collinear*
with beam A (poorly measured in A's forward lobe — see `§4.1`, the `(1−cosθ)²` axial
weight) is *transverse* to beam B (well measured in B's forward lobe, `sin²θ`
weight). The optimizer (`W4b`) searches over the **set** of channels, not a single
detector. Caveat: if two channels share a detector or LO their noise is correlated
and information does not simply add — flag that case; default channels are
independent.

---

## 5. Detection schemes and collection optics (model each)

- **Split / quadrant detection:** balanced difference of a split photodiode in the
  collected far field. Effective mode is the antisymmetric spatial mode; compute
  its overlap with `∂E_s/∂q` to get the achieved fraction of `F_q(Ω)`.
- **Balanced homodyne / heterodyne with external LO:** the LO spatial mode is a
  free function; the optimal LO ∝ `∂E_s/∂q` saturates `F_q(Ω)`. Report both the
  optimal-LO efficiency and the efficiency for a realistic Gaussian LO.
- **Self-homodyne / imaging (probe beam as its own LO):** the unscattered probe
  (532 nm imaging beam — **not** the trap) is the reference; scattered + probe
  interfere on a split/quadrant detector or camera. This is the likely default
  hardware. Model per beam and per direction: forward self-homodyne emphasizes the
  forward interference term (good for the two axes transverse to that beam, poor for
  the collinear axis); backward collection recovers the collinear axis. Quantify the
  shortfall vs the optimal LO, per direction, per beam.
- **Collection optics:** finite-NA cone (forward and/or backward), aplanatic lens
  apodization (`√cosθ` factor), aperture truncation, and mode-matching loss into a
  single-mode fiber if applicable. NA relates to cone half-angle by
  `NA = n_medium sin α` (vacuum here, so `NA = sin α`, capped below 1).

Each scheme is a map from `(E_s, ∂E_s/∂q, geometry)` to `(η_q, S_q^{imp})`.

---

## 6. Limits and internal consistency (used as validation gates)

- **Plane-wave limit of GLMT:** wide beam waist ⇒ BSCs → plane-wave weights ⇒
  §1 results. (`VALIDATION.md §3`.)
- **Rayleigh/dipole limit:** `x ≪ 1` ⇒ scattered field → dipole field; the
  information pattern → Tebbenjohanns 2019 dipole result. (M2 gate.)
- **Energy conservation / optical theorem:** `Q_ext` from forward amplitude equals
  the coefficient sum; absorbing-sphere `Q_ext > Q_sca`.
- **Reciprocity / symmetry:** patterns respect the polarization symmetry
  (`φ → φ+π` etc.).
- **Convergence:** results stable under `n_max → n_max + Δ` and under angular-grid
  refinement (Gauss–Legendre in `cosθ`).

---

## 7. Numerical regime and cost (drives architecture)

Validated in the study phase. **Detection is at `λ_det = 532 nm`**, so the operative
size parameters are large:

| a (µm) | x = 2πa/λ_det (532 nm) | n_max (Wiscombe) |
|---|---|---|
| 3   | 35.4  | 51  |
| 5   | 59.0  | 77  |
| 8   | 94.5  | 115 |
| 12  | 141.7 | 165 |
| 20  | 236.2 | 264 |

(For comparison, the same radii at the 1064 nm trap wavelength give half the `x`;
the trap value is irrelevant to detection but useful for GLMT trapping-field code.)

- Plane-wave Mie far-field: ~1 ms per evaluation at `x≈30`; scales with `n_max`.
- One `η(NA)` curve over a 600×360 grid: ~60 ms at small `x`. Grid-convergence check
  (angle-integrated `g` vs coefficient `g`) passes to `≤ 1e-10` up to `x=236` on a
  ~1400×420 grid — **large `x` needs a finer angular grid** (forward lobe width
  `~1/x ~ 4 mrad` at `x=236`); the study used 700–1400 `θ` nodes over `[3,20] µm`.
- **Cost driver at large `x` is GLMT translation-addition.** BSCs and the
  translation step per displacement dominate, and translation-addition at
  `n_max ~ 264` is heavy (translation coefficients scale steeply in order). Mitigate
  by: caching BSCs (they depend only on beam + `r_s`, reusable across DOF, direction,
  and the whole angular grid), vectorizing VSWF evaluation over the grid,
  precomputing angular functions once per `n_max`/grid, and — for the size sweep —
  sampling `a` coarsely at the top end (efficiencies are flat across 3–20 µm, `§8`).
- Optimization inner loop must avoid recomputing BSCs when only the collection
  geometry (a mask on a fixed field) or the channel's collection cone changes.

---

## 8. Golden values (regression targets)

`prototype/golden_values.json` holds validated plane-wave `Q_ext, Q_sca, Q_back, g`
(engine regression, wavelength-independent — they fix `x`), silica at nine `(m,x)`
points `x = 0.30 … 29.5` plus two sanity cases, all matching `miepython` to
`≤ 1.4e-12`.

`prototype/information_pattern_results.json` holds the operative info-pattern golden
set **at `λ_det = 532 nm`** (`n=1.4607`, vacuum) for `a = 3, 5, 8, 12, 20 µm`:
per-axis forward info fraction, peak polar angle, and `η` at NA 0.5/0.8/0.95 for
**both forward and backward** collection. (`information_pattern_results_1064.json` is
retained for reference / small-`x` dipole behavior.) Representative numbers to
sanity-check against:

- **Transverse (x), forward collection:** info forward fraction 0.90–0.94, peak
  ≈ 40–55°; `η_x(NA0.8, fwd)` ≈ 0.52 (a=3 µm) → 0.59 (a=20 µm). Backward is useless
  for transverse: `η_x(NA0.8, bwd)` ≈ 0.02.
- **Axial (z), backward collection is the story:** z-info forward fraction ≈ 0.21,
  peak ≈ 155–159° (backscatter). Forward collection is poor — `η_z(NA0.8, fwd)`
  ≈ 0.055–0.069 — but **backward** collection gives `η_z(NA0.8, bwd)` ≈ 0.65–0.69.
  With the backward port available (`MASTER_PLAN.md F2`), the collinear axis goes
  from ~6% to ~65–70% collected.
- **Size is a weak knob across 3–20 µm:** `η` varies smoothly and slightly (larger
  spheres marginally better for transverse). The strong Mie-resonance size
  dependence seen at sub-micron radii (e.g. the old a=2.5 µm forward dip at 1064 nm)
  has washed out at `x ≳ 35`.

These are plane-wave numbers in each beam's own frame. The GLMT engine shifts them;
its plane-wave limit must reproduce them, and the dipole limit (`x≪1`) must
reproduce Tebbenjohanns 2019.

---

## 9. Key references (seed; W3 expands and verifies in `LITERATURE.md`)

1. Tebbenjohanns, Frimmer, Novotny, "Optimal position detection of a dipolar
   scatterer in a focused field," **PRA 100, 043821 (2019)** (arXiv:1907.12838).
   *Dipole information pattern; backscatter-optimal axial detection; Heisenberg
   imprecision–backaction limit.* Primary theory anchor.
2. Maurer, González-Ballestero, Romero-Isart, "Quantum theory of light interaction
   with a Lorenz–Mie particle: optical detection and 3D ground-state cooling,"
   **PRA 108, 033714 (2023)**; Erratum **PRA 109, 019901 (2024)**. *Full-Mie
   generalization — the direct theoretical target of this project.*
3. Gonzalez-Ballestero, Aspelmeyer, Novotny, Quidant, Romero-Isart,
   "Levitodynamics," **Science 374, eabg3027 (2021)**. Field review / context.
4. Bohren & Huffman, *Absorption and Scattering of Light by Small Particles*
   (1983). Plane-wave Mie reference (§1 algorithm).
5. Gouesbet & Gréhan, *Generalized Lorenz–Mie Theories*, 2nd ed. (2017). GLMT/BSC.
6. Richards & Wolf, **Proc. R. Soc. A 253, 358 (1959)**. High-NA focal field.
7. Recent, on-topic (verify): arXiv:2409.00782 "Optimal displacement detection of
   arbitrarily-shaped levitated dielectric objects"; arXiv:2512.17894 "Visualizing
   detection efficiency in optomechanical scattering"; Pluchar et al.,
   "Imaging-based quantum optomechanics," **PRL 135, 023601 (2025)**.
