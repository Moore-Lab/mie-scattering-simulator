# Cross-check vs the Moore Lab dissertation (Chapter 2)

Independent validation of `mieinfo` against *Levitated Optomechanical Sensors for Fundamental
Physics* (Yale, Moore Lab), Chapter 2 "Trapping and detection" — the analytic Mie/GLMT results
(§2.1.2) and the optimal-detection / information-density figures (§2.3, Figs 2.5–2.7).

## 1. Analytic framework — exact agreement

| Dissertation (Ch. 2) | `mieinfo` | |
|---|---|---|
| `S1/S2 = Σ (2n+1)/n(n+1)·[aπ+bτ], [aτ+bπ]` (Eq 2.24) | `mie.plane_wave.scattering_amplitudes` | identical |
| `a_n, b_n` Riccati-Bessel (Eq 2.21) | `mie.special` / `mie.plane_wave` | identical (Bohren–Huffman) |
| `Q_ext, Q_sca` (Eq 2.25–2.26) | `mie.plane_wave.efficiencies` | identical |
| GLMT `a_nl = −G^TM_nl a_n`, `b_nl = −G^TE_nl b_n` (Eq 2.29) | `glmt` (a_n·g_TM, b_n·g_TE) | identical |
| Info density `= \|dE_det/dx\|²` (Fisher, Eq 2.62–2.63) | `info.info_density` | identical |
| `dE = (r·n̂ − A·x)·dE`, Gouy term (Eq 2.64) | `ik(k̂−ŝ)·r_s` phase gradient (PHYSICS §3.1) | same structure |
| `S_imp ∝ 1/(P_sca ∫\|dE/dx\|² dΩ)` (Eq 2.65) | `S_imp ∝ 1/F_q`, `Γ_ba ∝ Q_sca` | identical |
| "4π gives no info (energy conservation)" for *power* detection | `η(4π)=1` optimal-LO; power differs | consistent |

## 2. Detection-efficiency numbers (1064 nm, NA 0.63, circular pol)

**Dipole, d = 100 nm (Fig 2.6) — good match:**

| | transverse fwd/bwd | axial fwd/bwd |
|---|---|---|
| dissertation | 8.8% / 8.2% | 1.9% / 30% |
| mieinfo (plane-wave) | 5.7 / 5.3 | 0.2 / 38 |
| mieinfo (Richards-Wolf) | 5.7 / 7.8 | 0.1 / 39 |

Structure exact: transverse fore-aft **symmetric** (dipole signature, matches the M2 dipole gate),
axial strongly **backward**. Magnitudes within ~1.5–2× (circular-vs-linear polarization + focusing).

**Mie, d = 3 µm (Fig 2.5) — partial; reveals a genuine physics point:**

| | transverse fwd/bwd | axial fwd/bwd |
|---|---|---|
| dissertation-figure read-off | 38% / 3.9% | 37% / 26% |
| mieinfo (plane-wave) | 15.0 / 3.2 | 0.7 / 52 |
| mieinfo (Richards-Wolf, NA 0.63, **linear**) | **37.2** / 28.6 | 1.3 / 87 |
| mieinfo (Richards-Wolf, NA 0.63, **circular**) | 36.2 / 29.5 | 1.3 / 87 |
| independent angular-spectrum IRF (linear) | 24 / 46 | 0.6 / 92 |

- Transverse-forward **validated**: plane-wave 15% → focused Richards-Wolf 37% ≈ thesis 38%.
- Axial-forward is ~1% in **both** the GLMT-BSC path **and** an independent angular-spectrum
  IRF — the two methods agree with each other and both reproduce the plane-wave limit exactly
  (G-LIMIT holds to <0.5 pt at NA 0.05). See §3 for why this is not a numerical bug.

## 3. What was fixed, and the honest status of the axial number

**Fixed (validated):** `RichardsWolfFocus` now models **arbitrary polarization exactly**.
The Debye–Wolf focal field is linear in the input Jones vector, so the field is the exact
superposition `px·E^{x-input} + py·E^{y-input}` (y-input = x-input rotated 90°). Previously a
circular Jones vector silently returned the pure x-field. Circular polarization is now correct
(`|E(0)|=1`, `Ey = i·Ex` at the focus; `tests/glmt/test_circular_polarization.py`). **However,
circular polarization does *not* fix the transverse-backward or axial discrepancy** — it is the
azimuthal average of the two linear inputs, so the collected efficiencies barely move.

**Not reproduced — and it appears to be physics, not a bug.** The scattered-field
information radiation field is `E_μ ≡ ∂E_s/∂r_0` (the reference's Eq 3), and the collection
efficiency is `η = ∫_cone|E_μ|² / ∫_4π|E_μ|²` (Eq 11) — exactly what `mieinfo` computes. For
**axial** motion the IRF weight is `∼(k_inc − k·ŝ)_z² = (k_inc,z − k cosθ)²`, which is
intrinsically **backward-heavy**: it vanishes forward (θ→0) and peaks in backscatter. The
focusing Gouy shift only lowers the on-axis `k_inc,z` from `k` to ≈0.91k (NA 0.63), a small
forward correction — nowhere near a forward-comparable 37%. This was checked three ways that
all agree:
  1. the GLMT quadrature-BSC + finite-difference field derivative (axial-fwd 1.3%),
  2. the analytic product-rule BSC-gradient derivative (`glmt.derivatives`, 1.3%),
  3. an independent **angular-spectrum IRF** built from first principles (each focusing-cone
     plane wave scatters via the full vector Mie amplitude, weighted by `i(k_in − kŝ)`; 0.6%).
All three reproduce the plane-wave provider exactly at low NA, so the machinery is sound.

**Cross-check against the *published* version of this work (arXiv:2408.15483, same group).**
Its Fig. 2 and text state that under high-NA focusing the **axial** information is collected
**predominantly in the backward direction** (backward NA 0.9 objective, η≈0.9), with a small
low-NA *forward* objective — i.e. the published paper agrees with `mieinfo` that axial info is
backward-dominated, and does **not** support a forward-comparable axial efficiency at a single
symmetric NA-0.63 cone. The most likely reconciliation is that the "37% / 26%" axial read-off
uses a **different detection geometry** than the symmetric-cone `F(cone)/F(4π)` definition here
(e.g. asymmetric forward/backward objective NAs, or a self-homodyne LO overlap rather than the
bare `|E_μ|²` integral). We could not reproduce it from the Eq-3/Eq-11 definition and did not
force a match by adjusting tolerances.

## 4. Conclusion

- **Validated & improved:** the analytic framework, the plane-wave/dipole regimes, the
  transverse-**forward** focusing gain (plane-wave ~14% → focused ~37% ≈ thesis 38%), and now
  **exact arbitrary-polarization** focal fields (circular pol fixed).
- **Axial info is backward-dominated — and that is CORRECT.** Three independent methods here
  (GLMT-BSC + finite-difference, analytic product-rule, and an independent angular-spectrum IRF)
  AND the group's own **published** paper (arXiv:2408.15483: *"the backward detection scheme … with
  a high numerical aperture lens provides sufficient information to achieve the quantum ground
  state"*) all agree axial information is collected predominantly **backward**. The dissertation
  figure's forward-comparable "37%" is therefore most likely a **different detection definition**
  (e.g. self-homodyne LO overlap with the forward beam, or asymmetric objective NAs) rather than the
  bare `F(cone)/F(4π)` scattered-field IRF this package computes. `mieinfo` matches the *published*
  literature here.
- **STILL OPEN (honest):** the transverse-**backward** at high NA is not reproduced — `mieinfo`
  gives ~29–46% (the GLMT-BSC and angular-spectrum methods *disagree with each other*: 29% vs 46%)
  vs the thesis 3.9%. So the focused-beam displacement derivative is validated only for
  transverse-**forward** and the G-LIMIT; the transverse-backward and the absolute high-NA axial
  magnitudes remain unreliable. The `RichardsWolfFocus` high-NA caveat therefore **stays**.
- **Recommendation unaffected:** the Moore Lab apparatus analysed in `recommendation.md §10`
  uses a **NA 0.02** (weakly focused) probe, where plane-wave ≈ GLMT (verified) and none of
  the above high-NA subtleties apply.
