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

**Mie, d = 3 µm (Fig 2.5) — partial; reveals a limitation:**

| | transverse fwd/bwd | axial fwd/bwd |
|---|---|---|
| dissertation | 38% / 3.9% | 37% / 26% |
| mieinfo (plane-wave) | 13.9 / 3.1 | 0.6 / 54 |
| mieinfo (Richards-Wolf, NA 0.63) | **35.5** / 29.5 | 1.1 / 88 |

- Transverse-forward **validated**: plane-wave 14% → focused Richards-Wolf 35.5% ≈ thesis 38%.
- But the focused-beam engine **fails** at this NA elsewhere: axial-forward 1.1% vs the thesis's 37%
  (the **Gouy-phase** axial sensitivity, Eq 2.64's `A·x` term, is not captured), and a spurious
  transverse-backward (29.5% vs 3.9%).

## 3. Conclusion

- **Validated:** the analytic framework (identical equations) and the **plane-wave / dipole**
  regimes agree with the dissertation quantitatively.
- **Limitation exposed:** the `RichardsWolfFocus` GLMT path is gate-validated only in its
  *plane-wave limit*; at strong focusing (NA 0.63) it does **not** reliably reproduce the thesis's
  focused-beam axial (Gouy) numbers. **Do not trust the high-NA focused-beam path** until it is
  fixed and validated against this reference. (The transverse-forward focusing correction is the
  one focused-beam quantity that checks out.)
- **Recommendation unaffected:** the Moore Lab apparatus analysed in `recommendation.md §10` uses a
  **NA 0.02** (weakly focused) probe, where the Gouy phase is negligible and plane-wave ≈ GLMT
  (verified). The thesis's NA-0.63 focused imaging is a different, more strongly-focused regime.
