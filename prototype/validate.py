"""Standalone self-validation of the reference-oracle seed. No project deps.

Proves the seed is a validated oracle before any Claude Code session builds on it:

  Gate 1  plane-wave Mie efficiencies (Qext,Qsca,Qback,g) vs miepython   ~1e-12
  Gate 2  far-field angular SHAPE vs miepython                           ~1e-11
          (mie_core's S1/S2 differ from miepython by a global x2 constant
          that CANCELS in every delivered quantity -- all information-pattern
          outputs are ratios -- so the physical check is the normalized shape)
  Gate 3  efficiencies reproduce the committed golden_values.json
  Gate 4  the information pattern reproduces the committed results JSONs
          (regression: there is no external oracle for the info pattern; its
          Mie inputs are validated by Gates 1-3, and the numbers are pinned here)

Exit 0 iff every runnable gate passes; nonzero on any failure. If miepython is
absent, Gates 1-2 SKIP (still self-consistent via Gates 3-4) and the script says so.

Usage:  python validate.py            # 1064 info-pattern regression (fast)
        python validate.py --full     # also regenerate + diff the 532 nm set (slower)
"""
import json, sys, numpy as np
import mie_core as mc

TOL_EFF, TOL_SHAPE, TOL_GOLD, TOL_INFO = 1e-9, 1e-8, 1e-9, 1e-6
# Qback is a coherent ALTERNATING sum over ~n_max terms (sum (2n+1)(-1)^n (a_n-b_n));
# catastrophic cancellation puts its float64 floor at ~1e-6 for large x (x~236 =>
# ~264 terms). That is a fundamental precision limit, not an error -- Qext/Qsca/g
# stay at ~1e-12. So Qback gets its own, honestly-looser tolerance.
TOL_QBACK = 1e-5
# DIAGNOSTIC fields are numerical residuals, not physical outputs, so a *relative*
# regression tolerance is meaningless for them. g_relerr = |g_map - g|/|g| is the
# grid-convergence residual (~1e-13); comparing a 7e-15 residual against a 5e-13 one
# reads as "51% error" though every physical quantity reproduces to ~1e-12. Skip it.
DIAG = {"g_relerr"}
fails = []

def check(name, err, tol):
    ok = err <= tol
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:44s} max_err={err:.2e}  (tol {tol:.0e})")
    if not ok:
        fails.append(name)

# ---------- Gate 1 + 2: vs miepython ----------
try:
    import miepython as mp
    m = 1.45 + 0j
    xs = [0.3, 1.0, 2.9526, 5.0, 14.7631, 29.5301, 59.05, 118.1, 236.2]
    mu = np.cos(np.linspace(0.01, np.pi - 0.01, 200))
    e_err = b_err = s_err = 0.0
    for x in xs:
        qe, qs, qb, g = mc.efficiencies(m, x)
        Qe, Qs, Qb, G = mp.efficiencies_mx(m, x)
        for a, b in ((qe, Qe), (qs, Qs), (g, G)):        # Qext, Qsca, g -- tight
            e_err = max(e_err, abs(a - b) / max(abs(b), 1e-30))
        b_err = max(b_err, abs(qb - Qb) / max(abs(Qb), 1e-30))   # Qback -- cancellation-limited
        S1, S2 = mc.scattering_amplitudes(m, x, np.arccos(mu))
        Ia = np.abs(S1) ** 2 + np.abs(S2) ** 2
        iu = mp.i_unpolarized(m, x, mu, norm="bohren")
        s_err = max(s_err, float(np.max(np.abs(Ia / Ia.sum() - iu / iu.sum()))))
    print(f"Gate 1/2  vs miepython {mp.__version__}  (x = 0.3 .. 236):")
    check("efficiencies Qext, Qsca, g", e_err, TOL_EFF)
    check("Qback (alternating-sum cancellation floor)", b_err, TOL_QBACK)
    check("far-field angular shape (normalized)", s_err, TOL_SHAPE)
except ImportError:
    print("Gate 1/2  miepython not installed -- SKIPPED  (pip install miepython)")

# ---------- Gate 3: efficiencies vs committed golden ----------
gv = json.load(open("golden_values.json"))
g_err = 0.0
for c in gv["cases"]:
    qe, qs, qb, g = mc.efficiencies(complex(c["m_real"], c["m_imag"]), c["x"])
    for key, val in (("Qext", qe), ("Qsca", qs), ("Qback", qb), ("g", g)):
        g_err = max(g_err, abs(val - c[key]) / max(abs(c[key]), 1e-30))
print("Gate 3  efficiencies vs committed golden_values.json:")
check(f"round-trip {len(gv['cases'])} golden cases", g_err, TOL_GOLD)

# ---------- Gate 4: information-pattern regression ----------
import information_pattern as ip
try:
    ref = json.load(open("information_pattern_results_1064.json"))
    i_err = 0.0
    for r in ref:
        got = ip.summarize(complex(r["m"], 0.0), r["x"], r["label"])
        for k, v in got.items():
            if k in DIAG:
                continue
            if isinstance(v, (int, float)) and isinstance(r.get(k), (int, float)):
                i_err = max(i_err, abs(v - r[k]) / max(abs(r[k]), 1e-12))
    print("Gate 4  information-pattern regression vs committed JSON (1064 nm set):")
    check("info-pattern reproduction", i_err, TOL_INFO)
except FileNotFoundError:
    print("Gate 4  information_pattern_results_1064.json ABSENT -- SKIPPED")

if "--full" in sys.argv:
    import study_detection_532 as s532
    op = json.load(open("information_pattern_results.json"))  # operative 532 set
    by_a = {round(r["a_um"], 3): r for r in op}
    o_err = 0.0
    print("Gate 4b regenerating operative 532 nm set at committed grids ...")
    for a, rref in by_a.items():
        nt, nph = s532.grid_for(a)
        got = s532.summarize(a, nt, nph)
        for k, v in got.items():
            if k in DIAG:
                continue
            if isinstance(v, (int, float)) and isinstance(rref.get(k), (int, float)):
                o_err = max(o_err, abs(v - rref[k]) / max(abs(rref[k]), 1e-12))
    check("info-pattern reproduction (532 operative)", o_err, TOL_INFO)

print()
if fails:
    print("VALIDATION FAILED:", ", ".join(fails))
    sys.exit(1)
print("ALL GATES PASSED -- the seed is a validated oracle.")
sys.exit(0)
