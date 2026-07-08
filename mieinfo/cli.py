"""mieinfo CLI (INTERFACES.md §10): validate | optimize | compare | report.

    mieinfo validate [--full]   run the validation gates; nonzero exit on failure
    mieinfo optimize [config]   OptResult for a RunConfig (or the silica default) -> JSON + figures
    mieinfo compare             run all literature benchmarks -> comparison table + JSON
    mieinfo report              assemble the recommendation inputs (optimize + compare)

Time convention e^{-iωt}. The default configuration is the silica Coriolis setup
(a=5 µm, 532 nm, two probe beams z + x, forward + backward ports).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT / "data" / "literature" / "benchmarks.json"


# --------------------------------------------------------------------------- #
# validate
# --------------------------------------------------------------------------- #
def _run_gate(name: str, fn) -> bool:
    try:
        out = fn()
        ok = bool(out[0]) if isinstance(out, tuple) else bool(out)
        detail = ""
        if isinstance(out, tuple) and len(out) > 1 and not isinstance(out[1], dict):
            detail = f"(max err/margin {out[1]:.2e})"
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:44s} {detail}")
        return ok
    except Exception as exc:  # noqa: BLE001 — a gate that errors is a failed gate
        print(f"  [ERROR] {name:44s} {exc}")
        return False


def _validate(full: bool = False) -> int:
    from mieinfo.validation import golden, limits, convergence

    print("Validation gates (VALIDATION.md):")
    gates = [
        ("G-GOLD efficiencies", golden.check_efficiencies_golden),
        ("G-LIMIT energy / optical theorem", limits.check_energy_optical_theorem),
        ("G-LIMIT dipole information structure", limits.check_dipole_information_structure),
        ("G-LIMIT phi-parity", limits.check_symmetry_phi_parity),
        ("G-LIMIT x/y cross identity", limits.check_xy_cross_identity),
        ("G-CONV n_max", convergence.check_nmax_convergence),
        ("G-CONV angular grid", convergence.check_angular_grid_convergence),
    ]
    if full:
        gates.append(("G-GOLD information pattern (slow)", golden.check_information_pattern_golden))
    ok = all(_run_gate(n, f) for n, f in gates)
    print("ALL GATES PASSED" if ok else "VALIDATION FAILED")
    return 0 if ok else 1


# --------------------------------------------------------------------------- #
# shared: the default silica configuration
# --------------------------------------------------------------------------- #
def _silica_defaults(a_um: float = 5.0):
    from mieinfo.types import Sphere, Medium, AngularGrid
    from mieinfo.glmt.beam import PlaneWave
    from mieinfo.glmt.scatter import PlaneWaveProvider
    med = Medium(n=1.0, wavelength_vacuum_m=532e-9)
    sphere = Sphere(radius_m=a_um * 1e-6, m=1.46 + 0j)
    return PlaneWaveProvider(), sphere, PlaneWave(med), AngularGrid.full_sphere(400, 120)


def _optimize(config_path: str | None, out_dir: Path) -> int:
    from mieinfo.optimize import optimize_detection, Constraints
    provider, sphere, beam, grid = _silica_defaults()
    axes = {"x": np.array([1.0, 0, 0]), "y": np.array([0, 1.0, 0]), "z": np.array([0, 0, 1.0])}
    cons_single = dict(na_max=0.95, directions_allowed=("forward", "backward"),
                       schemes_allowed=("optimal", "self_homodyne"))

    result = {"config": config_path or "silica default (a=5um, 532nm)", "per_axis": {}}
    for name, n_hat in axes.items():
        r = optimize_detection(provider, sphere, beam, n_hat, Constraints(max_channels=1, **cons_single), grid)
        g, res = r.best
        result["per_axis"][name] = {"best_direction": g.direction, "best_NA": round(g.NA, 3),
                                    "best_scheme": g.lo_mode, "eta_q": round(res.eta_q, 4),
                                    "sensitivity": {k: float(v) for k, v in r.sensitivity.items()}}
        print(f"  n_hat={name}: {g.direction} NA={g.NA:.2f} {g.lo_mode}  eta_q={res.eta_q:.3f}")

    n_diag = np.array([1.0, 1, 1]) / np.sqrt(3)
    rmc = optimize_detection(provider, sphere, beam, n_diag,
                             Constraints(max_channels=2, beam_axes_lab=(axes["z"], axes["x"]), **cons_single), grid)
    mc = rmc.best_multichannel
    result["two_beam_(1,1,1)"] = {
        "channels": [f"{np.round(c.propagation_lab, 2).tolist()}:{c.geometry.direction}@NA{c.geometry.NA:.2f}"
                     for c in (rmc.best_channel_set or [])],
        "eta_q_total": round(mc.eta_q_total, 4), "fisher_total_rel": mc.fisher_total_rel}
    print(f"  two-beam n_hat=(1,1,1): {result['two_beam_(1,1,1)']['channels']}  "
          f"eta_q_total={mc.eta_q_total:.3f}")

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "optimize_result.json").write_text(json.dumps(result, indent=2))
    print(f"wrote {out_dir/'optimize_result.json'}")
    return 0


def _compare(out_dir: Path) -> int:
    from mieinfo.glmt.scatter import PlaneWaveProvider
    from mieinfo.literature.schema import load_benchmarks
    from mieinfo.literature.compare import compare_benchmark
    provider = PlaneWaveProvider()
    rows, npass = [], 0
    for b in load_benchmarks(str(BENCH)):
        r = compare_benchmark(provider, b)
        npass += int(r.within_tolerance)
        rows.append({"benchmark": r.benchmark_key, "predicted": r.predicted, "reported": r.reported,
                     "within_tolerance": r.within_tolerance, "note": r.discrepancy_note})
        print(f"  [{'PASS' if r.within_tolerance else 'note'}] {r.benchmark_key:38s} "
              f"pred={r.predicted:.3f} rep={r.reported:.3f}")
    print(f"G-LIT: {npass}/{len(rows)} within tolerance")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "comparison.json").write_text(json.dumps(rows, indent=2))
    print(f"wrote {out_dir/'comparison.json'}")
    return 0 if npass >= 2 else 1  # M4: >=2 independent results


def _report(out_dir: Path) -> int:
    print("== optimize ==");  rc1 = _optimize(None, out_dir)
    print("== compare ==");   rc2 = _compare(out_dir)
    print(f"\nSee docs/recommendation.md for the full synthesis (M5).")
    return rc1 or rc2


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="mieinfo", description="Information-pattern detection optimizer")
    sub = p.add_subparsers(dest="cmd", required=True)
    pv = sub.add_parser("validate"); pv.add_argument("--full", action="store_true")
    po = sub.add_parser("optimize"); po.add_argument("config", nargs="?", default=None)
    sub.add_parser("compare"); sub.add_parser("report")
    for pp in (po, sub.choices["compare"], sub.choices["report"]):
        pp.add_argument("--out", default=str(ROOT / "results"))
    args = p.parse_args(argv)

    if args.cmd == "validate":
        return _validate(args.full)
    out = Path(args.out)
    if args.cmd == "optimize":
        return _optimize(args.config, out)
    if args.cmd == "compare":
        return _compare(out)
    if args.cmd == "report":
        return _report(out)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
