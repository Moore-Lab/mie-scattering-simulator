"""Flask backend for the interactive information-pattern viewer.

Endpoints:
  GET /                 the viewer page
  GET /api/list         [{key, params, x}]  for every cached pattern
  GET /api/pattern?...  compute-or-load a pattern -> {theta, phi, r, params, x, cached}
"""
from __future__ import annotations

import json
import os

import numpy as np
from flask import Flask, jsonify, render_template, request
from scipy.ndimage import gaussian_filter1d

from ..detect.optics import CollectionGeometry, collection_efficiency
from ..glmt.beam import PlaneWave
from ..glmt.scatter import PlaneWaveProvider
from ..info.modes import information_pattern
from ..types import AngularGrid, Medium, Sphere

app = Flask(__name__)

CACHE = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "results", "patterns"))
os.makedirs(CACHE, exist_ok=True)

AXES = {"info_x": [1.0, 0, 0], "info_y": [0, 1.0, 0], "info_z": [0, 0, 1.0]}
QLABEL = {"intensity": "intensity |E_s|²", "info_x": "info: x-motion (transverse ∥pol)",
          "info_y": "info: y-motion (transverse ⊥pol)", "info_z": "info: z-motion (axial)"}


def _params(args) -> dict:
    return {
        "radius_um": round(float(args.get("radius_um", 5.0)), 4),
        "wavelength_nm": round(float(args.get("wavelength_nm", 532.0)), 2),
        "m": round(float(args.get("m", 1.46)), 4),
        "quantity": args.get("quantity", "info_z"),
        "ntheta": int(args.get("ntheta", 110)),
        "nphi": int(args.get("nphi", 110)),
    }


def _key(p: dict) -> str:
    return (f"a{p['radius_um']}um_l{p['wavelength_nm']}nm_m{p['m']}"
            f"_{p['quantity']}_{p['ntheta']}x{p['nphi']}")


def _compute(p: dict) -> dict:
    med = Medium(n=1.0, wavelength_vacuum_m=p["wavelength_nm"] * 1e-9)
    sph = Sphere(radius_m=p["radius_um"] * 1e-6, m=complex(p["m"], 0.0))
    grid = AngularGrid.full_sphere(p["ntheta"], p["nphi"])
    prov, beam = PlaneWaveProvider(), PlaneWave(med)
    if p["quantity"] == "intensity":
        dens = prov.field(grid, sph, beam, np.zeros(3)).intensity()
        eta = {}
    else:
        pat = information_pattern(prov, grid, sph, beam, np.zeros(3), np.array(AXES[p["quantity"]]))
        dens = pat.density
        eta = {d: round(collection_efficiency(pat, CollectionGeometry(d, 0.8)), 4)
               for d in ("forward", "backward")}
    dens = gaussian_filter1d(np.asarray(dens, float), 6, axis=0, mode="nearest")  # envelope for a clean surface
    r = dens / dens.max()
    return {"theta": grid.theta.tolist(), "phi": grid.phi.tolist(), "r": r.tolist(),
            "params": p, "label": QLABEL.get(p["quantity"], p["quantity"]),
            "x": round(float(med.k * sph.radius_m), 2), "eta_NA0.8": eta}


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/pattern")
def api_pattern():
    p = _params(request.args)
    if p["quantity"] not in QLABEL:
        return jsonify({"error": f"unknown quantity {p['quantity']}"}), 400
    path = os.path.join(CACHE, _key(p) + ".json")
    if os.path.exists(path):
        data = json.loads(open(path, encoding="utf-8").read())
        data["cached"] = True
        return jsonify(data)
    data = _compute(p)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f)
    # tiny sidecar for the list (avoids loading the big r array)
    with open(path + ".meta", "w", encoding="utf-8") as f:
        json.dump({"key": _key(p), "params": p, "x": data["x"], "label": data["label"]}, f)
    data["cached"] = False
    return jsonify(data)


@app.route("/api/list")
def api_list():
    items = []
    for f in sorted(os.listdir(CACHE)):
        if f.endswith(".json.meta"):
            try:
                items.append(json.loads(open(os.path.join(CACHE, f), encoding="utf-8").read()))
            except Exception:  # noqa: BLE001
                pass
    return jsonify(items)


def main(host: str = "127.0.0.1", port: int = 5000) -> None:
    print(f"mieinfo viewer -> http://{host}:{port}   (patterns cached in {CACHE})")
    app.run(host=host, port=port, debug=False)


if __name__ == "__main__":
    main()
