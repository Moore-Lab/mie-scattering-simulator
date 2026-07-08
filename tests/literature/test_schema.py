"""literature/schema.py — record validation + the complex-index JSON encoding
(INTERFACES.md §8, LITERATURE.md §3)."""
from __future__ import annotations

import json

from mieinfo.literature.schema import Benchmark, Experiment, load_experiments


def test_complex_index_parse_and_json_roundtrip():
    e = Experiment(key="k", reference="r", sphere_material="silica",
                   refractive_index="1.45+0.001j", provenance="Fig 2")
    assert e.refractive_index == complex(1.45, 0.001)
    d = json.loads(e.model_dump_json())
    assert d["refractive_index"] == [1.45, 0.001]                 # stored as [re, im]
    assert Experiment(**d).refractive_index == complex(1.45, 0.001)  # reloads


def test_complex_index_accepts_pair_and_scalar():
    assert Experiment(key="k", reference="r", sphere_material="s", provenance="p",
                      refractive_index=[1.46, 0.0]).refractive_index == complex(1.46, 0.0)
    assert Experiment(key="k", reference="r", sphere_material="s", provenance="p",
                      refractive_index=1.46).refractive_index == complex(1.46, 0.0)


def test_optional_fields_default_null():
    e = Experiment(key="k", reference="r", sphere_material="silica", provenance="p")
    assert e.sphere_radius_m is None
    assert e.refractive_index is None
    assert e.notes == ""


def test_load_experiments(tmp_path):
    p = tmp_path / "exp.json"
    p.write_text(json.dumps([{"key": "a", "reference": "r", "sphere_material": "silica",
                              "provenance": "p", "refractive_index": [1.46, 0.0]}]))
    exps = load_experiments(str(p))
    assert len(exps) == 1 and exps[0].refractive_index == complex(1.46, 0.0)


def test_benchmark_fields():
    b = Benchmark(experiment_key="a", run_config={"radius_m": 5e-6}, target_quantity="detection_efficiency",
                  target_value=0.52, target_tolerance=0.05, target_provenance="Table 1")
    assert b.target_value == 0.52 and b.run_config["radius_m"] == 5e-6
