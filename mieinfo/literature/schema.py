"""Literature schema (orchestrator-owned, frozen at M0 — INTERFACES.md §8).

Records validate against these models; W3 populates data/literature/*.json, W4c reads
them. Provenance is mandatory; unstated parameters are explicit null, never guessed
(LITERATURE.md §3). Complex refractive index is stored as [re, im] in JSON (pydantic has
no native complex); convention Im(m) >= 0 for absorption (PHYSICS.md §0).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, Any, Optional

from pydantic import BaseModel, BeforeValidator, ConfigDict, PlainSerializer


def _to_complex(v: Any) -> Optional[complex]:
    if v is None or isinstance(v, complex):
        return v
    if isinstance(v, (int, float)):
        return complex(v, 0.0)
    if isinstance(v, (list, tuple)) and len(v) == 2:
        return complex(float(v[0]), float(v[1]))
    if isinstance(v, str):
        return complex(v.replace(" ", "").replace("i", "j"))
    raise ValueError(f"cannot parse a complex refractive index from {v!r}")


# Stored/serialized as [re, im]; parsed from number | [re,im] | "n+kj".
Complex = Annotated[
    complex,
    BeforeValidator(_to_complex),
    PlainSerializer(lambda c: [c.real, c.imag], return_type=list, when_used="json"),
]


class Experiment(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)  # for the complex refractive_index

    key: str                                    # citation key, unique
    reference: str                              # full citation
    doi_or_arxiv: Optional[str] = None
    sphere_material: str
    sphere_radius_m: Optional[float] = None
    refractive_index: Optional[Complex] = None
    wavelength_m: Optional[float] = None
    collection_NA: Optional[float] = None
    collection_direction: Optional[str] = None  # forward/backward/split/cavity
    detection_scheme: Optional[str] = None       # homodyne/heterodyne/self-homodyne/split/imaging
    dof: Optional[str] = None                     # x/y/z/rotation/...
    reported_detection_efficiency: Optional[float] = None
    reported_imprecision: Optional[str] = None    # value + units as reported
    reported_backaction: Optional[str] = None
    notes: str = ""
    provenance: str                              # where in the paper (figure/eq/table)


class Benchmark(BaseModel):
    """An Experiment with enough parameters to reproduce numerically."""
    experiment_key: str
    run_config: dict                             # fully specifies a mieinfo RunConfig
    target_quantity: str                         # 'detection_efficiency' | 'imprecision_ratio' | ...
    target_value: float
    target_tolerance: float
    target_provenance: str


def load_experiments(path: str) -> list[Experiment]:
    return [Experiment(**d) for d in json.loads(Path(path).read_text())]


def load_benchmarks(path: str) -> list[Benchmark]:
    return [Benchmark(**d) for d in json.loads(Path(path).read_text())]
