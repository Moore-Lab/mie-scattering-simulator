# Interface Contracts

**These are frozen after M0.** Workers implement against them and code against
stubs of contracts they consume. A worker that needs a change files a
`CONTRACT-CHANGE` blocker (`CONVENTIONS.md §5`); it does not unilaterally edit a
shared signature. The orchestrator owns `mieinfo/types.py`, the `FieldProvider`
protocol in `mieinfo/glmt/scatter.py`, and `mieinfo/literature/schema.py`.

Signatures below are the contract. Types are illustrative Python; keep names,
argument order, return shapes, and units exactly. Units are SI unless a field name
says otherwise (`_um`, `_deg`). Complex fields use the `e^{-iωt}` convention
(`PHYSICS.md §0`).

---

## 1. Core value types — `mieinfo/types.py` (orchestrator-owned)

```python
@dataclass(frozen=True)
class Sphere:
    radius_m: float
    m: complex                 # relative refractive index n_particle/n_medium

@dataclass(frozen=True)
class Medium:
    n: float = 1.0                       # vacuum around the sphere
    wavelength_vacuum_m: float = 532e-9  # DETECTION default (imaging beam); trap Medium uses 1064e-9
    @property
    def k(self) -> float: ...  # 2*pi*n/wavelength_vacuum_m
    # Sphere.m is relative to the medium AT that beam's wavelength: silica
    # m ~ 1.4607 at 532 nm, ~1.4496 at 1064 nm (dispersive).

@dataclass(frozen=True)
class AngularGrid:
    """Gauss-Legendre in cos(theta), uniform in phi. Carries solid-angle weights."""
    theta: np.ndarray          # (Ntheta,)  polar nodes (rad), from leggauss on cos(theta)
    phi: np.ndarray            # (Nphi,)    uniform azimuth nodes (rad)
    w_solid: np.ndarray        # (Ntheta, Nphi) dOmega weights; sum ~ 4*pi over full sphere
    @classmethod
    def full_sphere(cls, ntheta: int, nphi: int) -> "AngularGrid": ...
    # theta/phi meshed as (Ntheta, Nphi) with indexing="ij" wherever fields are evaluated.

@dataclass(frozen=True)
class VectorField:
    """Far-field scattered E on an AngularGrid, spherical components, arbitrary common units."""
    grid: AngularGrid
    E_theta: np.ndarray        # (Ntheta, Nphi) complex
    E_phi: np.ndarray          # (Ntheta, Nphi) complex
    def intensity(self) -> np.ndarray: ...   # |E_theta|^2 + |E_phi|^2

@dataclass(frozen=True)
class FieldDerivative:
    """dE_s/dr_j for j in x,y,z on the same grid."""
    grid: AngularGrid
    dE_theta: np.ndarray       # (3, Ntheta, Nphi) complex, axis0 = (x,y,z)
    dE_phi:   np.ndarray       # (3, Ntheta, Nphi) complex
```

## 2. Incident beam — `mieinfo/glmt/beam.py` [W1b]

```python
class IncidentBeam(Protocol):
    """Abstracts plane wave / paraxial Gaussian / Richards-Wolf focus."""
    medium: Medium
    polarization: np.ndarray            # Jones vector in the (x,y) transverse plane
    def focal_field(self, xyz: np.ndarray) -> np.ndarray:
        """E(r) at cartesian points xyz (...,3) -> (...,3) complex. Reference for BSC quadrature."""
    def waist_m(self) -> float: ...      # characteristic transverse scale (inf for plane wave)

class PlaneWave(IncidentBeam): ...
class GaussianParaxial(IncidentBeam):   # waist_m, propagation +z
    def __init__(self, medium: Medium, waist_m: float, polarization=(1,0)): ...
class RichardsWolfFocus(IncidentBeam):  # high-NA aplanatic focus of a plane wave
    def __init__(self, medium: Medium, NA: float, filling_factor: float=1.0, polarization=(1,0)): ...

# Each beam carries its OWN Medium, hence its own wavelength. Detection beams use
# a 532 nm Medium; the 1064 nm trap Medium is separate. In a beam's own formulae
# k_hat = +z. A beam's LAB-frame orientation is a rotation applied by the channel
# (see section 6): the beam computes patterns in its frame; the channel rotates
# s_hat into the lab frame to place the collection cone and to combine channels.
def lab_from_beam_frame(propagation_lab: np.ndarray,
                        polarization_lab: np.ndarray) -> "np.ndarray":  # 3x3 rotation
    """Rotation taking a beam-frame direction to the lab frame, given the beam's
    lab propagation axis (maps +z) and polarization (maps +x)."""
```

## 3. Beam-shape coefficients — `mieinfo/glmt/bsc.py` [W1b]

```python
@dataclass(frozen=True)
class BSC:
    n_max: int
    g_tm: np.ndarray   # indexed [n, m] over n=1..n_max, m=-n..n (packed; document packing)
    g_te: np.ndarray

def bsc_quadrature(beam: IncidentBeam, sphere_center_m: np.ndarray, n_max: int,
                   surface_radius_m: float) -> BSC: ...   # reference (slow)
def bsc_localized(beam: GaussianParaxial, sphere_center_m: np.ndarray, n_max: int) -> BSC: ...
def bsc_angular_spectrum(beam: RichardsWolfFocus, sphere_center_m: np.ndarray, n_max: int) -> BSC: ...
# All three return the SAME packing; a plane wave through any path must give the
# plane-wave weights (VALIDATION.md §3).
```

## 4. FieldProvider — the cross-track seam — `mieinfo/glmt/scatter.py` (orchestrator-owned protocol)

```python
class FieldProvider(Protocol):
    """The one interface info/ and detect/ depend on. Satisfied by mie.plane_wave now,
    glmt.scatter later. Everything downstream is written against THIS, not an implementation."""
    def field(self, grid: AngularGrid, sphere: Sphere, beam: IncidentBeam,
              r_s_m: np.ndarray, n_max: int | None = None) -> VectorField: ...
    def q_sca(self, sphere: Sphere, beam: IncidentBeam, n_max: int | None = None) -> float: ...

# Concrete implementations:
class PlaneWaveProvider(FieldProvider): ...   # wraps mie.plane_wave; r_s via analytic phase (PHYSICS 3.1)
class GLMTProvider(FieldProvider): ...        # wraps glmt.scatter; r_s via translated BSCs (PHYSICS 3.2)

def field_derivative(provider: FieldProvider, grid: AngularGrid, sphere: Sphere,
                     beam: IncidentBeam, r_s_m: np.ndarray,
                     method: str = "analytic", n_max: int | None = None) -> FieldDerivative:
    """method='analytic' (glmt.derivatives) or 'finite_difference' (validation).
    PlaneWaveProvider supports the exact analytic phase-gradient derivative (PHYSICS 3.1)."""
```

Stub for W2: a `PlaneWaveProvider` returning the prototype field, plus a
`field_derivative(..., method='analytic')` using the §3.1 phase gradient. W2 needs
nothing from GLMT to build and fully test.

## 5. Fisher information — `mieinfo/info/` [W2a]

```python
def info_density(deriv_q: tuple[np.ndarray, np.ndarray], grid: AngularGrid) -> np.ndarray:
    """(dE_theta_q, dE_phi_q) -> dF_q/dOmega on grid (|dE/dq|^2), arbitrary common units."""

def fisher_total(density: np.ndarray, grid: AngularGrid) -> float: ...   # integral over 4*pi

def combine_direction(deriv: FieldDerivative, n_hat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """dE/dq for q = r_s . n_hat: linear combo of Cartesian derivatives (PHYSICS 4.4).
    Returns (dE_theta_q, dE_phi_q)."""

@dataclass(frozen=True)
class InformationPattern:
    grid: AngularGrid
    density: np.ndarray        # dF_q/dOmega
    n_hat: np.ndarray
    f_total: float
```

## 6. Detection — `mieinfo/detect/` [W2b]

```python
@dataclass(frozen=True)
class CollectionGeometry:
    direction: str             # 'forward' | 'backward' | 'both' | 'split'
    NA: float                  # numerical aperture n*sin(alpha)
    apodization: str = "aplanatic"   # 'aplanatic' (sqrt(cos)) | 'none'
    lo_mode: str = "optimal"   # 'optimal' | 'gaussian' | 'self_homodyne' | 'quadrant'

@dataclass(frozen=True)
class DetectionResult:
    geometry: CollectionGeometry
    eta_q: float               # collected info / total info, incl. LO overlap
    s_imp_rel: float           # imprecision PSD in units of 1/f_total (relative)
    gamma_ba_rel: float        # backaction rate ~ Q_sca (relative units)
    sql_distance: float        # >=1; 1 == at the SQL

def cone_mask(grid: AngularGrid, geometry: CollectionGeometry) -> np.ndarray: ...  # bool (Ntheta,Nphi)

def apply_scheme(field: VectorField, deriv_q: tuple[np.ndarray, np.ndarray],
                 pattern: InformationPattern, geometry: CollectionGeometry,
                 q_sca: float) -> DetectionResult:
    """The scheme->metrics map (PHYSICS 4.2-4.3, 5). LO overlap for homodyne/self-homodyne/quadrant."""

# ---- Multi-channel: a real setup has several beams/ports (PHYSICS 4.5) ----
@dataclass(frozen=True)
class DetectionChannel:
    """One probe beam + one collection port, oriented in the lab frame."""
    beam: IncidentBeam                 # carries its own Medium (wavelength)
    propagation_lab: np.ndarray        # unit vector: beam's +z in lab frame
    polarization_lab: np.ndarray       # unit vector: beam's +x in lab frame
    geometry: CollectionGeometry
    name: str = ""
    independent: bool = True           # False if it shares a detector/LO with another channel

@dataclass(frozen=True)
class MultiChannelResult:
    per_channel: list[DetectionResult]           # aligned with the input channels
    eta_q_total: float                           # combined collected-info fraction of sum_c F_c(4pi)
    fisher_total_rel: float                      # sum_c eta_c * F_c(4pi), relative units
    sql_distance: float

def evaluate_channels(provider: FieldProvider, sphere: Sphere,
                      channels: list[DetectionChannel], n_hat: np.ndarray,
                      grid: AngularGrid) -> MultiChannelResult:
    """Per PHYSICS 4.5: evaluate each channel in its beam frame (rotate n_hat and
    s_hat via lab_from_beam_frame), then combine. Independent channels: Fisher info
    ADDS. Channels flagged independent=False must be combined with their correlation,
    not summed -- raise if that path is unimplemented rather than over-adding."""
```

## 7. Optimization — `mieinfo/optimize/` [W4b]

```python
@dataclass(frozen=True)
class Constraints:
    na_max: float
    directions_allowed: tuple[str, ...]      # e.g. ('forward','backward')
    schemes_allowed: tuple[str, ...]
    radius_range_m: tuple[float, float] | None = None   # if sphere size is a design variable
    fixed_sphere: Sphere | None = None
    beam_axes_lab: tuple[np.ndarray, ...] | None = None  # candidate probe axes to place channels on
    max_channels: int = 1                    # >1 optimizes a SET of beams/ports (PHYSICS 4.5)

@dataclass(frozen=True)
class OptResult:
    ranked: list[tuple[CollectionGeometry, DetectionResult]]   # best first (single-channel view)
    sensitivity: dict[str, float]     # d eta / d param for radius, NA, waist, m
    best: tuple[CollectionGeometry, DetectionResult]
    best_channel_set: list[DetectionChannel] | None = None     # when max_channels > 1
    best_multichannel: MultiChannelResult | None = None

def optimize_detection(provider: FieldProvider, sphere: Sphere, beam: IncidentBeam,
                       n_hat: np.ndarray, constraints: Constraints,
                       grid: AngularGrid) -> OptResult:
    """Single-channel: rank CollectionGeometry for one beam. Multi-channel
    (max_channels>1): also search over sets of DetectionChannel drawn from
    beam_axes_lab x directions x schemes, maximizing eta_q_total (PHYSICS 4.5).
    For an arbitrary n_hat, each axis is transverse to some beams and axial to
    others -- the optimizer should discover that the best set covers n_hat's
    components across complementary beams."""
```

## 8. Literature — `mieinfo/literature/schema.py` (orchestrator-owned)

```python
class Experiment(BaseModel):
    key: str                      # citation key, unique
    reference: str                # full citation
    doi_or_arxiv: str | None
    sphere_material: str
    sphere_radius_m: float | None
    refractive_index: complex | None
    wavelength_m: float | None
    collection_NA: float | None
    collection_direction: str | None      # forward/backward/split/cavity
    detection_scheme: str | None          # homodyne/heterodyne/self-homodyne/split/imaging
    dof: str | None                        # x/y/z/rotation/...
    reported_detection_efficiency: float | None
    reported_imprecision: str | None       # value + units as reported
    reported_backaction: str | None
    notes: str = ""                        # optional free-text (constraints on reproduction)
    provenance: str                        # REQUIRED — where in the paper (figure/eq/table)

class Benchmark(BaseModel):
    """An Experiment with enough parameters to reproduce numerically."""
    experiment_key: str
    run_config: dict                       # fully specifies a mieinfo RunConfig
    target_quantity: str                   # 'detection_efficiency' | 'imprecision_ratio' | ...
    target_value: float
    target_tolerance: float
    target_provenance: str

def load_experiments(path: str) -> list[Experiment]: ...
def load_benchmarks(path: str) -> list[Benchmark]: ...
```

## 9. Comparison — `mieinfo/literature/compare.py` [W4c]

```python
@dataclass(frozen=True)
class ComparisonResult:
    benchmark_key: str
    predicted: float
    reported: float
    within_tolerance: bool
    discrepancy_note: str

def compare_benchmark(provider: FieldProvider, benchmark: Benchmark) -> ComparisonResult: ...
```

## 10. CLI / Make targets — `mieinfo/cli.py`

```
mieinfo validate           # run validation.* gates; nonzero exit on failure
mieinfo optimize <config>  # OptResult for a RunConfig; writes JSON + figures
mieinfo compare            # run all benchmarks -> comparison table
mieinfo report             # assemble docs/recommendation.md inputs
```

`Makefile`: `install`, `test`, `validate`, `optimize`, `compare`, `report`.
`make validate` is the integration gate the orchestrator runs before every merge.
