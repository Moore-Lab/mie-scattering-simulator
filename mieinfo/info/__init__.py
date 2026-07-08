"""Fisher information about a sensed parameter, distributed over scattering direction.

The information radiation pattern is dF_q/dΩ ∝ |∂E_s/∂q|² (PHYSICS.md §4.1) — NOT the
intensity |E_s|²; it is reweighted by the displacement derivative. Everything here is a
function of the passed FieldProvider + field_derivative, so the plane-wave→GLMT swap is
a one-line change at the call site. Time convention e^{-iωt}.
"""
