# Pressure and virtual-temperature draft

This isolated header implements both Native `calc_pressure.F90` routines in
source order, each with eight required arguments. Actual JAX supports both;
its extracted Exner helper, branch-safe operands and scan/vector restructuring
are not imported. Full source and caller outlines precede implementation.
Compilation and all numerical evidence are pending.

Promote `calc_pressure.hpp` and replace the direct eight-argument
`f2k::calculate_thvm` call in `advance_clubb_core_module.cpp` with the public C++
name at the same call site. Two existing C++ driver calls use the same argument
order and may replace their matching boundary as part of the reviewed promotion.
The corresponding generated whole-call interface can then be removed when no
caller remains. Header-only implementation needs no new compiled source.
`init_pressure` is also complete, but its current hydrostatic/input-package
caller remains Native until that owner is ported; do not mistake that boundary
for a newly integrated initialization path.

All constants (`Cp`, `Lv`, `ep1`, `ep2`, `one`, `grav`, `p0`, `kappa`,
`zero_threshold`) and lower grid/type/array owners already exist. The only lower
scientific call is the live C++ `zt2zm_api(nzm,nzt,ngrdcol,gr,thvm,zero_threshold)`.
There is no new bridge, selected constant, type, accessor, global cache, stats
backend or error object to integrate.

`calculate_thvm` takes rank-2 View inputs/output directly and launches in the
output View's execution space. This keeps host/strided driver fields on the
host and core fields on the core execution space. All inputs must be accessible
from that space; current direct C++ callers have matching residency. The exact
Native three-term association and all five field-specific inputs are retained.

`init_pressure` uses concrete core Views because its grid components already
live there. Its generic adapter stages host/strided callers once through the
existing shared array owner and restores all four caller output allocations.
Interpolation supplies owned per-call momentum theta_v scratch. Only grid
component Views and dimension/value scalars are captured by kernels.

Initialization preserves five ordered phases: surface values on momentum
storagelevel1, interpolation, first thermodynamic Exner/pressure, remaining
thermodynamic recurrence and independent remaining momentum updates. The
recurrence runs ascending k within each independent column. This follows
Native storage indices; it does not reinterpret the fixed first level through
`grid_dir`/physical boundary metadata. Momentum updates depend on the completed
thermodynamic value at k-1, never a preceding momentum Exner value.

The first and recursive thermodynamic branches compare against
`epsilon(core_rknd)*thvm(current thermo)`. The momentum branch uses
`epsilon(core_rknd)*thvm(previous thermo)`, exactly as Native; the generalized JAX
helper's second-input threshold would differ. Epsilon is materialized outside
device kernels from the actual `real_t` trait. `grav/Cp` remains a Native constant
parameter. Logs, real powers, divisions and subtraction/product association
remain ordinary expressions. No extra floor, clipping, reciprocal rewrite,
FMA or speculative rounding adaptation is present.

Native input readers also use legacy single-column Fortran sequence association
and materialized reference-theta expressions. Those retained Native caller
contracts are documented in the outline; they do not create a new rank-1 C++
API. Future C++ input ports should supply rank-2 Views with the Native logical
shape. The separately excluded adaptive-grid caller does not expand this task.

Full Native hydrostatic equations, routine descriptions, comments, units and
argument groups remain. The copied claim that thermodynamic storagelevel1 is
above the lower boundary is replaced by its actual storage-index meaning.
`SOURCE_AUDIT.json` matches all 13 ordered assignment trees and three threshold
predicates; `CLEANUP_AUDIT.json` records a separate ownership/format review.
Neither establishes execution or BFB results. CPU/GPU/model outputs, pressure
initialization, near-equal theta thresholds, host/strided callers and mathematical
library differences remain untested; root owns the maintained CPU gate.
