# Third-moment integration

The draft supplies the two existing native public APIs in `clubb_kokkos`:

| Public API | Arguments | Behavior |
| --- | --- | --- |
| advance_xp3 | 38 | Simplified rt/thl/scalar third-moment updates, then u/v diagnosis |
| diagnose_xp3 | 28 | ADG1 or modified diagnostic ansatz for all third moments |

Argument names/order and required status are unchanged. There are no optionals.
Inputs/outputs use the established execution-space `Real1D`, `Real2D`, `Real3D`,
`grid` and live `f2k::stats_type` owners. The existing live core calls already
use those exact View/type contracts. Field/scalar units are documented in the
positional interfaces; native third moments are updated in caller storage.

When root is ready, move source/header to the normal core directory, add the
source to the existing core library and include its header in the core caller.
Replace the current `f2k::advance_xp3` and `f2k::diagnose_xp3` calls in
advance_clubb_core_module.cpp:916/927, then remove only those two whole-routine
interfaces from the existing bridge manifest. The current caller selects advance
when `l_advance_xp3` and non-ADG1; otherwise it selects diagnose. This draft does
not edit that selection or any shared file.

## Actual lower APIs and constants

Every lower owner already exists live:

- `grid_class.hpp`: rank-2 `zm2zt_api`, `zt2zm_api` and `ddzm`. Existing generic
  interpolation overloads handle scalar sections while retaining memory space.
- `Skx_module.hpp`: `Skx_func` and `xp3_LG_2005_ansatz`, with caller-supplied
  execution-space outputs. Both use column-first rank-2 Views.
- `f2k/stats_bridge.hpp` and `f2k/procedure_bridge.hpp`: sampled `stats_update`,
  rank-2 `stats_begin_budget` and `stats_finalize_budget`. Array-valued budget
  and dissipation expressions are evaluated explicitly before this boundary.
- `array_types.hpp`: `to_host_space`, `to_execution_space` and
  `from_execution_space`. These remain the only staging owner.
- `constants_clubb.hpp`: rt/thl/w tolerances, `w_tol_sqd`, `zero_threshold`,
  `grav`, `zero`, `one`, `one_half` and `three`.
- `parameter_indices.hpp`: `ixp3_coef_base`, `ixp3_coef_slope` and existing
  tunable-parameter ordering. `model_flags.hpp` already exposes `iiPDF_ADG1`.

Private solve-type values 1/2/3 mirror native module-private constants and remain
local. The native `l_predict_xp3=false` and `C_xp3_dissipation=1` remain private
routine parameters. There are no missing bridges/constants or declaration-only
stubs in this folder. Private scalar helper names also occur in another native
module; anonymous C++ namespace scope preserves their native privacy.

## Necessary execution adaptations

The public interface stays concrete. The private simplified routine accepts
field subviews through templates, so scalar-section calls keep source order
without whole-field copies. The two private production/accumulation expressions
are device-callable definitions within the same translation unit as their cell
kernel. Native scientific derivation/discretization comments are retained.

Columns/levels are independent in this module; no vertical recurrence or band
matrix solve exists. Kokkos kernels flatten cells with adjacent columns first.
Native 1-based indices convert to zero-based View accesses. The private interior
remains exactly native levels `1:nzt-1`, with `kp1=min(k+1,nzt)`. The upper
boundary remains the last storage level `nzt` in either grid direction; it is
not replaced by a grid-direction boundary index.

Scalar thresholds are copied once per public parent through `to_host_space`
for typed scalar lower arguments. In diagnostic scalar loops, the established
array adapter aliases compatible output sections or allocates/restores a
compatible field for the typed ansatz call; input copying is disabled because
that argument is native intent(out). Kernels capture only component Views,
including the grid spacing View rather than the owning grid object.

All scratch is call-owned. Local fields returned by existing grid functions
retain the returned owned allocation directly; they do not allocate an unused
initial array. Other scratch is allocated without initialization when completely
written before reading. Native production/accumulation zero initialization and
non-ADG1 width initialization remain explicit. The old third-moment snapshot is
allocated only for the predictive/LMM branch that reads it.

The Fortran stats backend remains unchanged. Native empty scalar budget names,
predictive-only begin/finalize gates and all update ordering remain visible.
This native module has no fatal diagnostics, random draw, host I/O or retained
mutable state.

## Validation handoff

The independent source audit is recorded, but no build or numerical run was
performed. Root must compile and run the maintained whole-model harness after
promotion. CPU Fortran/CPU Kokkos at zero tolerance is the regular gate; GPU
comparisons are occasional checkpoints. Exercise both parent selections and
relevant scalar/PDF/LMM settings before making any broader behavior claim.
Keep the source's disabled predictive branch explicitly untested until a
configuration actually exercises it. No new isolated numerical suite is required.
