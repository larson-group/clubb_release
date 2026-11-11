# Advance-helper draft

[advance_helper_module.hpp](advance_helper_module.hpp) implements all 26 procedures
from [Native advance_helper_module.F90](../../../src/CLUBB_core/advance_helper_module.F90)
in Native definition order. It includes every scalar/array variant of the three
generics, optional boundary arguments, optional statistics and source-constant
false branches. No JAX omission is used to reduce the Native inventory.

The draft is source-audited and formatted. **Compilation, integration, execution
and numerical/BFB validation are pending.** It changes no production source,
CMake inputs, Native code or bridge generator.

## Promotion

Promote this header into `clubb_kokkos/src/CLUBB_core/` with its reached callers.
It is header-only: scalar smoothing definitions must remain visible to device
callers, and no new C++ compilation unit is needed. Replace reached `f2k::`
helper calls with their Native-named C++ equivalents. Remove obsolete bridge
contracts only after their last Native call sites are replaced.

The header already contains the exact `smooth_max_scalars` body from the
[partial PDF utility helper](../pdf_utilities/advance_helper_module.hpp), in its
canonical Native position. Promote this owner once; do not also copy the partial
helper into production. The full [PDF utility draft](../pdf_utilities/pdf_utilities.hpp)
already includes `advance_helper_module.hpp` by that canonical name. Its folder
was not changed during this task.

Dependencies are explicit:

- Production [grid_class.hpp](../../src/CLUBB_core/grid_class.hpp),
  [T_in_K_module.hpp](../../src/CLUBB_core/T_in_K_module.hpp),
  [array_types.hpp](../../src/CLUBB_core/array_types.hpp), constants and parameter
  indices already supply the grid/temperature operations and storage conversions.
- The [complete saturation draft](../saturation/saturation.hpp) supplies the
  Native-named rank-two `sat_mixrat_liq_api(nzt, ngrdcol, gr, p_in_Pa, T_in_K,
  saturation_formula)` call. Promote that dependency before compiling this
  header. No scalar GPU-to-Fortran saturation crossing is used.
- The existing [statistics API](../../src/CLUBB_core/stats_netcdf.hpp) supplies
  `stats_type`, sampled flags and `stats_update`. These are the only `f2k::`
  operations in the implementation; all occur on the host. No new physics
  bridge or model constant is required by this draft itself. Saturation's
  selector constants remain that draft's integration requirement.

This resolves the helper dependencies of the [mixing-length draft](../mixing_length/INTEGRATION.md):
all `smooth_min`/`smooth_max` ranks, `calc_Ri_zm` and `smooth_heaviside_peskin`.
Its separate `numerical_check::length_check` dependency remains explicit there.

Use the maintained model comparison after integration. The current acceptance
gate is zero-tolerance **CPU Kokkos versus CPU Fortran**. GPU comparisons are
occasional milestones; record their differences without adding speculative
special-function emulation or making GPU BFB a prerequisite for each promotion.

## Interfaces and ownership

There are 26 Native definitions, 11 generic dispatch overloads and 11 public
host/strided boundary overloads. The latter only use the shared array conversion
owner; they contain no alternate physics or new options. Core arrays are
column-first `LayoutLeft` execution-space Views. Smoothing results own their
storage and return to the array input's memory space; two-array overloads use
the first array's memory space. Conformability remains a caller contract.

LAPACK boundary indices retain Native **one-based values**, including optional
second-variable indices. Grid metadata is also Native one-based. View accesses
subtract one; tunable-parameter indices already use production zero-based
constants. Optional stats is `f2k::stats_type*`, defaulting to `nullptr`; a
present stats object is passed as `&stats` without transferring ownership.

Boundary helpers preserve lower1/upper1/lower2/upper2 overwrite order even for
repeated indices. The LHS work splits by band row, whose writes are independent;
the RHS updates stay in one work item. Native `ERROR STOP` messages become host
exceptions before matrix mutation. Generic adapters copy caller-owned InOut
storage back using `from_execution_space`.

Cell-independent work runs concurrently; each averaging-window search and its
accumulation retain Native per-cell order and grid direction. Search variables
retain Fortran's post-DO values when a search reaches the grid boundary. Profile
reductions use one work item with increasing-level accumulation and return one
host scalar. `vertical_integral` expands Native `SUM(field*rho_ds*dz)` into that
ordered loop; this has not yet established equality to a compiled Native SUM.
`calc_xpwp` writes only interior levels. Adapters retain caller-owned boundary
storage, and the rank-one launch preserves a zero-iteration Native DO on grids
with fewer than three momentum levels.

All scratch belongs to a call. Device bodies capture component Views or scalar
grid metadata, never owning grid/stats objects. In `pvertinterp`, unused Native
`l_error` becomes column-private; its bracket search, `k_upper=1` fallback and
absence of an error report are preserved. No new guard, clipping threshold,
undefined-input fallback, mutable global or derived-type registry is added.

## Source review and deliberate arithmetic choices

[SOURCE_OUTLINE.json](SOURCE_OUTLINE.json) records the full read-before-code
inventory, signatures, intent/optional/rank groups, branches and dependencies.
[SOURCE_AUDIT.json](SOURCE_AUDIT.json) records the separate manual cleanup/review;
[expression_audit.json](expression_audit.json) records contract, assignment-tree,
branch-predicate and ordered-call comparisons. Its source-only aid can be rerun
from the repository root:

```bash
python clubb_kokkos/drafts/advance_helper_module/source_audit.py
```

That command reads source/contracts and writes the draft audit record. It does
not compile or execute physics, and should not be installed as a production
module or expanded into an isolated oracle suite.

MIN/MAX ternaries select the second operand on ties, following the already
established Native convention in `BFB-003` of the
[BFB adaptation log](../../BFB_ADAPTATIONS.md). Integer squares/cubes are ordinary
products with Native grouping. The Richardson path deliberately preserves the
Native inverse multiplication rather than replacing it with division; that
choice and its BFB comment already exist in Native. No explicit FMA, forced
rounding, added floor or special-function workaround was introduced. Root owns
any extension of the permanent BFB log after integration/model evidence.

One pre-existing dormant Native defect is retained visibly: with
`l_Cx_fnc_Richardson_vert_avg` changed from its source `false` parameter to true,
`compute_Cx_fnc_Richardson` would compute a width average and then overwrite it
from uninitialized `Cx_fnc_Richardson_avg` (Native lines 911–923). The draft keeps
the false parameter and Native statement order and marks this exact block. This
port does not repair or claim support for enabling that broken Native branch.

Final owner cleanup preserves every executable token while restoring Native
continuation and actual intent groups in all26 definitions,24 prototypes and
seven science calls (stats remains InOut). Two matching boundary-adapter argument
lists were formatted the same way. Copied grid/stencil comments and one orphan
fragment were corrected; numerical behavior remains untested. Exact token
evidence: `/tmp/clubb_advance_helper_format_token_audit.json`.
