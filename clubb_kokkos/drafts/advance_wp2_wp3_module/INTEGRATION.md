# Vertical-velocity moment advance draft

This draft ports the complete supported
[`advance_wp2_wp3_module.F90`](../../../src/CLUBB_core/advance_wp2_wp3_module.F90)
source: the public `advance_wp2_wp3` and all 22 private routines, in native order.
The [public declaration](advance_wp2_wp3_module.hpp) retains all 68 required
argument names and positions. The [implementation](advance_wp2_wp3_module.cpp)
includes the native PDF choices, grid directions, stepping options, boundary
conditions, budgets, clipping, assertions and diagnostic dumps.

The files are isolated from the build. Compilation and numerical verification
are pending. Promote the reviewed source/header into the existing core owner,
establish the dependencies below, then replace the explicit parent f2k call in
`advance_clubb_core`. Remove the obsolete whole-parent bridge through its current
generator owner. Do not add the draft directory to production include paths.

## Lower calls and shared owners

Use the existing C++ grid operators, `Skx_func` and `band_solve`. The complete
[advance-helper draft](../advance_helper_module/advance_helper_module.hpp) supplies
`calc_wp3_on_wp2` and `vertical_avg`; promote those public APIs before this parent
is compiled. Native statistics remain under the existing `f2k::stats_type` owner.

The public native routines below remain explicit f2k calls. Their exact argument
order, kinds, ranks, dimensions, intents and optional status are recorded in
[DEPENDENCIES.json](DEPENDENCIES.json). Extend the existing bridge manifest and
generator; no new physics interface or global storage is required.

| Native module | Routine | Arguments | Optional argument |
| --- | --- | ---: | --- |
| diffusion | diffusion_zm_lhs | 10 | None |
| diffusion | diffusion_zt_lhs | 10 | None |
| mean_adv | term_ma_zm_lhs | 7 | None |
| mean_adv | term_ma_zt_lhs | 10 | None |
| fill_holes | fill_holes_vertical_api | 10 | None |
| fill_holes | fill_holes_wp2_from_horz_tke | 8 | None |
| clip_explicit | clip_variance | 9 | threshold_hi |
| clip_explicit | clip_skewness | 10 | None |
| stats_netcdf | stats_begin_budget_2d | 3 | None |
| stats_netcdf | stats_update_budget_2d | 3 | None |
| stats_netcdf | stats_finalize_budget_2d | 4 | l_count_sample |

`clip_variance` receives its optional upper threshold, `wp2_max`, at the native
call site. The budget wrappers retain the existing generic names
`stats_begin_budget`, `stats_update_budget` and `stats_finalize_budget`. The live
finalize wrapper currently exposes only three arguments; it must also support
the native optional `l_count_sample` so the two horizontal-variance budget calls
can pass `false`. Preserve absence at ordinary three-argument calls.

The existing `stats_update`, `var_on_stats_list` and
`clubb_at_least_debug_level_api` interfaces are reused. The rank-1 condition
diagnostic explicitly converts its View handle to `f2k::Real1D` to select the
existing rank-1 overload; the statistics owner stages the data.

Select these native compile-time values through the existing generated
`model_flags.hpp` owner before compiling, rather than introducing shadow flags:

- `l_explicit_turbulent_adv_wp3 = false`, native `model_flags.F90:120`.
- `l_force_descending_solves = false`, native `model_flags.F90:194`.

The private `clip_wp2 = 12`, `l_crank_nich_diff = false` and 2/3/5 diagonal
counts belong to this native module and remain private constants here.

## Layout, execution and diagnostics

- Physical fields retain `(column,level)` LayoutLeft storage. Local coefficient
  and matrix arrays keep native `(band,column,level)` shapes while lower calls
  use those contracts. The existing C++ band-solver boundary alone packs to
  `(column,band,row)` and restores the mutated LHS afterward. This is the
  existing layout adapter, with no scientific change.
- Forced descending solves reverse both matrix axes and the RHS, plus the old
  solution when needed. Owned reverse buffers preserve native whole-array
  assignment semantics. Only the final solution is reversed back, as in native.
- Native grid and loop metadata remain one-based within calculations; View
  accesses subtract one. Tunable-parameter constants already use the established
  zero-based C++ parameter axis. Mean-advection weight storage retains its native
  signed third-axis interpretation through the current grid and bridge owners.
- Independent column/level/band work uses Kokkos kernels. Local scratch belongs
  to each invocation. Component Views and scalar metadata are captured instead
  of owning UDTs. No f2k, statistics call, host read or print occurs in a kernel.
- C11/C16 assertions use per-column integer reductions and failure flags, then
  print failed cells in native level/column order. Both assertion checks run
  before the native all-column fatal return. NaN comparisons are unchanged.
- Native debug and fatal diagnostics stage only at their host boundary, retain
  field and label order, print arrays in Fortran linear order, and use
  `max_digits10`. Solver failure returns before prognostic solution assignment.
- The correlation-derived `wp2` floor evaluates each native MAX candidate once
  into a distinct const value. A local View saves pre-limit values for the debug
  warning, followed by the native ceiling at one. No scientific temporary is
  reused to hide another quantity.
- Native debug-level-3 `vertical_avg` passes rank-2 fields to an explicit rank-1
  dummy. The draft packs the first `nzm` column-fastest elements into rank-1
  Views, preserving native sequence association rather than substituting a
  per-column diagnostic.
- Array-valued statistics expressions are materialized in local Views directly
  at their native call sites. The backend retains budget ordering and the
  `l_count_sample=false` semantics described above.

## Supported scope and source notes

The two post-solve sponge blocks remain outside the established JAX/Kokkos
supported scope. The live C++ driver already rejects sponge-enabled cases;
retain that guard when promoting this module. The precise omitted boundary is
marked in the implementation. No sponge state, placeholder profile or new
global owner is introduced. Solver-method support remains the existing
`band_solve` policy.

Native `wp23_solve::wp3_pr3_lhs` omits INTENT but only reads that rank-3 field in
the sampled pressure budget. The draft's recorded parser contract is input;
neither native source nor the shared parser was modified. Source expressions
and branch asymmetries remain authoritative, including the distinct RHS TKE
coefficient indices.

The cleanup corrects unambiguous copied documentation errors: third-moment
units and density grid/units in `wp2_term_pr_dfsn_rhs`, the dimensional LHS
coefficient in `wp3_term_tp_lhs`, inverse-density/shear units and the `u'w'`
label in `wp3_term_pr_turb_rhs`, and two prose occurrences of `intent(in)`
where the native comment means “in”. Scientific expressions are unchanged.

No mismatch-driven arithmetic adjustment, explicit FMA or forced rounding has
been added. Integer powers and native MIN/MAX are ordinary language lowerings
with their expression trees recorded in the source audit. They still require
the normal numerical gate; this draft claims no BFB evidence or speedup.

## Review and validation handoff

[SOURCE_OUTLINE.json](SOURCE_OUTLINE.json) records the full-source read before
implementation. [SOURCE_AUDIT.json](SOURCE_AUDIT.json) distinguishes the final
source review from pending build/runtime evidence. The expression, parent
expression and call audit artifacts record 23 signatures, 233 explicitly
indexed scientific/budget assignment trees, 100 ordered lower calls with 634
arguments, and 171 kernel capture/local-variable checks. Whole-array, matrix,
diagnostic, warning and scope adapters received a separate source read. An
independent reviewer also read the complete source and target.

After integration, use the maintained short `multi_col_stats` comparison with
CPU Fortran and CPU Kokkos at zero tolerance. Longer duration and standard/all
statistics checks belong to larger integrations or targeted diagnoses. GPU
comparisons are occasional milestones; document differences rather than
silently changing comparison thresholds. Record any deliberate BFB arithmetic
adaptation in the existing [BFB log](../../BFB_ADAPTATIONS.md).
