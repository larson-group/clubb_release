# Mixing-length draft integration

This isolated draft translates the complete supported [mixing-length module](../../../src/CLUBB_core/mixing_length.F90). It is outside the build. Compilation and numerical verification are pending; the current source deliberately names helper APIs that have not been promoted yet.

## Start integration

Read [SOURCE_OUTLINE.md](SOURCE_OUTLINE.md) and [SOURCE_AUDIT.json](SOURCE_AUDIT.json), then promote [mixing_length.hpp](mixing_length.hpp) and [mixing_length.cpp](mixing_length.cpp) into `clubb_kokkos/src/CLUBB_core/` after their dependencies below are available. Add the source through the existing CMake owner. Replace the live `f2k::calc_Lscale` call with `calc_Lscale` and remove that whole-routine bridge at promotion time. The parent owns those shared changes; none were made in this draft.

The Native routine inventory and order are preserved: `calc_Lscale` (52 arguments), `set_Lscale_max` (5), private `compute_mixing_length` (20), `calc_Lscale_directly` (27), `diagnose_Lscale_from_tau` (39), private `compute_coherent_mixing_length` (23), and private scalar `sat_mixrat_liq_flatau` (2). The four Native public procedures are the four header APIs. JAX provides all three methods and seven routines; no method or option was excluded.

## Dependencies before compiling

The [saturation draft](../saturation/) is being prepared separately. `mixing_length.cpp` includes `saturation.hpp` and calls the Native-named device scalar `sat_mixrat_liq_api(real_t pressure, real_t temperature, int saturation_formula)` inside the legacy parcel trajectories. A host f2k call cannot run there. Promote the device-capable public saturation owner before enabling those branches; do not substitute the local Flatau helper, because method1 must honor the selected saturation formula.

Two initial first-layer saturation batches retain the explicit boundary `f2k::sat_mixrat_liq_api(nzt, ngrdcol, gr, p_in_Pa, tl_par_1, saturation_formula, start_index)`. `start_index` is a Native one-based thermodynamic level. The existing transient array overload drops that Native optional input, so it needs an overload selecting `sat_mixrat_liq_2D` with `start_index_in` supplied explicitly. Alternatively, once the saturation owner is promoted, replace those two callsites with its equivalent public View overload. Preserve the distinct Native scalar/2D behavior, including `I_sat_sphum`.

[BRIDGE_CONTRACTS.json](BRIDGE_CONTRACTS.json) records exact Native modules, source files, argument names/order, kinds, ranks, intents, optional status, and result ranks. The remaining calls use these public Native surfaces:

| Draft call | Native selected definition | Contract |
| --- | --- | --- |
| `f2k::calc_Ri_zm` | `advance_helper_module::calc_Ri_zm` | Seven required arguments; two rank2 inputs and rank2 output. |
| `f2k::length_check` | `numerical_check::length_check` | Five required arguments; three rank1 inputs and mutable `err_info`. |
| `f2k::smooth_max` | `smooth_max_sclr_idx` | `(nz, ngrdcol, scalar, rank2, scalar)` returns owned HostSpace rank2 values. |
| `f2k::smooth_max` | `smooth_max_array_scalar` | `(nz, ngrdcol, rank2, scalar, scalar)` returns owned HostSpace rank2 values. |
| `f2k::smooth_max` | `smooth_max_array_1D_scalar` | `(ngrdcol, rank1, scalar, scalar)` returns owned HostSpace rank1 values. |
| `f2k::smooth_min` | `smooth_min_sclr_idx` | `(nz, ngrdcol, scalar, rank2, scalar)` returns owned HostSpace rank2 values. |
| `f2k::smooth_min` | `smooth_min_array_scalar` | `(nz, ngrdcol, rank2, scalar, scalar)` returns owned HostSpace rank2 values. |
| `f2k::smooth_min` | `smooth_min_arrays` | `(nz, ngrdcol, rank2, rank2, scalar)` returns owned HostSpace rank2 values. |
| `f2k::smooth_heaviside_peskin` | `advance_helper_module::smooth_heaviside_peskin` | Four required arguments; rank2 input and result. |

The smooth-min/max selected definitions belong to `advance_helper_module`; its public generics retain the Native call names. None of these new contracts has optional arguments. Existing timer, debug metadata and stats boundaries remain in use. Helpers may instead be replaced with their public C++ implementations when their owners are promoted.

Extend the selected Native model constants in the existing generator owner with `Lscale_direct=1` and `Lscale_coherence=3`. `Lscale_tau=2`, the physical constants, and all mixing-length parameter indices already exist. Do not hard-code duplicate parameter-index constants in this draft.

## Execution and state ownership

Independent stages use column-first LayoutLeft execution-space Views, with adjacent cell threads accessing adjacent columns. Local scratch belongs to each call. Component View aliases keep the owning PDF/grid types out of kernel captures. Internal `i`, `j`, `k` and boundary metadata retain Native one-based values; only View accesses subtract one. Public C++ parameter indices are already zero based and are passed unchanged, including parameter-column subviews.

Parcel entrainment, energy exhaustion, quadratic stopping distances, endpoint extension and aggregation remain ordered. Coherent launches execute in Native CPU/OPENACC_BFB order within each column, so crossing/arrival additions need no atomics across independent columns. Both grid orientations and all three source overwrite choices remain present. The first low-coherence layer contributes before the early exit; downward arrival sums exclude the upper-boundary launch as Native does.

Error detection uses integer GPU reductions and copies failure diagnostics only when needed. Native debug gates, all-column or per-column status behavior, field-specific dumps, and return order remain visible. Diagnostic values use local streams with `max_digits10`; no shared stream formatting or global scratch is retained. Native `ERROR STOP` prints its source message before a C++ exception unwinds to the standalone error boundary.

## Arithmetic and BFB status

No experimentally justified FMA, forced reassociation or alternative math-library workaround has been introduced. Integer powers become explicit products; real powers use ordinary `Kokkos::pow`. Native `SQRT`, `EXP` and `SIGN` use the corresponding Kokkos operations; `copysign` preserves negative-zero SIGN behavior. MIN/MAX ternaries choose the second argument on equality, retaining the Native behavior already established in the grid phase. Cross-reference that existing evidence when promoting these semantic adaptations into the permanent BFB catalog.

The coherent cutoff is the unsuffixed Native `0.01`, so the draft uses `real_t(0.01f)` to preserve default-REAL rounding. Explicit `_core_rknd` literals retain the public `real_t` kind. The recorded Native build flags use `-g -O0 -Mbounds` without a default-real promotion flag. These are source-kind translations, not claims of newly demonstrated mismatch fixes.

The separate cleanup/source audit checked the full seven-routine inventory, signatures, call order, stats list, branch/formula/index/ownership details, diagnostics and formatting. **Build evidence: pending. Numerical evidence: pending.** No draft builds, tests or model runs were launched. The parent will integrate incrementally and use the maintained short `kokkos_vs_fortran` development loop, then the usual milestone comparisons. Passing has not been inferred from source review.

Owner format cleanup restores complete Native argument descriptions/units and
continuation/actual-intent groups in all seven definitions, three private
prototypes, four public declarations and35ordinary calls. Old undocumented
arguments in diagnose_Lscale_from_tau use the same named fields documented in
its parent. The l_avg_Lscale rationale is complete; index-convention boilerplate
is stated once, without inappropriate array notes in the scalar helper. Exact
executable-token equality against the stable pre-cleanup source is recorded in
`/tmp/clubb_mixing_format_token_audit.json`. This is source evidence, not runtime
validation; the distinct-agent full review is separate.
