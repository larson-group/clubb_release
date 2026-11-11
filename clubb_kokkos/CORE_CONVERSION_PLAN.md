# CLUBB core conversion plan

This records the supported C++/Kokkos core conversion completed at `ac37a3898`.
Its short CPU comparisons and full 19-case native-duration/native-timestep
matrix with standard stats pass strictly bit-for-bit (BFB). Its CUDA build and
default and active-hybrid GPU checkpoints complete. These results belong to
that core-phase executable pair, before the later driver ports.
The [driver conversion plan](DRIVER_CONVERSION_PLAN.md) tracks their current
completed integration and recorded numerical acceptance. The inventory below distinguishes
C++ bodies, retained Native services and excluded paths; source completion does
not establish BFB agreement on every branch.

## Scope and acceptance

- Use current [Fortran](../src/CLUBB_core/) as the source authority. Use the
  [JAX plan](../clubb_jax/JAX_CONVERSION_PLAN.md) and implementation to resolve
  supported-path decisions; retain Fortran numerical/error semantics when JAX differs.
- Port the supported timestep path and its required initialization/helpers.
  Keep source names, routine/argument/call order, comments and branch coverage.
  Every completed file needs a separate source audit and cleanup pass.
- Use the live call tree to select compatible groups of source-audited ports.
  Multiple files can be replaced and tested together, including in parallel.
  Accept the full group when it is CPU BFB; no separate model gate is required
  for each file. If a group fails, narrow the added files to identify the source
  before debugging its first CPU difference. Keep unported lower calls explicit;
  avoid dependency-first detours and broad isolated oracle suites.
- Draft and source-audit independent files in parallel under [drafts](drafts/)
  while builds and comparisons run. Keep shared integration changes owned by
  the integrating worker, and distinguish draft review from compiled/model
  evidence. Promote reviewed groups into the live call tree.
- The core phase retained Native statistics. The later
  [statistics API](src/CLUBB_core/stats_netcdf.hpp) now owns its C++ state and output.
  Optional BFB PDF matrix calls stay visibly in [f2k](src/f2k/); whole supported
  scientific routines are no longer forwarded.
- Keep column-first `LayoutLeft` device fields and call-owned scratch. Preserve
  vertical recurrences. Explain any necessary layout or language adaptation at
  its code boundary.
- Require CPU Kokkos BFB with the CPU Fortran reference in the recorded build
  before advancing live integration. Use GPU comparisons occasionally at
  milestones, recording the newly integrated piece where GPU results change,
  even slightly. Do not investigate GPU numerical causes unless requested.
  Diagnose CPU gate failures without changing tolerances. No
  optimized-performance or other-compiler claim follows from Debug validation.
- The agreed unsupported driver features are sponge damping, input-fields mode,
  the generalized-grid test, adaptive-grid methods, GFDL activation, BUGSrad and
  SILHS radiation. The fixed `l_add_dycore_grid` constructor is supported.
  The [solver wrapper](src/CLUBB_core/matrix_solver_wrapper.cpp) explicitly rejects
  LAPACK and BiCGStab solve methods. LAPACK condition diagnostics are also
  unsupported: requesting `rcond` returns fatal `err_info` before any solve.
- Radiation, microphysics, SILHS and other outer-driver packages were retained
  Native dependencies at the core-phase checkpoint. Their supported science
  now executes in C++; the later driver gate is tracked separately. No upstream
  runner selector integration is required; keep using `run_scm.py -exe`.

## Validation

Follow the [maintained comparison workflow](README.md#compare-with-fortran).
Use short multicolumn checks during development, native timestep/duration and
standard stats for final integration, and all stats occasionally after a large
piece or when budgets or diagnostics need coverage. Record deliberate arithmetic
choices made for BFB in the [adaptation log](BFB_ADAPTATIONS.md).
The regular numerical gate is CPU Fortran versus a separate CPU Kokkos build,
with zero tolerances. GPU comparisons are milestone evidence: record the newly
integrated piece where results changed, without tracing GPU numerical causes
unless requested. GPU BFB is not required before the next port.
Independent reviewed drafts may be integrated/tested in parallel frozen source
snapshots, each with separate builds and outputs. Review and combine passing
changes, then recheck the combined model with the same comparison.
The model comparison is the primary validation; small unit checks
cover specific edge cases it cannot exercise. CTest runs the existing checks documented
in the [README](README.md#derived-types-and-temporary-calls).

Before completing the whole phase, resolve every inventory row below, audit the
supported call graph for remaining scientific f2k delegations, run the full case
matrix and review every untested or deliberately excluded branch. A data-type
mirror alone does not complete a module with procedures.

## File inventory

“Ported/audited” describes that file’s supported routines, not every transitive
Native caller or every runtime branch. Retained backends/services and explicit
exclusions are recorded below. Passing groups remain scoped to their recorded
binaries and options; the final matrix below verifies the frozen core-phase
pair at `ac37a3898`, not the later driver implementation. Current caller and
service notes reflect the live source.

| Fortran file | C++ owner | State / remaining work |
| --- | --- | --- |
| [LY93_pdf.F90](../src/CLUBB_core/LY93_pdf.F90) | [LY93_pdf.hpp](src/CLUBB_core/LY93_pdf.hpp) | All three routines promoted and independently source-reviewed. Selector 6 requires input-fields mode, which the standalone rejects; algorithm runtime/BFB remains unverified. |
| [Nc_Ncn_eqns.F90](../src/CLUBB_core/Nc_Ncn_eqns.F90) | [Nc_Ncn_eqns.hpp](src/CLUBB_core/Nc_Ncn_eqns.hpp) | All six routines source-ported/audited, with device scalar and same-memory-space elemental View overloads. Metadata/setup isolated CPU cases and active cloudy all-stats pass; combined 19-case CPU/debug2 and four-case all-stats/debug2 checks pass. |
| [Skx_module.F90](../src/CLUBB_core/Skx_module.F90) | [Skx_module.cpp](src/CLUBB_core/Skx_module.cpp) | Ported/audited; direct or case comparisons pass. |
| [T_in_K_module.F90](../src/CLUBB_core/T_in_K_module.F90) | [T_in_K_module.hpp](src/CLUBB_core/T_in_K_module.hpp) | Ported/audited; direct or case comparisons pass. |
| [adg1_adg2_3d_luhar_pdf.F90](../src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90) | [adg1_adg2_3d_luhar_pdf.cpp](src/CLUBB_core/adg1_adg2_3d_luhar_pdf.cpp) | All nine routines promoted/source-reviewed, including independent review and device Luhar/root adaptations. ADG1 passes the 19 short native-timestep CPU cases. ADG2/Luhar selectors 2/3 remain runtime-gated by input-fields mode. |
| [advance_clubb_core_module.F90](../src/CLUBB_core/advance_clubb_core_module.F90) | [advance_clubb_core_module.cpp](src/CLUBB_core/advance_clubb_core_module.cpp) | Ported/audited; direct or case comparisons pass. |
| [advance_helper_module.F90](../src/CLUBB_core/advance_helper_module.F90) | [advance_helper_module.hpp](src/CLUBB_core/advance_helper_module.hpp) | All 26 routines/adapters promoted and independently source-reviewed. Reached stability/shear/splat/radiation and moment-parent calls pass recorded CPU groups, including targeted moist/shear/stability options. Other generic/reduction and compile-time paths retain their recorded runtime limits; profile SUM arithmetic is not universally verified. |
| [advance_windm_edsclrm_module.F90](../src/CLUBB_core/advance_windm_edsclrm_module.F90) | [advance_windm_edsclrm_module.cpp](src/CLUBB_core/advance_windm_edsclrm_module.cpp) | All 6 routines source-ported/audited; combined moment/lower-physics group passes19 short CPU cases and full 360-step/native 60s/standard-stats BOMEX. Predictor=false and passive-scalar paths also pass targeted lane checks. |
| [advance_wp2_wp3_module.F90](../src/CLUBB_core/advance_wp2_wp3_module.F90) | [advance_wp2_wp3_module.cpp](src/CLUBB_core/advance_wp2_wp3_module.cpp) | All 23 routines source-ported/audited; combined group passes19 short CPU cases and full 360-step BOMEX. No additional arithmetic adaptation; Native compile-time option coverage remains scoped to recorded lane checks. |
| [advance_xm_wpxp_module.F90](../src/CLUBB_core/advance_xm_wpxp_module.F90) | [advance_xm_wpxp_module.cpp](src/CLUBB_core/advance_xm_wpxp_module.cpp) | All 18 routines source-ported/audited; BFB-017 records the verified main-diagonal contraction adaptation. Combined group passes19 short CPU cases and full 360-step BOMEX. |
| [advance_xp2_xpyp_module.F90](../src/CLUBB_core/advance_xp2_xpyp_module.F90) | [advance_xp2_xpyp_module.cpp](src/CLUBB_core/advance_xp2_xpyp_module.cpp) | All 16 routines source-ported/audited; BFB-016 records the verified main-diagonal contraction adaptation. Combined group passes19 short CPU cases and full 360-step BOMEX. Multiple-LHS and compile-time-disabled paths are not runtime claims. |
| [advance_xp3_module.F90](../src/CLUBB_core/advance_xp3_module.F90) | [advance_xp3_module.cpp](src/CLUBB_core/advance_xp3_module.cpp) | All 5 routines independently source-audited/compiled; live diagnose path passes19 CPU cases and active alternate-PDF all-stats check. Compile-time-disabled advancement paths remain runtime-unexercised. |
| [array_index.F90](../src/CLUBB_core/array_index.F90) | [array_index.hpp](src/CLUBB_core/array_index.hpp) | Complete six-field scalar-index type mirror; no Native procedures exist in this module. Source field order and one-based values are retained, and generated UDT roundtrip checks pass. |
| [calc_pressure.F90](../src/CLUBB_core/calc_pressure.F90) | [calc_pressure.hpp](src/CLUBB_core/calc_pressure.hpp) | Both routines source-ported/audited; three calculate_thvm callers pass19 short CPU cases. init_pressure has no newly ported caller/runtime evidence. |
| [calc_roots.F90](../src/CLUBB_core/calc_roots.F90) | [calc_roots.hpp](src/CLUBB_core/calc_roots.hpp) | All three Native APIs and required device scalar overloads promoted and independently source-reviewed. Isolated CPU/CUDA compilation passes; numerical coverage follows actual PDF callers, without a blanket complex-root/BFB claim. |
| [calendar.F90](../src/CLUBB_core/calendar.F90) | [calendar.cpp](src/CLUBB_core/calendar.cpp) | All five routines and immutable tables source-ported/audited; CPU compilation and tiny century-leap, signed elapsed-time rollover and fatal-day checks pass. C++ radiation and input-date interpretation now call this owner. |
| [clip_explicit.F90](../src/CLUBB_core/clip_explicit.F90) | [clip_explicit.cpp](src/CLUBB_core/clip_explicit.cpp) | All 6 routines source-ported/audited and reached parent calls switched; combined 19 short CPU cases and full 360-step BOMEX are byte-identical. |
| [clubb_api_module.F90](../src/CLUBB_core/clubb_api_module.F90) | [clubb_api_module.hpp](src/CLUBB_core/clubb_api_module.hpp) | The 14 Native bodies in the supported standalone surface are represented, including four grid and two derived-parameter wrappers in Native relative order. Five external host/single-column or unused wrappers remain outside this -exe phase. Reached standalone calls pass the combined short CPU gate. |
| [clubb_precision.F90](../src/CLUBB_core/clubb_precision.F90) | — | Compile-time kinds are mapped by clubb_real.hpp (core float/double), optional BFB-mode precision-selected LAPACK bindings and double model-time fields. Statistical buffers and time use double precision in the permanent C++ stats owner. No runtime algorithm remains in this module. |
| [code_timer_module.F90](../src/CLUBB_core/code_timer_module.F90) | [code_timer_module.cpp](src/CLUBB_core/code_timer_module.cpp) | All 16 bodies are ported in Native order: five public APIs, CPU registry/stack/summary helpers and common path/failure helpers. Both process-CPU-time and direct GPTL backends retain time_clubb.py-compatible output. The agreed private per-thread bookkeeping adds no GPU fences. All five timer bridge entries are removed. |
| [constants_clubb.F90](../src/CLUBB_core/constants_clubb.F90) | [constants_clubb.hpp](src/CLUBB_core/constants_clubb.hpp) | Ported/audited; direct or case comparisons pass. |
| [corr_varnce_module.F90](../src/CLUBB_core/corr_varnce_module.F90) | [corr_varnce_module.cpp](src/CLUBB_core/corr_varnce_module.cpp) | All seven routines source-ported/audited; immutable Native default tables and public parameter types surround the unchanged generated metadata schema. Metadata/setup isolated CPU cases pass; combined 19-case CPU/debug2 and four-case all-stats/debug2 checks pass. |
| [diagnose_correlations_module.F90](../src/CLUBB_core/diagnose_correlations_module.F90) | [diagnose_correlations_module.hpp](src/CLUBB_core/diagnose_correlations_module.hpp) | All ten active routines source-ported/audited; l_calc_w_corr=true retains the JAX unsupported-path gate. Two diagnostic-option model probes are rejected by both models at the same symmetry assertion, not accepted trajectories. Normal group cases pass in isolation; combined 19-case CPU/debug2 and four-case all-stats/debug2 checks pass. |
| [diffusion.F90](../src/CLUBB_core/diffusion.F90) | [diffusion.cpp](src/CLUBB_core/diffusion.cpp) | Both routines source-ported/audited and reached parent calls switched; combined 19 short CPU cases and full 360-step BOMEX are byte-identical. |
| [endian.F90](../src/CLUBB_core/endian.F90) | — | Explicitly omit this unused I/O utility from the supported C++ call graph: no JAX mirror and no external src references to its module, generic, constants or private bodies. Native definitions remain unchanged. |
| [err_info_type_module.F90](../src/CLUBB_core/err_info_type_module.F90) | [err_info_type_module.cpp](src/CLUBB_core/err_info_type_module.cpp) | All 5 lifecycle procedures, including fixed-width header formatting, now execute in C++. Context and header allocations remain case-owned; I2/I3/F6.2 edge checks compare with Native separately from whole-model numerical acceptance. |
| [error_code.F90](../src/CLUBB_core/error_code.F90) | [error_code.cpp](src/CLUBB_core/error_code.cpp) / [header](src/CLUBB_core/error_code.hpp) | All three status constants, ANY check and signed debug-level accessors are C++. Private per-host-thread diagnostic metadata follows Native threadprivate state; statistics use the same C++ accessor, with no Native metadata synchronization. Generic debug bridges are removed. The unused process/thread-prefix global has no C++ copy. |
| [file_functions.F90](../src/CLUBB_core/file_functions.F90) | [file_functions.cpp](src/CLUBB_core/file_functions.cpp) | Both formatted readers now execute in C++, preserving physical-record READ behavior, repeated/null values and caller-owned output. Short cloudy M-PACE and focused record checks cover this later input port. |
| [fill_holes.F90](../src/CLUBB_core/fill_holes.F90) | [fill_holes.cpp](src/CLUBB_core/fill_holes.cpp) | All 10 routines source-ported with owner source/cleanup audit; integrated group passes19 short CPU cases, targeted active-scalar all-stats and full 360-step BOMEX. Later C++ microphysics callers exposed a padded-column `ANY` scan defect; logical rank1–3 indexing fixes it, with a failing-before/passing-after guard-column regression, DO/DS/RICO20-step standard checks and the later full six-case KK/SILHS pass. Not every fill strategy is runtime-covered. |
| [grid_adaptation_module.F90](../src/CLUBB_core/grid_adaptation_module.F90) | [grid_adaptation_module.hpp](src/CLUBB_core/grid_adaptation_module.hpp) | Bounded setup_gr_dycore promoted and independently source-reviewed; it now calls the C++ read_grid_heights. This fixed-grid constructor preserves the existing l_add_dycore_grid option. Broader adaptive-grid/generalization algorithms remain excluded. |
| [grid_class.F90](../src/CLUBB_core/grid_class.F90) | [grid_class.cpp](src/CLUBB_core/grid_class.cpp) | All 22 construction/cleanup/height/weight/operator routines and read_grid_heights are C++, including file/schema/monotonicity validation. Earlier 19 short cases/full BOMEX and isolated initialization grid types 1/2/3 pass; the later input port has focused grid/error checks and a byte-identical short CPU gate. |
| [hydromet_pdf_parameter_module.F90](../src/CLUBB_core/hydromet_pdf_parameter_module.F90) | [hydromet_pdf_parameter_module.cpp](src/CLUBB_core/hydromet_pdf_parameter_module.cpp) | All 3 lifecycle routines source-ported/audited; generated collection schema remains authoritative. Isolated and combined 19 short CPU cases and existing UDT checks pass. |
| [index_mapping.F90](../src/CLUBB_core/index_mapping.F90) | [index_mapping.cpp](src/CLUBB_core/index_mapping.cpp) | All 5 routines source-ported/audited; combined 19 short CPU cases and full 360-step BOMEX are byte-identical. |
| [interpolation.F90](../src/CLUBB_core/interpolation.F90) | [interpolation.hpp](src/CLUBB_core/interpolation.hpp) | All 8 routines ported/audited; direct native CPU/GPU checks, including isolated quintic-enabled bodies, pass. Shared debug metadata remains native. |
| [lapack_interfaces.F90](../src/CLUBB_core/lapack_interfaces.F90) | [lapack_interfaces.hpp](src/CLUBB_core/lapack_interfaces.hpp) | Four precision-selected external-library interfaces serve PDF matrices only with KOKKOS_BFB=ON; default operations are custom C++. Expert solver interfaces have been removed. Raw library ABI stays inside f2k; it does not delegate a CLUBB scientific wrapper. |
| [lapack_wrap.F90](../src/CLUBB_core/lapack_wrap.F90) | [lapack_wrap.cpp](src/CLUBB_core/lapack_wrap.cpp) | Retained expert interfaces reject LAPACK with fatal err_info and an explicit message. Solver and condition-diagnostic f2k calls have been removed; rejection checks replace the former numerical LAPACK comparisons. |
| [matrix_operations.F90](../src/CLUBB_core/matrix_operations.F90) | [matrix_operations.hpp](src/CLUBB_core/matrix_operations.hpp) | Five source-used routines are C++, including required Cholesky/eigenvalue/correlation diagnostics and set_lower_triangular_matrix for the C++ correlation-file reader. KOKKOS_BFB=ON retains the original external LAPACK ABI service; the default uses a JAX-style factorization and Jacobi eigenvalue diagnostics. Three unused utilities remain outside the supported C++ call graph. The original four-routine group passes its recorded isolated, combined 19-case CPU/debug2 and four-case all-stats/debug2 checks; the setter belongs to the later driver port. |
| [matrix_solver_wrapper.F90](../src/CLUBB_core/matrix_solver_wrapper.F90) | [matrix_solver_wrapper.cpp](src/CLUBB_core/matrix_solver_wrapper.cpp) | All 5 solves/3 NaN helpers ported/audited. Supported LU paths retain Native comparison coverage; LAPACK methods and rcond requests return fatal err_info before solving. All four moment parents call these C++ wrappers; prior combined 19 CPU cases/full BOMEX pass. |
| [mean_adv.F90](../src/CLUBB_core/mean_adv.F90) | [mean_adv.cpp](src/CLUBB_core/mean_adv.cpp) | Both routines source-ported/audited; combined 19 short CPU cases and full 360-step BOMEX are byte-identical. |
| [mixing_length.F90](../src/CLUBB_core/mixing_length.F90) | [mixing_length.cpp](src/CLUBB_core/mixing_length.cpp) | All 7 routines independently source-audited/compiled; method3 passes19 short CPU cases; methods1/2 each pass two cases. Required length_check and timers are C++; stats and debug metadata now use permanent C++ owners. |
| [model_flags.F90](../src/CLUBB_core/model_flags.F90) | [model_flags.cpp](src/CLUBB_core/model_flags.cpp) | All three default/initialization/printing procedures promoted/source-audited around the authoritative generated 66-field schema. Mutable flags are case-owned; selectors are immutable. Printing uses caller-owned C++ streams; the unused integer-unit bridge is removed. Isolated initialization cases pass; combined short CPU checks pass byte-identically. |
| [mono_flux_limiter.F90](../src/CLUBB_core/mono_flux_limiter.F90) | [mono_flux_limiter.cpp](src/CLUBB_core/mono_flux_limiter.cpp) | All 7 routines source-ported with owner source/cleanup audit; combined 19 short CPU cases and full 360-step BOMEX are byte-identical. Not every limiter/compile-time choice is runtime-covered. |
| [new_hybrid_pdf.F90](../src/CLUBB_core/new_hybrid_pdf.F90) | [new_hybrid_pdf.hpp](src/CLUBB_core/new_hybrid_pdf.hpp) | Nine helper routines promoted/source-audited with device scalar/profile APIs. Isolated default and active PDF7 CPU cases, including scalar/all-stats branches, are byte-identical. BFB-019/020/021 record the verified arithmetic adaptations; combined short CPU checks pass byte-identically. |
| [new_hybrid_pdf_main.F90](../src/CLUBB_core/new_hybrid_pdf_main.F90) | [new_hybrid_pdf_main.cpp](src/CLUBB_core/new_hybrid_pdf_main.cpp) | The full driver is promoted/source-audited. Isolated PDF7 multicolumn/all-stats and predictor=false scalar checks pass byte-identically; default ADG1 does not exercise this path. Combined short CPU checks pass byte-identically. |
| [new_pdf.F90](../src/CLUBB_core/new_pdf.F90) | [new_pdf.hpp](src/CLUBB_core/new_pdf.hpp) | All nine routines promoted and independently source-reviewed, including execution-space complex-root limits and all coefficient branches. Final grouping/comment cleanup preserves executable tokens. Selector 4 remains runtime-gated by input-fields mode. |
| [new_pdf_main.F90](../src/CLUBB_core/new_pdf_main.F90) | [new_pdf_main.cpp](src/CLUBB_core/new_pdf_main.cpp) | All four routines promoted and independently source-reviewed, including all three advection-flag branches, skew clipping and coefficient packing. Selector 4 remains runtime-gated; no algorithm runtime/BFB claim. |
| [new_tsdadg_pdf.F90](../src/CLUBB_core/new_tsdadg_pdf.F90) | [new_tsdadg_pdf.hpp](src/CLUBB_core/new_tsdadg_pdf.hpp) | All four routines promoted and independently source-reviewed, preserving setter ties, epsilon branches and ordered host warning payload. Selector 5 remains runtime-gated by input-fields mode. |
| [numerical_check.F90](../src/CLUBB_core/numerical_check.F90) | [numerical_check.cpp](src/CLUBB_core/numerical_check.cpp) | All 12 routines source-ported/audited; four reached whole-check bridges removed. Four short debug2 CPU cases pass byte-identically; combined lifecycle group passes all 19 short CPU cases with debug level 2. Parameter validation calls C++ check_parameters_api; shared debug-level metadata now uses the permanent C++ accessor. |
| [parameter_indices.F90](../src/CLUBB_core/parameter_indices.F90) | [parameter_indices.hpp](src/CLUBB_core/parameter_indices.hpp) | Ported/audited; direct or case comparisons pass. |
| [parameters_tunable.F90](../src/CLUBB_core/parameters_tunable.F90) | [parameters_tunable.hpp](src/CLUBB_core/parameters_tunable.hpp) | All eight supported default/catalog/check/derived/init/private pack/unpack bodies promoted/source-audited. C++ owns parameter values, calculations and initialized-array namelist input. Two JAX-omitted tuner readers remain excluded. Isolated initialization cases pass; combined short CPU checks pass byte-identically. |
| [pdf_closure_module.F90](../src/CLUBB_core/pdf_closure_module.F90) | [pdf_closure_module.cpp](src/CLUBB_core/pdf_closure_module.cpp) | All 18 source-used routines ported/audited; three uncalled private interpolation helpers excluded with call-graph evidence. CPU checks pass19 short cases, active cloudy all-stats cases and full 360-step BOMEX. All six alternate algorithm calls now use permanent C++ APIs. Default ADG1 skips fusion; recorded isolated PDF7 cases exercise it. Selectors 2–6 remain input-fields-gated, and combined short CPU checks pass byte-identically. |
| [pdf_parameter_module.F90](../src/CLUBB_core/pdf_parameter_module.F90) | [pdf_parameter_module.cpp](src/CLUBB_core/pdf_parameter_module.cpp) | All 7 JAX-supported lifecycle/pack/print routines source-ported/audited; isolated and combined 19 short CPU cases pass. Two Native-only copy helpers have no JAX body or callers and remain explicitly omitted. |
| [pdf_utilities.F90](../src/CLUBB_core/pdf_utilities.F90) | [pdf_utilities.hpp](src/CLUBB_core/pdf_utilities.hpp) | All 19 routines source-ported/audited; selected callers pass19 short CPU cases and four all-stats cases. Unused lognormal/double branches are not runtime claims. |
| [penta_bicgstab_solver.F90](../src/CLUBB_core/penta_bicgstab_solver.F90) | — | Explicitly excluded with JAX: the current C++ solver wrapper reports the unsupported method and sets fatal error status. No iterative solver port is claimed. |
| [penta_lu_solver.F90](../src/CLUBB_core/penta_lu_solver.F90) | [penta_lu_solver.cpp](src/CLUBB_core/penta_lu_solver.cpp) | Ported/audited; direct or case comparisons pass. |
| [pos_definite_module.F90](../src/CLUBB_core/pos_definite_module.F90) | [pos_definite_module.cpp](src/CLUBB_core/pos_definite_module.cpp) | Routine source-ported/audited; combined 19 short CPU cases, active-scalar all-stats and full 360-step BOMEX are byte-identical. |
| [precipitation_fraction.F90](../src/CLUBB_core/precipitation_fraction.F90) | [precipitation_fraction.cpp](src/CLUBB_core/precipitation_fraction.cpp) | All four routines source-ported/audited; isolated CPU group and active rainy/cloudy all-stats pass. Combined 19-case CPU/debug2 and four-case all-stats/debug2 checks pass. |
| [remapping_module.F90](../src/CLUBB_core/remapping_module.F90) | [remapping_module.cpp](src/CLUBB_core/remapping_module.cpp) | All 18 supported routines promoted and independently source-reviewed, including cached Ullrich operators, PPM boundary/limiter modes and diagnostics. Serial/CUDA compilation passes. C++ time-dependent dycore forcing now calls this owner; C++ stats output-grid remapping now uses this owner too. Core-phase tests did not exercise C++ remapping; the later 190-step active ATEX fixed-dycore comparison matches Native. Other remapping modes retain their explicit runtime limits. Source undefinedness is retained explicitly. |
| [saturation.F90](../src/CLUBB_core/saturation.F90) | [saturation.hpp](src/CLUBB_core/saturation.hpp) | All 10 routines source-ported/audited;19 short CPU cases and three additional liquid-selector BOMEX checks are byte-identical. Optional partial ranges and enabled specific humidity remain untested. |
| [setup_clubb_pdf_params.F90](../src/CLUBB_core/setup_clubb_pdf_params.F90) | [setup_clubb_pdf_params.cpp](src/CLUBB_core/setup_clubb_pdf_params.cpp) | All 18 routines source-ported/audited, reached through the C++ microphysics-preparation caller. BFB-018 records the verified emergency-variance contraction. Isolated 19-case CPU matrix and cloudy/SILHS all-stats pass; combined 19-case CPU/debug2 and four-case all-stats/debug2 checks pass. |
| [sfc_varnce_module.F90](../src/CLUBB_core/sfc_varnce_module.F90) | [sfc_varnce_module.cpp](src/CLUBB_core/sfc_varnce_module.cpp) | Full 33-argument parent independently source-audited; verified BFB-015 preserves rounded-quotient powers for heat, moisture and passive scalars. Combined19 short CPU cases, debug2/all-stats and full 360-step BOMEX are byte-identical; lower checker now C++. |
| [sigma_sqd_w_module.F90](../src/CLUBB_core/sigma_sqd_w_module.F90) | [sigma_sqd_w_module.cpp](src/CLUBB_core/sigma_sqd_w_module.cpp) | Ported/audited; direct or case comparisons pass. |
| [sponge_layer_damping.F90](../src/CLUBB_core/sponge_layer_damping.F90) | — | Sponge damping is rejected by the supported standalone driver. |
| [stats_clubb_utilities.F90](../src/CLUBB_core/stats_clubb_utilities.F90) | [stats_clubb_utilities.cpp](src/CLUBB_core/stats_clubb_utilities.cpp) | All three diagnostics routines ported; later stats milestone owns integrated acceptance. |
| [stats_netcdf.F90](../src/CLUBB_core/stats_netcdf.F90) | [stats_netcdf.cpp](src/CLUBB_core/stats_netcdf.cpp) | Complete lifecycle, registry, sampling, budgets, remapping and NetCDF writer port; later stats milestone owns integrated acceptance. |
| [tridiag_lu_solver.F90](../src/CLUBB_core/tridiag_lu_solver.F90) | [tridiag_lu_solver.cpp](src/CLUBB_core/tridiag_lu_solver.cpp) | Ported/audited; direct or case comparisons pass. |
| [turbulent_adv_pdf.F90](../src/CLUBB_core/turbulent_adv_pdf.F90) | [turbulent_adv_pdf.cpp](src/CLUBB_core/turbulent_adv_pdf.cpp) | All 4 routines source-ported/audited; combined 19 short CPU cases and full 360-step BOMEX are byte-identical. |

## Core-phase integration evidence (`ac37a3898`)

The results in this section predate the later driver-package ports. See the
[driver plan](DRIVER_CONVERSION_PLAN.md) for current integration checks and
the full driver-duration matrix and additional Morrison validation.

The four main moment parents and their compatible lower-physics group are live.
Their combined CPU comparison passes all 19 short native-timestep cases and the
full 360-step BOMEX case with standard stats. These Debug checks establish the
recorded numerical scope, without an optimized-performance claim.

Full numerical checks and the 15 type lifecycle routines pass the earlier
19-case CPU model gate with debug level 2 and four active all-stats cases, plus
existing UDT/fatal-return checks. The metadata/setup/precipitation group also
passes its combined matrix and cloudy/SILHS all-stats cases.

The PDF algorithms, parameter and flag initialization, grid construction, fixed
dycore-grid helper and remapping bodies are now promoted.
The isolated PDF and initialization groups pass their recorded short CPU cases;
PDF7 also passes active scalar/all-stats checks. The final combined CPU build
succeeds. Its five-step native-timestep matrix passes all 19 cases with
multicolumn stats and debug level 2. Default/PDF7 four-case
checks pass 8/8; COBRA all-stats predictor variants pass 3/3; fixed dycore-grid
BOMEX/ATEX checks pass 2/2 with five 30-second steps. All comparisons require file BFB. Existing CTest
checks pass 14/14 and focused Kokkos pytests pass 31/31.

The full 19-case native-duration/native-timestep matrix with standard stats
passes from a frozen executable pair: all saved files are byte-identical,
both models complete the independently expected 17,298 steps, and neither binary
changes during the run. The build uses NVHPC 26.9, double precision, Debug flags
and Serial Kokkos; this is not an optimized-performance claim. The existing
reset/reinitialization sequence also matches an ordinary two-step Native run,
including its observed step count. Frozen hashes, per-case counts and results
are under `build/validation/final_supported_core_cpu/`.

The complete CUDA executable builds, and
its four-case GPU checkpoint completes with no additional saved-output changes
from the preceding physics/metadata checkpoint. GPU versus CPU differences remain
recorded in the [BFB log](BFB_ADAPTATIONS.md#final-pdf-and-initialization-checkpoint).
An additional active Hybrid7 COBRA GPU checkpoint completes before and after the
integration. Its hybrid saved file changes while the default ADG1 file and both
Native reference files remain byte-identical; the
[record](BFB_ADAPTATIONS.md#active-hybrid7-gpu-checkpoint) identifies that boundary
without investigating GPU numerical causes.
The full matrix verifies the recorded case configurations; source-only and
untested branches remain distinguished in the inventory.

## Current Native dependencies

Statistics, debug levels, timers and input services now execute in C++.
The optional PDF matrix calls retain their external LAPACK ABI in
`KOKKOS_BFB=ON` mode. Native derived-type converters and numerical oracles are
test-only. The default production executable no longer links Native stats or
core scientific libraries; the local project still configures Fortran reference
targets. See the [statistics documentation](README.md#statistics) and the
[later driver milestones](DRIVER_CONVERSION_PLAN.md) for the current owners.

The core-phase acceptance evidence above remains valid for its recorded
binaries; it does not accept later driver or stats replacements. KK/Morrison
currently execute in host C++ with Kokkos View staging; GPU optimization is deferred.
