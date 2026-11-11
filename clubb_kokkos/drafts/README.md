# Port drafts and source reviews

These files hold source ports and review records before and after promotion.
CMake compiles the promoted copies in [src/CLUBB_core](../src/CLUBB_core/), not
this directory. Runtime evidence applies to the recorded live copies and options;
source review or a changed draft alone does not establish numerical agreement.
The supported source inventory is promoted, and the combined short CPU
comparisons pass strictly bit-for-bit (BFB). All 19 native-duration/native-timestep
CPU cases with standard stats also produce byte-identical saved files. The CUDA
executable builds and its default and active-hybrid GPU checkpoints complete; see the
[core inventory](../CORE_CONVERSION_PLAN.md) for accepted scope and retained
Native dependencies.

| Draft directory | Canonical source | Current work |
| --- | --- | --- |
| [stats_netcdf](stats_netcdf/) | [stats_netcdf.F90](../../src/CLUBB_core/stats_netcdf.F90) | All 40 Native routines promoted to the permanent C++ stats owner; source, cleanup and independent reviews retained here. See the later driver milestone for model acceptance. |
| [stats_clubb_utilities](stats_clubb_utilities/) | [stats_clubb_utilities.F90](../../src/CLUBB_core/stats_clubb_utilities.F90) | All three diagnostic accumulation routines promoted and independently reviewed. Host accumulation is the initial implementation; GPU tuning is deferred. |
| [adg1_adg2_3d_luhar_pdf](adg1_adg2_3d_luhar_pdf/) | [adg1_adg2_3d_luhar_pdf.F90](../../src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90) | All nine routines promoted and independently source-reviewed, including device Luhar/root adaptations. ADG1 passes the 19-case short CPU matrix; ADG2/Luhar selectors 2/3 remain input-fields-gated. |
| [mixing_length](mixing_length/) | [mixing_length.F90](../../src/CLUBB_core/mixing_length.F90) | Seven routines independently source-audited and promoted; all 19 short CPU cases pass. Methods 1 and 2 also pass two targeted cases each. |
| [pdf_utilities](pdf_utilities/) | [pdf_utilities.F90](../../src/CLUBB_core/pdf_utilities.F90) | All 19 routines promoted/source-audited;19 short CPU cases and four all-stats cases pass. GPU builds/runs; BOMEX/COBRA exact, DYCOMS/JUN25 small recorded differences. |
| [saturation](saturation/) | [saturation.F90](../../src/CLUBB_core/saturation.F90) | Ten routines promoted/source-audited;19 short CPU cases and three liquid-selector checks pass. The later GPU checkpoint records changed outputs for the combined moment/saturation/updraft group. |
| [advance_helper_module](advance_helper_module/) | [advance_helper_module.F90](../../src/CLUBB_core/advance_helper_module.F90) | All 26 procedures promoted and independently source-reviewed. Recorded stability/shear/splat/radiation and moment-parent CPU groups pass; other generic/reduction paths retain their runtime limits. |
| [calc_roots](calc_roots/) | [calc_roots.F90](../../src/CLUBB_core/calc_roots.F90) | All three Native APIs and two device scalar overloads promoted and independently source-reviewed. CPU/CUDA compilation passes; numerical evidence is limited to actual PDF callers. |
| [new_hybrid_pdf](new_hybrid_pdf/) | [new_hybrid_pdf_main.F90](../../src/CLUBB_core/new_hybrid_pdf_main.F90), [new_hybrid_pdf.F90](../../src/CLUBB_core/new_hybrid_pdf.F90) | All ten routines across both files promoted/source-audited; required constants and roots are integrated. Isolated default/PDF7 multicolumn and active scalar/all-stats checks pass. BFB-019/020/021 record verified arithmetic choices; combined short CPU checks pass byte-identically. |
| [new_pdf](new_pdf/) | [new_pdf_main.F90](../../src/CLUBB_core/new_pdf_main.F90), [new_pdf.F90](../../src/CLUBB_core/new_pdf.F90) | All 13 routines promoted and independently source-reviewed, including complex-root limits and all advection flags. Final comment/group cleanup preserves executable tokens. Selector 4 remains input-fields-gated; runtime/BFB is unverified. |
| [sfc_varnce_module](sfc_varnce_module/) | [sfc_varnce_module.F90](../../src/CLUBB_core/sfc_varnce_module.F90) | Complete independently reviewed parent promoted; verified BFB-015 covers rounded-quotient powers. Combined CPU cases, active all-stats/scalar checks and full BOMEX pass. |
| [advance_xm_wpxp_module](advance_xm_wpxp_module/) | [advance_xm_wpxp_module.F90](../../src/CLUBB_core/advance_xm_wpxp_module.F90) | All 18 routines source-audited and promoted; BFB-017 logs the CPU contraction adaptation. Combined 19 short CPU cases and full BOMEX pass. |
| [advance_wp2_wp3_module](advance_wp2_wp3_module/) | [advance_wp2_wp3_module.F90](../../src/CLUBB_core/advance_wp2_wp3_module.F90) | All 23 routines independently source-reviewed and promoted; combined 19 short CPU cases and full BOMEX pass. |
| [advance_xp2_xpyp_module](advance_xp2_xpyp_module/) | [advance_xp2_xpyp_module.F90](../../src/CLUBB_core/advance_xp2_xpyp_module.F90) | All 16 routines source-audited and promoted; BFB-016 logs the CPU contraction adaptation. Combined 19 short CPU cases and full BOMEX pass; multiple-LHS paths remain untested. |
| [advance_xp3_module](advance_xp3_module/) | [advance_xp3_module.F90](../../src/CLUBB_core/advance_xp3_module.F90) | All five routines independently source-reviewed and promoted; 19 short CPU cases and an active PDF7 all-stats check pass. Native compile-time-disabled advancement remains unexercised. |
| [advance_windm_edsclrm_module](advance_windm_edsclrm_module/) | [advance_windm_edsclrm_module.F90](../../src/CLUBB_core/advance_windm_edsclrm_module.F90) | All six routines source-audited and promoted; combined 19 short CPU cases/full BOMEX and targeted predictor=false scalar checks pass. |
| [pdf_closure_module](pdf_closure_module/) | [pdf_closure_module.F90](../../src/CLUBB_core/pdf_closure_module.F90) | All 18 supported routines promoted; three uncalled private interpolation helpers are excluded. Earlier CPU multicolumn/all-stats/full BOMEX checks pass. Recorded PDF7 cases exercise fusion; default ADG1 skips it. Combined short CPU checks pass byte-identically. |
| [clip_explicit](clip_explicit/) | [clip_explicit.F90](../../src/CLUBB_core/clip_explicit.F90) | All six routines/19 constants source-audited and promoted; combined 19 short CPU cases and full BOMEX pass. |
| [diffusion](diffusion/) | [diffusion.F90](../../src/CLUBB_core/diffusion.F90) | Both routines source-audited and promoted; combined 19 short CPU cases and full BOMEX pass. |
| [mean_adv](mean_adv/) | [mean_adv.F90](../../src/CLUBB_core/mean_adv.F90) | Both routines source-audited and promoted; combined 19 short CPU cases and full BOMEX pass. |
| [LY93_pdf](LY93_pdf/) | [LY93_pdf.F90](../../src/CLUBB_core/LY93_pdf.F90) | All three routines promoted, compiled and independently source-reviewed. Selector 6 remains input-fields-gated; algorithm runtime/BFB is unverified. |
| [new_tsdadg_pdf](new_tsdadg_pdf/) | [new_tsdadg_pdf.F90](../../src/CLUBB_core/new_tsdadg_pdf.F90) | All four routines promoted, compiled and independently source-reviewed. Selector 5 remains input-fields-gated; algorithm runtime/BFB is unverified. |
| [numerical_check](numerical_check/) | [numerical_check.F90](../../src/CLUBB_core/numerical_check.F90) | All 12 routines promoted/source-audited; earlier combined 19 short CPU cases and four active all-stats cases pass with debug level 2. Parameter validation is now C++; shared debug/header-format metadata remains Native. |
| [pdf_closure_module_updraft](pdf_closure_module_updraft/) | [pdf_closure_module.F90](../../src/CLUBB_core/pdf_closure_module.F90) | Cloudy-updraft helper promoted, with active nonzero CPU all-stats outputs exact. All18 source-used PDF routines are live; three uncalled private interpolation helpers are excluded with source-call-graph evidence. |
| [calc_pressure](calc_pressure/) | [calc_pressure.F90](../../src/CLUBB_core/calc_pressure.F90) | Both routines source-audited and promoted; three calculate_thvm callers pass 19 short CPU cases/full BOMEX. init_pressure has no newly ported caller. |
| [turbulent_adv_pdf](turbulent_adv_pdf/) | [turbulent_adv_pdf.F90](../../src/CLUBB_core/turbulent_adv_pdf.F90) | All four routines independently source-reviewed and promoted; combined 19 short CPU cases and full BOMEX pass. |
| [pos_definite_module](pos_definite_module/) | [pos_definite_module.F90](../../src/CLUBB_core/pos_definite_module.F90) | Routine independently source-reviewed and promoted; combined 19 short CPU cases/full BOMEX and active scalar all-stats checks pass. |
| [fill_holes](fill_holes/) | [fill_holes.F90](../../src/CLUBB_core/fill_holes.F90) | Ten supported routines source-audited and promoted; four JAX-omitted methods remain explicitly excluded. Combined CPU cases/full BOMEX pass; not every strategy is exercised. |
| [mono_flux_limiter](mono_flux_limiter/) | [mono_flux_limiter.F90](../../src/CLUBB_core/mono_flux_limiter.F90) | All seven routines source-audited and promoted; combined 19 short CPU cases and full BOMEX pass. |
| [index_mapping](index_mapping/) | [index_mapping.F90](../../src/CLUBB_core/index_mapping.F90) | All five routines source-audited and promoted; combined 19 short CPU cases and full BOMEX pass. |
| [corr_varnce_module](corr_varnce_module/) | [corr_varnce_module.F90](../../src/CLUBB_core/corr_varnce_module.F90) | All seven routines and immutable tables source-audited and promoted. Isolated metadata/setup CPU group passes; combined 19 CPU/debug2 cases and four all-stats/debug2 cases pass. |
| [matrix_operations](matrix_operations/) | [matrix_operations.F90](../../src/CLUBB_core/matrix_operations.F90) | Four approved routines and required LAPACK interfaces source-audited and promoted. Isolated group passes; combined 19 CPU/debug2 cases and four all-stats/debug2 cases pass. |
| [diagnose_correlations_module](diagnose_correlations_module/) | [diagnose_correlations_module.F90](../../src/CLUBB_core/diagnose_correlations_module.F90) | All ten active routines source-audited and promoted; unsupported l_calc_w_corr=true remains gated. Normal isolated group passes; two diagnose probes fail the same symmetry assertion in both models. Combined 19 CPU/debug2 cases and four all-stats/debug2 cases pass. |
| [err_info_type_module](err_info_type_module/) | [err_info_type_module.F90](../../src/CLUBB_core/err_info_type_module.F90) | All five lifecycle routines source-audited and promoted; exact Native formatting remains an explicit transient presence-overload call. Combined 19 CPU cases/four all-stats cases and UDT/fatal-return checks pass. |
| [hydromet_pdf_parameter_module](hydromet_pdf_parameter_module/) | [hydromet_pdf_parameter_module.F90](../../src/CLUBB_core/hydromet_pdf_parameter_module.F90) | All three lifecycle routines source-audited and promoted with the generated collection owner; combined 19 CPU cases/four all-stats cases and UDT checks pass. |
| [pdf_parameter_module](pdf_parameter_module/) | [pdf_parameter_module.F90](../../src/CLUBB_core/pdf_parameter_module.F90) | All seven JAX-supported lifecycle/pack/print routines source-audited and promoted; combined 19 CPU cases/four all-stats cases pass. Uncalled Native-only copy helpers are explicitly omitted. |
| [calendar](calendar/) | [calendar.F90](../../src/CLUBB_core/calendar.F90) | All five routines source-audited and promoted; CPU compilation and tiny rollover/century/error checks pass. Current model date callers remain Native. |
| [Nc_Ncn_eqns](Nc_Ncn_eqns/) | [Nc_Ncn_eqns.F90](../../src/CLUBB_core/Nc_Ncn_eqns.F90) | All six routines source-audited and promoted; isolated CPU group/cloudy all-stats pass, combined 19 CPU/debug2 cases and four all-stats/debug2 cases pass. |
| [precipitation_fraction](precipitation_fraction/) | [precipitation_fraction.F90](../../src/CLUBB_core/precipitation_fraction.F90) | All four routines source-audited and promoted; isolated CPU group/cloudy all-stats pass, combined 19 CPU/debug2 cases and four all-stats/debug2 cases pass. |
| [setup_clubb_pdf_params](setup_clubb_pdf_params/) | [setup_clubb_pdf_params.F90](../../src/CLUBB_core/setup_clubb_pdf_params.F90) | All 18 routines source-audited and promoted with the bounded outer wrapper. Verified BFB-018 and isolated CPU group/cloudy/SILHS all-stats pass; combined 19 CPU/debug2 cases and four all-stats/debug2 cases pass. |
| [grid_class](grid_class/) | [grid_class.F90](../../src/CLUBB_core/grid_class.F90) | Five remaining construction/height/cleanup/weight bodies promoted; existing 17 operators preserved. Grid-file reading remains Native I/O. Isolated grid types 1/2/3 pass; combined short CPU checks pass byte-identically. |
| [model_flags](model_flags/) | [model_flags.F90](../../src/CLUBB_core/model_flags.F90) | All three default/init/print procedures promoted around the generated schema. Case-owned flags and C++ stream output are permanent; integer-unit printing retains Native I/O. Isolated cases pass; combined short CPU checks pass byte-identically. |
| [parameters_tunable](parameters_tunable/) | [parameters_tunable.F90](../../src/CLUBB_core/parameters_tunable.F90) | All eight supported bodies promoted; Native namelist file I/O only, with required iunit/NEWUNIT adaptation documented. Two tuner readers remain excluded. Isolated cases pass; combined short CPU checks pass byte-identically. |
| [grid_adaptation_module](grid_adaptation_module/) | [grid_adaptation_module.F90](../../src/CLUBB_core/grid_adaptation_module.F90) | Bounded fixed setup_gr_dycore promoted and independently source-reviewed. Native grid-file reading remains; broader adaptive-grid/generalization algorithms remain excluded. |
| [remapping_module](remapping_module/) | [remapping_module.F90](../../src/CLUBB_core/remapping_module.F90) | All 18 supported routines promoted and independently source-reviewed. Serial/CUDA compilation passes; current Native outer/stat callers remain unchanged, so C++ remapping runtime/BFB is unverified. |

Each owner keeps its source outline, C++ files, source-audit record and
integration notes together. Notes identify required constants and lower calls,
unsupported or untested branches, and any deliberate BFB adaptations. Owners
edit only their draft directory; integration owns shared generators, CMake and
live call sites.

Promote compatible reviewed files together into `src/CLUBB_core/`, wire their
explicit lower boundaries and remove replaced whole-routine bridges. Record the
promoted group, build it and run the maintained Kokkos/Fortran comparison. A CPU
BFB group needs no separate per-file model checks; narrow the added files only
when the group fails. Use short `multi_col_stats` runs for the usual
cycle and broader stats/duration at milestones or for a specific diagnosis.
Draft work on other files can continue throughout that cycle. Independent CPU
integration/test lanes use frozen source copies and separate build/output trees.
After reviewing and combining passing deltas, recheck the combined model. Data-only
benchmark symlinks can be shared; source/runtime scripts must resolve within each lane.

Keep the [core inventory](../CORE_CONVERSION_PLAN.md) explicit about which
portions are integrated and verified. Record forced arithmetic choices in the
[BFB adaptation log](../BFB_ADAPTATIONS.md) when promoting them, including pending
verification. A draft's source audit does not establish BFB or GPU performance.


## Standalone driver groups

The supported driver groups are integrated in [src](../src/); their current
scope and numerical acceptance are recorded in the [driver plan](../DRIVER_CONVERSION_PLAN.md).
The JSON outlines and audits below describe the reviewed draft snapshots. Their
pending-build fields and historical hashes refer to those snapshots, rather than
the current implementation. Promoted draft source copies, one-off translation
scripts and obsolete integration instructions have been removed. Edit the
permanent source owners for further changes.

| Review records | Permanent source |
| --- | --- |
| [Benchmark outline](driver_benchmark_cases/SOURCE_OUTLINE.json) | [Benchmark_cases](../src/Benchmark_cases/) |
| [Input inventory](driver_input_fields/PACKAGE_OUTLINE.json) and [IOP outline](driver_input_fields/readiopdata_SOURCE_OUTLINE.json) | [Input_fields](../src/Input_fields/), [IOP wrapper](../src/readiopdata_module.cpp) and [IOP core](../src/readiopdata_core/) |
| [Radiation outline](driver_radiation_silhs/Radiation/SOURCE_OUTLINE.json) and [SILHS outline](driver_radiation_silhs/SILHS/SOURCE_OUTLINE.json) | [Radiation](../src/Radiation/) and [SILHS](../src/SILHS/) |
| [Microphysics outline](driver_microphysics/source_outline.json), [KK outline](driver_microphysics/KK_engine_outline.json) and [Morrison outline](driver_microphysics/Morrison_engine_outline.json) | [Microphys](../src/Microphys/), including both engine subdirectories |
