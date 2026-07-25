# E3SM Maint-3.2 Tunable Defaults Report

This report documents why each value in this directory was chosen. The goal is to make `input/parameter_and_flag_configs/e3sm_maint32/` auditable against E3SM maint-3.2, not to claim that every current CLUBB option existed in that E3SM branch. When a value comes from E3SM, the table names the source type: E3SM XML namelist default, E3SM CLUBB source constant, EAM wrapper behavior, or an older parameter name. When no maint-3.2 equivalent was found, the table says whether the current CLUBB default was retained as infrastructure or whether a current-only feature was shut off.

Primary sources checked:

- `/home/guntherhuebler/Code/E3SM/components/eam/bld/namelist_files/namelist_defaults_eam.xml`
- `/home/guntherhuebler/Code/E3SM/components/eam/bld/namelist_files/namelist_definition.xml`
- `/home/guntherhuebler/Code/E3SM/components/eam/src/physics/clubb/parameters_tunable.F90`
- `/home/guntherhuebler/Code/E3SM/components/eam/src/physics/clubb/model_flags.F90`
- `/home/guntherhuebler/Code/E3SM/components/eam/src/physics/cam/clubb_intr.F90`
- E3SM CLUBB source files where constants later became tunable in current CLUBB, including `advance_helper_module.F90`, `advance_clubb_core_module.F90`, `advance_xp2_xpyp_module.F90`, and `surface_varnce_module.F90`.

## `tunable_parameters.in`

| Name | Value | Rationale |
|---|---:|---|
| `C1` | `2.400000` | Found as `clubb_C1` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C1b` | `2.800000` | Found as `clubb_C1b` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C1c` | `0.750000` | Found as `clubb_C1c` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C2rt` | `1.750000` | Found as `clubb_C2rt` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C2thl` | `1.750000` | E3SM wrapper sets `C2thl = C2rt` when `clubb_C2rt` is specified and `clubb_C2thl` is not. |
| `C2rtthl` | `2.275000` | E3SM wrapper sets `C2rtthl = 1.3*C2rt` when `clubb_C2rt` is specified and `clubb_C2rtthl` is not. |
| `C4` | `5.200000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C_uu_shr` | `0.300000` | Current split of old E3SM pressure coefficient `C5`; E3SM `parameters_tunable.F90` default `C5 = 0.3`. |
| `C_uu_buoy` | `0.300000` | Current split of old E3SM pressure coefficient `C5`; E3SM `parameters_tunable.F90` default `C5 = 0.3`. |
| `C6rt` | `4.000000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C6rtb` | `7.500000` | Found as `clubb_C6rtb` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C6rtc` | `0.500000` | Found as `clubb_C6rtc` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C6thl` | `4.000000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C6thlb` | `7.500000` | Found as `clubb_C6thlb` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C6thlc` | `0.500000` | Found as `clubb_C6thlc` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C7` | `0.500000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C7b` | `0.500000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C7c` | `0.500000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C8` | `5.200000` | Found as generic `clubb_C8` with `phys="default"` in E3SM `namelist_defaults_eam.xml`; hgrid-specific `ne120np4` override was not used for the generic file. |
| `C8b` | `0.000000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C11` | `0.700000` | Found as `clubb_C11` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C11b` | `0.200000` | Found as `clubb_C11b` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C11c` | `0.850000` | Found as `clubb_C11c` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C12` | `1.000000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C14` | `2.500000` | Found as `clubb_C14` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `C_wp2_pr_dfsn` | `0.000000` | Current-only pressure term not found in E3SM maint-3.2; shut off. |
| `C_wp3_pr_tp` | `0.000000` | Current-only pressure term not found in E3SM maint-3.2; shut off. |
| `C_wp3_pr_turb` | `0.000000` | Current-only pressure term not found in E3SM maint-3.2; current source default is disabled. |
| `C_wp3_pr_dfsn` | `0.000000` | Current-only pressure term not found in E3SM maint-3.2; shut off. |
| `C_wp2_splat` | `0.00` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C6rt_Lscale0` | `14.0` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C6thl_Lscale0` | `14.0` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C7_Lscale0` | `0.850000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `wpxp_L_thresh` | `100.0` | Found as `clubb_wpxp_L_thresh` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `c_K` | `0.200000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `c_K1` | `0.750000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `nu1` | `20.00000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `c_K2` | `0.125000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `nu2` | `5.000000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `c_K6` | `0.375000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `nu6` | `5.000000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `c_K8` | `1.250000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `nu8` | `20.00000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `c_K9` | `0.250000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `nu9` | `20.00000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `nu10` | `0.000000` | Found as a source default constant in E3SM `parameters_tunable.F90`; this keeps eddy diffusion for extra scalar and wind fields disabled. |
| `c_K_hm` | `0.750000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `c_K_hmb` | `0.100000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `K_hm_min_coef` | `0.10000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `nu_hm` | `1.500000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `slope_coef_spread_DG_means_w` | `21.00000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `pdf_component_stdev_factor_w` | `6.500000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `coef_spread_DG_means_rt` | `0.800000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `coef_spread_DG_means_thl` | `0.800000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `beta` | `2.400000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `gamma_coef` | `0.120000` | Found as `clubb_gamma_coef` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `gamma_coefb` | `0.280000` | Found as `clubb_gamma_coefb` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `gamma_coefc` | `1.200000` | Found as `clubb_gamma_coefc` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `taumax` | `3600.000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `lmin_coef` | `0.100000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `omicron` | `0.800000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `zeta_vrnce_rat` | `0.000000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `upsilon_precip_frac_rat` | `0.900000` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `lambda0_stability_coef` | `0.03` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `mult_coef` | `1.000000` | Found as a source default constant in E3SM `parameters_tunable.F90` for the non-CLUBBND CAM path. |
| `mu` | `5.000E-4` | Found as `clubb_mu` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `Lscale_mu_coef` | `2.0` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `Lscale_pert_coef` | `0.1` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `alpha_corr` | `0.15` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `Skw_denom_coef` | `0.0` | Found as a conditional E3SM source default in `parameters_tunable.F90` for the CAM build path. |
| `c_K10` | `0.35` | Found as `clubb_c_K10` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `c_K10h` | `0.35` | Found as `clubb_c_K10h` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `thlp2_rad_coef` | `1.0` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `thlp2_rad_cloud_frac_thresh` | `0.1` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `up2_sfc_coef` | `2.0` | Found as old E3SM source parameter `up2_vp2_factor = 2.0`; current `up2_sfc_coef` is the same surface `up2` and `vp2` multiplier in the current formula. |
| `Skw_max_mag` | `4.5` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C_invrs_tau_bkgnd` | `1.0` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C_invrs_tau_sfc` | `0.1` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C_invrs_tau_shear` | `0.02` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C_invrs_tau_N2` | `0.1` | Found as a source default constant in E3SM `parameters_tunable.F90`. |
| `C_invrs_tau_N2_wp2` | `0.1` | Current split of old single E3SM `C_invrs_tau_N2`; set to the old E3SM coefficient because maint-3.2 did not have per-equation N2 coefficients. |
| `C_invrs_tau_N2_xp2` | `0.1` | Current split of old single E3SM `C_invrs_tau_N2`; set to the old E3SM coefficient because maint-3.2 did not have per-equation N2 coefficients. |
| `C_invrs_tau_N2_wpxp` | `0.1` | Current split of old single E3SM `C_invrs_tau_N2`; set to the old E3SM coefficient because maint-3.2 did not have per-equation N2 coefficients. |
| `C_invrs_tau_N2_clear_wp3` | `0.1` | Current split of old single E3SM `C_invrs_tau_N2`; set to the old E3SM coefficient because maint-3.2 did not have per-equation N2 coefficients. |
| `C_invrs_tau_wpxp_Ri` | `0.0` | Current-only Ri multiplier for `invrs_tau_wpxp`; not found in E3SM maint-3.2, so it is shut off. |
| `C_invrs_tau_wpxp_N2_thresh` | `3.3E-4` | Current-only threshold for the Ri multiplier path; retained current default because `C_invrs_tau_wpxp_Ri = 0.0` disables the multiplier. |
| `xp3_coef_base` | `0.25` | Current-only simple `xp3` tunability; no E3SM maint-3.2 equivalent found, so current source default was retained. |
| `xp3_coef_slope` | `0.01` | Current-only simple `xp3` tunability; no E3SM maint-3.2 equivalent found, so current source default was retained. |
| `altitude_threshold` | `100.0` | Found as a source default constant in E3SM `parameters_tunable.F90`; also exposed in E3SM as `clubb_altitude_thresh`. |
| `rtp2_clip_coef` | `0.5` | Found as a hard-coded constant in E3SM `advance_xp2_xpyp_module.F90`; now tunable in current CLUBB. |
| `Richardson_num_min` | `0.25` | Found as hard-coded `Richardson_num_min = one_fourth` in E3SM `advance_helper_module.F90`; now tunable in current CLUBB. |
| `Richardson_num_max` | `400.0` | Found as hard-coded `Richardson_num_max = 400` in E3SM `advance_helper_module.F90`; now tunable in current CLUBB. |
| `Cx_min` | `0.333333` | Found as hard-coded `Cx_min = one_third` in E3SM `advance_helper_module.F90`; now tunable in current CLUBB. |
| `Cx_max` | `0.95` | Found as hard-coded `Cx_max = 0.95` in E3SM `advance_helper_module.F90`; now tunable in current CLUBB. |
| `a3_coef_min` | `1.0` | Current tunable value matches current CLUBB source default; no separate E3SM maint-3.2 tunable was found. |
| `a_const` | `1.8` | Found as a hard-coded source constant in E3SM `surface_varnce_module.F90`; now tunable in current CLUBB. |
| `bv_efold` | `5.0` | Current-only mixed Brunt-Vaisala cloud-fraction control; no E3SM maint-3.2 equivalent found, so current source default was retained. |
| `wpxp_Ri_exp` | `.5` | Current-only exponent for the Ri multiplier path; retained current default because `C_invrs_tau_wpxp_Ri = 0.0` disables the multiplier. |
| `z_displace` | `20.0` | Found as hard-coded `z_displace = 20.0` in E3SM `advance_clubb_core_module.F90`; now tunable in current CLUBB. |

## `configurable_model_flags.in`: CLUBB flags and integer options

| Name | Value | Rationale |
|---|---:|---|
| `iiPDF_type` | `1` | Current CLUBB PDF-selection infrastructure; no E3SM maint-3.2 namelist equivalent found, so current standalone default was retained. |
| `ipdf_call_placement` | `2` | Found as `clubb_ipdf_call_placement` with `phys="default"` in E3SM `namelist_defaults_eam.xml`. |
| `penta_solve_method` | `2` | Current solver-selection infrastructure; no E3SM maint-3.2 namelist equivalent found, so current standalone default was retained. |
| `tridiag_solve_method` | `2` | Current solver-selection infrastructure; no E3SM maint-3.2 namelist equivalent found, so current standalone default was retained. |
| `saturation_formula` | `3` | E3SM `model_flags.F90` maps `flatau` to `saturation_flatau`, which is integer value `3`. |
| `grid_remap_method` | `1` | Current grid-remap infrastructure; no E3SM maint-3.2 namelist equivalent found, so current standalone default was retained. |
| `grid_adapt_in_time_method` | `0` | Current grid-adaptation infrastructure; no E3SM maint-3.2 namelist equivalent found, so current standalone default was retained. |
| `fill_holes_type` | `2` | Current hole-fill infrastructure; E3SM maint-3.2 only had older fixed hole-fill behavior, so current standalone default was retained. |
| `l_godunov_upwind_wpxp_ta` | `.false.` | Current-only Godunov variant not found in E3SM maint-3.2; shut off. |
| `l_upwind_xpyp_ta` | `.true.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_godunov_upwind_xpyp_ta` | `.false.` | Current-only Godunov variant not found in E3SM maint-3.2; shut off. |
| `l_upwind_xm_ma` | `.true.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_vert_avg_closure` | `.true.` | Found as `clubb_vert_avg_closure` default `.true.` in E3SM XML and as a source default in E3SM `model_flags.F90`. |
| `l_standard_term_ta` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_partial_upwind_wp3` | `.false.` | Current-only option not found in E3SM maint-3.2; shut off. |
| `l_use_cloud_cover` | `.true.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_rcm_supersat_adj` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_damp_wp3_Skw_squared` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_min_wp2_from_corr_wx` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_min_xp2_from_corr_wx` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_C2_cloud_frac` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_predict_upwp_vpwp` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_ho_nontrad_coriolis` | `.false.` | Current-only Coriolis option not found in E3SM maint-3.2; shut off. |
| `l_ho_trad_coriolis` | `.false.` | Current-only Coriolis option not found in E3SM maint-3.2; shut off. |
| `l_diag_Lscale_from_tau` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_stability_correct_tau_zm` | `.true.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_damp_wp2_using_em` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_use_C7_Richardson` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_use_precip_frac` | `.true.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_do_expldiff_rtm_thlm` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_use_C11_Richardson` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_use_shear_Richardson` | `.false.` | Current-only Richardson option not found in E3SM maint-3.2; shut off. |
| `l_prescribed_avg_deltaz` | `.false.` | Found as an E3SM non-GFDL source default in `parameters_tunable.F90`. |
| `l_diffuse_rtm_and_thlm` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`; E3SM wrapper only enables it when `clubb_stabcorrect` is true. |
| `l_stability_correct_Kh_N2_zm` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`; E3SM wrapper only enables it when `clubb_stabcorrect` is true. |
| `l_trapezoidal_rule_zt` | `.true.` | E3SM wrapper derives this from `clubb_vert_avg_closure = .true.`. |
| `l_trapezoidal_rule_zm` | `.true.` | E3SM wrapper derives this from `clubb_vert_avg_closure = .true.`. |
| `l_call_pdf_closure_twice` | `.true.` | E3SM wrapper derives this from `clubb_vert_avg_closure = .true.`. |
| `l_Lscale_plume_centered` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_brunt_vaisala_freq_moist` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_use_thvm_in_bv_freq` | `.false.` | Found as a source default constant in E3SM `model_flags.F90`. |
| `l_lmm_stepping` | `.false.` | Current-only time-stepping option not found in E3SM maint-3.2; shut off. |
| `l_e3sm_config` | `.false.` | Current-only compatibility switch not present in E3SM maint-3.2; shut off rather than assuming new behavior. |
| `l_vary_convect_depth` | `.false.` | Current-only option not found in E3SM maint-3.2; shut off. |
| `l_use_tke_in_wp3_pr_turb_term` | `.false.` | Current-only pressure-term option not found in E3SM maint-3.2; shut off. |
| `l_use_tke_in_wp2_wp3_K_dfsn` | `.false.` | Current-only diffusion option not found in E3SM maint-3.2; shut off. |
| `l_use_wp3_lim_with_smth_Heaviside` | `.false.` | Current-only limiter option not found in E3SM maint-3.2; shut off. |
| `l_smooth_Heaviside_tau_wpxp` | `.false.` | Current-only smooth Heaviside option not found in E3SM maint-3.2; shut off. |
| `l_modify_limiters_for_cnvg_test` | `.false.` | Current-only convergence-test option not found in E3SM maint-3.2; shut off. |
| `l_enable_relaxed_clipping` | `.false.` | Current-only clipping option not found in E3SM maint-3.2; shut off. |
| `l_linearize_pbl_winds` | `.false.` | Current-only PBL-wind option not found in E3SM maint-3.2; shut off. |
| `l_mono_flux_lim_thlm` | `.false.` | Current-only monotonic flux limiter not found in E3SM maint-3.2; shut off. |
| `l_mono_flux_lim_rtm` | `.false.` | Current-only monotonic flux limiter not found in E3SM maint-3.2; shut off. |
| `l_mono_flux_lim_um` | `.false.` | Current-only monotonic flux limiter not found in E3SM maint-3.2; shut off. |
| `l_mono_flux_lim_vm` | `.false.` | Current-only monotonic flux limiter not found in E3SM maint-3.2; shut off. |
| `l_mono_flux_lim_spikefix` | `.false.` | Current-only spike-fix option for monotonic flux limiting not found in E3SM maint-3.2; shut off. |
| `l_host_applies_sfc_fluxes` | `.false.` | Found as an EAM wrapper/source default in E3SM `clubb_intr.F90` and `model_flags.F90`. |
| `l_wp2_fill_holes_tke` | `.false.` | Current-only hole-fill option not found in E3SM maint-3.2; shut off. |
| `l_add_dycore_grid` | `.false.` | Current-only dycore-grid option not found in E3SM maint-3.2; shut off. |

## `configurable_model_flags.in`: SILHS flags and integer options

| Name | Value | Rationale |
|---|---:|---|
| `cluster_allocation_strategy` | `3` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `l_lh_importance_sampling` | `.true.` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `l_Lscale_vert_avg` | `.false.` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `l_lh_straight_mc` | `.false.` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `l_lh_clustered_sampling` | `.true.` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `l_rcm_in_cloud_k_lh_start` | `.true.` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `l_random_k_lh_start` | `.false.` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `l_max_overlap_in_cloud` | `.true.` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `l_lh_instant_var_covar_src` | `.true.` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `l_lh_limit_weights` | `.true.` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `l_lh_var_frac` | `.false.` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `l_lh_normalize_weights` | `.true.` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |

## `silhs_parameters.in`

| Name | Value | Rationale |
|---|---:|---|
| `importance_prob_thresh` | `1.0e-8` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
| `vert_decorr_coef` | `0.03` | No direct E3SM maint-3.2 CLUBB/SILHS namelist equivalent found; retained current standalone SILHS default. |
