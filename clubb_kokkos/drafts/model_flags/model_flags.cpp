// Set, assemble and print CLUBB's case-owned configuration flags.
//
// set_default_clubb_config_flags_api writes Native defaults into caller scalars;
// initialize_clubb_config_flags_type_api validates and assembles the existing
// configuration type. print_clubb_config_flags_api preserves Native row/value
// order for a C++ output stream, or retains public Native I/O for an integer unit.
// No mutable module-level flag state is owned here. Native namelist overrides
// follow default initialization at the existing driver I/O boundary.
#include "model_flags.hpp"
#include "f2k/procedure_bridge.hpp"

#include <ostream>
#include <sstream>
#include <stdexcept>

namespace clubb_kokkos {

// =============================================================================
void set_default_clubb_config_flags_api(int& iiPDF_type,                         // Out
                                        int& ipdf_call_placement,                // Out
                                        int& penta_solve_method,                 // Out
                                        int& tridiag_solve_method,               // Out
                                        int& saturation_formula,                 // Out
                                        int& Lscale_method,                      // Out
                                        int& grid_remap_method,                  // Out
                                        int& grid_adapt_in_time_method,          // Out
                                        int& fill_holes_type,                    // Out
                                        bool& l_use_precip_frac,                 // Out
                                        bool& l_predict_upwp_vpwp,               // Out
                                        bool& l_ho_nontrad_coriolis,             // Out
                                        bool& l_ho_trad_coriolis,                // Out
                                        bool& l_min_wp2_from_corr_wx,            // Out
                                        bool& l_min_xp2_from_corr_wx,            // Out
                                        bool& l_C2_cloud_frac,                   // Out
                                        bool& l_diffuse_rtm_and_thlm,            // Out
                                        bool& l_stability_correct_Kh_N2_zm,      // Out
                                        bool& l_calc_thlp2_rad,                  // Out
                                        bool& l_upwind_xpyp_ta,                  // Out
                                        bool& l_upwind_xm_ma,                    // Out
                                        bool& l_uv_nudge,                        // Out
                                        bool& l_rtm_nudge,                       // Out
                                        bool& l_vert_avg_closure,                // Out
                                        bool& l_trapezoidal_rule_zt,             // Out
                                        bool& l_trapezoidal_rule_zm,             // Out
                                        bool& l_call_pdf_closure_twice,          // Out
                                        bool& l_standard_term_ta,                // Out
                                        bool& l_partial_upwind_wp3,              // Out
                                        bool& l_godunov_upwind_wpxp_ta,          // Out
                                        bool& l_godunov_upwind_xpyp_ta,          // Out
                                        bool& l_use_cloud_cover,                 // Out
                                        bool& l_diagnose_correlations,           // Out
                                        bool& l_calc_w_corr,                     // Out
                                        bool& l_const_Nc_in_cloud,               // Out
                                        bool& l_fix_w_chi_eta_correlations,      // Out
                                        bool& l_stability_correct_tau_zm,        // Out
                                        bool& l_damp_wp2_using_em,               // Out
                                        bool& l_do_expldiff_rtm_thlm,            // Out
                                        bool& l_Lscale_plume_centered,           // Out
                                        bool& l_use_C7_Richardson,               // Out
                                        bool& l_use_C11_Richardson,              // Out
                                        bool& l_use_shear_Richardson,            // Out
                                        bool& l_brunt_vaisala_freq_moist,        // Out
                                        bool& l_use_thvm_in_bv_freq,             // Out
                                        bool& l_rcm_supersat_adj,                // Out
                                        bool& l_damp_wp3_Skw_squared,            // Out
                                        bool& l_prescribed_avg_deltaz,           // Out
                                        bool& l_lmm_stepping,                    // Out
                                        bool& l_e3sm_config,                     // Out
                                        bool& l_vary_convect_depth,              // Out
                                        bool& l_use_tke_in_wp3_pr_turb_term,     // Out
                                        bool& l_use_tke_in_wp2_wp3_K_dfsn,       // Out
                                        bool& l_use_wp3_lim_with_smth_Heaviside, // Out
                                        bool& l_smooth_Heaviside_tau_wpxp,       // Out
                                        bool& l_modify_limiters_for_cnvg_test,   // Out
                                        bool& l_enable_relaxed_clipping,         // Out
                                        bool& l_linearize_pbl_winds,             // Out
                                        bool& l_mono_flux_lim_thlm,              // Out
                                        bool& l_mono_flux_lim_rtm,               // Out
                                        bool& l_mono_flux_lim_um,                // Out
                                        bool& l_mono_flux_lim_vm,                // Out
                                        bool& l_mono_flux_lim_spikefix,          // Out
                                        bool& l_host_applies_sfc_fluxes,         // Out
                                        bool& l_wp2_fill_holes_tke,              // Out
                                        bool& l_add_dycore_grid                  // Out
) {
  // Description:
  //   Set all configurable CLUBB flags to their Native defaults.
  // References: None.
  // iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use for
  // the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component
  // PDF. [-].
  // ipdf_call_placement: Selected option for the placement of the call to CLUBB's PDF. [-].
  // penta_solve_method: Option to set the penta-diagonal matrix solving method [-].
  // tridiag_solve_method: Option to set the tri-diagonal matrix solving method [-].
  // saturation_formula: Integer that stores the saturation formula to be used [-].
  // Lscale_method: 1=legacy parcels, 2=tau, 3=coherent parcels [-].
  // grid_remap_method: Integer that stores what remapping technique should be used to remap from
  // one grid to another [-].
  // grid_adapt_in_time_method: Integer that stores how the grid should be adapted every timestep
  // or if the grid should not be adapted at all [-].
  // fill_holes_type: Option for which type of hole filler to use in the fill_holes_vertical
  // procedure [-].
  // l_use_precip_frac: Flag to use precipitation fraction in KK microphysics. The precipitation
  // fraction is automatically set to 1 when this flag is turned off. [-].
  // l_predict_upwp_vpwp: Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside the
  // advancement of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'> in subroutine
  // advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still approximated by eddy diffusivity when
  // <u> and <v> are advanced in subroutine advance_windm_edsclrm. [-].
  // l_ho_nontrad_coriolis: Flag to implement the nontraditional Coriolis terms in the prognostic
  // equations of <w'w'>, <u'w'>, and <u'u'>. [-].
  // l_ho_trad_coriolis: Flag to implement the traditional Coriolis terms in the prognostic
  // equations of <v'w'> and <u'w'>. [-].
  // l_min_wp2_from_corr_wx: Flag to base the threshold minimum value of wp2 on keeping the overall
  // correlation of w and x (w and rt, as well as w and theta-l) within the limits of
  // -max_mag_correlation_flux to max_mag_correlation_flux. [-].
  // l_min_xp2_from_corr_wx: Flag to base the threshold minimum value of xp2 (rtp2 and thlp2) on
  // keeping the overall correlation of w and x within the limits of -max_mag_correlation_flux to
  // max_mag_correlation_flux. [-].
  // l_C2_cloud_frac: Flag to use cloud fraction to adjust the value of the turbulent dissipation
  // coefficient, C2. [-].
  // l_diffuse_rtm_and_thlm: Diffuses rtm and thlm [-].
  // l_stability_correct_Kh_N2_zm: Divides Kh_N2_zm by a stability factor [-].
  // l_calc_thlp2_rad: Include the contribution of radiation to thlp2 [-].
  // l_upwind_xpyp_ta: This flag determines whether we want to use an upwind differencing
  // approximation rather than a centered differencing for turbulent or mean advection terms. It
  // affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, & sclrpthlp. [-].
  // l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
  // approximation rather than a centered differencing for turbulent or mean advection terms. It
  // affects rtm, thlm, sclrm, um and vm. [-].
  // l_uv_nudge: For wind speed nudging. [-].
  // l_rtm_nudge: For rtm nudging [-].
  // l_vert_avg_closure: Use 2 calls to pdf_closure and the trapezoidal rule to compute the
  // varibles that are output from high order closure [-].
  // l_trapezoidal_rule_zt: If true, the trapezoidal rule is called for the thermodynamic-level
  // variables output from pdf_closure. [-].
  // l_trapezoidal_rule_zm: If true, the trapezoidal rule is called for three momentum-level
  // variables - wpthvp, thlpthvp, and rtpthvp - output from pdf_closure. [-].
  // l_call_pdf_closure_twice: This logical flag determines whether or not to call subroutine
  // pdf_closure twice.  If true, pdf_closure is called first on thermodynamic levels and then on
  // momentum levels so that each variable is computed on its native level.  If false, pdf_closure
  // is only called on thermodynamic levels, and variables which belong on momentum levels are
  // interpolated. [-].
  // l_standard_term_ta: Use the standard discretization for the turbulent advection terms.
  // Setting to .false. means that a_1 and a_3 are pulled outside of the derivative in
  // advance_wp2_wp3_module.F90 and in advance_xp2_xpyp_module.F90. [-].
  // l_partial_upwind_wp3: Flag to use an "upwind" discretization rather than a centered
  // discretization for the portion of the wp3 turbulent advection term for ADG1 that is linearized
  // in terms of wp3<t+1>. (Requires ADG1 PDF and l_standard_term_ta). [-].
  // l_godunov_upwind_wpxp_ta: This flag determines whether we want to use an upwind differencing
  // approximation rather than a centered differencing for turbulent advection terms. It affects
  // wpxp only. [-].
  // l_godunov_upwind_xpyp_ta: This flag determines whether we want to use an upwind differencing
  // approximation rather than a centered differencing for turbulent advection terms. It affects
  // xpyp only. [-].
  // l_use_cloud_cover: Use cloud_cover and rcm_in_layer to help boost cloud_frac and rcm to help
  // increase cloudiness at coarser grid resolutions. [-].
  // l_diagnose_correlations: Diagnose correlations instead of using fixed ones [-].
  // l_calc_w_corr: Calculate the correlations between w and the hydrometeors [-].
  // l_const_Nc_in_cloud: Use a constant cloud droplet conc. within cloud (K&K) [-].
  // l_fix_w_chi_eta_correlations: Use a fixed correlation for s and t Mellor(chi/eta) [-].
  // l_stability_correct_tau_zm: Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability correction
  // [-].
  // l_damp_wp2_using_em: In wp2 equation, use a dissipation formula of -(2/3)*em/tau_zm, as in
  // Bougeault (1981) [-].
  // l_do_expldiff_rtm_thlm: Diffuse rtm and thlm explicitly [-].
  // l_Lscale_plume_centered: Alternate that uses the PDF to compute the perturbed values [-].
  // l_use_C7_Richardson: Parameterize C7 based on Richardson number [-].
  // l_use_C11_Richardson: Parameterize C11 and C16 based on Richardson number [-].
  // l_use_shear_Richardson: Use shear in the calculation of Richardson number [-].
  // l_brunt_vaisala_freq_moist: Use a different formula for the Brunt-Vaisala frequency in
  // saturated atmospheres (from Durran and Klemp, 1982) [-].
  // l_use_thvm_in_bv_freq: Use thvm in the calculation of Brunt-Vaisala frequency [-].
  // l_rcm_supersat_adj: Add excess supersaturated vapor to cloud water [-].
  // l_damp_wp3_Skw_squared: Set damping on wp3 to use Skw^2 rather than Skw^4 [-].
  // l_prescribed_avg_deltaz: If .true., avg_deltaz = deltaz [-].
  // l_lmm_stepping: Apply Linear Multistep Method (LMM) Stepping [-].
  // l_e3sm_config: Run model with E3SM settings [-].
  // l_vary_convect_depth: Flag used to calculate convective velocity using a variable estimate of
  // layer depth based on the depth over which wpthlp is positive near the ground when true More
  // information can be found by Looking at issue #905 on the clubb repo [-].
  // l_use_tke_in_wp3_pr_turb_term: Use TKE formulation for wp3 pr_turb term [-].
  // l_use_tke_in_wp2_wp3_K_dfsn: Use TKE in eddy diffusion for wp2 and wp3 [-].
  // l_use_wp3_lim_with_smth_Heaviside: Flag to activate mods on wp3 limiters for conv test [-].
  // l_smooth_Heaviside_tau_wpxp: Use smoothed Heaviside 'Peskin' function in the calculation of
  // H_invrs_tau_wpxp_N2 in src/CLUBB_core/mixing_length.F90 [-].
  // l_modify_limiters_for_cnvg_test: Flag to activate mods on limiters for conv test [-].
  // l_enable_relaxed_clipping: Flag to relax clipping on wpxp in xm_wpxp_clipping_and_stats [-].
  // l_linearize_pbl_winds: Code to linearize PBL winds [-].
  // l_mono_flux_lim_thlm: Flag to turn on monotonic flux limiter for thlm [-].
  // l_mono_flux_lim_rtm: Flag to turn on monotonic flux limiter for rtm [-].
  // l_mono_flux_lim_um: Flag to turn on monotonic flux limiter for um [-].
  // l_mono_flux_lim_vm: Flag to turn on monotonic flux limiter for vm [-].
  // l_mono_flux_lim_spikefix: Flag to implement monotonic flux limiter code that eliminates
  // spurious drying tendencies at model top [-].
  // l_host_applies_sfc_fluxes: Use to determine whether a host model has already applied the
  // surface flux, to avoid double counting. [-].
  // l_wp2_fill_holes_tke: Turn on additional hole-filling for wp2 that takes TKE from up2 and vp2,
  // if necessary [-].
  // l_add_dycore_grid: Turn on remapping from the dycore grid [-].

  // WARNING: THE DEFAULT VALUES OF THE FLAGS BELOW MAY BE OVERWRITTEN
  // BY NAMELIST VALUES FROM, E.G., configurable_clubb_flags_nl!!!

  iiPDF_type = iiPDF_ADG1;
  ipdf_call_placement = ipdf_post_advance_fields;
  penta_solve_method = penta_lu;
  tridiag_solve_method = tridiag_lu;
  saturation_formula = saturation_flatau;
  Lscale_method = Lscale_coherence;
  grid_remap_method = ppm_remap;
  grid_adapt_in_time_method = no_grid_adaptation;
  fill_holes_type = sliding_window;
  l_use_precip_frac = true;
  l_predict_upwp_vpwp = true;
  l_ho_nontrad_coriolis = false;
  l_ho_trad_coriolis = false;
  l_min_wp2_from_corr_wx = false;
  l_min_xp2_from_corr_wx = true;
  l_C2_cloud_frac = false;
  l_diffuse_rtm_and_thlm = false;
  l_stability_correct_Kh_N2_zm = false;
  l_calc_thlp2_rad = true;
  l_upwind_xpyp_ta = true;
  l_upwind_xm_ma = true;
  l_uv_nudge = false;
  l_rtm_nudge = false;
  l_vert_avg_closure = false;
  l_trapezoidal_rule_zt = false;
  l_trapezoidal_rule_zm = false;
  l_call_pdf_closure_twice = false;
  l_standard_term_ta = false;
  l_partial_upwind_wp3 = false;
  l_godunov_upwind_wpxp_ta = false;
  l_godunov_upwind_xpyp_ta = false;
  l_use_cloud_cover = false;
  l_diagnose_correlations = false;
  l_calc_w_corr = false;
  l_const_Nc_in_cloud = false;
  l_fix_w_chi_eta_correlations = true;
  l_stability_correct_tau_zm = false;
  l_damp_wp2_using_em = true;
  l_do_expldiff_rtm_thlm = false;
  l_Lscale_plume_centered = false;
  l_use_C7_Richardson = true;
  l_use_C11_Richardson = false;
  l_use_shear_Richardson = false;
  l_brunt_vaisala_freq_moist = false;
  l_use_thvm_in_bv_freq = false;
  l_rcm_supersat_adj = true;
  l_damp_wp3_Skw_squared = true;
  l_prescribed_avg_deltaz = false;
  l_lmm_stepping = false;
  l_e3sm_config = false;
  l_vary_convect_depth = false;
  l_use_tke_in_wp3_pr_turb_term = true;
  l_use_tke_in_wp2_wp3_K_dfsn = false;
  l_use_wp3_lim_with_smth_Heaviside = false;
  l_smooth_Heaviside_tau_wpxp = false;
  l_modify_limiters_for_cnvg_test = false;
  l_enable_relaxed_clipping = false;
  l_linearize_pbl_winds = false;
  l_mono_flux_lim_thlm = true;
  l_mono_flux_lim_rtm = true;
  l_mono_flux_lim_um = true;
  l_mono_flux_lim_vm = true;
  l_mono_flux_lim_spikefix = true;
  l_host_applies_sfc_fluxes = false;
  l_wp2_fill_holes_tke = true;
  l_add_dycore_grid = false;

  return;
}

// =============================================================================
void initialize_clubb_config_flags_type_api(int iiPDF_type,                             // In
                                            int ipdf_call_placement,                    // In
                                            int penta_solve_method,                     // In
                                            int tridiag_solve_method,                   // In
                                            int saturation_formula,                     // In
                                            int Lscale_method,                          // In
                                            int grid_remap_method,                      // In
                                            int grid_adapt_in_time_method,              // In
                                            int fill_holes_type,                        // In
                                            bool l_use_precip_frac,                     // In
                                            bool l_predict_upwp_vpwp,                   // In
                                            bool l_ho_nontrad_coriolis,                 // In
                                            bool l_ho_trad_coriolis,                    // In
                                            bool l_min_wp2_from_corr_wx,                // In
                                            bool l_min_xp2_from_corr_wx,                // In
                                            bool l_C2_cloud_frac,                       // In
                                            bool l_diffuse_rtm_and_thlm,                // In
                                            bool l_stability_correct_Kh_N2_zm,          // In
                                            bool l_calc_thlp2_rad,                      // In
                                            bool l_upwind_xpyp_ta,                      // In
                                            bool l_upwind_xm_ma,                        // In
                                            bool l_uv_nudge,                            // In
                                            bool l_rtm_nudge,                           // In
                                            bool l_vert_avg_closure,                    // In
                                            bool l_trapezoidal_rule_zt,                 // In
                                            bool l_trapezoidal_rule_zm,                 // In
                                            bool l_call_pdf_closure_twice,              // In
                                            bool l_standard_term_ta,                    // In
                                            bool l_partial_upwind_wp3,                  // In
                                            bool l_godunov_upwind_wpxp_ta,              // In
                                            bool l_godunov_upwind_xpyp_ta,              // In
                                            bool l_use_cloud_cover,                     // In
                                            bool l_diagnose_correlations,               // In
                                            bool l_calc_w_corr,                         // In
                                            bool l_const_Nc_in_cloud,                   // In
                                            bool l_fix_w_chi_eta_correlations,          // In
                                            bool l_stability_correct_tau_zm,            // In
                                            bool l_damp_wp2_using_em,                   // In
                                            bool l_do_expldiff_rtm_thlm,                // In
                                            bool l_Lscale_plume_centered,               // In
                                            bool l_use_C7_Richardson,                   // In
                                            bool l_use_C11_Richardson,                  // In
                                            bool l_use_shear_Richardson,                // In
                                            bool l_brunt_vaisala_freq_moist,            // In
                                            bool l_use_thvm_in_bv_freq,                 // In
                                            bool l_rcm_supersat_adj,                    // In
                                            bool l_damp_wp3_Skw_squared,                // In
                                            bool l_prescribed_avg_deltaz,               // In
                                            bool l_lmm_stepping,                        // In
                                            bool l_e3sm_config,                         // In
                                            bool l_vary_convect_depth,                  // In
                                            bool l_use_tke_in_wp3_pr_turb_term,         // In
                                            bool l_use_tke_in_wp2_wp3_K_dfsn,           // In
                                            bool l_use_wp3_lim_with_smth_Heaviside,     // In
                                            bool l_smooth_Heaviside_tau_wpxp,           // In
                                            bool l_modify_limiters_for_cnvg_test,       // In
                                            bool l_enable_relaxed_clipping,             // In
                                            bool l_linearize_pbl_winds,                 // In
                                            bool l_mono_flux_lim_thlm,                  // In
                                            bool l_mono_flux_lim_rtm,                   // In
                                            bool l_mono_flux_lim_um,                    // In
                                            bool l_mono_flux_lim_vm,                    // In
                                            bool l_mono_flux_lim_spikefix,              // In
                                            bool l_host_applies_sfc_fluxes,             // In
                                            bool l_wp2_fill_holes_tke,                  // In
                                            bool l_add_dycore_grid,                     // In
                                            clubb_config_flags_type& clubb_config_flags // Out
) {
  // Description:
  //   Validate selector bounds and initialize the caller-owned configuration type.
  // References: None.
  // iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use for
  // the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component
  // PDF. [-].
  // ipdf_call_placement: Selected option for the placement of the call to CLUBB's PDF. [-].
  // penta_solve_method: Option to set the penta-diagonal matrix solving method [-].
  // tridiag_solve_method: Option to set the tri-diagonal matrix solving method [-].
  // saturation_formula: Integer that stores the saturation formula to be used [-].
  // Lscale_method: 1=legacy parcels, 2=tau, 3=coherent parcels [-].
  // grid_remap_method: Integer that stores what remapping technique should be used to remap from
  // one grid to another (starts at 1, so 0 is an invalid option for this flag) [-].
  // grid_adapt_in_time_method: Integer that stores how the grid should be adapted every timestep
  // or if the grid should not be adapted at all [-].
  // fill_holes_type: Option for which type of hole filler to use in the fill_holes_vertical
  // procedure [-].
  // l_use_precip_frac: Flag to use precipitation fraction in KK microphysics. The precipitation
  // fraction is automatically set to 1 when this flag is turned off. [-].
  // l_predict_upwp_vpwp: Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside the
  // advancement of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>, and <w'sclr'> in subroutine
  // advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still approximated by eddy diffusivity when
  // <u> and <v> are advanced in subroutine advance_windm_edsclrm. [-].
  // l_ho_nontrad_coriolis: Flag to implement the nontraditional Coriolis terms in the prognostic
  // equations of <w'w'>, <u'w'>, and <u'u'>. [-].
  // l_ho_trad_coriolis: Flag to implement the traditional Coriolis terms in the prognostic
  // equations of <v'w'> and <u'w'>. [-].
  // l_min_wp2_from_corr_wx: Flag to base the threshold minimum value of wp2 on keeping the overall
  // correlation of w and x (w and rt, as well as w and theta-l) within the limits of
  // -max_mag_correlation_flux to max_mag_correlation_flux. [-].
  // l_min_xp2_from_corr_wx: Flag to base the threshold minimum value of xp2 (rtp2 and thlp2) on
  // keeping the overall correlation of w and x within the limits of -max_mag_correlation_flux to
  // max_mag_correlation_flux. [-].
  // l_C2_cloud_frac: Flag to use cloud fraction to adjust the value of the turbulent dissipation
  // coefficient, C2. [-].
  // l_diffuse_rtm_and_thlm: Diffuses rtm and thlm [-].
  // l_stability_correct_Kh_N2_zm: Divides Kh_N2_zm by a stability factor [-].
  // l_calc_thlp2_rad: Include the contribution of radiation to thlp2 [-].
  // l_upwind_xpyp_ta: This flag determines whether we want to use an upwind differencing
  // approximation rather than a centered differencing for turbulent or mean advection terms. It
  // affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, & sclrpthlp. [-].
  // l_upwind_xm_ma: This flag determines whether we want to use an upwind differencing
  // approximation rather than a centered differencing for turbulent or mean advection terms. It
  // affects rtm, thlm, sclrm, um and vm. [-].
  // l_uv_nudge: For wind speed nudging. [-].
  // l_rtm_nudge: For rtm nudging [-].
  // l_vert_avg_closure: Use 2 calls to pdf_closure and the trapezoidal rule to compute the
  // varibles that are output from high order closure [-].
  // l_trapezoidal_rule_zt: If true, the trapezoidal rule is called for the thermodynamic-level
  // variables output from pdf_closure. [-].
  // l_trapezoidal_rule_zm: If true, the trapezoidal rule is called for three momentum-level
  // variables - wpthvp, thlpthvp, and rtpthvp - output from pdf_closure. [-].
  // l_call_pdf_closure_twice: This logical flag determines whether or not to call subroutine
  // pdf_closure twice.  If true, pdf_closure is called first on thermodynamic levels and then on
  // momentum levels so that each variable is computed on its native level.  If false, pdf_closure
  // is only called on thermodynamic levels, and variables which belong on momentum levels are
  // interpolated. [-].
  // l_standard_term_ta: Use the standard discretization for the turbulent advection terms.
  // Setting to .false. means that a_1 and a_3 are pulled outside of the derivative in
  // advance_wp2_wp3_module.F90 and in advance_xp2_xpyp_module.F90. [-].
  // l_partial_upwind_wp3: Flag to use an "upwind" discretization rather than a centered
  // discretization for the portion of the wp3 turbulent advection term for ADG1 that is linearized
  // in terms of wp3<t+1>. (Requires ADG1 PDF and l_standard_term_ta). [-].
  // l_godunov_upwind_wpxp_ta: This flag determines whether we want to use an upwind differencing
  // approximation rather than a centered differencing for turbulent advection terms. It affects
  // wpxp only. [-].
  // l_godunov_upwind_xpyp_ta: This flag determines whether we want to use an upwind differencing
  // approximation rather than a centered differencing for turbulent advection terms. It affects
  // xpyp only. [-].
  // l_use_cloud_cover: Use cloud_cover and rcm_in_layer to help boost cloud_frac and rcm to help
  // increase cloudiness at coarser grid resolutions. [-].
  // l_diagnose_correlations: Diagnose correlations instead of using fixed ones [-].
  // l_calc_w_corr: Calculate the correlations between w and the hydrometeors [-].
  // l_const_Nc_in_cloud: Use a constant cloud droplet conc. within cloud (K&K) [-].
  // l_fix_w_chi_eta_correlations: Use a fixed correlation for s and t Mellor(chi/eta) [-].
  // l_stability_correct_tau_zm: Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability correction
  // [-].
  // l_damp_wp2_using_em: In wp2 equation, use a dissipation formula of -(2/3)*em/tau_zm, as in
  // Bougeault (1981) [-].
  // l_do_expldiff_rtm_thlm: Diffuse rtm and thlm explicitly [-].
  // l_Lscale_plume_centered: Alternate that uses the PDF to compute the perturbed values [-].
  // l_use_C7_Richardson: Parameterize C7 based on Richardson number [-].
  // l_use_C11_Richardson: Parameterize C11 and C16 based on Richardson number [-].
  // l_use_shear_Richardson: Use shear in the calculation of Richardson number [-].
  // l_brunt_vaisala_freq_moist: Use a different formula for the Brunt-Vaisala frequency in
  // saturated atmospheres (from Durran and Klemp, 1982) [-].
  // l_use_thvm_in_bv_freq: Use thvm in the calculation of Brunt-Vaisala frequency [-].
  // l_rcm_supersat_adj: Add excess supersaturated vapor to cloud water [-].
  // l_damp_wp3_Skw_squared: Set damping on wp3 to use Skw^2 rather than Skw^4 [-].
  // l_prescribed_avg_deltaz: If .true., avg_deltaz = deltaz [-].
  // l_lmm_stepping: Apply Linear Multistep Method (LMM) Stepping [-].
  // l_e3sm_config: Run model with E3SM settings [-].
  // l_vary_convect_depth: Flag used to calculate convective velocity using a variable estimate of
  // layer depth based on the depth over which wpthlp is positive near the ground when true More
  // information can be found by Looking at issue #905 on the clubb repo [-].
  // l_use_tke_in_wp3_pr_turb_term: Use TKE formulation for wp3 pr_turb term [-].
  // l_use_tke_in_wp2_wp3_K_dfsn: Use TKE in eddy diffusion for wp2 and wp3 [-].
  // l_use_wp3_lim_with_smth_Heaviside: Flag to activate mods on wp3 limiters for conv test [-].
  // l_smooth_Heaviside_tau_wpxp: Use smoothed Heaviside 'Peskin' function in the calculation of
  // H_invrs_tau_wpxp_N2 in src/CLUBB_core/mixing_length.F90 [-].
  // l_modify_limiters_for_cnvg_test: Flag to activate mods on limiters for conv test [-].
  // l_enable_relaxed_clipping: Flag to relax clipping on wpxp in xm_wpxp_clipping_and_stats [-].
  // l_linearize_pbl_winds: Code to linearize PBL winds [-].
  // l_mono_flux_lim_thlm: Flag to turn on monotonic flux limiter for thlm [-].
  // l_mono_flux_lim_rtm: Flag to turn on monotonic flux limiter for rtm [-].
  // l_mono_flux_lim_um: Flag to turn on monotonic flux limiter for um [-].
  // l_mono_flux_lim_vm: Flag to turn on monotonic flux limiter for vm [-].
  // l_mono_flux_lim_spikefix: Flag to implement monotonic flux limiter code that eliminates
  // spurious drying tendencies at model top [-].
  // l_host_applies_sfc_fluxes: Use to determine whether a host model has already applied the
  // surface flux, to avoid double counting. [-].
  // l_wp2_fill_holes_tke: Turn on additional hole-filling for wp2 that takes TKE from up2 and vp2,
  // if necessary [-].
  // l_add_dycore_grid: Turn on remapping from the dycore grid [-].
  // clubb_config_flags: Derived type holding all configurable CLUBB flags [-].

  if (grid_remap_method <= 0) {
    throw std::runtime_error(
        "Invalid value for the flag grid_remap_method. Should be greater or equal to 1.");
  }

  if (Lscale_method < Lscale_direct || Lscale_method > Lscale_coherence) {
    throw std::runtime_error("Invalid Lscale_method. Expected an integer from 1 through 3.");
  }

  clubb_config_flags.iiPDF_type = iiPDF_type;
  clubb_config_flags.ipdf_call_placement = ipdf_call_placement;
  clubb_config_flags.penta_solve_method = penta_solve_method;
  clubb_config_flags.tridiag_solve_method = tridiag_solve_method;
  clubb_config_flags.saturation_formula = saturation_formula;
  clubb_config_flags.Lscale_method = Lscale_method;
  clubb_config_flags.grid_remap_method = grid_remap_method;
  clubb_config_flags.grid_adapt_in_time_method = grid_adapt_in_time_method;
  clubb_config_flags.fill_holes_type = fill_holes_type;
  clubb_config_flags.l_use_precip_frac = l_use_precip_frac;
  clubb_config_flags.l_predict_upwp_vpwp = l_predict_upwp_vpwp;
  clubb_config_flags.l_ho_nontrad_coriolis = l_ho_nontrad_coriolis;
  clubb_config_flags.l_ho_trad_coriolis = l_ho_trad_coriolis;
  clubb_config_flags.l_min_wp2_from_corr_wx = l_min_wp2_from_corr_wx;
  clubb_config_flags.l_min_xp2_from_corr_wx = l_min_xp2_from_corr_wx;
  clubb_config_flags.l_C2_cloud_frac = l_C2_cloud_frac;
  clubb_config_flags.l_diffuse_rtm_and_thlm = l_diffuse_rtm_and_thlm;
  clubb_config_flags.l_stability_correct_Kh_N2_zm = l_stability_correct_Kh_N2_zm;
  clubb_config_flags.l_calc_thlp2_rad = l_calc_thlp2_rad;
  clubb_config_flags.l_upwind_xpyp_ta = l_upwind_xpyp_ta;
  clubb_config_flags.l_upwind_xm_ma = l_upwind_xm_ma;
  clubb_config_flags.l_uv_nudge = l_uv_nudge;
  clubb_config_flags.l_rtm_nudge = l_rtm_nudge;
  clubb_config_flags.l_vert_avg_closure = l_vert_avg_closure;
  clubb_config_flags.l_trapezoidal_rule_zt = l_trapezoidal_rule_zt;
  clubb_config_flags.l_trapezoidal_rule_zm = l_trapezoidal_rule_zm;
  clubb_config_flags.l_call_pdf_closure_twice = l_call_pdf_closure_twice;
  clubb_config_flags.l_standard_term_ta = l_standard_term_ta;
  clubb_config_flags.l_partial_upwind_wp3 = l_partial_upwind_wp3;
  clubb_config_flags.l_godunov_upwind_wpxp_ta = l_godunov_upwind_wpxp_ta;
  clubb_config_flags.l_godunov_upwind_xpyp_ta = l_godunov_upwind_xpyp_ta;
  clubb_config_flags.l_use_cloud_cover = l_use_cloud_cover;
  clubb_config_flags.l_diagnose_correlations = l_diagnose_correlations;
  clubb_config_flags.l_calc_w_corr = l_calc_w_corr;
  clubb_config_flags.l_const_Nc_in_cloud = l_const_Nc_in_cloud;
  clubb_config_flags.l_fix_w_chi_eta_correlations = l_fix_w_chi_eta_correlations;
  clubb_config_flags.l_stability_correct_tau_zm = l_stability_correct_tau_zm;
  clubb_config_flags.l_damp_wp2_using_em = l_damp_wp2_using_em;
  clubb_config_flags.l_do_expldiff_rtm_thlm = l_do_expldiff_rtm_thlm;
  clubb_config_flags.l_Lscale_plume_centered = l_Lscale_plume_centered;
  clubb_config_flags.l_use_C7_Richardson = l_use_C7_Richardson;
  clubb_config_flags.l_use_C11_Richardson = l_use_C11_Richardson;
  clubb_config_flags.l_use_shear_Richardson = l_use_shear_Richardson;
  clubb_config_flags.l_brunt_vaisala_freq_moist = l_brunt_vaisala_freq_moist;
  clubb_config_flags.l_use_thvm_in_bv_freq = l_use_thvm_in_bv_freq;
  clubb_config_flags.l_rcm_supersat_adj = l_rcm_supersat_adj;
  clubb_config_flags.l_damp_wp3_Skw_squared = l_damp_wp3_Skw_squared;
  clubb_config_flags.l_prescribed_avg_deltaz = l_prescribed_avg_deltaz;
  clubb_config_flags.l_lmm_stepping = l_lmm_stepping;
  clubb_config_flags.l_e3sm_config = l_e3sm_config;
  clubb_config_flags.l_vary_convect_depth = l_vary_convect_depth;
  clubb_config_flags.l_use_tke_in_wp3_pr_turb_term = l_use_tke_in_wp3_pr_turb_term;
  clubb_config_flags.l_use_tke_in_wp2_wp3_K_dfsn = l_use_tke_in_wp2_wp3_K_dfsn;
  clubb_config_flags.l_use_wp3_lim_with_smth_Heaviside = l_use_wp3_lim_with_smth_Heaviside;
  clubb_config_flags.l_smooth_Heaviside_tau_wpxp = l_smooth_Heaviside_tau_wpxp;
  clubb_config_flags.l_modify_limiters_for_cnvg_test = l_modify_limiters_for_cnvg_test;
  clubb_config_flags.l_enable_relaxed_clipping = l_enable_relaxed_clipping;
  clubb_config_flags.l_linearize_pbl_winds = l_linearize_pbl_winds;
  clubb_config_flags.l_mono_flux_lim_thlm = l_mono_flux_lim_thlm;
  clubb_config_flags.l_mono_flux_lim_rtm = l_mono_flux_lim_rtm;
  clubb_config_flags.l_mono_flux_lim_um = l_mono_flux_lim_um;
  clubb_config_flags.l_mono_flux_lim_vm = l_mono_flux_lim_vm;
  clubb_config_flags.l_mono_flux_lim_spikefix = l_mono_flux_lim_spikefix;
  clubb_config_flags.l_host_applies_sfc_fluxes = l_host_applies_sfc_fluxes;
  clubb_config_flags.l_wp2_fill_holes_tke = l_wp2_fill_holes_tke;
  clubb_config_flags.l_add_dycore_grid = l_add_dycore_grid;

  return;
}

// =============================================================================
void print_clubb_config_flags_api(int iunit, const clubb_config_flags_type& clubb_config_flags) {
  // Description:
  //   Print flags to a Native logical unit through the transient I/O service.
  // References: None.
  // iunit: The file to write to [-].
  // clubb_config_flags: Derived type holding all configurable CLUBB flags [-].

  // Preserve the original integer-unit public contract and exact Native formatting.
  f2k::print_clubb_config_flags_api(iunit, clubb_config_flags);
}

// =============================================================================
void print_clubb_config_flags_api(std::ostream& iunit,
                                  const clubb_config_flags_type& clubb_config_flags) {
  // Description:
  //   Print the Native flag rows to the existing C++ case-owned output stream.
  // iunit: Output stream replacing the Native logical-unit selector at this boundary.
  // clubb_config_flags: Caller-owned configurable CLUBB flags [-].
  // Keep integer/logical rendering local, independent of the caller's stream flags.
  // Native list-directed spacing is compiler-specific; row/value mapping is exact.
  std::ostringstream output;
  output << "iiPDF_type = " << clubb_config_flags.iiPDF_type << '\n';
  output << "ipdf_call_placement = " << clubb_config_flags.ipdf_call_placement << '\n';
  output << "penta_solve_method = " << clubb_config_flags.penta_solve_method << '\n';
  output << "tridiag_solve_method = " << clubb_config_flags.tridiag_solve_method << '\n';
  output << "Lscale_method = " << clubb_config_flags.Lscale_method << '\n';
  output << "grid_remap_method = " << clubb_config_flags.grid_remap_method << '\n';
  output << "grid_adapt_in_time_method = " << clubb_config_flags.grid_adapt_in_time_method << '\n';
  output << "l_use_precip_frac = " << (clubb_config_flags.l_use_precip_frac ? "T" : "F") << '\n';
  output << "l_predict_upwp_vpwp = " << (clubb_config_flags.l_predict_upwp_vpwp ? "T" : "F")
         << '\n';
  output << "l_ho_nontrad_coriolis = " << (clubb_config_flags.l_ho_nontrad_coriolis ? "T" : "F")
         << '\n';
  output << "l_ho_trad_coriolis = " << (clubb_config_flags.l_ho_trad_coriolis ? "T" : "F") << '\n';
  output << "l_min_wp2_from_corr_wx = " << (clubb_config_flags.l_min_wp2_from_corr_wx ? "T" : "F")
         << '\n';
  output << "l_min_xp2_from_corr_wx = " << (clubb_config_flags.l_min_xp2_from_corr_wx ? "T" : "F")
         << '\n';
  output << "l_C2_cloud_frac = " << (clubb_config_flags.l_C2_cloud_frac ? "T" : "F") << '\n';
  output << "l_diffuse_rtm_and_thlm = " << (clubb_config_flags.l_diffuse_rtm_and_thlm ? "T" : "F")
         << '\n';
  output << "l_stability_correct_Kh_N2_zm = "
         << (clubb_config_flags.l_stability_correct_Kh_N2_zm ? "T" : "F") << '\n';
  output << "l_calc_thlp2_rad = " << (clubb_config_flags.l_calc_thlp2_rad ? "T" : "F") << '\n';
  output << "l_upwind_xpyp_ta = " << (clubb_config_flags.l_upwind_xpyp_ta ? "T" : "F") << '\n';
  output << "l_upwind_xm_ma = " << (clubb_config_flags.l_upwind_xm_ma ? "T" : "F") << '\n';
  output << "l_uv_nudge = " << (clubb_config_flags.l_uv_nudge ? "T" : "F") << '\n';
  output << "l_rtm_nudge = " << (clubb_config_flags.l_rtm_nudge ? "T" : "F") << '\n';
  output << "l_vert_avg_closure = " << (clubb_config_flags.l_vert_avg_closure ? "T" : "F") << '\n';
  output << "l_trapezoidal_rule_zt = " << (clubb_config_flags.l_trapezoidal_rule_zt ? "T" : "F")
         << '\n';
  output << "l_trapezoidal_rule_zm = " << (clubb_config_flags.l_trapezoidal_rule_zm ? "T" : "F")
         << '\n';
  output << "l_call_pdf_closure_twice = "
         << (clubb_config_flags.l_call_pdf_closure_twice ? "T" : "F") << '\n';
  output << "l_standard_term_ta = " << (clubb_config_flags.l_standard_term_ta ? "T" : "F") << '\n';
  output << "l_partial_upwind_wp3 = " << (clubb_config_flags.l_partial_upwind_wp3 ? "T" : "F")
         << '\n';
  output << "l_godunov_upwind_wpxp_ta = "
         << (clubb_config_flags.l_godunov_upwind_wpxp_ta ? "T" : "F") << '\n';
  output << "l_godunov_upwind_xpyp_ta = "
         << (clubb_config_flags.l_godunov_upwind_xpyp_ta ? "T" : "F") << '\n';
  output << "l_use_cloud_cover = " << (clubb_config_flags.l_use_cloud_cover ? "T" : "F") << '\n';
  output << "l_diagnose_correlations = " << (clubb_config_flags.l_diagnose_correlations ? "T" : "F")
         << '\n';
  output << "l_calc_w_corr = " << (clubb_config_flags.l_calc_w_corr ? "T" : "F") << '\n';
  output << "l_const_Nc_in_cloud = " << (clubb_config_flags.l_const_Nc_in_cloud ? "T" : "F")
         << '\n';
  output << "l_fix_w_chi_eta_correlations = "
         << (clubb_config_flags.l_fix_w_chi_eta_correlations ? "T" : "F") << '\n';
  output << "l_stability_correct_tau_zm = "
         << (clubb_config_flags.l_stability_correct_tau_zm ? "T" : "F") << '\n';
  output << "l_damp_wp2_using_em = " << (clubb_config_flags.l_damp_wp2_using_em ? "T" : "F")
         << '\n';
  output << "l_do_expldiff_rtm_thlm = " << (clubb_config_flags.l_do_expldiff_rtm_thlm ? "T" : "F")
         << '\n';
  output << "l_Lscale_plume_centered = " << (clubb_config_flags.l_Lscale_plume_centered ? "T" : "F")
         << '\n';
  output << "l_use_C7_Richardson = " << (clubb_config_flags.l_use_C7_Richardson ? "T" : "F")
         << '\n';
  output << "l_use_C11_Richardson = " << (clubb_config_flags.l_use_C11_Richardson ? "T" : "F")
         << '\n';
  output << "l_use_shear_Richardson = " << (clubb_config_flags.l_use_shear_Richardson ? "T" : "F")
         << '\n';
  output << "l_brunt_vaisala_freq_moist = "
         << (clubb_config_flags.l_brunt_vaisala_freq_moist ? "T" : "F") << '\n';
  output << "l_use_thvm_in_bv_freq = " << (clubb_config_flags.l_use_thvm_in_bv_freq ? "T" : "F")
         << '\n';
  output << "l_rcm_supersat_adj = " << (clubb_config_flags.l_rcm_supersat_adj ? "T" : "F") << '\n';
  output << "l_damp_wp3_Skw_squared = " << (clubb_config_flags.l_damp_wp3_Skw_squared ? "T" : "F")
         << '\n';
  output << "l_prescribed_avg_deltaz = " << (clubb_config_flags.l_prescribed_avg_deltaz ? "T" : "F")
         << '\n';
  output << "l_lmm_stepping = " << (clubb_config_flags.l_lmm_stepping ? "T" : "F") << '\n';
  output << "l_e3sm_config = " << (clubb_config_flags.l_e3sm_config ? "T" : "F") << '\n';
  output << "l_vary_convect_depth = " << (clubb_config_flags.l_vary_convect_depth ? "T" : "F")
         << '\n';
  output << "l_use_tke_in_wp3_pr_turb_term = "
         << (clubb_config_flags.l_use_tke_in_wp3_pr_turb_term ? "T" : "F") << '\n';
  output << "l_use_tke_in_wp2_wp3_K_dfsn = "
         << (clubb_config_flags.l_use_tke_in_wp2_wp3_K_dfsn ? "T" : "F") << '\n';
  output << "l_use_wp3_lim_with_smth_Heaviside = "
         << (clubb_config_flags.l_use_wp3_lim_with_smth_Heaviside ? "T" : "F") << '\n';
  output << "l_smooth_Heaviside_tau_wpxp = "
         << (clubb_config_flags.l_smooth_Heaviside_tau_wpxp ? "T" : "F") << '\n';
  output << "l_modify_limiters_for_cnvg_test = "
         << (clubb_config_flags.l_modify_limiters_for_cnvg_test ? "T" : "F") << '\n';
  output << "l_enable_relaxed_clipping = "
         << (clubb_config_flags.l_enable_relaxed_clipping ? "T" : "F") << '\n';
  output << "l_linearize_pbl_winds = " << (clubb_config_flags.l_linearize_pbl_winds ? "T" : "F")
         << '\n';
  output << "l_mono_flux_lim_thlm = " << (clubb_config_flags.l_mono_flux_lim_thlm ? "T" : "F")
         << '\n';
  output << "l_mono_flux_lim_rtm = " << (clubb_config_flags.l_mono_flux_lim_rtm ? "T" : "F")
         << '\n';
  output << "l_mono_flux_lim_um = " << (clubb_config_flags.l_mono_flux_lim_vm ? "T" : "F") << '\n';
  output << "l_mono_flux_lim_vm = " << (clubb_config_flags.l_mono_flux_lim_um ? "T" : "F") << '\n';
  output << "l_mono_flux_lim_spikefix = "
         << (clubb_config_flags.l_mono_flux_lim_spikefix ? "T" : "F") << '\n';
  output << "l_host_applies_sfc_fluxes = "
         << (clubb_config_flags.l_host_applies_sfc_fluxes ? "T" : "F") << '\n';
  output << "l_wp2_fill_holes_tke = " << (clubb_config_flags.l_wp2_fill_holes_tke ? "T" : "F")
         << '\n';
  output << "l_add_dycore_grid = " << (clubb_config_flags.l_add_dycore_grid ? "T" : "F") << '\n';
  iunit << output.str();
}

} // namespace clubb_kokkos
