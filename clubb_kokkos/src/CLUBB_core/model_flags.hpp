// Define, initialize and print CLUBB's case-owned configuration flags.
//
// set_default_clubb_config_flags_api writes scalar defaults before namelist
// overrides. initialize_clubb_config_flags_type_api validates the remapping and
// mixing-length selectors, then assembles clubb_config_flags_type.
// print_clubb_config_flags_api writes the Native rows to a caller-owned C++ stream.
// Immutable selectors and the type schema are generated from model_flags.F90;
// each caller owns its mutable configuration. Native print row/value order is
// retained, including the source's swapped um/vm labels.
#pragma once

#include <Kokkos_Core.hpp>
#include <iosfwd>
#include <cstdint>
#include <string>
#include "clubb_real.hpp"

namespace clubb_kokkos {

// BEGIN GENERATED MODEL FLAGS TYPE
// Native compile-time options used by the ported core.
inline constexpr int order_xm_wpxp = 1;
inline constexpr int order_xp2_xpyp = 2;
inline constexpr int order_wp2_wp3 = 3;
inline constexpr int order_windm = 4;
inline constexpr int iiPDF_ADG1 = 1;
inline constexpr int iiPDF_ADG2 = 2;
inline constexpr int iiPDF_3D_Luhar = 3;
inline constexpr int iiPDF_new = 4;
inline constexpr int iiPDF_TSDADG = 5;
inline constexpr int iiPDF_LY93 = 6;
inline constexpr int iiPDF_new_hybrid = 7;
inline constexpr int saturation_bolton = 1;
inline constexpr int saturation_gfdl = 2;
inline constexpr int saturation_flatau = 3;
inline constexpr int saturation_lookup = 4;
inline constexpr bool I_sat_sphum = false;
inline constexpr int ipdf_pre_advance_fields = 1;
inline constexpr int ipdf_post_advance_fields = 2;
inline constexpr int ipdf_pre_post_advance_fields = 3;
inline constexpr bool l_advance_xp3 = false;
inline constexpr bool l_morr_xp2_mc = false;
inline constexpr bool l_evaporate_cold_rcm = false;
inline constexpr bool l_silhs_rad = false;
inline constexpr bool l_use_boussinesq = false;
inline constexpr bool l_test_grid_generalization = false;
inline constexpr bool l_explicit_turbulent_adv_wp3 = false;
inline constexpr bool l_upwind_Kh_dp_term = false;
inline constexpr int global_fill = 1;
inline constexpr int sliding_window = 2;
inline constexpr bool l_pos_def = false;
inline constexpr bool l_clip_turb_adv = false;
inline constexpr bool l_explicit_turbulent_adv_wpxp = false;
inline constexpr bool l_force_descending_solves = false;
inline constexpr bool l_explicit_turbulent_adv_xpyp = false;
inline constexpr bool l_gamma_Skw = true;
inline constexpr bool l_cubic_interp = false;
inline constexpr bool l_quintic_poly_interp = false;
inline constexpr int Lscale_direct = 1;
inline constexpr int Lscale_tau = 2;
inline constexpr int Lscale_coherence = 3;
inline constexpr int cons_ullrich_remap = 1;
inline constexpr int ppm_remap = 2;
inline constexpr int no_grid_adaptation = 0;
inline constexpr int lapack = 1;
inline constexpr int penta_lu = 2;
inline constexpr int tridiag_lu = 2;
inline constexpr int penta_bicgstab = 3;

// Mirrors src/CLUBB_core/model_flags.F90; arrays use Fortran-order LayoutLeft.
struct clubb_config_flags_type {
  int iiPDF_type{};
  int ipdf_call_placement{};
  int penta_solve_method{};
  int tridiag_solve_method{};
  int saturation_formula{};
  int Lscale_method{};
  int grid_remap_method{};
  int grid_adapt_in_time_method{};
  int fill_holes_type{};
  std::uint8_t l_use_precip_frac{};
  std::uint8_t l_predict_upwp_vpwp{};
  std::uint8_t l_ho_nontrad_coriolis{};
  std::uint8_t l_ho_trad_coriolis{};
  std::uint8_t l_min_wp2_from_corr_wx{};
  std::uint8_t l_min_xp2_from_corr_wx{};
  std::uint8_t l_C2_cloud_frac{};
  std::uint8_t l_diffuse_rtm_and_thlm{};
  std::uint8_t l_stability_correct_Kh_N2_zm{};
  std::uint8_t l_calc_thlp2_rad{};
  std::uint8_t l_upwind_xpyp_ta{};
  std::uint8_t l_upwind_xm_ma{};
  std::uint8_t l_uv_nudge{};
  std::uint8_t l_rtm_nudge{};
  std::uint8_t l_vert_avg_closure{};
  std::uint8_t l_trapezoidal_rule_zt{};
  std::uint8_t l_trapezoidal_rule_zm{};
  std::uint8_t l_call_pdf_closure_twice{};
  std::uint8_t l_standard_term_ta{};
  std::uint8_t l_partial_upwind_wp3{};
  std::uint8_t l_godunov_upwind_wpxp_ta{};
  std::uint8_t l_godunov_upwind_xpyp_ta{};
  std::uint8_t l_use_cloud_cover{};
  std::uint8_t l_diagnose_correlations{};
  std::uint8_t l_calc_w_corr{};
  std::uint8_t l_const_Nc_in_cloud{};
  std::uint8_t l_fix_w_chi_eta_correlations{};
  std::uint8_t l_stability_correct_tau_zm{};
  std::uint8_t l_damp_wp2_using_em{};
  std::uint8_t l_do_expldiff_rtm_thlm{};
  std::uint8_t l_Lscale_plume_centered{};
  std::uint8_t l_use_C7_Richardson{};
  std::uint8_t l_use_C11_Richardson{};
  std::uint8_t l_use_shear_Richardson{};
  std::uint8_t l_brunt_vaisala_freq_moist{};
  std::uint8_t l_use_thvm_in_bv_freq{};
  std::uint8_t l_rcm_supersat_adj{};
  std::uint8_t l_damp_wp3_Skw_squared{};
  std::uint8_t l_prescribed_avg_deltaz{};
  std::uint8_t l_lmm_stepping{};
  std::uint8_t l_e3sm_config{};
  std::uint8_t l_vary_convect_depth{};
  std::uint8_t l_use_tke_in_wp3_pr_turb_term{};
  std::uint8_t l_use_tke_in_wp2_wp3_K_dfsn{};
  std::uint8_t l_use_wp3_lim_with_smth_Heaviside{};
  std::uint8_t l_smooth_Heaviside_tau_wpxp{};
  std::uint8_t l_modify_limiters_for_cnvg_test{};
  std::uint8_t l_enable_relaxed_clipping{};
  std::uint8_t l_linearize_pbl_winds{};
  std::uint8_t l_mono_flux_lim_thlm{};
  std::uint8_t l_mono_flux_lim_rtm{};
  std::uint8_t l_mono_flux_lim_um{};
  std::uint8_t l_mono_flux_lim_vm{};
  std::uint8_t l_mono_flux_lim_spikefix{};
  std::uint8_t l_host_applies_sfc_fluxes{};
  std::uint8_t l_wp2_fill_holes_tke{};
  std::uint8_t l_add_dycore_grid{};
};

// END GENERATED MODEL FLAGS TYPE

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
);

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
);

// iunit is the C++ output stream replacing the Native logical-unit selector.
void print_clubb_config_flags_api(std::ostream& iunit,
                                  const clubb_config_flags_type& clubb_config_flags);

}  // namespace clubb_kokkos
