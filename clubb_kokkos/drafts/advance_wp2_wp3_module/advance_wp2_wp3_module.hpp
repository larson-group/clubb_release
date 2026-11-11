// Advance CLUBB vertical-velocity variance wp2 and third moment wp3.
//
// advance_wp2_wp3 assembles their coupled semi-implicit equations, invokes the
// existing band solver, applies Native hole-filling/clipping and records budgets.
// Physical fields use column-first LayoutLeft Views; local coefficient arrays
// retain Native band-first shapes until the explicit column-first solver boundary.
// Native index metadata stays one-based; parameter indices are already zero-based.
// Sponge-enabled cases remain outside the existing JAX/Kokkos driver scope.
// See INTEGRATION.md for dependencies and pending numerical/build evidence.
// References: Golaz et al.2002,JAS59,3540–3551,Eqs12/18; Equations for CLUBB sec6.
#pragma once

#include "array_types.hpp"
#include "grid_class.hpp"
#include "pdf_parameter_module.hpp"
#include "parameters_tunable.hpp"
#include "err_info_type_module.hpp"
#include "f2k/stats_bridge.hpp"

namespace clubb_kokkos {

void advance_wp2_wp3(int nzm, int nzt, int ngrdcol, const grid& gr, real_t dt,              // In
                     Real1D sfc_elevation, Real1D fcor_y, Real2D sigma_sqd_w, Real2D wm_zm, // In
                     Real2D wm_zt,                                                          // In
                     Real2D wpup2, Real2D wpvp2, Real2D wp2up2, Real2D wp2vp2, Real2D wp4,  // In
                     Real2D wpthvp, Real2D wp2thvp, Real2D wp2up, Real2D um, Real2D vm,     // In
                     Real2D upwp, Real2D vpwp,                                              // In
                     Real2D em, Real2D Kh_zm, Real2D Kh_zt, Real2D invrs_tau_C4_zm,         // In
                     Real2D invrs_tau_wp3_zt, Real2D invrs_tau_C1_zm,                       // In
                     Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zm,            // In
                     Real2D invrs_rho_ds_zt, Real2D thv_ds_zm,                              // In
                     Real2D thv_ds_zt, Real2D mixt_frac, Real2D Cx_fnc_Richardson,          // In
                     Real2D lhs_splat_wp2, Real2D lhs_splat_wp3,                            // In
                     const implicit_coefs_terms& pdf_implicit_coefs_terms,                  // In
                     Real2D wprtp, Real2D wpthlp, Real2D rtp2, Real2D thlp2,                // In
                     Real2D clubb_params, const nu_vertical_res_dep& nu_vert_res_dep,       // In
                     int iiPDF_type,                                                        // In
                     int penta_solve_method,                                                // In
                     int fill_holes_type,                                                   // In
                     bool l_min_wp2_from_corr_wx,                                           // In
                     bool l_upwind_xm_ma,                                                   // In
                     bool l_standard_term_ta,                                               // In
                     bool l_partial_upwind_wp3,                                             // In
                     bool l_damp_wp2_using_em,                                              // In
                     bool l_use_C11_Richardson,                                             // In
                     bool l_damp_wp3_Skw_squared,                                           // In
                     bool l_lmm_stepping,                                                   // In
                     bool l_use_tke_in_wp3_pr_turb_term,                                    // In
                     bool l_use_tke_in_wp2_wp3_K_dfsn,                                      // In
                     bool l_use_wp3_lim_with_smth_Heaviside,                                // In
                     bool l_wp2_fill_holes_tke,                                             // In
                     bool l_ho_nontrad_coriolis,                                            // In
                     bool l_implemented,                                                    // In
                     f2k::stats_type& stats,                                                // InOut
                     Real2D up2, Real2D vp2, Real2D wp2, Real2D wp3,                        // InOut
                     err_info_type& err_info);                                              // InOut
} // namespace clubb_kokkos
