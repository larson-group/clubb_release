// Coordinate CLUBB PDF closure on thermodynamic and momentum grid levels.
//
// pdf_closure selects the probability density function (PDF), diagnoses its
// component moments, and computes cloud, water and buoyancy moments on one grid.
// pdf_closure_driver coordinates both grids, optional vertical averaging,
// nudging, cloud-water adjustments and sampled output to C++ statistics.
//
// calc_wp4_pdf, calc_wp2xp_pdf, calc_wpxp2_pdf and calc_wpxpyp_pdf integrate
// fourth/third moments from the component means, variances and correlations.
// calc_w_up_in_cloud diagnoses cloudy updraft/downdraft velocities and fractions.
//
// Arrays use execution-space, column-first LayoutLeft Views with native field
// units. Grid boundary metadata retains its native one-based values. Scratch
// belongs to each call; err_info carries native fatal status to the caller.
// PDF algorithms and diagnostics call the permanent C++ owners.
#pragma once

#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "grid_class.hpp"
#include "pdf_parameter_module.hpp"
#include "err_info_type_module.hpp"
#include "stats_netcdf.hpp"

namespace clubb_kokkos {

void pdf_closure(ScratchWorkspace& workspace,                                               // InOut
                 int nz, int ngrdcol, int sclr_dim, Real1D sclr_tol, const grid& gr,        // In
                 int hydromet_dim, Real2D p_in_Pa, Real2D exner, Real2D thv_ds,             // In
                 Real2D wm, Real2D wp2, Real2D wp3,                                         // In
                 Real2D Skw, Real2D Skthl_in, Real2D Skrt_in, Real2D Sku_in, Real2D Skv_in, // In
                 Real2D rtm, Real2D rtp2, Real2D wprtp,                                     // In
                 Real2D thlm, Real2D thlp2, Real2D wpthlp,                                  // In
                 Real2D um, Real2D up2, Real2D upwp,                                        // In
                 Real2D vm, Real2D vp2, Real2D vpwp,                                        // In
                 Real2D rtpthlp,                                                            // In
                 Real3D sclrm, Real3D wpsclrp, Real3D sclrp2,                               // In
                 Real3D sclrprtp, Real3D sclrpthlp, Real3D Sksclr_in,                       // In
                 Real3D wphydrometp, Real3D wp2hmp,                                         // In
                 Real3D rtphmp, Real3D thlphmp,                                             // In
                 Real2D clubb_params, Real1D mixt_frac_max_mag,                             // In
                 int saturation_formula,                                                    // In
                 stats_type& stats,                                                         // InOut
                 int iiPDF_type,                                                            // In
                 Bool1D l_mix_rat_hm,                                                       // In
                 Real2D sigma_sqd_w,                                                        // InOut
                 pdf_parameter& pdf_params, implicit_coefs_terms& pdf_implicit_coefs_terms, // InOut
                 err_info_type& err_info,                                                   // InOut
                 Real2D wpup2, Real2D wpvp2,                                                // Out
                 Real2D wp2up2, Real2D wp2vp2, Real2D wp4,                                  // Out
                 Real2D wprtp2, Real2D wp2rtp,                                              // Out
                 Real2D wpthlp2, Real2D wp2thlp, Real2D wprtpthlp,                          // Out
                 Real2D cloud_frac, Real2D ice_supersat_frac,                               // Out
                 Real2D rcm, Real2D wpthvp, Real2D wp2thvp, Real2D wp2up, Real2D rtpthvp,   // Out
                 Real2D thlpthvp, Real2D wprcp, Real2D wp2rcp, Real2D rtprcp,               // Out
                 Real2D thlprcp, Real2D rcp2,                                               // Out
                 Real2D uprcp, Real2D vprcp,                                                // Out
                 Real2D w_up_in_cloud, Real2D w_down_in_cloud,                              // Out
                 Real2D cloudy_updraft_frac, Real2D cloudy_downdraft_frac,                  // Out
                 Real3D wpsclrprtp, Real3D wpsclrp2, Real3D sclrpthvp,                      // Out
                 Real3D wpsclrpthlp, Real3D sclrprcp, Real3D wp2sclrp,                      // Out
                 Real2D rc_coef);

void calc_wp4_pdf(int nz, int ngrdcol,                  // In
                  Real2D wm, Real2D w_1, Real2D w_2,    // In
                  Real2D varnce_w_1, Real2D varnce_w_2, // In
                  Real2D mixt_frac,                     // In
                  Real2D wp4);

void calc_wp2xp_pdf(int nz, int ngrdcol,                          // In
                    Real2D wm, Real2D xm, Real2D w_1, Real2D w_2, // In
                    Real2D x_1, Real2D x_2,                       // In
                    Real2D varnce_w_1, Real2D varnce_w_2,         // In
                    Real2D varnce_x_1, Real2D varnce_x_2,         // In
                    Real2D corr_w_x_1, Real2D corr_w_x_2,         // In
                    Real2D mixt_frac,                             // In
                    Real2D wp2xp);

void calc_wpxp2_pdf(int nz, int ngrdcol,                  // In
                    Real2D wm, Real2D xm, Real2D w_1,     // In
                    Real2D w_2, Real2D x_1, Real2D x_2,   // In
                    Real2D varnce_w_1, Real2D varnce_w_2, // In
                    Real2D varnce_x_1, Real2D varnce_x_2, // In
                    Real2D corr_w_x_1, Real2D corr_w_x_2, // In
                    Real2D mixt_frac,                     // In
                    Real2D wpxp2);

void calc_wpxpyp_pdf(int nz, int ngrdcol,                                     // In
                     Real2D wm, Real2D xm, Real2D ym, Real2D w_1, Real2D w_2, // In
                     Real2D x_1, Real2D x_2,                                  // In
                     Real2D y_1, Real2D y_2,                                  // In
                     Real2D varnce_w_1, Real2D varnce_w_2,                    // In
                     Real2D varnce_x_1, Real2D varnce_x_2,                    // In
                     Real2D varnce_y_1, Real2D varnce_y_2,                    // In
                     Real2D corr_w_x_1, Real2D corr_w_x_2,                    // In
                     Real2D corr_w_y_1, Real2D corr_w_y_2,                    // In
                     Real2D corr_x_y_1, Real2D corr_x_y_2,                    // In
                     Real2D mixt_frac,                                        // In
                     Real2D wpxpyp);

void calc_w_up_in_cloud(int nz, int ngrdcol,                                          // In
                        Real2D mixt_frac, Real2D cloud_frac_1, Real2D cloud_frac_2,   // In
                        Real2D w_1, Real2D w_2, Real2D varnce_w_1, Real2D varnce_w_2, // In
                        Real2D w_up_in_cloud, Real2D w_down_in_cloud,                 // Out
                        Real2D cloudy_updraft_frac, Real2D cloudy_downdraft_frac);    // Out

void pdf_closure_driver(ScratchWorkspace& workspace,                                        // InOut
                        const grid& gr, int nzm, int nzt, int ngrdcol,                      // In
                        real_t dt, int hydromet_dim, int sclr_dim, Real1D sclr_tol,         // In
                        Real2D wprtp, Real2D thlm, Real2D wpthlp, Real2D rtp2, Real2D rtp3, // In
                        Real2D thlp2, Real2D thlp3, Real2D rtpthlp, Real2D wp2,             // In
                        Real2D wp3, Real2D wm_zm, Real2D wm_zt,                             // In
                        Real2D um, Real2D up2, Real2D upwp, Real2D up3,                     // In
                        Real2D vm, Real2D vp2, Real2D vpwp, Real2D vp3,                     // In
                        Real2D p_in_Pa, Real2D exner,                                       // In
                        Real2D thv_ds_zm, Real2D thv_ds_zt, Real2D rtm_ref,                 // In
                        Real3D wphydrometp,                                                 // In
                        Real3D wp2hmp, Real3D rtphmp_zt, Real3D thlphmp_zt,                 // In
                        Real3D sclrm, Real3D wpsclrp, Real3D sclrp2,                        // In
                        Real3D sclrprtp, Real3D sclrpthlp, Real3D sclrp3,                   // In
                        Real1D p_sfc, bool l_samp_stats_in_pdf_call,                        // In
                        Real1D mixt_frac_max_mag, real_t ts_nudge,                          // In
                        real_t rtm_min, real_t rtm_nudge_max_altitude,                      // In
                        Real2D clubb_params,                                                // In
                        int iiPDF_type,                                                     // In
                        int saturation_formula,                                             // In
                        bool l_rtm_nudge,                                                   // In
                        bool l_trapezoidal_rule_zt,                                         // In
                        bool l_trapezoidal_rule_zm,                                         // In
                        bool l_call_pdf_closure_twice,                                      // In
                        bool l_use_cloud_cover,                                             // In
                        bool l_rcm_supersat_adj,                                            // In
                        Bool1D l_mix_rat_hm,                                                // In
                        stats_type& stats,                                                  // InOut
                        Real2D rtm, Real2D sigma_sqd_w,                                     // InOut
                        implicit_coefs_terms& pdf_implicit_coefs_terms,                     // InOut
                        pdf_parameter& pdf_params, pdf_parameter& pdf_params_zm,
                        err_info_type& err_info,                      // InOut
                        Real2D rcm, Real2D cloud_frac,                // Out
                        Real2D ice_supersat_frac, Real2D wprcp,       // Out
                        Real2D wpthvp, Real2D wp2thvp, Real2D wp2up,  // Out
                        Real2D rtpthvp, Real2D thlpthvp,              // Out
                        Real2D rcm_in_layer, Real2D cloud_cover,      // Out
                        Real2D rcp2_zt, Real2D thlprcp,               // Out
                        Real2D rc_coef_zm, Real3D sclrpthvp,          // Out
                        Real2D wpup2, Real2D wpvp2,                   // Out
                        Real2D wp2up2, Real2D wp2vp2, Real2D wp4,     // Out
                        Real2D wp2rtp, Real2D wprtp2, Real2D wp2thlp, // Out
                        Real2D wpthlp2, Real2D wprtpthlp,             // Out
                        Real2D uprcp, Real2D vprcp,                   // Out
                        Real2D w_up_in_cloud, Real2D w_down_in_cloud, // Out
                        Real2D cloudy_updraft_frac,                   // Out
                        Real2D cloudy_downdraft_frac,                 // Out
                        Real3D wp2sclrp, Real3D wpsclrp2,             // Out
                        Real3D wpsclrprtp, Real3D wpsclrpthlp);

} // namespace clubb_kokkos
