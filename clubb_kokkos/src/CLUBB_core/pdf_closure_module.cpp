// PDF closure orchestration and higher moments;
// see pdf_closure_module.hpp for the public interface.
#include "scratch_workspace.hpp"
#include "numerical_check.hpp"
#include "clip_explicit.hpp"
#include "pdf_closure_module.hpp"
#include "saturation.hpp"
#include "Skx_module.hpp"
#include "adg1_adg2_3d_luhar_pdf.hpp"
#include "new_pdf_main.hpp"
#include "new_hybrid_pdf_main.hpp"
#include "LY93_pdf.hpp"
#include "new_tsdadg_pdf.hpp"
#include "T_in_K_module.hpp"
#include "constants_clubb.hpp"
#include "error_code.hpp"
#include "model_flags.hpp"
#include "parameter_indices.hpp"
#include "stats_netcdf.hpp"
#include "pdf_utilities.hpp"
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace clubb_kokkos {

namespace {

void transform_pdf_chi_eta_component(ScratchWorkspace& workspace,                      // InOut
                                     int nz, int ngrdcol,                              // In
                                     Real2D tl, Real2D rsatl, Real2D rt, Real2D exner, // In
                                     Real2D varnce_thl, Real2D varnce_rt,              // In
                                     Real2D corr_rt_thl, Real2D chi,                   // In/Out
                                     Real2D crt, Real2D cthl,                          // Out
                                     Real2D stdev_chi, Real2D stdev_eta,               // Out
                                     Real2D covar_chi_eta,                             // Out
                                     Real2D corr_chi_eta);

void calc_wp2xp2_pdf(int nz, int ngrdcol,                  // In
                     Real2D wm, Real2D xm, Real2D w_1,     // In
                     Real2D w_2, Real2D x_1, Real2D x_2,   // In
                     Real2D varnce_w_1, Real2D varnce_w_2, // In
                     Real2D varnce_x_1, Real2D varnce_x_2, // In
                     Real2D corr_w_x_1, Real2D corr_w_x_2, // In
                     Real2D mixt_frac,                     // In
                     Real2D wp2xp2);

void calc_liquid_cloud_frac_component(int nz, int ngrdcol,               // In
                                      Real2D mean_chi, Real2D stdev_chi, // In
                                      Real2D cloud_frac, Real2D rc);

void calc_ice_cloud_frac_component(int nz, int ngrdcol,               // In
                                   Real2D mean_chi, Real2D stdev_chi, // In
                                   Real2D rc_in, Real2D cloud_frac,   // In
                                   Real2D p_in_Pa, Real2D tl,         // In
                                   Real2D rsatl, Real2D crt,          // In
                                   int saturation_formula,            // In
                                   Real2D ice_supersat_frac, Real2D rc);

void calc_xprcp_component(int nz, int ngrdcol, // In
                          Real2D wm, Real2D rtm, Real2D thlm, Real2D um, Real2D vm,
                          Real2D rcm,                                                  // In
                          Real2D w_i, Real2D rt_i,                                     // In
                          Real2D thl_i, Real2D u_i, Real2D v_i,                        // In
                          Real2D varnce_w_i, Real2D chi_i,                             // In
                          Real2D stdev_chi_i, Real2D stdev_eta_i,                      // In
                          Real2D corr_w_chi_i, Real2D corr_chi_eta_i,                  // In
                          Real2D crt_i, Real2D cthl_i,                                 // In
                          Real2D rc_i, Real2D cloud_frac_i, int iiPDF_type,            // In
                          Real2D wprcp_contrib_comp_i, Real2D wp2rcp_contrib_comp_i,   // Out
                          Real2D rtprcp_contrib_comp_i, Real2D thlprcp_contrib_comp_i, // Out
                          Real2D uprcp_contrib_comp_i, Real2D vprcp_contrib_comp_i);

void pdf_closure_driver_zm(ScratchWorkspace& workspace,                                     // InOut
                           const grid& gr, int nzm, int nzt, int ngrdcol, int hydromet_dim, // In
                           int sclr_dim, Real1D sclr_tol, Real1D p_sfc,                     // In
                           Real1D mixt_frac_max_mag, Real2D clubb_params,                   // In
                           int iiPDF_type, int saturation_formula,                          // In
                           Bool1D l_mix_rat_hm, Real2D p_in_Pa, Real2D thv_ds_zm,           // In
                           Real2D wm_zm, Real2D wp2, Real2D wp3_zm,                         // In
                           Real2D Skw_zm, Real2D Skthl_zm, Real2D Skrt_zm,                  // In
                           Real2D rtm, Real2D rtp2, Real2D wprtp,                           // In
                           Real2D thlm, Real2D thlp2, Real2D wpthlp,                        // In
                           Real2D um, Real2D up2, Real2D upwp, Real2D up3,                  // In
                           Real2D vm, Real2D vp2, Real2D vpwp, Real2D vp3,                  // In
                           Real2D rtpthlp, Real3D sclrm, Real3D wpsclrp, Real3D sclrp2,     // In
                           Real3D sclrprtp, Real3D sclrpthlp, Real3D sclrp3,                // In
                           Real3D wphydrometp,                                              // In
                           Real3D wp2hmp, Real3D rtphmp_zt, Real3D thlphmp_zt,              // In
                           stats_type& stats, Real2D sigma_sqd_w,
                           pdf_parameter& pdf_params_zm,                              // InOut
                           err_info_type& err_info,                                   // InOut
                           Real2D rtm_zm, Real2D thlm_zm,                             // Out
                           Real2D wp2up2, Real2D wp2vp2, Real2D wp4,                  // Out
                           Real2D wprtp2_zm, Real2D wp2rtp_zm,                        // Out
                           Real2D wpthlp2_zm, Real2D wp2thlp_zm, Real2D wprtpthlp_zm, // Out
                           Real2D cloud_frac_zm, Real2D ice_supersat_frac_zm,         // Out
                           Real2D rcm_zm, Real2D wpthvp, Real2D wp2thvp_zm, Real2D wp2up_zm, // Out
                           Real2D rtpthvp, Real2D thlpthvp, Real2D wprcp, Real2D wp2rcp_zm,  // Out
                           Real2D rtprcp, Real2D thlprcp, Real2D rcp2, Real2D uprcp,
                           Real2D vprcp,                                               // Out
                           Real3D wpsclrprtp_zm, Real3D wpsclrp2_zm, Real3D sclrpthvp, // Out
                           Real3D wpsclrpthlp_zm, Real3D sclrprcp, Real3D wp2sclrp_zm, // Out
                           Real2D rc_coef_zm);

void trapezoidal_rule_zt(int nzm, int nzt, int ngrdcol, int sclr_dim, const grid& gr,   // In
                         bool l_call_pdf_closure_twice,                                 // In
                         const stats_type& stats,                                       // In
                         Real2D wprtp2, Real2D wpthlp2,                                 // InOut
                         Real2D wprtpthlp, Real2D cloud_frac, Real2D ice_supersat_frac, // InOut
                         Real2D rcm, Real2D wp2thvp, Real2D wp2up, Real3D wpsclrprtp,
                         Real3D wpsclrp2,                                               // InOut
                         Real3D wpsclrpthlp,                                            // InOut
                         Real2D wprtp2_zm, Real2D wpthlp2_zm,                           // InOut
                         Real2D wprtpthlp_zm, Real2D cloud_frac_zm,                     // InOut
                         Real2D ice_supersat_frac_zm, Real2D rcm_zm, Real2D wp2thvp_zm, // InOut
                         Real2D wp2up_zm,                                               // InOut
                         Real3D wpsclrprtp_zm, Real3D wpsclrp2_zm, Real3D wpsclrpthlp_zm);

void trapezoidal_rule_zm(int nzm, int nzt, int ngrdcol, const grid& gr,           // In
                         Real2D wpthvp_zt, Real2D thlpthvp_zt, Real2D rtpthvp_zt, // In
                         Real2D wpthvp, Real2D thlpthvp, Real2D rtpthvp);

void calc_trapezoid_zt(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                       Real2D variable_zm,                            // In
                       Real2D variable_zt);

void calc_trapezoid_zm(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D variable_zt, // In
                       Real2D variable_zm);

void compute_cloud_cover(ScratchWorkspace& workspace,                                    // InOut
                         const grid& gr, int nzt, int ngrdcol,                           // In
                         const pdf_parameter& pdf_params, Real2D cloud_frac, Real2D rcm, // In
                         err_info_type& err_info,                                        // InOut
                         Real2D cloud_cover, Real2D rcm_in_layer);

} // namespace

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
                 Real2D rc_coef) {
  // Description:
  // Subroutine that computes pdf parameters analytically.
  //
  // Based of the original formulation, but with some tweaks
  // to remove some of the less realistic assumptions and
  // improve transport terms.

  //   Corrected version that should remove inconsistency

  // References:
  //   The shape of CLUBB's PDF is given by the expression in
  //   https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:clubb_pdf

  //   Eqn. 29, 30, 31, 32 & 33  on p. 3547 of
  //   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
  //   Method and Model Description'' Golaz, et al. (2002)
  //   JAS, Vol. 59, pp. 3540--3551.
  //----------------------------------------------------------------------

  //----------------------------- Local Variables -----------------------------
  auto workspace_scope = workspace.scope();
  Real2D u_1           = workspace_scope.view<Real2D>(
      "pdf_closure::u_1", ngrdcol, nz); // Mean of eastward wind (1st PDF component)         [m/s]
  Real2D u_2 = workspace_scope.view<Real2D>(
      "pdf_closure::u_2", ngrdcol, nz); // Mean of eastward wind (2nd PDF component)         [m/s]
  Real2D varnce_u_1 =
      workspace_scope.view<Real2D>("pdf_closure::varnce_u_1", ngrdcol,
                                   nz); // Variance of u (1st PDF component)             [m^2/s^2]
  Real2D varnce_u_2 =
      workspace_scope.view<Real2D>("pdf_closure::varnce_u_2", ngrdcol,
                                   nz); // Variance of u (2nd PDF component)             [m^2/s^2]
  Real2D v_1 = workspace_scope.view<Real2D>(
      "pdf_closure::v_1", ngrdcol, nz); // Mean of northward wind (1st PDF component)        [m/s]
  Real2D v_2 = workspace_scope.view<Real2D>(
      "pdf_closure::v_2", ngrdcol, nz); // Mean of northward wind (2nd PDF component)        [m/s]
  Real2D varnce_v_1 =
      workspace_scope.view<Real2D>("pdf_closure::varnce_v_1", ngrdcol,
                                   nz); // Variance of v (1st PDF component)             [m^2/s^2]
  Real2D varnce_v_2 =
      workspace_scope.view<Real2D>("pdf_closure::varnce_v_2", ngrdcol,
                                   nz); // Variance of v (2nd PDF component)             [m^2/s^2]
  Real2D alpha_u =
      workspace_scope.view<Real2D>("pdf_closure::alpha_u", ngrdcol,
                                   nz); // Factor relating to normalized variance for u        [-]
  Real2D alpha_v =
      workspace_scope.view<Real2D>("pdf_closure::alpha_v", ngrdcol,
                                   nz); // Factor relating to normalized variance for v        [-]
  Real2D corr_u_w_1 =
      workspace_scope.view<Real2D>("pdf_closure::corr_u_w_1", ngrdcol,
                                   nz); // Correlation of u and w   (1st PDF component)      [-]
  Real2D corr_u_w_2 =
      workspace_scope.view<Real2D>("pdf_closure::corr_u_w_2", ngrdcol,
                                   nz); // Correlation of u and w   (2nd PDF component)      [-]
  Real2D corr_v_w_1 =
      workspace_scope.view<Real2D>("pdf_closure::corr_v_w_1", ngrdcol,
                                   nz); // Correlation of v and w   (1st PDF component)      [-]
  Real2D corr_v_w_2 =
      workspace_scope.view<Real2D>("pdf_closure::corr_v_w_2", ngrdcol,
                                   nz); // Correlation of v and w   (2nd PDF component)      [-]
  Real3D sclr1 = workspace_scope.view<Real3D>("pdf_closure::sclr1", ngrdcol, nz, sclr_dim);
  Real3D sclr2 = workspace_scope.view<Real3D>("pdf_closure::sclr2", ngrdcol, nz, sclr_dim);
  Real3D varnce_sclr1 =
      workspace_scope.view<Real3D>("pdf_closure::varnce_sclr1", ngrdcol, nz, sclr_dim);
  Real3D varnce_sclr2 =
      workspace_scope.view<Real3D>("pdf_closure::varnce_sclr2", ngrdcol, nz, sclr_dim);
  Real3D alpha_sclr =
      workspace_scope.view<Real3D>("pdf_closure::alpha_sclr", ngrdcol, nz, sclr_dim);
  Real3D corr_sclr_thl_1 =
      workspace_scope.view<Real3D>("pdf_closure::corr_sclr_thl_1", ngrdcol, nz, sclr_dim);
  Real3D corr_sclr_thl_2 =
      workspace_scope.view<Real3D>("pdf_closure::corr_sclr_thl_2", ngrdcol, nz, sclr_dim);
  Real3D corr_sclr_rt_1 =
      workspace_scope.view<Real3D>("pdf_closure::corr_sclr_rt_1", ngrdcol, nz, sclr_dim);
  Real3D corr_sclr_rt_2 =
      workspace_scope.view<Real3D>("pdf_closure::corr_sclr_rt_2", ngrdcol, nz, sclr_dim);
  Real3D corr_w_sclr_1 =
      workspace_scope.view<Real3D>("pdf_closure::corr_w_sclr_1", ngrdcol, nz, sclr_dim);
  Real3D corr_w_sclr_2 =
      workspace_scope.view<Real3D>("pdf_closure::corr_w_sclr_2", ngrdcol, nz, sclr_dim);
  Real2D tl1      = workspace_scope.view<Real2D>("pdf_closure::tl1", ngrdcol, nz);
  Real2D tl2      = workspace_scope.view<Real2D>("pdf_closure::tl2", ngrdcol, nz);
  Real2D sqrt_wp2 = workspace_scope.view<Real2D>("pdf_closure::sqrt_wp2", ngrdcol,
                                                 nz); // Square root of wp2          [m/s]
  Real2D Skthl    = workspace_scope.view<Real2D>("pdf_closure::Skthl", ngrdcol,
                                                 nz); // Skewness of thl             [-]
  Real2D Skrt     = workspace_scope.view<Real2D>("pdf_closure::Skrt", ngrdcol,
                                                 nz); // Skewness of rt              [-]
  Real2D Sku      = workspace_scope.view<Real2D>("pdf_closure::Sku", ngrdcol,
                                                 nz); // Skewness of u               [-]
  Real2D Skv      = workspace_scope.view<Real2D>("pdf_closure::Skv", ngrdcol,
                                                 nz); // Skewness of v               [-]
  Real3D Sksclr   = workspace_scope.view<Real3D>("pdf_closure::Sksclr", ngrdcol, nz,
                                                 sclr_dim); // Skewness of rt              [-]
  Real2D wprcp_contrib_comp_1 =
      workspace_scope.view<Real2D>("pdf_closure::wprcp_contrib_comp_1", ngrdcol,
                                   nz); // <w'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
  Real2D wprcp_contrib_comp_2 =
      workspace_scope.view<Real2D>("pdf_closure::wprcp_contrib_comp_2", ngrdcol,
                                   nz); // <w'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]
  Real2D wp2rcp_contrib_comp_1 =
      workspace_scope.view<Real2D>("pdf_closure::wp2rcp_contrib_comp_1", ngrdcol,
                                   nz); // <w'^2rc'> contrib. (1st comp) [m^2/s^2(kg/kg)]
  Real2D wp2rcp_contrib_comp_2 =
      workspace_scope.view<Real2D>("pdf_closure::wp2rcp_contrib_comp_2", ngrdcol,
                                   nz); // <w'^2rc'> contrib. (2nd comp) [m^2/s^2(kg/kg)]
  Real2D rtprcp_contrib_comp_1 =
      workspace_scope.view<Real2D>("pdf_closure::rtprcp_contrib_comp_1", ngrdcol,
                                   nz); // <rt'rc'> contrib. (1st PDF comp.)  [kg^2/kg^2]
  Real2D rtprcp_contrib_comp_2 =
      workspace_scope.view<Real2D>("pdf_closure::rtprcp_contrib_comp_2", ngrdcol,
                                   nz); // <rt'rc'> contrib. (2nd PDF comp.)  [kg^2/kg^2]
  Real2D thlprcp_contrib_comp_1 =
      workspace_scope.view<Real2D>("pdf_closure::thlprcp_contrib_comp_1", ngrdcol,
                                   nz); // <thl'rc'> contrib. (1st PDF comp.)  [K(kg/kg)]
  Real2D thlprcp_contrib_comp_2 =
      workspace_scope.view<Real2D>("pdf_closure::thlprcp_contrib_comp_2", ngrdcol,
                                   nz); // <thl'rc'> contrib. (2nd PDF comp.)  [K(kg/kg)]
  Real2D uprcp_contrib_comp_1 =
      workspace_scope.view<Real2D>("pdf_closure::uprcp_contrib_comp_1", ngrdcol,
                                   nz); // <u'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
  Real2D uprcp_contrib_comp_2 =
      workspace_scope.view<Real2D>("pdf_closure::uprcp_contrib_comp_2", ngrdcol,
                                   nz); // <u'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]
  Real2D vprcp_contrib_comp_1 =
      workspace_scope.view<Real2D>("pdf_closure::vprcp_contrib_comp_1", ngrdcol,
                                   nz); // <v'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
  Real2D vprcp_contrib_comp_2 =
      workspace_scope.view<Real2D>("pdf_closure::vprcp_contrib_comp_2", ngrdcol,
                                   nz); // <v'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]
  Real2D rc_1_ice       = workspace_scope.view<Real2D>("pdf_closure::rc_1_ice", ngrdcol, nz);
  Real2D rc_2_ice       = workspace_scope.view<Real2D>("pdf_closure::rc_2_ice", ngrdcol, nz);
  Real2D wm_clubb_pdf   = workspace_scope.view<Real2D>("pdf_closure::wm_clubb_pdf", ngrdcol, nz);
  Real2D rtm_clubb_pdf  = workspace_scope.view<Real2D>("pdf_closure::rtm_clubb_pdf", ngrdcol, nz);
  Real2D thlm_clubb_pdf = workspace_scope.view<Real2D>("pdf_closure::thlm_clubb_pdf", ngrdcol, nz);
  Real2D wp2_clubb_pdf  = workspace_scope.view<Real2D>("pdf_closure::wp2_clubb_pdf", ngrdcol, nz);
  Real2D rtp2_clubb_pdf = workspace_scope.view<Real2D>("pdf_closure::rtp2_clubb_pdf", ngrdcol, nz);
  Real2D thlp2_clubb_pdf =
      workspace_scope.view<Real2D>("pdf_closure::thlp2_clubb_pdf", ngrdcol, nz);
  Real2D wp3_clubb_pdf  = workspace_scope.view<Real2D>("pdf_closure::wp3_clubb_pdf", ngrdcol, nz);
  Real2D rtp3_clubb_pdf = workspace_scope.view<Real2D>("pdf_closure::rtp3_clubb_pdf", ngrdcol, nz);
  Real2D thlp3_clubb_pdf =
      workspace_scope.view<Real2D>("pdf_closure::thlp3_clubb_pdf", ngrdcol, nz);
  Real2D Skw_clubb_pdf  = workspace_scope.view<Real2D>("pdf_closure::Skw_clubb_pdf", ngrdcol, nz);
  Real2D Skrt_clubb_pdf = workspace_scope.view<Real2D>("pdf_closure::Skrt_clubb_pdf", ngrdcol, nz);
  Real2D Skthl_clubb_pdf =
      workspace_scope.view<Real2D>("pdf_closure::Skthl_clubb_pdf", ngrdcol, nz);

  constexpr bool l_liq_ice_loading_test = false; // Native temporary loading-test flag.
  const bool l_scalar_calc              = sclr_dim > 0;

  // Capture only the PDF component Views used by each kernel. Native lower
  // calls fill these same Views without replacing their allocation owners.
  const auto pdf_params_alpha_thl           = pdf_params.alpha_thl;
  const auto pdf_params_alpha_rt            = pdf_params.alpha_rt;
  const auto pdf_params_corr_w_rt_1         = pdf_params.corr_w_rt_1;
  const auto pdf_params_corr_w_rt_2         = pdf_params.corr_w_rt_2;
  const auto pdf_params_corr_w_thl_1        = pdf_params.corr_w_thl_1;
  const auto pdf_params_corr_w_thl_2        = pdf_params.corr_w_thl_2;
  const auto pdf_params_thl_1               = pdf_params.thl_1;
  const auto pdf_params_thl_2               = pdf_params.thl_2;
  const auto pdf_params_mixt_frac           = pdf_params.mixt_frac;
  const auto pdf_params_ice_supersat_frac_1 = pdf_params.ice_supersat_frac_1;
  const auto pdf_params_ice_supersat_frac_2 = pdf_params.ice_supersat_frac_2;
  const auto pdf_params_cloud_frac_1        = pdf_params.cloud_frac_1;
  const auto pdf_params_cloud_frac_2        = pdf_params.cloud_frac_2;
  const auto pdf_params_rc_1                = pdf_params.rc_1;
  const auto pdf_params_rc_2                = pdf_params.rc_2;
  const auto pdf_params_corr_w_chi_1        = pdf_params.corr_w_chi_1;
  const auto pdf_params_corr_w_chi_2        = pdf_params.corr_w_chi_2;
  const auto pdf_params_corr_w_eta_1        = pdf_params.corr_w_eta_1;
  const auto pdf_params_corr_w_eta_2        = pdf_params.corr_w_eta_2;
  const auto pdf_params_varnce_rt_1         = pdf_params.varnce_rt_1;
  const auto pdf_params_varnce_rt_2         = pdf_params.varnce_rt_2;
  const auto pdf_params_varnce_thl_1        = pdf_params.varnce_thl_1;
  const auto pdf_params_varnce_thl_2        = pdf_params.varnce_thl_2;
  const auto pdf_params_crt_1               = pdf_params.crt_1;
  const auto pdf_params_crt_2               = pdf_params.crt_2;
  const auto pdf_params_cthl_1              = pdf_params.cthl_1;
  const auto pdf_params_cthl_2              = pdf_params.cthl_2;
  const auto pdf_params_chi_1               = pdf_params.chi_1;
  const auto pdf_params_stdev_chi_1         = pdf_params.stdev_chi_1;
  const auto pdf_params_chi_2               = pdf_params.chi_2;
  const auto pdf_params_stdev_chi_2         = pdf_params.stdev_chi_2;

  //----------------------------- Begin Code -----------------------------
  // Initialize default alpha coefficients and local skewness copies for
  // closures that may clip skewness without changing caller-owned inputs.
  if (iiPDF_type != iiPDF_ADG1 && iiPDF_type != iiPDF_ADG2) {
    Kokkos::parallel_for(
        "pdf_closure::default_alpha", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i                = cell % ngrdcol;
          const int k                = (cell / ngrdcol);
          pdf_params_alpha_thl(i, k) = one_half;
          pdf_params_alpha_rt(i, k)  = one_half;
        });

    Kokkos::parallel_for(
        "pdf_closure::copy_skewness", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = (cell / ngrdcol);
          // This allows for skewness to be clipped locally without passing the updated
          // value back out.
          Skrt(i, k)  = Skrt_in(i, k);
          Skthl(i, k) = Skthl_in(i, k);
          Sku(i, k)   = Sku_in(i, k);
          Skv(i, k)   = Skv_in(i, k);
        });

    Kokkos::parallel_for(
        "pdf_closure::scalar_defaults", ngrdcol * nz * sclr_dim, KOKKOS_LAMBDA(const int cell) {
          const int i        = cell % ngrdcol;
          const int k        = (cell / ngrdcol) % nz;
          const int sclr     = cell / (ngrdcol * nz);
          Sksclr(i, k, sclr) = Sksclr_in(i, k, sclr);
          if (l_scalar_calc) {
            alpha_sclr(i, k, sclr) = one_half;
          }
        });
  }

  // To avoid recomputing.
  Kokkos::parallel_for(
      "pdf_closure::sqrt_wp2", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = (cell / ngrdcol);
        // To avoid recomputing
        sqrt_wp2(i, k) = Kokkos::sqrt(wp2(i, k));
      });

  // Select the two-component PDF closure for w, rt, theta-l and passive
  // scalars. Each lower call returns component means, variances and weight.
  if (iiPDF_type == iiPDF_ADG1) {
    ADG1_pdf_driver(workspace, nz, ngrdcol, sclr_dim, sclr_tol, // In
                    wm, rtm, thlm, um, vm,                      // In
                    wp2, rtp2, thlp2, up2, vp2,                 // In
                    Skw, wprtp, wpthlp, upwp, vpwp, sqrt_wp2,   // In
                    sigma_sqd_w, Kokkos::subview(clubb_params, Kokkos::ALL, ibeta),
                    mixt_frac_max_mag,                                // In
                    sclrm, sclrp2, wpsclrp, l_scalar_calc,            // In
                    err_info,                                         // In/Out
                    pdf_params.w_1, pdf_params.w_2,                   // Out
                    pdf_params.rt_1, pdf_params.rt_2,                 // Out
                    pdf_params.thl_1, pdf_params.thl_2,               // Out
                    u_1, u_2, v_1, v_2,                               // Out
                    pdf_params.varnce_w_1, pdf_params.varnce_w_2,     // Out
                    pdf_params.varnce_rt_1, pdf_params.varnce_rt_2,   // Out
                    pdf_params.varnce_thl_1, pdf_params.varnce_thl_2, // Out
                    varnce_u_1, varnce_u_2,                           // Out
                    varnce_v_1, varnce_v_2,                           // Out
                    pdf_params.mixt_frac,                             // Out
                    pdf_params.alpha_rt, pdf_params.alpha_thl,        // Out
                    alpha_u, alpha_v,                                 // Out
                    sclr1, sclr2, varnce_sclr1,                       // Out
                    varnce_sclr2, alpha_sclr);                        // Out
  } else if (iiPDF_type == iiPDF_ADG2) {
    ADG2_pdf_driver(workspace, nz, ngrdcol, sclr_dim, sclr_tol, // In
                    wm, rtm, thlm, wp2, rtp2, thlp2,            // In
                    Skw, wprtp, wpthlp, sqrt_wp2,
                    Kokkos::subview(clubb_params, Kokkos::ALL, ibeta), // In
                    sclrm, sclrp2, wpsclrp, l_scalar_calc,             // In
                    err_info,                                          // In/Out
                    pdf_params.w_1, pdf_params.w_2,                    // Out
                    pdf_params.rt_1, pdf_params.rt_2,                  // Out
                    pdf_params.thl_1, pdf_params.thl_2,                // Out
                    pdf_params.varnce_w_1, pdf_params.varnce_w_2,      // Out
                    pdf_params.varnce_rt_1, pdf_params.varnce_rt_2,    // Out
                    pdf_params.varnce_thl_1, pdf_params.varnce_thl_2,  // Out
                    pdf_params.mixt_frac,                              // Out
                    pdf_params.alpha_rt, pdf_params.alpha_thl,         // Out
                    sigma_sqd_w, sclr1, sclr2,                         // Out
                    varnce_sclr1, varnce_sclr2, alpha_sclr);           // Out
  } else if (iiPDF_type == iiPDF_3D_Luhar) {
    for (int i = 0; i < ngrdcol; ++i) {
      Luhar_3D_pdf_driver(
          workspace, nz, Kokkos::subview(wm, i, Kokkos::ALL), Kokkos::subview(rtm, i, Kokkos::ALL),
          Kokkos::subview(thlm, i, Kokkos::ALL), Kokkos::subview(wp2, i, Kokkos::ALL),
          Kokkos::subview(rtp2, i, Kokkos::ALL), Kokkos::subview(thlp2, i, Kokkos::ALL), // In
          Kokkos::subview(Skw, i, Kokkos::ALL), Kokkos::subview(Skrt, i, Kokkos::ALL),
          Kokkos::subview(Skthl, i, Kokkos::ALL), Kokkos::subview(wprtp, i, Kokkos::ALL),
          Kokkos::subview(wpthlp, i, Kokkos::ALL), // In
          Kokkos::subview(pdf_params.w_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.w_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.rt_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.rt_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.thl_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.thl_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.varnce_w_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.varnce_w_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.varnce_rt_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.varnce_rt_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.varnce_thl_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.varnce_thl_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.mixt_frac, i, Kokkos::ALL));   // Out
    }
  } else if (iiPDF_type == iiPDF_new) {
    new_pdf_driver(workspace, nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2, // In
                   Skw,                                                     // In
                   wprtp, wpthlp, rtpthlp,                                  // In
                   clubb_params,                                            // In
                   Skrt, Skthl,                                             // In/Out
                   pdf_params.w_1, pdf_params.w_2,                          // Out
                   pdf_params.rt_1, pdf_params.rt_2,                        // Out
                   pdf_params.thl_1, pdf_params.thl_2,                      // Out
                   pdf_params.varnce_w_1, pdf_params.varnce_w_2,            // Out
                   pdf_params.varnce_rt_1, pdf_params.varnce_rt_2,          // Out
                   pdf_params.varnce_thl_1, pdf_params.varnce_thl_2,        // Out
                   pdf_params.mixt_frac,                                    // Out
                   pdf_implicit_coefs_terms);                               // Out
  } else if (iiPDF_type == iiPDF_TSDADG) {
    for (int i = 0; i < ngrdcol; ++i) {
      tsdadg_pdf_driver(
          workspace, nz, Kokkos::subview(wm, i, Kokkos::ALL), Kokkos::subview(rtm, i, Kokkos::ALL),
          Kokkos::subview(thlm, i, Kokkos::ALL), Kokkos::subview(wp2, i, Kokkos::ALL),
          Kokkos::subview(rtp2, i, Kokkos::ALL), Kokkos::subview(thlp2, i, Kokkos::ALL), // In
          Kokkos::subview(Skw, i, Kokkos::ALL), Kokkos::subview(Skrt, i, Kokkos::ALL),
          Kokkos::subview(Skthl, i, Kokkos::ALL), Kokkos::subview(wprtp, i, Kokkos::ALL),
          Kokkos::subview(wpthlp, i, Kokkos::ALL), // In
          Kokkos::subview(pdf_params.w_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.w_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.rt_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.rt_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.thl_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.thl_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.varnce_w_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.varnce_w_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.varnce_rt_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.varnce_rt_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.varnce_thl_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.varnce_thl_2, i, Kokkos::ALL), // Out
          Kokkos::subview(pdf_params.mixt_frac, i, Kokkos::ALL));   // Out
    }
  } else if (iiPDF_type == iiPDF_LY93) {
    for (int i = 0; i < ngrdcol; ++i) {
      LY93_driver(workspace, nz, Kokkos::subview(wm, i, Kokkos::ALL),
                  Kokkos::subview(rtm, i, Kokkos::ALL), Kokkos::subview(thlm, i, Kokkos::ALL),
                  Kokkos::subview(wp2, i, Kokkos::ALL), Kokkos::subview(rtp2, i, Kokkos::ALL), // In
                  Kokkos::subview(thlp2, i, Kokkos::ALL), Kokkos::subview(Skw, i, Kokkos::ALL),
                  Kokkos::subview(Skrt, i, Kokkos::ALL),
                  Kokkos::subview(Skthl, i, Kokkos::ALL), // In
                  Kokkos::subview(pdf_params.w_1, i, Kokkos::ALL),
                  Kokkos::subview(pdf_params.w_2, i, Kokkos::ALL), // Out
                  Kokkos::subview(pdf_params.rt_1, i, Kokkos::ALL),
                  Kokkos::subview(pdf_params.rt_2, i, Kokkos::ALL), // Out
                  Kokkos::subview(pdf_params.thl_1, i, Kokkos::ALL),
                  Kokkos::subview(pdf_params.thl_2, i, Kokkos::ALL), // Out
                  Kokkos::subview(pdf_params.varnce_w_1, i, Kokkos::ALL),
                  Kokkos::subview(pdf_params.varnce_w_2, i, Kokkos::ALL), // Out
                  Kokkos::subview(pdf_params.varnce_rt_1, i, Kokkos::ALL),
                  Kokkos::subview(pdf_params.varnce_rt_2, i, Kokkos::ALL), // Out
                  Kokkos::subview(pdf_params.varnce_thl_1, i, Kokkos::ALL),
                  Kokkos::subview(pdf_params.varnce_thl_2, i, Kokkos::ALL), // Out
                  Kokkos::subview(pdf_params.mixt_frac, i, Kokkos::ALL));   // Out
    }
  } else if (iiPDF_type == iiPDF_new_hybrid) {
    new_hybrid_pdf_driver(
        workspace, nz, ngrdcol, sclr_dim,                                          // In
        wm, rtm, thlm, um, vm,                                                     // In
        wp2, rtp2, thlp2, up2, vp2,                                                // In
        Skw, wprtp, wpthlp, upwp, vpwp,                                            // In
        sclrm, sclrp2, wpsclrp,                                                    // In
        clubb_params,                                                              // In
        Kokkos::subview(clubb_params, Kokkos::ALL, islope_coef_spread_DG_means_w), // In
        Kokkos::subview(clubb_params, Kokkos::ALL, ipdf_component_stdev_factor_w), // In
        Skrt, Skthl, Sku, Skv, Sksclr,                                             // I/O
        pdf_params.w_1, pdf_params.w_2,                                            // Out
        pdf_params.rt_1, pdf_params.rt_2,                                          // Out
        pdf_params.thl_1, pdf_params.thl_2,                                        // Out
        u_1, u_2, v_1, v_2,                                                        // Out
        pdf_params.varnce_w_1,                                                     // Out
        pdf_params.varnce_w_2,                                                     // Out
        pdf_params.varnce_rt_1,                                                    // Out
        pdf_params.varnce_rt_2,                                                    // Out
        pdf_params.varnce_thl_1,                                                   // Out
        pdf_params.varnce_thl_2,                                                   // Out
        varnce_u_1, varnce_u_2,                                                    // Out
        varnce_v_1, varnce_v_2,                                                    // Out
        sclr1, sclr2,                                                              // Out
        varnce_sclr1, varnce_sclr2,                                                // Out
        pdf_params.mixt_frac, sigma_sqd_w,                                         // Out
        pdf_implicit_coefs_terms);                                                 // Out
  }

  // Calculate the PDF component correlations of rt and thl.
  calc_comp_corrs_binormal(nz, ngrdcol,                                         // In
                           rtpthlp, rtm, thlm,                                  // In
                           pdf_params.rt_1, pdf_params.rt_2,                    // In
                           pdf_params.thl_1, pdf_params.thl_2,                  // In
                           pdf_params.varnce_rt_1, pdf_params.varnce_rt_2,      // In
                           pdf_params.varnce_thl_1, pdf_params.varnce_thl_2,    // In
                           pdf_params.mixt_frac,                                // In
                           pdf_params.corr_rt_thl_1, pdf_params.corr_rt_thl_2); // Out

  if (iiPDF_type == iiPDF_ADG1 || iiPDF_type == iiPDF_ADG2 || iiPDF_type == iiPDF_new_hybrid) {
    // These PDFs set component correlations involving w (including u/v) to zero.
    Kokkos::parallel_for(
        "pdf_closure::zero_w_correlations", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = (cell / ngrdcol);
          // These PDF types define corr_w_rt_1, corr_w_rt_2, corr_w_thl_1, and
          // corr_w_thl_2 to all have a value of 0, so skip the calculation.
          // The values of corr_u_w_1, corr_u_w_2, corr_v_w_1, and corr_v_w_2 are
          // all defined to be 0, as well.
          pdf_params_corr_w_rt_1(i, k)  = zero;
          pdf_params_corr_w_rt_2(i, k)  = zero;
          pdf_params_corr_w_thl_1(i, k) = zero;
          pdf_params_corr_w_thl_2(i, k) = zero;
          corr_u_w_1(i, k)              = zero;
          corr_u_w_2(i, k)              = zero;
          corr_v_w_1(i, k)              = zero;
          corr_v_w_2(i, k)              = zero;
        });

  } else {
    // Calculate component correlations of w with rt, then thl.
    calc_comp_corrs_binormal(nz, ngrdcol,                                       // In
                             wprtp, wm, rtm,                                    // In
                             pdf_params.w_1, pdf_params.w_2,                    // In
                             pdf_params.rt_1, pdf_params.rt_2,                  // In
                             pdf_params.varnce_w_1, pdf_params.varnce_w_2,      // In
                             pdf_params.varnce_rt_1, pdf_params.varnce_rt_2,    // In
                             pdf_params.mixt_frac,                              // In
                             pdf_params.corr_w_rt_1, pdf_params.corr_w_rt_2);   // Out
    calc_comp_corrs_binormal(nz, ngrdcol,                                       // In
                             wpthlp, wm, thlm,                                  // In
                             pdf_params.w_1, pdf_params.w_2,                    // In
                             pdf_params.thl_1, pdf_params.thl_2,                // In
                             pdf_params.varnce_w_1, pdf_params.varnce_w_2,      // In
                             pdf_params.varnce_thl_1, pdf_params.varnce_thl_2,  // In
                             pdf_params.mixt_frac,                              // In
                             pdf_params.corr_w_thl_1, pdf_params.corr_w_thl_2); // Out
  }

  if (l_scalar_calc) {
    // Calculate the PDF component correlations of each scalar with thl and rt.
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      calc_comp_corrs_binormal(
          nz, ngrdcol, // In
          Kokkos::subview(sclrpthlp, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr), thlm, // In
          Kokkos::subview(sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(sclr2, Kokkos::ALL, Kokkos::ALL, sclr), // In
          pdf_params.thl_1, pdf_params.thl_2,                     // In
          Kokkos::subview(varnce_sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(varnce_sclr2, Kokkos::ALL, Kokkos::ALL, sclr), // In
          pdf_params.varnce_thl_1, pdf_params.varnce_thl_2,              // In
          pdf_params.mixt_frac,                                          // In
          Kokkos::subview(corr_sclr_thl_1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(corr_sclr_thl_2, Kokkos::ALL, Kokkos::ALL, sclr)); // Out
    }
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      calc_comp_corrs_binormal(
          nz, ngrdcol, // In
          Kokkos::subview(sclrprtp, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr), rtm, // In
          Kokkos::subview(sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(sclr2, Kokkos::ALL, Kokkos::ALL, sclr), // In
          pdf_params.rt_1, pdf_params.rt_2,                       // In
          Kokkos::subview(varnce_sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(varnce_sclr2, Kokkos::ALL, Kokkos::ALL, sclr), // In
          pdf_params.varnce_rt_1, pdf_params.varnce_rt_2,                // In
          pdf_params.mixt_frac,                                          // In
          Kokkos::subview(corr_sclr_rt_1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(corr_sclr_rt_2, Kokkos::ALL, Kokkos::ALL, sclr)); // Out
    }
    if (iiPDF_type == iiPDF_ADG1 || iiPDF_type == iiPDF_ADG2 || iiPDF_type == iiPDF_new_hybrid) {
      // Component correlations involving w vanish for these PDF choices.
      Kokkos::parallel_for(
          "pdf_closure::zero_scalar_w_correlations", ngrdcol * nz * sclr_dim,
          KOKKOS_LAMBDA(const int cell) {
            const int i    = cell % ngrdcol;
            const int k    = (cell / ngrdcol) % nz;
            const int sclr = cell / (ngrdcol * nz);
            // These PDF types define all PDF component correlations involving w
            // to have a value of 0, so skip the calculation.
            corr_w_sclr_1(i, k, sclr) = zero;
            corr_w_sclr_2(i, k, sclr) = zero;
          });

    } else {
      for (int sclr = 0; sclr < sclr_dim; ++sclr) {
        calc_comp_corrs_binormal(
            nz, ngrdcol, // In
            Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr), wm,
            Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr), // In
            pdf_params.w_1, pdf_params.w_2,                         // In
            Kokkos::subview(sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
            Kokkos::subview(sclr2, Kokkos::ALL, Kokkos::ALL, sclr), // In
            pdf_params.varnce_w_1, pdf_params.varnce_w_2,           // In
            Kokkos::subview(varnce_sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
            Kokkos::subview(varnce_sclr2, Kokkos::ALL, Kokkos::ALL, sclr), // In
            pdf_params.mixt_frac,                                          // In
            Kokkos::subview(corr_w_sclr_1, Kokkos::ALL, Kokkos::ALL, sclr),
            Kokkos::subview(corr_w_sclr_2, Kokkos::ALL, Kokkos::ALL, sclr)); // Out
      }
    }
  }

  // Compute higher-order moments (these are interactive).
  calc_wp2xp_pdf(nz, ngrdcol, wm, rtm, pdf_params.w_1, pdf_params.w_2, pdf_params.rt_1,
                 pdf_params.rt_2, pdf_params.varnce_w_1, pdf_params.varnce_w_2,
                 pdf_params.varnce_rt_1, pdf_params.varnce_rt_2, pdf_params.corr_w_rt_1,
                 pdf_params.corr_w_rt_2, pdf_params.mixt_frac, wp2rtp);
  calc_wp2xp_pdf(nz, ngrdcol, wm, thlm, pdf_params.w_1, pdf_params.w_2, pdf_params.thl_1,
                 pdf_params.thl_2, pdf_params.varnce_w_1, pdf_params.varnce_w_2,
                 pdf_params.varnce_thl_1, pdf_params.varnce_thl_2, pdf_params.corr_w_thl_1,
                 pdf_params.corr_w_thl_2, pdf_params.mixt_frac, wp2thlp);
  calc_wp2xp_pdf(nz, ngrdcol, wm, um, pdf_params.w_1, pdf_params.w_2, u_1, u_2,
                 pdf_params.varnce_w_1, pdf_params.varnce_w_2, varnce_u_1, varnce_u_2, corr_u_w_1,
                 corr_u_w_2, pdf_params.mixt_frac, wp2up);

  // Compute higher-order moments (these may be interactive).
  calc_wpxp2_pdf(nz, ngrdcol, wm, um, pdf_params.w_1, pdf_params.w_2, u_1, u_2,
                 pdf_params.varnce_w_1, pdf_params.varnce_w_2, varnce_u_1, varnce_u_2, corr_u_w_1,
                 corr_u_w_2, pdf_params.mixt_frac, wpup2);
  calc_wpxp2_pdf(nz, ngrdcol, wm, vm, pdf_params.w_1, pdf_params.w_2, v_1, v_2,
                 pdf_params.varnce_w_1, pdf_params.varnce_w_2, varnce_v_1, varnce_v_2, corr_v_w_1,
                 corr_v_w_2, pdf_params.mixt_frac, wpvp2);
  calc_wp2xp2_pdf(nz, ngrdcol, wm, um, pdf_params.w_1, pdf_params.w_2, u_1, u_2,
                  pdf_params.varnce_w_1, pdf_params.varnce_w_2, varnce_u_1, varnce_u_2, corr_u_w_1,
                  corr_u_w_2, pdf_params.mixt_frac, wp2up2);
  calc_wp2xp2_pdf(nz, ngrdcol, wm, vm, pdf_params.w_1, pdf_params.w_2, v_1, v_2,
                  pdf_params.varnce_w_1, pdf_params.varnce_w_2, varnce_v_1, varnce_v_2, corr_v_w_1,
                  corr_v_w_2, pdf_params.mixt_frac, wp2vp2);
  calc_wp4_pdf(nz, ngrdcol, wm, pdf_params.w_1, pdf_params.w_2, pdf_params.varnce_w_1,
               pdf_params.varnce_w_2, pdf_params.mixt_frac, wp4);

  if (l_explicit_turbulent_adv_xpyp || var_on_stats_list(stats, "wprtp2")) {
    calc_wpxp2_pdf(nz, ngrdcol, wm, rtm, pdf_params.w_1, pdf_params.w_2, pdf_params.rt_1,
                   pdf_params.rt_2, pdf_params.varnce_w_1, pdf_params.varnce_w_2,
                   pdf_params.varnce_rt_1, pdf_params.varnce_rt_2, pdf_params.corr_w_rt_1,
                   pdf_params.corr_w_rt_2, pdf_params.mixt_frac, wprtp2);
  }
  if (l_explicit_turbulent_adv_xpyp || var_on_stats_list(stats, "wpthlp2")) {
    calc_wpxp2_pdf(nz, ngrdcol, wm, thlm, pdf_params.w_1, pdf_params.w_2, pdf_params.thl_1,
                   pdf_params.thl_2, pdf_params.varnce_w_1, pdf_params.varnce_w_2,
                   pdf_params.varnce_thl_1, pdf_params.varnce_thl_2, pdf_params.corr_w_thl_1,
                   pdf_params.corr_w_thl_2, pdf_params.mixt_frac, wpthlp2);
  }
  if (l_explicit_turbulent_adv_xpyp || var_on_stats_list(stats, "wprtpthlp")) {
    calc_wpxpyp_pdf(nz, ngrdcol, wm, rtm, thlm, pdf_params.w_1, pdf_params.w_2, pdf_params.rt_1,
                    pdf_params.rt_2, pdf_params.thl_1, pdf_params.thl_2, pdf_params.varnce_w_1,
                    pdf_params.varnce_w_2, pdf_params.varnce_rt_1, pdf_params.varnce_rt_2,
                    pdf_params.varnce_thl_1, pdf_params.varnce_thl_2, pdf_params.corr_w_rt_1,
                    pdf_params.corr_w_rt_2, pdf_params.corr_w_thl_1, pdf_params.corr_w_thl_2,
                    pdf_params.corr_rt_thl_1, pdf_params.corr_rt_thl_2, pdf_params.mixt_frac,
                    wprtpthlp);
  }

  // Scalar addition to higher-order moments.
  if (l_scalar_calc) {
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      calc_wp2xp_pdf(
          nz, ngrdcol, wm, Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.w_1,
          pdf_params.w_2, Kokkos::subview(sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(sclr2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.varnce_w_1,
          pdf_params.varnce_w_2, Kokkos::subview(varnce_sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(varnce_sclr2, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(corr_w_sclr_1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(corr_w_sclr_2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.mixt_frac,
          Kokkos::subview(wp2sclrp, Kokkos::ALL, Kokkos::ALL, sclr));
    }
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      calc_wpxp2_pdf(
          nz, ngrdcol, wm, Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.w_1,
          pdf_params.w_2, Kokkos::subview(sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(sclr2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.varnce_w_1,
          pdf_params.varnce_w_2, Kokkos::subview(varnce_sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(varnce_sclr2, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(corr_w_sclr_1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(corr_w_sclr_2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.mixt_frac,
          Kokkos::subview(wpsclrp2, Kokkos::ALL, Kokkos::ALL, sclr));
    }
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      calc_wpxpyp_pdf(
          nz, ngrdcol, wm, Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr), rtm,
          pdf_params.w_1, pdf_params.w_2, Kokkos::subview(sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(sclr2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.rt_1, pdf_params.rt_2,
          pdf_params.varnce_w_1, pdf_params.varnce_w_2,
          Kokkos::subview(varnce_sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(varnce_sclr2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.varnce_rt_1,
          pdf_params.varnce_rt_2, Kokkos::subview(corr_w_sclr_1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(corr_w_sclr_2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.corr_w_rt_1,
          pdf_params.corr_w_rt_2, Kokkos::subview(corr_sclr_rt_1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(corr_sclr_rt_2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.mixt_frac,
          Kokkos::subview(wpsclrprtp, Kokkos::ALL, Kokkos::ALL, sclr));
    }
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      calc_wpxpyp_pdf(
          nz, ngrdcol, wm, Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr), thlm,
          pdf_params.w_1, pdf_params.w_2, Kokkos::subview(sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(sclr2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.thl_1,
          pdf_params.thl_2, pdf_params.varnce_w_1, pdf_params.varnce_w_2,
          Kokkos::subview(varnce_sclr1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(varnce_sclr2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.varnce_thl_1,
          pdf_params.varnce_thl_2, Kokkos::subview(corr_w_sclr_1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(corr_w_sclr_2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.corr_w_thl_1,
          pdf_params.corr_w_thl_2, Kokkos::subview(corr_sclr_thl_1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(corr_sclr_thl_2, Kokkos::ALL, Kokkos::ALL, sclr), pdf_params.mixt_frac,
          Kokkos::subview(wpsclrpthlp, Kokkos::ALL, Kokkos::ALL, sclr));
    }
  }

  // Compute higher-order moments that include theta_v.
  // First compute component liquid-water temperatures (Sommeria & Deardorff
  // 1977, Eq.3); "1" and "2" denote the two Gaussians.
  Kokkos::parallel_for(
      "pdf_closure::component_temperatures", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = (cell / ngrdcol);
        // Compute higher order moments that include theta_v.
        // First compute some preliminary quantities.
        // "1" denotes first Gaussian; "2" denotes 2nd Gaussian
        // liq water temp (Sommeria & Deardorff 1977 (SD), eqn. 3)
        tl1(i, k) = pdf_params_thl_1(i, k) * exner(i, k);
        tl2(i, k) = pdf_params_thl_2(i, k) * exner(i, k);
      });

  sat_mixrat_liq_api(nz, ngrdcol, gr, p_in_Pa, tl1, saturation_formula, pdf_params.rsatl_1);
  sat_mixrat_liq_api(nz, ngrdcol, gr, p_in_Pa, tl2, saturation_formula, pdf_params.rsatl_2);

  const bool l_calc_ice_supersat_frac = true;
  transform_pdf_chi_eta_component(
      workspace, nz, ngrdcol, tl1, pdf_params.rsatl_1, pdf_params.rt_1, exner, // In
      pdf_params.varnce_thl_1, pdf_params.varnce_rt_1,                         // In
      pdf_params.corr_rt_thl_1, pdf_params.chi_1, // corr_rt_thl In; chi Out (native declaration)
      pdf_params.crt_1, pdf_params.cthl_1,        // Out
      pdf_params.stdev_chi_1, pdf_params.stdev_eta_1, // Out
      pdf_params.covar_chi_eta_1,                     // Out
      pdf_params.corr_chi_eta_1);                     // Out
  // Calculate cloud fraction component for PDF1.
  calc_liquid_cloud_frac_component(nz, ngrdcol, pdf_params.chi_1, pdf_params.stdev_chi_1, // In
                                   pdf_params.cloud_frac_1, pdf_params.rc_1);             // Out
  if (l_calc_ice_supersat_frac) {
    calc_ice_cloud_frac_component(nz, ngrdcol, pdf_params.chi_1, pdf_params.stdev_chi_1,
                                  pdf_params.rc_1, pdf_params.cloud_frac_1, p_in_Pa, tl1,
                                  pdf_params.rsatl_1, pdf_params.crt_1, saturation_formula,
                                  pdf_params.ice_supersat_frac_1, rc_1_ice);
  }
  transform_pdf_chi_eta_component(
      workspace, nz, ngrdcol, tl2, pdf_params.rsatl_2, pdf_params.rt_2, exner, // In
      pdf_params.varnce_thl_2, pdf_params.varnce_rt_2,                         // In
      pdf_params.corr_rt_thl_2, pdf_params.chi_2, // corr_rt_thl In; chi Out (native declaration)
      pdf_params.crt_2, pdf_params.cthl_2,        // Out
      pdf_params.stdev_chi_2, pdf_params.stdev_eta_2, // Out
      pdf_params.covar_chi_eta_2,                     // Out
      pdf_params.corr_chi_eta_2);                     // Out
  // Calculate cloud fraction component for PDF2.
  calc_liquid_cloud_frac_component(nz, ngrdcol, pdf_params.chi_2, pdf_params.stdev_chi_2, // In
                                   pdf_params.cloud_frac_2, pdf_params.rc_2);             // Out
  if (l_calc_ice_supersat_frac) {
    calc_ice_cloud_frac_component(nz, ngrdcol, pdf_params.chi_2, pdf_params.stdev_chi_2,
                                  pdf_params.rc_2, pdf_params.cloud_frac_2, p_in_Pa, tl2,
                                  pdf_params.rsatl_2, pdf_params.crt_2, saturation_formula,
                                  pdf_params.ice_supersat_frac_2, rc_2_ice);
    // Compute ice supersaturation fraction.
    Kokkos::parallel_for(
        "pdf_closure::ice_mixture", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = (cell / ngrdcol);
          // Compute ice cloud fraction, ice_supersat_frac
          ice_supersat_frac(i, k) =
              pdf_params_mixt_frac(i, k) * pdf_params_ice_supersat_frac_1(i, k) +
              (one - pdf_params_mixt_frac(i, k)) * pdf_params_ice_supersat_frac_2(i, k);
        });

  } else {
    // Ice supersaturation fraction would otherwise contain garbage values.
    Kokkos::deep_copy(ice_supersat_frac, zero);
    if (clubb_at_least_debug_level_api(1))
      std::cerr
          << "Warning: ice_supersat_frac has garbage values if do_liquid_only_in_clubb = .false.\n";
  }

  // Compute cloud fraction and mean cloud water mixing ratio.
  // Reference: https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:anl_int_cloud_terms
  Kokkos::parallel_for(
      "pdf_closure::cloud_mixture", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = (cell / ngrdcol);
        // Compute cloud fraction and mean cloud water mixing ratio.
        // Reference:
        // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:anl_int_cloud_terms
        // Native rounds the first component product and fuses the second;
        // implicit CUDA contraction selects the opposite product (BFB-011).
        cloud_frac(i, k) =
            Kokkos::fma(one - pdf_params_mixt_frac(i, k), pdf_params_cloud_frac_2(i, k),
                        pdf_params_mixt_frac(i, k) * pdf_params_cloud_frac_1(i, k));
        rcm(i, k) = Kokkos::fma(one - pdf_params_mixt_frac(i, k), pdf_params_rc_2(i, k),
                                pdf_params_mixt_frac(i, k) * pdf_params_rc_1(i, k));
        rcm(i, k) = ((zero_threshold) > (rcm(i, k)) ? (zero_threshold) : (rcm(i, k)));
      });

  if (iiPDF_type == iiPDF_ADG1 || iiPDF_type == iiPDF_ADG2 || iiPDF_type == iiPDF_new_hybrid) {
    // corr_w_rt and corr_w_thl are zero, so corr_w_chi and corr_w_eta are zero too.
    Kokkos::parallel_for(
        "pdf_closure::zero_chi_eta_correlations", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = (cell / ngrdcol);
          // corr_w_rt and corr_w_thl are zero for these pdf types so
          // corr_w_chi and corr_w_eta are zero as well
          pdf_params_corr_w_chi_1(i, k) = zero;
          pdf_params_corr_w_chi_2(i, k) = zero;
          pdf_params_corr_w_eta_1(i, k) = zero;
          pdf_params_corr_w_eta_2(i, k) = zero;
        });

  } else {
    // Independent correlations share one cell kernel, avoiding array-call scratch/copies.
    const auto pdf_params_stdev_eta_1 = pdf_params.stdev_eta_1;
    const auto pdf_params_stdev_eta_2 = pdf_params.stdev_eta_2;
    Kokkos::parallel_for(
        "pdf_closure::component_correlations", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i                    = cell % ngrdcol;
          const int k                    = cell / ngrdcol;
          const real_t sqrt_varnce_rt_1  = Kokkos::sqrt(pdf_params_varnce_rt_1(i, k));
          const real_t sqrt_varnce_rt_2  = Kokkos::sqrt(pdf_params_varnce_rt_2(i, k));
          const real_t sqrt_varnce_thl_1 = Kokkos::sqrt(pdf_params_varnce_thl_1(i, k));
          const real_t sqrt_varnce_thl_2 = Kokkos::sqrt(pdf_params_varnce_thl_2(i, k));
          // Correlation of w and chi for each component.
          pdf_params_corr_w_chi_1(i, k) =
              calc_corr_chi_x(pdf_params_crt_1(i, k), pdf_params_cthl_1(i, k), sqrt_varnce_rt_1,
                              sqrt_varnce_thl_1, pdf_params_stdev_chi_1(i, k),
                              pdf_params_corr_w_rt_1(i, k), pdf_params_corr_w_thl_1(i, k));
          pdf_params_corr_w_chi_2(i, k) =
              calc_corr_chi_x(pdf_params_crt_2(i, k), pdf_params_cthl_2(i, k), sqrt_varnce_rt_2,
                              sqrt_varnce_thl_2, pdf_params_stdev_chi_2(i, k),
                              pdf_params_corr_w_rt_2(i, k), pdf_params_corr_w_thl_2(i, k));
          // Correlation of w and eta for each component.
          pdf_params_corr_w_eta_1(i, k) =
              calc_corr_eta_x(pdf_params_crt_1(i, k), pdf_params_cthl_1(i, k), sqrt_varnce_rt_1,
                              sqrt_varnce_thl_1, pdf_params_stdev_eta_1(i, k),
                              pdf_params_corr_w_rt_1(i, k), pdf_params_corr_w_thl_1(i, k));
          pdf_params_corr_w_eta_2(i, k) =
              calc_corr_eta_x(pdf_params_crt_2(i, k), pdf_params_cthl_2(i, k), sqrt_varnce_rt_2,
                              sqrt_varnce_thl_2, pdf_params_stdev_eta_2(i, k),
                              pdf_params_corr_w_rt_2(i, k), pdf_params_corr_w_thl_2(i, k));
        });
  }

  // Moments depending on theta_v use its approximated, linearized equation:
  // theta_v = theta_l + ep1*thv_ds*r_t + [Lv/(Cp*exner) - ep2*thv_ds]*r_c.
  // Therefore the perturbation has the same coefficients on th_l', r_t', r_c'.
  // thv_ds is the reference value approximating theta_l.
  // Reference: https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:anl_int_buoy_terms
  // Calculate component1 and component2 contributions to liquid-water moments.
  calc_xprcp_component(nz, ngrdcol,                                          // In
                       wm, rtm, thlm, um, vm, rcm,                           // In
                       pdf_params.w_1, pdf_params.rt_1,                      // In
                       pdf_params.thl_1, u_1, v_1,                           // In
                       pdf_params.varnce_w_1, pdf_params.chi_1,              // In
                       pdf_params.stdev_chi_1, pdf_params.stdev_eta_1,       // In
                       pdf_params.corr_w_chi_1, pdf_params.corr_chi_eta_1,   // In
                       pdf_params.crt_1, pdf_params.cthl_1,                  // In
                       pdf_params.rc_1, pdf_params.cloud_frac_1, iiPDF_type, // In
                       wprcp_contrib_comp_1, wp2rcp_contrib_comp_1,          // Out
                       rtprcp_contrib_comp_1, thlprcp_contrib_comp_1,        // Out
                       uprcp_contrib_comp_1, vprcp_contrib_comp_1);          // Out
  calc_xprcp_component(nz, ngrdcol,                                          // In
                       wm, rtm, thlm, um, vm, rcm,                           // In
                       pdf_params.w_2, pdf_params.rt_2,                      // In
                       pdf_params.thl_2, u_2, v_2,                           // In
                       pdf_params.varnce_w_2, pdf_params.chi_2,              // In
                       pdf_params.stdev_chi_2, pdf_params.stdev_eta_2,       // In
                       pdf_params.corr_w_chi_2, pdf_params.corr_chi_eta_2,   // In
                       pdf_params.crt_2, pdf_params.cthl_2,                  // In
                       pdf_params.rc_2, pdf_params.cloud_frac_2, iiPDF_type, // In
                       wprcp_contrib_comp_2, wp2rcp_contrib_comp_2,          // Out
                       rtprcp_contrib_comp_2, thlprcp_contrib_comp_2,        // Out
                       uprcp_contrib_comp_2, vprcp_contrib_comp_2);          // Out
  Kokkos::parallel_for(
      "pdf_closure::water_moment_mixture", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = (cell / ngrdcol);
        // Calculate rc_coef, which is the coefficient on <x'rc'> in the <x'thv'> equation.
        rc_coef(i, k) = Lv / (exner(i, k) * Cp) - ep2 * thv_ds(i, k);
        // Calculate <w'rc'>, <w'^2 rc'>, <rt'rc'>, and <thl'rc'>.
        // Native rounds component1, then fuses component2 for each moment.
        // Keep the same contraction as the cloud/water means (BFB-012).
        wprcp(i, k)   = Kokkos::fma(one - pdf_params_mixt_frac(i, k), wprcp_contrib_comp_2(i, k),
                                    pdf_params_mixt_frac(i, k) * wprcp_contrib_comp_1(i, k));
        wp2rcp(i, k)  = Kokkos::fma(one - pdf_params_mixt_frac(i, k), wp2rcp_contrib_comp_2(i, k),
                                    pdf_params_mixt_frac(i, k) * wp2rcp_contrib_comp_1(i, k));
        rtprcp(i, k)  = Kokkos::fma(one - pdf_params_mixt_frac(i, k), rtprcp_contrib_comp_2(i, k),
                                    pdf_params_mixt_frac(i, k) * rtprcp_contrib_comp_1(i, k));
        thlprcp(i, k) = Kokkos::fma(one - pdf_params_mixt_frac(i, k), thlprcp_contrib_comp_2(i, k),
                                    pdf_params_mixt_frac(i, k) * thlprcp_contrib_comp_1(i, k));
        uprcp(i, k)   = Kokkos::fma(one - pdf_params_mixt_frac(i, k), uprcp_contrib_comp_2(i, k),
                                    pdf_params_mixt_frac(i, k) * uprcp_contrib_comp_1(i, k));
        vprcp(i, k)   = Kokkos::fma(one - pdf_params_mixt_frac(i, k), vprcp_contrib_comp_2(i, k),
                                    pdf_params_mixt_frac(i, k) * vprcp_contrib_comp_1(i, k));
      });

  // Calculate w'thv', w'^2thv', rt'thv' and thl'thv'.
  Kokkos::parallel_for(
      "pdf_closure::buoyancy_moments", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = (cell / ngrdcol);
        // Calculate <w'thv'>, <w'^2 thv'>, <rt'thv'>, and <thl'thv'>.
        wpthvp(i, k) =
            wpthlp(i, k) + ep1 * thv_ds(i, k) * wprtp(i, k) + rc_coef(i, k) * wprcp(i, k);
        wp2thvp(i, k) =
            wp2thlp(i, k) + ep1 * thv_ds(i, k) * wp2rtp(i, k) + rc_coef(i, k) * wp2rcp(i, k);
        rtpthvp(i, k) =
            rtpthlp(i, k) + ep1 * thv_ds(i, k) * rtp2(i, k) + rc_coef(i, k) * rtprcp(i, k);
        thlpthvp(i, k) =
            thlp2(i, k) + ep1 * thv_ds(i, k) * rtpthlp(i, k) + rc_coef(i, k) * thlprcp(i, k);
      });

  // Add precipitation loading in the x'thv' equation (native test flag is false).
  if (l_liq_ice_loading_test) {
    const auto l_mix_rat_hm_host = to_host_space(workspace, l_mix_rat_hm);
    for (int hm_idx = 0; hm_idx < hydromet_dim; ++hm_idx) {
      if (l_mix_rat_hm_host(hm_idx)) {
        Kokkos::parallel_for(
            "pdf_closure::precipitation_loading", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
              const int i    = cell % ngrdcol;
              const int k    = (cell / ngrdcol);
              wp2thvp(i, k)  = wp2thvp(i, k) - thv_ds(i, k) * wp2hmp(i, k, hm_idx);
              wpthvp(i, k)   = wpthvp(i, k) - thv_ds(i, k) * wphydrometp(i, k, hm_idx);
              thlpthvp(i, k) = thlpthvp(i, k) - thv_ds(i, k) * thlphmp(i, k, hm_idx);
              rtpthvp(i, k)  = rtpthvp(i, k) - thv_ds(i, k) * rtphmp(i, k, hm_idx);
            });
      }
    }
  }

  // Account for subplume correlation of scalar and theta_v.
  // Larson et al.(2002), "Small-scale...", Eqs.A13/A8 (their scalar is w).
  if (l_scalar_calc) {
    Kokkos::parallel_for(
        "pdf_closure::scalar_buoyancy", ngrdcol * nz * sclr_dim, KOKKOS_LAMBDA(const int cell) {
          const int i    = cell % ngrdcol;
          const int k    = (cell / ngrdcol) % nz;
          const int sclr = cell / (ngrdcol * nz);
          sclrprcp(i, k, sclr) =
              pdf_params_mixt_frac(i, k) *
                  ((sclr1(i, k, sclr) - sclrm(i, k, sclr)) * pdf_params_rc_1(i, k)) +
              (one - pdf_params_mixt_frac(i, k)) *
                  ((sclr2(i, k, sclr) - sclrm(i, k, sclr)) * pdf_params_rc_2(i, k)) +
              pdf_params_mixt_frac(i, k) * corr_sclr_rt_1(i, k, sclr) * pdf_params_crt_1(i, k) *
                  Kokkos::sqrt(varnce_sclr1(i, k, sclr) * pdf_params_varnce_rt_1(i, k)) *
                  pdf_params_cloud_frac_1(i, k) +
              (one - pdf_params_mixt_frac(i, k)) * corr_sclr_rt_2(i, k, sclr) *
                  pdf_params_crt_2(i, k) *
                  Kokkos::sqrt(varnce_sclr2(i, k, sclr) * pdf_params_varnce_rt_2(i, k)) *
                  pdf_params_cloud_frac_2(i, k) -
              pdf_params_mixt_frac(i, k) * corr_sclr_thl_1(i, k, sclr) * pdf_params_cthl_1(i, k) *
                  Kokkos::sqrt(varnce_sclr1(i, k, sclr) * pdf_params_varnce_thl_1(i, k)) *
                  pdf_params_cloud_frac_1(i, k) -
              (one - pdf_params_mixt_frac(i, k)) * corr_sclr_thl_2(i, k, sclr) *
                  pdf_params_cthl_2(i, k) *
                  Kokkos::sqrt(varnce_sclr2(i, k, sclr) * pdf_params_varnce_thl_2(i, k)) *
                  pdf_params_cloud_frac_2(i, k);
          sclrpthvp(i, k, sclr) = sclrpthlp(i, k, sclr) +
                                  ep1 * thv_ds(i, k) * sclrprtp(i, k, sclr) +
                                  rc_coef(i, k) * sclrprcp(i, k, sclr);
        });
  }

  // Liquid-water variance is needed only for statistics, except in CAM.
#ifndef CLUBB_CAM
  if (var_on_stats_list(stats, "rcp2")) {
#endif
    Kokkos::parallel_for(
        "pdf_closure::liquid_water_variance", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = (cell / ngrdcol);
          rcp2(i, k)  = pdf_params_mixt_frac(i, k) *
                            (pdf_params_chi_1(i, k) * pdf_params_rc_1(i, k) +
                             pdf_params_cloud_frac_1(i, k) *
                                 (pdf_params_stdev_chi_1(i, k) * pdf_params_stdev_chi_1(i, k))) +
                        (one - pdf_params_mixt_frac(i, k)) *
                            (pdf_params_chi_2(i, k) * pdf_params_rc_2(i, k) +
                             pdf_params_cloud_frac_2(i, k) *
                                 (pdf_params_stdev_chi_2(i, k) * pdf_params_stdev_chi_2(i, k))) -
                        (rcm(i, k) * rcm(i, k));
          rcp2(i, k)  = ((zero_threshold) > (rcp2(i, k)) ? (zero_threshold) : (rcp2(i, k)));
        });

#ifndef CLUBB_CAM
  }
#endif

  if ((iiPDF_type == iiPDF_ADG1 || iiPDF_type == iiPDF_ADG2 || iiPDF_type == iiPDF_new_hybrid) &&
      (var_on_stats_list(stats, "w_up_in_cloud") || var_on_stats_list(stats, "w_down_in_cloud"))) {
    calc_w_up_in_cloud(nz, ngrdcol,                                      // In
                       pdf_params.mixt_frac,                             // In
                       pdf_params.cloud_frac_1, pdf_params.cloud_frac_2, // In
                       pdf_params.w_1, pdf_params.w_2,                   // In
                       pdf_params.varnce_w_1, pdf_params.varnce_w_2,     // In
                       w_up_in_cloud, w_down_in_cloud,                   // Out
                       cloudy_updraft_frac, cloudy_downdraft_frac);      // Out
  } else {
    Kokkos::parallel_for(
        "pdf_closure::zero_cloudy_updraft", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i                 = cell % ngrdcol;
          const int k                 = (cell / ngrdcol);
          w_up_in_cloud(i, k)         = zero;
          w_down_in_cloud(i, k)       = zero;
          cloudy_updraft_frac(i, k)   = zero;
          cloudy_downdraft_frac(i, k) = zero;
        });
  }

#ifdef TUNER
  // Check the first ten levels of column1 for reasonable component
  // temperatures (190–1000K), as required by the native tuning build.
  {
    const auto previous_precision = std::cerr.precision();
    std::cerr.precision(std::numeric_limits<real_t>::max_digits10);
    const auto thl_1_host      = to_host_space(workspace, pdf_params.thl_1);
    const auto thl_2_host      = to_host_space(workspace, pdf_params.thl_2);
    const auto err_header_host = to_host_space(workspace, err_info.err_header);
    auto err_code_host         = to_host_space(workspace, err_info.err_code);
    for (int i = 0; i < 10 && i < thl_1_host.extent_int(1); ++i) {
      if (thl_1_host(0, i) < real_t(190)) {
        std::cerr.write(err_header_host.data(), err_header_host.extent(0)) << '\n';
        std::cerr << "Fatal error: pdf_params%thl_1 =" << thl_1_host(0, i)
                  << " < 190K at first grid column and grid level i = " << i + 1 << '\n';
        err_code_host(0) = clubb_fatal_error;
      }
      if (thl_2_host(0, i) < real_t(190)) {
        std::cerr.write(err_header_host.data(), err_header_host.extent(0)) << '\n';
        std::cerr << "Fatal error: pdf_params%thl_2 =" << thl_2_host(0, i)
                  << " < 190K at first grid column and grid level i = " << i + 1 << '\n';
        err_code_host(0) = clubb_fatal_error;
      }
      if (thl_1_host(0, i) > real_t(1000)) {
        std::cerr.write(err_header_host.data(), err_header_host.extent(0)) << '\n';
        std::cerr << "Fatal error: pdf_params%thl_1 =" << thl_1_host(0, i)
                  << " > 1000K at first grid column and grid level i = " << i + 1 << '\n';
        err_code_host(0) = clubb_fatal_error;
      }
      if (thl_2_host(0, i) > real_t(1000)) {
        std::cerr.write(err_header_host.data(), err_header_host.extent(0)) << '\n';
        std::cerr << "Fatal error: pdf_params%thl_2 =" << thl_2_host(0, i)
                  << " > 1000K at first grid column and grid level i = " << i + 1 << '\n';
        err_code_host(0) = clubb_fatal_error;
      }
    }
    Kokkos::deep_copy(err_info.err_code, err_code_host);
    std::cerr.precision(previous_precision);
  }
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error))
      return;
  }
#endif

  if (clubb_at_least_debug_level_api(2)) {
    for (int i = 0; i < ngrdcol; ++i) {
      pdf_closure_check(
          nz, sclr_dim, Kokkos::subview(wp4, i, Kokkos::ALL),
          Kokkos::subview(wprtp2, i, Kokkos::ALL), Kokkos::subview(wp2rtp, i, Kokkos::ALL),
          Kokkos::subview(wpthlp2, i, Kokkos::ALL), // intent(in)
          Kokkos::subview(wp2thlp, i, Kokkos::ALL), Kokkos::subview(cloud_frac, i, Kokkos::ALL),
          Kokkos::subview(rcm, i, Kokkos::ALL), Kokkos::subview(wpthvp, i, Kokkos::ALL),
          Kokkos::subview(wp2thvp, i, Kokkos::ALL), // intent(in)
          Kokkos::subview(wp2up, i, Kokkos::ALL),   // intent(in)
          Kokkos::subview(rtpthvp, i, Kokkos::ALL), Kokkos::subview(thlpthvp, i, Kokkos::ALL),
          Kokkos::subview(wprcp, i, Kokkos::ALL),
          Kokkos::subview(wp2rcp, i, Kokkos::ALL), // intent(in)
          Kokkos::subview(rtprcp, i, Kokkos::ALL), Kokkos::subview(thlprcp, i, Kokkos::ALL),
          Kokkos::subview(rcp2, i, Kokkos::ALL),
          Kokkos::subview(wprtpthlp, i, Kokkos::ALL), // intent(in)
          Kokkos::subview(pdf_params.crt_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.crt_2, i, Kokkos::ALL), // intent(in)
          Kokkos::subview(pdf_params.cthl_1, i, Kokkos::ALL),
          Kokkos::subview(pdf_params.cthl_2, i, Kokkos::ALL), // intent(in)
          pdf_params,                                         // intent(in)
          Kokkos::subview(sclrpthvp, i, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(sclrprcp, i, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(wpsclrp2, i, Kokkos::ALL, Kokkos::ALL), // intent(in)
          Kokkos::subview(wpsclrprtp, i, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(wpsclrpthlp, i, Kokkos::ALL, Kokkos::ALL),
          Kokkos::subview(wp2sclrp, i, Kokkos::ALL, Kokkos::ALL), // intent(in)
          stats,                                                  // intent(in)
          err_info);                                              // intent(inout)
      const auto err_code_host = to_host_space(workspace, err_info.err_code);
      if (err_code_host(i) == clubb_fatal_error) {
        const auto err_header_host = to_host_space(workspace, err_info.err_header);
        std::cerr.write(err_header_host.data() + i * err_header_host.extent(0),
                        err_header_host.extent(0))
            << '\n';
      }
    }
  }

  // Error reporting. Joshua Fasching, February2008.
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      const auto previous_precision = std::cerr.precision();
      std::cerr.precision(std::numeric_limits<real_t>::max_digits10);
      std::cerr << "Error in pdf_closure_new\n";
      std::cerr << "Intent(in)\n";
      {
        const auto p_in_Pa_host = to_host_space(workspace, p_in_Pa);
        std::cerr << "p_in_Pa =  ";
        for (std::size_t i1 = 0; i1 < p_in_Pa_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < p_in_Pa_host.extent(0); ++i0)
            std::cerr << p_in_Pa_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto exner_host = to_host_space(workspace, exner);
        std::cerr << "exner =  ";
        for (std::size_t i1 = 0; i1 < exner_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < exner_host.extent(0); ++i0)
            std::cerr << exner_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto thv_ds_host = to_host_space(workspace, thv_ds);
        std::cerr << "thv_ds =  ";
        for (std::size_t i1 = 0; i1 < thv_ds_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < thv_ds_host.extent(0); ++i0)
            std::cerr << thv_ds_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto wm_host = to_host_space(workspace, wm);
        std::cerr << "wm =  ";
        for (std::size_t i1 = 0; i1 < wm_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wm_host.extent(0); ++i0)
            std::cerr << wm_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto wp2_host = to_host_space(workspace, wp2);
        std::cerr << "wp2 =  ";
        for (std::size_t i1 = 0; i1 < wp2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wp2_host.extent(0); ++i0)
            std::cerr << wp2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto wp3_host = to_host_space(workspace, wp3);
        std::cerr << "wp3 =  ";
        for (std::size_t i1 = 0; i1 < wp3_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wp3_host.extent(0); ++i0)
            std::cerr << wp3_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto sigma_sqd_w_host = to_host_space(workspace, sigma_sqd_w);
        std::cerr << "sigma_sqd_w =  ";
        for (std::size_t i1 = 0; i1 < sigma_sqd_w_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < sigma_sqd_w_host.extent(0); ++i0)
            std::cerr << sigma_sqd_w_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto rtm_host = to_host_space(workspace, rtm);
        std::cerr << "rtm =  ";
        for (std::size_t i1 = 0; i1 < rtm_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < rtm_host.extent(0); ++i0)
            std::cerr << rtm_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto rtp2_host = to_host_space(workspace, rtp2);
        std::cerr << "rtp2 =  ";
        for (std::size_t i1 = 0; i1 < rtp2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < rtp2_host.extent(0); ++i0)
            std::cerr << rtp2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto wprtp_host = to_host_space(workspace, wprtp);
        std::cerr << "wprtp =  ";
        for (std::size_t i1 = 0; i1 < wprtp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wprtp_host.extent(0); ++i0)
            std::cerr << wprtp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto thlm_host = to_host_space(workspace, thlm);
        std::cerr << "thlm =  ";
        for (std::size_t i1 = 0; i1 < thlm_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < thlm_host.extent(0); ++i0)
            std::cerr << thlm_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto thlp2_host = to_host_space(workspace, thlp2);
        std::cerr << "thlp2 =  ";
        for (std::size_t i1 = 0; i1 < thlp2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < thlp2_host.extent(0); ++i0)
            std::cerr << thlp2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto wpthlp_host = to_host_space(workspace, wpthlp);
        std::cerr << "wpthlp =  ";
        for (std::size_t i1 = 0; i1 < wpthlp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wpthlp_host.extent(0); ++i0)
            std::cerr << wpthlp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto rtpthlp_host = to_host_space(workspace, rtpthlp);
        std::cerr << "rtpthlp =  ";
        for (std::size_t i1 = 0; i1 < rtpthlp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < rtpthlp_host.extent(0); ++i0)
            std::cerr << rtpthlp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      if (sclr_dim > 0) {
        {
          const auto sclrm_host = to_host_space(workspace, sclrm);
          std::cerr << "sclrm =  ";
          for (std::size_t i2 = 0; i2 < sclrm_host.extent(2); ++i2)
            for (std::size_t i1 = 0; i1 < sclrm_host.extent(1); ++i1)
              for (std::size_t i0 = 0; i0 < sclrm_host.extent(0); ++i0)
                std::cerr << sclrm_host(i0, i1, i2) << ' ';
          std::cerr << '\n';
        }
        {
          const auto wpsclrp_host = to_host_space(workspace, wpsclrp);
          std::cerr << "wpsclrp =  ";
          for (std::size_t i2 = 0; i2 < wpsclrp_host.extent(2); ++i2)
            for (std::size_t i1 = 0; i1 < wpsclrp_host.extent(1); ++i1)
              for (std::size_t i0 = 0; i0 < wpsclrp_host.extent(0); ++i0)
                std::cerr << wpsclrp_host(i0, i1, i2) << ' ';
          std::cerr << '\n';
        }
        {
          const auto sclrp2_host = to_host_space(workspace, sclrp2);
          std::cerr << "sclrp2 =  ";
          for (std::size_t i2 = 0; i2 < sclrp2_host.extent(2); ++i2)
            for (std::size_t i1 = 0; i1 < sclrp2_host.extent(1); ++i1)
              for (std::size_t i0 = 0; i0 < sclrp2_host.extent(0); ++i0)
                std::cerr << sclrp2_host(i0, i1, i2) << ' ';
          std::cerr << '\n';
        }
        {
          const auto sclrprtp_host = to_host_space(workspace, sclrprtp);
          std::cerr << "sclrprtp =  ";
          for (std::size_t i2 = 0; i2 < sclrprtp_host.extent(2); ++i2)
            for (std::size_t i1 = 0; i1 < sclrprtp_host.extent(1); ++i1)
              for (std::size_t i0 = 0; i0 < sclrprtp_host.extent(0); ++i0)
                std::cerr << sclrprtp_host(i0, i1, i2) << ' ';
          std::cerr << '\n';
        }
        {
          const auto sclrpthlp_host = to_host_space(workspace, sclrpthlp);
          std::cerr << "sclrpthlp =  ";
          for (std::size_t i2 = 0; i2 < sclrpthlp_host.extent(2); ++i2)
            for (std::size_t i1 = 0; i1 < sclrpthlp_host.extent(1); ++i1)
              for (std::size_t i0 = 0; i0 < sclrpthlp_host.extent(0); ++i0)
                std::cerr << sclrpthlp_host(i0, i1, i2) << ' ';
          std::cerr << '\n';
        }
      }
      std::cerr << "Intent(out)\n";
      {
        const auto wp4_host = to_host_space(workspace, wp4);
        std::cerr << "wp4 =  ";
        for (std::size_t i1 = 0; i1 < wp4_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wp4_host.extent(0); ++i0)
            std::cerr << wp4_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      if (l_explicit_turbulent_adv_xpyp || var_on_stats_list(stats, "wprtp2")) {
        {
          const auto wprtp2_host = to_host_space(workspace, wprtp2);
          std::cerr << "wprtp2 =  ";
          for (std::size_t i1 = 0; i1 < wprtp2_host.extent(1); ++i1)
            for (std::size_t i0 = 0; i0 < wprtp2_host.extent(0); ++i0)
              std::cerr << wprtp2_host(i0, i1) << ' ';
          std::cerr << '\n';
        }
      }
      {
        const auto wp2rtp_host = to_host_space(workspace, wp2rtp);
        std::cerr << "wp2rtp =  ";
        for (std::size_t i1 = 0; i1 < wp2rtp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wp2rtp_host.extent(0); ++i0)
            std::cerr << wp2rtp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      if (l_explicit_turbulent_adv_xpyp || var_on_stats_list(stats, "wpthlp2")) {
        {
          const auto wpthlp2_host = to_host_space(workspace, wpthlp2);
          std::cerr << "wpthlp2 =  ";
          for (std::size_t i1 = 0; i1 < wpthlp2_host.extent(1); ++i1)
            for (std::size_t i0 = 0; i0 < wpthlp2_host.extent(0); ++i0)
              std::cerr << wpthlp2_host(i0, i1) << ' ';
          std::cerr << '\n';
        }
      }
      {
        const auto cloud_frac_host = to_host_space(workspace, cloud_frac);
        std::cerr << "cloud_frac =  ";
        for (std::size_t i1 = 0; i1 < cloud_frac_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < cloud_frac_host.extent(0); ++i0)
            std::cerr << cloud_frac_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto ice_supersat_frac_host = to_host_space(workspace, ice_supersat_frac);
        std::cerr << "ice_supersat_frac =  ";
        for (std::size_t i1 = 0; i1 < ice_supersat_frac_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < ice_supersat_frac_host.extent(0); ++i0)
            std::cerr << ice_supersat_frac_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto rcm_host = to_host_space(workspace, rcm);
        std::cerr << "rcm =  ";
        for (std::size_t i1 = 0; i1 < rcm_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < rcm_host.extent(0); ++i0)
            std::cerr << rcm_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto wpthvp_host = to_host_space(workspace, wpthvp);
        std::cerr << "wpthvp =  ";
        for (std::size_t i1 = 0; i1 < wpthvp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wpthvp_host.extent(0); ++i0)
            std::cerr << wpthvp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto wp2thvp_host = to_host_space(workspace, wp2thvp);
        std::cerr << "wp2thvp =  ";
        for (std::size_t i1 = 0; i1 < wp2thvp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wp2thvp_host.extent(0); ++i0)
            std::cerr << wp2thvp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto wp2up_host = to_host_space(workspace, wp2up);
        std::cerr << "wp2up =  ";
        for (std::size_t i1 = 0; i1 < wp2up_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wp2up_host.extent(0); ++i0)
            std::cerr << wp2up_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto rtpthvp_host = to_host_space(workspace, rtpthvp);
        std::cerr << "rtpthvp =  ";
        for (std::size_t i1 = 0; i1 < rtpthvp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < rtpthvp_host.extent(0); ++i0)
            std::cerr << rtpthvp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto thlpthvp_host = to_host_space(workspace, thlpthvp);
        std::cerr << "thlpthvp =  ";
        for (std::size_t i1 = 0; i1 < thlpthvp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < thlpthvp_host.extent(0); ++i0)
            std::cerr << thlpthvp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto wprcp_host = to_host_space(workspace, wprcp);
        std::cerr << "wprcp =  ";
        for (std::size_t i1 = 0; i1 < wprcp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wprcp_host.extent(0); ++i0)
            std::cerr << wprcp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto wp2rcp_host = to_host_space(workspace, wp2rcp);
        std::cerr << "wp2rcp =  ";
        for (std::size_t i1 = 0; i1 < wp2rcp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wp2rcp_host.extent(0); ++i0)
            std::cerr << wp2rcp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto rtprcp_host = to_host_space(workspace, rtprcp);
        std::cerr << "rtprcp =  ";
        for (std::size_t i1 = 0; i1 < rtprcp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < rtprcp_host.extent(0); ++i0)
            std::cerr << rtprcp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto thlprcp_host = to_host_space(workspace, thlprcp);
        std::cerr << "thlprcp =  ";
        for (std::size_t i1 = 0; i1 < thlprcp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < thlprcp_host.extent(0); ++i0)
            std::cerr << thlprcp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
#ifndef CLUBB_CAM
      // if CLUBB is used in CAM we want this variable computed no matter what
      if (var_on_stats_list(stats, "rcp2")) {
#endif
        {
          const auto rcp2_host = to_host_space(workspace, rcp2);
          std::cerr << "rcp2 =  ";
          for (std::size_t i1 = 0; i1 < rcp2_host.extent(1); ++i1)
            for (std::size_t i0 = 0; i0 < rcp2_host.extent(0); ++i0)
              std::cerr << rcp2_host(i0, i1) << ' ';
          std::cerr << '\n';
        }
#ifndef CLUBB_CAM
      }
#endif
      if (l_explicit_turbulent_adv_xpyp || var_on_stats_list(stats, "wprtpthlp")) {
        {
          const auto wprtpthlp_host = to_host_space(workspace, wprtpthlp);
          std::cerr << "wprtpthlp =  ";
          for (std::size_t i1 = 0; i1 < wprtpthlp_host.extent(1); ++i1)
            for (std::size_t i0 = 0; i0 < wprtpthlp_host.extent(0); ++i0)
              std::cerr << wprtpthlp_host(i0, i1) << ' ';
          std::cerr << '\n';
        }
      }
      {
        const auto rcp2_host = to_host_space(workspace, rcp2);
        std::cerr << "rcp2 =  ";
        for (std::size_t i1 = 0; i1 < rcp2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < rcp2_host.extent(0); ++i0)
            std::cerr << rcp2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto wprtpthlp_host = to_host_space(workspace, wprtpthlp);
        std::cerr << "wprtpthlp =  ";
        for (std::size_t i1 = 0; i1 < wprtpthlp_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < wprtpthlp_host.extent(0); ++i0)
            std::cerr << wprtpthlp_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_w_1_host = to_host_space(workspace, pdf_params.w_1);
        std::cerr << "pdf_params%w_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_w_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_w_1_host.extent(0); ++i0)
            std::cerr << pdf_params_w_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_w_2_host = to_host_space(workspace, pdf_params.w_2);
        std::cerr << "pdf_params%w_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_w_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_w_2_host.extent(0); ++i0)
            std::cerr << pdf_params_w_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_varnce_w_1_host = to_host_space(workspace, pdf_params.varnce_w_1);
        std::cerr << "pdf_params%varnce_w_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_varnce_w_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_varnce_w_1_host.extent(0); ++i0)
            std::cerr << pdf_params_varnce_w_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_varnce_w_2_host = to_host_space(workspace, pdf_params.varnce_w_2);
        std::cerr << "pdf_params%varnce_w_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_varnce_w_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_varnce_w_2_host.extent(0); ++i0)
            std::cerr << pdf_params_varnce_w_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_rt_1_host = to_host_space(workspace, pdf_params.rt_1);
        std::cerr << "pdf_params%rt_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_rt_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_rt_1_host.extent(0); ++i0)
            std::cerr << pdf_params_rt_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_rt_2_host = to_host_space(workspace, pdf_params.rt_2);
        std::cerr << "pdf_params%rt_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_rt_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_rt_2_host.extent(0); ++i0)
            std::cerr << pdf_params_rt_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_varnce_rt_1_host = to_host_space(workspace, pdf_params.varnce_rt_1);
        std::cerr << "pdf_params%varnce_rt_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_varnce_rt_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_varnce_rt_1_host.extent(0); ++i0)
            std::cerr << pdf_params_varnce_rt_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_varnce_rt_2_host = to_host_space(workspace, pdf_params.varnce_rt_2);
        std::cerr << "pdf_params%varnce_rt_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_varnce_rt_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_varnce_rt_2_host.extent(0); ++i0)
            std::cerr << pdf_params_varnce_rt_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_thl_1_host = to_host_space(workspace, pdf_params.thl_1);
        std::cerr << "pdf_params%thl_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_thl_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_thl_1_host.extent(0); ++i0)
            std::cerr << pdf_params_thl_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_thl_2_host = to_host_space(workspace, pdf_params.thl_2);
        std::cerr << "pdf_params%thl_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_thl_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_thl_2_host.extent(0); ++i0)
            std::cerr << pdf_params_thl_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_varnce_thl_1_host = to_host_space(workspace, pdf_params.varnce_thl_1);
        std::cerr << "pdf_params%varnce_thl_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_varnce_thl_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_varnce_thl_1_host.extent(0); ++i0)
            std::cerr << pdf_params_varnce_thl_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_varnce_thl_2_host = to_host_space(workspace, pdf_params.varnce_thl_2);
        std::cerr << "pdf_params%varnce_thl_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_varnce_thl_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_varnce_thl_2_host.extent(0); ++i0)
            std::cerr << pdf_params_varnce_thl_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_w_rt_1_host = to_host_space(workspace, pdf_params.corr_w_rt_1);
        std::cerr << "pdf_params%corr_w_rt_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_w_rt_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_w_rt_1_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_w_rt_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_w_rt_2_host = to_host_space(workspace, pdf_params.corr_w_rt_2);
        std::cerr << "pdf_params%corr_w_rt_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_w_rt_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_w_rt_2_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_w_rt_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_w_thl_1_host = to_host_space(workspace, pdf_params.corr_w_thl_1);
        std::cerr << "pdf_params%corr_w_thl_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_w_thl_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_w_thl_1_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_w_thl_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_w_thl_2_host = to_host_space(workspace, pdf_params.corr_w_thl_2);
        std::cerr << "pdf_params%corr_w_thl_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_w_thl_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_w_thl_2_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_w_thl_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_rt_thl_1_host =
            to_host_space(workspace, pdf_params.corr_rt_thl_1);
        std::cerr << "pdf_params%corr_rt_thl_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_rt_thl_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_rt_thl_1_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_rt_thl_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_rt_thl_2_host =
            to_host_space(workspace, pdf_params.corr_rt_thl_2);
        std::cerr << "pdf_params%corr_rt_thl_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_rt_thl_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_rt_thl_2_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_rt_thl_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_alpha_thl_host = to_host_space(workspace, pdf_params.alpha_thl);
        std::cerr << "pdf_params%alpha_thl =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_alpha_thl_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_alpha_thl_host.extent(0); ++i0)
            std::cerr << pdf_params_alpha_thl_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_alpha_rt_host = to_host_space(workspace, pdf_params.alpha_rt);
        std::cerr << "pdf_params%alpha_rt =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_alpha_rt_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_alpha_rt_host.extent(0); ++i0)
            std::cerr << pdf_params_alpha_rt_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_crt_1_host = to_host_space(workspace, pdf_params.crt_1);
        std::cerr << "pdf_params%crt_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_crt_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_crt_1_host.extent(0); ++i0)
            std::cerr << pdf_params_crt_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_crt_2_host = to_host_space(workspace, pdf_params.crt_2);
        std::cerr << "pdf_params%crt_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_crt_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_crt_2_host.extent(0); ++i0)
            std::cerr << pdf_params_crt_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_cthl_1_host = to_host_space(workspace, pdf_params.cthl_1);
        std::cerr << "pdf_params%cthl_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_cthl_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_cthl_1_host.extent(0); ++i0)
            std::cerr << pdf_params_cthl_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_cthl_2_host = to_host_space(workspace, pdf_params.cthl_2);
        std::cerr << "pdf_params%cthl_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_cthl_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_cthl_2_host.extent(0); ++i0)
            std::cerr << pdf_params_cthl_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_chi_1_host = to_host_space(workspace, pdf_params.chi_1);
        std::cerr << "pdf_params%chi_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_chi_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_chi_1_host.extent(0); ++i0)
            std::cerr << pdf_params_chi_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_chi_2_host = to_host_space(workspace, pdf_params.chi_2);
        std::cerr << "pdf_params%chi_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_chi_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_chi_2_host.extent(0); ++i0)
            std::cerr << pdf_params_chi_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_stdev_chi_1_host = to_host_space(workspace, pdf_params.stdev_chi_1);
        std::cerr << "pdf_params%stdev_chi_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_stdev_chi_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_stdev_chi_1_host.extent(0); ++i0)
            std::cerr << pdf_params_stdev_chi_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_stdev_chi_2_host = to_host_space(workspace, pdf_params.stdev_chi_2);
        std::cerr << "pdf_params%stdev_chi_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_stdev_chi_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_stdev_chi_2_host.extent(0); ++i0)
            std::cerr << pdf_params_stdev_chi_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_stdev_eta_1_host = to_host_space(workspace, pdf_params.stdev_eta_1);
        std::cerr << "pdf_params%stdev_eta_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_stdev_eta_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_stdev_eta_1_host.extent(0); ++i0)
            std::cerr << pdf_params_stdev_eta_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_stdev_eta_2_host = to_host_space(workspace, pdf_params.stdev_eta_2);
        std::cerr << "pdf_params%stdev_eta_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_stdev_eta_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_stdev_eta_2_host.extent(0); ++i0)
            std::cerr << pdf_params_stdev_eta_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_covar_chi_eta_1_host =
            to_host_space(workspace, pdf_params.covar_chi_eta_1);
        std::cerr << "pdf_params%covar_chi_eta_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_covar_chi_eta_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_covar_chi_eta_1_host.extent(0); ++i0)
            std::cerr << pdf_params_covar_chi_eta_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_covar_chi_eta_2_host =
            to_host_space(workspace, pdf_params.covar_chi_eta_2);
        std::cerr << "pdf_params%covar_chi_eta_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_covar_chi_eta_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_covar_chi_eta_2_host.extent(0); ++i0)
            std::cerr << pdf_params_covar_chi_eta_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_w_chi_1_host = to_host_space(workspace, pdf_params.corr_w_chi_1);
        std::cerr << "pdf_params%corr_w_chi_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_w_chi_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_w_chi_1_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_w_chi_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_w_chi_2_host = to_host_space(workspace, pdf_params.corr_w_chi_2);
        std::cerr << "pdf_params%corr_w_chi_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_w_chi_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_w_chi_2_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_w_chi_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_w_eta_1_host = to_host_space(workspace, pdf_params.corr_w_eta_1);
        std::cerr << "pdf_params%corr_w_eta_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_w_eta_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_w_eta_1_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_w_eta_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_w_eta_2_host = to_host_space(workspace, pdf_params.corr_w_eta_2);
        std::cerr << "pdf_params%corr_w_eta_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_w_eta_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_w_eta_2_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_w_eta_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_chi_eta_1_host =
            to_host_space(workspace, pdf_params.corr_chi_eta_1);
        std::cerr << "pdf_params%corr_chi_eta_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_chi_eta_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_chi_eta_1_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_chi_eta_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_corr_chi_eta_2_host =
            to_host_space(workspace, pdf_params.corr_chi_eta_2);
        std::cerr << "pdf_params%corr_chi_eta_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_corr_chi_eta_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_corr_chi_eta_2_host.extent(0); ++i0)
            std::cerr << pdf_params_corr_chi_eta_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_rsatl_1_host = to_host_space(workspace, pdf_params.rsatl_1);
        std::cerr << "pdf_params%rsatl_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_rsatl_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_rsatl_1_host.extent(0); ++i0)
            std::cerr << pdf_params_rsatl_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_rsatl_2_host = to_host_space(workspace, pdf_params.rsatl_2);
        std::cerr << "pdf_params%rsatl_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_rsatl_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_rsatl_2_host.extent(0); ++i0)
            std::cerr << pdf_params_rsatl_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_rc_1_host = to_host_space(workspace, pdf_params.rc_1);
        std::cerr << "pdf_params%rc_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_rc_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_rc_1_host.extent(0); ++i0)
            std::cerr << pdf_params_rc_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_rc_2_host = to_host_space(workspace, pdf_params.rc_2);
        std::cerr << "pdf_params%rc_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_rc_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_rc_2_host.extent(0); ++i0)
            std::cerr << pdf_params_rc_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_cloud_frac_1_host = to_host_space(workspace, pdf_params.cloud_frac_1);
        std::cerr << "pdf_params%cloud_frac_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_cloud_frac_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_cloud_frac_1_host.extent(0); ++i0)
            std::cerr << pdf_params_cloud_frac_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_cloud_frac_2_host = to_host_space(workspace, pdf_params.cloud_frac_2);
        std::cerr << "pdf_params%cloud_frac_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_cloud_frac_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_cloud_frac_2_host.extent(0); ++i0)
            std::cerr << pdf_params_cloud_frac_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_mixt_frac_host = to_host_space(workspace, pdf_params.mixt_frac);
        std::cerr << "pdf_params%mixt_frac =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_mixt_frac_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_mixt_frac_host.extent(0); ++i0)
            std::cerr << pdf_params_mixt_frac_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_ice_supersat_frac_1_host =
            to_host_space(workspace, pdf_params.ice_supersat_frac_1);
        std::cerr << "pdf_params%ice_supersat_frac_1 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_ice_supersat_frac_1_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_ice_supersat_frac_1_host.extent(0); ++i0)
            std::cerr << pdf_params_ice_supersat_frac_1_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      {
        const auto pdf_params_ice_supersat_frac_2_host =
            to_host_space(workspace, pdf_params.ice_supersat_frac_2);
        std::cerr << "pdf_params%ice_supersat_frac_2 =  ";
        for (std::size_t i1 = 0; i1 < pdf_params_ice_supersat_frac_2_host.extent(1); ++i1)
          for (std::size_t i0 = 0; i0 < pdf_params_ice_supersat_frac_2_host.extent(0); ++i0)
            std::cerr << pdf_params_ice_supersat_frac_2_host(i0, i1) << ' ';
        std::cerr << '\n';
      }
      if (sclr_dim > 0) {
        {
          const auto sclrpthvp_host = to_host_space(workspace, sclrpthvp);
          std::cerr << "sclrpthvp =  ";
          for (std::size_t i2 = 0; i2 < sclrpthvp_host.extent(2); ++i2)
            for (std::size_t i1 = 0; i1 < sclrpthvp_host.extent(1); ++i1)
              for (std::size_t i0 = 0; i0 < sclrpthvp_host.extent(0); ++i0)
                std::cerr << sclrpthvp_host(i0, i1, i2) << ' ';
          std::cerr << '\n';
        }
        {
          const auto sclrprcp_host = to_host_space(workspace, sclrprcp);
          std::cerr << "sclrprcp =  ";
          for (std::size_t i2 = 0; i2 < sclrprcp_host.extent(2); ++i2)
            for (std::size_t i1 = 0; i1 < sclrprcp_host.extent(1); ++i1)
              for (std::size_t i0 = 0; i0 < sclrprcp_host.extent(0); ++i0)
                std::cerr << sclrprcp_host(i0, i1, i2) << ' ';
          std::cerr << '\n';
        }
        {
          const auto wpsclrp2_host = to_host_space(workspace, wpsclrp2);
          std::cerr << "wpsclrp2 =  ";
          for (std::size_t i2 = 0; i2 < wpsclrp2_host.extent(2); ++i2)
            for (std::size_t i1 = 0; i1 < wpsclrp2_host.extent(1); ++i1)
              for (std::size_t i0 = 0; i0 < wpsclrp2_host.extent(0); ++i0)
                std::cerr << wpsclrp2_host(i0, i1, i2) << ' ';
          std::cerr << '\n';
        }
        {
          const auto wpsclrprtp_host = to_host_space(workspace, wpsclrprtp);
          std::cerr << "wpsclrprtp =  ";
          for (std::size_t i2 = 0; i2 < wpsclrprtp_host.extent(2); ++i2)
            for (std::size_t i1 = 0; i1 < wpsclrprtp_host.extent(1); ++i1)
              for (std::size_t i0 = 0; i0 < wpsclrprtp_host.extent(0); ++i0)
                std::cerr << wpsclrprtp_host(i0, i1, i2) << ' ';
          std::cerr << '\n';
        }
        {
          const auto wpsclrpthlp_host = to_host_space(workspace, wpsclrpthlp);
          std::cerr << "wpsclrpthlp =  ";
          for (std::size_t i2 = 0; i2 < wpsclrpthlp_host.extent(2); ++i2)
            for (std::size_t i1 = 0; i1 < wpsclrpthlp_host.extent(1); ++i1)
              for (std::size_t i0 = 0; i0 < wpsclrpthlp_host.extent(0); ++i0)
                std::cerr << wpsclrpthlp_host(i0, i1, i2) << ' ';
          std::cerr << '\n';
        }
        {
          const auto wp2sclrp_host = to_host_space(workspace, wp2sclrp);
          std::cerr << "wp2sclrp =  ";
          for (std::size_t i2 = 0; i2 < wp2sclrp_host.extent(2); ++i2)
            for (std::size_t i1 = 0; i1 < wp2sclrp_host.extent(1); ++i1)
              for (std::size_t i0 = 0; i0 < wp2sclrp_host.extent(0); ++i0)
                std::cerr << wp2sclrp_host(i0, i1, i2) << ' ';
          std::cerr << '\n';
        }
      }

      std::cerr.precision(previous_precision);
      return;
    }

    // Error-check PDF parameters and moments for internal consistency.
    // Native Luhar diagnostics run serially on host mirrors to retain print
    // order. They do not change model state or close any prognostic equation.
    if (iiPDF_type == iiPDF_3D_Luhar) {
      const auto previous_precision = std::cerr.precision();
      std::cerr.precision(std::numeric_limits<real_t>::max_digits10);
      const auto Skrt_host                    = to_host_space(workspace, Skrt);
      const auto Skrt_clubb_pdf_host          = create_host_mirror(workspace, Skrt_clubb_pdf);
      const auto Skthl_host                   = to_host_space(workspace, Skthl);
      const auto Skthl_clubb_pdf_host         = create_host_mirror(workspace, Skthl_clubb_pdf);
      const auto Skw_host                     = to_host_space(workspace, Skw);
      const auto Skw_clubb_pdf_host           = create_host_mirror(workspace, Skw_clubb_pdf);
      const auto clubb_params_host            = to_host_space(workspace, clubb_params);
      const auto rtm_host                     = to_host_space(workspace, rtm);
      const auto rtm_clubb_pdf_host           = create_host_mirror(workspace, rtm_clubb_pdf);
      const auto rtp2_host                    = to_host_space(workspace, rtp2);
      const auto rtp2_clubb_pdf_host          = create_host_mirror(workspace, rtp2_clubb_pdf);
      const auto rtp3_clubb_pdf_host          = create_host_mirror(workspace, rtp3_clubb_pdf);
      const auto thlm_host                    = to_host_space(workspace, thlm);
      const auto thlm_clubb_pdf_host          = create_host_mirror(workspace, thlm_clubb_pdf);
      const auto thlp2_host                   = to_host_space(workspace, thlp2);
      const auto thlp2_clubb_pdf_host         = create_host_mirror(workspace, thlp2_clubb_pdf);
      const auto thlp3_clubb_pdf_host         = create_host_mirror(workspace, thlp3_clubb_pdf);
      const auto wm_host                      = to_host_space(workspace, wm);
      const auto wm_clubb_pdf_host            = create_host_mirror(workspace, wm_clubb_pdf);
      const auto wp2_host                     = to_host_space(workspace, wp2);
      const auto wp2_clubb_pdf_host           = create_host_mirror(workspace, wp2_clubb_pdf);
      const auto wp3_clubb_pdf_host           = create_host_mirror(workspace, wp3_clubb_pdf);
      const auto pdf_params_mixt_frac_host    = to_host_space(workspace, pdf_params.mixt_frac);
      const auto pdf_params_rt_1_host         = to_host_space(workspace, pdf_params.rt_1);
      const auto pdf_params_rt_2_host         = to_host_space(workspace, pdf_params.rt_2);
      const auto pdf_params_thl_1_host        = to_host_space(workspace, pdf_params.thl_1);
      const auto pdf_params_thl_2_host        = to_host_space(workspace, pdf_params.thl_2);
      const auto pdf_params_varnce_rt_1_host  = to_host_space(workspace, pdf_params.varnce_rt_1);
      const auto pdf_params_varnce_rt_2_host  = to_host_space(workspace, pdf_params.varnce_rt_2);
      const auto pdf_params_varnce_thl_1_host = to_host_space(workspace, pdf_params.varnce_thl_1);
      const auto pdf_params_varnce_thl_2_host = to_host_space(workspace, pdf_params.varnce_thl_2);
      const auto pdf_params_varnce_w_1_host   = to_host_space(workspace, pdf_params.varnce_w_1);
      const auto pdf_params_varnce_w_2_host   = to_host_space(workspace, pdf_params.varnce_w_2);
      const auto pdf_params_w_1_host          = to_host_space(workspace, pdf_params.w_1);
      const auto pdf_params_w_2_host          = to_host_space(workspace, pdf_params.w_2);
      real_t Skw_denom_coef;
      for (int i = 0; i < ngrdcol; ++i) {
        // Means
        for (int k = 0; k < nz; ++k)
          wm_clubb_pdf_host(i, k) =
              pdf_params_mixt_frac_host(i, k) * pdf_params_w_1_host(i, k) +
              (one - pdf_params_mixt_frac_host(i, k)) * pdf_params_w_2_host(i, k);
        for (int k = 0; k < nz; ++k) {
          if (Kokkos::abs((wm_clubb_pdf_host(i, k) - wm_host(i, k)) /
                          ((wm_host(i, k)) > (eps) ? (wm_host(i, k)) : (eps))) > real_t(.05)) {
            std::cerr << "wm error at thlm = " << ' ' << thlm_host(i, k) << ' '
                      << ((wm_clubb_pdf_host(i, k) - wm_host(i, k)) /
                          ((wm_host(i, k)) > (eps) ? (wm_host(i, k)) : (eps)))
                      << '\n';
          }
        }
        for (int k = 0; k < nz; ++k)
          rtm_clubb_pdf_host(i, k) =
              pdf_params_mixt_frac_host(i, k) * pdf_params_rt_1_host(i, k) +
              (one - pdf_params_mixt_frac_host(i, k)) * pdf_params_rt_2_host(i, k);
        for (int k = 0; k < nz; ++k) {
          if (Kokkos::abs((rtm_clubb_pdf_host(i, k) - rtm_host(i, k)) /
                          ((rtm_host(i, k)) > (eps) ? (rtm_host(i, k)) : (eps))) > real_t(.05)) {
            std::cerr << "rtm error at thlm = " << ' ' << thlm_host(i, k) << ' '
                      << ((rtm_clubb_pdf_host(i, k) - rtm_host(i, k)) /
                          ((rtm_host(i, k)) > (eps) ? (rtm_host(i, k)) : (eps)))
                      << '\n';
          }
        }
        for (int k = 0; k < nz; ++k)
          thlm_clubb_pdf_host(i, k) =
              pdf_params_mixt_frac_host(i, k) * pdf_params_thl_1_host(i, k) +
              (one - pdf_params_mixt_frac_host(i, k)) * pdf_params_thl_2_host(i, k);
        for (int k = 0; k < nz; ++k) {
          if (Kokkos::abs((thlm_clubb_pdf_host(i, k) - thlm_host(i, k)) / thlm_host(i, k)) >
              real_t(.05)) {
            std::cerr << "thlm error at thlm = " << ' ' << thlm_host(i, k) << ' '
                      << ((thlm_clubb_pdf_host(i, k) - thlm_host(i, k)) / thlm_host(i, k)) << '\n';
          }
        }
        // Variances
        for (int k = 0; k < nz; ++k)
          wp2_clubb_pdf_host(i, k) =
              pdf_params_mixt_frac_host(i, k) * (((pdf_params_w_1_host(i, k) - wm_host(i, k)) *
                                                  (pdf_params_w_1_host(i, k) - wm_host(i, k))) +
                                                 pdf_params_varnce_w_1_host(i, k)) +
              (one - pdf_params_mixt_frac_host(i, k)) *
                  (((pdf_params_w_2_host(i, k) - wm_host(i, k)) *
                    (pdf_params_w_2_host(i, k) - wm_host(i, k))) +
                   pdf_params_varnce_w_2_host(i, k));
        for (int k = 0; k < nz; ++k) {
          if (wp2_host(i, k) > (w_tol * w_tol)) {
            if (Kokkos::abs((wp2_clubb_pdf_host(i, k) - wp2_host(i, k)) / wp2_host(i, k)) >
                real_t(.05)) {
              std::cerr << "wp2 error at thlm = " << ' ' << thlm_host(i, k) << ' '
                        << ((wp2_clubb_pdf_host(i, k) - wp2_host(i, k)) / wp2_host(i, k)) << '\n';
            }
          }
        }
        for (int k = 0; k < nz; ++k)
          rtp2_clubb_pdf_host(i, k) =
              pdf_params_mixt_frac_host(i, k) * (((pdf_params_rt_1_host(i, k) - rtm_host(i, k)) *
                                                  (pdf_params_rt_1_host(i, k) - rtm_host(i, k))) +
                                                 pdf_params_varnce_rt_1_host(i, k)) +
              (one - pdf_params_mixt_frac_host(i, k)) *
                  (((pdf_params_rt_2_host(i, k) - rtm_host(i, k)) *
                    (pdf_params_rt_2_host(i, k) - rtm_host(i, k))) +
                   pdf_params_varnce_rt_2_host(i, k));
        for (int k = 0; k < nz; ++k) {
          if (rtp2_host(i, k) > (rt_tol * rt_tol)) {
            if (Kokkos::abs((rtp2_clubb_pdf_host(i, k) - rtp2_host(i, k)) / rtp2_host(i, k)) >
                real_t(.05)) {
              std::cerr << "rtp2 error at thlm = " << ' ' << thlm_host(i, k) << ' '
                        << "Error = " << ' '
                        << ((rtp2_clubb_pdf_host(i, k) - rtp2_host(i, k)) / rtp2_host(i, k))
                        << '\n';
            }
          }
        }
        for (int k = 0; k < nz; ++k)
          thlp2_clubb_pdf_host(i, k) =
              pdf_params_mixt_frac_host(i, k) * (((pdf_params_thl_1_host(i, k) - thlm_host(i, k)) *
                                                  (pdf_params_thl_1_host(i, k) - thlm_host(i, k))) +
                                                 pdf_params_varnce_thl_1_host(i, k)) +
              (one - pdf_params_mixt_frac_host(i, k)) *
                  (((pdf_params_thl_2_host(i, k) - thlm_host(i, k)) *
                    (pdf_params_thl_2_host(i, k) - thlm_host(i, k))) +
                   pdf_params_varnce_thl_2_host(i, k));
        for (int k = 0; k < nz; ++k) {
          if (thlp2_host(i, k) > (thl_tol * thl_tol)) {
            if (Kokkos::abs((thlp2_clubb_pdf_host(i, k) - thlp2_host(i, k)) / thlp2_host(i, k)) >
                real_t(.05)) {
              std::cerr << "thlp2 error at thlm = " << ' ' << thlm_host(i, k) << ' '
                        << "Error = " << ' '
                        << ((thlp2_clubb_pdf_host(i, k) - thlp2_host(i, k)) / thlp2_host(i, k))
                        << '\n';
            }
          }
        }
        // Third order moments
        for (int k = 0; k < nz; ++k)
          wp3_clubb_pdf_host(i, k) = pdf_params_mixt_frac_host(i, k) *
                                         (pdf_params_w_1_host(i, k) - wm_host(i, k)) *
                                         (((pdf_params_w_1_host(i, k) - wm_host(i, k)) *
                                           (pdf_params_w_1_host(i, k) - wm_host(i, k))) +
                                          three * pdf_params_varnce_w_1_host(i, k)) +
                                     (one - pdf_params_mixt_frac_host(i, k)) *
                                         (pdf_params_w_2_host(i, k) - wm_host(i, k)) *
                                         (((pdf_params_w_2_host(i, k) - wm_host(i, k)) *
                                           (pdf_params_w_2_host(i, k) - wm_host(i, k))) +
                                          three * pdf_params_varnce_w_2_host(i, k));
        for (int k = 0; k < nz; ++k)
          rtp3_clubb_pdf_host(i, k) = pdf_params_mixt_frac_host(i, k) *
                                          (pdf_params_rt_1_host(i, k) - rtm_host(i, k)) *
                                          (((pdf_params_rt_1_host(i, k) - rtm_host(i, k)) *
                                            (pdf_params_rt_1_host(i, k) - rtm_host(i, k))) +
                                           three * pdf_params_varnce_rt_1_host(i, k)) +
                                      (one - pdf_params_mixt_frac_host(i, k)) *
                                          (pdf_params_rt_2_host(i, k) - rtm_host(i, k)) *
                                          (((pdf_params_rt_2_host(i, k) - rtm_host(i, k)) *
                                            (pdf_params_rt_2_host(i, k) - rtm_host(i, k))) +
                                           three * pdf_params_varnce_rt_2_host(i, k));
        for (int k = 0; k < nz; ++k)
          thlp3_clubb_pdf_host(i, k) = pdf_params_mixt_frac_host(i, k) *
                                           (pdf_params_thl_1_host(i, k) - thlm_host(i, k)) *
                                           (((pdf_params_thl_1_host(i, k) - thlm_host(i, k)) *
                                             (pdf_params_thl_1_host(i, k) - thlm_host(i, k))) +
                                            three * pdf_params_varnce_thl_1_host(i, k)) +
                                       (one - pdf_params_mixt_frac_host(i, k)) *
                                           (pdf_params_thl_2_host(i, k) - thlm_host(i, k)) *
                                           (((pdf_params_thl_2_host(i, k) - thlm_host(i, k)) *
                                             (pdf_params_thl_2_host(i, k) - thlm_host(i, k))) +
                                            three * pdf_params_varnce_thl_2_host(i, k));
        // Skewness
        Skw_denom_coef = clubb_params_host(i, iSkw_denom_coef);
        for (int k = 0; k < nz; ++k)
          Skw_clubb_pdf_host(i, k) =
              wp3_clubb_pdf_host(i, k) /
              Kokkos::pow((wp2_clubb_pdf_host(i, k) + Skw_denom_coef * (w_tol * w_tol)),
                          real_t(1.5));
        for (int k = 0; k < nz; ++k) {
          if (Skw_host(i, k) > real_t(.05)) {
            if (Kokkos::abs((Skw_clubb_pdf_host(i, k) - Skw_host(i, k)) / Skw_host(i, k)) >
                real_t(.25)) {
              std::cerr << "Skw error at thlm = " << ' ' << thlm_host(i, k) << ' '
                        << "Error = " << ' '
                        << ((Skw_clubb_pdf_host(i, k) - Skw_host(i, k)) / Skw_host(i, k)) << ' '
                        << Skw_clubb_pdf_host(i, k) << ' ' << Skw_host(i, k) << '\n';
            }
          }
        }
        for (int k = 0; k < nz; ++k)
          Skrt_clubb_pdf_host(i, k) =
              rtp3_clubb_pdf_host(i, k) /
              Kokkos::pow((rtp2_clubb_pdf_host(i, k) + Skw_denom_coef * (rt_tol * rt_tol)),
                          real_t(1.5));
        for (int k = 0; k < nz; ++k) {
          if (Skrt_host(i, k) > real_t(.05)) {
            if (Kokkos::abs((Skrt_clubb_pdf_host(i, k) - Skrt_host(i, k)) / Skrt_host(i, k)) >
                real_t(.25)) {
              std::cerr << "Skrt error at thlm = " << ' ' << thlm_host(i, k) << ' '
                        << "Error = " << ' '
                        << ((Skrt_clubb_pdf_host(i, k) - Skrt_host(i, k)) / Skrt_host(i, k)) << ' '
                        << Skrt_clubb_pdf_host(i, k) << ' ' << Skrt_host(i, k) << '\n';
            }
          }
        }
        for (int k = 0; k < nz; ++k)
          Skthl_clubb_pdf_host(i, k) =
              thlp3_clubb_pdf_host(i, k) /
              Kokkos::pow((thlp2_clubb_pdf_host(i, k) + Skw_denom_coef * (thl_tol * thl_tol)),
                          real_t(1.5));
        for (int k = 0; k < nz; ++k) {
          if (Skthl_host(i, k) > real_t(.05)) {
            if (Kokkos::abs((Skthl_clubb_pdf_host(i, k) - Skthl_host(i, k)) / Skthl_host(i, k)) >
                real_t(.25)) {
              std::cerr << "Skthl error at thlm = " << ' ' << thlm_host(i, k) << ' '
                        << "Error = " << ' '
                        << ((Skthl_clubb_pdf_host(i, k) - Skthl_host(i, k)) / Skthl_host(i, k))
                        << ' ' << Skthl_clubb_pdf_host(i, k) << ' ' << Skthl_host(i, k) << '\n';
            }
          }
        }
      }
      std::cerr.precision(previous_precision);
    }
  }
}

namespace {

void transform_pdf_chi_eta_component(ScratchWorkspace& workspace,                      // InOut
                                     int nz, int ngrdcol,                              // In
                                     Real2D tl, Real2D rsatl, Real2D rt, Real2D exner, // In
                                     Real2D varnce_thl, Real2D varnce_rt,              // In
                                     Real2D corr_rt_thl, Real2D chi,                   // In/Out
                                     Real2D crt, Real2D cthl,                          // Out
                                     Real2D stdev_chi, Real2D stdev_eta,               // Out
                                     Real2D covar_chi_eta,                             // Out
                                     Real2D corr_chi_eta) {
  // Transform component rt/theta-l means and variances to chi/eta,
  // returning thermodynamic coefficients, spreads and their correlation.
  // Chi and eta have water-mixing-ratio units [kg/kg].
  auto workspace_scope = workspace.scope();
  Real2D denominator =
      workspace_scope.view<Real2D>("transform_pdf_chi_eta_component::denominator", ngrdcol, nz);
  constexpr real_t denom_thresh = chi_tol * eta_tol;
  constexpr real_t Cp_on_Lv     = Cp / Lv;
  Kokkos::parallel_for(
      "transform_pdf_chi_eta_component::transform", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = (cell / ngrdcol);
        real_t beta, invrs_beta_rsatl_p1;
        // SD's beta (eqn. 8)
        beta                = ep * (Lv * Lv) / (Rd * Cp * (tl(i, k) * tl(i, k)));
        invrs_beta_rsatl_p1 = one / (one + beta * rsatl(i, k));
        // s from Lewellen and Yoh 1993 (LY) eqn. 1
        chi(i, k) = (rt(i, k) - rsatl(i, k)) * invrs_beta_rsatl_p1;
        // For each normal distribution in the sum of two normal distributions,
        // s' = crt * rt'  +  cthl * thl';
        // therefore, x's' = crt * x'rt'  +  cthl * x'thl'.
        // Larson et al. May, 2001.
        crt(i, k)  = invrs_beta_rsatl_p1;
        cthl(i, k) = (one + beta * rt(i, k)) * (invrs_beta_rsatl_p1 * invrs_beta_rsatl_p1) *
                     Cp_on_Lv * beta * rsatl(i, k) * exner(i, k);
      });
  Kokkos::parallel_for(
      "transform_pdf_chi_eta_component::component_spreads", ngrdcol * nz,
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = (cell / ngrdcol);
        real_t varnce_rt_term, varnce_thl_term, corr_rt_thl_term, varnce_chi, varnce_eta;
        varnce_rt_term      = (crt(i, k) * crt(i, k)) * varnce_rt(i, k);
        varnce_thl_term     = (cthl(i, k) * cthl(i, k)) * varnce_thl(i, k);
        covar_chi_eta(i, k) = varnce_rt_term - varnce_thl_term;
    // Native rounds this product before the chi/eta sums. CUDA otherwise
    // contracts it across this assignment; __dmul_rn preserves that boundary
    // without a volatile memory spill (BFB-010 in BFB_ADAPTATIONS.md).
#if defined(__CUDA_ARCH__)
        corr_rt_thl_term = __dmul_rn(two * corr_rt_thl(i, k) * crt(i, k) * cthl(i, k),
                                     Kokkos::sqrt(varnce_rt(i, k) * varnce_thl(i, k)));
#else
        corr_rt_thl_term = two * corr_rt_thl(i, k) * crt(i, k) * cthl(i, k) *
                           Kokkos::sqrt(varnce_rt(i, k) * varnce_thl(i, k));
#endif
        varnce_chi        = varnce_rt_term - corr_rt_thl_term + varnce_thl_term;
        varnce_eta        = varnce_rt_term + corr_rt_thl_term + varnce_thl_term;
        stdev_chi(i, k)   = Kokkos::sqrt(varnce_chi);
        stdev_eta(i, k)   = Kokkos::sqrt(varnce_eta);
        denominator(i, k) = stdev_chi(i, k) * stdev_eta(i, k);
      });
  smooth_corr_quotient(ngrdcol, nz, covar_chi_eta, denominator, denom_thresh, corr_chi_eta);
}

} // namespace

void calc_wp4_pdf(int nz, int ngrdcol,                  // In
                  Real2D wm, Real2D w_1, Real2D w_2,    // In
                  Real2D varnce_w_1, Real2D varnce_w_2, // In
                  Real2D mixt_frac,                     // In
                  Real2D wp4) {
  // Description:
  // Calculates <w'^4> by integrating over the PDF of w.  The integral is:
  //
  // <w'^4> = INT(-inf:inf) ( w - <w> )^4 P(w) dw;
  //
  // where <w> is the overall mean of w and P(w) is a two-component normal
  // distribution of w.  The integrated equation is:
  //
  // <w'^4> = mixt_frac * ( 3 * sigma_w_1^4
  //                        + 6 * ( mu_w_1 - <w> )^2 * sigma_w_1^2
  //                        + ( mu_w_1 - <w> )^4 )
  //          + ( 1 - mixt_frac ) * ( 3 * sigma_w_2^4
  //                                  + 6 * ( mu_w_2 - <w> )^2 * sigma_w_2^2
  //                                  + ( mu_w_2 - <w> )^4 );
  //
  // where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
  // of w in the 2nd PDF component, sigma_w_1 is the standard deviation of w in
  // the 1st PDF component, sigma_w_2 is the standard deviation of w in the 2nd
  // PDF component, and mixt_frac is the mixture fraction, which is the weight
  // of the 1st PDF component.

  // References:
  //-----------------------------------------------------------------------

  // Arrays are indexed (column, level).
  // wm: Mean of w (overall) [m/s]
  // w_1: Mean of w (1st PDF component) [m/s]
  // w_2: Mean of w (2nd PDF component) [m/s]
  // varnce_w_1: Variance of w (1st PDF component) [m^2/s^2]
  // varnce_w_2: Variance of w (2nd PDF component) [m^2/s^2]
  // mixt_frac: Mixture fraction [-]
  // wp4: <w'^4> [m^4/s^4]

  // Calculate <w'^4> by integrating over the PDF.
  Kokkos::parallel_for(
      "calc_wp4_pdf", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        wp4(i, k) =
            mixt_frac(i, k) *
                (three * (varnce_w_1(i, k) * varnce_w_1(i, k)) +
                 six * (((w_1(i, k) - wm(i, k)) * (w_1(i, k) - wm(i, k)))) * varnce_w_1(i, k) +
                 (((w_1(i, k) - wm(i, k)) * (w_1(i, k) - wm(i, k))) *
                  ((w_1(i, k) - wm(i, k)) * (w_1(i, k) - wm(i, k))))) +
            (one - mixt_frac(i, k)) *
                (three * (varnce_w_2(i, k) * varnce_w_2(i, k)) +
                 six * (((w_2(i, k) - wm(i, k)) * (w_2(i, k) - wm(i, k)))) * varnce_w_2(i, k) +
                 (((w_2(i, k) - wm(i, k)) * (w_2(i, k) - wm(i, k))) *
                  ((w_2(i, k) - wm(i, k)) * (w_2(i, k) - wm(i, k)))));
      });
}

namespace {

void calc_wp2xp2_pdf(int nz, int ngrdcol,                  // In
                     Real2D wm, Real2D xm, Real2D w_1,     // In
                     Real2D w_2, Real2D x_1, Real2D x_2,   // In
                     Real2D varnce_w_1, Real2D varnce_w_2, // In
                     Real2D varnce_x_1, Real2D varnce_x_2, // In
                     Real2D corr_w_x_1, Real2D corr_w_x_2, // In
                     Real2D mixt_frac,                     // In
                     Real2D wp2xp2) {
  // Description:
  // Calculates <w'^2x'^2> by integrating over the PDF of w and x.  The
  // integral
  // is:
  //
  // <w'^2x'^2>
  // = INT(-inf:inf) INT(-inf:inf) ( w - <w> )^2 ( x - <x> )^2 P(w,x) dx dw;
  //
  // where <w> is the overall mean of w, <x> is the overall mean of x, and
  // P(w,x) is a two-component bivariate normal distribution of w and x.  The
  // integrated equation is:
  //
  // <w'^2x'^2>
  //   = mixt_frac
  //      * ( ( mu_w_1 - <w> )**2 * ( ( mu_x_1 - <x> )**2 + sigma_x_1^2 )
  //      + four * corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_x_1 - <x> ) * (
  //      mu_w_1 - <w> )
  //      + ( ( mu_x_1 - <x> )**2 + ( 1 + 2*corr_w_x_1**2 ) * sigma_x_1^2 ) *
  //      sigma_w_1^2 )
  //    + ( one - mixt_frac )
  //      * ( ( mu_w_2 - <w> )**2 * ( ( mu_x_2 - <x> )**2 + sigma_x_2^2 )
  //      + four * corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_x_2 - <x> ) * (
  //      mu_w_2 - <w> )
  //      + ( ( mu_x_2 - <x> )**2 + ( 1 + 2*corr_w_x_2**2 ) * sigma_x_2^2 ) *
  //      sigma_w_2^2 )
  //
  // where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
  // of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
  // component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
  // the standard deviation of w in the 1st PDF component, sigma_w_2 is the
  // standard deviation of w in the 2nd PDF component, sigma_x_1 is the
  // standard deviation of x in the 1st PDF component, sigma_x_2 is the
  // standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
  // correlation of w and x in the 1st PDF component, corr_w_x_2 is the
  // correlation of w and x in the 2nd PDF component, and mixt_frac is the
  // mixture fraction, which is the weight of the 1st PDF component.

  // References:
  //-----------------------------------------------------------------------
  Kokkos::parallel_for(
      "calc_wp2xp2_pdf::moment", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = (cell / ngrdcol);
        // Calculate <w'x'^2> by integrating over the PDF.
        // Native rounds each variance product before fusing the squared deviation.
        wp2xp2(i, k) =
            mixt_frac(i, k) *
                (((w_1(i, k) - wm(i, k)) * (w_1(i, k) - wm(i, k))) *
                     (((x_1(i, k) - xm(i, k)) * (x_1(i, k) - xm(i, k))) + varnce_x_1(i, k)) +
                 four * corr_w_x_1(i, k) * Kokkos::sqrt(varnce_w_1(i, k) * varnce_x_1(i, k)) *
                     (x_1(i, k) - xm(i, k)) * (w_1(i, k) - wm(i, k)) +
                 Kokkos::fma(x_1(i, k) - xm(i, k), x_1(i, k) - xm(i, k),
                             (1 + 2 * (corr_w_x_1(i, k) * corr_w_x_1(i, k))) * varnce_x_1(i, k)) *
                     varnce_w_1(i, k)) +
            (one - mixt_frac(i, k)) *
                (((w_2(i, k) - wm(i, k)) * (w_2(i, k) - wm(i, k))) *
                     (((x_2(i, k) - xm(i, k)) * (x_2(i, k) - xm(i, k))) + varnce_x_2(i, k)) +
                 four * corr_w_x_2(i, k) * Kokkos::sqrt(varnce_w_2(i, k) * varnce_x_2(i, k)) *
                     (x_2(i, k) - xm(i, k)) * (w_2(i, k) - wm(i, k)) +
                 Kokkos::fma(x_2(i, k) - xm(i, k), x_2(i, k) - xm(i, k),
                             (1 + 2 * (corr_w_x_2(i, k) * corr_w_x_2(i, k))) * varnce_x_2(i, k)) *
                     varnce_w_2(i, k));
      });
}

} // namespace

void calc_wp2xp_pdf(int nz, int ngrdcol,                          // In
                    Real2D wm, Real2D xm, Real2D w_1, Real2D w_2, // In
                    Real2D x_1, Real2D x_2,                       // In
                    Real2D varnce_w_1, Real2D varnce_w_2,         // In
                    Real2D varnce_x_1, Real2D varnce_x_2,         // In
                    Real2D corr_w_x_1, Real2D corr_w_x_2,         // In
                    Real2D mixt_frac,                             // In
                    Real2D wp2xp) {
  // Description:
  // Calculates <w'^2 x'> by integrating over the PDF of w and x.  The integral
  // is:
  //
  // <w'^2 x'>
  // = INT(-inf:inf) INT(-inf:inf) ( w - <w> )^2 ( x - <x> ) P(w,x) dx dw;
  //
  // where <w> is the overall mean of w, <x> is the overall mean of x, and
  // P(w,x) is a two-component bivariate normal distribution of w and x.  The
  // integrated equation is:
  //
  // <w'^2 x'>
  // = mixt_frac * ( ( mu_x_1 - <x> ) * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
  //                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
  //                   * ( mu_w_1 - <w> ) )
  //   + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )
  //                           * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 )
  //                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
  //                             * ( mu_w_2 - <w> ) );
  //
  // where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
  // of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
  // component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
  // the standard deviation of w in the 1st PDF component, sigma_w_2 is the
  // standard deviation of w in the 2nd PDF component, sigma_x_1 is the
  // standard deviation of x in the 1st PDF component, sigma_x_2 is the
  // standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
  // correlation of w and x in the 1st PDF component, corr_w_x_2 is the
  // correlation of w and x in the 2nd PDF component, and mixt_frac is the
  // mixture fraction, which is the weight of the 1st PDF component.

  // References:
  //-----------------------------------------------------------------------

  // Arrays are indexed (column, level).
  // wm: Mean of w (overall) [m/s]
  // xm: Mean of x (overall) [units vary]
  // w_1: Mean of w (1st PDF component) [m/s]
  // w_2: Mean of w (2nd PDF component) [m/s]
  // x_1: Mean of x (1st PDF component) [units vary]
  // x_2: Mean of x (2nd PDF component) [units vary]
  // varnce_w_1: Variance of w (1st PDF component) [m^2/s^2]
  // varnce_w_2: Variance of w (2nd PDF component) [m^2/s^2]
  // varnce_x_1: Variance of x (1st PDF component) [(units vary)^2]
  // varnce_x_2: Variance of x (2nd PDF component) [(units vary)^2]
  // corr_w_x_1: Correlation of w and x (1st PDF comp.) [-]
  // corr_w_x_2: Correlation of w and x (2nd PDF comp.) [-]
  // mixt_frac: Mixture fraction [-]
  // wp2xp: <w'^2 x'> [m^2/s^2 (units vary)]

  // Calculate <w'^2 x'> by integrating over the PDF.
  Kokkos::parallel_for(
      "calc_wp2xp_pdf", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        wp2xp(i, k) =
            mixt_frac(i, k) *
                ((((w_1(i, k) - wm(i, k)) * (w_1(i, k) - wm(i, k))) + varnce_w_1(i, k)) *
                     (x_1(i, k) - xm(i, k)) +
                 two * corr_w_x_1(i, k) * Kokkos::sqrt(varnce_w_1(i, k) * varnce_x_1(i, k)) *
                     (w_1(i, k) - wm(i, k))) +
            (one - mixt_frac(i, k)) *
                ((((w_2(i, k) - wm(i, k)) * (w_2(i, k) - wm(i, k))) + varnce_w_2(i, k)) *
                     (x_2(i, k) - xm(i, k)) +
                 two * corr_w_x_2(i, k) * Kokkos::sqrt(varnce_w_2(i, k) * varnce_x_2(i, k)) *
                     (w_2(i, k) - wm(i, k)));
      });
}

void calc_wpxp2_pdf(int nz, int ngrdcol,                  // In
                    Real2D wm, Real2D xm, Real2D w_1,     // In
                    Real2D w_2, Real2D x_1, Real2D x_2,   // In
                    Real2D varnce_w_1, Real2D varnce_w_2, // In
                    Real2D varnce_x_1, Real2D varnce_x_2, // In
                    Real2D corr_w_x_1, Real2D corr_w_x_2, // In
                    Real2D mixt_frac,                     // In
                    Real2D wpxp2) {
  // Description:
  // Calculates <w'x'^2> by integrating over the PDF of w and x.  The integral
  // is:
  //
  // <w'x'^2>
  // = INT(-inf:inf) INT(-inf:inf) ( w - <w> ) ( x - <x> )^2 P(w,x) dx dw;
  //
  // where <w> is the overall mean of w, <x> is the overall mean of x, and
  // P(w,x) is a two-component bivariate normal distribution of w and x.  The
  // integrated equation is:
  //
  // <w'x'^2>
  // = mixt_frac * ( ( mu_w_1 - <w> ) * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
  //                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
  //                   * ( mu_x_1 - <x> ) )
  //   + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )
  //                           * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 )
  //                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
  //                             * ( mu_x_2 - <x> ) );
  //
  // where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
  // of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
  // component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
  // the standard deviation of w in the 1st PDF component, sigma_w_2 is the
  // standard deviation of w in the 2nd PDF component, sigma_x_1 is the
  // standard deviation of x in the 1st PDF component, sigma_x_2 is the
  // standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
  // correlation of w and x in the 1st PDF component, corr_w_x_2 is the
  // correlation of w and x in the 2nd PDF component, and mixt_frac is the
  // mixture fraction, which is the weight of the 1st PDF component.

  // References:
  //-----------------------------------------------------------------------

  // Arrays are indexed (column, level).
  // wm: Mean of w (overall) [m/s]
  // xm: Mean of x (overall) [units vary]
  // w_1: Mean of w (1st PDF component) [m/s]
  // w_2: Mean of w (2nd PDF component) [m/s]
  // x_1: Mean of x (1st PDF component) [units vary]
  // x_2: Mean of x (2nd PDF component) [units vary]
  // varnce_w_1: Variance of w (1st PDF component) [m^2/s^2]
  // varnce_w_2: Variance of w (2nd PDF component) [m^2/s^2]
  // varnce_x_1: Variance of x (1st PDF component) [(units vary)^2]
  // varnce_x_2: Variance of x (2nd PDF component) [(units vary)^2]
  // corr_w_x_1: Correlation of w and x (1st PDF comp.) [-]
  // corr_w_x_2: Correlation of w and x (2nd PDF comp.) [-]
  // mixt_frac: Mixture fraction [-]
  // wpxp2: <w'x'^2> [m/s (units vary)^2]

  // Calculate <w'x'^2> by integrating over the PDF.
  Kokkos::parallel_for(
      "calc_wpxp2_pdf", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        wpxp2(i, k) =
            mixt_frac(i, k) *
                ((w_1(i, k) - wm(i, k)) *
                     (((x_1(i, k) - xm(i, k)) * (x_1(i, k) - xm(i, k))) + varnce_x_1(i, k)) +
                 two * corr_w_x_1(i, k) * Kokkos::sqrt(varnce_w_1(i, k) * varnce_x_1(i, k)) *
                     (x_1(i, k) - xm(i, k))) +
            (one - mixt_frac(i, k)) *
                ((w_2(i, k) - wm(i, k)) *
                     (((x_2(i, k) - xm(i, k)) * (x_2(i, k) - xm(i, k))) + varnce_x_2(i, k)) +
                 two * corr_w_x_2(i, k) * Kokkos::sqrt(varnce_w_2(i, k) * varnce_x_2(i, k)) *
                     (x_2(i, k) - xm(i, k)));
      });
}

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
                     Real2D wpxpyp) {
  // Description:
  // Calculates <w'x'y'> by integrating over the PDF of w, x, and y.  The
  // integral is:
  //
  // <w'x'y'>
  // = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
  //   ( w - <w> ) ( x - <x> ) ( y - <y> ) P(w,x,y) dy dx dw;
  //
  // where <w> is the overall mean of w, <x> is the overall mean of x, <y> is
  // the overall mean of y, and P(w,x,y) is a two-component trivariate normal
  // distribution of w, x, and y.  The integrated equation is:
  //
  // <w'x'y'>
  // = mixt_frac
  //   * ( ( mu_w_1 - <w> ) * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
  //       + corr_x_y_1 * sigma_x_1 * sigma_y_1 * ( mu_w_1 - <w> )
  //       + corr_w_y_1 * sigma_w_1 * sigma_y_1 * ( mu_x_1 - <x> )
  //       + corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_y_1 - <y> ) )
  //   + ( 1 - mixt_frac )
  //     * ( ( mu_w_2 - <w> ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
  //         + corr_x_y_2 * sigma_x_2 * sigma_y_2 * ( mu_w_2 - <w> )
  //         + corr_w_y_2 * sigma_w_2 * sigma_y_2 * ( mu_x_2 - <x> )
  //         + corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_y_2 - <y> ) );
  //
  // where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
  // of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
  // component, mu_x_2 is the mean of x in the 2nd PDF component, mu_y_1 is the
  // mean of y in the 1st PDF component, mu_y_2 is the mean of y in the 2nd PDF
  // component, sigma_w_1 is the standard deviation of w in the 1st PDF
  // component, sigma_w_2 is the standard deviation of w in the 2nd PDF
  // component, sigma_x_1 is the standard deviation of x in the 1st PDF
  // component, sigma_x_2 is the standard deviation of x in the 2nd PDF
  // component, sigma_y_1 is the standard deviation of y in the 1st PDF
  // component, sigma_y_2 is the standard deviation of y in the 2nd PDF
  // component, corr_w_x_1 is the correlation of w and x in the 1st PDF
  // component, corr_w_x_2 is the correlation of w and x in the 2nd PDF
  // component, corr_w_y_1 is the correlation of w and y in the 1st PDF
  // component, corr_w_y_2 is the correlation of w and y in the 2nd PDF
  // component, corr_x_y_1 is the correlation of x and y in the 1st PDF
  // component, corr_x_y_2 is the correlation of x and y in the 2nd PDF
  // component, and mixt_frac is the mixture fraction, which is the weight of
  // the 1st PDF component.

  // References:
  //-----------------------------------------------------------------------

  // Arrays are indexed (column, level).
  // wm: Mean of w (overall) [m/s]
  // xm: Mean of x (overall) [x units]
  // ym: Mean of y (overall) [y units]
  // w_1: Mean of w (1st PDF component) [m/s]
  // w_2: Mean of w (2nd PDF component) [m/s]
  // x_1: Mean of x (1st PDF component) [x units]
  // x_2: Mean of x (2nd PDF component) [x units]
  // y_1: Mean of y (1st PDF component) [y units]
  // y_2: Mean of y (2nd PDF component) [y units]
  // varnce_w_1: Variance of w (1st PDF component) [m^2/s^2]
  // varnce_w_2: Variance of w (2nd PDF component) [m^2/s^2]
  // varnce_x_1: Variance of x (1st PDF component) [(x units)^2]
  // varnce_x_2: Variance of x (2nd PDF component) [(x units)^2]
  // varnce_y_1: Variance of y (1st PDF component) [(y units)^2]
  // varnce_y_2: Variance of y (2nd PDF component) [(y units)^2]
  // corr_w_x_1: Correlation of w and x (1st PDF component) [-]
  // corr_w_x_2: Correlation of w and x (2nd PDF component) [-]
  // corr_w_y_1: Correlation of w and y (1st PDF component) [-]
  // corr_w_y_2: Correlation of w and y (2nd PDF component) [-]
  // corr_x_y_1: Correlation of x and y (1st PDF component) [-]
  // corr_x_y_2: Correlation of x and y (2nd PDF component) [-]
  // mixt_frac: Mixture fraction [-]
  // wpxpyp: <w'x'y'> [m/s (units vary)]

  // Calculate <w'x'y'> by integrating over the PDF.
  Kokkos::parallel_for(
      "calc_wpxpyp_pdf", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        wpxpyp(i, k) =
            mixt_frac(i, k) *
                ((w_1(i, k) - wm(i, k)) * (x_1(i, k) - xm(i, k)) * (y_1(i, k) - ym(i, k)) +
                 corr_x_y_1(i, k) * Kokkos::sqrt(varnce_x_1(i, k) * varnce_y_1(i, k)) *
                     (w_1(i, k) - wm(i, k)) +
                 corr_w_y_1(i, k) * Kokkos::sqrt(varnce_w_1(i, k) * varnce_y_1(i, k)) *
                     (x_1(i, k) - xm(i, k)) +
                 corr_w_x_1(i, k) * Kokkos::sqrt(varnce_w_1(i, k) * varnce_x_1(i, k)) *
                     (y_1(i, k) - ym(i, k))) +
            (one - mixt_frac(i, k)) *
                ((w_2(i, k) - wm(i, k)) * (x_2(i, k) - xm(i, k)) * (y_2(i, k) - ym(i, k)) +
                 corr_x_y_2(i, k) * Kokkos::sqrt(varnce_x_2(i, k) * varnce_y_2(i, k)) *
                     (w_2(i, k) - wm(i, k)) +
                 corr_w_y_2(i, k) * Kokkos::sqrt(varnce_w_2(i, k) * varnce_y_2(i, k)) *
                     (x_2(i, k) - xm(i, k)) +
                 corr_w_x_2(i, k) * Kokkos::sqrt(varnce_w_2(i, k) * varnce_x_2(i, k)) *
                     (y_2(i, k) - ym(i, k)));
      });
}

namespace {

void calc_liquid_cloud_frac_component(int nz, int ngrdcol,               // In
                                      Real2D mean_chi, Real2D stdev_chi, // In
                                      Real2D cloud_frac, Real2D rc) {
  // Description:
  // Calculates the PDF component cloud water mixing ratio, rc_i, and cloud
  // fraction, cloud_frac_i, for the ith PDF component.
  //
  // The equation for cloud water mixing ratio, rc, at any point is:
  //
  // rc = chi * H(chi);
  //
  // and the equation for cloud fraction at a point, fc, is:
  //
  // fc = H(chi);
  //
  // where where extended liquid water mixing ratio, chi, is equal to cloud
  // water mixing ratio, rc, when positive.  When the atmosphere is saturated
  // at this point, cloud water is found, and rc = chi, while fc = 1.
  // Otherwise, clear air is found at this point, and rc = fc = 0.
  //
  // The mean of rc and fc is calculated by integrating over the PDF, such
  // that:
  //
  // <rc> = INT(-inf:inf) chi * H(chi) * P(chi) dchi; and
  //
  // cloud_frac = <fc> = INT(-inf:inf) H(chi) * P(chi) dchi.
  //
  // This can be rewritten as:
  //
  // <rc> = INT(0:inf) chi * P(chi) dchi; and
  //
  // cloud_frac = <fc> = INT(0:inf) P(chi) dchi;
  //
  // and further rewritten as:
  //
  // <rc> = SUM(i=1,N) mixt_frac_i INT(0:inf) chi * P_i(chi) dchi; and
  //
  // cloud_frac = SUM(i=1,N) mixt_frac_i INT(0:inf) P_i(chi) dchi;
  //
  // where N is the number of PDF components.  The equation for mean rc in the
  // ith PDF component is:
  //
  // rc_i = INT(0:inf) chi * P_i(chi) dchi;
  //
  // and the equation for cloud fraction in the ith PDF component is:
  //
  // cloud_frac_i = INT(0:inf) P_i(chi) dchi.
  //
  // The component values are related to the overall values by:
  //
  // <rc> = SUM(i=1,N) mixt_frac_i * rc_i; and
  //
  // cloud_frac = SUM(i=1,N) mixt_frac_i * cloud_frac_i.

  // References:
  //----------------------------------------------------------------------
  constexpr real_t invrs_sqrt_2   = one / sqrt_2;
  constexpr real_t invrs_sqrt_2pi = one / sqrt_2pi;
  Kokkos::parallel_for(
      "calc_liquid_cloud_frac_component::cloud_fraction", ngrdcol * nz,
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = (cell / ngrdcol);
        real_t zeta;
        if ((Kokkos::abs(mean_chi(i, k)) <= eps && stdev_chi(i, k) <= chi_tol) ||
            (mean_chi(i, k) < -max_num_stdevs * stdev_chi(i, k))) {
          // The mean of chi is at saturation and does not vary in the ith PDF component
          cloud_frac(i, k) = zero;
          rc(i, k)         = zero;
        } else if (mean_chi(i, k) > max_num_stdevs * stdev_chi(i, k)) {
          // The mean of chi is multiple standard deviations above the saturation point.
          // Thus, all cloud in the ith PDF component.
          cloud_frac(i, k) = one;
          rc(i, k)         = mean_chi(i, k);
        } else {
          // The mean of chi is within max_num_stdevs of the saturation point.
          // Thus, layer is partly cloudy, requires calculation.
          zeta             = mean_chi(i, k) / stdev_chi(i, k);
          cloud_frac(i, k) = one_half * (one + Kokkos::erf(zeta * invrs_sqrt_2));
          rc(i, k) = mean_chi(i, k) * cloud_frac(i, k) +
                     stdev_chi(i, k) * Kokkos::exp(-one_half * (zeta * zeta)) * invrs_sqrt_2pi;
        }
      });
}

void calc_ice_cloud_frac_component(int nz, int ngrdcol,               // In
                                   Real2D mean_chi, Real2D stdev_chi, // In
                                   Real2D rc_in, Real2D cloud_frac,   // In
                                   Real2D p_in_Pa, Real2D tl,         // In
                                   Real2D rsatl, Real2D crt,          // In
                                   int saturation_formula,            // In
                                   Real2D ice_supersat_frac, Real2D rc) {
  // Description:
  //   A version of the cloud fraction calculation modified to work
  //   for layers that are potentially below freezing. If there are
  //   no below freezing levels, the ice_supersat_frac calculation is
  //   the same as cloud_frac.
  //
  //   For the below freezing levels, the saturation point will be
  //   non-zero, thus we need to calculate chi_at_ice_sat.
  //
  //   The description of the equations are located in the description
  //   of calc_liquid_cloud_frac_component.
  //----------------------------------------------------------------------
  constexpr real_t invrs_sqrt_2   = one / sqrt_2;
  constexpr real_t invrs_sqrt_2pi = one / sqrt_2pi;

  // Warm boxes use the liquid-cloud calculation. Integer OR reduction keeps
  // the native early return when all boxes are above freezing.
  int l_any_below_freezing = 0;
  Kokkos::parallel_reduce(
      "calc_ice_cloud_frac_component::freezing_gate", ngrdcol * nz,
      KOKKOS_LAMBDA(const int cell, int& any_below) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        if (tl(i, k) > T_freeze_K) {
          ice_supersat_frac(i, k) = cloud_frac(i, k);
          rc(i, k)                = rc_in(i, k);
        } else {
          any_below += 1;
        }
      },
      l_any_below_freezing);
  if (!l_any_below_freezing)
    return;

  Kokkos::parallel_for(
      "calc_ice_cloud_frac_component::ice_cloud_fraction", ngrdcol * nz,
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = (cell / ngrdcol);
        real_t zeta, chi_at_ice_sat;
        if (tl(i, k) <= T_freeze_K) {
          // Temperature is freezing, we must compute chi_at_ice_sat and
          // calculate the new cloud_frac_component
          // Calculate ice saturation only for the cold cell that consumes it.
          const real_t rsat_ice = sat_mixrat_ice(p_in_Pa(i, k), tl(i, k), saturation_formula);
          chi_at_ice_sat        = crt(i, k) * (rsat_ice - rsatl(i, k));
          if ((Kokkos::abs(mean_chi(i, k) - chi_at_ice_sat) <= eps && stdev_chi(i, k) <= chi_tol) ||
              (mean_chi(i, k) - chi_at_ice_sat < -max_num_stdevs * stdev_chi(i, k))) {
            // The mean of chi is at saturation and does not vary in the ith PDF component
            ice_supersat_frac(i, k) = zero;
            rc(i, k)                = zero;
          } else if (mean_chi(i, k) - chi_at_ice_sat > max_num_stdevs * stdev_chi(i, k)) {
            // The mean of chi is multiple standard deviations above the saturation point.
            // Thus, all cloud in the ith PDF component.
            ice_supersat_frac(i, k) = one;
            rc(i, k)                = mean_chi(i, k) - chi_at_ice_sat;
          } else {
            // The mean of chi is within max_num_stdevs of the saturation point.
            // Thus, layer is partly cloudy, requires calculation.
            zeta                    = (mean_chi(i, k) - chi_at_ice_sat) / stdev_chi(i, k);
            ice_supersat_frac(i, k) = one_half * (one + Kokkos::erf(zeta * invrs_sqrt_2));
            rc(i, k) = (mean_chi(i, k) - chi_at_ice_sat) * ice_supersat_frac(i, k) +
                       stdev_chi(i, k) * Kokkos::exp(-one_half * (zeta * zeta)) * invrs_sqrt_2pi;
          }
        }
      });
}

void calc_xprcp_component(int nz, int ngrdcol, // In
                          Real2D wm, Real2D rtm, Real2D thlm, Real2D um, Real2D vm,
                          Real2D rcm,                                                  // In
                          Real2D w_i, Real2D rt_i,                                     // In
                          Real2D thl_i, Real2D u_i, Real2D v_i,                        // In
                          Real2D varnce_w_i, Real2D chi_i,                             // In
                          Real2D stdev_chi_i, Real2D stdev_eta_i,                      // In
                          Real2D corr_w_chi_i, Real2D corr_chi_eta_i,                  // In
                          Real2D crt_i, Real2D cthl_i,                                 // In
                          Real2D rc_i, Real2D cloud_frac_i, int iiPDF_type,            // In
                          Real2D wprcp_contrib_comp_i, Real2D wp2rcp_contrib_comp_i,   // Out
                          Real2D rtprcp_contrib_comp_i, Real2D thlprcp_contrib_comp_i, // Out
                          Real2D uprcp_contrib_comp_i, Real2D vprcp_contrib_comp_i) {
  // Description:
  // Calculates the contribution to <w'rc'>, <w'^2 rc'>, <rt'rc'>, and
  // <thl'rc'> from the ith PDF component.
  //
  //
  // <w'rc'>
  // -------
  //
  // The value of <w'rc'> is calculated by integrating over the PDF:
  //
  // <w'rc'>
  // = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
  //   ( w - <w> ) ( rc - <rc> ) P(w,rt,thl) dthl drt dw;
  //
  // where <w> is the overall mean of w, <rc> is the overall mean of rc, and
  // P(w,rt,thl) is a two-component trivariate normal distribution of w, rt,
  // and thl.  This equation is rewritten as:
  //
  // <w'rc'>
  // = mixt_frac
  //   * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
  //     ( w - <w> ) ( rc - <rc> ) P_1(w,rt,thl) dthl drt dw
  //   + ( 1 - mixt_frac )
  //     * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
  //       ( w - <w> ) ( rc - <rc> ) P_2(w,rt,thl) dthl drt dw;
  //
  // where mixt_frac is the mixture fraction, which is the weight of the 1st
  // PDF component, and where P_1(w,rt,thl) and P_2(w,rt,thl) are the equations
  // for the trivariate normal PDF of w, rt, and thl in the 1st and 2nd PDF
  // components, respectively.  The contribution from the ith PDF component is:
  //
  // INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
  // ( w - <w> ) ( rc - <rc> ) P_i(w,rt,thl) dthl drt dw;
  //
  // where P_i(w,rt,thl) is the trivariate normal PDF of w, rt, and thl in the
  // ith PDF component.  The PDF undergoes a PDF transformation in each PDF
  // component, which is a change of variables and a translation, stretching,
  // and rotation of the axes.  The PDF becomes a trivariate normal PDF that is
  // written in terms of w, chi, and eta coordinates.  Cloud water mixing
  // ratio, rc, is written in terms of extended liquid water mixing ratio, chi,
  // such that:
  //
  // rc = chi H(chi);
  //
  // where H(chi) is the Heaviside step function.  The contribution from the
  // ith PDF component to <w'rc'> can be written as:
  //
  // INT(-inf:inf) INT(-inf:inf)
  // ( w - <w> ) ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw;
  //
  // where P_i(w,chi) is the bivariate normal PDF of w and chi in the ith PDF
  // component.  The solved equation for the <w'rc'> contribution from the ith
  // PDF component (wprcp_contrib_comp_i) is:
  //
  // wprcp_contrib_comp_i
  // = INT(-inf:inf) INT(-inf:inf)
  //   ( w - <w> ) ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw
  // = ( mu_w_i - <w> )
  //   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
  //       + 1/sqrt(2*pi) * sigma_chi_i
  //         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
  //   + corr_w_chi_i * sigma_w_i * sigma_chi_i
  //     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
  //
  // where mu_w_i is the mean of w in the ith PDF component, mu_chi_i is the
  // mean of chi in the ith PDF component, sigma_w_i is the standard deviation
  // of w in the ith PDF component, sigma_chi_i is the standard deviation of
  // chi in the ith PDF component, and corr_w_chi_i is the correlation of w and
  // chi in the ith PDF component.
  //
  // Special case:  sigma_chi_i = 0.
  //
  // In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
  // in the ith PDF component.  The equation becomes:
  //
  // wprcp_contrib_comp_i
  // = | ( mu_w_i - <w> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
  //   | ( mu_w_i - <w> ) * ( -<rc> ); when mu_chi_i <= 0.
  //
  //
  // <w'^2 rc'>
  // ----------
  //
  // The value of <w'^2 rc'> is calculated by integrating over the PDF:
  //
  // <w'^2 rc'>
  // = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
  //   ( w - <w> )^2 ( rc - <rc> ) P(w,rt,thl) dthl drt dw.
  //
  // This equation is rewritten as:
  //
  // <w'^2 rc'>
  // = mixt_frac
  //   * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
  //     ( w - <w> )^2 ( rc - <rc> ) P_1(w,rt,thl) dthl drt dw
  //   + ( 1 - mixt_frac )
  //     * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
  //       ( w - <w> )^2 ( rc - <rc> ) P_2(w,rt,thl) dthl drt dw.
  //
  // The contribution from the ith PDF component is:
  //
  // INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
  // ( w - <w> )^2 ( rc - <rc> ) P_i(w,rt,thl) dthl drt dw.
  //
  // The PDF undergoes a PDF transformation in each PDF component, and becomes
  // a trivariate normal PDF that is written in terms of w, chi, and eta
  // coordinates.  The contribution from the ith PDF component to <w'^2 rc'>
  // can be written as:
  //
  // INT(-inf:inf) INT(-inf:inf)
  // ( w - <w> )^2 ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw.
  //
  // The solved equation for the <w'^2 rc'> contribution from the ith PDF
  // component (wp2rcp_contrib_comp_i) is:
  //
  // wp2rcp_contrib_comp_i
  // = INT(-inf:inf) INT(-inf:inf)
  //   ( w - <w> )^2 ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw
  // = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 )
  //   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
  //       + 1/sqrt(2*pi) * sigma_chi_i
  //         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
  //   + ( mu_w_i - <w> ) * corr_w_chi_i * sigma_w_i * sigma_chi_i
  //     * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
  //   + 1/sqrt(2*pi) * corr_w_chi_i^2 * sigma_w_i^2 * sigma_chi_i
  //     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }.
  //
  // Special case:  sigma_chi_i = 0.
  //
  // In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
  // in the ith PDF component.  The equation becomes:
  //
  // wp2rcp_contrib_comp_i
  // = | ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( mu_chi_i - <rc> );
  //   |     when mu_chi_i > 0;
  //   | ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( -<rc> );
  //   |     when mu_chi_i <= 0.
  //
  //
  // <rt'rc'>
  // --------
  //
  // The value of <rt'rc'> is calculated by integrating over the PDF:
  //
  // <rt'rc'>
  // = INT(-inf:inf) INT(-inf:inf)
  //   ( rt - <rt> ) ( rc - <rc> ) P(rt,thl) dthl drt;
  //
  // where <rt> is the overall mean of rt, and where P(rt,thl) is a
  // two-component bivariate normal distribution of rt and thl.  This equation
  // is rewritten as:
  //
  // <rt'rc'>
  // = mixt_frac
  //   * INT(-inf:inf) INT(-inf:inf)
  //     ( rt - <rt> ) ( rc - <rc> ) P_1(rt,thl) dthl drt
  //   + ( 1 - mixt_frac )
  //     * INT(-inf:inf) INT(-inf:inf)
  //       ( rt - <rt> ) ( rc - <rc> ) P_2(rt,thl) dthl drt;
  //
  // where P_1(rt,thl) and P_2(rt,thl) are the equations for the bivariate
  // normal PDF of rt and thl in the 1st and 2nd PDF components, respectively.
  // The contribution from the ith PDF component is:
  //
  // INT(-inf:inf) INT(-inf:inf)
  // ( rt - <rt> ) ( rc - <rc> ) P_i(rt,thl) dthl drt;
  //
  // where P_i(rt,thl) is the bivariate normal PDF of rt and thl in the ith PDF
  // component.  The PDF undergoes a PDF transformation in each PDF component,
  // and becomes a bivariate normal PDF that is written in terms of chi and
  // eta coordinates.  Total water mixing ratio, rt, is rewritten in terms of
  // chi and eta by:
  //
  // rt = mu_rt_i
  //      + ( ( eta - mu_eta_i ) + ( chi - mu_chi_i ) ) / ( 2 * crt_i );
  //
  // where mu_rt_i is the mean of rt in the ith PDF component, mu_eta_i is the
  // mean of eta in the ith PDF component, and crt_i is a coefficient on rt in
  // the chi/eta transformation equations.  The contribution from the ith PDF
  // component to <rt'rc'> can be written as:
  //
  // INT(-inf:inf) INT(-inf:inf)
  // ( mu_rt_i - <rt> + ( eta - mu_eta_i ) / ( 2 * crt_i )
  //   + ( chi - mu_chi_i ) / ( 2 * crt_i ) )
  // * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi;
  //
  // where P_i(chi,eta) is the bivariate normal PDF of chi and eta in the ith
  // PDF component.  The solved equation for the <rt'rc'> contribution from the
  // ith PDF component (rtprcp_contrib_comp_i) is:
  //
  // rtprcp_contrib_comp_i
  // = INT(-inf:inf) INT(-inf:inf)
  //   ( mu_rt_i - <rt> + ( eta - mu_eta_i ) / ( 2 * crt_i )
  //     + ( chi - mu_chi_i ) / ( 2 * crt_i ) )
  //   * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi
  // = ( mu_rt_i - <rt> )
  //   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
  //       + 1/sqrt(2*pi) * sigma_chi_i
  //         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
  //   + ( corr_chi_eta_i * sigma_eta_i + sigma_chi_i ) / ( 2 * crt_i )
  //     * sigma_chi_i
  //     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
  //
  // where sigma_eta_i is the standard deviation of eta in the ith PDF
  // component and corr_chi_eta_i is the correlation of chi and eta in the ith
  // PDF component.
  //
  // Special case:  sigma_chi_i = 0.
  //
  // In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
  // in the ith PDF component.  The equation becomes:
  //
  // rtprcp_contrib_comp_i
  // = | ( mu_rt_i - <rt> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
  //   | ( mu_rt_i - <rt> ) * ( -<rc> ); when mu_chi_i <= 0.
  //
  //
  // <thl'rc'>
  // ---------
  //
  // The value of <thl'rc'> is calculated by integrating over the PDF:
  //
  // <thl'rc'>
  // = INT(-inf:inf) INT(-inf:inf)
  //   ( thl - <thl> ) ( rc - <rc> ) P(rt,thl) dthl drt;
  //
  // where <thl> is the overall mean of thl.  This equation is rewritten as:
  //
  // <thl'rc'>
  // = mixt_frac
  //   * INT(-inf:inf) INT(-inf:inf)
  //     ( thl - <thl> ) ( rc - <rc> ) P_1(rt,thl) dthl drt
  //   + ( 1 - mixt_frac )
  //     * INT(-inf:inf) INT(-inf:inf)
  //       ( thl - <thl> ) ( rc - <rc> ) P_2(rt,thl) dthl drt.
  //
  // The contribution from the ith PDF component is:
  //
  // INT(-inf:inf) INT(-inf:inf)
  // ( thl - <thl> ) ( rc - <rc> ) P_i(rt,thl) dthl drt.
  //
  // The PDF undergoes a PDF transformation in each PDF component, and becomes
  // a bivariate normal PDF that is written in terms of chi and eta
  // coordinates.  Liquid water potential temperature, thl, is rewritten in
  // terms of chi and eta by:
  //
  // thl = mu_thl_i
  //       + ( ( eta - mu_eta_i ) - ( chi - mu_chi_i ) ) / ( 2 * cthl_i );
  //
  // where mu_thl_i is the mean of thl in the ith PDF component and cthl_i is a
  // coefficient on thl in the chi/eta transformation equations.  The
  // contribution from the ith PDF component to <thl'rc'> can be written as:
  //
  // INT(-inf:inf) INT(-inf:inf)
  // ( mu_thl_i - <thl> + ( eta - mu_eta_i ) / ( 2 * cthl_i )
  //   - ( chi - mu_chi_i ) / ( 2 * cthl_i ) )
  // * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi.
  //
  // The solved equation for the <thl'rc'> contribution from the ith PDF
  // component (thlprcp_contrib_comp_i) is:
  //
  // thlprcp_contrib_comp_i
  // = INT(-inf:inf) INT(-inf:inf)
  //   ( mu_thl_i - <thl> + ( eta - mu_eta_i ) / ( 2 * cthl_i )
  //     - ( chi - mu_chi_i ) / ( 2 * cthl_i ) )
  //   * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi
  // = ( mu_thl_i - <thl> )
  //   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
  //       + 1/sqrt(2*pi) * sigma_chi_i
  //         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
  //   + ( corr_chi_eta_i * sigma_eta_i - sigma_chi_i ) / ( 2 * cthl_i )
  //     * sigma_chi_i
  //     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
  //
  // Special case:  sigma_chi_i = 0.
  //
  // In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
  // in the ith PDF component.  The equation becomes:
  //
  // thlprcp_contrib_comp_i
  // = | ( mu_thl_i - <thl> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
  //   | ( mu_thl_i - <thl> ) * ( -<rc> ); when mu_chi_i <= 0.
  //
  //
  // Use equations for PDF component cloud fraction cloud water mixing ratio
  // -----------------------------------------------------------------------
  //
  // The equation for cloud fraction in the ith PDF component, fc_i, is:
  //
  // fc_i = 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
  //
  // In the special case that sigma_chi_i = 0, the equation becomes:
  //
  // fc_i = | 1; when mu_chi_i > 0;
  //        | 0; when mu_chi_i <= 0.
  //
  // The equation for mean cloud water mixing ratio in the ith PDF component,
  // rc_i, is:
  //
  // rc_i
  // = mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
  //   + 1/sqrt(2*pi) * sigma_chi_i
  //     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }
  // = mu_chi_i * fc_i
  //   + 1/sqrt(2*pi) * sigma_chi_i
  //     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }.
  //
  // In the special case that sigma_chi_i = 0, the equation becomes:
  //
  // rc_i = | mu_chi_i; when mu_chi_i > 0;
  //        | 0; when mu_chi_i <= 0.
  //
  // The above equations can be substituted into the equations for
  // wprcp_contrib_comp_i, wp2rcp_contrib_comp_i, rtprcp_contrib_comp_i, and
  // thlprcp_contrib_comp_i.  The new equations are:
  //
  // wprcp_contrib_comp_i
  // = ( mu_w_i - <w> ) * ( rc_i - <rc> )
  //   + corr_w_chi_i * sigma_w_i * sigma_chi_i * fc_i;
  //
  // wp2rcp_contrib_comp_i
  // = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( rc_i - <rc> )
  //   + 2 * ( mu_w_i - <w> ) * corr_w_chi_i * sigma_w_i * sigma_chi_i * fc_i
  //   + 1/sqrt(2*pi) * corr_w_chi_i^2 * sigma_w_i^2 * sigma_chi_i
  //     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) };
  //
  // rtprcp_contrib_comp_i
  // = ( mu_rt_i - <rt> ) * ( rc_i - <rc> )
  //   + ( corr_chi_eta_i * sigma_eta_i + sigma_chi_i ) / ( 2 * crt_i )
  //     * sigma_chi_i * fc_i; and
  //
  // thlprcp_contrib_comp_i
  // = ( mu_thl_i - <thl> ) * ( rc_i - <rc> )
  //   + ( corr_chi_eta_i * sigma_eta_i - sigma_chi_i ) / ( 2 * cthl_i )
  //     * sigma_chi_i * fc_i.
  //
  // While the above equations reduce to their listed versions in the special
  // case that sigma_chi_i = 0, those versions are faster to calculate.  When
  // mu_chi_i > 0, they are:
  //
  // wprcp_contrib_comp_i = ( mu_w_i - <w> ) * ( mu_chi_i - <rc> );
  // wp2rcp_contrib_comp_i
  // = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( mu_chi_i - <rc> );
  // rtprcp_contrib_comp_i = ( mu_rt_i - <rt> ) * ( mu_chi_i - <rc> ); and
  // thlprcp_contrib_comp_i = ( mu_thl_i - <thl> ) * ( mu_chi_i - <rc> );
  //
  // and when mu_chi_i <= 0, they are:
  //
  // wprcp_contrib_comp_i = - ( mu_w_i - <w> ) * <rc>;
  // wp2rcp_contrib_comp_i = - ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * <rc>;
  // rtprcp_contrib_comp_i = - ( mu_rt_i - <rt> ) * <rc>; and
  // thlprcp_contrib_comp_i = - ( mu_thl_i - <thl> ) * <rc>.

  // References:
  //-----------------------------------------------------------------------
  // Changing these conditionals may result in inconsistencies with the
  // conditional statements located in calc_cloud_frac_component.
  Kokkos::parallel_for(
      "calc_xprcp_component::component_moments", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i                = cell % ngrdcol;
        const int k                = (cell / ngrdcol);
        wprcp_contrib_comp_i(i, k) = (w_i(i, k) - wm(i, k)) * (rc_i(i, k) - rcm(i, k));
        wp2rcp_contrib_comp_i(i, k) =
            (((w_i(i, k) - wm(i, k)) * (w_i(i, k) - wm(i, k))) + varnce_w_i(i, k)) *
            (rc_i(i, k) - rcm(i, k));
        rtprcp_contrib_comp_i(i, k) =
            (rt_i(i, k) - rtm(i, k)) * (rc_i(i, k) - rcm(i, k)) +
            (corr_chi_eta_i(i, k) * stdev_eta_i(i, k) + stdev_chi_i(i, k)) / (two * crt_i(i, k)) *
                stdev_chi_i(i, k) * cloud_frac_i(i, k);
        thlprcp_contrib_comp_i(i, k) =
            (thl_i(i, k) - thlm(i, k)) * (rc_i(i, k) - rcm(i, k)) +
            (corr_chi_eta_i(i, k) * stdev_eta_i(i, k) - stdev_chi_i(i, k)) / (two * cthl_i(i, k)) *
                stdev_chi_i(i, k) * cloud_frac_i(i, k);
        uprcp_contrib_comp_i(i, k) = (u_i(i, k) - um(i, k)) * (rc_i(i, k) - rcm(i, k));
        vprcp_contrib_comp_i(i, k) = (v_i(i, k) - vm(i, k)) * (rc_i(i, k) - rcm(i, k));
      });
  // Non-ADG components may have corr_w_chi_i /= 0. In principle,
  // uprcp_contrib_comp_i might also depend on corr_u_w_i here.
  if (!(iiPDF_type == iiPDF_ADG1 || iiPDF_type == iiPDF_ADG2 || iiPDF_type == iiPDF_new_hybrid)) {
    Kokkos::parallel_for(
        "calc_xprcp_component::within_component_correlation", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = (cell / ngrdcol);
          // Chi varies significantly in the ith PDF component (stdev_chi > chi_tol)
          // and there is some cloud (0 < cloud_frac <= 1)
          if (stdev_chi_i(i, k) > chi_tol && cloud_frac_i(i, k) > zero) {
            wprcp_contrib_comp_i(i, k) =
                wprcp_contrib_comp_i(i, k) + corr_w_chi_i(i, k) * Kokkos::sqrt(varnce_w_i(i, k)) *
                                                 stdev_chi_i(i, k) * cloud_frac_i(i, k);
            wp2rcp_contrib_comp_i(i, k) =
                wp2rcp_contrib_comp_i(i, k) +
                two * (w_i(i, k) - wm(i, k)) * corr_w_chi_i(i, k) * Kokkos::sqrt(varnce_w_i(i, k)) *
                    stdev_chi_i(i, k) * cloud_frac_i(i, k) +
                (corr_w_chi_i(i, k) * corr_w_chi_i(i, k)) * varnce_w_i(i, k) * stdev_chi_i(i, k) *
                    Kokkos::exp(-(chi_i(i, k) * chi_i(i, k)) /
                                (two * (stdev_chi_i(i, k) * stdev_chi_i(i, k)))) /
                    sqrt_2pi;
            // In principle, uprcp_contrib_comp_i might depend on corr_u_w_i here.
          }
        });
  }
}

} // namespace

void calc_w_up_in_cloud(int nz, int ngrdcol,                                          // In
                        Real2D mixt_frac, Real2D cloud_frac_1, Real2D cloud_frac_2,   // In
                        Real2D w_1, Real2D w_2, Real2D varnce_w_1, Real2D varnce_w_2, // In
                        Real2D w_up_in_cloud, Real2D w_down_in_cloud,                 // Out
                        Real2D cloudy_updraft_frac, Real2D cloudy_downdraft_frac) {   // Out
  // Description:
  // Subroutine that computes the mean cloudy updraft (and also calculates
  // the mean cloudy downdraft).
  //
  // In order to activate aerosol, we'd like to feed the activation scheme
  // a vertical velocity that's representative of cloudy updrafts. For skewed
  // layers, like cumulus layers, this might be an improvement over the square
  // root of wp2 that's currently used. At the same time, it would be simpler
  // and less expensive than feeding SILHS samples into the aerosol code
  // (see larson-group/e3sm#19 and larson-group/e3sm#26).
  //
  // The formulas are only valid for certain PDFs in CLUBB (ADG1, ADG2,
  // new hybrid), hence we omit calculation if another PDF type is used.
  //
  // References: https://www.overleaf.com/project/614a136d47846639af22ae34
  //----------------------------------------------------------------------

  // nz: Number of vertical levels
  // ngrdcol: Number of grid columns
  //----------- Input Variables -----------
  // mixt_frac: mixture fraction [-]
  // cloud_frac_1: cloud fraction (1st PDF component) [-]
  // cloud_frac_2: cloud fraction (2nd PDF component) [-]
  // w_1: upward velocity (1st PDF component) [m/s]
  // w_2: upward velocity (2nd PDF component) [m/s]
  // varnce_w_1: variance of w (1st PDF component) [m^2/s^2]
  // varnce_w_2: variance of w (2nd PDF component) [m^2/s^2]

  //----------- Output Variables -----------
  // w_up_in_cloud: mean cloudy updraft speed [m/s]
  // w_down_in_cloud: mean cloudy downdraft speed [m/s]
  // cloudy_updraft_frac: cloudy updraft fraction [-]
  // cloudy_downdraft_frac: cloudy downdraft fraction [-]

  Kokkos::parallel_for(
      "calc_w_up_in_cloud", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;

        //----------- Local Variables -----------
        real_t w_up_1, w_up_2;        // integral of w and Heaviside fnc, where w > 0
        real_t w_down_1, w_down_2;    // integral of w and Heaviside fnc, where w < 0
        real_t stdev_w_1, stdev_w_2;  // Standard deviation of w
        real_t ratio_w_1;             // mu_w_1 / (sqrt(2) * sigma_w_1)
        real_t ratio_w_2;             // mu_w_2 / (sqrt(2) * sigma_w_2)
        real_t erf_ratio_w_1;         // erf(ratio_w_1)
        real_t erf_ratio_w_2;         // erf(ratio_w_2)
        real_t exp_neg_ratio_w_1_sqd; // exp(-ratio_w_1^2)
        real_t exp_neg_ratio_w_2_sqd; // exp(-ratio_w_2^2)
        real_t updraft_frac_1;        // Fraction of 1st PDF comp. where w > 0
        real_t updraft_frac_2;        // Fraction of 2nd PDF comp. where w > 0
        real_t downdraft_frac_1;      // Fraction of 1st PDF comp. where w < 0
        real_t downdraft_frac_2;      // Fraction of 2nd PDF comp. where w < 0

        stdev_w_1 = Kokkos::sqrt(varnce_w_1(i, k));
        stdev_w_2 = Kokkos::sqrt(varnce_w_2(i, k));

        // Calculate quantities in the 1st PDF component.
        if (w_1(i, k) > max_num_stdevs * stdev_w_1) {

          // The mean of w in the 1st PDF component is more than
          // max_num_stdevs standard deviations above 0.
          // The entire 1st PDF component is found in an updraft (w > 0).
          w_up_1           = w_1(i, k);
          updraft_frac_1   = one;
          w_down_1         = zero;
          downdraft_frac_1 = zero;

        } else if (w_1(i, k) < -max_num_stdevs * stdev_w_1) {

          // The mean of w in the 1st PDF component is more than
          // max_num_stdevs standard deviations below 0.
          // The entire 1st PDF component is found in a downdraft (w < 0).
          w_up_1           = zero;
          updraft_frac_1   = zero;
          w_down_1         = w_1(i, k);
          downdraft_frac_1 = one;

        } else {

          // The 1st PDF component contains both updraft and downdraft.
          ratio_w_1             = w_1(i, k) / (sqrt_2 * (eps > stdev_w_1 ? eps : stdev_w_1));
          erf_ratio_w_1         = Kokkos::erf(ratio_w_1);
          exp_neg_ratio_w_1_sqd = Kokkos::exp(-(ratio_w_1 * ratio_w_1));

          w_up_1 = one_half * w_1(i, k) * (one + erf_ratio_w_1) +
                   (stdev_w_1 / sqrt_2pi) * exp_neg_ratio_w_1_sqd;

          updraft_frac_1 = one_half * (one + erf_ratio_w_1);

          w_down_1 = one_half * w_1(i, k) * (one - erf_ratio_w_1) -
                     (stdev_w_1 / sqrt_2pi) * exp_neg_ratio_w_1_sqd;

          // downdraft_frac_1 = one_half * (one - erf_ratio_w_1)
          downdraft_frac_1 = one - updraft_frac_1;
        }

        // Calculate quantities in the 2nd PDF component.
        if (w_2(i, k) > max_num_stdevs * stdev_w_2) {

          // The mean of w in the 2nd PDF component is more than
          // max_num_stdevs standard deviations above 0.
          // The entire 2nd PDF component is found in an updraft (w > 0).
          w_up_2           = w_2(i, k);
          updraft_frac_2   = one;
          w_down_2         = zero;
          downdraft_frac_2 = zero;

        } else if (w_2(i, k) < -max_num_stdevs * stdev_w_2) {

          // The mean of w in the 2nd PDF component is more than
          // max_num_stdevs standard deviations below 0.
          // The entire 2nd PDF component is found in a downdraft (w < 0).
          w_up_2           = zero;
          updraft_frac_2   = zero;
          w_down_2         = w_2(i, k);
          downdraft_frac_2 = one;

        } else {

          // The 2nd PDF component contains both updraft and downdraft.
          ratio_w_2             = w_2(i, k) / (sqrt_2 * (eps > stdev_w_2 ? eps : stdev_w_2));
          erf_ratio_w_2         = Kokkos::erf(ratio_w_2);
          exp_neg_ratio_w_2_sqd = Kokkos::exp(-(ratio_w_2 * ratio_w_2));

          w_up_2 = one_half * w_2(i, k) * (one + erf_ratio_w_2) +
                   (stdev_w_2 / sqrt_2pi) * exp_neg_ratio_w_2_sqd;

          updraft_frac_2 = one_half * (one + erf_ratio_w_2);

          w_down_2 = one_half * w_2(i, k) * (one - erf_ratio_w_2) -
                     (stdev_w_2 / sqrt_2pi) * exp_neg_ratio_w_2_sqd;

          // downdraft_frac_2 = one_half * (one - erf_ratio_w_2)
          downdraft_frac_2 = one - updraft_frac_2;
        }

        // Calculate the total cloudy updraft fraction.
        cloudy_updraft_frac(i, k) = mixt_frac(i, k) * cloud_frac_1(i, k) * updraft_frac_1 +
                                    (one - mixt_frac(i, k)) * cloud_frac_2(i, k) * updraft_frac_2;

        // Calculate the total cloudy downdraft fraction.
        cloudy_downdraft_frac(i, k) =
            mixt_frac(i, k) * cloud_frac_1(i, k) * downdraft_frac_1 +
            (one - mixt_frac(i, k)) * cloud_frac_2(i, k) * downdraft_frac_2;

        // Calculate the mean vertical velocity found in a cloudy updraft.
        w_up_in_cloud(i, k) = (mixt_frac(i, k) * cloud_frac_1(i, k) * w_up_1 +
                               (one - mixt_frac(i, k)) * cloud_frac_2(i, k) * w_up_2) /
                              (eps > cloudy_updraft_frac(i, k) ? eps : cloudy_updraft_frac(i, k));

        // Calculate the mean vertical velocity found in a cloudy downdraft.
        w_down_in_cloud(i, k) =
            (mixt_frac(i, k) * cloud_frac_1(i, k) * w_down_1 +
             (one - mixt_frac(i, k)) * cloud_frac_2(i, k) * w_down_2) /
            (eps > cloudy_downdraft_frac(i, k) ? eps : cloudy_downdraft_frac(i, k));
      });

  return;
}

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
                        Real3D wpsclrprtp, Real3D wpsclrpthlp) {
  // Interpolate moments to the PDF grids, apply the selected closure, then
  // return cloud and turbulent moments on their native grids. Sampled
  // diagnostics are sent to the C++ stats context supplied by the caller.

  //------------------------------- Local Variables -------------------------------
  auto workspace_scope = workspace.scope();
  Real2D wp2_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::wp2_zt", ngrdcol,
                                   nzt); // wp2 interpolated to thermodynamic levels   [m^2/s^2]
  Real2D rtp2_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::rtp2_zt", ngrdcol,
                                   nzt); // rtp2 interpolated to thermodynamic levels  [kg^2/kg^2]
  Real2D thlp2_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::thlp2_zt", ngrdcol,
                                   nzt); // thlp2 interpolated to thermodynamic levels [K^2]
  Real2D wprtp_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::wprtp_zt", ngrdcol,
                                   nzt); // wprtp interpolated to thermodynamic levels [m/s kg/kg]
  Real2D wpthlp_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::wpthlp_zt", ngrdcol,
                                   nzt); // wpthlp interpolated to thermodynamic levs. [m/s K]
  Real2D rtpthlp_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::rtpthlp_zt", ngrdcol,
                                   nzt); // rtpthlp interp. to thermodynamic levels    [kg/kg K]
  Real2D up2_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::up2_zt", ngrdcol,
                                   nzt); // up2 interpolated to thermodynamic levels   [m^2/s^2]
  Real2D vp2_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::vp2_zt", ngrdcol,
                                   nzt); // vp2 interpolated to thermodynamic levels   [m^2/s^2]
  Real2D upwp_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::upwp_zt", ngrdcol,
                                   nzt); // upwp interpolated to thermodynamic levels  [m^2/s^2]
  Real2D vpwp_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::vpwp_zt", ngrdcol,
                                   nzt); // vpwp interpolated to thermodynamic levels  [m^2/s^2]
  Real2D sigma_sqd_w_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::sigma_sqd_w_zt", ngrdcol,
                                   nzt); // PDF width parameter (thermodynamic levels) [-]
  Real2D Skw_zt = workspace_scope.view<Real2D>(
      "pdf_closure_driver::Skw_zt", ngrdcol, nzt); // Skewness of w on thermodynamic levels      [-]
  Real2D Skrt_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::Skrt_zt", ngrdcol,
                                   nzt); // Skewness of rt on thermodynamic levels     [-]
  Real2D Skthl_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::Skthl_zt", ngrdcol,
                                   nzt); // Skewness of thl on thermodynamic levels    [-]
  Real2D Sku_zt = workspace_scope.view<Real2D>(
      "pdf_closure_driver::Sku_zt", ngrdcol, nzt); // Skewness of u on thermodynamic levels      [-]
  Real2D Skv_zt = workspace_scope.view<Real2D>(
      "pdf_closure_driver::Skv_zt", ngrdcol, nzt); // Skewness of v on thermodynamic levels      [-]
  Real2D wp3_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::wp3_zm", ngrdcol,
                                   nzm); // wp3 interpolated to momentum levels        [m^3/s^3]
  Real2D rtp3_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::rtp3_zm", ngrdcol,
                                   nzm); // rtp3 interpolated to momentum levels       [kg^3/kg^3]
  Real2D thlp3_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::thlp3_zm", ngrdcol,
                                   nzm); // thlp3 interpolated to momentum levels      [K^3]
  Real2D Skw_zm = workspace_scope.view<Real2D>(
      "pdf_closure_driver::Skw_zm", ngrdcol, nzm); // Skewness of w on momentum levels           [-]
  Real2D Skrt_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::Skrt_zm", ngrdcol,
                                   nzm); // Skewness of rt on momentum levels          [-]
  Real2D Skthl_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::Skthl_zm", ngrdcol,
                                   nzm); // Skewness of thl on momentum levels         [-]
  Real3D wphydrometp_zt = workspace_scope.view<Real3D>(
      "pdf_closure_driver::wphydrometp_zt", ngrdcol, nzt,
      hydromet_dim); // Covariance of w and hm (on t-levs.) [(m/s) <hm units>]
  Real3D wpsclrp_zt =
      workspace_scope.view<Real3D>("pdf_closure_driver::wpsclrp_zt", ngrdcol, nzt,
                                   sclr_dim); // w' sclr' interpolated to thermo. levels
  Real3D sclrp2_zt =
      workspace_scope.view<Real3D>("pdf_closure_driver::sclrp2_zt", ngrdcol, nzt,
                                   sclr_dim); // sclr'^2 interpolated to thermo. levels
  Real3D sclrprtp_zt =
      workspace_scope.view<Real3D>("pdf_closure_driver::sclrprtp_zt", ngrdcol, nzt,
                                   sclr_dim); // sclr' r_t' interpolated to thermo. levels
  Real3D sclrpthlp_zt =
      workspace_scope.view<Real3D>("pdf_closure_driver::sclrpthlp_zt", ngrdcol, nzt,
                                   sclr_dim); // sclr' th_l' interpolated thermo. levels
  Real3D Sksclr_zt =
      workspace_scope.view<Real3D>("pdf_closure_driver::Sksclr_zt", ngrdcol, nzt,
                                   sclr_dim); // Skewness of sclr on thermodynamic levels      [-]
  Real2D wp2up2_zt = workspace_scope.view<Real2D>(
      "pdf_closure_driver::wp2up2_zt", ngrdcol, nzt); // w'^2u'^2 (on thermo. grid)       [m^4/s^4]
  Real2D wp2vp2_zt = workspace_scope.view<Real2D>(
      "pdf_closure_driver::wp2vp2_zt", ngrdcol, nzt); // w'^2v'^2 (on thermo. grid)       [m^4/s^4]
  Real2D wp4_zt = workspace_scope.view<Real2D>("pdf_closure_driver::wp4_zt", ngrdcol,
                                               nzt); // w'^4 (on thermo. grid)           [m^4/s^4]
  Real2D wpthvp_zt = workspace_scope.view<Real2D>(
      "pdf_closure_driver::wpthvp_zt", ngrdcol, nzt); // Buoyancy flux (on thermo. grid)  [(K m)/s]
  Real2D rtpthvp_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::rtpthvp_zt", ngrdcol,
                                   nzt); // r_t' th_v' (on thermo. grid)     [(kg K)/kg]
  Real2D thlpthvp_zt = workspace_scope.view<Real2D>("pdf_closure_driver::thlpthvp_zt", ngrdcol,
                                                    nzt); // th_l' th_v' (on thermo. grid)    [K^2]
  Real2D wprcp_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::wprcp_zt", ngrdcol,
                                   nzt); // w' r_c' (on thermo. grid)        [(m kg)/(s kg)]
  Real2D wp2rcp =
      workspace_scope.view<Real2D>("pdf_closure_driver::wp2rcp", ngrdcol,
                                   nzt); // w'^2 r_c' (on thermo. grid)     [(m^2 kg)/(s^2 kg)]
  Real2D rtprcp_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::rtprcp_zt", ngrdcol,
                                   nzt); // r_t' r_c' (on thermo. grid)      [(kg^2)/(kg^2)]
  Real2D thlprcp_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::thlprcp_zt", ngrdcol,
                                   nzt); // th_l' r_c' (on thermo. grid)     [(K kg)/kg]
  Real2D uprcp_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::uprcp_zt", ngrdcol,
                                   nzt); // u' r_c' (on thermo. grid)        [(m kg)/(s kg)]
  Real2D vprcp_zt =
      workspace_scope.view<Real2D>("pdf_closure_driver::vprcp_zt", ngrdcol,
                                   nzt); // v' r_c' (on thermo. grid)        [(m kg)/(s kg)]
  Real2D rc_coef =
      workspace_scope.view<Real2D>("pdf_closure_driver::rc_coef", ngrdcol,
                                   nzt); // Coefficient of X'r_c' (thermo. levs.)  [K/(kg/kg)]
  Real2D rcm_supersat_adj =
      workspace_scope.view<Real2D>("pdf_closure_driver::rcm_supersat_adj", ngrdcol,
                                   nzt); // Adjust. to rcm due to spurious supersat. [kg/kg]
  Real3D sclrpthvp_zt = workspace_scope.view<Real3D>("pdf_closure_driver::sclrpthvp_zt", ngrdcol,
                                                     nzt, sclr_dim); // sclr'th_v' (on thermo. grid)
  Real3D sclrprcp_zt = workspace_scope.view<Real3D>("pdf_closure_driver::sclrprcp_zt", ngrdcol, nzt,
                                                    sclr_dim); // sclr'rc' (on thermo. grid)
  Real2D wprtp2_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::wprtp2_zm", ngrdcol,
                                   nzm); // < w' r_t'^2 > on momentum levels      [m/s kg^2/kg^2]
  Real2D wp2rtp_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::wp2rtp_zm", ngrdcol,
                                   nzm); // < w'^2 r_t' > on momentum levels      [m^2/s^2 kg/kg]
  Real2D wpthlp2_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::wpthlp2_zm", ngrdcol,
                                   nzm); // < w' th_l'^2 > on momentum levels     [m/s K^2]
  Real2D wp2thlp_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::wp2thlp_zm", ngrdcol,
                                   nzm); // < w'^2 th_l' > on momentum levels     [m^2/s^2 K]
  Real2D wprtpthlp_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::wprtpthlp_zm", ngrdcol,
                                   nzm); // < w' r_t' th_l' > on momentum levels  [m/s kg/kg K]
  Real2D wp2thvp_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::wp2thvp_zm", ngrdcol,
                                   nzm); // < w'^2 th_v' > on momentum levels     [m^2/s^2 K]
  Real2D wp2up_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::wp2up_zm", ngrdcol,
                                   nzm); // < w'^2 u' > on momentum levels        [m^3/s^3]
  Real2D wp2rcp_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::wp2rcp_zm", ngrdcol,
                                   nzm); // < w'^2 r_c' > on momentum levles      [m^2/s^2 kg/kg]
  Real2D Skw_velocity =
      workspace_scope.view<Real2D>("pdf_closure_driver::Skw_velocity", ngrdcol,
                                   nzm); // Skewness velocity                        [m/s]
  Real2D cloud_frac_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::cloud_frac_zm", ngrdcol,
                                   nzm); // Cloud Fraction on momentum levels        [-]
  Real2D ice_supersat_frac_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::ice_supersat_frac_zm", ngrdcol,
                                   nzm); // Ice supersat. frac. on momentum levels   [-]
  Real2D rtm_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::rtm_zm", ngrdcol,
                                   nzm); // Total water mixing ratio at mom. levs.   [kg/kg]
  Real2D thlm_zm = workspace_scope.view<Real2D>(
      "pdf_closure_driver::thlm_zm", ngrdcol, nzm); // Liquid water pot. temp. at mom. levs.    [K]
  Real2D rcm_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver::rcm_zm", ngrdcol,
                                   nzm); // rcm at momentum levels                   [kg/kg]
  Real2D rtprcp =
      workspace_scope.view<Real2D>("pdf_closure_driver::rtprcp", ngrdcol,
                                   nzm); // < r_t' r_c' > (momentum levels)          [kg^2/kg^2]
  Real2D rcp2 =
      workspace_scope.view<Real2D>("pdf_closure_driver::rcp2", ngrdcol,
                                   nzm); // Variance of r_c (momentum levels)        [kg^2/kg^2]
  Real3D wpsclrprtp_zm = workspace_scope.view<Real3D>("pdf_closure_driver::wpsclrprtp_zm", ngrdcol,
                                                      nzm, sclr_dim); // w'sclr'rt' on momentum grid
  Real3D wpsclrp2_zm = workspace_scope.view<Real3D>("pdf_closure_driver::wpsclrp2_zm", ngrdcol, nzm,
                                                    sclr_dim); // w'sclr'^2 on momentum grid
  Real3D wpsclrpthlp_zm = workspace_scope.view<Real3D>(
      "pdf_closure_driver::wpsclrpthlp_zm", ngrdcol, nzm, sclr_dim); // w'sclr'thl' on momentum grid
  Real3D wp2sclrp_zm = workspace_scope.view<Real3D>("pdf_closure_driver::wp2sclrp_zm", ngrdcol, nzm,
                                                    sclr_dim); // w'^2 sclr' on momentum grid
  Real3D sclrprcp    = workspace_scope.view<Real3D>("pdf_closure_driver::sclrprcp", ngrdcol, nzm,
                                                    sclr_dim); // < sclr' r_c' > on momentum grid
  Real2D rsat = workspace_scope.view<Real2D>("pdf_closure_driver::rsat", ngrdcol,
                                             nzt); // Saturation mixing ratio from mean rt and thl.
  Real2D rel_humidity = workspace_scope.view<Real2D>(
      "pdf_closure_driver::rel_humidity", ngrdcol, nzt); // Relative humidity after PDF closure [-]
  Real2D T_in_K = workspace_scope.view<Real2D>("pdf_closure_driver::T_in_K", ngrdcol,
                                               nzt); // Absolute temperature [K]
  Real2D chi    = workspace_scope.view<Real2D>("pdf_closure_driver::chi", ngrdcol,
                                               nzt); // Mellor's 's'         [kg/kg]
  Real2D chip2  = workspace_scope.view<Real2D>("pdf_closure_driver::chip2", ngrdcol,
                                               nzt); // Variance of Mellor's 's' [kg/kg]

  const auto zt     = gr.zt;
  const int k_lb_zm = gr.k_lb_zm;
  const int k_ub_zm = gr.k_ub_zm;
  // Scalar tolerances are read by the host between independent species calls.
  const auto sclr_tol_host = to_host_space(workspace, sclr_tol);

  //-------------------------------------- Begin Code --------------------------------------

  //---------------------------------------------------------------------------
  // Interpolate wp3 to momentum levels, rtp3 and thlp3 to momentum levels
  // as needed, and
  // wp2, rtp2, thlp2, up2, vp2, and sclrp2 to thermodynamic levels, and then
  // compute Skw for both the momentum and thermodynamic grid levels. Skrt and
  // Skthl are computed for thermodynamic grid levels and for momentum levels
  // as needed. Sku, Skv, and Sksclr are computed only for thermodynamic levels here.
  //---------------------------------------------------------------------------
  zm2zt_api(nzm, nzt, ngrdcol, gr, wp2, wp2_zt, w_tol_sqd); // Positive definite
  Skx_func(nzt, ngrdcol, wp2_zt, wp3, w_tol, clubb_params, Skw_zt);

  zt2zm_api(nzm, nzt, ngrdcol, gr, wp3, wp3_zm);
  Skx_func(nzm, ngrdcol, wp2, wp3_zm, w_tol, clubb_params, Skw_zm);

  zm2zt_api(nzm, nzt, ngrdcol, gr, thlp2, thlp2_zt, thl_tol * thl_tol);
  Skx_func(nzt, ngrdcol, thlp2_zt, thlp3, thl_tol, clubb_params, Skthl_zt);

  zm2zt_api(nzm, nzt, ngrdcol, gr, rtp2, rtp2_zt, rt_tol * rt_tol);
  Skx_func(nzt, ngrdcol, rtp2_zt, rtp3, rt_tol, clubb_params, Skrt_zt);

  zm2zt_api(nzm, nzt, ngrdcol, gr, up2, up2_zt, w_tol_sqd);
  Skx_func(nzt, ngrdcol, up2_zt, up3, w_tol, clubb_params, Sku_zt);

  zm2zt_api(nzm, nzt, ngrdcol, gr, vp2, vp2_zt, w_tol_sqd);
  Skx_func(nzt, ngrdcol, vp2_zt, vp3, w_tol, clubb_params, Skv_zt);

  if (l_call_pdf_closure_twice || (stats.l_sample && l_samp_stats_in_pdf_call)) {
    zt2zm_api(nzm, nzt, ngrdcol, gr, thlp3, thlp3_zm);
    Skx_func(nzm, ngrdcol, thlp2, thlp3_zm, thl_tol, clubb_params, Skthl_zm);

    zt2zm_api(nzm, nzt, ngrdcol, gr, rtp3, rtp3_zm);
    Skx_func(nzm, ngrdcol, rtp2, rtp3_zm, rt_tol, clubb_params, Skrt_zm);
  }

  zm2zt_api(nzm, nzt, ngrdcol, gr, sigma_sqd_w, sigma_sqd_w_zt, zero_threshold);
  zm2zt_api(nzm, nzt, ngrdcol, gr, wprtp, wprtp_zt);
  zm2zt_api(nzm, nzt, ngrdcol, gr, wpthlp, wpthlp_zt);
  zm2zt_api(nzm, nzt, ngrdcol, gr, rtpthlp, rtpthlp_zt);
  zm2zt_api(nzm, nzt, ngrdcol, gr, upwp, upwp_zt);
  zm2zt_api(nzm, nzt, ngrdcol, gr, vpwp, vpwp_zt);

  for (int sclr = 0; sclr < sclr_dim; ++sclr) {
    // Put passive scalar input on the t grid for the PDF.
    zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr),
              Kokkos::subview(wpsclrp_zt, Kokkos::ALL, Kokkos::ALL, sclr));
    zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(sclrprtp, Kokkos::ALL, Kokkos::ALL, sclr),
              Kokkos::subview(sclrprtp_zt, Kokkos::ALL, Kokkos::ALL, sclr));
    zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(sclrpthlp, Kokkos::ALL, Kokkos::ALL, sclr),
              Kokkos::subview(sclrpthlp_zt, Kokkos::ALL, Kokkos::ALL, sclr));

    // Positive definite quantity.
    zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, sclr),
              Kokkos::subview(sclrp2_zt, Kokkos::ALL, Kokkos::ALL, sclr),
              sclr_tol_host(sclr) * sclr_tol_host(sclr));
    Skx_func(nzt, ngrdcol, Kokkos::subview(sclrp2_zt, Kokkos::ALL, Kokkos::ALL, sclr),
             Kokkos::subview(sclrp3, Kokkos::ALL, Kokkos::ALL, sclr), sclr_tol_host(sclr),
             clubb_params, Kokkos::subview(Sksclr_zt, Kokkos::ALL, Kokkos::ALL, sclr));
  }

  // Interpolate hydrometeor mixed moments to momentum levels.
  for (int j = 0; j < hydromet_dim; ++j)
    zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wphydrometp, Kokkos::ALL, Kokkos::ALL, j),
              Kokkos::subview(wphydrometp_zt, Kokkos::ALL, Kokkos::ALL, j));

  //----------------------------------------------------------------
  // Call closure scheme.
  //----------------------------------------------------------------
  pdf_closure(workspace, nzt, ngrdcol, sclr_dim, sclr_tol, gr, // In
              hydromet_dim, p_in_Pa, exner, thv_ds_zt,         // In
              wm_zt, wp2_zt, wp3,                              // In
              Skw_zt, Skthl_zt, Skrt_zt, Sku_zt, Skv_zt,       // In
              rtm, rtp2_zt, wprtp_zt,                          // In
              thlm, thlp2_zt, wpthlp_zt,                       // In
              um, up2_zt, upwp_zt,                             // In
              vm, vp2_zt, vpwp_zt,                             // In
              rtpthlp_zt,                                      // In
              sclrm, wpsclrp_zt, sclrp2_zt,                    // In
              sclrprtp_zt, sclrpthlp_zt, Sksclr_zt,            // In
              wphydrometp_zt, wp2hmp,                          // In
              rtphmp_zt, thlphmp_zt,                           // In
              clubb_params, mixt_frac_max_mag,                 // In
              saturation_formula,                              // In
              stats,                                           // InOut
              iiPDF_type,                                      // In
              l_mix_rat_hm,                                    // In
              sigma_sqd_w_zt,                                  // InOut
              pdf_params, pdf_implicit_coefs_terms,            // InOut
              err_info,                                        // InOut
              wpup2, wpvp2,                                    // Out
              wp2up2_zt, wp2vp2_zt, wp4_zt,                    // Out
              wprtp2, wp2rtp,                                  // Out
              wpthlp2, wp2thlp, wprtpthlp,                     // Out
              cloud_frac, ice_supersat_frac,                   // Out
              rcm, wpthvp_zt, wp2thvp, wp2up, rtpthvp_zt,      // Out
              thlpthvp_zt, wprcp_zt, wp2rcp, rtprcp_zt,        // Out
              thlprcp_zt, rcp2_zt,                             // Out
              uprcp_zt, vprcp_zt,                              // Out
              w_up_in_cloud, w_down_in_cloud,                  // Out
              cloudy_updraft_frac, cloudy_downdraft_frac,      // Out
              wpsclrprtp, wpsclrp2, sclrpthvp_zt,              // Out
              wpsclrpthlp, sclrprcp_zt, wp2sclrp,              // Out
              rc_coef);                                        // Out

  // Subroutine may produce NaN values, and if so, return.
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "After first call to pdf_closure in pdf_closure_driver\n";
      return;
    }
  }

  if (l_rtm_nudge) {
    // Nudge rtm to prevent excessive drying.
    Kokkos::parallel_for(
        "pdf_closure_driver::rtm_nudge", ngrdcol * nzt, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          if (rtm(i, k) < rtm_min && zt(i, k) < rtm_nudge_max_altitude)
            rtm(i, k) = rtm(i, k) + (rtm_ref(i, k) - rtm(i, k)) * (dt / ts_nudge);
        });
  }

  if (l_call_pdf_closure_twice) {
    pdf_closure_driver_zm(workspace, gr, nzm, nzt, ngrdcol, hydromet_dim, // In
                          sclr_dim, sclr_tol, p_sfc,                      // In
                          mixt_frac_max_mag, clubb_params,                // In
                          iiPDF_type, saturation_formula,                 // In
                          l_mix_rat_hm, p_in_Pa, thv_ds_zm,               // In
                          wm_zm, wp2, wp3_zm,                             // In
                          Skw_zm, Skthl_zm, Skrt_zm,                      // In
                          rtm, rtp2, wprtp,                               // In
                          thlm, thlp2, wpthlp,                            // In
                          um, up2, upwp, up3,                             // In
                          vm, vp2, vpwp, vp3,                             // In
                          rtpthlp, sclrm, wpsclrp, sclrp2,                // In
                          sclrprtp, sclrpthlp, sclrp3,                    // In
                          wphydrometp,                                    // In
                          wp2hmp, rtphmp_zt, thlphmp_zt,                  // In
                          stats, sigma_sqd_w, pdf_params_zm,              // InOut
                          err_info,                                       // InOut
                          rtm_zm, thlm_zm,                                // Out
                          wp2up2, wp2vp2, wp4,                            // Out
                          wprtp2_zm, wp2rtp_zm,                           // Out
                          wpthlp2_zm, wp2thlp_zm, wprtpthlp_zm,           // Out
                          cloud_frac_zm, ice_supersat_frac_zm,            // Out
                          rcm_zm, wpthvp, wp2thvp_zm, wp2up_zm,           // Out
                          rtpthvp, thlpthvp, wprcp, wp2rcp_zm,            // Out
                          rtprcp, thlprcp, rcp2, uprcp, vprcp,            // Out
                          wpsclrprtp_zm, wpsclrp2_zm, sclrpthvp,          // Out
                          wpsclrpthlp_zm, sclrprcp, wp2sclrp_zm,          // Out
                          rc_coef_zm);                                    // Out

    // Subroutine may produce NaN values, and if so, return.
    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        std::cerr << err_info.err_header_global << '\n';
        std::cerr << "After second call to pdf_closure in pdf_closure_driver\n";
        return;
      }
    }
  } else {
    // Interpolate momentum variables output from the first call to
    // pdf_closure back to momentum grid.
    zt2zm_api(nzm, nzt, ngrdcol, gr, wp4_zt, wp4, zero_threshold);
    zt2zm_api(nzm, nzt, ngrdcol, gr, wpthvp_zt, wpthvp);
    zt2zm_api(nzm, nzt, ngrdcol, gr, thlpthvp_zt, thlpthvp);
    zt2zm_api(nzm, nzt, ngrdcol, gr, rtpthvp_zt, rtpthvp);
    zt2zm_api(nzm, nzt, ngrdcol, gr, wprcp_zt, wprcp);
    zt2zm_api(nzm, nzt, ngrdcol, gr, rc_coef, rc_coef_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, rtprcp_zt, rtprcp);
    zt2zm_api(nzm, nzt, ngrdcol, gr, thlprcp_zt, thlprcp);
    zt2zm_api(nzm, nzt, ngrdcol, gr, uprcp_zt, uprcp);
    zt2zm_api(nzm, nzt, ngrdcol, gr, vprcp_zt, vprcp);
    zt2zm_api(nzm, nzt, ngrdcol, gr, wp2up2_zt, wp2up2);
    zt2zm_api(nzm, nzt, ngrdcol, gr, wp2vp2_zt, wp2vp2);

    Kokkos::parallel_for(
        "pdf_closure_driver::momentum_boundaries", ngrdcol, KOKKOS_LAMBDA(const int i) {
          // Since top momentum level is higher than top thermo level,
          // set variables at top momentum level to 0.
          wp4(i, k_ub_zm) = zero;
          // Set wp4 to 0 at the lowest momentum level (momentum level 1).
          wp4(i, k_lb_zm)        = zero;
          wpthvp(i, k_ub_zm)     = zero;
          thlpthvp(i, k_ub_zm)   = zero;
          rtpthvp(i, k_ub_zm)    = zero;
          wprcp(i, k_ub_zm)      = zero;
          rc_coef_zm(i, k_ub_zm) = zero;
          rtprcp(i, k_ub_zm)     = zero;
          thlprcp(i, k_ub_zm)    = zero;
          uprcp(i, k_ub_zm)      = zero;
          vprcp(i, k_ub_zm)      = zero;
          wp2up2(i, k_ub_zm)     = zero;
          wp2vp2(i, k_ub_zm)     = zero;
        });

    // Initialize variables to avoid uninitialized variables.
    Kokkos::parallel_for(
        "pdf_closure_driver::momentum_initialization", ngrdcol * nzm,
        KOKKOS_LAMBDA(const int cell) {
          const int i                = cell % ngrdcol;
          const int k                = cell / ngrdcol;
          cloud_frac_zm(i, k)        = zero;
          ice_supersat_frac_zm(i, k) = zero;
          rcm_zm(i, k)               = zero;
          rtm_zm(i, k)               = zero;
          thlm_zm(i, k)              = zero;
        });

    // Interpolate passive scalars back onto the m grid.
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      zt2zm_api(nzm, nzt, ngrdcol, gr,
                Kokkos::subview(sclrpthvp_zt, Kokkos::ALL, Kokkos::ALL, sclr),
                Kokkos::subview(sclrpthvp, Kokkos::ALL, Kokkos::ALL, sclr));
      zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(sclrprcp_zt, Kokkos::ALL, Kokkos::ALL, sclr),
                Kokkos::subview(sclrprcp, Kokkos::ALL, Kokkos::ALL, sclr));
      Kokkos::parallel_for(
          "pdf_closure_driver::scalar_boundaries", ngrdcol, KOKKOS_LAMBDA(const int i) {
            sclrpthvp(i, k_ub_zm, sclr) = zero;
            sclrprcp(i, k_ub_zm, sclr)  = zero;
          });
    }

    if (stats.l_sample && l_samp_stats_in_pdf_call && var_on_stats_list(stats, "rcp2")) {
      // Positive definite quantity.
      zt2zm_api(nzm, nzt, ngrdcol, gr, rcp2_zt, rcp2, zero_threshold);
      Kokkos::parallel_for(
          "pdf_closure_driver::rcp2_boundary", ngrdcol,
          KOKKOS_LAMBDA(const int i) { rcp2(i, k_ub_zm) = zero; });
    }
  }

  // If l_trapezoidal_rule_zt is true, call trapezoidal_rule_zt for
  // thermodynamic-level variables output from pdf_closure.
  // ldgrant June 2009
  if (l_trapezoidal_rule_zt) {
    trapezoidal_rule_zt(nzm, nzt, ngrdcol, sclr_dim, gr,             // In
                        l_call_pdf_closure_twice,                    // In
                        stats,                                       // In
                        wprtp2, wpthlp2,                             // InOut
                        wprtpthlp, cloud_frac, ice_supersat_frac,    // InOut
                        rcm, wp2thvp, wp2up, wpsclrprtp, wpsclrp2,   // InOut
                        wpsclrpthlp,                                 // InOut
                        wprtp2_zm, wpthlp2_zm,                       // InOut
                        wprtpthlp_zm, cloud_frac_zm,                 // InOut
                        ice_supersat_frac_zm, rcm_zm, wp2thvp_zm,    // InOut
                        wp2up_zm,                                    // InOut
                        wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm); // InOut
  } else {
    zt2zm_api(nzm, nzt, ngrdcol, gr, cloud_frac, cloud_frac_zm);
    // Since top momentum level is higher than top thermo. level,
    // set variables at top momentum level to 0.
    Kokkos::parallel_for(
        "pdf_closure_driver::cloud_fraction_boundary", ngrdcol,
        KOKKOS_LAMBDA(const int i) { cloud_frac_zm(i, k_ub_zm) = zero; });
  }

  // If l_trapezoidal_rule_zm is true, call trapezoidal_rule_zm for
  // the important momentum-level variables output from pdf_closure.
  // ldgrant Feb. 2010
  if (l_trapezoidal_rule_zm) {
    trapezoidal_rule_zm(nzm, nzt, ngrdcol, gr,              // In
                        wpthvp_zt, thlpthvp_zt, rtpthvp_zt, // In
                        wpthvp, thlpthvp, rtpthvp);         // InOut
  }

  // Vince Larson clipped rcm in order to prevent rvm < 0. 5 Apr 2008.
  // This code won't work unless rtm >= 0 !!!
  // We do not clip rcm_in_layer because rcm_in_layer only influences
  // radiation, and we do not want to bother recomputing it.
  // Code is duplicated from below to ensure that relative humidity
  // is calculated properly. 3 Sep 2009.
  clip_rcm(nzt, ngrdcol, rtm, "rtm < rcm after pdf_closure", rcm);

  // Compute variables cloud_cover and rcm_in_layer.
  // Added July 2009.
  compute_cloud_cover(workspace, gr, nzt, ngrdcol, pdf_params, cloud_frac, rcm, err_info,
                      cloud_cover, rcm_in_layer);

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "calling compute_cloud_cover in pdf_closure_driver\n";
      return;
    }
  }

  if (l_use_cloud_cover) {
    // Use cloud_cover and rcm_in_layer to help boost cloud_frac and rcm to help
    // increase cloudiness at coarser grid resolutions.
    Kokkos::parallel_for(
        "pdf_closure_driver::use_cloud_cover", ngrdcol * nzt, KOKKOS_LAMBDA(const int cell) {
          const int i      = cell % ngrdcol;
          const int k      = cell / ngrdcol;
          cloud_frac(i, k) = cloud_cover(i, k);
          rcm(i, k)        = rcm_in_layer(i, k);
        });
    // rcm_in_layer divides rcm by the vertical cloud fraction at a cloud top
    // or base, so it can exceed rtm. Clip rcm again so that rvm >= 0.
    clip_rcm(nzt, ngrdcol, rtm, "rtm < rcm_in_layer after pdf_closure", rcm);
  }

  Kokkos::parallel_for(
      "pdf_closure_driver::clip_cloud_fraction", ngrdcol * nzt, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        // Native MIN selects its second argument on a tie.
        // Clip cloud fraction here if it still exceeds 1.0 due to roundoff.
        cloud_frac(i, k) = real_t(1) < cloud_frac(i, k) ? real_t(1) : cloud_frac(i, k);
        // Ditto with ice cloud fraction.
        ice_supersat_frac(i, k) =
            real_t(1) < ice_supersat_frac(i, k) ? real_t(1) : ice_supersat_frac(i, k);
      });

  if (l_rcm_supersat_adj) {
    thlm2T_in_K_api(nzt, ngrdcol, thlm, exner, rcm, T_in_K);
    sat_mixrat_liq_api(nzt, ngrdcol, gr, p_in_Pa, T_in_K, saturation_formula, rsat);
    Kokkos::parallel_for(
        "pdf_closure_driver::relative_humidity", ngrdcol * nzt, KOKKOS_LAMBDA(const int cell) {
          const int i            = cell % ngrdcol;
          const int k            = cell / ngrdcol;
          rel_humidity(i, k)     = (rtm(i, k) - rcm(i, k)) / rsat(i, k);
          rcm_supersat_adj(i, k) = zero;
        });

    // +PAB mods, take remaining supersaturation that may exist
    // after CLUBB PDF call and add it to rcm. Supersaturation may exist
    // after PDF call due to issues with calling PDF on the thermo grid
    // and momentum grid and the interpolation between the two.
    // A nonzero integer reduction represents the native logical OR.
    int l_spur_supersat = 0;
    Kokkos::parallel_reduce(
        "pdf_closure_driver::supersaturation_adjustment", ngrdcol * nzt,
        KOKKOS_LAMBDA(const int cell, int& count) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          if (rel_humidity(i, k) > real_t(1)) {
            rcm_supersat_adj(i, k) = (rtm(i, k) - rcm(i, k)) - rsat(i, k);
            rcm(i, k)              = rcm(i, k) + rcm_supersat_adj(i, k);
            count += 1;
          }
        },
        l_spur_supersat);
    if (clubb_at_least_debug_level_api(1) && l_spur_supersat)
      std::cerr << "Warning: spurious supersaturation was removed after pdf_closure!\n";
  } else if (stats.l_sample && l_samp_stats_in_pdf_call) {
    Kokkos::deep_copy(rcm_supersat_adj, zero);
  }

  if (stats.l_sample && l_samp_stats_in_pdf_call) {
    stats_update("Skw_zt", Skw_zt, stats);
    stats_update("Skw_zm", Skw_zm, stats);
    stats_update("Skthl_zt", Skthl_zt, stats);
    stats_update("Skthl_zm", Skthl_zm, stats);
    stats_update("Skrt_zt", Skrt_zt, stats);
    stats_update("Skrt_zm", Skrt_zm, stats);
    stats_update("uprcp", uprcp, stats);
    stats_update("vprcp", vprcp, stats);

    // Compute skewness velocity for stats output purposes.
    if (var_on_stats_list(stats, "Skw_velocity")) {
      Kokkos::parallel_for(
          "pdf_closure_driver::skewness_velocity", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int i          = cell % ngrdcol;
            const int k          = cell / ngrdcol;
            const real_t wp2_max = wp2(i, k) > w_tol_sqd ? wp2(i, k) : w_tol_sqd;
            // NVHPC -O0 folds the native reciprocal times quotient into one
            // division, rounding the combined denominator before that division.
            Skw_velocity(i, k) = wp3_zm(i, k) / (wp2_max * (real_t(1) - sigma_sqd_w(i, k)));
          });
    }

    stats_update("rc_coef", rc_coef, stats);
    stats_update("wp2rcp", wp2rcp, stats);
    stats_update("rtprcp", rtprcp, stats);
    stats_update("rcp2", rcp2, stats);
    stats_update("Skw_velocity", Skw_velocity, stats);
    stats_update("cloud_frac", cloud_frac, stats);
    stats_update("ice_supersat_frac", ice_supersat_frac, stats);
    stats_update("rcm_in_layer", rcm_in_layer, stats);
    stats_update("cloud_cover", cloud_cover, stats);
    stats_update("sigma_sqd_w", sigma_sqd_w, stats);
    stats_update("wpthvp", wpthvp, stats);
    stats_update("wp2thvp", wp2thvp, stats);
    stats_update("wp2up", wp2up, stats);
    stats_update("rtpthvp", rtpthvp, stats);
    stats_update("thlpthvp", thlpthvp, stats);
    stats_update("wprcp", wprcp, stats);
    stats_update("rc_coef_zm", rc_coef_zm, stats);
    stats_update("thlprcp", thlprcp, stats);
    stats_update("wpup2", wpup2, stats);
    stats_update("wpvp2", wpvp2, stats);
    stats_update("wp2up2", wp2up2, stats);
    stats_update("wp2vp2", wp2vp2, stats);
    stats_update("wp4", wp4, stats);
    stats_update("wp2rtp", wp2rtp, stats);
    stats_update("wprtp2", wprtp2, stats);
    stats_update("wp2thlp", wp2thlp, stats);
    stats_update("wpthlp2", wpthlp2, stats);
    stats_update("wprtpthlp", wprtpthlp, stats);
    stats_update("w_up_in_cloud", w_up_in_cloud, stats);
    stats_update("w_down_in_cloud", w_down_in_cloud, stats);
    stats_update("cld_updr_frac", cloudy_updraft_frac, stats);
    stats_update("cld_downdr_frac", cloudy_downdraft_frac, stats);
    stats_update("cloud_frac_zm", cloud_frac_zm, stats);
    stats_update("ice_supersat_frac_zm", ice_supersat_frac_zm, stats);
    stats_update("rtm_zm", rtm_zm, stats);
    stats_update("thlm_zm", thlm_zm, stats);
    stats_update("rcm_zm", rcm_zm, stats);
    stats_update("rcm_supersat_adj", rcm_supersat_adj, stats);
    stats_update("mixt_frac", pdf_params.mixt_frac, stats);
    stats_update("w_1", pdf_params.w_1, stats);
    stats_update("w_2", pdf_params.w_2, stats);
    stats_update("varnce_w_1", pdf_params.varnce_w_1, stats);
    stats_update("varnce_w_2", pdf_params.varnce_w_2, stats);
    stats_update("thl_1", pdf_params.thl_1, stats);
    stats_update("thl_2", pdf_params.thl_2, stats);
    stats_update("varnce_thl_1", pdf_params.varnce_thl_1, stats);
    stats_update("varnce_thl_2", pdf_params.varnce_thl_2, stats);
    stats_update("rt_1", pdf_params.rt_1, stats);
    stats_update("rt_2", pdf_params.rt_2, stats);
    stats_update("varnce_rt_1", pdf_params.varnce_rt_1, stats);
    stats_update("varnce_rt_2", pdf_params.varnce_rt_2, stats);
    stats_update("rc_1", pdf_params.rc_1, stats);
    stats_update("rc_2", pdf_params.rc_2, stats);
    stats_update("rsatl_1", pdf_params.rsatl_1, stats);
    stats_update("rsatl_2", pdf_params.rsatl_2, stats);
    stats_update("cloud_frac_1", pdf_params.cloud_frac_1, stats);
    stats_update("cloud_frac_2", pdf_params.cloud_frac_2, stats);
    stats_update("chi_1", pdf_params.chi_1, stats);
    stats_update("chi_2", pdf_params.chi_2, stats);
    stats_update("stdev_chi_1", pdf_params.stdev_chi_1, stats);
    stats_update("stdev_chi_2", pdf_params.stdev_chi_2, stats);
    stats_update("stdev_eta_1", pdf_params.stdev_eta_1, stats);
    stats_update("stdev_eta_2", pdf_params.stdev_eta_2, stats);
    stats_update("covar_chi_eta_1", pdf_params.covar_chi_eta_1, stats);
    stats_update("covar_chi_eta_2", pdf_params.covar_chi_eta_2, stats);
    stats_update("corr_w_chi_1", pdf_params.corr_w_chi_1, stats);
    stats_update("corr_w_chi_2", pdf_params.corr_w_chi_2, stats);
    stats_update("corr_w_eta_1", pdf_params.corr_w_eta_1, stats);
    stats_update("corr_w_eta_2", pdf_params.corr_w_eta_2, stats);
    stats_update("corr_chi_eta_1", pdf_params.corr_chi_eta_1, stats);
    stats_update("corr_chi_eta_2", pdf_params.corr_chi_eta_2, stats);
    stats_update("corr_w_rt_1", pdf_params.corr_w_rt_1, stats);
    stats_update("corr_w_rt_2", pdf_params.corr_w_rt_2, stats);
    stats_update("corr_w_thl_1", pdf_params.corr_w_thl_1, stats);
    stats_update("corr_w_thl_2", pdf_params.corr_w_thl_2, stats);
    stats_update("corr_rt_thl_1", pdf_params.corr_rt_thl_1, stats);
    stats_update("corr_rt_thl_2", pdf_params.corr_rt_thl_2, stats);
    stats_update("crt_1", pdf_params.crt_1, stats);
    stats_update("crt_2", pdf_params.crt_2, stats);
    stats_update("cthl_1", pdf_params.cthl_1, stats);
    stats_update("cthl_2", pdf_params.cthl_2, stats);

    if (var_on_stats_list(stats, "chi") || var_on_stats_list(stats, "chip2")) {
      const auto mixt_frac = pdf_params.mixt_frac;
      const auto chi_1     = pdf_params.chi_1;
      const auto chi_2     = pdf_params.chi_2;
      Kokkos::parallel_for(
          "pdf_closure_driver::mean_chi", ngrdcol * nzt, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol;
            const int k = cell / ngrdcol;
            // Match native contraction: round the first component product,
            // then fuse the second product with the sum.
            chi(i, k) = Kokkos::fma(real_t(1) - mixt_frac(i, k), chi_2(i, k),
                                    mixt_frac(i, k) * chi_1(i, k));
          });
      stats_update("chi", chi, stats);
    }

    if (var_on_stats_list(stats, "chip2")) {
      compute_variance_binormal(chi, pdf_params.chi_1, pdf_params.chi_2, pdf_params.stdev_chi_1,
                                pdf_params.stdev_chi_2, pdf_params.mixt_frac, chip2);
      stats_update("chip2", chip2, stats);
    }

    if (l_call_pdf_closure_twice) {
      stats_update("w_1_zm", pdf_params_zm.w_1, stats);
      stats_update("w_2_zm", pdf_params_zm.w_2, stats);
      stats_update("varnce_w_1_zm", pdf_params_zm.varnce_w_1, stats);
      stats_update("varnce_w_2_zm", pdf_params_zm.varnce_w_2, stats);
      stats_update("mixt_frac_zm", pdf_params_zm.mixt_frac, stats);
    }

    if (sclr_dim > 0) {
      for (int sclr = 0; sclr < sclr_dim; ++sclr) {
        const std::string sclr_idx = std::to_string(sclr + 1);
        stats_update("sclr" + sclr_idx + "pthvp",
                     Kokkos::subview(sclrpthvp, Kokkos::ALL, Kokkos::ALL, sclr), stats);
        stats_update("sclr" + sclr_idx + "prcp",
                     Kokkos::subview(sclrprcp, Kokkos::ALL, Kokkos::ALL, sclr), stats);
        stats_update("wp2sclr" + sclr_idx + "p",
                     Kokkos::subview(wp2sclrp, Kokkos::ALL, Kokkos::ALL, sclr), stats);
        stats_update("wpsclr" + sclr_idx + "p2",
                     Kokkos::subview(wpsclrp2, Kokkos::ALL, Kokkos::ALL, sclr), stats);
        stats_update("wpsclr" + sclr_idx + "prtp",
                     Kokkos::subview(wpsclrprtp, Kokkos::ALL, Kokkos::ALL, sclr), stats);
        stats_update("wpsclr" + sclr_idx + "pthlp",
                     Kokkos::subview(wpsclrpthlp, Kokkos::ALL, Kokkos::ALL, sclr), stats);
      }
    }
  }
}

namespace {

void pdf_closure_driver_zm(ScratchWorkspace& workspace,                                     // InOut
                           const grid& gr, int nzm, int nzt, int ngrdcol, int hydromet_dim, // In
                           int sclr_dim, Real1D sclr_tol, Real1D p_sfc,                     // In
                           Real1D mixt_frac_max_mag, Real2D clubb_params,                   // In
                           int iiPDF_type, int saturation_formula,                          // In
                           Bool1D l_mix_rat_hm, Real2D p_in_Pa, Real2D thv_ds_zm,           // In
                           Real2D wm_zm, Real2D wp2, Real2D wp3_zm,                         // In
                           Real2D Skw_zm, Real2D Skthl_zm, Real2D Skrt_zm,                  // In
                           Real2D rtm, Real2D rtp2, Real2D wprtp,                           // In
                           Real2D thlm, Real2D thlp2, Real2D wpthlp,                        // In
                           Real2D um, Real2D up2, Real2D upwp, Real2D up3,                  // In
                           Real2D vm, Real2D vp2, Real2D vpwp, Real2D vp3,                  // In
                           Real2D rtpthlp, Real3D sclrm, Real3D wpsclrp, Real3D sclrp2,     // In
                           Real3D sclrprtp, Real3D sclrpthlp, Real3D sclrp3,                // In
                           Real3D wphydrometp,                                              // In
                           Real3D wp2hmp, Real3D rtphmp_zt, Real3D thlphmp_zt,              // In
                           stats_type& stats, Real2D sigma_sqd_w,
                           pdf_parameter& pdf_params_zm,                              // InOut
                           err_info_type& err_info,                                   // InOut
                           Real2D rtm_zm, Real2D thlm_zm,                             // Out
                           Real2D wp2up2, Real2D wp2vp2, Real2D wp4,                  // Out
                           Real2D wprtp2_zm, Real2D wp2rtp_zm,                        // Out
                           Real2D wpthlp2_zm, Real2D wp2thlp_zm, Real2D wprtpthlp_zm, // Out
                           Real2D cloud_frac_zm, Real2D ice_supersat_frac_zm,         // Out
                           Real2D rcm_zm, Real2D wpthvp, Real2D wp2thvp_zm, Real2D wp2up_zm, // Out
                           Real2D rtpthvp, Real2D thlpthvp, Real2D wprcp, Real2D wp2rcp_zm,  // Out
                           Real2D rtprcp, Real2D thlprcp, Real2D rcp2, Real2D uprcp,
                           Real2D vprcp,                                               // Out
                           Real3D wpsclrprtp_zm, Real3D wpsclrp2_zm, Real3D sclrpthvp, // Out
                           Real3D wpsclrpthlp_zm, Real3D sclrprcp, Real3D wp2sclrp_zm, // Out
                           Real2D rc_coef_zm) {
  // Call the closure on momentum levels, using interpolated means and mixed
  // moments. Return its momentum-grid outputs while discarding thermodynamic-
  // only diagnostics and the local implicit coefficients.

  //------------------------------- Local Variables -------------------------------
  auto workspace_scope = workspace.scope();
  Real2D p_in_Pa_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::p_in_Pa_zm", ngrdcol,
                                   nzm); // Pressure interpolated to momentum levels [Pa]
  Real2D exner_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::exner_zm", ngrdcol,
                                   nzm); // Exner interpolated to momentum levels    [-]
  Real2D um_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::um_zm", ngrdcol,
                                   nzm); // Grid-mean eastward wind on m-levs.       [m/s]
  Real2D vm_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::vm_zm", ngrdcol,
                                   nzm); // Grid-mean northward wind on m-levs.      [m/s]
  Real2D up3_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::up3_zm", ngrdcol,
                                   nzm); // u'^3 interpolated to momentum levels     [m^3/s^3]
  Real2D vp3_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::vp3_zm", ngrdcol,
                                   nzm); // v'^3 interpolated to momentum levels     [m^3/s^3]
  Real2D Sku_zm = workspace_scope.view<Real2D>("pdf_closure_driver_zm::Sku_zm", ngrdcol,
                                               nzm); // Skewness of u on momentum levels         [-]
  Real2D Skv_zm = workspace_scope.view<Real2D>("pdf_closure_driver_zm::Skv_zm", ngrdcol,
                                               nzm); // Skewness of v on momentum levels         [-]
  Real2D wpup2_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::wpup2_zm", ngrdcol,
                                   nzm); // w'u'^2 on momentum grid                  [m^3/s^3]
  Real2D wpvp2_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::wpvp2_zm", ngrdcol,
                                   nzm); // w'v'^2 on momentum grid                  [m^3/s^3]
  Real2D w_up_in_cloud_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::w_up_in_cloud_zm", ngrdcol,
                                   nzm); // Avg. cloudy updraft velocity; m-levs     [m/s]
  Real2D w_down_in_cloud_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::w_down_in_cloud_zm", ngrdcol,
                                   nzm); // Avg. cloudy downdraft velocity; m-levs   [m/s]
  Real2D cloudy_updraft_frac_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::cloudy_updraft_frac_zm", ngrdcol,
                                   nzm); // Cloudy updraft fraction; m-levs          [-]
  Real2D cloudy_downdraft_frac_zm =
      workspace_scope.view<Real2D>("pdf_closure_driver_zm::cloudy_downdraft_frac_zm", ngrdcol,
                                   nzm); // Cloudy downdraft fraction; m-levs        [-]
  Real3D wp2hmp_zm = workspace_scope.view<Real3D>(
      "pdf_closure_driver_zm::wp2hmp_zm", ngrdcol, nzm,
      hydromet_dim); // Moment <w'^2 hm'> on m-levs. [(m/s)^2 <hm units>]
  Real3D rtphmp = workspace_scope.view<Real3D>(
      "pdf_closure_driver_zm::rtphmp", ngrdcol, nzm,
      hydromet_dim); // Covariance of rt and hm          [(kg/kg) <hm units>]
  Real3D thlphmp = workspace_scope.view<Real3D>(
      "pdf_closure_driver_zm::thlphmp", ngrdcol, nzm,
      hydromet_dim); // Covariance of thl and hm               [K <hm units>]
  Real3D sclrm_zm =
      workspace_scope.view<Real3D>("pdf_closure_driver_zm::sclrm_zm", ngrdcol, nzm,
                                   sclr_dim); // Passive scalar mean on momentum grid [units vary]
  Real3D sclrp3_zm = workspace_scope.view<Real3D>(
      "pdf_closure_driver_zm::sclrp3_zm", ngrdcol, nzm,
      sclr_dim); // sclr'^3 interpolated to momentum levels [{units vary}^3]
  Real3D Sksclr_zm =
      workspace_scope.view<Real3D>("pdf_closure_driver_zm::Sksclr_zm", ngrdcol, nzm,
                                   sclr_dim); // Skewness of sclr on momentum levels     [-]
  implicit_coefs_terms pdf_implicit_coefs_terms_zm;

  const int k_lb_zm        = gr.k_lb_zm;
  const int k_ub_zm        = gr.k_ub_zm;
  const int k_ub_zt        = gr.k_ub_zt;
  const auto sclr_tol_host = to_host_space(workspace, sclr_tol);

  //-------------------------------------- Begin Code --------------------------------------

  // Call pdf_closure a second time on momentum levels, to output (rather than
  // interpolate) the variables which belong on the momentum levels.

  // Interpolate sclrm to the momentum level for use in the second call to
  // pdf_closure.
  for (int sclr = 0; sclr < sclr_dim; ++sclr) {
    // Clip if extrapolation causes sclrm_zm to be less than sclr_tol.
    zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr),
              Kokkos::subview(sclrm_zm, Kokkos::ALL, Kokkos::ALL, sclr), sclr_tol_host(sclr));
  }

  // Interpolate pressure, p_in_Pa, to momentum levels. Since the surface
  // (or model lower boundary) is located at momentum level k = 1, the
  // pressure there is p_sfc.
  zt2zm_api(nzm, nzt, ngrdcol, gr, p_in_Pa, p_in_Pa_zm);
  Kokkos::parallel_for(
      "pdf_closure_driver_zm::pressure_boundaries", ngrdcol, KOKKOS_LAMBDA(const int i) {
        p_in_Pa_zm(i, k_lb_zm) = p_sfc(i);
        // Clip pressure if the extrapolation leads to a negative pressure.
        p_in_Pa_zm(i, k_ub_zm) = p_in_Pa_zm(i, k_ub_zm) > real_t(0.5) * p_in_Pa(i, k_ub_zt)
                                     ? p_in_Pa_zm(i, k_ub_zm)
                                     : real_t(0.5) * p_in_Pa(i, k_ub_zt);
      });

  // Set exner at momentum levels, exner_zm, based on p_in_Pa_zm.
  Kokkos::parallel_for(
      "pdf_closure_driver_zm::exner", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
        const int i    = cell % ngrdcol;
        const int k    = cell / ngrdcol;
        exner_zm(i, k) = Kokkos::pow(p_in_Pa_zm(i, k) / p0, kappa);
      });

  // Clip if extrapolation at the top level causes rtm_zm to be < rt_tol.
  zt2zm_api(nzm, nzt, ngrdcol, gr, rtm, rtm_zm, rt_tol);
  // Clip if extrapolation at the top level causes thlm_zm to be < thl_tol.
  zt2zm_api(nzm, nzt, ngrdcol, gr, thlm, thlm_zm, thl_tol);

  // Interpolate hydrometeor mixed moments to momentum levels.
  for (int j = 0; j < hydromet_dim; ++j) {
    zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(rtphmp_zt, Kokkos::ALL, Kokkos::ALL, j),
              Kokkos::subview(rtphmp, Kokkos::ALL, Kokkos::ALL, j));
    zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(thlphmp_zt, Kokkos::ALL, Kokkos::ALL, j),
              Kokkos::subview(thlphmp, Kokkos::ALL, Kokkos::ALL, j));
    zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wp2hmp, Kokkos::ALL, Kokkos::ALL, j),
              Kokkos::subview(wp2hmp_zm, Kokkos::ALL, Kokkos::ALL, j));
  }

  zt2zm_api(nzm, nzt, ngrdcol, gr, um, um_zm);
  zt2zm_api(nzm, nzt, ngrdcol, gr, vm, vm_zm);
  zt2zm_api(nzm, nzt, ngrdcol, gr, up3, up3_zm);
  zt2zm_api(nzm, nzt, ngrdcol, gr, vp3, vp3_zm);

  Skx_func(nzm, ngrdcol, up2, up3_zm, w_tol, clubb_params, Sku_zm);
  Skx_func(nzm, ngrdcol, vp2, vp3_zm, w_tol, clubb_params, Skv_zm);

  for (int sclr = 0; sclr < sclr_dim; ++sclr) {
    zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(sclrp3, Kokkos::ALL, Kokkos::ALL, sclr),
              Kokkos::subview(sclrp3_zm, Kokkos::ALL, Kokkos::ALL, sclr));
    Skx_func(nzm, ngrdcol, Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, sclr),
             Kokkos::subview(sclrp3_zm, Kokkos::ALL, Kokkos::ALL, sclr), sclr_tol_host(sclr),
             clubb_params, Kokkos::subview(Sksclr_zm, Kokkos::ALL, Kokkos::ALL, sclr));
  }

  // pdf_implicit_coefs_terms is only used in the iiPDF_new and
  // iiPDF_new_hybrid closures. Only initialize our local _zm version
  // for those choices.
  if (iiPDF_type == iiPDF_new || iiPDF_type == iiPDF_new_hybrid) {
    init_pdf_implicit_coefs_terms_api(nzm, ngrdcol, sclr_dim, pdf_implicit_coefs_terms_zm);
  }

  // Call pdf_closure to output the variables which belong on the momentum grid.
  pdf_closure(workspace, nzm, ngrdcol, sclr_dim, sclr_tol, gr,  // In
              hydromet_dim, p_in_Pa_zm, exner_zm, thv_ds_zm,    // In
              wm_zm, wp2, wp3_zm,                               // In
              Skw_zm, Skthl_zm, Skrt_zm, Sku_zm, Skv_zm,        // In
              rtm_zm, rtp2, wprtp,                              // In
              thlm_zm, thlp2, wpthlp,                           // In
              um_zm, up2, upwp,                                 // In
              vm_zm, vp2, vpwp,                                 // In
              rtpthlp,                                          // In
              sclrm_zm, wpsclrp, sclrp2,                        // In
              sclrprtp, sclrpthlp, Sksclr_zm,                   // In
              wphydrometp, wp2hmp_zm,                           // In
              rtphmp, thlphmp,                                  // In
              clubb_params, mixt_frac_max_mag,                  // In
              saturation_formula,                               // In
              stats,                                            // InOut
              iiPDF_type,                                       // In
              l_mix_rat_hm,                                     // In
              sigma_sqd_w,                                      // InOut
              pdf_params_zm, pdf_implicit_coefs_terms_zm,       // InOut
              err_info,                                         // InOut
              wpup2_zm, wpvp2_zm,                               // Out
              wp2up2, wp2vp2, wp4,                              // Out
              wprtp2_zm, wp2rtp_zm,                             // Out
              wpthlp2_zm, wp2thlp_zm, wprtpthlp_zm,             // Out
              cloud_frac_zm, ice_supersat_frac_zm,              // Out
              rcm_zm, wpthvp, wp2thvp_zm, wp2up_zm, rtpthvp,    // Out
              thlpthvp, wprcp, wp2rcp_zm, rtprcp,               // Out
              thlprcp, rcp2,                                    // Out
              uprcp, vprcp,                                     // Out
              w_up_in_cloud_zm, w_down_in_cloud_zm,             // Out
              cloudy_updraft_frac_zm, cloudy_downdraft_frac_zm, // Out
              wpsclrprtp_zm, wpsclrp2_zm, sclrpthvp,            // Out
              wpsclrpthlp_zm, sclrprcp, wp2sclrp_zm,            // Out
              rc_coef_zm);                                      // Out
}

void trapezoidal_rule_zt(int nzm, int nzt, int ngrdcol, int sclr_dim, const grid& gr,   // In
                         bool l_call_pdf_closure_twice,                                 // In
                         const stats_type& stats,                                       // In
                         Real2D wprtp2, Real2D wpthlp2,                                 // InOut
                         Real2D wprtpthlp, Real2D cloud_frac, Real2D ice_supersat_frac, // InOut
                         Real2D rcm, Real2D wp2thvp, Real2D wp2up, Real3D wpsclrprtp,
                         Real3D wpsclrp2,                                               // InOut
                         Real3D wpsclrpthlp,                                            // InOut
                         Real2D wprtp2_zm, Real2D wpthlp2_zm,                           // InOut
                         Real2D wprtpthlp_zm, Real2D cloud_frac_zm,                     // InOut
                         Real2D ice_supersat_frac_zm, Real2D rcm_zm, Real2D wp2thvp_zm, // InOut
                         Real2D wp2up_zm,                                               // InOut
                         Real3D wpsclrprtp_zm, Real3D wpsclrp2_zm, Real3D wpsclrpthlp_zm) {
  // Take thermo-grid closure outputs and either interpolate them to momentum
  // levels or reuse the second closure call's values, then vertically average
  // them back onto the thermo grid with calc_trapezoid_zt.
  //
  // ldgrant June 2009
  // The final momentum-grid arguments are InOut: the two-call branch keeps
  // their existing values, while the single-call branch fills them here.

  const int k_ub_zm = gr.k_ub_zm;

  // Store components of pdf_params in the locally declared variables.
  // We only apply the trapezoidal rule to these when
  // l_apply_rule_to_pdf_params is true. Applying the rule to the final
  // result of pdf_closure rather than intermediate results can make the PDF
  // component and cloud membership inconsistent, breaking preferential
  // in-cloud Latin-hypercube sampling. -dschanen 13 Feb 2012

  // If l_call_pdf_closure_twice is true, the _zm variables already have
  // values from the second call to pdf_closure in advance_clubb_core.
  // Otherwise interpolate them to the _zm levels.
  if (!l_call_pdf_closure_twice) {
    // Interpolate thermodynamic variables to the momentum grid.
    zt2zm_api(nzm, nzt, ngrdcol, gr, wprtp2, wprtp2_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, wpthlp2, wpthlp2_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, wprtpthlp, wprtpthlp_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, cloud_frac, cloud_frac_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, ice_supersat_frac, ice_supersat_frac_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, rcm, rcm_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, wp2thvp, wp2thvp_zm);
    zt2zm_api(nzm, nzt, ngrdcol, gr, wp2up, wp2up_zm);

    // Since top momentum level is higher than top thermo level,
    // set variables at top momentum level to 0.
    Kokkos::parallel_for(
        "trapezoidal_rule_zt::momentum_boundaries", ngrdcol, KOKKOS_LAMBDA(const int i) {
          wprtp2_zm(i, k_ub_zm)            = zero;
          wpthlp2_zm(i, k_ub_zm)           = zero;
          wprtpthlp_zm(i, k_ub_zm)         = zero;
          cloud_frac_zm(i, k_ub_zm)        = zero;
          ice_supersat_frac_zm(i, k_ub_zm) = zero;
          rcm_zm(i, k_ub_zm)               = zero;
          wp2thvp_zm(i, k_ub_zm)           = zero;
          wp2up_zm(i, k_ub_zm)             = zero;
        });

    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wpsclrprtp, Kokkos::ALL, Kokkos::ALL, sclr),
                Kokkos::subview(wpsclrprtp_zm, Kokkos::ALL, Kokkos::ALL, sclr));
      zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wpsclrp2, Kokkos::ALL, Kokkos::ALL, sclr),
                Kokkos::subview(wpsclrp2_zm, Kokkos::ALL, Kokkos::ALL, sclr));
      zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wpsclrpthlp, Kokkos::ALL, Kokkos::ALL, sclr),
                Kokkos::subview(wpsclrpthlp_zm, Kokkos::ALL, Kokkos::ALL, sclr));
      Kokkos::parallel_for(
          "trapezoidal_rule_zt::scalar_boundaries", ngrdcol, KOKKOS_LAMBDA(const int i) {
            wpsclrprtp_zm(i, k_ub_zm, sclr)  = zero;
            wpsclrp2_zm(i, k_ub_zm, sclr)    = zero;
            wpsclrpthlp_zm(i, k_ub_zm, sclr) = zero;
          });
    }
  }

  if (stats.enabled) {
    // Use the trapezoidal rule to recompute the variables on the stats_zt level.
    if (var_on_stats_list(stats, "wprtp2"))
      calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, wprtp2_zm, wprtp2);
    if (var_on_stats_list(stats, "wpthlp2"))
      calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, wpthlp2_zm, wpthlp2);
    if (var_on_stats_list(stats, "wprtpthlp"))
      calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, wprtpthlp_zm, wprtpthlp);

    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      calc_trapezoid_zt(nzm, nzt, ngrdcol, gr,
                        Kokkos::subview(wpsclrprtp_zm, Kokkos::ALL, Kokkos::ALL, sclr),
                        Kokkos::subview(wpsclrprtp, Kokkos::ALL, Kokkos::ALL, sclr));
      calc_trapezoid_zt(nzm, nzt, ngrdcol, gr,
                        Kokkos::subview(wpsclrpthlp_zm, Kokkos::ALL, Kokkos::ALL, sclr),
                        Kokkos::subview(wpsclrpthlp, Kokkos::ALL, Kokkos::ALL, sclr));
      calc_trapezoid_zt(nzm, nzt, ngrdcol, gr,
                        Kokkos::subview(wpsclrp2_zm, Kokkos::ALL, Kokkos::ALL, sclr),
                        Kokkos::subview(wpsclrp2, Kokkos::ALL, Kokkos::ALL, sclr));
    }
  }

  calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, cloud_frac_zm, cloud_frac);
  calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, ice_supersat_frac_zm, ice_supersat_frac);
  calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, rcm_zm, rcm);
  calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, wp2thvp_zm, wp2thvp);
  calc_trapezoid_zt(nzm, nzt, ngrdcol, gr, wp2up_zm, wp2up);

  // End of trapezoidal rule.
}

void trapezoidal_rule_zm(int nzm, int nzt, int ngrdcol, const grid& gr,           // In
                         Real2D wpthvp_zt, Real2D thlpthvp_zt, Real2D rtpthvp_zt, // In
                         Real2D wpthvp, Real2D thlpthvp, Real2D rtpthvp) {
  // Recompute the important momentum-grid buoyancy moments wpthvp,
  // thlpthvp and rtpthvp by vertically averaging the first closure call's
  // thermodynamic outputs with the momentum values already computed.
  // ldgrant Feb. 2010

  calc_trapezoid_zm(nzm, nzt, ngrdcol, gr, wpthvp_zt, wpthvp);
  calc_trapezoid_zm(nzm, nzt, ngrdcol, gr, thlpthvp_zt, thlpthvp);
  calc_trapezoid_zm(nzm, nzt, ngrdcol, gr, rtpthvp_zt, rtpthvp);
}

void calc_trapezoid_zt(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                       Real2D variable_zm,                            // In
                       Real2D variable_zt) {
  // Recompute thermodynamic-grid closure outputs using the trapezoidal rule
  // between the adjacent momentum levels and the original thermo-grid value.
  // Heights are in metres; the result retains variable_zt's field units.
  // ldgrant June 2009

  const auto zm           = gr.zm;
  const auto zt           = gr.zt;
  const auto invrs_dzt    = gr.invrs_dzt;
  const real_t grid_dir   = gr.grid_dir;
  const int grid_dir_indx = gr.grid_dir_indx;
  const int k_lb_zt       = gr.k_lb_zt;
  const int nlevels       = (gr.k_ub_zt - gr.k_lb_zt) / grid_dir_indx + 1;

  Kokkos::parallel_for(
      "calc_trapezoid_zt", ngrdcol * nlevels, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = k_lb_zt + (cell / ngrdcol) * grid_dir_indx;
        int k_zm, k_zmp1;
        if (grid_dir_indx > 0) {
          // Ascending grid.
          k_zmp1 = k + 1;
          k_zm   = k;
        } else {
          // Descending grid.
          k_zmp1 = k;
          k_zm   = k + 1;
        }
        // Trapezoidal rule from calculus. NVHPC -O0 rounds the upper term,
        // then fuses the lower term's final inverse-spacing multiplication.
        variable_zt(i, k) =
            Kokkos::fma(one_half * (variable_zt(i, k) + variable_zm(i, k_zm)) *
                            (zt(i, k) - zm(i, k_zm)) * grid_dir,
                        invrs_dzt(i, k),
                        one_half * (variable_zm(i, k_zmp1) + variable_zt(i, k)) *
                            (zm(i, k_zmp1) - zt(i, k)) * grid_dir * invrs_dzt(i, k));
      });
}

void calc_trapezoid_zm(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D variable_zt, // In
                       Real2D variable_zm) {
  // Recompute interior momentum-grid closure outputs with the trapezoidal
  // rule, retaining both boundary values. Heights are in metres and outputs
  // retain variable_zm's field units.
  // ldgrant Feb. 2010

  const auto zm           = gr.zm;
  const auto zt           = gr.zt;
  const auto invrs_dzm    = gr.invrs_dzm;
  const real_t grid_dir   = gr.grid_dir;
  const int grid_dir_indx = gr.grid_dir_indx;
  const int k_lb_zm       = gr.k_lb_zm;
  const int nlevels       = nzm > 1 ? (gr.k_ub_zm - gr.k_lb_zm) / grid_dir_indx - 1 : 0;

  // Boundary conditions: trapezoidal rule not valid at top zm level, nzmax.
  // Trapezoidal rule also is not used at zm level 1.
  Kokkos::parallel_for(
      "calc_trapezoid_zm", ngrdcol * nlevels, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = k_lb_zm + (cell / ngrdcol + 1) * grid_dir_indx;
        int k_zt, k_ztm1;
        if (grid_dir_indx > 0) {
          // Ascending grid.
          k_zt   = k;
          k_ztm1 = k - 1;
        } else {
          // Descending grid.
          k_zt   = k - 1;
          k_ztm1 = k;
        }
        // Trapezoidal rule from calculus; preserve the same native fusion
        // as calc_trapezoid_zt above.
        variable_zm(i, k) = Kokkos::fma(one_half * (variable_zm(i, k) + variable_zt(i, k_ztm1)) *
                                            (zm(i, k) - zt(i, k_ztm1)) * grid_dir,
                                        invrs_dzm(i, k),
                                        one_half * (variable_zt(i, k_zt) + variable_zm(i, k)) *
                                            (zt(i, k_zt) - zm(i, k)) * grid_dir * invrs_dzm(i, k));
      });
}

void compute_cloud_cover(ScratchWorkspace& workspace,                                    // InOut
                         const grid& gr, int nzt, int ngrdcol,                           // In
                         const pdf_parameter& pdf_params, Real2D cloud_frac, Real2D rcm, // In
                         err_info_type& err_info,                                        // InOut
                         Real2D cloud_cover, Real2D rcm_in_layer) {
  // Compute cloud cover (the amount of sky covered by cloud) and liquid
  // water in the cloudy part of each grid layer from cloud fraction, mean
  // liquid water and the PDF's extended cloud water. Set a per-column fatal
  // error if the native exhaustive cloud-neighbor cases cannot be selected.
  //
  // Reference: G. L. Mellor (1977), "The Gaussian Cloud Model Relations",
  // JAS, Vol. 34, pp. 356–358. Added July 2009.

  //------------------------ Local variables ------------------------
  auto workspace_scope = workspace.scope();
  Real2D chi_mean =
      workspace_scope.view<Real2D>("compute_cloud_cover::chi_mean", ngrdcol,
                                   nzt); // Mean extended cloud water of the two Gaussians [kg/kg]
  Real2D vert_cloud_frac_upper =
      workspace_scope.view<Real2D>("compute_cloud_cover::vert_cloud_frac_upper", ngrdcol,
                                   nzt); // Fraction of cloud in top half of grid box
  Real2D vert_cloud_frac_lower =
      workspace_scope.view<Real2D>("compute_cloud_cover::vert_cloud_frac_lower", ngrdcol,
                                   nzt); // Fraction of cloud in bottom half of grid box
  Real2D vert_cloud_frac =
      workspace_scope.view<Real2D>("compute_cloud_cover::vert_cloud_frac", ngrdcol,
                                   nzt); // Fraction of cloud filling the grid box in the vertical

  const auto mixt_frac    = pdf_params.mixt_frac;
  const auto chi_1        = pdf_params.chi_1;
  const auto chi_2        = pdf_params.chi_2;
  const auto zt           = gr.zt;
  const auto zm           = gr.zm;
  const auto invrs_dzm    = gr.invrs_dzm;
  const auto err_code     = err_info.err_code;
  const int grid_dir_indx = gr.grid_dir_indx;
  const real_t grid_dir   = gr.grid_dir;
  const int k_lb_zt       = gr.k_lb_zt;
  const int k_ub_zt       = gr.k_ub_zt;
  const int nlevels       = (gr.k_ub_zt - gr.k_lb_zt) / grid_dir_indx;

  //------------------------ Begin code ------------------------
  Kokkos::parallel_for(
      "compute_cloud_cover::mean_chi", ngrdcol * nzt, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        // Preserve the same native contraction as the sampled mean chi.
        chi_mean(i, k) =
            Kokkos::fma(real_t(1) - mixt_frac(i, k), chi_2(i, k), mixt_frac(i, k) * chi_1(i, k));
      });

  Kokkos::parallel_for(
      "compute_cloud_cover::cloud_layers", ngrdcol * nlevels, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = k_lb_zt + (cell / ngrdcol) * grid_dir_indx;
        int km1, kp1;
        if (grid_dir_indx > 0) {
          km1 = k > 0 ? k - 1 : 0;
          kp1 = k + 1 < nzt ? k + 1 : nzt - 1;
        } else {
          km1 = k + 1 < nzt ? k + 1 : nzt - 1;
          kp1 = k > 0 ? k - 1 : 0;
        }

        if (rcm(i, k) < rc_tol) { // No cloud at this level.
          cloud_cover(i, k)  = cloud_frac(i, k);
          rcm_in_layer(i, k) = rcm(i, k);
        } else if (rcm(i, kp1) >= rc_tol && rcm(i, km1) >= rc_tol) {
          // There is cloud above and below,
          // so assume cloud fills grid box from top to bottom.
          cloud_cover(i, k)  = cloud_frac(i, k);
          rcm_in_layer(i, k) = rcm(i, k);
        } else if (rcm(i, kp1) < rc_tol || rcm(i, km1) < rc_tol) {
          // Cloud may fail to reach gridbox top or base or both.

          // First let cloud fill the entire grid box, then overwrite
          // vert_cloud_frac_upper/lower for cloud top, base or one-point cloud.
          vert_cloud_frac_upper(i, k) = real_t(0.5);
          vert_cloud_frac_lower(i, k) = real_t(0.5);
          int k_zmp1, k_zm;
          if (grid_dir_indx > 0) {
            // Ascending grid.
            k_zmp1 = k + 1;
            k_zm   = k;
          } else {
            // Descending grid.
            k_zmp1 = k;
            k_zm   = k + 1;
          }

          if (rcm(i, kp1) < rc_tol) { // Cloud top.
            vert_cloud_frac_upper(i, k) =
                ((real_t(0.5) / (grid_dir * invrs_dzm(i, k_zmp1))) / (zm(i, k_zmp1) - zt(i, k))) *
                (rcm(i, k) / (rcm(i, k) + Kokkos::abs(chi_mean(i, kp1))));
            vert_cloud_frac_upper(i, k) = real_t(0.5) < vert_cloud_frac_upper(i, k)
                                              ? real_t(0.5)
                                              : vert_cloud_frac_upper(i, k);

            // Make the transition in cloudiness more gradual than the min alone.
            vert_cloud_frac_upper(i, k) =
                vert_cloud_frac_upper(i, k) +
                ((rcm(i, kp1) / rc_tol) * (real_t(0.5) - vert_cloud_frac_upper(i, k)));
          } else {
            vert_cloud_frac_upper(i, k) = real_t(0.5);
          }

          if (rcm(i, km1) < rc_tol) { // Cloud base.
            vert_cloud_frac_lower(i, k) =
                ((real_t(0.5) / (grid_dir * invrs_dzm(i, k_zm))) / (zt(i, k) - zm(i, k_zm))) *
                (rcm(i, k) / (rcm(i, k) + Kokkos::abs(chi_mean(i, km1))));
            vert_cloud_frac_lower(i, k) = real_t(0.5) < vert_cloud_frac_lower(i, k)
                                              ? real_t(0.5)
                                              : vert_cloud_frac_lower(i, k);

            // Make the transition in cloudiness more gradual than the min alone.
            vert_cloud_frac_lower(i, k) =
                vert_cloud_frac_lower(i, k) +
                ((rcm(i, km1) / rc_tol) * (real_t(0.5) - vert_cloud_frac_lower(i, k)));
          } else {
            vert_cloud_frac_lower(i, k) = real_t(0.5);
          }

          vert_cloud_frac(i, k) = vert_cloud_frac_upper(i, k) + vert_cloud_frac_lower(i, k);
          const real_t limited_vertical_fraction =
              real_t(1) < vert_cloud_frac(i, k) ? real_t(1) : vert_cloud_frac(i, k);
          vert_cloud_frac(i, k) = cloud_frac(i, k) > limited_vertical_fraction
                                      ? cloud_frac(i, k)
                                      : limited_vertical_fraction;
          cloud_cover(i, k)     = cloud_frac(i, k) / vert_cloud_frac(i, k);
          rcm_in_layer(i, k)    = rcm(i, k) / vert_cloud_frac(i, k);
        } else {
          // This case should be literally unreachable: all possible options
          // are covered above. Preserve the native sentinel and per-column error.
          cloud_cover(i, k)  = unused_var;
          rcm_in_layer(i, k) = unused_var;
          // Different failing levels may report the same column concurrently.
          Kokkos::atomic_exchange(&err_code(i), clubb_fatal_error);
        }
      });

  Kokkos::parallel_for(
      "compute_cloud_cover::top_boundary", ngrdcol, KOKKOS_LAMBDA(const int i) {
        cloud_cover(i, k_ub_zt)  = cloud_frac(i, k_ub_zt);
        rcm_in_layer(i, k_ub_zt) = rcm(i, k_ub_zt);
      });

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      const auto mixt_frac_host  = to_host_space(workspace, pdf_params.mixt_frac);
      const auto chi_1_host      = to_host_space(workspace, pdf_params.chi_1);
      const auto chi_2_host      = to_host_space(workspace, pdf_params.chi_2);
      const auto cloud_frac_host = to_host_space(workspace, cloud_frac);
      const auto rcm_host        = to_host_space(workspace, rcm);

      // Retain distinguishing digits in the failure dump without changing
      // the formatting of later messages on the shared stream.
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "ERROR: compute_cloud_cover entered a conditional case it should not have \n";
      diagnostic << "pdf_params%mixt_frac = ";
      for (int k = 0; k < nzt; ++k)
        for (int i = 0; i < ngrdcol; ++i)
          diagnostic << mixt_frac_host(i, k) << ' ';
      diagnostic << '\n' << "pdf_params%chi_1 = ";
      for (int k = 0; k < nzt; ++k)
        for (int i = 0; i < ngrdcol; ++i)
          diagnostic << chi_1_host(i, k) << ' ';
      diagnostic << '\n' << "pdf_params%chi_2 = ";
      for (int k = 0; k < nzt; ++k)
        for (int i = 0; i < ngrdcol; ++i)
          diagnostic << chi_2_host(i, k) << ' ';
      diagnostic << '\n' << "cloud_frac = ";
      for (int k = 0; k < nzt; ++k)
        for (int i = 0; i < ngrdcol; ++i)
          diagnostic << cloud_frac_host(i, k) << ' ';
      diagnostic << '\n' << "rcm = ";
      for (int k = 0; k < nzt; ++k)
        for (int i = 0; i < ngrdcol; ++i)
          diagnostic << rcm_host(i, k) << ' ';
      diagnostic << '\n';
      std::cerr << diagnostic.str();
    }
  }
}

} // namespace

} // namespace clubb_kokkos
