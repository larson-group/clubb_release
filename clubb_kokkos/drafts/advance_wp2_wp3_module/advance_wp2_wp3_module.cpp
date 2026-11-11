// Advance the coupled vertical-velocity variance and third moment by one timestep.
//
// advance_wp2_wp3 assembles, solves and clips the Native semi-implicit wp2/wp3
// equations, including their sampled budgets. The public declaration documents
// physical field layout; INTEGRATION.md records lower-call dependencies and the
// existing sponge exclusion. Draft compilation and numerical checks are pending.
#include "advance_wp2_wp3_module.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "model_flags.hpp"
#include "error_code.hpp"
#include "Skx_module.hpp"
#include "advance_helper_module.hpp"
#include "matrix_solver_wrapper.hpp"
#include "f2k/procedure_bridge.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <utility>

namespace clubb_kokkos {
namespace {
// Private named constants avoid string comparisons. clip_wp2 must match the
// named constant in Native clip_explicit.
constexpr int clip_wp2 = 12;
// Set true for Crank-Nicholson diffusion, false for completely implicit diffusion.
// The latter stabilized deep convection with almost no boundary-layer effect
// in the Native comparison (Brian, 2008-01-04).
constexpr bool l_crank_nich_diff = false;
constexpr int ndiags2 = 2, ndiags3 = 3, ndiags5 = 5;

void wp23_solve(int nzm, int nzt, int ngrdcol, const grid& gr, real_t dt,              // In
                Real3D lhs, Real2D rhs,                                                // InOut
                Real3D lhs_ma_zm, Real2D lhs_dp1_wp2, Real3D lhs_diff_zm,              // In
                Real3D lhs_ta_wp2, Real2D lhs_pr1_wp2, Real2D lhs_pr1_wp3,             // In
                Real3D lhs_diff_zt, Real3D lhs_adv_tp_wp3, Real3D lhs_pr_tp_wp3,       // In
                Real3D wp3_pr3_lhs, Real3D lhs_ma_zt,                                  // In
                Real3D wp3_term_ta_lhs_result,                                         // In
                Real2D wm_zm, Real2D wm_zt,                                            // In
                Real1D sfc_elevation, Real2D C11_Skw_fnc,                              // In
                Real2D rho_ds_zm,                                                      // In
                Real2D upwp, Real2D vpwp, Real2D wprtp, Real2D wpthlp, Real2D rtp2,    // In
                Real2D thlp2,                                                          // In
                Real2D clubb_params,                                                   // In
                int penta_solve_method,                                                // In
                int fill_holes_type,                                                   // In
                bool l_min_wp2_from_corr_wx,                                           // In
                bool l_use_tke_in_wp2_wp3_K_dfsn,                                      // In
                bool l_use_wp3_lim_with_smth_Heaviside,                                // In
                bool l_wp2_fill_holes_tke,                                             // In
                bool l_implemented,                                                    // In
                f2k::stats_type& stats,                                                // InOut
                Real2D up2, Real2D vp2, Real2D wp2, Real2D wp3, Real2D wp3_zm,         // InOut
                Real2D wp2_zt, err_info_type& err_info);                               // InOut
void wp23_lhs(int nzm, int nzt, int ngrdcol, real_t dt,                                // In
              Real3D wp3_term_ta_lhs_result,                                           // In
              Real3D lhs_diff_zm, Real3D lhs_diff_zt, Real3D lhs_ma_zm,                // In
              Real3D lhs_ma_zt, Real3D lhs_ta_wp2,                                     // In
              Real3D lhs_tp_wp3,                                                       // In
              Real2D lhs_ac_pr2_wp2, Real2D lhs_ac_pr2_wp3, Real2D lhs_dp1_wp2,        // In
              Real2D lhs_pr1_wp3, Real2D lhs_pr1_wp2, Real2D lhs_splat_wp2,            // In
              Real2D lhs_splat_wp3,                                                    // In
              Real3D lhs);                                                             // Out
void wp23_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, real_t dt, Real1D fcor_y, // In
              Real3D wp3_term_ta_lhs_result,                                           // In
              Real3D lhs_diff_zm, Real3D lhs_diff_zt, Real3D lhs_diff_zm_crank,        // In
              Real3D lhs_diff_zt_crank,                                                // In
              Real3D lhs_tp_wp3, Real3D lhs_adv_tp_wp3, Real3D lhs_pr_tp_wp3,          // In
              Real3D lhs_ta_wp3, Real2D lhs_dp1_wp2, Real2D rhs_dp1_wp2, Real2D lhs_pr1_wp2, // In
              Real2D rhs_pr1_wp2, Real2D lhs_pr1_wp3, Real2D rhs_pr1_wp3,                    // In
              Real2D rhs_bp_pr2_wp2,                                                         // In
              Real2D rhs_pr_dfsn_wp2, Real2D rhs_bp1_pr2_wp3, Real2D rhs_pr3_wp2,            // In
              Real2D rhs_pr3_wp3,                                                            // In
              Real2D rhs_ta_wp3, Real2D rhs_pr_turb_wp3, Real2D rhs_pr_dfsn_wp3,             // In
              Real2D wp2, Real2D wp3, Real2D wpup2, Real2D wpvp2,                            // In
              Real2D wpthvp, Real2D wp2thvp, Real2D wp2up, Real2D up2, Real2D vp2,           // In
              Real2D upwp,                                                                   // In
              Real2D C11_Skw_fnc, Real2D thv_ds_zm, Real2D thv_ds_zt,                        // In
              Real2D lhs_splat_wp2, Real2D lhs_splat_wp3,                                    // In
              Real2D clubb_params,                                                           // In
              int iiPDF_type,                                                                // In
              bool l_use_tke_in_wp2_wp3_K_dfsn,                                              // In
              bool l_ho_nontrad_coriolis,                                                    // In
              f2k::stats_type& stats,                                    // InOut
              Real2D rhs);                                               // Out
void wp2_term_ta_lhs(int nzm, int nzt, int ngrdcol, const grid& gr,      // In
                     Real2D rho_ds_zt, Real2D invrs_rho_ds_zm,           // In
                     Real3D lhs_ta_wp2);                                 // Out
void wp2_terms_ac_pr2_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                          Real1D C_uu_shr, Real2D wm_zt,                 // In
                          Real2D lhs_ac_pr2_wp2);                        // Out
void wp2_term_dp1_lhs(int nzm, int ngrdcol, const grid& gr,              // In
                      Real2D C1_Skw_fnc, Real2D invrs_tau1m,             // In
                      Real2D lhs_dp1_wp2);                               // Out
void wp2_term_pr1_lhs(int nzm, int ngrdcol, const grid& gr, Real1D C4, Real2D invrs_tau_C4_zm, // In
                      Real2D lhs_pr1_wp2);                                          // Out
void wp2_terms_bp_pr2_rhs(int nzm, int ngrdcol, const grid& gr,                     // In
                          Real1D C_uu_buoy,                                         // In
                          Real2D thv_ds_zm, Real2D wpthvp,                          // In
                          Real2D rhs_bp_pr2_wp2);                                   // Out
void wp2_term_dp1_rhs(int nzm, int ngrdcol, const grid& gr, Real2D C1_Skw_fnc,      // In
                      Real2D invrs_tau1m, real_t threshold, Real2D up2, Real2D vp2, // In
                      bool l_damp_wp2_using_em,                                     // In
                      Real2D rhs_dp1_wp2);                                          // Out
void wp2_term_pr3_rhs(int nzm, int nzt, int ngrdcol, const grid& gr,                // In
                      Real1D C_uu_shr,                                              // In
                      Real1D C_uu_buoy,                                             // In
                      Real2D thv_ds_zm, Real2D wpthvp, Real2D upwp,                 // In
                      Real2D um, Real2D vpwp, Real2D vm,                            // In
                      Real2D rhs_pr3_wp2);                                          // Out
void wp2_term_pr1_rhs(int nzm, int ngrdcol, const grid& gr, Real1D C4,              // In
                      Real2D up2, Real2D vp2, Real2D invrs_tau_C4_zm,               // In
                      Real2D rhs_pr1_wp2);                                          // Out
void wp2_term_pr_dfsn_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real1D C_wp2_pr_dfsn, // In
                          Real2D rho_ds_zt, Real2D invrs_rho_ds_zm,                            // In
                          Real2D wpup2, Real2D wpvp2, Real2D wp3,                              // In
                          Real2D rhs_pr_dfsn_wp2);                                           // Out
void wp3_term_ta_new_pdf_lhs(int nzm, int nzt, int ngrdcol, const grid& gr,                  // In
                             Real2D coef_wp4_implicit,                                       // In
                             Real2D wp2, Real2D rho_ds_zm, Real2D invrs_rho_ds_zt,           // In
                             Real3D lhs_ta_wp3);                                             // Out
void wp3_term_ta_ADG1_lhs(int nzm, int nzt, int ngrdcol, const grid& gr,                     // In
                          Real2D wp2, Real2D a1_coef, Real2D a1_coef_zt,                     // In
                          Real2D a3_coef, Real2D a3_coef_zt,                                 // In
                          Real2D wp3_on_wp2, Real2D rho_ds_zm,                               // In
                          Real2D rho_ds_zt, Real2D invrs_rho_ds_zt,                          // In
                          bool l_standard_term_ta,                                           // In
                          bool l_partial_upwind_wp3,                                         // In
                          Real3D lhs_ta_wp3);                                                // Out
void wp3_term_tp_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real1D coef_wp3_tp,      // In
                     Real2D wp2, Real2D rho_ds_zm, Real2D invrs_rho_ds_zt,                   // In
                     Real3D lhs_tp_wp3);                                                     // Out
void wp3_terms_ac_pr2_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D C11_Skw_fnc, // In
                          Real2D wm_zm,                                                      // In
                          Real2D lhs_ac_pr2_wp3);                                            // Out
void wp3_term_pr1_lhs(int nzt, int ngrdcol, const grid& gr,                                  // In
                      Real1D C8, Real1D C8b,                                                 // In
                      Real2D invrs_tau_wp3_zt, Real2D Skw_zt,                                // In
                      bool l_damp_wp3_Skw_squared,                                           // In
                      Real2D lhs_pr1_wp3);                                                   // Out
void wp3_term_ta_explicit_rhs(int nzm, int nzt, int ngrdcol, const grid& gr,                 // In
                              Real2D wp4, Real2D rho_ds_zm, Real2D invrs_rho_ds_zt,          // In
                              Real2D rhs_ta_wp3);                                            // Out
void wp3_terms_bp1_pr2_rhs(int nzt, int ngrdcol, const grid& gr, Real2D C11_Skw_fnc,         // In
                           Real2D thv_ds_zt, Real2D wp2thvp,                                 // In
                           Real2D rhs_bp1_pr2_wp3);                                          // Out
void wp3_term_pr_turb_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real1D C_wp3_pr_turb, // In
                          Real2D Kh_zt, Real2D wpthvp, Real2D dum_dz, Real2D dvm_dz,           // In
                          Real2D upwp, Real2D vpwp,                                            // In
                          Real2D thv_ds_zt,                                                    // In
                          Real2D rho_ds_zm, Real2D invrs_rho_ds_zt,                            // In
                          Real2D em, Real2D wp2,                                               // In
                          Real2D rhs_pr_turb_wp3,              // Out
                          bool l_use_tke_in_wp3_pr_turb_term); // In
void wp3_term_pr_dfsn_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real1D C_wp3_pr_dfsn, // In
                          Real2D rho_ds_zm, Real2D invrs_rho_ds_zt,                            // In
                          Real2D wp2up2, Real2D wp2vp2, Real2D wp4,                            // In
                          Real2D up2, Real2D vp2, Real2D wp2,                                  // In
                          Real2D rhs_pr_dfsn_wp3);                        // Out
void wp3_term_pr1_rhs(int nzt, int ngrdcol, const grid& gr,               // In
                      Real1D C8, Real1D C8b,                              // In
                      Real2D invrs_tau_wp3_zt, Real2D Skw_zt, Real2D wp3, // In
                      bool l_damp_wp3_Skw_squared,                        // In
                      Real2D rhs_pr1_wp3);                                // Out
} // namespace

// =============================================================================
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
                     err_info_type& err_info) {                                             // InOut
  //
  // Description:
  // Advance w'^2 and w'^3 one timestep.
  //
  // References:
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wp2_wp3_eqns
  //
  // Eqn. 12 & 18 on p. 3545--3546 of
  // ``A PDF-Based Model for Boundary Layer Clouds. Part I:
  //   Method and Model Description'' Golaz, et al. (2002)
  // JAS, Vol. 59, pp. 3540--3551.
  //
  // See also
  // ``Equations for CLUBB'', Section 6:
  // /Implict solution for the vertical velocity moments/
  //------------------------------------------------------------------------

  // --------------------------- Input Variables ---------------------------
  // dt: Model timestep                            [s]
  // sfc_elevation: Elevation of ground level                 [m AMSL]
  // fcor_y: Nontraditional Coriolis parameter         [s^-1]
  // Meridional planetary vorticity. Proportional to cos(latitude)
  // sigma_sqd_w: sigma_sqd_w (momentum levels)             [-]
  // wm_zm: w wind component on momentum levels       [m/s]
  // wp2up2: w'^2u'^2 (momentum levels)                [m^4/s^4]
  // wp2vp2: w'^2v'^2 (momentum levels)                [m^4/s^4]
  // wp4: w'^4 (momentum levels)                    [m^4/s^4]
  // wpthvp: w'th_v' (momentum levels)                 [K m/s]
  // upwp: u'w' (momentum levels)                    [m^2/s^2]
  // vpwp: v'w' (momentum levels)                    [m^2/s^2]
  // em: Turbulence kinetic energy                 [m^2/s^2]
  // Kh_zm: Eddy diffusivity on momentum levels       [m^2/s]
  // invrs_tau_C4_zm: Inverse time-scale tau on momentum levels         [1/s]
  // invrs_tau_C1_zm: Inverse tau values used for the C1 (dp1) term in wp2 [1/s]
  // rho_ds_zm: Dry, static density on momentum levels    [kg/m^3]
  // invrs_rho_ds_zm: Inv. dry, static density @ momentum levs. [m^3/kg]
  // thv_ds_zm: Dry, base-state theta_v on momentum levs. [K]
  // wprtp: Flux of total water mixing ratio          [m/s kg/kg]
  // wpthlp: Flux of liquid water potential temp.      [m/s K]
  // rtp2: Variance of rt (overall)                  [kg^2/kg^2]
  // thlp2: Variance of thl (overall)                 [K^2]
  // Cx_fnc_Richardson: Cx_fnc from Richardson_num                [-]
  // lhs_splat_wp2: LHS coefficient of wp2 splatting term     [1/s]
  // wm_zt: w wind component on thermodynamic levels  [m/s]
  // wpup2: w'u'^2 (thermodynamic levels)             [m^3/s^3]
  // wpvp2: w'v'^2 (thermodynamic levels)             [m^3/s^3]
  // wp2thvp: w'^2th_v' (thermodynamic levels)          [K m^2/s^2]
  // wp2up: w'^2u' (thermodynamic levels)             [m^3/s^3]
  // um: u wind component (thermodynamic levels)   [m/s]
  // vm: v wind component (thermodynamic levels)   [m/s]
  // Kh_zt: Eddy diffusivity on thermodynamic levels  [m^2/s]
  // invrs_tau_wp3_zt: Inverse time-scale tau on thermodynamic levels    [1/s]
  // rho_ds_zt: Dry, static density on thermo. levels     [kg/m^3]
  // invrs_rho_ds_zt: Inv. dry, static density @ thermo. levs.  [m^3/kg]
  // thv_ds_zt: Dry, base-state theta_v on thermo. levs.  [K]
  // mixt_frac: Weight of 1st normal distribution         [-]
  // lhs_splat_wp3: LHS coefficient of wp3 splatting term     [1/s]
  // pdf_implicit_coefs_terms: Implicit coefs / explicit terms [units vary]
  // clubb_params: Array of CLUBB's tunable parameters    [units vary]
  // nu_vert_res_dep: Vertical resolution dependent nu values
  // iiPDF_type: Selected option for the two-component normal (double
  // Gaussian) PDF type to use for the w, rt, and theta-l (or
  // w, chi, and eta) portion of CLUBB's multivariate,
  // two-component PDF.
  // penta_solve_method: Method to solve then penta-diagonal system
  // fill_holes_type: Specifier for which hole filling method to use
  // l_min_wp2_from_corr_wx: Flag to base the threshold minimum value of wp2 on keeping the
  // overall correlation of w and x (w and rt, as well as w and
  // theta-l) within the limits of -max_mag_correlation_flux to
  // max_mag_correlation_flux.
  // l_upwind_xm_ma: This flag determines whether we want to use an upwind
  // differencing approximation rather than a centered differencing
  // for turbulent or mean advection terms. It affects rtm, thlm,
  // sclrm, um and vm.
  // l_standard_term_ta: Use the standard discretization for the turbulent advection
  // terms. Setting to .false. means that a_1 and a_3 are pulled
  // outside of the derivative in advance_wp2_wp3_module.F90 and in
  // advance_xp2_xpyp_module.F90.
  // l_partial_upwind_wp3: Flag to use an "upwind" discretization rather
  // than a centered discretization for the portion
  // of the wp3 turbulent advection term for ADG1
  // that is linearized in terms of wp3<t+1>.
  // (Requires ADG1 PDF and l_standard_term_ta).
  // l_damp_wp2_using_em: In the wp2 equation, use a dissipation formula of
  // -(2/3)*em/tau_zm,
  // as in Bougeault (1981)
  // l_use_C11_Richardson: Parameterize C16 based on Richardson number
  // l_damp_wp3_Skw_squared: Set damping on wp3 to use Skw^2 rather than Skw^4
  // l_lmm_stepping: Apply Linear Multistep Method (LMM) Stepping
  // l_use_tke_in_wp3_pr_turb_term: Use TKE formulation for wp3 pr_turb term
  // l_use_tke_in_wp2_wp3_K_dfsn: Use TKE in eddy diffusion for wp2 and wp3
  // l_use_wp3_lim_with_smth_Heaviside: Flag to activate mods on wp3 limiters for conv test
  // l_wp2_fill_holes_tke: Turn on additional hole-filling for wp2
  // that takes TKE from up2 and vp2, if necessary
  // l_ho_nontrad_coriolis: Flag to implement the nontraditional Coriolis terms in the
  // prognostic equations of <w'w'>, <u'w'>, and <u'u'>.
  // l_implemented: True if CLUBB is being implemented and run in a host model
  // --------------------- intent(inout) Variable ---------------------
  // --------------------------- Input/Output ---------------------------
  // up2: u'^2 (momentum levels)                    [m^2/s^2]
  // vp2: v'^2 (momentum levels)                    [m^2/s^2]
  // wp2: w'^2 (momentum levels)                    [m^2/s^2]
  // wp3: w'^3 (thermodynamic levels) [m^3/s^3]
  // err_info: err_info struct containing err_code and err_header

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D wp2_old("advance_wp2_wp3::wp2_old", ngrdcol, nzm); // w'^2 (momentum levels) [m^2/s^2]
  Real2D wp3_zm("advance_wp2_wp3::wp3_zm", ngrdcol,
                nzm); // w'^3 interpolated to momentum levels   [m^3/s^3]
  Real2D Skw_zm("advance_wp2_wp3::Skw_zm", ngrdcol, nzm); // Skewness of w on momentum levels [-]
  Real2D wp3_on_wp2("advance_wp2_wp3::wp3_on_wp2", ngrdcol,
                    nzm); // Smoothed version of wp3 / wp2         [m/s]
  Real2D wp2_smth("advance_wp2_wp3::wp2_smth", ngrdcol, nzm);
  Real2D em_smth("advance_wp2_wp3::em_smth", ngrdcol, nzm);
  Real2D wp3_old("advance_wp2_wp3::wp3_old", ngrdcol, nzt); // w'^3 (thermodynamic levels) [m^3/s^3]
  Real2D wp2_zt("advance_wp2_wp3::wp2_zt", ngrdcol,
                nzt); // w'^2 interpolated to thermodyamic levels [m^2/s^2]
  Real2D Skw_zt("advance_wp2_wp3::Skw_zt", ngrdcol,
                nzt); // Skewness of w on thermodynamic levels [-]
  Real2D wp3_on_wp2_zt("advance_wp2_wp3::wp3_on_wp2_zt", ngrdcol,
                       nzt); // Smoothed version of wp3 / wp2 on thermodynamic levels [m/s]
  Real2D Kw1("advance_wp2_wp3::Kw1", ngrdcol, nzt);
  Real2D Kw8("advance_wp2_wp3::Kw8", ngrdcol, nzm);
  Real2D C1_Skw_fnc("advance_wp2_wp3::C1_Skw_fnc", ngrdcol,
                    nzm); // C_1 parameter with Sk_w applied              [-]
  Real2D C11_Skw_fnc("advance_wp2_wp3::C11_Skw_fnc", ngrdcol,
                     nzt); // C_11 parameter with Sk_w applied             [-]
  Real2D C16_fnc("advance_wp2_wp3::C16_fnc", ngrdcol, nzt); // C_16 parameter [-]
  Real3D wp3_term_ta_lhs_result("advance_wp2_wp3::wp3_term_ta_lhs_result", ndiags5, ngrdcol, nzt);
  Real3D wp3_pr3_lhs("advance_wp2_wp3::wp3_pr3_lhs", ndiags5, ngrdcol, nzt);
  Real3D lhs_ta_wp2("advance_wp2_wp3::lhs_ta_wp2", ndiags2, ngrdcol,
                    nzm); // Turbulent advection terms for wp2
  Real3D lhs_tp_wp3("advance_wp2_wp3::lhs_tp_wp3", ndiags2, ngrdcol,
                    nzt); // Turbulent production terms of w'^3
  Real3D lhs_adv_tp_wp3("advance_wp2_wp3::lhs_adv_tp_wp3", ndiags2, ngrdcol,
                        nzt); // Turbulent production terms of w'^3 (for stats)
  Real3D lhs_pr_tp_wp3(
      "advance_wp2_wp3::lhs_pr_tp_wp3", ndiags2, ngrdcol,
      nzt); // Pressure scrambling terms for turbulent production of w'^3 (for stats)
  Real3D lhs_ta_wp3("advance_wp2_wp3::lhs_ta_wp3", ndiags2, ngrdcol,
                    nzt); // Turbulent advection terms for wp3
  Real2D lhs_dp1_wp2("advance_wp2_wp3::lhs_dp1_wp2", ngrdcol,
                     nzm); // wp2 "over-implicit" dissipation term
  Real2D rhs_dp1_wp2("advance_wp2_wp3::rhs_dp1_wp2", ngrdcol, nzm); // wp2 rhs dissipation term
  Real2D lhs_pr1_wp2("advance_wp2_wp3::lhs_pr1_wp2", ngrdcol,
                     nzm); // wp2 "over-implicit" pressure term 1
  Real2D rhs_pr1_wp2("advance_wp2_wp3::rhs_pr1_wp2", ngrdcol, nzm); // wp2 rhs pressure term 1
  Real2D rhs_bp_pr2_wp2("advance_wp2_wp3::rhs_bp_pr2_wp2", ngrdcol,
                        nzm); // wp2 bouyancy production and pressure term 2
  Real2D rhs_pr_dfsn_wp2("advance_wp2_wp3::rhs_pr_dfsn_wp2", ngrdcol,
                         nzm); // wp2 pressure diffusion term
  Real2D lhs_ac_pr2_wp2("advance_wp2_wp3::lhs_ac_pr2_wp2", ngrdcol,
                        nzm); // Accumulation terms of w'^2 and w'^2 pressure term 2
  Real2D rhs_pr3_wp2("advance_wp2_wp3::rhs_pr3_wp2", ngrdcol, nzm); // wp2 pressure term 3
  Real2D lhs_pr1_wp3("advance_wp2_wp3::lhs_pr1_wp3", ngrdcol,
                     nzt); // wp3 "over-implicit" pressure term 1
  Real2D rhs_pr1_wp3("advance_wp2_wp3::rhs_pr1_wp3", ngrdcol, nzt); // wp3 rhs pressure term 1
  Real2D rhs_bp1_pr2_wp3("advance_wp2_wp3::rhs_bp1_pr2_wp3", ngrdcol,
                         nzt); // wp3 bouyancy production 1 and pressure term 2
  Real2D rhs_pr3_wp3("advance_wp2_wp3::rhs_pr3_wp3", ngrdcol, nzt); // wp3 pressure term 3
  Real2D rhs_ta_wp3("advance_wp2_wp3::rhs_ta_wp3", ngrdcol, nzt);   // wp3 turbulent advection term
  Real2D rhs_pr_turb_wp3("advance_wp2_wp3::rhs_pr_turb_wp3", ngrdcol,
                         nzt); // wp3 pressure-turbulence correlation term !--EXPERIMENTAL--!
  Real2D rhs_pr_dfsn_wp3("advance_wp2_wp3::rhs_pr_dfsn_wp3", ngrdcol,
                         nzt); // wp3 pressure diffusion term
  Real2D lhs_ac_pr2_wp3("advance_wp2_wp3::lhs_ac_pr2_wp3", ngrdcol,
                        nzt); // Accumulation terms of w'^3 and w'^3 pressure term 2
  Real3D lhs_diff_zm("advance_wp2_wp3::lhs_diff_zm", ndiags3, ngrdcol,
                     nzm); // Completely implicit diffusion term for w'2
  Real3D lhs_diff_zm_crank("advance_wp2_wp3::lhs_diff_zm_crank", ndiags3, ngrdcol, nzm);
  Real3D lhs_ma_zm("advance_wp2_wp3::lhs_ma_zm", ndiags3, ngrdcol,
                   nzm); // Mean advection term for w'2
  Real3D lhs_diff_zt("advance_wp2_wp3::lhs_diff_zt", ndiags3, ngrdcol,
                     nzt); // Completely implicit diffusion term for w'3
  Real3D lhs_diff_zt_crank("advance_wp2_wp3::lhs_diff_zt_crank", ndiags3, ngrdcol, nzt);
  Real3D lhs_ma_zt("advance_wp2_wp3::lhs_ma_zt", ndiags3, ngrdcol,
                   nzt); // Mean advection term for w'3
  Real2D coef_wp4_implicit_zt("advance_wp2_wp3::coef_wp4_implicit_zt", ngrdcol,
                              nzt); // <w'^4>|_zt=coef_wp4_implicit_zt*<w'^2>|_zt^2 [-]
  Real2D coef_wp4_implicit("advance_wp2_wp3::coef_wp4_implicit", ngrdcol,
                           nzm); // <w'^4> = coef_wp4_implicit * <w'^2>^2        [-]
  Real2D a1_coef("advance_wp2_wp3::a1_coef", ngrdcol,
                 nzm); // a_1 (momentum levels); See eqn. 23 in `Equations for CLUBB' [-]
  Real2D a3_coef("advance_wp2_wp3::a3_coef", ngrdcol,
                 nzm); // a_3 (momentum levels); See eqn. 25 in `Equations for CLUBB' [-]
  Real2D a1_coef_zt("advance_wp2_wp3::a1_coef_zt", ngrdcol,
                    nzt); // a_1 interpolated to thermodynamic levels                    [-]
  Real2D a3_coef_zt("advance_wp2_wp3::a3_coef_zt", ngrdcol,
                    nzt); // a_3 interpolated to thermodynamic levels                    [-]
  Real2D dum_dz("advance_wp2_wp3::dum_dz", ngrdcol, nzm); // Vertical derivatives of um and vm
  Real2D dvm_dz("advance_wp2_wp3::dvm_dz", ngrdcol, nzm);
  Real3D lhs("advance_wp2_wp3::lhs", ndiags5, ngrdcol,
             ((2 * nzm) - 1)); // Implicit contributions to wp2/wp3 (band diag. matrix)
  Real2D rhs("advance_wp2_wp3::rhs", ngrdcol, ((2 * nzm) - 1)); // RHS of band matrix
  Real1D C_wp3_pr_tp("advance_wp2_wp3::C_wp3_pr_tp",
                     ngrdcol); // CLUBB tunable parameter C_wp3_pr_tp
  Real2D Kw1_zm("advance_wp2_wp3::Kw1_zm", ngrdcol,
                nzm); // Eddy diffusivity coefficient, momentum levels [m2/s]
  Real2D Kw8_zt("advance_wp2_wp3::Kw8_zt", ngrdcol,
                nzt); // Eddy diffusivity coefficient, thermo. levels [m2/s]
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto pdf_implicit_coefs_terms_coef_wp4_implicit =
      pdf_implicit_coefs_terms.coef_wp4_implicit;
  const auto nu_vert_res_dep_nu1 = nu_vert_res_dep.nu1;
  const auto nu_vert_res_dep_nu8 = nu_vert_res_dep.nu8;
  const auto gr_invrs_dzm = gr.invrs_dzm;
  const auto gr_weights_zm2zt = gr.weights_zm2zt;
  const auto gr_weights_zt2zm = gr.weights_zt2zm;
  const auto gr_invrs_dzt = gr.invrs_dzt;
  const auto gr_grid_dir = gr.grid_dir;
  const auto gr_zt = gr.zt;

  Kokkos::deep_copy(wp2_zt, zm2zt_api(nzm, nzt, ngrdcol, gr, wp2, w_tol_sqd));
  Kokkos::deep_copy(wp3_zm, zt2zm_api(nzm, nzt, ngrdcol, gr, wp3));
  Skx_func(nzt, ngrdcol, wp2_zt, wp3, w_tol, clubb_params, Skw_zt);
  Skx_func(nzm, ngrdcol, wp2, wp3_zm, w_tol, clubb_params, Skw_zm);
  calc_wp3_on_wp2(nzm, nzt, ngrdcol, gr, wp2, wp3, wp3_on_wp2, wp3_on_wp2_zt);
  // -----------------------------------------------------------------------
  // Compute the a3 coefficient (formula 25 in `Equations for CLUBB')
  // Note: a3 has been modified because the wp3 turbulent advection term is
  // now discretized on its own. This removes the "- 3" from the end.
  // a3_coef = 3.0_core_rknd * sigma_sqd_w*sigma_sqd_w  &
  // + 6.0_core_rknd*(1.0_core_rknd-sigma_sqd_w)*sigma_sqd_w  &
  // + (1.0_core_rknd-sigma_sqd_w)*(1.0_core_rknd-sigma_sqd_w)
  // This is a simplified version of the formula above.
  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "advance_wp2_wp3::cells_1", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          a3_coef(i - 1, k - 1) = ((-(real_t(2.)) * ((real_t(1.) - sigma_sqd_w(i - 1, k - 1)) *
                                                     (real_t(1.) - sigma_sqd_w(i - 1, k - 1)))) +
                                   real_t(3.0));
          // We found we obtain fewer spikes in wp3 when we clip a3 to be no greater
          // than -1.4 -dschanen 4 Jan 2011
          // a3_coef = max( a3_coef, -1.4_core_rknd ) ! Known magic number
          a3_coef(i - 1, k - 1) = (a3_coef(i - 1, k - 1) > clubb_params(i - 1, ia3_coef_min)
                                       ? a3_coef(i - 1, k - 1)
                                       : clubb_params(i - 1, ia3_coef_min));
        });
  }
  Kokkos::deep_copy(a3_coef_zt, zm2zt_api(nzm, nzt, ngrdcol, gr, a3_coef));
  if (stats.flags().l_sample) {
    f2k::stats_update("a3_coef_zt", a3_coef_zt, stats);
    f2k::stats_update("a3_coef", a3_coef, stats);
  }
  // Define tauw
  // tauw3t = tau_zt
  // / ( 1.
  // + 3.0_core_rknd * max(
  // min(1.-(mixt_frac-0.01_core_rknd)/(0.05_core_rknd-0.01_core_rknd)
  // ,1.)
  // ,0.)
  // + 3.0_core_rknd * max(
  // min(1.-(mixt_frac-0.99_core_rknd)/(0.95_core_rknd-0.99_core_rknd)
  // ,1.)
  // ,0.)
  // )
  // do k=1,gr%nz
  // Skw = abs( wp3(k)/max(wp2(k),1.e-8)**1.5_core_rknd )
  // Skw = min( 5.0_core_rknd, Skw )
  // tauw3t(k) = tau_zt(k) / ( 0.005_core_rknd*Skw**4 + one )
  // end do
  if ((l_crank_nich_diff && l_use_tke_in_wp2_wp3_K_dfsn)) {
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << err_info.err_header_global << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "The l_crank_nich_diff flag and l_use_tke_in_wp2_wp3_K_dfsn " << ' ';
      diagnostic << "flags cannot currently be used together." << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  // Vince Larson added code to make C11 function of Skw. 13 Mar 2005
  // If this code is used, C11 is no longer relevant, i.e. constants
  // are hardwired.
  if (l_use_C11_Richardson) {
    Kokkos::deep_copy(C11_Skw_fnc,
                      zm2zt_api(nzm, nzt, ngrdcol, gr, Cx_fnc_Richardson, zero_threshold));
  } else {
    {
      const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "advance_wp2_wp3::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
            real_t C11{};
            real_t C11b{};
            real_t C11c{};
            C11 = clubb_params(i - 1, iC11);
            C11b = clubb_params(i - 1, iC11b);
            C11c = clubb_params(i - 1, iC11c);
            // Calculate C_{1} and C_{11} as functions of skewness of w.
            // The if..then here is only for computational efficiency -dschanen 2 Sept 08
            if ((Kokkos::abs((C11 - C11b)) > ((Kokkos::abs((C11 + C11b)) * eps) / 2))) {
              C11_Skw_fnc(i - 1, k - 1) =
                  (C11b +
                   ((C11 - C11b) * Kokkos::exp((-(one_half) * ((Skw_zt(i - 1, k - 1) / C11c) *
                                                               (Skw_zt(i - 1, k - 1) / C11c))))));
            } else {
              C11_Skw_fnc(i - 1, k - 1) = C11b;
            }
          });
    }
  }
  // l_use_C11_Richardson
  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "advance_wp2_wp3::cells_3", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          real_t C1{};
          real_t C1b{};
          real_t C1c{};
          C1 = clubb_params(i - 1, iC1);
          C1b = clubb_params(i - 1, iC1b);
          C1c = clubb_params(i - 1, iC1c);
          // The if..then here is only for computational efficiency -dschanen 2 Sept 08
          if ((Kokkos::abs((C1 - C1b)) > ((Kokkos::abs((C1 + C1b)) * eps) / 2))) {
            C1_Skw_fnc(i - 1, k - 1) =
                (C1b + ((C1 - C1b) * Kokkos::exp((-(one_half) * ((Skw_zm(i - 1, k - 1) / C1c) *
                                                                 (Skw_zm(i - 1, k - 1) / C1c))))));
          } else {
            C1_Skw_fnc(i - 1, k - 1) = C1b;
          }
        });
  }
  if (l_damp_wp2_using_em) {
    // Insert 1/3 here to account for the fact that in the dissipation term,
    // (2/3)*em = (2/3)*(1/2)*(wp2+up2+vp2).  Then we can insert wp2, up2,
    // and vp2 directly into the dissipation subroutines without prefixing them by (1/3).
    {
      const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "advance_wp2_wp3::cells_4", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
            C1_Skw_fnc(i - 1, k - 1) = (one_third * C1_Skw_fnc(i - 1, k - 1));
          });
    }
  }
  // Set C16_fnc based on Richardson_num
  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "advance_wp2_wp3::cells_5", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          C16_fnc(i - 1, k - 1) = Cx_fnc_Richardson(i - 1, k - 1);
        });
  }
  if (f2k::clubb_at_least_debug_level_api(0)) {
    // Detect the Native per-cell assertions on device, then print failures in
    // Native level/column order. Only integer counts and failure flags reduce.
    {
      Kokkos::View<int**, Kokkos::LayoutLeft> invalid("advance_wp2_wp3::C11_Skw_fnc_invalid",
                                                      ngrdcol, nzt);
      const auto err_code = err_info.err_code;
      int count = 0;
      Kokkos::parallel_reduce(
          "advance_wp2_wp3::C11_Skw_fnc_check", ngrdcol,
          KOKKOS_LAMBDA(const int column, int& bad) {
            for (int level = 0; level < nzt; ++level) {
              invalid(column, level) = 0;
              if (C11_Skw_fnc(column, level) > one || C11_Skw_fnc(column, level) < real_t(0.0)) {
                invalid(column, level) = 1;
                err_code(column) = clubb_fatal_error;
                bad += 1;
              }
            }
          },
          count);
      if (count != 0) {
        const auto invalid_host = to_host_space(invalid);
        const auto header_host = to_host_space(err_info.err_header);
        for (int level = 0; level < nzt; ++level)
          for (int column = 0; column < ngrdcol; ++column) {
            if (invalid_host(column, level) != 0) {
              std::ostringstream diagnostic;
              diagnostic.write(header_host.data() + column * header_host.extent(0),
                               header_host.extent(0));
              diagnostic << '\n'
                         << "The C11_Skw_fnc is outside the valid range for this variable" << '\n';
              std::cerr << diagnostic.str();
            }
          }
      }
    }
    // Detect the Native per-cell assertions on device, then print failures in
    // Native level/column order. Only integer counts and failure flags reduce.
    {
      Kokkos::View<int**, Kokkos::LayoutLeft> invalid("advance_wp2_wp3::C16_fnc_invalid", ngrdcol,
                                                      nzt);
      const auto err_code = err_info.err_code;
      int count = 0;
      Kokkos::parallel_reduce(
          "advance_wp2_wp3::C16_fnc_check", ngrdcol,
          KOKKOS_LAMBDA(const int column, int& bad) {
            for (int level = 0; level < nzt; ++level) {
              invalid(column, level) = 0;
              if (C16_fnc(column, level) > one || C16_fnc(column, level) < real_t(0.0)) {
                invalid(column, level) = 1;
                err_code(column) = clubb_fatal_error;
                bad += 1;
              }
            }
          },
          count);
      if (count != 0) {
        const auto invalid_host = to_host_space(invalid);
        const auto header_host = to_host_space(err_info.err_header);
        for (int level = 0; level < nzt; ++level)
          for (int column = 0; column < ngrdcol; ++column) {
            if (invalid_host(column, level) != 0) {
              std::ostringstream diagnostic;
              diagnostic.write(header_host.data() + column * header_host.extent(0),
                               header_host.extent(0));
              diagnostic << '\n'
                         << "The C16_fnc is outside the valid range for this variable" << '\n';
              std::cerr << diagnostic.str();
            }
          }
      }
    }
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      return;
    }
  }
  if (stats.flags().l_sample) {
    f2k::stats_update("C11_Skw_fnc", C11_Skw_fnc, stats);
    f2k::stats_update("C1_Skw_fnc", C1_Skw_fnc, stats);
  }
  // Define the Coefficent of Eddy Diffusivity for the wp2 and wp3.
  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "advance_wp2_wp3::cells_6", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Kw1 is used for wp2, which is located on momentum levels.
          // Kw1 is located on thermodynamic levels.
          // Kw1 = c_K1 * Kh_zt ! c_K1 variable is commented out
          Kw1(i - 1, k - 1) = (clubb_params(i - 1, ic_K1) * Kh_zt(i - 1, k - 1));
        });
  }
  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "advance_wp2_wp3::cells_7", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Kw8 is used for wp3, which is located on thermodynamic levels.
          // Kw8 is located on momentum levels.
          // Note: Kw8 is usually defined to be 1/2 of Kh_zm.
          // Kw8 = c_K8 * Kh_zm ! c_K8 variable is commented out
          Kw8(i - 1, k - 1) = (clubb_params(i - 1, ic_K8) * Kh_zm(i - 1, k - 1));
        });
  }
  if (!(l_explicit_turbulent_adv_wp3)) {
    if (((iiPDF_type == iiPDF_new) || (iiPDF_type == iiPDF_new_hybrid))) {
      // Unpack coef_wp4_implicit from pdf_implicit_coefs_terms.
      // Since PDF parameters and the resulting implicit coefficients and
      // explicit terms are calculated on thermodynamic levels, the <w'^4>
      // implicit coefficient needs to be unpacked as coef_wp4_implicit_zt.
      {
        const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
        const int extent_1 =
            ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
        Kokkos::parallel_for(
            "advance_wp2_wp3::cells_8", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
              const int i = (1) + ((cell / (1)) % extent_1) * (1);
              const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
              coef_wp4_implicit_zt(i - 1, k - 1) =
                  pdf_implicit_coefs_terms_coef_wp4_implicit(i - 1, k - 1);
            });
      }
      // The values of <w'^4> are located on momentum levels.  Interpolate
      // coef_wp4_implicit_zt to momentum levels as coef_wp4_implicit.  The
      // discretization diagram is found in the description section of
      // function wp3_term_ta_new_pdf_lhs below.  These values are always
      // positive.
      Kokkos::deep_copy(coef_wp4_implicit,
                        zt2zm_api(nzm, nzt, ngrdcol, gr, coef_wp4_implicit_zt, zero_threshold));
      // Set the value of coef_wp4_implicit to 0 at the lower boundary and at
      // the upper boundary.  This sets the value of <w'^4> to 0 at the lower
      // and upper boundaries.
      Kokkos::deep_copy(Kokkos::subview(coef_wp4_implicit, Kokkos::ALL, 1 - 1), zero);
      Kokkos::deep_copy(Kokkos::subview(coef_wp4_implicit, Kokkos::ALL, nzm - 1), zero);
      if (stats.flags().l_sample) {
        f2k::stats_update("coef_wp4_implicit", coef_wp4_implicit, stats);
      }
    } else if ((iiPDF_type == iiPDF_ADG1)) {
      // Define a_1 and a_3 (both are located on momentum levels).
      // They are variables that are both functions of sigma_sqd_w (where
      // sigma_sqd_w is located on momentum levels).
      {
        const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
        const int extent_1 =
            ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
        Kokkos::parallel_for(
            "advance_wp2_wp3::cells_9", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
              const int i = (1) + ((cell / (1)) % extent_1) * (1);
              const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
              a1_coef(i - 1, k - 1) = (one / (one - sigma_sqd_w(i - 1, k - 1)));
            });
      }
      // Interpolate a_1 from momentum levels to thermodynamic levels.  This
      // will be used for the w'^3 turbulent advection (ta) term.
      Kokkos::deep_copy(a1_coef_zt, zm2zt_api(nzm, nzt, ngrdcol, gr, a1_coef, zero_threshold));
      // Positive def. quantity
    }
    // iiPDF_type
  }
  // .not. l_explicit_turbulent_adv_wp3
  // Not using pressure term, set to 0
  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "advance_wp2_wp3::cells_10", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          rhs_pr3_wp3(i - 1, k - 1) = zero;
        });
  }
  // Initiaize some terms to zero
  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    const int extent_2 =
        ((1) > 0 ? (1) <= (ndiags5) : (1) >= (ndiags5)) ? ((ndiags5) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "advance_wp2_wp3::cells_11", extent_0 * extent_1 * extent_2, KOKKOS_LAMBDA(const int cell) {
          const int b = (1) + ((cell / (1)) % extent_2) * (1);
          const int i = (1) + ((cell / (1 * extent_2)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_2 * extent_1)) % extent_0) * (1);
          wp3_term_ta_lhs_result(b - 1, i - 1, k - 1) = zero;
          wp3_pr3_lhs(b - 1, i - 1, k - 1) = zero;
        });
  }
  Kokkos::deep_copy(Kw1_zm, zt2zm_api(nzm, nzt, ngrdcol, gr, Kw1, zero));
  Kokkos::deep_copy(Kw8_zt, zm2zt_api(nzm, nzt, ngrdcol, gr, Kw8, zero));
  // Experimental bouyancy term
  // Experimental term from CLUBB TRAC ticket #411
  // Compute the vertical derivative of the u and v winds
  if (!(l_use_tke_in_wp3_pr_turb_term)) {
    Kokkos::deep_copy(dum_dz, ddzt(nzm, nzt, ngrdcol, gr, um));
    Kokkos::deep_copy(dvm_dz, ddzt(nzm, nzt, ngrdcol, gr, vm));
  }
  Kokkos::deep_copy(em_smth, zm2zt2zm(nzm, nzt, ngrdcol, gr, em));
  Kokkos::deep_copy(wp2_smth, zm2zt2zm(nzm, nzt, ngrdcol, gr, wp2));
  // Calculate term
  wp3_term_pr_turb_rhs(nzm, nzt, ngrdcol, gr,                                      // In
                       Kokkos::subview(clubb_params, Kokkos::ALL, iC_wp3_pr_turb), // In
                       Kh_zt, wpthvp,                                              // In
                       dum_dz, dvm_dz,                                             // In
                       upwp, vpwp,                                                 // In
                       thv_ds_zt,                                                  // In
                       rho_ds_zm, invrs_rho_ds_zt,                                 // In
                       em_smth, wp2_smth,                                          // In
                       rhs_pr_turb_wp3,                                            // Out
                       l_use_tke_in_wp3_pr_turb_term);                             // In

  wp3_term_pr_dfsn_rhs(nzm, nzt, ngrdcol, gr,                                      // In
                       Kokkos::subview(clubb_params, Kokkos::ALL, iC_wp3_pr_dfsn), // In
                       rho_ds_zm, invrs_rho_ds_zt,                                 // In
                       wp2up2, wp2vp2, wp4,                                        // In
                       up2, vp2, wp2,                                              // In
                       rhs_pr_dfsn_wp3);                                           // Out

  wp2_term_pr_dfsn_rhs(nzm, nzt, ngrdcol, gr,                                      // In
                       Kokkos::subview(clubb_params, Kokkos::ALL, iC_wp2_pr_dfsn), // In
                       rho_ds_zt, invrs_rho_ds_zm,                                 // In
                       wpup2, wpvp2, wp3,                                          // In
                       rhs_pr_dfsn_wp2);                                           // Out

  // This part handles the wp2 equation terms.
  f2k::diffusion_zm_lhs(nzm, nzt, ngrdcol, gr, Kw1, Kw1_zm, nu_vert_res_dep_nu1, // In
                        invrs_rho_ds_zm, rho_ds_zt,                              // In
                        lhs_diff_zm);                                            // Out

  // This part handles the wp3 equation terms.
  f2k::diffusion_zt_lhs(nzm, nzt, ngrdcol, gr, Kw8, Kw8_zt, nu_vert_res_dep_nu8, // In
                        invrs_rho_ds_zt, rho_ds_zm,                              // In
                        lhs_diff_zt);                                            // Out

  // Calculate RHS eddy diffusion terms for w'2 and w'3
  if (l_crank_nich_diff) {
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "advance_wp2_wp3::cells_12", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            lhs_diff_zm_crank(1 - 1, i - 1, k - 1) = (lhs_diff_zm(1 - 1, i - 1, k - 1) * one_half);
            lhs_diff_zm_crank(2 - 1, i - 1, k - 1) = (lhs_diff_zm(2 - 1, i - 1, k - 1) * one_half);
            lhs_diff_zm_crank(3 - 1, i - 1, k - 1) = (lhs_diff_zm(3 - 1, i - 1, k - 1) * one_half);
          });
    }
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "advance_wp2_wp3::cells_13", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            lhs_diff_zt_crank(1 - 1, i - 1, k - 1) =
                ((lhs_diff_zt(1 - 1, i - 1, k - 1) * clubb_params(i - 1, iC12)) * one_half);
            lhs_diff_zt_crank(2 - 1, i - 1, k - 1) =
                ((lhs_diff_zt(2 - 1, i - 1, k - 1) * clubb_params(i - 1, iC12)) * one_half);
            lhs_diff_zt_crank(3 - 1, i - 1, k - 1) =
                ((lhs_diff_zt(3 - 1, i - 1, k - 1) * clubb_params(i - 1, iC12)) * one_half);
          });
    }
  }
  // Calculate "over-implicit" pressure terms for w'2 and w'3
  wp2_term_pr1_rhs(nzm, ngrdcol, gr, Kokkos::subview(clubb_params, Kokkos::ALL, iC4), // In
                   up2, vp2, invrs_tau_C4_zm,                                         // In
                   rhs_pr1_wp2);                                                      // Out

  // Note:  An "over-implicit" weighted time step is applied to the  term.
  // A weighting factor of greater than 1 may be used to make the
  // term more numerically stable (see note below for w'^3 RHS
  // turbulent advection (ta) term).
  wp2_term_pr1_lhs(nzm, ngrdcol, gr,                                                 // In
                   Kokkos::subview(clubb_params, Kokkos::ALL, iC4), invrs_tau_C4_zm, // In
                   lhs_pr1_wp2);                                                     // Out

  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "advance_wp2_wp3::cells_14", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          C_wp3_pr_tp(i - 1) = one;
        });
  }
  // Calculate turbulent production terms of w'^3
  wp3_term_tp_lhs(nzm, nzt, ngrdcol, gr, C_wp3_pr_tp, // In
                  wp2, rho_ds_zm, invrs_rho_ds_zt,    // In
                  lhs_adv_tp_wp3);                    // Out

  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "advance_wp2_wp3::cells_15", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          C_wp3_pr_tp(i - 1) = -(clubb_params(i - 1, iC_wp3_pr_tp));
        });
  }
  // Calculate pressure damping of turbulent production of w'^3
  wp3_term_tp_lhs(nzm, nzt, ngrdcol, gr, C_wp3_pr_tp, // In
                  wp2, rho_ds_zm, invrs_rho_ds_zt,    // In
                  lhs_pr_tp_wp3);                     // Out

  // Sum contributions to turbulent production from standard term & damping
  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    const int extent_2 =
        ((1) > 0 ? (1) <= (ndiags2) : (1) >= (ndiags2)) ? ((ndiags2) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "advance_wp2_wp3::cells_16", extent_0 * extent_1 * extent_2, KOKKOS_LAMBDA(const int cell) {
          const int b = (1) + ((cell / (1)) % extent_2) * (1);
          const int i = (1) + ((cell / (1 * extent_2)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_2 * extent_1)) % extent_0) * (1);
          lhs_tp_wp3(b - 1, i - 1, k - 1) =
              (lhs_adv_tp_wp3(b - 1, i - 1, k - 1) + lhs_pr_tp_wp3(b - 1, i - 1, k - 1));
        });
  }
  // Calculate pressure terms 1 for w'^3
  wp3_term_pr1_lhs(nzt, ngrdcol, gr,                                 // In
                   Kokkos::subview(clubb_params, Kokkos::ALL, iC8),  // In
                   Kokkos::subview(clubb_params, Kokkos::ALL, iC8b), // In
                   invrs_tau_wp3_zt, Skw_zt,                         // In
                   l_damp_wp3_Skw_squared,                           // In
                   lhs_pr1_wp3);                                     // Out

  // Calculate dissipation terms 1 for w'^2
  wp2_term_dp1_lhs(nzm, ngrdcol, gr,            // In
                   C1_Skw_fnc, invrs_tau_C1_zm, // In
                   lhs_dp1_wp2);                // Out

  // Calculate buoyancy production of w'^2 and w'^2 pressure term 2
  wp2_terms_bp_pr2_rhs(nzm, ngrdcol, gr,                                       // In
                       Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_buoy), // In
                       thv_ds_zm, wpthvp,                                      // In
                       rhs_bp_pr2_wp2);                                        // Out

  // Calculate pressure terms 3 for w'^2
  wp2_term_pr3_rhs(nzm, nzt, ngrdcol, gr,                                  // In
                   Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_shr),  // In
                   Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_buoy), // In
                   thv_ds_zm, wpthvp, upwp,                                // In
                   um, vpwp, vm,                                           // In
                   rhs_pr3_wp2);                                           // Out

  // Calculate dissipation terms 1 for w'^2
  wp2_term_dp1_rhs(nzm, ngrdcol, gr, C1_Skw_fnc,         // In
                   invrs_tau_C1_zm, w_tol_sqd, up2, vp2, // In
                   l_damp_wp2_using_em,                  // In
                   rhs_dp1_wp2);                         // Out

  // Calculate buoyancy production of w'^3 and w'^3 pressure term 2
  wp3_terms_bp1_pr2_rhs(nzt, ngrdcol, gr, C11_Skw_fnc, // In
                        thv_ds_zt, wp2thvp,            // In
                        rhs_bp1_pr2_wp3);              // Out

  // Calculate pressure terms 1 for w'^3
  wp3_term_pr1_rhs(nzt, ngrdcol, gr,                                 // In
                   Kokkos::subview(clubb_params, Kokkos::ALL, iC8),  // In
                   Kokkos::subview(clubb_params, Kokkos::ALL, iC8b), // In
                   invrs_tau_wp3_zt, Skw_zt, wp3,                    // In
                   l_damp_wp3_Skw_squared,                           // In
                   rhs_pr1_wp3);                                     // Out

  if (l_explicit_turbulent_adv_wp3) {
    // The turbulent advection term is being solved explicitly.
    // The w'^3 turbulent advection term is being solved explicitly.
    // The turbulent advection stats code is still set up in two parts,
    // so call stat_begin_update_pt.  The implicit portion of the stat,
    // which has a value of 0, will still be called later.  Since
    // stat_begin_update_pt automatically subtracts the value sent in,
    // reverse the sign on the input value.
    wp3_term_ta_explicit_rhs(nzm, nzt, ngrdcol, gr,           // In
                             wp4, rho_ds_zm, invrs_rho_ds_zt, // In
                             rhs_ta_wp3);                     // Out

  } else {
    // The turbulent advection term is being solved implicitly.
    // The turbulent advection term is being solved implicitly. See note above
    if ((iiPDF_type == iiPDF_ADG1)) {
      // The ADG1 PDF is used.
      wp3_term_ta_ADG1_lhs(nzm, nzt, ngrdcol, gr,      // In
                           wp2, a1_coef, a1_coef_zt,   // In
                           a3_coef, a3_coef_zt,        // In
                           wp3_on_wp2, rho_ds_zm,      // In
                           rho_ds_zt, invrs_rho_ds_zt, // In
                           l_standard_term_ta,         // In
                           l_partial_upwind_wp3,       // In
                           wp3_term_ta_lhs_result);    // Out

    } else if (((iiPDF_type == iiPDF_new) || (iiPDF_type == iiPDF_new_hybrid))) {
      // The new PDF or the new hybrid PDF is used.
      // Calculate terms
      wp3_term_ta_new_pdf_lhs(nzm, nzt, ngrdcol, gr, coef_wp4_implicit, // In
                              wp2, rho_ds_zm, invrs_rho_ds_zt,          // In
                              lhs_ta_wp3);                              // Out

      // Add terms
      {
        const int extent_0 =
            ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
        const int extent_1 =
            ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
        Kokkos::parallel_for(
            "advance_wp2_wp3::cells_17", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
              const int i = (1) + ((cell / (1)) % extent_1) * (1);
              const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
              wp3_term_ta_lhs_result(2 - 1, i - 1, k - 1) = lhs_ta_wp3(1 - 1, i - 1, k - 1);
              wp3_term_ta_lhs_result(4 - 1, i - 1, k - 1) = lhs_ta_wp3(2 - 1, i - 1, k - 1);
            });
      }
    }
    // iiPDF_type
  }
  // l_explicit_turbulent_adv_wp3
  // Compute the explicit portion of the w'^2 and w'^3 equations.
  // Build the right-hand side vector.
  wp23_rhs(nzm, nzt, ngrdcol, gr, dt, fcor_y,                              // In
           wp3_term_ta_lhs_result,                                         // In
           lhs_diff_zm, lhs_diff_zt, lhs_diff_zm_crank, lhs_diff_zt_crank, // In
           lhs_tp_wp3, lhs_adv_tp_wp3, lhs_pr_tp_wp3,                      // In
           lhs_ta_wp3, lhs_dp1_wp2, rhs_dp1_wp2, lhs_pr1_wp2,              // In
           rhs_pr1_wp2, lhs_pr1_wp3, rhs_pr1_wp3, rhs_bp_pr2_wp2,          // In
           rhs_pr_dfsn_wp2, rhs_bp1_pr2_wp3, rhs_pr3_wp2, rhs_pr3_wp3,     // In
           rhs_ta_wp3, rhs_pr_turb_wp3, rhs_pr_dfsn_wp3,                   // In
           wp2, wp3, wpup2, wpvp2,                                         // In
           wpthvp, wp2thvp, wp2up, up2, vp2, upwp,                         // In
           C11_Skw_fnc, thv_ds_zm, thv_ds_zt,                              // In
           lhs_splat_wp2, lhs_splat_wp3,                                   // In
           clubb_params,                                                   // In
           iiPDF_type,                                                     // In
           l_use_tke_in_wp2_wp3_K_dfsn,                                    // In
           l_ho_nontrad_coriolis,                                          // In
           stats,                                                          // InOut
           rhs);                                                           // Out

  // Calculated mean advection term for w'2
  f2k::term_ma_zm_lhs(nzm, nzt, ngrdcol, wm_zm, // In
                      gr_invrs_dzm, gr_weights_zm2zt,
                      lhs_ma_zm); // Out

  // Calculated mean advection term for w'3
  f2k::term_ma_zt_lhs(nzm, nzt, ngrdcol, wm_zt, gr_weights_zt2zm, // In
                      gr_invrs_dzt, gr_invrs_dzm,                 // In
                      l_upwind_xm_ma, gr_grid_dir,                // In
                      lhs_ma_zt);                                 // Out

  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    const int extent_2 =
        ((1) > 0 ? (1) <= (ndiags3) : (1) >= (ndiags3)) ? ((ndiags3) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "advance_wp2_wp3::cells_18", extent_0 * extent_1 * extent_2, KOKKOS_LAMBDA(const int cell) {
          const int b = (1) + ((cell / (1)) % extent_2) * (1);
          const int i = (1) + ((cell / (1 * extent_2)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_2 * extent_1)) % extent_0) * (1);
          lhs_diff_zt(b - 1, i - 1, k - 1) =
              (lhs_diff_zt(b - 1, i - 1, k - 1) * clubb_params(i - 1, iC12));
        });
  }
  if (l_crank_nich_diff) {
    // Using a Crank-Nicholson time step for diffusion terms
    // Modify diffusion terms
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "advance_wp2_wp3::cells_19", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            lhs_diff_zm(1 - 1, i - 1, k - 1) = (lhs_diff_zm(1 - 1, i - 1, k - 1) * real_t(0.5));
            lhs_diff_zm(2 - 1, i - 1, k - 1) = (lhs_diff_zm(2 - 1, i - 1, k - 1) * real_t(0.5));
            lhs_diff_zm(3 - 1, i - 1, k - 1) = (lhs_diff_zm(3 - 1, i - 1, k - 1) * real_t(0.5));
          });
    }
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "advance_wp2_wp3::cells_20", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            lhs_diff_zt(1 - 1, i - 1, k - 1) = (lhs_diff_zt(1 - 1, i - 1, k - 1) * real_t(0.5));
            lhs_diff_zt(2 - 1, i - 1, k - 1) = (lhs_diff_zt(2 - 1, i - 1, k - 1) * real_t(0.5));
            lhs_diff_zt(3 - 1, i - 1, k - 1) = (lhs_diff_zt(3 - 1, i - 1, k - 1) * real_t(0.5));
          });
    }
  }
  // Calculate turbulent advection terms for wp2
  wp2_term_ta_lhs(nzm, nzt, ngrdcol, gr,      // In
                  rho_ds_zt, invrs_rho_ds_zm, // In
                  lhs_ta_wp2);                // Out

  // Calculate accumulation terms of w'^2 and w'^2 pressure term 2
  wp2_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr,                                        // In
                       Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_shr), wm_zt, // In
                       lhs_ac_pr2_wp2);                                              // Out

  // Calculate accumulation terms of w'^3 and w'^3 pressure terms 2
  wp3_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr, C11_Skw_fnc, wm_zm, // In
                       lhs_ac_pr2_wp3);                           // Out

  // Compute the implicit portion of the w'^2 and w'^3 equations.
  // Build the left-hand side matrix.
  wp23_lhs(nzm, nzt, ngrdcol, dt,                                  // In
           wp3_term_ta_lhs_result,                                 // In
           lhs_diff_zm, lhs_diff_zt, lhs_ma_zm,                    // In
           lhs_ma_zt, lhs_ta_wp2,                                  // In
           lhs_tp_wp3,                                             // In
           lhs_ac_pr2_wp2, lhs_ac_pr2_wp3, lhs_dp1_wp2,            // In
           lhs_pr1_wp3, lhs_pr1_wp2, lhs_splat_wp2, lhs_splat_wp3, // In
           lhs);                                                   // Out

  if (l_lmm_stepping) {
    {
      const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "advance_wp2_wp3::cells_21", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
            wp2_old(i - 1, k - 1) = wp2(i - 1, k - 1);
          });
    }
    {
      const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "advance_wp2_wp3::cells_22", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
            wp3_old(i - 1, k - 1) = wp3(i - 1, k - 1);
          });
    }
  }
  // l_lmm_stepping
  // Solve semi-implicitly
  wp23_solve(nzm, nzt, ngrdcol, gr, dt,                     // In
             lhs, rhs,                                      // InOut
             lhs_ma_zm, lhs_dp1_wp2, lhs_diff_zm,           // In
             lhs_ta_wp2, lhs_pr1_wp2, lhs_pr1_wp3,          // In
             lhs_diff_zt, lhs_adv_tp_wp3, lhs_pr_tp_wp3,    // In
             wp3_pr3_lhs, lhs_ma_zt,                        // In
             wp3_term_ta_lhs_result,                        // In
             wm_zm, wm_zt,                                  // In
             sfc_elevation, C11_Skw_fnc,                    // In
             rho_ds_zm,                                     // In
             upwp, vpwp, wprtp, wpthlp, rtp2, thlp2,        // In
             clubb_params,                                  // In
             penta_solve_method,                            // In
             fill_holes_type,                               // In
             l_min_wp2_from_corr_wx,                        // In
             l_use_tke_in_wp2_wp3_K_dfsn,                   // In
             l_use_wp3_lim_with_smth_Heaviside,             // In
             l_wp2_fill_holes_tke,                          // In
             l_implemented,                                 // In
             stats,                                         // InOut
             up2, vp2, wp2, wp3, wp3_zm, wp2_zt, err_info); // InOut
  if (l_lmm_stepping) {
    {
      const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "advance_wp2_wp3::cells_23", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
            wp2(i - 1, k - 1) = (one_half * (wp2_old(i - 1, k - 1) + wp2(i - 1, k - 1)));
          });
    }
    {
      const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "advance_wp2_wp3::cells_24", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
            wp3(i - 1, k - 1) = (one_half * (wp3_old(i - 1, k - 1) + wp3(i - 1, k - 1)));
          });
    }
  }
  // l_lmm_stepping
  // TODO: Native post-solve wp2/wp3 sponge blocks remain outside the supported
  // JAX/Kokkos scope. clubb_driver rejects sponge-enabled cases before allocation;
  // no placeholder flags, profiles or global sponge owner are introduced here.
  if (f2k::clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << err_info.err_header_global << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Error in advance_wp2_wp3" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "intent(in)" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "gr%zt(1,:) = " << ' ';
        {
          const auto field_host = to_host_space(gr_zt);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "dt = " << ' ';
        diagnostic << dt << ' ';
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "sfc_elevation = " << ' ';
        {
          const auto field_host = to_host_space(sfc_elevation);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "sigma_sqd_w = " << ' ';
        {
          const auto field_host = to_host_space(sigma_sqd_w);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wm_zm = " << ' ';
        {
          const auto field_host = to_host_space(wm_zm);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wm_zt = " << ' ';
        {
          const auto field_host = to_host_space(wm_zt);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wpup2 = " << ' ';
        {
          const auto field_host = to_host_space(wpup2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wpvp2 = " << ' ';
        {
          const auto field_host = to_host_space(wpvp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp2up2 = " << ' ';
        {
          const auto field_host = to_host_space(wp2up2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp2vp2 = " << ' ';
        {
          const auto field_host = to_host_space(wp2vp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp4 = " << ' ';
        {
          const auto field_host = to_host_space(wp4);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wpthvp = " << ' ';
        {
          const auto field_host = to_host_space(wpthvp);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp2thvp = " << ' ';
        {
          const auto field_host = to_host_space(wp2thvp);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp2up = " << ' ';
        {
          const auto field_host = to_host_space(wp2up);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "um = " << ' ';
        {
          const auto field_host = to_host_space(um);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "vm = " << ' ';
        {
          const auto field_host = to_host_space(vm);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "upwp = " << ' ';
        {
          const auto field_host = to_host_space(upwp);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "vpwp = " << ' ';
        {
          const auto field_host = to_host_space(vpwp);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "up2 = " << ' ';
        {
          const auto field_host = to_host_space(up2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "vp2 = " << ' ';
        {
          const auto field_host = to_host_space(vp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "em = " << ' ';
        {
          const auto field_host = to_host_space(em);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Kh_zm = " << ' ';
        {
          const auto field_host = to_host_space(Kh_zm);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Kh_zt = " << ' ';
        {
          const auto field_host = to_host_space(Kh_zt);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "invrs_tau_C4 zm = " << ' ';
        {
          const auto field_host = to_host_space(invrs_tau_C4_zm);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "invrs_tau_wp3_zt = " << ' ';
        {
          const auto field_host = to_host_space(invrs_tau_wp3_zt);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Skw_zm = " << ' ';
        {
          const auto field_host = to_host_space(Skw_zm);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Skw_zt = " << ' ';
        {
          const auto field_host = to_host_space(Skw_zt);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "mixt_frac = " << ' ';
        {
          const auto field_host = to_host_space(mixt_frac);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "a3_coef = " << ' ';
        {
          const auto field_host = to_host_space(a3_coef);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "a3_coef_zt = " << ' ';
        {
          const auto field_host = to_host_space(a3_coef_zt);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp3_on_wp2 = " << ' ';
        {
          const auto field_host = to_host_space(wp3_on_wp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "invrs_tau_C1_zm = " << ' ';
        {
          const auto field_host = to_host_space(invrs_tau_C1_zm);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "rho_ds_zm = " << ' ';
        {
          const auto field_host = to_host_space(rho_ds_zm);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "rho_ds_zt = " << ' ';
        {
          const auto field_host = to_host_space(rho_ds_zt);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "invrs_rho_ds_zm = " << ' ';
        {
          const auto field_host = to_host_space(invrs_rho_ds_zm);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "invrs_rho_ds_zt = " << ' ';
        {
          const auto field_host = to_host_space(invrs_rho_ds_zt);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "thv_ds_zm = " << ' ';
        {
          const auto field_host = to_host_space(thv_ds_zm);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "thv_ds_zt = " << ' ';
        {
          const auto field_host = to_host_space(thv_ds_zt);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Cx_fnc_Richardson = " << ' ';
        {
          const auto field_host = to_host_space(Cx_fnc_Richardson);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "lhs_splat_wp2 = " << ' ';
        {
          const auto field_host = to_host_space(lhs_splat_wp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "lhs_splat_wp3 = " << ' ';
        {
          const auto field_host = to_host_space(lhs_splat_wp3);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wprtp = " << ' ';
        {
          const auto field_host = to_host_space(wprtp);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wpthlp = " << ' ';
        {
          const auto field_host = to_host_space(wpthlp);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "rtp2 = " << ' ';
        {
          const auto field_host = to_host_space(rtp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "thlp2 = " << ' ';
        {
          const auto field_host = to_host_space(thlp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "pdf_implicit_coefs_terms%coef_wp4_implicit = " << ' ';
        {
          const auto field_host = to_host_space(pdf_implicit_coefs_terms_coef_wp4_implicit);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "intent(in/out)" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp2_zt = " << ' ';
        {
          const auto field_host = to_host_space(wp2_zt);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp3_zm = " << ' ';
        {
          const auto field_host = to_host_space(wp3_zm);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      if (l_lmm_stepping) {
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "wp2 (pre-solve) = " << ' ';
          {
            const auto field_host = to_host_space(wp2_old);
            for (std::size_t element = 0; element < field_host.size(); ++element)
              diagnostic << field_host.data()[element] << ' ';
          }
          diagnostic << "\n" << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp2 = " << ' ';
        {
          const auto field_host = to_host_space(wp2);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      if (l_lmm_stepping) {
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "wp3 (pre-solve) = " << ' ';
          {
            const auto field_host = to_host_space(wp3_old);
            for (std::size_t element = 0; element < field_host.size(); ++element)
              diagnostic << field_host.data()[element] << ' ';
          }
          diagnostic << "\n" << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp3 = " << ' ';
        {
          const auto field_host = to_host_space(wp3);
          for (std::size_t element = 0; element < field_host.size(); ++element)
            diagnostic << field_host.data()[element] << ' ';
        }
        diagnostic << "\n" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
    }
    // fatal error
  }
  if (stats.flags().l_sample) {
    f2k::stats_update("wp2_zt", wp2_zt, stats);
    f2k::stats_update("wp3_zm", wp3_zm, stats);
  }
  return;
}

namespace {
// =============================================================================
void wp23_solve(int nzm, int nzt, int ngrdcol, const grid& gr, real_t dt,           // In
                Real3D lhs, Real2D rhs,                                             // InOut
                Real3D lhs_ma_zm, Real2D lhs_dp1_wp2, Real3D lhs_diff_zm,           // In
                Real3D lhs_ta_wp2, Real2D lhs_pr1_wp2, Real2D lhs_pr1_wp3,          // In
                Real3D lhs_diff_zt, Real3D lhs_adv_tp_wp3, Real3D lhs_pr_tp_wp3,    // In
                Real3D wp3_pr3_lhs, Real3D lhs_ma_zt,                               // In
                Real3D wp3_term_ta_lhs_result,                                      // In
                Real2D wm_zm, Real2D wm_zt,                                         // In
                Real1D sfc_elevation, Real2D C11_Skw_fnc,                           // In
                Real2D rho_ds_zm,                                                   // In
                Real2D upwp, Real2D vpwp, Real2D wprtp, Real2D wpthlp, Real2D rtp2, // In
                Real2D thlp2,                                                       // In
                Real2D clubb_params,                                                // In
                int penta_solve_method,                                             // In
                int fill_holes_type,                                                // In
                bool l_min_wp2_from_corr_wx,                                        // In
                bool l_use_tke_in_wp2_wp3_K_dfsn,                                   // In
                bool l_use_wp3_lim_with_smth_Heaviside,                             // In
                bool l_wp2_fill_holes_tke,                                          // In
                bool l_implemented,                                                 // In
                f2k::stats_type& stats,                                             // InOut
                Real2D up2, Real2D vp2, Real2D wp2, Real2D wp3, Real2D wp3_zm,      // InOut
                Real2D wp2_zt, err_info_type& err_info) {                           // InOut
  //
  // Description:
  // Decompose, and back substitute the matrix for wp2/wp3
  //
  // References:
  // _Equations for CLUBB_ section 6.3
  //------------------------------------------------------------------------

  // ----------------------- Input Variables -----------------------
  // dt: Timestep                                  [s]
  // lhs: Implicit contributions to wp2/wp3 (band diag. matrix)
  // rhs: RHS of band matrix
  // lhs_dp1_wp2: wp2 "over-implicit" dissipation term
  // lhs_pr1_wp2: wp2 "over-implicit" pressure term 1
  // ----------------------- Input/Output Variables -----------------------
  // lhs_pr1_wp3: wp3 "over-implicit" pressure term 1
  // lhs_ta_wp2: Turbulent advection terms for wp2
  // lhs_adv_tp_wp3: Turbulent production terms of w'^3 (for stats)
  // lhs_pr_tp_wp3: Pressure scrambling terms for turbulent production of w'^3 (for stats)
  // lhs_diff_zm: Completely implicit diffusion term for w'2
  // lhs_ma_zm: Mean advection term for w'2
  // lhs_diff_zt: Completely implicit diffusion term for w'3
  // lhs_ma_zt: Mean advection term for w'3
  // sfc_elevation: Elevation of ground level                 [m AMSL]
  // wm_zm: w wind component on momentum levels       [m/s]
  // rho_ds_zm: Dry, static density on momentum levels    [kg/m^3]
  // wprtp: Flux of total water mixing ratio          [m/s kg/kg]
  // wpthlp: Flux of liquid water potential temp.      [m/s K]
  // upwp: Momentum flux of u                        [m^2/s^2]
  // vpwp: Momentum flux of v                        [m^2/s^2]
  // rtp2: Variance of rt (overall)                  [kg^2/kg^2]
  // thlp2: Variance of thl (overall)                 [K^2]
  // wm_zt: w wind component on thermodynamic levels  [m/s]
  // C11_Skw_fnc: C_11 parameter with Sk_w applied          [-]
  // clubb_params: Array of CLUBB's tunable parameters    [units vary]
  // penta_solve_method: Method to solve then penta-diagonal system
  // fill_holes_type: Specifier for which hole filling method to use
  // l_min_wp2_from_corr_wx: Flag to base the threshold minimum value of wp2 on keeping
  // the overall correlation of w and x (w and rt, as well as w
  // and theta-l) within the limits of -max_mag_correlation_flux
  // to max_mag_correlation_flux.
  // l_use_tke_in_wp2_wp3_K_dfsn: Use TKE in eddy diffusion for wp2 and wp3
  // l_use_wp3_lim_with_smth_Heaviside: Flag to activate mods on wp3 limiters for conv test
  // l_wp2_fill_holes_tke: Turn on additional hole-filling for wp2
  // that takes TKE from up2 and vp2, if necessary
  // l_implemented: True if CLUBB is being implemented and run in a host model
  // up2: u'^2 (momentum levels)                          [m^2/s^2]
  // vp2: v'^2 (momentum levels)                          [m^2/s^2]
  // wp2: w'^2 (momentum levels)                            [m^2/s^2]
  // wp3_zm: w'^3 interpolated to momentum levels              [m^3/s^3]
  // wp3: w'^3 (thermodynamic levels)                       [m^3/s^3]
  // wp2_zt: w'^2 interpolated to thermodyamic levels          [m^2/s^2]
  // err_info: err_info struct containing err_code and err_header

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D rhs_save("wp23_solve::rhs_save", ngrdcol, ((2 * nzm) - 1)); // Saved RHS of band matrix
  Real2D solut("wp23_solve::solut", ngrdcol, ((2 * nzm) - 1)); // Solution to band diagonal system.
  Real2D old_solut(
      "wp23_solve::old_solut", ngrdcol,
      ((2 * nzm) - 1)); // Old solution, used as an initial guess in the bicgstab method
  Real1D rcond("wp23_solve::rcond", ngrdcol); // Est. of the reciprocal of the condition #
  Real2D lhs_wp2_ac_term("wp23_solve::lhs_wp2_ac_term", ngrdcol,
                         nzm); // w'^2 term ac, used for stats
  Real2D lhs_wp2_pr2_term("wp23_solve::lhs_wp2_pr2_term", ngrdcol,
                          nzm); // w'^2 term pr2, used for stats
  Real2D wp2_min_array("wp23_solve::wp2_min_array", ngrdcol,
                       nzm); // Minimum values for wp2 [m^2/s^2]
  Real2D lhs_wp3_ac_term("wp23_solve::lhs_wp3_ac_term", ngrdcol,
                         nzt); // w'^3 term ac, used for stats
  Real2D lhs_wp3_pr2_term("wp23_solve::lhs_wp3_pr2_term", ngrdcol,
                          nzt); // w'^3 term pr2, used for stats
  Real2D C11_Skw_fnc_zeros("wp23_solve::C11_Skw_fnc_zeros", ngrdcol, nzt);
  Real2D C11_Skw_fnc_plus_one("wp23_solve::C11_Skw_fnc_plus_one", ngrdcol, nzt);
  Real2D stats_tmp_zm("wp23_solve::stats_tmp_zm", ngrdcol, nzm);
  Real2D stats_tmp_zt("wp23_solve::stats_tmp_zt", ngrdcol, nzt);
  Real1D C_uu_shr_zeros("wp23_solve::C_uu_shr_zeros", ngrdcol);
  Real1D C_uu_shr_plus_one("wp23_solve::C_uu_shr_plus_one", ngrdcol);
  real_t wp2_avg_before{};
  real_t wp2_avg_after{};
  int k{};
  int i{};
  bool l_clubb_at_least_debug_level_3{};
  Real2D wp2_min_before_limit("wp23_solve::wp2_min_before_limit", ngrdcol, nzm);
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_zm = gr.zm;
  const auto gr_zt = gr.zt;
  const auto gr_dzm = gr.dzm;
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;

  // Save the value of rhs, which will be overwritten with the solution as
  // part of the solving routine.
  {
    const int extent_0 = ((1) > 0 ? (1) <= (((2 * nzm) - 1)) : (1) >= (((2 * nzm) - 1)))
                             ? ((((2 * nzm) - 1)) - (1)) / (1) + 1
                             : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_solve::cells_1", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          rhs_save(i - 1, k - 1) = rhs(i - 1, k - 1);
        });
  }
  if ((penta_solve_method == penta_bicgstab)) {
    {
      const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
            old_solut(i - 1, ((2 * k) - 1) - 1) = wp2(i - 1, k - 1);
          });
    }
    {
      const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_3", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
            old_solut(i - 1, (2 * k) - 1) = wp3(i - 1, k - 1);
          });
    }
  }
  // Matrix solves are bit-different between ascending and descending.
  // This ensures matrices are solved in the same (descending) order,
  // which is useful for ensuring BFBness between grid modes
  if ((l_force_descending_solves && (gr_grid_dir_indx > 0))) {
    // We need to flip in both the vertical dimensions and the bands for the lhs
    {
      // Native whole-array reversal requires owned scratch before overwriting.
      const int matrix_rows = 2 * nzm - 1;
      Real3D reversed("wp23_solve::lhs_reversed", ndiags5, ngrdcol, matrix_rows);
      Kokkos::parallel_for(
          "wp23_solve::reverse_lhs", ndiags5 * ngrdcol * matrix_rows,
          KOKKOS_LAMBDA(const int cell) {
            const int band = cell % ndiags5;
            const int column = (cell / ndiags5) % ngrdcol;
            const int row = cell / (ndiags5 * ngrdcol);
            reversed(band, column, row) = lhs(ndiags5 - 1 - band, column, matrix_rows - 1 - row);
          });
      Kokkos::deep_copy(lhs, reversed);
    }
    {
      // Native whole-array reversal requires owned scratch before overwriting.
      const int matrix_rows = 2 * nzm - 1;
      Real2D reversed("wp23_solve::rhs_reversed", ngrdcol, matrix_rows);
      Kokkos::parallel_for(
          "wp23_solve::reverse_rhs", ngrdcol * matrix_rows, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int row = cell / ngrdcol;
            reversed(column, row) = rhs(column, matrix_rows - 1 - row);
          });
      Kokkos::deep_copy(rhs, reversed);
    }
    if ((penta_solve_method == penta_bicgstab)) {
      {
        // Native whole-array reversal requires owned scratch before overwriting.
        const int matrix_rows = 2 * nzm - 1;
        Real2D reversed("wp23_solve::old_solut_reversed", ngrdcol, matrix_rows);
        Kokkos::parallel_for(
            "wp23_solve::reverse_old_solut", ngrdcol * matrix_rows, KOKKOS_LAMBDA(const int cell) {
              const int column = cell % ngrdcol;
              const int row = cell / ngrdcol;
              reversed(column, row) = old_solut(column, matrix_rows - 1 - row);
            });
        Kokkos::deep_copy(old_solut, reversed);
      }
    }
  }
  // Solve the system with LAPACK
  if ((stats.flags().l_sample && f2k::var_on_stats_list(stats, "wp23_matrix_condt_num"))) {
    // Solve the system and return condition number
    // Note: When using lapack this can change the answer slightly
    {
      // Existing solver boundary: Native scratch stays (band,column,row),
      // while the C++ band_solve owner uses (column,band,row). Restore InOut lhs.
      const int matrix_rows = 2 * nzm - 1;
      Real3D lhs_solve("wp23_solve::lhs_solve", ngrdcol, ndiags5, matrix_rows);
      Kokkos::parallel_for(
          "wp23_solve::pack_lhs", ngrdcol * ndiags5 * matrix_rows, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int band = (cell / ngrdcol) % ndiags5;
            const int row = cell / (ngrdcol * ndiags5);
            lhs_solve(column, band, row) = lhs(band, column, row);
          });
      band_solve("wp2_wp3", penta_solve_method, ngrdcol, 2, 2, ((2 * nzm) - 1), l_implemented,
                 lhs_solve, rhs, err_info, solut, nullptr, &rcond);
      Kokkos::parallel_for(
          "wp23_solve::restore_lhs", ngrdcol * ndiags5 * matrix_rows,
          KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int band = (cell / ngrdcol) % ndiags5;
            const int row = cell / (ngrdcol * ndiags5);
            lhs(band, column, row) = lhs_solve(column, band, row);
          });
    }

    if (stats.flags().l_sample) {
      // Est. of the condition number of the w'^2/w^3 LHS matrix
      {
        // Materialize the Native array expression only at this stats boundary.
        Real1D rcond_stats_0("wp23_solve::rcond_stats_0", ngrdcol);
        Kokkos::parallel_for(
            "wp23_solve::rcond_stats_0", ngrdcol, KOKKOS_LAMBDA(const int cell) {
              const int i = cell % ngrdcol + 1;
              rcond_stats_0(i - 1) = (one / rcond(i - 1));
            });
        f2k::stats_update("wp23_matrix_condt_num", f2k::Real1D(rcond_stats_0), stats);
      }
    }
  } else {
    // Solve the system
    {
      // Existing solver boundary: Native scratch stays (band,column,row),
      // while the C++ band_solve owner uses (column,band,row). Restore InOut lhs.
      const int matrix_rows = 2 * nzm - 1;
      Real3D lhs_solve("wp23_solve::lhs_solve", ngrdcol, ndiags5, matrix_rows);
      Kokkos::parallel_for(
          "wp23_solve::pack_lhs", ngrdcol * ndiags5 * matrix_rows, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int band = (cell / ngrdcol) % ndiags5;
            const int row = cell / (ngrdcol * ndiags5);
            lhs_solve(column, band, row) = lhs(band, column, row);
          });
      band_solve("wp2_wp3", penta_solve_method, ngrdcol, 2, 2, ((2 * nzm) - 1), l_implemented,
                 lhs_solve, rhs, err_info, solut, &old_solut, nullptr);
      Kokkos::parallel_for(
          "wp23_solve::restore_lhs", ngrdcol * ndiags5 * matrix_rows,
          KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int band = (cell / ngrdcol) % ndiags5;
            const int row = cell / (ngrdcol * ndiags5);
            lhs(band, column, row) = lhs_solve(column, band, row);
          });
    }
  }
  // Flip the back to the ascending direction if we forced the solve
  // to be in descending mode
  if ((l_force_descending_solves && (gr_grid_dir_indx > 0))) {
    {
      // Native whole-array reversal requires owned scratch before overwriting.
      const int matrix_rows = 2 * nzm - 1;
      Real2D reversed("wp23_solve::solut_reversed", ngrdcol, matrix_rows);
      Kokkos::parallel_for(
          "wp23_solve::reverse_solut", ngrdcol * matrix_rows, KOKKOS_LAMBDA(const int cell) {
            const int column = cell % ngrdcol;
            const int row = cell / ngrdcol;
            reversed(column, row) = solut(column, matrix_rows - 1 - row);
          });
      Kokkos::deep_copy(solut, reversed);
    }
  }
  if (f2k::clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << err_info.err_header_global << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Error in wp23_solve calling band_solve for wp2_wp3" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp2 & wp3 LU decomp. failed" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp2 and wp3 LHS" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      for (k = 1; k <= ((nzm - 1)); k += 1) {
        for (i = 1; i <= (ngrdcol); i += 1) {
          {
            std::ostringstream diagnostic;
            diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            diagnostic << "grid col = " << ' ';
            diagnostic << i << ' ';
            diagnostic << "zm level = " << ' ';
            diagnostic << k << ' ';
            diagnostic << "height [m] = " << ' ';
            diagnostic << to_host_space(gr_zm)(i - 1, k - 1) << ' ';
            diagnostic << "LHS = " << ' ';
            {
              const auto field_host = to_host_space(
                  Kokkos::subview(lhs, std::make_pair(1 - 1, 5), i - 1, ((2 * k) - 1) - 1));
              for (std::size_t element = 0; element < field_host.size(); ++element)
                diagnostic << field_host.data()[element] << ' ';
            }
            std::cerr << diagnostic.str() << '\n';
          }
          {
            std::ostringstream diagnostic;
            diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            diagnostic << "grid col = " << ' ';
            diagnostic << i << ' ';
            diagnostic << "zt level = " << ' ';
            diagnostic << k << ' ';
            diagnostic << "height [m] = " << ' ';
            diagnostic << to_host_space(gr_zt)(i - 1, k - 1) << ' ';
            diagnostic << "LHS = " << ' ';
            {
              const auto field_host =
                  to_host_space(Kokkos::subview(lhs, std::make_pair(1 - 1, 5), i - 1, (2 * k) - 1));
              for (std::size_t element = 0; element < field_host.size(); ++element)
                diagnostic << field_host.data()[element] << ' ';
            }
            std::cerr << diagnostic.str() << '\n';
          }
        }
      }
      // k = 1, nzm-1
      for (i = 1; i <= (ngrdcol); i += 1) {
        k = nzm;
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "grid col = " << ' ';
          diagnostic << i << ' ';
          diagnostic << "zm level = " << ' ';
          diagnostic << k << ' ';
          diagnostic << "height [m] = " << ' ';
          diagnostic << to_host_space(gr_zm)(i - 1, k - 1) << ' ';
          diagnostic << "LHS = " << ' ';
          {
            const auto field_host = to_host_space(
                Kokkos::subview(lhs, std::make_pair(1 - 1, 5), i - 1, ((2 * k) - 1) - 1));
            for (std::size_t element = 0; element < field_host.size(); ++element)
              diagnostic << field_host.data()[element] << ' ';
          }
          std::cerr << diagnostic.str() << '\n';
        }
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "wp2 and wp3 RHS" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      for (k = 1; k <= ((nzm - 1)); k += 1) {
        for (i = 1; i <= (ngrdcol); i += 1) {
          {
            std::ostringstream diagnostic;
            diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            diagnostic << "grid col = " << ' ';
            diagnostic << i << ' ';
            diagnostic << "zm level = " << ' ';
            diagnostic << k << ' ';
            diagnostic << "height [m] = " << ' ';
            diagnostic << to_host_space(gr_zm)(i - 1, k - 1) << ' ';
            diagnostic << "RHS = " << ' ';
            diagnostic << to_host_space(rhs_save)(i - 1, ((2 * k) - 1) - 1) << ' ';
            std::cerr << diagnostic.str() << '\n';
          }
          {
            std::ostringstream diagnostic;
            diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            diagnostic << "grid col = " << ' ';
            diagnostic << i << ' ';
            diagnostic << "zt level = " << ' ';
            diagnostic << k << ' ';
            diagnostic << "height [m] = " << ' ';
            diagnostic << to_host_space(gr_zt)(i - 1, k - 1) << ' ';
            diagnostic << "RHS = " << ' ';
            diagnostic << to_host_space(rhs_save)(i - 1, (2 * k) - 1) << ' ';
            std::cerr << diagnostic.str() << '\n';
          }
        }
      }
      // k = 1, nzm-1
      for (i = 1; i <= (ngrdcol); i += 1) {
        k = nzm;
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "grid col = " << ' ';
          diagnostic << i << ' ';
          diagnostic << "zm level = " << ' ';
          diagnostic << k << ' ';
          diagnostic << "height [m] = " << ' ';
          diagnostic << to_host_space(gr_zm)(i - 1, k - 1) << ' ';
          diagnostic << "RHS = " << ' ';
          diagnostic << to_host_space(rhs_save)(i - 1, ((2 * k) - 1) - 1) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
      }
      return;
    }
  }
  // Copy result into output arrays and clip
  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_solve::cells_4", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          int k_wp2{};
          k_wp2 = ((2 * k) - 1);
          wp2(i - 1, k - 1) = solut(i - 1, k_wp2 - 1);
        });
  }
  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_solve::cells_5", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          int k_wp3{};
          k_wp3 = (2 * k);
          wp3(i - 1, k - 1) = solut(i - 1, k_wp3 - 1);
        });
  }
  if (stats.flags().l_sample) {
    Kokkos::deep_copy(C_uu_shr_zeros, zero);
    {
      Kokkos::parallel_for(
          "wp23_solve::array_1", ngrdcol, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            C_uu_shr_plus_one(i - 1) = (clubb_params(i - 1, iC_uu_shr) + one);
          });
    }
    Kokkos::deep_copy(C11_Skw_fnc_zeros, zero);
    {
      Kokkos::parallel_for(
          "wp23_solve::array_2", ngrdcol * nzt, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            C11_Skw_fnc_plus_one(i - 1, k - 1) = (C11_Skw_fnc(i - 1, k - 1) + one);
          });
    }
    // Note:  To find the contribution of w'^2 term ac, substitute 0 for the
    // C_uu_shr input to function wp2_terms_ac_pr2_lhs.
    wp2_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr, // In
                         C_uu_shr_zeros, wm_zt, // In
                         lhs_wp2_ac_term);      // Out

    // Note:  To find the contribution of w'^2 term pr2, add 1 to the
    // C_uu_shr input to function wp2_terms_ac_pr2_lhs.
    wp2_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr,    // In
                         C_uu_shr_plus_one, wm_zt, // In
                         lhs_wp2_pr2_term);        // Out

    // Note:  To find the contribution of w'^3 term ac, substitute 0 for the
    // C_ll skewness function input to function wp3_terms_ac_pr2_lhs.
    wp3_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr, C11_Skw_fnc_zeros, wm_zm, // In
                         lhs_wp3_ac_term);                                // Out

    // Note:  To find the contribution of w'^3 term pr2, add 1 to the
    // C_ll skewness function input to function wp3_terms_ac_pr2_lhs.
    wp3_terms_ac_pr2_lhs(nzm, nzt, ngrdcol, gr,       // In
                         C11_Skw_fnc_plus_one, wm_zm, // In
                         lhs_wp3_pr2_term);           // Out
  }
  if (stats.flags().l_sample) {
    // Finalize implicit contributions for wp2
    // Finalize the implicit pieces of wp2/wp3 budgets after the solve.
    // The companion explicit pieces were opened earlier with begin_budget.
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_6", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) =
                ((-(gamma_over_implicit_ts)*lhs_dp1_wp2(i - 1, k - 1)) * wp2(i - 1, k - 1));
          });
    }
    // w'^2 term dp1 has both implicit and explicit components;
    // call stat_end_update_pt.
    // w'^2 term dp1 has both implicit and explicit components.
    f2k::stats_finalize_budget("wp2_dp1", stats_tmp_zm, stats);
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_7", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            int km1{};
            int kp1{};
            km1 = ((k - gr_grid_dir_indx) > 1 ? (k - gr_grid_dir_indx) : 1);
            kp1 = ((k + gr_grid_dir_indx) < nzm ? (k + gr_grid_dir_indx) : nzm);
            stats_tmp_zm(i - 1, k - 1) =
                (((-(lhs_diff_zm((2 + gr_grid_dir_indx) - 1, i - 1, k - 1)) * wp2(i - 1, km1 - 1)) -
                  (lhs_diff_zm(2 - 1, i - 1, k - 1) * wp2(i - 1, k - 1))) -
                 (lhs_diff_zm((2 - gr_grid_dir_indx) - 1, i - 1, k - 1) * wp2(i - 1, kp1 - 1)));
          });
    }
    if ((l_crank_nich_diff || l_use_tke_in_wp2_wp3_K_dfsn)) {
      // dp2 has explicit and implicit pieces for CN/tke diffusion.
      // w'^2 term dp2 has both implicit and explicit components (if the
      // Crank-Nicholson scheme is selected or if l_use_tke_in_wp2_wp3_K_dfsn is true);
      // w'^2 term dp2 has both implicit and explicit components.
      f2k::stats_finalize_budget("wp2_dp2", stats_tmp_zm, stats);
    } else {
      // Otherwise dp2 is fully implicit.
      // If neither of these flags is true, then w'^2 term dp2 is
      // completely implicit; call stat_update_var_pt.
      // w'^2 term dp2 is completely implicit; call stat_update_var_pt.
      f2k::stats_update("wp2_dp2", stats_tmp_zm, stats);
    }
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_8", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) =
                ((-(lhs_ta_wp2(2 - 1, i - 1, k - 1)) * wp3(i - 1, (k - 1) - 1)) +
                 (-(lhs_ta_wp2(1 - 1, i - 1, k - 1)) * wp3(i - 1, k - 1)));
          });
    }
    // w'^2 term ta is completely implicit; call stat_update_var_pt.
    f2k::stats_update("wp2_ta", stats_tmp_zm, stats);
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_9", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            int km1{};
            int kp1{};
            km1 = ((k - gr_grid_dir_indx) > 1 ? (k - gr_grid_dir_indx) : 1);
            kp1 = ((k + gr_grid_dir_indx) < nzm ? (k + gr_grid_dir_indx) : nzm);
            stats_tmp_zm(i - 1, k - 1) =
                (((-(lhs_ma_zm((2 + gr_grid_dir_indx) - 1, i - 1, k - 1)) * wp2(i - 1, km1 - 1)) -
                  (lhs_ma_zm(2 - 1, i - 1, k - 1) * wp2(i - 1, k - 1))) -
                 (lhs_ma_zm((2 - gr_grid_dir_indx) - 1, i - 1, k - 1) * wp2(i - 1, kp1 - 1)));
          });
    }
    // w'^2 term ma is completely implicit; call stat_update_var_pt.
    f2k::stats_update("wp2_ma", stats_tmp_zm, stats);
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_10", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) = (-(lhs_wp2_ac_term(i - 1, k - 1)) * wp2(i - 1, k - 1));
          });
    }
    // w'^2 term ac is completely implicit; call stat_update_var_pt.
    f2k::stats_update("wp2_ac", stats_tmp_zm, stats);
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_11", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) =
                ((-(gamma_over_implicit_ts)*lhs_pr1_wp2(i - 1, k - 1)) * wp2(i - 1, k - 1));
          });
    }
    // w'^2 term pr1 has both implicit and explicit components;
    // call stat_end_update_pt.
    // w'^2 term pr1 has both implicit and explicit components.
    f2k::stats_finalize_budget("wp2_pr1", stats_tmp_zm, stats);
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_12", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) = (-(lhs_wp2_pr2_term(i - 1, k - 1)) * wp2(i - 1, k - 1));
          });
    }
    // w'^2 term pr2 has both implicit and explicit components;
    // call stat_end_update_pt.
    // w'^2 term pr2 has both implicit and explicit components.
    f2k::stats_finalize_budget("wp2_pr2", stats_tmp_zm, stats);
    // Finalize implicit contributions for wp3
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_13", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) =
                ((-(gamma_over_implicit_ts)*lhs_pr1_wp3(i - 1, k - 1)) * wp3(i - 1, k - 1));
          });
    }
    // w'^3 term pr1 has both implicit and explicit components;
    // call stat_end_update_pt
    // w'^3 term pr1 has both implicit and explicit components;
    // term more numerically stable (see note above for LHS turbulent
    // w'^3 term pr1 has both implicit and explicit components.
    f2k::stats_finalize_budget("wp3_pr1", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_14", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            int km1{};
            int kp1{};
            km1 = ((k - gr_grid_dir_indx) > 1 ? (k - gr_grid_dir_indx) : 1);
            kp1 = ((k + gr_grid_dir_indx) < nzt ? (k + gr_grid_dir_indx) : nzt);
            stats_tmp_zt(i - 1, k - 1) =
                (((-(lhs_diff_zt((2 + gr_grid_dir_indx) - 1, i - 1, k - 1)) * wp3(i - 1, km1 - 1)) -
                  (lhs_diff_zt(2 - 1, i - 1, k - 1) * wp3(i - 1, k - 1))) -
                 (lhs_diff_zt((2 - gr_grid_dir_indx) - 1, i - 1, k - 1) * wp3(i - 1, kp1 - 1)));
          });
    }
    if ((l_crank_nich_diff || l_use_tke_in_wp2_wp3_K_dfsn)) {
      // dp1 has explicit and implicit pieces for CN/tke diffusion.
      // w'^3 term dp1 has both implicit and explicit components (if the
      // Crank-Nicholson scheme is selected or l_use_tke_in_wp2_wp3_K_dfsn is true);
      // w'^3 term dp1 has both implicit and explicit components.
      f2k::stats_finalize_budget("wp3_dp1", stats_tmp_zt, stats);
    } else {
      // Otherwise dp1 is fully implicit.
      // If neither of these flags is true, then w'^3 term dp1 is
      // w'^3 term dp1 is completely implicit; call stat_update_var_pt.
      f2k::stats_update("wp3_dp1", stats_tmp_zt, stats);
    }
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_15", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) =
                ((((((-(gamma_over_implicit_ts)*wp3_term_ta_lhs_result(5 - 1, i - 1, k - 1)) *
                     wp3(i - 1, (k - 1) - 1)) -
                    ((gamma_over_implicit_ts * wp3_term_ta_lhs_result(4 - 1, i - 1, k - 1)) *
                     wp2(i - 1, k - 1))) -
                   ((gamma_over_implicit_ts * wp3_term_ta_lhs_result(3 - 1, i - 1, k - 1)) *
                    wp3(i - 1, k - 1))) -
                  ((gamma_over_implicit_ts * wp3_term_ta_lhs_result(2 - 1, i - 1, k - 1)) *
                   wp2(i - 1, (k + 1) - 1))) -
                 ((gamma_over_implicit_ts * wp3_term_ta_lhs_result(1 - 1, i - 1, k - 1)) *
                  wp3(i - 1, (k + 1) - 1)));
          });
    }
    // w'^3 term ta has both implicit and explicit components;
    // make the term more numerically stable (see note above for
    // RHS turbulent advection (ta) term).
    // Call stat_begin_update_pt.  Since stat_begin_update_pt
    // automatically subtracts the value sent in, reverse the sign
    // on the input value.
    // w'^3 term ta has both implicit and explicit components.
    f2k::stats_finalize_budget("wp3_ta", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_16", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) =
                (((-(gamma_over_implicit_ts)*lhs_adv_tp_wp3(2 - 1, i - 1, k - 1)) *
                  wp2(i - 1, k - 1)) -
                 ((gamma_over_implicit_ts * lhs_adv_tp_wp3(1 - 1, i - 1, k - 1)) *
                  wp2(i - 1, (k + 1) - 1)));
          });
    }
    // term more numerically stable (see note above for RHS turbulent
    // production (tp) term).  Call stat_begin_update_pt.  Since
    // w'^3 term tp has both implicit and explicit components;
    // w'^3 term tp has both implicit and explicit components.
    f2k::stats_finalize_budget("wp3_tp", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_17", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) =
                (((-(gamma_over_implicit_ts)*lhs_pr_tp_wp3(2 - 1, i - 1, k - 1)) *
                  wp2(i - 1, k - 1)) -
                 ((gamma_over_implicit_ts * lhs_pr_tp_wp3(1 - 1, i - 1, k - 1)) *
                  wp2(i - 1, (k + 1) - 1)));
          });
    }
    // w'^3 term pr_tp same as above tp term but opposite sign.
    f2k::stats_finalize_budget("wp3_pr_tp", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_18", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) =
                (((((-(wp3_pr3_lhs(5 - 1, i - 1, k - 1)) * wp3(i - 1, (k - 1) - 1)) -
                    (wp3_pr3_lhs(4 - 1, i - 1, k - 1) * wp2(i - 1, k - 1))) -
                   (wp3_pr3_lhs(3 - 1, i - 1, k - 1) * wp3(i - 1, k - 1))) -
                  (wp3_pr3_lhs(2 - 1, i - 1, k - 1) * wp2(i - 1, (k + 1) - 1))) -
                 (wp3_pr3_lhs(1 - 1, i - 1, k - 1) * wp3(i - 1, (k + 1) - 1)));
          });
    }
    // w'^3 pressure term 3 (pr3) explicit (rhs) contribution
    // w'^3 pressure term 3 (pr3) has both implicit and explicit components;
    // w'^3 pressure term 3 (pr3) has both implicit and explicit components.
    f2k::stats_finalize_budget("wp3_pr3", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_19", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            int km1{};
            int kp1{};
            km1 = ((k - gr_grid_dir_indx) > 1 ? (k - gr_grid_dir_indx) : 1);
            kp1 = ((k + gr_grid_dir_indx) < nzt ? (k + gr_grid_dir_indx) : nzt);
            stats_tmp_zt(i - 1, k - 1) =
                (((-(lhs_ma_zt((2 + gr_grid_dir_indx) - 1, i - 1, k - 1)) * wp3(i - 1, km1 - 1)) -
                  (lhs_ma_zt(2 - 1, i - 1, k - 1) * wp3(i - 1, k - 1))) -
                 (lhs_ma_zt((2 - gr_grid_dir_indx) - 1, i - 1, k - 1) * wp3(i - 1, kp1 - 1)));
          });
    }
    // w'^3 term ma is completely implicit; call stat_update_var_pt.
    f2k::stats_update("wp3_ma", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_20", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) = (-(lhs_wp3_ac_term(i - 1, k - 1)) * wp3(i - 1, k - 1));
          });
    }
    // w'^3 term ac is completely implicit; call stat_update_var_pt.
    f2k::stats_update("wp3_ac", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_21", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) = (-(lhs_wp3_pr2_term(i - 1, k - 1)) * wp3(i - 1, k - 1));
          });
    }
    // w'^3 term pr2 has both implicit and explicit components;
    // call stat_end_update_pt
    // w'^3 term pr2 has both implicit and explicit components.
    f2k::stats_finalize_budget("wp3_pr2", stats_tmp_zt, stats);
  }
  if (stats.flags().l_sample) {
    // Store previous value for effect of the positive definite scheme
    // Bracket the positive-definite/hole-filling adjustments so the pd
    // budget terms report only the limiter-induced change.
    {
      // Materialize the Native array expression only at this stats boundary.
      Real2D wp2_stats_3("wp23_solve::wp2_stats_3", ngrdcol, nzm);
      Kokkos::parallel_for(
          "wp23_solve::wp2_stats_3", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            wp2_stats_3(i - 1, k - 1) = (wp2(i - 1, k - 1) / dt);
          });
      f2k::stats_begin_budget("wp2_pd", wp2_stats_3, stats);
    }
    // up2_pd and vp2_pd are also modified in a call to pos_definite_variances
    // in advance_xp2_xpyp, so we need to use stat_modify here instead of stat_begin_update
    {
      // Materialize the Native array expression only at this stats boundary.
      Real2D up2_stats_4("wp23_solve::up2_stats_4", ngrdcol, nzm);
      Kokkos::parallel_for(
          "wp23_solve::up2_stats_4", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            up2_stats_4(i - 1, k - 1) = (up2(i - 1, k - 1) / dt);
          });
      f2k::stats_begin_budget("up2_pd", up2_stats_4, stats);
    }
    {
      // Materialize the Native array expression only at this stats boundary.
      Real2D vp2_stats_5("wp23_solve::vp2_stats_5", ngrdcol, nzm);
      Kokkos::parallel_for(
          "wp23_solve::vp2_stats_5", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            vp2_stats_5(i - 1, k - 1) = (vp2(i - 1, k - 1) / dt);
          });
      f2k::stats_begin_budget("vp2_pd", vp2_stats_5, stats);
    }
  }
  if ((fill_holes_type != 0)) {
    if (stats.flags().l_sample) {
      f2k::stats_update("wp2_before_hf", wp2, stats);
    }
    // Debugging to check if wp2 hole-filling is conservative:
    // Compute vertical_avg before and after hole-filling
    // Output warning in case the average changed from hole-filling
    if (f2k::clubb_at_least_debug_level_api(3)) {
      {
        // Native passes bulk rank2 arrays to an explicit rank1 dummy. Pack
        // its first nzm elements in Fortran column-fastest sequence association.
        Real1D density_window("wp23_solve::density_window", nzm);
        Real1D field_window("wp23_solve::field_window", nzm);
        Real1D thickness_window("wp23_solve::thickness_window", nzm);
        Kokkos::parallel_for(
            "wp23_solve::average_window", nzm, KOKKOS_LAMBDA(const int index) {
              const int column = index % ngrdcol;
              const int level = index / ngrdcol;
              density_window(index) = rho_ds_zm(column, level);
              field_window(index) = wp2(column, level);
              thickness_window(index) = gr_dzm(column, level);
            });
        wp2_avg_before = vertical_avg(nzm, density_window, field_window, thickness_window);
      }
    }
    // Use a simple hole filling algorithm
    // upper_hf_level = nzm-1 since we are filling the zm levels
    f2k::fill_holes_vertical_api(nzm, ngrdcol, w_tol_sqd, (gr_k_lb_zm + gr_grid_dir_indx),
                                 (gr_k_ub_zm - gr_grid_dir_indx), gr_dzm, rho_ds_zm,
                                 gr_grid_dir_indx, fill_holes_type, wp2);
    if (stats.flags().l_sample) {
      f2k::stats_update("wp2_after_hf", wp2, stats);
    }
    if (f2k::clubb_at_least_debug_level_api(3)) {
      {
        // Native passes bulk rank2 arrays to an explicit rank1 dummy. Pack
        // its first nzm elements in Fortran column-fastest sequence association.
        Real1D density_window("wp23_solve::density_window", nzm);
        Real1D field_window("wp23_solve::field_window", nzm);
        Real1D thickness_window("wp23_solve::thickness_window", nzm);
        Kokkos::parallel_for(
            "wp23_solve::average_window", nzm, KOKKOS_LAMBDA(const int index) {
              const int column = index % ngrdcol;
              const int level = index / ngrdcol;
              density_window(index) = rho_ds_zm(column, level);
              field_window(index) = wp2(column, level);
              thickness_window(index) = gr_dzm(column, level);
            });
        wp2_avg_after = vertical_avg(nzm, density_window, field_window, thickness_window);
      }
      if ((Kokkos::abs((wp2_avg_before - wp2_avg_after)) >
           (std::numeric_limits<real_t>::epsilon() * 1000))) {
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "Warning! Hole-filling on wp2 is not conservative! " << ' ';
          diagnostic << "The difference is " << ' ';
          diagnostic << (wp2_avg_after - wp2_avg_before) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
      }
    }
    if (l_wp2_fill_holes_tke) {
      f2k::fill_holes_wp2_from_horz_tke(nzm, ngrdcol, w_tol_sqd, 2, (nzm - 1), wp2, up2, vp2);
    }
  }
  // wp2 hole filling check
  if (stats.flags().l_sample) {
    // Store updated value for effect of the positive definite scheme
    {
      // Materialize the Native array expression only at this stats boundary.
      Real2D wp2_stats_6("wp23_solve::wp2_stats_6", ngrdcol, nzm);
      Kokkos::parallel_for(
          "wp23_solve::wp2_stats_6", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            wp2_stats_6(i - 1, k - 1) = (wp2(i - 1, k - 1) / dt);
          });
      f2k::stats_finalize_budget("wp2_pd", wp2_stats_6, stats);
    }
    {
      // Materialize the Native array expression only at this stats boundary.
      Real2D up2_stats_7("wp23_solve::up2_stats_7", ngrdcol, nzm);
      Kokkos::parallel_for(
          "wp23_solve::up2_stats_7", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            up2_stats_7(i - 1, k - 1) = (up2(i - 1, k - 1) / dt);
          });
      f2k::stats_finalize_budget("up2_pd", up2_stats_7, stats, false);
    }
    {
      // Materialize the Native array expression only at this stats boundary.
      Real2D vp2_stats_8("wp23_solve::vp2_stats_8", ngrdcol, nzm);
      Kokkos::parallel_for(
          "wp23_solve::vp2_stats_8", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            vp2_stats_8(i - 1, k - 1) = (vp2(i - 1, k - 1) / dt);
          });
      f2k::stats_finalize_budget("vp2_pd", vp2_stats_8, stats, false);
    }
  }
  // Clip <w'^2> at a minimum and maximum threshold.
  // We attempt to clip extreme values of wp2 to prevent a crash of the
  // type found on the Climate Process Team ticket #49.  Chris Golaz found that
  // instability caused by large wp2 in CLUBB led unrealistic results in AM3.
  // -dschanen 11 Apr 2011
  l_clubb_at_least_debug_level_3 = f2k::clubb_at_least_debug_level_api(3);
  // The value of <w'^2> is not allowed to become smaller than the threshold
  // value of w_tol^2.  Additionally, that threshold value may be boosted at
  // any grid level in order to keep the overall correlation of w and rt or
  // the overall correlation of w and theta-l between the values of
  // -max_mag_correlation_flux and max_mag_correlation_flux by boosting <w'^2>
  // rather than by limiting the magnitude of <w'rt'> or <w'thl'>.
  if (l_min_wp2_from_corr_wx) {
    // The overall correlation of w and rt is:
    // corr_w_rt = wprtp / ( sqrt( wp2 ) * sqrt( rtp2 ) );
    // and the overall correlation of w and thl is:
    // corr_w_thl = wpthlp / ( sqrt( wp2 ) * sqrt( thlp2 ) ).
    // Squaring both sides, the equations becomes:
    // corr_w_rt^2 = wprtp^2 / ( wp2 * rtp2 ); and
    // corr_w_thl^2 = wpthlp^2 / ( wp2 * thlp2 ).
    // Using max_mag_correlation_flux for the correlation and then solving for
    // the minimum of wp2, the equation becomes:
    // wp2|_min = max( wprtp^2 / ( rtp2 * max_mag_correlation_flux^2 ),
    // wpthlp^2 / ( thlp2 * max_mag_correlation_flux^2 ) ).
    // And similarly for upwp and vpwp
    {
      const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_22", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
            // Lower Native MAX in argument order, evaluating each velocity-variance
            // candidate once [m^2/s^2]; use distinct const values, not reused scratch.
            const real_t min_wp2_from_rt =
                ((wprtp(i - 1, k - 1) * wprtp(i - 1, k - 1)) /
                 (rtp2(i - 1, k - 1) * (max_mag_correlation_flux * max_mag_correlation_flux)));
            const real_t min_wp2_from_thl =
                ((wpthlp(i - 1, k - 1) * wpthlp(i - 1, k - 1)) /
                 (thlp2(i - 1, k - 1) * (max_mag_correlation_flux * max_mag_correlation_flux)));
            const real_t min_wp2_from_u =
                ((upwp(i - 1, k - 1) * upwp(i - 1, k - 1)) /
                 (up2(i - 1, k - 1) * (max_mag_correlation_flux * max_mag_correlation_flux)));
            const real_t min_wp2_from_v =
                ((vpwp(i - 1, k - 1) * vpwp(i - 1, k - 1)) /
                 (vp2(i - 1, k - 1) * (max_mag_correlation_flux * max_mag_correlation_flux)));
            const real_t min_wp2_after_1 =
                w_tol_sqd > min_wp2_from_rt ? w_tol_sqd : min_wp2_from_rt;
            const real_t min_wp2_after_2 =
                min_wp2_after_1 > min_wp2_from_thl ? min_wp2_after_1 : min_wp2_from_thl;
            const real_t min_wp2_after_3 =
                min_wp2_after_2 > min_wp2_from_u ? min_wp2_after_2 : min_wp2_from_u;
            wp2_min_array(i - 1, k - 1) =
                min_wp2_after_3 > min_wp2_from_v ? min_wp2_after_3 : min_wp2_from_v;
            wp2_min_before_limit(i - 1, k - 1) = wp2_min_array(i - 1, k - 1);
          });
    }
    // Native warnings describe the pre-limit values. Stage only for debug3,
    // print in source level/column order, then apply the unchanged one ceiling.
    if (l_clubb_at_least_debug_level_3) {
      const auto minimum_host = to_host_space(wp2_min_before_limit);
      const auto wprtp_host = to_host_space(wprtp);
      const auto rtp2_host = to_host_space(rtp2);
      const auto wpthlp_host = to_host_space(wpthlp);
      const auto thlp2_host = to_host_space(thlp2);
      const auto upwp_host = to_host_space(upwp);
      const auto up2_host = to_host_space(up2);
      const auto vpwp_host = to_host_space(vpwp);
      const auto vp2_host = to_host_space(vp2);
      for (int level = 0; level < nzm; ++level)
        for (int column = 0; column < ngrdcol; ++column) {
          if (minimum_host(column, level) > one) {
            std::ostringstream diagnostic;
            diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            diagnostic << "Warning: wp2_min_array(" << column + 1 << "," << level + 1
                       << ") = " << minimum_host(column, level)
                       << "> 1.0. The threshold value is limited to 1.0." << '\n';
            diagnostic << "Checking inputs:" << '\n';
            diagnostic << "wprtp = " << wprtp_host(column, level)
                       << ", rtp2 = " << rtp2_host(column, level) << '\n';
            diagnostic << "wpthlp = " << wpthlp_host(column, level)
                       << ", thlp2 = " << thlp2_host(column, level) << '\n';
            diagnostic << "upwp = " << upwp_host(column, level)
                       << ", up2 = " << up2_host(column, level) << '\n';
            diagnostic << "vpwp = " << vpwp_host(column, level)
                       << ", vp2 = " << vp2_host(column, level) << '\n';
            std::cerr << diagnostic.str();
          }
        }
    }
    Kokkos::parallel_for(
        "wp23_solve::limit_wp2_min", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
          const int column = cell % ngrdcol;
          const int level = cell / ngrdcol;
          wp2_min_array(column, level) =
              one < wp2_min_array(column, level) ? one : wp2_min_array(column, level);
        });
  } else {
    // Consider only the minimum tolerance threshold value for wp2.
    {
      const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_solve::cells_23", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
            wp2_min_array(i - 1, k - 1) = w_tol_sqd;
          });
    }
  }
  // l_min_wp2_from_corr_wx
  f2k::clip_variance(nzm, ngrdcol, gr, clip_wp2, dt, wp2_min_array, // In
                     stats,                                         // In
                     wp2,                                           // InOut
                     wp2_max);                                      // optional in

  // Interpolate w'^2 from momentum levels to thermodynamic levels.
  // This is used for the clipping of w'^3 according to the value
  // of Sk_w now that w'^2 and w'^3 have been advanced one timestep.
  Kokkos::deep_copy(wp2_zt, zm2zt_api(nzm, nzt, ngrdcol, gr, wp2, w_tol_sqd));
  // Positive definite quantity
  // Clip w'^3 by limiting skewness.
  f2k::clip_skewness(nzt, ngrdcol, gr, dt, sfc_elevation,                              // In
                     Kokkos::subview(clubb_params, Kokkos::ALL, iSkw_max_mag), wp2_zt, // In
                     l_use_wp3_lim_with_smth_Heaviside,                                // In
                     stats,                                                            // In
                     wp3);                                                             // InOut
  // Compute wp3_zm for output purposes
  Kokkos::deep_copy(wp3_zm, zt2zm_api(nzm, nzt, ngrdcol, gr, wp3));
  return;
}

// =============================================================================
void wp23_lhs(int nzm, int nzt, int ngrdcol, real_t dt,                         // In
              Real3D wp3_term_ta_lhs_result,                                    // In
              Real3D lhs_diff_zm, Real3D lhs_diff_zt, Real3D lhs_ma_zm,         // In
              Real3D lhs_ma_zt, Real3D lhs_ta_wp2,                              // In
              Real3D lhs_tp_wp3,                                                // In
              Real2D lhs_ac_pr2_wp2, Real2D lhs_ac_pr2_wp3, Real2D lhs_dp1_wp2, // In
              Real2D lhs_pr1_wp3, Real2D lhs_pr1_wp2, Real2D lhs_splat_wp2,     // In
              Real2D lhs_splat_wp3,                                             // In
              Real3D lhs) {                                                     // Out
  //
  // Description:
  // Compute LHS band diagonal matrix for w'^2 and w'^3.
  // This subroutine computes the implicit portion
  // of the w'^2 and w'^3 equations.
  //
  // Boundary conditions
  //
  //   Both wp2 and wp3 used fixed-point boundary conditions.
  //   Therefore, anything set in the above loop at both the upper
  //   and lower boundaries would be overwritten here.  However, the
  //   above loop does not extend to the boundary levels.  An array
  //   with a value of 1 at the main diagonal on the left-hand side
  //   and with values of 0 at all other diagonals on the left-hand
  //   side will preserve the right-hand side value at that level.
  //
  //      wp2(1)  wp3(1)  ...  wp3(nzt)  wp2(nzm)
  //     [  0.0     0.0          0.0       0.0  ]
  //     [  0.0     0.0          0.0       0.0  ]
  //     [  1.0     1.0   ...    1.0       1.0  ]
  //     [  0.0     0.0          0.0       0.0  ]
  //     [  0.0     0.0          0.0       0.0  ]
  //
  //
  //  WARNING: This subroutine has been optimized. Significant changes could
  //           noticeably  impact computational efficiency. See clubb:ticket:834
  //-------------------------------------------------------------------------------

  // ----------------------- Input Variables -----------------------
  // dt: Timestep length                            [s]
  // lhs_diff_zm: Completely implicit diffusion term for w'2
  // lhs_ma_zm: Mean advection term for w'2
  // lhs_diff_zt: Completely implicit diffusion term for w'3
  // lhs_ma_zt: Mean advection term for w'3
  // lhs_ta_wp2: Turbulent advection terms for wp2
  // lhs_tp_wp3: Turbulent production terms of w'^3
  // lhs_ac_pr2_wp2: Accumulation terms of w'^2 and w'^2 pressure term 2
  // lhs_dp1_wp2: Dissipation terms 1 for w'^2
  // lhs_pr1_wp2: Pressure term 1 for w'2
  // lhs_splat_wp2: LHS coefficient of wp2 splatting term  [1/s]
  // lhs_ac_pr2_wp3: Accumulation terms of w'^3 and w'^3 pressure term 2
  // lhs_pr1_wp3: Dissipation terms 1 for w'^3
  // lhs_splat_wp3: LHS coefficient of wp3 splatting term  [1/s]
  // ----------------------- Output Variable -----------------------
  // lhs: Implicit contributions to wp2/wp3 (band diag. matrix)

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  real_t invrs_dt{}; // Inverse of dt, 1/dt, used for computational efficiency

  // Calculate invrs_dt
  invrs_dt = (real_t(1.0) / dt);
  {
    const int extent_0 = ((1) > 0 ? (1) <= (((2 * nzm) - 1)) : (1) >= (((2 * nzm) - 1)))
                             ? ((((2 * nzm) - 1)) - (1)) / (1) + 1
                             : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    const int extent_2 =
        ((1) > 0 ? (1) <= (ndiags5) : (1) >= (ndiags5)) ? ((ndiags5) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_lhs::cells_1", extent_0 * extent_1 * extent_2, KOKKOS_LAMBDA(const int cell) {
          const int b = (1) + ((cell / (1)) % extent_2) * (1);
          const int i = (1) + ((cell / (1 * extent_2)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_2 * extent_1)) % extent_0) * (1);
          lhs(b - 1, i - 1, k - 1) = real_t(0.0);
        });
  }
  // Lower (upper) boundary for w'2 on ascending (descending) grid.
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ndiags5) : (1) >= (ndiags5)) ? ((ndiags5) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_lhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int b = (1) + ((cell / (1)) % extent_1) * (1);
          const int i = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          if ((b != 3)) {
            lhs(b - 1, i - 1, 1 - 1) = real_t(0.0);
          } else {
            lhs(b - 1, i - 1, 1 - 1) = real_t(1.0);
          }
        });
  }
  // Lower (upper) boundary for w'3 on ascending (descending) grid.
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ndiags5) : (1) >= (ndiags5)) ? ((ndiags5) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_lhs::cells_3", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int b = (1) + ((cell / (1)) % extent_1) * (1);
          const int i = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          if ((b != 3)) {
            lhs(b - 1, i - 1, 2 - 1) = real_t(0.0);
          } else {
            lhs(b - 1, i - 1, 2 - 1) = real_t(1.0);
          }
        });
  }
  // Combine terms to calculate non-boundary lhs values
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_lhs::cells_4", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          int k_wp2{};
          k_wp2 = ((2 * k) - 1);
          // ------ w'2 ------
          // LHS mean advection (ma) and diffusion (diff) terms
          lhs(1 - 1, i - 1, k_wp2 - 1) =
              ((lhs(1 - 1, i - 1, k_wp2 - 1) + lhs_ma_zm(1 - 1, i - 1, k - 1)) +
               lhs_diff_zm(1 - 1, i - 1, k - 1));
          // LHS turbulent advection (ta) term.
          lhs(2 - 1, i - 1, k_wp2 - 1) =
              (lhs(2 - 1, i - 1, k_wp2 - 1) + lhs_ta_wp2(1 - 1, i - 1, k - 1));
          // LHS mean advection (ma) and diffusion (diff) terms
          lhs(3 - 1, i - 1, k_wp2 - 1) =
              ((lhs(3 - 1, i - 1, k_wp2 - 1) + lhs_ma_zm(2 - 1, i - 1, k - 1)) +
               lhs_diff_zm(2 - 1, i - 1, k - 1));
          // LHS accumulation (ac) term and pressure term 2 (pr2).
          lhs(3 - 1, i - 1, k_wp2 - 1) =
              (lhs(3 - 1, i - 1, k_wp2 - 1) + lhs_ac_pr2_wp2(i - 1, k - 1));
          // LHS dissipation term 1 (dp1).
          // Note:  An "over-implicit" weighted time step is applied to this term.
          // A weighting factor of greater than 1 may be used to make the term
          // more numerically stable (see note below for w'^3 LHS turbulent
          // advection (ta) term).
          lhs(3 - 1, i - 1, k_wp2 - 1) =
              (lhs(3 - 1, i - 1, k_wp2 - 1) + (gamma_over_implicit_ts * lhs_dp1_wp2(i - 1, k - 1)));
          // LHS time tendency.
          lhs(3 - 1, i - 1, k_wp2 - 1) = (lhs(3 - 1, i - 1, k_wp2 - 1) + invrs_dt);
          // LHS turbulent advection (ta) term.
          lhs(4 - 1, i - 1, k_wp2 - 1) =
              (lhs(4 - 1, i - 1, k_wp2 - 1) + lhs_ta_wp2(2 - 1, i - 1, k - 1));
          // LHS mean advection (ma) and diffusion (diff) terms
          lhs(5 - 1, i - 1, k_wp2 - 1) =
              ((lhs(5 - 1, i - 1, k_wp2 - 1) + lhs_ma_zm(3 - 1, i - 1, k - 1)) +
               lhs_diff_zm(3 - 1, i - 1, k - 1));
        });
  }
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_lhs::cells_5", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          int k_wp3{};
          k_wp3 = (2 * k);
          // ------ w'3 ------
          // LHS mean advection (ma) and diffusion (diff) terms
          lhs(1 - 1, i - 1, k_wp3 - 1) =
              ((lhs(1 - 1, i - 1, k_wp3 - 1) + lhs_ma_zt(1 - 1, i - 1, k - 1)) +
               lhs_diff_zt(1 - 1, i - 1, k - 1));
          // LHS turbulent production (tp) term.
          // Note:  An "over-implicit" weighted time step is applied to this term.
          lhs(2 - 1, i - 1, k_wp3 - 1) =
              (lhs(2 - 1, i - 1, k_wp3 - 1) +
               (gamma_over_implicit_ts * lhs_tp_wp3(1 - 1, i - 1, k - 1)));
          // LHS mean advection (ma) and diffusion (diff) terms
          lhs(3 - 1, i - 1, k_wp3 - 1) =
              ((lhs(3 - 1, i - 1, k_wp3 - 1) + lhs_ma_zt(2 - 1, i - 1, k - 1)) +
               lhs_diff_zt(2 - 1, i - 1, k - 1));
          // LHS accumulation (ac) term and pressure term 2 (pr2).
          lhs(3 - 1, i - 1, k_wp3 - 1) =
              (lhs(3 - 1, i - 1, k_wp3 - 1) + lhs_ac_pr2_wp3(i - 1, k - 1));
          // LHS pressure term 1 (pr1).
          // Note:  An "over-implicit" weighted time step is applied to this term.
          lhs(3 - 1, i - 1, k_wp3 - 1) =
              (lhs(3 - 1, i - 1, k_wp3 - 1) + (gamma_over_implicit_ts * lhs_pr1_wp3(i - 1, k - 1)));
          // Add implicit splatting
          lhs(3 - 1, i - 1, k_wp3 - 1) =
              (lhs(3 - 1, i - 1, k_wp3 - 1) + lhs_splat_wp3(i - 1, k - 1));
          // LHS time tendency.
          lhs(3 - 1, i - 1, k_wp3 - 1) = (lhs(3 - 1, i - 1, k_wp3 - 1) + invrs_dt);
          // LHS turbulent production (tp) term.
          // Note:  An "over-implicit" weighted time step is applied to this term.
          lhs(4 - 1, i - 1, k_wp3 - 1) =
              (lhs(4 - 1, i - 1, k_wp3 - 1) +
               (gamma_over_implicit_ts * lhs_tp_wp3(2 - 1, i - 1, k - 1)));
          // LHS mean advection (ma) and diffusion (diff) terms
          lhs(5 - 1, i - 1, k_wp3 - 1) =
              ((lhs(5 - 1, i - 1, k_wp3 - 1) + lhs_ma_zt(3 - 1, i - 1, k - 1)) +
               lhs_diff_zt(3 - 1, i - 1, k - 1));
        });
  }
  // Upper (lower) boundary for w'3 on ascending (descending) grid.
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ndiags5) : (1) >= (ndiags5)) ? ((ndiags5) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_lhs::cells_6", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int b = (1) + ((cell / (1)) % extent_1) * (1);
          const int i = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          if ((b != 3)) {
            lhs(b - 1, i - 1, ((2 * nzm) - 2) - 1) = real_t(0.0);
          } else {
            lhs(b - 1, i - 1, ((2 * nzm) - 2) - 1) = real_t(1.0);
          }
        });
  }
  // Upper (lower) boundary for w'2 on ascending (desceding) grid.
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ndiags5) : (1) >= (ndiags5)) ? ((ndiags5) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_lhs::cells_7", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int b = (1) + ((cell / (1)) % extent_1) * (1);
          const int i = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          if ((b != 3)) {
            lhs(b - 1, i - 1, ((2 * nzm) - 1) - 1) = real_t(0.0);
          } else {
            lhs(b - 1, i - 1, ((2 * nzm) - 1) - 1) = real_t(1.0);
          }
        });
  }
  // LHS pressure term 1 (pr1) for wp2
  // Note:  An "over-implicit" weighted time step is applied to this term.
  // A weighting factor of greater than 1 may be used to make the term
  // more numerically stable (see note below for w'^3 LHS turbulent
  // advection (ta) term).
  // Reference:
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wp2_pr
  // Add terms to lhs
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_lhs::cells_8", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          int k_wp2{};
          k_wp2 = ((2 * k) - 1);
          lhs(3 - 1, i - 1, k_wp2 - 1) =
              (lhs(3 - 1, i - 1, k_wp2 - 1) + (gamma_over_implicit_ts * lhs_pr1_wp2(i - 1, k - 1)));
        });
  }
  // Add implicit splatting to wp2
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_lhs::cells_9", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          int k_wp2{};
          k_wp2 = ((2 * k) - 1);
          lhs(3 - 1, i - 1, k_wp2 - 1) =
              (lhs(3 - 1, i - 1, k_wp2 - 1) + lhs_splat_wp2(i - 1, k - 1));
        });
  }
  // LHS turbulent advection (ta) term for wp3
  if (!(l_explicit_turbulent_adv_wp3)) {
    // Note:  An "over-implicit" weighted time step is applied to this term.
    // The weight of the implicit portion of this term is controlled
    // by the factor gamma_over_implicit_ts (abbreviated "gamma" in
    // the expression below).  A factor is added to the right-hand
    // side of the equation in order to balance a weight that is not
    // equal to 1, such that:
    // -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
    // where X is the variable that is being solved for in a
    // predictive equation (w'^3 in this case), y(t) is the
    // linearized portion of the term that gets treated implicitly,
    // and RHS is the portion of the term that is always treated
    // explicitly (in the case of the w'^3 turbulent advection term,
    // RHS = 0).  A weight of greater than 1 can be applied to make
    // the term more numerically stable.
    // Add terms to lhs
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      const int extent_2 =
          ((1) > 0 ? (1) <= (ndiags5) : (1) >= (ndiags5)) ? ((ndiags5) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_lhs::cells_10", extent_0 * extent_1 * extent_2, KOKKOS_LAMBDA(const int cell) {
            const int b = (1) + ((cell / (1)) % extent_2) * (1);
            const int i = (1) + ((cell / (1 * extent_2)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_2 * extent_1)) % extent_0) * (1);
            int k_wp3{};
            k_wp3 = (2 * k);
            lhs(b - 1, i - 1, k_wp3 - 1) =
                (lhs(b - 1, i - 1, k_wp3 - 1) +
                 (gamma_over_implicit_ts * wp3_term_ta_lhs_result(b - 1, i - 1, k - 1)));
          });
    }
  }
  return;
}

// =============================================================================
void wp23_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, real_t dt, Real1D fcor_y,       // In
              Real3D wp3_term_ta_lhs_result,                                                 // In
              Real3D lhs_diff_zm, Real3D lhs_diff_zt, Real3D lhs_diff_zm_crank,              // In
              Real3D lhs_diff_zt_crank,                                                      // In
              Real3D lhs_tp_wp3, Real3D lhs_adv_tp_wp3, Real3D lhs_pr_tp_wp3,                // In
              Real3D lhs_ta_wp3, Real2D lhs_dp1_wp2, Real2D rhs_dp1_wp2, Real2D lhs_pr1_wp2, // In
              Real2D rhs_pr1_wp2, Real2D lhs_pr1_wp3, Real2D rhs_pr1_wp3,                    // In
              Real2D rhs_bp_pr2_wp2,                                                         // In
              Real2D rhs_pr_dfsn_wp2, Real2D rhs_bp1_pr2_wp3, Real2D rhs_pr3_wp2,            // In
              Real2D rhs_pr3_wp3,                                                            // In
              Real2D rhs_ta_wp3, Real2D rhs_pr_turb_wp3, Real2D rhs_pr_dfsn_wp3,             // In
              Real2D wp2, Real2D wp3, Real2D wpup2, Real2D wpvp2,                            // In
              Real2D wpthvp, Real2D wp2thvp, Real2D wp2up, Real2D up2, Real2D vp2,           // In
              Real2D upwp,                                                                   // In
              Real2D C11_Skw_fnc, Real2D thv_ds_zm, Real2D thv_ds_zt,                        // In
              Real2D lhs_splat_wp2, Real2D lhs_splat_wp3,                                    // In
              Real2D clubb_params,                                                           // In
              int iiPDF_type,                                                                // In
              bool l_use_tke_in_wp2_wp3_K_dfsn,                                              // In
              bool l_ho_nontrad_coriolis,                                                    // In
              f2k::stats_type& stats, // InOut
              Real2D rhs) {           // Out
  //
  // Description:
  //   Compute RHS vector for w'^2 and w'^3.
  //   This subroutine computes the explicit portion of
  //   the w'^2 and w'^3 equations.
  //
  //   Notes:
  //        For LHS turbulent advection (ta) term.
  //           An "over-implicit" weighted time step is applied to this term.
  //           The weight of the implicit portion of this term is controlled
  //           by the factor gamma_over_implicit_ts (abbreviated "gamma" in
  //           the expression below).  A factor is added to the right-hand
  //           side of the equation in order to balance a weight that is not
  //           equal to 1, such that:
  //                -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
  //           where X is the variable that is being solved for in a
  //           predictive equation (w'^3 in this case), y(t) is the
  //           linearized portion of the term that gets treated implicitly,
  //           and RHS is the portion of the term that is always treated
  //           explicitly (in the case of the w'^3 turbulent advection term,
  //           RHS = 0).  A weight of greater than 1 can be applied to make
  //           the term more numerically stable.
  //
  //
  //  WARNING: This subroutine has been optimized. Significant changes could
  //           noticeably  impact computational efficiency. See clubb:ticket:834
  //-------------------------------------------------------------------------------

  // --------------------- Input Variables ---------------------
  // dt: Timestep length                           [s]
  // fcor_y: Nontraditional Coriolis parameter         [s^-1]
  // Meridional planetary vorticity. Proportional to cos(latitude)
  // lhs_tp_wp3: Turbulent production terms of w'^3
  // lhs_adv_tp_wp3: Turbulent production terms of w'^3 (for stats)
  // lhs_pr_tp_wp3: Pressure scrambling terms for turbulent production of w'^3 (for stats)
  // lhs_ta_wp3: Turbulent advection terms for wp3
  // lhs_dp1_wp2: wp2 "over-implicit" dissipation term
  // rhs_dp1_wp2: wp2 rhs dissipation term
  // lhs_pr1_wp2: wp2 "over-implicit" pressure term 1
  // rhs_pr1_wp2: wp2 rhs pressure term 1
  // rhs_bp_pr2_wp2: wp2 bouyancy production and pressure term 2
  // rhs_pr_dfsn_wp2: wp2 pressure diffusion term
  // rhs_pr3_wp2: wp2 pressure term 3
  // lhs_pr1_wp3: wp3 "over-implicit" pressure term 1
  // rhs_pr1_wp3: wp3 rhs pressure term 1
  // rhs_bp1_pr2_wp3: wp3 bouyancy production 1 and pressure term 2
  // rhs_pr3_wp3: wp3 pressure term 3
  // rhs_ta_wp3: wp3 turbulent advection term
  // rhs_pr_turb_wp3: wp3 pressure-turbulence correlation term !--EXPERIMENTAL--!
  // rhs_pr_dfsn_wp3: wp3 pressure diffusion term
  // wp2: w'^2 (momentum levels)                    [m^2/s^2]
  // wpthvp: w'th_v' (momentum levels)                 [K m/s]
  // up2: u'^2 (momentum levels)                    [m^2/s^2]
  // vp2: v'^2 (momentum levels)                    [m^2/s^2]
  // upwp: u'w' (momentum levels)                    [m^2/s^2]
  // thv_ds_zm: Dry, base-state theta_v on momentum levs. [K]
  // lhs_splat_wp2: LHS coefficient of wp2 splatting term     [1/s]
  // wp3: w'^3 (thermodynamic levels)               [m^3/s^3]
  // wpup2: w'u'^2 (thermodynamic levels)             [m^3/s^3]
  // wpvp2: w'v'^2 (thermodynamic levels)             [m^3/s^3]
  // wp2thvp: w'^2th_v' (thermodynamic levels)        [K m^2/s^2]
  // wp2up: w'^2u' (thermodynamic levels)             [m^3/s^3]
  // C11_Skw_fnc: C_11 parameter with Sk_w applied          [-]
  // thv_ds_zt: Dry, base-state theta_v on thermo. levs.  [K]
  // lhs_splat_wp3: LHS coefficient of wp3 splatting term     [1/s]
  // clubb_params: Array of CLUBB's tunable parameters    [units vary]
  // iiPDF_type: Selected option for the two-component normal (double
  // Gaussian) PDF type to use for the w, rt, and theta-l (or
  // w, chi, and eta) portion of CLUBB's multivariate,
  // two-component PDF.
  // l_use_tke_in_wp2_wp3_K_dfsn: Use TKE in eddy diffusion for wp2 and wp3
  // l_ho_nontrad_coriolis: Flag to implement the nontraditional Coriolis terms in the
  // prognostic equations of <w'w'>, <u'w'>, and <u'u'>.
  // --------------------- Output Variable ---------------------
  // rhs: RHS of band matrix

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  int rhs_lb_idx_zm{}; // Index of (momentum) lower boundary level in rhs array
  int rhs_lb_idx_zt{}; // Index of (thermodynamic) lower boundary level in rhs array
  int rhs_ub_idx_zt{}; // Index of (thermodynamic) upper boundary level in rhs array
  int rhs_ub_idx_zm{}; // Index of (momentum) upper boundary level in rhs array
  Real2D rhs_bp_wp2("wp23_rhs::rhs_bp_wp2", ngrdcol, nzm);   // wp2 bouyancy production (stats only)
  Real2D rhs_pr2_wp2("wp23_rhs::rhs_pr2_wp2", ngrdcol, nzm); // wp2 pressure term 2 (stats only)
  Real2D rhs_bp1_wp3("wp23_rhs::rhs_bp1_wp3", ngrdcol,
                     nzm); // wp3 bouyancy production 1 (stats only)
  Real2D rhs_pr2_wp3("wp23_rhs::rhs_pr2_wp3", ngrdcol, nzm); // wp3 pressure term 2 (stats only)
  Real2D stats_tmp_zm("wp23_rhs::stats_tmp_zm", ngrdcol, nzm);
  Real2D stats_tmp_zt("wp23_rhs::stats_tmp_zt", ngrdcol, nzt);
  real_t invrs_dt{}; // Inverse of dt, 1/dt, used for computational efficiency
  Real2D C11_Skw_fnc_zeros("wp23_rhs::C11_Skw_fnc_zeros", ngrdcol, nzt);
  Real2D C11_Skw_fnc_plus_one("wp23_rhs::C11_Skw_fnc_plus_one", ngrdcol, nzt);
  Real1D C_uu_buoy_zeros("wp23_rhs::C_uu_buoy_zeros", ngrdcol);
  Real1D C_uu_buoy_plus_one("wp23_rhs::C_uu_buoy_plus_one", ngrdcol);
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_lb_zm = gr.k_lb_zm;

  // Calculate invers_dt
  invrs_dt = (real_t(1.0) / dt);
  // Initialize to zero
  {
    const int extent_0 = ((1) > 0 ? (1) <= (((2 * nzm) - 1)) : (1) >= (((2 * nzm) - 1)))
                             ? ((((2 * nzm) - 1)) - (1)) / (1) + 1
                             : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_rhs::cells_1", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          rhs(i - 1, k - 1) = real_t(0.0);
        });
  }
  // Experimental bouyancy term
  // Experimental term from CLUBB TRAC ticket #411
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_rhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          int k_wp3{};
          k_wp3 = (2 * k);
          rhs(i - 1, k_wp3 - 1) = (rhs_pr_turb_wp3(i - 1, k - 1) + rhs_pr_dfsn_wp3(i - 1, k - 1));
        });
  }
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_rhs::cells_3", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          int k_wp2{};
          k_wp2 = ((2 * k) - 1);
          rhs(i - 1, k_wp2 - 1) = (rhs(i - 1, k_wp2 - 1) + rhs_pr_dfsn_wp2(i - 1, k - 1));
        });
  }
  // These lines are for the diffusional term with a Crank-Nicholson
  // time step.  They are not used for completely implicit diffusion.
  if (l_crank_nich_diff) {
    // Add diffusion terms
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_4", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            int k_wp2{};
            k_wp2 = ((2 * k) - 1);
            rhs(i - 1, k_wp2 - 1) =
                (((rhs(i - 1, k_wp2 - 1) -
                   (lhs_diff_zm_crank((2 + gr_grid_dir_indx) - 1, i - 1, k - 1) *
                    wp2(i - 1, (k - gr_grid_dir_indx) - 1))) -
                  (lhs_diff_zm_crank(2 - 1, i - 1, k - 1) * wp2(i - 1, k - 1))) -
                 (lhs_diff_zm_crank((2 - gr_grid_dir_indx) - 1, i - 1, k - 1) *
                  wp2(i - 1, (k + gr_grid_dir_indx) - 1)));
          });
    }
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_5", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            int k_wp3{};
            k_wp3 = (2 * k);
            rhs(i - 1, k_wp3 - 1) =
                (((rhs(i - 1, k_wp3 - 1) -
                   (lhs_diff_zt_crank((2 + gr_grid_dir_indx) - 1, i - 1, k - 1) *
                    wp3(i - 1, (k - gr_grid_dir_indx) - 1))) -
                  (lhs_diff_zt_crank(2 - 1, i - 1, k - 1) * wp3(i - 1, k - 1))) -
                 (lhs_diff_zt_crank((2 - gr_grid_dir_indx) - 1, i - 1, k - 1) *
                  wp3(i - 1, (k + gr_grid_dir_indx) - 1)));
          });
    }
  }
  // This code block adds terms to the right-hand side so that TKE is being
  // used in eddy diffusion instead of just wp2 or wp3.  For example, in the
  // wp2 equation, if this flag is false, the eddy diffusion term would
  // normally be completely implicit (hence no right-hand side contribution),
  // and equal to +d/dz((K+nu)d/dz(wp2)).  With this flag set to true, the eddy
  // diffusion term will be +d/dz((K+nu)d/dz(up2+vp2+wp2)), but the up2 and vp2
  // parts are added on here as if they were right-hand side terms. For the wp3
  // equation, with this flag false, the eddy diffusion term is
  // +d/dz((K+nu)d/dz(wp3)), but with this flag true, it will be
  // +d/dz((K+nu)d/dz(wpup2+wpvp2+wp3)).
  if (l_use_tke_in_wp2_wp3_K_dfsn) {
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_6", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            int k_wp2{};
            k_wp2 = ((2 * k) - 1);
            rhs(i - 1, k_wp2 - 1) =
                (((rhs(i - 1, k_wp2 - 1) - (lhs_diff_zm((2 + gr_grid_dir_indx) - 1, i - 1, k - 1) *
                                            (up2(i - 1, (k - gr_grid_dir_indx) - 1) +
                                             vp2(i - 1, (k - gr_grid_dir_indx) - 1)))) -
                  (lhs_diff_zm(2 - 1, i - 1, k - 1) * (up2(i - 1, k - 1) + vp2(i - 1, k - 1)))) -
                 (lhs_diff_zm((2 - gr_grid_dir_indx) - 1, i - 1, k - 1) *
                  (up2(i - 1, (k + gr_grid_dir_indx) - 1) +
                   vp2(i - 1, (k + gr_grid_dir_indx) - 1))));
          });
    }
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_7", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            int k_wp3{};
            k_wp3 = (2 * k);
            rhs(i - 1, k_wp3 - 1) =
                (((rhs(i - 1, k_wp3 - 1) - (lhs_diff_zt((2 + gr_grid_dir_indx) - 1, i - 1, k - 1) *
                                            (wpup2(i - 1, (k - gr_grid_dir_indx) - 1) +
                                             wpvp2(i - 1, (k - gr_grid_dir_indx) - 1)))) -
                  (lhs_diff_zt(2 - 1, i - 1, k - 1) *
                   (wpup2(i - 1, k - 1) + wpvp2(i - 1, k - 1)))) -
                 (lhs_diff_zt((2 + gr_grid_dir_indx) - 1, i - 1, k - 1) *
                  (wpup2(i - 1, (k + gr_grid_dir_indx) - 1) +
                   wpvp2(i - 1, (k + gr_grid_dir_indx) - 1))));
          });
    }
  }
  // Add pressure terms and splat terms
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_rhs::cells_8", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          int k_wp2{};
          k_wp2 = ((2 * k) - 1);
          rhs(i - 1, k_wp2 - 1) = (rhs(i - 1, k_wp2 - 1) + rhs_pr1_wp2(i - 1, k - 1));
          rhs(i - 1, k_wp2 - 1) =
              (rhs(i - 1, k_wp2 - 1) + ((one - gamma_over_implicit_ts) *
                                        (-(lhs_pr1_wp2(i - 1, k - 1)) * wp2(i - 1, k - 1))));
          // Effect of vertical compression of eddies
          rhs(i - 1, k_wp2 - 1) = rhs(i - 1, k_wp2 - 1);
        });
  }
  if (l_ho_nontrad_coriolis) {
    // Add the nontraditional Coriolis term for wp2
    // Hing Ong, 19 July 2025
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_9", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            int k_wp2{};
            k_wp2 = ((2 * k) - 1);
            rhs(i - 1, k_wp2 - 1) =
                (rhs(i - 1, k_wp2 - 1) + ((two * fcor_y(i - 1)) * upwp(i - 1, k - 1)));
          });
    }
    // Add the nontraditional Coriolis term for wp3
    // Hing Ong, 1 Septempber 2025
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_10", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            int k_wp3{};
            k_wp3 = (2 * k);
            rhs(i - 1, k_wp3 - 1) =
                (rhs(i - 1, k_wp3 - 1) + ((three * fcor_y(i - 1)) * wp2up(i - 1, k - 1)));
          });
    }
  }
  // l_ho_nontrad_coriolis
  // Combine terms
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_rhs::cells_11", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          int k_wp3{};
          k_wp3 = (2 * k);
          // ------ Combine terms for 3rd moment of vertical velocity, <w'^3> ------ !
          // RHS time tendency.
          rhs(i - 1, k_wp3 - 1) = (rhs(i - 1, k_wp3 - 1) + (invrs_dt * wp3(i - 1, k - 1)));
          // RHS contribution from "over-implicit" turbulent production (tp) term.
          rhs(i - 1, k_wp3 - 1) =
              (rhs(i - 1, k_wp3 - 1) +
               ((one - gamma_over_implicit_ts) *
                ((-(lhs_tp_wp3(1 - 1, i - 1, k - 1)) * wp2(i - 1, (k + 1) - 1)) -
                 (lhs_tp_wp3(2 - 1, i - 1, k - 1) * wp2(i - 1, k - 1)))));
          // RHS buoyancy production (bp) term and pressure term 2 (pr2).
          rhs(i - 1, k_wp3 - 1) = (rhs(i - 1, k_wp3 - 1) + rhs_bp1_pr2_wp3(i - 1, k - 1));
          // RHS pressure term 1
          rhs(i - 1, k_wp3 - 1) = (rhs(i - 1, k_wp3 - 1) + rhs_pr1_wp3(i - 1, k - 1));
          // RHS "over implicit" pressure term 1 (pr1).
          rhs(i - 1, k_wp3 - 1) =
              (rhs(i - 1, k_wp3 - 1) + ((one - gamma_over_implicit_ts) *
                                        (-(lhs_pr1_wp3(i - 1, k - 1)) * wp3(i - 1, k - 1))));
        });
  }
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_rhs::cells_12", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          int k_wp2{};
          k_wp2 = ((2 * k) - 1);
          // ------ Combine terms for 2nd moment of vertical velocity, <w'^2> ------ !
          // RHS time tendency.
          rhs(i - 1, k_wp2 - 1) = (rhs(i - 1, k_wp2 - 1) + (invrs_dt * wp2(i - 1, k - 1)));
          // RHS buoyancy production (bp) term and pressure term 2 (pr2).
          rhs(i - 1, k_wp2 - 1) = (rhs(i - 1, k_wp2 - 1) + rhs_bp_pr2_wp2(i - 1, k - 1));
          // RHS pressure term 3 (pr3).
          rhs(i - 1, k_wp2 - 1) = (rhs(i - 1, k_wp2 - 1) + rhs_pr3_wp2(i - 1, k - 1));
          // RHS dissipation term 1 (dp1).
          rhs(i - 1, k_wp2 - 1) = (rhs(i - 1, k_wp2 - 1) + rhs_dp1_wp2(i - 1, k - 1));
          // RHS "over implicit" pressure term 1 (pr1).
          rhs(i - 1, k_wp2 - 1) =
              (rhs(i - 1, k_wp2 - 1) + ((one - gamma_over_implicit_ts) *
                                        (-(lhs_dp1_wp2(i - 1, k - 1)) * wp2(i - 1, k - 1))));
        });
  }
  if (l_explicit_turbulent_adv_wp3) {
    // The turbulent advection term is being solved explicitly.
    // Add RHS turbulent advection (ta) terms
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_13", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            int k_wp3{};
            k_wp3 = (2 * k);
            rhs(i - 1, k_wp3 - 1) = (rhs(i - 1, k_wp3 - 1) + rhs_ta_wp3(i - 1, k - 1));
          });
    }
  } else {
    // The turbulent advection term is being solved implicitly. See note above
    if ((iiPDF_type == iiPDF_ADG1)) {
      // The ADG1 PDF is used.
      // Add terms
      {
        const int extent_0 =
            ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
        const int extent_1 =
            ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
        Kokkos::parallel_for(
            "wp23_rhs::cells_14", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
              const int i = (1) + ((cell / (1)) % extent_1) * (1);
              const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
              int k_wp3{};
              k_wp3 = (2 * k);
              // Brian -- come up with better method for testing ascending vs descending
              if ((gr_grid_dir_indx > 0)) {
                // Ascending grid
                rhs(i - 1, k_wp3 - 1) =
                    (rhs(i - 1, k_wp3 - 1) +
                     ((one - gamma_over_implicit_ts) *
                      (((((-(wp3_term_ta_lhs_result(1 - 1, i - 1, k - 1)) *
                           wp3(i - 1, (k + 1) - 1)) -
                          (wp3_term_ta_lhs_result(2 - 1, i - 1, k - 1) * wp2(i - 1, (k + 1) - 1))) -
                         (wp3_term_ta_lhs_result(3 - 1, i - 1, k - 1) * wp3(i - 1, k - 1))) -
                        (wp3_term_ta_lhs_result(4 - 1, i - 1, k - 1) * wp2(i - 1, k - 1))) -
                       (wp3_term_ta_lhs_result(5 - 1, i - 1, k - 1) * wp3(i - 1, (k - 1) - 1)))));
              } else {
                // gr%grid_dir_indx < 0
                // Descending grid
                rhs(i - 1, k_wp3 - 1) =
                    (rhs(i - 1, k_wp3 - 1) +
                     ((one - gamma_over_implicit_ts) *
                      (((((-(wp3_term_ta_lhs_result(5 - 1, i - 1, k - 1)) *
                           wp3(i - 1, (k - 1) - 1)) -
                          (wp3_term_ta_lhs_result(4 - 1, i - 1, k - 1) * wp2(i - 1, k - 1))) -
                         (wp3_term_ta_lhs_result(3 - 1, i - 1, k - 1) * wp3(i - 1, k - 1))) -
                        (wp3_term_ta_lhs_result(2 - 1, i - 1, k - 1) * wp2(i - 1, (k + 1) - 1))) -
                       (wp3_term_ta_lhs_result(1 - 1, i - 1, k - 1) * wp3(i - 1, (k + 1) - 1)))));
              }
            });
      }
    } else if (((iiPDF_type == iiPDF_new) || (iiPDF_type == iiPDF_new_hybrid))) {
      // The new PDF or the new hybrid PDF is used.
      // Add terms
      {
        const int extent_0 =
            ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
        const int extent_1 =
            ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
        Kokkos::parallel_for(
            "wp23_rhs::cells_15", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
              const int i = (1) + ((cell / (1)) % extent_1) * (1);
              const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
              int k_wp3{};
              k_wp3 = (2 * k);
              rhs(i - 1, k_wp3 - 1) =
                  (rhs(i - 1, k_wp3 - 1) +
                   ((one - gamma_over_implicit_ts) *
                    ((-(lhs_ta_wp3(1 - 1, i - 1, k - 1)) * wp2(i - 1, (k + 1) - 1)) -
                     (lhs_ta_wp3(2 - 1, i - 1, k - 1) * wp2(i - 1, k - 1)))));
            });
      }
    }
    // iiPDF_type
  }
  // l_explicit_turbulent_adv_wp3
  // --------- Boundary Conditions ---------
  // Both wp2 and wp3 used fixed-point boundary conditions.
  // Therefore, anything set in the above loop at both the upper
  // and lower boundaries would be overwritten here.  However, the
  // above loop does not extend to the boundary levels.  An array
  // with a value of 1 at the main diagonal on the left-hand side
  // and with values of 0 at all other diagonals on the left-hand
  // side will preserve the right-hand side value at that level.
  // The value of w'^2 at the lower boundary will remain the same.
  // When the lower boundary is at the surface, the surface value of
  // w'^2 is set in subroutine calc_surface_varnce (surface_varnce_module.F).
  // The value of w'^3 at the lower boundary will be 0.
  // The value of w'^2 at the upper boundary will be set to the threshold
  // minimum value of w_tol_sqd.
  // The value of w'^3 at the upper boundary will be set to 0.
  if ((gr_grid_dir_indx > 0)) {
    // Ascending Grid
    rhs_lb_idx_zm = 1;
    rhs_lb_idx_zt = 2;
    rhs_ub_idx_zt = (2 * nzt);
    rhs_ub_idx_zm = ((2 * nzm) - 1);
  } else {
    // Descending Grid
    rhs_lb_idx_zm = ((2 * nzm) - 1);
    rhs_lb_idx_zt = (2 * nzt);
    rhs_ub_idx_zt = 2;
    rhs_ub_idx_zm = 1;
  }
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp23_rhs::cells_16", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Lower Boundary
          rhs(i - 1, rhs_lb_idx_zm - 1) = wp2(i - 1, gr_k_lb_zm - 1);
          rhs(i - 1, rhs_lb_idx_zt - 1) = real_t(0.0);
          // Upper Boundary
          rhs(i - 1, rhs_ub_idx_zt - 1) = real_t(0.0);
          rhs(i - 1, rhs_ub_idx_zm - 1) = w_tol_sqd;
        });
  }
  // --------- Statistics output ---------
  if (stats.flags().l_sample) {
    Kokkos::deep_copy(C_uu_buoy_zeros, zero);
    {
      Kokkos::parallel_for(
          "wp23_rhs::array_0", ngrdcol, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            C_uu_buoy_plus_one(i - 1) = (clubb_params(i - 1, iC_uu_buoy) + one);
          });
    }
    Kokkos::deep_copy(C11_Skw_fnc_zeros, zero);
    {
      Kokkos::parallel_for(
          "wp23_rhs::array_1", ngrdcol * nzt, KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol + 1;
            const int k = cell / ngrdcol + 1;
            C11_Skw_fnc_plus_one(i - 1, k - 1) = (C11_Skw_fnc(i - 1, k - 1) + one);
          });
    }
    // w'^2 term bp is completely explicit; call stat_update_var_pt.
    // Note:  To find the contribution of w'^2 term bp, substitute 0 for the
    // C_uu_buoy input to function wp2_terms_bp_pr2_rhs.
    wp2_terms_bp_pr2_rhs(nzm, ngrdcol, gr, C_uu_buoy_zeros, // In
                         thv_ds_zm, wpthvp,                 // In
                         rhs_bp_wp2);                       // Out

    // w'^2 term pr2 has both implicit and explicit components; call
    // stat_begin_update_pt.  Since stat_begin_update_pt automatically
    // subtracts the value sent in, reverse the sign on wp2_terms_bp_pr2_rhs.
    // Note:  To find the contribution of w'^2 term pr2, add 1 to the
    // C_uu_buoy input to function wp2_terms_bp_pr2_rhs.
    wp2_terms_bp_pr2_rhs(nzm, ngrdcol, gr, C_uu_buoy_plus_one, // In
                         thv_ds_zm, wpthvp,                    // In
                         rhs_pr2_wp2);                         // Out

    // w'^3 term bp is completely explicit; call stat_update_var_pt.
    // Note:  To find the contribution of w'^3 term bp, substitute 0 for the
    // C_11 skewness function input to function wp3_terms_bp1_pr2_rhs.
    wp3_terms_bp1_pr2_rhs(nzt, ngrdcol, gr, C11_Skw_fnc_zeros, // In
                          thv_ds_zt, wp2thvp,                  // In
                          rhs_bp1_wp3);                        // Out

    // w'^3 term pr2 has both implicit and explicit components; call
    // stat_begin_update_pt.  Since stat_begin_update_pt automatically
    // subtracts the value sent in, reverse the sign on wp3_terms_bp1_pr2_rhs.
    // Note:  To find the contribution of w'^3 term pr2, add 1 to the
    // C_11 skewness function input to function wp3_terms_bp1_pr2_rhs.
    wp3_terms_bp1_pr2_rhs(nzt, ngrdcol, gr, C11_Skw_fnc_plus_one, // In
                          thv_ds_zt, wp2thvp,                     // In
                          rhs_pr2_wp3);                           // Out
  }
  if (stats.flags().l_sample) {
    // ----------- w'2 -----------
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_17", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) = rhs_bp_wp2(i - 1, k - 1);
          });
    }
    f2k::stats_update("wp2_bp", stats_tmp_zm, stats);
    if (l_ho_nontrad_coriolis) {
      Kokkos::deep_copy(stats_tmp_zm, zero);
      {
        const int extent_0 =
            ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
        const int extent_1 =
            ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
        Kokkos::parallel_for(
            "wp23_rhs::cells_18", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
              const int i = (1) + ((cell / (1)) % extent_1) * (1);
              const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
              stats_tmp_zm(i - 1, k - 1) = ((two * fcor_y(i - 1)) * upwp(i - 1, k - 1));
            });
      }
      // w'^2 term nct is completely explicit; call stat_update_var_pt.
      // Hing Ong, 22 July 2025
      // Include effect of vertical compression of eddies in wp2 budget
      f2k::stats_update("wp2_nct", stats_tmp_zm, stats);
    }
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_19", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) = rhs_pr_dfsn_wp2(i - 1, k - 1);
          });
    }
    f2k::stats_update("wp2_pr_dfsn", stats_tmp_zm, stats);
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_20", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) = (-(lhs_splat_wp2(i - 1, k - 1)) * wp2(i - 1, k - 1));
          });
    }
    f2k::stats_update("wp2_splat", stats_tmp_zm, stats);
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_21", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) = rhs_pr3_wp2(i - 1, k - 1);
          });
    }
    // w'^2 term pr3 is completely explicit; call stat_update_var_pt.
    f2k::stats_update("wp2_pr3", stats_tmp_zm, stats);
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_22", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) = -(rhs_pr2_wp2(i - 1, k - 1));
          });
    }
    f2k::stats_begin_budget("wp2_pr2", stats_tmp_zm, stats);
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_23", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) = -(rhs_dp1_wp2(i - 1, k - 1));
          });
    }
    // w'^2 term dp1 has both implicit and explicit components; call
    // stat_begin_update_pt.  Since stat_begin_update_pt automatically
    // subtracts the value sent in, reverse the sign on wp2_term_dp1_rhs.
    f2k::stats_begin_budget("wp2_dp1", stats_tmp_zm, stats);
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_24", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) = ((one - gamma_over_implicit_ts) *
                                          (-(lhs_dp1_wp2(i - 1, k - 1)) * wp2(i - 1, k - 1)));
          });
    }
    f2k::stats_update_budget("wp2_dp1", stats_tmp_zm, stats);
    if ((l_crank_nich_diff || l_use_tke_in_wp2_wp3_K_dfsn)) {
      // w'^2 term dp2 and w'^3 term dp1 have both implicit and explicit
      // components (if the l_use_tke_in_wp2_wp3_K_dfsn flag is true;
      // call stat_begin_update_pt.
      // Crank-Nicholson scheme is selected); call stat_begin_update_pt.
      // Since stat_begin_update_pt automatically subtracts the value sent in,
      // reverse the sign on right-hand side diffusion component.  If
      // Crank-Nicholson diffusion is not selected, the stat_begin_update_pt
      // will not be called.
      Kokkos::deep_copy(stats_tmp_zm, zero);
      if (l_crank_nich_diff) {
        {
          const int extent_0 = ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1)))
                                   ? (((nzm - 1)) - (2)) / (1) + 1
                                   : 0;
          const int extent_1 =
              ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
          Kokkos::parallel_for(
              "wp23_rhs::cells_25", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
                const int i = (1) + ((cell / (1)) % extent_1) * (1);
                const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
                stats_tmp_zm(i - 1, k - 1) =
                    (((stats_tmp_zm(i - 1, k - 1) +
                       (lhs_diff_zm_crank((2 + gr_grid_dir_indx) - 1, i - 1, k - 1) *
                        wp2(i - 1, (k - gr_grid_dir_indx) - 1))) +
                      (lhs_diff_zm_crank(2 - 1, i - 1, k - 1) * wp2(i - 1, k - 1))) +
                     (lhs_diff_zm_crank((2 - gr_grid_dir_indx) - 1, i - 1, k - 1) *
                      wp2(i - 1, (k + gr_grid_dir_indx) - 1)));
              });
        }
      }
      if (l_use_tke_in_wp2_wp3_K_dfsn) {
        {
          const int extent_0 = ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1)))
                                   ? (((nzm - 1)) - (2)) / (1) + 1
                                   : 0;
          const int extent_1 =
              ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
          Kokkos::parallel_for(
              "wp23_rhs::cells_26", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
                const int i = (1) + ((cell / (1)) % extent_1) * (1);
                const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
                stats_tmp_zm(i - 1, k - 1) =
                    (((stats_tmp_zm(i - 1, k - 1) +
                       (lhs_diff_zm((2 + gr_grid_dir_indx) - 1, i - 1, k - 1) *
                        (up2(i - 1, (k - gr_grid_dir_indx) - 1) +
                         vp2(i - 1, (k - gr_grid_dir_indx) - 1)))) +
                      (lhs_diff_zm(2 - 1, i - 1, k - 1) *
                       (up2(i - 1, k - 1) + vp2(i - 1, k - 1)))) +
                     (lhs_diff_zm((2 - gr_grid_dir_indx) - 1, i - 1, k - 1) *
                      (up2(i - 1, (k + gr_grid_dir_indx) - 1) +
                       vp2(i - 1, (k + gr_grid_dir_indx) - 1))));
              });
        }
      }
      f2k::stats_begin_budget("wp2_dp2", stats_tmp_zm, stats);
    }
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_27", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) = -(rhs_pr1_wp2(i - 1, k - 1));
          });
    }
    f2k::stats_begin_budget("wp2_pr1", stats_tmp_zm, stats);
    Kokkos::deep_copy(stats_tmp_zm, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_28", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zm(i - 1, k - 1) = ((one - gamma_over_implicit_ts) *
                                          (-(lhs_pr1_wp2(i - 1, k - 1)) * wp2(i - 1, k - 1)));
          });
    }
    // w'^2 term pr1 has both implicit and explicit components; call
    // stat_begin_update_pt.  Since stat_begin_update_pt automatically
    // subtracts the value sent in, reverse the sign on wp2_term_pr1_rhs.
    // Note:  An "over-implicit" weighted time step is applied to this
    // term.  A weighting factor of greater than 1 may be used to
    // term more numerically stable (see note below for w'^3 LHS
    // make the term more numerically stable (see note below for
    // w'^3 RHS turbulent advection (ta) term).
    f2k::stats_update_budget("wp2_pr1", stats_tmp_zm, stats);
    // ----------- w'3 -----------
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_29", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) = rhs_bp1_wp3(i - 1, k - 1);
          });
    }
    f2k::stats_update("wp3_bp1", stats_tmp_zt, stats);
    if (l_ho_nontrad_coriolis) {
      Kokkos::deep_copy(stats_tmp_zt, zero);
      {
        const int extent_0 =
            ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
        const int extent_1 =
            ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
        Kokkos::parallel_for(
            "wp23_rhs::cells_30", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
              const int i = (1) + ((cell / (1)) % extent_1) * (1);
              const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
              stats_tmp_zt(i - 1, k - 1) = ((three * fcor_y(i - 1)) * wp2up(i - 1, k - 1));
            });
      }
      // w'^3 term nct is completely explicit; call stat_update_var_pt.
      // Hing Ong, 1 September 2025
      f2k::stats_update("wp3_nct", stats_tmp_zt, stats);
    }
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_31", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) = rhs_pr_turb_wp3(i - 1, k - 1);
          });
    }
    f2k::stats_update("wp3_pr_turb", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_32", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) = rhs_pr_dfsn_wp3(i - 1, k - 1);
          });
    }
    f2k::stats_update("wp3_pr_dfsn", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_33", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) = (-(lhs_splat_wp3(i - 1, k - 1)) * wp3(i - 1, k - 1));
          });
    }
    f2k::stats_update("wp3_splat", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_34", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) = -(rhs_pr2_wp3(i - 1, k - 1));
          });
    }
    f2k::stats_begin_budget("wp3_pr2", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_35", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) = -(rhs_pr1_wp3(i - 1, k - 1));
          });
    }
    // w'^3 term pr1 has both implicit and explicit components; call
    // stat_begin_update_pt.  Since stat_begin_update_pt automatically
    // subtracts the value sent in, reverse the sign on wp3_term_pr1_rhs.
    f2k::stats_begin_budget("wp3_pr1", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_36", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) = ((one - gamma_over_implicit_ts) *
                                          (-(lhs_pr1_wp3(i - 1, k - 1)) * wp3(i - 1, k - 1)));
          });
    }
    f2k::stats_update_budget("wp3_pr1", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    if (l_explicit_turbulent_adv_wp3) {
      {
        const int extent_0 =
            ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
        const int extent_1 =
            ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
        Kokkos::parallel_for(
            "wp23_rhs::cells_37", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
              const int i = (1) + ((cell / (1)) % extent_1) * (1);
              const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
              stats_tmp_zt(i - 1, k - 1) = -(rhs_ta_wp3(i - 1, k - 1));
            });
      }
    } else {
      if ((iiPDF_type == iiPDF_ADG1)) {
        {
          const int extent_0 = ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1)))
                                   ? (((nzt - 1)) - (2)) / (1) + 1
                                   : 0;
          const int extent_1 =
              ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
          Kokkos::parallel_for(
              "wp23_rhs::cells_38", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
                const int i = (1) + ((cell / (1)) % extent_1) * (1);
                const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
                stats_tmp_zt(i - 1, k - 1) =
                    (-((one - gamma_over_implicit_ts)) *
                     (((((-(wp3_term_ta_lhs_result(1 - 1, i - 1, k - 1)) *
                          wp3(i - 1, (k + 1) - 1)) -
                         (wp3_term_ta_lhs_result(2 - 1, i - 1, k - 1) * wp2(i - 1, (k + 1) - 1))) -
                        (wp3_term_ta_lhs_result(3 - 1, i - 1, k - 1) * wp3(i - 1, k - 1))) -
                       (wp3_term_ta_lhs_result(4 - 1, i - 1, k - 1) * wp2(i - 1, k - 1))) -
                      (wp3_term_ta_lhs_result(5 - 1, i - 1, k - 1) * wp3(i - 1, (k - 1) - 1))));
              });
        }
      } else if (((iiPDF_type == iiPDF_new) || (iiPDF_type == iiPDF_new_hybrid))) {
        {
          const int extent_0 = ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1)))
                                   ? (((nzt - 1)) - (2)) / (1) + 1
                                   : 0;
          const int extent_1 =
              ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
          Kokkos::parallel_for(
              "wp23_rhs::cells_39", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
                const int i = (1) + ((cell / (1)) % extent_1) * (1);
                const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
                stats_tmp_zt(i - 1, k - 1) =
                    (-((one - gamma_over_implicit_ts)) *
                     ((-(lhs_ta_wp3(1 - 1, i - 1, k - 1)) * wp2(i - 1, (k + 1) - 1)) -
                      (lhs_ta_wp3(2 - 1, i - 1, k - 1) * wp2(i - 1, k - 1))));
              });
        }
      }
    }
    f2k::stats_begin_budget("wp3_ta", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_40", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) =
                (-((one - gamma_over_implicit_ts)) *
                 ((-(lhs_adv_tp_wp3(1 - 1, i - 1, k - 1)) * wp2(i - 1, (k + 1) - 1)) -
                  (lhs_adv_tp_wp3(2 - 1, i - 1, k - 1) * wp2(i - 1, k - 1))));
          });
    }
    f2k::stats_begin_budget("wp3_tp", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_41", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) =
                (-((one - gamma_over_implicit_ts)) *
                 ((-(lhs_pr_tp_wp3(1 - 1, i - 1, k - 1)) * wp2(i - 1, (k + 1) - 1)) -
                  (lhs_pr_tp_wp3(2 - 1, i - 1, k - 1) * wp2(i - 1, k - 1))));
          });
    }
    f2k::stats_begin_budget("wp3_pr_tp", stats_tmp_zt, stats);
    Kokkos::deep_copy(stats_tmp_zt, zero);
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp23_rhs::cells_42", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            stats_tmp_zt(i - 1, k - 1) = rhs_pr3_wp3(i - 1, k - 1);
          });
    }
    f2k::stats_begin_budget("wp3_pr3", stats_tmp_zt, stats);
    if ((l_crank_nich_diff || l_use_tke_in_wp2_wp3_K_dfsn)) {
      Kokkos::deep_copy(stats_tmp_zt, zero);
      if (l_crank_nich_diff) {
        {
          const int extent_0 = ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1)))
                                   ? (((nzt - 1)) - (2)) / (1) + 1
                                   : 0;
          const int extent_1 =
              ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
          Kokkos::parallel_for(
              "wp23_rhs::cells_43", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
                const int i = (1) + ((cell / (1)) % extent_1) * (1);
                const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
                stats_tmp_zt(i - 1, k - 1) =
                    (((stats_tmp_zt(i - 1, k - 1) +
                       (lhs_diff_zt((2 + gr_grid_dir_indx) - 1, i - 1, k - 1) *
                        wp3(i - 1, (k - gr_grid_dir_indx) - 1))) +
                      (lhs_diff_zt(2 - 1, i - 1, k - 1) * wp3(i - 1, k - 1))) +
                     (lhs_diff_zt((2 - gr_grid_dir_indx) - 1, i - 1, k - 1) *
                      wp3(i - 1, (k + gr_grid_dir_indx) - 1)));
              });
        }
      }
      if (l_use_tke_in_wp2_wp3_K_dfsn) {
        {
          const int extent_0 = ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1)))
                                   ? (((nzt - 1)) - (2)) / (1) + 1
                                   : 0;
          const int extent_1 =
              ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
          Kokkos::parallel_for(
              "wp23_rhs::cells_44", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
                const int i = (1) + ((cell / (1)) % extent_1) * (1);
                const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
                stats_tmp_zt(i - 1, k - 1) =
                    (((stats_tmp_zt(i - 1, k - 1) +
                       (lhs_diff_zt((2 + gr_grid_dir_indx) - 1, i - 1, k - 1) *
                        (wpup2(i - 1, (k - gr_grid_dir_indx) - 1) +
                         wpvp2(i - 1, (k - gr_grid_dir_indx) - 1)))) +
                      (lhs_diff_zt(2 - 1, i - 1, k - 1) *
                       (wpup2(i - 1, k - 1) + wpvp2(i - 1, k - 1)))) +
                     (lhs_diff_zt((2 - gr_grid_dir_indx) - 1, i - 1, k - 1) *
                      (wpup2(i - 1, (k + gr_grid_dir_indx) - 1) +
                       wpvp2(i - 1, (k + gr_grid_dir_indx) - 1))));
              });
        }
      }
      f2k::stats_begin_budget("wp3_dp1", stats_tmp_zt, stats);
    }
  }
  return;
}

// =============================================================================
void wp2_term_ta_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                     Real2D rho_ds_zt, Real2D invrs_rho_ds_zm,      // In
                     Real3D lhs_ta_wp2) {                           // Out
  //
  // Description:
  // Turbulent advection term for w'^2:  implicit portion of the code.
  //
  // The d(w'^2)/dt equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * w'^3 )/dz.
  //
  // The term is solved for completely implicitly, such that:
  //
  // - (1/rho_ds) * d( rho_ds * w'^3(t+1) )/dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign
  //        is reversed and the leading "-" in front of the term is changed
  //        to a "+".
  //
  // The timestep index (t+1) means that the value of w'^3 being used is from
  // the next timestep, which is being advanced to in solving the d(w'^2)/dt
  // and d(w'^3)/dt equations.
  //
  // This term is discretized as follows:
  //
  // While the values of w'^2 are found on the momentum levels, the values of
  // w'^3 are found on the thermodynamic levels.  Additionally, the values of
  // rho_ds_zt are found on the thermodynamic levels, and the values of
  // invrs_rho_ds_zm are found on the momentum levels.  On the thermodynamic
  // levels, the values of rho_ds_zt are multiplied by the values of w'^3.  The
  // derivative of (rho_ds_zt * w'^3) is taken over the intermediate (central)
  // momentum level, where it is multiplied by invrs_rho_ds_zm, yielding the
  // desired results.
  //
  // -----rho_ds_zt----------wp3------------------------------ t(k)
  //
  // ========invrs_rho_ds_zm==========d(rho_ds*wp3)/dz======== m(k)
  //
  // -----rho_ds_zt----------wp3------------------------------ t(k-1)
  //
  // The vertical indices t(k), m(k), and t(k-1) correspond with altitudes
  // zt(k), zm(k), and zt(k-1), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )
  //
  // References:
  //-----------------------------------------------------------------------

  // Constant parameters
  // ------------------------ Input Variables ------------------------
  // rho_ds_zt: Dry, static density at thermodynamic levels  [kg/m^3]
  // invrs_rho_ds_zm: Inv. dry, static density at momentum levels  [m^3/kg]
  // ------------------------ Return Variable ------------------------
  // lhs_ta_wp2: LHS coefficient of wp2 turbulent advection  [1/m]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  constexpr int kp1_tdiag = 1; // Thermodynamic superdiagonal index.
  constexpr int k_tdiag = 2;   // Thermodynamic subdiagonal index.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;
  const auto gr_invrs_dzm = gr.invrs_dzm;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    const int extent_1 = ((1) > 0 ? (1) <= (2) : (1) >= (2)) ? ((2) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_ta_lhs::cells_1", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Set lower boundary to 0
          lhs_ta_wp2(i - 1, k - 1, gr_k_lb_zm - 1) = zero;
          // Set upper boundary to 0
          lhs_ta_wp2(i - 1, k - 1, gr_k_ub_zm - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_ta_lhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Thermodynamic superdiagonal: [ x wp3(k,<t+1>) ]
          lhs_ta_wp2(kp1_tdiag - 1, i - 1, k - 1) =
              ((+(invrs_rho_ds_zm(i - 1, k - 1)) * gr_invrs_dzm(i - 1, k - 1)) *
               rho_ds_zt(i - 1, k - 1));
          // Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
          lhs_ta_wp2(k_tdiag - 1, i - 1, k - 1) =
              ((-(invrs_rho_ds_zm(i - 1, k - 1)) * gr_invrs_dzm(i - 1, k - 1)) *
               rho_ds_zt(i - 1, (k - 1) - 1));
        });
  }
  return;
}

// =============================================================================
void wp2_terms_ac_pr2_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                          Real1D C_uu_shr, Real2D wm_zt,                 // In
                          Real2D lhs_ac_pr2_wp2) {                       // Out
  //
  // Description:
  // Accumulation of w'^2 and w'^2 pressure term 2:  implicit portion of the
  // code.
  //
  // The d(w'^2)/dt equation contains an accumulation term:
  //
  // - 2 w'^2 dw/dz;
  //
  // and pressure term 2:
  //
  // - C_5 ( -2 w'^2 dw/dz + 2 (g/th_0) w'th_v' ).
  //
  // The w'^2 accumulation term is completely implicit, while w'^2 pressure
  // term 2 has both implicit and explicit components.  The accumulation term
  // and the implicit portion of pressure term 2 are combined and solved
  // together as:
  //
  // + ( 1 - C_uu_shr ) ( -2 w'^2(t+1) dw/dz ).
  //
  // Note 1:  When the term is brought over to the left-hand side, the sign
  //          is reversed and the leading "-" in front of the "2" is changed
  //          to a "+".
  // Note 2:  We have broken C5 up into C_uu_shr for this term
  //          and C_uu_buoy for the buoyancy term.
  //
  // The timestep index (t+1) means that the value of w'^2 being used is from
  // the next timestep, which is being advanced to in solving the d(w'^2)/dt
  // equation.
  //
  // The terms are discretized as follows:
  //
  // The values of w'^2 are found on the momentum levels, while the values of
  // wm_zt (mean vertical velocity on thermodynamic levels) are found on the
  // thermodynamic levels.  The vertical derivative of wm_zt is taken over the
  // intermediate (central) momentum level.  It is then multiplied by w'^2
  // (implicitly calculated at timestep (t+1)) and the coefficients to yield
  // the desired results.
  //
  // -------wm_zt--------------------------------------------- t(k)
  //
  // ===============d(wm_zt)/dz============wp2================ m(k)
  //
  // -------wm_zt--------------------------------------------- t(k-1)
  //
  // The vertical indices t(k), m(k), and t(k-1) correspond with altitudes
  // zt(k), zm(k), and zt(k-1), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )
  //
  // References:
  //-----------------------------------------------------------------------

  // ------------------------ Input Variables ------------------------
  // wm_zt: w wind component at thermodynamic levels    [m/s]
  // C_uu_shr: Model parameter C_uu_shr                       [-]
  // ------------------------ Output Variable ------------------------
  // lhs_ac_pr2_wp2: LHS coefficient of wp2 ac and pr2 terms [1/s]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;
  const auto gr_invrs_dzm = gr.invrs_dzm;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_terms_ac_pr2_lhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set lower boundary to 0
          lhs_ac_pr2_wp2(i - 1, gr_k_lb_zm - 1) = zero;
          // Set upper boundary to 0
          lhs_ac_pr2_wp2(i - 1, gr_k_ub_zm - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_terms_ac_pr2_lhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Momentum main diagonal: [ x wp2(k,<t+1>) ]
          lhs_ac_pr2_wp2(i - 1, k - 1) =
              (((+((one - C_uu_shr(i - 1))) * two) * gr_invrs_dzm(i - 1, k - 1)) *
               (wm_zt(i - 1, k - 1) - wm_zt(i - 1, (k - 1) - 1)));
        });
  }
  return;
}

// =============================================================================
void wp2_term_dp1_lhs(int nzm, int ngrdcol, const grid& gr,  // In
                      Real2D C1_Skw_fnc, Real2D invrs_tau1m, // In
                      Real2D lhs_dp1_wp2) {                  // Out
  //
  // Description:
  // Dissipation term 1 for w'^2:  implicit portion of the code.
  //
  // The d(w'^2)/dt equation contains dissipation term 1:
  //
  // - ( C_1 / tau_1m ) w'^2.
  //
  // Since w'^2 has a minimum threshold, the term should be damped only to that
  // threshold.  The term becomes:
  //
  // - ( C_1 / tau_1m ) * ( w'^2 - threshold ).
  //
  // This term is broken into implicit and explicit portions.  The implicit
  // portion of this term is:
  //
  // - ( C_1 / tau_1m ) w'^2(t+1).
  //
  // Note:  When the implicit term is brought over to the left-hand side, the
  //        sign is reversed and the leading "-" in front of the term is
  //        changed to a "+".
  //
  // The timestep index (t+1) means that the value of w'^2 being used is from
  // the next timestep, which is being advanced to in solving the d(w'^2)/dt
  // equation.
  //
  // The values of w'^2 are found on the momentum levels.  The values of the
  // C_1 skewness function and time-scale tau1m are also found on the momentum
  // levels.
  //
  // References:
  //-----------------------------------------------------------------------

  // ------------------ Input Variables ------------------
  // C1_Skw_fnc: C_1 parameter with Sk_w applied    [-]
  // invrs_tau1m: Inverse time-scale tau at momentum levels  [1/s]
  // ------------------ Output Variable ------------------
  // lhs_dp1_wp2: LHS coefficient of wp2 dissipation term 1  [1/s]
  // ------------------ Local Variable ------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_dp1_lhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set lower boundary to 0
          lhs_dp1_wp2(i - 1, gr_k_lb_zm - 1) = zero;
          // Set upper boundary to 0
          lhs_dp1_wp2(i - 1, gr_k_ub_zm - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_dp1_lhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Momentum main diagonal: [ x wp2(k,<t+1>) ]
          lhs_dp1_wp2(i - 1, k - 1) = (+(C1_Skw_fnc(i - 1, k - 1)) * invrs_tau1m(i - 1, k - 1));
        });
  }
  return;
}

// =============================================================================
void wp2_term_pr1_lhs(int nzm, int ngrdcol, const grid& gr, Real1D C4, Real2D invrs_tau_C4_zm, // In
                      Real2D lhs_pr1_wp2) { // Out
  //
  // Description
  // Pressure term 1 for w'^2:  implicit portion of the code.
  //
  // The d(w'^2)/dt equation contains pressure term 1:
  //
  // - ( C_4 / tau_1m ) * ( w'^2 - (2/3)*em ),
  //
  // where em = (1/2) * ( w'^2 + u'^2 + v'^2 ).
  //
  // This simplifies to:
  //
  // - ( C_4 / tau_1m ) * (2/3) * w'^2
  //    + ( C_4 / tau_1m ) * (1/3) * ( u'^2 + v'^2 ).
  //
  // Pressure term 1 has both implicit and explicit components.  The implicit
  // portion is:
  //
  // - ( C_4 / tau_1m ) * (2/3) * w'^2(t+1);
  //
  // and is computed in this function.
  //
  // Note:  When the implicit term is brought over to the left-hand side, the
  //        sign is reversed and the leading "-" in front of the term is
  //        changed to a "+".
  //
  // The timestep index (t+1) means that the value of w'^2 being used is from
  // the next timestep, which is being advanced to in solving the d(w'^2)/dt
  // equation.
  //
  // The values of w'^2 are found on momentum levels, as are the values of
  // tau1m.
  //
  // References:
  //-----------------------------------------------------------------------

  // --------------------- Input Variables ---------------------
  // invrs_tau_C4_zm: Inverse time-scale tau at momentum levels  [1/s]
  // C4: Model parameter C_4                [-]
  // --------------------- Output Variable ---------------------
  // lhs_pr1_wp2: LHS coefficient of wp2 pressure term 1  [1/s]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_pr1_lhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set lower boundary to 0
          lhs_pr1_wp2(i - 1, gr_k_lb_zm - 1) = zero;
          // Set upper boundary to 0
          lhs_pr1_wp2(i - 1, gr_k_ub_zm - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_pr1_lhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Momentum main diagonal: [ x wp2(k,<t+1>) ]
          lhs_pr1_wp2(i - 1, k - 1) =
              (+(((two * C4(i - 1)) * invrs_tau_C4_zm(i - 1, k - 1))) / three);
        });
  }
  return;
}

// =============================================================================
void wp2_terms_bp_pr2_rhs(int nzm, int ngrdcol, const grid& gr, // In
                          Real1D C_uu_buoy,                     // In
                          Real2D thv_ds_zm, Real2D wpthvp,      // In
                          Real2D rhs_bp_pr2_wp2) {              // Out
  //
  // Description:
  // Buoyancy production of w'^2 and w'^2 pressure term 2:  explicit portion of
  // the code.
  //
  // The d(w'^2)/dt equation contains a buoyancy production term:
  //
  // + 2 (g/thv_ds) w'th_v';
  //
  // and pressure term 2:
  //
  // - C_5 ( -2 w'^2 dw/dz + 2 (g/thv_ds) w'th_v' ).
  //
  // The w'^2 buoyancy production term is completely explicit, while w'^2
  // pressure term 2 has both implicit and explicit components.  The buoyancy
  // production term and the explicit portion of pressure term 2 are combined
  // and solved together as:
  //
  // + ( 1 - C_uu_buoy ) ( 2 (g/thv_ds) w'th_v' ).
  //
  // Note:  We have broken C5 up into C_uu_shr for the accumulation term
  //        and C_uu_buoy for the buoyancy term.
  //
  // References:
  //-----------------------------------------------------------------------

  // ------------------ Input Variables ------------------
  // C_uu_buoy: Model parameter C_uu_buoy                         [-]
  // thv_ds_zm: Dry, base-state theta_v at momentum levels   [K]
  // wpthvp: w'th_v'                                      [K m/s]
  // ------------------ Output Variable ------------------
  // rhs_bp_pr2_wp2: RHS portion of wp2 from terms bp and pr2  [m^2/s^3]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_terms_bp_pr2_rhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          rhs_bp_pr2_wp2(i - 1, gr_k_lb_zm - 1) = zero;
          // Set upper boundary to 0
          rhs_bp_pr2_wp2(i - 1, gr_k_ub_zm - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_terms_bp_pr2_rhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          rhs_bp_pr2_wp2(i - 1, k - 1) =
              (((+((one - C_uu_buoy(i - 1))) * two) * (grav / thv_ds_zm(i - 1, k - 1))) *
               wpthvp(i - 1, k - 1));
        });
  }
  return;
}

// =============================================================================
void wp2_term_dp1_rhs(int nzm, int ngrdcol, const grid& gr, Real2D C1_Skw_fnc,      // In
                      Real2D invrs_tau1m, real_t threshold, Real2D up2, Real2D vp2, // In
                      bool l_damp_wp2_using_em,                                     // In
                      Real2D rhs_dp1_wp2) {                                         // Out
  //
  // Description:
  // When l_damp_wp2_using_em == .false., then
  // Dissipation term 1 for w'^2:  explicit portion of the code.
  //
  // The d(w'^2)/dt equation contains dissipation term 1:
  //
  // - ( C_1 / tau_1m ) w'^2.
  //
  // Since w'^2 has a minimum threshold, the term should be damped only to that
  // threshold.  The term becomes:
  //
  // - ( C_1 / tau_1m ) * ( w'^2 - threshold ).
  //
  // This term is broken into implicit and explicit portions.  The explicit
  // portion of this term is:
  //
  // + ( C_1 / tau_1m ) * threshold.
  //
  // The values of the C_1 skewness function, time-scale tau1m, and the
  // threshold are found on the momentum levels.
  //
  // if l_damp_wp2_using_em == .true., then
  // we damp wp2 using a more standard turbulence closure, -(2/3)*em/tau
  // This only works if C1=C14 and l_stability_correct_tau_zm =.false.
  // A factor of (1/3) is absorbed into C1.
  // The threshold is implicitly set to 0.
  //
  //
  // References:
  //-----------------------------------------------------------------------

  // -------------------- Input Variables --------------------
  // C1_Skw_fnc: C_1 parameter with Sk_w applied                  [-]
  // invrs_tau1m: Inverse time-scale tau at momentum levels        [1/s]
  // up2: Horizontal (east-west) velocity variance, u'^2   [m^2/s^2]
  // vp2: Horizontal (north-south) velocity variance, v'^2 [m^2/s^2]
  // threshold: Minimum allowable value of w'^2       [m^2/s^2]
  // l_damp_wp2_using_em: In the wp2 equation, use a dissipation formula of -(2/3)*em/tau_zm,
  // as in Bougeault (1981)
  // -------------------- Output Variable --------------------
  // rhs_dp1_wp2: RHS portion of wp2 from dissipation term 1  [m^2/s^3]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_dp1_rhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          rhs_dp1_wp2(i - 1, gr_k_lb_zm - 1) = zero;
          // Set upper boundary to 0
          rhs_dp1_wp2(i - 1, gr_k_ub_zm - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  if (l_damp_wp2_using_em) {
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp2_term_dp1_rhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            rhs_dp1_wp2(i - 1, k - 1) = (-((C1_Skw_fnc(i - 1, k - 1) * invrs_tau1m(i - 1, k - 1))) *
                                         (up2(i - 1, k - 1) + vp2(i - 1, k - 1)));
          });
    }
  } else {
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp2_term_dp1_rhs::cells_3", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            rhs_dp1_wp2(i - 1, k - 1) =
                (+((C1_Skw_fnc(i - 1, k - 1) * invrs_tau1m(i - 1, k - 1))) * threshold);
          });
    }
  }
  // l_damp_wp2_using_em
  return;
}

// =============================================================================
void wp2_term_pr3_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                      Real1D C_uu_shr,                               // In
                      Real1D C_uu_buoy,                              // In
                      Real2D thv_ds_zm, Real2D wpthvp, Real2D upwp,  // In
                      Real2D um, Real2D vpwp, Real2D vm,             // In
                      Real2D rhs_pr3_wp2) {                          // Out
  //
  // Description:
  // Pressure term 3 for w'^2:  explicit portion of the code.
  //
  // The d(w'^2)/dt equation contains pressure term 3:
  //
  // + (2/3) C_5 [ (g/thv_ds) w'th_v' - u'w' du/dz - v'w' dv/dz ].
  //
  // Note that below we have broken up C5 into C_uu_shr for shear terms and
  // C_uu_buoy for buoyancy terms.
  //
  // This term is solved for completely explicitly and is discretized as
  // follows:
  //
  // The values of w'th_v', u'w', and v'w' are found on the momentum levels,
  // whereas the values of um and vm are found on the thermodynamic levels.
  // Additionally, the values of thv_ds_zm are found on the momentum levels.
  // The derivatives of both um and vm are taken over the intermediate
  // (central) momentum level.  All the remaining mathematical operations take
  // place at the central momentum level, yielding the desired result.
  //
  // -----um--------------vm---------------------------------------- t(k)
  //
  // =upwp====d(um)/dz========d(vm)/dz==vpwp===thv_ds_zm==wpthvp==== m(k)
  //
  // -----um--------------vm---------------------------------------- t(k-1)
  //
  // The vertical indices t(k), m(k), and t(k-1) correspond with altitudes
  // zt(k), zm(k), and zt(k-1), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )
  //
  // References:
  //-----------------------------------------------------------------------

  // --------------------- Input Variables ---------------------
  // C_uu_shr: Model parameter                            [-]
  // C_uu_buoy: Model parameter                            [-]
  // thv_ds_zm: Dry, base-state theta_v at momentum level (k)  [K]
  // wpthvp: w'th_v'(k)                                     [K m/s]
  // upwp: u'w'(k)                                        [m^2/s^2]
  // vpwp: v'w'(k)                                        [m^2/s^2]
  // um: um(k)                                          [m/s]
  // vm: vm(k)                                          [m/s]
  // --------------------- Output Variable ---------------------
  // rhs_pr3_wp2: RHS portion of wp2 from pressure term 3  [m^2/s^3]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;
  const auto gr_invrs_dzm = gr.invrs_dzm;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_pr3_rhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          rhs_pr3_wp2(i - 1, gr_k_lb_zm - 1) = zero;
          // Set upper boundary to 0
          rhs_pr3_wp2(i - 1, gr_k_ub_zm - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_pr3_rhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Michael Falk, 2 August 2007
          // Use the following code for standard mixing, with c_k=0.548:
          rhs_pr3_wp2(i - 1, k - 1) =
              (+(two_thirds) *
               (((C_uu_buoy(i - 1) * (grav / thv_ds_zm(i - 1, k - 1))) * wpthvp(i - 1, k - 1)) +
                (C_uu_shr(i - 1) * (((-(upwp(i - 1, k - 1)) * gr_invrs_dzm(i - 1, k - 1)) *
                                     (um(i - 1, k - 1) - um(i - 1, (k - 1) - 1))) -
                                    ((vpwp(i - 1, k - 1) * gr_invrs_dzm(i - 1, k - 1)) *
                                     (vm(i - 1, k - 1) - vm(i - 1, (k - 1) - 1)))))));
          // Use the following code for alternate mixing, with c_k=0.1 or 0.2
          // = + two_thirds * C_uu_shr &
          // * ( ( grav / thv_ds_zm(k) ) * wpthvp(k) &
          // - 0. * upwp(k) * invrs_dzm(k) * ( um(k) - um(k-1) ) &
          // - 0. * vpwp(k) * invrs_dzm(k) * ( vm(k) - vm(k-1) ) &
          // )
          // eMFc
          // Added by dschanen for ticket #36
          // We have found that when shear generation is zero this term will only be
          // offset by hole-filling (wp2_pd) and reduces turbulence
          // unrealistically at lower altitudes to make up the difference.
          rhs_pr3_wp2(i - 1, k - 1) =
              (rhs_pr3_wp2(i - 1, k - 1) > zero_threshold ? rhs_pr3_wp2(i - 1, k - 1)
                                                          : zero_threshold);
        });
  }
  return;
}

// =============================================================================
void wp2_term_pr1_rhs(int nzm, int ngrdcol, const grid& gr, Real1D C4, // In
                      Real2D up2, Real2D vp2, Real2D invrs_tau_C4_zm,  // In
                      Real2D rhs_pr1_wp2) {                            // Out
  //
  // Description:
  // Pressure term 1 for w'^2:  explicit portion of the code.
  //
  // The d(w'^2)/dt equation contains pressure term 1:
  //
  // - ( C_4 / tau_1m ) * ( w'^2 - (2/3)*em );
  //
  // where em = (1/2) * ( w'^2 + u'^2 + v'^2 ).
  //
  // This simplifies to:
  //
  // - ( C_4 / tau_1m ) * (2/3) * w'^2
  //    + ( C_4 / tau_1m ) * (1/3) * ( u'^2 + v'^2 ).
  //
  // Pressure term 1 has both implicit and explicit components.
  // The explicit portion is:
  //
  // + ( C_4 / tau_1m ) * (1/3) * ( u'^2 + v'^2 );
  //
  // and is computed in this function.
  //
  // The values of u'^2 and v'^2 are found on momentum levels, as are the
  // values of tau1m.
  //
  // References:
  //-----------------------------------------------------------------------

  // ------------------------ Input Variables ------------------------
  // C4: Model parameter C_4                      [-]
  // up2: u'^2(k)                               [m^2/s^2]
  // vp2: v'^2(k)                               [m^2/s^2]
  // invrs_tau_C4_zm: Inverse time-scale tau at momentum levels [1/s]
  // ------------------------ Output Variable ------------------------
  // rhs_pr1_wp2: RHS portion of wp2 from pressure term 1  [m^2/s^3]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;

  // Set lower bounadry to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_pr1_rhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          rhs_pr1_wp2(i - 1, gr_k_lb_zm - 1) = zero;
          // Set upper boundary to 0
          rhs_pr1_wp2(i - 1, gr_k_ub_zm - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_pr1_rhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          rhs_pr1_wp2(i - 1, k - 1) = (+(((C4(i - 1) * (up2(i - 1, k - 1) + vp2(i - 1, k - 1))) *
                                          invrs_tau_C4_zm(i - 1, k - 1))) /
                                       three);
        });
  }
  return;
}

// =============================================================================
void wp2_term_pr_dfsn_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real1D C_wp2_pr_dfsn, // In
                          Real2D rho_ds_zt, Real2D invrs_rho_ds_zm,                            // In
                          Real2D wpup2, Real2D wpvp2, Real2D wp3,                              // In
                          Real2D rhs_pr_dfsn_wp2) { // Out
  //
  // Description:
  //
  // This term is intended to represent the "diffusion" part of the wp2
  // pressure correlation.  The total pressure diffusion term,
  //
  //   -1 / rho * ( d( <u_k'p'> )/dx_i + d( <u_i'p'> )/dx_k )
  //
  // becomes
  //
  //   -2 / rho * d( <w'p'> )/dz
  //
  // for the w'^2 equation.  The factor of 2 is replaced with a tunable
  // parameter, C_wp2_pr_dfsn, and p' is replaced with
  //
  //   p' ~ - rho * ( u_i*u_i - <u_i*u_i> ),
  //
  // following Lumley 1978.  The wp2 pressure diffusion term becomes
  //
  //   + C_wp2_pr_dfsn / rho * ( d( rho*<w'u_iu_i> )/dz )
  //
  // References:
  //   Lumley 1978, p. 170.  See eq. 6.47 and accompanying discussion.
  //-----------------------------------------------------------------------

  // ---------------------- Input Variables ----------------------
  // C_wp2_pr_dfsn: Model parameter C_wp2_pr_dfsn                [-]
  // invrs_rho_ds_zm: Inverse dry static density (momentum levels) [m^3/kg]
  // rho_ds_zt: Dry static density on thermo. levels    [kg/m^3]
  // wpup2: w'u'^2 on thermodynamic levels          [m^3/s^3]
  // wpvp2: w'v'^2 on thermodynamic levels          [m^3/s^3]
  // wp3: w'^3 on thermo levels                   [m^3/s^3]
  // ---------------------- Output Variable ----------------------
  // rhs_pr_dfsn_wp2: RHS portion of wp2 from pressure-diffusion correlation [m^2/s^3]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D wpuip2("wp2_term_pr_dfsn_rhs::wpuip2", ngrdcol,
                nzm); // 3rd-order moment sum <w'u_i'u_i'>     [m^3/s^3]
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_invrs_dzm = gr.invrs_dzm;
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_ub_zm = gr.k_ub_zm;

  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzt) : (1) >= (nzt)) ? ((nzt) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_pr_dfsn_rhs::cells_1", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          wpuip2(i - 1, k - 1) = ((wpup2(i - 1, k - 1) + wpvp2(i - 1, k - 1)) + wp3(i - 1, k - 1));
        });
  }
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzm - 1)) : (2) >= ((nzm - 1))) ? (((nzm - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_pr_dfsn_rhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          rhs_pr_dfsn_wp2(i - 1, k - 1) =
              (((+(C_wp2_pr_dfsn(i - 1)) * invrs_rho_ds_zm(i - 1, k - 1)) *
                gr_invrs_dzm(i - 1, k - 1)) *
               ((rho_ds_zt(i - 1, k - 1) * wpuip2(i - 1, k - 1)) -
                (rho_ds_zt(i - 1, (k - 1) - 1) * wpuip2(i - 1, (k - 1) - 1))));
        });
  }
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp2_term_pr_dfsn_rhs::cells_3", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set lower boundary condition
          rhs_pr_dfsn_wp2(i - 1, gr_k_lb_zm - 1) =
              rhs_pr_dfsn_wp2(i - 1, (gr_k_lb_zm + gr_grid_dir_indx) - 1);
          // Set upper boundary to 0
          rhs_pr_dfsn_wp2(i - 1, gr_k_ub_zm - 1) = zero;
        });
  }
  return;
}

// =============================================================================
void wp3_term_ta_new_pdf_lhs(int nzm, int nzt, int ngrdcol, const grid& gr,        // In
                             Real2D coef_wp4_implicit,                             // In
                             Real2D wp2, Real2D rho_ds_zm, Real2D invrs_rho_ds_zt, // In
                             Real3D lhs_ta_wp3) {                                  // Out
  //
  // Description:
  // Turbulent advection of <w'^3>:  implicit portion of the code.
  //
  // This implicit discretization is specifically for the new PDF.
  //
  // The d<w'^3>/dt equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * <w'^4> )/dz.
  //
  // A substitution, which is specific to the new PDF, is made in order to
  // close the turbulent advection term, such that:
  //
  // <w'^4> = coef_wp4_implicit * <w'^2>^2.
  //
  // The calculation of coef_wp4_implicit is detailed in function
  // calc_coef_wp4_implicit, which is found in module new_pdf in new_pdf.F90.
  //
  // The turbulent advection term is rewritten as:
  //
  // - (1/rho_ds) * d( rho_ds * coef_wp4_implicit * <w'^2>^2 )/dz.
  //
  // The <w'^2>^2 term is timestep split so that it can be expressed linearly
  // in terms of <w'^2> at the (t+1) timestep, such that:
  //
  // <w'^2>^2 = <w'^2>(t) * <w'^2>(t+1);
  //
  // which allows the turbulent advection term to be expressed implicitly as:
  //
  // - (1/rho_ds)
  //   * d( rho_ds * coef_wp4_implicit * <w'^2>(t) * <w'^2>(t+1) )/dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign is
  //        reversed and the leading "-" in front of all d[ ] / dz terms is
  //        changed to a "+".
  //
  // Timestep index (t) stands for the index of the current timestep, while
  // timestep index (t+1) stands for the index of the next timestep, which is
  // being advanced to in solving the d<w'^3>/dt and d<w'^2>/dt equations.
  //
  // The implicit discretization of this term is as follows:
  //
  // The values of <w'^3> are found on the thermodynamic levels, while the
  // values of <w'^2> are found on the momentum levels.  The values of
  // coef_wp4_implicit_zt are originally calculated by the PDF on the
  // thermodynamic levels.  They are interpolated to the intermediate momentum
  // levels as coef_wp4_implicit.  Additionally, the values of rho_ds_zm are
  // found on the momentum levels, and the values of invrs_rho_ds_zt are found
  // on the thermodynamic levels.  At the intermediate momentum levels, the
  // values of coef_wp4_implicit are multiplied by <w'^2>(t) * <w'^2>(t+1), and
  // the resulting product is also multiplied by rho_ds_zm.  This product is
  // referred to as G below.  Then, the derivative (d/dz) of that expression is
  // taken over the central thermodynamic level, where it is multiplied by
  // -invrs_rho_ds_zt.  This yields the desired result.  In this function,
  // the values of G are as follows:
  //
  // G = rho_ds_zm * coef_wp4_implicit * <w'^2>(t) * <w'^2>(t+1).
  //
  // -------coef_wp4_implicit_zt---------------------------------------- t(k+1)
  //
  // =======coef_wp4_implicit(interp)=======wp2=========rho_ds_zm======= m(k+1)
  //
  // -------coef_wp4_implicit_zt-----dG/dz-----invrs_rho_ds_zt----wp3--- t(k)
  //
  // =======coef_wp4_implicit(interp)=======wp2=========rho_ds_zm======= m(k)
  //
  // -------coef_wp4_implicit_zt---------------------------------------- t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is
  // used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
  //
  // References:
  //-----------------------------------------------------------------------

  // Constant parameters
  // ------------------------ Input Variables ------------------------
  // coef_wp4_implicit: <w'^4>=coef_wp4_implicit*<w'^2>^2; m-levs [-]
  // wp2: <w'^2>                                    [m^2/s^2]
  // rho_ds_zm: Dry, static density at momentum levels    [kg/m^3]
  // invrs_rho_ds_zt: Inv dry, static density at thermo levels  [m^3/kg]
  // ------------------------ Output Variable ------------------------
  // lhs_ta_wp3: LHS coefficient of wp3 turbulent advection  [m/s^2]
  // ------------------------ Local Variable ------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  constexpr int k_mdiag = 1;   // Momentum superdiagonal index.
  constexpr int km1_mdiag = 2; // Momentum subdiagonal index.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_invrs_dzt = gr.invrs_dzt;

  // Set term at lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_ta_new_pdf_lhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set term at lower boundary to 0
          // Lower the Native whole-band assignment without creating a device slice.
          for (int band = 0; band < ndiags2; ++band)
            lhs_ta_wp3(band, i - 1, gr_k_lb_zt - 1) = zero;
          // Set term at upper boundary to 0
          // Lower the Native whole-band assignment without creating a device slice.
          for (int band = 0; band < ndiags2; ++band)
            lhs_ta_wp3(band, i - 1, gr_k_ub_zt - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_ta_new_pdf_lhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Momentum superdiagonal: [ x wp2(k+1,<t+1>) ]
          lhs_ta_wp3(k_mdiag - 1, i - 1, k - 1) =
              ((((+(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                 rho_ds_zm(i - 1, (k + 1) - 1)) *
                coef_wp4_implicit(i - 1, (k + 1) - 1)) *
               wp2(i - 1, (k + 1) - 1));
          // Momentum subdiagonal: [ x wp2(k,<t+1>) ]
          lhs_ta_wp3(km1_mdiag - 1, i - 1, k - 1) =
              ((((-(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                 rho_ds_zm(i - 1, k - 1)) *
                coef_wp4_implicit(i - 1, k - 1)) *
               wp2(i - 1, k - 1));
        });
  }
  return;
}

// =============================================================================
void wp3_term_ta_ADG1_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, // In
                          Real2D wp2, Real2D a1_coef, Real2D a1_coef_zt, // In
                          Real2D a3_coef, Real2D a3_coef_zt,             // In
                          Real2D wp3_on_wp2, Real2D rho_ds_zm,           // In
                          Real2D rho_ds_zt, Real2D invrs_rho_ds_zt,      // In
                          bool l_standard_term_ta,                       // In
                          bool l_partial_upwind_wp3,                     // In
                          Real3D lhs_ta_wp3) {                           // Out
  //
  // Description:
  // Turbulent advection of w'^3:  implicit portion of the code.
  //
  // This implicit discretization is specifically for the ADG1 PDF.
  //
  // The d(w'^3)/dt equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * w'^4 )/dz.
  //
  // A substitution, which is specific to ADG1, is made in order to close the
  // turbulent advection term, such that:
  //
  // w'^4 = a_3 * (w'^2)^2  +  a_1 * ( (w'^3)^2 / w'^2 );
  //
  // where both a_1 and a_3 are variables that are functions of sigma_sqd_w,
  // such that:
  //
  // a_1 = 1 / (1 - sigma_sqd_w); and
  //
  // a_3 = 3*(sigma_sqd_w)^2 + 6*(1 - sigma_sqd_w)*sigma_sqd_w
  //       + (1 - sigma_sqd_w)^2.
  //
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:wp4_diagnosis
  //
  // The turbulent advection term is rewritten as:
  //
  // - (1/rho_ds) * d [ rho_ds * a_3 * (w'^2)^2 ] / dz
  // - (1/rho_ds) * d [ rho_ds * a_1 * ( (w'^3)^2 / w'^2 ) ] / dz.
  //
  // The (w'^2)^2 and (w'^3)^2 terms are both timestep split so that they can
  // be expressed linearly in terms of w'^2 and w'^3, respectively, at the
  // (t+1) timestep, such that:
  //
  // (w'^2)^2 = w'^2(t) * w'^2(t+1);
  // (w'^3)^2 = w'^3(t) * w'^3(t+1);
  //
  // which allows these terms to be expressed implicitly as:
  //
  // - (1/rho_ds) * d [ rho_ds * a_3 * w'^2(t) * w'^2(t+1) ] / dz
  // - (1/rho_ds) * d [ rho_ds * a_1 * w'^3(t) * w'^3(t+1) / w'^2(t) ] / dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign is
  //        reversed and the leading "-" in front of all d[ ] / dz terms is
  //        changed to a "+".
  //
  // Timestep index (t) stands for the index of the current timestep, while
  // timestep index (t+1) stands for the index of the next timestep, which is
  // being advanced to in solving the d(w'^3)/dt and d(w'^2)/dt equations.
  //
  // The implicit portion of these terms is discretized as follows:
  //
  // The values of w'^3 are found on the thermodynamic levels, while the values
  // of w'^2, a_1, and a_3 are found on the momentum levels.  Additionally, the
  // values of rho_ds_zm are found on the momentum levels, and the values of
  // invrs_rho_ds_zt are found on the thermodynamic levels.  The variable w'^3
  // is interpolated to the intermediate momentum levels.  The values of the
  // mathematical expressions (called F and G here) within the dF/dz and dG/dz
  // terms are computed on the momentum levels.  Then, the derivatives (d/dz)
  // of the expressions (F and G) are taken over the central thermodynamic
  // level, where dF/dz and dG/dz are multiplied by -invrs_rho_ds_zt.  This
  // yields the desired results.  In this function, the values of F and G are
  // as follows:
  //
  // F = rho_ds_zm * a_3(t) * w'^2(t) * w'^2(t+1); and
  //
  // G = rho_ds_zm * a_1(t) * w'^3(t) * w'^3(t+1) / w'^2(t).
  //
  // ------------------------------------------------wp3---------------- t(k+1)
  //
  // ===a3====wp2====rho_ds_zm====a1======================wp3(interp)=== m(k+1)
  //
  // -----------dF/dz----invrs_rho_ds_zt----dG/dz----wp3---------------- t(k)
  //
  // ===a3====wp2====rho_ds_zm====a1======================wp3(interp)=== m(k)
  //
  // ------------------------------------------------wp3---------------- t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is
  // used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
  //
  // References:
  //-----------------------------------------------------------------------

  // Constant parameters
  // ---------------------- Input Variables ----------------------
  // wp2: w'^2                                     [m^2/s^2]
  // a1_coef: a_1                                      [-]
  // a3_coef: a_3                                      [-]
  // wp3_on_wp2: w'^3 / w'^2 at momentum levels           [m/s]
  // rho_ds_zm: Dry, static density at momentum levels   [kg/m^3]
  // a1_coef_zt: a_1 interpolated to thermodynamic levels [-]
  // a3_coef_zt: a_3 interpolated to thermodynamic levels [-]
  // rho_ds_zt: Dry, static density at thermo. levels    [kg/m^3]
  // invrs_rho_ds_zt: Inv dry, static density at thermo levels [m^3/kg]
  // l_standard_term_ta: Use the standard discretization for the
  // turbulent advection terms.  Setting to .false.
  // means that a_1 and a_3 are pulled outside of the
  // derivative in advance_wp2_wp3_module.F90 and in
  // advance_xp2_xpyp_module.F90.
  // l_partial_upwind_wp3: Flag to use an "upwind" discretization rather
  // than a centered discretization for the portion
  // of the wp3 turbulent advection term for ADG1
  // that is linearized in terms of wp3<t+1>.
  // (Requires ADG1 PDF and l_standard_term_ta).
  // ---------------------- Output Variable ----------------------
  // lhs_ta_wp3: LHS coefficient of wp3 turbulent advection

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  constexpr int kp1_tdiag = 1; // Thermodynamic superdiagonal index.
  constexpr int k_mdiag = 2;   // Momentum superdiagonal index.
  constexpr int k_tdiag = 3;   // Thermodynamic main diagonal index.
  constexpr int km1_mdiag = 4; // Momentum subdiagonal index.
  constexpr int km1_tdiag = 5; // Thermodynamic subdiagonal index.
  constexpr int t_above = 1;   // Index for upper thermodynamic level grid weight.
  constexpr int t_below = 2;   // Index for lower thermodynamic level grid weight.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_invrs_dzt = gr.invrs_dzt;
  const auto gr_weights_zt2zm = gr.weights_zt2zm;
  const auto gr_grid_dir = gr.grid_dir;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    const int extent_1 = ((1) > 0 ? (1) <= (5) : (1) >= (5)) ? ((5) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_ta_ADG1_lhs::cells_1", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Set lower boundary to 0
          lhs_ta_wp3(i - 1, k - 1, gr_k_lb_zt - 1) = zero;
          // Set upper boundary to 0
          lhs_ta_wp3(i - 1, k - 1, gr_k_ub_zt - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  if (l_standard_term_ta) {
    // The turbulent advection term is discretized normally, in accordance
    // with the model equations found in the documentation and the description
    // listed above.
    if (!(l_partial_upwind_wp3)) {
      // All portions of the wp3 turbulent advection term for ADG1 use
      // centered discretization in accordance with description and diagram
      // shown above.
      {
        const int extent_0 =
            ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
        const int extent_1 =
            ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
        Kokkos::parallel_for(
            "wp3_term_ta_ADG1_lhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
              const int i = (1) + ((cell / (1)) % extent_1) * (1);
              const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
              // Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
              lhs_ta_wp3(kp1_tdiag - 1, i - 1, k - 1) =
                  (((((+(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                      rho_ds_zm(i - 1, (k + 1) - 1)) *
                     a1_coef(i - 1, (k + 1) - 1)) *
                    wp3_on_wp2(i - 1, (k + 1) - 1)) *
                   gr_weights_zt2zm(i - 1, (k + 1) - 1, t_above - 1));
              // Momentum superdiagonal: [ x wp2(k+1,<t+1>) ]
              lhs_ta_wp3(k_mdiag - 1, i - 1, k - 1) =
                  ((((+(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                     rho_ds_zm(i - 1, (k + 1) - 1)) *
                    a3_coef(i - 1, (k + 1) - 1)) *
                   wp2(i - 1, (k + 1) - 1));
              // Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
              lhs_ta_wp3(k_tdiag - 1, i - 1, k - 1) =
                  ((+(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                   ((((rho_ds_zm(i - 1, (k + 1) - 1) * a1_coef(i - 1, (k + 1) - 1)) *
                      wp3_on_wp2(i - 1, (k + 1) - 1)) *
                     gr_weights_zt2zm(i - 1, (k + 1) - 1, t_below - 1)) -
                    (((rho_ds_zm(i - 1, k - 1) * a1_coef(i - 1, k - 1)) *
                      wp3_on_wp2(i - 1, k - 1)) *
                     gr_weights_zt2zm(i - 1, k - 1, t_above - 1))));
              // Momentum subdiagonal: [ x wp2(k,<t+1>) ]
              lhs_ta_wp3(km1_mdiag - 1, i - 1, k - 1) =
                  ((((-(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                     rho_ds_zm(i - 1, k - 1)) *
                    a3_coef(i - 1, k - 1)) *
                   wp2(i - 1, k - 1));
              // Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
              lhs_ta_wp3(km1_tdiag - 1, i - 1, k - 1) =
                  (((((-(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                      rho_ds_zm(i - 1, k - 1)) *
                     a1_coef(i - 1, k - 1)) *
                    wp3_on_wp2(i - 1, k - 1)) *
                   gr_weights_zt2zm(i - 1, k - 1, t_below - 1));
            });
      }
    } else {
      // l_partial_upwind_wp3
      // Partial upwinding of the wp3 turbulent advection term, where the
      // portion of the wp3 turbulent advection term that is linearized in
      // terms of wp2<t+1> is still handled using centered discretization,
      // but the portion of the term that is linearized in terms of wp3<t+1>
      // is handled using an "upwind" discretization that also takes into
      // "winds" that converge or diverge around the central thermodynamic
      // grid level.  Provided by Chris Vogl and Shixuan Zhang.
      {
        const int extent_0 =
            ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
        const int extent_1 =
            ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
        Kokkos::parallel_for(
            "wp3_term_ta_ADG1_lhs::cells_3", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
              const int i = (1) + ((cell / (1)) % extent_1) * (1);
              const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
              // Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
              lhs_ta_wp3(kp1_tdiag - 1, i - 1, k - 1) =
                  ((((+(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                     rho_ds_zt(i - 1, (k + 1) - 1)) *
                    gr_grid_dir) *
                   (((gr_grid_dir * a1_coef(i - 1, (k + 1) - 1)) * wp3_on_wp2(i - 1, (k + 1) - 1)) <
                            zero
                        ? ((gr_grid_dir * a1_coef(i - 1, (k + 1) - 1)) *
                           wp3_on_wp2(i - 1, (k + 1) - 1))
                        : zero));
              // Momentum superdiagonal: [ x wp2(k+1,<t+1>) ]
              lhs_ta_wp3(k_mdiag - 1, i - 1, k - 1) =
                  ((((+(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                     rho_ds_zm(i - 1, (k + 1) - 1)) *
                    a3_coef(i - 1, (k + 1) - 1)) *
                   wp2(i - 1, (k + 1) - 1));
              // Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
              lhs_ta_wp3(k_tdiag - 1, i - 1, k - 1) =
                  ((((+(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                     rho_ds_zt(i - 1, k - 1)) *
                    gr_grid_dir) *
                   ((((gr_grid_dir * a1_coef(i - 1, (k + 1) - 1)) *
                      wp3_on_wp2(i - 1, (k + 1) - 1)) > zero
                         ? ((gr_grid_dir * a1_coef(i - 1, (k + 1) - 1)) *
                            wp3_on_wp2(i - 1, (k + 1) - 1))
                         : zero) -
                    (((gr_grid_dir * a1_coef(i - 1, k - 1)) * wp3_on_wp2(i - 1, k - 1)) < zero
                         ? ((gr_grid_dir * a1_coef(i - 1, k - 1)) * wp3_on_wp2(i - 1, k - 1))
                         : zero)));
              // Momentum subdiagonal: [ x wp2(k,<t+1>) ]
              lhs_ta_wp3(km1_mdiag - 1, i - 1, k - 1) =
                  ((((-(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                     rho_ds_zm(i - 1, k - 1)) *
                    a3_coef(i - 1, k - 1)) *
                   wp2(i - 1, k - 1));
              // Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
              lhs_ta_wp3(km1_tdiag - 1, i - 1, k - 1) =
                  ((((-(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                     rho_ds_zt(i - 1, (k - 1) - 1)) *
                    gr_grid_dir) *
                   (((gr_grid_dir * a1_coef(i - 1, k - 1)) * wp3_on_wp2(i - 1, k - 1)) > zero
                        ? ((gr_grid_dir * a1_coef(i - 1, k - 1)) * wp3_on_wp2(i - 1, k - 1))
                        : zero));
            });
      }
    }
    // .not. l_partial_upwind_wp3
  } else {
    // Alternate discretization for the turbulent advection term, which
    // contains the term:
    // - (1/rho_ds) * d [ rho_ds * a_1 * (w'^3)^2 / w'^2 ] / dz.  In order
    // to help stabilize w'^3, a_1 has been pulled outside of the derivative.
    // On the left-hand side of the equation, this effects the thermodynamic
    // superdiagonal (kp1_tdiag), the thermodynamic main diagonal (k_tdiag),
    // and the thermodynamic subdiagonal (km1_tdiag).
    // Additionally, the discretization of the turbulent advection term, which
    // contains the term:
    // - (1/rho_ds) * d [ rho_ds * a_3 * (w'^2)^2 ] / dz, has been altered to
    // pull a_3 outside of the derivative.  This was done in order to help
    // stabilize w'^3.  On the left-hand side of the equation, this effects
    // the momentum superdiagonal (k_mdiag) and the momentum subdiagonal
    // (km1_mdiag).
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp3_term_ta_ADG1_lhs::cells_4", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            // Thermodynamic superdiagonal: [ x wp3(k+1,<t+1>) ]
            lhs_ta_wp3(kp1_tdiag - 1, i - 1, k - 1) =
                (((((+(invrs_rho_ds_zt(i - 1, k - 1)) * a1_coef_zt(i - 1, k - 1)) *
                    gr_invrs_dzt(i - 1, k - 1)) *
                   rho_ds_zm(i - 1, (k + 1) - 1)) *
                  wp3_on_wp2(i - 1, (k + 1) - 1)) *
                 gr_weights_zt2zm(i - 1, (k + 1) - 1, t_above - 1));
            // Momentum superdiagonal: [ x wp2(k+1,<t+1>) ]
            lhs_ta_wp3(k_mdiag - 1, i - 1, k - 1) =
                ((((+(invrs_rho_ds_zt(i - 1, k - 1)) * a3_coef_zt(i - 1, k - 1)) *
                   gr_invrs_dzt(i - 1, k - 1)) *
                  rho_ds_zm(i - 1, (k + 1) - 1)) *
                 wp2(i - 1, (k + 1) - 1));
            // Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
            lhs_ta_wp3(k_tdiag - 1, i - 1, k - 1) =
                (((+(invrs_rho_ds_zt(i - 1, k - 1)) * a1_coef_zt(i - 1, k - 1)) *
                  gr_invrs_dzt(i - 1, k - 1)) *
                 (((rho_ds_zm(i - 1, (k + 1) - 1) * wp3_on_wp2(i - 1, (k + 1) - 1)) *
                   gr_weights_zt2zm(i - 1, (k + 1) - 1, t_below - 1)) -
                  ((rho_ds_zm(i - 1, k - 1) * wp3_on_wp2(i - 1, k - 1)) *
                   gr_weights_zt2zm(i - 1, k - 1, t_above - 1))));
            // Momentum subdiagonal: [ x wp2(k,<t+1>) ]
            lhs_ta_wp3(km1_mdiag - 1, i - 1, k - 1) =
                ((((-(invrs_rho_ds_zt(i - 1, k - 1)) * a3_coef_zt(i - 1, k - 1)) *
                   gr_invrs_dzt(i - 1, k - 1)) *
                  rho_ds_zm(i - 1, k - 1)) *
                 wp2(i - 1, k - 1));
            // Thermodynamic subdiagonal: [ x wp3(k-1,<t+1>) ]
            lhs_ta_wp3(km1_tdiag - 1, i - 1, k - 1) =
                (((((-(invrs_rho_ds_zt(i - 1, k - 1)) * a1_coef_zt(i - 1, k - 1)) *
                    gr_invrs_dzt(i - 1, k - 1)) *
                   rho_ds_zm(i - 1, k - 1)) *
                  wp3_on_wp2(i - 1, k - 1)) *
                 gr_weights_zt2zm(i - 1, k - 1, t_below - 1));
          });
    }
  }
  // l_standard_term_ta
  return;
}

// =============================================================================
void wp3_term_tp_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real1D coef_wp3_tp, // In
                     Real2D wp2, Real2D rho_ds_zm, Real2D invrs_rho_ds_zt,              // In
                     Real3D lhs_tp_wp3) {                                               // Out
  //
  // Description:
  // Turbulent production of w'^3:  implicit portion of the code.
  //
  // The d(w'^3)/dt equation contains a turbulent production term:
  //
  // + 3 * ( w'^2 / rho_ds ) * d( rho_ds * w'^2 )/dz.
  //
  // The turbulent production term is rewritten as:
  //
  // + 3 * ( w'^2 / rho_ds ) * d[ rho_ds * w'^2 ]/dz
  // = + (3/rho_ds) * d[ rho_ds * (w'^2)^2 ]/dz - (3/2) * d[ (w'^2)^2 ]/dz.
  //
  // The (w'^2)^2 terms are timestep split so that they can be expressed
  // linearly in terms of w'^2 at the (t+1) timestep, such that:
  //
  // (w'^2)^2 = w'^2(t) * w'^2(t+1).
  //
  // The term can now be expressed implicitly as:
  //
  // + (3/rho_ds) * d [ rho_ds * w'^2(t) * w'^2(t+1) ] / dz
  // - (3/2) * d [ w'^2(t) * w'^2(t+1) ] /dz.
  //
  // Note:  When the term is brought over to the left-hand side, the sign is
  //        reversed and the leading "-" in front of a d[ ] / dz term is
  //        changed to a "+".  Likewise, the leading "+" in front of a
  //        d[ ] / dz term is changed to a "-".
  //
  // Timestep index (t) stands for the index of the current timestep, while
  // timestep index (t+1) stands for the index of the next timestep, which is
  // being advanced to in solving the d(w'^3)/dt and d(w'^2)/dt equations.
  //
  // The implicit portion of these terms is discretized as follows:
  //
  // While the values of w'^3 are found on the thermodynamic levels, the values
  // of w'^2 are found on the momentum levels.  Additionally, the values of
  // rho_ds_zm are found on the momentum levels, and the values of
  // invrs_rho_ds_zt are found on the thermodynamic levels.  The values of the
  // mathematical expressions (called F and G below) within the dF/dz and dG/dz
  // terms are computed on the momentum levels.  Then, the derivatives (d/dz)
  // of the expressions (F and G) are taken over the central thermodynamic
  // level, where dF/dz and dG/dz are multiplied by -3 * invrs_rho_ds_zt and
  // 3/2, respectively, yielding the desired results.  In this function, the
  // values of F and G are as follows:
  //
  // F = rho_ds_zm * w'^2(t) * w'^2(t+1);
  //
  // G = w'^2(t) * w'^2(t+1).
  //
  // ====wp2=========rho_ds_zm========================================== m(k+1)
  //
  // -----------dF/dz----invrs_rho_ds_zt----dG/dz----wp3---------------- t(k)
  //
  // ====wp2=========rho_ds_zm========================================== m(k)
  //
  // The vertical indices m(k+1), t(k), and m(k) correspond with altitudes
  // zm(k+1), zt(k), and zm(k), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
  //
  // References:
  //-----------------------------------------------------------------------

  // Constant parameters
  // -------------------- Input Variables --------------------
  // coef_wp3_tp: Coefficient for tp pressure scrambling term   [-]
  // wp2: w'^2                                        [m^2/s^2]
  // rho_ds_zm: Dry, static density at momentum levels      [kg/m^3]
  // invrs_rho_ds_zt: Inv dry, static density at thermo levels    [m^3/kg]
  // -------------------- Output Variable --------------------
  // lhs_tp_wp3: LHS coefficient multiplying wp2 in the wp3 equation [m/s^2]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  constexpr int k_mdiag = 1;   // Momentum superdiagonal index.
  constexpr int km1_mdiag = 2; // Momentum subdiagonal index.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_invrs_dzt = gr.invrs_dzt;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_tp_lhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int k = (1) + ((cell / (1)) % extent_0) * (1);
          int i{};
          for (i = 1; i <= (2); i += 1) {
            // Set lower boundary to 0
            lhs_tp_wp3(i - 1, k - 1, gr_k_lb_zt - 1) = zero;
            // Set upper boundary to 0
            lhs_tp_wp3(i - 1, k - 1, gr_k_ub_zt - 1) = zero;
          }
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_tp_lhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Momentum superdiagonal: [ x wp2(k+1,<t+1>) ]
          lhs_tp_wp3(k_mdiag - 1, i - 1, k - 1) =
              ((((((-(coef_wp3_tp(i - 1)) * three) * invrs_rho_ds_zt(i - 1, k - 1)) *
                  gr_invrs_dzt(i - 1, k - 1)) *
                 rho_ds_zm(i - 1, (k + 1) - 1)) *
                wp2(i - 1, (k + 1) - 1)) +
               (((coef_wp3_tp(i - 1) * three_halves) * gr_invrs_dzt(i - 1, k - 1)) *
                wp2(i - 1, (k + 1) - 1)));
          // Momentum subdiagonal: [ x wp2(k,<t+1>) ]
          lhs_tp_wp3(km1_mdiag - 1, i - 1, k - 1) =
              ((((((+(coef_wp3_tp(i - 1)) * three) * invrs_rho_ds_zt(i - 1, k - 1)) *
                  gr_invrs_dzt(i - 1, k - 1)) *
                 rho_ds_zm(i - 1, k - 1)) *
                wp2(i - 1, k - 1)) -
               (((coef_wp3_tp(i - 1) * three_halves) * gr_invrs_dzt(i - 1, k - 1)) *
                wp2(i - 1, k - 1)));
        });
  }
  return;
}

// =============================================================================
void wp3_terms_ac_pr2_lhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real2D C11_Skw_fnc, // In
                          Real2D wm_zm,                                                      // In
                          Real2D lhs_ac_pr2_wp3) {                                           // Out
  //
  // Description:
  // Accumulation of w'^3 and w'^3 pressure term 2:  implicit portion of the
  // code.
  //
  // The d(w'^3)/dt equation contains an accumulation term:
  //
  // - 3 w'^3 dw/dz;
  //
  // and pressure term 2:
  //
  // - C_11 ( -3 w'^3 dw/dz + 3 (g/th_0) w'^2th_v' ).
  //
  // The w'^3 accumulation term is completely implicit, while w'^3 pressure
  // term 2 has both implicit and explicit components.  The accumulation term
  // and the implicit portion of pressure term 2 are combined and solved
  // together as:
  //
  // + ( 1 - C_11 ) ( -3 w'^3(t+1) dw/dz ).
  //
  // Note:  When the term is brought over to the left-hand side, the sign
  //        is reversed and the leading "-" in front of the "3" is changed
  //        to a "+".
  //
  // The timestep index (t+1) means that the value of w'^3 being used is from
  // the next timestep, which is being advanced to in solving the d(w'^3)/dt
  // equation.
  //
  // The terms are discretized as follows:
  //
  // The values of w'^3 are found on thermodynamic levels, while the values of
  // wm_zm (mean vertical velocity on momentum levels) are found on momentum
  // levels.  The vertical derivative of wm_zm is taken over the intermediate
  // (central) thermodynamic level.  It is then multiplied by w'^3 (implicitly
  // calculated at timestep (t+1)) and the coefficients to yield the desired
  // results.
  //
  // =======wm_zm============================================= m(k+1)
  //
  // ---------------d(wm_zm)/dz------------wp3---------------- t(k)
  //
  // =======wm_zm============================================= m(k)
  //
  // The vertical indices m(k+1), t(k), and m(k) correspond with altitudes
  // zm(k+1), zt(k), and zm(k), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
  //
  // References:
  //-----------------------------------------------------------------------

  // ------------------------ Input Variables ------------------------
  // C11_Skw_fnc: C_11 parameter with Sk_w applied       [-]
  // wm_zm: w wind component at momentum levels    [m/s]
  // ------------------------ Output Variable ------------------------
  // lhs_ac_pr2_wp3: LHS coefficient of wp3 from terms ac and pr2 [1/s]
  // ------------------------ Local variable ------------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_invrs_dzt = gr.invrs_dzt;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_terms_ac_pr2_lhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set lower boundary to 0
          lhs_ac_pr2_wp3(i - 1, gr_k_lb_zt - 1) = zero;
          // Set upper boundary to 0
          lhs_ac_pr2_wp3(i - 1, gr_k_ub_zt - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_terms_ac_pr2_lhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          // Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
          lhs_ac_pr2_wp3(i - 1, k - 1) =
              (((+((one - C11_Skw_fnc(i - 1, k - 1))) * three) * gr_invrs_dzt(i - 1, k - 1)) *
               (wm_zm(i - 1, (k + 1) - 1) - wm_zm(i - 1, k - 1)));
        });
  }
  return;
}

// =============================================================================
void wp3_term_pr1_lhs(int nzt, int ngrdcol, const grid& gr,   // In
                      Real1D C8, Real1D C8b,                  // In
                      Real2D invrs_tau_wp3_zt, Real2D Skw_zt, // In
                      bool l_damp_wp3_Skw_squared,            // In
                      Real2D lhs_pr1_wp3) {                   // Out
  //
  // Description:
  // Pressure term 1 for w'^3:  implicit portion of the code.
  //
  // Pressure term 1 is the term:
  //
  // - (C_8/tau_w3t) * ( C_8b * Sk_wt^2 + 1 ) * w'^3;
  //
  // where Sk_wt = w'^3 / (w'^2)^(3/2).
  //
  // This term needs to be linearized, so function L(w'^3) is defined to be
  // equal to this term (pressure term 1), such that:
  //
  // L(w'^3) = - (C_8/tau_w3t) * ( C_8b * (w'^3)^3 / (w'^2)^3 + w'^3 ).
  //
  // A Taylor Series expansion (truncated after the first derivative term) of
  // L(w'^3) around w'^3 = w'^3(t) is used to linearize pressure term 1.
  // Evaluating L(w'^3) at w'^3(t+1):
  //
  // L( w'^3(t+1) ) = L( w'^3(t) )
  //                  + ( d L(w'^3) / d w'^3 )|_(w'^3=w'^3(t))
  //                    * ( w'^3(t+1) - w'^3(t) ).
  //
  // After evaluating the expression above, the term has become linearized.  It
  // is broken down into implicit (LHS) and explicit (RHS) components.
  // The implicit portion is:
  //
  // - (C_8/tau_w3t) * ( 3 * C_8b * Sk_wt^2 + 1 ) * w'^3(t+1).
  //
  // Note:  When the term is brought over to the left-hand side, the sign
  //        is reversed and the leading "-" in front of the term is changed
  //        to a "+".
  //
  // Timestep index (t) stands for the index of the current timestep, while
  // timestep index (t+1) stands for the index of the next timestep, which is
  // being advanced to in solving the d(w'^3)/dt equation.
  //
  // The values of w'^3 are found on the thermodynamic levels, as are the
  // values of tau_w3t and Sk_wt (in Sk_wt, w'^3 is found on thermodynamic
  // levels and w'^2 is interpolated to thermodynamic levels).
  //
  // References:
  //-----------------------------------------------------------------------

  // ---------------------- Input Variables ----------------------
  // invrs_tau_wp3_zt: Inverse time-scale tau at thermodynamic levels  [1/s]
  // Skw_zt: Skewness of w at thermodynamic levels   [-]
  // C8: Model parameter C_8                     [-]
  // C8b: Model parameter C_8b                    [-]
  // l_damp_wp3_Skw_squared: Set damping on wp3 to use Skw^2 rather than Skw^4
  // ---------------------- Output Variable ----------------------
  // lhs_pr1_wp3: LHS coefficient of wp3 from pressure term 1  [1/s]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_pr1_lhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set lower boundary to 0
          lhs_pr1_wp3(i - 1, gr_k_lb_zt - 1) = zero;
          // Set upper boundary to 0
          lhs_pr1_wp3(i - 1, gr_k_ub_zt - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  if (l_damp_wp3_Skw_squared) {
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp3_term_pr1_lhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            // Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
            lhs_pr1_wp3(i - 1, k - 1) =
                (+((C8(i - 1) * invrs_tau_wp3_zt(i - 1, k - 1))) *
                 (((three * C8b(i - 1)) * (Skw_zt(i - 1, k - 1) * Skw_zt(i - 1, k - 1))) + one));
          });
    }
  } else {
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp3_term_pr1_lhs::cells_3", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            // Thermodynamic main diagonal: [ x wp3(k,<t+1>) ]
            lhs_pr1_wp3(i - 1, k - 1) =
                (+((C8(i - 1) * invrs_tau_wp3_zt(i - 1, k - 1))) *
                 (((five * C8b(i - 1)) * ((Skw_zt(i - 1, k - 1) * Skw_zt(i - 1, k - 1)) *
                                          (Skw_zt(i - 1, k - 1) * Skw_zt(i - 1, k - 1)))) +
                  one));
          });
    }
  }
  // l_damp_wp3_Skw_squared
  return;
}

// =============================================================================
void wp3_term_ta_explicit_rhs(int nzm, int nzt, int ngrdcol, const grid& gr,        // In
                              Real2D wp4, Real2D rho_ds_zm, Real2D invrs_rho_ds_zt, // In
                              Real2D rhs_ta_wp3) {                                  // Out
  //
  // Description:
  // Turbulent advection of <w'^3>:  explicit portion of the code.
  //
  // This explicit discretization works generally for any PDF.
  //
  // The d<w'^3>/dt equation contains a turbulent advection term:
  //
  // - (1/rho_ds) * d( rho_ds * <w'^4> )/dz.
  //
  // The value of <w'^4> is found by integrating over the PDF of w, as detailed
  // in function calc_wp4_pdf, which is found in module pdf_closure_module in
  // pdf_closure_module.F90.
  //
  // The explicit discretization of this term is as follows:
  //
  // The values of <w'^3> are found on the thermodynamic levels, while the
  // values of <w'^4> are found on the momentum levels.  The values of
  // <w'^4>|_zt are originally calculated by the PDF on the thermodynamic
  // levels.  They are interpolated to the intermediate momentum levels as
  // <w'^4>.  Additionally, the values of rho_ds_zm are found on the momentum
  // levels, and the values of invrs_rho_ds_zt are found on the thermodynamic
  // levels.  At the intermediate momentum levels, the values of <w'^4> are
  // multiplied by rho_ds_zm.  Then, the derivative (d/dz) of that expression
  // is taken over the central thermodynamic level, where it is multiplied by
  // -invrs_rho_ds_zt.  This yields the desired result.
  //
  // ---------wp4_zt---------------------------------------------------- t(k+1)
  //
  // =========wp4(interp)===========rho_ds_zm=========================== m(k+1)
  //
  // ---------wp4_zt-----d( rho_ds_zm * wp4 )/dz-----invrs_rho_ds_zt---- t(k)
  //
  // =========wp4(interp)===========rho_ds_zm=========================== m(k)
  //
  // ---------wp4_zt---------------------------------------------------- t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is
  // used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
  //
  // References:
  //-----------------------------------------------------------------------

  // ---------------------- Input Variables ----------------------
  // wp4: <w'^4>                                   [m^4/s^4]
  // rho_ds_zm: Dry, static density at momentum level    [kg/m^3]
  // invrs_rho_ds_zt: Inv dry, static density at thermo level  [m^3/kg]
  // ---------------------- Output Variable ----------------------
  // rhs_ta_wp3: Rate of change of wp3 from turbulent advection  [m^3/s^4]
  // ---------------------- Output variables ----------------------

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_invrs_dzt = gr.invrs_dzt;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_ta_explicit_rhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set lower boundary to 0
          rhs_ta_wp3(i - 1, gr_k_lb_zt - 1) = zero;
          // Set upper boundary to 0
          rhs_ta_wp3(i - 1, gr_k_ub_zt - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_ta_explicit_rhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          rhs_ta_wp3(i - 1, k - 1) =
              ((-(invrs_rho_ds_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
               ((rho_ds_zm(i - 1, (k + 1) - 1) * wp4(i - 1, (k + 1) - 1)) -
                (rho_ds_zm(i - 1, k - 1) * wp4(i - 1, k - 1))));
        });
  }
  return;
}

// =============================================================================
void wp3_terms_bp1_pr2_rhs(int nzt, int ngrdcol, const grid& gr, Real2D C11_Skw_fnc, // In
                           Real2D thv_ds_zt, Real2D wp2thvp,                         // In
                           Real2D rhs_bp1_pr2_wp3) {                                 // Out
  //
  // Description:
  // Buoyancy production of w'^3 and w'^3 pressure term 2:  explicit portion of
  // the code.
  //
  // The d(w'^3)/dt equation contains a buoyancy production term:
  //
  // + 3 (g/thv_ds) w'^2th_v';
  //
  // and pressure term 2:
  //
  // - C_11 ( -3 w'^3 dw/dz + 3 (g/thv_ds) w'^2th_v' ).
  //
  // The w'^3 buoyancy production term is completely explicit, while w'^3
  // pressure term 2 has both implicit and explicit components.  The buoyancy
  // production term and the explicit portion of pressure term 2 are combined
  // and solved together as:
  //
  // + ( 1 - C_ll ) ( 3 (g/thv_ds) w'^2th_v' ).
  //
  // References:
  //-----------------------------------------------------------------------

  // -------------------- Input Variables --------------------
  // C11_Skw_fnc: C_11 parameter with Sk_w applied         [-]
  // thv_ds_zt: Dry, base-state theta_v at thermo. levs  [K]
  // wp2thvp: w'^2 th_v'                               [K m^2/s^2]
  // -------------------- Output Variable --------------------
  // rhs_bp1_pr2_wp3: RHS portion of wp3 from terms bp1 and pr2 [m^3/s^4]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_terms_bp1_pr2_rhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set lower boundary to 0
          rhs_bp1_pr2_wp3(i - 1, gr_k_lb_zt - 1) = zero;
          // Set upper boundary to 0
          rhs_bp1_pr2_wp3(i - 1, gr_k_ub_zt - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_terms_bp1_pr2_rhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          rhs_bp1_pr2_wp3(i - 1, k - 1) =
              (((+((one - C11_Skw_fnc(i - 1, k - 1))) * three) * (grav / thv_ds_zt(i - 1, k - 1))) *
               wp2thvp(i - 1, k - 1));
        });
  }
  return;
}

// =============================================================================
void wp3_term_pr_turb_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real1D C_wp3_pr_turb, // In
                          Real2D Kh_zt, Real2D wpthvp, Real2D dum_dz, Real2D dvm_dz,           // In
                          Real2D upwp, Real2D vpwp,                                            // In
                          Real2D thv_ds_zt,                                                    // In
                          Real2D rho_ds_zm, Real2D invrs_rho_ds_zt,                            // In
                          Real2D em, Real2D wp2,                                               // In
                          Real2D rhs_pr_turb_wp3,               // Out
                          bool l_use_tke_in_wp3_pr_turb_term) { // In
  // Description:
  //   Experimental term from CLUBB TRAC ticket #411. The derivative here is of
  //   the form:
  //   - C_15 * Kh * ∂{ grav / thv_ds * [w'th_v'(k) - w'th_v'(k-1)]
  //                    -[ u'w'(k) * ∂u(k)/∂z - u'w'(k-1) * ∂u(k-1)/∂z ]
  //                    -[ v'w'(k) * ∂v(k)/∂z - v'w'(k-1) * ∂v(k-1)/∂z ]  }/∂z.
  //
  //   This does not appear in Andre et al. 1976 or Bougeault et al. 1981, but
  //   is based on experiments in matching LES data.
  //
  // References:
  //   None
  //-----------------------------------------------------------------------

  // --------------------- Input Variables ---------------------
  // C_wp3_pr_turb: Model parameter C_wp3_pr_turb                [-]
  // Kh_zt: Eddy-diffusivity on thermo. levels      [m^2/s]
  // thv_ds_zt: Dry, base-state theta_v at thermo. levs [K]
  // invrs_rho_ds_zt: Inverse dry static density on thermo. levs [m^3/kg]
  // wpthvp: w'th_v'                                 [K m/s]
  // dum_dz: derivative of u wind with respect to z  [s^-1]
  // dvm_dz: derivative of v wind with respect to z  [s^-1]
  // upwp: u'w'                                    [m^2/s^2]
  // vpwp: v'w'                                    [m^2/s^2]
  // rho_ds_zm: Dry static density on mom. levels       [kg/m^3]
  // wp2: w'^2                                    [m^2/s^2]
  // em: Turbulence kinetic energy               [m^2/s^2]
  // l_use_tke_in_wp3_pr_turb_term: Use TKE formulation for wp3 pr_turb term
  // Return Variable
  // rhs_pr_turb_wp3: RHS portion of wp3 from pressure-turbulence correlation [m^3/s^4]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_invrs_dzt = gr.invrs_dzt;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_pr_turb_rhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set lower boundary to 0
          rhs_pr_turb_wp3(i - 1, gr_k_lb_zt - 1) = zero;
          // Set upper boundary to 0
          rhs_pr_turb_wp3(i - 1, gr_k_ub_zt - 1) = zero;
        });
  }
  if (!(l_use_tke_in_wp3_pr_turb_term)) {
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp3_term_pr_turb_rhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            rhs_pr_turb_wp3(i - 1, k - 1) =
                (((-(C_wp3_pr_turb(i - 1)) * Kh_zt(i - 1, k - 1)) * gr_invrs_dzt(i - 1, k - 1)) *
                 ((((grav / thv_ds_zt(i - 1, k - 1)) *
                    (wpthvp(i - 1, (k + 1) - 1) - wpthvp(i - 1, k - 1))) -
                   ((upwp(i - 1, (k + 1) - 1) * dum_dz(i - 1, (k + 1) - 1)) -
                    (upwp(i - 1, k - 1) * dum_dz(i - 1, k - 1)))) -
                  ((vpwp(i - 1, (k + 1) - 1) * dvm_dz(i - 1, (k + 1) - 1)) -
                   (vpwp(i - 1, k - 1) * dvm_dz(i - 1, k - 1)))));
          });
    }
  } else {
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp3_term_pr_turb_rhs::cells_3", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            rhs_pr_turb_wp3(i - 1, k - 1) =
                (((-(C_wp3_pr_turb(i - 1)) * invrs_rho_ds_zt(i - 1, k - 1)) *
                  gr_invrs_dzt(i - 1, k - 1)) *
                 (((rho_ds_zm(i - 1, (k + 1) - 1) * wp2(i - 1, (k + 1) - 1)) *
                   em(i - 1, (k + 1) - 1)) -
                  ((rho_ds_zm(i - 1, k - 1) * wp2(i - 1, k - 1)) * em(i - 1, k - 1))));
          });
    }
  }
  return;
}

// =============================================================================
void wp3_term_pr_dfsn_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, Real1D C_wp3_pr_dfsn, // In
                          Real2D rho_ds_zm, Real2D invrs_rho_ds_zt,                            // In
                          Real2D wp2up2, Real2D wp2vp2, Real2D wp4,                            // In
                          Real2D up2, Real2D vp2, Real2D wp2,                                  // In
                          Real2D rhs_pr_dfsn_wp3) { // Out
  //
  // Description:
  //
  // This term is intended to represent the "diffusion" part of the total wp3
  // pressure correlation.  The total wp3 pressure term, -3w'^2/rho*dp'/dz, can be
  // split into
  //
  //   -3w'^2/rho*dp'/dz = + 3p'/rho*d(w'^2)/dz - 3/rho*d(w'^2p')/dz
  //
  // using the product rule.  The second term on the RHS we consider to be the
  // diffusion part, calculated by this subroutine.  We replace the factor of 3
  // with a tunable parameter, C_wp3_pr_dfsn, and we replace p' with
  //
  //   p' ~ - rho * ( u_i*u_i - <u_i*u_i> ),
  //
  // following Lumley 1978.  The wp3 pressure diffusion term then becomes
  //
  //   + C_wp3_pr_dfsn / rho * ( d( rho*( <w'^2u_i'u_i'> - <w'^2>*<u_i'u_i'> ) )/dz )
  //
  // References:
  //   Lumley 1978, p. 170.  See eq. 6.47 and accompanying discussion.
  //-----------------------------------------------------------------------

  // ---------------------- Input Variables ----------------------
  // C_wp3_pr_dfsn: Model parameter C_wp3_pr_dfsn              [-]
  // invrs_rho_ds_zt: Inverse dry static density (thermo levels) [m^3/kg]
  // rho_ds_zm: Dry static density on mom. levels          [kg/m^3]
  // wp2up2: w'^2u'^2 on momentum levels                [m^4/s^4]
  // wp2vp2: w'^2v'^2 on momentum levels                [m^4/s^4]
  // wp4: w'^4 on momentum levels                    [m^4/s^4]
  // up2: u'^2 on momentum levels                    [m^2/s^2]
  // vp2: v'^2 on momentum levels                    [m^2/s^2]
  // wp2: w'^2 on momentum levels                    [m^2/s^2]
  // ---------------------- Output Variable ----------------------
  // rhs_pr_dfsn_wp3: RHS portion of wp3 from pressure-diffusion correlation [m^3/s^4]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  Real2D wp2uip2("wp3_term_pr_dfsn_rhs::wp2uip2", ngrdcol,
                 nzm); // 4th-order moment sum <w'^2u_i'u_i'>     [m^4/s^4]
  Real2D wp2_uip2("wp3_term_pr_dfsn_rhs::wp2_uip2", ngrdcol,
                  nzm); // 2nd-order moment sum <w'^2>*<u_i'u_i'>  [m^4/s^4]
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;
  const auto gr_invrs_dzt = gr.invrs_dzt;

  {
    const int extent_0 = ((1) > 0 ? (1) <= (nzm) : (1) >= (nzm)) ? ((nzm) - (1)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_pr_dfsn_rhs::cells_1", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (1) + ((cell / (1 * extent_1)) % extent_0) * (1);
          wp2uip2(i - 1, k - 1) =
              ((wp2up2(i - 1, k - 1) + wp2vp2(i - 1, k - 1)) + wp4(i - 1, k - 1));
          wp2_uip2(i - 1, k - 1) =
              (((wp2(i - 1, k - 1) * up2(i - 1, k - 1)) + (wp2(i - 1, k - 1) * vp2(i - 1, k - 1))) +
               (wp2(i - 1, k - 1) * wp2(i - 1, k - 1)));
        });
  }
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_pr_dfsn_rhs::cells_2", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set lower boundary to 0
          rhs_pr_dfsn_wp3(i - 1, gr_k_lb_zt - 1) = zero;
          // Set upper boundary to 0
          rhs_pr_dfsn_wp3(i - 1, gr_k_ub_zt - 1) = zero;
        });
  }
  {
    const int extent_0 =
        ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
    const int extent_1 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_pr_dfsn_rhs::cells_3", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_1) * (1);
          const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
          rhs_pr_dfsn_wp3(i - 1, k - 1) =
              (((+(C_wp3_pr_dfsn(i - 1)) * invrs_rho_ds_zt(i - 1, k - 1)) *
                gr_invrs_dzt(i - 1, k - 1)) *
               ((rho_ds_zm(i - 1, (k + 1) - 1) *
                 (wp2uip2(i - 1, (k + 1) - 1) - wp2_uip2(i - 1, (k + 1) - 1))) -
                (rho_ds_zm(i - 1, k - 1) * (wp2uip2(i - 1, k - 1) - wp2_uip2(i - 1, k - 1)))));
        });
  }
  return;
}

// =============================================================================
void wp3_term_pr1_rhs(int nzt, int ngrdcol, const grid& gr,               // In
                      Real1D C8, Real1D C8b,                              // In
                      Real2D invrs_tau_wp3_zt, Real2D Skw_zt, Real2D wp3, // In
                      bool l_damp_wp3_Skw_squared,                        // In
                      Real2D rhs_pr1_wp3) {                               // Out
  //
  // Description:
  // Pressure term 1 for w'^3:  explicit portion of the code.
  //
  // Pressure term 1 is the term:
  //
  // - (C_8/tau_w3t) * ( C_8b * Sk_wt^2 + 1 ) * w'^3;
  //
  // where Sk_wt = w'^3 / (w'^2)^(3/2).
  //
  // This term needs to be linearized, so function L(w'^3) is defined to be
  // equal to this term (pressure term 1), such that:
  //
  // L(w'^3) = - (C_8/tau_w3t) * ( C_8b * (w'^3)^3 / (w'^2)^3 + w'^3 ).
  //
  // A Taylor Series expansion (truncated after the first derivative term) of
  // L(w'^3) around w'^3 = w'^3(t) is used to linearize pressure term 1.
  // Evaluating L(w'^3) at w'^3(t+1):
  //
  // L( w'^3(t+1) ) = L( w'^3(t) )
  //                  + ( d L(w'^3) / d w'^3 )|_(w'^3=w'^3(t))
  //                    * ( w'^3(t+1) - w'^3(t) ).
  //
  // After evaluating the expression above, the term has become linearized.  It
  // is broken down into implicit (LHS) and explicit (RHS) components.
  // The explicit portion is:
  //
  // + (C_8/tau_w3t) * ( 2 * C_8b * Sk_wt^2 + 1 ) * w'^3(t).
  //
  // Timestep index (t) stands for the index of the current timestep, while
  // timestep index (t+1) stands for the index of the next timestep, which is
  // being advanced to in solving the d(w'^3)/dt equation.
  //
  // The values of w'^3 are found on the thermodynamic levels, as are the
  // values of tau_w3t and Sk_wt (in Sk_wt, w'^3 is found on thermodynamic
  // levels and w'^2 is interpolated to thermodynamic levels).
  //
  // References:
  //-----------------------------------------------------------------------

  // --------------------- Input Variables ---------------------
  // C8: Model parameter C_8                        [-]
  // C8b: Model parameter C_8b                       [-]
  // invrs_tau_wp3_zt: Inverse time-scale tau at thermodynamic levels  [1/s]
  // Skw_zt: Skewness of w at thermodynamic levels      [-]
  // wp3: w'^3                                       [m^3/s^3]
  // l_damp_wp3_Skw_squared: Set damping on wp3 to use Skw^2 rather than Skw^4
  // --------------------- Output Variable ---------------------
  // rhs_pr1_wp3: RHS portion of wp3 from pressure term 1  [m^3/s^4]

  // Local arrays belong to this call. Internal indices retain Native one-based
  // values; View accesses subtract one. Parameter indices are already zero-based.
  // Capture component Views and native one-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt = gr.k_lb_zt;
  const auto gr_k_ub_zt = gr.k_ub_zt;

  // Set lower boundary to 0
  {
    const int extent_0 =
        ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
    Kokkos::parallel_for(
        "wp3_term_pr1_rhs::cells_1", extent_0, KOKKOS_LAMBDA(const int cell) {
          const int i = (1) + ((cell / (1)) % extent_0) * (1);
          // Set lower boundary to 0
          rhs_pr1_wp3(i - 1, gr_k_lb_zt - 1) = zero;
          // Set upper boundary to 0
          rhs_pr1_wp3(i - 1, gr_k_ub_zt - 1) = zero;
        });
  }
  // Calculate term at all interior grid levels.
  if (l_damp_wp3_Skw_squared) {
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp3_term_pr1_rhs::cells_2", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            rhs_pr1_wp3(i - 1, k - 1) =
                ((+((C8(i - 1) * invrs_tau_wp3_zt(i - 1, k - 1))) *
                  ((two * C8b(i - 1)) * (Skw_zt(i - 1, k - 1) * Skw_zt(i - 1, k - 1)))) *
                 wp3(i - 1, k - 1));
          });
    }
  } else {
    {
      const int extent_0 =
          ((1) > 0 ? (2) <= ((nzt - 1)) : (2) >= ((nzt - 1))) ? (((nzt - 1)) - (2)) / (1) + 1 : 0;
      const int extent_1 =
          ((1) > 0 ? (1) <= (ngrdcol) : (1) >= (ngrdcol)) ? ((ngrdcol) - (1)) / (1) + 1 : 0;
      Kokkos::parallel_for(
          "wp3_term_pr1_rhs::cells_3", extent_0 * extent_1, KOKKOS_LAMBDA(const int cell) {
            const int i = (1) + ((cell / (1)) % extent_1) * (1);
            const int k = (2) + ((cell / (1 * extent_1)) % extent_0) * (1);
            rhs_pr1_wp3(i - 1, k - 1) =
                ((+((C8(i - 1) * invrs_tau_wp3_zt(i - 1, k - 1))) *
                  ((four * C8b(i - 1)) * ((Skw_zt(i - 1, k - 1) * Skw_zt(i - 1, k - 1)) *
                                          (Skw_zt(i - 1, k - 1) * Skw_zt(i - 1, k - 1))))) *
                 wp3(i - 1, k - 1));
          });
    }
  }
  // l_damp_wp3_Skw_squared
  return;
}

} // namespace

} // namespace clubb_kokkos
