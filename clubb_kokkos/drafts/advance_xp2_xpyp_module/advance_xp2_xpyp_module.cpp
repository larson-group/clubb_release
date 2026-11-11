// Advance scalar variances/covariances and horizontal turbulent variances.
//
// advance_xp2_xpyp preserves Native advection/dissipation/pressure solves,
// clipping, scalar emulation and sampled Native budgets. update_xp2_mc computes
// evaporation contributions using a column-local ordered precipitation profile.
// Physics Views are column first; Native matrix arrays remain bands first until
// the existing matrix_solver_wrapper conversion boundary. All scratch is local.
// Sponge-enabled initialization remains excluded with JAX scope; see INTEGRATION.md.
// References: Golaz et al. (2002), JAS59,3540–3551, Eqs.13–15;
// https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:xpyp_eqns
#include "advance_xp2_xpyp_module.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "model_flags.hpp"
#include "clip_explicit.hpp" // Native selected constants pending at their owner.
#include "matrix_solver_wrapper.hpp"
#include "error_code.hpp"
#include "f2k/procedure_bridge.hpp"
#include "f2k_dependencies.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <algorithm>
#include <utility>
namespace clubb_kokkos {
namespace {
// Private named constants to avoid string comparisons.
constexpr int xp2_xpyp_rtp2       = 1;
constexpr int xp2_xpyp_thlp2      = 2;
constexpr int xp2_xpyp_rtpthlp    = 3;
constexpr int xp2_xpyp_up2_vp2    = 4;
constexpr int xp2_xpyp_vp2        = 6;
constexpr int xp2_xpyp_up2        = 5;
constexpr int xp2_xpyp_scalars    = 7;
constexpr int xp2_xpyp_sclrp2     = 8;
constexpr int xp2_xpyp_sclrprtp   = 9;
constexpr int xp2_xpyp_sclrpthlp  = 10;
constexpr int xp2_xpyp_single_lhs = 11;
constexpr int ndiags2             = 2;
constexpr int ndiags3             = 3;
constexpr int ndiags5             = 5;
template <class View> std::string dump_values(View values) {
  // Format Native arrays in first-index-fastest order at a host error boundary.
  auto host = to_host_space(values);
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  for (std::size_t n = 0; n < host.size(); ++n)
    out << host.data()[n] << ' ';
  return out.str();
}
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26>
void solve_xp2_xpyp_with_single_lhs(
    int nzm, int nzt, int ngrdcol, int sclr_dim, Array0 sclr_tol, const grid& gr,
    const sclr_idx_type& sclr_idx, Array1 C2x, Array2 invrs_tau_xp2_zm, Array3 rtm, Array4 thlm,
    Array5 wprtp, Array6 wpthlp, Array7 rtp2_forcing, Array8 thlp2_forcing, Array9 rtpthlp_forcing,
    Array10 sclrm, Array11 wpsclrp, Array12 lhs_ta, Array13 lhs_ma, Array14 lhs_diff,
    Array15 rhs_ta_wprtp2, Array16 rhs_ta_wpthlp2, Array17 rhs_ta_wprtpthlp,
    Array18 rhs_ta_wpsclrp2, Array19 rhs_ta_wprtpsclrp, Array20 rhs_ta_wpthlpsclrp, real_t dt,
    bool l_scalar_calc, bool l_lmm_stepping, bool l_implemented, int tridiag_solve_method,
    f2k::stats_type& stats, Array21 rtp2, Array22 thlp2, Array23 rtpthlp, Array24 sclrp2,
    Array25 sclrprtp, Array26 sclrpthlp, err_info_type& err_info);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34>
void solve_xp2_xpyp_with_multiple_lhs(
    int nzm, int nzt, int ngrdcol, int sclr_dim, Array0 sclr_tol, const grid& gr,
    const sclr_idx_type& sclr_idx, Array1 C2rt_1d, Array2 C2thl_1d, Array3 C2rtthl_1d,
    Array4 C2sclr_1d, Array5 invrs_tau_xp2_zm, Array6 rtm, Array7 thlm, Array8 wprtp, Array9 wpthlp,
    Array10 rtp2_forcing, Array11 thlp2_forcing, Array12 rtpthlp_forcing, Array13 sclrm,
    Array14 wpsclrp, Array15 lhs_ta_wprtp2, Array16 lhs_ta_wpthlp2, Array17 lhs_ta_wprtpthlp,
    Array18 lhs_ta_wpsclrp2, Array19 lhs_ta_wprtpsclrp, Array20 lhs_ta_wpthlpsclrp, Array21 lhs_ma,
    Array22 lhs_diff, Array23 rhs_ta_wprtp2, Array24 rhs_ta_wpthlp2, Array25 rhs_ta_wprtpthlp,
    Array26 rhs_ta_wpsclrp2, Array27 rhs_ta_wprtpsclrp, Array28 rhs_ta_wpthlpsclrp, real_t dt,
    int iiPDF_type, bool l_scalar_calc, bool l_lmm_stepping, bool l_implemented,
    int tridiag_solve_method, f2k::stats_type& stats, Array29 rtp2, Array30 thlp2, Array31 rtpthlp,
    Array32 sclrp2, Array33 sclrprtp, Array34 sclrpthlp, err_info_type& err_info);
template <class Array0, class Array1, class Array2, class Array3, class Array4>
void xp2_xpyp_lhs(int nzm, int ngrdcol, real_t dt, const grid& gr, Array0 lhs_ta, Array1 lhs_ma,
                  Array2 lhs_diff, Array3 lhs_dp1, Array4 lhs);
template <class Array0, class Array1, class Array2>
void xp2_xpyp_solve(int nzm, int ngrdcol, int solve_type, int nrhs, bool l_implemented,
                    const grid& gr, int tridiag_solve_method, f2k::stats_type& stats, Array0 rhs,
                    Array1 lhs, err_info_type& err_info, Array2 xapxbp);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
void stats_finalize_xp2_xpyp_terms(int nzm, int ngrdcol, int solve_type, const grid& gr,
                                   Array0 xapxbp, Array1 lhs_dp1, Array2 lhs_dp1_pr1,
                                   Array3 lhs_diff, Array4 lhs_ta, Array5 lhs_ma,
                                   f2k::stats_type& stats);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20>
void xp2_xpyp_uv_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, int solve_type, real_t dt,
                     Array0 wp2, Array1 wpthvp, Array2 invrs_tau_C4_zm, Array3 invrs_tau_C14_zm,
                     Array4 xam, Array5 xbm, Array6 wpxap, Array7 wpxbp, Array8 xap2, Array9 xbp2,
                     Array10 thv_ds_zm, Array11 C4, Array12 C_uu_shr, Array13 C_uu_buoy,
                     Array14 C14, Array15 lhs_splat_wp2, Array16 lhs_ta, Array17 rhs_ta,
                     Array18 lhs_dp1_C4, Array19 lhs_dp1_C14, f2k::stats_type& stats, Array20 rhs);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10>
void xp2_xpyp_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, int solve_type, real_t dt,
                  Array0 wpxap, Array1 wpxbp, Array2 xam, Array3 xbm, Array4 xapxbp,
                  Array5 xpyp_forcing, Array6 Cn, Array7 invrs_tau_zm, real_t threshold,
                  Array8 lhs_ta, Array9 rhs_ta, f2k::stats_type& stats, Array10 rhs);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43>
void calc_xp2_xpyp_ta_terms(
    int nzm, int nzt, int ngrdcol, int sclr_dim, const grid& gr, Array0 wprtp, Array1 wprtp2,
    Array2 wpthlp, Array3 wpthlp2, Array4 wprtpthlp, Array5 rtp2, Array6 thlp2, Array7 rtpthlp,
    Array8 upwp, Array9 vpwp, Array10 up2, Array11 vp2, Array12 wp2, Array13 wp2_zt,
    Array14 wpsclrp, Array15 wpsclrp2, Array16 wpsclrprtp, Array17 wpsclrpthlp, Array18 sclrp2,
    Array19 sclrprtp, Array20 sclrpthlp, Array21 rho_ds_zt, Array22 invrs_rho_ds_zm,
    Array23 rho_ds_zm, Array24 wp3_on_wp2, Array25 wp3_on_wp2_zt, Array26 sigma_sqd_w,
    const implicit_coefs_terms& pdf_implicit_coefs_terms, bool l_scalar_calc, Array27 beta,
    int iiPDF_type, bool l_upwind_xpyp_ta, bool l_godunov_upwind_xpyp_ta, f2k::stats_type& stats,
    Array28 lhs_ta_wprtp2, Array29 lhs_ta_wpthlp2, Array30 lhs_ta_wprtpthlp, Array31 lhs_ta_wpup2,
    Array32 lhs_ta_wpvp2, Array33 lhs_ta_wpsclrp2, Array34 lhs_ta_wprtpsclrp,
    Array35 lhs_ta_wpthlpsclrp, Array36 rhs_ta_wprtp2, Array37 rhs_ta_wpthlp2,
    Array38 rhs_ta_wprtpthlp, Array39 rhs_ta_wpup2, Array40 rhs_ta_wpvp2, Array41 rhs_ta_wpsclrp2,
    Array42 rhs_ta_wprtpsclrp, Array43 rhs_ta_wpthlpsclrp);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
void term_tp_rhs(int nzm, int nzt, int ngrdcol, Array0 xam, Array1 xbm, Array2 wpxbp, Array3 wpxap,
                 Array4 invrs_dzm, Array5 rhs);
template <class Array0, class Array1, class Array2>
void term_dp1_lhs(int nzm, int ngrdcol, const grid& gr, Array0 Cn, Array1 invrs_tau_zm, Array2 lhs);
template <class Array0, class Array1, class Array2>
void term_dp1_rhs(int nzm, int ngrdcol, Array0 Cn, Array1 invrs_tau_zm, real_t threshold,
                  Array2 rhs);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void term_pr1(int nzm, int ngrdcol, Array0 C4, Array1 C14, Array2 xbp2, Array3 wp2,
              Array4 invrs_tau_C4_zm, Array5 invrs_tau_C14_zm, Array6 rhs);
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8>
void term_pr2(int nzm, int nzt, int ngrdcol, const grid& gr, Array0 C_uu_shr, Array1 C_uu_buoy,
              Array2 thv_ds_zm, Array3 wpthvp, Array4 upwp, Array5 vpwp, Array6 um, Array7 vm,
              Array8 rhs_pr2);
template <class Array0, class Array1>
void pos_definite_variances(int nzm, int ngrdcol, const grid& gr, int solve_type,
                            int fill_holes_type, real_t dt, real_t tolerance, Array0 rho_ds_zm,
                            f2k::stats_type& stats, Array1 xp2_np1);
} // namespace

// nzm: Number of momentum vertical levels
// nzt: Number of thermodynamic vertical levels
// ngrdcol: Number of grid columns
// sclr_dim: Number of passive scalars
// sclr_tol: Threshold(s) on the passive scalars  [units vary]
// invrs_tau_xp2_zm: Inverse time-scale for xp2 on momentum levels [1/s]
// invrs_tau_C4_zm: Inverse time-scale for C4 terms on mom. levels [1/s]
// invrs_tau_C14_zm: Inverse time-scale for C14 terms on mom. levels [1/s]
// wm_zm: w-wind component on momentum levels   [m/s]
// rtm: Total water mixing ratio (t-levs)     [kg/kg]
// wprtp: <w'r_t'> (momentum levels)            [(m/s)(kg/kg)]
// thlm: Liquid potential temp. (t-levs)       [K]
// wpthlp: <w'th_l'> (momentum levels)           [(m K)/s]
// wpthvp: <w'th_v'> (momentum levels)           [(m K)/s]
// um: u wind (thermodynamic levels)         [m/s]
// vm: v wind (thermodynamic levels)         [m/s]
// wp2: <w'^2> (momentum levels)              [m^2/s^2]
// wp3: <w'^3> (thermodynamic levels)         [m^3/s^3]
// upwp: <u'w'> (momentum levels)              [m^2/s^2]
// vpwp: <v'w'> (momentum levels)              [m^2/s^2]
// sigma_sqd_w: sigma_sqd_w (momentum levels)         [-]
// wprtp2: <w'r_t'^2> (thermodynamic levels)     [m/s (kg/kg)^2]
// wpthlp2: <w'th_l'^2> (thermodynamic levels)    [m/s K^2]
// wprtpthlp: <w'r_t'th_l'> (thermodynamic levels)  [m/s (kg/kg) K]
// Kh_zt: Eddy diffusivity on thermo. levels    [m^2/s]
// rtp2_forcing: <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
// thlp2_forcing: <th_l'^2> forcing (momentum levels)   [K^2/s]
// rtpthlp_forcing: <r_t'th_l'> forcing (momentum levels) [(kg/kg)K/s]
// rho_ds_zm: Dry, static density on momentum levs. [kg/m^3]
// rho_ds_zt: Dry, static density on thermo. levels [kg/m^3]
// invrs_rho_ds_zm: Inv. dry, static density @ mom. levs. [m^3/kg]
// thv_ds_zm: Dry, base-state theta_v on mom. levs. [K]
// cloud_frac: Cloud fraction (thermodynamic levels) [-]
// pdf_implicit_coefs_terms: Implicit coefs / explicit terms [units vary]
// dt: Model timestep                                [s]
// fcor_y: Nontraditional Coriolis parameter         [s^-1] Meridional planetary vorticity.
// Proportional to cos(latitude)
// sclrm: Mean value; pass. scalar (t-levs.) [{sclr units}]
// wpsclrp: <w'sclr'> (momentum levels)        [m/s{sclr units}]
// wpsclrp2: <w'sclr'^2> (thermodynamic levels) [m/s{sclr units}^2]
// wpsclrprtp: <w'sclr'r_t'> (thermo. levels)     [m/s{sclr units)kg/kg]
// wpsclrpthlp: <w'sclr'th_l'> (thermo. levels)    [m/s{sclr units}K]
// lhs_splat_wp2: LHS coefficient of wp2 splatting term  [1/s]
// clubb_params: Array of CLUBB's tunable parameters    [units vary]
// nu_vert_res_dep: Vertical resolution dependent nu values
// iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use for
// the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component
// PDF.
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// fill_holes_type: Specifier for which hole filling method to use
// l_ho_nontrad_coriolis: Flag to implement the nontraditional Coriolis terms in the prognostic
// equations of <w'w'>, <u'w'>, and <u'u'>.
// l_min_xp2_from_corr_wx: Flag to base the threshold minimum value of xp2 (rtp2 and thlp2) on
// keeping the overall correlation of w and x within the limits of -max_mag_correlation_flux to
// max_mag_correlation_flux.
// l_C2_cloud_frac: Flag to use cloud fraction to adjust the value of the turbulent dissipation
// coefficient, C2.
// l_upwind_xpyp_ta: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, & sclrpthlp.
// l_godunov_upwind_xpyp_ta: This flag determines whether we want to use a Godunov-like upwind
// differencing approximation rather than a centered differencing for turbulent or mean advection
// terms. It affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, & sclrpthlp.
// l_lmm_stepping: Apply Linear Multistep Method (LMM) Stepping
// l_implemented: True if CLUBB is being implemented and run in a host model
// rtp2: <r_t'^2>                      [(kg/kg)^2]
// thlp2: <th_l'^2>                     [K^2]
// rtpthlp: <r_t'th_l'>                   [(kg K)/kg]
// up2: <u'^2>                        [m^2/s^2]
// vp2: <v'^2>                        [m^2/s^2]
// err_info: err_info struct containing err_code and err_header
void advance_xp2_xpyp(int nzm, int nzt, int ngrdcol, int sclr_dim, Real1D sclr_tol,
                      const grid& gr,                                  // In
                      const sclr_idx_type& sclr_idx,                   // In
                      Real2D invrs_tau_xp2_zm, Real2D invrs_tau_C4_zm, // In
                      Real2D invrs_tau_C14_zm, Real2D wm_zm,           // In
                      Real2D rtm, Real2D wprtp, Real2D thlm, Real2D wpthlp, Real2D wpthvp,
                      Real2D um,                                                         // In
                      Real2D vm,                                                         // In
                      Real2D wp2, Real2D wp3, Real2D upwp, Real2D vpwp,                  // In
                      Real2D sigma_sqd_w, Real2D wprtp2, Real2D wpthlp2,                 // In
                      Real2D wprtpthlp, Real2D Kh_zt, Real2D rtp2_forcing,               // In
                      Real2D thlp2_forcing, Real2D rtpthlp_forcing,                      // In
                      Real2D rho_ds_zm, Real2D rho_ds_zt, Real2D invrs_rho_ds_zm,        // In
                      Real2D thv_ds_zm, Real2D cloud_frac,                               // In
                      const implicit_coefs_terms& pdf_implicit_coefs_terms,              // In
                      real_t dt, Real1D fcor_y,                                          // In
                      Real3D sclrm, Real3D wpsclrp,                                      // In
                      Real3D wpsclrp2, Real3D wpsclrprtp, Real3D wpsclrpthlp,            // In
                      Real2D lhs_splat_wp2,                                              // In
                      Real2D clubb_params, const nu_vertical_res_dep& nu_vert_res_dep,   // In
                      int iiPDF_type,                                                    // In
                      int tridiag_solve_method,                                          // In
                      int fill_holes_type,                                               // In
                      bool l_ho_nontrad_coriolis,                                        // In
                      bool l_min_xp2_from_corr_wx,                                       // In
                      bool l_C2_cloud_frac,                                              // In
                      bool l_upwind_xpyp_ta,                                             // In
                      bool l_godunov_upwind_xpyp_ta,                                     // In
                      bool l_lmm_stepping,                                               // In
                      bool l_implemented,                                                // In
                      f2k::stats_type& stats,                                            // InOut
                      Real2D rtp2, Real2D thlp2, Real2D rtpthlp, Real2D up2, Real2D vp2, // InOut
                      Real3D sclrp2, Real3D sclrprtp, Real3D sclrpthlp,
                      err_info_type& err_info) // InOut
{
  // Description:
  // Prognose scalar variances, scalar covariances, and horizontal turbulence components.
  // References:
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:xpyp_eqns
  //
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:up2_vp2_eqns
  //
  //   Eqn. 13, 14, 15  on p. 3545 of
  //   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
  //     Method and Model Description'' Golaz, et al. (2002)
  //   JAS, Vol. 59, pp. 3540--3551.
  // See also:
  //   ``Equations for CLUBB'', Section 4:
  //   /Steady-state solution for the variances/
  //-----------------------------------------------------------------------
  //------------------------------ Input Variables ------------------------------
  // Passive scalar input
  //------------------- InOut Variables -------------------
  //------------------------------ Input/Output Variables ------------------------------
  // An attribute of (inout) is also needed to import the value of the variances
  // at the surface.  Brian.  12/18/05.
  // Passive scalar output
  //------------------------------ Local Variables ------------------------------
  // Eddy Diffusion for Variances and Covariances.
  // Turbulent advection terms
  // Implicit (LHS) turbulent advection terms
  // Explicit (RHS) turbulent advection terms
  // Implicit (LHS) turbulent advection terms for scalars
  // Explicit (RHS) turbulent advection terms for scalars
  // Constant parameters
  // Minimum value of cloud fraction to multiply C2 by to calculate the
  // adjusted value of C2.
  // Loop indices
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rtp2_old(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rtp2_old"), ngrdcol,
      nzm); // Saved value of <r_t'^2>         [(kg/kg)^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> thlp2_old(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::thlp2_old"), ngrdcol,
      nzm); // Saved value of <th_l'^2>        [K^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rtpthlp_old(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rtpthlp_old"), ngrdcol,
      nzm); // Saved value of <r_t'th_l'>      [(kg K)/kg]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> up2_old(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::up2_old"), ngrdcol,
      nzm); // Saved value of <u'^2>           [m^2/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vp2_old(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::vp2_old"), ngrdcol,
      nzm); // Saved value of <v'^2>           [m^2/s^2]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sclrp2_old(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::sclrp2_old"), ngrdcol, nzm,
      sclr_dim); // Saved value of <sclr'^2>     [units vary]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sclrprtp_old(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::sclrprtp_old"), ngrdcol,
      nzm, sclr_dim); // Saved value of <sclr'rt'>    [units vary]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sclrpthlp_old(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::sclrpthlp_old"), ngrdcol,
      nzm, sclr_dim); // Saved value of <sclr'thl'>   [units vary]
  real_t C2rt;        // CLUBB tunable parameter C2rt
  real_t C2thl;       // CLUBB tunable parameter C2thl
  real_t C2rtthl;     // CLUBB tunable parameter C2rtthl
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C2sclr_1d(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::C2sclr_1d"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C2rt_1d(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::C2rt_1d"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C2thl_1d(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::C2thl_1d"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C2rtthl_1d(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::C2rtthl_1d"), ngrdcol,
      nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C4_1d(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::C4_1d"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> C14_1d(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::C14_1d"), ngrdcol, nzm);
  real_t threshold; // Minimum value for variances                   [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> threshold_array(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::threshold_array"), ngrdcol,
      nzm); // Minimum value for variances [units vary]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs"), ndiags3, ngrdcol,
      nzm); // Tridiagonal matrix
  Kokkos::View<real_t***, Kokkos::LayoutLeft> uv_rhs(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::uv_rhs"), ngrdcol, nzm,
      2); // RHS vectors of tridiagonal system for up2/vp2
  Kokkos::View<real_t***, Kokkos::LayoutLeft> uv_solution(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::uv_solution"), ngrdcol,
      nzm, 2); // Solution to the tridiagonal system for up2/vp2
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Kw2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::Kw2"), ngrdcol,
      nzt); // For rtp2, thlp2, rtpthlp, and passive scalars  [m^2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Kw9(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::Kw9"), ngrdcol,
      nzt); // For up2 and vp2                                [m^2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Kw2_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::Kw2_zm"), ngrdcol,
      nzm); // Eddy diffusivity coefficient, momentum levels [m2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> Kw9_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::Kw9_zm"), ngrdcol,
      nzm); // Eddy diffusivity coefficient, momentum levels [m2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rtpthlp_chnge(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rtpthlp_chnge"), ngrdcol,
      nzm); // Net change in r_t'th_l' due to clipping [(kg/kg) K]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sclrprtp_chnge(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::sclrprtp_chnge"), ngrdcol,
      nzm, sclr_dim); // Net change in sclr'r_t' due to clipping  [units vary]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sclrpthlp_chnge(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::sclrpthlp_chnge"), ngrdcol,
      nzm, sclr_dim); // Net change in sclr'th_l' due to clipping [units vary]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ta_wprtp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_ta_wprtp2"), ndiags3,
      ngrdcol, nzm); // For <w'rt'^2>
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ta_wpthlp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_ta_wpthlp2"), ndiags3,
      ngrdcol, nzm); // For <w'thl'^2>
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ta_wprtpthlp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_ta_wprtpthlp"),
      ndiags3, ngrdcol, nzm); // For <w'rt'thl'>
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ta_wpup2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_ta_wpup2"), ndiags3,
      ngrdcol, nzm); // For <w'u'^2>
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ta_wpvp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_ta_wpvp2"), ndiags3,
      ngrdcol, nzm); // For <w'v'^2>
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_ta_wprtp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rhs_ta_wprtp2"), ngrdcol,
      nzm); // For <w'rt'^2>
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_ta_wpthlp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rhs_ta_wpthlp2"), ngrdcol,
      nzm); // For <w'thl'^2>
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_ta_wprtpthlp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rhs_ta_wprtpthlp"),
      ngrdcol, nzm); // For <w'rt'thl'>
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_ta_wpup2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rhs_ta_wpup2"), ngrdcol,
      nzm); // For <w'u'^2>
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_ta_wpvp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rhs_ta_wpvp2"), ngrdcol,
      nzm); // For <w'v'^2>
  Kokkos::View<real_t****, Kokkos::LayoutLeft> lhs_ta_wpsclrp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_ta_wpsclrp2"), ndiags3,
      ngrdcol, nzm, sclr_dim); // For <w'sclr'^2>
  Kokkos::View<real_t****, Kokkos::LayoutLeft> lhs_ta_wprtpsclrp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_ta_wprtpsclrp"),
      ndiags3, ngrdcol, nzm, sclr_dim); // For <w'rt'sclr'>
  Kokkos::View<real_t****, Kokkos::LayoutLeft> lhs_ta_wpthlpsclrp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_ta_wpthlpsclrp"),
      ndiags3, ngrdcol, nzm, sclr_dim); // For <w'thl'sclr'>
  Kokkos::View<real_t***, Kokkos::LayoutLeft> rhs_ta_wpsclrp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rhs_ta_wpsclrp2"), ngrdcol,
      nzm, sclr_dim); // For <w'sclr'^2>
  Kokkos::View<real_t***, Kokkos::LayoutLeft> rhs_ta_wprtpsclrp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rhs_ta_wprtpsclrp"),
      ngrdcol, nzm, sclr_dim); // For <w'sclr'rt'>
  Kokkos::View<real_t***, Kokkos::LayoutLeft> rhs_ta_wpthlpsclrp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rhs_ta_wpthlpsclrp"),
      ngrdcol, nzm, sclr_dim); // For <w'sclr'thl'>
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_diff(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_diff"), ndiags3,
      ngrdcol, nzm); // Diffusion contributions to lhs, dissipation term 2
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_diff_uv(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_diff_uv"), ndiags3,
      ngrdcol, nzm); // Diffusion contributions to lhs for <w'u'^2> and <w'v'^2>
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs_ma(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_ma"), ndiags3, ngrdcol,
      nzm); // Mean advection contributions to lhs
  bool l_scalar_calc;
  constexpr bool l_clip_large_rtp2     = true;
  constexpr real_t min_cloud_frac_mult = real_t(0.10);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> lhs_dp1(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_dp1"), ngrdcol,
      nzm); // LHS dissipation term 1
  Kokkos::View<real_t**, Kokkos::LayoutLeft> lhs_dp1_C14(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_dp1_C14"), ngrdcol,
      nzm); // LHS dissipation term 1, for up2 vp2 using C14
  Kokkos::View<real_t**, Kokkos::LayoutLeft> lhs_dp1_C4(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::lhs_dp1_C4"), ngrdcol,
      nzm); // LHS dissipation term 1, for up2 vp2 using C4
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rtm_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::rtm_zm"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> cloud_frac_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::cloud_frac_zm"), ngrdcol,
      nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_tmp_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_tmp_zm"), ngrdcol,
      nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wp3_on_wp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::wp3_on_wp2"), ngrdcol,
      nzm); // Smoothed version of <w'^3>/<w'^2> zm [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wp2_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::wp2_zt"), ngrdcol,
      nzt); // <w'^2> interpolated to thermo. levels [m^2/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wp3_on_wp2_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::wp3_on_wp2_zt"), ngrdcol,
      nzt); // Smoothed version of <w'^3>/<w'^2> zt [m/s]
  int sclr;
  int k;
  int i;
  bool l_single_solve_possible;

  //------------------------------ Begin Code ------------------------------

  f2k::calc_wp3_on_wp2(nzm, nzt, ngrdcol, gr, wp2, wp3, // In
                       wp3_on_wp2, wp3_on_wp2_zt);      // Out

  Kokkos::deep_copy(
      Kokkos::subview(wp2_zt, Kokkos::ALL, Kokkos::ALL),
      zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wp2, Kokkos::ALL, Kokkos::ALL), w_tol_sqd));

  // Use 3 different values of C2 for rtp2, thlp2, rtpthlp.
  if (l_C2_cloud_frac) {

    // Interpolate cloud_frac to momentum levels.
    Kokkos::deep_copy(cloud_frac_zm, zt2zm_api(nzm, nzt, ngrdcol, gr, cloud_frac));

    Kokkos::parallel_for(
        "advance_xp2_xpyp::loop1", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);

          if (cloud_frac_zm(i - 1, k - 1) >= cloud_frac_min) {
            C2rt_1d(i - 1, k - 1) =
                clubb_params(i - 1, iC2rt) * ((min_cloud_frac_mult) > (cloud_frac_zm(i - 1, k - 1))
                                                  ? (min_cloud_frac_mult)
                                                  : (cloud_frac_zm(i - 1, k - 1)));
            C2thl_1d(i - 1, k - 1) =
                clubb_params(i - 1, iC2thl) * ((min_cloud_frac_mult) > (cloud_frac_zm(i - 1, k - 1))
                                                   ? (min_cloud_frac_mult)
                                                   : (cloud_frac_zm(i - 1, k - 1)));
            C2rtthl_1d(i - 1, k - 1) = clubb_params(i - 1, iC2rtthl) *
                                       ((min_cloud_frac_mult) > (cloud_frac_zm(i - 1, k - 1))
                                            ? (min_cloud_frac_mult)
                                            : (cloud_frac_zm(i - 1, k - 1)));
          } else {
            C2rt_1d(i - 1, k - 1)    = clubb_params(i - 1, iC2rt);
            C2thl_1d(i - 1, k - 1)   = clubb_params(i - 1, iC2thl);
            C2rtthl_1d(i - 1, k - 1) = clubb_params(i - 1, iC2rtthl);
          }
        });
  } else {
    Kokkos::parallel_for(
        "advance_xp2_xpyp::loop2", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          C2rt_1d(i - 1, k - 1)    = clubb_params(i - 1, iC2rt);
          C2thl_1d(i - 1, k - 1)   = clubb_params(i - 1, iC2thl);
          C2rtthl_1d(i - 1, k - 1) = clubb_params(i - 1, iC2rtthl);
        });
  }

  Kokkos::parallel_for(
      "advance_xp2_xpyp::loop3", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzm);
        remaining /= (nzm);
        C2sclr_1d(i - 1, k - 1) = clubb_params(i - 1, iC2rt);
        C4_1d(i - 1, k - 1)     = two_thirds * clubb_params(i - 1, iC4);
        C14_1d(i - 1, k - 1)    = one_third * clubb_params(i - 1, iC14);
      });

  // Are we solving for passive scalars as well?
  if (sclr_dim > 0) {
    l_scalar_calc = true;
  } else {
    l_scalar_calc = false;
  }

  // Define the Coefficent of Eddy Diffusivity for the variances
  // and covariances.
  Kokkos::parallel_for(
      "advance_xp2_xpyp::loop4", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzt);
        remaining /= (nzt);
        // Kw2 is used for variances and covariances rtp2, thlp2, rtpthlp, and
        // passive scalars.  The variances and covariances are located on the
        // momentum levels.  Kw2 is located on the thermodynamic levels.
        // Kw2 = c_K2 * Kh_zt
        Kw2(i - 1, k - 1) = clubb_params(i - 1, ic_K2) * Kh_zt(i - 1, k - 1);

        // Kw9 is used for variances up2 and vp2.  The variances are located on
        // the momentum levels.  Kw9 is located on the thermodynamic levels.
        // Kw9 = c_K9 * Kh_zt
        Kw9(i - 1, k - 1) = clubb_params(i - 1, ic_K9) * Kh_zt(i - 1, k - 1);
      });

  Kokkos::deep_copy(Kw2_zm, zt2zm_api(nzm, nzt, ngrdcol, gr, Kw2, zero));
  Kokkos::deep_copy(Kw9_zm, zt2zm_api(nzm, nzt, ngrdcol, gr, Kw9, zero));

  if (l_lmm_stepping) {

    Kokkos::parallel_for(
        "advance_xp2_xpyp::loop5", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          thlp2_old(i - 1, k - 1)   = thlp2(i - 1, k - 1);
          rtp2_old(i - 1, k - 1)    = rtp2(i - 1, k - 1);
          rtpthlp_old(i - 1, k - 1) = rtpthlp(i - 1, k - 1);
        });

    if (sclr_dim > 0) {
      Kokkos::parallel_for(
          "advance_xp2_xpyp::loop6",
          std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            const int sclr = 1 + remaining % (sclr_dim);
            remaining /= (sclr_dim);
            sclrp2_old(i - 1, k - 1, sclr - 1)    = sclrp2(i - 1, k - 1, sclr - 1);
            sclrprtp_old(i - 1, k - 1, sclr - 1)  = sclrprtp(i - 1, k - 1, sclr - 1);
            sclrpthlp_old(i - 1, k - 1, sclr - 1) = sclrpthlp(i - 1, k - 1, sclr - 1);
          });
    }
  }

  // Calculate all the explicit and implicit turbulent advection terms
  calc_xp2_xpyp_ta_terms(
      nzm, nzt, ngrdcol, sclr_dim, gr, wprtp, wprtp2, wpthlp, wpthlp2, wprtpthlp, rtp2, // In
      thlp2, rtpthlp, upwp, vpwp, up2, vp2, wp2, wp2_zt, wpsclrp, wpsclrp2, wpsclrprtp, // In
      wpsclrpthlp, sclrp2, sclrprtp, sclrpthlp, rho_ds_zt, invrs_rho_ds_zm, rho_ds_zm,  // In
      wp3_on_wp2, wp3_on_wp2_zt, sigma_sqd_w, pdf_implicit_coefs_terms, l_scalar_calc,  // In
      Kokkos::subview(clubb_params, Kokkos::ALL, ibeta), iiPDF_type, l_upwind_xpyp_ta,  // In
      l_godunov_upwind_xpyp_ta,                                                         // In
      stats,                                                                            // InOut
      lhs_ta_wprtp2, lhs_ta_wpthlp2, lhs_ta_wprtpthlp, lhs_ta_wpup2, lhs_ta_wpvp2,      // Out
      lhs_ta_wpsclrp2, lhs_ta_wprtpsclrp, lhs_ta_wpthlpsclrp, rhs_ta_wprtp2,            // Out
      rhs_ta_wpthlp2, rhs_ta_wprtpthlp, rhs_ta_wpup2, rhs_ta_wpvp2, rhs_ta_wpsclrp2,    // Out
      rhs_ta_wprtpsclrp, rhs_ta_wpthlpsclrp);                                           // Out

  // Calculate LHS eddy diffusion term: dissipation term 2 (dp2). This is the
  // diffusion term for all LHS matrices except <w'u'^2> and <w'v'^2>
  f2k::diffusion_zm_lhs(nzm, nzt, ngrdcol, gr, Kw2, Kw2_zm, nu_vert_res_dep.nu2,
                        invrs_rho_ds_zm, // In
                        rho_ds_zt,       // In
                        lhs_diff);       // Out

  // Calculate LHS mean advection (ma) term, this term is equal for all LHS matrices
  f2k::term_ma_zm_lhs(nzm, nzt, ngrdcol, wm_zm, gr.invrs_dzm, gr.weights_zm2zt, // In
                      lhs_ma);                                                  // Out

  // Check if we can use a single LHS to solve
  if (l_explicit_turbulent_adv_xpyp || !l_explicit_turbulent_adv_xpyp && iiPDF_type == iiPDF_ADG1) {

    l_single_solve_possible = true;

    // Native logical AND reduction: reduce integer disagreement hits.
    int unequal_C2_columns = 0;
    Kokkos::parallel_reduce(
        "advance_xp2_xpyp::single_lhs", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& unequal) {
          const real_t C2rt    = clubb_params(column, iC2rt);
          const real_t C2thl   = clubb_params(column, iC2thl);
          const real_t C2rtthl = clubb_params(column, iC2rtthl);
          if (Kokkos::abs(C2rt - C2thl) > Kokkos::abs(C2rt + C2thl) / 2 * eps ||
              Kokkos::abs(C2rt - C2rtthl) > Kokkos::abs(C2rt + C2rtthl) / 2 * eps)
            ++unequal;
        },
        unequal_C2_columns);
    l_single_solve_possible = unequal_C2_columns == 0;

  } else {

    l_single_solve_possible = false;
  }

  if (l_single_solve_possible) {

    // All left hand side matricies are equal for rtp2, thlp2, rtpthlp, and scalars.
    // Thus only one solve is neccesary, using combined right hand sides
    solve_xp2_xpyp_with_single_lhs(
        nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, sclr_idx, C2rt_1d, invrs_tau_xp2_zm, // In
        rtm, thlm, wprtp, wpthlp, rtp2_forcing, thlp2_forcing, rtpthlp_forcing, sclrm,  // In
        wpsclrp, lhs_ta_wprtp2, lhs_ma, lhs_diff, rhs_ta_wprtp2, rhs_ta_wpthlp2,        // In
        rhs_ta_wprtpthlp, rhs_ta_wpsclrp2, rhs_ta_wprtpsclrp, rhs_ta_wpthlpsclrp, dt,   // In
        l_scalar_calc, l_lmm_stepping, l_implemented, tridiag_solve_method,             // In
        stats, rtp2, thlp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp, err_info);            // InOut
  } else {

    // Left hand sides are potentially different, this requires multiple solves
    solve_xp2_xpyp_with_multiple_lhs(
        nzm, nzt, ngrdcol, sclr_dim, sclr_tol, gr, sclr_idx, C2rt_1d, C2thl_1d,           // In
        C2rtthl_1d, C2sclr_1d, invrs_tau_xp2_zm, rtm, thlm, wprtp, wpthlp, rtp2_forcing,  // In
        thlp2_forcing, rtpthlp_forcing, sclrm, wpsclrp, lhs_ta_wprtp2, lhs_ta_wpthlp2,    // In
        lhs_ta_wprtpthlp, lhs_ta_wpsclrp2, lhs_ta_wprtpsclrp, lhs_ta_wpthlpsclrp, lhs_ma, // In
        lhs_diff, rhs_ta_wprtp2, rhs_ta_wpthlp2, rhs_ta_wprtpthlp, rhs_ta_wpsclrp2,       // In
        rhs_ta_wprtpsclrp, rhs_ta_wpthlpsclrp, dt, iiPDF_type, l_scalar_calc,             // In
        l_lmm_stepping, l_implemented, tridiag_solve_method,                              // In
        stats, rtp2, thlp2, rtpthlp, sclrp2, sclrprtp, sclrpthlp, err_info);              // InOut
  }

  //!!!!***** u'^2 / v'^2 *****!!!!!

  // Calculate LHS eddy diffusion term: dissipation term 2 (dp2), for <w'u'^2> and <w'v'^2>
  f2k::diffusion_zm_lhs(nzm, nzt, ngrdcol, gr, Kw9, Kw9_zm, nu_vert_res_dep.nu9,
                        invrs_rho_ds_zm, // In
                        rho_ds_zt,       // In
                        lhs_diff_uv);    // Out

  if (l_lmm_stepping) {
    Kokkos::parallel_for(
        "advance_xp2_xpyp::loop7", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          up2_old(i - 1, k - 1) = up2(i - 1, k - 1);
          vp2_old(i - 1, k - 1) = vp2(i - 1, k - 1);
        });
  }

  // Dissipation term 1 using C4
  term_dp1_lhs(nzm, ngrdcol, gr, C4_1d, invrs_tau_C4_zm, // In
               lhs_dp1_C4);                              // Out

  // Dissipation term 1 using C14
  term_dp1_lhs(nzm, ngrdcol, gr, C14_1d, invrs_tau_C14_zm, // In
               lhs_dp1_C14);                               // Out

  Kokkos::parallel_for(
      "advance_xp2_xpyp::loop8", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        lhs_dp1(i - 1, k - 1) =
            (lhs_dp1_C4(i - 1, k - 1) + lhs_dp1_C14(i - 1, k - 1)) * gamma_over_implicit_ts;
      });

  if (iiPDF_type == iiPDF_new_hybrid) {

    // Different LHS required for up2 and vp2.

    // Solve for up2

    // Implicit contributions to term up2
    xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wpup2, lhs_ma, lhs_diff_uv, lhs_dp1, // In
                 lhs);                                                             // Out

    // Explicit contributions to up2
    xp2_xpyp_uv_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_up2, dt, wp2, wpthvp, invrs_tau_C4_zm, // In
                    invrs_tau_C14_zm, um, vm, upwp, vpwp, up2, vp2, thv_ds_zm,             // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC4),                       // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_shr),                 // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_buoy),                // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC14), lhs_splat_wp2,
                    lhs_ta_wpup2,                                                // In
                    rhs_ta_wpup2, lhs_dp1_C4, lhs_dp1_C14,                       // In
                    stats,                                                       // InOut
                    Kokkos::subview(uv_rhs, Kokkos::ALL, Kokkos::ALL, (1 - 1))); // Out

    if (l_ho_nontrad_coriolis) {

      // Add the nontraditional Coriolis term
      // Hing Ong, 19 July 2025
      // Vertical compression of eddies.
      Kokkos::parallel_for(
          "advance_xp2_xpyp::loop9", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            uv_rhs(i - 1, k - 1, 1 - 1) =
                uv_rhs(i - 1, k - 1, 1 - 1) - two * fcor_y(i - 1) * upwp(i - 1, k - 1);
          });

      if (stats.flags().l_sample) {
        Kokkos::parallel_for(
            "advance_xp2_xpyp::array10",
            ((stats_tmp_zm.extent(0) - (1)) / (1) + 1) * ((stats_tmp_zm.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((stats_tmp_zm.extent(0) - (1)) / (1) + 1);
              remaining /= ((stats_tmp_zm.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((stats_tmp_zm.extent(1) - (1)) / (1) + 1);
              remaining /= ((stats_tmp_zm.extent(1) - (1)) / (1) + 1);
              stats_tmp_zm(slice_0, slice_1) = zero;
            });
        Kokkos::parallel_for(
            "advance_xp2_xpyp::loop11", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              stats_tmp_zm(i - 1, k - 1) = -two * fcor_y(i - 1) * upwp(i - 1, k - 1);
            });

        // u'u' term nct is completely explicit; call stat_update_var.
        // Hing Ong, 22 July 2025
        f2k::stats_update("up2_nct", stats_tmp_zm, // In
                          stats);                  // InOut
      }
    }

    // Solve the tridiagonal system
    xp2_xpyp_solve(nzm, ngrdcol, xp2_xpyp_up2_vp2, 1, l_implemented, gr, tridiag_solve_method, // In
                   stats, uv_rhs, lhs, err_info, // InOut
                   uv_solution);                 // Out

    Kokkos::parallel_for(
        "advance_xp2_xpyp::array12",
        ((up2.extent(0) - (1)) / (1) + 1) * ((up2.extent(1) - (1)) / (1) + 1),
        KOKKOS_LAMBDA(int cell) {
          int remaining     = cell;
          const int slice_0 = remaining % ((up2.extent(0) - (1)) / (1) + 1);
          remaining /= ((up2.extent(0) - (1)) / (1) + 1);
          const int slice_1 = remaining % ((up2.extent(1) - (1)) / (1) + 1);
          remaining /= ((up2.extent(1) - (1)) / (1) + 1);
          up2(slice_0, slice_1) =
              uv_solution((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (1 - 1));
        });

    if (l_lmm_stepping) {
      Kokkos::parallel_for(
          "advance_xp2_xpyp::array13",
          ((up2.extent(0) - (1)) / (1) + 1) * ((up2.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((up2.extent(0) - (1)) / (1) + 1);
            remaining /= ((up2.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((up2.extent(1) - (1)) / (1) + 1);
            remaining /= ((up2.extent(1) - (1)) / (1) + 1);
            up2(slice_0, slice_1) = one_half * (up2_old(slice_0, slice_1) + up2(slice_0, slice_1));
          });
    }

    if (stats.flags().l_sample) {
      {
        Real2D stats_lhs_14(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_lhs_14"),
            lhs_dp1_C14.extent(0), lhs_dp1_C14.extent(1));
        Kokkos::parallel_for(
            "advance_xp2_xpyp::stats_lhs_14", lhs_dp1_C14.size(), KOKKOS_LAMBDA(const int cell) {
              const int column = cell % lhs_dp1_C14.extent(0), level = cell / lhs_dp1_C14.extent(0);
              stats_lhs_14(column, level) = gamma_over_implicit_ts * lhs_dp1_C14(column, level);
            });
        Real2D stats_lhs_15(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_lhs_15"),
            lhs_dp1_C4.extent(0), lhs_dp1_C4.extent(1));
        Kokkos::parallel_for(
            "advance_xp2_xpyp::stats_lhs_15", lhs_dp1_C4.size(), KOKKOS_LAMBDA(const int cell) {
              const int column = cell % lhs_dp1_C4.extent(0), level = cell / lhs_dp1_C4.extent(0);
              stats_lhs_15(column, level) = gamma_over_implicit_ts * lhs_dp1_C4(column, level);
            });
        stats_finalize_xp2_xpyp_terms(nzm, ngrdcol, xp2_xpyp_up2, gr, up2, stats_lhs_14,
                                      stats_lhs_15, lhs_diff_uv, // In
                                      lhs_ta_wpup2, lhs_ma,      // In
                                      stats);                    // InOut
      }
    }

    // Implicit contributions to term vp2
    xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wpvp2, lhs_ma, lhs_diff_uv, lhs_dp1, // In
                 lhs);                                                             // Out

    // Explicit contributions to vp2
    xp2_xpyp_uv_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_vp2, dt, wp2, wpthvp, invrs_tau_C4_zm, // In
                    invrs_tau_C14_zm, vm, um, vpwp, upwp, vp2, up2, thv_ds_zm,             // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC4),                       // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_shr),                 // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_buoy),                // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC14), lhs_splat_wp2,
                    lhs_ta_wpvp2,                                                // In
                    rhs_ta_wpvp2, lhs_dp1_C4, lhs_dp1_C14,                       // In
                    stats,                                                       // InOut
                    Kokkos::subview(uv_rhs, Kokkos::ALL, Kokkos::ALL, (1 - 1))); // Out

    // Solve the tridiagonal system
    xp2_xpyp_solve(nzm, ngrdcol, xp2_xpyp_up2_vp2, 1, l_implemented, gr, tridiag_solve_method, // In
                   stats, uv_rhs, lhs, err_info, // InOut
                   uv_solution);                 // Out

    Kokkos::parallel_for(
        "advance_xp2_xpyp::array16",
        ((vp2.extent(0) - (1)) / (1) + 1) * ((vp2.extent(1) - (1)) / (1) + 1),
        KOKKOS_LAMBDA(int cell) {
          int remaining     = cell;
          const int slice_0 = remaining % ((vp2.extent(0) - (1)) / (1) + 1);
          remaining /= ((vp2.extent(0) - (1)) / (1) + 1);
          const int slice_1 = remaining % ((vp2.extent(1) - (1)) / (1) + 1);
          remaining /= ((vp2.extent(1) - (1)) / (1) + 1);
          vp2(slice_0, slice_1) =
              uv_solution((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (1 - 1));
        });

    if (l_lmm_stepping) {
      Kokkos::parallel_for(
          "advance_xp2_xpyp::array17",
          ((vp2.extent(0) - (1)) / (1) + 1) * ((vp2.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((vp2.extent(0) - (1)) / (1) + 1);
            remaining /= ((vp2.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((vp2.extent(1) - (1)) / (1) + 1);
            remaining /= ((vp2.extent(1) - (1)) / (1) + 1);
            vp2(slice_0, slice_1) = one_half * (vp2_old(slice_0, slice_1) + vp2(slice_0, slice_1));
          });
    }

    if (stats.flags().l_sample) {
      {
        Real2D stats_lhs_18(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_lhs_18"),
            lhs_dp1_C14.extent(0), lhs_dp1_C14.extent(1));
        Kokkos::parallel_for(
            "advance_xp2_xpyp::stats_lhs_18", lhs_dp1_C14.size(), KOKKOS_LAMBDA(const int cell) {
              const int column = cell % lhs_dp1_C14.extent(0), level = cell / lhs_dp1_C14.extent(0);
              stats_lhs_18(column, level) = gamma_over_implicit_ts * lhs_dp1_C14(column, level);
            });
        Real2D stats_lhs_19(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_lhs_19"),
            lhs_dp1_C4.extent(0), lhs_dp1_C4.extent(1));
        Kokkos::parallel_for(
            "advance_xp2_xpyp::stats_lhs_19", lhs_dp1_C4.size(), KOKKOS_LAMBDA(const int cell) {
              const int column = cell % lhs_dp1_C4.extent(0), level = cell / lhs_dp1_C4.extent(0);
              stats_lhs_19(column, level) = gamma_over_implicit_ts * lhs_dp1_C4(column, level);
            });
        stats_finalize_xp2_xpyp_terms(nzm, ngrdcol, xp2_xpyp_vp2, gr, vp2, stats_lhs_18,
                                      stats_lhs_19, lhs_diff_uv, // In
                                      lhs_ta_wpvp2, lhs_ma,      // In
                                      stats);                    // InOut
      }
    }

  } else {

    // ADG1 allows up2 and vp2 to use the same LHS.

    // Implicit contributions to term up2/vp2
    xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wpup2, lhs_ma, lhs_diff_uv, lhs_dp1, // In
                 lhs);                                                             // Out

    // Explicit contributions to up2
    xp2_xpyp_uv_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_up2, dt, wp2, wpthvp, invrs_tau_C4_zm, // In
                    invrs_tau_C14_zm, um, vm, upwp, vpwp, up2, vp2, thv_ds_zm,             // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC4),                       // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_shr),                 // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_buoy),                // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC14), lhs_splat_wp2,
                    lhs_ta_wpup2,                                                // In
                    rhs_ta_wpup2, lhs_dp1_C4, lhs_dp1_C14,                       // In
                    stats,                                                       // InOut
                    Kokkos::subview(uv_rhs, Kokkos::ALL, Kokkos::ALL, (1 - 1))); // Out

    if (l_ho_nontrad_coriolis) {

      // Add the nontraditional Coriolis term
      // Hing Ong, 19 July 2025
      // Vertical compression of eddies.
      Kokkos::parallel_for(
          "advance_xp2_xpyp::loop20", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            uv_rhs(i - 1, k - 1, 1 - 1) =
                uv_rhs(i - 1, k - 1, 1 - 1) - two * fcor_y(i - 1) * upwp(i - 1, k - 1);
          });

      if (stats.flags().l_sample) {
        Kokkos::parallel_for(
            "advance_xp2_xpyp::array21",
            ((stats_tmp_zm.extent(0) - (1)) / (1) + 1) * ((stats_tmp_zm.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((stats_tmp_zm.extent(0) - (1)) / (1) + 1);
              remaining /= ((stats_tmp_zm.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((stats_tmp_zm.extent(1) - (1)) / (1) + 1);
              remaining /= ((stats_tmp_zm.extent(1) - (1)) / (1) + 1);
              stats_tmp_zm(slice_0, slice_1) = zero;
            });
        Kokkos::parallel_for(
            "advance_xp2_xpyp::loop22", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              stats_tmp_zm(i - 1, k - 1) = -two * fcor_y(i - 1) * upwp(i - 1, k - 1);
            });

        // u'u' term nct is completely explicit; call stat_update_var.
        // Hing Ong, 22 July 2025
        f2k::stats_update("up2_nct", stats_tmp_zm, // In
                          stats);                  // InOut
      }
    }

    // Explicit contributions to vp2
    xp2_xpyp_uv_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_vp2, dt, wp2, wpthvp, invrs_tau_C4_zm, // In
                    invrs_tau_C14_zm, vm, um, vpwp, upwp, vp2, up2, thv_ds_zm,             // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC4),                       // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_shr),                 // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC_uu_buoy),                // In
                    Kokkos::subview(clubb_params, Kokkos::ALL, iC14), lhs_splat_wp2,
                    lhs_ta_wpup2,                                                // In
                    rhs_ta_wpvp2, lhs_dp1_C4, lhs_dp1_C14,                       // In
                    stats,                                                       // InOut
                    Kokkos::subview(uv_rhs, Kokkos::ALL, Kokkos::ALL, (2 - 1))); // Out

    // Solve the tridiagonal system
    xp2_xpyp_solve(nzm, ngrdcol, xp2_xpyp_up2_vp2, 2, l_implemented, gr, tridiag_solve_method, // In
                   stats, uv_rhs, lhs, err_info, // InOut
                   uv_solution);                 // Out

    Kokkos::parallel_for(
        "advance_xp2_xpyp::loop23", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          up2(i - 1, k - 1) = uv_solution(i - 1, k - 1, 1 - 1);
          vp2(i - 1, k - 1) = uv_solution(i - 1, k - 1, 2 - 1);
        });

    if (l_lmm_stepping) {
      Kokkos::parallel_for(
          "advance_xp2_xpyp::loop24", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            up2(i - 1, k - 1) = one_half * (up2_old(i - 1, k - 1) + up2(i - 1, k - 1));
            vp2(i - 1, k - 1) = one_half * (vp2_old(i - 1, k - 1) + vp2(i - 1, k - 1));
          });
    }

    if (stats.flags().l_sample) {
      {
        Real2D stats_lhs_25(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_lhs_25"),
            lhs_dp1_C14.extent(0), lhs_dp1_C14.extent(1));
        Kokkos::parallel_for(
            "advance_xp2_xpyp::stats_lhs_25", lhs_dp1_C14.size(), KOKKOS_LAMBDA(const int cell) {
              const int column = cell % lhs_dp1_C14.extent(0), level = cell / lhs_dp1_C14.extent(0);
              stats_lhs_25(column, level) = gamma_over_implicit_ts * lhs_dp1_C14(column, level);
            });
        Real2D stats_lhs_26(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_lhs_26"),
            lhs_dp1_C4.extent(0), lhs_dp1_C4.extent(1));
        Kokkos::parallel_for(
            "advance_xp2_xpyp::stats_lhs_26", lhs_dp1_C4.size(), KOKKOS_LAMBDA(const int cell) {
              const int column = cell % lhs_dp1_C4.extent(0), level = cell / lhs_dp1_C4.extent(0);
              stats_lhs_26(column, level) = gamma_over_implicit_ts * lhs_dp1_C4(column, level);
            });
        stats_finalize_xp2_xpyp_terms(nzm, ngrdcol, xp2_xpyp_up2, gr, up2, stats_lhs_25,
                                      stats_lhs_26, lhs_diff_uv, // In
                                      lhs_ta_wpup2, lhs_ma,      // In
                                      stats);                    // InOut
      }
      {
        Real2D stats_lhs_27(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_lhs_27"),
            lhs_dp1_C14.extent(0), lhs_dp1_C14.extent(1));
        Kokkos::parallel_for(
            "advance_xp2_xpyp::stats_lhs_27", lhs_dp1_C14.size(), KOKKOS_LAMBDA(const int cell) {
              const int column = cell % lhs_dp1_C14.extent(0), level = cell / lhs_dp1_C14.extent(0);
              stats_lhs_27(column, level) = gamma_over_implicit_ts * lhs_dp1_C14(column, level);
            });
        Real2D stats_lhs_28(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_lhs_28"),
            lhs_dp1_C4.extent(0), lhs_dp1_C4.extent(1));
        Kokkos::parallel_for(
            "advance_xp2_xpyp::stats_lhs_28", lhs_dp1_C4.size(), KOKKOS_LAMBDA(const int cell) {
              const int column = cell % lhs_dp1_C4.extent(0), level = cell / lhs_dp1_C4.extent(0);
              stats_lhs_28(column, level) = gamma_over_implicit_ts * lhs_dp1_C4(column, level);
            });
        stats_finalize_xp2_xpyp_terms(nzm, ngrdcol, xp2_xpyp_vp2, gr, vp2, stats_lhs_27,
                                      stats_lhs_28, lhs_diff_uv, // In
                                      lhs_ta_wpup2, lhs_ma,      // In
                                      stats);                    // InOut
      }
    }
  }

  // Apply the positive definite scheme to variances
  if (fill_holes_type != 0) {
    pos_definite_variances(nzm, ngrdcol, gr, xp2_xpyp_rtp2, fill_holes_type, dt,
                           (rt_tol * rt_tol), // In
                           rho_ds_zm,         // In
                           stats, rtp2);      // InOut

    pos_definite_variances(nzm, ngrdcol, gr, xp2_xpyp_thlp2, fill_holes_type, dt,
                           (thl_tol * thl_tol), // In
                           rho_ds_zm,           // In
                           stats, thlp2);       // InOut

    pos_definite_variances(nzm, ngrdcol, gr, xp2_xpyp_up2, fill_holes_type, dt, w_tol_sqd,
                           rho_ds_zm,   // In
                           stats, up2); // InOut

    pos_definite_variances(nzm, ngrdcol, gr, xp2_xpyp_vp2, fill_holes_type, dt, w_tol_sqd,
                           rho_ds_zm,   // In
                           stats, vp2); // InOut
  }

  // Clipping for r_t'^2

  // threshold = zero_threshold
  //
  // where ( wp2 >= w_tol_sqd ) &
  //    threshold = rt_tol*rt_tol

  // The value of rtp2 is not allowed to become smaller than the threshold
  // value of rt_tol^2.  Additionally, that threshold value may be boosted at
  // any grid level in order to keep the overall correlation of w and rt
  // between the values of -max_mag_correlation_flux and
  // max_mag_correlation_flux by boosting rtp2 rather than by limiting the
  // magnitude of wprtp.
  if (l_min_xp2_from_corr_wx) {

    // The overall correlation of w and rt is:
    //
    // corr_w_rt = wprtp / ( sqrt( wp2 ) * sqrt( rtp2 ) ).
    //
    // Squaring both sides, the equation becomes:
    //
    // corr_w_rt^2 = wprtp^2 / ( wp2 * rtp2 ).
    //
    // Using max_mag_correlation_flux for the correlation and then solving for
    // the minimum of rtp2, the equation becomes:
    //
    // rtp2|_min = wprtp^2 / ( wp2 * max_mag_correlation_flux^2 ).
    Kokkos::parallel_for(
        "advance_xp2_xpyp::loop29", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);

          threshold_array(i - 1, k - 1) =
              (((rt_tol * rt_tol)) >
                       ((wprtp(i - 1, k - 1) * wprtp(i - 1, k - 1)) /
                        (wp2(i - 1, k - 1) * (max_mag_correlation_flux * max_mag_correlation_flux)))
                   ? ((rt_tol * rt_tol))
                   : ((wprtp(i - 1, k - 1) * wprtp(i - 1, k - 1)) /
                      (wp2(i - 1, k - 1) * (max_mag_correlation_flux * max_mag_correlation_flux))));
        });

    f2k::clip_variance(nzm, ngrdcol, gr, xp2_xpyp_rtp2, dt, threshold_array, // In
                       stats, rtp2);                                         // InOut
  } else {

    // Consider only the minimum tolerance threshold value for rtp2.
    Kokkos::parallel_for(
        "advance_xp2_xpyp::loop30", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          threshold_array(i - 1, k - 1) = (rt_tol * rt_tol);
        });

    f2k::clip_variance(nzm, ngrdcol, gr, xp2_xpyp_rtp2, dt, threshold_array, // In
                       stats, rtp2);                                         // InOut
  }

  // Special clipping on the variance of rt to prevent a large variance at
  // higher altitudes.  This is done because we don't want the PDF to extend
  // into the negative, and found that for latin hypercube sampling a large
  // variance aloft leads to negative samples of total water.
  // -dschanen 8 Dec 2010
  if (l_clip_large_rtp2) {
    if (stats.flags().l_sample) {
      // This overwrites stats clipping data from clip_variance
      // Store previous value in order to calculate clipping
      {
        Real2D stats_values_31(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_values_31"),
            rtp2.extent(0), rtp2.extent(1));
        Kokkos::parallel_for(
            "advance_xp2_xpyp::stats_value", rtp2.size(), KOKKOS_LAMBDA(int cell) {
              const int i           = cell % rtp2.extent(0);
              const int k           = cell / rtp2.extent(0);
              stats_values_31(i, k) = rtp2(i, k) / dt;
            });
        f2k::stats_begin_budget("rtp2_cl", stats_values_31, stats);
      }
    }

    Kokkos::deep_copy(rtm_zm, zt2zm_api(nzm, nzt, ngrdcol, gr, rtm, zero_threshold));

    Kokkos::parallel_for(
        "advance_xp2_xpyp::loop32", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          real_t threshold;
          threshold = (((rt_tol * rt_tol)) > (clubb_params(i - 1, irtp2_clip_coef) *
                                              (rtm_zm(i - 1, k - 1) * rtm_zm(i - 1, k - 1)))
                           ? ((rt_tol * rt_tol))
                           : (clubb_params(i - 1, irtp2_clip_coef) *
                              (rtm_zm(i - 1, k - 1) * rtm_zm(i - 1, k - 1))));
          if (rtp2(i - 1, k - 1) > threshold) {
            rtp2(i - 1, k - 1) = threshold;
          }
        });

    if (stats.flags().l_sample) {
      // Store final value in order to calculate clipping
      {
        Real2D stats_values_33(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_values_33"),
            rtp2.extent(0), rtp2.extent(1));
        Kokkos::parallel_for(
            "advance_xp2_xpyp::stats_value", rtp2.size(), KOKKOS_LAMBDA(int cell) {
              const int i           = cell % rtp2.extent(0);
              const int k           = cell / rtp2.extent(0);
              stats_values_33(i, k) = rtp2(i, k) / dt;
            });
        f2k::stats_finalize_budget("rtp2_cl", stats_values_33, stats, false);
      }
    }
  }

  // Clipping for th_l'^2

  // threshold = zero_threshold
  //
  // where ( wp2 >= w_tol_sqd ) &
  //    threshold = thl_tol*thl_tol

  // The value of thlp2 is not allowed to become smaller than the threshold
  // value of thl_tol^2.  Additionally, that threshold value may be boosted at
  // any grid level in order to keep the overall correlation of w and theta-l
  // between the values of -max_mag_correlation_flux and
  // max_mag_correlation_flux by boosting thlp2 rather than by limiting the
  // magnitude of wpthlp.
  if (l_min_xp2_from_corr_wx) {

    // The overall correlation of w and theta-l is:
    //
    // corr_w_thl = wpthlp / ( sqrt( wp2 ) * sqrt( thlp2 ) ).
    //
    // Squaring both sides, the equation becomes:
    //
    // corr_w_thl^2 = wpthlp^2 / ( wp2 * thlp2 ).
    //
    // Using max_mag_correlation_flux for the correlation and then solving for
    // the minimum of thlp2, the equation becomes:
    //
    // thlp2|_min = wpthlp^2 / ( wp2 * max_mag_correlation_flux^2 ).
    Kokkos::parallel_for(
        "advance_xp2_xpyp::loop34", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          threshold_array(i - 1, k - 1) =
              (((thl_tol * thl_tol)) >
                       ((wpthlp(i - 1, k - 1) * wpthlp(i - 1, k - 1)) /
                        (wp2(i - 1, k - 1) * (max_mag_correlation_flux * max_mag_correlation_flux)))
                   ? ((thl_tol * thl_tol))
                   : ((wpthlp(i - 1, k - 1) * wpthlp(i - 1, k - 1)) /
                      (wp2(i - 1, k - 1) * (max_mag_correlation_flux * max_mag_correlation_flux))));
        });

    f2k::clip_variance(nzm, ngrdcol, gr, xp2_xpyp_thlp2, dt, threshold_array, // In
                       stats, thlp2);                                         // InOut

  } else {

    // Consider only the minimum tolerance threshold value for thlp2.
    Kokkos::parallel_for(
        "advance_xp2_xpyp::loop35", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          threshold_array(i - 1, k - 1) = (thl_tol * thl_tol);
        });

    f2k::clip_variance(nzm, ngrdcol, gr, xp2_xpyp_thlp2, dt, threshold_array, // In
                       stats, thlp2);                                         // InOut
  }

  // Clipping for u'^2

  // Clip negative values of up2
  Kokkos::parallel_for(
      "advance_xp2_xpyp::loop36", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzm);
        remaining /= (nzm);
        threshold_array(i - 1, k - 1) = w_tol_sqd;
      });

  f2k::clip_variance(nzm, ngrdcol, gr, xp2_xpyp_up2, dt, threshold_array, // In
                     stats, up2);                                         // InOut

  // Clip excessively large values of up2

  if (stats.flags().l_sample) {
    {
      Real2D stats_values_37(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_values_37"),
          up2.extent(0), up2.extent(1));
      Kokkos::parallel_for(
          "advance_xp2_xpyp::stats_value", up2.size(), KOKKOS_LAMBDA(int cell) {
            const int i           = cell % up2.extent(0);
            const int k           = cell / up2.extent(0);
            stats_values_37(i, k) = up2(i, k) / dt;
          });
      f2k::stats_begin_budget("up2_cl", stats_values_37, stats);
    }
  }

  Kokkos::parallel_for(
      "advance_xp2_xpyp::loop38", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzm);
        remaining /= (nzm);
        up2(i - 1, k - 1) =
            ((up2(i - 1, k - 1)) < (real_t(1000.)) ? (up2(i - 1, k - 1)) : (real_t(1000.)));
      });

  if (stats.flags().l_sample) {
    {
      Real2D stats_values_39(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_values_39"),
          up2.extent(0), up2.extent(1));
      Kokkos::parallel_for(
          "advance_xp2_xpyp::stats_value", up2.size(), KOKKOS_LAMBDA(int cell) {
            const int i           = cell % up2.extent(0);
            const int k           = cell / up2.extent(0);
            stats_values_39(i, k) = up2(i, k) / dt;
          });
      f2k::stats_finalize_budget("up2_cl", stats_values_39, stats, false);
    }
  }

  // Clipping for v'^2

  // Clip negative values of vp2
  Kokkos::parallel_for(
      "advance_xp2_xpyp::loop40", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzm);
        remaining /= (nzm);
        threshold_array(i - 1, k - 1) = w_tol_sqd;
      });

  f2k::clip_variance(nzm, ngrdcol, gr, xp2_xpyp_vp2, dt, threshold_array, // In
                     stats, vp2);                                         // InOut

  if (stats.flags().l_sample) {
    {
      Real2D stats_values_41(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_values_41"),
          vp2.extent(0), vp2.extent(1));
      Kokkos::parallel_for(
          "advance_xp2_xpyp::stats_value", vp2.size(), KOKKOS_LAMBDA(int cell) {
            const int i           = cell % vp2.extent(0);
            const int k           = cell / vp2.extent(0);
            stats_values_41(i, k) = vp2(i, k) / dt;
          });
      f2k::stats_begin_budget("vp2_cl", stats_values_41, stats);
    }
  }

  Kokkos::parallel_for(
      "advance_xp2_xpyp::loop42", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzm);
        remaining /= (nzm);
        vp2(i - 1, k - 1) =
            ((vp2(i - 1, k - 1)) < (real_t(1000.)) ? (vp2(i - 1, k - 1)) : (real_t(1000.)));
      });

  if (stats.flags().l_sample) {
    {
      Real2D stats_values_43(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_values_43"),
          vp2.extent(0), vp2.extent(1));
      Kokkos::parallel_for(
          "advance_xp2_xpyp::stats_value", vp2.size(), KOKKOS_LAMBDA(int cell) {
            const int i           = cell % vp2.extent(0);
            const int k           = cell / vp2.extent(0);
            stats_values_43(i, k) = vp2(i, k) / dt;
          });
      f2k::stats_finalize_budget("vp2_cl", stats_values_43, stats, false);
    }
  }

  // When selected, apply sponge damping after up2 and vp2 have been advanced.
  // Native up2_vp2_sponge_damp_settings%l_sponge_damping:
  // Sponge damping is excluded with the JAX scope. Integration must reject
  // sponge-enabled cases upstream; see INTEGRATION.md.

  // Clipping for r_t'th_l'
  // Clipping r_t'th_l' at each vertical level, based on the
  // correlation of r_t and th_l at each vertical level, such that:
  // corr_(r_t,th_l) = r_t'th_l' / [ sqrt(r_t'^2) * sqrt(th_l'^2) ];
  // -1 <= corr_(r_t,th_l) <= 1.
  // Since r_t'^2, th_l'^2, and r_t'th_l' are all computed in the
  // same place, clipping for r_t'th_l' only has to be done once.
  if (stats.flags().l_sample) {
    {
      Real2D stats_values_44(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_values_44"),
          rtpthlp.extent(0), rtpthlp.extent(1));
      Kokkos::parallel_for(
          "advance_xp2_xpyp::stats_value", rtpthlp.size(), KOKKOS_LAMBDA(int cell) {
            const int i           = cell % rtpthlp.extent(0);
            const int k           = cell / rtpthlp.extent(0);
            stats_values_44(i, k) = rtpthlp(i, k) / dt;
          });
      f2k::stats_begin_budget("rtpthlp_cl", stats_values_44, stats);
    }
  }
  f2k::clip_covar(nzm, ngrdcol, xp2_xpyp_rtpthlp, rtp2, thlp2, // In
                  rtpthlp,                                     // InOut
                  rtpthlp_chnge);                              // Out
  if (stats.flags().l_sample) {
    {
      Real2D stats_values_45(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp2_xpyp::stats_values_45"),
          rtpthlp.extent(0), rtpthlp.extent(1));
      Kokkos::parallel_for(
          "advance_xp2_xpyp::stats_value", rtpthlp.size(), KOKKOS_LAMBDA(int cell) {
            const int i           = cell % rtpthlp.extent(0);
            const int k           = cell / rtpthlp.extent(0);
            stats_values_45(i, k) = rtpthlp(i, k) / dt;
          });
      f2k::stats_finalize_budget("rtpthlp_cl", stats_values_45, stats);
    }
  }

  if (l_scalar_calc) {

    // Apply hole filling algorithm to the scalar variance terms
    if (fill_holes_type != 0) {
      for (sclr = 1; sclr <= sclr_dim; sclr += 1) {
        pos_definite_variances(
            nzm, ngrdcol, gr, xp2_xpyp_sclrp2, fill_holes_type, dt,                  // In
            (to_host_space(sclr_tol)(sclr - 1) * to_host_space(sclr_tol)(sclr - 1)), // In
            rho_ds_zm,                                                               // In
            stats, Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)));   // InOut

        if (sclr == sclr_idx.iisclr_rt) {
          // Here again, we do this kluge here to make sclr'rt' == rt'^2
          pos_definite_variances(
              nzm, ngrdcol, gr, xp2_xpyp_sclrprtp, fill_holes_type, dt,                // In
              (to_host_space(sclr_tol)(sclr - 1) * to_host_space(sclr_tol)(sclr - 1)), // In
              rho_ds_zm,                                                               // In
              stats, Kokkos::subview(sclrprtp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // InOut
        }
        if (sclr == sclr_idx.iisclr_thl) {
          // As with sclr'rt' above, but for sclr'thl'
          pos_definite_variances(
              nzm, ngrdcol, gr, xp2_xpyp_sclrpthlp, fill_holes_type, dt,                // In
              (to_host_space(sclr_tol)(sclr - 1) * to_host_space(sclr_tol)(sclr - 1)),  // In
              rho_ds_zm,                                                                // In
              stats, Kokkos::subview(sclrpthlp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // InOut
        }
      }
    }

    // Clipping for sclr'^2
    for (sclr = 1; sclr <= sclr_dim; sclr += 1) {

      Kokkos::parallel_for(
          "advance_xp2_xpyp::loop46", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            threshold_array(i - 1, k - 1) = (sclr_tol(sclr - 1) * sclr_tol(sclr - 1));
          });

      f2k::clip_variance(nzm, ngrdcol, gr, clip_sclrp2, dt, threshold_array, // In
                         stats,
                         Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // InOut
    }

    // Clipping for sclr'r_t'
    // Clipping sclr'r_t' at each vertical level, based on the
    // correlation of sclr and r_t at each vertical level, such that:
    // corr_(sclr,r_t) = sclr'r_t' / [ sqrt(sclr'^2) * sqrt(r_t'^2) ];
    // -1 <= corr_(sclr,r_t) <= 1.
    // Since sclr'^2, r_t'^2, and sclr'r_t' are all computed in the
    // same place, clipping for sclr'r_t' only has to be done once.
    for (sclr = 1; sclr <= sclr_dim; sclr += 1) {

      if (sclr == sclr_idx.iisclr_rt) {
        // Treat this like a variance if we're emulating rt
        Kokkos::parallel_for(
            "advance_xp2_xpyp::loop47", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              threshold_array(i - 1, k - 1) = sclr_tol(sclr - 1) * rt_tol;
            });

        f2k::clip_variance(
            nzm, ngrdcol, gr, clip_sclrprtp, dt, threshold_array,                    // In
            stats, Kokkos::subview(sclrprtp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // InOut
      } else {
        f2k::clip_covar(
            nzm, ngrdcol, clip_sclrprtp,                                            // In
            Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), rtp2,    // In
            Kokkos::subview(sclrprtp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),        // InOut
            Kokkos::subview(sclrprtp_chnge, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out
      }
    }

    // Clipping for sclr'th_l'
    // Clipping sclr'th_l' at each vertical level, based on the
    // correlation of sclr and th_l at each vertical level, such that:
    // corr_(sclr,th_l) = sclr'th_l' / [ sqrt(sclr'^2) * sqrt(th_l'^2) ];
    // -1 <= corr_(sclr,th_l) <= 1.
    // Since sclr'^2, th_l'^2, and sclr'th_l' are all computed in the
    // same place, clipping for sclr'th_l' only has to be done once.
    for (sclr = 1; sclr <= sclr_dim; sclr += 1) {
      if (sclr == sclr_idx.iisclr_thl) {
        // As above, but for thl
        Kokkos::parallel_for(
            "advance_xp2_xpyp::loop48", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              threshold_array(i - 1, k - 1) = sclr_tol(sclr - 1) * thl_tol;
            });

        f2k::clip_variance(
            nzm, ngrdcol, gr, clip_sclrpthlp, dt, threshold_array,                    // In
            stats, Kokkos::subview(sclrpthlp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // InOut
      } else {
        f2k::clip_covar(
            nzm, ngrdcol, clip_sclrpthlp,                                            // In
            Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), thlp2,    // In
            Kokkos::subview(sclrpthlp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),        // InOut
            Kokkos::subview(sclrpthlp_chnge, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out
      }
    }
  }

  if (f2k::clubb_at_least_debug_level_api(0)) {

    if (any_equal(err_info.err_code, clubb_fatal_error)) {

      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "Error in advance_xp2_xpyp" << "\n";
        std::cerr << message.str();
      }

      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "Intent(in)" << "\n";
        std::cerr << message.str();
      }

      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "invrs_tau_xp2_zm = " << " "
                << dump_values(Kokkos::subview(invrs_tau_xp2_zm, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "invrs_tau_C4_zm = " << " "
                << dump_values(Kokkos::subview(invrs_tau_C4_zm, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "invrs_tau_C14_zm = " << " "
                << dump_values(Kokkos::subview(invrs_tau_C14_zm, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "wm_zm = " << " "
                << dump_values(Kokkos::subview(wm_zm, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "rtm = " << " " << dump_values(Kokkos::subview(rtm, Kokkos::ALL, Kokkos::ALL))
                << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "wprtp = " << " "
                << dump_values(Kokkos::subview(wprtp, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "thlm = " << " " << dump_values(Kokkos::subview(thlm, Kokkos::ALL, Kokkos::ALL))
                << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "wpthlp = " << " "
                << dump_values(Kokkos::subview(wpthlp, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "wpthvp = " << " "
                << dump_values(Kokkos::subview(wpthvp, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "um = " << " " << dump_values(Kokkos::subview(um, Kokkos::ALL, Kokkos::ALL))
                << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "vm = " << " " << dump_values(Kokkos::subview(vm, Kokkos::ALL, Kokkos::ALL))
                << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "wp2 = " << " " << dump_values(Kokkos::subview(wp2, Kokkos::ALL, Kokkos::ALL))
                << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "wp3 = " << " " << dump_values(Kokkos::subview(wp3, Kokkos::ALL, Kokkos::ALL))
                << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "upwp = " << " " << dump_values(Kokkos::subview(upwp, Kokkos::ALL, Kokkos::ALL))
                << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "vpwp = " << " " << dump_values(Kokkos::subview(vpwp, Kokkos::ALL, Kokkos::ALL))
                << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "sigma_sqd_w = " << " "
                << dump_values(Kokkos::subview(sigma_sqd_w, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "Kh_zt = " << " "
                << dump_values(Kokkos::subview(Kh_zt, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "rtp2_forcing = " << " "
                << dump_values(Kokkos::subview(rtp2_forcing, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "thlp2_forcing = " << " "
                << dump_values(Kokkos::subview(thlp2_forcing, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "rtpthlp_forcing = " << " "
                << dump_values(Kokkos::subview(rtpthlp_forcing, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "rho_ds_zm = " << " "
                << dump_values(Kokkos::subview(rho_ds_zm, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "rho_ds_zt = " << " "
                << dump_values(Kokkos::subview(rho_ds_zt, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "invrs_rho_ds_zm = " << " "
                << dump_values(Kokkos::subview(invrs_rho_ds_zm, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "thv_ds_zm = " << " "
                << dump_values(Kokkos::subview(thv_ds_zm, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "wp2_zt = " << " "
                << dump_values(Kokkos::subview(wp2_zt, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }

      for (sclr = 1; sclr <= sclr_dim; sclr += 1) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "sclrm = " << " " << sclr << " "
                  << dump_values(Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)))
                  << "\n";
          std::cerr << message.str();
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "wpsclrp = " << " " << sclr << " "
                  << dump_values(Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)))
                  << "\n";
          std::cerr << message.str();
        }
      }

      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "Intent(In/Out)" << "\n";
        std::cerr << message.str();
      }

      if (l_lmm_stepping) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "rtp2 (pre-solve) = " << " "
                  << dump_values(Kokkos::subview(rtp2_old, Kokkos::ALL, Kokkos::ALL)) << "\n";
          std::cerr << message.str();
        }
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "rtp2 = " << " " << dump_values(Kokkos::subview(rtp2, Kokkos::ALL, Kokkos::ALL))
                << "\n";
        std::cerr << message.str();
      }
      if (l_lmm_stepping) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "thlp2 (pre-solve) = " << " "
                  << dump_values(Kokkos::subview(thlp2_old, Kokkos::ALL, Kokkos::ALL)) << "\n";
          std::cerr << message.str();
        }
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "thlp2 = " << " "
                << dump_values(Kokkos::subview(thlp2, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      if (l_lmm_stepping) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "rtpthlp (pre-solve) = " << " "
                  << dump_values(Kokkos::subview(rtpthlp_old, Kokkos::ALL, Kokkos::ALL)) << "\n";
          std::cerr << message.str();
        }
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "rtpthlp = " << " "
                << dump_values(Kokkos::subview(rtpthlp, Kokkos::ALL, Kokkos::ALL)) << "\n";
        std::cerr << message.str();
      }
      if (l_lmm_stepping) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "up2 (pre-solve) = " << " "
                  << dump_values(Kokkos::subview(up2_old, Kokkos::ALL, Kokkos::ALL)) << "\n";
          std::cerr << message.str();
        }
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "up2 = " << " " << dump_values(Kokkos::subview(up2, Kokkos::ALL, Kokkos::ALL))
                << "\n";
        std::cerr << message.str();
      }
      if (l_lmm_stepping) {
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "vp2 (pre-solve) = " << " "
                  << dump_values(Kokkos::subview(vp2_old, Kokkos::ALL, Kokkos::ALL)) << "\n";
          std::cerr << message.str();
        }
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "vp2 = " << " " << dump_values(Kokkos::subview(vp2, Kokkos::ALL, Kokkos::ALL))
                << "\n";
        std::cerr << message.str();
      }

      for (sclr = 1; sclr <= sclr_dim; sclr += 1) {
        if (l_lmm_stepping) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "sclrp2 (pre-solve) = " << " " << sclr << " "
                    << dump_values(
                           Kokkos::subview(sclrp2_old, Kokkos::ALL, Kokkos::ALL, (sclr - 1)))
                    << "\n";
            std::cerr << message.str();
          }
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "sclrp2 = " << " " << sclr << " "
                  << dump_values(Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)))
                  << "\n";
          std::cerr << message.str();
        }
        if (l_lmm_stepping) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "sclrprtp (pre-solve) = " << " " << sclr << " "
                    << dump_values(
                           Kokkos::subview(sclrprtp_old, Kokkos::ALL, Kokkos::ALL, (sclr - 1)))
                    << "\n";
            std::cerr << message.str();
          }
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "sclrprtp = " << " " << sclr << " "
                  << dump_values(Kokkos::subview(sclrprtp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)))
                  << "\n";
          std::cerr << message.str();
        }
        if (l_lmm_stepping) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "sclrthlp (pre-solve) = " << " " << sclr << " "
                    << dump_values(
                           Kokkos::subview(sclrpthlp_old, Kokkos::ALL, Kokkos::ALL, (sclr - 1)))
                    << "\n";
            std::cerr << message.str();
          }
        }
        {
          std::ostringstream message;
          message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          message << "sclrthlp = " << " " << sclr << " "
                  << dump_values(Kokkos::subview(sclrpthlp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)))
                  << "\n";
          std::cerr << message.str();
        }
      }
    }
  }

  return;
}

namespace {
// nzm: Number of momentum vertical levels
// nzt: Number of thermodynamic vertical levels
// ngrdcol: Number of grid columns
// sclr_dim: Number of passive scalars
// sclr_tol: Threshold(s) on the passive scalars  [units vary]
// invrs_tau_xp2_zm: Inverse time-scale for xp2 on momentum levels [1/s]
// rtm: Total water mixing ratio (t-levs)     [kg/kg]
// thlm: Liquid potential temp. (t-levs)       [K]
// wprtp: <w'r_t'> (momentum levels)            [(m/s)(kg/kg)]
// wpthlp: <w'th_l'> (momentum levels)           [(m K)/s]
// rtp2_forcing: <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
// thlp2_forcing: <th_l'^2> forcing (momentum levels)   [K^2/s]
// rtpthlp_forcing: <r_t'th_l'> forcing (momentum levels) [(kg/kg)K/s]
// sclrm: Mean value; pass. scalar (t-levs.) [{sclr units}]
// wpsclrp: <w'sclr'> (momentum levels)        [m/s{sclr units}]
// lhs_ma: Mean advection contributions to lhs
// lhs_diff: Diffusion contributions to lhs, dissipation term 2
// rhs_ta_wprtp2: For <w'rt'^2>
// rhs_ta_wpthlp2: For <w'thl'^2>
// rhs_ta_wprtpthlp: For <w'rt'thl'>
// rhs_ta_wpsclrp2: For <w'sclr'^2>
// rhs_ta_wprtpsclrp: For <w'sclr'rt'><w'sclr'^2><w'sclr'thl'>
// rhs_ta_wpthlpsclrp: For <w'sclr'thl'>
// dt: Model timestep                                [s]
// l_implemented: True if CLUBB is being implemented and run in a host model
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// rtp2: <r_t'^2>                      [(kg/kg)^2]
// thlp2: <th_l'^2>                     [K^2]
// rtpthlp: <r_t'th_l'>                   [(kg K)/kg]
// err_info: err_info struct containing err_code and err_header
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26>
void solve_xp2_xpyp_with_single_lhs(int nzm, int nzt, int ngrdcol, int sclr_dim, Array0 sclr_tol,
                                    const grid& gr,                // In
                                    const sclr_idx_type& sclr_idx, // In
                                    Array1 C2x, Array2 invrs_tau_xp2_zm, Array3 rtm, Array4 thlm,
                                    Array5 wprtp,                                             // In
                                    Array6 wpthlp, Array7 rtp2_forcing, Array8 thlp2_forcing, // In
                                    Array9 rtpthlp_forcing, Array10 sclrm, Array11 wpsclrp,   // In
                                    Array12 lhs_ta, Array13 lhs_ma, Array14 lhs_diff,         // In
                                    Array15 rhs_ta_wprtp2, Array16 rhs_ta_wpthlp2,            // In
                                    Array17 rhs_ta_wprtpthlp, Array18 rhs_ta_wpsclrp2,        // In
                                    Array19 rhs_ta_wprtpsclrp, Array20 rhs_ta_wpthlpsclrp,    // In
                                    real_t dt, bool l_scalar_calc, bool l_lmm_stepping,       // In
                                    bool l_implemented,                                       // In
                                    int tridiag_solve_method,                                 // In
                                    f2k::stats_type& stats,                       // InOut
                                    Array21 rtp2, Array22 thlp2, Array23 rtpthlp, // InOut
                                    Array24 sclrp2, Array25 sclrprtp, Array26 sclrpthlp,
                                    err_info_type& err_info) // InOut
{
  // Description:
  //     This subroutine generates a single lhs matrix and multiple rhs matricies, then
  //     solves the system. This should only be used when the lhs matrices for
  //     <w'r_t'>, <w'th_l'>, <w'rt'thl'> ,<w'sclr'^2>, <w'sclr'rt'>,  and <w'sclr'thl'>
  //     are identical. Otherwise multiple lhs matrices and multiple solves are required.
  //
  //     The ADG1 PDF must be in use for this subroutine to have the
  //     potential to be called.
  //----------------------------------------------------------------------------------
  // -------- Input Variables --------
  // Passive scalar input
  // -------- In/Out Variables --------
  // An attribute of (inout) is also needed to import the value of the variances
  // at the surface.  Brian.  12/18/05.
  // -------- Local Variables --------
  // LHS dissipation term 1 (dp1). An "over-implicit" weighted time
  // step is applied to this term (and to pressure term 1 for u'^2 and v'^2).
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:xp2_dp
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "solve_xp2_xpyp_with_single_lhs::lhs"),
      ndiags3, ngrdcol, nzm); // Tridiagonal matrix
  Kokkos::View<real_t***, Kokkos::LayoutLeft> rhs(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "solve_xp2_xpyp_with_single_lhs::rhs"),
      ngrdcol, nzm, 3 + 3 * sclr_dim);
  Kokkos::View<real_t***, Kokkos::LayoutLeft> solution(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "solve_xp2_xpyp_with_single_lhs::solution"),
      ngrdcol, nzm, 3 + 3 * sclr_dim);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sclrp2_forcing(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_single_lhs::sclrp2_forcing"),
      ngrdcol, nzm); // <sclr'^2> forcing (momentum levels)    [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sclrprtp_forcing(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_single_lhs::sclrprtp_forcing"),
      ngrdcol, nzm); // <sclr'r_t'> forcing (momentum levels)  [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sclrpthlp_forcing(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_single_lhs::sclrpthlp_forcing"),
      ngrdcol, nzm); // <sclr'th_l'> forcing (momentum levels) [units vary]
  real_t threshold;  // Minimum value for variances                   [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> lhs_dp1(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "solve_xp2_xpyp_with_single_lhs::lhs_dp1"),
      ngrdcol, nzm); // LHS dissipation term 1
  int sclr;
  int k;
  int i;

  // -------- Begin Code --------

  term_dp1_lhs(nzm, ngrdcol, gr, C2x, invrs_tau_xp2_zm, // In
               lhs_dp1);                                // Out

  Kokkos::parallel_for(
      "solve_xp2_xpyp_with_single_lhs::loop1", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        lhs_dp1(i - 1, k - 1) = lhs_dp1(i - 1, k - 1) * gamma_over_implicit_ts;
      });

  // Calculate lhs matrix
  xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta, lhs_ma, lhs_diff, lhs_dp1, // In
               lhs);                                                    // Out

  // Calculate rhs matricies
  xp2_xpyp_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_rtp2, dt, wprtp, wprtp, rtm, rtm, rtp2,        // In
               rtp2_forcing, C2x, invrs_tau_xp2_zm, (rt_tol * rt_tol), lhs_ta, rhs_ta_wprtp2, // In
               stats,                                                    // InOut
               Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, (1 - 1))); // Out

  xp2_xpyp_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_thlp2, dt, wpthlp, wpthlp, thlm, thlm, thlp2, // In
               thlp2_forcing, C2x, invrs_tau_xp2_zm, (thl_tol * thl_tol), lhs_ta,            // In
               rhs_ta_wpthlp2,                                                               // In
               stats,                                                    // InOut
               Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, (2 - 1))); // Out

  xp2_xpyp_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_rtpthlp, dt, wprtp, wpthlp, rtm, thlm, rtpthlp, // In
               rtpthlp_forcing, C2x, invrs_tau_xp2_zm, zero_threshold, lhs_ta,
               rhs_ta_wprtpthlp,                                         // In
               stats,                                                    // InOut
               Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, (3 - 1))); // Out

  if (l_scalar_calc) {

    for (sclr = 1; sclr <= sclr_dim; sclr += 1) {

      // Forcing for <sclr'^2>.
      Kokkos::parallel_for(
          "solve_xp2_xpyp_with_single_lhs::loop2", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            sclrp2_forcing(i - 1, k - 1) = zero;
          });

      //!!!!***** sclr'^2 *****!!!!!
      xp2_xpyp_rhs(
          nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrp2, dt,                                         // In
          Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),                     // In
          Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),                     // In
          Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),                       // In
          Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),                       // In
          Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), sclrp2_forcing, C2x, // In
          invrs_tau_xp2_zm,                                                                   // In
          (to_host_space(sclr_tol)(sclr - 1) * to_host_space(sclr_tol)(sclr - 1)), lhs_ta,    // In
          Kokkos::subview(rhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),             // In
          stats,                                                           // InOut
          Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, (3 + sclr - 1))); // Out

      //!!!!***** sclr'r_t' *****!!!!!
      if (sclr == sclr_idx.iisclr_rt) {
        // In this case we're trying to emulate rt'^2 with sclr'rt', so we
        // handle this as we would a variance, even though generally speaking
        // the scalar is not rt
        Kokkos::parallel_for(
            "solve_xp2_xpyp_with_single_lhs::loop3", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              sclrprtp_forcing(i - 1, k - 1) = rtp2_forcing(i - 1, k - 1);
            });

        threshold = (rt_tol * rt_tol);
      } else {
        Kokkos::parallel_for(
            "solve_xp2_xpyp_with_single_lhs::loop4", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              sclrprtp_forcing(i - 1, k - 1) = zero;
            });

        threshold = zero_threshold;
      }

      xp2_xpyp_rhs(
          nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrprtp, dt,                                      // In
          Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), wprtp,             // In
          Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), rtm,                 // In
          Kokkos::subview(sclrprtp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), sclrprtp_forcing, // In
          C2x, invrs_tau_xp2_zm, threshold, lhs_ta,                                          // In
          Kokkos::subview(rhs_ta_wprtpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),          // In
          stats,                                                                      // InOut
          Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, (3 + sclr + sclr_dim - 1))); // Out

      //!!!!***** sclr'th_l' *****!!!!!
      if (sclr == sclr_idx.iisclr_thl) {
        // In this case we're trying to emulate thl'^2 with sclr'thl', so we
        // handle this as we did with sclr_rt, above.
        Kokkos::parallel_for(
            "solve_xp2_xpyp_with_single_lhs::loop5", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              sclrpthlp_forcing(i - 1, k - 1) = thlp2_forcing(i - 1, k - 1);
            });

        threshold = (thl_tol * thl_tol);
      } else {
        Kokkos::parallel_for(
            "solve_xp2_xpyp_with_single_lhs::loop6", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              sclrpthlp_forcing(i - 1, k - 1) = zero;
            });

        threshold = zero_threshold;
      }

      xp2_xpyp_rhs(
          nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrpthlp, dt,                                       // In
          Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), wpthlp,              // In
          Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), thlm,                  // In
          Kokkos::subview(sclrpthlp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), sclrpthlp_forcing, // In
          C2x, invrs_tau_xp2_zm, threshold, lhs_ta,                                            // In
          Kokkos::subview(rhs_ta_wpthlpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),           // In
          stats,                                                                          // InOut
          Kokkos::subview(rhs, Kokkos::ALL, Kokkos::ALL, (3 + sclr + 2 * sclr_dim - 1))); // Out
    }
  }

  // Solve multiple rhs with single lhs
  xp2_xpyp_solve(nzm, ngrdcol, xp2_xpyp_single_lhs, 3 + 3 * sclr_dim, l_implemented, gr, // In
                 tridiag_solve_method,                                                   // In
                 stats, rhs, lhs, err_info,                                              // InOut
                 solution);                                                              // Out

  if (f2k::clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "calling xp2_xpyp_solve in solve_xp2_xpyp_with_single_lhs" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }

  if (l_lmm_stepping) {

    Kokkos::parallel_for(
        "solve_xp2_xpyp_with_single_lhs::loop7", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          rtp2(i - 1, k - 1)  = one_half * (rtp2(i - 1, k - 1) + solution(i - 1, k - 1, 1 - 1));
          thlp2(i - 1, k - 1) = one_half * (thlp2(i - 1, k - 1) + solution(i - 1, k - 1, 2 - 1));
          rtpthlp(i - 1, k - 1) =
              one_half * (rtpthlp(i - 1, k - 1) + solution(i - 1, k - 1, 3 - 1));
        });

    if (sclr_dim > 0) {
      Kokkos::parallel_for(
          "solve_xp2_xpyp_with_single_lhs::loop8",
          std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            const int sclr = 1 + remaining % (sclr_dim);
            remaining /= (sclr_dim);
            sclrp2(i - 1, k - 1, sclr - 1) =
                one_half * (sclrp2(i - 1, k - 1, sclr - 1) + solution(i - 1, k - 1, 3 + sclr - 1));
            sclrprtp(i - 1, k - 1, sclr - 1) =
                one_half * (sclrprtp(i - 1, k - 1, sclr - 1) +
                            solution(i - 1, k - 1, 3 + sclr_dim + sclr - 1));
            sclrpthlp(i - 1, k - 1, sclr - 1) =
                one_half * (sclrpthlp(i - 1, k - 1, sclr - 1) +
                            solution(i - 1, k - 1, 3 + 2 * sclr_dim + sclr - 1));
          });
    }

  } else {

    Kokkos::parallel_for(
        "solve_xp2_xpyp_with_single_lhs::loop9", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          rtp2(i - 1, k - 1)    = solution(i - 1, k - 1, 1 - 1);
          thlp2(i - 1, k - 1)   = solution(i - 1, k - 1, 2 - 1);
          rtpthlp(i - 1, k - 1) = solution(i - 1, k - 1, 3 - 1);
        });

    if (sclr_dim > 0) {
      Kokkos::parallel_for(
          "solve_xp2_xpyp_with_single_lhs::loop10",
          std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            const int sclr = 1 + remaining % (sclr_dim);
            remaining /= (sclr_dim);
            sclrp2(i - 1, k - 1, sclr - 1)    = solution(i - 1, k - 1, 3 + sclr - 1);
            sclrprtp(i - 1, k - 1, sclr - 1)  = solution(i - 1, k - 1, 3 + sclr_dim + sclr - 1);
            sclrpthlp(i - 1, k - 1, sclr - 1) = solution(i - 1, k - 1, 3 + 2 * sclr_dim + sclr - 1);
          });
    }
  }

  if (stats.flags().l_sample) {
    {

      stats_finalize_xp2_xpyp_terms(nzm, ngrdcol, xp2_xpyp_rtp2, gr, rtp2, lhs_dp1, lhs_dp1,
                                    lhs_diff, lhs_ta, // In
                                    lhs_ma,           // In
                                    stats);           // InOut
    }
    {

      stats_finalize_xp2_xpyp_terms(nzm, ngrdcol, xp2_xpyp_thlp2, gr, thlp2, lhs_dp1, lhs_dp1,
                                    lhs_diff, lhs_ta, // In
                                    lhs_ma,           // In
                                    stats);           // InOut
    }
    {

      stats_finalize_xp2_xpyp_terms(nzm, ngrdcol, xp2_xpyp_rtpthlp, gr, rtpthlp, lhs_dp1, lhs_dp1,
                                    lhs_diff, lhs_ta, // In
                                    lhs_ma,           // In
                                    stats);           // InOut
    }
  }

  return;
}

// nzm: Number of momentum vertical levels
// nzt: Number of thermodynamic vertical levels
// ngrdcol: Number of grid columns
// sclr_dim: Number of passive scalars
// sclr_tol: Threshold(s) on the passive scalars  [units vary]
// invrs_tau_xp2_zm: Inverse time-scale for xp2 on momentum levels [1/s]
// rtm: Total water mixing ratio (t-levs)     [kg/kg]
// thlm: Liquid potential temp. (t-levs)       [K]
// wprtp: <w'r_t'> (momentum levels)            [(m/s)(kg/kg)]
// wpthlp: <w'th_l'> (momentum levels)           [(m K)/s]
// rtp2_forcing: <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
// thlp2_forcing: <th_l'^2> forcing (momentum levels)   [K^2/s]
// rtpthlp_forcing: <r_t'th_l'> forcing (momentum levels) [(kg/kg)K/s]
// sclrm: Mean value; pass. scalar (t-levs.) [{sclr units}]
// wpsclrp: <w'sclr'> (momentum levels)        [m/s{sclr units}]
// lhs_ta_wprtp2: Turbulent advection term for <w'rt'^2>
// lhs_ta_wpthlp2: Turbulent advection term for <w'thl'^2>
// lhs_ta_wprtpthlp: Turbulent advection term for <w'rtp'thl'>
// lhs_ta_wpsclrp2: For <w'sclr'^2>
// lhs_ta_wprtpsclrp: For <w'rt'sclr'>
// lhs_ta_wpthlpsclrp: For <w'thl'sclr'>
// lhs_ma: Mean advection contributions to lhs
// lhs_diff: Diffusion contributions to lhs, dissipation term 2
// rhs_ta_wprtp2: For <w'rt'^2>
// rhs_ta_wpthlp2: For <w'thl'^2>
// rhs_ta_wprtpthlp: For <w'rt'thl'>
// rhs_ta_wpsclrp2: For <w'sclr'^2>
// rhs_ta_wprtpsclrp: For <w'sclr'rt'>
// rhs_ta_wpthlpsclrp: For <w'sclr'thl'>
// dt: Model timestep                                [s]
// iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use for
// the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component
// PDF.
// l_implemented: True if CLUBB is being implemented and run in a host model
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// rtp2: <r_t'^2>                      [(kg/kg)^2]
// thlp2: <th_l'^2>                     [K^2]
// rtpthlp: <r_t'th_l'>                   [(kg K)/kg]
// err_info: err_info struct containing err_code and err_header
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34>
void solve_xp2_xpyp_with_multiple_lhs(
    int nzm, int nzt, int ngrdcol, int sclr_dim, Array0 sclr_tol, const grid& gr,  // In
    const sclr_idx_type& sclr_idx,                                                 // In
    Array1 C2rt_1d, Array2 C2thl_1d, Array3 C2rtthl_1d, Array4 C2sclr_1d,          // In
    Array5 invrs_tau_xp2_zm, Array6 rtm, Array7 thlm, Array8 wprtp, Array9 wpthlp, // In
    Array10 rtp2_forcing, Array11 thlp2_forcing, Array12 rtpthlp_forcing,          // In
    Array13 sclrm, Array14 wpsclrp,                                                // In
    Array15 lhs_ta_wprtp2, Array16 lhs_ta_wpthlp2,                                 // In
    Array17 lhs_ta_wprtpthlp, Array18 lhs_ta_wpsclrp2,                             // In
    Array19 lhs_ta_wprtpsclrp, Array20 lhs_ta_wpthlpsclrp,                         // In
    Array21 lhs_ma, Array22 lhs_diff,                                              // In
    Array23 rhs_ta_wprtp2, Array24 rhs_ta_wpthlp2,                                 // In
    Array25 rhs_ta_wprtpthlp, Array26 rhs_ta_wpsclrp2,                             // In
    Array27 rhs_ta_wprtpsclrp, Array28 rhs_ta_wpthlpsclrp,                         // In
    real_t dt, int iiPDF_type, bool l_scalar_calc,                                 // In
    bool l_lmm_stepping,                                                           // In
    bool l_implemented,                                                            // In
    int tridiag_solve_method,                                                      // In
    f2k::stats_type& stats,                                                        // InOut
    Array29 rtp2, Array30 thlp2, Array31 rtpthlp,                                  // InOut
    Array32 sclrp2, Array33 sclrprtp, Array34 sclrpthlp, err_info_type& err_info)  // InOut
{
  // Description:
  //     This subroutine generates different lhs and rhs matrices to solve for.
  //
  //-----------------------------------------------------------------------
  //------------------------ Input Variables ------------------------
  // Passive scalar input
  //------------------------ In/Out Variables ------------------------
  // Input/Output variables
  // An attribute of (inout) is also needed to import the value of the variances
  // at the surface.  Brian.  12/18/05.
  //------------------------ Local Variables ------------------------
  Kokkos::View<real_t***, Kokkos::LayoutLeft> lhs(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "solve_xp2_xpyp_with_multiple_lhs::lhs"),
      ndiags3, ngrdcol, nzm); // Tridiagonal matrix
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "solve_xp2_xpyp_with_multiple_lhs::rhs"),
      ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sclrp2_forcing(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_multiple_lhs::sclrp2_forcing"),
      ngrdcol, nzm); // <sclr'^2> forcing (momentum levels)    [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sclrprtp_forcing(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_multiple_lhs::sclrprtp_forcing"),
      ngrdcol, nzm); // <sclr'r_t'> forcing (momentum levels)  [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sclrpthlp_forcing(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_multiple_lhs::sclrpthlp_forcing"),
      ngrdcol, nzm); // <sclr'th_l'> forcing (momentum levels) [units vary]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sclr_rhs(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "solve_xp2_xpyp_with_multiple_lhs::sclr_rhs"),
      ngrdcol, nzm, 3 * sclr_dim); // RHS vectors of tridiagonal system for the passive scalars
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sclr_solution(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_multiple_lhs::sclr_solution"),
      ngrdcol, nzm, 3 * sclr_dim); // Solution to tridiagonal system for the passive scalars
  real_t threshold;                // Minimum value for variances                   [units vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> lhs_dp1(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "solve_xp2_xpyp_with_multiple_lhs::lhs_dp1"),
      ngrdcol, nzm); // LHS dissipation term 1
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rtp2_solution(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_multiple_lhs::rtp2_solution"),
      ngrdcol, nzm); // <r_t'^2>                      [(kg/kg)^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> thlp2_solution(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_multiple_lhs::thlp2_solution"),
      ngrdcol, nzm); // <th_l'^2>                     [K^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rtpthlp_solution(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_multiple_lhs::rtpthlp_solution"),
      ngrdcol, nzm); // <r_t'th_l'>                   [(kg K)/kg]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sclrp2_solution(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_multiple_lhs::sclrp2_solution"),
      ngrdcol, nzm, sclr_dim);
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sclrprtp_solution(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_multiple_lhs::sclrprtp_solution"),
      ngrdcol, nzm, sclr_dim);
  Kokkos::View<real_t***, Kokkos::LayoutLeft> sclrpthlp_solution(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "solve_xp2_xpyp_with_multiple_lhs::sclrpthlp_solution"),
      ngrdcol, nzm, sclr_dim);
  int sclr;
  int k;
  int i;

  //------------------------ Begin Code ------------------------

  //!!!!***** r_t'^2 *****!!!!!
  term_dp1_lhs(nzm, ngrdcol, gr, C2rt_1d, invrs_tau_xp2_zm, // In
               lhs_dp1);                                    // Out

  Kokkos::parallel_for(
      "solve_xp2_xpyp_with_multiple_lhs::loop1", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        lhs_dp1(i - 1, k - 1) = lhs_dp1(i - 1, k - 1) * gamma_over_implicit_ts;
      });

  // Implicit contributions to term rtp2
  xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wprtp2, lhs_ma, lhs_diff, lhs_dp1, // In
               lhs);                                                           // Out

  xp2_xpyp_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_rtp2, dt, wprtp, wprtp, rtm, rtm, rtp2,    // In
               rtp2_forcing, C2rt_1d, invrs_tau_xp2_zm, (rt_tol * rt_tol), lhs_ta_wprtp2, // In
               rhs_ta_wprtp2,                                                             // In
               stats,                                                                     // InOut
               rhs);                                                                      // Out

  // Solve the tridiagonal system
  xp2_xpyp_solve(nzm, ngrdcol, xp2_xpyp_rtp2, 1, l_implemented, gr, tridiag_solve_method, // In
                 stats, rhs, lhs, err_info,                                               // InOut
                 rtp2_solution);                                                          // Out

  if (f2k::clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "calling xp2_xpyp_solve for rtp2 in solve_xp2_xpyp_with_multiple_lhs" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }

  if (l_lmm_stepping) {
    Kokkos::parallel_for(
        "solve_xp2_xpyp_with_multiple_lhs::loop2", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          rtp2(i - 1, k - 1) = one_half * (rtp2(i - 1, k - 1) + rtp2_solution(i - 1, k - 1));
        });
  } else {
    Kokkos::parallel_for(
        "solve_xp2_xpyp_with_multiple_lhs::loop3", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          rtp2(i - 1, k - 1) = rtp2_solution(i - 1, k - 1);
        });
  }

  if (stats.flags().l_sample) {
    {

      stats_finalize_xp2_xpyp_terms(nzm, ngrdcol, xp2_xpyp_rtp2, gr, rtp2, lhs_dp1, lhs_dp1,
                                    lhs_diff, lhs_ta_wprtp2, // In
                                    lhs_ma,                  // In
                                    stats);                  // InOut
    }
  }

  //!!!!***** th_l'^2 *****!!!!!
  term_dp1_lhs(nzm, ngrdcol, gr, C2thl_1d, invrs_tau_xp2_zm, // In
               lhs_dp1);                                     // Out

  Kokkos::parallel_for(
      "solve_xp2_xpyp_with_multiple_lhs::loop4", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        lhs_dp1(i - 1, k - 1) = lhs_dp1(i - 1, k - 1) * gamma_over_implicit_ts;
      });

  // Implicit contributions to term thlp2
  xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wpthlp2, lhs_ma, lhs_diff, lhs_dp1, // In
               lhs);                                                            // Out

  // Explicit contributions to thlp2
  xp2_xpyp_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_thlp2, dt, wpthlp, wpthlp, thlm, thlm, thlp2,   // In
               thlp2_forcing, C2thl_1d, invrs_tau_xp2_zm, (thl_tol * thl_tol), lhs_ta_wpthlp2, // In
               rhs_ta_wpthlp2,                                                                 // In
               stats, // InOut
               rhs);  // Out

  // Solve the tridiagonal system
  xp2_xpyp_solve(nzm, ngrdcol, xp2_xpyp_thlp2, 1, l_implemented, gr, tridiag_solve_method, // In
                 stats, rhs, lhs, err_info,                                                // InOut
                 thlp2_solution);                                                          // Out

  if (f2k::clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "calling xp2_xpyp_solve for thlp2 in solve_xp2_xpyp_with_multiple_lhs" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }

  if (l_lmm_stepping) {
    Kokkos::parallel_for(
        "solve_xp2_xpyp_with_multiple_lhs::loop5", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          thlp2(i - 1, k - 1) = one_half * (thlp2(i - 1, k - 1) + thlp2_solution(i - 1, k - 1));
        });
  } else {
    Kokkos::parallel_for(
        "solve_xp2_xpyp_with_multiple_lhs::loop6", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          thlp2(i - 1, k - 1) = thlp2_solution(i - 1, k - 1);
        });
  }

  if (stats.flags().l_sample) {
    {

      stats_finalize_xp2_xpyp_terms(nzm, ngrdcol, xp2_xpyp_thlp2, gr, thlp2, lhs_dp1, lhs_dp1,
                                    lhs_diff,               // In
                                    lhs_ta_wpthlp2, lhs_ma, // In
                                    stats);                 // InOut
    }
  }

  //!!!!***** r_t'th_l' *****!!!!!
  term_dp1_lhs(nzm, ngrdcol, gr, C2rtthl_1d, invrs_tau_xp2_zm, // In
               lhs_dp1);                                       // Out

  Kokkos::parallel_for(
      "solve_xp2_xpyp_with_multiple_lhs::loop7", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        lhs_dp1(i - 1, k - 1) = lhs_dp1(i - 1, k - 1) * gamma_over_implicit_ts;
      });

  // Implicit contributions to term rtpthlp
  xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, lhs_ta_wprtpthlp, lhs_ma, lhs_diff, lhs_dp1, // In
               lhs);                                                              // Out

  // Explicit contributions to rtpthlp
  xp2_xpyp_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_rtpthlp, dt, wprtp, wpthlp, rtm, thlm, rtpthlp, // In
               rtpthlp_forcing, C2rtthl_1d, invrs_tau_xp2_zm, zero_threshold,
               lhs_ta_wprtpthlp, // In
               rhs_ta_wprtpthlp, // In
               stats,            // InOut
               rhs);             // Out

  // Solve the tridiagonal system
  xp2_xpyp_solve(nzm, ngrdcol, xp2_xpyp_rtpthlp, 1, l_implemented, gr, tridiag_solve_method, // In
                 stats, rhs, lhs, err_info, // InOut
                 rtpthlp_solution);         // Out

  if (f2k::clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "calling xp2_xpyp_solve for rtpthlp in solve_xp2_xpyp_with_multiple_lhs" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }

  if (l_lmm_stepping) {
    Kokkos::parallel_for(
        "solve_xp2_xpyp_with_multiple_lhs::loop8", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          rtpthlp(i - 1, k - 1) =
              one_half * (rtpthlp(i - 1, k - 1) + rtpthlp_solution(i - 1, k - 1));
        });
  } else {
    Kokkos::parallel_for(
        "solve_xp2_xpyp_with_multiple_lhs::loop9", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          rtpthlp(i - 1, k - 1) = rtpthlp_solution(i - 1, k - 1);
        });
  }

  if (stats.flags().l_sample) {
    {

      stats_finalize_xp2_xpyp_terms(nzm, ngrdcol, xp2_xpyp_rtpthlp, gr, rtpthlp, lhs_dp1, lhs_dp1,
                                    lhs_diff,                 // In
                                    lhs_ta_wprtpthlp, lhs_ma, // In
                                    stats);                   // InOut
    }
  }

  if (l_scalar_calc) {

    term_dp1_lhs(nzm, ngrdcol, gr, C2sclr_1d, invrs_tau_xp2_zm, // In
                 lhs_dp1);                                      // Out

    Kokkos::parallel_for(
        "solve_xp2_xpyp_with_multiple_lhs::loop10", std::max(0, nzm) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          lhs_dp1(i - 1, k - 1) = lhs_dp1(i - 1, k - 1) * gamma_over_implicit_ts;
        });

    if (iiPDF_type != iiPDF_ADG1) {

      // Any PDF besides ADG1 is used

      for (sclr = 1; sclr <= sclr_dim; sclr += 1) {

        // Forcing for <sclr'^2>.
        Kokkos::parallel_for(
            "solve_xp2_xpyp_with_multiple_lhs::array11",
            ((sclrp2_forcing.extent(0) - (1)) / (1) + 1) *
                ((sclrp2_forcing.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((sclrp2_forcing.extent(0) - (1)) / (1) + 1);
              remaining /= ((sclrp2_forcing.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((sclrp2_forcing.extent(1) - (1)) / (1) + 1);
              remaining /= ((sclrp2_forcing.extent(1) - (1)) / (1) + 1);
              sclrp2_forcing((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = zero;
            });

        //!!!!***** sclr'^2 *****!!!!!
        xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, // In
                     Kokkos::subview(lhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     (sclr - 1)), // In
                     lhs_ma, lhs_diff, lhs_dp1,   // In
                     lhs);                        // Out

        xp2_xpyp_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrp2, dt,                              // In
                     Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),          // In
                     Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),          // In
                     Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),            // In
                     Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),            // In
                     Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),           // In
                     Kokkos::subview(sclrp2_forcing, Kokkos::ALL, Kokkos::ALL), C2sclr_1d,    // In
                     Kokkos::subview(invrs_tau_xp2_zm, Kokkos::ALL, Kokkos::ALL),             // In
                     (to_host_space(sclr_tol)(sclr - 1) * to_host_space(sclr_tol)(sclr - 1)), // In
                     Kokkos::subview(lhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     (sclr - 1)),                                            // In
                     Kokkos::subview(rhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), // In
                     stats, // InOut
                     rhs);  // Out

        // Solve the tridiagonal system
        xp2_xpyp_solve(
            nzm, ngrdcol, xp2_xpyp_scalars, 1, l_implemented, gr, tridiag_solve_method, // In
            stats, rhs, lhs, err_info,                                                  // InOut
            Kokkos::subview(sclrp2_solution, Kokkos::ALL, Kokkos::ALL, (sclr - 1)));    // Out

        if (f2k::clubb_at_least_debug_level_api(0)) {
          if (any_equal(err_info.err_code, clubb_fatal_error)) {
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << err_info.err_header_global << "\n";
              std::cerr << message.str();
            }
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << "calling xp2_xpyp_solve for sclrp2 with sclr index " << " " << sclr << " "
                      << " in solve_xp2_xpyp_with_multiple_lhs" << "\n";
              std::cerr << message.str();
            }
            return;
          }
        }

        if (l_lmm_stepping) {
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::array12",
              ((sclrp2.extent(0) - (1)) / (1) + 1) * ((sclrp2.extent(1) - (1)) / (1) + 1),
              KOKKOS_LAMBDA(int cell) {
                int remaining     = cell;
                const int slice_0 = remaining % ((sclrp2.extent(0) - (1)) / (1) + 1);
                remaining /= ((sclrp2.extent(0) - (1)) / (1) + 1);
                const int slice_1 = remaining % ((sclrp2.extent(1) - (1)) / (1) + 1);
                remaining /= ((sclrp2.extent(1) - (1)) / (1) + 1);
                sclrp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1)) =
                    one_half *
                    (sclrp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1)) +
                     sclrp2_solution((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1)));
              });
        } else {
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::array13",
              ((sclrp2.extent(0) - (1)) / (1) + 1) * ((sclrp2.extent(1) - (1)) / (1) + 1),
              KOKKOS_LAMBDA(int cell) {
                int remaining     = cell;
                const int slice_0 = remaining % ((sclrp2.extent(0) - (1)) / (1) + 1);
                remaining /= ((sclrp2.extent(0) - (1)) / (1) + 1);
                const int slice_1 = remaining % ((sclrp2.extent(1) - (1)) / (1) + 1);
                remaining /= ((sclrp2.extent(1) - (1)) / (1) + 1);
                sclrp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1)) =
                    sclrp2_solution((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1));
              });
        }

        //!!!!***** sclr'r_t' *****!!!!!
        if (sclr == sclr_idx.iisclr_rt) {
          // In this case we're trying to emulate rt'^2 with sclr'rt', so we
          // handle this as we would a variance, even though generally speaking
          // the scalar is not rt
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::array14",
              ((sclrprtp_forcing.extent(0) - (1)) / (1) + 1) *
                  ((sclrprtp_forcing.extent(1) - (1)) / (1) + 1),
              KOKKOS_LAMBDA(int cell) {
                int remaining     = cell;
                const int slice_0 = remaining % ((sclrprtp_forcing.extent(0) - (1)) / (1) + 1);
                remaining /= ((sclrprtp_forcing.extent(0) - (1)) / (1) + 1);
                const int slice_1 = remaining % ((sclrprtp_forcing.extent(1) - (1)) / (1) + 1);
                remaining /= ((sclrprtp_forcing.extent(1) - (1)) / (1) + 1);
                sclrprtp_forcing((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                    rtp2_forcing((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)));
              });
          threshold = (rt_tol * rt_tol);
        } else {
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::array15",
              ((sclrprtp_forcing.extent(0) - (1)) / (1) + 1) *
                  ((sclrprtp_forcing.extent(1) - (1)) / (1) + 1),
              KOKKOS_LAMBDA(int cell) {
                int remaining     = cell;
                const int slice_0 = remaining % ((sclrprtp_forcing.extent(0) - (1)) / (1) + 1);
                remaining /= ((sclrprtp_forcing.extent(0) - (1)) / (1) + 1);
                const int slice_1 = remaining % ((sclrprtp_forcing.extent(1) - (1)) / (1) + 1);
                remaining /= ((sclrprtp_forcing.extent(1) - (1)) / (1) + 1);
                sclrprtp_forcing((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = zero;
              });
          threshold = zero_threshold;
        }

        xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, // In
                     Kokkos::subview(lhs_ta_wprtpsclrp, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     (sclr - 1)), // In
                     lhs_ma, lhs_diff, lhs_dp1,   // In
                     lhs);                        // Out

        xp2_xpyp_rhs(nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrprtp, dt,                          // In
                     Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), wprtp, // In
                     Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), rtm,     // In
                     Kokkos::subview(sclrprtp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),
                     sclrprtp_forcing,                       // In
                     C2sclr_1d, invrs_tau_xp2_zm, threshold, // In
                     Kokkos::subview(lhs_ta_wprtpsclrp, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     (sclr - 1)),                                              // In
                     Kokkos::subview(rhs_ta_wprtpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), // In
                     stats, // InOut
                     rhs);  // Out

        // Solve the tridiagonal system
        xp2_xpyp_solve(
            nzm, ngrdcol, xp2_xpyp_scalars, 1, l_implemented, gr, tridiag_solve_method, // In
            stats, rhs, lhs, err_info,                                                  // InOut
            Kokkos::subview(sclrprtp_solution, Kokkos::ALL, Kokkos::ALL, (sclr - 1)));  // Out

        if (f2k::clubb_at_least_debug_level_api(0)) {
          if (any_equal(err_info.err_code, clubb_fatal_error)) {
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << err_info.err_header_global << "\n";
              std::cerr << message.str();
            }
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << "calling xp2_xpyp_solve for sclrprtp with sclr index " << " " << sclr
                      << " " << " in solve_xp2_xpyp_with_multiple_lhs" << "\n";
              std::cerr << message.str();
            }
            return;
          }
        }

        if (l_lmm_stepping) {
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::array16",
              ((sclrprtp.extent(0) - (1)) / (1) + 1) * ((sclrprtp.extent(1) - (1)) / (1) + 1),
              KOKKOS_LAMBDA(int cell) {
                int remaining     = cell;
                const int slice_0 = remaining % ((sclrprtp.extent(0) - (1)) / (1) + 1);
                remaining /= ((sclrprtp.extent(0) - (1)) / (1) + 1);
                const int slice_1 = remaining % ((sclrprtp.extent(1) - (1)) / (1) + 1);
                remaining /= ((sclrprtp.extent(1) - (1)) / (1) + 1);
                sclrprtp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1)) =
                    one_half *
                    (sclrprtp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1)) +
                     sclrprtp_solution((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)),
                                       (sclr - 1)));
              });
        } else {
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::array17",
              ((sclrprtp.extent(0) - (1)) / (1) + 1) * ((sclrprtp.extent(1) - (1)) / (1) + 1),
              KOKKOS_LAMBDA(int cell) {
                int remaining     = cell;
                const int slice_0 = remaining % ((sclrprtp.extent(0) - (1)) / (1) + 1);
                remaining /= ((sclrprtp.extent(0) - (1)) / (1) + 1);
                const int slice_1 = remaining % ((sclrprtp.extent(1) - (1)) / (1) + 1);
                remaining /= ((sclrprtp.extent(1) - (1)) / (1) + 1);
                sclrprtp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1)) =
                    sclrprtp_solution((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1));
              });
        }

        //!!!!***** sclr'th_l' *****!!!!!

        if (sclr == sclr_idx.iisclr_thl) {
          // In this case we're trying to emulate thl'^2 with sclr'thl', so we
          // handle this as we did with sclr_rt, above.
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::array18",
              ((sclrpthlp_forcing.extent(0) - (1)) / (1) + 1) *
                  ((sclrpthlp_forcing.extent(1) - (1)) / (1) + 1),
              KOKKOS_LAMBDA(int cell) {
                int remaining     = cell;
                const int slice_0 = remaining % ((sclrpthlp_forcing.extent(0) - (1)) / (1) + 1);
                remaining /= ((sclrpthlp_forcing.extent(0) - (1)) / (1) + 1);
                const int slice_1 = remaining % ((sclrpthlp_forcing.extent(1) - (1)) / (1) + 1);
                remaining /= ((sclrpthlp_forcing.extent(1) - (1)) / (1) + 1);
                sclrpthlp_forcing((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                    thlp2_forcing((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)));
              });
          threshold = (thl_tol * thl_tol);
        } else {
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::array19",
              ((sclrpthlp_forcing.extent(0) - (1)) / (1) + 1) *
                  ((sclrpthlp_forcing.extent(1) - (1)) / (1) + 1),
              KOKKOS_LAMBDA(int cell) {
                int remaining     = cell;
                const int slice_0 = remaining % ((sclrpthlp_forcing.extent(0) - (1)) / (1) + 1);
                remaining /= ((sclrpthlp_forcing.extent(0) - (1)) / (1) + 1);
                const int slice_1 = remaining % ((sclrpthlp_forcing.extent(1) - (1)) / (1) + 1);
                remaining /= ((sclrpthlp_forcing.extent(1) - (1)) / (1) + 1);
                sclrpthlp_forcing((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = zero;
              });
          threshold = zero_threshold;
        }

        xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, // In
                     Kokkos::subview(lhs_ta_wpthlpsclrp, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                                     (sclr - 1)), // In
                     lhs_ma, lhs_diff, lhs_dp1,   // In
                     lhs);                        // Out

        xp2_xpyp_rhs(
            nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrpthlp, dt,                          // In
            Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), wpthlp, // In
            Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), thlm,     // In
            Kokkos::subview(sclrpthlp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),
            sclrpthlp_forcing,                      // In
            C2sclr_1d, invrs_tau_xp2_zm, threshold, // In
            Kokkos::subview(lhs_ta_wpthlpsclrp, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                            (sclr - 1)),                                               // In
            Kokkos::subview(rhs_ta_wpthlpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), // In
            stats,                                                                     // InOut
            rhs);                                                                      // Out

        // Solve the tridiagonal system
        xp2_xpyp_solve(
            nzm, ngrdcol, xp2_xpyp_scalars, 1, l_implemented, gr, tridiag_solve_method, // In
            stats, rhs, lhs, err_info,                                                  // InOut
            Kokkos::subview(sclrpthlp_solution, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out

        if (f2k::clubb_at_least_debug_level_api(0)) {
          if (any_equal(err_info.err_code, clubb_fatal_error)) {
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << err_info.err_header_global << "\n";
              std::cerr << message.str();
            }
            {
              std::ostringstream message;
              message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
              message << "calling xp2_xpyp_solve for sclrpthlp with sclr index " << " " << sclr
                      << " " << " in solve_xp2_xpyp_with_multiple_lhs" << "\n";
              std::cerr << message.str();
            }
            return;
          }
        }

        if (l_lmm_stepping) {
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::array20",
              ((sclrpthlp.extent(0) - (1)) / (1) + 1) * ((sclrpthlp.extent(1) - (1)) / (1) + 1),
              KOKKOS_LAMBDA(int cell) {
                int remaining     = cell;
                const int slice_0 = remaining % ((sclrpthlp.extent(0) - (1)) / (1) + 1);
                remaining /= ((sclrpthlp.extent(0) - (1)) / (1) + 1);
                const int slice_1 = remaining % ((sclrpthlp.extent(1) - (1)) / (1) + 1);
                remaining /= ((sclrpthlp.extent(1) - (1)) / (1) + 1);
                sclrpthlp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1)) =
                    one_half *
                    (sclrpthlp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1)) +
                     sclrpthlp_solution((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)),
                                        (sclr - 1)));
              });
        } else {
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::array21",
              ((sclrpthlp.extent(0) - (1)) / (1) + 1) * ((sclrpthlp.extent(1) - (1)) / (1) + 1),
              KOKKOS_LAMBDA(int cell) {
                int remaining     = cell;
                const int slice_0 = remaining % ((sclrpthlp.extent(0) - (1)) / (1) + 1);
                remaining /= ((sclrpthlp.extent(0) - (1)) / (1) + 1);
                const int slice_1 = remaining % ((sclrpthlp.extent(1) - (1)) / (1) + 1);
                remaining /= ((sclrpthlp.extent(1) - (1)) / (1) + 1);
                sclrpthlp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1)) =
                    sclrpthlp_solution((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)),
                                       (sclr - 1));
              });
        }
      }

    } else {

      // The ADG1 PDF is used.

      // Implicit contributions to passive scalars

      //!!!!***** sclr'^2, sclr'r_t', sclr'th_l' *****!!!!!
      // Note:  For ADG1, the LHS arrays are the same for all scalar variables,
      //        and also for <sclr'^2>, <sclr'r_t'>, and <sclr'th_l'>.
      xp2_xpyp_lhs(
          nzm, ngrdcol, dt, gr,                                                             // In
          Kokkos::subview(lhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, (1 - 1)), // In
          lhs_ma, lhs_diff, lhs_dp1,                                                        // In
          lhs);                                                                             // Out

      // Explicit contributions to passive scalars
      for (sclr = 1; sclr <= sclr_dim; sclr += 1) {

        // Forcing for <sclr'^2>.
        Kokkos::parallel_for(
            "solve_xp2_xpyp_with_multiple_lhs::loop22", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              sclrp2_forcing(i - 1, k - 1) = zero;
            });

        //!!!!***** sclr'^2 *****!!!!!
        xp2_xpyp_rhs(
            nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrp2, dt,                                       // In
            Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),                   // In
            Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),                   // In
            Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),                     // In
            Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),                     // In
            Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), sclrp2_forcing,    // In
            C2sclr_1d, invrs_tau_xp2_zm,                                                      // In
            (to_host_space(sclr_tol)(sclr - 1) * to_host_space(sclr_tol)(sclr - 1)),          // In
            Kokkos::subview(lhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, (1 - 1)), // In
            Kokkos::subview(rhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),           // In
            stats,                                                            // InOut
            Kokkos::subview(sclr_rhs, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out

        //!!!!***** sclr'r_t' *****!!!!!
        if (sclr == sclr_idx.iisclr_rt) {
          // In this case we're trying to emulate rt'^2 with sclr'rt', so we
          // handle this as we would a variance, even though generally speaking
          // the scalar is not rt
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::loop23", std::max(0, nzm) * std::max(0, ngrdcol),
              KOKKOS_LAMBDA(const int cell) {
                int remaining = cell;
                const int i   = 1 + remaining % (ngrdcol);
                remaining /= (ngrdcol);
                const int k = 1 + remaining % (nzm);
                remaining /= (nzm);
                sclrprtp_forcing(i - 1, k - 1) = rtp2_forcing(i - 1, k - 1);
              });

          threshold = (rt_tol * rt_tol);
        } else {
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::loop24", std::max(0, nzm) * std::max(0, ngrdcol),
              KOKKOS_LAMBDA(const int cell) {
                int remaining = cell;
                const int i   = 1 + remaining % (ngrdcol);
                remaining /= (ngrdcol);
                const int k = 1 + remaining % (nzm);
                remaining /= (nzm);
                sclrprtp_forcing(i - 1, k - 1) = zero;
              });

          threshold = zero_threshold;
        }

        xp2_xpyp_rhs(
            nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrprtp, dt,                                      // In
            Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), wprtp,             // In
            Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), rtm,                 // In
            Kokkos::subview(sclrprtp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), sclrprtp_forcing, // In
            C2sclr_1d, invrs_tau_xp2_zm, threshold,                                            // In
            Kokkos::subview(lhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, (1 - 1)),  // In
            Kokkos::subview(rhs_ta_wprtpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),          // In
            stats,                                                                       // InOut
            Kokkos::subview(sclr_rhs, Kokkos::ALL, Kokkos::ALL, (sclr + sclr_dim - 1))); // Out

        //!!!!***** sclr'th_l' *****!!!!!

        if (sclr == sclr_idx.iisclr_thl) {
          // In this case we're trying to emulate thl'^2 with sclr'thl', so we
          // handle this as we did with sclr_rt, above.
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::loop25", std::max(0, nzm) * std::max(0, ngrdcol),
              KOKKOS_LAMBDA(const int cell) {
                int remaining = cell;
                const int i   = 1 + remaining % (ngrdcol);
                remaining /= (ngrdcol);
                const int k = 1 + remaining % (nzm);
                remaining /= (nzm);
                sclrpthlp_forcing(i - 1, k - 1) = thlp2_forcing(i - 1, k - 1);
              });

          threshold = (thl_tol * thl_tol);
        } else {
          Kokkos::parallel_for(
              "solve_xp2_xpyp_with_multiple_lhs::loop26", std::max(0, nzm) * std::max(0, ngrdcol),
              KOKKOS_LAMBDA(const int cell) {
                int remaining = cell;
                const int i   = 1 + remaining % (ngrdcol);
                remaining /= (ngrdcol);
                const int k = 1 + remaining % (nzm);
                remaining /= (nzm);
                sclrpthlp_forcing(i - 1, k - 1) = zero;
              });

          threshold = zero_threshold;
        }

        xp2_xpyp_rhs(
            nzm, nzt, ngrdcol, gr, xp2_xpyp_sclrpthlp, dt,                          // In
            Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), wpthlp, // In
            Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, (sclr - 1)), thlm,     // In
            Kokkos::subview(sclrpthlp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),
            sclrpthlp_forcing,                                                                // In
            C2sclr_1d, invrs_tau_xp2_zm, threshold,                                           // In
            Kokkos::subview(lhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, (1 - 1)), // In
            Kokkos::subview(rhs_ta_wpthlpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),        // In
            stats, // InOut
            Kokkos::subview(sclr_rhs, Kokkos::ALL, Kokkos::ALL, (sclr + 2 * sclr_dim - 1))); // Out
      }

      // Solve the tridiagonal system
      xp2_xpyp_solve(nzm, ngrdcol, xp2_xpyp_scalars, 3 * sclr_dim, l_implemented, gr, // In
                     tridiag_solve_method,                                            // In
                     stats, sclr_rhs, lhs, err_info,                                  // InOut
                     sclr_solution);                                                  // Out

      if (f2k::clubb_at_least_debug_level_api(0)) {
        if (any_equal(err_info.err_code, clubb_fatal_error)) {
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << err_info.err_header_global << "\n";
            std::cerr << message.str();
          }
          {
            std::ostringstream message;
            message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
            message << "calling xp2_xpyp_solve for sclr with sclr index " << " " << sclr << " "
                    << " in solve_xp2_xpyp_with_multiple_lhs" << "\n";
            std::cerr << message.str();
          }
          return;
        }
      }

      if (l_lmm_stepping) {

        Kokkos::parallel_for(
            "solve_xp2_xpyp_with_multiple_lhs::loop27",
            std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              const int sclr = 1 + remaining % (sclr_dim);
              remaining /= (sclr_dim);
              sclrp2(i - 1, k - 1, sclr - 1) = one_half * (sclrp2(i - 1, k - 1, sclr - 1) +
                                                           sclr_solution(i - 1, k - 1, sclr - 1));
              sclrprtp(i - 1, k - 1, sclr - 1) =
                  one_half * (sclrprtp(i - 1, k - 1, sclr - 1) +
                              sclr_solution(i - 1, k - 1, sclr_dim + sclr - 1));
              sclrpthlp(i - 1, k - 1, sclr - 1) =
                  one_half * (sclrpthlp(i - 1, k - 1, sclr - 1) +
                              sclr_solution(i - 1, k - 1, 2 * sclr_dim + sclr - 1));
            });

      } else {

        Kokkos::parallel_for(
            "solve_xp2_xpyp_with_multiple_lhs::loop28",
            std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              const int sclr = 1 + remaining % (sclr_dim);
              remaining /= (sclr_dim);
              sclrp2(i - 1, k - 1, sclr - 1)   = sclr_solution(i - 1, k - 1, sclr - 1);
              sclrprtp(i - 1, k - 1, sclr - 1) = sclr_solution(i - 1, k - 1, sclr_dim + sclr - 1);
              sclrpthlp(i - 1, k - 1, sclr - 1) =
                  sclr_solution(i - 1, k - 1, 2 * sclr_dim + sclr - 1);
            });
      }
    }
  }

  return;
}

// lhs_ta: Turbulent advection contributions to lhs
// lhs_ma: Mean advection contributions to lhs
// lhs_diff: Diffusion contributions to lhs, dissipation term 2
// lhs_dp1: LHS dissipation term 1
// lhs: Implicit contributions to the term
template <class Array0, class Array1, class Array2, class Array3, class Array4>
void xp2_xpyp_lhs(int nzm, int ngrdcol, real_t dt, const grid& gr,               // In
                  Array0 lhs_ta, Array1 lhs_ma, Array2 lhs_diff, Array3 lhs_dp1, // In
                  Array4 lhs)                                                    // Out
{
  // Description:
  //     Compute LHS tridiagonal matrix for a variance or covariance term
  //
  // Notes:
  //     For LHS turbulent advection terms:
  //         An "over-implicit" weighted time step is applied to this term.
  //         The weight of the implicit portion of this term is controlled by
  //         the factor gamma_over_implicit_ts (abbreviated "gamma" in the
  //         expression below).  A factor is added to the right-hand side of
  //         the equation in order to balance a weight that is not equal to 1,
  //         such that:
  //              -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
  //         where X is the variable that is being solved for in a predictive
  //         equation (<x'^2> or <x'y'> in this case), y(t) is the linearized
  //         portion of the term that gets treated implicitly, and RHS is the
  //         portion of the term that is always treated explicitly.  A weight
  //         of greater than 1 can be applied to make the term more
  //         numerically stable.
  //
  //----------------------------------------------------------------------------------
  //------------------- Input Variables -------------------
  //------------------- Output Variables -------------------
  //---------------- Local Variables -------------------
  // Array indices
  int i;
  int k;
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;

  //---------------- Begin Code -------------------

  // Combine all lhs terms into lhs, should be fully vectorized
  Kokkos::parallel_for(
      "xp2_xpyp_lhs::loop1", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        lhs(1 - 1, i - 1, k - 1) = lhs_diff(1 - 1, i - 1, k - 1) + lhs_ma(1 - 1, i - 1, k - 1) +
                                   lhs_ta(1 - 1, i - 1, k - 1) * gamma_over_implicit_ts;

        lhs(2 - 1, i - 1, k - 1) = lhs_diff(2 - 1, i - 1, k - 1) + lhs_ma(2 - 1, i - 1, k - 1) +
                                   lhs_ta(2 - 1, i - 1, k - 1) * gamma_over_implicit_ts +
                                   lhs_dp1(i - 1, k - 1);

        lhs(3 - 1, i - 1, k - 1) = lhs_diff(3 - 1, i - 1, k - 1) + lhs_ma(3 - 1, i - 1, k - 1) +
                                   lhs_ta(3 - 1, i - 1, k - 1) * gamma_over_implicit_ts;
      });

  // LHS time tendency.
  Kokkos::parallel_for(
      "xp2_xpyp_lhs::loop2", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        lhs(2 - 1, i - 1, k - 1) = lhs(2 - 1, i - 1, k - 1) + (one / dt);
      });

  // Boundary Conditions
  // These are set so that the surface_varnce value of the variances and
  // covariances can be used at the lowest boundary and the values of those
  // variables can be set to their respective threshold minimum values at the
  // top boundary.  Fixed-point boundary conditions are used for both the
  // variances and the covariances.

  Kokkos::parallel_for(
      "xp2_xpyp_lhs::loop3", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        // Set the lower boundaries for the first variable
        lhs(1 - 1, i - 1, gr_k_lb_zm - 1) = zero;
        lhs(2 - 1, i - 1, gr_k_lb_zm - 1) = one;
        lhs(3 - 1, i - 1, gr_k_lb_zm - 1) = zero;

        // Set the upper boundaries for the first variable
        lhs(1 - 1, i - 1, gr_k_ub_zm - 1) = zero;
        lhs(2 - 1, i - 1, gr_k_ub_zm - 1) = one;
        lhs(3 - 1, i - 1, gr_k_ub_zm - 1) = zero;
      });

  return;
}

// solve_type: Variable(s) description
// nrhs: Number of right hand side vectors
// l_implemented: True if CLUBB is being implemented and run in a host model
// tridiag_solve_method: Specifier for method to solve tridiagonal systems
// rhs: Explicit contributions to x variance/covariance term [units vary]
// lhs: Implicit contributions to x variance/covariance term [units vary]
// err_info: err_info struct containing err_code and err_header
// xapxbp: Computed value of the variable(s) at <t+1> [units vary]
template <class Array0, class Array1, class Array2>
void xp2_xpyp_solve(int nzm, int ngrdcol, int solve_type, int nrhs,  // In
                    bool l_implemented,                              // In
                    const grid& gr, int tridiag_solve_method,        // In
                    f2k::stats_type& stats,                          // InOut
                    Array0 rhs, Array1 lhs, err_info_type& err_info, // InOut
                    Array2 xapxbp)                                   // Out
{
  // Description:
  // Solve a tridiagonal system
  //
  // References:
  //   None
  //-----------------------------------------------------------------------
  // External
  // ---------------------- Input variables ----------------------
  //------------------- Input/Output variables ------------------
  // ---------------------- Input/Ouput variables ----------------------
  // ---------------------- Output Variables ----------------------
  // ---------------------- Local variables ----------------------
  Kokkos::View<real_t*, Kokkos::LayoutLeft> rcond(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_solve::rcond"),
      ngrdcol); // Est. of the reciprocal of the condition # on the matrix
  bool l_need_rcond;
  bool l_single_lhs_solve;
  std::string solve_type_str; // solve_type in string format for debug output purposes

  // --- Begin Code ---

  l_single_lhs_solve = false;
  l_need_rcond       = false;
  //------------------------------------------------------------------------
  // Note that these are diagnostics from inverting the matrix, not a budget

  switch (solve_type) {
  case xp2_xpyp_rtp2: {
    solve_type_str = "rtp2";
    l_need_rcond = stats.flags().l_sample && f2k::var_on_stats_list(stats, "rtp2_matrix_condt_num");
    break;
  }
  case xp2_xpyp_thlp2: {
    solve_type_str = "thlp2";
    l_need_rcond =
        stats.flags().l_sample && f2k::var_on_stats_list(stats, "thlp2_matrix_condt_num");
    break;
  }
  case xp2_xpyp_rtpthlp: {
    solve_type_str = "rtpthlp";
    l_need_rcond =
        stats.flags().l_sample && f2k::var_on_stats_list(stats, "rtpthlp_matrix_condt_num");
    break;
  }
  case xp2_xpyp_up2_vp2: {
    solve_type_str = "up2_vp2";
    l_need_rcond =
        stats.flags().l_sample && f2k::var_on_stats_list(stats, "up2_vp2_matrix_condt_num");
    break;
  }
  case xp2_xpyp_single_lhs: {
    l_single_lhs_solve = true;
    solve_type_str     = "xp2_xpyp_single_lhs";
    // In single solve tpype, condition number is either output for none or all
    // rtp2, thlp2, and rtpthlp together
    l_need_rcond =
        stats.flags().l_sample && (f2k::var_on_stats_list(stats, "rtp2_matrix_condt_num") ||
                                   f2k::var_on_stats_list(stats, "thlp2_matrix_condt_num") ||
                                   f2k::var_on_stats_list(stats, "rtpthlp_matrix_condt_num"));
    break;
  }
  default: {
    solve_type_str = "scalar";
    l_need_rcond   = false;
    break;
  }
  }

  // Matrix solves are bit-different between ascending and descending.
  // This ensures matrices are solved in the same (descending) order,
  // which is useful for ensuring BFBness between grid modes
  // Native lhs is (band,column,row); the existing solver owns (column,band,row).
  // Transpose only at this matrix API boundary; keep Native descending flips.
  const bool reverse = l_force_descending_solves && gr.grid_dir_indx > 0;
  Real3D lhs_solver("xp2_xpyp_solve::lhs_solver", ngrdcol, ndiags3, nzm);
  Real3D rhs_solver("xp2_xpyp_solve::rhs_solver", ngrdcol, nzm, nrhs);
  Real3D solution_solver("xp2_xpyp_solve::solution_solver", ngrdcol, nzm, nrhs);
  Kokkos::parallel_for(
      "xp2_xpyp_solve::matrix_order", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol, source_k = reverse ? nzm - 1 - k : k;
        for (int b = 0; b < ndiags3; ++b)
          lhs_solver(i, b, k) = lhs(reverse ? ndiags3 - 1 - b : b, i, source_k);
        for (int n = 0; n < nrhs; ++n) {
          if constexpr (decltype(rhs)::rank == 2)
            rhs_solver(i, k, n) = rhs(i, source_k);
          else
            rhs_solver(i, k, n) = rhs(i, source_k, n);
        }
      });

  if (l_need_rcond) {

    tridiag_solve(solve_type_str, tridiag_solve_method, ngrdcol, nzm, nrhs, l_implemented,
                  lhs_solver, rhs_solver, err_info, solution_solver, &rcond);
    // Restore Native InOut allocations and the unflipped result before diagnostics.
    Kokkos::parallel_for(
        "xp2_xpyp_solve::native_storage", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          for (int b = 0; b < ndiags3; ++b)
            lhs(b, i, k) = lhs_solver(i, b, k);
          for (int n = 0; n < nrhs; ++n) {
            if constexpr (decltype(rhs)::rank == 2)
              rhs(i, k) = rhs_solver(i, k, n);
            else
              rhs(i, k, n) = rhs_solver(i, k, n);
            if constexpr (decltype(xapxbp)::rank == 2)
              xapxbp(i, k) = solution_solver(i, k, n);
            else
              xapxbp(i, k, n) = solution_solver(i, k, n);
          }
        });

    if (stats.flags().l_sample) {
      if (l_single_lhs_solve) {
        // Single lhs solve including rtp2, thlp2, rtpthlp. Estimate for each.
        // Est. of the condition number of the variance LHS matrix
        {
          Real1D stats_values_1(
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_solve::stats_values_1"),
              rcond.extent(0));
          Kokkos::parallel_for(
              "xp2_xpyp_solve::stats_value", rcond.size(), KOKKOS_LAMBDA(int cell) {
                const int i       = cell;
                stats_values_1(i) = one / rcond(i);
              });
          f2k::stats_update("rtp2_matrix_condt_num", f2k::Real1D(stats_values_1), stats);
        }
        {
          Real1D stats_values_2(
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_solve::stats_values_2"),
              rcond.extent(0));
          Kokkos::parallel_for(
              "xp2_xpyp_solve::stats_value", rcond.size(), KOKKOS_LAMBDA(int cell) {
                const int i       = cell;
                stats_values_2(i) = one / rcond(i);
              });
          f2k::stats_update("thlp2_matrix_condt_num", f2k::Real1D(stats_values_2), stats);
        }
        {
          Real1D stats_values_3(
              Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_solve::stats_values_3"),
              rcond.extent(0));
          Kokkos::parallel_for(
              "xp2_xpyp_solve::stats_value", rcond.size(), KOKKOS_LAMBDA(int cell) {
                const int i       = cell;
                stats_values_3(i) = one / rcond(i);
              });
          f2k::stats_update("rtpthlp_matrix_condt_num", f2k::Real1D(stats_values_3), stats);
        }
      } else {
        switch (solve_type) {
        case xp2_xpyp_rtp2: {
          {
            Real1D stats_values_4(
                Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_solve::stats_values_4"),
                rcond.extent(0));
            Kokkos::parallel_for(
                "xp2_xpyp_solve::stats_value", rcond.size(), KOKKOS_LAMBDA(int cell) {
                  const int i       = cell;
                  stats_values_4(i) = one / rcond(i);
                });
            f2k::stats_update("rtp2_matrix_condt_num", f2k::Real1D(stats_values_4), stats);
          }
          break;
        }
        case xp2_xpyp_thlp2: {
          {
            Real1D stats_values_5(
                Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_solve::stats_values_5"),
                rcond.extent(0));
            Kokkos::parallel_for(
                "xp2_xpyp_solve::stats_value", rcond.size(), KOKKOS_LAMBDA(int cell) {
                  const int i       = cell;
                  stats_values_5(i) = one / rcond(i);
                });
            f2k::stats_update("thlp2_matrix_condt_num", f2k::Real1D(stats_values_5), stats);
          }
          break;
        }
        case xp2_xpyp_rtpthlp: {
          {
            Real1D stats_values_6(
                Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_solve::stats_values_6"),
                rcond.extent(0));
            Kokkos::parallel_for(
                "xp2_xpyp_solve::stats_value", rcond.size(), KOKKOS_LAMBDA(int cell) {
                  const int i       = cell;
                  stats_values_6(i) = one / rcond(i);
                });
            f2k::stats_update("rtpthlp_matrix_condt_num", f2k::Real1D(stats_values_6), stats);
          }
          break;
        }
        case xp2_xpyp_up2_vp2: {
          {
            Real1D stats_values_7(
                Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_solve::stats_values_7"),
                rcond.extent(0));
            Kokkos::parallel_for(
                "xp2_xpyp_solve::stats_value", rcond.size(), KOKKOS_LAMBDA(int cell) {
                  const int i       = cell;
                  stats_values_7(i) = one / rcond(i);
                });
            f2k::stats_update("up2_vp2_matrix_condt_num", f2k::Real1D(stats_values_7), stats);
          }
          break;
        }
        default: {
          // No condition number is setup for the passive scalars
          // No condition number output for passive scalars.
          break;
        }
        }
      }
    }

  } else {

    tridiag_solve(solve_type_str, tridiag_solve_method, ngrdcol, nzm, nrhs, l_implemented,
                  lhs_solver, rhs_solver, err_info, solution_solver, nullptr);
    // Restore Native InOut allocations and the unflipped result before diagnostics.
    Kokkos::parallel_for(
        "xp2_xpyp_solve::native_storage", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          for (int b = 0; b < ndiags3; ++b)
            lhs(b, i, k) = lhs_solver(i, b, k);
          for (int n = 0; n < nrhs; ++n) {
            if constexpr (decltype(rhs)::rank == 2)
              rhs(i, k) = rhs_solver(i, k, n);
            else
              rhs(i, k, n) = rhs_solver(i, k, n);
            if constexpr (decltype(xapxbp)::rank == 2)
              xapxbp(i, k) = solution_solver(i, k, n);
            else
              xapxbp(i, k, n) = solution_solver(i, k, n);
          }
        });
  }

  if (f2k::clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << err_info.err_header_global << "\n";
        std::cerr << message.str();
      }
      {
        std::ostringstream message;
        message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        message << "calling tridiag_solve in xp2_xpyp_solve" << "\n";
        std::cerr << message.str();
      }
      return;
    }
  }

  // Flip the back to the ascending direction if we forced the solve
  // to be in descending mode
  // Flip the computed answer only after the Native fatal-return check.
  if (l_force_descending_solves && gr.grid_dir_indx > 0) {
    Kokkos::parallel_for(
        "xp2_xpyp_solve::restore_direction", ngrdcol * nzm, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          for (int n = 0; n < nrhs; ++n) {
            if constexpr (decltype(xapxbp)::rank == 2)
              xapxbp(i, k) = solution_solver(i, nzm - 1 - k, n);
            else
              xapxbp(i, k, n) = solution_solver(i, nzm - 1 - k, n);
          }
        });
  }

  return;
}

template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
void stats_finalize_xp2_xpyp_terms(int nzm, int ngrdcol, int solve_type, const grid& gr,
                                   Array0 xapxbp, // In
                                   Array1 lhs_dp1, Array2 lhs_dp1_pr1, Array3 lhs_diff,
                                   Array4 lhs_ta,          // In
                                   Array5 lhs_ma,          // In
                                   f2k::stats_type& stats) // InOut
{
  // Finalize the field-specific implicit budget terms from the solved field and LHS components.
  // Inputs retain momentum levels/Native bands; mutate only the Native stats context.
  //----------------------- Input variables -----------------------
  //----------------------- InOut variables -----------------------
  //----------------------- Local Variables -----------------------
  // Budget indices
  std::string name_dp1;
  std::string name_dp2;
  std::string name_ta;
  std::string name_ma;
  std::string name_pr1;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_tmp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "stats_finalize_xp2_xpyp_terms::stats_tmp"),
      ngrdcol, nzm);
  int i;
  int k;
  int km1;
  int kp1;
  const auto gr_grid_dir_indx = gr.grid_dir_indx;

  //-----------------------Begin Code -----------------------

  if (!stats.flags().l_sample) {
    return;
  }

  // Finalize implicit contributions for r_t'^2, th_l'^2, r_t'th_l',
  // u'^2, and v'^2.

  switch (solve_type) {
  case xp2_xpyp_rtp2: {
    name_dp1 = "rtp2_dp1";
    name_dp2 = "rtp2_dp2";
    name_ta  = "rtp2_ta";
    name_ma  = "rtp2_ma";
    name_pr1 = "";
    break;
  }
  case xp2_xpyp_thlp2: {
    name_dp1 = "thlp2_dp1";
    name_dp2 = "thlp2_dp2";
    name_ta  = "thlp2_ta";
    name_ma  = "thlp2_ma";
    name_pr1 = "";
    break;
  }
  case xp2_xpyp_rtpthlp: {
    name_dp1 = "rtpthlp_dp1";
    name_dp2 = "rtpthlp_dp2";
    name_ta  = "rtpthlp_ta";
    name_ma  = "rtpthlp_ma";
    name_pr1 = "";
    break;
  }
  case xp2_xpyp_up2: {
    name_dp1 = "up2_dp1";
    name_dp2 = "up2_dp2";
    name_ta  = "up2_ta";
    name_ma  = "up2_ma";
    name_pr1 = "up2_pr1";
    break;
  }
  case xp2_xpyp_vp2: {
    name_dp1 = "vp2_dp1";
    name_dp2 = "vp2_dp2";
    name_ta  = "vp2_ta";
    name_ma  = "vp2_ma";
    name_pr1 = "vp2_pr1";
    break;
  }
  default: {
    return;
    break;
  }
  }

  Kokkos::parallel_for(
      "stats_finalize_xp2_xpyp_terms::array1",
      ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
      KOKKOS_LAMBDA(int cell) {
        int remaining     = cell;
        const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
        remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
        const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
        remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
        stats_tmp(slice_0, slice_1) = zero;
      });
  Kokkos::parallel_for(
      "stats_finalize_xp2_xpyp_terms::loop2", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        stats_tmp(i - 1, k - 1) = -lhs_dp1(i - 1, k - 1) * xapxbp(i - 1, k - 1);
      });

  // x'y' term dp1 has both implicit and explicit components;
  // call stat_end_update_pt.
  f2k::stats_finalize_budget(name_dp1, stats_tmp, // In
                             stats);              // InOut

  if (solve_type == xp2_xpyp_up2 || solve_type == xp2_xpyp_vp2) {
    Kokkos::parallel_for(
        "stats_finalize_xp2_xpyp_terms::array3",
        ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
        KOKKOS_LAMBDA(int cell) {
          int remaining     = cell;
          const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
          remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
          const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
          remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
          stats_tmp(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "stats_finalize_xp2_xpyp_terms::loop4", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (nzm - 2);
          remaining /= (nzm - 2);
          stats_tmp(i - 1, k - 1) = -lhs_dp1_pr1(i - 1, k - 1) * xapxbp(i - 1, k - 1);
        });

    // x'y' term pr1 has both implicit and explicit components;
    // call stat_end_update_pt.
    f2k::stats_finalize_budget(name_pr1, stats_tmp, // In
                               stats);              // InOut
  }

  Kokkos::parallel_for(
      "stats_finalize_xp2_xpyp_terms::array5",
      ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
      KOKKOS_LAMBDA(int cell) {
        int remaining     = cell;
        const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
        remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
        const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
        remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
        stats_tmp(slice_0, slice_1) = zero;
      });
  Kokkos::parallel_for(
      "stats_finalize_xp2_xpyp_terms::loop6", std::max(0, nzm - 2), KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int k   = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        int i;
        int km1;
        int kp1;
        km1 = ((k - gr_grid_dir_indx) > (1) ? (k - gr_grid_dir_indx) : (1));
        kp1 = ((k + gr_grid_dir_indx) < (nzm) ? (k + gr_grid_dir_indx) : (nzm));
        for (i = 1; i <= ngrdcol; i += 1) {
          stats_tmp(i - 1, k - 1) =
              (-gamma_over_implicit_ts * lhs_ta(2 + gr_grid_dir_indx - 1, i - 1, k - 1)) *
                  xapxbp(i - 1, km1 - 1) +
              (-gamma_over_implicit_ts * lhs_ta(2 - 1, i - 1, k - 1)) * xapxbp(i - 1, k - 1) +
              (-gamma_over_implicit_ts * lhs_ta(2 - gr_grid_dir_indx - 1, i - 1, k - 1)) *
                  xapxbp(i - 1, kp1 - 1);
        }
      });

  // x'y' term ta has both implicit and explicit components;
  f2k::stats_finalize_budget(name_ta, stats_tmp, // In
                             stats);             // InOut

  Kokkos::parallel_for(
      "stats_finalize_xp2_xpyp_terms::array7",
      ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
      KOKKOS_LAMBDA(int cell) {
        int remaining     = cell;
        const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
        remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
        const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
        remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
        stats_tmp(slice_0, slice_1) = zero;
      });
  Kokkos::parallel_for(
      "stats_finalize_xp2_xpyp_terms::loop8", std::max(0, nzm - 2), KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int k   = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        int i;
        int km1;
        int kp1;
        km1 = ((k - gr_grid_dir_indx) > (1) ? (k - gr_grid_dir_indx) : (1));
        kp1 = ((k + gr_grid_dir_indx) < (nzm) ? (k + gr_grid_dir_indx) : (nzm));
        for (i = 1; i <= ngrdcol; i += 1) {
          stats_tmp(i - 1, k - 1) =
              (-lhs_diff(2 + gr_grid_dir_indx - 1, i - 1, k - 1)) * xapxbp(i - 1, km1 - 1) +
              (-lhs_diff(2 - 1, i - 1, k - 1)) * xapxbp(i - 1, k - 1) +
              (-lhs_diff(2 - gr_grid_dir_indx - 1, i - 1, k - 1)) * xapxbp(i - 1, kp1 - 1);
        }
      });

  // x'y' term dp2 is completely implicit; call stat_update_var_pt.
  f2k::stats_update(name_dp2, stats_tmp, // In
                    stats);              // InOut

  Kokkos::parallel_for(
      "stats_finalize_xp2_xpyp_terms::array9",
      ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
      KOKKOS_LAMBDA(int cell) {
        int remaining     = cell;
        const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
        remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
        const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
        remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
        stats_tmp(slice_0, slice_1) = zero;
      });
  Kokkos::parallel_for(
      "stats_finalize_xp2_xpyp_terms::loop10", std::max(0, nzm - 2), KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int k   = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        int i;
        int km1;
        int kp1;
        km1 = ((k - gr_grid_dir_indx) > (1) ? (k - gr_grid_dir_indx) : (1));
        kp1 = ((k + gr_grid_dir_indx) < (nzm) ? (k + gr_grid_dir_indx) : (nzm));
        for (i = 1; i <= ngrdcol; i += 1) {
          stats_tmp(i - 1, k - 1) =
              (-lhs_ma(2 + gr_grid_dir_indx - 1, i - 1, k - 1)) * xapxbp(i - 1, km1 - 1) +
              (-lhs_ma(2 - 1, i - 1, k - 1)) * xapxbp(i - 1, k - 1) +
              (-lhs_ma(2 - gr_grid_dir_indx - 1, i - 1, k - 1)) * xapxbp(i - 1, kp1 - 1);
        }
      });

  // x'y' term ma is completely implicit; call stat_update_var_pt.
  f2k::stats_update(name_ma, stats_tmp, // In
                    stats);             // InOut

  return;
}

// wp2: w'^2 (momentum levels)                      [m^2/s^2]
// wpthvp: w'th_v' (momentum levels)                     [K m/s]
// invrs_tau_C4_zm: Inverse time-scale for C4 terms on mom. levels  [1/s]
// invrs_tau_C14_zm: Inverse time-scale for C14 terms on mom. levels [1/s]
// xam: x_am (thermodynamic levels)                     [m/s]
// xbm: x_bm (thermodynamic levels)                     [m/s]
// wpxap: w'x_a' (momentum levels)                    [m^2/s^2]
// wpxbp: w'x_b' (momentum levels)                    [m^2/s^2]
// xap2: x_a'^2 (momentum levels)                    [m^2/s^2]
// xbp2: x_b'^2 (momentum levels)                    [m^2/s^2]
// thv_ds_zm: Dry, base-state theta_v on momentum levs.         [K]
// C4: Model parameter C_4                         [-]
// C_uu_shr: Model parameter C_uu_shr                    [-]
// C_uu_buoy: Model parameter C_uu_buoy                   [-]
// C14: Model parameter C_{14}                      [-]
// lhs_splat_wp2: LHS coefficient of wp2 splatting term           [1/s]
// lhs_ta: LHS turbulent advection term
// rhs_ta: RHS turbulent advection terms
// lhs_dp1_C4: LHS dissipation term 1, for up2 vp2 using C14
// lhs_dp1_C14: LHS dissipation term 1, for up2 vp2 using C4
// rhs: Explicit contributions to x variance/covariance terms
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20>
void xp2_xpyp_uv_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, int solve_type, real_t dt, // In
                     Array0 wp2, Array1 wpthvp,                                                // In
                     Array2 invrs_tau_C4_zm, Array3 invrs_tau_C14_zm,                          // In
                     Array4 xam, Array5 xbm, Array6 wpxap, Array7 wpxbp, Array8 xap2,
                     Array9 xbp2,                                                        // In
                     Array10 thv_ds_zm, Array11 C4, Array12 C_uu_shr, Array13 C_uu_buoy, // In
                     Array14 C14, Array15 lhs_splat_wp2,                                 // In
                     Array16 lhs_ta, Array17 rhs_ta,                                     // In
                     Array18 lhs_dp1_C4, Array19 lhs_dp1_C14,                            // In
                     f2k::stats_type& stats,                                             // InOut
                     Array20 rhs)                                                        // Out
{
  // Description:
  //     Explicit contributions to u'^2 or v'^2
  //
  // Notes:
  //    For LHS turbulent advection (ta) term:
  //        An "over-implicit" weighted time step is applied to this term.
  //        The weight of the implicit portion of this term is controlled by
  //        the factor gamma_over_implicit_ts (abbreviated "gamma" in the
  //        expression below).  A factor is added to the right-hand side of
  //        the equation in order to balance a weight that is not equal to 1,
  //        such that:
  //             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
  //        where X is the variable that is being solved for in a predictive
  //        equation (x'^2 or x'y' in this case), y(t) is the linearized
  //        portion of the term that gets treated implicitly, and RHS is the
  //        portion of the term that is always treated explicitly.  A weight
  //        of greater than 1 can be applied to make the term more
  //        numerically stable.
  //
  //
  //   --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
  //   Significant changes to this routine may adversely affect computational speed
  //       - Gunther Huebler, Aug. 2018, clubb:ticket:834
  //----------------------------------------------------------------------------------
  //------------------- Input Variables -------------------
  // For "over-implicit" weighted time step.
  // This vector holds output from the LHS (implicit) portion of a term at a
  // given vertical level.  This output is weighted and applied to the RHS.
  // This is used if the implicit portion of the term is "over-implicit", which
  // means that the LHS contribution is given extra weight (>1) in order to
  // increase numerical stability.  A weighted factor must then be applied to
  // the RHS in order to balance the weight.
  //------------------- Output Variables -------------------
  //---------------- Local Variables -------------------
  // Array indices
  int i;
  int k;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_pr1(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_uv_rhs::rhs_pr1"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_pr2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_uv_rhs::rhs_pr2"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_pr1(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_uv_rhs::stats_pr1"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_pr2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_uv_rhs::stats_pr2"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_term_tp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_uv_rhs::rhs_term_tp"), ngrdcol,
      nzm);
  Kokkos::View<real_t*, Kokkos::LayoutLeft> C4_zeros(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_uv_rhs::C4_zeros"), ngrdcol);
  Kokkos::View<real_t*, Kokkos::LayoutLeft> C14_zeros(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_uv_rhs::C14_zeros"), ngrdcol);
  std::string name_ta;
  std::string name_tp;
  std::string name_dp1;
  std::string name_pr1;
  std::string name_pr2;
  std::string name_splat;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_zm_begin(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_uv_rhs::stats_zm_begin"), ngrdcol,
      nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_zm_mod(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_uv_rhs::stats_zm_mod"), ngrdcol,
      nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_tmp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_uv_rhs::stats_tmp"), ngrdcol, nzm);
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_lb_zm       = gr.k_lb_zm;
  const auto gr_k_ub_zm       = gr.k_ub_zm;

  //----------------------------- Begin Code ----------------------------------

  switch (solve_type) {
  case xp2_xpyp_vp2: {
    name_ta    = "vp2_ta";
    name_tp    = "vp2_tp";
    name_dp1   = "vp2_dp1";
    name_pr1   = "vp2_pr1";
    name_pr2   = "vp2_pr2";
    name_splat = "vp2_splat";
    break;
  }
  case xp2_xpyp_up2: {
    name_ta    = "up2_ta";
    name_tp    = "up2_tp";
    name_dp1   = "up2_dp1";
    name_pr1   = "up2_pr1";
    name_pr2   = "up2_pr2";
    name_splat = "up2_splat";
    break;
  }
  default: {
    name_ta    = "";
    name_tp    = "";
    name_dp1   = "";
    name_pr1   = "";
    name_pr2   = "";
    name_splat = "";
    break;
  }
  }

  // Vertical compression of eddies causes gustiness (increase in up2 and vp2)
  // Add half the contribution to up2 and half to vp2
  Kokkos::parallel_for(
      "xp2_xpyp_uv_rhs::loop1", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        rhs(i - 1, k - 1) =
            rhs_ta(i - 1, k - 1) + real_t(0.5) * lhs_splat_wp2(i - 1, k - 1) * wp2(i - 1, k - 1);
      });

  // Calculate RHS pressure term 1 (pr1).
  term_pr1(nzm, ngrdcol, C4, C14, xbp2, wp2, invrs_tau_C4_zm, invrs_tau_C14_zm, // In
           rhs_pr1);                                                            // Out

  // Calculate RHS pressure term 2 (pr2).
  term_pr2(nzm, nzt, ngrdcol, gr, C_uu_shr, C_uu_buoy, thv_ds_zm, wpthvp, wpxap, wpxbp, xam, // In
           xbm,                                                                              // In
           rhs_pr2);                                                                         // Out

  // RHS turbulent production (tp) term.
  // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:up2_pr
  term_tp_rhs(nzm, nzt, ngrdcol, xam, xam, wpxap, wpxap, gr.invrs_dzm, // In
              rhs_term_tp);                                            // Out

  // Finish RHS calc with vectorizable loop, functions are in source file and should
  // be inlined with an -O2 or above compiler optimization flag
  Kokkos::parallel_for(
      "xp2_xpyp_uv_rhs::loop2", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);

        rhs(i - 1, k - 1) =
            rhs(i - 1, k - 1) +
            (one - gamma_over_implicit_ts) * (-lhs_ta(2 - gr_grid_dir_indx - 1, i - 1, k - 1) *
                                                  xap2(i - 1, k + gr_grid_dir_indx - 1) -
                                              lhs_ta(2 - 1, i - 1, k - 1) * xap2(i - 1, k - 1) -
                                              lhs_ta(2 + gr_grid_dir_indx - 1, i - 1, k - 1) *
                                                  xap2(i - 1, k - gr_grid_dir_indx - 1));

        // RHS turbulent production (tp) term.
        // https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:up2_pr
        rhs(i - 1, k - 1) = rhs(i - 1, k - 1) + (one - C_uu_shr(i - 1)) * rhs_term_tp(i - 1, k - 1);

        // RHS pressure term 1 (pr1) (and dissipation term 1 (dp1)).
        rhs(i - 1, k - 1) = rhs(i - 1, k - 1) + rhs_pr1(i - 1, k - 1);

        // RHS contribution from "over-implicit" weighted time step
        // for LHS dissipation term 1 (dp1) and pressure term 1 (pr1).
        rhs(i - 1, k - 1) =
            rhs(i - 1, k - 1) + (one - gamma_over_implicit_ts) *
                                    (-lhs_dp1_C4(i - 1, k - 1) - lhs_dp1_C14(i - 1, k - 1)) *
                                    xap2(i - 1, k - 1);

        // RHS pressure term 2 (pr2).
        rhs(i - 1, k - 1) = rhs(i - 1, k - 1) + rhs_pr2(i - 1, k - 1);
      });

  // RHS time tendency.
  Kokkos::parallel_for(
      "xp2_xpyp_uv_rhs::loop3", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        rhs(i - 1, k - 1) = rhs(i - 1, k - 1) + one / dt * xap2(i - 1, k - 1);
      });

  if (stats.flags().l_sample) {

    Kokkos::parallel_for(
        "xp2_xpyp_uv_rhs::array4", ((C14_zeros.extent(0) - (1)) / (1) + 1),
        KOKKOS_LAMBDA(int cell) {
          int remaining     = cell;
          const int slice_0 = remaining % ((C14_zeros.extent(0) - (1)) / (1) + 1);
          remaining /= ((C14_zeros.extent(0) - (1)) / (1) + 1);
          C14_zeros(slice_0) = zero;
        });
    Kokkos::parallel_for(
        "xp2_xpyp_uv_rhs::array5", ((C4_zeros.extent(0) - (1)) / (1) + 1), KOKKOS_LAMBDA(int cell) {
          int remaining     = cell;
          const int slice_0 = remaining % ((C4_zeros.extent(0) - (1)) / (1) + 1);
          remaining /= ((C4_zeros.extent(0) - (1)) / (1) + 1);
          C4_zeros(slice_0) = zero;
        });

    term_pr1(nzm, ngrdcol, C4, C14_zeros, xbp2, wp2, invrs_tau_C4_zm, invrs_tau_C14_zm, // In
             stats_pr1);                                                                // Out

    term_pr1(nzm, ngrdcol, C4_zeros, C14, xbp2, wp2, invrs_tau_C4_zm, invrs_tau_C14_zm, // In
             stats_pr2);                                                                // Out

    if (f2k::var_on_stats_list(stats, name_ta)) {
      Kokkos::parallel_for(
          "xp2_xpyp_uv_rhs::array6",
          ((stats_zm_begin.extent(0) - (1)) / (1) + 1) *
              ((stats_zm_begin.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_zm_begin.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_zm_begin.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_zm_begin.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_zm_begin.extent(1) - (1)) / (1) + 1);
            stats_zm_begin(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_uv_rhs::array7",
          ((stats_zm_mod.extent(0) - (1)) / (1) + 1) * ((stats_zm_mod.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_zm_mod.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_zm_mod.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_zm_mod.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_zm_mod.extent(1) - (1)) / (1) + 1);
            stats_zm_mod(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_uv_rhs::loop8", std::max(0, nzm - 2) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 2 + remaining % (nzm - 2);
            remaining /= (nzm - 2);
            stats_zm_begin(i - 1, k - 1) = -rhs_ta(i - 1, k - 1);
            stats_zm_mod(i - 1, k - 1) =
                (one - gamma_over_implicit_ts) * (-lhs_ta(2 - gr_grid_dir_indx - 1, i - 1, k - 1) *
                                                      xap2(i - 1, k + gr_grid_dir_indx - 1) -
                                                  lhs_ta(2 - 1, i - 1, k - 1) * xap2(i - 1, k - 1) -
                                                  lhs_ta(2 + gr_grid_dir_indx - 1, i - 1, k - 1) *
                                                      xap2(i - 1, k - gr_grid_dir_indx - 1));
          });

      // ta has implicit and explicit pieces: begin with the rhs part and add
      // the over-implicit correction with update_budget.
      f2k::stats_begin_budget(name_ta, stats_zm_begin, // In
                              stats);                  // InOut
      f2k::stats_update_budget(name_ta, stats_zm_mod,  // In
                               stats);                 // InOut
    }

    if (f2k::var_on_stats_list(stats, name_pr1)) {
      Kokkos::parallel_for(
          "xp2_xpyp_uv_rhs::array9",
          ((stats_zm_begin.extent(0) - (1)) / (1) + 1) *
              ((stats_zm_begin.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_zm_begin.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_zm_begin.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_zm_begin.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_zm_begin.extent(1) - (1)) / (1) + 1);
            stats_zm_begin(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_uv_rhs::array10",
          ((stats_zm_mod.extent(0) - (1)) / (1) + 1) * ((stats_zm_mod.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_zm_mod.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_zm_mod.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_zm_mod.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_zm_mod.extent(1) - (1)) / (1) + 1);
            stats_zm_mod(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_uv_rhs::loop11", std::max(0, nzm - 2) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 2 + remaining % (nzm - 2);
            remaining /= (nzm - 2);
            stats_zm_begin(i - 1, k - 1) = -stats_pr1(i - 1, k - 1);
            stats_zm_mod(i - 1, k - 1) =
                (one - gamma_over_implicit_ts) * (-lhs_dp1_C4(i - 1, k - 1) * xap2(i - 1, k - 1));
          });

      f2k::stats_begin_budget(name_pr1, stats_zm_begin, // In
                              stats);                   // InOut
      f2k::stats_update_budget(name_pr1, stats_zm_mod,  // In
                               stats);                  // InOut
    }

    if (f2k::var_on_stats_list(stats, name_dp1)) {
      Kokkos::parallel_for(
          "xp2_xpyp_uv_rhs::array12",
          ((stats_zm_begin.extent(0) - (1)) / (1) + 1) *
              ((stats_zm_begin.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_zm_begin.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_zm_begin.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_zm_begin.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_zm_begin.extent(1) - (1)) / (1) + 1);
            stats_zm_begin(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_uv_rhs::array13",
          ((stats_zm_mod.extent(0) - (1)) / (1) + 1) * ((stats_zm_mod.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_zm_mod.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_zm_mod.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_zm_mod.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_zm_mod.extent(1) - (1)) / (1) + 1);
            stats_zm_mod(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_uv_rhs::loop14", std::max(0, nzm - 2) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 2 + remaining % (nzm - 2);
            remaining /= (nzm - 2);
            stats_zm_begin(i - 1, k - 1) = -stats_pr2(i - 1, k - 1);
            stats_zm_mod(i - 1, k - 1) =
                (one - gamma_over_implicit_ts) * (-lhs_dp1_C14(i - 1, k - 1) * xap2(i - 1, k - 1));
          });

      // dp1 has implicit and explicit pieces: begin with the rhs part and add
      // the over-implicit correction with update_budget.
      f2k::stats_begin_budget(name_dp1, stats_zm_begin, // In
                              stats);                   // InOut
      f2k::stats_update_budget(name_dp1, stats_zm_mod,  // In
                               stats);                  // InOut
    }

    Kokkos::parallel_for(
        "xp2_xpyp_uv_rhs::array15",
        ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
        KOKKOS_LAMBDA(int cell) {
          int remaining     = cell;
          const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
          remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
          const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
          remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
          stats_tmp(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xp2_xpyp_uv_rhs::loop16", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (nzm - 2);
          remaining /= (nzm - 2);
          stats_tmp(i - 1, k - 1) = rhs_pr2(i - 1, k - 1);
        });

    // pr2 is explicit in this branch.
    f2k::stats_update(name_pr2, stats_tmp, // In
                      stats);              // InOut

    Kokkos::parallel_for(
        "xp2_xpyp_uv_rhs::array17",
        ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
        KOKKOS_LAMBDA(int cell) {
          int remaining     = cell;
          const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
          remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
          const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
          remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
          stats_tmp(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xp2_xpyp_uv_rhs::loop18", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (nzm - 2);
          remaining /= (nzm - 2);
          stats_tmp(i - 1, k - 1) = (one - C_uu_shr(i - 1)) * rhs_term_tp(i - 1, k - 1);
        });

    // tp is explicit in this branch.
    f2k::stats_update(name_tp, stats_tmp, // In
                      stats);             // InOut

    Kokkos::parallel_for(
        "xp2_xpyp_uv_rhs::array19",
        ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
        KOKKOS_LAMBDA(int cell) {
          int remaining     = cell;
          const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
          remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
          const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
          remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
          stats_tmp(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xp2_xpyp_uv_rhs::loop20", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (nzm - 2);
          remaining /= (nzm - 2);
          stats_tmp(i - 1, k - 1) = real_t(0.5) * lhs_splat_wp2(i - 1, k - 1) * wp2(i - 1, k - 1);
        });

    // Vertical compression of eddies (splat) contribution.
    f2k::stats_update(name_splat, stats_tmp, // In
                      stats);                // InOut
  }

  // Boundary Conditions
  // These are set so that the surface_varnce value of u'^2 or v'^2 can be
  // used at the lowest boundary and the values of those variables can be
  // set to their respective threshold minimum values at the top boundary.
  // Fixed-point boundary conditions are used for the variances.

  Kokkos::parallel_for(
      "xp2_xpyp_uv_rhs::loop21", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        rhs(i - 1, gr_k_lb_zm - 1) = xap2(i - 1, gr_k_lb_zm - 1);
        // The value of u'^2 or v'^2 at the upper boundary will be set to the
        // threshold minimum value of w_tol_sqd.
        rhs(i - 1, gr_k_ub_zm - 1) = w_tol_sqd;
      });

  return;
}

// dt: Model timestep                              [s]
// wpxap: w'x_a' (momentum levels)     [m/s{x_a units}]
// wpxbp: w'x_b' (momentum levels)     [m/s{x_b units}]
// xam: x_am (thermodynamic levels)     [{x_a units}]
// xbm: x_bm (thermodynamic levels)     [{x_b units}]
// xapxbp: x_a'x_b' (m-levs)          [{x_a un}{x_b un}]
// xpyp_forcing: <x'y'> forcing (m-levs)      [{x un}{x un}/s]
// Cn: Coefficient C_n                           [-]
// invrs_tau_zm: Time-scale tau on momentum levels       [1/s]
// threshold: Smallest allowable mag. value for x_a'x_b'  [{x_am units} *{x_bm units}]
// lhs_ta: LHS turbulent advection (ta) term
// rhs_ta: RHS turbulent advection (ta) term
// rhs: Explicit contributions to x variance/covariance terms
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10>
void xp2_xpyp_rhs(int nzm, int nzt, int ngrdcol, const grid& gr, int solve_type, real_t dt, // In
                  Array0 wpxap, Array1 wpxbp,                                               // In
                  Array2 xam, Array3 xbm, Array4 xapxbp, Array5 xpyp_forcing,               // In
                  Array6 Cn, Array7 invrs_tau_zm, real_t threshold,                         // In
                  Array8 lhs_ta, Array9 rhs_ta,                                             // In
                  f2k::stats_type& stats,                                                   // InOut
                  Array10 rhs)                                                              // Out
{
  // Description:
  //   Explicit contributions to r_t'^2, th_l'^2, r_t'th_l', sclr'r_t',
  //   sclr'th_l', or sclr'^2.
  //
  //   Note:
  //     For LHS turbulent advection term:
  //        An "over-implicit" weighted time step is applied to this term.
  //        The weight of the implicit portion of this term is controlled by
  //        the factor gamma_over_implicit_ts (abbreviated "gamma" in the
  //        expression below).  A factor is added to the right-hand side of
  //        the equation in order to balance a weight that is not equal to 1,
  //        such that:
  //             -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
  //        where X is the variable that is being solved for in a predictive
  //        equation (x'^2 or x'y' in this case), y(t) is the linearized
  //        portion of the term that gets treated implicitly, and RHS is the
  //        portion of the term that is always treated explicitly.  A weight
  //        of greater than 1 can be applied to make the term more
  //        numerically stable.
  //
  //   --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
  //   Significant changes to this routine may adversely affect computational speed
  //       - Gunther Huebler, Aug. 2018, clubb:ticket:834
  //----------------------------------------------------------------------------------
  //------------------- Input Variables -------------------
  // For "over-implicit" weighted time step.
  // This vector holds output from the LHS (implicit) portion of a term at a
  // given vertical level.  This output is weighted and applied to the RHS.
  // This is used if the implicit portion of the term is "over-implicit", which
  // means that the LHS contribution is given extra weight (>1) in order to
  // increase numerical stability.  A weighted factor must then be applied to
  // the RHS in order to balance the weight.
  //------------------- Output Variables -------------------
  //---------------- Local Variables -------------------
  // Fractional amount of xp2 that microphysics (mc) is allowed to deplete,
  // where 0 <= mc_xp2_deplete_frac <= 1.  When mc_xp2_deplete_frac has a value
  // of 0, microphysics is not allowed to deplete xp2 at all, and when
  // mc_xp2_deplete_frac has a value of 1, microphysics is allowed to deplete
  // xp2 down to the threshold value of x_tol^2.  This is used only when the
  // l_clip_large_neg_mc flag is enabled.
  // Array indices
  real_t xp2_mc_limiter; // Largest allowable (negative) mc effect
  bool l_clip_large_neg_mc             = false;
  constexpr real_t mc_xp2_deplete_frac = real_t(0.75);
  int i;
  int k;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_term_tp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::rhs_term_tp"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rhs_term_dp1(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::rhs_term_dp1"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> lhs_term_dp1(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::lhs_term_dp1"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_tp1(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::stats_tp1"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_tp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::stats_tp2"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wpxp_zeros(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::wpxp_zeros"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> xm_zeros(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::xm_zeros"), ngrdcol, nzt);
  std::string name_ta;
  std::string name_tp;
  std::string name_tp1;
  std::string name_tp2;
  std::string name_dp1;
  std::string name_f;
  bool l_ta_stats;
  bool l_tp_stats;
  bool l_tp1_stats;
  bool l_tp2_stats;
  bool l_dp1_stats;
  bool l_f_stats;
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_ta_begin(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::stats_ta_begin"), ngrdcol,
      nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_dp1_begin(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::stats_dp1_begin"), ngrdcol,
      nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_ta_mod(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::stats_ta_mod"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_dp1_mod(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::stats_dp1_mod"), ngrdcol, nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> stats_tmp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "xp2_xpyp_rhs::stats_tmp"), ngrdcol, nzm);
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_lb_zm       = gr.k_lb_zm;
  const auto gr_k_ub_zm       = gr.k_ub_zm;

  //------------------------------ Begin Code ---------------------------------

  switch (solve_type) {
  // rtp2/thlp2 case (1 turbulent production term)
  case xp2_xpyp_rtp2: {
    name_ta  = "rtp2_ta";
    name_tp  = "rtp2_tp";
    name_tp1 = "";
    name_tp2 = "";
    name_dp1 = "rtp2_dp1";
    name_f   = "rtp2_forcing";
    break;
  }
  case xp2_xpyp_thlp2: {
    name_ta  = "thlp2_ta";
    name_tp  = "thlp2_tp";
    name_tp1 = "";
    name_tp2 = "";
    name_dp1 = "thlp2_dp1";
    name_f   = "thlp2_forcing";
    // rtpthlp case (2 turbulent production terms)
    break;
  }
  case xp2_xpyp_rtpthlp: {
    name_ta  = "rtpthlp_ta";
    name_tp  = "";
    name_tp1 = "rtpthlp_tp1";
    name_tp2 = "rtpthlp_tp2";
    name_dp1 = "rtpthlp_dp1";
    name_f   = "rtpthlp_forcing";
    break;
  }
  default: {
    name_ta  = "";
    name_tp  = "";
    name_tp1 = "";
    name_tp2 = "";
    name_dp1 = "";
    name_f   = "";
    break;
  }
  }

  // x'y' term pr2 is completely explicit; call stat_update_var_pt.
  // x'y' term tp is completely explicit; call stat_update_var_pt.
  // RHS turbulent production (tp) term.
  term_tp_rhs(nzm, nzt, ngrdcol, xam, xbm, wpxbp, wpxap, gr.invrs_dzm, // In
              rhs_term_tp);                                            // Out

  // RHS dissipation term 1 (dp1)
  term_dp1_rhs(nzm, ngrdcol, Cn, invrs_tau_zm, threshold, // In
               rhs_term_dp1);                             // Out

  // RHS contribution from "over-implicit" weighted time step
  // for LHS dissipation term 1 (dp1).
  term_dp1_lhs(nzm, ngrdcol, gr, Cn, invrs_tau_zm, // In
               lhs_term_dp1);                      // Out

  // Finish RHS calc with vectorizable loop, functions are in source file and should
  // be inlined with an -O2 or above compiler optimization flag
  Kokkos::parallel_for(
      "xp2_xpyp_rhs::loop1", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        rhs(i - 1, k - 1) =
            rhs_ta(i - 1, k - 1) +
            (one - gamma_over_implicit_ts) * (-lhs_ta(2 - gr_grid_dir_indx - 1, i - 1, k - 1) *
                                                  xapxbp(i - 1, k + gr_grid_dir_indx - 1) -
                                              lhs_ta(2 - 1, i - 1, k - 1) * xapxbp(i - 1, k - 1) -
                                              lhs_ta(2 + gr_grid_dir_indx - 1, i - 1, k - 1) *
                                                  xapxbp(i - 1, k - gr_grid_dir_indx - 1));

        // RHS turbulent production (tp) term.
        rhs(i - 1, k - 1) = rhs(i - 1, k - 1) + rhs_term_tp(i - 1, k - 1);

        // RHS dissipation term 1 (dp1)
        rhs(i - 1, k - 1) = rhs(i - 1, k - 1) + rhs_term_dp1(i - 1, k - 1);

        // RHS contribution from "over-implicit" weighted time step
        // for LHS dissipation term 1 (dp1).
        rhs(i - 1, k - 1) =
            rhs(i - 1, k - 1) +
            (one - gamma_over_implicit_ts) * (-lhs_term_dp1(i - 1, k - 1) * xapxbp(i - 1, k - 1));
      });

  // RHS <x'y'> forcing.
  // Note: <x'y'> forcing includes the effects of microphysics on <x'y'>.
  if (l_clip_large_neg_mc && (solve_type == xp2_xpyp_rtp2 || solve_type == xp2_xpyp_thlp2)) {
    Kokkos::parallel_for(
        "xp2_xpyp_rhs::loop2", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (nzm - 2);
          remaining /= (nzm - 2);
          real_t xp2_mc_limiter;

          // Limit the variance-depleting effects of excessively large
          // microphysics terms on rtp2 and thlp2 in order to reduce oscillations.
          if (rhs_term_tp(i - 1, k - 1) >= zero) {

            // Microphysics is allowed to deplete turbulent production and a
            // fraction of the variance, determined by mc_xp2_deplete_frac, down
            // to the threshold.
            xp2_mc_limiter = -mc_xp2_deplete_frac * (xapxbp(i - 1, k - 1) - threshold) / dt -
                             rhs_term_tp(i - 1, k - 1);
          } else {

            // Microphysics is allowed to deplete a fraction of the variance,
            // determined by mc_xp2_deplete_frac, down to the threshold.
            xp2_mc_limiter = -mc_xp2_deplete_frac * (xapxbp(i - 1, k - 1) - threshold) / dt;
          }

          // Note: the value of xp2_mc_limiter is always negative.
          rhs(i - 1, k - 1) = rhs(i - 1, k - 1) + ((xpyp_forcing(i - 1, k - 1)) > (xp2_mc_limiter)
                                                       ? (xpyp_forcing(i - 1, k - 1))
                                                       : (xp2_mc_limiter));
        });
  } else {
    Kokkos::parallel_for(
        "xp2_xpyp_rhs::loop3", std::max(0, nzm - 2) * std::max(0, ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int i   = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 2 + remaining % (nzm - 2);
          remaining /= (nzm - 2);
          rhs(i - 1, k - 1) = rhs(i - 1, k - 1) + xpyp_forcing(i - 1, k - 1);
        });
  }

  Kokkos::parallel_for(
      "xp2_xpyp_rhs::loop4", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        rhs(i - 1, k - 1) = rhs(i - 1, k - 1) + one / dt * xapxbp(i - 1, k - 1);
      });

  if (stats.flags().l_sample) {
    // Statistics: explicit contributions for up2 or vp2.
    // Statistics: explicit contributions for rtp2, thlp2, or rtpthlp.
    Kokkos::parallel_for(
        "xp2_xpyp_rhs::array5",
        ((xm_zeros.extent(0) - (1)) / (1) + 1) * ((xm_zeros.extent(1) - (1)) / (1) + 1),
        KOKKOS_LAMBDA(int cell) {
          int remaining     = cell;
          const int slice_0 = remaining % ((xm_zeros.extent(0) - (1)) / (1) + 1);
          remaining /= ((xm_zeros.extent(0) - (1)) / (1) + 1);
          const int slice_1 = remaining % ((xm_zeros.extent(1) - (1)) / (1) + 1);
          remaining /= ((xm_zeros.extent(1) - (1)) / (1) + 1);
          xm_zeros(slice_0, slice_1) = zero;
        });
    Kokkos::parallel_for(
        "xp2_xpyp_rhs::array6",
        ((wpxp_zeros.extent(0) - (1)) / (1) + 1) * ((wpxp_zeros.extent(1) - (1)) / (1) + 1),
        KOKKOS_LAMBDA(int cell) {
          int remaining     = cell;
          const int slice_0 = remaining % ((wpxp_zeros.extent(0) - (1)) / (1) + 1);
          remaining /= ((wpxp_zeros.extent(0) - (1)) / (1) + 1);
          const int slice_1 = remaining % ((wpxp_zeros.extent(1) - (1)) / (1) + 1);
          remaining /= ((wpxp_zeros.extent(1) - (1)) / (1) + 1);
          wpxp_zeros(slice_0, slice_1) = zero;
        });

    // Note:  To find the contribution of x'y' term tp1, substitute 0 for all
    //        the xam inputs and the wpxbp input to function term_tp.
    term_tp_rhs(nzm, nzt, ngrdcol, xm_zeros, xbm, wpxp_zeros, wpxap, gr.invrs_dzm, // In
                stats_tp1);                                                        // Out

    // Note:  To find the contribution of x'y' term tp2, substitute 0 for all
    //        the xbm inputs and the wpxap input to function term_tp.
    term_tp_rhs(nzm, nzt, ngrdcol, xam, xm_zeros, wpxbp, wpxp_zeros, gr.invrs_dzm, // In
                stats_tp2);                                                        // Out

    l_ta_stats  = name_ta.size() > 0 && f2k::var_on_stats_list(stats, name_ta);
    l_tp_stats  = name_tp.size() > 0 && f2k::var_on_stats_list(stats, name_tp);
    l_tp1_stats = name_tp1.size() > 0 && f2k::var_on_stats_list(stats, name_tp1);
    l_tp2_stats = name_tp2.size() > 0 && f2k::var_on_stats_list(stats, name_tp2);
    l_dp1_stats = name_dp1.size() > 0 && f2k::var_on_stats_list(stats, name_dp1);
    l_f_stats   = name_f.size() > 0 && f2k::var_on_stats_list(stats, name_f);

    if (l_ta_stats || l_dp1_stats) {
      Kokkos::parallel_for(
          "xp2_xpyp_rhs::array7",
          ((stats_ta_begin.extent(0) - (1)) / (1) + 1) *
              ((stats_ta_begin.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_ta_begin.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_ta_begin.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_ta_begin.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_ta_begin.extent(1) - (1)) / (1) + 1);
            stats_ta_begin(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_rhs::array8",
          ((stats_dp1_begin.extent(0) - (1)) / (1) + 1) *
              ((stats_dp1_begin.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_dp1_begin.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_dp1_begin.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_dp1_begin.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_dp1_begin.extent(1) - (1)) / (1) + 1);
            stats_dp1_begin(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_rhs::array9",
          ((stats_ta_mod.extent(0) - (1)) / (1) + 1) * ((stats_ta_mod.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_ta_mod.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_ta_mod.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_ta_mod.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_ta_mod.extent(1) - (1)) / (1) + 1);
            stats_ta_mod(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_rhs::array10",
          ((stats_dp1_mod.extent(0) - (1)) / (1) + 1) * ((stats_dp1_mod.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_dp1_mod.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_dp1_mod.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_dp1_mod.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_dp1_mod.extent(1) - (1)) / (1) + 1);
            stats_dp1_mod(slice_0, slice_1) = zero;
          });

      Kokkos::parallel_for(
          "xp2_xpyp_rhs::loop11", std::max(0, nzm - 2) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 2 + remaining % (nzm - 2);
            remaining /= (nzm - 2);
            if (l_ta_stats) {
              stats_ta_begin(i - 1, k - 1) = -rhs_ta(i - 1, k - 1);
              stats_ta_mod(i - 1, k - 1)   = (one - gamma_over_implicit_ts) *
                                             (-lhs_ta(2 - gr_grid_dir_indx - 1, i - 1, k - 1) *
                                                  xapxbp(i - 1, k + gr_grid_dir_indx - 1) -
                                              lhs_ta(2 - 1, i - 1, k - 1) * xapxbp(i - 1, k - 1) -
                                              lhs_ta(2 + gr_grid_dir_indx - 1, i - 1, k - 1) *
                                                  xapxbp(i - 1, k - gr_grid_dir_indx - 1));
            }
            if (l_dp1_stats) {
              stats_dp1_begin(i - 1, k - 1) = -rhs_term_dp1(i - 1, k - 1);
              stats_dp1_mod(i - 1, k - 1)   = (one - gamma_over_implicit_ts) *
                                              (-lhs_term_dp1(i - 1, k - 1) * xapxbp(i - 1, k - 1));
            }
          });
    }

    if (l_ta_stats) {
      // <x'y'> term ta has both implicit and explicit components; call
      // stat_begin_update_pt.  Since stat_begin_update_pt automatically
      // subtracts the value sent in, reverse the sign on term_ta_ADG1_rhs.
      // x'y' term ta has both implicit and explicit components; call
      // stat_begin_update_pt.  Since stat_begin_update_pt automatically
      // subtracts the value sent in, reverse the sign on term_ta_explicit_rhs.
      // ta has implicit and explicit pieces: begin with the rhs part and add
      // the over-implicit correction with update_budget.
      f2k::stats_begin_budget(name_ta, stats_ta_begin, // In
                              stats);                  // InOut
      f2k::stats_update_budget(name_ta, stats_ta_mod,  // In
                               stats);                 // InOut
    }
    if (l_dp1_stats) {
      // x'y' term dp1 has both implicit and explicit components; call
      // stat_begin_update_pt.  Since stat_begin_update_pt automatically
      // subtracts the value sent in, reverse the sign on term_dp1_rhs.
      // Note:  An "over-implicit" weighted time step is applied to this term.
      //        A weighting factor of greater than 1 may be used to make the
      //        term more numerically stable (see note above for RHS turbulent
      //        advection (ta) term).
      // dp1 has implicit and explicit pieces: begin with the rhs part and add
      // the over-implicit correction with update_budget.
      f2k::stats_begin_budget(name_dp1, stats_dp1_begin, // In
                              stats);                    // InOut
      f2k::stats_update_budget(name_dp1, stats_dp1_mod,  // In
                               stats);                   // InOut
    }
    if (l_tp_stats) {
      Kokkos::parallel_for(
          "xp2_xpyp_rhs::array12",
          ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
            stats_tmp(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_rhs::loop13", std::max(0, nzm - 2) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 2 + remaining % (nzm - 2);
            remaining /= (nzm - 2);
            stats_tmp(i - 1, k - 1) = rhs_term_tp(i - 1, k - 1);
          });

      // tp is explicit in this branch.
      f2k::stats_update(name_tp, stats_tmp, // In
                        stats);             // InOut
    }
    if (l_tp1_stats) {
      Kokkos::parallel_for(
          "xp2_xpyp_rhs::array14",
          ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
            stats_tmp(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_rhs::loop15", std::max(0, nzm - 2) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 2 + remaining % (nzm - 2);
            remaining /= (nzm - 2);
            stats_tmp(i - 1, k - 1) = stats_tp1(i - 1, k - 1);
          });

      // x'y' term tp1 is completely explicit; call stat_update_var_pt.
      // tp1 is explicit in this branch.
      f2k::stats_update(name_tp1, stats_tmp, // In
                        stats);              // InOut
    }
    if (l_tp2_stats) {
      Kokkos::parallel_for(
          "xp2_xpyp_rhs::array16",
          ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
            stats_tmp(slice_0, slice_1) = zero;
          });
      Kokkos::parallel_for(
          "xp2_xpyp_rhs::loop17", std::max(0, nzm - 2) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 2 + remaining % (nzm - 2);
            remaining /= (nzm - 2);
            stats_tmp(i - 1, k - 1) = stats_tp2(i - 1, k - 1);
          });

      // x'y' term tp2 is completely explicit; call stat_update_var_pt.
      // tp2 is explicit in this branch.
      f2k::stats_update(name_tp2, stats_tmp, // In
                        stats);              // InOut
    }
    if (l_f_stats) {
      Kokkos::parallel_for(
          "xp2_xpyp_rhs::array18",
          ((stats_tmp.extent(0) - (1)) / (1) + 1) * ((stats_tmp.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((stats_tmp.extent(0) - (1)) / (1) + 1);
            remaining /= ((stats_tmp.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((stats_tmp.extent(1) - (1)) / (1) + 1);
            remaining /= ((stats_tmp.extent(1) - (1)) / (1) + 1);
            stats_tmp(slice_0, slice_1) = zero;
          });
      if (l_clip_large_neg_mc && (solve_type == xp2_xpyp_rtp2 || solve_type == xp2_xpyp_thlp2)) {
        Kokkos::parallel_for(
            "xp2_xpyp_rhs::loop19", std::max(0, nzm - 2) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 2 + remaining % (nzm - 2);
              remaining /= (nzm - 2);
              stats_tmp(i - 1, k - 1) =
                  ((xpyp_forcing(i - 1, k - 1)) > (xp2_mc_limiter) ? (xpyp_forcing(i - 1, k - 1))
                                                                   : (xp2_mc_limiter));
            });
      } else {
        Kokkos::parallel_for(
            "xp2_xpyp_rhs::loop20", std::max(0, nzm - 2) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 2 + remaining % (nzm - 2);
              remaining /= (nzm - 2);
              stats_tmp(i - 1, k - 1) = xpyp_forcing(i - 1, k - 1);
            });
      }

      // x'y' forcing term is completely explicit; call stat_update_var_pt.
      f2k::stats_update(name_f, stats_tmp, // In
                        stats);            // InOut
    }
  }

  // Boundary Conditions
  // These are set so that the surface_varnce value of rtp2, thlp2, or rtpthlp
  // (or sclrp2, sclrprtp, or sclrpthlp) can be used at the lowest boundary and the
  // values of those variables can be set to their respective threshold minimum
  // values (which is 0 in the case of the covariances) at the top boundary.
  // Fixed-point boundary conditions are used for both the variances and the
  // covariances.

  // The value of the field at the upper boundary will be set to it's threshold
  // minimum value, as contained in the variable 'threshold'.
  Kokkos::parallel_for(
      "xp2_xpyp_rhs::loop21", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        rhs(i - 1, gr_k_lb_zm - 1) = xapxbp(i - 1, gr_k_lb_zm - 1);
        rhs(i - 1, gr_k_ub_zm - 1) = threshold;
      });

  return;
}

// nzm: Number of momentum vertical levels
// nzt: Number of thermodynamic vertical levels
// ngrdcol: Number of grid columns
// sclr_dim: Number of passive scalars
// wprtp: <w'r_t'> (momentum levels)            [(m/s)(kg/kg)]
// wprtp2: <w'r_t'^2> (thermodynamic levels)     [m/s (kg/kg)^2]
// wpthlp: <w'th_l'> (momentum levels)           [(m K)/s]
// wpthlp2: <w'th_l'^2> (thermodynamic levels)    [m/s K^2]
// wprtpthlp: <w'r_t'th_l'> (thermodynamic levels)  [m/s (kg/kg) K]
// rtp2: <r_t'^2>                              [(kg/kg)^2]
// thlp2: <th_l'^2>                             [K^2]
// rtpthlp: <r_t'th_l'>                           [(kg K)/kg]
// upwp: <u'w'> (momentum levels)              [m^2/s^2]
// vpwp: <v'w'> (momentum levels)              [m^2/s^2]
// up2: <u'^2> (momentum levels)              [m^s/s^2]
// vp2: <v'^2> (momentum levels)              [m^2/s^2]
// wp2: <w'^2> (momentum levels)              [m^2/s^2]
// wp2_zt: <w'^2> interpolated to thermo. levels [m^2/s^2]
// wpsclrp: <w'sclr'> (momentum levels)        [m/s{sclr units}]
// wpsclrp2: <w'sclr'^2> (thermodynamic levels) [m/s{sclr units}^2]
// wpsclrprtp: <w'sclr'r_t'> (thermo. levels)     [m/s{sclr units)kg/kg]
// wpsclrpthlp: <w'sclr'th_l'> (thermo. levels)    [m/s{sclr units}K]
// rho_ds_zt: Dry, static density on thermo. levels [kg/m^3]
// invrs_rho_ds_zm: Inv. dry, static density @ mom. levs. [m^3/kg]
// rho_ds_zm: Dry, static density on momentum levs. [kg/m^3]
// wp3_on_wp2: Smoothed version of <w'^3>/<w'^2> zm  [m/s]
// wp3_on_wp2_zt: Smoothed version of <w'^3>/<w'^2> zt  [m/s]
// sigma_sqd_w: sigma_sqd_w (momentum levels)         [-]
// pdf_implicit_coefs_terms: Implicit coefs / explicit terms [units vary]
// beta: CLUBB tunable parameter beta
// iiPDF_type: Selected option for the two-component normal (double Gaussian) PDF type to use for
// the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate, two-component
// PDF.
// l_upwind_xpyp_ta: This flag determines whether we want to use an upwind differencing
// approximation rather than a centered differencing for turbulent or mean advection terms. It
// affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, & sclrpthlp.
// l_godunov_upwind_xpyp_ta: This flag determines whether we want to use a Godunov-like upwind
// approximation rather than a centered differencing for  turbulent or mean advection terms. It
// affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, & sclrpthlp.
// lhs_ta_wprtp2: For <w'rt'^2>
// lhs_ta_wpthlp2: For <w'thl'^2>
// lhs_ta_wprtpthlp: For <w'rt'thl'>
// lhs_ta_wpup2: For <w'u'^2>
// lhs_ta_wpvp2: For <w'v'^2>
// lhs_ta_wpsclrp2: For <w'sclr'^2>
// lhs_ta_wprtpsclrp: For <w'rt'sclr'>
// lhs_ta_wpthlpsclrp: For <w'thl'sclr'>
// rhs_ta_wprtp2: For <w'rt'^2>
// rhs_ta_wpthlp2: For <w'thl'^2>
// rhs_ta_wprtpthlp: For <w'rt'thl'>
// rhs_ta_wpup2: For <w'u'^2>
// rhs_ta_wpvp2: For <w'v'^2>
// rhs_ta_wpsclrp2: For <w'sclr'^2>
// rhs_ta_wprtpsclrp: For <w'rt'sclr'>
// rhs_ta_wpthlpsclrp: For <w'thl'sclr'>
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43>
void calc_xp2_xpyp_ta_terms(
    int nzm, int nzt, int ngrdcol, int sclr_dim, const grid& gr,                  // In
    Array0 wprtp, Array1 wprtp2, Array2 wpthlp, Array3 wpthlp2, Array4 wprtpthlp, // In
    Array5 rtp2, Array6 thlp2, Array7 rtpthlp, Array8 upwp, Array9 vpwp,          // In
    Array10 up2, Array11 vp2, Array12 wp2,                                        // In
    Array13 wp2_zt, Array14 wpsclrp, Array15 wpsclrp2, Array16 wpsclrprtp,        // In
    Array17 wpsclrpthlp,                                                          // In
    Array18 sclrp2, Array19 sclrprtp, Array20 sclrpthlp,                          // In
    Array21 rho_ds_zt, Array22 invrs_rho_ds_zm, Array23 rho_ds_zm,                // In
    Array24 wp3_on_wp2, Array25 wp3_on_wp2_zt, Array26 sigma_sqd_w,               // In
    const implicit_coefs_terms& pdf_implicit_coefs_terms, bool l_scalar_calc,     // In
    Array27 beta, int iiPDF_type, bool l_upwind_xpyp_ta,                          // In
    bool l_godunov_upwind_xpyp_ta,                                                // In
    f2k::stats_type& stats,                                                       // InOut
    Array28 lhs_ta_wprtp2, Array29 lhs_ta_wpthlp2, Array30 lhs_ta_wprtpthlp,      // Out
    Array31 lhs_ta_wpup2, Array32 lhs_ta_wpvp2, Array33 lhs_ta_wpsclrp2,          // Out
    Array34 lhs_ta_wprtpsclrp, Array35 lhs_ta_wpthlpsclrp,                        // Out
    Array36 rhs_ta_wprtp2, Array37 rhs_ta_wpthlp2, Array38 rhs_ta_wprtpthlp,      // Out
    Array39 rhs_ta_wpup2, Array40 rhs_ta_wpvp2, Array41 rhs_ta_wpsclrp2,          // Out
    Array42 rhs_ta_wprtpsclrp, Array43 rhs_ta_wpthlpsclrp)                        // Out
{
  // Description:
  //   This procedure calculates all the turbulent advection terms needed by the
  //   various LHS and RHS matrices. In general, first the implicit coeficients are
  //   calculated and used to calculate the LHS turbulent advection terms, then the
  //   explicit terms are calculated and used to calculate the RHS turbulent advection
  //   terms.
  //
  //-------------------------------------------------------------------------------------------
  //------------------- Input Variables -------------------
  // Passive scalar output
  //------------------- Output Variables -------------------
  // Implicit (LHS) turbulent advection terms
  // Explicit (RHS) turbulent advection terms
  // Implicit (LHS) turbulent advection terms for scalars
  // Explicit (RHS) turbulent advection terms for scalars
  //------------------- Local Variable -------------------
  // Variables for turbulent advection of predictive variances and covariances.
  // <w'rt'^2> = coef_wprtp2_implicit * <rt'^2> + term_wprtp2_explicit
  // <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2> + term_wpthlp2_explicit
  // <w'rt'thl'> = coef_wprtpthlp_implicit*<rt'thl'> + term_wprtpthlp_explicit
  // CLUBB does not produce a PDF for horizontal wind components u and v.
  // However, turbulent advection of the variances of the horizontal wind
  // components, <u'^2> and <v'^2>, is still handled by equations of the form:
  // <w'u'^2> = coef_wpup2_implicit * <u'^2> + term_wpup2_explicit; and
  // <w'v'^2> = coef_wpvp2_implicit * <v'^2> + term_wpvp2_explicit.
  // Sign of turbulent velocity (used for "upwind" turbulent advection)
  // Sign of turbulent velocity (used for "upwind" turbulent advection in Godunov scheme)
  // <w'sclr'^2> = coef_wpsclrp2_implicit * <sclr'^2> + term_wpsclrp2_explicit
  // <w'rt'sclr'> = coef_wprtpsclrp_implicit * <sclr'rt'>
  //                + term_wprtpsclrp_explicit
  // <w'thl'sclr'> = coef_wpthlpsclrp_implicit * <sclr'thl'>
  //                 + term_wpthlpsclrp_explicit
  // Sign of turbulent velocity (used for "upwind" turbulent advection)
  // Sign of turbulent velocity (used for "upwind" turbulent advection in Godunov scheme)
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wprtp2_implicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wprtp2_implicit"),
      ngrdcol, nzt); // Coefficient that is multiplied by <rt'^2>  [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wprtp2_explicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wprtp2_explicit"),
      ngrdcol, nzt); // Term that is on the RHS           [m/s(kg/kg)^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wprtp2_implicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wprtp2_implicit_zm"),
      ngrdcol, nzm); // coef_wprtp2_implicit interp. to m-levs  [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wprtp2_explicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wprtp2_explicit_zm"),
      ngrdcol, nzm); // term_wprtp2_expl interp m-levs [m/s(kg/kg)^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpthlp2_implicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wpthlp2_implicit"),
      ngrdcol, nzt); // Coef. that is multiplied by <thl'^2>      [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpthlp2_explicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wpthlp2_explicit"),
      ngrdcol, nzt); // Term that is on the RHS               [m/s K^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpthlp2_implicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wpthlp2_implicit_zm"),
      ngrdcol, nzm); // coef_wpthlp2_implicit interp. m-levs   [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpthlp2_explicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wpthlp2_explicit_zm"),
      ngrdcol, nzm); // term_wpthlp2_expl interp to m-levs [m/s K^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wprtpthlp_implicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wprtpthlp_implicit"),
      ngrdcol, nzt); // Coef. that is multiplied by <rt'thl'>   [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wprtpthlp_explicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wprtpthlp_explicit"),
      ngrdcol, nzt); // Term that is on the RHS         [m/s(kg/kg)K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wprtpthlp_implicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wprtpthlp_implicit_zm"),
      ngrdcol, nzm); // coef_wprtpthlp_impl interp. m-levs   [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wprtpthlp_explicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wprtpthlp_explicit_zm"),
      ngrdcol, nzm); // term_wprtpthlp_expl intrp zm [m/s(kg/kg)K]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpup2_implicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wpup2_implicit"),
      ngrdcol, nzt); // Coef. that is multiplied by <u'^2>  [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpup2_explicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wpup2_explicit"),
      ngrdcol, nzt); // Term that is on the RHS (<u'^2>)    [m^3/s^3]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpvp2_implicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wpvp2_implicit"),
      ngrdcol, nzt); // Coef. that is multiplied by <v'^2>  [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpvp2_explicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wpvp2_explicit"),
      ngrdcol, nzt); // Term that is on the RHS (<v'^2>)    [m^3/s^3]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpup2_implicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wpup2_implicit_zm"),
      ngrdcol, nzm); // coef_wpup2__impl intrp m-levs   [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpup2_explicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wpup2_explicit_zm"),
      ngrdcol, nzm); // term_wpup2_expl interp. m-levs  [m^3/s^3]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpvp2_implicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wpvp2_implicit_zm"),
      ngrdcol, nzm); // coef_wpvp2_impl intrp m-levs    [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpvp2_explicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wpvp2_explicit_zm"),
      ngrdcol, nzm); // term_wpvp2_expl interp. m-levs  [m^3/s^3]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_rtp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::sgn_t_vel_rtp2"),
      ngrdcol, nzm); // Sign of the turbulent velocity for <rt'^2>      [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_thlp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::sgn_t_vel_thlp2"),
      ngrdcol, nzm); // Sign of the turbulent velocity for <thl'^2>     [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_rtpthlp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::sgn_t_vel_rtpthlp"),
      ngrdcol, nzm); // Sign of the turbulent velocity for <rt'thl'>    [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_up2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::sgn_t_vel_up2"),
      ngrdcol, nzm); // Sign of the turbulent vel. for <u'^2>           [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_vp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::sgn_t_vel_vp2"),
      ngrdcol, nzm); // Sign of the turbulent vel. for <v'^2>           [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_rtp2_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::sgn_t_vel_rtp2_zt"),
      ngrdcol, nzt); // Sign of the turbulent velocity for <rt'^2>      [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_thlp2_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::sgn_t_vel_thlp2_zt"),
      ngrdcol, nzt); // Sign of the turbulent velocity for <thl'^2>     [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_rtpthlp_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::sgn_t_vel_rtpthlp_zt"),
      ngrdcol, nzt); // Sign of the turbulent velocity for <rt'thl'>    [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpsclrp2_implicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wpsclrp2_implicit"),
      ngrdcol, nzt); // Coef. that is multiplied by <sclr'^2>    [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpsclrp2_explicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wpsclrp2_explicit"),
      ngrdcol, nzt); // Term that is on the RHS    [m/s(units vary)^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpsclrp2_implicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wpsclrp2_implicit_zm"),
      ngrdcol, nzm); // coef_wpsclrp2_impl interp zm   [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpsclrp2_explicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wpsclrp2_explicit_zm"),
      ngrdcol, nzm); // term_wpsclrp2_expl interp zm   [un vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wprtpsclrp_implicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wprtpsclrp_implicit"),
      ngrdcol, nzt); // Coef. that is multiplied by <sclr'rt'> [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wprtpsclrp_explicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wprtpsclrp_explicit"),
      ngrdcol, nzt); // Term that is on the RHS [m/s(kg/kg)(un. v.)]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wprtpsclrp_implicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wprtpsclrp_implicit_zm"),
      ngrdcol, nzm); // coef_wprtpsclrp_impl interp zm [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wprtpsclrp_explicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wprtpsclrp_explicit_zm"),
      ngrdcol, nzm); // term_wprtpsclrp_expl interp zm [un vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpthlpsclrp_implicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wpthlpsclrp_implicit"),
      ngrdcol, nzt); // Coef. that is mult. by <sclr'thl'>    [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpthlpsclrp_explicit(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wpthlpsclrp_explicit"),
      ngrdcol, nzt); // Term that is on the RHS  [(m/s)K(un. vary)]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> coef_wpthlpsclrp_implicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::coef_wpthlpsclrp_implicit_zm"),
      ngrdcol, nzm); // coef_wpsclrpthlp_impl intrp zm [m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> term_wpthlpsclrp_explicit_zm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::term_wpthlpsclrp_explicit_zm"),
      ngrdcol, nzm); // term_wpsclrpthlp_expl intrp zm [un vary]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_sclrp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::sgn_t_vel_sclrp2"),
      ngrdcol, nzm); // Sign of the turbulent velocity for <sclr'^2>  [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_sclrprtp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::sgn_t_vel_sclrprtp"),
      ngrdcol, nzm); // Sign of the turbulent velocity for <sclr'rt'> [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_sclrpthlp(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::sgn_t_vel_sclrpthlp"),
      ngrdcol, nzm); // Sign of the turbulent vel. for <sclr'thl'>    [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_sclrp2_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::sgn_t_vel_sclrp2_zt"),
      ngrdcol, nzt); // Sign of the turbulent velocity for <sclr'^2>  [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_sclrprtp_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::sgn_t_vel_sclrprtp_zt"),
      ngrdcol, nzt); // Sign of the turbulent velocity for <sclr'rt'> [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> sgn_t_vel_sclrpthlp_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing,
                         "calc_xp2_xpyp_ta_terms::sgn_t_vel_sclrpthlp_zt"),
      ngrdcol, nzt); // Sign of the turbulent vel. for <sclr'thl'>    [-]
  Kokkos::View<real_t***, Kokkos::LayoutLeft> wpsclrp_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::wpsclrp_zt"),
      ngrdcol, nzt, sclr_dim); // <w'sclr'> interp. to thermo. levels  [m/s {sclrm units}]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> a1_coef(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::a1_coef"), ngrdcol,
      nzm); // a_1 (momentum levels); See eqn. 24 in `Equations for CLUBB' [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wp_coef(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::wp_coef"), ngrdcol,
      nzm);
  Kokkos::View<real_t**, Kokkos::LayoutLeft> a1_coef_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::a1_coef_zt"),
      ngrdcol, nzt); // a_1 interpolated to thermodynamic levels      [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> upwp_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::upwp_zt"), ngrdcol,
      nzt); // <u'w'> interpolated to thermodynamic levels    [m^2/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> vpwp_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::vpwp_zt"), ngrdcol,
      nzt); // <v'w'> interpolated to thermodynamic levels    [m^2/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wprtp_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::wprtp_zt"), ngrdcol,
      nzt); // w'r_t' interpolated to thermodynamic levels   [(kg/kg) m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wpthlp_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::wpthlp_zt"), ngrdcol,
      nzt); // w'th_l' interpolated to thermodyamnic levels  [K m/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wp_coef_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_xp2_xpyp_ta_terms::wp_coef_zt"),
      ngrdcol, nzt);
  int sclr;
  int i;
  int b;
  int k; // Loop index
  const auto pdf_implicit_coefs_terms_coef_wprtp2_implicit =
      pdf_implicit_coefs_terms.coef_wprtp2_implicit;
  const auto pdf_implicit_coefs_terms_coef_wpthlp2_implicit =
      pdf_implicit_coefs_terms.coef_wpthlp2_implicit;
  const auto pdf_implicit_coefs_terms_coef_wprtpthlp_implicit =
      pdf_implicit_coefs_terms.coef_wprtpthlp_implicit;
  const auto pdf_implicit_coefs_terms_term_wprtpthlp_explicit =
      pdf_implicit_coefs_terms.term_wprtpthlp_explicit;
  const auto pdf_implicit_coefs_terms_term_wprtp2_explicit =
      pdf_implicit_coefs_terms.term_wprtp2_explicit;
  const auto pdf_implicit_coefs_terms_term_wpthlp2_explicit =
      pdf_implicit_coefs_terms.term_wpthlp2_explicit;
  const auto pdf_implicit_coefs_terms_coef_wpup2_implicit =
      pdf_implicit_coefs_terms.coef_wpup2_implicit;
  const auto pdf_implicit_coefs_terms_term_wpup2_explicit =
      pdf_implicit_coefs_terms.term_wpup2_explicit;
  const auto pdf_implicit_coefs_terms_coef_wpvp2_implicit =
      pdf_implicit_coefs_terms.coef_wpvp2_implicit;
  const auto pdf_implicit_coefs_terms_term_wpvp2_explicit =
      pdf_implicit_coefs_terms.term_wpvp2_explicit;
  const auto pdf_implicit_coefs_terms_coef_wpsclrp2_implicit =
      pdf_implicit_coefs_terms.coef_wpsclrp2_implicit;
  const auto pdf_implicit_coefs_terms_term_wpsclrp2_explicit =
      pdf_implicit_coefs_terms.term_wpsclrp2_explicit;
  const auto pdf_implicit_coefs_terms_coef_wprtpsclrp_implicit =
      pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit;
  const auto pdf_implicit_coefs_terms_term_wprtpsclrp_explicit =
      pdf_implicit_coefs_terms.term_wprtpsclrp_explicit;
  const auto pdf_implicit_coefs_terms_coef_wpthlpsclrp_implicit =
      pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit;
  const auto pdf_implicit_coefs_terms_term_wpthlpsclrp_explicit =
      pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit;

  //------------------- Begin Code -------------------

  // Define a_1 (located on momentum levels).
  // It is a variable that is a function of sigma_sqd_w (where sigma_sqd_w is
  // located on the momentum levels).  This will be used for the turbulent
  // advection (ta) terms for the ADG1 PDF.  This will also be used for the
  // turbulent advection of <u'^2> and <v'^2>, regardless of which PDF type or
  // turbulent advection option is used.
  Kokkos::parallel_for(
      "calc_xp2_xpyp_ta_terms::loop1", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzm);
        remaining /= (nzm);
        a1_coef(i - 1, k - 1) = one / (one - sigma_sqd_w(i - 1, k - 1));
      });

  // Interpolate a_1 from the momentum levels to the thermodynamic levels.
  Kokkos::deep_copy(Kokkos::subview(a1_coef_zt, Kokkos::ALL, Kokkos::ALL),
                    zm2zt_api(nzm, nzt, ngrdcol, gr,
                              Kokkos::subview(a1_coef, Kokkos::ALL, Kokkos::ALL), zero_threshold));

  Kokkos::parallel_for(
      "calc_xp2_xpyp_ta_terms::loop2", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzm);
        remaining /= (nzm);
        wp_coef(i - 1, k - 1) = (one - one_third * beta(i - 1)) *
                                (a1_coef(i - 1, k - 1) * a1_coef(i - 1, k - 1)) *
                                wp3_on_wp2(i - 1, k - 1) / wp2(i - 1, k - 1);
      });

  Kokkos::parallel_for(
      "calc_xp2_xpyp_ta_terms::loop3", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzt);
        remaining /= (nzt);
        wp_coef_zt(i - 1, k - 1) = (one - one_third * beta(i - 1)) *
                                   (a1_coef_zt(i - 1, k - 1) * a1_coef_zt(i - 1, k - 1)) *
                                   wp3_on_wp2_zt(i - 1, k - 1) / wp2_zt(i - 1, k - 1);
      });

  if (l_explicit_turbulent_adv_xpyp) {

    // The turbulent advection of <x'y'> is handled explicitly, the
    // terms are calculated only for the RHS matrices. The
    // term_wpxpyp_explicit terms are equal to <w'x'y'> as calculated using PDF
    // parameters, which are general for any PDF type. The values of
    // <w'x'y'> are calculated on thermodynamic levels.

    // These coefficients only need to be set if stats output is on

    if (stats.flags().l_sample) {
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::loop4", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzt);
            remaining /= (nzt);
            coef_wprtp2_implicit(i - 1, k - 1)    = zero;
            coef_wpthlp2_implicit(i - 1, k - 1)   = zero;
            coef_wprtpthlp_implicit(i - 1, k - 1) = zero;
          });
    }

    // The turbulent advection terms are handled entirely explicitly. Thus the LHS
    // terms can be set to zero.
    Kokkos::parallel_for(
        "calc_xp2_xpyp_ta_terms::loop5",
        std::max(0, nzm) * std::max(0, ngrdcol) * std::max(0, ndiags3),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int b   = 1 + remaining % (ndiags3);
          remaining /= (ndiags3);
          const int i = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          lhs_ta_wprtp2(b - 1, i - 1, k - 1)    = zero;
          lhs_ta_wpthlp2(b - 1, i - 1, k - 1)   = zero;
          lhs_ta_wprtpthlp(b - 1, i - 1, k - 1) = zero;
        });

    if (l_scalar_calc) {
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::loop6",
          std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol) * std::max(0, ndiags3),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int b   = 1 + remaining % (ndiags3);
            remaining /= (ndiags3);
            const int i = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            const int sclr = 1 + remaining % (sclr_dim);
            remaining /= (sclr_dim);
            lhs_ta_wpsclrp2(b - 1, i - 1, k - 1, sclr - 1)    = zero;
            lhs_ta_wprtpsclrp(b - 1, i - 1, k - 1, sclr - 1)  = zero;
            lhs_ta_wpthlpsclrp(b - 1, i - 1, k - 1, sclr - 1) = zero;
          });
    }

    // The termo-level terms only need to be set if we're not using l_upwind_xpyp_ta,
    // or if stats output is on

    if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::loop7", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzt);
            remaining /= (nzt);
            term_wprtp2_explicit(i - 1, k - 1)    = wprtp2(i - 1, k - 1);
            term_wpthlp2_explicit(i - 1, k - 1)   = wpthlp2(i - 1, k - 1);
            term_wprtpthlp_explicit(i - 1, k - 1) = wprtpthlp(i - 1, k - 1);
          });
    }

    // Interpolate wprtp2 to momentum levels, and calculate the sign of vertical velocity
    if (l_upwind_xpyp_ta) {
      Kokkos::deep_copy(
          Kokkos::subview(term_wprtp2_explicit_zm, Kokkos::ALL, Kokkos::ALL),
          zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wprtp2, Kokkos::ALL, Kokkos::ALL)));

      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::loop8", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            sgn_t_vel_rtp2(i - 1, k - 1) =
                Kokkos::copysign(one, term_wprtp2_explicit_zm(i - 1, k - 1) * rtp2(i - 1, k - 1));
          });
    }

    // Calculate the RHS turbulent advection term for <w'rt'^2>
    f2k::xpyp_term_ta_pdf_rhs(
        nzm, nzt, ngrdcol, gr, term_wprtp2_explicit, rho_ds_zt, rho_ds_zm,          // In
        invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtp2, term_wprtp2_explicit_zm, // In
        rhs_ta_wprtp2);                                                             // Out

    // Interpolate wpthlp2 to momentum levels, and calculate the sign of vertical velocity
    if (l_upwind_xpyp_ta) {
      Kokkos::deep_copy(
          Kokkos::subview(term_wpthlp2_explicit_zm, Kokkos::ALL, Kokkos::ALL),
          zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wpthlp2, Kokkos::ALL, Kokkos::ALL)));

      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::loop9", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            sgn_t_vel_thlp2(i - 1, k - 1) =
                Kokkos::copysign(one, term_wpthlp2_explicit_zm(i - 1, k - 1) * thlp2(i - 1, k - 1));
          });
    }

    // Calculate the RHS turbulent advection term for <w'thl'^2>
    f2k::xpyp_term_ta_pdf_rhs(
        nzm, nzt, ngrdcol, gr, term_wpthlp2_explicit, rho_ds_zt, rho_ds_zm,           // In
        invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_thlp2, term_wpthlp2_explicit_zm, // In
        rhs_ta_wpthlp2);                                                              // Out

    // Interpolate wprtpthlp to momentum levels, and calculate the sign of vertical velocity
    if (l_upwind_xpyp_ta) {
      Kokkos::deep_copy(
          Kokkos::subview(term_wprtpthlp_explicit_zm, Kokkos::ALL, Kokkos::ALL),
          zt2zm_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wprtpthlp, Kokkos::ALL, Kokkos::ALL)));

      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::loop10", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            sgn_t_vel_rtpthlp(i - 1, k - 1) = Kokkos::copysign(
                one, term_wprtpthlp_explicit_zm(i - 1, k - 1) * rtpthlp(i - 1, k - 1));
          });
    }

    // Calculate the RHS turbulent advection term for <w'rt'thl'>
    f2k::xpyp_term_ta_pdf_rhs(
        nzm, nzt, ngrdcol, gr, term_wprtpthlp_explicit, rho_ds_zt, rho_ds_zm,             // In
        invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, term_wprtpthlp_explicit_zm, // In
        rhs_ta_wprtpthlp);                                                                // Out

    if (l_scalar_calc) {

      for (sclr = 1; sclr <= sclr_dim; sclr += 1) {

        // Interpolate wpsclrp2 to momentum levels and calculate the sign of
        // vertical velocityif l_upwind_xpyp_ta, otherwise just use wpsclrp2
        if (l_upwind_xpyp_ta) {
          Kokkos::deep_copy(
              Kokkos::subview(term_wpsclrp2_explicit_zm, Kokkos::ALL, Kokkos::ALL),
              zt2zm_api(nzm, nzt, ngrdcol, gr,
                        Kokkos::subview(wpsclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1))));

          Kokkos::parallel_for(
              "calc_xp2_xpyp_ta_terms::loop11", std::max(0, nzm) * std::max(0, ngrdcol),
              KOKKOS_LAMBDA(const int cell) {
                int remaining = cell;
                const int i   = 1 + remaining % (ngrdcol);
                remaining /= (ngrdcol);
                const int k = 1 + remaining % (nzm);
                remaining /= (nzm);
                sgn_t_vel_sclrp2(i - 1, k - 1) = Kokkos::copysign(
                    one, term_wpsclrp2_explicit_zm(i - 1, k - 1) * sclrp2(i - 1, k - 1, sclr - 1));
              });
        } else {
          Kokkos::parallel_for(
              "calc_xp2_xpyp_ta_terms::loop12", std::max(0, nzt) * std::max(0, ngrdcol),
              KOKKOS_LAMBDA(const int cell) {
                int remaining = cell;
                const int i   = 1 + remaining % (ngrdcol);
                remaining /= (ngrdcol);
                const int k = 1 + remaining % (nzt);
                remaining /= (nzt);
                term_wpsclrp2_explicit(i - 1, k - 1) = wpsclrp2(i - 1, k - 1, sclr - 1);
              });
        }

        // Calculate the RHS turbulent advection term for <w'sclr'^2>
        f2k::xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wpsclrp2_explicit, rho_ds_zt, rho_ds_zm,            // In
            invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrp2, term_wpsclrp2_explicit_zm, // In
            Kokkos::subview(rhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)));        // Out
      }

      // Interpolate wpsclrprtp to momentum levels and calculate the sign of
      // vertical velocityif l_upwind_xpyp_ta, otherwise just use wpsclrprtp
      for (sclr = 1; sclr <= sclr_dim; sclr += 1) {
        if (l_upwind_xpyp_ta) {
          Kokkos::deep_copy(
              Kokkos::subview(term_wprtpsclrp_explicit_zm, Kokkos::ALL, Kokkos::ALL),
              zt2zm_api(nzm, nzt, ngrdcol, gr,
                        Kokkos::subview(wpsclrprtp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))));

          Kokkos::parallel_for(
              "calc_xp2_xpyp_ta_terms::loop13", std::max(0, nzm) * std::max(0, ngrdcol),
              KOKKOS_LAMBDA(const int cell) {
                int remaining = cell;
                const int i   = 1 + remaining % (ngrdcol);
                remaining /= (ngrdcol);
                const int k = 1 + remaining % (nzm);
                remaining /= (nzm);
                sgn_t_vel_sclrprtp(i - 1, k - 1) =
                    Kokkos::copysign(one, term_wprtpsclrp_explicit_zm(i - 1, k - 1) *
                                              sclrprtp(i - 1, k - 1, sclr - 1));
              });
        } else {
          Kokkos::parallel_for(
              "calc_xp2_xpyp_ta_terms::loop14", std::max(0, nzt) * std::max(0, ngrdcol),
              KOKKOS_LAMBDA(const int cell) {
                int remaining = cell;
                const int i   = 1 + remaining % (ngrdcol);
                remaining /= (ngrdcol);
                const int k = 1 + remaining % (nzt);
                remaining /= (nzt);
                term_wprtpsclrp_explicit(i - 1, k - 1) = wpsclrprtp(i - 1, k - 1, sclr - 1);
              });
        }

        // Calculate the RHS turbulent advection term for <w'sclr'rt'>
        f2k::xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wprtpsclrp_explicit, rho_ds_zt, rho_ds_zm,     // In
            invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrprtp,                     // In
            term_wprtpsclrp_explicit_zm,                                               // In
            Kokkos::subview(rhs_ta_wprtpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out
      }

      // Interpolate wpsclrpthlp to momentum levels and calculate the sign of
      // vertical velocityif l_upwind_xpyp_ta, otherwise just use wpsclrpthlp
      for (sclr = 1; sclr <= sclr_dim; sclr += 1) {
        if (l_upwind_xpyp_ta) {
          Kokkos::deep_copy(
              Kokkos::subview(term_wpthlpsclrp_explicit_zm, Kokkos::ALL, Kokkos::ALL),
              zt2zm_api(nzm, nzt, ngrdcol, gr,
                        Kokkos::subview(wpsclrpthlp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))));

          Kokkos::parallel_for(
              "calc_xp2_xpyp_ta_terms::loop15", std::max(0, nzm) * std::max(0, ngrdcol),
              KOKKOS_LAMBDA(const int cell) {
                int remaining = cell;
                const int i   = 1 + remaining % (ngrdcol);
                remaining /= (ngrdcol);
                const int k = 1 + remaining % (nzm);
                remaining /= (nzm);
                sgn_t_vel_sclrpthlp(i - 1, k - 1) =
                    Kokkos::copysign(one, term_wpthlpsclrp_explicit_zm(i - 1, k - 1) *
                                              sclrpthlp(i - 1, k - 1, sclr - 1));
              });
        } else {
          Kokkos::parallel_for(
              "calc_xp2_xpyp_ta_terms::loop16", std::max(0, nzt) * std::max(0, ngrdcol),
              KOKKOS_LAMBDA(const int cell) {
                int remaining = cell;
                const int i   = 1 + remaining % (ngrdcol);
                remaining /= (ngrdcol);
                const int k = 1 + remaining % (nzt);
                remaining /= (nzt);
                term_wpthlpsclrp_explicit(i - 1, k - 1) = wpsclrpthlp(i - 1, k - 1, sclr - 1);
              });
        }

        // Calculate the RHS turbulent advection term for <w'sclr'thl'>
        f2k::xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wpthlpsclrp_explicit, rho_ds_zt, rho_ds_zm,     // In
            invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrpthlp,                     // In
            term_wpthlpsclrp_explicit_zm,                                               // In
            Kokkos::subview(rhs_ta_wpthlpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out
      }
    }

  } else {

    // The turbulent advection of <x'y'> is handled implicitly or
    // semi-implicitly.

    if (iiPDF_type == iiPDF_ADG1) {

      // The ADG1 PDF is used.

      // For ADG1, the sign of the turbulent velocity is the sign of
      // <w'^3> / <w'^2>.  For simplicity, the sign of turbulent
      // velocity is set to wp3_on_wp2 for all terms.

      // The termodynamic grid level coefficients are only needed if l_upwind_xpyp_ta
      // is false, or if stats output is on

      if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop17", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzt);
              remaining /= (nzt);
              coef_wprtp2_implicit(i - 1, k - 1) =
                  one_third * beta(i - 1) * a1_coef_zt(i - 1, k - 1) * wp3_on_wp2_zt(i - 1, k - 1);
              coef_wpthlp2_implicit(i - 1, k - 1)   = coef_wprtp2_implicit(i - 1, k - 1);
              coef_wprtpthlp_implicit(i - 1, k - 1) = coef_wprtp2_implicit(i - 1, k - 1);
            });
      }

      // Calculate the momentum level coefficients and sign of vertical velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop18", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              coef_wprtp2_implicit_zm(i - 1, k - 1) =
                  one_third * beta(i - 1) * a1_coef(i - 1, k - 1) * wp3_on_wp2(i - 1, k - 1);
              sgn_t_vel_rtp2(i - 1, k - 1) = Kokkos::copysign(one, wp3_on_wp2(i - 1, k - 1));
            });
      }

      if (!l_godunov_upwind_xpyp_ta) {

        // Calculate the LHS turbulent advection term for <w'rt'^2>
        f2k::xpyp_term_ta_pdf_lhs(
            nzm, nzt, ngrdcol, gr, coef_wprtp2_implicit, rho_ds_zt, rho_ds_zm,          // In
            invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtp2, coef_wprtp2_implicit_zm, // In
            lhs_ta_wprtp2);                                                             // Out

      } else {

        // Godunov-like method for the vertical discretization of ta term
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop19", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzt);
              remaining /= (nzt);
              coef_wprtp2_implicit(i - 1, k - 1) =
                  one_third * beta(i - 1) * a1_coef_zt(i - 1, k - 1) * wp3_on_wp2_zt(i - 1, k - 1);
              coef_wpthlp2_implicit(i - 1, k - 1)   = coef_wprtp2_implicit(i - 1, k - 1);
              coef_wprtpthlp_implicit(i - 1, k - 1) = coef_wprtp2_implicit(i - 1, k - 1);
            });

        f2k::xpyp_term_ta_pdf_lhs_godunov(nzm, nzt, ngrdcol, gr, coef_wprtp2_implicit,
                                          invrs_rho_ds_zm, rho_ds_zm, // In
                                          lhs_ta_wprtp2);             // Out
      }

      // For ADG1, the LHS turbulent advection terms for
      // <w'rt'^2>, <w'thl'^2>, <w'rt'thl'>, and <w'sclr'x'> are all equal
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::loop20",
          std::max(0, nzm) * std::max(0, ngrdcol) * std::max(0, ndiags3),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int b   = 1 + remaining % (ndiags3);
            remaining /= (ndiags3);
            const int i = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            lhs_ta_wpthlp2(b - 1, i - 1, k - 1)   = lhs_ta_wprtp2(b - 1, i - 1, k - 1);
            lhs_ta_wprtpthlp(b - 1, i - 1, k - 1) = lhs_ta_wprtp2(b - 1, i - 1, k - 1);
          });

      if (l_scalar_calc) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop21",
            std::max(0, sclr_dim) * std::max(0, nzm) * std::max(0, ngrdcol) * std::max(0, ndiags3),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int b   = 1 + remaining % (ndiags3);
              remaining /= (ndiags3);
              const int i = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              const int sclr = 1 + remaining % (sclr_dim);
              remaining /= (sclr_dim);
              lhs_ta_wpsclrp2(b - 1, i - 1, k - 1, sclr - 1)   = lhs_ta_wprtp2(b - 1, i - 1, k - 1);
              lhs_ta_wprtpsclrp(b - 1, i - 1, k - 1, sclr - 1) = lhs_ta_wprtp2(b - 1, i - 1, k - 1);
              lhs_ta_wpthlpsclrp(b - 1, i - 1, k - 1, sclr - 1) =
                  lhs_ta_wprtp2(b - 1, i - 1, k - 1);
            });
      }

      // Explicit contributions

      // The termodynamic grid level term are only needed if l_upwind_xpyp_ta
      // is false, or if stats output is on

      if (!l_upwind_xpyp_ta || stats.flags().l_sample) {

        Kokkos::deep_copy(
            Kokkos::subview(wprtp_zt, Kokkos::ALL, Kokkos::ALL),
            zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wprtp, Kokkos::ALL, Kokkos::ALL)));
        Kokkos::deep_copy(
            Kokkos::subview(wpthlp_zt, Kokkos::ALL, Kokkos::ALL),
            zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wpthlp, Kokkos::ALL, Kokkos::ALL)));

        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop22", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzt);
              remaining /= (nzt);
              term_wprtp2_explicit(i - 1, k - 1) =
                  wp_coef_zt(i - 1, k - 1) * (wprtp_zt(i - 1, k - 1) * wprtp_zt(i - 1, k - 1));
              term_wpthlp2_explicit(i - 1, k - 1) =
                  wp_coef_zt(i - 1, k - 1) * (wpthlp_zt(i - 1, k - 1) * wpthlp_zt(i - 1, k - 1));
              term_wprtpthlp_explicit(i - 1, k - 1) =
                  wp_coef_zt(i - 1, k - 1) * wprtp_zt(i - 1, k - 1) * wpthlp_zt(i - 1, k - 1);
            });
      }

      // Calculate the momentum level terms and sign of vertical velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop23", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              term_wprtp2_explicit_zm(i - 1, k - 1) =
                  wp_coef(i - 1, k - 1) * (wprtp(i - 1, k - 1) * wprtp(i - 1, k - 1));
              sgn_t_vel_rtp2(i - 1, k - 1) = Kokkos::copysign(one, wp3_on_wp2(i - 1, k - 1));
            });
      }

      if (!l_godunov_upwind_xpyp_ta) {

        // Calculate the RHS turbulent advection term for <w'rt'^2>
        f2k::xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wprtp2_explicit, rho_ds_zt, rho_ds_zm,          // In
            invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtp2, term_wprtp2_explicit_zm, // In
            rhs_ta_wprtp2);                                                             // Out

      } else {

        // Using the godunov upwind scheme for the calculation of RHS turbulent
        // advection term for <w'rt'^2>. Here, we define the "wind" for godunov
        // scheme as ( one - one_third * beta ) * a1_coef_zt**2 * wp3_on_wp2_zt / wp2_zt,
        // and define the xpyp_term_ta_pdf_rhs_godunov subroutine in
        // turbulent_adv_pdf.F90 to process the calculation using godunov scheme
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop24", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              term_wprtp2_explicit_zm(i - 1, k - 1) = (wprtp(i - 1, k - 1) * wprtp(i - 1, k - 1));
            });

        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop25", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzt);
              remaining /= (nzt);
              sgn_t_vel_rtp2_zt(i - 1, k - 1) = Kokkos::copysign(one, wp_coef_zt(i - 1, k - 1));
            });

        f2k::xpyp_term_ta_pdf_rhs_godunov(nzm, nzt, ngrdcol, gr, term_wprtp2_explicit_zm,
                                          invrs_rho_ds_zm,              // In
                                          sgn_t_vel_rtp2_zt, rho_ds_zm, // In
                                          rhs_ta_wprtp2);               // Out
      }

      // Calculate the momentum level terms and sign of vertical velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop26", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              term_wpthlp2_explicit_zm(i - 1, k - 1) =
                  wp_coef(i - 1, k - 1) * (wpthlp(i - 1, k - 1) * wpthlp(i - 1, k - 1));
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop27", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzt);
              remaining /= (nzt);
              sgn_t_vel_thlp2(i - 1, k - 1) = Kokkos::copysign(one, wp3_on_wp2(i - 1, k - 1));
            });
      }

      if (!l_godunov_upwind_xpyp_ta) {

        // Calculate the RHS turbulent advection term for <w'thl'^2>
        f2k::xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wpthlp2_explicit, rho_ds_zt, rho_ds_zm,           // In
            invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_thlp2, term_wpthlp2_explicit_zm, // In
            rhs_ta_wpthlp2);                                                              // Out

      } else {

        // Using the godunov upwind scheme for the calculation of RHS
        // turbulent advection term for <w'thl'^2>.
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop28", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              term_wpthlp2_explicit_zm(i - 1, k - 1) =
                  (wpthlp(i - 1, k - 1) * wpthlp(i - 1, k - 1));
            });

        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop29", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzt);
              remaining /= (nzt);
              sgn_t_vel_thlp2_zt(i - 1, k - 1) = Kokkos::copysign(one, wp_coef_zt(i - 1, k - 1));
            });

        f2k::xpyp_term_ta_pdf_rhs_godunov(nzm, nzt, ngrdcol, gr, term_wpthlp2_explicit_zm,
                                          invrs_rho_ds_zm,               // In
                                          sgn_t_vel_thlp2_zt, rho_ds_zm, // In
                                          rhs_ta_wpthlp2);               // Out
      }

      // Calculate the momentum level terms and sign of vertical velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop30", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              term_wprtpthlp_explicit_zm(i - 1, k - 1) =
                  wp_coef(i - 1, k - 1) * wprtp(i - 1, k - 1) * wpthlp(i - 1, k - 1);
              sgn_t_vel_rtpthlp(i - 1, k - 1) = Kokkos::copysign(one, wp3_on_wp2(i - 1, k - 1));
            });
      }

      if (!l_godunov_upwind_xpyp_ta) {

        // Calculate the RHS turbulent advection term for <w'rt'thl'>
        f2k::xpyp_term_ta_pdf_rhs(
            nzm, nzt, ngrdcol, gr, term_wprtpthlp_explicit, rho_ds_zt, rho_ds_zm,             // In
            invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, term_wprtpthlp_explicit_zm, // In
            rhs_ta_wprtpthlp);                                                                // Out

      } else {

        // Using the godunov upwind scheme for the calculation of RHS
        // turbulent
        // advection term for <w'rt'thl'>.
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop31", std::max(0, nzm) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzm);
              remaining /= (nzm);
              term_wprtpthlp_explicit_zm(i - 1, k - 1) = wprtp(i - 1, k - 1) * wpthlp(i - 1, k - 1);
            });

        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::loop32", std::max(0, nzt) * std::max(0, ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              int remaining = cell;
              const int i   = 1 + remaining % (ngrdcol);
              remaining /= (ngrdcol);
              const int k = 1 + remaining % (nzt);
              remaining /= (nzt);
              sgn_t_vel_rtpthlp_zt(i - 1, k - 1) = Kokkos::copysign(one, wp_coef_zt(i - 1, k - 1));
            });

        f2k::xpyp_term_ta_pdf_rhs_godunov(nzm, nzt, ngrdcol, gr, term_wprtpthlp_explicit_zm,
                                          invrs_rho_ds_zm,                 // In
                                          sgn_t_vel_rtpthlp_zt, rho_ds_zm, // In
                                          rhs_ta_wprtpthlp);               // Out
      }

      if (l_scalar_calc) {

        // Explicit contributions to passive scalars

        // Interpolate wpsclrp to thermo levels if not using l_upwind_xpyp_ta
        if (!l_upwind_xpyp_ta) {
          for (sclr = 1; sclr <= sclr_dim; sclr += 1) {
            Kokkos::deep_copy(
                Kokkos::subview(wpsclrp_zt, Kokkos::ALL, Kokkos::ALL, (sclr - 1)),
                zm2zt_api(nzm, nzt, ngrdcol, gr,
                          Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))));
          }
        }

        for (sclr = 1; sclr <= sclr_dim; sclr += 1) {

          // Calculate the momentum level terms and sign of vertical velocity if
          // l_upwind_xpyp_ta is true, otherwise just calculate the thermo level terms
          if (l_upwind_xpyp_ta) {
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop33", std::max(0, nzm) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzm);
                  remaining /= (nzm);
                  term_wpsclrp2_explicit_zm(i - 1, k - 1) =
                      wp_coef(i - 1, k - 1) *
                      (wpsclrp(i - 1, k - 1, sclr - 1) * wpsclrp(i - 1, k - 1, sclr - 1));
                  sgn_t_vel_sclrp2(i - 1, k - 1) = Kokkos::copysign(one, wp3_on_wp2(i - 1, k - 1));
                });
          } else {
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop34", std::max(0, nzt) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzt);
                  remaining /= (nzt);
                  term_wpsclrp2_explicit(i - 1, k - 1) =
                      wp_coef_zt(i - 1, k - 1) *
                      (wpsclrp_zt(i - 1, k - 1, sclr - 1) * wpsclrp_zt(i - 1, k - 1, sclr - 1));
                });
          }

          if (!l_godunov_upwind_xpyp_ta) {

            // Calculate the RHS turbulent advection term for <w'sclr'^2>
            f2k::xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wpsclrp2_explicit, rho_ds_zt, rho_ds_zm, // In
                invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrp2,
                term_wpsclrp2_explicit_zm,                                               // In
                Kokkos::subview(rhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out

          } else {

            // Using the godunov upwind scheme for the calculation of RHS
            // turbulent
            // advection term for <w'sclr'^2>.
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop35", std::max(0, nzm) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzm);
                  remaining /= (nzm);
                  term_wpsclrp2_explicit_zm(i - 1, k - 1) =
                      (wpsclrp(i - 1, k - 1, sclr - 1) * wpsclrp(i - 1, k - 1, sclr - 1));
                });
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop36", std::max(0, nzt) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzt);
                  remaining /= (nzt);
                  sgn_t_vel_sclrp2_zt(i - 1, k - 1) =
                      Kokkos::copysign(one, wp_coef_zt(i - 1, k - 1));
                });

            f2k::xpyp_term_ta_pdf_rhs_godunov(
                nzm, nzt, ngrdcol, gr, term_wpsclrp2_explicit_zm, invrs_rho_ds_zm,       // In
                sgn_t_vel_sclrp2_zt, rho_ds_zm,                                          // In
                Kokkos::subview(rhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out
          }
        }

        // Calculate the momentum level terms and sign of vertical velocity if
        // l_upwind_xpyp_ta is true, otherwise just calculate the thermo level terms
        for (sclr = 1; sclr <= sclr_dim; sclr += 1) {

          if (l_upwind_xpyp_ta) {
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop37", std::max(0, nzm) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzm);
                  remaining /= (nzm);
                  term_wprtpsclrp_explicit_zm(i - 1, k - 1) =
                      wp_coef(i - 1, k - 1) * wpsclrp(i - 1, k - 1, sclr - 1) * wprtp(i - 1, k - 1);
                  sgn_t_vel_sclrprtp(i - 1, k - 1) =
                      Kokkos::copysign(one, wp3_on_wp2(i - 1, k - 1));
                });
          } else {
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop38", std::max(0, nzt) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzt);
                  remaining /= (nzt);
                  term_wprtpsclrp_explicit(i - 1, k - 1) = wp_coef_zt(i - 1, k - 1) *
                                                           wpsclrp_zt(i - 1, k - 1, sclr - 1) *
                                                           wprtp_zt(i - 1, k - 1);
                });
          }

          if (!l_godunov_upwind_xpyp_ta) {

            // Calculate the RHS turbulent advection term for <w'sclr'rt'>
            f2k::xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wprtpsclrp_explicit, rho_ds_zt, rho_ds_zm,     // In
                invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrprtp,                     // In
                term_wprtpsclrp_explicit_zm,                                               // In
                Kokkos::subview(rhs_ta_wprtpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out

          } else {

            // Using the godunov upwind scheme for the calculation of RHS
            // turbulent advection term for <w'sclr'rt'>.
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop39", std::max(0, nzm) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzm);
                  remaining /= (nzm);
                  term_wprtpsclrp_explicit_zm(i - 1, k - 1) =
                      wpsclrp(i - 1, k - 1, sclr - 1) * wprtp(i - 1, k - 1);
                });
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop40", std::max(0, nzt) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzt);
                  remaining /= (nzt);
                  sgn_t_vel_sclrprtp_zt(i - 1, k - 1) =
                      Kokkos::copysign(one, wp_coef_zt(i - 1, k - 1));
                });

            f2k::xpyp_term_ta_pdf_rhs_godunov(
                nzm, nzt, ngrdcol, gr, term_wprtpsclrp_explicit_zm, invrs_rho_ds_zm,       // In
                sgn_t_vel_sclrprtp_zt, rho_ds_zm,                                          // In
                Kokkos::subview(rhs_ta_wprtpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out
          }
        }

        // Calculate the momentum level terms and sign of vertical velocity if
        // l_upwind_xpyp_ta is true, otherwise just calculate the thermo level terms
        for (sclr = 1; sclr <= sclr_dim; sclr += 1) {

          if (l_upwind_xpyp_ta) {
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop41", std::max(0, nzm) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzm);
                  remaining /= (nzm);
                  term_wpthlpsclrp_explicit_zm(i - 1, k - 1) = wp_coef(i - 1, k - 1) *
                                                               wpsclrp(i - 1, k - 1, sclr - 1) *
                                                               wpthlp(i - 1, k - 1);
                  sgn_t_vel_sclrpthlp(i - 1, k - 1) =
                      Kokkos::copysign(one, wp3_on_wp2(i - 1, k - 1));
                });
          } else {
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop42", std::max(0, nzt) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzt);
                  remaining /= (nzt);
                  term_wpthlpsclrp_explicit(i - 1, k - 1) = wp_coef_zt(i - 1, k - 1) *
                                                            wpsclrp_zt(i - 1, k - 1, sclr - 1) *
                                                            wpthlp_zt(i - 1, k - 1);
                });
          }

          if (!l_godunov_upwind_xpyp_ta) {

            // Calculate the RHS turbulent advection term for <w'sclr'thl'>
            f2k::xpyp_term_ta_pdf_rhs(
                nzm, nzt, ngrdcol, gr, term_wpthlpsclrp_explicit, rho_ds_zt, rho_ds_zm,     // In
                invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrpthlp,                     // In
                term_wpthlpsclrp_explicit_zm,                                               // In
                Kokkos::subview(rhs_ta_wpthlpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out

          } else {

            // Using the godunov upwind scheme for the calculation of RHS
            // turbulent advection term for <w'sclr'thl'>.
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop43", std::max(0, nzm) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzm);
                  remaining /= (nzm);
                  term_wpthlpsclrp_explicit_zm(i - 1, k - 1) =
                      wpsclrp(i - 1, k - 1, sclr - 1) * wpthlp(i - 1, k - 1);
                });
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::loop44", std::max(0, nzt) * std::max(0, ngrdcol),
                KOKKOS_LAMBDA(const int cell) {
                  int remaining = cell;
                  const int i   = 1 + remaining % (ngrdcol);
                  remaining /= (ngrdcol);
                  const int k = 1 + remaining % (nzt);
                  remaining /= (nzt);
                  sgn_t_vel_sclrpthlp_zt(i - 1, k - 1) =
                      Kokkos::copysign(one, wp_coef_zt(i - 1, k - 1));
                });

            f2k::xpyp_term_ta_pdf_rhs_godunov(
                nzm, nzt, ngrdcol, gr, term_wpthlpsclrp_explicit_zm, invrs_rho_ds_zm,       // In
                sgn_t_vel_sclrpthlp_zt, rho_ds_zm,                                          // In
                Kokkos::subview(rhs_ta_wpthlpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out
          }
        }
      }

    } else if (iiPDF_type == iiPDF_new) {

      // The new PDF is used.

      // The termodynamic grid level coefficients are only needed if l_upwind_xpyp_ta
      // is false, or if stats output is on

      if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array45",
            ((coef_wprtp2_implicit.extent(0) - (1)) / (1) + 1) *
                ((coef_wprtp2_implicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((coef_wprtp2_implicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((coef_wprtp2_implicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((coef_wprtp2_implicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((coef_wprtp2_implicit.extent(1) - (1)) / (1) + 1);
              coef_wprtp2_implicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_coef_wprtp2_implicit((1 - 1 + slice_0 * (1)),
                                                                (1 - 1 + slice_1 * (1)));
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array46",
            ((coef_wpthlp2_implicit.extent(0) - (1)) / (1) + 1) *
                ((coef_wpthlp2_implicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((coef_wpthlp2_implicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((coef_wpthlp2_implicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((coef_wpthlp2_implicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((coef_wpthlp2_implicit.extent(1) - (1)) / (1) + 1);
              coef_wpthlp2_implicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_coef_wpthlp2_implicit((1 - 1 + slice_0 * (1)),
                                                                 (1 - 1 + slice_1 * (1)));
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array47",
            ((coef_wprtpthlp_implicit.extent(0) - (1)) / (1) + 1) *
                ((coef_wprtpthlp_implicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((coef_wprtpthlp_implicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((coef_wprtpthlp_implicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((coef_wprtpthlp_implicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((coef_wprtpthlp_implicit.extent(1) - (1)) / (1) + 1);
              coef_wprtpthlp_implicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_coef_wprtpthlp_implicit((1 - 1 + slice_0 * (1)),
                                                                   (1 - 1 + slice_1 * (1)));
            });
      }

      // Calculate the momentum level terms and sign of vertical velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::deep_copy(Kokkos::subview(coef_wprtp2_implicit_zm, Kokkos::ALL, Kokkos::ALL),
                          zt2zm_api(nzm, nzt, ngrdcol, gr,
                                    Kokkos::subview(pdf_implicit_coefs_terms.coef_wprtp2_implicit,
                                                    Kokkos::ALL, Kokkos::ALL)));
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array48",
            ((sgn_t_vel_rtp2.extent(0) - (1)) / (1) + 1) *
                ((sgn_t_vel_rtp2.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((sgn_t_vel_rtp2.extent(0) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_rtp2.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((sgn_t_vel_rtp2.extent(1) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_rtp2.extent(1) - (1)) / (1) + 1);
              sgn_t_vel_rtp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = Kokkos::copysign(
                  one, coef_wprtp2_implicit_zm((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) *
                           rtp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) *
                           rtp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))));
            });
      }

      // Calculate the LHS turbulent advection term for <w'rt'^2>
      f2k::xpyp_term_ta_pdf_lhs(
          nzm, nzt, ngrdcol, gr, coef_wprtp2_implicit, rho_ds_zt, rho_ds_zm,          // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtp2, coef_wprtp2_implicit_zm, // In
          lhs_ta_wprtp2);                                                             // Out

      // Calculate the momentum level terms and sign of vertical velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::deep_copy(Kokkos::subview(coef_wpthlp2_implicit_zm, Kokkos::ALL, Kokkos::ALL),
                          zt2zm_api(nzm, nzt, ngrdcol, gr,
                                    Kokkos::subview(pdf_implicit_coefs_terms.coef_wpthlp2_implicit,
                                                    Kokkos::ALL, Kokkos::ALL)));
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array49",
            ((sgn_t_vel_thlp2.extent(0) - (1)) / (1) + 1) *
                ((sgn_t_vel_thlp2.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((sgn_t_vel_thlp2.extent(0) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_thlp2.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((sgn_t_vel_thlp2.extent(1) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_thlp2.extent(1) - (1)) / (1) + 1);
              sgn_t_vel_thlp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = Kokkos::copysign(
                  one, coef_wpthlp2_implicit_zm((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) *
                           thlp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) *
                           thlp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))));
            });
      }

      // Calculate the LHS turbulent advection term for <w'thl'^2>
      f2k::xpyp_term_ta_pdf_lhs(
          nzm, nzt, ngrdcol, gr, coef_wpthlp2_implicit, rho_ds_zt, rho_ds_zm,           // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_thlp2, coef_wpthlp2_implicit_zm, // In
          lhs_ta_wpthlp2);                                                              // Out

      // Calculate the momentum level terms and sign of vertical velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::deep_copy(
            Kokkos::subview(coef_wprtpthlp_implicit_zm, Kokkos::ALL, Kokkos::ALL),
            zt2zm_api(nzm, nzt, ngrdcol, gr,
                      Kokkos::subview(pdf_implicit_coefs_terms.coef_wprtpthlp_implicit, Kokkos::ALL,
                                      Kokkos::ALL)));
        Kokkos::deep_copy(
            Kokkos::subview(term_wprtpthlp_explicit_zm, Kokkos::ALL, Kokkos::ALL),
            zt2zm_api(nzm, nzt, ngrdcol, gr,
                      Kokkos::subview(pdf_implicit_coefs_terms.term_wprtpthlp_explicit, Kokkos::ALL,
                                      Kokkos::ALL)));
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array50",
            ((sgn_t_vel_rtpthlp.extent(0) - (1)) / (1) + 1) *
                ((sgn_t_vel_rtpthlp.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((sgn_t_vel_rtpthlp.extent(0) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_rtpthlp.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((sgn_t_vel_rtpthlp.extent(1) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_rtpthlp.extent(1) - (1)) / (1) + 1);
              sgn_t_vel_rtpthlp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  Kokkos::copysign(one,
                                   (coef_wprtpthlp_implicit_zm((1 - 1 + slice_0 * (1)),
                                                               (1 - 1 + slice_1 * (1))) *
                                        rtpthlp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) +
                                    term_wprtpthlp_explicit_zm((1 - 1 + slice_0 * (1)),
                                                               (1 - 1 + slice_1 * (1)))) *
                                       rtpthlp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))));
            });
      }

      // Calculate the LHS turbulent advection term for <w'rt'thl'>
      f2k::xpyp_term_ta_pdf_lhs(
          nzm, nzt, ngrdcol, gr, coef_wprtpthlp_implicit, rho_ds_zt, rho_ds_zm,             // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, coef_wprtpthlp_implicit_zm, // In
          lhs_ta_wprtpthlp);                                                                // Out

      if (l_scalar_calc) {
        // The code for the scalar variables will be set up later.
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array51",
            ((lhs_ta_wpsclrp2.extent(0) - (1)) / (1) + 1) *
                ((lhs_ta_wpsclrp2.extent(1) - (1)) / (1) + 1) *
                ((lhs_ta_wpsclrp2.extent(2) - (1)) / (1) + 1) *
                ((lhs_ta_wpsclrp2.extent(3) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((lhs_ta_wpsclrp2.extent(0) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wpsclrp2.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((lhs_ta_wpsclrp2.extent(1) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wpsclrp2.extent(1) - (1)) / (1) + 1);
              const int slice_2 = remaining % ((lhs_ta_wpsclrp2.extent(2) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wpsclrp2.extent(2) - (1)) / (1) + 1);
              const int slice_3 = remaining % ((lhs_ta_wpsclrp2.extent(3) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wpsclrp2.extent(3) - (1)) / (1) + 1);
              lhs_ta_wpsclrp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)),
                              (1 - 1 + slice_2 * (1)), (1 - 1 + slice_3 * (1))) = zero;
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array52",
            ((lhs_ta_wprtpsclrp.extent(0) - (1)) / (1) + 1) *
                ((lhs_ta_wprtpsclrp.extent(1) - (1)) / (1) + 1) *
                ((lhs_ta_wprtpsclrp.extent(2) - (1)) / (1) + 1) *
                ((lhs_ta_wprtpsclrp.extent(3) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((lhs_ta_wprtpsclrp.extent(0) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wprtpsclrp.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((lhs_ta_wprtpsclrp.extent(1) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wprtpsclrp.extent(1) - (1)) / (1) + 1);
              const int slice_2 = remaining % ((lhs_ta_wprtpsclrp.extent(2) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wprtpsclrp.extent(2) - (1)) / (1) + 1);
              const int slice_3 = remaining % ((lhs_ta_wprtpsclrp.extent(3) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wprtpsclrp.extent(3) - (1)) / (1) + 1);
              lhs_ta_wprtpsclrp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)),
                                (1 - 1 + slice_2 * (1)), (1 - 1 + slice_3 * (1))) = zero;
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array53",
            ((lhs_ta_wpthlpsclrp.extent(0) - (1)) / (1) + 1) *
                ((lhs_ta_wpthlpsclrp.extent(1) - (1)) / (1) + 1) *
                ((lhs_ta_wpthlpsclrp.extent(2) - (1)) / (1) + 1) *
                ((lhs_ta_wpthlpsclrp.extent(3) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((lhs_ta_wpthlpsclrp.extent(0) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wpthlpsclrp.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((lhs_ta_wpthlpsclrp.extent(1) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wpthlpsclrp.extent(1) - (1)) / (1) + 1);
              const int slice_2 = remaining % ((lhs_ta_wpthlpsclrp.extent(2) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wpthlpsclrp.extent(2) - (1)) / (1) + 1);
              const int slice_3 = remaining % ((lhs_ta_wpthlpsclrp.extent(3) - (1)) / (1) + 1);
              remaining /= ((lhs_ta_wpthlpsclrp.extent(3) - (1)) / (1) + 1);
              lhs_ta_wpthlpsclrp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)),
                                 (1 - 1 + slice_2 * (1)), (1 - 1 + slice_3 * (1))) = zero;
            });
      }

      // The termodynamic grid level term are only needed if l_upwind_xpyp_ta
      // is false, or if stats output is on. The value of term_wprtp2_explicit_zm
      // and term_wpthlp2_explicit are always 0.

      if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array54",
            ((term_wprtp2_explicit.extent(0) - (1)) / (1) + 1) *
                ((term_wprtp2_explicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((term_wprtp2_explicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((term_wprtp2_explicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((term_wprtp2_explicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((term_wprtp2_explicit.extent(1) - (1)) / (1) + 1);
              term_wprtp2_explicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = zero;
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array55",
            ((term_wpthlp2_explicit.extent(0) - (1)) / (1) + 1) *
                ((term_wpthlp2_explicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((term_wpthlp2_explicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((term_wpthlp2_explicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((term_wpthlp2_explicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((term_wpthlp2_explicit.extent(1) - (1)) / (1) + 1);
              term_wpthlp2_explicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = zero;
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array56",
            ((term_wprtpthlp_explicit.extent(0) - (1)) / (1) + 1) *
                ((term_wprtpthlp_explicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((term_wprtpthlp_explicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((term_wprtpthlp_explicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((term_wprtpthlp_explicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((term_wprtpthlp_explicit.extent(1) - (1)) / (1) + 1);
              term_wprtpthlp_explicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_term_wprtpthlp_explicit((1 - 1 + slice_0 * (1)),
                                                                   (1 - 1 + slice_1 * (1)));
            });
      }

      // Calculate the RHS turbulent advection term for <w'rt'thl'>
      f2k::xpyp_term_ta_pdf_rhs(
          nzm, nzt, ngrdcol, gr, term_wprtpthlp_explicit, rho_ds_zt, rho_ds_zm,             // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, term_wprtpthlp_explicit_zm, // In
          rhs_ta_wprtpthlp);                                                                // Out

      // The <rt'^2> and <thl'^2> turbulent advection terms are entirely implicit, as
      // <w'rt'^2> = coef_wprtp2_implicit * <rt'^2>, and
      // <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2>.  So the values of these RHS
      // turbulent advection terms are always zero
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::array57",
          ((rhs_ta_wprtp2.extent(0) - (1)) / (1) + 1) * ((rhs_ta_wprtp2.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((rhs_ta_wprtp2.extent(0) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wprtp2.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((rhs_ta_wprtp2.extent(1) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wprtp2.extent(1) - (1)) / (1) + 1);
            rhs_ta_wprtp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = zero;
          });
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::array58",
          ((rhs_ta_wpthlp2.extent(0) - (1)) / (1) + 1) *
              ((rhs_ta_wpthlp2.extent(1) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((rhs_ta_wpthlp2.extent(0) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wpthlp2.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((rhs_ta_wpthlp2.extent(1) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wpthlp2.extent(1) - (1)) / (1) + 1);
            rhs_ta_wpthlp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = zero;
          });

      // The code for the scalar variables will be set up later.
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::array59",
          ((rhs_ta_wpsclrp2.extent(0) - (1)) / (1) + 1) *
              ((rhs_ta_wpsclrp2.extent(1) - (1)) / (1) + 1) *
              ((rhs_ta_wpsclrp2.extent(2) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((rhs_ta_wpsclrp2.extent(0) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wpsclrp2.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((rhs_ta_wpsclrp2.extent(1) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wpsclrp2.extent(1) - (1)) / (1) + 1);
            const int slice_2 = remaining % ((rhs_ta_wpsclrp2.extent(2) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wpsclrp2.extent(2) - (1)) / (1) + 1);
            rhs_ta_wpsclrp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)),
                            (1 - 1 + slice_2 * (1))) = zero;
          });
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::array60",
          ((rhs_ta_wprtpsclrp.extent(0) - (1)) / (1) + 1) *
              ((rhs_ta_wprtpsclrp.extent(1) - (1)) / (1) + 1) *
              ((rhs_ta_wprtpsclrp.extent(2) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((rhs_ta_wprtpsclrp.extent(0) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wprtpsclrp.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((rhs_ta_wprtpsclrp.extent(1) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wprtpsclrp.extent(1) - (1)) / (1) + 1);
            const int slice_2 = remaining % ((rhs_ta_wprtpsclrp.extent(2) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wprtpsclrp.extent(2) - (1)) / (1) + 1);
            rhs_ta_wprtpsclrp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)),
                              (1 - 1 + slice_2 * (1))) = zero;
          });
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::array61",
          ((rhs_ta_wpthlpsclrp.extent(0) - (1)) / (1) + 1) *
              ((rhs_ta_wpthlpsclrp.extent(1) - (1)) / (1) + 1) *
              ((rhs_ta_wpthlpsclrp.extent(2) - (1)) / (1) + 1),
          KOKKOS_LAMBDA(int cell) {
            int remaining     = cell;
            const int slice_0 = remaining % ((rhs_ta_wpthlpsclrp.extent(0) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wpthlpsclrp.extent(0) - (1)) / (1) + 1);
            const int slice_1 = remaining % ((rhs_ta_wpthlpsclrp.extent(1) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wpthlpsclrp.extent(1) - (1)) / (1) + 1);
            const int slice_2 = remaining % ((rhs_ta_wpthlpsclrp.extent(2) - (1)) / (1) + 1);
            remaining /= ((rhs_ta_wpthlpsclrp.extent(2) - (1)) / (1) + 1);
            rhs_ta_wpthlpsclrp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)),
                               (1 - 1 + slice_2 * (1))) = zero;
          });

    } else if (iiPDF_type == iiPDF_new_hybrid) {

      // The new hybrid PDF is used.

      if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array62",
            ((coef_wprtp2_implicit.extent(0) - (1)) / (1) + 1) *
                ((coef_wprtp2_implicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((coef_wprtp2_implicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((coef_wprtp2_implicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((coef_wprtp2_implicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((coef_wprtp2_implicit.extent(1) - (1)) / (1) + 1);
              coef_wprtp2_implicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_coef_wprtp2_implicit((1 - 1 + slice_0 * (1)),
                                                                (1 - 1 + slice_1 * (1)));
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array63",
            ((term_wprtp2_explicit.extent(0) - (1)) / (1) + 1) *
                ((term_wprtp2_explicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((term_wprtp2_explicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((term_wprtp2_explicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((term_wprtp2_explicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((term_wprtp2_explicit.extent(1) - (1)) / (1) + 1);
              term_wprtp2_explicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_term_wprtp2_explicit((1 - 1 + slice_0 * (1)),
                                                                (1 - 1 + slice_1 * (1)));
            });
      }

      // Calculate the momentum level terms and sign of turbulent velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::deep_copy(Kokkos::subview(coef_wprtp2_implicit_zm, Kokkos::ALL, Kokkos::ALL),
                          zt2zm_api(nzm, nzt, ngrdcol, gr,
                                    Kokkos::subview(pdf_implicit_coefs_terms.coef_wprtp2_implicit,
                                                    Kokkos::ALL, Kokkos::ALL)));
        Kokkos::deep_copy(Kokkos::subview(term_wprtp2_explicit_zm, Kokkos::ALL, Kokkos::ALL),
                          zt2zm_api(nzm, nzt, ngrdcol, gr,
                                    Kokkos::subview(pdf_implicit_coefs_terms.term_wprtp2_explicit,
                                                    Kokkos::ALL, Kokkos::ALL)));
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array64",
            ((sgn_t_vel_rtp2.extent(0) - (1)) / (1) + 1) *
                ((sgn_t_vel_rtp2.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((sgn_t_vel_rtp2.extent(0) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_rtp2.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((sgn_t_vel_rtp2.extent(1) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_rtp2.extent(1) - (1)) / (1) + 1);
              sgn_t_vel_rtp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = Kokkos::copysign(
                  one, (coef_wprtp2_implicit_zm((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) *
                            rtp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) +
                        term_wprtp2_explicit_zm((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)))) *
                           rtp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))));
            });
      }

      // Calculate the LHS turbulent advection term for <w'rt'^2>
      f2k::xpyp_term_ta_pdf_lhs(
          nzm, nzt, ngrdcol, gr, coef_wprtp2_implicit, rho_ds_zt, rho_ds_zm,          // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtp2, coef_wprtp2_implicit_zm, // In
          lhs_ta_wprtp2);                                                             // Out

      // Calculate the RHS turbulent advection term for <w'rt'^2>
      f2k::xpyp_term_ta_pdf_rhs(
          nzm, nzt, ngrdcol, gr, term_wprtp2_explicit, rho_ds_zt, rho_ds_zm,          // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtp2, term_wprtp2_explicit_zm, // In
          rhs_ta_wprtp2);                                                             // Out

      if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array65",
            ((coef_wpthlp2_implicit.extent(0) - (1)) / (1) + 1) *
                ((coef_wpthlp2_implicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((coef_wpthlp2_implicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((coef_wpthlp2_implicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((coef_wpthlp2_implicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((coef_wpthlp2_implicit.extent(1) - (1)) / (1) + 1);
              coef_wpthlp2_implicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_coef_wpthlp2_implicit((1 - 1 + slice_0 * (1)),
                                                                 (1 - 1 + slice_1 * (1)));
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array66",
            ((term_wpthlp2_explicit.extent(0) - (1)) / (1) + 1) *
                ((term_wpthlp2_explicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((term_wpthlp2_explicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((term_wpthlp2_explicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((term_wpthlp2_explicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((term_wpthlp2_explicit.extent(1) - (1)) / (1) + 1);
              term_wpthlp2_explicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_term_wpthlp2_explicit((1 - 1 + slice_0 * (1)),
                                                                 (1 - 1 + slice_1 * (1)));
            });
      }

      // Calculate the momentum level terms and sign of turbulent velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::deep_copy(Kokkos::subview(coef_wpthlp2_implicit_zm, Kokkos::ALL, Kokkos::ALL),
                          zt2zm_api(nzm, nzt, ngrdcol, gr,
                                    Kokkos::subview(pdf_implicit_coefs_terms.coef_wpthlp2_implicit,
                                                    Kokkos::ALL, Kokkos::ALL)));
        Kokkos::deep_copy(Kokkos::subview(term_wpthlp2_explicit_zm, Kokkos::ALL, Kokkos::ALL),
                          zt2zm_api(nzm, nzt, ngrdcol, gr,
                                    Kokkos::subview(pdf_implicit_coefs_terms.term_wpthlp2_explicit,
                                                    Kokkos::ALL, Kokkos::ALL)));
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array67",
            ((sgn_t_vel_thlp2.extent(0) - (1)) / (1) + 1) *
                ((sgn_t_vel_thlp2.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((sgn_t_vel_thlp2.extent(0) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_thlp2.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((sgn_t_vel_thlp2.extent(1) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_thlp2.extent(1) - (1)) / (1) + 1);
              sgn_t_vel_thlp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = Kokkos::copysign(
                  one,
                  (coef_wpthlp2_implicit_zm((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) *
                       thlp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) +
                   term_wpthlp2_explicit_zm((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)))) *
                      thlp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))));
            });
      }

      // Calculate the LHS turbulent advection term for <w'thl'^2>
      f2k::xpyp_term_ta_pdf_lhs(
          nzm, nzt, ngrdcol, gr, coef_wpthlp2_implicit, rho_ds_zt, rho_ds_zm,           // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_thlp2, coef_wpthlp2_implicit_zm, // In
          lhs_ta_wpthlp2);                                                              // Out

      // Calculate the RHS turbulent advection term for <w'thl'^2>
      f2k::xpyp_term_ta_pdf_rhs(
          nzm, nzt, ngrdcol, gr, term_wpthlp2_explicit, rho_ds_zt, rho_ds_zm,           // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_thlp2, term_wpthlp2_explicit_zm, // In
          rhs_ta_wpthlp2);                                                              // Out

      if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array68",
            ((coef_wprtpthlp_implicit.extent(0) - (1)) / (1) + 1) *
                ((coef_wprtpthlp_implicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((coef_wprtpthlp_implicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((coef_wprtpthlp_implicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((coef_wprtpthlp_implicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((coef_wprtpthlp_implicit.extent(1) - (1)) / (1) + 1);
              coef_wprtpthlp_implicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_coef_wprtpthlp_implicit((1 - 1 + slice_0 * (1)),
                                                                   (1 - 1 + slice_1 * (1)));
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array69",
            ((term_wprtpthlp_explicit.extent(0) - (1)) / (1) + 1) *
                ((term_wprtpthlp_explicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((term_wprtpthlp_explicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((term_wprtpthlp_explicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((term_wprtpthlp_explicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((term_wprtpthlp_explicit.extent(1) - (1)) / (1) + 1);
              term_wprtpthlp_explicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_term_wprtpthlp_explicit((1 - 1 + slice_0 * (1)),
                                                                   (1 - 1 + slice_1 * (1)));
            });
      }

      // Calculate the momentum level terms and sign of turbulent velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::deep_copy(
            Kokkos::subview(coef_wprtpthlp_implicit_zm, Kokkos::ALL, Kokkos::ALL),
            zt2zm_api(nzm, nzt, ngrdcol, gr,
                      Kokkos::subview(pdf_implicit_coefs_terms.coef_wprtpthlp_implicit, Kokkos::ALL,
                                      Kokkos::ALL)));
        Kokkos::deep_copy(
            Kokkos::subview(term_wprtpthlp_explicit_zm, Kokkos::ALL, Kokkos::ALL),
            zt2zm_api(nzm, nzt, ngrdcol, gr,
                      Kokkos::subview(pdf_implicit_coefs_terms.term_wprtpthlp_explicit, Kokkos::ALL,
                                      Kokkos::ALL)));
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array70",
            ((sgn_t_vel_rtpthlp.extent(0) - (1)) / (1) + 1) *
                ((sgn_t_vel_rtpthlp.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((sgn_t_vel_rtpthlp.extent(0) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_rtpthlp.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((sgn_t_vel_rtpthlp.extent(1) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_rtpthlp.extent(1) - (1)) / (1) + 1);
              sgn_t_vel_rtpthlp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  Kokkos::copysign(one,
                                   (coef_wprtpthlp_implicit_zm((1 - 1 + slice_0 * (1)),
                                                               (1 - 1 + slice_1 * (1))) *
                                        rtpthlp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) +
                                    term_wprtpthlp_explicit_zm((1 - 1 + slice_0 * (1)),
                                                               (1 - 1 + slice_1 * (1)))) *
                                       rtpthlp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))));
            });
      }

      // Calculate the LHS turbulent advection term for <w'rt'thl'>
      f2k::xpyp_term_ta_pdf_lhs(
          nzm, nzt, ngrdcol, gr, coef_wprtpthlp_implicit, rho_ds_zt, rho_ds_zm,             // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, coef_wprtpthlp_implicit_zm, // In
          lhs_ta_wprtpthlp);                                                                // Out

      // Calculate the RHS turbulent advection term for <w'rt'thl'>
      f2k::xpyp_term_ta_pdf_rhs(
          nzm, nzt, ngrdcol, gr, term_wprtpthlp_explicit, rho_ds_zt, rho_ds_zm,             // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_rtpthlp, term_wprtpthlp_explicit_zm, // In
          rhs_ta_wprtpthlp);                                                                // Out

      if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array71",
            ((coef_wpup2_implicit.extent(0) - (1)) / (1) + 1) *
                ((coef_wpup2_implicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((coef_wpup2_implicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((coef_wpup2_implicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((coef_wpup2_implicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((coef_wpup2_implicit.extent(1) - (1)) / (1) + 1);
              coef_wpup2_implicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_coef_wpup2_implicit((1 - 1 + slice_0 * (1)),
                                                               (1 - 1 + slice_1 * (1)));
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array72",
            ((term_wpup2_explicit.extent(0) - (1)) / (1) + 1) *
                ((term_wpup2_explicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((term_wpup2_explicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((term_wpup2_explicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((term_wpup2_explicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((term_wpup2_explicit.extent(1) - (1)) / (1) + 1);
              term_wpup2_explicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_term_wpup2_explicit((1 - 1 + slice_0 * (1)),
                                                               (1 - 1 + slice_1 * (1)));
            });
      }

      // Calculate the momentum level terms and sign of turbulent velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::deep_copy(Kokkos::subview(coef_wpup2_implicit_zm, Kokkos::ALL, Kokkos::ALL),
                          zt2zm_api(nzm, nzt, ngrdcol, gr,
                                    Kokkos::subview(pdf_implicit_coefs_terms.coef_wpup2_implicit,
                                                    Kokkos::ALL, Kokkos::ALL)));
        Kokkos::deep_copy(Kokkos::subview(term_wpup2_explicit_zm, Kokkos::ALL, Kokkos::ALL),
                          zt2zm_api(nzm, nzt, ngrdcol, gr,
                                    Kokkos::subview(pdf_implicit_coefs_terms.term_wpup2_explicit,
                                                    Kokkos::ALL, Kokkos::ALL)));
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array73",
            ((sgn_t_vel_up2.extent(0) - (1)) / (1) + 1) *
                ((sgn_t_vel_up2.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((sgn_t_vel_up2.extent(0) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_up2.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((sgn_t_vel_up2.extent(1) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_up2.extent(1) - (1)) / (1) + 1);
              sgn_t_vel_up2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = Kokkos::copysign(
                  one, (coef_wpup2_implicit_zm((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) *
                            up2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) +
                        term_wpup2_explicit_zm((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)))) *
                           up2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))));
            });
      }

      // Calculate the LHS turbulent advection term for <w'u'^2>
      f2k::xpyp_term_ta_pdf_lhs(
          nzm, nzt, ngrdcol, gr, coef_wpup2_implicit, rho_ds_zt, rho_ds_zm,         // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_up2, coef_wpup2_implicit_zm, // In
          lhs_ta_wpup2);                                                            // Out

      // Calculate the RHS turbulent advection term for <w'u'^2>
      f2k::xpyp_term_ta_pdf_rhs(
          nzm, nzt, ngrdcol, gr, term_wpup2_explicit, rho_ds_zt, rho_ds_zm,         // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_up2, term_wpup2_explicit_zm, // In
          rhs_ta_wpup2);                                                            // Out

      if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array74",
            ((coef_wpvp2_implicit.extent(0) - (1)) / (1) + 1) *
                ((coef_wpvp2_implicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((coef_wpvp2_implicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((coef_wpvp2_implicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((coef_wpvp2_implicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((coef_wpvp2_implicit.extent(1) - (1)) / (1) + 1);
              coef_wpvp2_implicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_coef_wpvp2_implicit((1 - 1 + slice_0 * (1)),
                                                               (1 - 1 + slice_1 * (1)));
            });
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array75",
            ((term_wpvp2_explicit.extent(0) - (1)) / (1) + 1) *
                ((term_wpvp2_explicit.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((term_wpvp2_explicit.extent(0) - (1)) / (1) + 1);
              remaining /= ((term_wpvp2_explicit.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((term_wpvp2_explicit.extent(1) - (1)) / (1) + 1);
              remaining /= ((term_wpvp2_explicit.extent(1) - (1)) / (1) + 1);
              term_wpvp2_explicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                  pdf_implicit_coefs_terms_term_wpvp2_explicit((1 - 1 + slice_0 * (1)),
                                                               (1 - 1 + slice_1 * (1)));
            });
      }

      // Calculate the momentum level terms and sign of turbulent velocity if
      // l_upwind_xpyp_ta is true
      if (l_upwind_xpyp_ta) {
        Kokkos::deep_copy(Kokkos::subview(coef_wpvp2_implicit_zm, Kokkos::ALL, Kokkos::ALL),
                          zt2zm_api(nzm, nzt, ngrdcol, gr,
                                    Kokkos::subview(pdf_implicit_coefs_terms.coef_wpvp2_implicit,
                                                    Kokkos::ALL, Kokkos::ALL)));
        Kokkos::deep_copy(Kokkos::subview(term_wpvp2_explicit_zm, Kokkos::ALL, Kokkos::ALL),
                          zt2zm_api(nzm, nzt, ngrdcol, gr,
                                    Kokkos::subview(pdf_implicit_coefs_terms.term_wpvp2_explicit,
                                                    Kokkos::ALL, Kokkos::ALL)));
        Kokkos::parallel_for(
            "calc_xp2_xpyp_ta_terms::array76",
            ((sgn_t_vel_vp2.extent(0) - (1)) / (1) + 1) *
                ((sgn_t_vel_vp2.extent(1) - (1)) / (1) + 1),
            KOKKOS_LAMBDA(int cell) {
              int remaining     = cell;
              const int slice_0 = remaining % ((sgn_t_vel_vp2.extent(0) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_vp2.extent(0) - (1)) / (1) + 1);
              const int slice_1 = remaining % ((sgn_t_vel_vp2.extent(1) - (1)) / (1) + 1);
              remaining /= ((sgn_t_vel_vp2.extent(1) - (1)) / (1) + 1);
              sgn_t_vel_vp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) = Kokkos::copysign(
                  one, (coef_wpvp2_implicit_zm((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) *
                            vp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) +
                        term_wpvp2_explicit_zm((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)))) *
                           vp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))));
            });
      }

      // Calculate the LHS turbulent advection term for <w'v'^2>
      f2k::xpyp_term_ta_pdf_lhs(
          nzm, nzt, ngrdcol, gr, coef_wpvp2_implicit, rho_ds_zt, rho_ds_zm,         // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_vp2, coef_wpvp2_implicit_zm, // In
          lhs_ta_wpvp2);                                                            // Out

      // Calculate the RHS turbulent advection term for <w'v'^2>
      f2k::xpyp_term_ta_pdf_rhs(
          nzm, nzt, ngrdcol, gr, term_wpvp2_explicit, rho_ds_zt, rho_ds_zm,         // In
          invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_vp2, term_wpvp2_explicit_zm, // In
          rhs_ta_wpvp2);                                                            // Out

      if (l_scalar_calc) {

        for (sclr = 1; sclr <= sclr_dim; sclr += 1) {
          if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::array77",
                ((coef_wpsclrp2_implicit.extent(0) - (1)) / (1) + 1) *
                    ((coef_wpsclrp2_implicit.extent(1) - (1)) / (1) + 1),
                KOKKOS_LAMBDA(int cell) {
                  int remaining = cell;
                  const int slice_0 =
                      remaining % ((coef_wpsclrp2_implicit.extent(0) - (1)) / (1) + 1);
                  remaining /= ((coef_wpsclrp2_implicit.extent(0) - (1)) / (1) + 1);
                  const int slice_1 =
                      remaining % ((coef_wpsclrp2_implicit.extent(1) - (1)) / (1) + 1);
                  remaining /= ((coef_wpsclrp2_implicit.extent(1) - (1)) / (1) + 1);
                  coef_wpsclrp2_implicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                      pdf_implicit_coefs_terms_coef_wpsclrp2_implicit(
                          (1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1));
                });
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::array78",
                ((term_wpsclrp2_explicit.extent(0) - (1)) / (1) + 1) *
                    ((term_wpsclrp2_explicit.extent(1) - (1)) / (1) + 1),
                KOKKOS_LAMBDA(int cell) {
                  int remaining = cell;
                  const int slice_0 =
                      remaining % ((term_wpsclrp2_explicit.extent(0) - (1)) / (1) + 1);
                  remaining /= ((term_wpsclrp2_explicit.extent(0) - (1)) / (1) + 1);
                  const int slice_1 =
                      remaining % ((term_wpsclrp2_explicit.extent(1) - (1)) / (1) + 1);
                  remaining /= ((term_wpsclrp2_explicit.extent(1) - (1)) / (1) + 1);
                  term_wpsclrp2_explicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                      pdf_implicit_coefs_terms_term_wpsclrp2_explicit(
                          (1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1));
                });
          }

          // Calculate the momentum level terms and sign of turbulent velocity if
          // l_upwind_xpyp_ta is true
          if (l_upwind_xpyp_ta) {
            Kokkos::deep_copy(
                Kokkos::subview(coef_wpsclrp2_implicit_zm, Kokkos::ALL, Kokkos::ALL),
                zt2zm_api(nzm, nzt, ngrdcol, gr,
                          Kokkos::subview(pdf_implicit_coefs_terms.coef_wpsclrp2_implicit,
                                          Kokkos::ALL, Kokkos::ALL, (sclr - 1))));
            Kokkos::deep_copy(
                Kokkos::subview(term_wpsclrp2_explicit_zm, Kokkos::ALL, Kokkos::ALL),
                zt2zm_api(nzm, nzt, ngrdcol, gr,
                          Kokkos::subview(pdf_implicit_coefs_terms.term_wpsclrp2_explicit,
                                          Kokkos::ALL, Kokkos::ALL, (sclr - 1))));
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::array79",
                ((sgn_t_vel_sclrp2.extent(0) - (1)) / (1) + 1) *
                    ((sgn_t_vel_sclrp2.extent(1) - (1)) / (1) + 1),
                KOKKOS_LAMBDA(int cell) {
                  int remaining     = cell;
                  const int slice_0 = remaining % ((sgn_t_vel_sclrp2.extent(0) - (1)) / (1) + 1);
                  remaining /= ((sgn_t_vel_sclrp2.extent(0) - (1)) / (1) + 1);
                  const int slice_1 = remaining % ((sgn_t_vel_sclrp2.extent(1) - (1)) / (1) + 1);
                  remaining /= ((sgn_t_vel_sclrp2.extent(1) - (1)) / (1) + 1);
                  sgn_t_vel_sclrp2((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                      Kokkos::copysign(one, (coef_wpsclrp2_implicit_zm((1 - 1 + slice_0 * (1)),
                                                                       (1 - 1 + slice_1 * (1))) *
                                                 sclrp2((1 - 1 + slice_0 * (1)),
                                                        (1 - 1 + slice_1 * (1)), (sclr - 1)) +
                                             term_wpsclrp2_explicit_zm((1 - 1 + slice_0 * (1)),
                                                                       (1 - 1 + slice_1 * (1)))) *
                                                sclrp2((1 - 1 + slice_0 * (1)),
                                                       (1 - 1 + slice_1 * (1)), (sclr - 1)));
                });
          }

          // Calculate the LHS turbulent advection term for <w'sclr'^2>
          f2k::xpyp_term_ta_pdf_lhs(
              nzm, nzt, ngrdcol, gr, coef_wpsclrp2_implicit, rho_ds_zt, rho_ds_zm,            // In
              invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrp2, coef_wpsclrp2_implicit_zm, // In
              Kokkos::subview(lhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL,
                              (sclr - 1))); // Out

          // Calculate the RHS turbulent advection term for <w'sclr'^2>
          f2k::xpyp_term_ta_pdf_rhs(
              nzm, nzt, ngrdcol, gr, term_wpsclrp2_explicit, rho_ds_zt, rho_ds_zm,            // In
              invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrp2, term_wpsclrp2_explicit_zm, // In
              Kokkos::subview(rhs_ta_wpsclrp2, Kokkos::ALL, Kokkos::ALL, (sclr - 1)));        // Out

          if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::array80",
                ((coef_wprtpsclrp_implicit.extent(0) - (1)) / (1) + 1) *
                    ((coef_wprtpsclrp_implicit.extent(1) - (1)) / (1) + 1),
                KOKKOS_LAMBDA(int cell) {
                  int remaining = cell;
                  const int slice_0 =
                      remaining % ((coef_wprtpsclrp_implicit.extent(0) - (1)) / (1) + 1);
                  remaining /= ((coef_wprtpsclrp_implicit.extent(0) - (1)) / (1) + 1);
                  const int slice_1 =
                      remaining % ((coef_wprtpsclrp_implicit.extent(1) - (1)) / (1) + 1);
                  remaining /= ((coef_wprtpsclrp_implicit.extent(1) - (1)) / (1) + 1);
                  coef_wprtpsclrp_implicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                      pdf_implicit_coefs_terms_coef_wprtpsclrp_implicit(
                          (1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1));
                });
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::array81",
                ((term_wprtpsclrp_explicit.extent(0) - (1)) / (1) + 1) *
                    ((term_wprtpsclrp_explicit.extent(1) - (1)) / (1) + 1),
                KOKKOS_LAMBDA(int cell) {
                  int remaining = cell;
                  const int slice_0 =
                      remaining % ((term_wprtpsclrp_explicit.extent(0) - (1)) / (1) + 1);
                  remaining /= ((term_wprtpsclrp_explicit.extent(0) - (1)) / (1) + 1);
                  const int slice_1 =
                      remaining % ((term_wprtpsclrp_explicit.extent(1) - (1)) / (1) + 1);
                  remaining /= ((term_wprtpsclrp_explicit.extent(1) - (1)) / (1) + 1);
                  term_wprtpsclrp_explicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                      pdf_implicit_coefs_terms_term_wprtpsclrp_explicit(
                          (1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1));
                });
          }

          // Calculate the momentum level terms and sign of turbulent velocity if
          // l_upwind_xpyp_ta is true
          if (l_upwind_xpyp_ta) {
            Kokkos::deep_copy(
                Kokkos::subview(coef_wprtpsclrp_implicit_zm, Kokkos::ALL, Kokkos::ALL),
                zt2zm_api(nzm, nzt, ngrdcol, gr,
                          Kokkos::subview(pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit,
                                          Kokkos::ALL, Kokkos::ALL, (sclr - 1))));
            Kokkos::deep_copy(
                Kokkos::subview(term_wprtpsclrp_explicit_zm, Kokkos::ALL, Kokkos::ALL),
                zt2zm_api(nzm, nzt, ngrdcol, gr,
                          Kokkos::subview(pdf_implicit_coefs_terms.term_wprtpsclrp_explicit,
                                          Kokkos::ALL, Kokkos::ALL, (sclr - 1))));
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::array82",
                ((sgn_t_vel_sclrprtp.extent(0) - (1)) / (1) + 1) *
                    ((sgn_t_vel_sclrprtp.extent(1) - (1)) / (1) + 1),
                KOKKOS_LAMBDA(int cell) {
                  int remaining     = cell;
                  const int slice_0 = remaining % ((sgn_t_vel_sclrprtp.extent(0) - (1)) / (1) + 1);
                  remaining /= ((sgn_t_vel_sclrprtp.extent(0) - (1)) / (1) + 1);
                  const int slice_1 = remaining % ((sgn_t_vel_sclrprtp.extent(1) - (1)) / (1) + 1);
                  remaining /= ((sgn_t_vel_sclrprtp.extent(1) - (1)) / (1) + 1);
                  sgn_t_vel_sclrprtp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                      Kokkos::copysign(one, (coef_wprtpsclrp_implicit_zm((1 - 1 + slice_0 * (1)),
                                                                         (1 - 1 + slice_1 * (1))) *
                                                 sclrprtp((1 - 1 + slice_0 * (1)),
                                                          (1 - 1 + slice_1 * (1)), (sclr - 1)) +
                                             term_wprtpsclrp_explicit_zm((1 - 1 + slice_0 * (1)),
                                                                         (1 - 1 + slice_1 * (1)))) *
                                                sclrprtp((1 - 1 + slice_0 * (1)),
                                                         (1 - 1 + slice_1 * (1)), (sclr - 1)));
                });
          }

          // Calculate the LHS turbulent advection term for <w'rt'sclr'>
          f2k::xpyp_term_ta_pdf_lhs(nzm, nzt, ngrdcol, gr, coef_wprtpsclrp_implicit, rho_ds_zt,
                                    rho_ds_zm,                                             // In
                                    invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrprtp, // In
                                    coef_wprtpsclrp_implicit_zm,                           // In
                                    Kokkos::subview(lhs_ta_wprtpsclrp, Kokkos::ALL, Kokkos::ALL,
                                                    Kokkos::ALL, (sclr - 1))); // Out

          // Calculate the RHS turbulent advection term for <w'rt'sclr'>
          f2k::xpyp_term_ta_pdf_rhs(
              nzm, nzt, ngrdcol, gr, term_wprtpsclrp_explicit, rho_ds_zt, rho_ds_zm,     // In
              invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrprtp,                     // In
              term_wprtpsclrp_explicit_zm,                                               // In
              Kokkos::subview(rhs_ta_wprtpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out

          if (!l_upwind_xpyp_ta || stats.flags().l_sample) {
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::array83",
                ((coef_wpthlpsclrp_implicit.extent(0) - (1)) / (1) + 1) *
                    ((coef_wpthlpsclrp_implicit.extent(1) - (1)) / (1) + 1),
                KOKKOS_LAMBDA(int cell) {
                  int remaining = cell;
                  const int slice_0 =
                      remaining % ((coef_wpthlpsclrp_implicit.extent(0) - (1)) / (1) + 1);
                  remaining /= ((coef_wpthlpsclrp_implicit.extent(0) - (1)) / (1) + 1);
                  const int slice_1 =
                      remaining % ((coef_wpthlpsclrp_implicit.extent(1) - (1)) / (1) + 1);
                  remaining /= ((coef_wpthlpsclrp_implicit.extent(1) - (1)) / (1) + 1);
                  coef_wpthlpsclrp_implicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                      pdf_implicit_coefs_terms_coef_wpthlpsclrp_implicit(
                          (1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1));
                });
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::array84",
                ((term_wpthlpsclrp_explicit.extent(0) - (1)) / (1) + 1) *
                    ((term_wpthlpsclrp_explicit.extent(1) - (1)) / (1) + 1),
                KOKKOS_LAMBDA(int cell) {
                  int remaining = cell;
                  const int slice_0 =
                      remaining % ((term_wpthlpsclrp_explicit.extent(0) - (1)) / (1) + 1);
                  remaining /= ((term_wpthlpsclrp_explicit.extent(0) - (1)) / (1) + 1);
                  const int slice_1 =
                      remaining % ((term_wpthlpsclrp_explicit.extent(1) - (1)) / (1) + 1);
                  remaining /= ((term_wpthlpsclrp_explicit.extent(1) - (1)) / (1) + 1);
                  term_wpthlpsclrp_explicit((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                      pdf_implicit_coefs_terms_term_wpthlpsclrp_explicit(
                          (1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1)), (sclr - 1));
                });
          }

          // Calculate the momentum level terms and sign of turbulent velocity if
          // l_upwind_xpyp_ta is true
          if (l_upwind_xpyp_ta) {
            Kokkos::deep_copy(
                Kokkos::subview(coef_wpthlpsclrp_implicit_zm, Kokkos::ALL, Kokkos::ALL),
                zt2zm_api(nzm, nzt, ngrdcol, gr,
                          Kokkos::subview(pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit,
                                          Kokkos::ALL, Kokkos::ALL, (sclr - 1))));
            Kokkos::deep_copy(
                Kokkos::subview(term_wpthlpsclrp_explicit_zm, Kokkos::ALL, Kokkos::ALL),
                zt2zm_api(nzm, nzt, ngrdcol, gr,
                          Kokkos::subview(pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit,
                                          Kokkos::ALL, Kokkos::ALL, (sclr - 1))));
            Kokkos::parallel_for(
                "calc_xp2_xpyp_ta_terms::array85",
                ((sgn_t_vel_sclrpthlp.extent(0) - (1)) / (1) + 1) *
                    ((sgn_t_vel_sclrpthlp.extent(1) - (1)) / (1) + 1),
                KOKKOS_LAMBDA(int cell) {
                  int remaining     = cell;
                  const int slice_0 = remaining % ((sgn_t_vel_sclrpthlp.extent(0) - (1)) / (1) + 1);
                  remaining /= ((sgn_t_vel_sclrpthlp.extent(0) - (1)) / (1) + 1);
                  const int slice_1 = remaining % ((sgn_t_vel_sclrpthlp.extent(1) - (1)) / (1) + 1);
                  remaining /= ((sgn_t_vel_sclrpthlp.extent(1) - (1)) / (1) + 1);
                  sgn_t_vel_sclrpthlp((1 - 1 + slice_0 * (1)), (1 - 1 + slice_1 * (1))) =
                      Kokkos::copysign(one,
                                       (coef_wpthlpsclrp_implicit_zm((1 - 1 + slice_0 * (1)),
                                                                     (1 - 1 + slice_1 * (1))) *
                                            sclrpthlp((1 - 1 + slice_0 * (1)),
                                                      (1 - 1 + slice_1 * (1)), (sclr - 1)) +
                                        term_wpthlpsclrp_explicit_zm((1 - 1 + slice_0 * (1)),
                                                                     (1 - 1 + slice_1 * (1)))) *
                                           sclrpthlp((1 - 1 + slice_0 * (1)),
                                                     (1 - 1 + slice_1 * (1)), (sclr - 1)));
                });
          }

          // Calculate the LHS turbulent advection term for <w'thl'sclr'>
          f2k::xpyp_term_ta_pdf_lhs(nzm, nzt, ngrdcol, gr, coef_wpthlpsclrp_implicit, rho_ds_zt,
                                    rho_ds_zm,                                              // In
                                    invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrpthlp, // In
                                    coef_wpthlpsclrp_implicit_zm,                           // In
                                    Kokkos::subview(lhs_ta_wpthlpsclrp, Kokkos::ALL, Kokkos::ALL,
                                                    Kokkos::ALL, (sclr - 1))); // Out

          // Calculate the RHS turbulent advection term for <w'thl'sclr'>
          f2k::xpyp_term_ta_pdf_rhs(
              nzm, nzt, ngrdcol, gr, term_wpthlpsclrp_explicit, rho_ds_zt, rho_ds_zm,     // In
              invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_sclrpthlp,                     // In
              term_wpthlpsclrp_explicit_zm,                                               // In
              Kokkos::subview(rhs_ta_wpthlpsclrp, Kokkos::ALL, Kokkos::ALL, (sclr - 1))); // Out
        }
      }
    }
  }

  if (iiPDF_type != iiPDF_new_hybrid) {

    // Set up the implicit coefficients and explicit terms for turbulent
    // advection of <u'^2> and <v'^2> following ADG1, then calculate the
    // turbulent advection terms.

    // CLUBB does not produce a PDF for horizontal wind components u and v.
    // However, the code for the ADG1 PDF is still used to handle the
    // turbulent advection for the variances of the horizontal wind
    // components, <u'^2> and <v'^2>.  The ADG1 code is used regardless of
    // which PDF type or turbulent advection option is used.  The implicit
    // coefficients and explicit terms are calculated on thermodynamic grid
    // levels.

    if (l_upwind_xpyp_ta) {
      // Calculate coef_wpup2_wpvp2_implicit, term_wpup2_explicit, and
      // term_wpvp2_explicit on momentum levels as
      // coef_wpup2_wpvp2_implicit_zm(i,:), term_wpup2_explicit_zm, and
      // term_wpvp2_explicit_zm, respectively.
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::loop86", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            coef_wpup2_implicit_zm(i - 1, k - 1) =
                one_third * beta(i - 1) * a1_coef(i - 1, k - 1) * wp3_on_wp2(i - 1, k - 1);
            coef_wpvp2_implicit_zm(i - 1, k - 1) = coef_wpup2_implicit_zm(i - 1, k - 1);
            term_wpup2_explicit_zm(i - 1, k - 1) =
                wp_coef(i - 1, k - 1) * (upwp(i - 1, k - 1) * upwp(i - 1, k - 1));
            term_wpvp2_explicit_zm(i - 1, k - 1) =
                wp_coef(i - 1, k - 1) * (vpwp(i - 1, k - 1) * vpwp(i - 1, k - 1));
          });

      // For ADG1, the sign of the turbulent velocity is the sign of
      // <w'^3> / <w'^2>.  For simplicity, the sign of turbulent velocity is
      // set to wp3_on_wp2.
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::loop87", std::max(0, nzm) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzm);
            remaining /= (nzm);
            sgn_t_vel_up2(i - 1, k - 1) = Kokkos::copysign(one, wp3_on_wp2(i - 1, k - 1));
            sgn_t_vel_vp2(i - 1, k - 1) = sgn_t_vel_up2(i - 1, k - 1);
          });
    } else {

      // Interpolate <u'w'> and <v'w'> from the momentum levels to the
      // thermodynamic levels.  These will be used for the turbulent
      // advection terms in each equation.
      Kokkos::deep_copy(
          Kokkos::subview(upwp_zt, Kokkos::ALL, Kokkos::ALL),
          zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(upwp, Kokkos::ALL, Kokkos::ALL)));
      Kokkos::deep_copy(
          Kokkos::subview(vpwp_zt, Kokkos::ALL, Kokkos::ALL),
          zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(vpwp, Kokkos::ALL, Kokkos::ALL)));

      // Implicit coefficient on <u'^2> or <v'^2> in <w'u'^2> or <w'v'^2>
      // equation.
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::loop88", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzt);
            remaining /= (nzt);
            coef_wpup2_implicit(i - 1, k - 1) =
                one_third * beta(i - 1) * a1_coef_zt(i - 1, k - 1) * wp3_on_wp2_zt(i - 1, k - 1);
            coef_wpvp2_implicit(i - 1, k - 1) = coef_wpup2_implicit(i - 1, k - 1);
          });

      // Explicit (RHS) term in <w'u'^2> equation.
      Kokkos::parallel_for(
          "calc_xp2_xpyp_ta_terms::loop89", std::max(0, nzt) * std::max(0, ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            int remaining = cell;
            const int i   = 1 + remaining % (ngrdcol);
            remaining /= (ngrdcol);
            const int k = 1 + remaining % (nzt);
            remaining /= (nzt);
            term_wpup2_explicit(i - 1, k - 1) =
                wp_coef_zt(i - 1, k - 1) * (upwp_zt(i - 1, k - 1) * upwp_zt(i - 1, k - 1));

            // Explicit (RHS) term in <w'v'^2> equation.
            term_wpvp2_explicit(i - 1, k - 1) =
                wp_coef_zt(i - 1, k - 1) * (vpwp_zt(i - 1, k - 1) * vpwp_zt(i - 1, k - 1));
          });
    }

    // Calculate the LHS turbulent advection term for <w'u'^2> and <w'v'^2>
    f2k::xpyp_term_ta_pdf_lhs(
        nzm, nzt, ngrdcol, gr, coef_wpup2_implicit, rho_ds_zt, rho_ds_zm,         // In
        invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_up2, coef_wpup2_implicit_zm, // In
        lhs_ta_wpup2);                                                            // Out

    // The same LHS is used for <w'u'^2> and <w'v'^2>
    Kokkos::parallel_for(
        "calc_xp2_xpyp_ta_terms::loop90",
        std::max(0, nzm) * std::max(0, ngrdcol) * std::max(0, ndiags3),
        KOKKOS_LAMBDA(const int cell) {
          int remaining = cell;
          const int b   = 1 + remaining % (ndiags3);
          remaining /= (ndiags3);
          const int i = 1 + remaining % (ngrdcol);
          remaining /= (ngrdcol);
          const int k = 1 + remaining % (nzm);
          remaining /= (nzm);
          lhs_ta_wpvp2(b - 1, i - 1, k - 1) = lhs_ta_wpup2(b - 1, i - 1, k - 1);
        });

    // Calculate the RHS turbulent advection term for <w'u'^2>
    f2k::xpyp_term_ta_pdf_rhs(
        nzm, nzt, ngrdcol, gr, term_wpup2_explicit, rho_ds_zt, rho_ds_zm,         // In
        invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_up2, term_wpup2_explicit_zm, // In
        rhs_ta_wpup2);                                                            // Out

    // Calculate the RHS turbulent advection term for <w'v'^2>
    f2k::xpyp_term_ta_pdf_rhs(
        nzm, nzt, ngrdcol, gr, term_wpvp2_explicit, rho_ds_zt, rho_ds_zm,         // In
        invrs_rho_ds_zm, l_upwind_xpyp_ta, sgn_t_vel_vp2, term_wpvp2_explicit_zm, // In
        rhs_ta_wpvp2);                                                            // Out
  }

  // Stats output for implicit coefficients and explicit terms of
  // <w'rt'^2>, <w'thl'^2>, and <w'rt'thl'> used in the calcualtion of
  // the turbulent advection terms

  if (stats.flags().l_sample) {
    f2k::stats_update("coef_wprtp2_implicit", coef_wprtp2_implicit,       // In
                      stats);                                             // InOut
    f2k::stats_update("term_wprtp2_explicit", term_wprtp2_explicit,       // In
                      stats);                                             // InOut
    f2k::stats_update("coef_wpthlp2_implicit", coef_wpthlp2_implicit,     // In
                      stats);                                             // InOut
    f2k::stats_update("term_wpthlp2_explicit", term_wpthlp2_explicit,     // In
                      stats);                                             // InOut
    f2k::stats_update("coef_wprtpthlp_implicit", coef_wprtpthlp_implicit, // In
                      stats);                                             // InOut
    f2k::stats_update("term_wprtpthlp_explicit", term_wprtpthlp_explicit, // In
                      stats);                                             // InOut
  }

  return;
}

// xam: x_am(k)                     [{x_am units}]
// xbm: x_bm(k)                     [{x_bm units}]
// wpxbp: w'x_b'(k)                   [m/s {x_bm units}]
// wpxap: w'x_a'(k)                   [m/s {x_am units}]
// invrs_dzm: Inverse of grid spacing (k) [1/m]
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5>
void term_tp_rhs(int nzm, int nzt, int ngrdcol, Array0 xam, Array1 xbm, // In
                 Array2 wpxbp, Array3 wpxap, Array4 invrs_dzm,          // In
                 Array5 rhs)                                            // Out
{
  // Description:
  // Turbulent production of x_a'x_b':  explicit portion of the code.
  //
  // The d(x_a'x_b')/dt equation contains a turbulent production term:
  //
  // - w'x_b' d(x_am)/dz - w'x_a' d(x_bm)/dz.
  //
  // This term is solved for completely explicitly and is discretized as
  // follows:
  //
  // The values of w'x_a' and w'x_b' are found on the momentum levels, whereas
  // the values of x_am and x_bm are found on the thermodynamic levels.  The
  // derivatives of both x_am and x_bm are taken over the intermediate
  // (central) momentum level.  All of the remaining mathematical operations
  // take place at the central momentum level, yielding the desired result.
  //
  // ---------xam--------------xbm---------------------------- t(k)
  //
  // ===wpxap======d(xam)/dz=========d(xbm)/dz===wpxbp======== m(k)
  //
  // ---------xam--------------xbm---------------------------- t(k-1)
  //
  // The vertical indices t(k), m(k), and t(k-1) correspond with altitudes
  // zt(k), zm(k), and zt(k-1), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )
  // References:
  //-----------------------------------------------------------------------
  //------------------------ Input variables ------------------------
  //------------------------ Return Variable ------------------------
  //------------------------ Local Variables ------------------------
  int i;
  int k;

  //------------------------ Begin Code ------------------------

  Kokkos::parallel_for(
      "term_tp_rhs::loop1", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        rhs(i - 1, k - 1) = -wpxbp(i - 1, k - 1) * invrs_dzm(i - 1, k - 1) *
                                (xam(i - 1, k - 1) - xam(i - 1, k - 1 - 1)) -
                            wpxap(i - 1, k - 1) * invrs_dzm(i - 1, k - 1) *
                                (xbm(i - 1, k - 1) - xbm(i - 1, k - 1 - 1));
      });

  return;
}

// Cn: Coefficient C_n                       [-]
// invrs_tau_zm: Inverse time-scale tau at momentum levels (k) [1/s]
template <class Array0, class Array1, class Array2>
void term_dp1_lhs(int nzm, int ngrdcol, const grid& gr, // In
                  Array0 Cn, Array1 invrs_tau_zm,       // In
                  Array2 lhs)                           // Out
{
  // Description:
  // Dissipation term 1 for x_a'x_b':  implicit portion of the code.
  //
  // The d(x_a'x_b')/dt equation contains dissipation term 1:
  //
  // - ( C_n / tau_zm ) x_a'x_b'.
  //
  // For cases where x_a'x_b' is a variance (in other words, where x_a and x_b
  // are the same variable), the term is damped to a certain positive
  // threshold, such that:
  //
  // - ( C_n / tau_zm ) * ( x_a'x_b' - threshold ).
  //
  // However, if x_a'x_b' is u'^2 or v'^2, damping to a minimum threshold value
  // is part of pressure term 1 and is handled as part of function 'term_pr1'.
  // Thus, for u'^2 and v'^2, function 'term_dp1_lhs' is called, but function
  // 'term_dp1_rhs' is not called, as function 'term_pr1' is called instead.
  //
  // For cases where x_a'x_b' is a covariance (in other words, where x_a and
  // x_b are different variables), threshold is set to 0, and the expression
  // reverts to the form found in the first equation.
  //
  // This term is broken into implicit and explicit portions.  The equations
  // for u'^2, v'^2, and any covariances only include the implicit portion.
  // The implicit portion of this term is:
  //
  // - ( C_n / tau_zm ) x_a'x_b'(t+1).
  //
  // Note:  When the implicit term is brought over to the left-hand side,
  //        the sign is reversed and the leading "-" in front of the term
  //        is changed to a "+".
  //
  // The timestep index (t+1) means that the value of x_a'x_b' being used is
  // from the next timestep, which is being advanced to in solving the
  // d(x_a'x_b')/dt equation.
  //
  // The values of x_a'x_b' are found on momentum levels.  The values of
  // time-scale tau_zm are also found on momentum levels.
  //
  // Note:  For equations that use pressure term 1 (such as the equations for
  //        u'^2 and v'^2), C_n = ( 2*C_4 + C_14 ) / 3; which combines the
  //        implicit contributions for dissipation term 1 and pressure term 1
  //        into one expression.  Otherwise, C_n = C_2.
  // References:
  //-----------------------------------------------------------------------
  //--------------------------- Input Variables ---------------------------
  //--------------------------- Return Variables ---------------------------
  //--------------------------- Local Variables ---------------------------
  int i;
  int k;
  const auto gr_k_lb_zm = gr.k_lb_zm;
  const auto gr_k_ub_zm = gr.k_ub_zm;

  //--------------------------- Begin Code ---------------------------

  Kokkos::parallel_for(
      "term_dp1_lhs::loop1", std::max(0, ngrdcol), KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        // Set lower boundary to 0
        lhs(i - 1, gr_k_lb_zm - 1) = zero;
        // Set upper boundary to 0
        lhs(i - 1, gr_k_ub_zm - 1) = zero;
      });

  Kokkos::parallel_for(
      "term_dp1_lhs::loop2", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        // Momentum main diagonal: [ x xapxbp(k,<t+1>) ]
        lhs(i - 1, k - 1) = +Cn(i - 1, k - 1) * invrs_tau_zm(i - 1, k - 1);
      });

  return;
}

// Cn: Coefficient C_n                               [-]
// invrs_tau_zm: Time-scale tau at momentum levels (k)         [1/s]
// threshold: Minimum allowable magnitude value of x_a'x_b' [units vary]
template <class Array0, class Array1, class Array2>
void term_dp1_rhs(int nzm, int ngrdcol,                             // In
                  Array0 Cn, Array1 invrs_tau_zm, real_t threshold, // In
                  Array2 rhs)                                       // Out
{
  // Description:
  // Dissipation term 1 for x_a'x_b':  explicit portion of the code.
  //
  // The d(x_a'x_b')/dt equation contains dissipation term 1:
  //
  // - ( C_n / tau_zm ) x_a'x_b'.
  //
  // For cases where x_a'x_b' is a variance (in other words, where x_a and x_b
  // are the same variable), the term is damped to a certain positive
  // threshold, such that:
  //
  // - ( C_n / tau_zm ) * ( x_a'x_b' - threshold ).
  //
  // However, if x_a'x_b' is u'^2 or v'^2, damping to a minimum threshold value
  // is part of pressure term 1 and is handled as part of function 'term_pr1'.
  // Thus, for u'^2 and v'^2, function 'term_dp1_lhs' is called, but function
  // 'term_dp1_rhs' is not called, as function 'term_pr1' is called instead.
  //
  // For cases where x_a'x_b' is a covariance (in other words, where x_a and
  // x_b are different variables), threshold is set to 0, and the expression
  // reverts to the form found in the first equation.
  //
  // This term is broken into implicit and explicit portions.  The equations
  // for u'^2, v'^2, and any covariances only include the implicit portion.
  // The explicit portion of this term is:
  //
  // + ( C_n / tau_zm ) * threshold.
  //
  // The values of time-scale tau_zm and the threshold are found on the
  // momentum levels.
  //
  // Note:  The equations that use pressure term 1 (such as the equations for
  //        u'^2 and v'^2) do not call this function.  Thus, within this
  //        function, C_n = C_2.
  // References:
  //-----------------------------------------------------------------------
  //--------------------------- Input Variables ---------------------------
  //--------------------------- Return Variable ---------------------------
  //--------------------------- Local Variables ---------------------------
  int i;
  int k;

  //--------------------------- Begin Code ---------------------------

  Kokkos::parallel_for(
      "term_dp1_rhs::loop1", std::max(0, nzm) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzm);
        remaining /= (nzm);
        rhs(i - 1, k - 1) = +Cn(i - 1, k - 1) * invrs_tau_zm(i - 1, k - 1) * threshold;
      });

  return;
}

// C4: Model parameter C_4                                 [-]
// C14: Model parameter C_14                                [-]
// xbp2: v'^2(k) (if solving for u'^2) or vice versa         [m^2/s^2]
// wp2: w'^2(k)                                             [m^2/s^2]
// invrs_tau_C4_zm: Time-scale tau for C4 terms at momentum levels (k)  [1/s]
// invrs_tau_C14_zm: Time-scale tau for C14 terms at momentum levels (k) [1/s]
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6>
void term_pr1(int nzm, int ngrdcol, Array0 C4, Array1 C14, Array2 xbp2,    // In
              Array3 wp2, Array4 invrs_tau_C4_zm, Array5 invrs_tau_C14_zm, // In
              Array6 rhs)                                                  // Out
{
  // Description:
  // Pressure term 1 for x_a'x_b':  explicit portion of the code.
  //
  // Note:  Pressure term 1 is only used when x_a'x_b' is either u'^2 or v'^2.
  //        For the following description, pressure term 2 for u'^2 is used as
  //        the example.  Pressure term 2 for v'^2 is the same as pressure
  //        term 2 for u'^2, except that the v'^2 and u'^2 variables are
  //        switched.
  //
  // The d(u'^2)/dt equation contains dissipation term 1:
  //
  // - ( C_4 / tau_zm ) * ( u'^2 - (2/3)*em );
  //
  // where em = (1/2) * ( u'^2 + v'^2 + w'^2 );
  //
  // and with the substitution applied, dissipation term 1 becomes:
  //
  // - ( C_4 / tau_zm ) * ( u'^2 - (1/3) * ( u'^2 + v'^2 + w'^2 ) ).
  //
  // The d(u'^2)/dt equation also contains pressure term 1:
  //
  // - (2/3) * epsilon;
  //
  // where epsilon = C_14 * ( em / tau_zm ).
  //
  // Additionally, since pressure term 1 is a damping term, em is damped only
  // to it's minimum threshold value, em_min, where:
  //
  // em_min = (1/2) * ( u'^2|_min + v'^2|_min + w'^2|_min )
  //      = (1/2) * ( w_tol^2 + w_tol^2 + w_tol^2 )
  //      = (3/2) * w_tol^2.
  //
  // With the damping threshold applied, epsilon becomes:
  //
  // epsilon = C_14 * ( ( em - em_min ) / tau_zm );
  //
  // and with all substitutions applied, pressure term 1 becomes:
  //
  // - (2/3) * ( C_14 / tau_zm )
  //         * [ (1/2) * ( u'^2 + v'^2 + w'^2 ) - (3/2) * w_tol^2 ].
  //
  // Dissipation term 1 and pressure term 1 are combined and simplify to:
  //
  // - [ ( 2*C_4 + C_14 ) / ( 3 * tau_zm ) ] * u'^2
  //    + [ ( C_4 - C_14 ) / ( 3 * tau_zm ) ] * ( v'^2 + w'^2 )
  //    + ( C_14 / tau_zm ) * w_tol^2.
  //
  // The combined term has both implicit and explicit components.
  // The implicit component is:
  //
  // - [ ( 2*C_4 + C_14 ) / ( 3 * tau_zm ) ] * u'^2(t+1).
  //
  // Note:  When the implicit term is brought over to the left-hand side,
  //        the sign is reversed and the leading "-" in front of the term
  //        is changed to a "+".
  //
  // Timestep index (t) stands for the index of the current timestep, while
  // timestep index (t+1) stands for the index of the next timestep, which is
  // being advanced to in solving the d(x_a'x_b')/dt equation.
  //
  // The implicit component of the combined dp1 and pr1 term is solved in
  // function "term_dp1_lhs" above, where "( 2*C_4 + C_14 ) / 3" is sent in
  // as "C_n".
  //
  // The explicit component of the combined dp1 and pr1 term is:
  //
  // + [ ( C_4 - C_14 ) / ( 3 * tau_zm ) ] * ( v'^2(t) + w'^2(t) )
  // + ( C_14 / tau_zm ) * w_tol^2;
  //
  // and is discretized as follows:
  //
  // The values for v'^2 and w'^2, as well as for tau_zm, are found on the
  // momentum levels.  The mathematical operations all take place on the
  // momentum levels, yielding the desired result.
  // References:
  //-----------------------------------------------------------------------
  //------------------------- Input Variables -------------------------
  // Return Variable
  //------------------------- Return Variables -------------------------
  //------------------------- Local Variables -------------------------
  int i;
  int k;

  //------------------------- Begin Code -------------------------

  Kokkos::parallel_for(
      "term_pr1::loop1", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        rhs(i - 1, k - 1) = +one_third * C4(i - 1) * (xbp2(i - 1, k - 1) + wp2(i - 1, k - 1)) *
                                invrs_tau_C4_zm(i - 1, k - 1) -
                            one_third * C14(i - 1) * (xbp2(i - 1, k - 1) + wp2(i - 1, k - 1)) *
                                invrs_tau_C14_zm(i - 1, k - 1) +
                            C14(i - 1) * invrs_tau_C14_zm(i - 1, k - 1) * w_tol_sqd;
      });

  return;
}

// C_uu_shr: Model parameter C_uu_shr                       [-]
// C_uu_buoy: Model parameter C_uu_buoy                      [-]
// thv_ds_zm: Dry, base-state theta_v at momentum level (k)  [K]
// wpthvp: w'th_v'(k)                                     [m/K/s]
// upwp: u'w'(k)                                        [m^2/s^2]
// vpwp: v'w'(k)                                        [m^2/s^2]
// um: mean zonal wind       [m/s]
// vm: mean meridional wind  [m/s]
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8>
void term_pr2(int nzm, int nzt, int ngrdcol, const grid& gr,                      // In
              Array0 C_uu_shr, Array1 C_uu_buoy, Array2 thv_ds_zm, Array3 wpthvp, // In
              Array4 upwp,                                                        // In
              Array5 vpwp, Array6 um, Array7 vm,                                  // In
              Array8 rhs_pr2)                                                     // Out
{
  // Description:
  // Pressure term 2 for x_a'x_b':  explicit portion of the code.
  //
  // Note:  Pressure term 2 is only used when x_a'x_b' is either u'^2 or v'^2.
  //        For the following description, pressure term 2 for u'^2 is used as
  //        the example.  Pressure term 2 for v'^2 is the exact same as
  //        pressure term 2 for u'^2.
  //
  // The d(u'^2)/dt equation contains pressure term 2:
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
  // References:
  //-----------------------------------------------------------------------
  // External
  // ------------ Input Variables ------------
  // Note: Entire arrays of um and vm are now required rather than um and vm
  // only at levels k and k+1.  The entire array is necessary when a vertical
  // average calculation of d(um)/dz and d(vm)/dz is used. --ldgrant March 2010
  // ------------ Output Variable ------------
  // ------------ Local Variables ------------
  int i;
  int k;
  const auto gr_invrs_dzm = gr.invrs_dzm;

  // calculation for d(um)/dz and d(vm)/dz

  // ------------ Begin code ------------

  // use original version of term_pr2

  // As applied to w'2
  Kokkos::parallel_for(
      "term_pr2::loop1", std::max(0, nzm - 2) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 2 + remaining % (nzm - 2);
        remaining /= (nzm - 2);
        rhs_pr2(i - 1, k - 1) =
            +two_thirds *
            (C_uu_buoy(i - 1) * (grav / thv_ds_zm(i - 1, k - 1)) * wpthvp(i - 1, k - 1) +
             C_uu_shr(i - 1) * (-upwp(i - 1, k - 1) * gr_invrs_dzm(i - 1, k - 1) *
                                    (um(i - 1, k - 1) - um(i - 1, k - 1 - 1)) -
                                vpwp(i - 1, k - 1) * gr_invrs_dzm(i - 1, k - 1) *
                                    (vm(i - 1, k - 1) - vm(i - 1, k - 1 - 1))));

        // Added by dschanen for ticket #36
        // We have found that when shear generation is zero this term will only be
        // offset by hole-filling (up2_pd/vp2_pd) and reduces turbulence
        // unrealistically at lower altitudes to make up the difference.
        rhs_pr2(i - 1, k - 1) =
            ((rhs_pr2(i - 1, k - 1)) > (zero_threshold) ? (rhs_pr2(i - 1, k - 1))
                                                        : (zero_threshold));
      });

  return;
}

// dt: Model timestep              [s]
// tolerance: Threshold for xp2_np1       [units vary]
// rho_ds_zm: Dry, static density on momentum levels         [kg/m^3]
// xp2_np1: Variance for <n+1>          [units vary]
template <class Array0, class Array1>
void pos_definite_variances(int nzm, int ngrdcol, const grid& gr, // In
                            int solve_type, int fill_holes_type,  // In
                            real_t dt, real_t tolerance,          // In
                            Array0 rho_ds_zm,                     // In
                            f2k::stats_type& stats,               // InOut
                            Array1 xp2_np1)                       // InOut
{
  // Description:
  // Use the hole filling code to make a variance term positive definite
  //-----------------------------------------------------------------------
  //------------------- Input variables ------------------
  //------------------- Local variables ------------------
  std::string name_pd;
  std::string name_hf_before;
  std::string name_hf_after;

  //------------------- Begin Code ------------------

  switch (solve_type) {
  case xp2_xpyp_rtp2: {
    name_pd        = "rtp2_pd";
    name_hf_before = "rtp2_before_hf";
    name_hf_after  = "rtp2_after_hf";
    break;
  }
  case xp2_xpyp_thlp2: {
    name_pd        = "thlp2_pd";
    name_hf_before = "thlp2_before_hf";
    name_hf_after  = "thlp2_after_hf";
    break;
  }
  case xp2_xpyp_up2: {
    name_pd        = "up2_pd";
    name_hf_before = "up2_before_hf";
    name_hf_after  = "up2_after_hf";
    break;
  }
  case xp2_xpyp_vp2: {
    name_pd        = "vp2_pd";
    name_hf_before = "vp2_before_hf";
    name_hf_after  = "vp2_after_hf";
    break;
  }
  default: {
    name_pd        = "";
    name_hf_before = "";
    name_hf_after  = "";
    break;
  }
  }

  if (stats.flags().l_sample) {
    // Store previous value for effect of the positive definite scheme
    if (name_pd.size() > 0) {
      {
        Real2D stats_values_1(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                                 "pos_definite_variances::stats_values_1"),
                              xp2_np1.extent(0), xp2_np1.extent(1));
        Kokkos::parallel_for(
            "pos_definite_variances::stats_value", xp2_np1.size(), KOKKOS_LAMBDA(int cell) {
              const int i          = cell % xp2_np1.extent(0);
              const int k          = cell / xp2_np1.extent(0);
              stats_values_1(i, k) = xp2_np1(i, k) / dt;
            });
        f2k::stats_begin_budget(name_pd, stats_values_1, stats);
      }
    }
  }

  // Call the hole-filling scheme.
  // The first pass-through should draw from only two levels on either side
  // of the hole.
  // upper_hf_level = nz-1 since we are filling the zm levels
  if (fill_holes_type != 0) {
    if (stats.flags().l_sample && name_hf_before.size() > 0) {
      f2k::stats_update(name_hf_before, xp2_np1, // In
                        stats);                  // InOut
    }

    f2k::fill_holes_vertical_api(nzm, ngrdcol, tolerance, gr.k_lb_zm + gr.grid_dir_indx, // In
                                 gr.k_ub_zm - gr.grid_dir_indx, gr.dzm, rho_ds_zm,
                                 gr.grid_dir_indx, // In
                                 fill_holes_type,  // In
                                 xp2_np1);         // InOut

    if (stats.flags().l_sample && name_hf_after.size() > 0) {
      f2k::stats_update(name_hf_after, xp2_np1, // In
                        stats);                 // InOut
    }
  }

  if (stats.flags().l_sample) {
    if (name_pd.size() > 0) {
      {
        Real2D stats_values_2(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                                 "pos_definite_variances::stats_values_2"),
                              xp2_np1.extent(0), xp2_np1.extent(1));
        Kokkos::parallel_for(
            "pos_definite_variances::stats_value", xp2_np1.size(), KOKKOS_LAMBDA(int cell) {
              const int i          = cell % xp2_np1.extent(0);
              const int k          = cell / xp2_np1.extent(0);
              stats_values_2(i, k) = xp2_np1(i, k) / dt;
            });
        f2k::stats_finalize_budget(name_pd, stats_values_2, stats);
      }
    }
  }

  return;
}

} // namespace

// nzm: Momentum points in the Vertical        [-]
// nzt: Thermodynamic points in the Vertical        [-]
// dt: Model timestep        [s]
// cloud_frac: Cloud fraction                        [-]
// rcm: Cloud water mixing ratio              [kg/kg]
// rvm: Vapor water mixing ratio              [kg/kg]
// thlm: Liquid potential temperature          [K]
// wm: Mean vertical velocity                [m/s]
// exner: Exner function                        [-]
// rrm_evap: Evaporation of rain                   [kg/kg/s] It is expected that this variable is
// negative, as that is the convention in Morrison microphysics
// pdf_params: PDF parameters
// rtp2_mc: Tendency of <rt'^2> due to evaporation   [(kg/kg)^2/s]
// thlp2_mc: Tendency of <thl'^2> due to evaporation  [K^2/s]
// wprtp_mc: Tendency of <w'rt'> due to evaporation   [m*(kg/kg)/s^2]
// wpthlp_mc: Tendency of <w'thl'> due to evaporation  [m*K/s^2]
// rtpthlp_mc: Tendency of <rt'thl'> due to evaporation [K*(kg/kg)/s]
void update_xp2_mc(const grid& gr, int nzm, int nzt, int ngrdcol, real_t dt,
                   Real2D cloud_frac,                                                         // In
                   Real2D rcm, Real2D rvm, Real2D thlm,                                       // In
                   Real2D wm, Real2D exner, Real2D rrm_evap, const pdf_parameter& pdf_params, // In
                   Real2D rtp2_mc, Real2D thlp2_mc, Real2D wprtp_mc, Real2D wpthlp_mc, // InOut
                   Real2D rtpthlp_mc)                                                  // InOut
{
  // Description:
  // This subroutine is for use when l_morr_xp2_mc = .true.
  // The effects of rain evaporation on rtp2 and thlp2 are included by
  // assuming rain falls through the moist (cold) portion of the pdf.
  // This is accomplished by defining a precip_fraction and assuming a double
  // delta shaped pdf, such that the evaporation makes the moist component
  // moister and the colder component colder. Calculations are done using
  // variables on the zt grid, and the outputs are on the zm grid --storer
  //  -------------------- Input Variables --------------------
  //  -------------------- Input/Output Variables --------------------
  //  -------------------- Local Variables --------------------
  Kokkos::View<real_t**, Kokkos::LayoutLeft> temp_rtp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "update_xp2_mc::temp_rtp2"), ngrdcol,
      nzt); // Used only to calculate rtp2_mc  [(kg/kg)^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> temp_thlp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "update_xp2_mc::temp_thlp2"), ngrdcol,
      nzt); // Used to calculate thlp2_mc      [K^2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> temp_wp2(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "update_xp2_mc::temp_wp2"), ngrdcol,
      nzt); // Used to calculate wpxp_mc       [m^2/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rtp2_mc_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "update_xp2_mc::rtp2_mc_zt"), ngrdcol,
      nzt); // Calculated on the zt grid             [(kg/kg)^2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> thlp2_mc_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "update_xp2_mc::thlp2_mc_zt"), ngrdcol,
      nzt); // Calculated on the zt grid             [(kg/kg)^2/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wprtp_mc_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "update_xp2_mc::wprtp_mc_zt"), ngrdcol,
      nzt); // Calculated on the zt grid             [m*(kg/kg)/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> wpthlp_mc_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "update_xp2_mc::wpthlp_mc_zt"), ngrdcol,
      nzt); // Calcualted on the zt grid             [m*K/s^2]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rtpthlp_mc_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "update_xp2_mc::rtpthlp_mc_zt"), ngrdcol,
      nzt); // Calculated on the zt grid             [K*(kg/kg)/s]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> precip_frac_double_delta(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "update_xp2_mc::precip_frac_double_delta"),
      ngrdcol, nzt); // Precipitation fraction for a double delta [-]
  Kokkos::View<real_t**, Kokkos::LayoutLeft> pf_const(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "update_xp2_mc::pf_const"), ngrdcol,
      nzt); // ( 1 - pf )/( pf )                    [-]
  int i;
  int k;
  const auto gr_k_ub_zt              = gr.k_ub_zt;
  const auto gr_grid_dir_indx        = gr.grid_dir_indx;
  const auto gr_k_lb_zt              = gr.k_lb_zt;
  const auto pdf_params_mixt_frac    = pdf_params.mixt_frac;
  const auto pdf_params_rt_1         = pdf_params.rt_1;
  const auto pdf_params_varnce_rt_1  = pdf_params.varnce_rt_1;
  const auto pdf_params_rt_2         = pdf_params.rt_2;
  const auto pdf_params_varnce_rt_2  = pdf_params.varnce_rt_2;
  const auto pdf_params_thl_1        = pdf_params.thl_1;
  const auto pdf_params_varnce_thl_1 = pdf_params.varnce_thl_1;
  const auto pdf_params_thl_2        = pdf_params.thl_2;
  const auto pdf_params_varnce_thl_2 = pdf_params.varnce_thl_2;
  const auto pdf_params_w_1          = pdf_params.w_1;
  const auto pdf_params_varnce_w_1   = pdf_params.varnce_w_1;
  const auto pdf_params_w_2          = pdf_params.w_2;
  const auto pdf_params_varnce_w_2   = pdf_params.varnce_w_2;

  // -------------------- Begin Code --------------------

  // Calculate precip_frac_double_delta
  Kokkos::parallel_for(
      "update_xp2_mc::array1", ((precip_frac_double_delta.extent(0) - (1)) / (1) + 1),
      KOKKOS_LAMBDA(int cell) {
        int remaining     = cell;
        const int slice_0 = remaining % ((precip_frac_double_delta.extent(0) - (1)) / (1) + 1);
        remaining /= ((precip_frac_double_delta.extent(0) - (1)) / (1) + 1);
        precip_frac_double_delta((1 - 1 + slice_0 * (1)), (gr_k_ub_zt - 1)) = real_t(0.0);
      });
  // Keep the Native downward recurrence serial within each independent column.
  Kokkos::parallel_for(
      "update_xp2_mc::precipitation_fraction", ngrdcol, KOKKOS_LAMBDA(const int column) {
        for (int k = gr_k_ub_zt - gr_grid_dir_indx;
             gr_grid_dir_indx > 0 ? k >= gr_k_lb_zt : k <= gr_k_lb_zt; k -= gr_grid_dir_indx) {
          if (cloud_frac(column, k - 1) > cloud_frac_min)
            precip_frac_double_delta(column, k - 1) = cloud_frac(column, k - 1);
          else
            precip_frac_double_delta(column, k - 1) =
                precip_frac_double_delta(column, k + gr_grid_dir_indx - 1);
        }
      });

  // pf_const is calculated so that when precip_frac_double_delta = 0, rtp2_mc and
  // thlp2_mc will both be zero.  This also avoids a divide by zero error
  Kokkos::parallel_for(
      "update_xp2_mc::loop2", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzt);
        remaining /= (nzt);
        if (precip_frac_double_delta(i - 1, k - 1) > cloud_frac_min) {
          pf_const(i - 1, k - 1) = (real_t(1.0) - precip_frac_double_delta(i - 1, k - 1)) /
                                   precip_frac_double_delta(i - 1, k - 1);
        } else {
          pf_const(i - 1, k - 1) = real_t(0.0);
        }
      });

  // Include effects of rain evaporation on rtp2
  Kokkos::parallel_for(
      "update_xp2_mc::loop3", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzt);
        remaining /= (nzt);
        temp_rtp2(i - 1, k - 1) =
            pdf_params_mixt_frac(i - 1, k - 1) *
                (((pdf_params_rt_1(i - 1, k - 1) - (rcm(i - 1, k - 1) + rvm(i - 1, k - 1))) *
                  (pdf_params_rt_1(i - 1, k - 1) - (rcm(i - 1, k - 1) + rvm(i - 1, k - 1)))) +
                 pdf_params_varnce_rt_1(i - 1, k - 1)) +
            (real_t(1.0) - pdf_params_mixt_frac(i - 1, k - 1)) *
                (((pdf_params_rt_2(i - 1, k - 1) - (rcm(i - 1, k - 1) + rvm(i - 1, k - 1))) *
                  (pdf_params_rt_2(i - 1, k - 1) - (rcm(i - 1, k - 1) + rvm(i - 1, k - 1)))) +
                 pdf_params_varnce_rt_2(i - 1, k - 1));

        rtp2_mc_zt(i - 1, k - 1) =
            (rrm_evap(i - 1, k - 1) * rrm_evap(i - 1, k - 1)) * pf_const(i - 1, k - 1) * dt +
            real_t(2.0) * Kokkos::abs(rrm_evap(i - 1, k - 1)) *
                Kokkos::sqrt(temp_rtp2(i - 1, k - 1) * pf_const(i - 1, k - 1));
        // use absolute value of evaporation, as evaporation will add
        // to rt_1
      });

  Kokkos::deep_copy(rtp2_mc, zt2zm_api(nzm, nzt, ngrdcol, gr, rtp2_mc_zt));

  // Include the effects of rain evaporation on thlp2
  Kokkos::parallel_for(
      "update_xp2_mc::loop4", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzt);
        remaining /= (nzt);
        temp_thlp2(i - 1, k - 1) = pdf_params_mixt_frac(i - 1, k - 1) *
                                       (((pdf_params_thl_1(i - 1, k - 1) - thlm(i - 1, k - 1)) *
                                         (pdf_params_thl_1(i - 1, k - 1) - thlm(i - 1, k - 1))) +
                                        pdf_params_varnce_thl_1(i - 1, k - 1)) +
                                   (real_t(1.0) - pdf_params_mixt_frac(i - 1, k - 1)) *
                                       (((pdf_params_thl_2(i - 1, k - 1) - thlm(i - 1, k - 1)) *
                                         (pdf_params_thl_2(i - 1, k - 1) - thlm(i - 1, k - 1))) +
                                        pdf_params_varnce_thl_2(i - 1, k - 1));

        thlp2_mc_zt(i - 1, k - 1) =
            ((rrm_evap(i - 1, k - 1) * Lv / (Cp * exner(i - 1, k - 1))) *
             (rrm_evap(i - 1, k - 1) * Lv / (Cp * exner(i - 1, k - 1)))) *
                pf_const(i - 1, k - 1) * dt +
            real_t(2.0) * Kokkos::abs(rrm_evap(i - 1, k - 1)) * Lv / (Cp * exner(i - 1, k - 1)) *
                Kokkos::sqrt(temp_thlp2(i - 1, k - 1) * pf_const(i - 1, k - 1));
      });

  Kokkos::deep_copy(thlp2_mc, zt2zm_api(nzm, nzt, ngrdcol, gr, thlp2_mc_zt));

  // Include effects of rain evaporation on other moments (wprtp, wpthlp, and
  // rtpthlp - added 07/13 rstorer
  Kokkos::parallel_for(
      "update_xp2_mc::loop5", std::max(0, nzt) * std::max(0, ngrdcol),
      KOKKOS_LAMBDA(const int cell) {
        int remaining = cell;
        const int i   = 1 + remaining % (ngrdcol);
        remaining /= (ngrdcol);
        const int k = 1 + remaining % (nzt);
        remaining /= (nzt);
        temp_wp2(i - 1, k - 1) = pdf_params_mixt_frac(i - 1, k - 1) *
                                     (((pdf_params_w_1(i - 1, k - 1) - wm(i - 1, k - 1)) *
                                       (pdf_params_w_1(i - 1, k - 1) - wm(i - 1, k - 1))) +
                                      pdf_params_varnce_w_1(i - 1, k - 1)) +
                                 (real_t(1.0) - pdf_params_mixt_frac(i - 1, k - 1)) *
                                     (((pdf_params_w_2(i - 1, k - 1) - wm(i - 1, k - 1)) *
                                       (pdf_params_w_2(i - 1, k - 1) - wm(i - 1, k - 1))) +
                                      pdf_params_varnce_w_2(i - 1, k - 1));

        wprtp_mc_zt(i - 1, k - 1) = Kokkos::abs(rrm_evap(i - 1, k - 1)) *
                                    Kokkos::sqrt(pf_const(i - 1, k - 1)) *
                                    Kokkos::sqrt(temp_wp2(i - 1, k - 1));

        wpthlp_mc_zt(i - 1, k - 1) =
            -real_t(1.0) * Lv / (Cp * exner(i - 1, k - 1)) * Kokkos::abs(rrm_evap(i - 1, k - 1)) *
            Kokkos::sqrt(pf_const(i - 1, k - 1)) * Kokkos::sqrt(temp_wp2(i - 1, k - 1));

        rtpthlp_mc_zt(i - 1, k - 1) =
            -real_t(1.0) * Kokkos::abs(rrm_evap(i - 1, k - 1)) *
                Kokkos::sqrt(pf_const(i - 1, k - 1)) *
                ((Lv / (Cp * exner(i - 1, k - 1))) * Kokkos::sqrt(temp_rtp2(i - 1, k - 1)) +
                 Kokkos::sqrt(temp_thlp2(i - 1, k - 1))) -
            (Lv / (Cp * exner(i - 1, k - 1))) * pf_const(i - 1, k - 1) *
                ((rrm_evap(i - 1, k - 1)) * (rrm_evap(i - 1, k - 1))) * dt;
      });

  Kokkos::deep_copy(wprtp_mc, zt2zm_api(nzm, nzt, ngrdcol, gr, wprtp_mc_zt));
  Kokkos::deep_copy(wpthlp_mc, zt2zm_api(nzm, nzt, ngrdcol, gr, wpthlp_mc_zt));
  Kokkos::deep_copy(rtpthlp_mc, zt2zm_api(nzm, nzt, ngrdcol, gr, rtpthlp_mc_zt));
}

} // namespace clubb_kokkos
