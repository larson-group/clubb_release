// Predict and diagnose third-order scalar and horizontal-velocity moments.
//
// advance_xp3 advances rt/thl/passive-scalar third moments with the Native
// simplified equation, then diagnoses u/v third moments. diagnose_xp3 applies
// the Larson–Golaz ansatz to all those fields with the selected PDF width.
// Inputs/outputs are execution-space LayoutLeft Views indexed (column,level)
// or (column,level,scalar). Native units are listed in the interfaces.
// Private production/accumulation expressions execute within independent cell
// kernels. Stats retain the existing Fortran backend; no state is kept here.
// Native predictive/LMM and steady-state branches are retained. The last
// storage level is the Native upper boundary in either grid direction.
// Reference: Larson and Golaz (2005) for the diagnostic ansatz.
#include "advance_xp3_module.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "model_flags.hpp"
#include "Skx_module.hpp"
#include "f2k/procedure_bridge.hpp"
#include <Kokkos_MathematicalFunctions.hpp>
#include <string>
namespace clubb_kokkos {
namespace {
constexpr int xp3_rtp3   = 1; // Named constant for solving rtp3.
constexpr int xp3_thlp3  = 2; // Named constant for solving thlp3.
constexpr int xp3_sclrp3 = 3; // Named constant for solving sclrp3.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8>
void advance_xp3_simplified(int nzm, int nzt, int ngrdcol, const grid& gr, int solve_type,
                            real_t dt, Array0 xm, Array1 xp2, Array2 wpxp, Array3 wpxp2,
                            Array4 rho_ds_zm, Array5 invrs_rho_ds_zt, Array6 invrs_tau_zt,
                            Array7 tau_max_zt, real_t x_tol, bool l_lmm_stepping,
                            f2k::stats_type& stats, Array8 xp3);
KOKKOS_INLINE_FUNCTION real_t term_tp_rhs(real_t xp2_zt, real_t wpxpp1, real_t wpxp,
                                          real_t rho_ds_zmp1, real_t rho_ds_zm,
                                          real_t invrs_rho_ds_zt, real_t invrs_dzt);
KOKKOS_INLINE_FUNCTION real_t term_ac_rhs(real_t xm_zmp1, real_t xm_zm, real_t wpxp2,
                                          real_t invrs_dzt);
} // namespace

void advance_xp3(
    int nzm,         // In: Number of momentum vertical levels
    int nzt,         // In: Number of thermodynamic vertical levels
    int ngrdcol,     // In: Number of grid columns
    int sclr_dim,    // In: Number of passive scalars
    Real1D sclr_tol, // In: Threshold(s) on the passive scalars  [units vary]
    const grid& gr,  // In
    real_t dt,       // In: Model timestep                            [s]
    Real2D rtm,      // In: Mean (overall) of rt (thermo. levels)                [kg/kg]
    Real2D thlm,     // In: Mean (overall) of thl (thermo. levels)               [K]
    Real2D rtp2,     // In: Variance (overall) of rt (m-levs.)                   [kg^2/kg^2]
    Real2D thlp2,    // In: Variance (overall) of thl (m-levs.)                  [K^2]
    Real2D wprtp,    // In: Turbulent flux of rt (momentum levs.)                [m/s kg/kg]
    Real2D wpthlp,   // In: Turbulent flux of thl (momentum levs.)               [m/s K]
    // In: <w'rt'^2> (thermodynamic levels)                     [m/s(kg/kg)^2]
    Real2D wprtp2,
    Real2D wpthlp2,   // In: <w'thl'^2> (thermodynamic levels)                    [m/s K^2]
    Real2D rho_ds_zm, // In: Dry, static density on momentum levels               [kg/m^3]
    // In: Inv. dry, static density at thermo. levels           [m^3/kg]
    Real2D invrs_rho_ds_zt,
    Real2D invrs_tau_zt,    // In: Inverse time-scale tau on thermodynamic levels       [1/s]
    Real2D tau_max_zt,      // In: Max. allowable eddy dissipation time scale on t-levs [s]
    Real3D sclrm,           // In: Mean (overall) of sclr (thermo. levels) [sclr units]
    Real3D sclrp2,          // In: Variance (overall) of sclr (m-levs.)    [(sclr units)^2]
    Real3D wpsclrp,         // In: Turbulent flux of sclr (momentum levs.) [m/s(sclr units)]
    Real3D wpsclrp2,        // In: <w'sclr'^2> (thermodynamic levels)      [m/s(sclr units)^2]
    Real2D wp2,             // In: w'^2 (momentum levels)                              [m^2/s^2]
    Real2D wp3,             // In: w'^3 on thermo. grid     [m^3/s^3]
    Real2D upwp,            // In: u'w' (momentum levels)                              [m^2/s^2]
    Real2D vpwp,            // In: v'w' (momentum levels)                              [m^2/s^2]
    Real2D up2,             // In: u'^2 (momentum levels)                              [m^2/s^2]
    Real2D vp2,             // In: v'^2 (momentum levels)                              [m^2/s^2]
    Real2D thvm,            // In: Virtual potential temperature on thermodynamic levels [K]
    Real2D clubb_params,    // In: Array of CLUBB's tunable parameters    [units vary]
    bool l_lmm_stepping,    // In: Apply Linear Multistep Method (LMM) Stepping
    f2k::stats_type& stats, // Inout
    Real2D rtp3,            // Inout: <rt'^3> (thermodynamic levels)     [kg^3/kg^3]
    Real2D thlp3,           // Inout: <thl'^3> (thermodynamic levels)    [K^3]
    Real3D sclrp3,          // Inout: <sclr'^3> (thermodynamic levels)    [(sclr units)^3]
    Real2D up3,             // Inout: u'^3 (thermodynamic levels)        [m^3/s^3]
    // Inout: v'^3 (thermodynamic levels)        [m^3/s^3]
    Real2D vp3) {
  // Description:
  // Advance <rt'^3>, <thl'^3>, and <sclr'^3> one model timestep using a
  // simplified form of the <x'^3> predictive equation.  The simplified <x'^3>
  // equation can either be advanced from its previous value or calculated
  // using a steady-state approximation.
  // References:
  //-----------------------------------------------------------------------
  //----------------------- Local Variables -----------------------
  Real2D thvm_zm;      // Virtual potential temperature on momentum levs. [K]
  Real2D ddzm_thvm_zm; // d(thvm_zm)/dz [K/m]
  Real2D brunt_vaisala_freq_sqd_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp3::brunt_vaisala_freq_sqd_zt"),
      ngrdcol, nzt); // Buoyancy frequency squared on t-levs. [s^-2]
  Real2D Skw_zt(Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp3::Skw_zt"), ngrdcol,
                nzt); // w skewness on thermodynamic levels [-]
  Real2D wp2_zt;      // w'^2 on thermo. grid [m^2/s^2]
  Real2D thlp2_zt;    // thl'^2 on thermo. grid [K^2]
  Real2D wpthlp_zt;   // w'thl' on thermo. grid [m K/s]
  Real2D wprtp_zt;    // w'rt' on thermo. grid [m kg/(kg s)]
  Real2D rtp2_zt;     // rt'^2 on thermo. grid [(kg/kg)^2]
  Real2D upwp_zt;     // u'w' on thermo. grid [m^2/s^2]
  Real2D vpwp_zt;     // v'w' on thermo. grid [m^2/s^2]
  Real2D up2_zt;      // u'^2 on thermo. grid [m^2/s^2]
  Real2D vp2_zt;      // v'^2 on thermo. grid [m^2/s^2]
  Real2D xp3_coef_fnc(Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp3::xp3_coef_fnc"),
                      ngrdcol, nzt); // Coefficient in simple xp3 equation [-]
  wp2_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, wp2, w_tol_sqd);

  // Advance <rt'^3> one model timestep or calculate <rt'^3> using a
  // steady-state approximation.
  advance_xp3_simplified(nzm, nzt, ngrdcol, gr, xp3_rtp3, dt, rtm, rtp2, wprtp, wprtp2, rho_ds_zm,
                         invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt, rt_tol, l_lmm_stepping, stats,
                         rtp3);

  // Advance <thl'^3> one model timestep or calculate <thl'^3> using a
  // steady-state approximation.
  advance_xp3_simplified(nzm, nzt, ngrdcol, gr, xp3_thlp3, dt, thlm, thlp2, wpthlp, wpthlp2,
                         rho_ds_zm, invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt, thl_tol,
                         l_lmm_stepping, stats, thlp3);

  // Advance <sclr'^3> one model timestep or calculate <sclr'^3> using a
  // steady-state approximation. Stage scalar thresholds at the typed scalar
  // boundary once; the source's per-scalar lower-call order is retained.
  auto sclr_tol_host = to_host_space(sclr_tol);
  for (int sclr = 0; sclr < sclr_dim; ++sclr) {
    advance_xp3_simplified(nzm, nzt, ngrdcol, gr, xp3_sclrp3, dt,
                           Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr),
                           Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, sclr),
                           Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr),
                           Kokkos::subview(wpsclrp2, Kokkos::ALL, Kokkos::ALL, sclr), rho_ds_zm,
                           invrs_rho_ds_zt, invrs_tau_zt, tau_max_zt, sclr_tol_host(sclr),
                           l_lmm_stepping, stats,
                           Kokkos::subview(sclrp3, Kokkos::ALL, Kokkos::ALL, sclr));
  }

  // Use a modified form of the Larson and Golaz (2005) ansatz for the
  // ADG1 PDF to calculate <u'^3> and <v'^3> for another type of PDF.
  Skx_func(nzt, ngrdcol, wp2_zt, wp3, w_tol, clubb_params, Skw_zt);
  wpthlp_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, wpthlp);
  wprtp_zt  = zm2zt_api(nzm, nzt, ngrdcol, gr, wprtp);
  // Positive def. quantity.
  thlp2_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, thlp2, thl_tol * thl_tol);
  // Positive def. quantity.
  rtp2_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, rtp2, rt_tol * rt_tol);
  upwp_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, upwp);
  vpwp_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, vpwp);
  // Positive def. quantity.
  up2_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, up2, w_tol_sqd);
  // Positive def. quantity.
  vp2_zt       = zm2zt_api(nzm, nzt, ngrdcol, gr, vp2, w_tol_sqd);
  thvm_zm      = zt2zm_api(nzm, nzt, ngrdcol, gr, thvm, zero_threshold);
  ddzm_thvm_zm = ddzm(nzm, nzt, ngrdcol, gr, thvm_zm);
  Kokkos::parallel_for(
      "advance_xp3::buoyancy_frequency", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        const real_t stability          = (grav / thvm(i, k)) * ddzm_thvm_zm(i, k);
        brunt_vaisala_freq_sqd_zt(i, k) = stability > zero ? stability : zero;
      });

  // The xp3_coef_fnc provides some extra tunability to the simple xp3
  // equation.
  // When xp3_coef_fnc goes to 0, the value of Skx goes to the smallest
  // magnitude permitted by the function. When xp3_coef_fnc goes to 1, the
  // magnitude of Skx becomes huge.
  Kokkos::parallel_for(
      "advance_xp3::xp3_coef_fnc", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        xp3_coef_fnc(i, k) = clubb_params(i, ixp3_coef_base) +
                             (one - clubb_params(i, ixp3_coef_slope)) *
                                 (one - Kokkos::exp(brunt_vaisala_freq_sqd_zt(i, k) /
                                                    clubb_params(i, ixp3_coef_slope)));
      });
  xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, upwp_zt, wp2_zt, up2_zt, xp3_coef_fnc, clubb_params,
                     w_tol, up3);
  xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, vpwp_zt, wp2_zt, vp2_zt, xp3_coef_fnc, clubb_params,
                     w_tol, vp3);

  if (stats.flags().l_sample) {
    f2k::stats_update("thlp2_zt", thlp2_zt, stats);
    f2k::stats_update("wpthlp_zt", wpthlp_zt, stats);
    f2k::stats_update("wprtp_zt", wprtp_zt, stats);
    f2k::stats_update("rtp2_zt", rtp2_zt, stats);
    f2k::stats_update("up2_zt", up2_zt, stats);
    f2k::stats_update("vp2_zt", vp2_zt, stats);
    f2k::stats_update("upwp_zt", upwp_zt, stats);
    f2k::stats_update("vpwp_zt", vpwp_zt, stats);
  }
}

void diagnose_xp3(int nzm,                // In: Number of momentum levels
                  int nzt,                // In: Number of thermodynamic levels
                  int ngrdcol,            // In: Number of grid columns
                  int sclr_dim,           // In: Number of passive scalars
                  Real1D sclr_tol,        // In: Thresholds for passive scalars [units vary]
                  const grid& gr,         // In: Grid structure
                  int iiPDF_type,         // In: PDF type selector
                  Real2D clubb_params,    // In: Array of CLUBB tunable parameters [units vary]
                  Real2D wp2,             // In: Variance of vertical velocity [m^2/s^2]
                  Real2D wp3,             // In: Third moment of vertical velocity [m^3/s^3]
                  Real2D thvm,            // In: Virtual potential temperature [K]
                  Real2D wprtp,           // In: Turbulent flux of rt [m kg/(kg s)]
                  Real2D wpthlp,          // In: Turbulent flux of thl [m K/s]
                  Real2D rtp2,            // In: Variance of rt [(kg/kg)^2]
                  Real2D thlp2,           // In: Variance of thl [K^2]
                  Real2D upwp,            // In: Turbulent flux of u [m^2/s^2]
                  Real2D vpwp,            // In: Turbulent flux of v [m^2/s^2]
                  Real2D up2,             // In: Variance of u wind [m^2/s^2]
                  Real2D vp2,             // In: Variance of v wind [m^2/s^2]
                  Real2D sigma_sqd_w,     // In: PDF width parameter [-]
                  Real3D wpsclrp,         // In: Turbulent flux of passive scalars [m sclr/s]
                  Real3D sclrp2,          // In: Variance of passive scalars [units vary]
                  f2k::stats_type& stats, // Inout: Statistics accumulator
                  Real2D rtp3,            // Inout: Third moment of rt [(kg/kg)^3]
                  Real2D thlp3,           // Inout: Third moment of thl [K^3]
                  Real2D up3,             // Inout: Third moment of u wind [m^3/s^3]
                  Real2D vp3,             // Inout: Third moment of v wind [m^3/s^3]
                  // Inout: Third moments of passive scalars [units vary]
                  Real3D sclrp3) {
  // Description:
  //   Diagnose third-order scalar and horizontal-velocity moments from the
  //   current second-order moments and PDF width information.
  //
  // References:
  //   Larson and Golaz (2005) for the diagnostic ansatz.
  //-----------------------------------------------------------------------
  //----------------------- Local Variables -----------------------
  Real2D thvm_zm;      // Virtual potential temperature on momentum levs. [K]
  Real2D ddzm_thvm_zm; // d(thvm_zm)/dz [K/m]
  Real2D brunt_vaisala_freq_sqd_zt(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "diagnose_xp3::brunt_vaisala_freq_sqd_zt"),
      ngrdcol, nzt); // Buoyancy frequency squared on t-levs. [s^-2]
  Real2D wp2_zt;     // Variance of vertical velocity on thermodynamic levels [m^2/s^2]
  Real2D Skw_zt(Kokkos::view_alloc(Kokkos::WithoutInitializing, "diagnose_xp3::Skw_zt"), ngrdcol,
                nzt);    // w skewness on thermodynamic levels [-]
  Real2D wpthlp_zt;      // w'thl' on thermo. grid [m K/s]
  Real2D wprtp_zt;       // w'rt' on thermo. grid [m kg/(kg s)]
  Real2D thlp2_zt;       // thl'^2 on thermo. grid [K^2]
  Real2D rtp2_zt;        // rt'^2 on thermo. grid [(kg/kg)^2]
  Real2D upwp_zt;        // u'w' on thermo. grid [m^2/s^2]
  Real2D vpwp_zt;        // v'w' on thermo. grid [m^2/s^2]
  Real2D up2_zt;         // u'^2 on thermo. grid [m^2/s^2]
  Real2D vp2_zt;         // v'^2 on thermo. grid [m^2/s^2]
  Real2D sigma_sqd_w_zt; // PDF width parameter (thermodynamic levels) [-]
  Real2D xp3_coef_fnc(Kokkos::view_alloc(Kokkos::WithoutInitializing, "diagnose_xp3::xp3_coef_fnc"),
                      ngrdcol, nzt); // Coefficient in simple xp3 equation [-]
  Real2D wpsclrp_zt;                 // Scalar flux on thermo. levels [un. vary]
  Real2D sclrp2_zt;                  // Scalar variance on thermo. levels [un. vary]
  //----------------------------- Begin Code ------------------------------
  wp2_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, wp2, w_tol_sqd);
  // The ADG1 PDF must use this option.
  Skx_func(nzt, ngrdcol, wp2_zt, wp3, w_tol, clubb_params, Skw_zt);
  wpthlp_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, wpthlp);
  wprtp_zt  = zm2zt_api(nzm, nzt, ngrdcol, gr, wprtp);
  // Positive def. quantity.
  thlp2_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, thlp2, thl_tol * thl_tol);
  // Positive def. quantity.
  rtp2_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, rtp2, rt_tol * rt_tol);
  upwp_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, upwp);
  vpwp_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, vpwp);
  // Positive def. quantity.
  up2_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, up2, w_tol_sqd);
  // Positive def. quantity.
  vp2_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, vp2, w_tol_sqd);
  // The typed ansatz API receives one threshold scalar per Native scalar call.
  auto sclr_tol_host = to_host_space(sclr_tol);

  if (iiPDF_type == iiPDF_ADG1) {
    // Use the Larson and Golaz (2005) ansatz for the ADG1 PDF to
    // calculate <rt'^3>, <thl'^3>, <u'^3>, <v'^3>, and <sclr'^3>.
    sigma_sqd_w_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, sigma_sqd_w, zero_threshold);
    xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, wpthlp_zt, wp2_zt, thlp2_zt, sigma_sqd_w_zt,
                       clubb_params, thl_tol, thlp3);
    xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, wprtp_zt, wp2_zt, rtp2_zt, sigma_sqd_w_zt,
                       clubb_params, rt_tol, rtp3);
    xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, upwp_zt, wp2_zt, up2_zt, sigma_sqd_w_zt, clubb_params,
                       w_tol, up3);
    xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, vpwp_zt, wp2_zt, vp2_zt, sigma_sqd_w_zt, clubb_params,
                       w_tol, vp3);
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      // Keep the Native ADG1 scalar-flux minimum as well as the variance minimum.
      wpsclrp_zt =
          zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr),
                    sclr_tol_host(sclr) * sclr_tol_host(sclr));
      sclrp2_zt =
          zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, sclr),
                    sclr_tol_host(sclr) * sclr_tol_host(sclr));
      // Existing array owner handles the typed ansatz output boundary. A compatible
      // contiguous scalar section aliases its allocation; otherwise restore it.
      auto sclrp3_section = Kokkos::subview(sclrp3, Kokkos::ALL, Kokkos::ALL, sclr);
      auto sclrp3_field   = to_execution_space(sclrp3_section, false);
      xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, wpsclrp_zt, wp2_zt, sclrp2_zt, sigma_sqd_w_zt,
                         clubb_params, sclr_tol_host(sclr), sclrp3_field);
      from_execution_space(sclrp3_section, sclrp3_field);
    }
  } else {
    // Use a modified form of the Larson and Golaz (2005) ansatz for the
    // ADG1 PDF to calculate <u'^3> and <v'^3> for another type of PDF.
    thvm_zm      = zt2zm_api(nzm, nzt, ngrdcol, gr, thvm, zero_threshold);
    ddzm_thvm_zm = ddzm(nzm, nzt, ngrdcol, gr, thvm_zm);
    Kokkos::parallel_for(
        "diagnose_xp3::buoyancy_frequency", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          const real_t stability          = (grav / thvm(i, k)) * ddzm_thvm_zm(i, k);
          brunt_vaisala_freq_sqd_zt(i, k) = stability > zero ? stability : zero;
        });
    // Initialize sigma_sqd_w_zt to zero so we don't break output.
    sigma_sqd_w_zt =
        Real2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "diagnose_xp3::sigma_sqd_w_zt"),
               ngrdcol, nzt);
    Kokkos::deep_copy(sigma_sqd_w_zt, zero);

    // The xp3_coef_fnc is used in place of sigma_sqd_w_zt when the
    // ADG1 PDF is not being used. The xp3_coef_fnc provides some extra
    // tunability to the simple xp3 equation.
    // When xp3_coef_fnc goes to 0, the value of Skx goes to the smallest
    // magnitude permitted by the function. When xp3_coef_fnc goes to 1,
    // the magnitude of Skx becomes huge.
    // The value of Skx becomes large near cloud top, where there is a
    // higher degree of static stability. The exp{ } portion of the
    // xp3_coef_fnc allows the xp3_coef_fnc to become larger in regions
    // of high static stability, producing larger magnitude values of Skx.
    Kokkos::parallel_for(
        "diagnose_xp3::xp3_coef_fnc", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          xp3_coef_fnc(i, k) = clubb_params(i, ixp3_coef_base) +
                               (one - clubb_params(i, ixp3_coef_slope)) *
                                   (one - Kokkos::exp(brunt_vaisala_freq_sqd_zt(i, k) /
                                                      clubb_params(i, ixp3_coef_slope)));
        });
    xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, wpthlp_zt, wp2_zt, thlp2_zt, xp3_coef_fnc,
                       clubb_params, thl_tol, thlp3);
    xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, wprtp_zt, wp2_zt, rtp2_zt, xp3_coef_fnc, clubb_params,
                       rt_tol, rtp3);
    xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, upwp_zt, wp2_zt, up2_zt, xp3_coef_fnc, clubb_params,
                       w_tol, up3);
    xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, vpwp_zt, wp2_zt, vp2_zt, xp3_coef_fnc, clubb_params,
                       w_tol, vp3);
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      // Native non-ADG1 scalar-flux interpolation has no minimum.
      wpsclrp_zt = zm2zt_api(nzm, nzt, ngrdcol, gr,
                             Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr));
      sclrp2_zt =
          zm2zt_api(nzm, nzt, ngrdcol, gr, Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, sclr),
                    sclr_tol_host(sclr) * sclr_tol_host(sclr));
      auto sclrp3_section = Kokkos::subview(sclrp3, Kokkos::ALL, Kokkos::ALL, sclr);
      auto sclrp3_field   = to_execution_space(sclrp3_section, false);
      xp3_LG_2005_ansatz(nzt, ngrdcol, Skw_zt, wpsclrp_zt, wp2_zt, sclrp2_zt, xp3_coef_fnc,
                         clubb_params, sclr_tol_host(sclr), sclrp3_field);
      from_execution_space(sclrp3_section, sclrp3_field);
    }
  }

  if (stats.flags().l_sample) {
    f2k::stats_update("thlp2_zt", thlp2_zt, stats);
    f2k::stats_update("wpthlp_zt", wpthlp_zt, stats);
    f2k::stats_update("wprtp_zt", wprtp_zt, stats);
    f2k::stats_update("rtp2_zt", rtp2_zt, stats);
    f2k::stats_update("up2_zt", up2_zt, stats);
    f2k::stats_update("vp2_zt", vp2_zt, stats);
    f2k::stats_update("upwp_zt", upwp_zt, stats);
    f2k::stats_update("vpwp_zt", vpwp_zt, stats);
  }
}

namespace {
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8>
void advance_xp3_simplified(
    int nzm,        // In
    int nzt,        // In
    int ngrdcol,    // In
    const grid& gr, // In
    int solve_type, // In: Flag for solving for rtp3, thlp3, or sclrp3
    real_t dt,      // In: Model timestep                            [s]
    Array0 xm,      // In: Mean (overall) of x (thermo. levels)                 [(x units)]
    Array1 xp2,     // In: Variance (overall) of x (m-levs.)       [(x units)^2]
    Array2 wpxp,    // In: Turbulent flux of x (momentum levs.)    [m/s(x units)]
    // In: <w'x'^2> (thermodynamic levels)                      [m/s(x units)^2]
    Array3 wpxp2,
    Array4 rho_ds_zm, // In: Dry, static density on momentum levels  [kg/m^3]
    // In: Inv. dry, static density at thermo. levels           [m^3/kg]
    Array5 invrs_rho_ds_zt,
    Array6 invrs_tau_zt,    // In: Inverse time-scale tau on thermodynamic levels       [1/s]
    Array7 tau_max_zt,      // In: Max. allowable eddy dissipation time scale on t-levs [s]
    real_t x_tol,           // In: Tolerance value of x                           [(x units)]
    bool l_lmm_stepping,    // In: Apply Linear Multistep Method (LMM) Stepping
    f2k::stats_type& stats, // Inout
    // Inout: <x'^3> (thermodynamic levels)    [(x units)^3]
    Array8 xp3) {
  // Description:
  // Predicts the value of <x'^3> using a simplified form of the <x'^3>
  // predictive equation.
  //
  // The full predictive equation for <x'^3>, where <x'^3> can be <rt'^3>,
  // <thl'^3>, or <sclr'^3>, is:
  //
  // d<x'^3>/dt = - <w> * d<x'^3>/dz
  //              - (1/rho_ds) * d( rho_ds * <w'x'^3> )/dz
  //              - 3 * <w'x'^2> * d<x>/dz
  //              + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz
  //              - ( C_xp3_dissipation / tau ) * <x'^3>
  //              + d ( ( K_xp3 + nu_xp3 ) * d<x'^3>/dz )/dz
  //              + 3 * < x'^2 (dx/dt)|_f' >;
  //
  // where (dx/dt)|_f is the "forcing" term, which may include effects such as
  // microphysical effects or radiative effects.  The tunable coefficients are
  // C_xp3_dissipation, K_xp3, and nu_xp3.  The terms are listed as follows:
  //
  // time tendency: d<x'^3>/dt;
  // mean advection: - <w> * d<x'^3>/dz;
  // turbulent advection: - (1/rho_ds) * d( rho_ds * <w'x'^3> )/dz;
  // accumulation: - 3 * <w'x'^2> * d<x>/dz;
  // turbulent production: + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz;
  // turbulent dissipation: - ( C_xp3_dissipation / tau ) * <x'^3>;
  // diffusion: + d ( ( K_xp3 + nu_xp3 ) * d<x'^3>/dz )/dz; and
  // microphysics/other forcing: + 3 * < x'^2 (dx/dt)|_f' >.
  //
  // The microphysics and turbulent advection terms are both found by
  // integration over the subgrid PDF.  This requires new integrated terms.
  // The turbulent advection term may need to be made semi-implicit in order
  // to aid model stability.  This may be difficult to do for <x'^3>.
  // Additionally, if it could be made semi-implicit, it involves a derivative
  // and would require a tridiagonal solver to include contributions from
  // <x'^3> on three grid levels.  While the microphysics term and turbulent
  // advection term are important contributors to <x'^3>, they are being
  // omitted because of the additional complications they bring.
  //
  // The mean advection and diffusion terms also would require a tridiagonal
  // solver in order to make the terms implicit because they involve
  // derivatives and values of <x'^3> on three grid levels.  While tridiagonal
  // solvers are not very computationally expensive, they are still more
  // expensive than a simplified one-line equation.  The mean advection and
  // diffusion terms are also rather small in magnitude, so they are also
  // being neglected.
  //
  // This leaves the following equation:
  //
  // d<x'^3>/dt = - 3 * <w'x'^2> * d<x>/dz
  //              + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz
  //              - ( C_xp3_dissipation / tau ) * <x'^3>;
  //
  // which is a balance of time-tendency, accumulation, turbulent production,
  // and turbulent dissipation.  This equation can be handled semi-implicitly
  // as:
  //
  // ( <x'^3>(t+1) - <x'^3>(t) ) / delta_t
  // = - 3 * <w'x'^2> * d<x>/dz
  //   + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz
  //   - ( C_xp3_dissipation / tau ) * <x'^3>(t+1);
  //
  // which can be rewritten as:
  //
  // ( 1 / delta_t + ( C_xp3_dissipation / tau ) ) * <x'^3>(t+1)
  // = ( <x'^3>(t) / delta_t )
  //   - 3 * <w'x'^2> * d<x>/dz
  //   + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz.
  //
  // The predictive equation can be solved for <x'^3> as:
  //
  // <x'^3>(t+1)
  // = ( ( <x'^3>(t) / delta_t )
  //     - 3 * <w'x'^2> * d<x>/dz
  //     + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz )
  //   / ( 1 / delta_t + ( C_xp3_dissipation / tau ) ).
  //
  // Alternatively, a steady-state approximation can be used, which
  // approximates d<x'^3>/dt = 0.  The equation becomes a balance of
  // accumulation, turbulent production, and turbulent dissipation, and is
  // written as:
  //
  // 0 = - 3 * <w'x'^2> * d<x>/dz
  //     + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz
  //     - ( C_xp3_dissipation / tau ) * <x'^3>.
  //
  // The equation can be solved for <x'^3> as:
  //
  // <x'^3>
  // = ( tau / C_xp3_dissipation )
  //   * ( - 3 * <w'x'^2> * d<x>/dz
  //       + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz ).
  //
  // When the flag l_predict_xp3 is enabled, the predictive version of <x'^3>
  // is used.  When the flag is turned off, the steady-state approximation is
  // used.
  // References:
  //-----------------------------------------------------------------------
  //----------------------- Local Variables -----------------------
  Real2D xp3_old; // Saved <x'^3> (thermodynamic levels)    [(x units)^3]
  Real2D xm_zm;   // Mean of x interpolated to momentum levels     [(x units)]
  Real2D xp2_zt;  // Variance of x interpolated to thermo. levels  [(x units)^2]
  Real2D term_tp(Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp3_simplified::term_tp"),
                 ngrdcol, nzt); // <x'^3> turbulent production term              [(x units)^3/s]
  Real2D term_ac(Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp3_simplified::term_ac"),
                 ngrdcol, nzt); // <x'^3> accumulation term                      [(x units)^3/s]
  // Coefficient in the <x'^3> turbulent dissipation term [-].
  constexpr real_t C_xp3_dissipation = real_t(1.0);
  // Flag to either predict <x'^3> or use steady-state approximation.
  constexpr bool l_predict_xp3 = false;
  std::string name_bt, name_tp, name_ac, name_dp;
  // The saved field is read only by the Native predictive/LMM branch.
  if (l_predict_xp3 && l_lmm_stepping)
    xp3_old =
        Real2D(Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp3_simplified::xp3_old"),
               ngrdcol, nzt);
  // Capture the grid component rather than its owning UDT in device kernels.
  const auto invrs_dzt = gr.invrs_dzt;
  //----------------------- Begin Code -----------------------
  switch (solve_type) {
  case xp3_rtp3:
    // Budget stats for rtp3.
    name_bt = "rtp3_bt";
    name_tp = "rtp3_tp";
    name_ac = "rtp3_ac";
    name_dp = "rtp3_dp";
    break;
  case xp3_thlp3:
    // Budget stats for thlp3.
    name_bt = "thlp3_bt";
    name_tp = "thlp3_tp";
    name_ac = "thlp3_ac";
    name_dp = "thlp3_dp";
    break;
  default:
    // Budgets aren't setup for the passive scalars.
    name_bt = "";
    name_tp = "";
    name_ac = "";
    name_dp = "";
    break;
  }
  if (l_predict_xp3) {
    if (stats.flags().l_sample) {
      Real2D stats_bt_before(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp3_simplified::bt_before"),
          ngrdcol, nzt);
      Kokkos::parallel_for(
          "advance_xp3_simplified::bt_before", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
            const int i = cell % ngrdcol, k = cell / ngrdcol;
            stats_bt_before(i, k) = xp3(i, k) / dt;
          });
      f2k::stats_begin_budget(name_bt, stats_bt_before, stats);
    }
  }
  // Initialize variables.
  Kokkos::deep_copy(term_tp, zero);
  Kokkos::deep_copy(term_ac, zero);
  // Interpolate <x> to momentum levels.
  xm_zm = zt2zm_api(nzm, nzt, ngrdcol, gr, xm, zero_threshold);
  // Interpolate <x'^2> to thermodynamic levels. Positive definite quantity.
  xp2_zt = zm2zt_api(nzm, nzt, ngrdcol, gr, xp2, x_tol * x_tol);

  // Independent cells retain Native k=1:nzt-1, with no vertical recurrence.
  Kokkos::parallel_for(
      "advance_xp3_simplified::terms", ngrdcol * (nzt > 0 ? nzt - 1 : 0), KOKKOS_LAMBDA(int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        // Define the km1 index.
        const int kp1 = k + 1 < nzt - 1 ? k + 1 : nzt - 1;
        // Calculate the <x'^3> turbulent production (tp) term.
        term_tp(i, k) = term_tp_rhs(xp2_zt(i, k), wpxp(i, kp1), wpxp(i, k), rho_ds_zm(i, kp1),
                                    rho_ds_zm(i, k), invrs_rho_ds_zt(i, k), invrs_dzt(i, k));
        // Calculate the <x'^3> accumulation (ac) term.
        term_ac(i, k) = term_ac_rhs(xm_zm(i, kp1), xm_zm(i, k), wpxp2(i, k), invrs_dzt(i, k));
        if (l_predict_xp3) {
          if (l_lmm_stepping)
            xp3_old(i, k) = xp3(i, k);
          // Advance <x'^3> one time step.
          xp3(i, k) = ((xp3(i, k) / dt) + term_tp(i, k) + term_ac(i, k)) /
                      ((one / dt) + (C_xp3_dissipation * invrs_tau_zt(i, k)));
          if (l_lmm_stepping)
            xp3(i, k) = one_half * (xp3_old(i, k) + xp3(i, k));
        } else {
          // Calculate <x'^3> using the steady-state approximation.
          const real_t reciprocal_tau = one / invrs_tau_zt(i, k);
          xp3(i, k) = (reciprocal_tau < tau_max_zt(i, k) ? reciprocal_tau : tau_max_zt(i, k)) *
                      one / C_xp3_dissipation * (term_tp(i, k) + term_ac(i, k));
        }
      });
  // Set Upper Boundary Condition: retain Native storage level nzt in both
  // grid directions, rather than substituting a grid boundary index.
  Kokkos::parallel_for(
      "advance_xp3_simplified::upper_boundary", ngrdcol,
      KOKKOS_LAMBDA(int i) { xp3(i, nzt - 1) = zero; });
  if (stats.flags().l_sample) {
    f2k::stats_update(name_tp, term_tp, stats);
    f2k::stats_update(name_ac, term_ac, stats);
    Real2D stats_dp(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp3_simplified::dissipation"),
        ngrdcol, nzt);
    Kokkos::parallel_for(
        "advance_xp3_simplified::dissipation", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          stats_dp(i, k) = -(C_xp3_dissipation * invrs_tau_zt(i, k)) * xp3(i, k);
        });
    f2k::stats_update(name_dp, stats_dp, stats);
    if (l_predict_xp3) {
      Real2D stats_bt_after(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "advance_xp3_simplified::bt_after"),
          ngrdcol, nzt);
      Kokkos::parallel_for(
          "advance_xp3_simplified::bt_after", ngrdcol * nzt, KOKKOS_LAMBDA(int cell) {
            const int i = cell % ngrdcol, k = cell / ngrdcol;
            stats_bt_after(i, k) = xp3(i, k) / dt;
          });
      f2k::stats_finalize_budget(name_bt, stats_bt_after, stats);
    }
  }
}

KOKKOS_INLINE_FUNCTION real_t
term_tp_rhs(real_t xp2_zt,          // In: <x'^2> interp. to thermo. level (k)     [(x units)^2]
            real_t wpxpp1,          // In: <w'x'> at momentum level (k+1)         [m/s(x units)]
            real_t wpxp,            // In: <w'x'> at momentum level (k)           [m/s(x units)]
            real_t rho_ds_zmp1,     // In: Dry, static density on momentum level (k+1)  [kg/m^3]
            real_t rho_ds_zm,       // In: Dry, static density on momentum level (k)    [kg/m^3]
            real_t invrs_rho_ds_zt, // In: Inv. dry, static density at thermo. lev. (k) [m^3/kg]
            // In: Inverse of grid spacing (k)                     [1/m]
            real_t invrs_dzt) {
  // Description:
  // Turbulent production of <x'^3>:  explicit portion of the code.
  //
  // The d<x'^3>/dt equation contains a turbulent production term:
  //
  // + 3 * ( <x'^2> / rho_ds ) * d( rho_ds * <w'x'> )/dz.
  //
  // The <x'^3> turbulent production term is completely explicit and is
  // discretized as follows:
  //
  // The values of <x'^3> are found on the thermodynamic levels, while the
  // values of <w'x'> and <x'^2> are found on the momentum levels.
  // Additionally, the values of rho_ds_zm are found on the momentum levels,
  // and the values of invrs_rho_ds_zt are found on the thermodynamic levels.
  // The values of <x'^2> are interpolated to the central thermodynamic level
  // as <x'^2>|_zt.  On the momentum levels, the values of <w'x'> are
  // multiplied by rho_ds_zm.  Then, the derivative (d/dz) of
  // rho_ds_zm * <w'x'> is taken over the central thermodynamic level.  At the
  // central thermodynamic level, the derivative is multiplied by
  // invrs_rho_ds_zt, and their product is also multiplied by 3 * <x'^2>|_zt,
  // yielding the desired results.
  //
  // =========wpxpp1=========rho_ds_zmp1===========xp2p1================ m(k+1)
  //
  // --xp3--d( rho_ds_zm * wpxp )/dz--invrs_rho_ds_zt--xp2_zt(interp.)-- t(k)
  //
  // =========wpxp===========rho_ds_zm=============xp2================== m(k)
  //
  // The vertical indices m(k+1), t(k), and m(k) correspond with altitudes
  // zm(k+1), zt(k), and zm(k), respectively.  The letter "t" is used for
  // thermodynamic levels and the letter "m" is used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
  // References:
  //-----------------------------------------------------------------------
  // The <x'^3> turbulent production term.
  return +three * xp2_zt * invrs_rho_ds_zt * invrs_dzt * (rho_ds_zmp1 * wpxpp1 - rho_ds_zm * wpxp);
}

KOKKOS_INLINE_FUNCTION real_t
term_ac_rhs(real_t xm_zmp1, // In: <x> interpolated to momentum level (k+1)  [(x units)]
            real_t xm_zm,   // In: <x> interpolated to momentum level (k)    [(x units)]
            real_t wpxp2,   // In: <w'x'^2> at thermodynamic level (k)       [m/s(x units)^2]
            // In: Inverse of grid spacing (k)               [1/m]
            real_t invrs_dzt) {
  // Description:
  // Accumulation of <x'^3>:  explicit portion of the code.
  //
  // The d<x'^3>/dt equation contains an accumulation term:
  //
  // - 3 * <w'x'^2> * d<x>/dz.
  //
  // The <x'^3> accumulation term is completely explicit and is discretized as
  // follows:
  //
  // The values of <x'^3>, <x>, and <w'x'^2> are found on thermodynamic levels.
  // The values of <x> are interpolated to the intermediate momentum levels as
  // <x>|_zm.  Then, the derivative (d/dz) of <x>|_zm is taken over the
  // central thermodynamic level, where it is multiplied by -3 * <w'x'^2>.
  //
  // ----------------------xmp1----------------------------------------- t(k+1)
  //
  // =========================xm_zmp1(interp.)========================== m(k+1)
  //
  // ----------xp3---------xm---------dxm_zm/dz---------wpxp2----------- t(k)
  //
  // =========================xm_zm(interp.)============================ m(k)
  //
  // ----------------------xmm1----------------------------------------- t(k-1)
  //
  // The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
  // with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
  // The letter "t" is used for thermodynamic levels and the letter "m" is
  // used for momentum levels.
  //
  // invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) )
  // References:
  //-----------------------------------------------------------------------
  // The <x'^3> accumulation term.
  return -three * wpxp2 * invrs_dzt * (xm_zmp1 - xm_zm);
}

} // namespace
} // namespace clubb_kokkos
