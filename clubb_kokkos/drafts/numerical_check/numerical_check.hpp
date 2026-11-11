// Validate model fields/settings while retaining Native diagnostics and status.
//
// Usage:
//   length_check: validate the three mixing-length profiles.
//   pdf_closure_check: check PDF moments and first-column PDF parameters.
//   parameterization_check: validate core inputs and reset beginning-step errors.
//   sfc_varnce_check: validate surface variances and scalar covariances.
//   rad_check: report negative radiation inputs.
//   invalid_model_arrays: report every non-finite selected model profile.
//   calculate_spurious_source: return the signed column conservation residual.
//   check_clubb_settings_api: validate parameters and mixed model settings.
//
// Profile APIs accept strided execution-space input Views without packing.
// NumericalNames stores Native character strings as (character,species).
// All scalar/profile arguments are required. Settings retains Native GPU compile
// gates and public lower boundaries; see INTEGRATION.md. This draft has no build
// or numerical evidence. The existing parameterization interface is unchanged.
#pragma once
#include "array_types.hpp"
#include "array_index.hpp"
#include "err_info_type_module.hpp"
#include "pdf_parameter_module.hpp"
#include "model_flags.hpp"
#include "f2k/stats_bridge.hpp"
#include <string>
namespace clubb_kokkos {
using NumericalProfile1D = Kokkos::View<const real_t*, Kokkos::LayoutStride>;
using NumericalProfile2D = Kokkos::View<const real_t**, Kokkos::LayoutStride>;
using NumericalNames     = Kokkos::View<const char**, Kokkos::LayoutStride>;
// Lscale: Mixing length                 [m]
// Lscale_up: Upward mixing length          [m]
// Lscale_down: Downward mixing length        [m]
// err_info: err_info struct containing err_code and err_header
void length_check(int nzt, NumericalProfile1D Lscale, NumericalProfile1D Lscale_up, // In
                  NumericalProfile1D Lscale_down,                                   // In
                  err_info_type& err_info);                                         // InOut

// wp4: w'^4                  [m^4/s^4]
// wprtp2: w' r_t'^2 [(m/s)(kg/kg)^2]
// wp2rtp: w'^2 r_t'             [(m^2 kg)/(s^2 kg)]
// wpthlp2: w' th_l'^2            [(m K^2)/s]
// wp2thlp: w'^2 th_l'            [(m^2 K)/s^2]
// cloud_frac: Cloud fraction        [-]
// rcm: Mean liquid water     [kg/kg]
// wpthvp: Buoyancy flux         [(K m)/s]
// wp2thvp: w'^2 th_v'            [(m^2 K)/s^2]
// wp2up: w'^2 u'               [m^3/s^3]
// rtpthvp: r_t' th_v'            [(kg K)/kg]
// thlpthvp: th_l' th_v'           [K^2]
// wprcp: w' r_c'               [(m kg)/(s kg)]
// wp2rcp: w'^2 r_c'             [(m^2 kg)/(s^2 kg)]
// rtprcp: r_t' r_c'             [(kg^2)/(kg^2)]
// thlprcp: th_l' r_c'            [(K kg)/kg]
// rcp2: r_c'^2                [(kg^2)/(kg^2)]
// wprtpthlp: w' r_t' th_l'         [(m kg K)/(s kg)]
// pdf_params: PDF parameters          [units vary]
// err_info: err_info struct containing err_code and err_header
void pdf_closure_check(
    int nz, int sclr_dim,                                                           // In
    NumericalProfile1D wp4, NumericalProfile1D wprtp2, NumericalProfile1D wp2rtp,   // In
    NumericalProfile1D wpthlp2,                                                     // In
    NumericalProfile1D wp2thlp, NumericalProfile1D cloud_frac,                      // In
    NumericalProfile1D rcm, NumericalProfile1D wpthvp, NumericalProfile1D wp2thvp,  // In
    NumericalProfile1D wp2up,                                                       // In
    NumericalProfile1D rtpthvp, NumericalProfile1D thlpthvp,                        // In
    NumericalProfile1D wprcp, NumericalProfile1D wp2rcp,                            // In
    NumericalProfile1D rtprcp, NumericalProfile1D thlprcp, NumericalProfile1D rcp2, // In
    NumericalProfile1D wprtpthlp,                                                   // In
    NumericalProfile1D crt_1, NumericalProfile1D crt_2, NumericalProfile1D cthl_1,  // In
    NumericalProfile1D cthl_2, const pdf_parameter& pdf_params,                     // In
    NumericalProfile2D sclrpthvp, NumericalProfile2D sclrprcp,                      // In
    NumericalProfile2D wpsclrp2,                                                    // In
    NumericalProfile2D wpsclrprtp, NumericalProfile2D wpsclrpthlp,                  // In
    NumericalProfile2D wp2sclrp,                                                    // In
    const f2k::stats_type& stats,                                                   // In
    err_info_type& err_info);                                                       // InOut

// clang-format off
void parameterization_check(
    int nzm, int nzt, int ngrdcol, int sclr_dim, int edsclr_dim,                              // In
    Real2D thlm_forcing, Real2D rtm_forcing, Real2D um_forcing,                               // In
    Real2D vm_forcing, Real2D wm_zm, Real2D wm_zt, Real2D p_in_Pa,                            // In
    Real2D rho_zm, Real2D rho, Real2D exner, Real2D rho_ds_zm,                                // In
    Real2D rho_ds_zt, Real2D invrs_rho_ds_zm, Real2D invrs_rho_ds_zt,                         // In
    Real2D thv_ds_zm, Real2D thv_ds_zt, Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D upwp_sfc, // In
    Real1D vpwp_sfc, Real1D p_sfc, Real2D um, Real2D upwp, Real2D vm, Real2D vpwp, Real2D up2,
    Real2D vp2,                                                                   // In
    Real2D rtm, Real2D wprtp, Real2D thlm, Real2D wpthlp, Real2D wp2, Real2D wp3, // In
    Real2D rtp2, Real2D thlp2, Real2D rtpthlp,                                    // In
    const std::string& prefix,                                                    // In
    Real2D wpsclrp_sfc, Real2D wpedsclrp_sfc, Real3D sclrm, Real3D wpsclrp,       // In
    Real3D sclrp2,                                                                // In
    Real3D sclrprtp, Real3D sclrpthlp, Real3D sclrm_forcing, Real3D edsclrm,      // In
    Real3D edsclrm_forcing,                                                       // In
    err_info_type& err_info);                                                     // InOut
// clang-format on

// wp2_sfc: Vertical velocity variance        [m^2/s^2]
// up2_sfc: u'^2                              [m^2/s^2]
// vp2_sfc: v'^2 [m^2/s^2]
// thlp2_sfc: thetal variance                   [K^2]
// rtp2_sfc: rt variance                       [(kg/kg)^2]
// rtpthlp_sfc: thetal rt covariance              [kg K/kg]
// sclrp2_sfc: Passive scalar variance                 [units^2]
// sclrprtp_sfc: Passive scalar r_t covariance           [units kg/kg]
// sclrpthlp_sfc: Passive scalar theta_l covariance       [units K]
// err_info: err_info struct containing err_code and err_header
void sfc_varnce_check(int sclr_dim, real_t wp2_sfc,      // In
                      real_t up2_sfc, real_t vp2_sfc,    // In
                      real_t thlp2_sfc, real_t rtp2_sfc, // In
                      real_t rtpthlp_sfc,                // In
                      Real1D sclrp2_sfc,                 // In
                      Real1D sclrprtp_sfc,               // In
                      Real1D sclrpthlp_sfc,              // In
                      err_info_type& err_info);          // InOut

// thlm: Liquid water potential temperature [K]
// rcm: Liquid Water Mixing Ratio            [kg/kg]
// rtm: Total Water Mixing Ratio             [kg/kg]
// rim: Ice Water Mixing Ratio               [kg/kg]
// cloud_frac: Cloud Fraction                       [-]
// p_in_Pa: Pressure                             [Pa]
// exner: Exner Function                       [-]
// rho_zm: Air Density                          [kg/m^3]
// err_info: err_info struct containing err_code and err_header
void rad_check(int ngrdcol, int nzm, int nzt, Real2D thlm, Real2D rcm, Real2D rtm, Real2D rim, // In
               Real2D cloud_frac, Real2D p_in_Pa, Real2D exner, Real2D rho_zm,                 // In
               err_info_type& err_info); // InOut

// hydromet_list: Hydrometeor names, ten characters per species
// um: eastward grid-mean wind comp. (thermo. levs.)  [m/s]
// vm: northward grid-mean wind comp. (thermo. levs.) [m/s]
// rtm: total water mixing ratio, r_t (thermo. levels) [kg/kg]
// wprtp: w' r_t' (momentum levels)                      [(kg/kg) m/s]
// thlm: liq. water pot. temp., th_l (thermo. levels)   [K]
// wpthlp: w'th_l' (momentum levels)                      [(m/s) K]
// rtp2: r_t'^2 (momentum levels)                       [(kg/kg)^2]
// thlp2: th_l'^2 (momentum levels)                      [K^2]
// rtpthlp: r_t'th_l' (momentum levels)                    [(kg/kg) K]
// wp2: w'^2 (momentum levels)                         [m^2/s^2]
// wp3: w'^3 (thermodynamic levels)                    [m^3/s^3]
// wp2thvp: < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
// wp2up: < w'^2 u' > (thermodynamic levels)             [m^3/s^3]
// rtpthvp: < r_t' th_v' > (momentum levels)               [kg/kg K]
// thlpthvp: < th_l' th_v' > (momentum levels)              [K^2]
// hydromet: Array of hydrometeors                          [units vary]
// sclrm: Passive scalar mean (thermo. levels)              [units vary]
// edsclrm: Eddy passive scalar grid-mean (thermo. levels)   [units vary]
bool invalid_model_arrays(int nzm, int nzt, int hydromet_dim, NumericalNames hydromet_list, // In
                          int sclr_dim, int edsclr_dim,                                     // In
                          NumericalProfile1D um, NumericalProfile1D vm,
                          NumericalProfile1D rtm, // In
                          NumericalProfile1D wprtp, NumericalProfile1D thlm,
                          NumericalProfile1D wpthlp, // In
                          NumericalProfile1D rtp2, NumericalProfile1D thlp2,
                          NumericalProfile1D rtpthlp,                              // In
                          NumericalProfile1D wp2, NumericalProfile1D wp3,          // In
                          NumericalProfile1D wp2thvp, NumericalProfile1D wp2up,    // In
                          NumericalProfile1D rtpthvp, NumericalProfile1D thlpthvp, // In
                          NumericalProfile2D hydromet, NumericalProfile2D sclrm,   // In
                          NumericalProfile2D edsclrm);                             // In

// Native scalar conservation residual; also callable within a device kernel.
// integral_after: Vertically-integrated quantity after dt time  [units vary]
// integral_before: Vertically-integrated quantity before dt time [units vary]
// flux_top: Total flux at the top of the domain           [units vary]
// flux_sfc: Total flux at the bottom of the domain        [units vary]
// integral_forcing: Vertically-integrated forcing                 [units vary]
// dt: Timestep size                                 [s]
KOKKOS_INLINE_FUNCTION real_t calculate_spurious_source(real_t integral_after,
                                                        real_t integral_before,             // In
                                                        real_t flux_top, real_t flux_sfc,   // In
                                                        real_t integral_forcing, real_t dt) // In
{
  //
  // Description:
  //   Checks whether there is conservation within the column and returns any
  //   imbalance as spurious_source where spurious_source is defined negative
  //   for a spurious sink.
  //
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  real_t spurious_source; // Signed column conservation residual [units vary]

  //--------------------------------------------------------------------

  // ---- Begin Code ----

  spurious_source =
      (integral_after - integral_before) / dt + flux_top - flux_sfc - integral_forcing;

  return spurious_source;
}

// ngrdcol: Number of grid columns
// edsclr_dim: Number of eddy-diffusivity passive scalars
// sclr_idx: Indices of the thlm/rtm/CO2 passive scalars
// params: Including C1, nu1, nu2, etc.
// l_implemented: True if CLUBB is being implemented and run in a host model
// l_input_fields: Flag for whether LES input fields are being used
// clubb_config_flags: Runtime model configuration flags
// err_info: err_info struct containing err_code and err_header
void check_clubb_settings_api(int ngrdcol,                                       // In
                              int edsclr_dim,                                    // In
                              const sclr_idx_type& sclr_idx,                     // In
                              Real2D params,                                     // In
                              bool l_implemented,                                // In
                              bool l_input_fields,                               // In
                              const clubb_config_flags_type& clubb_config_flags, // In
                              err_info_type& err_info);                          // InOut

// Actual surface callers already stage scalar diagnostics; adapt their scalar
// vectors through the existing array owner once before the Native-shaped check.
template <class A, class B, class C>
void sfc_varnce_check(int sclr_dim, real_t wp2_sfc, real_t up2_sfc, real_t vp2_sfc,
                      real_t thlp2_sfc, real_t rtp2_sfc, real_t rtpthlp_sfc, A sclrp2_sfc,
                      B sclrprtp_sfc, C sclrpthlp_sfc, err_info_type& err_info) {
  auto variance               = to_execution_space(sclrp2_sfc);
  auto water_covariance       = to_execution_space(sclrprtp_sfc);
  auto temperature_covariance = to_execution_space(sclrpthlp_sfc);
  sfc_varnce_check(sclr_dim, wp2_sfc, up2_sfc, vp2_sfc, thlp2_sfc, rtp2_sfc, rtpthlp_sfc, variance,
                   water_covariance, temperature_covariance, err_info);
}
} // namespace clubb_kokkos
