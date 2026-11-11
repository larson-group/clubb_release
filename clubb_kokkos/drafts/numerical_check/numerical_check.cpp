// Check model fields and settings with the complete Native diagnostic contracts.
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
// New profile diagnostics accept strided execution-space Views, so actual
// per-column callers need no bulk data packing. Native stats/debug metadata
// stays at explicit host boundaries. All-column fatal status and diagnostic
// order are preserved. The four existing live routines are copied verbatim.
#include "numerical_check.hpp"
#include "error_code.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "f2k/procedure_bridge.hpp"
#include "f2k/driver_io.hpp"
#include "f2k_dependencies.hpp"
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
namespace clubb_kokkos {
static void check_negative(Real2D var, const std::string& varname, const std::string& operation,
                           err_info_type& err_info);
template <class View>
static void check_nan_1d(View var, const std::string& varname, const std::string& operation,
                         err_info_type& err_info);
template <class View>
static void check_nan_2d(View var, const std::string& varname, const std::string& operation,
                         err_info_type& err_info);

// C++ dispatch expresses the native private generic check_nan interface.
template <class View>
static void check_nan(View var, const std::string& varname, const std::string& operation,
                      err_info_type& err_info) {
  // Select the native one-/two-dimensional finite-value check by View rank.
  static_assert(View::rank == 1 || View::rank == 2);
  if constexpr (View::rank == 1)
    check_nan_1d(var, varname, operation, err_info);
  else
    check_nan_2d(var, varname, operation, err_info);
}

static void check_nan_sclr(real_t var, const std::string& varname, const std::string& operation,
                           err_info_type& err_info);
static void check_nan(real_t var, const std::string& varname, const std::string& operation,
                      err_info_type& err_info) {
  // Native scalar generic dispatch; the existing View dispatcher is unchanged.
  check_nan_sclr(var, varname, operation, err_info);
}
template <class View> static bool any_nonfinite(View var) {
  // Device lowering of Native ANY(.NOT.IEEE_IS_FINITE(profile)); no host staging.
  static_assert(View::rank == 1);
  int invalid = 0;
  Kokkos::parallel_reduce(
      "invalid_model_arrays::finite", var.extent(0),
      KOKKOS_LAMBDA(const int level, int& count) { count += !Kokkos::isfinite(var(level)); },
      invalid);
  return invalid != 0;
}
template <class View> static std::string dump_values(View values) {
  // Complete Native array values are packed only at a host diagnostic boundary.
  const auto host = to_host_space(values);
  std::ostringstream out;
  out << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  for (std::size_t element = 0; element < host.size(); ++element)
    out << host.data()[element] << ' ';
  return out.str();
}
//=============================================================================
// Lscale: Mixing length                 [m]
// Lscale_up: Upward mixing length          [m]
// Lscale_down: Downward mixing length        [m]
// err_info: err_info struct containing err_code and err_header
void length_check(int nzt, NumericalProfile1D Lscale, NumericalProfile1D Lscale_up, // In
                  NumericalProfile1D Lscale_down,                                   // In
                  err_info_type& err_info)                                          // InOut
{
  //
  //        Description: This subroutine determines if any of the output
  //        variables for the length_new subroutine carry values that
  //        are NaNs.
  //
  //        Joshua Fasching February 2008
  //---------------------------------------------------------------------------------
  // Constant Parameters
  // Input Variables
  // Input/Output Variables
  const std::string proc_name = "compute_mixing_length";
  //-----------------------------------------------------------------------------

  check_nan(Lscale, "Lscale", proc_name,           // In
            err_info);                             // InOut
  check_nan(Lscale_up, "Lscale_up", proc_name,     // In
            err_info);                             // InOut
  check_nan(Lscale_down, "Lscale_down", proc_name, // In
            err_info);                             // InOut

  return;
}

//=============================================================================
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
    err_info_type& err_info)                                                        // InOut
{
  // Description: This subroutine determines if any of the output
  //   variables for the pdf_closure subroutine carry values that
  //   are NaNs.
  //
  // Joshua Fasching February 2008
  //---------------------------------------------------------------------------
  // Parameter Constants
  // Input Variables
  // Input (Optional passive scalar variables)
  // Input/Output variables
  // Local variables
  const std::string proc_name = "pdf_closure";
  int sclr; // Scalar loop index

  //-------------------------------------------------------------------------------

  // ---- Begin Code ----

  if (f2k::var_on_stats_list(stats, "wp4")) {
    check_nan(wp4, "wp4", proc_name, // In
              err_info);             // InOut
  }
  if (f2k::var_on_stats_list(stats, "wprtp2")) {
    check_nan(wprtp2, "wprtp2", proc_name, // In
              err_info);                   // InOut
  }
  check_nan(wp2rtp, "wp2rtp", proc_name, // In
            err_info);                   // InOut
  if (f2k::var_on_stats_list(stats, "wpthlp2")) {
    check_nan(wpthlp2, "wpthlp2", proc_name, // In
              err_info);                     // InOut
  }
  check_nan(wp2thlp, "wp2thlp", proc_name,       // In
            err_info);                           // InOut
  check_nan(cloud_frac, "cloud_frac", proc_name, // In
            err_info);                           // InOut
  check_nan(rcm, "rcm", proc_name,               // In
            err_info);                           // InOut
  check_nan(wpthvp, "wpthvp", proc_name,         // In
            err_info);                           // InOut
  check_nan(wp2thvp, "wp2thvp", proc_name,       // In
            err_info);                           // InOut
  check_nan(wp2up, "wp2up", proc_name,           // In
            err_info);                           // InOut
  check_nan(rtpthvp, "rtpthvp", proc_name,       // In
            err_info);                           // InOut
  check_nan(thlpthvp, "thlpthvp", proc_name,     // In
            err_info);                           // InOut
  check_nan(wprcp, "wprcp", proc_name,           // In
            err_info);                           // InOut
  check_nan(wp2rcp, "wp2rcp", proc_name,         // In
            err_info);                           // InOut
  check_nan(rtprcp, "rtprcp", proc_name,         // In
            err_info);                           // InOut
  check_nan(thlprcp, "thlprcp", proc_name,       // In
            err_info);                           // InOut
  if (f2k::var_on_stats_list(stats, "rcp2")) {
    check_nan(rcp2, "rcp2", proc_name, // In
              err_info);               // InOut
  }
  if (f2k::var_on_stats_list(stats, "wprtpthlp")) {
    check_nan(wprtpthlp, "wprtpthlp", proc_name, // In
              err_info);                         // InOut
  }
  check_nan(crt_1, "crt_1", proc_name,   // In
            err_info);                   // InOut
  check_nan(crt_2, "crt_2", proc_name,   // In
            err_info);                   // InOut
  check_nan(cthl_1, "cthl_1", proc_name, // In
            err_info);                   // InOut
  check_nan(cthl_2, "cthl_2", proc_name, // In
            err_info);                   // InOut
  // Check each PDF parameter at the grid level sent in.
  check_nan(Kokkos::subview(pdf_params.w_1, (1 - 1), Kokkos::ALL), "pdf_params%w_1(1,:)",   // In
            proc_name,                                                                      // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.w_2, (1 - 1), Kokkos::ALL), "pdf_params%w_2(1,:)",   // In
            proc_name,                                                                      // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.varnce_w_1, (1 - 1), Kokkos::ALL),                   // In
            "pdf_params%varnce_w_1(1,:)", proc_name,                                        // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.varnce_w_2, (1 - 1), Kokkos::ALL),                   // In
            "pdf_params%varnce_w_2(1,:)", proc_name,                                        // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.rt_1, (1 - 1), Kokkos::ALL), "pdf_params%rt_1(1,:)", // In
            proc_name,                                                                      // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.rt_2, (1 - 1), Kokkos::ALL), "pdf_params%rt_2(1,:)", // In
            proc_name,                                                                      // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.varnce_rt_1, (1 - 1), Kokkos::ALL),                  // In
            "pdf_params%varnce_rt_1(1,:)", proc_name,                                       // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.varnce_rt_2, (1 - 1), Kokkos::ALL),                  // In
            "pdf_params%varnce_rt_2(1,:)", proc_name,                                       // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.thl_1, (1 - 1), Kokkos::ALL), "pdf_params%thl_1(1,:)", // In
            proc_name,                                                                        // In
            err_info); // InOut
  check_nan(Kokkos::subview(pdf_params.thl_2, (1 - 1), Kokkos::ALL), "pdf_params%thl_2(1,:)", // In
            proc_name,                                                                        // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.varnce_thl_1, (1 - 1), Kokkos::ALL),                 // In
            "pdf_params%varnce_thl_1(1,:)", proc_name,                                      // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.varnce_thl_2, (1 - 1), Kokkos::ALL),                 // In
            "pdf_params%varnce_thl_2(1,:)", proc_name,                                      // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.mixt_frac, (1 - 1), Kokkos::ALL),                    // In
            "pdf_params%mixt_frac(1,:)", proc_name,                                         // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.corr_w_rt_1, (1 - 1), Kokkos::ALL),                  // In
            "pdf_params%corr_w_rt_1(1,:)", proc_name,                                       // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.corr_w_rt_2, (1 - 1), Kokkos::ALL),                  // In
            "pdf_params%corr_w_rt_2(1,:)", proc_name,                                       // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.corr_w_thl_1, (1 - 1), Kokkos::ALL),                 // In
            "pdf_params%corr_w_thl_1(1,:)", proc_name,                                      // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.corr_w_thl_2, (1 - 1), Kokkos::ALL),                 // In
            "pdf_params%corr_w_thl_2(1,:)", proc_name,                                      // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.corr_rt_thl_1, (1 - 1), Kokkos::ALL),                // In
            "pdf_params%corr_rt_thl_1(1,:)", proc_name,                                     // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.corr_rt_thl_2, (1 - 1), Kokkos::ALL),                // In
            "pdf_params%corr_rt_thl_2(1,:)", proc_name,                                     // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.rc_1, (1 - 1), Kokkos::ALL), "pdf_params%rc_1(1,:)", // In
            proc_name,                                                                      // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.rc_2, (1 - 1), Kokkos::ALL), "pdf_params%rc_2(1,:)", // In
            proc_name,                                                                      // In
            err_info);                                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.rsatl_1, (1 - 1), Kokkos::ALL),
            "pdf_params%rsatl_1(1,:)", // In
            proc_name,                 // In
            err_info);                 // InOut
  check_nan(Kokkos::subview(pdf_params.rsatl_2, (1 - 1), Kokkos::ALL),
            "pdf_params%rsatl_2(1,:)",                                      // In
            proc_name,                                                      // In
            err_info);                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.cloud_frac_1, (1 - 1), Kokkos::ALL), // In
            "pdf_params%cloud_frac_1(1,:)", proc_name,                      // In
            err_info);                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.cloud_frac_2, (1 - 1), Kokkos::ALL), // In
            "pdf_params%cloud_frac_2(1,:)", proc_name,                      // In
            err_info);                                                      // InOut
  check_nan(Kokkos::subview(pdf_params.chi_1, (1 - 1), Kokkos::ALL), "pdf_params%chi_1(1,:)", // In
            proc_name,                                                                        // In
            err_info); // InOut
  check_nan(Kokkos::subview(pdf_params.chi_2, (1 - 1), Kokkos::ALL), "pdf_params%chi_2(1,:)", // In
            proc_name,                                                                        // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.stdev_chi_1, (1 - 1), Kokkos::ALL),         // In
            "pdf_params%stdev_chi_1(1,:)", proc_name,                              // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.stdev_chi_2, (1 - 1), Kokkos::ALL),         // In
            "pdf_params%stdev_chi_2(1,:)", proc_name,                              // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.stdev_eta_1, (1 - 1), Kokkos::ALL),         // In
            "pdf_params%stdev_eta_1(1,:)", proc_name,                              // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.stdev_eta_2, (1 - 1), Kokkos::ALL),         // In
            "pdf_params%stdev_eta_2(1,:)", proc_name,                              // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.covar_chi_eta_1, (1 - 1), Kokkos::ALL),     // In
            "pdf_params%covar_chi_eta_1(1,:)", proc_name,                          // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.covar_chi_eta_2, (1 - 1), Kokkos::ALL),     // In
            "pdf_params%covar_chi_eta_2(1,:)", proc_name,                          // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.corr_w_chi_1, (1 - 1), Kokkos::ALL),        // In
            "pdf_params%corr_w_chi_1(1,:)", proc_name,                             // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.corr_w_chi_2, (1 - 1), Kokkos::ALL),        // In
            "pdf_params%corr_w_chi_2(1,:)", proc_name,                             // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.corr_w_eta_1, (1 - 1), Kokkos::ALL),        // In
            "pdf_params%corr_w_eta_1(1,:)", proc_name,                             // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.corr_w_eta_2, (1 - 1), Kokkos::ALL),        // In
            "pdf_params%corr_w_eta_2(1,:)", proc_name,                             // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.corr_chi_eta_1, (1 - 1), Kokkos::ALL),      // In
            "pdf_params%corr_chi_eta_1(1,:)", proc_name,                           // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.corr_chi_eta_2, (1 - 1), Kokkos::ALL),      // In
            "pdf_params%corr_chi_eta_2(1,:)", proc_name,                           // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.alpha_thl, (1 - 1), Kokkos::ALL),           // In
            "pdf_params%alpha_thl(1,:)", proc_name,                                // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.alpha_rt, (1 - 1), Kokkos::ALL),            // In
            "pdf_params%alpha_rt(1,:)", proc_name,                                 // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.ice_supersat_frac_1, (1 - 1), Kokkos::ALL), // In
            "pdf_params%ice_supersat_frac_1(1,:)", proc_name,                      // In
            err_info);                                                             // InOut
  check_nan(Kokkos::subview(pdf_params.ice_supersat_frac_2, (1 - 1), Kokkos::ALL), // In
            "pdf_params%ice_supersat_frac_2(1,:)", proc_name,                      // In
            err_info);                                                             // InOut

  if (sclr_dim > 0) {
    for (sclr = 1; sclr <= sclr_dim; sclr += 1) {
      check_nan(Kokkos::subview(sclrpthvp, Kokkos::ALL, (sclr - 1)), "sclrpthvp", proc_name, // In
                err_info);                                                                 // InOut
      check_nan(Kokkos::subview(sclrprcp, Kokkos::ALL, (sclr - 1)), "sclrprcp", proc_name, // In
                err_info);                                                                 // InOut
      check_nan(Kokkos::subview(wpsclrprtp, Kokkos::ALL, (sclr - 1)), "wpsclrprtp", proc_name, // In
                err_info);                                                                 // InOut
      check_nan(Kokkos::subview(wpsclrp2, Kokkos::ALL, (sclr - 1)), "wpsclrp2", proc_name, // In
                err_info);                                                                 // InOut
      check_nan(Kokkos::subview(wpsclrpthlp, Kokkos::ALL, (sclr - 1)), "wpsclrtlp", proc_name, // In
                err_info);                                                                 // InOut
      check_nan(Kokkos::subview(wp2sclrp, Kokkos::ALL, (sclr - 1)), "wp2sclrp", proc_name, // In
                err_info);                                                                 // InOut
      // sclr = 1, sclr_dim, 1
    }
  }

  return;
}

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
    err_info_type& err_info) {                                                    // InOut
  // Identify non-finite core fields, then check positive-definite inputs.
  // Preserve the native NaN early return before the beginning-of-step reset.
  // Joshua Fasching,February2008.
  // Native character(len=18) holds this procedure name without extra padding.
  const std::string proc_name = "advance_clubb_core";

  // -------- Input Nan Check ----------------------------------------------
  // Reductions inspect execution-space fields; host copies are needed only to
  // print bad coordinates, instead of unconditionally staging every input.
  check_nan(thlm_forcing, "thlm_forcing", prefix + proc_name, err_info);
  check_nan(rtm_forcing, "rtm_forcing", prefix + proc_name, err_info);
  check_nan(um_forcing, "um_forcing", prefix + proc_name, err_info);
  check_nan(vm_forcing, "vm_forcing", prefix + proc_name, err_info);

  check_nan(wm_zm, "wm_zm", prefix + proc_name, err_info);
  check_nan(wm_zt, "wm_zt", prefix + proc_name, err_info);
  check_nan(p_in_Pa, "p_in_Pa", prefix + proc_name, err_info);
  check_nan(rho_zm, "rho_zm", prefix + proc_name, err_info);
  check_nan(rho, "rho", prefix + proc_name, err_info);
  check_nan(exner, "exner", prefix + proc_name, err_info);
  check_nan(rho_ds_zm, "rho_ds_zm", prefix + proc_name, err_info);
  check_nan(rho_ds_zt, "rho_ds_zt", prefix + proc_name, err_info);
  check_nan(invrs_rho_ds_zm, "invrs_rho_ds_zm", prefix + proc_name, err_info);
  check_nan(invrs_rho_ds_zt, "invrs_rho_ds_zt", prefix + proc_name, err_info);
  check_nan(thv_ds_zm, "thv_ds_zm", prefix + proc_name, err_info);
  check_nan(thv_ds_zt, "thv_ds_zt", prefix + proc_name, err_info);

  check_nan(um, "um", prefix + proc_name, err_info);
  check_nan(upwp, "upwp", prefix + proc_name, err_info);
  check_nan(vm, "vm", prefix + proc_name, err_info);
  check_nan(vpwp, "vpwp", prefix + proc_name, err_info);
  check_nan(up2, "up2", prefix + proc_name, err_info);
  check_nan(vp2, "vp2", prefix + proc_name, err_info);
  check_nan(rtm, "rtm", prefix + proc_name, err_info);
  check_nan(wprtp, "wprtp", prefix + proc_name, err_info);
  check_nan(thlm, "thlm", prefix + proc_name, err_info);
  check_nan(wpthlp, "wpthlp", prefix + proc_name, err_info);
  check_nan(wp2, "wp2", prefix + proc_name, err_info);
  check_nan(wp3, "wp3", prefix + proc_name, err_info);
  check_nan(rtp2, "rtp2", prefix + proc_name, err_info);
  check_nan(thlp2, "thlp2", prefix + proc_name, err_info);
  check_nan(rtpthlp, "rtpthlp", prefix + proc_name, err_info);

  check_nan(wpthlp_sfc, "wpthlp_sfc", prefix + proc_name, err_info);
  check_nan(wprtp_sfc, "wprtp_sfc", prefix + proc_name, err_info);
  check_nan(upwp_sfc, "upwp_sfc", prefix + proc_name, err_info);
  check_nan(vpwp_sfc, "vpwp_sfc", prefix + proc_name, err_info);
  check_nan(p_sfc, "p_sfc", prefix + proc_name, err_info);

  for (int sclr = 0; sclr < sclr_dim; ++sclr) {
    check_nan(Kokkos::subview(sclrm_forcing, Kokkos::ALL, Kokkos::ALL, sclr), "sclrm_forcing",
              prefix + proc_name, err_info);
    check_nan(Kokkos::subview(wpsclrp_sfc, Kokkos::ALL, sclr), "wpsclrp_sfc", prefix + proc_name,
              err_info);
    check_nan(Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr), "sclrm", prefix + proc_name,
              err_info);
    check_nan(Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr), "wpsclrp",
              prefix + proc_name, err_info);
    check_nan(Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, sclr), "sclrp2", prefix + proc_name,
              err_info);
    check_nan(Kokkos::subview(sclrprtp, Kokkos::ALL, Kokkos::ALL, sclr), "sclrprtp",
              prefix + proc_name, err_info);
    check_nan(Kokkos::subview(sclrpthlp, Kokkos::ALL, Kokkos::ALL, sclr), "sclrpthlp",
              prefix + proc_name, err_info);
  }

  for (int edsclr = 0; edsclr < edsclr_dim; ++edsclr) {
    check_nan(Kokkos::subview(edsclrm_forcing, Kokkos::ALL, Kokkos::ALL, edsclr), "edsclrm_forcing",
              prefix + proc_name, err_info);
    check_nan(Kokkos::subview(wpedsclrp_sfc, Kokkos::ALL, edsclr), "wpedsclrp_sfc",
              prefix + proc_name, err_info);
    check_nan(Kokkos::subview(edsclrm, Kokkos::ALL, Kokkos::ALL, edsclr), "edsclrm",
              prefix + proc_name, err_info);
  }

  // ---------------------------------------------------------------------
  if (f2k::clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error))
      return;
  }

  check_negative(rtm, "rtm", prefix + proc_name, err_info);
  check_negative(p_in_Pa, "p_in_Pa", prefix + proc_name, err_info);
  check_negative(rho, "rho", prefix + proc_name, err_info);
  check_negative(rho_zm, "rho_zm", prefix + proc_name, err_info);
  check_negative(exner, "exner", prefix + proc_name, err_info);
  check_negative(rho_ds_zm, "rho_ds_zm", prefix + proc_name, err_info);
  check_negative(rho_ds_zt, "rho_ds_zt", prefix + proc_name, err_info);
  check_negative(invrs_rho_ds_zm, "invrs_rho_ds_zm", prefix + proc_name, err_info);
  check_negative(invrs_rho_ds_zt, "invrs_rho_ds_zt", prefix + proc_name, err_info);
  check_negative(thv_ds_zm, "thv_ds_zm", prefix + proc_name, err_info);
  check_negative(thv_ds_zt, "thv_ds_zt", prefix + proc_name, err_info);
  check_negative(up2, "up2", prefix + proc_name, err_info);
  check_negative(vp2, "vp2", prefix + proc_name, err_info);
  check_negative(wp2, "wp2", prefix + proc_name, err_info);
  check_negative(thlm, "thlm", prefix + proc_name, err_info);
  check_negative(rtp2, "rtp2", prefix + proc_name, err_info);
  check_negative(thlp2, "thlp2", prefix + proc_name, err_info);

  // Fortran character comparisons blank-pad their operands.
  const auto prefix_last = prefix.find_last_not_of(' ');
  const bool beginning   = prefix.compare(0, prefix_last == std::string::npos ? 0 : prefix_last + 1,
                                          "beginning of") == 0;
  if (any_equal(err_info.err_code, clubb_fatal_error) && beginning) {
    // Negative value generated by the host model, hence ignore that error.
    Kokkos::deep_copy(err_info.err_code, clubb_no_error);
  }

  // Check the first levels for cold temperatures (the native threshold is190K).
  const int first_levels = Kokkos::min(10, nzt);
  int cold               = 0;
  Kokkos::parallel_reduce(
      "parameterization_check::cold_levels", Kokkos::RangePolicy<>(0, ngrdcol * first_levels),
      KOKKOS_LAMBDA(const int cell, int& count) {
        count += thlm(cell % ngrdcol, cell / ngrdcol) < real_t(190.0);
      },
      cold);
  if (cold) {
    const auto host = to_host_space(thlm);
    for (int i = 0; i < ngrdcol; ++i)
      for (int k = 0; k < first_levels; ++k)
        if (host(i, k) < real_t(190.0)) {
          // Keep enough digits to distinguish a value just below190K without
          // changing the formatting of later messages on the shared stream.
          std::ostringstream message;
          message << "Liquid water potential temperature (thlm) < 190K "
                  << "at grid column i = " << i + 1 << " and grid level k = " << k + 1 << ": thlm("
                  << i + 1 << "," << k + 1
                  << ") = " << std::setprecision(std::numeric_limits<real_t>::max_digits10)
                  << host(i, k) << '\n';
          std::cerr << message.str();
        }
  }
}
// clang-format on

//=============================================================================
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
                      err_info_type& err_info)           // InOut
{
  //
  //       Description:This subroutine determines if any of the output
  //       variables for the calc_surface_varnce subroutine carry values that
  //       are nans.
  //
  //       Joshua Fasching February 2008
  //
  //
  //-----------------------------------------------------------------------
  // Constant Parameters
  // Name of the subroutine calling the check
  // Input Variables
  // Input/Output variables
  const std::string proc_name = "calc_surface_varnce";
  //-----------------------------------------------------------------------

  // ---- Begin Code ----

  check_nan(wp2_sfc, "wp2_sfc", proc_name,         // In
            err_info);                             // InOut
  check_nan(up2_sfc, "up2_sfc", proc_name,         // In
            err_info);                             // InOut
  check_nan(vp2_sfc, "vp2_sfc", proc_name,         // In
            err_info);                             // InOut
  check_nan(thlp2_sfc, "thlp2_sfc", proc_name,     // In
            err_info);                             // InOut
  check_nan(rtp2_sfc, "rtp2_sfc", proc_name,       // In
            err_info);                             // InOut
  check_nan(rtpthlp_sfc, "rtpthlp_sfc", proc_name, // In
            err_info);                             // InOut

  if (sclr_dim > 0) {
    check_nan(sclrp2_sfc, "sclrp2_sfc", proc_name, // In
              err_info);                           // InOut

    check_nan(sclrprtp_sfc, "sclrprtp_sfc", proc_name, // In
              err_info);                               // InOut

    check_nan(sclrpthlp_sfc, "sclrpthlp_sfc", proc_name, // In
              err_info);                                 // InOut
  }

  return;
}

//=============================================================================
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
               err_info_type& err_info) // InOut
{
  // Description:
  //   Checks radiation input variables. If they are < 0 it reports
  //   to the console.
  //------------------------------------------------------------------------
  // Constant Parameters
  // Input/Output variables
  // Input/Output variables
  // Local variables
  const std::string proc_name = "Before BUGSrad.";
  Kokkos::View<real_t**, Kokkos::LayoutLeft> rvm(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "rad_check::rvm"), ngrdcol,
      nzt); // Water vapor mixing ratio r_t-r_c [kg/kg]

  //-------------------------------------------------------------------------

  Kokkos::parallel_for(
      "rad_check::array1", ((rvm.extent(0) - (1)) / (1) + 1) * ((rvm.extent(1) - (1)) / (1) + 1),
      KOKKOS_LAMBDA(int cell) {
        int remaining     = cell;
        const int slice_0 = remaining % ((rvm.extent(0) - (1)) / (1) + 1);
        remaining /= ((rvm.extent(0) - (1)) / (1) + 1);
        const int slice_1 = remaining % ((rvm.extent(1) - (1)) / (1) + 1);
        remaining /= ((rvm.extent(1) - (1)) / (1) + 1);
        rvm(slice_0, slice_1) = rtm(slice_0, slice_1) - rcm(slice_0, slice_1);
      });

  check_negative(thlm, "thlm", proc_name,             // In
                 err_info);                           // InOut
  check_negative(rcm, "rcm", proc_name,               // In
                 err_info);                           // InOut
  check_negative(rtm, "rtm", proc_name,               // In
                 err_info);                           // InOut
  check_negative(rvm, "rvm", proc_name,               // In
                 err_info);                           // InOut
  check_negative(rim, "rim", proc_name,               // In
                 err_info);                           // InOut
  check_negative(cloud_frac, "cloud_frac", proc_name, // In
                 err_info);                           // InOut
  check_negative(p_in_Pa, "p_in_Pa", proc_name,       // In
                 err_info);                           // InOut
  check_negative(exner, "exner", proc_name,           // In
                 err_info);                           // InOut
  check_negative(rho_zm, "rho_zm", proc_name,         // In
                 err_info);                           // InOut

  return;
}

//=============================================================================
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
                          NumericalProfile2D edsclrm)                              // In
{
  //       Description:
  //       Checks for invalid floating point values in select model arrays.
  //       References:
  //       None
  //------------------------------------------------------------------------
  // Local Variables
  bool invalid_model_arrays;
  int i;      // Hydrometeor loop index
  int sclr;   // Passive scalar loop index
  int edsclr; // Eddy-diffusivity passive scalar loop index

  invalid_model_arrays = false;

  // Check whether any variable array contains a NaN for
  // um, vm, thlm, rtm, rtp2, thlp2, wprtp, wpthlp, rtpthlp,
  // wp2, & wp3.
  if (any_nonfinite(um)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in um model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "um= ", um
    invalid_model_arrays = true;
  }

  if (any_nonfinite(vm)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in vm model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "vm= ", vm
    invalid_model_arrays = true;
  }

  if (any_nonfinite(wp2)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in wp2 model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "wp2= ", wp2
    invalid_model_arrays = true;
  }

  if (any_nonfinite(wp3)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in wp3 model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "wp3= ", wp3
    invalid_model_arrays = true;
  }

  if (any_nonfinite(rtm)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in rtm model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "rtm= ", rtm
    invalid_model_arrays = true;
  }

  if (any_nonfinite(thlm)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in thlm model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "thlm= ", thlm
    invalid_model_arrays = true;
  }

  if (any_nonfinite(rtp2)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in rtp2 model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "rtp2= ", rtp2
    invalid_model_arrays = true;
  }

  if (any_nonfinite(thlp2)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in thlp2 model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "thlp2= ", thlp2
    invalid_model_arrays = true;
  }

  if (any_nonfinite(wprtp)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in wprtp model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "wprtp= ", wprtp
    invalid_model_arrays = true;
  }

  if (any_nonfinite(wpthlp)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in wpthlp model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "wpthlp= ", wpthlp
    invalid_model_arrays = true;
  }

  if (any_nonfinite(rtpthlp)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in rtpthlp model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "rtpthlp= ", rtpthlp
    invalid_model_arrays = true;
  }

  if (hydromet_dim > 0) {
    for (i = 1; i <= hydromet_dim; i += 1) {
      if (any_nonfinite(Kokkos::subview(hydromet, Kokkos::ALL, (i - 1)))) {
        { // Native character(len=10) list is stored as (character,species).
          const auto names_host = to_host_space(hydromet_list);
          std::string name;
          for (int character = 0; character < 10; ++character)
            name += names_host(character, i - 1);
          const auto last = name.find_last_not_of(' ');
          name.resize(last == std::string::npos ? 0 : last + 1);
          std::cerr << "NaN in a hydrometeor model array " + name << '\n';
        }
        //             write(fstderr,*) "hydromet= ", hydromet
        invalid_model_arrays = true;
        //             return
      }
    }
  }

  //       if ( any( .not. ieee_is_finite( wm_zt ) ) ) then
  //         write(fstderr,*) "NaN in wm_zt model array"
  //         write(fstderr,*) "wm_zt= ", wm_zt
  //         invalid_model_arrays = .true.
  //       end if

  if (any_nonfinite(wp2thvp)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in wp2thvp model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "wp2thvp = ", wp2thvp
    invalid_model_arrays = true;
  }

  if (any_nonfinite(wp2up)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in wp2up model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "wp2up = ", wp2up
    invalid_model_arrays = true;
  }

  if (any_nonfinite(rtpthvp)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in rtpthvp model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "rtpthvp = ", rtpthvp
    invalid_model_arrays = true;
  }

  if (any_nonfinite(thlpthvp)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "NaN in thlpthvp model array" << "\n";
      std::cerr << diagnostic_message.str();
    }
    //         write(fstderr,*) "thlpthvp = ", thlpthvp
    invalid_model_arrays = true;
  }

  for (sclr = 1; sclr <= sclr_dim; sclr += 1) {
    if (any_nonfinite(Kokkos::subview(sclrm, Kokkos::ALL, (sclr - 1)))) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "NaN in sclrm" << " " << sclr << " " << "model array" << "\n";
        std::cerr << diagnostic_message.str();
      }
      //           write(fstderr,'(a6,i2,a1)') "sclrm(", sclr, ")"
      //           write(fstderr,*) sclrm(:,sclr)
      invalid_model_arrays = true;
    }
  }

  for (edsclr = 1; edsclr <= edsclr_dim; edsclr += 1) {
    if (any_nonfinite(Kokkos::subview(edsclrm, Kokkos::ALL, (edsclr - 1)))) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "NaN in edsclrm" << " " << edsclr << " " << "model array" << "\n";
        std::cerr << diagnostic_message.str();
      }
      //           write(fstderr,'(a8,i2,a1)') "edsclrm(", edsclr, ")"
      //           write(fstderr,*) edsclrm(:,edsclr)
      invalid_model_arrays = true;
    }
  }

  return invalid_model_arrays;
}

// clang-format off
static void check_negative(Real2D var, const std::string& varname, const std::string& operation,
                           err_info_type& err_info) {
  // Report every negative value in native level/column order and mark all
  // columns fatal, as in the source's general-error assignment.
  const int ngrdcol = var.extent(0);
  int negative      = 0;
  Kokkos::parallel_reduce(
      "check_negative", Kokkos::RangePolicy<>(0, var.size()),
      KOKKOS_LAMBDA(const int cell, int& count) {
        count += var(cell % ngrdcol, cell / ngrdcol) < real_t(0.0);
      },
      negative);
  if (negative) {
    const auto host = to_host_space(var);
    for (int k = 0; k < host.extent_int(1); ++k)
      for (int i = 0; i < host.extent_int(0); ++i)
        if (host(i, k) < real_t(0.0)) {
          std::cerr << err_info.err_header_global << '\n';
          std::cerr << varname << " < 0 in " << operation << " at i = " << i + 1
                    << " and k = " << k + 1 << '\n';
        }
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
  }
}
// clang-format on

// clang-format off
template <class View>
static void check_nan_1d(View var, const std::string& varname, const std::string& operation,
                         err_info_type& err_info) {
  // Report a non-finite profile/column field and mark all columns fatal.
  int invalid = 0;
  Kokkos::parallel_reduce(
      "check_nan_1d", Kokkos::RangePolicy<>(0, var.extent(0)),
      KOKKOS_LAMBDA(const int i, int& count) { count += !Kokkos::isfinite(var(i)); }, invalid);
  if (invalid) {
    std::cerr << err_info.err_header_global << '\n';
    std::cerr << varname << " is NaN in " << operation << '\n';
    // Single-column error, but the column number is not available.
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
  }
}
// clang-format on

// clang-format off
template <class View>
static void check_nan_2d(View var, const std::string& varname, const std::string& operation,
                         err_info_type& err_info) {
  // Report a non-finite column/profile field and mark all columns fatal.
  const int ngrdcol = var.extent(0);
  int invalid       = 0;
  Kokkos::parallel_reduce(
      "check_nan_2d", Kokkos::RangePolicy<>(0, var.size()),
      KOKKOS_LAMBDA(const int cell, int& count) {
        count += !Kokkos::isfinite(var(cell % ngrdcol, cell / ngrdcol));
      },
      invalid);
  if (invalid) {
    std::cerr << err_info.err_header_global << '\n';
    std::cerr << varname << " is NaN in " << operation << '\n';
    // Single-column error, but the column number is not available.
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
  }
}
// clang-format on

//=============================================================================
// var: Variable being examined
// varname: Name of variable being examined
// operation: Procedure calling check_nan
// err_info: err_info struct containing err_code and err_header
static void check_nan_sclr(real_t var, const std::string& varname,
                           const std::string& operation, // In
                           err_info_type& err_info)      // InOut
{
  //
  // Description:
  //   Checks for a NaN in the scalar var then reports it.
  //
  //-----------------------------------------------------------------------
  // External
  // Input Variables
  // Input/Output variables

  //--------------------------------------------------------------------
  if (!Kokkos::isfinite(var)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << varname << " " << " is NaN in " << " " << operation << "\n";
      std::cerr << diagnostic_message.str();
    }
    // Single col error, but col nr is not available
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
  }

  return;
}

// Native calculate_spurious_source is header-defined for its device scalar API.

//=============================================================================
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
                              err_info_type& err_info)                           // InOut
{
  // Description:
  //   Subroutine to set up the model for execution.
  //
  // References:
  //   None
  //---------------------------------------------------------------------
  //---------------------- Input Variables ----------------------
  // Flags
  //---------------------- InOut Variables ----------------------

  //---------------------- Begin Code ----------------------

  f2k::initialize_error_headers();

  // This is the complete configuration-validation entry point.  The
  // parameter catalog and parameter-only relations belong to
  // parameters_tunable; flag and mixed consistency rules remain below.
  f2k::check_parameters_api(ngrdcol, params, // In
                            err_info);       // InOut
  // Reuse the existing Native metadata owner for source module sponge flags.
  const auto environment = f2k::driver_environment();

#ifdef CLUBB_GPU
  if (clubb_config_flags.penta_solve_method == lapack) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "WARNING: The penta-diagonal lapack solver is not GPU accelerated"
                         << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message
          << " Set penta_solve_method = 2, to use an accelerated penta-diagonal solver" << "\n";
      std::cerr << diagnostic_message.str();
    }
  }

  if (clubb_config_flags.tridiag_solve_method == lapack) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "WARNING: The tri-diagonal lapack solver is not GPU accelerated"
                         << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message
          << " Set tridiag_solve_method = 2, to use an accelerated tri-diagonal solver" << "\n";
      std::cerr << diagnostic_message.str();
    }
  }

  if (l_input_fields) {
    throw std::runtime_error("l_input_fields = .true. not usable when running on GPUs");
  }
#endif

  // Sanity check
  if (f2k::clubb_at_least_debug_level_api(0)) {

    int any_parameter_relation_1 = 0;
    Kokkos::parallel_reduce(
        "check_clubb_settings_api::relation1", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& hits) {
          if (Kokkos::abs(params((1 - 1 + column * (1)), iC1) -
                          params((1 - 1 + column * (1)), iC14)) >
              Kokkos::abs(params((1 - 1 + column * (1)), iC1) +
                          params((1 - 1 + column * (1)), iC14)) /
                  2 * eps)
            ++hits;
        },
        any_parameter_relation_1);
    if (clubb_config_flags.l_damp_wp2_using_em &&
        ((any_parameter_relation_1 != 0) || clubb_config_flags.l_stability_correct_tau_zm)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message
            << "l_damp_wp2_using_em = T requires C1=C14 and l_stability_correct_tau_zm = F" << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C1 = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC1)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C14 = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC14)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "l_stability_correct_tau_zm = " << " "
                           << bool(clubb_config_flags.l_stability_correct_tau_zm) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }

  // Sanity check for the saturation formula
  switch (clubb_config_flags.saturation_formula) {
  case saturation_bolton: {
    // Using the Bolton 1980 approximations for SVP over vapor/ice

    break;
  }
  case saturation_gfdl: {
    // Using the Flatau, et al. polynomial approximation for SVP over vapor/ice

    break;
  }
  case saturation_flatau: {
    // h1g, 2010-06-16
    // Using the GFDL SVP formula (Goff-Gratch)

    // Add new saturation formulas after this

    break;
  }
  case saturation_lookup: {
    // Using the lookup table

    break;
  }
  default: {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Unknown approx. of saturation vapor pressure: " << " "
                         << clubb_config_flags.saturation_formula << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
    break;
  }
  }

  // Check for the type of two component normal (double Gaussian) PDF being
  // used for w, rt, and theta-l (or w, chi, and eta).
  if (clubb_config_flags.iiPDF_type < iiPDF_ADG1 ||
      clubb_config_flags.iiPDF_type > iiPDF_new_hybrid) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Unknown type of double Gaussian PDF selected: " << " "
                         << clubb_config_flags.iiPDF_type << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "iiPDF_type = " << " " << clubb_config_flags.iiPDF_type << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
    // iiPDF_type < iiPDF_ADG1 or iiPDF_type > iiPDF_lY93
  }

  // The ADG2 and 3D Luhar PDFs can only be used as part of input fields.
  if (clubb_config_flags.iiPDF_type == iiPDF_ADG2) {
    if (!l_input_fields) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message
            << "The ADG2 PDF can only be used with input fields (l_input_fields = .true.)." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "iiPDF_type = " << " " << clubb_config_flags.iiPDF_type << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "l_input_fields = " << " " << l_input_fields << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
      // .not. l_input_fields
    }
    // iiPDF_type == iiPDF_ADG2
  }

  if (clubb_config_flags.iiPDF_type == iiPDF_3D_Luhar) {
    if (!l_input_fields) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message
            << "The 3D Luhar PDF can only be used with input fields (l_input_fields = .true.)."
            << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "iiPDF_type = " << " " << clubb_config_flags.iiPDF_type << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "l_input_fields = " << " " << l_input_fields << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
      // .not. l_input_fields
    }
    // iiPDF_type == iiPDF_3D_Luhar
  }

  // This also currently applies to the new PDF until it has been fully
  // implemented.
  if (clubb_config_flags.iiPDF_type == iiPDF_new) {
    if (!l_input_fields) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message
            << "The new PDF can only be used with input fields (l_input_fields = .true.)." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "iiPDF_type = " << " " << clubb_config_flags.iiPDF_type << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "l_input_fields = " << " " << l_input_fields << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
      // .not. l_input_fields
    }
    // iiPDF_type == iiPDF_new
  }

  // This also currently applies to the TSDADG PDF until it has been fully
  // implemented.
  if (clubb_config_flags.iiPDF_type == iiPDF_TSDADG) {
    if (!l_input_fields) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message
            << "The new TSDADG PDF can only be used with input fields (l_input_fields = .true.)."
            << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "iiPDF_type = " << " " << clubb_config_flags.iiPDF_type << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "l_input_fields = " << " " << l_input_fields << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
      // .not. l_input_fields
    }
    // iiPDF_type == iiPDF_TSDADG
  }

  // This also applies to Lewellen and Yoh (1993).
  if (clubb_config_flags.iiPDF_type == iiPDF_LY93) {
    if (!l_input_fields) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "The Lewellen and Yoh PDF can only be used with input fields "
                              "(l_input_fields = .true.)."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "iiPDF_type = " << " " << clubb_config_flags.iiPDF_type << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "l_input_fields = " << " " << l_input_fields << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
      // .not. l_input_fields
    }
    // iiPDF_type == iiPDF_LY93
  }

  // Check the option for the placement of the call to CLUBB's PDF.
  if (clubb_config_flags.ipdf_call_placement < ipdf_pre_advance_fields ||
      clubb_config_flags.ipdf_call_placement > ipdf_pre_post_advance_fields) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Invalid option selected for ipdf_call_placement: " << " "
                         << clubb_config_flags.ipdf_call_placement << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }

  // The l_predict_upwp_vpwp flag requires that the ADG1 PDF is used
  // implicitly in subroutine advance_xm_wpxp.
  if (clubb_config_flags.l_predict_upwp_vpwp) {

    // When l_predict_upwp_vpwp is enabled, the
    // l_explicit_turbulent_adv_wpxp flag must be turned off.
    // Otherwise, explicit turbulent advection would require PDF parameters
    // for u and v to be calculated in PDF closure.  These would be needed
    // to calculate integrated fields such as wp2up, etc.
    if (l_explicit_turbulent_adv_wpxp) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "The l_explicit_turbulent_adv_wpxp option is not currently set up "
                              "for use with the l_predict_upwp_vpwp code."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
      // l_explicit_turbulent_adv_wpxp
    }

    // When l_predict_upwp_vpwp is enabled, the PDF type must be set to
    // the ADG1 PDF or the new hybrid PDF.  The other PDFs are not currently
    // set up to calculate variables needed for implicit or semi-implicit
    // turbulent advection, such as coef_wp2up_implicit, etc.
    if ((clubb_config_flags.iiPDF_type != iiPDF_ADG1) &&
        (clubb_config_flags.iiPDF_type != iiPDF_new_hybrid)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Currently, only the ADG1 PDF and the new hybrid PDF are set up for "
                              "use with the l_predict_upwp_vpwp code."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
      // iiPDF_type /= iiPDF_ADG1
    }

    // l_predict_upwp_vpwp
  }

  // The l_do_expldiff_rtm_thlm flag explicitly diffuses copies of thlm and
  // rtm that advance_windm_edsclrm stores in the eddy scalars
  // sclr_idx%iiedsclr_thl and sclr_idx%iiedsclr_rt.  With fewer than two
  // eddy scalars there is nowhere to store the copies, so the flag has no
  // effect.
  if (clubb_config_flags.l_do_expldiff_rtm_thlm && edsclr_dim > 1) {

    // The copies of thlm and rtm need distinct, valid eddy-scalar slots.
    if (sclr_idx.iiedsclr_thl < 1 || sclr_idx.iiedsclr_thl > edsclr_dim ||
        sclr_idx.iiedsclr_rt < 1 || sclr_idx.iiedsclr_rt > edsclr_dim ||
        sclr_idx.iiedsclr_thl == sclr_idx.iiedsclr_rt) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "l_do_expldiff_rtm_thlm = T requires distinct iiedsclr_thl and "
                              "iiedsclr_rt in 1:edsclr_dim."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "iiedsclr_thl = " << " " << sclr_idx.iiedsclr_thl << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "iiedsclr_rt = " << " " << sclr_idx.iiedsclr_rt << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "edsclr_dim = " << " " << edsclr_dim << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
      // invalid iiedsclr_thl or iiedsclr_rt
    }

    // l_do_expldiff_rtm_thlm
  }

  // The flags l_min_xp2_from_corr_wx and l_enable_relaxed_clipping must
  // have opposite values.
  if ((clubb_config_flags.l_min_xp2_from_corr_wx) &&
      (clubb_config_flags.l_enable_relaxed_clipping)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message
          << "Invalid configuration: l_min_xp2_from_corr_wx = T and l_enable_relaxed_clipping = T"
          << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "They must have opposite values" << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  } else if ((!clubb_config_flags.l_min_xp2_from_corr_wx) &&
             (!clubb_config_flags.l_enable_relaxed_clipping)) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message
          << "Invalid configuration: l_min_xp2_from_corr_wx = F and l_enable_relaxed_clipping = F"
          << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "They must have opposite values" << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    // err_info%err_code = clubb_fatal_error
    // return
  }

  // Checking for the code that orders CLUBB's advance_ subroutines
  if (order_xm_wpxp < 1 || order_xm_wpxp > 4) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "The variable order_xm_wpxp must have a value between 1 and 4" << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_xm_wpxp = " << " " << order_xm_wpxp << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  } else if (order_xm_wpxp == order_wp2_wp3 || order_xm_wpxp == order_xp2_xpyp ||
             order_xm_wpxp == order_windm) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "The variable order_xm_wpxp has the same value as another order_ "
                            "variable.  Please give each order index a unique value."
                         << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_xm_wpxp = " << " " << order_xm_wpxp << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_wp2_wp3 = " << " " << order_wp2_wp3 << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_xp2_xpyp = " << " " << order_xp2_xpyp << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_windm = " << " " << order_windm << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }

  if (order_wp2_wp3 < 1 || order_wp2_wp3 > 4) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "The variable order_wp2_wp3 must have a value between 1 and 4" << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_wp2_wp3 = " << " " << order_wp2_wp3 << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  } else if (order_wp2_wp3 == order_xm_wpxp || order_wp2_wp3 == order_xp2_xpyp ||
             order_wp2_wp3 == order_windm) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "The variable order_wp2_wp3 has the same value as another order_ "
                            "variable.  Please give each order index a unique value."
                         << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_wp2_wp3 = " << " " << order_wp2_wp3 << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_xm_wpxp = " << " " << order_xm_wpxp << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_xp2_xpyp = " << " " << order_xp2_xpyp << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_windm = " << " " << order_windm << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }

  if (order_xp2_xpyp < 1 || order_xp2_xpyp > 4) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "The variable order_xp2_xpyp must have a value between 1 and 4" << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_xp2_xpyp = " << " " << order_xp2_xpyp << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  } else if (order_xp2_xpyp == order_wp2_wp3 || order_xp2_xpyp == order_xm_wpxp ||
             order_xp2_xpyp == order_windm) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "The variable order_xp2_xpyp has the same value as another order_ "
                            "variable.  Please give each order index a unique value."
                         << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_xp2_xpyp = " << " " << order_xp2_xpyp << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_wp2_wp3 = " << " " << order_wp2_wp3 << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_xm_wpxp = " << " " << order_xm_wpxp << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_windm = " << " " << order_windm << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }

  if (order_windm < 1 || order_windm > 4) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "The variable order_windm must have a value between 1 and 4" << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_windm = " << " " << order_windm << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  } else if (order_windm == order_wp2_wp3 || order_windm == order_xp2_xpyp ||
             order_windm == order_xm_wpxp) {
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << err_info.err_header_global << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "The variable order_windm has the same value as another order_ "
                            "variable.  Please give each order index a unique value."
                         << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_windm = " << " " << order_windm << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_wp2_wp3 = " << " " << order_wp2_wp3 << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_xp2_xpyp = " << " " << order_xp2_xpyp << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "order_xm_wpxp = " << " " << order_xm_wpxp << "\n";
      std::cerr << diagnostic_message.str();
    }
    {
      std::ostringstream diagnostic_message;
      diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
      std::cerr << diagnostic_message.str();
    }
    // General error -> set all entries to clubb_fatal_error
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }

  // The tau-based method tunes C_invrs_tau parameters instead of these Cx values.
  // relevant Cx tunable parameters are all set to a value of 1 (as
  // you're supposed to tune the C_invrs_tau_ parameters instead).
  if (clubb_config_flags.Lscale_method == Lscale_tau) {

    // Note: someday when we can successfully run with all these parameters
    // having a value of 1, the "Warning" messages should be removed and the
    // "Fatal error" messages should be uncommented.

    // C1 must have a value of 1
    int any_parameter_relation_2 = 0;
    Kokkos::parallel_reduce(
        "check_clubb_settings_api::relation2", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& hits) {
          if (params((1 - 1 + column * (1)), iC1) > one ||
              params((1 - 1 + column * (1)), iC1) < one)
            ++hits;
        },
        any_parameter_relation_2);
    if ((any_parameter_relation_2 != 0)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "When Lscale_method=2, C1 must have a value of 1." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C1 = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC1)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Warning in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // write(fstderr,*) "Fatal error in check_clubb_settings_api"
      //  General error -> set all entries to clubb_fatal_error
      // err_info%err_code = clubb_fatal_error
      //  C1 check
    }

    // C1b must have a value of 1
    int any_parameter_relation_3 = 0;
    Kokkos::parallel_reduce(
        "check_clubb_settings_api::relation3", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& hits) {
          if (params((1 - 1 + column * (1)), iC1b) > one ||
              params((1 - 1 + column * (1)), iC1b) < one)
            ++hits;
        },
        any_parameter_relation_3);
    if ((any_parameter_relation_3 != 0)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "When Lscale_method=2, C1b must have a value of 1." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C1b = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC1b)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Warning in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // write(fstderr,*) "Fatal error in check_clubb_settings_api"
      //  General error -> set all entries to clubb_fatal_error
      // err_info%err_code = clubb_fatal_error
      //  C1b check
    }

    // C2rt must have a value of 1
    int any_parameter_relation_4 = 0;
    Kokkos::parallel_reduce(
        "check_clubb_settings_api::relation4", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& hits) {
          if (params((1 - 1 + column * (1)), iC2rt) > one ||
              params((1 - 1 + column * (1)), iC2rt) < one)
            ++hits;
        },
        any_parameter_relation_4);
    if ((any_parameter_relation_4 != 0)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "When Lscale_method=2, C2rt must have a value of 1." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C2rt = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC2rt)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Warning in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // write(fstderr,*) "Fatal error in check_clubb_settings_api"
      //  General error -> set all entries to clubb_fatal_error
      // err_info%err_code = clubb_fatal_error
      //  C2rt check
    }

    // C2thl must have a value of 1
    int any_parameter_relation_5 = 0;
    Kokkos::parallel_reduce(
        "check_clubb_settings_api::relation5", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& hits) {
          if (params((1 - 1 + column * (1)), iC2thl) > one ||
              params((1 - 1 + column * (1)), iC2thl) < one)
            ++hits;
        },
        any_parameter_relation_5);
    if ((any_parameter_relation_5 != 0)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "When Lscale_method=2, C2thl must have a value of 1." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C2thl = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC2thl)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Warning in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // write(fstderr,*) "Fatal error in check_clubb_settings_api"
      //  General error -> set all entries to clubb_fatal_error
      // err_info%err_code = clubb_fatal_error
      //  C2thl check
    }

    // C2rtthl must have a value of 1
    int any_parameter_relation_6 = 0;
    Kokkos::parallel_reduce(
        "check_clubb_settings_api::relation6", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& hits) {
          if (params((1 - 1 + column * (1)), iC2rtthl) > one ||
              params((1 - 1 + column * (1)), iC2rtthl) < one)
            ++hits;
        },
        any_parameter_relation_6);
    if ((any_parameter_relation_6 != 0)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "When Lscale_method=2, C2rtthl must have a value of 1." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C2rtthl = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC2rtthl)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Warning in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // write(fstderr,*) "Fatal error in check_clubb_settings_api"
      //  General error -> set all entries to clubb_fatal_error
      // err_info%err_code = clubb_fatal_error
      //  C2rtthl check
    }

    // C6rt must have a value of 1
    int any_parameter_relation_7 = 0;
    Kokkos::parallel_reduce(
        "check_clubb_settings_api::relation7", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& hits) {
          if (params((1 - 1 + column * (1)), iC6rt) > one ||
              params((1 - 1 + column * (1)), iC6rt) < one)
            ++hits;
        },
        any_parameter_relation_7);
    if ((any_parameter_relation_7 != 0)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "When Lscale_method=2, C6rt must have a value of 1." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C6rt = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC6rt)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Warning in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // write(fstderr,*) "Fatal error in check_clubb_settings_api"
      //  General error -> set all entries to clubb_fatal_error
      // err_info%err_code = clubb_fatal_error
      //  C6rt check
    }

    // C6rtb must have a value of 1
    int any_parameter_relation_8 = 0;
    Kokkos::parallel_reduce(
        "check_clubb_settings_api::relation8", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& hits) {
          if (params((1 - 1 + column * (1)), iC6rtb) > one ||
              params((1 - 1 + column * (1)), iC6rtb) < one)
            ++hits;
        },
        any_parameter_relation_8);
    if ((any_parameter_relation_8 != 0)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "When Lscale_method=2, C6rtb must have a value of 1." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C6rtb = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC6rtb)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Warning in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // write(fstderr,*) "Fatal error in check_clubb_settings_api"
      //  General error -> set all entries to clubb_fatal_error
      // err_info%err_code = clubb_fatal_error
      //  C6rtb check
    }

    // C6thl must have a value of 1
    int any_parameter_relation_9 = 0;
    Kokkos::parallel_reduce(
        "check_clubb_settings_api::relation9", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& hits) {
          if (params((1 - 1 + column * (1)), iC6thl) > one ||
              params((1 - 1 + column * (1)), iC6thl) < one)
            ++hits;
        },
        any_parameter_relation_9);
    if ((any_parameter_relation_9 != 0)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "When Lscale_method=2, C6thl must have a value of 1." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C6thl = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC6thl)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Warning in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // write(fstderr,*) "Fatal error in check_clubb_settings_api"
      //  General error -> set all entries to clubb_fatal_error
      // err_info%err_code = clubb_fatal_error
      //  C6thl check
    }

    // C6thlb must have a value of 1
    int any_parameter_relation_10 = 0;
    Kokkos::parallel_reduce(
        "check_clubb_settings_api::relation10", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& hits) {
          if (params((1 - 1 + column * (1)), iC6thlb) > one ||
              params((1 - 1 + column * (1)), iC6thlb) < one)
            ++hits;
        },
        any_parameter_relation_10);
    if ((any_parameter_relation_10 != 0)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "When Lscale_method=2, C6thlb must have a value of 1." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C6thlb = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC6thlb)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Warning in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // write(fstderr,*) "Fatal error in check_clubb_settings_api"
      //  General error -> set all entries to clubb_fatal_error
      // err_info%err_code = clubb_fatal_error
      //  C6thlb check
    }

    // C14 must have a value of 1
    int any_parameter_relation_11 = 0;
    Kokkos::parallel_reduce(
        "check_clubb_settings_api::relation11", ngrdcol,
        KOKKOS_LAMBDA(const int column, int& hits) {
          if (params((1 - 1 + column * (1)), iC14) > one ||
              params((1 - 1 + column * (1)), iC14) < one)
            ++hits;
        },
        any_parameter_relation_11);
    if ((any_parameter_relation_11 != 0)) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "When Lscale_method=2, C14 must have a value of 1." << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "C14 = " << " "
                           << dump_values(Kokkos::subview(params, Kokkos::ALL, iC14)) << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Warning in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // write(fstderr,*) "Fatal error in check_clubb_settings_api"
      //  General error -> set all entries to clubb_fatal_error
      // err_info%err_code = clubb_fatal_error
      //  C14 check
    }

    // tau-based mixing length
  }

  if (l_implemented) {

    if (clubb_config_flags.l_rtm_nudge) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "l_rtm_nudge must be set to .false. when l_implemented = .true."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }

    if (clubb_config_flags.l_uv_nudge) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "l_rtm_nudge must be set to .false. when l_implemented = .true."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }

    if (environment.thlm_sponge_damp_settings__l_sponge_damping) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "thlm_sponge_damp_settings%l_sponge_damping must be set to .false. "
                              "when  l_implemented = .true."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }

    if (environment.rtm_sponge_damp_settings__l_sponge_damping) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "rtm_sponge_damp_settings%l_sponge_damping must be set to .false. "
                              "when  l_implemented = .true."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }

    if (environment.uv_sponge_damp_settings__l_sponge_damping) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "uv_sponge_damp_settings%l_sponge_damping must be set to .false. "
                              "when  l_implemented = .true."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }

    if (environment.wp2_sponge_damp_settings__l_sponge_damping) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "wp2_sponge_damp_settings%l_sponge_damping must be set to .false. "
                              "when  l_implemented = .true."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }

    if (environment.wp3_sponge_damp_settings__l_sponge_damping) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "wp3_sponge_damp_settings%l_sponge_damping must be set to .false. "
                              "when  l_implemented = .true."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }

    if (environment.up2_vp2_sponge_damp_settings__l_sponge_damping) {
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << err_info.err_header_global << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "up2_vp2_sponge_damp_settings%l_sponge_damping must be set to "
                              ".false. when  l_implemented = .true."
                           << "\n";
        std::cerr << diagnostic_message.str();
      }
      {
        std::ostringstream diagnostic_message;
        diagnostic_message << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic_message << "Fatal error in check_clubb_settings_api" << "\n";
        std::cerr << diagnostic_message.str();
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    }
  }

  return;
}

} // namespace clubb_kokkos
