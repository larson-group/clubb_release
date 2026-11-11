// Core diagnostic sequence; see numerical_check.hpp for its public contract.
#include "numerical_check.hpp"
#include "error_code.hpp"
#include "f2k/procedure_bridge.hpp"

#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>

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

// Other native procedures preceding this scoped port remain native until their
// live callers are converted. Definitions below retain the source file's order.
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

// -----------------------------------------------------------------------
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

// -----------------------------------------------------------------------
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

// -----------------------------------------------------------------------
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

} // namespace clubb_kokkos
