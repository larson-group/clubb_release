// Mixing-length port; the Native source is src/CLUBB_core/mixing_length.F90.
// Public APIs/field units are documented in mixing_length.hpp. Scratch is per call.
// Internal column/level indices and grid endpoints are zero-based;
// parameter indices already use the public C++ zero-based constants.
#include "scratch_workspace.hpp"
#include "code_timer_module.hpp"
#include "mixing_length.hpp"
#include "constants_clubb.hpp"
#include "parameter_indices.hpp"
#include "model_flags.hpp"
#include "error_code.hpp"
#include "stats_netcdf.hpp"
#include "advance_helper_module.hpp"
#include "numerical_check.hpp"
#include "saturation.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace clubb_kokkos {

namespace {
void compute_mixing_length(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol,
                           const grid& gr, Real2D thvm,
                           Real2D thlm,                                                        // In
                           Real2D rtm, Real2D em, Real1D Lscale_max, Real2D p_in_Pa,           // In
                           Real2D exner, Real2D thv_ds, Real1D Lscale_mu_entrain, Real1D lmin, // In
                           int saturation_formula,                                             // In
                           bool l_implemented,                                                 // In
                           err_info_type& err_info,                            // InOut
                           Real2D Lscale, Real2D Lscale_up, Real2D Lscale_down // Out
);
void compute_coherent_mixing_length(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol,
                                    const grid& gr, Real2D thvm,
                                    Real2D thlm, // In
                                    Real2D rtm, const pdf_parameter& pdf_params, Real2D em,
                                    Real1D Lscale_max,                                     // In
                                    Real2D p_in_Pa,                                        // In
                                    Real2D exner, Real2D thv_ds, Real1D Lscale_mu_entrain, // In
                                    Real1D Lscale_tke_multiplier,                          // In
                                    Real1D Lscale_pdf_launch_frac,                         // In
                                    Real1D Lscale_up_weight, Real1D lmin,                  // In
                                    bool l_implemented,                                    // In
                                    err_info_type& err_info,                               // InOut
                                    Real2D Lscale, Real2D Lscale_up, Real2D Lscale_down    // Out
);
KOKKOS_INLINE_FUNCTION real_t sat_mixrat_liq_flatau(real_t pressure, real_t temperature // In
);
} // namespace

void calc_Lscale(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, const grid& gr,
                 bool l_implemented,                                                        // In
                 Real1D host_dx, Real1D host_dy,                                            // In
                 Real2D p_in_Pa, Real2D exner, Real2D rtm, Real2D thlm, Real2D thvm,        // In
                 Real2D thlp2, Real2D rtp2, Real2D rtpthlp,                                 // In
                 const pdf_parameter& pdf_params, Real2D em, Real2D thv_ds_zt, Real1D lmin, // In
                 Real1D upwp_sfc, Real1D vpwp_sfc, Real2D ddzt_umvm_sqd,                    // In
                 Real2D ice_supersat_frac,                                                  // In
                 real_t ufmin, real_t tau_const, Real1D sfc_elevation, Real2D clubb_params, // In
                 int saturation_formula, bool l_Lscale_plume_centered,                      // In
                 int Lscale_method,                                                         // In
                 bool l_e3sm_config,                                                        // In
                 bool l_smooth_Heaviside_tau_wpxp,                                          // In
                 bool l_modify_limiters_for_cnvg_test,                                      // In
                 bool l_use_invrs_tau_N2_iso,                                               // In
                 Real2D brunt_vaisala_freq_sqd_smth,                                        // In
                 stats_type& stats, err_info_type& err_info,                                // InOut
                 Real2D invrs_tau_zt, Real2D invrs_tau_zm, Real2D invrs_tau_xp2_zm,         // Out
                 Real2D invrs_tau_wp3_zt,                                                   // Out
                 Real2D invrs_tau_C1_zm, Real2D invrs_tau_C4_zm,                            // Out
                 Real2D invrs_tau_C6_zm, Real2D invrs_tau_C14_zm,                           // Out
                 Real2D tau_max_zm, Real2D tau_max_zt, Real2D tau_zm,                       // Out
                 Real2D Lscale, Real2D Lscale_zm, Real2D Lscale_up, Real2D Lscale_down      // Out
) {
  // Description:
  // Compute CLUBB's mixing length and dissipation time scales, including
  // Richardson-number-dependent damping terms and Lscale-related stats.
  // Read thermodynamic fields on nzt levels, moments on nzm levels, and
  // per-column tunable parameters; Lscale_method selects methods 1--3.
  // Return lengths [m], time scales [s], and inverse times [1/s], and update
  // stats and err_info. Grid and PDF state must be initialized by the caller.
  // References:
  // Guo et al. (2021, JAMES) for the diagnosed Lscale-from-tau option.
  // -----------------------------------------------------------------------

  // Arguments (Native order)
  // nzm: Number of momentum levels
  // nzt: Number of thermodynamic levels
  // ngrdcol: Number of grid columns
  // gr: Grid structure
  // l_implemented: True when current configuration is implemented
  // host_dx: Host-model grid spacing in x direction [m]
  // host_dy: Host-model grid spacing in y direction [m]
  // p_in_Pa: Pressure [Pa]
  // exner: Exner function [-]
  // rtm: Total water mixing ratio [kg/kg]
  // thlm: Liquid water potential temperature [K]
  // thvm: Virtual potential temperature [K]
  // thlp2: Variance of thl [K^2]
  // rtp2: Variance of rt [(kg/kg)^2]
  // rtpthlp: Covariance of rt and thl [K kg/kg]
  // pdf_params: PDF parameters
  // em: Turbulent kinetic energy [m^2/s^2]
  // thv_ds_zt: Reference-state virtual potential temperature [K]
  // lmin: Minimum mixing length [m]
  // upwp_sfc: Surface value of u'w' [m^2/s^2]
  // vpwp_sfc: Surface value of v'w' [m^2/s^2]
  // ddzt_umvm_sqd: Squared vertical shear of horizontal wind [s^-2]
  // ice_supersat_frac: Ice supersaturation fraction [-]
  // ufmin: Minimum velocity scale [m/s]
  // tau_const: Constant eddy dissipation time scale [s]
  // sfc_elevation: Surface elevation [m]
  // clubb_params: Tunable CLUBB parameters
  // saturation_formula: Saturation vapor pressure formula selector
  // l_Lscale_plume_centered: Whether to use plume-centered Lscale calculation
  // Lscale_method: Mixing-length formulation
  // l_e3sm_config: Whether E3SM configuration options are active
  // l_smooth_Heaviside_tau_wpxp: Whether to smooth Heaviside functions in wpxp tau
  // l_modify_limiters_for_cnvg_test: Whether to modify limiters for convergence tests
  // l_use_invrs_tau_N2_iso: Whether to use isotropic N2 inverse tau
  // brunt_vaisala_freq_sqd_smth: Smoothed Brunt-Vaisala frequency squared [s^-2]
  // stats: Statistics accumulator
  // err_info: Error information
  // invrs_tau_zt: Inverse time scale on thermodynamic levels [1/s]
  // invrs_tau_zm: Inverse time scale on momentum levels [1/s]
  // invrs_tau_xp2_zm: Inverse time scale for xp2 terms on momentum levels [1/s]
  // invrs_tau_wp3_zt: Inverse time scale for wp3 terms on thermodynamic levels [1/s]
  // invrs_tau_C1_zm: Inverse time scale used for C1 term on momentum levels [1/s]
  // invrs_tau_C4_zm: Inverse time scale used for C4 terms on momentum levels [1/s]
  // invrs_tau_C6_zm: Inverse time scale used for C6 term on momentum levels [1/s]
  // invrs_tau_C14_zm: Inverse time scale used for C14 terms on momentum levels [1/s]
  // tau_max_zm: Maximum allowable eddy dissipation time scale on momentum levels [s]
  // tau_max_zt: Maximum allowable eddy dissipation time scale on thermodynamic levels [s]
  // tau_zm: Eddy dissipation time scale on momentum levels [s]
  // Lscale: Turbulent mixing length on thermodynamic levels [m]
  // Lscale_zm: Turbulent mixing length on momentum levels [m]
  // Lscale_up: Upward component of turbulent mixing length [m]
  // Lscale_down: Downward component of turbulent mixing length [m]

  auto workspace_scope = workspace.scope();
  Real1D Lscale_max =
      workspace_scope.view<Real1D>(Kokkos::WithoutInitializing, "calc_Lscale::Lscale_max",
                                   ngrdcol); // Maximum allowable mixing length [m]
  constexpr bool l_smooth_min_max = false;
  Real2D em_zt      = (Lscale_method != Lscale_tau)
                          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                         "calc_Lscale::em_zt", ngrdcol, nzt)
                          : Real2D{}; // TKE interpolated to thermodynamic levels [m^2/s^2]
  Real2D thlp2_zt   = (Lscale_method != Lscale_tau)
                          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                         "calc_Lscale::thlp2_zt", ngrdcol, nzt)
                          : Real2D{}; // Variance of thl on thermodynamic levels [K^2]
  Real2D rtp2_zt    = (Lscale_method != Lscale_tau)
                          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                                         "calc_Lscale::rtp2_zt", ngrdcol, nzt)
                          : Real2D{}; // Variance of rt on thermodynamic levels [(kg/kg)^2]
  Real2D rtpthlp_zt = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "calc_Lscale::rtpthlp_zt", ngrdcol,
      nzt); // Covariance of rt and thl on thermodynamic levels [K kg/kg]
  Real2D tau_zt =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "calc_Lscale::tau_zt", ngrdcol,
                                   nzt); // Eddy dissipation time scale on thermodynamic levels [s]
  Real2D invrs_tau_sfc =
      (Lscale_method == Lscale_tau)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "calc_Lscale::invrs_tau_sfc",
                                         ngrdcol, nzm)
          : Real2D{}; // Surface contribution to inverse tau [1/s]
  Real2D invrs_tau_no_N2_zm =
      (Lscale_method == Lscale_tau)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "calc_Lscale::invrs_tau_no_N2_zm", ngrdcol, nzm)
          : Real2D{}; // Inverse tau without N2 contribution on momentum levels [1/s]
  Real2D invrs_tau_bkgnd =
      (Lscale_method == Lscale_tau)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "calc_Lscale::invrs_tau_bkgnd", ngrdcol, nzm)
          : Real2D{}; // Background inverse tau [1/s]
  Real2D invrs_tau_N2_iso =
      (Lscale_method == Lscale_tau)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "calc_Lscale::invrs_tau_N2_iso", ngrdcol, nzm)
          : Real2D{}; // Isotropic N2 contribution to inverse tau [1/s]
  Real2D invrs_tau_wp2_zm = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "calc_Lscale::invrs_tau_wp2_zm", ngrdcol,
      nzm); // Inverse time scale for wp2 terms on momentum levels [1/s]
  Real2D invrs_tau_wp3_zm = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "calc_Lscale::invrs_tau_wp3_zm", ngrdcol,
      nzm); // Inverse time scale for wp3 terms on momentum levels [1/s]
  Real2D invrs_tau_wpxp_zm = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "calc_Lscale::invrs_tau_wpxp_zm", ngrdcol,
      nzm); // Inverse time scale for wpxp terms on momentum levels [1/s]
  Real2D invrs_tau_shear =
      (Lscale_method == Lscale_tau)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "calc_Lscale::invrs_tau_shear", ngrdcol, nzm)
          : Real2D{}; // Shear contribution to inverse tau [1/s]
  Real2D Ri_zm =
      workspace_scope.view<Real2D>(Kokkos::WithoutInitializing, "calc_Lscale::Ri_zm", ngrdcol,
                                   nzm); // Gradient Richardson number on momentum levels [-]
  Real2D ddzt_umvm_sqd_clipped =
      (l_smooth_min_max)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "calc_Lscale::ddzt_umvm_sqd_clipped", ngrdcol, nzm)
          : Real2D{}; // Clipped squared vertical shear of horizontal wind [s^-2]
  Real2D brunt_vaisala_freq_clipped =
      (l_smooth_min_max)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "calc_Lscale::brunt_vaisala_freq_clipped", ngrdcol, nzm)
          : Real2D{}; // Clipped Brunt-Vaisala frequency squared [s^-2]

  timer_start("calc_Lscale" // In
  );
  // Determine the maximum allowable value for Lscale (in meters).
  set_Lscale_max(ngrdcol, l_implemented, host_dx, host_dy, // In
                 Lscale_max                                // Out
  );
  // Interpolate rtpthlp here so the diagnostic stats value is generated at
  // the same point as the Lscale calculation that consumes it.
  zm2zt_api(nzm, nzt, ngrdcol, gr, rtpthlp, rtpthlp_zt);
  if (stats.l_sample) {
    stats_update("rtpthlp_zt", rtpthlp_zt, // In
                 stats                     // InOut
    );
  }
  // Calculate Richardson number Ri_zm.
  if (l_modify_limiters_for_cnvg_test) {
    calc_Ri_zm(nzm, ngrdcol, brunt_vaisala_freq_sqd_smth, ddzt_umvm_sqd, // In
               real_t(0.0), real_t(1.0e-12),                             // In
               Ri_zm                                                     // Out
    );
    zm2zt2zm(workspace, nzm, nzt, ngrdcol, gr, Ri_zm, Ri_zm);
  } else {
    if (l_smooth_min_max) {
      smooth_max(nzm, ngrdcol, real_t(1.0e-7), brunt_vaisala_freq_sqd_smth,
                 (real_t(1.0e-4) * min_max_smth_mag), brunt_vaisala_freq_clipped);
      smooth_max(nzm, ngrdcol, ddzt_umvm_sqd, real_t(1.0e-7), (real_t(1.0e-6) * min_max_smth_mag),
                 ddzt_umvm_sqd_clipped);
      calc_Ri_zm(nzm, ngrdcol, brunt_vaisala_freq_clipped,        // In
                 ddzt_umvm_sqd_clipped, real_t(0.0), real_t(0.0), // In
                 Ri_zm                                            // Out
      );
    } else {
      calc_Ri_zm(nzm, ngrdcol, brunt_vaisala_freq_sqd_smth, ddzt_umvm_sqd, // In
                 real_t(1.0e-7), real_t(1.0e-7),                           // In
                 Ri_zm                                                     // Out
      );
    }
  }
  if (stats.l_sample) {
    stats_update("Ri_zm", Ri_zm, // In
                 stats           // InOut
    );
  }
  if ((Lscale_method != Lscale_tau)) {
    // Compute Lscale first using the buoyant parcel calculation.
    // Interpolate thlp2 and rtp2 to thermodynamic levels.
    zm2zt_api(nzm, nzt, ngrdcol, gr, thlp2, thlp2_zt, (thl_tol * thl_tol));
    // Positive def. quantity
    zm2zt_api(nzm, nzt, ngrdcol, gr, rtp2, rtp2_zt, (rt_tol * rt_tol));
    // Positive def. quantity
    calc_Lscale_directly(workspace, ngrdcol, nzm, nzt, gr,          // In
                         l_implemented, p_in_Pa,                    // In
                         exner, rtm, thlm, thvm,                    // In
                         rtp2_zt, thlp2_zt, rtpthlp_zt, pdf_params, // In
                         em, thv_ds_zt, Lscale_max, lmin,           // In
                         clubb_params,                              // In
                         saturation_formula,                        // In
                         l_Lscale_plume_centered,                   // In
                         Lscale_method,                             // In
                         stats, err_info,                           // InOut
                         Lscale, Lscale_up, Lscale_down             // Out
    );
    if (clubb_at_least_debug_level_api(0)) {
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
          diagnostic << "Error calling calc_Lscale_directly in calc_Lscale" << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        return;
      }
    }
    // Calculate CLUBB's turbulent eddy-turnover time scale as CLUBB's
    // length scale divided by a velocity scale.
    zm2zt_api(nzm, nzt, ngrdcol, gr, em, em_zt, em_min);
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "calc_Lscale::cells_0", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i  = (cell % ngrdcol);
            const int k  = (cell / ngrdcol);
            tau_zt(i, k) = ((Lscale(i, k) / Kokkos::sqrt(em_zt(i, k))) < clubb_params(i, itaumax)
                                ? (Lscale(i, k) / Kokkos::sqrt(em_zt(i, k)))
                                : clubb_params(i, itaumax));
          });
    }
    zt2zm_api(nzm, nzt, ngrdcol, gr, Lscale, Lscale_zm, zero_threshold);
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "calc_Lscale::cells_1", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            tau_zm(i, k) =
                ((Lscale_zm(i, k) / Kokkos::sqrt((em_min > em(i, k) ? em_min : em(i, k)))) <
                         clubb_params(i, itaumax)
                     ? (Lscale_zm(i, k) / Kokkos::sqrt((em_min > em(i, k) ? em_min : em(i, k))))
                     : clubb_params(i, itaumax));
          });
    }
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "calc_Lscale::cells_2", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i             = (cell % ngrdcol);
            const int k             = (cell / ngrdcol);
            invrs_tau_zm(i, k)      = (one / tau_zm(i, k));
            invrs_tau_wp2_zm(i, k)  = invrs_tau_zm(i, k);
            invrs_tau_xp2_zm(i, k)  = invrs_tau_zm(i, k);
            invrs_tau_wpxp_zm(i, k) = invrs_tau_zm(i, k);
            invrs_tau_wp3_zm(i, k)  = invrs_tau_zm(i, k);
            tau_max_zm(i, k)        = clubb_params(i, itaumax);
          });
    }
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "calc_Lscale::cells_3", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i            = (cell % ngrdcol);
            const int k            = (cell / ngrdcol);
            invrs_tau_zt(i, k)     = (one / tau_zt(i, k));
            invrs_tau_wp3_zt(i, k) = invrs_tau_zt(i, k);
            tau_max_zt(i, k)       = clubb_params(i, itaumax);
          });
    }
  } else {
    // Diagnose simple tau and Lscale.
    diagnose_Lscale_from_tau(workspace, nzm, nzt, ngrdcol, gr,                      // In
                             upwp_sfc, vpwp_sfc, ddzt_umvm_sqd,                     // In
                             ice_supersat_frac,                                     // In
                             em,                                                    // In
                             ufmin, tau_const,                                      // In
                             sfc_elevation, Lscale_max,                             // In
                             clubb_params,                                          // In
                             stats,                                                 // InOut
                             l_e3sm_config,                                         // In
                             l_smooth_Heaviside_tau_wpxp,                           // In
                             brunt_vaisala_freq_sqd_smth, Ri_zm,                    // In
                             err_info,                                              // InOut
                             invrs_tau_zt, invrs_tau_zm,                            // Out
                             invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd,    // Out
                             invrs_tau_shear, invrs_tau_N2_iso,                     // Out
                             invrs_tau_wp2_zm, invrs_tau_xp2_zm,                    // Out
                             invrs_tau_wp3_zm, invrs_tau_wp3_zt, invrs_tau_wpxp_zm, // Out
                             tau_max_zm, tau_max_zt, tau_zm, tau_zt,                // Out
                             Lscale, Lscale_up, Lscale_down                         // Out
    );
    if (clubb_at_least_debug_level_api(0)) {
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
          diagnostic << "Error calling diagnose_Lscale_from_tau in calc_Lscale" << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        return;
      }
    }
  }
  zt2zm_api(nzm, nzt, ngrdcol, gr, Lscale, Lscale_zm);
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "calc_Lscale::cells_4", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i            = (cell % ngrdcol);
          const int k            = (cell / ngrdcol);
          invrs_tau_C6_zm(i, k)  = invrs_tau_wpxp_zm(i, k);
          invrs_tau_C1_zm(i, k)  = invrs_tau_wp2_zm(i, k);
          invrs_tau_C14_zm(i, k) = invrs_tau_wp2_zm(i, k);
        });
  }

  if (((Lscale_method != Lscale_tau) && l_use_invrs_tau_N2_iso)) {
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "Error! l_use_invrs_tau_N2_iso is not used when Lscale_method is not the tau "
                    "method.If you want to use Lscale code, go to file "
                    "src/CLUBB_core/mixing_length.F90 and change l_use_invrs_tau_N2_iso to false"
                 << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    throw std::runtime_error("Native mixing_length ERROR STOP");
  }
  if (!l_use_invrs_tau_N2_iso) {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "calc_Lscale::cells_6", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i           = (cell % ngrdcol);
            const int k           = (cell / ngrdcol);
            invrs_tau_C4_zm(i, k) = invrs_tau_wp2_zm(i, k);
          });
    }
  } else {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "calc_Lscale::cells_7", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i           = (cell % ngrdcol);
            const int k           = (cell / ngrdcol);
            invrs_tau_C4_zm(i, k) = invrs_tau_N2_iso(i, k);
          });
    }
  }
  if (stats.l_sample) {
    stats_update("Lscale", Lscale, // In
                 stats             // InOut
    );
    stats_update("Lscale_up", Lscale_up, // In
                 stats                   // InOut
    );
    stats_update("Lscale_down", Lscale_down, // In
                 stats                       // InOut
    );
    stats_update("tau_zm", tau_zm, // In
                 stats             // InOut
    );
    stats_update("invrs_tau_zm", invrs_tau_zm, // In
                 stats                         // InOut
    );
    stats_update("invrs_tau_xp2_zm", invrs_tau_xp2_zm, // In
                 stats                                 // InOut
    );
    stats_update("invrs_tau_wp2_zm", invrs_tau_wp2_zm, // In
                 stats                                 // InOut
    );
    stats_update("invrs_tau_wpxp_zm", invrs_tau_wpxp_zm, // In
                 stats                                   // InOut
    );
    stats_update("invrs_tau_wp3_zm", invrs_tau_wp3_zm, // In
                 stats                                 // InOut
    );
    stats_update("tau_zt", tau_zt, // In
                 stats             // InOut
    );
    if ((Lscale_method == Lscale_tau)) {
      stats_update("invrs_tau_no_N2_zm", invrs_tau_no_N2_zm, // In
                   stats                                     // InOut
      );
      stats_update("invrs_tau_bkgnd", invrs_tau_bkgnd, // In
                   stats                               // InOut
      );
      stats_update("invrs_tau_sfc", invrs_tau_sfc, // In
                   stats                           // InOut
      );
      stats_update("invrs_tau_shear", invrs_tau_shear, // In
                   stats                               // InOut
      );
    }
  }
  timer_stop("calc_Lscale" // In
  );
  return;
}

void set_Lscale_max(int ngrdcol, bool l_implemented, Real1D host_dx, Real1D host_dy, // In
                    Real1D Lscale_max                                                // Out
) {
  // Description:
  // This subroutine sets the value of Lscale_max, which is the maximum
  // allowable value of Lscale.  For standard CLUBB, it is set to a very large
  // value so that Lscale will not be limited.  However, when CLUBB is running
  // as part of a host model, the value of Lscale_max is dependent on the size
  // of the host model's horizontal grid spacing.  The smaller the host model's
  // horizontal grid spacing, the smaller the value of Lscale_max.  When Lscale
  // is limited to a small value, the value of time-scale Tau is reduced, which
  // in turn produces greater damping on CLUBB's turbulent parameters.  This
  // is the desired effect on turbulent parameters for a host model with small
  // horizontal grid spacing, for small areas usually contain much less
  // variation in meteorological quantities than large areas.
  // References:
  // None
  // -----------------------------------------------------------------------

  // Arguments (Native order)
  // ngrdcol: Number of grid columns
  // l_implemented: True if CLUBB is being implemented and run in a host model
  // host_dx: Host model's east-west horizontal grid spacing     [m]
  // host_dy: Host model's north-south horizontal grid spacing   [m]
  // Lscale_max: Maximum allowable value for Lscale   [m]

  // Determine the maximum allowable value for Lscale (in meters).
  if (l_implemented) {
    Kokkos::parallel_for(
        "set_Lscale_max::columns_0", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i   = column;
          Lscale_max(i) = (real_t(0.25) * (host_dx(i) < host_dy(i) ? host_dx(i) : host_dy(i)));
        });
  } else {
    Kokkos::parallel_for(
        "set_Lscale_max::columns_1", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i   = column;
          Lscale_max(i) = real_t(1.0e5);
        });
  }
  return;
}

namespace {

void compute_mixing_length(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol,
                           const grid& gr, Real2D thvm,
                           Real2D thlm,                                                        // In
                           Real2D rtm, Real2D em, Real1D Lscale_max, Real2D p_in_Pa,           // In
                           Real2D exner, Real2D thv_ds, Real1D Lscale_mu_entrain, Real1D lmin, // In
                           int saturation_formula,                                             // In
                           bool l_implemented,                                                 // In
                           err_info_type& err_info,                            // InOut
                           Real2D Lscale, Real2D Lscale_up, Real2D Lscale_down // Out
) {
  // Description:
  // Larson's 5th moist, nonlocal length scale
  // References:
  // Section 3b ( /Eddy length formulation/ ) of
  // ``A PDF-Based Model for Boundary Layer Clouds. Part I:
  // Method and Model Description'' Golaz, et al. (2002)
  // JAS, Vol. 59, pp. 3540--3551.
  // Notes:
  // The equation for the rate of change of theta_l and r_t of the parcel with
  // respect to height, due to entrainment, is:
  // d(thl_par)/dz = - Lscale_mu_entrain * ( thl_parcel - thl_environment );
  // d(rt_par)/dz = - Lscale_mu_entrain * ( rt_parcel - rt_environment );
  // where Lscale_mu_entrain is the entrainment rate,
  // such that:
  // Lscale_mu_entrain = (1/m)*(dm/dz);
  // where m is the mass of the parcel.  The value of Lscale_mu_entrain is set to be a
  // constant.
  // The differential equations are solved for given the boundary condition
  // and given the fact that the value of thl_environment and rt_environment
  // are treated as changing linearly for a parcel of air from one grid level
  // to the next.
  // For the special case where entrainment rate, Lscale_mu_entrain, is set to 0,
  // thl_parcel and rt_parcel remain constant
  // The equation for Lscale_up is:
  // INT(z_i:z_i+Lscale_up) g * ( thv_par - thvm ) / thvm dz = -em(z_i);
  // and for Lscale_down
  // INT(z_i-Lscale_down:z_i) g * ( thv_par - thvm ) / thvm dz = em(z_i);
  // where thv_par is theta_v of the parcel, thvm is the mean
  // environmental value of theta_v, z_i is the altitude that the parcel
  // started from, and em is the mean value of TKE at
  // altitude z_i (which gives the parcel its initial boost)
  // The increment of CAPE (convective air potential energy) for any two
  // successive vertical levels is:
  // Upwards:
  // CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
  // Downwards:
  // CAPE_incr = INT(z_(-1):z_0) g * ( thv_par - thvm ) / thvm dz
  // Thus, the derivative of CAPE with respect to height is:
  // dCAPE/dz = g * ( thv_par - thvm ) / thvm.
  // A purely trapezoidal rule is used between levels, and is considered
  // to vary linearly at all altitudes.  Thus, dCAPE/dz is considered to be
  // of the form:  A * (z-zo) + dCAPE/dz|_(z_0),
  // where A = ( dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) ) / ( z_1 - z_0 )
  // The integral is evaluated to find the CAPE increment between two
  // successive vertical levels.  The result either adds to or depletes
  // from the total amount of energy that keeps the parcel ascending/descending.
  // IMPORTANT NOTE:
  // This subroutine has been optimized by adding precalculations, rearranging
  // equations to avoid divides, and modifying the algorithm entirely.
  // -Gunther Huebler
  // The algorithm previously used looped over every grid level, following a
  // a parcel up from its initial grid level to its max. The very nature of this
  // algorithm is an N^2
  // --------------------------------------------------------------------------------
  // Lscale_mu_entrain = (1/M) dM/dz > 0.  Lscale_mu_entrain=0 for no entrainment.
  // Siebesma recommends Lscale_mu_entrain=2e-3, although most schemes use Lscale_mu_entrain=1e-4
  // When Lscale_mu_entrain was fixed, we used the value Lscale_mu_entrain = 6.e-4

  // Arguments (Native order)
  // nzm: Number of momentum levels
  // nzt: Number of thermodynamic levels
  // ngrdcol: Number of grid columns
  // gr: Grid heights, interpolation weights, and orientation
  // thvm: Virtual potential temp. on themodynamic level  [K]
  // thlm: Liquid potential temp. on themodynamic level   [K]
  // rtm: Total water mixing ratio on themodynamic level [kg/kg]
  // em: em = 3/2 * w'^2; on momentum level             [m^2/s^2]
  // Lscale_max: Maximum allowable value for Lscale             [m]
  // p_in_Pa: Pressure on thermodynamic level                [Pa]
  // exner: Exner function on thermodynamic level          [-]
  // thv_ds: Dry, base-state theta_v on thermodynamic level [K]
  // Lscale_mu_entrain: Fractional entrainment rate per unit altitude [1/m]
  // lmin: CLUBB tunable parameter lmin
  // saturation_formula: Integer that stores the saturation formula to be used
  // l_implemented: True if CLUBB is being implemented and run in a host model
  // err_info: err_info struct containing err_code and err_header
  // Lscale: Mixing length      [m]
  // Lscale_up: Mixing length up   [m]
  // Lscale_down: Mixing length down [m]

  auto workspace_scope                 = workspace.scope();
  constexpr real_t zlmin               = real_t(0.1);  // Minimum value for Lscale [m]
  constexpr real_t Lscale_sfclyr_depth = real_t(500.); // [m]
  int i{-1};
  int start_index{-1};
  Real2D grav_on_thvm =
      workspace_scope.view<Real2D>("compute_mixing_length::grav_on_thvm", ngrdcol, nzt);
  Real2D Lv_coef = workspace_scope.view<Real2D>("compute_mixing_length::Lv_coef", ngrdcol, nzt);
  Real2D thl_par_j_precalc =
      workspace_scope.view<Real2D>("compute_mixing_length::thl_par_j_precalc", ngrdcol, nzt);
  Real2D rt_par_j_precalc =
      workspace_scope.view<Real2D>("compute_mixing_length::rt_par_j_precalc", ngrdcol, nzt);
  Real2D tl_par_1 = workspace_scope.view<Real2D>("compute_mixing_length::tl_par_1", ngrdcol, nzt);
  Real2D rt_par_1 = workspace_scope.view<Real2D>("compute_mixing_length::rt_par_1", ngrdcol, nzt);
  Real2D rsatl_par_1 =
      workspace_scope.view<Real2D>("compute_mixing_length::rsatl_par_1", ngrdcol, nzt);
  Real2D thl_par_1 = workspace_scope.view<Real2D>("compute_mixing_length::thl_par_1", ngrdcol, nzt);
  Real2D dCAPE_dz_1 =
      workspace_scope.view<Real2D>("compute_mixing_length::dCAPE_dz_1", ngrdcol, nzt);
  Real2D s_par_1  = workspace_scope.view<Real2D>("compute_mixing_length::s_par_1", ngrdcol, nzt);
  Real2D rc_par_1 = workspace_scope.view<Real2D>("compute_mixing_length::rc_par_1", ngrdcol, nzt);
  Real2D CAPE_incr_1 =
      workspace_scope.view<Real2D>("compute_mixing_length::CAPE_incr_1", ngrdcol, nzt);
  Real2D thv_par_1 = workspace_scope.view<Real2D>("compute_mixing_length::thv_par_1", ngrdcol, nzt);
  Real2D exp_mu_dzm =
      workspace_scope.view<Real2D>("compute_mixing_length::exp_mu_dzm", ngrdcol, nzm);
  Real2D invrs_dzm_on_mu =
      workspace_scope.view<Real2D>("compute_mixing_length::invrs_dzm_on_mu", ngrdcol, nzm);
  Real2D entrain_coef =
      workspace_scope.view<Real2D>("compute_mixing_length::entrain_coef", ngrdcol, nzm);
  Real2D tke_i = workspace_scope.view<Real2D>("compute_mixing_length::tke_i", ngrdcol, nzt);
  real_t Lv2_coef{};
  real_t invrs_Lscale_sfclyr_depth{};
  // Capture component Views and zero-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_grid_dir      = gr.grid_dir;
  const auto gr_dzm           = gr.dzm;
  const auto gr_invrs_dzm     = gr.invrs_dzm;
  const auto gr_k_lb_zt       = (gr.k_lb_zt);
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_ub_zt       = (gr.k_ub_zt);
  const auto gr_zt            = gr.zt;
  const auto gr_zm            = gr.zm;
  const auto gr_k_lb_zm       = (gr.k_lb_zm);

  // Native per-column fatal detection runs on device; print only failures on host.
  Kokkos::View<int*, Kokkos::LayoutLeft> bad_columns =
      workspace_scope.view<Kokkos::View<int*, Kokkos::LayoutLeft>>(
          "compute_mixing_length::bad_columns", ngrdcol);
  const auto err_code_check = err_info.err_code;
  int bad_count             = 0;
  Kokkos::parallel_reduce(
      "compute_mixing_length::input_check", ngrdcol,
      KOKKOS_LAMBDA(const int column, int& count) {
        bad_columns(column) = 0;
        if (Kokkos::abs(Lscale_mu_entrain(column)) < eps) {
          bad_columns(column)    = 1;
          err_code_check(column) = clubb_fatal_error;
          count += 1;
        }
      },
      bad_count);
  if (bad_count != 0) {
    const auto bad_host    = to_host_space(workspace, bad_columns);
    const auto header_host = to_host_space(workspace, err_info.err_header);
    const auto mu_host     = to_host_space(workspace, Lscale_mu_entrain);
    for (int column = 0; column < ngrdcol; ++column) {
      if (Kokkos::abs(mu_host(column)) < eps) {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic.write(header_host.data() + (column * header_host.extent(0)),
                         header_host.extent(0));
        diagnostic << '\n';
        diagnostic << "Lscale_mu_entrain = " << mu_host(column) << '\n';
        std::cerr << diagnostic.str();
      }
    }
  }
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Entrainment rate Lscale_mu_entrain cannot be 0" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Fatal error in subroutine compute_mixing_length" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      return;
    }
  }
  // Calculate initial turbulent kinetic energy for each grid level
  zm2zt_api(nzm, nzt, ngrdcol, gr, em, tke_i);
  // Initialize arrays and precalculate values for computational efficiency
  {
    const int level_count = nzt;
    Kokkos::parallel_for(
        "compute_mixing_length::cells_0", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          // Initialize up and down arrays
          Lscale_up(i, k)   = zlmin;
          Lscale_down(i, k) = zlmin;
          // Precalculate values to avoid unnecessary calculations later
          grav_on_thvm(i, k) = (grav / thvm(i, k));
          Lv_coef(i, k)      = ((Lv / (exner(i, k) * Cp)) - (ep2 * thv_ds(i, k)));
        });
  }
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "compute_mixing_length::cells_1", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          // Precalculate values to avoid unnecessary calculations later
          exp_mu_dzm(i, k) = Kokkos::exp(((-(Lscale_mu_entrain(i)) * gr_grid_dir) * gr_dzm(i, k)));
          invrs_dzm_on_mu(i, k) = ((gr_grid_dir * gr_invrs_dzm(i, k)) / Lscale_mu_entrain(i));
          entrain_coef(i, k)    = ((one - exp_mu_dzm(i, k)) * invrs_dzm_on_mu(i, k));
        });
  }
  // Precalculations of single values to avoid unnecessary calculations later
  Lv2_coef                  = ((ep * (Lv * Lv)) / (Rd * Cp));
  invrs_Lscale_sfclyr_depth = (one / Lscale_sfclyr_depth);
  // ---------------- Upwards Length Scale Calculation ----------------
  // Precalculate values for upward Lscale, these are useful only if a parcel can rise
  // more than one level. They are used in the equations that calculate thl and rt
  // recursively for a parcel as it ascends
  {
    const int level_count =
        (gr_grid_dir_indx > 0 ? gr_grid_dir_indx + gr_k_lb_zt < -gr_grid_dir_indx + gr_k_ub_zt + 1
                              : gr_grid_dir_indx + gr_k_lb_zt >= -gr_grid_dir_indx + gr_k_ub_zt)
            ? (((-2 * gr_grid_dir_indx - gr_k_lb_zt + gr_k_ub_zt) / gr_grid_dir_indx) + 1)
            : 0;
    Kokkos::parallel_for(
        "compute_mixing_length::cells_2", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int j = (gr_grid_dir_indx * (cell / ngrdcol) + gr_grid_dir_indx + gr_k_lb_zt);
          int j_zm{-1};
          if ((gr_grid_dir_indx > 0)) {
            // Ascending grid
            j_zm = j;
          } else {
            // gr%grid_dir_indx < 0
            // Descending grid
            j_zm = (j + 1);
          }
          thl_par_j_precalc(i, j) =
              ((thlm(i, j) - (thlm(i, (-gr_grid_dir_indx + j)) * exp_mu_dzm(i, j_zm))) -
               ((thlm(i, j) - thlm(i, (-gr_grid_dir_indx + j))) * entrain_coef(i, j_zm)));
          rt_par_j_precalc(i, j) =
              ((rtm(i, j) - (rtm(i, (-gr_grid_dir_indx + j)) * exp_mu_dzm(i, j_zm))) -
               ((rtm(i, j) - rtm(i, (-gr_grid_dir_indx + j))) * entrain_coef(i, j_zm)));
        });
  }
  // Calculate the initial change in TKE for each level. This is done for computational
  // efficiency, it helps because there will be at least one calculation for each grid level,
  // meaning the first one can be done for every grid level and therefore the calculations can
  // be vectorized, clubb:ticket:834. After the initial calculation however, it is uncertain
  // how many more iterations should be done for each individual grid level, and calculating
  // one change in TKE for each level until all are exhausted will result in many unnessary
  // and expensive calculations.
  // Calculate initial thl, tl, and rt for parcels at each grid level
  Kokkos::deep_copy(thl_par_1, zero);
  Kokkos::deep_copy(tl_par_1, zero);
  Kokkos::deep_copy(rt_par_1, zero);
  {
    const int level_count =
        (gr_grid_dir_indx > 0 ? gr_grid_dir_indx + gr_k_lb_zt < -gr_grid_dir_indx + gr_k_ub_zt + 1
                              : gr_grid_dir_indx + gr_k_lb_zt >= -gr_grid_dir_indx + gr_k_ub_zt)
            ? (((-2 * gr_grid_dir_indx - gr_k_lb_zt + gr_k_ub_zt) / gr_grid_dir_indx) + 1)
            : 0;
    Kokkos::parallel_for(
        "compute_mixing_length::cells_3", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int j = (gr_grid_dir_indx * (cell / ngrdcol) + gr_grid_dir_indx + gr_k_lb_zt);
          int j_zm{-1};
          if ((gr_grid_dir_indx > 0)) {
            // Ascending grid
            j_zm = j;
          } else {
            // gr%grid_dir_indx < 0
            // Descending grid
            j_zm = (j + 1);
          }
          thl_par_1(i, j) = (thlm(i, j) - ((thlm(i, j) - thlm(i, (-gr_grid_dir_indx + j))) *
                                           entrain_coef(i, j_zm)));
          tl_par_1(i, j)  = (thl_par_1(i, j) * exner(i, j));
          rt_par_1(i, j) =
              (rtm(i, j) - ((rtm(i, j) - rtm(i, (-gr_grid_dir_indx + j))) * entrain_coef(i, j_zm)));
        });
  }
  // Caclculate initial rsatl for parcels at each grid level
  // The entire pressure and temperature arrays are passed as
  // argument and the sub-arrays are choosen using
  // start_index. This workaround is used to solve
  // subarray issues with OpenACC.
  // rsatl_par_1(i,3:) = sat_mixrat_liq_acc( nz-2, ngrdcol, p_in_Pa(i,3:), tl_par_1(i,3:) )
  // since subarray 3:, the start_index is 3 and it is an optional argument
  start_index = (gr_grid_dir_indx + gr_k_lb_zt);
  sat_mixrat_liq_api(nzt, ngrdcol, gr, p_in_Pa, tl_par_1, saturation_formula, rsatl_par_1,
                     start_index);
  // Calculate initial dCAPE_dz and CAPE_incr for parcels at each grid level
  Kokkos::parallel_for(
      "compute_mixing_length::columns_4", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column;
        int j{-1};
        int k{-1};
        int j_zm{-1};
        int kp1_zm{-1};
        real_t tke{};
        real_t CAPE_incr{};
        real_t dCAPE_dz_j{};
        real_t dCAPE_dz_j_minus_1{};
        real_t thl_par_j{};
        real_t rt_par_j{};
        real_t rc_par_j{};
        real_t thv_par_j{};
        real_t tl_par_j{};
        real_t rsatl_par_j{};
        real_t s_par_j{};
        real_t Lscale_up_max_alt{};
        real_t tl_par_j_sqd{};
        real_t invrs_dCAPE_diff{};
        for (j = (gr_grid_dir_indx + gr_k_lb_zt);
             (gr_grid_dir_indx > 0 ? j < gr_k_ub_zt + 1 : j >= gr_k_ub_zt); j += gr_grid_dir_indx) {
          if ((gr_grid_dir_indx > 0)) {
            // Ascending grid
            j_zm = j;
          } else {
            // gr%grid_dir_indx < 0
            // Descending grid
            j_zm = (j + 1);
          }
          tl_par_j_sqd = (tl_par_1(i, j) * tl_par_1(i, j));
          // s from Lewellen and Yoh 1993 (LY) eqn. 1
          // s = ( rt - rsatl ) / ( 1 + beta * rsatl )
          // and SD's beta (eqn. 8),
          // beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
          // Simplified by multiplying top and bottom by tl^2 to avoid a divide and precalculating
          // ep * Lv**2 / ( Rd * cp )
          s_par_1(i, j)  = (((rt_par_1(i, j) - rsatl_par_1(i, j)) * tl_par_j_sqd) /
                            (tl_par_j_sqd + (Lv2_coef * rsatl_par_1(i, j))));
          rc_par_1(i, j) = (s_par_1(i, j) > zero_threshold ? s_par_1(i, j) : zero_threshold);
          // theta_v of entraining parcel at grid level j
          thv_par_1(i, j) = ((thl_par_1(i, j) + ((ep1 * thv_ds(i, j)) * rt_par_1(i, j))) +
                             (Lv_coef(i, j) * rc_par_1(i, j)));
          // dCAPE/dz = g * ( thv_par - thvm ) / thvm.
          dCAPE_dz_1(i, j) = (grav_on_thvm(i, j) * (thv_par_1(i, j) - thvm(i, j)));
          // CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
          // Trapezoidal estimate between grid levels, dCAPE at z_0 = 0 for this initial calculation
          CAPE_incr_1(i, j) = (((one_half * dCAPE_dz_1(i, j)) * gr_grid_dir) * gr_dzm(i, j_zm));
        }
        // Calculate Lscale_up for each grid level. If the TKE from a parcel has not been completely
        // exhausted by the initial change then continue the exhaustion calculations here for a
        // single
        // grid level at a time until the TKE is exhausted.
        Lscale_up_max_alt = zero;
        // Set initial max value for Lscale_up to 0
        for (k = gr_k_lb_zt; (gr_grid_dir_indx > 0 ? k < -2 * gr_grid_dir_indx + gr_k_ub_zt + 1
                                                   : k >= -2 * gr_grid_dir_indx + gr_k_ub_zt);
             k += gr_grid_dir_indx) {
          // If the initial turbulent kinetic energy (tke) has not been exhausted for this grid
          // level
          if (((tke_i(i, k) + CAPE_incr_1(i, (gr_grid_dir_indx + k))) > zero)) {
            // Calculate new TKE for parcel
            tke = (tke_i(i, k) + CAPE_incr_1(i, (gr_grid_dir_indx + k)));
            // Set j to 2 levels above current Lscale_up level, this is because we've already
            // determined that the parcel can rise at least 1 full level
            j = (2 * gr_grid_dir_indx + k);
            // Set initial thl, rt, and dCAPE_dz to the values found by the intial calculations
            thl_par_j          = thl_par_1(i, (gr_grid_dir_indx + k));
            rt_par_j           = rt_par_1(i, (gr_grid_dir_indx + k));
            dCAPE_dz_j_minus_1 = dCAPE_dz_1(i, (gr_grid_dir_indx + k));
            // Continue change in TKE calculations until it is exhausted or the max grid
            // level has been reached. j is the next grid level above the level that can
            // be reached for a parcel starting at level k. If TKE is exhausted in this loop
            // that means the parcel starting at k cannot reach level j, but has reached j-1
            // Update for generalized gridding:
            // Ascending grid: do while ( j < nzt ) -- j is increasing
            // Descending grid: do while ( j > 1 ) -- j is decreasing
            // Either way, j is going upwards in altitude.
            // Reconciling the grids:
            // Step 1: Use gr%k_ub_zt for the upper boundary:
            // Ascending grid: do while ( j < gr%k_ub_zt ) -- j is increasing
            // Descending grid: do while ( j > gr%k_ub_zt ) -- j is decreasing
            // Step 2: Reconciling the greater than / less than:
            // Ascending grid: do while ( j < nzt ) -- j is increasing
            // Descending grid: do while ( -j < -1 ) -- j is increasing
            // Generalized Expression:
            // do while ( gr%grid_dir_indx * j < gr%grid_dir_indx * gr%k_ub_zt );
            // where for an ascending grid, gr%grid_dir_indx = 1 and gr%k_ub_zt = nzt;
            // and for a descending grid, gr%grid_dir_indx = -1 and gr%k_ub_zt = 1.
            while ((gr_grid_dir_indx * j + gr_grid_dir_indx <
                    gr_grid_dir_indx * gr_k_ub_zt + gr_grid_dir_indx)) {
              if ((gr_grid_dir_indx > 0)) {
                // Ascending grid
                j_zm = j;
              } else {
                // gr%grid_dir_indx < 0
                // Descending grid
                j_zm = (j + 1);
              }
              // thl, rt of parcel are conserved except for entrainment
              // The values of thl_env and rt_env are treated as changing linearly for a parcel
              // of air ascending from level j-1 to level j
              // theta_l of the parcel starting at grid level k, and currenly
              // at grid level j
              // d(thl_par)/dz = - Lscale_mu_entrain * ( thl_par - thl_env )
              thl_par_j = (thl_par_j_precalc(i, j) + (thl_par_j * exp_mu_dzm(i, j_zm)));
              // r_t of the parcel starting at grid level k, and currenly
              // at grid level j
              // d(rt_par)/dz = - Lscale_mu_entrain * ( rt_par - rt_env )
              rt_par_j = (rt_par_j_precalc(i, j) + (rt_par_j * exp_mu_dzm(i, j_zm)));
              // Include effects of latent heating on Lscale_up 6/12/00
              // Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
              // Probably should use properties of bump 1 in Gaussian, not mean!!!
              tl_par_j     = (thl_par_j * exner(i, j));
              rsatl_par_j  = sat_mixrat_liq_api(p_in_Pa(i, j), tl_par_j, saturation_formula);
              tl_par_j_sqd = (tl_par_j * tl_par_j);
              // s from Lewellen and Yoh 1993 (LY) eqn. 1
              // s = ( rt - rsatl ) / ( 1 + beta * rsatl )
              // and SD's beta (eqn. 8),
              // beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
              // Simplified by multiplying top and bottom by tl^2 to avoid a
              // divide and precalculating ep * Lv**2 / ( Rd * cp )
              s_par_j  = (((rt_par_j - rsatl_par_j) * tl_par_j_sqd) /
                          (tl_par_j_sqd + (Lv2_coef * rsatl_par_j)));
              rc_par_j = (s_par_j > zero_threshold ? s_par_j : zero_threshold);
              // theta_v of entraining parcel at grid level j
              thv_par_j =
                  ((thl_par_j + ((ep1 * thv_ds(i, j)) * rt_par_j)) + (Lv_coef(i, j) * rc_par_j));
              // dCAPE/dz = g * ( thv_par - thvm ) / thvm.
              dCAPE_dz_j = (grav_on_thvm(i, j) * (thv_par_j - thvm(i, j)));
              // CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
              // Trapezoidal estimate between grid levels j and j-1
              CAPE_incr = (((one_half * (dCAPE_dz_j + dCAPE_dz_j_minus_1)) * gr_grid_dir) *
                           gr_dzm(i, j_zm));
              // Exit loop early if tke has been exhaused between level j and j+1
              if (((tke + CAPE_incr) <= zero)) {
                break;
              }
              // Save previous dCAPE value for next cycle
              dCAPE_dz_j_minus_1 = dCAPE_dz_j;
              // Caclulate new TKE and increment j
              tke = (tke + CAPE_incr);
              j   = (gr_grid_dir_indx + j);
            }
            // Add full grid level thickness for each grid level that was passed without the TKE
            // being exhausted, difference between starting level (k) and last level passed (j-1)
            Lscale_up(i, k) = ((Lscale_up(i, k) + gr_zt(i, (-gr_grid_dir_indx + j))) - gr_zt(i, k));
            if ((gr_grid_dir_indx * j + gr_grid_dir_indx <
                 gr_grid_dir_indx * gr_k_ub_zt + gr_grid_dir_indx)) {
              // Loop terminated early, meaning TKE was completely exhaused at grid level j.
              // Add the thickness z - z_0 (where z_0 < z <= z_1) to Lscale_up.
              if (((Kokkos::abs((dCAPE_dz_j - dCAPE_dz_j_minus_1)) * 2) <=
                   (Kokkos::abs((dCAPE_dz_j + dCAPE_dz_j_minus_1)) * eps))) {
                // Special case where dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) = 0
                // Find the remaining distance z - z_0 that it takes to
                // exhaust the remaining TKE
                Lscale_up(i, k) = (Lscale_up(i, k) + (-(tke) / dCAPE_dz_j));
              } else {
                if ((gr_grid_dir_indx > 0)) {
                  // Ascending grid
                  j_zm = j;
                } else {
                  // gr%grid_dir_indx < 0
                  // Descending grid
                  j_zm = (j + 1);
                }
                // Case used for most scenarios where dCAPE/dz|_(z_1) /= dCAPE/dz|_(z_0)
                // Find the remaining distance z - z_0 that it takes to exhaust the
                // remaining TKE (tke_i), using the quadratic formula (only the
                // negative (-) root works in this scenario).
                invrs_dCAPE_diff = (one / (dCAPE_dz_j - dCAPE_dz_j_minus_1));
                Lscale_up(i, k) =
                    ((Lscale_up(i, k) -
                      (((dCAPE_dz_j_minus_1 * invrs_dCAPE_diff) * gr_grid_dir) * gr_dzm(i, j_zm))) -
                     (((Kokkos::sqrt(((dCAPE_dz_j_minus_1 * dCAPE_dz_j_minus_1) -
                                      ((((two * tke) * gr_grid_dir) * gr_invrs_dzm(i, j_zm)) *
                                       (dCAPE_dz_j - dCAPE_dz_j_minus_1)))) *
                        invrs_dCAPE_diff) *
                       gr_grid_dir) *
                      gr_dzm(i, j_zm)));
              }
            }
          } else {
            // TKE for parcel at level (k) was exhaused before one full grid level
            if ((gr_grid_dir_indx > 0)) {
              // Ascending grid
              kp1_zm = (k + 1);
            } else {
              // gr%grid_dir_indx < 0
              // Descending grid
              kp1_zm = k;
            }
            // Find the remaining distance z - z_0 that it takes to exhaust the
            // remaining TKE (tke_i), using the quadratic formula. Simplified
            // since dCAPE_dz_j_minus_1 = 0.0
            Lscale_up(i, k) =
                (Lscale_up(i, k) -
                 (Kokkos::sqrt(((((-(two)*tke_i(i, k)) * gr_grid_dir) * gr_dzm(i, kp1_zm)) *
                                dCAPE_dz_1(i, (gr_grid_dir_indx + k)))) /
                  dCAPE_dz_1(i, (gr_grid_dir_indx + k))));
          }
          // If a parcel at a previous grid level can rise past the parcel at this grid level
          // then this one should also be able to rise up to that height. This feature insures
          // that the profile of Lscale_up will be smooth, thus reducing numerical instability.
          if (((gr_zt(i, k) + Lscale_up(i, k)) < Lscale_up_max_alt)) {
            // A lower starting parcel can ascend higher than this one, set height to the max
            // that any lower starting parcel can ascend to
            Lscale_up(i, k) = (Lscale_up_max_alt - gr_zt(i, k));
          } else {
            // This parcel can ascend higher than any below it, save final height
            Lscale_up_max_alt = (Lscale_up(i, k) + gr_zt(i, k));
          }
        }
      });
  // ---------------- Downwards Length Scale Calculation ----------------
  // Precalculate values for downward Lscale, these are useful only if a parcel can descend
  // more than one level. They are used in the equations that calculate thl and rt
  // recursively for a parcel as it descends
  {
    const int level_count =
        (gr_grid_dir_indx > 0 ? gr_k_lb_zt < -gr_grid_dir_indx + gr_k_ub_zt + 1
                              : gr_k_lb_zt >= -gr_grid_dir_indx + gr_k_ub_zt)
            ? (((-gr_grid_dir_indx - gr_k_lb_zt + gr_k_ub_zt) / gr_grid_dir_indx) + 1)
            : 0;
    Kokkos::parallel_for(
        "compute_mixing_length::cells_5", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int j = (gr_grid_dir_indx * (cell / ngrdcol) + gr_k_lb_zt);
          int jp1_zm{-1};
          if ((gr_grid_dir_indx > 0)) {
            // Ascending grid
            jp1_zm = (j + 1);
          } else {
            // gr%grid_dir_indx < 0
            // Descending grid
            jp1_zm = j;
          }
          thl_par_j_precalc(i, j) =
              ((thlm(i, j) - (thlm(i, (gr_grid_dir_indx + j)) * exp_mu_dzm(i, jp1_zm))) -
               ((thlm(i, j) - thlm(i, (gr_grid_dir_indx + j))) * entrain_coef(i, jp1_zm)));
          rt_par_j_precalc(i, j) =
              ((rtm(i, j) - (rtm(i, (gr_grid_dir_indx + j)) * exp_mu_dzm(i, jp1_zm))) -
               ((rtm(i, j) - rtm(i, (gr_grid_dir_indx + j))) * entrain_coef(i, jp1_zm)));
        });
  }
  // Calculate the initial change in TKE for each level. This is done for computational
  // efficiency, it helps because there will be at least one calculation for each grid level,
  // meaning the first one can be done for every grid level and therefore the calculations can
  // be vectorized, clubb:ticket:834. After the initial calculation however, it is uncertain
  // how many more iterations should be done for each individual grid level, and calculating
  // one change in TKE for each level until all are exhausted will result in many unnessary
  // and expensive calculations.
  // Calculate initial thl, tl, and rt for parcels at each grid level
  {
    const int level_count =
        (gr_grid_dir_indx > 0 ? gr_k_lb_zt < -gr_grid_dir_indx + gr_k_ub_zt + 1
                              : gr_k_lb_zt >= -gr_grid_dir_indx + gr_k_ub_zt)
            ? (((-gr_grid_dir_indx - gr_k_lb_zt + gr_k_ub_zt) / gr_grid_dir_indx) + 1)
            : 0;
    Kokkos::parallel_for(
        "compute_mixing_length::cells_6", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int j = (gr_grid_dir_indx * (cell / ngrdcol) + gr_k_lb_zt);
          int jp1_zm{-1};
          if ((gr_grid_dir_indx > 0)) {
            // Ascending grid
            jp1_zm = (j + 1);
          } else {
            // gr%grid_dir_indx < 0
            // Descending grid
            jp1_zm = j;
          }
          thl_par_1(i, j) = (thlm(i, j) - ((thlm(i, j) - thlm(i, (gr_grid_dir_indx + j))) *
                                           entrain_coef(i, jp1_zm)));
          tl_par_1(i, j)  = (thl_par_1(i, j) * exner(i, j));
          rt_par_1(i, j)  = (rtm(i, j) - ((rtm(i, j) - rtm(i, (gr_grid_dir_indx + j))) *
                                          entrain_coef(i, jp1_zm)));
        });
  }
  // Caclculate initial rsatl for parcels at each grid level, this function is elemental
  // The entire pressure and temperature arrays are passed as
  // argument and the sub-arrays are choosen using
  // start_index. This workaround is used to solve
  // subarray issues with OpenACC.
  // rsatl_par_1(i,2:) = sat_mixrat_liq_acc( nz-1, p_in_Pa(i,2:), tl_par_1(i,2:) )
  // since subarray 2:, the start_index is 2 and it is an optional argument
  start_index = gr_k_lb_zt;
  sat_mixrat_liq_api(nzt, ngrdcol, gr, p_in_Pa, tl_par_1, saturation_formula, rsatl_par_1,
                     start_index);
  // Calculate initial dCAPE_dz and CAPE_incr for parcels at each grid level
  Kokkos::parallel_for(
      "compute_mixing_length::columns_7", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column;
        int j{-1};
        int k{-1};
        int jp1_zm{-1};
        int k_zm{-1};
        real_t tke{};
        real_t CAPE_incr{};
        real_t dCAPE_dz_j{};
        real_t dCAPE_dz_j_plus_1{};
        real_t thl_par_j{};
        real_t rt_par_j{};
        real_t rc_par_j{};
        real_t thv_par_j{};
        real_t tl_par_j{};
        real_t rsatl_par_j{};
        real_t s_par_j{};
        real_t Lscale_down_min_alt{};
        real_t tl_par_j_sqd{};
        real_t invrs_dCAPE_diff{};
        for (j = gr_k_lb_zt; (gr_grid_dir_indx > 0 ? j < -gr_grid_dir_indx + gr_k_ub_zt + 1
                                                   : j >= -gr_grid_dir_indx + gr_k_ub_zt);
             j += gr_grid_dir_indx) {
          if ((gr_grid_dir_indx > 0)) {
            // Ascending grid
            jp1_zm = (j + 1);
          } else {
            // gr%grid_dir_indx < 0
            // Descending grid
            jp1_zm = j;
          }
          tl_par_j_sqd = (tl_par_1(i, j) * tl_par_1(i, j));
          // s from Lewellen and Yoh 1993 (LY) eqn. 1
          // s = ( rt - rsatl ) / ( 1 + beta * rsatl )
          // and SD's beta (eqn. 8),
          // beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
          // Simplified by multiplying top and bottom by tl^2 to avoid a divide and precalculating
          // ep * Lv**2 / ( Rd * cp )
          s_par_1(i, j)  = (((rt_par_1(i, j) - rsatl_par_1(i, j)) * tl_par_j_sqd) /
                            (tl_par_j_sqd + (Lv2_coef * rsatl_par_1(i, j))));
          rc_par_1(i, j) = (s_par_1(i, j) > zero_threshold ? s_par_1(i, j) : zero_threshold);
          // theta_v of entraining parcel at grid level j
          thv_par_1(i, j) = ((thl_par_1(i, j) + ((ep1 * thv_ds(i, j)) * rt_par_1(i, j))) +
                             (Lv_coef(i, j) * rc_par_1(i, j)));
          // dCAPE/dz = g * ( thv_par - thvm ) / thvm.
          dCAPE_dz_1(i, j) = (grav_on_thvm(i, j) * (thv_par_1(i, j) - thvm(i, j)));
          // CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
          // Trapezoidal estimate between grid levels, dCAPE at z_0 = 0 for this initial calculation
          CAPE_incr_1(i, j) = (((one_half * dCAPE_dz_1(i, j)) * gr_grid_dir) * gr_dzm(i, jp1_zm));
        }
        // Calculate Lscale_down for each grid level. If the TKE from a parcel has not been
        // completely
        // exhausted by the initial change then continue the exhaustion calculations here for a
        // single
        // grid level at a time until the TKE is exhausted.
        Lscale_down_min_alt = gr_zt(i, gr_k_ub_zt);
        // Set initial min value for Lscale_down to max zt
        for (k = gr_k_ub_zt; (-1 * gr_grid_dir_indx > 0 ? k < gr_grid_dir_indx + gr_k_lb_zt + 1
                                                        : k >= gr_grid_dir_indx + gr_k_lb_zt);
             k += (-1 * gr_grid_dir_indx)) {
          // If the initial turbulent kinetic energy (tke) has not been exhausted for this grid
          // level
          if (((tke_i(i, k) - CAPE_incr_1(i, (-gr_grid_dir_indx + k))) > zero)) {
            // Calculate new TKE for parcel
            tke = (tke_i(i, k) - CAPE_incr_1(i, (-gr_grid_dir_indx + k)));
            // Set j to 2 levels below current Lscale_down level, this is because we've already
            // determined that the parcel can descend at least 1 full level
            j = (-2 * gr_grid_dir_indx + k);
            // Set initial thl, rt, and dCAPE_dz to the values found by the intial calculations
            thl_par_j         = thl_par_1(i, (-gr_grid_dir_indx + k));
            rt_par_j          = rt_par_1(i, (-gr_grid_dir_indx + k));
            dCAPE_dz_j_plus_1 = dCAPE_dz_1(i, (-gr_grid_dir_indx + k));
            // Continue change in TKE calculations until it is exhausted or the min grid
            // level has been reached. j is the next grid level below the level that can
            // be reached for a parcel starting at level k. If TKE is exhausted in this loop
            // that means the parcel starting at k cannot sink to level j, but can sink to j+1
            // Update for generalized gridding:
            // Ascending grid: do while ( j >= 1 ) -- j is decreasing
            // Descending grid: do while ( j <= nzt ) -- j is increasing
            // Either way, j is going downwards in altitude.
            // Reconciling the grids:
            // Step 1: Use gr%k_lb_zt for the lower boundary:
            // Ascending grid: do while ( j >= gr%k_lb_zt ) -- j is decreasing
            // Descending grid: do while ( j <= gr%k_lb_zt ) -- j is increasing
            // Step 2: Reconciling the greater than / less than:
            // Ascending grid: do while ( j >= 1 ) -- j is decreasing
            // Descending grid: do while ( -j >= -nzt ) -- j is decreasing
            // Generalized Expression:
            // do while ( gr%grid_dir_indx * j >= gr%grid_dir_indx * gr%k_lb_zt );
            // where for an ascending grid, gr%grid_dir_indx = 1 and gr%k_lb_zt = 1;
            // and for a descending grid, gr%grid_dir_indx = -1 and gr%k_lb_zt = nzt.
            while ((gr_grid_dir_indx * j + gr_grid_dir_indx >=
                    gr_grid_dir_indx * gr_k_lb_zt + gr_grid_dir_indx)) {
              if ((gr_grid_dir_indx > 0)) {
                // Ascending grid
                jp1_zm = (j + 1);
              } else {
                // gr%grid_dir_indx < 0
                // Descending grid
                jp1_zm = j;
              }
              // thl, rt of parcel are conserved except for entrainment
              // The values of thl_env and rt_env are treated as changing linearly for a parcel
              // of air descending from level j to level j-1
              // theta_l of the parcel starting at grid level k, and currenly
              // at grid level j
              // d(thl_par)/dz = - Lscale_mu_entrain * ( thl_par - thl_env )
              thl_par_j = (thl_par_j_precalc(i, j) + (thl_par_j * exp_mu_dzm(i, jp1_zm)));
              // r_t of the parcel starting at grid level k, and currenly
              // at grid level j
              // d(rt_par)/dz = - Lscale_mu_entrain * ( rt_par - rt_env )
              rt_par_j = (rt_par_j_precalc(i, j) + (rt_par_j * exp_mu_dzm(i, jp1_zm)));
              // Include effects of latent heating on Lscale_up 6/12/00
              // Use thermodynamic formula of Bougeault 1981 JAS Vol. 38, 2416
              // Probably should use properties of bump 1 in Gaussian, not mean!!!
              tl_par_j     = (thl_par_j * exner(i, j));
              rsatl_par_j  = sat_mixrat_liq_api(p_in_Pa(i, j), tl_par_j, saturation_formula);
              tl_par_j_sqd = (tl_par_j * tl_par_j);
              // s from Lewellen and Yoh 1993 (LY) eqn. 1
              // s = ( rt - rsatl ) / ( 1 + beta * rsatl )
              // and SD's beta (eqn. 8),
              // beta = ep * ( Lv / ( Rd * tl ) ) * ( Lv / ( cp * tl ) )
              // Simplified by multiplying top and bottom by tl^2 to avoid a
              // divide and precalculating ep * Lv**2 / ( Rd * cp )
              s_par_j  = (((rt_par_j - rsatl_par_j) * tl_par_j_sqd) /
                          (tl_par_j_sqd + (Lv2_coef * rsatl_par_j)));
              rc_par_j = (s_par_j > zero_threshold ? s_par_j : zero_threshold);
              // theta_v of entraining parcel at grid level j
              thv_par_j =
                  ((thl_par_j + ((ep1 * thv_ds(i, j)) * rt_par_j)) + (Lv_coef(i, j) * rc_par_j));
              // dCAPE/dz = g * ( thv_par - thvm ) / thvm.
              dCAPE_dz_j = (grav_on_thvm(i, j) * (thv_par_j - thvm(i, j)));
              // CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz
              // Trapezoidal estimate between grid levels j+1 and j
              CAPE_incr = (((one_half * (dCAPE_dz_j + dCAPE_dz_j_plus_1)) * gr_grid_dir) *
                           gr_dzm(i, jp1_zm));
              // Exit loop early if tke has been exhaused between level j+1 and j
              if (((tke - CAPE_incr) <= zero)) {
                break;
              }
              // Save previous dCAPE value for next cycle
              dCAPE_dz_j_plus_1 = dCAPE_dz_j;
              // Caclulate new TKE and increment j
              tke = (tke - CAPE_incr);
              j   = (-gr_grid_dir_indx + j);
            }
            // Add full grid level thickness for each grid level that was passed without the TKE
            // being exhausted, difference between starting level (k) and last level passed (j+1)
            Lscale_down(i, k) =
                ((Lscale_down(i, k) + gr_zt(i, k)) - gr_zt(i, (gr_grid_dir_indx + j)));
            if ((gr_grid_dir_indx * j + gr_grid_dir_indx >=
                 gr_grid_dir_indx * gr_k_lb_zt + gr_grid_dir_indx)) {
              // Loop terminated early, meaning TKE was completely exhaused at grid level j.
              // Add the thickness z - z_0 (where z_0 < z <= z_1) to Lscale_up.
              if (((Kokkos::abs((dCAPE_dz_j - dCAPE_dz_j_plus_1)) * 2) <=
                   (Kokkos::abs((dCAPE_dz_j + dCAPE_dz_j_plus_1)) * eps))) {
                // Special case where dCAPE/dz|_(z_(-1)) - dCAPE/dz|_(z_0) = 0
                // Find the remaining distance z_0 - z that it takes to
                // exhaust the remaining TKE
                Lscale_down(i, k) = (Lscale_down(i, k) + (tke / dCAPE_dz_j));
              } else {
                if ((gr_grid_dir_indx > 0)) {
                  // Ascending grid
                  jp1_zm = (j + 1);
                } else {
                  // gr%grid_dir_indx < 0
                  // Descending grid
                  jp1_zm = j;
                }
                // Case used for most scenarios where dCAPE/dz|_(z_(-1)) /= dCAPE/dz|_(z_0)
                // Find the remaining distance z_0 - z that it takes to exhaust the
                // remaining TKE (tke_i), using the quadratic formula (only the
                // negative (-) root works in this scenario) -- however, the
                // negative (-) root is divided by another negative (-) factor,
                // which results in an overall plus (+) sign in front of the
                // square root term in the equation below).
                invrs_dCAPE_diff = (one / (dCAPE_dz_j - dCAPE_dz_j_plus_1));
                Lscale_down(i, k) =
                    ((Lscale_down(i, k) - (((dCAPE_dz_j_plus_1 * invrs_dCAPE_diff) * gr_grid_dir) *
                                           gr_dzm(i, jp1_zm))) +
                     (((Kokkos::sqrt(((dCAPE_dz_j_plus_1 * dCAPE_dz_j_plus_1) +
                                      ((((two * tke) * gr_grid_dir) * gr_invrs_dzm(i, jp1_zm)) *
                                       (dCAPE_dz_j - dCAPE_dz_j_plus_1)))) *
                        invrs_dCAPE_diff) *
                       gr_grid_dir) *
                      gr_dzm(i, jp1_zm)));
              }
            }
          } else {
            // TKE for parcel at level (k) was exhaused before one full grid level
            if ((gr_grid_dir_indx > 0)) {
              // Ascending grid
              k_zm = k;
            } else {
              // gr%grid_dir_indx < 0
              // Descending grid
              k_zm = (k + 1);
            }
            // Find the remaining distance z_0 - z that it takes to exhaust the
            // remaining TKE (tke_i), using the quadratic formula. Simplified
            // since dCAPE_dz_j_plus_1 = 0.0
            Lscale_down(i, k) =
                (Lscale_down(i, k) +
                 (Kokkos::sqrt(((((two * tke_i(i, k)) * gr_grid_dir) * gr_dzm(i, k_zm)) *
                                dCAPE_dz_1(i, (-gr_grid_dir_indx + k)))) /
                  dCAPE_dz_1(i, (-gr_grid_dir_indx + k))));
          }
          // If a parcel at a previous grid level can descend past the parcel at this grid level
          // then this one should also be able to descend down to that height. This feature insures
          // that the profile of Lscale_down will be smooth, thus reducing numerical instability.
          if (((gr_zt(i, k) - Lscale_down(i, k)) > Lscale_down_min_alt)) {
            Lscale_down(i, k) = (gr_zt(i, k) - Lscale_down_min_alt);
          } else {
            Lscale_down_min_alt = (gr_zt(i, k) - Lscale_down(i, k));
          }
        }
      });
  // ---------------- Final Lscale Calculation ----------------
  Kokkos::parallel_for(
      "compute_mixing_length::columns_8", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column;
        int k{-1};
        real_t lminh{};
        for (k = 0; k < nzt; k += 1) {
          // Make lminh a linear function starting at value lmin at the bottom
          // and going to zero at 500 meters in altitude.
          if (l_implemented) {
            // Within a host model, increase mixing length in 500 m layer above *ground*
            lminh = (((zero_threshold > (Lscale_sfclyr_depth - (gr_zt(i, k) - gr_zm(i, gr_k_lb_zm)))
                           ? zero_threshold
                           : (Lscale_sfclyr_depth - (gr_zt(i, k) - gr_zm(i, gr_k_lb_zm)))) *
                      lmin(i)) *
                     invrs_Lscale_sfclyr_depth);
          } else {
            // In standalone mode, increase mixing length in 500 m layer above *mean sea level*
            lminh = (((zero_threshold > (Lscale_sfclyr_depth - gr_zt(i, k))
                           ? zero_threshold
                           : (Lscale_sfclyr_depth - gr_zt(i, k))) *
                      lmin(i)) *
                     invrs_Lscale_sfclyr_depth);
          }
          Lscale_up(i, k)   = (lminh > Lscale_up(i, k) ? lminh : Lscale_up(i, k));
          Lscale_down(i, k) = (lminh > Lscale_down(i, k) ? lminh : Lscale_down(i, k));
          // When L is large, turbulence is weakly damped
          // When L is small, turbulence is strongly damped
          // Use a geometric mean to determine final Lscale so that L tends to become small
          // if either Lscale_up or Lscale_down becomes small.
          Lscale(i, k) = Kokkos::sqrt((Lscale_up(i, k) * Lscale_down(i, k)));
        }
        // Set the value of Lscale at the upper and lower boundaries.
        Lscale(i, gr_k_ub_zt) = Lscale(i, (-gr_grid_dir_indx + gr_k_ub_zt));
        // Vince Larson limited Lscale to allow host
        // model to take over deep convection.  13 Feb 2008.
        for (int level = 0; level < nzt; ++level)
          Lscale(i, level) = (Lscale(i, level) < Lscale_max(i) ? Lscale(i, level) : Lscale_max(i));
      });
  // Ensure that no Lscale values are NaN
  if (clubb_at_least_debug_level_api(1)) {
    for (i = 0; i < ngrdcol; i += 1) {
      length_check(nzt, Kokkos::subview(Lscale, i, Kokkos::ALL), // In
                   Kokkos::subview(Lscale_up, i, Kokkos::ALL),   // In
                   Kokkos::subview(Lscale_down, i, Kokkos::ALL), // In
                   err_info                                      // InOut
      );
      if ((to_host_space(workspace, err_info.err_code)(i) == clubb_fatal_error)) {
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          const auto header_host = to_host_space(workspace, err_info.err_header);
          diagnostic.write(header_host.data() + (header_host.extent(0) * i), header_host.extent(0));
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "Errors in compute_mixing_length subroutine in grid column " << ' ';
          diagnostic << ((i + 1)) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "Intent(in)" << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "thvm = " << ' ';
          const auto field_host = to_host_space(workspace, Kokkos::subview(thvm, i, Kokkos::ALL));
          for (std::size_t element = 0; element < field_host.extent(0); ++element)
            diagnostic << field_host(element) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "thlm = " << ' ';
          const auto field_host = to_host_space(workspace, Kokkos::subview(thlm, i, Kokkos::ALL));
          for (std::size_t element = 0; element < field_host.extent(0); ++element)
            diagnostic << field_host(element) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "rtm = " << ' ';
          const auto field_host = to_host_space(workspace, Kokkos::subview(rtm, i, Kokkos::ALL));
          for (std::size_t element = 0; element < field_host.extent(0); ++element)
            diagnostic << field_host(element) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "em = " << ' ';
          const auto field_host = to_host_space(workspace, Kokkos::subview(em, i, Kokkos::ALL));
          for (std::size_t element = 0; element < field_host.extent(0); ++element)
            diagnostic << field_host(element) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "exner = " << ' ';
          const auto field_host = to_host_space(workspace, Kokkos::subview(exner, i, Kokkos::ALL));
          for (std::size_t element = 0; element < field_host.extent(0); ++element)
            diagnostic << field_host(element) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "p_in_Pa = " << ' ';
          const auto field_host =
              to_host_space(workspace, Kokkos::subview(p_in_Pa, i, Kokkos::ALL));
          for (std::size_t element = 0; element < field_host.extent(0); ++element)
            diagnostic << field_host(element) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "thv_ds = " << ' ';
          const auto field_host = to_host_space(workspace, Kokkos::subview(thv_ds, i, Kokkos::ALL));
          for (std::size_t element = 0; element < field_host.extent(0); ++element)
            diagnostic << field_host(element) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "Intent(out)" << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "Lscale = " << ' ';
          const auto field_host = to_host_space(workspace, Kokkos::subview(Lscale, i, Kokkos::ALL));
          for (std::size_t element = 0; element < field_host.extent(0); ++element)
            diagnostic << field_host(element) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "Lscale_up = " << ' ';
          const auto field_host =
              to_host_space(workspace, Kokkos::subview(Lscale_up, i, Kokkos::ALL));
          for (std::size_t element = 0; element < field_host.extent(0); ++element)
            diagnostic << field_host(element) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
        {
          std::ostringstream diagnostic;
          diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
          diagnostic << "Lscale_down = " << ' ';
          const auto field_host =
              to_host_space(workspace, Kokkos::subview(Lscale_down, i, Kokkos::ALL));
          for (std::size_t element = 0; element < field_host.extent(0); ++element)
            diagnostic << field_host(element) << ' ';
          std::cerr << diagnostic.str() << '\n';
        }
      }
      // Fatal error
    }
  }
  return;
}

} // namespace

void calc_Lscale_directly(ScratchWorkspace& workspace, int ngrdcol, int nzm, int nzt,
                          const grid& gr,                                               // In
                          bool l_implemented, Real2D p_in_Pa, Real2D exner, Real2D rtm, // In
                          Real2D thlm, Real2D thvm, Real2D rtp2_zt, Real2D thlp2_zt,
                          Real2D rtpthlp_zt,                                            // In
                          const pdf_parameter& pdf_params, Real2D em, Real2D thv_ds_zt, // In
                          Real1D Lscale_max, Real1D lmin,                               // In
                          Real2D clubb_params,                                          // In
                          int saturation_formula,                                       // In
                          bool l_Lscale_plume_centered,                                 // In
                          int Lscale_method,                                            // In
                          stats_type& stats, err_info_type& err_info,                   // InOut
                          Real2D Lscale, Real2D Lscale_up, Real2D Lscale_down           // Out
) {
  // Description:
  // Select the legacy or coherent parcel calculation (methods 1 or 3)
  // using environmental profiles, momentum-level TKE, PDF properties, and
  // per-column parameters. Return thermodynamic-level combined, upward,
  // and downward lengths [m], updating stats and err_info.
  // If source constant l_avg_Lscale is enabled, compute perturbed launches
  // before the mean launch so the directional outputs describe the mean.
  // -----------------------------------------------------------------------

  // Arguments (Native order)
  // ngrdcol: Number of grid columns
  // nzm: Number of momentum levels
  // nzt: Number of thermodynamic levels
  // gr: Grid structure
  // l_implemented: True if CLUBB is being implemented and run in a host model
  // p_in_Pa: Air pressure (thermodynamic levels)       [Pa]
  // exner: Exner function [-]
  // rtm: Mean total-water mixing ratio [kg/kg]
  // thlm: Mean liquid-water potential temperature [K]
  // thvm: Mean virtual potential temperature [K]
  // rtp2_zt: Total-water variance [(kg/kg)^2]
  // thlp2_zt: Liquid-water potential-temperature variance [K^2]
  // rtpthlp_zt: Total-water/temperature covariance [K kg/kg]
  // pdf_params: PDF Parameters  [units vary]
  // em: Turbulent kinetic energy [m^2/s^2]
  // thv_ds_zt: Reference-state virtual potential temperature [K]
  // Lscale_max: Per-column cap on combined mixing length [m]
  // lmin: Surface value of the height-dependent length floor [m]
  // clubb_params: Array of CLUBB's tunable parameters    [units vary]
  // saturation_formula: Integer that stores the saturation formula to be used
  // l_Lscale_plume_centered: Alternate that uses the PDF to compute the perturbed values
  // Lscale_method: Mixing-length formulation
  // stats: Statistics accumulator
  // err_info: err_info struct containing err_code and err_header
  // Lscale: Mixing length      [m]
  // Lscale_up: Mixing length up   [m]
  // Lscale_down: Mixing length down [m]

  // Lscale is calculated in subroutine compute_mixing_length. If l_avg_Lscale
  // is true, compute_mixing_length is called two additional times with
  // perturbed values of rtm and thlm. An average Lscale from the three calls
  // is then calculated. This reduces temporal noise in RICO, BOMEX, LBA and
  // other cases, as described in Native's parameter declaration.
  auto workspace_scope        = workspace.scope();
  constexpr bool l_avg_Lscale = false;
  Real2D sign_rtpthlp_zt =
      workspace_scope.view<Real2D>("calc_Lscale_directly::sign_rtpthlp_zt", ngrdcol,
                                   nzt); // Sign of the covariance rtpthlp       [-]
  Real2D Lscale_pert_1 = workspace_scope.view<Real2D>(
      "calc_Lscale_directly::Lscale_pert_1", ngrdcol, nzt); // For avg. calculation of Lscale  [m]
  Real2D Lscale_pert_2 =
      workspace_scope.view<Real2D>("calc_Lscale_directly::Lscale_pert_2", ngrdcol, nzt);
  Real2D thlm_pert_1 = workspace_scope.view<Real2D>("calc_Lscale_directly::thlm_pert_1", ngrdcol,
                                                    nzt); // For avg. calculation of Lscale  [K]
  Real2D thlm_pert_2 =
      workspace_scope.view<Real2D>("calc_Lscale_directly::thlm_pert_2", ngrdcol, nzt);
  Real2D rtm_pert_1 = workspace_scope.view<Real2D>("calc_Lscale_directly::rtm_pert_1", ngrdcol,
                                                   nzt); // For avg. calculation of Lscale  [kg/kg]
  Real2D rtm_pert_2 =
      workspace_scope.view<Real2D>("calc_Lscale_directly::rtm_pert_2", ngrdcol, nzt);
  Real2D thlm_pert_pos_rt = workspace_scope.view<Real2D>(
      "calc_Lscale_directly::thlm_pert_pos_rt", ngrdcol, nzt); // For avg. calculation of Lscale [K]
  Real2D thlm_pert_neg_rt =
      workspace_scope.view<Real2D>("calc_Lscale_directly::thlm_pert_neg_rt", ngrdcol, nzt);
  Real2D rtm_pert_pos_rt =
      workspace_scope.view<Real2D>("calc_Lscale_directly::rtm_pert_pos_rt", ngrdcol,
                                   nzt); // For avg. calculation of Lscale [kg/kg]
  Real2D rtm_pert_neg_rt =
      workspace_scope.view<Real2D>("calc_Lscale_directly::rtm_pert_neg_rt", ngrdcol, nzt);
  Real1D mu_pert_1      = workspace_scope.view<Real1D>("calc_Lscale_directly::mu_pert_1", ngrdcol);
  Real1D mu_pert_2      = workspace_scope.view<Real1D>("calc_Lscale_directly::mu_pert_2", ngrdcol);
  Real1D mu_pert_pos_rt = workspace_scope.view<Real1D>("calc_Lscale_directly::mu_pert_pos_rt",
                                                       ngrdcol); // For l_Lscale_plume_centered
  Real1D mu_pert_neg_rt =
      workspace_scope.view<Real1D>("calc_Lscale_directly::mu_pert_neg_rt", ngrdcol);
  // Capture component Views and zero-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto pdf_params_rt_1         = pdf_params.rt_1;
  const auto pdf_params_rt_2         = pdf_params.rt_2;
  const auto pdf_params_varnce_rt_1  = pdf_params.varnce_rt_1;
  const auto pdf_params_thl_1        = pdf_params.thl_1;
  const auto pdf_params_varnce_thl_1 = pdf_params.varnce_thl_1;
  const auto pdf_params_thl_2        = pdf_params.thl_2;
  const auto pdf_params_varnce_thl_2 = pdf_params.varnce_thl_2;
  const auto pdf_params_varnce_rt_2  = pdf_params.varnce_rt_2;

  if ((l_avg_Lscale && (Lscale_method != Lscale_direct))) {
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << err_info.err_header_global << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "l_avg_Lscale requires Lscale_method=1" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "Fatal error in calc_Lscale_directly" << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  if (clubb_at_least_debug_level_api(0)) {
    if ((l_Lscale_plume_centered && !l_avg_Lscale)) {
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << err_info.err_header_global << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "l_Lscale_plume_centered requires l_avg_Lscale" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Fatal error in calc_Lscale_directly" << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      // General error -> set all entries to clubb_fatal_error
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
  if ((l_avg_Lscale && !l_Lscale_plume_centered)) {
    // Call compute length two additional times with perturbed values
    // of rtm and thlm so that an average value of Lscale may be calculated.
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "calc_Lscale_directly::cells_0", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i           = (cell % ngrdcol);
            const int k           = (cell / ngrdcol);
            sign_rtpthlp_zt(i, k) = Kokkos::copysign(Kokkos::abs(one), rtpthlp_zt(i, k));
          });
    }
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "calc_Lscale_directly::cells_1", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            rtm_pert_1(i, k) =
                (rtm(i, k) +
                 (clubb_params(i, iLscale_pert_coef) *
                  Kokkos::sqrt(
                      (rtp2_zt(i, k) > (rt_tol * rt_tol) ? rtp2_zt(i, k) : (rt_tol * rt_tol)))));
          });
    }
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "calc_Lscale_directly::cells_2", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            thlm_pert_1(i, k) =
                (thlm(i, k) +
                 ((sign_rtpthlp_zt(i, k) * clubb_params(i, iLscale_pert_coef)) *
                  Kokkos::sqrt((thlp2_zt(i, k) > (thl_tol * thl_tol) ? thlp2_zt(i, k)
                                                                     : (thl_tol * thl_tol)))));
          });
    }
    Kokkos::parallel_for(
        "calc_Lscale_directly::columns_3", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i  = column;
          mu_pert_1(i) = (clubb_params(i, iLscale_mu_entrain) / clubb_params(i, iLscale_mu_coef));
        });
    compute_mixing_length(workspace, nzm, nzt, ngrdcol, gr, thvm, thlm_pert_1, // In
                          rtm_pert_1, em, Lscale_max, p_in_Pa,                 // In
                          exner, thv_ds_zt, mu_pert_1, lmin,                   // In
                          saturation_formula, l_implemented,                   // In
                          err_info,                                            // InOut
                          Lscale_pert_1, Lscale_up, Lscale_down                // Out
    );
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "calc_Lscale_directly::cells_4", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            rtm_pert_2(i, k) =
                (rtm(i, k) -
                 (clubb_params(i, iLscale_pert_coef) *
                  Kokkos::sqrt(
                      (rtp2_zt(i, k) > (rt_tol * rt_tol) ? rtp2_zt(i, k) : (rt_tol * rt_tol)))));
          });
    }
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "calc_Lscale_directly::cells_5", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            thlm_pert_2(i, k) =
                (thlm(i, k) -
                 ((sign_rtpthlp_zt(i, k) * clubb_params(i, iLscale_pert_coef)) *
                  Kokkos::sqrt((thlp2_zt(i, k) > (thl_tol * thl_tol) ? thlp2_zt(i, k)
                                                                     : (thl_tol * thl_tol)))));
          });
    }
    Kokkos::parallel_for(
        "calc_Lscale_directly::columns_6", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i  = column;
          mu_pert_2(i) = (clubb_params(i, iLscale_mu_entrain) * clubb_params(i, iLscale_mu_coef));
        });
    compute_mixing_length(workspace, nzm, nzt, ngrdcol, gr, thvm, thlm_pert_2, // In
                          rtm_pert_2, em, Lscale_max, p_in_Pa,                 // In
                          exner, thv_ds_zt, mu_pert_2, lmin,                   // In
                          saturation_formula, l_implemented,                   // In
                          err_info,                                            // InOut
                          Lscale_pert_2, Lscale_up, Lscale_down                // Out
    );
  } else if ((l_avg_Lscale && l_Lscale_plume_centered)) {
    // Take the values of thl and rt based one 1st or 2nd plume
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "calc_Lscale_directly::cells_7", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i           = (cell % ngrdcol);
            const int k           = (cell / ngrdcol);
            sign_rtpthlp_zt(i, k) = Kokkos::copysign(Kokkos::abs(one), rtpthlp_zt(i, k));
          });
    }
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "calc_Lscale_directly::cells_8", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            real_t Lscale_pert_coef{};
            Lscale_pert_coef = clubb_params(i, iLscale_pert_coef);
            if ((pdf_params_rt_1(i, k) > pdf_params_rt_2(i, k))) {
              rtm_pert_pos_rt(i, k) =
                  (pdf_params_rt_1(i, k) +
                   (Lscale_pert_coef *
                    Kokkos::sqrt((pdf_params_varnce_rt_1(i, k) > (rt_tol * rt_tol)
                                      ? pdf_params_varnce_rt_1(i, k)
                                      : (rt_tol * rt_tol)))));
              thlm_pert_pos_rt(i, k) =
                  (pdf_params_thl_1(i, k) +
                   ((sign_rtpthlp_zt(i, k) * Lscale_pert_coef) *
                    Kokkos::sqrt((pdf_params_varnce_thl_1(i, k) > (thl_tol * thl_tol)
                                      ? pdf_params_varnce_thl_1(i, k)
                                      : (thl_tol * thl_tol)))));
              thlm_pert_neg_rt(i, k) =
                  (pdf_params_thl_2(i, k) -
                   ((sign_rtpthlp_zt(i, k) * Lscale_pert_coef) *
                    Kokkos::sqrt((pdf_params_varnce_thl_2(i, k) > (thl_tol * thl_tol)
                                      ? pdf_params_varnce_thl_2(i, k)
                                      : (thl_tol * thl_tol)))));
              rtm_pert_neg_rt(i, k) =
                  (pdf_params_rt_2(i, k) -
                   (Lscale_pert_coef *
                    Kokkos::sqrt((pdf_params_varnce_rt_2(i, k) > (rt_tol * rt_tol)
                                      ? pdf_params_varnce_rt_2(i, k)
                                      : (rt_tol * rt_tol)))));
              // Lscale_weight = pdf_params%mixt_frac(i,k)
            } else {
              rtm_pert_pos_rt(i, k) =
                  (pdf_params_rt_2(i, k) +
                   (Lscale_pert_coef *
                    Kokkos::sqrt((pdf_params_varnce_rt_2(i, k) > (rt_tol * rt_tol)
                                      ? pdf_params_varnce_rt_2(i, k)
                                      : (rt_tol * rt_tol)))));
              thlm_pert_pos_rt(i, k) =
                  (pdf_params_thl_2(i, k) +
                   ((sign_rtpthlp_zt(i, k) * Lscale_pert_coef) *
                    Kokkos::sqrt((pdf_params_varnce_thl_2(i, k) > (thl_tol * thl_tol)
                                      ? pdf_params_varnce_thl_2(i, k)
                                      : (thl_tol * thl_tol)))));
              thlm_pert_neg_rt(i, k) =
                  (pdf_params_thl_1(i, k) -
                   ((sign_rtpthlp_zt(i, k) * Lscale_pert_coef) *
                    Kokkos::sqrt((pdf_params_varnce_thl_1(i, k) > (thl_tol * thl_tol)
                                      ? pdf_params_varnce_thl_1(i, k)
                                      : (thl_tol * thl_tol)))));
              rtm_pert_neg_rt(i, k) =
                  (pdf_params_rt_1(i, k) -
                   (Lscale_pert_coef *
                    Kokkos::sqrt((pdf_params_varnce_rt_1(i, k) > (rt_tol * rt_tol)
                                      ? pdf_params_varnce_rt_1(i, k)
                                      : (rt_tol * rt_tol)))));
              // Lscale_weight = 1.0_core_rknd - pdf_params%mixt_frac(i,k)
            }
          });
    }
    Kokkos::parallel_for(
        "calc_Lscale_directly::columns_9", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          mu_pert_pos_rt(i) =
              (clubb_params(i, iLscale_mu_entrain) / clubb_params(i, iLscale_mu_coef));
          mu_pert_neg_rt(i) =
              (clubb_params(i, iLscale_mu_entrain) * clubb_params(i, iLscale_mu_coef));
        });
    // Call length with perturbed values of thl and rt
    compute_mixing_length(workspace, nzm, nzt, ngrdcol, gr, thvm, thlm_pert_pos_rt, // In
                          rtm_pert_pos_rt, em, Lscale_max, p_in_Pa,                 // In
                          exner, thv_ds_zt, mu_pert_pos_rt, lmin,                   // In
                          saturation_formula, l_implemented,                        // In
                          err_info,                                                 // InOut
                          Lscale_pert_1, Lscale_up, Lscale_down                     // Out
    );
    compute_mixing_length(workspace, nzm, nzt, ngrdcol, gr, thvm, thlm_pert_neg_rt, // In
                          rtm_pert_neg_rt, em, Lscale_max, p_in_Pa,                 // In
                          exner, thv_ds_zt, mu_pert_neg_rt, lmin,                   // In
                          saturation_formula, l_implemented,                        // In
                          err_info,                                                 // InOut
                          Lscale_pert_2, Lscale_up, Lscale_down                     // Out
    );
  } else {
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "calc_Lscale_directly::cells_10", (level_count * ngrdcol), KOKKOS_LAMBDA(const int cell) {
            const int i         = (cell % ngrdcol);
            const int k         = (cell / ngrdcol);
            Lscale_pert_1(i, k) = unused_var;
            // Undefined
            Lscale_pert_2(i, k) = unused_var;
            // Undefined
          });
    }
  }
  // l_avg_Lscale
  if (stats.l_sample) {
    stats_update("Lscale_pert_1", Lscale_pert_1, // In
                 stats                           // InOut
    );
    stats_update("Lscale_pert_2", Lscale_pert_2, // In
                 stats                           // InOut
    );
  }
  // ********** NOTE: **********
  // This call to compute_mixing_length must be last.  Otherwise, the values
  // of
  // Lscale_up and Lscale_down in stats will be based on perturbation length
  // scales
  // rather than the mean length scale.
  // Diagnose CLUBB's turbulent mixing length scale.
  switch (Lscale_method) {
  case Lscale_coherence: {
    compute_coherent_mixing_length(
        workspace, nzm, nzt, ngrdcol, gr, thvm, thlm,                        // In
        rtm, pdf_params, em, Lscale_max, p_in_Pa,                            // In
        exner, thv_ds_zt,                                                    // In
        Kokkos::subview(clubb_params, Kokkos::ALL, iLscale_mu_entrain),      // In
        Kokkos::subview(clubb_params, Kokkos::ALL, iLscale_tke_multiplier),  // In
        Kokkos::subview(clubb_params, Kokkos::ALL, iLscale_pdf_launch_frac), // In
        Kokkos::subview(clubb_params, Kokkos::ALL, iLscale_up_weight),       // In
        lmin, l_implemented,                                                 // In
        err_info,                                                            // InOut
        Lscale, Lscale_up, Lscale_down                                       // Out
    );
    break;
  }
  case Lscale_direct: {
    compute_mixing_length(workspace, nzm, nzt, ngrdcol, gr, thvm, thlm, // In
                          rtm, em, Lscale_max, p_in_Pa,                 // In
                          exner, thv_ds_zt,                             // In
                          Kokkos::subview(clubb_params, Kokkos::ALL, iLscale_mu_entrain),
                          lmin,                              // In
                          saturation_formula, l_implemented, // In
                          err_info,                          // InOut
                          Lscale, Lscale_up, Lscale_down     // Out
    );
    break;
  }
  default: {
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << err_info.err_header_global << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    {
      std::ostringstream diagnostic;
      diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
      diagnostic << "calc_Lscale_directly cannot diagnose Lscale_method=" << ' ';
      diagnostic << Lscale_method << ' ';
      std::cerr << diagnostic.str() << '\n';
    }
    Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
    return;
  }
  }
  if (l_avg_Lscale) {
    if (l_Lscale_plume_centered) {
      // Weighted average of mean, pert_1, & pert_2
      // Lscale = 0.5_core_rknd * ( Lscale + Lscale_weight*Lscale_pert_1 &
      // + (1.0_core_rknd-Lscale_weight)*Lscale_pert_2
      // )
      // Weighted average of just the perturbed values
      // Lscale = Lscale_weight*Lscale_pert_1 +
      // (1.0_core_rknd-Lscale_weight)*Lscale_pert_2
      // Un-weighted average of just the perturbed values
      {
        const int level_count = nzt;
        Kokkos::parallel_for(
            "calc_Lscale_directly::cells_11", (level_count * ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i  = (cell % ngrdcol);
              const int k  = (cell / ngrdcol);
              Lscale(i, k) = (one_half * (Lscale_pert_1(i, k) + Lscale_pert_2(i, k)));
            });
      }
    } else {
      {
        const int level_count = nzt;
        Kokkos::parallel_for(
            "calc_Lscale_directly::cells_12", (level_count * ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i = (cell % ngrdcol);
              const int k = (cell / ngrdcol);
              Lscale(i, k) =
                  (one_third * ((Lscale(i, k) + Lscale_pert_1(i, k)) + Lscale_pert_2(i, k)));
            });
      }
    }
  }
  return;
}

void diagnose_Lscale_from_tau(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol,
                              const grid& gr,                                         // In
                              Real1D upwp_sfc, Real1D vpwp_sfc, Real2D ddzt_umvm_sqd, // In
                              Real2D ice_supersat_frac,                               // In
                              Real2D em,                                              // In
                              real_t ufmin, real_t tau_const,                         // In
                              Real1D sfc_elevation, Real1D Lscale_max,                // In
                              Real2D clubb_params,                                    // In
                              stats_type& stats,                                      // InOut
                              bool l_e3sm_config,                                     // In
                              bool l_smooth_Heaviside_tau_wpxp,                       // In
                              Real2D brunt_vaisala_freq_sqd_smth, Real2D Ri_zm,       // In
                              err_info_type& err_info,                                // InOut
                              Real2D invrs_tau_zt, Real2D invrs_tau_zm,               // Out
                              Real2D invrs_tau_sfc, Real2D invrs_tau_no_N2_zm,
                              Real2D invrs_tau_bkgnd,                           // Out
                              Real2D invrs_tau_shear, Real2D invrs_tau_N2_iso,  // Out
                              Real2D invrs_tau_wp2_zm, Real2D invrs_tau_xp2_zm, // Out
                              Real2D invrs_tau_wp3_zm, Real2D invrs_tau_wp3_zt,
                              Real2D invrs_tau_wpxp_zm, // Out
                              Real2D tau_max_zm, Real2D tau_max_zt, Real2D tau_zm,
                              Real2D tau_zt,                                      // Out
                              Real2D Lscale, Real2D Lscale_up, Real2D Lscale_down // Out
) {
  // Description:
  // Diagnose inverse damping times and mixing length from surface stress,
  // momentum-level TKE, shear, and stability profiles (method 2). Return the
  // separate inverse-time contributions [1/s], bounded time scales [s], and
  // thermodynamic-level lengths [m]; update statistics and error information.
  // References:
  // Guo et al.(2021, JAMES)
  // -----------------------------------------------------------------------------------------------

  // Arguments (Native order)
  // nzm: Number of momentum levels
  // nzt: Number of thermodynamic levels
  // ngrdcol: Number of grid columns
  // gr: Grid structure
  // upwp_sfc: Surface value of u'w' [m^2/s^2]
  // vpwp_sfc: Surface value of v'w' [m^2/s^2]
  // ddzt_umvm_sqd: Squared vertical shear of horizontal wind [s^-2]
  // ice_supersat_frac: Ice supersaturation fraction [-]
  // em: Turbulent kinetic energy [m^2/s^2]
  // ufmin: Minimum velocity scale [m/s]
  // tau_const: Constant eddy dissipation time scale [s]
  // sfc_elevation: Surface elevation [m]
  // Lscale_max: Maximum allowable mixing length [m]
  // clubb_params: Array of CLUBB's tunable parameters    [units vary]
  // stats: Statistics accumulator
  // l_e3sm_config: Whether E3SM configuration options are active
  // l_smooth_Heaviside_tau_wpxp: Use the smoothed Heaviside 'Peskin' function to compute
  // invrs_tau_wpxp_zm
  // brunt_vaisala_freq_sqd_smth: smoothed Buoyancy frequency squared, N^2     [s^-2]
  // Ri_zm: Richardson number
  // err_info: err_info struct containing err_code and err_header
  // invrs_tau_zt: Inverse time scale on thermodynamic levels [1/s]
  // invrs_tau_zm: Inverse time scale on momentum levels [1/s]
  // invrs_tau_sfc: Surface contribution to inverse tau [1/s]
  // invrs_tau_no_N2_zm: Inverse tau without N2 contribution on momentum levels [1/s]
  // invrs_tau_bkgnd: Background inverse tau [1/s]
  // invrs_tau_shear: Shear contribution to inverse tau [1/s]
  // invrs_tau_N2_iso: Isotropic N2 contribution to inverse tau [1/s]
  // invrs_tau_wp2_zm: Inverse time scale for wp2 terms on momentum levels [1/s]
  // invrs_tau_xp2_zm: Inverse time scale for xp2 terms on momentum levels [1/s]
  // invrs_tau_wp3_zm: Inverse time scale for wp3 terms on momentum levels [1/s]
  // invrs_tau_wp3_zt: Inverse time scale for wp3 terms on thermodynamic levels [1/s]
  // invrs_tau_wpxp_zm: Inverse time scale for wpxp terms on momentum levels [1/s]
  // tau_max_zm: Maximum allowable eddy dissipation time scale on momentum levels [s]
  // tau_max_zt: Maximum allowable eddy dissipation time scale on thermodynamic levels [s]
  // tau_zm: Eddy dissipation time scale on momentum levels [s]
  // tau_zt: Eddy dissipation time scale on thermodynamic levels [s]
  // Lscale: Turbulent mixing length on thermodynamic levels [m]
  // Lscale_up: Upward component of turbulent mixing length [m]
  // Lscale_down: Downward component of turbulent mixing length [m]

  auto workspace_scope  = workspace.scope();
  Real2D brunt_freq_pos = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "diagnose_Lscale_from_tau::brunt_freq_pos", ngrdcol, nzm);
  Real2D brunt_freq_out_cloud = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "diagnose_Lscale_from_tau::brunt_freq_out_cloud", ngrdcol, nzm);
  Real2D bvf_thresh =
      (l_smooth_Heaviside_tau_wpxp)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "diagnose_Lscale_from_tau::bvf_thresh", ngrdcol, nzm)
          : Real2D{}; // temporatory array
  Real2D H_invrs_tau_wpxp_N2 = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "diagnose_Lscale_from_tau::H_invrs_tau_wpxp_N2", ngrdcol,
      nzm); // Heaviside function for clippings of invrs_tau_wpxp_N2
  Real1D ustar = workspace_scope.view<Real1D>(Kokkos::WithoutInitializing,
                                              "diagnose_Lscale_from_tau::ustar", ngrdcol);
  constexpr real_t heaviside_smth_range =
      real_t(1.0e-0); // range where Heaviside function is smoothed
  constexpr bool l_smooth_min_max = false;
  Real2D norm_ddzt_umvm           = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "diagnose_Lscale_from_tau::norm_ddzt_umvm", ngrdcol, nzm);
  Real2D smooth_norm_ddzt_umvm = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "diagnose_Lscale_from_tau::smooth_norm_ddzt_umvm", ngrdcol, nzm);
  Real2D brunt_vaisala_freq_clipped =
      (l_smooth_min_max)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "diagnose_Lscale_from_tau::brunt_vaisala_freq_clipped",
                                         ngrdcol, nzm)
          : Real2D{};
  Real2D ice_supersat_frac_zm = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "diagnose_Lscale_from_tau::ice_supersat_frac_zm", ngrdcol, nzm);
  Real2D invrs_tau_shear_smooth = workspace_scope.view<Real2D>(
      Kokkos::WithoutInitializing, "diagnose_Lscale_from_tau::invrs_tau_shear_smooth", ngrdcol,
      nzm);
  Real2D Ri_zm_smooth =
      (l_smooth_min_max)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "diagnose_Lscale_from_tau::Ri_zm_smooth", ngrdcol, nzm)
          : Real2D{};
  Real2D em_clipped =
      (l_smooth_min_max)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "diagnose_Lscale_from_tau::em_clipped", ngrdcol, nzm)
          : Real2D{};
  Real2D tau_zm_unclipped =
      (l_smooth_min_max)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "diagnose_Lscale_from_tau::tau_zm_unclipped", ngrdcol, nzm)
          : Real2D{};
  Real2D tmp_calc =
      (l_smooth_min_max)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "diagnose_Lscale_from_tau::tmp_calc", ngrdcol, nzm)
          : Real2D{};
  Real2D tmp_calc_max =
      (l_smooth_min_max)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "diagnose_Lscale_from_tau::tmp_calc_max", ngrdcol, nzm)
          : Real2D{};
  Real2D tmp_calc_min_max =
      (l_smooth_min_max)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "diagnose_Lscale_from_tau::tmp_calc_min_max", ngrdcol, nzm)
          : Real2D{};
  Real2D em_zt = workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                              "diagnose_Lscale_from_tau::em_zt", ngrdcol, nzt);
  Real2D tau_zt_unclipped =
      (l_smooth_min_max)
          ? workspace_scope.view<Real2D>(Kokkos::WithoutInitializing,
                                         "diagnose_Lscale_from_tau::tau_zt_unclipped", ngrdcol, nzt)
          : Real2D{};
  Real1D tmp_calc_ngrdcol = workspace_scope.view<Real1D>(
      Kokkos::WithoutInitializing, "diagnose_Lscale_from_tau::tmp_calc_ngrdcol", ngrdcol);
  // Capture component Views and zero-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_zm = gr.zm;

  zm2zt_api(nzm, nzt, ngrdcol, gr, em, em_zt, em_min);
  const auto zm_error     = gr.zm;
  const int k_lb_zm_error = (gr.k_lb_zm);
  // Native per-column fatal detection runs on device; print only failures on host.
  Kokkos::View<int*, Kokkos::LayoutLeft> bad_columns =
      workspace_scope.view<Kokkos::View<int*, Kokkos::LayoutLeft>>(
          Kokkos::WithoutInitializing, "diagnose_Lscale_from_tau::bad_columns", ngrdcol);
  const auto err_code_check = err_info.err_code;
  int bad_count             = 0;
  Kokkos::parallel_reduce(
      "diagnose_Lscale_from_tau::input_check", ngrdcol,
      KOKKOS_LAMBDA(const int column, int& count) {
        bad_columns(column) = 0;
        if (zm_error(column, k_lb_zm_error) - sfc_elevation(column) +
                clubb_params(column, iz_displace) <
            eps) {
          bad_columns(column)    = 1;
          err_code_check(column) = clubb_fatal_error;
          count += 1;
        }
      },
      bad_count);
  if (bad_count != 0) {
    auto diagnostic_scope  = workspace.scope();
    const auto bad_host    = to_host_space(workspace, bad_columns);
    const auto header_host = to_host_space(workspace, err_info.err_header);
    for (int column = 0; column < ngrdcol; ++column) {
      if (bad_host(column) != 0) {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic.write(header_host.data() + (column * header_host.extent(0)),
                         header_host.extent(0));
        diagnostic << '\n';
        std::cerr << diagnostic.str();
      }
    }
  }
  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      {
        std::ostringstream diagnostic;
        diagnostic << std::setprecision(std::numeric_limits<real_t>::max_digits10);
        diagnostic << "Lowest zm grid level is below ground in CLUBB." << ' ';
        std::cerr << diagnostic.str() << '\n';
      }
      return;
    }
  }
  Kokkos::parallel_for(
      "diagnose_Lscale_from_tau::columns_0", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column;
        tmp_calc_ngrdcol(i) =
            Kokkos::pow(((upwp_sfc(i) * upwp_sfc(i)) + (vpwp_sfc(i) * vpwp_sfc(i))), one_fourth);
      });
  if (l_smooth_min_max) {
    // tmp_calc_ngrdcol used here because openacc has an issue with
    // a variable being both input and output of a function
    smooth_max(ngrdcol, tmp_calc_ngrdcol, ufmin, min_max_smth_mag, ustar);
  } else {
    Kokkos::parallel_for(
        "diagnose_Lscale_from_tau::columns_1", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          ustar(i)    = (tmp_calc_ngrdcol(i) > ufmin ? tmp_calc_ngrdcol(i) : ufmin);
        });
  }
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "diagnose_Lscale_from_tau::cells_2", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i           = (cell % ngrdcol);
          const int k           = (cell / ngrdcol);
          invrs_tau_bkgnd(i, k) = (clubb_params(i, iC_invrs_tau_bkgnd) / tau_const);
          norm_ddzt_umvm(i, k)  = Kokkos::sqrt(ddzt_umvm_sqd(i, k));
        });
  }

  zm2zt2zm(workspace, nzm, nzt, ngrdcol, gr, norm_ddzt_umvm, smooth_norm_ddzt_umvm);
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "diagnose_Lscale_from_tau::cells_4", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          invrs_tau_shear_smooth(i, k) =
              (clubb_params(i, iC_invrs_tau_shear) * smooth_norm_ddzt_umvm(i, k));
        });
  }
  // Enforce that invrs_tau_shear is positive
  smooth_max(nzm, ngrdcol, invrs_tau_shear_smooth, zero_threshold, min_max_smth_mag,
             invrs_tau_shear);
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "diagnose_Lscale_from_tau::cells_5", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i         = (cell % ngrdcol);
          const int k         = (cell / ngrdcol);
          invrs_tau_sfc(i, k) = ((clubb_params(i, iC_invrs_tau_sfc) * (ustar(i) / vonk)) /
                                 ((gr_zm(i, k) - sfc_elevation(i)) + clubb_params(i, iz_displace)));
          // C_invrs_tau_sfc * ( wp2 / vonk /ustar ) / ( gr%zm(1,:) -sfc_elevation + z_displace )
        });
  }
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "diagnose_Lscale_from_tau::cells_6", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          invrs_tau_no_N2_zm(i, k) =
              ((invrs_tau_bkgnd(i, k) + invrs_tau_sfc(i, k)) + invrs_tau_shear(i, k));
        });
  }
  if (l_smooth_min_max) {
    smooth_max(nzm, ngrdcol, zero_threshold, brunt_vaisala_freq_sqd_smth,
               (real_t(1.0e-4) * min_max_smth_mag), brunt_vaisala_freq_clipped);
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_7", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i          = (cell % ngrdcol);
            const int k          = (cell / ngrdcol);
            brunt_freq_pos(i, k) = Kokkos::sqrt(brunt_vaisala_freq_clipped(i, k));
          });
    }
  } else {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_8", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i          = (cell % ngrdcol);
            const int k          = (cell / ngrdcol);
            brunt_freq_pos(i, k) = Kokkos::sqrt((zero_threshold > brunt_vaisala_freq_sqd_smth(i, k)
                                                     ? zero_threshold
                                                     : brunt_vaisala_freq_sqd_smth(i, k)));
          });
    }
  }
  zt2zm_api(nzm, nzt, ngrdcol, gr, ice_supersat_frac, ice_supersat_frac_zm, zero_threshold);
  if (l_smooth_min_max) {
    // roll this back as well once checks have passed
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_9", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i    = (cell % ngrdcol);
            const int k    = (cell / ngrdcol);
            tmp_calc(i, k) = (one - (ice_supersat_frac_zm(i, k) / real_t(0.001)));
          });
    }
    smooth_max(nzm, ngrdcol, zero_threshold, tmp_calc, min_max_smth_mag, tmp_calc_max);
    smooth_min(nzm, ngrdcol, one, tmp_calc_max, min_max_smth_mag, tmp_calc_min_max);
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_10", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i                = (cell % ngrdcol);
            const int k                = (cell / ngrdcol);
            brunt_freq_out_cloud(i, k) = (brunt_freq_pos(i, k) * tmp_calc_min_max(i, k));
          });
    }
  } else {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_11", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            brunt_freq_out_cloud(i, k) =
                (brunt_freq_pos(i, k) *
                 (one < (zero_threshold > (one - (ice_supersat_frac_zm(i, k) / real_t(0.001)))
                             ? zero_threshold
                             : (one - (ice_supersat_frac_zm(i, k) / real_t(0.001))))
                      ? one
                      : (zero_threshold > (one - (ice_supersat_frac_zm(i, k) / real_t(0.001)))
                             ? zero_threshold
                             : (one - (ice_supersat_frac_zm(i, k) / real_t(0.001))))));
          });
    }
  }
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "diagnose_Lscale_from_tau::cells_12", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          if ((gr_zm(i, k) < clubb_params(i, ialtitude_threshold))) {
            brunt_freq_out_cloud(i, k) = zero;
          }
        });
  }
  // Write both bv extra terms to invrs_taus to disk
  if (stats.l_sample) {
    stats_update("bv_freq_pos", brunt_freq_pos, // In
                 stats                          // InOut
    );
    stats_update("bv_freq_out_cloud", brunt_freq_out_cloud, // In
                 stats                                      // InOut
    );
  }
  // This time scale is used optionally for the return-to-isotropy term. It
  // omits invrs_tau_sfc based on the rationale that the isotropization
  // rate shouldn't be enhanced near the ground.
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "diagnose_Lscale_from_tau::cells_13", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          real_t C_invrs_tau_N2{};
          real_t C_invrs_tau_N2_wp2{};
          C_invrs_tau_N2         = clubb_params(i, iC_invrs_tau_N2);
          C_invrs_tau_N2_wp2     = clubb_params(i, iC_invrs_tau_N2_wp2);
          invrs_tau_N2_iso(i, k) = ((invrs_tau_bkgnd(i, k) + invrs_tau_shear(i, k)) +
                                    (C_invrs_tau_N2_wp2 * brunt_freq_pos(i, k)));
          invrs_tau_wp2_zm(i, k) =
              ((invrs_tau_no_N2_zm(i, k) + (C_invrs_tau_N2 * brunt_freq_pos(i, k))) +
               (C_invrs_tau_N2_wp2 * brunt_freq_out_cloud(i, k)));
          invrs_tau_zm(i, k) = (invrs_tau_no_N2_zm(i, k) + (C_invrs_tau_N2 * brunt_freq_pos(i, k)));
        });
  }
  if (l_e3sm_config) {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_14", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i        = (cell % ngrdcol);
            const int k        = (cell / ngrdcol);
            invrs_tau_zm(i, k) = (one_half * invrs_tau_zm(i, k));
          });
    }
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_15", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            invrs_tau_xp2_zm(i, k) =
                ((invrs_tau_no_N2_zm(i, k) +
                  (clubb_params(i, iC_invrs_tau_N2_xp2) * brunt_freq_pos(i, k))) + // 0
                 (((clubb_params(i, iC_invrs_tau_sfc) * two) * Kokkos::sqrt(em(i, k))) /
                  ((gr_zm(i, k) - sfc_elevation(i)) + clubb_params(i, iz_displace)))); // small
            // small
            // 0
          });
    }
    if (l_smooth_min_max) {
      smooth_max(nzm, ngrdcol, real_t(1.0e-7), brunt_vaisala_freq_sqd_smth,
                 (real_t(1.0e-4) * min_max_smth_mag), brunt_vaisala_freq_clipped);
      {
        const int level_count = nzm;
        Kokkos::parallel_for(
            "diagnose_Lscale_from_tau::cells_16", (level_count * ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i = (cell % ngrdcol);
              const int k = (cell / ngrdcol);
              tmp_calc(i, k) =
                  Kokkos::sqrt((ddzt_umvm_sqd(i, k) / brunt_vaisala_freq_clipped(i, k)));
            });
      }
      smooth_max(nzm, ngrdcol, tmp_calc, real_t(0.3), (real_t(0.3) * min_max_smth_mag),
                 tmp_calc_max);
      smooth_min(nzm, ngrdcol, tmp_calc_max, one, min_max_smth_mag, tmp_calc_min_max);
      {
        const int level_count = nzm;
        Kokkos::parallel_for(
            "diagnose_Lscale_from_tau::cells_17", (level_count * ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i            = (cell % ngrdcol);
              const int k            = (cell / ngrdcol);
              invrs_tau_xp2_zm(i, k) = (tmp_calc_min_max(i, k) * invrs_tau_xp2_zm(i, k));
            });
      }
    } else {
      {
        const int level_count = nzm;
        Kokkos::parallel_for(
            "diagnose_Lscale_from_tau::cells_18", (level_count * ngrdcol),
            KOKKOS_LAMBDA(const int cell) {
              const int i = (cell % ngrdcol);
              const int k = (cell / ngrdcol);
              invrs_tau_xp2_zm(i, k) =
                  (((Kokkos::sqrt((ddzt_umvm_sqd(i, k) /
                                   (real_t(1.0e-7) > brunt_vaisala_freq_sqd_smth(i, k)
                                        ? real_t(1.0e-7)
                                        : brunt_vaisala_freq_sqd_smth(i, k)))) > real_t(0.3)
                         ? Kokkos::sqrt((ddzt_umvm_sqd(i, k) /
                                         (real_t(1.0e-7) > brunt_vaisala_freq_sqd_smth(i, k)
                                              ? real_t(1.0e-7)
                                              : brunt_vaisala_freq_sqd_smth(i, k))))
                         : real_t(0.3)) < one
                        ? (Kokkos::sqrt((ddzt_umvm_sqd(i, k) /
                                         (real_t(1.0e-7) > brunt_vaisala_freq_sqd_smth(i, k)
                                              ? real_t(1.0e-7)
                                              : brunt_vaisala_freq_sqd_smth(i, k)))) > real_t(0.3)
                               ? Kokkos::sqrt((ddzt_umvm_sqd(i, k) /
                                               (real_t(1.0e-7) > brunt_vaisala_freq_sqd_smth(i, k)
                                                    ? real_t(1.0e-7)
                                                    : brunt_vaisala_freq_sqd_smth(i, k))))
                               : real_t(0.3))
                        : one) *
                   invrs_tau_xp2_zm(i, k));
            });
      }
    }
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_19", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            invrs_tau_wpxp_zm(i, k) =
                ((two * invrs_tau_zm(i, k)) +
                 (clubb_params(i, iC_invrs_tau_N2_wpxp) * brunt_freq_out_cloud(i, k)));
          });
    }
  } else {
    // l_e3sm_config = false
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_20", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            invrs_tau_xp2_zm(i, k) =
                ((invrs_tau_no_N2_zm(i, k) +
                  (clubb_params(i, iC_invrs_tau_N2) * brunt_freq_pos(i, k))) +
                 (clubb_params(i, iC_invrs_tau_N2_xp2) * brunt_freq_out_cloud(i, k)));
          });
    }
    // The earlier interpolation remains valid: its source and grid are unchanged.
    // !$acc parallel loop gang vector collapse(2) default(present)
    // do k = 1, nzm
    // do i = 1, ngrdcol
    // if ( ice_supersat_frac_zm(i,k) <= 0.01_core_rknd &
    // .and. invrs_tau_xp2_zm(i,k)  >= 0.003_core_rknd ) then
    // invrs_tau_xp2_zm(i,k) = 0.003_core_rknd
    // end if
    // end do
    // end do
    // !$acc end parallel loop
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_21", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            invrs_tau_wpxp_zm(i, k) =
                ((invrs_tau_no_N2_zm(i, k) +
                  (clubb_params(i, iC_invrs_tau_N2) * brunt_freq_pos(i, k))) +
                 (clubb_params(i, iC_invrs_tau_N2_wpxp) * brunt_freq_out_cloud(i, k)));
          });
    }
  }
  // l_e3sm_config
  if (l_smooth_Heaviside_tau_wpxp) {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_22", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i      = (cell % ngrdcol);
            const int k      = (cell / ngrdcol);
            bvf_thresh(i, k) = ((brunt_vaisala_freq_sqd_smth(i, k) /
                                 clubb_params(i, iC_invrs_tau_wpxp_N2_thresh)) -
                                one);
          });
    }
    smooth_heaviside_peskin(nzm, ngrdcol, bvf_thresh, heaviside_smth_range, H_invrs_tau_wpxp_N2);
  } else {
    // l_smooth_Heaviside_tau_wpxp = .false.
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_23", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            if ((brunt_vaisala_freq_sqd_smth(i, k) >
                 clubb_params(i, iC_invrs_tau_wpxp_N2_thresh))) {
              H_invrs_tau_wpxp_N2(i, k) = one;
            } else {
              H_invrs_tau_wpxp_N2(i, k) = zero;
            }
          });
    }
  }
  // l_smooth_Heaviside_tau_wpxp
  if (l_smooth_min_max) {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_24", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i    = (cell % ngrdcol);
            const int k    = (cell / ngrdcol);
            tmp_calc(i, k) = (clubb_params(i, iC_invrs_tau_wpxp_Ri) *
                              Kokkos::pow(Ri_zm(i, k), clubb_params(i, iwpxp_Ri_exp)));
          });
    }
    smooth_min(nzm, ngrdcol, tmp_calc, real_t(12.0), (real_t(12.0) * min_max_smth_mag),
               Ri_zm_smooth);
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_25", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            if ((gr_zm(i, k) > clubb_params(i, ialtitude_threshold))) {
              invrs_tau_wpxp_zm(i, k) = (invrs_tau_wpxp_zm(i, k) *
                                         (one + (H_invrs_tau_wpxp_N2(i, k) * Ri_zm_smooth(i, k))));
            }
          });
    }
  } else {
    // l_smooth_min_max
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_26", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            if ((gr_zm(i, k) > clubb_params(i, ialtitude_threshold))) {
              invrs_tau_wpxp_zm(i, k) =
                  (invrs_tau_wpxp_zm(i, k) *
                   (one + (H_invrs_tau_wpxp_N2(i, k) *
                           ((clubb_params(i, iC_invrs_tau_wpxp_Ri) *
                             Kokkos::pow((Ri_zm(i, k) > zero ? Ri_zm(i, k) : zero),
                                         clubb_params(i, iwpxp_Ri_exp))) < real_t(12.0)
                                ? (clubb_params(i, iC_invrs_tau_wpxp_Ri) *
                                   Kokkos::pow((Ri_zm(i, k) > zero ? Ri_zm(i, k) : zero),
                                               clubb_params(i, iwpxp_Ri_exp)))
                                : real_t(12.0)))));
            }
          });
    }
  }
  {
    const int level_count = nzm;
    Kokkos::parallel_for(
        "diagnose_Lscale_from_tau::cells_27", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          invrs_tau_wp3_zm(i, k) =
              (invrs_tau_wp2_zm(i, k) +
               (clubb_params(i, iC_invrs_tau_N2_clear_wp3) * brunt_freq_out_cloud(i, k)));
        });
  }
  // Calculate the maximum allowable value of time-scale tau,
  // which depends of the value of Lscale_max.
  if (l_smooth_min_max) {
    smooth_max(nzm, ngrdcol, em, em_min, min_max_smth_mag, em_clipped);
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_28", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i      = (cell % ngrdcol);
            const int k      = (cell / ngrdcol);
            tau_max_zt(i, k) = (Lscale_max(i) / Kokkos::sqrt(em_zt(i, k)));
          });
    }
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_29", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i      = (cell % ngrdcol);
            const int k      = (cell / ngrdcol);
            tau_max_zm(i, k) = (Lscale_max(i) / Kokkos::sqrt(em_clipped(i, k)));
          });
    }
  } else {
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_30", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i      = (cell % ngrdcol);
            const int k      = (cell / ngrdcol);
            tau_max_zt(i, k) = (Lscale_max(i) / Kokkos::sqrt(em_zt(i, k)));
          });
    }
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_31", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            tau_max_zm(i, k) =
                (Lscale_max(i) / Kokkos::sqrt((em(i, k) > em_min ? em(i, k) : em_min)));
          });
    }
  }
  if (l_smooth_min_max) {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_32", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i            = (cell % ngrdcol);
            const int k            = (cell / ngrdcol);
            tau_zm_unclipped(i, k) = (one / invrs_tau_zm(i, k));
          });
    }
    smooth_min(nzm, ngrdcol, tau_zm_unclipped, tau_max_zm, (real_t(1.0e3) * min_max_smth_mag),
               tau_zm);
    zm2zt_api(nzm, nzt, ngrdcol, gr, tau_zm, tau_zt_unclipped);
    smooth_min(nzt, ngrdcol, tau_zt_unclipped, tau_max_zt, (real_t(1.0e3) * min_max_smth_mag),
               tau_zt);
  } else {
    {
      const int level_count = nzm;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_33", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i = (cell % ngrdcol);
            const int k = (cell / ngrdcol);
            tau_zm(i, k) =
                ((one / invrs_tau_zm(i, k)) < tau_max_zm(i, k) ? (one / invrs_tau_zm(i, k))
                                                               : tau_max_zm(i, k));
          });
    }
    zm2zt_api(nzm, nzt, ngrdcol, gr, tau_zm, tau_zt);
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "diagnose_Lscale_from_tau::cells_34", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i  = (cell % ngrdcol);
            const int k  = (cell / ngrdcol);
            tau_zt(i, k) = (tau_zt(i, k) < tau_max_zt(i, k) ? tau_zt(i, k) : tau_max_zt(i, k));
          });
    }
  }
  zm2zt_api(nzm, nzt, ngrdcol, gr, invrs_tau_zm, invrs_tau_zt);
  zm2zt_api(nzm, nzt, ngrdcol, gr, invrs_tau_wp3_zm, invrs_tau_wp3_zt);
  {
    const int level_count = nzt;
    Kokkos::parallel_for(
        "diagnose_Lscale_from_tau::cells_35", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i  = (cell % ngrdcol);
          const int k  = (cell / ngrdcol);
          Lscale(i, k) = (tau_zt(i, k) * Kokkos::sqrt(em_zt(i, k)));
          // Lscale_up and Lscale_down aren't calculated with this option.
          // They are set to 0 for stats output.
          Lscale_up(i, k)   = zero;
          Lscale_down(i, k) = zero;
        });
  }
  return;
}

namespace {

void compute_coherent_mixing_length(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol,
                                    const grid& gr, Real2D thvm,
                                    Real2D thlm, // In
                                    Real2D rtm, const pdf_parameter& pdf_params, Real2D em,
                                    Real1D Lscale_max,                                     // In
                                    Real2D p_in_Pa,                                        // In
                                    Real2D exner, Real2D thv_ds, Real1D Lscale_mu_entrain, // In
                                    Real1D Lscale_tke_multiplier,                          // In
                                    Real1D Lscale_pdf_launch_frac,                         // In
                                    Real1D Lscale_up_weight, Real1D lmin,                  // In
                                    bool l_implemented,                                    // In
                                    err_info_type& err_info,                               // InOut
                                    Real2D Lscale, Real2D Lscale_up, Real2D Lscale_down    // Out
) {
  // Description:
  // Diagnose mixing length by weighting entraining moist-parcel travel with
  // the fraction of assumed launch energies able to cross each energy barrier.
  // Read thermodynamic-level environmental and PDF properties, momentum-level
  // TKE, positive entrainment rates [1/m], and per-column launch/mixing weights.
  // Return combined, upward, and downward lengths on thermodynamic levels [m];
  // err_info is passed through unchanged.
  // Analytic layer updates advance each trajectory in path order. With launch
  // energy scale E and running adverse-work barrier H [m^2/s^2], endpoint
  // coherence is (E/(E+H))^2. Integrate this weight assuming H is linear across
  // the layer, include the first layer at mean coherence <= 0.01, then stop.
  // A source constant selects local, maximum-endpoint, or arrival-weighted
  // directional lengths before floors, geometric mixing, and the final cap.
  // Parcel saturation always uses the local Flatau polynomial.
  // References:
  // doc/coherence_mixing_length.md: survival law and aggregation equations.
  // Golaz et al. (2002), J. Atmos. Sci., 59, 3540--3551: parcel background.
  // -----------------------------------------------------------------------

  // Arguments (Native order)
  // nzm: Number of momentum levels
  // nzt: Number of thermodynamic levels
  // ngrdcol: Number of grid columns
  // gr: Grid structure
  // thvm: Environmental virtual potential temperature [K]
  // thlm: Environmental liquid-water potential temperature [K]
  // rtm: Environmental total-water mixing ratio [kg/kg]
  // pdf_params: PDF parameters used for directional launch properties
  // em: Turbulent kinetic energy on momentum levels [m^2/s^2]
  // Lscale_max: Maximum allowable mixing length [m]
  // p_in_Pa: Pressure [Pa]
  // exner: Exner function [-]
  // thv_ds: Reference-state virtual potential temperature [K]
  // Lscale_mu_entrain: Fractional entrainment rate [1/m]
  // Lscale_tke_multiplier: Mean launch energy divided by TKE [-]
  // Lscale_pdf_launch_frac: PDF departure retained in launch properties [-]
  // Lscale_up_weight: Upward weight in directional geometric mixing [-]
  // lmin: Surface minimum mixing length [m]
  // l_implemented: True when CLUBB is run within a host model
  // err_info: Error information
  // Lscale: Turbulent mixing length [m]
  // Lscale_up: Upward component of turbulent mixing length [m]
  // Lscale_down: Downward component of turbulent mixing length [m]

  auto workspace_scope                      = workspace.scope();
  constexpr int no_overwrite                = 0; // Retain each locally launched reach
  constexpr int legacy_max_overwrite        = 1; // Extend the most distant parcel endpoint
  constexpr int coherence_average_overwrite = 2; // Average remaining reach by arrival coherence
  constexpr int overwrite_mode         = coherence_average_overwrite; // Selected overwrite method
  constexpr real_t zlmin               = real_t(0.1);  // Directional floor after aggregation [m]
  constexpr real_t surface_layer_depth = real_t(500.); // Depth of surface floor [m]
  constexpr real_t invs_sfclyr_depth =
      (one / surface_layer_depth); // Inverse surface-floor depth [1/m]
  Real2D coherence_energy =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::coherence_energy", ngrdcol,
                                   nzt); // Coherent parcel energy scale [m^2/s^2]
  Real2D crossed_length_up =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::crossed_length_up", ngrdcol,
                                   nzt); // Upward coherence integral crossing each layer [m]
  Real2D crossed_length_down =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::crossed_length_down", ngrdcol,
                                   nzt); // Downward coherence integral crossing each layer [m]
  Real2D arrival_coherence_up =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::arrival_coherence_up", ngrdcol,
                                   nzt); // Upward parcel coherence at each target [-]
  Real2D arrival_coherence_down =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::arrival_coherence_down",
                                   ngrdcol, nzt); // Downward parcel coherence at each target [-]
  Real2D layer_retention =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::layer_retention", ngrdcol,
                                   nzt); // Incoming parcel-property retention factor [-]
  Real2D thl_entrainment =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::thl_entrainment", ngrdcol,
                                   nzt); // Environmental contribution to parcel thl [K]
  Real2D rt_entrainment =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::rt_entrainment", ngrdcol,
                                   nzt); // Environmental contribution to parcel rt [kg/kg]
  Real2D thl_launch =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::thl_launch", ngrdcol,
                                   nzt); // Directional parcel launch thl [K]
  Real2D rt_launch =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::rt_launch", ngrdcol,
                                   nzt); // Directional parcel launch rt [kg/kg]
  Real2D launch_buoyancy =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::launch_buoyancy", ngrdcol,
                                   nzt); // Directional parcel launch buoyancy [m/s^2]
  Real2D grav_coef =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::grav_coef", ngrdcol,
                                   nzt); // Gravity divided by environmental thv [m/s^2/K]
  Real2D vapor_coef =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::vapor_coef", ngrdcol,
                                   nzt); // Water-vapor contribution to parcel thv [K]
  Real2D latent_coef =
      workspace_scope.view<Real2D>("compute_coherent_mixing_length::latent_coef", ngrdcol,
                                   nzt); // Liquid-water contribution to parcel thv [K]
  real_t liquid_coef{};                  // Constant in the saturation-excess formula [K^2]
  // Capture component Views and zero-based grid metadata; each kernel
  // uses only the aliases it references, keeping owning UDTs on the host.
  const auto gr_k_lb_zt       = (gr.k_lb_zt);
  const auto gr_grid_dir_indx = gr.grid_dir_indx;
  const auto gr_k_ub_zt       = (gr.k_ub_zt);
  const auto gr_zt            = gr.zt;
  const auto pdf_params_thl_1 = pdf_params.thl_1;
  const auto pdf_params_rt_1  = pdf_params.rt_1;
  const auto pdf_params_thl_2 = pdf_params.thl_2;
  const auto pdf_params_rt_2  = pdf_params.rt_2;
  const auto gr_zm            = gr.zm;
  const auto gr_k_lb_zm       = (gr.k_lb_zm);

  // Interpolate launch TKE to thermodynamic levels.
  zm2zt_api(nzm, nzt, ngrdcol, gr, em, coherence_energy);
  liquid_coef = ((ep * (Lv * Lv)) / (Rd * Cp));
  // Precompute the inverse-square energy scale and buoyancy coefficients.
  {
    const int level_count = nzt;
    Kokkos::parallel_for(
        "compute_coherent_mixing_length::cells_0", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          // This is the value we use for the mean of the assumed parcel energy distribution.
          // It is a somewhat arbitrary choice, so it was given a scaling parameter
          coherence_energy(i, k) = (coherence_energy(i, k) * Lscale_tke_multiplier(i));
          // These coefficients reconstruct parcel thv and scale its buoyancy.
          grav_coef(i, k)   = (grav / thvm(i, k));
          vapor_coef(i, k)  = (ep1 * thv_ds(i, k));
          latent_coef(i, k) = ((Lv / (exner(i, k) * Cp)) - (ep2 * thv_ds(i, k)));
        });
  }
  // The coherence-average mode needs directional crossing sums. At each target,
  // its own local launch contributes unit arrival coherence, even at zero TKE.
  if ((overwrite_mode == coherence_average_overwrite)) {
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "compute_coherent_mixing_length::cells_1", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i                  = (cell % ngrdcol);
            const int k                  = (cell / ngrdcol);
            crossed_length_up(i, k)      = zero;
            crossed_length_down(i, k)    = zero;
            arrival_coherence_up(i, k)   = one;
            arrival_coherence_down(i, k) = one;
          });
    }
  }
  // -------------------- Upward calculation ---------------------
  // Precompute the upward entrainment update for every layer.
  {
    const int level_count =
        (gr_grid_dir_indx > 0 ? gr_grid_dir_indx + gr_k_lb_zt < gr_k_ub_zt + 1
                              : gr_grid_dir_indx + gr_k_lb_zt >= gr_k_ub_zt)
            ? (((-gr_grid_dir_indx - gr_k_lb_zt + gr_k_ub_zt) / gr_grid_dir_indx) + 1)
            : 0;
    Kokkos::parallel_for(
        "compute_coherent_mixing_length::cells_2", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int j = (gr_grid_dir_indx * (cell / ngrdcol) + gr_grid_dir_indx + gr_k_lb_zt);
          int prev_j{-1};
          real_t layer_depth{};
          real_t retention{};
          real_t response{};
          prev_j      = (-gr_grid_dir_indx + j);
          layer_depth = Kokkos::abs((gr_zt(i, j) - gr_zt(i, prev_j)));
          // Integrate entrainment analytically through a layer-linear environment.
          retention             = Kokkos::exp((-(Lscale_mu_entrain(i)) * layer_depth));
          response              = ((one - retention) / (Lscale_mu_entrain(i) * layer_depth));
          layer_retention(i, j) = retention;
          thl_entrainment(i, j) = ((thlm(i, j) - (retention * thlm(i, prev_j))) -
                                   (response * (thlm(i, j) - thlm(i, prev_j))));
          rt_entrainment(i, j)  = ((rtm(i, j) - (retention * rtm(i, prev_j))) -
                                   (response * (rtm(i, j) - rtm(i, prev_j))));
        });
  }
  // Diagnose every upward launch state and buoyancy in parallel.
  {
    const int level_count = nzt;
    Kokkos::parallel_for(
        "compute_coherent_mixing_length::cells_3", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          real_t liquid_temperature{};
          real_t saturation_rt{};
          real_t rc_parcel{};
          // Blend toward PDF component 1 to define the upward launch state.
          // Using Lscale_pdf_launch_frac > 0 is experimental, and seems to introduce
          // a tight feedback loop that can cause instability
          thl_launch(i, k) =
              (thlm(i, k) + (Lscale_pdf_launch_frac(i) * (pdf_params_thl_1(i, k) - thlm(i, k))));
          rt_launch(i, k) =
              (rtm(i, k) + (Lscale_pdf_launch_frac(i) * (pdf_params_rt_1(i, k) - rtm(i, k))));
          // Saturation-adjust the launch state before diagnosing its buoyancy.
          liquid_temperature = (thl_launch(i, k) * exner(i, k));
          saturation_rt      = sat_mixrat_liq_flatau(p_in_Pa(i, k), liquid_temperature);
          liquid_temperature = (liquid_temperature * liquid_temperature);
          rc_parcel = ((((rt_launch(i, k) - saturation_rt) * liquid_temperature) /
                        (liquid_temperature + (liquid_coef * saturation_rt))) > zero_threshold
                           ? (((rt_launch(i, k) - saturation_rt) * liquid_temperature) /
                              (liquid_temperature + (liquid_coef * saturation_rt)))
                           : zero_threshold);
          launch_buoyancy(i, k) =
              ((grav_coef(i, k) * ((thl_launch(i, k) + (vapor_coef(i, k) * rt_launch(i, k))) +
                                   (latent_coef(i, k) * rc_parcel))) -
               grav);
        });
  }
  // Compute the coherence-weighted upward length for each launch level.
  // Visit launches in a fixed order within each column for reproducible sums.
  Kokkos::parallel_for(
      "compute_coherent_mixing_length::columns_4", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column;
        int j{-1};
        int k{-1};
        real_t layer_depth{};
        real_t thl_parcel{};
        real_t rt_parcel{};
        real_t liquid_temperature{};
        real_t saturation_rt{};
        real_t rc_parcel{};
        real_t oriented_buoyancy{};
        real_t prev_buoyancy{};
        real_t prev_work{};
        real_t max_barrier{};
        real_t prev_barrier{};
        real_t layer_coherence{};
        real_t endpoint_coherence{};
        real_t layer_length{};
        real_t Lscale_tmp{};
        for (k = gr_k_lb_zt; (gr_grid_dir_indx > 0 ? k < gr_k_ub_zt + 1 : k >= gr_k_ub_zt);
             k += gr_grid_dir_indx) {
          Lscale_tmp = zero;
          // Initialize the parcel state and path-integral trackers at launch.
          thl_parcel    = thl_launch(i, k);
          rt_parcel     = rt_launch(i, k);
          prev_buoyancy = launch_buoyancy(i, k);
          prev_work     = zero;
          prev_barrier  = zero;
          // The upper boundary is outside the active mixing-length domain.
          for (j = (gr_grid_dir_indx + k);
               (gr_grid_dir_indx > 0 ? j < -gr_grid_dir_indx + gr_k_ub_zt + 1
                                     : j >= -gr_grid_dir_indx + gr_k_ub_zt);
               j += gr_grid_dir_indx) {
            layer_depth = Kokkos::abs((gr_zt(i, j) - gr_zt(i, (-gr_grid_dir_indx + j))));
            // Entrain layer-mean air into the parcel as it rises.
            thl_parcel = (thl_entrainment(i, j) + (layer_retention(i, j) * thl_parcel));
            rt_parcel  = (rt_entrainment(i, j) + (layer_retention(i, j) * rt_parcel));
            // Diagnose cloud liquid water for the parcel buoyancy calculation.
            liquid_temperature = (thl_parcel * exner(i, j));
            saturation_rt      = sat_mixrat_liq_flatau(p_in_Pa(i, j), liquid_temperature);
            liquid_temperature = (liquid_temperature * liquid_temperature);
            rc_parcel          = (((rt_parcel - saturation_rt) * liquid_temperature) /
                                  (liquid_temperature + (liquid_coef * saturation_rt)));
            rc_parcel          = (rc_parcel > zero_threshold ? rc_parcel : zero_threshold);
            // Diagnose upward-oriented buoyancy: B = g * ( thv_parcel / thvm - 1 ).
            oriented_buoyancy = ((grav_coef(i, j) * ((thl_parcel + (vapor_coef(i, j) * rt_parcel)) +
                                                     (latent_coef(i, j) * rc_parcel))) -
                                 grav);
            // Integrate buoyancy over the layer to update cumulative parcel work.
            prev_work =
                (prev_work + ((one_half * (prev_buoyancy + oriented_buoyancy)) * layer_depth));
            // The barrier is the greatest energy deficit encountered along the path.
            max_barrier = (prev_barrier > -(prev_work) ? prev_barrier : -(prev_work));
            // Integrate survival exactly for a barrier linear between these endpoints.
            layer_coherence = ((coherence_energy(i, k) * coherence_energy(i, k)) /
                               ((coherence_energy(i, k) + prev_barrier) *
                                (coherence_energy(i, k) + max_barrier)));
            layer_length    = (layer_coherence * layer_depth);
            if ((overwrite_mode == coherence_average_overwrite)) {
              endpoint_coherence = ((coherence_energy(i, k) * coherence_energy(i, k)) /
                                    ((coherence_energy(i, k) + max_barrier) *
                                     (coherence_energy(i, k) + max_barrier)));
              // These sums reconstruct remaining reach without storing trajectories.
              crossed_length_up(i, j)    = (crossed_length_up(i, j) + layer_length);
              arrival_coherence_up(i, j) = (arrival_coherence_up(i, j) + endpoint_coherence);
            }
            // Accumulate coherence-weighted distance. Coherence is one without a barrier
            // and cannot increase because the running barrier cannot decrease.
            Lscale_tmp = (Lscale_tmp + layer_length);
            // Include this layer, then neglect the remaining low-coherence tail.
            // This has almost no effect on the result, and reduces the runtime cost
            if ((layer_coherence <= real_t(0.01f))) {
              break;
            }
            prev_buoyancy = oriented_buoyancy;
            prev_barrier  = max_barrier;
          }
          Lscale_up(i, k) = Lscale_tmp;
        }
      });
  // ------------------- Downward calculation --------------------
  // Reuse the same scratch arrays for the downward layer updates.
  {
    const int level_count =
        (-1 * gr_grid_dir_indx > 0 ? -gr_grid_dir_indx + gr_k_ub_zt < gr_k_lb_zt + 1
                                   : -gr_grid_dir_indx + gr_k_ub_zt >= gr_k_lb_zt)
            ? (((gr_grid_dir_indx + gr_k_lb_zt - gr_k_ub_zt) / (-1 * gr_grid_dir_indx)) + 1)
            : 0;
    Kokkos::parallel_for(
        "compute_coherent_mixing_length::cells_5", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int j = (-gr_grid_dir_indx * (cell / ngrdcol) - gr_grid_dir_indx + gr_k_ub_zt);
          int prev_j{-1};
          real_t layer_depth{};
          real_t retention{};
          real_t response{};
          prev_j      = (gr_grid_dir_indx + j);
          layer_depth = Kokkos::abs((gr_zt(i, j) - gr_zt(i, prev_j)));
          // Integrate entrainment analytically through a layer-linear environment.
          retention             = Kokkos::exp((-(Lscale_mu_entrain(i)) * layer_depth));
          response              = ((one - retention) / (Lscale_mu_entrain(i) * layer_depth));
          layer_retention(i, j) = retention;
          thl_entrainment(i, j) = ((thlm(i, j) - (retention * thlm(i, prev_j))) -
                                   (response * (thlm(i, j) - thlm(i, prev_j))));
          rt_entrainment(i, j)  = ((rtm(i, j) - (retention * rtm(i, prev_j))) -
                                   (response * (rtm(i, j) - rtm(i, prev_j))));
        });
  }
  // Diagnose every downward launch state and buoyancy in parallel.
  {
    const int level_count = nzt;
    Kokkos::parallel_for(
        "compute_coherent_mixing_length::cells_6", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (cell / ngrdcol);
          real_t liquid_temperature{};
          real_t saturation_rt{};
          real_t rc_parcel{};
          // Blend toward PDF component 2 to define the downward launch state.
          // Using Lscale_pdf_launch_frac > 0 is experimental, and seems to introduce
          // a tight feedback loop that can cause instability
          thl_launch(i, k) =
              (thlm(i, k) + (Lscale_pdf_launch_frac(i) * (pdf_params_thl_2(i, k) - thlm(i, k))));
          rt_launch(i, k) =
              (rtm(i, k) + (Lscale_pdf_launch_frac(i) * (pdf_params_rt_2(i, k) - rtm(i, k))));
          // Saturation-adjust the launch state before diagnosing its buoyancy.
          liquid_temperature = (thl_launch(i, k) * exner(i, k));
          saturation_rt      = sat_mixrat_liq_flatau(p_in_Pa(i, k), liquid_temperature);
          liquid_temperature = (liquid_temperature * liquid_temperature);
          rc_parcel = ((((rt_launch(i, k) - saturation_rt) * liquid_temperature) /
                        (liquid_temperature + (liquid_coef * saturation_rt))) > zero_threshold
                           ? (((rt_launch(i, k) - saturation_rt) * liquid_temperature) /
                              (liquid_temperature + (liquid_coef * saturation_rt)))
                           : zero_threshold);
          launch_buoyancy(i, k) =
              (grav -
               (grav_coef(i, k) * ((thl_launch(i, k) + (vapor_coef(i, k) * rt_launch(i, k))) +
                                   (latent_coef(i, k) * rc_parcel))));
        });
  }
  // Compute the coherence-weighted downward length for each launch level.
  // Match the upward launch order so each column accumulates reproducibly.
  Kokkos::parallel_for(
      "compute_coherent_mixing_length::columns_7", ngrdcol, KOKKOS_LAMBDA(const int column) {
        const int i = column;
        int j{-1};
        int k{-1};
        real_t layer_depth{};
        real_t thl_parcel{};
        real_t rt_parcel{};
        real_t liquid_temperature{};
        real_t saturation_rt{};
        real_t rc_parcel{};
        real_t oriented_buoyancy{};
        real_t prev_buoyancy{};
        real_t prev_work{};
        real_t max_barrier{};
        real_t prev_barrier{};
        real_t layer_coherence{};
        real_t endpoint_coherence{};
        real_t layer_length{};
        real_t Lscale_tmp{};
        for (k = gr_k_lb_zt; (gr_grid_dir_indx > 0 ? k < gr_k_ub_zt + 1 : k >= gr_k_ub_zt);
             k += gr_grid_dir_indx) {
          Lscale_tmp = zero;
          // Initialize the parcel state and path-integral trackers at launch.
          thl_parcel    = thl_launch(i, k);
          rt_parcel     = rt_launch(i, k);
          prev_buoyancy = launch_buoyancy(i, k);
          prev_work     = zero;
          prev_barrier  = zero;
          // Track the parcel from its launch level down to the lower boundary.
          for (j = (-gr_grid_dir_indx + k);
               (-1 * gr_grid_dir_indx > 0 ? j < gr_k_lb_zt + 1 : j >= gr_k_lb_zt);
               j += (-1 * gr_grid_dir_indx)) {
            layer_depth = Kokkos::abs((gr_zt(i, j) - gr_zt(i, (gr_grid_dir_indx + j))));
            // Entrain layer-mean air into the parcel as it descends.
            thl_parcel = (thl_entrainment(i, j) + (layer_retention(i, j) * thl_parcel));
            rt_parcel  = (rt_entrainment(i, j) + (layer_retention(i, j) * rt_parcel));
            // Diagnose cloud liquid water for the parcel buoyancy calculation.
            liquid_temperature = (thl_parcel * exner(i, j));
            saturation_rt      = sat_mixrat_liq_flatau(p_in_Pa(i, j), liquid_temperature);
            liquid_temperature = (liquid_temperature * liquid_temperature);
            rc_parcel = ((((rt_parcel - saturation_rt) * liquid_temperature) /
                          (liquid_temperature + (liquid_coef * saturation_rt))) > zero_threshold
                             ? (((rt_parcel - saturation_rt) * liquid_temperature) /
                                (liquid_temperature + (liquid_coef * saturation_rt)))
                             : zero_threshold);
            // Reverse the physical buoyancy sign to orient it along the downward path.
            oriented_buoyancy =
                (grav - (grav_coef(i, j) * ((thl_parcel + (vapor_coef(i, j) * rt_parcel)) +
                                            (latent_coef(i, j) * rc_parcel))));
            // Integrate oriented buoyancy to update cumulative parcel work.
            prev_work =
                (prev_work + ((one_half * (prev_buoyancy + oriented_buoyancy)) * layer_depth));
            // The barrier is the greatest energy deficit encountered along the path.
            max_barrier = (prev_barrier > -(prev_work) ? prev_barrier : -(prev_work));
            // Integrate survival for a barrier linear between the layer endpoints.
            layer_coherence = ((coherence_energy(i, k) * coherence_energy(i, k)) /
                               ((coherence_energy(i, k) + prev_barrier) *
                                (coherence_energy(i, k) + max_barrier)));
            layer_length    = (layer_coherence * layer_depth);
            // Accumulate coherence-weighted distance before testing the tail cutoff.
            if (((overwrite_mode == coherence_average_overwrite) && (k != gr_k_ub_zt))) {
              endpoint_coherence           = ((coherence_energy(i, k) * coherence_energy(i, k)) /
                                              ((coherence_energy(i, k) + max_barrier) *
                                               (coherence_energy(i, k) + max_barrier)));
              crossed_length_down(i, j)    = (crossed_length_down(i, j) + layer_length);
              arrival_coherence_down(i, j) = (arrival_coherence_down(i, j) + endpoint_coherence);
            }
            Lscale_tmp = (Lscale_tmp + layer_length);
            // Include this layer, then neglect the remaining low-coherence tail.
            // This has almost no effect on the result, and reduces the runtime cost
            if ((layer_coherence <= real_t(0.01f))) {
              break;
            }
            prev_buoyancy = oriented_buoyancy;
            prev_barrier  = max_barrier;
          }
          Lscale_down(i, k) = Lscale_tmp;
        }
      });
  // -------------------- Aggregation ---------------------
  // Do some final smoothing on the result
  if ((overwrite_mode == legacy_max_overwrite)) {
    // The legacy_max_overwrite method replaces Lscale values with the
    // maximum of all parcels that crossed it. For example, if a parcel
    // that starts on k=10 rises to k=20, but one from k=11 only rises
    // to k=15, this overwriting method will boost the k=11 result as
    // if it rose to k=20.
    // This was created primarily to boost the final Lscale result,
    // but also has a smoothing effect.
    Kokkos::parallel_for(
        "compute_coherent_mixing_length::columns_8", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          int step{};
          int up_k{-1};
          int down_k{-1};
          real_t maximum_endpoint{};
          real_t minimum_endpoint{};
          maximum_endpoint =
              (gr_zt(i, gr_k_lb_zt) +
               (zlmin > Lscale_up(i, gr_k_lb_zt) ? zlmin : Lscale_up(i, gr_k_lb_zt)));
          minimum_endpoint =
              (gr_zt(i, gr_k_ub_zt) -
               (zlmin > Lscale_down(i, gr_k_ub_zt) ? zlmin : Lscale_down(i, gr_k_ub_zt)));
          for (step = 0; step < nzt - 1; step += 1) {
            up_k               = (gr_grid_dir_indx * step + gr_k_lb_zt);
            down_k             = (-gr_grid_dir_indx * step - gr_grid_dir_indx + gr_k_ub_zt);
            Lscale_up(i, up_k) = (zlmin > Lscale_up(i, up_k) ? zlmin : Lscale_up(i, up_k));
            Lscale_down(i, down_k) =
                (zlmin > Lscale_down(i, down_k) ? zlmin : Lscale_down(i, down_k));
            // Match the legacy active bounds at the surface and upper boundary.
            if (((step > 0) && (step < nzt - 2))) {
              Lscale_up(i, up_k) = (Lscale_up(i, up_k) > (maximum_endpoint - gr_zt(i, up_k))
                                        ? Lscale_up(i, up_k)
                                        : (maximum_endpoint - gr_zt(i, up_k)));
              maximum_endpoint   = (maximum_endpoint > (gr_zt(i, up_k) + Lscale_up(i, up_k))
                                        ? maximum_endpoint
                                        : (gr_zt(i, up_k) + Lscale_up(i, up_k)));
            }
            if ((step < nzt - 2)) {
              Lscale_down(i, down_k) =
                  (Lscale_down(i, down_k) > (gr_zt(i, down_k) - minimum_endpoint)
                       ? Lscale_down(i, down_k)
                       : (gr_zt(i, down_k) - minimum_endpoint));
              minimum_endpoint = (minimum_endpoint < (gr_zt(i, down_k) - Lscale_down(i, down_k))
                                      ? minimum_endpoint
                                      : (gr_zt(i, down_k) - Lscale_down(i, down_k)));
            }
          }
        });
  } else if ((overwrite_mode == coherence_average_overwrite)) {
    // The coherence_average_overwrite averages how far arriving parcels can still
    // travel, weighted by their coherence at the level.
    // For an upward example on an evenly spaced grid, suppose the mixing lengths
    // from k=10,11,12,13 would reach k=20, and those from k=14,15 would reach k=18.
    // If all arriving parcels have coherence 1 and none from below k=10 arrive,
    // the averaged lengths at k=14,15 correspond to reaching
    // k=(4/5)*20+(1/5)*18 and k=(4/6)*20+(2/6)*18 respectively.
    Kokkos::parallel_for(
        "compute_coherent_mixing_length::columns_9", ngrdcol, KOKKOS_LAMBDA(const int column) {
          const int i = column;
          int step{};
          int up_k{-1};
          int down_k{-1};
          real_t remaining_length_up{};
          real_t remaining_length_down{};
          remaining_length_up   = Lscale_up(i, gr_k_lb_zt);
          remaining_length_down = Lscale_down(i, (-gr_grid_dir_indx + gr_k_ub_zt));
          for (step = 0; step < nzt - 1; step += 1) {
            up_k   = (gr_grid_dir_indx * step + gr_k_lb_zt);
            down_k = (-gr_grid_dir_indx * step - gr_grid_dir_indx + gr_k_ub_zt);
            if ((step > 0)) {
              remaining_length_up =
                  ((remaining_length_up - crossed_length_up(i, up_k)) + Lscale_up(i, up_k));
              remaining_length_down = ((remaining_length_down - crossed_length_down(i, down_k)) +
                                       Lscale_down(i, down_k));
            }
            Lscale_up(i, up_k) = (zlmin > (remaining_length_up / arrival_coherence_up(i, up_k))
                                      ? zlmin
                                      : (remaining_length_up / arrival_coherence_up(i, up_k)));
            Lscale_down(i, down_k) =
                (zlmin > (remaining_length_down / arrival_coherence_down(i, down_k))
                     ? zlmin
                     : (remaining_length_down / arrival_coherence_down(i, down_k)));
          }
        });
  } else if ((overwrite_mode == no_overwrite)) {
    // No overwrite, just ensure zlmin is respected
    {
      const int level_count = nzt;
      Kokkos::parallel_for(
          "compute_coherent_mixing_length::cells_10", (level_count * ngrdcol),
          KOKKOS_LAMBDA(const int cell) {
            const int i       = (cell % ngrdcol);
            const int k       = (cell / ngrdcol);
            Lscale_up(i, k)   = (zlmin > Lscale_up(i, k) ? zlmin : Lscale_up(i, k));
            Lscale_down(i, k) = (zlmin > Lscale_down(i, k) ? zlmin : Lscale_down(i, k));
          });
    }
  }
  // Apply surface constraints and combine the two directional lengths.
  {
    const int level_count =
        (gr_grid_dir_indx > 0 ? gr_k_lb_zt < -gr_grid_dir_indx + gr_k_ub_zt + 1
                              : gr_k_lb_zt >= -gr_grid_dir_indx + gr_k_ub_zt)
            ? (((-gr_grid_dir_indx - gr_k_lb_zt + gr_k_ub_zt) / gr_grid_dir_indx) + 1)
            : 0;
    Kokkos::parallel_for(
        "compute_coherent_mixing_length::cells_11", (level_count * ngrdcol),
        KOKKOS_LAMBDA(const int cell) {
          const int i = (cell % ngrdcol);
          const int k = (gr_grid_dir_indx * (cell / ngrdcol) + gr_k_lb_zt);
          real_t surface_floor{};
          // The minimum length decays linearly over the lowest 500 m.
          if (l_implemented) {
            surface_floor =
                (((surface_layer_depth - (gr_zt(i, k) - gr_zm(i, gr_k_lb_zm))) * lmin(i)) *
                 invs_sfclyr_depth);
          } else {
            surface_floor = (((surface_layer_depth - gr_zt(i, k)) * lmin(i)) * invs_sfclyr_depth);
          }
          // Apply the same surface floor before combining the directional lengths.
          Lscale_up(i, k) = (surface_floor > Lscale_up(i, k) ? surface_floor : Lscale_up(i, k));
          Lscale_down(i, k) =
              (surface_floor > Lscale_down(i, k) ? surface_floor : Lscale_down(i, k));
          // Weighted geometric mean between Lscale_up and Lscale_down.
          Lscale(i, k) = (Kokkos::pow(Lscale_up(i, k), Lscale_up_weight(i)) *
                          Kokkos::pow(Lscale_down(i, k), (one - Lscale_up_weight(i))));
          // Enforce the configured upper bound on the combined mixing length.
          Lscale(i, k) = (Lscale(i, k) < Lscale_max(i) ? Lscale(i, k) : Lscale_max(i));
          // The highest active level also supplies the zero-gradient lid value.
          if ((k == -gr_grid_dir_indx + gr_k_ub_zt)) {
            Lscale(i, gr_k_ub_zt) = Lscale(i, k);
          }
        });
  }
  return;
}

} // namespace

namespace {

KOKKOS_INLINE_FUNCTION real_t sat_mixrat_liq_flatau(real_t pressure, real_t temperature // In
) {
  // Description:
  // Compute liquid saturation mixing ratio using only the Flatau polynomial.
  // This compact specialization is colocated with the mixing-length closure
  // so the compiler can inline it into parcel trajectories.
  // It is also intentionally the serial method as opposed to the
  // complex pipeline / multi-fma unit friendlier version. This was
  // to minimize register usage in the complex loop it is used in.
  // See saturation.F90 for more details.
  // Note:
  // It is usually frowned upon to make custom saturation equations like
  // this because it can cause different parts of code to detect cloud
  // transitions differently, leading to inconsistent calculations.
  // We believe that there is no similar concern in this instance
  // because this is used only to compute the behavior of "virtual" parcels
  // in the above methods. Virtual in the sense that their lifetime is
  // only inside the Lscale calculation, which is already somewhat
  // arbitrary, and none of that internal information escapes.
  // Inputs and output:
  // pressure [Pa] and temperature [K] are scalar inputs; saturation_rt [kg/kg]
  // is returned. Evaluate the polynomial no lower than -85 deg C. If pressure
  // minus saturation vapor pressure is less than 1 Pa, return ep (the ratio
  // of dry-air to water-vapor gas constants) instead of dividing.
  // References:
  // Flatau et al. (1992), J. Appl. Meteor., 31, 1507--1513.
  // saturation.F90: shared saturation formulas and polynomial evaluations.
  // -----------------------------------------------------------------------

  // Arguments (Native order)
  // pressure: Pressure [Pa]
  // temperature: Temperature [K]

  real_t saturation_rt{};                                // Saturation mixing ratio [kg/kg]
  constexpr real_t minimum_temperature = -(real_t(85.)); // Polynomial lower bound [deg C]
  real_t temperature_celsius{};                          // Clipped temperature [deg C]
  real_t saturation_pressure{};                          // Saturation vapor pressure [Pa]

  temperature_celsius =
      ((temperature - T_freeze_K) > minimum_temperature ? (temperature - T_freeze_K)
                                                        : minimum_temperature);
  saturation_pressure =
      (real_t(100.) *
       (real_t(6.11583699) +
        (temperature_celsius *
         (real_t(0.444606896) +
          (temperature_celsius *
           (real_t(0.143177157e-1) +
            (temperature_celsius *
             (real_t(0.264224321e-3) +
              (temperature_celsius *
               (real_t(0.299291081e-5) +
                (temperature_celsius *
                 (real_t(0.203154182e-7) +
                  (temperature_celsius *
                   (real_t(0.702620698e-10) +
                    (temperature_celsius *
                     (real_t(0.379534310e-13) -
                      (temperature_celsius * real_t(0.321582393e-15))))))))))))))))));
  if (((pressure - saturation_pressure) < real_t(1.))) {
    saturation_rt = ep;
  } else {
    saturation_rt = ((ep * saturation_pressure) / (pressure - saturation_pressure));
  }
  return saturation_rt;
}

} // namespace

} // namespace clubb_kokkos
