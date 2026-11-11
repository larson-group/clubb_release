// Radiation orchestration and diagnostics; see the matching public interface.
#include "scratch_workspace.hpp"
#include "radiation_module.hpp"
#include "cos_solar_zen_module.hpp"
#include "rad_lwsw_module.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "CLUBB_core/stats_netcdf.hpp"

#include <iostream>
#include <stdexcept>

namespace clubb_kokkos {
void advance_clubb_radiation(
    ScratchWorkspace& workspace,                                                    // InOut
    const grid& gr, int ngrdcol, int hydromet_dim, int pdf_dim, int lh_num_samples, // In
    bool l_rad_itime, real_t dt_main, int day, int month, int year,                 // In
    real_t lat_vals, real_t lon_vals, double time_current, double time_initial,     // In
    Real2D rho, Real2D rho_zm, Real2D p_in_Pa, Real2D exner,                        // In
    Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D p_sfc,                              // In
    Real2D cloud_frac, Real2D ice_supersat_frac,                                    // In
    Real2D thlm, Real2D rtm, Real2D rcm,                                            // In
    Real4D X_nl_all_levs, Real3D lh_rt_clipped, Real3D lh_thl_clipped, Real3D lh_rc_clipped, // In
    Real3D lh_sample_point_weights, Real3D hydromet, const hm_metadata_type& hm_metadata,    // In
    stats_type& stats, err_info_type& err_info,                         // InOut
    Real1D deep_soil_T_in_K, Real1D sfc_soil_T_in_K, Real1D veg_T_in_K, // InOut
    Real2D radht,                                                       // InOut
    const parameters_radiation_state& parameters_radiation,             // Native module In
    radiation_variables_state& radiation_variables,                     // Native module InOut
    const simple_rad_state& simple_rad,                                 // Native module In
    const soil_vegetation_state& soil_vegetation                        // Native module In
) {
  // Advance the active scheme when needed and update radiation stats on each
  // sampling timestep. Retained fluxes are also inputs to soil/vegetation.
  auto workspace_scope    = workspace.scope();
  const auto Frad         = radiation_variables.Frad;
  const auto Frad_SW_up   = radiation_variables.Frad_SW_up;
  const auto Frad_LW_up   = radiation_variables.Frad_LW_up;
  const auto Frad_SW_down = radiation_variables.Frad_SW_down;
  const auto Frad_LW_down = radiation_variables.Frad_LW_down;

  // Soil/vegetation is kept here because the only present user is GABLS3,
  // which also uses radiation. The call consumes previous-step surface fluxes.
  if (soil_vegetation.l_soil_veg) {
    advance_soil_veg(workspace, ngrdcol, dt_main, Kokkos::subview(rho_zm, Kokkos::ALL, 0), // In
                     Kokkos::subview(Frad_SW_up, Kokkos::ALL, 0),                          // In
                     Kokkos::subview(Frad_SW_down, Kokkos::ALL, 0),                        // In
                     Kokkos::subview(Frad_LW_down, Kokkos::ALL, 0),                        // In
                     wpthlp_sfc, wprtp_sfc, p_sfc,                                         // In
                     stats, deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K);                // InOut
  }

  if (l_rad_itime) {
    // Snow/ice have been updated by microphysics, but thlm/rtm have not. Moving
    // radiation before microphysics would change that Native ordering.
    // Native l_silhs_rad remains rejected at initialization, as in actual JAX.
    // Its explicit public entry below also rejects it; no Native computation
    // is used as a fallback for an unported radiation call.
    radiation_driver(workspace, gr, time_current, time_initial, hydromet_dim, ngrdcol, // In
                     day, month, year, lat_vals, lon_vals,                             // In
                     rho, rho_zm, p_in_Pa, exner, cloud_frac, ice_supersat_frac,       // In
                     thlm, rtm, rcm, hydromet, hm_metadata,                            // In
                     stats, err_info,                                                  // InOut
                     radht, Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down,  // Out
                     parameters_radiation, radiation_variables, simple_rad);

    if (clubb_at_least_debug_level_api(0)) {
      if (any_equal(err_info.err_code, clubb_fatal_error)) {
        std::cerr << err_info.err_header_global << '\n'
                  << "Fatal error in radiation, check your parameter values and timestep\n";
        return;
      }
    }
  }

  // Update on every sampled timestep, even if radiation was not recomputed.
  if (stats.l_sample) {
    update_radiation_variables(ngrdcol, gr.nzm, gr.nzt, radht, Frad, Frad_SW_up, Frad_LW_up, // In
                               Frad_SW_down, Frad_LW_down,                                   // In
                               stats, // InOut
                               parameters_radiation, radiation_variables);
  }
}

void radiation_driver(ScratchWorkspace& workspace, const grid& gr, double time_current,
                      double time_initial,
                      int hydromet_dim,                                                // In
                      int ngrdcol,                                                     // In
                      int day, int month, int year,                                    // In
                      real_t lat_vals, real_t lon_vals,                                // In
                      Real2D rho, Real2D rho_zm, Real2D p_in_Pa,                       // In
                      Real2D exner, Real2D cloud_frac, Real2D ice_supersat_frac,       // In
                      Real2D thlm, Real2D rtm, Real2D rcm, Real3D hydromet,            // In
                      const hm_metadata_type& hm_metadata,                             // In
                      stats_type& stats,                                               // InOut
                      err_info_type& err_info,                                         // InOut
                      Real2D radht, Real2D Frad, Real2D Frad_SW_up, Real2D Frad_LW_up, // Out
                      Real2D Frad_SW_down, Real2D Frad_LW_down,                        // Out
                      const parameters_radiation_state& parameters_radiation, // Native module In
                      radiation_variables_state& radiation_variables,         // Native module InOut
                      const simple_rad_state& simple_rad                      // Native module In
) {
  // Compute a radiation tendency, keeping initialization, solar geometry,
  // scheme selection, component diagnostics and fatal handling in source order.
  auto workspace_scope = workspace.scope();
  const int nzm = gr.nzm, nzt = gr.nzt;
  const auto radht_LW = radiation_variables.radht_LW;
  const auto radht_SW = radiation_variables.radht_SW;
  const auto Frad_SW  = radiation_variables.Frad_SW;
  const auto Frad_LW  = radiation_variables.Frad_LW;
  real_t Fs0, amu0_core_rknd;
  double amu0;
  int amu0_index;
  constexpr bool l_center = true; // Centered interpolation in sunray_sw.

  // Native TRIM removes trailing spaces only.
  std::string rad_scheme = parameters_radiation.rad_scheme;
  const auto last        = rad_scheme.find_last_not_of(' ');
  if (last == std::string::npos)
    rad_scheme.clear();
  else
    rad_scheme.erase(last + 1);

  // Initialize only the six Native outputs here. Component fields are reset
  // inside the simplified scheme, not unconditionally as in the JAX adaptation.
  Kokkos::parallel_for(
      "radiation_driver::zero_fluxes", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol, k = cell / ngrdcol;
        Frad(i, k)         = real_t(0);
        Frad_SW_up(i, k)   = real_t(0);
        Frad_LW_up(i, k)   = real_t(0);
        Frad_SW_down(i, k) = real_t(0);
        Frad_LW_down(i, k) = real_t(0);
      });
  Kokkos::deep_copy(radht, real_t(0));

  if (parameters_radiation.l_sw_radiation) {
    if (parameters_radiation.l_fix_cos_solar_zen) {
      if (parameters_radiation.nparam > 1) {
        using HostTable = Kokkos::View<const real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
        const HostTable times(parameters_radiation.cos_solar_zen_times.data(),
                              parameters_radiation.nparam);
        // Existing binary_search returns the zero-based Native upper bracket,
        // including its special first-interval behavior; do not use searchsorted.
        amu0_index = binary_search(parameters_radiation.nparam, times, real_t(time_current));
      } else {
        amu0_index = 0;
      }
      if (amu0_index != -1) {
        amu0 = double(parameters_radiation.cos_solar_zen_values[amu0_index]);
      } else {
        std::cerr << "Time not found in cos_solar_zen_times\n";
        throw std::runtime_error("Critical error.");
      }
    } else {
      amu0 = cos_solar_zen(day, month, year, time_current, lat_vals, lon_vals);
    }
  } else {
    amu0 = 0.0;
  }

  if (rad_scheme == "bugsrad") {
    // User-approved engine exclusion. Case initialization rejects this scheme.
    throw std::runtime_error("BUGSrad is outside the supported C++ radiation scope");
  } else if (rad_scheme == "simplified") {
    Kokkos::parallel_for(
        "radiation_driver::zero_component_fluxes", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          Frad_SW(i, k) = real_t(0);
          Frad_LW(i, k) = real_t(0);
        });
    Kokkos::parallel_for(
        "radiation_driver::zero_component_heating", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          radht_SW(i, k) = real_t(0);
          radht_LW(i, k) = real_t(0);
        });

    if (parameters_radiation.l_sw_radiation && amu0 > 0.0) {
      amu0_core_rknd = real_t(amu0);
      if (parameters_radiation.nparam > 1) {
        using HostTable = Kokkos::View<const real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
        const HostTable cosines(parameters_radiation.cos_solar_zen_values.data(),
                                parameters_radiation.nparam);
        const HostTable fluxes(parameters_radiation.Fs_values.data(), parameters_radiation.nparam);
        lin_interpolate_on_grid_api(parameters_radiation.nparam, cosines, fluxes, amu0_core_rknd,
                                    Fs0);
      } else {
        Fs0 = parameters_radiation.Fs_values[0];
      }
      sunray_sw(workspace, ngrdcol, nzt, rcm, rho, amu0_core_rknd, gr.dzt,                // In
                gr.zm, gr.zt,                                                             // In
                parameters_radiation.eff_drop_radius, real_t(parameters_radiation.alvdr), // In
                parameters_radiation.gc, Fs0, parameters_radiation.omega, l_center,       // In
                Frad_SW);                                                                 // Out
      const Real2D radht_SW_ddzm = ddzm(nzm, nzt, ngrdcol, gr, Frad_SW);
      Kokkos::parallel_for(
          "radiation_driver::SW_heating", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
          KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol, k = cell / ngrdcol;
            radht_SW(i, k) = -radht_SW_ddzm(i, k) / (rho(i, k) * Cp);
          });
    }

    clubb_kokkos::simple_rad(workspace, gr, ngrdcol, rho, rho_zm, rtm, rcm, exner, // In
                             stats,                                                // InOut
                             Frad_LW, radht_LW,                                    // Out
                             parameters_radiation);
    Kokkos::parallel_for(
        "radiation_driver::total_flux", Kokkos::RangePolicy<>(0, ngrdcol * nzm),
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          Frad(i, k) = Frad_SW(i, k) + Frad_LW(i, k);
        });
    Kokkos::parallel_for(
        "radiation_driver::total_heating", Kokkos::RangePolicy<>(0, ngrdcol * nzt),
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol, k = cell / ngrdcol;
          radht(i, k) = radht_SW(i, k) + radht_LW(i, k);
        });
  } else if (rad_scheme == "simplified_bomex") {
    simple_rad_bomex(workspace, gr, ngrdcol, radht);
    Kokkos::deep_copy(Frad, real_t(0));
  } else if (rad_scheme == "lba") {
    simple_rad_lba(workspace, gr, ngrdcol, time_current, time_initial, radht, simple_rad);
    Kokkos::deep_copy(Frad, real_t(0));
  } else {
    // Native normally skips this dispatch for rad_scheme='none'. Retain its
    // fatal default rather than the extra 'none' acceptance in actual JAX.
    std::cerr << "Undefined value for namelist variable rad_scheme: " << rad_scheme << '\n';
    throw std::runtime_error("Fatal error encountered in radiation_driver.");
  }

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      std::cerr << err_info.err_header_global << '\n' << "Fatal error in radiation_driver:\n";
      return;
    }
  }
}

void update_radiation_variables(
    int ngrdcol, int nzm, int nzt, Real2D radht, Real2D Frad, Real2D Frad_SW_up, // In
    Real2D Frad_LW_up,                                                           // In
    Real2D Frad_SW_down, Real2D Frad_LW_down,                                    // In
    stats_type& stats,                                                           // InOut
    const parameters_radiation_state& parameters_radiation,                      // Native module In
    const radiation_variables_state& radiation_variables                         // Native module In
) {
  // Update sampled radiation variables in the exact Native name/call order.
  if (!stats.l_sample)
    return;
  stats_update("Frad", Frad, stats);
  stats_update("radht", radht, stats);
  stats_update("Frad_SW_up", Frad_SW_up, stats);
  stats_update("Frad_LW_up", Frad_LW_up, stats);
  stats_update("Frad_SW_down", Frad_SW_down, stats);
  stats_update("Frad_LW_down", Frad_LW_down, stats);

  std::string rad_scheme = parameters_radiation.rad_scheme;
  const auto last        = rad_scheme.find_last_not_of(' ');
  if (last == std::string::npos)
    rad_scheme.clear();
  else
    rad_scheme.erase(last + 1);
  if (rad_scheme == "simplified" || rad_scheme == "bugsrad") {
    stats_update("radht_LW", radiation_variables.radht_LW, stats);
    stats_update("radht_SW", radiation_variables.radht_SW, stats);
    stats_update("Frad_SW", radiation_variables.Frad_SW, stats);
    stats_update("Frad_LW", radiation_variables.Frad_LW, stats);
  }
  // The l_output_rad_files/BUGS-only extended diagnostics are excluded with the
  // BUGSrad engine and its existing initialization gate, not Native fallbacks.
}

void silhs_radiation_driver(const grid& gr, int ngrdcol, int nzm, int nzt, int lh_num_samples,
                            int pdf_dim,                                           // In
                            int hydromet_dim, const hm_metadata_type& hm_metadata, // In
                            int day, int month, int year,                          // In
                            real_t lat_vals, real_t lon_vals,                      // In
                            double time_current, double time_initial, Real2D rho, Real2D rho_zm,
                            Real2D p_in_Pa,                                                    // In
                            Real2D exner,                                                      // In
                            Real2D cloud_frac, Real2D ice_supersat_frac, Real4D X_nl_all_levs, // In
                            Real3D lh_rt_clipped, Real3D lh_thl_clipped, Real3D lh_rc_clipped, // In
                            Real3D lh_sample_point_weights, Real3D hydromet,                   // In
                            err_info_type& err_info, // InOut
                            Real2D radht, Real2D Frad, Real2D Frad_SW_up, Real2D Frad_LW_up,
                            Real2D Frad_SW_down, // Out
                            Real2D Frad_LW_down  // Out
) {
  // Explicit existing JAX/Kokkos exclusion, retained per the requested scope.
  // Native also multiplies momentum-level sample fluxes (nzm) by thermodynamic
  // weights (nzt), a nonconforming product; no silent numerical repair is made.
  throw std::runtime_error("SILHS radiation is outside the port scope");
}

} // namespace clubb_kokkos
