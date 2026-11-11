// Radiation configuration owned by one CLUBB case.
//
// This replaces the mutable module fields in parameters_radiation.F90. Call
// init_radiation before evaluating a radiation scheme. Fields retain Native
// names, declaration order and precision; lookup tables have 20 entries.
// BUGSrad configuration is retained for shared sounding/I/O contracts, although
// the BUGSrad engine remains outside the supported radiation calculation.
#pragma once

#include "clubb_real.hpp"

#include <array>
#include <string>

namespace clubb_kokkos {

struct parameters_radiation_state {
  std::string rad_scheme; // Radiation scheme name.

  double sol_const{};     // Solar constant [W/m^2].
  real_t radiation_top{}; // Top of the atmosphere supplied to radiation [m].

  // Albedo values; alvdr is also used by the simplified schemes.
  double alvdr{}; // Visible direct surface albedo [-].
  double alndr{}; // Near-IR direct surface albedo [-].
  double alvdf{}; // Visible diffuse surface albedo [-].
  double alndf{}; // Near-IR diffuse surface albedo [-].

  // Long-wave constants for simplified radiation.
  real_t kappa{}; // Duynkerke equation 5 absorption constant [m^2/kg].
  real_t F0{};    // Cloud-top heating coefficient [W/m^2].
  real_t F1{};    // Cloud-base heating coefficient [W/m^2].

  // Short-wave constants.
  real_t eff_drop_radius{}; // Effective droplet radius [m].
  real_t gc{};              // Asymmetry parameter [-].
  real_t omega{};           // Single-scattering albedo [-].
  double slr{};             // Fraction of daylight [-].

  std::array<real_t, 20> Fs_values{};            // Incident shortwave flux lookup [W/m^2].
  std::array<real_t, 20> cos_solar_zen_times{};  // Lookup times [s].
  std::array<real_t, 20> cos_solar_zen_values{}; // Cosines of solar zenith angle [-].

  bool l_fix_cos_solar_zen{}; // Use the prescribed solar-zenith lookup.
  bool l_sw_radiation{};      // Enable shortwave radiation.
  bool l_rad_above_cloud{};   // Include the DYCOMS-II RF02 Heaviside contribution.

  // Native init_radiation assigns this only while counting non-sentinel table
  // entries. An entirely absent table leaves Native nparam undefined; zero is
  // the explicit initial case state, not a numerical correction to valid data.
  int nparam{}; // Number of supplied solar-zenith lookup entries [-].

  bool l_use_default_std_atmosphere{}; // Use the 1976 U.S. standard profile.
};

// Common initialization is defined in Native bugsrad_driver.F90, but it also
// supplies every non-BUGS scheme. Defaults, counting, namelist input and
// setup-text output all execute in C++.
void init_radiation(int iunit, const std::string& namelist_file,       // In
                    const std::string& case_info_file,                 // In
                    bool& l_calc_thlp2_rad,                            // InOut
                    parameters_radiation_state& parameters_radiation); // Native module InOut

} // namespace clubb_kokkos
