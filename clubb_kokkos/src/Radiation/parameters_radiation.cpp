// Common radiation initialization from bugsrad_driver::init_radiation.
#include "CLUBB_core/error_code.hpp"
#include "parameters_radiation.hpp"
#include "CLUBB_core/constants_clubb.hpp"
#include "Radiation/radiation_io.hpp"

#include <Kokkos_Core.hpp>

namespace clubb_kokkos {

void init_radiation(int iunit, const std::string& namelist_file,     // In
                    const std::string& case_info_file,               // In
                    bool& l_calc_thlp2_rad,                          // InOut
                    parameters_radiation_state& parameters_radiation // Native module InOut
) {
  // Set common defaults, read radiation_setting, print its metadata and count
  // prescribed solar-zenith values. No BUGSrad engine computation occurs here.
  // The explicit case owner replaces Native parameters_radiation module data.
  constexpr bool l_write_to_file = true;

  auto& rad_scheme                   = parameters_radiation.rad_scheme;
  auto& sol_const                    = parameters_radiation.sol_const;
  auto& radiation_top                = parameters_radiation.radiation_top;
  auto& alvdr                        = parameters_radiation.alvdr;
  auto& alndr                        = parameters_radiation.alndr;
  auto& alvdf                        = parameters_radiation.alvdf;
  auto& alndf                        = parameters_radiation.alndf;
  auto& F0                           = parameters_radiation.F0;
  auto& F1                           = parameters_radiation.F1;
  auto& kappa                        = parameters_radiation.kappa;
  auto& gc                           = parameters_radiation.gc;
  auto& omega                        = parameters_radiation.omega;
  auto& slr                          = parameters_radiation.slr;
  auto& l_rad_above_cloud            = parameters_radiation.l_rad_above_cloud;
  auto& l_sw_radiation               = parameters_radiation.l_sw_radiation;
  auto& l_use_default_std_atmosphere = parameters_radiation.l_use_default_std_atmosphere;
  auto& l_fix_cos_solar_zen          = parameters_radiation.l_fix_cos_solar_zen;
  auto& Fs_values                    = parameters_radiation.Fs_values;
  auto& cos_solar_zen_values         = parameters_radiation.cos_solar_zen_values;
  auto& cos_solar_zen_times          = parameters_radiation.cos_solar_zen_times;
  auto& eff_drop_radius              = parameters_radiation.eff_drop_radius;
  auto& nparam                       = parameters_radiation.nparam;

  // ---- Begin Code ----
  // Set default values, then read the namelist.
  rad_scheme = "none";

  // BUGSrad parameters are retained as configuration shared with sounding.
  sol_const     = 1367.0; // Solar constant [W/m^2].
  alvdf         = 0.1;
  alndr         = 0.1;
  alndf         = 0.1;
  radiation_top = real_t(50000); // Top of the 1976 U.S. standard profile [m].
  alvdr         = 0.1;           // Visible direct albedo used by both schemes.

  F0                           = real_t(100.0);
  F1                           = real_t(20.0);
  kappa                        = real_t(119.0);
  gc                           = real_t(0.86);
  omega                        = real_t(0.9965);
  slr                          = 1.0;
  l_rad_above_cloud            = false;
  l_sw_radiation               = false;
  l_use_default_std_atmosphere = true;
  l_fix_cos_solar_zen          = false;

  // Incoming beam-direction insolation [W/m^2] and prescribed cosine/time lists.
  Fs_values.fill(real_t(0));
  cos_solar_zen_values.fill(real_t(-999));
  cos_solar_zen_times.fill(real_t(-999));
  eff_drop_radius = real_t(1.e-5);

  read_radiation_setting(iunit, namelist_file,                    // In
                         parameters_radiation, l_calc_thlp2_rad); // InOut

  if (clubb_at_least_debug_level_api(1)) {
    // Native text_writer writes to stdout and appends radiation_setting to the
    // setup file created by the driver, in the original field/record order.
    if (l_write_to_file) {
      append_radiation_setting(iunit, case_info_file, parameters_radiation);
    }
  }

  // The loop intentionally does not reset nparam: Native leaves its old value
  // if no entries precede the sentinel. New case state starts at the documented0.
  for (int k = 0; k < int(cos_solar_zen_values.size()); ++k) {
    if (Kokkos::abs(cos_solar_zen_values[k] - real_t(-999)) <
        Kokkos::abs(cos_solar_zen_values[k] + real_t(-999)) / 2 * eps) {
      break;
    } else {
      nparam = k + 1;
    }
  }
}

} // namespace clubb_kokkos
