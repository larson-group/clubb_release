// Bind and append the common Native radiation namelist settings.
//
// read_radiation_setting supplies the 22 Native members to the shared namelist
// parser, preserving incoming defaults and scalar kinds. The lookup arrays are
// borrowed as HostSpace Views over the existing 20-entry case-owned tables.
// append_radiation_setting retains Native text order and double-to-core print
// conversions through the shared text writer. Defaults, scheme selection and
// lookup counting remain in init_radiation.
#include "Radiation/radiation_io.hpp"

#include "namelist.hpp"
#include "text_writer.hpp"

#include <Kokkos_Core.hpp>
#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace clubb_kokkos {

void read_radiation_setting(int iunit, const std::string& namelist_file,
                            parameters_radiation_state& parameters_radiation,
                            bool& l_calc_thlp2_rad) {
  // Description:
  // Read &radiation_setting from initialized caller values. Native source:
  // bugsrad_driver::init_radiation; this common configuration also supplies
  // every non-BUGS radiation scheme.
  // A call-local C++ stream replaces the Native logical unit argument.
  (void)iunit;
  using HostList = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace,
                                Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  namelist::Group group("radiation_setting");
  group.bind("rad_scheme", parameters_radiation.rad_scheme, 20);
  group.bind("sol_const", parameters_radiation.sol_const);
  group.bind("alvdr", parameters_radiation.alvdr);
  group.bind("alvdf", parameters_radiation.alvdf);
  group.bind("alndr", parameters_radiation.alndr);
  group.bind("alndf", parameters_radiation.alndf);
  group.bind("kappa", parameters_radiation.kappa);
  group.bind("F0", parameters_radiation.F0);
  group.bind("F1", parameters_radiation.F1);
  group.bind("eff_drop_radius", parameters_radiation.eff_drop_radius);
  group.bind("gc", parameters_radiation.gc);
  group.bind("omega", parameters_radiation.omega);
  group.bind("Fs_values", HostList(parameters_radiation.Fs_values.data(), 20));
  group.bind("cos_solar_zen_values",
             HostList(parameters_radiation.cos_solar_zen_values.data(), 20));
  group.bind("cos_solar_zen_times", HostList(parameters_radiation.cos_solar_zen_times.data(), 20));
  group.bind("radiation_top", parameters_radiation.radiation_top);
  group.bind("l_fix_cos_solar_zen", parameters_radiation.l_fix_cos_solar_zen);
  group.bind("l_sw_radiation", parameters_radiation.l_sw_radiation);
  group.bind("slr", parameters_radiation.slr);
  group.bind("l_rad_above_cloud", parameters_radiation.l_rad_above_cloud);
  group.bind("l_use_default_std_atmosphere", parameters_radiation.l_use_default_std_atmosphere);
  group.bind("l_calc_thlp2_rad", l_calc_thlp2_rad);
  if (group.read_file(namelist_file) != 0)
    throw std::runtime_error(group.error());
}

void append_radiation_setting(int iunit, const std::string& case_info_file,
                              const parameters_radiation_state& parameters_radiation) {
  // Description:
  // Append the Native radiation setup records and emit the same values to
  // stdout. Each complete 20-entry table is printed in its Native order.
  // The Native STATUS='old' contract requires a previously created setup file.
  (void)iunit;
  if (!std::filesystem::exists(case_info_file))
    throw std::runtime_error("Radiation setup file does not exist: " + case_info_file);
  std::ofstream case_log(case_info_file, std::ios::app);
  if (!case_log)
    throw std::runtime_error("Cannot append radiation settings: " + case_info_file);

  write_text("--------------------------------------------------", case_log);
  write_text("&radiation_setting", case_log);
  write_text("--------------------------------------------------", case_log);

  write_text("rad_scheme = " + parameters_radiation.rad_scheme, case_log);
  write_text("sol_const = ", static_cast<real_t>(parameters_radiation.sol_const), case_log);
  write_text("alvdr = ", static_cast<real_t>(parameters_radiation.alvdr), case_log);
  write_text("alvdf = ", static_cast<real_t>(parameters_radiation.alvdf), case_log);
  write_text("alndr = ", static_cast<real_t>(parameters_radiation.alndr), case_log);
  write_text("alndf = ", static_cast<real_t>(parameters_radiation.alndf), case_log);
  write_text("radiation_top = ", parameters_radiation.radiation_top, case_log);
  write_text("F0 = ", parameters_radiation.F0, case_log);
  write_text("F1 = ", parameters_radiation.F1, case_log);
  write_text("kappa = ", parameters_radiation.kappa, case_log);
  write_text("gc = ", parameters_radiation.gc, case_log);
  write_text("omega = ", parameters_radiation.omega, case_log);
  write_text("slr = ", static_cast<real_t>(parameters_radiation.slr), case_log);
  write_text("l_rad_above_cloud = ", parameters_radiation.l_rad_above_cloud, case_log);
  write_text("l_sw_radiation = ", parameters_radiation.l_sw_radiation, case_log);
  write_text("l_fix_cos_solar_zen = ", parameters_radiation.l_fix_cos_solar_zen, case_log);
  write_text("l_use_default_std_atmosphere = ", parameters_radiation.l_use_default_std_atmosphere,
             case_log);
  write_text("Fs_values = ", parameters_radiation.Fs_values, case_log);
  write_text("cos_solar_zen_values = ", parameters_radiation.cos_solar_zen_values, case_log);
  write_text("cos_solar_zen_times = ", parameters_radiation.cos_solar_zen_times, case_log);
  write_text("eff_drop_radius = ", parameters_radiation.eff_drop_radius, case_log);
}

} // namespace clubb_kokkos
