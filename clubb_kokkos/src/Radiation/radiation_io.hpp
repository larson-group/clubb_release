// Read and report one case's radiation settings in C++.
//
// read_radiation_setting updates initialized caller values from the Native
// &radiation_setting group. append_radiation_setting writes the Native ordered
// setup records to stdout and the existing case-info file. The three lookup
// arrays have 20 entries; case ownership and initialization stay in
// parameters_radiation_state/init_radiation. Neither operation retains state.
#pragma once

#include "Radiation/parameters_radiation.hpp"

namespace clubb_kokkos {

void read_radiation_setting(int iunit, const std::string& namelist_file,      // In
                            parameters_radiation_state& parameters_radiation, // InOut
                            bool& l_calc_thlp2_rad);                          // InOut

void append_radiation_setting(int iunit, const std::string& case_info_file,            // In
                              const parameters_radiation_state& parameters_radiation); // In

} // namespace clubb_kokkos
