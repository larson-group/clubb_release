// Read and interpret SAM-like case input tables.
//
// The records retain Native30-character name fields and array dimension order.
// values owns a Kokkos View; formatted readers truncate metadata at30 characters
// and fill call-local HostSpace buffers before transferring the scientific data.
#pragma once

#include "CLUBB_core/array_types.hpp"
#include <string>
#include <optional>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace clubb_kokkos {

// Derived type for representing a rank1 variable read from an input file.
struct one_dim_read_var {
  std::string name;     // Name of the variable, Native character(len=30).
  std::string dim_name; // Name of the dimension along which the variable varies.
  Real1D values;        // Values of that variable.
};

// Derived type for representing a rank2 variable read from an input file.
struct two_dim_read_var {
  std::string name;      // Name of the variable, Native character(len=30).
  std::string dim1_name; // Name of one of the dimensions.
  std::string dim2_name; // Name of the other dimension.
  // Native stores values(inner row,outer row): typically(height,time), even
  // though dim1_name is Time[s]. Retain that source order and LayoutLeft.
  Real2D values;
};

// Source procedure order; Native logical-unit iunit remains required. Each C++
// reader owns its input records and file handle, with no persistent unit registry.
void read_two_dim_file(int iunit, int nCol, const std::string& filename,            // In
                       std::vector<two_dim_read_var>& read_vars,                    // Out
                       one_dim_read_var& other_dim);                                // Out
void read_one_dim_file(int iunit, int nCol, const std::string& filename,            // In
                       std::vector<one_dim_read_var>& read_vars);                   // Out
void fill_blanks_one_dim_vars(int num_vars,                                         // In
                              std::vector<one_dim_read_var>& one_dim_vars);         // InOut
void fill_blanks_two_dim_vars(int num_vars, const one_dim_read_var& other_dim,      // In
                              std::vector<two_dim_read_var>& two_dim_vars);         // InOut
void deallocate_one_dim_vars(int num_vars,                                          // In
                             std::vector<one_dim_read_var>& one_dim_vars);          // InOut
void deallocate_two_dim_vars(int num_vars,                                          // In
                             std::vector<two_dim_read_var>& two_dim_vars,           // InOut
                             one_dim_read_var& other_dim);                          // InOut
Real2D read_x_table(int nvar, int xdim, int ydim, const std::string& target_name,   // In
                    const std::vector<two_dim_read_var>& retVars);                  // In
Real1D read_x_profile(int nvar, int dim_size, const std::string& target_name,       // In
                      const std::vector<one_dim_read_var>& retVars,                 // In
                      const std::optional<std::string>& input_file = std::nullopt); // In
// Return a zero-based variable slot, or -1 when the variable is absent.
int get_target_index(int nvar, const std::string& target_name,      // In
                     const std::vector<one_dim_read_var>& retVars); // In
namespace input_reader_detail {
inline std::string input_text_trim(const std::string& value) {
  const auto end = value.find_last_not_of(' ');
  return value.substr(0, end == std::string::npos ? 0 : (end + 1));
}
} // namespace input_reader_detail

// Copy the defined source prefix into existing storage; leave its tail unchanged.
template <class Output, std::enable_if_t<Kokkos::is_view<Output>::value, int> = 0>
void read_x_profile(int nvar, int dim_size, const std::string& target_name,      // In
                    const std::vector<one_dim_read_var>& retVars,                // In
                    Output x,                                                    // Out
                    const std::optional<std::string>& input_file = std::nullopt) // In
{
  (void)dim_size;
  const int i = get_target_index(nvar, target_name, retVars);
  if (i >= 0) {
    const auto values = retVars[i].values;
    const int count   = values.extent(0);
    Kokkos::deep_copy(Kokkos::subview(x, Kokkos::pair<int, int>(0, count)), values);
  } else {
    if (input_file)
      std::cerr << input_reader_detail::input_text_trim(target_name)
                << " could not be found. Check the file " << *input_file << '\n';
    else
      std::cerr << input_reader_detail::input_text_trim(target_name)
                << " could not be found. Check your sounding.in file.\n";
    throw std::runtime_error("Fatal error in read_x_profile");
  }
}

// Copy the selected table into existing storage; retain Native lookup/errors.
template <class Output>
void read_x_table(int nvar, int xdim, int ydim, const std::string& target_name,
                  const std::vector<two_dim_read_var>& retVars, Output x) {
  (void)xdim;
  (void)ydim;
  bool l_found = false;
  int i        = 0;
  while (i < nvar && !l_found) {
    if (input_reader_detail::input_text_trim(retVars[i].name) ==
        input_reader_detail::input_text_trim(target_name)) {
      l_found = true;
      Kokkos::deep_copy(x, retVars[i].values);
    }
    ++i;
  }
  if (!l_found) {
    std::cerr << input_reader_detail::input_text_trim(target_name) << " could not be found.\n";
    throw std::runtime_error("Fatal error in function read_x_table");
  }
}

int count_columns(int iunit, const std::string& filename); // In

} // namespace clubb_kokkos
