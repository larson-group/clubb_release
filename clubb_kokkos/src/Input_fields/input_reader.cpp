// SAM-like case input table reading, lookup and blank filling.
#include "input_reader.hpp"
#include "input_names.hpp"
#include "CLUBB_core/interpolation.hpp"
#include "CLUBB_core/file_functions.hpp"
#include <array>
#include <iostream>
#include <stdexcept>

namespace clubb_kokkos {

namespace {
constexpr real_t blank_value = real_t(-999.9);
// Fortran character comparisons blank-pad the shorter operand; input metadata
// has the Native30-character limit. Public C++ names drop trailing padding, preserving
// leading blanks so direct std::string comparisons retain Native semantics.
using input_reader_detail::input_text_trim;
bool input_text_equal(const std::string& a, const std::string& b) {
  return input_text_trim(a) == input_text_trim(b);
}
} // namespace

template <class Grid, class Var>
static Real1D linear_fill_blanks(int dim_grid, Grid grid, Var var, real_t default_value);

//=============================================================================
void read_two_dim_file(int iunit, int nCol, const std::string& filename, // In
                       std::vector<two_dim_read_var>& read_vars,         // Out
                       one_dim_read_var& other_dim)                      // Out
{
  //
  // Description: This subroutine reads from a file containing data that varies
  //              in two dimensions. These dimensions are typically height
  //              and time.
  //
  //-----------------------------------------------------------------------------------------------

  using file_functions_detail::formatted_file;
  using file_functions_detail::read_item;
  (void)iunit;
  read_vars.assign(nCol, {});
  other_dim = {};
  std::vector<std::string> names(nCol); // Native character(len=30) variable names.
  int nRowI = 0, nRowO = 0;             // Inner and outer row counts.
  int k, j, i;
  bool isComment = true;
  std::string tmpline;
  std::vector<real_t> tmp(nCol);
  int input_status;

  // First run through, take names and determine how large the data file is.
  formatted_file file(input_text_trim(filename));
  // Skip all the comments at the top of the file. Native examines the raw
  // first200 characters, including an exclamation mark inside quoted text.
  while (isComment) {
    if (file.read_record(tmpline) != 0)
      throw std::runtime_error("Error reading input file " + filename);
    isComment = tmpline.substr(0, 200).find('!') != std::string::npos;
  }
  // Go back to the line that wasn't a comment.
  file.backspace();
  std::vector<read_item> name_items;
  for (auto& name : names)
    name_items.emplace_back(name, 30);
  if (file.read(name_items) != 0)
    throw std::runtime_error("Error reading input file " + filename);

  while (true) {
    std::array<read_item, 2> header{read_item(tmp[0]), read_item(nRowI)};
    input_status = file.read(header);
    // If input_status shows an end of data, then exit the loop.
    if (input_status < 0) {
      break;
    } else if (input_status > 0) {
      std::cerr << "Error reading data from file: " << input_text_trim(filename) << '\n';
      throw std::runtime_error("Fatal error input_reader");
    }
    if (nRowI < 1)
      throw std::runtime_error(
          "Number of elements must be an integer and greater than zero in two-dim input file.");
    for (k = 0; k < nRowI; ++k) {
      std::vector<read_item> row;
      for (auto& value : tmp)
        row.emplace_back(value);
      if (file.read(row) != 0)
        throw std::runtime_error("Error reading input file " + filename);
    }
    ++nRowO;
  }
  // Native has no defined inner extent if the file has no data blocks.
  if (nRowO == 0)
    throw std::runtime_error("No data rows in two-dim input file " + filename);
  for (i = 0; i < nRowO; ++i) {
    file.backspace();
    for (j = 0; j < nRowI; ++j)
      file.backspace();
  }
  file.backspace();

  // Store the names into the structure and allocate accordingly.
  using HostValues = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  std::vector<HostValues> values;
  for (k = 0; k < nCol; ++k) {
    read_vars[k].name      = input_text_trim(names[k]);
    read_vars[k].dim1_name = input_text_trim(time_name);
    read_vars[k].dim2_name = input_text_trim(names[0]);
    values.emplace_back("read_two_dim_file::values", nRowI, nRowO);
  }
  other_dim.name     = input_text_trim(time_name);
  other_dim.dim_name = input_text_trim(time_name);
  Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> other_values(
      "read_two_dim_file::other_dim", nRowO);

  // Read in the data again to the newly allocated arrays.
  for (k = 0; k < nRowO; ++k) {
    std::array<read_item, 1> coordinate{read_item(other_values(k))};
    if (file.read(coordinate) != 0)
      throw std::runtime_error("Error reading input file " + filename);
    for (j = 0; j < nRowI; ++j) {
      std::vector<read_item> row;
      for (i = 0; i < nCol; ++i)
        row.emplace_back(values[i](j, k));
      if (file.read(row) != 0)
        throw std::runtime_error("Error reading input file " + filename);
    }
  }
  // Keep Native(inner row,outer row), commonly(height,time), in the public Views.
  for (k = 0; k < nCol; ++k)
    read_vars[k].values = to_execution_space(values[k]);
  other_dim.values = to_execution_space(other_values);
}

//=============================================================================
void read_one_dim_file(int iunit, int nCol, const std::string& filename, // In
                       std::vector<one_dim_read_var>& read_vars)         // Out
{
  //
  // Description:
  //   This subroutine reads from a file containing data that varies
  //   in one dimension. The dimension is typically time.
  //
  // References:
  //   None
  //----------------------------------------------------------------------------------------------

  using file_functions_detail::formatted_file;
  using file_functions_detail::read_item;
  (void)iunit;
  read_vars.assign(nCol, {});
  std::vector<std::string> names(nCol); // Native character(len=30) variable names.
  std::string tmpline;
  int nRow = 0, k, j;
  std::vector<real_t> tmp(nCol);
  bool isComment = true;
  int input_status;

  // First run through, take names and determine how large the data file is.
  formatted_file file(input_text_trim(filename));
  // Skip all the comments at the top of the file.
  while (isComment) {
    if (file.read_record(tmpline) != 0)
      throw std::runtime_error("Error reading input file " + filename);
    isComment = tmpline.substr(0, 200).find('!') != std::string::npos;
  }
  // Go back to the line that wasn't a comment.
  file.backspace();
  std::vector<read_item> name_items;
  for (auto& name : names)
    name_items.emplace_back(name, 30);
  if (file.read(name_items) != 0)
    throw std::runtime_error("Error reading input file " + filename);
  // Count up the number of rows. Native defers a positive READ status to the
  // second pass; it still counts that row before eventually reporting failure.
  while (true) {
    std::vector<read_item> row;
    for (auto& value : tmp)
      row.emplace_back(value);
    input_status = file.read(row);
    if (input_status < 0)
      break;
    ++nRow;
  }
  // Rewind that many rows, first backing over the endfile record.
  for (int record = 0; record < nRow + 1; ++record)
    file.backspace();

  // Store the names into the structure and allocate accordingly.
  using HostValues = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  std::vector<HostValues> values;
  for (k = 0; k < nCol; ++k) {
    read_vars[k].name     = input_text_trim(names[k]);
    read_vars[k].dim_name = input_text_trim(names[0]);
    values.emplace_back("read_one_dim_file::values", nRow);
  }
  // Read in the data again to the newly allocated arrays.
  for (k = 0; k < nRow; ++k) {
    std::vector<read_item> row;
    for (j = 0; j < nCol; ++j)
      row.emplace_back(values[j](k));
    if (file.read(row) != 0)
      throw std::runtime_error("Error reading input file " + filename);
  }
  for (k = 0; k < nCol; ++k)
    read_vars[k].values = to_execution_space(values[k]);
}

//=============================================================================
void fill_blanks_one_dim_vars(int num_vars,                                // In
                              std::vector<one_dim_read_var>& one_dim_vars) // InOut
{
  //
  // Description:
  //   This subroutine fills in the blank spots (signified by constant blank_value)
  //   with values linearly interpolated using the first element of the array as a
  //   guide.
  //
  // References:
  //   None
  //----------------------------------------------------------------------------------------------

  for (int i = 0; i < num_vars; ++i)
    one_dim_vars[i].values =
        linear_fill_blanks(static_cast<int>(one_dim_vars[i].values.extent(0)), // In
                           one_dim_vars[0].values, one_dim_vars[i].values,     // In
                           real_t(0.0));                                       // In
}

//=============================================================================
void fill_blanks_two_dim_vars(int num_vars, const one_dim_read_var& other_dim, // In
                              std::vector<two_dim_read_var>& two_dim_vars)     // InOut
{
  //
  // Description:
  //   This subroutine fills in the blank spots (signified by the
  //   constant blank_value with values linearly interpolated using the first
  //   element of the array and the values in the other_dim argument as a guide.
  //
  //   This is a two step process. First we assume that the other_dim values
  //   have no holes, but there are blanks for that variable across that
  //   dimension. Then we fill holes across the dimension whose values are first
  //   in the array of two_dim_vars.
  //
  //   Ex. Time is the 'other_dim' and Height in meters is the first element in
  //   two_dim_vars.
  //
  // References:
  //   None
  //----------------------------------------------------------------------------------------------

  const int dim_size = two_dim_vars[0].values.extent(0);
  const int other_dim_size = other_dim.values.extent(0);
  for (int i = 1; i < num_vars; ++i) {
    // Interpolate along main dim.
    for (int j = 0; j < other_dim_size; ++j) {
      auto grid    = Kokkos::subview(two_dim_vars[0].values, Kokkos::ALL, j);
      auto var     = Kokkos::subview(two_dim_vars[i].values, Kokkos::ALL, j);
      auto var_out = linear_fill_blanks(dim_size, grid, var, blank_value);
      Kokkos::deep_copy(var, var_out);
    }
    // Interpolate along other dim. These are the source's strided row sections;
    // no array slice is constructed inside a device kernel.
    for (int j = 0; j < dim_size; ++j) {
      auto var     = Kokkos::subview(two_dim_vars[i].values, j, Kokkos::ALL);
      auto var_out = linear_fill_blanks(other_dim_size, other_dim.values, var, blank_value);
      Kokkos::deep_copy(var, var_out);
    }
  }
}

//=============================================================================
template <class Grid, class Var>
static Real1D linear_fill_blanks(int dim_grid, Grid grid, Var var, real_t default_value) {
  //
  // Description:
  //   This function fills blanks in array var using the grid
  //   as a guide. Blank values in var are signified by being
  //   less than or equal to the constant blank_value.
  //
  // References:
  //   None
  //-----------------------------------------------------------------------------------------------

  Real1D var_out("linear_fill_blanks::var_out", dim_grid);
  Real1D temp_grid("linear_fill_blanks::temp_grid", dim_grid);
  Real1D temp_var("linear_fill_blanks::temp_var", dim_grid);
  Kokkos::View<int*> counts("linear_fill_blanks::counts", 2);
  // Essentially this code leverages the previously written zlinterp function.
  // A smaller temporary grid and var variable are being created to pass to
  // zlinterp. zlinterp then performs the work of taking the temporary var
  // array and interpolating it to the actual grid array.
  Kokkos::parallel_for(
      "linear_fill_blanks::compact", 1, KOKKOS_LAMBDA(int) {
        int amt = 0;
        bool reversed = false;
        for (int i = 0; i < dim_grid; ++i) {
          if (var(i) > blank_value) {
            ++amt;
            temp_var((amt - 1))  = var(i);
            temp_grid((amt - 1)) = grid(i);
          }
          if (i > 0) {
            if (grid(i) < grid((i - 1)))
              reversed = true;
          }
        }
        counts(0) = amt;
        counts(1) = reversed;
      });
  // The Native control decision needs only the two integer compaction results.
  const auto host_counts = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, counts);
  const int amt = host_counts(0);
  const bool reversed = host_counts(1);
  if (amt == 0) {
    Kokkos::deep_copy(var_out, default_value);
  } else if (amt < dim_grid) {
    if (reversed) {
      Real1D negative_grid("linear_fill_blanks::negative_grid", dim_grid);
      Real1D negative_temp_grid("linear_fill_blanks::negative_temp_grid", amt);
      Kokkos::parallel_for(
          "linear_fill_blanks::negate", dim_grid, KOKKOS_LAMBDA(int i) {
            negative_grid(i) = -grid(i);
            if (i < amt)
              negative_temp_grid(i) = -temp_grid(i);
          });
      var_out = zlinterp_fnc(dim_grid, amt, negative_grid, negative_temp_grid, temp_var);
    } else {
      // dim_src=amt restricts both scratch inputs to their defined Native prefix.
      var_out = zlinterp_fnc(dim_grid, amt, grid, temp_grid, temp_var);
    }
  } else {
    Kokkos::deep_copy(var_out, var);
  }
  return var_out;
}

//=============================================================================
void deallocate_one_dim_vars(int num_vars,                                // In
                             std::vector<one_dim_read_var>& one_dim_vars) // InOut
{
  //
  //  Description:
  //    This subroutine deallocates the pointer stored in
  //    one_dim_vars%value for the whole array.
  //
  //------------------------------------------------------------------------------

  for (int i = 0; i < num_vars; ++i) {
    if (one_dim_vars[i].values.is_allocated())
      one_dim_vars[i].values = Real1D{};
  }
}

//=============================================================================
void deallocate_two_dim_vars(int num_vars,                                // In
                             std::vector<two_dim_read_var>& two_dim_vars, // InOut
                             one_dim_read_var& other_dim)                 // InOut
{
  //
  // Description:
  //   This subroutine deallocates the pointer stored in
  //   two_dim_vars%value for the whole array
  //
  // References:
  //   None
  //----------------------------------------------------------------------------------------------

  for (int i = 0; i < num_vars; ++i) {
    if (two_dim_vars[i].values.is_allocated())
      two_dim_vars[i].values = Real2D{};
  }
  if (other_dim.values.is_allocated())
    other_dim.values = Real1D{};
}

//=============================================================================
Real2D read_x_table(int nvar, int xdim, int ydim, const std::string& target_name, // In
                    const std::vector<two_dim_read_var>& retVars)                 // In
{
  //
  // Description:
  //   Searches for the variable specified by target_name in the
  //   collection of retVars. If the function finds the variable then it returns
  //   it. If it does not the program using this function will exit gracefully
  //   with a warning message.
  //
  // References:
  //   None
  //-----------------------------------------------------------------------------------------------

  Real2D x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "read_x_table::x"), xdim, ydim);
  read_x_table(nvar, xdim, ydim, target_name, retVars, x);
  return x;
}

//=============================================================================
Real1D read_x_profile(int nvar, int dim_size, const std::string& target_name, // In
                      const std::vector<one_dim_read_var>& retVars,           // In
                      const std::optional<std::string>& input_file)           // In
{
  //
  //  Description:
  //    Searches for the variable specified by target_name in the
  //    collection of retVars. If the function finds the variable then it returns
  //    it. If it does not the program using this function will exit gracefully
  //    with a warning message.
  //
  // Modified by Cavyn, June 2010
  //----------------------------------------------------------------------------------------------

  // Native writes only1:size(source); sounding callers use that defined prefix.
  // Do not invent values for the unwritten tail of its dim_size result.
  Real1D x(Kokkos::view_alloc(Kokkos::WithoutInitializing, "read_x_profile::x"), dim_size);
  read_x_profile(nvar, dim_size, target_name, retVars, x, input_file);
  return x;
}

//=============================================================================
int get_target_index(int nvar, const std::string& target_name,     // In
                     const std::vector<one_dim_read_var>& retVars) // In
{
  //
  // Description:
  //   Returns the index of the variable specified by target_name in the
  //   collection of retVars. Returns -1 if variable does not exist in retVars
  //
  // References:
  //   None
  //
  // Created by Cavyn, July 2010
  //----------------------------------------------------------------------------------------------

  for (int i = 0; i < nvar; ++i) {
    if (input_text_equal(retVars[i].name, target_name))
      return i;
  }
  return -1;
}

//=============================================================================
int count_columns(int iunit, const std::string& filename) // In
{
  // Description:
  //   This function counts the number of columns in a file, assuming that the
  //   first line of the file contains only column headers. (Comments are OK)

  // References:
  //   None

  // Created by Cavyn, July 2010
  //-----------------------------------------------------------------------------

  using file_functions_detail::formatted_file;
  using file_functions_detail::read_item;
  (void)iunit;
  int i, status_var;
  std::string tmp;
  std::array<std::string, 50> colArray{}; // Native maximum of50 column headers.
  bool isComment = true;
  formatted_file file(input_text_trim(filename));
  // Skip all the comments at the top of the file.
  while (isComment) {
    if (file.read_record(tmp) != 0)
      throw std::runtime_error("Error reading input file " + filename);
    isComment = tmp.substr(0, 200).find('!') != std::string::npos;
  }
  // Go back to the line that wasn't a comment.
  file.backspace();
  int nCols  = 0;
  status_var = file.read_record(tmp);
  if (status_var == 0) {
    // Move all words from the Native200-character internal record into an array.
    formatted_file record(std::vector<std::string>{tmp.substr(0, 200)});
    std::vector<read_item> headers;
    for (auto& name : colArray)
      headers.emplace_back(name, 200);
    status_var = record.read(headers);
  } else if (status_var > 0) {
    throw std::runtime_error(
        "Fatal error reading data in time_dependent_input function count_columns");
  }
  for (i = 0; i < static_cast<int>(colArray.size()); ++i)
    if (!input_text_trim(colArray[i]).empty())
      ++nCols;
  return nCols;
}

} // namespace clubb_kokkos
