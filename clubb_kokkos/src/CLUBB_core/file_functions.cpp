// Numeric record readers and the shared call-local formatted-file adapter.
#include "file_functions.hpp"
#include "fortran_io.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace clubb_kokkos {
namespace file_functions_detail {
formatted_file::formatted_file(const std::string& filename) {
  // Open an existing file and retain its records for Native BACKSPACE semantics.
  // Fortran file names ignore trailing character padding, retaining leading spaces.
  const auto end         = filename.find_last_not_of(' ');
  const std::string path = filename.substr(0, end == std::string::npos ? 0 : end + 1);
  std::ifstream input(path);
  if (!input)
    throw std::ios_base::failure("Error opening input file " + path);
  std::string record;
  while (std::getline(input, record)) {
    if (!record.empty() && record.back() == '\r')
      record.pop_back();
    records_.push_back(record);
  }
  if (input.bad())
    throw std::runtime_error("Error reading input file " + path);
}
formatted_file::formatted_file(std::vector<std::string> records) : records_(std::move(records)) {
  // An internal character READ has one caller-supplied record rather than a file.
}
int formatted_file::read_record(std::string& record) {
  // A formatted A read consumes one record; an EOF READ positions past ENDFILE.
  if (position_ == records_.size()) {
    after_endfile_ = true;
    return -1;
  }
  record         = records_[position_++];
  after_endfile_ = false;
  return 0;
}
int formatted_file::read(std::span<read_item> items) {
  // Complete one list-directed READ, continuing short records and discarding
  // the rest of the final record. Nulls retain destinations; slash ends the list.
  using fortran_io::Kind;
  std::string record;
  int status = read_record(record);
  if (status != 0 || items.empty())
    return status;
  const bool character_values = std::all_of(items.begin(), items.end(), [](const auto& item) {
    return std::holds_alternative<std::string*>(item.destination);
  });
  std::string source          = record;
  std::size_t token_index = 0, item_index = 0;
  bool expecting_value = true;
  while (true) {
    bool unterminated = false;
    const auto tokens = fortran_io::tokenize(source, character_values, &unterminated);
    while (token_index < tokens.size()) {
      const auto& token = tokens[token_index];
      if (token.kind == Kind::end_record) {
        // A leading comma on the next physical record denotes a new null item.
        expecting_value = true;
        ++token_index;
        continue;
      }
      if (token.kind == Kind::slash)
        return 0;
      if (token.kind == Kind::comma) {
        if (expecting_value && ++item_index == items.size())
          return 0;
        expecting_value = true;
        ++token_index;
        continue;
      }
      if (token.kind != Kind::value)
        return 1;
      int repeats                    = 1;
      bool null_value                = false;
      const fortran_io::Token* value = &token;
      if (token_index + 1 < tokens.size() && tokens[token_index + 1].kind == Kind::star &&
          token.end == tokens[token_index + 1].begin) {
        const auto star_end = tokens[token_index + 1].end;
        if (unterminated && token_index + 2 == tokens.size() && star_end < source.size() &&
            (source[star_end] == '\'' || source[star_end] == '"'))
          break;
        if (token.quoted || token.text.empty() ||
            !std::all_of(token.text.begin(), token.text.end(),
                         [](char c) { return c >= '0' && c <= '9'; }) ||
            !fortran_io::parse_value(token.text, repeats) || repeats < 1)
          return 1;
        token_index += 2;
        null_value = token_index == tokens.size() || tokens[token_index].begin != star_end ||
                     tokens[token_index].kind == Kind::comma ||
                     tokens[token_index].kind == Kind::end_record ||
                     tokens[token_index].kind == Kind::slash;
        if (!null_value) {
          if (tokens[token_index].kind != Kind::value)
            return 1;
          value = &tokens[token_index++];
        }
      } else {
        ++token_index;
      }
      for (int repeat = 0; repeat < repeats; ++repeat) {
        if (!null_value) {
          auto& item        = items[item_index];
          const bool parsed = std::visit(
              [&](auto* destination) {
                using T = std::remove_reference_t<decltype(*destination)>;
                if constexpr (std::is_same_v<T, std::string>) {
                  std::string text;
                  if (!fortran_io::parse_value(value->text, text))
                    return false;
                  *destination = text.substr(0, item.character_length);
                  return true;
                } else {
                  return !value->quoted && fortran_io::parse_value(value->text, *destination);
                }
              },
              item.destination);
          if (!parsed)
            return 1;
        }
        if (++item_index == items.size())
          return 0;
      }
      expecting_value = false;
    }
    // A required quoted value may continue across records; the shared lexer
    // returns only its complete prefix until the closing delimiter is available.
    status = read_record(record);
    if (status != 0)
      return status;
    source += '\n';
    source += record;
  }
}
void formatted_file::backspace() {
  // Native first backs over ENDFILE, then over one physical record per call.
  if (after_endfile_) {
    after_endfile_ = false;
  } else if (position_ != 0) {
    --position_;
  }
}
} // namespace file_functions_detail

//===============================================================================
void file_read_1d(int file_unit, const std::string& path_and_filename,                     // In
                  int num_datapts, int entries_per_line,                                   // In
                  Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> variable) { // Out
  // Description:
  //   Read a declared number of data points, taking entries_per_line values
  //   from each row before a possible final partial row. Values follow row order:
  //   row1,column1; row1,column2; then row2,column1, and so on.
  //
  //   For example,18 points with4 entries per row use four complete records and
  //   a final two-value record. Extra fields in each completed READ are discarded.
  //
  //            i = 1     i = 2     i = 3     i = 4
  //          ---------------------------------------
  //   k = 1 |   DP1       DP2       DP3       DP4
  //   k = 2 |   DP5       DP6       DP7       DP8
  //   k = 3 |   DP9       DP10      DP11      DP12
  //   k = 4 |   DP13      DP14      DP15      DP16
  //   k = 5 |   DP17      DP18
  //------------------------------------------------------------------------
  (void)file_unit;
  // Native serializes logical-unit access under OpenMP. Each C++ call owns its
  // file, so concurrent readers need neither a shared unit nor a global registry.
  auto file = [&]() {
    try {
      return file_functions_detail::formatted_file(path_and_filename);
    } catch (const std::ios_base::failure&) {
      std::cerr << "CLUBB encountered an error trying to open " << path_and_filename << '\n';
      throw std::runtime_error("Error opening forcings file");
    }
  }();
  // Michael Falk wrote this format for mpace_a: each full row has a declared
  // number of values; the final row contains only the remaining values (24 September 2007).
  int k, i;
  std::vector<file_functions_detail::read_item> items;
  items.reserve(entries_per_line);
  for (k = 0; k < num_datapts / entries_per_line; ++k) {
    items.clear();
    for (i = 0; i < entries_per_line; ++i)
      items.emplace_back(variable(k * entries_per_line + i));
    if (file.read(items) != 0)
      throw std::runtime_error("Error reading input file " + path_and_filename);
  }
  if (num_datapts % entries_per_line != 0) {
    k = num_datapts / entries_per_line;
    items.clear();
    for (i = 0; i < num_datapts % entries_per_line; ++i)
      items.emplace_back(variable(k * entries_per_line + i));
    if (file.read(items) != 0)
      throw std::runtime_error("Error reading input file " + path_and_filename);
  }
}

//===============================================================================
void file_read_2d(int device, const std::string& file_path,                                 // In
                  int file_dimension1, int file_dimension2, int file_per_line,              // In
                  Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> variable) { // Out
  // Description:
  //   Read file_dimension2 values on each of file_dimension1 vertical levels.
  //   Full rows have file_per_line entries; the last row on a level holds its
  //   remaining entries before the next level starts a new full row.
  //
  //   Michael Falk introduced this mpace_a format on 24 September 2007.
  //-------------------------------------------------------------------------------
  (void)device;
  Kokkos::deep_copy(Kokkos::subview(variable, std::make_pair(0, file_dimension1),
                                    std::make_pair(0, file_dimension2)),
                    real_t(-999)); // Initialize to nonsense values.
  file_functions_detail::formatted_file file(file_path);
  int i, j, k;
  std::vector<file_functions_detail::read_item> items;
  items.reserve(file_per_line);
  for (k = 0; k < file_dimension1; ++k) {
    for (j = 0; j < file_dimension2 / file_per_line; ++j) {
      items.clear();
      for (i = 0; i < file_per_line; ++i)
        items.emplace_back(variable(k, j * file_per_line + i));
      if (file.read(items) != 0)
        throw std::runtime_error("Error reading input file " + file_path);
    }
    items.clear();
    for (i = 0; i < file_dimension2 % file_per_line; ++i)
      items.emplace_back(variable(k, j * file_per_line + i));
    // Native performs this READ even when the implied-DO has no entries.
    if (file.read(items) != 0)
      throw std::runtime_error("Error reading input file " + file_path);
  }
}
} // namespace clubb_kokkos
