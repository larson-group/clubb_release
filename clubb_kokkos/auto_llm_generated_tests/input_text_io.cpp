// Check text-record semantics and grid errors outside ordinary model inputs.
#include "scratch_workspace_owner.hpp"
#include "CLUBB_core/file_functions.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/error_code.hpp"
#include "Input_fields/input_reader.hpp"
#include "Input_fields/sounding.hpp"
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
  // Report the input contract that failed, including small error-path checks.
  if (!condition)
    throw std::runtime_error(message);
}
void write(const std::filesystem::path& path, const std::string& contents) {
  // Create only the fixture records consumed by the production reader.
  std::ofstream file(path);
  file << contents;
  require(static_cast<bool>(file), "Could not write text-reader fixture");
}
void check(clubb_kokkos::ScratchWorkspace& workspace, const std::filesystem::path& directory) {
  // Use hand-written records and expected values, not a copy of parser logic.
  auto workspace_scope = workspace.scope();
  using namespace clubb_kokkos;
  using Host1D       = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  using Host2D       = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  const auto numeric = directory / "numeric.in";
  write(numeric, "1D0,2*2d0,4,999\n5 6\n7d0\n");
  Host2D padded("padded", 3, 7);
  Kokkos::deep_copy(padded, real_t(-77));
  auto values = Kokkos::subview(padded, 1, Kokkos::ALL);
  file_read_1d(27, numeric.string(), 7, 4, values);
  const std::array<real_t, 7> expected{1, 2, 2, 4, 5, 6, 7};
  for (int i = 0; i < 7; ++i) {
    require(values(i) == expected[i], "1D repeat/D-exponent/record continuation changed");
    require(padded(0, i) == -77 && padded(2, i) == -77, "Strided output touched a guard column");
  }

  const auto levels = directory / "levels.in";
  write(levels, "1 2 'ignored unfinished quote\n3 4\nignored empty READ record\n"
                "5 6\n7 8\nignored empty READ record\n");
  Host2D field("field", 2, 4);
  file_read_2d(27, levels.string(), 2, 4, 2, field);
  for (int j = 0; j < 4; ++j) {
    require(field(0, j) == j + 1 && field(1, j) == j + 5,
            "2D per-level full/empty READ record semantics changed");
  }
  using file_functions_detail::formatted_file;
  using file_functions_detail::read_item;
  std::array<real_t, 5> nullable{-9, -9, -9, -9, -9};
  std::array<read_item, 5> targets{nullable[0], nullable[1], nullable[2], nullable[3], nullable[4]};
  formatted_file nulls(std::vector<std::string>{",2*3,,/"});
  require(nulls.read(targets) == 0 && nullable == std::array<real_t, 5>{-9, 3, 3, -9, -9},
          "Null/repeated/slash input no longer retains untouched destinations");

  std::array<real_t, 3> spaced_values{-9, -9, -9};
  std::array<read_item, 3> spaced_targets{spaced_values[0], spaced_values[1], spaced_values[2]};
  formatted_file spaced_nulls(std::vector<std::string>{"2* 3"});
  require(spaced_nulls.read(spaced_targets) == 0 &&
              spaced_values == std::array<real_t, 3>{-9, -9, 3},
          "Whitespace after a repeat star must produce repeated null values");
  spaced_values.fill(-9);
  formatted_file separated_star(std::vector<std::string>{"2 *3"});
  require(separated_star.read(spaced_targets) > 0 && spaced_values[0] == 2 &&
              spaced_values[1] == -9,
          "Whitespace before a star must not form a repetition");
  spaced_values.fill(-9);
  formatted_file leading_record_comma(std::vector<std::string>{"1", ",2"});
  require(leading_record_comma.read(spaced_targets) == 0 &&
              spaced_values == std::array<real_t, 3>{1, -9, 2},
          "A leading comma on a new record must denote a null item");

  std::array<std::string, 3> names;
  std::array<read_item, 3> name_targets{names[0], names[1], names[2]};
  formatted_file character_stars(std::vector<std::string>{"a*b 2*xy"});
  require(character_stars.read(name_targets) == 0 &&
              names == std::array<std::string, 3>{"a*b", "xy", "xy"},
          "A character-value internal star became repeat syntax");

  std::array<std::string, 2> continued_names;
  std::array<read_item, 2> continued_targets{continued_names[0], continued_names[1]};
  formatted_file continued_quotes(std::vector<std::string>{"'ab", "cd' done"});
  require(continued_quotes.read(continued_targets) == 0 && continued_names[0] == "abcd" &&
              continued_names[1] == "done",
          "Quoted text no longer continues across records");

  const auto one = directory / "one.in";
  write(one, "! leading comment\nignored 'a!b' header\n"
             "z[m] 'a/b' '123456789012345678901234567890EXTRA'\n1D0 2d0 3\n4 5 6\n");
  require(count_columns(27, one.string()) == 3, "Quoted header column count changed");
  std::vector<one_dim_read_var> profiles;
  read_one_dim_file(27, 3, one.string(), profiles);
  require(profiles.size() == 3 && profiles[1].name == "a/b" &&
              profiles[2].name == "123456789012345678901234567890" &&
              profiles[1].dim_name == "z[m]",
          "Header quote/truncation/dimension metadata changed");
  auto profile = to_host_space(profiles[2].values);
  require(profile.extent(0) == 2 && profile(0) == 3 && profile(1) == 6,
          "One-dimensional table values or row count changed");

  const auto two = directory / "two.in";
  write(two, "z[m] 'a/b'\n0D0 2\n1 10\n2 20\n6d1 2\n1 30\n2 40\n");
  std::vector<two_dim_read_var> tables;
  one_dim_read_var coordinate;
  read_two_dim_file(27, 2, two.string(), tables, coordinate);
  auto table = to_host_space(tables[1].values);
  auto time  = to_host_space(coordinate.values);
  require(table.extent(0) == 2 && table.extent(1) == 2 && table(1, 0) == 20 && table(0, 1) == 30 &&
              time(0) == 0 && time(1) == 60 && tables[1].dim1_name == "Time[s]" &&
              tables[1].dim2_name == "z[m]",
          "Two-dimensional table inner/outer shape, values or metadata changed");
  write(two, "z[m] 'a/b'\n0 2.0\n1 10\n2 20\n");
  read_two_dim_file(27, 2, two.string(), tables, coordinate);
  require(tables[0].values.extent(0) == 2,
          "Native decimal-to-integer READ conversion changed for a block size");
  write(two, "z[m] 'a/b'\n0 bad\n1 10\n2 20\n");
  bool rejected = false;
  try {
    read_two_dim_file(27, 2, two.string(), tables, coordinate);
  } catch (const std::runtime_error& error) {
    rejected = std::string(error.what()) == "Fatal error input_reader";
  }
  require(rejected, "A malformed two-dimensional block size was accepted");

  err_info_type error;
  error.err_code          = decltype(error.err_code)("errors", 3);
  error.err_header_global = "test grid";
  Host1D momentum("momentum", 4), thermodynamic("thermodynamic", 3);
  const auto heights = directory / "heights.in";
  write(heights, "1D0 999\n2d0\n4\n");
  Kokkos::deep_copy(error.err_code, 0);
  read_grid_heights(4, 2, "", heights.string(), 27, momentum, thermodynamic, error);
  require(thermodynamic(0) == 1 && thermodynamic(2) == 4 && momentum(0) == 0 && momentum(3) == 0,
          "Stretched-grid scalar record reads or unused output changed");
  Kokkos::deep_copy(Kokkos::subview(error.err_code, 1), clubb_fatal_error);
  read_grid_heights(4, 2, "", heights.string(), 27, momentum, thermodynamic, error);
  auto codes = to_host_space(error.err_code);
  require(codes(0) == 0 && codes(1) == clubb_fatal_error && codes(2) == 0 && thermodynamic(0) == 0,
          "Preexisting thermodynamic-grid fatal handling changed");
  write(heights, "0\n1\n1\n4\n");
  Kokkos::deep_copy(error.err_code, 0);
  read_grid_heights(4, 3, heights.string(), "", 27, momentum, thermodynamic, error);
  codes = to_host_space(error.err_code);
  require(codes(0) == clubb_fatal_error && codes(1) == clubb_fatal_error &&
              codes(2) == clubb_fatal_error,
          "A nonincreasing grid did not fail every column");
  write(heights, "1\nnot_numeric\n4\n");
  Kokkos::deep_copy(error.err_code, 0);
  read_grid_heights(4, 2, "", heights.string(), 27, momentum, thermodynamic, error);
  codes = to_host_space(error.err_code);
  require(codes(0) == clubb_fatal_error && thermodynamic(0) == 0,
          "A malformed grid value no longer reports a fatal input error");

  const auto generic_profile = directory / "profile.in";
  write(generic_profile, "&profile nlevels=2, z=0,200, var=10,30 /\n");
  grid gr;
  gr.nzt = 2;
  gr.zt  = decltype(gr.zt)("profile heights", 1, 2);
  Host2D profile_heights("profile heights host", 1, 2);
  profile_heights(0, 0) = 50;
  profile_heights(0, 1) = 150;
  Kokkos::deep_copy(gr.zt, profile_heights);
  Real1D interpolated("interpolated profile", 2);
  read_profile(workspace, gr, generic_profile.string(), interpolated);
  const auto interpolated_host = to_host_space(interpolated);
  require(interpolated_host(0) == 15 && interpolated_host(1) == 25,
          "Generic profile namelist input or interpolation changed");
  write(generic_profile, "&profile z=0,200, var=10,30 /\n");
  rejected = false;
  try {
    read_profile(workspace, gr, generic_profile.string(), interpolated);
  } catch (const std::runtime_error&) {
    rejected = true;
  }
  require(rejected, "A generic profile with no level count was accepted");
}
} // namespace
int main(int argc, char** argv) {
  // Fixtures stay in a caller-selected scratch directory; no case data is modified.
  if (argc != 2) {
    std::cerr << "usage: input_text_io scratch_directory\n";
    return 2;
  }
  Kokkos::initialize(argc, argv);
  int result = 0;
  try {
    clubb_kokkos::test::ScratchWorkspaceOwner scratch_owner;
    auto& workspace = scratch_owner.workspace;
    std::filesystem::create_directories(argv[1]);
    check(workspace, argv[1]);
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    result = 1;
  }
  Kokkos::finalize();
  return result;
}
