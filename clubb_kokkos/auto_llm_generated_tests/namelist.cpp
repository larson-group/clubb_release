// Check typed namelist syntax, Native defaults and logical View indexing.
// Real model comparisons own numerical acceptance. These small inputs exercise
// null/repeated values, sections, fixed strings, derived components and errors.
#include "namelist.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[]) {
  // Use one call-local group per read, as the permanent driver binders do.
  std::filesystem::path directory;
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    using namespace clubb_kokkos;
    directory = std::filesystem::temp_directory_path() /
                ("clubb_namelist_" +
                 std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(directory);
    const auto filename = directory / "input.nml";
    Kokkos::View<double*, Kokkos::HostSpace> values("values", 6);
    Kokkos::View<double**, Kokkos::LayoutLeft, Kokkos::HostSpace> owner("matrix_owner", 4, 3);
    auto matrix  = Kokkos::subview(owner, Kokkos::pair<int, int>(0, 2), Kokkos::ALL);
    int integer  = 7;
    bool logical = false, component_flag = false;
    double component  = 0;
    std::string label = "initial";
    std::string error;
    const auto read = [&](const std::string& source) {
      std::ofstream(filename) << source;
      namelist::Group group("test");
      group.bind("a", values);
      group.bind("matrix", matrix);
      group.bind("i", integer);
      group.bind("flag", logical);
      group.bind("label", label, 8);
      group.bind("object%x", component);
      group.bind("object%l", component_flag);
      const int status = group.read_file(filename.string());
      error            = group.error();
      return status;
    };
    Kokkos::deep_copy(values, 9.0);
    Kokkos::deep_copy(owner, -7.0);
    if (read("&OTHER ignored='x!y' /\n&TeSt A(2)=1D0,2*2.5,,6, FLAG=.TRUE., I=+3, "
             "LABEL='a!b''c', OBJECT%X=1.5-2, OBJECT%L=F / trailing ignored") != 0)
      throw std::runtime_error(error);
    const double expected[6] = {9, 1, 2.5, 2.5, 9, 6};
    for (int i = 0; i < 6; ++i)
      if (values(i) != expected[i])
        throw std::runtime_error("null/repeated/indexed assignment lost a default");
    if (!logical || integer != 3 || label != "a!b'c" || component != 0.015 || component_flag)
      throw std::runtime_error("scalar/quoted/derived component assignment mismatch");

    if (read("$test a(2:4)=3*4, matrix(:,2)=11,12, label='0123456789' $end") != 0)
      throw std::runtime_error(error);
    if (values(1) != 4 || values(2) != 4 || values(3) != 4 || label != "01234567")
      throw std::runtime_error("array section or fixed string width mismatch");
    for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 4; ++i) {
        const double expected_cell = j == 1 && i < 2 ? 11.0 + i : -7.0;
        if (owner(i, j) != expected_cell)
          throw std::runtime_error("namelist array section touched a guard column");
      }
    if (read("&test label(3:5)='X', object=2.5,.true. /") != 0)
      throw std::runtime_error(error);
    if (label != "01X  567" || component != 2.5 || !component_flag)
      throw std::runtime_error("character substring or aggregate object mismatch");

    Kokkos::deep_copy(values, 3.0);
    if (read("&test a=2*,1, label='abc   ' /") != 0 || values(0) != 3 || values(1) != 3 ||
        values(2) != 1 || label != "abc")
      throw std::runtime_error("null repetition or right-trim mismatch: " + error);
    if (read("&test i=1.5 /") != 0 || integer != 1 || read("&test i=-1.9 /") != 0 ||
        integer != -1 || read("&test i=3D1 /") != 0 || integer != 30)
      throw std::runtime_error("Native namelist integer truncation mismatch");
    if (read("&other i=2 /") >= 0 || integer != 30)
      throw std::runtime_error("missing group did not preserve initialized values");
    if (read("It's a model input\n&test i=3 /\n") != 0 || integer != 3 ||
        read("prefix &test i=4 /") >= 0 || read("& test i=4 /") >= 0)
      throw std::runtime_error("Native group search must ignore unrelated header text");
    const auto previous_label = label;
    if (read("&test label(5:4)='' &endjunk") != 0 || label != previous_label)
      throw std::runtime_error("Native empty substring or end marker mismatch");
    for (const std::string source :
         {"&test unknown=1 /", "&test a(7)=1 /", "&test i='1' /", "&test label=hello /",
          "&test a(6:2:-2)=3*4 /", "&test a=7*1 /", "&test label='unfinished"})
      if (read(source) <= 0 || error.empty())
        throw std::runtime_error("malformed namelist was accepted: " + source);

    std::ofstream(filename) << "&test i=8 /\n";
    namelist::Group padded_filename("test");
    padded_filename.bind("i", integer);
    if (padded_filename.read_file(filename.string() + "   ") != 0 || integer != 8)
      throw std::runtime_error("Native file-name padding was not ignored");

    std::filesystem::remove_all(directory);
    std::cout << "Namelist defaults, indices, sections, text and error checks passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    if (!directory.empty())
      std::filesystem::remove_all(directory);
    return 1;
  }
}
