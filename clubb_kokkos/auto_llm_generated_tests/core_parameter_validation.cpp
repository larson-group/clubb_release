// Compare parameter-validation error ownership with the actual Native routine.
// Five columns cover a hard-bound failure, retained fatal/generalized statuses,
// a later-column C6 equality failure and Native's unordered NaN comparisons.
#include "CLUBB_core/parameters_tunable.hpp"

#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

extern "C" void native_check_parameters_edge_c(int ngrdcol, int nparams_in,
                                                const clubb_kokkos::real_t* parameters,
                                                int* codes);

int main(int argc, char* argv[]) {
  // Empty-file initialization supplies real permanent defaults. Faults are
  // column-local; comparison with Native must also preserve pre-existing codes.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    using namespace clubb_kokkos;
    using Host2D = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    constexpr int ngrdcol = 5;
    Host2D parameters("parameter_validation::parameters", ngrdcol, nparams);
    init_clubb_params_api(ngrdcol, 0, std::string{}, parameters);

    parameters(0, iC8) = real_t(-1); // Native's C8 lower hard bound is zero.
    parameters(3, iC6rt) = real_t(3); // Positive/in-bounds, unequal to default C6thl=2.
    parameters(4, iC8) = std::numeric_limits<real_t>::quiet_NaN();
    parameters(4, iC6rt) = std::numeric_limits<real_t>::quiet_NaN();
    const std::array<int, ngrdcol> initial = {
        clubb_no_error, clubb_fatal_error, clubb_generalized_grd_test_err,
        clubb_no_error, clubb_no_error};
    const std::array<int, ngrdcol> expected = {
        clubb_fatal_error, clubb_fatal_error, clubb_generalized_grd_test_err,
        clubb_fatal_error, clubb_no_error};
    auto native_codes = initial;
    native_check_parameters_edge_c(ngrdcol, nparams, parameters.data(), native_codes.data());

    err_info_type err_info{};
    err_info.err_code = decltype(err_info.err_code)("parameter_validation::codes", ngrdcol);
    err_info.err_header =
        decltype(err_info.err_header)("parameter_validation::headers", 200, ngrdcol);
    auto codes = Kokkos::create_mirror_view(err_info.err_code);
    auto headers = Kokkos::create_mirror_view(err_info.err_header);
    Kokkos::deep_copy(headers, ' ');
    for (int i = 0; i < ngrdcol; ++i) {
      codes(i) = initial[i];
      const std::string header = "parameter validation column " + std::to_string(i + 1);
      for (std::size_t j = 0; j < header.size(); ++j)
        headers(j, i) = header[j];
    }
    Kokkos::deep_copy(err_info.err_code, codes);
    Kokkos::deep_copy(err_info.err_header, headers);
    check_parameters_api(ngrdcol, parameters, err_info);

    codes = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, err_info.err_code);
    for (int i = 0; i < ngrdcol; ++i) {
      if (native_codes[i] != expected[i] || codes(i) != native_codes[i])
        throw std::runtime_error("parameter status mismatch in column " + std::to_string(i + 1) +
                                 ": expected " + std::to_string(expected[i]) + ", Native " +
                                 std::to_string(native_codes[i]) + ", C++ " +
                                 std::to_string(codes(i)));
    }
    std::cout << "Parameter validator preserves Native five-column statuses: 99 99 50 99 0\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
