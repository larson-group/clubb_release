#include "clubb_driver.hpp"
#include "CLUBB_core/error_code.hpp"

#include <Kokkos_Core.hpp>

#include <iostream>
#include <stdexcept>

namespace {

void check_error(const clubb_kokkos::err_info_type& err_info) {
  // Fail after the call that reports a fatal error, before a later reset hides it.
  if (!err_info.err_code.is_allocated()) throw std::runtime_error("missing error codes");
  const auto codes =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, err_info.err_code);
  for (std::size_t column = 0; column < codes.extent(0); ++column) {
    if (codes(column) == clubb_kokkos::clubb_fatal_error) {
      throw std::runtime_error("fatal error in driver sequence");
    }
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  // Exercise parameter snapshots/overrides, reinitialization, a suppressed
  // advance, reset and successive windows. Compare with a normal two-step run.
  if (argc != 2) return 1;
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    clubb_kokkos::ClubbDriver driver;
    clubb_kokkos::err_info_type err_info;

    driver.init_clubb_case(argv[1], err_info);
    check_error(err_info);
    clubb_kokkos::DriverArray2 parameters, independent_parameters;
    driver.get_clubb_params_all(parameters);
    const auto first_parameter = parameters(0, 0);
    parameters(0, 0) += 1;
    driver.get_clubb_params_all(independent_parameters);
    if (independent_parameters(0, 0) != first_parameter)
      throw std::runtime_error("parameter getter returned an alias");
    parameters(0, 0) = first_parameter;
    {
      clubb_kokkos::stats_type snapshot;
      driver.get_stats(snapshot);
      if (!snapshot.enabled)
        throw std::runtime_error("driver stats snapshot is not initialized");
    }
    driver.clean_up_clubb(err_info);
    check_error(err_info);
    driver.init_clubb_case(argv[1], err_info, &parameters);
    check_error(err_info);

    int total_param_sets = 0;
    int batch_size = 0;
    driver.get_runtime_batch_config(total_param_sets, batch_size);
    if (total_param_sets != batch_size) {
      throw std::runtime_error("driver sequence requires one runtime batch");
    }

    driver.advance_clubb_to_end(false, err_info, true, 1, 1);
    check_error(err_info);
    driver.set_case_initial_conditions(err_info, &parameters);
    check_error(err_info);
    // Print completed windows for the model comparison's step-count check.
    driver.advance_clubb_to_end(true, err_info, false, std::nullopt, 1);
    check_error(err_info);
    driver.advance_clubb_to_end(true, err_info, false, 2, std::nullopt);
    check_error(err_info);
    driver.clean_up_clubb(err_info);
    check_error(err_info);
    clubb_kokkos::cleanup_err_info_api(err_info);
    return 6;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
