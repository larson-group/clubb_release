// -----------------------------------------------------------------------------
// Run complete CLUBB cases from aggregate namelists.
//
// Usage:
//   clubb_standalone [namelist_filename] runs once (default: clubb.in).
//   clubb_standalone -worker reads one namelist filename per stdin line until EOF.
//
// Understanding:
//   This standalone preserves the native one-shot and worker interface.
//   run_case owns err_info and its checks, as in the Fortran standalone.
//   clubb_driver.cpp owns the scientific stages and statistics in C++.
//   A worker acknowledges completion after those checks.
// -----------------------------------------------------------------------------

#include <Kokkos_Core.hpp>

#include <iostream>
#include <stdexcept>
#include <string>

#include "CLUBB_core/err_info_type_module.hpp"
#include "CLUBB_core/error_code.hpp"
#include "clubb_driver.hpp"
#include "driver_io.hpp"

namespace {

int worker_main();
void run_case(const std::string& namelist_filename);

} // namespace

int main(int argc, char* argv[]) {
  // Select one case or worker input, initialize Kokkos for the returned Views,
  // and return CLUBB's success code after all requested cases finish.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);

    // Run one namelist or read filenames from stdin, matching the native entry.
    std::string arg = "clubb.in";
    if (argc >= 2 && argv[1][0] != '\0') {
      arg = argv[1];
    }

    if (arg == "-worker") {
      return worker_main();
    }

    run_case(arg);
    return 6;
  } catch (const std::exception& error) {
    std::cerr << "ERROR STOP " << error.what() << '\n';
    return 1;
  }
}

namespace {

// -----------------------------------------------------------------------------
int worker_main() {
  // Read aggregate namelist paths until EOF. Report each completion after the
  // C++ driver closes output and the standalone checks and releases its error status.
  std::cout << "CLUBB_WORKER_READY" << std::endl;

  std::string namelist_filename;
  while (std::getline(std::cin, namelist_filename)) {
    // Fortran's fixed-length read and trim discard trailing blanks.
    while (!namelist_filename.empty() &&
           (namelist_filename.back() == ' ' || namelist_filename.back() == '\r')) {
      namelist_filename.pop_back();
    }
    if (namelist_filename.empty()) {
      std::cerr << "Worker namelist filename is empty\n";
      return 1;
    }

    run_case(namelist_filename);
    std::cout << "CLUBB_WORKER_DONE" << std::endl;
  }

  if (std::cin.bad()) {
    std::cerr << "Unable to read worker namelist filename\n";
    return 1;
  }

  return 6;
}

// -----------------------------------------------------------------------------
void run_case(const std::string& namelist_filename) {
  // Run one aggregate namelist through the C++ driver and check its returned
  // error codes before releasing the caller-owned arrays.
  clubb_kokkos::err_info_type err_info;

  // The driver owns the model lifecycle; this scope owns its output status.
  clubb_kokkos::ClubbDriver driver;
  driver.run_clubb(namelist_filename, true, err_info);
  if (!err_info.err_code.is_allocated()) {
    throw std::runtime_error("CLUBB driver returned no error codes");
  }
  // Keep the native standalone's fatal and generalized-grid checks in order.
  if (!clubb_kokkos::l_test_grid_generalization) {
    if (clubb_kokkos::any_equal(err_info.err_code, clubb_kokkos::clubb_fatal_error)) {
      throw std::runtime_error("Fatal error in clubb, check your parameter values and timestep");
    }

    std::cerr << " Program exited normally\n";
  } else {
    if (clubb_kokkos::any_equal(err_info.err_code, clubb_kokkos::clubb_generalized_grd_test_err)) {
      throw std::runtime_error("Error in generalized grid test; check the error messages");
    }

    std::cerr << " Generalized grid test passed\n";
  }

  clubb_kokkos::cleanup_err_info_api(err_info);
}

} // namespace
