// Verify fatal LAPACK stubs preserve inputs and clear uncomputed outputs.
// Both expert entry points must report the unsupported feature and mark every
// column fatal, including columns that already carry a nonfatal status.
#include "CLUBB_core/lapack_wrap.hpp"
#include "CLUBB_core/error_code.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[]) {
  // Exercise the public band/tridiagonal interfaces with multiple columns and
  // RHSs; sentinel inputs expose any accidental factorization or solve.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    using namespace clubb_kokkos;
    constexpr int nc = 3, ndim = 4, nrhs = 2;
    for (bool use_band : {false, true}) {
      Real3D lhs("lhs", nc, use_band ? 5 : 3, ndim);
      Real3D rhs("rhs", nc, ndim, nrhs), soln("soln", nc, ndim, nrhs);
      Real1D rcond("rcond", nc);
      Kokkos::deep_copy(lhs, real_t(42));
      Kokkos::deep_copy(rhs, real_t(17));
      Kokkos::deep_copy(soln, real_t(9876.5));
      Kokkos::deep_copy(rcond, real_t(1234));
      err_info_type err_info;
      err_info.err_code = decltype(err_info.err_code)("codes", nc);
      Kokkos::deep_copy(err_info.err_code, 50);
      err_info.err_header_global = "LAPACK rejection test";
      std::ostringstream diagnostics;
      auto* previous = std::cerr.rdbuf(diagnostics.rdbuf());
      try {
        if (use_band)
          lapack_band_solvex("test system", 2, 2, ndim, nrhs, nc, // In
                             lhs, rhs, err_info,                  // InOut
                             soln, rcond);                        // Out
        else
          lapack_tridiag_solvex("test system", ndim, nrhs, nc, // In
                                lhs, rhs, err_info,            // InOut
                                soln, rcond);                  // Out
      } catch (...) {
        std::cerr.rdbuf(previous);
        throw;
      }
      std::cerr.rdbuf(previous);
      const std::string message = diagnostics.str();
      if (message.find("LAPACK is not supported in Kokkos.") == std::string::npos ||
          message.find(err_info.err_header_global) == std::string::npos ||
          message.find("test system") == std::string::npos)
        throw std::runtime_error("missing LAPACK fatal diagnostic or caller context");
      const auto codes_h = to_host_space(err_info.err_code);
      const auto lhs_h = to_host_space(lhs);
      const auto rhs_h = to_host_space(rhs);
      const auto soln_h = to_host_space(soln);
      const auto rcond_h = to_host_space(rcond);
      for (int i = 0; i < nc; ++i)
        if (codes_h(i) != clubb_fatal_error || rcond_h(i) != real_t(0))
          throw std::runtime_error("LAPACK stub did not return all-column fatal status");
      for (std::size_t cell = 0; cell < lhs_h.size(); ++cell)
        if (lhs_h.data()[cell] != real_t(42))
          throw std::runtime_error("LAPACK stub changed lhs");
      for (std::size_t cell = 0; cell < rhs_h.size(); ++cell)
        if (rhs_h.data()[cell] != real_t(17) || soln_h.data()[cell] != real_t(0))
          throw std::runtime_error("LAPACK stub changed rhs or returned a solution");
    }
    std::cout << "Both LAPACK entry points return fatal err_info without solving\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
