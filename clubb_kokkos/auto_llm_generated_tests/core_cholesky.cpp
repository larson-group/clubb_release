// Check portable Cholesky, scaling, retry and failure-only eigenvalue diagnostics.
// Analytic small matrices exercise lower-triangle ownership and reconstruction;
// the portable kernel also runs on the configured Kokkos execution space.
#include "scratch_workspace_owner.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/matrix_operations.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {
using namespace clubb_kokkos;
using HostMatrix = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;

void check_factor(HostMatrix input, Real2D factor, Real1D scaling, bool scaled,
                  real_t shift = real_t(0)) {
  // Reconstruct the symmetric input (plus a retry's diagonal shift) from L*L^T.
  // Upper input entries deliberately contain NaNs and must never participate.
  const auto factor_h = to_host_space(factor);
  const auto scaling_h = to_host_space(scaling);
  const int ndim = input.extent(0);
  const real_t tolerance = real_t(128 * ndim) * std::numeric_limits<real_t>::epsilon();
  for (int j = 0; j < ndim; ++j)
    for (int i = j; i < ndim; ++i) {
      real_t reconstructed = real_t(0);
      for (int p = 0; p <= j; ++p)
        reconstructed += factor_h(i, p) * factor_h(j, p);
      const real_t expected = input(i, j) * (scaled ? scaling_h(i) * scaling_h(j) : real_t(1)) +
                              (i == j ? shift : real_t(0));
      if (!std::isfinite(reconstructed) ||
          Kokkos::abs(reconstructed - expected) > tolerance * Kokkos::max(real_t(1), Kokkos::abs(expected))) {
        std::ostringstream message;
        message << "Cholesky reconstruction failed: " << factor.label() << " n=" << ndim
                << " (" << i << ',' << j << ") got=" << reconstructed << " expected=" << expected;
        throw std::runtime_error(message.str());
      }
      if (i > j && !std::isnan(factor_h(j, i)))
        throw std::runtime_error("Cholesky changed unused upper triangle");
    }
}
} // namespace

int main(int argc, char* argv[]) {
  // Check both public wrapper modes against reconstruction invariants, then
  // exercise non-positive/nonfinite pivots inside an actual Kokkos kernel.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    clubb_kokkos::test::ScratchWorkspaceOwner scratch_owner;
    auto& workspace = scratch_owner.workspace;
    using namespace clubb_kokkos;
    set_clubb_debug_level_api(2);
    for (int ndim : {1, 4, 8}) {
      HostMatrix lower("known_factor", ndim, ndim), input("input", ndim, ndim);
      for (int j = 0; j < ndim; ++j)
        for (int i = j; i < ndim; ++i)
          lower(i, j) = i == j ? real_t(i + 2) : real_t(i + j + 1) / real_t(8);
      for (int j = 0; j < ndim; ++j)
        for (int i = 0; i < ndim; ++i) {
          input(i, j) = std::numeric_limits<real_t>::quiet_NaN();
          if (i >= j) {
            input(i, j) = real_t(0);
            for (int p = 0; p <= j; ++p)
              input(i, j) += lower(i, p) * lower(j, p);
          }
        }
      Real2D kernel_factor("portable_factor", ndim, ndim);
      Kokkos::deep_copy(kernel_factor, input);
      Kokkos::View<int*> info("factor_status", 1);
      Kokkos::parallel_for("portable_cholesky_test", 1, KOKKOS_LAMBDA(int) {
        info(0) = matrix_operations_detail::portable_cholesky(ndim, kernel_factor);
      });
      if (to_host_space(info)(0) != 0)
        throw std::runtime_error("portable factor rejected a positive definite matrix");
      Real1D scaling("scaling", ndim);
      check_factor(input, kernel_factor, scaling, false);
      Real2D factor("factor", ndim, ndim);
      HostMatrix before("original_input", ndim, ndim);
      Kokkos::deep_copy(before, input);
      bool scaled;
      Cholesky_factor(workspace, ndim, input, scaling, factor, scaled);
      check_factor(input, factor, scaling, scaled);
      for (int j = 0; j < ndim; ++j)
        for (int i = 0; i < ndim; ++i)
          if (std::isnan(before(i, j)) ? !std::isnan(input(i, j)) : before(i, j) != input(i, j))
            throw std::runtime_error("Cholesky changed its input");
    }

    HostMatrix input("scaling_and_retry_input", 2, 2);
    Real1D scaling("scaling", 2);
    Real2D factor("factor", 2, 2);
    input(0, 0) = real_t(1.e-6);
    input(1, 0) = real_t(0.001);
    input(0, 1) = std::numeric_limits<real_t>::quiet_NaN();
    input(1, 1) = real_t(4);
    bool scaled;
    Cholesky_factor(workspace, 2, input, scaling, factor, scaled);
    if (!scaled)
      throw std::runtime_error("widely differing diagonals were not scaled");
    check_factor(input, factor, scaling, scaled);

    input(0, 0) = input(1, 1) = real_t(1);
    input(1, 0) = real_t(1.05);
    Cholesky_factor(workspace, 2, input, scaling, factor, scaled);
    if (scaled)
      throw std::runtime_error("correlation matrix was unexpectedly scaled");
    check_factor(input, factor, scaling, scaled, real_t(0.1));
    Real1D eigenvalues("eigenvalues", 2);
    matrix_operations_detail::Symm_matrix_eigenvalues(workspace, 2, input, eigenvalues);
    const auto eigenvalues_h = to_host_space(eigenvalues);
    const real_t tolerance = real_t(128) * std::numeric_limits<real_t>::epsilon();
    if (Kokkos::abs(eigenvalues_h(0) - real_t(-0.05)) > tolerance ||
        Kokkos::abs(eigenvalues_h(1) - real_t(2.05)) > tolerance)
      throw std::runtime_error("symmetric eigenvalue diagnostic is incorrect");
    if (input(0, 0) != real_t(1) || input(1, 0) != real_t(1.05) || !std::isnan(input(0, 1)))
      throw std::runtime_error("retry or eigenvalue diagnostic changed its input");

    for (real_t correlation : {real_t(0), real_t(0.25), real_t(-0.75)}) {
      // Equicorrelation has one eigenvalue 1+2*r and two eigenvalues 1-r.
      // Include a repeated spectrum and an indefinite matrix without retries.
      HostMatrix symmetric("equicorrelation", 3, 3);
      for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i)
          symmetric(i, j) = i == j ? real_t(1) : i > j ? correlation :
                                              std::numeric_limits<real_t>::quiet_NaN();
      Real1D spectrum("spectrum", 3);
      matrix_operations_detail::Symm_matrix_eigenvalues(workspace, 3, symmetric, spectrum);
      const auto spectrum_h = to_host_space(spectrum);
      real_t expected[3] = {real_t(1) + real_t(2) * correlation,
                            real_t(1) - correlation, real_t(1) - correlation};
      std::sort(expected, expected + 3);
      for (int i = 0; i < 3; ++i)
        if (!std::isfinite(spectrum_h(i)) || Kokkos::abs(spectrum_h(i) - expected[i]) > tolerance)
          throw std::runtime_error("repeated/indefinite eigenvalue diagnostic failed");
    }

    for (real_t diagonal : {real_t(0), real_t(-1), std::numeric_limits<real_t>::infinity(),
                            std::numeric_limits<real_t>::quiet_NaN()}) {
      Real2D invalid("invalid_factor", 1, 1);
      Kokkos::deep_copy(invalid, diagonal);
      Kokkos::View<int*> info("invalid_status", 1);
      Kokkos::parallel_for("portable_cholesky_invalid_test", 1, KOKKOS_LAMBDA(int) {
        info(0) = matrix_operations_detail::portable_cholesky(1, invalid);
      });
      if (to_host_space(info)(0) != 1)
        throw std::runtime_error("portable factor accepted an invalid pivot");
    }
#ifdef KOKKOS_BFB
    std::cout << "Cholesky wrapper: Fortran/f2k; portable kernel and edge checks passed\n";
#else
    std::cout << "Cholesky wrapper: custom; portable kernel and edge checks passed\n";
#endif
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
