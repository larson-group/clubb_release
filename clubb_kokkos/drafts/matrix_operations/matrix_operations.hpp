// Factor covariance matrices and mirror their lower triangles.
//
// Cholesky_factor retains Native LAPACK scaling, ten total factorization
// attempts, retry modification and complete diagnostic paths. The source-used
// covariance-to-correlation and eigenvalue helpers support those diagnostics.
// mirror_lower_triangular_matrix copies lower to upper without changing the
// diagonal. Matrices retain first-index-fast LayoutLeft; all scratch is local.
// Cholesky stages contiguous host buffers at the existing LAPACK owner; other
// kernels retain their caller's memory space. See INTEGRATION.md for scope.
#pragma once
#include "array_types.hpp"
#include "constants_clubb.hpp"
#include "lapack_interfaces.hpp"
#include "lapack_dependencies.hpp" // Pending declarations in the existing LAPACK owner.
#include "f2k/procedure_bridge.hpp"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace clubb_kokkos {
namespace matrix_operations_detail {
template <class Input, class Eigenvalues>
void Symm_matrix_eigenvalues(int ndim, Input a_input, Eigenvalues a_eigenvalues);
} // namespace matrix_operations_detail

//-----------------------------------------------------------------------
template <class Covariance, class Correlation>
void symm_covar_matrix_2_corr_matrix(int ndim, Covariance covar, // In
                                     Correlation corr)           // Out
{
  // Description:
  //   Convert a matrix of covariances in to a matrix of correlations.
  //   This only does the computation the lower triangular portion of the
  //   matrix.
  // References:
  //   None
  //-----------------------------------------------------------------------

  // covar has covariance units; corr is dimensionless. Both Views must be
  // accessible from the output execution space, as at the Native matrix call.
  Kokkos::deep_copy(corr, real_t(0.)); // Initialize to0.
  using Policy = Kokkos::RangePolicy<typename Correlation::execution_space>;
  Kokkos::parallel_for(
      "symm_covar_matrix_2_corr_matrix", Policy(0, ndim * ndim), KOKKOS_LAMBDA(int cell) {
        const int i = cell % ndim + 1;
        const int j = cell / ndim + 1;
        if (j <= i)
          corr(i - 1, j - 1) =
              covar(i - 1, j - 1) / Kokkos::sqrt(covar(i - 1, i - 1) * covar(j - 1, j - 1));
      });
}

//-------------------------------------------------------------------------------
template <class Input, class Scaling, class Factor>
void Cholesky_factor(int ndim, Input a_input,                              // In
                     Scaling a_scaling, Factor a_Cholesky, bool& l_scaled) // Out
{
  //  Description:
  //    Create a Cholesky factorization of a_input.
  //    If the factorization fails we use a modified a_input matrix and attempt
  //    to factorize again.
  //
  //  References:
  //    <http://www.netlib.org/lapack/explore-html/a00868.html> dpotrf
  //    <http://www.netlib.org/lapack/explore-html/a00860.html> dpoequ
  //    <http://www.netlib.org/lapack/explore-html/a00753.html> dlaqsy
  //-------------------------------------------------------------------------------

  constexpr int itermax = 10;            // Maximum total factorization attempts.
  constexpr real_t d_coef = real_t(0.1); // Failed-decomposition diagonal coefficient.
  using HostMatrix = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  using HostVector = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  // Host copies are call-owned LAPACK buffers; no scientific wrapper is bridged.
  const auto a_input_host = to_host_space(a_input);
  HostVector a_scaling_host("Cholesky_factor::a_scaling", ndim);
  HostMatrix a_Cholesky_host("Cholesky_factor::a_Cholesky", ndim, ndim);
  HostVector a_eigenvalues("Cholesky_factor::a_eigenvalues", ndim);
  HostMatrix a_corr("Cholesky_factor::a_corr", ndim, ndim);
  HostMatrix a_scaled("Cholesky_factor::a_scaled", ndim, ndim);
  real_t tau, d_smallest, amax, scond;
  int info;
  char equed;

  Kokkos::deep_copy(a_scaled, a_input_host); // Copy input array into local scaled array.
  equed = 'N';

  // Compute scaling for a_input, using Lapack routine spoequ for single precision,
  // or dpoequ for double precision.
  lapack_poequ(ndim, a_input_host.data(), ndim,           // In
               a_scaling_host.data(), scond, amax, info); // Out
  if (info == 0) {
    // Apply scaling to a_input, using Lapack routine slaqsy for single precision,
    // or dlaqsy for double precision.
    lapack_laqsy("Lower", ndim, a_scaled.data(), ndim, a_scaling_host.data(), scond, amax, equed);
  }

  // Determine if scaling was necessary.
  if (equed == 'Y') {
    l_scaled = true;
    Kokkos::deep_copy(a_Cholesky_host, a_scaled);
  } else {
    l_scaled = false;
    Kokkos::deep_copy(a_Cholesky_host, a_input_host);
  }

  for (int iter = 1; iter <= itermax; ++iter) {
    // Lapack Cholesky factorization, spotrf for single or dpotrf for double precision.
    lapack_potrf("Lower", ndim, a_Cholesky_host.data(), ndim, info);
    if (info < 0) {
      std::cerr << "Cholesky_factor "
                   " illegal value for argument "
                << -info << '\n';
      // Native ERROR STOP terminates this path; throw at the existing C++ fatal boundary.
      throw std::runtime_error("Cholesky_factor: illegal LAPACK argument");
    } else if (info == 0) {
      // Success!
      if (f2k::clubb_at_least_debug_level_api(1) && iter > 1) {
        std::ostringstream message;
        message << "a_factored (worked)=\n" << std::setprecision(3);
        for (int i = 1; i <= ndim; ++i) {
          for (int j = 1; j <= i; ++j)
            message << std::setw(10) << a_Cholesky_host(i - 1, j - 1);
          message << '\n';
        }
        std::cerr << message.str();
      }
      break;
    } else {
      if (f2k::clubb_at_least_debug_level_api(1)) {
        // This shouldn't happen now that the s and t Mellor(chi/eta) elements have been
        // modified to never be perfectly correlated, but it's here just in case.
        // -dschanen 10 Sept 2010
        std::ostringstream message;
        message << "Cholesky_factor: leading minor of order " << info
                << " is not positive definite.\nfactorization failed.\na_input=\n"
                << std::setprecision(3);
        for (int i = 1; i <= ndim; ++i) {
          for (int j = 1; j <= i; ++j)
            message << std::setw(10) << a_input_host(i - 1, j - 1);
          message << '\n';
        }
        message << "a_Cholesky=\n";
        for (int i = 1; i <= ndim; ++i) {
          for (int j = 1; j <= i; ++j)
            message << std::setw(10) << a_Cholesky_host(i - 1, j - 1);
          message << '\n';
        }
        std::cerr << message.str();
      }

      if (f2k::clubb_at_least_debug_level_api(2)) {
        matrix_operations_detail::Symm_matrix_eigenvalues(ndim, a_input_host, // In
                                                          a_eigenvalues);     // Out
        std::ostringstream message;
        message << "a_eigenvalues=\n" << std::setprecision(3);
        for (int i = 1; i <= ndim; ++i)
          message << std::setw(10) << a_eigenvalues(i - 1);
        message << '\n';
        symm_covar_matrix_2_corr_matrix(ndim, a_input_host, // In
                                        a_corr);            // Out
        message << "a_correlations=\n";
        for (int i = 1; i <= ndim; ++i) {
          for (int j = 1; j <= i; ++j)
            message << std::setw(10) << a_corr(i - 1, j - 1);
          message << '\n';
        }
        std::cerr << message.str();
      }

      if (iter == itermax) {
        std::cerr << "iteration =" << iter << "itermax =" << itermax << '\n';
        std::cerr << "Fatal error in Cholesky_factor\n";
      } else if (f2k::clubb_at_least_debug_level_api(1)) {
        // Adding a STOP statement to prevent this problem from slipping under
        // the rug. Native here only prints, then continues its modification.
        std::cerr << "Fatal error in Cholesky_factor\n";
        std::cerr << "Attempting to modify matrix to allow factorization.\n";
      }

      if (l_scaled) {
        Kokkos::deep_copy(a_Cholesky_host, a_scaled);
      } else {
        Kokkos::deep_copy(a_Cholesky_host, a_input_host);
      }
      // The number used for tau here is case specific to the Sigma covariance
      // matrix in the latin hypercube code and is not at all general.
      // Tau should be a number that is small relative to the other diagonal
      // elements of the matrix to keep the error caused by modifying 'a' low.
      // -dschanen 30 Aug 2010
      d_smallest = a_Cholesky_host(1 - 1, 1 - 1);
      for (int i = 2; i <= ndim; ++i) {
        if (d_smallest > a_Cholesky_host(i - 1, i - 1))
          d_smallest = a_Cholesky_host(i - 1, i - 1);
      }
      // Use the smallest element * d_coef * iteration.
      tau = d_smallest * d_coef * real_t(iter);
      for (int i = 1; i <= ndim; ++i) {
        for (int j = 1; j <= ndim; ++j) {
          if (i == j) {
            a_Cholesky_host(i - 1, j - 1) = a_Cholesky_host(i - 1, j - 1) + tau;
          } else {
            a_Cholesky_host(i - 1, j - 1) = a_Cholesky_host(i - 1, j - 1);
          }
        }
      }

      if (f2k::clubb_at_least_debug_level_api(2)) {
        matrix_operations_detail::Symm_matrix_eigenvalues(ndim, a_Cholesky_host, // In
                                                          a_eigenvalues);        // Out
        std::ostringstream message;
        message << "a_modified eigenvalues=\n" << std::scientific << std::setprecision(3);
        for (int i = 1; i <= ndim; ++i)
          message << std::setw(10) << a_eigenvalues(i - 1);
        message << '\n';
        std::cerr << message.str();
      }
    }
  }
  from_execution_space(a_scaling, a_scaling_host);
  from_execution_space(a_Cholesky, a_Cholesky_host);
}

//----------------------------------------------------------------------
namespace matrix_operations_detail {
template <class Input, class Eigenvalues>
void Symm_matrix_eigenvalues(int ndim, Input a_input,   // In
                             Eigenvalues a_eigenvalues) // Out
{
  //   Description:
  //     Computes the eigevalues of a_input
  //
  //   References:
  //     None
  //-----------------------------------------------------------------------

  constexpr int lwork = 180; // Native optimum obtained for n=5,31 Aug 2010.
  using HostVector = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  auto a_scratch = to_host_space(a_input); // Native value copy, never an input alias.
  HostVector a_eigenvalues_host("Symm_matrix_eigenvalues::output", ndim);
  HostVector work("Symm_matrix_eigenvalues::work", lwork);
  int info;
  // Lapack routine for eigenvalues and, optionally, eigenvectors: ssyev for
  // single precision or dsyev for double precision.
  lapack_syev("No eigenvectors", "Lower", ndim, a_scratch.data(), ndim, a_eigenvalues_host.data(),
              work.data(), lwork, info);
  if (info < 0) {
    std::cerr << "Symm_matrix_eigenvalues:"
                 " illegal value for argument "
              << -info << '\n';
  } else if (info == 0) {
    // Success!
  } else {
    std::cerr << "Symm_matrix_eigenvalues: Algorithm failed to converge.\n";
  }
  from_execution_space(a_eigenvalues, a_eigenvalues_host);
}
} // namespace matrix_operations_detail

//-----------------------------------------------------------------------
template <class Matrix>
void mirror_lower_triangular_matrix(int nvars,     // In
                                    Matrix matrix) // InOut
{
  // Description:
  //   Mirrors the elements of a lower triangular matrix to the upper
  //   triangle so that it is symmetric.

  // References:
  //   None
  //-----------------------------------------------------------------------

  // nvars is the square-matrix extent; first index is contiguous. Upper writes
  // read only the unchanged lower triangle, so all elements are independent.
  if (nvars > 1) {
    using Policy = Kokkos::RangePolicy<typename Matrix::execution_space>;
    Kokkos::parallel_for(
        "mirror_lower_triangular_matrix", Policy(0, nvars * nvars), KOKKOS_LAMBDA(int cell) {
          const int row = cell % nvars + 1;
          const int col = cell / nvars + 1;
          if (col >= 2 && row < col)
            matrix(row - 1, col - 1) = matrix(col - 1, row - 1);
        });
  }
}

} // namespace clubb_kokkos
