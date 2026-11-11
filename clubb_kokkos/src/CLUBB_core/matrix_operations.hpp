// Factor covariance matrices and mirror their lower triangles.
//
// Cholesky_factor uses a JAX-style column-wise factorization with Native
// scaling, ten total attempts and retry diagnostics. KOKKOS_BFB selects the
// original Fortran/f2k scaling, factorization and eigenvalue calls. The
// covariance-to-correlation and eigenvalue helpers support those diagnostics.
// mirror_lower_triangular_matrix copies lower to upper without changing the
// diagonal. set_lower_triangular_matrix writes one lower-triangle entry when
// both zero-based indices are nonnegative. Matrices retain first-index-fast
// LayoutLeft; temporary storage belongs to the active case workspace.
// Cholesky retains scoped host staging in both modes; its portable kernel
// also accepts device Views. Eigenvalues use a Jacobi diagnostic in custom mode.
// Other kernels retain their caller's memory space. Factorization diagnostics
// print locally in C++.
#pragma once
#include "error_code.hpp"
#include "array_types.hpp"
#include "constants_clubb.hpp"
#ifdef KOKKOS_BFB
#include "lapack_interfaces.hpp"
#endif
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace clubb_kokkos {
namespace matrix_operations_detail {
template <class Matrix> KOKKOS_INLINE_FUNCTION int portable_cholesky(int ndim, Matrix matrix) {
  // Factor the lower triangle in place using JAX's portable column algorithm.
  // Upper entries are preserved: Native callers may supply only the lower
  // triangle, so averaging with the upper triangle would change the contract.
  // Return the first non-positive/nonfinite pivot for the existing retry path.
  for (int column = 0; column < ndim; ++column) {
    real_t product = real_t(0);
    for (int p = 0; p < column; ++p)
      product += matrix(column, p) * matrix(column, p);

    const real_t diagonal = matrix(column, column) - product;
    if (!(diagonal > real_t(0)) || !Kokkos::isfinite(diagonal))
      return (column + 1);

    matrix(column, column) = Kokkos::sqrt(diagonal);
    for (int row = (column + 1); row < ndim; ++row) {
      product = real_t(0);
      for (int p = 0; p < column; ++p)
        product += matrix(row, p) * matrix(column, p);

      matrix(row, column) = (matrix(row, column) - product) / matrix(column, column);
    }
  }
  return 0;
}

template <class Input, class Eigenvalues>
void Symm_matrix_eigenvalues(ScratchWorkspace& workspace, int ndim, Input a_input,
                             Eigenvalues a_eigenvalues);
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
      "symm_covar_matrix_2_corr_matrix", Policy(0, (ndim * ndim)), KOKKOS_LAMBDA(int cell) {
        const int i = (cell % ndim);
        const int j = (cell / ndim);
        if (j < i + 1)
          corr(i, j) = covar(i, j) / Kokkos::sqrt(covar(i, i) * covar(j, j));
      });
}

//-------------------------------------------------------------------------------
template <class Input, class Scaling, class Factor>
void Cholesky_factor(ScratchWorkspace& workspace,                          // InOut
                     int ndim, Input a_input,                              // In
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

  auto workspace_scope    = workspace.scope();
  constexpr int itermax   = 10;          // Maximum total factorization attempts.
  constexpr real_t d_coef = real_t(0.1); // Failed-decomposition diagonal coefficient.
  using HostMatrix        = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  using HostVector        = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  // Scoped host staging preserves the current retry/diagnostic coordinator.
  // Reuse buffers across matrices; failure-only diagnostics reserve storage lazily.
  const auto a_input_host    = to_host_space(workspace, a_input);
  HostVector a_scaling_host  = workspace_scope.view<HostVector>("Cholesky_factor::a_scaling", ndim);
  HostMatrix a_Cholesky_host = workspace_scope.view<HostMatrix>(
      Kokkos::WithoutInitializing, "Cholesky_factor::a_Cholesky", ndim, ndim);
  HostVector a_eigenvalues;
  HostMatrix a_corr;
  HostMatrix a_scaled;
  real_t tau, d_smallest, amax, scond;
  int info;
  char equed;

  // Scale the factor input directly; retain a scaled retry copy only when used.
  Kokkos::deep_copy(a_Cholesky_host, a_input_host);
  equed = 'N';

#ifdef KOKKOS_BFB
  // Compute scaling for a_input, using Lapack routine spoequ for single precision,
  // or dpoequ for double precision.
  lapack_poequ(ndim, a_input_host.data(), ndim,           // In
               a_scaling_host.data(), scond, amax, info); // Out
  if (info == 0) {
    // Apply scaling to a_input, using Lapack routine slaqsy for single precision,
    // or dlaqsy for double precision.
    lapack_laqsy("Lower", ndim, a_Cholesky_host.data(), ndim, a_scaling_host.data(), scond, amax,
                 equed);
  }

#else
  // Compute Native diagonal scaling and apply it only when necessary.
  info  = 0;
  scond = real_t(1);
  amax  = real_t(0);
  if (ndim > 0) {
    real_t smallest = a_input_host(0, 0);
    amax            = smallest;
    for (int i = 0; i < ndim; ++i) {
      a_scaling_host(i) = a_input_host(i, i);
      smallest          = Kokkos::min(smallest, a_scaling_host(i));
      amax              = Kokkos::max(amax, a_scaling_host(i));
      if (info == 0 && !(a_scaling_host(i) > real_t(0)))
        info = i + 1;
    }

    if (info == 0) {
      for (int i = 0; i < ndim; ++i)
        a_scaling_host(i) = real_t(1) / Kokkos::sqrt(a_scaling_host(i));

      scond                      = Kokkos::sqrt(smallest) / Kokkos::sqrt(amax);
      constexpr real_t threshold = real_t(0.1); // Native LAQSY equilibration threshold.
      constexpr real_t small =
          std::numeric_limits<real_t>::min() / std::numeric_limits<real_t>::epsilon();
      constexpr real_t large = real_t(1) / small;
      if (scond < threshold || amax < small || amax > large) {
        equed = 'Y';
        for (int j = 0; j < ndim; ++j)
          for (int i = j; i < ndim; ++i)
            a_Cholesky_host(i, j) = a_scaling_host(j) * a_scaling_host(i) * a_Cholesky_host(i, j);
      }
    }
  }
#endif

  // Determine if scaling was necessary.
  l_scaled = equed == 'Y';
  if (l_scaled) {
    a_scaled = workspace_scope.view<HostMatrix>(Kokkos::WithoutInitializing,
                                                "Cholesky_factor::a_scaled", ndim, ndim);
    Kokkos::deep_copy(a_scaled, a_Cholesky_host);
  }

  for (int iter = 1; iter < itermax + 1; ++iter) {
#ifdef KOKKOS_BFB
    // Retain the original library's arithmetic for strict CPU comparisons.
    lapack_potrf("Lower", ndim, a_Cholesky_host.data(), ndim, info);
#else
    info = matrix_operations_detail::portable_cholesky(ndim, a_Cholesky_host);
#endif
    if (info < 0) {
      std::cerr << "Cholesky_factor "
                   " illegal value for argument "
                << (-1 * info) << '\n';
      // Native ERROR STOP terminates this path; throw at the existing C++ fatal boundary.
      throw std::runtime_error("Cholesky_factor: illegal LAPACK argument");
    } else if (info == 0) {
      // Success!
      if (clubb_at_least_debug_level_api(1) && iter > 1) {
        std::ostringstream message;
        message << "a_factored (worked)=\n" << std::setprecision(3);
        for (int i = 0; i < ndim; ++i) {
          for (int j = 0; j < i + 1; ++j)
            message << std::setw(10) << a_Cholesky_host(i, j);
          message << '\n';
        }
        std::cerr << message.str();
      }
      break;
    } else {
      if (clubb_at_least_debug_level_api(1)) {
        // This shouldn't happen now that the s and t Mellor(chi/eta) elements have been
        // modified to never be perfectly correlated, but it's here just in case.
        // -dschanen 10 Sept 2010
        std::ostringstream message;
        message << "Cholesky_factor: leading minor of order " << info
                << " is not positive definite.\nfactorization failed.\na_input=\n"
                << std::setprecision(3);
        for (int i = 0; i < ndim; ++i) {
          for (int j = 0; j < i + 1; ++j)
            message << std::setw(10) << a_input_host(i, j);
          message << '\n';
        }
        message << "a_Cholesky=\n";
        for (int i = 0; i < ndim; ++i) {
          for (int j = 0; j < i + 1; ++j)
            message << std::setw(10) << a_Cholesky_host(i, j);
          message << '\n';
        }
        std::cerr << message.str();
      }

      if (clubb_at_least_debug_level_api(2)) {
        if (!a_eigenvalues.is_allocated())
          a_eigenvalues = workspace_scope.view<HostVector>(Kokkos::WithoutInitializing,
                                                           "Cholesky_factor::a_eigenvalues", ndim);
        if (!a_corr.is_allocated())
          a_corr = workspace_scope.view<HostMatrix>(Kokkos::WithoutInitializing,
                                                    "Cholesky_factor::a_corr", ndim, ndim);
        matrix_operations_detail::Symm_matrix_eigenvalues(workspace, ndim, a_input_host, // In
                                                          a_eigenvalues);                // Out
        std::ostringstream message;
        message << "a_eigenvalues=\n" << std::setprecision(3);
        for (int i = 0; i < ndim; ++i)
          message << std::setw(10) << a_eigenvalues(i);
        message << '\n';
        symm_covar_matrix_2_corr_matrix(ndim, a_input_host, // In
                                        a_corr);            // Out
        message << "a_correlations=\n";
        for (int i = 0; i < ndim; ++i) {
          for (int j = 0; j < i + 1; ++j)
            message << std::setw(10) << a_corr(i, j);
          message << '\n';
        }
        std::cerr << message.str();
      }

      if (iter == itermax) {
        std::cerr << "iteration =" << iter << "itermax =" << itermax << '\n';
        std::cerr << "Fatal error in Cholesky_factor\n";
      } else if (clubb_at_least_debug_level_api(1)) {
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
      d_smallest = a_Cholesky_host(0, 0);
      for (int i = 1; i < ndim; ++i) {
        if (d_smallest > a_Cholesky_host(i, i))
          d_smallest = a_Cholesky_host(i, i);
      }
      // Use the smallest element * d_coef * iteration.
      tau = d_smallest * d_coef * real_t(iter);
      for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) {
          if (i == j) {
            a_Cholesky_host(i, j) = a_Cholesky_host(i, j) + tau;
          } else {
            a_Cholesky_host(i, j) = a_Cholesky_host(i, j);
          }
        }
      }

      if (clubb_at_least_debug_level_api(2)) {
        matrix_operations_detail::Symm_matrix_eigenvalues(workspace, ndim, a_Cholesky_host, // In
                                                          a_eigenvalues);                   // Out
        std::ostringstream message;
        message << "a_modified eigenvalues=\n" << std::scientific << std::setprecision(3);
        for (int i = 0; i < ndim; ++i)
          message << std::setw(10) << a_eigenvalues(i);
        message << '\n';
        std::cerr << message.str();
      }
    }
  }
  from_execution_space(workspace, a_scaling, a_scaling_host);
  from_execution_space(workspace, a_Cholesky, a_Cholesky_host);
}

//----------------------------------------------------------------------
namespace matrix_operations_detail {
template <class Input, class Eigenvalues>
void Symm_matrix_eigenvalues(ScratchWorkspace& workspace, // InOut
                             int ndim, Input a_input,     // In
                             Eigenvalues a_eigenvalues)   // Out
{
  //   Description:
  //     Compute ascending eigenvalues from the lower triangle of a_input.
  //     The input is preserved; this routine supplies failure diagnostics.
  //
  //   References:
  //     None
  //-----------------------------------------------------------------------

  using HostVector     = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  auto workspace_scope = workspace.scope();
  auto a_scratch = to_host_space(workspace, a_input); // Native value copy, never an input alias.
  HostVector a_eigenvalues_host =
      workspace_scope.view<HostVector>("Symm_matrix_eigenvalues::output", ndim);
  int info;
#ifdef KOKKOS_BFB
  constexpr int lwork = 180; // Native optimum obtained for n=5,31 Aug 2010.
  HostVector work     = workspace_scope.view<HostVector>("Symm_matrix_eigenvalues::work", lwork);
  // Lapack routine for eigenvalues and, optionally, eigenvectors: ssyev for
  // single precision or dsyev for double precision.
  lapack_syev("No eigenvectors", "Lower", ndim, a_scratch.data(), ndim, a_eigenvalues_host.data(),
              work.data(), lwork, info);
#else
  // Cyclic Jacobi rotations preserve the failure-only eigenvalue diagnostic
  // without a LAPACK dependency. Normalize first to avoid squared large values.
  real_t scale = real_t(0);
  info         = 0;
  for (int j = 0; j < ndim; ++j)
    for (int i = j; i < ndim; ++i) {
      if (!std::isfinite(a_scratch(i, j)))
        info = 1;
      scale = Kokkos::max(scale, Kokkos::abs(a_scratch(i, j)));
    }

  if (info == 0 && scale > real_t(0)) {
    for (int j = 0; j < ndim; ++j)
      for (int i = j; i < ndim; ++i)
        a_scratch(j, i) = a_scratch(i, j) = a_scratch(i, j) / scale;

    constexpr int max_sweeps = 50; // Complete cyclic sweeps of the small PDF matrix.
    info                     = 1;
    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
      real_t diagonal_max = real_t(0), off_diagonal_max = real_t(0);
      for (int p = 0; p < ndim; ++p) {
        diagonal_max = Kokkos::max(diagonal_max, Kokkos::abs(a_scratch(p, p)));
        for (int q = p + 1; q < ndim; ++q)
          off_diagonal_max = Kokkos::max(off_diagonal_max, Kokkos::abs(a_scratch(q, p)));
      }
      const real_t tolerance =
          std::numeric_limits<real_t>::epsilon() * Kokkos::max(real_t(1), diagonal_max);
      if (off_diagonal_max <= tolerance) {
        info = 0;
        break;
      }

      for (int p = 0; p < ndim; ++p)
        for (int q = p + 1; q < ndim; ++q) {
          const real_t entry = a_scratch(q, p);
          if (Kokkos::abs(entry) <= tolerance)
            continue;

          const real_t tau = (a_scratch(q, q) - a_scratch(p, p)) / (real_t(2) * entry);
          const real_t tangent =
              std::copysign(real_t(1), tau) / (Kokkos::abs(tau) + std::hypot(real_t(1), tau));
          const real_t cosine = real_t(1) / Kokkos::sqrt(real_t(1) + tangent * tangent);
          const real_t sine   = tangent * cosine;
          a_scratch(p, p) -= tangent * entry;
          a_scratch(q, q) += tangent * entry;
          a_scratch(q, p) = a_scratch(p, q) = real_t(0);
          for (int i = 0; i < ndim; ++i) {
            if (i == p || i == q)
              continue;

            const real_t ip = a_scratch(i, p), iq = a_scratch(i, q);
            a_scratch(i, p) = a_scratch(p, i) = cosine * ip - sine * iq;
            a_scratch(i, q) = a_scratch(q, i) = sine * ip + cosine * iq;
          }
        }
    }
  }
  for (int i = 0; i < ndim; ++i)
    a_eigenvalues_host(i) =
        info == 0 ? a_scratch(i, i) * scale : std::numeric_limits<real_t>::quiet_NaN();
  if (info == 0 && ndim > 1)
    std::sort(a_eigenvalues_host.data(), a_eigenvalues_host.data() + ndim);
#endif
  if (info < 0) {
    std::cerr << "Symm_matrix_eigenvalues:"
                 " illegal value for argument "
              << (-1 * info) << '\n';
  } else if (info == 0) {
    // Success!
  } else {
    std::cerr << "Symm_matrix_eigenvalues: Algorithm failed to converge.\n";
  }
  from_execution_space(workspace, a_eigenvalues, a_eigenvalues_host);
}
} // namespace matrix_operations_detail

//=============================================================================
template <class Matrix>
void set_lower_triangular_matrix(int pdf_dim, int index1, int index2, real_t xpyp, // In
                                 Matrix matrix)                                    // InOut
{
  // Description:
  // Set a value in the lower triangular portion of a square matrix. index1/2
  // are zero-based variate slots; either negative index skips it.
  // xpyp is usually a correlation or covariance, with the caller's units.
  (void)pdf_dim; // Number of variates; the owning View supplies its extent.

  // Reverse these to set values of the upper triangular matrix.
  const int i = Kokkos::max(index1, index2);
  const int j = Kokkos::min(index1, index2);
  if (i > -1 && j > -1) {
    Kokkos::parallel_for(
        "set_lower_triangular_matrix", Kokkos::RangePolicy<typename Matrix::execution_space>(0, 1),
        KOKKOS_LAMBDA(int) { matrix(i, j) = xpyp; });
  }
}

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
        "mirror_lower_triangular_matrix", Policy(0, (nvars * nvars)), KOKKOS_LAMBDA(int cell) {
          const int row = (cell % nvars);
          const int col = (cell / nvars);
          if (col >= 1 && row < col)
            matrix(row, col) = matrix(col, row);
        });
  }
}

} // namespace clubb_kokkos
