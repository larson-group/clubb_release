// Diagnose and factor hydrometeor correlation matrices for SILHS.
//
// diagnose_correlations fills supported prescribed-w correlations; scalar
// calc_mean/calc_varnce/calc_w_corr retain device-callable mixture formulas.
// calc_cholesky_corr_mtx_approx, setup_corr_cholesky_mtx and
// cholesky_to_corr_mtx_approx construct the Native approximate factor/matrix.
// iiPDF_w remains one-based. Fields use first-index-fast LayoutLeft, with
// call-local scratch and ordered products/sums. Generic adapters retain caller
// sections through the existing array owner. True l_calc_w_corr is explicitly
// excluded, as in JAX, because the Native working matrix is uninitialized.
// See INTEGRATION.md for the precise scope and pending numerical evidence.
#pragma once
#include "array_types.hpp"
#include "constants_clubb.hpp"
#include "f2k/procedure_bridge.hpp"
#include <stdexcept>

namespace clubb_kokkos {
namespace diagnose_correlations_detail {
inline void diagnose_corr(int n_variables, Real1D sqrt_sigma2_on_mu2_ip,
                          Real2D corr_matrix_prescribed, Real2D corr_matrix_approx);
inline void corr_array_assertion_checks(int n_variables, Real2D corr_array);
inline void rearrange_corr_array(int pdf_dim, int iiPDF_w, Real2D corr_array,
                                 Real2D corr_array_swapped);
} // namespace diagnose_correlations_detail

//-----------------------------------------------------------------------
inline void diagnose_correlations(int pdf_dim, int iiPDF_w, // In
                                  Real2D corr_array_pre,    // In
                                  bool l_calc_w_corr,       // In
                                  Real2D corr_array)        // Out
{
  // Description:
  //   This subroutine diagnoses the correlation matrix in order to feed it
  //   into SILHS microphysics.

  // References:
  //   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
  //   (see CLUBB Trac ticket#514)
  //-----------------------------------------------------------------------

  // TODO: Native has no implementation that initializes corr_array_swapped
  // when this flag is true (its approx_w_corr routines are commented out).
  // Preserve the explicit JAX-supported exclusion instead of inventing data.
  if (l_calc_w_corr)
    throw std::runtime_error("diagnose_correlations: l_calc_w_corr=true is unsupported; "
                             "Native working matrix is uninitialized");
  Real2D corr_array_pre_swapped("diagnose_correlations::pre_swapped", pdf_dim, pdf_dim);
  Real2D corr_array_swapped("diagnose_correlations::swapped", pdf_dim, pdf_dim);
  // We actually don't need this right now: sigma_x^2/mu_x^2 (ith component), ip [-].
  Real1D sigma2_on_mu2_ip_array("diagnose_correlations::sigma2_on_mu2_ip", pdf_dim);
  // Initialize sigma2_on_mu2_ip_array.
  Kokkos::deep_copy(sigma2_on_mu2_ip_array, zero);

  // Swap the w-correlations to the first row for the prescribed correlations.
  diagnose_correlations_detail::rearrange_corr_array(pdf_dim, iiPDF_w,        // In
                                                     corr_array_pre,          // In
                                                     corr_array_pre_swapped); // Out
  // Diagnose correlations.
  if (!l_calc_w_corr)
    Kokkos::deep_copy(corr_array_swapped, corr_array_pre_swapped);
  // Materialize the Native whole-array SQRT argument without a host crossing.
  Real1D sqrt_sigma2_on_mu2_ip("diagnose_correlations::sqrt_sigma2_on_mu2_ip", pdf_dim);
  Kokkos::parallel_for(
      "diagnose_correlations::sqrt_argument", pdf_dim,
      KOKKOS_LAMBDA(int i) { sqrt_sigma2_on_mu2_ip(i) = Kokkos::sqrt(sigma2_on_mu2_ip_array(i)); });
  diagnose_correlations_detail::diagnose_corr(pdf_dim, sqrt_sigma2_on_mu2_ip, // In
                                              corr_array_pre_swapped,         // In
                                              corr_array_swapped);            // InOut
  // Swap rows back.
  diagnose_correlations_detail::rearrange_corr_array(pdf_dim, iiPDF_w,   // In
                                                     corr_array_swapped, // In
                                                     corr_array);        // Out
}

//-----------------------------------------------------------------------
namespace diagnose_correlations_detail {
inline void diagnose_corr(int n_variables, Real1D sqrt_sigma2_on_mu2_ip, // In
                          Real2D corr_matrix_prescribed,                 // In
                          Real2D corr_matrix_approx)                     // InOut
{
  // Description:
  //   This subroutine diagnoses the correlation matrix for each timestep.

  // References:
  //   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
  //   (see CLUBB Trac ticket#514)
  //-----------------------------------------------------------------------

  Real1D s_1j("diagnose_corr::s_1j", n_variables); // sqrt(1-c_1j^2) [-].
  // Native references this input only in its permanently disabled print block.
  (void)sqrt_sigma2_on_mu2_ip;
  // Calculate all square roots. Preserve Native SQRT without JAX gradient floors.
  Kokkos::parallel_for(
      "diagnose_corr::s_1j", n_variables, KOKKOS_LAMBDA(int level) {
        const int i = level + 1;
        s_1j(i - 1) = Kokkos::sqrt(real_t(1.) - corr_matrix_approx(i - 1, 1 - 1) *
                                                    corr_matrix_approx(i - 1, 1 - 1));
      });
  // Diagnose the missing correlations (Native calls this the upper triangle).
  Kokkos::parallel_for(
      "diagnose_corr::missing", n_variables * n_variables, KOKKOS_LAMBDA(int cell) {
        const int i = cell % n_variables + 1;
        const int j = cell / n_variables + 1;
        if (j >= 2 && j <= n_variables - 1 && i >= j + 1) {
          // Formula (16) in Larson et al. (2011), disabled in Native:
          // f_ij = alpha_corr * sqrt_sigma2_on_mu2_ip(i) * sqrt_sigma2_on_mu2_ip(j)
          //        * sign(1.0_core_rknd,corr_matrix_approx(1,i)*corr_matrix_approx(1,j)).
          // If the predicting c1i's are small then cij will be closer to the prescribed value. If
          // the c1i's are bigger, then cij will be closer to formula (15) from the reference paper.
          // See clubb:ticket:514:comment:61 for details. The alternative blend is also disabled:
          // f_ij = (1-abs(corr_matrix_approx(1,i)*corr_matrix_approx(1,j)))
          //          *corr_matrix_prescribed(i,j)
          //        +abs(corr_matrix_approx(1,i)*corr_matrix_approx(1,j))*f_ij_o.
          real_t f_ij = corr_matrix_prescribed(i - 1, j - 1);
          // Make sure -1 < f_ij < 1.
          if (f_ij < -max_mag_correlation)
            f_ij = -max_mag_correlation;
          else if (f_ij > max_mag_correlation)
            f_ij = max_mag_correlation;
          // Formula(15) in Larson et al.(2011).
          corr_matrix_approx(i - 1, j - 1) =
              corr_matrix_approx(i - 1, 1 - 1) * corr_matrix_approx(j - 1, 1 - 1) +
              f_ij * s_1j(i - 1) * s_1j(j - 1);
        }
      });
}
} // namespace diagnose_correlations_detail

//-----------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
real_t calc_w_corr(real_t wpxp, real_t stdev_w, real_t stdev_x, real_t w_tol, real_t x_tol) {
  // Description:
  // Compute the correlations of w with the hydrometeors.

  // References:
  // clubb:ticket:514
  //-----------------------------------------------------------------------

  // stdev_w/w_tol [m/s]; stdev_x/x_tol and covariance units follow x.
  real_t calc_w_corr = wpxp / (Kokkos::max(stdev_x, x_tol) * Kokkos::max(stdev_w, w_tol));
  // Make sure the correlation is in [-1,1].
  if (calc_w_corr < -max_mag_correlation)
    calc_w_corr = -max_mag_correlation;
  else if (calc_w_corr > max_mag_correlation)
    calc_w_corr = max_mag_correlation;
  return calc_w_corr;
}

//-----------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
real_t calc_varnce(real_t mixt_frac, real_t x1, real_t x2, real_t xm, real_t x1p2, real_t x2p2) {
  // Description:
  // Calculate the variance xp2 from the components x1, x2.

  // References:
  //   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02,
  //   page 3535
  //-----------------------------------------------------------------------

  // mixt_frac is dimensionless; x1/x2/xm carry field units, variances squared units.
  return mixt_frac * ((x1 - xm) * (x1 - xm) + x1p2) +
         (real_t(1.0) - mixt_frac) * ((x2 - xm) * (x2 - xm) + x2p2);
}

//-----------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
real_t calc_mean(real_t mixt_frac, real_t x1, real_t x2) {
  // Description:
  // Calculate the mean xm from the components x1, x2.

  // References:
  //   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02,
  //   page 3535
  //-----------------------------------------------------------------------

  // mixt_frac is dimensionless; component means and result carry field units.
  return mixt_frac * x1 + (real_t(1.0) - mixt_frac) * x2;
}

//-----------------------------------------------------------------------
inline void setup_corr_cholesky_mtx(int n_variables, Real2D corr_matrix,
                                    Real2D corr_cholesky_mtx_t);
inline void cholesky_to_corr_mtx_approx(int n_variables, Real2D corr_cholesky_mtx_t,
                                        Real2D corr_matrix_approx);
inline void calc_cholesky_corr_mtx_approx(int n_variables, int iiPDF_w,                     // In
                                          Real2D corr_matrix,                               // In
                                          Real2D corr_cholesky_mtx, Real2D corr_mtx_approx) // Out
{
  // Description:
  //   This subroutine calculates the transposed correlation cholesky matrix
  //   from the correlation matrix
  //
  // References:
  //   1 Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
  //   2 CLUBB Trac ticket#514
  //-----------------------------------------------------------------------

  // All matrices are dimensionless. iiPDF_w retains Native one-based metadata.
  Real2D corr_cholesky_mtx_swap("calc_cholesky_corr_mtx_approx::factor", n_variables, n_variables);
  Real2D corr_mtx_approx_swap("calc_cholesky_corr_mtx_approx::approx", n_variables, n_variables);
  Real2D corr_mtx_swap("calc_cholesky_corr_mtx_approx::input", n_variables, n_variables);
  diagnose_correlations_detail::rearrange_corr_array(n_variables, iiPDF_w,                 // In
                                                     corr_matrix,                          // In
                                                     corr_mtx_swap);                       // Out
  setup_corr_cholesky_mtx(n_variables, corr_mtx_swap,                                      // In
                          corr_cholesky_mtx_swap);                                         // Out
  diagnose_correlations_detail::rearrange_corr_array(n_variables, iiPDF_w,                 // In
                                                     corr_cholesky_mtx_swap,               // In
                                                     corr_cholesky_mtx);                   // Out
  cholesky_to_corr_mtx_approx(n_variables, corr_cholesky_mtx_swap,                         // In
                              corr_mtx_approx_swap);                                       // Out
  diagnose_correlations_detail::rearrange_corr_array(n_variables, iiPDF_w,                 // In
                                                     corr_mtx_approx_swap,                 // In
                                                     corr_mtx_approx);                     // Out
  diagnose_correlations_detail::corr_array_assertion_checks(n_variables, corr_mtx_approx); // In
  // Set Native (j,i) slots j<i to zero for conformity (actual upper triangle).
  Kokkos::parallel_for(
      "calc_cholesky_corr_mtx_approx::zero_triangle", n_variables * n_variables,
      KOKKOS_LAMBDA(int cell) {
        const int i = cell / n_variables + 1;
        const int j = cell % n_variables + 1;
        if (i >= 2 && j < i)
          corr_mtx_approx(j - 1, i - 1) = zero;
      });
}

//-----------------------------------------------------------------------
inline void setup_corr_cholesky_mtx(int n_variables, Real2D corr_matrix, // In
                                    Real2D corr_cholesky_mtx_t)          // Out
{
  // Description:
  //   This subroutine calculates the transposed correlation cholesky matrix
  //   from the correlation matrix
  //
  // References:
  //   1 Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
  //   2 CLUBB Trac ticket#514
  //-----------------------------------------------------------------------

  Real2D s("setup_corr_cholesky_mtx::s", n_variables, n_variables); // sqrt(1-c(i,j)^2) [-].
  // Calculate all necessary square roots.
  Kokkos::parallel_for(
      "setup_corr_cholesky_mtx::square_roots", n_variables * n_variables, KOKKOS_LAMBDA(int cell) {
        const int i = cell / n_variables + 1;
        const int j = cell % n_variables + 1;
        if (i <= n_variables - 1 && j >= i + 1)
          s(j - 1, i - 1) =
              Kokkos::sqrt(real_t(1.) - corr_matrix(j - 1, i - 1) * corr_matrix(j - 1, i - 1));
      });
  // Calculate transposed correlation Cholesky matrix; reference1 formula10.
  // Initialize matrix to zero.
  Kokkos::deep_copy(corr_cholesky_mtx_t, zero);
  // Initialize Native (j,i) slots j>=i and diagonal to one.
  Kokkos::parallel_for(
      "setup_corr_cholesky_mtx::ones", n_variables * n_variables, KOKKOS_LAMBDA(int cell) {
        const int i = cell / n_variables + 1;
        const int j = cell % n_variables + 1;
        if (j >= i)
          corr_cholesky_mtx_t(j - 1, i - 1) = one;
      });
  // Set diagonal elements. Each diagonal product retains ascending Native i order.
  Kokkos::parallel_for(
      "setup_corr_cholesky_mtx::diagonal", n_variables, KOKKOS_LAMBDA(int level) {
        const int j = level + 1;
        if (j >= 2) {
          for (int i = 1; i <= j - 1; ++i)
            corr_cholesky_mtx_t(j - 1, j - 1) = corr_cholesky_mtx_t(j - 1, j - 1) * s(j - 1, i - 1);
        }
      });
  // Set first row, stored at Native(j,1) as in the source.
  Kokkos::parallel_for(
      "setup_corr_cholesky_mtx::first_row", n_variables, KOKKOS_LAMBDA(int level) {
        const int j = level + 1;
        if (j >= 2)
          corr_cholesky_mtx_t(j - 1, 1 - 1) = corr_matrix(j - 1, 1 - 1);
      });
  // Set upper triangle (Native wording); preserve its actual(j,i),j>i slots.
  Kokkos::parallel_for(
      "setup_corr_cholesky_mtx::triangle", n_variables * n_variables, KOKKOS_LAMBDA(int cell) {
        const int i = cell / n_variables + 1;
        const int j = cell % n_variables + 1;
        if (i >= 2 && i <= n_variables - 1 && j >= i + 1) {
          for (int k = 1; k <= i - 1; ++k)
            corr_cholesky_mtx_t(j - 1, i - 1) = corr_cholesky_mtx_t(j - 1, i - 1) * s(j - 1, k - 1);
          corr_cholesky_mtx_t(j - 1, i - 1) =
              corr_cholesky_mtx_t(j - 1, i - 1) * corr_matrix(j - 1, i - 1);
        }
      });
}

//-----------------------------------------------------------------------
inline void cholesky_to_corr_mtx_approx(int n_variables, Real2D corr_cholesky_mtx_t, // In
                                        Real2D corr_matrix_approx)                   // Out
{
  // Description:
  //   This subroutine approximates the correlation matrix from the correlation
  //   cholesky matrix
  //
  // References:
  //   1 Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
  //   2 CLUBB Trac ticket#514
  //-----------------------------------------------------------------------

  // Approximate the correlation matrix; reference1 formula(8). Materialize
  // Native MATMUL(L,TRANSPOSE(L)) with one ordered contraction per output cell;
  // no parallel floating reduction. Intrinsic/compiler contraction is unverified.
  Kokkos::parallel_for(
      "cholesky_to_corr_mtx_approx", n_variables * n_variables, KOKKOS_LAMBDA(int cell) {
        const int i = cell % n_variables + 1;
        const int j = cell / n_variables + 1;
        real_t correlation_sum = real_t(0.);
        for (int k = 1; k <= n_variables; ++k)
          correlation_sum = correlation_sum +
                            corr_cholesky_mtx_t(i - 1, k - 1) * corr_cholesky_mtx_t(j - 1, k - 1);
        corr_matrix_approx(i - 1, j - 1) = correlation_sum;
      });
}

//-----------------------------------------------------------------------
namespace diagnose_correlations_detail {
inline void corr_array_assertion_checks(int n_variables, Real2D corr_array) // In
{
  // Description:
  //   This subroutine does the assertion checks for the corr_array.

  // References:
  //
  //
  //-----------------------------------------------------------------------

  constexpr real_t tol = real_t(1.e-6); // Maximum diagonal deviation from one.
  // Native ERROR STOP becomes the existing C++ fatal exception boundary.
  if (f2k::clubb_at_least_debug_level_api(1)) {
    const auto host = to_host_space(corr_array);
    for (int i = 1; i <= n_variables - 1; ++i) {
      for (int j = i + 1; j <= n_variables; ++j) {
        // Check upper and lower values against correlation boundaries.
        if (host(i - 1, j - 1) < -max_mag_correlation || host(i - 1, j - 1) > max_mag_correlation ||
            host(j - 1, i - 1) < -max_mag_correlation || host(j - 1, i - 1) > max_mag_correlation)
          throw std::runtime_error("Error: A value in the correlation matrix is out of range.");
      }
    }
  }
  if (f2k::clubb_at_least_debug_level_api(2)) {
    const auto host = to_host_space(corr_array);
    for (int i = 1; i <= n_variables; ++i) {
      // Check diagonal elements are one up to a tolerance.
      if (host(i - 1, i - 1) > one + tol || host(i - 1, i - 1) < one - tol)
        throw std::runtime_error(
            "Error: Diagonal element(s) of the correlation matrix are unequal to one.");
    }
  }
}

//-----------------------------------------------------------------------
inline void rearrange_corr_array(int pdf_dim, int iiPDF_w,  // In
                                 Real2D corr_array,         // In
                                 Real2D corr_array_swapped) // Out
{
  // Description:
  //   This subroutine swaps the w-correlations to the first row if the input
  //   matrix is in the same order as the *_corr_array_cloud.in files. It swaps
  //   the rows back to the order of the *_corr_array_cloud.in files if the
  //   input matrix is already swapped (first row w-correlations).
  //
  // References:
  //
  //-----------------------------------------------------------------------

  Real1D swap_array("rearrange_corr_array::swap_array", pdf_dim);
  // Swap w-correlations to the first row. Native metadata is one-based.
  Kokkos::deep_copy(corr_array_swapped, corr_array);
  Kokkos::deep_copy(swap_array, Kokkos::subview(corr_array_swapped, Kokkos::ALL, 0));
  // Keep each Native array assignment as a separate ordered stage. The first
  // reversed slice overlaps only its unchanged final cell; no cross-cell race.
  Kokkos::parallel_for(
      "rearrange_corr_array::first_reversed", iiPDF_w, KOKKOS_LAMBDA(int cell) {
        const int i = cell + 1;
        corr_array_swapped(i - 1, 1 - 1) = corr_array_swapped(iiPDF_w - 1, iiPDF_w - i + 1 - 1);
      });
  Kokkos::parallel_for(
      "rearrange_corr_array::first_tail", pdf_dim - iiPDF_w, KOKKOS_LAMBDA(int cell) {
        const int i = iiPDF_w + cell + 1;
        corr_array_swapped(i - 1, 1 - 1) = corr_array_swapped(i - 1, iiPDF_w - 1);
      });
  Kokkos::parallel_for(
      "rearrange_corr_array::w_reversed", iiPDF_w, KOKKOS_LAMBDA(int cell) {
        const int j = cell + 1;
        corr_array_swapped(iiPDF_w - 1, j - 1) = swap_array(iiPDF_w - j + 1 - 1);
      });
  Kokkos::parallel_for(
      "rearrange_corr_array::w_tail", pdf_dim - iiPDF_w, KOKKOS_LAMBDA(int cell) {
        const int i = iiPDF_w + cell + 1;
        corr_array_swapped(i - 1, iiPDF_w - 1) = swap_array(i - 1);
      });
}
} // namespace diagnose_correlations_detail

// Existing public array boundary retains host/strided caller sections.
template <class Input, class Output>
void diagnose_correlations(int pdf_dim, int iiPDF_w, Input corr_array_pre, bool l_calc_w_corr,
                           Output corr_array) {
  // Stage only this public call; private calls retain execution-space scratch.
  auto pre_exec = to_execution_space(corr_array_pre);
  auto out_exec = to_execution_space(corr_array, false);
  diagnose_correlations(pdf_dim, iiPDF_w, pre_exec, l_calc_w_corr, out_exec);
  from_execution_space(corr_array, out_exec);
}
template <class Input, class Factor, class Approximation>
void calc_cholesky_corr_mtx_approx(int n_variables, int iiPDF_w, Input corr_matrix,
                                   Factor corr_cholesky_mtx, Approximation corr_mtx_approx) {
  // Keep all output allocations owned by their original caller.
  auto matrix_exec = to_execution_space(corr_matrix);
  auto factor_exec = to_execution_space(corr_cholesky_mtx, false);
  auto approx_exec = to_execution_space(corr_mtx_approx, false);
  calc_cholesky_corr_mtx_approx(n_variables, iiPDF_w, matrix_exec, factor_exec, approx_exec);
  from_execution_space(corr_cholesky_mtx, factor_exec);
  from_execution_space(corr_mtx_approx, approx_exec);
}
template <class Input, class Output>
void setup_corr_cholesky_mtx(int n_variables, Input corr_matrix, Output corr_cholesky_mtx_t) {
  // Reuse the array owner rather than duplicating host/strided matrix conversion.
  auto matrix_exec = to_execution_space(corr_matrix);
  auto factor_exec = to_execution_space(corr_cholesky_mtx_t, false);
  setup_corr_cholesky_mtx(n_variables, matrix_exec, factor_exec);
  from_execution_space(corr_cholesky_mtx_t, factor_exec);
}
template <class Factor, class Output>
void cholesky_to_corr_mtx_approx(int n_variables, Factor corr_cholesky_mtx_t,
                                 Output corr_matrix_approx) {
  // Preserve the existing public array boundary and caller-owned output.
  auto factor_exec = to_execution_space(corr_cholesky_mtx_t);
  auto approx_exec = to_execution_space(corr_matrix_approx, false);
  cholesky_to_corr_mtx_approx(n_variables, factor_exec, approx_exec);
  from_execution_space(corr_matrix_approx, approx_exec);
}
} // namespace clubb_kokkos
