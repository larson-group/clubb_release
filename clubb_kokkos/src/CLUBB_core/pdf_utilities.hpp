// Convert PDF component moments and correlations between normal/lognormal variables.
//
// Lognormal conversions: mean_L2N/mean_L2N_dp and stdev_L2N/stdev_L2N_dp give
// moments of ln(x); corr_NL2NN/corr_NL2NN_dp and corr_LL2NN/corr_LL2NN_dp map
// correlations into normal space; corr_NN2NL and corr_NN2LL perform the inverse.
// Binormal moments: compute_mean_binormal, compute_variance_binormal and
// calc_comp_corrs_binormal combine the two normal PDF components.
// Coordinate conversions: calc_corr_chi_x/calc_corr_eta_x and their inverses
// calc_corr_rt_x/calc_corr_thl_x map correlations between thermodynamic and PDF
// coordinates. calc_xp2 combines lognormal/zero component variances;
// smooth_corr_quotient regularizes a correlation denominator smoothly.
//
// This mirrors all 19 public routines in src/CLUBB_core/pdf_utilities.F90.
// Scalar elemental forms are device callable; native array procedures use
// execution-space Views indexed (column, level), with LayoutLeft. The three
// rank-two elemental overloads at the end serve existing PDF closure callers.
// Inputs are required and conformable, as in Native. Array overloads can fill
// caller-owned results; returning overloads allocate independent owned storage.
// Distinct _dp behavior is retained even when real_t is already double.
// Scientific references are retained beside the corresponding procedures.
#pragma once

#include "array_types.hpp"
#include "constants_clubb.hpp"
#include "advance_helper_module.hpp"
#include <Kokkos_NumericTraits.hpp>

namespace clubb_kokkos {

KOKKOS_INLINE_FUNCTION real_t smooth_corr_quotient(real_t numerator, real_t denominator,
                                                   real_t denom_thresh) {
  // Regularize a same-unit numerator/denominator into a dimensionless correlation.
  const real_t smoothing = Kokkos::min(real_t(min_max_smth_mag), denom_thresh);
  real_t tmp_denom =
      smooth_max(Kokkos::abs(numerator) / max_mag_correlation, denominator, smoothing);
  tmp_denom = smooth_max(tmp_denom, denom_thresh, smoothing);
  return numerator / tmp_denom;
}

// Forward declaration permits the native calc_comp_corrs_binormal call order.
inline void smooth_corr_quotient(int ngrdcol, int nz, Real2D numerator, Real2D denominator,
                                 real_t denom_thresh, Real2D quotient);

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t mean_L2N(real_t mu_x, real_t sigma2_on_mu2) {
  // Description:
  // For a lognormally-distributed variable x, this function finds the mean of
  // ln x (mu_x_n) for the ith component of the PDF, given the mean of x (mu_x)
  // and the variance of x (sigma_sqd_x) for the ith component of the PDF. The
  // value ln x is distributed normally when x is distributed lognormally.

  // References:
  //  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
  //    Marcel Dekker, 401 pp.
  //  -- App. B.
  //-----------------------------------------------------------------------

  // Input Variables
  // mu_x: Mean of x (ith PDF component)                     [-]
  // sigma2_on_mu2: Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)    [-]

  // Return Variable
  real_t mu_x_n;

  // Find the mean of ln x for the ith component of the PDF.
  // The max( mu_x / sqrt( 1 + sigma_x^2 / mu_x^2 ), tiny( mu_x ) ) statement
  // is used to prevent taking ln 0, which will produce a result of -infinity.
  // This would happen when mu_x is 0.  However, this code should not be
  // entered when mu_x has a value of 0.
  mu_x_n = Kokkos::log(
      Kokkos::max(mu_x / Kokkos::sqrt(one + sigma2_on_mu2), real_t(Kokkos::norm_min_v<real_t>)));

  return mu_x_n;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION double mean_L2N_dp(double mu_x, double sigma2_on_mu2) {
  // Description:
  // For a lognormally-distributed variable x, this function finds the mean of
  // ln x (mu_x_n) for the ith component of the PDF, given the mean of x (mu_x)
  // and the variance of x (sigma_sqd_x) for the ith component of the PDF. The
  // value ln x is distributed normally when x is distributed lognormally.
  // This function uses double precision variables.

  // References:
  //  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
  //    Marcel Dekker, 401 pp.
  //  -- App. B.
  //-----------------------------------------------------------------------

  // Input Variables
  // mu_x: Mean of x (ith PDF component)                     [-]
  // sigma2_on_mu2: Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)    [-]

  // Return Variable
  double mu_x_n;

  // Find the mean of ln x for the ith component of the PDF.
  mu_x_n = Kokkos::log(mu_x / Kokkos::sqrt(one_dp + sigma2_on_mu2));

  return mu_x_n;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t stdev_L2N(real_t sigma2_on_mu2) {
  // Description:
  // For a lognormally-distributed variable x, this function finds the standard
  // deviation of ln x (sigma_x_n) for the ith component of the PDF, given the
  // mean of x (mu_x) and the variance of x (sigma_sqd_x) for the ith component
  // of the PDF.  The value ln x is distributed normally when x is distributed
  // lognormally.

  // References:
  //  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
  //    Marcel Dekker, 401 pp.
  //  -- App. B.
  //-----------------------------------------------------------------------

  // Input Variables
  // sigma2_on_mu2: Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)    [-]

  // Return Variable
  real_t sigma_x_n;

  // Find the standard deviation of ln x for the ith component of the PDF.
  sigma_x_n = Kokkos::sqrt(Kokkos::log(one + sigma2_on_mu2));

  return sigma_x_n;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION double stdev_L2N_dp(double sigma2_on_mu2) {
  // Description:
  // For a lognormally-distributed variable x, this function finds the standard
  // deviation of ln x (sigma_x_n) for the ith component of the PDF, given the
  // mean of x (mu_x) and the variance of x (sigma_sqd_x) for the ith component
  // of the PDF.  The value ln x is distributed normally when x is distributed
  // lognormally.
  // This function uses double precision variables.

  // References:
  //  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
  //    Marcel Dekker, 401 pp.
  //  -- App. B.
  //-----------------------------------------------------------------------

  // Input Variables
  // sigma2_on_mu2: Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)    [-]

  // Return Variable
  double sigma_x_n;

  // Find the standard deviation of ln x for the ith component of the PDF.
  sigma_x_n = Kokkos::sqrt(Kokkos::log(one_dp + sigma2_on_mu2));

  return sigma_x_n;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t corr_NL2NN(real_t corr_x_y, real_t sigma_y_n,
                                         real_t y_sigma2_on_mu2) {
  // Description:
  // For a normally-distributed variable x and a lognormally-distributed
  // variable y, this function finds the correlation of x and ln y (corr_x_y_n)
  // for the ith component of the PDF, given the correlation of x and y
  // (corr_x_y) and the standard deviation of ln y (sigma_y_n) for the ith
  // component of the PDF.  The value ln y is distributed normally when y is
  // distributed lognormally.

  // References:
  //  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
  //    Marcel Dekker, 401 pp.
  //  -- Eq. B-1.
  //-----------------------------------------------------------------------

  // Input Variables
  // corr_x_y: Correlation of x and y (ith PDF component)       [-]
  // sigma_y_n: Standard deviation of ln y (ith PDF component)   [-]
  // y_sigma2_on_mu2: Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

  // Return Variable
  real_t corr_x_y_n;

  //-------------------------------------- Begin Code --------------------------------------

  // Find the correlation of x and ln y for the ith component of the PDF.
  // When sigma_y = 0 and mu_y > 0, y_sigma2_on_mu2 = 0.  This results in
  // sigma_y_n = 0.  The resulting corr_x_y_n is undefined.  However, the
  // divide-by-zero problem needs to be addressed in the code.
  if (sigma_y_n > zero) {
    corr_x_y_n = corr_x_y * Kokkos::sqrt(y_sigma2_on_mu2) / sigma_y_n;
  } else {
    // The value of sqrt( y_sigma2_on_mu2 ) / sigma_y_n can be rewritten as:
    // sqrt( y_sigma2_on_mu2 ) / sqrt( ln( 1 + y_sigma2_on_mu2 ) ).
    // This can be further rewritten as:
    // sqrt( y_sigma2_on_mu2 / ln( 1 + y_sigma2_on_mu2 ) ),
    // which has a limit of 1 as y_sigma2_on_mu2 approaches 0 from the right.
    // When sigma_y_n = 0, the value of corr_x_y_n is undefined, so set it
    // to corr_x_y.
    corr_x_y_n = corr_x_y;
  }

  // Clip the magnitude of the correlation of x and ln y in the ith PDF
  // component, just in case the correlation (ith PDF component) of x and y and
  // the standard deviation (ith PDF component) of ln y are inconsistent,
  // resulting in an unrealizable value for corr_x_y_n.
  if (corr_x_y_n > max_mag_correlation) {
    corr_x_y_n = max_mag_correlation;
  } else if (corr_x_y_n < -max_mag_correlation) {
    corr_x_y_n = -max_mag_correlation;
  }

  return corr_x_y_n;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION double corr_NL2NN_dp(double corr_x_y, double sigma_y_n,
                                            double y_sigma2_on_mu2) {
  // Description:
  // For a normally-distributed variable x and a lognormally-distributed
  // variable y, this function finds the correlation of x and ln y (corr_x_y_n)
  // for the ith component of the PDF, given the correlation of x and y
  // (corr_x_y) and the standard deviation of ln y (sigma_y_n) for the ith
  // component of the PDF.  The value ln y is distributed normally when y is
  // distributed lognormally.
  // This function uses double precision variables.

  // References:
  //  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
  //    Marcel Dekker, 401 pp.
  //  -- Eq. B-1.
  //-----------------------------------------------------------------------

  // Input Variables
  // corr_x_y: Correlation of x and y (ith PDF component)       [-]
  // sigma_y_n: Standard deviation of ln y (ith PDF component)   [-]
  // y_sigma2_on_mu2: Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

  // Return Variable
  double corr_x_y_n;

  //-------------------------------------- Begin Code --------------------------------------

  // Find the correlation of x and ln y for the ith component of the PDF.
  // When sigma_y = 0 and mu_y > 0, y_sigma2_on_mu2 = 0.  This results in
  // sigma_y_n = 0.  The resulting corr_x_y_n is undefined.  However, the
  // divide-by-zero problem needs to be addressed in the code.
  if (sigma_y_n > zero_dp) {
    corr_x_y_n = corr_x_y * Kokkos::sqrt(y_sigma2_on_mu2) / sigma_y_n;
  } else {
    // The value of sqrt( y_sigma2_on_mu2 ) / sigma_y_n can be rewritten as:
    // sqrt( y_sigma2_on_mu2 ) / sqrt( ln( 1 + y_sigma2_on_mu2 ) ).
    // This can be further rewritten as:
    // sqrt( y_sigma2_on_mu2 / ln( 1 + y_sigma2_on_mu2 ) ),
    // which has a limit of 1 as y_sigma2_on_mu2 approaches 0 from the right.
    // When sigma_y_n = 0, the value of corr_x_y_n is undefined, so set it
    // to corr_x_y.
    corr_x_y_n = corr_x_y;
  }

  // Clip the magnitude of the correlation of x and ln y in the ith PDF
  // component, just in case the correlation (ith PDF component) of x and y and
  // the standard deviation (ith PDF component) of ln y are inconsistent,
  // resulting in an unrealizable value for corr_x_y_n.
  if (corr_x_y_n > static_cast<double>(max_mag_correlation)) {
    corr_x_y_n = static_cast<double>(max_mag_correlation);
  } else if (corr_x_y_n < -static_cast<double>(max_mag_correlation)) {
    corr_x_y_n = -static_cast<double>(max_mag_correlation);
  }

  return corr_x_y_n;
}

// =============================================================================
inline void corr_NN2NL(int nz, int ngrdcol, Real2D corr_x_y_n, Real2D sigma_y_n,
                       Real2D y_sigma2_on_mu2, // In
                       Real2D corr_x_y) {      // Out
  // Description:
  // For a normally-distributed variable x and a lognormally-distributed
  // variable y, this function finds the correlation of x and y (corr_x_y) for
  // the ith component of the PDF, given the correlation of x and ln y
  // (corr_x_y_n) and the standard deviation of ln y (sigma_y_n) for the ith
  // component of the PDF.  The value ln y is distributed normally when y is
  // distributed lognormally.

  // References:
  //  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
  //    Marcel Dekker, 401 pp.
  //  -- Eq. B-1.
  //-----------------------------------------------------------------------

  // Input Variables
  // nz: Number of vertical levels
  // ngrdcol: Number of grid columns
  // corr_x_y_n: Correlation of x and ln y (ith PDF component)    [-]
  // sigma_y_n: Standard deviation of ln y (ith PDF component)   [-]
  // y_sigma2_on_mu2: Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]
  // corr_x_y: Correlation of x and y (ith PDF component) [-]

  //-------------------------------------- Begin Code --------------------------------------
  // Find the correlation of x and y for the ith component of the PDF.
  // When sigma_y=0 and mu_y>0, y_sigma2_on_mu2=0 and sigma_y_n=0.
  // The correlations are then undefined; preserve the native limiting value.
  Kokkos::parallel_for(
      "corr_NN2NL", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int j = cell % ngrdcol;
        const int k = cell / ngrdcol;
        if (sigma_y_n(j, k) > zero) {
          // This ratio is already positive here. The tiny floor prevents an
          // underflowed denominator from causing division by zero.
          corr_x_y(j, k) =
              corr_x_y_n(j, k) * sigma_y_n(j, k) /
              Kokkos::sqrt(Kokkos::max(y_sigma2_on_mu2(j, k), real_t(Kokkos::norm_min_v<real_t>)));
        } else {
          // sqrt(log(1+r)/r) tends to1 as r approaches0 from the right.
          corr_x_y(j, k) = corr_x_y_n(j, k);
        }
        // Preserve the native ordered upper/lower clipping checks.
        if (corr_x_y(j, k) > max_mag_correlation) {
          corr_x_y(j, k) = max_mag_correlation;
        } else if (corr_x_y(j, k) < -max_mag_correlation) {
          corr_x_y(j, k) = -max_mag_correlation;
        }
      });
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t corr_LL2NN(real_t corr_x_y, real_t sigma_x_n, real_t sigma_y_n, // In
                                         real_t x_sigma2_on_mu2, real_t y_sigma2_on_mu2) {
  // Description:
  // For lognormally-distributed variables x and y, this function finds the
  // correlation of ln x and ln y (corr_x_y_n) for the ith component of the
  // PDF, given the correlation of x and y (corr_x_y), the standard deviation
  // of ln x (sigma_x_n), and the standard deviation of ln y (sigma_y_n) for
  // the ith component of the PDF.  The value of ln x (or ln y) is distributed
  // normally when x (or y) is distributed lognormally.

  // References:
  //  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
  //    Marcel Dekker, 401 pp.
  //  -- Eq. C-3.
  //-----------------------------------------------------------------------

  // Input Variables
  // corr_x_y: Correlation of x and y (ith PDF component)       [-]
  // sigma_x_n: Standard deviation of ln x (ith PDF component)   [-]
  // sigma_y_n: Standard deviation of ln y (ith PDF component)   [-]
  // x_sigma2_on_mu2: Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)   [-]
  // y_sigma2_on_mu2: Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

  // Return Variable
  real_t corr_x_y_n;

  // Local Variable
  real_t log_arg; // Input into the natural logarithm [-]

  //-------------------------------------- Begin Code --------------------------------------

  // Find the correlation of ln x and ln y for the ith component of the PDF.
  // When sigma_x = 0 and mu_x > 0, x_sigma2_on_mu2 = 0.  This results in
  // sigma_x_n = 0.  The resulting corr_x_y_n is undefined.  The same holds
  // true when sigma_y = 0 and mu_y > 0.  However, the divide-by-zero problem
  // needs to be addressed in the code.
  if (sigma_x_n > zero && sigma_y_n > zero) {
    //      corr_x_y_n = log( one + corr_x_y * sqrt( exp( sigma_x_n**2 ) - one ) &
    //                                       * sqrt( exp( sigma_y_n**2 ) - one )  ) &
    //                   / ( sigma_x_n * sigma_y_n )
    log_arg    = one + corr_x_y * Kokkos::sqrt(x_sigma2_on_mu2 * y_sigma2_on_mu2);
    corr_x_y_n = Kokkos::log(log_arg) / (sigma_x_n * sigma_y_n);
  } else {
    // The value of corr_x_y_n is undefined, so set it to corr_x_y.
    corr_x_y_n = corr_x_y;
  }

  // Clip the magnitude of the correlation of ln x and ln y in the ith PDF
  // component, just in case the correlation (ith PDF component) of x and y,
  // the standard deviation (ith PDF component) of ln x, and the standard
  // deviation (ith PDF component) of ln y are inconsistent, resulting in an
  // unrealizable value for corr_x_y_n.
  if (corr_x_y_n > max_mag_correlation) {
    corr_x_y_n = max_mag_correlation;
  } else if (corr_x_y_n < -max_mag_correlation) {
    corr_x_y_n = -max_mag_correlation;
  }

  return corr_x_y_n;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION double corr_LL2NN_dp(double corr_x_y, double sigma_x_n,
                                            double sigma_y_n, // In
                                            double x_sigma2_on_mu2, double y_sigma2_on_mu2) {
  // Description:
  // For lognormally-distributed variables x and y, this function finds the
  // correlation of ln x and ln y (corr_x_y_n) for the ith component of the
  // PDF, given the correlation of x and y (corr_x_y), the standard deviation
  // of ln x (sigma_x_n), and the standard deviation of ln y (sigma_y_n) for
  // the ith component of the PDF.  The value of ln x (or ln y) is distributed
  // normally when x (or y) is distributed lognormally.
  // This function uses double precision variables.

  // References:
  //  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
  //    Marcel Dekker, 401 pp.
  //  -- Eq. C-3.
  //-----------------------------------------------------------------------

  // Input Variables
  // corr_x_y: Correlation of x and y (ith PDF component)       [-]
  // sigma_x_n: Standard deviation of ln x (ith PDF component)   [-]
  // sigma_y_n: Standard deviation of ln y (ith PDF component)   [-]
  // x_sigma2_on_mu2: Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)   [-]
  // y_sigma2_on_mu2: Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]

  // Return Variable
  double corr_x_y_n;

  //-------------------------------------- Begin Code --------------------------------------

  // Find the correlation of ln x and ln y for the ith component of the PDF.
  // When sigma_x = 0 and mu_x > 0, x_sigma2_on_mu2 = 0.  This results in
  // sigma_x_n = 0.  The resulting corr_x_y_n is undefined.  The same holds
  // true when sigma_y = 0 and mu_y > 0.  However, the divide-by-zero problem
  // needs to be addressed in the code.
  if (sigma_x_n > zero_dp && sigma_y_n > zero_dp) {
    corr_x_y_n = Kokkos::log(one_dp + corr_x_y * Kokkos::sqrt(x_sigma2_on_mu2 * y_sigma2_on_mu2)) /
                 (sigma_x_n * sigma_y_n);
  } else {
    // The value of corr_x_y_n is undefined, so set it to corr_x_y.
    corr_x_y_n = corr_x_y;
  }

  // Clip the magnitude of the correlation of ln x and ln y in the ith PDF
  // component, just in case the correlation (ith PDF component) of x and y,
  // the standard deviation (ith PDF component) of ln x, and the standard
  // deviation (ith PDF component) of ln y are inconsistent, resulting in an
  // unrealizable value for corr_x_y_n.
  if (corr_x_y_n > static_cast<double>(max_mag_correlation)) {
    corr_x_y_n = static_cast<double>(max_mag_correlation);
  } else if (corr_x_y_n < -static_cast<double>(max_mag_correlation)) {
    corr_x_y_n = -static_cast<double>(max_mag_correlation);
  }

  return corr_x_y_n;
}

// =============================================================================
inline void corr_NN2LL(int nz, int ngrdcol, Real2D corr_x_y_n, Real2D sigma_x_n,
                       Real2D sigma_y_n,                               // In
                       Real2D x_sigma2_on_mu2, Real2D y_sigma2_on_mu2, // In
                       Real2D corr_x_y) {                              // Out
  // Description:
  // For lognormally-distributed variables x and y, this function finds the
  // correlation of x and y (corr_x_y) for the ith component of the PDF, given
  // the correlation of ln x and ln y (corr_x_y_n), the standard deviation of
  // ln x (sigma_x_n), and the standard deviation of ln y (sigma_y_n) for
  // the ith component of the PDF.  The value of ln x (or ln y) is distributed
  // normally when x (or y) is distributed lognormally.

  // References:
  //  Garvey, P. R., 2000: Probability methods for cost uncertainty analysis.
  //    Marcel Dekker, 401 pp.
  //  -- Eq. C-3.
  //-----------------------------------------------------------------------

  // Input Variables
  // nz: Number of vertical levels
  // ngrdcol: Number of grid columns
  // corr_x_y_n: Correlation of ln x and ln y (ith PDF component) [-]
  // sigma_x_n: Standard deviation of ln x (ith PDF component)   [-]
  // sigma_y_n: Standard deviation of ln y (ith PDF component)   [-]
  // x_sigma2_on_mu2: Ratio:  sigma_x^2 / mu_x^2 (ith PDF component)   [-]
  // y_sigma2_on_mu2: Ratio:  sigma_y^2 / mu_y^2 (ith PDF component)   [-]
  // corr_x_y: Correlation of x and y (ith PDF component)  [-]

  //-------------------------------------- Begin Code --------------------------------------
  // Either zero transformed standard deviation makes both correlations
  // undefined. Return the supplied correlation in that native limiting case.
  Kokkos::parallel_for(
      "corr_NN2LL", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int j = cell % ngrdcol;
        const int k = cell / ngrdcol;
        if (sigma_x_n(j, k) > zero && sigma_y_n(j, k) > zero) {
          // Native alternate form:
          // (exp(sigma_x_n*sigma_y_n*corr_x_y_n)-one)
          //   /(sqrt(exp(sigma_x_n**2)-one)*sqrt(exp(sigma_y_n**2)-one)).
          corr_x_y(j, k) =
              (Kokkos::exp(sigma_x_n(j, k) * sigma_y_n(j, k) * corr_x_y_n(j, k)) - one) /
              Kokkos::sqrt(x_sigma2_on_mu2(j, k) * y_sigma2_on_mu2(j, k));
        } else {
          corr_x_y(j, k) = corr_x_y_n(j, k);
        }
        // Clip inconsistent component moments to a realizable correlation.
        if (corr_x_y(j, k) > max_mag_correlation) {
          corr_x_y(j, k) = max_mag_correlation;
        } else if (corr_x_y(j, k) < -max_mag_correlation) {
          corr_x_y(j, k) = -max_mag_correlation;
        }
      });
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t compute_mean_binormal(real_t mu_x_1, real_t mu_x_2,
                                                    real_t mixt_frac) {
  // Description:
  // Computes the overall grid-box mean of a binormal distribution from the
  // mean of each component

  // References:
  //   None
  //-----------------------------------------------------------------------

  // Input Variables
  // mu_x_1: First PDF component mean of 'x'                       [?]
  // mu_x_2: Second PDF component mean of 'x'                      [?]
  // mixt_frac: Weight of the first PDF component                     [-]

  // Return Variable
  real_t xm;

  //----- Begin Code -----
  xm = mixt_frac * mu_x_1 + (one - mixt_frac) * mu_x_2;

  return xm;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t compute_variance_binormal(real_t xm, real_t mu_x_1,
                                                        real_t mu_x_2,                      // In
                                                        real_t stdev_x_1, real_t stdev_x_2, // In
                                                        real_t mixt_frac) {
  // Description:
  // Computes the overall grid-box variance of a binormal distribution from the
  // variance of each component.

  // References:
  //   None
  //-----------------------------------------------------------------------

  // Input Variables
  // xm: Overall mean of 'x'                                   [?]
  // mu_x_1: First PDF component mean of 'x'                       [?]
  // mu_x_2: Second PDF component mean of 'x'                      [?]
  // stdev_x_1: Standard deviation of 'x' in the first PDF component  [?]
  // stdev_x_2: Standard deviation of 'x' in the second PDF component [?]
  // mixt_frac: Weight of the first PDF component                     [-]

  // Return Variable
  real_t xp2;

  //----- Begin Code -----
  xp2 = mixt_frac * (((mu_x_1 - xm) * (mu_x_1 - xm)) + (stdev_x_1 * stdev_x_1)) +
        (one - mixt_frac) * (((mu_x_2 - xm) * (mu_x_2 - xm)) + (stdev_x_2 * stdev_x_2));

  return xp2;
}

// =============================================================================
inline void calc_comp_corrs_binormal(int nz, int ngrdcol,               // In
                                     Real2D xpyp, Real2D xm, Real2D ym, // In
                                     Real2D mu_x_1, Real2D mu_x_2,      // In
                                     Real2D mu_y_1, Real2D mu_y_2,      // In
                                     Real2D sigma_x_1_sqd,              // In
                                     Real2D sigma_x_2_sqd,              // In
                                     Real2D sigma_y_1_sqd,              // In
                                     Real2D sigma_y_2_sqd,              // In
                                     Real2D mixt_frac,                  // In
                                     Real2D corr_x_y_1,                 // Out
                                     Real2D corr_x_y_2) {               // Out
  // Description:
  // Calculates the PDF component correlations of variables x and y, where
  // x and y are both distributed as two-component normals (or binormals).
  // The PDF component correlations are set equal to each other.
  //
  // The overall covariance of x and y, <x'y'>, can be expressed in terms of
  // PDF parameters by integrating over the PDF:
  //
  // <x'y'> = INT(-inf:inf) INT(-inf:inf) ( x - <x> ) ( y - <y> ) P(x,y) dy dx;
  //
  // where <x> is the overall mean of x, <y> is the overall mean of y, and
  // P(x,y) is the equation for the two-component normal PDF of x and y.
  //
  // The integral is evaluated, and the equation for <x'y'> is:
  //
  // <x'y'> = mixt_frac * ( ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
  //                        + corr_x_y_1 * sigma_x_1 * sigma_y_1 )
  //          + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
  //                                  + corr_x_y_2 * sigma_x_2 * sigma_y_2 );
  //
  // where mu_x_1 is the mean of x in the 1st PDF component, mu_x_2 is the mean
  // of x in the 2nd PDF component, mu_y_1 is the mean of y in the 1st PDF
  // component, mu_y_2 is the mean of y in the 2nd PDF component, sigma_x_1 is
  // the standard deviation of x in the 1st PDF component, sigma_x_2 is the
  // standard deviation of x in the 2nd PDF component, sigma_y_1 is the
  // standard deviation of y in the 1st PDF component, sigma_y_2 is the
  // standard deviation of y in the 2nd PDF component, corr_x_y_1 is the
  // correlation of x and y in the 1st PDF component, corr_x_y_2 is the
  // correlation of x and y in the 2nd PDF component, and mixt_frac is the
  // mixture fraction (weight of the 1st PDF component).
  //
  // This equation can be rewritten as:
  //
  // <x'y'> = mixt_frac * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
  //          + mixt_frac * corr_x_y_1 * sigma_x_1 * sigma_y_1
  //          + ( 1 - mixt_frac ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
  //          + ( 1 - mixt_frac ) * corr_x_y_2 * sigma_x_2 * sigma_y_2.
  //
  // Setting the two PDF component correlations equal to each other
  // (corr_x_y_1 = corr_x_y_2), the equation can be solved for the PDF
  // component correlations:
  //
  // corr_x_y_1 = corr_x_y_2
  // = ( <x'y'> - mixt_frac * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
  //            - ( 1 - mixt_frac ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> ) )
  //   / ( mixt_frac * sigma_x_1 * sigma_y_1
  //       + ( 1 - mixt_frac ) * sigma_x_2 * sigma_y_2 );
  //
  // where -1 <= corr_x_y_1 = corr_x_y_2 <= 1.
  //
  // When sigma_x_1 * sigma_y_1 = 0 and sigma_x_2 * sigma_y_2 = 0, at least one
  // of x or y are constant within each PDF component, and both PDF component
  // correlations are undefined.

  // References:
  //-----------------------------------------------------------------------

  // Input/Output Variables
  // nz: Number of vertical levels
  // ngrdcol: Number of grid columns
  // xpyp: Covariance of x and y (overall)    [(x units)(y units)]
  // xm: Mean of x (overall)                [x units]
  // ym: Mean of y (overall)                [y units]
  // mu_x_1: Mean of x (1st PDF component)      [x units]
  // mu_x_2: Mean of x (2nd PDF component)      [x units]
  // mu_y_1: Mean of y (1st PDF component)      [y units]
  // mu_y_2: Mean of y (2nd PDF component)      [y units]
  // sigma_x_1_sqd: Variance of x (1st PDF component)  [(x units)^2]
  // sigma_x_2_sqd: Variance of x (2nd PDF component)  [(x units)^2]
  // sigma_y_1_sqd: Variance of y (1st PDF component)  [(y units)^2]
  // sigma_y_2_sqd: Variance of y (2nd PDF component)  [(y units)^2]
  // mixt_frac: Mixture fraction                   [-]
  // corr_x_y_1: Correlation of x and y (1st PDF component)    [-]
  // corr_x_y_2: Correlation of x and y (2nd PDF component)    [-]

  // Both components use the same pointwise correlation. Keep its numerator,
  // denominator and regularization local instead of separate arrays/kernels.
  Kokkos::parallel_for(
      "calc_comp_corrs_binormal::components", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        const real_t numerator =
            xpyp(i, k) - mixt_frac(i, k) * (mu_x_1(i, k) - xm(i, k)) * (mu_y_1(i, k) - ym(i, k)) -
            (one - mixt_frac(i, k)) * (mu_x_2(i, k) - xm(i, k)) * (mu_y_2(i, k) - ym(i, k));
        const real_t denominator =
            mixt_frac(i, k) * Kokkos::sqrt(sigma_x_1_sqd(i, k) * sigma_y_1_sqd(i, k)) +
            (one - mixt_frac(i, k)) * Kokkos::sqrt(sigma_x_2_sqd(i, k) * sigma_y_2_sqd(i, k));
        const real_t correlation = smooth_corr_quotient(numerator, denominator, eps);
        corr_x_y_1(i, k)         = correlation;
        // Set corr_x_y_2 equal to corr_x_y_1.
        corr_x_y_2(i, k) = correlation;
      });
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t calc_corr_chi_x(real_t crt_i, real_t cthl_i,           // In
                                              real_t sigma_rt_i, real_t sigma_thl_i, // In
                                              real_t sigma_chi_i,                    // In
                                              real_t corr_rt_x_i, real_t corr_thl_x_i) {
  // Description:
  // This function calculates the correlation of extended liquid water mixing
  // ratio, chi (old s), and a generic variable x, within the ith component of
  // the PDF.  The variable chi can be split into mean and turbulent
  // components, such that:
  //
  // chi = <chi> + chi';
  //
  // where < > denotes a mean field an ' denotes a turbulent component.
  //
  // The linearized equation for chi' is given in Larson et al. (2001), where
  // within the ith component of the PDF:
  //
  // chi_(i)' = Coef_rt(i) * r_t(i)' - Coef_thl(i) * th_l(i)'.
  //
  // The equation for chi' can be multiplied by x'.  The equation becomes:
  //
  // chi'x'_(i) = Coef_rt(i) * r_t'x'_(i) - Coef_thl(i) * th_l'x'_(i).
  //
  // Averaging both sides, the covariance <chi'x'> is given by the equation:
  //
  // <chi'x'_(i)> = Coef_rt(i) * <r_t'x'_(i)> - Coef_thl(i) * <th_l'x'_(i)>.
  //
  // This equation can be rewritten as:
  //
  // sigma_chi(i) * sigma_x(i) * corr_chi_x(i)
  //   = Coef_rt(i) * sigma_rt(i) * sigma_x(i) * corr_rt_x(i)
  //     - Coef_thl(i) * sigma_thl(i) * sigma_x(i) * corr_thl_x(i).
  //
  // This equation can be solved for corr_chi_x(i):
  //
  // corr_chi_x(i)
  // = Coef_rt(i) * ( sigma_rt(i) / sigma_chi(i) ) * corr_rt_x(i)
  //   - Coef_thl(i) * ( sigma_thl(i) / sigma_chi(i) ) * corr_thl_x(i).
  //
  // The correlation of chi and x within the ith component of the PDF is
  // calculated.

  // References:
  // Eq. (13) and Eq. (14) of Larson, V. E., R. Wood, P. R. Field, J.-C. Golaz,
  // T. H. Vonder Haar, W. R. Cotton, 2001:  Systematic Biases in the
  // Microphysics and Thermodynamics of Numerical Models That Ignore
  // Subgrid-Scale Variability. J. Atmos. Sci., 58, 1117--1128,
  // doi:https://doi.org/10.1175/1520-0469(2001)058%3C1117:SBITMA%3E2.0.CO;2.
  //
  // Eq. (A29) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //-----------------------------------------------------------------------

  // Input Variables
  // crt_i: Coefficient of r_t for chi (old s) (ith PDF comp.)   [-]
  // cthl_i: Coefficient of th_l for chi (ith PDF comp.)  [(kg/kg)/K]
  // sigma_rt_i: Standard deviation of r_t (ith PDF component)    [kg/kg]
  // sigma_thl_i: Standard deviation of th_l (ith PDF component)       [K]
  // sigma_chi_i: Standard deviation of chi (ith PDF component)    [kg/kg]
  // corr_rt_x_i: Correlation of r_t and x (ith PDF component)         [-]
  // corr_thl_x_i: Correlation of th_l and x (ith PDF component)        [-]

  // Return Variable
  real_t corr_chi_x_i;

  // Calculate the correlation of chi and x in the ith PDF component.
  if (sigma_chi_i > chi_tol) {

    corr_chi_x_i = crt_i * (sigma_rt_i / sigma_chi_i) * corr_rt_x_i -
                   cthl_i * (sigma_thl_i / sigma_chi_i) * corr_thl_x_i;

  } else {

    // The standard deviation of chi in the ith PDF component is 0.  This
    // means that chi is constant within the ith PDF component, and the ith
    // PDF component covariance of chi and x is also 0.  The correlation of
    // chi and x is undefined in the ith PDF component, so a value of 0 will
    // be used.
    corr_chi_x_i = zero;
  }

  // Clip the magnitude of the correlation of chi and x in the ith PDF
  // component, just in case the correlations and standard deviations used in
  // calculating it are inconsistent, resulting in an unrealizable value for
  // corr_chi_x_i.
  if (corr_chi_x_i > max_mag_correlation) {
    corr_chi_x_i = max_mag_correlation;
  } else if (corr_chi_x_i < -max_mag_correlation) {
    corr_chi_x_i = -max_mag_correlation;
  }

  return corr_chi_x_i;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t calc_corr_eta_x(real_t crt_i, real_t cthl_i,            // In
                                              real_t sigma_rt_i, real_t sigma_thl_i,  // In
                                              real_t sigma_eta_i, real_t corr_rt_x_i, // In
                                              real_t corr_thl_x_i) {
  // Description:
  // This function calculates the correlation of the variable that is
  // orthogonal to extended liquid water mixing ratio in a PDF transformation,
  // eta (old t), and a generic variable x, within the ith component of
  // the PDF.

  // References:
  // Eq. (A30) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //-----------------------------------------------------------------------

  // Input Variables
  // crt_i: Coefficient of r_t for chi (old s) (ith PDF comp.)   [-]
  // cthl_i: Coefficient of th_l for chi (ith PDF comp.)  [(kg/kg)/K]
  // sigma_rt_i: Standard deviation of r_t (ith PDF component)    [kg/kg]
  // sigma_thl_i: Standard deviation of th_l (ith PDF component)       [K]
  // sigma_eta_i: Standard deviation of eta (ith PDF component)    [kg/kg]
  // corr_rt_x_i: Correlation of r_t and x (ith PDF component)         [-]
  // corr_thl_x_i: Correlation of th_l and x (ith PDF component)        [-]

  // Return Variable
  real_t corr_eta_x_i;

  // Calculate the correlation of eta and x in the ith PDF component.
  if (sigma_eta_i > eta_tol) {

    corr_eta_x_i = crt_i * (sigma_rt_i / sigma_eta_i) * corr_rt_x_i +
                   cthl_i * (sigma_thl_i / sigma_eta_i) * corr_thl_x_i;

  } else {

    // The standard deviation of eta in the ith PDF component is 0.  This
    // means that eta is constant within the ith PDF component, and the ith
    // PDF component covariance of eta and x is also 0.  The correlation of
    // eta and x is undefined in the ith PDF component, so a value of 0 will
    // be used.
    corr_eta_x_i = zero;
  }

  // Clip the magnitude of the correlation of eta and x in the ith PDF
  // component, just in case the correlations and standard deviations used in
  // calculating it are inconsistent, resulting in an unrealizable value for
  // corr_eta_x_i.
  if (corr_eta_x_i > max_mag_correlation) {
    corr_eta_x_i = max_mag_correlation;
  } else if (corr_eta_x_i < -max_mag_correlation) {
    corr_eta_x_i = -max_mag_correlation;
  }

  return corr_eta_x_i;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t calc_corr_rt_x(real_t crt_i, real_t sigma_rt_i,
                                             real_t sigma_chi_i, // In
                                             real_t sigma_eta_i, real_t corr_chi_x_i,
                                             real_t corr_eta_x_i) {
  // Description:
  // This function calculates the correlation of rt and x based on the
  // correlation of chi and x and the correlation of eta and x.

  // References:
  //-----------------------------------------------------------------------

  // Input Variables
  // crt_i: Coef. of r_t in chi/eta eqns. (ith PDF component)    [-]
  // sigma_rt_i: Standard deviation of r_t (ith PDF component)    [kg/kg]
  // sigma_chi_i: Standard deviation of chi (ith PDF component)    [kg/kg]
  // sigma_eta_i: Standard deviation of eta (ith PDF component)    [kg/kg]
  // corr_chi_x_i: Correlation of chi and x (ith PDF component)         [-]
  // corr_eta_x_i: Correlation of eta and x (ith PDF component)         [-]

  // Return Variable
  real_t corr_rt_x_i;

  // Calculate the correlation of rt and x in the ith PDF component.
  if (sigma_rt_i > rt_tol) {

    corr_rt_x_i =
        (sigma_eta_i * corr_eta_x_i + sigma_chi_i * corr_chi_x_i) / (two * crt_i * sigma_rt_i);

  } else {

    // The standard deviation of rt in the ith PDF component is 0.  This means
    // that rt is constant within the ith PDF component, and the ith PDF
    // component covariance of rt and x is also 0.  The correlation of rt and
    // x is undefined in the ith PDF component, so a value of 0 will be used.
    corr_rt_x_i = zero;
  }

  // Clip the magnitude of the correlation of rt and x in the ith PDF
  // component, just in case the correlations and standard deviations used in
  // calculating it are inconsistent, resulting in an unrealizable value for
  // corr_rt_x_i.
  if (corr_rt_x_i > max_mag_correlation) {
    corr_rt_x_i = max_mag_correlation;
  } else if (corr_rt_x_i < -max_mag_correlation) {
    corr_rt_x_i = -max_mag_correlation;
  }

  return corr_rt_x_i;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t calc_corr_thl_x(real_t cthl_i, real_t sigma_thl_i,
                                              real_t sigma_chi_i,                      // In
                                              real_t sigma_eta_i, real_t corr_chi_x_i, // In
                                              real_t corr_eta_x_i) {
  // Description:
  // This function calculates the correlation of thl and x based on the
  // correlation of chi and x and the correlation of eta and x.

  // References:
  //-----------------------------------------------------------------------

  // Input Variables
  // cthl_i: Coef. of thl:  chi/eta eqns. (ith PDF comp.) [(kg/kg)/K]
  // sigma_thl_i: Standard deviation of thl (ith PDF component)        [K]
  // sigma_chi_i: Standard deviation of chi (ith PDF component)    [kg/kg]
  // sigma_eta_i: Standard deviation of eta (ith PDF component)    [kg/kg]
  // corr_chi_x_i: Correlation of chi and x (ith PDF component)         [-]
  // corr_eta_x_i: Correlation of eta and x (ith PDF component)         [-]

  // Return Variable
  real_t corr_thl_x_i;

  // Calculate the correlation of thl and x in the ith PDF component.
  if (sigma_thl_i > thl_tol) {

    corr_thl_x_i =
        (sigma_eta_i * corr_eta_x_i - sigma_chi_i * corr_chi_x_i) / (two * cthl_i * sigma_thl_i);

  } else {

    // The standard deviation of thl in the ith PDF component is 0.  This
    // means that thl is constant within the ith PDF component, and the ith
    // PDF component covariance of thl and x is also 0.  The correlation of
    // thl and x is undefined in the ith PDF component, so a value of 0 will
    // be used.
    corr_thl_x_i = zero;
  }

  // Clip the magnitude of the correlation of thl and x in the ith PDF
  // component, just in case the correlations and standard deviations used in
  // calculating it are inconsistent, resulting in an unrealizable value for
  // corr_thl_x_i.
  if (corr_thl_x_i > max_mag_correlation) {
    corr_thl_x_i = max_mag_correlation;
  } else if (corr_thl_x_i < -max_mag_correlation) {
    corr_thl_x_i = -max_mag_correlation;
  }

  return corr_thl_x_i;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t calc_xp2(real_t mu_x_1, real_t mu_x_2,                       // In
                                       real_t mu_x_1_n, real_t mu_x_2_n,                   // In
                                       real_t sigma_x_1, real_t sigma_x_2,                 // In
                                       real_t sigma_x_1_n, real_t sigma_x_2_n,             // In
                                       real_t mixt_frac, real_t x_frac_1, real_t x_frac_2, // In
                                       real_t x_mean) {
  // Description:
  // Calculates the overall variance of x, <x'^2>, where the distribution of x
  // is a combination of a lognormal distribution and/or 0 in each PDF
  // component.  The fraction of each component where x is lognormally
  // distributed (amd greater than 0) is x_frac_i (x_frac_1 and x_frac_2 for
  // PDF components 1 and 2, respectively).  The fraction of each component
  // where x has a value of 0 is ( 1 - x_frac_i ).  This function should be
  // called to calculate the total variance for x when <x'^2> is not provided
  // by a predictive (or other) equation.
  //
  // This function is used to calculate the overall variance for rain water
  // mixing ratio, <r_r'^2>, and the overall variance for rain drop
  // concentration, <N_r'^2>.  The ratio of variance to mean-value-squared is
  // specified for the in-precip values of r_r and N_r within each PDF
  // component, allowing for the calculation of sigma_rr_i and sigma_Nr_i,
  // as well as sigma_rr_i_n and sigma_Nr_i_n.

  // References:
  //-----------------------------------------------------------------------

  // Input Variables
  // mu_x_1: Mean of x (1st PDF comp.) in x_frac                  [-]
  // mu_x_2: Mean of x (2nd PDF comp.) in x_frac                  [-]
  // mu_x_1_n: Mean of ln x (1st PDF comp.) in x_frac               [-]
  // mu_x_2_n: Mean of ln x (2nd PDF comp.) in x_frac               [-]
  // sigma_x_1: Standard deviation of x (1st PDF comp.) in x_frac    [-]
  // sigma_x_2: Standard deviation of x (2nd PDF comp.) in x_frac    [-]
  // sigma_x_1_n: Standard deviation of ln x (1st PDF comp.) in x_frac [-]
  // sigma_x_2_n: Standard deviation of ln x (2nd PDF comp.) in x_frac [-]
  // mixt_frac: Mixture fraction                                     [-]
  // x_frac_1: Fraction: x distributed lognormally (1st PDF comp.)  [-]
  // x_frac_2: Fraction: x distributed lognormally (2nd PDF comp.)  [-]
  // x_mean: Overall mean value of x                              [-]

  // Return Variable
  real_t xp2;

  // Calculate overall variance of x, <x'^2>.
  if (Kokkos::abs(sigma_x_1) < eps && Kokkos::abs(sigma_x_2) < eps) {

    // The value of x is constant within both PDF components.
    xp2 = (mixt_frac * x_frac_1 * (mu_x_1 * mu_x_1) +
           (one - mixt_frac) * x_frac_2 * (mu_x_2 * mu_x_2)) -
          (x_mean * x_mean);

  } else if (Kokkos::abs(sigma_x_1) < eps) {

    // The value of x is constant within the 1st PDF component.
    xp2 = (mixt_frac * x_frac_1 * (mu_x_1 * mu_x_1) +
           (one - mixt_frac) * x_frac_2 *
               Kokkos::exp(two * mu_x_2_n + two * (sigma_x_2_n * sigma_x_2_n))) -
          (x_mean * x_mean);

  } else if (Kokkos::abs(sigma_x_2) < eps) {

    // The value of x is constant within the 2nd PDF component.
    xp2 = (mixt_frac * x_frac_1 * Kokkos::exp(two * mu_x_1_n + two * (sigma_x_1_n * sigma_x_1_n)) +
           (one - mixt_frac) * x_frac_2 * (mu_x_2 * mu_x_2)) -
          (x_mean * x_mean);

  } else {

    // The value of x varies within both PDF component.
    xp2 = (mixt_frac * x_frac_1 * Kokkos::exp(two * mu_x_1_n + two * (sigma_x_1_n * sigma_x_1_n)) +
           (one - mixt_frac) * x_frac_2 *
               Kokkos::exp(two * mu_x_2_n + two * (sigma_x_2_n * sigma_x_2_n))) -
          (x_mean * x_mean);
  }

  // As a check, prevent negative values for hydrometeor variances due to
  // numerical loss of precision error.
  if (xp2 < zero) {
    xp2 = zero;
  }

  return xp2;
}

// =============================================================================
inline void smooth_corr_quotient(int ngrdcol, int nz, Real2D numerator, Real2D denominator, // In
                                 real_t denom_thresh,                                       // In
                                 Real2D quotient) {                                         // Out
  // Bound the correlation smoothly by first raising its denominator toward
  // abs(numerator)/max_mag_correlation, then toward denom_thresh. numerator
  // and denominator have the same units; quotient is dimensionless.
  Kokkos::parallel_for(
      "smooth_corr_quotient", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        quotient(i, k) = smooth_corr_quotient(numerator(i, k), denominator(i, k), denom_thresh);
      });
}

// =============================================================================
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Output>
inline void compute_variance_binormal(Field0 xm, Field1 mu_x_1, Field2 mu_x_2, Field3 stdev_x_1,
                                      Field4 stdev_x_2, Field5 mixt_frac, Output xp2) {
  // Fill the caller's conformable (column, level) result with the native
  // scalar formula. Aligned input/output mappings may alias pointwise.
  static_assert(Output::rank == 2);
  static_assert(std::is_same_v<typename Output::non_const_value_type, real_t>);
  static_assert(Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                           typename Output::memory_space>::accessible);
  static_assert(Field0::rank == 2 && Field1::rank == 2 && Field2::rank == 2 && Field3::rank == 2 &&
                Field4::rank == 2 && Field5::rank == 2);
  static_assert(Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                           typename Field0::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field1::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field2::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field3::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field4::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field5::memory_space>::accessible,
                "PDF utility arrays must be execution-space Views");
  const int ngrdcol = xm.extent_int(0);
  const int nz      = xm.extent_int(1);
  Kokkos::parallel_for(
      "compute_variance_binormal::columns", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        xp2(i, k) = compute_variance_binormal(xm(i, k), mu_x_1(i, k), mu_x_2(i, k), stdev_x_1(i, k),
                                              stdev_x_2(i, k), mixt_frac(i, k));
      });
}

template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5>
inline Real2D compute_variance_binormal(Field0 xm, Field1 mu_x_1, Field2 mu_x_2, Field3 stdev_x_1,
                                        Field4 stdev_x_2, Field5 mixt_frac) {
  // Return independent owned storage through the same direct-output kernel.
  Real2D xp2(Kokkos::view_alloc(Kokkos::WithoutInitializing, "compute_variance_binormal::result"),
             xm.extent(0), xm.extent(1));
  compute_variance_binormal(xm, mu_x_1, mu_x_2, stdev_x_1, stdev_x_2, mixt_frac, xp2);
  return xp2;
}

// =============================================================================
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Output>
inline void calc_corr_chi_x(Field0 crt_i, Field1 cthl_i, Field2 sigma_rt_i, Field3 sigma_thl_i,
                            Field4 sigma_chi_i, Field5 corr_rt_x_i, Field6 corr_thl_x_i,
                            Output corr_chi_x_i) {
  // Fill the caller's conformable (column, level) result with the native
  // scalar formula. Aligned input/output mappings may alias pointwise.
  static_assert(Output::rank == 2);
  static_assert(std::is_same_v<typename Output::non_const_value_type, real_t>);
  static_assert(Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                           typename Output::memory_space>::accessible);
  static_assert(Field0::rank == 2 && Field1::rank == 2 && Field2::rank == 2 && Field3::rank == 2 &&
                Field4::rank == 2 && Field5::rank == 2 && Field6::rank == 2);
  static_assert(Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                           typename Field0::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field1::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field2::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field3::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field4::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field5::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field6::memory_space>::accessible,
                "PDF utility arrays must be execution-space Views");
  const int ngrdcol = crt_i.extent_int(0);
  const int nz      = crt_i.extent_int(1);
  Kokkos::parallel_for(
      "calc_corr_chi_x::columns", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        corr_chi_x_i(i, k) =
            calc_corr_chi_x(crt_i(i, k), cthl_i(i, k), sigma_rt_i(i, k), sigma_thl_i(i, k),
                            sigma_chi_i(i, k), corr_rt_x_i(i, k), corr_thl_x_i(i, k));
      });
}

template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6>
inline Real2D calc_corr_chi_x(Field0 crt_i, Field1 cthl_i, Field2 sigma_rt_i, Field3 sigma_thl_i,
                              Field4 sigma_chi_i, Field5 corr_rt_x_i, Field6 corr_thl_x_i) {
  // Return independent owned storage through the same direct-output kernel.
  Real2D corr_chi_x_i(Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_corr_chi_x::result"),
                      crt_i.extent(0), crt_i.extent(1));
  calc_corr_chi_x(crt_i, cthl_i, sigma_rt_i, sigma_thl_i, sigma_chi_i, corr_rt_x_i, corr_thl_x_i,
                  corr_chi_x_i);
  return corr_chi_x_i;
}

// =============================================================================
template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6, class Output>
inline void calc_corr_eta_x(Field0 crt_i, Field1 cthl_i, Field2 sigma_rt_i, Field3 sigma_thl_i,
                            Field4 sigma_eta_i, Field5 corr_rt_x_i, Field6 corr_thl_x_i,
                            Output corr_eta_x_i) {
  // Fill the caller's conformable (column, level) result with the native
  // scalar formula. Aligned input/output mappings may alias pointwise.
  static_assert(Output::rank == 2);
  static_assert(std::is_same_v<typename Output::non_const_value_type, real_t>);
  static_assert(Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                           typename Output::memory_space>::accessible);
  static_assert(Field0::rank == 2 && Field1::rank == 2 && Field2::rank == 2 && Field3::rank == 2 &&
                Field4::rank == 2 && Field5::rank == 2 && Field6::rank == 2);
  static_assert(Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                           typename Field0::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field1::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field2::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field3::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field4::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field5::memory_space>::accessible &&
                    Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                               typename Field6::memory_space>::accessible,
                "PDF utility arrays must be execution-space Views");
  const int ngrdcol = crt_i.extent_int(0);
  const int nz      = crt_i.extent_int(1);
  Kokkos::parallel_for(
      "calc_corr_eta_x::columns", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        corr_eta_x_i(i, k) =
            calc_corr_eta_x(crt_i(i, k), cthl_i(i, k), sigma_rt_i(i, k), sigma_thl_i(i, k),
                            sigma_eta_i(i, k), corr_rt_x_i(i, k), corr_thl_x_i(i, k));
      });
}

template <class Field0, class Field1, class Field2, class Field3, class Field4, class Field5,
          class Field6>
inline Real2D calc_corr_eta_x(Field0 crt_i, Field1 cthl_i, Field2 sigma_rt_i, Field3 sigma_thl_i,
                              Field4 sigma_eta_i, Field5 corr_rt_x_i, Field6 corr_thl_x_i) {
  // Return independent owned storage through the same direct-output kernel.
  Real2D corr_eta_x_i(Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_corr_eta_x::result"),
                      crt_i.extent(0), crt_i.extent(1));
  calc_corr_eta_x(crt_i, cthl_i, sigma_rt_i, sigma_thl_i, sigma_eta_i, corr_rt_x_i, corr_thl_x_i,
                  corr_eta_x_i);
  return corr_eta_x_i;
}

} // namespace clubb_kokkos
