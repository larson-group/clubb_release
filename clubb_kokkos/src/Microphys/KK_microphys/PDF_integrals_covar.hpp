// Native PDF_integrals_covar numerical formulas, evaluated in source arithmetic order.
// Appended caller-owned KK parameters replace implicit module state.
#pragma once
#include "Microphys/KK_microphys/engine_contracts.hpp"
#include "Microphys/KK_microphys/KK_utilities.hpp"
namespace clubb_kokkos::kk_engine {
inline double quadrivar_NNLL_covar(double mu_x1, double mu_x2, double mu_x3_n, double mu_x4_n,
                                   double sigma_x1, double sigma_x2, double sigma_x3_n,
                                   double sigma_x4_n, double rho_x1x2, double rho_x1x3_n,
                                   double rho_x1x4_n, double rho_x2x3_n, double rho_x2x4_n,
                                   double rho_x3x4_n, double x1_mean,
                                   double x2_alpha_x3_beta_x4_gamma_mean, double alpha_exp,
                                   double beta_exp, double gamma_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J1) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S9) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variable
  double quadrivar_NNLL_covar;
  double s_cc;

  s_cc = (((mu_x2 / sigma_x2) + ((rho_x2x3_n * sigma_x3_n) * beta_exp)) +
          ((rho_x2x4_n * sigma_x4_n) * gamma_exp));

  // BFB-024: Native forms the lognormal variance/cross group separately from
  // the mean pair. Preserve its fused correlation differences and -O0 sum order.
  const double variance_x3 =
      ((one_half_dp * Kokkos::fma(-rho_x2x3_n, rho_x2x3_n, one_dp)) *
       (sigma_x3_n * sigma_x3_n)) * (beta_exp * beta_exp);
  const double variance_x4_prefix =
      (one_half_dp * Kokkos::fma(-rho_x2x4_n, rho_x2x4_n, one_dp)) *
      (sigma_x4_n * sigma_x4_n);
  const double variance_pair =
      Kokkos::fma(variance_x4_prefix, gamma_exp * gamma_exp, variance_x3);
  const double cross_prefix =
      ((Kokkos::fma(-rho_x2x3_n, rho_x2x4_n, rho_x3x4_n) * sigma_x3_n) * beta_exp) * sigma_x4_n;
  const double variance_cross = Kokkos::fma(cross_prefix, gamma_exp, variance_pair);
  const double mean_pair = Kokkos::fma(mu_x4_n, gamma_exp, mu_x3_n * beta_exp);
  const double lognormal_exponent = variance_cross + mean_pair;

  quadrivar_NNLL_covar =
      ((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * Kokkos::pow((-sigma_x2), alpha_exp)) *
          Kokkos::exp(lognormal_exponent)) *
         Kokkos::exp((((one_fourth_dp * (s_cc * s_cc)) - ((mu_x2 / sigma_x2) * s_cc)) +
                      (one_half_dp * ((mu_x2 * mu_x2) / (sigma_x2 * sigma_x2)))))) *
        (((((-rho_x1x2) * sigma_x1) * parabolic::gamma((alpha_exp + two_dp))) *
          Dv_fnc((-(alpha_exp + two_dp)), s_cc, kk)) +
         ((((((mu_x1 - x1_mean) - (((mu_x2 / sigma_x2) * rho_x1x2) * sigma_x1)) +
             ((((rho_x1x3_n - (rho_x1x2 * rho_x2x3_n)) * sigma_x1) * sigma_x3_n) * beta_exp)) +
            ((((rho_x1x4_n - (rho_x1x2 * rho_x2x4_n)) * sigma_x1) * sigma_x4_n) * gamma_exp)) *
           parabolic::gamma((alpha_exp + one_dp))) *
          Dv_fnc((-(alpha_exp + one_dp)), s_cc, kk)))) -
       (x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean)));

  return quadrivar_NNLL_covar;
}
inline double quadrivar_NNLL_covar_const_x1(double mu_x1, double mu_x2, double mu_x3_n,
                                            double mu_x4_n, double sigma_x2, double sigma_x3_n,
                                            double sigma_x4_n, double rho_x2x3_n, double rho_x2x4_n,
                                            double rho_x3x4_n, double x1_mean,
                                            double x2_alpha_x3_beta_x4_gamma_mean, double alpha_exp,
                                            double beta_exp, double gamma_exp,
                                            const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J2) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S10) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variable
  double quadrivar_NNLL_covar_const_x1;
  double s_cc;

  s_cc = (((mu_x2 / sigma_x2) + ((rho_x2x3_n * sigma_x3_n) * beta_exp)) +
          ((rho_x2x4_n * sigma_x4_n) * gamma_exp));

  // BFB-028: Native forms the lognormal variance/cross group separately from
  // the mean pair. Preserve its fused correlation differences and -O0 sum order.
  const double variance_x3 =
      ((one_half_dp * Kokkos::fma(-rho_x2x3_n, rho_x2x3_n, one_dp)) *
       (sigma_x3_n * sigma_x3_n)) * (beta_exp * beta_exp);
  const double variance_x4_prefix =
      (one_half_dp * Kokkos::fma(-rho_x2x4_n, rho_x2x4_n, one_dp)) *
      (sigma_x4_n * sigma_x4_n);
  const double variance_pair =
      Kokkos::fma(variance_x4_prefix, gamma_exp * gamma_exp, variance_x3);
  const double cross_prefix =
      ((Kokkos::fma(-rho_x2x3_n, rho_x2x4_n, rho_x3x4_n) * sigma_x3_n) * beta_exp) * sigma_x4_n;
  const double variance_cross = Kokkos::fma(cross_prefix, gamma_exp, variance_pair);
  const double mean_pair = Kokkos::fma(mu_x4_n, gamma_exp, mu_x3_n * beta_exp);
  const double lognormal_exponent = variance_cross + mean_pair;

  quadrivar_NNLL_covar_const_x1 =
      ((((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * (mu_x1 - x1_mean)) *
            Kokkos::pow((-sigma_x2), alpha_exp)) *
           Kokkos::exp(lognormal_exponent)) *
          Kokkos::exp((((one_fourth_dp * (s_cc * s_cc)) - ((mu_x2 / sigma_x2) * s_cc)) +
                       (one_half_dp * ((mu_x2 * mu_x2) / (sigma_x2 * sigma_x2)))))) *
         parabolic::gamma((alpha_exp + one_dp))) *
        Dv_fnc((-(alpha_exp + one_dp)), s_cc, kk)) -
       (x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean)));

  return quadrivar_NNLL_covar_const_x1;
}
inline double quadrivar_NNLL_covar_const_x2(double mu_x1, double mu_x2, double mu_x3_n,
                                            double mu_x4_n, double sigma_x1, double sigma_x3_n,
                                            double sigma_x4_n, double rho_x1x3_n, double rho_x1x4_n,
                                            double rho_x3x4_n, double x1_mean,
                                            double x2_alpha_x3_beta_x4_gamma_mean, double alpha_exp,
                                            double beta_exp, double gamma_exp,
                                            const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J3) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S11) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double quadrivar_NNLL_covar_const_x2;

  if ((mu_x2 <= zero_dp)) {

    // BFB-028: Native contracts variance/cross separately from the mean pair.
    const double variance_x3 =
        (one_half_dp * (sigma_x3_n * sigma_x3_n)) * (beta_exp * beta_exp);
    const double variance_pair =
        Kokkos::fma(one_half_dp * (sigma_x4_n * sigma_x4_n), gamma_exp * gamma_exp, variance_x3);
    const double cross_prefix = ((rho_x3x4_n * sigma_x3_n) * beta_exp) * sigma_x4_n;
    const double variance_cross = Kokkos::fma(cross_prefix, gamma_exp, variance_pair);
    const double mean_pair = Kokkos::fma(mu_x4_n, gamma_exp, mu_x3_n * beta_exp);
    const double lognormal_exponent = variance_cross + mean_pair;

    quadrivar_NNLL_covar_const_x2 =
        (((Kokkos::pow(mu_x2, alpha_exp) *
           (((mu_x1 - x1_mean) + (((rho_x1x3_n * sigma_x1) * sigma_x3_n) * beta_exp)) +
            (((rho_x1x4_n * sigma_x1) * sigma_x4_n) * gamma_exp))) *
          Kokkos::exp(lognormal_exponent)) -
         (x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean)));

  } else {

    quadrivar_NNLL_covar_const_x2 = ((-x2_alpha_x3_beta_x4_gamma_mean) * (mu_x1 - x1_mean));
  }

  return quadrivar_NNLL_covar_const_x2;
}
inline double quadrivar_NNLL_covar_const_x3(double mu_x1, double mu_x2, double mu_x3,
                                            double mu_x4_n, double sigma_x1, double sigma_x2,
                                            double sigma_x4_n, double rho_x1x2, double rho_x1x4_n,
                                            double rho_x2x4_n, double x1_mean,
                                            double x2_alpha_x3_beta_x4_gamma_mean, double alpha_exp,
                                            double beta_exp, double gamma_exp,
                                            const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J4) and Eq. (J5) of Griffin, B. M., 2016:  Improving the
  // Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
  // Effects Using a Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S12) and Eq. (S13) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of Parameterizing microphysical effects on variances and
  // covariances of moisture and heat content using a multivariate probability
  // density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
  // 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variable
  double quadrivar_NNLL_covar_const_x3;
  double s_cc;

  s_cc = ((mu_x2 / sigma_x2) + ((rho_x2x4_n * sigma_x4_n) * gamma_exp));

  quadrivar_NNLL_covar_const_x3 =
      ((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * Kokkos::pow((-sigma_x2), alpha_exp)) *
          Kokkos::pow(mu_x3, beta_exp)) *
         Kokkos::exp((((mu_x4_n * gamma_exp) +
                       ((one_half_dp * (sigma_x4_n * sigma_x4_n)) * (gamma_exp * gamma_exp))) -
                      (one_fourth_dp * (s_cc * s_cc))))) *
        (((((-rho_x1x2) * sigma_x1) * parabolic::gamma((alpha_exp + two_dp))) *
          Dv_fnc((-(alpha_exp + two_dp)), s_cc, kk)) +
         (((((mu_x1 - x1_mean) - (((mu_x2 / sigma_x2) * rho_x1x2) * sigma_x1)) +
            ((((rho_x1x4_n - (rho_x1x2 * rho_x2x4_n)) * sigma_x1) * sigma_x4_n) * gamma_exp)) *
           parabolic::gamma((alpha_exp + one_dp))) *
          Dv_fnc((-(alpha_exp + one_dp)), s_cc, kk)))) -
       (x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean)));

  return quadrivar_NNLL_covar_const_x3;
}
inline double quadrivar_NNLL_covar_const_x1x2(double mu_x1, double mu_x2, double mu_x3_n,
                                              double mu_x4_n, double sigma_x3_n, double sigma_x4_n,
                                              double rho_x3x4_n, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J6) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S14) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double quadrivar_NNLL_covar_const_x1x2;

  if ((mu_x2 <= zero_dp)) {

    // BFB-028: Native contracts variance/cross separately from the mean pair.
    const double variance_x3 =
        (one_half_dp * (sigma_x3_n * sigma_x3_n)) * (beta_exp * beta_exp);
    const double variance_pair =
        Kokkos::fma(one_half_dp * (sigma_x4_n * sigma_x4_n), gamma_exp * gamma_exp, variance_x3);
    const double cross_prefix = ((rho_x3x4_n * sigma_x3_n) * beta_exp) * sigma_x4_n;
    const double variance_cross = Kokkos::fma(cross_prefix, gamma_exp, variance_pair);
    const double mean_pair = Kokkos::fma(mu_x4_n, gamma_exp, mu_x3_n * beta_exp);
    const double lognormal_exponent = variance_cross + mean_pair;

    quadrivar_NNLL_covar_const_x1x2 =
        ((((mu_x1 - x1_mean) * Kokkos::pow(mu_x2, alpha_exp)) *
          Kokkos::exp(lognormal_exponent)) -
         (x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean)));

  } else {

    quadrivar_NNLL_covar_const_x1x2 = ((-x2_alpha_x3_beta_x4_gamma_mean) * (mu_x1 - x1_mean));
  }

  return quadrivar_NNLL_covar_const_x1x2;
}
inline double quadrivar_NNLL_covar_const_x1x3(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4_n, double sigma_x2, double sigma_x4_n,
                                              double rho_x2x4_n, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J7) and Eq. (J8) of Griffin, B. M., 2016:  Improving the
  // Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
  // Effects Using a Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S15) and Eq. (S16) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of Parameterizing microphysical effects on variances and
  // covariances of moisture and heat content using a multivariate probability
  // density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
  // 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variable
  double quadrivar_NNLL_covar_const_x1x3;
  double s_cc;

  s_cc = ((mu_x2 / sigma_x2) + ((rho_x2x4_n * sigma_x4_n) * gamma_exp));

  quadrivar_NNLL_covar_const_x1x3 =
      ((((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * (mu_x1 - x1_mean)) *
            Kokkos::pow((-sigma_x2), alpha_exp)) *
           Kokkos::pow(mu_x3, beta_exp)) *
          Kokkos::exp((((mu_x4_n * gamma_exp) +
                        ((one_half_dp * (sigma_x4_n * sigma_x4_n)) * (gamma_exp * gamma_exp))) -
                       (one_fourth_dp * (s_cc * s_cc))))) *
         parabolic::gamma((alpha_exp + one_dp))) *
        Dv_fnc((-(alpha_exp + one_dp)), s_cc, kk)) -
       (x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean)));

  return quadrivar_NNLL_covar_const_x1x3;
}
inline double quadrivar_NNLL_covar_const_x2x3(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4_n, double sigma_x1, double sigma_x4_n,
                                              double rho_x1x4_n, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J9) and Eq. (J10) of Griffin, B. M., 2016:  Improving the
  // Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
  // Effects Using a Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S17) and Eq. (S18) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of Parameterizing microphysical effects on variances and
  // covariances of moisture and heat content using a multivariate probability
  // density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
  // 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double quadrivar_NNLL_covar_const_x2x3;

  if ((mu_x2 <= zero_dp)) {

    quadrivar_NNLL_covar_const_x2x3 =
        ((((Kokkos::pow(mu_x2, alpha_exp) * Kokkos::pow(mu_x3, beta_exp)) *
           ((mu_x1 - x1_mean) + (((rho_x1x4_n * sigma_x1) * sigma_x4_n) * gamma_exp))) *
          Kokkos::exp(((mu_x4_n * gamma_exp) +
                       ((one_half_dp * (sigma_x4_n * sigma_x4_n)) * (gamma_exp * gamma_exp))))) -
         (x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean)));

  } else {

    quadrivar_NNLL_covar_const_x2x3 = ((-x2_alpha_x3_beta_x4_gamma_mean) * (mu_x1 - x1_mean));
  }

  return quadrivar_NNLL_covar_const_x2x3;
}
inline double quadrivar_NNLL_covar_const_x3x4(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4, double sigma_x1, double sigma_x2,
                                              double rho_x1x2, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J11) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S19) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double quadrivar_NNLL_covar_const_x3x4;

  quadrivar_NNLL_covar_const_x3x4 =
      (((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * Kokkos::pow((-sigma_x2), alpha_exp)) *
           Kokkos::pow(mu_x3, beta_exp)) *
          Kokkos::pow(mu_x4, gamma_exp)) *
         Kokkos::exp(((-one_fourth_dp) * ((mu_x2 * mu_x2) / (sigma_x2 * sigma_x2))))) *
        (((((-rho_x1x2) * sigma_x1) * parabolic::gamma((alpha_exp + two_dp))) *
          Dv_fnc((-(alpha_exp + two_dp)), (mu_x2 / sigma_x2), kk)) +
         ((((mu_x1 - x1_mean) - (((mu_x2 / sigma_x2) * rho_x1x2) * sigma_x1)) *
           parabolic::gamma((alpha_exp + one_dp))) *
          Dv_fnc((-(alpha_exp + one_dp)), (mu_x2 / sigma_x2), kk)))) -
       (x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean)));

  return quadrivar_NNLL_covar_const_x3x4;
}
inline double quadrivar_NNLL_covar_cst_x1x2x3(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4_n, double sigma_x4_n, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J12) and Eq. (J13) of Griffin, B. M., 2016:  Improving the
  // Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
  // Effects Using a Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S20) and Eq. (S21) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of Parameterizing microphysical effects on variances and
  // covariances of moisture and heat content using a multivariate probability
  // density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
  // 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double quadrivar_NNLL_covar_cst_x1x2x3;

  if ((mu_x2 <= zero_dp)) {

    quadrivar_NNLL_covar_cst_x1x2x3 =
        (((((mu_x1 - x1_mean) * Kokkos::pow(mu_x2, alpha_exp)) * Kokkos::pow(mu_x3, beta_exp)) *
          Kokkos::exp(((mu_x4_n * gamma_exp) +
                       ((one_half_dp * (sigma_x4_n * sigma_x4_n)) * (gamma_exp * gamma_exp))))) -
         (x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean)));

  } else {

    quadrivar_NNLL_covar_cst_x1x2x3 = ((-x2_alpha_x3_beta_x4_gamma_mean) * (mu_x1 - x1_mean));
  }

  return quadrivar_NNLL_covar_cst_x1x2x3;
}
inline double quadrivar_NNLL_covar_cst_x1x3x4(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4, double sigma_x2, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J14) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S22) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double quadrivar_NNLL_covar_cst_x1x3x4;

  quadrivar_NNLL_covar_cst_x1x3x4 =
      (((((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * (mu_x1 - x1_mean)) *
             Kokkos::pow((-sigma_x2), alpha_exp)) *
            Kokkos::pow(mu_x3, beta_exp)) *
           Kokkos::pow(mu_x4, gamma_exp)) *
          Kokkos::exp(((-one_fourth_dp) * ((mu_x2 * mu_x2) / (sigma_x2 * sigma_x2))))) *
         parabolic::gamma((alpha_exp + one_dp))) *
        Dv_fnc((-(alpha_exp + one_dp)), (mu_x2 / sigma_x2), kk)) -
       (x2_alpha_x3_beta_x4_gamma_mean * (mu_x1 - x1_mean)));

  return quadrivar_NNLL_covar_cst_x1x3x4;
}
inline double quadrivar_NNLL_covar_cst_x2x3x4(double mu_x1, double mu_x2, double mu_x3,
                                              double mu_x4, double x1_mean,
                                              double x2_alpha_x3_beta_x4_gamma_mean,
                                              double alpha_exp, double beta_exp, double gamma_exp,
                                              const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J15) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S23) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double quadrivar_NNLL_covar_cst_x2x3x4;

  if ((mu_x2 <= zero_dp)) {

    quadrivar_NNLL_covar_cst_x2x3x4 =
        ((mu_x1 - x1_mean) * (((Kokkos::pow(mu_x2, alpha_exp) * Kokkos::pow(mu_x3, beta_exp)) *
                               Kokkos::pow(mu_x4, gamma_exp)) -
                              x2_alpha_x3_beta_x4_gamma_mean));

  } else {

    quadrivar_NNLL_covar_cst_x2x3x4 = ((-x2_alpha_x3_beta_x4_gamma_mean) * (mu_x1 - x1_mean));
  }

  return quadrivar_NNLL_covar_cst_x2x3x4;
}
inline double quadrivar_NNLL_covar_const_all(double mu_x1, double mu_x2, double mu_x3, double mu_x4,
                                             double x1_mean, double x2_alpha_x3_beta_x4_gamma_mean,
                                             double alpha_exp, double beta_exp, double gamma_exp,
                                             const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J16) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S24) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double quadrivar_NNLL_covar_const_all;

  if ((mu_x2 <= zero_dp)) {

    quadrivar_NNLL_covar_const_all =
        ((mu_x1 - x1_mean) * (((Kokkos::pow(mu_x2, alpha_exp) * Kokkos::pow(mu_x3, beta_exp)) *
                               Kokkos::pow(mu_x4, gamma_exp)) -
                              x2_alpha_x3_beta_x4_gamma_mean));

  } else {

    quadrivar_NNLL_covar_const_all = ((-x2_alpha_x3_beta_x4_gamma_mean) * (mu_x1 - x1_mean));
  }

  return quadrivar_NNLL_covar_const_all;
}
inline double trivar_NNL_covar(double mu_x1, double mu_x2, double mu_x3_n, double sigma_x1,
                               double sigma_x2, double sigma_x3_n, double rho_x1x2,
                               double rho_x1x3_n, double rho_x2x3_n, double x1_mean,
                               double x2_alpha_x3_beta_mean, double alpha_exp, double beta_exp,
                               const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J17) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S25) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variable
  double trivar_NNL_covar;
  double s_c;

  s_c = ((mu_x2 / sigma_x2) + ((rho_x2x3_n * sigma_x3_n) * beta_exp));

  trivar_NNL_covar =
      (((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * Kokkos::pow(sigma_x2, alpha_exp)) *
         Kokkos::exp((((mu_x3_n * beta_exp) +
                       ((one_half_dp * (sigma_x3_n * sigma_x3_n)) * (beta_exp * beta_exp))) -
                      (one_fourth_dp * (s_c * s_c))))) *
        ((((rho_x1x2 * sigma_x1) * parabolic::gamma((alpha_exp + two_dp))) *
          Dv_fnc((-(alpha_exp + two_dp)), (-s_c), kk)) +
         (((((mu_x1 - x1_mean) - (((mu_x2 / sigma_x2) * rho_x1x2) * sigma_x1)) +
            ((((rho_x1x3_n - (rho_x1x2 * rho_x2x3_n)) * sigma_x1) * sigma_x3_n) * beta_exp)) *
           parabolic::gamma((alpha_exp + one_dp))) *
          Dv_fnc((-(alpha_exp + one_dp)), (-s_c), kk)))) -
       (x2_alpha_x3_beta_mean * (mu_x1 - x1_mean)));

  return trivar_NNL_covar;
}
inline double trivar_NNL_covar_const_x1(double mu_x1, double mu_x2, double mu_x3_n, double sigma_x2,
                                        double sigma_x3_n, double rho_x2x3_n, double x1_mean,
                                        double x2_alpha_x3_beta_mean, double alpha_exp,
                                        double beta_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J18) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S26) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variable
  double trivar_NNL_covar_const_x1;
  double s_c;

  s_c = ((mu_x2 / sigma_x2) + ((rho_x2x3_n * sigma_x3_n) * beta_exp));

  trivar_NNL_covar_const_x1 =
      (((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * (mu_x1 - x1_mean)) *
           Kokkos::pow(sigma_x2, alpha_exp)) *
          Kokkos::exp((((mu_x3_n * beta_exp) +
                        ((one_half_dp * (sigma_x3_n * sigma_x3_n)) * (beta_exp * beta_exp))) -
                       (one_fourth_dp * (s_c * s_c))))) *
         parabolic::gamma((alpha_exp + one_dp))) *
        Dv_fnc((-(alpha_exp + one_dp)), (-s_c), kk)) -
       (x2_alpha_x3_beta_mean * (mu_x1 - x1_mean)));

  return trivar_NNL_covar_const_x1;
}
inline double trivar_NNL_covar_const_x2(double mu_x1, double mu_x2, double mu_x3_n, double sigma_x1,
                                        double sigma_x3_n, double rho_x1x3_n, double x1_mean,
                                        double x2_alpha_x3_beta_mean, double alpha_exp,
                                        double beta_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J19) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S27) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double trivar_NNL_covar_const_x2;

  if ((mu_x2 >= zero_dp)) {

    trivar_NNL_covar_const_x2 =
        (((Kokkos::pow(mu_x2, alpha_exp) *
           ((mu_x1 - x1_mean) + (((rho_x1x3_n * sigma_x1) * sigma_x3_n) * beta_exp))) *
          Kokkos::exp(((mu_x3_n * beta_exp) +
                       ((one_half_dp * (sigma_x3_n * sigma_x3_n)) * (beta_exp * beta_exp))))) -
         (x2_alpha_x3_beta_mean * (mu_x1 - x1_mean)));

  } else {

    trivar_NNL_covar_const_x2 = ((-x2_alpha_x3_beta_mean) * (mu_x1 - x1_mean));
  }

  return trivar_NNL_covar_const_x2;
}
inline double trivar_NNL_covar_const_x3(double mu_x1, double mu_x2, double mu_x3, double sigma_x1,
                                        double sigma_x2, double rho_x1x2, double x1_mean,
                                        double x2_alpha_x3_beta_mean, double alpha_exp,
                                        double beta_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J20) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S28) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double trivar_NNL_covar_const_x3;

  trivar_NNL_covar_const_x3 =
      ((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * Kokkos::pow(sigma_x2, alpha_exp)) *
          Kokkos::pow(mu_x3, beta_exp)) *
         Kokkos::exp(((-one_fourth_dp) * ((mu_x2 * mu_x2) / (sigma_x2 * sigma_x2))))) *
        ((((rho_x1x2 * sigma_x1) * parabolic::gamma((alpha_exp + two_dp))) *
          Dv_fnc((-(alpha_exp + two_dp)), (-(mu_x2 / sigma_x2)), kk)) +
         ((((mu_x1 - x1_mean) - (((mu_x2 / sigma_x2) * rho_x1x2) * sigma_x1)) *
           parabolic::gamma((alpha_exp + one_dp))) *
          Dv_fnc((-(alpha_exp + one_dp)), (-(mu_x2 / sigma_x2)), kk)))) -
       (x2_alpha_x3_beta_mean * (mu_x1 - x1_mean)));

  return trivar_NNL_covar_const_x3;
}
inline double trivar_NNL_covar_const_x1x2(double mu_x1, double mu_x2, double mu_x3_n,
                                          double sigma_x3_n, double x1_mean,
                                          double x2_alpha_x3_beta_mean, double alpha_exp,
                                          double beta_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J21) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S29) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double trivar_NNL_covar_const_x1x2;

  if ((mu_x2 >= zero_dp)) {

    trivar_NNL_covar_const_x1x2 =
        (((Kokkos::pow(mu_x2, alpha_exp) * (mu_x1 - x1_mean)) *
          Kokkos::exp(((mu_x3_n * beta_exp) +
                       ((one_half_dp * (sigma_x3_n * sigma_x3_n)) * (beta_exp * beta_exp))))) -
         (x2_alpha_x3_beta_mean * (mu_x1 - x1_mean)));

  } else {

    trivar_NNL_covar_const_x1x2 = ((-x2_alpha_x3_beta_mean) * (mu_x1 - x1_mean));
  }

  return trivar_NNL_covar_const_x1x2;
}
inline double trivar_NNL_covar_const_x1x3(double mu_x1, double mu_x2, double mu_x3, double sigma_x2,
                                          double x1_mean, double x2_alpha_x3_beta_mean,
                                          double alpha_exp, double beta_exp,
                                          const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J22) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S30) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double trivar_NNL_covar_const_x1x3;

  trivar_NNL_covar_const_x1x3 =
      ((((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * (mu_x1 - x1_mean)) *
            Kokkos::pow(sigma_x2, alpha_exp)) *
           Kokkos::pow(mu_x3, beta_exp)) *
          Kokkos::exp(((-one_fourth_dp) * ((mu_x2 * mu_x2) / (sigma_x2 * sigma_x2))))) *
         parabolic::gamma((alpha_exp + one_dp))) *
        Dv_fnc((-(alpha_exp + one_dp)), (-(mu_x2 / sigma_x2)), kk)) -
       (x2_alpha_x3_beta_mean * (mu_x1 - x1_mean)));

  return trivar_NNL_covar_const_x1x3;
}
inline double trivar_NNL_covar_const_x2x3(double mu_x1, double mu_x2, double mu_x3, double x1_mean,
                                          double x2_alpha_x3_beta_mean, double alpha_exp,
                                          double beta_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J23) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S31) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double trivar_NNL_covar_const_x2x3;

  if ((mu_x2 >= zero_dp)) {

    trivar_NNL_covar_const_x2x3 =
        ((mu_x1 - x1_mean) *
         ((Kokkos::pow(mu_x2, alpha_exp) * Kokkos::pow(mu_x3, beta_exp)) - x2_alpha_x3_beta_mean));

  } else {

    trivar_NNL_covar_const_x2x3 = ((-x2_alpha_x3_beta_mean) * (mu_x1 - x1_mean));
  }

  return trivar_NNL_covar_const_x2x3;
}
inline double trivar_NNL_covar_const_all(double mu_x1, double mu_x2, double mu_x3, double x1_mean,
                                         double x2_alpha_x3_beta_mean, double alpha_exp,
                                         double beta_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (J24) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S32) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double trivar_NNL_covar_const_all;

  if ((mu_x2 >= zero_dp)) {

    trivar_NNL_covar_const_all =
        ((mu_x1 - x1_mean) *
         ((Kokkos::pow(mu_x2, alpha_exp) * Kokkos::pow(mu_x3, beta_exp)) - x2_alpha_x3_beta_mean));

  } else {

    trivar_NNL_covar_const_all = ((-x2_alpha_x3_beta_mean) * (mu_x1 - x1_mean));
  }

  return trivar_NNL_covar_const_all;
}
} // namespace clubb_kokkos::kk_engine
