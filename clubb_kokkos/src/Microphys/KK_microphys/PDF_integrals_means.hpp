// Native PDF_integrals_means numerical formulas, evaluated in source arithmetic order.
// Appended caller-owned KK parameters replace implicit module state.
#pragma once
#include "Microphys/KK_microphys/engine_contracts.hpp"
#include "Microphys/KK_microphys/KK_utilities.hpp"
namespace clubb_kokkos::kk_engine {
inline double trivar_NLL_mean(double mu_x1, double mu_x2_n, double mu_x3_n, double sigma_x1,
                              double sigma_x2_n, double sigma_x3_n, double rho_x1x2_n,
                              double rho_x1x3_n, double rho_x2x3_n, double alpha_exp,
                              double beta_exp, double gamma_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (60) and Eq. (61) of Larson, V. E. and B. M. Griffin, 2013:  Analytic
  // upscaling of a local microphysics scheme. Part I: Derivation.
  // Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
  // doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Eq. (C24) and Eq. (J25) of Griffin, B. M., 2016:  Improving the
  // Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
  // Effects Using a Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S24) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Eq. (S33) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variable
  double trivar_NLL_mean;
  double s_cc;

  s_cc = (((mu_x1 / sigma_x1) + ((rho_x1x2_n * sigma_x2_n) * beta_exp)) +
          ((rho_x1x3_n * sigma_x3_n) * gamma_exp));

  trivar_NLL_mean =
      (((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * Kokkos::pow((-sigma_x1), alpha_exp)) *
           Kokkos::exp(((mu_x2_n * beta_exp) + (mu_x3_n * gamma_exp)))) *
          Kokkos::exp(
              (one_half_dp *
               (((((one_dp - (rho_x1x2_n * rho_x1x2_n)) * (sigma_x2_n * sigma_x2_n)) *
                  (beta_exp * beta_exp)) +
                 (((one_dp - (rho_x1x3_n * rho_x1x3_n)) * (sigma_x3_n * sigma_x3_n)) *
                  (gamma_exp * gamma_exp))) +
                (((((two_dp * (rho_x2x3_n - (rho_x1x2_n * rho_x1x3_n))) * sigma_x2_n) * beta_exp) *
                  sigma_x3_n) *
                 gamma_exp))))) *
         Kokkos::exp((((one_fourth_dp * (s_cc * s_cc)) - ((mu_x1 / sigma_x1) * s_cc)) +
                      (one_half_dp * ((mu_x1 * mu_x1) / (sigma_x1 * sigma_x1)))))) *
        parabolic::gamma((alpha_exp + one_dp))) *
       Dv_fnc((-(alpha_exp + one_dp)), s_cc, kk));

  return trivar_NLL_mean;
}
inline double trivar_NLL_mean_const_x1(double mu_x1, double mu_x2_n, double mu_x3_n,
                                       double sigma_x2_n, double sigma_x3_n, double rho_x2x3_n,
                                       double alpha_exp, double beta_exp, double gamma_exp,
                                       const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (62) and Eq. (63) of Larson, V. E. and B. M. Griffin, 2013:  Analytic
  // upscaling of a local microphysics scheme. Part I: Derivation.
  // Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
  // doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Eq. (C29) and Eq. (J26) of Griffin, B. M., 2016:  Improving the
  // Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
  // Effects Using a Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S29) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Eq. (S34) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double trivar_NLL_mean_const_x1;

  if ((mu_x1 <= zero_dp)) {

    // BFB-026: Preserve Native's separately contracted variance/cross and mean groups.
    const double variance_x2 =
        (one_half_dp * (sigma_x2_n * sigma_x2_n)) * (beta_exp * beta_exp);
    const double variance_pair =
        Kokkos::fma(one_half_dp * (sigma_x3_n * sigma_x3_n), gamma_exp * gamma_exp, variance_x2);
    const double cross_prefix = ((rho_x2x3_n * sigma_x2_n) * beta_exp) * sigma_x3_n;
    const double variance_cross = Kokkos::fma(cross_prefix, gamma_exp, variance_pair);
    const double mean_pair = Kokkos::fma(mu_x3_n, gamma_exp, mu_x2_n * beta_exp);
    trivar_NLL_mean_const_x1 =
        Kokkos::pow(mu_x1, alpha_exp) * Kokkos::exp(variance_cross + mean_pair);

  } else {

    trivar_NLL_mean_const_x1 = zero_dp;
  }

  return trivar_NLL_mean_const_x1;
}
inline double trivar_NLL_mean_const_x2(double mu_x1, double mu_x2, double mu_x3_n, double sigma_x1,
                                       double sigma_x3_n, double rho_x1x3_n, double alpha_exp,
                                       double beta_exp, double gamma_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (C26), Eq. (C27), Eq. (J27), and Eq. (J28) of Griffin, B. M., 2016:
  // Improving the Subgrid-Scale Representation of Hydrometeors and
  // Microphysical Feedback Effects Using a Multivariate PDF.  Doctoral
  // dissertation, University of Wisconsin -- Milwaukee, Milwaukee, WI,
  // Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S26) and Eq. (S27) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of A new subgrid-scale representation of hydrometeor fields
  // using a multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Eq. (S35) and Eq. (S36) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of Parameterizing microphysical effects on variances and
  // covariances of moisture and heat content using a multivariate probability
  // density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
  // 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variable
  double trivar_NLL_mean_const_x2;
  double s_cc;

  s_cc = ((mu_x1 / sigma_x1) + ((rho_x1x3_n * sigma_x3_n) * gamma_exp));

  trivar_NLL_mean_const_x2 =
      ((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * Kokkos::pow((-sigma_x1), alpha_exp)) *
          Kokkos::pow(mu_x2, beta_exp)) *
         Kokkos::exp((((mu_x3_n * gamma_exp) +
                       ((one_half_dp * (sigma_x3_n * sigma_x3_n)) * (gamma_exp * gamma_exp))) -
                      (one_fourth_dp * (s_cc * s_cc))))) *
        parabolic::gamma((alpha_exp + one_dp))) *
       Dv_fnc((-(alpha_exp + one_dp)), s_cc, kk));

  return trivar_NLL_mean_const_x2;
}
inline double trivar_NLL_mean_const_x1x2(double mu_x1, double mu_x2, double mu_x3_n,
                                         double sigma_x3_n, double alpha_exp, double beta_exp,
                                         double gamma_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (C30), Eq. (C31), Eq. (J29), and Eq. (J30) of Griffin, B. M., 2016:
  // Improving the Subgrid-Scale Representation of Hydrometeors and
  // Microphysical Feedback Effects Using a Multivariate PDF.  Doctoral
  // dissertation, University of Wisconsin -- Milwaukee, Milwaukee, WI,
  // Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S30) and Eq. (S31) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of A new subgrid-scale representation of hydrometeor fields
  // using a multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Eq. (S37) and Eq. (S38) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of Parameterizing microphysical effects on variances and
  // covariances of moisture and heat content using a multivariate probability
  // density function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9,
  // 11, doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double trivar_NLL_mean_const_x1x2;

  if ((mu_x1 <= zero_dp)) {

    trivar_NLL_mean_const_x1x2 =
        ((Kokkos::pow(mu_x1, alpha_exp) * Kokkos::pow(mu_x2, beta_exp)) *
         Kokkos::exp(((mu_x3_n * gamma_exp) +
                      ((one_half_dp * (sigma_x3_n * sigma_x3_n)) * (gamma_exp * gamma_exp)))));

  } else {

    trivar_NLL_mean_const_x1x2 = zero_dp;
  }

  return trivar_NLL_mean_const_x1x2;
}
inline double trivar_NLL_mean_const_x2x3(double mu_x1, double mu_x2, double mu_x3, double sigma_x1,
                                         double alpha_exp, double beta_exp, double gamma_exp,
                                         const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (C28) and Eq. (J31) of Griffin, B. M., 2016:  Improving the
  // Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
  // Effects Using a Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S28) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Eq. (S39) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double trivar_NLL_mean_const_x2x3;

  trivar_NLL_mean_const_x2x3 =
      (((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * Kokkos::pow((-sigma_x1), alpha_exp)) *
           Kokkos::pow(mu_x2, beta_exp)) *
          Kokkos::pow(mu_x3, gamma_exp)) *
         Kokkos::exp(((-one_fourth_dp) * ((mu_x1 * mu_x1) / (sigma_x1 * sigma_x1))))) *
        parabolic::gamma((alpha_exp + one_dp))) *
       Dv_fnc((-(alpha_exp + one_dp)), (mu_x1 / sigma_x1), kk));

  return trivar_NLL_mean_const_x2x3;
}
inline double trivar_NLL_mean_const_all(double mu_x1, double mu_x2, double mu_x3, double alpha_exp,
                                        double beta_exp, double gamma_exp,
                                        const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (C32) and Eq. (J32) of Griffin, B. M., 2016:  Improving the
  // Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
  // Effects Using a Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S32) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Eq. (S40) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double trivar_NLL_mean_const_all;

  if ((mu_x1 <= zero_dp)) {

    trivar_NLL_mean_const_all = ((Kokkos::pow(mu_x1, alpha_exp) * Kokkos::pow(mu_x2, beta_exp)) *
                                 Kokkos::pow(mu_x3, gamma_exp));

  } else {

    trivar_NLL_mean_const_all = zero_dp;
  }

  return trivar_NLL_mean_const_all;
}
inline double bivar_NL_mean(double mu_x1, double mu_x2_n, double sigma_x1, double sigma_x2_n,
                            double rho_x1x2_n, double alpha_exp, double beta_exp,
                            const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (33), Eq. (34), Eq. (41), and Eq. (42) of Larson, V. E. and
  // B. M. Griffin, 2013:  Analytic upscaling of a local microphysics scheme.
  // Part I: Derivation.  Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
  // doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Eq. (C4), Eq. (C14), and Eq. (J33) of Griffin, B. M., 2016:  Improving
  // the Subgrid-Scale Representation of Hydrometeors and Microphysical
  // Feedback Effects Using a Multivariate PDF.  Doctoral dissertation,
  // University of Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp.,
  // URL http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S4) and Eq. (S14) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of A new subgrid-scale representation of hydrometeor fields
  // using a multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Eq. (S41) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variable
  double bivar_NL_mean;
  double s_c;

  s_c = ((mu_x1 / sigma_x1) + ((rho_x1x2_n * sigma_x2_n) * beta_exp));

  bivar_NL_mean =
      (((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * Kokkos::pow(sigma_x1, alpha_exp)) *
         Kokkos::exp((((mu_x2_n * beta_exp) +
                       ((one_half_dp * (sigma_x2_n * sigma_x2_n)) * (beta_exp * beta_exp))) -
                      (one_fourth_dp * (s_c * s_c))))) *
        parabolic::gamma((alpha_exp + one_dp))) *
       Dv_fnc((-(alpha_exp + one_dp)), (-s_c), kk));

  return bivar_NL_mean;
}
inline double bivar_NL_mean_const_x1(double mu_x1, double mu_x2_n, double sigma_x2_n,
                                     double alpha_exp, double beta_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (35), Eq. (36), Eq. (43), and Eq. (44) of Larson, V. E. and
  // B. M. Griffin, 2013:  Analytic upscaling of a local microphysics scheme.
  // Part I: Derivation.  Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
  // doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Eq. (C9), Eq. (C16), and Eq. (J34) of Griffin, B. M., 2016:  Improving
  // the Subgrid-Scale Representation of Hydrometeors and Microphysical
  // Feedback Effects Using a Multivariate PDF.  Doctoral dissertation,
  // University of Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp.,
  // URL http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S9) and Eq. (S16) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of A new subgrid-scale representation of hydrometeor fields
  // using a multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Eq. (S42) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double bivar_NL_mean_const_x1;

  if ((mu_x1 >= zero_dp)) {

    bivar_NL_mean_const_x1 =
        (Kokkos::pow(mu_x1, alpha_exp) *
         Kokkos::exp(((mu_x2_n * beta_exp) +
                      ((one_half_dp * (sigma_x2_n * sigma_x2_n)) * (beta_exp * beta_exp)))));

  } else {

    bivar_NL_mean_const_x1 = zero_dp;
  }

  return bivar_NL_mean_const_x1;
}
inline double bivar_NL_mean_const_x2(double mu_x1, double mu_x2, double sigma_x1, double alpha_exp,
                                     double beta_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (C8), Eq. (C15), and Eq. (J35) of Griffin, B. M., 2016:  Improving
  // the Subgrid-Scale Representation of Hydrometeors and Microphysical
  // Feedback Effects Using a Multivariate PDF.  Doctoral dissertation,
  // University of Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp.,
  // URL http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S8) and Eq. (S15) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of A new subgrid-scale representation of hydrometeor fields
  // using a multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Eq. (S43) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double bivar_NL_mean_const_x2;

  bivar_NL_mean_const_x2 =
      ((((((one_dp / Kokkos::sqrt((two_dp * pi_dp))) * Kokkos::pow(sigma_x1, alpha_exp)) *
          Kokkos::pow(mu_x2, beta_exp)) *
         Kokkos::exp(((-one_fourth_dp) * ((mu_x1 * mu_x1) / (sigma_x1 * sigma_x1))))) *
        parabolic::gamma((alpha_exp + one_dp))) *
       Dv_fnc((-(alpha_exp + one_dp)), (-(mu_x1 / sigma_x1)), kk));

  return bivar_NL_mean_const_x2;
}
inline double bivar_NL_mean_const_all(double mu_x1, double mu_x2, double alpha_exp, double beta_exp,
                                      const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (C10), Eq. (C17), and Eq. (J36) of Griffin, B. M., 2016:  Improving
  // the Subgrid-Scale Representation of Hydrometeors and Microphysical
  // Feedback Effects Using a Multivariate PDF.  Doctoral dissertation,
  // University of Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp.,
  // URL http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S10) and Eq. (S17) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of A new subgrid-scale representation of hydrometeor fields
  // using a multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Eq. (S44) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double bivar_NL_mean_const_all;

  if ((mu_x1 >= zero_dp)) {

    bivar_NL_mean_const_all = (Kokkos::pow(mu_x1, alpha_exp) * Kokkos::pow(mu_x2, beta_exp));

  } else {

    bivar_NL_mean_const_all = zero_dp;
  }

  return bivar_NL_mean_const_all;
}
inline double bivar_LL_mean(double mu_x1_n, double mu_x2_n, double sigma_x1_n, double sigma_x2_n,
                            double rho_x1x2_n, double alpha_exp, double beta_exp,
                            const parameters_KK& kk) {
  // Description:
  // Calculates the analytic solution to the integral:
  //
  // INT(0:INF) INT(0:INF) x1^alpha x2^beta P_LL(i)(x1,x2) dx2 dx1;
  //
  // where x1 and x2 are both variables that have a lognormal marginal
  // distribution in the ith component of the PDF.
  //
  // The analytic solution to the integral is:
  //
  // INT(0:INF) INT(0:INF) x1^alpha x2^beta P_LL(i)(x1,x2) dx2 dx1
  // = exp{ mu_x1_n * alpha + mu_x2_n * beta
  //        + 1/2 * sigma_x1_n^2 * alpha_exp^2
  //        + 1/2 * sigma_x2_n^2 * beta_exp^2
  //        + rho_x1x2_n * sigma_x1_n * alpha * sigma_x2_n * beta }
  // References:
  // Eq. (26) of Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of
  // a local microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol.
  // Soc., 139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Eq. (C39) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S39) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double bivar_LL_mean;

  // BFB-023: Native groups variance/cross terms separately from the mean pair
  // before exp. Explicit FMAs retain its -O0 NVHPC contraction and addition order.
  const double variance_x1 =
      (one_half_dp * (sigma_x1_n * sigma_x1_n)) * (alpha_exp * alpha_exp);
  const double variance_pair =
      Kokkos::fma(one_half_dp * (sigma_x2_n * sigma_x2_n), beta_exp * beta_exp, variance_x1);
  const double variance_cross =
      Kokkos::fma(((rho_x1x2_n * sigma_x1_n) * alpha_exp) * sigma_x2_n, beta_exp, variance_pair);
  const double mean_pair = Kokkos::fma(mu_x2_n, beta_exp, mu_x1_n * alpha_exp);
  bivar_LL_mean = Kokkos::exp(variance_cross + mean_pair);

  return bivar_LL_mean;
}
inline double bivar_LL_mean_const_x1(double mu_x1, double mu_x2_n, double sigma_x2_n,
                                     double alpha_exp, double beta_exp, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (C40) and Eq. (C41) of Griffin, B. M., 2016:  Improving the
  // Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
  // Effects Using a Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S40) and Eq. (S41) of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of A new subgrid-scale representation of hydrometeor fields
  // using a multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double bivar_LL_mean_const_x1;

  bivar_LL_mean_const_x1 =
      (Kokkos::pow(mu_x1, alpha_exp) *
       Kokkos::exp(((mu_x2_n * beta_exp) +
                    ((one_half_dp * (sigma_x2_n * sigma_x2_n)) * (beta_exp * beta_exp)))));

  return bivar_LL_mean_const_x1;
}
inline double bivar_LL_mean_const_all(double mu_x1, double mu_x2, double alpha_exp, double beta_exp,
                                      const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (C42) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (S42) of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double bivar_LL_mean_const_all;

  bivar_LL_mean_const_all = (Kokkos::pow(mu_x1, alpha_exp) * Kokkos::pow(mu_x2, beta_exp));

  return bivar_LL_mean_const_all;
}
} // namespace clubb_kokkos::kk_engine
