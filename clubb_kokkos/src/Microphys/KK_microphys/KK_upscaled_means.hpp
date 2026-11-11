// Native KK_upscaled_means numerical formulas, evaluated in source arithmetic order.
// Appended caller-owned KK parameters replace implicit module state.
#pragma once
#include "Microphys/KK_microphys/engine_contracts.hpp"
#include "Microphys/KK_microphys/KK_utilities.hpp"
namespace clubb_kokkos::kk_engine {
inline void KK_upscaled_means_driver(
    real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1,
    real_t mu_Nr_2, real_t mu_Ncn_1, real_t mu_Ncn_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
    real_t mu_Nr_1_n, real_t mu_Nr_2_n, real_t mu_Ncn_1_n, real_t mu_Ncn_2_n, real_t sigma_chi_1,
    real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2,
    real_t sigma_Ncn_1, real_t sigma_Ncn_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n,
    real_t sigma_Nr_1_n, real_t sigma_Nr_2_n, real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n,
    real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n, real_t corr_chi_Nr_1_n, real_t corr_chi_Nr_2_n,
    real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
    real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2, real_t KK_evap_coef,
    real_t KK_auto_coef, real_t KK_accr_coef, real_t KK_mvr_coef, real_t& KK_evap_tndcy,
    real_t& KK_auto_tndcy, real_t& KK_accr_tndcy, real_t& KK_mean_vol_rad,
    const parameters_KK& kk) {
  // Description:
  // References:
  // Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of a local
  //    microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol. Soc.,
  //    139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Griffin, B. M. and V. E. Larson, 2016:  Supplement of A new subgrid-scale
  //    representation of hydrometeor fields using a multivariate PDF.
  //    Geosci. Model Dev., 9, 6,
  //    doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Output Variables

  //!! Calculate the upscaled KK evaporation tendency.
  KK_evap_tndcy = KK_evap_upscaled_mean(
      mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n,
      mu_Nr_2_n, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2,
      sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n,
      corr_chi_Nr_1_n, corr_chi_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_evap_coef, mixt_frac,
      precip_frac_1, precip_frac_2, kk);

  //!! Calculate the upscaled KK autoconversion tendency.
  KK_auto_tndcy = KK_auto_upscaled_mean(mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n,
                                        mu_Ncn_2_n, sigma_chi_1, sigma_chi_2, sigma_Ncn_1,
                                        sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, corr_chi_Ncn_1_n,
                                        corr_chi_Ncn_2_n, KK_auto_coef, mixt_frac, kk);

  //!! Calculate the upscaled KK accretion tendency.
  KK_accr_tndcy = KK_accr_upscaled_mean(
      mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_chi_1, sigma_chi_2,
      sigma_rr_1, sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n,
      KK_accr_coef, mixt_frac, precip_frac_1, precip_frac_2, kk);

  //!! Calculate the upscaled KK rain drop mean volume radius.
  KK_mean_vol_rad = KK_mvr_upscaled_mean(
      mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1,
      sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
      corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mvr_coef, mixt_frac, precip_frac_1, precip_frac_2, kk);

  return;
}
inline real_t KK_evap_upscaled_mean(
    real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1,
    real_t mu_Nr_2, real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n, real_t mu_Nr_2_n,
    real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1,
    real_t sigma_Nr_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t sigma_Nr_1_n,
    real_t sigma_Nr_2_n, real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n, real_t corr_chi_Nr_1_n,
    real_t corr_chi_Nr_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t KK_evap_coef,
    real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2, const parameters_KK& kk) {
  // Description:
  // This function calculates the mean value of the upscaled KK rain water
  // evaporation tendency.
  // References:
  // Section 7 of Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of
  // a local microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol.
  // Soc., 139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Section C.3 of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Section S1.3 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t KK_evap_upscaled_mean;
  real_t alpha_exp;
  real_t beta_exp;
  real_t gamma_exp;

  // Values of the KK exponents.
  alpha_exp = KK_evap_Supersat_exp;
  beta_exp  = KK_evap_rr_exp;
  gamma_exp = KK_evap_Nr_exp;

  // Calculate the mean KK evaporation tendency.
  KK_evap_upscaled_mean =
      (KK_evap_coef *
       (((mixt_frac * precip_frac_1) *
         trivar_NLL_mean_eq(mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_chi_1,
                            sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_chi_rr_1_n,
                            corr_chi_Nr_1_n, corr_rr_Nr_1_n, alpha_exp, beta_exp, gamma_exp, kk)) +
        (((one - mixt_frac) * precip_frac_2) *
         trivar_NLL_mean_eq(mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_chi_2,
                            sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_chi_rr_2_n,
                            corr_chi_Nr_2_n, corr_rr_Nr_2_n, alpha_exp, beta_exp, gamma_exp, kk))));

  return KK_evap_upscaled_mean;
}
inline real_t KK_auto_upscaled_mean(real_t mu_chi_1, real_t mu_chi_2, real_t mu_Ncn_1,
                                    real_t mu_Ncn_2, real_t mu_Ncn_1_n, real_t mu_Ncn_2_n,
                                    real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_Ncn_1,
                                    real_t sigma_Ncn_2, real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n,
                                    real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n,
                                    real_t KK_auto_coef, real_t mixt_frac,
                                    const parameters_KK& kk) {
  // Description:
  // This function calculates the mean value of the upscaled KK rain water
  // autoconversion tendency.
  //
  // Note:  Cloud droplet concentration is based on cloud nuclei concentration.
  //        Cloud nuclei concentration, N_cn, has replaced cloud droplet
  //        concentration, N_c, in CLUBB's PDF.
  // References:
  // Section 5 of Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of
  // a local microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol.
  // Soc., 139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Section C.2 of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Section S1.2 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t KK_auto_upscaled_mean;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_auto_rc_exp;
  beta_exp  = KK_auto_Nc_exp;

  // Calculate the mean KK autoconversion tendency.
  KK_auto_upscaled_mean =
      (KK_auto_coef *
       ((mixt_frac * bivar_NL_mean_eq(mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_chi_1, sigma_Ncn_1,
                                      sigma_Ncn_1_n, corr_chi_Ncn_1_n, Nc_tol, alpha_exp, beta_exp,
                                      kk)) +
        ((one - mixt_frac) * bivar_NL_mean_eq(mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_chi_2,
                                              sigma_Ncn_2, sigma_Ncn_2_n, corr_chi_Ncn_2_n, Nc_tol,
                                              alpha_exp, beta_exp, kk))));

  return KK_auto_upscaled_mean;
}
inline real_t KK_accr_upscaled_mean(real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1,
                                    real_t mu_rr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
                                    real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_rr_1,
                                    real_t sigma_rr_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n,
                                    real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n,
                                    real_t KK_accr_coef, real_t mixt_frac, real_t precip_frac_1,
                                    real_t precip_frac_2, const parameters_KK& kk) {
  // Description:
  // This function calculates the mean value of the upscaled KK rain water
  // accretion tendency.
  // References:
  // Section 6 of Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of
  // a local microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol.
  // Soc., 139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Section C.1 of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Section S1.1 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t KK_accr_upscaled_mean;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_accr_rc_exp;
  beta_exp  = KK_accr_rr_exp;

  // Calculate the mean KK accretion tendency.
  KK_accr_upscaled_mean =
      (KK_accr_coef *
       (((mixt_frac * precip_frac_1) * bivar_NL_mean_eq(mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_chi_1,
                                                        sigma_rr_1, sigma_rr_1_n, corr_chi_rr_1_n,
                                                        rr_tol, alpha_exp, beta_exp, kk)) +
        (((one - mixt_frac) * precip_frac_2) *
         bivar_NL_mean_eq(mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_chi_2, sigma_rr_2, sigma_rr_2_n,
                          corr_chi_rr_2_n, rr_tol, alpha_exp, beta_exp, kk))));

  return KK_accr_upscaled_mean;
}
inline real_t KK_mvr_upscaled_mean(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                                   real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                                   real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                                   real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                                   real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                                   real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t KK_mvr_coef,
                                   real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2,
                                   const parameters_KK& kk) {
  // Description:
  // This function calculates the mean value of the upscaled KK rain drop mean
  // volume radius.
  // References:
  // Section 4 of Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of
  // a local microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol.
  // Soc., 139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Section C.4 of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Section S1.4 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t KK_mvr_upscaled_mean;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_mvr_rr_exp;
  beta_exp  = KK_mvr_Nr_exp;

  // Calculate the KK mean volume radius of rain drops.
  KK_mvr_upscaled_mean =
      (KK_mvr_coef *
       (((mixt_frac * precip_frac_1) *
         bivar_LL_mean_eq(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_rr_1, sigma_Nr_1,
                          sigma_rr_1_n, sigma_Nr_1_n, corr_rr_Nr_1_n, alpha_exp, beta_exp, kk)) +
        (((one - mixt_frac) * precip_frac_2) *
         bivar_LL_mean_eq(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_rr_2, sigma_Nr_2,
                          sigma_rr_2_n, sigma_Nr_2_n, corr_rr_Nr_2_n, alpha_exp, beta_exp, kk))));

  return KK_mvr_upscaled_mean;
}
inline real_t trivar_NLL_mean_eq(real_t mu_chi_i, real_t mu_rr_i, real_t mu_Nr_i, real_t mu_rr_i_n,
                                 real_t mu_Nr_i_n, real_t sigma_chi_i, real_t sigma_rr_i,
                                 real_t sigma_Nr_i, real_t sigma_rr_i_n, real_t sigma_Nr_i_n,
                                 real_t corr_chi_rr_i_n, real_t corr_chi_Nr_i_n,
                                 real_t corr_rr_Nr_i_n, real_t alpha_exp_in, real_t beta_exp_in,
                                 real_t gamma_exp_in, const parameters_KK& kk) {
  // Description:
  // This function calculates the contribution by the ith PDF component to the
  // expression < x1^alpha x2^beta x3^gamma >, where x1 = chi, x2 = r_r, and
  // x3 = N_r.  The total value of KK mean evaporation tendency is given by:
  //
  // < KK_evap >
  //    = KK_evap_coef
  //      * ( mixt_frac < chi^alpha r_r^beta N_r^gamma (1) >
  //          + ( 1 - mixt_frac ) < chi^alpha r_r^beta N_r^gamma (2) > ).
  //
  // One of eight functions is called, based on whether any one, any two, all
  // three, or none of x1 (chi), x2 (r_r), and x3 (N_r) vary.  Each one of
  // these eight functions are the result of an evaluated integral based on the
  // specific situation.
  // References:
  // Section 7 of Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of
  // a local microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol.
  // Soc., 139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Section C.3 and Section J.3 of Griffin, B. M., 2016:  Improving the
  // Subgrid-Scale Representation of Hydrometeors and Microphysical Feedback
  // Effects Using a Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Section S1.3 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Section S7 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t trivar_NLL_mean_eq;
  double mu_x1;
  double mu_x2;
  double mu_x3;
  double mu_x2_n;
  double mu_x3_n;
  double sigma_x1;
  double sigma_x2;
  double sigma_x3;
  double sigma_x2_n;
  double sigma_x3_n;
  double rho_x1x2_n;
  double rho_x1x3_n;
  double rho_x2x3_n;
  double alpha_exp;
  double beta_exp;
  double gamma_exp;
  double x1_tol;
  double x2_tol;
  double x3_tol;
  double s_cc;

  // Means for the ith PDF component.
  mu_x1 = static_cast<double>(mu_chi_i);
  if ((beta_exp_in >= zero)) {
    mu_x2 = static_cast<double>(mu_rr_i);
  } else {
    mu_x2 = static_cast<double>(((mu_rr_i) > (rr_tol) ? (mu_rr_i) : (rr_tol)));
  }
  if ((gamma_exp_in >= zero)) {
    mu_x3 = static_cast<double>(mu_Nr_i);
  } else {
    mu_x3 = static_cast<double>(((mu_Nr_i) > (Nr_tol) ? (mu_Nr_i) : (Nr_tol)));
  }
  mu_x2_n = static_cast<double>(mu_rr_i_n);
  mu_x3_n = static_cast<double>(mu_Nr_i_n);

  // Standard deviations for the ith PDF component.
  sigma_x1   = static_cast<double>(sigma_chi_i);
  sigma_x2   = static_cast<double>(sigma_rr_i);
  sigma_x3   = static_cast<double>(sigma_Nr_i);
  sigma_x2_n = static_cast<double>(sigma_rr_i_n);
  sigma_x3_n = static_cast<double>(sigma_Nr_i_n);

  // Correlations for the ith PDF component.
  rho_x1x2_n = static_cast<double>(corr_chi_rr_i_n);
  rho_x1x3_n = static_cast<double>(corr_chi_Nr_i_n);
  rho_x2x3_n = static_cast<double>(corr_rr_Nr_i_n);

  // Exponents.
  alpha_exp = static_cast<double>(alpha_exp_in);
  beta_exp  = static_cast<double>(beta_exp_in);
  gamma_exp = static_cast<double>(gamma_exp_in);

  // Tolerance values.
  // When the standard deviation of a variable is below the tolerance values,
  // it is considered to be zero, and the variable is considered to have a
  // constant value.
  x1_tol = static_cast<double>(chi_tol);
  x2_tol = static_cast<double>(rr_tol);
  x3_tol = static_cast<double>(Nr_tol);

  // Determine the value of the parabolic cylinder function input value, s_cc.
  // The value s_cc is being fed into the parabolic cylinder function.  When
  // the value of s_cc is too large in magnitude (depending on the order of the
  // parabolic cylinder function), overflow occurs, and the output of the
  // parabolic cylinder function is +/-Inf.  This is primarily due to a large
  // ratio of mu_x1 to sigma_x1.  When the value of s_cc is very large, the
  // distribution of x1 is basically a spike near the mean, so x1 is treated as
  // a constant.
  if ((sigma_x1 > x1_tol)) {
    s_cc = (((mu_x1 / sigma_x1) + ((rho_x1x2_n * sigma_x2_n) * beta_exp)) +
            ((rho_x1x3_n * sigma_x3_n) * gamma_exp));
  } else {
    // Note:  s_cc is +inf when mu_x1 > 0 and sigma_x1 = 0, and s_cc is -inf
    //        when mu_x1 < 0 and sigma_x1 = 0.  Furthermore, s_cc is undefined
    //        when mu_x1 = 0 and sigma_x1 = 0.  However, within the context of
    //        this particular function, only the absolute value of s_cc is
    //        relevant, and furthermore the absolute value of s_cc is only
    //        relevant when sigma_x1 > 0.  Therefore, this statement only
    //        serves as divide-by-zero and compiler warning prevention.
    s_cc = std::numeric_limits<double>::max();
  }

  // Based on the value of sigma_x1 (including the value of s_cc compared to
  // parab_cyl_max_input), find the correct form of the trivariate equation to
  // use.

  if (((((sigma_x1 <= x1_tol) || (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input))) &&
        (sigma_x2 <= x2_tol)) &&
       (sigma_x3 <= x3_tol))) {

    // The ith PDF component variance of each of chi, r_r, and N_r is 0.
    trivar_NLL_mean_eq = static_cast<real_t>(
        trivar_NLL_mean_const_all(mu_x1, mu_x2, mu_x3, alpha_exp, beta_exp, gamma_exp, kk));

  } else {
    if ((((sigma_x1 <= x1_tol) || (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input))) &&
         (sigma_x2 <= x2_tol))) {

      // The ith PDF component variance of both chi and r_r is 0.
      trivar_NLL_mean_eq = static_cast<real_t>(trivar_NLL_mean_const_x1x2(
          mu_x1, mu_x2, mu_x3_n, sigma_x3_n, alpha_exp, beta_exp, gamma_exp, kk));

    } else {
      if ((((sigma_x1 <= x1_tol) ||
            (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input))) &&
           (sigma_x3 <= x3_tol))) {

        // The ith PDF component variance of both chi and N_r is 0.
        trivar_NLL_mean_eq = static_cast<real_t>(trivar_NLL_mean_const_x1x2(
            mu_x1, mu_x3, mu_x2_n, sigma_x2_n, alpha_exp, gamma_exp, beta_exp, kk));

      } else {
        if (((sigma_x2 <= x2_tol) && (sigma_x3 <= x3_tol))) {

          // The ith PDF component variance of both r_r and N_r is 0.
          trivar_NLL_mean_eq = static_cast<real_t>(trivar_NLL_mean_const_x2x3(
              mu_x1, mu_x2, mu_x3, sigma_x1, alpha_exp, beta_exp, gamma_exp, kk));

        } else {
          if (((sigma_x1 <= x1_tol) ||
               (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input)))) {

            // The ith PDF component variance of chi is 0.
            trivar_NLL_mean_eq = static_cast<real_t>(
                trivar_NLL_mean_const_x1(mu_x1, mu_x2_n, mu_x3_n, sigma_x2_n, sigma_x3_n,
                                         rho_x2x3_n, alpha_exp, beta_exp, gamma_exp, kk));

          } else {
            if ((sigma_x2 <= x2_tol)) {

              // The ith PDF component variance of r_r is 0.
              trivar_NLL_mean_eq = static_cast<real_t>(
                  trivar_NLL_mean_const_x2(mu_x1, mu_x2, mu_x3_n, sigma_x1, sigma_x3_n, rho_x1x3_n,
                                           alpha_exp, beta_exp, gamma_exp, kk));

            } else {
              if ((sigma_x3 <= x3_tol)) {

                // The ith PDF component variance of N_r is 0.
                trivar_NLL_mean_eq = static_cast<real_t>(
                    trivar_NLL_mean_const_x2(mu_x1, mu_x3, mu_x2_n, sigma_x1, sigma_x2_n,
                                             rho_x1x2_n, alpha_exp, gamma_exp, beta_exp, kk));

              } else {

                // All fields vary in the ith PDF component.
                trivar_NLL_mean_eq = static_cast<real_t>(trivar_NLL_mean(
                    mu_x1, mu_x2_n, mu_x3_n, sigma_x1, sigma_x2_n, sigma_x3_n, rho_x1x2_n,
                    rho_x1x3_n, rho_x2x3_n, alpha_exp, beta_exp, gamma_exp, kk));
              }
            }
          }
        }
      }
    }
  }

  return trivar_NLL_mean_eq;
}
inline real_t bivar_NL_mean_eq(real_t mu_chi_i, real_t mu_y_i, real_t mu_y_i_n, real_t sigma_chi_i,
                               real_t sigma_y_i, real_t sigma_y_i_n, real_t corr_chi_y_i_n,
                               real_t y_tol, real_t alpha_exp_in, real_t beta_exp_in,
                               const parameters_KK& kk) {
  // Description:
  // This function calculates the contribution by the ith PDF component to the
  // expression < x1^alpha x2^beta >, where x1 = chi and x2 = N_cn or r_r,
  // depending on whether this function is being called for autoconversion or
  // accretion, respectively.  The total value of KK mean microphysics tendency
  // is given by:
  //
  // < KK_mc > = KK_mc_coef
  //             * ( mixt_frac < chi^alpha y^beta (1) >
  //                 + ( 1 - mixt_frac ) < chi^alpha y^beta (2) > );
  //
  // where y stands for either N_cn or r_r.  One of four functions is called,
  // based on whether either one, both, or none of x1 (s) and x2 (N_cn or r_r)
  // vary.  Each one of these four functions are the result of an evaluated
  // integral based on the specific situation.
  // References:
  // Section 5 and Section 6 of Larson, V. E. and B. M. Griffin, 2013:
  // Analytic upscaling of a local microphysics scheme. Part I: Derivation.
  // Q. J. Roy. Meteorol. Soc., 139, 670, 46--57,
  // doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Section C.1, Section C.2, and Section J.4 of Griffin, B. M., 2016:
  // Improving the Subgrid-Scale Representation of Hydrometeors and
  // Microphysical Feedback Effects Using a Multivariate PDF.  Doctoral
  // dissertation, University of Wisconsin -- Milwaukee, Milwaukee, WI,
  // Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Section S1.1 and Section S1.2 of Griffin, B. M. and V. E. Larson, 2016:
  // Supplement of A new subgrid-scale representation of hydrometeor fields
  // using a multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //
  // Section S8 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t bivar_NL_mean_eq;
  double mu_x1;
  double mu_x2;
  double mu_x2_n;
  double sigma_x1;
  double sigma_x2;
  double sigma_x2_n;
  double rho_x1x2_n;
  double alpha_exp;
  double beta_exp;
  double x1_tol;
  double x2_tol;
  double s_c;

  // Means for the ith PDF component.
  mu_x1 = static_cast<double>(mu_chi_i);
  if ((beta_exp_in >= zero)) {
    mu_x2 = static_cast<double>(mu_y_i);
  } else {
    mu_x2 = static_cast<double>(((mu_y_i) > (y_tol) ? (mu_y_i) : (y_tol)));
  }
  mu_x2_n = static_cast<double>(mu_y_i_n);

  // Standard deviations for the ith PDF component.
  sigma_x1   = static_cast<double>(sigma_chi_i);
  sigma_x2   = static_cast<double>(sigma_y_i);
  sigma_x2_n = static_cast<double>(sigma_y_i_n);

  // Correlations for the ith PDF component.
  rho_x1x2_n = static_cast<double>(corr_chi_y_i_n);

  // Exponents.
  alpha_exp = static_cast<double>(alpha_exp_in);
  beta_exp  = static_cast<double>(beta_exp_in);

  // Tolerance values.
  // When the standard deviation of a variable is below the tolerance values,
  // it is considered to be zero, and the variable is considered to have a
  // constant value.
  x1_tol = static_cast<double>(chi_tol);
  x2_tol = static_cast<double>(y_tol);

  // Determine the value of the parabolic cylinder function input value, s_c.
  // The value s_c is being fed into the parabolic cylinder function.  When
  // the value of s_c is too large in magnitude (depending on the order of the
  // parabolic cylinder function), overflow occurs, and the output of the
  // parabolic cylinder function is +/-Inf.  This is primarily due to a large
  // ratio of mu_x1 to sigma_x1.  When the value of s_c is very large, the
  // distribution of x1 is basically a spike near the mean, so x1 is treated as
  // a constant.
  if ((sigma_x1 > x1_tol)) {
    s_c = ((mu_x1 / sigma_x1) + ((rho_x1x2_n * sigma_x2_n) * beta_exp));
  } else {
    // Note:  s_c is +inf when mu_x1 > 0 and sigma_x1 = 0, and s_c is -inf
    //        when mu_x1 < 0 and sigma_x1 = 0.  Furthermore, s_c is undefined
    //        when mu_x1 = 0 and sigma_x1 = 0.  However, within the context of
    //        this particular function, only the absolute value of s_c is
    //        relevant, and furthermore the absolute value of s_c is only
    //        relevant when sigma_x1 > 0.  Therefore, this statement only
    //        serves as divide-by-zero and compiler warning prevention.
    s_c = std::numeric_limits<double>::max();
  }

  // Based on the values of sigma_x1 (including the value of s_c compared to
  // parab_cyl_max_input) and sigma_x2, find the correct form of the bivariate
  // equation to use.

  if ((((sigma_x1 <= x1_tol) || (Kokkos::abs(s_c) > static_cast<double>(parab_cyl_max_input))) &&
       (sigma_x2 <= x2_tol))) {

    // The ith PDF component variance of both chi and y (r_r or N_cn) is 0.
    bivar_NL_mean_eq =
        static_cast<real_t>(bivar_NL_mean_const_all(mu_x1, mu_x2, alpha_exp, beta_exp, kk));

  } else {
    if (((sigma_x1 <= x1_tol) || (Kokkos::abs(s_c) > static_cast<double>(parab_cyl_max_input)))) {

      // The ith PDF component variance of chi is 0.
      bivar_NL_mean_eq = static_cast<real_t>(
          bivar_NL_mean_const_x1(mu_x1, mu_x2_n, sigma_x2_n, alpha_exp, beta_exp, kk));

    } else {
      if ((sigma_x2 <= x2_tol)) {

        // The ith PDF component variance of y (r_r or N_cn) is 0.
        bivar_NL_mean_eq = static_cast<real_t>(
            bivar_NL_mean_const_x2(mu_x1, mu_x2, sigma_x1, alpha_exp, beta_exp, kk));

      } else {

        // All fields vary in the ith PDF component.
        bivar_NL_mean_eq = static_cast<real_t>(bivar_NL_mean(mu_x1, mu_x2_n, sigma_x1, sigma_x2_n,
                                                             rho_x1x2_n, alpha_exp, beta_exp, kk));
      }
    }
  }

  return bivar_NL_mean_eq;
}
inline real_t bivar_LL_mean_eq(real_t mu_rr_i, real_t mu_Nr_i, real_t mu_rr_i_n, real_t mu_Nr_i_n,
                               real_t sigma_rr_i, real_t sigma_Nr_i, real_t sigma_rr_i_n,
                               real_t sigma_Nr_i_n, real_t corr_rr_Nr_i_n, real_t alpha_exp_in,
                               real_t beta_exp_in, const parameters_KK& kk) {
  // Description:
  // This function calculates the expression < x1^alpha x2^beta >, where
  // x1 = r_r and x2 = N_r.  The value of mean volume radius is given by:
  //
  // < KK_mvr > = KK_mvr_coef
  //              * ( mixt_frac < r_r^alpha N_r^beta (1) >
  //                  + ( 1 - mixt_frac ) < r_r^alpha N_r^beta (2) > ).
  //
  // One of four functions is called, based on whether either, both, or none of
  // x1 (r_r) and x2 (N_r) vary.  Each one of these four functions are the
  // result of an evaluated integral based on the specific situation.
  // References:
  // Section 4 of Larson, V. E. and B. M. Griffin, 2013:  Analytic upscaling of
  // a local microphysics scheme. Part I: Derivation.  Q. J. Roy. Meteorol.
  // Soc., 139, 670, 46--57, doi:http://dx.doi.org/10.1002/qj.1967.
  //
  // Section C.4 of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Section S1.4 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // A new subgrid-scale representation of hydrometeor fields using a
  // multivariate PDF.  Geosci. Model Dev., 9, 6,
  // doi:http://dx.doi.org/10.5194/gmd-9-2031-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t bivar_LL_mean_eq;
  double mu_x1;
  double mu_x2;
  double mu_x1_n;
  double mu_x2_n;
  double sigma_x1;
  double sigma_x2;
  double sigma_x1_n;
  double sigma_x2_n;
  double rho_x1x2_n;
  double alpha_exp;
  double beta_exp;
  double x1_tol;
  double x2_tol;

  // Means for the ith PDF component.
  if ((alpha_exp_in >= zero)) {
    mu_x1 = static_cast<double>(mu_rr_i);
  } else {
    mu_x1 = static_cast<double>(((mu_rr_i) > (rr_tol) ? (mu_rr_i) : (rr_tol)));
  }
  if ((beta_exp_in >= zero)) {
    mu_x2 = static_cast<double>(mu_Nr_i);
  } else {
    mu_x2 = static_cast<double>(((mu_Nr_i) > (Nr_tol) ? (mu_Nr_i) : (Nr_tol)));
  }
  mu_x1_n = static_cast<double>(mu_rr_i_n);
  mu_x2_n = static_cast<double>(mu_Nr_i_n);

  // Standard deviations for the ith PDF component.
  sigma_x1   = static_cast<double>(sigma_rr_i);
  sigma_x2   = static_cast<double>(sigma_Nr_i);
  sigma_x1_n = static_cast<double>(sigma_rr_i_n);
  sigma_x2_n = static_cast<double>(sigma_Nr_i_n);

  // Correlations for the ith PDF component.
  rho_x1x2_n = static_cast<double>(corr_rr_Nr_i_n);

  // Exponents.
  alpha_exp = static_cast<double>(alpha_exp_in);
  beta_exp  = static_cast<double>(beta_exp_in);

  // Tolerance values.
  // When the standard deviation of a variable is below the tolerance values,
  // it is considered to be zero, and the variable is considered to have a
  // constant value.
  x1_tol = static_cast<double>(rr_tol);
  x2_tol = static_cast<double>(Nr_tol);

  // Calculate the mean of the bivariate lognormal equation.
  if (((sigma_x1 <= x1_tol) && (sigma_x2 <= x2_tol))) {

    // The ith PDF component variance of both r_r and N_r is 0.
    bivar_LL_mean_eq =
        static_cast<real_t>(bivar_LL_mean_const_all(mu_x1, mu_x2, alpha_exp, beta_exp, kk));

  } else {
    if ((sigma_x1 <= x1_tol)) {

      // The ith PDF component variance of r_r is 0.
      bivar_LL_mean_eq = static_cast<real_t>(
          bivar_LL_mean_const_x1(mu_x1, mu_x2_n, sigma_x2_n, alpha_exp, beta_exp, kk));

    } else {
      if ((sigma_x2 <= x2_tol)) {

        // The ith PDF component variance of N_r is 0.
        bivar_LL_mean_eq = static_cast<real_t>(
            bivar_LL_mean_const_x1(mu_x2, mu_x1_n, sigma_x1_n, beta_exp, alpha_exp, kk));

      } else {

        // All fields vary in the ith PDF component.
        bivar_LL_mean_eq = static_cast<real_t>(bivar_LL_mean(
            mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, rho_x1x2_n, alpha_exp, beta_exp, kk));
      }
    }
  }

  return bivar_LL_mean_eq;
}
} // namespace clubb_kokkos::kk_engine
