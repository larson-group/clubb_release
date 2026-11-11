// Native KK_upscaled_covariances numerical formulas, evaluated in source arithmetic order.
// Appended caller-owned KK parameters replace implicit module state.
#pragma once
#include "Microphys/KK_microphys/engine_contracts.hpp"
#include "Microphys/KK_microphys/KK_utilities.hpp"
namespace clubb_kokkos::kk_engine {
inline void KK_upscaled_covar_driver(
    real_t w_mean, real_t rtm, real_t thlm, real_t exner, real_t rrainm, real_t Nrm, real_t mu_w_1,
    real_t mu_w_2, real_t mu_chi_1, real_t mu_chi_2, real_t mu_eta_1, real_t mu_eta_2,
    real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_Ncn_1,
    real_t mu_Ncn_2, real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n, real_t mu_Nr_2_n,
    real_t mu_Ncn_1_n, real_t mu_Ncn_2_n, real_t sigma_w_1, real_t sigma_w_2, real_t sigma_chi_1,
    real_t sigma_chi_2, real_t sigma_eta_1, real_t sigma_eta_2, real_t sigma_rr_1,
    real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_Ncn_1, real_t sigma_Ncn_2,
    real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
    real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n, real_t corr_w_chi_1, real_t corr_w_chi_2,
    real_t corr_w_rr_1_n, real_t corr_w_rr_2_n, real_t corr_w_Nr_1_n, real_t corr_w_Nr_2_n,
    real_t corr_w_Ncn_1_n, real_t corr_w_Ncn_2_n, real_t corr_chi_eta_1, real_t corr_chi_eta_2,
    real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n, real_t corr_chi_Nr_1_n, real_t corr_chi_Nr_2_n,
    real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n, real_t corr_eta_rr_1_n,
    real_t corr_eta_rr_2_n, real_t corr_eta_Nr_1_n, real_t corr_eta_Nr_2_n, real_t corr_eta_Ncn_1_n,
    real_t corr_eta_Ncn_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t mixt_frac,
    real_t precip_frac_1, real_t precip_frac_2, real_t KK_evap_coef, real_t KK_auto_coef,
    real_t KK_accr_coef, real_t KK_evap_tndcy, real_t KK_auto_tndcy, real_t KK_accr_tndcy,
    real_t mu_rt_1, real_t mu_rt_2, real_t mu_thl_1, real_t mu_thl_2, real_t crt1, real_t crt2,
    real_t cthl1, real_t cthl2, real_t& wprtp_mc_src_tndcy, real_t& wpthlp_mc_src_tndcy,
    real_t& rtp2_mc_src_tndcy, real_t& thlp2_mc_src_tndcy, real_t& rtpthlp_mc_src_tndcy,
    real_t& w_KK_evap_covar, real_t& rt_KK_evap_covar, real_t& thl_KK_evap_covar,
    real_t& w_KK_auto_covar, real_t& rt_KK_auto_covar, real_t& thl_KK_auto_covar,
    real_t& w_KK_accr_covar, real_t& rt_KK_accr_covar, real_t& thl_KK_accr_covar,
    const parameters_KK& kk) {
  // Description:
  // References:
  // Griffin, B. M., 2016:  Improving the Subgrid-Scale Representation of
  //    Hydrometeors and Microphysical Feedback Effects Using a Multivariate
  //    PDF.  Doctoral dissertation, University of Wisconsin -- Milwaukee,
  //    Milwaukee, WI, Paper 1144, 165 pp., URL
  //    http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Griffin, B. M. and V. E. Larson, 2016:  Parameterizing microphysical
  //    effects on variances and covariances of moisture and heat content using
  //    a multivariate probability density function: a study with CLUBB (tag
  //    MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  //    doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  //-----------------------------------------------------------------------
  // Input Variables
  // Output Variables

  // Calculate the covariance of vertical velocity and KK evaporation tendency.
  if (((rrainm > rr_tol) && (Nrm > Nr_tol))) {

    w_KK_evap_covar = covar_x_KK_evap(
        mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n,
        mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_w_1, sigma_w_2, sigma_chi_1, sigma_chi_2, sigma_rr_1,
        sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
        corr_w_chi_1, corr_w_chi_2, corr_w_rr_1_n, corr_w_rr_2_n, corr_w_Nr_1_n, corr_w_Nr_2_n,
        corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n, corr_chi_Nr_2_n, corr_rr_Nr_1_n,
        corr_rr_Nr_2_n, w_mean, KK_evap_tndcy, KK_evap_coef, w_tol, mixt_frac, precip_frac_1,
        precip_frac_2, kk);

  } else {

    w_KK_evap_covar = zero;
  }

  // Calculate the covariance of total water mixing ratio and KK evaporation
  // tendency.
  if (((rrainm > rr_tol) && (Nrm > Nr_tol))) {

    rt_KK_evap_covar = covar_rt_KK_evap(
        mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n,
        mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2,
        sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n,
        sigma_Nr_2_n, corr_chi_eta_1, corr_chi_eta_2, corr_eta_rr_1_n, corr_eta_rr_2_n,
        corr_eta_Nr_1_n, corr_eta_Nr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n,
        corr_chi_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, mixt_frac, precip_frac_1, precip_frac_2,
        rtm, mu_rt_1, mu_rt_2, KK_evap_tndcy, KK_evap_coef, eta_tol, crt1, crt2, kk);

  } else {

    rt_KK_evap_covar = zero;
  }

  // Calculate the covariance of liquid water potential temperature and
  // KK evaporation tendency.
  if (((rrainm > rr_tol) && (Nrm > Nr_tol))) {

    thl_KK_evap_covar = covar_thl_KK_evap(
        mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n,
        mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2,
        sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n,
        sigma_Nr_2_n, corr_chi_eta_1, corr_chi_eta_2, corr_eta_rr_1_n, corr_eta_rr_2_n,
        corr_eta_Nr_1_n, corr_eta_Nr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n, corr_chi_Nr_1_n,
        corr_chi_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, mixt_frac, precip_frac_1, precip_frac_2,
        thlm, mu_thl_1, mu_thl_2, KK_evap_tndcy, KK_evap_coef, eta_tol, cthl1, cthl2, kk);

  } else {

    thl_KK_evap_covar = zero;
  }

  // Calculate the covariance of vertical velocity and KK autoconversion
  // tendency.
  w_KK_auto_covar = covar_x_KK_auto(
      mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n, sigma_w_1,
      sigma_w_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n,
      corr_w_chi_1, corr_w_chi_2, corr_w_Ncn_1_n, corr_w_Ncn_2_n, corr_chi_Ncn_1_n,
      corr_chi_Ncn_2_n, w_mean, KK_auto_tndcy, KK_auto_coef, w_tol, mixt_frac, kk);

  // Calculate the covariance of total water mixing ratio and KK autoconversion
  // tendency.
  rt_KK_auto_covar = covar_rt_KK_auto(
      mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n,
      sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2, sigma_Ncn_1_n,
      sigma_Ncn_2_n, corr_chi_eta_1, corr_chi_eta_2, corr_eta_Ncn_1_n, corr_eta_Ncn_2_n,
      corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, rtm, mu_rt_1, mu_rt_2, KK_auto_tndcy, KK_auto_coef,
      eta_tol, crt1, crt2, mixt_frac, kk);

  // Calculate the covariance of liquid water potential temperature and
  // KK autoconversion tendency.
  thl_KK_auto_covar = covar_thl_KK_auto(
      mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, mu_Ncn_1_n, mu_Ncn_2_n,
      sigma_eta_1, sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_Ncn_1, sigma_Ncn_2, sigma_Ncn_1_n,
      sigma_Ncn_2_n, corr_chi_eta_1, corr_chi_eta_2, corr_eta_Ncn_1_n, corr_eta_Ncn_2_n,
      corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, thlm, mu_thl_1, mu_thl_2, KK_auto_tndcy, KK_auto_coef,
      eta_tol, cthl1, cthl2, mixt_frac, kk);

  // Calculate the covariance of vertical velocity and KK accretion tendency.
  if ((rrainm > rr_tol)) {

    w_KK_accr_covar = covar_x_KK_accr(
        mu_w_1, mu_w_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_w_1,
        sigma_w_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n,
        corr_w_chi_1, corr_w_chi_2, corr_w_rr_1_n, corr_w_rr_2_n, corr_chi_rr_1_n, corr_chi_rr_2_n,
        w_mean, KK_accr_tndcy, KK_accr_coef, w_tol, mixt_frac, precip_frac_1, precip_frac_2, kk);

  } else {

    w_KK_accr_covar = zero;
  }

  // Calculate the covariance of total water mixing ratio and KK accretion
  // tendency.
  if ((rrainm > rr_tol)) {

    rt_KK_accr_covar = covar_rt_KK_accr(
        mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_eta_1,
        sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n,
        corr_chi_eta_1, corr_chi_eta_2, corr_eta_rr_1_n, corr_eta_rr_2_n, corr_chi_rr_1_n,
        corr_chi_rr_2_n, rtm, mu_rt_1, mu_rt_2, KK_accr_tndcy, KK_accr_coef, eta_tol, crt1, crt2,
        mixt_frac, precip_frac_1, precip_frac_2, kk);

  } else {

    rt_KK_accr_covar = zero;
  }

  // Calculate the covariance of liquid water potential temperature and
  // KK accretion tendency.
  if ((rrainm > rr_tol)) {

    thl_KK_accr_covar = covar_thl_KK_accr(
        mu_eta_1, mu_eta_2, mu_chi_1, mu_chi_2, mu_rr_1, mu_rr_2, mu_rr_1_n, mu_rr_2_n, sigma_eta_1,
        sigma_eta_2, sigma_chi_1, sigma_chi_2, sigma_rr_1, sigma_rr_2, sigma_rr_1_n, sigma_rr_2_n,
        corr_chi_eta_1, corr_chi_eta_2, corr_eta_rr_1_n, corr_eta_rr_2_n, corr_chi_rr_1_n,
        corr_chi_rr_2_n, thlm, mu_thl_1, mu_thl_2, KK_accr_tndcy, KK_accr_coef, eta_tol, cthl1,
        cthl2, mixt_frac, precip_frac_1, precip_frac_2, kk);

  } else {

    thl_KK_accr_covar = zero;
  }

  // Calculate the microphysics tendency for <w'r_t'>.
  wprtp_mc_src_tndcy = (-((w_KK_auto_covar + w_KK_accr_covar) + w_KK_evap_covar));

  // Calculate the microphysics tendency for <w'th_l'>.
  wpthlp_mc_src_tndcy =
      ((Lv / (Cp * exner)) * ((w_KK_auto_covar + w_KK_accr_covar) + w_KK_evap_covar));

  // Calculate the microphysics tendency for <r_t'^2>.
  rtp2_mc_src_tndcy = ((-two) * ((rt_KK_auto_covar + rt_KK_accr_covar) + rt_KK_evap_covar));

  // Calculate the microphysics tendency for <th_l'^2>.
  thlp2_mc_src_tndcy =
      ((two * (Lv / (Cp * exner))) * ((thl_KK_auto_covar + thl_KK_accr_covar) + thl_KK_evap_covar));

  // Calculate the microphysics tendency for <r_t'th_l'>.
  rtpthlp_mc_src_tndcy =
      (((Lv / (Cp * exner)) * ((rt_KK_auto_covar + rt_KK_accr_covar) + rt_KK_evap_covar)) -
       ((thl_KK_auto_covar + thl_KK_accr_covar) + thl_KK_evap_covar));

  return;
}
inline real_t covar_x_KK_evap(
    real_t mu_x_1, real_t mu_x_2, real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1, real_t mu_rr_2,
    real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
    real_t mu_Nr_2_n, real_t sigma_x_1, real_t sigma_x_2, real_t sigma_chi_1, real_t sigma_chi_2,
    real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
    real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n, real_t corr_x_chi_1,
    real_t corr_x_chi_2, real_t corr_x_rr_1_n, real_t corr_x_rr_2_n, real_t corr_x_Nr_1_n,
    real_t corr_x_Nr_2_n, real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n, real_t corr_chi_Nr_1_n,
    real_t corr_chi_Nr_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t x_mean,
    real_t KK_evap_tndcy, real_t KK_evap_coef, real_t x_tol, real_t mixt_frac, real_t precip_frac_1,
    real_t precip_frac_2, const parameters_KK& kk) {
  // Description:
  // This function calculates the covariance of x and KK evaporation
  // tendency, which can be written as < x'((dr_r/dt)_KKevap)' >, or more
  // simply as < x'KK_evap' >.
  // References:
  // Eq. (E21) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (A28) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
  // microphysical effects on variances and covariances of moisture and heat
  // content using a multivariate probability density function: a study with
  // CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_x_KK_evap;
  real_t alpha_exp;
  real_t beta_exp;
  real_t gamma_exp;

  // Values of the KK exponents.
  alpha_exp = KK_evap_Supersat_exp;
  beta_exp  = KK_evap_rr_exp;
  gamma_exp = KK_evap_Nr_exp;

  // Calculate the covariance of x and KK evaporation tendency.
  covar_x_KK_evap =
      ((mixt_frac *
        (((KK_evap_coef * precip_frac_1) *
          quadrivar_NNLL_covar_eq(
              mu_x_1, mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_x_1, sigma_chi_1,
              sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_x_chi_1, corr_x_rr_1_n,
              corr_x_Nr_1_n, corr_chi_rr_1_n, corr_chi_Nr_1_n, corr_rr_Nr_1_n, x_mean,
              KK_evap_tndcy, KK_evap_coef, x_tol, alpha_exp, beta_exp, gamma_exp, kk)) -
         (((one - precip_frac_1) * (mu_x_1 - x_mean)) * KK_evap_tndcy))) +
       ((one - mixt_frac) *
        (((KK_evap_coef * precip_frac_2) *
          quadrivar_NNLL_covar_eq(
              mu_x_2, mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_x_2, sigma_chi_2,
              sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_x_chi_2, corr_x_rr_2_n,
              corr_x_Nr_2_n, corr_chi_rr_2_n, corr_chi_Nr_2_n, corr_rr_Nr_2_n, x_mean,
              KK_evap_tndcy, KK_evap_coef, x_tol, alpha_exp, beta_exp, gamma_exp, kk)) -
         (((one - precip_frac_2) * (mu_x_2 - x_mean)) * KK_evap_tndcy))));

  return covar_x_KK_evap;
}
inline real_t covar_rt_KK_evap(
    real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1,
    real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
    real_t mu_Nr_1_n, real_t mu_Nr_2_n, real_t sigma_eta_1, real_t sigma_eta_2, real_t sigma_chi_1,
    real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2,
    real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
    real_t corr_eta_chi_1, real_t corr_eta_chi_2, real_t corr_eta_rr_1_n, real_t corr_eta_rr_2_n,
    real_t corr_eta_Nr_1_n, real_t corr_eta_Nr_2_n, real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n,
    real_t corr_chi_Nr_1_n, real_t corr_chi_Nr_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
    real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2, real_t rtm, real_t mu_rt_1,
    real_t mu_rt_2, real_t KK_evap_tndcy, real_t KK_evap_coef, real_t eta_tol, real_t crt1,
    real_t crt2, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (E24) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (A31) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
  // microphysical effects on variances and covariances of moisture and heat
  // content using a multivariate probability density function: a study with
  // CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_rt_KK_evap;
  real_t alpha_exp;
  real_t beta_exp;
  real_t gamma_exp;
  real_t comp_1_contrib;
  real_t comp_2_contrib;

  // Values of the KK exponents.
  alpha_exp = KK_evap_Supersat_exp;
  beta_exp  = KK_evap_rr_exp;
  gamma_exp = KK_evap_Nr_exp;

  // Calculate the contribution from PDF component 1 to the covariance of
  // r_t and KK evaporation tendency.
  comp_1_contrib =
      ((KK_evap_coef * precip_frac_1) *
       (((one / (two * crt1)) *
         (quadrivar_NNLL_covar_eq(
              mu_eta_1, mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_eta_1, sigma_chi_1,
              sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_eta_chi_1, corr_eta_rr_1_n,
              corr_eta_Nr_1_n, corr_chi_rr_1_n, corr_chi_Nr_1_n, corr_rr_Nr_1_n, mu_eta_1,
              KK_evap_tndcy, KK_evap_coef, eta_tol, alpha_exp, beta_exp, gamma_exp, kk) +
          trivar_NLL_mean_eq(mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_chi_1,
                             sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_chi_rr_1_n,
                             corr_chi_Nr_1_n, corr_rr_Nr_1_n, (alpha_exp + one), beta_exp,
                             gamma_exp, kk))) +
        (((mu_rt_1 - rtm) - (mu_chi_1 / (two * crt1))) *
         trivar_NLL_mean_eq(mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_chi_1,
                            sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_chi_rr_1_n,
                            corr_chi_Nr_1_n, corr_rr_Nr_1_n, alpha_exp, beta_exp, gamma_exp, kk))));

  // Calculate the contribution from PDF component 2 to the covariance of
  // r_t and KK evaporation tendency.
  comp_2_contrib =
      ((KK_evap_coef * precip_frac_2) *
       (((one / (two * crt2)) *
         (quadrivar_NNLL_covar_eq(
              mu_eta_2, mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_eta_2, sigma_chi_2,
              sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_eta_chi_2, corr_eta_rr_2_n,
              corr_eta_Nr_2_n, corr_chi_rr_2_n, corr_chi_Nr_2_n, corr_rr_Nr_2_n, mu_eta_2,
              KK_evap_tndcy, KK_evap_coef, eta_tol, alpha_exp, beta_exp, gamma_exp, kk) +
          trivar_NLL_mean_eq(mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_chi_2,
                             sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_chi_rr_2_n,
                             corr_chi_Nr_2_n, corr_rr_Nr_2_n, (alpha_exp + one), beta_exp,
                             gamma_exp, kk))) +
        (((mu_rt_2 - rtm) - (mu_chi_2 / (two * crt2))) *
         trivar_NLL_mean_eq(mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_chi_2,
                            sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_chi_rr_2_n,
                            corr_chi_Nr_2_n, corr_rr_Nr_2_n, alpha_exp, beta_exp, gamma_exp, kk))));

  // Calculate the covariance of r_t and KK evaporation tendency.
  covar_rt_KK_evap = ((mixt_frac * comp_1_contrib) + ((one - mixt_frac) * comp_2_contrib));

  return covar_rt_KK_evap;
}
inline real_t covar_thl_KK_evap(
    real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2, real_t mu_rr_1,
    real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
    real_t mu_Nr_1_n, real_t mu_Nr_2_n, real_t sigma_eta_1, real_t sigma_eta_2, real_t sigma_chi_1,
    real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2,
    real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
    real_t corr_eta_chi_1, real_t corr_eta_chi_2, real_t corr_eta_rr_1_n, real_t corr_eta_rr_2_n,
    real_t corr_eta_Nr_1_n, real_t corr_eta_Nr_2_n, real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n,
    real_t corr_chi_Nr_1_n, real_t corr_chi_Nr_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
    real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2, real_t thlm, real_t mu_thl_1,
    real_t mu_thl_2, real_t KK_evap_tndcy, real_t KK_evap_coef, real_t eta_tol, real_t cthl1,
    real_t cthl2, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (E27) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (A34) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
  // microphysical effects on variances and covariances of moisture and heat
  // content using a multivariate probability density function: a study with
  // CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_thl_KK_evap;
  real_t alpha_exp;
  real_t beta_exp;
  real_t gamma_exp;
  real_t comp_1_contrib;
  real_t comp_2_contrib;

  // Values of the KK exponents.
  alpha_exp = KK_evap_Supersat_exp;
  beta_exp  = KK_evap_rr_exp;
  gamma_exp = KK_evap_Nr_exp;

  // Calculate the contribution from PDF component 1 to the covariance of
  // th_l and KK evaporation tendency.
  comp_1_contrib =
      ((KK_evap_coef * precip_frac_1) *
       (((one / (two * cthl1)) *
         (quadrivar_NNLL_covar_eq(
              mu_eta_1, mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_eta_1, sigma_chi_1,
              sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_eta_chi_1, corr_eta_rr_1_n,
              corr_eta_Nr_1_n, corr_chi_rr_1_n, corr_chi_Nr_1_n, corr_rr_Nr_1_n, mu_eta_1,
              KK_evap_tndcy, KK_evap_coef, eta_tol, alpha_exp, beta_exp, gamma_exp, kk) -
          trivar_NLL_mean_eq(mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_chi_1,
                             sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_chi_rr_1_n,
                             corr_chi_Nr_1_n, corr_rr_Nr_1_n, (alpha_exp + one), beta_exp,
                             gamma_exp, kk))) +
        (((mu_thl_1 - thlm) + (mu_chi_1 / (two * cthl1))) *
         trivar_NLL_mean_eq(mu_chi_1, mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_chi_1,
                            sigma_rr_1, sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_chi_rr_1_n,
                            corr_chi_Nr_1_n, corr_rr_Nr_1_n, alpha_exp, beta_exp, gamma_exp, kk))));

  // Calculate the contribution from PDF component 2 to the covariance of
  // th_l and KK evaporation tendency.
  comp_2_contrib =
      ((KK_evap_coef * precip_frac_2) *
       (((one / (two * cthl2)) *
         (quadrivar_NNLL_covar_eq(
              mu_eta_2, mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_eta_2, sigma_chi_2,
              sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_eta_chi_2, corr_eta_rr_2_n,
              corr_eta_Nr_2_n, corr_chi_rr_2_n, corr_chi_Nr_2_n, corr_rr_Nr_2_n, mu_eta_2,
              KK_evap_tndcy, KK_evap_coef, eta_tol, alpha_exp, beta_exp, gamma_exp, kk) -
          trivar_NLL_mean_eq(mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_chi_2,
                             sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_chi_rr_2_n,
                             corr_chi_Nr_2_n, corr_rr_Nr_2_n, (alpha_exp + one), beta_exp,
                             gamma_exp, kk))) +
        (((mu_thl_2 - thlm) + (mu_chi_2 / (two * cthl2))) *
         trivar_NLL_mean_eq(mu_chi_2, mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_chi_2,
                            sigma_rr_2, sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_chi_rr_2_n,
                            corr_chi_Nr_2_n, corr_rr_Nr_2_n, alpha_exp, beta_exp, gamma_exp, kk))));

  // Calculate the covariance of th_l and KK evaporation tendency.
  covar_thl_KK_evap = ((mixt_frac * comp_1_contrib) + ((one - mixt_frac) * comp_2_contrib));

  return covar_thl_KK_evap;
}
inline real_t covar_x_KK_auto(real_t mu_x_1, real_t mu_x_2, real_t mu_chi_1, real_t mu_chi_2,
                              real_t mu_Ncn_1, real_t mu_Ncn_2, real_t mu_Ncn_1_n,
                              real_t mu_Ncn_2_n, real_t sigma_x_1, real_t sigma_x_2,
                              real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_Ncn_1,
                              real_t sigma_Ncn_2, real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n,
                              real_t corr_x_chi_1, real_t corr_x_chi_2, real_t corr_x_Ncn_1_n,
                              real_t corr_x_Ncn_2_n, real_t corr_chi_Ncn_1_n,
                              real_t corr_chi_Ncn_2_n, real_t x_mean, real_t KK_auto_tndcy,
                              real_t KK_auto_coef, real_t x_tol, real_t mixt_frac,
                              const parameters_KK& kk) {
  // Description:
  // This function calculates the covariance of x and KK autoconversion
  // tendency, which can be written as < x'((dr_r/dt)_KKauto)' >, or more
  // simply as < x'KK_auto' >.
  // References:
  // Eq. (E3) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (A9) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
  // microphysical effects on variances and covariances of moisture and heat
  // content using a multivariate probability density function: a study with
  // CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_x_KK_auto;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_auto_rc_exp;
  beta_exp  = KK_auto_Nc_exp;

  // Calculate the covariance of x and KK autoconversion tendency.
  covar_x_KK_auto =
      (KK_auto_coef *
       ((mixt_frac * trivar_NNL_covar_eq(mu_x_1, mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_x_1,
                                         sigma_chi_1, sigma_Ncn_1, sigma_Ncn_1_n, corr_x_chi_1,
                                         corr_x_Ncn_1_n, corr_chi_Ncn_1_n, x_mean, KK_auto_tndcy,
                                         KK_auto_coef, x_tol, Ncn_tol, alpha_exp, beta_exp, kk)) +
        ((one - mixt_frac) * trivar_NNL_covar_eq(mu_x_2, mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_x_2,
                                                 sigma_chi_2, sigma_Ncn_2, sigma_Ncn_2_n,
                                                 corr_x_chi_2, corr_x_Ncn_2_n, corr_chi_Ncn_2_n,
                                                 x_mean, KK_auto_tndcy, KK_auto_coef, x_tol,
                                                 Ncn_tol, alpha_exp, beta_exp, kk))));

  return covar_x_KK_auto;
}
inline real_t covar_rt_KK_auto(real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2,
                               real_t mu_Ncn_1, real_t mu_Ncn_2, real_t mu_Ncn_1_n,
                               real_t mu_Ncn_2_n, real_t sigma_eta_1, real_t sigma_eta_2,
                               real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_Ncn_1,
                               real_t sigma_Ncn_2, real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n,
                               real_t corr_eta_chi_1, real_t corr_eta_chi_2,
                               real_t corr_eta_Ncn_1_n, real_t corr_eta_Ncn_2_n,
                               real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n, real_t rtm,
                               real_t mu_rt_1, real_t mu_rt_2, real_t KK_auto_tndcy,
                               real_t KK_auto_coef, real_t eta_tol, real_t crt1, real_t crt2,
                               real_t mixt_frac, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (E6) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (A12) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
  // microphysical effects on variances and covariances of moisture and heat
  // content using a multivariate probability density function: a study with
  // CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_rt_KK_auto;
  real_t alpha_exp;
  real_t beta_exp;
  real_t comp_1_contrib;
  real_t comp_2_contrib;

  // Values of the KK exponents.
  alpha_exp = KK_auto_rc_exp;
  beta_exp  = KK_auto_Nc_exp;

  // Calculate the contribution from PDF component 1 to the covariance of
  // r_t and KK autoconversion tendency.
  comp_1_contrib =
      (KK_auto_coef *
       (((one / (two * crt1)) *
         (trivar_NNL_covar_eq(mu_eta_1, mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_eta_1, sigma_chi_1,
                              sigma_Ncn_1, sigma_Ncn_1_n, corr_eta_chi_1, corr_eta_Ncn_1_n,
                              corr_chi_Ncn_1_n, mu_eta_1, KK_auto_tndcy, KK_auto_coef, eta_tol,
                              Ncn_tol, alpha_exp, beta_exp, kk) +
          bivar_NL_mean_eq(mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_chi_1, sigma_Ncn_1, sigma_Ncn_1_n,
                           corr_chi_Ncn_1_n, Nc_tol, (alpha_exp + one), beta_exp, kk))) +
        (((mu_rt_1 - rtm) - (mu_chi_1 / (two * crt1))) *
         bivar_NL_mean_eq(mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_chi_1, sigma_Ncn_1, sigma_Ncn_1_n,
                          corr_chi_Ncn_1_n, Nc_tol, alpha_exp, beta_exp, kk))));

  // Calculate the contribution from PDF component 2 to the covariance of
  // r_t and KK autoconversion tendency.
  comp_2_contrib =
      (KK_auto_coef *
       (((one / (two * crt2)) *
         (trivar_NNL_covar_eq(mu_eta_2, mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_eta_2, sigma_chi_2,
                              sigma_Ncn_2, sigma_Ncn_2_n, corr_eta_chi_2, corr_eta_Ncn_2_n,
                              corr_chi_Ncn_2_n, mu_eta_2, KK_auto_tndcy, KK_auto_coef, eta_tol,
                              Ncn_tol, alpha_exp, beta_exp, kk) +
          bivar_NL_mean_eq(mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_chi_2, sigma_Ncn_2, sigma_Ncn_2_n,
                           corr_chi_Ncn_2_n, Nc_tol, (alpha_exp + one), beta_exp, kk))) +
        (((mu_rt_2 - rtm) - (mu_chi_2 / (two * crt2))) *
         bivar_NL_mean_eq(mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_chi_2, sigma_Ncn_2, sigma_Ncn_2_n,
                          corr_chi_Ncn_2_n, Nc_tol, alpha_exp, beta_exp, kk))));

  // Calculate the covariance of r_t and KK autoconversion tendency.
  covar_rt_KK_auto = ((mixt_frac * comp_1_contrib) + ((one - mixt_frac) * comp_2_contrib));

  return covar_rt_KK_auto;
}
inline real_t covar_thl_KK_auto(real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2,
                                real_t mu_Ncn_1, real_t mu_Ncn_2, real_t mu_Ncn_1_n,
                                real_t mu_Ncn_2_n, real_t sigma_eta_1, real_t sigma_eta_2,
                                real_t sigma_chi_1, real_t sigma_chi_2, real_t sigma_Ncn_1,
                                real_t sigma_Ncn_2, real_t sigma_Ncn_1_n, real_t sigma_Ncn_2_n,
                                real_t corr_eta_chi_1, real_t corr_eta_chi_2,
                                real_t corr_eta_Ncn_1_n, real_t corr_eta_Ncn_2_n,
                                real_t corr_chi_Ncn_1_n, real_t corr_chi_Ncn_2_n, real_t thlm,
                                real_t mu_thl_1, real_t mu_thl_2, real_t KK_auto_tndcy,
                                real_t KK_auto_coef, real_t eta_tol, real_t cthl1, real_t cthl2,
                                real_t mixt_frac, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (E9) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (A15) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
  // microphysical effects on variances and covariances of moisture and heat
  // content using a multivariate probability density function: a study with
  // CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_thl_KK_auto;
  real_t alpha_exp;
  real_t beta_exp;
  real_t comp_1_contrib;
  real_t comp_2_contrib;

  // Values of the KK exponents.
  alpha_exp = KK_auto_rc_exp;
  beta_exp  = KK_auto_Nc_exp;

  // Calculate the contribution from PDF component 1 to the covariance of
  // th_l and KK autoconversion tendency.
  comp_1_contrib =
      (KK_auto_coef *
       (((one / (two * cthl1)) *
         (trivar_NNL_covar_eq(mu_eta_1, mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_eta_1, sigma_chi_1,
                              sigma_Ncn_1, sigma_Ncn_1_n, corr_eta_chi_1, corr_eta_Ncn_1_n,
                              corr_chi_Ncn_1_n, mu_eta_1, KK_auto_tndcy, KK_auto_coef, eta_tol,
                              Ncn_tol, alpha_exp, beta_exp, kk) -
          bivar_NL_mean_eq(mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_chi_1, sigma_Ncn_1, sigma_Ncn_1_n,
                           corr_chi_Ncn_1_n, Nc_tol, (alpha_exp + one), beta_exp, kk))) +
        (((mu_thl_1 - thlm) + (mu_chi_1 / (two * cthl1))) *
         bivar_NL_mean_eq(mu_chi_1, mu_Ncn_1, mu_Ncn_1_n, sigma_chi_1, sigma_Ncn_1, sigma_Ncn_1_n,
                          corr_chi_Ncn_1_n, Nc_tol, alpha_exp, beta_exp, kk))));

  // Calculate the contribution from PDF component 2 to the covariance of
  // th_l and KK autoconversion tendency.
  comp_2_contrib =
      (KK_auto_coef *
       (((one / (two * cthl2)) *
         (trivar_NNL_covar_eq(mu_eta_2, mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_eta_2, sigma_chi_2,
                              sigma_Ncn_2, sigma_Ncn_2_n, corr_eta_chi_2, corr_eta_Ncn_2_n,
                              corr_chi_Ncn_2_n, mu_eta_2, KK_auto_tndcy, KK_auto_coef, eta_tol,
                              Ncn_tol, alpha_exp, beta_exp, kk) -
          bivar_NL_mean_eq(mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_chi_2, sigma_Ncn_2, sigma_Ncn_2_n,
                           corr_chi_Ncn_2_n, Nc_tol, (alpha_exp + one), beta_exp, kk))) +
        (((mu_thl_2 - thlm) + (mu_chi_2 / (two * cthl2))) *
         bivar_NL_mean_eq(mu_chi_2, mu_Ncn_2, mu_Ncn_2_n, sigma_chi_2, sigma_Ncn_2, sigma_Ncn_2_n,
                          corr_chi_Ncn_2_n, Nc_tol, alpha_exp, beta_exp, kk))));

  // Calculate the covariance of th_l and KK autoconversion tendency.
  covar_thl_KK_auto = ((mixt_frac * comp_1_contrib) + ((one - mixt_frac) * comp_2_contrib));

  return covar_thl_KK_auto;
}
inline real_t covar_x_KK_accr(real_t mu_x_1, real_t mu_x_2, real_t mu_chi_1, real_t mu_chi_2,
                              real_t mu_rr_1, real_t mu_rr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
                              real_t sigma_x_1, real_t sigma_x_2, real_t sigma_chi_1,
                              real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2,
                              real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t corr_x_chi_1,
                              real_t corr_x_chi_2, real_t corr_x_rr_1_n, real_t corr_x_rr_2_n,
                              real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n, real_t x_mean,
                              real_t KK_accr_tndcy, real_t KK_accr_coef, real_t x_tol,
                              real_t mixt_frac, real_t precip_frac_1, real_t precip_frac_2,
                              const parameters_KK& kk) {
  // Description:
  // This function calculates the covariance of x and KK accretion
  // tendency, which can be written as < x'((dr_r/dt)_KKaccr)' >, or more
  // simply as < x'KK_accr' >.
  // References:
  // Eq. (E12) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (A18) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
  // microphysical effects on variances and covariances of moisture and heat
  // content using a multivariate probability density function: a study with
  // CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_x_KK_accr;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_accr_rc_exp;
  beta_exp  = KK_accr_rr_exp;

  // Calculate the covariance of x and KK accretion tendency.
  covar_x_KK_accr =
      ((mixt_frac * (((KK_accr_coef * precip_frac_1) *
                      trivar_NNL_covar_eq(mu_x_1, mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_x_1,
                                          sigma_chi_1, sigma_rr_1, sigma_rr_1_n, corr_x_chi_1,
                                          corr_x_rr_1_n, corr_chi_rr_1_n, x_mean, KK_accr_tndcy,
                                          KK_accr_coef, x_tol, rr_tol, alpha_exp, beta_exp, kk)) -
                     (((one - precip_frac_1) * (mu_x_1 - x_mean)) * KK_accr_tndcy))) +
       ((one - mixt_frac) *
        (((KK_accr_coef * precip_frac_2) *
          trivar_NNL_covar_eq(mu_x_2, mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_x_2, sigma_chi_2,
                              sigma_rr_2, sigma_rr_2_n, corr_x_chi_2, corr_x_rr_2_n,
                              corr_chi_rr_2_n, x_mean, KK_accr_tndcy, KK_accr_coef, x_tol, rr_tol,
                              alpha_exp, beta_exp, kk)) -
         (((one - precip_frac_2) * (mu_x_2 - x_mean)) * KK_accr_tndcy))));

  return covar_x_KK_accr;
}
inline real_t covar_rt_KK_accr(real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2,
                               real_t mu_rr_1, real_t mu_rr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
                               real_t sigma_eta_1, real_t sigma_eta_2, real_t sigma_chi_1,
                               real_t sigma_chi_2, real_t sigma_rr_1, real_t sigma_rr_2,
                               real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t corr_eta_chi_1,
                               real_t corr_eta_chi_2, real_t corr_eta_rr_1_n,
                               real_t corr_eta_rr_2_n, real_t corr_chi_rr_1_n,
                               real_t corr_chi_rr_2_n, real_t rtm, real_t mu_rt_1, real_t mu_rt_2,
                               real_t KK_accr_tndcy, real_t KK_accr_coef, real_t eta_tol,
                               real_t crt1, real_t crt2, real_t mixt_frac, real_t precip_frac_1,
                               real_t precip_frac_2, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (E15) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (A21) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
  // microphysical effects on variances and covariances of moisture and heat
  // content using a multivariate probability density function: a study with
  // CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_rt_KK_accr;
  real_t alpha_exp;
  real_t beta_exp;
  real_t comp_1_contrib;
  real_t comp_2_contrib;

  // Values of the KK exponents.
  alpha_exp = KK_accr_rc_exp;
  beta_exp  = KK_accr_rr_exp;

  // Calculate the contribution from PDF component 1 to the covariance of
  // r_t and KK accretion tendency.
  comp_1_contrib =
      ((KK_accr_coef * precip_frac_1) *
       (((one / (two * crt1)) *
         (trivar_NNL_covar_eq(mu_eta_1, mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_eta_1, sigma_chi_1,
                              sigma_rr_1, sigma_rr_1_n, corr_eta_chi_1, corr_eta_rr_1_n,
                              corr_chi_rr_1_n, mu_eta_1, KK_accr_tndcy, KK_accr_coef, eta_tol,
                              rr_tol, alpha_exp, beta_exp, kk) +
          bivar_NL_mean_eq(mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_chi_1, sigma_rr_1, sigma_rr_1_n,
                           corr_chi_rr_1_n, rr_tol, (alpha_exp + one), beta_exp, kk))) +
        (((mu_rt_1 - rtm) - (mu_chi_1 / (two * crt1))) *
         bivar_NL_mean_eq(mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_chi_1, sigma_rr_1, sigma_rr_1_n,
                          corr_chi_rr_1_n, rr_tol, alpha_exp, beta_exp, kk))));

  // Calculate the contribution from PDF component 2 to the covariance of
  // r_t and KK accretion tendency.
  comp_2_contrib =
      ((KK_accr_coef * precip_frac_2) *
       (((one / (two * crt2)) *
         (trivar_NNL_covar_eq(mu_eta_2, mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_eta_2, sigma_chi_2,
                              sigma_rr_2, sigma_rr_2_n, corr_eta_chi_2, corr_eta_rr_2_n,
                              corr_chi_rr_2_n, mu_eta_2, KK_accr_tndcy, KK_accr_coef, eta_tol,
                              rr_tol, alpha_exp, beta_exp, kk) +
          bivar_NL_mean_eq(mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_chi_2, sigma_rr_2, sigma_rr_2_n,
                           corr_chi_rr_2_n, rr_tol, (alpha_exp + one), beta_exp, kk))) +
        (((mu_rt_2 - rtm) - (mu_chi_2 / (two * crt2))) *
         bivar_NL_mean_eq(mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_chi_2, sigma_rr_2, sigma_rr_2_n,
                          corr_chi_rr_2_n, rr_tol, alpha_exp, beta_exp, kk))));

  // Calculate the covariance of r_t and KK accretion tendency.
  covar_rt_KK_accr = ((mixt_frac * comp_1_contrib) + ((one - mixt_frac) * comp_2_contrib));

  return covar_rt_KK_accr;
}
inline real_t
covar_thl_KK_accr(real_t mu_eta_1, real_t mu_eta_2, real_t mu_chi_1, real_t mu_chi_2,
                  real_t mu_rr_1, real_t mu_rr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
                  real_t sigma_eta_1, real_t sigma_eta_2, real_t sigma_chi_1, real_t sigma_chi_2,
                  real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n,
                  real_t corr_eta_chi_1, real_t corr_eta_chi_2, real_t corr_eta_rr_1_n,
                  real_t corr_eta_rr_2_n, real_t corr_chi_rr_1_n, real_t corr_chi_rr_2_n,
                  real_t thlm, real_t mu_thl_1, real_t mu_thl_2, real_t KK_accr_tndcy,
                  real_t KK_accr_coef, real_t eta_tol, real_t cthl1, real_t cthl2, real_t mixt_frac,
                  real_t precip_frac_1, real_t precip_frac_2, const parameters_KK& kk) {
  // Description:
  // References:
  // Eq. (E18) of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Eq. (A24) of Griffin, B. M. and V. E. Larson, 2016:  Parameterizing
  // microphysical effects on variances and covariances of moisture and heat
  // content using a multivariate probability density function: a study with
  // CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11, 4273--4295,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_thl_KK_accr;
  real_t alpha_exp;
  real_t beta_exp;
  real_t comp_1_contrib;
  real_t comp_2_contrib;

  // Values of the KK exponents.
  alpha_exp = KK_accr_rc_exp;
  beta_exp  = KK_accr_rr_exp;

  // Calculate the contribution from PDF component 1 to the covariance of
  // th_l and KK evaporation tendency.
  comp_1_contrib =
      ((KK_accr_coef * precip_frac_1) *
       (((one / (two * cthl1)) *
         (trivar_NNL_covar_eq(mu_eta_1, mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_eta_1, sigma_chi_1,
                              sigma_rr_1, sigma_rr_1_n, corr_eta_chi_1, corr_eta_rr_1_n,
                              corr_chi_rr_1_n, mu_eta_1, KK_accr_tndcy, KK_accr_coef, eta_tol,
                              rr_tol, alpha_exp, beta_exp, kk) -
          bivar_NL_mean_eq(mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_chi_1, sigma_rr_1, sigma_rr_1_n,
                           corr_chi_rr_1_n, rr_tol, (alpha_exp + one), beta_exp, kk))) +
        (((mu_thl_1 - thlm) + (mu_chi_1 / (two * cthl1))) *
         bivar_NL_mean_eq(mu_chi_1, mu_rr_1, mu_rr_1_n, sigma_chi_1, sigma_rr_1, sigma_rr_1_n,
                          corr_chi_rr_1_n, rr_tol, alpha_exp, beta_exp, kk))));

  // Calculate the contribution from PDF component 2 to the covariance of
  // th_l and KK evaporation tendency.
  comp_2_contrib =
      ((KK_accr_coef * precip_frac_2) *
       (((one / (two * cthl2)) *
         (trivar_NNL_covar_eq(mu_eta_2, mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_eta_2, sigma_chi_2,
                              sigma_rr_2, sigma_rr_2_n, corr_eta_chi_2, corr_eta_rr_2_n,
                              corr_chi_rr_2_n, mu_eta_2, KK_accr_tndcy, KK_accr_coef, eta_tol,
                              rr_tol, alpha_exp, beta_exp, kk) -
          bivar_NL_mean_eq(mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_chi_2, sigma_rr_2, sigma_rr_2_n,
                           corr_chi_rr_2_n, rr_tol, (alpha_exp + one), beta_exp, kk))) +
        (((mu_thl_2 - thlm) + (mu_chi_2 / (two * cthl2))) *
         bivar_NL_mean_eq(mu_chi_2, mu_rr_2, mu_rr_2_n, sigma_chi_2, sigma_rr_2, sigma_rr_2_n,
                          corr_chi_rr_2_n, rr_tol, alpha_exp, beta_exp, kk))));

  // Calculate the covariance of th_l and KK accretion tendency.
  covar_thl_KK_accr = ((mixt_frac * comp_1_contrib) + ((one - mixt_frac) * comp_2_contrib));

  return covar_thl_KK_accr;
}
inline real_t quadrivar_NNLL_covar_eq(
    real_t mu_x_i, real_t mu_chi_i, real_t mu_rr_i, real_t mu_Nr_i, real_t mu_rr_i_n,
    real_t mu_Nr_i_n, real_t sigma_x_i, real_t sigma_chi_i, real_t sigma_rr_i, real_t sigma_Nr_i,
    real_t sigma_rr_i_n, real_t sigma_Nr_i_n, real_t corr_x_chi_i, real_t corr_x_rr_i_n,
    real_t corr_x_Nr_i_n, real_t corr_chi_rr_i_n, real_t corr_chi_Nr_i_n, real_t corr_rr_Nr_i_n,
    real_t x_mean, real_t mc_tndcy_mean, real_t mc_coef, real_t x_tol, real_t alpha_exp_in,
    real_t beta_exp_in, real_t gamma_exp_in, const parameters_KK& kk) {
  // Description:
  // This function calculates the contribution by the ith PDF component to the
  // expression < y1'y2'_(i) >, where y1 = x1 ( = x, which is w or t), and
  // where y2 = x2^alpha x3^beta x4^gamma ( = chi^alpha r_r^beta N_r^gamma,
  // which also equals KK_evap_tndcy / KK_evap_coef).  The value of covariance
  // of x and KK evaporation tendency is:
  //
  // < x'KK_evap' > = KK_evap_coef
  //                  * ( mixt_frac < y1'y2'_(1) >
  //                      + ( 1 - mixt_frac ) < y1'y2'_(2) > ).
  //
  // One of four functions are called, based on whether x1 and/or x2 (x and/or
  // chi) vary.  Each one of these four functions is the result of an evaluated
  // integral based on the specific situation.
  // References:
  // Section J.1 of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Section S5 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t quadrivar_NNLL_covar_eq;
  double mu_x1;
  double mu_x2;
  double mu_x3;
  double mu_x4;
  double mu_x3_n;
  double mu_x4_n;
  double sigma_x1;
  double sigma_x2;
  double sigma_x3;
  double sigma_x4;
  double sigma_x3_n;
  double sigma_x4_n;
  double rho_x1x2;
  double rho_x1x3_n;
  double rho_x1x4_n;
  double rho_x2x3_n;
  double rho_x2x4_n;
  double rho_x3x4_n;
  double x1_mean;
  double x2_alpha_x3_beta_x4_gamma_mean;
  double alpha_exp;
  double beta_exp;
  double gamma_exp;
  double x1_tol;
  double x2_tol;
  double x3_tol;
  double x4_tol;
  double s_cc;

  // Means for the ith PDF component.
  mu_x1 = static_cast<double>(mu_x_i);
  mu_x2 = static_cast<double>(mu_chi_i);
  if ((beta_exp_in >= zero)) {
    mu_x3 = static_cast<double>(mu_rr_i);
  } else {
    mu_x3 = static_cast<double>(((mu_rr_i) > (rr_tol) ? (mu_rr_i) : (rr_tol)));
  }
  if ((gamma_exp_in >= zero)) {
    mu_x4 = static_cast<double>(mu_Nr_i);
  } else {
    mu_x4 = static_cast<double>(((mu_Nr_i) > (Nr_tol) ? (mu_Nr_i) : (Nr_tol)));
  }
  mu_x3_n = static_cast<double>(mu_rr_i_n);
  mu_x4_n = static_cast<double>(mu_Nr_i_n);

  // Standard deviations for the ith PDF component.
  sigma_x1   = static_cast<double>(sigma_x_i);
  sigma_x2   = static_cast<double>(sigma_chi_i);
  sigma_x3   = static_cast<double>(sigma_rr_i);
  sigma_x4   = static_cast<double>(sigma_Nr_i);
  sigma_x3_n = static_cast<double>(sigma_rr_i_n);
  sigma_x4_n = static_cast<double>(sigma_Nr_i_n);

  // Correlations for the ith PDF component.
  rho_x1x2   = static_cast<double>(corr_x_chi_i);
  rho_x1x3_n = static_cast<double>(corr_x_rr_i_n);
  rho_x1x4_n = static_cast<double>(corr_x_Nr_i_n);
  rho_x2x3_n = static_cast<double>(corr_chi_rr_i_n);
  rho_x2x4_n = static_cast<double>(corr_chi_Nr_i_n);
  rho_x3x4_n = static_cast<double>(corr_rr_Nr_i_n);

  // Overall means.
  x1_mean                        = static_cast<double>(x_mean);
  x2_alpha_x3_beta_x4_gamma_mean = static_cast<double>((mc_tndcy_mean / mc_coef));

  // Exponents.
  alpha_exp = static_cast<double>(alpha_exp_in);
  beta_exp  = static_cast<double>(beta_exp_in);
  gamma_exp = static_cast<double>(gamma_exp_in);

  // Tolerance values.
  // When the standard deviation of a variable is below the tolerance values,
  // it is considered to be zero, and the variable is considered to have a
  // constant value.
  x1_tol = static_cast<double>(x_tol);
  x2_tol = static_cast<double>(chi_tol);
  x3_tol = static_cast<double>(rr_tol);
  x4_tol = static_cast<double>(Nr_tol);

  // Determine the value of the parabolic cylinder function input value, s_cc.
  // The value s_cc is being fed into the parabolic cylinder function.  When
  // the value of s_cc is too large in magnitude (depending on the order of the
  // parabolic cylinder function), overflow occurs, and the output of the
  // parabolic cylinder function is +/-Inf.  This is primarily due to a large
  // ratio of mu_x2 to sigma_x2.  When the value of s_cc is very large, the
  // distribution of x2 is basically a spike near the mean, so x2 is treated as
  // a constant.
  if ((sigma_x2 > x2_tol)) {
    s_cc = (((mu_x2 / sigma_x2) + ((rho_x2x3_n * sigma_x3_n) * beta_exp)) +
            ((rho_x2x4_n * sigma_x4_n) * gamma_exp));
  } else {
    // Note:  s_cc is +inf when mu_x2 > 0 and sigma_x2 = 0, and s_cc is -inf
    //        when mu_x2 < 0 and sigma_x2 = 0.  Furthermore, s_cc is undefined
    //        when mu_x2 = 0 and sigma_x2 = 0.  However, within the context of
    //        this particular function, only the absolute value of s_cc is
    //        relevant, and furthermore the absolute value of s_cc is only
    //        relevant when sigma_x2 > 0.  Therefore, this statement only
    //        serves as divide-by-zero and compiler warning prevention.
    s_cc = std::numeric_limits<double>::max();
  }

  // Based on the values of sigma_x1, sigma_x2 (including the value of s_cc
  // compared to parab_cyl_max_input), sigma_x3, and sigma_x4, find the correct
  // form of the quadrivariate equation to use.

  if (((((sigma_x1 <= x1_tol) &&
         ((sigma_x2 <= x2_tol) ||
          (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input)))) &&
        (sigma_x3 <= x3_tol)) &&
       (sigma_x4 <= x4_tol))) {

    // The ith PDF component variances of x (w or eta), chi, r_r (in-precip),
    // and N_r (in-precip) are all 0.
    quadrivar_NNLL_covar_eq = static_cast<real_t>(quadrivar_NNLL_covar_const_all(
        mu_x1, mu_x2, mu_x3, mu_x4, x1_mean, x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp,
        gamma_exp, kk));

  } else {
    if ((((sigma_x1 <= x1_tol) &&
          ((sigma_x2 <= x2_tol) ||
           (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input)))) &&
         (sigma_x3 <= x3_tol))) {

      // The ith PDF component variances of x (w or eta), chi, and r_r
      // (in-precip) are 0.
      quadrivar_NNLL_covar_eq = static_cast<real_t>(quadrivar_NNLL_covar_cst_x1x2x3(
          mu_x1, mu_x2, mu_x3, mu_x4_n, sigma_x4_n, x1_mean, x2_alpha_x3_beta_x4_gamma_mean,
          alpha_exp, beta_exp, gamma_exp, kk));

    } else {
      if ((((sigma_x1 <= x1_tol) &&
            ((sigma_x2 <= x2_tol) ||
             (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input)))) &&
           (sigma_x4 <= x4_tol))) {

        // The ith PDF component variances of x (w or eta), chi, and N_r
        // (in-precip) are 0.
        quadrivar_NNLL_covar_eq = static_cast<real_t>(quadrivar_NNLL_covar_cst_x1x2x3(
            mu_x1, mu_x2, mu_x4, mu_x3_n, sigma_x3_n, x1_mean, x2_alpha_x3_beta_x4_gamma_mean,
            alpha_exp, gamma_exp, beta_exp, kk));

      } else {
        if ((((sigma_x1 <= x1_tol) && (sigma_x3 <= x3_tol)) && (sigma_x4 <= x4_tol))) {

          // The ith PDF component variances of x (w or eta), r_r (in-precip), and
          // N_r (in-precip) are 0.
          quadrivar_NNLL_covar_eq = static_cast<real_t>(quadrivar_NNLL_covar_cst_x1x3x4(
              mu_x1, mu_x2, mu_x3, mu_x4, sigma_x2, x1_mean, x2_alpha_x3_beta_x4_gamma_mean,
              alpha_exp, beta_exp, gamma_exp, kk));

        } else {
          if (((((sigma_x2 <= x2_tol) ||
                 (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input))) &&
                (sigma_x3 <= x3_tol)) &&
               (sigma_x4 <= x4_tol))) {

            // The ith PDF component variances of chi, r_r (in-precip), and N_r
            // (in-precip) are 0.
            quadrivar_NNLL_covar_eq = static_cast<real_t>(quadrivar_NNLL_covar_cst_x2x3x4(
                mu_x1, mu_x2, mu_x3, mu_x4, x1_mean, x2_alpha_x3_beta_x4_gamma_mean, alpha_exp,
                beta_exp, gamma_exp, kk));

          } else {
            if (((sigma_x1 <= x1_tol) &&
                 ((sigma_x2 <= x2_tol) ||
                  (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input))))) {

              // The ith PDF component variances of both x (w or eta) and chi are 0.
              quadrivar_NNLL_covar_eq = static_cast<real_t>(quadrivar_NNLL_covar_const_x1x2(
                  mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x3_n, sigma_x4_n, rho_x3x4_n, x1_mean,
                  x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp, kk));

            } else {
              if (((sigma_x1 <= x1_tol) && (sigma_x3 <= x3_tol))) {

                // The ith PDF component variances of both x (w or eta) and r_r
                // (in-precip) are 0.
                quadrivar_NNLL_covar_eq = static_cast<real_t>(quadrivar_NNLL_covar_const_x1x3(
                    mu_x1, mu_x2, mu_x3, mu_x4_n, sigma_x2, sigma_x4_n, rho_x2x4_n, x1_mean,
                    x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp, kk));

              } else {
                if (((sigma_x1 <= x1_tol) && (sigma_x4 <= x4_tol))) {

                  // The ith PDF component variances of both x (w or eta) and N_r
                  // (in-precip) are 0.
                  quadrivar_NNLL_covar_eq = static_cast<real_t>(quadrivar_NNLL_covar_const_x1x3(
                      mu_x1, mu_x2, mu_x4, mu_x3_n, sigma_x2, sigma_x3_n, rho_x2x3_n, x1_mean,
                      x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, gamma_exp, beta_exp, kk));

                } else {
                  if ((((sigma_x2 <= x2_tol) ||
                        (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input))) &&
                       (sigma_x3 <= x3_tol))) {

                    // The ith PDF component variances of both chi and r_r (in-precip) are 0.
                    quadrivar_NNLL_covar_eq = static_cast<real_t>(quadrivar_NNLL_covar_const_x2x3(
                        mu_x1, mu_x2, mu_x3, mu_x4_n, sigma_x1, sigma_x4_n, rho_x1x4_n, x1_mean,
                        x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp, kk));

                  } else {
                    if ((((sigma_x2 <= x2_tol) ||
                          (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input))) &&
                         (sigma_x4 <= x4_tol))) {

                      // The ith PDF component variances of both chi and N_r (in-precip) are 0.
                      quadrivar_NNLL_covar_eq = static_cast<real_t>(quadrivar_NNLL_covar_const_x2x3(
                          mu_x1, mu_x2, mu_x4, mu_x3_n, sigma_x1, sigma_x3_n, rho_x1x3_n, x1_mean,
                          x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, gamma_exp, beta_exp, kk));

                    } else {
                      if (((sigma_x3 <= x3_tol) && (sigma_x4 <= x4_tol))) {

                        // The ith PDF component variances of both r_r (in-precip) and N_r
                        // (in-precip) are 0.
                        quadrivar_NNLL_covar_eq =
                            static_cast<real_t>(quadrivar_NNLL_covar_const_x3x4(
                                mu_x1, mu_x2, mu_x3, mu_x4, sigma_x1, sigma_x2, rho_x1x2, x1_mean,
                                x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp,
                                kk));

                      } else {
                        if ((sigma_x1 <= x1_tol)) {

                          // The ith PDF component variance of x (w or eta) is 0.
                          quadrivar_NNLL_covar_eq =
                              static_cast<real_t>(quadrivar_NNLL_covar_const_x1(
                                  mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x2, sigma_x3_n, sigma_x4_n,
                                  rho_x2x3_n, rho_x2x4_n, rho_x3x4_n, x1_mean,
                                  x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp,
                                  kk));

                        } else {
                          if (((sigma_x2 <= x2_tol) ||
                               (Kokkos::abs(s_cc) > static_cast<double>(parab_cyl_max_input)))) {

                            // The ith PDF component variance of chi is 0.
                            quadrivar_NNLL_covar_eq =
                                static_cast<real_t>(quadrivar_NNLL_covar_const_x2(
                                    mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x1, sigma_x3_n,
                                    sigma_x4_n, rho_x1x3_n, rho_x1x4_n, rho_x3x4_n, x1_mean,
                                    x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp, gamma_exp,
                                    kk));

                          } else {
                            if ((sigma_x3 <= x3_tol)) {

                              // The ith PDF component variance of r_r (in-precip) is 0.
                              quadrivar_NNLL_covar_eq =
                                  static_cast<real_t>(quadrivar_NNLL_covar_const_x3(
                                      mu_x1, mu_x2, mu_x3, mu_x4_n, sigma_x1, sigma_x2, sigma_x4_n,
                                      rho_x1x2, rho_x1x4_n, rho_x2x4_n, x1_mean,
                                      x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, beta_exp,
                                      gamma_exp, kk));

                            } else {
                              if ((sigma_x4 <= x4_tol)) {

                                // The ith PDF component variance of N_r (in-precip) is 0.
                                quadrivar_NNLL_covar_eq =
                                    static_cast<real_t>(quadrivar_NNLL_covar_const_x3(
                                        mu_x1, mu_x2, mu_x4, mu_x3_n, sigma_x1, sigma_x2,
                                        sigma_x3_n, rho_x1x2, rho_x1x3_n, rho_x2x3_n, x1_mean,
                                        x2_alpha_x3_beta_x4_gamma_mean, alpha_exp, gamma_exp,
                                        beta_exp, kk));

                              } else {

                                // This is the complete value of the quadrivariate.
                                // All fields vary in the ith PDF component.
                                quadrivar_NNLL_covar_eq = static_cast<real_t>(quadrivar_NNLL_covar(
                                    mu_x1, mu_x2, mu_x3_n, mu_x4_n, sigma_x1, sigma_x2, sigma_x3_n,
                                    sigma_x4_n, rho_x1x2, rho_x1x3_n, rho_x1x4_n, rho_x2x3_n,
                                    rho_x2x4_n, rho_x3x4_n, x1_mean, x2_alpha_x3_beta_x4_gamma_mean,
                                    alpha_exp, beta_exp, gamma_exp, kk));
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return quadrivar_NNLL_covar_eq;
}
inline real_t trivar_NNL_covar_eq(real_t mu_x_i, real_t mu_chi_i, real_t mu_y_i, real_t mu_y_i_n,
                                  real_t sigma_x_i, real_t sigma_chi_i, real_t sigma_y_i,
                                  real_t sigma_y_i_n, real_t corr_x_chi_i, real_t corr_x_y_i_n,
                                  real_t corr_chi_y_i_n, real_t x_mean, real_t mc_tndcy_mean,
                                  real_t mc_coef, real_t x_tol, real_t y_tol, real_t alpha_exp_in,
                                  real_t beta_exp_in, const parameters_KK& kk) {
  // Description:
  // This function calculates the contribution by the ith PDF component to the
  // expression < y1'y2'_(i) >, where y1 = x1 ( = x, which is w or t), and
  // where y2 = x2^alpha x3^beta ( = chi^alpha y^beta, where y is N_cn or r_r
  // for autoconversion or accretion, respectively, and which also equals
  // KK_auto_tndcy / KK_auto_coef or KK_accr_tndcy / KK_accr_coef,
  // respectively).  The value of covariance of x and the KK microphysics
  // tendency is:
  //
  // < x'KK_mc' > = KK_mc_coef
  //                * ( mixt_frac < y1'y2'_(1) >
  //                    + ( 1 - mixt_frac ) < y1'y2'_(2) > ).
  //
  // One of four functions are called, based on whether x1 and/or x2 (x and/or
  // chi) vary.  Each one of these four functions is the result of an evaluated
  // integral based on the specific situation.
  // References:
  // Section J.2 of Griffin, B. M., 2016:  Improving the Subgrid-Scale
  // Representation of Hydrometeors and Microphysical Feedback Effects Using a
  // Multivariate PDF.  Doctoral dissertation, University of
  // Wisconsin -- Milwaukee, Milwaukee, WI, Paper 1144, 165 pp., URL
  // http://dc.uwm.edu/cgi/viewcontent.cgi?article=2149&context=etd.
  //
  // Section S6 of Griffin, B. M. and V. E. Larson, 2016:  Supplement of
  // Parameterizing microphysical effects on variances and covariances of
  // moisture and heat content using a multivariate probability density
  // function: a study with CLUBB (tag MVCS).  Geosci. Model Dev., 9, 11,
  // doi:http://dx.doi.org/10.5194/gmd-9-4273-2016-supplement.
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t trivar_NNL_covar_eq;
  double mu_x1;
  double mu_x2;
  double mu_x3;
  double mu_x3_n;
  double sigma_x1;
  double sigma_x2;
  double sigma_x3;
  double sigma_x3_n;
  double rho_x1x2;
  double rho_x1x3_n;
  double rho_x2x3_n;
  double x1_mean;
  double x2_alpha_x3_beta_mean;
  double alpha_exp;
  double beta_exp;
  double x1_tol;
  double x2_tol;
  double x3_tol;
  double s_c;

  // Means for the ith PDF component.
  mu_x1 = static_cast<double>(mu_x_i);
  mu_x2 = static_cast<double>(mu_chi_i);
  if ((beta_exp_in >= zero)) {
    mu_x3 = static_cast<double>(mu_y_i);
  } else {
    mu_x3 = static_cast<double>(((mu_y_i) > (y_tol) ? (mu_y_i) : (y_tol)));
  }
  mu_x3_n = static_cast<double>(mu_y_i_n);

  // Standard deviations for the ith PDF component.
  sigma_x1   = static_cast<double>(sigma_x_i);
  sigma_x2   = static_cast<double>(sigma_chi_i);
  sigma_x3   = static_cast<double>(sigma_y_i);
  sigma_x3_n = static_cast<double>(sigma_y_i_n);

  // Correlations for the ith PDF component.
  rho_x1x2   = static_cast<double>(corr_x_chi_i);
  rho_x1x3_n = static_cast<double>(corr_x_y_i_n);
  rho_x2x3_n = static_cast<double>(corr_chi_y_i_n);

  // Overall means.
  x1_mean               = static_cast<double>(x_mean);
  x2_alpha_x3_beta_mean = static_cast<double>((mc_tndcy_mean / mc_coef));

  // Exponents.
  alpha_exp = static_cast<double>(alpha_exp_in);
  beta_exp  = static_cast<double>(beta_exp_in);

  // Tolerance values.
  // When the standard deviation of a variable is below the tolerance values,
  // it is considered to be zero, and the variable is considered to have a
  // constant value.
  x1_tol = static_cast<double>(x_tol);
  x2_tol = static_cast<double>(chi_tol);
  x3_tol = static_cast<double>(y_tol);

  // Determine the value of the parabolic cylinder function input value, s_c.
  // The value s_c is being fed into the parabolic cylinder function.  When
  // the value of s_c is too large in magnitude (depending on the order of the
  // parabolic cylinder function), overflow occurs, and the output of the
  // parabolic cylinder function is +/-Inf.  This is primarily due to a large
  // ratio of mu_x2 to sigma_x2.  When the value of s_c is very large, the
  // distribution of x2 is basically a spike near the mean, so x2 is treated as
  // a constant.
  if ((sigma_x2 > x2_tol)) {
    s_c = ((mu_x2 / sigma_x2) + ((rho_x2x3_n * sigma_x3_n) * beta_exp));
  } else {
    // Note:  s_c is +inf when mu_x2 > 0 and sigma_x2 = 0, and s_c is -inf
    //        when mu_x2 < 0 and sigma_x2 = 0.  Furthermore, s_c is undefined
    //        when mu_x2 = 0 and sigma_x2 = 0.  However, within the context of
    //        this particular function, only the absolute value of s_c is
    //        relevant, and furthermore the absolute value of s_c is only
    //        relevant when sigma_x2 > 0.  Therefore, this statement only
    //        serves as divide-by-zero and compiler warning prevention.
    s_c = std::numeric_limits<double>::max();
  }

  // Based on the values of sigma_x1, sigma_x2 (including the value of s_c
  // compared to parab_cyl_max_input), and sigma_x3, find the correct form of
  // the trivariate equation to use.

  if ((((sigma_x1 <= x1_tol) &&
        ((sigma_x2 <= x2_tol) || (Kokkos::abs(s_c) > static_cast<double>(parab_cyl_max_input)))) &&
       (sigma_x3 <= x3_tol))) {

    // The ith PDF component variances of x (w or eta), chi, and y (N_cn or
    // r_r (in-precip)) are all 0.
    trivar_NNL_covar_eq = static_cast<real_t>(trivar_NNL_covar_const_all(
        mu_x1, mu_x2, mu_x3, x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp, kk));

  } else {
    if (((sigma_x1 <= x1_tol) &&
         ((sigma_x2 <= x2_tol) || (Kokkos::abs(s_c) > static_cast<double>(parab_cyl_max_input))))) {

      // The ith PDF component variances of both x (w or eta) and chi are 0.
      trivar_NNL_covar_eq = static_cast<real_t>(
          trivar_NNL_covar_const_x1x2(mu_x1, mu_x2, mu_x3_n, sigma_x3_n, x1_mean,
                                      x2_alpha_x3_beta_mean, alpha_exp, beta_exp, kk));

    } else {
      if (((sigma_x1 <= x1_tol) && (sigma_x3 <= x3_tol))) {

        // The ith PDF component variances of both x (w or eta) and y (N_cn or
        // r_r (in-precip)) are 0.
        trivar_NNL_covar_eq = static_cast<real_t>(
            trivar_NNL_covar_const_x1x3(mu_x1, mu_x2, mu_x3, sigma_x2, x1_mean,
                                        x2_alpha_x3_beta_mean, alpha_exp, beta_exp, kk));

      } else {
        if ((((sigma_x2 <= x2_tol) ||
              (Kokkos::abs(s_c) > static_cast<double>(parab_cyl_max_input))) &&
             (sigma_x3 <= x3_tol))) {

          // The ith PDF component variances of both chi and y (N_cn or
          // r_r (in-precip)) are 0.
          trivar_NNL_covar_eq = static_cast<real_t>(trivar_NNL_covar_const_x2x3(
              mu_x1, mu_x2, mu_x3, x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp, kk));

        } else {
          if ((sigma_x1 <= x1_tol)) {

            // The ith PDF component variance of x (w or eta) is 0.
            trivar_NNL_covar_eq = static_cast<real_t>(
                trivar_NNL_covar_const_x1(mu_x1, mu_x2, mu_x3_n, sigma_x2, sigma_x3_n, rho_x2x3_n,
                                          x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp, kk));

          } else {
            if (((sigma_x2 <= x2_tol) ||
                 (Kokkos::abs(s_c) > static_cast<double>(parab_cyl_max_input)))) {

              // The ith PDF component variance of chi is 0.
              trivar_NNL_covar_eq = static_cast<real_t>(trivar_NNL_covar_const_x2(
                  mu_x1, mu_x2, mu_x3_n, sigma_x1, sigma_x3_n, rho_x1x3_n, x1_mean,
                  x2_alpha_x3_beta_mean, alpha_exp, beta_exp, kk));

            } else {
              if ((sigma_x3 <= x3_tol)) {

                // The ith PDF component variance of y (N_cn or r_r (in-precip)) is 0.
                trivar_NNL_covar_eq = static_cast<real_t>(trivar_NNL_covar_const_x3(
                    mu_x1, mu_x2, mu_x3, sigma_x1, sigma_x2, rho_x1x2, x1_mean,
                    x2_alpha_x3_beta_mean, alpha_exp, beta_exp, kk));

              } else {

                // This is the complete value of the trivariate.
                // All fields vary in the ith PDF component.
                trivar_NNL_covar_eq = static_cast<real_t>(trivar_NNL_covar(
                    mu_x1, mu_x2, mu_x3_n, sigma_x1, sigma_x2, sigma_x3_n, rho_x1x2, rho_x1x3_n,
                    rho_x2x3_n, x1_mean, x2_alpha_x3_beta_mean, alpha_exp, beta_exp, kk));
              }
            }
          }
        }
      }
    }
  }

  return trivar_NNL_covar_eq;
}
} // namespace clubb_kokkos::kk_engine
