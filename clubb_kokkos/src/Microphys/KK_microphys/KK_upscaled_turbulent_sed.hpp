// Native KK_upscaled_turbulent_sed numerical formulas, evaluated in source arithmetic order.
// Appended caller-owned KK parameters replace implicit module state.
#pragma once
#include "Microphys/KK_microphys/engine_contracts.hpp"
#include "Microphys/KK_microphys/KK_utilities.hpp"
namespace clubb_kokkos::kk_engine {
inline void KK_sed_vel_covars(real_t rrm, real_t rr_1, real_t rr_2, real_t Nrm, real_t Nr_1,
                              real_t Nr_2, real_t KK_mean_vol_rad, real_t mu_rr_1, real_t mu_rr_2,
                              real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n, real_t mu_rr_2_n,
                              real_t mu_Nr_1_n, real_t mu_Nr_2_n, real_t sigma_rr_1,
                              real_t sigma_rr_2, real_t sigma_Nr_1, real_t sigma_Nr_2,
                              real_t sigma_rr_1_n, real_t sigma_rr_2_n, real_t sigma_Nr_1_n,
                              real_t sigma_Nr_2_n, real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
                              real_t KK_mvr_coef, real_t mixt_frac, real_t precip_frac_1,
                              real_t precip_frac_2, real_t& Vrrprrp_impc, real_t& Vrrprrp_expc,
                              real_t& VNrpNrp_impc, real_t& VNrpNrp_expc, real_t& rr_KK_mvr_covar,
                              real_t& Nr_KK_mvr_covar, const parameters_KK& kk) {
  // Description:
  // The covariance of sedimentation velocity (of rain water mixing ratio) and
  // rain water mixing ratio, < V_rr'r_r' > is given by the equation:
  //
  // < V_rr'r_r' > = 0.012 ( 10^6 microm/m ) < r_r'R_vr' >.
  //
  // Likewise, the covariance of sedimentation velocity (of rain drop
  // concentration) and rain drop concentration, < V_Nr'N_r' > is given by the
  // equation:
  //
  // < V_Nr'N_r' > = 0.007 ( 10^6 microm/m ) < N_r'R_vr' >.
  //
  // The covariances < r_r'R_vr' > and < N_r'R_vr' > are calculated in this
  // subroutine.
  //
  // It is advantageous to write the < V_rr'r_r' > term and the < V_Nr'N_r' >
  // term in terms of < r_r > and < N_r >, respectively, because < V_rr'r_r' >
  // and < V_Nr'N_r' > are found in the turbulent sedimentation term as part
  // of the predictive equations for < r_r > and < N_r >, respectively.  This
  // allows the turbulent sedimentation term to be handled semi-implicitly in
  // the < r_r > and < N_r > predictive equations, increasing the numerical
  // stability of the solutions.  However, the option to handle these terms
  // in the traditional, completely explicit and less computationally-expensive
  // manner is included.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Output Variables
  // Local Variables
  real_t rr_KK_mvr_covar_coefA;
  real_t rr_KK_mvr_covar_termB;
  real_t Nr_KK_mvr_covar_coefA;
  real_t Nr_KK_mvr_covar_termB;
  constexpr bool l_semi_imp_turbulent_sed = true;

  // Calculate the covariance of rain drop mean volume radius and r_r,
  // < R_vr'r_r' >.
  if (l_semi_imp_turbulent_sed) {

    // Turbulent sedimentation will be handled semi-implicitly in the
    // hydrometeor predictive equation set.

    rr_KK_mvr_covar_coefA = covar_rr_KK_mvr_coefA(
        mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1,
        sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
        corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mean_vol_rad, KK_mvr_coef, kk);

    rr_KK_mvr_covar_termB = covar_rr_KK_mvr_termB(
        rr_1, rr_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
        sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n,
        sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mvr_coef, mixt_frac, kk);

  } else {

    // Turbulent sedimentation will be handled completely explicitly in the
    // hydrometeor predictive equation set.

    rr_KK_mvr_covar_coefA = zero;

    rr_KK_mvr_covar_termB = covar_rr_KK_mvr(
        mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1,
        sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
        corr_rr_Nr_1_n, corr_rr_Nr_2_n, rrm, KK_mean_vol_rad, KK_mvr_coef, mixt_frac, precip_frac_1,
        precip_frac_2, kk);
  }

  // Calculate the covariance of rain drop mean volume radius and N_r,
  // < R_vr'N_r' >.
  if (l_semi_imp_turbulent_sed) {

    // Turbulent sedimentation will be handled semi-implicitly in the
    // hydrometeor predictive equation set.

    Nr_KK_mvr_covar_coefA = covar_Nr_KK_mvr_coefA(
        mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1,
        sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
        corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mean_vol_rad, KK_mvr_coef, kk);

    Nr_KK_mvr_covar_termB = covar_Nr_KK_mvr_termB(
        Nr_1, Nr_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n,
        sigma_rr_1, sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n,
        sigma_Nr_2_n, corr_rr_Nr_1_n, corr_rr_Nr_2_n, KK_mvr_coef, mixt_frac, kk);

  } else {

    // Turbulent sedimentation will be handled completely explicitly in the
    // hydrometeor predictive equation set.

    Nr_KK_mvr_covar_coefA = zero;

    Nr_KK_mvr_covar_termB = covar_Nr_KK_mvr(
        mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, mu_Nr_2_n, sigma_rr_1,
        sigma_rr_2, sigma_Nr_1, sigma_Nr_2, sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, sigma_Nr_2_n,
        corr_rr_Nr_1_n, corr_rr_Nr_2_n, Nrm, KK_mean_vol_rad, KK_mvr_coef, mixt_frac, precip_frac_1,
        precip_frac_2, kk);
  }

  // Terms used to calculate the covariance of V_rr and r_r, < V_rr'r_r' >.
  Vrrprrp_impc = (((-real_t(0.012)) * micron_per_m) * rr_KK_mvr_covar_coefA);
  Vrrprrp_expc = (((-real_t(0.012)) * micron_per_m) * rr_KK_mvr_covar_termB);

  // Terms used to calculate the covariance of V_Nr and N_r, < V_Nr'N_r' >.
  VNrpNrp_impc = (((-real_t(0.007)) * micron_per_m) * Nr_KK_mvr_covar_coefA);
  VNrpNrp_expc = (((-real_t(0.007)) * micron_per_m) * Nr_KK_mvr_covar_termB);

  rr_KK_mvr_covar = ((rr_KK_mvr_covar_coefA * rrm) + rr_KK_mvr_covar_termB);
  Nr_KK_mvr_covar = ((Nr_KK_mvr_covar_coefA * Nrm) + Nr_KK_mvr_covar_termB);

  return;
}
inline real_t covar_rr_KK_mvr_coefA(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                                    real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                                    real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                                    real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                                    real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                                    real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
                                    real_t KK_mean_vol_rad, real_t KK_mvr_coef,
                                    const parameters_KK& kk) {
  // Description:
  // This function partially calculates the covariance of r_r and KK mean
  // volume radius of rain drops (R_vr), which can be written as < r_r'R_vr' >.
  //
  // The equation for the covariance of r_r and KK rain drop mean volume radius
  // is:
  //
  // < r_r'R_vr' >
  // = KK_mvr_coef
  //   ( a f_p(1)
  //     exp{ mu_rr_1_n ( alpha + 1 ) + mu_Nr_1_n beta
  //          + (1/2) sigma_rr_1_n^2 ( alpha + 1 )^2
  //          + (1/2) sigma_Nr_1_n^2 beta^2
  //          + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
  //   + ( 1 - a ) f_p(2)
  //     exp{ mu_rr_2_n ( alpha + 1 ) + mu_Nr_2_n beta
  //          + (1/2) sigma_rr_2_n^2 ( alpha + 1 )^2
  //          + (1/2) sigma_Nr_2_n^2 beta^2
  //          + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
  //   )
  //   - < R_vr > < r_r >;
  //
  // where "a" is the mixture fraction (weight of the 1st PDF component), and
  // f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
  // components, respectively.  The equation for < r_r'R_vr' > needs to be put
  // in the form:  < r_r'R_vr' > = coef_A < r_r > + term_B.  Since < r_r > can
  // be written as:
  //
  // < r_r > = a f_p(1) exp{ mu_rr_1_n + (1/2) sigma_rr_1_n^2 }
  //           + ( 1 - a ) f_p(2) exp{ mu_rr_2_n + (1/2) sigma_rr_2_n^2 }
  //         = a f_p(1) mu_rr_1 + ( 1 - a ) f_p(2) mu_rr_2
  //         = a rr_1 + ( 1 - a ) rr_2;
  //
  // the covariance < r_r'R_vr' > can be rewritten in terms of < r_r >:
  //
  // < r_r'R_vr' >
  // = ( KK_mvr_coef
  //     ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_1_n^2 beta^2
  //            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
  //     + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_2_n^2 beta^2
  //            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
  //     )
  //     - < R_vr > ) < r_r >
  //   - KK_mvr_coef
  //     ( a rr_1
  //       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_2_n^2 beta^2
  //            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
  //     + ( 1 - a ) rr_2
  //       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_1_n^2 beta^2
  //            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
  //     ).
  //
  // This function calculates the coefficient A, which is given by:
  //
  // coef_A
  // = KK_mvr_coef
  //   ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //          + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
  //          + (1/2) sigma_Nr_1_n^2 beta^2
  //          + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
  //   + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //          + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
  //          + (1/2) sigma_Nr_2_n^2 beta^2
  //          + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
  //   )
  //   - < R_vr >.
  //
  // The covariance of sedimentation velocity (of rain water mixing ratio) and
  // rain water mixing ratio, < V_rr'r_r' > is given by the equation:
  //
  // < V_rr'r_r' > = 0.012 ( 10^6 microm/m ) < r_r'R_vr' >.
  //
  // It is advantageous to write the < V_rr'r_r' > term in terms of < r_r >
  // because < V_rr'r_r' > is found in the turbulent sedimentation term as part
  // of the predictive equation for < r_r >.  This allows the turbulent
  // sedimentation term to be handled semi-implicitly in the < r_r > predictive
  // equation, increasing the numerical stability of the solution.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_rr_KK_mvr_coefA;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_mvr_rr_exp;
  beta_exp  = KK_mvr_Nr_exp;

  // Calculate the coefficient A in the < r_r'R_vr' > covariance equation,
  // where < r_r'R_vr' > = coef_A <r_r> + term_B, and:
  //
  // coef_A
  // = KK_mvr_coef
  //   ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //          + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
  //          + (1/2) sigma_Nr_1_n^2 beta^2
  //          + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
  //   + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //          + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
  //          + (1/2) sigma_Nr_2_n^2 beta^2
  //          + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
  //   )
  //   - < R_vr >.
  covar_rr_KK_mvr_coefA =
      ((KK_mvr_coef * (bivar_LL_covar_partial_rr(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_rr_1,
                                                 sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n,
                                                 corr_rr_Nr_1_n, alpha_exp, beta_exp, kk) +
                       bivar_LL_covar_partial_rr(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_rr_2,
                                                 sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n,
                                                 corr_rr_Nr_2_n, alpha_exp, beta_exp, kk))) -
       KK_mean_vol_rad);

  return covar_rr_KK_mvr_coefA;
}
inline real_t covar_rr_KK_mvr_termB(real_t rr_1, real_t rr_2, real_t mu_rr_1, real_t mu_rr_2,
                                    real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n,
                                    real_t mu_rr_2_n, real_t mu_Nr_1_n, real_t mu_Nr_2_n,
                                    real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1,
                                    real_t sigma_Nr_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n,
                                    real_t sigma_Nr_1_n, real_t sigma_Nr_2_n, real_t corr_rr_Nr_1_n,
                                    real_t corr_rr_Nr_2_n, real_t KK_mvr_coef, real_t mixt_frac,
                                    const parameters_KK& kk) {
  // Description:
  // This function partially calculates the covariance of r_r and KK mean
  // volume radius of rain drops (R_vr), which can be written as < r_r'R_vr' >.
  //
  // The equation for the covariance of r_r and KK rain drop mean volume radius
  // is:
  //
  // < r_r'R_vr' >
  // = KK_mvr_coef
  //   ( a f_p(1)
  //     exp{ mu_rr_1_n ( alpha + 1 ) + mu_Nr_1_n beta
  //          + (1/2) sigma_rr_1_n^2 ( alpha + 1 )^2
  //          + (1/2) sigma_Nr_1_n^2 beta^2
  //          + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
  //   + ( 1 - a ) f_p(2)
  //     exp{ mu_rr_2_n ( alpha + 1 ) + mu_Nr_2_n beta
  //          + (1/2) sigma_rr_2_n^2 ( alpha + 1 )^2
  //          + (1/2) sigma_Nr_2_n^2 beta^2
  //          + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
  //   )
  //   - < R_vr > < r_r >;
  //
  // where "a" is the mixture fraction (weight of the 1st PDF component), and
  // f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
  // components, respectively.  The equation for < r_r'R_vr' > needs to be put
  // in the form:  < r_r'R_vr' > = coef_A < r_r > + term_B.  Since < r_r > can
  // be written as:
  //
  // < r_r > = a f_p(1) exp{ mu_rr_1_n + (1/2) sigma_rr_1_n^2 }
  //           + ( 1 - a ) f_p(2) exp{ mu_rr_2_n + (1/2) sigma_rr_2_n^2 }
  //         = a f_p(1) mu_rr_1 + ( 1 - a ) f_p(2) mu_rr_2
  //         = a rr_1 + ( 1 - a ) rr_2;
  //
  // the covariance < r_r'R_vr' > can be rewritten in terms of < r_r >:
  //
  // < r_r'R_vr' >
  // = ( KK_mvr_coef
  //     ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_1_n^2 beta^2
  //            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
  //     + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_2_n^2 beta^2
  //            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
  //     )
  //     - < R_vr > ) < r_r >
  //   - KK_mvr_coef
  //     ( a rr_1
  //       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_2_n^2 beta^2
  //            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
  //     + ( 1 - a ) rr_2
  //       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_1_n^2 beta^2
  //            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
  //     ).
  //
  // This function calculates term B, which is given by:
  //
  // term_B
  // = - KK_mvr_coef
  //     ( a rr_1
  //       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_2_n^2 beta^2
  //            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
  //     + ( 1 - a ) rr_2
  //       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_1_n^2 beta^2
  //            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
  //     ).
  //
  // The covariance of sedimentation velocity (of rain water mixing ratio) and
  // rain water mixing ratio, < V_rr'r_r' > is given by the equation:
  //
  // < V_rr'r_r' > = 0.012 ( 10^6 microm/m ) < r_r'R_vr' >.
  //
  // It is advantageous to write the < V_rr'r_r' > term in terms of < r_r >
  // because < V_rr'r_r' > is found in the turbulent sedimentation term as part
  // of the predictive equation for < r_r >.  This allows the turbulent
  // sedimentation term to be handled semi-implicitly in the < r_r > predictive
  // equation, increasing the numerical stability of the solution.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_rr_KK_mvr_termB;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_mvr_rr_exp;
  beta_exp  = KK_mvr_Nr_exp;

  // Calculate term B in the < r_r'R_vr' > covariance equation, where
  // < r_r'R_vr' > = coef_A <r_r> + term_B, and:
  //
  // term_B
  // = - KK_mvr_coef
  //     ( a rr_1
  //       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_2_n^2 beta^2
  //            + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
  //     + ( 1 - a ) rr_2
  //       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 ( alpha^2 + 2 alpha )
  //            + (1/2) sigma_Nr_1_n^2 beta^2
  //            + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
  //     ).
  covar_rr_KK_mvr_termB =
      ((-KK_mvr_coef) *
       (((mixt_frac * rr_1) * bivar_LL_covar_partial_rr(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n,
                                                        sigma_rr_2, sigma_Nr_2, sigma_rr_2_n,
                                                        sigma_Nr_2_n, corr_rr_Nr_2_n, alpha_exp,
                                                        beta_exp, kk)) +
        (((one - mixt_frac) * rr_2) *
         bivar_LL_covar_partial_rr(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_rr_1, sigma_Nr_1,
                                   sigma_rr_1_n, sigma_Nr_1_n, corr_rr_Nr_1_n, alpha_exp, beta_exp,
                                   kk))));

  return covar_rr_KK_mvr_termB;
}
inline real_t covar_rr_KK_mvr(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                              real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                              real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                              real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                              real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                              real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t rrm,
                              real_t KK_mean_vol_rad, real_t KK_mvr_coef, real_t mixt_frac,
                              real_t precip_frac_1, real_t precip_frac_2, const parameters_KK& kk) {
  // Description:
  // This function calculates the covariance of r_r and KK mean volume radius
  // of rain drops (R_vr), which can be written as < r_r'R_vr' >.
  //
  // The equation for the covariance of r_r and KK rain drop mean volume radius
  // is:
  //
  // < r_r'R_vr' >
  // = KK_mvr_coef
  //   ( a f_p(1)
  //     exp{ mu_rr_1_n ( alpha + 1 ) + mu_Nr_1_n beta
  //          + (1/2) sigma_rr_1_n^2 ( alpha + 1 )^2
  //          + (1/2) sigma_Nr_1_n^2 beta^2
  //          + corr_rr_Nr_1_n sigma_rr_1_n ( alpha + 1 ) sigma_Nr_1_n beta }
  //   + ( 1 - a ) f_p(2)
  //     exp{ mu_rr_2_n ( alpha + 1 ) + mu_Nr_2_n beta
  //          + (1/2) sigma_rr_2_n^2 ( alpha + 1 )^2
  //          + (1/2) sigma_Nr_2_n^2 beta^2
  //          + corr_rr_Nr_2_n sigma_rr_2_n ( alpha + 1 ) sigma_Nr_2_n beta }
  //   )
  //   - < R_vr > < r_r >;
  //
  // where "a" is the mixture fraction (weight of the 1st PDF component), and
  // f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
  // components, respectively.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_rr_KK_mvr;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_mvr_rr_exp;
  beta_exp  = KK_mvr_Nr_exp;

  // Calculate the covariance of r_r and KK mean volume radius of rain drops.
  covar_rr_KK_mvr =
      ((KK_mvr_coef * (((mixt_frac * precip_frac_1) *
                        bivar_LL_mean_eq(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_rr_1,
                                         sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_rr_Nr_1_n,
                                         (alpha_exp + one), beta_exp, kk)) +
                       (((one - mixt_frac) * precip_frac_2) *
                        bivar_LL_mean_eq(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_rr_2,
                                         sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_rr_Nr_2_n,
                                         (alpha_exp + one), beta_exp, kk)))) -
       (rrm * KK_mean_vol_rad));

  return covar_rr_KK_mvr;
}
inline real_t covar_Nr_KK_mvr_coefA(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                                    real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                                    real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                                    real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                                    real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                                    real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n,
                                    real_t KK_mean_vol_rad, real_t KK_mvr_coef,
                                    const parameters_KK& kk) {
  // Description:
  // This function partially calculates the covariance of N_r and KK mean
  // volume radius of rain drops (R_vr), which can be written as < N_r'R_vr' >.
  //
  // The equation for the covariance of N_r and KK rain drop mean volume radius
  // is:
  //
  // < N_r'R_vr' >
  // = KK_mvr_coef
  //   ( a f_p(1)
  //     exp{ mu_rr_1_n alpha + mu_Nr_1_n ( beta + 1 )
  //          + (1/2) sigma_rr_1_n^2 alpha^2
  //          + (1/2) sigma_Nr_1_n^2 ( beta + 1 )^2
  //          + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
  //   + ( 1 - a ) f_p(2)
  //     exp{ mu_rr_2_n alpha + mu_Nr_2_n ( beta + 1 )
  //          + (1/2) sigma_rr_2_n^2 alpha^2
  //          + (1/2) sigma_Nr_2_n^2 ( beta + 1 )^2
  //          + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
  //   )
  //   - < R_vr > < N_r >;
  //
  // where "a" is the mixture fraction (weight of the 1st PDF component), and
  // f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
  // components, respectively.  The equation for < N_r'R_vr' > needs to be put
  // in the form:  < N_r'R_vr' > = coef_A < N_r > + term_B.  Since < N_r > can
  // be written as:
  //
  // < N_r > = a f_p(1) exp{ mu_Nr_1_n + (1/2) sigma_Nr_1_n^2 }
  //           + ( 1 - a ) f_p(2) exp{ mu_Nr_2_n + (1/2) sigma_Nr_2_n^2 }
  //         = a f_p(1) mu_Nr_1 + ( 1 - a ) f_p(2) mu_Nr_2
  //         = a Nr_1 + ( 1 - a ) Nr_2;
  //
  // the covariance < N_r'R_vr' > can be rewritten in terms of < N_r >:
  //
  // < N_r'R_vr' >
  // = ( KK_mvr_coef
  //     ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 alpha^2
  //            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
  //     + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 alpha^2
  //            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
  //     )
  //     - < R_vr > ) < N_r >
  //   - KK_mvr_coef
  //     ( a Nr_1
  //       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 alpha^2
  //            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
  //     + ( 1 - a ) Nr_2
  //       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 alpha^2
  //            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
  //     ).
  //
  // This function calculates the coefficient A, which is given by:
  //
  // coef_A
  // = KK_mvr_coef
  //   ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //          + (1/2) sigma_rr_1_n^2 alpha^2
  //          + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
  //          + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
  //   + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //          + (1/2) sigma_rr_2_n^2 alpha^2
  //          + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
  //          + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
  //   )
  //   - < R_vr >.
  //
  // The covariance of sedimentation velocity (of rain drop concentration) and
  // rain drop concentration, < V_Nr'N_r' > is given by the equation:
  //
  // < V_Nr'N_r' > = 0.007 ( 10^6 microm/m ) < N_r'R_vr' >.
  //
  // It is advantageous to write the < V_Nr'N_r' > term in terms of < N_r >
  // because < V_Nr'N_r' > is found in the turbulent sedimentation term as part
  // of the predictive equation for < N_r >.  This allows the turbulent
  // sedimentation term to be handled semi-implicitly in the < N_r > predictive
  // equation, increasing the numerical stability of the solution.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_Nr_KK_mvr_coefA;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_mvr_rr_exp;
  beta_exp  = KK_mvr_Nr_exp;

  // Calculate the coefficient A in the < N_r'R_vr' > covariance equation,
  // where < N_r'R_vr' > = coef_A <N_r> + term_B, and:
  //
  // coef_A
  // = KK_mvr_coef
  //   ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //          + (1/2) sigma_rr_1_n^2 alpha^2
  //          + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
  //          + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
  //   + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //          + (1/2) sigma_rr_2_n^2 alpha^2
  //          + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
  //          + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
  //   )
  //   - < R_vr >.
  covar_Nr_KK_mvr_coefA =
      ((KK_mvr_coef * (bivar_LL_covar_partial_Nr(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_rr_1,
                                                 sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n,
                                                 corr_rr_Nr_1_n, alpha_exp, beta_exp, kk) +
                       bivar_LL_covar_partial_Nr(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_rr_2,
                                                 sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n,
                                                 corr_rr_Nr_2_n, alpha_exp, beta_exp, kk))) -
       KK_mean_vol_rad);

  return covar_Nr_KK_mvr_coefA;
}
inline real_t covar_Nr_KK_mvr_termB(real_t Nr_1, real_t Nr_2, real_t mu_rr_1, real_t mu_rr_2,
                                    real_t mu_Nr_1, real_t mu_Nr_2, real_t mu_rr_1_n,
                                    real_t mu_rr_2_n, real_t mu_Nr_1_n, real_t mu_Nr_2_n,
                                    real_t sigma_rr_1, real_t sigma_rr_2, real_t sigma_Nr_1,
                                    real_t sigma_Nr_2, real_t sigma_rr_1_n, real_t sigma_rr_2_n,
                                    real_t sigma_Nr_1_n, real_t sigma_Nr_2_n, real_t corr_rr_Nr_1_n,
                                    real_t corr_rr_Nr_2_n, real_t KK_mvr_coef, real_t mixt_frac,
                                    const parameters_KK& kk) {
  // Description:
  // This function partially calculates the covariance of N_r and KK mean
  // volume radius of rain drops (R_vr), which can be written as < N_r'R_vr' >.
  //
  // The equation for the covariance of N_r and KK rain drop mean volume radius
  // is:
  //
  // < N_r'R_vr' >
  // = KK_mvr_coef
  //   ( a f_p(1)
  //     exp{ mu_rr_1_n alpha + mu_Nr_1_n ( beta + 1 )
  //          + (1/2) sigma_rr_1_n^2 alpha^2
  //          + (1/2) sigma_Nr_1_n^2 ( beta + 1 )^2
  //          + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
  //   + ( 1 - a ) f_p(2)
  //     exp{ mu_rr_2_n alpha + mu_Nr_2_n ( beta + 1 )
  //          + (1/2) sigma_rr_2_n^2 alpha^2
  //          + (1/2) sigma_Nr_2_n^2 ( beta + 1 )^2
  //          + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
  //   )
  //   - < R_vr > < N_r >;
  //
  // where "a" is the mixture fraction (weight of the 1st PDF component), and
  // f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
  // components, respectively.  The equation for < N_r'R_vr' > needs to be put
  // in the form:  < N_r'R_vr' > = coef_A < N_r > + term_B.  Since < N_r > can
  // be written as:
  //
  // < N_r > = a f_p(1) exp{ mu_Nr_1_n + (1/2) sigma_Nr_1_n^2 }
  //           + ( 1 - a ) f_p(2) exp{ mu_Nr_2_n + (1/2) sigma_Nr_2_n^2 }
  //         = a f_p(1) mu_Nr_1 + ( 1 - a ) f_p(2) mu_Nr_2
  //         = a Nr_1 + ( 1 - a ) Nr_2;
  //
  // the covariance < N_r'R_vr' > can be rewritten in terms of < N_r >:
  //
  // < N_r'R_vr' >
  // = ( KK_mvr_coef
  //     ( exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 alpha^2
  //            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
  //     + exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 alpha^2
  //            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
  //     )
  //     - < R_vr > ) < N_r >
  //   - KK_mvr_coef
  //     ( a Nr_1
  //       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 alpha^2
  //            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
  //     + ( 1 - a ) Nr_2
  //       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 alpha^2
  //            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
  //     ).
  //
  // This function calculates term B, which is given by:
  //
  // term_B
  // = - KK_mvr_coef
  //     ( a Nr_1
  //       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 alpha^2
  //            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
  //     + ( 1 - a ) Nr_2
  //       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 alpha^2
  //            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
  //     ).
  //
  // The covariance of sedimentation velocity (of rain drop concentration) and
  // rain drop concentration, < V_Nr'N_r' > is given by the equation:
  //
  // < V_Nr'N_r' > = 0.007 ( 10^6 microm/m ) < N_r'R_vr' >.
  //
  // It is advantageous to write the < V_Nr'N_r' > term in terms of < N_r >
  // because < V_Nr'N_r' > is found in the turbulent sedimentation term as part
  // of the predictive equation for < N_r >.  This allows the turbulent
  // sedimentation term to be handled semi-implicitly in the < N_r > predictive
  // equation, increasing the numerical stability of the solution.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_Nr_KK_mvr_termB;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_mvr_rr_exp;
  beta_exp  = KK_mvr_Nr_exp;

  // Calculate term B in the < N_r'R_vr' > covariance equation, where
  // < N_r'R_vr' > = coef_A <N_r> + term_B, and:
  //
  // term_B
  // = - KK_mvr_coef
  //     ( a Nr_1
  //       exp{ mu_rr_2_n alpha + mu_Nr_2_n beta
  //            + (1/2) sigma_rr_2_n^2 alpha^2
  //            + (1/2) sigma_Nr_2_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
  //     + ( 1 - a ) Nr_2
  //       exp{ mu_rr_1_n alpha + mu_Nr_1_n beta
  //            + (1/2) sigma_rr_1_n^2 alpha^2
  //            + (1/2) sigma_Nr_1_n^2 ( beta^2 + 2 beta )
  //            + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
  //     ).
  covar_Nr_KK_mvr_termB =
      ((-KK_mvr_coef) *
       (((mixt_frac * Nr_1) * bivar_LL_covar_partial_Nr(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n,
                                                        sigma_rr_2, sigma_Nr_2, sigma_rr_2_n,
                                                        sigma_Nr_2_n, corr_rr_Nr_2_n, alpha_exp,
                                                        beta_exp, kk)) +
        (((one - mixt_frac) * Nr_2) *
         bivar_LL_covar_partial_Nr(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_rr_1, sigma_Nr_1,
                                   sigma_rr_1_n, sigma_Nr_1_n, corr_rr_Nr_1_n, alpha_exp, beta_exp,
                                   kk))));

  return covar_Nr_KK_mvr_termB;
}
inline real_t covar_Nr_KK_mvr(real_t mu_rr_1, real_t mu_rr_2, real_t mu_Nr_1, real_t mu_Nr_2,
                              real_t mu_rr_1_n, real_t mu_rr_2_n, real_t mu_Nr_1_n,
                              real_t mu_Nr_2_n, real_t sigma_rr_1, real_t sigma_rr_2,
                              real_t sigma_Nr_1, real_t sigma_Nr_2, real_t sigma_rr_1_n,
                              real_t sigma_rr_2_n, real_t sigma_Nr_1_n, real_t sigma_Nr_2_n,
                              real_t corr_rr_Nr_1_n, real_t corr_rr_Nr_2_n, real_t Nrm,
                              real_t KK_mean_vol_rad, real_t KK_mvr_coef, real_t mixt_frac,
                              real_t precip_frac_1, real_t precip_frac_2, const parameters_KK& kk) {
  // Description:
  // This function calculates the covariance of N_r and KK mean volume radius
  // of rain drops (R_vr), which can be written as < N_r'R_vr' >.
  //
  // The equation for the covariance of N_r and KK rain drop mean volume radius
  // is:
  //
  // < N_r'R_vr' >
  // = KK_mvr_coef
  //   ( a f_p(1)
  //     exp{ mu_rr_1_n alpha + mu_Nr_1_n ( beta + 1 )
  //          + (1/2) sigma_rr_1_n^2 alpha^2
  //          + (1/2) sigma_Nr_1_n^2 ( beta + 1 )^2
  //          + corr_rr_Nr_1_n sigma_rr_1_n alpha sigma_Nr_1_n ( beta + 1 ) }
  //   + ( 1 - a ) f_p(2)
  //     exp{ mu_rr_2_n alpha + mu_Nr_2_n ( beta + 1 )
  //          + (1/2) sigma_rr_2_n^2 alpha^2
  //          + (1/2) sigma_Nr_2_n^2 ( beta + 1 )^2
  //          + corr_rr_Nr_2_n sigma_rr_2_n alpha sigma_Nr_2_n ( beta + 1 ) }
  //   )
  //   - < R_vr > < N_r >;
  //
  // where "a" is the mixture fraction (weight of the 1st PDF component), and
  // f_p(1) and f_p(2) are the precipitation fractions in the 1st and 2nd PDF
  // components, respectively.
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t covar_Nr_KK_mvr;
  real_t alpha_exp;
  real_t beta_exp;

  // Values of the KK exponents.
  alpha_exp = KK_mvr_rr_exp;
  beta_exp  = KK_mvr_Nr_exp;

  // Calculate the covariance of N_r and KK mean volume radius of rain drops.
  covar_Nr_KK_mvr =
      ((KK_mvr_coef * (((mixt_frac * precip_frac_1) *
                        bivar_LL_mean_eq(mu_rr_1, mu_Nr_1, mu_rr_1_n, mu_Nr_1_n, sigma_rr_1,
                                         sigma_Nr_1, sigma_rr_1_n, sigma_Nr_1_n, corr_rr_Nr_1_n,
                                         alpha_exp, (beta_exp + one), kk)) +
                       (((one - mixt_frac) * precip_frac_2) *
                        bivar_LL_mean_eq(mu_rr_2, mu_Nr_2, mu_rr_2_n, mu_Nr_2_n, sigma_rr_2,
                                         sigma_Nr_2, sigma_rr_2_n, sigma_Nr_2_n, corr_rr_Nr_2_n,
                                         alpha_exp, (beta_exp + one), kk)))) -
       (Nrm * KK_mean_vol_rad));

  return covar_Nr_KK_mvr;
}
inline real_t bivar_LL_covar_partial_rr(real_t mu_rr_i, real_t mu_Nr_i, real_t mu_rr_i_n,
                                        real_t mu_Nr_i_n, real_t sigma_rr_i, real_t sigma_Nr_i,
                                        real_t sigma_rr_i_n, real_t sigma_Nr_i_n,
                                        real_t corr_rr_Nr_i_n, real_t alpha_exp_in,
                                        real_t beta_exp_in, const parameters_KK& kk) {
  // Description:
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t bivar_LL_covar_partial_rr;
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

  // Calculate (partially) the covariance of the bivariate lognormal equation.
  if (((sigma_x1 <= x1_tol) && (sigma_x2 <= x2_tol))) {

    // The ith PDF component variance of both r_r and N_r is 0.
    bivar_LL_covar_partial_rr =
        static_cast<real_t>(bivar_LL_mean_const_all(mu_x1, mu_x2, alpha_exp, beta_exp, kk));

  } else {
    if ((sigma_x1 <= x1_tol)) {

      // The ith PDF component variance of r_r is 0.
      bivar_LL_covar_partial_rr = static_cast<real_t>(
          bivar_LL_mean_const_x1(mu_x1, mu_x2_n, sigma_x2_n, alpha_exp, beta_exp, kk));

    } else {
      if ((sigma_x2 <= x2_tol)) {

        // The ith PDF component variance of N_r is 0.
        bivar_LL_covar_partial_rr = static_cast<real_t>(
            bivar_LL_covar_const_x2_partial(mu_x1_n, mu_x2, sigma_x1_n, alpha_exp, beta_exp, kk));

      } else {

        // All fields vary in the ith PDF component.
        bivar_LL_covar_partial_rr = static_cast<real_t>(bivar_LL_covar_partial(
            mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, rho_x1x2_n, alpha_exp, beta_exp, kk));
      }
    }
  }

  return bivar_LL_covar_partial_rr;
}
inline real_t bivar_LL_covar_partial_Nr(real_t mu_rr_i, real_t mu_Nr_i, real_t mu_rr_i_n,
                                        real_t mu_Nr_i_n, real_t sigma_rr_i, real_t sigma_Nr_i,
                                        real_t sigma_rr_i_n, real_t sigma_Nr_i_n,
                                        real_t corr_rr_Nr_i_n, real_t alpha_exp_in,
                                        real_t beta_exp_in, const parameters_KK& kk) {
  // Description:
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  // Local Variables
  real_t bivar_LL_covar_partial_Nr;
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

  // Normally, x1 is used for r_r and x2 is used for N_r.  Here, x1 will be
  // used for N_r and x2 will be used for r_r.  This means that the exponent
  // on N_r, beta_exp_in, will become alpha_exp, and the exponent on r_r,
  // alpha_exp_in, will become beta_exp.

  // Means for the ith PDF component.
  if ((beta_exp_in >= zero)) {
    mu_x1 = static_cast<double>(mu_Nr_i);
  } else {
    mu_x1 = static_cast<double>(((mu_Nr_i) > (Nr_tol) ? (mu_Nr_i) : (Nr_tol)));
  }
  if ((alpha_exp_in >= zero)) {
    mu_x2 = static_cast<double>(mu_rr_i);
  } else {
    mu_x2 = static_cast<double>(((mu_rr_i) > (rr_tol) ? (mu_rr_i) : (rr_tol)));
  }
  mu_x1_n = static_cast<double>(mu_Nr_i_n);
  mu_x2_n = static_cast<double>(mu_rr_i_n);

  // Standard deviations for the ith PDF component.
  sigma_x1   = static_cast<double>(sigma_Nr_i);
  sigma_x2   = static_cast<double>(sigma_rr_i);
  sigma_x1_n = static_cast<double>(sigma_Nr_i_n);
  sigma_x2_n = static_cast<double>(sigma_rr_i_n);

  // Correlations for the ith PDF component.
  rho_x1x2_n = static_cast<double>(corr_rr_Nr_i_n);

  // Exponents.
  alpha_exp = static_cast<double>(beta_exp_in);
  beta_exp  = static_cast<double>(alpha_exp_in);

  // Tolerance values.
  // When the standard deviation of a variable is below the tolerance values,
  // it is considered to be zero, and the variable is considered to have a
  // constant value.
  x1_tol = static_cast<double>(Nr_tol);
  x2_tol = static_cast<double>(rr_tol);

  // Calculate (partially) the covariance of the bivariate lognormal equation.
  if (((sigma_x1 <= x1_tol) && (sigma_x2 <= x2_tol))) {

    // The ith PDF component variance of both r_r and N_r is 0.
    bivar_LL_covar_partial_Nr =
        static_cast<real_t>(bivar_LL_mean_const_all(mu_x1, mu_x2, alpha_exp, beta_exp, kk));

  } else {
    if ((sigma_x1 <= x1_tol)) {

      // The ith PDF component variance of N_r is 0.
      bivar_LL_covar_partial_Nr = static_cast<real_t>(
          bivar_LL_mean_const_x1(mu_x1, mu_x2_n, sigma_x2_n, alpha_exp, beta_exp, kk));

    } else {
      if ((sigma_x2 <= x2_tol)) {

        // The ith PDF component variance of r_r is 0.
        bivar_LL_covar_partial_Nr = static_cast<real_t>(
            bivar_LL_covar_const_x2_partial(mu_x1_n, mu_x2, sigma_x1_n, alpha_exp, beta_exp, kk));

      } else {

        // All fields vary in the ith PDF component.
        bivar_LL_covar_partial_Nr = static_cast<real_t>(bivar_LL_covar_partial(
            mu_x1_n, mu_x2_n, sigma_x1_n, sigma_x2_n, rho_x1x2_n, alpha_exp, beta_exp, kk));
      }
    }
  }

  return bivar_LL_covar_partial_Nr;
}
inline double bivar_LL_covar_partial(double mu_x1_n, double mu_x2_n, double sigma_x1_n,
                                     double sigma_x2_n, double rho_x1x2_n, double alpha_exp,
                                     double beta_exp, const parameters_KK& kk) {
  // Description:
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double bivar_LL_covar_partial;

  // BFB-025: Preserve Native's alpha polynomial FMA and separate variance/cross
  // and mean groups before exp, including the swapped Nr orientation.
  const double variance_x1 =
      (one_half_dp * (sigma_x1_n * sigma_x1_n)) *
      Kokkos::fma(two_dp, alpha_exp, alpha_exp * alpha_exp);
  const double variance_pair =
      Kokkos::fma(one_half_dp * (sigma_x2_n * sigma_x2_n), beta_exp * beta_exp, variance_x1);
  const double cross_prefix =
      ((rho_x1x2_n * sigma_x1_n) * (alpha_exp + one_dp)) * sigma_x2_n;
  const double variance_cross = Kokkos::fma(cross_prefix, beta_exp, variance_pair);
  const double mean_pair = Kokkos::fma(mu_x2_n, beta_exp, mu_x1_n * alpha_exp);
  bivar_LL_covar_partial = Kokkos::exp(variance_cross + mean_pair);

  return bivar_LL_covar_partial;
}
inline double bivar_LL_covar_const_x2_partial(double mu_x1_n, double mu_x2, double sigma_x1_n,
                                              double alpha_exp, double beta_exp,
                                              const parameters_KK& kk) {
  // Description:
  // References:
  //-----------------------------------------------------------------------
  // Input Variables
  // Return Variable
  double bivar_LL_covar_const_x2_partial;

  bivar_LL_covar_const_x2_partial =
      (Kokkos::pow(mu_x2, beta_exp) *
       Kokkos::exp(((mu_x1_n * alpha_exp) + ((one_half_dp * (sigma_x1_n * sigma_x1_n)) *
                                             ((alpha_exp * alpha_exp) + (two_dp * alpha_exp))))));

  return bivar_LL_covar_const_x2_partial;
}
} // namespace clubb_kokkos::kk_engine
