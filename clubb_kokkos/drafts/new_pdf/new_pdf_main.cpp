// Setter/responder PDF orchestration; see the public header and integration notes.
#include "new_pdf_main.hpp"
#include "new_pdf.hpp"
#include "model_flags.hpp"
#include "parameter_indices.hpp"
namespace clubb_kokkos {
namespace {

template <class xmView, class xp2View, class sgn_wpxpView, class mixt_fracView,
          class exp_factor_xView, class max_Skx2_pos_Skx_sgn_wpxpView,
          class max_Skx2_neg_Skx_sgn_wpxpView, class SkxView, class mu_x_1View, class mu_x_2View,
          class sigma_x_1_sqdView, class sigma_x_2_sqdView, class coef_sigma_x_1_sqdView,
          class coef_sigma_x_2_sqdView, class F_xView, class min_F_xView, class max_F_xView>
void calc_responder_var(int nz, xmView xm, xp2View xp2, sgn_wpxpView sgn_wpxp,
                        mixt_fracView mixt_frac, real_t coef_spread_DG_means_x,
                        exp_factor_xView exp_factor_x,
                        max_Skx2_pos_Skx_sgn_wpxpView max_Skx2_pos_Skx_sgn_wpxp,
                        max_Skx2_neg_Skx_sgn_wpxpView max_Skx2_neg_Skx_sgn_wpxp, SkxView Skx,
                        mu_x_1View mu_x_1, mu_x_2View mu_x_2, sigma_x_1_sqdView sigma_x_1_sqd,
                        sigma_x_2_sqdView sigma_x_2_sqd, coef_sigma_x_1_sqdView coef_sigma_x_1_sqd,
                        coef_sigma_x_2_sqdView coef_sigma_x_2_sqd, F_xView F_x, min_F_xView min_F_x,
                        max_F_xView max_F_x);
template <class SkxView, class F_xView, class zeta_xView, class min_F_xView, class max_F_xView>
void calc_F_x_zeta_x_setter(int nz, SkxView Skx, real_t slope_coef_spread_DG_means_x,
                            real_t pdf_component_stdev_factor_x, real_t lambda_, F_xView F_x,
                            zeta_xView zeta_x, min_F_xView min_F_x, max_F_xView max_F_x);
KOKKOS_INLINE_FUNCTION void calc_F_x_zeta_x_setter(real_t Skx, real_t slope_coef_spread_DG_means_x,
                                                   real_t pdf_component_stdev_factor_x,
                                                   real_t lambda_, real_t& F_x, real_t& zeta_x,
                                                   real_t& min_F_x, real_t& max_F_x);
template <class exp_factor_xView, class min_F_xView, class max_F_xView>
auto calc_F_x_responder(int nz, real_t coef_spread_DG_means_x, exp_factor_xView exp_factor_x,
                        min_F_xView min_F_x, max_F_xView max_F_x)
    -> Kokkos::View<real_t*, Kokkos::LayoutLeft, typename exp_factor_xView::memory_space>;
KOKKOS_INLINE_FUNCTION real_t calc_F_x_responder(real_t coef_spread_DG_means_x, real_t exp_factor_x,
                                                 real_t min_F_x, real_t max_F_x);
} // namespace

void new_pdf_driver(int nz, int ngrdcol, Real2D wm, Real2D rtm, Real2D thlm, Real2D wp2,
                    Real2D rtp2, Real2D thlp2, Real2D Skw, Real2D wprtp, Real2D wpthlp,
                    Real2D rtpthlp, Real2D clubb_params, Real2D Skrt, Real2D Skthl, Real2D mu_w_1,
                    Real2D mu_w_2, Real2D mu_rt_1, Real2D mu_rt_2, Real2D mu_thl_1, Real2D mu_thl_2,
                    Real2D sigma_w_1_sqd, Real2D sigma_w_2_sqd, Real2D sigma_rt_1_sqd,
                    Real2D sigma_rt_2_sqd, Real2D sigma_thl_1_sqd, Real2D sigma_thl_2_sqd,
                    Real2D mixt_frac, implicit_coefs_terms& pdf_implicit_coefs_terms) {
  // Description:
  // Selects which variable is used to set the mixture fraction for the PDF
  // ("the setter") and which variables are handled after that mixture fraction
  // has been set ("the responders").  Traditionally, w has been used to set
  // the PDF.

  // References:
  //-----------------------------------------------------------------------

  // Arguments
  // wm: Mean of w (overall)                 [m/s]
  // rtm: Mean of rt (overall)                [kg/kg]
  // thlm: Mean of thl (overall)               [K]
  // wp2: Variance of w (overall)             [m^2/s^2]
  // rtp2: Variance of rt (overall)            [kg^2/kg^2]
  // thlp2: Variance of thl (overall)           [K^2]
  // Skw: Skewness of w (overall)             [-]
  // wprtp: Covariance of w and rt (overall)    [(m/s)kg/kg]
  // wpthlp: Covariance of w and thl (overall)   [(m/s)K]
  // rtpthlp: Covariance of rt and thl (overall)  [(kg/kg)K]
  // clubb_params: Array of CLUBB's tunable parameters    [units vary]
  // Skrt: Skewness of rt (overall)            [-]
  // Skthl: Skewness of thl (overall)           [-]
  // mu_w_1: Mean of w (1st PDF component)        [m/s]
  // mu_w_2: Mean of w (2nd PDF component)        [m/s]
  // mu_rt_1: Mean of rt (1st PDF component)       [kg/kg]
  // mu_rt_2: Mean of rt (2nd PDF component)       [kg/kg]
  // mu_thl_1: Mean of thl (1st PDF component)      [K]
  // mu_thl_2: Mean of thl (2nd PDF component)      [K]
  // sigma_w_1_sqd: Variance of w (1st PDF component)    [m^2/s^2]
  // sigma_w_2_sqd: Variance of w (2nd PDF component)    [m^2/s^2]
  // sigma_rt_1_sqd: Variance of rt (1st PDF component)   [kg^2/kg^2]
  // sigma_rt_2_sqd: Variance of rt (2nd PDF component)   [kg^2/kg^2]
  // sigma_thl_1_sqd: Variance of thl (1st PDF component)  [K^2]
  // sigma_thl_2_sqd: Variance of thl (2nd PDF component)  [K^2]
  // mixt_frac: Mixture fraction                     [-]
  // pdf_implicit_coefs_terms: Implicit coefs / explicit terms [units vary]

  //------------------------ Local Variables ------------------------
  // Output only for recording statistics.
  // Output functionality deprecated, see https://github.com/larson-group/clubb/issues/1176
  Real2D F_w("new_pdf_driver::F_w", ngrdcol,
             nz); // Parameter for the spread of the PDF component means of w    [-]
  Real2D F_rt("new_pdf_driver::F_rt", ngrdcol,
              nz); // Parameter for the spread of the PDF component means of rt   [-]
  Real2D F_thl("new_pdf_driver::F_thl", ngrdcol,
               nz); // Parameter for the spread of the PDF component means of thl  [-]
  Real2D min_F_w("new_pdf_driver::min_F_w", ngrdcol,
                 nz); // Minimum allowable value of parameter F_w      [-]
  Real2D max_F_w("new_pdf_driver::max_F_w", ngrdcol,
                 nz); // Maximum allowable value of parameter F_w      [-]
  Real2D min_F_rt("new_pdf_driver::min_F_rt", ngrdcol,
                  nz); // Minimum allowable value of parameter F_rt     [-]
  Real2D max_F_rt("new_pdf_driver::max_F_rt", ngrdcol,
                  nz); // Maximum allowable value of parameter F_rt     [-]
  Real2D min_F_thl("new_pdf_driver::min_F_thl", ngrdcol,
                   nz); // Minimum allowable value of parameter F_thl    [-]
  Real2D max_F_thl("new_pdf_driver::max_F_thl", ngrdcol,
                   nz); // Maximum allowable value of parameter F_thl    [-]
  Real1D sigma_w_1("new_pdf_driver::sigma_w_1",
                   nz); // Standard deviation of w (1st PDF component)      [m/s]
  Real1D sigma_w_2("new_pdf_driver::sigma_w_2",
                   nz); // Standard deviation of w (2nd PDF component)      [m/s]
  Real1D sgn_wprtp("new_pdf_driver::sgn_wprtp",
                   nz); // Sign of the covariance of w and rt (overall)     [-]
  Real1D sgn_wpthlp("new_pdf_driver::sgn_wpthlp",
                    nz); // Sign of the covariance of w and thl (overall)    [-]
  Real1D sgn_wp2("new_pdf_driver::sgn_wp2",
                 nz); // Sign of the variance of w (overall); always pos. [-]
  Real1D coef_sigma_w_1_sqd("new_pdf_driver::coef_sigma_w_1_sqd",
                            nz); // sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>    [-]
  Real1D coef_sigma_w_2_sqd("new_pdf_driver::coef_sigma_w_2_sqd",
                            nz); // sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>    [-]
  Real1D coef_sigma_rt_1_sqd("new_pdf_driver::coef_sigma_rt_1_sqd",
                             nz); // sigma_rt_1^2 = coef_sigma_rt_1_sqd * <rt'^2> [-]
  Real1D coef_sigma_rt_2_sqd("new_pdf_driver::coef_sigma_rt_2_sqd",
                             nz); // sigma_rt_2^2 = coef_sigma_rt_2_sqd * <rt'^2> [-]
  Real1D coef_sigma_thl_1_sqd("new_pdf_driver::coef_sigma_thl_1_sqd",
                              nz); // sigma_thl_1^2=coef_sigma_thl_1_sqd*<thl'^2>  [-]
  Real1D coef_sigma_thl_2_sqd("new_pdf_driver::coef_sigma_thl_2_sqd",
                              nz); // sigma_thl_2^2=coef_sigma_thl_2_sqd*<thl'^2>  [-]
  Real1D max_Skx2_pos_Skx_sgn_wpxp("new_pdf_driver::max_Skx2_pos_Skx_sgn_wpxp",
                                   nz); // Maximum Skx^2 when Skx*sgn(<w'x'>) >= 0 [-]
  Real1D max_Skx2_neg_Skx_sgn_wpxp("new_pdf_driver::max_Skx2_neg_Skx_sgn_wpxp",
                                   nz);        // Maximum Skx^2 when Skx*sgn(<w'x'>) < 0  [-]
  Real1D zeta_w("new_pdf_driver::zeta_w", nz); // Parameter for the PDF component variances of w [-]
  real_t lambda_w; // Param. that increases or decreases Skw dependence  [-]
  Real1D exp_factor_rt("new_pdf_driver::exp_factor_rt",
                       nz); // Factor of the form 1 - exp{} that reduces F_rt   [-]
  Real1D exp_factor_thl("new_pdf_driver::exp_factor_thl",
                        nz); // Don't reduce F_thl by exp_factor_thl        [-]
  Real1D adj_corr_rt_thl("new_pdf_driver::adj_corr_rt_thl",
                         nz); // Adjusted (overall) correlation of rt and theta-l [-]
  Real1D coef_wp4_implicit("new_pdf_driver::coef_wp4_implicit",
                           nz); // <w'^4> = coef_wp4_implicit * <w'^2>^2       [-]
  Real1D coef_wprtp2_implicit("new_pdf_driver::coef_wprtp2_implicit",
                              nz); // <w'rt'^2> = coef_wprtp2_implicit*<rt'^2>  [m/s]
  Real1D coef_wpthlp2_implicit("new_pdf_driver::coef_wpthlp2_implicit",
                               nz); // <w'thl'^2>=coef_wpthlp2_implicit*<thl'^2> [m/s]
  Real1D coef_wp2rtp_implicit("new_pdf_driver::coef_wp2rtp_implicit",
                              nz); // Coefficient that is multiplied by <w'rt'>  [m/s]
  Real1D term_wp2rtp_explicit("new_pdf_driver::term_wp2rtp_explicit",
                              nz); // Term that is on the RHS          [m^2/s^2 kg/kg]
  Real1D coef_wp2thlp_implicit("new_pdf_driver::coef_wp2thlp_implicit",
                               nz); // Coef. that is multiplied by <w'thl'>      [m/s]
  Real1D term_wp2thlp_explicit("new_pdf_driver::term_wp2thlp_explicit",
                               nz); // Term that is on the RHS             [m^2/s^2 K]
  Real1D coef_wprtpthlp_implicit("new_pdf_driver::coef_wprtpthlp_implicit",
                                 nz); // Coef. that is multiplied by <rt'thl'>   [m/s]
  Real1D term_wprtpthlp_explicit("new_pdf_driver::term_wprtpthlp_explicit",
                                 nz); // Term that is on the RHS         [m/s(kg/kg)K]

  // Native profile helpers receive scalar tunables. Mirror once at this host boundary.
  const auto clubb_params_host = to_host_space(clubb_params);
  // Capture component Views rather than the owning derived type.
  const auto pdf_coef_wp4_implicit       = pdf_implicit_coefs_terms.coef_wp4_implicit;
  const auto pdf_coef_wprtp2_implicit    = pdf_implicit_coefs_terms.coef_wprtp2_implicit;
  const auto pdf_coef_wpthlp2_implicit   = pdf_implicit_coefs_terms.coef_wpthlp2_implicit;
  const auto pdf_coef_wp2rtp_implicit    = pdf_implicit_coefs_terms.coef_wp2rtp_implicit;
  const auto pdf_term_wp2rtp_explicit    = pdf_implicit_coefs_terms.term_wp2rtp_explicit;
  const auto pdf_coef_wp2thlp_implicit   = pdf_implicit_coefs_terms.coef_wp2thlp_implicit;
  const auto pdf_term_wp2thlp_explicit   = pdf_implicit_coefs_terms.term_wp2thlp_explicit;
  const auto pdf_coef_wprtpthlp_implicit = pdf_implicit_coefs_terms.coef_wprtpthlp_implicit;
  const auto pdf_term_wprtpthlp_explicit = pdf_implicit_coefs_terms.term_wprtpthlp_explicit;

  // ------------------------- Begin Code -------------------------

  for (int i = 0; i < ngrdcol; ++i) {

    // Calculate sgn( <w'rt'> ).
    Kokkos::parallel_for(
        "new_pdf_driver::profile_1", nz, KOKKOS_LAMBDA(const int k) {
          if (wprtp(i, k) >= zero) {
            sgn_wprtp(k) = one;
          } else {
            sgn_wprtp(k) = -one;
          }
        });

    // Calculate sgn( <w'thl'> ).
    Kokkos::parallel_for(
        "new_pdf_driver::profile_2", nz, KOKKOS_LAMBDA(const int k) {
          if (wpthlp(i, k) >= zero) {
            sgn_wpthlp(k) = one;
          } else {
            sgn_wpthlp(k) = -one;
          }
        });

    // Sign of the variance of w (overall), which is always positive.
    Kokkos::parallel_for(
        "new_pdf_driver::profile_3", nz, KOKKOS_LAMBDA(const int k) { sgn_wp2(k) = one; });

    lambda_w = real_t(0.5);

    // Calculate the adjusted (overall) correlation of rt and theta-l, and the
    // value of exp_factor_rt.
    Kokkos::parallel_for(
        "new_pdf_driver::profile_4", nz, KOKKOS_LAMBDA(const int k) {
          if (rtp2(i, k) >= ((rt_tol) * (rt_tol)) && thlp2(i, k) >= ((thl_tol) * (thl_tol))) {
            adj_corr_rt_thl(k) = rtpthlp(i, k) / Kokkos::sqrt(rtp2(i, k) * thlp2(i, k)) *
                                 sgn_wprtp(k) * sgn_wpthlp(k);
            adj_corr_rt_thl(k) =
                ((((adj_corr_rt_thl(k)) > (-max_mag_correlation) ? (adj_corr_rt_thl(k))
                                                                 : (-max_mag_correlation))) <
                         (max_mag_correlation)
                     ? (((adj_corr_rt_thl(k)) > (-max_mag_correlation) ? (adj_corr_rt_thl(k))
                                                                       : (-max_mag_correlation)))
                     : (max_mag_correlation));
            exp_factor_rt(k) =
                one - Kokkos::exp(-real_t(0.2) * Kokkos::pow((adj_corr_rt_thl(k) + one), 5));
          } else {
            adj_corr_rt_thl(k) = zero;
            exp_factor_rt(k)   = one;
          }
        });

    // The value of F_thl is not reduced by exp_factor_thl.
    Kokkos::parallel_for(
        "new_pdf_driver::profile_5", nz, KOKKOS_LAMBDA(const int k) { exp_factor_thl(k) = one; });

    // Vertical velocity, w, will always be the setter variable.
    calc_F_x_zeta_x_setter(nz, Kokkos::subview(Skw, i, Kokkos::ALL()),
                           clubb_params_host(i, islope_coef_spread_DG_means_w),
                           clubb_params_host(i, ipdf_component_stdev_factor_w), lambda_w,
                           Kokkos::subview(F_w, i, Kokkos::ALL()), zeta_w,
                           Kokkos::subview(min_F_w, i, Kokkos::ALL()),
                           Kokkos::subview(max_F_w, i, Kokkos::ALL()));

    // Calculate the PDF parameters, including mixture fraction, for the
    // setter variable, w.
    calc_setter_var_params(
        nz, Kokkos::subview(wm, i, Kokkos::ALL()), Kokkos::subview(wp2, i, Kokkos::ALL()),
        Kokkos::subview(Skw, i, Kokkos::ALL()), sgn_wp2, Kokkos::subview(F_w, i, Kokkos::ALL()),
        zeta_w, Kokkos::subview(mu_w_1, i, Kokkos::ALL()),
        Kokkos::subview(mu_w_2, i, Kokkos::ALL()), sigma_w_1, sigma_w_2,
        Kokkos::subview(mixt_frac, i, Kokkos::ALL()), coef_sigma_w_1_sqd, coef_sigma_w_2_sqd);

    Kokkos::parallel_for(
        "new_pdf_driver::profile_6", nz, KOKKOS_LAMBDA(const int k) {
          sigma_w_1_sqd(i, k) = ((sigma_w_1(k)) * (sigma_w_1(k)));
          sigma_w_2_sqd(i, k) = ((sigma_w_2(k)) * (sigma_w_2(k)));
        });

    // Calculate the upper limit on the magnitude of skewness for responding
    // variables.
    Kokkos::parallel_for(
        "new_pdf_driver::profile_7", nz, KOKKOS_LAMBDA(const int k) {
          max_Skx2_pos_Skx_sgn_wpxp(k) = four *
                                         (((one - mixt_frac(i, k))) * ((one - mixt_frac(i, k)))) /
                                         (mixt_frac(i, k) * (two - mixt_frac(i, k)));
        });

    Kokkos::parallel_for(
        "new_pdf_driver::profile_8", nz, KOKKOS_LAMBDA(const int k) {
          max_Skx2_neg_Skx_sgn_wpxp(k) = four * ((mixt_frac(i, k)) * (mixt_frac(i, k))) /
                                         (one - ((mixt_frac(i, k)) * (mixt_frac(i, k))));
        });

    // Calculate the PDF parameters for responder variable rt.
    calc_responder_var(
        nz, Kokkos::subview(rtm, i, Kokkos::ALL()), Kokkos::subview(rtp2, i, Kokkos::ALL()),
        sgn_wprtp, Kokkos::subview(mixt_frac, i, Kokkos::ALL()),
        clubb_params_host(i, icoef_spread_DG_means_rt), exp_factor_rt, max_Skx2_pos_Skx_sgn_wpxp,
        max_Skx2_neg_Skx_sgn_wpxp, Kokkos::subview(Skrt, i, Kokkos::ALL()),
        Kokkos::subview(mu_rt_1, i, Kokkos::ALL()), Kokkos::subview(mu_rt_2, i, Kokkos::ALL()),
        Kokkos::subview(sigma_rt_1_sqd, i, Kokkos::ALL()),
        Kokkos::subview(sigma_rt_2_sqd, i, Kokkos::ALL()), coef_sigma_rt_1_sqd, coef_sigma_rt_2_sqd,
        Kokkos::subview(F_rt, i, Kokkos::ALL()), Kokkos::subview(min_F_rt, i, Kokkos::ALL()),
        Kokkos::subview(max_F_rt, i, Kokkos::ALL()));

    // Calculate the PDF parameters for responder variable thl.
    calc_responder_var(
        nz, Kokkos::subview(thlm, i, Kokkos::ALL()), Kokkos::subview(thlp2, i, Kokkos::ALL()),
        sgn_wpthlp, Kokkos::subview(mixt_frac, i, Kokkos::ALL()),
        clubb_params_host(i, icoef_spread_DG_means_thl), exp_factor_thl, max_Skx2_pos_Skx_sgn_wpxp,
        max_Skx2_neg_Skx_sgn_wpxp, Kokkos::subview(Skthl, i, Kokkos::ALL()),
        Kokkos::subview(mu_thl_1, i, Kokkos::ALL()), Kokkos::subview(mu_thl_2, i, Kokkos::ALL()),
        Kokkos::subview(sigma_thl_1_sqd, i, Kokkos::ALL()),
        Kokkos::subview(sigma_thl_2_sqd, i, Kokkos::ALL()), coef_sigma_thl_1_sqd,
        coef_sigma_thl_2_sqd, Kokkos::subview(F_thl, i, Kokkos::ALL()),
        Kokkos::subview(min_F_thl, i, Kokkos::ALL()), Kokkos::subview(max_F_thl, i, Kokkos::ALL()));

    if (!l_explicit_turbulent_adv_wp3) {

      // Turbulent advection of <w'^3> is being handled implicitly.

      // <w'^4> = coef_wp4_implicit * <w'^2>^2.
      coef_wp4_implicit = calc_coef_wp4_implicit(nz, Kokkos::subview(mixt_frac, i, Kokkos::ALL()),
                                                 Kokkos::subview(F_w, i, Kokkos::ALL()),
                                                 coef_sigma_w_1_sqd, coef_sigma_w_2_sqd);

    } else {

      // Turbulent advection of <w'^3> is being handled explicitly.
      Kokkos::parallel_for(
          "new_pdf_driver::profile_9", nz,
          KOKKOS_LAMBDA(const int k) { coef_wp4_implicit(k) = zero; });
    }

    if (!l_explicit_turbulent_adv_xpyp) {

      // Turbulent advection of <rt'^2> and <thl'^2> is being handled
      // implicitly.  Turbulent advection of <rt'thl'> is being handled
      // semi-implicitly.

      // <w'rt'^2> = coef_wprtp2_implicit * <rt'^2>
      coef_wprtp2_implicit = calc_coef_wpxp2_implicit(
          nz, Kokkos::subview(wp2, i, Kokkos::ALL()), Kokkos::subview(rtp2, i, Kokkos::ALL()),
          Kokkos::subview(wprtp, i, Kokkos::ALL()), sgn_wprtp,
          Kokkos::subview(mixt_frac, i, Kokkos::ALL()), Kokkos::subview(F_w, i, Kokkos::ALL()),
          Kokkos::subview(F_rt, i, Kokkos::ALL()), coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
          coef_sigma_rt_1_sqd, coef_sigma_rt_2_sqd);

      // <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2>
      coef_wpthlp2_implicit = calc_coef_wpxp2_implicit(
          nz, Kokkos::subview(wp2, i, Kokkos::ALL()), Kokkos::subview(thlp2, i, Kokkos::ALL()),
          Kokkos::subview(wpthlp, i, Kokkos::ALL()), sgn_wpthlp,
          Kokkos::subview(mixt_frac, i, Kokkos::ALL()), Kokkos::subview(F_w, i, Kokkos::ALL()),
          Kokkos::subview(F_thl, i, Kokkos::ALL()), coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
          coef_sigma_thl_1_sqd, coef_sigma_thl_2_sqd);

      // <w'rt'thl'> = coef_wprtpthlp_implicit * <rt'thl'>
      //               + term_wprtpthlp_explicit
      calc_coefs_wpxpyp_semiimpl(
          nz, Kokkos::subview(wp2, i, Kokkos::ALL()), Kokkos::subview(rtp2, i, Kokkos::ALL()),
          Kokkos::subview(thlp2, i, Kokkos::ALL()), Kokkos::subview(wprtp, i, Kokkos::ALL()),
          Kokkos::subview(wpthlp, i, Kokkos::ALL()), sgn_wprtp, sgn_wpthlp,
          Kokkos::subview(mixt_frac, i, Kokkos::ALL()), Kokkos::subview(F_w, i, Kokkos::ALL()),
          Kokkos::subview(F_rt, i, Kokkos::ALL()), Kokkos::subview(F_thl, i, Kokkos::ALL()),
          coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_rt_1_sqd, coef_sigma_rt_2_sqd,
          coef_sigma_thl_1_sqd, coef_sigma_thl_2_sqd, coef_wprtpthlp_implicit,
          term_wprtpthlp_explicit);

    } else {

      // Turbulent advection of <rt'^2>, <thl'^2>, and <rt'thl'> is being
      // handled explicitly.
      Kokkos::parallel_for(
          "new_pdf_driver::profile_10", nz, KOKKOS_LAMBDA(const int k) {
            coef_wprtp2_implicit(k)    = zero;
            coef_wpthlp2_implicit(k)   = zero;
            coef_wprtpthlp_implicit(k) = zero;
            term_wprtpthlp_explicit(k) = zero;
          });
    }

    if (!l_explicit_turbulent_adv_wpxp) {

      // Turbulent advection of <w'rt'> and <w'thl'> is being handled
      // semi-implicitly.

      // <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
      calc_coefs_wp2xp_semiimpl(
          nz, Kokkos::subview(wp2, i, Kokkos::ALL()), Kokkos::subview(rtp2, i, Kokkos::ALL()),
          sgn_wprtp, Kokkos::subview(mixt_frac, i, Kokkos::ALL()),
          Kokkos::subview(F_w, i, Kokkos::ALL()), Kokkos::subview(F_rt, i, Kokkos::ALL()),
          coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_rt_1_sqd, coef_sigma_rt_2_sqd,
          coef_wp2rtp_implicit, term_wp2rtp_explicit);

      // <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
      calc_coefs_wp2xp_semiimpl(
          nz, Kokkos::subview(wp2, i, Kokkos::ALL()), Kokkos::subview(thlp2, i, Kokkos::ALL()),
          sgn_wpthlp, Kokkos::subview(mixt_frac, i, Kokkos::ALL()),
          Kokkos::subview(F_w, i, Kokkos::ALL()), Kokkos::subview(F_thl, i, Kokkos::ALL()),
          coef_sigma_w_1_sqd, coef_sigma_w_2_sqd, coef_sigma_thl_1_sqd, coef_sigma_thl_2_sqd,
          coef_wp2thlp_implicit, term_wp2thlp_explicit);

    } else {

      // Turbulent advection of <w'rt'> and <w'thl'> is being handled
      // explicitly.
      Kokkos::parallel_for(
          "new_pdf_driver::profile_11", nz, KOKKOS_LAMBDA(const int k) {
            coef_wp2rtp_implicit(k)  = zero;
            term_wp2rtp_explicit(k)  = zero;
            coef_wp2thlp_implicit(k) = zero;
            term_wp2thlp_explicit(k) = zero;
          });
    }

    // Pack the implicit coefficients and explicit terms into a single type
    // variable for output.
    Kokkos::parallel_for(
        "new_pdf_driver::profile_12", nz, KOKKOS_LAMBDA(const int k) {
          pdf_coef_wp4_implicit(i, k)       = coef_wp4_implicit(k);
          pdf_coef_wprtp2_implicit(i, k)    = coef_wprtp2_implicit(k);
          pdf_coef_wpthlp2_implicit(i, k)   = coef_wpthlp2_implicit(k);
          pdf_coef_wp2rtp_implicit(i, k)    = coef_wp2rtp_implicit(k);
          pdf_term_wp2rtp_explicit(i, k)    = term_wp2rtp_explicit(k);
          pdf_coef_wp2thlp_implicit(i, k)   = coef_wp2thlp_implicit(k);
          pdf_term_wp2thlp_explicit(i, k)   = term_wp2thlp_explicit(k);
          pdf_coef_wprtpthlp_implicit(i, k) = coef_wprtpthlp_implicit(k);
          pdf_term_wprtpthlp_explicit(i, k) = term_wprtpthlp_explicit(k);
        });
  }

  return;
}

namespace {
template <class xmView, class xp2View, class sgn_wpxpView, class mixt_fracView,
          class exp_factor_xView, class max_Skx2_pos_Skx_sgn_wpxpView,
          class max_Skx2_neg_Skx_sgn_wpxpView, class SkxView, class mu_x_1View, class mu_x_2View,
          class sigma_x_1_sqdView, class sigma_x_2_sqdView, class coef_sigma_x_1_sqdView,
          class coef_sigma_x_2_sqdView, class F_xView, class min_F_xView, class max_F_xView>
void calc_responder_var(int nz, xmView xm, xp2View xp2, sgn_wpxpView sgn_wpxp,
                        mixt_fracView mixt_frac, real_t coef_spread_DG_means_x,
                        exp_factor_xView exp_factor_x,
                        max_Skx2_pos_Skx_sgn_wpxpView max_Skx2_pos_Skx_sgn_wpxp,
                        max_Skx2_neg_Skx_sgn_wpxpView max_Skx2_neg_Skx_sgn_wpxp, SkxView Skx,
                        mu_x_1View mu_x_1, mu_x_2View mu_x_2, sigma_x_1_sqdView sigma_x_1_sqd,
                        sigma_x_2_sqdView sigma_x_2_sqd, coef_sigma_x_1_sqdView coef_sigma_x_1_sqd,
                        coef_sigma_x_2_sqdView coef_sigma_x_2_sqd, F_xView F_x, min_F_xView min_F_x,
                        max_F_xView max_F_x) {
  // Description:
  // This is the sub-driver for a responder variable.  The upper limits of the
  // magnitude of Skx are calculated, and Skx is clipped when its magnitude
  // exceeds the upper limits.  The limits of the F_x parameter are calculated,
  // and the value of F_x is set within those limits.  Then, the PDF parameters
  // for responder variable x are calculated.

  // References:
  //-----------------------------------------------------------------------

  // Native profile order is retained across the two complex-root solves.
  auto xm_device                        = to_execution_space(xm);
  auto xp2_device                       = to_execution_space(xp2);
  auto sgn_wpxp_device                  = to_execution_space(sgn_wpxp);
  auto mixt_frac_device                 = to_execution_space(mixt_frac);
  auto exp_factor_x_device              = to_execution_space(exp_factor_x);
  auto max_Skx2_pos_Skx_sgn_wpxp_device = to_execution_space(max_Skx2_pos_Skx_sgn_wpxp);
  auto max_Skx2_neg_Skx_sgn_wpxp_device = to_execution_space(max_Skx2_neg_Skx_sgn_wpxp);
  auto Skx_device                       = to_execution_space(Skx);
  auto mu_x_1_device                    = to_execution_space(mu_x_1, false);
  auto mu_x_2_device                    = to_execution_space(mu_x_2, false);
  auto sigma_x_1_sqd_device             = to_execution_space(sigma_x_1_sqd, false);
  auto sigma_x_2_sqd_device             = to_execution_space(sigma_x_2_sqd, false);
  auto coef_sigma_x_1_sqd_device        = to_execution_space(coef_sigma_x_1_sqd, false);
  auto coef_sigma_x_2_sqd_device        = to_execution_space(coef_sigma_x_2_sqd, false);
  auto F_x_device                       = to_execution_space(F_x, false);
  auto min_F_x_device                   = to_execution_space(min_F_x, false);
  auto max_F_x_device                   = to_execution_space(max_F_x, false);

  // Calculate the upper limit of the magnitude of Skx.
  Kokkos::parallel_for(
      "calc_responder_var::clip_Skx", nz, KOKKOS_LAMBDA(const int k) {
        if (Skx_device(k) * sgn_wpxp_device(k) >= zero) {
          if (((Skx_device(k)) * (Skx_device(k))) >= max_Skx2_pos_Skx_sgn_wpxp_device(k)) {
            if (Skx_device(k) >= zero) {
              Skx_device(k) = Kokkos::sqrt(real_t(0.99) * max_Skx2_pos_Skx_sgn_wpxp_device(k));
            } else {
              Skx_device(k) = -Kokkos::sqrt(real_t(0.99) * max_Skx2_pos_Skx_sgn_wpxp_device(k));
            }
          }
        } else {
          if (((Skx_device(k)) * (Skx_device(k))) >= max_Skx2_neg_Skx_sgn_wpxp_device(k)) {
            if (Skx_device(k) >= zero) {
              Skx_device(k) = Kokkos::sqrt(real_t(0.99) * max_Skx2_neg_Skx_sgn_wpxp_device(k));
            } else {
              Skx_device(k) = -Kokkos::sqrt(real_t(0.99) * max_Skx2_neg_Skx_sgn_wpxp_device(k));
            }
          }
        }
      });

  calc_limits_F_x_responder(nz, mixt_frac_device, Skx_device, sgn_wpxp_device,
                            max_Skx2_pos_Skx_sgn_wpxp_device, max_Skx2_neg_Skx_sgn_wpxp_device,
                            min_F_x_device, max_F_x_device);

  // F_x must have a value between min_F_x and max_F_x.
  F_x_device = calc_F_x_responder(nz, coef_spread_DG_means_x, exp_factor_x_device, min_F_x_device,
                                  max_F_x_device);

  calc_responder_params(nz, xm_device, xp2_device, Skx_device, sgn_wpxp_device, F_x_device,
                        mixt_frac_device, mu_x_1_device, mu_x_2_device, sigma_x_1_sqd_device,
                        sigma_x_2_sqd_device, coef_sigma_x_1_sqd_device, coef_sigma_x_2_sqd_device);

  from_execution_space(Skx, Skx_device);
  from_execution_space(mu_x_1, mu_x_1_device);
  from_execution_space(mu_x_2, mu_x_2_device);
  from_execution_space(sigma_x_1_sqd, sigma_x_1_sqd_device);
  from_execution_space(sigma_x_2_sqd, sigma_x_2_sqd_device);
  from_execution_space(coef_sigma_x_1_sqd, coef_sigma_x_1_sqd_device);
  from_execution_space(coef_sigma_x_2_sqd, coef_sigma_x_2_sqd_device);
  from_execution_space(F_x, F_x_device);
  from_execution_space(min_F_x, min_F_x_device);
  from_execution_space(max_F_x, max_F_x_device);
}

KOKKOS_INLINE_FUNCTION void calc_F_x_zeta_x_setter(real_t Skx, real_t slope_coef_spread_DG_means_x,
                                                   real_t pdf_component_stdev_factor_x,
                                                   real_t lambda_, real_t& F_x, real_t& zeta_x,
                                                   real_t& min_F_x, real_t& max_F_x) {
  // Description:
  // Calculates the values of F_x and zeta_x for the setter variable (which is
  // the variable that sets the mixture fraction).
  //
  // The value of F_x is calculated between 0 (min_F_x) and 1 (max_F_x).  The
  // equation is:
  //
  // F_x = max_F_x + ( min_F_x - max_F_x )
  //                 * exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x };
  //
  // which reduces to:
  //
  // F_x = 1 - exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x };
  //
  // where lambda > 0 and slope_coef_spread_DG_means_x > 0.  As |Skx| goes
  // toward 0, the value of F_x goes toward 0, and as |Skx| becomes large, the
  // value of F_x goes toward 1.  When slope_coef_spread_DG_means_x is small,
  // the value of F_x tends toward 1, and when slope_coef_spread_DG_means_x is
  // large, the value of F_x tends toward 0.  When lambda is small, the value
  // of F_x is less dependent on Skx, and when lambda is large, the value of
  // F_x is more dependent on Skx.
  //
  // Mathematically, this equation will always produce a value of F_x that
  // falls between min_F_x and max_F_x.  However, in order to prevent a value
  // of F_x from being calculated outside the bounds of min_F_x and max_F_x
  // owing to numerical underflow or loss of precision, this equation can be
  // rewritten as:
  //
  // F_x
  // = min_F_x * exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x }
  //   + max_F_x * ( 1 - exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x } ).
  //
  // The value of zeta_x used to adjust the PDF component standard devations:
  //
  // 1 + zeta_x = ( mixt_frac * sigma_x_1^2 )
  //              / ( ( 1 - mixt_frac ) * sigma_x_2^2 );
  //
  // where zeta_x > -1.  The sign of zeta_x is used to easily determine if
  // mixt_frac * sigma_x_1^2 is greater than ( 1 - mixt_frac ) * sigma_x_2^2
  // (when zeta_x is positive), mixt_frac * sigma_x_1^2 is less than
  // ( 1 - mixt_frac ) * sigma_x_2^2 (when zeta_x is negative), or
  // mixt_frac * sigma_x_1^2 is equal to ( 1 - mixt_frac ) * sigma_x_2^2 (when
  // zeta_x is 0).
  //
  // In order to allow for a tunable parameter that is the pure ratio of
  // mixt_frac * sigma_x_1^2 to ( 1 - mixt_frac ) * sigma_x_2^2, zeta_x is
  // related to the parameter pdf_component_stdev_factor_x, where:
  //
  // 1 + zeta_x = pdf_component_stdev_factor_x.

  // References:
  //-----------------------------------------------------------------------

  // Independent Native profile cell; nz is represented by the enclosing View extent.
  // Skx: Skewness of x (overall)              [-]
  // slope_coef_spread_DG_means_x: Slope coef: spread PDF comp. means x [-]
  // pdf_component_stdev_factor_x: Param.: PDF comp. standard devs.; x  [-]
  // lambda: Param. for Skx dependence            [-]
  // F_x: Parameter for the spread of the PDF component means of x  [-]
  // zeta_x: Parameter for the PDF component variances of x            [-]
  // min_F_x: Minimum allowable value of parameter F_x                  [-]
  // max_F_x: Maximum allowable value of parameter F_x                  [-]

  real_t exp_Skx_interp_factor; // Function to interp. between min. and max.   [-]

  // Set min_F_x to 0 and max_F_x to 1 for the setter variable.
  if (Kokkos::abs(Skx) > zero) {
    min_F_x = real_t(1.0e-3);
  } else {
    min_F_x = zero;
  }
  max_F_x = one;

  // F_x must have a value between min_F_x and max_F_x.
  exp_Skx_interp_factor =
      Kokkos::exp(-Kokkos::pow(Kokkos::abs(Skx), lambda_) / slope_coef_spread_DG_means_x);

  F_x = min_F_x * exp_Skx_interp_factor + max_F_x * (one - exp_Skx_interp_factor);

  // The value of zeta_x must be greater than -1.
  zeta_x = pdf_component_stdev_factor_x - one;

  return;
}

template <class SkxView, class F_xView, class zeta_xView, class min_F_xView, class max_F_xView>
void calc_F_x_zeta_x_setter(int nz, SkxView Skx, real_t slope_coef_spread_DG_means_x,
                            real_t pdf_component_stdev_factor_x, real_t lambda_, F_xView F_x,
                            zeta_xView zeta_x, min_F_xView min_F_x, max_F_xView max_F_x) {
  // Description:
  // Calculates the values of F_x and zeta_x for the setter variable (which is
  // the variable that sets the mixture fraction).
  // Profile adaptation of the device scalar formula and Native method above.

  // Preserve Native profile argument order; stage Views at the public memory boundary.
  auto Skx_device     = to_execution_space(Skx);
  auto F_x_device     = to_execution_space(F_x, false);
  auto zeta_x_device  = to_execution_space(zeta_x, false);
  auto min_F_x_device = to_execution_space(min_F_x, false);
  auto max_F_x_device = to_execution_space(max_F_x, false);
  Kokkos::parallel_for(
      "calc_F_x_zeta_x_setter", nz, KOKKOS_LAMBDA(const int k) {
        calc_F_x_zeta_x_setter(Skx_device(k), slope_coef_spread_DG_means_x,
                               pdf_component_stdev_factor_x, lambda_, F_x_device(k),
                               zeta_x_device(k), min_F_x_device(k), max_F_x_device(k));
      });
  from_execution_space(F_x, F_x_device);
  from_execution_space(zeta_x, zeta_x_device);
  from_execution_space(min_F_x, min_F_x_device);
  from_execution_space(max_F_x, max_F_x_device);
}

KOKKOS_INLINE_FUNCTION real_t calc_F_x_responder(real_t coef_spread_DG_means_x, real_t exp_factor_x,
                                                 real_t min_F_x, real_t max_F_x) {
  // Description:
  // Calculates the value of F_x as a tunable function between min_F_x and
  // max_F_x.
  //
  // The value of F_x is calculated between min_F_x and max_F_x.  The equation
  // is:
  //
  // F_x = min_F_x + ( max_F_x - min_F_x )
  //                 * coef_spread_DG_means_x * exp_factor_x;
  //
  // where 0 <= coef_spread_DG_means_x <= 1.  As coef_spread_DG_means_x
  // goes toward 0, the value of F_x goes toward min_F_x, and as
  // coef_spread_DG_means_x goes toward 1, the value of F_x goes toward
  // max_F_x.  The exp_factor_x is a factor of the form 1 - exp{ }.  The range
  // of values of exp_factor_x is 0 <= exp_factor_x <= 1.  Here, exp_factor_x
  // is used to reduce the value of F_x under special conditions.
  //
  // The equation for exp_factor_x depends on which responder variable is being
  // solved for.  For rt, the F_rt equation is:
  //
  // F_rt = min_F_rt + ( max_F_rt - min_F_rt )
  //                   * coef_spread_DG_means_rt * exp_factor_rt;
  //
  // where exp_factor_rt is given by:
  //
  // exp_factor_rt = 1 - exp{ -0.2 * ( adj_corr_rt_thl + 1 )^5 };
  //
  // and where the adjusted (overall) correlation of rt and theta-l, denoted
  // adj_corr_rt_thl, is given by:
  //
  // adj_corr_rt_thl = <rt'thl'> / ( sqrt( <rt'^2> ) * sqrt( <thl'^2> ) )
  //                   * sgn( <w'rt'> ) * sgn( <w'thl'> ).
  //
  // The range of values of the adjusted (overall) correlation of rt and
  // theta-l is -1 <= adj_corr_rt_thl <= 1.  The adjusted (overall) correlation
  // of rt and theta-l has the same absolute magnitude as the (overall)
  // correlation of rt and theta-l, but the sign of the adjusted correlation is
  // dependent on the agreement between the signs of <w'rt'> and <w'thl'>.
  // When all three covariances are in sign agreement, which is when:
  //
  // sgn( <w'rt'> ) * sgn( <w'thl'> ) * sgn( <rt'thl'> ) = 1,
  //
  // the value of adj_corr_rt_thl is positive.  However, when the three
  // covariances aren't consistent in sign, which is when:
  //
  // sgn( <w'rt'> ) * sgn( <w'thl'> ) * sgn( <rt'thl'> ) = -1,
  //
  // the value of adj_corr_rt_thl is negative.
  //
  // For theta-l, the F_thl equation is:
  //
  // F_thl = min_F_thl + ( max_F_thl - min_F_thl ) * coef_spread_DG_means_thl;
  //
  // where exp_factor_thl is set to 1 (1 - exp{-inf} = 1) because reducing
  // F_thl through the use of exp_factor_thl is not desired for theta-l.
  //
  // The direction of the two PDF component means of x (mu_x_1 and mu_x_2) with
  // respect to each other and the overall mean of x (<x>) is given by
  // sgn( <w'x'> ).  When sgn( <w'x'> ) = 1, ( mu_x_1 - <x> ) >= 0 and
  // ( mu_x_2 - <x> ) <= 0.  When sgn( <w'x'> ) = -1, ( mu_x_1 - <x> ) <= 0 and
  // ( mu_x_2 - <x> ) >= 0.  This helps to promote realizability of the PDF
  // component correlations (corr_w_x_1 and corr_w_x_2).
  //
  // The realizability of the PDF component correlations corr_w_rt_1,
  // corr_w_rt_2, corr_w_thl_1, and corr_w_thl_2 is promoted by using
  // sgn( <w'rt'> ) and sgn( <w'thl'> ), respectively to choose the direction
  // of the PDF component means of rt and theta-l.  However, when the three
  // covariances (<w'rt'>, <w'thl'>, and <rt'thl'>) aren't consistent in sign
  // (as shown above), it becomes difficult to keep the PDF component
  // correlations corr_rt_thl_1 and corr_rt_thl_2 in the realizable range.
  // However, in this situation, the above equation for exp_factor_rt brings
  // F_rt closer to min_F_rt, which promotes realizability for corr_rt_thl_1
  // and corr_rt_thl_2.
  //
  // Mathematically, this equation will always produce a value of F_x that
  // falls between min_F_x and max_F_x.  However, in order to prevent a value
  // of F_x from being calculated outside the bounds of min_F_x and max_F_x
  // owing to numerical underflow or loss of precision, this equation can be
  // rewritten as:
  //
  // F_x = min_F_x * ( 1 - coef_spread_DG_means_x * exp_factor_x )
  //       + max_F_x * coef_spread_DG_means_x * exp_factor_x.

  // References:
  //-----------------------------------------------------------------------

  // Independent Native profile cell; nz is represented by the enclosing View extent.
  // coef_spread_DG_means_x: Coef.: spread betw. PDF comp. means of x   [-]
  // exp_factor_x: Factor of the form 1 - exp{}; reduces F_x  [-]
  // min_F_x: Minimum allowable value of parameter F_x   [-]
  // max_F_x: Maximum allowable value of parameter F_x   [-]

  real_t F_x; // Parameter for the spread between the PDF component means of x [-]

  // F_x must have a value between min_F_x and max_F_x.
  F_x = min_F_x * (one - coef_spread_DG_means_x * exp_factor_x) +
        max_F_x * coef_spread_DG_means_x * exp_factor_x;

  return F_x;
}

template <class exp_factor_xView, class min_F_xView, class max_F_xView>
auto calc_F_x_responder(int nz, real_t coef_spread_DG_means_x, exp_factor_xView exp_factor_x,
                        min_F_xView min_F_x, max_F_xView max_F_x)
    -> Kokkos::View<real_t*, Kokkos::LayoutLeft, typename exp_factor_xView::memory_space> {
  // Description:
  // Calculates the value of F_x as a tunable function between min_F_x and
  // max_F_x.
  // Profile adaptation of the device scalar formula and Native method above.

  // Preserve Native profile argument order; stage Views at the public memory boundary.
  auto exp_factor_x_device = to_execution_space(exp_factor_x);
  auto min_F_x_device      = to_execution_space(min_F_x);
  auto max_F_x_device      = to_execution_space(max_F_x);
  Real1D result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_F_x_responder"), nz);
  Kokkos::parallel_for(
      "calc_F_x_responder", nz, KOKKOS_LAMBDA(const int k) {
        result(k) = calc_F_x_responder(coef_spread_DG_means_x, exp_factor_x_device(k),
                                       min_F_x_device(k), max_F_x_device(k));
      });
  if constexpr (std::is_same_v<typename exp_factor_xView::memory_space,
                               typename Real1D::memory_space>)
    return result;
  else {
    Kokkos::View<real_t*, Kokkos::LayoutLeft, typename exp_factor_xView::memory_space> output(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "calc_F_x_responder"), nz);
    from_execution_space(output, result);
    return output;
  }
}

} // namespace
} // namespace clubb_kokkos
