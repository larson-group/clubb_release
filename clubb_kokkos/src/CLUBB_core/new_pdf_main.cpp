// Setter/responder PDF orchestration; see the public header for field ownership.
#include "scratch_workspace.hpp"
#include "new_pdf_main.hpp"
#include "new_pdf.hpp"
#include "model_flags.hpp"
#include "parameter_indices.hpp"
namespace clubb_kokkos {
namespace {

KOKKOS_INLINE_FUNCTION void
calc_responder_var(real_t xm, real_t xp2, real_t sgn_wpxp, real_t mixt_frac,           // In
                   real_t coef_spread_DG_means_x, real_t exp_factor_x,                 // In
                   real_t max_Skx2_pos_Skx_sgn_wpxp, real_t max_Skx2_neg_Skx_sgn_wpxp, // In
                   real_t& Skx,                                                        // InOut
                   real_t& mu_x_1, real_t& mu_x_2,                                     // Out
                   real_t& sigma_x_1_sqd, real_t& sigma_x_2_sqd,                       // Out
                   real_t& coef_sigma_x_1_sqd, real_t& coef_sigma_x_2_sqd,             // Out
                   real_t& F_x, real_t& min_F_x, real_t& max_F_x                       // Out
);
KOKKOS_INLINE_FUNCTION void calc_F_x_zeta_x_setter(real_t Skx,                          // In
                                                   real_t slope_coef_spread_DG_means_x, // In
                                                   real_t pdf_component_stdev_factor_x, // In
                                                   real_t lambda_,                      // In
                                                   real_t& F_x, real_t& zeta_x,         // Out
                                                   real_t& min_F_x, real_t& max_F_x     // Out
);
KOKKOS_INLINE_FUNCTION real_t calc_F_x_responder(real_t coef_spread_DG_means_x, // In
                                                 real_t exp_factor_x,           // In
                                                 real_t min_F_x, real_t max_F_x // In
);
} // namespace

void new_pdf_driver(ScratchWorkspace& workspace, int nz, int ngrdcol, Real2D wm, Real2D rtm,
                    Real2D thlm, Real2D wp2,                        // In
                    Real2D rtp2, Real2D thlp2,                      // In
                    Real2D Skw,                                     // In
                    Real2D wprtp, Real2D wpthlp, Real2D rtpthlp,    // In
                    Real2D clubb_params,                            // In
                    Real2D Skrt, Real2D Skthl,                      // InOut
                    Real2D mu_w_1, Real2D mu_w_2,                   // Out
                    Real2D mu_rt_1, Real2D mu_rt_2,                 // Out
                    Real2D mu_thl_1, Real2D mu_thl_2,               // Out
                    Real2D sigma_w_1_sqd, Real2D sigma_w_2_sqd,     // Out
                    Real2D sigma_rt_1_sqd, Real2D sigma_rt_2_sqd,   // Out
                    Real2D sigma_thl_1_sqd, Real2D sigma_thl_2_sqd, // Out
                    Real2D mixt_frac,                               // Out
                    implicit_coefs_terms& pdf_implicit_coefs_terms  // InOut
) {
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
  auto workspace_scope = workspace.scope();
  // Keep only quantities consumed by later kernels in column-first scratch.
  // Setter bounds, standard deviations and responder bounds remain cell local.
  Real2D F_w = workspace_scope.view<Real2D>(
      "new_pdf_driver::F_w", ngrdcol,
      nz); // Parameter for the spread of the PDF component means of w    [-]
  Real2D F_rt = workspace_scope.view<Real2D>(
      "new_pdf_driver::F_rt", ngrdcol,
      nz); // Parameter for the spread of the PDF component means of rt   [-]
  Real2D F_thl = workspace_scope.view<Real2D>(
      "new_pdf_driver::F_thl", ngrdcol,
      nz); // Parameter for the spread of the PDF component means of thl  [-]
  Real2D sgn_wprtp =
      workspace_scope.view<Real2D>("new_pdf_driver::sgn_wprtp", ngrdcol,
                                   nz); // Sign of the covariance of w and rt (overall)     [-]
  Real2D sgn_wpthlp =
      workspace_scope.view<Real2D>("new_pdf_driver::sgn_wpthlp", ngrdcol,
                                   nz); // Sign of the covariance of w and thl (overall)    [-]
  Real2D exp_factor_rt =
      workspace_scope.view<Real2D>("new_pdf_driver::exp_factor_rt", ngrdcol,
                                   nz); // Factor of the form 1 - exp{} that reduces F_rt   [-]
  Real2D coef_sigma_w_1_sqd =
      workspace_scope.view<Real2D>("new_pdf_driver::coef_sigma_w_1_sqd", ngrdcol,
                                   nz); // sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>    [-]
  Real2D coef_sigma_w_2_sqd =
      workspace_scope.view<Real2D>("new_pdf_driver::coef_sigma_w_2_sqd", ngrdcol,
                                   nz); // sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>    [-]
  Real2D coef_sigma_rt_1_sqd =
      workspace_scope.view<Real2D>("new_pdf_driver::coef_sigma_rt_1_sqd", ngrdcol,
                                   nz); // sigma_rt_1^2 = coef_sigma_rt_1_sqd * <rt'^2> [-]
  Real2D coef_sigma_rt_2_sqd =
      workspace_scope.view<Real2D>("new_pdf_driver::coef_sigma_rt_2_sqd", ngrdcol,
                                   nz); // sigma_rt_2^2 = coef_sigma_rt_2_sqd * <rt'^2> [-]
  Real2D coef_sigma_thl_1_sqd =
      workspace_scope.view<Real2D>("new_pdf_driver::coef_sigma_thl_1_sqd", ngrdcol,
                                   nz); // sigma_thl_1^2=coef_sigma_thl_1_sqd*<thl'^2>  [-]
  Real2D coef_sigma_thl_2_sqd =
      workspace_scope.view<Real2D>("new_pdf_driver::coef_sigma_thl_2_sqd", ngrdcol,
                                   nz); // sigma_thl_2^2=coef_sigma_thl_2_sqd*<thl'^2>  [-]
  Real2D max_Skx2_pos_Skx_sgn_wpxp =
      workspace_scope.view<Real2D>("new_pdf_driver::max_Skx2_pos_Skx_sgn_wpxp", ngrdcol,
                                   nz); // Maximum Skx^2 when Skx*sgn(<w'x'>) >= 0 [-]
  Real2D max_Skx2_neg_Skx_sgn_wpxp =
      workspace_scope.view<Real2D>("new_pdf_driver::max_Skx2_neg_Skx_sgn_wpxp", ngrdcol,
                                   nz); // Maximum Skx^2 when Skx*sgn(<w'x'>) < 0  [-]

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

  // Fuse the covariance signs and adjusted correlation with their exponential.
  Kokkos::parallel_for(
      "new_pdf_driver::responder_preparation", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        real_t adj_corr_rt_thl;

        if (wprtp(i, k) >= zero) {
          sgn_wprtp(i, k) = one;
        } else {
          sgn_wprtp(i, k) = -one;
        }

        if (wpthlp(i, k) >= zero) {
          sgn_wpthlp(i, k) = one;
        } else {
          sgn_wpthlp(i, k) = -one;
        }

        if (rtp2(i, k) >= ((rt_tol) * (rt_tol)) && thlp2(i, k) >= ((thl_tol) * (thl_tol))) {
          adj_corr_rt_thl = rtpthlp(i, k) / Kokkos::sqrt(rtp2(i, k) * thlp2(i, k)) *
                            sgn_wprtp(i, k) * sgn_wpthlp(i, k);
          adj_corr_rt_thl =
              ((((adj_corr_rt_thl) > (-max_mag_correlation) ? (adj_corr_rt_thl)
                                                            : (-max_mag_correlation))) <
                       (max_mag_correlation)
                   ? (((adj_corr_rt_thl) > (-max_mag_correlation) ? (adj_corr_rt_thl)
                                                                  : (-max_mag_correlation)))
                   : (max_mag_correlation));
          exp_factor_rt(i, k) =
              one - Kokkos::exp(-real_t(0.2) * Kokkos::pow((adj_corr_rt_thl + one), 5));
        } else {
          adj_corr_rt_thl     = zero;
          exp_factor_rt(i, k) = one;
        }
      });

  // Vertical velocity sets the mixture fraction. Fuse its cell-local bounds,
  // component standard deviations and responding-variable skewness limits.
  Kokkos::parallel_for(
      "new_pdf_driver::setter", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        real_t zeta_w, min_F_w, max_F_w, sigma_w_1, sigma_w_2;
        calc_F_x_zeta_x_setter(Skw(i, k), clubb_params(i, islope_coef_spread_DG_means_w),
                               clubb_params(i, ipdf_component_stdev_factor_w), real_t(0.5),
                               F_w(i, k), zeta_w, min_F_w, max_F_w);
        // The sign of the overall w variance is always positive.
        calc_setter_var_params(wm(i, k), wp2(i, k), Skw(i, k), one, F_w(i, k), zeta_w, mu_w_1(i, k),
                               mu_w_2(i, k), sigma_w_1, sigma_w_2, mixt_frac(i, k),
                               coef_sigma_w_1_sqd(i, k), coef_sigma_w_2_sqd(i, k));
        sigma_w_1_sqd(i, k) = ((sigma_w_1) * (sigma_w_1));
        sigma_w_2_sqd(i, k) = ((sigma_w_2) * (sigma_w_2));

        max_Skx2_pos_Skx_sgn_wpxp(i, k) = four *
                                          (((one - mixt_frac(i, k))) * ((one - mixt_frac(i, k)))) /
                                          (mixt_frac(i, k) * (two - mixt_frac(i, k)));

        max_Skx2_neg_Skx_sgn_wpxp(i, k) = four * ((mixt_frac(i, k)) * (mixt_frac(i, k))) /
                                          (one - ((mixt_frac(i, k)) * (mixt_frac(i, k))));
      });

  // Calculate the PDF parameters for responder variable rt.
  Kokkos::parallel_for(
      "new_pdf_driver::calc_responder_var", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        real_t min_F_rt;
        real_t max_F_rt;
        calc_responder_var(rtm(i, k),                                 // In
                           rtp2(i, k), sgn_wprtp(i, k),               // In
                           mixt_frac(i, k),                           // In
                           clubb_params(i, icoef_spread_DG_means_rt), // In
                           exp_factor_rt(i, k),                       // In
                           max_Skx2_pos_Skx_sgn_wpxp(i, k),           // In
                           max_Skx2_neg_Skx_sgn_wpxp(i, k),           // In
                           Skrt(i, k),                                // InOut
                           mu_rt_1(i, k),                             // Out
                           mu_rt_2(i, k),                             // Out
                           sigma_rt_1_sqd(i, k),                      // Out
                           sigma_rt_2_sqd(i, k),                      // Out
                           coef_sigma_rt_1_sqd(i, k),                 // Out
                           coef_sigma_rt_2_sqd(i, k),                 // Out
                           F_rt(i, k),                                // Out
                           min_F_rt,                                  // Out
                           max_F_rt                                   // Out
        );
      });

  // Calculate the PDF parameters for responder variable thl.
  Kokkos::parallel_for(
      "new_pdf_driver::calc_responder_var", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        real_t min_F_thl;
        real_t max_F_thl;
        calc_responder_var(thlm(i, k),                                 // In
                           thlp2(i, k), sgn_wpthlp(i, k),              // In
                           mixt_frac(i, k),                            // In
                           clubb_params(i, icoef_spread_DG_means_thl), // In
                           one,                                        // In
                           max_Skx2_pos_Skx_sgn_wpxp(i, k),            // In
                           max_Skx2_neg_Skx_sgn_wpxp(i, k),            // In
                           Skthl(i, k),                                // InOut
                           mu_thl_1(i, k),                             // Out
                           mu_thl_2(i, k),                             // Out
                           sigma_thl_1_sqd(i, k),                      // Out
                           sigma_thl_2_sqd(i, k),                      // Out
                           coef_sigma_thl_1_sqd(i, k),                 // Out
                           coef_sigma_thl_2_sqd(i, k),                 // Out
                           F_thl(i, k),                                // Out
                           min_F_thl,                                  // Out
                           max_F_thl                                   // Out
        );
      });

  if (!l_explicit_turbulent_adv_wp3) {

    // Turbulent advection of <w'^3> is being handled implicitly.

    // <w'^4> = coef_wp4_implicit * <w'^2>^2.
    Kokkos::parallel_for(
        "new_pdf_driver::calc_coef_wp4_implicit", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i                 = cell % ngrdcol;
          const int k                 = cell / ngrdcol;
          pdf_coef_wp4_implicit(i, k) = calc_coef_wp4_implicit( // In
              mixt_frac(i, k),                                  // In
              F_w(i, k),                                        // In
              coef_sigma_w_1_sqd(i, k),                         // In
              coef_sigma_w_2_sqd(i, k)                          // In
          );
        });

  } else {

    // Turbulent advection of <w'^3> is being handled explicitly.
    Kokkos::parallel_for(
        "new_pdf_driver::profile_9", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i                 = cell % ngrdcol;
          const int k                 = cell / ngrdcol;
          pdf_coef_wp4_implicit(i, k) = zero;
        });
  }

  if (!l_explicit_turbulent_adv_xpyp) {

    // Turbulent advection of <rt'^2> and <thl'^2> is being handled
    // implicitly.  Turbulent advection of <rt'thl'> is being handled
    // semi-implicitly.

    // <w'rt'^2> = coef_wprtp2_implicit * <rt'^2>
    Kokkos::parallel_for(
        "new_pdf_driver::calc_coef_wpxp2_implicit", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i                    = cell % ngrdcol;
          const int k                    = cell / ngrdcol;
          pdf_coef_wprtp2_implicit(i, k) = calc_coef_wpxp2_implicit( // In
              wp2(i, k), rtp2(i, k),                                 // In
              wprtp(i, k), sgn_wprtp(i, k),                          // In
              mixt_frac(i, k),                                       // In
              F_w(i, k), F_rt(i, k),                                 // In
              coef_sigma_w_1_sqd(i, k),                              // In
              coef_sigma_w_2_sqd(i, k),                              // In
              coef_sigma_rt_1_sqd(i, k),                             // In
              coef_sigma_rt_2_sqd(i, k)                              // In
          );
        });

    // <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2>
    Kokkos::parallel_for(
        "new_pdf_driver::calc_coef_wpxp2_implicit", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i                     = cell % ngrdcol;
          const int k                     = cell / ngrdcol;
          pdf_coef_wpthlp2_implicit(i, k) = calc_coef_wpxp2_implicit( // In
              wp2(i, k), thlp2(i, k),                                 // In
              wpthlp(i, k), sgn_wpthlp(i, k),                         // In
              mixt_frac(i, k),                                        // In
              F_w(i, k), F_thl(i, k),                                 // In
              coef_sigma_w_1_sqd(i, k),                               // In
              coef_sigma_w_2_sqd(i, k),                               // In
              coef_sigma_thl_1_sqd(i, k),                             // In
              coef_sigma_thl_2_sqd(i, k)                              // In
          );
        });

    // <w'rt'thl'> = coef_wprtpthlp_implicit * <rt'thl'>
    //               + term_wprtpthlp_explicit
    Kokkos::parallel_for(
        "new_pdf_driver::calc_coefs_wpxpyp_semiimpl", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          calc_coefs_wpxpyp_semiimpl(wp2(i, k),                                       // In
                                     rtp2(i, k), thlp2(i, k),                         // In
                                     wprtp(i, k),                                     // In
                                     wpthlp(i, k), sgn_wprtp(i, k), sgn_wpthlp(i, k), // In
                                     mixt_frac(i, k),                                 // In
                                     F_w(i, k), F_rt(i, k),                           // In
                                     F_thl(i, k),                                     // In
                                     coef_sigma_w_1_sqd(i, k),                        // In
                                     coef_sigma_w_2_sqd(i, k),                        // In
                                     coef_sigma_rt_1_sqd(i, k),                       // In
                                     coef_sigma_rt_2_sqd(i, k),                       // In
                                     coef_sigma_thl_1_sqd(i, k),                      // In
                                     coef_sigma_thl_2_sqd(i, k),                      // In
                                     pdf_coef_wprtpthlp_implicit(i, k),               // Out
                                     pdf_term_wprtpthlp_explicit(i, k)                // Out
          );
        });

  } else {

    // Turbulent advection of <rt'^2>, <thl'^2>, and <rt'thl'> is being
    // handled explicitly.
    Kokkos::parallel_for(
        "new_pdf_driver::profile_10", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i                       = cell % ngrdcol;
          const int k                       = cell / ngrdcol;
          pdf_coef_wprtp2_implicit(i, k)    = zero;
          pdf_coef_wpthlp2_implicit(i, k)   = zero;
          pdf_coef_wprtpthlp_implicit(i, k) = zero;
          pdf_term_wprtpthlp_explicit(i, k) = zero;
        });
  }

  if (!l_explicit_turbulent_adv_wpxp) {

    // Turbulent advection of <w'rt'> and <w'thl'> is being handled
    // semi-implicitly.

    // <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
    Kokkos::parallel_for(
        "new_pdf_driver::calc_coefs_wp2xp_semiimpl", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          calc_coefs_wp2xp_semiimpl(wp2(i, k),                      // In
                                    rtp2(i, k), sgn_wprtp(i, k),    // In
                                    mixt_frac(i, k),                // In
                                    F_w(i, k),                      // In
                                    F_rt(i, k),                     // In
                                    coef_sigma_w_1_sqd(i, k),       // In
                                    coef_sigma_w_2_sqd(i, k),       // In
                                    coef_sigma_rt_1_sqd(i, k),      // In
                                    coef_sigma_rt_2_sqd(i, k),      // In
                                    pdf_coef_wp2rtp_implicit(i, k), // Out
                                    pdf_term_wp2rtp_explicit(i, k)  // Out
          );
        });

    // <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
    Kokkos::parallel_for(
        "new_pdf_driver::calc_coefs_wp2xp_semiimpl", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          calc_coefs_wp2xp_semiimpl(wp2(i, k),                       // In
                                    thlp2(i, k), sgn_wpthlp(i, k),   // In
                                    mixt_frac(i, k),                 // In
                                    F_w(i, k),                       // In
                                    F_thl(i, k),                     // In
                                    coef_sigma_w_1_sqd(i, k),        // In
                                    coef_sigma_w_2_sqd(i, k),        // In
                                    coef_sigma_thl_1_sqd(i, k),      // In
                                    coef_sigma_thl_2_sqd(i, k),      // In
                                    pdf_coef_wp2thlp_implicit(i, k), // Out
                                    pdf_term_wp2thlp_explicit(i, k)  // Out
          );
        });

  } else {

    // Turbulent advection of <w'rt'> and <w'thl'> is being handled
    // explicitly.
    Kokkos::parallel_for(
        "new_pdf_driver::profile_11", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i                     = cell % ngrdcol;
          const int k                     = cell / ngrdcol;
          pdf_coef_wp2rtp_implicit(i, k)  = zero;
          pdf_term_wp2rtp_explicit(i, k)  = zero;
          pdf_coef_wp2thlp_implicit(i, k) = zero;
          pdf_term_wp2thlp_explicit(i, k) = zero;
        });
  }

  return;
}

namespace {
KOKKOS_INLINE_FUNCTION void
calc_responder_var(real_t xm, real_t xp2, real_t sgn_wpxp, real_t mixt_frac,           // In
                   real_t coef_spread_DG_means_x, real_t exp_factor_x,                 // In
                   real_t max_Skx2_pos_Skx_sgn_wpxp, real_t max_Skx2_neg_Skx_sgn_wpxp, // In
                   real_t& Skx,                                                        // InOut
                   real_t& mu_x_1, real_t& mu_x_2,                                     // Out
                   real_t& sigma_x_1_sqd, real_t& sigma_x_2_sqd,                       // Out
                   real_t& coef_sigma_x_1_sqd, real_t& coef_sigma_x_2_sqd,             // Out
                   real_t& F_x, real_t& min_F_x, real_t& max_F_x                       // Out
) {
  // Description:
  // This is the sub-driver for a responder variable.  The upper limits of the
  // magnitude of Skx are calculated, and Skx is clipped when its magnitude
  // exceeds the upper limits.  The limits of the F_x parameter are calculated,
  // and the value of F_x is set within those limits.  Then, the PDF parameters
  // for responder variable x are calculated.

  // References:
  //-----------------------------------------------------------------------

  if (Skx * sgn_wpxp >= zero) {
    if (((Skx) * (Skx)) >= max_Skx2_pos_Skx_sgn_wpxp) {
      if (Skx >= zero) {
        Skx = Kokkos::sqrt(real_t(0.99) * max_Skx2_pos_Skx_sgn_wpxp);
      } else {
        Skx = -Kokkos::sqrt(real_t(0.99) * max_Skx2_pos_Skx_sgn_wpxp);
      }
    }
  } else {
    if (((Skx) * (Skx)) >= max_Skx2_neg_Skx_sgn_wpxp) {
      if (Skx >= zero) {
        Skx = Kokkos::sqrt(real_t(0.99) * max_Skx2_neg_Skx_sgn_wpxp);
      } else {
        Skx = -Kokkos::sqrt(real_t(0.99) * max_Skx2_neg_Skx_sgn_wpxp);
      }
    }
  }
  calc_limits_F_x_responder(mixt_frac, Skx, sgn_wpxp,                             // In
                            max_Skx2_pos_Skx_sgn_wpxp, max_Skx2_neg_Skx_sgn_wpxp, // In
                            min_F_x, max_F_x);                                    // Out
  // F_x must have a value between min_F_x and max_F_x.
  F_x = calc_F_x_responder(coef_spread_DG_means_x, exp_factor_x, min_F_x, max_F_x);
  calc_responder_params(xm, xp2, Skx, sgn_wpxp, F_x, mixt_frac,       // In
                        mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, // Out
                        coef_sigma_x_1_sqd, coef_sigma_x_2_sqd);      // Out
}

KOKKOS_INLINE_FUNCTION void calc_F_x_zeta_x_setter(real_t Skx,                          // In
                                                   real_t slope_coef_spread_DG_means_x, // In
                                                   real_t pdf_component_stdev_factor_x, // In
                                                   real_t lambda_,                      // In
                                                   real_t& F_x, real_t& zeta_x,         // Out
                                                   real_t& min_F_x, real_t& max_F_x     // Out
) {
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

KOKKOS_INLINE_FUNCTION real_t calc_F_x_responder(real_t coef_spread_DG_means_x, // In
                                                 real_t exp_factor_x,           // In
                                                 real_t min_F_x, real_t max_F_x // In
) {
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

} // namespace
} // namespace clubb_kokkos
