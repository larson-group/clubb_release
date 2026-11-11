// Hybrid PDF orchestration; see the public header for field and ownership details.
#include "scratch_workspace.hpp"
#include "new_hybrid_pdf_main.hpp"
#include "new_hybrid_pdf.hpp"
#include "Skx_module.hpp"
#include "constants_clubb.hpp"
#include "model_flags.hpp"
#include <iostream>
#include <stdexcept>
namespace clubb_kokkos {
namespace {

KOKKOS_INLINE_FUNCTION void calc_responder_driver(real_t xm, real_t xp2, real_t wpxp,
                                                  real_t wp2,                     // In
                                                  real_t mixt_frac, real_t F_w,   // In
                                                  real_t& Skx,                    // InOut
                                                  real_t& mu_x_1, real_t& mu_x_2, // Out
                                                  real_t& sigma_x_1_sqd,          // Out
                                                  real_t& sigma_x_2_sqd,          // Out
                                                  real_t& coef_sigma_x_1_sqd,     // Out
                                                  real_t& coef_sigma_x_2_sqd);
KOKKOS_INLINE_FUNCTION void
calc_F_w_zeta_w(real_t Skw, real_t wprtp, real_t wpthlp, real_t upwp, real_t vpwp, // In
                real_t wp2, real_t rtp2, real_t thlp2, real_t up2, real_t vp2,     // In
                real_t gamma_Skw_fnc,                                              // In
                real_t slope_coef_spread_DG_means_w,                               // In
                real_t pdf_component_stdev_factor_w,                               // In
                real_t max_corr_w_sclr_sqd,                                        // In
                real_t& F_w, real_t& zeta_w, real_t& min_F_w, real_t& max_F_w);
} // namespace

void new_hybrid_pdf_driver(ScratchWorkspace& workspace,                                   // InOut
                           int nz, int ngrdcol, int sclr_dim,                             // In
                           Real2D wm, Real2D rtm, Real2D thlm, Real2D um, Real2D vm,      // In
                           Real2D wp2, Real2D rtp2, Real2D thlp2, Real2D up2, Real2D vp2, // In
                           Real2D Skw, Real2D wprtp, Real2D wpthlp, Real2D upwp, Real2D vpwp, // In
                           Real3D sclrm, Real3D sclrp2, Real3D wpsclrp,                       // In
                           Real2D clubb_params,                                               // In
                           Real1D slope_coef_spread_DG_means_w,                               // In
                           Real1D pdf_component_stdev_factor_w,                               // In
                           Real2D Skrt, Real2D Skthl, Real2D Sku, Real2D Skv,
                           Real3D Sksclr,                                              // InOut
                           Real2D mu_w_1, Real2D mu_w_2,                               // Out
                           Real2D mu_rt_1, Real2D mu_rt_2,                             // Out
                           Real2D mu_thl_1, Real2D mu_thl_2,                           // Out
                           Real2D mu_u_1, Real2D mu_u_2, Real2D mu_v_1, Real2D mu_v_2, // Out
                           Real2D sigma_w_1_sqd, Real2D sigma_w_2_sqd,                 // Out
                           Real2D sigma_rt_1_sqd, Real2D sigma_rt_2_sqd,               // Out
                           Real2D sigma_thl_1_sqd, Real2D sigma_thl_2_sqd,             // Out
                           Real2D sigma_u_1_sqd, Real2D sigma_u_2_sqd,                 // Out
                           Real2D sigma_v_1_sqd, Real2D sigma_v_2_sqd,                 // Out
                           Real3D mu_sclr_1, Real3D mu_sclr_2,                         // Out
                           Real3D sigma_sclr_1_sqd, Real3D sigma_sclr_2_sqd,           // Out
                           Real2D mixt_frac, Real2D sigma_sqd_w,                       // Out
                           implicit_coefs_terms& pdf_implicit_coefs_terms) {
  // Description:
  // Calculate the PDF parameters for w (including mixture fraction), rt,
  // theta-l, and optionally, u, v, and passive scalar variables.

  // References:
  //-----------------------------------------------------------------------

  // Arguments
  // wm: Mean of w (overall)                 [m/s]
  // rtm: Mean of rt (overall)                [kg/kg]
  // thlm: Mean of thl (overall)               [K]
  // um: Mean of u (overall)                 [m/s]
  // vm: Mean of v (overall)                 [m/s]
  // wp2: Variance of w (overall)             [m^2/s^2]
  // rtp2: Variance of rt (overall)            [kg^2/kg^2]
  // thlp2: Variance of thl (overall)           [K^2]
  // up2: Variance of u (overall)             [m^2/s^2]
  // vp2: Variance of v (overall)             [m^2/s^2]
  // Skw: Skewness of w (overall)             [-]
  // wprtp: Covariance of w and rt (overall)    [(m/s)kg/kg]
  // wpthlp: Covariance of w and thl (overall)   [(m/s)K]
  // upwp: Covariance of u and w (overall)     [m^2/s^2]
  // vpwp: Covariance of v and w (overall)     [m^2/s^2]
  // sclrm: Mean of sclr (overall)                [units vary]
  // sclrp2: Variance of sclr (overall)            [(units vary)^2]
  // wpsclrp: Covariance of w and sclr (overall)    [(m/s)(units vary)]
  // clubb_params: Array of CLUBB tunable parameters [units vary]
  // Skrt: Skewness of rt (overall)            [-]
  // Skthl: Skewness of thl (overall)           [-]
  // Sku: Skewness of u (overall)             [-]
  // Skv: Skewness of v (overall)             [-]
  // Sksclr: Skewness of sclr (overall)         [-]
  // mu_w_1: Mean of w (1st PDF component)        [m/s]
  // mu_w_2: Mean of w (2nd PDF component)        [m/s]
  // mu_rt_1: Mean of rt (1st PDF component)       [kg/kg]
  // mu_rt_2: Mean of rt (2nd PDF component)       [kg/kg]
  // mu_thl_1: Mean of thl (1st PDF component)      [K]
  // mu_thl_2: Mean of thl (2nd PDF component)      [K]
  // mu_u_1: Mean of u (1st PDF component)        [m/s]
  // mu_u_2: Mean of u (2nd PDF component)        [m/s]
  // mu_v_1: Mean of v (1st PDF component)        [m/s]
  // mu_v_2: Mean of v (2nd PDF component)        [m/s]
  // sigma_w_1_sqd: Variance of w (1st PDF component)    [m^2/s^2]
  // sigma_w_2_sqd: Variance of w (2nd PDF component)    [m^2/s^2]
  // sigma_rt_1_sqd: Variance of rt (1st PDF component)   [kg^2/kg^2]
  // sigma_rt_2_sqd: Variance of rt (2nd PDF component)   [kg^2/kg^2]
  // sigma_thl_1_sqd: Variance of thl (1st PDF component)  [K^2]
  // sigma_thl_2_sqd: Variance of thl (2nd PDF component)  [K^2]
  // sigma_u_1_sqd: Variance of u (1st PDF component)    [m^2/s^2]
  // sigma_u_2_sqd: Variance of u (2nd PDF component)    [m^2/s^2]
  // sigma_v_1_sqd: Variance of v (1st PDF component)    [m^2/s^2]
  // sigma_v_2_sqd: Variance of v (2nd PDF component)    [m^2/s^2]
  // mu_sclr_1: Mean of sclr (1st PDF component)      [units vary]
  // mu_sclr_2: Mean of sclr (2nd PDF component)      [units vary]
  // sigma_sclr_1_sqd: Variance of sclr (1st PDF component)  [(un. vary)^2]
  // sigma_sclr_2_sqd: Variance of sclr (2nd PDF component)  [(un. vary)^2]
  // mixt_frac: Mixture fraction                     [-]
  // pdf_implicit_coefs_terms: Implicit coefs / explicit terms [units vary]

  //------------------------ Local Variables ------------------------
  auto workspace_scope = workspace.scope();
  Real2D F_w           = workspace_scope.view<Real2D>(
      "new_hybrid_pdf_driver::F_w", ngrdcol,
      nz); // Parameter for the spread of the PDF component means of w  [-]
  Real2D gamma_Skw_fnc = workspace_scope.view<Real2D>(
      "new_hybrid_pdf_driver::gamma_Skw_fnc", ngrdcol, nz); // Value of parameter gamma from Skw [-]
  Real2D coef_sigma_w_1_sqd =
      workspace_scope.view<Real2D>("new_hybrid_pdf_driver::coef_sigma_w_1_sqd", ngrdcol, nz); // Out
  Real2D coef_sigma_w_2_sqd =
      workspace_scope.view<Real2D>("new_hybrid_pdf_driver::coef_sigma_w_2_sqd", ngrdcol,
                                   nz); // sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>    [-]
  Real2D coef_sigma_rt_1_sqd =
      workspace_scope.view<Real2D>("new_hybrid_pdf_driver::coef_sigma_rt_1_sqd", ngrdcol, nz); // In
  Real2D coef_sigma_rt_2_sqd =
      workspace_scope.view<Real2D>("new_hybrid_pdf_driver::coef_sigma_rt_2_sqd", ngrdcol, nz); // In
  Real2D coef_sigma_thl_1_sqd = workspace_scope.view<Real2D>(
      "new_hybrid_pdf_driver::coef_sigma_thl_1_sqd", ngrdcol, nz); // In
  Real2D coef_sigma_thl_2_sqd = workspace_scope.view<Real2D>(
      "new_hybrid_pdf_driver::coef_sigma_thl_2_sqd", ngrdcol, nz); // In
  Real2D coef_sigma_u_1_sqd =
      workspace_scope.view<Real2D>("new_hybrid_pdf_driver::coef_sigma_u_1_sqd", ngrdcol, nz); // In
  Real2D coef_sigma_u_2_sqd =
      workspace_scope.view<Real2D>("new_hybrid_pdf_driver::coef_sigma_u_2_sqd", ngrdcol, nz); // In
  Real2D coef_sigma_v_1_sqd =
      workspace_scope.view<Real2D>("new_hybrid_pdf_driver::coef_sigma_v_1_sqd", ngrdcol, nz); // In
  Real2D coef_sigma_v_2_sqd =
      workspace_scope.view<Real2D>("new_hybrid_pdf_driver::coef_sigma_v_2_sqd", ngrdcol, nz); // In
  Real3D coef_sigma_sclr_1_sqd =
      workspace_scope.view<Real3D>("new_hybrid_pdf_driver::coef_sigma_sclr_1_sqd", ngrdcol, nz,
                                   sclr_dim); // Coefficient that is multiplied by <sclr'^2> [-]
  Real3D coef_sigma_sclr_2_sqd =
      workspace_scope.view<Real3D>("new_hybrid_pdf_driver::coef_sigma_sclr_2_sqd", ngrdcol, nz,
                                   sclr_dim); // Coefficient that is multiplied by <sclr'^2> [-]

  // Capture individual component Views, retaining the caller-owned UDT allocations.
  const auto pdf_coef_wp4_implicit         = pdf_implicit_coefs_terms.coef_wp4_implicit;
  const auto pdf_coef_wp2rtp_implicit      = pdf_implicit_coefs_terms.coef_wp2rtp_implicit;
  const auto pdf_term_wp2rtp_explicit      = pdf_implicit_coefs_terms.term_wp2rtp_explicit;
  const auto pdf_coef_wp2thlp_implicit     = pdf_implicit_coefs_terms.coef_wp2thlp_implicit;
  const auto pdf_term_wp2thlp_explicit     = pdf_implicit_coefs_terms.term_wp2thlp_explicit;
  const auto pdf_coef_wp2up_implicit       = pdf_implicit_coefs_terms.coef_wp2up_implicit;
  const auto pdf_term_wp2up_explicit       = pdf_implicit_coefs_terms.term_wp2up_explicit;
  const auto pdf_coef_wp2vp_implicit       = pdf_implicit_coefs_terms.coef_wp2vp_implicit;
  const auto pdf_term_wp2vp_explicit       = pdf_implicit_coefs_terms.term_wp2vp_explicit;
  const auto pdf_coef_wprtp2_implicit      = pdf_implicit_coefs_terms.coef_wprtp2_implicit;
  const auto pdf_term_wprtp2_explicit      = pdf_implicit_coefs_terms.term_wprtp2_explicit;
  const auto pdf_coef_wpthlp2_implicit     = pdf_implicit_coefs_terms.coef_wpthlp2_implicit;
  const auto pdf_term_wpthlp2_explicit     = pdf_implicit_coefs_terms.term_wpthlp2_explicit;
  const auto pdf_coef_wprtpthlp_implicit   = pdf_implicit_coefs_terms.coef_wprtpthlp_implicit;
  const auto pdf_term_wprtpthlp_explicit   = pdf_implicit_coefs_terms.term_wprtpthlp_explicit;
  const auto pdf_coef_wpup2_implicit       = pdf_implicit_coefs_terms.coef_wpup2_implicit;
  const auto pdf_term_wpup2_explicit       = pdf_implicit_coefs_terms.term_wpup2_explicit;
  const auto pdf_coef_wpvp2_implicit       = pdf_implicit_coefs_terms.coef_wpvp2_implicit;
  const auto pdf_term_wpvp2_explicit       = pdf_implicit_coefs_terms.term_wpvp2_explicit;
  const auto pdf_coef_wp2sclrp_implicit    = pdf_implicit_coefs_terms.coef_wp2sclrp_implicit;
  const auto pdf_term_wp2sclrp_explicit    = pdf_implicit_coefs_terms.term_wp2sclrp_explicit;
  const auto pdf_coef_wpsclrp2_implicit    = pdf_implicit_coefs_terms.coef_wpsclrp2_implicit;
  const auto pdf_term_wpsclrp2_explicit    = pdf_implicit_coefs_terms.term_wpsclrp2_explicit;
  const auto pdf_coef_wprtpsclrp_implicit  = pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit;
  const auto pdf_term_wprtpsclrp_explicit  = pdf_implicit_coefs_terms.term_wprtpsclrp_explicit;
  const auto pdf_coef_wpthlpsclrp_implicit = pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit;
  const auto pdf_term_wpthlpsclrp_explicit = pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit;

  //------------------------ Begin Code ------------------------
  // Independent columns and levels share each kernel; scalar species remain
  // ordered within each cell, and the Native PDF phases retain their order.

  compute_gamma_Skw(nz, ngrdcol, Skw, clubb_params, gamma_Skw_fnc);

  // Fuse setter preparation with its cell-local PDF calculation. The scalar
  // correlation fold keeps the Native species order and comparison on ties.
  Kokkos::parallel_for(
      "new_hybrid_pdf_driver::setter", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i                = cell % ngrdcol;
        const int k                = cell / ngrdcol;
        real_t max_corr_w_sclr_sqd = zero;
        for (int sclr = 0; sclr < sclr_dim; ++sclr) {
          if (wp2(i, k) * sclrp2(i, k, sclr) > zero) {
            const real_t corr_w_sclr_sqd =
                ((wpsclrp(i, k, sclr)) * (wpsclrp(i, k, sclr))) / (wp2(i, k) * sclrp2(i, k, sclr));
            max_corr_w_sclr_sqd =
                corr_w_sclr_sqd > max_corr_w_sclr_sqd ? corr_w_sclr_sqd : max_corr_w_sclr_sqd;
          }
        }
        real_t zeta_w, min_F_w, max_F_w;
        // Vertical velocity sets the mixture fraction.
        calc_F_w_zeta_w(Skw(i, k), wprtp(i, k), wpthlp(i, k), upwp(i, k), vpwp(i, k), wp2(i, k),
                        rtp2(i, k), thlp2(i, k), up2(i, k), vp2(i, k), gamma_Skw_fnc(i, k),
                        slope_coef_spread_DG_means_w(i), pdf_component_stdev_factor_w(i),
                        max_corr_w_sclr_sqd, F_w(i, k), zeta_w, min_F_w, max_F_w);
        real_t sigma_w_1, sigma_w_2;
        calculate_w_params(wm(i, k), wp2(i, k), Skw(i, k), F_w(i, k), zeta_w, mu_w_1(i, k),
                           mu_w_2(i, k), sigma_w_1, sigma_w_2, mixt_frac(i, k),
                           coef_sigma_w_1_sqd(i, k), coef_sigma_w_2_sqd(i, k));
        sigma_w_1_sqd(i, k) = ((sigma_w_1) * (sigma_w_1));
        sigma_w_2_sqd(i, k) = ((sigma_w_2) * (sigma_w_2));
      });

  // Check the whole batch before any responder work, matching the Native fatal phase.
  int invalid_mixture = 0;
  Kokkos::parallel_reduce(
      "new_hybrid_pdf_driver::mixture_check", ngrdcol * nz,
      KOKKOS_LAMBDA(const int cell, int& count) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        if (mixt_frac(i, k) < zero)
          ++count;
      },
      invalid_mixture);
  if (invalid_mixture) {
    std::cerr << "Mixture fraction cannot be calculated.\n";
    std::cerr << "The value of F_w must be greater than 0 when | Skw | > 0.\n";
    // Native ERROR STOP terminates this driver; the host exception reaches the runner.
    throw std::runtime_error("new_hybrid_pdf_driver: invalid mixture fraction");
  }

  // Calculate the PDF parameters for responder variable rt.
  Kokkos::parallel_for(
      "new_hybrid_pdf_driver::calc_responder_driver_rtm", ngrdcol * nz,
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        calc_responder_driver(rtm(i, k), rtp2(i, k), wprtp(i, k), wp2(i, k), mixt_frac(i, k),
                              F_w(i, k), Skrt(i, k), mu_rt_1(i, k), mu_rt_2(i, k),
                              sigma_rt_1_sqd(i, k), sigma_rt_2_sqd(i, k), coef_sigma_rt_1_sqd(i, k),
                              coef_sigma_rt_2_sqd(i, k));
      });

  // Calculate the PDF parameters for responder variable thl.
  Kokkos::parallel_for(
      "new_hybrid_pdf_driver::calc_responder_driver_thlm", ngrdcol * nz,
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        calc_responder_driver(thlm(i, k), thlp2(i, k), wpthlp(i, k), wp2(i, k), mixt_frac(i, k),
                              F_w(i, k), Skthl(i, k), mu_thl_1(i, k), mu_thl_2(i, k),
                              sigma_thl_1_sqd(i, k), sigma_thl_2_sqd(i, k),
                              coef_sigma_thl_1_sqd(i, k), coef_sigma_thl_2_sqd(i, k));
      });

  // Calculate the PDF parameters for responder variable u.
  Kokkos::parallel_for(
      "new_hybrid_pdf_driver::calc_responder_driver_um", ngrdcol * nz,
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        calc_responder_driver(um(i, k), up2(i, k), upwp(i, k), wp2(i, k), mixt_frac(i, k),
                              F_w(i, k), Sku(i, k), mu_u_1(i, k), mu_u_2(i, k), sigma_u_1_sqd(i, k),
                              sigma_u_2_sqd(i, k), coef_sigma_u_1_sqd(i, k),
                              coef_sigma_u_2_sqd(i, k));
      });

  // Calculate the PDF parameters for responder variable v.
  Kokkos::parallel_for(
      "new_hybrid_pdf_driver::calc_responder_driver_vm", ngrdcol * nz,
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        calc_responder_driver(vm(i, k), vp2(i, k), vpwp(i, k), wp2(i, k), mixt_frac(i, k),
                              F_w(i, k), Skv(i, k), mu_v_1(i, k), mu_v_2(i, k), sigma_v_1_sqd(i, k),
                              sigma_v_2_sqd(i, k), coef_sigma_v_1_sqd(i, k),
                              coef_sigma_v_2_sqd(i, k));
      });

  // Calculate the PDF parameters for responder variable sclr.
  if (sclr_dim > 0) {

    Kokkos::parallel_for(
        "new_hybrid_pdf_driver::responding_scalars", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          for (int sclr = 0; sclr < sclr_dim; ++sclr) {

            calc_responder_driver(
                sclrm(i, k, sclr), sclrp2(i, k, sclr), wpsclrp(i, k, sclr), wp2(i, k),
                mixt_frac(i, k), F_w(i, k), Sksclr(i, k, sclr), mu_sclr_1(i, k, sclr),
                mu_sclr_2(i, k, sclr), sigma_sclr_1_sqd(i, k, sclr), sigma_sclr_2_sqd(i, k, sclr),
                coef_sigma_sclr_1_sqd(i, k, sclr), coef_sigma_sclr_2_sqd(i, k, sclr));
          }
        });
  }

  if (!l_explicit_turbulent_adv_wp3) {

    // Turbulent advection of <w'^3> is being handled implicitly.

    // <w'^4> = coef_wp4_implicit * <w'^2>^2.
    Kokkos::parallel_for(
        "new_hybrid_pdf_driver::set_coef_wp4_implicit", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell) {
          const int i                 = cell % ngrdcol;
          const int k                 = cell / ngrdcol;
          pdf_coef_wp4_implicit(i, k) = calculate_coef_wp4_implicit(
              mixt_frac(i, k), F_w(i, k), coef_sigma_w_1_sqd(i, k), coef_sigma_w_2_sqd(i, k));
        });

  } else {

    // Turbulent advection of <w'^3> is being handled explicitly.
    Kokkos::parallel_for(
        "new_hybrid_pdf_driver::set_coef_wp4_implicit", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell) {
          const int i                 = cell % ngrdcol;
          const int k                 = cell / ngrdcol;
          pdf_coef_wp4_implicit(i, k) = zero;
        });
  }

  if (!l_explicit_turbulent_adv_wpxp) {

    // Turbulent advection of <w'rt'>, <w'thl'>, <u'w'>, <v'w'>, and <w'sclr'>
    // are all being handled implicitly.

    Kokkos::parallel_for(
        "new_hybrid_pdf_driver::wp2xp_implicit", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          // <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'>

          pdf_coef_wp2rtp_implicit(i, k) =
              calc_coef_wp2xp_implicit(wp2(i, k), mixt_frac(i, k), F_w(i, k),
                                       coef_sigma_w_1_sqd(i, k), coef_sigma_w_2_sqd(i, k));

          // <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'>;
          // <w'^2 u'> = coef_wp2up_implicit * <u'w'>; and
          // <w'^2 v'> = coef_wp2vp_implicit * <v'w'>;
          // where each coef_wp2xp_implicit is the same as coef_wp2rtp_implicit.

          pdf_coef_wp2thlp_implicit(i, k) = pdf_coef_wp2rtp_implicit(i, k);
          pdf_coef_wp2up_implicit(i, k)   = pdf_coef_wp2rtp_implicit(i, k);
          pdf_coef_wp2vp_implicit(i, k)   = pdf_coef_wp2rtp_implicit(i, k);

          // <w'^2 sclr'> = coef_wp2sclrp_implicit * <w'sclr'>;
          // where each coef_wp2xp_implicit is the same as coef_wp2rtp_implicit.
          if (sclr_dim > 0) {

            for (int sclr = 0; sclr < sclr_dim; ++sclr) {
              pdf_coef_wp2sclrp_implicit(i, k, sclr) = pdf_coef_wp2rtp_implicit(i, k);
            }
          }
        });

  } else {

    // Turbulent advection of <w'rt'>, <w'thl'>, <u'w'>, <v'w'>, and <w'sclr'>
    // are all being handled explicitly.
    Kokkos::parallel_for(
        "new_hybrid_pdf_driver::set_coef_wp2rtp_implicit", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell) {
          const int i                     = cell % ngrdcol;
          const int k                     = cell / ngrdcol;
          pdf_coef_wp2rtp_implicit(i, k)  = zero;
          pdf_coef_wp2thlp_implicit(i, k) = zero;
          pdf_coef_wp2up_implicit(i, k)   = zero;
          pdf_coef_wp2vp_implicit(i, k)   = zero;
        });
    if (sclr_dim > 0) {
      Kokkos::parallel_for(
          "new_hybrid_pdf_driver::set_coef_wp2sclrp_implicit", ngrdcol * nz * sclr_dim,
          KOKKOS_LAMBDA(const int cell) {
            const int i                            = cell % ngrdcol;
            const int k                            = (cell / ngrdcol) % nz;
            const int sclr                         = cell / (ngrdcol * nz);
            pdf_coef_wp2sclrp_implicit(i, k, sclr) = zero;
          });
    }
  }

  if (!l_explicit_turbulent_adv_xpyp) {

    // Turbulent advection of <rt'^2>, <thl'^2>, <rt'thl'>, <u'^2>, <v'^2>,
    // <sclr'^2>, <sclr'rt'>, and <sclr'thl'> are all being handled
    // semi-implicitly.

    // <w'rt'^2> = coef_wprtp2_implicit * <rt'^2> + term_wprtp2_explicit
    Kokkos::parallel_for(
        "new_hybrid_pdf_driver::calc_coefs_wpxp2_semiimpl_wp2", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          calc_coefs_wpxp2_semiimpl(wp2(i, k), wprtp(i, k), mixt_frac(i, k), F_w(i, k),
                                    coef_sigma_rt_1_sqd(i, k), coef_sigma_rt_2_sqd(i, k),
                                    pdf_coef_wprtp2_implicit(i, k), pdf_term_wprtp2_explicit(i, k));
        });

    // <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2> + term_wprtp2_explicit
    Kokkos::parallel_for(
        "new_hybrid_pdf_driver::calc_coefs_wpxp2_semiimpl_wp2", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          calc_coefs_wpxp2_semiimpl(wp2(i, k), wpthlp(i, k), mixt_frac(i, k), F_w(i, k),
                                    coef_sigma_thl_1_sqd(i, k), coef_sigma_thl_2_sqd(i, k),
                                    pdf_coef_wpthlp2_implicit(i, k),
                                    pdf_term_wpthlp2_explicit(i, k));
        });

    // <w'rt'thl'> = coef_wprtpthlp_implicit * <rt'thl'>
    //               + term_wprtpthlp_explicit
    Kokkos::parallel_for(
        "new_hybrid_pdf_driver::calc_coefs_wpxpyp_semiimpl_wp2", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          calc_coefs_wpxpyp_semiimpl(wp2(i, k), wprtp(i, k), wpthlp(i, k), mixt_frac(i, k),
                                     F_w(i, k), coef_sigma_rt_1_sqd(i, k),
                                     coef_sigma_rt_2_sqd(i, k), coef_sigma_thl_1_sqd(i, k),
                                     coef_sigma_thl_2_sqd(i, k), pdf_coef_wprtpthlp_implicit(i, k),
                                     pdf_term_wprtpthlp_explicit(i, k));
        });

    // <w'u'^2> = coef_wpup2_implicit * <u'^2> + term_wpup2_explicit
    Kokkos::parallel_for(
        "new_hybrid_pdf_driver::calc_coefs_wpxp2_semiimpl_wp2", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          calc_coefs_wpxp2_semiimpl(wp2(i, k), upwp(i, k), mixt_frac(i, k), F_w(i, k),
                                    coef_sigma_u_1_sqd(i, k), coef_sigma_u_2_sqd(i, k),
                                    pdf_coef_wpup2_implicit(i, k), pdf_term_wpup2_explicit(i, k));
        });

    // <w'v'^2> = coef_wpvp2_implicit * <v'^2> + term_wpvp2_explicit
    Kokkos::parallel_for(
        "new_hybrid_pdf_driver::calc_coefs_wpxp2_semiimpl_wp2", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          calc_coefs_wpxp2_semiimpl(wp2(i, k), vpwp(i, k), mixt_frac(i, k), F_w(i, k),
                                    coef_sigma_v_1_sqd(i, k), coef_sigma_v_2_sqd(i, k),
                                    pdf_coef_wpvp2_implicit(i, k), pdf_term_wpvp2_explicit(i, k));
        });

    if (sclr_dim > 0) {

      Kokkos::parallel_for(
          "new_hybrid_pdf_driver::scalar_implicit_terms", ngrdcol * nz,
          KOKKOS_LAMBDA(const int cell) {
            const int i = cell % ngrdcol;
            const int k = cell / ngrdcol;
            for (int sclr = 0; sclr < sclr_dim; ++sclr) {

              // <w'sclr'^2> = coef_wpsclrp2_implicit * <sclr'^2>
              //               + term_wpsclrp2_explicit

              calc_coefs_wpxp2_semiimpl(
                  wp2(i, k), wpsclrp(i, k, sclr), mixt_frac(i, k), F_w(i, k),
                  coef_sigma_sclr_1_sqd(i, k, sclr), coef_sigma_sclr_2_sqd(i, k, sclr),
                  pdf_coef_wpsclrp2_implicit(i, k, sclr), pdf_term_wpsclrp2_explicit(i, k, sclr));

              // <w'rt'sclr'> = coef_wprtpsclrp_implicit * <sclr'rt'>
              //                + term_wprtpsclrp_explicit

              calc_coefs_wpxpyp_semiimpl(
                  wp2(i, k), wprtp(i, k), wpsclrp(i, k, sclr), mixt_frac(i, k), F_w(i, k),
                  coef_sigma_rt_1_sqd(i, k), coef_sigma_rt_2_sqd(i, k),
                  coef_sigma_sclr_1_sqd(i, k, sclr), coef_sigma_sclr_2_sqd(i, k, sclr),
                  pdf_coef_wprtpsclrp_implicit(i, k, sclr),
                  pdf_term_wprtpsclrp_explicit(i, k, sclr));

              // <w'thl'sclr'> = coef_wpthlpsclrp_implicit * <sclr'thl'>
              //                 + term_wpthlpsclrp_explicit

              calc_coefs_wpxpyp_semiimpl(
                  wp2(i, k), wpthlp(i, k), wpsclrp(i, k, sclr), mixt_frac(i, k), F_w(i, k),
                  coef_sigma_thl_1_sqd(i, k), coef_sigma_thl_2_sqd(i, k),
                  coef_sigma_sclr_1_sqd(i, k, sclr), coef_sigma_sclr_2_sqd(i, k, sclr),
                  pdf_coef_wpthlpsclrp_implicit(i, k, sclr),
                  pdf_term_wpthlpsclrp_explicit(i, k, sclr));
            }
          });
    }

  } else {

    // Turbulent advection of <rt'^2>, <thl'^2>, <rt'thl'>, <u'^2>, <v'^2>,
    // <sclr'^2>, <sclr'rt'>, and <sclr'thl'> are being handled explicitly.
    Kokkos::parallel_for(
        "new_hybrid_pdf_driver::set_coef_wprtp2_implicit", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell) {
          const int i                       = cell % ngrdcol;
          const int k                       = cell / ngrdcol;
          pdf_coef_wprtp2_implicit(i, k)    = zero;
          pdf_term_wprtp2_explicit(i, k)    = zero;
          pdf_coef_wpthlp2_implicit(i, k)   = zero;
          pdf_term_wpthlp2_explicit(i, k)   = zero;
          pdf_coef_wprtpthlp_implicit(i, k) = zero;
          pdf_term_wprtpthlp_explicit(i, k) = zero;
          pdf_coef_wpup2_implicit(i, k)     = zero;
          pdf_term_wpup2_explicit(i, k)     = zero;
          pdf_coef_wpvp2_implicit(i, k)     = zero;
          pdf_term_wpvp2_explicit(i, k)     = zero;
        });
    if (sclr_dim > 0) {
      Kokkos::parallel_for(
          "new_hybrid_pdf_driver::set_coef_wpsclrp2_implicit", ngrdcol * nz * sclr_dim,
          KOKKOS_LAMBDA(const int cell) {
            const int i                               = cell % ngrdcol;
            const int k                               = (cell / ngrdcol) % nz;
            const int sclr                            = cell / (ngrdcol * nz);
            pdf_coef_wpsclrp2_implicit(i, k, sclr)    = zero;
            pdf_term_wpsclrp2_explicit(i, k, sclr)    = zero;
            pdf_coef_wprtpsclrp_implicit(i, k, sclr)  = zero;
            pdf_term_wprtpsclrp_explicit(i, k, sclr)  = zero;
            pdf_coef_wpthlpsclrp_implicit(i, k, sclr) = zero;
            pdf_term_wpthlpsclrp_explicit(i, k, sclr) = zero;
          });
    }
  }

  // The calculation of skewness of rt, thl, u, v, and scalars is hard-wired
  // for use with the ADG1 code, which contains the variable sigma_sqd_w.
  // In order to use an equivalent expression for these skewnesses using the
  // new hybrid PDF (without doing more recoding), set the value of
  // sigma_sqd_w to 1 - F_w.
  // The hybrid closure has no explicit wp2*x contribution. Write the
  // constant directly, rather than allocating and copying zero profiles.
  Kokkos::parallel_for(
      "new_hybrid_pdf_driver::final_outputs", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i                     = cell % ngrdcol;
        const int k                     = cell / ngrdcol;
        sigma_sqd_w(i, k)               = one - F_w(i, k);
        pdf_term_wp2rtp_explicit(i, k)  = zero;
        pdf_term_wp2thlp_explicit(i, k) = zero;
        pdf_term_wp2up_explicit(i, k)   = zero;
        pdf_term_wp2vp_explicit(i, k)   = zero;
        for (int sclr = 0; sclr < sclr_dim; ++sclr)
          pdf_term_wp2sclrp_explicit(i, k, sclr) = zero;
      });

  return;
}

namespace {
// =============================================================================
KOKKOS_INLINE_FUNCTION void calc_responder_driver(real_t xm, real_t xp2, real_t wpxp,
                                                  real_t wp2,                     // In
                                                  real_t mixt_frac, real_t F_w,   // In
                                                  real_t& Skx,                    // InOut
                                                  real_t& mu_x_1, real_t& mu_x_2, // Out
                                                  real_t& sigma_x_1_sqd,          // Out
                                                  real_t& sigma_x_2_sqd,          // Out
                                                  real_t& coef_sigma_x_1_sqd,     // Out
                                                  real_t& coef_sigma_x_2_sqd) {
  // Description:
  // This is the sub-driver for a responder variable.  The limits of the range
  // of values of Skx that the PDF is able to represent are calculated, and
  // Skx is clipped when it exceeds the upper or lower limit.  Then, the PDF
  // parameters for responder variable x are calculated.

  // References:
  //-----------------------------------------------------------------------

  // Arguments retain Native order; scalar reference outputs implement elemental Out/InOut.
  // xm: Mean of x (overall)                          [units vary]
  // xp2: Variance of x (overall)                  [(units vary)^2]
  // wpxp: Covariance of w and x (overall)        [m/s (units vary)]
  // wp2: Variance of w (overall)                         [m^2/s^2]
  // mixt_frac: Mixture fraction                                      [-]
  // F_w: Parameter for the spread of the PDF comp. means of w  [-]
  // Skx: Skewness of x (overall)              [-]
  // mu_x_1: Mean of x (1st PDF component)        [units vary]
  // mu_x_2: Mean of x (2nd PDF component)        [units vary]
  // sigma_x_1_sqd: Variance of x (1st PDF component)    [(units vary)^2]
  // sigma_x_2_sqd: Variance of x (2nd PDF component)    [(units vary)^2]
  // coef_sigma_x_1_sqd: sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>    [-]
  // coef_sigma_x_2_sqd: sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>    [-]

  real_t corr_w_x;
  real_t min_Skx;
  real_t max_Skx;

  // Calculate the limits of Skx.
  if (F_w > zero) {

    // Calculate the overall correlation of w and x.
    if (wp2 * xp2 > zero) {

      // The overall correlation of w and x is:
      //
      // corr_w_x = <w'x'> / sqrt( <w'^2> * <x'^2> ).
      corr_w_x = wpxp / Kokkos::sqrt(wp2 * xp2);

    } else {

      // When <w'^2> = 0 or <x'^2> = 0, <w'x'> = 0.  The correlation of w
      // and x is undefined, however, since <w'x'> = 0, Skx = 0.  Setting
      // corr_w_x = 0 in this scenario will set max_Skx = min_Skx = 0.
      corr_w_x = zero;
    }

    if (wpxp >= zero) {

      // BFB-021: Native rounds this quotient before the corr^3/F_w^1.5 factors.
      const real_t skewness_factor =
          (one + mixt_frac) / Kokkos::sqrt(mixt_frac * (one - mixt_frac));
      min_Skx = skewness_factor * ((corr_w_x) * (corr_w_x) * (corr_w_x)) /
                    Kokkos::pow(F_w, three_halves) -
                Kokkos::sqrt(mixt_frac / (one - mixt_frac)) * three * corr_w_x / Kokkos::sqrt(F_w);

      max_Skx = (mixt_frac - two) / Kokkos::sqrt(mixt_frac * (one - mixt_frac)) *
                    ((corr_w_x) * (corr_w_x) * (corr_w_x)) / Kokkos::pow(F_w, three_halves) +
                Kokkos::sqrt((one - mixt_frac) / mixt_frac) * three * corr_w_x / Kokkos::sqrt(F_w);

    } else {

      min_Skx = (mixt_frac - two) / Kokkos::sqrt(mixt_frac * (one - mixt_frac)) *
                    ((corr_w_x) * (corr_w_x) * (corr_w_x)) / Kokkos::pow(F_w, three_halves) +
                Kokkos::sqrt((one - mixt_frac) / mixt_frac) * three * corr_w_x / Kokkos::sqrt(F_w);

      // BFB-021: Native rounds this quotient before the corr^3/F_w^1.5 factors.
      const real_t skewness_factor =
          (one + mixt_frac) / Kokkos::sqrt(mixt_frac * (one - mixt_frac));
      max_Skx = skewness_factor * ((corr_w_x) * (corr_w_x) * (corr_w_x)) /
                    Kokkos::pow(F_w, three_halves) -
                Kokkos::sqrt(mixt_frac / (one - mixt_frac)) * three * corr_w_x / Kokkos::sqrt(F_w);
    }

  } else {

    // When F_w = 0, <w'x'> = 0, and Skx = 0.
    min_Skx = zero;
    max_Skx = zero;
  }

  // Limit Skx so that min_Skx <= Skx <= max_Skx.
  if (Skx > max_Skx) {
    Skx = max_Skx;
  } else if (Skx < min_Skx) {
    Skx = min_Skx;
  }

  // Calculate the PDF parameters for responder variable x.
  calculate_responder_params(xm, xp2, Skx, wpxp, wp2, F_w, mixt_frac, mu_x_1, mu_x_2, sigma_x_1_sqd,
                             sigma_x_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd);

  return;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION void
calc_F_w_zeta_w(real_t Skw, real_t wprtp, real_t wpthlp, real_t upwp, real_t vpwp, // In
                real_t wp2, real_t rtp2, real_t thlp2, real_t up2, real_t vp2,     // In
                real_t gamma_Skw_fnc,                                              // In
                real_t slope_coef_spread_DG_means_w,                               // In
                real_t pdf_component_stdev_factor_w,                               // In
                real_t max_corr_w_sclr_sqd,                                        // In
                real_t& F_w, real_t& zeta_w, real_t& min_F_w, real_t& max_F_w) {
  // Description:
  // Calculates the values of F_w and zeta_w for w, which is the setter
  // variable (which is the variable that sets the mixture fraction).
  //
  // Based purely on the PDF of w and not considering restrictions caused by
  // other variables in the multivariate PDF, the value of F_w is calculated
  // between 0 (min_F_w) and 1 (max_F_w).  However, other variables in the PDF
  // place restrictions on the minimum value of F_w.  The range of acceptible
  // values for F_w is given by:
  //
  // max( <w'x'>^2 / ( <w'^2> * <x'^2> ), for all variables x ) <= F_w <= 1.
  //
  // Additionally, the value of F_w should still increase with an increasing
  // magnitude of Skw.  The value of F_w can be 0 only when Skw = 0 (as well as
  // all <w'x'> = 0, as well).  The F_w skewness equation used is:
  //
  // F_w = max_F_w + ( min_F_w - max_F_w )
  //                 * exp{ -|Skw|^lambda / slope_coef_spread_DG_means_w };
  //
  // where lambda > 0 and slope_coef_spread_DG_means_w > 0.  As |Skw| goes
  // toward 0, the value of F_w goes toward min_F_w, and as |Skw| becomes
  // large, the value of F_w goes toward max_F_w.  When the value of
  // slope_coef_spread_DG_means_w is small, the value of F_w tends toward
  // max_F_w, and when the value of slope_coef_spread_DG_means_w is large, the
  // value of F_w tends toward min_F_w.  When lambda is small, the value of F_w
  // is less dependent on Skw, and when lambda is large, the value of F_w is
  // more dependent on Skw.
  //
  // Mathematically, this equation will always produce a value of F_w that
  // falls between min_F_w and max_F_w.  However, in order to prevent a value
  // of F_w from being calculated outside the bounds of min_F_w and max_F_w
  // owing to numerical underflow or loss of precision, this equation can be
  // rewritten as:
  //
  // F_w
  // = min_F_w * exp{ -|Skw|^lambda / slope_coef_spread_DG_means_w }
  //   + max_F_w * ( 1 - exp{ -|Skw|^lambda / slope_coef_spread_DG_means_w } ).
  //
  // The value of zeta_w used to adjust the PDF component standard devations:
  //
  // 1 + zeta_w = ( mixt_frac * sigma_w_1^2 )
  //              / ( ( 1 - mixt_frac ) * sigma_w_2^2 );
  //
  // where zeta_w > -1.  The sign of zeta_w is used to easily determine if
  // mixt_frac * sigma_w_1^2 is greater than ( 1 - mixt_frac ) * sigma_w_2^2
  // (when zeta_w is positive), mixt_frac * sigma_w_1^2 is less than
  // ( 1 - mixt_frac ) * sigma_w_2^2 (when zeta_w is negative), or
  // mixt_frac * sigma_w_1^2 is equal to ( 1 - mixt_frac ) * sigma_w_2^2 (when
  // zeta_w is 0).
  //
  // In order to allow for a tunable parameter that is the pure ratio of
  // mixt_frac * sigma_w_1^2 to ( 1 - mixt_frac ) * sigma_w_2^2, zeta_w is
  // related to the parameter pdf_component_stdev_factor_w, where:
  //
  // 1 + zeta_w = pdf_component_stdev_factor_w.

  // References:
  //-----------------------------------------------------------------------

  // Arguments retain Native order; scalar reference outputs implement elemental Out/InOut.
  // Skw: Skewness of w (overall)              [-]
  // wprtp: Covariance (overall) of w and rt     [m/s (kg/kg)]
  // wpthlp: Covariance (overall) of w and thl    [m/s K]
  // upwp: Covariance (overall) of u and w      [m^2/s^2]
  // vpwp: Covariance (overall) of v and w      [m^2/s^2]
  // wp2: Variance (overall) of w              [m^2/s^2]
  // rtp2: Variance (overall) of rt             [(kg/kg)^2]
  // thlp2: Variance (overall) of thl            [K^2]
  // up2: Variance (overall) of u              [m^2/s^2]
  // vp2: Variance (overall) of v              [m^2/s^2]
  // gamma_Skw_fnc: Value of parameter gamma from tunable Skw function  [-]
  // max_corr_w_sclr_sqd: Max value of wpsclrp^2 / ( wp2 * sclrp2 )     [-]
  // F_w: Parameter for the spread of the PDF component means of w  [-]
  // zeta_w: Parameter for the PDF component variances of w            [-]
  // min_F_w: Minimum allowable value of parameter F_w                  [-]
  // max_F_w: Maximum allowable value of parameter F_w                  [-]

  (void)slope_coef_spread_DG_means_w; // Used only by the inactive exponential formulation.
  real_t corr_w_rt_sqd;
  real_t corr_w_thl_sqd;
  real_t corr_u_w_sqd;
  real_t corr_v_w_sqd;
  real_t max_corr_w_x_sqd;
  real_t zeta_w_star;

  // Find the maximum value of <w'x'>^2 / ( <w'^2> * <x'^2> ) for all
  // variables x that are Double Gaussian PDF responder variables.  This
  // includes rt and theta-l.  When l_predict_upwp_vpwp is enabled, u and v are
  // also calculated as part of the PDF, and they are included as well.
  // Additionally, when sclr_dim > 0, passive scalars (sclr) are also included.

  // Calculate the squared value of the correlation of w and rt.
  if (wp2 * rtp2 > zero) {
    corr_w_rt_sqd = ((wprtp) * (wprtp)) / (wp2 * rtp2);
  } else {
    corr_w_rt_sqd = zero;
  }

  // Calculate the squared value of the correlation of w and thl.
  if (wp2 * thlp2 > zero) {
    corr_w_thl_sqd = ((wpthlp) * (wpthlp)) / (wp2 * thlp2);
  } else {
    corr_w_thl_sqd = zero;
  }

  // Calculate the squared value of the correlation of u and w.
  if (up2 * wp2 > zero) {
    corr_u_w_sqd = ((upwp) * (upwp)) / (up2 * wp2);
  } else {
    corr_u_w_sqd = zero;
  }

  // Calculate the squared value of the correlation of v and w.
  if (vp2 * wp2 > zero) {
    corr_v_w_sqd = ((vpwp) * (vpwp)) / (vp2 * wp2);
  } else {
    corr_v_w_sqd = zero;
  }

  // Calculate the maximum value of corr_w_rt^2, corr_w_thl^2, corr_u_w^2, and
  // corr_v_w^2 in the calculation of max(corr_w_x^2).  Include the maximum
  // value of corr_w_sclr^2 (for all sclrs) when scalars are found.
  // Fold the Native variadic MAX in order, selecting its second operand on ties.
  max_corr_w_x_sqd = ((corr_w_rt_sqd) > (corr_w_thl_sqd) ? (corr_w_rt_sqd) : (corr_w_thl_sqd));
  max_corr_w_x_sqd = ((max_corr_w_x_sqd) > (corr_u_w_sqd) ? (max_corr_w_x_sqd) : (corr_u_w_sqd));
  max_corr_w_x_sqd = ((max_corr_w_x_sqd) > (corr_v_w_sqd) ? (max_corr_w_x_sqd) : (corr_v_w_sqd));
  max_corr_w_x_sqd =
      ((max_corr_w_x_sqd) > (max_corr_w_sclr_sqd) ? (max_corr_w_x_sqd) : (max_corr_w_sclr_sqd));

  max_corr_w_x_sqd =
      ((max_corr_w_x_sqd) < (((max_mag_correlation_flux) * (max_mag_correlation_flux)))
           ? (max_corr_w_x_sqd)
           : (((max_mag_correlation_flux) * (max_mag_correlation_flux))));

  // Calculate min_F_w and set max_F_w to 1.
  if (Kokkos::abs(Skw) > zero) {
    min_F_w = ((max_corr_w_x_sqd) > (real_t(1.0e-3)) ? (max_corr_w_x_sqd) : (real_t(1.0e-3)));
  } else {
    min_F_w = max_corr_w_x_sqd;
  }
  max_F_w = one;

  // F_w must have a value between min_F_w and max_F_w.
  // The exponential belongs only to the inactive formulation below; do not evaluate it.
  //    F_w = min_F_w * exp_Skw_interp_factor &
  //          + max_F_w * ( one - exp_Skw_interp_factor )

  // For now, use a formulation similar to what is used for ADG1.
  // Tunable parameter gamma.
  // When gamma goes to 0, the standard deviations of w in each PDF component
  // become small, and the spread between the two PDF component means of w
  // becomes large.  F_w goes to min_F_w.
  // When gamma goes to 1, the standard deviations of w in each PDF component
  // become large, and the spread between the two PDF component means of w
  // becomes small.  F_w goes to max_F_w.
  F_w = max_F_w - gamma_Skw_fnc * (max_F_w - min_F_w);

  // The value of zeta_w must be greater than -1.
  zeta_w_star = pdf_component_stdev_factor_w - one;

  // Make the PDF of w symmetric.  In other words, the PDF at a value of
  // positive skewness will look like a mirror image of the PDF at the
  // opposite value of negative skewness.
  if (Skw >= zero) {
    zeta_w = zeta_w_star;
  } else {
    zeta_w = -zeta_w_star / (zeta_w_star + one);
  }

  return;
}

} // namespace
} // namespace clubb_kokkos
