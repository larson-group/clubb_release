// ADG1, ADG2 and 3-D Luhar component PDF closures; see the public header.
#include "scratch_workspace.hpp"
#include "adg1_adg2_3d_luhar_pdf.hpp"
#include "constants_clubb.hpp"
#include "error_code.hpp"
#include "calc_roots.hpp"
#include <iostream>
#include <Kokkos_NumericTraits.hpp>
namespace clubb_kokkos {
namespace adg_pdf_detail {
KOKKOS_INLINE_FUNCTION void calc_Luhar_params(real_t Skx, real_t wpxp,      // In
                                              real_t xp2, real_t x_tol_sqd, // In
                                              real_t& mixt_frac, real_t& big_m, real_t& small_m);
KOKKOS_INLINE_FUNCTION void close_Luhar_pdf(real_t xm, real_t xp2, real_t mixt_frac,       // In
                                            real_t small_m, real_t wpxp, real_t x_tol_sqd, // In
                                            real_t& sigma_sqd_x_1, real_t& sigma_sqd_x_2,  // Out
                                            real_t& varnce_x_1, real_t& varnce_x_2,        // Out
                                            real_t& x_1_n, real_t& x_2_n, real_t& x_1, real_t& x_2);
} // namespace adg_pdf_detail
namespace {
void ADG1_ADG2_responder_params(int nz, int ngrdcol,                                       // In
                                Real2D xm, Real2D xp2, Real2D wp2, Real2D sqrt_wp2,        // In
                                Real2D wpxp, Real2D w_1_n, Real2D w_2_n, Real2D mixt_frac, // In
                                Real2D sigma_sqd_w, Real1D beta,                           // In
                                real_t x_tol,                                              // In
                                err_info_type& err_info,                                   // InOut
                                Real2D x_1, Real2D x_2, Real2D varnce_x_1,                 // Out
                                Real2D varnce_x_2, Real2D alpha_x);

KOKKOS_INLINE_FUNCTION void backsolve_Luhar_params(real_t Sk_max, real_t Skx,          // In
                                                   real_t big_m_max, real_t mixt_frac, // In
                                                   real_t& big_m_x, real_t& small_m_x);

KOKKOS_INLINE_FUNCTION real_t max_cubic_root(const Kokkos::Array<real_t, 1>& a_coef,
                                             const Kokkos::Array<real_t, 1>& b_coef,
                                             const Kokkos::Array<real_t, 1>& c_coef,
                                             const Kokkos::Array<real_t, 1>& d_coef);
} // namespace

void ADG1_pdf_driver(ScratchWorkspace& workspace,                                   // InOut
                     int nz, int ngrdcol, int sclr_dim, Real1D sclr_tol,            // In
                     Real2D wm, Real2D rtm, Real2D thlm, Real2D um, Real2D vm,      // In
                     Real2D wp2, Real2D rtp2, Real2D thlp2, Real2D up2, Real2D vp2, // In
                     Real2D Skw, Real2D wprtp, Real2D wpthlp, Real2D upwp, Real2D vpwp,
                     Real2D sqrt_wp2,                                                 // In
                     Real2D sigma_sqd_w, Real1D beta, Real1D mixt_frac_max_mag,       // In
                     Real3D sclrm, Real3D sclrp2, Real3D wpsclrp, bool l_scalar_calc, // In
                     err_info_type& err_info,                                         // InOut
                     Real2D w_1, Real2D w_2,                                          // Out
                     Real2D rt_1, Real2D rt_2,                                        // Out
                     Real2D thl_1, Real2D thl_2,                                      // Out
                     Real2D u_1, Real2D u_2, Real2D v_1, Real2D v_2,                  // Out
                     Real2D varnce_w_1, Real2D varnce_w_2,                            // Out
                     Real2D varnce_rt_1, Real2D varnce_rt_2,                          // Out
                     Real2D varnce_thl_1, Real2D varnce_thl_2,                        // Out
                     Real2D varnce_u_1, Real2D varnce_u_2,                            // Out
                     Real2D varnce_v_1, Real2D varnce_v_2,                            // Out
                     Real2D mixt_frac,                                                // Out
                     Real2D alpha_rt, Real2D alpha_thl,                               // Out
                     Real2D alpha_u, Real2D alpha_v,                                  // Out
                     Real3D sclr_1, Real3D sclr_2, Real3D varnce_sclr_1,              // Out
                     Real3D varnce_sclr_2, Real3D alpha_sclr) {
  // Calculate ADG1 component parameters from column/level moments, using
  // the vertical-velocity closure followed by rt, thl, u, v and scalar
  // responders. Return component means/variances and update native error state.

  // Local Variables
  auto workspace_scope = workspace.scope();
  Real2D w_1_n         = workspace_scope.view<Real2D>("ADG1_pdf_driver::w_1_n", ngrdcol,
                                                      nz); // Normalized component1 w mean [-]
  Real2D w_2_n         = workspace_scope.view<Real2D>("ADG1_pdf_driver::w_2_n", ngrdcol,
                                                      nz); // Normalized component2 w mean [-]

  // Calculate mixture fraction and PDF component means and variances of w.
  ADG1_w_closure(nz, ngrdcol, wm, wp2, Skw, sigma_sqd_w, // In
                 sqrt_wp2, mixt_frac_max_mag,            // In
                 err_info,                               // InOut
                 w_1, w_2, w_1_n, w_2_n,                 // Out
                 varnce_w_1, varnce_w_2, mixt_frac);     // Out

  if (clubb_at_least_debug_level_api(0)) {
    if (any_equal(err_info.err_code, clubb_fatal_error)) {
      return;
    }
  }

  // Calculate the PDF component means and variances of rt.
  ADG1_ADG2_responder_params(nz, ngrdcol,                    // In
                             rtm, rtp2, wp2, sqrt_wp2,       // In
                             wprtp, w_1_n, w_2_n, mixt_frac, // In
                             sigma_sqd_w, beta,              // In
                             rt_tol,                         // In
                             err_info,                       // InOut
                             rt_1, rt_2, varnce_rt_1,        // Out
                             varnce_rt_2, alpha_rt);         // Out

  // Calculate the PDF component means and variances of thl.
  ADG1_ADG2_responder_params(nz, ngrdcol,                     // In
                             thlm, thlp2, wp2, sqrt_wp2,      // In
                             wpthlp, w_1_n, w_2_n, mixt_frac, // In
                             sigma_sqd_w, beta,               // In
                             thl_tol,                         // In
                             err_info,                        // InOut
                             thl_1, thl_2, varnce_thl_1,      // Out
                             varnce_thl_2, alpha_thl);        // Out

  // Calculate the PDF component means and variances of u wind.
  ADG1_ADG2_responder_params(nz, ngrdcol,                   // In
                             um, up2, wp2, sqrt_wp2,        // In
                             upwp, w_1_n, w_2_n, mixt_frac, // In
                             sigma_sqd_w, beta,             // In
                             thl_tol,                       // In
                             err_info,                      // InOut
                             u_1, u_2, varnce_u_1,          // Out
                             varnce_u_2, alpha_u);          // Out

  // Calculate the PDF component means and variances of v wind.
  ADG1_ADG2_responder_params(nz, ngrdcol,                   // In
                             vm, vp2, wp2, sqrt_wp2,        // In
                             vpwp, w_1_n, w_2_n, mixt_frac, // In
                             sigma_sqd_w, beta,             // In
                             thl_tol,                       // In
                             err_info,                      // InOut
                             v_1, v_2, varnce_v_1,          // Out
                             varnce_v_2, alpha_v);          // Out

  // Calculate the PDF component means and variances of passive scalars.
  if (l_scalar_calc) {
    // Per-species tolerances cross from the execution space only at this
    // host dispatch boundary; scalar fields remain column-first device Views.
    const auto sclr_tol_host = to_host_space(workspace, sclr_tol);
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      ADG1_ADG2_responder_params(
          nz, ngrdcol, // In
          Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, sclr),                 // In
          wp2, sqrt_wp2, Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr), // In
          w_1_n, w_2_n, mixt_frac,                                                 // In
          sigma_sqd_w, beta,                                                       // In
          sclr_tol_host(sclr),                                                     // In
          err_info,                                                                // InOut
          Kokkos::subview(sclr_1, Kokkos::ALL, Kokkos::ALL, sclr),
          Kokkos::subview(sclr_2, Kokkos::ALL, Kokkos::ALL, sclr),        // Out
          Kokkos::subview(varnce_sclr_1, Kokkos::ALL, Kokkos::ALL, sclr), // Out
          Kokkos::subview(varnce_sclr_2, Kokkos::ALL, Kokkos::ALL, sclr), // Out
          Kokkos::subview(alpha_sclr, Kokkos::ALL, Kokkos::ALL, sclr));   // Out
    }
  }
}

void ADG2_pdf_driver(ScratchWorkspace& workspace,                        // InOut
                     int nz, int ngrdcol, int sclr_dim, Real1D sclr_tol, // In
                     Real2D wm, Real2D rtm, Real2D thlm, Real2D wp2, Real2D rtp2,
                     Real2D thlp2,                                                          // In
                     Real2D Skw, Real2D wprtp, Real2D wpthlp, Real2D sqrt_wp2, Real1D beta, // In
                     Real3D sclrm, Real3D sclrp2, Real3D wpsclrp, bool l_scalar_calc,       // In
                     err_info_type& err_info,                                               // InOut
                     Real2D w_1, Real2D w_2,                                                // Out
                     Real2D rt_1, Real2D rt_2,                                              // Out
                     Real2D thl_1, Real2D thl_2,                                            // Out
                     Real2D varnce_w_1, Real2D varnce_w_2,                                  // Out
                     Real2D varnce_rt_1, Real2D varnce_rt_2,                                // Out
                     Real2D varnce_thl_1, Real2D varnce_thl_2,                              // Out
                     Real2D mixt_frac,                                                      // Out
                     Real2D alpha_rt, Real2D alpha_thl,                                     // Out
                     Real2D sigma_sqd_w, Real3D sclr_1, Real3D sclr_2,                      // Out
                     Real3D varnce_sclr_1, Real3D varnce_sclr_2, Real3D alpha_sclr) {
  // Calculates the PDF parameters using the Analytic Double Gaussian 2 (ADG2)
  // PDF closure.

  // References:
  //-----------------------------------------------------------------------

  auto workspace_scope     = workspace.scope();
  Real2D w_1_n             = workspace_scope.view<Real2D>("ADG2_pdf_driver::w_1_n", ngrdcol, nz);
  Real2D w_2_n             = workspace_scope.view<Real2D>("ADG2_pdf_driver::w_2_n", ngrdcol, nz);
  const auto sclr_tol_host = to_host_space(workspace, sclr_tol);
  // Calculate the mixture fraction and the PDF component means and variances
  // of w. Independent cells retain Native calc/close/width order while sharing
  // work-item-local Luhar parameters; no host profile staging is needed.
  // Reproduce ADG2_w_closure using separate functions.
  Kokkos::parallel_for(
      "ADG2_pdf_driver::w_closure", ngrdcol * nz, KOKKOS_LAMBDA(const int index) {
        const int i = index % ngrdcol;
        const int k = index / ngrdcol;
        real_t small_m_w, big_m_w, sigma_sqd_w_1, sigma_sqd_w_2;
        adg_pdf_detail::calc_Luhar_params(Skw(i, k), wp2(i, k), wp2(i, k), w_tol_sqd,
                                          mixt_frac(i, k), big_m_w, small_m_w);
        adg_pdf_detail::close_Luhar_pdf(wm(i, k), wp2(i, k), mixt_frac(i, k), small_m_w, wp2(i, k),
                                        w_tol_sqd, sigma_sqd_w_1, sigma_sqd_w_2, varnce_w_1(i, k),
                                        varnce_w_2(i, k), w_1_n(i, k), w_2_n(i, k), w_1(i, k),
                                        w_2(i, k));
        // Overwrite sigma_sqd_w for consistency with ADG1.
        const real_t width = one / (one + small_m_w * small_m_w);
        sigma_sqd_w(i, k)  = width < real_t(0.99) ? width : real_t(0.99);
      });
  // Calculate the PDF component means and variances of rt.
  ADG1_ADG2_responder_params(nz, ngrdcol, rtm, rtp2, wp2, sqrt_wp2, wprtp, w_1_n, w_2_n, mixt_frac,
                             sigma_sqd_w, beta, rt_tol, err_info, rt_1, rt_2, varnce_rt_1,
                             varnce_rt_2, alpha_rt);
  // Calculate the PDF component means and variances of thl.
  ADG1_ADG2_responder_params(nz, ngrdcol, thlm, thlp2, wp2, sqrt_wp2, wpthlp, w_1_n, w_2_n,
                             mixt_frac, sigma_sqd_w, beta, thl_tol, err_info, thl_1, thl_2,
                             varnce_thl_1, varnce_thl_2, alpha_thl);
  // Calculate the PDF component means and variances of passive scalars.
  if (l_scalar_calc) {
    for (int sclr = 0; sclr < sclr_dim; ++sclr) {
      ADG1_ADG2_responder_params(nz, ngrdcol,
                                 Kokkos::subview(sclrm, Kokkos::ALL, Kokkos::ALL, sclr),
                                 Kokkos::subview(sclrp2, Kokkos::ALL, Kokkos::ALL, sclr), wp2,
                                 sqrt_wp2, Kokkos::subview(wpsclrp, Kokkos::ALL, Kokkos::ALL, sclr),
                                 w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta, sclr_tol_host(sclr),
                                 err_info, Kokkos::subview(sclr_1, Kokkos::ALL, Kokkos::ALL, sclr),
                                 Kokkos::subview(sclr_2, Kokkos::ALL, Kokkos::ALL, sclr),
                                 Kokkos::subview(varnce_sclr_1, Kokkos::ALL, Kokkos::ALL, sclr),
                                 Kokkos::subview(varnce_sclr_2, Kokkos::ALL, Kokkos::ALL, sclr),
                                 Kokkos::subview(alpha_sclr, Kokkos::ALL, Kokkos::ALL, sclr));
    }
  }
}

void Luhar_3D_pdf_driver(int nz, Real1D wm, Real1D rtm, Real1D thlm, Real1D wp2, Real1D rtp2,
                         Real1D thlp2, // In
                         Real1D Skw, Real1D Skrt, Real1D Skthl, Real1D wprtp,
                         Real1D wpthlp,                            // In
                         Real1D w_1, Real1D w_2,                   // Out
                         Real1D rt_1, Real1D rt_2,                 // Out
                         Real1D thl_1, Real1D thl_2,               // Out
                         Real1D varnce_w_1, Real1D varnce_w_2,     // Out
                         Real1D varnce_rt_1, Real1D varnce_rt_2,   // Out
                         Real1D varnce_thl_1, Real1D varnce_thl_2, // Out
                         Real1D mixt_frac) {
  // Calculates the PDF parameters using the 3D Luhar PDF closure.

  // References:
  //-----------------------------------------------------------------------

  // Native one-element profile calls become device scalar calls, retaining
  // each level's setter/responder order and the exact skewness tie branches.
  Kokkos::parallel_for(
      "Luhar_3D_pdf_driver", nz, KOKKOS_LAMBDA(const int k) {
        real_t w_1_n,        // Normalized mean of w, component1 [-]
            w_2_n,           // Normalized mean of w, component2 [-]
            rt_1_n,          // Normalized mean of rt, component1 [-]
            rt_2_n,          // Normalized mean of rt, component2 [-]
            thl_1_n,         // Normalized mean of thl, component1 [-]
            thl_2_n,         // Normalized mean of thl, component2 [-]
            small_m_w,       // Luhar m for w [-]
            big_m_w,         // Luhar M for w [-]
            small_m_rt,      // Luhar m for rt [-]
            big_m_rt,        // Luhar M for rt [-]
            small_m_thl,     // Luhar m for thl [-]
            big_m_thl,       // Luhar M for thl [-]
            sigma_sqd_w_1,   // Normalized width of w, component1 [-]
            sigma_sqd_w_2,   // Normalized width of w, component2 [-]
            sigma_sqd_rt_1,  // Normalized width of rt, component1 [-]
            sigma_sqd_rt_2,  // Normalized width of rt, component2 [-]
            sigma_sqd_thl_1, // Normalized width of thl, component1 [-]
            sigma_sqd_thl_2; // Normalized width of thl, component2 [-]
        if ((Kokkos::abs(Skw(k)) >= Kokkos::abs(Skthl(k))) &&
            (Kokkos::abs(Skw(k)) >= Kokkos::abs(Skrt(k)))) {
          // w has the greatest magnitude of skewness.
          // Solve for the w PDF
          adg_pdf_detail::calc_Luhar_params(Skw(k), wp2(k), wp2(k), w_tol_sqd, mixt_frac(k),
                                            big_m_w, small_m_w);
          adg_pdf_detail::close_Luhar_pdf(wm(k), wp2(k), mixt_frac(k), small_m_w, wp2(k), w_tol_sqd,
                                          sigma_sqd_w_1, sigma_sqd_w_2, varnce_w_1(k),
                                          varnce_w_2(k), w_1_n, w_2_n, w_1(k), w_2(k));
          // Solve for the thl PDF
          backsolve_Luhar_params(Skw(k), Skthl(k), big_m_w, mixt_frac(k), big_m_thl, small_m_thl);
          adg_pdf_detail::close_Luhar_pdf(thlm(k), thlp2(k), mixt_frac(k), small_m_thl, wpthlp(k),
                                          (thl_tol * thl_tol), sigma_sqd_thl_1, sigma_sqd_thl_2,
                                          varnce_thl_1(k), varnce_thl_2(k), thl_1_n, thl_2_n,
                                          thl_1(k), thl_2(k));
          // Solve for the rt PDF
          backsolve_Luhar_params(Skw(k), Skrt(k), big_m_w, mixt_frac(k), big_m_rt, small_m_rt);
          adg_pdf_detail::close_Luhar_pdf(rtm(k), rtp2(k), mixt_frac(k), small_m_rt, wprtp(k),
                                          (rt_tol * rt_tol), sigma_sqd_rt_1, sigma_sqd_rt_2,
                                          varnce_rt_1(k), varnce_rt_2(k), rt_1_n, rt_2_n, rt_1(k),
                                          rt_2(k));
        } else if ((Kokkos::abs(Skthl(k)) > Kokkos::abs(Skw(k))) &&
                   (Kokkos::abs(Skthl(k)) >= Kokkos::abs(Skrt(k)))) {
          // theta-l has the greatest magnitude of skewness.
          // Solve for the thl PDF
          adg_pdf_detail::calc_Luhar_params(Skthl(k), wpthlp(k), thlp2(k), (thl_tol * thl_tol),
                                            mixt_frac(k), big_m_thl, small_m_thl);
          // Solve for the thl PDF
          adg_pdf_detail::close_Luhar_pdf(thlm(k), thlp2(k), mixt_frac(k), small_m_thl, wpthlp(k),
                                          (thl_tol * thl_tol), sigma_sqd_thl_1, sigma_sqd_thl_2,
                                          varnce_thl_1(k), varnce_thl_2(k), thl_1_n, thl_2_n,
                                          thl_1(k), thl_2(k));
          // Solve for the w PDF
          backsolve_Luhar_params(Skthl(k), Skw(k), big_m_thl, mixt_frac(k), big_m_w, small_m_w);
          adg_pdf_detail::close_Luhar_pdf(wm(k), wp2(k), mixt_frac(k), small_m_w, wp2(k), w_tol_sqd,
                                          sigma_sqd_w_1, sigma_sqd_w_2, varnce_w_1(k),
                                          varnce_w_2(k), w_1_n, w_2_n, w_1(k), w_2(k));
          // Solve for the rt PDF
          backsolve_Luhar_params(Skthl(k), Skrt(k), big_m_thl, mixt_frac(k), big_m_rt, small_m_rt);
          adg_pdf_detail::close_Luhar_pdf(rtm(k), rtp2(k), mixt_frac(k), small_m_rt, wprtp(k),
                                          (rt_tol * rt_tol), sigma_sqd_rt_1, sigma_sqd_rt_2,
                                          varnce_rt_1(k), varnce_rt_2(k), rt_1_n, rt_2_n, rt_1(k),
                                          rt_2(k));
        } else {
          // rt has the greatest magnitude of skewness.
          // Solve for the rt PDF
          adg_pdf_detail::calc_Luhar_params(Skrt(k), wprtp(k), rtp2(k), (rt_tol * rt_tol),
                                            mixt_frac(k), big_m_rt, small_m_rt);
          // Solve for the rt PDF
          adg_pdf_detail::close_Luhar_pdf(rtm(k), rtp2(k), mixt_frac(k), small_m_rt, wprtp(k),
                                          (rt_tol * rt_tol), sigma_sqd_rt_1, sigma_sqd_rt_2,
                                          varnce_rt_1(k), varnce_rt_2(k), rt_1_n, rt_2_n, rt_1(k),
                                          rt_2(k));
          // Solve for the w PDF
          backsolve_Luhar_params(Skrt(k), Skw(k), big_m_rt, mixt_frac(k), big_m_w, small_m_w);
          adg_pdf_detail::close_Luhar_pdf(wm(k), wp2(k), mixt_frac(k), small_m_w, wp2(k), w_tol_sqd,
                                          sigma_sqd_w_1, sigma_sqd_w_2, varnce_w_1(k),
                                          varnce_w_2(k), w_1_n, w_2_n, w_1(k), w_2(k));
          // Solve for the thl PDF
          backsolve_Luhar_params(Skrt(k), Skthl(k), big_m_rt, mixt_frac(k), big_m_thl, small_m_thl);
          adg_pdf_detail::close_Luhar_pdf(thlm(k), thlp2(k), mixt_frac(k), small_m_thl, wpthlp(k),
                                          (thl_tol * thl_tol), sigma_sqd_thl_1, sigma_sqd_thl_2,
                                          varnce_thl_1(k), varnce_thl_2(k), thl_1_n, thl_2_n,
                                          thl_1(k), thl_2(k));
        }
      });
}

void ADG1_w_closure(int nz, int ngrdcol, Real2D wm, Real2D wp2, Real2D Skw,
                    Real2D sigma_sqd_w,                                 // In
                    Real2D sqrt_wp2, Real1D mixt_frac_max_mag,          // In
                    err_info_type& err_info,                            // InOut
                    Real2D w_1, Real2D w_2, Real2D w_1_n, Real2D w_2_n, // Out
                    Real2D varnce_w_1, Real2D varnce_w_2, Real2D mixt_frac) {
  // Description:
  // Calculates the mixture fraction, the PDF component means of w, and the PDF
  // component variances of w for the Analytic Double Gaussian 1 (ADG1)
  // closure.  It sets the widths of both w Gaussians to be the same, and
  // furthermore, relates them to the overall variance of w (<w'^2>) by a
  // parameter, sigma_sqd_w.  The equation is:
  //
  // sigma_w_1^2 = sigma_w_2^2 = sigma_sqd_w * <w'^2>;
  //
  // where sigma_w_1^2 is the variance of w in the 1st PDF component,
  // sigma_w_2^2 is the variance of w in the 2nd PDF component, and parameter
  // sigma_sqd_w must have a value within the range 0 <= sigma_sqd_w <= 1.
  //
  // References:
  // Golaz, J-C., V. E. Larson, and W. R. Cotton, 2002a: A PDF-based model for
  // boundary layer clouds. Part I: Method and model description. J. Atmos.
  // Sci., 59, 3540–3551.
  //
  // Vincent E. Larson and Jean-Christophe Golaz, 2005: Using Probability
  // Density Functions to Derive Consistent Closure Relationships among
  // Higher-Order Moments. Mon. Wea. Rev., 133, 1023–1042.
  //-----------------------------------------------------------------------

  if (clubb_at_least_debug_level_api(2)) {
    int below_tolerance = 0;
    Kokkos::parallel_reduce(
        "ADG1_w_closure::variance_check", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell, int& count) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          if (wp2(i, k) < w_tol_sqd)
            ++count;
        },
        below_tolerance);
    if (below_tolerance) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "Error in ADG1_w_closure: wp2 must not be less than the w tolerance value.\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }
  Kokkos::parallel_for(
      "ADG1_w_closure::component_parameters", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        // Width (standard deviation) parameters are non-zero
        // The variable "mixt_frac" is the weight of the 1st PDF component.  The
        // weight of the 2nd PDF component is "1-mixt_frac".  If there isn't any
        // skewness of w (Sk_w = 0 because w'^3 = 0), mixt_frac = 0.5, and both
        // PDF components are equally weighted.  If there is positive skewness of
        // w (Sk_w > 0 because w'^3 > 0), 0 < mixt_frac < 0.5, and the 2nd PDF
        // component has greater weight than does the 1st PDF component.  If there
        // is negative skewness of w (Sk_w < 0 because w'^3 < 0),
        // 0.5 < mixt_frac < 1, and the 1st PDF component has greater weight than
        // does the 2nd PDF component.
        if (Kokkos::abs(Skw(i, k)) <= real_t(1.0e-5)) {
          mixt_frac(i, k) = one_half;
        } else {
          mixt_frac(i, k) =
              one_half * (one - Skw(i, k) / Kokkos::sqrt(four * ((one - sigma_sqd_w(i, k)) *
                                                                 (one - sigma_sqd_w(i, k)) *
                                                                 (one - sigma_sqd_w(i, k))) +
                                                         (Skw(i, k) * Skw(i, k))));
        }
        // Clip mixt_frac, and 1 - mixt_frac, to avoid dividing by a small number.
        // Formula for mixt_frac_max_mag =
        // 1 - ( 1/2 * ( 1 - Skw_max
        //                   / sqrt( 4*( 1 - sigma_sqd_w )^3 + Skw_max^2 ) ) ),
        // where sigma_sqd_w is fixed at 0.4 for this calculation.
        mixt_frac(i, k) = ((((mixt_frac(i, k)) > (one - mixt_frac_max_mag(i))
                                 ? (mixt_frac(i, k))
                                 : (one - mixt_frac_max_mag(i)))) < (mixt_frac_max_mag(i))
                               ? (((mixt_frac(i, k)) > (one - mixt_frac_max_mag(i))
                                       ? (mixt_frac(i, k))
                                       : (one - mixt_frac_max_mag(i))))
                               : (mixt_frac_max_mag(i)));
        // The normalized mean of w for Gaussian "plume" 1 is w_1_n.  It's value
        // will always be greater than 0.  As an example, a value of 1.0 would
        // indicate that the actual mean of w for Gaussian "plume" 1 is found 1.0
        // standard deviation above the overall mean for w.
        w_1_n(i, k) =
            Kokkos::sqrt(((one - mixt_frac(i, k)) / mixt_frac(i, k)) * (one - sigma_sqd_w(i, k)));
        // The normalized mean of w for Gaussian "plume" 2 is w_2_n.  It's value
        // will always be less than 0.  As an example, a value of -0.5 would
        // indicate that the actual mean of w for Gaussian "plume" 2 is found 0.5
        // standard deviations below the overall mean for w.
        w_2_n(i, k) =
            -Kokkos::sqrt((mixt_frac(i, k) / (one - mixt_frac(i, k))) * (one - sigma_sqd_w(i, k)));
        // The mean of w for Gaussian "plume" 1 is w_1.
        w_1(i, k) = wm(i, k) + sqrt_wp2(i, k) * w_1_n(i, k);
        // The mean of w for Gaussian "plume" 2 is w_2.
        w_2(i, k) = wm(i, k) + sqrt_wp2(i, k) * w_2_n(i, k);
        // The variance of w for Gaussian "plume" 1 for varnce_w_1.
        varnce_w_1(i, k) = sigma_sqd_w(i, k) * wp2(i, k);
        // The variance of w for Gaussian "plume" 2 for varnce_w_2.
        // The variance in both Gaussian "plumes" is defined to be the same.
        varnce_w_2(i, k) = sigma_sqd_w(i, k) * wp2(i, k);
      });
}

void calc_Luhar_params(int nz, Real1D Skx, Real1D wpxp, // In
                       Real1D xp2, real_t x_tol_sqd,    // In
                       Real1D mixt_frac, Real1D big_m, Real1D small_m) {
  // Description:
  // For the Luhar closure, this subroutine takes Skx (and Skw) as input and
  // outputs the mixture fraction, big_m, and small_m. This code was written
  // using the equations and nomenclature of Larson et al. (2002) Appendix
  // section e.
  //
  // The relationship between skewness of x (Skx), mixture fraction (a), and
  // Luhar's small m (m) is given by:
  //
  // Skx^2 = ( m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3 )
  //         * ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).
  //
  // Luhar's large M (M) is used to more easily express the factor involving
  // the m's:
  //
  // M = ( m^2 + 1 )^3 / ( m^2 * ( m^2 + 3 )^2 ).
  //
  // The equation involving skewness of x becomes:
  //
  // Skx^2 = ( 1 / M ) * ( 1 - 2*a )^2 / ( a * ( 1 - a ) );
  //
  // or:
  //
  // M * Skx^2 = ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).
  //
  // This equation can be rewritten as:
  //
  // ( a * ( 1 - a ) ) * M * Skx^2 = ( 1 - 2*a )^2;
  //
  // as well as:
  //
  // ( a - a^2 ) * M * Skx^2 = 1 - 4*a + 4*a^2;
  //
  // and eventually as:
  //
  // ( 4 + M * Skx^2 ) * a^2 - ( 4 + M * Skx^2 ) * a + 1 = 0.
  //
  // Solving the quadratic equation for a:
  //
  // a = (1/2) * ( 1 +- Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).
  //
  // Since by definition, mu_w_1 >= mu_w_2, a < 0.5 when Skw > 0, the equation
  // for mixture fraction is:
  //
  // a = (1/2) * ( 1 - Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).
  //
  // For 3-D Luhar, the variable (w, rt, or theta-l) with the greatest
  // magnitude of skewness is used to calculate mixture fraction.  Since it is
  // desirable to still have a < 0.5 when Skw > 0 and a > 0.5 when Skw < 0, the
  // sign function is used.  The value of Skx is replaced by:
  //
  // Skx|_adj = sgn( <w'x'> ) * Skx;
  //
  // where
  //
  // sgn( <w'x'> ) = | 1 when <w'x'> >= 0
  //                 | -1 when <w'x'> < 0.
  //
  // Since Skx|_adj^2 = ( sgn( <w'x'> ) * Skx )^2 = ( sgn( <w'x'> ) )^2 * Skx^2
  // = Skx^2, the equation for mixture fraction is:
  //
  // a = (1/2) * ( 1 - sgn( <w'x'> ) * Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).
  //
  // When using the ADG2 closure or when using the 3-D Luhar closure when the
  // variable with the greatest magnitude of skewness is w, Skw = Skx and
  // sgn( <w'^2> ) is always equal to 1, reducing the equation to its previous
  // form.

  // References:
  // Vincent E. Larson, Jean-Christophe Golaz, and William R. Cotton, 2002:
  // Small-Scale and Mesoscale Variability in Cloudy Boundary Layers: Joint
  // Probability Density Functions. J. Atmos. Sci., 59, 3519–3539.
  //-----------------------------------------------------------------------

  Kokkos::parallel_for(
      "calc_Luhar_params", nz, KOKKOS_LAMBDA(const int k) {
        adg_pdf_detail::calc_Luhar_params(Skx(k), wpxp(k), xp2(k), x_tol_sqd, mixt_frac(k),
                                          big_m(k), small_m(k));
      });
}

namespace adg_pdf_detail {
KOKKOS_INLINE_FUNCTION void calc_Luhar_params(real_t Skx, real_t wpxp,      // In
                                              real_t xp2, real_t x_tol_sqd, // In
                                              real_t& mixt_frac, real_t& big_m, real_t& small_m) {
  // Device scalar form of the Native independent-level equations.
  real_t small_m_sqd,       // Square of Luhar's small m [-]
      sgn_wpxp;             // Sign of covariance wpxp [-]
  real_t small_m_candidate; // Evaluate the native MAX operand once.
  if (xp2 > x_tol_sqd) {
    // Width (standard deviation) parameters are non-zero
    // Calculate Luhar's m (small m).
    // If Skx is very small, then small_m will tend to zero which risks
    // divide-by-zero.  To ameliorate this problem, we enforce abs( x_1_n )
    // and abs( x_2_n ) > 0.05.
    // Note:  Luhar's small_m (m) is the only tunable parameter in the Luhar
    //        closure, so this equation can be changed.  However, the value
    //        of m should go toward 0 as Skx goes toward 0 so that the double
    //        Gaussian reduces to a single Gaussian when the distribution is
    //        unskewed.
    small_m_candidate = two_thirds * Kokkos::pow(Kokkos::abs(Skx), one_third);
    small_m = ((small_m_candidate) > (real_t(0.05)) ? (small_m_candidate) : (real_t(0.05)));
    // Calculate m^2.
    small_m_sqd = (small_m * small_m);
    // Calculate Luhar's M (big M).
    big_m = ((one + small_m_sqd) * (one + small_m_sqd) * (one + small_m_sqd)) /
            (((three + small_m_sqd) * (three + small_m_sqd)) * small_m_sqd);
    // Calculate sgn( <w'x'> ).
    if (wpxp >= zero) {
      sgn_wpxp = one;
    } else {
      sgn_wpxp = -one;
    }
    // Calculate mixture fraction.
    mixt_frac =
        one_half * (one - sgn_wpxp * Skx * Kokkos::sqrt(one / ((four / big_m) + (Skx * Skx))));
  } else {
    // Variable x doesn't vary.
    mixt_frac = one_half;
    // For output purposes, set small_m and big_m to 0.
    small_m = zero;
    big_m   = zero;
  }
}
} // namespace adg_pdf_detail

void close_Luhar_pdf(int nz, Real1D xm, Real1D xp2, Real1D mixt_frac, // In
                     Real1D small_m, Real1D wpxp, real_t x_tol_sqd,   // In
                     Real1D sigma_sqd_x_1, Real1D sigma_sqd_x_2,      // Out
                     Real1D varnce_x_1, Real1D varnce_x_2,            // Out
                     Real1D x_1_n, Real1D x_2_n, Real1D x_1, Real1D x_2) {
  // Description:
  // For the Luhar closure, this subroutine takes Skx, xm, xp2, and mixt_frac,
  // big_m, and small_m (calculated in calc_Luhar_params) as input and outputs
  // the PDF component means and variances of a variable x in the joint-PDF
  // according to Luhar et al. (1996).  This code was written using the
  // equations and nomenclature of Larson et al. (2002) Appendix section e.

  // References:
  //    Vincent E. Larson, Jean-Christophe Golaz, and William R. Cotton, 2002:
  //    Small-Scale and Mesoscale Variability in Cloudy Boundary Layers: Joint
  //    Probability Density Functions. J. Atmos. Sci., 59, 3519–3539.
  //-----------------------------------------------------------------------

  Kokkos::parallel_for(
      "close_Luhar_pdf", nz, KOKKOS_LAMBDA(const int k) {
        adg_pdf_detail::close_Luhar_pdf(xm(k), xp2(k), mixt_frac(k), small_m(k), wpxp(k), x_tol_sqd,
                                        sigma_sqd_x_1(k), sigma_sqd_x_2(k), varnce_x_1(k),
                                        varnce_x_2(k), x_1_n(k), x_2_n(k), x_1(k), x_2(k));
      });
}

namespace adg_pdf_detail {
KOKKOS_INLINE_FUNCTION void close_Luhar_pdf(real_t xm, real_t xp2, real_t mixt_frac,       // In
                                            real_t small_m, real_t wpxp, real_t x_tol_sqd, // In
                                            real_t& sigma_sqd_x_1, real_t& sigma_sqd_x_2,  // Out
                                            real_t& varnce_x_1, real_t& varnce_x_2,        // Out
                                            real_t& x_1_n, real_t& x_2_n, real_t& x_1,
                                            real_t& x_2) {
  // Device scalar form of the Native independent-level equations.
  real_t sqrt_xp2, // Overall standard deviation [units of x]
      sgn_wpxp;    // Sign of covariance wpxp [-]
  // Definedness deviation: Native initializes sign only for xp2>tolerance
  // then reads it below tolerance too. Match JAX's defined covariance sign;
  // physical means/variances retain Native formulas. No Native BFB claim
  // applies to its undefined nonvarying normalized means (Native1022–1060).
  sgn_wpxp = wpxp >= zero ? one : -one;
  if (xp2 > x_tol_sqd) {
    // Width (standard deviation) parameters are non-zero
    // Calculate sgn( <w'x'> ).
    // Sign is initialized for both branches immediately above.
    // Calculate the square root of the overall variance of x.
    sqrt_xp2 = Kokkos::sqrt(xp2);
    // Normalized width parameter of x in the 1st PDF component.
    sigma_sqd_x_1 = (one - mixt_frac) / (mixt_frac * (one + (small_m * small_m)));
    // The variance of x in the 1st PDF component.
    varnce_x_1 = sigma_sqd_x_1 * xp2;
    // Normalized width parameter of x in the 2nd PDF component.
    sigma_sqd_x_2 = mixt_frac / ((one - mixt_frac) * (one + (small_m * small_m)));
    // The variance of x in the 2nd PDF component.
    varnce_x_2 = sigma_sqd_x_2 * xp2;
    // Normalized mean of x in the 1st PDF component.
    x_1_n = sgn_wpxp * small_m * Kokkos::sqrt(sigma_sqd_x_1);
    // Normalized mean of x in the 2nd PDF component.
    x_2_n = -sgn_wpxp * small_m * Kokkos::sqrt(sigma_sqd_x_2);
    // The mean of x in the 1st PDF component.
    x_1 = xm + sqrt_xp2 * x_1_n;
    // The mean of x in the 2nd PDF component.
    x_2 = xm + sqrt_xp2 * x_2_n;
  } else {
    // Variable x doesn't vary.
    sigma_sqd_x_1 = (one / (one + (small_m * small_m)));
    sigma_sqd_x_2 = (one / (one + (small_m * small_m)));
    varnce_x_1    = zero;
    varnce_x_2    = zero;
    x_1_n         = sgn_wpxp * small_m * Kokkos::sqrt(sigma_sqd_x_1);
    x_2_n         = -sgn_wpxp * small_m * Kokkos::sqrt(sigma_sqd_x_2);
    x_1           = xm;
    x_2           = xm;
  }
}
} // namespace adg_pdf_detail

namespace {

void ADG1_ADG2_responder_params(int nz, int ngrdcol,                                       // In
                                Real2D xm, Real2D xp2, Real2D wp2, Real2D sqrt_wp2,        // In
                                Real2D wpxp, Real2D w_1_n, Real2D w_2_n, Real2D mixt_frac, // In
                                Real2D sigma_sqd_w, Real1D beta,                           // In
                                real_t x_tol,                                              // In
                                err_info_type& err_info,                                   // InOut
                                Real2D x_1, Real2D x_2, Real2D varnce_x_1,                 // Out
                                Real2D varnce_x_2, Real2D alpha_x) {
  // Description:
  // Calculates the PDF component means and PDF component variances for thl,
  // rt, or sclr when either the ADG1 PDF or ADG2 PDF are used.
  //
  // The normalized variance for thl, rt, and sclr for "plume" 1 is:
  //
  // { 1 - [1/(1-sigma_sqd_w)]*[ <w'x'>^2 / (<w'^2> * <x'^2>) ] / mixt_frac }
  // * { (1/3)*beta + mixt_frac*( 1 - (2/3)*beta ) };
  //
  // where "x" stands for thl, rt, or sclr; "mixt_frac" is the weight of
  // Gaussian "plume" 1, and 0 <= beta <= 3.
  //
  // The factor { (1/3)*beta + mixt_frac*( 1 - (2/3)*beta ) } does not depend
  // on which varable "x" stands for.  The factor is multiplied by 2 and
  // defined as width_factor_1.
  //
  // The factor
  // { 1 - [1/(1-sigma_sqd_w)]*[ (w'x')^2 / (w'^2 * x'^2) ] / mixt_frac }
  // depends on which variable "x" stands for.  It is multiplied by 1/2 and
  // defined as alpha_x, where "x" stands for thl, rt, or sclr.

  // Vince Larson added a dimensionless factor so that the
  // width of plumes in theta_l, rt can vary.
  // beta is a constant defined in CLUBB's tunable parameters.
  //   Set 0<beta<3.
  // beta=1.5_core_rknd recovers Chris Golaz' simplified formula.
  // 3 Nov 2003

  // References:

  //-------------------------- Begin Code --------------------------
  if (clubb_at_least_debug_level_api(2)) {
    // Integer ANY reduction leaves valid data on the execution space.
    int below_tolerance = 0;
    Kokkos::parallel_reduce(
        "ADG1_ADG2_responder_params::variance_check", ngrdcol * nz,
        KOKKOS_LAMBDA(const int cell, int& count) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          if (xp2(i, k) < x_tol * x_tol)
            ++count;
        },
        below_tolerance);
    if (below_tolerance) {
      std::cerr << err_info.err_header_global << '\n';
      std::cerr << "Error in ADG1_ADG2_responder_params: xp2 must not be less "
                   "than its tolerance value.\n";
      Kokkos::deep_copy(err_info.err_code, clubb_fatal_error);
      return;
    }
  }

  Kokkos::parallel_for(
      "ADG1_ADG2_responder_params::component_parameters", ngrdcol * nz,
      KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        real_t width_factor_1; // Width factor relating to PDF component1 [-]

        // This if-statement has been commented out in favor of the error check above.
        // These fields should never be below tolerance, and these tolerances are too large
        // to make this if-case a proper asymptote of the computation in the else-case.
        // See https://github.com/larson-group/clubb/pull/1276 for more detail.
        // if (wp2(i,k) < w_tol_sqd || xp2(i,k) < x_tol*x_tol) {
        //   // Vertical velocity doesn't vary. Variable x is treated as a single
        //   // Gaussian in this situation.
        //   x_1(i,k) = xm(i,k);
        //   x_2(i,k) = xm(i,k);
        //   varnce_x_1(i,k) = xp2(i,k);
        //   varnce_x_2(i,k) = xp2(i,k);
        //   alpha_x(i,k) = one_half;
        // } else {

        x_1(i, k) = xm(i, k) - wpxp(i, k) / (sqrt_wp2(i, k) * w_2_n(i, k));
        x_2(i, k) = xm(i, k) - wpxp(i, k) / (sqrt_wp2(i, k) * w_1_n(i, k));

        alpha_x(i, k) = one_half * (one - wpxp(i, k) * wpxp(i, k) /
                                              ((one - sigma_sqd_w(i, k)) * wp2(i, k) * xp2(i, k)));

        // Native MIN/MAX select their second argument on an equality tie.
        alpha_x(i, k) = alpha_x(i, k) < one ? alpha_x(i, k) : one;
        alpha_x(i, k) = alpha_x(i, k) > zero_threshold ? alpha_x(i, k) : zero_threshold;

        // Native rounds the second term, then fuses two_thirds*beta into the sum.
        width_factor_1 =
            Kokkos::fma(two_thirds, beta(i), two * mixt_frac(i, k) * (one - two_thirds * beta(i)));

        // Vince Larson multiplied original expressions by width_factor_1,2
        // to generalize scalar skewnesses. 05Nov2003.
        varnce_x_1(i, k) = width_factor_1 * xp2(i, k) * alpha_x(i, k) / mixt_frac(i, k);
        varnce_x_2(i, k) =
            (two - width_factor_1) * xp2(i, k) * alpha_x(i, k) / (one - mixt_frac(i, k));

        // } // Former single-Gaussian branch remains inactive.
      });
}

KOKKOS_INLINE_FUNCTION void backsolve_Luhar_params(real_t Sk_max, real_t Skx,          // In
                                                   real_t big_m_max, real_t mixt_frac, // In
                                                   real_t& big_m_x, real_t& small_m_x) {
  // Description:
  // This subroutine calculates Luhar's big_m and small_m for the variate 'x'
  // consistent with the mixture fraction of the variate with the largest
  // skewness.
  //
  // The relationship between skewness of x (Skx), mixture fraction (a), and
  // Luhar's small m (m) is given by:
  //
  // Skx^2 = ( m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3 )
  //         * ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).
  //
  // Moving the factor involving mixture fraction to the right-hand side:
  //
  // ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2
  // = m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3.
  //
  // This can be rewritten as:
  //
  // ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2
  // = ( m^6 + 6*m^4 + 9*m^2 ) / ( m^6 + 3*m^4 + 3*m^2 + 1 ).
  //
  // Setting alpha = ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2, the
  // equation can be rewritten as:
  //
  // ( m^6 + 3*m^4 + 3*m^2 + 1 ) * alpha = m^6 + 6*m^4 + 9*m^2.
  //
  // This can be rearranged and rewritten as:
  //
  // ( alpha - 1 ) * m^6 + ( 3 * alpha - 6 ) * m^4
  // + ( 3 * alpha - 9 ) * m^2 + alpha = 0.
  //
  // This can be rewritten again as:
  //
  // ( alpha - 1 ) * (m^2)^3 + ( 3 * alpha - 6 ) * (m^2)^2
  // + ( 3 * alpha - 9 ) * (m^2) + alpha = 0.
  //
  // The goal is to solve for m^2, and then take the square root of m^2 to
  // solve for m.  This can be accomplished by using the cubic formula (with
  // the l_use_cubic_backsolve option), or else by a quadratic approximation.

  // References:
  //-----------------------------------------------------------------------

  real_t alpha, a, b, c, alpha_upr, alpha_low, discrim;
  real_t small_m_root; // Own the native function-result temporary: call once.
  Kokkos::Array<real_t, 1> a_coef_in;
  Kokkos::Array<real_t, 1> b_coef_in;
  Kokkos::Array<real_t, 1> c_coef_in;
  Kokkos::Array<real_t, 1> d_coef_in;
  constexpr bool l_use_cubic_backsolve = true;
  if (l_use_cubic_backsolve) {
    if (Kokkos::abs(mixt_frac - one_half) < real_t(0.001)) {
      // When mixture fraction = 0.5 (based on the variable with the largest
      // magnitude of skewness), all variables must have a skewness of 0.
      // Set m to the minimum threshold of 0.05.
      small_m_x = real_t(0.05);
      // Calculate the corresponding value of big_m_x.
      big_m_x = ((one + (small_m_x * small_m_x)) * (one + (small_m_x * small_m_x)) *
                 (one + (small_m_x * small_m_x))) /
                (((three + (small_m_x * small_m_x)) * (three + (small_m_x * small_m_x))) *
                 (small_m_x * small_m_x));
    } else if (Kokkos::abs(Skx) < eps) {
      // Mixture fraction /= 0.5 because the variable with the largest
      // magnitude of skewness has a skewness /= 0.  However, variable x has
      // a skewness of 0.  In order to reproduce the correct skewness for
      // variable x, set m to 0 (regardless of minimum thresholds used in
      // other parts of the code).
      small_m_x = zero;
      // The value of big_m_x should be inf.  Set it to huge.  This is not
      // used in any calculation, anyway.
      // The Kokkos trait supplies HUGE without a host-only function call.
      big_m_x = Kokkos::finite_max_v<real_t>;
    } else {
      // Backsolve for m, given mixt_frac and Skx.
      // alpha = 1/M is given by:
      // [ mixt_frac * ( 1 - mixt_frac ) / ( 1 - 2 * mixt_frac )^2 ] * Skx^2.
      alpha =
          (mixt_frac * (one - mixt_frac) / ((one - two * mixt_frac) * (one - two * mixt_frac))) *
          (Skx * Skx);
      // Calculate big_m_x.
      big_m_x = one / alpha;
      // Solve the cubic equation for m^2:
      // ( alpha - 1 ) * (m^2)^3 + ( 3 * alpha - 6 ) * (m^2)^2
      // + ( 3 * alpha - 9 ) * (m^2) + alpha = 0.
      // The largest root is preferred.
      a_coef_in[0] = alpha - one;
      b_coef_in[0] = three * alpha - six;
      c_coef_in[0] = three * alpha - nine;
      d_coef_in[0] = alpha;
      small_m_root = max_cubic_root(a_coef_in, b_coef_in, c_coef_in, d_coef_in);
      small_m_x    = Kokkos::sqrt(((small_m_root) > ((real_t(0.05) * real_t(0.05)))
                                       ? (small_m_root)
                                       : ((real_t(0.05) * real_t(0.05)))));
    }
  } else {
    alpha =
        ((Skx * Skx) / ((((Sk_max * Sk_max)) > (eps) ? ((Sk_max * Sk_max)) : (eps)) * big_m_max));
    // This limit keeps the discriminant >= 0
    alpha_upr = two * Kokkos::sqrt(real_t(13.0)) - five;
    alpha_low = eps;
    // For this approximation, alpha must be less than 2*sqrt(13) - 5 to get a
    // real ans.
    alpha = ((alpha) < (alpha_upr) ? (alpha) : (alpha_upr));
    // For testing, eliminate possibility of divide by zero
    alpha = ((alpha) > (alpha_low) ? (alpha) : (alpha_low));
    // Use a piece-wise approximation
    if (alpha < one) {
      a         = ((three * alpha - six) > (eps) ? (three * alpha - six) : (eps));
      b         = three * alpha - nine;
      c         = alpha;
      discrim   = (b * b) - four * a * c;
      small_m_x = Kokkos::sqrt((-b - Kokkos::sqrt(discrim)) / (two * a));
    } else {
      // For this approximation, alpha must be less than 2*sqrt(13) - 5
      // to get a real ans.
      alpha     = ((alpha) < (two) ? (alpha) : (two));
      a         = ((six * alpha - nine) > (eps) ? (six * alpha - nine) : (eps));
      b         = -six;
      c         = two * alpha - one;
      discrim   = (b * b) - four * a * c;
      small_m_x = Kokkos::sqrt((-b - Kokkos::sqrt(discrim)) / (two * a));
    }
    // Clip consistently with subroutine calc_Luhar_params
    small_m_x = ((real_t(5e-2)) > (small_m_x) ? (real_t(5e-2)) : (small_m_x));
    big_m_x   = one / alpha;
  }
}

KOKKOS_INLINE_FUNCTION real_t max_cubic_root(const Kokkos::Array<real_t, 1>& a_coef,
                                             const Kokkos::Array<real_t, 1>& b_coef,
                                             const Kokkos::Array<real_t, 1>& c_coef,
                                             const Kokkos::Array<real_t, 1>& d_coef) {
  // Description:
  // Calculates the largest root that results from solving a cubic equation of
  // the form a*x^3 + b*x^2 + c*x + d = 0.
  //
  // This is done to backsolve for m^2 for the 3-D Luhar closure, given the
  // values of mixt_frac and Skx.

  // References:
  //-----------------------------------------------------------------------

  Kokkos::Array<complex_t, 3> cubic_roots;
  Kokkos::Array<complex_t, 2> quadratic_roots;
  real_t max_root, a_coef_thresh, b_coef_thresh;
  // Calculate a minimum threshold for |a| to call this a cubic equation.
  a_coef_thresh =
      real_t(0.001) *
      ((((Kokkos::abs(b_coef[0])) > (Kokkos::abs(c_coef[0])) ? (Kokkos::abs(b_coef[0]))
                                                             : (Kokkos::abs(c_coef[0])))) >
               (Kokkos::abs(d_coef[0]))
           ? (((Kokkos::abs(b_coef[0])) > (Kokkos::abs(c_coef[0])) ? (Kokkos::abs(b_coef[0]))
                                                                   : (Kokkos::abs(c_coef[0]))))
           : (Kokkos::abs(d_coef[0])));
  // Calculate a minimum threshold for |b| to call this a quadratic equation.
  // This only matters when |a| <= a_coef_thresh.
  b_coef_thresh = real_t(0.001) * ((Kokkos::abs(c_coef[0])) > (Kokkos::abs(d_coef[0]))
                                       ? (Kokkos::abs(c_coef[0]))
                                       : (Kokkos::abs(d_coef[0])));
  if (Kokkos::abs(a_coef[0]) > a_coef_thresh) {
    // The equation is a cubic equation.
    cubic_roots = cubic_solve(a_coef[0], b_coef[0], c_coef[0], d_coef[0]);
    if (Kokkos::abs(cubic_roots[1].imag()) < eps && Kokkos::abs(cubic_roots[2].imag()) < eps) {
      // Find the maximum root of the three roots.
      max_root =
          ((((cubic_roots[0].real()) > (cubic_roots[1].real()) ? (cubic_roots[0].real())
                                                               : (cubic_roots[1].real()))) >
                   (cubic_roots[2].real())
               ? (((cubic_roots[0].real()) > (cubic_roots[1].real()) ? (cubic_roots[0].real())
                                                                     : (cubic_roots[1].real())))
               : (cubic_roots[2].real()));
    } else {
      max_root = cubic_roots[0].real();
    }
  } else if (Kokkos::abs(b_coef[0]) > b_coef_thresh) {
    // The equation is a quadratic equation, since a = 0, but b /= 0.
    // This should very rarely occur for 3-D Luhar.  When it does, the result
    // will always be two real-valued roots.
    quadratic_roots = quadratic_solve(b_coef[0], c_coef[0], d_coef[0]);
    // Find the maximum root of the two roots.
    max_root =
        ((quadratic_roots[0].real()) > (quadratic_roots[1].real()) ? (quadratic_roots[0].real())
                                                                   : (quadratic_roots[1].real()));
  } else {
    // The equation is a linear equation.
    // This won't happen for 3-D Luhar.
    max_root = -d_coef[0] / c_coef[0];
  }
  return max_root;
}

} // namespace
} // namespace clubb_kokkos
