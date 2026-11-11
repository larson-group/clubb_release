// ADG1 closure orchestration and responder parameters; see the public header.
#include "adg1_adg2_3d_luhar_pdf.hpp"
#include "constants_clubb.hpp"
#include "error_code.hpp"
#include "f2k/procedure_bridge.hpp"
#include <iostream>

namespace clubb_kokkos {

namespace {

void ADG1_ADG2_responder_params(int nz, int ngrdcol,                                       // In
                                Real2D xm, Real2D xp2, Real2D wp2, Real2D sqrt_wp2,        // In
                                Real2D wpxp, Real2D w_1_n, Real2D w_2_n, Real2D mixt_frac, // In
                                Real2D sigma_sqd_w, Real1D beta,                           // In
                                real_t x_tol,                                              // In
                                err_info_type& err_info,                                   // InOut
                                Real2D x_1, Real2D x_2, Real2D varnce_x_1,                 // Out
                                Real2D varnce_x_2, Real2D alpha_x);

} // namespace

void ADG1_pdf_driver(int nz, int ngrdcol, int sclr_dim, Real1D sclr_tol,            // In
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
  Real2D w_1_n("ADG1_pdf_driver::w_1_n", ngrdcol, nz); // Normalized component1 w mean [-]
  Real2D w_2_n("ADG1_pdf_driver::w_2_n", ngrdcol, nz); // Normalized component2 w mean [-]

  // Calculate mixture fraction and PDF component means and variances of w.
  // TODO: Replace this native lower call when ADG1_w_closure is ported.
  f2k::ADG1_w_closure(nz, ngrdcol, wm, wp2, Skw, sigma_sqd_w, // In
                      sqrt_wp2, mixt_frac_max_mag,            // In
                      err_info,                               // InOut
                      w_1, w_2, w_1_n, w_2_n,                 // Out
                      varnce_w_1, varnce_w_2, mixt_frac);     // Out

  if (f2k::clubb_at_least_debug_level_api(0)) {
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
    const auto sclr_tol_host = to_host_space(sclr_tol);
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
  if (f2k::clubb_at_least_debug_level_api(2)) {
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

        width_factor_1 =
            two_thirds * beta(i) + two * mixt_frac(i, k) * (one - two_thirds * beta(i));

        // Vince Larson multiplied original expressions by width_factor_1,2
        // to generalize scalar skewnesses. 05Nov2003.
        varnce_x_1(i, k) = width_factor_1 * xp2(i, k) * alpha_x(i, k) / mixt_frac(i, k);
        varnce_x_2(i, k) =
            (two - width_factor_1) * xp2(i, k) * alpha_x(i, k) / (one - mixt_frac(i, k));

        // } // Former single-Gaussian branch remains inactive.
      });
}

} // namespace

} // namespace clubb_kokkos
