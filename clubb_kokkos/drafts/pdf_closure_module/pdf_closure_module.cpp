// Higher moments integrated over the two-component normal PDF.
#include "pdf_closure_module.hpp"
#include "constants_clubb.hpp"

namespace clubb_kokkos {

void calc_wp4_pdf(int nz, int ngrdcol,                  // In
                  Real2D wm, Real2D w_1, Real2D w_2,    // In
                  Real2D varnce_w_1, Real2D varnce_w_2, // In
                  Real2D mixt_frac,                     // In
                  Real2D wp4) {
  // Description:
  // Calculates <w'^4> by integrating over the PDF of w.  The integral is:
  //
  // <w'^4> = INT(-inf:inf) ( w - <w> )^4 P(w) dw;
  //
  // where <w> is the overall mean of w and P(w) is a two-component normal
  // distribution of w.  The integrated equation is:
  //
  // <w'^4> = mixt_frac * ( 3 * sigma_w_1^4
  //                        + 6 * ( mu_w_1 - <w> )^2 * sigma_w_1^2
  //                        + ( mu_w_1 - <w> )^4 )
  //          + ( 1 - mixt_frac ) * ( 3 * sigma_w_2^4
  //                                  + 6 * ( mu_w_2 - <w> )^2 * sigma_w_2^2
  //                                  + ( mu_w_2 - <w> )^4 );
  //
  // where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
  // of w in the 2nd PDF component, sigma_w_1 is the standard deviation of w in
  // the 1st PDF component, sigma_w_2 is the standard deviation of w in the 2nd
  // PDF component, and mixt_frac is the mixture fraction, which is the weight
  // of the 1st PDF component.

  // References:
  //-----------------------------------------------------------------------

  // Arrays are indexed (column, level).
  // wm: Mean of w (overall) [m/s]
  // w_1: Mean of w (1st PDF component) [m/s]
  // w_2: Mean of w (2nd PDF component) [m/s]
  // varnce_w_1: Variance of w (1st PDF component) [m^2/s^2]
  // varnce_w_2: Variance of w (2nd PDF component) [m^2/s^2]
  // mixt_frac: Mixture fraction [-]
  // wp4: <w'^4> [m^4/s^4]

  // Calculate <w'^4> by integrating over the PDF.
  Kokkos::parallel_for(
      "calc_wp4_pdf", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        wp4(i, k) =
            mixt_frac(i, k) *
                (three * (varnce_w_1(i, k) * varnce_w_1(i, k)) +
                 six * (((w_1(i, k) - wm(i, k)) * (w_1(i, k) - wm(i, k)))) * varnce_w_1(i, k) +
                 (((w_1(i, k) - wm(i, k)) * (w_1(i, k) - wm(i, k))) *
                  ((w_1(i, k) - wm(i, k)) * (w_1(i, k) - wm(i, k))))) +
            (one - mixt_frac(i, k)) *
                (three * (varnce_w_2(i, k) * varnce_w_2(i, k)) +
                 six * (((w_2(i, k) - wm(i, k)) * (w_2(i, k) - wm(i, k)))) * varnce_w_2(i, k) +
                 (((w_2(i, k) - wm(i, k)) * (w_2(i, k) - wm(i, k))) *
                  ((w_2(i, k) - wm(i, k)) * (w_2(i, k) - wm(i, k)))));
      });
}

void calc_wp2xp_pdf(int nz, int ngrdcol,                          // In
                    Real2D wm, Real2D xm, Real2D w_1, Real2D w_2, // In
                    Real2D x_1, Real2D x_2,                       // In
                    Real2D varnce_w_1, Real2D varnce_w_2,         // In
                    Real2D varnce_x_1, Real2D varnce_x_2,         // In
                    Real2D corr_w_x_1, Real2D corr_w_x_2,         // In
                    Real2D mixt_frac,                             // In
                    Real2D wp2xp) {
  // Description:
  // Calculates <w'^2 x'> by integrating over the PDF of w and x.  The integral
  // is:
  //
  // <w'^2 x'>
  // = INT(-inf:inf) INT(-inf:inf) ( w - <w> )^2 ( x - <x> ) P(w,x) dx dw;
  //
  // where <w> is the overall mean of w, <x> is the overall mean of x, and
  // P(w,x) is a two-component bivariate normal distribution of w and x.  The
  // integrated equation is:
  //
  // <w'^2 x'>
  // = mixt_frac * ( ( mu_x_1 - <x> ) * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
  //                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
  //                   * ( mu_w_1 - <w> ) )
  //   + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )
  //                           * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 )
  //                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
  //                             * ( mu_w_2 - <w> ) );
  //
  // where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
  // of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
  // component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
  // the standard deviation of w in the 1st PDF component, sigma_w_2 is the
  // standard deviation of w in the 2nd PDF component, sigma_x_1 is the
  // standard deviation of x in the 1st PDF component, sigma_x_2 is the
  // standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
  // correlation of w and x in the 1st PDF component, corr_w_x_2 is the
  // correlation of w and x in the 2nd PDF component, and mixt_frac is the
  // mixture fraction, which is the weight of the 1st PDF component.

  // References:
  //-----------------------------------------------------------------------

  // Arrays are indexed (column, level).
  // wm: Mean of w (overall) [m/s]
  // xm: Mean of x (overall) [units vary]
  // w_1: Mean of w (1st PDF component) [m/s]
  // w_2: Mean of w (2nd PDF component) [m/s]
  // x_1: Mean of x (1st PDF component) [units vary]
  // x_2: Mean of x (2nd PDF component) [units vary]
  // varnce_w_1: Variance of w (1st PDF component) [m^2/s^2]
  // varnce_w_2: Variance of w (2nd PDF component) [m^2/s^2]
  // varnce_x_1: Variance of x (1st PDF component) [(units vary)^2]
  // varnce_x_2: Variance of x (2nd PDF component) [(units vary)^2]
  // corr_w_x_1: Correlation of w and x (1st PDF comp.) [-]
  // corr_w_x_2: Correlation of w and x (2nd PDF comp.) [-]
  // mixt_frac: Mixture fraction [-]
  // wp2xp: <w'^2 x'> [m^2/s^2 (units vary)]

  // Calculate <w'^2 x'> by integrating over the PDF.
  Kokkos::parallel_for(
      "calc_wp2xp_pdf", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        wp2xp(i, k) =
            mixt_frac(i, k) *
                ((((w_1(i, k) - wm(i, k)) * (w_1(i, k) - wm(i, k))) + varnce_w_1(i, k)) *
                     (x_1(i, k) - xm(i, k)) +
                 two * corr_w_x_1(i, k) * Kokkos::sqrt(varnce_w_1(i, k) * varnce_x_1(i, k)) *
                     (w_1(i, k) - wm(i, k))) +
            (one - mixt_frac(i, k)) *
                ((((w_2(i, k) - wm(i, k)) * (w_2(i, k) - wm(i, k))) + varnce_w_2(i, k)) *
                     (x_2(i, k) - xm(i, k)) +
                 two * corr_w_x_2(i, k) * Kokkos::sqrt(varnce_w_2(i, k) * varnce_x_2(i, k)) *
                     (w_2(i, k) - wm(i, k)));
      });
}

void calc_wpxp2_pdf(int nz, int ngrdcol,                  // In
                    Real2D wm, Real2D xm, Real2D w_1,     // In
                    Real2D w_2, Real2D x_1, Real2D x_2,   // In
                    Real2D varnce_w_1, Real2D varnce_w_2, // In
                    Real2D varnce_x_1, Real2D varnce_x_2, // In
                    Real2D corr_w_x_1, Real2D corr_w_x_2, // In
                    Real2D mixt_frac,                     // In
                    Real2D wpxp2) {
  // Description:
  // Calculates <w'x'^2> by integrating over the PDF of w and x.  The integral
  // is:
  //
  // <w'x'^2>
  // = INT(-inf:inf) INT(-inf:inf) ( w - <w> ) ( x - <x> )^2 P(w,x) dx dw;
  //
  // where <w> is the overall mean of w, <x> is the overall mean of x, and
  // P(w,x) is a two-component bivariate normal distribution of w and x.  The
  // integrated equation is:
  //
  // <w'x'^2>
  // = mixt_frac * ( ( mu_w_1 - <w> ) * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
  //                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
  //                   * ( mu_x_1 - <x> ) )
  //   + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )
  //                           * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 )
  //                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
  //                             * ( mu_x_2 - <x> ) );
  //
  // where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
  // of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
  // component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
  // the standard deviation of w in the 1st PDF component, sigma_w_2 is the
  // standard deviation of w in the 2nd PDF component, sigma_x_1 is the
  // standard deviation of x in the 1st PDF component, sigma_x_2 is the
  // standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
  // correlation of w and x in the 1st PDF component, corr_w_x_2 is the
  // correlation of w and x in the 2nd PDF component, and mixt_frac is the
  // mixture fraction, which is the weight of the 1st PDF component.

  // References:
  //-----------------------------------------------------------------------

  // Arrays are indexed (column, level).
  // wm: Mean of w (overall) [m/s]
  // xm: Mean of x (overall) [units vary]
  // w_1: Mean of w (1st PDF component) [m/s]
  // w_2: Mean of w (2nd PDF component) [m/s]
  // x_1: Mean of x (1st PDF component) [units vary]
  // x_2: Mean of x (2nd PDF component) [units vary]
  // varnce_w_1: Variance of w (1st PDF component) [m^2/s^2]
  // varnce_w_2: Variance of w (2nd PDF component) [m^2/s^2]
  // varnce_x_1: Variance of x (1st PDF component) [(units vary)^2]
  // varnce_x_2: Variance of x (2nd PDF component) [(units vary)^2]
  // corr_w_x_1: Correlation of w and x (1st PDF comp.) [-]
  // corr_w_x_2: Correlation of w and x (2nd PDF comp.) [-]
  // mixt_frac: Mixture fraction [-]
  // wpxp2: <w'x'^2> [m/s (units vary)^2]

  // Calculate <w'x'^2> by integrating over the PDF.
  Kokkos::parallel_for(
      "calc_wpxp2_pdf", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        wpxp2(i, k) =
            mixt_frac(i, k) *
                ((w_1(i, k) - wm(i, k)) *
                     (((x_1(i, k) - xm(i, k)) * (x_1(i, k) - xm(i, k))) + varnce_x_1(i, k)) +
                 two * corr_w_x_1(i, k) * Kokkos::sqrt(varnce_w_1(i, k) * varnce_x_1(i, k)) *
                     (x_1(i, k) - xm(i, k))) +
            (one - mixt_frac(i, k)) *
                ((w_2(i, k) - wm(i, k)) *
                     (((x_2(i, k) - xm(i, k)) * (x_2(i, k) - xm(i, k))) + varnce_x_2(i, k)) +
                 two * corr_w_x_2(i, k) * Kokkos::sqrt(varnce_w_2(i, k) * varnce_x_2(i, k)) *
                     (x_2(i, k) - xm(i, k)));
      });
}

void calc_wpxpyp_pdf(int nz, int ngrdcol,                                     // In
                     Real2D wm, Real2D xm, Real2D ym, Real2D w_1, Real2D w_2, // In
                     Real2D x_1, Real2D x_2,                                  // In
                     Real2D y_1, Real2D y_2,                                  // In
                     Real2D varnce_w_1, Real2D varnce_w_2,                    // In
                     Real2D varnce_x_1, Real2D varnce_x_2,                    // In
                     Real2D varnce_y_1, Real2D varnce_y_2,                    // In
                     Real2D corr_w_x_1, Real2D corr_w_x_2,                    // In
                     Real2D corr_w_y_1, Real2D corr_w_y_2,                    // In
                     Real2D corr_x_y_1, Real2D corr_x_y_2,                    // In
                     Real2D mixt_frac,                                        // In
                     Real2D wpxpyp) {
  // Description:
  // Calculates <w'x'y'> by integrating over the PDF of w, x, and y.  The
  // integral is:
  //
  // <w'x'y'>
  // = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
  //   ( w - <w> ) ( x - <x> ) ( y - <y> ) P(w,x,y) dy dx dw;
  //
  // where <w> is the overall mean of w, <x> is the overall mean of x, <y> is
  // the overall mean of y, and P(w,x,y) is a two-component trivariate normal
  // distribution of w, x, and y.  The integrated equation is:
  //
  // <w'x'y'>
  // = mixt_frac
  //   * ( ( mu_w_1 - <w> ) * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
  //       + corr_x_y_1 * sigma_x_1 * sigma_y_1 * ( mu_w_1 - <w> )
  //       + corr_w_y_1 * sigma_w_1 * sigma_y_1 * ( mu_x_1 - <x> )
  //       + corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_y_1 - <y> ) )
  //   + ( 1 - mixt_frac )
  //     * ( ( mu_w_2 - <w> ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
  //         + corr_x_y_2 * sigma_x_2 * sigma_y_2 * ( mu_w_2 - <w> )
  //         + corr_w_y_2 * sigma_w_2 * sigma_y_2 * ( mu_x_2 - <x> )
  //         + corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_y_2 - <y> ) );
  //
  // where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
  // of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
  // component, mu_x_2 is the mean of x in the 2nd PDF component, mu_y_1 is the
  // mean of y in the 1st PDF component, mu_y_2 is the mean of y in the 2nd PDF
  // component, sigma_w_1 is the standard deviation of w in the 1st PDF
  // component, sigma_w_2 is the standard deviation of w in the 2nd PDF
  // component, sigma_x_1 is the standard deviation of x in the 1st PDF
  // component, sigma_x_2 is the standard deviation of x in the 2nd PDF
  // component, sigma_y_1 is the standard deviation of y in the 1st PDF
  // component, sigma_y_2 is the standard deviation of y in the 2nd PDF
  // component, corr_w_x_1 is the correlation of w and x in the 1st PDF
  // component, corr_w_x_2 is the correlation of w and x in the 2nd PDF
  // component, corr_w_y_1 is the correlation of w and y in the 1st PDF
  // component, corr_w_y_2 is the correlation of w and y in the 2nd PDF
  // component, corr_x_y_1 is the correlation of x and y in the 1st PDF
  // component, corr_x_y_2 is the correlation of x and y in the 2nd PDF
  // component, and mixt_frac is the mixture fraction, which is the weight of
  // the 1st PDF component.

  // References:
  //-----------------------------------------------------------------------

  // Arrays are indexed (column, level).
  // wm: Mean of w (overall) [m/s]
  // xm: Mean of x (overall) [x units]
  // ym: Mean of y (overall) [y units]
  // w_1: Mean of w (1st PDF component) [m/s]
  // w_2: Mean of w (2nd PDF component) [m/s]
  // x_1: Mean of x (1st PDF component) [x units]
  // x_2: Mean of x (2nd PDF component) [x units]
  // y_1: Mean of y (1st PDF component) [y units]
  // y_2: Mean of y (2nd PDF component) [y units]
  // varnce_w_1: Variance of w (1st PDF component) [m^2/s^2]
  // varnce_w_2: Variance of w (2nd PDF component) [m^2/s^2]
  // varnce_x_1: Variance of x (1st PDF component) [(x units)^2]
  // varnce_x_2: Variance of x (2nd PDF component) [(x units)^2]
  // varnce_y_1: Variance of y (1st PDF component) [(y units)^2]
  // varnce_y_2: Variance of y (2nd PDF component) [(y units)^2]
  // corr_w_x_1: Correlation of w and x (1st PDF component) [-]
  // corr_w_x_2: Correlation of w and x (2nd PDF component) [-]
  // corr_w_y_1: Correlation of w and y (1st PDF component) [-]
  // corr_w_y_2: Correlation of w and y (2nd PDF component) [-]
  // corr_x_y_1: Correlation of x and y (1st PDF component) [-]
  // corr_x_y_2: Correlation of x and y (2nd PDF component) [-]
  // mixt_frac: Mixture fraction [-]
  // wpxpyp: <w'x'y'> [m/s (units vary)]

  // Calculate <w'x'y'> by integrating over the PDF.
  Kokkos::parallel_for(
      "calc_wpxpyp_pdf", ngrdcol * nz, KOKKOS_LAMBDA(const int cell) {
        const int i = cell % ngrdcol;
        const int k = cell / ngrdcol;
        wpxpyp(i, k) =
            mixt_frac(i, k) *
                ((w_1(i, k) - wm(i, k)) * (x_1(i, k) - xm(i, k)) * (y_1(i, k) - ym(i, k)) +
                 corr_x_y_1(i, k) * Kokkos::sqrt(varnce_x_1(i, k) * varnce_y_1(i, k)) *
                     (w_1(i, k) - wm(i, k)) +
                 corr_w_y_1(i, k) * Kokkos::sqrt(varnce_w_1(i, k) * varnce_y_1(i, k)) *
                     (x_1(i, k) - xm(i, k)) +
                 corr_w_x_1(i, k) * Kokkos::sqrt(varnce_w_1(i, k) * varnce_x_1(i, k)) *
                     (y_1(i, k) - ym(i, k))) +
            (one - mixt_frac(i, k)) *
                ((w_2(i, k) - wm(i, k)) * (x_2(i, k) - xm(i, k)) * (y_2(i, k) - ym(i, k)) +
                 corr_x_y_2(i, k) * Kokkos::sqrt(varnce_x_2(i, k) * varnce_y_2(i, k)) *
                     (w_2(i, k) - wm(i, k)) +
                 corr_w_y_2(i, k) * Kokkos::sqrt(varnce_w_2(i, k) * varnce_y_2(i, k)) *
                     (x_2(i, k) - xm(i, k)) +
                 corr_w_x_2(i, k) * Kokkos::sqrt(varnce_w_2(i, k) * varnce_x_2(i, k)) *
                     (y_2(i, k) - ym(i, k)));
      });
}

} // namespace clubb_kokkos
