// Integrate fourth/third moments over the two-component normal PDF.
// Source-only draft: calc_wp4_pdf, calc_wp2xp_pdf, calc_wpxp2_pdf and calc_wpxpyp_pdf.
// Column-first View arguments retain Native field units; not compiled or tested.
#pragma once
#include "array_types.hpp"

namespace clubb_kokkos {

void calc_wp4_pdf(int nz, int ngrdcol,                  // In
                  Real2D wm, Real2D w_1, Real2D w_2,    // In
                  Real2D varnce_w_1, Real2D varnce_w_2, // In
                  Real2D mixt_frac,                     // In
                  Real2D wp4);

void calc_wp2xp_pdf(int nz, int ngrdcol,                          // In
                    Real2D wm, Real2D xm, Real2D w_1, Real2D w_2, // In
                    Real2D x_1, Real2D x_2,                       // In
                    Real2D varnce_w_1, Real2D varnce_w_2,         // In
                    Real2D varnce_x_1, Real2D varnce_x_2,         // In
                    Real2D corr_w_x_1, Real2D corr_w_x_2,         // In
                    Real2D mixt_frac,                             // In
                    Real2D wp2xp);

void calc_wpxp2_pdf(int nz, int ngrdcol,                  // In
                    Real2D wm, Real2D xm, Real2D w_1,     // In
                    Real2D w_2, Real2D x_1, Real2D x_2,   // In
                    Real2D varnce_w_1, Real2D varnce_w_2, // In
                    Real2D varnce_x_1, Real2D varnce_x_2, // In
                    Real2D corr_w_x_1, Real2D corr_w_x_2, // In
                    Real2D mixt_frac,                     // In
                    Real2D wpxp2);

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
                     Real2D wpxpyp);

} // namespace clubb_kokkos
