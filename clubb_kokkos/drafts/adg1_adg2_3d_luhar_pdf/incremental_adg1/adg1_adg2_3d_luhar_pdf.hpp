// Compute ADG1 PDF component means and variances from grid-box moments.
//
// ADG1_pdf_driver first closes vertical velocity, then computes responder
// parameters for total water, liquid-water potential temperature, horizontal
// winds and optional passive scalars. Component means retain each field's
// units; variances have squared units and alpha/mixture factors are dimensionless.
// All profile Views use execution-space, column-first LayoutLeft storage;
// scalar fields add a final species dimension. The caller owns outputs and
// err_info. Vertical-velocity closure remains an explicit temporary f2k call.
// Source: src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90:ADG1_pdf_driver.
#pragma once

#include "array_types.hpp"
#include "err_info_type_module.hpp"

namespace clubb_kokkos {

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
                     Real3D varnce_sclr_2, Real3D alpha_sclr);

} // namespace clubb_kokkos
