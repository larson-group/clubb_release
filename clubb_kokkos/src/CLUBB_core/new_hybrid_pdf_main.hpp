// Close hybrid PDF parameters and implicit moment terms for column/level fields.
//
// new_hybrid_pdf_driver computes the vertical-velocity mixture, then rt, thl,
// horizontal-wind and passive-scalar component means/variances. Responder
// skewnesses may be clipped. The caller-owned implicit_coefs_terms receives
// coefficients/explicit terms selected by the three Native advection flags.
// Independent columns/levels share each execution-space kernel; species remain
// ordered within each cell. Model Views use column-first LayoutLeft, and the
// explicit workspace owns scoped rank-2/rank-3 scratch for the batch.
#pragma once
#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "pdf_parameter_module.hpp"
namespace clubb_kokkos {

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
                           implicit_coefs_terms& pdf_implicit_coefs_terms);

} // namespace clubb_kokkos
