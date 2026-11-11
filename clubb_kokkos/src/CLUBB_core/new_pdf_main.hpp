// Original two-component PDF driver for w, rt and thl fields.
//
// new_pdf_driver computes the setter mixture and component moments, then closes
// responder rt/thl parameters and writes implicit coefficients and explicit terms
// directly into the caller-owned implicit_coefs_terms fields. Skrt and Skthl may
// be clipped in place. Independent columns/levels share each execution-space
// kernel; the explicit workspace owns scoped column-first LayoutLeft scratch.
// All three Native advection flags are retained. Responder limits use the
// execution-space cubic root solver in calc_roots.hpp.
#pragma once
#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "pdf_parameter_module.hpp"
namespace clubb_kokkos {

void new_pdf_driver(ScratchWorkspace& workspace,                                         // InOut
                    int nz, int ngrdcol, Real2D wm, Real2D rtm, Real2D thlm, Real2D wp2, // In
                    Real2D rtp2, Real2D thlp2,                                           // In
                    Real2D Skw,                                                          // In
                    Real2D wprtp, Real2D wpthlp, Real2D rtpthlp,                         // In
                    Real2D clubb_params,                                                 // In
                    Real2D Skrt, Real2D Skthl,                                           // InOut
                    Real2D mu_w_1, Real2D mu_w_2,                                        // Out
                    Real2D mu_rt_1, Real2D mu_rt_2,                                      // Out
                    Real2D mu_thl_1, Real2D mu_thl_2,                                    // Out
                    Real2D sigma_w_1_sqd, Real2D sigma_w_2_sqd,                          // Out
                    Real2D sigma_rt_1_sqd, Real2D sigma_rt_2_sqd,                        // Out
                    Real2D sigma_thl_1_sqd, Real2D sigma_thl_2_sqd,                      // Out
                    Real2D mixt_frac,                                                    // Out
                    implicit_coefs_terms& pdf_implicit_coefs_terms                       // InOut
);

} // namespace clubb_kokkos
