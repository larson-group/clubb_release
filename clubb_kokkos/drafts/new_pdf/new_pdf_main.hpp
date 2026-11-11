// Original two-component PDF driver for w, rt and thl fields.
//
// new_pdf_driver computes the setter mixture and component moments, then closes
// responder rt/thl parameters and packs implicit-coefficient/explicit-term fields.
// Skrt and Skthl may be clipped in place. Fields are column-first LayoutLeft;
// scratch belongs to the call and the caller owns implicit_coefs_terms storage.
// All three Native advection flags are retained. The only transient lower science
// dependency is the two cubic-root solves inside each responder's limit call.
// Isolated source draft: not compiled, integrated or numerically tested.
#pragma once
#include "array_types.hpp"
#include "pdf_parameter_module.hpp"
namespace clubb_kokkos {

void new_pdf_driver(int nz, int ngrdcol, Real2D wm, Real2D rtm, Real2D thlm, Real2D wp2,
                    Real2D rtp2, Real2D thlp2, Real2D Skw, Real2D wprtp, Real2D wpthlp,
                    Real2D rtpthlp, Real2D clubb_params, Real2D Skrt, Real2D Skthl, Real2D mu_w_1,
                    Real2D mu_w_2, Real2D mu_rt_1, Real2D mu_rt_2, Real2D mu_thl_1, Real2D mu_thl_2,
                    Real2D sigma_w_1_sqd, Real2D sigma_w_2_sqd, Real2D sigma_rt_1_sqd,
                    Real2D sigma_rt_2_sqd, Real2D sigma_thl_1_sqd, Real2D sigma_thl_2_sqd,
                    Real2D mixt_frac, implicit_coefs_terms& pdf_implicit_coefs_terms);

} // namespace clubb_kokkos
