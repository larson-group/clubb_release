// Diagnose cloudy updraft and downdraft velocities from the two PDF components.
//
// calc_w_up_in_cloud returns two mean vertical velocities [m/s] and their
// cloudy fractions [-]. Its seven input fields and four caller-owned output
// fields use (column, level) LayoutLeft Views. The caller selects the supported
// ADG1, ADG2 or new-hybrid PDF; this routine has no state or statistics updates.
// Isolated source draft for merging into pdf_closure_module; runtime unverified.
#pragma once

#include "array_types.hpp"

namespace clubb_kokkos {

void calc_w_up_in_cloud(int nz, int ngrdcol,                                          // In
                        Real2D mixt_frac, Real2D cloud_frac_1, Real2D cloud_frac_2,   // In
                        Real2D w_1, Real2D w_2, Real2D varnce_w_1, Real2D varnce_w_2, // In
                        Real2D w_up_in_cloud, Real2D w_down_in_cloud,                 // Out
                        Real2D cloudy_updraft_frac, Real2D cloudy_downdraft_frac);    // Out

} // namespace clubb_kokkos
