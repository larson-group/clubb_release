// Compute the PDF width subject to the responder-variable realizability limits.
//
// compute_sigma_sqd_w uses momentum-level variances and fluxes, interpolates
// vertical-velocity skewness from thermodynamic levels, then smooths the width
// vertically. Arrays are indexed (column, level) in the execution space.
// Grid interpolation uses the C++ grid operators in the same execution space.
// Reference: Equations for CLUBB, Eq. 22.
#pragma once

#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "grid_class.hpp"

namespace clubb_kokkos {

void compute_sigma_sqd_w(ScratchWorkspace& workspace,                       // InOut
                         int nzm, int nzt, int ngrdcol, const grid& gr,     // In
                         Real2D wp3, Real2D wp2, Real2D thlp2, Real2D rtp2, // In
                         Real2D up2, Real2D vp2, Real2D wpthlp, Real2D wprtp, Real2D upwp,
                         Real2D vpwp,              // In
                         Real2D clubb_params,      // In
                         bool l_predict_upwp_vpwp, // In
                         Real2D sigma_sqd_w);      // Out

} // namespace clubb_kokkos
