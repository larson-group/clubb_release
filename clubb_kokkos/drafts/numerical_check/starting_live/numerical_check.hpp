// Check the CLUBB core's inputs and outputs for non-finite or negative fields.
//
// parameterization_check follows the native input-diagnostic sequence and updates
// caller-owned err_info. Any non-finite field marks all columns fatal. Negative
// positive-definite fields also mark all columns fatal, but the native beginning-
// of-step call resets those host-provided negative errors after reporting them.
// A cold-temperature warning covers the first ten thermodynamic levels (<190K).
// Arrays are column-first execution-space Views; field units follow core inputs.
// Other numerical_check routines remain explicit native calls until reached.
#pragma once

#include "array_types.hpp"
#include "err_info_type_module.hpp"
#include <string>

namespace clubb_kokkos {

void parameterization_check(
    int nzm, int nzt, int ngrdcol, int sclr_dim, int edsclr_dim,                              // In
    Real2D thlm_forcing, Real2D rtm_forcing, Real2D um_forcing,                               // In
    Real2D vm_forcing, Real2D wm_zm, Real2D wm_zt, Real2D p_in_Pa,                            // In
    Real2D rho_zm, Real2D rho, Real2D exner, Real2D rho_ds_zm,                                // In
    Real2D rho_ds_zt, Real2D invrs_rho_ds_zm, Real2D invrs_rho_ds_zt,                         // In
    Real2D thv_ds_zm, Real2D thv_ds_zt, Real1D wpthlp_sfc, Real1D wprtp_sfc, Real1D upwp_sfc, // In
    Real1D vpwp_sfc, Real1D p_sfc, Real2D um, Real2D upwp, Real2D vm, Real2D vpwp, Real2D up2,
    Real2D vp2,                                                                   // In
    Real2D rtm, Real2D wprtp, Real2D thlm, Real2D wpthlp, Real2D wp2, Real2D wp3, // In
    Real2D rtp2, Real2D thlp2, Real2D rtpthlp,                                    // In
    const std::string& prefix,                                                    // In
    Real2D wpsclrp_sfc, Real2D wpedsclrp_sfc, Real3D sclrm, Real3D wpsclrp,       // In
    Real3D sclrp2,                                                                // In
    Real3D sclrprtp, Real3D sclrpthlp, Real3D sclrm_forcing, Real3D edsclrm,      // In
    Real3D edsclrm_forcing,                                                       // In
    err_info_type& err_info);                                                     // InOut

} // namespace clubb_kokkos
