// Diagnose surface velocity/scalar variances and covariances for CLUBB.
//
// calc_sfc_varnce updates only the Native lower momentum level using surface
// fluxes, mean winds, parcel length and per-column parameters. Passive scalar
// emulation, optional convective-depth diagnosis and the source-constant-false
// Andre1978 branch are retained. Fields use column-first LayoutLeft Views;
// sclr_idx/grid metadata retain Native one-based values. Stats remains Native.
// Reference: Andre et al.1978, JAS35,1861–1883.
#pragma once

#include "scratch_workspace.hpp"
#include "array_types.hpp"
#include "array_index.hpp"
#include "grid_class.hpp"
#include "err_info_type_module.hpp"
#include "stats_netcdf.hpp"

namespace clubb_kokkos {

void calc_sfc_varnce(ScratchWorkspace& workspace, int nzm, int nzt, int ngrdcol, int sclr_dim,
                     const sclr_idx_type& sclr_idx,                   // In
                     const grid& gr, real_t dt, Real1D sfc_elevation, // In
                     Real1D upwp_sfc, Real1D vpwp_sfc, Real2D wpthlp,
                     Real1D wprtp_sfc, // In
                     Real2D um, Real2D vm, Real2D Lscale_up,
                     Real2D wpsclrp_sfc,                               // In
                     Real2D lhs_splat_wp2, Real2D tau_zm,              // In
                     bool l_vary_convect_depth,                        // In
                     real_t T0,                                        // In
                     Real1D up2_sfc_coef,                              // In
                     Real1D a_const,                                   // In
                     stats_type& stats,                                // InOut
                     Real2D wp2, Real2D up2, Real2D vp2,               // InOut
                     Real2D thlp2, Real2D rtp2, Real2D rtpthlp,        // InOut
                     Real3D sclrp2, Real3D sclrprtp, Real3D sclrpthlp, // InOut
                     err_info_type& err_info)                          // InOut
    ;
} // namespace clubb_kokkos
