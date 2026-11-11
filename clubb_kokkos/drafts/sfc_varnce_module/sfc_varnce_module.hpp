// Diagnose surface velocity/scalar variances and covariances for CLUBB.
//
// calc_sfc_varnce updates only the Native lower momentum level using surface
// fluxes, mean winds, parcel length and per-column parameters. Passive scalar
// emulation, optional convective-depth diagnosis and the source-constant-false
// Andre1978 branch are retained. Fields use column-first LayoutLeft Views;
// sclr_idx/grid metadata retain Native one-based values. Stats remains Native.
// See INTEGRATION.md for the temporary public diagnostic checker boundary.
// Reference: Andre et al.1978, JAS35,1861–1883.
#pragma once

#include "array_types.hpp"
#include "array_index.hpp"
#include "grid_class.hpp"
#include "err_info_type_module.hpp"
#include "f2k/stats_bridge.hpp"

namespace clubb_kokkos {

void calc_sfc_varnce(int nzm, int nzt, int ngrdcol, int sclr_dim, const sclr_idx_type& sclr_idx,
                     const grid& gr, real_t dt, Real1D sfc_elevation, Real1D upwp_sfc,
                     Real1D vpwp_sfc, Real2D wpthlp, Real1D wprtp_sfc, Real2D um, Real2D vm,
                     Real2D Lscale_up, Real2D wpsclrp_sfc, Real2D lhs_splat_wp2, Real2D tau_zm,
                     bool l_vary_convect_depth, real_t T0, Real1D up2_sfc_coef, Real1D a_const,
                     f2k::stats_type& stats, Real2D wp2, Real2D up2, Real2D vp2, Real2D thlp2,
                     Real2D rtp2, Real2D rtpthlp, Real3D sclrp2, Real3D sclrprtp, Real3D sclrpthlp,
                     err_info_type& err_info);
} // namespace clubb_kokkos
