// Advance CLUBB's prognostic moments and diagnose its PDF for one timestep.
//
// advance_clubb_core follows src/CLUBB_core/advance_clubb_core_module.F90.
// Pass column-first device Views, initialized grid/configuration/PDF values and
// preallocated outputs. Scratch Views live only for this call. Statistics use
// the case's native f2k context until the stats backend is ported.
// The host-facing clubb_api_module.hpp stages driver arrays when needed.
// The first array index is the column; the second is a thermodynamic (nzt)
// or momentum (nzm) level, and a third index selects scalar/species fields.
// dt is in seconds. Input forcing fields use the corresponding state units
// per second. Outputs include diffusivities [m^2/s], mixing length [m],
// inverse dissipation time [1/s] and cloud/PDF diagnostics.
//
// References: Golaz et al. (2002), JAS 59, 3540-3551;
// https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:overview_clubb
#pragma once
#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/model_flags.hpp"
#include "CLUBB_core/parameters_tunable.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
#include "CLUBB_core/err_info_type_module.hpp"
#include "stats_netcdf.hpp"

namespace clubb_kokkos {
void advance_clubb_core(
    ScratchWorkspace& workspace,                                                        // InOut
    const grid& gr, int nzm, int nzt, int ngrdcol,                                      // In
    bool l_implemented, real_t dt, const Real1D& fcor, const Real1D& fcor_y,            // In
    const Real1D& sfc_elevation,                                                        // In
    int hydromet_dim,                                                                   // In
    int sclr_dim, const Real1D& sclr_tol, int edsclr_dim,                               // In
    const sclr_idx_type& sclr_idx,                                                      // In
    const Real2D& thlm_forcing, const Real2D& rtm_forcing, const Real2D& um_forcing,    // In
    const Real2D& vm_forcing,                                                           // In
    const Real3D& sclrm_forcing, const Real3D& edsclrm_forcing,                         // In
    const Real2D& wprtp_forcing,                                                        // In
    const Real2D& wpthlp_forcing, const Real2D& rtp2_forcing,                           // In
    const Real2D& thlp2_forcing,                                                        // In
    const Real2D& rtpthlp_forcing, const Real2D& wm_zm, const Real2D& wm_zt,            // In
    const Real1D& wpthlp_sfc, const Real1D& wprtp_sfc, const Real1D& upwp_sfc,          // In
    const Real1D& vpwp_sfc, const Real1D& p_sfc,                                        // In
    const Real2D& wpsclrp_sfc, const Real2D& wpedsclrp_sfc,                             // In
    const Real1D& upwp_sfc_pert, const Real1D& vpwp_sfc_pert,                           // In
    const Real2D& rtm_ref, const Real2D& thlm_ref, const Real2D& um_ref,                // In
    const Real2D& vm_ref, const Real2D& ug, const Real2D& vg,                           // In
    Real2D p_in_Pa,                                                                     // InOut
    const Real2D& rho_zm, const Real2D& rho,                                            // In
    Real2D exner,                                                                       // InOut
    const Real2D& rho_ds_zm, const Real2D& rho_ds_zt, const Real2D& invrs_rho_ds_zm,    // In
    const Real2D& invrs_rho_ds_zt, const Real2D& thv_ds_zm, const Real2D& thv_ds_zt,    // In
    const Bool1D& l_mix_rat_hm,                                                         // In
    const Real2D& rfrzm,                                                                // In
    const Real3D& wphydrometp, const Real3D& wp2hmp, const Real3D& rtphmp_zt,           // In
    const Real3D& thlphmp_zt,                                                           // In
    const Real1D& host_dx, const Real1D& host_dy,                                       // In
    const Real2D& clubb_params, const nu_vertical_res_dep& nu_vert_res_dep,             // In
    const Real1D& lmin,                                                                 // In
    const Real1D& mixt_frac_max_mag, real_t T0, real_t ts_nudge,                        // In
    real_t rtm_min, real_t rtm_nudge_max_altitude,                                      // In
    const clubb_config_flags_type& clubb_config_flags,                                  // In
    stats_type& stats,                                                                  // InOut
    Real2D um, Real2D vm, Real2D upwp, Real2D vpwp, Real2D up2, Real2D vp2, Real2D up3, // InOut
    Real2D vp3,                                                                         // InOut
    Real2D thlm, Real2D rtm, Real2D wprtp, Real2D wpthlp,                               // InOut
    Real2D wp2, Real2D wp3, Real2D rtp2, Real2D rtp3, Real2D thlp2, Real2D thlp3,       // InOut
    Real2D rtpthlp,                                                                     // InOut
    Real3D sclrm,                                                                       // InOut
    Real3D sclrp2, Real3D sclrp3, Real3D sclrprtp, Real3D sclrpthlp,                    // InOut
    Real3D wpsclrp, Real3D edsclrm,                                                     // InOut
    Real2D rcm, Real2D cloud_frac,                                                      // InOut
    Real2D wpthvp, Real2D wp2thvp, Real2D wp2up, Real2D rtpthvp, Real2D thlpthvp,       // InOut
    Real3D sclrpthvp,                                                                   // InOut
    Real2D wp2rtp, Real2D wp2thlp, Real2D uprcp, Real2D vprcp, Real2D rc_coef_zm,       // InOut
    Real2D wp4,                                                                         // InOut
    Real2D wpup2, Real2D wpvp2, Real2D wp2up2, Real2D wp2vp2, Real2D ice_supersat_frac, // InOut
    Real2D um_pert, Real2D vm_pert, Real2D upwp_pert, Real2D vpwp_pert,                 // InOut
    pdf_parameter& pdf_params, pdf_parameter& pdf_params_zm,                            // InOut
    implicit_coefs_terms& pdf_implicit_coefs_terms,                                     // InOut
    err_info_type& err_info,                                                            // InOut
    Real2D Kh_zm, Real2D Kh_zt,                                                         // Out
    Real2D thlprcp, Real2D wprcp, Real2D w_up_in_cloud, Real2D w_down_in_cloud,         // Out
    Real2D cloudy_updraft_frac, Real2D cloudy_downdraft_frac,                           // Out
    Real2D rcm_in_layer, Real2D cloud_cover, Real2D invrs_tau_zm,                       // Out
    Real2D Lscale);                                                                     // Out
} // namespace clubb_kokkos
