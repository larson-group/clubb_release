// Accumulate CLUBB state, hydrometeor and SILHS diagnostics for sampled timesteps.
//
// stats_accumulate records model profiles, column integrals, vertical averages
// and conservation residuals. stats_accumulate_hydromet_api records hydrometeor
// profiles and water paths. stats_accumulate_lh_tend records sampled microphysics
// tendencies and Kessler estimates. Calls preserve Native sampling and registry
// guards; stats owns its buffers and output state.
//
// Scientific fields use Native column-first LayoutLeft order. The implementation
// works on call-local HostSpace snapshots; generic overloads pack host, strided
// or default-execution fields through the existing array owner only when sampled.
// Existing C++ scientific kernels remain explicit lower calls. Grid/species
// metadata retains Native one-based indices; stats column arguments are zero-based.
#pragma once

#include "scratch_workspace.hpp"
#include "CLUBB_core/array_types.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/stats_netcdf.hpp"

namespace clubb_kokkos {

void stats_accumulate(ScratchWorkspace& workspace,                                 // InOut
                      int nzm, int nzt, int ngrdcol, int sclr_dim, int edsclr_dim, // In
                      const grid& gr, real_t dt,                                   // In
                      bool l_implemented, bool l_host_applies_sfc_fluxes,          // In
                      bool l_stability_correct_tau_zm,                             // In
                      StatsHostReal2D clubb_params,                                // In
                      StatsHostReal2D um, StatsHostReal2D vm, StatsHostReal2D upwp,
                      StatsHostReal2D vpwp, StatsHostReal2D up2, StatsHostReal2D vp2, // In
                      StatsHostReal2D thlm, StatsHostReal2D rtm, StatsHostReal2D thlm_before,
                      StatsHostReal2D rtm_before, StatsHostReal2D thlm_forcing,
                      StatsHostReal2D rtm_forcing, // In
                      StatsHostReal1D wpthlp_sfc, StatsHostReal1D wprtp_sfc, StatsHostReal2D wprtp,
                      StatsHostReal2D wpthlp, // In
                      StatsHostReal2D wp2, StatsHostReal2D wp3, StatsHostReal2D rtp2,
                      StatsHostReal2D rtp3, StatsHostReal2D thlp2, StatsHostReal2D thlp3, // In
                      StatsHostReal2D rtpthlp,                                            // In
                      StatsHostReal2D p_in_Pa, StatsHostReal2D exner, StatsHostReal2D rho,
                      StatsHostReal2D rho_zm, // In
                      StatsHostReal2D rho_ds_zm, StatsHostReal2D rho_ds_zt,
                      StatsHostReal2D thv_ds_zm,
                      StatsHostReal2D thv_ds_zt,                                           // In
                      StatsHostReal2D wm_zt, StatsHostReal2D wm_zm, StatsHostReal2D rcm,   // In
                      StatsHostReal2D cloud_frac,                                          // In
                      StatsHostReal2D thvm, StatsHostReal2D ug, StatsHostReal2D vg,        // In
                      StatsHostReal2D ddzt_umvm_sqd, StatsHostReal2D stability_correction, // In
                      StatsHostReal2D Kh_zt,                                               // In
                      StatsHostReal2D rsat,                                                // In
                      StatsHostReal2D Kh_zm,                                               // In
                      StatsHostReal2D em,                                                  // In
                      StatsHostReal3D sclrm, StatsHostReal3D sclrp2,                       // In
                      StatsHostReal3D sclrprtp, StatsHostReal3D sclrpthlp,
                      StatsHostReal3D sclrm_forcing, // In
                      StatsHostReal3D wpsclrp, StatsHostReal3D wpedsclrp,
                      StatsHostReal3D edsclrm,         // In
                      StatsHostReal3D edsclrm_forcing, // In
                      int saturation_formula,          // In
                      stats_type& stats                // InOut
);

void stats_accumulate_hydromet_api(ScratchWorkspace& workspace,                   // InOut
                                   const grid& gr, int ngrdcol, int hydromet_dim, // In
                                   const hm_metadata_type& hm_metadata,
                                   StatsHostReal3D hydromet,  // In
                                   StatsHostReal2D rho_ds_zt, // In
                                   stats_type& stats          // InOut
);

void stats_accumulate_lh_tend(ScratchWorkspace& workspace, // InOut
                              const grid& gr, int ngrdcol, int hydromet_dim,
                              const hm_metadata_type& hm_metadata,                       // In
                              StatsHostReal3D lh_hydromet_mc, StatsHostReal2D lh_Ncm_mc, // In
                              StatsHostReal2D lh_thlm_mc, StatsHostReal2D lh_rvm_mc,
                              StatsHostReal2D lh_rcm_mc, // In
                              StatsHostReal2D lh_AKm, StatsHostReal2D AKm, StatsHostReal2D AKstd,
                              StatsHostReal2D AKstd_cld, // In
                              StatsHostReal2D lh_rcm_avg, StatsHostReal2D AKm_rcm,
                              StatsHostReal2D AKm_rcc, // In
                              stats_type& stats        // InOut
);

// Host I/O boundaries pack arbitrary caller Views without changing the Native
// leading arguments. No output field or caller allocation is replaced.
template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11,
          class Array12, class Array13, class Array14, class Array15, class Array16, class Array17,
          class Array18, class Array19, class Array20, class Array21, class Array22, class Array23,
          class Array24, class Array25, class Array26, class Array27, class Array28, class Array29,
          class Array30, class Array31, class Array32, class Array33, class Array34, class Array35,
          class Array36, class Array37, class Array38, class Array39, class Array40, class Array41,
          class Array42, class Array43, class Array44, class Array45, class Array46, class Array47,
          class Array48, class Array49, class Array50, class Array51, class Array52, class Array53>
void stats_accumulate(ScratchWorkspace& workspace,                                 // InOut
                      int nzm, int nzt, int ngrdcol, int sclr_dim, int edsclr_dim, // In
                      const grid& gr, real_t dt,                                   // In
                      bool l_implemented, bool l_host_applies_sfc_fluxes,          // In
                      bool l_stability_correct_tau_zm,                             // In
                      Array0 clubb_params,                                         // In
                      Array1 um, Array2 vm, Array3 upwp, Array4 vpwp, Array5 up2, Array6 vp2, // In
                      Array7 thlm, Array8 rtm, Array9 thlm_before, Array10 rtm_before,
                      Array11 thlm_forcing, Array12 rtm_forcing,                            // In
                      Array13 wpthlp_sfc, Array14 wprtp_sfc, Array15 wprtp, Array16 wpthlp, // In
                      Array17 wp2, Array18 wp3, Array19 rtp2, Array20 rtp3, Array21 thlp2,
                      Array22 thlp3,                                               // In
                      Array23 rtpthlp,                                             // In
                      Array24 p_in_Pa, Array25 exner, Array26 rho, Array27 rho_zm, // In
                      Array28 rho_ds_zm, Array29 rho_ds_zt, Array30 thv_ds_zm,
                      Array31 thv_ds_zt,                                          // In
                      Array32 wm_zt, Array33 wm_zm, Array34 rcm,                  // In
                      Array35 cloud_frac,                                         // In
                      Array36 thvm, Array37 ug, Array38 vg,                       // In
                      Array39 ddzt_umvm_sqd, Array40 stability_correction,        // In
                      Array41 Kh_zt,                                              // In
                      Array42 rsat,                                               // In
                      Array43 Kh_zm,                                              // In
                      Array44 em,                                                 // In
                      Array45 sclrm, Array46 sclrp2,                              // In
                      Array47 sclrprtp, Array48 sclrpthlp, Array49 sclrm_forcing, // In
                      Array50 wpsclrp, Array51 wpedsclrp, Array52 edsclrm,        // In
                      Array53 edsclrm_forcing,                                    // In
                      int saturation_formula,                                     // In
                      stats_type& stats                                           // InOut
) {
  // Stage scientific inputs only for a Native sampled timestep.
  auto workspace_scope = workspace.scope();
  if (!stats.l_sample)
    return;
  static_assert(Array0::rank == 2);
  static_assert(Array1::rank == 2);
  static_assert(Array2::rank == 2);
  static_assert(Array3::rank == 2);
  static_assert(Array4::rank == 2);
  static_assert(Array5::rank == 2);
  static_assert(Array6::rank == 2);
  static_assert(Array7::rank == 2);
  static_assert(Array8::rank == 2);
  static_assert(Array9::rank == 2);
  static_assert(Array10::rank == 2);
  static_assert(Array11::rank == 2);
  static_assert(Array12::rank == 2);
  static_assert(Array13::rank == 1);
  static_assert(Array14::rank == 1);
  static_assert(Array15::rank == 2);
  static_assert(Array16::rank == 2);
  static_assert(Array17::rank == 2);
  static_assert(Array18::rank == 2);
  static_assert(Array19::rank == 2);
  static_assert(Array20::rank == 2);
  static_assert(Array21::rank == 2);
  static_assert(Array22::rank == 2);
  static_assert(Array23::rank == 2);
  static_assert(Array24::rank == 2);
  static_assert(Array25::rank == 2);
  static_assert(Array26::rank == 2);
  static_assert(Array27::rank == 2);
  static_assert(Array28::rank == 2);
  static_assert(Array29::rank == 2);
  static_assert(Array30::rank == 2);
  static_assert(Array31::rank == 2);
  static_assert(Array32::rank == 2);
  static_assert(Array33::rank == 2);
  static_assert(Array34::rank == 2);
  static_assert(Array35::rank == 2);
  static_assert(Array36::rank == 2);
  static_assert(Array37::rank == 2);
  static_assert(Array38::rank == 2);
  static_assert(Array39::rank == 2);
  static_assert(Array40::rank == 2);
  static_assert(Array41::rank == 2);
  static_assert(Array42::rank == 2);
  static_assert(Array43::rank == 2);
  static_assert(Array44::rank == 2);
  static_assert(Array45::rank == 3);
  static_assert(Array46::rank == 3);
  static_assert(Array47::rank == 3);
  static_assert(Array48::rank == 3);
  static_assert(Array49::rank == 3);
  static_assert(Array50::rank == 3);
  static_assert(Array51::rank == 3);
  static_assert(Array52::rank == 3);
  static_assert(Array53::rank == 3);
  stats_accumulate(
      workspace,

      nzm, nzt, ngrdcol, sclr_dim, edsclr_dim, gr, dt, l_implemented, l_host_applies_sfc_fluxes,
      l_stability_correct_tau_zm, to_host_space(workspace, clubb_params),
      to_host_space(workspace, um), to_host_space(workspace, vm), to_host_space(workspace, upwp),
      to_host_space(workspace, vpwp), to_host_space(workspace, up2), to_host_space(workspace, vp2),
      to_host_space(workspace, thlm), to_host_space(workspace, rtm),
      to_host_space(workspace, thlm_before), to_host_space(workspace, rtm_before),
      to_host_space(workspace, thlm_forcing), to_host_space(workspace, rtm_forcing),
      to_host_space(workspace, wpthlp_sfc), to_host_space(workspace, wprtp_sfc),
      to_host_space(workspace, wprtp), to_host_space(workspace, wpthlp),
      to_host_space(workspace, wp2), to_host_space(workspace, wp3), to_host_space(workspace, rtp2),
      to_host_space(workspace, rtp3), to_host_space(workspace, thlp2),
      to_host_space(workspace, thlp3), to_host_space(workspace, rtpthlp),
      to_host_space(workspace, p_in_Pa), to_host_space(workspace, exner),
      to_host_space(workspace, rho), to_host_space(workspace, rho_zm),
      to_host_space(workspace, rho_ds_zm), to_host_space(workspace, rho_ds_zt),
      to_host_space(workspace, thv_ds_zm), to_host_space(workspace, thv_ds_zt),
      to_host_space(workspace, wm_zt), to_host_space(workspace, wm_zm),
      to_host_space(workspace, rcm), to_host_space(workspace, cloud_frac),
      to_host_space(workspace, thvm), to_host_space(workspace, ug), to_host_space(workspace, vg),
      to_host_space(workspace, ddzt_umvm_sqd), to_host_space(workspace, stability_correction),
      to_host_space(workspace, Kh_zt), to_host_space(workspace, rsat),
      to_host_space(workspace, Kh_zm), to_host_space(workspace, em),
      to_host_space(workspace, sclrm), to_host_space(workspace, sclrp2),
      to_host_space(workspace, sclrprtp), to_host_space(workspace, sclrpthlp),
      to_host_space(workspace, sclrm_forcing), to_host_space(workspace, wpsclrp),
      to_host_space(workspace, wpedsclrp), to_host_space(workspace, edsclrm),
      to_host_space(workspace, edsclrm_forcing), saturation_formula, stats);
}

template <class Array0, class Array1>
void stats_accumulate_hydromet_api(ScratchWorkspace& workspace,                          // InOut
                                   const grid& gr, int ngrdcol, int hydromet_dim,        // In
                                   const hm_metadata_type& hm_metadata, Array0 hydromet, // In
                                   Array1 rho_ds_zt,                                     // In
                                   stats_type& stats                                     // InOut
) {
  // Stage scientific inputs only for a Native sampled timestep.
  auto workspace_scope = workspace.scope();
  if (!stats.l_sample)
    return;
  static_assert(Array0::rank == 3);
  static_assert(Array1::rank == 2);
  stats_accumulate_hydromet_api(workspace, gr, ngrdcol, hydromet_dim, hm_metadata,
                                to_host_space(workspace, hydromet),
                                to_host_space(workspace, rho_ds_zt), stats);
}

template <class Array0, class Array1, class Array2, class Array3, class Array4, class Array5,
          class Array6, class Array7, class Array8, class Array9, class Array10, class Array11>
void stats_accumulate_lh_tend(ScratchWorkspace& workspace, // InOut
                              const grid& gr, int ngrdcol, int hydromet_dim,
                              const hm_metadata_type& hm_metadata,                       // In
                              Array0 lh_hydromet_mc, Array1 lh_Ncm_mc,                   // In
                              Array2 lh_thlm_mc, Array3 lh_rvm_mc, Array4 lh_rcm_mc,     // In
                              Array5 lh_AKm, Array6 AKm, Array7 AKstd, Array8 AKstd_cld, // In
                              Array9 lh_rcm_avg, Array10 AKm_rcm, Array11 AKm_rcc,       // In
                              stats_type& stats                                          // InOut
) {
  // Stage scientific inputs only for a Native sampled timestep.
  if (!stats.l_sample)
    return;
  static_assert(Array0::rank == 3);
  static_assert(Array1::rank == 2);
  static_assert(Array2::rank == 2);
  static_assert(Array3::rank == 2);
  static_assert(Array4::rank == 2);
  static_assert(Array5::rank == 2);
  static_assert(Array6::rank == 2);
  static_assert(Array7::rank == 2);
  static_assert(Array8::rank == 2);
  static_assert(Array9::rank == 2);
  static_assert(Array10::rank == 2);
  static_assert(Array11::rank == 2);
  auto workspace_scope = workspace.scope();
  stats_accumulate_lh_tend(
      workspace, gr, ngrdcol, hydromet_dim, hm_metadata, to_host_space(workspace, lh_hydromet_mc),
      to_host_space(workspace, lh_Ncm_mc), to_host_space(workspace, lh_thlm_mc),
      to_host_space(workspace, lh_rvm_mc), to_host_space(workspace, lh_rcm_mc),
      to_host_space(workspace, lh_AKm), to_host_space(workspace, AKm),
      to_host_space(workspace, AKstd), to_host_space(workspace, AKstd_cld),
      to_host_space(workspace, lh_rcm_avg), to_host_space(workspace, AKm_rcm),
      to_host_space(workspace, AKm_rcc), stats);
}

} // namespace clubb_kokkos
