// Verify the core's beginning-of-step fatal return through the host/device API.
// A NaN forcing must reach the core input diagnostic and return err_info before
// surface fluxes or any lower solver can change the prognostic fields.
#include "scratch_workspace_owner.hpp"
#include "CLUBB_core/error_code.hpp"
#include "CLUBB_core/clubb_api_module.hpp"
#include "CLUBB_core/parameter_indices.hpp"
#include <iostream>
#include <limits>
#include <stdexcept>

int main(int argc, char* argv[]) {
  // Two columns and zero scalar species exercise status conversion without
  // requiring initialized grid/PDF solvers: the diagnostic must return first.
  try {
    Kokkos::ScopeGuard kokkos(argc, argv);
    clubb_kokkos::test::ScratchWorkspaceOwner scratch_owner;
    auto& workspace = scratch_owner.workspace;
    using namespace clubb_kokkos;
    using HostReal1D  = Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    using HostReal2D  = Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    using HostReal3D  = Kokkos::View<real_t***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    using HostBool1D  = Kokkos::View<std::uint8_t*, Kokkos::LayoutLeft, Kokkos::HostSpace>;
    grid gr{};
    int nzm     = 3;
    int nzt     = 2;
    int ngrdcol = 2;
    bool l_implemented{};
    real_t dt = real_t(60);
    HostReal1D fcor("fcor", ngrdcol);
    HostReal1D fcor_y("fcor_y", ngrdcol);
    HostReal1D sfc_elevation("sfc_elevation", ngrdcol);
    int hydromet_dim = 0;
    int sclr_dim     = 0;
    HostReal1D sclr_tol("sclr_tol", sclr_dim);
    int edsclr_dim = 0;
    sclr_idx_type sclr_idx{};
    HostReal2D thlm_forcing("thlm_forcing", ngrdcol, nzt);
    HostReal2D rtm_forcing("rtm_forcing", ngrdcol, nzt);
    HostReal2D um_forcing("um_forcing", ngrdcol, nzt);
    HostReal2D vm_forcing("vm_forcing", ngrdcol, nzt);
    HostReal3D sclrm_forcing("sclrm_forcing", ngrdcol, nzt, sclr_dim);
    HostReal3D edsclrm_forcing("edsclrm_forcing", ngrdcol, nzt, edsclr_dim);
    HostReal2D wprtp_forcing("wprtp_forcing", ngrdcol, nzm);
    HostReal2D wpthlp_forcing("wpthlp_forcing", ngrdcol, nzm);
    HostReal2D rtp2_forcing("rtp2_forcing", ngrdcol, nzm);
    HostReal2D thlp2_forcing("thlp2_forcing", ngrdcol, nzm);
    HostReal2D rtpthlp_forcing("rtpthlp_forcing", ngrdcol, nzm);
    HostReal2D wm_zm("wm_zm", ngrdcol, nzm);
    HostReal2D wm_zt("wm_zt", ngrdcol, nzt);
    HostReal1D wpthlp_sfc("wpthlp_sfc", ngrdcol);
    HostReal1D wprtp_sfc("wprtp_sfc", ngrdcol);
    HostReal1D upwp_sfc("upwp_sfc", ngrdcol);
    HostReal1D vpwp_sfc("vpwp_sfc", ngrdcol);
    HostReal1D p_sfc("p_sfc", ngrdcol);
    HostReal2D wpsclrp_sfc("wpsclrp_sfc", ngrdcol, sclr_dim);
    HostReal2D wpedsclrp_sfc("wpedsclrp_sfc", ngrdcol, edsclr_dim);
    HostReal1D upwp_sfc_pert("upwp_sfc_pert", ngrdcol);
    HostReal1D vpwp_sfc_pert("vpwp_sfc_pert", ngrdcol);
    HostReal2D rtm_ref("rtm_ref", ngrdcol, nzt);
    HostReal2D thlm_ref("thlm_ref", ngrdcol, nzt);
    HostReal2D um_ref("um_ref", ngrdcol, nzt);
    HostReal2D vm_ref("vm_ref", ngrdcol, nzt);
    HostReal2D ug("ug", ngrdcol, nzt);
    HostReal2D vg("vg", ngrdcol, nzt);
    HostReal2D p_in_Pa("p_in_Pa", ngrdcol, nzt);
    HostReal2D rho_zm("rho_zm", ngrdcol, nzm);
    HostReal2D rho("rho", ngrdcol, nzt);
    HostReal2D exner("exner", ngrdcol, nzt);
    HostReal2D rho_ds_zm("rho_ds_zm", ngrdcol, nzm);
    HostReal2D rho_ds_zt("rho_ds_zt", ngrdcol, nzt);
    HostReal2D invrs_rho_ds_zm("invrs_rho_ds_zm", ngrdcol, nzm);
    HostReal2D invrs_rho_ds_zt("invrs_rho_ds_zt", ngrdcol, nzt);
    HostReal2D thv_ds_zm("thv_ds_zm", ngrdcol, nzm);
    HostReal2D thv_ds_zt("thv_ds_zt", ngrdcol, nzt);
    HostBool1D l_mix_rat_hm("l_mix_rat_hm", hydromet_dim);
    HostReal2D rfrzm("rfrzm", ngrdcol, nzt);
    HostReal3D wphydrometp("wphydrometp", ngrdcol, nzm, hydromet_dim);
    HostReal3D wp2hmp("wp2hmp", ngrdcol, nzt, hydromet_dim);
    HostReal3D rtphmp_zt("rtphmp_zt", ngrdcol, nzt, hydromet_dim);
    HostReal3D thlphmp_zt("thlphmp_zt", ngrdcol, nzt, hydromet_dim);
    HostReal1D host_dx("host_dx", ngrdcol);
    HostReal1D host_dy("host_dy", ngrdcol);
    HostReal2D clubb_params("clubb_params", ngrdcol, nparams);
    nu_vertical_res_dep nu_vert_res_dep{};
    HostReal1D lmin("lmin", ngrdcol);
    HostReal1D mixt_frac_max_mag("mixt_frac_max_mag", ngrdcol);
    real_t T0{};
    real_t ts_nudge{};
    real_t rtm_min{};
    real_t rtm_nudge_max_altitude{};
    clubb_config_flags_type clubb_config_flags{};
    clubb_kokkos::stats_type stats{};
    HostReal2D um("um", ngrdcol, nzt);
    HostReal2D vm("vm", ngrdcol, nzt);
    HostReal2D upwp("upwp", ngrdcol, nzm);
    HostReal2D vpwp("vpwp", ngrdcol, nzm);
    HostReal2D up2("up2", ngrdcol, nzm);
    HostReal2D vp2("vp2", ngrdcol, nzm);
    HostReal2D up3("up3", ngrdcol, nzt);
    HostReal2D vp3("vp3", ngrdcol, nzt);
    HostReal2D thlm("thlm", ngrdcol, nzt);
    HostReal2D rtm("rtm", ngrdcol, nzt);
    HostReal2D wprtp("wprtp", ngrdcol, nzm);
    HostReal2D wpthlp("wpthlp", ngrdcol, nzm);
    HostReal2D wp2("wp2", ngrdcol, nzm);
    HostReal2D wp3("wp3", ngrdcol, nzt);
    HostReal2D rtp2("rtp2", ngrdcol, nzm);
    HostReal2D rtp3("rtp3", ngrdcol, nzt);
    HostReal2D thlp2("thlp2", ngrdcol, nzm);
    HostReal2D thlp3("thlp3", ngrdcol, nzt);
    HostReal2D rtpthlp("rtpthlp", ngrdcol, nzm);
    HostReal3D sclrm("sclrm", ngrdcol, nzt, sclr_dim);
    HostReal3D sclrp2("sclrp2", ngrdcol, nzm, sclr_dim);
    HostReal3D sclrp3("sclrp3", ngrdcol, nzt, sclr_dim);
    HostReal3D sclrprtp("sclrprtp", ngrdcol, nzm, sclr_dim);
    HostReal3D sclrpthlp("sclrpthlp", ngrdcol, nzm, sclr_dim);
    HostReal3D wpsclrp("wpsclrp", ngrdcol, nzm, sclr_dim);
    HostReal3D edsclrm("edsclrm", ngrdcol, nzt, edsclr_dim);
    err_info_type err_info_api{};
    HostReal2D rcm("rcm", ngrdcol, nzt);
    HostReal2D cloud_frac("cloud_frac", ngrdcol, nzt);
    HostReal2D wpthvp("wpthvp", ngrdcol, nzm);
    HostReal2D wp2thvp("wp2thvp", ngrdcol, nzt);
    HostReal2D wp2up("wp2up", ngrdcol, nzt);
    HostReal2D rtpthvp("rtpthvp", ngrdcol, nzm);
    HostReal2D thlpthvp("thlpthvp", ngrdcol, nzm);
    HostReal3D sclrpthvp("sclrpthvp", ngrdcol, nzm, sclr_dim);
    HostReal2D wp2rtp("wp2rtp", ngrdcol, nzt);
    HostReal2D wp2thlp("wp2thlp", ngrdcol, nzt);
    HostReal2D uprcp("uprcp", ngrdcol, nzm);
    HostReal2D vprcp("vprcp", ngrdcol, nzm);
    HostReal2D rc_coef_zm("rc_coef_zm", ngrdcol, nzm);
    HostReal2D wp4("wp4", ngrdcol, nzm);
    HostReal2D wpup2("wpup2", ngrdcol, nzt);
    HostReal2D wpvp2("wpvp2", ngrdcol, nzt);
    HostReal2D wp2up2("wp2up2", ngrdcol, nzm);
    HostReal2D wp2vp2("wp2vp2", ngrdcol, nzm);
    HostReal2D ice_supersat_frac("ice_supersat_frac", ngrdcol, nzt);
    HostReal2D um_pert("um_pert", ngrdcol, nzt);
    HostReal2D vm_pert("vm_pert", ngrdcol, nzt);
    HostReal2D upwp_pert("upwp_pert", ngrdcol, nzm);
    HostReal2D vpwp_pert("vpwp_pert", ngrdcol, nzm);
    pdf_parameter pdf_params{};
    pdf_parameter pdf_params_zm{};
    implicit_coefs_terms pdf_implicit_coefs_terms{};
    HostReal2D Kh_zm("Kh_zm", ngrdcol, nzm);
    HostReal2D Kh_zt("Kh_zt", ngrdcol, nzt);
    HostReal2D thlprcp("thlprcp", ngrdcol, nzm);
    HostReal2D wprcp("wprcp", ngrdcol, nzm);
    HostReal2D w_up_in_cloud("w_up_in_cloud", ngrdcol, nzt);
    HostReal2D w_down_in_cloud("w_down_in_cloud", ngrdcol, nzt);
    HostReal2D cloudy_updraft_frac("cloudy_updraft_frac", ngrdcol, nzt);
    HostReal2D cloudy_downdraft_frac("cloudy_downdraft_frac", ngrdcol, nzt);
    HostReal2D rcm_in_layer("rcm_in_layer", ngrdcol, nzt);
    HostReal2D cloud_cover("cloud_cover", ngrdcol, nzt);
    HostReal2D invrs_tau_zm("invrs_tau_zm", ngrdcol, nzm);
    HostReal2D Lscale("Lscale", ngrdcol, nzt);

    set_clubb_debug_level_api(2);
    init_default_err_info_api(workspace, ngrdcol, err_info_api);
    // Retain the native fixed character width while identifying this check.
    err_info_api.err_header_global.replace(0, 17, "core_error_return");
    const auto original_header = err_info_api.err_header_global;
    Kokkos::deep_copy(um, real_t(23));
    Kokkos::deep_copy(wpthlp, real_t(17));
    thlm_forcing(0, 0) = std::numeric_limits<real_t>::quiet_NaN();

    advance_clubb_core_api(
        workspace, gr, nzm, nzt, ngrdcol, l_implemented, dt, fcor, fcor_y, sfc_elevation,
        hydromet_dim, sclr_dim, sclr_tol, edsclr_dim, sclr_idx, thlm_forcing, rtm_forcing,
        um_forcing, vm_forcing, sclrm_forcing, edsclrm_forcing, wprtp_forcing, wpthlp_forcing,
        rtp2_forcing, thlp2_forcing, rtpthlp_forcing, wm_zm, wm_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc,
        vpwp_sfc, p_sfc, wpsclrp_sfc, wpedsclrp_sfc, upwp_sfc_pert, vpwp_sfc_pert, rtm_ref,
        thlm_ref, um_ref, vm_ref, ug, vg, p_in_Pa, rho_zm, rho, exner, rho_ds_zm, rho_ds_zt,
        invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, l_mix_rat_hm, rfrzm, wphydrometp,
        wp2hmp, rtphmp_zt, thlphmp_zt, host_dx, host_dy, clubb_params, nu_vert_res_dep, lmin,
        mixt_frac_max_mag, T0, ts_nudge, rtm_min, rtm_nudge_max_altitude, clubb_config_flags, stats,
        um, vm, upwp, vpwp, up2, vp2, up3, vp3, thlm, rtm, wprtp, wpthlp, wp2, wp3, rtp2, rtp3,
        thlp2, thlp3, rtpthlp, sclrm, sclrp2, sclrp3, sclrprtp, sclrpthlp, wpsclrp, edsclrm,
        err_info_api, rcm, cloud_frac, wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, sclrpthvp, wp2rtp,
        wp2thlp, uprcp, vprcp, rc_coef_zm, wp4, wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac,
        um_pert, vm_pert, upwp_pert, vpwp_pert, pdf_params, pdf_params_zm, pdf_implicit_coefs_terms,
        Kh_zm, Kh_zt, thlprcp, wprcp, w_up_in_cloud, w_down_in_cloud, cloudy_updraft_frac,
        cloudy_downdraft_frac, rcm_in_layer, cloud_cover, invrs_tau_zm, Lscale);

    auto codes = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, err_info_api.err_code);
    for (int i = 0; i < ngrdcol; ++i) {
      if (codes(i) != clubb_fatal_error)
        throw std::runtime_error("diagnostic status did not return through the core API");
      for (int k = 0; k < nzt; ++k)
        if (um(i, k) != real_t(23))
          throw std::runtime_error("core continued after the initial fatal diagnostic");
      for (int k = 0; k < nzm; ++k)
        if (wpthlp(i, k) != real_t(17))
          throw std::runtime_error("surface fluxes changed after the initial fatal diagnostic");
    }
    if (err_info_api.err_header_global != original_header)
      throw std::runtime_error("core diagnostic changed the caller's header");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
