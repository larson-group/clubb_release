"""User-facing wrappers for routines from CLUBB_core/stats_clubb_utilities.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.pdf_params import pdf_parameter
from clubb_python.derived_types.pdf_params_converter import (
    set_fortran_pdf_params,
    set_fortran_pdf_params_zm,
)


def stats_accumulate(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, edsclr_dim: int,
    invrs_dzm, zt, dzm, dzt, dt: float,
    k_lb_zm: int, k_ub_zm: int, l_implemented: bool, l_host_applies_sfc_fluxes: bool,
    l_stability_correct_tau_zm: bool, clubb_params,
    um, vm, upwp, vpwp, up2, vp2,
    thlm, rtm, thlm_before, rtm_before, thlm_forcing, rtm_forcing,
    wpthlp_sfc, wprtp_sfc, wprtp, wpthlp,
    wp2, wp3, rtp2, rtp3, thlp2, thlp3,
    rtpthlp,
    wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp,
    p_in_pa, exner, rho, rho_zm,
    rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt,
    wm_zt, wm_zm, rcm, wprcp, rc_coef,
    rc_coef_zm,
    rcm_zm, rtm_zm, thlm_zm, cloud_frac,
    ice_supersat_frac,
    cloud_frac_zm, ice_supersat_frac_zm, rcm_in_layer,
    cloud_cover, rcm_supersat_adj, sigma_sqd_w,
    thvm, ug, vg, lscale, wpthlp2, wp2thlp,
    wprtp2, wp2rtp,
    lscale_up, lscale_down, kh_zt, wp2rcp,
    wprtpthlp, rsat, wpup2, wpvp2,
    wp2up2, wp2vp2, wp4,
    tau_zm, kh_zm, thlprcp,
    rtprcp, rcp2, em, wp3_on_wp2, wp3_on_wp2_zt, skw_velocity,
    ddzt_umvm_sqd, stability_correction, a3_coef,
    w_up_in_cloud, w_down_in_cloud,
    cloudy_updraft_frac, cloudy_downdraft_frac,
    pdf_params: pdf_parameter,
    pdf_params_zm: pdf_parameter,
    sclrm, sclrp2,
    sclrprtp, sclrpthlp, sclrm_forcing, sclrpthvp,
    wpsclrp, sclrprcp, wp2sclrp, wpsclrp2,
    wpsclrprtp,
    wpsclrpthlp, wpedsclrp, edsclrm,
    edsclrm_forcing,
    saturation_formula: int,
    l_call_pdf_closure_twice: bool,
):
    """Run stats_accumulate with strict direct argument mapping to the wrapper."""
    set_fortran_grid(gr)
    set_fortran_pdf_params(pdf_params)
    set_fortran_pdf_params_zm(pdf_params_zm)

    clubb_f2py.f2py_stats_accumulate(
        int(sclr_dim), int(edsclr_dim),
        f_arr(invrs_dzm), f_arr(zt), f_arr(dzm), f_arr(dzt), float(dt),
        int(k_lb_zm), int(k_ub_zm), bool(l_implemented), bool(l_host_applies_sfc_fluxes),
        bool(l_stability_correct_tau_zm), f_arr(clubb_params),
        f_arr(um), f_arr(vm), f_arr(upwp), f_arr(vpwp), f_arr(up2), f_arr(vp2),
        f_arr(thlm), f_arr(rtm), f_arr(thlm_before), f_arr(rtm_before),
        f_arr(thlm_forcing), f_arr(rtm_forcing),
        f_arr(wpthlp_sfc), f_arr(wprtp_sfc), f_arr(wprtp), f_arr(wpthlp),
        f_arr(wp2), f_arr(wp3), f_arr(rtp2), f_arr(rtp3), f_arr(thlp2), f_arr(thlp3),
        f_arr(rtpthlp),
        f_arr(wpthvp), f_arr(wp2thvp), f_arr(wp2up), f_arr(rtpthvp), f_arr(thlpthvp),
        f_arr(p_in_pa), f_arr(exner), f_arr(rho), f_arr(rho_zm),
        f_arr(rho_ds_zm), f_arr(rho_ds_zt), f_arr(thv_ds_zm), f_arr(thv_ds_zt),
        f_arr(wm_zt), f_arr(wm_zm), f_arr(rcm), f_arr(wprcp), f_arr(rc_coef),
        f_arr(rc_coef_zm),
        f_arr(rcm_zm), f_arr(rtm_zm), f_arr(thlm_zm), f_arr(cloud_frac),
        f_arr(ice_supersat_frac),
        f_arr(cloud_frac_zm), f_arr(ice_supersat_frac_zm), f_arr(rcm_in_layer),
        f_arr(cloud_cover), f_arr(rcm_supersat_adj), f_arr(sigma_sqd_w),
        f_arr(thvm), f_arr(ug), f_arr(vg), f_arr(lscale), f_arr(wpthlp2), f_arr(wp2thlp),
        f_arr(wprtp2), f_arr(wp2rtp),
        f_arr(lscale_up), f_arr(lscale_down), f_arr(kh_zt), f_arr(wp2rcp),
        f_arr(wprtpthlp), f_arr(rsat), f_arr(wpup2), f_arr(wpvp2),
        f_arr(wp2up2), f_arr(wp2vp2), f_arr(wp4),
        f_arr(tau_zm), f_arr(kh_zm), f_arr(thlprcp),
        f_arr(rtprcp), f_arr(rcp2), f_arr(em), f_arr(wp3_on_wp2),
        f_arr(wp3_on_wp2_zt), f_arr(skw_velocity),
        f_arr(ddzt_umvm_sqd), f_arr(stability_correction), f_arr(a3_coef),
        f_arr(w_up_in_cloud), f_arr(w_down_in_cloud),
        f_arr(cloudy_updraft_frac), f_arr(cloudy_downdraft_frac),
        f_arr(sclrm), f_arr(sclrp2),
        f_arr(sclrprtp), f_arr(sclrpthlp), f_arr(sclrm_forcing), f_arr(sclrpthvp),
        f_arr(wpsclrp), f_arr(sclrprcp), f_arr(wp2sclrp), f_arr(wpsclrp2),
        f_arr(wpsclrprtp),
        f_arr(wpsclrpthlp), f_arr(wpedsclrp), f_arr(edsclrm),
        f_arr(edsclrm_forcing),
        int(saturation_formula),
        bool(l_call_pdf_closure_twice),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
