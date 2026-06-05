"""User-facing wrappers for routines from CLUBB_core/stats_clubb_utilities.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def stats_accumulate(
    nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, edsclr_dim: int, gr: Grid,
    dt: float,
    l_implemented: bool, l_host_applies_sfc_fluxes: bool,
    l_stability_correct_tau_zm: bool, clubb_params,
    um, vm, upwp, vpwp, up2, vp2,
    thlm, rtm, thlm_before, rtm_before, thlm_forcing, rtm_forcing,
    wpthlp_sfc, wprtp_sfc, wprtp, wpthlp,
    wp2, wp3, rtp2, rtp3, thlp2, thlp3,
    rtpthlp,
    p_in_pa, exner, rho, rho_zm,
    rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt,
    wm_zt, wm_zm, rcm, cloud_frac,
    thvm, ug, vg, ddzt_umvm_sqd, stability_correction, kh_zt, rsat, kh_zm, em,
    sclrm, sclrp2,
    sclrprtp, sclrpthlp, sclrm_forcing,
    wpsclrp, wpedsclrp, edsclrm,
    edsclrm_forcing,
    saturation_formula: int,
):
    """Run stats_accumulate with strict direct argument mapping to Fortran."""
    set_fortran_grid(gr)

    clubb_f2py.f2py_stats_accumulate(
        int(sclr_dim), int(edsclr_dim), float(dt),
        bool(l_implemented), bool(l_host_applies_sfc_fluxes),
        bool(l_stability_correct_tau_zm), f_arr(clubb_params),
        f_arr(um), f_arr(vm), f_arr(upwp), f_arr(vpwp), f_arr(up2), f_arr(vp2),
        f_arr(thlm), f_arr(rtm), f_arr(thlm_before), f_arr(rtm_before),
        f_arr(thlm_forcing), f_arr(rtm_forcing),
        f_arr(wpthlp_sfc), f_arr(wprtp_sfc), f_arr(wprtp), f_arr(wpthlp),
        f_arr(wp2), f_arr(wp3), f_arr(rtp2), f_arr(rtp3), f_arr(thlp2), f_arr(thlp3),
        f_arr(rtpthlp),
        f_arr(p_in_pa), f_arr(exner), f_arr(rho), f_arr(rho_zm),
        f_arr(rho_ds_zm), f_arr(rho_ds_zt), f_arr(thv_ds_zm), f_arr(thv_ds_zt),
        f_arr(wm_zt), f_arr(wm_zm), f_arr(rcm), f_arr(cloud_frac),
        f_arr(thvm), f_arr(ug), f_arr(vg),
        f_arr(ddzt_umvm_sqd), f_arr(stability_correction),
        f_arr(kh_zt), f_arr(rsat), f_arr(kh_zm), f_arr(em),
        f_arr(sclrm), f_arr(sclrp2),
        f_arr(sclrprtp), f_arr(sclrpthlp), f_arr(sclrm_forcing),
        f_arr(wpsclrp), f_arr(wpedsclrp), f_arr(edsclrm),
        f_arr(edsclrm_forcing),
        int(saturation_formula),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
