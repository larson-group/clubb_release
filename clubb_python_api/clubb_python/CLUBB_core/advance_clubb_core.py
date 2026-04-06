"""User-facing wrappers for routines from CLUBB_core/advance_clubb_core_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.sclr_idx import SclrIdx
from clubb_python.derived_types.sclr_idx_converter import set_fortran_sclr_idx
from clubb_python.derived_types.config_flags import ConfigFlags
from clubb_python.derived_types.config_flags_converter import set_fortran_config_flags
from clubb_python.derived_types.nu_vert_res_dep import NuVertResDep
from clubb_python.derived_types.nu_vert_res_dep_converter import set_fortran_nu_vert_res_dep
from clubb_python.derived_types.pdf_params import pdf_parameter
from clubb_python.derived_types.pdf_params_converter import (
    get_fortran_implicit_coefs,
    get_fortran_pdf_params,
    get_fortran_pdf_params_zm,
    set_fortran_pdf_params,
    set_fortran_pdf_params_zm,
)
from clubb_python.derived_types.pdf_params import implicit_coefs_terms
from clubb_python.derived_types.pdf_params_converter import set_fortran_implicit_coefs
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def advance_clubb_core(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, l_implemented: bool, dt: float,
    fcor, fcor_y, sfc_elevation,
    hydromet_dim: int, sclr_dim: int, edsclr_dim: int,
    sclr_tol,
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing,
    sclrm_forcing, edsclrm_forcing,
    wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, rtpthlp_forcing, wm_zm, wm_zt,
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, p_sfc,
    wpsclrp_sfc, wpedsclrp_sfc,
    upwp_sfc_pert, vpwp_sfc_pert,
    rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg,
    p_in_Pa, rho_zm, rho, exner,
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt,
    l_mix_rat_hm,
    rfrzm,
    wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt,
    host_dx, host_dy,
    clubb_params, lmin: float, mixt_frac_max_mag: float,
    t0_val: float, ts_nudge: float,
    rtm_min: float, rtm_nudge_max_altitude: float,
    um, vm, upwp, vpwp, up2, vp2, up3, vp3,
    thlm, rtm, wprtp, wpthlp,
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp,
    sclrm,
    sclrp2, sclrp3, sclrprtp, sclrpthlp,
    wpsclrp, edsclrm,
    rcm, cloud_frac,
    wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp,
    sclrpthvp,
    wp2rtp, wp2thlp, uprcp, vprcp, rc_coef_zm, wp4,
    wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac,
    um_pert, vm_pert, upwp_pert, vpwp_pert,
    sclr_idx: SclrIdx,
    config_flags: ConfigFlags,
    nu_vert_res_dep: NuVertResDep,
    pdf_params: pdf_parameter,
    pdf_params_zm: pdf_parameter,
    pdf_implicit_coefs_terms: implicit_coefs_terms,
    err_info: ErrInfo,
):
    """Advance CLUBB core one timestep.

    Hidden Fortran-derived-type data is passed through `derived_type_storage`,
    but the Python return follows the compressed source Fortran
    `intent(inout/out)` order: numeric `inout`, refreshed hidden UDTs, then
    numeric `out`.
    """
    set_fortran_grid(gr)
    set_fortran_sclr_idx(sclr_idx)
    set_fortran_config_flags(config_flags)
    set_fortran_nu_vert_res_dep(nu_vert_res_dep)
    set_fortran_pdf_params(pdf_params)
    set_fortran_pdf_params_zm(pdf_params_zm)
    set_fortran_implicit_coefs(pdf_implicit_coefs_terms)
    set_fortran_err_info(err_info)

    result = clubb_f2py.f2py_advance_clubb_core(
        l_implemented=bool(l_implemented),
        dt=dt,
        fcor=f_arr(fcor),
        fcor_y=f_arr(fcor_y),
        sfc_elevation=f_arr(sfc_elevation),
        hydromet_dim=int(hydromet_dim),
        edsclr_dim=int(edsclr_dim),
        sclr_tol=f_arr(sclr_tol),
        thlm_forcing=f_arr(thlm_forcing),
        rtm_forcing=f_arr(rtm_forcing),
        um_forcing=f_arr(um_forcing),
        vm_forcing=f_arr(vm_forcing),
        sclrm_forcing=f_arr(sclrm_forcing),
        edsclrm_forcing=f_arr(edsclrm_forcing),
        wprtp_forcing=f_arr(wprtp_forcing),
        wpthlp_forcing=f_arr(wpthlp_forcing),
        rtp2_forcing=f_arr(rtp2_forcing),
        thlp2_forcing=f_arr(thlp2_forcing),
        rtpthlp_forcing=f_arr(rtpthlp_forcing),
        wm_zm=f_arr(wm_zm),
        wm_zt=f_arr(wm_zt),
        wpthlp_sfc=f_arr(wpthlp_sfc),
        wprtp_sfc=f_arr(wprtp_sfc),
        upwp_sfc=f_arr(upwp_sfc),
        vpwp_sfc=f_arr(vpwp_sfc),
        p_sfc=f_arr(p_sfc),
        wpsclrp_sfc=f_arr(wpsclrp_sfc),
        wpedsclrp_sfc=f_arr(wpedsclrp_sfc),
        upwp_sfc_pert=f_arr(upwp_sfc_pert),
        vpwp_sfc_pert=f_arr(vpwp_sfc_pert),
        rtm_ref=f_arr(rtm_ref),
        thlm_ref=f_arr(thlm_ref),
        um_ref=f_arr(um_ref),
        vm_ref=f_arr(vm_ref),
        ug=f_arr(ug),
        vg=f_arr(vg),
        rho_zm=f_arr(rho_zm),
        rho=f_arr(rho),
        rho_ds_zm=f_arr(rho_ds_zm),
        rho_ds_zt=f_arr(rho_ds_zt),
        invrs_rho_ds_zm=f_arr(invrs_rho_ds_zm),
        invrs_rho_ds_zt=f_arr(invrs_rho_ds_zt),
        thv_ds_zm=f_arr(thv_ds_zm),
        thv_ds_zt=f_arr(thv_ds_zt),
        l_mix_rat_hm=np.asfortranarray(l_mix_rat_hm, dtype=np.bool_),
        rfrzm=f_arr(rfrzm),
        wphydrometp=f_arr(wphydrometp),
        wp2hmp=f_arr(wp2hmp),
        rtphmp_zt=f_arr(rtphmp_zt),
        thlphmp_zt=f_arr(thlphmp_zt),
        host_dx=f_arr(host_dx),
        host_dy=f_arr(host_dy),
        clubb_params=f_arr(clubb_params),
        lmin=lmin,
        mixt_frac_max_mag=mixt_frac_max_mag,
        t0_val=t0_val,
        ts_nudge=ts_nudge,
        rtm_min=rtm_min,
        rtm_nudge_max_altitude=rtm_nudge_max_altitude,
        um=f_arr(um),
        vm=f_arr(vm),
        up3=f_arr(up3),
        vp3=f_arr(vp3),
        thlm=f_arr(thlm),
        rtm=f_arr(rtm),
        rtp3=f_arr(rtp3),
        thlp3=f_arr(thlp3),
        wp3=f_arr(wp3),
        upwp=f_arr(upwp),
        vpwp=f_arr(vpwp),
        up2=f_arr(up2),
        vp2=f_arr(vp2),
        wprtp=f_arr(wprtp),
        wpthlp=f_arr(wpthlp),
        rtp2=f_arr(rtp2),
        thlp2=f_arr(thlp2),
        rtpthlp=f_arr(rtpthlp),
        wp2=f_arr(wp2),
        sclrm=f_arr(sclrm),
        sclrp3=f_arr(sclrp3),
        wpsclrp=f_arr(wpsclrp),
        sclrp2=f_arr(sclrp2),
        sclrprtp=f_arr(sclrprtp),
        sclrpthlp=f_arr(sclrpthlp),
        p_in_pa=f_arr(p_in_Pa),
        exner=f_arr(exner),
        rcm=f_arr(rcm),
        cloud_frac=f_arr(cloud_frac),
        wp2thvp=f_arr(wp2thvp),
        wp2up=f_arr(wp2up),
        wpthvp=f_arr(wpthvp),
        rtpthvp=f_arr(rtpthvp),
        thlpthvp=f_arr(thlpthvp),
        sclrpthvp=f_arr(sclrpthvp),
        wp2rtp=f_arr(wp2rtp),
        wp2thlp=f_arr(wp2thlp),
        wpup2=f_arr(wpup2),
        wpvp2=f_arr(wpvp2),
        ice_supersat_frac=f_arr(ice_supersat_frac),
        uprcp=f_arr(uprcp),
        vprcp=f_arr(vprcp),
        rc_coef_zm=f_arr(rc_coef_zm),
        wp4=f_arr(wp4),
        wp2up2=f_arr(wp2up2),
        wp2vp2=f_arr(wp2vp2),
        um_pert=f_arr(um_pert),
        vm_pert=f_arr(vm_pert),
        upwp_pert=f_arr(upwp_pert),
        vpwp_pert=f_arr(vpwp_pert),
        edsclrm=f_arr(edsclrm),
        nzm=int(nzm),
        nzt=int(nzt),
        ngrdcol=int(ngrdcol),
        sclr_dim=int(sclr_dim),
    )
    return (
        *result[:51],
        get_fortran_pdf_params(),
        get_fortran_pdf_params_zm(),
        get_fortran_implicit_coefs(),
        get_fortran_err_info(),
        *result[51:],
    )
