"""User-facing wrappers for routines from CLUBB_core/numerical_check.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.config_flags import ConfigFlags
from clubb_python.derived_types.config_flags_converter import set_fortran_config_flags
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info
from clubb_python.derived_types.pdf_params import pdf_parameter
from clubb_python.derived_types.pdf_params_converter import set_fortran_pdf_params


def calculate_spurious_source(
    integral_after: float, integral_before: float,
    flux_top: float, flux_sfc: float, integral_forcing: float, dt: float,
):
    """Compute column conservation residual (negative means spurious sink)."""
    return float(clubb_f2py.f2py_calculate_spurious_source(
        float(integral_after), float(integral_before), float(flux_top),
        float(flux_sfc), float(integral_forcing), float(dt)))


def check_clubb_settings(
    ngrdcol: int,
    params,
    config_flags: ConfigFlags,
    err_info: ErrInfo,
    l_implemented: bool,
    l_input_fields: bool,
):
    """Validate CLUBB configuration settings."""
    set_fortran_config_flags(config_flags)
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_check_clubb_settings(
        f_arr(params), l_implemented, l_input_fields, ngrdcol=int(ngrdcol))
    return get_fortran_err_info()


def sfc_varnce_check(
    sclr_dim: int, wp2_sfc: float, up2_sfc: float, vp2_sfc: float, thlp2_sfc: float,
    rtp2_sfc: float, rtpthlp_sfc: float, sclrp2_sfc, sclrprtp_sfc, sclrpthlp_sfc,
    err_info: ErrInfo,
):
    """Run numerical surface-variance NaN checks against stored err_info."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_sfc_varnce_check(
        int(sclr_dim), float(wp2_sfc), float(up2_sfc), float(vp2_sfc), float(thlp2_sfc),
        float(rtp2_sfc), float(rtpthlp_sfc), f_arr(sclrp2_sfc), f_arr(sclrprtp_sfc), f_arr(sclrpthlp_sfc),
    )
    return get_fortran_err_info()


def parameterization_check(
    nzm: int, nzt: int, sclr_dim: int, edsclr_dim: int,
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, wm_zm, wm_zt, p_in_pa,
    rho_zm, rho, exner, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
    thv_ds_zm, thv_ds_zt, wpthlp_sfc: float, wprtp_sfc: float, upwp_sfc: float,
    vpwp_sfc: float, p_sfc: float, um, upwp, vm, vpwp, up2, vp2, rtm, wprtp,
    thlm, wpthlp, wp2, wp3, rtp2, thlp2, rtpthlp, prefix: str, wpsclrp_sfc,
    wpedsclrp_sfc, sclrm, wpsclrp, sclrp2, sclrprtp, sclrpthlp, sclrm_forcing,
    edsclrm, edsclrm_forcing, err_info: ErrInfo,
):
    """Run numerical parameterization checks."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_parameterization_check(
        nzm=int(nzm), nzt=int(nzt), sclr_dim=int(sclr_dim), edsclr_dim=int(edsclr_dim),
        thlm_forcing=f_arr(thlm_forcing), rtm_forcing=f_arr(rtm_forcing),
        um_forcing=f_arr(um_forcing), vm_forcing=f_arr(vm_forcing),
        wm_zm=f_arr(wm_zm), wm_zt=f_arr(wm_zt), p_in_pa=f_arr(p_in_pa),
        rho_zm=f_arr(rho_zm), rho=f_arr(rho), exner=f_arr(exner),
        rho_ds_zm=f_arr(rho_ds_zm), rho_ds_zt=f_arr(rho_ds_zt),
        invrs_rho_ds_zm=f_arr(invrs_rho_ds_zm), invrs_rho_ds_zt=f_arr(invrs_rho_ds_zt),
        thv_ds_zm=f_arr(thv_ds_zm), thv_ds_zt=f_arr(thv_ds_zt),
        wpthlp_sfc=float(wpthlp_sfc), wprtp_sfc=float(wprtp_sfc),
        upwp_sfc=float(upwp_sfc), vpwp_sfc=float(vpwp_sfc), p_sfc=float(p_sfc),
        um=f_arr(um), upwp=f_arr(upwp), vm=f_arr(vm), vpwp=f_arr(vpwp), up2=f_arr(up2), vp2=f_arr(vp2),
        rtm=f_arr(rtm), wprtp=f_arr(wprtp), thlm=f_arr(thlm), wpthlp=f_arr(wpthlp),
        wp2=f_arr(wp2), wp3=f_arr(wp3), rtp2=f_arr(rtp2), thlp2=f_arr(thlp2), rtpthlp=f_arr(rtpthlp),
        prefix=str(prefix),
        wpsclrp_sfc=f_arr(wpsclrp_sfc), wpedsclrp_sfc=f_arr(wpedsclrp_sfc),
        sclrm=f_arr(sclrm), wpsclrp=f_arr(wpsclrp), sclrp2=f_arr(sclrp2), sclrprtp=f_arr(sclrprtp),
        sclrpthlp=f_arr(sclrpthlp), sclrm_forcing=f_arr(sclrm_forcing),
        edsclrm=f_arr(edsclrm), edsclrm_forcing=f_arr(edsclrm_forcing),
    )
    return get_fortran_err_info()


def pdf_closure_check(
    nz: int, sclr_dim: int,
    wp4, wprtp2, wp2rtp, wpthlp2, wp2thlp, cloud_frac, rcm, wpthvp, wp2thvp, wp2up,
    rtpthvp, thlpthvp, wprcp, wp2rcp, rtprcp, thlprcp, rcp2, wprtpthlp, crt_1, crt_2,
    cthl_1, cthl_2, sclrpthvp, sclrprcp, wpsclrp2, wpsclrprtp, wpsclrpthlp, wp2sclrp,
    pdf_params: pdf_parameter,
    err_info: ErrInfo,
):
    """Run pdf_closure_check using pdf_params/stats/err_info data in module storage."""
    set_fortran_pdf_params(pdf_params)
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_pdf_closure_check(
        int(sclr_dim), f_arr(wp4), f_arr(wprtp2), f_arr(wp2rtp), f_arr(wpthlp2), f_arr(wp2thlp), f_arr(cloud_frac),
        f_arr(rcm), f_arr(wpthvp), f_arr(wp2thvp), f_arr(wp2up), f_arr(rtpthvp), f_arr(thlpthvp),
        f_arr(wprcp), f_arr(wp2rcp), f_arr(rtprcp), f_arr(thlprcp), f_arr(rcp2), f_arr(wprtpthlp),
        f_arr(crt_1), f_arr(crt_2), f_arr(cthl_1), f_arr(cthl_2), f_arr(sclrpthvp), f_arr(sclrprcp),
        f_arr(wpsclrp2), f_arr(wpsclrprtp), f_arr(wpsclrpthlp), f_arr(wp2sclrp),
        nz=int(nz))
    return get_fortran_err_info()


def length_check(nzt: int, lscale, lscale_up, lscale_down, err_info: ErrInfo):
    """Run the mixing-length NaN check against err_info data in module storage."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_length_check(f_arr(lscale), f_arr(lscale_up), f_arr(lscale_down), nzt=int(nzt))
    return get_fortran_err_info()
