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


def _parameterization_profile_arr(value, ngrdcol: int):
    arr = np.asarray(value)
    if arr.ndim == 1:
        arr = arr.reshape((1, arr.shape[0]))
    elif arr.ndim == 2 and arr.shape[0] != int(ngrdcol):
        arr = arr.reshape((1, arr.shape[0], arr.shape[1]))
    return f_arr(arr)


def _parameterization_sfc_arr(value):
    arr = np.asarray(value)
    if arr.ndim == 0:
        arr = arr.reshape((1,))
    return f_arr(arr)


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
    nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, edsclr_dim: int,
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
    profile_arr = lambda value: _parameterization_profile_arr(value, ngrdcol)
    sfc_arr = _parameterization_sfc_arr

    clubb_f2py.f2py_parameterization_check(
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
        sclr_dim=int(sclr_dim), edsclr_dim=int(edsclr_dim),
        thlm_forcing=profile_arr(thlm_forcing), rtm_forcing=profile_arr(rtm_forcing),
        um_forcing=profile_arr(um_forcing), vm_forcing=profile_arr(vm_forcing),
        wm_zm=profile_arr(wm_zm), wm_zt=profile_arr(wm_zt), p_in_pa=profile_arr(p_in_pa),
        rho_zm=profile_arr(rho_zm), rho=profile_arr(rho), exner=profile_arr(exner),
        rho_ds_zm=profile_arr(rho_ds_zm), rho_ds_zt=profile_arr(rho_ds_zt),
        invrs_rho_ds_zm=profile_arr(invrs_rho_ds_zm), invrs_rho_ds_zt=profile_arr(invrs_rho_ds_zt),
        thv_ds_zm=profile_arr(thv_ds_zm), thv_ds_zt=profile_arr(thv_ds_zt),
        wpthlp_sfc=sfc_arr(wpthlp_sfc), wprtp_sfc=sfc_arr(wprtp_sfc),
        upwp_sfc=sfc_arr(upwp_sfc), vpwp_sfc=sfc_arr(vpwp_sfc), p_sfc=sfc_arr(p_sfc),
        um=profile_arr(um), upwp=profile_arr(upwp), vm=profile_arr(vm), vpwp=profile_arr(vpwp),
        up2=profile_arr(up2), vp2=profile_arr(vp2),
        rtm=profile_arr(rtm), wprtp=profile_arr(wprtp), thlm=profile_arr(thlm), wpthlp=profile_arr(wpthlp),
        wp2=profile_arr(wp2), wp3=profile_arr(wp3), rtp2=profile_arr(rtp2),
        thlp2=profile_arr(thlp2), rtpthlp=profile_arr(rtpthlp),
        prefix=str(prefix),
        wpsclrp_sfc=profile_arr(wpsclrp_sfc), wpedsclrp_sfc=profile_arr(wpedsclrp_sfc),
        sclrm=profile_arr(sclrm), wpsclrp=profile_arr(wpsclrp), sclrp2=profile_arr(sclrp2),
        sclrprtp=profile_arr(sclrprtp), sclrpthlp=profile_arr(sclrpthlp),
        sclrm_forcing=profile_arr(sclrm_forcing),
        edsclrm=profile_arr(edsclrm), edsclrm_forcing=profile_arr(edsclrm_forcing),
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
