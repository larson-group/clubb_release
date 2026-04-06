"""User-facing wrappers for routines from CLUBB_core/mixing_length.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.pdf_params import pdf_parameter
from clubb_python.derived_types.pdf_params_converter import set_fortran_pdf_params
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def set_lscale_max(ngrdcol: int, l_implemented: bool, host_dx, host_dy):
    """Compute maximum allowed Lscale from host horizontal grid spacing."""
    return clubb_f2py.f2py_set_lscale_max(
        l_implemented, f_arr(host_dx), f_arr(host_dy), ngrdcol=int(ngrdcol)
    )


def diagnose_lscale_from_tau(
    gr: Grid,
    nzm: int,
    nzt: int,
    ngrdcol: int,
    upwp_sfc, vpwp_sfc, ddzt_umvm_sqd,
    ice_supersat_frac, em, sqrt_em_zt,
    ufmin: float, tau_const: float,
    sfc_elevation, lscale_max, clubb_params,
    l_e3sm_config: bool,
    l_smooth_heaviside_tau_wpxp: bool,
    brunt_vaisala_freq_sqd_smth, ri_zm,
    err_info: ErrInfo,
):
    """Diagnose tau and Lscale profiles from shear/buoyancy inputs."""
    set_fortran_grid(gr)
    set_fortran_err_info(err_info)
    result = clubb_f2py.f2py_diagnose_lscale_from_tau(
        f_arr(upwp_sfc), f_arr(vpwp_sfc), f_arr(ddzt_umvm_sqd),
        f_arr(ice_supersat_frac), f_arr(em), f_arr(sqrt_em_zt),
        float(ufmin), float(tau_const),
        f_arr(sfc_elevation), f_arr(lscale_max), f_arr(clubb_params),
        l_e3sm_config, l_smooth_heaviside_tau_wpxp,
        f_arr(brunt_vaisala_freq_sqd_smth), f_arr(ri_zm),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
    return get_fortran_err_info(), *result


def calc_lscale_directly(
    gr: Grid,
    ngrdcol: int,
    nzm: int,
    nzt: int,
    l_implemented: bool,
    p_in_pa, exner, rtm, thlm, thvm, newmu,
    rtp2_zt, thlp2_zt, rtpthlp_zt, em, thv_ds_zt, lscale_max, lmin: float,
    clubb_params, saturation_formula: int,
    l_lscale_plume_centered: bool,
    pdf_params: pdf_parameter,
    err_info: ErrInfo,
):
    """Diagnose Lscale directly from thermodynamic profiles and PDF data."""
    set_fortran_grid(gr)
    set_fortran_pdf_params(pdf_params)
    set_fortran_err_info(err_info)
    result = clubb_f2py.f2py_calc_lscale_directly(
        l_implemented,
        f_arr(p_in_pa), f_arr(exner), f_arr(rtm), f_arr(thlm), f_arr(thvm), f_arr(newmu),
        f_arr(rtp2_zt), f_arr(thlp2_zt), f_arr(rtpthlp_zt), f_arr(em), f_arr(thv_ds_zt),
        f_arr(lscale_max), float(lmin), f_arr(clubb_params), int(saturation_formula),
        l_lscale_plume_centered,
        ngrdcol=int(ngrdcol), nzm=int(nzm), nzt=int(nzt),
    )
    return get_fortran_err_info(), *result
