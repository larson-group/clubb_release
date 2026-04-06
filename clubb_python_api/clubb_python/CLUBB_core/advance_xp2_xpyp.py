"""User-facing wrappers for routines from CLUBB_core/advance_xp2_xpyp_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.sclr_idx import SclrIdx
from clubb_python.derived_types.sclr_idx_converter import set_fortran_sclr_idx
from clubb_python.derived_types.nu_vert_res_dep import NuVertResDep
from clubb_python.derived_types.nu_vert_res_dep_converter import set_fortran_nu_vert_res_dep
from clubb_python.derived_types.pdf_params import implicit_coefs_terms, pdf_parameter
from clubb_python.derived_types.pdf_params_converter import set_fortran_implicit_coefs, set_fortran_pdf_params
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def advance_xp2_xpyp(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, sclr_tol,
    invrs_tau_xp2_zm, invrs_tau_c4_zm, invrs_tau_c14_zm, wm_zm,
    rtm, wprtp, thlm, wpthlp, wpthvp, um, vm,
    wp2, wp2_zt, wp3, upwp, vpwp,
    sigma_sqd_w, wprtp2, wpthlp2, wprtpthlp, kh_zt,
    rtp2_forcing, thlp2_forcing, rtpthlp_forcing,
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, thv_ds_zm, cloud_frac,
    wp3_on_wp2, wp3_on_wp2_zt, dt: float, fcor_y,
    sclrm, wpsclrp, wpsclrp2, wpsclrprtp, wpsclrpthlp, lhs_splat_wp2,
    clubb_params, iipdf_type: int, tridiag_solve_method: int, fill_holes_type: int,
    l_predict_upwp_vpwp: bool, l_ho_nontrad_coriolis: bool, l_min_xp2_from_corr_wx: bool,
    l_c2_cloud_frac: bool, l_upwind_xpyp_ta: bool, l_godunov_upwind_xpyp_ta: bool,
    l_lmm_stepping: bool,
    rtp2, thlp2, rtpthlp, up2, vp2, sclrp2, sclrprtp, sclrpthlp,
    sclr_idx: SclrIdx,
    nu_vert_res_dep: NuVertResDep,
    pdf_implicit_coefs_terms: implicit_coefs_terms,
    err_info: ErrInfo,
):
    """Advance scalar variances/covariances and horizontal turbulence components.

    Pushes required Python UDT mirrors into `derived_type_storage`
    before the Fortran call.
    """
    set_fortran_grid(gr)
    set_fortran_sclr_idx(sclr_idx)
    set_fortran_nu_vert_res_dep(nu_vert_res_dep)
    set_fortran_implicit_coefs(pdf_implicit_coefs_terms)
    set_fortran_err_info(err_info)

    result = clubb_f2py.f2py_advance_xp2_xpyp(
        int(sclr_dim), f_arr(sclr_tol),
        f_arr(invrs_tau_xp2_zm), f_arr(invrs_tau_c4_zm), f_arr(invrs_tau_c14_zm), f_arr(wm_zm),
        f_arr(rtm), f_arr(wprtp), f_arr(thlm), f_arr(wpthlp), f_arr(wpthvp), f_arr(um), f_arr(vm),
        f_arr(wp2), f_arr(wp2_zt), f_arr(wp3), f_arr(upwp), f_arr(vpwp),
        f_arr(sigma_sqd_w), f_arr(wprtp2), f_arr(wpthlp2), f_arr(wprtpthlp), f_arr(kh_zt),
        f_arr(rtp2_forcing), f_arr(thlp2_forcing), f_arr(rtpthlp_forcing),
        f_arr(rho_ds_zm), f_arr(rho_ds_zt), f_arr(invrs_rho_ds_zm), f_arr(thv_ds_zm), f_arr(cloud_frac),
        f_arr(wp3_on_wp2), f_arr(wp3_on_wp2_zt), float(dt), f_arr(fcor_y),
        f_arr(sclrm), f_arr(wpsclrp), f_arr(wpsclrp2), f_arr(wpsclrprtp), f_arr(wpsclrpthlp), f_arr(lhs_splat_wp2),
        f_arr(clubb_params), int(iipdf_type), int(tridiag_solve_method), int(fill_holes_type),
        l_predict_upwp_vpwp, l_ho_nontrad_coriolis, l_min_xp2_from_corr_wx,
        l_c2_cloud_frac, l_upwind_xpyp_ta, l_godunov_upwind_xpyp_ta,
        l_lmm_stepping,
        f_arr(rtp2), f_arr(thlp2), f_arr(rtpthlp), f_arr(up2), f_arr(vp2),
        f_arr(sclrp2), f_arr(sclrprtp), f_arr(sclrpthlp),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
    return *result, get_fortran_err_info()


def update_xp2_mc(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, dt: float,
    cloud_frac, rcm, rvm, thlm, wm, exner, rrm_evap,
    pdf_params: pdf_parameter,
    rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc, rtpthlp_mc,
):
    """Update evaporation-driven covariance tendencies using grid/pdf data in module storage."""
    set_fortran_grid(gr)
    set_fortran_pdf_params(pdf_params)
    return clubb_f2py.f2py_update_xp2_mc(
        float(dt), f_arr(cloud_frac), f_arr(rcm), f_arr(rvm), f_arr(thlm), f_arr(wm), f_arr(exner), f_arr(rrm_evap),
        f_arr(rtp2_mc), f_arr(thlp2_mc), f_arr(wprtp_mc), f_arr(wpthlp_mc), f_arr(rtpthlp_mc),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
