"""User-facing wrappers for routines from CLUBB_core/advance_wp2_wp3_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.nu_vert_res_dep import NuVertResDep
from clubb_python.derived_types.nu_vert_res_dep_converter import set_fortran_nu_vert_res_dep
from clubb_python.derived_types.pdf_params import implicit_coefs_terms
from clubb_python.derived_types.pdf_params_converter import set_fortran_implicit_coefs
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def advance_wp2_wp3(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, dt: float,
    sfc_elevation, fcor_y, sigma_sqd_w, wm_zm, wm_zt,
    a3_coef, a3_coef_zt, wp3_on_wp2,
    wpup2, wpvp2, wp2up2, wp2vp2, wp4,
    wpthvp, wp2thvp, wp2up, um, vm, upwp, vpwp,
    em, kh_zm, kh_zt, invrs_tau_c4_zm, invrs_tau_wp3_zt, invrs_tau_c1_zm,
    skw_zm, skw_zt, rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,
    thv_ds_zm, thv_ds_zt, mixt_frac, cx_fnc_richardson,
    lhs_splat_wp2, lhs_splat_wp3,
    wprtp, wpthlp, rtp2, thlp2, clubb_params,
    iipdf_type: int, penta_solve_method: int, fill_holes_type: int,
    l_min_wp2_from_corr_wx: bool, l_upwind_xm_ma: bool, l_tke_aniso: bool,
    l_standard_term_ta: bool, l_partial_upwind_wp3: bool, l_damp_wp2_using_em: bool,
    l_use_c11_richardson: bool, l_damp_wp3_skw_squared: bool, l_lmm_stepping: bool,
    l_use_tke_in_wp3_pr_turb_term: bool, l_use_tke_in_wp2_wp3_k_dfsn: bool,
    l_use_wp3_lim_with_smth_heaviside: bool, l_wp2_fill_holes_tke: bool,
    l_ho_nontrad_coriolis: bool,
    up2, vp2, wp2, wp3, wp3_zm, wp2_zt,
    nu_vert_res_dep: NuVertResDep,
    pdf_implicit_coefs_terms: implicit_coefs_terms,
    err_info: ErrInfo,
):
    """Advance w'^2 and w'^3 one model timestep.

    Pushes required Python UDT mirrors into `derived_type_storage`
    before the Fortran call.
    """
    set_fortran_grid(gr)
    set_fortran_nu_vert_res_dep(nu_vert_res_dep)
    set_fortran_implicit_coefs(pdf_implicit_coefs_terms)
    set_fortran_err_info(err_info)

    result = clubb_f2py.f2py_advance_wp2_wp3(
        float(dt),
        f_arr(sfc_elevation), f_arr(fcor_y), f_arr(sigma_sqd_w), f_arr(wm_zm), f_arr(wm_zt),
        f_arr(a3_coef), f_arr(a3_coef_zt), f_arr(wp3_on_wp2),
        f_arr(wpup2), f_arr(wpvp2), f_arr(wp2up2), f_arr(wp2vp2), f_arr(wp4),
        f_arr(wpthvp), f_arr(wp2thvp), f_arr(wp2up), f_arr(um), f_arr(vm), f_arr(upwp), f_arr(vpwp),
        f_arr(em), f_arr(kh_zm), f_arr(kh_zt), f_arr(invrs_tau_c4_zm), f_arr(invrs_tau_wp3_zt),
        f_arr(invrs_tau_c1_zm),
        f_arr(skw_zm), f_arr(skw_zt), f_arr(rho_ds_zm), f_arr(rho_ds_zt),
        f_arr(invrs_rho_ds_zm), f_arr(invrs_rho_ds_zt), f_arr(thv_ds_zm), f_arr(thv_ds_zt),
        f_arr(mixt_frac), f_arr(cx_fnc_richardson), f_arr(lhs_splat_wp2), f_arr(lhs_splat_wp3),
        f_arr(wprtp), f_arr(wpthlp), f_arr(rtp2), f_arr(thlp2), f_arr(clubb_params),
        int(iipdf_type), int(penta_solve_method), int(fill_holes_type),
        l_min_wp2_from_corr_wx, l_upwind_xm_ma, l_tke_aniso,
        l_standard_term_ta, l_partial_upwind_wp3, l_damp_wp2_using_em,
        l_use_c11_richardson, l_damp_wp3_skw_squared, l_lmm_stepping,
        l_use_tke_in_wp3_pr_turb_term, l_use_tke_in_wp2_wp3_k_dfsn,
        l_use_wp3_lim_with_smth_heaviside, l_wp2_fill_holes_tke,
        l_ho_nontrad_coriolis,
        f_arr(up2), f_arr(vp2), f_arr(wp2), f_arr(wp3), f_arr(wp3_zm), f_arr(wp2_zt),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
    return *result, get_fortran_err_info()
