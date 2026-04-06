"""User-facing wrappers for routines from CLUBB_core/advance_xm_wpxp_module.F90."""

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


def advance_xm_wpxp(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, sclr_tol, dt: float,
    sigma_sqd_w, wm_zm, wm_zt, wp2, lscale_zm,
    wp3_on_wp2, wp3_on_wp2_zt, kh_zt, kh_zm,
    stability_correction,
    invrs_tau_c6_zm, tau_max_zm, skw_zm, wp2rtp, rtpthvp,
    rtm_forcing, wprtp_forcing, rtm_ref, wp2thlp,
    thlpthvp, thlm_forcing, wpthlp_forcing, thlm_ref,
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2,
    w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm,
    l_implemented: bool, em, wp2sclrp, sclrpthvp, sclrm_forcing, sclrp2, cx_fnc_richardson,
    um_forcing, vm_forcing, ug, vg, wpthvp,
    fcor, fcor_y, um_ref, vm_ref, up2, vp2, uprcp, vprcp, rc_coef_zm,
    clubb_params, ts_nudge: float,
    iipdf_type: int, penta_solve_method: int, tridiag_solve_method: int, fill_holes_type: int,
    l_predict_upwp_vpwp: bool, l_ho_nontrad_coriolis: bool, l_ho_trad_coriolis: bool,
    l_diffuse_rtm_and_thlm: bool, l_stability_correct_kh_n2_zm: bool,
    l_godunov_upwind_wpxp_ta: bool, l_upwind_xm_ma: bool, l_uv_nudge: bool,
    l_tke_aniso: bool, l_diag_lscale_from_tau: bool, l_use_c7_richardson: bool,
    l_lmm_stepping: bool, l_enable_relaxed_clipping: bool, l_linearize_pbl_winds: bool,
    l_mono_flux_lim_thlm: bool, l_mono_flux_lim_rtm: bool, l_mono_flux_lim_um: bool,
    l_mono_flux_lim_vm: bool, l_mono_flux_lim_spikefix: bool,
    order_xm_wpxp: int, order_xp2_xpyp: int, order_wp2_wp3: int,
    rtm, wprtp, thlm, wpthlp, sclrm, wpsclrp, um, upwp, vm, vpwp,
    um_pert, vm_pert, upwp_pert, vpwp_pert,
    nu_vert_res_dep: NuVertResDep,
    pdf_implicit_coefs_terms: implicit_coefs_terms,
    err_info: ErrInfo,
):
    """Advance mean fields and turbulent fluxes one model timestep.

    Pushes required Python UDT mirrors into `derived_type_storage`
    before the Fortran call.
    """
    set_fortran_grid(gr)
    set_fortran_nu_vert_res_dep(nu_vert_res_dep)
    set_fortran_implicit_coefs(pdf_implicit_coefs_terms)
    set_fortran_err_info(err_info)

    result = clubb_f2py.f2py_advance_xm_wpxp(
        int(sclr_dim), f_arr(sclr_tol), float(dt),
        f_arr(sigma_sqd_w), f_arr(wm_zm), f_arr(wm_zt), f_arr(wp2), f_arr(lscale_zm),
        f_arr(wp3_on_wp2), f_arr(wp3_on_wp2_zt), f_arr(kh_zt), f_arr(kh_zm),
        f_arr(stability_correction),
        f_arr(invrs_tau_c6_zm), f_arr(tau_max_zm), f_arr(skw_zm), f_arr(wp2rtp), f_arr(rtpthvp),
        f_arr(rtm_forcing), f_arr(wprtp_forcing), f_arr(rtm_ref), f_arr(wp2thlp),
        f_arr(thlpthvp), f_arr(thlm_forcing), f_arr(wpthlp_forcing), f_arr(thlm_ref),
        f_arr(rho_ds_zm), f_arr(rho_ds_zt), f_arr(invrs_rho_ds_zm), f_arr(invrs_rho_ds_zt),
        f_arr(thv_ds_zm), f_arr(rtp2), f_arr(thlp2),
        f_arr(w_1_zm), f_arr(w_2_zm), f_arr(varnce_w_1_zm), f_arr(varnce_w_2_zm), f_arr(mixt_frac_zm),
        l_implemented, f_arr(em), f_arr(wp2sclrp), f_arr(sclrpthvp), f_arr(sclrm_forcing),
        f_arr(sclrp2), f_arr(cx_fnc_richardson),
        f_arr(um_forcing), f_arr(vm_forcing), f_arr(ug), f_arr(vg), f_arr(wpthvp),
        f_arr(fcor), f_arr(fcor_y), f_arr(um_ref), f_arr(vm_ref), f_arr(up2), f_arr(vp2),
        f_arr(uprcp), f_arr(vprcp), f_arr(rc_coef_zm),
        f_arr(clubb_params), float(ts_nudge),
        int(iipdf_type), int(penta_solve_method), int(tridiag_solve_method), int(fill_holes_type),
        l_predict_upwp_vpwp, l_ho_nontrad_coriolis, l_ho_trad_coriolis,
        l_diffuse_rtm_and_thlm, l_stability_correct_kh_n2_zm,
        l_godunov_upwind_wpxp_ta, l_upwind_xm_ma, l_uv_nudge,
        l_tke_aniso, l_diag_lscale_from_tau, l_use_c7_richardson,
        l_lmm_stepping, l_enable_relaxed_clipping, l_linearize_pbl_winds,
        l_mono_flux_lim_thlm, l_mono_flux_lim_rtm, l_mono_flux_lim_um,
        l_mono_flux_lim_vm, l_mono_flux_lim_spikefix,
        int(order_xm_wpxp), int(order_xp2_xpyp), int(order_wp2_wp3),
        f_arr(rtm), f_arr(wprtp), f_arr(thlm), f_arr(wpthlp), f_arr(sclrm), f_arr(wpsclrp),
        f_arr(um), f_arr(upwp), f_arr(vm), f_arr(vpwp),
        f_arr(um_pert), f_arr(vm_pert), f_arr(upwp_pert), f_arr(vpwp_pert),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
    return *result, get_fortran_err_info()
