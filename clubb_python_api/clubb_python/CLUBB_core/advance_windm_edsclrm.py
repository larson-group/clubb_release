"""User-facing wrappers for routines from CLUBB_core/advance_windm_edsclrm_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.nu_vert_res_dep import NuVertResDep
from clubb_python.derived_types.nu_vert_res_dep_converter import set_fortran_nu_vert_res_dep
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def advance_windm_edsclrm(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, edsclr_dim: int, dt: float,
    wm_zt, km_zm, kmh_zm,
    ug, vg, um_ref, vm_ref,
    wp2, up2, vp2, um_forcing, vm_forcing,
    edsclrm_forcing,
    rho_ds_zm, invrs_rho_ds_zt,
    fcor, l_implemented: bool, ts_nudge: float, tridiag_solve_method: int,
    l_predict_upwp_vpwp: bool, l_upwind_xm_ma: bool, l_uv_nudge: bool,
    l_tke_aniso: bool, l_lmm_stepping: bool, l_linearize_pbl_winds: bool,
    order_xp2_xpyp: int, order_wp2_wp3: int, order_windm: int,
    um, vm, edsclrm, upwp, vpwp, wpedsclrp,
    um_pert, vm_pert, upwp_pert, vpwp_pert,
    nu_vert_res_dep: NuVertResDep,
    err_info: ErrInfo,
):
    """Advance mean winds and eddy scalars one model timestep.

    Pushes required Python UDT mirrors into `derived_type_storage`
    before the Fortran call.
    """
    set_fortran_grid(gr)
    set_fortran_nu_vert_res_dep(nu_vert_res_dep)
    set_fortran_err_info(err_info)

    result = clubb_f2py.f2py_advance_windm_edsclrm(
        int(edsclr_dim), float(dt),
        f_arr(wm_zt), f_arr(km_zm), f_arr(kmh_zm),
        f_arr(ug), f_arr(vg), f_arr(um_ref), f_arr(vm_ref),
        f_arr(wp2), f_arr(up2), f_arr(vp2), f_arr(um_forcing), f_arr(vm_forcing),
        f_arr(edsclrm_forcing),
        f_arr(rho_ds_zm), f_arr(invrs_rho_ds_zt),
        f_arr(fcor), l_implemented, float(ts_nudge), int(tridiag_solve_method),
        l_predict_upwp_vpwp, l_upwind_xm_ma, l_uv_nudge,
        l_tke_aniso, l_lmm_stepping, l_linearize_pbl_winds,
        int(order_xp2_xpyp), int(order_wp2_wp3), int(order_windm),
        f_arr(um), f_arr(vm), f_arr(edsclrm), f_arr(upwp), f_arr(vpwp), f_arr(wpedsclrp),
        f_arr(um_pert), f_arr(vm_pert), f_arr(upwp_pert), f_arr(vpwp_pert),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
    return *result, get_fortran_err_info()
