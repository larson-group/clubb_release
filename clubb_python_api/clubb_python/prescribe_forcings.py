"""User-facing wrappers for forcings routines."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid
from clubb_python.derived_types.sclr_idx import SclrIdx
from clubb_python.derived_types.sclr_idx_converter import set_fortran_sclr_idx
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def prescribe_forcings(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, sclr_dim: int, edsclr_dim: int,
    sclr_idx: SclrIdx, runtype: str, sfctype: int,
    time_current: float, time_initial: float, dt: float,
    um, vm, thlm, p_in_Pa, exner, rho, rho_zm, thvm, veg_t_in_k=None,
    l_modify_bc_for_cnvg_test: bool = False,
    saturation_formula: int = 0, stats=None, l_add_dycore_grid: bool = False, grid_remap_method: int = 0,
    total_idx_rho_lin_spline: int = 0, rho_lin_spline_vals=None, rho_lin_spline_levels=None, gr_dycore: Grid | None = None,
    rtm=None, wm_zm=None, wm_zt=None, ug=None, vg=None, um_ref=None, vm_ref=None,
    thlm_forcing=None, rtm_forcing=None, um_forcing=None, vm_forcing=None,
    wprtp_forcing=None, wpthlp_forcing=None, rtp2_forcing=None, thlp2_forcing=None, rtpthlp_forcing=None,
    wpsclrp=None, sclrm_forcing=None, edsclrm_forcing=None,
    wpthlp_sfc=None, wprtp_sfc=None, upwp_sfc=None, vpwp_sfc=None,
    T_sfc=None, p_sfc=None, sens_ht: float = 0.0, latent_ht: float = 0.0,
    wpsclrp_sfc=None, wpedsclrp_sfc=None,
    err_info: ErrInfo | None = None,
):
    """Call Fortran prescribe_forcings through the F2PY wrapper."""
    if err_info is None:
        raise ValueError("prescribe_forcings requires err_info.")
    zt_in = gr.zt
    l_t_dependent = False
    l_ignore_forcings = False
    l_input_xpwp_sfc = False
    grid_adapt_in_time_method = 0
    set_fortran_grid(gr)
    set_fortran_sclr_idx(sclr_idx)
    set_fortran_err_info(err_info)

    result = clubb_f2py.f2py_prescribe_forcings(
        int(sclr_dim), int(edsclr_dim), str(runtype), int(sfctype),
        float(time_current), float(time_initial), float(dt),
        f_arr(um), f_arr(vm), f_arr(thlm), f_arr(p_in_Pa), f_arr(exner), f_arr(rho), f_arr(rho_zm), f_arr(thvm), f_arr(zt_in),
        l_t_dependent, l_ignore_forcings, l_input_xpwp_sfc,
        l_modify_bc_for_cnvg_test,
        int(saturation_formula), l_add_dycore_grid, int(grid_remap_method),
        int(grid_adapt_in_time_method),
        f_arr(rtm), f_arr(wm_zm), f_arr(wm_zt), f_arr(ug), f_arr(vg), f_arr(um_ref), f_arr(vm_ref),
        f_arr(thlm_forcing), f_arr(rtm_forcing), f_arr(um_forcing), f_arr(vm_forcing),
        f_arr(wprtp_forcing), f_arr(wpthlp_forcing), f_arr(rtp2_forcing), f_arr(thlp2_forcing),
        f_arr(rtpthlp_forcing),
        f_arr(wpsclrp), f_arr(sclrm_forcing), f_arr(edsclrm_forcing),
        f_arr(wpthlp_sfc), f_arr(wprtp_sfc), f_arr(upwp_sfc), f_arr(vpwp_sfc),
        f_arr(T_sfc), f_arr(p_sfc), float(sens_ht), float(latent_ht),
        f_arr(wpsclrp_sfc), f_arr(wpedsclrp_sfc),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
    return *result, get_fortran_err_info()
