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
    runtype: str, sfctype: int,
    time_current: float, time_initial: float, dt: float,
    um, vm, thlm, p_in_Pa, exner, rho, rho_zm, thvm, zt_in,
    l_t_dependent: bool, l_ignore_forcings: bool, l_input_xpwp_sfc: bool,
    l_modify_bc_for_cnvg_test: bool,
    saturation_formula: int, l_add_dycore_grid: bool, grid_remap_method: int,
    grid_adapt_in_time_method: int,
    rtm, wm_zm, wm_zt, ug, vg, um_ref, vm_ref,
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing,
    wprtp_forcing, wpthlp_forcing, rtp2_forcing, thlp2_forcing, rtpthlp_forcing,
    wpsclrp, sclrm_forcing, edsclrm_forcing,
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
    T_sfc, p_sfc, sens_ht: float, latent_ht: float,
    wpsclrp_sfc, wpedsclrp_sfc,
    sclr_idx: SclrIdx,
    err_info: ErrInfo,
):
    """Call Fortran prescribe_forcings through the F2PY wrapper."""
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
