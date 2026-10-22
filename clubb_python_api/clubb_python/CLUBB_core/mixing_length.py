"""Call the Fortran mixing-length and turbulence-time routines from Python.

Usage:
  - ``calc_lscale`` returns the complete length/time diagnosis for one timestep.
  - ``set_lscale_max`` returns a standalone or host-spacing length cap.
  - ``diagnose_lscale_from_tau`` computes method 2 and its inverse-time terms.
  - ``calc_lscale_directly`` selects parcel method 1 or coherent method 3.

Understanding:
  These wrappers require the compiled Python API. They copy the supplied grid,
  PDF parameters when needed, and error state into shared Fortran storage,
  then return a new ErrInfo followed by numeric outputs in Fortran order.
  The caller must retain that returned error state. Fortran updates the stored
  statistics context; statistics are not an extra Python return value.
  Fields use (ngrdcol,nzt) thermodynamic or (ngrdcol,nzm) momentum axes.
  Lengths are in m, times s, inverse times 1/s, pressure Pa, temperature K,
  water kg/kg, and TKE m^2/s^2. lscale_method is an integer 1/2/3, not a boolean.

References:
  ``src/CLUBB_core/mixing_length.F90`` is the implementation authority;
  ``doc/coherence_mixing_length.md`` describes coherent method 3.
"""

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
    """Return per-column length caps [m] from host_dx/host_dy [m], or 100000 m in
    standalone mode.
    """
    return clubb_f2py.f2py_set_lscale_max(
        l_implemented, f_arr(host_dx), f_arr(host_dy), ngrdcol=int(ngrdcol)
    )


def calc_lscale(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr: Grid,
    l_implemented: bool,
    host_dx, host_dy,
    p_in_pa, exner, rtm, thlm, thvm,
    thlp2, rtp2, rtpthlp,
    pdf_params: pdf_parameter,
    em, thv_ds_zt, lmin: float,
    upwp_sfc, vpwp_sfc, ddzt_umvm_sqd, ice_supersat_frac,
    ufmin: float, tau_const: float, sfc_elevation, clubb_params,
    saturation_formula: int,
    l_lscale_plume_centered: bool,
    lscale_method: int,
    l_e3sm_config: bool,
    l_smooth_heaviside_tau_wpxp: bool,
    l_modify_limiters_for_cnvg_test: bool,
    l_use_invrs_tau_n2_iso: bool,
    brunt_vaisala_freq_sqd_smth,
    err_info: ErrInfo,
):
    """Return mixing lengths and damping times for lscale_method 1, 2, or 3.

    Read initialized grid/PDF state, environmental profiles, moments, and
    per-column parameters. Return ErrInfo followed by inverse times, bounded
    times, and combined/directional lengths in the Fortran output order.
    Update the shared Fortran grid, PDF, error, and statistics contexts.
    """
    set_fortran_grid(gr)
    set_fortran_pdf_params(pdf_params)
    set_fortran_err_info(err_info)
    result = clubb_f2py.f2py_calc_lscale(
        bool(l_implemented), f_arr(host_dx), f_arr(host_dy),
        f_arr(p_in_pa), f_arr(exner), f_arr(rtm), f_arr(thlm), f_arr(thvm),
        f_arr(thlp2), f_arr(rtp2), f_arr(rtpthlp), f_arr(em),
        f_arr(thv_ds_zt), float(lmin),
        f_arr(upwp_sfc), f_arr(vpwp_sfc), f_arr(ddzt_umvm_sqd), f_arr(ice_supersat_frac),
        float(ufmin), float(tau_const), f_arr(sfc_elevation), f_arr(clubb_params),
        int(saturation_formula), bool(l_lscale_plume_centered), int(lscale_method),
        bool(l_e3sm_config), bool(l_smooth_heaviside_tau_wpxp),
        bool(l_modify_limiters_for_cnvg_test), bool(l_use_invrs_tau_n2_iso),
        f_arr(brunt_vaisala_freq_sqd_smth),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
    return get_fortran_err_info(), *result


def diagnose_lscale_from_tau(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr: Grid,
    upwp_sfc, vpwp_sfc, ddzt_umvm_sqd,
    ice_supersat_frac, em,
    ufmin: float, tau_const: float,
    sfc_elevation, lscale_max, clubb_params,
    l_e3sm_config: bool,
    l_smooth_heaviside_tau_wpxp: bool,
    brunt_vaisala_freq_sqd_smth, ri_zm,
    err_info: ErrInfo,
):
    """Return method-2 inverse-time contributions [1/s], times [s], and lengths [m].

    Use surface stress, TKE, shear, and stability on the supplied grid.
    Return ErrInfo first, followed by Fortran numeric outputs; directional
    lengths are zero placeholders. Update shared grid/error/statistics state.
    """
    set_fortran_grid(gr)
    set_fortran_err_info(err_info)
    em_zt = clubb_f2py.f2py_zm2zt_2d(
        int(nzt), f_arr(em), nzm=int(nzm), ngrdcol=int(ngrdcol))
    sqrt_em_zt = np.sqrt(np.maximum(em_zt, 0.0))
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
    ngrdcol: int,
    nzm: int,
    nzt: int,
    gr: Grid,
    l_implemented: bool,
    p_in_pa, exner, rtm, thlm, thvm,
    rtp2_zt, thlp2_zt, rtpthlp_zt, pdf_params: pdf_parameter,
    em, thv_ds_zt, lscale_max, lmin: float,
    clubb_params, saturation_formula: int,
    l_lscale_plume_centered: bool,
    lscale_method: int,
    err_info: ErrInfo,
):
    """Return ErrInfo and combined/upward/downward parcel lengths [m].

    lscale_method selects legacy method 1 or coherent method 3 from
    thermodynamic profiles, momentum-level TKE, PDF properties, and parameters.
    All lengths have (ngrdcol,nzt) shape. Update shared Fortran grid, PDF,
    error, and statistics contexts.
    """
    set_fortran_grid(gr)
    set_fortran_pdf_params(pdf_params)
    set_fortran_err_info(err_info)
    result = clubb_f2py.f2py_calc_lscale_directly(
        l_implemented,
        f_arr(p_in_pa), f_arr(exner), f_arr(rtm), f_arr(thlm), f_arr(thvm),
        f_arr(rtp2_zt), f_arr(thlp2_zt), f_arr(rtpthlp_zt),
        f_arr(em), f_arr(thv_ds_zt),
        f_arr(lscale_max), float(lmin), f_arr(clubb_params), int(saturation_formula),
        l_lscale_plume_centered, int(lscale_method),
        ngrdcol=int(ngrdcol), nzm=int(nzm), nzt=int(nzt),
    )
    return get_fortran_err_info(), *result
