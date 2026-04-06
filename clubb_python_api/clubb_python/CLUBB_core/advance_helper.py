"""User-facing wrappers for routines from CLUBB_core/advance_helper_module.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def calculate_thlp2_rad(
    gr: Grid, ngrdcol: int, nzm: int, nzt: int, rcm, thlprcp, radht, clubb_params, thlp2_forcing
):
    """Calculate radiation contribution to thlp2 forcing."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_calculate_thlp2_rad(
        f_arr(rcm), f_arr(thlprcp), f_arr(radht), f_arr(clubb_params), f_arr(thlp2_forcing),
        ngrdcol=int(ngrdcol), nzm=int(nzm), nzt=int(nzt))


def calc_ri_zm(nzm: int, ngrdcol: int, bv_freq_sqd, shear, lim_bv: float, lim_shear: float):
    """Compute Richardson number on momentum levels."""
    return clubb_f2py.f2py_calc_ri_zm(
        f_arr(bv_freq_sqd), f_arr(shear), lim_bv, lim_shear,
        nzm=int(nzm), ngrdcol=int(ngrdcol))


def vertical_avg(total_idx: int, rho_ds, field, dz) -> float:
    """Compute density-weighted vertical average of a 1D profile."""
    return float(clubb_f2py.f2py_vertical_avg(
        f_arr(rho_ds), f_arr(field), f_arr(dz), total_idx=int(total_idx)))


def vertical_integral(total_idx: int, rho_ds, field, dz) -> float:
    """Compute density-weighted vertical integral of a 1D profile."""
    return float(clubb_f2py.f2py_vertical_integral(
        f_arr(rho_ds), f_arr(field), f_arr(dz), total_idx=int(total_idx)))


def pvertinterp(gr: Grid, nzt: int, ngrdcol: int, p_mid, p_out: float, input_var):
    """Interpolate a profile from model pressure levels to target pressure."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_pvertinterp(
        f_arr(p_mid), p_out, f_arr(input_var), nzt=int(nzt), ngrdcol=int(ngrdcol))


def smooth_heaviside_peskin(nz: int, ngrdcol: int, input, smth_range: float):
    """Compute Peskin-smoothed Heaviside transform elementwise."""
    return clubb_f2py.f2py_smooth_heaviside_peskin(
        f_arr(input), float(smth_range), nz=int(nz), ngrdcol=int(ngrdcol))


def smooth_max(nz: int, ngrdcol: int, input_var1, input_var2, smth_coef: float):
    """Compute the smoothed max for supported 2D scalar/array combinations."""
    if np.isscalar(input_var1):
        return clubb_f2py.f2py_smooth_max_scalar_array(
            float(input_var1), f_arr(input_var2), float(smth_coef), nz=int(nz), ngrdcol=int(ngrdcol)
        )
    if np.isscalar(input_var2):
        return clubb_f2py.f2py_smooth_max_array_scalar(
            f_arr(input_var1), float(input_var2), float(smth_coef), nz=int(nz), ngrdcol=int(ngrdcol)
        )
    raise TypeError("smooth_max currently supports scalar/2D-array and 2D-array/scalar inputs only.")


def smooth_min(nz: int, ngrdcol: int, input_var1, input_var2, smth_coef: float):
    """Compute the smoothed min for supported 2D scalar/array combinations."""
    if np.isscalar(input_var1):
        return clubb_f2py.f2py_smooth_min_scalar_array(
            float(input_var1), f_arr(input_var2), float(smth_coef), nz=int(nz), ngrdcol=int(ngrdcol)
        )
    if np.isscalar(input_var2):
        return clubb_f2py.f2py_smooth_min_array_scalar(
            f_arr(input_var1), float(input_var2), float(smth_coef), nz=int(nz), ngrdcol=int(ngrdcol)
        )
    raise TypeError("smooth_min currently supports scalar/2D-array and 2D-array/scalar inputs only.")


def calc_xpwp(gr: Grid, nzm: int, nzt: int, ngrdcol: int, km_zm, xm):
    """Compute x'w' from grid spacing, diffusivity, and a thermo-level field."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_calc_xpwp_2d(
        f_arr(km_zm), f_arr(xm), nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol))


def lscale_width_vert_avg(gr: Grid, nzm: int, ngrdcol: int, smth_type: int, var_profile, lscale_zm,
                          rho_ds_zm, var_below_ground_value: float):
    """Compute running vertical average over half-width Lscale_zm."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_lscale_width_vert_avg(
        int(smth_type), f_arr(var_profile), f_arr(lscale_zm), f_arr(rho_ds_zm),
        float(var_below_ground_value), nzm=int(nzm), ngrdcol=int(ngrdcol))


def wp2_term_splat_lhs(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, c_wp2_splat, brunt_vaisala_freq_sqd_splat
):
    """Compute splatting lhs contribution for wp2 equation."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_wp2_term_splat_lhs(
        int(nzt), f_arr(c_wp2_splat), f_arr(brunt_vaisala_freq_sqd_splat),
        nzm=int(nzm), ngrdcol=int(ngrdcol))


def wp3_term_splat_lhs(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, c_wp2_splat, brunt_vaisala_freq_sqd_splat
):
    """Compute splatting lhs contribution for wp3 equation."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_wp3_term_splat_lhs(
        int(nzt), f_arr(c_wp2_splat), f_arr(brunt_vaisala_freq_sqd_splat),
        nzm=int(nzm), ngrdcol=int(ngrdcol))


def calc_brunt_vaisala_freq_sqd(
    gr: Grid,
    nzm: int,
    nzt: int,
    ngrdcol: int,
    thlm, exner, rtm, rcm, p_in_Pa, thvm, ice_supersat_frac,
    saturation_formula: int,
    l_brunt_vaisala_freq_moist: bool,
    l_use_thvm_in_bv_freq: bool,
    l_modify_limiters_for_cnvg_test: bool,
    bv_efold, T0: float,
):
    """Compute Brunt-Vaisala frequency squared on momentum levels."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_calc_brunt_vaisala_freq_sqd(
        int(nzm),
        f_arr(thlm), f_arr(exner), f_arr(rtm), f_arr(rcm),
        f_arr(p_in_Pa), f_arr(thvm), f_arr(ice_supersat_frac),
        saturation_formula,
        l_brunt_vaisala_freq_moist,
        l_use_thvm_in_bv_freq,
        l_modify_limiters_for_cnvg_test,
        f_arr(bv_efold), T0, nzt=int(nzt), ngrdcol=int(ngrdcol))


def compute_cx_fnc_richardson(
    gr: Grid,
    nzm: int,
    nzt: int,
    ngrdcol: int,
    Lscale_zm, ddzt_umvm_sqd, rho_ds_zm,
    brunt_vaisala_freq_sqd, brunt_vaisala_freq_sqd_mixed,
    clubb_params,
    l_use_shear_Richardson: bool,
    l_modify_limiters_for_cnvg_test: bool,
):
    """Compute Cx as function of Richardson number."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_compute_cx_fnc_richardson(
        int(nzt),
        f_arr(Lscale_zm), f_arr(ddzt_umvm_sqd), f_arr(rho_ds_zm),
        f_arr(brunt_vaisala_freq_sqd), f_arr(brunt_vaisala_freq_sqd_mixed),
        f_arr(clubb_params),
        l_use_shear_Richardson,
        l_modify_limiters_for_cnvg_test, nzm=int(nzm), ngrdcol=int(ngrdcol))


def calc_stability_correction(
    gr: Grid,
    nzm: int,
    nzt: int,
    ngrdcol: int,
    thlm, Lscale_zm, em, exner, rtm, rcm,
    p_in_Pa, thvm, ice_supersat_frac,
    lambda0_stability_coef, bv_efold, T0: float,
    saturation_formula: int,
    l_brunt_vaisala_freq_moist: bool,
    l_use_thvm_in_bv_freq: bool,
    l_modify_limiters_for_cnvg_test: bool,
):
    """Compute stability correction factor for turbulence."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_calc_stability_correction(
        f_arr(thlm), f_arr(Lscale_zm), f_arr(em), f_arr(exner),
        f_arr(rtm), f_arr(rcm), f_arr(p_in_Pa), f_arr(thvm),
        f_arr(ice_supersat_frac),
        f_arr(lambda0_stability_coef), f_arr(bv_efold), T0,
        saturation_formula,
        l_brunt_vaisala_freq_moist,
        l_use_thvm_in_bv_freq,
        l_modify_limiters_for_cnvg_test,
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol))
