"""User-facing wrappers for routines from CLUBB_core/turbulent_adv_pdf.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def xpyp_term_ta_pdf_lhs(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, coef_wpxpyp_implicit, rho_ds_zt, rho_ds_zm,
    invrs_rho_ds_zm, l_upwind_xpyp_turbulent_adv: bool, sgn_turbulent_vel, coef_wpxpyp_implicit_zm,
):
    """Compute the implicit turbulent-advection lhs coefficients for xpyp-family equations."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_xpyp_term_ta_pdf_lhs(
        f_arr(coef_wpxpyp_implicit), f_arr(rho_ds_zt), f_arr(rho_ds_zm), f_arr(invrs_rho_ds_zm),
        l_upwind_xpyp_turbulent_adv, f_arr(sgn_turbulent_vel), f_arr(coef_wpxpyp_implicit_zm),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )


def xpyp_term_ta_pdf_lhs_godunov(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, coef_wpxpyp_implicit, invrs_rho_ds_zm, rho_ds_zm
):
    """Compute the Godunov-form implicit turbulent-advection lhs coefficients."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_xpyp_term_ta_pdf_lhs_godunov(
        f_arr(coef_wpxpyp_implicit), f_arr(invrs_rho_ds_zm), f_arr(rho_ds_zm),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )


def xpyp_term_ta_pdf_rhs(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, term_wpxpyp_explicit, rho_ds_zt, rho_ds_zm,
    invrs_rho_ds_zm, l_upwind_xpyp_turbulent_adv: bool, sgn_turbulent_vel, term_wpxpyp_explicit_zm,
):
    """Compute the explicit turbulent-advection rhs term for xpyp-family equations."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_xpyp_term_ta_pdf_rhs(
        f_arr(term_wpxpyp_explicit), f_arr(rho_ds_zt), f_arr(rho_ds_zm), f_arr(invrs_rho_ds_zm),
        l_upwind_xpyp_turbulent_adv, f_arr(sgn_turbulent_vel), f_arr(term_wpxpyp_explicit_zm),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )


def xpyp_term_ta_pdf_rhs_godunov(
    gr: Grid, nzm: int, nzt: int, ngrdcol: int, term_wpxpyp_explicit_zm, invrs_rho_ds_zm,
    sgn_turbulent_vel, rho_ds_zm,
):
    """Compute the Godunov-form explicit turbulent-advection rhs term."""
    set_fortran_grid(gr)
    return clubb_f2py.f2py_xpyp_term_ta_pdf_rhs_godunov(
        f_arr(term_wpxpyp_explicit_zm), f_arr(invrs_rho_ds_zm), f_arr(sgn_turbulent_vel), f_arr(rho_ds_zm),
        nzm=int(nzm), nzt=int(nzt), ngrdcol=int(ngrdcol),
    )
