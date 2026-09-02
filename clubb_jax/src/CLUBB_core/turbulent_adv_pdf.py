"""JAX port of `src/CLUBB_core/turbulent_adv_pdf.F90`.

Description:
  Calculates the turbulent advection term in the predictive equation for a
  variance or covariance where turbulent advection is calculated by integrating
  over the PDF. This includes <w'thl'>, <w'rt'>, <rt'^2>, <thl'^2>, and
  <rt'thl'>, as well as passive scalar fields <w'sclr'>, <sclr'^2>,
  <sclr'rt'>, and <sclr'thl'>. CLUBB does not produce a PDF for horizontal wind
  components u and v. However, the horizontal wind variances <u'^2> and <v'^2>
  still use this code, as well.

The routines are written generally in terms of <x'y'> and
coef_wpxpyp_implicit or term_wpxpyp_explicit, but also apply to <x'^2> and
<w'x'> with their corresponding PDF coefficients and explicit terms.

Porting deviations:
  Fortran intent(out) arrays (`lhs_ta` and `rhs_ta`) are returned directly.
  The JAX LHS routines return stacked diagonal arrays in the same
  [superdiagonal, main diagonal, subdiagonal] order used by the Fortran
  constants. Fortran 1-based grid loops are represented by 0-based slices.
"""

from __future__ import annotations

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp


@partial(
    jax.jit,
    static_argnames=("nzm", "nzt", "ngrdcol", "l_upwind_xpyp_turbulent_adv"),
)
def xpyp_term_ta_pdf_lhs(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr,
    coef_wpxpyp_implicit,
    rho_ds_zt,
    rho_ds_zm,
    invrs_rho_ds_zm,
    l_upwind_xpyp_turbulent_adv: bool,
    sgn_turbulent_vel,
    coef_wpxpyp_implicit_zm,
):
    """Turbulent advection of <w'x'>, <x'^2>, and <x'y'>: implicit portion.

    The <x'y'> d/dt equation contains a turbulent advection term:

      - (1/rho_ds) * d( rho_ds * <w'x'y'> )/dz.

    The value of <w'x'y'> is written in terms of PDF parameters in the form:

      <w'x'y'> = coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit.

    The variable <x'y'> is evaluated at the (t+1) timestep, which allows the
    turbulent advection term to be expressed semi-implicitly as:

      - (1/rho_ds) * d( rho_ds * coef_wpxpyp_implicit * <x'y'>(t+1) )/dz
      - (1/rho_ds) * d( rho_ds * term_wpxpyp_explicit )/dz.

    The implicit portion is:

      - (1/rho_ds) * d( rho_ds * coef_wpxpyp_implicit * <x'y'>(t+1) )/dz.

    Note: When the implicit term is brought over to the left-hand side, the
    sign is reversed and the leading "-" in front of all implicit d[ ] / dz
    terms is changed to a "+".

    Centered Discretization:

    The values of <x'y'> are found on the momentum levels, while the values of
    coef_wpxpyp_implicit are found on thermodynamic levels.  The values of
    <x'y'> are interpolated to the intermediate thermodynamic levels as
    <x'y'>|_zt.  At the thermodynamic levels, coef_wpxpyp_implicit is
    multiplied by <x'y'>|_zt and by rho_ds_zt.  Then the derivative of that
    expression is taken over the central momentum level, where it is multiplied
    by -invrs_rho_ds_zm.

    "Upwind" Discretization:

    The values of coef_wpxpyp_implicit are interpolated to the intermediate
    momentum levels as coef_wpxpyp_implicit_zm.  The sign of the turbulent
    velocity is found on the central momentum level.  When the sign of the
    turbulent velocity is positive, the derivative involves the central
    momentum level and the momentum level immediately below it.  When the sign
    is negative, the derivative involves the central momentum level and the
    momentum level immediately above it.

    JAX adaptation: the Fortran out-argument `lhs_ta` is returned as an array
    with diagonal order [superdiagonal, main diagonal, subdiagonal].  Fortran
    1-based grid indices are represented by 0-based Python slices.
    """
    # Set lower boundary array to 0.
    if not l_upwind_xpyp_turbulent_adv:
        # Centered discretization.
        fac = invrs_rho_ds_zm[:, 1:-1] * gr.invrs_dzm[:, 1:-1]
        rho_coef_k = rho_ds_zt[:, 1:] * coef_wpxpyp_implicit[:, 1:]
        rho_coef_km1 = rho_ds_zt[:, :-1] * coef_wpxpyp_implicit[:, :-1]

        # Momentum superdiagonal: [ x xpyp(k+1,<t+1>) ]
        super_int = fac * rho_coef_k * gr.weights_zm2zt[:, 1:, 0]

        # Momentum main diagonal: [ x xpyp(k,<t+1>) ]
        main_int = fac * (
            rho_coef_k * gr.weights_zm2zt[:, 1:, 1]
            - rho_coef_km1 * gr.weights_zm2zt[:, :-1, 0]
        )

        # Momentum subdiagonal: [ x xpyp(k-1,<t+1>) ]
        sub_int = -fac * rho_coef_km1 * gr.weights_zm2zt[:, :-1, 1]
    else:
        # "Upwind" discretization.
        irho = invrs_rho_ds_zm[:, 1:-1]
        grid_sgn = gr.grid_dir * sgn_turbulent_vel[:, 1:-1]

        # The "wind" is blowing upward for an ascending grid or the "wind" is
        # blowing downward for a descending grid.
        main_up = (
            irho * gr.invrs_dzt[:, :-1]
            * rho_ds_zm[:, 1:-1] * coef_wpxpyp_implicit_zm[:, 1:-1]
        )
        sub_up = (
            -irho * gr.invrs_dzt[:, :-1]
            * rho_ds_zm[:, :-2] * coef_wpxpyp_implicit_zm[:, :-2]
        )

        # The "wind" is blowing downward for an ascending grid or the "wind" is
        # blowing upward for a descending grid.
        super_down = (
            irho * gr.invrs_dzt[:, 1:]
            * rho_ds_zm[:, 2:] * coef_wpxpyp_implicit_zm[:, 2:]
        )
        main_down = (
            -irho * gr.invrs_dzt[:, 1:]
            * rho_ds_zm[:, 1:-1] * coef_wpxpyp_implicit_zm[:, 1:-1]
        )
        zeros_int = jnp.zeros((ngrdcol, nzm - 2), dtype=jnp.float64)

        super_int = jnp.where(grid_sgn > 0.0, zeros_int, super_down)
        main_int = jnp.where(grid_sgn > 0.0, main_up, main_down)
        sub_int = jnp.where(grid_sgn > 0.0, sub_up, zeros_int)

    # Set upper boundary array to 0.
    zeros_bnd = jnp.zeros((ngrdcol, 1), dtype=jnp.float64)
    superdiag = jnp.concatenate([zeros_bnd, super_int, zeros_bnd], axis=1)
    maindiag = jnp.concatenate([zeros_bnd, main_int, zeros_bnd], axis=1)
    subdiag = jnp.concatenate([zeros_bnd, sub_int, zeros_bnd], axis=1)
    return jnp.stack([superdiag, maindiag, subdiag], axis=0)


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def xpyp_term_ta_pdf_lhs_godunov(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr,
    coef_wpxpyp_implicit,
    invrs_rho_ds_zm,
    rho_ds_zm,
):
    """Godunov-form implicit turbulent-advection lhs coefficients.

    This subroutine is a revised version of xpyp_term_ta_pdf_lhs_all. The
    revisions are maded to use the  Godunov-like upwind scheme for the
    vertical discretization of the turbulent advection term. This subroutine
    returns an array of 3 dimensional arrays, one for every grid level.

    Notes:
    This subroutine exists for testing of the Godunov-like upwind scheme.  This
    subroutine does not handle boundary conditions and sets them to 0.
    """
    # Set lower boundary array to 0.
    # Godunov-like upwind discretization.
    gd = gr.grid_dir
    zero = jnp.zeros((ngrdcol, nzm - 2), dtype=jnp.float64)
    coef_k = coef_wpxpyp_implicit[:, 1:]
    coef_km1 = coef_wpxpyp_implicit[:, :-1]
    fac = invrs_rho_ds_zm[:, 1:-1] * gr.invrs_dzm[:, 1:-1] * gd

    # Momentum superdiagonal: [ x xpyp(k+1,<t+1>) ]
    super_int = fac * rho_ds_zm[:, 2:] * jnp.minimum(zero, gd * coef_k)

    # Momentum main diagonal: [ x xpyp(k,<t+1>) ]
    main_int = (
        fac * rho_ds_zm[:, 1:-1]
        * (jnp.maximum(zero, gd * coef_k) - jnp.minimum(zero, gd * coef_km1))
    )

    # Momentum subdiagonal: [ x xpyp(k-1,<t+1>) ]
    sub_int = -fac * rho_ds_zm[:, :-2] * jnp.maximum(zero, gd * coef_km1)

    # Set upper boundary array to 0.
    zeros_bnd = jnp.zeros((ngrdcol, 1), dtype=jnp.float64)
    superdiag = jnp.concatenate([zeros_bnd, super_int, zeros_bnd], axis=1)
    maindiag = jnp.concatenate([zeros_bnd, main_int, zeros_bnd], axis=1)
    subdiag = jnp.concatenate([zeros_bnd, sub_int, zeros_bnd], axis=1)
    return jnp.stack([superdiag, maindiag, subdiag], axis=0)


@partial(
    jax.jit,
    static_argnames=("nzm", "nzt", "ngrdcol", "l_upwind_xpyp_turbulent_adv"),
)
def xpyp_term_ta_pdf_rhs(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr,
    term_wpxpyp_explicit,
    rho_ds_zt,
    rho_ds_zm,
    invrs_rho_ds_zm,
    l_upwind_xpyp_turbulent_adv: bool,
    sgn_turbulent_vel,
    term_wpxpyp_explicit_zm,
):
    """Turbulent advection of <w'x'>, <x'^2>, and <x'y'>: explicit portion.

    The explicit portion of the <x'y'> turbulent advection term is:

      - (1/rho_ds) * d( rho_ds * term_wpxpyp_explicit )/dz.

    Centered Discretization:

    The values of term_wpxpyp_explicit are found on thermodynamic levels.  At
    the thermodynamic levels, term_wpxpyp_explicit is multiplied by rho_ds_zt.
    Then the derivative of that expression is taken over the central momentum
    level, where it is multiplied by -invrs_rho_ds_zm.

    "Upwind" Discretization:

    The values of term_wpxpyp_explicit are interpolated to the intermediate
    momentum levels as term_wpxpyp_explicit_zm.  When the sign of the turbulent
    velocity is positive, the derivative involves the central momentum level
    and the momentum level immediately below it.  When the sign is negative,
    the derivative involves the central momentum level and the momentum level
    immediately above it.

    JAX adaptation: the Fortran out-argument `rhs_ta` is returned.  Fortran
    1-based grid indices are represented by 0-based Python slices.
    """
    # Set lower boundary value to 0.
    if not l_upwind_xpyp_turbulent_adv:
        # Centered discretization.
        rhs_int = (
            -invrs_rho_ds_zm[:, 1:-1] * gr.invrs_dzm[:, 1:-1]
            * (
                rho_ds_zt[:, 1:] * term_wpxpyp_explicit[:, 1:]
                - rho_ds_zt[:, :-1] * term_wpxpyp_explicit[:, :-1]
            )
        )
    else:
        # "Upwind" discretization.
        grid_sgn = gr.grid_dir * sgn_turbulent_vel[:, 1:-1]

        # The "wind" is blowing upward for an ascending grid or the "wind" is
        # blowing downward for a descending grid.
        rhs_up = (
            -invrs_rho_ds_zm[:, 1:-1] * gr.invrs_dzt[:, :-1]
            * (
                rho_ds_zm[:, 1:-1] * term_wpxpyp_explicit_zm[:, 1:-1]
                - rho_ds_zm[:, :-2] * term_wpxpyp_explicit_zm[:, :-2]
            )
        )

        # The "wind" is blowing downward for an ascending grid or the "wind" is
        # blowing upward for a descending grid.
        rhs_down = (
            -invrs_rho_ds_zm[:, 1:-1] * gr.invrs_dzt[:, 1:]
            * (
                rho_ds_zm[:, 2:] * term_wpxpyp_explicit_zm[:, 2:]
                - rho_ds_zm[:, 1:-1] * term_wpxpyp_explicit_zm[:, 1:-1]
            )
        )
        rhs_int = jnp.where(grid_sgn > 0.0, rhs_up, rhs_down)

    # Set upper boundary value to 0.
    zeros_bnd = jnp.zeros((ngrdcol, 1), dtype=jnp.float64)
    return jnp.concatenate([zeros_bnd, rhs_int, zeros_bnd], axis=1)


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def xpyp_term_ta_pdf_rhs_godunov(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr,
    term_wpxpyp_explicit_zm,
    invrs_rho_ds_zm,
    sgn_turbulent_vel,
    rho_ds_zm,
):
    """Godunov-form explicit turbulent-advection rhs term.

    This subroutine intends to add godunov upwind difference scheme based
    on xpyp_term_ta_pdf_rhs.  The revisions are maded to use the Godunov-like
    upwind scheme for the vertical discretization.
    This subroutine returns an array of values
    for every grid level.

    Optional Arguements:
    The optional arguements can be used to override the default indices.
    from_level - low index, default 2
    to level - high index, default gr%nzm-1

    Notes:
    This subroutine exists for testing of Godunov-like upwind scheme.
    THIS SUBROUTINE DOES NOT HANDLE BOUNDARY CONDITIONS AND SETS THEM TO 0
    """
    # Set lower boundary value to 0.
    gd = gr.grid_dir
    zero = jnp.zeros((ngrdcol, nzm - 2), dtype=jnp.float64)
    sgn_k = sgn_turbulent_vel[:, 1:]
    sgn_km1 = sgn_turbulent_vel[:, :-1]
    term_k = term_wpxpyp_explicit_zm[:, 1:-1]

    rhs_int = (
        -invrs_rho_ds_zm[:, 1:-1] * gr.invrs_dzm[:, 1:-1] * gd
        * (
            jnp.minimum(zero, gd * sgn_k)
            * rho_ds_zm[:, 2:] * term_wpxpyp_explicit_zm[:, 2:]
            + jnp.maximum(zero, gd * sgn_k) * rho_ds_zm[:, 1:-1] * term_k
            - jnp.minimum(zero, gd * sgn_km1) * rho_ds_zm[:, 1:-1] * term_k
            - jnp.maximum(zero, gd * sgn_km1)
            * rho_ds_zm[:, :-2] * term_wpxpyp_explicit_zm[:, :-2]
        )
    )

    # Set upper boundary value to 0.
    zeros_bnd = jnp.zeros((ngrdcol, 1), dtype=jnp.float64)
    return jnp.concatenate([zeros_bnd, rhs_int, zeros_bnd], axis=1)


__all__ = [
    "xpyp_term_ta_pdf_lhs",
    "xpyp_term_ta_pdf_lhs_godunov",
    "xpyp_term_ta_pdf_rhs",
    "xpyp_term_ta_pdf_rhs_godunov",
]
