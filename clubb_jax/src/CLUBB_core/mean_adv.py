"""JAX port of ``src/CLUBB_core/mean_adv.F90``.

Description:
Module mean_adv computes the mean advection terms for all of the time-tendency
(prognostic) equations in the CLUBB parameterization.  All of the mean
advection terms are solved for completely implicitly, and therefore become part
of the left-hand side of their respective equations.

Function term_ma_zt_lhs handles the mean advection terms for the variables
located at thermodynamic grid levels.  These variables are: rtm, thlm, wp3,
all hydrometeor species, and sclrm.

Function term_ma_zm_lhs handles the mean advection terms for the variables
located at momentum grid levels.  The variables are: wprtp, wpthlp, wp2, rtp2,
thlp2, rtpthlp, up2, vp2, wpsclrp, sclrprtp, sclrpthlp, and sclrp2.

Porting deviations:
  * Fortran fills the ``lhs_ma`` out-argument; JAX returns the array.
  * Fortran uses 1-based band indices ``kp1_*diag``, ``k_*diag``, and
    ``km1_*diag``.  JAX stacks the same diagonals in Python order
    ``[superdiagonal, main diagonal, subdiagonal]``.
  * Fortran loops over columns and levels; JAX uses array slices and
    ``jnp.where`` for the upwind branch.
"""

from __future__ import annotations

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp


@partial(
    jax.jit,
    static_argnames=("nzm", "nzt", "ngrdcol", "l_upwind_xm_ma"),
)
def term_ma_zt_lhs(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    wm_zt,
    weights_zt2zm,
    invrs_dzt,
    invrs_dzm,
    l_upwind_xm_ma: bool,
    grid_dir: float,
):
    """Mean advection of var_zt: implicit portion of the code.

    The variable "var_zt" stands for a variable that is located at
    thermodynamic grid levels.

    The d(var_zt)/dt equation contains a mean advection term:

      - w * d(var_zt)/dz.

    This term is solved for completely implicitly, such that:

      - w * d( var_zt(t+1) )/dz.

    Note: When the term is brought over to the left-hand side, the sign is
    reversed and the leading "-" in front of the term is changed to a "+".

    The timestep index (t+1) means that the value of var_zt being used is from
    the next timestep, which is being advanced to in solving the d(var_zt)/dt
    equation.

    This term is discretized as follows:

    The values of var_zt are found on the thermodynamic levels, as are the
    values of wm_zt (mean vertical velocity on thermodynamic levels).  The
    variable var_zt is interpolated to the intermediate momentum levels.  The
    derivative of the interpolated values is taken over the central
    thermodynamic level.  The derivative is multiplied by wm_zt at the central
    thermodynamic level to get the desired result.

      -----var_zt---------------------------------------------- t(k+1)

      =================var_zt(interp)========================== m(k+1)

      -----var_zt---------------------d(var_zt)/dz-----wm_zt--- t(k)

      =================var_zt(interp)========================== m(k)

      -----var_zt---------------------------------------------- t(k-1)

    The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
    with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
    The letter "t" is used for thermodynamic levels and the letter "m" is used
    for momentum levels.

      invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )

    Special discretization for upper and lower boundary levels:

    Zero derivative method: the derivative d(var_zt)/dz is set to 0 over the
    model top and the model bottom.

    This method corresponds with the "zero-flux" boundary condition option for
    eddy diffusion, where d(var_zt)/dz is set to 0 across the upper boundary.

    In order to discretize the upper boundary condition, consider a new level
    outside the model (thermodynamic level gr%nzt+1) just above the upper
    boundary level (thermodynamic level gr%nzt).  The value of var_zt at the
    level just outside the model is defined to be the same as the value of
    var_zt at thermodynamic level gr%nzt.  Therefore, the value of
    d(var_zt)/dz between the level just outside the model and the uppermost
    thermodynamic level is 0, staying consistent with the zero-flux boundary
    condition option for the eddy diffusion portion of the code.  Therefore,
    the value of var_zt at momentum level gr%nzm, which is the upper boundary
    of the model, would be the same as the value of var_zt at the uppermost
    thermodynamic level.

    The values of var_zt are found on the thermodynamic levels, as are the
    values of wm_zt (mean vertical velocity on the thermodynamic levels).  The
    variable var_zt is interpolated to momentum level gr%nzm-1, based on the
    values of var_zt at thermodynamic levels gr%nzt and gr%nzt-1.  The value
    of var_zt at momentum level gr%nzm is set equal to the value of var_zt at
    thermodynamic level gr%nzt, as described above.  The derivative of the set
    and interpolated values, d(var_zt)/dz, is taken over the central
    thermodynamic level.  The derivative is multiplied by wm_zt at the central
    thermodynamic level to get the desired result.

    For the following diagram, k = gr%nzm, which is the uppermost level of the
    model:

      --[var_zt(k+1) = var_zt(k)]----(level outside model)----- t(k+1)

      ==[var_zt(top) = var_zt(k)]===[d(var_zt)/dz|_(top) = 0]== m(k+1) Boundary

      -----var_zt(k)------------------d(var_zt)/dz-----wm_zt--- t(k)

      =================var_zt(interp)========================== m(k)

      -----var_zt(k-1)----------------------------------------- t(k-1)

    where (top) stands for the grid index of momentum level k = gr%nzm, which
    is the upper boundary of the model.

    An analogous discretization occurs at the lower boundary.

    JAX adaptation: the Fortran out-argument `lhs_ma` is returned as an array
    with diagonal order [superdiagonal, main diagonal, subdiagonal].  Fortran
    1-based grid indices are represented by 0-based Python slices.
    """
    # -------------------------- Begin Code --------------------------

    if (not l_upwind_xm_ma):  # Use centered differencing
        # Most of the interior model; normal conditions.
        fac = wm_zt[:, 1:-1] * invrs_dzt[:, 1:-1]

        # Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
        super_int = fac * weights_zt2zm[:, 2:-1, 0]

        # Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
        main_int = fac * (
            weights_zt2zm[:, 2:-1, 1]
            - weights_zt2zm[:, 1:-2, 0]
        )

        # Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
        sub_int = -fac * weights_zt2zm[:, 1:-2, 1]

        # Upper Boundary for an Ascending Grid
        # or Lower Boundary for a Descending Grid.

        # Special discretization for zero derivative method, where the
        # derivative d(var_zt)/dz over the model boundary is set to 0.
        fac_top = wm_zt[:, -1:] * invrs_dzt[:, -1:]
        super_top = jnp.zeros((ngrdcol, 1), dtype=jnp.float64)
        main_top = fac_top * (1.0 - weights_zt2zm[:, nzm - 2:nzm - 1, 0])
        sub_top = -fac_top * weights_zt2zm[:, nzm - 2:nzm - 1, 1]

        # Lower Boundary for an Ascending Grid
        # or Upper Boundary for a Descending Grid.

        # Special discretization for zero derivative method, where the
        # derivative d(var_zt)/dz over the model boundary is set to 0.
        fac_bot = wm_zt[:, :1] * invrs_dzt[:, :1]
        super_bot = fac_bot * weights_zt2zm[:, 1:2, 0]
        main_bot = -fac_bot * (1.0 - weights_zt2zm[:, 1:2, 1])
        sub_bot = jnp.zeros((ngrdcol, 1), dtype=jnp.float64)

        superdiag = jnp.concatenate([super_bot, super_int, super_top], axis=1)
        maindiag = jnp.concatenate([main_bot, main_int, main_top], axis=1)
        subdiag = jnp.concatenate([sub_bot, sub_int, sub_top], axis=1)
        return jnp.stack([superdiag, maindiag, subdiag], axis=0)

    # l_upwind_xm_ma == .true.; use "upwind" differencing

    # Most of the interior model; normal conditions.
    wm_int = wm_zt[:, 1:-1]
    grid_wm_int = grid_dir * wm_int
    invrs_dzm_k = invrs_dzm[:, 1:-2]
    invrs_dzm_kp1 = invrs_dzm[:, 2:-1]

    # Mean wind is in an upward direction and an ascending grid is used
    # or mean wind is in a downward direction and a descending grid
    # is used.
    #
    # Otherwise, mean wind is in a downward direction and an ascending grid is
    # used or mean wind is in an upward direction and a descending grid
    # is used.
    super_int = jnp.where(grid_wm_int >= 0.0, 0.0, wm_int * invrs_dzm_kp1)
    main_int = jnp.where(
        grid_wm_int >= 0.0,
        wm_int * invrs_dzm_k,
        -wm_int * invrs_dzm_kp1,
    )
    sub_int = jnp.where(grid_wm_int >= 0.0, -wm_int * invrs_dzm_k, 0.0)

    # Upper Boundary for an Ascending Grid
    # or Lower Boundary for a Descending Grid.
    wm_top = wm_zt[:, -1:]
    grid_wm_top = grid_dir * wm_top
    invrs_dzm_top = invrs_dzm[:, nzm - 2:nzm - 1]
    super_top = jnp.zeros((ngrdcol, 1), dtype=jnp.float64)
    main_top = jnp.where(grid_wm_top >= 0.0, wm_top * invrs_dzm_top, 0.0)
    sub_top = jnp.where(grid_wm_top >= 0.0, -wm_top * invrs_dzm_top, 0.0)

    # Lower Boundary for an Ascending Grid
    # or Upper Boundary for a Descending Grid.
    wm_bot = wm_zt[:, :1]
    grid_wm_bot = grid_dir * wm_bot
    invrs_dzm_bot = invrs_dzm[:, 1:2]
    super_bot = jnp.where(grid_wm_bot >= 0.0, 0.0, wm_bot * invrs_dzm_bot)
    main_bot = jnp.where(grid_wm_bot >= 0.0, 0.0, -wm_bot * invrs_dzm_bot)
    sub_bot = jnp.zeros((ngrdcol, 1), dtype=jnp.float64)

    superdiag = jnp.concatenate([super_bot, super_int, super_top], axis=1)
    maindiag = jnp.concatenate([main_bot, main_int, main_top], axis=1)
    subdiag = jnp.concatenate([sub_bot, sub_int, sub_top], axis=1)
    return jnp.stack([superdiag, maindiag, subdiag], axis=0)


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def term_ma_zm_lhs(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    wm_zm,
    invrs_dzm,
    weights_zm2zt,
):
    """Mean advection of var_zm: implicit portion of the code.

    The variable "var_zm" stands for a variable that is located at momentum
    grid levels.

    The d(var_zm)/dt equation contains a mean advection term:

      - w * d(var_zm)/dz.

    This term is solved for completely implicitly, such that:

      - w * d( var_zm(t+1) )/dz.

    Note: When the term is brought over to the left-hand side, the sign is
    reversed and the leading "-" in front of the term is changed to a "+".

    The timestep index (t+1) means that the value of var_zm being used is from
    the next timestep, which is being advanced to in solving the d(var_zm)/dt
    equation.

    This term is discretized as follows:

    The values of var_zm are found on the momentum levels, as are the values of
    wm_zm (mean vertical velocity on momentum levels).  The variable var_zm is
    interpolated to the intermediate thermodynamic levels.  The derivative of
    the interpolated values is taken over the central momentum level.  The
    derivative is multiplied by wm_zm at the central momentum level to get the
    desired result.

      =====var_zm============================================== m(k+1)

      -----------------var_zm(interp)-------------------------- t(k)

      =====var_zm=====================d(var_zm)/dz=====wm_zm=== m(k)

      -----------------var_zm(interp)-------------------------- t(k-1)

      =====var_zm============================================== m(k-1)

    The vertical indices m(k+1), t(k), m(k), t(k-1), and m(k-1) correspond
    with altitudes zm(k+1), zt(k), zm(k), zt(k-1), and zm(k-1), respectively.
    The letter "t" is used for thermodynamic levels and the letter "m" is used
    for momentum levels.

      invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )

    JAX adaptation: the Fortran out-argument `lhs_ma` is returned as an array
    with diagonal order [superdiagonal, main diagonal, subdiagonal].  Fortran
    1-based grid indices are represented by 0-based Python slices.
    """
    # -------------------------- Begin Code --------------------------

    # Set boundary array to 0
    # This is the Lower Boundary for an Ascending Grid
    # or the Upper Boundary for a Descending Grid.
    # Most of the interior model; normal conditions.
    fac = wm_zm[:, 1:-1] * invrs_dzm[:, 1:-1]

    # Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
    super_int = fac * weights_zm2zt[:, 1:, 0]

    # Momentum main diagonal: [ x var_zm(k,<t+1>) ]
    main_int = fac * (
        weights_zm2zt[:, 1:, 1]
        - weights_zm2zt[:, :-1, 0]
    )

    # Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
    sub_int = -fac * weights_zm2zt[:, :-1, 1]

    # Set boundary array to 0
    # This is the Upper Boundary for an Ascending Grid
    # or the Lower Boundary for a Descending Grid.
    zeros_bnd = jnp.zeros((ngrdcol, 1), dtype=jnp.float64)
    superdiag = jnp.concatenate([zeros_bnd, super_int, zeros_bnd], axis=1)
    maindiag = jnp.concatenate([zeros_bnd, main_int, zeros_bnd], axis=1)
    subdiag = jnp.concatenate([zeros_bnd, sub_int, zeros_bnd], axis=1)
    return jnp.stack([superdiag, maindiag, subdiag], axis=0)


__all__ = [
    "term_ma_zt_lhs",
    "term_ma_zm_lhs",
]
