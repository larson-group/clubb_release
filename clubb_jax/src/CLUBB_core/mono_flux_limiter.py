"""JAX port of `src/CLUBB_core/mono_flux_limiter.F90`.

Description:
  Limits the value of w'x' and corrects the value of xm when the xm turbulent
  advection term is not monotonic.  A monotonic turbulent advection scheme
  will not create new extrema for variable x, based only on turbulent
  advection (not considering mean advection and xm forcings).

  A monotonic turbulent advection scheme does not allow new extrema for variable
  x to be created by turbulent advection.  When w'x' would move xm outside the
  allowable range implied by nearby previous-time values, forcings, and mean
  advection, w'x' is limited and xm is re-solved.

Porting deviations:
  `calc_turb_adv_range` keeps the Fortran routine boundary, but its ascending
  and descending grid loops are split into adjacent private helpers so the
  branch-specific `lax.while_loop` carries stay small and static.
  The Fortran constant-thickness branch is represented only by the disabled
  `l_constant_thickness = False` guard because the Fortran parameter is
  `.false.`.
"""

from __future__ import annotations

from functools import partial

import jax
import sys
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level

import jax.numpy as jnp
from jax.scipy.special import erf

from clubb_jax.src.CLUBB_core.constants_clubb import (
    eps,
    one,
    unused_var,
    zero,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.model_flags import l_force_descending_solves
from clubb_jax.src.CLUBB_core.grid_class import zm2zt, zt2zm
from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs
from clubb_jax.src.CLUBB_core.matrix_solver_wrapper import tridiag_solve


# Private named constants to avoid string comparisons.
# NOTE: These values must match the values for xm_wpxp_thlm
# and xm_wpxp_rtm given in advance_xm_wpxp_module!
mono_flux_thlm = 1
mono_flux_rtm = 2
mono_flux_um = 4
mono_flux_vm = 5

ndiags3 = 3


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "solve_type",
        "tridiag_solve_method",
        "l_implemented",
        "l_upwind_xm_ma",
        "l_mono_flux_lim_spikefix",
    ),
)
def monotonic_turbulent_flux_limit(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr,
    solve_type: int,
    dt: float,
    xm_old,
    xp2,
    wm_zt,
    xm_forcing,
    rho_ds_zm,
    rho_ds_zt,
    invrs_rho_ds_zm,
    invrs_rho_ds_zt,
    xp2_threshold: float,
    xm_tol: float,
    l_implemented: bool,
    low_lev_effect,
    high_lev_effect,
    tridiag_solve_method: int,
    l_upwind_xm_ma: bool,
    l_mono_flux_lim_spikefix: bool,
    stats,
    xm,
    wpxp,
    err_info,
):
    r"""Limits the value of w'x' and corrects the value of xm when the xm turbulent
    advection term is not monotonic.  A monotonic turbulent advection scheme
    will not create new extrema for variable x, based only on turbulent
    advection (not considering mean advection and xm forcings).

    Montonic turbulent advection
    ----------------------------

    A monotonic turbulent advection scheme does not allow new extrema for
    variable x to be created (by means of turbulent advection).  In a
    monotonic turbulent advection scheme, when only the effects of turbulent
    advection are considered (neglecting forcings and mean advection), the
    value of variable x at a given point should not increase above the
    greatest value of variable x at nearby points, nor decrease below the
    smallest value of variable x at nearby points.  Nearby points are points
    that are close enough to the given point so that the value of variable x
    at the given point is effected by the values of variable x at the nearby
    points by means of transfer by turbulent winds during a time step.  Again,
    a monotonic scheme insures that advection only transfers around values of
    variable x and does not create new extrema for variable x.  A monotonic
    turbulent advection scheme is useful because the turbulent advection term
    (w'x') may go numerically unstable, resulting in large instabilities in
    the mean field (xm).  A monotonic turbulent advection scheme will limit
    the change in xm, and also in w'x'.

    The following example illustrates the concept of monotonic turbulent
    advection.  Three successive vertical grid levels are shown (k-1, k, and
    k+1).  Three point values of theta-l are listed at every vertical grid
    level.  All three vertical levels have a mean theta-l (thlm) of 288.0 K.
    A circulation is occuring (in the direction of the arrows) in the vertical
    (w wind component) and in the horizontal (u and/or v wind components),
    such that the mean value of vertical velocity (wmm) is 0, but there is a
    turbulent component such that w'^2 > 0.

    level = k+1 || --- 287.0 K --- 288.0 K --- 289.0 K --- || thlm = 288.0
                ||      / \--------------------->|         ||
                ||       |                       |         || wmm = 0; wp2 > 0
                ||       |<---------------------\ /        ||
    level = k   || --- 288.0 K --- 288.0 K --- 288.0 K --- || thlm = 288.0
                ||       |<---------------------/ \        ||
                ||       |                       |         || wmm = 0; wp2 > 0
                ||      \ /--------------------->|         ||
    level = k-1 || --- 287.5 K --- 288.0 K --- 288.5 K --- || thlm = 288.0

    Neglecting any contributions from thlm forcings (effects of radiation,
    microphysics, large-scale horizontal advection, etc.), the values of
    theta-l as shown will be altered by only turbulent advection.  As a side
    note, the contribution of mean advection will be 0 since wmm = 0.  The
    diagram shows that the value of theta-l at the point on the right at level
    k will increase.  However, the values of theta-l at the other two points
    at level k will remain the same.  Thus, the value of thlm at level k will
    become greater than 288.0 K.  In the same manner, the values of thlm at
    the other two vertical levels (k-1 and k+1) will become smaller than
    288.0 K.  However, the monotonic turbulent advection scheme insures that
    any theta-l point value cannot become smaller than the smallest theta-l
    point value (287.0 K) or larger than the largest theta-l point value
    (289.0 K).  Since all theta-l point values must fall between 287.0 K and
    289.0 K, the level averages of theta-l (thlm) must fall between 287.0 K
    and 289.0 K.  Thus, any values of the turbulent flux, w'th_l', that would
    cause thlm to rise above 289.0 K or fall below 287.0 K, not considering
    the effect of other terms on thlm (such as forcings), are faulty and need
    to be limited appropriately.  The values of thlm also need to be corrected
    appropriately.

    Formula for the limitation of w'x' and xm
    -----------------------------------------

    The equation for change in the mean field, xm, over time is:

    d(xm)/dt = -w*d(xm)/dz - (1/rho_ds) * d( rho_ds * w'x' )/dz + xm_forcing;

    where w*d(xm)/dz is the mean advection component,
    (1/rho_ds) * d( rho_ds * w'x' )/dz is the turbulent advection component,
    and xm_forcing is the xm forcing component.  The d(xm)/dt time tendency
    component is discretized as:

    xm(k,<t+1>)/dt = xm(k,<t>)/dt - w*d(xm)/dz
                     - (1/rho_ds) * d( rho_ds * w'x' )/dz + xm_forcing.

    The value of xm after it has been advanced to timestep (t+1) must be in an
    appropriate range based on the values of xm at timestep (t), the amount of
    xm forcings applied over the ensuing time step, and the amount of mean
    advection applied over the ensuing time step.  This is exactly the same
    thing as saying that the value of xm(k,<t+1>), with the contribution of
    turbulent advection included, must fall into a certain range based on the
    value of xm(k,<t+1>) without the contribution of the turbulent advection
    component over the last time step.  The following inequality is used to
    limit the value of xm(k,<t+1>):

    MIN{ xm(k-1,<t>) + dt*xm_forcing(k-1) - dt*wm_zt(k-1)*d(xm)/dz|_(k-1)
            - x_max_dev_low(k-1,<t>),
         xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
            - x_max_dev_low(k,<t>),
         xm(k+1,<t>) + dt*xm_forcing(k+1) - dt*wm_zt(k+1)*d(xm)/dz|_(k+1)
            - x_max_dev_low(k+1,<t>) }
    <= xm(k,<t+1>) <=
    MAX{ xm(k-1,<t>) + dt*xm_forcing(k-1) - dt*wm_zt(k-1)*d(xm)/dz|_(k-1)
            + x_max_dev_high(k-1,<t>),
         xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
            + x_max_dev_high(k,<t>),
         xm(k+1,<t>) + dt*xm_forcing(k+1) - dt*wm_zt(k+1)*d(xm)/dz|_(k+1)
            + x_max_dev_high(k+1,<t>) };

    where x_max_dev_low is the absolute value of the deviation from the mean
    of the smallest point value of variable x at the given vertical level and
    timestep; and where x_max_dev_high is the deviation from the mean of the
    largest point value of variable x at the given vertical level and
    timestep.  For example, at vertical level (k+1) and timestep (t):

    x_max_dev_low(k+1,<t>)  = | MIN( x(k+1,<t>) ) - xm(k+1,<t>) |;
    x_max_dev_high(k+1,<t>) = MAX( x(k+1,<t>) ) - xm(k+1,<t>).

    The inequality shown above only takes into account values from the central
    level, one-level-below the central level, and one-level-above the central
    level.  This is the minimal amount of vertical levels that can have their
    values taken into consideration.  Any vertical level that can have it's
    properties advect to the given level during the course of a single time
    step can be taken into consideration.  However, only three levels will be
    considered in this example for the sake of simplicity.

    The inequality will be written in more simple terms:

    xm_lower_lim_allowable(k) <= xm(k,<t+1>) <= xm_upper_lim_allowable(k).

    The inequality can now be related to the turbulent flux, w'x'(k,<t+1>),
    through a substitution that is made for xm(k,<t+1>), such that:

    xm(k,<t+1>) = xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
                  - dt * (1/rho_ds) * d( rho_ds * w'x' )/dz|_(k).

    The inequality becomes:

    xm_lower_lim_allowable(k)
    <=
       xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
       - dt * (1/rho_ds) * d( rho_ds * w'x' )/dz|_(k)
    <=
    xm_upper_lim_allowable(k).

    The inequality is rearranged, and the turbulent advection term,
    d(w'x')/dz, is discretized:

    xm_lower_lim_allowable(k)
    - [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k) ]
    <=
       - dt * (1/rho_ds_zt(k))
              * invrs_dzt(k)
                * [   rho_ds_zm(k+1) * w'x'(k+1,<t+1>)
                    - rho_ds_zm(k) * w'x'(k,<t+1>) ]
    <=
    xm_upper_lim_allowable(k)
    - [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k) ];

    where invrs_dzt(k) = 1 / ( zm(k+1) - zm(k) ).

    Multiplying the inequality by -rho_ds_zt(k)/(dt*invrs_dzt(k)):

    rho_ds_zt(k)/(dt*invrs_dzt(k))
    * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
        - xm_lower_lim_allowable(k) ]
    >=
       rho_ds_zm(k+1) * w'x'(k+1,<t+1>) - rho_ds_zm(k) * w'x'(k,<t+1>)
    >=
    rho_ds_zt(k)/(dt*invrs_dzt(k))
    * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
        - xm_upper_lim_allowable(k) ].

    Note:  The inequality symbols have been flipped due to multiplication
           involving a (-) sign.

    Adding rho_ds_zm(k) * w'x'(k,<t+1>) to the inequality:

    rho_ds_zt(k)/(dt*invrs_dzt(k))
    * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
        - xm_lower_lim_allowable(k) ]
    + rho_ds_zm(k) * w'x'(k,<t+1>)
    >= rho_ds_zm(k+1) * w'x'(k+1,<t+1>) >=
    rho_ds_zt(k)/(dt*invrs_dzt(k))
    * [ xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
        - xm_upper_lim_allowable(k) ]
    + rho_ds_zm(k) * w'x'(k,<t+1>).

    The inequality is then rearranged to be based around w'x'(k+1,<t+1>):

    (1/rho_ds_zm(k+1))
    * [ rho_ds_zt(k)/(dt*invrs_dzt(k))
        * { xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
            - xm_lower_lim_allowable(k) }
        + rho_ds_zm(k) * w'x'(k,<t+1>) ]
    >=   w'x'(k+1,<t+1>)   >=
    (1/rho_ds_zm(k+1))
    * [ rho_ds_zt(k)/(dt*invrs_dzt(k))
        * { xm(k,<t>) + dt*xm_forcing(k) - dt*wm_zt(k)*d(xm)/dz|_(k)
            - xm_upper_lim_allowable(k) }
        + rho_ds_zm(k) * w'x'(k,<t+1>) ].

    The values of w'x' are found on the momentum levels, while the values of
    xm are found on the thermodynamic levels.  Additionally, the values of
    rho_ds_zm are found on the momentum levels, and the values of rho_ds_zt
    are found on the thermodynamic levels.  The inequality is applied to
    w'x'(k+1,<t+1>) from vertical levels (k+1=) 2 through the second-highest
    level, (k+1=) gr.nzm-1.  The value of w'x' at level 1 is a set surface
    (or lowest level) flux.  The value of w'x' at the highest level is also a
    set value, and therefore is not altered.

    Approximating maximum and minimum values of x at any given vertical level
    -------------------------------------------------------------------------

    The CLUBB code provides means, variances, and covariances for certain
    variables at all vertical levels.  However, there is no way to find the
    maximum or minimum point value of any variable on any vertical level.
    Without that information, x_max_dev_low and x_max_dev_high can't be found,
    and the inequality above is useless.  However, there is a way to
    approximate the maximum and minimum point values at any given vertical
    level.  The maximum and minimum point values can be approximated through
    the use of the variance, x'^2.

    Just as the mean value of x, which is xm, and the turbulent flux of x,
    which is w'x', are known, so is the variance of x, which is x'^2.  The
    standard deviation of x is the square root of the variance of x.  The
    distribution of x along the horizontal plane (at vertical level k) is
    approximated to be the sum of two normal (or Gaussian) distributions.
    Most of the values in a normal distribution are found within 2 standard
    deviations from the mean.  Thus, the maximum point value of x along the
    horizontal plance at any vertical level can be approximated as:
    xm + 2*sqrt(x'^2).  Likewise, the minimum value of x along the horizontal
    plane at any vertical level can be approximated as:  xm - 2*sqrt(x'^2).

    The values of x'^2 are found on the momentum levels.  The values of xm
    are found on the thermodynamic levels.  Thus, the values of x'^2 are
    interpolated to the thermodynamic levels in order to find the maximum
    and minimum point values of variable x.

    The one downfall of this method is that instabilities can arise in the
    model where unphysically large values of x'^2 are produced.  Thus, this
    allows for an unphysically large deviation of xm from its values at the
    previous time step due to turbulent advection.  Thus, for purposes of
    determining the maximum and minimum point values of x, a upper limit
    is placed on x'^2, in order to limit the standard deviation of x.  This
    limit is only applied in this subroutine, and is not applied to x'^2
    elsewhere in the model code.

    The vertical scan over levels uses jax.lax.scan.

    Returns:
        (xm, wpxp, err_info, stats)

    Fortran declarations (monotonic_turbulent_flux_limit in src/CLUBB_core/mono_flux_limiter.F90):

        !----------------------- Input Variables -----------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        type (grid), intent(in) :: gr

        integer, intent(in) ::  &
          solve_type  ! Variables being solved for.

        real( kind = core_rknd ), intent(in) :: &
          dt          ! Model timestep length                           [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          xm_old,          & ! xm at previous time step (thermo. levs.) [units vary]
          wm_zt,           & ! w wind component on thermodynamic levels [m/s]
          xm_forcing,      & ! xm forcings (thermodynamic levels)       [units vary]
          rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
          invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          xp2,             & ! x'^2 (momentum levels)                   [units vary]
          rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
          invrs_rho_ds_zm    ! Inv. dry, static density @ moment. levs. [m^3/kg]

        real( kind = core_rknd ), intent(in) :: &
          xp2_threshold, &   ! Lower limit of x'^2                      [units vary]
          xm_tol             ! Lower limit of maxdev                    [units vary]

        logical, intent(in) :: &
          l_implemented   ! True if CLUBB is being implemented and run in a host model

        integer, dimension(ngrdcol,nzt), intent(in) :: &
          low_lev_effect, & ! Index of lowest level that has an effect (for lev. k)
          high_lev_effect   ! Index of highest level that has an effect (for lev. k)

        integer, intent(in) :: &
          tridiag_solve_method  ! Specifier for method to solve tridiagonal systems

        logical, intent(in) :: &
          l_upwind_xm_ma, & ! This flag determines whether we want to use an upwind differencing
                            ! approximation rather than a centered differencing for turbulent or
                            ! mean advection terms. It affects rtm, thlm, sclrm, um and vm.
          l_mono_flux_lim_spikefix ! Flag to implement monotonic flux limiter code that
                                   ! eliminates spurious drying tendencies at model top

        !----------------------- Input/Output Variables -----------------------
        !------------------------- Inout Variables -------------------------
        type(stats_type), intent(inout) :: &
          stats

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
          xm          ! xm at current time step (thermodynamic levels)  [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
          wpxp        ! w'x' (momentum levels)                          [units vary]

        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
    """
    low_lev_effect = jnp.asarray(low_lev_effect, dtype=jnp.int32)
    high_lev_effect = jnp.asarray(high_lev_effect, dtype=jnp.int32)

    name_wpxp_mfl = ""
    name_xm_mfl = ""
    name_xm_enter_mfl = ""
    name_xm_old = ""
    name_wpxp_enter_mfl = ""
    name_xm_without_ta = ""
    name_xm_mfl_min = ""
    name_xm_mfl_max = ""
    name_wpxp_mfl_min = ""
    name_wpxp_mfl_max = ""
    name_xm_exit_mfl = ""
    name_wpxp_exit_mfl = ""

    if solve_type == mono_flux_rtm:
        max_xp2 = 5.0e-6
        name_wpxp_mfl = "wprtp_mfl"
        name_xm_mfl = "rtm_mfl"
        name_xm_enter_mfl = "rtm_enter_mfl"
        name_xm_old = "rtm_old"
        name_wpxp_enter_mfl = "wprtp_enter_mfl"
        name_xm_without_ta = "rtm_without_ta"
        name_xm_mfl_min = "rtm_mfl_min"
        name_xm_mfl_max = "rtm_mfl_max"
        name_wpxp_mfl_min = "wprtp_mfl_min"
        name_wpxp_mfl_max = "wprtp_mfl_max"
        name_xm_exit_mfl = "rtm_exit_mfl"
        name_wpxp_exit_mfl = "wprtp_exit_mfl"
    elif solve_type == mono_flux_thlm:
        max_xp2 = 5.0
        name_wpxp_mfl = "wpthlp_mfl"
        name_xm_mfl = "thlm_mfl"
        name_xm_enter_mfl = "thlm_enter_mfl"
        name_xm_old = "thlm_old"
        name_wpxp_enter_mfl = "wpthlp_enter_mfl"
        name_xm_without_ta = "thlm_without_ta"
        name_xm_mfl_min = "thlm_mfl_min"
        name_xm_mfl_max = "thlm_mfl_max"
        name_wpxp_mfl_min = "wpthlp_mfl_min"
        name_wpxp_mfl_max = "wpthlp_mfl_max"
        name_xm_exit_mfl = "thlm_exit_mfl"
        name_wpxp_exit_mfl = "wpthlp_exit_mfl"
    elif solve_type == mono_flux_um:
        max_xp2 = 10.0
        name_wpxp_mfl = "upwp_mfl"
        name_xm_mfl = "um_mfl"
    elif solve_type == mono_flux_vm:
        max_xp2 = 10.0
        name_wpxp_mfl = "vpwp_mfl"
        name_xm_mfl = "vm_mfl"
    else:
        max_xp2 = 5.0

    if stats.l_sample:
        stats = stats.begin_budget(name_wpxp_mfl, wpxp / dt)
        stats = stats.begin_budget(name_xm_mfl, xm / dt)
        stats = stats.update(name_xm_enter_mfl, xm)
        stats = stats.update(name_xm_old, xm_old)
        stats = stats.update(name_wpxp_enter_mfl, wpxp)

    invrs_dt = one / dt

    # Interpolate x'^2 to thermodynamic levels.
    xp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, xp2)

    # Place an upper limit on xp2_zt.
    xp2_zt = jnp.minimum(jnp.maximum(xp2_zt, xp2_threshold), max_xp2)

    # Store the value of xm as it enters the mfl.
    xm_enter_mfl = xm

    # Find the maximum and minimum usable values of variable x at each
    # vertical level.
    max_dev = jnp.maximum(2.0 * jnp.sqrt(xp2_zt), xm_tol)

    # Shut off to avoid using new, possibly corrupt mean advection term.
    xm_without_ta = xm_old + dt * xm_forcing

    # Since variable x must be one of theta_l, r_t, or a scalar, all of which
    # are positive definite quantities, the value must be >= 0.  Variable x
    # may also be u or v, in which case it is not positive definite.
    is_uv = solve_type in (mono_flux_um, mono_flux_vm)
    min_x_allowable_lev = jnp.where(
        is_uv,
        xm_without_ta - max_dev,
        jnp.maximum(xm_without_ta - max_dev, zero_threshold),
    )
    max_x_allowable_lev = xm_without_ta + max_dev

    # Find the smallest and largest value of all relevant levels for variable x.
    low_lev = jnp.minimum(low_lev_effect, high_lev_effect)
    high_lev = jnp.maximum(low_lev_effect, high_lev_effect)
    low_lev = jnp.maximum(low_lev, 0)
    high_lev = jnp.minimum(high_lev, nzt - 1)
    level_idx = jnp.arange(nzt)
    window = (
        (level_idx[None, None, :] >= low_lev[:, :, None])
        & (level_idx[None, None, :] <= high_lev[:, :, None])
    )
    min_x_allowable = jnp.min(
        jnp.where(window, min_x_allowable_lev[:, None, :], jnp.inf),
        axis=2,
    )
    max_x_allowable = jnp.max(
        jnp.where(window, max_x_allowable_lev[:, None, :], -jnp.inf),
        axis=2,
    )

    wpxp_thresh_term_zt = (
        invrs_dt * gr.grid_dir * gr.dzt * (xm_without_ta - min_x_allowable)
    )
    wpxp_mfl_max_term_zt = rho_ds_zt * wpxp_thresh_term_zt
    wpxp_mfl_min_term_zt = (
        rho_ds_zt * invrs_dt * gr.grid_dir * gr.dzt * (xm_without_ta - max_x_allowable)
    )

    # Interpolate wpxp_thresh_term_zt to momentum levels.
    wpxp_thresh_term = zt2zm(nzm, nzt, ngrdcol, gr, wpxp_thresh_term_zt)

    if gr.grid_dir_indx > 0:
        scan_levels = jnp.arange(1, nzm - 1)
        initial_neighbor = wpxp[:, gr.k_lb_zm]
    else:
        scan_levels = jnp.arange(nzm - 2, 0, -1)
        initial_neighbor = wpxp[:, gr.k_lb_zm]

    spikefix_rtm = bool(l_mono_flux_lim_spikefix and solve_type == mono_flux_rtm)

    def clip_step(carry, k):
        previous_wpxp = carry
        k_zt = jnp.where(gr.grid_dir_indx > 0, k - 1, k)
        previous_k = k - gr.grid_dir_indx

        previous_thresh = jnp.take(wpxp_thresh_term, previous_k, axis=1)
        previous_rho = jnp.take(rho_ds_zm, previous_k, axis=1)
        invrs_rho = jnp.take(invrs_rho_ds_zm, k, axis=1)
        max_term = jnp.take(wpxp_mfl_max_term_zt, k_zt, axis=1)
        min_term = jnp.take(wpxp_mfl_min_term_zt, k_zt, axis=1)
        current_wpxp = jnp.take(wpxp, k, axis=1)

        spikefix_cond = (
            (jnp.abs(previous_wpxp) > previous_thresh)
            & (previous_wpxp < zero)
        )
        wpxp_mfl_max_raw = invrs_rho * (max_term + previous_rho * previous_wpxp)
        wpxp_mfl_max_k = jnp.where(spikefix_cond & spikefix_rtm, zero, wpxp_mfl_max_raw)

        over_max = current_wpxp > wpxp_mfl_max_k
        wpxp_mfl_min_k = invrs_rho * (min_term + previous_rho * previous_wpxp)
        under_min = (~over_max) & (current_wpxp < wpxp_mfl_min_k)

        adjusted_wpxp = jnp.where(
            over_max,
            wpxp_mfl_max_k,
            jnp.where(under_min, wpxp_mfl_min_k, current_wpxp),
        )
        wpxp_net_adjust_k = jnp.where(
            over_max,
            wpxp_mfl_max_k - current_wpxp,
            jnp.where(under_min, wpxp_mfl_min_k - current_wpxp, zero),
        )
        wpxp_mfl_min_stat = jnp.where(over_max, unused_var, wpxp_mfl_min_k)
        l_xm_adjustment_needed = jnp.abs(wpxp_net_adjust_k) > eps
        return adjusted_wpxp, (
            k,
            adjusted_wpxp,
            wpxp_net_adjust_k,
            wpxp_mfl_min_stat,
            wpxp_mfl_max_k,
            l_xm_adjustment_needed,
        )

    _carry, (
        k_scan,
        wpxp_scan,
        wpxp_net_adjust_scan,
        wpxp_mfl_min_scan,
        wpxp_mfl_max_scan,
        l_xm_adjustment_scan,
    ) = jax.lax.scan(clip_step, initial_neighbor, scan_levels)

    wpxp = wpxp.at[:, k_scan].set(wpxp_scan.T)
    wpxp_net_adjust = jnp.zeros((ngrdcol, nzm))
    wpxp_mfl_min = jnp.zeros((ngrdcol, nzm))
    wpxp_mfl_max = jnp.zeros((ngrdcol, nzm))
    wpxp_net_adjust = wpxp_net_adjust.at[:, k_scan].set(wpxp_net_adjust_scan.T)
    wpxp_mfl_min = wpxp_mfl_min.at[:, k_scan].set(wpxp_mfl_min_scan.T)
    wpxp_mfl_max = wpxp_mfl_max.at[:, k_scan].set(wpxp_mfl_max_scan.T)
    l_xm_adjustment_needed = jnp.any(l_xm_adjustment_scan, axis=0)
    l_any_adjustment_needed = jnp.any(l_xm_adjustment_needed)

    def adjust_xm(args):
        xm_in, err_info_in = args
        lhs_mfl_xm = mfl_xm_lhs(
            nzm,
            nzt,
            ngrdcol,
            dt,
            gr.weights_zt2zm,
            gr.invrs_dzt,
            gr.invrs_dzm,
            wm_zt,
            l_implemented,
            l_upwind_xm_ma,
            gr.grid_dir,
        )
        rhs_mfl_xm = mfl_xm_rhs(
            nzm,
            nzt,
            ngrdcol,
            dt,
            xm_old,
            wpxp,
            xm_forcing,
            gr.invrs_dzt,
            rho_ds_zm,
            invrs_rho_ds_zt,
        )
        xm_mfl, err_info_out = mfl_xm_solve(
            nzt,
            ngrdcol,
            gr,
            solve_type,
            tridiag_solve_method,
            l_implemented,
            lhs_mfl_xm,
            rhs_mfl_xm,
            err_info_in,
        )
        if clubb_at_least_debug_level(0):
            jax.lax.cond(
                err_info_out.any_fatal(),
                lambda _: jax.debug.callback(
                    lambda info: print(
                        info.global_error_header_host(),
                        "Error in monotonic_turbulent_flux_limit after calling mfl_xm_solve",
                        sep="\n",
                        file=sys.stderr,
                    ),
                    err_info_out, ordered=True,
                ),
                lambda _: None, operand=None,
            )
        xm_out = jnp.where(l_xm_adjustment_needed[:, None], xm_mfl, xm_in)

        # Ensure there are no spikes at the top of the domain.
        dz = gr.zm[:, gr.k_ub_zm] - gr.zm[:, gr.k_ub_zm - gr.grid_dir_indx]
        xm_density_weighted = (
            rho_ds_zt[:, gr.k_ub_zt]
            * (xm_out[:, gr.k_ub_zt] - xm_enter_mfl[:, gr.k_ub_zt])
            * dz
        )
        if gr.grid_dir_indx > 0:
            integral_mask = jnp.arange(nzt) < gr.k_ub_zt
        else:
            integral_mask = jnp.arange(nzt) > gr.k_ub_zt
        xm_vert_integral = jnp.sum(
            jnp.where(
                integral_mask[None, :],
                rho_ds_zt * xm_out * gr.grid_dir * gr.dzt,
                zero,
            ),
            axis=1,
        )
        l_top_spike = (
            jnp.abs(xm_out[:, gr.k_ub_zt] - xm_enter_mfl[:, gr.k_ub_zt])
            > 10.0 * xm_tol
        )
        l_small_integral = jnp.abs(xm_vert_integral) < eps
        xm_top_removed = xm_out.at[:, gr.k_ub_zt].set(xm_enter_mfl[:, gr.k_ub_zt])
        xm_adj_coef = xm_density_weighted / jnp.where(
            jnp.abs(xm_vert_integral) > zero,
            xm_vert_integral,
            one,
        )
        def print_spike_warnings(small, low):
            import numpy as np
            for _ in np.flatnonzero(small):
                print(
                    'Vertical integral of xm is zero;mfl will remove spike at top of domain,but it will no'
                    't conserve xm.',
                    file=sys.stderr,
                )
            for _ in np.flatnonzero(low):
                print("xm_adj_coef in mfl less than -0.99, mx_adj_coef set to -0.99", file=sys.stderr)
        small_warning = l_top_spike & l_small_integral
        low_warning = l_top_spike & ~l_small_integral & (xm_adj_coef < -0.99)
        jax.lax.cond(
            jnp.any(small_warning | low_warning),
            lambda _: jax.debug.callback(print_spike_warnings, small_warning, low_warning, ordered=True),
            lambda _: None, operand=None,
        )
        xm_adj_coef = jnp.maximum(xm_adj_coef, -0.99)
        xm_scaled = (xm_out * (one + xm_adj_coef)[:, None]).at[:, gr.k_ub_zt].set(
            xm_enter_mfl[:, gr.k_ub_zt]
        )
        xm_out = jnp.where(
            l_top_spike[:, None],
            jnp.where(l_small_integral[:, None], xm_top_removed, xm_scaled),
            xm_out,
        )
        return xm_out, err_info_out

    xm, err_info = jax.lax.cond(
        l_any_adjustment_needed,
        adjust_xm,
        lambda args: args,
        (xm, err_info),
    )

    if stats.l_sample:
        stats = stats.update(name_xm_without_ta, xm_without_ta)
        stats = stats.update(name_xm_mfl_min, min_x_allowable)
        stats = stats.update(name_xm_mfl_max, max_x_allowable)
        stats = stats.update(name_wpxp_mfl_min, wpxp_mfl_min)
        stats = stats.update(name_wpxp_mfl_max, wpxp_mfl_max)
        stats = stats.finalize_budget(name_wpxp_mfl, wpxp / dt)
        stats = stats.finalize_budget(name_xm_mfl, xm / dt)
        stats = stats.update(name_xm_exit_mfl, xm)
        stats = stats.update(name_wpxp_exit_mfl, wpxp)

    return xm, wpxp, err_info, stats


@partial(
    jax.jit,
    static_argnames=("nzm", "nzt", "ngrdcol", "l_implemented", "l_upwind_xm_ma"),
)
def mfl_xm_lhs(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    dt: float,
    weights_zt2zm,
    invrs_dzt,
    invrs_dzm,
    wm_zt,
    l_implemented: bool,
    l_upwind_xm_ma: bool,
    grid_dir: float,
):
    """This subroutine is part of the process of re-solving for xm at timestep
    index (t+1).  This is done because the original solving process produced
    values outside of what is deemed acceptable by the monotonic flux limiter.
    Unlike the original formulation for advancing xm one timestep, which
    combines w'x' and xm in a band-diagonal solver, this formulation uses a
    tridiagonal solver to solve for only the value of xm(t+1), for w'x'(t+1)
    is known.

    Subroutine mfl_xm_lhs sets up the left-hand side of the matrix equation.

    Returns:
        lhs

    Fortran declarations (mfl_xm_lhs in src/CLUBB_core/mono_flux_limiter.F90):

        !---------------------------- Input Variables ----------------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        real( kind = core_rknd ), intent(in) :: &
          dt     ! Model timestep length                      [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          wm_zt         ! w wind component on thermodynamic levels   [m/s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          invrs_dzt

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          invrs_dzm

        real( kind = core_rknd ), dimension(ngrdcol,nzm,t_above:t_below), intent(in) :: &
          weights_zt2zm

        logical, intent(in) :: &
          l_implemented   ! True if CLUBB is being implemented and run in a host model

        logical, intent(in) :: &
          l_upwind_xm_ma ! This flag determines whether we want to use an upwind differencing
                         ! approximation rather than a centered differencing for turbulent or
                         ! mean advection terms. It affects rtm, thlm, sclrm, um and vm.

        real( kind = core_rknd ), intent(in) :: &
          grid_dir    ! Grid direction multiplier (+1 = ascending; -1 = descending)

        !---------------------------- Output Variables ----------------------------
        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzt), intent(out) :: &
          lhs    ! Left hand side of tridiagonal matrix
    """
    if not l_implemented:
        lhs = term_ma_zt_lhs(
            nzm,
            nzt,
            ngrdcol,
            wm_zt,
            weights_zt2zm,
            invrs_dzt,
            invrs_dzm,
            l_upwind_xm_ma,
            grid_dir,
        )
    else:
        lhs = jnp.zeros((ndiags3, ngrdcol, nzt))

    # LHS xm time tendency.
    lhs = lhs.at[1, :, :].add(one / dt)
    return lhs


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def mfl_xm_rhs(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    dt: float,
    xm_old,
    wpxp,
    xm_forcing,
    invrs_dzt,
    rho_ds_zm,
    invrs_rho_ds_zt,
):
    """This subroutine is part of the process of re-solving for xm at timestep
    index (t+1).  This is done because the original solving process produced
    values outside of what is deemed acceptable by the monotonic flux limiter.
    Unlike the original formulation for advancing xm one timestep, which
    combines w'x' and xm in a band-diagonal solver, this formulation uses a
    tridiagonal solver to solve for only the value of xm(t+1), for w'x'(t+1)
    is known.

    Subroutine mfl_xm_rhs sets up the right-hand side of the matrix equation.

    Returns:
        rhs

    Fortran declarations (mfl_xm_rhs in src/CLUBB_core/mono_flux_limiter.F90):

        !---------------------------- Input Variables ----------------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        real( kind = core_rknd ), intent(in) :: &
          dt                 ! Model timestep length                    [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          invrs_dzt,       & ! The inverse spacing between momentum grid levels;
                             ! centered over thermodynamic grid levels.
          xm_old,          & ! xm; timestep (t) (thermodynamic levels)  [units vary]
          xm_forcing,      & ! xm forcings (thermodynamic levels)       [units vary]
          invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) ::  &
          wpxp,            & ! w'x'; timestep (t+1); limited (m-levs.)  [units vary]
          rho_ds_zm          ! Dry, static density on momentum levels   [kg/m^3]

        !---------------------------- Output Variable ----------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
          rhs         ! Right hand side of tridiagonal matrix equation
    """
    del nzm, nzt, ngrdcol
    invrs_dt = one / dt
    rhs = (
        xm_old * invrs_dt
        - invrs_rho_ds_zt
        * invrs_dzt
        * (
            rho_ds_zm[:, 1:] * wpxp[:, 1:]
            - rho_ds_zm[:, :-1] * wpxp[:, :-1]
        )
        + xm_forcing
    )
    return rhs


@partial(
    jax.jit,
    static_argnames=("nzt", "ngrdcol", "solve_type", "tridiag_solve_method", "l_implemented"),
)
def mfl_xm_solve(
    nzt: int,
    ngrdcol: int,
    gr,
    solve_type: int,
    tridiag_solve_method: int,
    l_implemented: bool,
    lhs,
    rhs,
    err_info,
):
    """This subroutine is part of the process of re-solving for xm at timestep
    index (t+1).  This is done because the original solving process produced
    values outside of what is deemed acceptable by the monotonic flux limiter.
    Unlike the original formulation for advancing xm one timestep, which
    combines w'x' and xm in a band-diagonal solver, this formulation uses a
    tridiagonal solver to solve for only the value of xm(t+1), for w'x'(t+1)
    is known.

    Subroutine mfl_xm_solve solves the tridiagonal matrix equation for xm at
    timestep index (t+1).

    Returns:
        (xm, err_info)

    Fortran declarations (mfl_xm_solve in src/CLUBB_core/mono_flux_limiter.F90):

        !---------------------------- Input Variables ----------------------------
        integer, intent(in) :: &
          nzt, &
          ngrdcol

        type( grid ), intent(in) :: &
          gr

        integer, intent(in) :: &
          solve_type  ! Variables being solved for.

        integer, intent(in) :: &
          tridiag_solve_method  ! Specifier for method to solve tridiagonal systems

        logical, intent(in) :: &
          l_implemented  ! True if CLUBB is being implemented and run in a host model

        !---------------------------- InOut Variables ----------------------------
        real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzt), intent(inout) :: &
          lhs  ! Left hand side of tridiagonal matrix

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
          rhs  ! Right hand side of tridiagonal matrix equation

        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header

        !---------------------------- Output Variables ----------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
          xm   ! Value of variable being solved for at timestep (t+1)   [units vary]
    """
    if solve_type == mono_flux_rtm:
        solve_type_str = "rtm"
    elif solve_type == mono_flux_thlm:
        solve_type_str = "thlm"
    else:
        solve_type_str = "scalars"

    if l_force_descending_solves and gr.grid_dir_indx > 0:
        # Matrix solves are bit-different between ascending and descending.
        # This ensures matrices are solved in the same (descending) order,
        # which is useful for ensuring BFBness between grid modes
        # We need to flip in both the vertical dimensions and the bands for the lhs
        lhs = lhs[::-1, :, ::-1]
        rhs = rhs[:, ::-1]

    # Solve for xm at timestep index (t+1) using the tridiagonal solver.
    err_info, xm, _rcond = tridiag_solve(
        solve_type_str,
        tridiag_solve_method,
        ngrdcol,
        nzt,
        lhs,
        rhs,
        err_info,
        l_implemented=l_implemented,
    )

    if l_force_descending_solves and gr.grid_dir_indx > 0:
        # Flip the back to the ascending direction if we forced the solve
        # to be in descending mode
        xm = xm[:, ::-1]

    if clubb_at_least_debug_level(0):
        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(
                lambda info: print(
                    info.global_error_header_host(),
                    "calling tridiag_solve in mfl_xm_solve",
                    sep="\n",
                    file=sys.stderr,
                ),
                err_info, ordered=True,
            ),
            lambda _: None, operand=None,
        )
    return xm, err_info


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def calc_turb_adv_range(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr,
    dt: float,
    w_1_zm,
    w_2_zm,
    varnce_w_1_zm,
    varnce_w_2_zm,
    mixt_frac_zm,
    stats,
):
    """Calculates the lowermost and uppermost thermodynamic grid levels that can
    effect the base (or central) thermodynamic level through the effects of
    turbulent advection over the course of one time step.  This is used as
    part of the monotonic turbulent advection scheme.

    One method is to use the vertical velocity at each level to determine the
    amount of time that it takes to travel across that particular grid level.
    The method is to keep on advancing one grid level until either (a) the
    total sum of time taken reaches or exceeds the model time step length,
    (b) the top or bottom of the model is reached, or (c) a level is reached
    where the vertical velocity component (with turbulence included) is
    oriented completely opposite of the direction of travel towards the base
    (or central) thermodynamic level.  An example of situation (c) would be,
    while starting from a higher altitude and searching downward for all
    upward vertical velocity components, encountering a strong downdraft
    where the vertical velocity at every single point is oriented downward.
    Such a situation would occur when the mean vertical velocity (wm_zm)
    exceeds any turbulent component (w') that would be oriented upwards.

    Another method is to simply set the thickness (in meters) of the layer
    that turbulent advection is allowed to act over, for purposes of the
    monotonic turbulent advection scheme.  The lowermost and uppermost
    grid level that can effect the base (or central) thermodynamic level
    is computed based on the thickness and altitude of each level.

    Fortran keeps both grid-direction loops inside this routine; the JAX
    function delegates them to private helpers so that the lax.while_loop
    carries stay small and static.

    Returns:
        (low_lev_effect, high_lev_effect, stats)

    Fortran declarations (calc_turb_adv_range in src/CLUBB_core/mono_flux_limiter.F90):

        !------------------------- Input Variables -------------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        type (grid), intent(in) :: gr

        real( kind = core_rknd ), intent(in) :: &
          dt ! Model timestep length                       [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          w_1_zm,        & ! Mean w (1st PDF component)                   [m/s]
          w_2_zm,        & ! Mean w (2nd PDF component)                   [m/s]
          varnce_w_1_zm, & ! Variance of w (1st PDF component)            [m^2/s^2]
          varnce_w_2_zm, & ! Variance of w (2nd PDF component)            [m^2/s^2]
          mixt_frac_zm     ! Weight of 1st PDF component (Sk_w dependent) [-]

        type(stats_type), intent(inout) :: &
          stats

        !------------------------- Output Variables -------------------------
        integer, dimension(ngrdcol,nzt), intent(out) :: &
          low_lev_effect, & ! Index of lowest level that has an effect (for lev. k)
          high_lev_effect   ! Index of highest level that has an effect (for lev. k)
    """
    # Toggle constant or variable thickness. The Fortran parameter is .false.
    l_constant_thickness = False
    del l_constant_thickness

    invrs_dt = one / dt
    w_min = gr.grid_dir * gr.dzm * invrs_dt

    # Find the average upwards vertical velocity and the average downwards
    # vertical velocity.
    # Note:  A level that has all vertical wind moving downwards will have a
    #        vert_vel_up value that is 0, and vice versa.
    mean_w_down_zm, mean_w_up_zm, stats = mean_vert_vel_up_down(
        nzm,
        ngrdcol,
        w_1_zm,
        w_2_zm,
        varnce_w_1_zm,
        varnce_w_2_zm,
        mixt_frac_zm,
        zero,
        w_min,
        stats,
    )

    if gr.grid_dir_indx > 0:
        low_lev_effect, high_lev_effect = _calc_turb_adv_range_ascending(
            nzt, ngrdcol, gr, dt, mean_w_down_zm, mean_w_up_zm,
        )
    else:
        low_lev_effect, high_lev_effect = _calc_turb_adv_range_descending(
            nzt, ngrdcol, gr, dt, mean_w_down_zm, mean_w_up_zm,
        )

    return low_lev_effect, high_lev_effect, stats


# Target-only split of the two grid-direction branches inside Fortran
# `calc_turb_adv_range`; keep these helpers adjacent to that routine.
def _calc_turb_adv_range_ascending(nzt, ngrdcol, gr, dt, mean_w_down_zm, mean_w_up_zm):
    """Ascending-grid branch of calc_turb_adv_range.

    Porting deviation: this helper is extracted from the active non-constant
    thickness branch of Fortran `calc_turb_adv_range` solely to keep JAX loop
    carries small and static.
    """
    low_lev_effect = jnp.zeros((ngrdcol, nzt), dtype=jnp.int32)
    high_lev_effect = jnp.zeros((ngrdcol, nzt), dtype=jnp.int32)

    def loop_body(k, effects):
        low_lev_effect, high_lev_effect = effects
        # Compute the number of levels that effect the central thermodynamic
        # level through upwards motion (traveling from lower thermodynamic
        # levels to reach the central thermodynamic level).
        start_low = jnp.full((ngrdcol,), k - 1, dtype=jnp.int32)

        def low_cond(carry):
            j, _dt_all, _low, done = carry
            return jnp.any((j >= 0) & (~done))

        def low_body(carry):
            j, dt_all, low, done = carry
            active = (j >= 0) & (~done)
            j_adj = j + 1
            # Continue if there is some component of upwards vertical velocity.
            vert_vel_up = jnp.take_along_axis(mean_w_up_zm, j_adj[:, None], axis=1)[:, 0]
            dzm = jnp.take_along_axis(gr.dzm, j_adj[:, None], axis=1)[:, 0]
            has_up = vert_vel_up > zero
            # Compute the amount of time it takes to travel one grid level
            # upwards:  delta_t = delta_z / vert_vel_up_zm.
            dt_one = gr.grid_dir * dzm / jnp.where(has_up, vert_vel_up, one)
            # Total time elapsed for crossing all grid levels that have been
            # passed, thus far.
            dt_next = jnp.where(active & has_up, dt_all + dt_one, dt_all)
            # Stop if has taken more than one model time step (overall) to
            # travel the entire extent of the current vertical grid level.
            reached = active & has_up & (dt_next >= dt)
            # Stop if there isn't a component of upwards vertical velocity.
            no_up = active & (~has_up)
            low_next = jnp.where(active, j, low)
            low_next = jnp.where(no_up, j + gr.grid_dir_indx, low_next)
            done_next = done | reached | no_up
            j_next = jnp.where(active & (~done_next), j - 1, j)
            return j_next, dt_next, low_next, done_next

        _j, _dt_all, low_k, _done = jax.lax.while_loop(
            low_cond,
            low_body,
            (
                start_low,
                jnp.zeros((ngrdcol,)),
                start_low,
                jnp.zeros((ngrdcol,), dtype=bool),
            ),
        )
        low_lev_effect = low_lev_effect.at[:, k].set(low_k)

        # Compute the number of levels that effect the central thermodynamic
        # level through downwards motion (traveling from higher thermodynamic
        # levels to reach the central thermodynamic level).
        start_high = jnp.full((ngrdcol,), k + 1, dtype=jnp.int32)

        def high_cond(carry):
            j, _dt_all, _high, done = carry
            return jnp.any((j <= nzt - 1) & (~done))

        def high_body(carry):
            j, dt_all, high, done = carry
            active = (j <= nzt - 1) & (~done)
            j_adj = j
            # Continue if there is some component of downwards vertical velocity.
            vert_vel_down = jnp.take_along_axis(mean_w_down_zm, j_adj[:, None], axis=1)[:, 0]
            dzm = jnp.take_along_axis(gr.dzm, j_adj[:, None], axis=1)[:, 0]
            has_down = vert_vel_down < zero
            # Compute the amount of time it takes to travel one grid level
            # downwards:  delta_t = - delta_z / vert_vel_down_zm.
            # Note:  There is a (-) sign in front of delta_z because the
            #        distance traveled is downwards.  Since vert_vel_down
            #        has a negative value, dt_one_grid_lev will be a
            #        positive value.
            dt_one = -gr.grid_dir * dzm / jnp.where(has_down, vert_vel_down, -one)
            # Total time elapsed for crossing all grid levels that have been
            # passed, thus far.
            dt_next = jnp.where(active & has_down, dt_all + dt_one, dt_all)
            # Stop if has taken more than one model time step (overall) to
            # travel the entire extent of the current vertical grid level.
            reached = active & has_down & (dt_next >= dt)
            # Stop if there isn't a component of downwards vertical velocity.
            no_down = active & (~has_down)
            high_next = jnp.where(active, j, high)
            high_next = jnp.where(no_down, j - gr.grid_dir_indx, high_next)
            done_next = done | reached | no_down
            j_next = jnp.where(active & (~done_next), j + 1, j)
            return j_next, dt_next, high_next, done_next

        _j, _dt_all, high_k, _done = jax.lax.while_loop(
            high_cond,
            high_body,
            (
                start_high,
                jnp.zeros((ngrdcol,)),
                start_high,
                jnp.zeros((ngrdcol,), dtype=bool),
            ),
        )
        high_lev_effect = high_lev_effect.at[:, k].set(high_k)
        return low_lev_effect, high_lev_effect

    low_lev_effect, high_lev_effect = jax.lax.fori_loop(
        1, nzt - 2, loop_body, (low_lev_effect, high_lev_effect)
    )

    # Information for levels gr%k_lb_zt, gr%k_ub_zt-gr%grid_dir_indx, and
    # gr%k_ub_zt is not needed. However, set the values at these levels for
    # purposes of not having odd values in the arrays.
    low_lev_effect = low_lev_effect.at[:, gr.k_lb_zt].set(gr.k_lb_zt)
    high_lev_effect = high_lev_effect.at[:, gr.k_lb_zt].set(gr.k_lb_zt)
    low_lev_effect = low_lev_effect.at[:, gr.k_ub_zt - gr.grid_dir_indx].set(
        gr.k_ub_zt - gr.grid_dir_indx
    )
    high_lev_effect = high_lev_effect.at[:, gr.k_ub_zt - gr.grid_dir_indx].set(gr.k_ub_zt)
    low_lev_effect = low_lev_effect.at[:, gr.k_ub_zt].set(gr.k_ub_zt)
    high_lev_effect = high_lev_effect.at[:, gr.k_ub_zt].set(gr.k_ub_zt)

    return low_lev_effect, high_lev_effect


def _calc_turb_adv_range_descending(nzt, ngrdcol, gr, dt, mean_w_down_zm, mean_w_up_zm):
    """Descending-grid branch of calc_turb_adv_range.

    Porting deviation: this helper is extracted from the active non-constant
    thickness branch of Fortran `calc_turb_adv_range` solely to keep JAX loop
    carries small and static.
    """
    low_lev_effect = jnp.zeros((ngrdcol, nzt), dtype=jnp.int32)
    high_lev_effect = jnp.zeros((ngrdcol, nzt), dtype=jnp.int32)

    def loop_body(i, effects):
        low_lev_effect, high_lev_effect = effects
        # Map the forward lax loop counter to Fortran's descending level order.
        k = nzt - 2 - i
        # Compute the number of levels that effect the central thermodynamic
        # level through upwards motion (traveling from lower thermodynamic
        # levels to reach the central thermodynamic level).
        start_low = jnp.full((ngrdcol,), k - gr.grid_dir_indx, dtype=jnp.int32)

        def low_cond(carry):
            j, _dt_all, _low, done = carry
            return jnp.any((j <= gr.k_lb_zt) & (~done))

        def low_body(carry):
            j, dt_all, low, done = carry
            active = (j <= gr.k_lb_zt) & (~done)
            j_adj = j
            # Continue if there is some component of upwards vertical velocity.
            vert_vel_up = jnp.take_along_axis(mean_w_up_zm, j_adj[:, None], axis=1)[:, 0]
            dzm = jnp.take_along_axis(gr.dzm, j_adj[:, None], axis=1)[:, 0]
            has_up = vert_vel_up > zero
            # Compute the amount of time it takes to travel one grid level
            # upwards:  delta_t = delta_z / vert_vel_up_zm.
            dt_one = gr.grid_dir * dzm / jnp.where(has_up, vert_vel_up, one)
            dt_next = jnp.where(active & has_up, dt_all + dt_one, dt_all)
            reached = active & has_up & (dt_next >= dt)
            no_up = active & (~has_up)
            low_next = jnp.where(active, j, low)
            low_next = jnp.where(no_up, j + gr.grid_dir_indx, low_next)
            done_next = done | reached | no_up
            j_next = jnp.where(active & (~done_next), j - gr.grid_dir_indx, j)
            return j_next, dt_next, low_next, done_next

        _j, _dt_all, low_k, _done = jax.lax.while_loop(
            low_cond,
            low_body,
            (
                start_low,
                jnp.zeros((ngrdcol,)),
                start_low,
                jnp.zeros((ngrdcol,), dtype=bool),
            ),
        )
        low_lev_effect = low_lev_effect.at[:, k].set(low_k)

        # Compute the number of levels that effect the central thermodynamic
        # level through downwards motion (traveling from higher thermodynamic
        # levels to reach the central thermodynamic level).
        start_high = jnp.full((ngrdcol,), k + gr.grid_dir_indx, dtype=jnp.int32)

        def high_cond(carry):
            j, _dt_all, _high, done = carry
            return jnp.any((j >= gr.k_ub_zt) & (~done))

        def high_body(carry):
            j, dt_all, high, done = carry
            active = (j >= gr.k_ub_zt) & (~done)
            j_adj = j + 1
            # Continue if there is some component of downwards vertical velocity.
            vert_vel_down = jnp.take_along_axis(mean_w_down_zm, j_adj[:, None], axis=1)[:, 0]
            dzm = jnp.take_along_axis(gr.dzm, j_adj[:, None], axis=1)[:, 0]
            has_down = vert_vel_down < zero
            # Compute the amount of time it takes to travel one grid level
            # downwards:  delta_t = - delta_z / vert_vel_down_zm.
            # Note:  There is a (-) sign in front of delta_z because the
            #        distance traveled is downwards.  Since vert_vel_down
            #        has a negative value, dt_one_grid_lev will be a
            #        positive value.
            dt_one = -gr.grid_dir * dzm / jnp.where(has_down, vert_vel_down, -one)
            dt_next = jnp.where(active & has_down, dt_all + dt_one, dt_all)
            reached = active & has_down & (dt_next >= dt)
            no_down = active & (~has_down)
            high_next = jnp.where(active, j, high)
            high_next = jnp.where(no_down, j - gr.grid_dir_indx, high_next)
            done_next = done | reached | no_down
            j_next = jnp.where(active & (~done_next), j + gr.grid_dir_indx, j)
            return j_next, dt_next, high_next, done_next

        _j, _dt_all, high_k, _done = jax.lax.while_loop(
            high_cond,
            high_body,
            (
                start_high,
                jnp.zeros((ngrdcol,)),
                start_high,
                jnp.zeros((ngrdcol,), dtype=bool),
            ),
        )
        high_lev_effect = high_lev_effect.at[:, k].set(high_k)
        return low_lev_effect, high_lev_effect

    low_lev_effect, high_lev_effect = jax.lax.fori_loop(
        0, nzt - 3, loop_body, (low_lev_effect, high_lev_effect)
    )

    low_lev_effect = low_lev_effect.at[:, gr.k_lb_zt].set(gr.k_lb_zt)
    high_lev_effect = high_lev_effect.at[:, gr.k_lb_zt].set(gr.k_lb_zt)
    low_lev_effect = low_lev_effect.at[:, gr.k_ub_zt - gr.grid_dir_indx].set(
        gr.k_ub_zt - gr.grid_dir_indx
    )
    high_lev_effect = high_lev_effect.at[:, gr.k_ub_zt - gr.grid_dir_indx].set(gr.k_ub_zt)
    low_lev_effect = low_lev_effect.at[:, gr.k_ub_zt].set(gr.k_ub_zt)
    high_lev_effect = high_lev_effect.at[:, gr.k_ub_zt].set(gr.k_ub_zt)
    return low_lev_effect, high_lev_effect


@partial(jax.jit, static_argnames=("nzm", "ngrdcol"))
def mean_vert_vel_up_down(
    nzm: int,
    ngrdcol: int,
    w_1_zm,
    w_2_zm,
    varnce_w_1_zm,
    varnce_w_2_zm,
    mixt_frac_zm,
    w_ref: float,
    w_min,
    stats,
):
    """Description
    The values of vertical velocity, along a horizontal plane at any given
    vertical level, are not allowed by CLUBB to be uniform.  In other words,
    there must be some variance in vertical velocity.  This subroutine
    calculates the mean of all values of vertical velocity, at any given
    vertical level, that are greater than a certain reference velocity.  This
    subroutine also calculates the mean of all values of vertical velocity, at
    any given vertical level, that are less than a certain reference velocity.
    The reference velocity is usually 0 m/s, in which case this subroutine
    calculates the average positive (upward) velocity and the average negative
    (downward) velocity.  However, the reference velocity may be other values,
    such as wm_zm, which is the overall mean vertical velocity.  If the
    reference velocity is wm_zm, this subroutine calculates the average of all
    values of w that are on the positive ("upward") side of the mean and the
    average of all values of w that are on the negative ("downward") side of
    the mean.  These mean positive and negative vertical velocities are useful
    in determining how long, on average, it takes a parcel of air, being
    driven by subgrid updrafts or downdrafts, to traverse the length of the
    vertical grid level.

    Method
    ------

    The CLUBB model uses a joint PDF of vertical velocity, liquid water
    potential temperature, and total water mixing ratio to determine subgrid
    variability.

    The values of vertical velocity, w, along an undefined horizontal plane
    at any vertical level, are considered to approximately follow a
    distribution that is a mixture of two normal (or Gaussian) distributions.
    The values of w that are a part of the 1st normal distribution are
    referred to as w_1, and the values of w that are part of the 2nd normal
    distribution are referred to as w_2.  Note that these distributions
    overlap, and there are many values of w that are found in both w_1 and w_2.

    The probability density function (PDF) for w, P(w), is:

    P(w) = mixt_frac*P(w_1) + (1-mixt_frac)*P(w_2);

    where "mixt_frac" is the weight of the 1st normal distribution, and P(w_1) and
    P(w_2) are the equations for the 1st and 2nd normal distributions,
    respectively:

    P(w_1) = 1 / ( sigma_w_1 * sqrt(2*PI) )
            * EXP[ -(w_1-mu_w_1)^2 / (2*sigma_w_1^2) ]; and

    P(w_2) = 1 / ( sigma_w_2 * sqrt(2*PI) )
            * EXP[ -(w_2-mu_w_2)^2 / (2*sigma_w_2^2) ].

    The mean of the 1st normal distribution is mu_w_1, and the standard
    deviation of the 1st normal distribution is sigma_w_1.  The mean of the
    2nd normal distribution is mu_w_2, and the standard deviation of the 2nd
    normal distribution is sigma_w_2.

    The average value of w, distributed according to the probability
    distribution, between limits alpha and beta, is:

    <w|_(alpha:beta)> = INT(alpha:beta) w P(w) dw.

    The average value of w over a certain domain is used to determine the
    average positive and negative (as compared to the reference velocity)
    values of w at any vertical level.

    Average Negative Vertical Velocity
    ----------------------------------

    The average of all values of w in the distribution that are below the
    reference velocity, w|_ref, is the mean value of w over the domain
    -inf <= w <= w|_ref, such that:

    <w|_(-inf:w|_ref)> = INT(-inf:w|_ref) w P(w) dw.
                       = mixt_frac * INT(-inf:w|_ref) w_1 P(w_1) dw_1
                         + (1-mixt_frac) * INT(-inf:w|_ref) w_2 P(w_2) dw_2.

    For each normal distribution in the mixture of normal distribution, i
    (where "i" can be 1 or 2):

    INT(-inf:w|_ref) wi P(wi) dwi =
      - ( sigma_wi / sqrt(2*PI) ) * EXP[ -(w|_ref-mu_wi)^2 / (2*sigma_wi^2) ]
      + mu_wi * (1/2)*[ 1 + erf( (w|_ref-mu_wi) / (sqrt(2)*sigma_wi) ) ];

    where mu_wi is the mean of w for the ith normal distribution, sigma_wi is
    the standard deviations of w for the ith normal distribution, and erf( )
    is the error function.

    The mean of all values of w <= w|_ref is:

    <w|_(-inf:w|_ref)> =
    mixt_frac * { - ( sigma_w_1 / sqrt(2*PI) )
                    * EXP[ -(w|_ref-mu_w_1)^2 / (2*sigma_w_1^2) ]
                  + mu_w_1 * (1/2)
                    *[1 + erf( (w|_ref-mu_w_1) / (sqrt(2)*sigma_w_1) )] }
    + (1-mixt_frac) * { - ( sigma_w_2 / sqrt(2*PI) )
                          * EXP[ -(w|_ref-mu_w_2)^2 / (2*sigma_w_2^2) ]
                        + mu_w_2 * (1/2)
                          *[1 + erf( (w|_ref-mu_w_2) / (sqrt(2)*sigma_w_2) )] }.

    Average Positive Vertical Velocity
    ----------------------------------

    The average of all values of w in the distribution that are above the
    reference velocity, w|_ref, is the mean value of w over the domain
    w|_ref <= w <= inf, such that:

    <w|_(w|_ref:inf)> = INT(w|_ref:inf) w P(w) dw.
                      = mixt_frac * INT(w|_ref:inf) w_1 P(w_1) dw_1
                        + (1-mixt_frac) * INT(w|_ref:inf) w_2 P(w_2) dw_2.

    For each normal distribution in the mixture of normal distribution, i
    (where "i" can be 1 or 2):

    INT(w|_ref:inf) wi P(wi) dwi =
        ( sigma_wi / sqrt(2*PI) ) * EXP[ -(w|_ref-mu_wi)^2 / (2*sigma_wi^2) ]
      + mu_wi * (1/2)*[ 1 - erf( (w|_ref-mu_wi) / (sqrt(2)*sigma_wi) ) ];

    where mu_wi is the mean of w for the ith normal distribution, sigma_wi is
    the standard deviations of w for the ith normal distribution, and erf( )
    is the error function.

    The mean of all values of w >= w|_ref is:

    <w|_(w|_ref:inf)> =
    mixt_frac * {   ( sigma_w_1 / sqrt(2*PI) )
                   * EXP[ -(w|_ref-mu_w_1)^2 / (2*sigma_w_1^2) ]
                  + mu_w_1 * (1/2)
                    *[1 - erf( (w|_ref-mu_w_1) / (sqrt(2)*sigma_w_1) )] }
    + (1-mixt_frac) * {   ( sigma_w_2 / sqrt(2*PI) )
                         * EXP[ -(w|_ref-mu_w_2)^2 / (2*sigma_w_2^2) ]
                        + mu_w_2 * (1/2)
                          *[1 - erf( (w|_ref-mu_w_2) / (sqrt(2)*sigma_w_2) )] }.

    Special Limitations:
    --------------------

    A normal distribution has a domain from -inf to inf.  However, the mixture
    of normal distributions is an approximation of the distribution of values
    of w along a horizontal plane at any given vertical level.  Vertical
    velocity, w, has absolute minimum and maximum values (that cannot be
    predicted by the PDF).  The absolute maximum and minimum for each normal
    distribution is most likely found within 2 or 3 standard deviations of the
    mean for the relevant normal distribution.  In other words, for each
    normal distribution in the mixture of normal distributions, all the values
    of w are found within 2 or 3 standard deviations on both sides of the
    mean.  Therefore, if one (or both) of the normal distributions has a mean
    that is more than 3 standard deviations away from the reference velocity,
    then that entire w distribution is found on ONE side of the reference
    velocity.

    Therefore:

    a) where mu_wi + 3*sigma_wi <= w|_ref:

          The entire ith normal distribution of w is on the negative side of
          w|_ref; and

          INT(-inf:w|_ref) wi P(wi) dwi = mu_wi; and
          INT(inf:w|_ref) wi P(wi) dwi = 0.

    b) where mu_wi - 3*sigma_wi >= w|_ref:

          The entire ith normal distribution of w is on the positive side of
          w|_ref; and

          INT(-inf:w|_ref) wi P(wi) dwi = 0; and
          INT(inf:w|_ref) wi P(wi) dwi = mu_wi.

    Notes: A value of 3 standard deviations above and below the mean of the
           ith normal distribution was chosen for the approximate maximum and
           minimum values of the ith normal distribution because 99.7% of
           values in a normal distribution are found within 3 standard
           deviations from the mean (compared to 95.4% for 2 standard
           deviations).  The value of 3 standard deviations provides for a
           reasonable estimate of the absolute maximum and minimum of w, while
           covering a great majority of the normal distribution.

           In addition to approximating the up and down components of w
           by checking if the pdfs are greater than 3 standard deviations
           from the mean, there is now a case to approximate when w is
           too small in general. The input array, w_min, contains the
           minimum values of vertical velocity that would be required
           at a given grid level for that grid box to be able to affect
           the adjacent levels. If the magnitude of w at a given level
           is less than 3 standard deviations below w_min for that level,
           then there is no significant portion of the air from that grid
           box that is capable of interacting with the next level, and
           the upward and downward components for that pdf are set to 0.

    Returns:
        (mean_w_down, mean_w_up, stats)

    Fortran declarations (mean_vert_vel_up_down in src/CLUBB_core/mono_flux_limiter.F90):

        !------------------------- Input Variables -------------------------
        integer, intent(in) :: &
          nzm, &
          ngrdcol

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          w_1_zm,        & ! Mean w (1st PDF component)                   [m/s]
          w_2_zm,        & ! Mean w (2nd PDF component)                   [m/s]
          varnce_w_1_zm, & ! Variance of w (1st PDF component)            [m^2/s^2]
          varnce_w_2_zm, & ! Variance of w (2nd PDF component)            [m^2/s^2]
          mixt_frac_zm,  & ! Weight of 1st PDF component (Sk_w dependent) [-]
          w_min            ! Minimum velocity to affect adjacent level    [m/s]

        real( kind = core_rknd ), intent(in) :: &
          w_ref          ! Reference velocity, w|_ref (normally = 0)   [m/s]

        type(stats_type), intent(inout) :: &
          stats

        !------------------------- Output Variables -------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
          mean_w_down, & ! Overall mean w (<= w|_ref)                  [m/s]
          mean_w_up      ! Overall mean w (>= w|_ref)                  [m/s]
    """
    mean_w_down_1st, mean_w_up_1st = calc_mean_w_up_down_component(
        nzm, ngrdcol, w_1_zm, varnce_w_1_zm, w_ref, w_min,
    )
    mean_w_down_2nd, mean_w_up_2nd = calc_mean_w_up_down_component(
        nzm, ngrdcol, w_2_zm, varnce_w_2_zm, w_ref, w_min,
    )

    # Overall mean of downwards w.
    mean_w_down = (
        mixt_frac_zm * mean_w_down_1st
        + (one - mixt_frac_zm) * mean_w_down_2nd
    )

    # Overall mean of upwards w.
    mean_w_up = (
        mixt_frac_zm * mean_w_up_1st
        + (one - mixt_frac_zm) * mean_w_up_2nd
    )

    if stats.l_sample:
        stats = stats.update("mean_w_up", mean_w_up)
        stats = stats.update("mean_w_down", mean_w_down)

    return mean_w_down, mean_w_up, stats


@partial(jax.jit, static_argnames=("nzm", "ngrdcol"))
def calc_mean_w_up_down_component(
    nzm: int,
    ngrdcol: int,
    w_i_zm,
    varnce_w_i_zm,
    w_ref: float,
    w_min,
):
    """This procedure is used to split the PDF component of
                  vertical velocity into upward and downward components.

                  The method used is described in the description of
                  mean_vert_vel_up_down, which calls this function.

    The Fortran routine has an optional MKL vectorized path; the JAX function
    uses array expressions directly and needs no MKL packing and unpacking.

    Returns:
        (mean_w_down_i, mean_w_up_i)

    Fortran declarations (calc_mean_w_up_down_component in src/CLUBB_core/mono_flux_limiter.F90):

        integer, intent(in) :: &
          nzm, &
          ngrdcol

        !------------------------- Input Variables -------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) ::  &
          w_i_zm,        & ! Mean of w                                   [m/s]
          varnce_w_i_zm, & ! Variance of w                               [m^2/s^2]
          w_min            ! Min. velocity required to affect adj. level [m/s]

        real( kind = core_rknd ), intent(in) ::  &
          w_ref         ! Reference velocity, w|_ref (normally = 0)   [m/s]

        !------------------------- Output Variables -------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) ::  &
          mean_w_down_i, & ! Mean w (<= w|_ref) from normal distribution [m/s]
          mean_w_up_i      ! Mean w (>= w|_ref) from normal distribution [m/s]
    """
    del ngrdcol

    invrs_sqrt_2pi = one / jnp.sqrt(2.0 * jnp.pi)
    sqrt_2 = jnp.sqrt(2.0)

    sigma_w_i_zm = jnp.sqrt(jnp.maximum(varnce_w_i_zm, zero))
    sigma_safe = jnp.where(sigma_w_i_zm > zero, sigma_w_i_zm, one)

    # The entire normal is too weak to affect adjacent grid levels
    # w is considered to be 0 in both up and down directions
    too_weak = jnp.abs(w_i_zm) + 3.0 * sigma_w_i_zm <= w_min
    # The entire normal is on the negative side of w|_ref.
    all_down = (~too_weak) & (w_i_zm + 3.0 * sigma_w_i_zm <= w_ref)
    # The entire normal is on the positive side of w|_ref.
    all_up = (~too_weak) & (~all_down) & (w_i_zm - 3.0 * sigma_w_i_zm >= w_ref)

    # The normal has significant values on both sides of w_ref.
    exp_cache = jnp.exp(-((w_ref - w_i_zm) ** 2) / (2.0 * sigma_safe ** 2))
    erf_cache = erf((w_ref - w_i_zm) / (sqrt_2 * sigma_safe))

    mean_w_down_mixed = (
        -sigma_w_i_zm * invrs_sqrt_2pi * exp_cache
        + w_i_zm * 0.5 * (one + erf_cache)
    )
    mean_w_up_mixed = (
        sigma_w_i_zm * invrs_sqrt_2pi * exp_cache
        + w_i_zm * 0.5 * (one - erf_cache)
    )

    mean_w_down_i = jnp.where(
        too_weak,
        zero,
        jnp.where(all_down, w_i_zm, jnp.where(all_up, zero, mean_w_down_mixed)),
    )
    mean_w_up_i = jnp.where(
        too_weak,
        zero,
        jnp.where(all_down, zero, jnp.where(all_up, w_i_zm, mean_w_up_mixed)),
    )

    # Upper and lower levels are not used, set to 0 to be safe and avoid NaN problems.
    boundary = (jnp.arange(nzm)[None, :] == 0) | (jnp.arange(nzm)[None, :] == nzm - 1)
    mean_w_down_i = jnp.where(boundary, zero, mean_w_down_i)
    mean_w_up_i = jnp.where(boundary, zero, mean_w_up_i)
    return mean_w_down_i, mean_w_up_i
