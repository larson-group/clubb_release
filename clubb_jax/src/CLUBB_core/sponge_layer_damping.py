"""Batched sponge-damping kernels from sponge_layer_damping.F90.

Time scales and layer depths follow each column's grid heights.
Each damping call applies its own column's top-down layer exit. Standalone
configuration guards continue to reject sponge damping until driver integration.
"""
from __future__ import annotations

from typing import NamedTuple

import jax.numpy as jnp


class sponge_damp_settings(NamedTuple):
    """Time-scale bounds in seconds and sponge depth as a domain-height fraction."""
    tau_sponge_damp_min: float
    tau_sponge_damp_max: float
    sponge_damp_depth: float
    l_sponge_damping: bool = False


class sponge_damp_profile(NamedTuple):
    """Column/level time scales and column layer depths; None means uninitialized."""
    tau_sponge_damp: object
    sponge_layer_depth: object


def sponge_damp_xm(ngrdcol, nzm, nzt, gr, dt, zt, zm, xm_ref, xm, damping_profile):
    """Damps each column of a mean field toward its reference profile.
    The damping profile must be initialized before calling this function;
    otherwise a stop is issued.

    The layer boundary follows each column's current heights.

    Returns:
        xm_p

    Fortran declarations (sponge_damp_xm in src/CLUBB_core/sponge_layer_damping.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          ngrdcol, & ! Number of columns [-]
          nzm,     & ! Number of momentum levels [-]
          nzt        ! Number of thermodynamic levels [-]

        type( grid ), intent(in) :: &
          gr    ! Grid bounds and vertical direction

        real( kind = core_rknd ), intent(in) :: &
          dt    ! Model Timestep  [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          zt    ! Thermodynamic grid heights [m]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          zm    ! Momentum grid heights [m]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          xm_ref, & ! Reference mean field [units vary]
          xm        ! Mean field before damping [units vary]

        type(sponge_damp_profile), intent(in) :: &
          damping_profile    ! Column time scales and layer depths

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
          xm_p   ! Damped value of xm  [units_vary]
    """
    if damping_profile.tau_sponge_damp is None:
        raise ValueError("tau_sponge_damp in sponge_damp_xm used before initialization")
    levels = gr.k_ub_zt - gr.grid_dir_indx * jnp.arange(nzt)
    in_sponge = (zm[:, gr.k_ub_zm, None] - zt[:, levels]) < damping_profile.sponge_layer_depth[:, None]
    l_damp = jnp.cumsum(~in_sponge, axis=1) == 0
    dt_on_tau = dt / damping_profile.tau_sponge_damp[:, levels]
    damped = (xm[:, levels] + dt_on_tau * xm_ref[:, levels]) / (1.0 + dt_on_tau)
    xm_p = xm.at[:, levels].set(jnp.where(l_damp, damped, xm[:, levels]))
    return xm_p


def sponge_damp_xp2(ngrdcol, nzm, gr, dt, zm, xp2, x_tol_sqd, damping_profile):
    """Calculates "sponge"-layer damping of the variance xp2 in each column.

    Sponge damping on a local value of x is given by the equation:

    x_d = x - ( delta_t / tau ) * ( x - <x> ),

    where x is the local value prior to damping, x_d is the damped local value
    of x, <x> is the grid-level mean value of x, delta_t is the model time
    step duration, and tau is the damping time scale.  Since delta_t / tau has
    the same value everywhere at a grid level, the grid-level mean of x does
    not change as a result of damping.

    Subtracting <x> from both sides:

    x_d - <x> = ( x - <x> ) - ( delta_t / tau ) * ( x - <x> ),

    which results in:

    x_d - <x> = ( 1 - delta_t / tau ) * ( x - <x> ).

    Squaring both sides:

    ( x_d - <x> )^2 = ( 1 - delta_t / tau )^2 * ( x - <x> )^2.

    After averaging both sides, the damped value of xp2 is:

    < x_d'^2 > = ( 1 - delta_t / tau )^2 * < x'^2 >.

    Any sponge damping is applied to (predictive) xp2 after the field has been
    advanced in time.  This allows sponge damping to be applied in an implicit
    manner.  The damped value of xp2 must also be limited at a minimum value
    of x_tol^2.

    Returns:
        xp2_damped

    Fortran declarations (sponge_damp_xp2 in src/CLUBB_core/sponge_layer_damping.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          ngrdcol, & ! Number of columns [-]
          nzm        ! Number of momentum levels [-]

        type( grid ), intent(in) :: &
          gr    ! Grid bounds and vertical direction

        real( kind = core_rknd ), intent(in) :: &
          dt    ! Model Timestep  [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          zm,  & ! Momentum grid heights [m]
          xp2    ! Variance before damping [units vary]

        real( kind = core_rknd ), intent(in) :: &
          x_tol_sqd    ! Minimum variance [units vary]

        type(sponge_damp_profile), intent(in) :: &
          damping_profile    ! Column time scales and layer depths

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
          xp2_damped    ! Variance of x, <x'^2>, after damping   [units vary]
    """
    if damping_profile.tau_sponge_damp is None:
        raise ValueError("tau_sponge_damp in sponge_damp_xp2 used before initialization")
    levels = gr.k_ub_zm - gr.grid_dir_indx * jnp.arange(nzm)
    in_sponge = (zm[:, gr.k_ub_zm, None] - zm[:, levels]) < damping_profile.sponge_layer_depth[:, None]
    l_damp = jnp.cumsum(~in_sponge, axis=1) == 0
    dt_on_tau = dt / damping_profile.tau_sponge_damp[:, levels]
    damped = jnp.maximum((1.0 - dt_on_tau)**2 * xp2[:, levels], x_tol_sqd)
    xp2_damped = xp2.at[:, levels].set(jnp.where(l_damp, damped, xp2[:, levels]))
    return xp2_damped


def sponge_damp_xp3(ngrdcol, nzm, nzt, gr, dt, z, zm, xp3, damping_profile):
    """Calculates "sponge"-layer damping of the third moment xp3 in each column.

    Sponge damping on a local value of x is given by the equation:

    x_d = x - ( delta_t / tau ) * ( x - <x> ),

    where x is the local value prior to damping, x_d is the damped local value
    of x, <x> is the grid-level mean value of x, delta_t is the model time
    step duration, and tau is the damping time scale.  Since delta_t / tau has
    the same value everywhere at a grid level, the grid-level mean of x does
    not change as a result of damping.

    Subtracting <x> from both sides:

    x_d - <x> = ( x - <x> ) - ( delta_t / tau ) * ( x - <x> ),

    which results in:

    x_d - <x> = ( 1 - delta_t / tau ) * ( x - <x> ).

    Taking both sides to the third power:

    ( x_d - <x> )^3 = ( 1 - delta_t / tau )^3 * ( x - <x> )^3.

    After averaging both sides, the damped value of xp3 is:

    < x_d'^3 > = ( 1 - delta_t / tau )^3 * < x'^3 >.

    Any sponge damping is applied to (predictive) xp3 after the field has been
    advanced in time.  This allows sponge damping to be applied in an implicit
    manner.

    Returns:
        xp3_damped

    Fortran declarations (sponge_damp_xp3 in src/CLUBB_core/sponge_layer_damping.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          ngrdcol, & ! Number of columns [-]
          nzm,     & ! Number of momentum levels [-]
          nzt        ! Number of thermodynamic levels [-]

        type( grid ), intent(in) :: &
          gr    ! Grid bounds and vertical direction

        real( kind = core_rknd ), intent(in) :: &
          dt    ! Model Timestep  [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          z    ! Thermodynamic grid heights [m]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          zm    ! Momentum grid heights [m]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          xp3    ! Third moment before damping [units vary]

        type(sponge_damp_profile), intent(in) :: &
          damping_profile    ! Column time scales and layer depths

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
          xp3_damped    ! <x'^3> after damping   [units vary]
    """
    if damping_profile.tau_sponge_damp is None:
        raise ValueError("tau_sponge_damp in sponge_damp_xp3 used before initialization")
    levels = gr.k_ub_zt - gr.grid_dir_indx * jnp.arange(nzt)
    in_sponge = (zm[:, gr.k_ub_zm, None] - z[:, levels]) < damping_profile.sponge_layer_depth[:, None]
    l_damp = jnp.cumsum(~in_sponge, axis=1) == 0
    dt_on_tau = dt / damping_profile.tau_sponge_damp[:, levels]
    damped = (1.0 - dt_on_tau)**3 * xp3[:, levels]
    xp3_damped = xp3.at[:, levels].set(jnp.where(l_damp, damped, xp3[:, levels]))
    return xp3_damped


def initialize_tau_sponge_damp_api(ngrdcol, nz, gr, dt, z, settings):
    """Initialize the damping time scale independently in each column.  The time scale
    attains its maximum value, tau_sponge_damp_max, at the bottom of the
    "sponge" damping layer, which results in minimal damping.  Likewise, the
    time scale attains its minimum value, tau_sponge_damp_min, at the top of
    the model, which results in maximum damping.  At levels in-between the top
    of the model and the base of the sponge damping layer, the value of
    tau_sponge_damp is in-between tau_sponge_damp_min and tau_sponge_damp_max,
    as calculated by an interpolation formula.

    Returns:
        damping_profile

    Fortran declarations (initialize_tau_sponge_damp_api in src/CLUBB_core/sponge_layer_damping.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          ngrdcol, & ! Number of columns [-]
          nz         ! Number of vertical levels [-]

        type (grid), intent(in) :: gr    ! Model grid

        real( kind = core_rknd ), intent(in) :: &
          dt    ! Model Timestep    [s]

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          z     ! Height of model grid levels    [m]

        type(sponge_damp_settings), intent(in) :: &
          settings    ! Shared damping configuration

        !--------------------------- Output Variables ---------------------------
        type(sponge_damp_profile), intent(out) :: &
          damping_profile    ! Column time scales and layer depths
    """
    if settings.tau_sponge_damp_min < 2.0 * dt:
        raise ValueError("tau_sponge_damp_min is too small (< 2*dt)")
    zm_top = jnp.asarray(gr.zm[:, gr.k_ub_zm])
    sponge_layer_depth = settings.sponge_damp_depth * zm_top
    levels = gr.k_ub_zm - gr.grid_dir_indx * jnp.arange(nz)
    if nz != gr.nzm:
        levels = gr.k_ub_zt - gr.grid_dir_indx * jnp.arange(nz)
    distance = zm_top[:, None] - jnp.asarray(z)[:, levels]
    l_damp = jnp.cumsum(~(distance < sponge_layer_depth[:, None]), axis=1) == 0
    safe_depth = jnp.where(l_damp, sponge_layer_depth[:, None], 1.0)
    exponent = jnp.where(l_damp, distance, 0.0) / safe_depth
    tau_in_layer = settings.tau_sponge_damp_min * (
        settings.tau_sponge_damp_max / settings.tau_sponge_damp_min
    )**exponent
    tau = jnp.full((ngrdcol, nz), settings.tau_sponge_damp_max)
    tau = tau.at[:, levels].set(jnp.where(l_damp, tau_in_layer, settings.tau_sponge_damp_max))
    damping_profile = sponge_damp_profile(tau, sponge_layer_depth)
    return damping_profile


def finalize_tau_sponge_damp_api(damping_profile):
    """Frees memory allocated in initialize_tau_sponge_damp_api

    Fortran declarations (finalize_tau_sponge_damp_api in src/CLUBB_core/sponge_layer_damping.F90):

        ! Input/Output Variable(s)
        type(sponge_damp_profile), intent(inout) :: &
          damping_profile ! Information for damping the profile
    """
    return damping_profile._replace(tau_sponge_damp=None, sponge_layer_depth=None)
