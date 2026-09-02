"""Contains subroutines for the GABLS3 LES.

References:
    http://www4.ncsu.edu/~sbasu5/GABLS3/
"""

from __future__ import annotations

from functools import partial

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import ep1
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor
from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.sfc_flux import compute_momentum_flux


# -----------------------------------------------------------------------
def gabls3_night_sfclyr(
    ngrdcol,
    time,
    um_sfc,
    vm_sfc,
    thlm_sfc,
    rtm_sfc,
    lowest_level,
) -> tuple:
    """This subroutine computes surface fluxes of horizontal momentum,
    heat and moisture according to GCSS ATEX specifications.

    References:
        http://www4.ncsu.edu/~sbasu5/GABLS3/
    """
    z0 = 0.15

    if not time_dependent_input.l_t_dependent:
        return None

    # Use time_select to determine the time indexes before and after 'time', as well as
    # the time fraction necessary for linear_interp_factor
    before_time, after_time, time_frac = time_dependent_input.time_select(
        time,
        time_dependent_input.time_sfc_given.size,
        time_dependent_input.time_sfc_given,
    )
    ts = linear_interp_factor(
        time_frac,
        time_dependent_input.thlm_sfc_given[after_time],
        time_dependent_input.thlm_sfc_given[before_time],
    )
    qs = linear_interp_factor(
        time_frac,
        time_dependent_input.rtm_sfc_given[after_time],
        time_dependent_input.rtm_sfc_given[before_time],
    )

    # Compute heat and moisture fluxes
    # TODO(port-mirror): vmap expresses the source's column loop around the
    # scalar landflx routine without changing landflx's source-level interface.
    wpthlp_sfc, wprtp_sfc, ubar, ustar = jax.vmap(
        landflx,
        in_axes=(0, None, 0, None, 0, 0, 0, None),
    )(thlm_sfc, ts, rtm_sfc, qs, um_sfc, vm_sfc, lowest_level, z0)

    if time_dependent_input.l_input_xpwp_sfc:
        # Feed in momentum fluxes
        upwp_sfc_interp = linear_interp_factor(
            time_frac,
            time_dependent_input.upwp_sfc_given[after_time],
            time_dependent_input.upwp_sfc_given[before_time],
        )
        vpwp_sfc_interp = linear_interp_factor(
            time_frac,
            time_dependent_input.vpwp_sfc_given[after_time],
            time_dependent_input.vpwp_sfc_given[before_time],
        )
        upwp_sfc = jnp.full(
            ngrdcol,
            upwp_sfc_interp,
        )
        vpwp_sfc = jnp.full(
            ngrdcol,
            vpwp_sfc_interp,
        )
    else:
        # Compute momentum fluxes
        upwp_sfc, vpwp_sfc = compute_momentum_flux(
            ngrdcol, um_sfc, vm_sfc, ubar, ustar,
        )

    return upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, ustar


@partial(jax.jit)
def psi_h(x, xlmo):
    return (-5.0 * x) / xlmo


@partial(jax.jit)
def gm1(x):
    return (1.0 - 15.0 * x) ** 0.25


@partial(jax.jit)
def gh1(x):
    return jnp.sqrt(1.0 - 9.0 * x) / 0.74


@partial(jax.jit)
def fm1(x):
    # The source converts through default real inside alog. Preserve that
    # precision boundary explicitly before returning to the core dtype.
    core_dtype = jnp.result_type(x)
    x_real = jnp.asarray(x, dtype=jnp.float32)
    x_squared_real = jnp.asarray(x * x, dtype=jnp.float32)
    # XLA's float32 log may differ from host ALOG by one ulp. Evaluate in
    # float64, then explicitly round through the source's default real.
    log_one = jnp.log(
        ((1.0 + x_real) / 2.0).astype(jnp.float64)
    ).astype(jnp.float32).astype(core_dtype)
    log_two = jnp.log(
        ((1.0 + x_squared_real) / 2.0).astype(jnp.float64)
    ).astype(jnp.float32).astype(core_dtype)
    pii = jnp.arccos(jnp.asarray(-1.0, dtype=core_dtype)) / 2.0
    return 2.0 * log_one + log_two - 2.0 * jnp.arctan(x) + pii


@partial(jax.jit)
def fh1(x):
    # The source converts through default real inside alog. Preserve that
    # precision boundary explicitly before returning to the core dtype.
    core_dtype = jnp.result_type(x)
    x_real = jnp.asarray(x, dtype=jnp.float32)
    log_arg = (1.0 + jnp.float32(0.74) * x_real) / 2.0
    return 2.0 * jnp.log(log_arg.astype(jnp.float64)).astype(jnp.float32).astype(core_dtype)


@partial(jax.jit)
def landflx(th, ts, qh, qs, uh, vh, h, z0):
    """landflx.F90 from SAM 6.7.5."""
    core_dtype = jnp.result_type(th, ts, qh, qs, uh, vh, h, z0)
    zody_arg = (
        jnp.asarray(h, dtype=jnp.float32) / jnp.asarray(z0, dtype=jnp.float32)
    )
    zody = jnp.log(zody_arg.astype(jnp.float64)).astype(jnp.float32).astype(core_dtype)

    vel = jnp.sqrt(jnp.maximum(0.5, uh ** 2 + vh ** 2))
    r = (
        9.81 / ts
        * (th * (1.0 + ep1 * qh) - ts * (1.0 + ep1 * qs))
        * h
        / vel ** 2
    )

    # TODO(port-mirror): lax.cond requires callable branches; these two local
    # functions preserve the source r < 0 branch without becoming module APIs.
    def unstable_branch(_):
        xsi = jnp.asarray(0.0, dtype=core_dtype)

        xm = gm1(xsi)
        xh = gh1(xsi)
        fm = zody - fm1(xm)
        fh = 0.74 * (zody - fh1(xh))
        xsi = r / fh * fm ** 2

        xsi = -jnp.abs(xsi)
        xm = gm1(xsi)
        xh = gh1(xsi)
        fm = zody - fm1(xm)
        fh = 0.74 * (zody - fh1(xh))
        xsi = r / fh * fm ** 2

        xsi = -jnp.abs(xsi)
        xm = gm1(xsi)
        xh = gh1(xsi)
        fm = zody - fm1(xm)
        fh = 0.74 * (zody - fh1(xh))
        xsi = r / fh * fm ** 2
        return xsi, fm, fh

    def stable_branch(_):
        a = 4.8 * 4.8 * r - 1.00 * 6.35
        b = (2.0 * r * 4.8 - 1.00) * zody
        c = r * zody ** 2
        d = jnp.sqrt(b * b - 4.0 * a * c)
        xsi1 = (-b + d) / a / 2.0
        xsi2 = (-b - d) / a / 2.0
        xsi = jnp.maximum(
            jnp.asarray(xsi1, dtype=jnp.float32),
            jnp.asarray(xsi2, dtype=jnp.float32),
        ).astype(core_dtype)
        fm = zody + 4.8 * xsi
        fh = zody + 7.8 * xsi
        return xsi, fm, fh

    xsi, fm, fh = jax.lax.cond(r < 0.0, unstable_branch, stable_branch, None)

    vel = jnp.sqrt(uh ** 2 + vh ** 2)

    # Modification for GABLS3_night
    # Specification states how to compute these
    # Joshua Fasching January 2009
    ustar = 0.4 / fm * vel

    xsi = jax.lax.cond(
        xsi >= 0.0,
        lambda value: jnp.maximum(1.0e-5, value),
        lambda value: jnp.minimum(-1.0e-5, value),
        xsi,
    )
    xlmo = h / xsi
    #------ Modification of GABLS3_night
    surface_log_arg = jnp.asarray(h, dtype=jnp.float32) / jnp.float32(0.25)
    surface_log = jnp.log(
        surface_log_arg.astype(jnp.float64)
    ).astype(jnp.float32).astype(core_dtype)
    denominator = surface_log - psi_h(h, xlmo) + psi_h(0.25, xlmo)
    shf = 0.4 * ustar * (ts - th) / denominator
    lhf = 0.4 * ustar * (qs - qh) / denominator
    #-----------------------------
    return shf, lhf, vel, ustar


__all__ = [
    "gabls3_night_sfclyr",
]
