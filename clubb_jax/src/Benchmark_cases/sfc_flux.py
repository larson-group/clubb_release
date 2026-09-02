"""This module contains generalized subroutines for determining surface
fluxes.

References:
    None.
"""

from __future__ import annotations

from functools import partial

import jax
import jax.numpy as jnp
from jax import lax

from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv
from clubb_jax.src.CLUBB_core.interpolation import linear_interp_factor


# ===============================================================================
@partial(jax.jit, static_argnames=("ngrdcol",))
def compute_momentum_flux(ngrdcol, um_sfc, vm_sfc, ubar, ustar):
    """This subroutine computes the momentum fluxes upwp_sfc and vpwp_sfc.

    References:
        None.
    """

    # ---- Begin Code ----
    # Compute momentum fluxes
    upwp_sfc = -um_sfc * ustar**2 / ubar
    vpwp_sfc = -vm_sfc * ustar**2 / ubar
    return upwp_sfc, vpwp_sfc


# ===============================================================================
@partial(jax.jit, static_argnames=("ngrdcol",))
def compute_ubar(ngrdcol, um_sfc, vm_sfc):
    """This function determines the value of ubar based on the momentum at
    the surface.

    References:
        None.
    """

    # Constant paramter(s)
    ubmin = 0.25

    # ---- Begin Code ----
    return jnp.maximum(ubmin, jnp.sqrt(um_sfc**2 + vm_sfc**2))


# ===============================================================================
@partial(jax.jit, static_argnames=("ntimes",))
def compute_ht_mostr_flux(time_in, ntimes):
    """Compute heat and moisture fluxes from ARM data in (W/m2).

    References:
        None.
    """
    from clubb_jax.src.Benchmark_cases import time_dependent_input

    time_sfc_given = time_dependent_input.time_sfc_given
    sens_ht_given = time_dependent_input.sens_ht_given
    latent_ht_given = time_dependent_input.latent_ht_given
    # ---------------------BEGIN CODE-------------------------
    # Default Initialization
    heat_flx = jnp.asarray(0.0, dtype=sens_ht_given.dtype)
    moisture_flx = jnp.asarray(0.0, dtype=latent_ht_given.dtype)
    time = jnp.asarray(time_in, dtype=sens_ht_given.dtype)

    # Compute heat and moisture fluxes from ARM data in (W/m2)
    # TODO(port-mirror): These local callables only express the source IF and
    # DO WHILE blocks through JAX control flow; remove them if JAX gains direct
    # traced Python control flow.
    def before_first(_):
        return sens_ht_given[0], latent_ht_given[0]

    def after_last(_):
        return sens_ht_given[ntimes - 1], latent_ht_given[ntimes - 1]

    def between_times(_):
        def loop_condition(carry):
            return carry[0] <= ntimes - 2

        def loop_body(carry):
            i1, heat_flx, moisture_flx = carry
            i2 = i1 + 1
            in_interval = (time >= time_sfc_given[i1]) & (time < time_sfc_given[i2])
            time_frac = (time - time_sfc_given[i1]) / (
                time_sfc_given[i2] - time_sfc_given[i1]
            )
            heat_flx = jnp.where(
                in_interval,
                linear_interp_factor(time_frac, sens_ht_given[i2], sens_ht_given[i1]),
                heat_flx,
            )
            moisture_flx = jnp.where(
                in_interval,
                linear_interp_factor(
                    time_frac,
                    latent_ht_given[i2],
                    latent_ht_given[i1],
                ),
                moisture_flx,
            )
            return i2, heat_flx, moisture_flx

        _, heat_flx_result, moisture_flx_result = lax.while_loop(
            loop_condition,
            loop_body,
            (0, heat_flx, moisture_flx),
        )
        return heat_flx_result, moisture_flx_result

    return lax.cond(
        time <= time_sfc_given[0],
        before_first,
        lambda _: lax.cond(
            time >= time_sfc_given[ntimes - 1],
            after_last,
            between_times,
            operand=None,
        ),
        operand=None,
    )


# ===============================================================================
@partial(jax.jit, static_argnames=("ngrdcol",))
def compute_wpthlp_sfc(ngrdcol, Cd, ubar, thlm_sfc, T_sfc, exner_sfc):
    """This function determins the surface flux of heat.

    References:
        None.
    """

    # ---- Begin Code ----
    return -Cd * ubar * (thlm_sfc - T_sfc / exner_sfc)


# ===============================================================================
@partial(jax.jit, static_argnames=("ngrdcol",))
def compute_wprtp_sfc(ngrdcol, Cd, ubar, rtm_sfc, adjustment):
    """This function determines the surface flux of moisture.

    References:
        None.
    """

    # ---- Begin Code ----
    return -Cd * ubar * (rtm_sfc - adjustment)


# ===============================================================================
def set_sclr_sfc_rtm_thlm(
    ngrdcol,
    sclr_dim,
    edsclr_dim,
    sclr_idx,
    wpthlp_sfc,
    wprtp_sfc,
):
    """This function determines the surface flux of moisture.

    References:
        None.
    """
    # --------------------- Begin Code ---------------------
    wpsclrp_sfc = jnp.zeros((ngrdcol, sclr_dim), dtype=wpthlp_sfc.dtype)
    wpedsclrp_sfc = jnp.zeros((ngrdcol, edsclr_dim), dtype=wpthlp_sfc.dtype)

    if sclr_dim > 0:
        # Let passive scalars be equal to rt and theta_l for now
        if sclr_idx.iisclr_thl > 0:
            wpsclrp_sfc = wpsclrp_sfc.at[:, sclr_idx.iisclr_thl - 1].set(wpthlp_sfc)
        if sclr_idx.iisclr_rt > 0:
            wpsclrp_sfc = wpsclrp_sfc.at[:, sclr_idx.iisclr_rt - 1].set(wprtp_sfc)

        if sclr_idx.iiedsclr_thl > 0:
            wpedsclrp_sfc = wpedsclrp_sfc.at[:, sclr_idx.iiedsclr_thl - 1].set(wpthlp_sfc)
        if sclr_idx.iiedsclr_rt > 0:
            wpedsclrp_sfc = wpedsclrp_sfc.at[:, sclr_idx.iiedsclr_rt - 1].set(wprtp_sfc)

    if edsclr_dim > 0:
        wpedsclrp_sfc = jnp.zeros((ngrdcol, edsclr_dim), dtype=wpthlp_sfc.dtype)

    return wpsclrp_sfc, wpedsclrp_sfc


# ===============================================================================
@partial(jax.jit)
def convert_sens_ht_to_km_s(sens_ht, rho_sfc):
    """This function converts sensible heat flux in W/m^2 to
    natural units of k m/s for the wpthlp_sfc variable.
    """
    # --------------------BEGIN CODE-----------------------
    return sens_ht / (rho_sfc * Cp)


# ===============================================================================
@partial(jax.jit)
def convert_latent_ht_to_m_s(latent_ht, rho_sfc):
    """This function converts latent heat flux in W/m^2 to
    natural units of m/s for the wprtp_sfc variable.
    """
    # --------------------BEGIN CODE-----------------------
    return latent_ht / (rho_sfc * Lv)


__all__ = [
    "compute_momentum_flux",
    "compute_ubar",
    "compute_ht_mostr_flux",
    "compute_wprtp_sfc",
    "compute_wpthlp_sfc",
    "set_sclr_sfc_rtm_thlm",
    "convert_sens_ht_to_km_s",
    "convert_latent_ht_to_m_s",
]
