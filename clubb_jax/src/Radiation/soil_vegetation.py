"""JAX port of ``Radiation/soil_vegetation.F90``.

Description:
    This subroutine updates the surface and soil temp, soil heat flux,
    while assuming that net radiation and turbulent heat fluxes are
    already available from another subroutine.
"""

from functools import partial

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv, kappa, p0, stefan_boltzmann

configure_jax_precision()


@partial(jax.jit, static_argnames=("ngrdcol",))
def advance_soil_veg(
    ngrdcol,
    dt,  # Current model timestep (Must be < 60s) [s]
    rho_sfc,  # Air density at the surface   [kg/m^3]
    Frad_SW_up_sfc,  # SW Net                       [W/m^3]
    Frad_SW_down_sfc,  # SW radiative upwelling flux  [W/m^2]
    Frad_LW_down_sfc,  # LW downwelling flux          [W/m^2]
    wpthlp_sfc, wprtp_sfc, p_sfc, stats,
    deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K,
):
    """Advance the interactive soil and vegetation temperatures one timestep.

    The surface temperature (sfc_soil_T_in_K) is calculated  from
    the surface energy budget.

    ************************************** 2
              *                   *
              *                   *
            --*         +         *--  1
              *                   *
    *** surface ************+************* 0
                       sfc_soil_T_in_K              l (evel)

    Heat conduction (in a homogeneous medium) can be described by
    the equation:

             d(sfc_soil_T_in_K)/dt=ks*d(d(sfc_soil_T_in_K)/dz)/dz

    In which ks is the soil thermal diffusity.
    We consider a semi half infinite medium, initially at the
    constant temperature deep_soil_T_in_K. if we vary the surface temperature
    sfc_soil_T_in_K sinussodally in time we can deduce a relation between the
    surface temperature, the soil heat flux (soil_heat_flux)  and deep_soil_T_in_K:

             d(sfc_soil_T_in_K)/dt=c1*soil_heat_flux - c2*(sfc_soil_T_in_K-deep_soil_T_in_K)

    However in reality the temperature deep_soil_T_in_K also varies in time, it
    may be calculated from:

             d(deep_soil_T_in_K)/dt= c2*soil_heat_flux

    The equations given above are analogous to those used by
    Deardorff  (1978).

    Reference:
    Duynkerke, Peter G.  "Radiation Fog: A Comparison of Model
    Simulation with Detailed Observations"
    (February 1991) _Monthly Weather Review_ Vol 119, p. 324-341.

    This subroutine does not produce any output variables. Instead the module
    variables listed below are updated.

    veg_T_in_K_in, &            ! Temperature of vegetation layer [K]
    sfc_soil_T_in_K_in, &       ! Temperature of surface soil layer [K]
    deep_soil_T_in_K_in         ! Temperature of deep soil layer [K]

    Returns Fortran inout values in order: ``stats``, ``deep_soil_T_in_K``,
    ``sfc_soil_T_in_K``, and ``veg_T_in_K``.
    """
    # External

    # Input variables

    # Local variables
    # --------------------------------- Begin Code ---------------------------------

    # ----------------------------
    #  Soil parameters
    # ---------------------------
    # soil heat capacity              [Jg/K]
    cs = 2.00e3  # cs
    # soil density                    [g/m^3]
    rs = 1.00e3  # ps
    # soil heat diffusivity           [m^2/s]
    ks = 2.00e-7  # as
    d1 = jnp.sqrt(ks * 3600.0 * 24.0)  # Known magic number
    c1 = 2.0 * jnp.sqrt(jnp.pi) / (rs * cs * d1)  # coefficient in force restore 1
    c2 = (  # coefficient in force restore 2
        2.0 * jnp.pi / (3600.0 * 24.0)  # Omega - known magic number
    )
    # coefficient in force restore 3
    c3 = jnp.sqrt(jnp.pi * 2.0) / (
        jnp.exp(jnp.pi / 4.0) * rs * cs * jnp.sqrt(ks * 3600.0 * 24.0 * 365.0)
    )  # Known magic number

    if stats.l_sample:
        stats = stats.update("veg_T_in_K", veg_T_in_K)
        stats = stats.update("sfc_soil_T_in_K", sfc_soil_T_in_K)
        stats = stats.update("deep_soil_T_in_K", deep_soil_T_in_K)

    Frad_LW_up_sfc = stefan_boltzmann * veg_T_in_K ** 4  # LW upwelling flux [W/m2]
    # Turbulent Flux of equivalent potential temperature   [K]
    wpthep = wpthlp_sfc + (Lv / Cp) * (p0 / p_sfc) ** kappa * wprtp_sfc

    # Calculate net radiation minus turbulent heat flux
    veg_heat_flux = (
        Frad_LW_down_sfc - Frad_LW_up_sfc - wpthep * rho_sfc * Cp
        + (Frad_SW_down_sfc - Frad_SW_up_sfc)
    )

    # Calculate soil heat flux
    # Duynkerke (1991) used a coefficient of 3.0 W/m^2*K, not 10.0 W/m^2*K
    #
    # Equation 19 p.328
    soil_heat_flux = (  # Soil Heat flux [W/m^2]
        10.0 * (veg_T_in_K - sfc_soil_T_in_K)
        + 0.05 * Frad_SW_down_sfc  # Known magic number
    )

    # Update surf veg temp
    veg_T_in_K = (
        veg_T_in_K
        + dt * 5.0e-5 * (veg_heat_flux - soil_heat_flux)  # Known magic number
    )

    # Update soil temp
    sfc_soil_T_in_K = sfc_soil_T_in_K + dt * (
        c1 * soil_heat_flux - c2 * (sfc_soil_T_in_K - deep_soil_T_in_K)
    )

    # Update deep soil temp
    deep_soil_T_in_K = deep_soil_T_in_K + dt * c3 * soil_heat_flux

    if stats.l_sample:
        # The source makes one scalar update per column. JaxStats accepts the
        # equivalent all-column surface vector in one call.
        stats = stats.update("soil_heat_flux", soil_heat_flux)

    return stats, deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K


def initialize_soil_veg(ngrdcol):
    """Sets some default values for the soil scheme.

    References:
        None

    Returns source outputs ``deep_soil_T_in_K``, ``sfc_soil_T_in_K``, and
    ``veg_T_in_K`` in that order.
    """
    # ---- Begin Code ----
    # These default values are the values for gabls3
    deep_soil_T_in_K = jnp.full((ngrdcol,), 288.58, dtype=jnp.float64)
    sfc_soil_T_in_K = jnp.full((ngrdcol,), 300.0, dtype=jnp.float64)
    veg_T_in_K = jnp.full((ngrdcol,), 300.0, dtype=jnp.float64)
    # Disable this for most cases
    # JAX adaptation: the source resets module ``l_soil_veg`` here. Its
    # immutable equivalent is initialized later from the case namelist.
    return deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K


__all__ = ["advance_soil_veg", "initialize_soil_veg"]
