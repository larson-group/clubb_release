"""JAX port of ``Radiation/cos_solar_zen_module.F90``.

Description:
    A function based on coefficients from Liou and the Clayson and
    Curry formula.  Approximates the cosine of the solar zenith
    angle anywhere in the world based on current Greenwich mean
    time and the latitude and longitude.

References:
    Clayson and Curry formula from C. A. Clayson and J. A. Curry ,
    J. Geophys.
    Res. Vol. 101, No. C12, Pages 28515-28528, 15 Dec. 1996.
    Liou ``An Introduction to Atmospheric Radiation''
      Table 2.2 and Eqn. 2.2.10
"""

from functools import partial

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
from clubb_jax.src.CLUBB_core.constants_clubb import radians_per_deg_dp, sec_per_hr
from clubb_jax.src.CLUBB_core.calendar import (
    compute_current_date_api,
    gregorian2julian_day,
    leap_year,
)

configure_jax_precision()

# Input Variables
@partial(jax.jit, static_argnames=("day", "month", "year"))
def cos_solar_zen(
    day,  # Day of month at model start
    month,  # Month of year at model start
    year,  # Year at model start
    current_time,  # Current time since start date [s]
    lat_in_degrees,  # Latitude       [degrees_N]
    lon_in_degrees,  # Longitude      [degrees_E]
):
    """Return the source function result ``cos_solar_zen``."""
    # External

    # Constant Parameters

    # Liou's coefficients
    c0 = 0.006918  # [-]
    c1 = -0.399912  # [-]
    c2 = -0.006758  # [-]
    c3 = -0.002697  # [-]
    d1 = 0.070257  # [-]
    d2 = 0.000907  # [-]
    d3 = 0.000148  # [-]

    # Local Variables
    present_day, present_month, present_year, present_time = compute_current_date_api(
        day, month, year, current_time
    )

    jul_day = gregorian2julian_day(present_day, present_month, present_year)
    days_in_year = jnp.where(leap_year(present_year), 366.0, 365.0)

    # Determine the number of hours
    hour = present_time / sec_per_hr
    t = 2.0 * jnp.pi * (jul_day - 1) / days_in_year

    delta = (
        c0
        + c1 * jnp.cos(t) + d1 * jnp.sin(t)
        + c2 * jnp.cos(2.0 * t) + d2 * jnp.sin(2.0 * t)
        + c3 * jnp.cos(3.0 * t) + d3 * jnp.sin(3.0 * t)
    )

    # The angle  longang  is equivalent to the
    # hour angle in the formula for cosZ .
    # References: Source file zenith.f
    #   from http://magic.gfdi.fsu.edu/seaflux/DIURNAL/README.txt
    #   Clayson and Curry formula from C. A. Clayson and J. A. Curry ,
    #   J. Geophys.
    #   Res. Vol. 101, No. C12, Pages 28515-28528, 15 Dec. 1996 .

    #   June 6, 2006

    h = jnp.floor(hour).astype(jnp.int32)
    l_first_half_day = (h >= 0) & (h <= 11)
    l_second_half_day = (h >= 12) & (h <= 23)
    # JAX adaptation: ``compute_current_date_api`` normalizes finite model time
    # to [0, 86400), so the source fatal default is unreachable. Preserve an
    # invalid value as NaN rather than silently selecting the second branch.
    zln = jnp.where(
        l_first_half_day,
        180.0 - hour * 15.0,  # Known magic number
        jnp.where(
            l_second_half_day,
            540.0 - hour * 15.0,  # Known magic number
            jnp.nan,
        ),
    )

    longang = jnp.abs(lon_in_degrees - zln) * radians_per_deg_dp
    latang = lat_in_degrees * radians_per_deg_dp

    # Return Variable
    # Cosine of the solar zenith angle (sometimes denoted amu0).
    cos_solar_zen_result = (
        jnp.sin(latang) * jnp.sin(delta)
        + jnp.cos(latang) * jnp.cos(delta) * jnp.cos(longang)
    )

    return cos_solar_zen_result


__all__ = ["cos_solar_zen"]
