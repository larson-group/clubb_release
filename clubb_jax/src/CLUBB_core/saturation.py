"""JAX implementations of selected routines from `src/CLUBB_core/saturation.F90`.

Description:
  Contains functions that compute saturation with respect
  to liquid or ice.

Porting deviations:
- Fortran exposes generic interfaces with scalar and 2D implementations.
  JAX uses broadcasting functions for both cases.
- `sat_mixrat_liq` mirrors `I_sat_sphum = .false.` because the JAX runtime does
  not currently expose GFDL specific-humidity saturation semantics.
- `saturation_lookup` depends on the Fortran lookup table and is intentionally
  not ported here.
- `sat_mixrat_ice` currently mirrors the Flatau ice polynomial path used by the
  JAX closure callers. The Fortran routine also has Bolton and GFDL branches.
- `rcm_sat_adj` uses fixed-count bisection with `jax.lax.scan` rather than the
  Fortran tolerance/`itermax` loop so it remains JIT compatible.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import (
    Cp,
    Lv,
    T_freeze_K,
    ep,
)
from clubb_jax.src.CLUBB_core.model_flags import (
    saturation_bolton,
    saturation_flatau,
    saturation_gfdl,
    saturation_lookup,
)

SATURATION_FLATAU = saturation_flatau
SATURATION_BOLTON = saturation_bolton
SATURATION_GFDL = saturation_gfdl
SATURATION_LOOKUP = saturation_lookup

_FLATAU_ICE_MIN_T_C = -90.0
_FLATAU_ICE_A = (
    100.0 * 6.09868993,
    100.0 * 0.499320233,
    100.0 * 0.184672631e-1,
    100.0 * 0.402737184e-3,
    100.0 * 0.565392987e-5,
    100.0 * 0.521693933e-7,
    100.0 * 0.307839583e-9,
    100.0 * 0.105785160e-11,
    100.0 * 0.161444444e-14,
)


def sat_vapor_press_liq_flatau(T_in_K):
    """Computes SVP for water vapor.

    References:
      ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
        and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
        pp. 1507--1513
    """
    # Determine deg K - 273.15
    T_in_C = jnp.maximum(T_in_K - T_freeze_K, -85.0)

    # Since this approximation is only good out to -85 degrees Celsius we
    # truncate the result here (Flatau, et al. 1992)
    T_in_C_sqd = T_in_C ** 2

    # Polynomial approx. (Flatau, et al. 1992)

    # Factoring the polynomial above and changing it into this form allows the cpu
    # to complete the calculations out of order. This is because modern cpus can complete
    # multiple instructions at once if they do not depend on eachother, in the above case
    # each instruction relies on the result of the last. In this version however, the terms
    # in the parentheses could potentially be calculated in parallel by different execution
    # units in the cpu, then only when those terms are being multiplied together do the
    # instructions need to be done one at a time. See clubb issue 834 for more info.
    #   - Gunther Huebler, Aug 2018
    return (
        -3.21582393e-14
        * (T_in_C - 646.5835252598777)
        * (T_in_C + 90.72381630364440)
        * (T_in_C_sqd + 111.0976961559954 * T_in_C + 6459.629194243118)
        * (T_in_C_sqd + 152.3131930092453 * T_in_C + 6499.774954705265)
        * (T_in_C_sqd + 174.4279584934021 * T_in_C + 7721.679732114084)
    )


def sat_vapor_press_liq_bolton(T_in_K):
    """Computes SVP for water vapor.

    References:
      Bolton 1980
    """
    # (Bolton 1980) approx.
    # Generally this more computationally expensive than the Flatau polnomial expansion
    return 611.2 * jnp.exp(17.67 * (T_in_K - T_freeze_K) / (T_in_K - 29.65))


def sat_vapor_press_liq_gfdl(T_in_K):
    """Compute saturation vapor pressure with respect to liquid.

    Description:
      copy from "GFDL polysvp.F90"
       Compute saturation vapor pressure with respect to liquid  by using
      function from Goff and Gratch (1946)

       Polysvp returned in units of pa.
       T_in_K  is input in units of K.

    JAX note:
      This helper is shared with `sat_mixrat_liq`; it keeps the lower
      temperature bound used in the Fortran `sat_mixrat_liq_*` GFDL branch.
    """
    # Since the Goff-Gratch approximation is valid only down to -70 degrees Celsius,
    #   we threshold the temperature.  This will yield a minimal saturation at
    #   cold temperatures.
    T_in_K_clipped = jnp.maximum(173.15, T_in_K)

    # Goff Gratch equation, uncertain below -70 C
    return (
        10.0 ** (
            -7.90298 * (373.16 / T_in_K_clipped - 1.0)
            + 5.02808 * jnp.log10(373.16 / T_in_K_clipped)
            - 1.3816e-7
            * (10.0 ** (11.344 * (1.0 - T_in_K_clipped / 373.16)) - 1.0)
            + 8.1328e-3
            * (10.0 ** (-3.49149 * (373.16 / T_in_K_clipped - 1.0)) - 1.0)
            + jnp.log10(1013.246)
        )
        * 100.0
    )


def sat_vapor_press_liq(T_in_K, saturation_formula: int):
    """Computes SVP for water vapor.

    Description:
      Computes SVP for water vapor. Calls one of the other functions
      that calculate an approximation to SVP.

    References:
      None
    """
    T_in_K = jnp.asarray(T_in_K, dtype=jnp.float64)
    # Saturation Vapor Pressure, esat, can be found to be approximated
    # in many different ways.
    if saturation_formula == saturation_flatau:
        # Using the Flatau, et al. polynomial approximation for SVP over vapor
        return sat_vapor_press_liq_flatau(T_in_K)
    if saturation_formula == saturation_bolton:
        # Using the Bolton 1980 approximations for SVP over vapor
        return sat_vapor_press_liq_bolton(T_in_K)
    if saturation_formula == saturation_gfdl:
        # Using GFDL polynomial approximation for SVP with respect to liquid
        return sat_vapor_press_liq_gfdl(T_in_K)
    if saturation_formula == saturation_lookup:
        # Use the lookup table to determine the saturation vapor pressure.
        raise ValueError("saturation_lookup is not ported to JAX yet.")
    # Undefined approximation
    raise ValueError(f"Unsupported saturation_formula={saturation_formula}")


def sat_mixrat_liq(p_in_Pa, T_in_K, saturation_formula: int):
    """Used to compute the saturation mixing ratio of liquid water.

    References:
      Formula from Emanuel 1994, 4.4.14
    """
    p_in_Pa = jnp.asarray(p_in_Pa, dtype=jnp.float64)
    T_in_K = jnp.asarray(T_in_K, dtype=jnp.float64)
    # Calculate the SVP for water vapor.
    esat = sat_vapor_press_liq(T_in_K, saturation_formula)

    # If esat exceeds the air pressure, then assume esat~=0.5*pressure
    #   and set rsat = ep = 0.622
    denom = p_in_Pa - esat
    safe = denom >= 1.0
    denom_safe = jnp.where(safe, denom, 1.0)
    # Formula for Saturation Mixing Ratio:
    #
    # rs = (epsilon) * [ esat / ( p - esat ) ];
    # where epsilon = R_d / R_v
    return jnp.where(safe, ep * esat / denom_safe, ep)


def sat_mixrat_ice(p_in_Pa, T_in_K):
    """Used to compute the saturation mixing ratio of ice.

    References:
      Formula from Emanuel 1994, 4.4.15
    """
    p_in_Pa = jnp.asarray(p_in_Pa, dtype=jnp.float64)
    T_in_K = jnp.asarray(T_in_K, dtype=jnp.float64)
    # Determine deg K - 273.15
    T_in_C = jnp.maximum(T_in_K - T_freeze_K, _FLATAU_ICE_MIN_T_C)
    a = _FLATAU_ICE_A
    # Using the Flatau, et al. polynomial approximation for SVP over ice
    # Since this approximation is only good out to -90 degrees Celsius we
    # truncate the result here (Flatau, et al. 1992)
    # Polynomial approx. (Flatau, et al. 1992)
    esat_ice = (
        a[0]
        + T_in_C
        * (
            a[1]
            + T_in_C
            * (
                a[2]
                + T_in_C
                * (
                    a[3]
                    + T_in_C
                    * (
                        a[4]
                        + T_in_C
                        * (
                            a[5]
                            + T_in_C
                            * (
                                a[6]
                                + T_in_C * (a[7] + T_in_C * a[8])
                            )
                        )
                    )
                )
            )
        )
    )
    # If esat_ice exceeds the air pressure, then assume esat_ice~=0.5*pressure
    #   and set rsat = ep = 0.622
    denom = p_in_Pa - esat_ice
    safe = denom >= 1.0
    denom_safe = jnp.where(safe, denom, 1.0)
    # Formula for Saturation Mixing Ratio:
    #
    # rs = (epsilon) * [ esat / ( p - esat ) ];
    # where epsilon = R_d / R_v
    return jnp.where(safe, ep * esat_ice / denom_safe, ep)


def rcm_sat_adj(thlm, rtm, p_in_Pa, exner, saturation_formula: int):
    """Find rcm from an initial profile that has saturation at some point.

    Description:

      This function uses an iterative method to find the value of rcm
      from an initial profile that has saturation at some point.

    References:
      None
    """
    thlm = jnp.asarray(thlm, dtype=jnp.float64)
    rtm = jnp.asarray(rtm, dtype=jnp.float64)
    p_in_Pa = jnp.asarray(p_in_Pa, dtype=jnp.float64)
    exner = jnp.asarray(exner, dtype=jnp.float64)

    tolerance = 0.001
    zero_threshold = 0.0

    theta = thlm
    too_low = jnp.zeros_like(theta)
    too_high = jnp.zeros_like(theta)
    done = jnp.zeros_like(theta, dtype=bool)

    def step(carry, iteration):
        theta, too_low, too_high, done = carry
        rsat = sat_mixrat_liq(p_in_Pa, theta * exner, saturation_formula)
        answer = theta - (Lv / (Cp * exner)) * jnp.maximum(
            rtm - rsat,
            zero_threshold,
        )
        diff = answer - thlm
        converged = jnp.abs(diff) <= tolerance
        active = ~done

        set_high = active & (~converged) & (diff > tolerance)
        set_low = active & (~converged) & ((-diff) > tolerance)
        too_high_new = jnp.where(set_high, theta, too_high)
        too_low_new = jnp.where(set_low, theta, too_low)

        too_high_new = jnp.where(
            active & (~converged) & (iteration == 0),
            theta + 20.0,
            too_high_new,
        )
        theta_new = 0.5 * (too_low_new + too_high_new)

        done_new = done | (active & converged)
        theta = jnp.where(active & (~converged), theta_new, theta)
        too_low = jnp.where(active, too_low_new, too_low)
        too_high = jnp.where(active, too_high_new, too_high)
        return (theta, too_low, too_high, done_new), None

    (theta, _too_low, _too_high, _done), _ = jax.lax.scan(
        step,
        (theta, too_low, too_high, done),
        jnp.arange(64),
    )
    rsat = sat_mixrat_liq(p_in_Pa, theta * exner, saturation_formula)
    return jnp.maximum(rtm - rsat, zero_threshold)
