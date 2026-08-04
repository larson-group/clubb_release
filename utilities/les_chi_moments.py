"""Derive linearized LES chi moments using standalone CLUBB thermodynamics."""

from __future__ import annotations

import numpy as np


# These match standalone CLUBB's constants and its default Flatau liquid
# saturation transform in saturation.F90 and pdf_closure_module.F90.
CP_DRY_AIR = 1004.67
LATENT_HEAT_VAPORIZATION = 2.5e6
RD_DRY_AIR = 287.04
RV_WATER_VAPOR = 461.5
REFERENCE_PRESSURE_PA = 1.0e5
FREEZING_TEMPERATURE_K = 273.15
EPSILON = RD_DRY_AIR / RV_WATER_VAPOR
KAPPA = RD_DRY_AIR / CP_DRY_AIR


def flatau_saturation_mixing_ratio(pressure_pa, temperature_k):
    """Reproduce CLUBB's default liquid saturation-mixing-ratio formula."""
    pressure_pa = np.asarray(pressure_pa, dtype=float)
    temperature_k = np.asarray(temperature_k, dtype=float)
    temperature_c = np.maximum(temperature_k - FREEZING_TEMPERATURE_K, -85.0)
    temperature_c_squared = np.square(temperature_c)
    saturation_vapor_pressure = (
        -3.21582393e-14
        * (temperature_c - 646.5835252598777)
        * (temperature_c + 90.72381630364440)
        * (temperature_c_squared + 111.0976961559954 * temperature_c + 6459.629194243118)
        * (temperature_c_squared + 152.3131930092453 * temperature_c + 6499.774954705265)
        * (temperature_c_squared + 174.4279584934021 * temperature_c + 7721.679732114084)
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        mixing_ratio = EPSILON * saturation_vapor_pressure / (
            pressure_pa - saturation_vapor_pressure
        )
    return np.where(pressure_pa - saturation_vapor_pressure < 1.0, EPSILON, mixing_ratio)


def derive_chi_moments(
    mean_rt,
    mean_thl,
    pressure_pa,
    var_rt,
    var_thl,
    covar_rt_thl,
    covar_w_rt,
    covar_w_thl,
    var_w,
):
    """Return mean/variance/covariance for chi'=crt*rt'-cthl*thl'."""
    arrays = np.broadcast_arrays(
        *[
            np.asarray(value, dtype=float)
            for value in (
                mean_rt,
                mean_thl,
                pressure_pa,
                var_rt,
                var_thl,
                covar_rt_thl,
                covar_w_rt,
                covar_w_thl,
                var_w,
            )
        ]
    )
    (
        mean_rt,
        mean_thl,
        pressure_pa,
        var_rt,
        var_thl,
        covar_rt_thl,
        covar_w_rt,
        covar_w_thl,
        var_w,
    ) = arrays
    valid = np.logical_and.reduce(
        [np.isfinite(value) for value in arrays]
        + [pressure_pa > 0.0, mean_thl > 0.0]
    )

    with np.errstate(all="ignore"):
        exner = np.power(pressure_pa / REFERENCE_PRESSURE_PA, KAPPA)
        liquid_water_temperature = mean_thl * exner
        rsatl = flatau_saturation_mixing_ratio(pressure_pa, liquid_water_temperature)
        beta = (
            EPSILON
            * LATENT_HEAT_VAPORIZATION**2
            / (RD_DRY_AIR * CP_DRY_AIR * np.square(liquid_water_temperature))
        )
        inverse_denominator = 1.0 / (1.0 + beta * rsatl)
        mean_chi = (mean_rt - rsatl) * inverse_denominator
        coef_rt = inverse_denominator
        coef_thl = (
            (1.0 + beta * mean_rt)
            * np.square(inverse_denominator)
            * (CP_DRY_AIR / LATENT_HEAT_VAPORIZATION)
            * beta
            * rsatl
            * exner
        )
        var_chi_raw = (
            np.square(coef_rt) * var_rt
            - 2.0 * coef_rt * coef_thl * covar_rt_thl
            + np.square(coef_thl) * var_thl
        )
        covar_w_chi_raw = coef_rt * covar_w_rt - coef_thl * covar_w_thl
        var_chi = np.maximum(var_chi_raw, 0.0)
        nonnegative_var_w = np.maximum(var_w, 0.0)
        covariance_limit = np.sqrt(var_chi * nonnegative_var_w)
        covar_w_chi = np.clip(covar_w_chi_raw, -covariance_limit, covariance_limit)

    finite_results = (
        np.isfinite(mean_chi)
        & np.isfinite(var_chi)
        & np.isfinite(covar_w_chi)
        & np.isfinite(coef_rt)
        & np.isfinite(coef_thl)
        & np.isfinite(rsatl)
    )
    valid &= finite_results
    adjusted = valid & (
        ~np.isclose(var_chi, var_chi_raw, rtol=1.0e-10, atol=1.0e-20)
        | ~np.isclose(covar_w_chi, covar_w_chi_raw, rtol=1.0e-10, atol=1.0e-20)
    )

    def masked(values):
        return np.where(valid, values, np.nan)

    return {
        "mean_chi": masked(mean_chi),
        "var_chi": masked(var_chi),
        "covar_w_chi": masked(covar_w_chi),
        "coef_rt": masked(coef_rt),
        "coef_thl": masked(coef_thl),
        "rsatl": masked(rsatl),
        "adjusted": adjusted,
    }
