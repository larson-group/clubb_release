"""Kessler sampling diagnostics from est_kessler_microphys_module.F90."""

import jax
import jax.numpy as jnp
from jax.scipy.special import erf
from clubb_jax.src.CLUBB_core.constants_clubb import chi_tol
from clubb_jax.src.CLUBB_core.advance_helper_module import sqrt_clipped


# -----------------------------------------------------------------------------
def est_kessler_microphys_api(
    nzt, num_samples, pdf_dim, ngrdcol,             # In
    X_nl_all_levs, pdf_params, rcm, cloud_frac,     # In
    X_mixt_comp_all_levs, lh_sample_point_weights,  # In
    l_lh_importance_sampling,                       # In
):
    """This subroutine computes microphysical grid box averages of the
    Kessler autoconversion scheme, using both Latin hypercube sampling
    and analytic integration, given a Latin Hypercube sample.

    X_nl_all_levs has shape (ngrdcol, num_samples, nzt, pdf_dim).
    Return lh_AKm, AKm, AKstd, AKstd_cld, AKm_rcm, AKm_rcc [kg/kg/s],
    followed by lh_rcm_avg [kg/kg], all on thermodynamic levels, and the
    per-column source ERROR STOP status for the host's ErrInfo.

    Arguments:
        nzt: Number of vertical levels
        num_samples: Number of sample points
        pdf_dim: Number of variates
        ngrdcol: Number of model columns
        X_nl_all_levs: Sample that is transformed ultimately to normal-lognormal
        pdf_params: PDF parameters [units vary]
        rcm: Liquid water mixing ratio [kg/kg]
        cloud_frac: Cloud fraction [-]
        X_mixt_comp_all_levs: Whether we're in mixture component 1 or 2
        lh_sample_point_weights: Weight for cloud weighted sampling
        l_lh_importance_sampling: Do importance sampling (SILHS) [-]
    """
    # -------------------------------------------------------------------------
    # Call Kessler autoconversion using the Latin-hypercube sample, then compute
    # its grid-box averages analytically. The first PDF variate is chi.
    # -------------------------------------------------------------------------
    rcm_sample = jnp.maximum(X_nl_all_levs[..., 0], 0.0)

    # Idealized autoconversion: A_K = coeff * max(rc - r_crit, 0).
    # The source uses coeff = 1.e-3 [1/s] and r_crit = 0.2e-3 [kg/kg].
    # The source uses unit component cloud fractions for these estimates.
    cloud_frac_1 = cloud_frac_2 = jnp.ones_like(cloud_frac)
    lh_AKm, l_error_AKm = calc_estimate(
        num_samples, pdf_params.mixt_frac,             # In
        cloud_frac_1, cloud_frac_2, rcm_sample,         # In
        X_mixt_comp_all_levs, lh_sample_point_weights,  # In
        l_lh_importance_sampling,                       # In
        1.0e-3, 0.2e-3,                                # In
    )

    # Monte Carlo estimate of liquid water, for comparison with analytic rcm.
    lh_rcm_avg, l_error_rcm = calc_estimate(
        num_samples, pdf_params.mixt_frac,             # In
        cloud_frac_1, cloud_frac_2, rcm_sample,         # In
        X_mixt_comp_all_levs, lh_sample_point_weights,  # In
        l_lh_importance_sampling,                       # In
        1.0, 0.0,                                      # In
    )

    # Exact Kessler autoconversion [kg/kg/s] in the two PDF components.
    r_crit = 0.2e-3
    K_one = 1.0e-3
    chi_n_1_crit = (pdf_params.chi_1 - r_crit) / jnp.maximum(pdf_params.stdev_chi_1, chi_tol)
    cloud_frac_1_crit = 0.5 * (1.0 + erf(chi_n_1_crit / jnp.sqrt(2.0)))
    AK1 = K_one * (
        (pdf_params.chi_1 - r_crit) * cloud_frac_1_crit
        + pdf_params.stdev_chi_1 * jnp.exp(-0.5 * chi_n_1_crit**2) / jnp.sqrt(2.0 * jnp.pi)
    )
    chi_n_2_crit = (pdf_params.chi_2 - r_crit) / jnp.maximum(pdf_params.stdev_chi_2, chi_tol)
    cloud_frac_2_crit = 0.5 * (1.0 + erf(chi_n_2_crit / jnp.sqrt(2.0)))
    AK2 = K_one * (
        (pdf_params.chi_2 - r_crit) * cloud_frac_2_crit
        + pdf_params.stdev_chi_2 * jnp.exp(-0.5 * chi_n_2_crit**2) / jnp.sqrt(2.0 * jnp.pi)
    )
    AKm = pdf_params.mixt_frac * AK1 + (1.0 - pdf_params.mixt_frac) * AK2

    # Exact Kessler standard deviation [kg/kg/s]. The component variances can
    # be slightly negative from roundoff, so the source clips them at zero.
    AK1var = jnp.maximum(
        0.0,
        K_one * (pdf_params.chi_1 - r_crit) * AK1
        + K_one * K_one * pdf_params.stdev_chi_1**2 * cloud_frac_1_crit
        - AK1**2,
    )
    AK2var = jnp.maximum(
        0.0,
        K_one * (pdf_params.chi_2 - r_crit) * AK2
        + K_one * K_one * pdf_params.stdev_chi_2**2 * cloud_frac_2_crit
        - AK2**2,
    )

    # Grid-box average standard deviation.
    AKstd = sqrt_clipped(
        pdf_params.mixt_frac * ((AK1 - AKm) ** 2 + AK1var)
        + (1.0 - pdf_params.mixt_frac) * ((AK2 - AKm) ** 2 + AK2var)
    )

    # Within-cloud standard deviation. JAX masks the denominator before the
    # calculation so the inactive cloud-free branch cannot divide by zero.
    safe_cf = jnp.where(cloud_frac > 0.0, cloud_frac, 1.0)
    AKstd_cld = jnp.where(
        cloud_frac > 0.0,
        sqrt_clipped(
            (
                pdf_params.mixt_frac * (AK1**2 + AK1var)
                + (1.0 - pdf_params.mixt_frac) * (AK2**2 + AK2var)
            )
            / safe_cf
            - (AKm / safe_cf) ** 2
        ),
        0.0,
    )

    # Local Kessler autoconversion using grid-box-average liquid water.
    AKm_rcm = K_one * jnp.maximum(0.0, rcm - r_crit)

    # Local Kessler autoconversion using within-cloud liquid water. The source
    # 0.001 cloud-fraction threshold avoids NaNs for very small cloud fractions
    # (D. Schanen, 3 June 2009).
    AKm_rcc = jnp.where(
        cloud_frac > 0.001,
        cloud_frac * K_one * jnp.maximum(0.0, rcm / safe_cf - r_crit),
        0.0,
    )
    l_error = jnp.any(l_error_AKm | l_error_rcm, axis=1)
    return lh_AKm, AKm, AKstd, AKstd_cld, AKm_rcm, AKm_rcc, lh_rcm_avg, l_error


# -----------------------------------------------------------------------------
def calc_estimate(
    num_samples, mixt_frac,                        # In
    cloud_frac_1, cloud_frac_2, rc,                # In
    X_mixt_comp_one_lev, lh_sample_point_weights,  # In
    l_lh_importance_sampling,                      # In
    coeff, r_crit,                                 # In
):
    """Compute the Monte Carlo estimate of idealized autoconversion.

    Source l_cloud_weighted_averaging is a local False constant. Both component
    sums are combined and divided by the total sample count. Return the estimate
    and local source ERROR STOP status so traced callers can update ErrInfo.

    Arguments:
        num_samples: Number of calls to microphysics (normally=2)
        mixt_frac: Mixture fraction of Gaussians
        cloud_frac_1: Cloud fraction associated w/ 1st, 2nd mixture component
        cloud_frac_2: Cloud fraction associated w/ 1st, 2nd mixture component
        rc: Sample specific liquid-water content (when positive) [kg/kg]
        X_mixt_comp_one_lev: Whether we're in the first or second mixture component
        lh_sample_point_weights: Weight for cloud weighted sampling
        l_lh_importance_sampling: Do importance sampling (SILHS)
        coeff: Autoconversion coefficient [1/s]
        r_crit: Cloud-water threshold for autoconversion [kg/kg]
    """
    # Handle possible errors in the ranges of mixt_frac and cloud_frac_1/2.
    l_error = (
        (mixt_frac > 1.0) | (mixt_frac < 0.0)
        | (cloud_frac_1 > 1.0) | (cloud_frac_1 < 0.0)
        | (cloud_frac_2 > 1.0) | (cloud_frac_2 < 0.0)
    )
    if num_samples == 0:
        # Source ERROR STOP: no sample points in calc_estimate.
        return jnp.zeros_like(mixt_frac), jnp.ones_like(mixt_frac, dtype=bool)

    # Initialize autoconversion in each mixture component. A scan preserves
    # the source sample-loop addition order; weights apply only for SILHS.
    weight = lh_sample_point_weights if l_lh_importance_sampling else 1.0
    sample_estimate = coeff * jnp.maximum(0.0, rc - r_crit) * weight

    # TODO(port-mirror): this scan body expresses the source sequential sample
    # loop; retain it until JAX can compile that loop without a callable body.
    def accumulate(est_m, sample):
        est_m1, est_m2 = est_m
        value, component = sample
        first = component == 1
        return (
            est_m1 + jnp.where(first, value, 0.0),
            est_m2 + jnp.where(first, 0.0, value),
        ), None

    (est_m1, est_m2), _ = jax.lax.scan(
        accumulate,
        (jnp.zeros_like(mixt_frac), jnp.zeros_like(mixt_frac)),
        (jnp.moveaxis(sample_estimate, 1, 0), jnp.moveaxis(X_mixt_comp_one_lev, 1, 0)),
    )
    return (est_m1 + est_m2) / num_samples, l_error
