"""Uniform-to-normal/lognormal transformations from transform_to_pdf_module.F90.

Fortran Sigma_Cholesky/std_normal argument layouts are retained; physical sample
arrays use (column, sample, level, variable). Output arguments become returns.
"""

import jax
import jax.numpy as jnp
from jax.scipy.special import erfc


# -----------------------------------------------------------------------------
def transform_uniform_samples_to_pdf(
    nzt, ngrdcol, num_samples, pdf_dim, d_uniform_extra,  # In
    hm_metadata,                                          # In
    Sigma_Cholesky1, Sigma_Cholesky2,                     # In
    mu1, mu2, X_mixt_comp_all_levs,                       # In
    X_u_all_levs, cloud_frac,                             # In
    l_in_precip_all_levs,                                 # In
):
    """Transform uniform samples to samples from CLUBB's PDF.

    Physical sample arrays use (column, sample, level, variable); the source
    standard-normal and Cholesky layouts are retained by the lower routines.

    Arguments:
        nzt: Number of vertical grid levels
        ngrdcol: Number of grid columns
        num_samples: Number of subcolumn samples
        pdf_dim: `d' Number of variates (normally 3 + microphysics specific variables)
        d_uniform_extra: Number of variates included in uniform sample only (often 2)
        hm_metadata: Hydrometeor/PDF variable index metadata
        Sigma_Cholesky1: Correlations Cholesky matrix, 1st component [-]
        Sigma_Cholesky2: Correlations Cholesky matrix, 2nd component [-]
        mu1: Means of the hydrometeors,(chi, eta, w, <hydrometeors>), 1st component [units
            vary]
        mu2: Means of the hydrometeors,(chi, eta, w, <hydrometeors>), 2nd component [units
            vary]
        X_mixt_comp_all_levs: Whether we're in the first or second mixture component
        X_u_all_levs: Sample drawn from uniform distribution from a particular grid level
        cloud_frac: Cloud fraction [-]
        l_in_precip_all_levs: Whether we are in precipitation (T/F)
    """
    # -------------------------------------------------------------------------
    # Generate sample points for a microphysics/radiation scheme.
    # -------------------------------------------------------------------------
    # From the Latin-hypercube sample, generate a standard-normal sample.
    std_normal = cdfnorminv(pdf_dim, nzt, ngrdcol, num_samples, X_u_all_levs[..., :pdf_dim])

    # Compute the nonstandard normal sample from Sigma's Cholesky factor,
    # the standard-normal sample, and the component means.
    X_nl_all_levs = multiply_Cholesky(
        nzt, ngrdcol, num_samples, pdf_dim, std_normal,  # In
        Sigma_Cholesky1, Sigma_Cholesky2,                # In
        mu1, mu2, X_mixt_comp_all_levs,                  # In
    )

    # Determine lognormal variables: chi, eta and w are normal; Ncn and the
    # hydrometeors that follow them are lognormal. Convert the latter with exp.
    p = max(hm_metadata.iiPDF_chi, hm_metadata.iiPDF_eta, hm_metadata.iiPDF_w) + 1
    X_nl_all_levs = X_nl_all_levs.at[..., p:].set(jnp.exp(X_nl_all_levs[..., p:]))

    # Zero precipitation hydrometeors outside precipitation; retain Ncn.
    p = hm_metadata.iiPDF_Ncn + 1
    X_nl_all_levs = X_nl_all_levs.at[..., p:].set(
        jnp.where(l_in_precip_all_levs[..., None], X_nl_all_levs[..., p:], 0.0)
    )

    # Clip extreme sample chi values. PDF closure clips component cloud fraction
    # under extreme conditions; force chi to the same saturated/unsaturated side.
    eps = jnp.finfo(jnp.float64).eps
    chi = X_nl_all_levs[..., hm_metadata.iiPDF_chi]
    chi = jnp.where(
        cloud_frac < eps,
        jnp.minimum(chi, 0.0),
        jnp.where(cloud_frac > 1.0 - eps, jnp.maximum(chi, eps), chi),
    )
    return X_nl_all_levs.at[..., hm_metadata.iiPDF_chi].set(chi)


# -----------------------------------------------------------------------------
def cdfnorminv(
    pdf_dim, nzt, ngrdcol, num_samples, X_u_all_levs,  # In
):
    """This function computes the inverse of the cumulative normal distribution function.
    The return value is the lower tail quantile for the standard normal distribution.
    This is equivalent to SQRT(2) * ERFINV(2*P-1), but is designed for computational
    efficiency on GPUs, however it also has a signficant performance boost when run
    on CPUs compared to the previously used ltqnorm. The GPU based performance mainly
    comes from the reduction of the chance for warp divergence.
    THIS FUNCTION ONLY HAS SINGLE PRECISION ACCURACY, BUT ACCEPTS DOUBLE PRECISION ARGUMENTS

    Arguments:
        pdf_dim: `d' Number of variates (normally 3 + microphysics specific variables)
        nzt: Number of vertical grid levels
        ngrdcol: Number of grid columns
        num_samples: Number of subcolumn samples
        X_u_all_levs: Uniform variates, arranged as (column, sample, level, variate)
    """
    # Coefficients of Mike Giles's inverse-erf approximation, in source
    # Horner order: a for the central region and b for the tails. The Fortran
    # literals have default-real precision before assignment to core_rknd;
    # preserve that rounding before promoting them to the sample dtype.
    a = jnp.array(
        [
            2.81022636e-8,
            3.43273939e-7,
            -3.5233877e-6,
            -4.39150654e-6,
            0.00021858087,
            -0.00125372503,
            -0.00417768164,
            0.246640727,
            1.50140941,
        ],
        dtype=jnp.float32,
    ).astype(X_u_all_levs.dtype)
    b = jnp.array(
        [
            -0.000200214257,
            0.000100950558,
            0.00134934322,
            -0.00367342844,
            0.00573950773,
            -0.0076224613,
            0.00943887047,
            1.00167406,
            2.83297682,
        ],
        dtype=jnp.float32,
    ).astype(X_u_all_levs.dtype)
    x = 2.0 * X_u_all_levs[..., :pdf_dim] - 1.0
    w = -jnp.log((1.0 - x) * (1.0 + x))

    # Central-region shift, or tail-region square-root shift.
    central = w < 5.0
    w = jnp.where(central, w - 2.5, jnp.sqrt(jnp.where(central, 1.0, w)) - 3.0)
    std_normal = jnp.sqrt(2.0) * x * jnp.where(central, jnp.polyval(a, w), jnp.polyval(b, w))

    # Preserve source std_normal(variable, column, level, sample) layout.
    return jnp.transpose(std_normal, (3, 0, 2, 1))


# -----------------------------------------------------------------------------
def ltqnorm(p_core_rknd):
    """This function is ported to Fortran from the same function written in Matlab,
    see the following description of this function.  Hongli Jiang, 2/17/2004
    Converted to double precision by Vince Larson 2/22/2004;
    this improves results for input values of p near 1.
    LTQNORM Lower tail quantile for standard normal distribution.
    Z = LTQNORM(P) returns the lower tail quantile for the standard normal
    distribution function.  I.e., it returns the Z satisfying Pr{X < Z} = P,
    where X has a standard normal distribution.
    LTQNORM(P) is the same as SQRT(2) * ERFINV(2*P-1), but the former returns a
    more accurate value when P is close to zero.
    The algorithm uses a minimax approximation by rational functions and the
    result has a relative error less than 1.15e-9. A last refinement by
    Halley's rational method is applied to achieve full machine precision.
    Author:      Peter J. Acklam
    Time-stamp:  2003-04-23 08:26:51 +0200
    E-mail:      pjacklam@online.no
    URL:         http://home.online.no/~pjacklam
    """
    # Coefficients in the rational approximations (Acklam's a, b, c, d arrays).
    a = jnp.array(
        [
            -39.69683028665376,
            220.9460984245205,
            -275.9285104469687,
            138.3577518672690,
            -30.66479806614716,
            2.506628277459239,
        ]
    )
    b = jnp.array(
        [
            -54.47609879822406,
            161.5858368580409,
            -155.6989798598866,
            66.80131188771972,
            -13.28068155288572,
            1.0,
        ]
    )
    c = jnp.array(
        [
            -0.007784894002430293,
            -0.3223964580411365,
            -2.400758277161838,
            -2.549732539343734,
            4.374664141464968,
            2.938163982698783,
        ]
    )
    d = jnp.array(
        [
            0.007784695709041462,
            0.3224671290700398,
            2.445134137142996,
            3.754408661907416,
            1.0,
        ]
    )
    p = jnp.asarray(p_core_rknd, dtype=jnp.float64)

    # JAX adaptation: mask the input before logarithms so invalid/end-point
    # values do not contaminate the inactive approximation branches.
    valid = (p > 0.0) & (p < 1.0)
    safe_p = jnp.where(valid, p, 0.5)

    # Rational approximations for the lower and upper regions share |q|.
    q = jnp.sqrt(-2.0 * jnp.log(jnp.minimum(safe_p, 1.0 - safe_p)))
    tail = jnp.polyval(c, q) / jnp.polyval(d, q)

    # Rational approximation for the central region; break points are 0.02425
    # and 1 - 0.02425, as in the source.
    q = safe_p - 0.5
    r = q * q
    z = jnp.where(
        safe_p < 0.02425,
        tail,
        jnp.where(safe_p > 1.0 - 0.02425, -tail, jnp.polyval(a, r) * q / jnp.polyval(b, r)),
    )

    # One iteration of Halley's rational method improves the relative error of
    # the approximation to full machine precision (Eric Raut, 23 August 2014).
    e = 0.5 * erfc(-z / jnp.sqrt(2.0)) - safe_p
    u = e * jnp.sqrt(2.0 * jnp.pi) * jnp.exp(z * z / 2.0)
    z = z - u / (1.0 + z * u / 2.0)

    # Source end-point/error results, expressed without dividing by zero.
    return jnp.where(valid, z, jnp.where(p == 0.0, -jnp.inf, jnp.where(p == 1.0, jnp.inf, jnp.nan)))


# -----------------------------------------------------------------------------
def multiply_Cholesky(
    nzt, ngrdcol, num_samples, pdf_dim, std_normal,  # In
    Sigma_Cholesky1, Sigma_Cholesky2,                # In
    mu1, mu2, X_mixt_comp_all_levs,                  # In
):
    """Computes X_nl_all_levs from the Cholesky factorization of Sigma,
    std_normal, and mu.
    X_nl_all_levs = Sigma_Cholesky * std_normal + mu.

    Arguments:
        nzt: Number of vertical grid levels
        ngrdcol: Number of grid columns
        num_samples: Number of samples
        pdf_dim: Number of variates (normally=5)
        std_normal: vector of d-variate standard normal distribution [-]
        Sigma_Cholesky1: Cholesky factorization of the Sigma matrix, 1st component [units
            vary]
        Sigma_Cholesky2: Cholesky factorization of the Sigma matrix, 2nd component [units
            vary]
        mu1: d-dimensional column vector of means of Gaussian, 1st component [units vary]
        mu2: d-dimensional column vector of means of Gaussian, 2nd component [units vary]
        X_mixt_comp_all_levs: Whether we're in the first or second mixture component
    """
    # Select each sample's mixture component and initialize with its mean.
    first = X_mixt_comp_all_levs == 1
    X_nl_all_levs = jnp.where(first[..., None], mu1[:, None], mu2[:, None])

    # Compute Sigma_Cholesky * std_normal. Sequential accumulation over j
    # preserves the source triangular sum; entries j > p contribute zero.
    def step(j, values):
        coeff = jnp.where(
            first[..., None], Sigma_Cholesky1[j, :, None], Sigma_Cholesky2[j, :, None]
        )
        coeff = jnp.where(jnp.arange(pdf_dim) >= j, coeff, 0.0)
        return values + coeff * jnp.swapaxes(std_normal[j], 1, 2)[..., None]

    return jax.lax.fori_loop(0, pdf_dim, step, X_nl_all_levs)


# -----------------------------------------------------------------------------
def chi_eta_2_rtthl(
    nzt, ngrdcol, num_samples,  # In
    rt_1, thl_1,                # In
    rt_2, thl_2,                # In
    crt_1, cthl_1,              # In
    crt_2, cthl_2,              # In
    mu_chi_1, mu_chi_2,         # In
    chi, eta,                   # In
    X_mixt_comp_all_levs,       # In
):
    """Converts from chi(s), eta(t) variables to rt, thl.  Also sets a limit on the value
    of cthl_1 and cthl_2 to prevent extreme values of temperature.

    Arguments:
        nzt: Vertical grid levels
        ngrdcol: Columns
        num_samples: Number of subcolumn samples
        rt_1: n dimensional column vector of rt [kg/kg]
        thl_1: n dimensional column vector of thetal [K]
        rt_2: n dimensional column vector of rt [kg/kg]
        thl_2: n dimensional column vector of thetal [K]
        crt_1: Constants from plumes 1 & 2 of rt
        cthl_1: Constants from plumes 1 & 2 of thetal
        crt_2: Constants from plumes 1 & 2 of rt
        cthl_2: Constants from plumes 1 & 2 of thetal
        mu_chi_1: Mean for chi_1 and chi_2 [kg/kg]
        mu_chi_2: Mean for chi_1 and chi_2 [kg/kg]
        chi: [kg/kg]
        eta: [-]
        X_mixt_comp_all_levs: Whether we're in the first or second mixture component
    """
    thl_dev_lim = 5.0  # Source local limit on temperature deviations [K].

    # Select component 1 or 2 for each sample before applying its transform.
    first = X_mixt_comp_all_levs == 1
    rt = jnp.where(first, rt_1[:, None], rt_2[:, None])
    thl = jnp.where(first, thl_1[:, None], thl_2[:, None])
    crt = jnp.where(first, crt_1[:, None], crt_2[:, None])
    cthl = jnp.where(first, cthl_1[:, None], cthl_2[:, None])
    mu_chi = jnp.where(first, mu_chi_1[:, None], mu_chi_2[:, None])
    lh_rt = rt + (0.5 / crt) * (chi - mu_chi) + (0.5 / crt) * eta

    # Limit the quantity by which temperature can vary [K].
    lh_dev_thl_lim = (-0.5 / cthl) * (chi - mu_chi) + (0.5 / cthl) * eta
    lh_dev_thl_lim = jnp.maximum(jnp.minimum(lh_dev_thl_lim, thl_dev_lim), -thl_dev_lim)
    return lh_rt, thl + lh_dev_thl_lim
