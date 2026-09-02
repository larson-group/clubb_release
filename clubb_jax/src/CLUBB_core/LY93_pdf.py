"""JAX implementations of selected routines from ``LY93_pdf.F90``.

Description:
The multivariate, two-component PDF of Lewellen and Yoh (1993).

References:
Lewellen, W. S. and Yoh, S., 1993.  Binormal Model of Ensemble Partial
Cloudiness.  J. Atmos. Sci., 50, 9, 1228--1237.

Porting deviations:
The Fortran routines take an explicit ``nz`` and loop over 1-D vertical
profiles; the JAX routines broadcast over the provided array shapes.
Fortran mutates output arrays, while these Python functions return tuples.
``calc_mixt_frac_LY93`` uses a fixed-length ``lax.scan`` with frozen converged
points instead of an unbounded ``do`` loop so it remains JIT compatible.
``_scbrt`` is a local JAX helper for the real cube root used in the LY93 formula.
"""
from jax import lax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()


# Grad-safe sqrt(max(x,0)); the Fortran uses sqrt directly.
from clubb_jax.src.CLUBB_core.pdf_utilities import _safe_sqrt as _ssqrt


def _scbrt(x):
    """Real cube root of x>=0 with a finite gradient at x=0 (cbrt'(0)=inf otherwise)."""
    xp = jnp.where(x > 0.0, x, 1.0)
    return jnp.where(x > 0.0, jnp.cbrt(xp), 0.0)


def calc_params_LY93(xm, xp2, Skx, mixt_frac):
    """Calculate the LY93 PDF parameters for one variable.

    Description:
    Calculates the PDF component means and PDF component variances for
    variable x according to Lewellen and Yoh.

    References:
    Eq. (14), Eq. (15), Eq. (16), Eq. (17), and Eq. (18) of
    Lewellen, W. S. and Yoh, S., 1993.  Binormal Model of Ensemble Partial
    Cloudiness.  J. Atmos. Sci., 50, 9, 1228--1237.
    """
    xm = jnp.asarray(xm, dtype=jnp.float64)
    xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    Skx = jnp.asarray(Skx, dtype=jnp.float64)
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64)
    omf = 1.0 - mf

    # Find the sign of Skx
    sgn = jnp.where(Skx >= 0.0, 1.0, -1.0)

    # Calculate B_x, the LY function for the spread of the PDF component means.
    B_x = sgn * _ssqrt(xp2) * _scbrt(jnp.abs(Skx) / omf)

    # Calculate the mean of x in the 1st PDF component.
    mu_x_1 = xm - B_x * omf

    # Calculate the mean of x in the 2nd PDF component.
    mu_x_2 = xm + B_x * mf

    # Calculate the variance of x in the 1st PDF component.
    sigma_x_1_sqd = xp2 - B_x ** 2 * omf * (1.0 + mf + mf ** 2) / (3.0 * mf)

    # Calculate the variance of x in the 2nd PDF component.
    sigma_x_2_sqd = xp2 + B_x ** 2 * omf ** 2 / 3.0
    return mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd


def calc_mixt_frac_LY93(Sk_max, itermax=60):
    """Calculate the LY93 mixture fraction.

    Description:
    Calculates mixture fraction iteratively according to Lewellen and Yoh.

    References:
    Eq. (21) of Lewellen, W. S. and Yoh, S., 1993.  Binormal Model of Ensemble
    Partial Cloudiness.  J. Atmos. Sci., 50, 9, 1228--1237.
    """
    Sk_max = jnp.asarray(Sk_max, dtype=jnp.float64)

    # Tolerance for mixture fraction in solver [-]
    tol = 1.0e-4
    use_iter = Sk_max > 0.84

    def step(carry, _):
        mf, low, high, done = carry
        expr = mf ** 6 - Sk_max ** 2 * (1.0 - mf)
        hit = jnp.abs(expr) < tol
        new_high = jnp.where(expr > 0.0, mf, high)
        new_low = jnp.where(expr < 0.0, mf, low)
        new_done = done | hit

        # Mixture fraction has been solved for within the specificed
        # tolerance.
        new_mf = jnp.where(new_done, mf, 0.5 * (new_low + new_high))
        new_low = jnp.where(done, low, new_low)
        new_high = jnp.where(done, high, new_high)
        return (new_mf, new_low, new_high, new_done), None

    half = 0.5 * jnp.ones_like(Sk_max)
    one_ = jnp.ones_like(Sk_max)
    init = (0.75 * jnp.ones_like(Sk_max), half, one_, jnp.zeros_like(Sk_max, dtype=bool))
    (mf_final, _, _, _), _ = lax.scan(step, init, None, length=itermax)
    return jnp.where(use_iter, mf_final, 0.75)


def LY93_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl):
    """Calculate the LY93 PDF parameters for w, rt, and theta-l.

    Description:
    Calculates the mixture fraction and the PDF component means and PDF
    component variances of w, rt, and theta-l following Lewellen and Yoh.

    References:
    Lewellen, W. S. and Yoh, S., 1993.  Binormal Model of Ensemble Partial
    Cloudiness.  J. Atmos. Sci., 50, 9, 1228--1237.
    """
    Skw = jnp.asarray(Skw, dtype=jnp.float64)
    Skrt = jnp.asarray(Skrt, dtype=jnp.float64)
    Skthl = jnp.asarray(Skthl, dtype=jnp.float64)

    # Find the maximum of the magnitudes of skewness.
    Sk_max = jnp.maximum(jnp.maximum(jnp.abs(Skw), jnp.abs(Skrt)), jnp.abs(Skthl))

    # Calculate mixture fraction.
    mixt_frac = calc_mixt_frac_LY93(Sk_max)

    # Calculate the PDF parameters for w.
    mu_w_1, mu_w_2, sig_w_1, sig_w_2 = calc_params_LY93(wm, wp2, Skw, mixt_frac)

    # Calculate the PDF parameters for rt.
    mu_rt_1, mu_rt_2, sig_rt_1, sig_rt_2 = calc_params_LY93(rtm, rtp2, Skrt, mixt_frac)

    # Calculate the PDF parameters for thl.
    mu_thl_1, mu_thl_2, sig_thl_1, sig_thl_2 = calc_params_LY93(thlm, thlp2, Skthl, mixt_frac)
    return (mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2,
            sig_w_1, sig_w_2, sig_rt_1, sig_rt_2, sig_thl_1, sig_thl_2, mixt_frac)
