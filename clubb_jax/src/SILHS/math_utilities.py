"""Weighted sample moments from math_utilities.F90."""

import jax
import jax.numpy as jnp


# -----------------------------------------------------------------------------
def compute_sample_mean(
    n_levels, n_samples, ngrdcol,  # In
    weight, x_sample,              # In
):
    """Find the mean of a set of sample points in every model column.

    Arguments:
        n_levels: Number of sample levels
        n_samples: Number of sample points
        ngrdcol: Number of model columns
        weight: Weights for individual points of the vector
        x_sample: Collection of sample points [units vary]
    """
    # Preserve sample accumulation order for each column; batch columns/levels
    # in JAX rather than using the source's contiguous innermost column loop.
    # weight and x_sample: (ngrdcol, n_samples, n_levels) [units vary].
    mean, _ = jax.lax.scan(
        lambda mean, xs: (mean + xs[0] * xs[1], None),
        jnp.zeros_like(x_sample[:, 0]),
        (jnp.swapaxes(weight, 0, 1), jnp.swapaxes(x_sample, 0, 1)),
    )
    return mean / n_samples


# -----------------------------------------------------------------------------
def compute_sample_variance(
    n_levels, n_samples, ngrdcol,  # In
    x_sample, weight, x_mean,      # In
):
    """Compute the variance of a set of sample points in every model column.

    Arguments:
        n_levels: Number of sample levels in the mean / variance
        n_samples: Number of sample points to compute the variance of
        ngrdcol: Number of model columns
        x_sample: Collection of sample points [units vary]
        weight: Coefficient to weight the nth sample point by [-]
        x_mean: Mean sample points [units vary]
    """
    # Weighted squared departures from the supplied mean, divided by n_samples.
    # This is the source population moment, without an n_samples - 1 correction.
    return compute_sample_mean(
        n_levels, n_samples, ngrdcol,               # In
        weight, (x_sample - x_mean[:, None]) ** 2,  # In
    )


# -----------------------------------------------------------------------------
def compute_sample_covariance(
    n_levels, n_samples, ngrdcol,  # In
    weight, x_sample, x_mean,      # In
    y_sample, y_mean,              # In
):
    """Compute the covariance of a set of sample points of 2 variables
    in every model column.

    Arguments:
        n_levels: Number of sample levels in the mean / variance
        n_samples: Number of sample points to compute the covariance of
        ngrdcol: Number of model columns
        weight: Coefficient to weight the nth sample point by [-]
        x_sample: Collection of sample points [units vary]
        x_mean: Mean sample points [units vary]
        y_sample: Collection of sample points [units vary]
        y_mean: Mean sample points [units vary]
    """
    # x_mean and y_mean have shape (ngrdcol, n_levels). The means may come from
    # the analytic PDF rather than from the finite sample.
    return compute_sample_mean(
        n_levels, n_samples, ngrdcol,                                         # In
        weight, (x_sample - x_mean[:, None]) * (y_sample - y_mean[:, None]),  # In
    )


# -----------------------------------------------------------------------------
def rand_integer_in_range(low, high, key):
    """Return a uniformly distributed integer in the inclusive range [low, high].

    JAX adaptation: use an explicit native key instead of the source MT95
    integer draw and MOD reduction; randint's exclusive upper bound is high + 1.

    Arguments:
        low: Lowest possible returned value
        high: Highest possible returned value
        key: Native JAX random key for these draws
    """
    return jax.random.randint(key, jnp.shape(low), low, high + 1)
