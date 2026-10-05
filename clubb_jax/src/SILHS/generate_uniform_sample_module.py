"""Uniform Latin-hypercube samples from generate_uniform_sample_module.F90.

User-requested adaptation: native JAX random operations replace MT95/CLUBB
random machinery. Keys and permutation storage are explicit inputs/returns.
"""

import jax
import jax.numpy as jnp
from clubb_jax.src.SILHS.latin_hypercube_arrays import LatinHypercubeArrays


# -----------------------------------------------------------------------------
def rand_uniform_real(key, shape=()):
    """Generate uniformly distributed real numbers in the open interval (0, 1).

    JAX adaptation: key and shape replace the implicit MT95 generator. The
    open interval prevents infinite values in inverse-CDF transformations.
    """
    # The source clips a rounded value of one when generator/model kinds differ.
    # Native draws exclude one; the positive lower bound also excludes zero.
    return jax.random.uniform(
        key, shape, dtype=jnp.float64, minval=jnp.finfo(jnp.float64).eps, maxval=1.0
    )


# -----------------------------------------------------------------------------
def generate_uniform_lh_sample(
    iter, num_samples, sequence_length, n_vars,  # In
    l_lh_deterministic_test,                     # In
    key,                                         # In
    sampling_state,                              # InOut
):
    """Generates a matrix X that contains a Latin Hypercube sample.
    The sample is uniformly distributed.

    iter is the model iteration; num_samples is n; sequence_length is n_t;
    n_vars is the number of uniform variates. X_u_one_lev has shape
    (num_samples, n_vars), with one n_vars-dimensional sample per row.

    References:
        Art B. Owen (2003), "Quasi-Monte Carlo Sampling," SIGGRAPH 2003.
        https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:lh_algorithm

    JAX adaptation: initialization allocates one_height_time_matrix; updated
    permutation/prior-iteration storage is returned with X_u_one_lev.

    Arguments:
        iter: Model iteration number
        num_samples: `n' Number of samples generated
        sequence_length: `n_t' Number of timesteps before the permutation repeats
        n_vars: Number of uniform variables to generate
        l_lh_deterministic_test: Ordered strata and repeating draws for repeatable testing
        key: Native JAX random key for these draws
        sampling_state: Case-owned permutation and prior iteration; updated state is returned
    """
    nt_repeat = num_samples * sequence_length

    # Latin hypercube sample generation: generate an nt_repeat x n_vars array
    # of random integers when this iteration starts a new sequence.
    i_rmd = (iter - 1) % sequence_length
    perm_key, draw_key = jax.random.split(key)
    # Repeatable testing only; numerically acceptable sampling is not guaranteed.
    # Use strata 0, ..., nt_repeat-1 for every variate. Refresh at the start
    # of each sequence and reuse the first num_samples rows between refreshes.
    one_height_time_matrix = jax.lax.cond(
        i_rmd == 0,
        lambda _: (
            jnp.broadcast_to(jnp.arange(nt_repeat, dtype=jnp.int32)[:, None], (nt_repeat, n_vars))
            if l_lh_deterministic_test else permute_height_time(nt_repeat, n_vars, perm_key)
        ),
        lambda _: sampling_state.one_height_time_matrix,
        operand=None,
    )

    # Choose values using the permuted vector and uniform random jitter.
    # Preserve Fortran's first-num_samples row selection between reshuffles.
    if l_lh_deterministic_test:
        # Repeatable testing only; numerically acceptable sampling is not guaranteed.
        # Cycle through offsets (0.125, 0.625, 0.375, 0.875) using the sum of
        # zero-based sample, variate and timestep indices modulo four.
        # Each point is (stored stratum + offset)/nt_repeat.
        test_uniform_draws = jnp.array((0.125, 0.625, 0.375, 0.875))
        j = jnp.arange(num_samples, dtype=jnp.int32)[:, None]
        k = jnp.arange(n_vars, dtype=jnp.int32)[None, :]
        i_draw = (j + k + (iter - 1)) % 4
        X_u_one_lev = (1.0 / nt_repeat) * (
            one_height_time_matrix[:num_samples] + test_uniform_draws[i_draw]
        )
    else:
        X_u_one_lev = choose_permuted_random(nt_repeat, one_height_time_matrix[:num_samples], draw_key)

    # Check that the iteration number increments correctly. In the source,
    # allocation initializes prior_iter; here its zero value marks the first call.
    # Match the source diagnostic's prior-iteration update rules. Duplicate
    # calls in a multi-column sequence warn but do not alter sampling results.
    prior_iter = sampling_state.prior_iter
    if sequence_length > 1:

        def warn(_):
            jax.debug.print(
                "The iteration number in latin_hypercube_driver is not incrementing properly: "
                "prior={p}, current={i}",
                p=prior_iter,
                i=iter,
            )
            return None

        jax.lax.cond(
            (prior_iter != 0) & (prior_iter != iter - 1),
            warn,
            lambda _: None,
            operand=None,
        )
        prior_iter = jnp.where((prior_iter == 0) | (prior_iter == iter - 1), iter, prior_iter)
    else:
        prior_iter = jnp.where(prior_iter == 0, iter, prior_iter)

    return X_u_one_lev, LatinHypercubeArrays(
        one_height_time_matrix, jnp.asarray(prior_iter, dtype=jnp.int32)
    )


# -----------------------------------------------------------------------------
def choose_permuted_random(nt_repeat, p_matrix_element, key):
    """Chooses a permuted random, using native JAX jitter in each stratum.

    Arguments:
        nt_repeat: Number of samples before the sequence repeats
        p_matrix_element: Permuted integer
        key: Native JAX random key for these draws
    """
    rand = rand_uniform_real(key, jnp.shape(p_matrix_element))
    return (1.0 / nt_repeat) * (p_matrix_element + rand)


# -----------------------------------------------------------------------------
def permute_height_time(nt_repeat, n_vars, key):
    """Generates a matrix one_height_time_matrix, which is a nt_repeat x n_vars
    matrix whose 1st dimension is random permutations of the integer sequence
    (0,...,nt_repeat-1).

    Arguments:
        nt_repeat: Total number of sample points before sequence repeats.
        n_vars: The number of variates in the uniform sample
        key: Native JAX random key for these draws
    """
    # Choose an integer Latin-hypercube permutation for each uniform variate.
    return jax.vmap(lambda k: rand_permute(nt_repeat, k))(jax.random.split(key, n_vars)).T


# -----------------------------------------------------------------------------
def rand_permute(n, key):
    """Generate a vector containing 0, ..., n - 1 in random order.

    References:
        Art Owen, "Quasi-Monte Carlo Sampling," Section 1.3, following
        Luc Devroye, "Non-Uniform Random Variate Generation" (1986).

    JAX adaptation: permutation uses the supplied key; it does not reseed an
    implicit generator or port the source MT95 shuffle.

    Arguments:
        n: Number of elements to permute
        key: Native JAX random key for these draws
    """
    return jax.random.permutation(key, jnp.arange(n, dtype=jnp.int32))
