"""Native-random SILHS invariants and deterministic Fortran formula checks."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from scipy.special import ndtri
from clubb_jax.src.SILHS.latin_hypercube_arrays import LatinHypercubeArrays
from clubb_jax.src.SILHS.generate_uniform_sample_module import (
    generate_uniform_lh_sample,
)
from clubb_jax.src.SILHS.silhs_importance_sample_module import (
    define_importance_categories,
    compute_category_real_probs,
    compute_category_sample_weights,
    eight_cluster_allocation,
    four_cluster_no_precip,
    two_cluster_cp_nocp,
    limit_category_weights,
    importance_sampling_driver,
)
from clubb_jax.src.SILHS.transform_to_pdf_module import (
    cdfnorminv,
    ltqnorm,
    multiply_Cholesky,
)
from clubb_jax.src.SILHS.latin_hypercube_driver_module import compute_arb_overlap
from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision

configure_jax_precision()


def test_stratification_reproducibility_and_sequence_reuse():
    n, seq, d = 16, 3, 5
    initial = LatinHypercubeArrays(jnp.zeros((n * seq, d), jnp.int32), jnp.array(0))
    run = jax.jit(lambda it, key, state: generate_uniform_lh_sample(
        it, n, seq, d,  # In
        False,         # In
        key,            # In
        state,          # InOut
    ))
    x, state = run(1, jax.random.PRNGKey(41), initial)
    np.testing.assert_array_equal(jnp.floor(x * n * seq), state.one_height_time_matrix[:n])
    np.testing.assert_array_equal(
        jnp.sort(state.one_height_time_matrix, axis=0),
        jnp.broadcast_to(jnp.arange(n * seq)[:, None], (n * seq, d)),
    )
    same, _ = run(1, jax.random.PRNGKey(41), initial)
    np.testing.assert_array_equal(x, same)
    y, reused = run(2, jax.random.PRNGKey(42), state)
    np.testing.assert_array_equal(reused.one_height_time_matrix, state.one_height_time_matrix)
    assert not np.array_equal(x, y)
    _, new = run(4, jax.random.PRNGKey(43), reused)
    assert not np.array_equal(new.one_height_time_matrix, state.one_height_time_matrix)
    one = LatinHypercubeArrays(jnp.zeros((n, d), jnp.int32), jnp.array(0))
    x, _ = generate_uniform_lh_sample(
        1, n, 1, d,             # In
        False,                 # In
        jax.random.PRNGKey(1),  # In
        one,                    # InOut
    )
    np.testing.assert_array_equal(
        jnp.sort(jnp.floor(x * n), axis=0),
        jnp.broadcast_to(jnp.arange(n)[:, None], (n, d)),
    )
    assert np.all((np.asarray(x) > 0) & (np.asarray(x) < 1))


@pytest.mark.parametrize("sequence_length", [1, 3])
def test_deterministic_strata_cycle_and_sequence_reuse(sequence_length):
    n, d = 8, 5
    nt_repeat = n * sequence_length
    initial = LatinHypercubeArrays(jnp.zeros((nt_repeat, d), jnp.int32), jnp.array(0))
    run = jax.jit(lambda it, key, state: generate_uniform_lh_sample(
        it, n, sequence_length, d,  # In
        True,                      # In
        key,                        # In
        state,                      # InOut
    ))
    expected_matrix = np.broadcast_to(np.arange(nt_repeat)[:, None], (nt_repeat, d))
    expected_offsets = np.tile(np.array([
        [0.125, 0.625, 0.375, 0.875, 0.125],
        [0.625, 0.375, 0.875, 0.125, 0.625],
        [0.375, 0.875, 0.125, 0.625, 0.375],
        [0.875, 0.125, 0.625, 0.375, 0.875],
    ]), (2, 1))
    state = initial
    for iteration in range(1, sequence_length + 2):
        x, state = run(iteration, jax.random.PRNGKey(iteration), state)
        expected = (1.0 / nt_repeat) * (
            expected_matrix[:n] + np.roll(expected_offsets, -(iteration - 1), axis=0)
        )
        np.testing.assert_array_equal(x, expected)
        np.testing.assert_array_equal(np.floor(np.asarray(x) * nt_repeat), expected_matrix[:n])
        np.testing.assert_array_equal(state.one_height_time_matrix, expected_matrix)
        assert int(state.prior_iter) == (iteration if sequence_length > 1 else 1)
    np.testing.assert_array_equal(initial.one_height_time_matrix, 0)


def test_deterministic_overlap_pool_is_independent_of_seed():
    from clubb_jax.src.SILHS.latin_hypercube_driver_module import generate_random_pool

    run = jax.jit(lambda seed: generate_random_pool(
        4, 2, 5, 8, 2,  # In
        seed, None,       # In
        True,            # In
    ))
    pool = run(1)
    assert pool.shape == (2, 8, 4, 7)
    np.testing.assert_array_equal(run(99), pool)
    cycle = np.array([0.125, 0.625, 0.375, 0.875])
    np.testing.assert_array_equal(np.unique(pool), np.sort(cycle))
    np.testing.assert_array_equal(pool[0, :, 0, 0], np.tile(cycle, 2))
    np.testing.assert_array_equal(pool[0, 0, :, 0], cycle)
    np.testing.assert_array_equal(pool[0, 0, 0, :], np.tile(cycle, 2)[:7])
    np.testing.assert_array_equal(pool[1, :, :, :], np.roll(pool[0, :, :, :], -1, axis=0))


@pytest.mark.parametrize(
    "allocation",
    [eight_cluster_allocation, four_cluster_no_precip, two_cluster_cp_nocp],
)
@pytest.mark.parametrize("variance", [False, True])
@pytest.mark.parametrize(
    "fractions",
    [(0.2, 0.6, 0.3, 0.4, 0.7), (0.0, 0.0, 1.0, 0.0, 0.0), (1.0, 1.0, 0.0, 1.0, 1.0)],
)
def test_category_probabilities_weights_and_empty_clusters(allocation, variance, fractions):
    cat = define_importance_categories()
    real = compute_category_real_probs(cat, *fractions)
    prescribed = allocation(cat, real, variance)
    prescribed, error = limit_category_weights(real, prescribed)
    assert not bool(error)
    weights = compute_category_sample_weights(real, prescribed)
    np.testing.assert_allclose(jnp.sum(real), 1.0, atol=1.0e-15)
    np.testing.assert_allclose(jnp.sum(prescribed), 1.0, atol=1.0e-15)
    np.testing.assert_allclose(jnp.sum(weights * prescribed), 1.0, atol=1.0e-15)
    assert np.isfinite(weights).all()
    assert np.all(np.asarray(prescribed) >= -1.0e-15)
    assert np.all(np.asarray(weights)[np.asarray(real) >= 1.0e-8] <= 2.0 + 1.0e-14)


@pytest.mark.parametrize("strategy", [1, 2, 3])
@pytest.mark.parametrize("clustered", [False, True])
def test_importance_samples_membership_and_weighted_estimate(strategy, clustered, monkeypatch):
    from clubb_jax.src.CLUBB_core import error_code

    monkeypatch.setattr(error_code, "_debug_level", 2)
    # Unnormalized importance weights should recover the cloud/component/precip
    # probabilities. Thousands of stratified points avoid a flaky small-N oracle.
    n = 8192
    key = jax.random.PRNGKey(11)
    base = jax.random.uniform(key, (n, 3), dtype=jnp.float64)
    run = jax.jit(
        lambda x: importance_sampling_driver(
            n,                          # In
            0.2, 0.6,                   # In
            0.3,                        # In
            0.4, 0.7,                   # In
            strategy, clustered,        # In
            True, False, False,         # In
            x[:, 0], x[:, 1], x[:, 2],  # InOut
            jax.random.PRNGKey(12),     # In
        )
    )
    chi, comp, prec, w, error = run(base)
    assert not bool(error)
    cf = jnp.where(comp < 0.3, 0.2, 0.6)
    pf = jnp.where(comp < 0.3, 0.4, 0.7)
    np.testing.assert_allclose(jnp.mean(w * (chi >= 1.0 - cf)), 0.3 * 0.2 + 0.7 * 0.6, atol=0.015)
    np.testing.assert_allclose(jnp.mean(w * (comp < 0.3)), 0.3, atol=0.015)
    np.testing.assert_allclose(jnp.mean(w * (prec < pf)), 0.3 * 0.4 + 0.7 * 0.7, atol=0.015)


def test_inverse_normal_and_cholesky_reference():
    u = jnp.array([3.0e-8, 0.001, 0.1, 0.5, 0.9, 0.999, 1.0 - 3.0e-8]).reshape(1, 1, 1, 7)
    np.testing.assert_allclose(
        cdfnorminv(
            7, 1, 1, 1, u,  # In
        ).reshape(-1),
        ndtri(np.asarray(u).reshape(-1)),
        atol=3.0e-7,
    )
    np.testing.assert_allclose(ltqnorm(u), ndtri(u), atol=2.0e-10)
    normals = jnp.array([[[[1.0, 2.0]]], [[[3.0, 4.0]]]])
    sigma = jnp.array([[[[2.0, 0.5]]], [[[0.0, 3.0]]]])
    mu = jnp.array([[[10.0, 20.0]]])
    comp = jnp.array([[[1], [2]]])
    result = multiply_Cholesky(
        1, 1, 2, 2, normals,  # In
        sigma, 2.0 * sigma,   # In
        mu, -mu, comp,        # In
    )
    np.testing.assert_allclose(result, [[[[12.0, 29.5]], [[-2.0, 6.0]]]])


def test_vertical_overlap_matches_scalar_recurrence_and_maximum_overlap():
    rng = np.random.default_rng(2)
    ncol, ns, nz, d = 2, 3, 6, 4
    start = np.array([1, 4])
    corr = rng.uniform(0.1, 1.0, (ncol, nz))
    pool = rng.uniform(size=(ncol, ns, nz, d))
    x = np.zeros_like(pool)
    x[np.arange(ncol), :, start, :] = rng.uniform(0.1, 0.9, (ncol, ns, d))
    expected = x.copy()
    for i in range(ncol):
        for direction in (1, -1):
            point = x[i, :, start[i], :].copy()
            for k in range(start[i] + direction, nz if direction == 1 else -1, direction):
                width = 1.0 - corr[i, k]
                point = point - width + 2.0 * width * pool[i, :, k, :]
                point = np.where(
                    point > 1.0 - 3.0e-8,
                    2.0 - point - 6.0e-8,
                    np.where(point < 3.0e-8, -point + 6.0e-8, point),
                )
                expected[i, :, k, :] = point
    run = jax.jit(
        lambda c, p, v: compute_arb_overlap(
            nz, ncol, ns, d - 2, 2,  # In
            jnp.array(start), c, p,  # In
            v,                       # InOut
        )
    )
    np.testing.assert_allclose(run(corr, pool, x), expected, atol=3.0e-16)
    maximum = run(jnp.ones_like(corr), pool, x)
    np.testing.assert_allclose(
        maximum,
        np.broadcast_to(x[np.arange(ncol), :, start, :][:, :, None, :], pool.shape),
    )


@pytest.mark.parametrize("instantaneous", [False, True])
def test_microphysical_variance_sources_include_optional_dt_terms(instantaneous):
    from types import SimpleNamespace
    from clubb_jax.src.SILHS.lh_microphys_var_covar_module import (
        lh_microphys_var_covar_driver_api,
    )

    z = jnp.zeros((1, 2))
    pdf = SimpleNamespace(
        mixt_frac=0.5 * jnp.ones_like(z), rt_1=z, rt_2=z, thl_1=z, thl_2=z, w_1=z, w_2=z
    )
    rt = jnp.array([[[-1.0, -1.0], [1.0, 1.0]]])
    weights = jnp.ones_like(rt)
    run = jax.jit(
        lambda x: lh_microphys_var_covar_driver_api(
            2, 2, 1, 0.5, weights,      # In
            pdf, x, x, 4.0 * x,         # In
            0.0 * x, 2.0 * x, 3.0 * x,  # In
            instantaneous,              # In
        )
    )
    output = run(rt)
    expected = (4.0, 6.0, 8.0, 12.0, 5.0) if instantaneous else (6.0, 10.5, 8.0, 12.0, 8.0)
    for value, ref in zip(output, expected):
        np.testing.assert_array_equal(value, ref)


def test_category_rms_uses_weighted_second_moment_and_real_probability():
    from types import SimpleNamespace
    from clubb_jax.src.Microphys.silhs_category_variance_module import (
        silhs_sample_category_variance,
    )

    metadata = SimpleNamespace(iiPDF_chi=0, iiPDF_rr=1)
    half = jnp.full((1, 1), 0.5)
    pdf = SimpleNamespace(cloud_frac_1=half, cloud_frac_2=half, mixt_frac=half)
    precip = SimpleNamespace(precip_frac_1=half, precip_frac_2=half)
    x = jnp.array([[[[1.0, 1.0]], [[-1.0, 0.0]]]])
    comp = jnp.array([[[1], [2]]])
    values = jnp.array([[[2.0], [3.0]]])
    run = jax.jit(
        lambda samples: silhs_sample_category_variance(
            1, 1, 2, 2, x,                     # In
            comp, samples,                     # In
            jnp.ones((1, 2, 1)), pdf, precip,  # In
            metadata,                          # In
        )
    )
    np.testing.assert_array_equal(run(values), [[[4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0]]])


@pytest.mark.parametrize("all_levels", [False, True])
@pytest.mark.parametrize("straight_mc", [False, True])
def test_uniform_driver_options_and_shared_permutation_state(all_levels, straight_mc):
    from types import SimpleNamespace
    from clubb_jax.src.SILHS.latin_hypercube_driver_module import (
        generate_all_uniform_samples,
    )

    ncol, nz, n, dim = 2, 4, 8, 4
    one = jnp.ones((ncol, nz))
    precip = SimpleNamespace(precip_frac_1=0.4 * one, precip_frac_2=0.7 * one)
    pool = jax.random.uniform(jax.random.PRNGKey(1), (ncol, n, nz, dim + 2), dtype=jnp.float64)
    initial = LatinHypercubeArrays(jnp.zeros((n, dim + 2), jnp.int32), jnp.array(0, jnp.int32))
    run = jax.jit(
        lambda state: generate_all_uniform_samples(
            1, dim, 2, n, 1,                         # In
            nz, ncol, jnp.array([1, 2]), one, pool,  # In
            0,                                       # In
            0.2 * one,                               # In
            0.6 * one,                               # In
            0.3 * one, precip,                       # In
            3,                                       # In
            True,                                    # In
            straight_mc,                             # In
            True,                                    # In
            True,                                    # In
            False,                                   # In
            True,                                    # In
            False,                                   # In
            all_levels,                              # In
            jax.random.PRNGKey(2),                   # In
            state,                                   # InOut
        )
    )
    x, w, state, error = run(initial)
    assert x.shape == pool.shape and w.shape == pool.shape[:3]
    assert not bool(error)
    assert np.all((np.asarray(x) > 0) & (np.asarray(x) < 1))
    np.testing.assert_allclose(jnp.sum(w, axis=1), n, atol=1.0e-13)
    if straight_mc:
        np.testing.assert_array_equal(state.one_height_time_matrix, initial.one_height_time_matrix)
        if all_levels:
            np.testing.assert_array_equal(x, pool)
        else:
            # Fortran still applies maximum vertical overlap to straight MC.
            starts = pool[jnp.arange(ncol), :, jnp.array([1, 2]), :]
            np.testing.assert_array_equal(x, jnp.broadcast_to(starts[:, :, None, :], x.shape))
    elif not all_levels:
        np.testing.assert_array_equal(x, jnp.broadcast_to(x[:, :, :1, :], x.shape))


@pytest.mark.parametrize("straight_mc", [False, True])
def test_old_cloud_weighted_sample_count_check_precedes_dispatch(monkeypatch, straight_mc):
    import inspect
    from clubb_jax.src.SILHS import latin_hypercube_driver_module as driver

    monkeypatch.setattr(driver, "l_lh_old_cloud_weighted", True)
    # The source rejects odd sample counts before reading either branch's
    # arrays, including when straight MC bypasses the old weighted sampler.
    args = {name: None for name in inspect.signature(driver.generate_all_uniform_samples).parameters}
    args.update(num_samples=3, l_lh_straight_mc=straight_mc)
    with pytest.raises(ValueError, match="even sample count"):
        driver.generate_all_uniform_samples(**args)


def test_invalid_category_draw_is_reported_in_jit(monkeypatch, capsys):
    from clubb_jax.src.CLUBB_core import error_code
    from clubb_jax.src.SILHS.silhs_importance_sample_module import (
        pick_sample_categories,
    )

    monkeypatch.setattr(error_code, "_debug_level", 0)
    result = jax.jit(lambda r: pick_sample_categories(3, jnp.full(8, 0.125), r))(
        jnp.array([-0.1, 0.5, 1.1])
    )
    jax.block_until_ready(result)
    jax.effects_barrier()
    np.testing.assert_array_equal(result, [-1, 4, -1])
    assert "Invalid rand_vect number" in capsys.readouterr().out


@pytest.mark.parametrize("direction", [1, -1])
def test_random_starting_level_bounds_and_no_cloud_fallback(direction):
    from types import SimpleNamespace
    from clubb_jax.src.SILHS.latin_hypercube_driver_module import compute_k_lh_start

    ncol = 8192
    gr = SimpleNamespace(grid_dir_indx=direction)
    rc = jnp.broadcast_to(jnp.array([0.0, 0.08, 0.06, 0.03, 0.0]), (ncol, 5))
    cf = jnp.broadcast_to(jnp.array([0.0, 0.9, 0.6, 0.1, 0.0]), (ncol, 5))
    pdf = SimpleNamespace(mixt_frac=jnp.full_like(cf, 0.5), cloud_frac_1=cf, cloud_frac_2=cf)
    run = jax.jit(lambda r: compute_k_lh_start(
        gr, 5, ncol, r, pdf,  # In
        True,                 # In
        True,                 # In
        17,                   # In
    ))
    levels = run(rc)
    assert int(jnp.min(levels)) == 1
    assert int(jnp.max(levels)) == 3
    np.testing.assert_allclose(jnp.mean(levels), 2.0, atol=0.025)
    np.testing.assert_array_equal(run(jnp.zeros_like(rc)), 1 if direction > 0 else 3)
