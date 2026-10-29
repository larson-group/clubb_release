"""Check shared Dv routing, degenerate PDFs, and differentiation of KK covariances."""
import gc

import jax
import jax.numpy as jnp
import numpy as np
import pytest


from clubb_jax.src.Microphys.KK_microphys import KK_upscaled_covariances as kk

jax.config.update("jax_enable_x64", True)


@pytest.fixture(autouse=True)
def release_compilations():
    yield
    jax.clear_caches()
    gc.collect()


def _inputs(degenerate=False):
    """Distinct PDF components, both signs of chi, and all 16 variance masks."""
    index = np.arange(16).reshape(2, 8)
    full = lambda value: jnp.full(index.shape, value, dtype=jnp.float64)
    values = dict(w_mean=full(.1), rtm=full(.01), thlm=full(300.), exner=full(.95),
                  mixt_frac=full(.37), precip_frac_1=full(.6), precip_frac_2=full(.8),
                  KK_evap_coef=full(.2), KK_auto_coef=full(1350.), KK_accr_coef=full(67.),
                  KK_evap_tndcy=full(-1e-7), KK_auto_tndcy=full(2e-7), KK_accr_tndcy=full(3e-7),
                  crt1=full(.6), crt2=full(.7), cthl1=full(.4), cthl2=full(.5))
    for i in (1, 2):
        factor = 1.0 + .2 * i
        for name, mean in (('w', .3), ('eta', .2), ('rr', 1e-4), ('Nr', 1e6),
                           ('Ncn', 1e8), ('rt', .011), ('thl', 300.1)):
            values[f'mu_{name}_{i}'] = full(mean * factor)
        for name, mean in (('rr', -9.), ('Nr', 14.), ('Ncn', 18.)):
            values[f'mu_{name}_{i}_n'] = full(mean + .1 * i)
            values[f'sigma_{name}_{i}_n'] = full(.3 * factor)
        for name, sigma, bit in (('w', .4, 1), ('eta', .3, 1), ('chi', .2, 2),
                                 ('rr', 4e-5, 4), ('Nr', 5e5, 8), ('Ncn', 5e7, 4)):
            small = kk._CHI_TOL if name == 'chi' and i == 2 else 0.0
            values[f'sigma_{name}_{i}'] = jnp.where(
                (index & bit) != 0, small, sigma * factor) if degenerate else full(sigma * factor)
        ratio = jnp.where(index % 2 == 0, -.8, .7) * (1.0 if i == 1 else -1.3)
        values[f'mu_chi_{i}'] = ratio * jnp.maximum(values[f'sigma_chi_{i}'], kk._CHI_TOL)
        for pair in ('w_chi', 'chi_eta'):
            values[f'corr_{pair}_{i}'] = full(.12 * factor)
        for pair in ('w_rr', 'w_Nr', 'w_Ncn', 'chi_rr', 'chi_Nr', 'chi_Ncn',
                     'eta_rr', 'eta_Nr', 'eta_Ncn', 'rr_Nr'):
            values[f'corr_{pair}_{i}_n'] = full(.08 * factor)
    return values


def _without_batch(monkeypatch, inputs):
    # Exercise the original per-integral Dv calls as an independent routing oracle.
    with monkeypatch.context() as patch:
        patch.setattr(kk, 'batch_covariance_dv', lambda *args: (None, None, None))
        return jax.device_get(kk.KK_upscaled_covar_driver(**inputs))


@pytest.mark.parametrize('degenerate', [False, True])
def test_driver_matches_unbatched_integrals(monkeypatch, degenerate):
    inputs = _inputs(degenerate)
    reference = _without_batch(monkeypatch, inputs)
    compiled = jax.jit(kk.KK_upscaled_covar_driver)
    actual = compiled(**inputs)
    for got, expected in zip(actual, reference, strict=True):
        assert np.isfinite(got).all() and np.isfinite(expected).all()
        np.testing.assert_allclose(got, expected, rtol=2e-12, atol=1e-20)
    # More than dry/zero tendencies must be exercised in the normal PDF case.
    if not degenerate:
        assert all(np.max(np.abs(value)) > 1e-15 for value in reference)


def test_nested_jit_gradient_matches_unbatched_finite_difference(monkeypatch):
    inputs = _inputs()
    # Follow chi into the batch: a detached/cached Dv table would give a wrong derivative.
    scales = jnp.asarray([max(float(np.max(np.abs(x))), 1e-20)
                          for x in _without_batch(monkeypatch, inputs)])

    def loss(scale, inputs, scales):
        updated = dict(inputs, mu_chi_1=inputs['mu_chi_1'] * scale)
        outputs = kk.KK_upscaled_covar_driver(**updated)
        return sum(jnp.sum(value / scales[i]) for i, value in enumerate(outputs))

    derivative = jax.jit(jax.grad(loss))(1.0, inputs, scales)
    with monkeypatch.context() as patch:
        patch.setattr(kk, 'batch_covariance_dv', lambda *args: (None, None, None))
        epsilon = 1e-5
        finite_difference = (loss(1.0 + epsilon, inputs, scales)
                             - loss(1.0 - epsilon, inputs, scales)) / (2 * epsilon)
    assert np.isfinite(derivative) and np.isfinite(finite_difference)
    np.testing.assert_allclose(derivative, finite_difference, rtol=2e-6, atol=1e-7)
