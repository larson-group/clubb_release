"""Check that microphysics uses each column's parameters and vertical metrics.

The PDF preparation check isolates its parameter boundary with a small setup
double. The hole-filling check exercises the actual driver, comparing batched
columns with independent calls and checking density-weighted mass conservation.
Both checks cover eager/JIT execution and derivatives across these boundaries.
References: pdf_hydromet_microphys_wrapper.F90 and fill_holes.F90.
"""

from types import SimpleNamespace

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.auto_llm_generated_tests.pytests.microphysics_test_inputs import initialize
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_driver_api
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.Microphys import pdf_hydromet_microphys_wrapper as wrapper

jax.config.update("jax_enable_x64", True)


@pytest.mark.parametrize("compiled", [False, True])
def test_pdf_preparation_preserves_per_column_parameters(monkeypatch, compiled):
    initialize(microphys_scheme="khairoutdinov_kogan")
    gr = SimpleNamespace(ngrdcol=2, nzt=3, nzm=4)
    stats = JaxStats.empty(l_sample=False, names=(), ncol=2, max_nlev=4)
    flags = SimpleNamespace(
        iiPDF_type=1, l_use_precip_frac=True, l_diagnose_correlations=False,
        l_calc_w_corr=False, l_const_Nc_in_cloud=True,
        l_fix_w_chi_eta_correlations=True,
    )

    def setup(*args):
        params = args[16]
        return (
            args[-4],
            jnp.repeat(params[..., None], 2, axis=2),
            *(jnp.zeros((2, 3, 4)),) * 4,
            *(jnp.zeros((2, 3, 4, 4)),) * 4,
            args[-3], args[-2], args[-1],
        )

    def mixed(*args):
        return *(jnp.zeros((2, 3, 2)),) * 3, args[-1]

    monkeypatch.setattr(wrapper, "setup_pdf_parameters_api", setup)
    monkeypatch.setattr(wrapper, "hydrometeor_mixed_moments", mixed)

    def run(params):
        return wrapper.pdf_hydromet_microphys_prep.__wrapped__(
            gr, 2, 4, 2, 0, 0.0,
            None, None, None, None, None, None, jnp.ones((2, 3, 2)), None,
            None, None, None, None, params, flags, None, False,
            stats, ErrInfo.initialized(2), None, None,
        )[2][..., 0]

    run = jax.jit(run) if compiled else run
    params = jnp.arange(1.0, 9.0).reshape(2, 4)
    result = run(params)
    assert result.shape == params.shape
    np.testing.assert_array_equal(result, params)
    tangent = jnp.array([[0.0] * 4, [0.1, 0.2, 0.3, 0.4]])
    np.testing.assert_array_equal(jax.jvp(run, (params,), (tangent,))[1], tangent)
    weights = jnp.arange(2.0, 10.0).reshape(2, 4)
    np.testing.assert_array_equal(jax.grad(lambda p: jnp.sum(run(p) * weights))(params), weights)


@pytest.mark.parametrize("compiled", [False, True])
@pytest.mark.parametrize("fill_type", [1, 2])
def test_hydrometeor_fill_uses_each_columns_metric(compiled, fill_type):
    _, _, metadata, *_ = initialize(microphys_scheme="khairoutdinov_kogan")
    nzt = 20
    dz = jnp.stack((jnp.linspace(20.0, 80.0, nzt), jnp.linspace(110.0, 35.0, nzt)))
    rho = jnp.stack((jnp.linspace(1.2, 0.8, nzt), jnp.linspace(1.1, 0.7, nzt)))
    mass = jnp.full((2, nzt), 2e-4).at[:, 9].set(-1e-4)
    hydromet = jnp.stack((mass, jnp.full_like(mass, 1e6)), axis=2)

    def run(hm, metrics, density):
        ncol = hm.shape[0]
        zero = jnp.zeros((ncol, nzt))
        stats = JaxStats.empty(l_sample=False, names=(), ncol=ncol, max_nlev=nzt)
        return fill_holes_driver_api(
            SimpleNamespace(dzt=metrics), ncol, nzt, 10.0, 2,
            metadata, False, density, jnp.ones_like(zero), fill_type,
            stats, zero, zero, hm,
        )[-1]

    run = jax.jit(run) if compiled else run
    result = run(hydromet, dz, rho)
    independent = jnp.concatenate([
        run(hydromet[i:i+1], dz[i:i+1], rho[i:i+1]) for i in range(2)
    ])
    assert result.shape == independent.shape == hydromet.shape
    np.testing.assert_allclose(result, independent, rtol=1e-13, atol=1e-15)
    weights = rho * dz
    np.testing.assert_allclose(
        jnp.sum(result[..., 0] * weights, axis=1),
        jnp.sum(mass * weights, axis=1), rtol=1e-13,
    )
    assert np.min(result[..., 0]) >= 0.0
    gradient = jax.grad(lambda hm: jnp.sum(run(hm, dz, rho)[..., 0] * weights))(hydromet)
    np.testing.assert_allclose(gradient[..., 0], weights, rtol=1e-13, atol=1e-13)
    np.testing.assert_array_equal(gradient[..., 1], 0.0)
