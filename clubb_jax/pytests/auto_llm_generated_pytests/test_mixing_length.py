"""Check legacy parcel stopping, grid orientation, errors, and derivatives.

Analytic first-layer roots check energy exhaustion. Legacy and tau checks
compare forward/reverse derivatives and finite differences at smooth inputs,
including zero-energy finiteness and runtime parameter reuse.
"""
import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core import error_code
from clubb_jax.src.CLUBB_core.constants_clubb import ep1, grav
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.mixing_length import (
    calc_Lscale, compute_mixing_length,
)
from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params, PNAME_IDX
from clubb_jax.src.CLUBB_core.pdf_params import init_pdf_params


@pytest.fixture(autouse=True)
def restore_debug():
    """Run each test with debug -1 and restore the previous global level after it completes."""
    previous = error_code._debug_level
    error_code.set_debug_level(-1)
    yield
    error_code.set_debug_level(previous)


def fixture():
    """Return a two-column moist/stratified grid, thermodynamics, PDF departures, and
    dynamic test inputs.
    """
    gr = setup_grid(ngrdcol=2, deltaz=80., zm_init=0., zm_top=960.)
    thlm = 298. + .003 * gr.zt + jnp.array([[0.], [.3]])
    rt = jnp.broadcast_to(jnp.array([[.008], [.018]]), thlm.shape)
    exner = jnp.exp(-3e-5 * gr.zt)
    pdf = init_pdf_params(gr.nzt, 2).replace(
        thl_1=thlm + .7, thl_2=thlm - .4, rt_1=rt + .001, rt_2=rt - .0005)
    cp = jnp.asarray(init_clubb_params(2, ''))
    inputs = dict(thlm=thlm, rtm=rt, em=jnp.full((2, gr.nzm), .4), clubb_params=cp)
    return gr, exner, 1e5 * exner**3.5, thlm * (1. + ep1 * rt), pdf, inputs


def test_legacy_first_layer_root_and_debug_predicate():
    """Compare low-energy dry parcel lengths with analytic first-layer roots while debug
    checks are enabled.
    """
    gr = setup_grid(ngrdcol=2, deltaz=80., zm_init=0., zm_top=960.)
    thl = 298. + .003 * gr.zt
    mu = jnp.array([.001, .002])
    energy = .0001
    error_code.set_debug_level(2)
    actual = jax.jit(lambda e: compute_mixing_length(
        gr.nzm, gr.nzt, 2, gr, thl, thl, jnp.zeros_like(thl), e,
        jnp.full((2,), 1e5), jnp.full_like(thl, 1e5), jnp.ones_like(thl), thl,
        mu, 0., 3, False, ErrInfo.initialized(2)))(jnp.full((2, gr.nzm), energy))
    a = (1.-np.exp(-np.asarray(mu)*80.))/(np.asarray(mu)*80.)
    up, down = np.full(thl.shape, .1), np.full(thl.shape, .1)
    du = -grav/np.asarray(thl[:, 1:-1]) * .003*80.*a[:, None]
    dd = grav/np.asarray(thl[:, :-1]) * .003*80.*a[:, None]
    up[:, :-2] += np.sqrt(-2.*energy*80.*du)/(-du)
    down[:, 1:] += np.sqrt(2.*energy*80.*dd)/dd
    combined = np.sqrt(up*down)
    combined[:, -1] = combined[:, -2]
    assert not actual[0].is_fatal()
    for result, reference in zip(actual[1:], [combined, up, down]):
        assert result.shape == reference.shape
        np.testing.assert_allclose(result, reference, rtol=2e-10, atol=1e-10)


@pytest.mark.parametrize('method', [1, 2])
def test_methods_gradients_and_runtime_parameter_reuse(method):
    """Check each method's gradients and compilation reuse on two moist columns.

    Require finite zero-energy derivatives, JVP/reverse/finite-difference
    agreement at smooth inputs, and runtime parameter reuse.
    """
    gr, exner, pressure, thvm, pdf, x = fixture()
    stats = JaxStats.empty(l_sample=False, names=(), ncol=2, max_nlev=gr.nzm)
    def calculate(x):
        """Return full outputs for the selected method using the dynamic test input dictionary."""
        return calc_Lscale(gr.nzm, gr.nzt, 2, gr, False, jnp.full((2,), 1e5),
            jnp.full((2,), 1e5), pressure, exner, x['rtm'], x['thlm'], thvm,
            jnp.full((2, gr.nzm), .1), jnp.full((2, gr.nzm), 1e-6),
            jnp.zeros((2, gr.nzm)), pdf, x['em'], x['thlm'], 0.,
            jnp.full((2,), .02), jnp.full((2,), .01), jnp.full((2, gr.nzm), 1e-5),
            jnp.zeros_like(thvm), .1, 300., jnp.zeros((2,)), x['clubb_params'],
            3, False, method, False, True, False, False,
            jnp.full((2, gr.nzm), 1e-5), stats, ErrInfo.initialized(2))
    def loss(x):
        """Return a vertically weighted sum of combined lengths and the calculation error state."""
        result = calculate(x)
        return jnp.sum(result[12] * jnp.linspace(.5, 1.5, gr.nzt)), result[0]
    compiled = jax.jit(jax.value_and_grad(loss, has_aux=True))
    (_, error), gradient = compiled(x)
    assert not error.is_fatal()
    assert all(np.isfinite(g).all() for g in gradient.values())
    direction = dict(thlm=jnp.cos(gr.zt)*.01, rtm=jnp.ones_like(x['rtm'])*1e-6,
                     em=jnp.ones_like(x['em'])*.01, clubb_params=jnp.zeros_like(x['clubb_params']))
    direction['clubb_params'] = direction['clubb_params'].at[:, PNAME_IDX['Lscale_mu_entrain']].set(1e-5)
    scalar = lambda values: loss(values)[0]
    _, tangent = jax.jit(lambda x, d: jax.jvp(scalar, (x,), (d,)))(x, direction)
    reverse = sum(jnp.vdot(gradient[k], direction[k]) for k in x)
    np.testing.assert_allclose(tangent, reverse, rtol=2e-9, atol=1e-8)
    h = 1e-4
    plus = {k: x[k]+h*direction[k] for k in x}
    minus = {k: x[k]-h*direction[k] for k in x}
    evaluate = jax.jit(scalar)
    fd = (evaluate(plus)-evaluate(minus))/(2*h)
    np.testing.assert_allclose(reverse, fd, rtol=3e-5, atol=1e-6)
    cache_size = compiled._cache_size()
    compiled(plus)
    assert compiled._cache_size() == cache_size
    # Exact zero launch energy is a clipping boundary. Its chosen derivative
    # must remain finite even though inactive path formulas contain divisions.
    at_zero = dict(x, em=jnp.zeros_like(x['em']))
    (_, error), gradient = compiled(at_zero)
    assert not error.is_fatal()
    assert all(np.isfinite(g).all() for g in gradient.values())


@pytest.mark.parametrize('method', [1])
def test_parcel_orientation_stretched_grid_host_floor_and_cap(method):
    """Compare physical lengths after native grid reversal, enforce host caps, and check
    parameter-column isolation.
    """
    from clubb_jax.src.CLUBB_core.grid_class import GRID_TYPE_STRETCHED_ZM
    from clubb_jax.src.CLUBB_core.mixing_length import calc_Lscale_directly
    heights = np.array([[100., 150., 230., 340., 500., 720., 1000.],
                        [100., 160., 250., 380., 550., 760., 1000.]])
    gr = setup_grid(ngrdcol=2, deltaz=80., zm_init=100., zm_top=1000.,
                    grid_type=GRID_TYPE_STRETCHED_ZM, momentum_heights=heights)
    thl = 298. + .003*gr.zt
    rt = jnp.full_like(thl, .008)
    pdf = init_pdf_params(gr.nzt, 2).replace(thl_1=thl+.5, thl_2=thl-.3, rt_1=rt+.001, rt_2=rt-.001)
    cp = jnp.asarray(init_clubb_params(2, ''))
    em = jnp.broadcast_to(jnp.array([[.01], [1.]]), (2, gr.nzm))
    stats = JaxStats.empty(l_sample=False, names=(), ncol=2, max_nlev=gr.nzt)

    def calculate(grid, thl, rt, em, pdf, cp):
        """Return legacy outputs on the supplied native grid with host floors
        and per-column caps.
        """
        exner = jnp.exp(-3e-5*grid.zt)
        return calc_Lscale_directly(2, gr.nzm, gr.nzt, grid, True,
            1e5*exner**3.5, exner, rt, thl, thl*(1.+ep1*rt),
            jnp.full_like(rt, 1e-6), jnp.full_like(thl, .1), jnp.zeros_like(rt),
            pdf, em, thl, jnp.array([15., 70.]), 20., cp, 3, False,
            method, stats, ErrInfo.initialized(2))

    ordinary = jax.jit(calculate)(gr, thl, rt, em, pdf, cp)
    # Grid construction currently gates descending standalone runs. Construct
    # the physically identical native layout here to check the kernel boundary.
    reverse = gr.replace(
        **{name: -getattr(gr, name)[:, ::-1] for name in ('dzm', 'dzt', 'invrs_dzm', 'invrs_dzt')},
        zm=gr.zm[:, ::-1], zt=gr.zt[:, ::-1],
        weights_zm2zt=gr.weights_zm2zt[:, ::-1, ::-1],
        weights_zt2zm=gr.weights_zt2zm[:, ::-1, ::-1],
        grid_dir=-1., grid_dir_indx=-1, k_lb_zm=gr.nzm-1, k_ub_zm=0,
        k_lb_zt=gr.nzt-1, k_ub_zt=0)
    pdf_reverse = pdf.replace(**{name: getattr(pdf, name)[:, ::-1]
                                for name in ('thl_1', 'thl_2', 'rt_1', 'rt_2')})
    flipped = jax.jit(calculate)(reverse, thl[:, ::-1], rt[:, ::-1], em[:, ::-1], pdf_reverse, cp)
    assert not ordinary[0].is_fatal() and not flipped[0].is_fatal()
    for result, reversed_result in zip(ordinary[1:4], flipped[1:4]):
        assert result.shape == reversed_result.shape == (2, gr.nzt)
        np.testing.assert_allclose(result, reversed_result[:, ::-1], rtol=1e-11, atol=1e-9)
    assert np.all(np.asarray(ordinary[1]) <= np.array([[15.], [70.]]))
    # Altering one column's runtime parameters cannot change its neighbor.
    changed_cp = cp.at[0, PNAME_IDX['Lscale_mu_entrain']].set(.01)
    changed = jax.jit(calculate)(gr, thl, rt, em, pdf, changed_cp)
    for result, neighbor in zip(changed[1:4], ordinary[1:4]):
        np.testing.assert_array_equal(result[1], neighbor[1])


def test_legacy_invalid_entrainment_preserves_healthy_column(capfd):
    """Inject zero entrainment in one column and require its fatal diagnostic plus finite,
    unchanged neighboring results.
    """
    gr, exner, pressure, thvm, _, x = fixture()
    def calculate(mu):
        """Return legacy parcel lengths/errors for the supplied two-column entrainment vector."""
        return compute_mixing_length(gr.nzm, gr.nzt, 2, gr, thvm, x['thlm'],
            x['rtm'], x['em'], jnp.full((2,), 1e5), pressure, exner, x['thlm'],
            mu, 20., 3, False, ErrInfo.initialized(2))
    run = jax.jit(calculate)
    healthy = run(jnp.array([.001, .002]))
    invalid = run(jnp.array([0., .002]))
    np.testing.assert_array_equal(invalid[0].fatal_mask(), [True, False])
    jax.effects_barrier()
    assert 'Entrainment rate Lscale_mu_entrain cannot be 0' in capfd.readouterr().err
    for result, reference in zip(invalid[1:], healthy[1:]):
        assert np.isfinite(result).all()
        np.testing.assert_array_equal(result[1], reference[1])


def test_cpu_launch_groups_match_single_columns_and_transpose():
    """Compare partial legacy launch groups with singleton columns, fatal masks, and
    forward/reverse/finite-difference derivatives.
    """
    # Nine columns and ten/eleven active launch levels exercise partial final
    # groups in both directions. Mixed entrainment values expose cross-column
    # leakage; column 0 checks fatal-status preservation through grouped scans.
    gr = setup_grid(ngrdcol=9, deltaz=80., zm_init=0., zm_top=960.)
    mu = jnp.linspace(.0007, .002, 9).at[0].set(0.)
    def calculate(gr, mu):
        """Return legacy outputs on the supplied grid for a mixed valid/invalid
        entrainment vector.
        """
        thl = 298. + .003*gr.zt + .02*jnp.sin(gr.zt)
        rt = jnp.full_like(thl, .008)
        exner = jnp.exp(-3e-5*gr.zt)
        return compute_mixing_length(gr.nzm, gr.nzt, gr.ngrdcol, gr,
            thl*(1.+ep1*rt), thl, rt, jnp.full((gr.ngrdcol, gr.nzm), .4),
            jnp.full((gr.ngrdcol,), 1e5), 1e5*exner**3.5, exner, thl,
            mu, 0., 3, False, ErrInfo.initialized(gr.ngrdcol))
    run = jax.jit(calculate)
    batch = run(gr, mu)
    np.testing.assert_array_equal(batch[0].fatal_mask(), [True]+[False]*8)
    for i in range(9):
        one_grid = jax.tree.map(lambda a: a[i:i+1], gr).replace(ngrdcol=1)
        single = run(one_grid, mu[i:i+1])
        for result, reference in zip(batch[1:], single[1:]):
            assert reference.shape == (1, gr.nzt)
            np.testing.assert_allclose(result[i:i+1], reference, rtol=2e-11, atol=2e-10)
    scalar = lambda m: jnp.sum(calculate(gr, m)[1])
    gradient = jax.jit(jax.grad(scalar))(mu)
    assert np.isfinite(gradient).all()
    direction = jnp.cos(jnp.arange(9))*1e-5
    direction = direction.at[0].set(0.)
    _, tangent = jax.jit(lambda m, d: jax.jvp(scalar, (m,), (d,)))(mu, direction)
    reverse = jnp.vdot(gradient, direction)
    np.testing.assert_allclose(reverse, tangent, rtol=1e-9, atol=1e-9)
    evaluate = jax.jit(scalar)
    h = 1e-4
    finite_difference = (evaluate(mu+h*direction)-evaluate(mu-h*direction))/(2*h)
    np.testing.assert_allclose(reverse, finite_difference, rtol=2e-5, atol=1e-7)
