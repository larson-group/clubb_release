"""Check JAX parcel stopping, coherent tails, and derivative behavior.

Run these generated checks through the repository owner:
    bash human_reviewed_tests/run_pytests.sh -jax -k mixing_length

The coherent reference uses NumPy parcel equations with scalar launch/target
loops and actual early exits, sharing only the native grid interpolation.
The legacy first-layer reference is dry and linearly stratified, with an
analytic stopping distance. Other checks cover all three runtime methods,
zero energy, finite/JVP/reverse/finite-difference derivatives, compilation
reuse, grid reversal, column isolation, host caps, and incomplete groups.
Each test documents its setup and assertions below; fixtures restore debug
state. No Fortran/F2PY runtime is imported. Whole-driver/Fortran comparisons
and real-case gradients belong to the separate integration_tests/ workflows.
"""
import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core import error_code
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv, Rd, ep, ep1, ep2, grav, zero_threshold
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.grid_class import setup_grid, zm2zt
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.mixing_length import (
    calc_Lscale, compute_mixing_length, compute_coherent_mixing_length, sat_mixrat_liq_flatau,
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


def flatau_numpy(pressure, temperature):
    """Return a NumPy polynomial saturation reference [kg/kg] for pressure [Pa] and
    temperature [K].
    """
    # Use NumPy's polynomial evaluator as the independent saturation reference.
    coefficients = [6.11583699, .444606896, .0143177157, .000264224321,
                    .00000299291081, .0000000203154182, .0000000000702620698,
                    .0000000000000379534310, -.000000000000000321582393]
    t = np.maximum(np.asarray(temperature) - 273.15, -85.)
    vapor_pressure = 100. * np.polynomial.polynomial.polyval(t, coefficients)
    return ep * vapor_pressure / np.where(pressure - vapor_pressure < 1., vapor_pressure,
                                         pressure - vapor_pressure)


def test_local_flatau_temperature_floor_and_pressure_cap():
    """Compare clipped/capped saturation values with NumPy and require finite, zero
    derivatives on the flat branches.
    """
    temperature = jnp.array([100., 187., 250., 273.15, 300., 350.])
    pressure = jnp.array([1e5, 1e5, 1e5, 1e5, 1e5, 10.])
    actual = jax.jit(sat_mixrat_liq_flatau)(pressure, temperature)
    expected = flatau_numpy(np.asarray(pressure), np.asarray(temperature))
    assert actual.shape == expected.shape
    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=1e-15)
    gradient = jax.jit(jax.grad(lambda t: jnp.sum(sat_mixrat_liq_flatau(pressure, t))))(temperature)
    assert np.isfinite(gradient).all()
    np.testing.assert_array_equal(gradient[jnp.array([0, 1, 5])], 0.)


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
    cp = jnp.asarray(init_clubb_params(2, '')).at[:, PNAME_IDX['Lscale_pdf_launch_frac']].set(.35)
    cp = cp.at[:, PNAME_IDX['Lscale_up_weight']].set(jnp.array([.25, .8]))
    inputs = dict(thlm=thlm, rtm=rt, em=jnp.full((2, gr.nzm), .4), clubb_params=cp)
    return gr, exner, 1e5 * exner**3.5, thlm * (1. + ep1 * rt), pdf, inputs


def coherent_numpy(gr, thvm, thlm, rt, pdf, em, pressure, exner, mu, multiplier, frac, weight):
    """Return serial coherent lengths and the number of actual tail cutoffs.

    Use ascending grid fields with the same units as the production routine;
    mu, multiplier, fraction, and weight are per-column vectors. Reconstruct
    arrival-weighted reaches without the production scan/grouping control flow.
    The surface floor and cap are omitted; tests choose a zero floor and a
    nonbinding cap for these comparisons.
    """
    z = np.asarray(gr.zt)
    ncol, nlev = z.shape
    energy = np.maximum(np.asarray(zm2zt(gr.nzm, gr.nzt, ncol, gr, em)), 0.) * multiplier[:, None]
    directional = []
    cutoff_count = 0
    for sign, component in [(1, 1), (-1, 2)]:
        raw = np.zeros_like(z)
        crossed = np.zeros_like(z)
        arrival = np.ones_like(z)
        for i in range(ncol):
            def buoyancy(thl_parcel, rt_parcel, j):
                """Return path-oriented moist buoyancy [m/s^2] at target j for one column."""
                t = thl_parcel * exner[i, j]
                rs = flatau_numpy(pressure[i, j], t)
                rc = max((rt_parcel - rs) * t*t / (t*t + ep*Lv*Lv/(Rd*Cp)*rs), zero_threshold)
                thv = thl_parcel + ep1 * thlm[i, j] * rt_parcel + (Lv/(exner[i, j]*Cp) - ep2*thlm[i, j])*rc
                return sign * (grav/thvm[i, j]*thv - grav)
            for k in range(nlev):
                if energy[i, k] <= 0.:
                    continue
                thl_parcel = thlm[i, k] + frac[i] * (np.asarray(getattr(pdf, f'thl_{component}'))[i, k] - thlm[i, k])
                rt_parcel = rt[i, k] + frac[i] * (np.asarray(getattr(pdf, f'rt_{component}'))[i, k] - rt[i, k])
                previous_buoyancy = buoyancy(thl_parcel, rt_parcel, k)
                work = barrier = 0.
                targets = range(k+1, nlev-1) if sign == 1 else range(k-1, -1, -1)
                for j in targets:
                    previous = j-sign
                    depth = abs(z[i, j] - z[i, previous])
                    retention = np.exp(-mu[i]*depth)
                    response = (1.-retention)/(mu[i]*depth)
                    thl_parcel = thlm[i, j] + retention*(thl_parcel-thlm[i, previous]) - response*(thlm[i, j]-thlm[i, previous])
                    rt_parcel = rt[i, j] + retention*(rt_parcel-rt[i, previous]) - response*(rt[i, j]-rt[i, previous])
                    current_buoyancy = buoyancy(thl_parcel, rt_parcel, j)
                    work += .5*(previous_buoyancy + current_buoyancy)*depth
                    next_barrier = max(barrier, -work)
                    coherence = energy[i, k]**2 / ((energy[i, k]+barrier)*(energy[i, k]+next_barrier))
                    length = coherence*depth
                    raw[i, k] += length
                    if sign == 1 or k != nlev-1:
                        crossed[i, j] += length
                        arrival[i, j] += energy[i, k]**2/(energy[i, k]+next_barrier)**2
                    if coherence <= .01:
                        cutoff_count += 1
                        break
                    previous_buoyancy, barrier = current_buoyancy, next_barrier
            order = range(nlev-1) if sign == 1 else range(nlev-2, -1, -1)
            remaining = raw[i, 0 if sign == 1 else nlev-2]
            for position, k in enumerate(order):
                if position:
                    remaining += raw[i, k] - crossed[i, k]
                raw[i, k] = max(.1, remaining/arrival[i, k])
        directional.append(raw)
    up, down = directional
    result = up[:, :-1]**weight[:, None] * down[:, :-1]**(1.-weight[:, None])
    return (np.concatenate((result, result[:, -1:]), axis=1), up, down), cutoff_count


@pytest.mark.parametrize('energy', [.4, .0001, 0.])
def test_coherent_matches_serial_early_exit_oracle(energy):
    """Compare both moist columns with the serial oracle at ordinary, cutoff-triggering,
    and zero launch energies.
    """
    gr, exner, pressure, thvm, pdf, x = fixture()
    em = jnp.full_like(x['em'], energy)
    mu, multiplier = jnp.array([.001, .002]), jnp.array([1., .7])
    frac, weight = jnp.array([.35, 1.]), jnp.array([0., 1.])
    actual = jax.jit(lambda e: compute_coherent_mixing_length(
        gr.nzm, gr.nzt, 2, gr, thvm, x['thlm'], x['rtm'], pdf, e,
        jnp.full((2,), 1e5), pressure, exner, x['thlm'], mu, multiplier,
        frac, weight, 0., False, ErrInfo.initialized(2)))(em)
    expected, cutoffs = coherent_numpy(gr, np.asarray(thvm), np.asarray(x['thlm']),
        np.asarray(x['rtm']), pdf, em, np.asarray(pressure), np.asarray(exner),
        np.asarray(mu), np.asarray(multiplier), np.asarray(frac), np.asarray(weight))
    assert not actual[0].is_fatal()
    for result, reference in zip(actual[1:], expected):
        assert result.shape == reference.shape == (2, gr.nzt)
        np.testing.assert_allclose(result, reference, rtol=2e-10, atol=2e-9)
    if energy == .0001:
        assert cutoffs > 0  # This oracle actually exercises the masked tail.


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


@pytest.mark.parametrize('method', [1, 2, 3])
def test_three_methods_gradients_and_runtime_parameter_reuse(method):
    """Check each method's gradients and compilation reuse on two moist columns.

    Require finite zero-energy derivatives, JVP/reverse/finite-difference
    agreement at smooth inputs, and nonzero coherent-parameter sensitivities.
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
    for name, amount in [('Lscale_mu_entrain', 1e-5), ('Lscale_tke_multiplier', .02),
                         ('Lscale_pdf_launch_frac', .01), ('Lscale_up_weight', .01)]:
        direction['clubb_params'] = direction['clubb_params'].at[:, PNAME_IDX[name]].set(amount)
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
    if method == 3:
        for name in ('Lscale_mu_entrain', 'Lscale_tke_multiplier', 'Lscale_pdf_launch_frac', 'Lscale_up_weight'):
            assert np.any(np.asarray(gradient['clubb_params'][:, PNAME_IDX[name]]) != 0.), name

    # Exact zero launch energy is a clipping boundary. Its chosen derivative
    # must remain finite even though inactive path formulas contain divisions.
    at_zero = dict(x, em=jnp.zeros_like(x['em']))
    (_, error), gradient = compiled(at_zero)
    assert not error.is_fatal()
    assert all(np.isfinite(g).all() for g in gradient.values())


@pytest.mark.parametrize('method', [1, 3])
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
    cp = jnp.asarray(init_clubb_params(2, '')).at[:, PNAME_IDX['Lscale_pdf_launch_frac']].set(.3)
    em = jnp.broadcast_to(jnp.array([[.01], [1.]]), (2, gr.nzm))
    stats = JaxStats.empty(l_sample=False, names=(), ncol=2, max_nlev=gr.nzt)

    def calculate(grid, thl, rt, em, pdf, cp):
        """Return direct/coherent outputs on the supplied native grid with host floors
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


@pytest.mark.parametrize('descending', [False, True])
def test_coherent_partial_launch_groups_preserve_columns_and_derivatives(descending):
    """Compare a nine-column coherent batch with serial and singleton references.

    Use partial launch/offset blocks, mixed moisture and energy, real early
    cutoffs, both grid orientations, column-isolated changes, and JVP/reverse/
    finite-difference checks to detect padding or masking errors.
    """
    # A partial final launch group and offset block must not contribute padded
    # parcels. Mixed moisture/energy columns include zero energy, early cutoffs,
    # and long trajectories; the serial NumPy oracle uses actual early exits.
    ncol = 9
    gr = setup_grid(ngrdcol=ncol, deltaz=80., zm_init=0., zm_top=2960.)
    thlm = 298. + .003*gr.zt + jnp.arange(ncol)[:, None]*.03
    rt = jnp.broadcast_to(jnp.linspace(.008, .018, ncol)[:, None], thlm.shape)
    em = jnp.broadcast_to(jnp.linspace(.1, .8, ncol)[:, None], (ncol, gr.nzm))
    em = em.at[0].set(0.).at[1].set(.0001)
    cp = jnp.asarray(init_clubb_params(ncol, ''))
    cp = cp.at[:, PNAME_IDX['Lscale_mu_entrain']].set(jnp.linspace(.0007, .002, ncol))
    cp = cp.at[:, PNAME_IDX['Lscale_pdf_launch_frac']].set(jnp.linspace(0., 1., ncol))
    cp = cp.at[:, PNAME_IDX['Lscale_up_weight']].set(jnp.linspace(0., 1., ncol))
    exner = jnp.exp(-3e-5*gr.zt)
    pdf = init_pdf_params(gr.nzt, ncol).replace(
        thl_1=thlm+.7, thl_2=thlm-.4, rt_1=rt+.001, rt_2=rt-.0005)
    expected, cutoffs = coherent_numpy(gr, np.asarray(thlm*(1.+ep1*rt)),
        np.asarray(thlm), np.asarray(rt), pdf, np.asarray(em),
        np.asarray(1e5*exner**3.5), np.asarray(exner),
        np.asarray(cp[:, PNAME_IDX['Lscale_mu_entrain']]),
        np.asarray(cp[:, PNAME_IDX['Lscale_tke_multiplier']]),
        np.asarray(cp[:, PNAME_IDX['Lscale_pdf_launch_frac']]),
        np.asarray(cp[:, PNAME_IDX['Lscale_up_weight']]))
    assert cutoffs > 0
    if descending:
        gr = gr.replace(
            **{name: -getattr(gr, name)[:, ::-1]
               for name in ('dzm', 'dzt', 'invrs_dzm', 'invrs_dzt')},
            zm=gr.zm[:, ::-1], zt=gr.zt[:, ::-1],
            weights_zm2zt=gr.weights_zm2zt[:, ::-1, ::-1],
            weights_zt2zm=gr.weights_zt2zm[:, ::-1, ::-1],
            grid_dir=-1., grid_dir_indx=-1, k_lb_zm=gr.nzm-1, k_ub_zm=0,
            k_lb_zt=gr.nzt-1, k_ub_zt=0)
        thlm, rt, em = thlm[:, ::-1], rt[:, ::-1], em[:, ::-1]
        expected = [a[:, ::-1] for a in expected]
    x = dict(thlm=thlm, rt=rt, em=em, params=cp)

    def calculate(grid, values):
        """Return coherent outputs from dynamic state/parameters on the supplied native grid."""
        ng = grid.ngrdcol
        thl, rt = values['thlm'], values['rt']
        exner = jnp.exp(-3e-5*grid.zt)
        params = values['params']
        pdf = init_pdf_params(grid.nzt, ng).replace(
            thl_1=thl+.7, thl_2=thl-.4, rt_1=rt+.001, rt_2=rt-.0005)
        return compute_coherent_mixing_length(grid.nzm, grid.nzt, ng, grid,
            thl*(1.+ep1*rt), thl, rt, pdf, values['em'], jnp.full((ng,), 1e5),
            1e5*exner**3.5, exner, thl,
            params[:, PNAME_IDX['Lscale_mu_entrain']],
            params[:, PNAME_IDX['Lscale_tke_multiplier']],
            params[:, PNAME_IDX['Lscale_pdf_launch_frac']],
            params[:, PNAME_IDX['Lscale_up_weight']], 0., False, ErrInfo.initialized(ng))

    run = jax.jit(calculate)
    batch = run(gr, x)
    assert not batch[0].is_fatal()
    for actual, reference in zip(batch[1:], expected):
        assert actual.shape == reference.shape == (ncol, gr.nzt)
        np.testing.assert_allclose(actual, reference, rtol=2e-10, atol=2e-9)
    for i in range(ncol):
        one_grid = jax.tree.map(lambda a: a[i:i+1], gr).replace(ngrdcol=1)
        single = run(one_grid, {k: a[i:i+1] for k, a in x.items()})
        for actual, reference in zip(batch[1:], single[1:]):
            assert reference.shape == (1, gr.nzt)
            np.testing.assert_allclose(actual[i:i+1], reference, rtol=2e-10, atol=2e-9)
    changed = run(gr, dict(x, params=cp.at[4, PNAME_IDX['Lscale_mu_entrain']].set(.003)))
    for actual, reference in zip(changed[1:], batch[1:]):
        np.testing.assert_array_equal(np.delete(actual, 4, axis=0), np.delete(reference, 4, axis=0))

    loss = lambda values: jnp.sum(calculate(gr, values)[1] * jnp.linspace(.5, 1.5, gr.nzt))
    gradient = jax.jit(jax.grad(loss))(x)
    assert all(np.isfinite(g).all() for g in gradient.values())
    direction = dict(thlm=jnp.cos(gr.zt)*.01, rt=jnp.ones_like(rt)*1e-6,
                     em=jnp.ones_like(em).at[0].set(0.)*.001, params=jnp.zeros_like(cp))
    direction['params'] = direction['params'].at[:, PNAME_IDX['Lscale_mu_entrain']].set(1e-5)
    _, tangent = jax.jit(lambda x, d: jax.jvp(loss, (x,), (d,)))(x, direction)
    reverse = sum(jnp.vdot(gradient[k], direction[k]) for k in x)
    np.testing.assert_allclose(tangent, reverse, rtol=2e-9, atol=1e-8)
    evaluate = jax.jit(loss)
    h = 1e-4
    plus = {k: x[k]+h*direction[k] for k in x}
    minus = {k: x[k]-h*direction[k] for k in x}
    finite_difference = (evaluate(plus)-evaluate(minus))/(2*h)
    np.testing.assert_allclose(reverse, finite_difference, rtol=5e-5, atol=1e-6)
