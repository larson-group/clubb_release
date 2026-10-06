"""Check source-derived reset, parameter and absolute-window contracts.

The lifecycle follows src/clubb_driver.F90 (init_clubb_case,
set_case_initial_conditions and advance_clubb_to_end); derived parameters
follow src/CLUBB_core/parameters_tunable.F90. These focused checks cover state
snapshots, parameter refresh, batch selection and window validation. The real
reset-and-rerun scenario adapted from src/clubb_driver_test.F90 lives in
clubb_jax/src/clubb_driver_test.py; additional window/batch checks live in
clubb_jax/tests/run_driver_extensions_test.py.
"""
from types import SimpleNamespace
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src import advance_clubb_to_end as driver
from clubb_jax.src import clubb_case_initalization as initialization
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core import error_code
from clubb_jax.src.CLUBB_core.parameters_tunable import (
    PNAME_IDX, calc_derived_params, init_clubb_params,
)
from clubb_jax.src.CLUBB_core.stats_netcdf import StatsWriter
from clubb_jax.src.Input_fields import namelist
from utilities.create_case_namelist import (
    create_case_namelist_file, prune_clubb_stats_namelist, set_stats_string,
)


def _grid():
    return SimpleNamespace(
        nzt=3, nzm=4, k_lb_zt=0, k_ub_zt=2, k_lb_zm=0, k_ub_zm=3,
        zt=jnp.array([[10., 30., 50.], [20., 100., 180.]]),
        zm=jnp.array([[0., 20., 40., 60.], [0., 80., 160., 240.]]),
    )


def _writer(tmp_path, *, total_columns=2):
    tmp_path.mkdir(parents=True, exist_ok=True)
    registry = tmp_path / 'stats.in'
    registry.write_text(
        '&clubb_stats_nl\n entry(1)="thlm | zt | K | Potential temperature"\n/\n',
    )
    return StatsWriter(
        registry_path=str(registry), output_path='', nzt=3, nzm=4, ngrdcol=2,
        zt=np.array([10., 30., 50.]), zm=np.array([0., 20., 40., 60.]),
        stats_tsamp=60., stats_tout=120., dt_main=60.,
        day=1, month=1, year=2000, time_initial=0., ncol_total=total_columns,
    )


@pytest.fixture
def lifecycle_state(tmp_path, monkeypatch):
    monkeypatch.setattr(error_code, '_debug_level', -1)
    params = jnp.asarray(init_clubb_params(4, ''))
    params = params.at[:, PNAME_IDX['C8']].set(jnp.array([.2, .4, .6, .8]))
    state = dict(
        ngrdcol=2, total_param_sets=4, clubb_params_all=params,
        clubb_params=params[:2], gr=_grid(),
        cfg=dict(grid_type=1, deltaz_nl=20.),
        flags=SimpleNamespace(l_prescribed_avg_deltaz=False),
        stats_writer=_writer(tmp_path, total_columns=4), err_info=ErrInfo.initialized(2),
        iinit=4, thlm=jnp.ones((2, 3)),
        sampling_state=SimpleNamespace(prior_iter=3, perm=jnp.arange(4)),
    )
    state['_initial_state'] = dict(state)
    return state


@pytest.mark.parametrize('assignment,expected', [
    ('C8=0.25', [.25, .5, .5, .5]),
    ('C8=0.25,0.75', [.25, .75, .5, .5]),
    ('C8=,0.75,,0.125', [.5, .75, .5, .125]),
    ('C8(2)=0.25', [.5, .25, .5, .5]),
    ('C8(2:3)=0.25,0.75', [.5, .25, .75, .5]),
    ('C8(1:4:2)=0.25,0.75', [.25, .5, .75, .5]),
    ('C8(4:1:-1)=0.2,0.4,0.6,0.8', [.8, .6, .4, .2]),
    ('C8(4:1:-2)=0.25,0.75', [.5, .75, .5, .25]),
    ('C8 ( 4 : 1 : -1 )=0.2,,0.6', [.5, .6, .5, .2]),
])
def test_parameter_namelist_assigns_source_selected_elements(tmp_path, monkeypatch, assignment, expected):
    # Exercise the parser used in the provisioned runtime without f90nml.
    monkeypatch.setattr(namelist, 'f90nml', None)
    path = tmp_path / 'parameters.in'
    path.write_text(f'&clubb_params_nl\n {assignment}\n/\n')
    params = init_clubb_params(4, str(path))
    np.testing.assert_array_equal(params[:, PNAME_IDX['C8']], expected)


def test_parameter_defaults_and_assignment_bounds(tmp_path):
    defaults = init_clubb_params(2, '')
    np.testing.assert_array_equal(defaults[0], defaults[1])
    assert defaults[0, PNAME_IDX['C8']] == .5
    path = tmp_path / 'overflow.in'
    path.write_text('&clubb_params_nl\n C8=.2,.4,.6\n/\n')
    with pytest.raises(ValueError, match='exceeds ngrdcol'):
        init_clubb_params(2, str(path))


@pytest.mark.parametrize('section,message', [
    ('1:4:0', 'stride must not be zero'),
    (':4', 'Unsupported parameter section'),
    ('0:3', 'exceeds ngrdcol'),
])
def test_parameter_sections_reject_invalid_or_unhandled_syntax(tmp_path, monkeypatch, section, message):
    monkeypatch.setattr(namelist, 'f90nml', None)
    path = tmp_path / 'parameters.in'
    path.write_text(f'&clubb_params_nl\n C8({section})=.2\n/\n')
    with pytest.raises(ValueError, match=message):
        init_clubb_params(4, str(path))


@pytest.mark.parametrize('grid_type,prescribed,expected_spacing', [
    (1, False, [20., 40.]),
    (2, False, [20., 80.]),
    (3, False, [20., 80.]),
    (3, True, [20., 40.]),
])
def test_derived_parameters_follow_source_grid_branches(grid_type, prescribed, expected_spacing):
    params = jnp.asarray(init_clubb_params(2, ''))
    params = params.at[:, PNAME_IDX['mult_coef']].set(jnp.array([.25, .75]))
    params = params.at[:, PNAME_IDX['lmin_coef']].set(jnp.array([.5, .6]))
    params = params.at[:, PNAME_IDX['Skw_max_mag']].set(jnp.array([8., 4.]))
    nu, lmin, mixture_bound = calc_derived_params(
        _grid(), 2, grid_type, jnp.array([20., 40.]), params, prescribed,
    )
    spacing = np.asarray(expected_spacing)
    factor = np.where(spacing > 40., 1. + np.array([.25, .75]) * np.log(spacing / 40.), 1.)
    for name in ('nu1', 'nu2', 'nu6', 'nu8', 'nu9', 'nu10', 'nu_hm'):
        np.testing.assert_allclose(getattr(nu, name), params[:, PNAME_IDX[name]] * factor, rtol=1.e-14)
    assert float(lmin) == 24.
    np.testing.assert_allclose(mixture_bound, 1. - .5 * (1. - 4. / np.sqrt(4. * (1. - .4)**3 + 4.**2)))


def test_derived_parameters_compile_and_differentiate():
    params = jnp.asarray(init_clubb_params(2, ''))

    def value(parameters):
        nu, lmin, _ = calc_derived_params(
            _grid(), 2, 3, jnp.array([20., 40.]), parameters, False,
        )
        return jnp.sum(nu.nu1) + lmin

    gradient = jax.jit(jax.grad(value))(params)
    assert np.isfinite(gradient).all()
    np.testing.assert_allclose(gradient[:, PNAME_IDX['nu1']], [1., 1. + .5 * np.log(2.)])
    np.testing.assert_array_equal(gradient[:, PNAME_IDX['lmin_coef']], [0., 40.])


def test_reset_restores_fields_sampling_state_and_restart_clock(lifecycle_state):
    state = lifecycle_state
    initial = state['_initial_state']
    state['thlm'] = state['thlm'] + 20.
    state['sampling_state'] = SimpleNamespace(prior_iter=100, perm=jnp.zeros(4))
    state['iinit'] = 101
    state['err_info'] = state['err_info'].set_fatal()
    state['_jax_stats'] = object()
    initialization.set_case_initial_conditions(state)
    np.testing.assert_array_equal(state['thlm'], initial['thlm'])
    assert state['sampling_state'] is initial['sampling_state']
    assert state['iinit'] == 4
    assert not state['err_info'].is_fatal()
    assert '_jax_stats' not in state


def test_reset_retains_current_parameters_and_recomputes_derived_state(lifecycle_state):
    state = lifecycle_state
    changed = state['clubb_params'].at[:, PNAME_IDX['nu1']].set(3.)
    changed = changed.at[:, PNAME_IDX['lmin_coef']].set(.6)
    initialization.set_case_initial_conditions(state, changed)
    state['thlm'] = state['thlm'] + 1.
    initialization.set_case_initial_conditions(state)
    np.testing.assert_array_equal(state['clubb_params'], changed)
    np.testing.assert_array_equal(state['nu_vert_res_dep'].nu1, [3., 3.])
    assert float(state['lmin']) == 24.


def test_batch_selection_and_repeated_batch_cycle(lifecycle_state):
    state = lifecycle_state
    writer = state['stats_writer']
    for batch_num, offset, expected in [(1, 0, [.2, .4]), (2, 2, [.6, .8]), (1, 0, [.2, .4])]:
        initialization.set_case_initial_conditions(state, batch_num=batch_num)
        np.testing.assert_array_equal(state['clubb_params'][:, PNAME_IDX['C8']], expected)
        assert writer.active_batch_offset == offset
        assert writer._time_index == 0
        assert not writer.l_sample


@pytest.mark.parametrize('batch_num', [0, 3])
def test_reset_rejects_invalid_batch_number(lifecycle_state, batch_num):
    with pytest.raises(ValueError, match='batch_num'):
        initialization.set_case_initial_conditions(lifecycle_state, batch_num=batch_num)


def test_explicit_parameters_take_precedence_and_keep_runtime_shape(lifecycle_state):
    state = lifecycle_state
    changed = state['clubb_params'].at[:, PNAME_IDX['C8']].set(.9)
    initialization.set_case_initial_conditions(state, changed, batch_num=2)
    np.testing.assert_array_equal(state['clubb_params'], changed)
    with pytest.raises(ValueError, match='shape'):
        initialization.set_case_initial_conditions(state, changed[:1])


@pytest.mark.parametrize('debug_level,raises', [(-1, False), (0, False), (1, True), (2, True)])
def test_reset_validates_replacement_parameters_at_source_debug_level(
    lifecycle_state, monkeypatch, debug_level, raises,
):
    monkeypatch.setattr(error_code, '_debug_level', debug_level)
    state = lifecycle_state
    changed = state['clubb_params'].at[:, PNAME_IDX['C_uu_shr']].set(1.5)
    if raises:
        with pytest.raises(RuntimeError, match='check_clubb_settings'):
            initialization.set_case_initial_conditions(state, changed)
    else:
        initialization.set_case_initial_conditions(state, changed)


def _driver_state(tmp_path):
    state = dict(
        dt_main=60., dt_rad=120., time_initial=0., time_final=360.,
        iinit=1, ifinal=6, l_stats=True, stats_writer=_writer(tmp_path),
        ngrdcol=2, nzt=3, nzm=4, l_calc_thlp2_rad=False,
        err_info=ErrInfo.initialized(2),
    )
    for name in (
        'thlm', 'rtm', 'rcm', 'exner', 'thv_ds_zt', 'thlm_forcing',
        'rtm_forcing', 'radht', 'rcm_mc', 'rvm_mc', 'thlm_mc',
        'wprtp_forcing', 'wpthlp_forcing', 'rtp2_forcing', 'thlp2_forcing',
        'rtpthlp_forcing', 'wprtp_mc', 'wpthlp_mc', 'rtp2_mc', 'thlp2_mc', 'rtpthlp_mc',
    ):
        state[name] = jnp.zeros((2, 3))
    return state


@pytest.fixture
def patched_physics(monkeypatch):
    monkeypatch.setattr(error_code, '_debug_level', -1)
    calls = []
    monkeypatch.setattr(driver, 'calculate_thvm', lambda **kwargs: kwargs['thlm'])
    monkeypatch.setattr(driver, '_prescribe_forcings', lambda state, time: None)
    monkeypatch.setattr(driver, '_advance_radiation', lambda **kwargs: None)

    def core(state):
        state['thlm'] = state['thlm'] + 1.
        state['_jax_stats'] = state['_jax_stats'].update('thlm', state['thlm'])

    def microphysics(state, itime, time_current, l_rad_itime):
        calls.append((itime, time_current, l_rad_itime))

    monkeypatch.setattr(driver, '_advance_clubb_core', core)
    monkeypatch.setattr(driver, '_advance_microphysics', microphysics)
    return calls


@pytest.mark.parametrize('iinit,options,expected', [
    (3, {}, [3, 4, 5, 6]),
    (3, {'max_steps': 2}, [3, 4]),
    (3, {'itime_end': 4}, [3, 4]),
    (3, {'itime_start': 4, 'itime_end': 5}, [4, 5]),
    (1, {'itime_start': 3, 'max_steps': 2}, [3, 4]),
    (1, {'itime_start': 5, 'itime_end': 4}, []),
])
def test_absolute_driver_window_indices(tmp_path, patched_physics, iinit, options, expected):
    state = _driver_state(tmp_path)
    state['iinit'] = iinit
    driver.advance_clubb_to_end(state, l_stdout=False, **options)
    assert patched_physics == [(itime, (itime - 1) * 60., itime % 2 == 0 or itime == 1) for itime in expected]


def test_successive_windows_preserve_model_and_stats_state(tmp_path, patched_physics):
    whole = _driver_state(tmp_path / 'whole')
    windows = _driver_state(tmp_path / 'windows')
    driver.advance_clubb_to_end(whole, l_stdout=False, itime_end=4)
    driver.advance_clubb_to_end(windows, l_stdout=False, itime_end=2)
    driver.advance_clubb_to_end(windows, l_stdout=False, itime_start=3, itime_end=4)
    np.testing.assert_array_equal(windows['thlm'], whole['thlm'])
    for first, second in zip(windows['_jax_stats'].buffers, whole['_jax_stats'].buffers):
        np.testing.assert_array_equal(first, second)
    for first, second in zip(windows['_jax_stats'].nsamples, whole['_jax_stats'].nsamples):
        np.testing.assert_array_equal(first, second)
    assert windows['stats_writer']._time_index == whole['stats_writer']._time_index == 2


def test_source_stats_suppression_preserves_layout_for_later_window(tmp_path, patched_physics):
    # The native optional argument suppresses this advance, not the configured
    # writer. A later window must retain its registry and sample only new steps.
    state = _driver_state(tmp_path)
    driver.advance_clubb_to_end(state, False, True, itime_end=2)
    assert state['l_stats']
    assert state['stats_writer']._time_index == 0
    assert state['_jax_stats'].names == state['stats_writer'].get_jax_layout().names
    assert all(not np.asarray(count).any() for count in state['_jax_stats'].nsamples)
    driver.advance_clubb_to_end(state, False, itime_start=3, itime_end=4)
    assert state['stats_writer']._time_index == 1
    np.testing.assert_array_equal(state['thlm'], np.full((2, 3), 4.))


@pytest.mark.parametrize('write_netcdf', [False, True])
def test_batched_silhs_source_output_guard(tmp_path, monkeypatch, write_netcdf):
    # Initialization only: the source rejects any enabled SILHS NetCDF output
    # in batch mode, including when optional sample fields are switched off.
    root = Path(initialization.__file__).resolve().parents[2]
    monkeypatch.setattr(error_code, '_debug_level', error_code._debug_level)
    path = create_case_namelist_file(
        'rico_silhs', tmp_path, multicol='4', batch_size=2,
        stats=str(root / 'input/stats/standard_stats.in'), debug='-1', max_iters=1,
    )
    text = prune_clubb_stats_namelist(path.read_text(), ['thlm'])
    text = set_stats_string(text, 'stats_output_filename', 'case_stats.nc' if write_netcdf else '')
    path.write_text(text)
    try:
        if write_netcdf:
            with pytest.raises(ValueError, match='Batch-mode stats NetCDF output'):
                initialization.init_clubb_case(str(path))
        else:
            state = initialization.init_clubb_case(str(path))
            try:
                assert state['ngrdcol'] == 2
                assert state['total_param_sets'] == 4
                assert state['stats_writer'].enabled
                assert state['stats_writer']._ncid is None
            finally:
                initialization.clean_up_clubb(state)
    finally:
        jax.clear_caches()


@pytest.mark.parametrize('name', ['C88', 'unknown(2)'])
def test_parameter_namelist_rejects_unknown_names(tmp_path, name):
    path = tmp_path / 'parameters.in'
    path.write_text(f'&clubb_params_nl\n {name}=0.25\n/\n')
    with pytest.raises(ValueError, match='Unknown CLUBB parameter'):
        init_clubb_params(2, str(path))
