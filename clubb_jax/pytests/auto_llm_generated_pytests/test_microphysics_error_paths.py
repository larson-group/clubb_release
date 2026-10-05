"""Source fatal-return boundaries, including whole-batch error handling.

Tests with lightweight argument doubles call the undecorated routine body;
the outer jax.jit(run) still checks each error path under compilation. Real
standalone integrations exercise the decorated entry points with model pytrees.
"""
from types import SimpleNamespace
import inspect
import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.Microphys import advance_microphys_module as transport
from clubb_jax.src.Microphys import parameters_microphys as parameters
from clubb_jax.pytests.microphysics_test_inputs import initialize


def test_namelist_cannot_override_enumerations():
    initialize(lh_microphys_disabled=1, morrison_lognormal=99)
    assert parameters.lh_microphys_disabled == 3
    assert parameters.morrison_lognormal == 2
    with pytest.raises(ValueError, match="SILHS"):
        initialize(
            lh_microphys_disabled=1, lh_microphys_type="interactive", microphys_scheme="none"
        )


@pytest.mark.parametrize('ncol', [1, 3])
@pytest.mark.parametrize('compiled', [False, True])
def test_failed_solve_preserves_statistics(monkeypatch, ncol, compiled):
    from clubb_jax.src.CLUBB_core import matrix_solver_wrapper
    gr = SimpleNamespace(ngrdcol=ncol, nzt=4)
    stats = JaxStats.empty(l_sample=True, names=('rrm_ma', 'rrm_sd', 'rrm_ta'),
                           ncol=ncol, max_nlev=4)
    stats = stats.begin_budget('rrm_ta', jnp.ones((ncol, 4)))
    def failed_solver(solve_type, method, ncol, nzt, lhs, rhs, err):
        return err.set_fatal(jnp.arange(ncol) == ncol-1), -jnp.ones_like(rhs), None
    monkeypatch.setattr(matrix_solver_wrapper, 'tridiag_solve', failed_solver)
    bands = jnp.ones((3, ncol, 4))
    def run(rhs):
        return transport.microphys_solve(gr, ncol, 'rrm', True, bands, bands,
            bands, bands, jnp.ones_like(rhs), 2, stats, bands, rhs,
            jnp.ones_like(rhs), ErrInfo.initialized(ncol))
    result = (jax.jit(run) if compiled else run)(jnp.ones((ncol, 4)))
    assert result[-1].is_fatal()
    np.testing.assert_array_equal(result[3], -np.ones((ncol, 4)))
    for before, after in zip(jax.tree_util.tree_leaves(stats),
                             jax.tree_util.tree_leaves(result[0])):
        np.testing.assert_array_equal(before, after)


@pytest.mark.parametrize('compiled', [False, True])
def test_failed_cloud_number_solve_skips_clipping_and_final_stats(monkeypatch, compiled):
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    initialize(microphys_scheme='morrison', l_predict_Nc=True,
               specify_aerosol='morrison_no_aerosol', l_in_cloud_Nc_diff=False)
    gr = setup_grid(2, 100., 0., 600.)
    one = jnp.ones((2, gr.nzt)); zero = jnp.zeros_like(one)
    stats = JaxStats.empty(l_sample=True, names=('Ncm_cl', 'wpNcp'),
        grids=('zt', 'zm'), ncol=2, max_nlev=gr.nzm,
        grid_nlev=(gr.nzt, gr.nzm, 1, gr.nzt, 1, gr.nzt, gr.nzm))
    def failed_solve(gr, ncol, solve_type, l_sed, lhs_ta, lhs_ma,
                     sed_turb_lhs, sed_diff_lhs, cloud_frac, method,
                     stats, lhs, rhs, hmm, err):
        return stats, lhs, rhs, -jnp.ones_like(hmm), err.set_fatal(jnp.array([False, True]))
    monkeypatch.setattr(transport, 'microphys_solve', failed_solve)
    def run(source):
        return transport.advance_Ncm(gr, 2, 10., zero, .5*one,
            jnp.zeros((2, gr.nzm)), zero, jnp.ones((2, gr.nzm)), one, one,
            source, SimpleNamespace(nu_hm=jnp.zeros(2)), False, 2, stats,
            5.e7*one, 1.e8*one, ErrInfo.initialized(2))
    result = (jax.jit(run) if compiled else run)(zero)
    assert result[3].is_fatal()
    np.testing.assert_array_equal(result[1], -one)
    np.testing.assert_array_equal(result[2], -2*one)
    for counts in result[0].nsamples:
        assert np.all(np.asarray(counts) == 0)
    np.testing.assert_array_equal(result[4], 0.)


def test_fatal_host_dump_emits_fields_in_source_order(capsys):
    kwargs = {name: jnp.array([1.]) for name in inspect.signature(transport.write_adv_micro_errors).parameters}
    kwargs.update(gr=None, ngrdcol=1, hydromet_dim=2,
                  nu_vert_res_dep=SimpleNamespace(nu_hm=jnp.array([0.])),
                  err_info=ErrInfo.initialized(1))
    transport.write_adv_micro_errors(**kwargs)
    assert capsys.readouterr().err == ''
    kwargs['err_info'] = kwargs['err_info'].set_fatal()
    transport.write_adv_micro_errors(**kwargs)
    output = capsys.readouterr().err
    assert 'Error in advance_microphys' in output
    assert output.index('Intent(in)') < output.index('Intent(inout)') < output.index('Intent(out)')
    assert output.index('hydromet_mc =') < output.index('hydromet =') < output.index('wpNcp =')


@pytest.mark.parametrize("compiled", [False, True])
def test_pdf_fatal_skips_mixed_moments_and_statistics(monkeypatch, compiled):
    from clubb_jax.src.Microphys import pdf_hydromet_microphys_wrapper as wrapper

    initialize(microphys_scheme="khairoutdinov_kogan")
    gr = SimpleNamespace(ngrdcol=2, nzt=3, nzm=4)
    stats = JaxStats.empty(l_sample=True, names=("rtprrp",), ncol=2, max_nlev=3)

    def setup(*args):
        hydromet_pdf_params = args[-2]
        err = args[-4].set_fatal(jnp.array([False, True]))
        means = (jnp.ones((2, 3, 4)),) * 4
        matrices = (jnp.ones((2, 3, 4, 4)),) * 4
        return (
            err,
            jnp.ones((2, 4, 2)),
            *means,
            *matrices,
            args[-3],
            hydromet_pdf_params,
            args[-1],
        )

    def mixed(*args):
        values = jnp.ones((2, 3, 2))
        return values, values, values, args[-1].update("rtprrp", values[..., 0])

    monkeypatch.setattr(wrapper, "setup_pdf_parameters_api", setup)
    monkeypatch.setattr(wrapper, "hydrometeor_mixed_moments", mixed)
    flags = SimpleNamespace(
        iiPDF_type=1,
        l_use_precip_frac=True,
        l_diagnose_correlations=False,
        l_calc_w_corr=False,
        l_const_Nc_in_cloud=True,
        l_fix_w_chi_eta_correlations=True,
    )

    def run(hydromet):
        return wrapper.pdf_hydromet_microphys_prep.__wrapped__(
            gr, 2, 4, 2,                       # In
            0, 0.0,                            # In
            None, None, None,                  # In
            None, None, None, hydromet, None,  # In
            None, None,                        # In
            None, None, jnp.ones((2, 4)),      # In
            flags, None,                       # In
            False,                             # In
            stats,                             # InOut
            ErrInfo.initialized(2),            # InOut
            None,                              # InOut
            None,                              # InOut
        )

    result = (jax.jit(run) if compiled else run)(jnp.ones((2, 3, 2)))
    assert result[1].is_fatal()
    for value in result[12:15]:
        np.testing.assert_array_equal(value, 0.0)
    assert np.all(np.asarray(result[0].nsamples[0]) == 0)


@pytest.mark.parametrize('compiled', [False, True])
def test_failed_hydrometeor_transport_does_not_advance_cloud_number(monkeypatch, compiled):
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    _, _, metadata, *_ = initialize(microphys_scheme='khairoutdinov_kogan')
    gr = setup_grid(2, 100., 0., 600.)
    one = jnp.ones((2, gr.nzt)); zero = jnp.zeros_like(one)
    hydromet = jnp.ones((2, gr.nzt, 2))
    flux = jnp.zeros((2, gr.nzm, 2))
    stats = JaxStats.empty(l_sample=True, names=('Ncm', 'Nc_in_cloud'), ncol=2, max_nlev=gr.nzt)
    monkeypatch.setattr(transport, 'calculate_K_hm', lambda *args: flux)
    def failed_hydro(*args):
        # Preserve the strict production return contract; only transport fails.
        return (args[19], 2*hydromet, hydromet, flux, zero, zero,
                args[-1].set_fatal(jnp.array([False, True])), flux, flux, flux, hydromet)
    monkeypatch.setattr(transport, 'advance_hydrometeor', failed_hydro)
    def run(ncm):
        return transport.advance_microphys.__wrapped__(
            gr, 2, 10., 10., 2, metadata,
            zero, jnp.ones((2, gr.nzm)), one, one, jnp.ones((2, gr.nzm)),
            zero, .5*one, jnp.ones((2, gr.nzm)), jnp.zeros((2, gr.nzm)),
            jnp.ones((2, gr.nzm)), one, one, hydromet, zero, jnp.ones((2, gr.nzm)),
            hydromet, hydromet, None, SimpleNamespace(nu_hm=jnp.zeros(2)),
            2, 1, False, stats, hydromet, hydromet, flux, flux, ncm,
            1.e8*one, zero, zero, ErrInfo.initialized(2),
        )
    result = (jax.jit(run) if compiled else run)(3.e7*one)
    assert result[9].is_fatal()
    np.testing.assert_array_equal(result[1], 2*hydromet)
    np.testing.assert_array_equal(result[5], 3.e7*one)
    assert np.all(np.asarray(result[0].nsamples[0]) == 0)


@pytest.mark.parametrize("compiled", [False, True])
def test_silhs_fatal_skips_clipping_and_statistics(monkeypatch, compiled):
    from clubb_jax.src.Microphys import pdf_hydromet_microphys_wrapper as wrapper
    from clubb_jax.src.SILHS import silhs_api_module as api
    from clubb_jax.src.SILHS import latin_hypercube_driver_module as driver
    from clubb_jax.src.SILHS.latin_hypercube_arrays import LatinHypercubeArrays

    initialize(
        microphys_scheme="khairoutdinov_kogan", lh_microphys_type="interactive", l_local_kk=True
    )
    gr = SimpleNamespace(ngrdcol=2, nzt=3, nzm=4, dzt=jnp.ones((2, 3)))
    ns = parameters.lh_num_samples
    stats = JaxStats.empty(l_sample=True, names=("lh_rtm",), ncol=2, max_nlev=3)
    sampling_state = LatinHypercubeArrays(jnp.zeros((ns, 6), dtype=jnp.int32), jnp.int32(0))

    def setup(*args):
        return (
            args[-4],
            jnp.ones((2, 4, 2)),
            *(jnp.ones((2, 3, 4)),) * 4,
            *(jnp.ones((2, 3, 4, 4)),) * 4,
            args[-3],
            args[-2],
            args[-1],
        )

    def mixed(*args):
        return *(jnp.zeros((2, 3, 2)),) * 3, args[-1]

    def generate(*args):
        return (
            args[-3].set_fatal(jnp.array([False, True])),
            jnp.ones((2, ns, 3, 4)),
            jnp.zeros((2, ns, 3), dtype=jnp.int32),
            jnp.ones((2, ns, 3)),
            args[-2],
            args[-1],
        )

    def clip(*args):
        return 2 * args[7], *(jnp.ones((2, ns, 3)),) * 5

    def accumulate(*args):
        return args[-1].update("lh_rtm", jnp.ones((2, 3)))

    monkeypatch.setattr(wrapper, "setup_pdf_parameters_api", setup)
    monkeypatch.setattr(wrapper, "hydrometeor_mixed_moments", mixed)
    monkeypatch.setattr(api, "generate_silhs_sample_api", generate)
    monkeypatch.setattr(api, "clip_transform_silhs_output_api", clip)
    monkeypatch.setattr(driver, "stats_accumulate_lh_api", accumulate)
    flags = SimpleNamespace(
        iiPDF_type=1,
        l_use_precip_frac=True,
        l_diagnose_correlations=False,
        l_calc_w_corr=False,
        l_const_Nc_in_cloud=True,
        l_fix_w_chi_eta_correlations=True,
    )

    def run(hydromet):
        return wrapper.pdf_hydromet_microphys_prep.__wrapped__(
            gr, 2, 4, 2,                       # In
            1, 0.0,                            # In
            None, None, None,                  # In
            None, None, None, hydromet, None,  # In
            None, None,                        # In
            None, None, jnp.ones((2, 4)),      # In
            flags, None,                       # In
            False,                             # In
            stats,                             # InOut
            ErrInfo.initialized(2),            # InOut
            None,                              # InOut
            sampling_state,                    # InOut
        )

    result = (jax.jit(run) if compiled else run)(jnp.ones((2, 3, 2)))
    assert result[1].is_fatal()
    np.testing.assert_array_equal(result[15], 1.0)
    for value in result[18:23]:
        np.testing.assert_array_equal(value, 0.0)
    np.testing.assert_array_equal(result[0].nsamples[0], 0)
