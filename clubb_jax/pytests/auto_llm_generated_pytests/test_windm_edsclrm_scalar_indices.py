"""Explicit diffusion must respect scalar identities, independent of slot order."""
import numpy as np
import pytest
import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core import ErrInfo, NuVertResDep, SclrIdx
from clubb_jax.src.CLUBB_core.advance_windm_edsclrm_module import advance_windm_edsclrm
from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.model_flags import get_default_config_flags
from clubb_jax.src.CLUBB_core.numerical_check import check_clubb_settings
from clubb_jax.src.CLUBB_core.parameter_indices import nparams
from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params


def _args(thl_index, rt_index, enabled=True):
    gr = setup_grid(ngrdcol=2, deltaz=100., zm_init=0., zm_top=500., grid_type=1)
    nc, nt, nm = gr.ngrdcol, gr.nzt, gr.nzm
    zt = lambda v: jnp.full((nc, nt), v, dtype=jnp.float64)
    zm = lambda v: jnp.full((nc, nm), v, dtype=jnp.float64)
    return dict(
        nzm=nm, nzt=nt, ngrdcol=nc, edsclr_dim=4, gr=gr,
        sclr_idx=SclrIdx(-1, -1, -1, rt_index, thl_index, -1), dt=60.,
        wm_zt=zt(0), Kh_zm=zm(0), clubb_params=jnp.ones((nc, nparams)),
        ug=zt(1), vg=zt(0), um_ref=zt(1), vm_ref=zt(0),
        wp2=zm(.01), up2=zm(.01), vp2=zm(.01), um_forcing=zt(0), vm_forcing=zt(0),
        edsclrm_forcing=jnp.broadcast_to(jnp.array([.01, .02, .03, .04]), (nc, nt, 4)),
        p_in_Pa=jnp.broadcast_to(jnp.linspace(100000., 60000., nt), (nc, nt)),
        rho_ds_zm=zm(1), rho_ds_zt=zt(1), invrs_rho_ds_zt=zt(1),
        fcor=jnp.zeros(nc), l_implemented=True,
        nu_vert_res_dep=NuVertResDep(nm, *(jnp.zeros(nc) for _ in range(7))),
        ts_nudge=3600., tridiag_solve_method=2, l_predict_upwp_vpwp=False,
        l_upwind_xm_ma=True, l_uv_nudge=False, l_lmm_stepping=False,
        l_linearize_pbl_winds=False, l_do_expldiff_rtm_thlm=enabled, fill_holes_type=0,
        upwp_cl_num=0, vpwp_cl_num=0,
        stats=JaxStats.empty(l_sample=False, names=(), ncol=nc, max_nlev=nm),
        um=zt(1), vm=zt(0), thlm=zt(300), rtm=zt(.01),
        edsclrm=jnp.full((nc, nt, 4), 5.),
        upwp=zm(0), vpwp=zm(0), wpedsclrp=jnp.zeros((nc, nm, 4)),
        um_pert=zt(0), vm_pert=zt(0), upwp_pert=zm(0), vpwp_pert=zm(0),
        err_info=ErrInfo.initialized(nc),
    )


def _forcing_without_copies(args, thl_index, rt_index):
    forcing = np.array(args["edsclrm_forcing"])
    forcing[:, :, thl_index - 1] = 0.
    forcing[:, :, rt_index - 1] = 0.
    return forcing


@pytest.mark.parametrize("indices", [(1, 2), (2, 1), (1, 3), (4, 2)])
def test_named_slots_preserve_other_tracer_forcing(indices, capsys):
    args = _args(*indices)
    args["edsclrm_forcing"] = jnp.asarray(_forcing_without_copies(args, *indices))
    expected = np.array(args["edsclrm"] + args["dt"] * args["edsclrm_forcing"])
    expected[:, :, indices[0] - 1] = args["thlm"]
    expected[:, :, indices[1] - 1] = args["rtm"]
    out = advance_windm_edsclrm(**args)
    jax.block_until_ready(out)
    jax.effects_barrier()
    assert capsys.readouterr().err == ""
    assert not out[-2].is_fatal_host()
    np.testing.assert_allclose(out[4], 300., rtol=0, atol=1e-12)
    np.testing.assert_allclose(out[5], .01, rtol=0, atol=1e-14)
    np.testing.assert_allclose(out[6], expected, rtol=0, atol=1e-12)


@pytest.mark.parametrize("forced_copy", ["thl", "rt"])
def test_forcing_on_copies_is_fatal(forced_copy, capsys):
    """A forced copy prints its source diagnostic and flags only column two."""
    thl_index, rt_index = 1, 3
    args = _args(thl_index, rt_index)
    forcing = _forcing_without_copies(args, thl_index, rt_index)
    index = thl_index if forced_copy == "thl" else rt_index
    forcing[1, 0, index - 1] = 1.e-5  # second column, lowest level only
    args["edsclrm_forcing"] = jnp.asarray(forcing)
    err_info = advance_windm_edsclrm(**args)[-2]
    jax.block_until_ready(err_info)
    jax.effects_barrier()
    np.testing.assert_array_equal(np.asarray(err_info.fatal_mask()), [False, True])
    text = capsys.readouterr().err
    assert "requires zero edsclrm_forcing" in text
    assert "Fatal error in CLUBB - MPI Process" in text


@pytest.mark.parametrize("flux_copy", ["thl", "rt"])
def test_surface_flux_on_copies_is_fatal(flux_copy, capsys):
    """A copy's surface flux prints its source diagnostic and flags column two."""
    thl_index, rt_index = 1, 3
    args = _args(thl_index, rt_index)
    args["edsclrm_forcing"] = jnp.asarray(_forcing_without_copies(args, thl_index, rt_index))
    index = thl_index if flux_copy == "thl" else rt_index
    gr = args["gr"]
    wpedsclrp = np.zeros((gr.ngrdcol, gr.nzm, 4))
    wpedsclrp[:, gr.k_lb_zm, 1] = 1.e-3  # a surface flux on another tracer is fine
    wpedsclrp[1, gr.k_lb_zm, index - 1] = 1.e-5  # second column only
    args["wpedsclrp"] = jnp.asarray(wpedsclrp)
    err_info = advance_windm_edsclrm(**args)[-2]
    jax.block_until_ready(err_info)
    jax.effects_barrier()
    np.testing.assert_array_equal(np.asarray(err_info.fatal_mask()), [False, True])
    text = capsys.readouterr().err
    assert "requires zero wpedsclrp_sfc" in text
    assert "Fatal error in CLUBB - MPI Process" in text


def _check_settings(edsclr_dim, thl_index, rt_index, enabled=True):
    flags = get_default_config_flags()._replace(l_do_expldiff_rtm_thlm=enabled)
    return check_clubb_settings(
        err_info=ErrInfo.initialized(1), ngrdcol=1, edsclr_dim=edsclr_dim,
        sclr_idx=SclrIdx(-1, -1, -1, rt_index, thl_index, -1),
        params=init_clubb_params(1, ""), config_flags=flags,
        l_implemented=False, l_input_fields=False,
    )


@pytest.mark.parametrize("edsclr_dim,indices", [
    (4, (-1, 2)), (4, (1, 0)), (4, (1, 5)), (4, (2, 2)),
])
def test_settings_check_rejects_invalid_eddy_scalars(edsclr_dim, indices):
    assert _check_settings(edsclr_dim, *indices).is_fatal_host()


# With fewer than two eddy scalars the flag has no effect, so it is accepted.
@pytest.mark.parametrize("edsclr_dim,indices,enabled", [
    (2, (1, 2), True), (4, (4, 2), True), (0, (-1, -1), False),
    (0, (-1, -1), True), (1, (1, 1), True),
])
def test_settings_check_accepts_valid_or_disabled_eddy_scalars(edsclr_dim, indices, enabled):
    assert not _check_settings(edsclr_dim, *indices, enabled=enabled).is_fatal_host()


def test_disabled_option_keeps_all_forcing():
    args = _args(-1, -1, enabled=False)
    out = advance_windm_edsclrm(**args)
    assert not out[-2].is_fatal_host()
    np.testing.assert_allclose(
        out[6], args["edsclrm"] + args["dt"] * args["edsclrm_forcing"], rtol=0, atol=1e-12,
    )
