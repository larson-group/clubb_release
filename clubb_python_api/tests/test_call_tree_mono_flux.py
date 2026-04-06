"""Test wrapper for mono_flux_limiter::calc_turb_adv_range."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def _make_grid(ngrdcol=1, nzmax=11, dz=300.0):
    """Create a standard test grid."""
    clubb_api.init_err_info(ngrdcol)
    flags = clubb_api.get_default_config_flags()
    clubb_api.init_config_flags(flags)

    zm = np.arange(nzmax, dtype=np.float64) * dz
    zt = 0.5 * (zm[:-1] + zm[1:])

    gr, _err_info = clubb_api.setup_grid(
        nzmax=nzmax, ngrdcol=ngrdcol,
        sfc_elevation=np.zeros(ngrdcol, dtype=np.float64),
        l_implemented=True, l_ascending_grid=True,
        grid_type=1,
        deltaz=np.full(ngrdcol, dz, dtype=np.float64),
        zm_init=np.zeros(ngrdcol, dtype=np.float64),
        zm_top=np.full(ngrdcol, zm[-1], dtype=np.float64),
        momentum_heights=np.tile(zm, (ngrdcol, 1)),
        thermodynamic_heights=np.tile(zt, (ngrdcol, 1)),
        err_info=ErrInfo(ngrdcol=ngrdcol),
    )
    return gr


def test_calc_turb_adv_range_basic_bounds():
    """calc_turb_adv_range should return in-domain low/high bounds for each level."""
    gr = _make_grid()

    # Mildly asymmetric two-component vertical velocity distribution.
    w_1 = np.full((gr.ngrdcol, gr.nzm), 0.4, dtype=np.float64)
    w_2 = np.full((gr.ngrdcol, gr.nzm), -0.3, dtype=np.float64)
    var_1 = np.full((gr.ngrdcol, gr.nzm), 0.2, dtype=np.float64)
    var_2 = np.full((gr.ngrdcol, gr.nzm), 0.2, dtype=np.float64)
    mixt = np.full((gr.ngrdcol, gr.nzm), 0.6, dtype=np.float64)

    low, high = clubb_api.calc_turb_adv_range(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol, dt=60.0,
        w_1_zm=w_1, w_2_zm=w_2,
        varnce_w_1_zm=var_1, varnce_w_2_zm=var_2,
        mixt_frac_zm=mixt,
    )

    assert low.shape == (gr.ngrdcol, gr.nzt)
    assert high.shape == (gr.ngrdcol, gr.nzt)

    for i in range(gr.ngrdcol):
        for k0 in range(gr.nzt):
            assert 0 <= low[i, k0] < gr.nzt
            assert 0 <= high[i, k0] < gr.nzt
            assert low[i, k0] <= k0 <= high[i, k0]


def test_monotonic_turbulent_flux_limit_preserves_zero_flux_state(tmp_path):
    """With zero flux and forcing, the monotonic limiter should leave xm and wpxp unchanged."""
    gr = _make_grid()
    ngrdcol, nzm, nzt = gr.ngrdcol, gr.nzm, gr.nzt
    err_info = ErrInfo(ngrdcol=ngrdcol)
    flags = clubb_api.get_default_config_flags()

    registry = str(Path(__file__).resolve().parent / "test_stats_registry.in")
    clubb_params = clubb_api.init_clubb_params(ngrdcol, iunit=10, filename="")

    clubb_api.init_stats(
        registry_path=registry,
        output_path=str(tmp_path / "mono_flux_stats.nc"),
        ncol=ngrdcol,
        stats_tsamp=60.0,
        stats_tout=60.0,
        dt_main=60.0,
        day_in=1,
        month_in=1,
        year_in=2000,
        time_initial=0.0,
        nzt=nzt,
        zt=np.asfortranarray(gr.zt),
        nzm=nzm,
        zm=np.asfortranarray(gr.zm),
        err_info=err_info,
        sclr_dim=0,
        edsclr_dim=0,
        clubb_params=clubb_params,
        param_names=clubb_api.get_param_names(),
    )

    try:
        xm = np.full((ngrdcol, nzt), 300.0, dtype=np.float64, order="F")
        wpxp = np.zeros((ngrdcol, nzm), dtype=np.float64, order="F")
        low = np.tile(np.arange(nzt, dtype=np.int32), (ngrdcol, 1))
        high = np.tile(np.arange(nzt, dtype=np.int32), (ngrdcol, 1))

        xm_out, wpxp_out = clubb_api.monotonic_turbulent_flux_limit(
            gr=gr,
            nzm=nzm,
            nzt=nzt,
            ngrdcol=ngrdcol,
            solve_type=1,
            dt=60.0,
            xm_old=xm,
            xp2=np.full((ngrdcol, nzm), 1.0e-6, dtype=np.float64, order="F"),
            wm_zt=np.zeros((ngrdcol, nzt), dtype=np.float64, order="F"),
            xm_forcing=np.zeros((ngrdcol, nzt), dtype=np.float64, order="F"),
            rho_ds_zm=np.ones((ngrdcol, nzm), dtype=np.float64, order="F"),
            rho_ds_zt=np.ones((ngrdcol, nzt), dtype=np.float64, order="F"),
            invrs_rho_ds_zm=np.ones((ngrdcol, nzm), dtype=np.float64, order="F"),
            invrs_rho_ds_zt=np.ones((ngrdcol, nzt), dtype=np.float64, order="F"),
            xp2_threshold=1.0e-12,
            xm_tol=1.0e-12,
            l_implemented=True,
            low_lev_effect=low,
            high_lev_effect=high,
            tridiag_solve_method=flags.tridiag_solve_method,
            l_upwind_xm_ma=flags.l_upwind_xm_ma,
            l_mono_flux_lim_spikefix=flags.l_mono_flux_lim_spikefix,
            xm=xm.copy(order="F"),
            wpxp=wpxp.copy(order="F"),
        )
    finally:
        clubb_api.finalize_stats(err_info=err_info)

    np.testing.assert_allclose(xm_out, xm)
    np.testing.assert_allclose(wpxp_out, wpxp)
