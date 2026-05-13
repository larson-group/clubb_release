"""Test wrapper coverage for advance_xp3 call-tree routine."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def _setup_grid_and_stats(tmp_path: Path, ngrdcol: int = 1):
    """Initialize grid, config, and stats state used by advance_xp3."""
    nzmax = 5
    dz = 100.0

    clubb_api.init_err_info(ngrdcol)
    clubb_api.init_config_flags(clubb_api.get_default_config_flags())
    err_info = ErrInfo(ngrdcol=ngrdcol)

    zm = np.arange(nzmax, dtype=np.float64) * dz
    zt = 0.5 * (zm[:-1] + zm[1:])

    gr, _err_info = clubb_api.setup_grid(
        nzmax=nzmax,
        ngrdcol=ngrdcol,
        sfc_elevation=np.zeros(ngrdcol, dtype=np.float64),
        l_implemented=True,
        l_ascending_grid=True,
        grid_type=1,
        deltaz=np.full(ngrdcol, dz, dtype=np.float64),
        zm_init=np.zeros(ngrdcol, dtype=np.float64),
        zm_top=np.full(ngrdcol, zm[-1], dtype=np.float64),
        momentum_heights=np.tile(zm, (ngrdcol, 1)),
        thermodynamic_heights=np.tile(zt, (ngrdcol, 1)),
        err_info=err_info,
    )

    clubb_params = clubb_api.init_clubb_params(ngrdcol, iunit=10, filename="")
    clubb_api.init_stats(
        registry_path="tests/test_stats_registry.in",
        output_path=str(tmp_path / "xp3_stats.nc"),
        ncol=ngrdcol,
        stats_tsamp=60.0,
        stats_tout=60.0,
        dt_main=60.0,
        day_in=1,
        month_in=1,
        year_in=2000,
        time_initial=0.0,
        nzt=gr.nzt,
        zt=np.asfortranarray(gr.zt),
        nzm=gr.nzm,
        zm=np.asfortranarray(gr.zm),
        err_info=err_info,
        sclr_dim=1,
        edsclr_dim=0,
        clubb_params=clubb_params,
        param_names=clubb_api.get_param_names(),
    )

    return gr, clubb_params, err_info


def _make_xp3_args(gr):
    """Construct physically reasonable inputs for advance_xp3."""
    ngrdcol, nzt, nzm = gr.ngrdcol, gr.nzt, gr.nzm

    def full(shape, val):
        return np.full(shape, val, dtype=np.float64, order="F")

    args = {
        "gr": gr,
        "nzm": nzm,
        "nzt": nzt,
        "ngrdcol": ngrdcol,
        "sclr_dim": 1,
        "sclr_tol": np.array([1.0e-8], dtype=np.float64, order="F"),
        "dt": 60.0,
        "rtm": full((ngrdcol, nzt), 0.01),
        "thlm": full((ngrdcol, nzt), 300.0),
        "rtp2": full((ngrdcol, nzm), 1.0e-6),
        "thlp2": full((ngrdcol, nzm), 0.01),
        "wprtp": full((ngrdcol, nzm), 1.0e-4),
        "wpthlp": full((ngrdcol, nzm), 1.0e-3),
        "wprtp2": full((ngrdcol, nzt), 1.0e-6),
        "wpthlp2": full((ngrdcol, nzt), 1.0e-4),
        "rho_ds_zm": full((ngrdcol, nzm), 1.0),
        "invrs_rho_ds_zt": full((ngrdcol, nzt), 1.0),
        "invrs_tau_zt": full((ngrdcol, nzt), 1.0 / 300.0),
        "tau_max_zt": full((ngrdcol, nzt), 300.0),
        "sclrm": full((ngrdcol, nzt, 1), 0.0),
        "sclrp2": full((ngrdcol, nzm, 1), 1.0e-6),
        "wpsclrp": full((ngrdcol, nzm, 1), 0.0),
        "wpsclrp2": full((ngrdcol, nzt, 1), 0.0),
        "wp2": full((ngrdcol, nzm), 0.1),
        "wp3": full((ngrdcol, nzt), 0.0),
        "upwp": full((ngrdcol, nzm), 0.0),
        "vpwp": full((ngrdcol, nzm), 0.0),
        "up2": full((ngrdcol, nzm), 0.1),
        "vp2": full((ngrdcol, nzm), 0.1),
        "thvm": full((ngrdcol, nzt), 300.0),
        "sigma_sqd_w": full((ngrdcol, nzm), 0.1),
        "clubb_params": clubb_api.init_clubb_params(ngrdcol, iunit=10, filename=""),
        "l_lmm_stepping": False,
        "rtp3": full((ngrdcol, nzt), 0.0),
        "thlp3": full((ngrdcol, nzt), 0.0),
        "sclrp3": full((ngrdcol, nzt, 1), 0.0),
        "up3": full((ngrdcol, nzt), 0.0),
        "vp3": full((ngrdcol, nzt), 0.0),
    }
    return args


def test_advance_xp3_returns_finite_arrays(tmp_path):
    """advance_xp3 should run and return finite outputs."""
    gr, clubb_params, err_info = _setup_grid_and_stats(tmp_path, ngrdcol=1)
    args = _make_xp3_args(gr)
    args["clubb_params"] = clubb_params

    try:
        out = clubb_api.advance_xp3(**args)
    finally:
        clubb_api.finalize_stats(err_info=err_info)

    assert isinstance(out, tuple)
    assert len(out) == 5
    for arr in out:
        assert np.all(np.isfinite(arr))


def test_advance_xp3_updates_match_return_values(tmp_path):
    """Returned arrays should match the in-place updated input buffers."""
    gr, clubb_params, err_info = _setup_grid_and_stats(tmp_path, ngrdcol=1)
    args = _make_xp3_args(gr)
    args["clubb_params"] = clubb_params

    rtp3_in = args["rtp3"]
    thlp3_in = args["thlp3"]
    sclrp3_in = args["sclrp3"]
    up3_in = args["up3"]
    vp3_in = args["vp3"]

    try:
        rtp3_out, thlp3_out, sclrp3_out, up3_out, vp3_out = clubb_api.advance_xp3(**args)
    finally:
        clubb_api.finalize_stats(err_info=err_info)

    np.testing.assert_allclose(rtp3_out, rtp3_in)
    np.testing.assert_allclose(thlp3_out, thlp3_in)
    np.testing.assert_allclose(sclrp3_out, sclrp3_in)
    np.testing.assert_allclose(up3_out, up3_in)
    np.testing.assert_allclose(vp3_out, vp3_in)


def test_compute_xp3_returns_finite_arrays(tmp_path):
    """compute_xp3 should run and return finite diagnosed outputs."""
    gr, clubb_params, err_info = _setup_grid_and_stats(tmp_path, ngrdcol=1)
    args = _make_xp3_args(gr)
    args["clubb_params"] = clubb_params

    try:
        out = clubb_api.compute_xp3(
            gr=args["gr"], nzm=args["nzm"], nzt=args["nzt"], ngrdcol=args["ngrdcol"],
            sclr_dim=args["sclr_dim"], sclr_tol=args["sclr_tol"], iipdf_type=1,
            clubb_params=args["clubb_params"], wp2=args["wp2"], wp3=args["wp3"],
            thvm=args["thvm"], wprtp=args["wprtp"], wpthlp=args["wpthlp"],
            rtp2=args["rtp2"], thlp2=args["thlp2"], upwp=args["upwp"], vpwp=args["vpwp"],
            up2=args["up2"], vp2=args["vp2"], sigma_sqd_w=args["sigma_sqd_w"],
            wpsclrp=args["wpsclrp"], sclrp2=args["sclrp2"], rtp3=args["rtp3"],
            thlp3=args["thlp3"], up3=args["up3"], vp3=args["vp3"], sclrp3=args["sclrp3"],
        )
    finally:
        clubb_api.finalize_stats(err_info=err_info)

    assert isinstance(out, tuple)
    assert len(out) == 5
    for arr in out:
        assert np.all(np.isfinite(arr))
