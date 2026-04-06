"""Test wrapper coverage for advance_xp2_xpyp call-tree routine."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.sclr_idx import SclrIdx


def _setup_env(tmp_path: Path, ngrdcol: int = 1):
    """Initialize grid and persistent UDT state used by advance_xp2_xpyp."""
    nzmax = 5
    dz = 100.0

    clubb_api.init_err_info(ngrdcol)
    flags = clubb_api.get_default_config_flags()
    clubb_api.init_config_flags(flags)

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
        err_info=ErrInfo(ngrdcol=ngrdcol),
    )

    clubb_params = clubb_api.init_clubb_params(ngrdcol, iunit=10, filename="")
    nu_vert_res_dep, _, _ = clubb_api.calc_derrived_params(
        gr=gr,
        ngrdcol=ngrdcol,
        grid_type=1,
        deltaz=np.full(ngrdcol, dz, dtype=np.float64),
        clubb_params=clubb_params,
        nu_vert_res_dep=None,
        l_prescribed_avg_deltaz=False,
    )
    pdf_implicit_coefs_terms = clubb_api.init_pdf_implicit(gr.nzm, ngrdcol, sclr_dim=1)
    sclr_idx = SclrIdx(0, 0, 0, 0, 0, 0)
    clubb_api.set_sclr_idx(*sclr_idx)
    err_info = ErrInfo(ngrdcol=ngrdcol)

    err_info = clubb_api.init_stats(
        registry_path="tests/test_stats_registry.in",
        output_path=str(tmp_path / "xp2_xpyp_stats.nc"),
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

    return gr, flags, clubb_params, sclr_idx, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info


def _make_args(gr, flags, clubb_params, sclr_idx, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info):
    """Construct stable, finite test inputs for advance_xp2_xpyp."""
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
        "invrs_tau_xp2_zm": full((ngrdcol, nzm), 1.0 / 300.0),
        "invrs_tau_c4_zm": full((ngrdcol, nzm), 1.0 / 300.0),
        "invrs_tau_c14_zm": full((ngrdcol, nzm), 1.0 / 300.0),
        "wm_zm": full((ngrdcol, nzm), 0.0),
        "rtm": full((ngrdcol, nzt), 0.01),
        "wprtp": full((ngrdcol, nzm), 0.0),
        "thlm": full((ngrdcol, nzt), 300.0),
        "wpthlp": full((ngrdcol, nzm), 0.0),
        "wpthvp": full((ngrdcol, nzm), 0.0),
        "um": full((ngrdcol, nzt), 1.0),
        "vm": full((ngrdcol, nzt), 0.0),
        "wp2": full((ngrdcol, nzm), 0.01),
        "wp2_zt": full((ngrdcol, nzt), 0.01),
        "wp3": full((ngrdcol, nzt), 0.0),
        "upwp": full((ngrdcol, nzm), 0.0),
        "vpwp": full((ngrdcol, nzm), 0.0),
        "sigma_sqd_w": full((ngrdcol, nzm), 0.5),
        "wprtp2": full((ngrdcol, nzt), 0.0),
        "wpthlp2": full((ngrdcol, nzt), 0.0),
        "wprtpthlp": full((ngrdcol, nzt), 0.0),
        "kh_zt": full((ngrdcol, nzt), 1.0),
        "rtp2_forcing": full((ngrdcol, nzm), 0.0),
        "thlp2_forcing": full((ngrdcol, nzm), 0.0),
        "rtpthlp_forcing": full((ngrdcol, nzm), 0.0),
        "rho_ds_zm": full((ngrdcol, nzm), 1.0),
        "rho_ds_zt": full((ngrdcol, nzt), 1.0),
        "invrs_rho_ds_zm": full((ngrdcol, nzm), 1.0),
        "thv_ds_zm": full((ngrdcol, nzm), 300.0),
        "cloud_frac": full((ngrdcol, nzt), 0.0),
        "wp3_on_wp2": full((ngrdcol, nzm), 0.0),
        "wp3_on_wp2_zt": full((ngrdcol, nzt), 0.0),
        "dt": 60.0,
        "fcor_y": full((ngrdcol,), 0.0),
        "sclrm": full((ngrdcol, nzt, 1), 0.0),
        "wpsclrp": full((ngrdcol, nzm, 1), 0.0),
        "wpsclrp2": full((ngrdcol, nzt, 1), 0.0),
        "wpsclrprtp": full((ngrdcol, nzt, 1), 0.0),
        "wpsclrpthlp": full((ngrdcol, nzt, 1), 0.0),
        "lhs_splat_wp2": full((ngrdcol, nzm), 0.0),
        "clubb_params": clubb_params,
        "iipdf_type": int(flags.iiPDF_type),
        "tridiag_solve_method": int(flags.tridiag_solve_method),
        "fill_holes_type": int(flags.fill_holes_type),
        "l_predict_upwp_vpwp": bool(flags.l_predict_upwp_vpwp),
        "l_ho_nontrad_coriolis": bool(flags.l_ho_nontrad_coriolis),
        "l_min_xp2_from_corr_wx": bool(flags.l_min_xp2_from_corr_wx),
        "l_c2_cloud_frac": bool(flags.l_C2_cloud_frac),
        "l_upwind_xpyp_ta": bool(flags.l_upwind_xpyp_ta),
        "l_godunov_upwind_xpyp_ta": bool(flags.l_godunov_upwind_xpyp_ta),
        "l_lmm_stepping": bool(flags.l_lmm_stepping),
        "rtp2": full((ngrdcol, nzm), 1.0e-6),
        "thlp2": full((ngrdcol, nzm), 0.01),
        "rtpthlp": full((ngrdcol, nzm), 0.0),
        "up2": full((ngrdcol, nzm), 0.01),
        "vp2": full((ngrdcol, nzm), 0.01),
        "sclrp2": full((ngrdcol, nzm, 1), 1.0e-6),
        "sclrprtp": full((ngrdcol, nzm, 1), 0.0),
        "sclrpthlp": full((ngrdcol, nzm, 1), 0.0),
        "sclr_idx": sclr_idx,
        "nu_vert_res_dep": nu_vert_res_dep,
        "pdf_implicit_coefs_terms": pdf_implicit_coefs_terms,
        "err_info": err_info,
    }
    return args


def test_advance_xp2_xpyp_returns_finite_arrays(tmp_path):
    """advance_xp2_xpyp should run and return finite output arrays."""
    gr, flags, clubb_params, sclr_idx, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info = _setup_env(tmp_path, ngrdcol=1)
    args = _make_args(gr, flags, clubb_params, sclr_idx, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info)

    try:
        out = clubb_api.advance_xp2_xpyp(**args)
    finally:
        err_info = clubb_api.finalize_stats(err_info=err_info)

    assert isinstance(out, tuple)
    assert len(out) == 9
    for arr in out[:-1]:
        assert np.all(np.isfinite(arr))
    assert isinstance(out[-1], ErrInfo)


def test_advance_xp2_xpyp_updates_match_return_values(tmp_path):
    """Returned arrays should match in-place updated input buffers."""
    gr, flags, clubb_params, sclr_idx, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info = _setup_env(tmp_path, ngrdcol=1)
    args = _make_args(gr, flags, clubb_params, sclr_idx, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info)

    rtp2_in = args["rtp2"]
    thlp2_in = args["thlp2"]
    rtpthlp_in = args["rtpthlp"]
    up2_in = args["up2"]
    vp2_in = args["vp2"]
    sclrp2_in = args["sclrp2"]
    sclrprtp_in = args["sclrprtp"]
    sclrpthlp_in = args["sclrpthlp"]

    try:
        (
            rtp2_out, thlp2_out, rtpthlp_out, up2_out, vp2_out,
            sclrp2_out, sclrprtp_out, sclrpthlp_out, err_info_out,
        ) = clubb_api.advance_xp2_xpyp(**args)
    finally:
        err_info = clubb_api.finalize_stats(err_info=err_info)

    np.testing.assert_allclose(rtp2_out, rtp2_in)
    np.testing.assert_allclose(thlp2_out, thlp2_in)
    np.testing.assert_allclose(rtpthlp_out, rtpthlp_in)
    np.testing.assert_allclose(up2_out, up2_in)
    np.testing.assert_allclose(vp2_out, vp2_in)
    np.testing.assert_allclose(sclrp2_out, sclrp2_in)
    np.testing.assert_allclose(sclrprtp_out, sclrprtp_in)
    np.testing.assert_allclose(sclrpthlp_out, sclrpthlp_in)
    assert isinstance(err_info_out, ErrInfo)


def test_update_xp2_mc_returns_finite_arrays(tmp_path):
    """update_xp2_mc should run and return finite zm-grid tendency arrays."""
    gr, flags, clubb_params, sclr_idx, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info = _setup_env(tmp_path, ngrdcol=1)
    ngrdcol, nzt, nzm = gr.ngrdcol, gr.nzt, gr.nzm
    pdf_params = clubb_api.init_pdf_params(nzt, ngrdcol)

    def full(shape, val):
        return np.full(shape, val, dtype=np.float64, order="F")

    try:
        out = clubb_api.update_xp2_mc(
            gr=gr,
            nzm=nzm,
            nzt=nzt,
            ngrdcol=ngrdcol,
            dt=60.0,
            cloud_frac=full((ngrdcol, nzt), 0.2),
            rcm=full((ngrdcol, nzt), 1.0e-4),
            rvm=full((ngrdcol, nzt), 0.01),
            thlm=full((ngrdcol, nzt), 300.0),
            wm=full((ngrdcol, nzt), 0.0),
            exner=full((ngrdcol, nzt), 1.0),
            rrm_evap=full((ngrdcol, nzt), -1.0e-7),
            pdf_params=pdf_params,
            rtp2_mc=full((ngrdcol, nzm), 0.0),
            thlp2_mc=full((ngrdcol, nzm), 0.0),
            wprtp_mc=full((ngrdcol, nzm), 0.0),
            wpthlp_mc=full((ngrdcol, nzm), 0.0),
            rtpthlp_mc=full((ngrdcol, nzm), 0.0),
        )
    finally:
        err_info = clubb_api.finalize_stats(err_info=err_info)

    assert isinstance(out, tuple)
    assert len(out) == 5
    for arr in out:
        assert arr.shape == (ngrdcol, nzm)
        assert np.all(np.isfinite(arr))
