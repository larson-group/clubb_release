"""Test wrapper coverage for advance_wp2_wp3 call-tree routine."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def _setup_env(tmp_path: Path, ngrdcol: int = 1):
    """Initialize grid and persistent UDT state used by advance_wp2_wp3."""
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
    pdf_implicit_coefs_terms = clubb_api.init_pdf_implicit(gr.nzm, ngrdcol, sclr_dim=0)
    err_info = ErrInfo(ngrdcol=ngrdcol)

    err_info = clubb_api.init_stats(
        registry_path="tests/test_stats_registry.in",
        output_path=str(tmp_path / "wp2_wp3_stats.nc"),
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
        sclr_dim=0,
        edsclr_dim=0,
        clubb_params=clubb_params,
        param_names=clubb_api.get_param_names(),
    )

    return gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info


def _make_args(gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info):
    """Construct stable, finite test inputs for advance_wp2_wp3."""
    ngrdcol, nzm, nzt = gr.ngrdcol, gr.nzm, gr.nzt

    def full(shape, val):
        return np.full(shape, val, dtype=np.float64, order="F")

    args = {
        "gr": gr,
        "nzm": nzm,
        "nzt": nzt,
        "ngrdcol": ngrdcol,
        "dt": 60.0,
        "sfc_elevation": full((ngrdcol,), 0.0),
        "fcor_y": full((ngrdcol,), 0.0),
        "sigma_sqd_w": full((ngrdcol, nzm), 0.5),
        "wm_zm": full((ngrdcol, nzm), 0.0),
        "wm_zt": full((ngrdcol, nzt), 0.0),
        "a3_coef": full((ngrdcol, nzm), 0.1),
        "a3_coef_zt": full((ngrdcol, nzt), 0.1),
        "wp3_on_wp2": full((ngrdcol, nzm), 0.0),
        "wpup2": full((ngrdcol, nzt), 0.0),
        "wpvp2": full((ngrdcol, nzt), 0.0),
        "wp2up2": full((ngrdcol, nzm), 0.0),
        "wp2vp2": full((ngrdcol, nzm), 0.0),
        "wp4": full((ngrdcol, nzm), 1.0e-4),
        "wpthvp": full((ngrdcol, nzm), 0.0),
        "wp2thvp": full((ngrdcol, nzt), 0.0),
        "wp2up": full((ngrdcol, nzt), 0.0),
        "um": full((ngrdcol, nzt), 1.0),
        "vm": full((ngrdcol, nzt), 0.0),
        "upwp": full((ngrdcol, nzm), 0.0),
        "vpwp": full((ngrdcol, nzm), 0.0),
        "em": full((ngrdcol, nzm), 0.1),
        "kh_zm": full((ngrdcol, nzm), 1.0),
        "kh_zt": full((ngrdcol, nzt), 1.0),
        "invrs_tau_c4_zm": full((ngrdcol, nzm), 1.0 / 300.0),
        "invrs_tau_wp3_zt": full((ngrdcol, nzt), 1.0 / 300.0),
        "invrs_tau_c1_zm": full((ngrdcol, nzm), 1.0 / 300.0),
        "skw_zm": full((ngrdcol, nzm), 0.0),
        "skw_zt": full((ngrdcol, nzt), 0.0),
        "rho_ds_zm": full((ngrdcol, nzm), 1.0),
        "rho_ds_zt": full((ngrdcol, nzt), 1.0),
        "invrs_rho_ds_zm": full((ngrdcol, nzm), 1.0),
        "invrs_rho_ds_zt": full((ngrdcol, nzt), 1.0),
        "thv_ds_zm": full((ngrdcol, nzm), 300.0),
        "thv_ds_zt": full((ngrdcol, nzt), 300.0),
        "mixt_frac": full((ngrdcol, nzt), 0.5),
        "cx_fnc_richardson": full((ngrdcol, nzm), 1.0),
        "lhs_splat_wp2": full((ngrdcol, nzm), 0.0),
        "lhs_splat_wp3": full((ngrdcol, nzt), 0.0),
        "wprtp": full((ngrdcol, nzm), 0.0),
        "wpthlp": full((ngrdcol, nzm), 0.0),
        "rtp2": full((ngrdcol, nzm), 1.0e-6),
        "thlp2": full((ngrdcol, nzm), 0.01),
        "clubb_params": clubb_params,
        "iipdf_type": int(flags.iiPDF_type),
        "penta_solve_method": int(flags.penta_solve_method),
        "fill_holes_type": int(flags.fill_holes_type),
        "l_min_wp2_from_corr_wx": bool(flags.l_min_wp2_from_corr_wx),
        "l_upwind_xm_ma": bool(flags.l_upwind_xm_ma),
        "l_tke_aniso": bool(flags.l_tke_aniso),
        "l_standard_term_ta": bool(flags.l_standard_term_ta),
        "l_partial_upwind_wp3": bool(flags.l_partial_upwind_wp3),
        "l_damp_wp2_using_em": bool(flags.l_damp_wp2_using_em),
        "l_use_c11_richardson": bool(flags.l_use_C11_Richardson),
        "l_damp_wp3_skw_squared": bool(flags.l_damp_wp3_Skw_squared),
        "l_lmm_stepping": bool(flags.l_lmm_stepping),
        "l_use_tke_in_wp3_pr_turb_term": bool(flags.l_use_tke_in_wp3_pr_turb_term),
        "l_use_tke_in_wp2_wp3_k_dfsn": False,
        "l_use_wp3_lim_with_smth_heaviside": bool(flags.l_use_wp3_lim_with_smth_Heaviside),
        "l_wp2_fill_holes_tke": bool(flags.l_wp2_fill_holes_tke),
        "l_ho_nontrad_coriolis": bool(flags.l_ho_nontrad_coriolis),
        "up2": full((ngrdcol, nzm), 0.01),
        "vp2": full((ngrdcol, nzm), 0.01),
        "wp2": full((ngrdcol, nzm), 0.01),
        "wp3": full((ngrdcol, nzt), 0.0),
        "wp3_zm": full((ngrdcol, nzm), 0.0),
        "wp2_zt": full((ngrdcol, nzt), 0.01),
        "nu_vert_res_dep": nu_vert_res_dep,
        "pdf_implicit_coefs_terms": pdf_implicit_coefs_terms,
        "err_info": err_info,
    }
    return args


def test_advance_wp2_wp3_returns_finite_arrays(tmp_path):
    """advance_wp2_wp3 should run and return finite output arrays."""
    gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info = _setup_env(tmp_path, ngrdcol=1)
    args = _make_args(gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info)

    try:
        out = clubb_api.advance_wp2_wp3(**args)
    finally:
        err_info = clubb_api.finalize_stats(err_info=err_info)

    assert isinstance(out, tuple)
    assert len(out) == 7
    for arr in out[:-1]:
        assert np.all(np.isfinite(arr))
    assert isinstance(out[-1], ErrInfo)


def test_advance_wp2_wp3_updates_match_return_values(tmp_path):
    """Returned arrays should match in-place updated input buffers."""
    gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info = _setup_env(tmp_path, ngrdcol=1)
    args = _make_args(gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info)

    up2_in = args["up2"]
    vp2_in = args["vp2"]
    wp2_in = args["wp2"]
    wp3_in = args["wp3"]
    wp3_zm_in = args["wp3_zm"]
    wp2_zt_in = args["wp2_zt"]

    try:
        up2_out, vp2_out, wp2_out, wp3_out, wp3_zm_out, wp2_zt_out, err_info_out = clubb_api.advance_wp2_wp3(**args)
    finally:
        err_info = clubb_api.finalize_stats(err_info=err_info)

    np.testing.assert_allclose(up2_out, up2_in)
    np.testing.assert_allclose(vp2_out, vp2_in)
    np.testing.assert_allclose(wp2_out, wp2_in)
    np.testing.assert_allclose(wp3_out, wp3_in)
    np.testing.assert_allclose(wp3_zm_out, wp3_zm_in)
    np.testing.assert_allclose(wp2_zt_out, wp2_zt_in)
    assert isinstance(err_info_out, ErrInfo)
