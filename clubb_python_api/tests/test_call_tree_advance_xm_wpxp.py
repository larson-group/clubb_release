"""Test wrapper coverage for advance_xm_wpxp call-tree routine."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def _setup_env(tmp_path: Path, ngrdcol: int = 1):
    """Initialize grid and persistent UDT state used by advance_xm_wpxp."""
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
    err_info = ErrInfo(ngrdcol=ngrdcol)
    clubb_api.set_sclr_idx(0, 0, 0, 0, 0, 0)

    err_info = clubb_api.init_stats(
        registry_path="tests/test_stats_registry.in",
        output_path=str(tmp_path / "xm_wpxp_stats.nc"),
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

    return gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info


def _make_args(gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info):
    """Construct stable, finite test inputs for advance_xm_wpxp."""
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
        "sigma_sqd_w": full((ngrdcol, nzm), 0.5),
        "wm_zm": full((ngrdcol, nzm), 0.0),
        "wm_zt": full((ngrdcol, nzt), 0.0),
        "wp2": full((ngrdcol, nzm), 0.01),
        "lscale_zm": full((ngrdcol, nzm), 50.0),
        "wp3_on_wp2": full((ngrdcol, nzm), 0.0),
        "wp3_on_wp2_zt": full((ngrdcol, nzt), 0.0),
        "kh_zt": full((ngrdcol, nzt), 1.0),
        "kh_zm": full((ngrdcol, nzm), 1.0),
        "stability_correction": full((ngrdcol, nzm), 1.0),
        "invrs_tau_c6_zm": full((ngrdcol, nzm), 1.0 / 300.0),
        "tau_max_zm": full((ngrdcol, nzm), 300.0),
        "skw_zm": full((ngrdcol, nzm), 0.0),
        "wp2rtp": full((ngrdcol, nzt), 0.0),
        "rtpthvp": full((ngrdcol, nzm), 0.0),
        "rtm_forcing": full((ngrdcol, nzt), 0.0),
        "wprtp_forcing": full((ngrdcol, nzm), 0.0),
        "rtm_ref": full((ngrdcol, nzt), 0.01),
        "wp2thlp": full((ngrdcol, nzt), 0.0),
        "thlpthvp": full((ngrdcol, nzm), 0.0),
        "thlm_forcing": full((ngrdcol, nzt), 0.0),
        "wpthlp_forcing": full((ngrdcol, nzm), 0.0),
        "thlm_ref": full((ngrdcol, nzt), 300.0),
        "rho_ds_zm": full((ngrdcol, nzm), 1.0),
        "rho_ds_zt": full((ngrdcol, nzt), 1.0),
        "invrs_rho_ds_zm": full((ngrdcol, nzm), 1.0),
        "invrs_rho_ds_zt": full((ngrdcol, nzt), 1.0),
        "thv_ds_zm": full((ngrdcol, nzm), 300.0),
        "rtp2": full((ngrdcol, nzm), 1.0e-6),
        "thlp2": full((ngrdcol, nzm), 0.01),
        "w_1_zm": full((ngrdcol, nzm), 0.0),
        "w_2_zm": full((ngrdcol, nzm), 0.0),
        "varnce_w_1_zm": full((ngrdcol, nzm), 0.01),
        "varnce_w_2_zm": full((ngrdcol, nzm), 0.01),
        "mixt_frac_zm": full((ngrdcol, nzm), 0.5),
        "l_implemented": True,
        "em": full((ngrdcol, nzm), 0.05),
        "wp2sclrp": full((ngrdcol, nzt, 1), 0.0),
        "sclrpthvp": full((ngrdcol, nzm, 1), 0.0),
        "sclrm_forcing": full((ngrdcol, nzt, 1), 0.0),
        "sclrp2": full((ngrdcol, nzm, 1), 1.0e-6),
        "cx_fnc_richardson": full((ngrdcol, nzm), 1.0),
        "um_forcing": full((ngrdcol, nzt), 0.0),
        "vm_forcing": full((ngrdcol, nzt), 0.0),
        "ug": full((ngrdcol, nzt), 1.0),
        "vg": full((ngrdcol, nzt), 0.0),
        "wpthvp": full((ngrdcol, nzm), 0.0),
        "fcor": full((ngrdcol,), 1.0e-4),
        "fcor_y": full((ngrdcol,), 0.0),
        "um_ref": full((ngrdcol, nzt), 1.0),
        "vm_ref": full((ngrdcol, nzt), 0.0),
        "up2": full((ngrdcol, nzm), 0.01),
        "vp2": full((ngrdcol, nzm), 0.01),
        "uprcp": full((ngrdcol, nzm), 0.0),
        "vprcp": full((ngrdcol, nzm), 0.0),
        "rc_coef_zm": full((ngrdcol, nzm), 0.0),
        "clubb_params": clubb_params,
        "ts_nudge": 3600.0,
        "iipdf_type": int(flags.iiPDF_type),
        "penta_solve_method": int(flags.penta_solve_method),
        "tridiag_solve_method": int(flags.tridiag_solve_method),
        "fill_holes_type": int(flags.fill_holes_type),
        "l_predict_upwp_vpwp": bool(flags.l_predict_upwp_vpwp),
        "l_ho_nontrad_coriolis": bool(flags.l_ho_nontrad_coriolis),
        "l_ho_trad_coriolis": bool(flags.l_ho_trad_coriolis),
        "l_diffuse_rtm_and_thlm": bool(flags.l_diffuse_rtm_and_thlm),
        "l_stability_correct_kh_n2_zm": bool(flags.l_stability_correct_Kh_N2_zm),
        "l_godunov_upwind_wpxp_ta": bool(flags.l_godunov_upwind_wpxp_ta),
        "l_upwind_xm_ma": bool(flags.l_upwind_xm_ma),
        "l_uv_nudge": bool(flags.l_uv_nudge),
        "l_tke_aniso": bool(flags.l_tke_aniso),
        "l_diag_lscale_from_tau": bool(flags.l_diag_Lscale_from_tau),
        "l_use_c7_richardson": bool(flags.l_use_C7_Richardson),
        "l_lmm_stepping": bool(flags.l_lmm_stepping),
        "l_enable_relaxed_clipping": bool(flags.l_enable_relaxed_clipping),
        "l_linearize_pbl_winds": bool(flags.l_linearize_pbl_winds),
        "l_mono_flux_lim_thlm": bool(flags.l_mono_flux_lim_thlm),
        "l_mono_flux_lim_rtm": bool(flags.l_mono_flux_lim_rtm),
        "l_mono_flux_lim_um": bool(flags.l_mono_flux_lim_um),
        "l_mono_flux_lim_vm": bool(flags.l_mono_flux_lim_vm),
        "l_mono_flux_lim_spikefix": bool(flags.l_mono_flux_lim_spikefix),
        "order_xm_wpxp": 1,
        "order_xp2_xpyp": 2,
        "order_wp2_wp3": 3,
        "rtm": full((ngrdcol, nzt), 0.01),
        "wprtp": full((ngrdcol, nzm), 0.0),
        "thlm": full((ngrdcol, nzt), 300.0),
        "wpthlp": full((ngrdcol, nzm), 0.0),
        "sclrm": full((ngrdcol, nzt, 1), 0.0),
        "wpsclrp": full((ngrdcol, nzm, 1), 0.0),
        "um": full((ngrdcol, nzt), 1.0),
        "upwp": full((ngrdcol, nzm), 0.0),
        "vm": full((ngrdcol, nzt), 0.0),
        "vpwp": full((ngrdcol, nzm), 0.0),
        "um_pert": full((ngrdcol, nzt), 0.0),
        "vm_pert": full((ngrdcol, nzt), 0.0),
        "upwp_pert": full((ngrdcol, nzm), 0.0),
        "vpwp_pert": full((ngrdcol, nzm), 0.0),
        "nu_vert_res_dep": nu_vert_res_dep,
        "pdf_implicit_coefs_terms": pdf_implicit_coefs_terms,
        "err_info": err_info,
    }
    return args


def test_advance_xm_wpxp_returns_finite_arrays(tmp_path):
    """advance_xm_wpxp should run and return finite output arrays."""
    gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info = _setup_env(tmp_path, ngrdcol=1)
    args = _make_args(gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info)

    try:
        out = clubb_api.advance_xm_wpxp(**args)
    finally:
        err_info = clubb_api.finalize_stats(err_info=err_info)

    assert isinstance(out, tuple)
    assert len(out) == 15
    for arr in out[:-1]:
        assert np.all(np.isfinite(arr))
    assert isinstance(out[-1], ErrInfo)


def test_advance_xm_wpxp_updates_match_return_values(tmp_path):
    """Returned arrays should match in-place updated input buffers."""
    gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info = _setup_env(tmp_path, ngrdcol=1)
    args = _make_args(gr, flags, clubb_params, nu_vert_res_dep, pdf_implicit_coefs_terms, err_info)

    rtm_in = args["rtm"]
    wprtp_in = args["wprtp"]
    thlm_in = args["thlm"]
    wpthlp_in = args["wpthlp"]
    sclrm_in = args["sclrm"]
    wpsclrp_in = args["wpsclrp"]
    um_in = args["um"]
    upwp_in = args["upwp"]
    vm_in = args["vm"]
    vpwp_in = args["vpwp"]
    um_pert_in = args["um_pert"]
    vm_pert_in = args["vm_pert"]
    upwp_pert_in = args["upwp_pert"]
    vpwp_pert_in = args["vpwp_pert"]

    try:
        (
            rtm_out, wprtp_out, thlm_out, wpthlp_out,
            sclrm_out, wpsclrp_out, um_out, upwp_out, vm_out, vpwp_out,
            um_pert_out, vm_pert_out, upwp_pert_out, vpwp_pert_out, err_info_out,
        ) = clubb_api.advance_xm_wpxp(**args)
    finally:
        err_info = clubb_api.finalize_stats(err_info=err_info)

    np.testing.assert_allclose(rtm_out, rtm_in)
    np.testing.assert_allclose(wprtp_out, wprtp_in)
    np.testing.assert_allclose(thlm_out, thlm_in)
    np.testing.assert_allclose(wpthlp_out, wpthlp_in)
    np.testing.assert_allclose(sclrm_out, sclrm_in)
    np.testing.assert_allclose(wpsclrp_out, wpsclrp_in)
    np.testing.assert_allclose(um_out, um_in)
    np.testing.assert_allclose(upwp_out, upwp_in)
    np.testing.assert_allclose(vm_out, vm_in)
    np.testing.assert_allclose(vpwp_out, vpwp_in)
    np.testing.assert_allclose(um_pert_out, um_pert_in)
    np.testing.assert_allclose(vm_pert_out, vm_pert_in)
    np.testing.assert_allclose(upwp_pert_out, upwp_pert_in)
    np.testing.assert_allclose(vpwp_pert_out, vpwp_pert_in)
    assert isinstance(err_info_out, ErrInfo)
