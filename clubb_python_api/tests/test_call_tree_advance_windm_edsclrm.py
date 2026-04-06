"""Test wrapper coverage for advance_windm_edsclrm call-tree routine."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def _setup_env(tmp_path: Path, ngrdcol: int = 1):
    """Initialize grid and persistent UDT state used by advance_windm_edsclrm."""
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
    err_info = ErrInfo(ngrdcol=ngrdcol)

    err_info = clubb_api.init_stats(
        registry_path="tests/test_stats_registry.in",
        output_path=str(tmp_path / "windm_edsclrm_stats.nc"),
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
        edsclr_dim=1,
        clubb_params=clubb_params,
        param_names=clubb_api.get_param_names(),
    )

    return gr, flags, nu_vert_res_dep, err_info


def _make_args(gr, flags, nu_vert_res_dep, err_info):
    """Construct stable, finite test inputs for advance_windm_edsclrm."""
    ngrdcol, nzt, nzm = gr.ngrdcol, gr.nzt, gr.nzm

    def full(shape, val):
        return np.full(shape, val, dtype=np.float64, order="F")

    args = {
        "gr": gr,
        "nzm": nzm,
        "nzt": nzt,
        "ngrdcol": ngrdcol,
        "edsclr_dim": 1,
        "dt": 60.0,
        "wm_zt": full((ngrdcol, nzt), 0.0),
        "km_zm": full((ngrdcol, nzm), 1.0),
        "kmh_zm": full((ngrdcol, nzm), 1.0),
        "ug": full((ngrdcol, nzt), 1.0),
        "vg": full((ngrdcol, nzt), 0.0),
        "um_ref": full((ngrdcol, nzt), 1.0),
        "vm_ref": full((ngrdcol, nzt), 0.0),
        "wp2": full((ngrdcol, nzm), 0.01),
        "up2": full((ngrdcol, nzm), 0.01),
        "vp2": full((ngrdcol, nzm), 0.01),
        "um_forcing": full((ngrdcol, nzt), 0.0),
        "vm_forcing": full((ngrdcol, nzt), 0.0),
        "edsclrm_forcing": full((ngrdcol, nzt, 1), 0.0),
        "rho_ds_zm": full((ngrdcol, nzm), 1.0),
        "invrs_rho_ds_zt": full((ngrdcol, nzt), 1.0),
        "fcor": full((ngrdcol,), 1.0e-4),
        "l_implemented": True,
        "ts_nudge": 3600.0,
        "tridiag_solve_method": int(flags.tridiag_solve_method),
        "l_predict_upwp_vpwp": False,
        "l_upwind_xm_ma": bool(flags.l_upwind_xm_ma),
        "l_uv_nudge": bool(flags.l_uv_nudge),
        "l_tke_aniso": bool(flags.l_tke_aniso),
        "l_lmm_stepping": bool(flags.l_lmm_stepping),
        "l_linearize_pbl_winds": bool(flags.l_linearize_pbl_winds),
        "order_xp2_xpyp": 2,
        "order_wp2_wp3": 3,
        "order_windm": 4,
        "um": full((ngrdcol, nzt), 1.0),
        "vm": full((ngrdcol, nzt), 0.0),
        "edsclrm": full((ngrdcol, nzt, 1), 0.0),
        "upwp": full((ngrdcol, nzm), 0.0),
        "vpwp": full((ngrdcol, nzm), 0.0),
        "wpedsclrp": full((ngrdcol, nzm, 1), 0.0),
        "um_pert": full((ngrdcol, nzt), 0.0),
        "vm_pert": full((ngrdcol, nzt), 0.0),
        "upwp_pert": full((ngrdcol, nzm), 0.0),
        "vpwp_pert": full((ngrdcol, nzm), 0.0),
        "nu_vert_res_dep": nu_vert_res_dep,
        "err_info": err_info,
    }
    return args


def test_advance_windm_edsclrm_returns_finite_arrays(tmp_path):
    """advance_windm_edsclrm should run and return finite output arrays."""
    gr, flags, nu_vert_res_dep, err_info = _setup_env(tmp_path, ngrdcol=1)
    args = _make_args(gr, flags, nu_vert_res_dep, err_info)

    try:
        out = clubb_api.advance_windm_edsclrm(**args)
    finally:
        err_info = clubb_api.finalize_stats(err_info=err_info)

    assert isinstance(out, tuple)
    assert len(out) == 11
    for arr in out[:-1]:
        assert np.all(np.isfinite(arr))
    assert isinstance(out[-1], ErrInfo)


def test_advance_windm_edsclrm_updates_match_return_values(tmp_path):
    """Returned arrays should match in-place updated input buffers."""
    gr, flags, nu_vert_res_dep, err_info = _setup_env(tmp_path, ngrdcol=1)
    args = _make_args(gr, flags, nu_vert_res_dep, err_info)

    um_in = args["um"]
    vm_in = args["vm"]
    edsclrm_in = args["edsclrm"]
    upwp_in = args["upwp"]
    vpwp_in = args["vpwp"]
    wpedsclrp_in = args["wpedsclrp"]
    um_pert_in = args["um_pert"]
    vm_pert_in = args["vm_pert"]
    upwp_pert_in = args["upwp_pert"]
    vpwp_pert_in = args["vpwp_pert"]

    try:
        (
            um_out, vm_out, edsclrm_out, upwp_out, vpwp_out,
            wpedsclrp_out, um_pert_out, vm_pert_out, upwp_pert_out, vpwp_pert_out, err_info_out,
        ) = clubb_api.advance_windm_edsclrm(**args)
    finally:
        err_info = clubb_api.finalize_stats(err_info=err_info)

    np.testing.assert_allclose(um_out, um_in)
    np.testing.assert_allclose(vm_out, vm_in)
    np.testing.assert_allclose(edsclrm_out, edsclrm_in)
    np.testing.assert_allclose(upwp_out, upwp_in)
    np.testing.assert_allclose(vpwp_out, vpwp_in)
    np.testing.assert_allclose(wpedsclrp_out, wpedsclrp_in)
    np.testing.assert_allclose(um_pert_out, um_pert_in)
    np.testing.assert_allclose(vm_pert_out, vm_pert_in)
    np.testing.assert_allclose(upwp_pert_out, upwp_pert_in)
    np.testing.assert_allclose(vpwp_pert_out, vpwp_pert_in)
    assert isinstance(err_info_out, ErrInfo)
