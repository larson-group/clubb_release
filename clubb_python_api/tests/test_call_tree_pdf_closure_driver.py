"""Test wrapper coverage for pdf_closure_driver call-tree routine."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def _setup_grid(ngrdcol: int = 1):
    """Initialize grid and err_info state used by pdf_closure_driver wrapper."""
    nzmax = 5
    dz = 100.0

    clubb_api.init_err_info(ngrdcol)

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
    return gr


def _make_driver_inputs(gr, sclr_dim: int, hydromet_dim: int):
    """Create synthetic state arrays for pdf_closure_driver."""
    ncol, nzm, nzt = gr.ngrdcol, gr.nzm, gr.nzt

    def fill(shape, value):
        return np.full(shape, value, dtype=np.float64, order="F")

    return {
        "wprtp": fill((ncol, nzm), 2.0e-3),
        "thlm": fill((ncol, nzt), 300.0),
        "wpthlp": fill((ncol, nzm), 1.0e-1),
        "rtp2": fill((ncol, nzm), 1.0e-4),
        "rtp3": fill((ncol, nzt), 1.0e-6),
        "thlp2": fill((ncol, nzm), 5.0e-1),
        "thlp3": fill((ncol, nzt), 1.0e-3),
        "rtpthlp": fill((ncol, nzm), 2.0e-3),
        "wp2": fill((ncol, nzm), 2.0e-1),
        "wp3": fill((ncol, nzt), 1.0e-2),
        "wm_zm": np.zeros((ncol, nzm), dtype=np.float64, order="F"),
        "wm_zt": np.zeros((ncol, nzt), dtype=np.float64, order="F"),
        "um": np.zeros((ncol, nzt), dtype=np.float64, order="F"),
        "up2": fill((ncol, nzm), 2.0e-1),
        "upwp": fill((ncol, nzm), 1.0e-2),
        "up3": fill((ncol, nzt), 1.0e-2),
        "vm": np.zeros((ncol, nzt), dtype=np.float64, order="F"),
        "vp2": fill((ncol, nzm), 2.0e-1),
        "vpwp": fill((ncol, nzm), -1.0e-2),
        "vp3": fill((ncol, nzt), 1.0e-2),
        "p_in_pa": fill((ncol, nzt), 9.0e4),
        "exner": fill((ncol, nzt), 1.0),
        "thv_ds_zm": fill((ncol, nzm), 300.0),
        "thv_ds_zt": fill((ncol, nzt), 300.0),
        "rtm_ref": fill((ncol, nzt), 1.0e-2),
        "wphydrometp": np.zeros((ncol, nzm, hydromet_dim), dtype=np.float64, order="F"),
        "wp2hmp": np.zeros((ncol, nzt, hydromet_dim), dtype=np.float64, order="F"),
        "rtphmp_zt": np.zeros((ncol, nzt, hydromet_dim), dtype=np.float64, order="F"),
        "thlphmp_zt": np.zeros((ncol, nzt, hydromet_dim), dtype=np.float64, order="F"),
        "sclrm": np.zeros((ncol, nzt, sclr_dim), dtype=np.float64, order="F"),
        "wpsclrp": fill((ncol, nzm, sclr_dim), 2.0e-3),
        "sclrp2": fill((ncol, nzm, sclr_dim), 1.0e-4),
        "sclrprtp": fill((ncol, nzm, sclr_dim), 5.0e-5),
        "sclrpthlp": fill((ncol, nzm, sclr_dim), 2.0e-3),
        "sclrp3": fill((ncol, nzt, sclr_dim), 1.0e-6),
        "p_sfc": fill((ncol,), 1.0e5),
        "rtm": fill((ncol, nzt), 1.0e-2),
    }


def _driver_controls(gr, sclr_dim: int, hydromet_dim: int):
    """Build explicit control inputs for strict direct pdf_closure_driver mapping."""
    return {
        "sclr_tol": np.asfortranarray(np.full((sclr_dim,), 1.0e-8, dtype=np.float64)),
        "l_mix_rat_hm": np.zeros((hydromet_dim,), dtype=bool),
        "clubb_params": np.asfortranarray(clubb_api.init_clubb_params(gr.ngrdcol, iunit=10, filename=""), dtype=np.float64),
        "mixt_frac_max_mag": 0.99,
        "ts_nudge": 3600.0,
        "rtm_min": 1.0e-8,
        "rtm_nudge_max_altitude": 0.0,
        "l_samp_stats_in_pdf_call": False,
        "l_rtm_nudge": False,
        "l_trapezoidal_rule_zt": False,
        "l_trapezoidal_rule_zm": False,
        "l_use_cloud_cover": False,
        "l_rcm_supersat_adj": False,
    }


def test_pdf_closure_driver_returns_finite_arrays():
    """pdf_closure_driver should return finite cloud/rcm/wp4 fields."""
    gr = _setup_grid(ngrdcol=1)
    flags = clubb_api.get_default_config_flags()
    inputs = _make_driver_inputs(gr, sclr_dim=1, hydromet_dim=1)
    controls = _driver_controls(gr, sclr_dim=1, hydromet_dim=1)
    pdf_params = clubb_api.init_pdf_params(gr.nzt, gr.ngrdcol)
    pdf_params_zm = clubb_api.init_pdf_params_zm(gr.nzm, gr.ngrdcol)
    implicit = clubb_api.init_pdf_implicit(gr.nzt, gr.ngrdcol, sclr_dim=1)
    err_info = ErrInfo(ngrdcol=gr.ngrdcol)
    clubb_api.reset_err_code()

    out = clubb_api.pdf_closure_driver(
        gr=gr,
        nzm=gr.nzm,
        nzt=gr.nzt,
        ngrdcol=gr.ngrdcol,
        dt=60.0,
        sclr_dim=1,
        hydromet_dim=1,
        saturation_formula=int(flags.saturation_formula),
        iiPDF_type=int(flags.iiPDF_type),
        l_predict_upwp_vpwp=bool(flags.l_predict_upwp_vpwp),
        l_call_pdf_closure_twice=bool(flags.l_call_pdf_closure_twice),
        pdf_params=pdf_params,
        pdf_params_zm=pdf_params_zm,
        pdf_implicit_coefs_terms=implicit,
        err_info=err_info,
        **controls,
        **inputs,
    )
    rtm_out = out[0]
    _implicit_out, _pdf_params_out, _pdf_params_zm_out, _err_info_out = out[1:5]
    rcm, cloud_frac, wp4 = out[5], out[6], out[26]

    assert cloud_frac.shape == (gr.ngrdcol, gr.nzt)
    assert rcm.shape == (gr.ngrdcol, gr.nzt)
    assert wp4.shape == (gr.ngrdcol, gr.nzm)

    assert np.all(np.isfinite(cloud_frac))
    assert np.all(np.isfinite(rcm))
    assert np.all(np.isfinite(wp4))


def test_pdf_closure_driver_is_repeatable():
    """Repeated calls with identical setup should produce identical outputs."""
    gr = _setup_grid(ngrdcol=1)
    flags = clubb_api.get_default_config_flags()
    inputs1 = _make_driver_inputs(gr, sclr_dim=1, hydromet_dim=1)
    inputs2 = _make_driver_inputs(gr, sclr_dim=1, hydromet_dim=1)
    controls = _driver_controls(gr, sclr_dim=1, hydromet_dim=1)
    pdf_params = clubb_api.init_pdf_params(gr.nzt, gr.ngrdcol)
    pdf_params_zm = clubb_api.init_pdf_params_zm(gr.nzm, gr.ngrdcol)
    implicit = clubb_api.init_pdf_implicit(gr.nzt, gr.ngrdcol, sclr_dim=1)
    err_info = ErrInfo(ngrdcol=gr.ngrdcol)

    clubb_api.reset_err_code()
    out1 = clubb_api.pdf_closure_driver(
        gr=gr,
        nzm=gr.nzm,
        nzt=gr.nzt,
        ngrdcol=gr.ngrdcol,
        dt=60.0,
        sclr_dim=1,
        hydromet_dim=1,
        saturation_formula=int(flags.saturation_formula),
        iiPDF_type=int(flags.iiPDF_type),
        l_predict_upwp_vpwp=bool(flags.l_predict_upwp_vpwp),
        l_call_pdf_closure_twice=bool(flags.l_call_pdf_closure_twice),
        pdf_params=pdf_params,
        pdf_params_zm=pdf_params_zm,
        pdf_implicit_coefs_terms=implicit,
        err_info=err_info,
        **controls,
        **inputs1,
    )
    result1 = out1

    clubb_api.reset_err_code()
    out2 = clubb_api.pdf_closure_driver(
        gr=gr,
        nzm=gr.nzm,
        nzt=gr.nzt,
        ngrdcol=gr.ngrdcol,
        dt=60.0,
        sclr_dim=1,
        hydromet_dim=1,
        saturation_formula=int(flags.saturation_formula),
        iiPDF_type=int(flags.iiPDF_type),
        l_predict_upwp_vpwp=bool(flags.l_predict_upwp_vpwp),
        l_call_pdf_closure_twice=bool(flags.l_call_pdf_closure_twice),
        pdf_params=pdf_params,
        pdf_params_zm=pdf_params_zm,
        pdf_implicit_coefs_terms=implicit,
        err_info=err_info,
        **controls,
        **inputs2,
    )
    result2 = out2

    np.testing.assert_allclose(result1[5], result2[5])
    np.testing.assert_allclose(result1[6], result2[6])
    np.testing.assert_allclose(result1[26], result2[26])


def test_pdf_closure_driver_allows_zero_dims():
    """Wrapper should support dim=0 calls (with internal size-1 buffers)."""
    gr = _setup_grid(ngrdcol=1)
    flags = clubb_api.get_default_config_flags()
    inputs = _make_driver_inputs(gr, sclr_dim=0, hydromet_dim=0)
    controls = _driver_controls(gr, sclr_dim=0, hydromet_dim=0)
    pdf_params = clubb_api.init_pdf_params(gr.nzt, gr.ngrdcol)
    pdf_params_zm = clubb_api.init_pdf_params_zm(gr.nzm, gr.ngrdcol)
    implicit = clubb_api.init_pdf_implicit(gr.nzt, gr.ngrdcol, sclr_dim=0)
    err_info = ErrInfo(ngrdcol=gr.ngrdcol)
    clubb_api.reset_err_code()

    out = clubb_api.pdf_closure_driver(
        gr=gr,
        nzm=gr.nzm,
        nzt=gr.nzt,
        ngrdcol=gr.ngrdcol,
        dt=60.0,
        sclr_dim=0,
        hydromet_dim=0,
        saturation_formula=int(flags.saturation_formula),
        iiPDF_type=int(flags.iiPDF_type),
        l_predict_upwp_vpwp=bool(flags.l_predict_upwp_vpwp),
        l_call_pdf_closure_twice=bool(flags.l_call_pdf_closure_twice),
        pdf_params=pdf_params,
        pdf_params_zm=pdf_params_zm,
        pdf_implicit_coefs_terms=implicit,
        err_info=err_info,
        **controls,
        **inputs,
    )
    rtm_out = out[0]
    _implicit_out, _pdf_params_out, _pdf_params_zm_out, _err_info_out = out[1:5]
    rcm, cloud_frac, wp4 = out[5], out[6], out[26]

    assert cloud_frac.shape == (gr.ngrdcol, gr.nzt)
    assert rcm.shape == (gr.ngrdcol, gr.nzt)
    assert wp4.shape == (gr.ngrdcol, gr.nzm)
    assert np.all(np.isfinite(cloud_frac))
    assert np.all(np.isfinite(rcm))
    assert np.all(np.isfinite(wp4))


def test_pdf_closure_driver_can_return_pdf_udt_state():
    """pdf_closure_driver should support push/pull of PDF UDT globals."""
    gr = _setup_grid(ngrdcol=1)
    flags = clubb_api.get_default_config_flags()
    inputs = _make_driver_inputs(gr, sclr_dim=1, hydromet_dim=1)
    controls = _driver_controls(gr, sclr_dim=1, hydromet_dim=1)
    clubb_api.reset_err_code()

    pdf_params = clubb_api.init_pdf_params(gr.nzt, gr.ngrdcol)
    pdf_params_zm = clubb_api.init_pdf_params_zm(gr.nzm, gr.ngrdcol)
    implicit = clubb_api.init_pdf_implicit(gr.nzt, gr.ngrdcol, sclr_dim=1)
    err_info = ErrInfo(ngrdcol=gr.ngrdcol)

    out = clubb_api.pdf_closure_driver(
        gr=gr,
        nzm=gr.nzm,
        nzt=gr.nzt,
        ngrdcol=gr.ngrdcol,
        dt=60.0,
        sclr_dim=1,
        hydromet_dim=1,
        saturation_formula=int(flags.saturation_formula),
        iiPDF_type=int(flags.iiPDF_type),
        l_predict_upwp_vpwp=bool(flags.l_predict_upwp_vpwp),
        l_call_pdf_closure_twice=bool(flags.l_call_pdf_closure_twice),
        **controls,
        pdf_params=pdf_params,
        pdf_params_zm=pdf_params_zm,
        pdf_implicit_coefs_terms=implicit,
        err_info=err_info,
        **inputs,
    )
    rtm_out = out[0]
    implicit_out, pdf_params_out, pdf_params_zm_out, err_info_out = out[1:5]
    rcm, cloud_frac, wp4 = out[5], out[6], out[26]

    assert cloud_frac.shape == (gr.ngrdcol, gr.nzt)
    assert rcm.shape == (gr.ngrdcol, gr.nzt)
    assert wp4.shape == (gr.ngrdcol, gr.nzm)

    assert pdf_params_out.ngrdcol == gr.ngrdcol
    assert pdf_params_out.nz == gr.nzt
    assert pdf_params_out.w_1.shape == (gr.ngrdcol, gr.nzt)
    assert pdf_params_out.mixt_frac.shape == (gr.ngrdcol, gr.nzt)
    assert np.any(np.abs(pdf_params_out.w_1) > 0.0)

    assert pdf_params_zm_out.ngrdcol == gr.ngrdcol
    assert pdf_params_zm_out.nz == gr.nzm
    assert pdf_params_zm_out.w_1.shape == (gr.ngrdcol, gr.nzm)

    assert implicit_out.ngrdcol == gr.ngrdcol
    assert implicit_out.nz == gr.nzt
    assert implicit_out.sclr_dim == 1
    assert implicit_out.coef_wp4_implicit.shape == (gr.ngrdcol, gr.nzt)
    assert implicit_out.coef_wp2sclrp_implicit is not None
    assert implicit_out.coef_wp2sclrp_implicit.shape == (gr.ngrdcol, gr.nzt, 1)
    assert err_info_out.ngrdcol == gr.ngrdcol
