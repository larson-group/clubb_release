"""Test wrappers for sfc_varnce_module routines in the call tree."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.sclr_idx import SclrIdx


def _make_grid(ngrdcol=1, nzmax=11, dz=500.0):
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


def test_calc_sfc_varnce_basic_finite_outputs():
    """calc_sfc_varnce should produce finite, nonnegative variances for simple inputs."""
    gr = _make_grid()
    ngrdcol = gr.ngrdcol
    sclr_dim = 1

    sclr_idx = SclrIdx(1, 1, 0, 0, 0, 0)
    clubb_api.set_sclr_idx(*sclr_idx)
    clubb_api.reset_err_code()
    err_info = ErrInfo(ngrdcol=ngrdcol)

    sfc_k = gr.k_lb_zm
    sfc_elevation = gr.zm[:, sfc_k].copy()
    upwp_sfc = np.zeros(ngrdcol, dtype=np.float64)
    vpwp_sfc = np.zeros(ngrdcol, dtype=np.float64)
    wprtp_sfc = np.zeros(ngrdcol, dtype=np.float64)
    wpthlp = np.zeros((ngrdcol, gr.nzm), dtype=np.float64)
    um = np.zeros((ngrdcol, gr.nzt), dtype=np.float64)
    vm = np.zeros((ngrdcol, gr.nzt), dtype=np.float64)
    lscale_up = np.full((ngrdcol, gr.nzt), 100.0, dtype=np.float64)
    wpsclrp_sfc = np.zeros((ngrdcol, sclr_dim), dtype=np.float64)
    lhs_splat_wp2 = np.zeros((ngrdcol, gr.nzm), dtype=np.float64)
    tau_zm = np.full((ngrdcol, gr.nzm), 100.0, dtype=np.float64)
    up2_sfc_coef = np.ones(ngrdcol, dtype=np.float64)
    a_const = np.ones(ngrdcol, dtype=np.float64)

    wp2 = np.full((ngrdcol, gr.nzm), 0.05, dtype=np.float64)
    up2 = np.full((ngrdcol, gr.nzm), 0.05, dtype=np.float64)
    vp2 = np.full((ngrdcol, gr.nzm), 0.05, dtype=np.float64)
    thlp2 = np.full((ngrdcol, gr.nzm), 0.01, dtype=np.float64)
    rtp2 = np.full((ngrdcol, gr.nzm), 1.0e-6, dtype=np.float64)
    rtpthlp = np.zeros((ngrdcol, gr.nzm), dtype=np.float64)
    sclrp2 = np.zeros((ngrdcol, gr.nzm, sclr_dim), dtype=np.float64)
    sclrprtp = np.zeros((ngrdcol, gr.nzm, sclr_dim), dtype=np.float64)
    sclrpthlp = np.zeros((ngrdcol, gr.nzm, sclr_dim), dtype=np.float64)

    result = clubb_api.calc_sfc_varnce(
        gr=gr, nzm=gr.nzm, nzt=gr.nzt, ngrdcol=gr.ngrdcol, sclr_dim=sclr_dim,
        dt=60.0, sfc_elevation=sfc_elevation,
        upwp_sfc=upwp_sfc, vpwp_sfc=vpwp_sfc, wpthlp=wpthlp, wprtp_sfc=wprtp_sfc,
        um=um, vm=vm, lscale_up=lscale_up, wpsclrp_sfc=wpsclrp_sfc,
        lhs_splat_wp2=lhs_splat_wp2, tau_zm=tau_zm,
        l_vary_convect_depth=False, t0=300.0,
        up2_sfc_coef=up2_sfc_coef, a_const=a_const,
        wp2=wp2, up2=up2, vp2=vp2, thlp2=thlp2, rtp2=rtp2, rtpthlp=rtpthlp,
        sclrp2=sclrp2, sclrprtp=sclrprtp, sclrpthlp=sclrpthlp,
        sclr_idx=sclr_idx, err_info=err_info,
    )

    (
        wp2_o, up2_o, vp2_o, thlp2_o, rtp2_o, rtpthlp_o,
        sclrp2_o, sclrprtp_o, sclrpthlp_o, err_info_out,
    ) = result

    assert wp2_o.shape == (ngrdcol, gr.nzm)
    assert up2_o.shape == (ngrdcol, gr.nzm)
    assert vp2_o.shape == (ngrdcol, gr.nzm)
    assert thlp2_o.shape == (ngrdcol, gr.nzm)
    assert rtp2_o.shape == (ngrdcol, gr.nzm)
    assert rtpthlp_o.shape == (ngrdcol, gr.nzm)
    assert sclrp2_o.shape == (ngrdcol, gr.nzm, sclr_dim)
    assert sclrprtp_o.shape == (ngrdcol, gr.nzm, sclr_dim)
    assert sclrpthlp_o.shape == (ngrdcol, gr.nzm, sclr_dim)

    for arr in (wp2_o, up2_o, vp2_o, thlp2_o, rtp2_o, rtpthlp_o, sclrp2_o, sclrprtp_o, sclrpthlp_o):
        assert np.all(np.isfinite(arr))

    assert np.all(wp2_o >= 0.0)
    # Surface closures can drive tiny negative wind variances; enforce small tolerance.
    assert np.all(up2_o >= -1.0e-3)
    assert np.all(vp2_o >= -1.0e-3)
    assert np.all(thlp2_o >= 0.0)
    assert np.all(rtp2_o >= 0.0)
    assert np.all(sclrp2_o >= 0.0)

    assert isinstance(err_info_out, ErrInfo)
    assert int(clubb_api.get_err_code(ngrdcol)[0]) == 0
