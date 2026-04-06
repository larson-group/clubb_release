"""Test wrapper coverage for advance_clubb_core call-tree entrypoint."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.pdf_params import implicit_coefs_terms, pdf_parameter
from clubb_python.derived_types.sclr_idx import SclrIdx


def _make_advance_inputs(ngrdcol: int):
    """Build a minimal, physically sane argument set for advance_clubb_core."""
    nzmax = 5
    nzt = nzmax - 1
    nzm = nzmax
    dz = 100.0

    clubb_api.init_err_info(ngrdcol)
    config_flags = clubb_api.get_default_config_flags()
    clubb_api.init_config_flags(config_flags)

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
    nu_vert_res_dep, lmin, mixt_frac_max_mag = clubb_api.calc_derrived_params(
        gr=gr,
        ngrdcol=ngrdcol,
        grid_type=1,
        deltaz=np.full(ngrdcol, dz, dtype=np.float64),
        clubb_params=clubb_params,
        nu_vert_res_dep=None,
        l_prescribed_avg_deltaz=False,
    )

    pdf_params = clubb_api.init_pdf_params(nzt, ngrdcol)
    pdf_params_zm = clubb_api.init_pdf_params_zm(nzm, ngrdcol)
    pdf_implicit_coefs_terms = clubb_api.init_pdf_implicit(nzm, ngrdcol, sclr_dim=1)
    sclr_idx = SclrIdx(0, 0, 0, 0, 0, 0)
    clubb_api.set_sclr_idx(*sclr_idx)
    err_info = ErrInfo(ngrdcol=ngrdcol)

    def full(shape, val):
        return np.full(shape, val, dtype=np.float64, order="F")

    shzt = (ngrdcol, nzt)
    shzm = (ngrdcol, nzm)
    shzts = (ngrdcol, nzt, 1)
    shzms = (ngrdcol, nzm, 1)

    return {
        "gr": gr,
        "nzm": nzm,
        "nzt": nzt,
        "ngrdcol": ngrdcol,
        "dt": 60.0,
        "l_implemented": True,
        "fcor": full((ngrdcol,), 1.0e-4),
        "fcor_y": full((ngrdcol,), 0.0),
        "sfc_elevation": full((ngrdcol,), 0.0),
        "hydromet_dim": 1,
        "sclr_dim": 1,
        "edsclr_dim": 1,
        "sclr_tol": np.array([1.0e-8], dtype=np.float64, order="F"),
        "thlm_forcing": full(shzt, 0.0),
        "rtm_forcing": full(shzt, 0.0),
        "um_forcing": full(shzt, 0.0),
        "vm_forcing": full(shzt, 0.0),
        "wm_zt": full(shzt, 0.0),
        "rho": full(shzt, 1.0),
        "rho_ds_zt": full(shzt, 1.0),
        "invrs_rho_ds_zt": full(shzt, 1.0),
        "thv_ds_zt": full(shzt, 300.0),
        "rfrzm": full(shzt, 0.0),
        "wprtp_forcing": full(shzm, 0.0),
        "wpthlp_forcing": full(shzm, 0.0),
        "rtp2_forcing": full(shzm, 0.0),
        "thlp2_forcing": full(shzm, 0.0),
        "rtpthlp_forcing": full(shzm, 0.0),
        "wm_zm": full(shzm, 0.0),
        "rho_zm": full(shzm, 1.0),
        "rho_ds_zm": full(shzm, 1.0),
        "invrs_rho_ds_zm": full(shzm, 1.0),
        "thv_ds_zm": full(shzm, 300.0),
        "wpthlp_sfc": full((ngrdcol,), 0.0),
        "wprtp_sfc": full((ngrdcol,), 0.0),
        "upwp_sfc": full((ngrdcol,), 0.0),
        "vpwp_sfc": full((ngrdcol,), 0.0),
        "p_sfc": full((ngrdcol,), 101325.0),
        "upwp_sfc_pert": full((ngrdcol,), 0.0),
        "vpwp_sfc_pert": full((ngrdcol,), 0.0),
        "rtm_ref": full(shzt, 0.01),
        "thlm_ref": full(shzt, 300.0),
        "um_ref": full(shzt, 1.0),
        "vm_ref": full(shzt, 0.0),
        "ug": full(shzt, 1.0),
        "vg": full(shzt, 0.0),
        "host_dx": full((ngrdcol,), 1000.0),
        "host_dy": full((ngrdcol,), 1000.0),
        "clubb_params": clubb_params,
        "lmin": lmin,
        "mixt_frac_max_mag": mixt_frac_max_mag,
        "t0_val": 300.0,
        "ts_nudge": 0.0,
        "rtm_min": 1.0e-8,
        "rtm_nudge_max_altitude": 5000.0,
        "l_mix_rat_hm": np.array([True], dtype=bool),
        "wphydrometp": full((ngrdcol, nzm, 1), 0.0),
        "wp2hmp": full((ngrdcol, nzt, 1), 0.0),
        "rtphmp_zt": full((ngrdcol, nzt, 1), 0.0),
        "thlphmp_zt": full((ngrdcol, nzt, 1), 0.0),
        "sclrm_forcing": full(shzts, 0.0),
        "edsclrm_forcing": full(shzts, 0.0),
        "wpsclrp_sfc": full((ngrdcol, 1), 0.0),
        "wpedsclrp_sfc": full((ngrdcol, 1), 0.0),
        "um": full(shzt, 1.0),
        "vm": full(shzt, 0.0),
        "up3": full(shzt, 0.0),
        "vp3": full(shzt, 0.0),
        "rtm": full(shzt, 0.01),
        "thlm": full(shzt, 300.0),
        "rtp3": full(shzt, 0.0),
        "thlp3": full(shzt, 0.0),
        "wp3": full(shzt, 0.0),
        "p_in_Pa": full(shzt, 90000.0),
        "exner": full(shzt, 1.0),
        "rcm": full(shzt, 0.0),
        "cloud_frac": full(shzt, 0.0),
        "wp2thvp": full(shzt, 0.0),
        "wp2up": full(shzt, 0.0),
        "wp2rtp": full(shzt, 0.0),
        "wp2thlp": full(shzt, 0.0),
        "wpup2": full(shzt, 0.0),
        "wpvp2": full(shzt, 0.0),
        "ice_supersat_frac": full(shzt, 0.0),
        "um_pert": full(shzt, 0.0),
        "vm_pert": full(shzt, 0.0),
        "upwp": full(shzm, 0.0),
        "vpwp": full(shzm, 0.0),
        "up2": full(shzm, 0.01),
        "vp2": full(shzm, 0.01),
        "wprtp": full(shzm, 0.0),
        "wpthlp": full(shzm, 0.0),
        "rtp2": full(shzm, 1.0e-6),
        "thlp2": full(shzm, 0.01),
        "rtpthlp": full(shzm, 0.0),
        "wp2": full(shzm, 0.01),
        "wpthvp": full(shzm, 0.0),
        "rtpthvp": full(shzm, 0.0),
        "thlpthvp": full(shzm, 0.0),
        "uprcp": full(shzm, 0.0),
        "vprcp": full(shzm, 0.0),
        "rc_coef_zm": full(shzm, 0.0),
        "wp4": full(shzm, 0.0),
        "wp2up2": full(shzm, 0.0),
        "wp2vp2": full(shzm, 0.0),
        "upwp_pert": full(shzm, 0.0),
        "vpwp_pert": full(shzm, 0.0),
        "sclrm": full(shzts, 0.0),
        "sclrp3": full(shzts, 0.0),
        "wpsclrp": full(shzms, 0.0),
        "sclrp2": full(shzms, 0.0),
        "sclrprtp": full(shzms, 0.0),
        "sclrpthlp": full(shzms, 0.0),
        "sclrpthvp": full(shzms, 0.0),
        "edsclrm": full(shzts, 0.0),
        "sclr_idx": sclr_idx,
        "config_flags": config_flags,
        "nu_vert_res_dep": nu_vert_res_dep,
        "pdf_params": pdf_params,
        "pdf_params_zm": pdf_params_zm,
        "pdf_implicit_coefs_terms": pdf_implicit_coefs_terms,
        "err_info": err_info,
    }


def test_advance_clubb_core_returns_finite_arrays():
    """advance_clubb_core should run and return finite arrays for a minimal input state."""
    args = _make_advance_inputs(ngrdcol=1)
    out = clubb_api.advance_clubb_core(**args)

    assert isinstance(out, tuple)
    assert len(out) > 10
    for arr in out[:51]:
        assert np.all(np.isfinite(arr))
    for arr in out[55:]:
        assert np.all(np.isfinite(arr))
    assert all(
        isinstance(udt, expected_type)
        for udt, expected_type in zip(
            out[51:55],
            (pdf_parameter, pdf_parameter, implicit_coefs_terms, ErrInfo),
        )
    )


def test_advance_clubb_core_repeatability_for_identical_inputs():
    """Calling advance_clubb_core twice with identical inputs should yield identical outputs."""
    args1 = _make_advance_inputs(ngrdcol=1)
    args2 = _make_advance_inputs(ngrdcol=1)

    out1 = clubb_api.advance_clubb_core(**args1)
    out2 = clubb_api.advance_clubb_core(**args2)

    assert len(out1) == len(out2)
    for arr1, arr2 in zip(out1[:51], out2[:51]):
        np.testing.assert_allclose(arr1, arr2)
    for arr1, arr2 in zip(out1[55:], out2[55:]):
        np.testing.assert_allclose(arr1, arr2)
