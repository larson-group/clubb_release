"""Test wrapper for precipitation_fraction::precip_fraction."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def _make_grid(ngrdcol=1, nzmax=6, dz=250.0):
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


def test_precip_fraction_returns_zero_when_no_hydrometeors(tmp_path):
    """With zero hydrometeors everywhere, precip fractions should remain zero."""
    gr = _make_grid()
    ngrdcol, nzt = gr.ngrdcol, gr.nzt
    hydromet_dim = 1
    err_info = ErrInfo(ngrdcol=ngrdcol)

    registry = str(Path(__file__).resolve().parent / "test_stats_registry.in")
    clubb_params = clubb_api.init_clubb_params(ngrdcol, iunit=10, filename="")

    clubb_api.init_stats(
        registry_path=registry,
        output_path=str(tmp_path / "precip_frac_stats.nc"),
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
        nzm=gr.nzm,
        zm=np.asfortranarray(gr.zm),
        err_info=err_info,
        sclr_dim=0,
        edsclr_dim=0,
        clubb_params=clubb_params,
        param_names=clubb_api.get_param_names(),
    )

    try:
        zeros_2d = np.zeros((ngrdcol, nzt), dtype=np.float64, order="F")
        pulled_err, precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol = clubb_api.precip_fraction(
            gr=gr,
            nzt=nzt,
            ngrdcol=ngrdcol,
            hydromet_dim=hydromet_dim,
            hydromet=np.zeros((ngrdcol, nzt, hydromet_dim), dtype=np.float64, order="F"),
            cloud_frac=np.full((ngrdcol, nzt), 0.3, dtype=np.float64, order="F"),
            cloud_frac_1=zeros_2d,
            l_mix_rat_hm=np.array([True], dtype=bool),
            l_frozen_hm=np.array([False], dtype=bool),
            hydromet_tol=np.array([1.0e-8], dtype=np.float64),
            cloud_frac_2=zeros_2d,
            ice_supersat_frac=zeros_2d,
            ice_supersat_frac_1=zeros_2d,
            ice_supersat_frac_2=zeros_2d,
            mixt_frac=np.full((ngrdcol, nzt), 0.5, dtype=np.float64, order="F"),
            clubb_params=clubb_params[0, :],
            err_info=err_info,
        )
    finally:
        clubb_api.finalize_stats(err_info=err_info)

    assert precip_frac.shape == (ngrdcol, nzt)
    assert precip_frac_1.shape == (ngrdcol, nzt)
    assert precip_frac_2.shape == (ngrdcol, nzt)
    assert precip_frac_tol.shape == (ngrdcol,)
    assert np.allclose(precip_frac, 0.0)
    assert np.allclose(precip_frac_1, 0.0)
    assert np.allclose(precip_frac_2, 0.0)
    assert np.all(precip_frac_tol > 0.0)
    assert np.all(pulled_err.err_code == 0)
