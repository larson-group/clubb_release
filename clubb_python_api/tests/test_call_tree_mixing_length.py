"""Test wrappers for additional mixing-length and PDF implicit routines."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.pdf_params_converter import (
    get_fortran_implicit_coefs,
    set_fortran_implicit_coefs,
)
from clubb_python.derived_types.pdf_params import (
    unpack_implicit_coefs,
    pack_implicit_coefs_2d,
    pack_implicit_coefs_3d,
)


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


@pytest.fixture(scope="module")
def gr():
    """Shared test grid."""
    return _make_grid()


def test_zero_pdf_implicit_coefs_terms_resets_state():
    """zero_pdf_implicit_coefs_terms should zero all packed implicit arrays."""
    ic = clubb_api.init_pdf_implicit(nz=10, ngrdcol=1, sclr_dim=2)
    data_2d = np.asfortranarray(np.full_like(pack_implicit_coefs_2d(ic), 7.0))
    data_3d = pack_implicit_coefs_3d(ic)
    assert data_3d is not None
    data_3d = np.asfortranarray(np.full_like(data_3d, -3.0))
    ic = unpack_implicit_coefs(data_2d=data_2d, data_3d=data_3d)
    set_fortran_implicit_coefs(ic)

    clubb_api.zero_pdf_implicit_coefs_terms()
    got = get_fortran_implicit_coefs()

    np.testing.assert_allclose(pack_implicit_coefs_2d(got), 0.0)
    got_3d = pack_implicit_coefs_3d(got)
    assert got_3d is not None
    np.testing.assert_allclose(got_3d, 0.0)


def test_diagnose_lscale_from_tau_basic(gr):
    """diagnose_lscale_from_tau should return finite arrays with expected shapes."""
    ngrdcol, nzm, nzt = gr.ngrdcol, gr.nzm, gr.nzt
    clubb_params = clubb_api.init_clubb_params(ngrdcol, iunit=10, filename="")
    err_info = ErrInfo(ngrdcol=ngrdcol)

    out = clubb_api.diagnose_lscale_from_tau(
        gr,
        nzm=nzm,
        nzt=nzt,
        ngrdcol=ngrdcol,
        upwp_sfc=np.full(ngrdcol, 0.02, dtype=np.float64),
        vpwp_sfc=np.full(ngrdcol, 0.01, dtype=np.float64),
        ddzt_umvm_sqd=np.full((ngrdcol, nzm), 1.0e-6, dtype=np.float64),
        ice_supersat_frac=np.zeros((ngrdcol, nzt), dtype=np.float64),
        em=np.full((ngrdcol, nzm), 0.1, dtype=np.float64),
        sqrt_em_zt=np.full((ngrdcol, nzt), np.sqrt(0.1), dtype=np.float64),
        ufmin=0.1,
        tau_const=300.0,
        sfc_elevation=np.full(ngrdcol, -10.0, dtype=np.float64),
        lscale_max=np.full(ngrdcol, 1.0e5, dtype=np.float64),
        clubb_params=clubb_params,
        l_e3sm_config=False,
        l_smooth_heaviside_tau_wpxp=True,
        brunt_vaisala_freq_sqd_smth=np.full((ngrdcol, nzm), 1.0e-6, dtype=np.float64),
        ri_zm=np.full((ngrdcol, nzm), 0.1, dtype=np.float64),
        err_info=err_info,
    )

    assert len(out) == 20
    assert isinstance(out[0], ErrInfo)
    arrays = out[1:]
    nzt_idx = {0, 10, 13, 15, 16, 17, 18}
    for idx, arr in enumerate(arrays):
        expected_shape = (ngrdcol, nzt) if idx in nzt_idx else (ngrdcol, nzm)
        assert arr.shape == expected_shape
        assert np.all(np.isfinite(arr))

    err = clubb_api.get_err_code(ngrdcol)
    assert np.all(err == 0), f"diagnose_lscale_from_tau set err_code={err}"


def test_calc_lscale_directly_basic(gr):
    """calc_lscale_directly should return finite Lscale arrays with nzt shape."""
    ngrdcol, nzm, nzt = gr.ngrdcol, gr.nzm, gr.nzt
    clubb_params = clubb_api.init_clubb_params(ngrdcol, iunit=10, filename="")
    flags = clubb_api.get_default_config_flags()

    # Ensure the underlying stored PDF state is allocated, even though the
    # current Fortran branch for l_avg_Lscale=.false. does not consume it.
    pdf_params = clubb_api.init_pdf_params(nzt, ngrdcol)
    err_info = ErrInfo(ngrdcol=ngrdcol)

    err_info_out, lscale, lscale_up, lscale_down = clubb_api.calc_lscale_directly(
        gr,
        ngrdcol=ngrdcol,
        nzm=nzm,
        nzt=nzt,
        l_implemented=True,
        p_in_pa=np.full((ngrdcol, nzt), 100000.0, dtype=np.float64),
        exner=np.full((ngrdcol, nzt), 1.0, dtype=np.float64),
        rtm=np.full((ngrdcol, nzt), 0.01, dtype=np.float64),
        thlm=np.full((ngrdcol, nzt), 300.0, dtype=np.float64),
        thvm=np.full((ngrdcol, nzt), 300.0, dtype=np.float64),
        newmu=np.full(ngrdcol, 1.0e-3, dtype=np.float64),
        rtp2_zt=np.full((ngrdcol, nzt), 1.0e-6, dtype=np.float64),
        thlp2_zt=np.full((ngrdcol, nzt), 1.0e-4, dtype=np.float64),
        rtpthlp_zt=np.zeros((ngrdcol, nzt), dtype=np.float64),
        em=np.full((ngrdcol, nzm), 0.2, dtype=np.float64),
        thv_ds_zt=np.full((ngrdcol, nzt), 300.0, dtype=np.float64),
        lscale_max=np.full(ngrdcol, 1.0e5, dtype=np.float64),
        lmin=1.0,
        clubb_params=clubb_params,
        saturation_formula=int(flags.saturation_formula),
        l_lscale_plume_centered=False,
        pdf_params=pdf_params,
        err_info=err_info,
    )

    assert lscale.shape == (ngrdcol, nzt)
    assert lscale_up.shape == (ngrdcol, nzt)
    assert lscale_down.shape == (ngrdcol, nzt)
    assert np.all(np.isfinite(lscale))
    assert np.all(np.isfinite(lscale_up))
    assert np.all(np.isfinite(lscale_down))
    assert isinstance(err_info_out, ErrInfo)

    err = clubb_api.get_err_code(ngrdcol)
    assert np.all(err == 0), f"calc_lscale_directly set err_code={err}"
