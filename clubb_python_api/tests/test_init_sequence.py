"""Test the CLUBB initialization sequence through the Python interface.

Validates that each init routine can be called without errors and that
the full init sequence produces a physically reasonable grid and parameter set.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python.derived_types.config_flags import ConfigFlags, CONFIG_FLAG_NAMES
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.pdf_params import implicit_coefs_terms
from clubb_python.derived_types.pdf_params import pdf_parameter
from clubb_python.derived_types.pdf_params import (
    pack_pdf_params,
    pack_implicit_coefs_2d,
    pack_implicit_coefs_3d,
)
from clubb_python.clubb_api import (
    get_default_config_flags, init_config_flags,
    init_err_info, get_err_code,
    setup_grid,
    init_clubb_params, get_param_names, calc_derrived_params,
    init_pdf_params, init_pdf_params_zm, init_pdf_implicit, init_pdf_implicit_coefs_terms,
)


class TestConfigFlags:
    """Test config flag get/set round-trip."""

    def test_default_flags_returns_namedtuple(self):
        """get_default_config_flags returns a ConfigFlags with 67 fields."""
        flags = get_default_config_flags()
        assert isinstance(flags, ConfigFlags)
        assert len(flags) == 67

    def test_default_integer_flags_positive(self):
        """Integer enum flags should be positive (valid enum values)."""
        flags = get_default_config_flags()
        assert flags.iiPDF_type > 0
        assert flags.saturation_formula > 0

    def test_default_logical_flags_are_bool(self):
        """Logical flags should be Python bools."""
        flags = get_default_config_flags()
        assert isinstance(flags.l_use_precip_frac, bool)
        assert isinstance(flags.l_tke_aniso, bool)

    def test_known_defaults(self):
        """Some known defaults from model_flags.F90."""
        flags = get_default_config_flags()
        assert flags.l_use_precip_frac is True
        assert flags.l_predict_upwp_vpwp is True
        assert flags.l_ho_nontrad_coriolis is False
        assert flags.l_tke_aniso is True
        assert flags.l_damp_wp2_using_em is True
        assert flags.l_mono_flux_lim_thlm is True
        assert flags.l_add_dycore_grid is False

    def test_init_config_flags_roundtrip(self):
        """init_config_flags stores flags that were obtained from defaults."""
        flags = get_default_config_flags()
        # Should not raise
        init_config_flags(flags)

    def test_modified_flags(self):
        """Can modify a flag and store it."""
        flags = get_default_config_flags()
        modified = flags._replace(l_tke_aniso=False)
        assert modified.l_tke_aniso is False
        init_config_flags(modified)


class TestErrInfo:
    """Test error info initialization."""

    def test_init_err_info(self):
        """init_err_info should set all error codes to 0 (no error)."""
        ngrdcol = 3
        init_err_info(ngrdcol)
        err = get_err_code(ngrdcol)
        assert err.shape == (ngrdcol,)
        assert np.all(err == 0), f"Expected all zeros, got {err}"


class TestSetupGrid:
    """Test grid setup through the Fortran setup_grid routine."""

    def test_evenly_spaced_ascending(self):
        """Setup an evenly-spaced ascending grid."""
        ngrdcol = 1
        nzmax = 11  # momentum levels
        dz = 500.0

        init_err_info(ngrdcol)

        # Momentum levels: 0, 500, 1000, ..., 5000
        zm = np.arange(nzmax, dtype=np.float64) * dz
        # Thermodynamic levels: midpoints
        zt = 0.5 * (zm[:-1] + zm[1:])

        momentum_heights = zm.reshape(1, -1)
        thermodynamic_heights = zt.reshape(1, -1)

        gr, err_info = setup_grid(
            nzmax=nzmax, ngrdcol=ngrdcol,
            sfc_elevation=np.zeros(ngrdcol, dtype=np.float64),
            l_implemented=True, l_ascending_grid=True,
            grid_type=1,  # evenly spaced
            deltaz=np.full(ngrdcol, dz, dtype=np.float64),
            zm_init=np.zeros(ngrdcol, dtype=np.float64),
            zm_top=np.full(ngrdcol, zm[-1], dtype=np.float64),
            momentum_heights=momentum_heights,
            thermodynamic_heights=thermodynamic_heights,
        err_info=ErrInfo(ngrdcol=ngrdcol),
        )

        # No errors
        err = get_err_code(ngrdcol)
        assert np.all(err == 0), f"Grid setup failed with err_code={err}"

        # setup_grid now returns Grid directly
        assert gr.nzm == nzmax
        assert gr.nzt == nzmax - 1
        assert gr.ngrdcol == ngrdcol

        # Grid should be ascending
        assert gr.zm[0, -1] > gr.zm[0, 0]

    def test_multi_column_grid(self):
        """Setup grid with 2 columns."""
        ngrdcol = 2
        nzmax = 6
        dz = 1000.0

        init_err_info(ngrdcol)

        zm = np.arange(nzmax, dtype=np.float64) * dz
        zt = 0.5 * (zm[:-1] + zm[1:])

        momentum_heights = np.tile(zm, (ngrdcol, 1))
        thermodynamic_heights = np.tile(zt, (ngrdcol, 1))

        gr, err_info = setup_grid(
            nzmax=nzmax, ngrdcol=ngrdcol,
            sfc_elevation=np.zeros(ngrdcol, dtype=np.float64),
            l_implemented=True, l_ascending_grid=True,
            grid_type=1,
            deltaz=np.full(ngrdcol, dz, dtype=np.float64),
            zm_init=np.zeros(ngrdcol, dtype=np.float64),
            zm_top=np.full(ngrdcol, zm[-1], dtype=np.float64),
            momentum_heights=momentum_heights,
            thermodynamic_heights=thermodynamic_heights,
        err_info=ErrInfo(ngrdcol=ngrdcol),
        )

        err = get_err_code(ngrdcol)
        assert np.all(err == 0)
        assert gr.ngrdcol == ngrdcol


class TestClubbParams:
    """Test parameter initialization."""

    def test_default_params_shape(self):
        """init_clubb_params returns (ngrdcol, nparams) array."""
        ngrdcol = 1
        params = init_clubb_params(ngrdcol, iunit=10, filename="")
        assert params.shape == (ngrdcol, len(get_param_names()))

    def test_default_params_finite(self):
        """All default parameters should be finite."""
        params = init_clubb_params(1, iunit=10, filename="")
        assert np.all(np.isfinite(params))

    def test_multi_column_params(self):
        """Each column gets the same defaults."""
        ngrdcol = 3
        params = init_clubb_params(ngrdcol, iunit=10, filename="")
        assert params.shape == (ngrdcol, len(get_param_names()))
        # All columns should be identical for defaults
        np.testing.assert_array_equal(params[0], params[1])
        np.testing.assert_array_equal(params[0], params[2])


class TestCalcDerrivedParams:
    """Test derived parameter calculation."""

    def test_produces_finite_outputs(self):
        """calc_derrived_params should return finite lmin and mixt_frac_max_mag."""
        ngrdcol = 1
        nzmax = 11
        dz = 500.0

        # Full init sequence up to calc_derrived_params
        init_err_info(ngrdcol)

        zm = np.arange(nzmax, dtype=np.float64) * dz
        zt = 0.5 * (zm[:-1] + zm[1:])

        gr, err_info = setup_grid(
            nzmax=nzmax, ngrdcol=ngrdcol,
            sfc_elevation=np.zeros(ngrdcol, dtype=np.float64),
            l_implemented=True, l_ascending_grid=True,
            grid_type=1,
            deltaz=np.full(ngrdcol, dz, dtype=np.float64),
            zm_init=np.zeros(ngrdcol, dtype=np.float64),
            zm_top=np.full(ngrdcol, zm[-1], dtype=np.float64),
            momentum_heights=zm.reshape(1, -1),
            thermodynamic_heights=zt.reshape(1, -1),
        err_info=ErrInfo(ngrdcol=ngrdcol),
        )

        params = init_clubb_params(ngrdcol, iunit=10, filename="")
        deltaz_arr = np.full(ngrdcol, dz, dtype=np.float64)
        nu_vert_res_dep, lmin, mixt_frac_max_mag = calc_derrived_params(
            gr=gr,
            ngrdcol=ngrdcol,
            grid_type=1, deltaz=deltaz_arr,
            clubb_params=params, nu_vert_res_dep=None, l_prescribed_avg_deltaz=False,
        )

        assert np.isfinite(lmin), f"lmin not finite: {lmin}"
        assert np.isfinite(mixt_frac_max_mag), f"mixt_frac_max_mag not finite: {mixt_frac_max_mag}"
        assert nu_vert_res_dep.nzm == gr.ngrdcol
        assert lmin > 0, f"lmin should be positive, got {lmin}"


class TestPdfInit:
    """Test PDF parameter initialization."""

    def test_init_pdf_params(self):
        """init_pdf_params should return pulled pdf_parameter state."""
        p = init_pdf_params(10, 1)
        assert isinstance(p, pdf_parameter)
        assert p.ngrdcol == 1
        assert p.nz == 10
        packed = pack_pdf_params(p)
        assert packed.shape == (1, 10, 47)
        np.testing.assert_allclose(packed, 0.0)

    def test_init_pdf_params_zm(self):
        """init_pdf_params_zm should return pulled pdf_parameter state."""
        p = init_pdf_params_zm(10, 1)
        assert isinstance(p, pdf_parameter)
        assert p.ngrdcol == 1
        assert p.nz == 10
        packed = pack_pdf_params(p)
        assert packed.shape == (1, 10, 47)
        np.testing.assert_allclose(packed, 0.0)

    def test_init_pdf_implicit(self):
        """init_pdf_implicit should return pulled implicit-coefs state."""
        ic = init_pdf_implicit(10, 1, 0)
        assert isinstance(ic, implicit_coefs_terms)
        assert ic.ngrdcol == 1
        assert ic.nz == 10
        assert ic.sclr_dim == 0
        data_2d = pack_implicit_coefs_2d(ic)
        assert data_2d.shape == (1, 10, 19)
        np.testing.assert_allclose(data_2d, 0.0)
        assert pack_implicit_coefs_3d(ic) is None

    def test_init_pdf_implicit_with_scalars(self):
        """init_pdf_implicit with sclr_dim > 0 should return 3D packed state."""
        ic = init_pdf_implicit(10, 1, 2)
        assert isinstance(ic, implicit_coefs_terms)
        assert ic.ngrdcol == 1
        assert ic.nz == 10
        assert ic.sclr_dim == 2
        data_2d = pack_implicit_coefs_2d(ic)
        data_3d = pack_implicit_coefs_3d(ic)
        assert data_2d.shape == (1, 10, 19)
        assert data_3d is not None
        assert data_3d.shape == (1, 10, 2, 8)
        np.testing.assert_allclose(data_2d, 0.0)
        np.testing.assert_allclose(data_3d, 0.0)

    def test_init_pdf_implicit_coefs_terms_alias(self):
        """Alias should call the same Fortran-backed implicit-coefs initializer."""
        ic = init_pdf_implicit_coefs_terms(10, 1, 0)
        assert isinstance(ic, implicit_coefs_terms)
        assert ic.ngrdcol == 1
        assert ic.nz == 10
        assert ic.sclr_dim == 0


class TestFullInitSequence:
    """Test the complete initialization sequence end-to-end."""

    def test_full_init_no_errors(self):
        """Run the full init sequence and verify no fatal errors."""
        ngrdcol = 1
        nzmax = 11
        nzt = nzmax - 1
        dz = 500.0

        # 1. Init error info
        init_err_info(ngrdcol)

        # 2. Config flags (defaults)
        flags = get_default_config_flags()
        init_config_flags(flags)

        # 3. Setup grid (now returns Grid directly)
        zm = np.arange(nzmax, dtype=np.float64) * dz
        zt = 0.5 * (zm[:-1] + zm[1:])
        gr, err_info = setup_grid(
            nzmax=nzmax, ngrdcol=ngrdcol,
            sfc_elevation=np.zeros(ngrdcol, dtype=np.float64),
            l_implemented=True, l_ascending_grid=True,
            grid_type=1,
            deltaz=np.full(ngrdcol, dz, dtype=np.float64),
            zm_init=np.zeros(ngrdcol, dtype=np.float64),
            zm_top=np.full(ngrdcol, zm[-1], dtype=np.float64),
            momentum_heights=zm.reshape(1, -1),
            thermodynamic_heights=zt.reshape(1, -1),
        err_info=ErrInfo(ngrdcol=ngrdcol),
        )
        err = get_err_code(ngrdcol)
        assert np.all(err == 0), f"Grid setup error: {err}"

        # 4. Init params + derived params
        params = init_clubb_params(ngrdcol, iunit=10, filename="")
        nu_vert_res_dep, lmin, mfmm = calc_derrived_params(
            gr=gr, ngrdcol=ngrdcol, grid_type=1,
            deltaz=np.full(ngrdcol, dz, dtype=np.float64),
            clubb_params=params,
            nu_vert_res_dep=None,
            l_prescribed_avg_deltaz=False,
        )
        assert lmin > 0
        assert nu_vert_res_dep.nzm == gr.ngrdcol

        # 5. Init PDF structures
        init_pdf_params(nzt, ngrdcol)
        init_pdf_params_zm(nzt, ngrdcol)
        init_pdf_implicit(nzt, ngrdcol, sclr_dim=0)

        # 6. Final check: no errors
        err = get_err_code(ngrdcol)
        assert np.all(err == 0), f"Post-init error: {err}"

        # 7. Grid should be physically valid
        assert gr.zm[0, -1] > gr.zm[0, 0]
        assert np.all(gr.dzm > 0)
        assert np.all(gr.dzt > 0)
