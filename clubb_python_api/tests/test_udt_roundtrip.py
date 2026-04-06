"""Test push/pull roundtrips for all CLUBB UDT types.

Validates that Python NamedTuples can be pushed to Fortran module state
and pulled back with identical values.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python.derived_types import Grid, SclrIdx, NuVertResDep, pdf_parameter, implicit_coefs_terms
from clubb_python.clubb_api import init_pdf_params, init_pdf_params_zm, init_pdf_implicit
from clubb_python.derived_types import (
    set_fortran_grid, get_fortran_grid,
    set_fortran_sclr_idx, get_fortran_sclr_idx,
    set_fortran_nu_vert_res_dep, get_fortran_nu_vert_res_dep,
    get_fortran_pdf_params, set_fortran_pdf_params,
    get_fortran_pdf_params_zm, set_fortran_pdf_params_zm,
    get_fortran_implicit_coefs, set_fortran_implicit_coefs,
)
from clubb_python.derived_types.pdf_params import (
    pack_pdf_params,
    unpack_pdf_params,
    pack_implicit_coefs_2d,
    pack_implicit_coefs_3d,
    unpack_implicit_coefs,
)


def _make_test_grid(ngrdcol=1, nzm=11, nzt=10):
    """Construct a simple test grid with known values."""
    zm_1d = np.linspace(0.0, 1000.0, nzm)
    zt_1d = np.linspace(50.0, 950.0, nzt)

    zm = np.tile(zm_1d, (ngrdcol, 1)).astype(np.float64)
    zt = np.tile(zt_1d, (ngrdcol, 1)).astype(np.float64)

    dzm = np.full((ngrdcol, nzm), 100.0, dtype=np.float64)
    dzt = np.full((ngrdcol, nzt), 100.0, dtype=np.float64)
    invrs_dzm = np.full((ngrdcol, nzm), 0.01, dtype=np.float64)
    invrs_dzt = np.full((ngrdcol, nzt), 0.01, dtype=np.float64)

    weights_zt2zm = np.full((ngrdcol, nzm, 2), 0.5, dtype=np.float64)
    weights_zm2zt = np.full((ngrdcol, nzt, 2), 0.5, dtype=np.float64)

    return Grid(
        nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
        zm=zm, zt=zt, dzm=dzm, dzt=dzt,
        invrs_dzm=invrs_dzm, invrs_dzt=invrs_dzt,
        weights_zt2zm=weights_zt2zm, weights_zm2zt=weights_zm2zt,
        k_lb_zm=0, k_ub_zm=nzm - 1, k_lb_zt=0, k_ub_zt=nzt - 1,
        grid_dir_indx=1, grid_dir=1.0,
    )


class TestGrid:
    """Test Grid push/pull roundtrip."""

    def test_single_column_roundtrip(self):
        """Single-column grid should roundtrip exactly."""
        gr = _make_test_grid(ngrdcol=1, nzm=11, nzt=10)
        set_fortran_grid(gr)
        got = get_fortran_grid()

        for field in (
            "zm", "zt", "dzm", "dzt", "invrs_dzm", "invrs_dzt",
            "weights_zt2zm", "weights_zm2zt",
        ):
            np.testing.assert_array_equal(getattr(got, field), getattr(gr, field))

        assert got.k_lb_zm == gr.k_lb_zm
        assert got.k_ub_zm == gr.k_ub_zm
        assert got.k_lb_zt == gr.k_lb_zt
        assert got.k_ub_zt == gr.k_ub_zt
        assert got.grid_dir_indx == gr.grid_dir_indx
        assert got.grid_dir == gr.grid_dir

    def test_multi_column_roundtrip(self):
        """Multi-column grid should roundtrip exactly."""
        ngrdcol, nzm, nzt = 3, 6, 5
        rng = np.random.default_rng(123)

        gr = Grid(
            nzm=nzm, nzt=nzt, ngrdcol=ngrdcol,
            zm=rng.uniform(0, 2000, (ngrdcol, nzm)),
            zt=rng.uniform(0, 2000, (ngrdcol, nzt)),
            dzm=rng.uniform(50, 200, (ngrdcol, nzm)),
            dzt=rng.uniform(50, 200, (ngrdcol, nzt)),
            invrs_dzm=np.zeros((ngrdcol, nzm), dtype=np.float64),
            invrs_dzt=np.zeros((ngrdcol, nzt), dtype=np.float64),
            weights_zt2zm=rng.uniform(0, 1, (ngrdcol, nzm, 2)),
            weights_zm2zt=rng.uniform(0, 1, (ngrdcol, nzt, 2)),
            k_lb_zm=1, k_ub_zm=4, k_lb_zt=0, k_ub_zt=3,
            grid_dir_indx=-1, grid_dir=-1.0,
        )
        gr = gr._replace(invrs_dzm=1.0 / gr.dzm, invrs_dzt=1.0 / gr.dzt)

        set_fortran_grid(gr)
        got = get_fortran_grid()

        for field in (
            "zm", "zt", "dzm", "dzt", "invrs_dzm", "invrs_dzt",
            "weights_zt2zm", "weights_zm2zt",
        ):
            np.testing.assert_array_equal(getattr(got, field), getattr(gr, field))

        assert got.k_lb_zm == 1
        assert got.k_ub_zm == 4
        assert got.grid_dir_indx == -1
        assert got.grid_dir == -1.0

    def test_overwrite(self):
        """Second Grid push should replace the first."""
        set_fortran_grid(_make_test_grid(ngrdcol=1, nzm=5, nzt=4))
        gr = _make_test_grid(ngrdcol=2, nzm=8, nzt=7)
        set_fortran_grid(gr)

        got = get_fortran_grid()
        np.testing.assert_array_equal(got.zm, gr.zm)
        np.testing.assert_array_equal(got.zt, gr.zt)
        assert got.k_lb_zm == gr.k_lb_zm


class TestSclrIdx:
    """Test SclrIdx push/pull roundtrip."""

    def test_roundtrip(self):
        """Push arbitrary values, pull, compare."""
        idx = SclrIdx(iisclr_rt=1, iisclr_thl=2, iisclr_CO2=3,
                      iiedsclr_rt=4, iiedsclr_thl=5, iiedsclr_CO2=6)
        set_fortran_sclr_idx(idx)
        got = get_fortran_sclr_idx()
        assert got == idx

    def test_zero_values(self):
        """Push zero values roundtrip."""
        idx = SclrIdx(0, 0, 0, 0, 0, 0)
        set_fortran_sclr_idx(idx)
        got = get_fortran_sclr_idx()
        assert got == idx

    def test_overwrite(self):
        """Second SclrIdx push should replace the first."""
        set_fortran_sclr_idx(SclrIdx(10, 20, 30, 40, 50, 60))
        idx = SclrIdx(1, 2, 3, 4, 5, 6)
        set_fortran_sclr_idx(idx)
        got = get_fortran_sclr_idx()
        assert got == idx


class TestNuVertResDep:
    """Test NuVertResDep push/pull roundtrip."""

    def test_roundtrip(self):
        """Push random arrays, pull, compare."""
        nzm = 11
        rng = np.random.default_rng(42)
        nu = NuVertResDep(
            nzm=nzm,
            nu1=rng.random(nzm),
            nu2=rng.random(nzm),
            nu6=rng.random(nzm),
            nu8=rng.random(nzm),
            nu9=rng.random(nzm),
            nu10=rng.random(nzm),
            nu_hm=rng.random(nzm),
        )
        set_fortran_nu_vert_res_dep(nu)
        got = get_fortran_nu_vert_res_dep(nzm)
        assert got.nzm == nzm
        for name in ['nu1', 'nu2', 'nu6', 'nu8', 'nu9', 'nu10', 'nu_hm']:
            np.testing.assert_array_equal(getattr(got, name), getattr(nu, name))

    def test_overwrite_different_size(self):
        """Second NuVertResDep push should replace the first allocation."""
        set_fortran_nu_vert_res_dep(
            NuVertResDep(
                nzm=5,
                nu1=np.ones(5, dtype=np.float64),
                nu2=np.ones(5, dtype=np.float64),
                nu6=np.ones(5, dtype=np.float64),
                nu8=np.ones(5, dtype=np.float64),
                nu9=np.ones(5, dtype=np.float64),
                nu10=np.ones(5, dtype=np.float64),
                nu_hm=np.ones(5, dtype=np.float64),
            )
        )
        nu = NuVertResDep(
            nzm=20,
            nu1=np.full(20, 99.0, dtype=np.float64),
            nu2=np.full(20, 99.0, dtype=np.float64),
            nu6=np.full(20, 99.0, dtype=np.float64),
            nu8=np.full(20, 99.0, dtype=np.float64),
            nu9=np.full(20, 99.0, dtype=np.float64),
            nu10=np.full(20, 99.0, dtype=np.float64),
            nu_hm=np.full(20, 99.0, dtype=np.float64),
        )
        set_fortran_nu_vert_res_dep(nu)
        got = get_fortran_nu_vert_res_dep(20)
        np.testing.assert_array_equal(got.nu1, nu.nu1)


class TestPdfParams:
    """Test pdf_parameter packed push/pull roundtrip."""

    def test_init_pull_zeros(self):
        """After init, all pdf_params should be zero."""
        ngrdcol, nz = 1, 10
        p = init_pdf_params(nz, ngrdcol)
        assert p.ngrdcol == ngrdcol
        assert p.nz == nz
        packed = pack_pdf_params(p)
        assert packed.shape == (ngrdcol, nz, 47)
        np.testing.assert_array_equal(packed, 0.0)

    def test_roundtrip(self):
        """Modify packed values, push, pull, verify."""
        ngrdcol, nz = 2, 8
        init_pdf_params(nz, ngrdcol)

        rng = np.random.default_rng(123)
        data = rng.random((ngrdcol, nz, 47))
        p = unpack_pdf_params(data)
        set_fortran_pdf_params(p)

        got = get_fortran_pdf_params()
        got_packed = pack_pdf_params(got)
        assert got_packed.shape == (ngrdcol, nz, 47)
        np.testing.assert_allclose(got_packed, data, atol=1e-14)

    def test_zm_roundtrip(self):
        """Same roundtrip for pdf_params_zm."""
        ngrdcol, nz = 1, 5
        init_pdf_params_zm(nz, ngrdcol)

        rng = np.random.default_rng(456)
        data = rng.random((ngrdcol, nz, 47))
        p = unpack_pdf_params(data)
        set_fortran_pdf_params_zm(p)

        got = get_fortran_pdf_params_zm()
        got_packed = pack_pdf_params(got)
        assert got_packed.shape == (ngrdcol, nz, 47)
        np.testing.assert_allclose(got_packed, data, atol=1e-14)


class TestImplicitCoefs:
    """Test implicit_coefs_terms packed push/pull roundtrip."""

    def test_init_pull_zeros_no_sclr(self):
        """After init with sclr_dim=0, 2D data should be zero, no 3D data."""
        ngrdcol, nz = 1, 10
        ic = init_pdf_implicit(nz, ngrdcol, sclr_dim=0)
        assert ic.ngrdcol == ngrdcol
        assert ic.nz == nz
        assert ic.sclr_dim == 0
        data_2d = pack_implicit_coefs_2d(ic)
        assert data_2d.shape == (ngrdcol, nz, 19)
        np.testing.assert_array_equal(data_2d, 0.0)
        assert pack_implicit_coefs_3d(ic) is None

    def test_2d_roundtrip(self):
        """Modify 2D packed values, push, pull, verify."""
        ngrdcol, nz = 2, 6
        ic0 = init_pdf_implicit(nz, ngrdcol, sclr_dim=0)

        rng = np.random.default_rng(789)
        data_2d = rng.random((ngrdcol, nz, 19))
        ic = unpack_implicit_coefs(data_2d=data_2d, data_3d=None)
        set_fortran_implicit_coefs(ic)

        got = get_fortran_implicit_coefs()
        got_2d = pack_implicit_coefs_2d(got)
        assert got_2d.shape == (ngrdcol, nz, 19)
        np.testing.assert_allclose(got_2d, data_2d, atol=1e-14)

    def test_3d_roundtrip(self):
        """Roundtrip with sclr_dim > 0 for 3D scalar fields."""
        ngrdcol, nz, sclr_dim = 1, 8, 2
        ic0 = init_pdf_implicit(nz, ngrdcol, sclr_dim=sclr_dim)

        rng = np.random.default_rng(101112)
        data_2d = rng.random((ngrdcol, nz, 19))
        data_3d = rng.random((ngrdcol, nz, sclr_dim, 8))
        ic = unpack_implicit_coefs(data_2d=data_2d, data_3d=data_3d)
        set_fortran_implicit_coefs(ic)

        got = get_fortran_implicit_coefs()
        got_2d = pack_implicit_coefs_2d(got)
        got_3d = pack_implicit_coefs_3d(got)
        assert got_2d.shape == (ngrdcol, nz, 19)
        assert got_3d is not None
        assert got_3d.shape == (ngrdcol, nz, sclr_dim, 8)
        np.testing.assert_allclose(got_2d, data_2d, atol=1e-14)
        np.testing.assert_allclose(got_3d, data_3d, atol=1e-14)
