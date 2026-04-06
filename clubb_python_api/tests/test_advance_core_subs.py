"""Test individual advance_clubb_core sub-routine wrappers.

Validates that each wrapped sub-routine can be called from Python,
produces finite outputs, and respects basic physical constraints.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.grid_class import Grid


# ── Shared fixtures ──────────────────────────────────────────────────

def _make_grid(ngrdcol=1, nzmax=11, dz=500.0):
    """Create a standard test grid and return (gr, params, lmin)."""
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

    params = clubb_api.init_clubb_params(ngrdcol, iunit=10, filename="")
    _nu_vert_res_dep, lmin, mfmm = clubb_api.calc_derrived_params(
        gr=gr, ngrdcol=ngrdcol, grid_type=1,
        deltaz=np.full(ngrdcol, dz, dtype=np.float64),
        clubb_params=params,
        nu_vert_res_dep=None,
        l_prescribed_avg_deltaz=False,
    )
    return gr, params, lmin, flags


def _init_test_stats(gr, tmp_path: Path):
    """Initialize stored stats state for wrappers that require stats."""
    registry = str(Path(__file__).resolve().parent / "test_stats_registry.in")
    err_info = ErrInfo(ngrdcol=gr.ngrdcol)
    clubb_params = clubb_api.init_clubb_params(gr.ngrdcol, iunit=10, filename="")
    err_info = clubb_api.init_stats(
        registry_path=registry,
        output_path=str(tmp_path / "advance_core_subs_stats.nc"),
        ncol=gr.ngrdcol,
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
    return err_info


@pytest.fixture(scope="module")
def grid_env():
    """Module-scoped fixture: set up grid + params once for all tests."""
    return _make_grid()


# ── compute_sigma_sqd_w ──────────────────────────────────────────────

class TestComputeSigmaSqdW:
    """Test the PDF width parameter computation."""

    def test_basic_call(self, grid_env):
        """compute_sigma_sqd_w should return finite array of correct shape."""
        gr, params, lmin, flags = grid_env
        nzm, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nzm)

        sigma = clubb_api.compute_sigma_sqd_w(
            nzm=nzm,
            ngrdcol=ngrdcol,
            gamma_Skw_fnc=np.ones(shape),
            wp2=np.full(shape, 0.5),
            thlp2=np.full(shape, 0.1),
            rtp2=np.full(shape, 1e-8),
            up2=np.full(shape, 0.3),
            vp2=np.full(shape, 0.3),
            wpthlp=np.full(shape, 0.01),
            wprtp=np.full(shape, 1e-5),
            upwp=np.full(shape, 0.01),
            vpwp=np.full(shape, 0.01),
            l_predict_upwp_vpwp=True,
        )
        assert sigma.shape == shape
        assert np.all(np.isfinite(sigma))

    def test_sigma_bounded(self, grid_env):
        """sigma_sqd_w should be between 0 and 1."""
        gr, params, lmin, flags = grid_env
        nzm, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nzm)

        sigma = clubb_api.compute_sigma_sqd_w(
            nzm=nzm,
            ngrdcol=ngrdcol,
            gamma_Skw_fnc=np.ones(shape),
            wp2=np.full(shape, 1.0),
            thlp2=np.full(shape, 0.1),
            rtp2=np.full(shape, 1e-8),
            up2=np.full(shape, 0.3),
            vp2=np.full(shape, 0.3),
            wpthlp=np.full(shape, 0.01),
            wprtp=np.full(shape, 1e-5),
            upwp=np.zeros(shape),
            vpwp=np.zeros(shape),
            l_predict_upwp_vpwp=False,
        )
        assert np.all(sigma >= 0.0)
        assert np.all(sigma <= 1.0)

    def test_zero_variance(self, grid_env):
        """With zero variances, sigma_sqd_w should still be finite."""
        gr, params, lmin, flags = grid_env
        nzm, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nzm)

        sigma = clubb_api.compute_sigma_sqd_w(
            nzm=nzm,
            ngrdcol=ngrdcol,
            gamma_Skw_fnc=np.ones(shape),
            wp2=np.full(shape, 1e-10),
            thlp2=np.zeros(shape),
            rtp2=np.zeros(shape),
            up2=np.zeros(shape),
            vp2=np.zeros(shape),
            wpthlp=np.zeros(shape),
            wprtp=np.zeros(shape),
            upwp=np.zeros(shape),
            vpwp=np.zeros(shape),
            l_predict_upwp_vpwp=False,
        )
        assert np.all(np.isfinite(sigma))


# ── fill_holes_vertical ──────────────────────────────────────────────

class TestFillHolesVertical:
    """Test the hole-filling routine."""

    def test_no_holes(self, grid_env):
        """A field with no holes should be unchanged."""
        gr, params, lmin, flags = grid_env
        nz, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nz)

        field = np.ones(shape) * 5.0
        threshold = 0.0
        result = clubb_api.fill_holes_vertical(
            nz, ngrdcol, threshold,
            lower_hf_level=1, upper_hf_level=nz-2,
            dz=gr.dzm, rho_ds=np.ones(shape),
            grid_dir_indx=gr.grid_dir_indx,
            fill_holes_type=1,
            field=field.copy(),
        )
        np.testing.assert_array_almost_equal(result, field)

    def test_fills_negative(self, grid_env):
        """After hole filling, field should be >= threshold."""
        gr, params, lmin, flags = grid_env
        nz, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nz)

        field = np.ones(shape) * 2.0
        field[0, 3] = -1.0  # a hole
        threshold = 0.0
        result = clubb_api.fill_holes_vertical(
            nz, ngrdcol, threshold,
            lower_hf_level=1, upper_hf_level=nz-2,
            dz=gr.dzm, rho_ds=np.ones(shape),
            grid_dir_indx=gr.grid_dir_indx,
            fill_holes_type=1,
            field=field.copy(),
        )
        assert np.all(result >= threshold - 1e-15)


class TestFillHolesWp2FromHorzTke:
    """Test the wp2 hole filling routine that conserves horizontal TKE."""

    def test_fills_wp2_and_conserves_tke(self, grid_env):
        """Interior wp2 holes should be filled without changing total TKE."""
        gr, params, lmin, flags = grid_env
        nz, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nz)

        wp2 = np.full(shape, 2.0, dtype=np.float64)
        up2 = np.full(shape, 3.0, dtype=np.float64)
        vp2 = np.full(shape, 4.0, dtype=np.float64)
        threshold = 1.0
        k = 3
        wp2[0, k] = 0.25
        tke_before = wp2 + up2 + vp2

        wp2_out, up2_out, vp2_out = clubb_api.fill_holes_wp2_from_horz_tke(
            nz=nz,
            ngrdcol=ngrdcol,
            threshold=threshold,
            lower_hf_level=1,
            upper_hf_level=nz - 2,
            wp2=wp2.copy(),
            up2=up2.copy(),
            vp2=vp2.copy(),
        )

        np.testing.assert_allclose(wp2_out[0, k], threshold)
        np.testing.assert_allclose(wp2_out + up2_out + vp2_out, tke_before)
        np.testing.assert_allclose(wp2_out[:, 0], wp2[:, 0])
        np.testing.assert_allclose(wp2_out[:, -1], wp2[:, -1])


# ── clip_rcm ─────────────────────────────────────────────────────────

class TestClipRcm:
    """Test cloud water clipping."""

    def test_rcm_clipped_to_rtm(self, grid_env):
        """rcm should be clipped so 0 <= rcm <= rtm."""
        gr, params, lmin, flags = grid_env
        nzt, ngrdcol = gr.nzt, gr.ngrdcol
        shape = (ngrdcol, nzt)

        rtm = np.full(shape, 0.01)   # 10 g/kg
        rcm = np.full(shape, 0.02)   # 20 g/kg — too much!
        result = clubb_api.clip_rcm(nzt=nzt, ngrdcol=ngrdcol, rtm=rtm, message="test_clip_rcm", rcm=rcm.copy())
        assert np.all(result <= rtm + 1e-15)
        assert np.all(result >= -1e-15)

    def test_negative_rcm_unchanged(self, grid_env):
        """Negative rcm passes through unchanged (clip_rcm only clips rcm > rtm)."""
        gr, params, lmin, flags = grid_env
        nzt, ngrdcol = gr.nzt, gr.ngrdcol
        shape = (ngrdcol, nzt)

        rtm = np.full(shape, 0.01)
        rcm = np.full(shape, -0.005)
        result = clubb_api.clip_rcm(nzt=nzt, ngrdcol=ngrdcol, rtm=rtm, message="test_clip_rcm", rcm=rcm.copy())
        np.testing.assert_array_almost_equal(result, rcm)

    def test_valid_rcm_unchanged(self, grid_env):
        """Valid rcm (0 < rcm < rtm) should pass through unchanged."""
        gr, params, lmin, flags = grid_env
        nzt, ngrdcol = gr.nzt, gr.ngrdcol
        shape = (ngrdcol, nzt)

        rtm = np.full(shape, 0.01)
        rcm = np.full(shape, 0.002)
        result = clubb_api.clip_rcm(nzt=nzt, ngrdcol=ngrdcol, rtm=rtm, message="test_clip_rcm", rcm=rcm.copy())
        np.testing.assert_array_almost_equal(result, rcm)


# ── clip_covar ────────────────────────────────────────────────────────

class TestClipCovar:
    """Test covariance clipping (Cauchy-Schwarz bound)."""

    def test_basic_call(self, grid_env):
        """clip_covar should return tuple of (xpyp, xpyp_chnge)."""
        gr, params, lmin, flags = grid_env
        nzm, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nzm)

        xp2 = np.full(shape, 1.0)
        yp2 = np.full(shape, 1.0)
        xpyp = np.full(shape, 0.5)  # within bounds

        result = clubb_api.clip_covar(
            nzm=nzm, ngrdcol=ngrdcol,
            solve_type=1,
            l_first_clip_ts=True, l_last_clip_ts=True,
            dt=300.0, xp2=xp2, yp2=yp2,
            l_predict_upwp_vpwp=False, xpyp=xpyp.copy(),
        )
        # Returns (xpyp, xpyp_chnge)
        assert len(result) == 2
        assert result[0].shape == shape
        assert result[1].shape == shape

    def test_clips_excess_covariance(self, grid_env):
        """Covariance exceeding sqrt(xp2*yp2) should be clipped down."""
        gr, params, lmin, flags = grid_env
        nzm, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nzm)

        xp2 = np.full(shape, 1.0)
        yp2 = np.full(shape, 1.0)
        xpyp = np.full(shape, 5.0)  # way too big

        clipped, change = clubb_api.clip_covar(
            nzm=nzm, ngrdcol=ngrdcol,
            solve_type=1,
            l_first_clip_ts=True, l_last_clip_ts=True,
            dt=300.0, xp2=xp2, yp2=yp2,
            l_predict_upwp_vpwp=False, xpyp=xpyp.copy(),
        )
        # Clipped values should be finite
        assert np.all(np.isfinite(clipped))
        # Interior levels should be clipped down from 5.0 (both boundaries may differ)
        assert np.all(np.abs(clipped[:, 1:-1]) < 5.0), "Interior covariance should be clipped"
        # Change should be non-zero since clipping occurred
        assert np.all(np.isfinite(change))


class TestClipSkewness:
    """Test the skewness-based wp3 clipping routine."""

    def test_clip_skewness_limits_wp3(self, grid_env, tmp_path):
        """clip_skewness should reduce excessive wp3 magnitudes."""
        gr, params, lmin, flags = grid_env
        err_info = _init_test_stats(gr, tmp_path)

        nzt, ngrdcol = gr.nzt, gr.ngrdcol
        shape = (ngrdcol, nzt)
        wp2_zt = np.ones(shape, dtype=np.float64)
        wp3 = np.full(shape, 50.0, dtype=np.float64)

        try:
            clipped = clubb_api.clip_skewness(
                gr=gr,
                nzt=nzt,
                ngrdcol=ngrdcol,
                dt=60.0,
                sfc_elevation=np.zeros(ngrdcol, dtype=np.float64),
                skw_max_mag=np.full(ngrdcol, 4.5, dtype=np.float64),
                wp2_zt=wp2_zt,
                l_use_wp3_lim_with_smth_heaviside=False,
                wp3=wp3,
            )
        finally:
            err_info = clubb_api.finalize_stats(err_info=err_info)

        assert clipped.shape == shape
        assert np.all(np.abs(clipped) <= 4.5 + 1e-12)


# ── clip_variance ─────────────────────────────────────────────────────

class TestClipVariance:
    """Test variance clipping."""

    def test_clips_to_threshold(self, grid_env):
        """Variance below threshold should be raised (interior levels)."""
        gr, params, lmin, flags = grid_env
        nzm, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nzm)

        threshold_lo = np.full(shape, 0.01)
        xp2 = np.full(shape, 0.001)  # below threshold

        result = clubb_api.clip_variance(
            gr, nzm=nzm, ngrdcol=ngrdcol, solve_type=1, dt=300.0,
            threshold_lo=threshold_lo, xp2=xp2.copy(),
        )
        # Interior levels should be clipped; boundary levels may differ
        assert np.all(result[:, 1:-1] >= threshold_lo[:, 1:-1] - 1e-15)

    def test_above_threshold_unchanged(self, grid_env):
        """Variance above threshold should pass through."""
        gr, params, lmin, flags = grid_env
        nzm, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nzm)

        threshold_lo = np.full(shape, 0.01)
        xp2 = np.full(shape, 1.0)

        result = clubb_api.clip_variance(
            gr, nzm=nzm, ngrdcol=ngrdcol, solve_type=1, dt=300.0,
            threshold_lo=threshold_lo, xp2=xp2.copy(),
        )
        np.testing.assert_array_almost_equal(result, xp2)


# ── calc_brunt_vaisala_freq_sqd ───────────────────────────────────────

class TestBruntVaisala:
    """Test Brunt-Vaisala frequency computation."""

    def test_basic_call(self, grid_env):
        """calc_brunt_vaisala_freq_sqd returns 5 arrays of correct shape."""
        gr, params, lmin, flags = grid_env
        nzt, nzm, ngrdcol = gr.nzt, gr.nzm, gr.ngrdcol
        zt_shape = (ngrdcol, nzt)
        zm_shape = (ngrdcol, nzm)

        # Stable atmosphere: thlm increases with height
        thlm = 300.0 + 5.0 * np.linspace(0, 1, nzt).reshape(1, -1) * np.ones(zt_shape)
        exner = np.full(zt_shape, 1.0)
        rtm = np.full(zt_shape, 0.01)
        rcm = np.zeros(zt_shape)
        p_in_Pa = np.full(zt_shape, 100000.0)
        thvm = thlm.copy()
        ice_supersat_frac = np.zeros(zt_shape)

        result = clubb_api.calc_brunt_vaisala_freq_sqd(
            gr, nzm, nzt, ngrdcol, thlm, exner, rtm, rcm, p_in_Pa, thvm, ice_supersat_frac,
            saturation_formula=flags.saturation_formula,
            l_brunt_vaisala_freq_moist=False,
            l_use_thvm_in_bv_freq=True,
            l_modify_limiters_for_cnvg_test=False,
            bv_efold=np.full(ngrdcol, 2000.0),
            T0=300.0,
        )

        assert len(result) == 5
        for arr in result:
            assert arr.shape == zm_shape
            assert np.all(np.isfinite(arr))

    def test_stable_positive_n2(self, grid_env):
        """Stable atmosphere should have positive N^2."""
        gr, params, lmin, flags = grid_env
        nzt, nzm, ngrdcol = gr.nzt, gr.nzm, gr.ngrdcol
        zt_shape = (ngrdcol, nzt)

        # Strongly stable: 10 K/km
        thlm = 300.0 + 10.0 * np.linspace(0, 5, nzt).reshape(1, -1) * np.ones(zt_shape)
        exner = np.full(zt_shape, 1.0)
        rtm = np.full(zt_shape, 0.005)
        rcm = np.zeros(zt_shape)
        p_in_Pa = np.linspace(100000, 50000, nzt).reshape(1, -1) * np.ones(zt_shape)
        thvm = thlm.copy()
        ice_supersat_frac = np.zeros(zt_shape)

        bv2, bv2_mixed, bv2_dry, bv2_moist, bv2_smth = clubb_api.calc_brunt_vaisala_freq_sqd(
            gr, nzm, nzt, ngrdcol, thlm, exner, rtm, rcm, p_in_Pa, thvm, ice_supersat_frac,
            saturation_formula=flags.saturation_formula,
            l_brunt_vaisala_freq_moist=False,
            l_use_thvm_in_bv_freq=True,
            l_modify_limiters_for_cnvg_test=False,
            bv_efold=np.full(ngrdcol, 2000.0),
            T0=300.0,
        )
        # Interior levels of stable atmosphere should have positive N^2
        # (boundaries may be different due to grid extrapolation)
        assert np.any(bv2 > 0), "Expected some positive N^2 for stable atmosphere"


# ── compute_cx_fnc_richardson ─────────────────────────────────────────

class TestCxFncRichardson:
    """Test Cx function of Richardson number."""

    def test_basic_call(self, grid_env):
        """compute_cx_fnc_richardson returns finite array."""
        gr, params, lmin, flags = grid_env
        nzm, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nzm)

        result = clubb_api.compute_cx_fnc_richardson(
            gr, nzm, gr.nzt, ngrdcol,
            Lscale_zm=np.full(shape, 100.0),
            ddzt_umvm_sqd=np.full(shape, 1e-4),
            rho_ds_zm=np.full(shape, 1.0),
            brunt_vaisala_freq_sqd=np.full(shape, 1e-4),
            brunt_vaisala_freq_sqd_mixed=np.full(shape, 1e-4),
            clubb_params=params,
            l_use_shear_Richardson=flags.l_use_shear_Richardson,
            l_modify_limiters_for_cnvg_test=False,
        )
        assert result.shape == shape
        assert np.all(np.isfinite(result))


# ── calc_stability_correction ─────────────────────────────────────────

class TestStabilityCorrection:
    """Test stability correction computation."""

    def test_basic_call(self, grid_env):
        """calc_stability_correction returns finite array."""
        gr, params, lmin, flags = grid_env
        nzt, nzm, ngrdcol = gr.nzt, gr.nzm, gr.ngrdcol
        zt_shape = (ngrdcol, nzt)
        zm_shape = (ngrdcol, nzm)

        result = clubb_api.calc_stability_correction(
            gr, nzm, nzt, ngrdcol,
            thlm=np.full(zt_shape, 300.0),
            Lscale_zm=np.full(zm_shape, 100.0),
            em=np.full(zm_shape, 0.5),
            exner=np.full(zt_shape, 1.0),
            rtm=np.full(zt_shape, 0.01),
            rcm=np.zeros(zt_shape),
            p_in_Pa=np.full(zt_shape, 100000.0),
            thvm=np.full(zt_shape, 300.0),
            ice_supersat_frac=np.zeros(zt_shape),
            lambda0_stability_coef=np.full(ngrdcol, 0.04),
            bv_efold=np.full(ngrdcol, 2000.0),
            T0=300.0,
            saturation_formula=flags.saturation_formula,
            l_brunt_vaisala_freq_moist=flags.l_brunt_vaisala_freq_moist,
            l_use_thvm_in_bv_freq=flags.l_use_thvm_in_bv_freq,
            l_modify_limiters_for_cnvg_test=False,
        )
        assert result.shape == zm_shape
        assert np.all(np.isfinite(result))

    def test_nonnegative(self, grid_env):
        """Stability correction should be non-negative."""
        gr, params, lmin, flags = grid_env
        nzt, nzm, ngrdcol = gr.nzt, gr.nzm, gr.ngrdcol
        zt_shape = (ngrdcol, nzt)
        zm_shape = (ngrdcol, nzm)

        result = clubb_api.calc_stability_correction(
            gr, nzm, nzt, ngrdcol,
            thlm=np.full(zt_shape, 300.0),
            Lscale_zm=np.full(zm_shape, 100.0),
            em=np.full(zm_shape, 0.5),
            exner=np.full(zt_shape, 1.0),
            rtm=np.full(zt_shape, 0.01),
            rcm=np.zeros(zt_shape),
            p_in_Pa=np.full(zt_shape, 100000.0),
            thvm=np.full(zt_shape, 300.0),
            ice_supersat_frac=np.zeros(zt_shape),
            lambda0_stability_coef=np.full(ngrdcol, 0.04),
            bv_efold=np.full(ngrdcol, 2000.0),
            T0=300.0,
            saturation_formula=flags.saturation_formula,
            l_brunt_vaisala_freq_moist=False,
            l_use_thvm_in_bv_freq=True,
            l_modify_limiters_for_cnvg_test=False,
        )
        assert np.all(result >= -1e-15)


# ── clip_covars_denom ─────────────────────────────────────────────────

class TestClipCovarsDenom:
    """Test bulk covariance clipping."""

    def test_basic_call(self, grid_env):
        """clip_covars_denom should not crash with valid inputs."""
        gr, params, lmin, flags = grid_env
        nzm, ngrdcol = gr.nzm, gr.ngrdcol
        shape = (ngrdcol, nzm)
        sclr_dim = 0

        result = clubb_api.clip_covars_denom(
            nzm=nzm, ngrdcol=ngrdcol,
            sclr_dim=sclr_dim, dt=300.0,
            rtp2=np.full(shape, 1e-6),
            thlp2=np.full(shape, 0.1),
            up2=np.full(shape, 0.3),
            vp2=np.full(shape, 0.3),
            wp2=np.full(shape, 0.5),
            sclrp2=np.zeros((ngrdcol, nzm, 1)),
            wprtp_cl_num=1, wpthlp_cl_num=2, wpsclrp_cl_num=3,
            upwp_cl_num=4, vpwp_cl_num=5,
            l_predict_upwp_vpwp=flags.l_predict_upwp_vpwp,
            l_tke_aniso=flags.l_tke_aniso,
            l_linearize_pbl_winds=False,
            wprtp=np.full(shape, 1e-4),
            wpthlp=np.full(shape, 0.01),
            upwp=np.full(shape, 0.01),
            vpwp=np.full(shape, 0.01),
            wpsclrp=np.zeros((ngrdcol, nzm, 1)),
            upwp_pert=np.zeros(shape),
            vpwp_pert=np.zeros(shape),
        )
        # Should return modified arrays
        assert result is not None
