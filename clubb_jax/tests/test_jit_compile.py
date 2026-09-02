"""Compile-check the currently jitted CLUBB-JAX helper routines.

A ``@jax.jit`` decorator only wraps a function.  The useful check is forcing a
representative call through ``lower().compile()``, which verifies that JAX can
stage and compile that exact static-argument branch.

Run directly:

  /path/to/.venv-jax/bin/python clubb_jax/tests/test_jit_compile.py

Or with pytest:

  /path/to/.venv-jax/bin/python -m pytest clubb_jax/tests/test_jit_compile.py
"""

from __future__ import annotations

import math
import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import jax

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core import clip_explicit, model_flags, parameter_indices
from clubb_jax.src.CLUBB_core.advance_helper_module import (
    Lscale_width_vert_avg,
    calc_Ri_zm,
    calc_brunt_vaisala_freq_sqd,
    calc_ddzt_umvm_sqd,
    calc_stability_correction,
    calc_wp3_on_wp2,
    calc_xpwp,
    calculate_thlp2_rad,
    compute_Cx_fnc_Richardson,
    smooth_heaviside_peskin,
    smooth_max,
    smooth_min,
    wp23_term_splat_lhs,
)
from clubb_jax.src.CLUBB_core.advance_xp2_xpyp_module import (
    advance_xp2_xpyp,
    calc_xp2_xpyp_ta_terms,
    pos_definite_variances,
    solve_xp2_xpyp_with_multiple_lhs,
    solve_xp2_xpyp_with_single_lhs,
    stats_finalize_xp2_xpyp_terms,
    term_dp1_lhs,
    term_dp1_rhs,
    term_pr1,
    term_pr2,
    term_tp_rhs,
    xp2_xpyp_lhs,
    xp2_xpyp_rhs,
    xp2_xpyp_single_lhs_valid,
    xp2_xpyp_solve,
    xp2_xpyp_uv_rhs,
    xp2_xpyp_rtp2,
    xp2_xpyp_thlp2,
    xp2_xpyp_rtpthlp,
    xp2_xpyp_up2,
    xp2_xpyp_up2_vp2,
    xp2_xpyp_single_lhs,
)
from clubb_jax.src.CLUBB_core.clip_explicit import (
    clip_covar,
    clip_rcm,
    clip_skewness_core,
)
from clubb_jax.src.CLUBB_core.diffusion import diffusion_zm_lhs, diffusion_zt_lhs
from clubb_jax.src.CLUBB_core.fill_holes import (
    fill_holes_global,
    fill_holes_sliding_window,
    fill_holes_vertical,
    fill_holes_wp2_from_horz_tke,
)
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zm_lhs, term_ma_zt_lhs
from clubb_jax.src.CLUBB_core.penta_lu_solver import (
    penta_lu_solve_multiple_rhs_lhs,
    penta_lu_solve_single_rhs_multiple_lhs,
)
from clubb_jax.src.CLUBB_core.tridiag_lu_solver import (
    tridiag_lu_solve_multiple_rhs_lhs,
    tridiag_lu_solve_single_rhs_lhs,
    tridiag_lu_solve_single_rhs_multiple_lhs,
)
from clubb_jax.src.CLUBB_core.turbulent_adv_pdf import (
    xpyp_term_ta_pdf_lhs,
    xpyp_term_ta_pdf_lhs_godunov,
    xpyp_term_ta_pdf_rhs,
    xpyp_term_ta_pdf_rhs_godunov,
)
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.grid_class import Grid
from clubb_jax.src.CLUBB_core.nu_vert_res_dep import NuVertResDep
from clubb_jax.src.CLUBB_core.pdf_params import implicit_coefs_terms
from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx


NGRDCOL = 2
NZM = 6
NZT = 5
NRHS = 3


def _compile(function, *args, **kwargs):
    """Force JAX compilation for one explicit call signature."""
    function.lower(*args, **kwargs).compile()


def _array(shape: tuple[int, ...], base: float = 1.0, step: float = 0.01):
    values = jnp.arange(math.prod(shape), dtype=jnp.float64).reshape(shape)
    return base + step * values


def _zt_array(base: float = 1.0):
    return _array((NGRDCOL, NZT), base)


def _zm_array(base: float = 1.0):
    return _array((NGRDCOL, NZM), base)


def _nu_array(base: float = 0.01):
    return _array((NGRDCOL,), base)


def _sample_grid() -> Grid:
    """Small shape-correct grid fixture for compile-only tests."""
    zm = jnp.broadcast_to(
        jnp.linspace(0.0, 500.0, NZM, dtype=jnp.float64),
        (NGRDCOL, NZM),
    )
    zt = jnp.broadcast_to(
        jnp.linspace(50.0, 450.0, NZT, dtype=jnp.float64),
        (NGRDCOL, NZT),
    )
    dzm = jnp.full((NGRDCOL, NZM), 100.0, dtype=jnp.float64)
    dzt = jnp.full((NGRDCOL, NZT), 100.0, dtype=jnp.float64)
    weights_zt2zm = jnp.full((NGRDCOL, NZM, 2), 0.5, dtype=jnp.float64)
    weights_zm2zt = jnp.full((NGRDCOL, NZT, 2), 0.5, dtype=jnp.float64)

    return Grid(
        nzm=NZM,
        nzt=NZT,
        ngrdcol=NGRDCOL,
        zm=zm,
        zt=zt,
        dzm=dzm,
        dzt=dzt,
        invrs_dzm=1.0 / dzm,
        invrs_dzt=1.0 / dzt,
        weights_zt2zm=weights_zt2zm,
        weights_zm2zt=weights_zm2zt,
        k_lb_zm=1,
        k_ub_zm=NZM,
        k_lb_zt=1,
        k_ub_zt=NZT,
        grid_dir_indx=1,
        grid_dir=1.0,
    )


def _clubb_params():
    params = jnp.ones((NGRDCOL, parameter_indices.nparams), dtype=jnp.float64)
    params = params.at[:, parameter_indices.iRichardson_num_min].set(-1.0)
    params = params.at[:, parameter_indices.iRichardson_num_max].set(1.0)
    params = params.at[:, parameter_indices.iCx_min].set(0.1)
    params = params.at[:, parameter_indices.iCx_max].set(1.0)
    params = params.at[:, parameter_indices.ithlp2_rad_coef].set(0.2)
    return params


def _tridiag_lhs_single():
    lhs = jnp.zeros((3, NZM), dtype=jnp.float64)
    lhs = lhs.at[0].set(0.1)
    lhs = lhs.at[1].set(4.0)
    lhs = lhs.at[2].set(0.1)
    return lhs


def _tridiag_lhs_multi():
    lhs = jnp.zeros((3, NGRDCOL, NZM), dtype=jnp.float64)
    lhs = lhs.at[0].set(0.1)
    lhs = lhs.at[1].set(4.0)
    lhs = lhs.at[2].set(0.1)
    return lhs


def _penta_lhs_multi():
    lhs = jnp.zeros((5, NGRDCOL, NZM), dtype=jnp.float64)
    lhs = lhs.at[0].set(0.01)
    lhs = lhs.at[1].set(0.1)
    lhs = lhs.at[2].set(4.0)
    lhs = lhs.at[3].set(0.1)
    lhs = lhs.at[4].set(0.01)
    return lhs


def _implicit_coefs_terms():
    scalar_zt = _zt_array(1.0)
    scalar_fields = [_array((NGRDCOL, NZT, 1), 1.0)] * 8
    return implicit_coefs_terms(
        NGRDCOL,
        NZT,
        1,
        *([scalar_zt] * 19),
        *scalar_fields,
    )


def _stats(l_sample: bool = False):
    names = (
        "rtp2_ta",
        "rtp2_tp",
        "rtp2_dp1",
        "rtp2_forcing",
        "thlp2_ta",
        "thlp2_tp",
        "thlp2_dp1",
        "thlp2_forcing",
        "rtpthlp_ta",
        "rtpthlp_tp1",
        "rtpthlp_tp2",
        "rtpthlp_dp1",
        "rtpthlp_forcing",
        "up2_ta",
        "up2_tp",
        "up2_dp1",
        "up2_pr1",
        "up2_pr2",
        "up2_splat",
        "vp2_ta",
        "vp2_tp",
        "vp2_dp1",
        "vp2_pr1",
        "vp2_pr2",
        "vp2_splat",
    )
    return JaxStats.empty(
        l_sample=l_sample,
        names=names,
        ncol=NGRDCOL,
        max_nlev=NZM,
    )


def _nu_vert_res_dep():
    return NuVertResDep(NZM, *([_nu_array(0.01)] * 7))


def _lhs3(base: float = 0.0):
    return _array((3, NGRDCOL, NZM), base)


def _advance_xp2_xpyp_args():
    gr = _sample_grid()
    sclr_dim = 1
    sclr_idx = SclrIdx(1, 0, 0, 0, 0, 0)
    zm = _zm_array(1.0)
    zt = _zt_array(1.0)
    sclr_zm = _array((NGRDCOL, NZM, sclr_dim), 1.0)
    sclr_zt = _array((NGRDCOL, NZT, sclr_dim), 1.0)

    return (
        NZM,
        NZT,
        NGRDCOL,
        sclr_dim,
        jnp.array([1.0e-4], dtype=jnp.float64),
        gr,
        sclr_idx,
        _zm_array(0.1),
        _zm_array(0.1),
        _zm_array(0.1),
        _zm_array(0.01),
        _zt_array(0.01),
        _zm_array(0.01),
        _zt_array(300.0),
        _zm_array(0.01),
        _zm_array(0.01),
        _zt_array(1.0),
        _zt_array(1.0),
        _zm_array(1.0),
        _zt_array(0.1),
        _zm_array(0.01),
        _zm_array(0.01),
        _zm_array(0.1),
        _zt_array(0.01),
        _zt_array(0.01),
        _zt_array(0.01),
        _zt_array(0.1),
        _zm_array(0.0),
        _zm_array(0.0),
        _zm_array(0.0),
        _zm_array(1.0),
        _zt_array(1.0),
        1.0 / _zm_array(1.0),
        _zm_array(300.0),
        _zt_array(0.0),
        _implicit_coefs_terms(),
        1.0,
        jnp.array([1.0e-4, 1.0e-4], dtype=jnp.float64),
        sclr_zt,
        sclr_zm,
        sclr_zt,
        sclr_zt,
        sclr_zt,
        _zm_array(0.01),
        _clubb_params(),
        _nu_vert_res_dep(),
        model_flags.iiPDF_ADG1,
        model_flags.tridiag_lu,
        0,
        False,
        False,
        False,
        False,
        False,
        False,
        True,
        _stats(l_sample=False),
        _zm_array(1.0),
        _zm_array(1.0),
        _zm_array(0.01),
        _zm_array(1.0),
        _zm_array(1.0),
        sclr_zm,
        sclr_zm,
        sclr_zm,
        ErrInfo.initialized(NGRDCOL),
    )


def test_derived_type_pytree_metadata_is_static():
    gr = _sample_grid()
    sclr_idx = SclrIdx(1, 2, 3, 4, 5, 6)
    nu = _nu_vert_res_dep()
    coefs = _implicit_coefs_terms()
    err_info = ErrInfo.initialized(NGRDCOL)

    assert len(jax.tree_util.tree_leaves(gr)) == 8
    assert len(jax.tree_util.tree_leaves(sclr_idx)) == 0
    assert len(jax.tree_util.tree_leaves(nu)) == 7
    assert len(jax.tree_util.tree_leaves(coefs)) == 27
    assert len(jax.tree_util.tree_leaves(err_info)) == 2


def test_xp2_xpyp_single_lhs_valid_compiles():
    _compile(
        xp2_xpyp_single_lhs_valid,
        _clubb_params(),
        model_flags.iiPDF_ADG1,
        False,
    )


def test_term_tp_rhs_compiles():
    _compile(
        term_tp_rhs,
        NZM,
        NZT,
        NGRDCOL,
        _zt_array(1.0),
        _zt_array(2.0),
        _zm_array(0.1),
        _zm_array(0.2),
        _sample_grid().invrs_dzm,
    )


def test_term_dp1_lhs_compiles():
    _compile(
        term_dp1_lhs,
        NZM,
        NGRDCOL,
        _sample_grid(),
        _zm_array(1.0),
        _zm_array(0.1),
    )


def test_term_dp1_rhs_compiles():
    _compile(term_dp1_rhs, NZM, NGRDCOL, _zm_array(1.0), _zm_array(0.1), 1.0e-8)


def test_term_pr1_compiles():
    _compile(
        term_pr1,
        NZM,
        NGRDCOL,
        jnp.array([1.0, 1.0], dtype=jnp.float64),
        jnp.array([1.0, 1.0], dtype=jnp.float64),
        _zm_array(1.0),
        _zm_array(1.0),
        _zm_array(0.1),
        _zm_array(0.1),
    )


def test_term_pr2_compiles():
    _compile(
        term_pr2,
        NZM,
        NZT,
        NGRDCOL,
        _sample_grid(),
        jnp.array([1.0, 1.0], dtype=jnp.float64),
        jnp.array([1.0, 1.0], dtype=jnp.float64),
        _zm_array(300.0),
        _zm_array(0.01),
        _zm_array(0.01),
        _zm_array(0.01),
        _zt_array(1.0),
        _zt_array(1.0),
    )


def test_xp2_xpyp_lhs_compiles():
    _compile(
        xp2_xpyp_lhs,
        NZM,
        NGRDCOL,
        1.0,
        _sample_grid(),
        _lhs3(0.0),
        _lhs3(0.0),
        _lhs3(0.0),
        _zm_array(0.1),
    )


def test_xp2_xpyp_solve_compiles():
    _compile(
        xp2_xpyp_solve,
        NZM,
        NGRDCOL,
        xp2_xpyp_single_lhs,
        True,
        _sample_grid(),
        model_flags.tridiag_lu,
        _stats(l_sample=False),
        _array((NGRDCOL, NZM, 3), 1.0),
        _tridiag_lhs_multi(),
        ErrInfo.initialized(NGRDCOL),
    )


def test_stats_finalize_xp2_xpyp_terms_compiles():
    _compile(
        stats_finalize_xp2_xpyp_terms,
        NZM,
        NGRDCOL,
        xp2_xpyp_rtp2,
        _sample_grid(),
        _zm_array(1.0),
        _zm_array(0.1),
        _zm_array(0.1),
        _lhs3(0.0),
        _lhs3(0.0),
        _lhs3(0.0),
        _stats(l_sample=False),
    )


def test_xp2_xpyp_uv_rhs_compiles():
    _compile(
        xp2_xpyp_uv_rhs,
        NZM,
        NZT,
        NGRDCOL,
        _sample_grid(),
        xp2_xpyp_up2,
        1.0,
        _zm_array(1.0),
        _zm_array(0.01),
        _zm_array(0.1),
        _zm_array(0.1),
        _zt_array(1.0),
        _zt_array(1.0),
        _zm_array(0.01),
        _zm_array(0.01),
        _zm_array(1.0),
        _zm_array(1.0),
        _zm_array(300.0),
        jnp.array([1.0, 1.0], dtype=jnp.float64),
        jnp.array([1.0, 1.0], dtype=jnp.float64),
        jnp.array([1.0, 1.0], dtype=jnp.float64),
        jnp.array([1.0, 1.0], dtype=jnp.float64),
        _zm_array(0.01),
        _lhs3(0.0),
        _zm_array(0.0),
        _zm_array(0.1),
        _zm_array(0.1),
        _stats(l_sample=False),
    )


def test_xp2_xpyp_rhs_compiles():
    _compile(
        xp2_xpyp_rhs,
        NZM,
        NZT,
        NGRDCOL,
        _sample_grid(),
        xp2_xpyp_rtp2,
        1.0,
        _zm_array(0.01),
        _zm_array(0.01),
        _zt_array(1.0),
        _zt_array(1.0),
        _zm_array(1.0),
        _zm_array(0.0),
        _zm_array(1.0),
        _zm_array(0.1),
        1.0e-8,
        _lhs3(0.0),
        _zm_array(0.0),
        _stats(l_sample=False),
    )


def test_calc_xp2_xpyp_ta_terms_compiles():
    _compile(
        calc_xp2_xpyp_ta_terms,
        NZM,
        NZT,
        NGRDCOL,
        1,
        _sample_grid(),
        _zm_array(0.01),
        _zt_array(0.01),
        _zm_array(0.01),
        _zt_array(0.01),
        _zt_array(0.01),
        _zm_array(1.0),
        _zm_array(1.0),
        _zm_array(0.01),
        _zm_array(0.01),
        _zm_array(0.01),
        _zm_array(1.0),
        _zm_array(1.0),
        _zm_array(1.0),
        _zt_array(1.0),
        _array((NGRDCOL, NZM, 1), 0.01),
        _array((NGRDCOL, NZT, 1), 0.01),
        _array((NGRDCOL, NZT, 1), 0.01),
        _array((NGRDCOL, NZT, 1), 0.01),
        _array((NGRDCOL, NZM, 1), 1.0),
        _array((NGRDCOL, NZM, 1), 0.01),
        _array((NGRDCOL, NZM, 1), 0.01),
        _zt_array(1.0),
        1.0 / _zm_array(1.0),
        _zm_array(1.0),
        _zm_array(0.1),
        _zt_array(0.1),
        _zm_array(0.1),
        _implicit_coefs_terms(),
        True,
        jnp.array([1.0, 1.0], dtype=jnp.float64),
        model_flags.iiPDF_ADG1,
        False,
        False,
        _stats(l_sample=False),
    )


def test_pos_definite_variances_compiles():
    _compile(
        pos_definite_variances,
        NZM,
        NGRDCOL,
        _sample_grid(),
        xp2_xpyp_rtp2,
        0,
        1.0,
        1.0e-8,
        _zm_array(1.0),
        _stats(l_sample=False),
        _zm_array(1.0),
    )


def test_solve_xp2_xpyp_with_single_lhs_compiles():
    _compile(
        solve_xp2_xpyp_with_single_lhs,
        NZM,
        NZT,
        NGRDCOL,
        1,
        jnp.array([1.0e-4], dtype=jnp.float64),
        _sample_grid(),
        SclrIdx(1, 0, 0, 0, 0, 0),
        _zm_array(1.0),
        _zm_array(0.1),
        _zt_array(0.01),
        _zt_array(300.0),
        _zm_array(0.01),
        _zm_array(0.01),
        _zm_array(0.0),
        _zm_array(0.0),
        _zm_array(0.0),
        _array((NGRDCOL, NZT, 1), 1.0),
        _array((NGRDCOL, NZM, 1), 0.01),
        _lhs3(0.0),
        _lhs3(0.0),
        _lhs3(0.0),
        _zm_array(0.0),
        _zm_array(0.0),
        _zm_array(0.0),
        _array((NGRDCOL, NZM, 1), 0.0),
        _array((NGRDCOL, NZM, 1), 0.0),
        _array((NGRDCOL, NZM, 1), 0.0),
        1.0,
        True,
        False,
        True,
        model_flags.tridiag_lu,
        _stats(l_sample=False),
        _zm_array(1.0),
        _zm_array(1.0),
        _zm_array(0.01),
        _array((NGRDCOL, NZM, 1), 1.0),
        _array((NGRDCOL, NZM, 1), 0.01),
        _array((NGRDCOL, NZM, 1), 0.01),
        ErrInfo.initialized(NGRDCOL),
    )


def test_solve_xp2_xpyp_with_multiple_lhs_compiles():
    _compile(
        solve_xp2_xpyp_with_multiple_lhs,
        NZM,
        NZT,
        NGRDCOL,
        1,
        jnp.array([1.0e-4], dtype=jnp.float64),
        _sample_grid(),
        SclrIdx(1, 0, 0, 0, 0, 0),
        _zm_array(1.0),
        _zm_array(1.0),
        _zm_array(1.0),
        _zm_array(1.0),
        _zm_array(0.1),
        _zt_array(0.01),
        _zt_array(300.0),
        _zm_array(0.01),
        _zm_array(0.01),
        _zm_array(0.0),
        _zm_array(0.0),
        _zm_array(0.0),
        _array((NGRDCOL, NZT, 1), 1.0),
        _array((NGRDCOL, NZM, 1), 0.01),
        _lhs3(0.0),
        _lhs3(0.0),
        _lhs3(0.0),
        _array((3, NGRDCOL, NZM, 1), 0.0),
        _array((3, NGRDCOL, NZM, 1), 0.0),
        _array((3, NGRDCOL, NZM, 1), 0.0),
        _lhs3(0.0),
        _lhs3(0.0),
        _zm_array(0.0),
        _zm_array(0.0),
        _zm_array(0.0),
        _array((NGRDCOL, NZM, 1), 0.0),
        _array((NGRDCOL, NZM, 1), 0.0),
        _array((NGRDCOL, NZM, 1), 0.0),
        1.0,
        model_flags.iiPDF_ADG1,
        True,
        False,
        True,
        model_flags.tridiag_lu,
        _stats(l_sample=False),
        _zm_array(1.0),
        _zm_array(1.0),
        _zm_array(0.01),
        _array((NGRDCOL, NZM, 1), 1.0),
        _array((NGRDCOL, NZM, 1), 0.01),
        _array((NGRDCOL, NZM, 1), 0.01),
        ErrInfo.initialized(NGRDCOL),
    )


def test_advance_xp2_xpyp_compiles():
    _compile(advance_xp2_xpyp, *_advance_xp2_xpyp_args())


def test_calc_ddzt_umvm_sqd_compiles():
    gr = _sample_grid()
    _compile(calc_ddzt_umvm_sqd, NZM, NZT, NGRDCOL, gr, _zt_array(1.0), _zt_array(2.0))


def test_calc_wp3_on_wp2_compiles():
    gr = _sample_grid()
    _compile(calc_wp3_on_wp2, NZM, NZT, NGRDCOL, gr, _zm_array(1.0), _zt_array(0.2))


def test_calc_stability_correction_compiles():
    _compile(
        calc_stability_correction,
        NZM,
        NGRDCOL,
        _zm_array(0.02),
        _zm_array(100.0),
        _zm_array(0.5),
        jnp.array([0.1, 0.2], dtype=jnp.float64),
    )


def test_calc_brunt_vaisala_freq_sqd_compiles():
    gr = _sample_grid()
    _compile(
        calc_brunt_vaisala_freq_sqd,
        NZM,
        NZT,
        NGRDCOL,
        gr,
        _zt_array(300.0),
        _zt_array(0.9),
        _zt_array(0.01),
        _zt_array(0.001),
        _zt_array(90000.0),
        _zt_array(301.0),
        _zt_array(0.5),
        model_flags.saturation_bolton,
        True,
        True,
        False,
        jnp.array([0.1, 0.2], dtype=jnp.float64),
        300.0,
    )


def test_calc_Ri_zm_compiles():
    _compile(calc_Ri_zm, NZM, NGRDCOL, _zm_array(0.02), _zm_array(0.5), 1.0e-7, 1.0e-7)


def test_compute_Cx_fnc_Richardson_compiles():
    gr = _sample_grid()
    _compile(
        compute_Cx_fnc_Richardson,
        NZM,
        NZT,
        NGRDCOL,
        gr,
        _zm_array(100.0),
        _zm_array(0.5),
        _zm_array(1.0),
        _zm_array(0.02),
        _zm_array(0.02),
        _clubb_params(),
        True,
        False,
    )


def test_Lscale_width_vert_avg_compiles():
    gr = _sample_grid()
    _compile(
        Lscale_width_vert_avg,
        NZM,
        NGRDCOL,
        gr,
        2,
        _zm_array(1.0),
        _zm_array(100.0),
        _zm_array(1.0),
        0.01,
    )


def test_wp23_term_splat_lhs_compiles():
    gr = _sample_grid()
    _compile(
        wp23_term_splat_lhs,
        NZM,
        NZT,
        NGRDCOL,
        gr,
        jnp.array([0.2, 0.3], dtype=jnp.float64),
        _zm_array(0.02),
        _zm_array(100.0),
        _zm_array(1.0),
    )


def test_smooth_min_compiles():
    _compile(smooth_min, NZM, NGRDCOL, _zm_array(1.0), _zm_array(2.0), 0.1)


def test_smooth_max_compiles():
    _compile(smooth_max, NZM, NGRDCOL, _zm_array(1.0), _zm_array(2.0), 0.1)


def test_smooth_heaviside_peskin_compiles():
    _compile(smooth_heaviside_peskin, NZM, NGRDCOL, _zm_array(0.0), 0.6)


def test_calc_xpwp_compiles():
    gr = _sample_grid()
    _compile(calc_xpwp, NZM, NZT, NGRDCOL, gr, _zm_array(0.1), _zt_array(1.0))


def test_calculate_thlp2_rad_compiles():
    gr = _sample_grid()
    _compile(
        calculate_thlp2_rad,
        NGRDCOL,
        NZM,
        NZT,
        gr,
        _zt_array(0.001),
        _zm_array(0.001),
        _zt_array(0.01),
        _clubb_params(),
        _zm_array(0.0),
    )


def test_clip_covar_compiles():
    _compile(
        clip_covar,
        NZM,
        NGRDCOL,
        clip_explicit.clip_wprtp,
        _zm_array(1.0),
        _zm_array(1.5),
        _zm_array(0.1),
    )


def test_clip_skewness_core_compiles():
    gr = _sample_grid()
    _compile(
        clip_skewness_core,
        NZT,
        NGRDCOL,
        gr,
        jnp.array([0.0, 10.0], dtype=jnp.float64),
        jnp.array([4.5, 4.5], dtype=jnp.float64),
        _zt_array(1.0),
        True,
        _zt_array(0.1),
    )


def test_clip_rcm_compiles():
    _compile(clip_rcm, NZT, NGRDCOL, _zt_array(0.01), "jit-compile-test", _zt_array(0.001))


def test_diffusion_zt_lhs_compiles():
    gr = _sample_grid()
    _compile(
        diffusion_zt_lhs,
        NZM,
        NZT,
        NGRDCOL,
        gr,
        _zm_array(0.1),
        _zt_array(0.1),
        jnp.array([0.01, 0.02], dtype=jnp.float64),
        1.0 / _zt_array(1.0),
        _zm_array(1.0),
    )


def test_diffusion_zm_lhs_compiles():
    gr = _sample_grid()
    _compile(
        diffusion_zm_lhs,
        NZM,
        NZT,
        NGRDCOL,
        gr,
        _zt_array(0.1),
        _zm_array(0.1),
        jnp.array([0.01, 0.02], dtype=jnp.float64),
        1.0 / _zm_array(1.0),
        _zt_array(1.0),
    )


def test_fill_holes_global_compiles():
    gr = _sample_grid()
    _compile(fill_holes_global, NZM, NGRDCOL, 0.2, 1, 4, gr.dzm, _zm_array(1.0), _zm_array(0.1))


def test_fill_holes_sliding_window_compiles():
    gr = _sample_grid()
    _compile(
        fill_holes_sliding_window,
        NZM,
        NGRDCOL,
        0.2,
        0,
        5,
        gr.dzm,
        _zm_array(1.0),
        _zm_array(0.1),
    )


def test_fill_holes_vertical_compiles():
    gr = _sample_grid()
    _compile(
        fill_holes_vertical,
        NZM,
        NGRDCOL,
        0.2,
        1,
        4,
        gr.dzm,
        _zm_array(1.0),
        1,
        model_flags.global_fill,
        _zm_array(0.1),
    )
    _compile(
        fill_holes_vertical,
        NZM,
        NGRDCOL,
        0.2,
        0,
        5,
        gr.dzm,
        _zm_array(1.0),
        1,
        model_flags.sliding_window,
        _zm_array(0.1),
    )


def test_fill_holes_wp2_from_horz_tke_compiles():
    _compile(
        fill_holes_wp2_from_horz_tke,
        NZM,
        NGRDCOL,
        0.2,
        1,
        4,
        _zm_array(0.1),
        _zm_array(0.5),
        _zm_array(0.6),
    )


def test_term_ma_zt_lhs_compiles():
    gr = _sample_grid()
    _compile(
        term_ma_zt_lhs,
        NZM,
        NZT,
        NGRDCOL,
        _zt_array(0.1),
        gr.weights_zt2zm,
        gr.invrs_dzt,
        gr.invrs_dzm,
        False,
        1.0,
    )
    _compile(
        term_ma_zt_lhs,
        NZM,
        NZT,
        NGRDCOL,
        _zt_array(0.1),
        gr.weights_zt2zm,
        gr.invrs_dzt,
        gr.invrs_dzm,
        True,
        1.0,
    )


def test_term_ma_zm_lhs_compiles():
    gr = _sample_grid()
    _compile(term_ma_zm_lhs, NZM, NZT, NGRDCOL, _zm_array(0.1), gr.invrs_dzm, gr.weights_zm2zt)


def test_penta_lu_solve_single_rhs_multiple_lhs_compiles():
    _compile(penta_lu_solve_single_rhs_multiple_lhs, NZM, NGRDCOL, _penta_lhs_multi(), _zm_array(1.0))


def test_penta_lu_solve_multiple_rhs_lhs_compiles():
    _compile(
        penta_lu_solve_multiple_rhs_lhs,
        NZM,
        NRHS,
        NGRDCOL,
        _penta_lhs_multi(),
        _array((NGRDCOL, NZM, NRHS), 1.0),
    )


def test_tridiag_lu_solve_single_rhs_lhs_compiles():
    _compile(tridiag_lu_solve_single_rhs_lhs, NZM, _tridiag_lhs_single(), _array((NZM,), 1.0))


def test_tridiag_lu_solve_single_rhs_multiple_lhs_compiles():
    _compile(
        tridiag_lu_solve_single_rhs_multiple_lhs,
        NZM,
        NGRDCOL,
        _tridiag_lhs_multi(),
        _zm_array(1.0),
    )


def test_tridiag_lu_solve_multiple_rhs_lhs_compiles():
    _compile(
        tridiag_lu_solve_multiple_rhs_lhs,
        NZM,
        NRHS,
        NGRDCOL,
        _tridiag_lhs_multi(),
        _array((NGRDCOL, NZM, NRHS), 1.0),
    )


def test_xpyp_term_ta_pdf_lhs_compiles():
    gr = _sample_grid()
    _compile(
        xpyp_term_ta_pdf_lhs,
        NZM,
        NZT,
        NGRDCOL,
        gr,
        _zt_array(1.0),
        _zt_array(1.0),
        _zm_array(1.0),
        1.0 / _zm_array(1.0),
        False,
        _zm_array(1.0),
        _zm_array(1.0),
    )
    _compile(
        xpyp_term_ta_pdf_lhs,
        NZM,
        NZT,
        NGRDCOL,
        gr,
        _zt_array(1.0),
        _zt_array(1.0),
        _zm_array(1.0),
        1.0 / _zm_array(1.0),
        True,
        _zm_array(1.0),
        _zm_array(1.0),
    )


def test_xpyp_term_ta_pdf_lhs_godunov_compiles():
    gr = _sample_grid()
    _compile(xpyp_term_ta_pdf_lhs_godunov, NZM, NZT, NGRDCOL, gr, _zt_array(1.0), 1.0 / _zm_array(1.0), _zm_array(1.0))


def test_xpyp_term_ta_pdf_rhs_compiles():
    gr = _sample_grid()
    _compile(
        xpyp_term_ta_pdf_rhs,
        NZM,
        NZT,
        NGRDCOL,
        gr,
        _zt_array(1.0),
        _zt_array(1.0),
        _zm_array(1.0),
        1.0 / _zm_array(1.0),
        False,
        _zm_array(1.0),
        _zm_array(1.0),
    )
    _compile(
        xpyp_term_ta_pdf_rhs,
        NZM,
        NZT,
        NGRDCOL,
        gr,
        _zt_array(1.0),
        _zt_array(1.0),
        _zm_array(1.0),
        1.0 / _zm_array(1.0),
        True,
        _zm_array(1.0),
        _zm_array(1.0),
    )


def test_xpyp_term_ta_pdf_rhs_godunov_compiles():
    gr = _sample_grid()
    _compile(
        xpyp_term_ta_pdf_rhs_godunov,
        NZM,
        NZT,
        NGRDCOL,
        gr,
        _zm_array(1.0),
        1.0 / _zm_array(1.0),
        _zt_array(1.0),
        _zm_array(1.0),
    )


def main() -> int:
    test_derived_type_pytree_metadata_is_static()
    test_xp2_xpyp_single_lhs_valid_compiles()
    test_term_tp_rhs_compiles()
    test_term_dp1_lhs_compiles()
    test_term_dp1_rhs_compiles()
    test_term_pr1_compiles()
    test_term_pr2_compiles()
    test_xp2_xpyp_lhs_compiles()
    test_xp2_xpyp_solve_compiles()
    test_stats_finalize_xp2_xpyp_terms_compiles()
    test_xp2_xpyp_uv_rhs_compiles()
    test_xp2_xpyp_rhs_compiles()
    test_calc_xp2_xpyp_ta_terms_compiles()
    test_pos_definite_variances_compiles()
    test_solve_xp2_xpyp_with_single_lhs_compiles()
    test_solve_xp2_xpyp_with_multiple_lhs_compiles()
    test_advance_xp2_xpyp_compiles()
    test_calc_ddzt_umvm_sqd_compiles()
    test_calc_wp3_on_wp2_compiles()
    test_calc_stability_correction_compiles()
    test_calc_brunt_vaisala_freq_sqd_compiles()
    test_calc_Ri_zm_compiles()
    test_compute_Cx_fnc_Richardson_compiles()
    test_Lscale_width_vert_avg_compiles()
    test_wp23_term_splat_lhs_compiles()
    test_smooth_min_compiles()
    test_smooth_max_compiles()
    test_smooth_heaviside_peskin_compiles()
    test_calc_xpwp_compiles()
    test_calculate_thlp2_rad_compiles()
    test_clip_covar_compiles()
    test_clip_skewness_core_compiles()
    test_clip_rcm_compiles()
    test_diffusion_zt_lhs_compiles()
    test_diffusion_zm_lhs_compiles()
    test_fill_holes_global_compiles()
    test_fill_holes_sliding_window_compiles()
    test_fill_holes_vertical_compiles()
    test_fill_holes_wp2_from_horz_tke_compiles()
    test_term_ma_zt_lhs_compiles()
    test_term_ma_zm_lhs_compiles()
    test_penta_lu_solve_single_rhs_multiple_lhs_compiles()
    test_penta_lu_solve_multiple_rhs_lhs_compiles()
    test_tridiag_lu_solve_single_rhs_lhs_compiles()
    test_tridiag_lu_solve_single_rhs_multiple_lhs_compiles()
    test_tridiag_lu_solve_multiple_rhs_lhs_compiles()
    test_xpyp_term_ta_pdf_lhs_compiles()
    test_xpyp_term_ta_pdf_lhs_godunov_compiles()
    test_xpyp_term_ta_pdf_rhs_compiles()
    test_xpyp_term_ta_pdf_rhs_godunov_compiles()
    print("Compiled explicit representative JIT call signatures.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
