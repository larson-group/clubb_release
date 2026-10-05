"""Check static pytree metadata and compiled single-LHS validity."""

from __future__ import annotations
from utilities.output_paths import REPO_ROOT as _REPO_ROOT

import math

REPO_ROOT = _REPO_ROOT


import jax

jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core import model_flags, parameter_indices
from clubb_jax.src.CLUBB_core.advance_xp2_xpyp_module import xp2_xpyp_single_lhs_valid
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


def _nu_vert_res_dep():
    return NuVertResDep(NZM, *([_nu_array(0.01)] * 7))


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
