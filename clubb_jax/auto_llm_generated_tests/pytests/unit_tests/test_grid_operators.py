#!/usr/bin/env python3
"""validate the JAX vertical grid operators."""


import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.grid_class import setup_grid

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def test_descending_grid_rejected():
    """`setup_grid` fail-loud rejects `l_ascending_grid=False`: the JAX grid operators (zt2zm boundary handling) are
    ascending-only, so a descending grid would silently mis-compute the two boundary levels (verified iter 483 — interior
    interp is correct, both boundaries diverge ~0.5·Δfield from the Fortran). No case uses a descending grid; the guard
    makes the unsupported case explicit instead of a silent footgun. (iter 483)"""
    try:
        setup_grid(ngrdcol=2, deltaz=_DZ, zm_init=0.0, zm_top=_ZTOP, grid_type=1, l_ascending_grid=False)
    except ValueError as e:
        assert "ascending" in str(e).lower(), f"unexpected rejection message: {e}"
        print("  setup_grid fail-loud rejects descending grid (l_ascending_grid=False)  PASS")
        return
    raise AssertionError("setup_grid did NOT reject l_ascending_grid=False (descending grid)")
