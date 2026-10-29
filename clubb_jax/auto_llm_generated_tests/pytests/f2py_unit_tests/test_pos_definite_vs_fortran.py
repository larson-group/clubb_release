"""Check the JAX pos_definite_adj port (Smolarkiewicz 1989 flux limiter)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.pos_definite_module as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_rel_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.pos_definite_module import pos_definite_adj

_NG, _DZ, _ZTOP, _DT = 2, 40.0, 1200.0, 60.0


def _make_inputs(jgr):
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rng = np.random.default_rng(20260603)
    # Flux on zm, field on zt. Make column 0 develop a below-zero field (exercise the limiter branches),
    # leave column 1 all-positive (exercise the "no negative -> *_pd = 0" branch).
    flux_np1 = rng.standard_normal((ng, nzm)) * 1e-3
    field_n = np.abs(rng.standard_normal((ng, nzt))) * 1e-2 + 1e-3
    field_np1 = field_n + rng.standard_normal((ng, nzt)) * 5e-3
    field_np1[0, nzt // 3] = -2e-3      # a genuine negative in column 0
    field_np1[1, :] = np.abs(field_np1[1, :]) + 1e-4   # column 1 strictly positive
    return (np.asfortranarray(field_np1), np.asfortranarray(flux_np1), np.asfortranarray(field_n))


def test_f2py_oracle():
    """pos_definite_adj (flux-conservative positive-definite limiter) matches Fortran in all four outputs."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    field_np1, flux_np1, field_n = _make_inputs(jgr)
    f_out = fortran.pos_definite_adj(nzm, nzt, ng, fgr, _DT, field_np1.copy(), flux_np1.copy(), field_n.copy())
    j_out = pos_definite_adj(nzm, nzt, ng, jgr, _DT, field_np1, flux_np1, field_n)
    for name, jx, fx in zip(("field_np1", "flux_np1", "field_pd", "flux_pd"), j_out, f_out):
        jx, fx = np.asarray(jx), np.asarray(fx)
        assert jx.shape == fx.shape, f"{name}: JAX shape {jx.shape} != Fortran {fx.shape}"
        rel = max_rel_diff("pos_definite", jx, fx, 1e-30)
        assert rel < 1e-11, f"{name} Fortran mismatch rel {rel:.2e}"
