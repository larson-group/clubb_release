"""Check the JAX compute_sigma_sqd_w port (sigma_sqd_w_module.F90)
against Fortran."""

from input.scripts.create_case_namelist import create_case_namelist_file

import os

import numpy as np

import clubb_python.CLUBB_core.sigma_sqd_w_module as fortran

from utilities.output_paths import REPO_ROOT

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.sigma_sqd_w_module import compute_sigma_sqd_w

_TUNABLE_PARAMS = str(create_case_namelist_file("bomex"))


_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


def _fields(nzm, nzt, seed):
    rng = np.random.default_rng(seed)
    sh = (_NG, nzm)
    return dict(wp3=rng.uniform(-2.0, 2.0, (_NG, nzt)),
                wp2=rng.uniform(0.05, 2.0, sh), thlp2=rng.uniform(0.01, 1.0, sh),
                rtp2=rng.uniform(1e-8, 1e-6, sh), up2=rng.uniform(0.05, 2.0, sh),
                vp2=rng.uniform(0.05, 2.0, sh), wpthlp=rng.uniform(-0.5, 0.5, sh),
                wprtp=rng.uniform(-1e-4, 1e-4, sh), upwp=rng.uniform(-0.5, 0.5, sh),
                vpwp=rng.uniform(-0.5, 0.5, sh))


def test_f2py_oracle():
    """compute_sigma_sqd_w (including its internal Skw and gamma(Skw)) matches Fortran, with
    l_predict_upwp_vpwp on and off."""
    from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    params = np.asarray(init_clubb_params(ng, filename=_TUNABLE_PARAMS))
    order = ('wp3', 'wp2', 'thlp2', 'rtp2', 'up2', 'vp2', 'wpthlp', 'wprtp', 'upwp', 'vpwp')
    worst = 0.0
    for seed in (5, 6):
        f = _fields(nzm, nzt, seed)
        for l_pred in (False, True):
            ref = fortran.compute_sigma_sqd_w(nzm, nzt, ng, fgr, *[f[k] for k in order], params, l_pred)
            got = compute_sigma_sqd_w(nzm, nzt, ng, jgr, *[f[k] for k in order], params, l_pred)
            worst = max(worst, max_abs_diff(f"compute_sigma_sqd_w(l_predict_upwp_vpwp={l_pred})", got, ref))
    assert worst < 1e-11, f"compute_sigma_sqd_w Fortran mismatch {worst:.2e}"
