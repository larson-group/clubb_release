"""Check advance_xp3_module.py:advance_xp3 against Fortran advance_xp3.

advance_xp3_simplified is private in the Fortran module, so it is checked
through advance_xp3, which calls it for rtp3, thlp3 and each scalar.
"""

import os

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.advance_xp3 as fortran

from utilities.output_paths import REPO_ROOT

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff, no_stats, setup_matching_grids
from clubb_jax.src.CLUBB_core.advance_xp3_module import advance_xp3
from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params

_TUNABLE_PARAMS = os.path.join(REPO_ROOT, "input", "parameter_and_flag_configs", "default",
                               "tunable_parameters.in")
_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
_DT = 60.0
_SCLR_TOL = np.array([1e-8, 1e-2])
_OUTPUTS = ("rtp3", "thlp3", "sclrp3", "up3", "vp3")
# Argument order shared by the JAX and Fortran advance_xp3 after dt.
_FIELDS = ("rtm", "thlm", "rtp2", "thlp2", "wprtp", "wpthlp", "wprtp2", "wpthlp2",
           "rho_ds_zm", "invrs_rho_ds_zt", "invrs_tau_zt", "tau_max_zt",
           "sclrm", "sclrp2", "wpsclrp", "wpsclrp2",
           "wp2", "wp3", "upwp", "vpwp", "up2", "vp2", "thvm")


def _inputs(rng, ng, nzm, nzt, sclr_dim):
    zt_frac = np.linspace(0.0, 1.0, nzt)
    return dict(
        rtm=rng.uniform(1e-3, 2e-2, (ng, nzt)),
        thlm=290.0 + 10.0 * zt_frac + rng.uniform(-0.5, 0.5, (ng, nzt)),
        rtp2=rng.uniform(1e-8, 1e-5, (ng, nzm)),
        thlp2=rng.uniform(1e-3, 1.0, (ng, nzm)),
        wprtp=rng.uniform(-1e-3, 1e-3, (ng, nzm)),
        wpthlp=rng.uniform(-0.5, 0.5, (ng, nzm)),
        wprtp2=rng.uniform(-1e-6, 1e-6, (ng, nzt)),
        wpthlp2=rng.uniform(-0.1, 0.1, (ng, nzt)),
        rho_ds_zm=rng.uniform(0.5, 1.2, (ng, nzm)),
        invrs_rho_ds_zt=rng.uniform(0.8, 2.0, (ng, nzt)),
        invrs_tau_zt=rng.uniform(1e-3, 1e-2, (ng, nzt)),
        tau_max_zt=rng.uniform(300.0, 900.0, (ng, nzt)),
        sclrm=rng.uniform(0.0, 5.0, (ng, nzt, sclr_dim)),
        sclrp2=rng.uniform(1e-3, 1.0, (ng, nzm, sclr_dim)),
        wpsclrp=rng.uniform(-0.5, 0.5, (ng, nzm, sclr_dim)),
        wpsclrp2=rng.uniform(-0.1, 0.1, (ng, nzt, sclr_dim)),
        wp2=rng.uniform(1e-3, 2.0, (ng, nzm)),
        wp3=rng.uniform(-1.0, 1.0, (ng, nzt)),
        upwp=rng.uniform(-0.5, 0.5, (ng, nzm)),
        vpwp=rng.uniform(-0.5, 0.5, (ng, nzm)),
        up2=rng.uniform(1e-3, 2.0, (ng, nzm)),
        vp2=rng.uniform(1e-3, 2.0, (ng, nzm)),
        # Stably stratified, so the up3/vp3 xp3_coef_fnc sees a realistic Brunt-Vaisala frequency.
        thvm=300.0 + 6.0 * zt_frac + rng.uniform(-0.2, 0.2, (ng, nzt)),
    )


def test_f2py_oracle():
    """advance_xp3 (steady-state rtp3/thlp3/sclrp3 and diagnosed up3/vp3) matches Fortran."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    sclr_dim = _SCLR_TOL.size
    clubb_params = np.asarray(init_clubb_params(ng, filename=_TUNABLE_PARAMS))
    rng = np.random.default_rng(12)
    for l_lmm_stepping in (False, True):
        worst = {}
        for _ in range(6):
            a = _inputs(rng, ng, nzm, nzt, sclr_dim)
            fields = [a[k] for k in _FIELDS]
            xp3_in = [rng.uniform(-1.0, 1.0, (ng, nzt)) for _ in range(2)] + [
                rng.uniform(-1.0, 1.0, (ng, nzt, sclr_dim))] + [
                rng.uniform(-1.0, 1.0, (ng, nzt)) for _ in range(2)]
            *j_out, _stats = advance_xp3(
                nzm, nzt, ng, sclr_dim, jnp.asarray(_SCLR_TOL), jgr, _DT,
                *(jnp.asarray(x) for x in fields), jnp.asarray(clubb_params), l_lmm_stepping,
                no_stats(ng), *(jnp.asarray(x) for x in xp3_in))
            f_out = fortran.advance_xp3(
                nzm, nzt, ng, sclr_dim, _SCLR_TOL, fgr, _DT, *fields, clubb_params, l_lmm_stepping,
                *(x.copy() for x in xp3_in))
            assert len(j_out) == len(f_out) == len(_OUTPUTS), (
                f"output counts: JAX {len(j_out)}, Fortran {len(f_out)}, expected {len(_OUTPUTS)}")
            for name, j, f in zip(_OUTPUTS, j_out, f_out, strict=True):
                worst[name] = max(worst.get(name, 0.0), max_abs_diff(name, j, f))
        bad = {k: v for k, v in worst.items() if not v < 1e-12}
        assert not bad, f"advance_xp3 (l_lmm_stepping={l_lmm_stepping}) Fortran mismatch: {bad}"
