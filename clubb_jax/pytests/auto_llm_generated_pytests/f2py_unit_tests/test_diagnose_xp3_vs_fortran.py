"""Check advance_xp3_module.py:diagnose_xp3 against Fortran diagnose_xp3."""

from input.scripts.create_case_namelist import create_case_namelist_file

import os

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.advance_xp3 as fortran

from utilities.output_paths import REPO_ROOT

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff, no_stats, setup_matching_grids
from clubb_jax.src.CLUBB_core.advance_xp3_module import diagnose_xp3
from clubb_jax.src.CLUBB_core.model_flags import iiPDF_ADG1, iiPDF_new
from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params

_TUNABLE_PARAMS = str(create_case_namelist_file("bomex"))
_NG, _DZ, _ZTOP = 2, 40.0, 1200.0
_SCLR_TOL = np.array([1e-8, 1e-2])
_OUTPUTS = ("rtp3", "thlp3", "up3", "vp3", "sclrp3")


def _inputs(rng, ng, nzm, nzt, sclr_dim):
    zt_frac = np.linspace(0.0, 1.0, nzt)
    return dict(
        wp2=rng.uniform(1e-3, 2.0, (ng, nzm)),
        wp3=rng.uniform(-1.0, 1.0, (ng, nzt)),
        # Stably stratified virtual potential temperature, so the non-ADG1
        # branch sees a realistic Brunt-Vaisala frequency.
        thvm=300.0 + 6.0 * zt_frac + rng.uniform(-0.2, 0.2, (ng, nzt)),
        wprtp=rng.uniform(-1e-3, 1e-3, (ng, nzm)),
        wpthlp=rng.uniform(-0.5, 0.5, (ng, nzm)),
        rtp2=rng.uniform(1e-8, 1e-5, (ng, nzm)),
        thlp2=rng.uniform(1e-3, 1.0, (ng, nzm)),
        upwp=rng.uniform(-0.5, 0.5, (ng, nzm)),
        vpwp=rng.uniform(-0.5, 0.5, (ng, nzm)),
        up2=rng.uniform(1e-3, 2.0, (ng, nzm)),
        vp2=rng.uniform(1e-3, 2.0, (ng, nzm)),
        sigma_sqd_w=rng.uniform(0.0, 0.5, (ng, nzm)),
        wpsclrp=rng.uniform(-0.5, 0.5, (ng, nzm, sclr_dim)),
        sclrp2=rng.uniform(1e-3, 1.0, (ng, nzm, sclr_dim)),
    )


def test_f2py_oracle():
    """diagnose_xp3 matches Fortran for the ADG1 ansatz and the xp3_coef_fnc branch, including scalars."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    sclr_dim = _SCLR_TOL.size
    clubb_params = np.asarray(init_clubb_params(ng, filename=_TUNABLE_PARAMS))
    rng = np.random.default_rng(5)
    for iipdf_type in (iiPDF_ADG1, iiPDF_new):
        worst = {}
        for _ in range(6):
            a = _inputs(rng, ng, nzm, nzt, sclr_dim)
            moments = [a[k] for k in ("wp2", "wp3", "thvm", "wprtp", "wpthlp", "rtp2", "thlp2",
                                      "upwp", "vpwp", "up2", "vp2", "sigma_sqd_w", "wpsclrp", "sclrp2")]
            xp3_in = [np.zeros((ng, nzt))] * 4 + [np.zeros((ng, nzt, sclr_dim))]
            *j_out, _stats = diagnose_xp3(
                nzm, nzt, ng, sclr_dim, jnp.asarray(_SCLR_TOL), jgr, iipdf_type, jnp.asarray(clubb_params),
                *(jnp.asarray(x) for x in moments), no_stats(ng), *(jnp.asarray(x) for x in xp3_in))
            f_out = fortran.diagnose_xp3(
                nzm, nzt, ng, sclr_dim, _SCLR_TOL, fgr, iipdf_type, clubb_params,
                *moments, *(x.copy() for x in xp3_in))
            assert len(j_out) == len(f_out) == len(_OUTPUTS), (
                f"output counts: JAX {len(j_out)}, Fortran {len(f_out)}, expected {len(_OUTPUTS)}")
            for name, j, f in zip(_OUTPUTS, j_out, f_out, strict=True):
                worst[name] = max(worst.get(name, 0.0), max_abs_diff(name, j, f))
        bad = {k: v for k, v in worst.items() if not v < 1e-12}
        assert not bad, f"diagnose_xp3 (iiPDF_type={iipdf_type}) Fortran mismatch: {bad}"
