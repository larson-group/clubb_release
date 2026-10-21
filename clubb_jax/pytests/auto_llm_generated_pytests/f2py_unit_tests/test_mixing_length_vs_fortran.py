"""Check mixing_length.F90 ports against Fortran."""

import os

import numpy as np

import clubb_python.CLUBB_core.mixing_length as fortran
from clubb_python.CLUBB_core.pdf_parameter_module import init_pdf_params as fortran_init_pdf_params
from clubb_python.derived_types.err_info import ErrInfo as FortranErrInfo

from utilities.output_paths import REPO_ROOT as _REPO_ROOT

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff, max_rel_diff, no_stats, setup_matching_grids
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.mixing_length import (
    set_Lscale_max,
    calc_Lscale_directly,
    diagnose_Lscale_from_tau,
)
from clubb_jax.src.CLUBB_core.model_flags import Lscale_direct
from clubb_jax.src.CLUBB_core.parameter_indices import iLscale_mu_entrain
from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params
from clubb_jax.src.CLUBB_core.pdf_params import init_pdf_params
from clubb_jax.src.CLUBB_core.saturation import SATURATION_FLATAU

_ROOT = str(_REPO_ROOT)
_TUNABLE_PARAMS = os.path.join(_ROOT, "input", "parameter_and_flag_configs", "default", "tunable_parameters.in")

_NG, _DZ, _ZTOP = 2, 40.0, 1600.0


def test_set_Lscale_max_f2py():
    """set_Lscale_max matches Fortran for host-model and standalone configurations."""
    rng = np.random.default_rng(1)
    worst = 0.0
    for l_impl in (False, True):
        for _ in range(10):
            dx = rng.uniform(100.0, 5e5, 3); dy = rng.uniform(100.0, 5e5, 3)
            worst = max(worst, max_abs_diff("set_Lscale_max", set_Lscale_max(3, l_impl, dx, dy),
                                            fortran.set_lscale_max(3, l_impl, dx, dy)))
    assert worst < 1e-12, f"set_Lscale_max Fortran mismatch {worst:.2e}"


def test_calc_Lscale_directly_f2py():
    """calc_Lscale_directly (parcel-ascent Lscale, Lscale_up, Lscale_down) matches Fortran.

    The parcel iteration accumulates ~1e-13 relative rounding differences, so a relative bound is used."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    worst_rel = 0.0
    for seed in range(6):
        rng = np.random.default_rng(seed)
        thlm = np.sort(rng.uniform(290.0, 330.0, (ng, nzt)), axis=1); exner = rng.uniform(0.8, 1.0, (ng, nzt))
        rtm = rng.uniform(2e-3, 1.4e-2, (ng, nzt)); thvm = thlm * (1.0 + 0.6 * rtm)
        em = rng.uniform(1e-2, 1.0, (ng, nzm)); p = np.sort(rng.uniform(6e4, 1e5, (ng, nzt)), axis=1)[:, ::-1].copy()
        thv_ds = thvm.copy(); Lmax = np.full(ng, 1e5); lmin = 20.0
        cp = np.array(init_clubb_params(ng, filename=_TUNABLE_PARAMS))
        cp[:, iLscale_mu_entrain] = rng.uniform(0.4, 1.0, ng)
        z = np.zeros((ng, nzt))
        common = (False, p, exner, rtm, thlm, thvm, z, z, z)
        tail = (em, thv_ds, Lmax, lmin, cp, SATURATION_FLATAU, False, Lscale_direct)
        fr = fortran.calc_lscale_directly(ng, nzm, nzt, fgr, *common, fortran_init_pdf_params(nzt, ng), *tail,
                                          FortranErrInfo(ngrdcol=ng))
        jr = calc_Lscale_directly(ng, nzm, nzt, jgr, *common, init_pdf_params(nzt, ng), *tail,
                                  no_stats(ng), ErrInfo(ngrdcol=ng))
        assert not fr[0].is_fatal() and not jr[0].is_fatal()
        for name, a, b in zip(("Lscale", "Lscale_up", "Lscale_down"), jr[1:4], fr[1:4]):
            worst_rel = max(worst_rel, max_rel_diff(name, a, b, 1.0))
    assert worst_rel < 1e-9, f"calc_Lscale_directly Fortran rel mismatch {worst_rel:.2e}"


def test_diagnose_Lscale_from_tau_f2py():
    """diagnose_Lscale_from_tau (the l_diag_Lscale_from_tau path) matches Fortran in all tau/Lscale outputs.

    No validated case enables l_diag_Lscale_from_tau, so this is the routine's only Fortran comparison."""
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    cp = np.asarray(init_clubb_params(ng, filename=_TUNABLE_PARAMS))   # realistic defaults (well-conditioned tau)
    worst = 0.0
    for seed in range(4):
        rng = np.random.default_rng(seed)
        upwp = rng.uniform(-0.1, 0.1, ng); vpwp = rng.uniform(-0.1, 0.1, ng)
        ddzt = rng.uniform(0.0, 1e-3, (ng, nzm)); isf = rng.uniform(0.0, 0.3, (ng, nzt))
        em = rng.uniform(1e-2, 1.0, (ng, nzm))
        sfce = np.zeros(ng); Lmax = np.full(ng, 1e5)
        Ri = rng.uniform(-1.0, 5.0, (ng, nzm)); bvs = rng.uniform(-1e-3, 1e-2, (ng, nzm))
        head = (upwp, vpwp, ddzt, isf, em, 0.01, 900.0, sfce, Lmax, cp)
        for l_smth in (False, True):
            fr = fortran.diagnose_lscale_from_tau(nzm, nzt, ng, fgr, *head, False, l_smth, bvs, Ri,
                                                  FortranErrInfo(ngrdcol=ng))
            jr = diagnose_Lscale_from_tau(nzm, nzt, ng, jgr, *head, no_stats(ng), False, l_smth, bvs, Ri,
                                          ErrInfo(ngrdcol=ng))
            n_out = len(fr) - 1
            assert n_out == 19, f"unexpected Fortran output count {n_out}"
            for k in range(1, n_out + 1):
                worst = max(worst, max_rel_diff(f"diagnose_Lscale_from_tau output {k}", jr[k], fr[k], 1e-3))
    assert worst < 1e-9, f"diagnose_Lscale_from_tau Fortran rel mismatch {worst:.2e}"
