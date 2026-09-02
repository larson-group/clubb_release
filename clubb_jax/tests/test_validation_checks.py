#!/usr/bin/env python3
"""test_validation_checks.py — validate the ported CLUBB check routines.

Covers the in-scope NaN/validation check routines (corr_varnce_module.F90 / numerical_check.F90). Their f2py
wrappers set err_code on `stored_err_info` and `return` (no error-stop, and the f2py does NOT expose err_code), so
there is no observable f2py oracle — validation is by the Fortran source logic (behavioral) plus an f2py no-crash
cross-check on the valid case.
  assert_corr_symmetric: True iff symmetric within 1e-6 AND unit diagonal within eps (1e-10).
  sfc_varnce_check:       True iff every surface variance/covariance is finite.
  length_check:           True iff Lscale/Lscale_up/Lscale_down are all finite.
  pdf_closure_check:      True iff every pdf_closure output AND pdf_params component is finite.
  rad_check:              True iff every radiation input (incl. rvm=rtm-rcm) is non-negative.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for p in (_ROOT, _ROOT + "/clubb_python_api"):
    if p not in sys.path:
        sys.path.append(p)

import numpy as np

from clubb_jax.src.CLUBB_core.corr_varnce_module import assert_corr_symmetric
from clubb_jax.src.CLUBB_core.numerical_check import (
    sfc_varnce_check, length_check, pdf_closure_check, rad_check, invalid_model_arrays,
)
from clubb_jax.src.CLUBB_core.pdf_params import init_pdf_params


def test_assert_corr_symmetric_behavior():
    n = 5
    rng = np.random.default_rng(2)
    lower = rng.uniform(-0.5, 0.5, (n, n))
    sym = np.tril(lower, -1) + np.tril(lower, -1).T + np.eye(n)   # symmetric, unit diagonal
    assert assert_corr_symmetric(sym) is True, "valid symmetric unit-diagonal matrix rejected"
    # Asymmetric (perturb one off-diagonal beyond tol)
    asym = sym.copy(); asym[0, 1] += 1e-3
    assert assert_corr_symmetric(asym) is False, "asymmetric matrix accepted"
    # Within tolerance (perturb below tol) → still accepted
    near = sym.copy(); near[0, 1] += 5e-7; near[1, 0] += 0.0
    assert assert_corr_symmetric(near) is True, "within-tol asymmetry wrongly rejected"
    # Non-unit diagonal
    nd = sym.copy(); nd[2, 2] = 0.9
    assert assert_corr_symmetric(nd) is False, "non-unit-diagonal matrix accepted"
    print("  assert_corr_symmetric: symmetric+unit-diag PASS / asym + non-unit-diag FAIL  PASS")


def test_assert_corr_symmetric_f2py_nocrash():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py assert_corr_symmetric (no-crash valid case): SKIP ({type(e).__name__})")
        return
    n = 6
    rng = np.random.default_rng(7)
    lower = rng.uniform(-0.4, 0.4, (n, n))
    sym = np.tril(lower, -1) + np.tril(lower, -1).T + np.eye(n)
    # The f2py sets err_code on stored_err_info and returns (no error-stop). A valid matrix must not raise.
    clubb_f2py.f2py_assert_corr_symmetric(np.asfortranarray(sym))
    assert assert_corr_symmetric(sym) is True
    print("  f2py assert_corr_symmetric: valid matrix runs without error + JAX agrees  PASS")


def test_sfc_varnce_check_behavior():
    ng = 1
    ok = dict(wp2_sfc=0.3, up2_sfc=0.2, vp2_sfc=0.2, thlp2_sfc=0.5, rtp2_sfc=1e-6, rtpthlp_sfc=1e-4)
    assert sfc_varnce_check(0, **ok) is True, "all-finite surface variances rejected"
    bad = dict(ok); bad['thlp2_sfc'] = np.nan
    assert sfc_varnce_check(0, **bad) is False, "NaN thlp2_sfc accepted"
    inf = dict(ok); inf['wp2_sfc'] = np.inf
    assert sfc_varnce_check(0, **inf) is False, "Inf wp2_sfc accepted"
    # scalar fields finite but a NaN passive scalar with sclr_dim>0
    assert sfc_varnce_check(1, sclrp2_sfc=np.array([np.nan]), sclrprtp_sfc=np.array([0.0]),
                            sclrpthlp_sfc=np.array([0.0]), **ok) is False, "NaN scalar variance accepted"
    print("  sfc_varnce_check: all-finite PASS / NaN+Inf FAIL (incl. passive scalars)  PASS")


def test_length_check_behavior():
    g = np.ones(8)
    assert length_check(g, g, g) is True, "all-finite mixing lengths rejected"
    assert length_check(g, np.array([np.nan] + [1.0] * 7), g) is False, "NaN Lscale_up accepted"
    assert length_check(g, g, np.array([np.inf] + [1.0] * 7)) is False, "Inf Lscale_down accepted"
    print("  length_check: all-finite PASS / NaN+Inf FAIL  PASS")


def test_pdf_closure_check_behavior():
    pp = init_pdf_params(10, 1)   # nz=10, ngrdcol=1 → all-zero (finite) pdf_params components
    closure = dict(wp4=np.zeros(10), wp2rtp=np.zeros(10), rcm=np.zeros(10), crt_1=np.ones(10))
    assert pdf_closure_check(closure, pp, sclr_dim=0) is True, "all-finite pdf_closure outputs rejected"
    bad = dict(closure); bad['rcm'] = np.array([np.nan] + [0.0] * 9)
    assert pdf_closure_check(bad, pp, sclr_dim=0) is False, "NaN closure output accepted"
    # NaN inside a pdf_params component must also fail
    pp_bad = pp._replace(w_1=np.array([[np.nan] + [0.0] * 9]))
    assert pdf_closure_check(closure, pp_bad, sclr_dim=0) is False, "NaN pdf_params%w_1 accepted"
    # NaN passive-scalar array with sclr_dim>0
    assert pdf_closure_check(closure, pp, sclr_dim=1,
                             sclr_fields=dict(sclrpthvp=np.array([np.nan]))) is False, "NaN scalar accepted"
    print("  pdf_closure_check: all-finite PASS / NaN closure+pdf_params+scalar FAIL  PASS")


def test_rad_check_behavior():
    n = 6
    ok = dict(thlm=300 * np.ones(n), rcm=np.zeros(n), rtm=1e-3 * np.ones(n), rim=np.zeros(n),
              cloud_frac=0.5 * np.ones(n), p_in_Pa=1e5 * np.ones(n), exner=np.ones(n), rho_zm=np.ones(n))
    assert rad_check(**ok) is True, "all-nonnegative radiation inputs rejected"
    neg = dict(ok); neg['rim'] = np.array([-1e-6] + [0.0] * (n - 1))
    assert rad_check(**neg) is False, "negative rim accepted"
    # rvm = rtm - rcm < 0 must also fail (rcm > rtm)
    rvm_neg = dict(ok); rvm_neg['rcm'] = 2e-3 * np.ones(n)
    assert rad_check(**rvm_neg) is False, "negative rvm (rcm>rtm) accepted"
    print("  rad_check: all-nonnegative PASS / negative input+rvm FAIL  PASS")


def test_invalid_model_arrays_behavior():
    n = 8
    ok = {k: np.ones(n) for k in ('um', 'vm', 'rtm', 'wprtp', 'thlm', 'wpthlp', 'rtp2', 'thlp2',
                                  'rtpthlp', 'wp2', 'wp3', 'wp2thvp', 'wp2up', 'rtpthvp', 'thlpthvp')}
    # NOTE inverse polarity: returns True if INVALID (matching the Fortran name + driver `if invalid(...)`).
    assert invalid_model_arrays(**ok) is False, "all-finite arrays flagged invalid"
    bad = dict(ok); bad['wp3'] = np.array([np.inf] + [1.0] * (n - 1))
    assert invalid_model_arrays(**bad) is True, "Inf in wp3 not flagged"
    hm = np.ones((n, 2)); hm[2, 1] = np.nan
    assert invalid_model_arrays(**ok, hydromet=hm, hydromet_list=['rrm', 'Nrm']) is True, "NaN hydromet not flagged"
    assert invalid_model_arrays(**ok, sclrm=np.ones((n, 1)), edsclrm=np.ones((n, 1))) is False
    print("  invalid_model_arrays: clean->False / Inf wp3 + NaN hydromet->True (inverse polarity)  PASS")


def main():
    print("test_validation_checks:")
    for t in (test_assert_corr_symmetric_behavior, test_assert_corr_symmetric_f2py_nocrash,
              test_sfc_varnce_check_behavior, test_length_check_behavior,
              test_pdf_closure_check_behavior, test_rad_check_behavior,
              test_invalid_model_arrays_behavior):
        t()
    print("All validation-check ports PASSED")


if __name__ == "__main__":
    main()
