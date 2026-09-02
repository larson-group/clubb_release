#!/usr/bin/env python3
"""test_hydromet_pdf_parameter.py — validate the hydromet-PDF parameter containers.

This Fortran module is pure data-container init (no physics / no f2py oracle): the only "correct" behavior is
all-zero fields with the right shapes and dims metadata. Tests assert exactly that, plus the zero round-trip.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np

from clubb_jax.src.CLUBB_core.hydromet_pdf_parameter_module import (
    MAX_HYDROMET_DIM, init_hydromet_pdf_params, init_precip_fracs, zero_precip_fracs)


def test_hydromet_pdf_params_zero():
    p = init_hydromet_pdf_params()
    vecs = [p.hm_1, p.hm_2, p.mu_hm_1, p.mu_hm_2, p.sigma_hm_1, p.sigma_hm_2,
            p.corr_w_hm_1, p.corr_w_hm_2, p.corr_chi_hm_1, p.corr_chi_hm_2,
            p.corr_eta_hm_1, p.corr_eta_hm_2]
    for v in vecs:
        assert v.shape == (MAX_HYDROMET_DIM,) and np.all(np.asarray(v) == 0.0)
    for m in (p.corr_hmx_hmy_1, p.corr_hmx_hmy_2):
        assert m.shape == (MAX_HYDROMET_DIM, MAX_HYDROMET_DIM) and np.all(np.asarray(m) == 0.0)
    assert p.mu_Ncn_1 == 0.0 and p.mu_Ncn_2 == 0.0 and p.sigma_Ncn_1 == 0.0 and p.sigma_Ncn_2 == 0.0
    print(f"  init_hydromet_pdf_params: all {len(vecs)} vecs + 2 mats + 4 scalars zero, shapes OK  PASS")


def test_precip_fracs_init_and_zero():
    nzt, ngrdcol = 37, 3
    pf = init_precip_fracs(nzt, ngrdcol)
    assert pf.ngrdcol == ngrdcol and pf.nzt == nzt
    for f in (pf.precip_frac, pf.precip_frac_1, pf.precip_frac_2):
        assert f.shape == (ngrdcol, nzt) and np.all(np.asarray(f) == 0.0)
    # Round-trip: zero a non-zero copy.
    import dataclasses
    import jax.numpy as jnp
    dirty = dataclasses.replace(pf, precip_frac=jnp.ones((ngrdcol, nzt)))
    cleaned = zero_precip_fracs(dirty)
    assert cleaned.ngrdcol == ngrdcol and cleaned.nzt == nzt
    assert np.all(np.asarray(cleaned.precip_frac) == 0.0), "zero_precip_fracs failed to zero"
    print(f"  init_precip_fracs/zero_precip_fracs: shape ({ngrdcol},{nzt}), dims + zero round-trip OK  PASS")


def main():
    print("test_hydromet_pdf_parameter:")
    for t in (test_hydromet_pdf_params_zero, test_precip_fracs_init_and_zero):
        t()
    print("All hydromet_pdf_parameter checks PASSED")


if __name__ == "__main__":
    main()
