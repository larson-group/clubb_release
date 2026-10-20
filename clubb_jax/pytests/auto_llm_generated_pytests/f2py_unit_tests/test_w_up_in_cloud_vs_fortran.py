"""Check the JAX calc_w_up_in_cloud port against Fortran."""

import numpy as np

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.pdf_closure_module import calc_w_up_in_cloud

NG, NZ = 2, 8


def _fields(seed):
    rng = np.random.default_rng(seed)
    a = rng.uniform(0.2, 0.8, (NG, NZ))
    cf1 = rng.uniform(0.0, 1.0, (NG, NZ)); cf2 = rng.uniform(0.0, 1.0, (NG, NZ))
    w1 = rng.uniform(-3, 3, (NG, NZ)); w2 = rng.uniform(-3, 3, (NG, NZ))
    v1 = rng.uniform(0.01, 1.0, (NG, NZ)); v2 = rng.uniform(0.01, 1.0, (NG, NZ))
    return a, cf1, cf2, w1, w2, v1, v2


def test_f2py_oracle():
    import clubb_f2py
    worst = 0.0
    for seed in (11, 22, 33):
        a, cf1, cf2, w1, w2, v1, v2 = _fields(seed)
        f = clubb_f2py.f2py_calc_w_up_in_cloud(a, cf1, cf2, w1, w2, v1, v2)
        g = calc_w_up_in_cloud(a, cf1, cf2, w1, w2, v1, v2)
        for fi, gi in zip(f, g):
            worst = max(worst, max_abs_diff("w_up_in_cloud", np.asarray(gi), np.asarray(fi)))
    # ~1e-11 residual is the jax.scipy.special.erf vs Fortran intrinsic erf implementation difference
    # (amplified by the 1/cloudy_updraft_frac division), not a transcription error — closed-form faithful.
    assert worst < 1e-9, f"f2py mismatch {worst:.2e}"
