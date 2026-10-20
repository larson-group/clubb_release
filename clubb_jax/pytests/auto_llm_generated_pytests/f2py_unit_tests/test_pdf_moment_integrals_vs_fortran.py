"""Check the JAX binormal/trinormal PDF moment integrals against Fortran."""

import numpy as np

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.pdf_closure_module import (
    calc_wp4_pdf,
    calc_wp2xp_pdf,
    calc_wpxp2_pdf,
    calc_wpxpyp_pdf,
)

NG, NZ = 2, 5


def _binormal_params(rng, shape):
    return dict(w_1=rng.uniform(-2, 2, shape), w_2=rng.uniform(-2, 2, shape),
                x_1=rng.uniform(-2, 2, shape), x_2=rng.uniform(-2, 2, shape),
                varnce_w_1=rng.uniform(0.1, 1.5, shape), varnce_w_2=rng.uniform(0.1, 1.5, shape),
                varnce_x_1=rng.uniform(0.1, 1.5, shape), varnce_x_2=rng.uniform(0.1, 1.5, shape),
                corr_w_x_1=rng.uniform(-0.8, 0.8, shape), corr_w_x_2=rng.uniform(-0.8, 0.8, shape),
                mixt_frac=rng.uniform(0.3, 0.7, shape))


def test_f2py_oracle():
    import clubb_f2py
    rng = np.random.default_rng(3)
    sh = (NG, NZ)
    worst = 0.0
    for _ in range(15):
        p = _binormal_params(rng, sh)
        wm = p['mixt_frac'] * p['w_1'] + (1 - p['mixt_frac']) * p['w_2']
        xm = p['mixt_frac'] * p['x_1'] + (1 - p['mixt_frac']) * p['x_2']
        # wp4
        r = np.asarray(clubb_f2py.f2py_calc_wp4_pdf(wm, p['w_1'], p['w_2'], p['varnce_w_1'], p['varnce_w_2'], p['mixt_frac']))
        g = np.asarray(calc_wp4_pdf(wm, p['w_1'], p['w_2'], p['varnce_w_1'], p['varnce_w_2'], p['mixt_frac']))
        worst = max(worst, max_abs_diff("pdf_moment_integrals", g, r))
        # wp2xp
        args = (wm, xm, p['w_1'], p['w_2'], p['x_1'], p['x_2'], p['varnce_w_1'], p['varnce_w_2'],
                p['varnce_x_1'], p['varnce_x_2'], p['corr_w_x_1'], p['corr_w_x_2'], p['mixt_frac'])
        worst = max(worst, max_abs_diff("pdf_moment_integrals", np.asarray(calc_wp2xp_pdf(*args)), np.asarray(clubb_f2py.f2py_calc_wp2xp_pdf(*args))))
        worst = max(worst, max_abs_diff("pdf_moment_integrals", np.asarray(calc_wpxp2_pdf(*args)), np.asarray(clubb_f2py.f2py_calc_wpxp2_pdf(*args))))
        # wpxpyp (reuse x params as y for simplicity, distinct corrs)
        ym = xm
        wargs = (wm, xm, ym, p['w_1'], p['w_2'], p['x_1'], p['x_2'], p['x_1'], p['x_2'],
                 p['varnce_w_1'], p['varnce_w_2'], p['varnce_x_1'], p['varnce_x_2'], p['varnce_x_1'], p['varnce_x_2'],
                 p['corr_w_x_1'], p['corr_w_x_2'], p['corr_w_x_1'], p['corr_w_x_2'], p['corr_w_x_1'], p['corr_w_x_2'],
                 p['mixt_frac'])
        worst = max(worst, max_abs_diff("pdf_moment_integrals", np.asarray(calc_wpxpyp_pdf(*wargs)), np.asarray(clubb_f2py.f2py_calc_wpxpyp_pdf(*wargs))))
    assert worst < 1e-11, f"pdf-moment f2py mismatch {worst:.2e}"
