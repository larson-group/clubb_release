"""Check the JAX calc_comp_corrs_binormal + smooth_corr_quotient ports
against Fortran."""

import jax.numpy as jnp
import numpy as np

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.pdf_utilities import calc_comp_corrs_binormal, smooth_corr_quotient

NG, NZ = 2, 6


def _fields(seed):
    rng = np.random.default_rng(seed)
    mu_x_1 = rng.uniform(-2, 2, (NG, NZ)); mu_x_2 = rng.uniform(-2, 2, (NG, NZ))
    mu_y_1 = rng.uniform(-2, 2, (NG, NZ)); mu_y_2 = rng.uniform(-2, 2, (NG, NZ))
    sx1 = rng.uniform(0.1, 2.0, (NG, NZ)); sx2 = rng.uniform(0.1, 2.0, (NG, NZ))
    sy1 = rng.uniform(0.1, 2.0, (NG, NZ)); sy2 = rng.uniform(0.1, 2.0, (NG, NZ))
    a = rng.uniform(0.2, 0.8, (NG, NZ))
    xm = a * mu_x_1 + (1 - a) * mu_x_2
    ym = a * mu_y_1 + (1 - a) * mu_y_2
    xpyp = rng.uniform(-1, 1, (NG, NZ))
    return dict(xpyp=xpyp, xm=xm, ym=ym, mu_x_1=mu_x_1, mu_x_2=mu_x_2, mu_y_1=mu_y_1, mu_y_2=mu_y_2,
                sigma_x_1_sqd=sx1, sigma_x_2_sqd=sx2, sigma_y_1_sqd=sy1, sigma_y_2_sqd=sy2, mixt_frac=a)


def test_f2py_oracle():
    import clubb_f2py
    worst = 0.0
    for seed in (11, 22, 33):
        f = _fields(seed)
        r1, r2 = clubb_f2py.f2py_calc_comp_corrs_binormal(
            f['xpyp'], f['xm'], f['ym'], f['mu_x_1'], f['mu_x_2'], f['mu_y_1'], f['mu_y_2'],
            f['sigma_x_1_sqd'], f['sigma_x_2_sqd'], f['sigma_y_1_sqd'], f['sigma_y_2_sqd'], f['mixt_frac'])
        g1, g2 = calc_comp_corrs_binormal(**f)
        d = max(max_abs_diff("calc_comp_corrs_binormal", np.asarray(g1), np.asarray(r1)),
                max_abs_diff("calc_comp_corrs_binormal", np.asarray(g2), np.asarray(r2)))
        worst = max(worst, d)
        assert d < 1e-12, f"f2py mismatch seed {seed}: {d:.2e}"


def test_smooth_corr_quotient_f2py():
    """smooth_corr_quotient (pdf_utilities — the smoothly |corr|<=max_mag clamped numerator/denominator) was only
    checked by the bound property above; validate it directly against the f2py oracle on general num/denom."""
    import clubb_f2py
    rng = np.random.default_rng(4)
    worst = 0.0
    for _ in range(40):
        ng, nz = 2, 8
        num = rng.uniform(-1e3, 1e3, (ng, nz)); den = rng.uniform(1e-4, 1e2, (ng, nz)); thr = 1e-10
        j = np.asarray(smooth_corr_quotient(jnp.asarray(num), jnp.asarray(den), thr))
        f = np.asarray(clubb_f2py.f2py_smooth_corr_quotient(num, den, thr))
        worst = max(worst, float(max_abs_diff("calc_comp_corrs_binormal", j, f)))
    assert worst < 1e-12, f"smooth_corr_quotient f2py mismatch {worst:.2e}"
