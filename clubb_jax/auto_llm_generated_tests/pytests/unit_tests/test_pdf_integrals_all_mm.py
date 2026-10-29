#!/usr/bin/env python3
"""Check mixed-moment PDF probability limits and normal central moments."""
import math
import pytest
from clubb_jax.src.Microphys.KK_microphys import parabolic_cylinder


import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.Microphys.KK_microphys.PDF_integrals_all_MM import trivar_NNL_MM, quadrivar_NNLL_MM


def _phi(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


# Default series accuracy is 1e-4; high accuracy preserves the strict analytic check.
@pytest.mark.parametrize('high_accuracy,tolerance', [(False, 1e-4), (True, 1e-9)])
def test_base_case_probability_mass(monkeypatch, high_accuracy, tolerance):
    """a=b=0 → Phi(mu_x2/sigma_x2), independent of alpha/beta/the other moments."""
    # The series stopping accuracy is selected globally and captured when JIT traces.
    try:
        with monkeypatch.context() as context:
            context.setattr(parabolic_cylinder, 'l_high_accuracy_parab_cyl_fnc', high_accuracy)
            parabolic_cylinder.dv_parabolic_cylinder.clear_cache()
            for mu_x2, sigma_x2 in ((0.3, 0.5), (-0.2, 0.8), (1.0, 0.4), (0.0, 1.0)):
                got = float(trivar_NNL_MM(
                    mu_x1=0.5, mu_x2=mu_x2, mu_x3_n=-0.4, sigma_x1=0.6, sigma_x2=sigma_x2, sigma_x3_n=0.5,
                    rho_x1x2=0.3, rho_x1x3_n=0.2, rho_x2x3_n=-0.1, x1_mean=0.4, x2_alpha_x3_beta_mean=1.3,
                    alpha_exp=1.5, beta_exp=0.5, a_exp=0, b_exp=0))
                ref = _phi(mu_x2 / sigma_x2)
                rel = abs(got - ref) / (abs(ref) + 1e-30)
                assert rel < tolerance, f"base case mu_x2={mu_x2}: rel {rel:.2e} (got {got}, Phi {ref})"
            print(f"  probability mass: high_accuracy={high_accuracy}, tolerance={tolerance:g}  PASS")
    finally:
        parabolic_cylinder.dv_parabolic_cylinder.clear_cache()



def _normal_central(mu, sigma, c, n):
    """n-th central moment E[(X-c)^n] for X~N(mu,sigma) (independent: raw double-factorial expansion)."""
    cm = lambda k: 0.0 if k % 2 else sigma ** k * math.factorial(k) / (2 ** (k // 2) * math.factorial(k // 2))
    return sum(math.comb(n, k) * (mu - c) ** (n - k) * cm(k) for k in range(n + 1))


def _Q():
    return dict(mu_x1=0.5, mu_x2=-0.3, mu_x3_n=-0.4, mu_x4_n=-0.6, s_x1=0.6, s_x2=0.5, s_x3n=0.45, s_x4n=0.5,
                r12=0.25, r13n=0.2, r14n=-0.1, r23n=-0.15, r24n=0.1, r34n=0.3,
                x1_mean=0.4, M=1.3, alpha=1.5, beta=0.5, gamma=0.4)


def _call_quad(q, a, b, **over):
    p = {**q, **over}
    return float(quadrivar_NNLL_MM(
        p['mu_x1'], p['mu_x2'], p['mu_x3_n'], p['mu_x4_n'], p['s_x1'], p['s_x2'], p['s_x3n'], p['s_x4n'],
        p['r12'], p['r13n'], p['r14n'], p['r23n'], p['r24n'], p['r34n'], p['x1_mean'], p['M'],
        p['alpha'], p['beta'], p['gamma'], a, b))


def test_quadrivar_base_and_psum():
    q = _Q()
    # (1) a=b=0 -> Phi(-mu_x2/sigma_x2) (x2<0 mass), independent of all else.
    for mx2, sx2 in ((-0.3, 0.5), (0.2, 0.8), (0.0, 1.0)):
        got = _call_quad(q, 0, 0, mu_x2=mx2, s_x2=sx2)
        ref = _phi(-mx2 / sx2)
        assert abs(got - ref) / (abs(ref) + 1e-30) < 1e-9, f"quad base mu_x2={mx2}"
    # (2) a=2,b=0,rho12=0 -> (2nd central moment of x1 about x1_mean) * Phi(-mu_x2/sigma_x2)
    for a in (1, 2, 3):
        got = _call_quad(q, a, 0, r12=0.0)
        ref = _normal_central(q['mu_x1'], q['s_x1'], q['x1_mean'], a) * _phi(-q['mu_x2'] / q['s_x2'])
        assert abs(got - ref) / (abs(ref) + 1e-30) < 1e-9, f"quad p-sum a={a}"
    print("  quadrivar_NNLL_MM: a=b=0 -> Phi(-mu_x2/sigma_x2) + (a,b=0,rho12=0) -> central-moment*Phi: <1e-9  PASS")
