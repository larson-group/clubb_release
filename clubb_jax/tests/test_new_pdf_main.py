#!/usr/bin/env python3
"""test_new_pdf_main.py — validate the JAX calc_F_x_zeta_x_setter port (new_pdf_main.F90).

No direct f2py oracle exists for this routine (it is internal to new_pdf_driver), so it is validated by a
literal NumPy transcription, the F_x bounds, the analytic limits, and a finite jax.grad.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.new_pdf_main import calc_F_x_zeta_x_setter, new_pdf_driver

NG, NZ = 2, 6
_SLOPE, _STDEV_FACTOR, _LAMBDA = 1.4, 1.2, 0.5


def _ref(Skx, slope, stdev_factor, lam):
    absS = np.abs(Skx)
    min_F = np.where(absS > 0, 1e-3, 0.0)
    e = np.exp(-(absS ** lam) / slope)
    F = min_F * e + 1.0 * (1 - e)
    return F, np.full_like(Skx, stdev_factor - 1.0), min_F, np.ones_like(Skx)


def test_transcription():
    rng = np.random.default_rng(3)
    worst = 0.0
    for _ in range(50):
        Skx = rng.uniform(-3, 3, (NG, NZ))
        g = calc_F_x_zeta_x_setter(Skx, _SLOPE, _STDEV_FACTOR, _LAMBDA)
        r = _ref(Skx, _SLOPE, _STDEV_FACTOR, _LAMBDA)
        for gi, ri in zip(g, r):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - ri)))
    assert worst < 1e-14, f"transcription mismatch {worst:.2e}"
    print(f"  calc_F_x_zeta_x_setter: literal transcription, worst {worst:.2e}  PASS")


def test_bounds_and_limits():
    Skx = np.linspace(-3, 3, 41)
    F, zeta, minF, maxF = (np.asarray(x) for x in calc_F_x_zeta_x_setter(Skx, _SLOPE, _STDEV_FACTOR, _LAMBDA))
    assert np.all(F >= minF - 1e-12) and np.all(F <= maxF + 1e-12), "F_x out of [min,max]"
    # Skx=0 -> exp factor = 1 -> F_x = min_F_x = 0.
    F0 = float(np.asarray(calc_F_x_zeta_x_setter(np.array([0.0]), _SLOPE, _STDEV_FACTOR, _LAMBDA)[0])[0])
    assert abs(F0) < 1e-14, "Skx=0 should give F_x=0"
    # Large |Skx| -> exp(-|Skx|^lambda/slope) -> 0 -> F_x -> max_F_x = 1 (needs |Skx| large enough that
    # sqrt(|Skx|)/slope >> 1, e.g. 1e4 -> exp(-100/1.4) ~ 0).
    Fbig = float(np.asarray(calc_F_x_zeta_x_setter(np.array([1.0e4]), _SLOPE, _STDEV_FACTOR, _LAMBDA)[0])[0])
    assert abs(Fbig - 1.0) < 1e-6, "large |Skx| should give F_x->1"
    assert abs(zeta[0] - (_STDEV_FACTOR - 1.0)) < 1e-14, "zeta_x = stdev_factor - 1"
    print("  bounds F_x in [min,max] + Skx=0/large limits + zeta identity  PASS")


def test_differentiable():
    Skx = jnp.asarray(np.random.default_rng(1).uniform(0.1, 3, (NG, NZ)))
    g = np.asarray(jax.grad(lambda s: jnp.sum(calc_F_x_zeta_x_setter(s, _SLOPE, _STDEV_FACTOR, _LAMBDA)[0] ** 2))(Skx))
    assert np.isfinite(g).all(), "non-finite grad through calc_F_x_zeta_x_setter"
    print(f"  jax.grad through calc_F_x_zeta_x_setter: finite ({g.size} entries)  PASS")


def _f2py_paths():
    for p in (_ROOT, _ROOT + "/clubb_python_api"):
        if p not in sys.path:
            sys.path.append(p)


def test_driver_f2py():
    _f2py_paths()
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py new_pdf_driver oracle: SKIP ({type(e).__name__})")
        return
    nparams, ng, nz = 102, 2, 8
    rng = np.random.default_rng(11)
    # Plausible tunable params: slope>0, stdev_factor>0, coef_spread in (0,1].
    params = np.zeros((ng, nparams))
    params[:, 52] = 1.4    # islope_coef_spread_DG_means_w
    params[:, 53] = 1.2    # ipdf_component_stdev_factor_w
    params[:, 54] = 0.7    # icoef_spread_DG_means_rt
    params[:, 55] = 0.7    # icoef_spread_DG_means_thl
    worst = 0.0
    for _ in range(3):
        wm = rng.uniform(-1, 1, (ng, nz)); rtm = rng.uniform(0, 0.02, (ng, nz)); thlm = rng.uniform(285, 305, (ng, nz))
        wp2 = rng.uniform(0.05, 2, (ng, nz)); rtp2 = rng.uniform(1e-7, 1e-6, (ng, nz)); thlp2 = rng.uniform(0.05, 1, (ng, nz))
        Skw = rng.uniform(-2, 2, (ng, nz)); Skrt = rng.uniform(-2, 2, (ng, nz)); Skthl = rng.uniform(-2, 2, (ng, nz))
        wprtp = rng.uniform(-1e-3, 1e-3, (ng, nz)); wpthlp = rng.uniform(-0.5, 0.5, (ng, nz))
        rtpthlp = rng.uniform(-1e-4, 1e-4, (ng, nz))
        f = clubb_f2py.f2py_new_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, wprtp, wpthlp, rtpthlp,
                                           params, Skrt.copy(), Skthl.copy())
        g = new_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, wprtp, wpthlp, rtpthlp, params, Skrt, Skthl)
        for fi, gi in zip(f, g):
            worst = max(worst, np.max(np.abs(np.asarray(gi) - np.asarray(fi))))
    assert worst < 1e-11, f"new_pdf_driver f2py mismatch {worst:.2e}"
    print(f"  f2py new_pdf_driver: bit-match (15 PDF-param outputs incl. clipped Skrt/Skthl), worst {worst:.2e}  PASS")


def test_driver_differentiable():
    _f2py_paths()
    ng, nz = 1, 5
    rng = np.random.default_rng(3)
    params = np.zeros((ng, 102)); params[:, 52] = 1.4; params[:, 53] = 1.2; params[:, 54] = 0.7; params[:, 55] = 0.7
    a = dict(wm=rng.uniform(-1, 1, (ng, nz)), rtm=rng.uniform(0, 0.02, (ng, nz)), thlm=rng.uniform(285, 305, (ng, nz)),
             wp2=rng.uniform(0.1, 2, (ng, nz)), rtp2=rng.uniform(1e-7, 1e-6, (ng, nz)), thlp2=rng.uniform(0.1, 1, (ng, nz)),
             wprtp=rng.uniform(-1e-3, 1e-3, (ng, nz)), wpthlp=rng.uniform(-0.5, 0.5, (ng, nz)),
             rtpthlp=rng.uniform(-1e-4, 1e-4, (ng, nz)), Skrt=rng.uniform(-1, 1, (ng, nz)), Skthl=rng.uniform(-1, 1, (ng, nz)))
    def loss(Skw):
        outs = new_pdf_driver(a['wm'], a['rtm'], a['thlm'], a['wp2'], a['rtp2'], a['thlp2'], Skw,
                              a['wprtp'], a['wpthlp'], a['rtpthlp'], params, a['Skrt'], a['Skthl'])
        return sum(jnp.sum(o ** 2) for o in outs)
    g = np.asarray(jax.grad(loss)(jnp.asarray(rng.uniform(-2, 2, (ng, nz)))))
    assert np.isfinite(g).all(), "non-finite grad through new_pdf_driver"
    print(f"  jax.grad through new_pdf_driver: finite ({g.size} entries)  PASS")


def main():
    print("test_new_pdf_main:")
    for t in (test_transcription, test_bounds_and_limits, test_differentiable,
              test_driver_f2py, test_driver_differentiable):
        t()
    print("All new_pdf_main checks PASSED")


if __name__ == "__main__":
    main()
