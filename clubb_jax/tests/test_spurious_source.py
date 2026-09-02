#!/usr/bin/env python3
"""test_spurious_source.py — validate the JAX calculate_spurious_source port (advance_clubb_core_module).

spurious_source = (integral_after - integral_before)/dt + flux_top - flux_sfc - integral_forcing — the
conservation-error budget diagnostic. Ported but lacking a dedicated f2py test. Oracles:
  1. f2py bit-shadow vs f2py_calculate_spurious_source. SKIPs if clubb_f2py is unbuilt.
  2. Conservation identity: when the column-integral change exactly equals the boundary fluxes plus forcing,
     the spurious source is zero.
  3. A finite jax.grad.
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
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.numerical_check import calculate_spurious_source


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py calculate_spurious_source oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(4)
    worst = 0.0
    dt = 60.0
    for _ in range(200):
        ia, ib, ft, fs, ifc = rng.uniform(-10, 10, 5)
        ref = float(clubb_f2py.f2py_calculate_spurious_source(ia, ib, ft, fs, ifc, dt))
        got = float(calculate_spurious_source(ia, ib, ft, fs, ifc, dt))
        worst = max(worst, abs(got - ref))
    assert worst < 1e-12, f"spurious_source f2py mismatch {worst:.2e}"
    print(f"  f2py calculate_spurious_source: bit-match over 200 configs, worst {worst:.2e}  PASS")


def test_conservation_identity():
    dt = 30.0
    ib, ft, fs, ifc = 5.0, 0.3, 0.1, 0.02
    # Construct integral_after so the budget closes exactly -> spurious source = 0.
    ia = ib + dt * (ifc - ft + fs)
    s = float(calculate_spurious_source(ia, ib, ft, fs, ifc, dt))
    assert abs(s) < 1e-12, f"closed budget should give zero spurious source, got {s}"
    # A perturbation of integral_after by eps -> spurious source = eps/dt.
    s2 = float(calculate_spurious_source(ia + 0.5, ib, ft, fs, ifc, dt))
    assert abs(s2 - 0.5 / dt) < 1e-12, "spurious source sensitivity to integral_after"
    print("  conservation identity: closed budget -> 0; d/d(integral_after) = 1/dt  PASS")


def test_differentiable():
    def loss(args):
        ia, ib, ft, fs, ifc = args
        return calculate_spurious_source(ia, ib, ft, fs, ifc, 60.0) ** 2
    g = np.asarray(jax.grad(loss)(jnp.array([3.0, 2.0, 0.1, 0.05, 0.01])))
    assert np.isfinite(g).all(), "non-finite grad through calculate_spurious_source"
    print(f"  jax.grad through calculate_spurious_source: finite ({g.size} entries)  PASS")


def main():
    print("test_spurious_source:")
    for t in (test_f2py_oracle, test_conservation_identity, test_differentiable):
        t()
    print("All spurious_source checks PASSED")


if __name__ == "__main__":
    main()
