#!/usr/bin/env python3
"""validate the JAX calculate_spurious_source port (advance_clubb_core_module)."""

from clubb_jax.auto_llm_generated_tests.pytests.array_inputs import call_jax


from clubb_jax.src.CLUBB_core.numerical_check import calculate_spurious_source


def test_conservation_identity():
    dt = 30.0
    ib, ft, fs, ifc = 5.0, 0.3, 0.1, 0.02
    # Construct integral_after so the budget closes exactly -> spurious source = 0.
    ia = ib + dt * (ifc - ft + fs)
    s = float(call_jax(calculate_spurious_source, ia, ib, ft, fs, ifc, dt))
    assert abs(s) < 1e-12, f"closed budget should give zero spurious source, got {s}"
    # A perturbation of integral_after by eps -> spurious source = eps/dt.
    s2 = float(call_jax(calculate_spurious_source, ia + 0.5, ib, ft, fs, ifc, dt))
    assert abs(s2 - 0.5 / dt) < 1e-12, "spurious source sensitivity to integral_after"
    print("  conservation identity: closed budget -> 0; d/d(integral_after) = 1/dt  PASS")
