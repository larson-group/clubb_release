#!/usr/bin/env python3
"""test_time_select.py — validate `time_select`, the time-bracket selector used by the supported per-case
`<case>_read_t_dependent` readers (jun25/nov11) and generic time-dependent forcing.

`time_select(time, times)` (time_dependent_input.py:26 ↔ time_dependent_input.F90:1125) returns the (before, after,
frac) bracket for linear time interpolation. It is the single mechanism underneath all time-dependent surface-flux /
forcing interpolation, but was only exercised transitively (via apply_time_dependent_forcings), never directly — so
its edge cases were unpinned. This test pins:
  1. Strictly-interior times: the 0-based (before, after, frac) triple matches a literal transcription of the Fortran
     do-loop bracket search EXACTLY.
  2. Interpolation EQUIVALENCE at ALL query times incl. exact interior nodes: (1-frac)*v[before] + frac*v[after]
     equals the Fortran result. (At an exact interior node the JAX returns the UPPER bracket with frac=0 while the
     Fortran transcription returns the LOWER bracket with frac=1 — different triple, identical interpolated value;
     this test documents that intentional, benign difference rather than asserting triple-equality there.)
  3. Endpoints (time==times[0] / times[-1]) and out-of-range rejection (the Fortran `error stop`; JAX raises).
Oracle-independent (the Fortran bracket logic is a few lines, transcribed directly); never SKIPs. (iter 532)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np

from clubb_jax.src.Benchmark_cases.time_dependent_input import time_select


def _ref_time_select(time, times):
    """Literal transcription of time_dependent_input.F90:time_select (0-based). Picks the FIRST bracket [k,k+1]
    with times[k] <= time <= times[k+1]; raises on out-of-range (the Fortran `error stop`)."""
    n = len(times)
    if time < times[0] or time > times[-1]:
        raise ValueError("out of range")
    before = None
    for k in range(n - 1):
        assert times[k] < times[k + 1], "times not strictly increasing"
        if before is None and times[k] <= time <= times[k + 1]:
            before = k
    after = before + 1
    frac = (time - times[before]) / (times[after] - times[before])
    return before, after, float(frac)


_TIMES = np.array([0.0, 100.0, 250.0, 250.0 + 300.0, 900.0])  # non-uniform, strictly increasing
_VALS = np.array([2.0, 5.0, -1.0, 3.5, 0.25])                  # a flux profile sampled at _TIMES


def _interp(before, after, frac, vals):
    return (1.0 - frac) * vals[before] + frac * vals[after]


def test_interior_triples_match_fortran():
    worst = 0.0
    for t in (10.0, 50.0, 99.999, 100.001, 175.0, 400.0, 549.0, 700.0, 899.0):
        jb, ja, jf = time_select(float(t), len(_TIMES), _TIMES)
        rb, ra, rf = _ref_time_select(float(t), _TIMES)
        assert (jb, ja) == (rb, ra), f"t={t}: JAX bracket ({jb},{ja}) != Fortran ({rb},{ra})"
        worst = max(worst, abs(jf - rf))
    assert worst < 1e-12, f"interior frac mismatch {worst:.2e}"
    print(f"  interior times: (before,after,frac) match the Fortran bracket loop (worst Δfrac {worst:.1e})  PASS")


def test_interpolation_equivalent_everywhere():
    """The value that actually matters — interpolation result — matches the Fortran at ALL query times, including
    the exact interior nodes where the (before,after,frac) triple differs but the interpolated value does not."""
    worst = 0.0
    queries = list(_TIMES) + [10.0, 175.0, 400.0, 700.0, 899.999]   # incl. every exact node
    for t in queries:
        jval = _interp(*time_select(float(t), len(_TIMES), _TIMES), _VALS)
        rval = _interp(*_ref_time_select(float(t), _TIMES), _VALS)
        worst = max(worst, abs(jval - rval))
    assert worst < 1e-12, f"interpolated-value mismatch {worst:.2e}"
    # And exact-node queries must reproduce the node value exactly.
    for i, t in enumerate(_TIMES):
        assert abs(_interp(*time_select(float(t), len(_TIMES), _TIMES), _VALS) - _VALS[i]) < 1e-12, f"node {i} not exact"
    print(f"  interpolation equivalent to Fortran at all times incl. exact nodes (worst {worst:.1e})  PASS")


def test_endpoints_and_out_of_range():
    b0, a0, f0 = time_select(float(_TIMES[0]), len(_TIMES), _TIMES)
    assert (b0, a0, f0) == (0, 1, 0.0), f"first endpoint triple {b0, a0, f0}"
    bN, aN, fN = time_select(float(_TIMES[-1]), len(_TIMES), _TIMES)
    assert (bN, aN) == (len(_TIMES) - 2, len(_TIMES) - 1) and abs(fN - 1.0) < 1e-12, f"last endpoint {bN, aN, fN}"
    for bad in (_TIMES[0] - 1.0, _TIMES[-1] + 1.0):
        try:
            time_select(float(bad), len(_TIMES), _TIMES)
        except ValueError:
            pass
        else:
            raise AssertionError(f"out-of-range time {bad} was NOT rejected (mirrors the Fortran error stop)")
    print("  endpoints (frac 0 / 1) + out-of-range rejection (mirrors Fortran error stop)  PASS")


def main():
    print("test_time_select:")
    test_interior_triples_match_fortran()
    test_interpolation_equivalent_everywhere()
    test_endpoints_and_out_of_range()
    print("All time_select checks PASSED")


if __name__ == "__main__":
    main()
