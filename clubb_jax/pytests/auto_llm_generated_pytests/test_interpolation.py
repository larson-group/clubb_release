#!/usr/bin/env python3
"""validate the JAX interpolation.py port (lin_interpolate_two_points, mono_cubic_interp)."""


import numpy as np
import jax
jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.interpolation import lin_interpolate_two_points, mono_cubic_interp, linear_interp_factor, zlinterp_fnc, lin_interp_between_grids


def _binary_search(array, var):
    """Literal transcription of interpolation.F90:binary_search (1-based index of the >= bracket, -1 if out)."""
    n = len(array)
    low, high = 2, n
    if var < array[0] or var > array[n - 1] or n < 2:
        return -1
    if array[0] <= var <= array[1]:
        return 2
    while low <= high:
        i = (low + high) // 2
        if array[i - 2] < var <= array[i - 1]:
            return i
        elif var < array[i - 1]:
            high = i - 1
        else:
            low = i + 1
    return -1


def _zlinterp_ref(grid_out, grid_src, var_src):
    """Literal Fortran zlinterp_fnc (binary_search + lin_interpolate_two_points, zero outside range)."""
    out = np.zeros(len(grid_out))
    for kint, go in enumerate(grid_out):
        if go < grid_src[0]:
            continue
        k = _binary_search(grid_src, go)
        if k == -1:
            break
        km1 = max(1, k - 1)
        out[kint] = ((go - grid_src[km1 - 1]) / (grid_src[k - 1] - grid_src[km1 - 1])
                     * (var_src[k - 1] - var_src[km1 - 1]) + var_src[km1 - 1])
    return out

# Branch configurations: (km1, k00, kp1, kp2) exercising km1==k00 / kp1==kp2 / interior / extrapolate.
_CONFIGS = [(0, 0, 1, 2), (0, 1, 2, 2), (0, 1, 2, 3), (2, 1, 2, 3)]
_Z = (0.0, 100.0, 250.0, 450.0)        # zm1, z00, zp1, zp2 (monotone increasing)
_F = (1.0, 2.5, 3.2, 3.9)              # fm1, f00, fp1, fp2 (monotone increasing)


def test_lin_interp_identity():
    val = float(lin_interpolate_two_points(150.0, 200.0, 100.0, 5.0, 1.0))
    assert abs(val - ((150.0 - 100.0) / (200.0 - 100.0) * (5.0 - 1.0) + 1.0)) < 1e-14
    # Endpoints reproduce the known values.
    assert abs(float(lin_interpolate_two_points(100.0, 200.0, 100.0, 5.0, 1.0)) - 1.0) < 1e-14
    assert abs(float(lin_interpolate_two_points(200.0, 200.0, 100.0, 5.0, 1.0)) - 5.0) < 1e-14
    print("  lin_interpolate_two_points: closed-form + endpoints  PASS")


def test_monotonicity():
    # Steffen's method keeps the interpolant within [f00, fp1] for monotone data between z00 and zp1.
    zm1, z00, zp1, zp2 = _Z
    fm1, f00, fp1, fp2 = _F
    for km1, k00, kp1, kp2 in ((0, 1, 2, 3),):
        for z_in in np.linspace(z00, zp1, 21):
            v = float(mono_cubic_interp(z_in, km1, k00, kp1, kp2, zm1, z00, zp1, zp2, fm1, f00, fp1, fp2))
            assert f00 - 1e-12 <= v <= fp1 + 1e-12, f"non-monotone at z={z_in}: {v}"
    print("  Steffen monotonicity: interpolant stays within [f00, fp1]  PASS")


def test_linear_interp_factor():
    assert abs(float(linear_interp_factor(0.25, 8.0, 4.0)) - (0.25 * (8.0 - 4.0) + 4.0)) < 1e-14
    assert abs(float(linear_interp_factor(0.0, 8.0, 4.0)) - 4.0) < 1e-14
    assert abs(float(linear_interp_factor(1.0, 8.0, 4.0)) - 8.0) < 1e-14
    print("  linear_interp_factor: closed-form + endpoints  PASS")


def test_zlinterp():
    rng = np.random.default_rng(13)
    grid_src = np.sort(rng.uniform(0.0, 10000.0, 30))
    var_src = rng.standard_normal(30)
    grid_out = np.sort(rng.uniform(-500.0, 11000.0, 50))   # straddles both ends -> zero-fill exercised
    got = np.asarray(zlinterp_fnc(grid_out, grid_src, var_src))
    ref = _zlinterp_ref(grid_out, grid_src, var_src)
    assert np.max(np.abs(got - ref)) < 1e-12, f"zlinterp mismatch {np.max(np.abs(got-ref)):.2e}"
    # Zero-fill below/above the source range.
    assert got[grid_out < grid_src[0]].tolist() == [0.0] * int((grid_out < grid_src[0]).sum())
    assert np.all(got[grid_out > grid_src[-1]] == 0.0)
    print("  zlinterp_fnc: matches literal binary_search+lin_interp, zero-fill outside range  PASS")


def _lin_interp_between_grids_ref(interp_alt, cur_alt, cur_val, tol=1e-6):
    """Literal transcription of interpolation.F90:lin_interp_between_grids (per-point search + clamp)."""
    n = len(cur_alt)
    out = np.empty(len(interp_alt))
    for ii, x in enumerate(interp_alt):
        done = False
        k = 0
        while (not done) and k < n:
            if abs(x - cur_alt[k]) < tol:
                out[ii] = cur_val[k]; done = True
            elif x < cur_alt[k]:
                if k > 0:
                    out[ii] = float(lin_interpolate_two_points(x, cur_alt[k], cur_alt[k - 1],
                                                               cur_val[k], cur_val[k - 1]))
                else:
                    out[ii] = cur_val[0]
                done = True
            k += 1
        if not done and k >= n:
            out[ii] = cur_val[-1]
    return out


def test_lin_interp_between_grids():
    rng = np.random.default_rng(21)
    worst = 0.0
    for _ in range(100):
        cur_alt = np.unique(np.sort(rng.uniform(0.0, 5000.0, rng.integers(4, 40))))
        if len(cur_alt) < 2:
            continue
        cur_val = rng.standard_normal(len(cur_alt))
        tgt = np.concatenate([rng.uniform(-200.0, 5200.0, 30), cur_alt])  # in-range + clamp + exact-match
        got = np.asarray(lin_interp_between_grids(tgt, cur_alt, cur_val))
        ref = _lin_interp_between_grids_ref(tgt, cur_alt, cur_val)
        worst = max(worst, float(np.max(np.abs(got - ref))))
    assert worst < 1e-12, f"lin_interp_between_grids mismatch {worst:.2e}"
    print(f"  lin_interp_between_grids: matches literal Fortran loop (clamp+exact-match), worst {worst:.1e}  PASS")
