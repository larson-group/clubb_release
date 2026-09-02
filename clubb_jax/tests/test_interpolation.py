#!/usr/bin/env python3
"""test_interpolation.py — validate the JAX interpolation.py port (lin_interpolate_two_points, mono_cubic_interp).

Oracles:
  1. f2py bit-shadow: clubb_f2py.f2py_lin_interpolate_two_points and f2py_mono_cubic_interp on the same args
     (the compiled default uses the Steffen cubic, l_quintic_poly_interp=False). SKIPs if clubb_f2py is unbuilt.
  2. Closed-form linear interpolation identity.
  3. Steffen monotonicity: the cubic stays within [f00, fp1] for monotone data on [z00, zp1].
  4. A finite jax.grad.
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

from clubb_jax.src.CLUBB_core.interpolation import (
    lin_interpolate_two_points, mono_cubic_interp, linear_interp_factor, zlinterp_fnc,
    lin_interp_between_grids, plinterp_fnc)


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


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py interpolation oracle: SKIP ({type(e).__name__})")
        return
    # lin_interpolate_two_points
    r = float(clubb_f2py.f2py_lin_interpolate_two_points(150.0, 200.0, 100.0, 5.0, 1.0))
    j = float(lin_interpolate_two_points(150.0, 200.0, 100.0, 5.0, 1.0))
    assert abs(j - r) < 1e-13, f"lin_interp f2py mismatch {abs(j-r):.2e}"
    # mono_cubic_interp over all branch configs + a few interpolation altitudes.
    worst = 0.0
    zm1, z00, zp1, zp2 = _Z
    fm1, f00, fp1, fp2 = _F
    for (km1, k00, kp1, kp2) in _CONFIGS:
        for z_in in (120.0, 175.0, 240.0):
            ref = float(clubb_f2py.f2py_mono_cubic_interp(
                z_in, km1, k00, kp1, kp2, zm1, z00, zp1, zp2, fm1, f00, fp1, fp2))
            got = float(mono_cubic_interp(
                z_in, km1, k00, kp1, kp2, zm1, z00, zp1, zp2, fm1, f00, fp1, fp2))
            worst = max(worst, abs(got - ref))
    assert worst < 1e-11, f"mono_cubic_interp f2py mismatch {worst:.2e}"
    # zlinterp_fnc (vertical linear interp + zero-fill outside range) — was only checked vs a literal
    # re-implementation (_zlinterp_ref, which could share a bug with the JAX); the f2py oracle is the independent
    # Fortran ground truth (iter 410). Exercise general source/target grids incl. out-of-range (zero-extrap).
    rng = np.random.default_rng(3)
    zworst = 0.0
    for _ in range(20):
        ns, no = int(rng.integers(4, 12)), int(rng.integers(4, 15))
        grid_src = np.sort(rng.uniform(0.0, 3000.0, ns))
        var_src = rng.uniform(-5.0, 5.0, ns)
        grid_out = np.sort(rng.uniform(-200.0, 3200.0, no))   # spans beyond src range → out-of-range zero-extrap
        gj = np.asarray(zlinterp_fnc(grid_out, grid_src, var_src))
        gf = np.asarray(clubb_f2py.f2py_zlinterp_fnc(grid_out, grid_src, var_src))
        zworst = max(zworst, float(np.max(np.abs(gj - gf))))
    assert zworst < 1e-12, f"zlinterp_fnc f2py mismatch {zworst:.2e}"
    print(f"  f2py lin + mono_cubic_interp + zlinterp_fnc (20 cases, incl. out-of-range): bit-match, "
          f"worst {max(worst, zworst):.2e}  PASS")


def test_monotonicity():
    # Steffen's method keeps the interpolant within [f00, fp1] for monotone data between z00 and zp1.
    zm1, z00, zp1, zp2 = _Z
    fm1, f00, fp1, fp2 = _F
    for km1, k00, kp1, kp2 in ((0, 1, 2, 3),):
        for z_in in np.linspace(z00, zp1, 21):
            v = float(mono_cubic_interp(z_in, km1, k00, kp1, kp2, zm1, z00, zp1, zp2, fm1, f00, fp1, fp2))
            assert f00 - 1e-12 <= v <= fp1 + 1e-12, f"non-monotone at z={z_in}: {v}"
    print("  Steffen monotonicity: interpolant stays within [f00, fp1]  PASS")


def test_differentiable():
    zm1, z00, zp1, zp2 = _Z
    def loss(f):
        fm1, f00, fp1, fp2 = f
        return mono_cubic_interp(175.0, 0, 1, 2, 3, zm1, z00, zp1, zp2, fm1, f00, fp1, fp2)
    g = np.asarray(jax.grad(loss)(jnp.array(_F)))
    assert np.isfinite(g).all(), "non-finite grad through mono_cubic_interp"
    print(f"  jax.grad through mono_cubic_interp: finite ({g.size} entries)  PASS")


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


def test_zlinterp_differentiable():
    grid_src = np.linspace(0.0, 1000.0, 12)
    grid_out = np.linspace(50.0, 950.0, 20)
    def loss(v):
        return jnp.sum(zlinterp_fnc(grid_out, grid_src, v) ** 2)
    g = np.asarray(jax.grad(loss)(jnp.asarray(np.random.default_rng(1).standard_normal(12))))
    assert np.isfinite(g).all(), "non-finite grad through zlinterp_fnc"
    print(f"  jax.grad through zlinterp_fnc: finite ({g.size} entries)  PASS")


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


def test_plinterp():
    """plinterp_fnc (interpolation.F90:plinterp_fnc) does linear interp in PRESSURE coords — it is exactly
    zlinterp_fnc with BOTH grids negated (so a decreasing pressure grid becomes increasing, with the same
    zero-fill-outside-range semantics). It is a ported mirror with NO live-path caller and NO f2py wrapper of
    its own, so it was previously unvalidated. Pin it three ways: (1) the negated-grid identity vs the literal
    _zlinterp_ref, (2) an f2py bit-shadow via f2py_zlinterp_fnc on the negated grids (the documented identity,
    SKIP if unbuilt), (3) decreasing-pressure inputs (the actual use case). (iter 525)"""
    rng = np.random.default_rng(525)
    # Pressure-like decreasing source grid (Pa) + an output grid spanning in-range + below + above.
    worst_ref = 0.0
    f2py = None
    try:
        import clubb_f2py
        f2py = clubb_f2py
    except Exception:
        pass
    worst_f2py = 0.0
    for _ in range(40):
        n_src = rng.integers(4, 20)
        grid_src = np.sort(rng.uniform(2.0e4, 1.0e5, n_src))[::-1].copy()   # decreasing pressure [Pa]
        var_src = rng.standard_normal(n_src)
        # monotonically-decreasing output pressure grid (the real use case; spans in-range + below + above).
        # zlinterp/plinterp assume a monotone target grid (they break on the first out-of-range-above point).
        grid_out = np.sort(rng.uniform(1.5e4, 1.05e5, rng.integers(4, 25)))[::-1].copy()
        got = np.asarray(plinterp_fnc(grid_out, grid_src, var_src))
        # (1) negated-grid identity vs the literal zlinterp reference
        ref = _zlinterp_ref(-grid_out, -grid_src, var_src)
        worst_ref = max(worst_ref, float(np.max(np.abs(got - ref))))
        # (2) f2py bit-shadow on the negated grids
        if f2py is not None:
            gf = np.asarray(f2py.f2py_zlinterp_fnc(-grid_out, -grid_src, var_src))
            worst_f2py = max(worst_f2py, float(np.max(np.abs(got - gf))))
    assert worst_ref < 1e-12, f"plinterp_fnc vs negated-grid zlinterp reference mismatch {worst_ref:.2e}"
    if f2py is not None:
        assert worst_f2py < 1e-12, f"plinterp_fnc f2py bit-shadow mismatch {worst_f2py:.2e}"
        print(f"  plinterp_fnc: negated-grid identity + f2py_zlinterp bit-shadow (40 decreasing-P cases), "
              f"worst {max(worst_ref, worst_f2py):.1e}  PASS")
    else:
        print(f"  plinterp_fnc: negated-grid identity vs literal zlinterp ref (40 cases), worst {worst_ref:.1e} "
              "(f2py SKIP)  PASS")


def test_plinterp_differentiable():
    grid_out = jnp.asarray([9.5e4, 8.0e4, 6.0e4, 4.0e4])
    grid_src = jnp.asarray([1.0e5, 8.5e4, 7.0e4, 5.0e4, 3.0e4])    # decreasing pressure
    v0 = jnp.asarray([1.0, 2.0, 1.5, 0.5, -0.5])
    g = jax.grad(lambda v: jnp.sum(plinterp_fnc(grid_out, grid_src, v) ** 2))(v0)
    assert np.isfinite(np.asarray(g)).all(), "non-finite grad through plinterp_fnc"
    print(f"  jax.grad through plinterp_fnc: finite ({g.size} entries)  PASS")


def main():
    print("test_interpolation:")
    for t in (test_lin_interp_identity, test_f2py_oracle, test_monotonicity, test_differentiable,
              test_linear_interp_factor, test_zlinterp, test_zlinterp_differentiable,
              test_lin_interp_between_grids, test_plinterp, test_plinterp_differentiable):
        t()
    print("All interpolation checks PASSED")


if __name__ == "__main__":
    main()
