"""Check the JAX interpolation.py ports (lin_interpolate_two_points, mono_cubic_interp,
zlinterp_fnc, plinterp_fnc) against Fortran."""

import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.interpolation import (
    lin_interpolate_two_points,
    mono_cubic_interp,
    zlinterp_fnc,
    plinterp_fnc,
)


_CONFIGS = [(0, 0, 1, 2), (0, 1, 2, 2), (0, 1, 2, 3), (2, 1, 2, 3)]
_Z = (0.0, 100.0, 250.0, 450.0)        # zm1, z00, zp1, zp2 (monotone increasing)
_F = (1.0, 2.5, 3.2, 3.9)              # fm1, f00, fp1, fp2 (monotone increasing)


def test_f2py_oracle():
    import clubb_f2py
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
            worst = max(worst, max_abs_diff("mono_cubic_interp", got, ref))
    assert worst < 1e-11, f"mono_cubic_interp f2py mismatch {worst:.2e}"
    # zlinterp_fnc (vertical linear interp + zero-fill outside range) on general source/target grids,
    # including targets outside the source range.
    rng = np.random.default_rng(3)
    zworst = 0.0
    for _ in range(20):
        ns, no = int(rng.integers(4, 12)), int(rng.integers(4, 15))
        grid_src = np.sort(rng.uniform(0.0, 3000.0, ns))
        var_src = rng.uniform(-5.0, 5.0, ns)
        grid_out = np.sort(rng.uniform(-200.0, 3200.0, no))   # spans beyond src range → out-of-range zero-extrap
        gj = np.asarray(zlinterp_fnc(grid_out, grid_src, var_src))
        gf = np.asarray(clubb_f2py.f2py_zlinterp_fnc(grid_out, grid_src, var_src))
        zworst = max(zworst, float(max_abs_diff("interpolation", gj, gf)))
    assert zworst < 1e-12, f"zlinterp_fnc f2py mismatch {zworst:.2e}"


def test_plinterp_f2py():
    """plinterp_fnc (linear interpolation in pressure coordinates, zero outside the source range) matches Fortran."""
    import clubb_f2py
    rng = np.random.default_rng(525)
    worst = 0.0
    for _ in range(40):
        n_src = rng.integers(4, 20)
        grid_src = np.sort(rng.uniform(2.0e4, 1.0e5, n_src))[::-1].copy()   # decreasing pressure [Pa]
        var_src = rng.standard_normal(n_src)
        # Decreasing target pressures spanning in-range, below and above the source grid.
        grid_out = np.sort(rng.uniform(1.5e4, 1.05e5, rng.integers(4, 25)))[::-1].copy()
        got = np.asarray(plinterp_fnc(grid_out, grid_src, var_src))
        ref = np.asarray(clubb_f2py.f2py_plinterp_fnc(grid_out, grid_src, var_src))
        worst = max(worst, float(max_abs_diff("plinterp_fnc", got, ref)))
    assert worst < 1e-12, f"plinterp_fnc f2py mismatch {worst:.2e}"
