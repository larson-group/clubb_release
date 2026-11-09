"""Check the JAX hole fillers against the Fortran routines: fill_holes_vertical for every
fill_holes_type, and fill_holes_wp2_from_horz_tke.

Run with `bash human_reviewed_tests/run_pytests.sh -jax_f2py_human_reviewed`,
which puts the compiled CLUBB Python API (clubb_f2py) on the import path.
"""

import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_vertical, fill_holes_wp2_from_horz_tke


NZ = 40
THRESHOLD = 1e-3
FILL_HOLES_TYPES = [1, 2, 3, 4, 5, 6]


def _hole_scenario(name):
    """Return one column (NZ levels) for a named hole-filling scenario.

    Levels 15-20 hold the holes. With rho_ds*dz between about 18 and 84, the holes in
    "needs_wide_range" need more mass than the 10 levels around them provide but
    less than the whole range, and those in "insufficient_mass" need more than
    the whole range provides.
    """
    k = np.arange(NZ)
    if name == "no_holes":
        return 0.5 + 0.2 * np.sin(k)
    if name == "local_holes":
        field = 0.5 + 0.2 * np.sin(k)
        field[[15, 18]] = -0.01
        return field
    if name == "needs_wide_range":
        field = np.full(NZ, 0.2)
        field[15:21] = -0.5
        return field
    if name == "insufficient_mass":
        field = np.full(NZ, 0.02)
        field[15:21] = -0.5
        return field
    raise ValueError(name)


def _run_both(field, rho, dz, grid_dir, fill_holes_type):
    """Run JAX and Fortran fill_holes_vertical on 0-based levels 1..NZ-1 (or NZ-2..0 for grid_dir=-1)."""
    import clubb_f2py
    ng = field.shape[0]
    lower, upper = (1, NZ - 1) if grid_dir == 1 else (NZ - 2, 0)
    j = np.asarray(fill_holes_vertical(
        NZ, ng, THRESHOLD, lower, upper, jnp.asarray(dz), jnp.asarray(rho), grid_dir, fill_holes_type,
        jnp.asarray(field),
    ))
    f = np.asarray(clubb_f2py.f2py_fill_holes_vertical(
        THRESHOLD, lower + 1, upper + 1, dz.copy(), rho.copy(), grid_dir, fill_holes_type, field.copy(),
    ))
    return j, f, slice(min(lower, upper), max(lower, upper) + 1)


@pytest.mark.parametrize("fill_holes_type", FILL_HOLES_TYPES)
@pytest.mark.parametrize("grid_dir", [1, -1])
@pytest.mark.parametrize("scenario", ["no_holes", "local_holes", "needs_wide_range", "insufficient_mass"])
def test_fill_holes_vertical_scenarios_f2py(scenario, grid_dir, fill_holes_type):
    """fill_holes_vertical (fill_holes.F90) for every fill_holes_type on four constructed columns:
    - no_holes: nothing below threshold; the field must come back unchanged.
    - local_holes: two shallow holes that the levels next to them can fill.
    - needs_wide_range: a block of deep holes that only a range wider than 10 levels, or the global fill, can fill.
    - insufficient_mass: a block of deep holes the whole range cannot fill (weighted average below threshold).
    The column is mirrored for grid_dir=-1. A second column holds a no-hole field, so filling one column must
    leave the other alone. Checks that JAX matches Fortran to 1e-12, that the mass sum(field*rho*dz) over the
    fill range is conserved, that fillable holes end at or above threshold and that insufficient mass leaves
    holes."""
    column = _hole_scenario(scenario)
    rho_column = np.linspace(0.6, 1.2, NZ)
    dz_column = np.linspace(30.0, 70.0, NZ)
    if grid_dir == -1:
        column, rho_column, dz_column = column[::-1], rho_column[::-1], dz_column[::-1]
    field = np.stack([column, _hole_scenario("no_holes")])
    rho = np.stack([rho_column, rho_column])
    dz = np.stack([dz_column, dz_column])

    j, f, hf_range = _run_both(field, rho, dz, grid_dir, fill_holes_type)

    assert max_abs_diff("fill_holes_vertical", j, f) < 1e-12
    mass_before = np.sum((field * rho * dz)[:, hf_range], axis=1)
    mass_after = np.sum((j * rho * dz)[:, hf_range], axis=1)
    np.testing.assert_allclose(mass_after, mass_before, rtol=1e-12)
    np.testing.assert_array_equal(j[1], field[1])
    if scenario == "no_holes":
        np.testing.assert_array_equal(j, field)
    elif scenario == "insufficient_mass":
        assert np.min(j[0, hf_range]) < THRESHOLD
    else:
        assert np.min(j[0, hf_range]) >= THRESHOLD - 1e-12


@pytest.mark.parametrize("fill_holes_type", FILL_HOLES_TYPES)
def test_fill_holes_vertical_random_f2py(fill_holes_type):
    """fill_holes_vertical for every fill_holes_type on 20 random two-column fields per grid direction, compared
    with the f2py oracle to 1e-12. Each field has scattered shallow holes and a random block of 8 deep holes, so
    hole positions relative to the windows vary between draws."""
    rng = np.random.default_rng(3)
    worst = 0.0
    for grid_dir in (1, -1):
        for _ in range(20):
            field = rng.uniform(-0.05, 1.0, (2, NZ))
            k_hole = rng.integers(1, NZ - 12)
            field[:, k_hole:k_hole + 8] = rng.uniform(-2.0, -0.5, (2, 8))
            rho = rng.uniform(0.3, 1.2, (2, NZ)); dz = rng.uniform(20.0, 120.0, (2, NZ))
            j, f, _ = _run_both(field, rho, dz, grid_dir, fill_holes_type)
            worst = max(worst, float(max_abs_diff("fill_holes_vertical", j, f)))
    assert worst < 1e-12, f"fill_holes_vertical type {fill_holes_type} f2py mismatch {worst:.2e}"


def test_fill_holes_wp2_from_horz_tke_f2py():
    """fill_holes_wp2_from_horz_tke (fill_holes.F90) — the TKE-conserving wp2 hole-fill: where wp2 < threshold and
    up2/vp2 have spare TKE, it borrows from up2/vp2 proportionally to raise wp2. Unlike fill_holes_vertical this is
    a direct proportional borrow (NOT the FP-floor mass_frac redistribution), so it bit-shadows the f2py oracle
    cleanly. The JAX lower_k/upper_k are 0-based (Fortran lower/upper_hf_level = +1)."""
    import clubb_f2py
    rng = np.random.default_rng(4)
    worst = 0.0
    for _ in range(20):
        ng, nzm = 2, 12
        wp2 = rng.uniform(-0.1, 1.0, (ng, nzm))   # negatives / sub-threshold -> holes to fill
        up2 = rng.uniform(0.05, 1.0, (ng, nzm)); vp2 = rng.uniform(0.05, 1.0, (ng, nzm)); thr = 1e-3
        j = fill_holes_wp2_from_horz_tke(
            nzm, ng, thr, 0, nzm - 1, jnp.asarray(wp2), jnp.asarray(up2), jnp.asarray(vp2),
        )
        f = clubb_f2py.f2py_fill_holes_wp2_from_horz_tke(thr, 1, nzm, np.asfortranarray(wp2.copy()),
                                                         np.asfortranarray(up2.copy()), np.asfortranarray(vp2.copy()))
        for k, out_name in enumerate(("wp2", "up2", "vp2")):
            worst = max(worst, float(max_abs_diff(f"fill_holes_wp2_from_horz_tke {out_name}",
                                                  np.asarray(j[k]), np.asarray(f[k]))))
    assert worst < 1e-12, f"fill_holes_wp2_from_horz_tke f2py mismatch {worst:.2e}"
