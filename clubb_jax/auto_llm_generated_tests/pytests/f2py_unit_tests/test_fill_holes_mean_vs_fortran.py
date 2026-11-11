"""Check the MEAN-FIELD vertical hole-filler (the rtm_cl / thlm_cl mechanism)
against Fortran."""

import jax.numpy as jnp
import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_vertical, fill_holes_wp2_from_horz_tke


def test_fill_holes_vertical_f2py():
    """fill_holes_vertical (fill_holes.F90) — the mass-conserving vertical hole-filler, validated vs the f2py oracle
    for BOTH fill_holes_type=1 (single global pass) and 2 (sliding 3-pt window + global fallback). The JAX lower_k/
    upper_k are 0-based (Fortran lower/upper_hf_level = +1). Requires float64 (jax_enable_x64) to match."""
    import clubb_f2py
    rng = np.random.default_rng(0)
    worst = 0.0
    for fht in (1, 2):
        for _ in range(20):
            ng, nz = 2, 16
            field = rng.uniform(-0.05, 1.0, (ng, nz))   # negatives = holes to fill
            rho = rng.uniform(0.3, 1.2, (ng, nz)); dz = rng.uniform(20.0, 120.0, (ng, nz)); thr = 1e-3
            j = np.asarray(fill_holes_vertical(
                nz, ng, thr, 0, nz - 1, jnp.asarray(dz), jnp.asarray(rho), 1, fht, jnp.asarray(field),
            ))
            f = np.asarray(clubb_f2py.f2py_fill_holes_vertical(thr, 1, nz, dz.copy(), rho.copy(), 1, fht, field.copy()))
            worst = max(worst, float(max_abs_diff("fill_holes_mean", j, f)))
    assert worst < 1e-12, f"fill_holes_vertical f2py mismatch {worst:.2e}"


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
        for k in range(3):   # wp2, up2, vp2
            worst = max(worst, float(max_abs_diff("fill_holes_mean", np.asarray(j[k]), np.asarray(f[k]))))
    assert worst < 1e-12, f"fill_holes_wp2_from_horz_tke f2py mismatch {worst:.2e}"
