"""Check the cos_solar_zen_module.F90 port (cosine of the solar zenith angle)
against Fortran."""

import numpy as np

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.Radiation.cos_solar_zen_module import cos_solar_zen


def test_f2py_oracle():
    import clubb_f2py
    rng = np.random.default_rng(3)
    worst = 0.0
    n = 0
    for _ in range(60):
        d = int(rng.integers(1, 28)); m = int(rng.integers(1, 13)); y = int(rng.choice([2000, 2007, 2008, 2023]))
        t = float(rng.uniform(0.0, 3 * 86400.0))   # spans multiple days (exercises the date rollover)
        lat = float(rng.uniform(-80.0, 80.0)); lon = float(rng.uniform(-180.0, 180.0))
        j = float(cos_solar_zen(d, m, y, t, lat, lon))
        f = float(clubb_f2py.f2py_cos_solar_zen(d, m, y, t, lat, lon))
        worst = max(worst, max_abs_diff("cos_solar_zen", j, f)); n += 1
    assert worst < 1e-12, f"cos_solar_zen f2py mismatch {worst:.2e}"
