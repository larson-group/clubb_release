#!/usr/bin/env python3
"""test_cos_solar_zen.py — validate the cos_solar_zen_module.F90 port (cosine of the solar zenith angle).

cos_solar_zen builds on the calendar routines (compute_current_date_api / gregorian2julian_day / leap_year, validated
in test_calendar) and adds the solar-declination + hour-angle geometry. test_bugsrad bypasses it with the fixed
cos_solar_zen path (l_fix_cos_solar_zen), so this is its first dedicated oracle test. Oracles:
  1. f2py bit-shadow vs f2py_cos_solar_zen across dates / latitudes / longitudes / multi-day elapsed times. SKIPs
     if clubb_f2py is unbuilt.
  2. Physical bounds: cos(zenith) in [-1, 1].
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

import numpy as np
from clubb_jax.src.Radiation.cos_solar_zen_module import cos_solar_zen

_DATES = [(21, 3, 2008), (21, 6, 2008), (21, 12, 2007), (1, 1, 2000), (15, 7, 2023)]


def test_bounds():
    """cos(zenith) is a cosine, so it must lie in [-1, 1] for every date/time/latitude (negative = sun below
    the horizon — cos_solar_zen returns the raw cosine, not the night-clamped value)."""
    lo, hi = 2.0, -2.0
    for d, m, y in _DATES:
        for t in (0.0, 6 * 3600.0, 12 * 3600.0, 18 * 3600.0):
            for lat in (-60.0, 0.0, 45.0, 80.0):
                cz = float(cos_solar_zen(d, m, y, t, lat, 0.0))
                lo, hi = min(lo, cz), max(hi, cz)
    assert -1.0 - 1e-12 <= lo and hi <= 1.0 + 1e-12, f"cos_zen out of [-1,1]: [{lo}, {hi}]"
    print(f"  cos_solar_zen physical bounds: cos(zenith) in [{lo:.3f}, {hi:.3f}] subset of [-1, 1]  PASS")


def test_f2py_oracle():
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py cos_solar_zen oracle: SKIP ({type(e).__name__})")
        return
    rng = np.random.default_rng(3)
    worst = 0.0
    n = 0
    for _ in range(60):
        d = int(rng.integers(1, 28)); m = int(rng.integers(1, 13)); y = int(rng.choice([2000, 2007, 2008, 2023]))
        t = float(rng.uniform(0.0, 3 * 86400.0))   # spans multiple days (exercises the date rollover)
        lat = float(rng.uniform(-80.0, 80.0)); lon = float(rng.uniform(-180.0, 180.0))
        j = float(cos_solar_zen(d, m, y, t, lat, lon))
        f = float(clubb_f2py.f2py_cos_solar_zen(d, m, y, t, lat, lon))
        worst = max(worst, abs(j - f)); n += 1
    assert worst < 1e-12, f"cos_solar_zen f2py mismatch {worst:.2e}"
    print(f"  cos_solar_zen vs f2py oracle ({n} cases, dates/lat/lon/multi-day t): bit-match, worst {worst:.2e}  PASS")


def main():
    print("test_cos_solar_zen:")
    test_bounds()
    test_f2py_oracle()
    print("All cos_solar_zen checks PASSED")


if __name__ == "__main__":
    main()
