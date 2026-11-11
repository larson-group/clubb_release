#!/usr/bin/env python3
"""validate the cos_solar_zen_module.F90 port (cosine of the solar zenith angle)."""


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
