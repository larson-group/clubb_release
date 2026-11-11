#!/usr/bin/env python3
"""validate the calendar.F90 port (Fliegel & van Flandern Julian-Day-Number routines)."""


from clubb_jax.src.CLUBB_core.calendar import (
    compute_current_date_api, gregorian2julian_date, julian2gregorian_date,
    gregorian2julian_day, leap_year,
)

_DATES = [(1, 1, 1990), (29, 2, 2000), (28, 2, 1900), (31, 12, 2024),
          (15, 7, 2008), (1, 3, 2000), (30, 6, 2023), (1, 1, 2008)]


def _old_compute(day, month, year, current_time_s):
    """Independent month-walking reference for compute_current_date."""
    total_days = int(current_time_s // 86400)
    time_in_day = current_time_s - total_days * 86400.0
    d, m, y = day, month, year
    remaining = total_days
    while remaining > 0:
        dim = [0, 31, 29 if leap_year(y) else 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        days_left = dim[m] - d + 1
        if remaining < days_left:
            d += remaining; remaining = 0
        else:
            remaining -= days_left; m += 1
            if m > 12: m = 1; y += 1
            d = 1
    return d, m, y, time_in_day


def test_jdn_anchors_and_roundtrip():
    # Published anchors (proleptic Gregorian JDN at noon).
    assert gregorian2julian_date(1, 1, 2000) == 2451545, gregorian2julian_date(1, 1, 2000)
    assert gregorian2julian_date(1, 1, 1970) == 2440588, gregorian2julian_date(1, 1, 1970)
    # Exact inverse over a range of dates.
    for d, m, y in _DATES:
        jd = gregorian2julian_date(d, m, y)
        assert julian2gregorian_date(jd) == (d, m, y), f"roundtrip {d}/{m}/{y} -> {jd} -> {julian2gregorian_date(jd)}"
    print("  gregorian2julian_date/julian2gregorian_date: anchors + roundtrip  PASS")


def test_compute_current_date_api_matches_reference():
    times = [0.0, 100.0, 86400.0, 86400 * 45 + 3600.5, 86400 * 400 + 12345.0,
             86400 * 1000.0, 86400 * 366 + 0.25]
    n = 0
    for d, m, y in _DATES:
        for secs in times:
            n += 1
            assert compute_current_date_api(d, m, y, secs) == _old_compute(d, m, y, secs), \
                f"compute_current_date_api mismatch at {(d, m, y, secs)}"
    print(f"  compute_current_date_api: matches month-walking reference over {n} cases  PASS")


def test_gregorian2julian_day():
    assert gregorian2julian_day(1, 1, 2008) == 1
    assert gregorian2julian_day(31, 12, 2008) == 366   # leap year
    assert gregorian2julian_day(1, 3, 2001) == 60      # non-leap: 31+28+1
    print("  gregorian2julian_day: day-of-year  PASS")
