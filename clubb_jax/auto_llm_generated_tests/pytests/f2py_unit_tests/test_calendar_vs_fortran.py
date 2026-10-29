"""Check the calendar.F90 port (Fliegel & van Flandern Julian-Day-Number routines)
against Fortran."""

from clubb_jax.src.CLUBB_core.calendar import (
    compute_current_date_api,
    gregorian2julian_date,
    julian2gregorian_date,
    gregorian2julian_day,
    leap_year,
)

_DATES = [(1, 1, 1990), (29, 2, 2000), (28, 2, 1900), (31, 12, 2024),
          (15, 7, 2008), (1, 3, 2000), (30, 6, 2023), (1, 1, 2008)]


def test_f2py_oracle():
    """Validate the calendar ports against the actual Fortran oracle (f2py exposes compute_current_date,
    gregorian2julian_day, julian2gregorian_date, leap_year) — the independent ground truth."""
    import clubb_f2py
    times = [0.0, 100.0, 86400.0, 86400 * 45 + 3600.5, 86400 * 400 + 12345.0, 86400 * 366 + 0.25]
    n = 0
    for d, m, y in _DATES:
        assert bool(leap_year(y)) == bool(clubb_f2py.f2py_leap_year(y)), f"leap_year {y}"
        assert int(gregorian2julian_day(d, m, y)) == int(clubb_f2py.f2py_gregorian2julian_day(d, m, y)), f"g2jd {(d, m, y)}"
        jd = gregorian2julian_date(d, m, y)
        assert tuple(julian2gregorian_date(jd)) == tuple(int(v) for v in clubb_f2py.f2py_julian2gregorian_date(jd)), f"j2g {jd}"
        for secs in times:
            n += 1
            jx = compute_current_date_api(d, m, y, secs)
            fp = clubb_f2py.f2py_compute_current_date(d, m, y, secs)
            assert (jx[0], jx[1], jx[2]) == (int(fp[0]), int(fp[1]), int(fp[2])) and abs(jx[3] - fp[3]) <= 1e-9, \
                f"compute_current_date {(d, m, y, secs)}: jax={jx} f2py={tuple(fp)}"
