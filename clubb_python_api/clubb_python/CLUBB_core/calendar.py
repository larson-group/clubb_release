"""User-facing wrappers for routines from CLUBB_core/calendar.F90."""

import clubb_f2py


def julian2gregorian_date(julian_date: int):
    """Convert a Julian date to Gregorian day, month, year."""
    return clubb_f2py.f2py_julian2gregorian_date(int(julian_date))


def compute_current_date(previous_day: int, previous_month: int, previous_year: int,
                         seconds_since_previous_date: float):
    """Advance Gregorian date/time by elapsed seconds."""
    return clubb_f2py.f2py_compute_current_date(
        int(previous_day),
        int(previous_month),
        int(previous_year),
        float(seconds_since_previous_date),
    )


def gregorian2julian_day(day: int, month: int, year: int) -> int:
    """Convert Gregorian day/month/year to Julian day-of-year."""
    return int(clubb_f2py.f2py_gregorian2julian_day(int(day), int(month), int(year)))


def leap_year(year: int) -> bool:
    """Return True when the Gregorian year is a leap year."""
    return bool(clubb_f2py.f2py_leap_year(int(year)))
