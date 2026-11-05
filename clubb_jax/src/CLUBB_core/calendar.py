"""JAX-compatible port of ``src/CLUBB_core/calendar.F90``.

Description:
    Gregorian and Julian calendar conversion routines used by solar-zenith
    radiation. Fortran subroutine outputs return tuples in source output order.

JAX adaptation:
    Integer calculations remain JAX arrays so calendar rollover works while
    elapsed model time is traced.
"""

import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import sec_per_day

month_names = ("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
days_per_month = jnp.asarray((31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))


def _itrunc_div(a, b):
    """JAX expression for the source integer division, which truncates toward zero."""
    return jnp.trunc(jnp.asarray(a) / jnp.asarray(b)).astype(jnp.int32)


def gregorian2julian_date(day, month, year):
    """Computes the Julian Date (gregorian2julian), or the number of days since
    1 January 4713 BC, given a Gregorian Calender date (day, month, year).

    References:
    Fliegel, H. F. and van Flandern, T. C.,
    Communications of the ACM, Vol. 11, No. 10 (October, 1968)

    Returns:
        gregorian2julian_date

    Fortran declarations (gregorian2julian_date in src/CLUBB_core/calendar.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          day,        & ! Gregorian Calendar Day for given Month        [dd]
          month,      & ! Gregorian Calendar Month for given Year       [mm]
          year          ! Gregorian Calendar Year                       [yyyy]
    """
    # ---- Begin Code ----
    I = jnp.asarray(year, dtype=jnp.int32)
    J = jnp.asarray(month, dtype=jnp.int32)
    K = jnp.asarray(day, dtype=jnp.int32)
    gregorian2julian_date = (
        K - 32075
        + _itrunc_div(1461 * (I + 4800 + _itrunc_div(J - 14, 12)), 4)
        + _itrunc_div(367 * (J - 2 - _itrunc_div(J - 14, 12) * 12), 12)
        - _itrunc_div(3 * _itrunc_div(I + 4900 + _itrunc_div(J - 14, 12), 100), 4)
    )
    return gregorian2julian_date


def julian2gregorian_date(julian_date):
    """Computes the Gregorina Calendar date (day, month, year)
    given the Julian date (julian_date).

    References:
    Fliegel, H. F. and van Flandern, T. C.,
    Communications of the ACM, Vol. 11, No. 10 (October, 1968)
    http://portal.acm.org/citation.cfm?id=364097

    Returns:
        (K, J, I)

    Fortran declarations (julian2gregorian_date in src/CLUBB_core/calendar.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: julian_date ! Julian date being converted from

        !--------------------------- Output Variables ---------------------------
        integer, intent(out)::  &
          day,     & ! Gregorian calender day for given Month       [dd]
          month,   & ! Gregorian calender month for given Year      [mm]
          year       ! Gregorian calender year                      [yyyy]
    """
    # ---- Begin Code ----
    L = jnp.asarray(julian_date, dtype=jnp.int32) + 68569
    N = _itrunc_div(4 * L, 146097)
    L = L - _itrunc_div(146097 * N + 3, 4)
    I = _itrunc_div(4000 * (L + 1), 1461001)
    L = L - _itrunc_div(1461 * I, 4) + 31
    J = _itrunc_div(80 * L, 2447)
    K = L - _itrunc_div(2447 * J, 80)
    L = _itrunc_div(J, 11)
    J = J + 2 - 12 * L
    I = 100 * (N - 49) + I + L
    return K, J, I


def leap_year(year):
    """Determines if the given year is a leap year.

    Returns:
        leap_year

    Fortran declarations (leap_year in src/CLUBB_core/calendar.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: year ! Gregorian Calendar Year [yyyy]
    """
    # ---- Begin Code ----
    year = jnp.asarray(year, dtype=jnp.int32)
    leap_year = (year % 4 == 0) & ~((year % 100 == 0) & (year % 400 != 0))
    return leap_year


def compute_current_date_api(
    previous_day, previous_month, previous_year,
    seconds_since_previous_date,
):
    """Computes the current Gregorian date from a previous date and
    the seconds that have transpired since that date.

    Returns:
        (current_day, current_month, current_year, seconds_since_current_date)

    Fortran declarations (compute_current_date_api in src/CLUBB_core/calendar.F90):

        !--------------------------- Input Variables ---------------------------

        ! Previous date
        integer, intent(in) :: &
          previous_day,    & ! Day of the month      [dd]
          previous_month,  & ! Month of the year     [mm]
          previous_year      ! Year                  [yyyy]

        real(kind=time_precision), intent(in) :: &
          seconds_since_previous_date ! [s]

        !--------------------------- Output Variables ---------------------------

        ! Current date
        integer, intent(out) :: &
          current_day,     & ! Day of the month      [dd]
          current_month,   & ! Month of the year     [mm]
          current_year       ! Year                  [yyyy]

        real(kind=time_precision), intent(out) :: &
          seconds_since_current_date
    """
    # ---- Begin Code ----

    # Using Julian dates, add the number of days that the model has run.
    days_since_1jan4713bc = gregorian2julian_date(
        previous_day, previous_month, previous_year,
    )
    days_since_start = jnp.floor(
        seconds_since_previous_date / sec_per_day,
    ).astype(jnp.int32)
    days_since_1jan4713bc = days_since_1jan4713bc + days_since_start
    seconds_since_current_date = (
        seconds_since_previous_date - days_since_start * sec_per_day
    )
    current_day, current_month, current_year = julian2gregorian_date(days_since_1jan4713bc)
    return current_day, current_month, current_year, seconds_since_current_date


def gregorian2julian_day(day, month, year):
    """This subroutine determines the Julian day (1-366)
    for a given Gregorian calendar date(e.g. July 1, 2008).

    Returns:
        julian_day

    Fortran declarations (gregorian2julian_day in src/CLUBB_core/calendar.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
         day,             & ! Day of the Month      [dd]
         month,           & ! Month of the Year     [mm]
         year               ! Year                  [yyyy]
    """
    # ---- Begin Code ----
    month = jnp.asarray(month, dtype=jnp.int32)
    # Add the days from the preceding months.
    julian_day = jnp.asarray(day, dtype=jnp.int32) + jnp.sum(
        jnp.where(jnp.arange(1, 13) < month, days_per_month, 0)
    )
    # Kluge for leap years: March 1 and later require the extra day.
    julian_day = julian_day + jnp.where(leap_year(year) & (month > 2), 1, 0)
    # TODO(port-mirror): source uses ``error stop`` for an invalid Julian day.
    # The radiation configuration boundary validates external dates before JIT;
    # dates produced by ``compute_current_date_api`` are valid by construction.
    return julian_day
