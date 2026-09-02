"""JAX-compatible port of ``src/CLUBB_core/calendar.F90``.

Description:
    Gregorian and Julian calendar conversion routines used by solar-zenith
    radiation. Fortran subroutine outputs return tuples in source output order.

JAX adaptation:
    Integer calculations remain JAX arrays so calendar rollover works while
    elapsed model time is traced.
"""

import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
from clubb_jax.src.CLUBB_core.constants_clubb import sec_per_day

configure_jax_precision()

month_names = ("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
days_per_month = jnp.asarray((31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))


def _itrunc_div(a, b):
    """JAX expression for the source integer division, which truncates toward zero."""
    return jnp.trunc(jnp.asarray(a) / jnp.asarray(b)).astype(jnp.int32)


def gregorian2julian_date(day, month, year):
    """Compute the Julian Date from a Gregorian date.

    Reference: Fliegel and van Flandern (1968).
    """
    # ---- Begin Code ----
    I = jnp.asarray(year, dtype=jnp.int32)
    J = jnp.asarray(month, dtype=jnp.int32)
    K = jnp.asarray(day, dtype=jnp.int32)
    return (
        K - 32075
        + _itrunc_div(1461 * (I + 4800 + _itrunc_div(J - 14, 12)), 4)
        + _itrunc_div(367 * (J - 2 - _itrunc_div(J - 14, 12) * 12), 12)
        - _itrunc_div(3 * _itrunc_div(I + 4900 + _itrunc_div(J - 14, 12), 100), 4)
    )


def julian2gregorian_date(julian_date):
    """Compute the Gregorian calendar date from a Julian date.

    Returns Fortran outputs ``day, month, year`` in that order.
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
    """Determine whether the given year is a leap year."""
    # ---- Begin Code ----
    year = jnp.asarray(year, dtype=jnp.int32)
    return (year % 4 == 0) & ~((year % 100 == 0) & (year % 400 != 0))


def compute_current_date_api(
    previous_day, previous_month, previous_year,
    seconds_since_previous_date,
):
    """Compute the current Gregorian date from a previous date and elapsed seconds.

    Returns Fortran outputs ``current_day, current_month, current_year,
    seconds_since_current_date`` in that order.
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
    """Determine the Julian day (1-366) for a given Gregorian calendar date."""
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


__all__ = [
    "month_names", "days_per_month", "gregorian2julian_date",
    "julian2gregorian_date", "leap_year", "compute_current_date_api",
    "gregorian2julian_day",
]
