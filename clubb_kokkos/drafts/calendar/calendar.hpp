// Convert Gregorian dates, Julian dates and elapsed model time.
//
// Public routines:
//   gregorian2julian_date: Gregorian day/month/year to an integer Julian date.
//   julian2gregorian_date: Julian date to Gregorian day/month/year.
//   leap_year: test the Gregorian leap-year rule.
//   compute_current_date_api: advance a date by elapsed seconds, returning date and time.
//   gregorian2julian_day: return the one-based day of year, retaining the Native fatal check.
//
// Months are one-based. Integer divisions retain Native truncation toward zero;
// elapsed days use FLOOR, including negative elapsed times. Time arguments use
// double, matching Native time_precision independently of the core real kind.
// The two public month tables are immutable. Callers supply valid Native integer
// dates and month indices; no additional date-validation policy is introduced.
// Reference: Fliegel and van Flandern, Communications of the ACM 11(10), 1968.
#pragma once

#include <array>
#include <string_view>

namespace clubb_kokkos {

// Three-letter month abbreviations in Gregorian month order.
inline constexpr std::array<std::string_view, 12> month_names{
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"};

// Days in each month of a non-leap year [-].
inline constexpr std::array<int, 12> days_per_month{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

int gregorian2julian_date(int day, int month, int year); // In

void julian2gregorian_date(int julian_date,                // In
                           int& day, int& month, int& year // Out
);

bool leap_year(int year); // In

void compute_current_date_api(int previous_day, int previous_month, // In
                              int previous_year,                    // In
                              double seconds_since_previous_date,   // In
                              int& current_day, int& current_month, // Out
                              int& current_year,                    // Out
                              double& seconds_since_current_date    // Out
);

int gregorian2julian_day(int day, int month, int year); // In

} // namespace clubb_kokkos
