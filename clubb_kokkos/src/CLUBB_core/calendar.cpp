// Calendar conversions and elapsed-time rollover for CLUBB.
//
// gregorian2julian_date/julian2gregorian_date convert complete calendar dates;
// leap_year and gregorian2julian_day provide Gregorian leap/day-of-year rules.
// compute_current_date_api returns the date and remaining seconds after elapsed
// model time. Public tables and calling contracts are declared in calendar.hpp.
// Reference: Fliegel and van Flandern, Communications of the ACM 11(10), 1968.
#include "calendar.hpp"
#include "constants_clubb.hpp"

#include <cmath>
#include <stdexcept>

namespace clubb_kokkos {

// -----------------------------------------------------------------------------
int gregorian2julian_date(int day, int month, int year) { // In
  // Description:
  //   Compute the Julian date, or days since 1 January 4713 BC, from a
  //   Gregorian day, month and year.
  // Reference:
  //   Fliegel, H. F. and van Flandern, T. C., Communications of the ACM,
  //   Vol. 11, No. 10 (October, 1968).
  // day: Gregorian calendar day for the given month [dd].
  // month: Gregorian calendar month for the given year [mm].
  // year: Gregorian calendar year [yyyy].
  int I; // Gregorian year used by the integer conversion [yyyy].
  int J; // Gregorian month used by the integer conversion [mm].
  int K; // Gregorian day used by the integer conversion [dd].

  I = year;
  J = month;
  K = day;

  // C++ integer division, like Native INTEGER division, truncates toward zero.
  return K - 32075 + 1461 * (I + 4800 + (J - 14) / 12) / 4 +
         367 * (J - 2 - (J - 14) / 12 * 12) / 12 - 3 * ((I + 4900 + (J - 14) / 12) / 100) / 4;
}

// -----------------------------------------------------------------------------
void julian2gregorian_date(int julian_date,                // In
                           int& day, int& month, int& year // Out
) {
  // Description:
  //   Compute the Gregorian calendar day, month and year from a Julian date.
  // Reference:
  //   Fliegel, H. F. and van Flandern, T. C., Communications of the ACM,
  //   Vol. 11, No. 10 (October, 1968).
  //   http://portal.acm.org/citation.cfm?id=364097
  // julian_date: Integer Julian date being converted [-].
  // day: Gregorian calendar day for the given month [dd].
  // month: Gregorian calendar month for the given year [mm].
  // year: Gregorian calendar year [yyyy].
  int I, J, K, N, L; // Integer working values in the date-conversion algorithm [-].

  L = julian_date + 68569;      // Known magic number.
  N = 4 * L / 146097;           // Known magic number.
  L = L - (146097 * N + 3) / 4; // Known magic number.
  I = 4000 * (L + 1) / 1461001; // Known magic number.
  L = L - 1461 * I / 4 + 31;    // Known magic number.
  J = 80 * L / 2447;            // Known magic number.
  K = L - 2447 * J / 80;        // Known magic number.
  L = J / 11;                   // Known magic number.
  J = J + 2 - 12 * L;           // Known magic number.
  I = 100 * (N - 49) + I + L;   // Known magic number.

  year  = I;
  month = J;
  day   = K;
}

// -----------------------------------------------------------------------------
bool leap_year(int year) { // In
  // Description:
  //   Determine whether the given Gregorian calendar year is a leap year.
  // References: None.
  // year: Gregorian calendar year [yyyy].
  return year % 4 == 0 && !(year % 100 == 0 && year % 400 != 0);
}

// -----------------------------------------------------------------------------
void compute_current_date_api(int previous_day, int previous_month, // In
                              int previous_year,                    // In
                              double seconds_since_previous_date,   // In
                              int& current_day, int& current_month, // Out
                              int& current_year,                    // Out
                              double& seconds_since_current_date    // Out
) {
  // Description:
  //   Compute the current Gregorian date from a previous date and elapsed
  //   seconds, returning seconds since the current day's start.
  // References: None.
  // previous_day: Previous day of the month [dd].
  // previous_month: Previous month of the year [mm].
  // previous_year: Previous year [yyyy].
  // seconds_since_previous_date: Elapsed time from the previous date [s].
  // current_day: Current day of the month [dd].
  // current_month: Current month of the year [mm].
  // current_year: Current year [yyyy].
  // seconds_since_current_date: Remaining seconds since the current date [s].
  int days_since_1jan4713bc; // Julian date of the starting/current date [days].
  int days_since_start;      // Whole elapsed days, rounded toward negative infinity [days].

  // Using Julian dates, add the days that the model has been running.
  // Determine the Julian date of the starting Gregorian date.
  days_since_1jan4713bc = gregorian2julian_date(previous_day,                 // In
                                                previous_month, previous_year // In
  );

  // Determine the number of elapsed days since the start date.
  days_since_start =
      static_cast<int>(std::floor(seconds_since_previous_date / static_cast<double>(sec_per_day)));

  // Set days_since_1jan4713bc to the present Julian date.
  days_since_1jan4713bc = days_since_1jan4713bc + days_since_start;

  // Set present time to seconds since the Julian date.
  seconds_since_current_date =
      seconds_since_previous_date -
      (static_cast<double>(days_since_start) * static_cast<double>(sec_per_day));

  julian2gregorian_date(days_since_1jan4713bc,                   // In
                        current_day, current_month, current_year // Out
  );
}

// -----------------------------------------------------------------------------
int gregorian2julian_day(int day, int month, int year) { // In
  // Description:
  //   Determine the Julian day of year (1–366) for a Gregorian calendar date.
  // References: None.
  // day: Gregorian day of the month [dd].
  // month: Gregorian month of the year [mm].
  // year: Gregorian calendar year [yyyy].
  int previous_month_days = 0; // Native SUM of the preceding month lengths [days].

  // Add the days from the previous months in Native storage order.
  for (int m = 0; m < month - 1; ++m)
    previous_month_days += days_per_month[m];

  int gregorian2julian_day = day + previous_month_days; // Native function result [days].

  // Leap-year correction: 29 February needs no bump; 1 March and later do.
  if (leap_year(year) && month > 2) {
    gregorian2julian_day = gregorian2julian_day + 1;
  }

  if ((leap_year(year) && gregorian2julian_day > 366) ||
      (!leap_year(year) && gregorian2julian_day > 365)) {
    throw std::runtime_error("Problem with Julian day conversion in gregorian2julian_day.");
  }

  return gregorian2julian_day;
}

} // namespace clubb_kokkos
