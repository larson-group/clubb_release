// Check calendar boundaries not reached by the retained Native I/O callers.
#include "CLUBB_core/calendar.hpp"
#include <iostream>
#include <stdexcept>

int main() {
  // Exercise century leap rules, both rollover directions and the Native
  // upper-day error branch without changing the model comparison settings.
  try {
    using namespace clubb_kokkos;
    if (leap_year(1900) || !leap_year(2000) || !leap_year(2024) || leap_year(2025))
      throw std::runtime_error("Gregorian century leap rule changed");
    if (gregorian2julian_date(1, 1, 2000) != 2451545 || gregorian2julian_day(1, 3, 2000) != 61 ||
        gregorian2julian_day(1, 3, 1900) != 60)
      throw std::runtime_error("Julian date/day conversion changed");

    int day, month, year;
    double seconds;
    julian2gregorian_date(2451545, day, month, year);
    if (day != 1 || month != 1 || year != 2000)
      throw std::runtime_error("Julian date inverse changed");
    compute_current_date_api(1, 1, 2001, -0.25, day, month, year, seconds);
    if (day != 31 || month != 12 || year != 2000 || seconds != 86399.75)
      throw std::runtime_error("Negative elapsed time must floor to the previous day");
    compute_current_date_api(28, 2, 2000, 86400.5, day, month, year, seconds);
    if (day != 29 || month != 2 || year != 2000 || seconds != 0.5)
      throw std::runtime_error("Leap-day rollover changed");

    bool rejected = false;
    try {
      gregorian2julian_day(367, 1, 2001);
    } catch (const std::runtime_error& error) {
      rejected = std::string(error.what()) ==
                 "Problem with Julian day conversion in gregorian2julian_day.";
    }
    if (!rejected)
      throw std::runtime_error("Native upper-day error branch changed");
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
