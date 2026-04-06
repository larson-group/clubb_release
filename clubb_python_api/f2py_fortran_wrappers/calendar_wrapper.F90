! calendar_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module calendar

subroutine f2py_julian2gregorian_date(julian_date, day, month, year)

  use calendar, only: julian2gregorian_date

  implicit none

  integer, intent(in) :: julian_date
  integer, intent(out) :: day, month, year

  call julian2gregorian_date(julian_date, day, month, year)

end subroutine f2py_julian2gregorian_date

subroutine f2py_gregorian2julian_day(day, month, year, julian_day)

  use calendar, only: gregorian2julian_day

  implicit none

  integer, intent(in) :: day, month, year
  integer, intent(out) :: julian_day

  julian_day = gregorian2julian_day(day, month, year)

end subroutine f2py_gregorian2julian_day


subroutine f2py_leap_year(year, is_leap_year)

  use calendar, only: leap_year

  implicit none

  integer, intent(in) :: year
  logical, intent(out) :: is_leap_year

  is_leap_year = leap_year(year)

end subroutine f2py_leap_year

subroutine f2py_compute_current_date(previous_day, previous_month, previous_year, &
    seconds_since_previous_date, current_day, current_month, current_year, &
    seconds_since_current_date)

  use clubb_precision, only: time_precision
  use calendar, only: compute_current_date_api

  implicit none

  integer, intent(in) :: previous_day, previous_month, previous_year
  real(time_precision), intent(in) :: seconds_since_previous_date
  integer, intent(out) :: current_day, current_month, current_year
  real(time_precision), intent(out) :: seconds_since_current_date

  call compute_current_date_api(previous_day, previous_month, previous_year, &
    seconds_since_previous_date, current_day, current_month, current_year, &
    seconds_since_current_date)

end subroutine f2py_compute_current_date
