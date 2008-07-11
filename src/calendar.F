!$Id: calendar.F,v 1.5 2008-07-11 15:59:38 faschinj Exp $
      module calendar

        implicit none

        public :: gregorian2julian_date, julian2gregorian_date, 
     .            leap_year, compute_current_date,
     .            gregorian2julian_day

        private ! Default Scope

        ! Constant Parameters

        ! Month Abbreviations
        character(len=3), dimension(12), public, parameter ::
     .  month = (/'JAN','FEB','MAR','APR','MAY','JUN',
     .            'JUL','AUG','SEP','OCT','NOV','DEC'/)

        ! Number of days per month (Jan..Dec) for a non leap year
        integer, public, dimension(12), parameter ::
     .  days_per_month = (/31, 28, 31, 30, 31, 30,
     .                     31, 31, 30, 31, 30, 31/)
  
        contains
!-----------------------------------------------------------------------        
        integer function gregorian2julian_date(year, month, day)
!
!       Description:
!       Computes the Julian Date (gregorian2julian), or the number of days since 
!       1 January 4713 BC, given a Gregorian Calender date (year, month, day).
!
!       Reference: 
!       Fliegel, H. F. and van Flandern, T. C., 
!       Communications of the ACM, Vol. 11, No. 10 (October, 1968)
!----------------------------------------------------------------------
        implicit none
        
        ! Output Variables
        INTEGER,intent(in) :: 
     .  YEAR,   ! Gregorian Calendar Year
     .  MONTH,  ! Gregorian Calendar Month
     .  DAY     ! Gregorian Calendar Day
        
        ! Local Variables
        integer :: I,J,K

        I= YEAR
        J= MONTH
        K= DAY

        gregorian2julian_date = K-32075+1461*(I+4800+(J-14)/12)/4+367*
     .         (J-2-(J-14)/12*12)/12-3*((I+4900+(J-14)/12)/100)/4

        return
        end function gregorian2julian_date
        
!------------------------------------------------------------------
        subroutine julian2gregorian_date(julian_date, year, month, day)
!
!       Description:
!       Computes the Gregorina Calendar date (year, month, day) 
!       given the Julian date (jd).
!
!       Reference: 
!       Fliegel, H. F. and van Flandern, T. C., 
!       Communications of the ACM, Vol. 11, No. 10 (October, 1968)
!------------------------------------------------------------------
        implicit none

        ! Input Variable
        integer,intent(in) :: julian_date ! Julian date being converted from

        ! Output Variable
        integer,intent(out):: 
     .  year,   ! Gregorian calender year       [####]
     .  month,  ! Gregorian calender month      [##]
     .  day     ! Gregorian calender day        [##]
        
        ! Local Variables
        integer :: i, j, k, n, l

        L= julian_date+68569
        N= 4*L/146097
        L= L-(146097*N+3)/4
        I= 4000*(L+1)/1461001
        L= L-1461*I/4+31
        J= 80*L/2447
        K= L-2447*J/80
        L= J/11
        J= J+2-12*L
        I= 100*(N-49)+I+L

        YEAR= I
        MONTH= J
        DAY= K

        RETURN
        
        end subroutine julian2gregorian_date
!-----------------------------------------------------------------------------
        logical function leap_year( year )
!
!       Description: Determines if the given year is a leap year.
!-----------------------------------------------------------------------------
        implicit none

        integer, intent(in) :: year ! Gregorian Calendar Year [####]

        leap_year = ( (mod(year,4) == 0) .and.
     .       (.not.(  mod(year,100) == 0 .and. mod(year,400) /= 0 ) ) )
  
        end function leap_year

!----------------------------------------------------------------------------
        subroutine compute_current_date( previous_day, previous_month,
     .                           previous_year, 
     .                           seconds_since_previous_date,
     .                           current_day, current_month,
     .                           current_year,
     .                           seconds_since_current_date )
!
!       Description: Compuntes the current Gregorian date from a previous date and
!       the seconds that have transpired since that date.
!
!----------------------------------------------------------------------------
        use stats_precision, only:
     .  time_precision  ! Variable(s)

        use constants, only:
     .  sec_per_day     ! Variable(s)

        implicit none

        integer, intent(in) ::
     .  previous_day,   ! Day of the month      [##]
     .  previous_month, ! Month of the year     [##]
     .  previous_year   ! Year                  [####]

        real(kind=time_precision), intent(in) ::
     .  seconds_since_previous_date ! [s]

        integer, intent(out) ::
     .  current_day, 
     .  current_month,
     .  current_year

        real(kind=time_precision), intent(out) ::
     .  seconds_since_current_date

        integer ::
     .  days_since_1jan4713bc,
     .  days_since_start

        ! Using Julian dates we are able to add the days that the model
        ! has been running

        ! Determine the Julian Date of the starting date, 
        !    written in Gregorian (year, month, day) form
        days_since_1jan4713bc = gregorian2julian_date( previous_year, 
     .                                    previous_month, previous_day )

        ! Determine the amount of days that have passed since start date
        days_since_start = 
     .        floor( seconds_since_previous_date / sec_per_day )

        ! Set days_since_1jan4713 to the present Julian date
        days_since_1jan4713bc = days_since_1jan4713bc + days_since_start

        ! Set Present time to be seconds since the Julian date
        seconds_since_current_date = 
     .  seconds_since_previous_date - ( days_since_start * sec_per_day )

        call julian2gregorian_date
     .         ( days_since_1jan4713bc,
     .           current_year, current_month, current_day)

      
        end subroutine compute_current_date       

!-------------------------------------------------------------------------------------
        subroutine gregorian2julian_day( day, month, year, julian_day )
!
!       Description: This subroutine determines the Julian day (1-366)
!       for a given Gregorian calendar date.
!
!-------------------------------------------------------------------------------------
        
        implicit none

        ! Input Variable(s)

        integer, intent(in) ::
     .   day,            ! Day of the Month      [##]
     .   month,          ! Month of the Year     [##]
     .   year            ! Year                  [####]

        ! Output Variable(s)

        integer, intent(out) ::
     .  julian_day      ! Julian Day (1-366)


        ! Local Variable
        integer :: j

        julian_day = day

        ! Add the days from the previous months
        do j = 1, month-1, 1
          julian_day = julian_day + days_per_month(j)
        end do

        ! Kluge for a leap year
        ! If the date were 29 Feb 2000 this would not increment julian_day
        ! However 01 March 2000 would need the 1 day bump
        if ( leap_year(year) .and. month > 2 ) then
            julian_day = julian_day + 1
        end if

        if ( ( leap_year( year ) .and. julian_day > 366 ) .or.
     .       ( .not. leap_year( year ) .and. julian_day > 365 ) ) then
                stop "Problem with Julian day conversion."
        endif

        end subroutine gregorian2julian_day

      end module calendar
