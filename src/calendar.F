        module calendar

        implicit none

        public :: gregorian2julian, julian2gregorian, leap_year

        private

        contains
!-----------------------------------------------------------------------        
        integer function gregorian2julian (year, month, day)
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
        INTEGER,intent(in) :: YEAR, MONTH, DAY
        
        ! Local Variables
        integer :: I,J,K

        I= YEAR
        J= MONTH
        K= DAY

        gregorian2julian = K-32075+1461*(I+4800+(J-14)/12)/4+367*
     .         (J-2-(J-14)/12*12)/12-3*((I+4900+(J-14)/12)/100)/4

        return
        end function gregorian2julian
        
!------------------------------------------------------------------
        subroutine julian2gregorian(julian_date, year, month, day)
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
     .  year, ! Gregorian calender year
     .  month, ! Gregorian calender month
     .  day ! Gregorian calender day
        
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
        
        end subroutine julian2gregorian
!-----------------------------------------------------------------------------
        logical function leap_year( year )
!
!       Description: Determines if the given year is a leap year.
!-----------------------------------------------------------------------------
        implicit none

        integer, intent(in) :: year

        leap_year = ( (mod(year,4) == 0) .and.
     .       (.not.(  mod(year,100) == 0 .and. mod(year,400) /= 0 ) ) )
  
        end function leap_year
        
        end module calendar
