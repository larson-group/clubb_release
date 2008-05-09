! $Id: cos_solar_zen.F,v 1.4 2008-05-09 22:25:24 faschinj Exp $
!-----------------------------------------------------------------------
        module cos_solar_zen_mod

        private

        public :: cos_solar_zen

        private :: gregorian2julian, julian2gregorian
        
        contains
        double precision function cos_solar_zen
     .           ( day, month, year, current_time, lat_in_degrees,
     .             lon_in_degrees )
!       Description:
!       A function based on coefficients from Liou and the Clayson and 
!       Curry formula.  Approximates the cosine of the solar zenith
!       angle anywhere in the world based on current Greenwich mean 
!       time and the latitude and longitude.

!       References:
!       Clayson and Curry formula from C. A. Clayson and J. A. Curry , 
!         J. Geophys.
!         Res. Vol. 101, No. C12, Pages 28515-28528, 15 Dec. 1996.
!       Liou ``An Introduction to Atmospheric Radiation'' 
!         Table 2.2 and Eqn. 2.2.10
!-----------------------------------------------------------------------
        use constants, only: pi_dp, fstderr, sec_per_day, sec_per_hr
        use stats_prec, only: time_prec

        implicit none

        ! External
        intrinsic :: sin, cos, mod, abs, int

        ! Constant Parameters
        integer, dimension(12), parameter ::
     .  ndays = (/31, 28, 31, 30, 31, 30,
     .            31, 31, 30, 31, 30, 31/)  ! Days per month (Jan..Dec)

        ! Liou's coefficients
        double precision, parameter ::
     .  c0 =  0.006918,  ! [-]
     .  c1 = -0.399912,  ! [-]
     .  c2 = -0.006758,  ! [-]
     .  c3 = -0.002697,  ! [-]
     .  d1 =  0.070257,  ! [-]
     .  d2 =  0.000907,  ! [-]
     .  d3 =  0.000148   ! [-]

        ! Input Variables
        integer, intent(in) :: 
     .  day,   ! Day of month at model start
     .  month, ! Month of year at model start
     .  year   ! Year at model start

        real(kind=time_prec), intent(in) ::
     .  current_time   ! Current time since start date [s]

        real, intent(in) ::
     .  lat_in_degrees,! Latitude       [degrees_N]
     .  lon_in_degrees ! Longitude      [degrees_E]

!        ! Output Variables
!        double precision ::
!     .  amu0 ! Cosine of the solar zenith angle

        ! Local Variables
        double precision ::
     .  t, h, 
     .  delta,
     .  zln,
     .  longang,
     .  latang,
     .  hour, 
     .  present_time

        integer ::
     .  jul_day, j, days_in_year,
     .  present_year, present_month, present_day, days_since_start,
     .  days_since_1jan4713bc  
        
        ! Using Julian dates we are able to add the days that the model
        ! has been running
        
        ! Determine the Julian Date of the starting date, 
        !    written in Gregorian (year, month, day) form
        days_since_1jan4713bc = gregorian2julian( year, month, day )
        
        ! Determine the amount of days that have passed since start date
        days_since_start = int( current_time / sec_per_day )

        ! Set days_since_1jan4713 to the present Julian date
        days_since_1jan4713bc = days_since_1jan4713bc + days_since_start

        ! Set Present time to be seconds since the Julian date
        present_time = current_time - ( days_since_start * sec_per_day )

        call julian2gregorian 
     .         ( days_since_1jan4713bc,                  
     .           present_year, present_month, present_day)
        
        ! Now the current Gregorian Date of the model should now be in
        ! the variables present_year, present_month, and present_day
        ! Joshua Fasching May 2008
        
        ! A version of Dr. Golaz's leap year code (from outputgrads)
        !   **** Turn into subroutine; put in calendar module; use here and in outputgrads ***
        if ( (mod(present_year,4) == 0) .and.
     .       (.not.(  mod(present_year,100) == 0 
     .                .and. mod(present_year,400) /= 0 
     .              ) 
     .       )
     .      )  then 
                 days_in_year = 366
        else
                 days_in_year = 365
        end if

        jul_day = present_day

        ! Add the days from the previous months
        do j = 1, present_month-1, 1
          jul_day = jul_day + ndays(j)
        end do
        
        ! Kluge for a leap year
        if ( days_in_year == 366 .and. present_month > 2 ) then
            jul_day = jul_day + 1
        end if

        !delta_in_degrees = delta*(180/pi_dp)

        ! Compute hour angle (old code)
        ! h = 2*pi_dp*t_since_noon/86400

        ! Determine the number of hours
        hour = present_time / sec_per_hr

        t = 2*pi_dp*(jul_day-1)/days_in_year

        delta = c0 
     .        + c1*cos( t ) + d1*sin( t )
     .        + c2*cos( 2*t ) + d2*sin( 2*t )
     .        + c3*cos( 3*t ) + d3*sin( 3*t )

        ! Just in case
        if ( jul_day > days_in_year ) then 
          stop "Problem with days solar zenith code."
        end if

        ! The angle  longang  is equivalent to the
        ! hour angle in the formula for cosZ .
        ! References: Source file zenith.f
        !   from http://magic.gfdi.fsu.edu/seaflux/DIURNAL/README.txt
        !   Clayson and Curry formula from C. A. Clayson and J. A. Curry , 
        !   J. Geophys.
        !   Res. Vol. 101, No. C12, Pages 28515-28528, 15 Dec. 1996 .

        !   June 6, 2006

        select case( int( hour ) )
        case( 0:11 )
          zln = 180.00 - hour*15.00
        case( 12:23 )
          zln = 540.00 - hour*15.00
        case default
          write(unit=fstderr,fmt=*) "Hour=", hour
          stop " > 24 hours in cosine solar zenith code"
        end select

        longang = abs( lon_in_degrees - zln ) * pi_dp/180.0
        latang  = lat_in_degrees * pi_dp/180.0


        ! Cosine of the solar zenith angle (sometimes denoted amu0).
        cos_solar_zen = sin(latang)*sin(delta)
     .       + cos(latang)*cos(delta)*cos(longang)

        !write(*,'(a,f15.6)') "cosine solar zenith", cos_solar_zen !%% debug

        return
        end function cos_solar_zen
        
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
        
        end module cos_solar_zen_mod
