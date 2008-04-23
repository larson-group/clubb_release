! $Id: cos_solar_zen.F,v 1.2 2008-04-23 17:52:52 faschinj Exp $
!-----------------------------------------------------------------------
        module cos_solar_zen_mod

        private

        public :: cos_solar_zen
 
        contains
        function cos_solar_zen
     .           ( day, month, year, current_time, lat_in_degrees,
     .             lon_in_degrees ) result( amu0 )
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
        use constants, only: pi_dp, fstderr

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
     .  day,   ! Day of month
     .  month, ! Month of year
     .  year   ! Current year

        real, intent(in) ::
     .  current_time,  ! Current time   [s]
     .  lat_in_degrees,! Latitude       [degrees_N]
     .  lon_in_degrees ! Longitude      [degrees_E]

        ! Output Variables
        double precision ::
     .  amu0 ! Cosine of the solar zenith angle

        ! Local Variables
        double precision ::
     .  t, h,
     .  delta,
     .  zln,
     .  longang,
     .  latang,
     .  hour 

        integer ::
     .  jul_day, j, daysinyear

        ! A version of Dr. Golaz's leap year code (from outputgrads)
        if ( (mod(year,4) == 0) .and. 
     .       (.not.(  mod(year,100) == 0 .and. mod(year,400) /= 0 ) )
     .     ) then 
          daysinyear = 366
        else
          daysinyear = 365
        end if

        jul_day = day

        ! Add the days from the previous months
        do j = 1, month-1, 1
          jul_day = jul_day + ndays(j)
        end do

        ! Kluge for a leap year
        if ( daysinyear == 366 .and. month > 2 ) jul_day = jul_day + 1

        t = 2*pi_dp*(jul_day-1)/daysinyear

        delta = c0 
     .        + c1*cos( t ) + d1*sin( t )
     .        + c2*cos( 2*t ) + d2*sin( 2*t )
     .        + c3*cos( 3*t ) + d3*sin( 3*t )

        !delta_in_degrees = delta*(180/pi_dp)

        ! Compute hour angle (old code)
        ! h = 2*pi_dp*t_since_noon/86400

        hour = current_time / 3600.0

        ! Fix the time if we've been running for more than a day
        j = int( hour / 24.0 )
        hour = hour - ( j * 24.0 )
        jul_day = jul_day + j

        ! Just in case
        if ( jul_day > daysinyear ) then 
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


        ! Cosine of the solar zenith angle.
        amu0 = sin(latang)*sin(delta)
     .       + cos(latang)*cos(delta)*cos(longang)

        !write(*,'(a,f15.6)') "cosine solar zenith", amu0 !%% debug

        return
        end function cos_solar_zen

        end module cos_solar_zen_mod
