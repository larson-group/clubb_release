! $Id: cos_solar_zen.F90,v 1.1 2008-07-22 16:04:12 faschinj Exp $
!-----------------------------------------------------------------------
        module cos_solar_zen_mod

        implicit none

        public :: cos_solar_zen

        private ! Default Scope

        contains

        double precision function cos_solar_zen & 
                 ( day, month, year, current_time, lat_in_degrees, & 
                   lon_in_degrees )
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

        use constants, only: pi_dp, fstderr, sec_per_day, sec_per_hr ! Variable(s)

        use stats_precision, only: time_precision ! Variable(s)

        use calendar, only:  & 
            gregorian2julian_day, compute_current_date, leap_year ! Procedure(s)
        
        implicit none

        ! External
        intrinsic :: sin, cos, mod, abs, int

        ! Constant Parameters

        ! Liou's coefficients
        double precision, parameter :: & 
        c0 =  0.006918,   & ! [-]
        c1 = -0.399912,   & ! [-]
        c2 = -0.006758,   & ! [-]
        c3 = -0.002697,   & ! [-]
        d1 =  0.070257,   & ! [-]
        d2 =  0.000907,   & ! [-]
        d3 =  0.000148   ! [-]

        ! Input Variables
        integer, intent(in) ::  & 
        day,    & ! Day of month at model start
        month,  & ! Month of year at model start
        year   ! Year at model start

        real(kind=time_precision), intent(in) :: & 
        current_time   ! Current time since start date [s]

        real, intent(in) :: & 
        lat_in_degrees, & ! Latitude       [degrees_N]
        lon_in_degrees ! Longitude      [degrees_E]

!        ! Output Variables
!        double precision ::
!     .  amu0 ! Cosine of the solar zenith angle

        ! Local Variables
        double precision :: & 
        t,  & 
!     .  h,
        delta, & 
        zln, & 
        longang, & 
        latang, & 
        hour,  & 
        present_time

        integer :: & 
        jul_day, days_in_year, & 
        present_year, present_month, present_day 
        
        call compute_current_date( day, month, & 
                                   year, & 
                                   current_time, & 
                                   present_day, present_month, & 
                                   present_year, & 
                                   present_time )

        jul_day = gregorian2julian_day( present_day, present_month, & 
                                   present_year )

        if( leap_year( present_year ) )then
                days_in_year = 366
        else
                days_in_year = 365
        endif

        !delta_in_degrees = delta*(180/pi_dp)

        ! Compute hour angle (old code)
        ! h = 2*pi_dp*t_since_noon/86400

        ! Determine the number of hours
        hour = present_time / sec_per_hr

        t = 2*pi_dp*(jul_day-1)/days_in_year

        delta = c0  & 
              + c1*cos( t ) + d1*sin( t ) & 
              + c2*cos( 2*t ) + d2*sin( 2*t ) & 
              + c3*cos( 3*t ) + d3*sin( 3*t )

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
          zln = 0.0
          write(unit=fstderr,fmt=*) "Hour=", hour
          stop " > 24 hours in cosine solar zenith code"
        end select

        longang = abs( lon_in_degrees - zln ) * pi_dp/180.0
        latang  = lat_in_degrees * pi_dp/180.0


        ! Cosine of the solar zenith angle (sometimes denoted amu0).
        cos_solar_zen = sin(latang)*sin(delta) & 
             + cos(latang)*cos(delta)*cos(longang)

        !write(*,'(a,f15.6)') "cosine solar zenith", cos_solar_zen !%% debug

        return

        end function cos_solar_zen
                
        end module cos_solar_zen_mod
