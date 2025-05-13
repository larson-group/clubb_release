! $Id$
!-----------------------------------------------------------------------
module cos_solar_zen_module

  use clubb_precision, only: &
    dp ! double precision

  implicit none

  public :: cos_solar_zen

  private ! Default Scope

  contains

!-----------------------------------------------------------------------
  function cos_solar_zen & 
           ( day, month, year, current_time, lat_in_degrees, & 
             lon_in_degrees )

! Description:
!   A function based on coefficients from Liou and the Clayson and
!   Curry formula.  Approximates the cosine of the solar zenith
!   angle anywhere in the world based on current Greenwich mean
!   time and the latitude and longitude.

! References:
!   Clayson and Curry formula from C. A. Clayson and J. A. Curry ,
!   J. Geophys.
!   Res. Vol. 101, No. C12, Pages 28515-28528, 15 Dec. 1996.
!   Liou ``An Introduction to Atmospheric Radiation''
!     Table 2.2 and Eqn. 2.2.10
!-----------------------------------------------------------------------

    use constants_clubb, only: pi_dp, fstderr, sec_per_hr,&
                               radians_per_deg_dp ! Variable(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use calendar, only: &
        gregorian2julian_day, compute_current_date_api, leap_year ! Procedure(s)

    implicit none

    ! External
    intrinsic :: sin, cos, mod, abs, int

    ! Constant Parameters

    ! Liou's coefficients
    real( kind = dp ), parameter :: &
      c0 =  0.006918_dp,   & ! [-]
      c1 = -0.399912_dp,   & ! [-]
      c2 = -0.006758_dp,   & ! [-]
      c3 = -0.002697_dp,   & ! [-]
      d1 =  0.070257_dp,   & ! [-]
      d2 =  0.000907_dp,   & ! [-]
      d3 =  0.000148_dp      ! [-]

    ! Input Variables
    integer, intent(in) :: &
      day,    & ! Day of month at model start
      month,  & ! Month of year at model start
      year      ! Year at model start

    real(kind=time_precision), intent(in) :: &
      current_time   ! Current time since start date [s]

    real( kind = core_rknd ), intent(in) :: &
      lat_in_degrees, & ! Latitude       [degrees_N]
      lon_in_degrees    ! Longitude      [degrees_E]

    ! Return Variable
    real( kind = dp ) :: &
      cos_solar_zen

    ! Local Variables
    real( kind = dp ) :: &
      t,  &
      delta, &
      zln, &
      longang, &
      latang, &
      hour

    real( kind = time_precision ) :: &
      present_time

    integer :: &
      jul_day, days_in_year, &
      present_year, present_month, present_day

    call compute_current_date_api( day, month, &
                                   year, &
                                   current_time, &
                                   present_day, present_month, &
                                   present_year, &
                                   present_time )

    jul_day = gregorian2julian_day( present_day, present_month, &
                               present_year )

    if ( leap_year( present_year ) ) then
      days_in_year = 366
    else
      days_in_year = 365
    end if

    !delta_in_degrees = delta*(180/pi_dp)

    ! Compute hour angle (old code)
    ! h = 2*pi_dp*t_since_noon/86400

    ! Determine the number of hours
    hour = present_time / real(sec_per_hr,kind=time_precision)

    t = 2._dp*pi_dp*real(jul_day-1, kind=dp) / real(days_in_year, kind=dp)

    delta = c0  & 
          + c1*cos( t ) + d1*sin( t ) & 
          + c2*cos( 2._dp*t ) + d2*sin( 2._dp*t ) & 
          + c3*cos( 3._dp*t ) + d3*sin( 3._dp*t )

! The angle  longang  is equivalent to the
! hour angle in the formula for cosZ .
! References: Source file zenith.f
!   from http://magic.gfdi.fsu.edu/seaflux/DIURNAL/README.txt
!   Clayson and Curry formula from C. A. Clayson and J. A. Curry ,
!   J. Geophys.
!   Res. Vol. 101, No. C12, Pages 28515-28528, 15 Dec. 1996 .

!   June 6, 2006

    select case ( int( hour ) )
    case ( 0:11 )
      zln = 180.00_dp - hour*15.00_dp ! Known magic number
    case ( 12:23 )
      zln = 540.00_dp - hour*15.00_dp ! Known magic number
    case default
      zln = 0.0_dp
      write(unit=fstderr,fmt=*) "Hour=", hour
      error stop " > 24 hours in cosine solar zenith code"
    end select

    longang = abs( real(lon_in_degrees, kind=dp) - zln ) * radians_per_deg_dp
    latang  = real(lat_in_degrees, kind=dp) * radians_per_deg_dp


    ! Cosine of the solar zenith angle (sometimes denoted amu0).
    cos_solar_zen = sin( latang ) * sin( delta ) & 
                  + cos( latang ) * cos( delta ) * cos( longang )

    !write(*,'(a,f15.6)') "cosine solar zenith angle=", cos_solar_zen !%% debug

    return

  end function cos_solar_zen

end module cos_solar_zen_module
