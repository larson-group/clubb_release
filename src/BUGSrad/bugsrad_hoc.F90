!-----------------------------------------------------------------------
! $Id: bugsrad_hoc.F90,v 1.15 2007-06-19 21:04:59 dschanen Exp $

subroutine bugsrad_hoc( alt, nz, lat_in_degrees, lon_in_degrees, &
                        day, month, year, time,                  &
                        thlm, rcm, rtm, rrm, rim,                & 
                        cf, pinpa, exner, rhom, Tsfc,            &
                        radht, Frad, thlm_forcing )
! Description:
! Does the necessary operations to interface the HOC model with
! the bugsrad subprogram.

! References:
! Stevens, et al., (2001) _Journal of Atmospheric Science_, Vol 58, p.3391-3409

! Contact for information on BUGSrad (other than this routine)
!   Norm Wood <norm@atmos.colostate.edu>

! All code external to this based on the BUGSrad source from 2004/7/10
!-----------------------------------------------------------------------

  use constants
#ifdef STATS
use hoc_stats, only: zt, zm, lstats_samp, &
    iFrad_SW, iFrad_LW, iradht_SW, iradht_LW
#endif

  implicit none

! External
  double precision, external :: cos_solar_zen 

  intrinsic :: dble, real

! Constant parameters
  integer, parameter :: &
  nlen = 1, &   ! Length of the total domain
  slen = 1      ! Length of the sub domain

! Number of levels to take from U.S. Std. Atmos tables
  integer, parameter :: std_atmos_buffer = 10 

! Number of levels to interpolate from the bottom of std_atmos to the top
! of the HOC profile, hopefully enough to eliminate cooling spikes, etc. 
  integer, parameter :: lin_int_buffer = 20

! The sum of the above to two buffers
  integer, parameter :: buffer = lin_int_buffer + std_atmos_buffer

  integer, parameter :: std_atmos_dim = 25

! Parameters from U.S. Standard Atmosphere, 1976;  Starting at 1 km altitude
! Pressure in millibars
  double precision, parameter, dimension(std_atmos_dim) ::      &
  std_pinmb = (/8.986e2, 7.950e2, 7.012e2, 6.166e2, 5.405e2,    &
                4.722e2, 4.111e2, 3.565e2, 3.080e2, 2.650e2,    & 
                2.270e2, 1.940e2, 1.658e2, 1.417e2, 1.211e2,    &
                1.035e2, 8.850e1, 7.565e1, 6.467e1, 5.529e1,    &
                4.729e1, 4.047e1, 3.467e1, 2.972e1, 2.546e1/)

! Temperature in degrees Kelvin
  double precision, parameter, dimension(std_atmos_dim) :: &
  std_tempk = (/281.6, 275.1, 268.7, 262.2, 255.7,         &
                249.2, 242.7, 236.2, 229.7, 223.2,         &
                216.8, 216.6, 216.6, 216.6, 216.6,         &
                216.6, 216.6, 216.6, 216.6, 216.6,         &
                217.6, 218.6, 219.6, 220.6, 221.6/)

! Specific Humidity ( Water Vapor / Density )
  double precision, parameter, dimension(std_atmos_dim) ::                &
  std_sp_hmdty = (/0.378e-02, 0.288e-02, 0.198e-02, 0.134e-02, 0.869e-03, & 
                   0.576e-03, 0.356e-03, 0.228e-03, 0.985e-04, 0.435e-04, &
                   0.225e-04, 0.119e-04, 0.675e-05, 0.369e-05, 0.370e-05, &
                   0.366e-05, 0.365e-05, 0.362e-05, 0.423e-05, 0.495e-05, &
                   0.634e-05, 0.806e-05, 0.104e-04, 0.130e-04, 0.165e-04/)

! Ozone ( O_3 / Density )
  double precision, parameter, dimension(std_atmos_dim) ::             &
  std_o3l   = (/0.486e-07, 0.536e-07, 0.550e-07, 0.561e-07, 0.611e-07, &
                0.682e-07, 0.814e-07, 0.989e-07, 0.152e-06, 0.218e-06, &
                0.356e-06, 0.513e-06, 0.638e-06, 0.834e-06, 0.108e-05, &
                0.138e-05, 0.197e-05, 0.263e-05, 0.337e-05, 0.427e-05, &
                0.502e-05, 0.605e-05, 0.691e-05, 0.767e-05, 0.848e-05/)

! Input Variables
  real, intent(in) :: &
  alt,           &! Maximum altitude in the model domain [m]
  lat_in_degrees,&! Latitude                             [Degrees North]
  lon_in_degrees,&! Longitude                            [Degrees East]
  time            ! Model time                           [s]
  
  integer, intent(in) :: &
  nz,              & ! Vertical extent;  i.e. nnzp in the grid class
  day, month, year   ! Time of year

  real, intent(in), dimension(nz) :: &
  thlm,  & ! Liquid potential temp.     [K]
  rcm,   & ! Liquid water mixing ratio  [kg/kg]
  rrm,   & ! Rain water mixing ratio    [kg/kg]
  rim,   & ! Ice water mixing ratio     [kg/kg]
  rtm,   & ! Total water mixing ratio   [kg/kg]
  rhom,  & ! Density                    [kg/m^3]
  cf,    & ! Cloud fraction             [%]
  pinpa, & ! Pressure                   [Pa]
  exner    ! Exner function             [-]

  real, intent(in) :: Tsfc ! Theta at the surface [K]

! Input/Output Variables
  real, intent(inout), dimension(nz) :: &
  thlm_forcing ! Theta_l LS tendency [K/s]

! Output Variables
  real, intent(out), dimension(nz) :: &
  Frad, & ! Total radiative flux       [W/m^2]
  radht   ! Total heating rate         [K/s]

! Local Variables
  real, dimension(nz) :: &
  Frad_SW, & ! SW radiative flux          [W/m^2]
  Frad_LW, & ! LW radiative flux          [W/m^2]
  radht_SW,& ! SW heating rate            [K/s]
  radht_LW   ! LW heating rate            [K/s]

! Altered 3 Oct 2005 to be buffer levels higher
  double precision, dimension(nlen,(nz-1)+buffer) :: &
  tempk,& ! Temperature            [K]
  rcil, & ! Ice mixing ratio       [kg/kg]
  o3l     ! Ozone mixing ratio     [kg/kg]

  double precision, dimension(nlen,(nz-1)+buffer+1) :: &
  Frad_uLW, & ! LW upwelling flux         [W/m^2]
  Frad_dLW, & ! LW downwelling flux       [W/m^2]
  Frad_uSW, & ! SW upwelling flux         [W/m^2]
  Frad_dSW    ! SW downwelling flux       [W/m^2]

  double precision, dimension(nlen,(nz-1)+buffer) :: &
  sp_humidity, & ! Specific humidity      [kg/kg]
  pinmb          ! Pressure in millibars  [hPa]

! Pressure in millibars for layers (calculated as an average of pinmb)
  double precision, dimension(nlen,(nz-1)+buffer+1) :: &
  playerinmb ! [hPa]

  double precision, dimension(nlen,(nz-1)+buffer) :: &
  dpl, &          ! Difference in pressure levels       [hPa]
  rrm2, rcm2, cf2 ! Two-dimensional copies of the input parameters

  double precision, dimension(nlen,(nz-1)+buffer+1) :: &
  radht_SW2,&! SW Radiative heating rate        [W/m^2]
  radht_LW2  ! LW Radiative heating rate        [W/m^2]

  double precision, dimension(nlen) :: &
  alvdr,&! Visible direct surface albedo        [-]
  alvdf,&! Visible diffuse surface albedo       [-]
  alndr,&! Near-IR direct surface albedo        [-]
  alndf  ! Near-IR diffuse surface albedo       [-]

  double precision, dimension(nlen) :: &
  slr, & ! Fraction of daylight  
  ts,  & ! Surface temperature [K]
  amu0   ! Cosine of the solar zenith angle

  double precision :: z1_fact, z2_fact ! Temp storage

  integer :: i, j, z, z1, z2  ! Loop indices

!-----------------------------------------------------------------------

! amu0 = 0.4329 ! Nov 11 Altocu value
! Calculated value
  amu0 = cos_solar_zen( day, month, year, time, lat_in_degrees, lon_in_degrees )

! Convert to millibars
  pinmb(1,1:(nz-1))    = dble( pinpa(2:nz) / 100.0 ) ! t grid in HOC

  playerinmb(1,2:nz-1) = dble( (pinmb(1,1:(nz-2)) + pinmb(1,2:nz-1))/ 2 ) ! m grid in HOC
  playerinmb(1,1)      = ( dble( pinpa(1) / 100.0 ) + pinmb(1,2) ) / 2
  playerinmb(1,nz)     = 2 * playerinmb(1,nz-1) - playerinmb(1,nz-2)


! Convert theta_l to temperature
!   kappa: Dry air gas constant / Dry air specific heat at p
!   Lv:    Latent heat of vaporization
  tempk(1,1:(nz-1)) = thlm(2:nz) * ( 1000.0d0 / pinmb(1,1:(nz-1)) )**(-kappa) &
                    + Lv*rcm(2:nz) / Cp

! Derive Specific humidity from rc & rt.
  sp_humidity(1,1:(nz-1)) = dble( rtm(2:nz) - rcm(2:nz) ) / dble( 1+rtm(2:nz) )

! Setup miscellaneous variables

! Albedo values
  alvdr = 0.1d0
  alvdf = 0.1d0
  alndr = 0.1d0
  alndf = 0.1d0

  slr  = 1.0d0 ! Fraction of daylight

! Changed for MPACE inter comparison -dschanen 3 Nov 2006
! rcil(1,1:(nz-1)+buffer) = 0.0d0 ! Assume no ice for HOC
  

! Ozone = 5.4e-5 g/m^3 from U.S. Standard Atmosphere, 1976. 
!   Convert from g to kg.
  o3l(1,1:(nz-1)) = dble( ( 5.4e-5 / rhom(1:(nz-1)) ) * 0.001 )

! Convert and transpose as needed
  rcil(1,buffer+1:(nz-1)+buffer)  = flip( dble( rim(2:nz) ), nz-1 )
  rrm2(1,buffer+1:(nz-1)+buffer)  = flip( dble( rrm(2:nz) ), nz-1 )
  rcm2(1,buffer+1:(nz-1)+buffer)  = flip( dble( rcm(2:nz) ), nz-1 )
  cf2(1,buffer+1:(nz-1)+buffer)   = flip( dble( cf(2:nz) ), nz-1 ) 

  tempk(1,buffer+1:(nz-1)+buffer) = flip( tempk(1,1:(nz-1)), nz-1 )

  sp_humidity(1,buffer+1:(nz-1)+buffer) = flip( sp_humidity(1,1:(nz-1)), nz-1 )

  pinmb(1,(buffer+1):(nz-1+buffer))        = flip( pinmb(1,1:(nz-1)), nz-1 )
  playerinmb(1,(buffer+1):(nz-1+buffer+1)) = flip( playerinmb(1,1:nz), nz )

  o3l(1,buffer+1:(nz-1)+buffer) = flip( o3l(1,1:(nz-1)), nz-1 )

! Assume these are all zero above the HOC profile
  rrm2(1,1:buffer) = 0.0d0
  rcil(1,1:buffer) = 0.0d0
  rcm2(1,1:buffer) = 0.0d0
  cf2(1,1:buffer)  = 0.0d0

! Determine the altitude of the HOC model compared to the standard atmospheric 
! profile.  Then take the values at altitude j km on the table and use 
! std_atmos_buffer kilometers of it as the top of the profile fed into BUGsrad.

  j = 1 ! initial altitude
  do while ( real( j * 1000 ) < alt )
    j = j + 1
    if ( (j + std_atmos_buffer ) > std_atmos_dim ) then
      stop "bugsrad_hoc: cannot handle this altitude" ! exceeds a 25 km altitude
    endif
  end do

  ! Add the standard atmospheric profile above the linear interpolation
  do i = 1, std_atmos_buffer, 1
    tempk(1,i)       = std_tempk((std_atmos_buffer+j)-(i-1)) 
    sp_humidity(1,i) = std_sp_hmdty((std_atmos_buffer+j)-(i-1))
    o3l(1,i)         = std_o3l((std_atmos_buffer+j)-(i-1)) 
    pinmb(1,i)       = std_pinmb((std_atmos_buffer+j)-(i-1)) 
  end do

! Do a linear interpolation to produce the levels between the standard 
! atmospheric levels and the HOC levels;  
! These levels should number the lin_int_buffer parameter
  z1 = buffer + 1
  z2 = std_atmos_buffer
  do z = buffer, std_atmos_buffer+1, -1
    z1_fact = dble(z2 - z) / dble(z2 - z1) 
    z2_fact = dble(z - z1) / dble(z2 - z1) 

    tempk(1,z) = z1_fact * tempk(1,z1) + z2_fact * tempk(1,z2)

    sp_humidity(1,z) = z1_fact * sp_humidity(1,z1) + z2_fact * sp_humidity(1,z2)

    o3l(1,z) = z1_fact * o3l(1,z1) + z2_fact * o3l(1,z2)

    pinmb(1,z) = z1_fact * pinmb(1,z1) + z2_fact * pinmb(1,z2)
  end do

  playerinmb(1,2:buffer+1) = ( pinmb(1,1:buffer) + pinmb(1,2:buffer+1) ) / 2
  playerinmb(1,1) = 2 * playerinmb(1,2) - playerinmb(1,3)

! Calculate the difference in pressure layers (including buffer levels)
  do i = 1, (nz-1)+buffer
    dpl(1,i) = playerinmb(1,i+1) - playerinmb(1,i)
  end do

  ts(1) = tempk(1,(nz-1)+buffer)

!  Write a profile for Kurt's driver program for debugging purposes
! open(10, file="profile.dat")
! write(10,'(2i4,a10)') nlen, (nz-1)+buffer, "TROPICAL"
! do i=1, (nz-1)+buffer
!   write(10,'(i4,9f12.6)') i, pinmb(1,i), playerinmb(1,i),tempk(1,i),         &
!   sp_humidity(1,i), 100000.0*o3l(1,i), rcm2(1,i), rcil(1,i),cf2(1,i),dpl(1,i)
! enddo
! write(10,'(a4,a12,3f12.6)') "","", playerinmb(1,nz+buffer), ts(1), amu0(1)
! close(10)
! pause

  call bugs_rad( nlen, slen, (nz-1)+buffer, playerinmb,          &
                 pinmb, dpl, tempk, sp_humidity,                 &
                 rcm2, rcil, rrm2, o3l,                          &
                 ts, amu0, slr, alvdf,                           &
                 alndf, alvdr, alndr, dble( sol_const ),         &
                 dble( grav ), dble( Cp ), radht_SW2, radht_LW2, &
                 Frad_dSW, Frad_uSW, Frad_dLW, Frad_uLW,         &
                 cf2 )

! Michael pointed out that this was a temperature tendency, not a theta_l
! tendency.  The 2nd line should fix both.  -dschanen 28 July 2006
  radht_SW(2:nz) = real( flip( radht_SW2(1,buffer+1:nz+buffer), nz-1 ) ) &
                   * ( 1.0 / exner(2:nz) )

  radht_LW(2:nz) = real( flip( radht_LW2(1,buffer+1:nz+buffer), nz-1 ) ) &
                   * ( 1.0 / exner(2:nz) )

  ! No radiative heating below ground
  radht_SW(1) = 0.0
  radht_LW(1) = 0.0

  radht = radht_SW + radht_LW 

! These are on the m grid, and require no adjusting
  Frad_SW(1:nz) = real( flip(  Frad_uSW(1,buffer+1:nz+buffer) &
                             - Frad_dSW(1,buffer+1:nz+buffer), nz ) )
  Frad_LW(1:nz) = real( flip(  Frad_uLW(1,buffer+1:nz+buffer) &
                             - Frad_dLW(1,buffer+1:nz+buffer), nz ) )

  Frad(1:nz) = Frad_SW(1:nz) + Frad_LW(1:nz)

  thlm_forcing(1:nz) = thlm_forcing(1:nz) + radht(1:nz)

#ifdef STATS
        if ( lstats_samp ) then

          if ( iradht_LW > 0 ) then
            zt%x(:,iradht_LW) = zt%x(:,iradht_LW) + radht_LW
            zt%n(:,iradht_LW) = zt%n(:,iradht_LW) + 1
          end if
          if ( iradht_SW > 0 ) then
            zt%x(:,iradht_SW) = zt%x(:,iradht_SW) + radht_SW
            zt%n(:,iradht_SW) = zt%n(:,iradht_SW) + 1
          end if

          if ( iFrad_SW > 0 ) then
            zm%x(:,iFrad_SW) = zm%x(:,iFrad_SW) + Frad_SW
            zm%n(:,iFrad_SW) = zm%n(:,iFrad_SW) + 1
          end if
          if ( iFrad_LW > 0 ) then
            zm%x(:,iFrad_LW) = zm%x(:,iFrad_LW) + Frad_LW
            zm%n(:,iFrad_LW) = zm%n(:,iFrad_LW) + 1
          end if

        end if
#endif /*STATS*/

  return
!-----------------------------------------------------------------------
  contains

!-----------------------------------------------------------------------
! Flips a single dimension array (i.e. a vector), so the first element 
! becomes the last and vice versa for the whole column.  This is a 
! necessary part of the code because BUGSrad and HOC store altitudes in
! reverse order
!-----------------------------------------------------------------------
  function flip( x, xdim )
    implicit none

!   Input
    integer, intent(in) :: xdim

    double precision, dimension(xdim), intent(in) :: x

!   Output
    double precision, dimension(xdim) :: flip

!   Internal
    double precision, dimension(xdim) :: tmp
    integer indx

    do indx = 1, xdim, 1
      tmp(indx) = x((xdim+1) - (indx))
    end do

    flip = tmp

    return
  end function flip
!-----------------------------------------------------------------------

end subroutine bugsrad_hoc
