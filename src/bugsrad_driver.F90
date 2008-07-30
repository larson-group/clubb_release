!-----------------------------------------------------------------------
! $Id: bugsrad_driver.F90,v 1.3 2008-07-30 21:03:08 faschinj Exp $
module bugsrad_hoc_mod

implicit none

public :: bugsrad_hoc

private ! Default Scope

contains

subroutine bugsrad_hoc &
           ( alt, nz, lin_int_buffer,        &
             lat_in_degrees, lon_in_degrees, &
             day, month, year, time,         &
             thlm, rcm, rtm, rsnwm, rim,     & 
             cf, p_in_Pa, p_in_Pam, exner, rhom, &
             radht, Frad, thlm_forcing )
! Description:
! Does the necessary operations to interface the HOC model with
! the bugsrad subprogram.

! Grid Layout:
!
! /////////////// Layers from U.S. Standard Atmosphere   ///////////////
! ///////////////       Dimension: std_atmos_buffer      ///////////////
!                                  .
!                                  .
! ---------------         Interpolated Layers            --------------- 
! ---------------         Dimension: lin_int_buffer      --------------- 
!                                  .
!                                  .
! ///////////////          Top of HOC Grid               ///////////////
! ///////////////          Dimension: nz                 ///////////////

! References:
! Stevens, et al., (2001) _Journal of Atmospheric Science_, Vol 58, p.3391-3409
! McClatchey, et al., (1972) _Environmental Research Papers_, No. 411, p.94

! Contact for information on BUGSrad (other than this routine)
!   Norm Wood <norm@atmos.colostate.edu>

! All code external to this based on the BUGSrad source from 2004/7/10
!-----------------------------------------------------------------------

  use constants, only: fstderr, sol_const, grav, Cp ! Variable(s)
  
  use std_atmosphere_mod, only: std_atmos_dim, std_alt, std_pinmb, & ! Variable(s)
      std_T_in_K, std_sp_hmdty, std_o3l
      
  use stats_precision, only: time_precision ! Variable(s)

  use cos_solar_zen_mod, only: cos_solar_zen ! Procedure(s)

  use T_in_K_mod, only: thlm2T_in_K ! Procedure(s)

  use error_code, only: clubb_at_debug_level ! Procedure(s)

 
  use stats_type, only: stat_update_var ! Procedure(s)

  use stats_variables, only: zt, zm, l_stats_samp, & ! Variable(s)
    iFrad_SW, iFrad_LW, iradht_SW, iradht_LW

  implicit none

  intrinsic :: dble, real

! Constant parameters
  integer, parameter :: &
  nlen = 1, &   ! Length of the total domain
  slen = 1      ! Length of the sub domain

! Number of levels to take from U.S. Std. Atmos tables
  integer, parameter :: std_atmos_buffer = 10 ! For typical cases

! Input Variables
  real, intent(in) :: &
  lat_in_degrees,&! Latitude   [Degrees North]
  lon_in_degrees  ! Longitude  [Degrees East]

  real(kind=time_precision), intent(in) :: &
  time ! Model time [s]
  
  integer, intent(in) :: &
  nz,              & ! Vertical extent;  i.e. nnzp in the grid class
  day, month, year   ! Date of model start

! Number of levels to interpolate from the bottom of std_atmos to the top
! of the HOC profile, hopefully enough to eliminate cooling spikes, etc. 
  integer, intent(in) :: lin_int_buffer

  real, intent(in), dimension(nz) :: &
  alt,   & ! Altitudes of the model     [m]
  thlm,  & ! Liquid potential temp.     [K]
  rcm,   & ! Liquid water mixing ratio  [kg/kg]
  rsnwm, & ! Snow water mixing ratio    [kg/kg]
  rim,   & ! Ice water mixing ratio     [kg/kg]
  rtm,   & ! Total water mixing ratio   [kg/kg]
  rhom,  & ! Density                    [kg/m^3]
  cf,    & ! Cloud fraction             [%]
  p_in_Pa, & ! Pressure on the t grid     [Pa]
  p_in_Pam,& ! Pressure on the m grid     [Pa]
  exner    ! Exner function             [-]

! Input/Output Variables
  real, intent(inout), dimension(nz) :: &
  thlm_forcing ! Theta_l LS tendency [K/s]

! Output Variables
  real, intent(out), dimension(nz) :: &
  Frad, & ! Total radiative flux       [W/m^2]
  radht   ! Total heating rate         [K/s]

! Local Variables
  real, dimension(nz) :: &
  Frad_SW, & ! SW radiative flux       [W/m^2]
  Frad_LW, & ! LW radiative flux       [W/m^2]
  radht_SW,& ! SW heating rate         [K/s]
  radht_LW   ! LW heating rate         [K/s]

! Altered 3 Oct 2005 to be buffer levels higher
  double precision, dimension(nlen,(nz-1)+lin_int_buffer+std_atmos_buffer) :: &
  T_in_K,& ! Temperature            [K]
  rcil, & ! Ice mixing ratio       [kg/kg]
  o3l     ! Ozone mixing ratio     [kg/kg]

  double precision, dimension(nlen,(nz-1)+lin_int_buffer+std_atmos_buffer+1) :: &
  Frad_uLW, & ! LW upwelling flux         [W/m^2]
  Frad_dLW, & ! LW downwelling flux       [W/m^2]
  Frad_uSW, & ! SW upwelling flux         [W/m^2]
  Frad_dSW    ! SW downwelling flux       [W/m^2]

  double precision, dimension(nlen,(nz-1)+lin_int_buffer+std_atmos_buffer) :: &
  sp_humidity, & ! Specific humidity      [kg/kg]
  pinmb          ! Pressure in millibars  [hPa]

! Pressure in millibars for layers (calculated as an average of pinmb)
  double precision, dimension(nlen,(nz-1)+lin_int_buffer+std_atmos_buffer+1) :: &
  playerinmb ! [hPa]

  double precision, dimension(nlen,(nz-1)+lin_int_buffer+std_atmos_buffer) :: &
  dpl, &          ! Difference in pressure levels       [hPa]
  rsnwm2, rcm2, cf2 ! Two-dimensional copies of the input parameters

  double precision, dimension(nlen,(nz-1)+lin_int_buffer+std_atmos_buffer+1) :: &
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

  double precision :: z1_fact, z2_fact, tmp ! Temp storage

  integer :: i, j, z, z1, z2  ! Loop indices

  integer :: buffer ! The sum of the two buffers
!-----------------------------------------------------------------------

  buffer = lin_int_buffer + std_atmos_buffer

! amu0 = 0.4329 ! Nov 11 Altocu value

! Calculated value of cosine of the solar zenith angle
  amu0 = cos_solar_zen( day, month, year, time, lat_in_degrees, lon_in_degrees )

! Convert to millibars
  pinmb(1,1:(nz-1))  = dble( p_in_Pa(2:nz) / 100.0 ) ! t grid in HOC

  playerinmb(1,1:nz) = dble( p_in_Pam / 100.0 ) ! m grid in HOC


! Convert theta_l to temperature
  
  T_in_K(1,1:(nz-1)) = thlm2T_in_K( thlm(2:nz), exner(2:nz), rcm(2:nz) )
  
! Derive Specific humidity from rc & rt.
  do z = 2, nz
    if ( rtm(z) < rcm(z) ) then
      sp_humidity(1,z-1) = 0.0d0
      if ( clubb_at_debug_level(1) ) then
      write(fstderr,*) "rvm < 0 at ", z, " before BUGSrad, specific humidity set to 0."
      endif
    else
      sp_humidity(1,z-1) &
        = dble( rtm(z) - rcm(z) ) / dble( 1.0+rtm(z) )
     end if
   end do

! Setup miscellaneous variables

! Albedo values
  alvdr = 0.1d0
  alvdf = 0.1d0
  alndr = 0.1d0
  alndf = 0.1d0

  slr  = 1.0d0 ! Fraction of daylight

! Ozone at < 1 km = 5.4e-5 g/m^3 from U.S. Standard Atmosphere, 1976. 
!   Convert from g to kg.
  o3l(1,1:(nz-1)) = dble( ( 5.4e-5 / rhom(1:(nz-1)) ) * 0.001 )

! Convert and transpose as needed
  rcil(1,buffer+1:(nz-1)+buffer)  = flip( dble( rim(2:nz) ), nz-1 )
  rsnwm2(1,buffer+1:(nz-1)+buffer)= flip( dble( rsnwm(2:nz) ), nz-1 )
  rcm2(1,buffer+1:(nz-1)+buffer)  = flip( dble( rcm(2:nz) ), nz-1 )
  cf2(1,buffer+1:(nz-1)+buffer)   = flip( dble( cf(2:nz) ), nz-1 ) 

  T_in_K(1,buffer+1:(nz-1)+buffer) = flip( T_in_K(1,1:(nz-1)), nz-1 )

  sp_humidity(1,buffer+1:(nz-1)+buffer) = flip( sp_humidity(1,1:(nz-1)), nz-1 )

  pinmb(1,(buffer+1):(nz-1+buffer))        = flip( pinmb(1,1:(nz-1)), nz-1 )
  playerinmb(1,(buffer+1):(nz-1+buffer+1)) = flip( playerinmb(1,1:nz), nz )

  o3l(1,buffer+1:(nz-1)+buffer) = flip( o3l(1,1:(nz-1)), nz-1 )

! Assume these are all zero above the HOC profile
  rsnwm2(1,1:buffer) = 0.0d0
  rcil(1,1:buffer)   = 0.0d0
  rcm2(1,1:buffer)   = 0.0d0
  cf2(1,1:buffer)    = 0.0d0

! Determine the altitude of the HOC model compared to the standard atmospheric 
! profile.  Then take the values at altitude j km on the table and use 
! std_atmos_buffer kilometers of it as the top of the profile fed into BUGsrad.

! e.g.
! if the HOC maximum altitude = 3200 m, then
! the nearest layer above in the std_atmos arrays is 4000 m, so
! the array used in BUGSrad will be:
! HOC layers up to 3200 m (grid spacing determined by model) +
! lin_int_buffer layers between 3200 m and 4000 m (variable grid spacing) +
! std_atmos_buffer layers from 4000 m to 14000 m  (1 km grid spacing)

  j = 1 ! initial altitude
  do while ( ( std_alt(j) ) < alt(nz) )
    j = j + 1
    if ( (j + std_atmos_buffer ) > std_atmos_dim ) then
      write(fstderr,*) "j = ", j, "alt = ", alt(nz), " m"
      stop "bugsrad_hoc: cannot handle this altitude" ! exceeds a 50 km altitude
    end if
  end do

  ! Add the standard atmospheric profile above the linear interpolation
  do i = 1, std_atmos_buffer, 1
    T_in_K(1,i)       = std_T_in_K((std_atmos_buffer+j)-(i-1)) 
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
    z1_fact = dble( z2 - z ) / dble( z2 - z1 ) 
    z2_fact = dble( z - z1 ) / dble( z2 - z1 ) 

    T_in_K(1,z) = z1_fact * T_in_K(1,z1) + z2_fact * T_in_K(1,z2)

    sp_humidity(1,z) = z1_fact * sp_humidity(1,z1) + z2_fact * sp_humidity(1,z2)

    o3l(1,z) = z1_fact * o3l(1,z1) + z2_fact * o3l(1,z2)

    pinmb(1,z) = z1_fact * pinmb(1,z1) + z2_fact * pinmb(1,z2)
  end do

  ! Do a linear interpolation to find playerinmb.  Since this interpolation
  ! occurs at levels above the top of the CLUBB model, the CLUBB zt2zm function 
  ! or CLUBB weighted averages do not apply.  The variable playerinmb is being
  ! defined on momentum levels above the top of the CLUBB model, which are 
  ! being defined here at points half-way inbetween the thermodynamic levels
  ! above the top of the CLUBB model.  Brian Griffin; May 13, 2008.
  playerinmb(1,2:buffer+1) = ( pinmb(1,1:buffer) + pinmb(1,2:buffer+1) ) / 2.

  ! Do a linear extension to find playerinmb at the uppermost standard 
  ! atmosphere momentum level.  The grid is evenly-spaced at these points.
  ! Brian Griffin; May 13, 2008.
  tmp = 2. * playerinmb(1,2) - playerinmb(1,3)
  if ( tmp > 0. ) then
    playerinmb(1,1) = tmp
  else ! Assuming a linear extension didn't work
    playerinmb(1,1) = .5 * playerinmb(1,2)
  end if

! Calculate the difference in pressure layers (including buffer levels)
  do i = 1, (nz-1)+buffer
    dpl(1,i) = playerinmb(1,i+1) - playerinmb(1,i)
  end do

  ts(1) = T_in_K(1,(nz-1)+buffer)

!  Write a profile for Kurt's driver program for debugging purposes
! open(10, file="profile.dat")
! write(10,'(2i4,a10)') nlen, (nz-1)+buffer, "TROPICAL"
! do i=1, (nz-1)+buffer
!   write(10,'(i4,9f12.6)') i, pinmb(1,i), playerinmb(1,i),T_in_K(1,i),         &
!   sp_humidity(1,i), 100000.0*o3l(1,i), rcm2(1,i), rcil(1,i),cf2(1,i),dpl(1,i)
! enddo
! write(10,'(a4,a12,3f12.6)') "","", playerinmb(1,nz+buffer), ts(1), amu0(1)
! close(10)
! pause
 
  call bugs_rad( nlen, slen, (nz-1)+buffer, playerinmb,          &
                 pinmb, dpl, T_in_K, sp_humidity,                 &
                 rcm2, rcil, rsnwm2, o3l,                        &
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

  if ( l_stats_samp ) then
 
    call stat_update_var( iradht_LW, radht_LW, zt )   

    call stat_update_var( iradht_SW, radht_SW, zt )

    call stat_update_var( iFrad_SW, Frad_SW, zm )

    call stat_update_var( iFrad_LW, Frad_LW, zm )

  end if ! lstats_samp
 

  return
end subroutine bugsrad_hoc
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
function flip( x, xdim )
! Description:
! Flips a single dimension array (i.e. a vector), so the first element 
! becomes the last and vice versa for the whole column.  This is a 
! necessary part of the code because BUGSrad and HOC store altitudes in
! reverse order
!-----------------------------------------------------------------------
  implicit none

  ! Input
  integer, intent(in) :: xdim

  double precision, dimension(xdim), intent(in) :: x

  ! Output
    double precision, dimension(xdim) :: flip

  ! Internal
  double precision, dimension(xdim) :: tmp
  integer :: indx

  do indx = 1, xdim, 1
    tmp(indx) = x((xdim+1) - (indx))
  end do

  flip = tmp

  return
end function flip
!-----------------------------------------------------------------------

end module bugsrad_hoc_mod
