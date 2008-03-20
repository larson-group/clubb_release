!-----------------------------------------------------------------------
! $Id: bugsrad_hoc.F90,v 1.19 2008-03-20 20:20:33 dschanen Exp $

subroutine bugsrad_hoc( alt, nz, lat_in_degrees, lon_in_degrees, &
                        day, month, year, time,                  &
                        thlm, rcm, rtm, rrm, rim,                & 
                        cf, pinpa, exner, rhom,                  &
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

  use constants
#ifdef STATS
use stats_hoc, only: zt, zm, lstats_samp, &
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
  integer, parameter :: std_atmos_buffer = 10 ! For typical cases

! Number of levels to interpolate from the bottom of std_atmos to the top
! of the HOC profile, hopefully enough to eliminate cooling spikes, etc. 
  integer, parameter :: lin_int_buffer = 20

! The sum of the above to two buffers
  integer, parameter :: buffer = lin_int_buffer + std_atmos_buffer

  integer, parameter :: std_atmos_dim = 50

! Parameters from U.S. Standard Atmosphere, 1976;  Starting at 1 km altitude

! Altitude in meters
  double precision, parameter, dimension(std_atmos_dim) :: &
  std_alt = (/  &
    1000.00,       2000.00,       3000.00,       4000.00,       5000.00,  &
    6000.00,       7000.00,       8000.00,       9000.00,       10000.0,  &
    11000.0,       12000.0,       13000.0,       14000.0,       15000.0,  &
    16000.0,       17000.0,       18000.0,       19000.0,       20000.0,  &
    21000.0,       22000.0,       23000.0,       24000.0,       25000.0,  &
    26000.0,       27000.0,       28000.0,       29000.0,       30000.0,  &
    31000.0,       32000.0,       33000.0,       34000.0,       35000.0,  &
    36000.0,       37000.0,       38000.0,       39000.0,       40000.0,  &
    41000.0,       42000.0,       43000.0,       44000.0,       45000.0,  &
    46000.0,       47000.0,       48000.0,       49000.0,       50000.0   /)

! Pressure in millibars
  double precision, parameter, dimension(std_atmos_dim) :: &
  std_pinmb = (/  &
    898.600,       795.000,       701.200,       616.600,       540.500,  &
    472.200,       411.100,       356.500,       308.000,       265.000,  &
    227.000,       194.000,       165.800,       141.700,       121.100,  &
    103.500,       88.5000,       75.6500,       64.6700,       55.2900,  &
    47.2900,       40.4700,       34.6700,       29.7200,       25.4600,  &
    22.7620,       20.0640,       17.3660,       14.6680,       11.9700,  &
    10.7252,       9.48040,       8.23560,       6.99080,       5.74600,  &
    5.17100,       4.59600,       4.02100,       3.44600,       2.87100,  &
    2.59500,       2.31900,       2.04300,       1.76700,       1.49100,  &
    1.35236,       1.21372,       1.07508,      0.936440,      0.797800   /)

! Temperature in degrees Kelvin
  double precision, parameter, dimension(std_atmos_dim) :: &
  std_tempk = (/  &
    281.600,       275.100,       268.700,       262.200,       255.700,  &
    249.200,       242.700,       236.200,       229.700,       223.200,  &
    216.800,       216.600,       216.600,       216.600,       216.600,  &
    216.600,       216.600,       216.600,       216.600,       216.600,  &
    217.600,       218.600,       219.600,       220.600,       221.600,  &
    222.580,       223.560,       224.540,       225.520,       226.500,  &
    228.500,       230.500,       232.500,       234.500,       236.500,  &
    239.280,       242.060,       244.840,       247.620,       250.400,  &
    253.160,       255.920,       258.680,       261.440,       264.200,  &
    265.480,       266.760,       268.040,       269.320,       270.600   /)

! Specific Humidity ( Water Vapor / Density )
  double precision, parameter, dimension(std_atmos_dim) :: &
  std_sp_hmdty = (/ &
   0.378038E-02,  0.287984E-02,  0.197954E-02,  0.134261E-02,  0.869093E-03, &
   0.575670E-03,  0.355932E-03,  0.228224E-03,  0.984800E-04,  0.435308E-04, &
   0.224781E-04,  0.118628E-04,  0.675169E-05,  0.368583E-05,  0.369610E-05, &
   0.366366E-05,  0.365425E-05,  0.361842E-05,  0.423077E-05,  0.494882E-05, &
   0.633914E-05,  0.806077E-05,  0.103636E-04,  0.129953E-04,  0.164671E-04, &
   0.173018E-04,  0.181366E-04,  0.189714E-04,  0.198062E-04,  0.206410E-04, & 
   0.202939E-04,  0.199469E-04,  0.195999E-04,  0.192529E-04,  0.189058E-04, &
   0.184780E-04,  0.180502E-04,  0.176224E-04,  0.171946E-04,  0.167668E-04, &
   0.166688E-04,  0.165707E-04,  0.164727E-04,  0.163747E-04,  0.162767E-04, &
   0.153583E-04,  0.144398E-04,  0.135214E-04,  0.126030E-04,  0.116845E-04  /)


! Ozone ( O_3 / Density )
  double precision, parameter, dimension(std_atmos_dim) :: & 
  std_o3l= (/ & 
   0.486049E-07,  0.536246E-07,  0.549874E-07,  0.561455E-07,  0.611081E-07, &
   0.681715E-07,  0.813559E-07,  0.988969E-07,  0.152002E-06,  0.217654E-06, &
   0.356360E-06,  0.512985E-06,  0.637659E-06,  0.833699E-06,  0.107803E-05, &
   0.138138E-05,  0.196767E-05,  0.263158E-05,  0.336538E-05,  0.427398E-05, &
   0.501849E-05,  0.604557E-05,  0.690909E-05,  0.766937E-05,  0.848303E-05, &
   0.895916E-05,  0.943528E-05,  0.991141E-05,  0.103875E-04,  0.108637E-04, &
   0.112905E-04,  0.117173E-04,  0.121441E-04,  0.125709E-04,  0.129978E-04, &
   0.128507E-04,  0.127036E-04,  0.125565E-04,  0.124094E-04,  0.122623E-04, &
   0.115392E-04,  0.108162E-04,  0.100931E-04,  0.937005E-05,  0.864700E-05, &
   0.769657E-05,  0.674614E-05,  0.579570E-05,  0.484527E-05,  0.389484E-05  /)

! Input Variables
  real, intent(in) :: &
  lat_in_degrees,&! Latitude   [Degrees North]
  lon_in_degrees,&! Longitude  [Degrees East]
  time            ! Model time [s]
  
  integer, intent(in) :: &
  nz,              & ! Vertical extent;  i.e. nnzp in the grid class
  day, month, year   ! Time of year

  real, intent(in), dimension(nz) :: &
  alt,   & ! Altitudes of the model     [m]
  thlm,  & ! Liquid potential temp.     [K]
  rcm,   & ! Liquid water mixing ratio  [kg/kg]
  rrm,   & ! Rain water mixing ratio    [kg/kg]
  rim,   & ! Ice water mixing ratio     [kg/kg]
  rtm,   & ! Total water mixing ratio   [kg/kg]
  rhom,  & ! Density                    [kg/m^3]
  cf,    & ! Cloud fraction             [%]
  pinpa, & ! Pressure                   [Pa]
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

  double precision :: z1_fact, z2_fact, tmp ! Temp storage

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
  do z = 2, nz
    if ( rtm(z) < rcm(z) ) then
      sp_humidity(1,z-1) = 0.0d0
      write(fstderr,*) "rvm < 0 at ", z, " before BUGSrad, specific humidity set to 0."
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
    z1_fact = dble( z2 - z ) / dble( z2 - z1 ) 
    z2_fact = dble( z - z1 ) / dble( z2 - z1 ) 

    tempk(1,z) = z1_fact * tempk(1,z1) + z2_fact * tempk(1,z2)

    sp_humidity(1,z) = z1_fact * sp_humidity(1,z1) + z2_fact * sp_humidity(1,z2)

    o3l(1,z) = z1_fact * o3l(1,z1) + z2_fact * o3l(1,z2)

    pinmb(1,z) = z1_fact * pinmb(1,z1) + z2_fact * pinmb(1,z2)
  end do

  playerinmb(1,2:buffer+1) = ( pinmb(1,1:buffer) + pinmb(1,2:buffer+1) ) / 2.

  tmp = 2. * playerinmb(1,2) - playerinmb(1,3)
  if ( tmp > 0. ) then
    playerinmb(1,1) = tmp
  else ! Assuming a linear extrapolation didn't work
    playerinmb(1,1) = .5 * playerinmb(1,2)
  end if

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

  end if ! lstats_samp
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
