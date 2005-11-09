!-----------------------------------------------------------------------
! $Id: bugsrad_hoc.F90,v 1.2 2005-11-09 22:42:02 dschanen Exp $
! SUBROUTINE bugsrad_hoc
! Does the necessary operations to interface the HOC model with
! the bugsrad subprogram.
!-----------------------------------------------------------------------
subroutine bugsrad_hoc( alt, nz, thlm, rcm, rtm, rrm, cf, pinpa, rhom, Tsfc,   &
                        radht, radht_SW, radht_LW, Frad, Frad_SW, Frad_LW,     &
                        thlm_forcing )
  use constants

  implicit none

! Parameters
  integer, parameter :: nlen = 1 ! length of the total domain
  integer, parameter :: slen = 1 ! length of the sub domain

! Number of levels to take from U.S. Std. Atmos tables
  integer, parameter :: std_atmos_buffer = 10 

! Number of levels to interpolate from the bottom of std_atmos to the top
! of the HOC profile, hopefully enough to elminate cooling spikes, etc. 
  integer, parameter :: lin_int_buffer = 20

! The sum of the above to two
  integer, parameter :: buffer = lin_int_buffer + std_atmos_buffer

  integer, parameter :: std_atmos_dim = 25

! Parameters from U.S. Standard Atmosphere, 1976;  Starting at 1 km altitude
! Pressure in millibars
  double precision, parameter, dimension(std_atmos_dim) ::                     &
  std_pinmb = (/8.986e2, 7.950e2, 7.012e2, 6.166e2, 5.405e2,                   &
                4.722e2, 4.111e2, 3.565e2, 3.080e2, 2.650e2,                   &
                2.270e2, 1.940e2, 1.658e2, 1.417e2, 1.211e2,                   &
                1.035e2, 8.850e1, 7.565e1, 6.467e1, 5.529e1,                   &
                4.729e1, 4.047e1, 3.467e1, 2.972e1, 2.546e1/)

! Temperature in degrees Kelvin
  double precision, parameter, dimension(std_atmos_dim) ::                     &
  std_tempk = (/281.6, 275.1, 268.7, 262.2, 255.7,                             &
                249.2, 242.7, 236.2, 229.7, 223.2,                             &
                216.8, 216.6, 216.6, 216.6, 216.6,                             &
                216.6, 216.6, 216.6, 216.6, 216.6,                             &
                217.6, 218.6, 219.6, 220.6, 221.6/)

! Specific Humidity ( Water Vapor / Density )
  double precision, parameter, dimension(std_atmos_dim) ::                     &
  std_sp_hmdty = (/0.378e-02, 0.288e-02, 0.198e-02, 0.134e-02, 0.869e-03,      &
                   0.576e-03, 0.356e-03, 0.228e-03, 0.985e-04, 0.435e-04,      &
                   0.225e-04, 0.119e-04, 0.675e-05, 0.369e-05, 0.370e-05,      &
                   0.366e-05, 0.365e-05, 0.362e-05, 0.423e-05, 0.495e-05,      &
                   0.634e-05, 0.806e-05, 0.104e-04, 0.130e-04, 0.165e-04/)

! Ozone ( O_3 / Density )
  double precision, parameter, dimension(std_atmos_dim) ::                     &
  std_o3l   = (/0.486e-07, 0.536e-07, 0.550e-07, 0.561e-07, 0.611e-07,         &
                0.682e-07, 0.814e-07, 0.989e-07, 0.152e-06, 0.218e-06,         &
                0.356e-06, 0.513e-06, 0.638e-06, 0.834e-06, 0.108e-05,         &
                0.138e-05, 0.197e-05, 0.263e-05, 0.337e-05, 0.427e-05,         &
                0.502e-05, 0.605e-05, 0.691e-05, 0.767e-05, 0.848e-05/)

! Input
  real, intent(in)    :: alt ! maximum altitude in the model domain (kilometers)
  integer, intent(in) :: nz  ! nnzp in the grid class
  
  real, intent(in), dimension(nz) :: thlm  ! Liquid potential temperature (K)
  real, intent(in), dimension(nz) :: rcm   ! Liquid water mixing ratio (kg/kg)
  real, intent(in), dimension(nz) :: rrm   ! Rain water mixing ratio (kg/kg)
  real, intent(in), dimension(nz) :: rtm   ! Total water mixing ratio (kg/kg)
  real, intent(in), dimension(nz) :: rhom  ! Rho m
  real, intent(in), dimension(nz) :: cf    ! Cloud fraction (%)
  real, intent(in), dimension(nz) :: pinpa ! Pressure in pascals

  real, intent(in) :: Tsfc     ! Temperature in K (surface)

! Output
  real, intent(out), dimension(nz) :: Frad    ! Total radiative flux
  real, intent(out), dimension(nz) :: Frad_SW ! SW radiative flux
  real, intent(out), dimension(nz) :: Frad_LW ! LW radiative flux

  real, intent(out), dimension(nz) :: radht    ! Total heating rate
  real, intent(out), dimension(nz) :: radht_SW ! SW heating rate
  real, intent(out), dimension(nz) :: radht_LW ! LW heating rate

  real, intent(inout), dimension(nz) :: thlm_forcing

! Internal
! Altered 3 Oct 2005 to be buffer levels higher
  double precision, dimension(nlen,(nz-1)+buffer) :: tempk! Temperature in K
  double precision, dimension(nlen,(nz-1)+buffer) :: rcil ! Ice mixing ratio (kg/kg)
  double precision, dimension(nlen,(nz-1)+buffer) :: o3l  ! Ozone mixing ratio (kg/kg)

  double precision, dimension(nlen,(nz-1)+buffer+1) :: Frad_uLW ! LW upwelling flux (in W/m^2)
  double precision, dimension(nlen,(nz-1)+buffer+1) :: Frad_dLW ! LW downwelling flux (in W/m^2)
  double precision, dimension(nlen,(nz-1)+buffer+1) :: Frad_uSW ! SW upwelling flux (in W/m^2)
  double precision, dimension(nlen,(nz-1)+buffer+1) :: Frad_dSW ! SW downwelling flux (in W/m^2)

  double precision, dimension(nlen,(nz-1)+buffer)   :: sp_humidity ! Specific humidity (kg/kg)
  double precision, dimension(nlen,(nz-1)+buffer)   :: pinmb       ! Pressure in millibars
  double precision, dimension(nlen,(nz-1)+buffer+1) :: playerinmb  ! Pressure in millibars for layers (calculated as an average of pinmb)
  double precision, dimension(nlen,(nz-1)+buffer)   :: dpl         ! Difference in pressure levels (mb)

  double precision, dimension(nlen,(nz-1)+buffer) :: rrm2, rcm2, cf2 ! Two-dimensional copies of the input parameters

  double precision, dimension(nlen,(nz-1)+buffer+1) :: radht_SW2    ! SW Radiative heating rate
  double precision, dimension(nlen,(nz-1)+buffer+1) :: radht_LW2    ! LW Radiative heating rate

  double precision alvdr(nlen) ! Visible direct surface albedo  
  double precision alvdf(nlen) ! Visible diffuse surface albedo  
  double precision alndr(nlen) ! Near-IR direct surface albedo  
  double precision alndf(nlen) ! Near-IR diffuse surface albedo  

  double precision slr(nlen)  ! Fraction of daylight  
  double precision ts(nlen)   ! Surface temperature in K

  integer i, j, z, z1, z2  ! loop indices
  double precision z1_fact, z2_fact

! Convert to millibars
  pinmb(1,1:(nz-1))    = dble( pinpa(2:nz) / 100.0 ) ! t grid in HOC

  playerinmb(1,2:nz-1) = dble( (pinmb(1,1:(nz-2)) + pinmb(1,2:nz-1))/ 2 ) ! m grid in HOC
  playerinmb(1,1)      = ( dble( pinpa(1) / 100.0 ) + pinmb(1,2) ) / 2
  playerinmb(1,nz)     = 2 * playerinmb(1,nz-1) - playerinmb(1,nz-2)


! Convert theta_l to temperature
!   kappa: Dry air gas constant / Dry air specific heat at p
!   Lv:    Latent heat of vaporization
  tempk(1,1:(nz-1)) = thlm(2:nz) * ( 1000.0d0 / pinmb(1,1:(nz-1)) )**(-kappa)  &
                    + Lv*rcm(2:nz) / Cp

! Derive Specific humidity from rc & rt.
  sp_humidity(1,1:(nz-1)) = dble( rtm(2:nz) - rcm(2:nz) ) / dble( 1+rtm(2:nz) )

! Setup misc. variables
  alvdr = 0.1d0
  alvdf = 0.1d0
  alndr = 0.1d0
  alndf = 0.1d0

  slr  = 1.0d0 ! fraction of daylight
 
  rcil(1,1:(nz-1)+buffer) = 0.0d0 ! Assume no ice for HOC

! Ozone = 5.4e-5 g/m^3 from U.S. Standard Atmosphere, 1976. 
!   Convert from g to kg.
  o3l(1,1:(nz-1)) = dble( ( 5.4e-5 / rhom(1:(nz-1)) ) * 0.001 )

! Convert and transpose as needed
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
! write(10,'(a4,a12,3f12.6)') "","", playerinmb(1,nz+buffer), ts(1), amu0
! close(10)
! pause

  call bugs_rad( nlen, slen, (nz-1)+buffer, playerinmb, pinmb, dpl, tempk,     &
                 sp_humidity, rcm2, rcil, rrm2, o3l, ts, dble( amu0 ),         &
                 slr, alvdf, alndf, alvdr, alndr,                              &
                 dble(sol_const), dble(grav), dble(Cp),                        &
                 radht_SW2, radht_LW2,                                         &
                 Frad_dSW, Frad_uSW, Frad_dLW, Frad_uLW, cf2 )

  radht_SW(2:nz) = real( flip( radht_SW2(1,buffer+1:nz+buffer), nz-1 ) )
  radht_LW(2:nz) = real( flip( radht_LW2(1,buffer+1:nz+buffer), nz-1 ) )

! No radiative heating below ground
  radht_SW(1) = 0.0
  radht_LW(1) = 0.0

  radht = radht_SW + radht_LW

! These are on the m grid, and require no adjusting
  Frad_SW(1:nz) = real( flip(  Frad_uSW(1,buffer+1:nz+buffer)                  &
                             - Frad_dSW(1,buffer+1:nz+buffer), nz ) )
  Frad_LW(1:nz) = real( flip(  Frad_uLW(1,buffer+1:nz+buffer)                  &
                             - Frad_dLW(1,buffer+1:nz+buffer), nz ) )

  Frad(1:nz) = Frad_SW(1:nz) + Frad_LW(1:nz)

  thlm_forcing(1:nz) = thlm_forcing(1:nz) + radht(1:nz)

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
    enddo

    flip = tmp

    return
  end function flip
!-----------------------------------------------------------------------

end subroutine bugsrad_hoc
!-----------------------------------------------------------------------
