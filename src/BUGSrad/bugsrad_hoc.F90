#define BFFR 20 /* additional z-levels to eliminate the cooling spike */

!-----------------------------------------------------------------------
! $Id: bugsrad_hoc.F90,v 1.1 2005-10-27 20:06:50 dschanen Exp $
! SUBROUTINE bugsrad_hoc
! Does the necessary operations to interface the HOC model with
! the bugsrad subprogram.  
!-----------------------------------------------------------------------
subroutine bugsrad_hoc( nz, thlm, rcm, rtm, rrm, cf, pinpa, rhom, Tsfc,        &
                        radht, radht_SW, radht_LW, Frad, Frad_SW, Frad_LW,     &
                        thlm_forcing )
  use constants

  implicit none

! Parameters
  integer, parameter :: nlen = 1  ! Single column

! Input
  integer, intent(in) :: nz
  
  real, intent(in), dimension(nz) :: thlm  ! Liquid potential temperature (K)
  real, intent(in), dimension(nz) :: rcm   ! Liquid water mixing ratio
  real, intent(in), dimension(nz) :: rrm   ! Rain water mixing ratio
  real, intent(in), dimension(nz) :: rtm   ! Total water mixing ratio
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
! Altered 3 Oct 2005 to be BFFR levels higher
  double precision, dimension(nlen,(nz-1)+BFFR) :: tempk! Temperature in K
  double precision, dimension(nlen,(nz-1)+BFFR) :: rcil ! Ice mixing ratio (kg/kg)
  double precision, dimension(nlen,(nz-1)+BFFR) :: o3l  ! Ozone mixing ratio (kg/kg)

  double precision, dimension(nlen,(nz-1)+BFFR+1) :: Frad_uLW ! LW upwelling flux (in W/m^2)
  double precision, dimension(nlen,(nz-1)+BFFR+1) :: Frad_dLW ! LW downwelling flux (in W/m^2)
  double precision, dimension(nlen,(nz-1)+BFFR+1) :: Frad_uSW ! SW upwelling flux (in W/m^2)
  double precision, dimension(nlen,(nz-1)+BFFR+1) :: Frad_dSW ! SW downwelling flux (in W/m^2)

  double precision, dimension(nlen,(nz-1)+BFFR)   :: sp_humidity ! Specific humidity (kg/kg)
  double precision, dimension(nlen,(nz-1)+BFFR)   :: pinmb       ! Pressure in millibars
  double precision, dimension(nlen,(nz-1)+BFFR+1) :: playerinmb  ! Pressure in millibars for layers (calculated as an average of pinmb)
  double precision, dimension(nlen,(nz-1)+BFFR)   :: dpl         ! Difference in pressure levels (mb)

  double precision, dimension(nlen,(nz-1)+BFFR) :: rrm2, rcm2, cf2 ! Two-dimensional copies of the input parameters

  double precision, dimension(nlen,(nz-1)+BFFR+1) :: radht_SW2    ! SW Radiative heating rate
  double precision, dimension(nlen,(nz-1)+BFFR+1) :: radht_LW2    ! LW Radiative heating rate

  double precision alvdr(nlen) ! Visible direct surface albedo  
  double precision alvdf(nlen) ! Visible diffuse surface albedo  
  double precision alndr(nlen) ! Near-IR direct surface albedo  
  double precision alndf(nlen) ! Near-IR diffuse surface albedo  

  double precision slr(nlen)  ! Fraction of daylight  
  double precision ts(nlen)   ! Surface temperature in K

  integer i, j, l   ! loop indices

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
 
  rcil(1,1:(nz-1)+BFFR) = 0.0d0 ! Assume no ice for HOC

! Ozone = 5.4e-5 g/m^3 from U.S. Standard Atmosphere, 1976. 
!   Convert from g to kg.
  o3l(1,1:(nz-1)) = dble( ( 5.4e-5 / rhom(1:(nz-1)) ) * 0.001 )

! Convert and transpose as needed
  rrm2(1,BFFR+1:(nz-1)+BFFR)  = flip( dble( rrm(2:nz) ), nz-1 )
  rcm2(1,BFFR+1:(nz-1)+BFFR)  = flip( dble( rcm(2:nz) ), nz-1 )
  cf2(1,BFFR+1:(nz-1)+BFFR)   = flip( dble( cf(2:nz) ), nz-1 ) 

  tempk(1,BFFR+1:(nz-1)+BFFR) = flip( tempk(1,1:(nz-1)), nz-1 )

  sp_humidity(1,BFFR+1:(nz-1)+BFFR) = flip( sp_humidity(1,1:(nz-1)), nz-1 )

  pinmb(1,(BFFR+1):(nz-1+BFFR))        = flip( pinmb(1,1:(nz-1)), nz-1 )
  playerinmb(1,(BFFR+1):(nz-1+BFFR+1)) = flip( playerinmb(1,1:nz), nz )

  o3l(1,BFFR+1:(nz-1)+BFFR) = flip( o3l(1,1:(nz-1)), nz-1 )

! Add BFFR levels
  rrm2(1,1:BFFR) = rrm2(1,BFFR+1)
  rcm2(1,1:BFFR) = rcm2(1,BFFR+1)
  cf2(1,1:BFFR)  = cf2(1,BFFR+1)

  tempk(1,1:BFFR) = tempk(1,BFFR+1)

  sp_humidity(1,1:BFFR) = sp_humidity(1,BFFR+1)

  o3l(1,1:BFFR) = o3l(1,BFFR+1)

! Linear interpolation to generate pressure
  do i = BFFR, 1, -1 ! Buffer...1
    pinmb(1,i) = 2 * pinmb(1,i+1) - pinmb(1,i+2)
    playerinmb(1,i) = 2 * playerinmb(1,i+1) - playerinmb(1,i+2)
  enddo


! Calculate the difference in pressure layers (including buffer levels)
  do i = 1, (nz-1)+BFFR
    dpl(1,i) = playerinmb(1,i+1) - playerinmb(1,i)
  enddo

!  ts(1) = dble( Tsfc )
  ts(1) = tempk(1,(nz-1)+BFFR)

!  Write a profile for Kurt's driver program
!  open(10, file="profile.dat")
!  write(10,'(2i4,a10)') nlen, (nz-1)+BFFR, "TROPICAL"
!  do i=1, (nz-1)+BFFR
!    write(10,'(i4,9f12.6)') i, pinmb(1,i), playerinmb(1,i),tempk(1,i),         &
!    sp_humidity(1,i), 100000.0*o3l(1,i), rcm2(1,i), rcil(1,i),cf2(1,i),dpl(1,i)
!  enddo
!  write(10,'(a4,a12,3f12.6)') "","", playerinmb(1,nz+BFFR), ts(1), amu0
!  close(10)

  call bugs_rad( nlen, nlen, (nz-1)+BFFR, playerinmb, pinmb, dpl, tempk,       &
                 sp_humidity, rcm2, rcil, rrm2, o3l, ts, dble( amu0 ),         &
                 slr, alvdf, alndf, alvdr, alndr,                              &
                 dble(sol_const), dble(grav), dble(Cp),                        &
                 radht_SW2, radht_LW2,                                         &
                 Frad_dSW, Frad_uSW, Frad_dLW, Frad_uLW, cf2 )

  radht_SW(2:nz) = real( flip( radht_SW2(1,BFFR+1:nz+BFFR), nz-1 ) )
  radht_LW(2:nz) = real( flip( radht_LW2(1,BFFR+1:nz+BFFR), nz-1 ) )

! No radiative heating below ground
  radht_SW(1) = 0.0
  radht_LW(1) = 0.0

  radht    = radht_SW + radht_LW

! These are on the m grid, and require no adjusting
  Frad_SW(1:nz) = real( flip(  Frad_uSW(1,BFFR+1:nz+BFFR)                      &
                             - Frad_dSW(1,BFFR+1:nz+BFFR), nz ) )
  Frad_LW(1:nz) = real( flip(  Frad_uLW(1,BFFR+1:nz+BFFR)                      &
                             - Frad_dLW(1,BFFR+1:nz+BFFR), nz ) )

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
