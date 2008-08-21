! $Id: diag_variables.F90,v 1.15 2008-08-17 16:10:40 griffinb Exp $
module diagnostic_variables

! This module contains definitions of all diagnostic
! arrays used in the single column model, as well as subroutines
! to allocate, deallocate and initialize them.

! Note that while these are all same dimension, there is a
! thermodynamic and momentum grid and they have different levels
!-----------------------------------------------------------------------

implicit none

private ! Set default scope

public :: setup_diagnostic_variables, & 
          cleanup_diagnostic_variables


! Diagnostic variables

real, target, allocatable, dimension(:), public :: & 
  sigma_sqd_w_zt,     & ! PDF width parameter: t point          [-]
  Skw_zm,    & ! Skw on moment. grid                   [-]
  Skw_zt,    & ! Skw on thermo. grid                   [-]
  ug,        & ! u geostrophic wind                    [m/s]
  vg,        & ! v geostrophic wind                    [m/s]
  um_ref,    & ! Initial u wind; Michael Falk,         [m/s]
  vm_ref,    & ! Initial v wind; Michael Falk,         [m/s]
  thvm,      & ! Virtual potential Temperature         [K]
  shear        ! Wind shear production

!$omp   threadprivate(sigma_sqd_w_zt, Skw_zm, Skw_zt, ug, vg)
!$omp   threadprivate(thvm, shear)
!$omp   threadprivate(um_ref, vm_ref)

real, target, allocatable, dimension(:), public :: & 
  rsat ! Saturation mixing ratio  ! Brian
       !$omp   threadprivate(rsat)

real, target, allocatable, dimension(:), public :: & 
  Frad,  & ! Radiative flux (momentum point)
  radht    ! SW + LW heating rate

!$omp   threadprivate(Frad, radht)

! Second order moments
real, target, allocatable, dimension(:), public :: & 
  wprcp,    & ! w'rc'                [m kg/s kg]
  thlprcp,  & ! thl'rc'              [K kg/kg]
  rtprcp,   & ! rt'rc'               [kg^2/kg^2]
  rcp2        ! rc'^2                [kg^2/kg^2]

!$omp   threadprivate(wprcp, thlprcp, rtprcp, rcp2)

! Third order moments
real, target, allocatable, dimension(:), public :: & 
  wpthlp2,   & ! w'thl'^2    [m K^2/s]
  wp2thlp,   & ! w'^2 thl'   [m^2 K/s^2]
  wprtp2,    & ! w'rt'^2     [m kg^2/kg^2]
  wp2rtp,    & ! w'^2rt'     [m^2 kg/kg]
  wprtpthlp, & ! w'rt'thl'   [m kg K/kg s]
  wp2rcp       ! w'^2 rc'    [m^2 kg/kg s^2]

!$omp   threadprivate(wpthlp2, wp2thlp, wprtp2, wp2rtp )
!$omp   threadprivate(wprtpthlp, wp2rcp)

! Fourth order moments
real, target, allocatable, dimension(:), public :: & 
  wp4 ! w'^4      [m^4/s^4]

!$omp   threadprivate(wp4)

! Buoyancy related moments
real, target, allocatable, dimension(:), public :: & 
  wpthvp,   & ! w'thv'       [m K/s]
  rtpthvp,  & ! rt' thv'     [kg K/kg]
  thlpthvp, & ! thl'thv'     [K^2] 
  wp2thvp     ! w'^2 thv'    [m^2 K/s^2]

!$omp   threadprivate(wpthvp, rtpthvp, thlpthvp, wp2thvp)

real, target, allocatable, dimension(:), public :: & 
  Kh_zt,  & ! Eddy diffusivity: zt grid        [m^2/s]
  Kh_zm     ! Eddy diffusivity: zm grid        [m^2/s]

!$omp   threadprivate(Kh_zt, Kh_zm)

! Mixing lengths
real, target, allocatable, dimension(:), public :: & 
  Lscale, Lscale_up, Lscale_down ! [m]

!$omp   threadprivate(Lscale, Lscale_up, Lscale_down)

real, target, allocatable, dimension(:), public :: & 
  em,   & ! em               [m^2/s^2]
  tau_zt  ! Dissipation time [s]

!$omp   threadprivate(em, tau_zt)

! hydrometeors variable array
real, allocatable, dimension(:,:), public :: hydromet
! When running with COAMPS microphysics this contains:
! 1 rrainm      Rain water mixing ratio               [kg/kg]
! 2 Nrm      Rain drop number concentration        [num/kg]
! 3 rsnow    Snow water mixing ratio               [kg/kg]
! 4 rice     Ice water mixing ratio                [kg/kg]
! 5 rgraupel Graupel water mixing ratio            [kg/kg]
!$omp   threadprivate(hydromet)

real, target, allocatable, dimension(:), public :: & 
  Ncm,   & ! Cloud droplet number concentration      [num/kg]
  Ncnm,  & ! Cloud nuclei number concentration       [num/m^3]
  Nim      ! Ice nuclei number concentration         [num/m^3]
!$omp   threadprivate(Ncm, Ncnm, Nim)


! Surface data
real, public  :: ustar ! Average value of friction velocity [m/s]

!$omp   threadprivate(ustar)

! Passive scalar variables

real, target, allocatable, dimension(:,:), public :: & 
  wpedsclrp   ! w'edsclr'

real, target, allocatable, dimension(:,:), public :: & 
  sclrpthvp,   & ! sclr'th_v'
  sclrprtp,    & ! sclr'rt'
  sclrp2,      & ! sclr'^2
  sclrpthlp,   & ! sclr'th_l'
  sclrprcp,    & ! sclr'rc'
  wp2sclrp,    & ! w'^2 sclr'
  wpsclrp2,    & ! w'sclr'^2
  wpsclrprtp,  & ! w'sclr'rt'
  wpsclrpthlp    ! w'sclr'thl'

!$omp   threadprivate(wpedsclrp)
!$omp   threadprivate(sclrpthvp, sclrprtp, sclrp2, sclrpthlp)
!$omp   threadprivate(wp2sclrp, wpsclrprtp, wpsclrpthlp)

! Interpolated variables for tuning
! 
real, target, allocatable, dimension(:), public :: & 
  wp2_zt,     & ! w'^2 on thermo. grid
  thlp2_zt,   & ! thl'^2 on thermo. grid
  wpthlp_zt,  & ! w'thl' on thermo. grid
  wprtp_zt,   & ! w'rt' on thermo. grid
  rtp2_zt,    & ! rt'^2 on therm. grid
  rtpthlp_zt    ! rt'thl' on thermo. grid

!$omp   threadprivate(wp2_zt, thlp2_zt, wpthlp_zt, wprtp_zt) 
!$omp   threadprivate(rtp2_zt, rtpthlp_zt)

! 

! Variables needed for the pdf closure scheme
!
!       pdf_parms contains the parameters of the pdf:
!
!        pdf_parms(:,1) = w1
!        pdf_parms(:,2) = w2
!        pdf_parms(:,3) = sw1
!        pdf_parms(:,4) = sw2
!        pdf_parms(:,5) = rt1
!        pdf_parms(:,6) = rt2
!        pdf_parms(:,7) = srt1
!        pdf_parms(:,8) = srt2
!        pdf_parms(:,9) = thl1
!        pdf_parms(:,10) = thl2
!        pdf_parms(:,11) = sthl1
!        pdf_parms(:,12) = sthl2
!        pdf_parms(:,13) = a
!        pdf_parms(:,14) = rc1
!        pdf_parms(:,15) = rc2
!        pdf_parms(:,16) = rsl1
!        pdf_parms(:,17) = rsl2
!        pdf_parms(:,18) = R1
!        pdf_parms(:,19) = R2
!        pdf_parms(:,20) = s1
!        pdf_parms(:,21) = s2
!        pdf_parms(:,22) = ss1
!        pdf_parms(:,23) = ss2
!        pdf_parms(:,24) = rrtthl

real, target, allocatable, dimension(:,:), public :: & 
  pdf_parms

!$omp   threadprivate(pdf_parms)

! Latin Hypercube arrays.  Vince Larson 22 May 2005
real, target, allocatable, dimension(:), public :: & 
  AKm_est,   & ! Kessler ac estimate                 [kg/kg]
  AKm,       & ! Exact Kessler ac                    [kg/kg]
  AKstd,     & ! St dev of exact Kessler ac          [???]
  AKstd_cld, & ! Stdev of exact w/in cloud ac        [???]
  rcm_est,   & ! Monte Carlo rcm estimate            [kg/kg]
  AKm_rcm,   & ! Kessler ac based on rcm             [???]
  AKm_rcc      ! Kessler ac based on rcm/cf          [???]

!$omp   threadprivate(AKm_est, AKm, AKstd, AKstd_cld, rcm_est, AKm_rcm)
!$omp   threadprivate(AKm_rcc)

contains

!-----------------------------------------------------------------------
!  Allocates and initializes prognostic scalar and array variables 
!  for the HOC model code
!-----------------------------------------------------------------------
subroutine setup_diagnostic_variables( nzmax )
use model_flags, only:  & 
    l_LH_on ! Variable(s)

use constants, only:  & 
    emin ! Variables

use parameters, only: & 
    hydromet_dim, & 
    sclr_dim

implicit none

! Constant Parameters 
integer, parameter :: pdf_dimension = 26

! Input Variables
integer, intent(in) :: nzmax

! Local Variables
integer :: i

!$omp   parallel
!   --- Allocation ---

! Diagnostic variables

allocate( sigma_sqd_w_zt(1:nzmax) )  ! PDF width parameter: t point
allocate( Skw_zm(1:nzmax) )          ! Skw
allocate( Skw_zt(1:nzmax) )          ! Skw
allocate( ug(1:nzmax) )              ! u geostrophic wind
allocate( vg(1:nzmax) )              ! v geostrophic wind
allocate( um_ref(1:nzmax) )          ! Reference u wind for nudging; Michael Falk, 17 Oct 2007
allocate( vm_ref(1:nzmax) )          ! Reference v wind for nudging; Michael Falk, 17 Oct 2007
allocate( thvm(1:nzmax) )            ! Virtual potential temperature

allocate( rsat(1:nzmax) )       ! Saturation mixing ratio  ! Brian

allocate( Frad(1:nzmax) )      ! radiative flux (momentum point)

allocate( radht(1:nzmax) )     ! SW + LW heating rate

allocate( shear(1:nzmax) )     ! wind shear production

! Second order moments

allocate( wprcp(1:nzmax) )     ! w'rc'
allocate( thlprcp(1:nzmax) )   ! thl'rc'
allocate( rtprcp(1:nzmax) )    ! rt'rc'
allocate( rcp2(1:nzmax) )      ! rc'^2

! Third order moments

allocate( wpthlp2(1:nzmax) )   ! w'thl'^2
allocate( wp2thlp(1:nzmax) )   ! w'^2thl'
allocate( wprtp2(1:nzmax) )    ! w'rt'^2
allocate( wp2rtp(1:nzmax) )    ! w'^2rt'
allocate( wprtpthlp(1:nzmax) ) ! w'rt'thl'
allocate( wp2rcp(1:nzmax) )    ! w'^2rc'

! Fourth order moments

allocate( wp4(1:nzmax) )

! Buoyancy related moments

allocate( wpthvp(1:nzmax) )
allocate( rtpthvp(1:nzmax) )
allocate( thlpthvp(1:nzmax) )
allocate( wp2thvp(1:nzmax) )

allocate( Kh_zt(1:nzmax) )
allocate( Kh_zm(1:nzmax) )

allocate( em(1:nzmax) )
allocate( Lscale(1:nzmax) )
allocate( Lscale_up(1:nzmax) )
allocate( Lscale_down(1:nzmax) )
allocate( tau_zt(1:nzmax) )
!       allocate( tau_zm(1:nzmax) )

 
! Tuning Variables
allocate( wp2_zt(1:nzmax) )     ! w'^2 on thermo. grid
allocate( thlp2_zt(1:nzmax) )   ! thl'^2 on thermo. grid
allocate( wpthlp_zt(1:nzmax) )  ! w'thl' on thermo. grid
allocate( wprtp_zt(1:nzmax) )   ! w'rt' on thermo. grid
allocate( rtp2_zt(1:nzmax) )    ! rt'^2 on thermo. grid
allocate( rtpthlp_zt(1:nzmax) ) ! rt'thl' on thermo. grid
 

! Array fpr pdf closure scheme

allocate( pdf_parms(1:nzmax,1:pdf_dimension) )

allocate( Ncm(1:nzmax) )
allocate( Ncnm(1:nzmax) )
allocate( Nim(1:nzmax) )
allocate( hydromet(1:nzmax,1:hydromet_dim) ) ! All hydrometeor fields


! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
!       if ( l_LH_on ) then
  allocate( AKm_est(1:nzmax) )    ! Kessler ac estimate
  allocate( AKm(1:nzmax) )        ! Exact Kessler ac
  allocate( AKstd(1:nzmax) )      ! St dev of exact Kessler ac
  allocate( AKstd_cld(1:nzmax) )  ! St dev of exact w/in cloud Kessler ac
  allocate( rcm_est(1:nzmax) )      ! Monte Carlo rcm estimate
  allocate( AKm_rcm(1:nzmax) )      ! Kessler ac based on rcm
  allocate( AKm_rcc(1:nzmax) )      ! Kessler ac based on rcm/cf
!       end if ! l_LH_on
! End of variables for Latin hypercube.

! Variables for new mixing scheme
allocate( sclrprtp(1:nzmax, 1:sclr_dim) )
allocate( sclrp2(1:nzmax, 1:sclr_dim) )
allocate( sclrpthvp(1:nzmax, 1:sclr_dim) )
allocate( sclrpthlp(1:nzmax, 1:sclr_dim) )
allocate( sclrprcp(1:nzmax, 1:sclr_dim) )

allocate( wp2sclrp(1:nzmax, 1:sclr_dim) )
allocate( wpsclrp2(1:nzmax, 1:sclr_dim) )
allocate( wpsclrprtp(1:nzmax, 1:sclr_dim) )
allocate( wpsclrpthlp(1:nzmax, 1:sclr_dim) )

! Eddy Diff. Scalars
allocate( wpedsclrp(1:nzmax, 1:sclr_dim) )

!   --- Initializaton ---

! Diagnostic variables

sigma_sqd_w_zt = 0.0      ! PDF width parameter: t point

Skw_zm = 0.0
Skw_zt = 0.0

ug  = 0.0      ! u geostrophic wind
vg  = 0.0      ! v geostrophic wind
um_ref   = 0.0 !
vm_ref   = 0.0 !
 
thvm = 0.0  ! Virtual potential temperature
rsat  = 0.0  ! Saturation mixing ratio  ! Brian

radht = 0.0 ! Heating rate
Frad  = 0.0 ! Radiative flux

shear = 0.0    ! Wind shear production

! Second order moments
wprcp   = 0.0
thlprcp = 0.0
rtprcp  = 0.0
rcp2    = 0.0

! Third order moments
wpthlp2   = 0.0
wp2thlp   = 0.0
wprtp2    = 0.0
wp2rtp    = 0.0
wp2rcp    = 0.0
wprtpthlp = 0.0

! Fourth order moments
wp4 = 0.0

! Buoyancy related moments
wpthvp   = 0.0
rtpthvp  = 0.0
thlpthvp = 0.0
wp2thvp  = 0.0

! Eddy diffusivity
Kh_zt      = 0.0
Kh_zm      = 0.0

em       = emin

! Length scale
Lscale   = 0.0
Lscale_up      = 0.0
Lscale_down    = 0.0

! Dissipation time
tau_zt     = 0.0

! Hydrometer types
Ncm(1:nzmax)  = 0.0
Ncnm(1:nzmax) = 0.0
Nim(1:nzmax)  = 0.0

do i = 1, hydromet_dim, 1
  hydromet(1:nzmax,i) = 0.0
end do

! Array for pdf closure scheme
pdf_parms(1:nzmax,:) = 0.0

! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
if ( l_LH_on ) then
  AKm_est   = 0.0  ! Kessler ac estimate
  AKm       = 0.0  ! Exact Kessler ac
  AKstd     = 0.0  ! St dev of exact Kessler ac
  AKstd_cld = 0.0  ! St dev of exact w/in cloud Kessler ac
  rcm_est   = 0.0  ! Monte Carlo rcm estimate
  AKm_rcm   = 0.0  ! Kessler ac based on rcm
  AKm_rcc   = 0.0  ! Kessler ac based on rcm/cf
end if ! l_LH_on

! Passive scalars
if ( sclr_dim > 0 ) then
  sclrprtp(:,:)      = 0.0
  sclrp2(:,:)        = 0.0
  sclrpthvp(:,:)     = 0.0
  sclrpthlp(:,:)     = 0.0
  sclrprcp(:,:)      = 0.0

  wp2sclrp(:,:)      = 0.0
  wpsclrp2(:,:)      = 0.0
  wpsclrprtp(:,:)    = 0.0
  wpsclrpthlp(:,:)   = 0.0

  wpedsclrp(:,:)     = 0.0
end if
!$omp   end parallel

return
end subroutine setup_diagnostic_variables

!------------------------------------------------------------------------
subroutine cleanup_diagnostic_variables( )

!       Description:
!       Subroutine to deallocate variables defined in module global
!------------------------------------------------------------------------
use model_flags, only: & 
    l_LH_on ! Variable(s)

implicit none

!$omp   parallel

! --- Deallocate --- 

deallocate( sigma_sqd_w_zt )       ! PDF width parameter: t point
deallocate( Skw_zm )
deallocate( Skw_zt )
deallocate( ug )        ! u geostrophic wind
deallocate( vg )        ! v geostrophic wind
deallocate( um_ref )    ! u initial
deallocate( vm_ref )    ! v initial
 
deallocate( thvm )      ! virtual potential temperature
deallocate( rsat )      ! saturation mixing ratio  ! Brian

deallocate( Frad )      ! radiative flux (momentum point)

deallocate( radht )     ! SW + LW heating rate

deallocate( shear )     ! wind shear production

! Second order moments

deallocate( wprcp )     ! w'rc'
deallocate( thlprcp )   ! thl'rc'
deallocate( rtprcp )    ! rt'rc'
deallocate( rcp2 )      ! rc'^2

! Third order moments

deallocate( wpthlp2 )   ! w'thl'^2
deallocate( wp2thlp )   ! w'^2thl'
deallocate( wprtp2 )    ! w'rt'^2
deallocate( wp2rtp )    ! w'^2rt'
deallocate( wprtpthlp ) ! w'rt'thl'
deallocate( wp2rcp )    ! w'^2rc'

! Fourth order moments

deallocate( wp4 )

! Buoyancy related moments

deallocate( wpthvp )
deallocate( rtpthvp )
deallocate( thlpthvp )
deallocate( wp2thvp )

deallocate( Kh_zt )
deallocate( Kh_zm )

deallocate( em )
deallocate( Lscale )
deallocate( Lscale_up )
deallocate( Lscale_down )
deallocate( tau_zt )

! Cloud water variables

deallocate( Ncm )
deallocate( Ncnm )
deallocate( Nim )

deallocate( hydromet )  ! Hydrometeor fields

 
! Interpolated variables for tuning
deallocate( wp2_zt )     ! w'^2 on t
deallocate( thlp2_zt )   ! th_l'^2 on t 
deallocate( wpthlp_zt )  ! w'th_l' on t
deallocate( wprtp_zt )   ! w'rt' on t
deallocate( rtp2_zt )    ! rt'^2 on t
deallocate( rtpthlp_zt ) ! rt'th_l' on t
 


! Array for pdf closure scheme

deallocate( pdf_parms )

! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
!       if ( l_LH_on ) then
  deallocate( AKm_est )   ! Kessler ac estimate
  deallocate( AKm )       ! Exact Kessler ac
  deallocate( AKstd )     ! St dev of exact Kessler ac
  deallocate( AKstd_cld ) ! St dev of exact w/in cloud Kessler ac
  deallocate( rcm_est )   ! Monte Carlo rcm estimate
  deallocate( AKm_rcm )   ! Kessler ac based on rcm
  deallocate( AKm_rcc )   ! Kessler ac based on rcm/cf
!       end if ! l_LH_on

! Passive scalars
deallocate( sclrprtp )
deallocate( sclrp2 )
deallocate( sclrpthvp )
deallocate( sclrpthlp )
deallocate( sclrprcp )

deallocate( wp2sclrp )
deallocate( wpsclrp2 )
deallocate( wpsclrprtp )
deallocate( wpsclrpthlp )

deallocate( wpedsclrp )
!$omp   end parallel

return
end subroutine cleanup_diagnostic_variables

end module diagnostic_variables
