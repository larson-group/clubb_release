!-----------------------------------------------------------------------
! $Id: prog_variables.F,v 1.17 2008-03-25 04:58:22 dschanen Exp $
        module prognostic_variables

!       This module contains definitions of all prognostic
!       arrays used in the single column model, as well as subroutines
!       to allocate, deallocate and initialize them.

!       Note that while these are all same dimension, there is a
!       thermodynamic grid and a momentum grid, and the grids have 
!       different points.
!-----------------------------------------------------------------------
        use constants, only: sclr_dim

        implicit none

        ! Prognostic variables
        real, target, allocatable, dimension(:) ::
     .  um,     ! u wind                        [m/s]
     .  vm,     ! v wind                        [m/s]
     .  upwp,   ! vertical u momentum flux      [m^2/s^2]
     .  vpwp,   ! vertical v momentum flux      [m^2/s^2]
     .  up2,    ! u'^2                          [m^2/s^2]
     .  vp2,    ! v'^2                          [m^2/s^2]
     .  thlm,   ! liquid potential temperature  [K]
     .  rtm,    ! total water mixing ratio      [kg/kg]
     .  wprtp,  ! w'rt'                         [m kg/s kg]
     .  wpthlp, ! w'thl'                        [m K /s]
     .  wp2,    ! w'^2                          [m^2/s^2]
     .  wp3,    ! w'^3                          [m^3/s^3]
     .  rtp2,   ! rt'^2                         [kg/kg]
     .  thlp2,  ! thl'^2                        [K^2]
     .  rtpthlp ! rt'thl'                       [kg/kg K]

!$omp   threadprivate(um, vm, upwp, vpwp, up2, vp2)
!$omp   threadprivate(thlm, rtm, wprtp, wpthlp, wp2)
!$omp   threadprivate(wp3, rtp2, thlp2, rtpthlp)

        real, target, allocatable, dimension(:) ::
     .  p,           ! Pressure (Pa) on thermodynamic points    [Pa]
     .  exner,       ! Exner function = ( p / p0 ) ** kappa     [-]
     .  rhot,        ! Density                                  [kg/m^3]
     .  rhom,        ! Density                                  [kg/m^3]
     .  thlm_forcing,! thlm large-scale forcing                 [K/s]
     .  rtm_forcing  ! rtm large-scale forcing                  [kg/kg/s]

!$omp   threadprivate(p, exner, rhot, rhom, thlm_forcing, rtm_forcing)

        ! Imposed large scale w
        real, target, allocatable, dimension(:) ::
     .  wmm,! w on momentum levels              [m/s]
     .  wmt ! w on thermodynamic levels         [m/s]

        ! PDF width parameter: momentum levels
        real, target, allocatable, dimension(:) :: Scm  ! [-]

        ! Mixing Lengths
        real, target, allocatable, dimension(:) :: taum ! [s]

!$omp   threadprivate(wmm, wmt, Scm, taum)

        ! Cloud water variables
        real, target, allocatable, dimension(:) ::
     .  rcm,  ! Cloud water mixing ratio                [kg/kg]
     .  Ncm,  ! Cloud droplet number concentration      [num/kg]
     .  Ncnm, ! Cloud nuclei number concentration       [num/m^3]
     .  Nim,  ! Ice nuclei number concentration         [num/m^3]
     .  cf    ! Cloud fraction                          [%]

!$omp   threadprivate(rcm, Ncm, Ncnm, Nim, cf)

        ! Prognostic drizzle variable array
        real, target, allocatable, dimension(:,:) :: hydromet
        ! Contains:
        ! 1 rrm      Rain water mixing ratio               [kg/kg]
        ! 2 Nrm      Rain drop number concentration        [num/kg]
        ! 3 rsnow    Snow water mixing ratio               [kg/kg]
        ! 4 rice     Ice water mixing ratio                [kg/kg]
        ! 5 rgraupel Graupel water mixing ratio            [kg/kg]

!$omp   threadprivate(hydromet)

        ! Surface fluxes
        real :: 
     .  wpthlp_sfc,        ! w'thl'      [m K/s]
     .  wprtp_sfc,         ! w'rt'       [m kg/kg s]
     .  upwp_sfc, vpwp_sfc ! u'w' & v'w' [m^2/s^2]

!$omp   threadprivate(wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc)

#ifdef SCALARS 
        ! Surface fluxes for passive scalars
        real, dimension(sclr_dim) ::
     .  wpsclrp_sfc,    ! w'sclr' at surface    [units m/s]
     .  wpedsclrp_sfc   ! w'edsclr' at surface  [units m/s]

!$omp   threadprivate(wpsclrp_sfc, wpedsclrp_sfc)
#endif /* SCALARS */

        ! More surface data
        real :: 
     .  Tsfc, ! surface temperature     [K]
     .  psfc, ! surface pressure        [Pa]
     .  SE,   ! sensible heat flux      [K/s]
     .  LE    ! latent heat flux        [1/s]

!$omp   threadprivate(Tsfc, psfc, SE, LE)

#ifdef SCALARS

        ! Passive scalars 
        real, target, allocatable, dimension(:,:) ::
     .  sclrm,         ! Mean passive scalars           [units vary]
     .  sclrm_forcing, ! Scalars' forcing               [units/s]
     .  edsclrm        ! Mean eddy-diffusivity scalars  [units vary]

!$omp   threadprivate(sclrm, sclrm_forcing, edsclrm)

#endif /*SCALARS*/

        contains
!-----------------------------------------------------------------------
        subroutine setup_prognostic_variables( nzmax )

!       Description:
!       Allocates and Initializes prognostic scalar and array variables 
!       for the HOC model code

!       References:
!       None
!-----------------------------------------------------------------------
        use constants, only: sclr_dim, emin, hydromet_dim
        use model_flags

        implicit none

        integer, intent(in) :: nzmax
!$omp   parallel
!   --- Allocation ---

! Prognostic variables

        allocate( um(1:nzmax) )        ! u wind
        allocate( vm(1:nzmax) )        ! v wind

        allocate( upwp(1:nzmax) )      ! vertical u momentum flux
        allocate( vpwp(1:nzmax) )      ! vertical v momentum flux

        allocate( up2(1:nzmax) )
        allocate( vp2(1:nzmax) )

        allocate( thlm(1:nzmax) )      ! liquid potential temperature
        allocate( rtm(1:nzmax) )       ! total water mixing ratio
        allocate( wprtp(1:nzmax) )     ! w'rt'
        allocate( wpthlp(1:nzmax) )    ! w'thl'
        allocate( wp2(1:nzmax) )       ! w'^2
        allocate( wp3(1:nzmax) )       ! w'^3
        allocate( rtp2(1:nzmax) )      ! rt'^2
        allocate( thlp2(1:nzmax) )     ! thl'^2
        allocate( rtpthlp(1:nzmax) )   ! rt'thlp'

        allocate( p(1:nzmax) )         ! pressure (pascals)
        allocate( exner(1:nzmax) )     ! exner
        allocate( rhot(1:nzmax) )      ! density: t points
        allocate( rhom(1:nzmax) )      ! density: m points

        allocate( thlm_forcing(1:nzmax) ) ! thlm ls forcing
        allocate( rtm_forcing(1:nzmax) )  ! rtm ls forcing

! Imposed large scale w

        allocate( wmm(1:nzmax) )       ! momentum levels
        allocate( wmt(1:nzmax) )       ! thermodynamic levels

! PDF width parameter: momentum levels

        allocate( Scm(1:nzmax) ) 

! Mixing lengths

        allocate( taum(1:nzmax) ) 

        ! Cloud water variables

        allocate( rcm(1:nzmax) )
        allocate( Ncm(1:nzmax) )
        allocate( Ncnm(1:nzmax) )
        allocate( Nim(1:nzmax) )
        allocate( cf(1:nzmax) )

        allocate( hydromet(1:nzmax,1:hydromet_dim) ) ! All hydrometeor fields

#ifdef SCALARS

        ! Passive scalar variables 

        allocate( sclrm(1:nzmax, 1:sclr_dim) )
        allocate( sclrm_forcing(1:nzmax, 1:sclr_dim) )

        allocate( edsclrm(1:nzmax, 1:sclr_dim) )

#endif /*SCALARS*/

!--------- Set initial values for array variables ---------

        ! Prognostic variables

        um(1:nzmax)      = 0.0     ! u wind
        vm (1:nzmax)     = 0.0     ! v wind

        upwp(1:nzmax)    = 0.0     ! vertical u momentum flux
        vpwp(1:nzmax)    = 0.0     ! vertical v momentum flux

        up2(1:nzmax)     = 2./3. * emin ! u'^2
        vp2(1:nzmax)     = 2./3. * emin ! v'^2
        wp2(1:nzmax)     = 2./3. * emin ! w'^2

        thlm(1:nzmax)    = 0.0     ! liquid potential temperature
        rtm(1:nzmax)     = 0.0     ! total water mixing ratio
        wprtp(1:nzmax)   = 0.0     ! w'rt'
        wpthlp(1:nzmax)  = 0.0     ! w'thl'
        wp3(1:nzmax)     = 0.0     ! w'^3
        rtp2(1:nzmax)    = 0.0     ! rt'^2
        thlp2(1:nzmax)   = 0.0     ! thl'^2
        rtpthlp(1:nzmax) = 0.0     ! rt'thl'
 
        p(1:nzmax)     = 0.0    ! pressure (Pa)
        exner(1:nzmax) = 0.0    ! exner
        rhot(1:nzmax)  = 0.0    ! density on thermo. levels
        rhom(1:nzmax)  = 0.0    ! density on moment. levels

        thlm_forcing(1:nzmax) = 0.0     ! thlm large-scale forcing
        rtm_forcing(1:nzmax)  = 0.0     ! rtm large-scale forcing

        ! Imposed large scale w

        wmm(1:nzmax) = 0.0      ! Momentum levels
        wmt(1:nzmax) = 0.0      ! Thermodynamic levels

        ! PDF width parameter: momentum levels

        Scm(1:nzmax)  = 0.0 

        ! Mixing lengths

        taum(1:nzmax) = 0.0

        ! Cloud water variables

        rcm(1:nzmax)  = 0.0
        Ncm(1:nzmax)  = 0.0
        Ncnm(1:nzmax) = 0.0
        Nim(1:nzmax)  = 0.0
        cf(1:nzmax)   = 0.0

        ! Hydrometer types
        hydromet(1:nzmax,1:hydromet_dim) = 0.0

        ! Surface fluxes
        wpthlp_sfc = 0.0
        wprtp_sfc  = 0.0
        upwp_sfc   = 0.0
        vpwp_sfc   = 0.0

#ifdef SCALARS
        ! Passive scalars 
        wpsclrp_sfc(1:sclr_dim)   = 0.0
        wpedsclrp_sfc(1:sclr_dim) = 0.0

        sclrm(1:nzmax,1:sclr_dim)         = 0.0
        sclrm_forcing(1:nzmax,1:sclr_dim) = 0.0

        edsclrm(1:nzmax,1:sclr_dim)       = 0.0
#endif /*SCALARS*/
!$omp   end parallel

        return
        end subroutine setup_prognostic_variables
!-----------------------------------------------------------------------
        subroutine cleanup_prognostic_variables
        implicit none

        ! Prognostic variables
!$omp   parallel

        deallocate( um )        ! u wind
        deallocate( vm )        ! v wind

        deallocate( upwp )      ! vertical u momentum flux
        deallocate( vpwp )      ! vertical v momentum flux

        deallocate( up2, vp2 )

        deallocate( thlm )      ! liquid potential temperature
        deallocate( rtm )       ! total water mixing ratio
        deallocate( wprtp )     ! w'rt'
        deallocate( wpthlp )    ! w'thl'
        deallocate( wp2 )       ! w'^2
        deallocate( wp3 )       ! w'^3
        deallocate( rtp2 )      ! rt'^2
        deallocate( thlp2 )     ! thl'^2
        deallocate( rtpthlp )   ! rt'thl'

        deallocate( p )         ! pressure
        deallocate( exner )     ! exner
        deallocate( rhot )      ! density: t points
        deallocate( rhom )      ! density: m points

        deallocate( thlm_forcing )
        deallocate( rtm_forcing )

        ! Imposed large scale w

        deallocate( wmm )       ! momentum levels
        deallocate( wmt )       ! thermodynamic levels

        ! PDF width parameter

        deallocate( Scm )

        ! Mixing lengths

        deallocate( taum )

        ! Cloud water variables

        deallocate( rcm )
        deallocate( Ncm )
        deallocate( Ncnm )
        deallocate( Nim )
        deallocate( cf )

        ! Prognostic drizzle variables

        deallocate( hydromet )  ! Hydrometeor fields

#ifdef SCALARS
        ! Passive scalars

        deallocate( sclrm )
        deallocate( sclrm_forcing )

        deallocate( edsclrm )

#endif /*SCALARS*/
!$omp   end parallel

        return
        end subroutine cleanup_prognostic_variables

        end module prognostic_variables
