!-----------------------------------------------------------------------
! $Id: prog_variables.F,v 1.9 2006-08-04 21:08:56 dschanen Exp $
        module prognostic_variables

! This module contains definitions of all prognostic
! arrays used in the single column model, as well as subroutines
! to allocate, deallocate and initialize them.

! Note that while these are all same dimension, there is a
! thermodynamic and momentum grid and they have different levels
!-----------------------------------------------------------------------
        implicit none

! Prognostic variables

        real, target, allocatable :: um(:)      ! u wind
        real, target, allocatable :: vm(:)      ! v wind

        real, allocatable, target :: upwp(:)      ! vertical u momentum flux
        real, allocatable, target :: vpwp(:)      ! vertical v momentum flux

        real, target, allocatable :: thlm(:)    ! liquid potential temperature
        real, target, allocatable :: rtm(:)     ! total water mixing ratio
        real, target, allocatable :: wprtp(:)   ! w'rt'
        real, target, allocatable :: wpthlp(:)  ! w'thl'
        real, target, allocatable :: wp2(:)     ! w'^2
        real, target, allocatable :: wp3(:)     ! w'^3
        real, target, allocatable :: rtp2(:)    ! rt'^2
        real, target, allocatable :: thlp2(:)   ! thl'^2
        real, target, allocatable :: rtpthlp(:) ! rt'thl'

        real, target, allocatable :: p(:)
                         ! pressure (Pa) on thermodynamic points
        real, target, allocatable :: exner(:)        
                         ! exner = ( p / p0 ) ** kappa

        real, target, allocatable :: rhot(:) ! density (kg/m^3)
        real, target, allocatable :: rhom(:) ! density (kg/m^3)

        real, target, allocatable :: thlm_forcing(:) ! thlm ls forcing
        real, target, allocatable :: rtm_forcing(:)  ! rtm ls forcing

! Imposed large scale w

        real, target, allocatable :: wmm(:) ! momentum levels
        real, target, allocatable :: wmt(:) ! thermodynamic levels

! PDF width parameter: momentum levels

        real, target, allocatable :: Scm(:)   

! Mixing Lengths

        real, target, allocatable :: taum(:)

! Cloud water variables

        real, allocatable, target :: rcm(:)
        real, allocatable, target :: Ncm(:)
        real, allocatable, target :: cf(:)

! Prognostic drizzle variables

        real, allocatable, target :: rrm(:)  ! rain water mixing ratio  ! Brian
        real, allocatable, target :: Nrm(:)  ! rain droplet number conc.! Brian

! Surface fluxes

        real :: wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc

! More surface data

        real :: Tsfc             ! surface temperature (K)
        real :: psfc             ! surface pressure (Pa)
        real :: SE               ! sensible heat flux (K/s)
        real :: LE               ! latent heat flux (1/s)

#ifdef SCALARS
! New mixing scheme variables

        real, target, allocatable :: sclrm(:,:)          ! passive scalars
        real, target, allocatable :: sclrm_forcing(:,:)  ! scalars' forcing

        real, target, allocatable :: edsclrm(:,:)  ! eddy-diff. scalar mean
        real, target, allocatable :: edsclrmt(:,:) ! edsclrm tendency

#endif /*SCALARS*/

        contains
!-----------------------------------------------------------------------
!  Allocates and Initializes prognostic scalar and array variables 
!  for the HOC model code
!-----------------------------------------------------------------------
        subroutine setup_prognostic_variables( nzmax )
        use constants, only: sclr_dim, emin
        use model_flags

        implicit none

        integer, intent(in) :: nzmax

!   --- Allocation ---

! Prognostic variables

        allocate( um(1:nzmax) )        ! u wind
        allocate( vm(1:nzmax) )        ! v wind

        allocate( upwp(1:nzmax) )      ! vertical u momentum flux
        allocate( vpwp(1:nzmax) )      ! vertical v momentum flux

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
        allocate( cf(1:nzmax) )

! Prognostic drizzle variables

        allocate( rrm(1:nzmax) )       ! rain water mixing ratio  ! Brian
        allocate( Nrm(1:nzmax) )       ! rain droplet number conc.  ! Brian

#ifdef SCALARS

! Variables for new mixing scheme

        allocate( sclrm(1:nzmax, 1:sclr_dim) )
        allocate( sclrm_forcing(1:nzmax, 1:sclr_dim) )

        allocate( edsclrm(1:nzmax, 1:sclr_dim) )
        allocate( edsclrmt(1:nzmax, 1:sclr_dim) )


#endif /*SCALARS*/

!   --- Initializaton ---

! Prognostic variables

        um           = 0.0     ! u wind
        vm           = 0.0     ! v wind

        upwp         = 0.0     ! vertical u momentum flux
        vpwp         = 0.0     ! vertical v momentum flux

        thlm         = 0.0     ! liquid potential temperature
        rtm          = 0.0     ! total water mixing ratio
        wprtp        = 0.0     ! w'rt'
        wpthlp       = 0.0     ! w'thl'
        wp2          = 2./3. * emin  ! w'^2
        wp3          = 0.0     ! w'^3
        rtp2         = 0.0     ! rt'^2
        thlp2        = 0.0     ! thl'^2
        rtpthlp      = 0.0     ! rt'thl'
 
        p     = 0.0    ! pressure (Pa)
        exner = 0.0    ! exner
        rhot  = 0.0    ! density on thermo. levels
        rhom  = 0.0    ! density on moment. levels

        thlm_forcing = 0.0     ! thlm large scale forcing
        rtm_forcing  = 0.0     ! rtm large scale forcing

! Imposed large scale w

        wmm = 0.0      ! momentum levels
        wmt = 0.0      ! thermodynamic levels

! PDF width parameter: momentum levels

        Scm  = 0.0 

! Mixing lengths

        taum = 0.0

! Cloud water variables

        rcm = 0.0
        Ncm = 0.0
        cf  = 0.0

! Prognostic drizzle variables
        
        ! These must always be allocated, but not necessarily initialized
        if ( kk_rain ) then
          rrm          = 0.0  ! rain water mixing ratio  ! Brian
          Nrm          = 0.0  ! rain droplet number conc.! Brian
        end if ! kk_rain

#ifdef SCALARS
! New mixing scheme Variables
        sclrm(:,:)         = 0.0
        sclrm_forcing(:,:) = 0.0

        edsclrm(:,:)       = 0.0
        edsclrmt(:,:)      = 0.0
#endif /*SCALARS*/


        return
        end subroutine setup_prognostic_variables
!-----------------------------------------------------------------------
        subroutine cleanup_prognostic_variables
        implicit none

! Prognostic variables

        deallocate( um )        ! u wind
        deallocate( vm )        ! v wind

        deallocate( upwp )      ! vertical u momentum flux
        deallocate( vpwp )      ! vertical v momentum flux

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
        deallocate( cf )

! Prognostic drizzle variables

        deallocate( Nrm )       ! rain droplet number conc.! Brian
        deallocate( rrm )       ! rain water mixing ratio  ! Brian


#ifdef SCALARS
! New mixing scheme structures

        deallocate( sclrm )
        deallocate( sclrm_forcing )

        deallocate( edsclrm )
        deallocate( edsclrmt )

#endif /*SCALARS*/

        return
        end subroutine cleanup_prognostic_variables

        end module prognostic_variables
