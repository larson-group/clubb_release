!-----------------------------------------------------------------------
! $Id: prog_variables.F,v 1.3 2006-04-28 20:57:04 dschanen Exp $
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
        real, target, allocatable :: rtp2(:)    ! rt'^2
        real, target, allocatable :: thlp2(:)   ! thl'^2
        real, target, allocatable :: rtpthlp(:) ! rt'thl'
        real, target, allocatable :: wprtp(:)   ! w'rt'
        real, target, allocatable :: wpthlp(:)  ! w'thl'
        real, target, allocatable :: wp2(:)     ! w'^2
        real, target, allocatable :: wp3(:)     ! w'^3

        real, target, allocatable :: p(:)
                         ! pressure (Pa) on thermodynamic points
        real, target, allocatable :: exner(:)        
                         ! exner = ( p / p0 ) ** kappa

        real, target, allocatable :: rhot(:) ! density (kg/m^3)
        real, target, allocatable :: rhom(:) ! density (kg/m^3)

! Tendency arrays for prognostic variables

        real, allocatable :: umt(:)             ! u wind
        real, allocatable :: vmt(:)             ! v wind

        real, target, allocatable :: thlm_forcing(:) ! thlm ls forcing
        real, target, allocatable :: rtm_forcing(:)  ! rtm ls forcing

! PDF width parameter: momentum levels

        real, target, allocatable :: Scm(:)   

! Mixing Lengths

        real, target, allocatable :: taum(:)

! Imposed large scale w: m point

        real, target, allocatable :: wmm(:)    

! Prognostic drizzle variables

        real, allocatable, target :: rrm(:)  ! rain water mixing ratio  ! Brian
        real, allocatable, target :: Nrm(:)  ! rain droplet number conc.! Brian

! Surface data

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
        use constants, only: sclrm_dimension, emin
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
        allocate( wp2(1:nzmax) )       ! w'^2
        allocate( rtp2(1:nzmax) )      ! rt'^2
        allocate( thlp2(1:nzmax) )     ! thl'^2
        allocate( wprtp(1:nzmax) )     ! w'rt'
        allocate( wpthlp(1:nzmax) )    ! w'thl'
        allocate( rtpthlp(1:nzmax) )   ! rt'thlp'
        allocate( wp3(1:nzmax) )       ! w'^3

        allocate( p(1:nzmax) )         ! pressure (pascals)
        allocate( exner(1:nzmax) )     ! exner
        allocate( rhot(1:nzmax) )      ! density: t points
        allocate( rhom(1:nzmax) )      ! density: m points

! Tendency arrays for prognostic variables

        allocate( umt(1:nzmax) )       ! u wind
        allocate( vmt(1:nzmax) )       ! v wind

        allocate( thlm_forcing(1:nzmax) ) ! thlm ls forcing
        allocate( rtm_forcing(1:nzmax) )  ! rtm ls forcing

! PDF width parameter: momentum levels

        allocate( Scm(1:nzmax) ) 

! Mixing lengths

        allocate( taum(1:nzmax) ) 

! Imposed large scale w: m point

        allocate( wmm(1:nzmax) )       

! Prognostic drizzle variables

        allocate( rrm(1:nzmax) )       ! rain water mixing ratio  ! Brian
        allocate( Nrm(1:nzmax) )       ! rain droplet number conc.  ! Brian

#ifdef SCALARS

! Variables for new mixing scheme

        allocate( sclrm(1:nzmax, 1:sclrm_dimension) )
        allocate( sclrm_forcing(1:nzmax, 1:sclrm_dimension) )

        allocate( edsclrm(1:nzmax, 1:sclrm_dimension) )
        allocate( edsclrmt(1:nzmax, 1:sclrm_dimension) )


#endif /*SCALARS*/

!   --- Initializaton ---

! Prognostic variables

        um           = 0.0     ! u wind
        vm           = 0.0     ! v wind
        upwp         = 0.0     ! vertical u momentum flux
        vpwp         = 0.0     ! vertical v momentum flux
        thlm_forcing = 0.0     ! thlm large scale forcing
        rtm_forcing  = 0.0     ! rtm large scale forcing
        thlm         = 0.0     ! liquid potential temperature
        rtm          = 0.0     ! total water mixing ratio
        wp2          = 2./3. * emin  ! w'^2
        rtp2         = 0.0     ! rt'^2
        thlp2        = 0.0     ! thl'^2
        wprtp        = 0.0     ! w'rt'
        wpthlp       = 0.0     ! w'thl'
        rtpthlp      = 0.0     ! rt'thl'
        wp3          = 0.0     ! w'^3
 
        p     = 0.0    ! pressure (Pa)
        exner = 0.0    ! exner
        rhot  = 0.0    ! density on thermo. levels
        rhom  = 0.0    ! density on moment. levels



! Tendency arrays for prognostic variables

        umt = 0.0  ! u wind
        vmt = 0.0  ! v wind

! PDF width parameter: momentum levels

        Scm  = 0.0 

! Mixing lengths

        taum = 0.0

! Imposed large scale w: m point

        wmm = 0.0      
 
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

        deallocate( thlm_forcing )
        deallocate( rtm_forcing )
        deallocate( thlm )      ! liquid potential temperature
        deallocate( rtm )       ! total water mixing ratio
        deallocate( wp2 )       ! w'^2
        deallocate( rtp2 )      ! rt'^2
        deallocate( thlp2 )     ! thl'^2
        deallocate( wprtp )     ! w'rt'
        deallocate( wpthlp )    ! w'thl'
        deallocate( rtpthlp )   ! rt'thl'
        deallocate( wp3 )       ! w'^3

        deallocate( p )         ! pressure
        deallocate( exner )     ! exner
        deallocate( rhot )      ! density: t points
        deallocate( rhom )      ! density: m points

! PDF width parameter

        deallocate( Scm )

! Mixing lengths

        deallocate( taum )

! Imposed large scale w: m point

        deallocate( wmm )       

! Tendency arrays for prognostic variables

        deallocate( umt )       ! u wind
        deallocate( vmt )       ! v wind

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
