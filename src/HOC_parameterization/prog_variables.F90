!-----------------------------------------------------------------------
! $Id: prog_variables.F90,v 1.8 2008-08-06 21:38:59 faschinj Exp $
        module prognostic_variables

!       This module contains definitions of all prognostic
!       arrays used in the single column model, as well as subroutines
!       to allocate, deallocate and initialize them.

!       Note that while these are all same dimension, there is a
!       thermodynamic grid and a momentum grid, and the grids have 
!       different points.
!-----------------------------------------------------------------------

        implicit none

        private ! Set Default Scoping

        public ::  & 
        setup_prognostic_variables,  & 
        cleanup_prognostic_variables

        ! Prognostic variables
        real, target, allocatable, dimension(:), public :: & 
        um,      & ! u wind                        [m/s]
        vm,      & ! v wind                        [m/s]
        upwp,    & ! vertical u momentum flux      [m^2/s^2]
        vpwp,    & ! vertical v momentum flux      [m^2/s^2]
        up2,     & ! u'^2                          [m^2/s^2]
        vp2,     & ! v'^2                          [m^2/s^2]
        thlm,    & ! liquid potential temperature  [K]
        rtm,     & ! total water mixing ratio      [kg/kg]
        wprtp,   & ! w'rt'                         [m kg/s kg]
        wpthlp,  & ! w'thl'                        [m K /s]
        wp2,     & ! w'^2                          [m^2/s^2]
        wp3,     & ! w'^3                          [m^3/s^3]
        rtp2,    & ! rt'^2                         [kg/kg]
        thlp2,   & ! thl'^2                        [K^2]
        rtpthlp ! rt'thl'                       [kg/kg K]

!$omp   threadprivate(um, vm, upwp, vpwp, up2, vp2)
!$omp   threadprivate(thlm, rtm, wprtp, wpthlp, wp2)
!$omp   threadprivate(wp3, rtp2, thlp2, rtpthlp)

        real, target, allocatable, dimension(:), public :: & 
        p_in_Pa,      & ! Pressure (Pa) on thermodynamic points    [Pa]
        exner,        & ! Exner function = ( p / p0 ) ** kappa     [-]
        rho,         & ! Density                                  [kg/m^3]
        rho_zm,         & ! Density                                  [kg/m^3]
        thlm_forcing, & ! thlm large-scale forcing                 [K/s]
        rtm_forcing  ! rtm large-scale forcing                  [kg/kg/s]

!$omp   threadprivate(p, exner, rho, rho_zm, thlm_forcing, rtm_forcing)

        ! Imposed large scale w
        real, target, allocatable, dimension(:), public :: & 
        wm_zm, & ! w on momentum levels              [m/s]
        wm_zt ! w on thermodynamic levels         [m/s]

        ! PDF width parameter: momentum levels
        real, target, allocatable, dimension(:), public :: sigma_sqd_w  ! [-]

        ! Mixing Lengths
        real, target, allocatable, dimension(:), public :: tau_zm ! [s]

!$omp   threadprivate(wm_zm, wm_zt, sigma_sqd_w, tau_zm)

        ! Cloud water variables
        real, target, allocatable, dimension(:), public :: & 
        rcm,   & ! Cloud water mixing ratio                [kg/kg]
        cf    ! Cloud fraction                          [%]

!$omp   threadprivate(rcm, cf)

        ! Surface fluxes
        real, public ::  & 
        wpthlp_sfc,         & ! w'thl'      [m K/s]
        wprtp_sfc,          & ! w'rt'       [m kg/kg s]
        upwp_sfc, vpwp_sfc ! u'w' & v'w' [m^2/s^2]

!$omp   threadprivate(wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc)

        ! Surface fluxes for passive scalars
        real, dimension(:), allocatable, public :: & 
        wpsclrp_sfc,     & ! w'sclr' at surface    [units m/s]
        wpedsclrp_sfc   ! w'edsclr' at surface  [units m/s]

!$omp   threadprivate(wpsclrp_sfc, wpedsclrp_sfc)

        ! More surface data
        real, public ::  & 
        Tsfc,  & ! surface temperature     [K]
        psfc,  & ! surface pressure        [Pa]
        SE,    & ! sensible heat flux      [K/s]
        LE    ! latent heat flux        [1/s]

!$omp   threadprivate(Tsfc, psfc, SE, LE)

        ! Passive scalars 
        real, target, allocatable, dimension(:,:), public :: & 
        sclrm,          & ! Mean passive scalars           [units vary]
        sclrm_forcing,  & ! Scalars' forcing               [units/s]
        edsclrm,        & ! Mean eddy-diffusivity scalars  [units vary]
        wpsclrp        ! w'sclr'                        [units vary m/s]

!$omp   threadprivate(sclrm, sclrm_forcing, edsclrm, wpsclrp)

        contains
!-----------------------------------------------------------------------
        subroutine setup_prognostic_variables( nzmax )

!       Description:
!       Allocates and Initializes prognostic scalar and array variables 
!       for the HOC model code

!       References:
!       None
!-----------------------------------------------------------------------
        use constants, only:  & 
            emin, &
            rttol, &
            thltol

        use parameters, only: & 
            sclr_dim ! Variable(s) 

        implicit none

        integer, intent(in) :: nzmax

        integer :: i
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

        allocate( p_in_Pa(1:nzmax) )         ! pressure (pascals)
        allocate( exner(1:nzmax) )     ! exner function
        allocate( rho(1:nzmax) )      ! density: t points
        allocate( rho_zm(1:nzmax) )      ! density: m points

        allocate( thlm_forcing(1:nzmax) ) ! thlm ls forcing
        allocate( rtm_forcing(1:nzmax) )  ! rtm ls forcing

        ! Imposed large scale w

        allocate( wm_zm(1:nzmax) )       ! momentum levels
        allocate( wm_zt(1:nzmax) )       ! thermodynamic levels

        ! PDF width parameter: momentum levels

        allocate( sigma_sqd_w(1:nzmax) ) 

        ! Mixing lengths

        allocate( tau_zm(1:nzmax) ) 

        ! Cloud water variables

        allocate( rcm(1:nzmax) )
        allocate( cf(1:nzmax) )

        ! Passive scalar variables 
        ! Note that sclr_dim can be 0
        allocate( wpsclrp_sfc(1:sclr_dim), wpedsclrp_sfc(1:sclr_dim) )
        allocate( sclrm(1:nzmax, 1:sclr_dim) )
        allocate( sclrm_forcing(1:nzmax, 1:sclr_dim) )

        allocate( edsclrm(1:nzmax, 1:sclr_dim) )
        allocate( wpsclrp(1:nzmax, 1:sclr_dim) )


!--------- Set initial values for array variables ---------

        ! Prognostic variables

        um(1:nzmax)      = 0.0     ! u wind
        vm (1:nzmax)     = 0.0     ! v wind

        upwp(1:nzmax)    = 0.0     ! vertical u momentum flux
        vpwp(1:nzmax)    = 0.0     ! vertical v momentum flux

        up2(1:nzmax)     = 2./3. * emin ! u'^2
        vp2(1:nzmax)     = 2./3. * emin ! v'^2
        wp2(1:nzmax)     = 2./3. * emin ! w'^2

        thlm(1:nzmax)    = 0.0         ! liquid potential temperature
        rtm(1:nzmax)     = 0.0         ! total water mixing ratio
        wprtp(1:nzmax)   = 0.0         ! w'rt'
        wpthlp(1:nzmax)  = 0.0         ! w'thl'
        wp3(1:nzmax)     = 0.0         ! w'^3
        rtp2(1:nzmax)    = rttol**2    ! rt'^2
        thlp2(1:nzmax)   = thltol**2   ! thl'^2
        rtpthlp(1:nzmax) = 0.0         ! rt'thl'
 
        p_in_Pa(1:nzmax)= 0.0    ! pressure (Pa)
        exner(1:nzmax) = 0.0    ! exner
        rho(1:nzmax)  = 0.0    ! density on thermo. levels
        rho_zm(1:nzmax)  = 0.0    ! density on moment. levels

        thlm_forcing(1:nzmax) = 0.0     ! thlm large-scale forcing
        rtm_forcing(1:nzmax)  = 0.0     ! rtm large-scale forcing

        ! Imposed large scale w

        wm_zm(1:nzmax) = 0.0      ! Momentum levels
        wm_zt(1:nzmax) = 0.0      ! Thermodynamic levels

        ! PDF width parameter: momentum levels

        sigma_sqd_w(1:nzmax)  = 0.0 

        ! Mixing lengths

        tau_zm(1:nzmax) = 0.0

        ! Cloud water variables

        rcm(1:nzmax)  = 0.0
        cf(1:nzmax)   = 0.0


        ! Surface fluxes
        wpthlp_sfc = 0.0
        wprtp_sfc  = 0.0
        upwp_sfc   = 0.0
        vpwp_sfc   = 0.0

        ! Passive scalars 
        do i = 1, sclr_dim, 1
          wpsclrp_sfc(i)   = 0.0
          wpedsclrp_sfc(i) = 0.0

          sclrm(1:nzmax,i)         = 0.0
          sclrm_forcing(1:nzmax,i) = 0.0

          edsclrm(1:nzmax,i) = 0.0
          wpsclrp(1:nzmax,i) = 0.0
        end do

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

        deallocate( p_in_Pa )   ! pressure
        deallocate( exner )     ! exner
        deallocate( rho )      ! density: t points
        deallocate( rho_zm )      ! density: m points

        deallocate( thlm_forcing )
        deallocate( rtm_forcing )

        ! Imposed large scale w

        deallocate( wm_zm )       ! momentum levels
        deallocate( wm_zt )       ! thermodynamic levels

        ! PDF width parameter

        deallocate( sigma_sqd_w )

        ! Mixing lengths

        deallocate( tau_zm )

        ! Cloud water variables

        deallocate( rcm )
        deallocate( cf )


        ! Passive scalars
        deallocate( wpsclrp_sfc, wpedsclrp_sfc )
        deallocate( sclrm )
        deallocate( sclrm_forcing )

        deallocate( edsclrm )
        deallocate( wpsclrp )

!$omp   end parallel

        return
        end subroutine cleanup_prognostic_variables

        end module prognostic_variables
