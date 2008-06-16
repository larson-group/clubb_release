#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: fire.F,v 1.11 2008-06-16 22:26:11 faschinj Exp $
        module fire

!       Description:
!       Contains subroutines for the GCSS FIRE case.
!----------------------------------------------------------------------
        public :: fire_tndcy, sfc_momentum_fluxes, sfc_thermo_fluxes

        private ! Default Scope

        contains

!----------------------------------------------------------------------
        subroutine fire_tndcy
     .             ( time, rhot, rcm, exner, 
     .               wmt, wmm, Frad, radht, 
     .               thlm_forcing, rtm_forcing,
     .               sclrm_forcing )
!       Description:
!       Subroutine to large-scale subsidence for FIRE case. Calls
!       cloud_rad for computing radiation

!       References:
!       None
!----------------------------------------------------------------------

        use constants, only: sclr_dim ! Variable(s)

        use model_flags, only: lbugsrad  ! Variable(s)

        use grid_class, only: gr ! Variable(s)

        use grid_class, only: zt2zm ! Procedure(s)

        use atex_cloud_rad, only: cloud_rad ! Procedure(s)

        use stats_prec, only: time_prec ! Variable(s)

#ifdef STATS
        use stats_type, only: stat_update_var ! Procedure(s)

        use stats_hoc, only: zt, iradht_LW, lstats_samp ! Variable(s)
#endif

        implicit none

        ! Input Variables
        real(kind=time_prec), intent(in) :: time ! Model time [s]

        real, intent(in), dimension(gr%nnzp) ::
     .  rhot, ! Density                         [kg/m^3]
     .  rcm,  ! Liquid water mixing ratio       [kg/kg]
     .  exner ! Exner function                  [-]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) ::
     .  wmt,         ! w wind on thermodynamic grid     [m/s]
     .  wmm,         ! w wind on momentum grid          [m/s]
     .  Frad,        ! Radiative flux                   [W/m^2]
     .  radht,       ! Radiative heating rate           [K/s]
     .  thlm_forcing,! Liquid water potential temperature tendency [K/s]
     .  rtm_forcing  ! Total water mixing ratio tendency [kg/kg/s]

        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclrm_forcing ! Passive scalar tendency [units vary]

!       Internal variables

        integer :: k

!       Large-scale subsidence

        do k = 2, gr%nnzp

           if ( gr%zt(k) >= 0. .and. gr%zt(k) < 1500. ) then
              wmt(k) = - 5.e-6 * gr%zt(k)
           end if

        end do

        ! Boundary condition.
        wmt(1) = 0.0        ! Below surface

        wmm = zt2zm( wmt )

        ! Boundary condition.
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
        ! Radiative theta-l tendency is computed interactively elsewhere

        thlm_forcing = 0.0

        ! Large scale advective moisture tendency

        rtm_forcing = 0.0

        ! Use cloud_rad to compute radiation
        if ( .not. lbugsrad ) then
          call cloud_rad( rhot, rcm, exner, Frad, radht, thlm_forcing )
        end if

#ifdef STATS
        if ( .not. lbugsrad .and. lstats_samp ) then
          call stat_update_var( iradht_LW, radht, zt )
        end if
#endif
        ! Testing for scalars
        if ( present( sclrm_forcing ) ) then
          sclrm_forcing(1:gr%nnzp,SCLR_THETA) = thlm_forcing(1:gr%nnzp)
          sclrm_forcing(1:gr%nnzp,SCLR_RT)    = rtm_forcing(1:gr%nnzp)
        end if

        return
        end subroutine fire_tndcy

!------------------------------------------------------------------------
        subroutine sfc_momentum_fluxes( u, v, upwp_sfc, vpwp_sfc, 
     .                                  ustar )

!       Description:
!       This subroutine computes surface momentum fluxes using aerodynamic
!       formulas.

!       References:
!       None
!------------------------------------------------------------------------

        implicit none

        ! External
        intrinsic :: sqrt

        ! Constant parameter
!        real, intent(out) :: 
!     .  ustar = 0.3

        ! Input variables
        real, intent(in) :: 
     .  u, ! u wind first level above ground    [m/s]
     .  v  ! v wind first level above ground    [m/s]

        ! Output Variables
        real, intent(out) :: 
     .  upwp_sfc,! sfc u momentum flux (m^2/s^2)
     .  vpwp_sfc,! sfc v momentum flux (m^2/s^2)
     .  ustar    ! surface friction velocity [m/s]

        ! Local Variables
        real :: M ! total wind speed above ground

        ! Declare the value of ustar
        ustar = 0.3

        ! Computes fluxes

        M = sqrt( u*u + v*v )
        upwp_sfc = - ustar*ustar * u / M
        vpwp_sfc = - ustar*ustar * v / M

        return
        end subroutine sfc_momentum_fluxes

!------------------------------------------------------------------------
        subroutine sfc_thermo_fluxes( u, v, Tsfc, psfc, thlair, rtair,
     .                                wpthlp_sfc, wprtp_sfc,
     .                                sclrair, wpsclrp_sfc,
     .                                wpedsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of heat and moisture 
!       using aerodynamic formulas.

!       References:
!       None
!------------------------------------------------------------------------

        use constants, only: sclr_dim, kappa, p0 ! Variable(s)

        use saturation, only: rsat ! Procedure(s)

        implicit none
        
        ! External
        intrinsic :: present, sqrt

        

        ! Parameter
        real, parameter :: C = 1.3e-3

        ! Input Variables
        real, intent(in) :: 
     .  u,      ! u wind                        [m/s]
     .  v,      ! u wind                        [m/s]
     .  Tsfc,   ! Surface temperature           [K]
     .  psfc,   ! Surface pressure              [Pa]
     .  thlair, ! theta_l at first model layer  [K]
     .  rtair   ! rt at first model layer       [kg/kg]


        ! Output Variables
        real, intent(out) :: 
     .  wpthlp_sfc,! surface thetal flux        [K m/s]
     .  wprtp_sfc  ! surface moisture flux      [kg/kg m/s]

        ! Input Variables (optional)
        real, optional, intent(in), dimension(sclr_dim) :: 
     .  sclrair ! Passive scalars at first model layer  [units vary]

        ! Output Variables (optional) 
        real, optional, intent(out), dimension(sclr_dim) :: 
     .  wpsclrp_sfc,   ! scalar surface flux            [units m/s]
     .  wpedsclrp_sfc  ! eddy-scalar surface flux       [units m/s]

        ! Local Variables
        real :: M  ! Total wind speed above ground

        ! Compute fluxes
        M = sqrt( u*u + v*v )
        wpthlp_sfc = -C * M * ( thlair - Tsfc * (psfc/p0)**kappa )
        wprtp_sfc  = -C * M * ( rtair - rsat( psfc, Tsfc ) )

        ! Compute passive scalar as == to wprtp & wpthlp for now
        if ( present( wpsclrp_sfc ) ) then
          wpsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpsclrp_sfc(SCLR_RT)    = wprtp_sfc
        endif

        if ( present( wpedsclrp_sfc ) ) then
          wpedsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpedsclrp_sfc(SCLR_RT)    = wprtp_sfc
        endif

        return
        end subroutine sfc_thermo_fluxes

        end module fire
