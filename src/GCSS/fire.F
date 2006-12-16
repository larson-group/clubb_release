#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: fire.F,v 1.1 2006-12-16 00:32:45 dschanen Exp $
        module fire

!       Description:
!       Contains subroutines for the GCSS FIRE case.
!----------------------------------------------------------------------
        public :: fire_tndcy

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

        use constants, only: sclr_dim
        use model_flags, only: lbugsrad
        use grid_class
        use atex_cloud_rad

        implicit none

        ! Input Variables
        real, intent(in) :: time ! Model time [s]

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

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
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

        if ( present( sclrm_forcing ) ) then
          sclrm_forcing(1:gr%nnzp,SCLR_THETA) = thlm_forcing(1:gr%nnzp)
          sclrm_forcing(1:gr%nnzp,SCLR_RT)    = rtm_forcing(1:gr%nnzp)
        end if

        return
        end subroutine fire_tndcy

        end module fire
