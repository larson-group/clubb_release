#define SCLR_THETA 1
#define SCLR_RT 1
!----------------------------------------------------------------------
!$Id: atex.F,v 1.2 2006-12-20 03:16:08 dschanen Exp $
        module atex

!       Description:
!       Contains subroutines for the GCSS ATEX case.
!----------------------------------------------------------------------
        public :: atex_tndcy, atex_sfclyr

        contains

!----------------------------------------------------------------------
        subroutine atex_tndcy( time, time_initial, rtm, 
     .                         rhot, rcm, exner, wmt, wmm, Frad, radht,
     .                         thlm_forcing, rtm_forcing, isValid,
     .                         sclrm_forcing )
!       Description:
!       Subroutine to set theta-l and water tendencies for ATEX case

!       References:

!----------------------------------------------------------------------
        use constants, only: sclr_dim, fstderr
        use model_flags, only: lbugsrad
        use grid_class
        use atex_cloud_rad
#ifdef STATS
        use statistics
#endif

        implicit none

        ! Input Variables
        real, intent(in) :: 
     .  time,        ! Current time     [s]
     .  time_initial ! Initial time     [s]

        real, intent(in), dimension(gr%nnzp) ::
     .  rtm,  ! Total water mixing ratio        [kg/kg]
     .  rhot, ! Density                         [kg/m^3]
     .  rcm,  ! Liquid water mixing ratio       [kg/kg]
     .  exner ! Exner function                  [-]

        ! Input/output
        logical, intent(inout) :: isValid ! Diagnostic 

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) ::
     .  wmt,         ! w wind on thermodynamic grid                [m/s]
     .  wmm,         ! w wind on momentum grid                     [m/s]
     .  Frad,        ! Radiative flux                              [W/m^2]
     .  radht,       ! Radiative heating rate                      [K/s]
     .  thlm_forcing,! Liquid water potential temperature tendency [K/s]
     .  rtm_forcing  ! Total water mixing ratio tendency           [kg/kg/s]


        ! Output (optional)
        real, optional, intent(out), dimension(gr%nnzp, sclr_dim) ::
     .  sclrm_forcing ! Passive scalar tendency         [units/s]


        ! Internal variables
        integer :: i
        real :: zi

        ! Forcings are applied only after t = 5400 s
        wmt = 0.
        wmm = 0.

        thlm_forcing = 0.
        rtm_forcing  = 0.

        if ( time >= time_initial + 5400.0 ) then

        !  Identify height of 6.5 g/kg moisture level

           i = 2
           do while ( i <= gr%nnzp .and. rtm(i) > 6.5e-3 )
              i = i + 1
           end do
           if ( i == gr%nnzp+1 .or. i == 2 ) then
             write(fstderr,*) "Identification of 6.5 g/kg level failed"
             write(fstderr,*) "Subroutine: atex_tndcy. File: atex.F"
             write(fstderr,*) "i = ", i
             write(fstderr,*) "rtm(i) = ",rtm(i)
             isValid = .false.
             return
           end if
           zi = gr%zt(i-1)

!          Large scale subsidence

           do i = 2, gr%nnzp

              if ( gr%zt(i) > 0. .and. gr%zt(i) <= zi ) then
                 wmt(i) 
     .             = -0.0065 * gr%zt(i)/zi
              else if ( gr%zt(i) > zi .and. gr%zt(i) <= zi+300. ) then
                 wmt(i)
     .             = - 0.0065 * ( 1. - (gr%zt(i)-zi)/300. )
              else
                 wmt(i) = 0.
              end if

           end do

           wmm = zt2zm( wmt )

           ! Boundary conditions.
           wmt(1) = 0.0        ! Below surface
           wmm(1) = 0.0        ! At surface
           wmm(gr%nnzp) = 0.0  ! Model top
        
           ! Theta-l tendency

           do i = 2, gr%nnzp

              if ( gr%zt(i) > 0. .and. gr%zt(i) < zi ) then
                 thlm_forcing(i) 
     .             = -1.1575e-5 * ( 3. - gr%zt(i)/zi )
              else if ( gr%zt(i) > zi .and. gr%zt(i) <= zi+300. ) then
                 thlm_forcing(i)
     .             = -2.315e-5 * ( 1. - (gr%zt(i)-zi)/300. )
              else
                 thlm_forcing(i) = 0.0
              end if

           end do

           ! Moisture tendency
           do i = 2, gr%nnzp

              if ( gr%zt(i) > 0. .and. gr%zt(i) < zi ) then
                 rtm_forcing(i) = -1.58e-8 * ( 1. - gr%zt(i)/zi )  ! Brian
              else
                 rtm_forcing(i) = 0.0       ! Brian
              end if

           end do

           ! Boundary conditions
           thlm_forcing(1) = 0.0  ! Below surface
           rtm_forcing(1)  = 0.0  ! Below surface

        end if ! time >= time_initial + 5400.0

        ! Use cloud_rad() to compute radiation
        if ( .not. lbugsrad ) then
          call cloud_rad( rhot, rcm, exner, Frad, radht, thlm_forcing )
        end if

#ifdef STATS
        if ( .not. lbugsrad .and. lstats_samp ) then
          if ( iradht_LW > 0 ) then
            zt%x(:,iradht_LW) = zt%x(:,iradht_LW) + radht
            zt%n(:,iradht_LW) = zt%n(:,iradht_LW) + 1
          end if
        end if
#endif
        if ( present( sclrm_forcing ) ) then
          sclrm_forcing(1:gr%nnzp,SCLR_THETA) = thlm_forcing(1:gr%nnzp)
          sclrm_forcing(1:gr%nnzp,SCLR_RT)    = rtm_forcing(1:gr%nnzp)
        end if

        return
        end subroutine atex_tndcy

!----------------------------------------------------------------------
        subroutine atex_sfclyr( up, vp, thp, rtp, uw, vw, wtp, wrp,
     .                          sclrm_sfc, wpsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS ATEX specifications

!       References:

!----------------------------------------------------------------------

        use constants

        implicit none

        ! Constants

        real, parameter :: 
     .  ubmin = 0.25,
     .  ustar = 0.3,
     .  C_10 = 0.0013,
     .  SST = 298 

        ! Input variables
        real, intent(in) :: up, vp, thp, rtp

        ! Input variables (optional)
        real, intent(in), optional, dimension(sclr_dim) :: sclrm_sfc

        ! Output variables
        real, intent(out) :: uw, vw, wtp, wrp

        ! Output variables (optional)

        real, optional, intent(out) :: wpsclrp_sfc(sclr_dim)

        ! Local Variables
        real :: ubar

        ! Compute heat and moisture fluxes

        ubar = max( ubmin, sqrt( up*up + vp*vp ) )

        wtp  = -C_10 * ubar * ( thp - SST * (1000./1015.)**kappa )
        wrp  = -C_10 * ubar * ( rtp - 0.0198293 )

        ! Optional passive scalar
        if ( present( wpsclrp_sfc ) .and. present( sclrm_sfc ) ) then
          wpsclrp_sfc(SCLR_THETA) = -C_10 * ubar * 
     .           ( sclrm_sfc(SCLR_THETA) - SST * (1000./1015.)**kappa )

          wpsclrp_sfc(SCLR_RT) = -C_10 * ubar * 
     .                     ( sclrm_sfc(SCLR_RT) - 0.0198293 )
        end if

        ! Compute momentum fluxes

        uw = -up * ustar*ustar / ubar
        vw = -vp * ustar*ustar / ubar

        return
        end subroutine atex_sfclyr

        end module atex
