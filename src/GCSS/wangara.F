#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: wangara.F,v 1.1 2006-12-16 00:32:45 dschanen Exp $
        module wangara

        public :: wangara_tndcy, wangard_sfclyr

        contains
!----------------------------------------------------------------------
        subroutine wangara_tndcy( time, wmt, wmm, 
     .                            thlm_forcing, rtm_forcing,
     .                            sclr_tndcy )
!       Description:
!       Subroutine to set theta and water tendencies for Wangara case
!       References;
!       None
!----------------------------------------------------------------------
        use grid_class
        use constants, only: sclr_dim

        implicit none

        ! Input Variables
        real, intent(in) :: time ! Model time [s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) ::
     .  wmt,         ! w wind on thermodynamic grid                [m/s]
     .  wmm,         ! w wind on momentum grid                     [m/s]
     .  thlm_forcing,! Liquid water potential temperature tendency [K/s]
     .  rtm_forcing  ! Total water mixing ratio tendency           [kg/kg/s]

        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclr_tndcy ! Passive scalar tendency [units vary]

        ! No large-scale subsidence for now
        wmt = 0.0
        wmm = 0.0

        ! No large-scale water tendency or cooling

        rtm_forcing  = 0.0
        thlm_forcing = 0.0

        if ( present( sclr_tndcy ) ) then
          sclr_tndcy(1:gr%nnzp,SCLR_THETA) = thlm_forcing(1:gr%nnzp)
          sclr_tndcy(1:gr%nnzp,SCLR_RT)    = rtm_forcing(1:gr%nnzp)
        end if

        return
        end subroutine wangara_tndcy

!----------------------------------------------------------------------
        subroutine wangara_sfclyr( time, up, vp, uw, vw, wtp, wrp,
     .                             wpsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture for Wangara day 33
!----------------------------------------------------------------------

        use constants, only: sclr_dim, pi, fstderr

        implicit none

        intrinsic :: mod, max, cos, sqrt, present

        ! Constants
        real, parameter :: ubmin = 0.25
        real, parameter :: ustar = 0.13

        ! Input variables

        real, intent(in) :: time
        real, intent(in) :: up, vp

        ! Output variables

        real, intent(out) :: uw, vw, wtp, wrp
        
        real, optional, intent(out), dimension(sclr_dim) :: 
     .  wpsclrp_sfc

!       Internal variables

        real :: ubar
        real :: time_utc, time_est

        ! Compute UTC time of the day in seconds

        time_utc = mod( time, 86400. )

        ! Now convert UTC time to Australia EST (local time)

        time_est = mod( time_utc + 36000., 86400. )

        if ( time_est < 27000 .or. time_est > 63000 ) then
           write(fstderr,*) "wangara_sfclyr: error local time must"
     .       //" be between 730 and 1730."
           write(fstderr,*) 'time_est = ',time_est
           stop
        end if

        ! Compute heat and moisture fluxes

        wtp = 0.18
     .        * cos( (time_est-45000.0)/36000.0 * pi )
        wrp = 1.3e-4 * wtp

        ! Compute new mixing variables
        if ( present(wpsclrp_sfc) ) wpsclrp_sfc(SCLR_THETA) = 0.18
     .        * cos( (time_est-45000.0)/36000.0 * pi )

        if ( present(wpsclrp_sfc) ) wpsclrp_sfc(SCLR_RT) = 1.3e-4*wtp

        ! Compute momentum fluxes

        ubar = max(ubmin, sqrt(up*up + vp*vp))

        uw = -up * ustar*ustar / ubar
        vw = -vp * ustar*ustar / ubar

        return
        end subroutine wangara_sfclyr

!----------------------------------------------------------------------

        end module wangara
