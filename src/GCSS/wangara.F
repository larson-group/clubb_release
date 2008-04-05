#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: wangara.F,v 1.6 2008-04-05 19:54:20 griffinb Exp $
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
        use stats_prec, only: time_prec

        implicit none

        ! Input Variables
        real(kind=time_prec), intent(in) :: time ! Model time [s]

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
        subroutine wangara_sfclyr( time, um_sfc, vm_sfc, 
     .                             upwp_sfc, vpwp_sfc, 
     .                             wpthlp_sfc, wprtp_sfc, ustar,
     .                             wpsclrp_sfc, wpedsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture for Wangara day 33
!----------------------------------------------------------------------

        use constants, only: sclr_dim, pi, fstderr, sec_per_day
        use stats_prec, only: time_prec

        implicit none

        intrinsic :: mod, max, cos, sqrt, present

        ! Constants
        real, parameter :: 
     .  ubmin = 0.25
!     .  ustar = 0.13

        ! Input variables
        real(kind=time_prec), intent(in) :: 
     .  time    ! Current time  [s]

        real, intent(in) :: 
     .  um_sfc, ! um(2)         [m/s]
     .  vm_sfc  ! vm(2)         [m/s]

        ! Output variables
        real, intent(out) :: 
     .  upwp_sfc,    ! u'w' at (1)      [m^2/s^2]
     .  vpwp_sfc,    ! v'w'at (1)       [m^2/s^2]
     .  wpthlp_sfc,  ! w'th_l' at (1)   [(m K)/s]  
     .  wprtp_sfc,   ! w'r_t'(1) at (1) [(m kg)/(s kg)]
     .  ustar        ! surface friction velocity [m/s]

        ! Output variables (optional)
        real, intent(out), optional, dimension(sclr_dim) :: 
     .  wpsclrp_sfc,  ! Passive scalar surface flux      [units m/s]
     .  wpedsclrp_sfc ! Passive eddy-scalar surface flux [units m/s]

        ! Local variables

        real :: ubar 
        real(kind=time_prec) :: time_utc, time_est


        ! Declare the value of ustar.
        ustar = 0.13

        ! Compute UTC time of the day in seconds

        time_utc = mod( time, sec_per_day )

        ! Now convert UTC time to Australia EST (local time)

        time_est = mod( time_utc + 36000., sec_per_day )

        if ( time_est < 27000 .or. time_est > 63000 ) then
           write(fstderr,*) "wangara_sfclyr: error local time must"
     .       //" be between 730 and 1730."
           write(fstderr,*) 'time_est = ',time_est
           stop
        end if

        ! Compute heat and moisture fluxes

        wpthlp_sfc = 0.18 * cos( (time_est-45000.0)/36000.0 * pi )
        wprtp_sfc  = 1.3e-4 * wpthlp_sfc

        ! Let the passive scalars be equal to rt and theta_l for now
        if ( present( wpsclrp_sfc ) ) then
          wpsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpsclrp_sfc(SCLR_RT)    = wprtp_sfc
        end if
        if ( present( wpedsclrp_sfc  ) ) then
          wpedsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpedsclrp_sfc(SCLR_RT)    = wprtp_sfc
        end if

        ! Compute momentum fluxes

        ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

        upwp_sfc = -um_sfc * ustar**2 / ubar
        vpwp_sfc = -vm_sfc * ustar**2 / ubar

        return
        end subroutine wangara_sfclyr

!----------------------------------------------------------------------

        end module wangara
