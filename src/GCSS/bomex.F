#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: bomex.F,v 1.5 2007-12-15 23:58:04 griffinb Exp $
        module bomex

!       Description:
!       Contains subroutines for the GCSS BOMEX case.
!----------------------------------------------------------------------
        public :: bomex_tndcy, bomex_sfclyr

        contains

!----------------------------------------------------------------------
        subroutine bomex_tndcy( time, wmt, wmm, radht,
     .                          thlm_forcing, rtm_forcing,
     .                          sclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for BOMEX case

!       References:
!       <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!----------------------------------------------------------------------

        use grid_class
        use model_flags, only: lbugsrad
        use constants, only: sclr_dim

        implicit none

        ! Input
        real, intent(in)  :: time ! Current model time [s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) ::
     .  wmt,          ! w wind on thermodynamic grid                 [m/s]
     .  wmm,          ! w wind on momentum grid                      [m/s]
     .  radht,        ! Radiative heating rate                       [K/s]
     .  thlm_forcing, ! Liquid water potential temperature tendency  [K/s]
     .  rtm_forcing   ! Total water mixing ratio tendency            [kg/kg/s]

        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclrm_forcing ! Passive scalar forcing [units vary]

        ! Local Variables
        integer :: k

        ! Large scale subsidence
        do k = 2, gr%nnzp, 1

           if ( gr%zt(k) >= 0. .and. gr%zt(k) < 1500. ) then
              wmt(k) = - ( 0.0065 / 1500. ) * gr%zt(k)
           else if ( gr%zt(k) >= 1500. .and. gr%zt(k) < 2100. ) then
              wmt(k)
     .          = - 0.0065 
     .            + 0.0065 * ( gr%zt(k) - 1500. ) / ( 2100. - 1500. )
           else
              wmt(k) = 0.
           end if

        end do ! k=2..gr%nnzp

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
        if ( .not. lbugsrad ) then

        ! Radiative theta-l tendency
          do k = 2, gr%nnzp

            if ( gr%zt(k) >= 0. .and. gr%zt(k) < 1500. ) then
              radht(k) = -2.315e-5
            else if ( gr%zt(k) >= 1500. .and. gr%zt(k) < 2500. ) then
              radht(k)
     .          = - 2.315e-5 
     .            + 2.315e-5 
     .              * ( gr%zt(k) - 1500. ) / ( 2500. - 1500. )
            else
              radht(k) = 0.
            end if

          end do

          thlm_forcing = radht
        else ! Compute radht interactively with BUGSrad

          thlm_forcing = 0.0

        end if ! ~lbugsrad

        ! Large scale advective moisture tendency
        do k = 2, gr%nnzp

           if ( gr%zt(k) >= 0. .and. gr%zt(k) < 300. ) then
              rtm_forcing(k) = -1.2e-8
           else if ( gr%zt(k) >= 300. .and. gr%zt(k) < 500. ) then
              rtm_forcing(k) 
     .          = - 1.2e-8 
     .              * ( 1. - ( gr%zt(k) - 300. )/( 500. - 300. ) )
           else
              rtm_forcing(k) = 0.
           end if

        end do


        ! Boundary conditions
        thlm_forcing(1) = 0.0  ! Below surface
        rtm_forcing(1)  = 0.0  ! Below surface

        if ( present( sclrm_forcing ) ) then 
          sclrm_forcing(1:gr%nnzp,SCLR_THETA) = thlm_forcing(1:gr%nnzp)
          sclrm_forcing(1:gr%nnzp,SCLR_RT)    = rtm_forcing(1:gr%nnzp)
        end if

        return
        end subroutine bomex_tndcy

!----------------------------------------------------------------------
        subroutine bomex_sfclyr( um_sfc, vm_sfc, 
     .                           upwp_sfc, vpwp_sfc,
     .                           wpthlp_sfc, wprtp_sfc, ustar,
     .                           wpsclrp_sfc, wpedsclrp_sfc )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS BOMEX specifications

!       References:
!----------------------------------------------------------------------

        use constants, only: sclr_dim

        implicit none

        ! Constant Parameters
        real, parameter :: 
     .  ubmin = 0.25
!     .  ustar = 0.28

        ! Input Variables
        real, intent(in) :: 
     .  um_sfc, ! um(2) [m/s]
     .  vm_sfc  ! vm(2) [m/s]

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

        ! Declare the value of ustar.
        ustar = 0.28

        ! Compute heat and moisture fluxes

        wpthlp_sfc = 8.e-3
        wprtp_sfc  = 5.2e-5

        ! Let passive scalars be equal to rt and theta_l for now 
        if ( present( wpsclrp_sfc ) ) then 
          wpsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpsclrp_sfc(SCLR_RT)    = wprtp_sfc
        end if

        if ( present( wpsclrp_sfc ) ) then 
          wpedsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpedsclrp_sfc(SCLR_RT)    = wprtp_sfc
        end if

        ! Compute momentum fluxes

        ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

        upwp_sfc = -um_sfc * ustar**2 / ubar
        vpwp_sfc = -vm_sfc * ustar**2 / ubar

        return
        end subroutine bomex_sfclyr

        end module bomex
