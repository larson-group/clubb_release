#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: astex.F,v 1.2 2006-12-20 03:16:08 dschanen Exp $
        module astex

!       Description:
!       Contains subroutines for the ASTEX KK case.
!----------------------------------------------------------------------

        public :: astex_tndcy, astex_sfclyr

        private

        contains

!----------------------------------------------------------------------
        subroutine astex_tndcy( time, wmt, wmm, 
     .                          thlm_forcing, rtm_forcing,
     .                          sclrm_forcing )

!       Description:
!       Subroutine to set theta and water tendencies for ASTEX KK case
!       References:
!----------------------------------------------------------------------

        use constants, only: sclr_dim
        use grid_class

        implicit none

        ! Input Variables
        real, intent(in) :: time ! Current time         [s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) :: 
     .  wmt,          ! w wind on the thermodynamic grid        [m/s]
     .  wmm,          ! w wind on the momentum grid             [m/s]
     .  thlm_forcing, ! Liquid potential temperature tendency   [K/s]
     .  rtm_forcing   ! Total water mixing ratio tendency       [kg/kg/s]

        real, intent(out), optional, dimension(gr%nnzp,sclr_dim) :: 
     .  sclrm_forcing ! Passive scalar forcing  [units/s]

        ! Local variables

        integer :: i

        ! Large-scale subsidence

        do i=2,gr%nnzp

           wmt(i) = - 5.e-6 * gr%zt(i)

        end do

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
        ! Radiative theta-l tendency

        thlm_forcing = 0.0

        ! Large scale advective moisture tendency

        rtm_forcing = 0.0

        ! Passive scalar testing

        if ( present( sclrm_forcing ) ) sclrm_forcing(:,:) = 0.0

        return
        end subroutine astex_tndcy

!----------------------------------------------------------------------
        subroutine astex_sfclyr( rho0, up, vp, upwp_sfc, vpwp_sfc, 
     .                           wpthlp_sfc, wprtp_sfc, wpsclrp_sfc )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to ASTEX with Khairoutdinov and Kogan
!       alteration.

!       References:
!----------------------------------------------------------------------

        use constants

        implicit none

        ! Input variables

        real, intent(in) :: 
     .  rho0,    ! Density at (1)          [kg/m^3]
     .  up, vp   ! um(2), vm(2) (not used) [m/s]

        ! Output variables

        real, intent(out) :: 
     .  upwp_sfc,    ! u'w' at (1)      [m^2/s^2]
     .  vpwp_sfc,    ! v'w'at (1)       [m^2/s^2]
     .  wpthlp_sfc,  ! w'th_l' at (1)   [(m K)/s]  
     .  wprtp_sfc    ! w'r_t'(1) at (1) [(m kg)/(s kg)]

        ! Output variables (optional) 

        real, intent(out), optional, dimension(sclr_dim) :: 
     .  wpsclrp_sfc  ! w' scalar at (1) [units vary]

        ! Local variables

        real :: HF, LF

        ! Compute heat and moisture fluxes

        HF = 10.0
        LF = 25.0

        wpthlp_sfc = HF/( rho0*Cp )
        wprtp_sfc  = LF/( rho0*Lv )

        ! Compute momentum fluxes

        upwp_sfc = 0.09
        vpwp_sfc = 0.09

        ! Test scalars
        if ( present( wpsclrp_sfc ) ) then
          wpsclrp_sfc(SCLR_THETA) = HF/( rho0*Cp )
          wpsclrp_sfc(SCLR_RT)    = LF/( rho0*Lv )
        end if

        return
        end subroutine astex_sfclyr

        end module astex
