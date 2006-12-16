#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: astex.F,v 1.1 2006-12-16 00:32:45 dschanen Exp $
        module astex

!       Description:
!       Contains subroutines for the ASTEX KK case.
!----------------------------------------------------------------------

        public :: astex_tndcy, astex_sfclyr

        private

        contains

!----------------------------------------------------------------------
        subroutine astex_tndcy( time, wmt, wmm, theta_tndcy, rt_tndcy,
     .                          sclr_tndcy )

!       Description:
!       Subroutine to set theta and water tendencies for ASTEX KK case
!       References:
!----------------------------------------------------------------------

        use constants, only: sclr_dim
        use grid_class

        implicit none

!       Input
        real, intent(in) :: time

!       Output
        real, intent(out) :: wmt(gr%nnzp)
        real, intent(out) :: wmm(gr%nnzp)
        real, intent(out) :: theta_tndcy(gr%nnzp)
        real, intent(out) :: rt_tndcy(gr%nnzp)

        real, optional, 
     .        intent(out) :: sclr_tndcy(gr%nnzp, sclr_dim)

!       Internal variables

        integer i

!       Large-scale subsidence

        do i=2,gr%nnzp

           wmt(i) = - 5.e-6 * gr%zt(i)

        end do

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
!       Radiative theta-l tendency

        theta_tndcy = 0.

!       Large scale advective moisture tendency

        rt_tndcy = 0.

        if (present( sclr_tndcy )) sclr_tndcy = 0.0

        return
        end subroutine astex_tndcy

!----------------------------------------------------------------------
        subroutine astex_sfclyr( rho0, up, vp, uw, vw, wtp, wrp,
     .                           wpsclrp_sfc )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to ASTEX with Khairoutdinov and Kogan
!       alteration.

!       References:
!----------------------------------------------------------------------

        use constants

        implicit none

!       Input variables

        real, intent(in)  :: rho0, up, vp

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp

!       Output variables (optional) 

        real, optional, intent(out) :: wpsclrp_sfc(sclr_dim)

!       Internal variables

        real HF, LF

! Compute heat and moisture fluxes

        HF = 10.0
        LF = 25.0

        wtp = HF/(rho0*Cp)
        wrp = LF/(rho0*Lv)

! Compute momentum fluxes

        uw = 0.09
        vw = 0.09

        if (present( wpsclrp_sfc )) 
     .    wpsclrp_sfc(SCLR_THETA) = HF/(rho0*Cp)
        if (present( wpsclrp_sfc )) wpsclrp_sfc(SCLR_RT) = LF/(rho0*Lv)

        return
        end subroutine astex_sfclyr

        end module astex
