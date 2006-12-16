#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: arm.F,v 1.1 2006-12-16 00:32:45 dschanen Exp $
        module arm

!       Description:
!       Contains subroutines for the GCSS ARM case.
!----------------------------------------------------------------------

        public :: arm_tndcy, arm_sfclyr

        private

        contains

!----------------------------------------------------------------------
        subroutine arm_tndcy( time, theta_tndcy, rt_tndcy,
     .                        sclr_tndcy )
!       Description:
!       Subroutine to set theta and water tendencies for ARM case

!       References:
!       None
!----------------------------------------------------------------------

        use grid_class
        use constants, only: sclr_dim
        use model_flags, only: lbugsrad

        implicit none

        ! External
        intrinsic ::
     .  int, min, max, present

        ! Constant Parameters
        real, parameter, dimension(6) :: 
     .  atheta = (/ 0.000, 0.000,  0.000, -0.080, -0.160, -0.160/),
     .  rtheta = (/-0.125, 0.000,  0.000,  0.000,  0.000, -0.100/),
     .  art    = (/ 0.080, 0.020, -0.040, -0.100, -0.160, -0.300/)

        ! Input Variables
        real, intent(in) :: time ! Model time [s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) :: 
     .  theta_tndcy, ! Liquid water potential temperature tendency [K/s]
     .  rt_tndcy     ! Total water mixing ratio tendency [kg/kg/s]

        ! Output (optional) Variables
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) :: 
     .  sclr_tndcy ! Passive scalar tendency [units vary]

        ! Local variables
        integer :: i, i1, i2 ! Loop indices
        real :: 
     .  a, b,      ! [-]
     .  true_time, ! [s]
     .  theta_tmp, ! [K/s]
     .  rt_tmp     ! [kg/kg/s]

!-----------------------------------------------------------------------

        true_time = time

        ! Interpolate in time to get theta and rt tendency

        i1 = int( ( true_time - 41400. ) / 10800. ) + 1
        i1 = min( max( i1, 1 ), 5 )
        i2 = i1 + 1

        if (i1 < 5) then
           a = ( true_time - (41400. + 10800. * (i1-1)) ) / 10800.
        else
           a = ( true_time - (41400. + 10800. * (i1-1)) ) / 9000.
        end if

        if ( .not. lbugsrad ) then
          theta_tmp = ( 1. - a ) * ( atheta(i1) + rtheta(i1) )
     .              + a * ( atheta(i2) + rtheta(i2) )
        else ! factor in radiation later
          theta_tmp = ( 1. - a ) * ( atheta(i1) + 0.0 )
     .              + a * ( atheta(i2) + 0.0 )
        end if

        rt_tmp = ( 1. - a ) * art(i1) + a * art(i2)

        ! Convert to the right units

        theta_tmp = theta_tmp / 3600.
        rt_tmp    = rt_tmp / ( 3600. * 1000. )

        ! Interpolate with respect to height

        do i = 2, gr%nnzp
          if ( gr%zt(i) >= 0. .and. gr%zt(i) < 1000. ) then
            rt_tndcy(i)    = rt_tmp
            theta_tndcy(i) = theta_tmp
          else if ( gr%zt(i) >= 1000. .and. gr%zt(i) < 3000. ) then
            b              = 1. - ( gr%zt(i) - 1000. ) / 2000.
            rt_tndcy(i)    = b * rt_tmp
            theta_tndcy(i) = b * theta_tmp
          else
            rt_tndcy(i)    = 0.0
            theta_tndcy(i) = 0.0
          end if
        end do

        if ( present( sclr_tndcy ) ) then
          sclr_tndcy(1:gr%nnzp,SCLR_THETA) = theta_tndcy
          sclr_tndcy(1:gr%nnzp,SCLR_RT)    = rt_tndcy
        end if

        return
        end subroutine arm_tndcy
!----------------------------------------------------------------------
        subroutine arm_sfclyr( time, z, dn0, thp, up, vp, uw, vw, 
     .                         wtp, wrp, ustar,
     .                         wpsclrp_sfc )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS ARM specifications
!----------------------------------------------------------------------

        use constants, only: Cp, Lv, grav, sclr_dim

        implicit none

        ! External
        real, external :: diag_ustar

        ! Parameter Constants

        real, parameter ::
     .  ubmin = 0.25

        ! Input variables

        real, intent(in) :: time
        real, intent(in) :: z
        real, intent(in) :: dn0
        real, intent(in) :: thp
        real, intent(in) :: up,vp

        real, intent(inout) :: ustar

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp

!       Output variables (optional)
        real, optional, intent(out) :: wpsclrp_sfc(sclr_dim)

!       Internal variables

        integer i, j
        real usfc, vsfc, ubar
        real true_time
        real heat_flx, moisture_flx
        real heat_flx2, moisture_flx2
        real bflx

!       ARM roughness height
        real, parameter :: z0 = 0.035  ! momentum roughness height

! Compute heat and moisture fluxes from ARM data in (W/m2)

        true_time = time
        call arm_sfcflx( true_time, heat_flx, moisture_flx )

! Compute momentum fluxes

!       Convert heat_flx and moisture_flx to natural units
        heat_flx2     = heat_flx / ( Cp * dn0 )    ! (K m/s)
        moisture_flx2 = moisture_flx / ( Lv * dn0 )! (m/s)

!       Heat flux in units of (m2/s3) (needed by diag_ustar)
        bflx = grav/thp * heat_flx2

!       Sfc winds
        usfc = up
        vsfc = vp
        ubar = max( ubmin, sqrt( usfc ** 2 + vsfc ** 2 ) )

!       Compute ustar
        ustar = diag_ustar( z, bflx, ubar, z0 )

!       Assign fluxes

        uw  = -usfc/ubar * ustar * ustar
        vw  = -vsfc/ubar * ustar * ustar
        wtp = heat_flx2
        wrp = moisture_flx2

        if (present( wpsclrp_sfc )) wpsclrp_sfc(SCLR_THETA) = wtp
        if (present( wpsclrp_sfc )) wpsclrp_sfc(SCLR_RT) = wrp

        return
        end subroutine arm_sfclyr

!------------------------------------------------------------------------
        subroutine arm_sfcflx( time, heat_flx, moisture_flx )

!       Description:
!       This subroutine computes surface heat and moisture for a specific time
!       according to GCSS ARM specifications. Flux returned are in (W/m2)
!------------------------------------------------------------------------
        implicit none

        ! Parameter constants
        integer, parameter :: ntimes = 7

        real, parameter, dimension(ntimes) :: 
     .  times = (/ 41400., 55800., 64800., 68400.,
     .             77400., 86400., 93600. /),
        ! H and LE specifications
     .  H  = (/-30,  90, 140, 140, 100, -10, -10/),
     .  LE = (/  5, 250, 450, 500, 420, 180,   0/)

        ! Input variable
        real, intent(in) :: time !  Current time [s]

        ! Output variables
        real, intent(out) :: heat_flx, moisture_flx

        ! Local variables
        integer :: i1, i2
        real :: a

        if ( time <= times(1) ) then
           heat_flx     = H(1)
           moisture_flx = LE(1)
        else if ( time >= times(ntimes) ) then
           heat_flx     = H(ntimes)
           moisture_flx = LE(ntimes)
        else
           i1 = 1
           do while ( i1 <= ntimes-1 )
              i2 = i1 + 1
              if ( time >= times(i1) .and. time < times(i2) ) then
                 a            = (time-times(i1))/(times(i2)-times(i1))
                 heat_flx     = ( 1. - a ) * H(i1) + a * H(i2)
                 moisture_flx = ( 1. - a ) * LE(i1) + a * LE(i2)
                 i1           = ntimes
              end if
              i1 = i2
           end do
        end if ! time <= times(1)

        return
        end subroutine arm_sfcflx

        end module arm
