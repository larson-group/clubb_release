#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: gcss.F,v 1.25 2006-12-14 21:57:52 mjfalk Exp $ 
        module gcss

!       Description:
!       This file contains extra subroutines to simulate the GCSS ARM, 
!       GCSS BOMEX, GCSS ATEX, and GCSS FIRE intercomparison cases 
!       as well as Wangara day 33.

!       References:
!       <http://www.gewex.org/gcss.html>
!----------------------------------------------------------------------
        use model_flags

        implicit none

        private ! Default to private

        public ::
     .  cloud_rad, ! Used by both FIRE and ATEX
     .  atex_tndcy, 
     .  atex_sfclyr, 
     .  fire_tndcy,
     .  wangara_tndcy, 
     .  wangara_sfclyr,
     .  astex_tndcy, 
     .  astex_sfclyr,
     .  arm_tndcy, 
     .  arm_sfclyr,
     .  bomex_tndcy, 
     .  bomex_sfclyr,
     .  dycoms2_rf01_tndcy, 
     .  dycoms2_rf01_sfclyr,
     .  dycoms2_rf02_tndcy, 
     .  dycoms2_rf02_sfclyr,
     .  nov11_altocu_tndcy, ! Note that neither altocu case 
     .  jun25_altocu_tndcy, !   starts at the surface
     .  mpace_tndcy,        ! Michael Falk, 1 Nov 2006
     .  mpace_sfclyr,       ! Michael Falk, 1 Nov 2006
     .  cobra_tndcy,
     .  cobra_sfclyr,
     .  rico_tndcy,        ! Michael Falk, 13 Dec 2006
     .  rico_sfclyr        ! Michael Falk, 13 Dec 2006

        private ::
     .  diag_ustar, 
     .  arm_sfcflx, 
     .  Diff_denom,
     .  altocu_icedf

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
        subroutine bomex_tndcy( time, wmt, wmm, theta_tndcy, rt_tndcy,
     .                          sclr_tndcy )
!       Description:
!       Subroutine to set theta and water tendencies for BOMEX case

!       References:
!       None
!----------------------------------------------------------------------

        use grid_class
        use constants, only: sclr_dim

        implicit none

        ! Input
        real, intent(in)  :: time ! Model time [s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) ::
     .  wmt,         ! w wind on thermodynamic grid [m/s]
     .  wmm,         ! w wind on momentum grid [m/s]
     .  theta_tndcy, ! Liquid water potential temperature tendency [K/s]
     .  rt_tndcy     ! Total water mixing ratio tendency [kg/kg/s]

        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclr_tndcy ! Passive scalar forcing [units vary]

        ! Local Variables

        integer :: i

!       Large scale subsidence

        do i = 2, gr%nnzp, 1

           if ( gr%zt(i) >= 0. .and. gr%zt(i) < 1500. ) then
              wmt(i) = - ( 0.0065 / 1500. ) * gr%zt(i)
           else if ( gr%zt(i) >= 1500. .and. gr%zt(i) < 2100. ) then
              wmt(i)
     .          = - 0.0065 
     .            + 0.0065 * ( gr%zt(i) - 1500. ) / ( 2100. - 1500. )
           else
              wmt(i) = 0.
           end if

        end do

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
        ! Radiative theta-l tendency

        do i = 2, gr%nnzp

           if ( gr%zt(i) >= 0. .and. gr%zt(i) < 1500. ) then
              theta_tndcy(i) = -2.315e-5
           else if ( gr%zt(i) >= 1500. .and. gr%zt(i) < 2500. ) then
              theta_tndcy(i)
     .          = - 2.315e-5 
     .            + 2.315e-5 
     .              * ( gr%zt(i) - 1500. ) / ( 2500. - 1500. )
           else
              theta_tndcy(i) = 0.
           end if

        end do

!       Large scale advective moisture tendency

        do i = 2, gr%nnzp

           if ( gr%zt(i) >= 0. .and. gr%zt(i) < 300. ) then
              rt_tndcy(i) = -1.2e-8
           else if ( gr%zt(i) >= 300. .and. gr%zt(i) < 500. ) then
              rt_tndcy(i) 
     .          = - 1.2e-8 
     .              * ( 1. - ( gr%zt(i) - 300. )/( 500. - 300. ) )
           else
              rt_tndcy(i) = 0.
           end if

        end do

        ! Boundary conditions
        theta_tndcy(1) = 0.0  ! Below surface
        rt_tndcy(1)    = 0.0  ! Below surface

        if ( present( sclr_tndcy ) ) then 
          sclr_tndcy(1:gr%nnzp,SCLR_THETA) = theta_tndcy
          sclr_tndcy(1:gr%nnzp,SCLR_RT)    = rt_tndcy
        end if

        return
        end subroutine bomex_tndcy

!----------------------------------------------------------------------
        subroutine fire_tndcy( time, wmt, wmm, theta_tndcy, rt_tndcy,
     .                         sclr_tndcy )
!       Description:
!       Subroutine to large-scale subsidence for FIRE case

!       References:
!       None
!----------------------------------------------------------------------

        use constants, only: sclr_dim
        use grid_class

        implicit none

        ! Input
        real, intent(in) :: time ! Model time [s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) ::
     .  wmt,         ! w wind on thermodynamic grid [m/s]
     .  wmm,         ! w wind on momentum grid [m/s]
     .  theta_tndcy, ! Liquid water potential temperature tendency [K/s]
     .  rt_tndcy     ! Total water mixing ratio tendency [kg/kg/s]

        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclr_tndcy ! Passive scalar tendency [units vary]

!       Internal variables

        integer :: i

!       Large-scale subsidence

        do i = 2, gr%nnzp

           if ( gr%zt(i) >= 0. .and. gr%zt(i) < 1500. ) then
              wmt(i) = - 5.e-6 * gr%zt(i)
           end if

        end do

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
!       Radiative theta-l tendency is computed interactively elsewhere

        theta_tndcy = 0.0

!       Large scale advective moisture tendency

        rt_tndcy    = 0.0

        if ( present( sclr_tndcy ) ) then
          sclr_tndcy(1:gr%nnzp,SCLR_THETA) = theta_tndcy(1:gr%nnzp)
          sclr_tndcy(1:gr%nnzp,SCLR_RT)    = rt_tndcy(1:gr%nnzp)
        end if

        return
        end subroutine fire_tndcy

!----------------------------------------------------------------------
        subroutine wangara_tndcy( time, wmt, wmm, theta_tndcy, 
     .                            rt_tndcy,
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
     .  wmt,         ! w wind on thermodynamic grid [m/s]
     .  wmm,         ! w wind on momentum grid [m/s]
     .  theta_tndcy, ! Liquid water potential temperature tendency [K/s]
     .  rt_tndcy     ! Total water mixing ratio tendency [kg/kg/s]

        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclr_tndcy ! Passive scalar tendency [units vary]

!       No large-scale subsidence for now
        wmt = 0.0
        wmm = 0.0

!       No large-scale water tendency or cooling

        rt_tndcy    = 0.0
        theta_tndcy = 0.0

        if ( present( sclr_tndcy ) ) then
          sclr_tndcy(1:gr%nnzp,SCLR_THETA) = theta_tndcy(1:gr%nnzp)
          sclr_tndcy(1:gr%nnzp,SCLR_RT)    = rt_tndcy(1:gr%nnzp)
        end if

        return
        end subroutine wangara_tndcy

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for ATEX case
!
        subroutine atex_tndcy( time, time_initial, rtm, wmt, wmm,
     .                         theta_tndcy, rt_tndcy, isValid,
     .                         sclr_tndcy )

        use constants, only: sclr_dim
        use grid_class

        implicit none

!       Input
        real, intent(in) :: time
        real, intent(in) :: time_initial
        real, intent(in) :: rtm(gr%nnzp)

!       Output
        real, intent(out) :: wmt(gr%nnzp)
        real, intent(out) :: wmm(gr%nnzp)
        real, intent(out) :: theta_tndcy(gr%nnzp)
        real, intent(out) :: rt_tndcy(gr%nnzp)

!       Output (optional)
        real, optional, intent(out), dimension(gr%nnzp, sclr_dim)
     .                    :: sclr_tndcy

!       Misc
        logical, intent(inout) :: isValid

!       Internal variables

        integer i
        real zi

!       Forcings are applied only after t = 5400 s

        wmt         = 0.
        wmm         = 0.
        theta_tndcy = 0.
        rt_tndcy    = 0.

        if ( time >= time_initial + 5400.0 ) then

!          Identify height of 6.5 g/kg moisture level

           i = 2
           do while ( i <= gr%nnzp .and. rtm(i) > 6.5e-3 )
              i = i + 1
           end do
           if ( i == gr%nnzp+1 .or. i == 2 ) then
              write(*,*) 'Identification of 6.5 g/kg level failed'
              write(*,*) 'Subroutine: atex_tndcy. File: gcss.f'
              write(*,*) 'i = ',i
              write(*,*) 'rtm(i) = ',rtm(i)
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
                 theta_tndcy(i) 
     .             = -1.1575e-5 * ( 3. - gr%zt(i)/zi )
              else if ( gr%zt(i) > zi .and. gr%zt(i) <= zi+300. ) then
                 theta_tndcy(i)
     .             = -2.315e-5 * ( 1. - (gr%zt(i)-zi)/300. )
              else
                 theta_tndcy(i) = 0.0
              end if

           end do

           ! Moisture tendency

           do i = 2, gr%nnzp

              if ( gr%zt(i) > 0. .and. gr%zt(i) < zi ) then
                 rt_tndcy(i) = -1.58e-8 * ( 1. - gr%zt(i)/zi )  ! Brian
              else
                 rt_tndcy(i) = 0.0       ! Brian
              end if

           end do

           ! Boundary conditions
           theta_tndcy(1) = 0.0  ! Below surface
           rt_tndcy(1)    = 0.0  ! Below surface

        end if ! time >= time_initial + 5400.0

        if ( present( sclr_tndcy ) ) then
          sclr_tndcy(1:gr%nnzp,SCLR_THETA) = theta_tndcy(1:gr%nnzp)
          sclr_tndcy(1:gr%nnzp,SCLR_RT)    = rt_tndcy(1:gr%nnzp)
        end if

        return
        end subroutine atex_tndcy

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to GCSS ARM specifications
!
        subroutine arm_sfclyr( time, z, dn0, thp, up, vp, uw, vw, 
     .                         wtp, wrp, ustar,
     .                         wpsclrp_sfc )

        use constants, only: Cp, Lv, grav, sclr_dim
        implicit none

        ! External
        ! In module calls:
        ! real diag_ustar

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

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to GCSS BOMEX specifications
!
        subroutine bomex_sfclyr( up, vp, uw, vw, wtp, wrp,
     .                           wpsclrp_sfc )

        use constants, only: sclr_dim

        implicit none

!       Constant

        real, parameter :: ubmin = 0.25
        real, parameter :: ustar = 0.28

!       Input variables

        real, intent(in) :: up, vp

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp

!       Output variables (optional)

        real, optional, intent(out) :: wpsclrp_sfc(sclr_dim)

!       Internal variables

        real ubar

! Compute heat and moisture fluxes

        wtp = 8.e-3
        wrp = 5.2e-5
        
        if ( present(wpsclrp_sfc) ) then 
          wpsclrp_sfc(SCLR_THETA) = 8.e-3
          wpsclrp_sfc(SCLR_RT)    = 5.2e-5
        end if

! Compute momentum fluxes

        ubar = max( ubmin, sqrt(up*up + vp*vp) )

        uw = -up * ustar*ustar / ubar
        vw = -vp * ustar*ustar / ubar

        return
        end subroutine bomex_sfclyr

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to GCSS ATEX specifications
!
        subroutine atex_sfclyr( up, vp, thp, rtp, uw, vw, wtp, wrp,
     .                          sclrm_sfc, wpsclrp_sfc )

        use constants

        implicit none

!       Constants

        real, parameter :: ubmin = 0.25
        real, parameter :: ustar = 0.3
        real, parameter :: C_10 = 0.0013
        real, parameter :: SST = 298

!       Input variables

        real, intent(in) :: up,vp,thp,rtp

!       Input variables (optional)

        real, optional, intent(in) :: sclrm_sfc(sclr_dim)

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp

!       Output variables (optional)

        real, optional, intent(out) :: wpsclrp_sfc(sclr_dim)

!       Internal variables

        real ubar

! Compute heat and moisture fluxes

        ubar = max(ubmin, sqrt(up*up + vp*vp))

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

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture for Wangara day 33
!
        subroutine wangara_sfclyr( time, up, vp, uw, vw, wtp, wrp,
     .                             wpsclrp_sfc )
        use constants, only: sclr_dim, pi

        implicit none

!       Constants

        real, parameter :: ubmin = 0.25
        real, parameter :: ustar = 0.13

!       Input variables

        real, intent(in) :: time
        real, intent(in) :: up, vp

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp
        
        real, optional, intent(out), dimension(sclr_dim) :: 
     .        wpsclrp_sfc

!       Internal variables

        real ubar
        real time_utc, time_est

! Compute UTC time of the day in seconds

        time_utc = mod( time, 86400. )

! Now convert UTC time to Australia EST (local time)

        time_est = mod( time_utc + 36000., 86400. )

        if ( time_est < 27000 .or. time_est > 63000 ) then
           write(*,*)
     .   'wangara_sfclyr: error local time must be between 730 and 1730'
           write(*,*) 'time_est = ',time_est
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

!------------------------------------------------------------------------
! This subroutine computes surface heat and moisture for a specific time
! according to GCSS ARM specifications. Flux returned are in (W/m2)
!

        subroutine arm_sfcflx( time, heat_flx, moisture_flx )
        implicit none

!       Input variable
        real time       ! time in seconds

!       Output variables
        real heat_flx, moisture_flx

!       Internal variables
        integer i1,i2
        real a

        integer ntimes
        parameter ( ntimes = 7 )
        real times(ntimes),H(ntimes), LE(ntimes)
        data times / 41400., 55800., 64800., 68400.,
     +               77400., 86400., 93600. /

!       H and LE specifications
        data H  / -30, 90, 140, 140, 100, -10, -10 /
        data LE / 5, 250, 450, 500, 420, 180, 0 /

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
        end if

        return
        end subroutine arm_sfcflx

! ----------------------------------------------------------------------
!
! DISCLAIMER : this code appears to be correct but has not been
!              very thouroughly tested. If you do notice any
!              anomalous behaviour then please contact Andy and/or
!              Bjorn
!
! Function diag_ustar:  returns value of ustar using the below
! similarity functions and a specified buoyancy flux (bflx) given in
! kinematic units
!
! phi_m (zeta > 0) =  (1 + am * zeta)
! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
!
! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
!
! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
! Layer, in Workshop on Micormeteorology, pages 67-100.
!
! Code writen March, 1999 by Bjorn Stevens
!
      real function diag_ustar( z, bflx, wnd, z0 ) 

      use constants, only: eps, grav, vonk, pi

      implicit none

      real, parameter      :: am   =  4.8   !   "          "         "
      real, parameter      :: bm   = 19.3   !   "          "         "

      real, intent (in)    :: z             ! height where u locates
      real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
      real, intent (in)    :: wnd           ! wind speed at z
      real, intent (in)    :: z0            ! momentum roughness height

      integer :: iterate
      real    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

      lnz   = log( z / z0 )
      klnz  = vonk/lnz
      c1    = pi / 2.0 - 3.0*log( 2.0 )

      ustar =  wnd*klnz
!      if (bflx /= 0.0) then
      if (abs(bflx) > 1.e-6) then
!      if (abs(bflx) > 1.e-4) then
        do iterate=1,4
!          lmo   = -bflx * vonk/(ustar**3 + eps)
          lmo   = -ustar**3 / ( vonk * bflx )
          zeta  = z/lmo
          if (zeta > 0.) then
            ustar =  vonk*wnd  /(lnz + am*zeta)
          else
            x     = sqrt( sqrt( 1.0 - bm*zeta ) )
            psi1  = 2.*log( 1.0+x ) + log( 1.0+x*x ) - 2.*atan( x ) + c1
            ustar = wnd*vonk/(lnz - psi1)
          end if
        end do
      end if

      diag_ustar = ustar

      return
      end function diag_ustar

!----------------------------------------------------------------------
! Subroutine to compute cloud IR radiation using a simple scheme
! based on LWP
!
! Based on GCSS ATEX intercomparison case
!----------------------------------------------------------------------

        subroutine cloud_rad( rhot, rcm, exner, Frad, radht, 
     .                        thlm_forcing, sclr_rate )

        use grid_class
        use constants

        implicit none

!       Input

        real, intent(in) :: rhot(gr%nnzp)       ! density (thermo point)
        real, intent(in) :: rcm(gr%nnzp)        ! liquid water mixing ratio
        real, intent(in) :: exner(gr%nnzp)      ! exner function

!       Output

        real, intent(out) :: Frad(gr%nnzp)          ! IR radiative flux
        real, intent(out) :: radht(gr%nnzp)         ! Radiative heating rate

        real, intent(inout) :: thlm_forcing(gr%nnzp)! Radht + 

!       Output (optional)

        real, optional, intent(out) :: sclr_rate(gr%nnzp)   ! Radiative heating rate ( scalar(1) )

!       Local

        real, dimension(1:gr%nnzp):: LWP        ! Liquid Water Path

        real F0,k
        parameter (F0 = 74., k = 130.)

        integer i, im1, ip1

        if ( .not. lbugsrad ) then
!       Compute liquid water path from top of the model
!       We define liquid water path on momentum levels

          LWP(gr%nnzp) = 0.
          do i = gr%nnzp-1, 1, -1
            LWP(i) = LWP(i+1) + rhot(i+1) * rcm(i+1) / gr%dzt(i+1)
          end do

!         Compute IR radiative flux

          do i = 1, gr%nnzp, 1
            Frad(i) = F0 * EXP( -k * 1.0 * LWP(i) )
          end do

!         Compute IR heating rate

          radht          = ( -1.0/(Cp*rhot) ) 
     .                   * ddzm( Frad ) * ( 1.0 / exner )
          radht(1)       = 0.
          radht(gr%nnzp) = 0.

          ! note that for ATEX after 90 minutes, advect and clear air
          ! radiation must be added in from atex_tndcy
          thlm_forcing(:) = thlm_forcing(:) + radht(:)
          if ( present( sclr_rate ) ) then
            sclr_rate(:) = thlm_forcing(:) ! Just theta
          end if

        end if ! ~ lbugsrad

        return
        end subroutine cloud_rad

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for DYCOMS RF01 case
!
        subroutine dycoms2_rf01_tndcy( time, rhot, rhom, rtm, rcm,
     .                                 exner, wmt, wmm, Frad, radht, 
     .                                 theta_tndcy, rt_tndcy, isValid,
     .                                 sclr, sclr_tndcy )

        use grid_class
        use constants
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

!       Input
        real, intent(in) :: time
        real, intent(in) :: rhot(gr%nnzp)
        real, intent(in) :: rhom(gr%nnzp)
        real, intent(in) :: rtm(gr%nnzp)
        real, intent(in) :: rcm(gr%nnzp)
        real, intent(in) :: exner(gr%nnzp)

        real, optional, intent(in) :: sclr(gr%nnzp, sclr_dim)

!       Output
        real, intent(out) :: wmt(gr%nnzp)
        real, intent(out) :: wmm(gr%nnzp)

        real, intent(out) :: Frad(gr%nnzp)
        real, intent(out) :: radht(gr%nnzp)
        real, intent(out) :: theta_tndcy(gr%nnzp)

        real, intent(out) :: rt_tndcy(gr%nnzp)

        real, optional, 
     .  intent(out) :: sclr_tndcy(gr%nnzp, sclr_dim)

!       Misc
        logical, intent(inout) :: isValid

!       Parameter

        real, parameter :: lsdiv =  3.75e-6
        real, parameter :: F0 = 70.0, F1 = 22.0, kay = 85.0

!       Internal variables

        integer i
        real zi
        real lwp(gr%nnzp)

        wmt         = 0.
        wmm         = 0.
        theta_tndcy = 0.
        rt_tndcy    = 0.

!       Identify height of 8.0 g/kg moisture level

        i = 2
        do while ( i <= gr%nnzp .and. rtm(i) > 8.0e-3 )
           i = i + 1
        end do
        if ( i == gr%nnzp+1 .or. i == 2 ) then
           write(*,*) 'Identification of 8.0 g/kg level failed'
           write(*,*) 'Subroutine: dycoms2_rf01_tndcy. File: gcss.f'
           write(*,*) 'i = ',i
           write(*,*) 'rtm(i) = ',rtm(i)
           isValid = .false.
           return
        end if
        zi = (gr%zt(i)-gr%zt(i-1))/(rtm(i)-rtm(i-1))*(8.0e-3-rtm(i-1))
     .     + gr%zt(i-1)
!        x_sfc(1,izi) = zi
#ifdef STATS
        if ( lstats_samp .and. izi > 0 ) then
          sfc%x(1,izi) = sfc%x(1,izi) + zi
          sfc%n(1,izi) = sfc%n(1,izi) + 1
        end if
#endif /*STATS*/

!       Large scale subsidence

        do i=2,gr%nnzp
           wmt(i) = - lsdiv * gr%zt(i)
        end do

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
!       theta-l radiative tendency

        if ( .not. lbugsrad ) then
!         Compute liquid water path from top of the model
!         We define liquid water path on momentum levels

          lwp(gr%nnzp) = 0.0
          do i = gr%nnzp-1, 1, -1
            lwp(i) = lwp(i+1) + rhot(i+1) * rcm(i+1) / gr%dzt(i+1)
          end do
!         x_sfc(1,ilwp) = lwp(1)

!         Compute IR radiative flux

          do i = 1, gr%nnzp, 1
            Frad(i) = F0 * EXP( -kay * lwp(i) )
     .              + F1 * EXP( -kay * (lwp(1)-lwp(i)) )
            if ( zi > 0 .and. gr%zm(i) > zi ) then
              Frad(i) = Frad(i)
     .                + rhom(i) * cp * lsdiv
     .                  * ( 0.25*(gr%zm(i)-zi)**(4.0/3.0)
     .                      + zi*(gr%zm(i)-zi)**(1.0/3.0) )
             end if
          end do

!         Compute IR heating rate

          radht          = ( -1.0/(Cp*rhot) ) * ddzm( Frad )
     .                   * ( 1.0 / exner )
          radht(1)       = 0.
          radht(gr%nnzp) = 0.

        end if ! ~ lbugsrad

!       Add heating rate to theta-l forcing

        if ( .not. lbugsrad ) theta_tndcy = theta_tndcy + radht
        
        if (present( sclr_tndcy )) sclr_tndcy(:,SCLR_THETA)= theta_tndcy
        if (present( sclr_tndcy )) sclr_tndcy(:,SCLR_RT)   = rt_tndcy

        return
        end subroutine dycoms2_rf01_tndcy

!----------------------------------------------------------------------

        SUBROUTINE dycoms2_rf02_tndcy
     .             ( time, time_initial, rhot, rhom, rtm, rcm, exner, 
     .               wmt, wmm, theta_tndcy, rt_tndcy, 
     .               Frad, radht, Ncm, Ncnm, isValid, 
     .               sclr_tndcy )

        USE grid_class
        USE constants
#ifdef STATS
        USE statistics
#endif /*STATS*/

        ! Input Variables
        real, intent(in) ::
     .  time,         ! Current time    [s]
     .  time_initial  ! Initial time    [s]

        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rhot
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rhom
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rtm
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rcm
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: exner

        ! Output
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: wmt
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: wmm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: theta_tndcy
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: rt_tndcy
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Frad
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: radht
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Ncnm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Ncm

        ! Optional output
        REAL, OPTIONAL, 
     .  DIMENSION(gr%nnzp, sclr_dim), INTENT(OUT):: sclr_tndcy

        ! Misc
        LOGICAL, INTENT(INOUT) :: isValid

        ! Local
        REAL, DIMENSION(1:gr%nnzp):: LWP
        REAL, DIMENSION(1:gr%nnzp):: Heaviside

        REAL:: z_i

        REAL, PARAMETER :: ls_div = 3.75e-6
        REAL, PARAMETER :: kap    = 85.0  ! m^2/kg
        REAL, PARAMETER :: F0     = 70.0  ! W/m^2
        REAL, PARAMETER :: F1     = 22.0  ! W/m^2

        INTEGER:: k

        ! Large-scale subsidence

        DO k = 2, gr%nnzp, 1
          wmt(k) = -ls_div * gr%zt(k)
        END DO

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
        IF ( .not. lbugsrad ) THEN
        ! Radiation

!         Compute liquid water path from top of the model
!         We define liquid water path on momentum levels

          LWP(gr%nnzp) = 0.0
          DO k = gr%nnzp-1, 1, -1
            LWP(k) = LWP(k+1) + rhot(k+1)*rcm(k+1)/gr%dzt(k+1)
          END DO

!         Find the height of the isotherm rtm = 8.0 g/kg.

          k = 2
          DO WHILE ( k <= gr%nnzp .AND. rtm(k) > 8.0e-3 )
            k = k + 1
          END DO
          IF ( k == gr%nnzp+1 .or. k == 2 ) THEN
            write(*,*) 'Identification of 8.0 g/kg level failed'
            write(*,*) 'Subroutine: dycoms2_rf02_tndcy. File: gcss.f'
            write(*,*) 'k = ', k
            write(*,*) 'rtm(k) = ', rtm(k)
            isValid = .false.
            return
          END IF
          z_i = ( (gr%zt(k)-gr%zt(k-1))/(rtm(k)-rtm(k-1)) )
     .        * (8.0e-3-rtm(k-1)) + gr%zt(k-1)

!         Compute the Heaviside step function for z - z_i.

          DO k = 1, gr%nnzp, 1
            IF ( gr%zm(k) - z_i  <  0.0 ) THEN
              Heaviside(k) = 0.0
            ELSE IF ( gr%zm(k) - z_i  ==  0.0 ) THEN
              Heaviside(k) = 0.5
            ELSE IF ( gr%zm(k) - z_i  >  0.0 ) THEN
              Heaviside(k) = 1.0
            END IF
          END DO

!         Compute radiative flux profile (Frad).
!         Radiative flux is defined on momentum levels.

          DO k = 1, gr%nnzp, 1

            Frad(k) = F0 * EXP( -kap * LWP(k) )
     .              + F1 * EXP( -kap * (LWP(1) - LWP(k)) )

            IF ( Heaviside(k) > 0.0 ) THEN
              Frad(k) = Frad(k)
     .                + rhom(k) * Cp * ls_div * Heaviside(k)
     .                  * ( 0.25 * ((gr%zm(k)-z_i)**(4.0/3.0))
     .                + z_i * ((gr%zm(k)-z_i)**(1.0/3.0)) )
            END IF

          END DO 

        ! Compute the radiative heating rate.
        ! The radiative heating rate is defined on thermodynamic levels.

          DO k = 2, gr%nnzp, 1
            radht(k) = ( 1.0 / exner(k) ) * ( -1.0/(Cp*rhot(k)) )
     .               * ( Frad(k) - Frad(k-1) ) * gr%dzt(k)
          END DO
          radht(1) = radht(2)
        END IF ! ~ lbugsrad

        ! Enter the final theta-l and rtm tendencies

        IF ( .not. lbugsrad ) THEN
          theta_tndcy(1:gr%nnzp) = radht(1:gr%nnzp)
        ELSE
          theta_tndcy(1:gr%nnzp) = 0.0
        END IF 
        rt_tndcy(1:gr%nnzp) = 0.0

        ! Update surface statistics
#ifdef STATS
        if ( lstats_samp .and. izi > 0 ) then
          sfc%x(1,izi) = sfc%x(1,izi) + z_i
          sfc%n(1,izi) = sfc%n(1,izi) + 1
        end if
#endif /*STATS*/

        ! The following lines of code specify cloud droplet
        ! concentration (Ncm).  The cloud droplet concentration has
        ! been moved here instead of being stated in Subroutine rain
        ! for the following reasons:
        !    a) The effects of cloud droplet sedimentation can be computed
        !       without having to call the precipitation scheme.
        !    b) Ncm tends to be a case-specific parameter.  Therefore, it
        !       is appropriate to declare in the same place as other
        !       case-specific parameters.
        !
        ! Someday, we could move the setting of Ncm to pdf_closure_new
        ! for the following reasons:
        !    a) The cloud water mixing ratio (rcm) is computed using the
        !       PDF scheme.  Perhaps someday Ncm can also be computed by
        !       the same scheme.
        !    b) It seems more appropriate to declare Ncm in the same place
        !       where rcm is computed.
        !
        ! Since cloud base (zb) is determined by the mixing ratio 1.0E-6,
        ! so will cloud droplet number concentration (Ncm).

        if ( lcoamps_micro .and. time == time_initial ) then
          ! Taken from COAMPS subroutine ncn_init()
          Ncnm(1:gr%nnzp) = 55000000.0 / rhot(1:gr%nnzp)

        else
          ! K & K or no micro
          DO k = 1, gr%nnzp, 1
            IF ( rcm(k) > 1.0E-6 ) THEN
              Ncm(k) = 55000000.0
            ELSE
              Ncm(k) = 0.0
            END IF
          END DO

        end if

        if ( present( sclr_tndcy ) ) 
     .    sclr_tndcy(:,SCLR_THETA) = theta_tndcy 
        if ( present( sclr_tndcy ) ) sclr_tndcy(:,SCLR_RT) = rt_tndcy 

        RETURN
        END SUBROUTINE dycoms2_rf02_tndcy

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to GCSS DYCOMS specifications
!
        subroutine dycoms2_rf01_sfclyr( sfctype, Tsfc, psfc, exnersfc,
     .                                  up, vp, thp, rtp, rhop,
     .                                  uw, vw, wtp, wrp,
     .                                  sclrm_sfc, wpsclrp_sfc )

        use constants
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

!       External

        real, external :: rsat

!       Constants

        real, parameter :: ubmin = 0.25
        real, parameter :: ustar = 0.25
        real, parameter :: Cd    = 0.0011

!       Input variables

        integer, intent(in) :: sfctype
        real, intent(in)    :: Tsfc, psfc, exnersfc, up, vp
        real, intent(in)    :: thp, rtp, rhop

        real, optional, intent(in) :: sclrm_sfc(sclr_dim) 

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp

        real, optional, intent(out) :: wpsclrp_sfc(sclr_dim) 

!       Internal variables

        real ubar

! Compute heat and moisture fluxes

        ubar = max(ubmin,sqrt(up*up + vp*vp))

! Compute momentum fluxes

        uw = -up * ustar*ustar / ubar
        vw = -vp * ustar*ustar / ubar

        if ( sfctype == 0 ) then

          wtp =  15.0 / ( rhop * cp )
          wrp = 115.0 / ( rhop * lv )

          if ( present( sclrm_sfc ).and.present( wpsclrp_sfc ) ) then
            wpsclrp_sfc(SCLR_THETA) = wtp
            wpsclrp_sfc(SCLR_RT)    = wrp
          end if
        else if ( sfctype == 1 ) then

          wtp = -Cd * ubar * ( thp - Tsfc/exnersfc )
          wrp = -Cd * ubar * ( rtp - rsat(psfc,Tsfc) )        

          if ( present( sclrm_sfc ).and.present( wpsclrp_sfc ) ) then
            wpsclrp_sfc(SCLR_THETA) = -Cd * ubar * 
     .                        ( sclrm_sfc(SCLR_THETA) - Tsfc/exnersfc )
            wpsclrp_sfc(SCLR_RT)    = -Cd * ubar * 
     .                        ( sclrm_sfc(SCLR_RT) - rsat(psfc, Tsfc) )
          end if
        else

          write(*,*) 'Invalid sfctype value = ', sfctype
          stop

        end if
#ifdef STATS
! Save ustar for statistics
!        x_sfc(1,iustar) = ustar
        if ( lstats_samp .and. iustar > 0 ) then
          sfc%x(1,iustar) = sfc%x(1,iustar) + ustar
          sfc%n(1,iustar) = sfc%n(1,iustar) + 1
        end if
#endif /*STATS*/

        return
        end subroutine dycoms2_rf01_sfclyr

!----------------------------------------------------------------------

        SUBROUTINE dycoms2_rf02_sfclyr( u_sfc, v_sfc, uw, vw, wtp, wrp,
     .                                  ustar, sclr_sfc, wpsclr_sfc )

        USE constants

! Input
        REAL, INTENT(IN):: u_sfc
        REAL, INTENT(IN):: v_sfc

! Input (optional)
        REAL, OPTIONAL, INTENT(IN):: sclr_sfc(sclr_dim)

! Output
        REAL, INTENT(OUT):: uw
        REAL, INTENT(OUT):: vw
        REAL, INTENT(OUT):: wtp
        REAL, INTENT(OUT):: wrp
        REAL, INTENT(OUT):: ustar

! Output (optional)
        REAL, OPTIONAL, INTENT(OUT):: wpsclr_sfc(sclr_dim)

! Internal
        REAL:: wind_sfc

        REAL, PARAMETER:: SH = 16.0
        REAL, PARAMETER:: LH = 93.0

        ustar = 0.25

        wind_sfc = SQRT( u_sfc**2.0 + v_sfc**2.0 )

        IF (wind_sfc > 0.0) THEN

           uw = -u_sfc * (ustar**2.0) / wind_sfc
           vw = -v_sfc * (ustar**2.0) / wind_sfc

        ELSE

           uw = 0.0
           vw = 0.0

        END IF

        wtp = SH / (1.21 * Cp)
        wrp = LH / (1.21 * Lv)

        if (present( wpsclr_sfc )) then
          wpsclr_sfc(SCLR_THETA) = SH / (1.21 * Cp)

          wpsclr_sfc(SCLR_RT) = LH / (1.21 * Lv)
        end if

        RETURN
        END SUBROUTINE dycoms2_rf02_sfclyr

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for ASTEX KK case
!
        subroutine astex_tndcy( time, wmt, wmm, theta_tndcy, rt_tndcy,
     .                          sclr_tndcy )

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
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to ASTEX with Khairoutdinov and Kogan
! alteration.

        subroutine astex_sfclyr( rho0, up, vp, uw, vw, wtp, wrp,
     .                           wpsclrp_sfc )

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

!-----------------------------------------------------------------------

        subroutine nov11_altocu_tndcy
     .             ( time, time_initial, dt, thlm,
     .               rcm, p, exner, rhot, rtm, wmt,
     .               wmm, thl_tndcy, rt_tndcy,
     .               Frad, radht, sclr_tndcy )

        use grid_class
        use constants
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c INPUT VARIABLES                                                             c
c time         :  time of simulation since start                 Unit: s      c
c time_initial :  initial time of simulation                     Unit: s      c
c dt           :  length of model time step                      Unit: s      c
c thlm         :  profile of liquid water potential temperature  Unit: K      c
c rcm          :  profile of liquid (cloud) water mixing ratio   Unit: kg/kg  c
c p            :  basic state pressure                           Unit: hPa    c
c exner        :  basic state Exner function                     Unit: NONE   c
c rhot         :  density of reference state at thermodynamic levels          c
c                                                                Unit: kg/m^3 c
c                                                                             c
c OUTPUT VARIABLES                                                            c
c rtm          :  profile of total water mixing ratio            Unit: kg/kg  c
c wmt          :  profile of vertical ascent/descent at thermodynamic levels  c
c                                                                Unit: m/s    c
c wmm          :  profile of vertical ascent/descent at momentum levels       c
c                                                                Unit: m/s    c
c                                                                             c
c thl_tndcy    :  change in liq. water potential temperature due to           c
c                 radiative heating and ice diffusion            Unit: K      c
c rt_tndcy     :  change in total water due to ice diffusion     Unit: kg/kg  c
c                                                                             c
c Frad         :  profile of total radiative flux (lw + sw)      Unit: W/m^2  c
c                                                                             c
c radht        :  profile of total radiative heating (lw + sw)   Unit: K/s    c
c                                                                             c
c sclr_tndcy   :                                                 Unit: ??     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !!!!!!!!!!!!!!!!!!!!!
        !! Input variables !!
        !!!!!!!!!!!!!!!!!!!!!
        REAL, INTENT(IN):: time
        REAL, INTENT(IN):: time_initial
        REAL, INTENT(IN):: dt
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: thlm
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rcm
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: p
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: exner
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rhot

        !!!!!!!!!!!!!!!!!!!!!!
        !! Output variables !!
        !!!!!!!!!!!!!!!!!!!!!!

        REAL, DIMENSION(1:gr%nnzp), INTENT(INOUT):: rtm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: wmt
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: wmm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: thl_tndcy
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: rt_tndcy

        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Frad

        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: radht

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Output variables (optional) !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, OPTIONAL, DIMENSION(gr%nnzp, sclr_dim), 
     .                  INTENT(OUT) :: sclr_tndcy

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c LOCAL VARIABLES                                                             c
c dqc_dt_icedfs : change in liquid water over time, due to ice diffusional    c
c                 growth.                                       Unit: 1/s     c
c                                                                             c
c LWP           : liquid water path from domain top.            Unit: kg/m^2  c
c                                                                             c
c TEMPORARY ARRAYS USED FOR COAMPS RADIATIVE SCHEME                           c
c (SEE COMMENTS BELOW FOR SCHEME DESCRIPTION)                                 c
c thm_rad       : temporary "flipped" array of potential temperature          c
c                                                               Unit: K       c
c rcm_rad       : temporary "flipped" array of liquid water mixing ratio      c
c                                                               Unit: kg/kg   c
c rhot_rad      : temporary "flipped" array of air density      Unit: kg/m^3  c
c dsigm         : temporary "flipped" array of grid spacing     Unit: m       c
c                                                                             c
c frad_out      : temporary "flipped" array of total radiative flux           c
c                 from rad_lwsw                                 Unit: W/m^2   c
c frad_lw_out   : temporary "flipped" array of LW radiative flux              c
c                 from rad_lwsw                                 Unit: W/m^2   c
c frad_sw_out   : temporary "flipped" array of SW radiative flux              c
c                 from rad_lwsw                                 Unit: W/m^2   c
c                                                                             c
c radhtk        : temporary "flipped" array of total radiative heating        c
c                 from rad_lwsw                                 Unit: K/s     c
c radht_lw_out  : temporary "flipped" array of LW radiative heating           c
c                 from rad_lwsw                                 Unit: K/s     c
c radht_sw_out  : temporary "flipped" array of SW radiative heating           c
c                 from rad_lwsw                                 Unit: K/s     c
c                                                                             c
c INTERPOLATION ARRAYS AND CONSTANTS                                          c
c zsubs         : heights at which wmt data is supplied                       c
c                 (used for subsidence interpolation)           Unit: m       c
c tsubs         : times after initialization at which wmt data is supplied    c
c                 (NOT USED IN NOV.11 CASE)                     Unit: s       c
c wtX(Y)        : vertical velocity specified at height Y and time X          c
c                 (ONLY wt1 IS USED IN NOV.11 CASE)             Unit: m/s     c
c w1-w2         : vertical velocity before (w1) and after (w2)                c
c                 the current time at the specified level                     c
c                 (NOT USED IN NOV.11 CASE)                     Unit: m/s     c
c                                                                             c
c ADDITIONAL PARAMETERS FOR NOV.11 SUBSIDENCE                                 c
c subs_on       : logical variable tells us whether to turn subsidence on     c
c                                                               Unit: NONE    c
c wmax          : defines value of maximum subsidence in profile              c
c                                                               Unit: cm/s    c
c zi            : defines approx. height of inversion within cloud            c
c                 (subsidence is equal to wmax at this height)  Unit: m       c
c dazi          : defines height above inversion                              c
c                 (above this height, subsidence linearly tapers off to zero) c
c                                                               Unit: m       c
c dbzi          : defines height above inversion                              c
c                 (below this height, subsidence linearly tapers off to zero) c
c                                                               Unit: m       c
c dac           : defines height above cloud                                  c
c                 (at / above this height, we have NO subsidence)             c
c                                                               Unit: m       c
c dbc           : defines height below cloud                                  c
c                 (at / below this height, we have NO subsidence)             c
c                                                               Unit: m       c
c                                                                             c
c RADIATION PARAMETERS                                                        c
c sw_on         : logical variable passed to radiation scheme                 c
c                 - is SW radiation on?                         Unit: NONE    c
c lw_on         : logical variable passed to radiation scheme                 c
c                 - is LW radiation on?                         Unit: NONE    c
c center        : use centered differencing (as opposed to a one-sided        c
c                 forward difference) in radiation code         Unit: NONE    c
c                                                                             c
c xi_abs        : cosine of the solar zenith angle              Unit: NONE    c
c F0            : coefficient for cloud top heating (see Stevens)             c
c                                                               Unit: W/m^2   c
c F1            : coefficient for cloud base heating (see Stevens)            c
c                                                               Unit: W/m^2   c
c kap           : "a constant" according to Duynkerke eqn. 5, where his       c
c                 value is 130 m^2/kg.                          Unit: m^2/kg  c
c radius        : effective droplet radius                      Unit: m       c
c AA            : albedo -- sea surface, according to Lenderink.              c
c                                                               Unit: NONE    c
c gc            : asymmetry parameter, "g" in Duynkerke.        Unit: NONE    c
c Fs0           : The incident incoming SW insolation at cloud top in the     c
c                 direction of the incoming beam (not the vertical).          c
c                                                               Unit: W/m^2   c
c omega         : single-scattering albedo                      Unit: NONE    c
c                                                                             c
c SOLAR ZENITH ANGLE PARAMETERS (NOT USED IN NOV.11 CASE)                     c
c c0            : coefficient for calculation of declination angle from       c
c                 Liou Table 2.2 and Eqn. 2.2.10                Unit: NONE    c
c c1            : same as above                                 Unit: NONE    c
c c2            : same as above                                 Unit: NONE    c
c c3            : same as above                                 Unit: NONE    c
c d1            : Same as above                                 Unit: NONE    c
c d2            : Same as above                                 Unit: NONE    c
c d3            : Same as above                                 Unit: NONE    c
c                                                                             c
c sda_t         : Linear function of day of the year.                         c
c                 sda_t=0 January 1 and sda_t -> 2*pi December 31.            c
c sda_delta     : Solar declination angle function from Liou 2.2.10           c
c sda_h         : Hour angle (Angle through which the earth must to turn      c
c                 to put sun directly overhead on a point's meridian)         c
c                 (Angle between current time and solar noon)                 c
c t_since_noon  : Number of seconds since noon (after noon > 0)               c
c julday        : Julian day of the year (January 1=1; December 31=365)       c
c                                                                             c
c start_time_until_noon: number of seconds between start time and solar noon  c
c                                                               Unit: s       c
c                                                                             c
c Fs0 INTERPOLATION PARAMETERS                                                c
c nparam        : Number of Fs0 values in the list.                           c
c xilist        : Values of cosine of solar zenith angle corresponding to     c
c                 the values in Fslist                                        c
c Fslist        : Values of Fs0 corresponding to the values in xilist.        c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Local ice diffusional arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, DIMENSION(1:gr%nnzp):: dqc_dt_icedfs

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Local radiation arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, DIMENSION(1:gr%nnzp) :: 
     .  Frad_LW, ! Long wave radiative flux     [W/m^2]
     .  Frad_SW, ! Short wave radiative flux    [W/m^2]
     .  radht_LW,! Long wave radiative heating  [K/s]
     .  radht_SW ! Short wave radiative heating [K/s]
        REAL, DIMENSION(1:gr%nnzp):: LWP
        REAL, DIMENSION(1:gr%nnzp):: thm_rad
        REAL, DIMENSION(1:gr%nnzp):: rcm_rad
        REAL, DIMENSION(1:gr%nnzp):: rhot_rad
        REAL, DIMENSION(1:gr%nnzp):: dsigm

        REAL, DIMENSION(1:gr%nnzp):: frad_out
        REAL, DIMENSION(1:gr%nnzp):: frad_lw_out
        REAL, DIMENSION(1:gr%nnzp):: frad_sw_out

        REAL, DIMENSION(1:gr%nnzp):: radhtk
        REAL, DIMENSION(1:gr%nnzp):: radht_lw_out
        REAL, DIMENSION(1:gr%nnzp):: radht_sw_out

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Working arrays for subsidence interpolation !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real zsubs(7)
        real wt1(7)       ! ONLY wt1 IS NEEDED FOR NOV.11 CASE

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Subsidence constant and variables (for Nov.11 case only) !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL:: wmax, zi, dazi, dbzi, dbc, dac

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! LW Radiative constants !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, PARAMETER:: F0   = 104.0  ! W/m^2
        REAL, PARAMETER:: F1   = 62.0   ! W/m^2
        REAL, PARAMETER:: kap  = 94.2 ! m^2/kg

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Working arrays for SW radiation interpolation !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER, PARAMETER :: nparam = 2
        REAL, DIMENSION(1:nparam):: xilist, Fslist

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! SW Radiative constants !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, PARAMETER:: radius = 1.0e-5
        REAL, PARAMETER:: A      = 0.1
        REAL, PARAMETER:: gc     = 0.86
        REAL, PARAMETER:: omega  = .9965

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Additional SW radiative variables !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL:: xi_abs, Fs0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variables used for calculating solar zenith angle                 !!
        !! (Currently not used in Nov.11 case)                               !!
        !!                                                                   !!
        !! NOTE BY ADAM SMITH, 13 April 2006                                 !!
        !! In COAMPS, these constants are called c0, c1, c2, c3, d1, d2,     !!
        !! and d3.  Since HOC uses these variables elsewhere, I have renamed !!
        !! them here.                                                        !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL c_0, c_1, c_2, c_3, d_1, d_2, d_3
        REAL sda_t, sda_delta, sda_h, t_since_noon, julday
        REAL start_time_until_noon

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Toggle for implementing differencing method in interpolations !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LOGICAL :: center 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Toggles for activating/deactivating forcings                  !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LOGICAL :: icedfs_on, subs_on, lw_on, sw_on

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variable used for working within vertical arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER:: k

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! END OF VARIABLE DECLARATION !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Toggles for activating/deactivating forcings                            cc
cc To turn off a specific forcing, set the corresponding toggle to .FALSE. cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         icedfs_on = .TRUE.
         subs_on   = .TRUE.
         lw_on     = .TRUE.
         sw_on     = .TRUE.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Toggle for centered/forward differencing (in interpolations) cc
cc To use centered differencing, set the toggle to .TRUE.       cc
cc To use forward differencing, set the toggle to .FALSE.       cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         center    = .TRUE.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Calculation of Solar Declination Angle                               c
c                                                                      c
c This section of code calculates the solar declination angle based    c
c on the simulaiton's time, date, and latitude.  These parameters are  c
c applied in the scheme, and a value of xi_abs is calculated.          c
c                                                                      c
c If xi_abs <= 0.0 at any time step, we shut off shortwave radiation   c
c for that time step (by setting sw_on = .FALSE.).  If xi_abs > 0.0,   c
c we reactivate solar radiation by setting sw_on = .TRUE.              c
c                                                                      c
c This version is based on solar_declin.m, written by Dr. Vince Larson c
c and based on Liou 2.2.10.                                            c
c                                                                      c
c FOR NOV.11 CASE                                                      c
c ---------------                                                      c
c The Nov.11 case uses a constant value for the cosine of the solar    c
c zenith angle.  Because the simulation is only 4 hours long, and it   c
c starts at 18Z (11am local time), we have rather constant sunlight    c
c through the entire simulation period.  As a result, we will use a    c
c constant value of xi_abs = 0.4329.  For now, the following           c
c calculation of xi_abs will be commented out, then followed by the    c
c manual declaration of xi_abs as a constant.  We may later decide to  c
c use the code later to calculate a more accurate xi_abs.              c
c                                                                      c
c Comment by Adam Smith (ajsmith4) on 26 June 2006                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      julday = 315
c      start_time_until_noon = 7916
c      t_since_noon   = time - start_time_until_noon
c      sda_t = 2*pi*(julday-1)/365
c
c      c0 = 0.006918
c      c1 = -0.399912
c      c2 = -0.006758
c      c3 = -0.002697
c      d1 = 0.070257
c      d2 = 0.000907
c      d3 = 0.000148
c
c      sda_delta = c0 + c1*cos(sda_t) + d1*sin(sda_t) +
c     .            c2*cos(2*sda_t) + d2*sin(2*sda_t) +
c     .            c3*cos(3*sda_t) + d3*sin(3*sda_t)
c
c      sda_h = 2*pi*t_since_noon/86400
c
c      xi_abs = sin(rlat*pi/180) * sin(sda_delta) +
c     .         cos(rlat*pi/180) * cos(sda_delta) * cos(sda_h)
c
c      xi_abs = max(xi_abs,0.)

       xi_abs = .4329

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Modification by Adam Smith 26 June 2006                                   !!
!! It is difficult to remember to set xi_abs = 0 when we want to shut off    !!
!! solar radiation.  If sw_on = .FALSE. above, we will automatically set     !!
!! xi_abs to 0 to avoid confusion or errors.                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (.not. sw_on) then
        xi_abs = 0.
      end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End of ajsmith4's Modification !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (xi_abs .eq. 0.) then
        sw_on = .FALSE.
      else
        sw_on = .TRUE.
      end if


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Fs0 Interpolation Parameters-- these also from Kurt Kotenberg's      c
c BUGSrad output.  Fs0 changes somewhat over the range of solar zenith c
c angles, and we obtained these values by solving                      c
c Fs0 = F_vertical / xi_abs .                                          c
c                                                                      c
c The linearinterpolation function returns Fs0.                        c
c                                                                      c
c FOR NOV.11 CASE                                                      c
c ---------------                                                      c
c As explained above, the solar declination angle is presumed to be    c
c constant for Nov.11 cases.  Because of this, we also use a constant  c
c value of Fs0 for solar radiation.                                    c
c                                                                      c
c Because the Jun.25 case uses a linearinterpolation to calculate      c
c Fs0, we will duplicate that scheme here to keep all of our group's   c
c code consistent.  This means that Fslist will be a 2D array with the c
c same Fs0 at all heights.  Therefore, Fs0 will be constant over the   c
c entire model run.                                                    c
c                                                                      c
c Comment by Adam Smith on 26 June 2006                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      xilist(1) = 0.
      xilist(2) = 1.

      Fslist(1) = 1212.75
      Fslist(2) = 1212.75

      call linearinterpolation( nparam, xilist, Fslist, xi_abs, Fs0 )


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Subsidence Parameters                                                c
c                                                                      c
c FOR NOV.11 CASE                                                      c
c ---------------                                                      c
c The Nov.11 case uses a constant subsidence profile, initiated after  c
c 1 hour of model runtime.  This initial hour is used to "spinup" the  c
c simulation and produce a realistic cloud.                            c
c                                                                      c
c In jun25.F (in mjfalk's /coamps/mod/consolidated5 directory on       c
c condella), subsidence varies over time.  As a result, he uses a      c
c number of arrays defining subsidence profiles for different times.   c
c Since Nov.11's subsidence does not vary with time, we only need one  c
c of these arrays, and the rest have been removed.  The array listing  c
c different times for subsidence has been removed as well.             c
c                                                                      c
c Comment by Adam Smith on 26 June 2006                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subsidence Parameters !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        wmax =  -0.03
          zi = 1400.0
        dazi =  400.0
        dbzi = 1000.0
         dbc =  300.0
         dac =  200.0

        zsubs(1) = 0.
        zsubs(2) = zi-dbzi-dbc
        zsubs(3) = zi-dbzi
        zsubs(4) = zi
        zsubs(5) = zi+dazi
        zsubs(6) = zi+dazi+dac
        zsubs(7) = 2400.

          wt1(1) = 0.
          wt1(2) = 0.
          wt1(3) = wmax
          wt1(4) = wmax
          wt1(5) = wmax
          wt1(6) = 0.
          wt1(7) = 0.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!! SPECIAL NOV.11 CONDITION FOR TOTAL WATER ABOVE CLOUD          !!
!! One hour after the initial time, the total water above cloud  !!
!! is adjusted to be 0.89 of what it previously was.             !!
!!                                                               !!
!! The conditional statement here is set so that if the timestep !!
!! is such that there is no timestep at exactly 3600.0 seconds,  !!
!! then the operation still happnens at the first timestep and   !!
!! only the first timestep after 3600.0 seconds.                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF ( time >= time_initial + 3600.0    .AND.
     .       time <  time_initial + 3600.0 + dt ) THEN

           DO k = 1, gr%nnzp, 1
              IF ( gr%zt(k) > 1555.0 ) THEN
                 rtm(k) = 0.89 * rtm(k)
              END IF
           ENDDO

        END IF


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c SPECIAL METHOD USED TO CALCULATE RADIATION                           c
c Grid descriptions by Adam Smith, 27 June 2006                        c
c                                                                      c
c In order to verify our HOC simulations are working properly, we      c
c have first developed a series of 3D simulations using the COAMPS-LES c
c model.  This large-eddy simulation (LES) simulation uses specific    c
c methods to calculate radiation, subsidence, and other microphysical  c
c processes.  To make the two models simluate clouds as closely as     c
c possible, we use the same radiation scheme in both models.           c
c                                                                      c
c In COAMPS-LES, we use a separate subroutine, rad_lwsw, to implement  c
c all radiation code.  This allows the subroutine to be duplicated     c
c exactly in many different models.  However, the subroutine uses the  c
c COAMPS vertical grid.  Therefore, for HOC to implement this code     c
c correctly, we must modify some of our variable profiles before       c
c calling the radiation subroutine.                                    c
c                                                                      c
c The following diagram describes the differences in model grids:      c
c                                                                      c
c       COAMPS-LES                                    HOC              c
c                                                                      c
c  k= 1   (m) ----------    <MODEL TOP>    ---------- k=kk+1      (m)  c
c                                                                      c
c  k= 1   (t) ----------                   ---------- k=kk+1      (t)  c
c                                                                      c
c  k= 2   (m) ----------                   ---------- k=kk        (m)  c
c                                                                      c
c  k= 2   (t) ----------                   ---------- k=kk        (t)  c
c                                                                      c
c            .                  .                     .                c
c            .                  .                     .                c
c            .                  .                     .                c
c                                                                      c
c  k=kk-1 (m) ----------  m = mass levels  ---------- k=3         (m)  c
c                                                                      c
c  k=kk-1 (t) ----------  t = thermo       ---------- k=3         (t)  c
c                                  levels                              c
c  k=kk   (m) ----------                   ---------- k=2         (m)  c
c                                                                      c
c  k=kk   (t) ----------  kk = number of   ---------- k=2         (t)  c
c                              vertical                                c
c  k=kk+1 (m) ----------       heights     ---------- k=1         (m)  c
c                                                                      c
c //////////////////////// MODEL SURFACE ///////////////////////////// c
c                                          ---------- k=1         (t)  c
c                                                                      c
c                                                                      c
c The major difference in the grids is that HOC uses an additional     c
c thermodynamic level below the model "surface".  This means that all  c
c HOC thermodynamic heights are shifted down one vertical level, and   c
c HOC also has one fewer momentum level than COAMPS.  Therefore, we    c
c use one additional vertical level in HOC, to make sure that the      c
c vertical domain matches in both models.                              c
c                                                                      c
c Also, the COAMPS grid indices are numbered from the top of the model c
c downward, while the HOC grid indices are numbered from the bottom    c
c up.  Therefore, since we are using a COAMPS radiation scheme, we     c
c flip moisture and temperature profiles that are passed into the      c
c rad_lwsw subroutine.  The rad scheme will produce results in using   c
c the COAMPS grid scheme, so all radiation output will be flipped      c
c back to the HOC grid before being applied to the model.              c
c                                                                      c
c Finally, since the COAMPS scheme does not have a gridpoint below     c
c model surface, we add that point to all radiative output files once  c
c they are converted back to HOC setup.  This allows all averages and  c
c calculations to be done correctly.                                   c
c                                                                      c
c                                                                      c
c Computation of radiative fluxes on staggered grid                    c
c Comments by Michael Falk, 16 February 2005.                          c
c                                                                      c
c Frad (and its components Frad_LW and Frad_SW) should be computed on  c
c w points, not on mass points, which is apparent from its formulation c
c and from its location in stats_sw instead of stats_sm.  The grid     c
c looks like this:                                                     c
c                                                                      c
c                                                                      c
c -----Frad----------------------------------    k = 1  (w level)      c
c     /    \            |-dwm                                          c
c -LWP------radht----------------------------    k = 1  (mass level)   c
c     \    /            |                                        c
c -----Frad----------------------------------    k = 2  (w level)      c
c     /    \                                                           c
c -LWP------radht----------------------------    k = 2  (mass level)   c
c     \    /                                                           c
c -----Frad----------------------------------    k = 3  (w level)      c
c     /    \                                                           c
c -LWP------radht----------------------------    k = 3  (mass level)   c
c                                                                      c
c If you consider Frad to take place on mass levels, then computing    c
c LWP is a forward difference and is only first-order accurate, while  c
c if Frad computed in between LWP levels, it is a centered difference  c
c which is second-order accurate.                                      c
c                                                                      c
c The coding implementation requires that Frad depend on LWP(k) and    c
c LWP(k-1) since the w level for a given k is at a higher altitude     c
c than the mass level.  radht, back on mass levels, depends on Frad(k) c
c and Frad(k+1).                                                       c
c                                                                      c
c Additionally, these computations assume that the distance between    c
c mass levels (dsigma) is constant, and that the w levels (spaced by   c
c dsigmw) always fall exactly halfway in between the mass levels.  If  c
c this is not the case, consider dwm to be the distance between a w    c
c level and the mass level below it, and dmw to be the distance        c
c between a mass level and the w level below it.  Then, the            c
c formulation for Frad_LW, for instance, would use a weighted average: c
c                                                                      c
c (dwm/(dwm+dmw)) * lwp(k) + (dmw/(dwm+dmw)) * lwp(k-1)                c
c which, for dwm always == dmw, reduces to                             c
c (1/2) * (lwp(k)) + (1/2) * (lwp(k-1))                                c
c which is identical to the current formulation.                       c
c ((lwp(k)+lwp(k-1))/2)                                                c
c                                                                      c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        !! We only implement this section if we choose not to use the !!
        !! BUGSRAD radiation scheme.                                  !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF ( .not. lbugsrad ) THEN

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! This code transforms these profiles from HOC grid to COAMPS !!
        !! grid.  The COAMPS-grid profiles are then passed to rad_lwsw !!
        !! for implementation.                                         !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k = 1, gr%nnzp
            thm_rad(k) = thlm(gr%nnzp-k+1) 
     .                 + Lv/(Cp*exner(k))*rcm(gr%nnzp-k+1)
            rcm_rad(k) = rcm(gr%nnzp-k+1)
            rhot_rad(k) = rhot(gr%nnzp-k+1)
            dsigm(k) = 1.0 / gr%dzt(gr%nnzp-k+1)
          end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Calling the radiation subroutine, which uses the COAMPS    !!
        !! grid method.  All input and output profiles use the COAMPS !!
        !! grid setup.                                                !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call rad_lwsw(thm_rad, rcm_rad, rhot_rad, dsigm,
     .                  Frad_out, Frad_LW_out, Frad_SW_out,
     .                  radhtk, radht_LW_out, radht_SW_out,
     .                  gr%nnzp-1, center,
     .                  xi_abs, F0, F1, kap, radius, A, gc, Fs0, omega,
     .                  sw_on, lw_on)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! This code transforms the radiation results back into HOC !!
        !! grid setup.  These Frad and radht arrays are actually    !!
        !! applied to the HOC model.                                !!
        !!                                                          !!
        !! The radht results are initially calculated in terms of   !!
        !! standard temperature (T).  However, HOC calculates       !!
        !! temperature in terms of potential temperature (theta).   !!
        !! Therefore, we multiply all radht results by (1.0/exner)  !!
        !! to convert from T to theta.                              !! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k = 1, gr%nnzp-1
            Frad(k)     = Frad_out(gr%nnzp-k+1)
            Frad_LW(k)  = Frad_LW_out(gr%nnzp-k+1)
            Frad_SW(k)  = Frad_SW_out(gr%nnzp-k+1)

            radht(k)    = ( 1.0/exner(k) ) * radhtk(gr%nnzp-k+1)
            radht_LW(k) = ( 1.0/exner(k) ) * radht_LW_out(gr%nnzp-k+1)
            radht_SW(k) = ( 1.0/exner(k) ) * radht_SW_out(gr%nnzp-k+1)
          end do

          Frad(1) = Frad(2)
          Frad_LW(1) = Frad_LW(2)
          Frad_SW(1) = Frad_SW(2)

          radht(1) = radht(2)
          radht_LW(1) = radht_LW(2)
          radht_SW(1) = radht_SW(2)

        END IF ! ~ lbugsrad

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Compute the loss of total water due to diffusional       !!
        !! growth of ice.  This is defined on thermodynamic levels. !!
        !!                                                          !!
        !! NOTE FOR JUN.25 CASE                                     !!
        !! --------------------                                     !!
        !! Ice diffusion is generally not used in the Jun.25 case.  !!
        !! The code is still available here if needed.              !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (time >= time_initial + 3600.0 .and. icedfs_on) THEN

           CALL altocu_icedf( dt, thlm, rcm, p, exner, rhot,
     .                              dqc_dt_icedfs )

        ELSE

           DO k = 1, gr%nnzp, 1
              dqc_dt_icedfs(k) = 0.0
           ENDDO

        END IF


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                    c
c Using linear interpolation scheme to interpolate subsidence        c
c                                                                    c
c FOR NOV.11 CASE                                                    c
c ---------------                                                    c
c As mentioned above, we want to implement a constant subsidence     c
c profile throughout the entire simulation (except for the 1 hour    c
c initial "spinup" period).  Because we do not have variations in    c
c subsidence over time, all tsubs sections (used in Jun.25 case)     c
c have been removed below.  Only one loop remains below, which will  c
c implement the interpolation subroutine, then return the w_ls       c
c value.  Unlike jun25, we use this w_ls value directly instead of   c
c using interpolation to calculate a value between time steps.       c
c                                                                    c
c DIAGRAM OF NOV. 11 SUBSIDENCE PROFILE                              c
c -------------------------------------                              c
c       |      <- No subsidence in this region                       c
c       |                                                            c
c       ------------ Height = (zi + dazi + dac)                      c
c        \                                                           c
c         \    <- Subsidence tapers linearly in this region          c
c          \                                                         c
c           -------- Height = (zi + dazi)                            c
c           |                                                        c
c           |  <- Subsidence equals wmax in this region              c
c           |                                                        c
c           -------- Height = zi                                     c
c           |                                                        c
c           -------- Height = (zi - dbzi)                            c
c          /                                                         c
c         /    <- Subsidence tapers linearly in this region          c
c        /                                                           c
c       ------------ Height = (zi - dbzi - dbc)                      c
c       |                                                            c
c       |      <- No subsidence in this region                       c
c                                                                    c
c                                                                    c
c Comment by Adam Smith on 26 June 2006                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do k=2,gr%nnzp
          if ( (time >= time_initial + 3600.0) .and. subs_on ) then
            call linearinterpolation(7,zsubs,wt1,gr%zt(k),wmt(k))
          else
c           If time is not yet one hour, we have no subsidence
            wmt(k) = 0.0
          end if

          wmt(1) = wmt(2)
        end do

        wmm = zt2zm(wmt)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Enter the final theta-l and rtm tendencies !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO k = 1, gr%nnzp, 1
           IF ( .not. lbugsrad ) THEN
             thl_tndcy(k) = radht(k)
     .                    - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           ELSE
             thl_tndcy(k) = - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           END IF
           rt_tndcy(k) = dqc_dt_icedfs(k)
        ENDDO

        IF ( present( sclr_tndcy ) ) THEN
          sclr_tndcy(:,SCLR_THETA) = thl_tndcy
          sclr_tndcy(:,SCLR_RT)    = rt_tndcy
        END IF

#ifdef STATS
        if ( .not.lbugsrad .and. lstats_samp ) then

          if ( iradht_LW > 0 ) then
            zt%x(:,iradht_LW) = zt%x(:,iradht_LW) + radht_LW
            zt%n(:,iradht_LW) = zt%n(:,iradht_LW) + 1
          end if
          if ( iradht_SW > 0 ) then
            zt%x(:,iradht_SW) = zt%x(:,iradht_SW) + radht_SW
            zt%n(:,iradht_SW) = zt%n(:,iradht_SW) + 1
          end if

          if ( iFrad_SW > 0 ) then
            zm%x(:,iFrad_SW) = zm%x(:,iFrad_SW) + Frad_SW
            zm%n(:,iFrad_SW) = zm%n(:,iFrad_SW) + 1
          end if
          if ( iFrad_LW > 0 ) then
            zm%x(:,iFrad_LW) = zm%x(:,iFrad_LW) + Frad_LW
            zm%n(:,iFrad_LW) = zm%n(:,iFrad_LW) + 1
          end if

        end if
#endif /*STATS*/

        RETURN
        END SUBROUTINE nov11_altocu_tndcy

!-----------------------------------------------------------------------

        ! This subroutine is based on a COAMPS subroutine (nov11_icedfs)
        ! written by Adam Smith and Vince Larson to calculate the
        ! diffusional growth of ice for the Nov. 11 altocumulus case.

        !---------------Brian's comment--------------------------------------!
        ! This code does not use actual microphysics.  Diffusional growth of !
        ! ice is supposed to be the growth of ice due to diffusion of water  !
        ! vapor.  Liquid water is not involved in diffusional growth.        !
        ! However, in mixed phase clouds (both ice and liquid water), most   !
        ! of the water vapor condenses onto the liquid droplets due to the   !
        ! fact that they have so much more available surface area.  This     !
        ! brings the amount of water vapor in the atmosphere to the          !
        ! saturation level with respect to liquid water.  However, since the !
        ! saturation vapor pressure with respect to ice is less than the     !
        ! saturation vapor pressure with respect to liquid water, a          !
        ! saturated atmosphere with respect to liquid water is still         !
        ! supersaturated with respect to ice.  As a result, ice still grows  !
        ! due to diffusion.  When this happens, the environmental vapor      !
        ! pressure drops to the point of saturation with respect to ice.     !
        ! This leaves the atmosphere subsaturated with respect to liquid     !
        ! water.  As a result, some of the liquid water evaporates until     !
        ! the atmosphere becomes saturated with respect to liquid water      !
        ! again.  The process then repeats itself.  As a result, the ice     !
        ! essentially grows at the expense of the liquid water.  This is     !
        ! why the diffusional growth of ice is being deducted from liquid    !
        ! water in this subroutine.
        !--------------------------------------------------------------------!

        SUBROUTINE altocu_icedf( dt, thlm, rcm, p, exner, rhot,
     .                           dqc_dt_icedfs )

        USE grid_class
        USE constants
        USE statistics

        implicit none

        external rsat
        real rsat

        ! Input variables
        REAL, INTENT(IN):: dt
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: thlm
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rcm
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: p
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: exner
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rhot

        ! Output variables
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: dqc_dt_icedfs

        ! Local variables
        REAL, DIMENSION(1:gr%nnzp):: m
        REAL, DIMENSION(1:gr%nnzp):: thm
        REAL, DIMENSION(1:gr%nnzp):: TT
        REAL, DIMENSION(1:gr%nnzp):: rs
        REAL, DIMENSION(1:gr%nnzp):: es
        REAL, DIMENSION(1:gr%nnzp):: ei
        REAL, DIMENSION(1:gr%nnzp):: Si
        REAL, DIMENSION(1:gr%nnzp):: Denom
        REAL, DIMENSION(1:gr%nnzp):: dm
        REAL, DIMENSION(1:gr%nnzp):: diam
        REAL, DIMENSION(1:gr%nnzp):: u_T_cm

        REAL:: a, b, k_u, q, n

        ! For function call
!        REAL:: Diff_denom

        REAL, PARAMETER:: Lf = 3.33e5
        ! Vince Larson avgd legs 2 and 7 (Fleishauer et al)  21 Jan 2005
        REAL, PARAMETER:: Ni = 2000

        INTEGER:: k

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                     !
        ! Coefficients for mass-diameter relationship, Mitchell (1996)        !
        ! mass = a (diameter/(1 meter))^b,  [a] = kg, [b] = []                !
        !                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        a = 2.05e-3
        b = 1.8

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                     !
        !  Coefficients for mass-diameter relationship, Kajikawa (1989)       !
        !  mass = a (diam/(1m))^b,  [a] = kg, [b] = []                        !
        !                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       a = 2.50e-4
!       b = 1.4

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                     !
        !  Coefficients for fallspeed-diameter relationship, Mitchell (1996)  !
        !  u_T = k_u rho^{-q} (diameter/(1 meter))^n,                         !
        !       [k_u] = m/s, [q] = [], [n] = [], [rho] = kg m^{-3}            !
        !                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        k_u = 55.
        q   = 0.17
        n   = 0.70

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                     !
        !  Coefficients for fallspeed-diameter relationship, Kajikawa (1989)  !
        !  u_T = k_u rho^{-q} (diam/(1m))^n,  [k_u] = m/s, [q] = [], [n] = [] !
        !       [rho] = kg m^{-3}                                             !
        !                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       k_u = 0.438
!       q = 0.0
!       n = 0.0742



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                     !
        !  Initialize ice particle mass                                       !
        !                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DO k = 1, gr%nnzp, 1
           m(k) = 1.0e-11
        ENDDO

        DO k = gr%nnzp, 2, -1

          !  Compute Temperature [K]
          thm(k) = thlm(k) + ( Lv / (Cp*exner(k)) ) * rcm(k)
          TT(k) = thm(k) * exner(k)

          ! Check whether we're in cloud and below freezing.
          ! Note:  A value of 1.0E-5 kg/kg is used as a threshold value
          ! for rcm because the HOC model shows a small amount of liquid
          ! water all the way to the model top, which messes with the
          ! ice diffusion calculations.
          IF ( rcm(k) >= 1.0E-5 .AND. TT(k) < 273.15 ) THEN

          ! Find saturation mixing ratio over vapor [kg kg^{-1}]
             rs(k) = rsat( p(k), TT(k) )

          ! Saturation vapor pressure over liquid in Pa
             es(k) = ( rs(k)*p(k) ) / ( ep + rs(k) )

          ! Saturation vapor pressure over ice in Pa, Eq. 2.15 Rogers and Yau
             ei(k) = es(k) / EXP( ( Lf/(Rv*273.16) )
     .                           * ( 273.16/TT(k) - 1.0 ) )

          ! Saturation ratio in a liquid-saturated cloud, p. 158 Rogers and Yau
          !---------------Brian's comment--------------------------------------!
          ! The actual formula is:  Si = e/ei = (e/es)*(es/ei) = S*(es/ei)     !
          ! It is assumed that any supersaturation forms liquid water and that !
          ! the atmosphere is then saturated with respect to liquid water.     !
          ! Therefore, S = 1.0, allowing Si = es/ei.                           !
          !--------------------------------------------------------------------!
             Si(k) = es(k)/ei(k)

          ! Denominator of diffusional growth equation, 9.4 of Rogers and Yau
             Denom(k) = Diff_denom( TT(k), p(k), ei(k) )

          ! Change in mass of a single ice crystal, m,
          ! as it falls a distance gr%dzt in meters

          !---------------Brian's comment--------------------------------------!
          ! dm/dt = 4*pi*C*(Si-1)/Denom; Rogers and Yau, Eq. 9.4.              !
          ! For plate-type ice crystals, C = 2r/pi (Rogers and Yau, p. 159).   !
          ! Since 2r = D, C = D/pi, and the equation becomes:                  !
          ! dm/dt = 4*D*(Si-1)/Denom.                                          !
          ! The mass-diameter relationship for an ice crystal is:              !
          ! D = (m/a)^(1/b); Rogers and Yau, Eq. 9.7.  This means:             !
          ! dm/dt = [(m/a)^(1/b)]*4*(Si-1)/Denom;                              !
          ! Dividing by rho yields the change in mixing ratio over time        !
          ! for an individual crystal.  Multiplying that by the ice crystal    !
          ! concentration yields the overall change in mixing ratio over time. !
          !--------------------------------------------------------------------!
             dqc_dt_icedfs(k) = - (Ni/rhot(k))
     .          * ( 4 * (Si(k) - 1) / Denom(k) )
     .          * (m(k)/a)**(1/b)

             ! Ensure that liquid is not over-depleted
             IF ( rcm(k) + dqc_dt_icedfs(k)*dt < 0.0 ) THEN
                dqc_dt_icedfs(k) = -rcm(k)/dt
             END IF

             !---------------Brian's comment-----------------------------------!
             ! dm = (dm/dt)*(dt/dz)*dz                                         !
             ! dm = (dm/dt)*(1/u_T)*dz                                         !
             !-----------------------------------------------------------------!
             dm(k) = ( 4 * (Si(k) - 1) / Denom(k) )
     .          * (k_u**(-1.0)) * (rhot(k)**q)
     .          * ( (m(k)/a)**((1.0-n)/b) )
     .          * (1.0/gr%dzm(k-1))
             m(k-1) = m(k) + dm(k)

             ! Diameter of ice crystal in meters.
             diam(k) = (m(k)/a)**(1/b)

             ! Fallspeed of ice crystal in cm/s.
             u_T_cm(k) = 100. * k_u *
     .                   ((m(k)/a)**(n/b)) * (rhot(k)**(-q))

          ELSE   ! There's no liquid and/or ice present; assume no ice growth

             m(k-1) = m(k)
             dqc_dt_icedfs(k) = 0.0
             diam(k)          = 0.0  ! Set zero to remind that we don't grow ice
             u_T_cm(k)        = 0.0  ! Set zero to remind that we don't grow ice

          END IF

        ENDDO

! Michael Falk added boundary condx, 31 July 2006

        m(1) = m(2)
        dqc_dt_icedfs(1) = dqc_dt_icedfs(2)
        diam(1) = diam(2)
        u_T_cm(1) = u_T_cm(2)

! eMFc

!
#ifdef STATS
        if ( lstats_samp ) then
!       diam(:) ! Icedfs diameter; Michael Falk, 1 Nov 2006
!       m(:)    ! Icedfs mass; Michael Falk, 1 Nov 2006
!       dqc_dt_icedfs(:) ! Icedfs change in liquid; Michael Falk, 1 Nov 2006
!       u_T_cm(:)        ! Icedfs fallspeed (cm/s); Michael Falk, 1 Nov 2006

          if ( idqc_dt_icedfs > 0 ) then
            zt%x(:,idqc_dt_icedfs) = zt%x(:,idqc_dt_icedfs)
     .                               + dqc_dt_icedfs
            zt%n(:,idqc_dt_icedfs) = zt%n(:,idqc_dt_icedfs) + 1
          end if
          if ( idiam > 0 ) then
            zt%x(:,idiam) = zt%x(:,idiam) + diam
            zt%n(:,idiam) = zt%n(:,idiam) + 1
          end if
          if ( im > 0 ) then
            zt%x(:,im) = zt%x(:,im) + m
            zt%n(:,im) = zt%n(:,im) + 1
          end if
          if ( iu_T_cm > 0 ) then
            zt%x(:,iu_T_cm) = zt%x(:,iu_T_cm) + u_T_cm
            zt%n(:,iu_T_cm) = zt%n(:,iu_T_cm) + 1
          end if
        end if
#endif /*STATS*/

        RETURN
        END SUBROUTINE altocu_icedf

!-----------------------------------------------------------------------
        subroutine jun25_altocu_tndcy
     .             ( time, time_initial, rlat, dt, thlm, 
     .               rcm, p, exner, rhot, rtm, wmt,
     .               wmm, thl_tndcy, rt_tndcy, 
     .               Frad, radht, 
     .               sclr_tndcy ) ! Optional
!       Description:
!       References:
!-----------------------------------------------------------------------
        use grid_class
        use constants
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

        ! Input variables
        real, intent(in) :: 
     .  time,         ! Time of simulation since start        [s]
     .  time_initial, ! Initial time of simulation            [s]
     .  rlat,         ! Reference latitude should be 37.6     [Degrees North]
     .  dt            ! Length of model timestep              [s]

        real, dimension(gr%nnzp), intent(in) :: 
     .  thlm,  ! Liquid potential temperature           [K]
     .  rcm,   ! Liquid water mixing ratio              [kg/kg]
     .  p,     ! Basic state pressure                   [Pa]
     .  exner, ! Exner function                         [-]
     .  rhot   ! Density of reference state on t grid   [kg/m^3]

        ! Input/Output Variables
        real, dimension(gr%nnzp), intent(inout) :: 
     .  rtm    ! Total water mixing ratio               [kg/kg]

        ! Output variables
        real, dimension(gr%nnzp), intent(inout) :: 
     .  wmt,      ! Vertical ascent/descent on therm. grid      [m/s]
     .  wmm,      ! Vertical ascent/descent on moment. grid     [m/s]
     .  thl_tndcy ! Change in liq. water potential temperature 
                  ! due to radiative heating and ice diffusion  [K/s]

        real, dimension(gr%nnzp), intent(inout) :: 
     .  rt_tndcy, ! Change in total water due to ice diffusion  [kg/kg/s]
     .  Frad,     ! Total radiative flux (LW + SW)              [W/m^2]
     .  radht     ! Total radiative heating (LW +SW)            [K/s]

        ! Output variables (optional)
        real, optional, dimension(gr%nnzp,sclr_dim),intent(out) :: 
     .  sclr_tndcy ! Large-scale tendency for passive scalars

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c LOCAL VARIABLES                                                             c
c dqc_dt_icedfs : change in liquid water over time, due to ice diffusional    c
c                 growth.                                       Unit: 1/s     c
c                                                                             c
c LWP           : liquid water path from domain top.            Unit: kg/m^2  c
c                                                                             c
c TEMPORARY ARRAYS USED FOR COAMPS RADIATIVE SCHEME                           c
c (SEE COMMENTS BELOW FOR SCHEME DESCRIPTION)                                 c
c thm_rad       : temporary "flipped" array of potential temperature          c
c                                                               Unit: K       c
c rcm_rad       : temporary "flipped" array of liquid water mixing ratio      c
c                                                               Unit: kg/kg   c
c rhot_rad      : temporary "flipped" array of air density      Unit: kg/m^3  c
c dsigm         : temporary "flipped" array of grid spacing     Unit: m       c
c                                                                             c
c frad_out      : temporary "flipped" array of total radiative flux           c
c                 from rad_lwsw                                 Unit: W/m^2   c
c frad_lw_out   : temporary "flipped" array of LW radiative flux              c
c                 from rad_lwsw                                 Unit: W/m^2   c
c frad_sw_out   : temporary "flipped" array of SW radiative flux              c
c                 from rad_lwsw                                 Unit: W/m^2   c
c                                                                             c
c radhtk        : temporary "flipped" array of total radiative heating        c
c                 from rad_lwsw                                 Unit: K/s     c
c radht_lw_out  : temporary "flipped" array of LW radiative heating           c
c                 from rad_lwsw                                 Unit: K/s     c
c radht_sw_out  : temporary "flipped" array of SW radiative heating           c
c                 from rad_lwsw                                 Unit: K/s     c
c                                                                             c
c INTERPOLATION ARRAYS AND CONSTANTS                                          c
c zsubs         : heights at which wmt data is supplied                       c
c                 (used for subsidence interpolation)           Unit: m       c
c tsubs         : times after initialization at which wmt data is supplied    c
c                 (NOT USED IN NOV.11 CASE)                     Unit: s       c
c wtX(Y)        : vertical velocity specified at height Y and time X          c
c                 (ONLY wt1 IS USED IN NOV.11 CASE)             Unit: m/s     c
c w1-w2         : vertical velocity before (w1) and after (w2)                c
c                 the current time at the specified level                     c
c                 (NOT USED IN NOV.11 CASE)                     Unit: m/s     c
c                                                                             c
c ADDITIONAL PARAMETERS FOR NOV.11 SUBSIDENCE (NOT FOR JUN.25 CASE)           c
c subs_on       : logical variable tells us whether to turn subsidence on     c
c                                                               Unit: NONE    c
c wmax          : defines value of maximum subsidence in profile              c
c                                                               Unit: cm/s    c
c zi            : defines approx. height of inversion within cloud            c
c                 (subsidence is equal to wmax at this height)  Unit: m       c
c dazi          : defines height above inversion                              c
c                 (above this height, subsidence linearly tapers off to zero) c
c                                                               Unit: m       c
c dbzi          : defines height above inversion                              c
c                 (below this height, subsidence linearly tapers off to zero) c
c                                                               Unit: m       c
c dac           : defines height above cloud                                  c
c                 (at / above this height, we have NO subsidence)             c
c                                                               Unit: m       c
c dbc           : defines height below cloud                                  c
c                 (at / below this height, we have NO subsidence)             c
c                                                               Unit: m       c
c                                                                             c
c RADIATION PARAMETERS                                                        c
c sw_on         : logical variable passed to radiation scheme                 c
c                 - is SW radiation on?                         Unit: NONE    c
c lw_on         : logical variable passed to radiation scheme                 c
c                 - is LW radiation on?                         Unit: NONE    c
c center        : use centered differencing (as opposed to a one-sided        c
c                 forward difference) in radiation code         Unit: NONE    c
c                                                                             c
c xi_abs        : cosine of the solar zenith angle              Unit: NONE    c
c F0            : coefficient for cloud top heating (see Stevens)             c
c                                                               Unit: W/m^2   c
c F1            : coefficient for cloud base heating (see Stevens)            c
c                                                               Unit: W/m^2   c
c kap           : "a constant" according to Duynkerke eqn. 5, where his       c
c                 value is 130 m^2/kg.                          Unit: m^2/kg  c
c radius        : effective droplet radius                      Unit: m       c
c AA            : albedo -- sea surface, according to Lenderink.              c
c                                                               Unit: NONE    c
c gc            : asymmetry parameter, "g" in Duynkerke.        Unit: NONE    c
c Fs0           : The incident incoming SW insolation at cloud top in the     c
c                 direction of the incoming beam (not the vertical).          c
c                                                               Unit: W/m^2   c
c omega         : single-scattering albedo                      Unit: NONE    c
c                                                                             c
c SOLAR ZENITH ANGLE PARAMETERS (NOT USED IN NOV.11 CASE)                     c
c c0            : coefficient for calculation of declination angle from       c
c                 Liou Table 2.2 and Eqn. 2.2.10                Unit: NONE    c
c c1            : same as above                                 Unit: NONE    c
c c2            : same as above                                 Unit: NONE    c
c c3            : same as above                                 Unit: NONE    c
c d1            : Same as above                                 Unit: NONE    c
c d2            : Same as above                                 Unit: NONE    c
c d3            : Same as above                                 Unit: NONE    c
c                                                                             c
c sda_t         : Linear function of day of the year.                         c
c                 sda_t=0 January 1 and sda_t -> 2*pi December 31.            c
c sda_delta     : Solar declination angle function from Liou 2.2.10           c
c sda_h         : Hour angle (Angle through which the earth must to turn      c
c                 to put sun directly overhead on a point's meridian)         c
c                 (Angle between current time and solar noon)                 c
c t_since_noon  : Number of seconds since noon (after noon > 0)               c
c julday        : Julian day of the year (January 1=1; December 31=365)       c
c                                                                             c
c start_time_until_noon: number of seconds between start time and solar noon  c
c                                                               Unit: s       c
c                                                                             c
c Fs0 INTERPOLATION PARAMETERS                                                c
c nparam        : Number of Fs0 values in the list.                           c
c xilist        : Values of cosine of solar zenith angle corresponding to     c
c                 the values in Fslist                                        c
c Fslist        : Values of Fs0 corresponding to the values in xilist.        c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Local ice diffusional arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, DIMENSION(gr%nnzp):: dqc_dt_icedfs

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Local radiation arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real, dimension(gr%nnzp) :: 
     .  Frad_LW, ! Long wave radiative flux     [W/m^2]
     .  Frad_SW, ! Short wave radiative flux    [W/m^2]
     .  radht_LW,! Long wave radiative heating  [K/s]
     .  radht_SW ! Short wave radiative heating [K/s]

        REAL, DIMENSION(gr%nnzp):: LWP
        REAL, DIMENSION(gr%nnzp):: thm_rad
        REAL, DIMENSION(gr%nnzp):: rcm_rad
        REAL, DIMENSION(gr%nnzp):: rhot_rad
        REAL, DIMENSION(gr%nnzp):: dsigm

        REAL, DIMENSION(gr%nnzp):: frad_out
        REAL, DIMENSION(gr%nnzp):: frad_lw_out
        REAL, DIMENSION(gr%nnzp):: frad_sw_out

        REAL, DIMENSION(gr%nnzp):: radhtk
        REAL, DIMENSION(gr%nnzp):: radht_lw_out
        REAL, DIMENSION(gr%nnzp):: radht_sw_out

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Working arrays for subsidence interpolation !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL zsubs(5)
        REAL tsubs(6)
        REAL wt1(5), wt2(5), wt3(5), wt4(5), wt5(5), wt6(5)
        REAL w1(gr%nnzp), w2(gr%nnzp)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! LW Radiative constants !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, PARAMETER :: F0   = 107.0  ! W/m^2
        REAL, PARAMETER :: F1   = 61.0   ! W/m^2
        REAL, PARAMETER :: kap  = 100.0  ! m^2/kg

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Working arrays for SW radiation interpolation !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER, PARAMETER :: nparam = 12
        REAL, DIMENSION(nparam) :: xilist, Fslist

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! SW Radiative constants !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, PARAMETER :: radius = 1.0e-5
        REAL, PARAMETER :: AA     = 0.1
        REAL, PARAMETER :: gc     = 0.85
        REAL, PARAMETER :: omega  = 0.992

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Additional SW radiative variables !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL xi_abs, Fs0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variables used for calculating solar zenith angle                 !!
        !! (Currently not used in Nov.11 case)                               !!
        !!                                                                   !!
        !! NOTE BY ADAM SMITH, 13 April 2006                                 !!
        !! In COAMPS, these constants are called c0, c1, c2, c3, d1, d2,     !!
        !! and d3.  Since HOC uses these variables names elsewhere, I have   !!
        !! renamed them here.                                                !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL c_0, c_1, c_2, c_3, d_1, d_2, d_3
        REAL sda_t, sda_delta, sda_h, t_since_noon, julday
        REAL start_time_until_noon

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Toggle for implementing differencing method in interpolations !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LOGICAL :: center

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Toggles for activating/deactivating forcings                  !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LOGICAL :: icedfs_on, subs_on, lw_on, sw_on

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variable used for working within vertical arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER:: k

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! END OF VARIABLE DECLARATION !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Toggles for activating/deactivating forcings                            cc
cc To turn off a specific forcing, set the corresponding toggle to .FALSE. cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         icedfs_on = .FALSE.
         subs_on   = .TRUE.
         lw_on     = .TRUE.
         sw_on     = .TRUE.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc Toggle for centered/forward differencing (in interpolations) cc
cc To use centered differencing, set the toggle to .TRUE.       cc
cc To use forward differencing, set the toggle to .FALSE.       cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         center    = .TRUE.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Calculation of Solar Declination Angle                               c
c                                                                      c
c This section of code calculates the solar declination angle based    c
c on the simulaiton's time, date, and latitude.  These parameters are  c
c applied in the scheme, and a value of xi_abs is calculated.          c
c                                                                      c
c If xi_abs <= 0.0 at any time step, we shut off shortwave radiation   c
c for that time step (by setting sw_on = .FALSE.).  If xi_abs > 0.0,   c
c we reactivate solar radiation by setting sw_on = .TRUE.              c
c                                                                      c
c This version is based on solar_declin.m, written by Dr. Vince Larson c
c and based on Liou 2.2.10.                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        julday = 176
        start_time_until_noon = 30540
        t_since_noon   = time - start_time_until_noon
        sda_t = 2.0*pi*(julday-1)/365

        c_0 = 0.006918
        c_1 = -0.399912
        c_2 = -0.006758
        c_3 = -0.002697
        d_1 = 0.070257
        d_2 = 0.000907
        d_3 = 0.000148

        sda_delta = c_0 + c_1*cos(sda_t)   + d_1*sin(sda_t)   +
     .                    c_2*cos(2*sda_t) + d_2*sin(2*sda_t) +
     .                    c_3*cos(3*sda_t) + d_3*sin(3*sda_t)

        sda_h = 2.0*pi*t_since_noon/86400.0

        xi_abs = sin(rlat*pi/180) * sin(sda_delta) +
     .           cos(rlat*pi/180) * cos(sda_delta) * cos(sda_h)

        xi_abs = max(xi_abs,0.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Modification by Adam Smith 26 June 2006                                   !!
!! It is difficult to remember to set xi_abs = 0 when we want to shut off    !!
!! solar radiation.  If sw_on = .FALSE. above, we will automatically set     !!
!! xi_abs to 0 to avoid confusion or errors.                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( .not. sw_on ) then
        xi_abs = 0.
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End of ajsmith4's Modification !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (xi_abs == 0.) then
        sw_on = .FALSE.
      else
        sw_on = .TRUE.
      end if


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Fs0 Interpolation Parameters-- these also from Kurt Kotenberg's      c
c BUGSrad output.  Fs0 changes somewhat over the range of solar zenith c
c angles, and we obtained these values by solving                      c
c Fs0 = F_vertical / xi_abs .                                          c
c                                                                      c
c The linearinterpolation function returns Fs0.                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      xilist(1) = 0.0
      xilist(2) = 0.01
      xilist(3) = 0.1
      xilist(4) = 0.2
      xilist(5) = 0.3
      xilist(6) = 0.4
      xilist(7) = 0.5
      xilist(8) = 0.6
      xilist(9) = 0.7
      xilist(10) = 0.8
      xilist(11) = 0.9
      xilist(12) = 1.0

      Fslist(1)  = 0.0
      Fslist(2)  = 715.86
      Fslist(3)  = 1073.577
      Fslist(4)  = 1165.0905
      Fslist(5)  = 1204.7033
      Fslist(6)  = 1227.6898
      Fslist(7)  = 1243.1772
      Fslist(8)  = 1254.5893
      Fslist(9)  = 1263.5491
      Fslist(10) = 1270.8668
      Fslist(11) = 1277.0474
      Fslist(12) = 1282.3994

      call linearinterpolation( nparam, xilist, Fslist, xi_abs, Fs0 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subsidence Parameters !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        zsubs(1) = 0
        zsubs(2) = 360
        zsubs(3) = 1090
        zsubs(4) = 1890
        zsubs(5) = 2500

        tsubs(1) = 0
        tsubs(2) = 10800
        tsubs(3) = 28800
        tsubs(4) = 36000
        tsubs(5) = 36000
        tsubs(6) = 36000

        wt1(1) = 0.
        wt1(2) = .004
        wt1(3) = .004
        wt1(4) = .004
        wt1(5) = 0.

        wt2(1) = 0.
        wt2(2) = .004
        wt2(3) = .004
        wt2(4) = .004
        wt2(5) = 0.

        wt3(1) = 0.
        wt3(2) = -.003
        wt3(3) = -.003
        wt3(4) = -.003
        wt3(5) = 0.

        wt4(1) = 0.
        wt4(2) = -.003
        wt4(3) = -.003
        wt4(4) = -.003
        wt4(5) = 0.

        wt5(1) = 0.
        wt5(2) = -.003
        wt5(3) = -.003
        wt5(4) = -.003
        wt5(5) = 0.

        wt6(1) = 0.
        wt6(2) = -.003
        wt6(3) = -.003
        wt6(4) = -.003
        wt6(5) = 0.


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c SPECIAL METHOD USED TO CALCULATE RADIATION                           c
c Grid descriptions by Adam Smith, 27 June 2006                        c
c                                                                      c
c In order to verify our HOC simulations are working properly, we      c
c have first developed a series of 3D simulations using the COAMPS-LES c
c model.  This large-eddy simulation (LES) simulation uses specific    c
c methods to calculate radiation, subsidence, and other microphysical  c
c processes.  To make the two models simluate clouds as closely as     c
c possible, we use the same radiation scheme in both models.           c
c                                                                      c
c In COAMPS-LES, we use a separate subroutine, rad_lwsw, to implement  c
c all radiation code.  This allows the subroutine to be duplicated     c
c exactly in many different models.  However, the subroutine uses the  c
c COAMPS vertical grid.  Therefore, for HOC to implement this code     c
c correctly, we must modify some of our variable profiles before       c
c calling the radiation subroutine.                                    c
c                                                                      c
c The following diagram describes the differences in model grids:      c
c                                                                      c
c       COAMPS-LES                                    HOC              c
c                                                                      c
c  k= 1   (m) ----------    <MODEL TOP>    ---------- k=kk+1      (m)  c
c                                                                      c
c  k= 1   (t) ----------                   ---------- k=kk+1      (t)  c
c                                                                      c
c  k= 2   (m) ----------                   ---------- k=kk        (m)  c
c                                                                      c
c  k= 2   (t) ----------                   ---------- k=kk        (t)  c
c                                                                      c
c            .                  .                     .                c
c            .                  .                     .                c
c            .                  .                     .                c
c                                                                      c
c  k=kk-1 (m) ----------  m = mass levels  ---------- k=3         (m)  c
c                                                                      c
c  k=kk-1 (t) ----------  t = thermo       ---------- k=3         (t)  c
c                                  levels                              c
c  k=kk   (m) ----------                   ---------- k=2         (m)  c
c                                                                      c
c  k=kk   (t) ----------  kk = number of   ---------- k=2         (t)  c
c                              vertical                                c
c  k=kk+1 (m) ----------       heights     ---------- k=1         (m)  c
c                                                                      c
c //////////////////////// MODEL SURFACE ///////////////////////////// c
c                                          ---------- k=1         (t)  c
c                                                                      c
c                                                                      c
c The major difference in the grids is that HOC uses an additional     c
c thermodynamic level below the model "surface".  This means that all  c
c HOC thermodynamic heights are shifted down one vertical level, and   c
c HOC also has one fewer momentum level than COAMPS.  Therefore, we    c
c use one additional vertical level in HOC, to make sure that the      c
c vertical domain matches in both models.                              c
c                                                                      c
c Also, the COAMPS grid indices are numbered from the top of the model c
c downward, while the HOC grid indices are numbered from the bottom    c
c up.  Therefore, since we are using a COAMPS radiation scheme, we     c
c flip moisture and temperature profiles that are passed into the      c
c rad_lwsw subroutine.  The rad scheme will produce results in using   c
c the COAMPS grid scheme, so all radiation output will be flipped      c
c back to the HOC grid before being applied to the model.              c
c                                                                      c
c Finally, since the COAMPS scheme does not have a gridpoint below     c
c model surface, we add that point to all radiative output files once  c
c they are converted back to HOC setup.  This allows all averages and  c
c calculations to be done correctly.                                   c
c                                                                      c
c                                                                      c
c Computation of radiative fluxes on staggered grid                    c
c Comments by Michael Falk, 16 February 2005.                          c
c                                                                      c
c Frad (and its components Frad_LW and Frad_SW) should be computed on  c
c w points, not on mass points, which is apparent from its formulation c
c and from its location in stats_sw instead of stats_sm.  The grid     c
c looks like this:                                                     c
c                                                                      c
c                                                                      c
c -----Frad----------------------------------    k = 1  (w level)      c
c     /    \            |-dwm                                          c
c -LWP------radht----------------------------    k = 1  (mass level)   c
c     \    /            |-dmw                                          c
c -----Frad----------------------------------    k = 2  (w level)      c
c     /    \                                                           c
c -LWP------radht----------------------------    k = 2  (mass level)   c
c     \    /                                                           c
c -----Frad----------------------------------    k = 3  (w level)      c
c     /    \                                                           c
c -LWP------radht----------------------------    k = 3  (mass level)   c
c                                                                      c
c If you consider Frad to take place on mass levels, then computing    c
c LWP is a forward difference and is only first-order accurate, while  c
c if Frad computed in between LWP levels, it is a centered difference  c
c which is second-order accurate.                                      c
c                                                                      c
c The coding implementation requires that Frad depend on LWP(k) and    c
c LWP(k-1) since the w level for a given k is at a higher altitude     c
c than the mass level.  radht, back on mass levels, depends on Frad(k) c
c and Frad(k+1).                                                       c
c                                                                      c
c Additionally, these computations assume that the distance between    c
c mass levels (dsigma) is constant, and that the w levels (spaced by   c
c dsigmw) always fall exactly halfway in between the mass levels.  If  c
c this is not the case, consider dwm to be the distance between a w    c
c level and the mass level below it, and dmw to be the distance        c
c between a mass level and the w level below it.  Then, the            c
c formulation for Frad_LW, for instance, would use a weighted average: c
c                                                                      c
c (dwm/(dwm+dmw)) * lwp(k) + (dmw/(dwm+dmw)) * lwp(k-1)                c
c which, for dwm always == dmw, reduces to                             c
c (1/2) * (lwp(k)) + (1/2) * (lwp(k-1))                                c
c which is identical to the current formulation.                       c
c ((lwp(k)+lwp(k-1))/2)                                                c
c                                                                      c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        !! We only implement this section if we choose not to use the !!
        !! BUGSRAD radiation scheme.                                  !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IF ( .not. lbugsrad ) THEN

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! This code transforms these profiles from HOC grid to COAMPS !!
        !! grid.  The COAMPS-grid profiles are then passed to rad_lwsw !!
        !! for implementation.                                         !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k = 1, gr%nnzp
            thm_rad(k) = thlm(gr%nnzp-k+1) 
     .                 + Lv/(Cp*exner(k))*rcm(gr%nnzp-k+1)
            rcm_rad(k) = rcm(gr%nnzp-k+1)
            rhot_rad(k) = rhot(gr%nnzp-k+1)
            dsigm(k) = 1.0 / gr%dzt(gr%nnzp-k+1)
          end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Calling the radiation subroutine, which uses the COAMPS    !!
        !! grid method.  All input and output profiles use the COAMPS !!
        !! grid setup.                                                !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call rad_lwsw(thm_rad, rcm_rad, rhot_rad, dsigm,
     .                  Frad_out, Frad_LW_out, Frad_SW_out,
     .                  radhtk, radht_LW_out, radht_SW_out,
     .                  gr%nnzp-1, center,
     .                  xi_abs, F0, F1, kap, radius, AA, gc, Fs0, omega,
     .                  sw_on, lw_on)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! This code transforms the radiation results back into HOC !!
        !! grid setup.  These Frad and radht arrays are actually    !!
        !! applied to the HOC model.                                !!
        !!                                                          !!
        !! The radht results are initially calculated in terms of   !!
        !! standard temperature (T).  However, HOC calculates       !!
        !! temperature in terms of potential temperature (theta).   !!
        !! Therefore, we multiply all radht results by (1.0/exner)  !!
        !! to convert from T to theta.                              !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do k = 1, gr%nnzp-1
            Frad(k)     = Frad_out(gr%nnzp-k+1)
            Frad_LW(k)  = Frad_LW_out(gr%nnzp-k+1)
            Frad_SW(k)  = Frad_SW_out(gr%nnzp-k+1)

            radht(k)    = ( 1.0/exner(k) ) * radhtk(gr%nnzp-k+1)
            radht_LW(k) = ( 1.0/exner(k) ) * radht_LW_out(gr%nnzp-k+1)
            radht_SW(k) = ( 1.0/exner(k) ) * radht_SW_out(gr%nnzp-k+1)
          end do

          Frad(1) = Frad(2)
          Frad_LW(1) = Frad_LW(2)
          Frad_SW(1) = Frad_SW(2)

          radht(1) = radht(2)
          radht_LW(1) = radht_LW(2)
          radht_SW(1) = radht_SW(2)

        END IF ! ~ lbugsrad


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Compute the loss of total water due to diffusional       !!
        !! growth of ice.  This is defined on thermodynamic levels. !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (time >= time_initial + 3600.0 .and. icedfs_on) THEN

           CALL altocu_icedf( dt, thlm, rcm, p, exner, rhot,
     .                        dqc_dt_icedfs )

        ELSE

           DO k = 1, gr%nnzp, 1
              dqc_dt_icedfs(k) = 0.0
           ENDDO

        END IF


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Using linear interpolation to calculate subsidence !!
        !! Original code by Michael Falk                      !!
        !! Added for Jun.25 case by Adam Smith, 13 April 2006 !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (time.lt.tsubs(1)) then
        do k=1,gr%nnzp
            call linearinterpolation(5,zsubs,wt1,gr%zt(k),wmt(k))
        end do
 
        else if (time < tsubs(2)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt1,gr%zt(k),w1(k))
          call linearinterpolation(5,zsubs,wt2,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(1))/(tsubs(2)-tsubs(1))*(w2(k)-w1(k))+w1(k)
        end do
 
        else if (time < tsubs(3)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt2,gr%zt(k),w1(k))
          call linearinterpolation(5,zsubs,wt3,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(2))/(tsubs(3)-tsubs(2))*(w2(k)-w1(k))+w1(k)
        end do
 
        else if (time < tsubs(4)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt3,gr%zt(k),w1(k))
          call linearinterpolation(5,zsubs,wt4,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(3))/(tsubs(4)-tsubs(3))*(w2(k)-w1(k))+w1(k)
        end do
 
        else if (time < tsubs(5)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt4,gr%zt(k),w1(k))
          call linearinterpolation(5,zsubs,wt5,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(4))/(tsubs(5)-tsubs(4))*(w2(k)-w1(k))+w1(k)
        end do
 
        else if (time < tsubs(6)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt5,gr%zt(k),w1(k))
          call linearinterpolation(5,zsubs,wt6,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(5))/(tsubs(6)-tsubs(5))*(w2(k)-w1(k))+w1(k)
        end do
 
        else if (time >= tsubs(6)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt6,gr%zt(k),wmt(k))
        end do
        end if

        wmt(1) = wmt(2)

        wmm = zt2zm(wmt)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Enter the final theta-l and rtm tendencies !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO k = 1, gr%nnzp, 1
           IF ( .not. lbugsrad ) THEN
             thl_tndcy(k) = radht(k)
     .                    - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           ELSE
             thl_tndcy(k) = - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           END IF
           rt_tndcy(k) = dqc_dt_icedfs(k)
        ENDDO

        if ( present( sclr_tndcy ) ) then
          sclr_tndcy(:,SCLR_THETA) = thl_tndcy
          sclr_tndcy(:,SCLR_RT)    = rt_tndcy
        end if 

#ifdef STATS
        if ( .not.lbugsrad .and. lstats_samp ) then

          if ( iradht_LW > 0 ) then
            zt%x(:,iradht_LW) = zt%x(:,iradht_LW) + radht_LW
            zt%n(:,iradht_LW) = zt%n(:,iradht_LW) + 1
          end if
          if ( iradht_SW > 0 ) then
            zt%x(:,iradht_SW) = zt%x(:,iradht_SW) + radht_SW
            zt%n(:,iradht_SW) = zt%n(:,iradht_SW) + 1
          end if

          if ( iFrad_SW > 0 ) then
            zm%x(:,iFrad_SW) = zm%x(:,iFrad_SW) + Frad_SW
            zm%n(:,iFrad_SW) = zm%n(:,iFrad_SW) + 1
          end if
          if ( iFrad_LW > 0 ) then
            zm%x(:,iFrad_LW) = zm%x(:,iFrad_LW) + Frad_LW
            zm%n(:,iFrad_LW) = zm%n(:,iFrad_LW) + 1
          end if

        end if
#endif /*STATS*/

        RETURN
        END SUBROUTINE jun25_altocu_tndcy

!----------------------------------------------------------------------
        subroutine mpace_tndcy
     .  ( time, time_initial, dt, rlat, thlm,
     .    exner, rhot, rtm, p, thvm, rcm,
     .    wmt, wmm, theta_tndcy, rt_tndcy,
     .    Ncnm, Ncm, Frad, radht,
     .    sclr_tndcy )

!        Description:
!          Subroutine to large-scale subsidence for MPACE case (Michael
!          Falk, 21 July 2006).  Added ice and radiation based on Adam
!          Smith Nov 11 case, 27 July 2006.  Comments and documentation
!          added 31 July 2006.
!
!        References:
!          Liou, Wallace and Hobbs, Shettle and Weinman
!-----------------------------------------------------------------------

        use constants
        use grid_class
#ifdef STATS
        use statistics
#endif /*STATS*/
        implicit none

        ! Local constants, subsidence
        real, parameter ::
     .  grav0 = 9.8,     ! m/s
     .  D     = 5.8e-6,  ! 1/s
     .  psfc  = 101000., ! Pa
     .  p00   = 100000., ! Pa; for theta calcs
     .  pinv  = 85000.   ! Pa; ditto

        ! Local constants, LW radiation (from DYCOMS II-RF01)
        real, parameter ::
     .  F0  = 70.0,
     .  F1  = 22.0,
     .  kap = 85.0

        ! Local constants, SW radiation (Shettle and Weinman)
        real, parameter ::
     .  Fs0    = 1212.75,
     .  radius = 1.0e-5,
     .  A      = 0.1,
     .  gc     = 0.86,
     .  omega  = 0.9965
!    .  rlat = 71.75

        ! Local constants, SW radiation (Liou solar angle scheme)
        real, parameter ::
     .  c_0 = 0.006918,
     .  c_1 = -0.399912,
     .  c_2 = -0.006758,
     .  c_3 = -0.002697,
     .  d_1 = 0.070257,
     .  d_2 = 0.000907,
     .  d_3 = 0.000148

        ! Input Variables
        real, intent(in) :: 
     .  time,         ! Current time of simulation      [s]
     .  time_initial, ! Initial time of simulation      [s]
     .  dt,           ! Current length of timestep      [s]
     .  rlat          ! Latitude                        [Degrees North]

        real, dimension(gr%nnzp), intent(in) ::
     .  thlm,  ! Liquid water potential temperature     [K]
     .  exner, ! Exner function                         [-]
     .  rhot,  ! Density of air                         [kg/m^3]
     .  rtm,   ! Total water mixing ratio               [kg/kg]
     .  p,     ! Pressure                               [Pa]
     .  thvm,  ! Virtual potential temperature          [K]
     .  rcm    ! Cloud water mixing ratio               [kg/kg]


        ! Output Variables
        real, dimension(gr%nnzp), intent(out) :: 
     .  wmt,         ! Large-scale vertical motion on t grid   [m/s]
     .  wmm,         ! Large-scale vertical motion on m grid   [m/s]
     .  theta_tndcy, ! Large-scale thlm tendency               [K/s]
     .  rt_tndcy,    ! Large-scale rtm tendency                [kg/kg/s]
     .  Frad,        ! Total radiative flux                    [W/m^2]
     .  Ncm,         ! Cloud droplet number concentration      [count/m^3]
     .  Ncnm,        ! Cloud nuclei number concentration       [count/m^3]
     .  radht        ! dT/dt, then d Theta/dt, due to rad.     [K/s]


        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclr_tndcy ! Passive scalar LS tendency            [units/s]


        ! Local Variables, radiation scheme
        real, dimension(gr%nnzp) :: 
     .  radht_LW,! dT/dt, then d Theta/dt, due to LW rad.  [K/s]
     .  radht_SW,! dT/dt, then d Theta/dt, due to SW rad.  [K/s]
     .  Frad_LW, ! Longwave radiative flux                 [W/m^2]
     .  Frad_SW  ! Shortwave radiative flux                [W/m^2]


        ! Local Variables, general
        integer :: i, k ! Loop indices


        ! Local Variables, subsidence scheme
        real ::
     .  velocity_omega


        ! Local Variables, radiation scheme
        real ::
     .  xi_abs,
     .  sda_t,
     .  sda_delta,
     .  sda_h,
     .  t_since_noon,
     .  julday,
     .  start_time_until_noon,
     .  t_tendency

        real, dimension(gr%nnzp) ::
     .  radht_theta,
     .  radht_LW_theta,
     .  radht_SW_theta,
     .  LWP,
     .  thm_rad,
     .  rcm_rad,
     .  rhot_rad,
     .  dsigm,
     .  frad_out,
     .  frad_lw_out,
     .  frad_sw_out,
     .  radhtk,
     .  radht_lw_out,
     .  radht_sw_out

        ! Local variables, simplified ice (icedfs) microphysics
        real, dimension(gr%nnzp) :: dqc_dt_icedfs

        ! Local variables, on/off switches for individual schemes
        logical :: 
     .  subs_on,
     .  lw_on,
     .  sw_on,
     .  center


        ! Set which schemes to use
        subs_on         = .TRUE.
        lw_on           = .TRUE.
        sw_on           = .TRUE.
        center          = .TRUE.

        ! Compute vertical motion
        do i=2,gr%nnzp
          velocity_omega = min( D*(psfc-p(i)), D*(psfc-pinv) )
          wmt(i) = -velocity_omega * Rd * thvm(i) / p(i) / grav0
        end do
        wmm = zt2zm( wmt )


        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        

        ! Compute large-scale tendencies
        do i=1,gr%nnzp
         t_tendency = min( -4.,-15.*(1.-((psfc-p(i))/21818.)) ) ! K/day
         theta_tndcy(i) = (t_tendency * ((psfc/p(i)) ** (Rd/Cp))) 
     .                    / 86400. ! K/s
         rt_tndcy(i)  = min( 0.164,-3*(1-((psfc-p(i))/15171.)) ) / 
     .                 1000. / 86400. ! g/kg/day -> kg/kg/s
        end do

        ! Compute radiation
        julday = 282
        start_time_until_noon = 18840 + 61200
        t_since_noon   = time - start_time_until_noon
        sda_t = 2*3.14*(julday-1)/365

        sda_delta = c_0 + c_1*cos(sda_t) + d_1*sin(sda_t) +
     .              c_2*cos(2*sda_t) + d_2*sin(2*sda_t) +
     .              c_3*cos(3*sda_t) + d_3*sin(3*sda_t)

        sda_h = 2*3.14*t_since_noon/86400

        xi_abs = sin(rlat*3.14/180) * sin(sda_delta) +
     .           cos(rlat*3.14/180) * cos(sda_delta) * cos(sda_h)

        xi_abs = max(xi_abs,0.)

        if (xi_abs == 0.) then
          sw_on = .FALSE.
        end if

        if (.not. sw_on) then
          xi_abs = 0.
        end if

        if ( .not. lbugsrad ) then
          do k = 1, gr%nnzp
            thm_rad(k)  = thlm(gr%nnzp-k+1) + Lv/Cp*rcm(gr%nnzp-k+1)
            rcm_rad(k)  = rcm(gr%nnzp-k+1)
            rhot_rad(k) = rhot(gr%nnzp-k+1)
            dsigm(k)    = 1.0 / gr%dzt(gr%nnzp-k+1)
          end do

          call rad_lwsw(thm_rad, rcm_rad, rhot_rad, dsigm,
     .                  Frad_out, Frad_LW_out, Frad_SW_out,
     .                  radhtk, radht_LW_out, radht_SW_out,
     .                  gr%nnzp-1, center,
     .                  xi_abs, F0, F1, kap, radius, A, gc, Fs0, omega,
     .                  sw_on, lw_on)

          do k = 1, gr%nnzp-1
            Frad(k)     = Frad_out(gr%nnzp-k+1)
            Frad_LW(k)  = Frad_LW_out(gr%nnzp-k+1)
            Frad_SW(k)  = Frad_SW_out(gr%nnzp-k+1)

            radht(k)    = radhtk(gr%nnzp-k+1)
            radht_LW(k) = radht_LW_out(gr%nnzp-k+1)
            radht_SW(k) = radht_SW_out(gr%nnzp-k+1)

            radht_theta(k)    = radht(k) * ((p00/p(k))**(Rd/Cp))
            radht_LW_theta(k) = radht_LW(k) * ((p00/p(k))**(Rd/Cp))
            radht_SW_theta(k) = radht_SW(k) * ((p00/p(k))**(Rd/Cp))
          end do ! k

          Frad(1)    = Frad(2)
          Frad_LW(1) = Frad_LW(2)
          Frad_SW(1) = Frad_SW(2)
          radht_theta(1)    = radht_theta(2)
          radht_LW_theta(1) = radht_LW_theta(2)
          radht_SW_theta(1) = radht_SW_theta(2)

          Frad(gr%nnzp)    = Frad(gr%nnzp-1)
          Frad_LW(gr%nnzp) = Frad_LW(gr%nnzp-1)
          Frad_SW(gr%nnzp) = Frad_SW(gr%nnzp-1)
          radht_theta(gr%nnzp)    = radht_theta(gr%nnzp-1)
          radht_LW_theta(gr%nnzp) = radht_LW_theta(gr%nnzp-1)
          radht_SW_theta(gr%nnzp) = radht_SW_theta(gr%nnzp-1)

          radht(1:gr%nnzp)    = radht_theta(1:gr%nnzp)
          radht_LW(1:gr%nnzp) = radht_LW_theta(1:gr%nnzp)
          radht_SW(1:gr%nnzp) = radht_SW_theta(1:gr%nnzp)

          do k = 1, gr%nnzp
            theta_tndcy(k) = theta_tndcy(k) + radht_theta(k)
          end do

        end if ! ~ lbugsrad

#ifdef STATS
        if ( .not.lbugsrad .and. lstats_samp ) then

          if ( iradht_LW > 0 ) then
            zt%x(:,iradht_LW) = zt%x(:,iradht_LW) + radht_LW
            zt%n(:,iradht_LW) = zt%n(:,iradht_LW) + 1
          end if
          if ( iradht_SW > 0 ) then
            zt%x(:,iradht_SW) = zt%x(:,iradht_SW) + radht_SW
            zt%n(:,iradht_SW) = zt%n(:,iradht_SW) + 1
          end if

          if ( iFrad_SW > 0 ) then
            zm%x(:,iFrad_SW) = zm%x(:,iFrad_SW) + Frad_SW
            zm%n(:,iFrad_SW) = zm%n(:,iFrad_SW) + 1
          end if
          if ( iFrad_LW > 0 ) then
            zm%x(:,iFrad_LW) = zm%x(:,iFrad_LW) + Frad_LW
            zm%n(:,iFrad_LW) = zm%n(:,iFrad_LW) + 1
          end if

        end if
#endif /*STATS*/


        ! Compute icedfs microphysics
        if ( licedfs ) then
          call altocu_icedf( dt, thlm, rcm, p, exner, rhot,
     .                       dqc_dt_icedfs )
          do k=1,gr%nnzp
            theta_tndcy(k) = theta_tndcy(k)
     .                       - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
            rt_tndcy(k) = rt_tndcy(k) + dqc_dt_icedfs(k)
          end do
        end if

        ! Initialize Ncnm on first timestep
        if ( lcoamps_micro .and. time == time_initial ) then
          Ncnm(1:gr%nnzp) 
     .    = 30.0 * (1.0 + exp(-gr%zt(1:gr%nnzp)/2000.0)) * 1.e6

        else if ( kk_rain ) then
          do k=1, gr%nnzp, 1
            if ( rcm(k) >= 1.e-6 ) then
              Ncm(k) = 30.0 * (1.0 + exp(-gr%zt(k)/2000.0)) * 1.e6
     .                 * rhot(k) 
            end if
          end do
        end if


        ! Set scalar tendencies
        if ( present( sclr_tndcy ) ) then
          sclr_tndcy(:,SCLR_THETA) = theta_tndcy
          sclr_tndcy(:,SCLR_RT)    = rt_tndcy
        end if

        return
        end subroutine mpace_tndcy
!----------------------------------------------------------------------

!----------------------------------------------------------------------
        subroutine mpace_sfclyr( rho0, up, vp, uw, vw,
     .                          wtp, wrp, wpsclrp_sfc )

!----------------------------------------------------------------------
!        Description:
!          Surface forcing subroutine for MPACE case.  Written July-
!          November 2006 by Michael Falk.
!
!        References:
!          MPACE specification, arm.gov
!-----------------------------------------------------------------------

        use constants
        implicit none

        ! Input Variables
        real, intent(in)  ::
     .   rho0,
     .   up,
     .   vp

        ! Output Variables
        real, intent(out) ::
     .    uw,
     .    vw,
     .    wtp,
     .    wrp

        ! Output Variables (optional) 
        real, optional, intent(out) ::
     .    wpsclrp_sfc(sclr_dim)

        ! Internal Variables
        real ::
     .    ubar

        real, parameter ::
     .    ubmin = 0.25,
     .    ustar = 0.25,
     .    HF    = 136.5, ! Sensible Heat Flux is 136.5 W m^-2, Specified by MPACE specification
     .    LF    = 107.7  ! Latent Heat Flux is 107.7 W m^-2, Specified by MPACE specification

!-----------------------------------------------------------------------

        ! Compute heat and moisture fluxes
        wtp = HF/(rho0*Cp)
        wrp = LF/(rho0*Lv)

        ! Compute momentum fluxes
        ubar = max(ubmin,sqrt(up*up + vp*vp))
        uw = -up * ustar*ustar / ubar
        vw = -vp * ustar*ustar / ubar

        ! Compute scalar fluxes
        if (present( wpsclrp_sfc )) wpsclrp_sfc(SCLR_THETA) = wtp
        if (present( wpsclrp_sfc )) wpsclrp_sfc(SCLR_RT) = wrp

        return
        end subroutine mpace_sfclyr

!----------------------------------------------------------------------

        subroutine cobra_tndcy( time, wmt, wmm, theta_tndcy, 
     .                          rt_tndcy,
     .                          sclr_tndcy )
!       Description:
!       Subroutine to set theta and water tendencies for COBRA CO2 case
!       References;
!       None
        use grid_class
        use constants, only: sclr_dim

        implicit none

        ! Input Variables
        real, intent(in) :: time ! Model time [s]

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) ::
     .  wmt,         ! w wind on thermodynamic grid [m/s]
     .  wmm,         ! w wind on momentum grid [m/s]
     .  theta_tndcy, ! Liquid water potential temperature tendency [K/s]
     .  rt_tndcy     ! Total water mixing ratio tendency [kg/kg/s]

        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclr_tndcy ! Passive scalar tendency [units vary]

!       No large-scale subsidence for now
        wmt = 0.0
        wmm = 0.0

!       No large-scale water tendency or cooling

        rt_tndcy    = 0.0
        theta_tndcy = 0.0

        if ( present( sclr_tndcy ) ) then
!          sclr_tndcy(1:gr%nnzp,SCLR_THETA) = theta_tndcy(1:gr%nnzp)
          sclr_tndcy(1:gr%nnzp,SCLR_THETA) = 0.0
          sclr_tndcy(1:gr%nnzp,SCLR_RT)    = rt_tndcy(1:gr%nnzp)
        end if

        return
        end subroutine cobra_tndcy

!-----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to the format used for the GCSS ARM case.
!
! The data has been altered so it can be used for the COBRA CO2 case.

        subroutine cobra_sfclyr( time, z, dn0, thp, up, vp,
     .                           uw, vw, wtp, wrp, ustar,
     .                           wpsclrp_sfc, wpedsclrp_sfc )

        use constants, only: Cp, Lv, grav, sclr_dim
        implicit none

        ! External
        ! In module calls:
        ! real diag_ustar

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
        real, optional, intent(out) :: wpedsclrp_sfc(sclr_dim)

!       Internal variables

        integer i, j
        real usfc, vsfc, ubar
        real true_time
        real heat_flx, moisture_flx
        real heat_flx2, moisture_flx2
        real bflx

        real CO2_flx, CO2_flx2
        real, parameter :: M_da = 0.02897  ! Molecular weight of dry air.

!       COBRA roughness height
!        real, parameter :: z0 = 0.035  ! ARM momentum roughness height
        real, parameter :: z0 = 1.75   ! momentum roughness height

! Compute heat and moisture fluxes from ARM data in (W/m2)

        true_time = time
        call cobra_sfcflx( true_time, heat_flx, moisture_flx
#ifdef SCALARS
     .                    ,CO2_flx
#endif /* SCALARS */
     .                   )

! Compute momentum fluxes

!       Convert heat_flx and moisture_flx to natural units
        heat_flx2     = heat_flx / ( Cp * dn0 )    ! (K m/s)
        moisture_flx2 = moisture_flx / ( Lv * dn0 )! (m/s)

#ifdef SCALARS
!       Convert CO2 surface flux to natural units.
!       The CO2 flux has been given in units of:  umol/(m^2 s).
!       umol stands for micromoles.  The CO2 concentration in
!       this code is in units of ppmv, which is also the molar
!       mixing ratio times 10^6.
!       The units are:  10^6 * [ mol (CO2) / mol (dry air) ].
!       w'CO2' = (Flux) * [ M (dry air) / rho (dry air) ];
!       where M is the molecular weight of dry air.
        CO2_flx2 = CO2_flx * ( M_da / dn0 )
#endif /* SCALARS */

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

!        if (present( wpsclrp_sfc ))
!     .                 wpsclrp_sfc(SCLR_THETA) = wtp
        if (present( wpsclrp_sfc ))
     .                 wpsclrp_sfc(SCLR_THETA) = CO2_flx2
        if (present( wpsclrp_sfc ))
     .                 wpsclrp_sfc(SCLR_RT) = wrp

!        if (present( wpedsclrp_sfc ))
!     .                 wpedsclrp_sfc(SCLR_THETA) = wtp
        if (present( wpedsclrp_sfc ))
     .                 wpedsclrp_sfc(SCLR_THETA) = CO2_flx2
        if (present( wpedsclrp_sfc ))
     .                 wpedsclrp_sfc(SCLR_RT) = wrp

        return
        end subroutine cobra_sfclyr

!-----------------------------------------------------------------------
! This subroutine computes surface heat and moisture for a specific time
! according to the format used for the GCSS ARM case.
! Fluxes returned are in (W/m2).
!
! The data has been altered so it can be used for the COBRA CO2 case.

        subroutine cobra_sfcflx( time, heat_flx, moisture_flx,
     .                           CO2_flx )

        implicit none

!       Input variable
        real, intent(in) :: time       ! time in seconds

!       Output variables
        real, intent(out) :: heat_flx, moisture_flx

!       Output variables (optional)
        real, optional, intent(out) :: CO2_flx

!       Internal variables
        integer i1,i2
        real a

        integer ntimes
        parameter ( ntimes = 6 )
        real times(ntimes), H(ntimes), LE(ntimes)
#ifdef SCALARS
        real CO2(ntimes)
#endif /* SCALARS */
        ! Best fit line
!        data times / 57600.0, 64800.0, 73800.0, 75600.0 /
        ! Every hour (3 half-hour point weighted average)
        data times / 57600.0, 61200.0, 64800.0, 68400.0,
     .               72000.0, 75600.0 /

!       H and LE specifications
        ! Every hour (3 half-hour point weighted average)
        data H   / 415.3, 388.3, 404.0, 326.8, 
     .             236.9, 220.1 /
        data LE  / 121.2, 119.1, 122.9, 102.7,
     .             104.8,  87.1 /

#ifdef SCALARS
!       CO2 specifications
        ! Every hour (3 half-hour point weighted average)
        data CO2 / -12.4, -11.3, -10.4,  -8.0,
     .              -7.2,  -7.1 /
#endif /* SCALARS */

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
#ifdef SCALARS
                 CO2_flx      = ( 1. - a ) * CO2(i1) + a * CO2(i2)
#endif /* SCALARS */
                 i1           = ntimes
              end if
              i1 = i2
           end do
        end if

        return
        end subroutine cobra_sfcflx





!----------------------------------------------------------------------
        subroutine rico_tndcy
     .  ( dt, exner, wmt, wmm, theta_tndcy,
     .    rt_tndcy, radht,
     .    sclr_tndcy )

!        Description:
!          Subroutine to apply case-specific forcings to RICO case
!          (Michael Falk, 13 Dec 2006).
!
!        References:
!-----------------------------------------------------------------------

        use constants
        use grid_class
#ifdef STATS
        use statistics
#endif /*STATS*/
        implicit none

        ! Input Variables
        real, intent(in) :: 
     .  dt            ! Current length of timestep      [s]

        real, dimension(gr%nnzp), intent(in) ::
     .  exner        ! Exner function                         [-]

        ! Output Variables
        real, dimension(gr%nnzp), intent(out) :: 
     .  wmt,         ! Large-scale vertical motion on t grid   [m/s]
     .  wmm,         ! Large-scale vertical motion on m grid   [m/s]
     .  theta_tndcy, ! Large-scale thlm tendency               [K/s]
     .  rt_tndcy,    ! Large-scale rtm tendency                [kg/kg/s]
     .  radht        ! dT/dt, then d Theta/dt, due to rad.     [K/s]

        ! Output Variables (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclr_dim) ::
     .  sclr_tndcy ! Passive scalar LS tendency            [units/s]

        ! Local Variables, general
        integer :: i, k ! Loop indices
        real    :: t_tendency ! Temperature (not potential temperature) tendency

        ! Compute vertical motion
        do k=1,gr%nnzp
          if (gr%zt(k) < 2100) then
            wmt(k) = -(0.005 / 2100) * gr%zt(k)
          else if (gr%zt(k) < 3000) then
            wmt(k) = -0.005 + (0.005 / (3000 - 2100))
     .               * (gr%zt(k) - 2100)
          else
            wmt(k) = 0.
          end if
        end do
        wmm = zt2zm( wmt )

        ! Boundary conditions on vertical motion.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        

        ! Compute large-scale horizontal temperature advection
        do k=1,gr%nnzp
          if (gr%zt(k) < 2100) then
            t_tendency = - 0.8 / 86400  ! Units [K/s]
          else if (gr%zt(k) < 3000) then
            t_tendency = - 0.8 / 86400
     .                   + (0.8 / 86400) / (3000 - 2100)
     .                   * (gr%zt(k) - 2100)  ! Units [K/s]
          else
            t_tendency = 0.  ! Units [K/s]
          end if
          ! Convert to units of [K/s] but potential T instead of T
!          theta_tndcy(k) = (t_tendency * ((psfc/p(k)) ** (Rd/Cp)))
          theta_tndcy(k) = (t_tendency / exner(k))
        end do


        ! Compute large-scale horizontal moisture advection [g/kg/s]
        do k=1,gr%nnzp
          if (gr%zt(k) < 550) then
            rt_tndcy(k) = - 1.0 / 86400  ! Units [g/kg/s]
          else if (gr%zt(k) < 1500) then
            rt_tndcy(k) = - 1.0 / 86400
     .                    + ((0.0 + 1.0) / 86400)
     .                    / (1500 - 550)
     .                   * (gr%zt(k) - 550)  ! Units [g/kg/s]
          else
            rt_tndcy(k) = 0.  ! Units [K/s]
          end if
          rt_tndcy(k) = rt_tndcy(k) / 1000.  ! Converts [g/kg/s] to [kg/kg/s]
        end do


        ! Compute radiative cooling
        do k=1,gr%nnzp
          if (gr%zt(k) < 2100) then
            radht(k) = - (2 / 86400)
     .              + ((-1.3 + 2) / 86400) / 2100 * gr%zt(k)  ! Units [K/s]
          else if (gr%zt(k) < 3000) then
            radht(k) = - (1.3 / 86400)
     .              + (1.3 / 86400) / (3000 - 2100) * (gr%zt(k) - 2100)  ! Units [K/s]
          else
            radht(k) = 0.
          end if
!          radht(k) = (radht(k) * ((psfc/p(k)) ** (Rd/Cp)))
          radht(k) = (radht(k) / exner(k))
          theta_tndcy(k) = theta_tndcy(k) + radht(k)
        end do


        ! Set scalar tendencies
        if ( present( sclr_tndcy ) ) then
          sclr_tndcy(:,SCLR_THETA) = theta_tndcy
          sclr_tndcy(:,SCLR_RT)    = rt_tndcy
        end if

        return
        end subroutine rico_tndcy
!----------------------------------------------------------------------




!----------------------------------------------------------------------
        subroutine rico_sfclyr( up, vp, thp, rtp,
     .                          lowestlevel, sst, psfc,
     .                          uw, vw, wtp, wrp,
     .                          sclrm_sfc, wpsclrp_sfc )

!----------------------------------------------------------------------
!        Description:
!          Surface forcing subroutine for MPACE case.  Written July-
!          November 2006 by Michael Falk.
!
!        References:
!          MPACE specification, arm.gov
!-----------------------------------------------------------------------

        use constants
        implicit none

!       External
        real, external :: rsat

!       Constants
        real, parameter :: ubmin = 0.25    ! I don't know where this comes from.  All of the other cases use it, though.
        real, parameter :: ustar = 0.3     ! Defined by ATEX specification
        real, parameter :: C_10  = 0.0013  ! Defined by ATEX specification
        real, parameter :: z0    = 0.00015 ! Roughness length, defined by ATEX specification

!       Internal variables
        real ubar,Cz

!       Input variables
        real, intent(in) :: up,vp,thp,rtp,lowestlevel,sst,psfc

!       Input variables (optional)
        real, optional, intent(in) :: sclrm_sfc(sclr_dim)

!       Output variables
        real, intent(out) :: uw, vw, wtp, wrp

!       Output variables (optional)
        real, optional, intent(out) :: wpsclrp_sfc(sclr_dim)

! Compute heat and moisture fluxes

        ubar = max(ubmin, sqrt(up*up + vp*vp))
        Cz   = C_10 * ((log(10/z0))/(log(lowestlevel/z0))) * 
     .         ((log(10/z0))/(log(lowestlevel/z0)))              ! Modification in case lowest model level isn't at 10 m, from ATEX specification
        wtp  = -Cz * ubar * ( thp - sst * (p0/psfc)**kappa )     ! K * m/s
        wrp  = -Cz * ubar * ( rtp - rsat(psfc,sst) )             ! kg/kg * m/s

! Optional passive scalar
        if ( present( wpsclrp_sfc ) .and. present( sclrm_sfc ) ) then
          wpsclrp_sfc(SCLR_THETA) = -C_10 * ubar * 
     .           ( sclrm_sfc(SCLR_THETA) - sst * (p0/psfc)**kappa )

          wpsclrp_sfc(SCLR_RT) = -C_10 * ubar * 
     .                     ( sclrm_sfc(SCLR_RT) - rsat(psfc,sst) )
        end if

! Compute momentum fluxes

        uw = -up * ustar*ustar / ubar
        vw = -vp * ustar*ustar / ubar

        ! Compute scalar fluxes
        if (present( wpsclrp_sfc )) wpsclrp_sfc(SCLR_THETA) = wtp
        if (present( wpsclrp_sfc )) wpsclrp_sfc(SCLR_RT) = wrp

        return
        end subroutine rico_sfclyr

!----------------------------------------------------------------------




!-----------------------------------------------------------------------

        FUNCTION Diff_denom(Temp, press, ei)

        USE constants

        IMPLICIT NONE

        REAL, INTENT(IN):: Temp, press, ei
        REAL:: Diff_denom

        REAL:: Ka, Dv
        REAL:: Fk, Fd
        REAL:: Celsius

!        REAL, PARAMETER:: Ls = 2.834e6

        Celsius = Temp - 273.16

        Ka = (5.69 + 0.017*Celsius)*0.00001  ! Ka in cal./(cm.*sec.*C)
        Ka = 4.1868*100.0*Ka  ! Ka in J./(m.*sec.*K)

        Dv = 0.221*((Temp/273.16)**1.94)*(101325.0/press)
                                ! Dv in (cm.^2)/sec.  ! .221 is correct.
        Dv = Dv/10000.0  ! Dv in (m.^2)/sec.

        Fk = (Ls/(Rv*Temp) - 1.0) * Ls/(Ka*Temp)
        Fd = (Rv*Temp)/(Dv*ei)

        Diff_denom = Fk + Fd

        RETURN
        END FUNCTION Diff_denom

      end module gcss
!------------------------------------------------------------------------
