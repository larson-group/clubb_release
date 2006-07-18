!----------------------------------------------------------------------
! $Id: gcss.F,v 1.5 2006-07-18 16:54:40 dschanen Exp $ 

! This file contains extra subroutines to simulate the GCSS ARM and
! GCSS BOMEX, GCSS ATEX, GCSS FIRE intercomparison cases as well 
! as Wangara day 33
!----------------------------------------------------------------------
#define SCLR_THETA 1
#define SCLR_RT 2
      module gcss
      use model_flags

      implicit none

      public cloud_rad, atex_tndcy, atex_sfclyr, fire_tndcy
      public wangara_tndcy, wangara_sfclyr
      public astex_tndcy, astex_sfclyr
      public arm_tndcy, arm_sfclyr
      public bomex_tndcy, bomex_sfclyr
      public dycoms2_rf01_tndcy, dycoms2_rf01_sfclyr 
      public dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr
      public nov11_altocu_tndcy 

      private diag_ustar, arm_sfcflx, Diff_denom
      private altocu_icedf

      contains
!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for ARM case
!----------------------------------------------------------------------
        subroutine arm_tndcy( time, theta_tndcy, rt_tndcy,
     .                        sclr_tndcy )

        use grid_class
        use constants, only: sclrm_dimension

        implicit none

!       Input
        real, intent(in)  :: time

!       Output
        real, intent(out) :: theta_tndcy(gr%nnzp)
        real, intent(out) :: rt_tndcy(gr%nnzp)

!       Output (optional)
        real, optional, intent(out), 
     .  dimension(gr%nnzp,sclrm_dimension) :: sclr_tndcy

!       Internal variables
        integer i, i1, i2
        real a, b
        real true_time
        real theta_tmp, rt_tmp

        real atheta(6), rtheta(6), art(6)
        data atheta / 0.000, 0.000, 0.000, -0.080, -0.160, -0.160 /
        data rtheta / -0.125, 0.000, 0.000, 0.000, 0.000, -0.100 /
        data art    / 0.080, 0.020, -0.040, -0.100, -0.160, -0.300 /

        true_time = time

!       Interpolate in time to get theta and rt tendency

        i1 = int( ( true_time - 41400. ) / 10800. ) + 1
        i1 = min( max( i1, 1 ), 5 )
        i2 = i1 + 1

        if (i1 < 5) then
           a = ( true_time - (41400. + 10800. * (i1-1)) ) / 10800.
        else
           a = ( true_time - (41400. + 10800. * (i1-1)) ) / 9000.
        endif

        if ( .not. lbugsrad ) then
          theta_tmp = ( 1. - a ) * ( atheta(i1) + rtheta(i1) )
     .              + a * ( atheta(i2) + rtheta(i2) )
        else ! factor in radiation later
          theta_tmp = ( 1. - a ) * ( atheta(i1) + 0.0 )
     .              + a * ( atheta(i2) + 0.0 )
        endif

        rt_tmp    = ( 1. - a ) * art(i1) + a * art(i2)

!       Convert to the right units

        theta_tmp = theta_tmp / 3600.
        rt_tmp    = rt_tmp / ( 3600. * 1000. )

!       Interpolate with respect to height

        do i=2, gr%nnzp
          if ( gr%zt(i) >= 0. .and. gr%zt(i) < 1000. ) then
            rt_tndcy(i)    = rt_tmp
            theta_tndcy(i) = theta_tmp
          elseif ( gr%zt(i) >= 1000. .and. gr%zt(i) < 3000. ) then
            b              = 1. - ( gr%zt(i) - 1000. ) / 2000.
            rt_tndcy(i)    = b * rt_tmp
            theta_tndcy(i) = b * theta_tmp
          else
            rt_tndcy(i)    = 0.
            theta_tndcy(i) = 0.
          endif
        enddo

        if (present( sclr_tndcy )) sclr_tndcy(:,SCLR_THETA)= theta_tndcy
        if (present( sclr_tndcy )) sclr_tndcy(:,SCLR_RT)   = rt_tndcy

        return
        end subroutine arm_tndcy

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for BOMEX case
!----------------------------------------------------------------------
        subroutine bomex_tndcy( time, wmt, wmm, theta_tndcy, rt_tndcy,
     .                          sclr_tndcy )

        use grid_class
        use constants, only: sclrm_dimension

        implicit none

!       Input
        real, intent(in)  :: time

!       Output
        real, intent(out) :: wmt(gr%nnzp)
        real, intent(out) :: wmm(gr%nnzp)
        real, intent(out) :: theta_tndcy(gr%nnzp)
        real, intent(out) :: rt_tndcy(gr%nnzp)

!       Output (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclrm_dimension)
     .                    :: sclr_tndcy

!       Internal variables

        integer i

!       Large scale subsidence

        do i=2, gr%nnzp

           if ( gr%zt(i) >= 0. .and. gr%zt(i) < 1500. ) then
              wmt(i) = - ( 0.0065 / 1500. ) * gr%zt(i)
           elseif ( gr%zt(i) >= 1500. .and. gr%zt(i) < 2100. ) then
              wmt(i)
     .          = - 0.0065 
     .            + 0.0065 * ( gr%zt(i) - 1500. ) / ( 2100. - 1500. )
           else
              wmt(i) = 0.
           endif

        enddo

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
!       Radiative theta-l tendency

        do i=2,gr%nnzp

           if ( gr%zt(i) >= 0. .and. gr%zt(i) < 1500. ) then
              theta_tndcy(i) = -2.315e-5
           elseif ( gr%zt(i) >= 1500. .and. gr%zt(i) < 2500. ) then
              theta_tndcy(i)
     .          = - 2.315e-5 
     .            + 2.315e-5 
     .              * ( gr%zt(i) - 1500. ) / ( 2500. - 1500. )
           else
              theta_tndcy(i) = 0.
           endif

        enddo

!       Large scale advective moisture tendency

        do i=2,gr%nnzp

           if ( gr%zt(i) >= 0. .and. gr%zt(i) < 300. ) then
              rt_tndcy(i) = -1.2e-8
           elseif ( gr%zt(i) >= 300. .and. gr%zt(i) < 500. ) then
              rt_tndcy(i) 
     .          = - 1.2e-8 
     .              * ( 1. - ( gr%zt(i) - 300. )/( 500. - 300. ) )
           else
              rt_tndcy(i) = 0.
           endif

        enddo

        ! Boundary conditions
        theta_tndcy(1) = 0.0  ! Below surface
        rt_tndcy(1)    = 0.0  ! Below surface

        if (present( sclr_tndcy )) sclr_tndcy(:,SCLR_THETA)=theta_tndcy
        if (present( sclr_tndcy )) sclr_tndcy(:,SCLR_RT) = rt_tndcy

        return
        end subroutine bomex_tndcy

!----------------------------------------------------------------------
! Subroutine to large-scale subsidence for FIRE case
!----------------------------------------------------------------------
        subroutine fire_tndcy( time, wmt, wmm, theta_tndcy, rt_tndcy,
     .                         sclr_tndcy )

        use constants, only: sclrm_dimension
        use grid_class

        implicit none

!       Input
        real, intent(in) :: time

!       Output
        real, intent(out) :: wmt(gr%nnzp)
        real, intent(out) :: wmm(gr%nnzp)
        real, intent(out) :: theta_tndcy(gr%nnzp)
        real, intent(out) :: rt_tndcy(gr%nnzp)

!       Output (optional)
        real, optional, intent(out), dimension(gr%nnzp,sclrm_dimension)
     .                    :: sclr_tndcy

!       Internal variables

        integer i

!       Large-scale subsidence

        do i=2,gr%nnzp

           if ( gr%zt(i) >= 0. .and. gr%zt(i) < 1500. ) then
              wmt(i) = - 5.e-6 * gr%zt(i)
           endif

        enddo

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
!       Radiative theta-l tendency is computed interactively elsewhere

        theta_tndcy = 0.

!       Large scale advective moisture tendency

        rt_tndcy    = 0.

        if ( present( sclr_tndcy ) ) 
     .    sclr_tndcy(:,SCLR_THETA) = theta_tndcy(:)
        if ( present( sclr_tndcy ) ) sclr_tndcy(:,SCLR_RT) = rt_tndcy(:)

        return
        end subroutine fire_tndcy

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for Wangara case
!
        subroutine wangara_tndcy( time, wmt, wmm, theta_tndcy, 
     .                            rt_tndcy,
     .                            sclr_tndcy )
        use grid_class
        use constants, only: sclrm_dimension

        implicit none

!       Input
        real, intent(in) :: time

!       Output
        real, dimension(gr%nnzp), intent(out) :: wmt
        real, dimension(gr%nnzp), intent(out) :: wmm
        real, dimension(gr%nnzp), intent(out) :: theta_tndcy
        real, dimension(gr%nnzp), intent(out) :: rt_tndcy

!       Optional Output
        real, optional, dimension(gr%nnzp, sclrm_dimension), 
     .        intent(out) :: sclr_tndcy

!       No large-scale subsidence for now
        wmt = 0.
        wmm = 0.

!       No large-scale water tendency or cooling

        rt_tndcy    = 0.
        theta_tndcy = 0.

!       Scalar r forcing is just set to be = to rtm_forcing for now
        if ( present( sclr_tndcy ) ) 
     .    sclr_tndcy(:,SCLR_THETA) = theta_tndcy
        if ( present( sclr_tndcy ) ) sclr_tndcy(:,SCLR_RT) = rt_tndcy

        return
        end subroutine wangara_tndcy

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for ATEX case
!
        subroutine atex_tndcy( time, time_initial, rtm, wmt, wmm,
     .                         theta_tndcy, rt_tndcy, isValid,
     .                         sclr_tndcy )

        use constants, only: sclrm_dimension
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
        real, optional, intent(out), dimension(gr%nnzp, sclrm_dimension)
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

           do i=2,gr%nnzp

              if ( gr%zt(i) > 0. .and. gr%zt(i) <= zi ) then
                 wmt(i) 
     .             = -0.0065 * gr%zt(i)/zi
              elseif ( gr%zt(i) > zi .and. gr%zt(i) <= zi+300. ) then
                 wmt(i)
     .             = - 0.0065 * ( 1. - (gr%zt(i)-zi)/300. )
              else
                 wmt(i) = 0.
              endif

           end do

           wmm = zt2zm( wmt )

           ! Boundary conditions.
           wmt(1) = 0.0        ! Below surface
           wmm(1) = 0.0        ! At surface
           wmm(gr%nnzp) = 0.0  ! Model top
        
!          theta-l tendency

           do i=2,gr%nnzp

              if ( gr%zt(i) > 0. .and. gr%zt(i) < zi ) then
                 theta_tndcy(i) 
     .             = -1.1575e-5 * ( 3. - gr%zt(i)/zi )
              elseif ( gr%zt(i) > zi .and. gr%zt(i) <= zi+300. ) then
                 theta_tndcy(i)
     .             = -2.315e-5 * ( 1. - (gr%zt(i)-zi)/300. )
              else
                 theta_tndcy(i) = 0.
              endif

           end do

!          moisture tendency

           do i=2,gr%nnzp

              if ( gr%zt(i) > 0. .and. gr%zt(i) < zi ) then
                 rt_tndcy(i) = -1.58e-8 * ( 1. - gr%zt(i)/zi )  ! Brian
              else
                 rt_tndcy(i) = 0.0       ! Brian
              endif

           end do

           ! Boundary conditions
           theta_tndcy(1) = 0.0  ! Below surface
           rt_tndcy(1)    = 0.0  ! Below surface

        end if

        if ( present( sclr_tndcy ) ) 
     .    sclr_tndcy(:,SCLR_THETA) = theta_tndcy
        if ( present( sclr_tndcy ) ) sclr_tndcy(:,SCLR_RT) = rt_tndcy

        return
        end subroutine atex_tndcy

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to GCSS ARM specifications
!
        subroutine arm_sfclyr( time, z, dn0, thp, up, vp, uw, vw, 
     .                         wtp, wrp, ustar,
     .                         wpsclrp_sfc )

        use constants
        implicit none

!       External (not needed with module)

!        real diag_ustar
!        external diag_ustar

!       Constants (now in constants.F)

!        real Cp, Lv
!        parameter ( Cp = 1004.67 )
!        parameter ( Lv = 2.5e6 )
!        real g
!        parameter ( g = 9.81 )
        real ubmin
        parameter (ubmin = 0.25)

!       Input variables

        real, intent(in) :: time
        real, intent(in) :: z
        real, intent(in) :: dn0
        real, intent(in) :: thp
        real, intent(in) :: up,vp

        real, intent(inout) :: ustar

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp

!       Output variables (optional)
        real, optional, intent(out) :: wpsclrp_sfc(sclrm_dimension)

!       Internal variables

        integer i, j
        real usfc, vsfc, ubar
        real true_time
        real heat_flx, moisture_flx
        real heat_flx2, moisture_flx2
        real bflx

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
        ustar = diag_ustar( z, bflx, ubar )

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

        use constants, only: sclrm_dimension

        implicit none

!       Constant

        real, parameter :: ubmin = 0.25
        real, parameter :: ustar = 0.28

!       Input variables

        real, intent(in) :: up, vp

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp

!       Output variables (optional)

        real, optional, intent(out) :: wpsclrp_sfc(sclrm_dimension)

!       Internal variables

        real ubar

! Compute heat and moisture fluxes

        wtp = 8.e-3
        wrp = 5.2e-5
        
        if ( present(wpsclrp_sfc) ) then 
          wpsclrp_sfc(SCLR_THETA) = 8.e-3
          wpsclrp_sfc(SCLR_RT)    = 5.2e-5
        endif

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

        real, optional, intent(in) :: sclrm_sfc(sclrm_dimension)

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp

!       Output variables (optional)

        real, optional, intent(out) :: wpsclrp_sfc(sclrm_dimension)

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
        endif

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
        use constants, only: sclrm_dimension
        implicit none

!       Constants

        real, parameter :: ubmin = 0.25
        real, parameter :: ustar = 0.13

!       Input variables

        real, intent(in) :: time
        real, intent(in) :: up, vp

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp
        
        real, optional, intent(out), dimension(sclrm_dimension) :: 
     .        wpsclrp_sfc

!       Internal variables

        real ubar
        real time_utc,time_est

! Compute UTC time of the day in seconds

        time_utc = mod( time, 86400. )

! Now convert UTC time to Australia EST (local time)

        time_est = mod( time_utc + 36000., 86400. )

        if ( time_est < 27000 .or. time_est > 63000 ) then
           write(*,*)
     .   'wangara_sfclyr: error local time must be between 730 and 1730'
           write(*,*) 'time_est = ',time_est
           stop
        endif

! Compute heat and moisture fluxes

        wtp = 0.18
     .        * cos( (time_est-45000.)/36000. * 3.141592654 )
        wrp = 1.3e-4 * wtp

! Compute new mixing variables
        if ( present(wpsclrp_sfc) ) wpsclrp_sfc(SCLR_THETA) = 0.18
     .        * cos( (time_est-45000.)/36000. * 3.141592654 )

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
              endif
              i1 = i2
           enddo
        endif

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
      real function diag_ustar( z, bflx, wnd ) 

      use constants, only: eps, grav, vonk

      implicit none

!     Added to constants.F
!      real, parameter      :: vonk =  0.4   ! von Karmans constant
!      real, parameter      :: g    = 9.81   ! gravitational acceleration
!      real, parameter      :: eps  = 1.e-10 ! non-zero, small number
      real, parameter      :: z0   = 0.035  ! momentum roughness height
      real, parameter      :: am   =  4.8   !   "          "         "
      real, parameter      :: bm   = 19.3   !   "          "         "

      real, intent (in)    :: z             ! height where u locates
      real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
      real, intent (in)    :: wnd           ! wind speed at z

      integer :: iterate
      real    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

      lnz   = log(z/z0)
      klnz  = vonk/lnz
      c1    = 3.14159/2. - 3.*log(2.)

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
            psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
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
          enddo

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
          endif

        endif ! ~ lbugsrad

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
#endif/*STATS*/

        implicit none

!       Input
        real, intent(in) :: time
        real, intent(in) :: rhot(gr%nnzp)
        real, intent(in) :: rhom(gr%nnzp)
        real, intent(in) :: rtm(gr%nnzp)
        real, intent(in) :: rcm(gr%nnzp)
        real, intent(in) :: exner(gr%nnzp)

        real, optional, intent(in) :: sclr(gr%nnzp, sclrm_dimension)

!       Output
        real, intent(out) :: wmt(gr%nnzp)
        real, intent(out) :: wmm(gr%nnzp)

        real, intent(out) :: Frad(gr%nnzp)
        real, intent(out) :: radht(gr%nnzp)
        real, intent(out) :: theta_tndcy(gr%nnzp)

        real, intent(out) :: rt_tndcy(gr%nnzp)

        real, optional, 
     .  intent(out) :: sclr_tndcy(gr%nnzp, sclrm_dimension)

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
#endif STATS

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
          enddo
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
             endif
          enddo

!         Compute IR heating rate

          radht          = ( -1.0/(Cp*rhot) ) * ddzm( Frad )
     .                   * ( 1.0 / exner )
          radht(1)       = 0.
          radht(gr%nnzp) = 0.

        endif ! ~ lbugsrad

!       Add heating rate to theta-l forcing

        if ( .not. lbugsrad ) theta_tndcy = theta_tndcy + radht
        
        if (present( sclr_tndcy )) sclr_tndcy(:,SCLR_THETA)= theta_tndcy
        if (present( sclr_tndcy )) sclr_tndcy(:,SCLR_RT)   = rt_tndcy

        return
        end subroutine dycoms2_rf01_tndcy

!----------------------------------------------------------------------

        SUBROUTINE dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, exner, 
     .                                 wmt, wmm, theta_tndcy, rt_tndcy, 
     .                                 Frad, radht, Ncm,
     .                                 isValid,
     .                                 sclr_tndcy )

        USE grid_class
        USE constants
#ifdef STATS
        USE statistics
#endif

        ! Input
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
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Ncm

        ! Optional output
        REAL, OPTIONAL, 
     .  DIMENSION(gr%nnzp, sclrm_dimension), INTENT(OUT):: sclr_tndcy

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
        ENDDO

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
          ENDDO

!         Find the height of the isotherm rtm = 8.0 g/kg.

          k = 2
          DO WHILE ( k <= gr%nnzp .AND. rtm(k) > 8.0e-3 )
            k = k + 1
          ENDDO
          IF ( k == gr%nnzp+1 .or. k == 2 ) THEN
            write(*,*) 'Identification of 8.0 g/kg level failed'
            write(*,*) 'Subroutine: dycoms2_rf02_tndcy. File: gcss.f'
            write(*,*) 'k = ', k
            write(*,*) 'rtm(k) = ', rtm(k)
            isValid = .false.
            return
          ENDIF
          z_i = ( (gr%zt(k)-gr%zt(k-1))/(rtm(k)-rtm(k-1)) )
     .        * (8.0e-3-rtm(k-1)) + gr%zt(k-1)

!         Compute the Heaviside step function for z - z_i.

          DO k = 1, gr%nnzp, 1
            IF ( gr%zm(k) - z_i  <  0.0 ) THEN
              Heaviside(k) = 0.0
            ELSEIF ( gr%zm(k) - z_i  ==  0.0 ) THEN
              Heaviside(k) = 0.5
            ELSEIF ( gr%zm(k) - z_i  >  0.0 ) THEN
              Heaviside(k) = 1.0
            ENDIF
          ENDDO

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
            ENDIF

          ENDDO 

        ! Compute the radiative heating rate.
        ! The radiative heating rate is defined on thermodynamic levels.

          DO k = 2, gr%nnzp, 1
            radht(k) = ( 1.0 / exner(k) ) * ( -1.0/(Cp*rhot(k)) )
     .               * ( Frad(k) - Frad(k-1) ) * gr%dzt(k)
          ENDDO
          radht(1) = radht(2)
        ENDIF ! ~ lbugsrad

        ! Enter the final theta-l and rtm tendencies

        IF ( .not. lbugsrad ) THEN
          theta_tndcy(1:gr%nnzp) = radht(1:gr%nnzp)
        ELSE
          theta_tndcy(1:gr%nnzp) = 0.0
        ENDIF 
        rt_tndcy(1:gr%nnzp) = 0.0

        ! Update surface statistics
#ifdef STATS
        if ( lstats_samp .and. izi > 0 ) then
          sfc%x(1,izi) = sfc%x(1,izi) + z_i
          sfc%n(1,izi) = sfc%n(1,izi) + 1
        endif
#endif
        if ( present( sclr_tndcy ) ) 
     .    sclr_tndcy(:,SCLR_THETA) = theta_tndcy 
        if ( present( sclr_tndcy ) ) sclr_tndcy(:,SCLR_RT) = rt_tndcy 


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
        DO k = 1, gr%nnzp, 1
           IF ( rcm(k) > 1.0E-6 ) THEN
              Ncm(k) = 55000000.0
           ELSE
              Ncm(k) = 0.0
           ENDIF
        ENDDO


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
#endif STATS

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

        real, optional, intent(in) :: sclrm_sfc(sclrm_dimension) 

!       Output variables

        real, intent(out) :: uw, vw, wtp, wrp

        real, optional, intent(out) :: wpsclrp_sfc(sclrm_dimension) 

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
          endif
        elseif ( sfctype == 1 ) then

          wtp = -Cd * ubar * ( thp - Tsfc/exnersfc )
          wrp = -Cd * ubar * ( rtp - rsat(psfc,Tsfc) )        

          if ( present( sclrm_sfc ).and.present( wpsclrp_sfc ) ) then
            wpsclrp_sfc(SCLR_THETA) = -Cd * ubar * 
     .                        ( sclrm_sfc(SCLR_THETA) - Tsfc/exnersfc )
            wpsclrp_sfc(SCLR_RT)    = -Cd * ubar * 
     .                        ( sclrm_sfc(SCLR_RT) - rsat(psfc, Tsfc) )
          endif
        else

          write(*,*) 'Invalid sfctype value = ', sfctype
          stop

        endif
#ifdef STATS
! Save ustar for statistics
!        x_sfc(1,iustar) = ustar
        if ( lstats_samp .and. iustar > 0 ) then
          sfc%x(1,iustar) = sfc%x(1,iustar) + ustar
          sfc%n(1,iustar) = sfc%n(1,iustar) + 1
        end if
#endif STATS

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
        REAL, OPTIONAL, INTENT(IN):: sclr_sfc(sclrm_dimension)

! Output
        REAL, INTENT(OUT):: uw
        REAL, INTENT(OUT):: vw
        REAL, INTENT(OUT):: wtp
        REAL, INTENT(OUT):: wrp
        REAL, INTENT(OUT):: ustar

! Output (optional)
        REAL, OPTIONAL, INTENT(OUT):: wpsclr_sfc(sclrm_dimension)

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

        ENDIF

        wtp = SH / (1.21 * Cp)
        wrp = LH / (1.21 * Lv)

        if (present( wpsclr_sfc )) then
          wpsclr_sfc(SCLR_THETA) = SH / (1.21 * Cp)

          wpsclr_sfc(SCLR_RT) = LH / (1.21 * Lv)
        endif

        RETURN
        END SUBROUTINE dycoms2_rf02_sfclyr

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for ASTEX KK case
!
        subroutine astex_tndcy( time, wmt, wmm, theta_tndcy, rt_tndcy,
     .                          sclr_tndcy )

        use constants, only: sclrm_dimension
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
     .        intent(out) :: sclr_tndcy(gr%nnzp, sclrm_dimension)

!       Internal variables

        integer i

!       Large-scale subsidence

        do i=2,gr%nnzp

           wmt(i) = - 5.e-6 * gr%zt(i)

        enddo

        wmm = zt2zm( wmt )

        ! Boundary conditions.
        wmt(1) = 0.0        ! Below surface
        wmm(1) = 0.0        ! At surface
        wmm(gr%nnzp) = 0.0  ! Model top
        
!       Radiative theta-l tendency

        theta_tndcy = 0.

!       Large scale advective moisture tendency

        rt_tndcy = 0.

        if (present( sclr_tndcy )) sclr_tndcy = 0.

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

        real, optional, intent(out) :: wpsclrp_sfc(sclrm_dimension)

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

        SUBROUTINE nov11_altocu_tndcy( time, time_initial, dt, thlm,
     .                                 rcm, p, exner, rhot, rtm, wmt,
     .                                 wmm, thl_tndcy, rt_tndcy,
     .                                 Frad, Frad_LW, Frad_SW,
     .                                 radht, radht_LW, radht_SW,
     .                                 sclr_tndcy )

        USE grid_class
        USE constants

        implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c INPUT VARIABLES                                                                            c
c time         :  time of simulation since start                               Unit: s       c
c time_initial :  initial time of simulation                                   Unit: s       c
c dt           :                                                               Unit: s       c
c thlm         :  profile of liquid water potential temperature                Unit: K       c
c rcm          :  profile of liquid (cloud) water mixing ratio                 Unit: kg/kg   c
c p            :  basic state pressure                                         Unit: hPa     c
c exner        :  basic state Exner function                                   Unit: NONE    c
c rhot         :  density of reference state at thermodynamic levels           Unit: kg/m^3  c
c                                                                                            c
c OUTPUT VARIABLES                                                                           c
c rtm          :  profile of total water mixing ratio                          Unit: kg/kg   c
c wmt          :  profile of vertical ascent/descent at thermodynamic levels   Unit: m/s     c
c wmm          :  profile of vertical ascent/descent at momentum levels        Unit: m/s     c
c                                                                                            c
c thl_tndcy    :  ??                                                           Unit: K ??    c
c rt_tndcy     :  ??                                                           Unit: kg/kg   c
c                                                                                            c
c Frad         :  profile of total radiative flux (lw + sw)                    Unit: W/m^2   c
c Frad_LW      :  profile of longwave radiative flux                           Unit: W/m^2   c
c Frad_SW      :  profile of shortwave radiative flux                          Unit: W/m^2   c
c                                                                                            c
c radht        :  profile of total radiative heating (lw + sw)                 Unit: K/s     c
c radht_LW     :  profile of longwave radiative heating                        Unit: K/s     c
c radht_SW     :  profile of shortwave radiative heating                       Unit: K/s     c
c                                                                                            c
c sclr_tndcy   :                                                               Unit: ??      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Frad_LW
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Frad_SW

        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: radht
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: radht_LW
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: radht_SW

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Output variables (optional) !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, OPTIONAL, DIMENSION(gr%nnzp, sclrm_dimension), 
     .                  INTENT(OUT) :: sclr_tndcy

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c LOCAL VARIABLES                                                                                     c
c dqc_dt_icedfs : change in liquid water over time, due to ice diffusional growth.     Unit: 1/s      c
c                                                                                                     c
c LWP           : liquid water path from domain top.                                   Unit: kg/m^2   c
c                                                                                                     c
c TEMPORARY ARRAYS USED FOR COAMPS RADIATIVE SCHEME                                                   c
c (SEE COMMENTS BELOW FOR SCHEME DESCRIPTION)                                                         c
c thm_rad       : temporary "flipped" array of potential temperature                   Unit: K        c
c rcm_rad       : temporary "flipped" array of liquid water mixing ratio               Unit: kg/kg    c
c rhot_rad      : temporary "flipped" array of air density                             Unit: kg/m^3   c
c dsigm         : temporary "flipped" array of grid spacing                            Unit: m        c
c                                                                                                     c
c frad_out      : temporary "flipped" array of total radiative flux from rad_lwsw      Unit: W/m^2    c
c frad_lw_out   : temporary "flipped" array of LW radiative flux from rad_lwsw         Unit: W/m^2    c
c frad_sw_out   : temporary "flipped" array of SW radiative flux from rad_lwsw         Unit: W/m^2    c
c                                                                                                     c
c radhtk        : temporary "flipped" array of total radiative heating from rad_lwsw   Unit: K/s      c
c radht_lw_out  : temporary "flipped" array of LW radiative heating from rad_lwsw      Unit: K/s      c
c radht_sw_out  : temporary "flipped" array of SW radiative heating from rad_lwsw      Unit: K/s      c
c                                                                                                     c
c INTERPOLATION ARRAYS AND CONSTANTS                                                                  c
c zsubs         : heights at which wmt data is supplied (used for subs interpolation)  Unit: m        c
c tsubs         : times after initialization at which wmt data is supplied             Unit: s        c
c                 (NOT USED IN NOV.11 CASE)                                                           c
c wtX(Y)        : vertical velocity specified at height Y and time X                   Unit: m/s      c
c               : (ONLY wt1 IS USED IN NOV.11 CASE)                                                   c
c w1-w2         : vertical velocity before (w1) and after (w2)                         Unit: m/s      c
c                 the current time at the specified level                                             c
c                 (NOT USED IN NOV.11 CASE)                                                           c
c                                                                                                     c
c ADDITIONAL PARAMETERS FOR NOV.11 SUBSIDENCE                                                         c
c subs_on       : logical variable tells us whether to turn subsidence on              Unit: NONE     c
c wmax          : defines value of maximum subsidence in profile                       Unit: cm/s     c
c zi            : defines approx. height of inversion within cloud                     Unit: m        c
c                 (subsidence is equal to wmax at this height)                                        c
c dazi          : defines height above inversion                                       Unit: m        c
c                 (above this height, subsidence linearly tapers off to zero)                         c
c dbzi          : defines height above inversion                                       Unit: m        c
c                 (below this height, subsidence linearly tapers off to zero)                         c
c dac           : defines height above cloud                                           Unit: m        c
c                 (at / above this height, we have NO subsidence)                      Unit: m        c
c dbc           : defines height below cloud                                           Unit: m        c
c                 (at / below this height, we have NO subsidence)                                     c
c                                                                                                     c
c RADIATION PARAMETERS                                                                                c
c sw_on         : logical variable passed to radiation scheme - is SW radiation on?    Unit: NONE     c
c lw_on         : logical variable passed to radiation scheme - is LW radiation on?    Unit: NONE     c
c center        : use centered differencing (as opposed to a one-sided forward         Unit: NONE     c
c                 difference) in radiation code                                                       c
c                                                                                                     c
c xi_abs        : cosine of the solar zenith angle                                     Unit: NONE     c
c F0            : coefficient for cloud top heating (see Stevens)                      Unit: W/m^2    c
c F1            : coefficient for cloud base heating (see Stevens)                     Unit: W/m^2    c
c kap           : "a constant" according to Duynkerke eqn. 5, where his value          Unit: m^2/kg   c
c                 is 130 m^2/kg.                                                                      c
c radius        : effective droplet radius                                             Unit: m        c
c AA            : albedo -- sea surface, according to Lenderink.                       Unit: NONE     c
c gc            : asymmetry parameter, "g" in Duynkerke.                               Unit: NONE     c
c Fs0           : The incident incoming SW insolation at cloud top in the              Unit: W/m^2    c
c                 direction of the incoming beam (not the vertical).                                  c
c omega         : single-scattering albedo                                             Unit: NONE     c
c                                                                                                     c
c SOLAR ZENITH ANGLE PARAMETERS (NOT USED IN NOV.11 CASE)                                             c
c c0            : coefficient for calculation of declination angle from                Unit: NONE     c
c                 Liou Table 2.2 and Eqn. 2.2.10                                                      c
c c1            : same as above                                                        Unit: NONE     c
c c2            : same as above                                                        Unit: NONE     c
c c3            : same as above                                                        Unit: NONE     c
c d1            : Same as above                                                                       c
c d2            : Same as above                                                                       c
c d3            : Same as above                                                                       c
c                                                                                                     c
c sda_t         : Linear function of day of the year.                                                 c
c                 sda_t=0 January 1 and sda_t -> 2*pi December 31.                                    c
c sda_delta     : Solar declination angle function from Liou 2.2.10                                   c
c sda_h         : Hour angle (Angle through which the earth must to turn                              c
c               to put sun directly overhead on a point's meridian)                                   c
c               (Angle between current time and solar noon)                                           c
c t_since_noon  : Number of seconds since noon (after noon > 0)                                       c
c julday        : Julian day of the year (January 1=1; December 31=365)                               c
c                                                                                                     c
c start_time_until_noon: number of seconds between start time and solar noon                          c
c                                                                                                     c
c Fs0 INTERPOLATION PARAMETERS                                                                        c
c nparam        : Number of Fs0 values in the list.                                                   c
c xilist        : Values of cosine of solar zenith angle corresponding to the values in Fslist        c
c Fslist        : Values of Fs0 corresponding to the values in xilist.                                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Local ice diffusional arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, DIMENSION(1:gr%nnzp):: dqc_dt_icedfs

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Local radiation arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variables used for calculating solar zenith angle                     !!
        !! (Currently not used in Nov.11 case)                                   !!
        !!                                                                       !!
        !! NOTE BY ADAM SMITH, 13 April 2006                                     !!
        !! In COAMPS, these constants are called c0, c1, c2, c3, d1, d2, and d3. !!
        !! Since HOC uses these variables elsewhere, I have renamed them here.   !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL c_0, c_1, c_2, c_3, d_1, d_2, d_3
        REAL sda_t, sda_delta, sda_h, t_since_noon, julday
        REAL start_time_until_noon

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Toggle for implementing differencing method in interpolations !!
        !! .FALSE. = forward differencing                                !!
        !! .TRUE.  = centered differencing                               !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LOGICAL :: center 

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Toggles for activating/deactivating forcings                            !!
        !! To turn off a specific forcing, set the corresponding toggle to .FALSE. !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
c      sda_t = 2*3.14*(julday-1)/365
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
c      sda_h = 2*3.14*t_since_noon/86400
c
c      xi_abs = sin(rlat*3.14/180) * sin(sda_delta) +
c     .         cos(rlat*3.14/180) * cos(sda_delta) * cos(sda_h)
c
c      xi_abs = max(xi_abs,0.)
c
c      print*, 'xi_abs=', xi_abs
c
       xi_abs = .4329

c Modification by Adam Smith 26 June 2006
c It is difficult to remember to set xi_abs = 0 when we want to shut off solar
c radiation.  If sw_on = .FALSE. above, we will automatically set xi_abs to 0
c to avoid confusion or errors.
      if ( .NOT. sw_on ) then
        xi_abs = 0.
      endif
c End of ajsmith4's Modification

      if (xi_abs .eq. 0.) then
        sw_on = .FALSE.
      else
        sw_on = .TRUE.
      endif


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

      call linearinterpolation(nparam,xilist,Fslist,xi_abs,Fs0)


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

c Nov. 11 subsidence parameters
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


        ! Special Nov. 11 condition for total water above cloud.
        ! One hour after the initial time, the total water above cloud
        ! is adjusted to be 0.89 of what it previously was.

        ! The conditional statement here is set so that if the timestep
        ! is such that there is no timestep at exactly 3600.0 seconds,
        ! then the operation still happnens at the first timestep and
        ! only the first timestep after 3600.0 seconds.
        IF ( time >= time_initial + 3600.0    .AND.
     .       time <  time_initial + 3600.0 + dt ) THEN

           DO k = 1, gr%nnzp, 1
              IF ( gr%zt(k) > 1555.0 ) THEN
                 rtm(k) = 0.89 * rtm(k)
              ENDIF
           ENDDO

        ENDIF


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Calling rad_lwsw to calculate LW and SW radiation
        ! More comments and description to follow
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

        ! We only implement this section if we choose not to use the
        ! BUGSRAD radiation scheme.
        IF ( .not. lbugsrad ) THEN

          ! This code transforms these profiles from HOC grid to COAMPS grid.
          ! The COAMPS-grid profiles are then passed to rad_lwsw for implementation.
          do k = 1, gr%nnzp
            thm_rad(k) = thlm(gr%nnzp-k+1) + Lv/Cp*rcm(gr%nnzp-k+1)
            rcm_rad(k) = rcm(gr%nnzp-k+1)
            rhot_rad(k) = rhot(gr%nnzp-k+1)
            dsigm(k) = 1.0 / gr%dzt(gr%nnzp-k+1)
          enddo

          ! Calling the radiation subroutine, which uses the COAMPS grid method.
          ! All input and output profiles use the COAMPS grid setup.
          call rad_lwsw(thm_rad, rcm_rad, rhot_rad, dsigm,
     .                  Frad_out, Frad_LW_out, Frad_SW_out,
     .                  radhtk, radht_LW_out, radht_SW_out,
     .                  gr%nnzp-1, center,
     .                  xi_abs, F0, F1, kap, radius, A, gc, Fs0, omega,
     .                  sw_on, lw_on)

          ! This code transforms the radiation results back into HOC grid setup.
          ! These Frad and radht arrays are actually applied to the HOC model.
          do k = 1, gr%nnzp-1
            Frad(k)     = Frad_out(gr%nnzp-k+1)
            Frad_LW(k)  = Frad_LW_out(gr%nnzp-k+1)
            Frad_SW(k)  = Frad_SW_out(gr%nnzp-k+1)

            radht(k)    = radhtk(gr%nnzp-k+1)
            radht_LW(k) = radht_LW_out(gr%nnzp-k+1)
            radht_SW(k) = radht_SW_out(gr%nnzp-k+1)
          enddo

          Frad(1) = Frad(2)
          Frad_LW(1) = Frad_LW(2)
          Frad_SW(1) = Frad_SW(2)

          radht(1) = radht(2)
          radht_LW(1) = radht_LW(2)
          radht_SW(1) = radht_SW(2)

        ENDIF ! ~ lbugsrad


        ! Compute the loss of total water due to diffusional growth of ice.
        ! This is defined on thermodynamic levels.

        IF (time >= time_initial + 3600.0) THEN

           CALL altocu_icedf( dt, thlm, rcm, p, exner, rhot,
     .                              dqc_dt_icedfs )

        ELSE

           DO k = 1, gr%nnzp, 1
              dqc_dt_icedfs(k) = 0.0
           ENDDO

        ENDIF


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Using linear interpolation scheme to interpolate subsidence          c
c                                                                      c
c FOR NOV.11 CASE                                                      c
c ---------------                                                      c
c As mentioned above, we want to implement a constant subsidence       c
c profile throughout the entire simulation (except for the 1 hour      c
c initial "spinup" period).  Because we do not have variations in      c
c subsidence over time, all tsubs sections (used in Jun.25 case) have  c
c been removed below.  Only one loop remains below, which will         c
c implement the interpolation subroutine, then return the w_ls value.  c
c Unlike jun25, we use this w_ls value directly instead of using       c
c interpolation to calculate a value between time steps.               c
c                                                                      c
c DIAGRAM OF NOV. 11 SUBSIDENCE PROFILE                                c
c -------------------------------------                                c
c       |      <- No subsidence in this region                         c
c       |                                                              c
c       ------------ Height = (zi + dazi + dac)                        c
c        \                                                             c
c         \    <- Subsidence tapers linearly in this region            c
c          \                                                           c
c           -------- Height = (zi + dazi)                              c
c           |                                                          c
c           |  <- Subsidence equals wmax in this region                c
c           |                                                          c
c           -------- Height = zi                                       c
c           |                                                          c
c           -------- Height = (zi - dbzi)                              c
c          /                                                           c
c         /    <- Subsidence tapers linearly in this region            c
c        /                                                             c
c       ------------ Height = (zi - dbzi - dbc)                        c
c       |                                                              c
c       |      <- No subsidence in this region                         c
c                                                                      c
c                                                                      c
c Comment by Adam Smith on 26 June 2006                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do k=2,gr%nnzp
          if ( (time >= time_initial + 3600.0) .and. subs_on ) then
            call linearinterpolation(7,zsubs,wt1,gr%zt(k),wmt(k))
          else
c           If time is not yet one hour, we have no subsidence
            wmt(k) = 0.0
          endif

          wmt(1) = wmt(2)
        enddo

        wmm = zt2zm(wmt)


        ! Enter the final theta-l and rtm tendencies

        DO k = 1, gr%nnzp, 1
           IF ( .not. lbugsrad ) THEN
             thl_tndcy(k) = radht(k)
     .                    - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           ELSE
             thl_tndcy(k) = - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           ENDIF
           rt_tndcy(k) = dqc_dt_icedfs(k)
        ENDDO

        IF ( present( sclr_tndcy ) ) THEN
          sclr_tndcy(:,SCLR_THETA) = thl_tndcy
          sclr_tndcy(:,SCLR_RT)    = rt_tndcy
        ENDIF

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
     .                                 dqc_dt_icedfs )

        USE grid_class
        USE constants

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
             ! Therefore, S = 1.0, allowing Si = es/ei.
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
             ENDIF

             !---------------Brian's comment--------------------------------------!
             ! dm = (dm/dt)*(dt/dz)*dz                                            !
             ! dm = (dm/dt)*(1/u_T)*dz                                            !
             !--------------------------------------------------------------------!
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

          ENDIF

        ENDDO

        RETURN
        END SUBROUTINE altocu_icedf

!-----------------------------------------------------------------------

        SUBROUTINE jun25_altocu_tndcy( time, time_initial, dt, thlm,
     .                                 rcm, p, exner, rhot, rtm, wmt,
     .                                 wmm, thl_tndcy, rt_tndcy, Frad,
     .                                 Frad_lw, Frad_sw,
     .                                 radht, radht_LW, radht_SW,
     .                                 sclr_tndcy )

        USE grid_class
        USE constants

        implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c INPUT VARIABLES                                                                            c
c time         :  time of simulation since start                               Unit: s       c
c time_initial :  initial time of simulation                                   Unit: s       c
c dt           :                                                               Unit: s       c
c thlm         :  profile of liquid water potential temperature                Unit: K       c
c rcm          :  profile of liquid (cloud) water mixing ratio                 Unit: kg/kg   c
c p            :  basic state pressure                                         Unit: hPa     c
c exner        :  basic state Exner function                                   Unit: NONE    c
c rhot         :  density of reference state at thermodynamic levels           Unit: kg/m^3  c
c                                                                                            c
c OUTPUT VARIABLES                                                                           c
c rtm          :  profile of total water mixing ratio                          Unit: kg/kg   c
c wmt          :  profile of vertical ascent/descent at thermodynamic levels   Unit: m/s     c
c wmm          :  profile of vertical ascent/descent at momentum levels        Unit: m/s     c
c                                                                                            c
c thl_tndcy    :  theta_l forcing                                              Unit: K/s     c
c rt_tndcy     :  total liquid water tendency                                  Unit: kg/kg/s c
c                                                                                            c
c Frad         :  profile of total radiative flux (lw + sw)                    Unit: W/m^2   c
c Frad_LW      :  profile of longwave radiative flux                           Unit: W/m^2   c
c Frad_SW      :  profile of shortwave radiative flux                          Unit: W/m^2   c
c                                                                                            c
c radht        :  profile of total radiative heating (lw + sw)                 Unit: K/s     c
c radht_LW     :  profile of longwave radiative heating                        Unit: K/s     c
c radht_SW     :  profile of shortwave radiative heating                       Unit: K/s     c
c                                                                                            c
c sclr_tndcy   :  passive scalar tendency                                      Unit: ??      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Frad_LW
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Frad_SW

        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: radht
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: radht_LW
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: radht_SW

        ! Output variables (optional)
        REAL, OPTIONAL, DIMENSION(gr%nnzp, sclrm_dimension),
     .                  INTENT(OUT) :: sclr_tndcy

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Local ice diffusional arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL, DIMENSION(1:gr%nnzp):: dqc_dt_icedfs

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Local radiation arrays !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        REAL, PARAMETER :: omega  = .992
        REAL, PARAMETER :: rlat = 37.6

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Additional SW radiative variables !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL xi_abs, Fs0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Variables used for calculating solar zenith angle                     !!
        !! (Currently not used in Nov.11 case)                                   !!
        !!                                                                       !!
        !! NOTE BY ADAM SMITH, 13 April 2006                                     !!
        !! In COAMPS, these constants are called c0, c1, c2, c3, d1, d2, and d3. !!
        !! Since HOC uses these variables elsewhere, I have renamed them here.   !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL c_0, c_1, c_2, c_3, d_1, d_2, d_3
        REAL sda_t, sda_delta, sda_h, t_since_noon, julday
        REAL start_time_until_noon

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Toggle for implementing differencing method in interpolations !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LOGICAL :: center

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Toggles for activating/deactivating forcings                            !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        julday = 176
        start_time_until_noon = 30540
        t_since_noon   = time - start_time_until_noon
        sda_t = 2*3.14*(julday-1)/365

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

        sda_h = 2*3.14*t_since_noon/86400

        xi_abs = sin(rlat*3.14/180) * sin(sda_delta) +
     .           cos(rlat*3.14/180) * cos(sda_delta) * cos(sda_h)


c        print*, 'xi_abs(before) = ', xi_abs
c        print*, 'sw_on(before) = ', sw_on

        xi_abs = max(xi_abs,0.)

c        print*, 'xi_abs(after) = ', xi_abs
c        print*, 'sw_on(after) = ', sw_on

c Modification by Adam Smith 15 June 2006
c It is difficult to remember to set xi_abs = 0 when we want to shut off solar
c radiation.  If sw_on = .FALSE. above, we will automatically set xi_abs to 0
c to avoid confusion or errors.

        if ( .NOT. sw_on ) then
          xi_abs = 0.
        endif

c        print*, 'xi_abs =', xi_abs

c End of Adam's Modification

        if (xi_abs .eq. 0.) then
          sw_on = .FALSE.
        else
          sw_on = .TRUE.
        endif

c        print*, 'xi_abs(end) = ', xi_abs
c        print*, 'sw_on(end) = ', sw_on
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

        call linearinterpolation(nparam,xilist,Fslist,xi_abs,Fs0)

ccccccccccccccccccccccccc
c                       c
c Subsidence Parameters c
c                       c
ccccccccccccccccccccccccc

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

        ! We only implement this section if we choose not to use the
        ! BUGSRAD radiation scheme.

        IF ( .not. lbugsrad ) THEN

          ! This code transforms these profiles from HOC grid to COAMPS grid.
          ! The COAMPS-grid profiles are then passed to rad_lwsw for implementation.
          do k = 1, gr%nnzp
            thm_rad(k) = thlm(gr%nnzp-k+1) + Lv/Cp*rcm(gr%nnzp-k+1)
            rcm_rad(k) = rcm(gr%nnzp-k+1)
            rhot_rad(k) = rhot(gr%nnzp-k+1)
            dsigm(k) = 1.0 / gr%dzt(gr%nnzp-k+1)
          enddo

          ! Calling the radiation subroutine, which uses the COAMPS grid method.
          ! All input and output profiles use the COAMPS grid setup.
          call rad_lwsw(thm_rad, rcm_rad, rhot_rad, dsigm,
     .                  Frad_out, Frad_LW_out, Frad_SW_out,
     .                  radhtk, radht_LW_out, radht_SW_out,
     .                  gr%nnzp-1, center,
     .                  xi_abs, F0, F1, kap, radius, AA, gc, Fs0, omega,
     .                  sw_on, lw_on)

c          do k = 1, gr%nnzp-1
c          print*,'time = ',time, ', k = ',k, ', Frad_SW_out = ',
c     .         Frad_SW_out(k)
c          enddo

          ! This code transforms the radiation results back into HOC grid setup.
          ! These Frad and radht arrays are actually applied to the HOC model.
          do k = 1, gr%nnzp-1
            Frad(k)     = Frad_out(gr%nnzp-k+1)
            Frad_LW(k)  = Frad_LW_out(gr%nnzp-k+1)
            Frad_SW(k)  = Frad_SW_out(gr%nnzp-k+1)

            radht(k)    = radhtk(gr%nnzp-k+1)
            radht_LW(k) = radht_LW_out(gr%nnzp-k+1)
            radht_SW(k) = radht_SW_out(gr%nnzp-k+1)
          enddo


          Frad(1) = Frad(2)
          Frad_LW(1) = Frad_LW(2)
          Frad_SW(1) = Frad_SW(2)

          radht(1) = radht(2)
          radht_LW(1) = radht_LW(2)
          radht_SW(1) = radht_SW(2)


c          ! "Flipping" input variables for COAMPS rad variable scheme
c          DO k = 1, gr%nnzp
c            rcm_rad(k+1)  = rcm(k)
c            rhot_rad(gr%nnzp-k+1) = rhot(k)
c            dsigm(gr%nnzp-k+1)   = gr%dzt(k)
c          ENDDO
c
c          call rad_lwsw(rcm_rad, rhot_rad, deltaz, gr%nnzp, center,
c     .                  xi_abs, F0, F1, kap, radius, AA, gc, Fs0, omega,
c     .                  sw_on, lw_on)
c
c          print*, 'hello3'
c
c          ! "Flipping" output variables to return to HOC variable scheme
c          DO k = 1, gr%nnzp
cc            frad(k)     = 0
cc            frad_lw(k)  = 0
cc            frad_sw(k)  = 0
cc
cc            radht(k)    = 0
cc            radht_lw(k) = 0
cc            radht_sw(k) = 0
c
c            frad(gr%nnzp-k+1)     = frad(k)
c            frad_lw(gr%nnzp-k+1)  = frad_lw(k)
c            frad_sw(gr%nnzp-k+1)  = frad_sw(k)
c
c            radht(gr%nnzp-k+1)    = rcm(k)
c            radht_lw(gr%nnzp-k+1) = rhot(k)
c            radht_sw(gr%nnzp-k+1) = gr%dzt(k)
c
c          ENDDO
c
c          print*, 'hello4'

        ENDIF ! ~ lbugsrad

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                               c
c Using linear interpolation scheme to interpolate subsidence   c
c Original code by Michael Falk                                 c
c Added for Jun.25 HOC simulation by Adam Smith, 13 April 2006  c
c                                                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c        do k=2,gr%nnzp
c          wmt(k) = 0
c        enddo
 
        if (time.lt.tsubs(1)) then
        do k=1,gr%nnzp
            call linearinterpolation(5,zsubs,wt1,gr%zt(k),wmt(k))
        enddo
 
        elseif (time.lt.tsubs(2)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt1,gr%zt(k),w1(k))
          call linearinterpolation(5,zsubs,wt2,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(1))/(tsubs(2)-tsubs(1))*(w2(k)-w1(k))+w1(k)
        enddo
 
        elseif (time.lt.tsubs(3)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt2,gr%zt(k),w1(k))
          call linearinterpolation(5,zsubs,wt3,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(2))/(tsubs(3)-tsubs(2))*(w2(k)-w1(k))+w1(k)
        enddo
 
        elseif (time.lt.tsubs(4)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt3,gr%zt(k),w1(k))
          call linearinterpolation(5,zsubs,wt4,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(3))/(tsubs(4)-tsubs(3))*(w2(k)-w1(k))+w1(k)
        enddo
 
        elseif (time.lt.tsubs(5)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt4,gr%zt(k),w1(k))
          call linearinterpolation(5,zsubs,wt5,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(4))/(tsubs(5)-tsubs(4))*(w2(k)-w1(k))+w1(k)
        enddo
 
        elseif (time.lt.tsubs(6)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt5,gr%zt(k),w1(k))
          call linearinterpolation(5,zsubs,wt6,gr%zt(k),w2(k))
        wmt(k) = (time-tsubs(5))/(tsubs(6)-tsubs(5))*(w2(k)-w1(k))+w1(k)
        enddo
 
        elseif (time.ge.tsubs(6)) then
        do k=2,gr%nnzp
          call linearinterpolation(5,zsubs,wt6,gr%zt(k),wmt(k))
        enddo
        endif

        wmt(1) = wmt(2)

        wmm = zt2zm(wmt)


        ! We don't use ice diffusion in the Jun.25 case, so we set
        ! dqc_dt_icedfs = 0 to avoid confusion and errors.
        DO k = 1, gr%nnzp, 1
            dqc_dt_icedfs(k) = 0.0
        ENDDO


        DO k = 1, gr%nnzp, 1
           IF ( .not. lbugsrad ) THEN
             thl_tndcy(k) = radht(k)
     .                    - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           ELSE
             thl_tndcy(k) = - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           ENDIF
           rt_tndcy(k) = dqc_dt_icedfs(k)
        ENDDO

        IF ( present( sclr_tndcy ) ) THEN
          sclr_tndcy(:,SCLR_THETA) = thl_tndcy
          sclr_tndcy(:,SCLR_RT)    = rt_tndcy
        ENDIF

        RETURN
        END SUBROUTINE jun25_altocu_tndcy


!-----------------------------------------------------------------------

        FUNCTION Diff_denom(Temp, press, ei)

        USE constants

        IMPLICIT NONE

        REAL, INTENT(IN):: Temp, press, ei
        REAL:: Diff_denom

        REAL:: Ka, Dv
        REAL:: Fk, Fd
        REAL:: Celsius

        REAL, PARAMETER:: Ls = 2.834e6

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
