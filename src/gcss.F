!----------------------------------------------------------------------
! $Id: gcss.F,v 1.8 2005-06-23 20:07:53 dschanen Exp $ 

! This file contains extra subroutines to simulate the GCSS ARM and
! GCSS BOMEX, GCSS ATEX, GCSS FIRE intercomparison cases as well 
! as Wangara day 33
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for ARM case
!----------------------------------------------------------------------
        subroutine arm_tndcy( time, theta_tndcy, rt_tndcy )

        use grid_class

        implicit none

!       Input
        real time

!       Output
        real theta_tndcy(gr%nnzp)
        real rt_tndcy(gr%nnzp)

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

        theta_tmp = ( 1. - a ) * ( atheta(i1) + rtheta(i1) )
     +                + a * ( atheta(i2) + rtheta(i2) )

        rt_tmp    = ( 1. - a ) * art(i1) + a * art(i2)

!       Convert to the right units

        theta_tmp = theta_tmp / 3600.
        rt_tmp    = rt_tmp / ( 3600. * 1000. )

!       Interpolate with respect to height

        do i=2, gr%nnzp

           if ( gr%zt(i) >= 0. .and. gr%zt(i) < 1000. ) then
              theta_tndcy(i) = theta_tmp
              rt_tndcy(i)    = rt_tmp
           elseif ( gr%zt(i) >= 1000. .and. gr%zt(i) < 3000. ) then
              b              = 1. - ( gr%zt(i) - 1000. ) / 2000.
              theta_tndcy(i) = b * theta_tmp
              rt_tndcy(i)    = b * rt_tmp
           else
              theta_tndcy(i) = 0.
              rt_tndcy(i)    = 0.
           endif

        enddo

        return
        end subroutine arm_tndcy

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for BOMEX case
!----------------------------------------------------------------------
        subroutine bomex_tndcy( time, wmt, wmm, theta_tndcy, rt_tndcy )

        use grid_class

        implicit none

!       Input
        real time

!       Output
        real wmt(gr%nnzp)
        real wmm(gr%nnzp)
        real theta_tndcy(gr%nnzp)
        real rt_tndcy(gr%nnzp)

!       Internal variables

        integer i

!       Large scale subsidence

        do i=2,gr%nnzp

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
        wmt(1) = 0.
        wmt(2) = 0.

        wmm = zt2zm( wmt )
        wmm(gr%nnzp) = 0.
        
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

        return
        end subroutine bomex_tndcy

!----------------------------------------------------------------------
! Subroutine to large-scale subsidence for FIRE case
!----------------------------------------------------------------------
        subroutine fire_tndcy( time, wmt, wmm, theta_tndcy, rt_tndcy )

        use grid_class

        implicit none

!       Input
        real time

!       Output
        real wmt(gr%nnzp)
        real wmm(gr%nnzp)
        real theta_tndcy(gr%nnzp)
        real rt_tndcy(gr%nnzp)

!       Internal variables

        integer i

!       Large-scale subsidence

        do i=2,gr%nnzp

           if ( gr%zt(i) >= 0. .and. gr%zt(i) < 1500. ) then
              wmt(i) = - 5.e-6 * gr%zt(i)
           endif

        enddo
        wmt(1) = 0.
        wmt(2) = 0.

        wmm = zt2zm( wmt )

        wmm(gr%nnzp) = 0.
        
!       Radiative theta-l tendency is computed interactively elsewhere

        theta_tndcy = 0.

!       Large scale advective moisture tendency

        rt_tndcy    = 0.

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
        real, optional, dimension(sclrm_dimension, gr%nnzp), 
     .        intent(out) :: sclr_tndcy

!       No large-scale subsidence for now
        wmt = 0.
        wmm = 0.

!       No large-scale water tendency or cooling

        rt_tndcy    = 0.
        theta_tndcy = 0.

!       Scalar r forcing is just set to be rtm_forcing for now
        if ( present(sclr_tndcy) ) sclr_tndcy = 0.

        return
        end subroutine wangara_tndcy

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for ATEX case
!
        subroutine atex_tndcy( time, time_initial, rtm, wmt, wmm,
     .                         theta_tndcy, rt_tndcy, isValid )

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
           wmt(1) = 0.
           wmt(2) = 0.

           wmm = zt2zm( wmt )
           wmm(gr%nnzp) = 0.
        
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

        end if

        return
        end subroutine atex_tndcy

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to GCSS ARM specifications
!
        subroutine arm_sfclyr( time, z, dn0, thp, up, vp, uw, vw, 
     .                         wtp, wrp, ustar )

        use constants, only: grav
        implicit none

!       External

        real diag_ustar
        external diag_ustar

!       Constant

        real ubmin
        parameter (ubmin = 0.25)
        real Cp, Lv
        parameter ( Cp = 1004.67 )
        parameter ( Lv = 2.5e6 )
!        real g
!        parameter ( g = 9.81 )

!       Input variables

        real time
        real z
        real dn0
        real thp
        real up,vp
        real ustar

!       Output variables

        real uw,vw,wtp,wrp

!       Internal variables

        integer i,j
        real usfc,vsfc,ubar
        real true_time
        real heat_flx,moisture_flx
        real heat_flx2,moisture_flx2
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

        return

        end subroutine arm_sfclyr

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to GCSS BOMEX specifications
!
        subroutine bomex_sfclyr( up, vp, uw, vw, wtp, wrp )

        implicit none

!       Constant

        real, parameter :: ubmin = 0.25
        real, parameter :: ustar = 0.28

!       Input variables

        real, intent(in) :: up,vp

!       Output variables

        real, intent(out) :: uw,vw,wtp,wrp

!       Internal variables

        real ubar

! Compute heat and moisture fluxes

        wtp = 8.e-3
        wrp = 5.2e-5

! Compute momentum fluxes

        ubar = max(ubmin,sqrt(up*up + vp*vp))

        uw = -up * ustar*ustar / ubar
        vw = -vp * ustar*ustar / ubar

        return
        end subroutine bomex_sfclyr

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to GCSS ATEX specifications
!
        subroutine atex_sfclyr( up, vp, thp, rtp, uw, vw, wtp, wrp )

        use constants

        implicit none

!       Constants

        real, parameter :: ubmin = 0.25
        real, parameter :: ustar = 0.3
        real, parameter :: C_10 = 0.0013
        real, parameter :: SST = 298

!       Input variables

        real, intent(in) :: up,vp,thp,rtp

!       Output variables

        real, intent(out) :: uw,vw,wtp,wrp

!       Internal variables

        real ubar

! Compute heat and moisture fluxes

        ubar = max(ubmin,sqrt(up*up + vp*vp))

        wtp  = -C_10 * ubar * ( thp - SST * (1000./1015.)**kappa )
        wrp  = -C_10 * ubar * ( rtp - 0.0198293 )

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

! Compute new mixing variable as wrp
        if ( present(wpsclrp_sfc) ) wpsclrp_sfc(1) = 1.3e-4 * wtp

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
!
        subroutine cloud_rad( rhot, rcm, Frad, htrt )

        use grid_class
        use constants

        implicit none

!       Input

        real, intent(in) :: rhot(gr%nnzp)       ! density (thermo point)
        real, intent(in) :: rcm(gr%nnzp)        ! liquid water mixing ratio

!        Output

        real, intent(out) :: Frad(gr%nnzp)        ! IR radiative flux
        real, intent(out) :: htrt(gr%nnzp)        ! Radiative heating rate

        real F0,k
        parameter (F0 = 74., k = 130.)

        integer i,im1,ip1

!       Compute liquid water path from top of the model
!       We define liquid water path on momentum levels

        Frad(gr%nnzp) = 0.
        do i=gr%nnzp-1,2,-1
           Frad(i) = Frad(i+1) + rhot(i+1) * rcm(i+1) / gr%dzt(i+1)
        enddo
        Frad(1) = Frad(2)

!       Compute IR radiative flux

        do i=1,gr%nnzp
           Frad(i) = F0 * exp( -k * 1.0 * Frad(i) )
        enddo

!       Compute IR heating rate

        htrt          = htrt - 1. / ( Cp * rhot ) * ddzm( Frad )
        htrt(1)       = 0.
        htrt(gr%nnzp) = 0.

        return
        end subroutine cloud_rad

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for DYCOMS RF01 case
!
        subroutine dycoms2_rf01_tndcy( time, rhot, rhom, rtm, rcm,
     .                                 wmt, wmm, Frad, theta_tndcy, 
     .                                 rt_tndcy, isValid )

        use grid_class
        use constants
#ifdef STATS
        use statistics
#endif STATS

        implicit none

!       Input
        real, intent(in) :: time
        real, intent(in) :: rhot(gr%nnzp)
        real, intent(in) :: rhom(gr%nnzp)
        real, intent(in) :: rtm(gr%nnzp)
        real, intent(in) :: rcm(gr%nnzp)

!       Output
        real, intent(out) :: wmt(gr%nnzp)
        real, intent(out) :: wmm(gr%nnzp)
        real, intent(out) :: Frad(gr%nnzp)
        real, intent(out) :: theta_tndcy(gr%nnzp)
        real, intent(out) :: rt_tndcy(gr%nnzp)

!       Misc
        logical, intent(inout) :: isValid

!       Parameter

        real, parameter :: lsdiv =  3.75e-6
        real, parameter :: F0 = 70.0, F1 = 22.0, kay = 85.0

!       Internal variables

        integer i
        real zi
        real lwp(gr%nnzp), htrt(gr%nnzp)

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
           write(*,*) 'Subroutine: dycoms_tndcy. File: gcss.f'
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
        wmt(1) = 0.
        wmt(2) = 0.

        wmm = zt2zm( wmt )
        wmm(gr%nnzp) = 0.
        
!       theta-l radiative tendency

!       Compute liquid water path from top of the model
!       We define liquid water path on momentum levels

        lwp(gr%nnzp) = 0.0
        do i=gr%nnzp-1,2,-1
           lwp(i) = lwp(i+1) + rhot(i+1) * rcm(i+1) / gr%dzt(i+1)
        enddo
        lwp(1) = lwp(2)
!        x_sfc(1,ilwp) = lwp(2)

!       Compute IR radiative flux

        do i=1,gr%nnzp
           Frad(i) = F0 * exp( -kay *lwp(i) )
     .             + F1 * exp( -kay * (lwp(1)-lwp(i)) )
           if ( zi > 0 .and. gr%zm(i) > zi ) then
             Frad(i) = Frad(i)
     .               + rhom(i) * cp * lsdiv
     .                 * ( 0.25*(gr%zm(i)-zi)**(4.0/3.0)
     .                     + zi*(gr%zm(i)-zi)**(1.0/3.0) )
           endif
        enddo

!       Compute IR heating rate

        htrt          = - 1. / ( Cp * rhot ) * ddzm( Frad )
        htrt(1)       = 0.
        htrt(gr%nnzp) = 0.

!       Add heating rate to theta-l forcing

        theta_tndcy = theta_tndcy + htrt
        
        return
        end subroutine dycoms2_rf01_tndcy

!----------------------------------------------------------------------

        SUBROUTINE dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, wmt, wmm,
     .                                 theta_tndcy, rt_tndcy, isValid )

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

        ! Output
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: wmt
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: wmm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: theta_tndcy
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: rt_tndcy

        ! Misc
        LOGICAL, INTENT(INOUT) :: isValid

        ! Local
        REAL, DIMENSION(1:gr%nnzp):: LWP
        REAL, DIMENSION(1:gr%nnzp):: Heaviside
        REAL, DIMENSION(1:gr%nnzp):: F_rad
        REAL, DIMENSION(1:gr%nnzp):: htrt

        REAL:: z_i

        REAL, PARAMETER :: ls_div = 3.75e-6
        REAL, PARAMETER :: kap    = 85.0 ! m^2/kg
        REAL, PARAMETER :: F0     = 70.0  ! W/m^2
        REAL, PARAMETER :: F1     = 22.0  ! W/m^2

        INTEGER:: k

        ! Large-scale subsidence

        DO k = 2, gr%nnzp, 1
           wmt(k) = -ls_div * gr%zt(k)
        ENDDO
        wmt(1) = 0.0
        wmt(2) = 0.0

        wmm = zt2zm( wmt )
        wmm(gr%nnzp) = 0.0

        ! Radiation

!       Compute liquid water path from top of the model
!       We define liquid water path on momentum levels

        LWP(gr%nnzp) = 0.0
        DO k = gr%nnzp-1, 2, -1
           LWP(k) = LWP(k+1) + rhot(k+1)*rcm(k+1)/gr%dzt(k+1)
        ENDDO
        LWP(1) = LWP(2)

!       Find the height of the isotherm rtm = 8.0 g/kg.

        k = 2
        DO WHILE ( k <= gr%nnzp .AND. rtm(k) > 8.0e-3 )
           k = k + 1
        ENDDO
        IF ( k == gr%nnzp+1 .or. k == 2 ) THEN
           write(*,*) 'Identification of 8.0 g/kg level failed'
           write(*,*) 'Subroutine: dycoms2_tndcy. File: gcss.f'
           write(*,*) 'k = ', k
           write(*,*) 'rtm(k) = ', rtm(k)
           isValid = .false.
           return
        ENDIF
        z_i = (gr%zt(k)-gr%zt(k-1))/(rtm(k)-rtm(k-1))*(8.0e-3-rtm(k-1))
     .       + gr%zt(k-1)

!       Compute the Heaviside step function for z - z_i.

        DO k = 1, gr%nnzp, 1
           IF ( gr%zm(k) - z_i  <  0.0 ) THEN
              Heaviside(k) = 0.0
           ELSEIF ( gr%zm(k) - z_i  ==  0.0 ) THEN
              Heaviside(k) = 0.5
           ELSEIF ( gr%zm(k) - z_i  >  0.0 ) THEN
              Heaviside(k) = 1.0
           ENDIF
        ENDDO

!       Compute radiative flux profile (F_rad).
!       Radiative flux is defined on momentum levels.

        DO k = 1, gr%nnzp, 1

           F_rad(k) =  F0 * EXP( -kap * LWP(k) )
     .               + F1 * EXP( -kap * (LWP(1) - LWP(k)) )

           IF ( Heaviside(k) > 0.0 ) THEN
              F_rad(k) = F_rad(k)
     .                  + rhom(k) * Cp * ls_div * Heaviside(k)
     .                  * ( 0.25 * ((gr%zm(k)-z_i)**(4.0/3.0))
     .                     + z_i * ((gr%zm(k)-z_i)**(1.0/3.0)) )
          ENDIF

        ENDDO 

        ! Compute the radiative heating rate.
        ! The radiative heating rate is defined on thermodynamic levels.

        DO k = 2, gr%nnzp, 1
           htrt(k) = ( -1.0/(Cp*rhot(k)) ) 
     .              * ( F_rad(k) - F_rad(k-1) ) * gr%dzt(k)
        ENDDO
        htrt(1) = htrt(2)

        ! Enter the final theta-l and rtm tendencies

        DO k = 1, gr%nnzp, 1
           theta_tndcy(k) = htrt(k)
           rt_tndcy(k) = 0.0
        ENDDO 

        ! Update surface statistics
#ifdef STATS
        if ( lstats_samp .and. izi > 0 ) then
          sfc%x(1,izi) = sfc%x(1,izi) + z_i
          sfc%n(1,izi) = sfc%n(1,izi) + 1
        endif
#endif

        END SUBROUTINE dycoms2_rf02_tndcy

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to GCSS DYCOMS specifications
!
        subroutine dycoms2_rf01_sfclyr( sfctype, Tsfc, psfc, exnersfc,
     .                                  up, vp, thp, rtp, rhop,
     .                                  uw, vw, wtp, wrp )

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

!       Output variables

        real, intent(out) :: uw,vw,wtp,wrp

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

        elseif ( sfctype == 1 ) then

          wtp = -Cd * ubar * ( thp - Tsfc/exnersfc )
          wrp = -Cd * ubar * ( rtp - rsat(psfc,Tsfc) )        

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

        SUBROUTINE dycoms2_rf02_sfclyr( u_sfc, v_sfc, thl_sfc, rt_sfc,
     .                                  p_sfc, uw, vw, wtp, wrp )

        USE constants

        REAL, INTENT(IN):: u_sfc
        REAL, INTENT(IN):: v_sfc
        REAL, INTENT(IN):: thl_sfc
        REAL, INTENT(IN):: rt_sfc
        REAL, INTENT(IN):: p_sfc

        REAL, INTENT(OUT):: uw
        REAL, INTENT(OUT):: vw
        REAL, INTENT(OUT):: wtp
        REAL, INTENT(OUT):: wrp

        REAL:: wind_sfc, Cd

        REAL, PARAMETER:: ustar = 0.28

        wind_sfc = SQRT( u_sfc**2.0 + v_sfc**2.0 )

        IF (wind_sfc > 0.0) THEN

           ! Calculation of drag coefficient.
           Cd = (ustar/wind_sfc)**2.0

        ELSE

           ! If there is no wind, there is no drag coefficient.
           Cd = 0.0

        ENDIF

        uw = -Cd * wind_sfc * u_sfc
        vw = -Cd * wind_sfc * v_sfc

        wtp = -Cd * wind_sfc * (thl_sfc - 292.0*(p0/p_sfc)**kappa)
        wrp = -Cd * wind_sfc * (rt_sfc - 0.013581)

        END SUBROUTINE dycoms2_rf02_sfclyr

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for ASTEX KK case
!
        subroutine astex_tndcy( time, wmt, wmm, theta_tndcy, rt_tndcy )

        use grid_class

        implicit none

!       Input
        real time

!       Output
        real wmt(gr%nnzp)
        real wmm(gr%nnzp)
        real theta_tndcy(gr%nnzp)
        real rt_tndcy(gr%nnzp)

!       Internal variables

        integer i

!       Large-scale subsidence

        do i=2,gr%nnzp

           wmt(i) = - 5.e-6 * gr%zt(i)

        enddo

        wmt(1) = 0.
        wmt(2) = 0.

        wmm = zt2zm( wmt )
        wmm(gr%nnzp) = 0.
        
!       Radiative theta-l tendency

        theta_tndcy = 0.

!       Large scale advective moisture tendency

        rt_tndcy = 0.

        return
        end subroutine astex_tndcy

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture according to ASTEX with Khairoutdinov and Kogan
! alteration.

        subroutine astex_sfclyr( rho0, up, vp, uw, vw, wtp, wrp )

        use constants

        implicit none

!       Input variables

        real rho0,up,vp

!       Output variables

        real uw,vw,wtp,wrp

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

        return
        end subroutine astex_sfclyr

!-----------------------------------------------------------------------

        SUBROUTINE nov11_altocu_tndcy( time, time_initial, dt, thlm,
     .                                 rcm, p, exner, rhot, wmt, wmm,
     .                                 thl_tndcy, rt_tndcy, Frad,
     .                                 Frad_lw, Frad_sw )

        USE grid_class
        USE constants

        implicit none

        ! Input variables
        REAL, INTENT(IN):: time
        REAL, INTENT(IN):: time_initial
        REAL, INTENT(IN):: dt
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: thlm
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rcm
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: p
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: exner
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rhot

        ! Output variables
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: wmt
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: wmm
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: thl_tndcy
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: rt_tndcy
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Frad
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Frad_lw
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Frad_sw

        ! Local variables
        REAL, DIMENSION(1:gr%nnzp):: LWP
        REAL, DIMENSION(1:gr%nnzp):: htrt
        REAL, DIMENSION(1:gr%nnzp):: dqc_dt_icedfs

        REAL, PARAMETER:: wmax = -3.0E-2
        REAL, PARAMETER:: kap = 100.0 ! m^2/kg
        REAL, PARAMETER:: F0 = 107.0  ! W/m^2
        REAL, PARAMETER:: F1 = 61.0   ! W/m^2

        INTEGER:: k


        ! Apply subsidence based on height after 3600.0 seconds into
        ! the model run.

        DO k = 1, gr%nnzp, 1

           IF (time >= time_initial + 3600.0) THEN
              IF (gr%zt(k) >= 0.0 .AND. gr%zt(k) <= 100.0) THEN
                 wmt(k) = 0.0
              ELSEIF (gr%zt(k) > 100.0 .AND. gr%zt(k) <= 400.0) THEN
                 wmt(k) = ( (wmax - 0.0)/(400.0 - 100.0) )
     .                   * (gr%zt(k) - 100.0) + 0.0
              ELSEIF (gr%zt(k) > 400.0 .AND. gr%zt(k) <= 1400.0) THEN
                 wmt(k) = wmax
              ELSEIF (gr%zt(k) > 1400.0 .AND. gr%zt(k) <= 1800.0) THEN
                 wmt(k) = wmax
              ELSEIF (gr%zt(k) > 1800.0 .AND. gr%zt(k) <= 2000.0) THEN
                 wmt(k) = ( (0.0 - wmax)/(2000.0 - 1800.0) )
     .                   * (gr%zt(k) - 1800.0) + wmax
              ELSEIF (gr%zt(k) > 2000.0) THEN
                 wmt(k) = 0.0
              ENDIF
           ELSE
              wmt(k) = 0.0
           ENDIF

        ENDDO

        wmm = zt2zm(wmt)


        ! Long-wave Radiation

!       Compute liquid water path from top of the model
!       We define liquid water path on momentum levels

        LWP(gr%nnzp) = 0.0
        DO k = gr%nnzp-1, 1, -1
           LWP(k) = LWP(k+1) + rhot(k+1)*rcm(k+1)/gr%dzt(k+1)
        ENDDO

!       Compute long-wave radiative flux profile (Frad_lw).
!       Radiative flux is defined on momentum levels.

        DO k = 1, gr%nnzp, 1

           Frad_lw(k) =  F0 * EXP( -kap * LWP(k) )
     .                 + F1 * EXP( -kap * (LWP(1) - LWP(k)) )

        ENDDO


        ! Short-wave Radiation
        ! Radiative flux is defined on momentum levels.

        CALL nov11_altocu_sw(rhot, rcm, Frad_sw)


        ! Total Radiation
        ! Radiative flux is defined on momentum levels.

        DO k = 1, gr%nnzp, 1
           Frad(k) = Frad_lw(k) + Frad_sw(k)
        ENDDO


        ! Compute the radiative heating rate.
        ! The radiative heating rate is defined on thermodynamic levels.

        DO k = 2, gr%nnzp, 1
           htrt(k) = ( -1.0/(Cp*rhot(k)) )
     .              * ( Frad(k) - Frad(k-1) ) * gr%dzt(k)
        ENDDO
        htrt(1) = htrt(2)


        ! Compute the loss of total water due to diffusional growth of ice.
        ! This is defined on thermodynamic levels.

        IF (time >= time_initial + 3600.0) THEN

           CALL nov11_altocu_icedf( dt, thlm, rcm, p, exner, rhot,
     .                              dqc_dt_icedfs )

        ELSE

           DO k = 1, gr%nnzp, 1
              dqc_dt_icedfs(k) = 0.0
           ENDDO

        ENDIF


        ! Enter the final theta-l and rtm tendencies

        DO k = 1, gr%nnzp, 1
           thl_tndcy(k) = htrt(k)
     .                   - ( Lv/(Cp*exner(k)) ) * dqc_dt_icedfs(k)
           rt_tndcy(k) = dqc_dt_icedfs(k)
        ENDDO


        END SUBROUTINE nov11_altocu_tndcy

!-----------------------------------------------------------------------

        ! This subroutine is based on a COAMPS subroutine (sunray_sw)
        ! written by Michael Falk to calculate short-wave radiation
        ! for the Nov. 11 altocumulus case.

        SUBROUTINE nov11_altocu_sw( rhot, rcm, Frad_sw )

        USE grid_class

        implicit none

        ! Input variables
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rhot
        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: rcm

        ! Output variables
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: Frad_sw

        ! Local variables
        REAL, DIMENSION(1:gr%nnzp):: tau
        REAL, DIMENSION(1:gr%nnzp):: taude
        REAL, DIMENSION(1:gr%nnzp):: F_diff
        REAL, DIMENSION(1:gr%nnzp):: F_dir

        REAL:: taupath, tauc, t1, t2, t3, cc1, cc2, omegade,
     .         x1, x2, x3, rk, rk2, xi_abs2, rp, alpha_2, beta_2, rtt,
     .         exmu0, expk, exmk, xp23p, xm23p, ap23b, am23b, taucde

        REAL:: ff, gcde

        REAL, PARAMETER:: radius = 1.0e-5
        REAL, PARAMETER:: AA     = 0.1
        REAL, PARAMETER:: gc     = 0.86
        REAL, PARAMETER:: omega  = .9965
        REAL, PARAMETER:: xi_abs = .4329
        REAL, PARAMETER:: Fs0    = 1212.75
        REAL, PARAMETER:: rho_l  = 1000.0

        INTEGER:: k


        ff = gc*gc
        gcde = gc/(1.0+gc)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                     !
        !  Computation of tau and omega variables                             !
        !                                                                     !
        !                                                                     !
        ! tauc   : column total optical depth                   Unit: none    !
        ! taupath: column total Delta-Eddington optical depth.  Unit: none    !
        ! omega  : single-scattering albedo                     Unit: none    !
        ! omegade: D-E omega-- from Duynkerke eqn.18.           Unit: none    !
        ! taucde : D-E tauc -- from Duynkerke eqn.18.           Unit: none    !
        ! taude  : D-E tau  -- from Duynkerke eqn.18.           Unit: none    !
        !                                                                     !
        ! Comments by Michael Falk, 26 January 2005                           !
        !                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! tau is defined on momentum levels.  Brian

        tauc         = 0.0
        tau(gr%nnzp) = 0.0      ! optical depth at top of atm. is 0.  Brian.
        do k = gr%nnzp-1, 1, -1
          tau(k) = 1.5 * (rcm(k+1)*rhot(k+1)/gr%dzt(k+1))
     .           / radius / rho_l
          tauc   = tauc + tau(k)
        enddo

        omegade = (1.0-ff)*omega/(1.0-omega*ff)
        taucde  = (1.0-omega*ff)*tauc

        do k = 1, gr%nnzp, 1
          taude(k) = (1.0-omega*ff)*tau(k)
        enddo


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                     !
        !  Computation of constants for radiative transfer equations          !
        !                                                                     !
        !  These variables come from Duynkerke eqn.20, which, with slight     !
        !  modifications, matches Shettle and Weinman between eqns.12b and 13.!
        !  Duynkerke uses slightly different formulations than Shettle and    !
        !  Weinman:                                                           !
        !                                                                     !
        ! F0(S&W)    = F0(Duynkerke)*pi.                                      !
        ! alpha(S&W) = alpha(Duynkerke)*F0(S&W).                              !
        ! beta(S&W)  = beta(Duynkerke)*F0(S&W).                               !
        ! c1(S&W)    = c1(Duynkerke)*F0(S&W).                                 !
        ! c2(S&W)    = c2(Duynkerke)*F0(S&W).                                 !
        ! x3(S&W)    = x3(Duynkerke)*F0(S&W).                                 !
        ! y3(S&W)    = y3(Duynkerke)*F0(S&W).                                 !
        !                                                                     !
        ! F0 is divided out of each term in several equations, and then       !
        ! reintroduced when the actual radiative flux is computed.  The       !
        ! computations here follow Lenderink's original sunray_sw code which  !
        ! uses the Duynkerke formulations for these variables.                !
        !                                                                     !
        ! x1     : term 1 in k equation                         Unit: none    !
        ! x2     : term 2 in k equation                         Unit: none    !
        ! rk     : k equation                                   Unit: none    !
        ! rk2    : k equation squared                           Unit: none    !
        ! x3     : term in denominator of alpha and beta        Unit: none    !
        ! rp     : p equation                                   Unit: none    !
        ! alpha  : alpha equation                               Unit: none    !
        ! beta   : beta equation                                Unit: none    !
        !                                                                     !
        ! The following variables are used by Lenderink to solve for          !
        ! Duynkerke's parameters C1 and C2.  They are all dimensionless.      !
        !                                                                     !
        ! rtt    : 2/3.                                                       !
        ! exmu0  : exponential term used in S&W eqn.14-- originally from      !
        !          eqn 1 (also Salby 9.35) in the source function for SW      !
        !          radiation                                                  !
        ! expk   : one of the coefficients of C2 on the left side of Shettle  !
        !          and Weinman eqn.14, originally from eqn.12a and eqn.12b    !
        ! exmk   : reciprocal of expk, one of the coefficients of C1 on the   !
        !          left side of Shettle and Weinman eqn.14, originally from   !
        !          eqn.12a and eqn.12b                                        !
        ! xp23p  : coefficient of C1 - left side of Shettle and Weinman eqn.13!
        ! xm23p  : coefficient of C2 - left side of Shettle and Weinman eqn.13!
        ! ap23b  : right side of Shettle and Weinman eqn.13                   !
        ! t1     : the other coefficient of C1 on left side of Shettle and    !
        !          Weinman eqn.14                                             !
        ! t2     : the other coefficient of C2 on left side of Shettle and    !
        !          Weinman eqn.14                                             !
        ! t3     : the coefficient of exmu0 on the right side of Shettle and  !
        !          Weinman eqn.14                                             !
        !                                                                     !
        ! Comments by Michael Falk, 26 January 2005                           !
        !                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        x1 = 1.0-omegade*gcde
        x2 = 1.0-omegade
        rk = sqrt(3.0*x2*x1)
        xi_abs2 = xi_abs*xi_abs
        rk2 = rk*rk
        x3  = 4.0*(1.0-rk2*xi_abs2)
        rp  = sqrt(3.0*x2/x1)
        alpha_2 = 3.0*omegade*xi_abs2*(1.0+gcde*x2)/x3
        beta_2  = 3.0*omegade*xi_abs*(1.0+3.0*gcde*xi_abs2*x2)/x3

        rtt   = 2.0/3.0
        exmu0 = exp(-taucde/xi_abs)
        expk  = exp(rk*taucde)
        exmk  = 1.0/expk
        xp23p = 1.0+rtt*rp
        xm23p = 1.0-rtt*rp
        ap23b = alpha_2+rtt*beta_2

        t1 = 1.-AA-rtt*(1.0+AA)*rp
        t2 = 1.-AA+rtt*(1.0+AA)*rp
        t3 = (1.-AA)*alpha_2-rtt*(1.+AA)*beta_2+AA*xi_abs


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                     !
        ! Shettle and Weinman 13 and 14, adapted into Duynkerke, give two     !
        ! equations and two unknowns which can be solved by linear combination!
        ! (CC2) and then substitution (CC1).                                    !
        !                                                                     !
        ! Comments by Michael Falk, 26 January 2005                           !
        !                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        cc2 = (xp23p*t3*exmu0-t1*ap23b*exmk)
     .     / (xp23p*t2*expk-xm23p*t1*exmk)
        cc1 = (ap23b-cc2*xm23p)/xp23p


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                     !
        !  Computation of diffuse and direct shortwave flux                   !
        !                                                                     !
        ! F_diff is the first term in Duynkerke eqn.19.  The F0 and pi        !
        ! constants which are divided out in the Duynkerke formulation are    !
        ! reintroduced here.                                                  !
        !                                                                     !
        ! Duynkerke eqn.19's F_diff term comes from Shettle and Weinman eqn.8,!
        ! where F_diff = F(upward)-F(downward).  Then F_diff = (-4/3)*pi*I1,  !
        ! where I1 is solved in Shettle and Weinman eqn.12b.  Capital P in    !
        ! Shettle and Weinman eqn.12b should actually be a lowercase p.       !
        !                                                                     !
        ! F_dir is the second term in Duynkerke eqn.19.                       !
        !                                                                     !
        ! The negative sign for F_diff and F_dir is related to the definition !
        ! of which direction is a positive flux.                              !
        !                                                                     !
        ! Comments by Michael Falk, 26 January 2005                           !
        !                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Since tau is on a momentum level, taude has to be on a momentum
        ! level.  taupath is a running total of taude from the top of the model
        ! on down.  It is also located on a momentum level.  Brian

        taupath = 0.

        DO k = gr%nnzp, 1, -1

           taupath = taupath + taude(k)

           F_diff(k) = (-4.0/3.0) * Fs0
     .               * (
     .                  rp *
     .                      (
     .                         cc1*exp(-rk*taupath)
     .                       - cc2*exp(rk*taupath)
     .                      )
     .                  - beta_2*exp(-taupath/xi_abs)
     .                 )

           F_dir(k) = -Fs0*xi_abs*exp(-taupath/xi_abs)

           Frad_sw(k) = F_diff(k) + F_dir(k)

        ENDDO


        END SUBROUTINE nov11_altocu_sw

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

        SUBROUTINE nov11_altocu_icedf( dt, thlm, rcm, p, exner, rhot,
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
        REAL:: Diff_denom

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

          ! Check whether we're in cloud and below freezing
          IF ( rcm(k) >= 0.0 .AND. TT(k) < 273.15 ) THEN

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
     .          * gr%dzt(k)
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


        END SUBROUTINE nov11_altocu_icedf

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

        END FUNCTION Diff_denom

!------------------------------------------------------------------------
