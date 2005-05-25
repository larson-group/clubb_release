! This file contains extra subroutine to simulate the GCSS ARM and
! GCSS BOMEX, GCSS ATEX, GCSS FIRE intercomparison cases as well 
! as Wangara day 33

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for ARM case
!
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
        real a,b
        real true_time
        real theta_tmp,rt_tmp

        real atheta(6), rtheta(6), art(6)
        data atheta / 0.000, 0.000, 0.000, -0.080, -0.160, -0.160 /
        data rtheta / -0.125, 0.000, 0.000, 0.000, 0.000, -0.100 /
        data art    / 0.080, 0.020, -0.040, -0.100, -0.160, -0.300 /

        true_time = time

!       Interpolate in time to get theta and rt tendency

        i1 = int( ( true_time - 41400. ) / 10800. ) + 1
        i1 = min( max( i1, 1 ), 5 )
        i2 = i1 + 1

        if (i1.lt.5) then
           a = ( true_time - (41400. + 10800. * (i1-1)) ) / 10800.
        else
           a = ( true_time - (41400. + 10800. * (i1-1)) ) / 9000.
        endif

        theta_tmp = ( 1. - a ) * ( atheta(i1) + rtheta(i1) )
     +                + a * ( atheta(i2) + rtheta(i2) )

        rt_tmp = ( 1. - a ) * art(i1) + a * art(i2)

!       Convert to the right units

        theta_tmp = theta_tmp / 3600.
        rt_tmp = rt_tmp / ( 3600. * 1000. )

!       Interpolate with respect to height

        do i=2,gr%nnzp

           if ( gr%zt(i).ge.0. .and. gr%zt(i).lt.1000. ) then
              theta_tndcy(i) = theta_tmp
              rt_tndcy(i) = rt_tmp
           elseif ( gr%zt(i).ge.1000. .and. gr%zt(i).lt.3000. ) then
              b = 1. - ( gr%zt(i) - 1000. ) / 2000.
              theta_tndcy(i) = b * theta_tmp
              rt_tndcy(i) = b * rt_tmp
           else
              theta_tndcy(i) = 0.
              rt_tndcy(i) = 0.
           endif

        enddo

        return
        end subroutine arm_tndcy

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for BOMEX case
!
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

           if ( gr%zt(i).ge.0. .and. gr%zt(i).lt.1500. ) then
              wmt(i) = - ( 0.0065 / 1500. ) * gr%zt(i)
           elseif ( gr%zt(i).ge.1500. .and. gr%zt(i).lt.2100. ) then
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

           if ( gr%zt(i).ge.0. .and. gr%zt(i).lt.1500. ) then
              theta_tndcy(i) = -2.315e-5
           elseif ( gr%zt(i).ge.1500. .and. gr%zt(i).lt.2500. ) then
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

           if ( gr%zt(i).ge.0. .and. gr%zt(i).lt.300. ) then
              rt_tndcy(i) = -1.2e-8
           elseif ( gr%zt(i).ge.300. .and. gr%zt(i).lt.500. ) then
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
!
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

           if ( gr%zt(i).ge.0. .and. gr%zt(i).lt.1500. ) then
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

        rt_tndcy = 0.

        return
        end subroutine fire_tndcy

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for Wangara case
!
        subroutine wangara_tndcy( time, wmt, wmm, theta_tndcy, 
     .                            rt_tndcy )

        use grid_class

        implicit none

!       Input
        real time

!       Output
        real wmt(gr%nnzp)
        real wmm(gr%nnzp)
        real theta_tndcy(gr%nnzp)
        real rt_tndcy(gr%nnzp)

!       No large-scale subsidence for now

        wmt = 0.
        wmm = 0.

!       No large-scale water tendendy or cooling

        rt_tndcy = 0.
        theta_tndcy = 0.

        return
        end subroutine wangara_tndcy

!----------------------------------------------------------------------
! Subroutine to set theta and water tendencies for ATEX case
!
        subroutine atex_tndcy(time, rtm, wmt, wmm, theta_tndcy, 
     .                        rt_tndcy, isValid)

        use grid_class

        implicit none

!       Input
        real, intent(in) :: time
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

        wmt = 0.
        wmm = 0.
        theta_tndcy = 0.
        rt_tndcy = 0.

!        if ( time >= 5400. ) then
        if ( time >= 48600. ) then   ! Brian

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
        call arm_sfcflx(true_time, heat_flx, moisture_flx)

! Compute momentum fluxes

!       Convert heat_flx and moisture_flx to natural units
        heat_flx2 = heat_flx / ( Cp * dn0 )        ! (K m/s)
        moisture_flx2 = moisture_flx / ( Lv * dn0 )! (m/s)

!       Heat flux in units of (m2/s3) (needed by diag_ustar)
        bflx = grav/thp * heat_flx2

!       Sfc winds
        usfc = up
        vsfc = vp
        ubar = max(ubmin,sqrt(usfc ** 2 + vsfc ** 2))

!       Compute ustar
        ustar = diag_ustar(z,bflx,ubar)

!       Assign fluxes

        uw = -usfc/ubar * ustar * ustar
        vw = -vsfc/ubar * ustar * ustar
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

        real ubmin
        parameter (ubmin = 0.25)
        real ustar
        parameter (ustar = 0.28)

!       Input variables

        real up,vp

!       Output variables

        real uw,vw,wtp,wrp

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

        wtp = -C_10 * ubar * ( thp - SST * (1000./1015.)**kappa )
        wrp = -C_10 * ubar * ( rtp - 0.0198293 )

! Compute momentum fluxes

        uw = -up * ustar*ustar / ubar
        vw = -vp * ustar*ustar / ubar

        return
        end subroutine atex_sfclyr

!----------------------------------------------------------------------
! This subroutine computes surface fluxes of horizontal momentum,
! heat and moisture for Wangara day 33
!
        subroutine wangara_sfclyr( time, up, vp, uw, vw, wtp, wrp )

        implicit none

!       Constants

        real ubmin
        parameter (ubmin = 0.25)
        real ustar
        parameter (ustar = 0.13)

!       Input variables

        real time
        real up,vp

!       Output variables

        real uw,vw,wtp,wrp

!       Internal variables

        real ubar
        real time_utc,time_est

! Compute UTC time of the day in seconds

        time_utc = mod( time, 86400. )

! Now convert UTC time to Australia EST (local time)

        time_est = mod( time_utc + 36000., 86400. )

        if ( time_est.lt.27000 .or. time_est.gt.63000 ) then
           write(*,*)
     .   'wangara_sfclyr: error local time must be between 730 and 1730'
           write(*,*) 'time_est = ',time_est
           stop
        endif

! Compute heat and moisture fluxes

        wtp = 0.18
     .        * cos( (time_est-45000.)/36000. * 3.141592654 )
        wrp = 1.3e-4 * wtp

! Compute momentum fluxes

        ubar = max(ubmin,sqrt(up*up + vp*vp))

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

        if ( time .le. times(1) ) then
           heat_flx = H(1)
           moisture_flx = LE(1)
        else if ( time .ge. times(ntimes) ) then
           heat_flx = H(ntimes)
           moisture_flx = LE(ntimes)
        else
           i1 = 1
           do while ( i1 .le. ntimes-1 )
              i2 = i1 + 1
              if ( time.ge.times(i1) .and. time.lt.times(i2) ) then
                 a = (time-times(i1))/(times(i2)-times(i1))
                 heat_flx = ( 1. - a ) * H(i1) + a * H(i2)
                 moisture_flx = ( 1. - a ) * LE(i1) + a * LE(i2)
                 i1 = ntimes
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

        htrt = htrt - 1. / ( Cp * rhot ) * ddzm( Frad )
        htrt(1) = 0.
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

        wmt = 0.
        wmm = 0.
        theta_tndcy = 0.
        rt_tndcy = 0.

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

        htrt = - 1. / ( Cp * rhot ) * ddzm( Frad )
        htrt(1) = 0.
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

        REAL, PARAMETER:: ls_div = 3.75e-6
        REAL, PARAMETER:: kap = 85.0 ! m^2/kg
        REAL, PARAMETER:: F0 = 70.0  ! W/m^2
        REAL, PARAMETER:: F1 = 22.0  ! W/m^2

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
        real, parameter :: Cd = 0.0011

!       Input variables

        integer, intent(in) :: sfctype
        real, intent(in)    :: Tsfc,psfc,exnersfc,up,vp,thp,rtp,rhop

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
