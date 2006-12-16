#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: dycoms2_rf01.F,v 1.1 2006-12-16 00:32:45 dschanen Exp $
        module dycoms2_rf01

!       Description:
!       Contains subroutines for the DYCOMS II RF01 case.
!----------------------------------------------------------------------
        public :: dycoms2_rf01_tndcy, dycoms2_rf01_sfclyr

        contains

!----------------------------------------------------------------------
        subroutine dycoms2_rf01_tndcy( time, rhot, rhom, rtm, rcm,
     .                                 exner, wmt, wmm, Frad, radht, 
     .                               thlm_forcing, rtm_forcing, isValid,
     .                                 sclrm, sclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for DYCOMS RF01 case.

!       References:
!----------------------------------------------------------------------

        use grid_class
        use constants
        use model_flags, only: lbugsrad
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

        intrinsic :: exp, sqrt

        ! Constant Parameters
        real, parameter :: 
     .  lsdiv =  3.75e-6,
     .  F0 = 70.0, F1 = 22.0, 
     .  kay = 85.0

        ! Input Variables
        real, intent(in) :: time ! Current time         [s]

        real, dimension(gr%nnzp), intent(in) :: 
     .  rhom, ! Density on moment. grid         [kg/m^3]
     .  rhot, ! Density on thermo. grid         [kg/m^3] 
     .  rtm,  ! Total water mixing ratio        [kg/kg]
     .  rcm,  ! Cloud water mixing ratio        [kg/kg]
     .  exner ! Exner function                  [-]

        ! Input/Output Variables
        logical, intent(inout) :: isValid

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) :: 
     .  wmt,          ! w wind on thermodynamic grid                 [m/s]
     .  wmm,          ! w wind on momentum grid                      [m/s]
     .  thlm_forcing, ! Liquid water potential temperature tendency  [K/s]
     .  rtm_forcing,  ! Total water mixing ratio tendency            [kg/kg/s]
     .  radht,        ! Radiative heating rate                       [K/s]
     .  Frad          ! Radiative flux                               [W/m^2]

        ! Input (optional)
        real, optional, intent(in) :: 
     .  sclrm(gr%nnzp, sclr_dim)

        ! Output (optional)
        real, optional, intent(out), dimension(gr%nnzp, sclr_dim) :: 
     .  sclrm_forcing

        ! Internal variables
        real, dimension(gr%nnzp) :: lwp

        integer :: i, zi

        wmt          = 0.
        wmm          = 0.
        thlm_forcing = 0.
        rtm_forcing  = 0.

        ! Identify height of 8.0 g/kg moisture level

        i = 2
        do while ( i <= gr%nnzp .and. rtm(i) > 8.0e-3 )
           i = i + 1
        end do
        if ( i == gr%nnzp+1 .or. i == 2 ) then
          write(fstderr,*) "Identification of 8.0 g/kg level failed"
          write(fstderr,*) "Subroutine: dycoms2_rf01_tndcy.  "
     .      //" File: dycoms2_rf01.F"
          write(fstderr,*) "i = ",i
          write(fstderr,*) "rtm(i) = ",rtm(i)
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
        
        ! Theta-l radiative tendency

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

        ! Compute IR heating rate

          radht          = ( -1.0/(Cp*rhot) ) * ddzm( Frad )
     .                   * ( 1.0 / exner )
          radht(1)       = 0.
          radht(gr%nnzp) = 0.

        end if ! ~ lbugsrad

        ! Add heating rate to theta-l forcing

        if ( .not. lbugsrad ) thlm_forcing = thlm_forcing + radht
        
        if (present( sclrm_forcing )) sclrm_forcing(:,SCLR_THETA)
     .  = thlm_forcing
        if (present( sclrm_forcing )) sclrm_forcing(:,SCLR_RT)
     .  = rtm_forcing

        return
        end subroutine dycoms2_rf01_tndcy
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

        end module dycoms2_rf01
