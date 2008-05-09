#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: dycoms2_rf01.F,v 1.13 2008-05-09 22:39:04 faschinj Exp $
        module dycoms2_rf01

!       Description:
!       Contains subroutines for the DYCOMS II RF01 case.
!----------------------------------------------------------------------
        public :: dycoms2_rf01_tndcy, dycoms2_rf01_sfclyr

        contains

!----------------------------------------------------------------------
        subroutine dycoms2_rf01_tndcy
     .             ( time, rhot, rhom, rtm, rcm, exner,
     .               wmt, wmm, Frad, radht, thlm_forcing,
     .               rtm_forcing, err_code,
     .               sclrm, sclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for DYCOMS RF01 case.

!       References:
!----------------------------------------------------------------------

        use grid_class
        use constants
        use model_flags, only: lbugsrad
        use stats_prec, only: time_prec
        use error_code, only: clubb_rtm_level_not_found
#ifdef STATS
        use stats_hoc
#endif /*STATS*/

        implicit none

        intrinsic :: exp, sqrt, present

        ! Constant Parameters
        real, parameter :: 
     .  lsdiv =  3.75e-6,
     .  F0 = 70.0, F1 = 22.0, 
     .  kay = 85.0

        ! Input Variables
        real(kind=time_prec), intent(in) :: time ! Model time [s]

        real, dimension(gr%nnzp), intent(in) :: 
     .  rhom, ! Density on moment. grid         [kg/m^3]
     .  rhot, ! Density on thermo. grid         [kg/m^3] 
     .  rtm,  ! Total water mixing ratio        [kg/kg]
     .  rcm,  ! Cloud water mixing ratio        [kg/kg]
     .  exner ! Exner function                  [-]

        ! Input/Output Variables
        integer, intent(inout) :: err_code

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
          err_code = clubb_rtm_level_not_found
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

          ! Compute liquid water path from top of the model
          ! We define liquid water path on momentum levels

          lwp(gr%nnzp) = 0.0
          do i = gr%nnzp-1, 1, -1
            lwp(i) = lwp(i+1) + rhot(i+1) * rcm(i+1) / gr%dzt(i+1)
          end do
!         x_sfc(1,ilwp) = lwp(1)

          ! Compute IR radiative flux

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

#ifdef STATS
          if ( lstats_samp .and. iradht_LW > 0 ) then
            zt%x(:,iradht_LW) = zt%x(:,iradht_LW) + radht
            zt%n(:,iradht_LW) = zt%n(:,iradht_LW) + 1
          end if
#endif
        end if ! ~ lbugsrad

        ! Add heating rate to theta-l forcing

        if ( .not. lbugsrad ) thlm_forcing = thlm_forcing + radht
        
        if (present( sclrm_forcing )) sclrm_forcing(:,SCLR_THETA)
     .  = thlm_forcing(:)
        if (present( sclrm_forcing )) sclrm_forcing(:,SCLR_RT)
     .  = rtm_forcing(:)

        return
        end subroutine dycoms2_rf01_tndcy

!----------------------------------------------------------------------
        subroutine dycoms2_rf01_sfclyr
     .             ( sfctype, Tsfc, psfc, 
     .               exner_sfc, um_sfc, vm_sfc, 
     .               thlm_sfc, rtm_sfc, 
     .               rhom_sfc, upwp_sfc, vpwp_sfc,
     .               wpthlp_sfc, wprtp_sfc, ustar,
     .               sclrm_sfc, wpsclrp_sfc, wpedsclrp_sfc )
!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS DYCOMS II RF 01 specifications

!       References:
!----------------------------------------------------------------------

        use constants
        use saturation, only: rsat

#ifdef STATS
        use stats_hoc
#endif /*STATS*/

        implicit none

        intrinsic :: max, sqrt, present

        ! Constant parameters
        real, parameter :: 
     .  ubmin = 0.25,
!     .  ustar = 0.25,
     .  Cd    = 0.0011

        ! Input variables
        integer, intent(in) :: sfctype

        real, intent(in) :: 
     .  Tsfc,      ! Surface temperature                [K]
     .  psfc,      ! Surface pressure                   [Pa]
     .  exner_sfc, ! Exner function at the surface      [-]
     .  um_sfc,    ! u wind at the surface              [m/s]
     .  vm_sfc,    ! v wind at the surface              [m/s]
     .  thlm_sfc,  ! theta_l at the surface             [K]
     .  rtm_sfc,   ! r_t at the surface                 [kg/kg]
     .  rhom_sfc   ! Density at the surface             [kg/m^3]

        ! Optional input variables
        real, intent(in), optional, dimension(sclr_dim) :: 
     .  sclrm_sfc  ! Passive scalar at the surface      [units vary]

        ! Output variables
        real, intent(out) :: 
     .  upwp_sfc,   ! u' w' at the surface              [m^2/s^2]
     .  vpwp_sfc,   ! v' w' at the surface              [m^2/s^2]
     .  wpthlp_sfc, ! w' thl' at the surface            [m K/s]
     .  wprtp_sfc,  ! w' rt'  at the surface            [m kg/kg]
     .  ustar       ! surface friction velocity         [m/s]

        real, intent(out), optional, dimension(sclr_dim) :: 
     .  wpsclrp_sfc,  ! w' sclr' at the surface         [m units/s]
     .  wpedsclrp_sfc ! w' edsclr' at the surface       [m units/s]

        ! Local variables

        real :: ubar

        ! Declare the value of ustar.
        ustar = 0.25

        ! Compute heat and moisture fluxes

        ubar = max( ubmin, sqrt( um_sfc**2 + vm_sfc**2 ) )

        ! Compute momentum fluxes

        upwp_sfc = -um_sfc * ustar**2 / ubar
        vpwp_sfc = -vm_sfc * ustar**2 / ubar

        if ( sfctype == 0 ) then

          wpthlp_sfc =  15.0 / ( rhom_sfc * Cp )
          wprtp_sfc  = 115.0 / ( rhom_sfc * Lv )

        else if ( sfctype == 1 ) then

          wpthlp_sfc = -Cd * ubar * ( thlm_sfc - Tsfc/exner_sfc )
          wprtp_sfc  = -Cd * ubar * ( rtm_sfc  - rsat( psfc, Tsfc ) )

        else  ! Undefined value for sfctype

          write(fstderr,*) "Invalid sfctype value = ", sfctype
          stop

        end if

        ! Let passive scalars be rt and theta for now
        if ( present( wpsclrp_sfc ) ) then
          wpsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpsclrp_sfc(SCLR_RT)    = wprtp_sfc
        end if

        if ( present( wpedsclrp_sfc ) ) then
          wpedsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpedsclrp_sfc(SCLR_RT)    = wprtp_sfc
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
