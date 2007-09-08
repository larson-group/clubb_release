#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: dycoms2_rf02.F,v 1.7 2007-09-08 18:14:51 griffinb Exp $
        module dycoms2_rf02

!       Description:
!       Contains subroutines for the DYCOMS II RF02 case.
!----------------------------------------------------------------------

        public :: dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr

        contains

!----------------------------------------------------------------------
        SUBROUTINE dycoms2_rf02_tndcy
     .             ( time, time_initial, rhot, rhom, rtm, rcm, exner, 
     .               wmt, wmm, thlm_forcing, rtm_forcing, 
     .               Frad, radht, Ncm, Ncnm, isValid, 
     .               sclrm_forcing )
!       Description:
!       Compute wm, thlm_ls, rtm_ls, radiative heating rate, and cloud
!       droplet number concentration as needed.

!       References:
!       ``Dynamics and Chemistry of Marine Stratocumulus -- DYCOMS-II''
!       Stevens, Bjorn, et al., (2003) 
!       Bull. Amer. Meteorol. Soc., 84, 579-593.
!----------------------------------------------------------------------

        USE grid_class
        USE constants
        USE model_flags, only: lbugsrad, lcoamps_micro, kk_rain
#ifdef STATS
        USE hoc_stats
#endif /*STATS*/

        ! Constant parameters
        real, parameter :: 
     .  ls_div = 3.75e-6,
     .  kap    = 85.0,    ! [m^2/kg]
     .  F0     = 70.0,    ! [W/m^2]
     .  F1     = 22.0     ! [W/m^2]

        ! Input Variables
        real, intent(in) ::
     .  time,         ! Current time    [s]
     .  time_initial  ! Initial time    [s]

        real, intent(in), dimension(gr%nnzp) ::
     .  rhot,          ! Density on thermo. grid        [kg/m^3] 
     .  rhom,          ! Density on moment. grid        [kg/m^3]
     .  rtm,           ! Total water mixing ratio       [kg/kg]
     .  rcm,           ! Cloud water mixing ratio       [kg/kg]
     .  exner          ! Exner function.                [-]

        ! Input/Output Variables
        logical, intent(inout) :: isValid

        ! Output Variables
        real, intent(out), dimension(gr%nnzp) :: 
     .  wmt,           ! wm on thermodynamic grid       [m/s]
     .  wmm,           ! wm on momentum grid            [m/s]
     .  thlm_forcing,  ! theta_l forcing                [K/s]
     .  rtm_forcing,   ! r_t forcing                    [(kg/kg)/s] 
     .  Frad,          ! Radiative flux                 [W/m^2]
     .  radht,         ! Radiative heating rate         [K/s] 
     .  Ncm,           ! Cloud droplet number conc.     [#/kg]
     .  Ncnm           ! Cloud nuclei number conc.      [#/m^3]

        ! Optional Output Variables
        real, intent(out), optional, dimension(gr%nnzp,sclr_dim) ::
     .  sclrm_forcing  ! Passive scalar tendency        [units/s]

        ! Local Variables
        real, dimension(gr%nnzp) :: 
     .  LWP,       ! Liquid water path
     .  Heaviside

        real :: z_i

        integer :: k  ! Loop index

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

         ! Compute liquid water path from top of the model
         ! We define liquid water path on momentum levels

          LWP(gr%nnzp) = 0.0
          DO k = gr%nnzp-1, 1, -1
            LWP(k) = LWP(k+1) + rhot(k+1)*rcm(k+1)/gr%dzt(k+1)
          END DO  ! k = gr%nnzp..1

!         Find the height of the isotherm rtm = 8.0 g/kg.

          k = 2
          DO WHILE ( k <= gr%nnzp .AND. rtm(k) > 8.0e-3 )
            k = k + 1
          END DO
          IF ( k == gr%nnzp+1 .or. k == 2 ) THEN
            write(fstderr,*) "Identification of 8.0 g/kg level failed"
            write(fstderr,*) "Subroutine: dycoms2_rf02_tndcy. "
     .        // "File: dycoms2_rf02.F"
            write(fstderr,*) "k = ", k
            write(fstderr,*) "rtm(k) = ", rtm(k)
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

#ifdef STATS
          if ( lstats_samp .and. iradht_LW > 0 ) then
            zt%x(:,iradht_LW) = zt%x(:,iradht_LW) + radht
            zt%n(:,iradht_LW) = zt%n(:,iradht_LW) + 1
          end if
#endif
        END IF ! ~ lbugsrad

        ! Enter the final theta-l and rtm tendencies

        IF ( .not. lbugsrad ) THEN
          thlm_forcing(1:gr%nnzp) = radht(1:gr%nnzp)
        ELSE
          thlm_forcing(1:gr%nnzp) = 0.0
        END IF

        rtm_forcing(1:gr%nnzp) = 0.0

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
        ! Since cloud base (zb) is determined by the mixing ratio rc_tol,
        ! so will cloud droplet number concentration (Ncm).

        if ( lcoamps_micro .and. time == time_initial ) then

          ! Taken from COAMPS subroutine ncn_init()
          Ncnm(1:gr%nnzp) = 55000000.0 / rhot(1:gr%nnzp)

        else

          ! K & K or no microphysics
          DO k = 1, gr%nnzp, 1
            IF ( rcm(k) > rc_tol ) THEN
              ! The specified cloud droplet concentration is 55 cm^-3, which
              ! is then converted to m^-3, and then divided by rhot to get the
              ! concentration in units of kg^-1.
              Ncm(k) = 55000000.0 / rhot(k)
            ELSE
              Ncm(k) = 0.0
            END IF
          END DO ! k=1..gr%nnzp

        end if

        if ( present( sclrm_forcing ) )  then
          sclrm_forcing(:,SCLR_THETA) = thlm_forcing(:) 
          sclrm_forcing(:,SCLR_RT)    = rtm_forcing(:)
        end if

        RETURN
        END SUBROUTINE dycoms2_rf02_tndcy


!----------------------------------------------------------------------

        subroutine dycoms2_rf02_sfclyr( um_sfc, vm_sfc, 
     .                                  upwp_sfc, vpwp_sfc, 
     .                                  wpthlp_sfc, wprtp_sfc, ustar,
     .                                  wpsclrp_sfc, wpedsclrp_sfc )

        use constants

        ! External
        intrinsic :: sqrt, present

        ! Constant parameters
        real, parameter :: 
     .  SH = 16.0,
     .  LH = 93.0

        ! Input Variables
        real, intent(in) :: 
     .  um_sfc, ! um(2) [m/s]
     .  vm_sfc  ! vm(2) [m/s]

        ! Output
        real, intent(out) :: 
     .  upwp_sfc,    ! u'w' at (1)      [m^2/s^2]
     .  vpwp_sfc,    ! v'w'at (1)       [m^2/s^2]
     .  wpthlp_sfc,  ! w'th_l' at (1)   [(m K)/s]  
     .  wprtp_sfc,   ! w'r_t'(1) at (1) [(m kg)/(s kg)]
     .  ustar

        ! Output Variables (optional)
        real, intent(out), optional, dimension(sclr_dim) :: 
     .  wpsclrp_sfc,   ! w' scalar at surface [units m/s]
     .  wpedsclrp_sfc  ! w' scalar at surface [units m/s]

        ! Local Variables
        real :: wind_sfc  ! ? [m^2/s^2]?

        ustar = 0.25

        wind_sfc = sqrt( um_sfc**2 + vm_sfc**2 )

        if (wind_sfc > 0.0) then
          upwp_sfc = -um_sfc * ( ustar**2 ) / wind_sfc
          vpwp_sfc = -vm_sfc * ( ustar**2 ) / wind_sfc

        else
           upwp_sfc = 0.0
           vpwp_sfc = 0.0

        end if  ! wind_sfc > 0

        wpthlp_sfc = SH / (1.21 * Cp)
        wprtp_sfc  = LH / (1.21 * Lv)

        ! Let passive scalars be equal to rt and theta_l for now
        if ( present( wpsclrp_sfc ) ) then
          wpsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpsclrp_sfc(SCLR_RT)    = wprtp_sfc
        end if

        if ( present( wpedsclrp_sfc ) ) then
          wpedsclrp_sfc(SCLR_THETA) = wpthlp_sfc
          wpedsclrp_sfc(SCLR_RT)    = wprtp_sfc
        end if

        return
        end subroutine dycoms2_rf02_sfclyr

        end module dycoms2_rf02
