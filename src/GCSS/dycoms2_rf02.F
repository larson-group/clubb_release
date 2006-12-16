#define SCLR_THETA 1
#define SCLR_RT 2
!----------------------------------------------------------------------
! $Id: dycoms2_rf02.F,v 1.1 2006-12-16 00:32:45 dschanen Exp $
        module dycoms2_rf02

!       Description:
!       Contains subroutines for the DYCOMS II RF02 case.
!----------------------------------------------------------------------

        public :: dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr

        contains

!----------------------------------------------------------------------
        SUBROUTINE dycoms2_rf02_tndcy
     .             ( time, time_initial, rhot, rhom, rtm, rcm, exner, 
     .               wmt, wmm, theta_tndcy, rt_tndcy, 
     .               Frad, radht, Ncm, Ncnm, isValid, 
     .               sclr_tndcy )
!       Description:

!       References:
!----------------------------------------------------------------------

        USE grid_class
        USE constants
        USE model_flags, only: lbugsrad, lcoamps_micro, kk_rain
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

        if ( PRESENT( wpsclr_sfc ) ) then
          wpsclr_sfc(SCLR_THETA) = SH / (1.21 * Cp)

          wpsclr_sfc(SCLR_RT) = LH / (1.21 * Lv)
        end if

        RETURN
        END SUBROUTINE dycoms2_rf02_sfclyr

        end module dycoms2_rf02
