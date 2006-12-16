        module icedfs

        public  :: altocu_icedf
        private :: Diff_denom

        contains
!-----------------------------------------------------------------------
        SUBROUTINE altocu_icedf( dt, thlm, rcm, p, exner, rhot,
     .                           dqc_dt_icedfs )
        ! Description:
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

        USE grid_class
        USE constants
        USE statistics

        implicit none

        ! For function call
        REAL, EXTERNAL :: rsat

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
        FUNCTION Diff_denom( Temp, press, ei )

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

        end module icedfs
