! $Id: ice_dfsn_mod.F,v 1.3 2008-04-05 19:54:32 griffinb Exp $        

        module ice_dfsn_mod

        public  :: ice_dfsn
        private :: Diff_denom

        contains
!-----------------------------------------------------------------------
        SUBROUTINE ice_dfsn( dt, thlm, rcm, press, exner, rhot,
     .                       rcm_icedfsn )
        ! Description:
        ! This subroutine is based on a COAMPS subroutine (nov11_icedfs)
        ! written by Adam Smith and Vince Larson to calculate the
        ! depletion of cloud water by the diffusional growth of ice.

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
        
        !     References:
        !       Section 4.2 of Larson et al. (2006), "What determines altocumulus
        !         dissipation time?", J. Geophys. Res., Vol. 111, D19207.
        !       Mitchell, D. L. (1996), "Use of mass- and area- ...", J. Atmos. Sci.
        !         Vol. 53, 1710--1723.        
        !       Rogers and Yau (1989), "A Short Course in Cloud Physics", 3rd. Ed. 
        !       Fleishauer et al. (2002), "Observed microphysical structure of
        !         midlevel, mixed-phase clouds", J. Atmos. Sci., Vol. 59, 
        !         pp. 1779--1804.


        USE grid_class
        USE constants
        USE stats_prec, only: time_prec
        USE saturation, only: rsat
#ifdef STATS
        use stats_hoc, only: zt, lstats_samp,
     .      ircm_icedfs, idiam, imass_ice_cryst, iu_T_cm
#endif

        implicit none

        ! Input variables
        REAL(KIND=time_prec), INTENT(IN):: 
     .  dt      ! Model timestep                                     [s]

        REAL, DIMENSION(1:gr%nnzp), INTENT(IN):: 
     .  thlm,   ! Mean liquid water potential temperature            [K]
     .  rcm,    ! Cloud water mixing ratio                           [kg kg^{-1}]
     .  press,  ! Air pressure                                       [Pa]
     .  exner,  ! Exner function                                     []
     .  rhot    ! Air density on thermodynamic grid                  [kg m^{-3}]

        ! Output variables
        REAL, DIMENSION(1:gr%nnzp), INTENT(OUT):: 
     .  rcm_icedfsn   ! Time tendency of rcm due to ice diffusional growth  
     .                !                                               [kg kg^{-1} s^{-1}]

        ! Local variables
        REAL, DIMENSION(1:gr%nnzp):: 
     .  mass_ice_cryst, ! Mass of a single ice crystal      [kg]
     .  theta,           ! Potential temperature                       [K]
     .  temp_in_K,       ! Temperature                                 [K]
     .  r_s,             ! Saturation mixing ratio over vapor          [kg kg^{-1}] 
     .  e_s,             ! Saturation vapor pressure over liquid       [Pa]
     .  e_i,             ! Saturation vapor pressure over ice          [Pa]
     .  S_i,             ! Ratio of saturation w.r.t. liquid to that for ice []
     .  Denom,           ! Denominator of diffusional growth equation  [m s kg^{-1}] 
     .  dmass_ice_cryst, ! Change in ice mass over vertical grid box   [kg m^{-1}]
     .  diam,            ! Diameter of ice crystal                     [m]
     .  u_T_cm           ! Fallspeed of ice crystal in cm/s            [cm s^{-1}]

        REAL:: 
     .  a_coef,    ! Pre-factor for mass-diameter relationship, Mitchell (1996) [kg] 
     .  b_expn,    ! Exponential for mass-diameter relationship, Mitchell (1996) []
     .  k_u_coef,  ! Pre-factor for fallspeed-diameter formula                  [m s^{-1}]
     .  q_expn,    ! Exponential of density in fallspeed-diameter formula       []   
     .  n_expn     ! Exponential of diameter in fallspeed-diameter formula      []
        ! Number of ice crystals per unit volume of air    [m^{-3}]
        ! Vince Larson avgd legs 2 and 7 (Fleishauer et al)  21 Jan 2005
        REAL, PARAMETER:: N_i = 2000

        INTEGER:: k

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                                                                     !
        ! Coefficients for mass-diameter relationship, Mitchell (1996)        !
        ! mass = a (diameter/(1 meter))^b,  [a] = kg, [b] = []                !
        !                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        a_coef = 2.05e-3
        b_expn = 1.8

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

        k_u_coef = 55.
        q_expn   = 0.17
        n_expn   = 0.70

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
           mass_ice_cryst(k) = 1.0e-11
        END DO

        DO k = gr%nnzp, 2, -1

          !  Compute Temperature [K]
          theta(k) = thlm(k) + ( Lv / (Cp*exner(k)) ) * rcm(k)
          temp_in_K(k) = theta(k) * exner(k)

          ! Check whether we're in cloud and below freezing.
          ! Note:  A value of 1.0E-5 kg/kg is used as a threshold value
          ! for rcm because the HOC model shows a small amount of liquid
          ! water all the way to the model top, which messes with the
          ! ice diffusion calculations.
          IF ( rcm(k) >= 1.0E-5 .AND. temp_in_K(k) < 273.15 ) THEN

          ! Find saturation mixing ratio over vapor [kg kg^{-1}]
             r_s(k) = rsat( press(k), temp_in_K(k) )

          ! Saturation vapor pressure over liquid in Pa
             e_s(k) = ( r_s(k)*press(k) ) / ( ep + r_s(k) )

          ! Saturation vapor pressure over ice in Pa, Eq. 2.15 Rogers and Yau
             e_i(k) = e_s(k) / EXP( ( Lf/(Rv*273.16) )
     .                           * ( 273.16/temp_in_K(k) - 1.0 ) )

          ! Saturation ratio in a liquid-saturated cloud, p. 158 Rogers and Yau
          !---------------Brian's comment--------------------------------------!
          ! The actual formula is:  Si = e/ei = (e/es)*(es/ei) = S*(es/ei)     !
          ! It is assumed that any supersaturation forms liquid water and that !
          ! the atmosphere is then saturated with respect to liquid water.     !
          ! Therefore, S = 1.0, allowing Si = es/ei.                           !
          !--------------------------------------------------------------------!
             S_i(k) = e_s(k)/e_i(k)

          ! Denominator of diffusional growth equation, 9.4 of Rogers and Yau
             Denom(k) = Diff_denom( temp_in_K(k), press(k), e_i(k) )

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
             rcm_icedfsn(k) = - (N_i/rhot(k))
     .          * ( 4 * (S_i(k) - 1) / Denom(k) )
     .          * (mass_ice_cryst(k)/a_coef)**(1/b_expn)

             ! Ensure that liquid is not over-depleted
             IF ( rcm(k) + rcm_icedfsn(k)*dt < 0.0 ) THEN
               rcm_icedfsn(k) = -rcm(k)/dt
             END IF

             !---------------Brian's comment-----------------------------------!
             ! dm = (dm/dt)*(dt/dz)*dz                                         !
             ! dm = (dm/dt)*(1/u_T)*dz                                         !
             !-----------------------------------------------------------------!
             dmass_ice_cryst(k) = ( 4 * (S_i(k) - 1) / Denom(k) )
     .          * (k_u_coef**(-1.0)) * ( rhot(k)**q_expn )
     .          * ( (mass_ice_cryst(k)/a_coef)**((1.0-n_expn)/b_expn) )
     .          * (1.0/gr%dzm(k-1))
             mass_ice_cryst(k-1) = mass_ice_cryst(k) 
     .                                    + dmass_ice_cryst(k)

             ! Diameter of ice crystal in meters.
             diam(k) = (mass_ice_cryst(k)/a_coef)**(1/b_expn)

             ! Fallspeed of ice crystal in cm/s.
             u_T_cm(k) = 100. * k_u_coef *
     .                   ((mass_ice_cryst(k)/a_coef)**(n_expn/b_expn)) 
     .                         * (rhot(k)**(-q_expn))

          ELSE   ! There's no liquid and/or ice present; assume no ice growth

             mass_ice_cryst(k-1) = mass_ice_cryst(k)
             rcm_icedfsn(k) = 0.0
             diam(k)        = 0.0  ! Set zero to remind that we don't grow ice
             u_T_cm(k)      = 0.0  ! Set zero to remind that we don't grow ice

          END IF

        END DO

! Michael Falk added boundary condx, 31 July 2006

        mass_ice_cryst(1) = mass_ice_cryst(2)
        rcm_icedfsn(1)    = rcm_icedfsn(2)
        diam(1)           = diam(2)
        u_T_cm(1)         = u_T_cm(2)

! eMFc

!
#ifdef STATS
        if ( lstats_samp ) then
!       diam(:) ! Icedfs diameter; Michael Falk, 1 Nov 2006
!       m(:)    ! Icedfs mass; Michael Falk, 1 Nov 2006
!       dqc_dt_icedfs(:) ! Icedfs change in liquid; Michael Falk, 1 Nov 2006
!       u_T_cm(:)        ! Icedfs fallspeed (cm/s); Michael Falk, 1 Nov 2006

          if ( ircm_icedfs > 0 ) then
            zt%x(:,ircm_icedfs) = zt%x(:,ircm_icedfs)
     .                          + rcm_icedfsn
            zt%n(:,ircm_icedfs) = zt%n(:,ircm_icedfs) + 1
          end if
          if ( idiam > 0 ) then
            zt%x(:,idiam) = zt%x(:,idiam) + diam
            zt%n(:,idiam) = zt%n(:,idiam) + 1
          end if
          if ( imass_ice_cryst > 0 ) then
            zt%x(:,imass_ice_cryst) = 
     .               zt%x(:,imass_ice_cryst) + mass_ice_cryst
            zt%n(:,imass_ice_cryst) = zt%n(:,imass_ice_cryst) + 1
          end if
          if ( iu_T_cm > 0 ) then
            zt%x(:,iu_T_cm) = zt%x(:,iu_T_cm) + u_T_cm
            zt%n(:,iu_T_cm) = zt%n(:,iu_T_cm) + 1
          end if
        end if
#endif /*STATS*/

        RETURN
        END SUBROUTINE ice_dfsn

!-----------------------------------------------------------------------

        FUNCTION Diff_denom( temp_in_K, press, e_i )

        USE constants

        IMPLICIT NONE

        ! Compute denominator of diffusional growth equation

        ! Reference:  Eqn. 9.4 of Rogers and Yau (1989), "A Short Course on Cloud Physics"

        REAL, INTENT(IN) :: 
     .   temp_in_K,   ! Temperature                               [K]
     .   press,       ! Air pressure                              [Pa]
     .   e_i          ! Vapor pressure over ice                   [Pa]

        REAL :: 
     .  Diff_denom   ! Denominator of diffusional growth equation  [m s kg^{-1}]

        REAL:: Ka, Dv
        REAL:: Fk, Fd
        REAL:: Celsius

!        REAL, PARAMETER:: Ls = 2.834e6

        Celsius = temp_in_K - 273.16

        Ka = (5.69 + 0.017*Celsius)*0.00001  ! Ka in cal./(cm.*sec.*C)
        Ka = 4.1868*100.0*Ka  ! Ka in J./(m.*sec.*K)

        Dv = 0.221 * ( (temp_in_K/273.16)**1.94 ) * (101325.0/press)
                                ! Dv in (cm.^2)/sec.  ! .221 is correct.
        Dv = Dv/10000.0  ! Dv in (m.^2)/sec.

        Fk = ( Ls/(Rv*temp_in_K) - 1.0 ) * Ls / (Ka*temp_in_K)
        Fd = (Rv*temp_in_K) / (Dv*e_i)

        Diff_denom = Fk + Fd

        RETURN
        END FUNCTION Diff_denom

        end module ice_dfsn_mod
