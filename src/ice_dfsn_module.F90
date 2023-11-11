! $Id$

module ice_dfsn_module

  implicit none

  public  :: ice_dfsn
  private :: Diff_denom

  private ! Default Scope

  contains
!-------------------------------------------------------------------------------
  subroutine ice_dfsn( gr, dt, thlm, rcm, exner, p_in_Pa, rho, & 
                       stats_metadata, &
                       stats_zt, &
                       rcm_icedfsn, thlm_icedfsn )
! Description:
!   This subroutine is based on a COAMPS subroutine (nov11_icedfs)
!   written by Adam Smith and Vince Larson to calculate the
!   depletion of cloud water by the diffusional growth of ice.
!
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
!-------------------------------------------------------------------------------

! References:
!   Section 4.2 of Larson et al. (2006), "What determines altocumulus
!     dissipation time?", J. Geophys. Res., Vol. 111, D19207.
!
!   Mitchell, D. L. (1996), "Use of mass- and area- ...", J. Atmos. Sci.
!     Vol. 53, 1710--1723.
!
!   Rogers and Yau (1989), "A Short Course in Cloud Physics", 3rd. Ed.
!
!   Fleishauer et al. (2002), "Observed microphysical structure of
!     midlevel, mixed-phase clouds", J. Atmos. Sci., Vol. 59,
!     pp. 1779--1804.
!-------------------------------------------------------------------------------

    use grid_class, only: grid

    use constants_clubb, only: & 
        Cp,  & ! Constant(s)
        Lv, & 
        ep, & 
        Rv, & 
        Lf,&
        T_freeze_K, &
        cm_per_m

    use clubb_precision, only:  & 
        core_rknd ! Variable(s)

    use saturation, only:  & 
        sat_mixrat_liq ! Procedure(s)

    use T_in_K_module, only: &
        thlm2T_in_K ! Procedure(s)

    use stats_type_utilities, only: & 
        stat_update_var

    use stats_variables, only: &
        stats_metadata_type

    use stats_type, only: &
        stats ! Type

    implicit none

    ! Constant Parameters
    ! Number of ice crystals per unit volume of air    [m^{-3}]
    ! Vince Larson avgd legs 2 and 7 (Fleishauer et al)  21 Jan 2005
    real(kind = core_rknd), parameter :: N_i = 2000._core_rknd

    !---------------------- Input variables ----------------------
    type(grid), target, intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in)::  & 
      dt      ! Model timestep                                     [s]

    real(kind = core_rknd), dimension(gr%nz), intent(in)::  & 
      thlm,    & ! Liquid potential temperature         [K]
      rcm,     & ! Cloud water mixing ratio             [kg kg^{-1}]
      exner,   & ! Exner function                       [-]
      p_in_Pa, & ! Air pressure                         [Pa]
      rho        ! Air density on thermodynamic grid    [kg m^{-3}]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !---------------------- InOut variables ----------------------
    type(stats), target, intent(inout) :: &
      stats_zt

    !---------------------- Output variables ----------------------
    real(kind = core_rknd), dimension(gr%nz), intent(out)::  & 
      rcm_icedfsn, & ! Time tendency of rcm due to ice diffusional growth  [kg kg^{-1} s^{-1}]
      thlm_icedfsn   ! Time tendency of thlm due to ice diffusional growth [K/s]

    ! Local variables
    real(kind = core_rknd), dimension(gr%nz)::  & 
      T_in_K,           & ! Absolute temperature                        [K]
      mass_ice_cryst,   & ! Mass of a single ice crystal                [kg]
      r_s,              & ! Saturation mixing ratio over vapor          [kg kg^{-1}] 
      e_s,              & ! Saturation vapor pressure over liquid       [Pa]
      e_i,              & ! Saturation vapor pressure over ice          [Pa]
      S_i,              & ! Ratio of saturation w.r.t. liquid to that for ice []
      Denom,            & ! Denominator of diffusional growth equation  [m s kg^{-1}] 
      dmass_ice_cryst,  & ! Change in ice mass over vertical grid box   [kg m^{-1}]
      diam,             & ! Diameter of ice crystal                     [m]
      u_T_cm              ! Fallspeed of ice crystal in cm/s            [cm s^{-1}]

    REAL(KIND=core_rknd)::  & 
      a_coef,     & ! Pre-factor for mass-diameter relationship, Mitchell (1996) [kg] 
      b_expn,     & ! Exponential for mass-diameter relationship, Mitchell (1996) []
      k_u_coef,   & ! Pre-factor for fallspeed-diameter formula                  [m s^{-1}]
      q_expn,     & ! Exponential of density in fallspeed-diameter formula       []   
      n_expn        ! Exponential of diameter in fallspeed-diameter formula      []

    integer :: k

    !---------------------- Begin Code ----------------------

    ! Determine absolute temperature
    T_in_K = thlm2T_in_K( gr%nz, thlm, exner, rcm )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                     !
    ! Coefficients for mass-diameter relationship, Mitchell (1996)        !
    ! mass = a (diameter/(1 meter))^b,  [a] = kg, [b] = []                !
    !                                                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    a_coef = 2.05e-3_core_rknd
    b_expn = 1.8_core_rknd

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                     !
    !  Coefficients for mass-diameter relationship, Kajikawa (1989)       !
    !  mass = a (diam/(1m))^b,  [a] = kg, [b] = []                        !
    !                                                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       a = 2.50e-4_core_rknd
!       b = 1.4_core_rknd

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                     !
    !  Coefficients for fallspeed-diameter relationship, Mitchell (1996)  !
    !  u_T = k_u rho^{-q} (diameter/(1 meter))^n,                         !
    !       [k_u] = m/s, [q] = [], [n] = [], [rho] = kg m^{-3}            !
    !                                                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    k_u_coef = 55._core_rknd
    q_expn   = 0.17_core_rknd
    n_expn   = 0.70_core_rknd

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                     !
    !  Coefficients for fallspeed-diameter relationship, Kajikawa (1989)  !
    !  u_T = k_u rho^{-q} (diam/(1m))^n,  [k_u] = m/s, [q] = [], [n] = [] !
    !       [rho] = kg m^{-3}                                             !
    !                                                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       k_u = 0.438_core_rknd
!       q = 0.0_core_rknd
!       n = 0.0742_core_rknd



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                     !
    !  Initialize ice particle mass                                       !
    !                                                                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do k = 1, gr%nz, 1
      mass_ice_cryst(k) = 1.0e-11_core_rknd
    end do

    do k = gr%nz, 2, -1

      ! Check whether we're in cloud and below freezing.
      ! Note:  A value of 1.0E-5 kg/kg is used as a threshold value
      ! for rcm because the CLUBB model shows a small amount of liquid
      ! water all the way to the model top, which messes with the
      ! ice diffusion calculations.
      IF ( rcm(k) >= 1.0E-5_core_rknd .AND. T_in_K(k) < T_freeze_K ) THEN

        ! Find saturation mixing ratio over vapor [kg kg^{-1}]
        r_s(k) = sat_mixrat_liq( p_in_Pa(k), T_in_K(k) )

        ! Saturation vapor pressure over liquid in Pa
        e_s(k) = ( r_s(k)*p_in_Pa(k) ) / ( ep + r_s(k) )

        ! Saturation vapor pressure over ice in Pa, Eq. 2.15 Rogers and Yau
        e_i(k) = e_s(k) / EXP( ( Lf/(Rv*T_freeze_K) ) & 
                            * ( T_freeze_K/T_in_K(k) - 1.0_core_rknd ) )

        ! Saturation ratio in a liquid-saturated cloud, p. 158 Rogers and Yau
        !---------------Brian's comment--------------------------------------!
        ! The actual formula is:  Si = e/ei = (e/esat)*(esat/ei) = S*(esat/ei)     !
        ! It is assumed that any supersaturation forms liquid water and that !
        ! the atmosphere is then saturated with respect to liquid water.     !
        ! Therefore, S = 1.0, allowing Si = esat/ei.                           !
        !--------------------------------------------------------------------!
        S_i(k) = e_s(k)/e_i(k)

        ! Denominator of diffusional growth equation, 9.4 of Rogers and Yau
        Denom(k) = Diff_denom( T_in_K(k), p_in_Pa(k), e_i(k) )

        ! Change in mass of a single ice crystal, m,
        ! as it falls a distance gr%invrs_dzt(1,:) in meters

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
        rcm_icedfsn(k) = - (N_i/rho(k)) & 
           * ( 4._core_rknd * (S_i(k) - 1._core_rknd) / Denom(k) ) & 
           * (mass_ice_cryst(k)/a_coef)**(1._core_rknd/b_expn)

        ! Ensure that liquid is not over-depleted
        IF ( rcm(k) + rcm_icedfsn(k)*dt < 0.0_core_rknd ) THEN
          rcm_icedfsn(k) = -rcm(k)/dt
        end if

        !---------------Brian's comment-----------------------------------!
        ! dm = (dm/dt)*(dt/dz)*dz                                         !
        ! dm = (dm/dt)*(1/u_T)*dz                                         !
        !-----------------------------------------------------------------!
        dmass_ice_cryst(k) = ( 4._core_rknd * (S_i(k) - 1._core_rknd) / Denom(k) ) & 
           * (1.0_core_rknd/k_u_coef) * ( rho(k)**q_expn ) & 
           * ( (mass_ice_cryst(k)/a_coef)**((1.0_core_rknd-n_expn)/b_expn) ) & 
           * (1.0_core_rknd/gr%invrs_dzm(1,k-1))
        mass_ice_cryst(k-1) = mass_ice_cryst(k)  & 
                                     + dmass_ice_cryst(k)

        ! Diameter of ice crystal in meters.
        diam(k) = (mass_ice_cryst(k)/a_coef)**(1._core_rknd/b_expn)

        ! Fallspeed of ice crystal in cm/s.
        u_T_cm(k) = cm_per_m * k_u_coef * & 
                    ((mass_ice_cryst(k)/a_coef)**(n_expn/b_expn))  & 
                          * (rho(k)**(-q_expn))

      else   ! There's no liquid and/or ice present; assume no ice growth

        mass_ice_cryst(k-1) = mass_ice_cryst(k)
        rcm_icedfsn(k) = 0.0_core_rknd
        diam(k)        = 0.0_core_rknd  ! Set zero to remind that we don't grow ice
        u_T_cm(k)      = 0.0_core_rknd  ! Set zero to remind that we don't grow ice

      end if

    end do ! k = gr%nz, 2, -1

    ! Michael Falk added boundary condx, 31 July 2006

    mass_ice_cryst(1) = mass_ice_cryst(2)
    rcm_icedfsn(1)    = rcm_icedfsn(2)
    diam(1)           = diam(2)
    u_T_cm(1)         = u_T_cm(2)

    ! eMFc

!
    if ( stats_metadata%l_stats_samp ) then

!       diam(:) ! Icedfs diameter; Michael Falk, 1 Nov 2006
!       m(:)    ! Icedfs mass; Michael Falk, 1 Nov 2006
!       dqc_dt_icedfs(:) ! Icedfs change in liquid; Michael Falk, 1 Nov 2006
!       u_T_cm(:)        ! Icedfs fallspeed (cm/s); Michael Falk, 1 Nov 2006
      call stat_update_var( stats_metadata%ircm_icedfs, rcm_icedfsn, stats_zt )

      call stat_update_var( stats_metadata%idiam, diam, stats_zt )

      call stat_update_var( stats_metadata%imass_ice_cryst, mass_ice_cryst, stats_zt )

      call stat_update_var( stats_metadata%iu_T_cm, u_T_cm, stats_zt )

    end if

    ! Determine time tendency of liquid potential temperature
    thlm_icedfsn(1:gr%nz) = - ( Lv/(Cp*exner(1:gr%nz)) ) * rcm_icedfsn(1:gr%nz)

    return
  end subroutine ice_dfsn

  !-----------------------------------------------------------------------------
  function Diff_denom( T_in_K, p_in_Pa, e_i )

    use constants_clubb, only: & 
        Ls,  & ! Constant(s)
        Rv, &
        T_freeze_K

    USE clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Description:
    !   Compute denominator of diffusional growth equation
    !
    ! References:
    !   Eqn. 9.4 of Rogers and Yau (1989), "A Short Course on Cloud Physics"
    !
    !-----------------------------------------------------------------------------

! Reference:  Eqn. 9.4 of Rogers and Yau (1989), "A Short Course on Cloud Physics"
    ! Constant Parameters
!   real, parameter :: Ls = 2.834e6

    real(kind=core_rknd), intent(in) ::  & 
     T_in_K,       & ! Temperature                               [K]
     p_in_Pa,        & ! Air pressure                              [Pa]
     e_i          ! Vapor pressure over ice                   [Pa]

    real(kind=core_rknd) ::  & 
      Diff_denom, & ! Denominator of diffusional growth equation  [m s kg^{-1}]
      Ka, Dv, &
      Fk, Fd, &
      Celsius

    ! ---- Begin Code ----

    Celsius = T_in_K - T_freeze_K

    Ka = (5.69_core_rknd + 0.017_core_rknd*Celsius)*0.00001_core_rknd  ! Ka in cal./(cm.*sec.*C)
    Ka = 4.1868_core_rknd*100.0_core_rknd*Ka  ! Ka in J./(m.*sec.*K)

    Dv = 0.221_core_rknd * ( (T_in_K/T_freeze_K)**1.94_core_rknd ) * &
          (101325.0_core_rknd/p_in_Pa)
    ! Dv in (cm.^2)/sec.  ! .221 is correct.
    Dv = Dv/10000.0_core_rknd  ! Dv in (m.^2)/sec.

    Fk = ( Ls/(Rv*T_in_K) - 1.0_core_rknd ) * Ls / (Ka*T_in_K)
    Fd = (Rv*T_in_K) / (Dv*e_i)

    Diff_denom = Fk + Fd

    return
  end function Diff_denom

end module ice_dfsn_module
