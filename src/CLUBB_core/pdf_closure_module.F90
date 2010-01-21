! $Id$
module pdf_closure_module

  implicit none

  public :: pdf_closure

  private ! Set Default Scope

  contains 
!------------------------------------------------------------------------
  subroutine pdf_closure &
             ( p_in_Pa, exner, thv_ds, wm,       &
               wp2, wp3, sigma_sqd_w,            &
               Skw, rtm, rtp2,                   &
               wprtp, thlm, thlp2,               &
               wpthlp, rtpthlp, sclrm,           &
               wpsclrp, sclrp2, sclrprtp,        &
               sclrpthlp, level,                 &
               wp4, wprtp2, wp2rtp,              &
               wpthlp2, wp2thlp, wprtpthlp,      &
               cloud_frac, rcm, wpthvp,          &
               wp2thvp, rtpthvp, thlpthvp,       &
               wprcp, wp2rcp, rtprcp,            &
               thlprcp, rcp2, pdf_params,        &
               err_code,                         &
               wpsclrprtp, wpsclrp2, sclrpthvp,  &
               wpsclrpthlp, sclrprcp, wp2sclrp )

!       Description:
!       Subroutine that computes pdf parameters analytically.

!       Based of the original formulation, but with some tweaks
!       to remove some of the less realistic assumptions and
!       improve transport terms.

!       Corrected version that should remove inconsistency

!       References:
!       Eqn. 29, 30, 31, 32 & 33  on p. 3547 of 
!       ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!         Method and Model Description'' Golaz, et al. (2002)
!       JAS, Vol. 59, pp. 3540--3551.
!------------------------------------------------------------------------

    use constants, only: & 
      ! Constants
      sqrt_2pi,      & ! sqrt(2*pi)
      sqrt_2,        & ! sqrt(2)
      pi,            & ! The ratio of radii to their circumference
      Cp,            & ! Dry air specific heat at constant p [J/kg/K]
      Lv,            & ! Latent heat of vaporization         [J/kg]
      Rd,            & ! Dry air gas constant                [J/kg/K]
      Rv,            & ! Water vapor gas constant            [J/kg/K]
      ep,            & ! Rd / Rv;     ep  = 0.622            [-]
      ep1,           & ! (1.0-ep)/ep; ep1 = 0.61             [-]
      ep2,           & ! 1.0/ep;      ep2 = 1.61             [-]
      wtol_sqd,      & ! Tolerance for w'^2                  [m^2/s^2]
      rttol,         & ! Tolerance for r_t                   [kg/kg]
      thltol,        & ! Tolerance for th_l                  [K]
      s_mellor_tol,  & ! Tolerance for pdf parameter s       [kg/kg]
      fstderr,       &
      zero_threshold  

    use parameters_model, only: &
      sclrtol,   & ! Array of passive scalar tolerances  [units vary]
      sclr_dim,  & ! Number of passive scalar variables
      a_max_mag    ! Maximum values for PDF parameter 'a'

    use parameters_tunable, only: & 
      beta  ! Variable(s)
            ! Plume widths for th_l and r_t [-]

    use variables_prognostic_module, only:  &
        pdf_parameter  ! type

    use anl_erf, only:  & 
      erf ! Procedure(s)
          ! The error function

    use numerical_check, only:  & 
      pdf_closure_check ! Procedure(s)

    use saturation, only:  & 
      sat_mixrat_liq ! Procedure(s)

    use error_code, only:  & 
      clubb_var_equals_NaN,  & ! Variable(s)
      clubb_at_least_debug_level ! Procedure(s)

    use stats_variables, only: &
      iwp4,       & ! Variables
      ircp2,      &
      iwprtp2,    &
      iwprtpthlp, &
      iwpthlp2

    implicit none

    intrinsic :: sqrt, exp, min, max, abs

    ! Input Variables
    real, intent(in) ::  & 
      p_in_Pa,     & ! Pressure                                   [Pa]
      exner,       & ! Exner function                             [-]
      thv_ds,      & ! Dry, base-state theta_v (ref. th_l here)   [K]
      wm,          & ! mean w-wind component (vertical velocity)  [m/s] 
      wp2,         & ! w'^2                                       [m^2/s^2] 
      wp3,         & ! w'^3                                       [m^3/s^3]
      sigma_sqd_w, & ! Width of individual w plumes               [-]
      Skw,         & ! Skewness of w                              [-]
      rtm,         & ! Mean total water mixing ratio              [kg/kg]
      rtp2,        & ! r_t'^2                                     [(kg/kg)^2]
      wprtp,       & ! w'r_t'                                     [(kg/kg)(m/s)]
      thlm,        & ! Mean liquid water potential temperature    [K]
      thlp2,       & ! th_l'^2                                    [K^2]
      wpthlp,      & ! w'th_l'                                    [K(m/s)]
      rtpthlp        ! r_t'th_l'                                  [K(kg/kg)]

    real, dimension(sclr_dim), intent(in) ::  & 
      sclrm,       & ! Mean passive scalar        [units vary]
      wpsclrp,     & ! w' sclr'                   [units vary]
      sclrp2,      & ! sclr'^2                    [units vary]
      sclrprtp,    & ! sclr' r_t'                 [units vary]
      sclrpthlp      ! sclr' th_l'                [units vary]

    integer, intent(in) ::  &
      level  ! Thermodynamic level for which calculations are taking place.

    ! Output Variables

    real, intent(out) ::  & 
      wp4,         & ! w'^4                  [m^4/s^4]
      wprtp2,      & ! w' r_t'               [(m kg)/(s kg)]
      wp2rtp,      & ! w'^2 r_t'             [(m^2 kg)/(s^2 kg)]
      wpthlp2,     & ! w' th_l'^2            [(m K^2)/s]
      wp2thlp,     & ! w'^2 th_l'            [(m^2 K)/s^2]
      cloud_frac,  & ! Cloud fraction        [-]
      rcm,         & ! Mean liquid water     [kg/kg]
      wpthvp,      & ! Buoyancy flux         [(K m)/s] 
      wp2thvp,     & ! w'^2 th_v'            [(m^2 K)/s^2]
      rtpthvp,     & ! r_t' th_v'            [(kg K)/kg]
      thlpthvp,    & ! th_l' th_v'           [K^2]
      wprcp,       & ! w' r_c'               [(m kg)/(s kg)]
      wp2rcp,      & ! w'^2 r_c'             [(m^2 kg)/(s^2 kg)]
      rtprcp,      & ! r_t' r_c'             [(kg^2)/(kg^2)]
      thlprcp,     & ! th_l' r_c'            [(K kg)/kg]
      rcp2,        & ! r_c'^2                [(kg^2)/(kg^2)]
      wprtpthlp      ! w' r_t' th_l'         [(m kg K)/(s kg)]

    type(pdf_parameter), intent(out) :: & 
      pdf_params     ! pdf paramters         [units vary]

    integer, intent(out) :: & 
      err_code       ! Are the outputs usable numbers?

    ! Output (passive scalar variables)

    real, intent(out), dimension(sclr_dim) ::  & 
      sclrpthvp, & 
      sclrprcp, & 
      wpsclrp2, & 
      wpsclrprtp, & 
      wpsclrpthlp, & 
      wp2sclrp

    ! Local Variables 

    real ::  & 
      w1_n, w2_n 
!     thl1_n, thl2_n, 
!     rt1_n, rt2_n

    ! Variables that are stored in derived data type pdf_params.
    real ::  &
      w1,          & ! Mean of w for 1st normal distribution                 [m/s]
      w2,          & ! Mean of w for 2nd normal distribution                 [m/s]
      varnce_w1,   & ! Variance of w for 1st normal distribution         [m^2/s^2]
      varnce_w2,   & ! Variance of w for 2nd normal distribution         [m^2/s^2]
      rt1,         & ! Mean of r_t for 1st normal distribution             [kg/kg]
      rt2,         & ! Mean of r_t for 2nd normal distribution             [kg/kg]
      varnce_rt1,  & ! Variance of r_t for 1st normal distribution     [kg^2/kg^2]
      varnce_rt2,  & ! Variance of r_t for 2nd normal distribution     [kg^2/kg^2]
      crt1,        & ! Coefficient for s'                                      [-]
      crt2,        & ! Coefficient for s'                                      [-]
      cthl1,       & ! Coefficient for s'                                    [1/K]
      cthl2,       & ! Coefficient for s'                                    [1/K]
      thl1,        & ! Mean of th_l for 1st normal distribution                [K]
      thl2,        & ! Mean of th_l for 2nd normal distribution                [K]
      varnce_thl1, & ! Variance of th_l for 1st normal distribution          [K^2]
      varnce_thl2, & ! Variance of th_l for 2nd normal distribution          [K^2]
      a,           & ! Weight of 1st normal distribution (Sk_w dependent)      [-]
      rc1,         & ! Mean of r_c for 1st normal distribution             [kg/kg]
      rc2,         & ! Mean of r_c for 2nd normal distribution             [kg/kg]
      rsl1,        & ! Mean of r_sl for 1st normal distribution            [kg/kg]
      rsl2,        & ! Mean of r_sl for 2nd normal distribution            [kg/kg]
      cloud_frac1, & ! Cloud fraction for 1st normal distribution              [-]
      cloud_frac2, & ! Cloud fraction for 2nd normal distribution              [-]
      s1,          & ! Mean of s for 1st normal distribution               [kg/kg]
      s2,          & ! Mean of s for 2nd normal distribution               [kg/kg]
      stdev_s1,    & ! Standard deviation of s for 1st normal distribution [kg/kg]
      stdev_s2,    & ! Standard deviation of s for 2nd normal distribution [kg/kg]
      rrtthl,      & ! Within-a-normal (sub-plume) correlation of r_t and th_l [-]
      alpha_thl,   & ! Factor relating to normalized variance for th_l         [-]
      alpha_rt       ! Factor relating to normalized variance for r_t          [-]

                   ! Note:  alpha coefficients = 0.5 * ( 1 - correlations^2 ).
                   !        These are used to calculate the scalar widths
                   !        varnce_thl1, varnce_thl2, varnce_rt1, and varnce_rt2 as in Eq. (34) of
                   !        Larson and Golaz (2005)

    ! Passive scalar local variables

    real, dimension(sclr_dim) ::  & 
      sclr1, sclr2,  &
      varnce_sclr1, varnce_sclr2, & 
      alpha_sclr,  & 
      rsclrthl, rsclrrt
!     sclr1_n, sclr2_n,
 
    logical :: l_scalar_calc ! True if sclr_dim > 0

    ! Quantities needed to predict higher order moments
    real ::  & 
      tl1, tl2,  & 
      beta1, beta2,  & 
      zeta1, zeta2

    real :: sqrt_wp2

    ! Thermodynamic quantity

    real :: rc_coef

    ! variables for a generalization of Chris Golaz' closure
    ! varies width of plumes in theta_l, rt
    real :: width_factor_1, width_factor_2

    integer :: i   ! Index

!------------------------ Code Begins ----------------------------------

    ! Check whether the passive scalars are present.

    if ( sclr_dim > 0 ) then
      l_scalar_calc = .true.
    else
      l_scalar_calc = .false.
    end if

    ! If there is no velocity, then use single delta fnc. as pdf
    ! Otherwise width parameters (e.g. varnce_w1, varnce_w2, etc.) are non-zero.
    if ( wp2 <= wtol_sqd )  then

      a         = 0.5
      w1        = wm
      w2        = wm
      varnce_w1 = 0.
      varnce_w2 = 0.
      rt1       = rtm
      rt2       = rtm
      alpha_rt  = 0.5
      varnce_rt1      = 0.
      varnce_rt2      = 0.
      thl1      = thlm
      thl2      = thlm
      alpha_thl = 0.5
      varnce_thl1     = 0.
      varnce_thl2     = 0.
      rrtthl    = 0.

      if ( l_scalar_calc ) then
        do i = 1, sclr_dim, 1
          sclr1(i)        = sclrm(i)
          sclr2(i)        = sclrm(i)
          varnce_sclr1(i) = 0.0
          varnce_sclr2(i) = 0.0
          alpha_sclr(i)   = 0.5
          rsclrrt(i)      = 0.0
          rsclrthl(i)     = 0.0
        end do ! 1..sclr_dim
      end if

    else ! Width (standard deviation) parameters are non-zero

      ! The variable "a" is the weight of Gaussian "plume" 1.  The weight of
      ! Gaussian "plume" 2 is "1-a".  If there isn't any skewness of w
      ! (Sk_w = 0 because w'^3 = 0), a = 0.5, and both Gaussian "plumes" are
      ! equally weighted.  If there is positive skewness of w (Sk_w > 0 because
      ! w'^3 > 0), 0 < a < 0.5, and Gaussian "plume" 2 has greater weight than
      ! does Gaussian "plume" 1.  If there is negative skewness of w (Sk_w < 0
      ! because w'^3 < 0), 0.5 < a < 1, and Gaussian "plume" 1 has greater
      ! weight than does Gaussian "plume" 2.
      if ( abs( Skw ) <= 1e-5 ) then
        a = 0.5
      else
        a = 0.5 * ( 1.0 - Skw/sqrt( 4.0*( 1.0 - sigma_sqd_w )**3 + Skw**2 ) )
      endif

      ! Determine sqrt( wp2 ) here to avoid re-computing it
      sqrt_wp2 = sqrt( wp2 )

      ! Clip a, 1-a, to avoid dividing by zero
      ! Formula for a_max_mag =
      ! 1 - ( 1/2 * ( 1 - Skw_max/sqrt( 4*( 1 - sigma_sqd_w )^3 + Skw_max^2 ) ) )
      ! Where sigma_sqd_w is fixed at 0.4
      a = min( max( a, 1.0-a_max_mag ), a_max_mag )

      ! The normalized mean of w for Gaussian "plume" 1 is w1_n.  It's value
      ! will always be greater than 0.  As an example, a value of 1.0 would
      ! indicate that the actual mean of w for Gaussian "plume" 1 is found
      ! 1.0 standard deviation above the overall mean for w.
      w1_n = sqrt( ( (1.-a)/a )*(1.-sigma_sqd_w) )
      ! The normalized mean of w for Gaussian "plume" 2 is w2_n.  It's value
      ! will always be less than 0.  As an example, a value of -0.5 would
      ! indicate that the actual mean of w for Gaussian "plume" 2 is found
      ! 0.5 standard deviations below the overall mean for w.
      w2_n = -sqrt( ( a/(1.-a) )*(1.-sigma_sqd_w) )
      ! The mean of w for Gaussian "plume" 1 is w1.
      w1   = wm + sqrt_wp2*w1_n
      ! The mean of w for Gaussian "plume" 2 is w2.
      w2   = wm + sqrt_wp2*w2_n

      ! The variance of w for Gaussian "plume" 1 for varnce_w1.
      varnce_w1  = sigma_sqd_w*wp2
      ! The variance of w for Gaussian "plume" 2 for varnce_w2.
      ! The variance in both Gaussian "plumes" is defined to be the same.
      varnce_w2  = sigma_sqd_w*wp2


      ! The normalized variance for thl, rt, and sclr for "plume" 1 is:
      !
      ! { 1 - [1/(1-sigma_sqd_w)]*[ (w'x')^2 / (w'^2 * x'^2) ] / a }
      ! * { (1/3)*beta + a*( 1 - (2/3)*beta ) };
      !
      ! where "x" stands for thl, rt, or sclr; "a" is the weight of Gaussian
      ! "plume" 1, and 0 <= beta <= 3.
      !
      ! The factor { (1/3)*beta + a*( 1 - (2/3)*beta ) } does not depend on
      ! which varable "x" stands for.  The factor is multiplied by 2 and defined
      ! as width_factor_1.
      !
      ! The factor { 1 - [1/(1-sigma_sqd_w)]*[ (w'x')^2 / (w'^2 * x'^2) ] / a }
      ! depends on which variable "x" stands for.  It is multiplied by 0.5 and
      ! defined as alpha_x, where "x" stands for thl, rt, or sclr.

      ! Vince Larson added a dimensionless factor so that the
      ! width of plumes in theta_l, rt can vary.  
      ! beta is a constant defined in module parameters_tunable  
      !   Set 0<beta<3.
      ! beta=1.5 recovers Chris Golaz' simplified formula.
      ! 3 Nov 2003

      width_factor_1 = ( 2.0/3.0 )*beta + 2.0*a*( 1.0 - ( 2.0/3.0 )*beta )
      width_factor_2 = 2.0 - width_factor_1

      if ( thlp2 <= thltol**2 ) then
        thl1      = thlm
        thl2      = thlm
        varnce_thl1     = 0.0
        varnce_thl2     = 0.0
        alpha_thl = 0.5
      else
!       thl1_n = - (wpthlp/(sqrt( wp2 )*sqrt( thlp2 )))/w2_n
!       thl2_n = - (wpthlp/(sqrt( wp2 )*sqrt( thlp2 )))/w1_n

        thl1 = thlm - ( wpthlp/sqrt_wp2 )/w2_n
        thl2 = thlm - ( wpthlp/sqrt_wp2 )/w1_n

        alpha_thl = 0.5 * ( 1.0 - wpthlp*wpthlp /((1.0-sigma_sqd_w)*wp2*thlp2) )

        alpha_thl = max( min( alpha_thl, 1.0 ), zero_threshold )

        ! Vince Larson multiplied original expressions by width_factor_1,2
        !   to generalize scalar skewnesses.  05 Nov 03
        varnce_thl1 = ( alpha_thl / a * thlp2 ) * width_factor_1
        varnce_thl2 = ( alpha_thl / (1.-a) * thlp2 ) * width_factor_2

      end if ! thlp2 <= thltol**2

      if ( rtp2 <= rttol**2 ) then
        rt1      = rtm
        rt2      = rtm
        varnce_rt1     = 0.0
        varnce_rt2     = 0.0
        alpha_rt = 0.5
      else
!       rt1_n = -( wprtp / ( sqrt( wp2 )*sqrt( rtp2 ) ) ) / w2_n
!       rt2_n = -( wprtp / ( sqrt( wp2 )*sqrt( rtp2 ) ) ) / w1_n

        rt1 = rtm - ( wprtp / sqrt_wp2 ) / w2_n
        rt2 = rtm - ( wprtp / sqrt_wp2 ) / w1_n

        alpha_rt = 0.5 * ( 1.0 - wprtp*wprtp /((1.0-sigma_sqd_w)*wp2*rtp2) )

        alpha_rt = max( min( alpha_rt, 1.0 ), zero_threshold )

      ! Vince Larson multiplied original expressions by width_factor_1,2
      !   to generalize scalar skewnesses.  05 Nov 03
        varnce_rt1 = ( alpha_rt / a * rtp2 ) * width_factor_1
        varnce_rt2 = ( alpha_rt / (1.-a) * rtp2 ) * width_factor_2

      end if ! rtp2 <= rttol**2 

      ! Compute pdf parameters for passive scalars
      if ( l_scalar_calc ) then
        do i = 1, sclr_dim
          if ( sclrp2(i) <= sclrtol(i)**2 ) then
            ! Set plume sclr for plume 1,2 to the mean
            sclr1(i)= sclrm(i)
            sclr2(i)= sclrm(i)
            ! Set the variance to zero
            varnce_sclr1(i) = 0.0
            varnce_sclr2(i) = 0.0

            alpha_sclr(i) = 0.5
          else
!           sclr1_n(i) = - ( wpsclrp(i) / (sqrt( wp2 ) &
!                        * sqrt( sclrp2(i) )) )/w2_n
!           sclr2_n(i) = - ( wpsclrp(i) / (sqrt( wp2 ) &
!                        * sqrt( sclrp2(i) )) )/w1_n
 
            sclr1(i) = sclrm(i)  & 
                     - ( wpsclrp(i) / sqrt_wp2 ) / w2_n
            sclr2(i) = sclrm(i)  & 
                     - ( wpsclrp(i) / sqrt_wp2 ) / w1_n

            alpha_sclr(i) = 0.5 * ( 1.0 - wpsclrp(i)*wpsclrp(i) & 
                    / ((1.0-sigma_sqd_w)*wp2*sclrp2(i)) )

            alpha_sclr(i) = max( min( alpha_sclr(i), 1.0 ), zero_threshold )

            ! Vince Larson multiplied original expressions by width_factor_1,2
            !  to generalize scalar skewnesses.  05 Nov 03
            varnce_sclr1(i) = ( alpha_sclr(i) / a * sclrp2(i) ) * width_factor_1
            varnce_sclr2(i) = ( alpha_sclr(i) / (1.-a) * sclrp2(i) ) * width_factor_2
          end if ! sclrp2(i) <= sclrtol(i)**2
        end do ! i=1, sclr_dim
      end if ! l_scalar_calc

      ! We include sub-plume correlation with coeff rrtthl.

      if ( varnce_rt1*varnce_thl1 > 0 .and. varnce_rt2*varnce_thl2 > 0) then
        rrtthl = ( rtpthlp - a * ( rt1-rtm ) * ( thl1-thlm ) & 
                   - (1.-a) * ( rt2-rtm ) * ( thl2-thlm ) ) & 
                / ( a*sqrt( varnce_rt1*varnce_thl1 ) + (1.-a)*sqrt( varnce_rt2*varnce_thl2 ) )
        if ( rrtthl < -1.0 ) then
          rrtthl = -1.0
        end if
        if ( rrtthl > 1.0 ) then
          rrtthl = 1.0
        end if
      else
        rrtthl = 0.0
      end if ! varnce_rt1*varnce_thl1 > 0 .and. varnce_rt2*varnce_thl2 > 0

      ! Sub-plume correlation, rsclrthl, between passive scalar and theta_l.
      if ( l_scalar_calc ) then
        do i=1, sclr_dim
          if ( varnce_sclr1(i)*varnce_thl1 > 0. .and. varnce_sclr2(i)*varnce_thl2 > 0. ) then
            rsclrthl(i) = ( sclrpthlp(i)  & 
            - a * ( sclr1(i)-sclrm(i) ) * ( thl1-thlm ) & 
            - (1.-a) * ( sclr2(i)-sclrm(i) ) * ( thl2-thlm ) ) & 
                / ( a*sqrt( varnce_sclr1(i)*varnce_thl1 )  & 
                         + (1.-a)*sqrt( varnce_sclr2(i)*varnce_thl2 ) )
            if ( rsclrthl(i) < -1.0 ) then
              rsclrthl(i) = -1.0
            end if
            if ( rsclrthl(i) > 1.0 ) then
              rsclrthl(i) = 1.0
            end if
          else
            rsclrthl(i) = 0.0
          end if

          ! Sub-plume correlation, rsclrrt, between passive scalar
          !   and total water.

          if ( varnce_sclr1(i)*varnce_rt1 > 0 .and. varnce_sclr2(i)*varnce_rt2 > 0 ) then
            rsclrrt(i) = ( sclrprtp(i) - a * ( sclr1(i)-sclrm(i) ) * ( rt1-rtm )&
                         - (1.-a) * ( sclr2(i)-sclrm(i) ) * ( rt2-rtm ) ) & 
             / ( a*sqrt( varnce_sclr1(i)*varnce_rt1 ) + (1.-a)*sqrt( varnce_sclr2(i)*varnce_rt2 ) )
            if ( rsclrrt(i) < -1.0 ) then
              rsclrrt(i) = -1.0
            end if
            if ( rsclrrt(i) > 1.0 ) then
              rsclrrt(i) = 1.0
            end if
          else
            rsclrrt(i) = 0.0
          end if
        end do ! i=1, sclr_dim
      end if ! l_scalar_calc 

    end if  ! Widths non-zero

    ! Compute higher order moments (these are interactive)
    wp2rtp  = a * ( (w1-wm)**2+varnce_w1 ) * ( rt1-rtm ) & 
            + (1.-a) * ( (w2-wm)**2+varnce_w2 ) * ( rt2-rtm )

    wp2thlp = a * ( (w1-wm)**2+varnce_w1 ) * ( thl1-thlm ) & 
            + (1.-a) * ( (w2-wm)**2+varnce_w2 ) * ( thl2-thlm )

    ! Compute higher order moments (these are non-interactive diagnostics)
    if ( iwp4 > 0 ) then
      wp4 = a * ( 3.*varnce_w1**2 + 6.*((w1-wm)**2)*varnce_w1 + (w1-wm)**4 ) & 
          + (1.-a) * ( 3.*varnce_w2**2 + 6.*((w2-wm)**2)*varnce_w2 + (w2-wm)**4 )
    end if

    if ( iwprtp2 > 0 ) then
      wprtp2  = a * ( w1-wm )*( (rt1-rtm)**2 + varnce_rt1 )  & 
              + (1.-a) * ( w2-wm )*( (rt2-rtm)**2 + varnce_rt2)
    end if

    if ( iwpthlp2 > 0 ) then
      wpthlp2 = a * ( w1-wm )*( (thl1-thlm)**2 + varnce_thl1 )  & 
              + (1.-a) * ( w2-wm )*( (thl2-thlm)**2+varnce_thl2 )
    end if

    if ( iwprtpthlp > 0 ) then
      wprtpthlp = a * ( w1-wm )*( (rt1-rtm)*(thl1-thlm)  & 
                + rrtthl*sqrt( varnce_rt1*varnce_thl1 ) ) & 
                + ( 1.-a ) * ( w2-wm )*( (rt2-rtm)*(thl2-thlm) & 
                + rrtthl*sqrt( varnce_rt2*varnce_thl2 ) )
    end if

    ! Scalar Addition to higher order moments
    if ( l_scalar_calc ) then
      do i=1, sclr_dim

        wp2sclrp(i)  = a * ( (w1-wm)**2+varnce_w1 )*( sclr1(i)-sclrm(i) ) & 
                     + (1.-a) * ( (w2-wm)**2+varnce_w2 ) * ( sclr2(i)-sclrm(i) )

        wpsclrp2(i) = a * ( w1-wm ) * ( (sclr1(i)-sclrm(i))**2 + varnce_sclr1(i) )  & 
                    + (1.-a) * ( w2-wm ) * ( (sclr2(i)-sclrm(i))**2 + varnce_sclr2(i) )

        wpsclrprtp(i) = a * ( w1-wm ) * ( ( rt1-rtm )*( sclr1(i)-sclrm(i) )  & 
          + rsclrrt(i)*sqrt( varnce_rt1*varnce_sclr1(i) ) ) &
          + ( 1.-a )*( w2-wm ) *  &
            ( ( rt2-rtm )*( sclr2(i)-sclrm(i) ) + rsclrrt(i)*sqrt( varnce_rt2*varnce_sclr2(i) ) )

        wpsclrpthlp(i) = a * ( w1-wm ) * ( ( sclr1(i)-sclrm(i) )*( thl1-thlm )  & 
          + rsclrthl(i)*sqrt( varnce_sclr1(i)*varnce_thl1 ) ) & 
          + ( 1.-a ) * ( w2-wm ) * &
            ( ( sclr2(i)-sclrm(i) )*( thl2-thlm ) &
              + rsclrthl(i)*sqrt( varnce_sclr2(i)*varnce_thl2 ) )

      end do ! i=1, sclr_dim
    end if ! l_scalar_calc

    ! Compute higher order moments that include theta_v.

    ! First compute some preliminary quantities.
    ! "1" denotes first Gaussian; "2" denotes 2nd Gaussian
    ! liq water temp (Sommeria & Deardorff 1977 (SD), eqn. 3)

    tl1  = thl1*exner
    tl2  = thl2*exner

    rsl1 = sat_mixrat_liq( p_in_Pa, tl1 )
    rsl2 = sat_mixrat_liq( p_in_Pa, tl2 )

    ! SD's beta (eqn. 8)
    beta1 = ep * ( Lv/(Rd*tl1) ) * ( Lv/(Cp*tl1) )
    beta2 = ep * ( Lv/(Rd*tl2) ) * ( Lv/(Cp*tl2) )

    ! s from Lewellen and Yoh 1993 (LY) eqn. 1
    s1 = ( rt1 - rsl1 ) / ( 1. + beta1 * rsl1 )
    s2 = ( rt2 - rsl2 ) / ( 1. + beta2 * rsl2 )

    ! Coefficients for s'
    ! For each normal distribution in the sum of two normal distributions,
    ! s' = crt * rt'  +  cthl * thl';
    ! therefore, x's' = crt * x'rt'  +  cthl * x'thl'.
    ! Larson et al. May, 2001.

    crt1  = 1./( 1. + beta1*rsl1)
    crt2  = 1./( 1. + beta2*rsl2)

    cthl1 = ( (1. + beta1 * rt1) / ( 1. + beta1*rsl1)**2 ) & 
             * ( Cp/Lv ) * beta1 * rsl1 * exner
    cthl2 = ( (1. + beta2 * rt2) / ( 1. + beta2*rsl2 )**2 ) & 
             * ( Cp/Lv ) * beta2 * rsl2 * exner

    ! Standard deviation of s
    ! include subplume correlation of qt, thl
    ! Because of round-off error, 
    ! stdev_s1 (and probably stdev_s2) can become negative when rrtthl=1
    ! One could also write this as a squared term
    ! plus a postive correction; this might be a neater format

    stdev_s1 = sqrt( max( zero_threshold, ( varnce_rt1*crt1**2 + varnce_thl1*cthl1**2  &
        - 2.0*rrtthl*crt1*sqrt( varnce_rt1*varnce_thl1 )*cthl1 )  & 
               ) &  ! max
          ) ! sqrt
    stdev_s2 = sqrt( max( zero_threshold, ( varnce_rt2*crt2**2 + varnce_thl2*cthl2**2 & 
        - 2.0*rrtthl*crt2*sqrt( varnce_rt2*varnce_thl2 )*cthl2 )  & 
               )  &  ! max
          ) ! sqrt

!   stdev_s1 = sqrt( (sqrt(varnce_rt1)*crt1 - sqrt(varnce_thl1)*cthl1)**2 &
!                + (1.-rrtthl)*2.*crt1*sqrt(varnce_rt1)*cthl1*sqrt(varnce_thl1)  )
!   stdev_s2 = sqrt( (sqrt(varnce_rt2)*crt2 - sqrt(varnce_thl2)*cthl2)**2 &
!                + (1.-rrtthl)*2.*crt2*sqrt(varnce_rt2)*cthl2*sqrt(varnce_thl2)  )


    ! We need to introduce a threshold value for the variance of s

    if ( stdev_s1 > s_mellor_tol ) then
      zeta1 = s1/stdev_s1
      cloud_frac1  = 0.5*( 1. + erf( zeta1/sqrt_2 )  )
      rc1          = s1*cloud_frac1+stdev_s1*exp( -0.5*zeta1**2 )/( sqrt_2pi )
    else
      if ( s1 < 0.0 ) then
        cloud_frac1  = 0.0
        rc1          = 0.0
      else
        cloud_frac1  = 1.0
        rc1          = s1
      end if ! s1 < 0
    end if ! stdev_s1 > s_mellor_tol

    if ( stdev_s2 > s_mellor_tol ) then
      zeta2       = s2/stdev_s2
      cloud_frac2 = 0.5*( 1. + erf( zeta2/sqrt_2 ) )
      rc2         = s2*cloud_frac2+stdev_s2*exp( -0.5*zeta2**2 )/( sqrt_2pi )
    else
      if ( s2 < 0.0 ) then
        cloud_frac2  = 0.0
        rc2          = 0.0
      else
        cloud_frac2  = 1.0
        rc2          = s2
      end if ! s2 < 0
    end if ! stdev_s2 > s_mellor_tol

    ! Compute moments that depend on theta_v
    !
    ! The moments that depend on th_v' are calculated based on an approximated
    ! and linearized form of the theta_v equation:
    !
    ! theta_v = theta_l + { (R_v/R_d) - 1 } * thv_ds * r_t
    !                   + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c;
    !
    ! and therefore:
    !
    ! th_v' = th_l' + { (R_v/R_d) - 1 } * thv_ds * r_t'
    !               + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c';
    !
    ! where thv_ds is used as a reference value to approximate theta_l.

    rc_coef = Lv / (exner*Cp) - ep2 * thv_ds

    wp2rcp = a * ((w1-wm)**2 + varnce_w1)*rc1 + (1.-a) * ((w2-wm)**2 + varnce_w2)*rc2 & 
           - wp2 * (a*rc1+(1.-a)*rc2)

    wp2thvp = wp2thlp + ep1*thv_ds*wp2rtp + rc_coef*wp2rcp

    wprcp = a * (w1-wm)*rc1 + (1.-a) * (w2-wm)*rc2

    wpthvp = wpthlp + ep1*thv_ds*wprtp + rc_coef*wprcp

    ! Account for subplume correlation in qt-thl
    thlprcp  = a * ( (thl1-thlm)*rc1 - (cthl1*varnce_thl1)*cloud_frac1 ) & 
             + (1.-a) * ( (thl2-thlm)*rc2 - (cthl2*varnce_thl2)*cloud_frac2 ) & 
             + a*rrtthl*crt1*sqrt( varnce_rt1*varnce_thl1 )*cloud_frac1 & 
             + (1.-a)*rrtthl*crt2*sqrt( varnce_rt2*varnce_thl2 )*cloud_frac2
    thlpthvp = thlp2 + ep1*thv_ds*rtpthlp + rc_coef*thlprcp

    ! Account for subplume correlation in qt-thl
    rtprcp = a * ( (rt1-rtm)*rc1 + (crt1*varnce_rt1)*cloud_frac1 ) & 
           + (1.-a) * ( (rt2-rtm)*rc2 + (crt2*varnce_rt2)*cloud_frac2 ) & 
           - a*rrtthl*cthl1*sqrt( varnce_rt1*varnce_thl1 )*cloud_frac1 & 
           - (1.-a)*rrtthl*cthl2*sqrt( varnce_rt2*varnce_thl2 )*cloud_frac2

    rtpthvp  = rtpthlp + ep1*thv_ds*rtp2 + rc_coef*rtprcp

    ! Account for subplume correlation between scalar, theta_v.
    ! See Eqs. A13, A8 from Larson et al. (2002) ``Small-scale...''
    !  where the ``scalar'' in this paper is w.
    if ( l_scalar_calc ) then
      do i=1, sclr_dim
        sclrprcp(i) &
        = a * ( ( sclr1(i)-sclrm(i) ) * rc1 ) + (1.-a) * ( ( sclr2(i)-sclrm(i) ) * rc2 ) & 
        + a*rsclrrt(i) * crt1  * sqrt( varnce_sclr1(i) * varnce_rt1 ) * cloud_frac1 & 
        + (1.-a) * rsclrrt(i) * crt2  * sqrt( varnce_sclr2(i) * varnce_rt2 ) * cloud_frac2 & 
        - a * rsclrthl(i) * cthl1  * sqrt( varnce_sclr1(i) * varnce_thl1 ) * cloud_frac1 & 
        - (1.-a) * rsclrthl(i) * cthl2  * sqrt( varnce_sclr2(i) * varnce_thl2 ) * cloud_frac2

        sclrpthvp(i) = sclrpthlp(i) + ep1*thv_ds*sclrprtp(i) + rc_coef*sclrprcp(i)
      end do ! i=1, sclr_dim
    end if ! l_scalar_calc

    ! Compute mean cloud fraction and cloud water

    cloud_frac = a * cloud_frac1 + (1.-a) * cloud_frac2
    rcm        = a * rc1         + (1.-a) * rc2

    ! Note: Brian added the following lines to ensure that there
    ! are never any negative liquid water values (or any negative
    ! cloud fraction values, for that matter).  According to
    ! Vince Larson, the analytic formula should not produce any
    ! negative results, but such computer-induced errors such as
    ! round-off error may produce such a value.  This has been
    ! corrected because Brian found a small negative value of
    ! rcm in the first timestep of the FIRE case.

    cloud_frac  = max( zero_threshold, cloud_frac )
    rcm = max( zero_threshold, rcm )

    ! Compute variance of liquid water mixing ratio.
    ! This is not needed for closure.  Statistical Analysis only.
    if ( ircp2 > 0 ) then

      rcp2 = a * ( s1*rc1 + cloud_frac1*stdev_s1**2 ) + ( 1.-a ) &
            * ( s2*rc2 + cloud_frac2*stdev_s2**2 ) - rcm**2
      rcp2 = max( zero_threshold, rcp2 )

    end if


    ! Save PDF parameters
    pdf_params%w1(level)          = w1
    pdf_params%w2(level)          = w2
    pdf_params%varnce_w1(level)   = varnce_w1
    pdf_params%varnce_w2(level)   = varnce_w2
    pdf_params%rt1(level)         = rt1
    pdf_params%rt2(level)         = rt2
    pdf_params%varnce_rt1(level)  = varnce_rt1
    pdf_params%varnce_rt2(level)        = varnce_rt2
    pdf_params%crt1(level)        = crt1
    pdf_params%crt2(level)        = crt2
    pdf_params%cthl1(level)       = cthl1
    pdf_params%cthl2(level)       = cthl2
    pdf_params%thl1(level)        = thl1
    pdf_params%thl2(level)        = thl2
    pdf_params%varnce_thl1(level) = varnce_thl1
    pdf_params%varnce_thl2(level) = varnce_thl2
    pdf_params%a(level)           = a
    pdf_params%rc1(level)         = rc1
    pdf_params%rc2(level)         = rc2
    pdf_params%rsl1(level)        = rsl1
    pdf_params%rsl2(level)        = rsl2
    pdf_params%cloud_frac1(level) = cloud_frac1
    pdf_params%cloud_frac2(level) = cloud_frac2
    pdf_params%s1(level)          = s1
    pdf_params%s2(level)          = s2
    pdf_params%stdev_s1(level)         = stdev_s1
    pdf_params%stdev_s2(level)         = stdev_s2
    pdf_params%rrtthl(level)      = rrtthl
    pdf_params%alpha_thl(level)   = alpha_thl
    pdf_params%alpha_rt(level)    = alpha_rt


    if ( clubb_at_least_debug_level( 2 ) ) then 
      call pdf_closure_check & 
           ( wp4, wprtp2, wp2rtp, wpthlp2, & 
             wp2thlp, cloud_frac, rcm, wpthvp, wp2thvp, & 
             rtpthvp, thlpthvp, wprcp, wp2rcp, & 
             rtprcp, thlprcp, rcp2, wprtpthlp, & 
             crt1, crt2, cthl1, cthl2, pdf_params, level, &
             err_code, & 
             sclrpthvp, sclrprcp, wpsclrp2, & 
             wpsclrprtp, wpsclrpthlp, wp2sclrp )

      ! Error Reporting
      ! Joshua Fasching February 2008
        
      if ( err_code == clubb_var_equals_NaN ) then
        
        write(fstderr,*) "Error in pdf_closure_new" 
          
        write(fstderr,*) "Intent(in)"
            
        write(fstderr,*) "p_in_Pa = ", p_in_Pa 
        write(fstderr,*) "exner = ", exner
        write(fstderr,*) "thv_ds = ", thv_ds
        write(fstderr,*) "wm = ", wm 
        write(fstderr,*) "wp2 = ", wp2 
        write(fstderr,*) "wp3 = ", wp3
        write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w
        write(fstderr,*) "rtm = ", rtm
        write(fstderr,*) "rtp2 = ", rtp2
        write(fstderr,*) "wprtp = ", wprtp
        write(fstderr,*) "thlm = ", thlm
        write(fstderr,*) "thlp2 = ", thlp2
        write(fstderr,*) "wpthlp = ", wpthlp
        write(fstderr,*) "rtpthlp = ", rtpthlp

        if ( sclr_dim > 0 ) then
          write(fstderr,*) "sclrm = ", sclrm
          write(fstderr,*) "wpsclrp = ", wpsclrp
          write(fstderr,*) "sclrp2 = ", sclrp2
          write(fstderr,*) "sclrprtp = ", sclrprtp
          write(fstderr,*) "sclrpthlp = ", sclrpthlp
        end if

        write(fstderr,*) "level = ", level

        write(fstderr,*) "Intent(out)"
  
        write(fstderr,*) "wp4 = ", wp4
        write(fstderr,*) "wprtp2 = ", wprtp2
        write(fstderr,*) "wp2rtp = ", wp2rtp
        write(fstderr,*) "wpthlp2 = ", wpthlp2
        write(fstderr,*) "cloud_frac = ", cloud_frac
        write(fstderr,*) "rcm = ", rcm
        write(fstderr,*) "wpthvp = ", wpthvp
        write(fstderr,*) "wp2thvp = ", wp2thvp
        write(fstderr,*) "rtpthvp = ", rtpthvp
        write(fstderr,*) "thlpthvp = ", thlpthvp
        write(fstderr,*) "wprcp = ", wprcp
        write(fstderr,*) "wp2rcp = ", wp2rcp
        write(fstderr,*) "rtprcp = ", rtprcp
        write(fstderr,*) "thlprcp = ", thlprcp
        write(fstderr,*) "rcp2 = ", rcp2
        write(fstderr,*) "wprtpthlp = ", wprtpthlp
        write(fstderr,*) "crt1 = ", crt1
        write(fstderr,*) "crt2 = ", crt2
        write(fstderr,*) "cthl1 = ", cthl1
        write(fstderr,*) "cthl2 = ", cthl2
        write(fstderr,*) "pdf_params%w1 = ", pdf_params%w1(level)
        write(fstderr,*) "pdf_params%w2 = ", pdf_params%w2(level)
        write(fstderr,*) "pdf_params%varnce_w1 = ", pdf_params%varnce_w1(level)
        write(fstderr,*) "pdf_params%varnce_w2 = ", pdf_params%varnce_w2(level)
        write(fstderr,*) "pdf_params%rt1 = ", pdf_params%rt1(level)
        write(fstderr,*) "pdf_params%rt2 = ", pdf_params%rt2(level)
        write(fstderr,*) "pdf_params%varnce_rt1 = ", pdf_params%varnce_rt1(level)
        write(fstderr,*) "pdf_params%varnce_rt2 = ", pdf_params%varnce_rt2(level)
        write(fstderr,*) "pdf_params%thl1 = ", pdf_params%thl1(level)
        write(fstderr,*) "pdf_params%thl2 = ", pdf_params%thl2(level)
        write(fstderr,*) "pdf_params%varnce_thl1 = ", pdf_params%varnce_thl1(level)
        write(fstderr,*) "pdf_params%varnce_thl2 = ", pdf_params%varnce_thl2(level)
        write(fstderr,*) "pdf_params%a = ", pdf_params%a(level)
        write(fstderr,*) "pdf_params%rrtthl = ", pdf_params%rrtthl(level)
        write(fstderr,*) "pdf_params%rc1 = ", pdf_params%rc1(level)
        write(fstderr,*) "pdf_params%rc2 = ", pdf_params%rc2(level)
        write(fstderr,*) "pdf_params%rsl1 = ", pdf_params%rsl1(level)
        write(fstderr,*) "pdf_params%rsl2 = ", pdf_params%rsl2(level)
        write(fstderr,*) "pdf_params%cloud_frac1 = ", pdf_params%cloud_frac1(level)
        write(fstderr,*) "pdf_params%cloud_frac2 = ", pdf_params%cloud_frac2(level)
        write(fstderr,*) "pdf_params%s1 = ", pdf_params%s1(level)
        write(fstderr,*) "pdf_params%s2 = ", pdf_params%s2(level)
        write(fstderr,*) "pdf_params%stdev_s1 = ", pdf_params%stdev_s1(level)
        write(fstderr,*) "pdf_params%stdev_s2 = ", pdf_params%stdev_s2(level)
        write(fstderr,*) "pdf_params%alpha_thl = ", pdf_params%alpha_thl(level)
        write(fstderr,*) "pdf_params%alpha_rt = ", pdf_params%alpha_rt(level)

        if ( sclr_dim > 0 )then
          write(fstderr,*) "sclrpthvp = ", sclrpthvp
          write(fstderr,*) "sclrprcp = ", sclrprcp
          write(fstderr,*) "wpsclrp2 = ", wpsclrp2
          write(fstderr,*) "wpsclrprtp = ", wpsclrprtp
          write(fstderr,*) "wpsclrpthlp = ", wpsclrpthlp
          write(fstderr,*) "wp2sclrp = ", wp2sclrp
        end if
                
      end if ! err_code == clubb_var_equals_NaN
  
    end if ! clubb_at_least_debug_level

    return
  end subroutine pdf_closure

end module pdf_closure_module
