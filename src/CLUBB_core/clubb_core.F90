!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------
module clubb_core

! Description:
!   The module containing the `core' of the CLUBB model.

! References:
!   None
!-----------------------------------------------------------------------

  implicit none

  public ::  & 
    setup_clubb_core, & 
    advance_clubb_core, & 
    cleanup_clubb_core

  private ! Default Scope

  contains

  !-----------------------------------------------------------------------
  subroutine advance_clubb_core & 
             ( l_implemented, dt, fcor, & 
               thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & 
               sclrm_forcing, edsclrm_forcing, wm_zm, wm_zt, &
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, & 
               wpsclrp_sfc, wpedsclrp_sfc, &
               p_in_Pa, rho_zm, rho, exner, & 
               um, vm, upwp, vpwp, up2, vp2, & 
               thlm, rtm, wprtp, wpthlp, wpthvp, &
               Kh_zt, wp2, wp3, & 
               rtp2, thlp2, rtpthlp, & 
               sigma_sqd_w, tau_zm, rcm, cf, & 
               sclrm, sclrp2, sclrprtp, sclrpthlp, &
               wpsclrp, edsclrm, pdf_params, &
               err_code ) 

    ! Description:
    !   Subroutine to advance the model one timestep

    ! References:
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !     Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    !-----------------------------------------------------------------------

    ! Modules to be included

    use constants, only: & 
      wtol,  & ! Variable(s)
      emin, & 
      thltol, & 
      rttol, &
      wtol_sqd, &
      ep2, & 
      Cp, & 
      Lv, & 
      ep1, & 
      eps, &
      fstderr, &
      zero_threshold

    use parameters_tunable, only: & 
      gamma_coefc,  & ! Variable(s)
      gamma_coefb, & 
      gamma_coef, & 
      taumax, & 
      c_K    

    use parameters_model, only: &
      T0, &
      sclr_dim, &
      edsclr_dim

    use model_flags, only: & 
      l_tke_aniso, & 
      l_gamma_Skw

    use grid_class, only: & 
      gr,  & ! Variable(s)
      zm2zt,  & ! Procedure(s)
      zt2zm, & 
      ddzm

    use numerical_check, only: & 
      parameterization_check ! Procedure(s)

    use variables_diagnostic_module, only: & 
      Skw_zt,  & ! Varible(s)
      Skw_zm, & 
      sigma_sqd_w_zt, & 
      wp4, & 
      thlpthvp, & 
      rtpthvp, & 
      wprcp, & 
      rtprcp, & 
      thlprcp, & 
      rcp2, & 
      rsat, & 
      shear, & 
      wprtp2, & 
      wp2rtp, & 
      wpthlp2, & 
      wp2thlp, & 
      wprtpthlp, & 
      wp2thvp, & 
      wp2rcp, & 
      thvm, & 
      em, & 
      Lscale, & 
      tau_zt, & 
      Kh_zm, & 
      vg, & 
      ug, & 
      um_ref, & 
      vm_ref, & 
      wp2_zt, & 
      thlp2_zt, & 
      wpthlp_zt, & 
      wprtp_zt, & 
      rtp2_zt, & 
      rtpthlp_zt, &
      rtm_ref, &
      thlm_ref

    use variables_diagnostic_module, only: & 
      wpedsclrp, & 
      sclrpthvp,   & ! sclr'th_v'
      sclrprcp,    & ! sclr'rc'
      wp2sclrp,    & ! w'^2 sclr'
      wpsclrp2,    & ! w'sclr'^2
      wpsclrprtp,  & ! w'sclr'rt'
      wpsclrpthlp    ! w'sclr'thl'

    use variables_prognostic_module, only: &
      pdf_parameter

    use advance_xm_wpxp_module, only: & 
      ! Variable(s) 
      advance_xm_wpxp          ! Compute mean/flux terms

    use advance_xp2_xpyp_module, only: & 
      ! Variable(s) 
      advance_xp2_xpyp     ! Computes variance terms

    use surface_var, only:  & 
      sfc_var ! Procedure

    use pdf_closure_module, only:  & 
      ! Procedure 
      pdf_closure     ! Prob. density function

    use mixing_length, only: & 
      compute_length ! Procedure

    use advance_windm_edsclrm_module, only:  & 
      advance_windm_edsclrm  ! Procedure(s)

    use saturation, only:  & 
      ! Procedure
      sat_mixrat_liq ! Saturation mixing ratio

    use advance_wp2_wp3_module, only:  & 
      advance_wp2_wp3 ! Procedure

    use stats_precision, only:  & 
      time_precision ! Variable(s)

    use error_code, only :  & 
      clubb_var_equals_NaN, & ! Variable(s)
      clubb_var_out_of_bounds, & 
      lapack_error,  & ! Procedure(s)
      clubb_at_least_debug_level

    use Skw_module, only:  & 
      Skw_func ! Procedure

    use clip_explicit, only: & 
      clip_covariances_denom ! Procedure(s)

    use T_in_K_mod, only: &
      ! Read values from namelist
      thlm2T_in_K ! Procedure

    use stats_subs, only: & 
      stats_accumulate ! Procedure

    use stats_type, only: & 
      stat_update_var_pt ! Procedure(s)

    use stats_variables, only: &
      ishear,     & ! Variables
      ircp2,      &
      iwp4,       &
      irsat,      &
      iwprtp_zt,  &
      iwpthlp_zt

    implicit none

    ! External
    intrinsic :: sqrt, min, max, exp, mod

    ! Input Variables
    logical, intent(in) ::  & 
      l_implemented ! Is this part of a larger host model (T/F) ?

    ! Note on dt, dmain, and dtclosure: since being moved out of
    ! hoc.F, all subroutines within advance_clubb_core now use
    ! dt for time dependent calculations.  The old dt is noted in
    ! each section of the code -dschanen 20 April 2006
    real(kind=time_precision), intent(in) ::  & 
      dt            ! Current timestep size    [s]

    real, intent(in) ::  & 
      fcor          ! Coriolis forcing         [s^-1]

    real, intent(in), dimension(gr%nnzp) ::  & 
      thlm_forcing,   & ! theta_l forcing.        [K/s]
      rtm_forcing,    & ! r_t forcing.            [(kg/kg)/s]
      um_forcing,     & ! u wind forcing          [m/s/s]
      vm_forcing,     & ! v wind forcing          [m/s/s]
      wm_zm,          & ! wm on moment. grid.     [m/s]
      wm_zt,          & ! wm on thermo. grid.     [m/s]
      p_in_Pa,        & ! Pressure.               [Pa] 
      rho_zm,         & ! Density on moment. grid [kg/m^3]
      rho,            & ! Density on thermo. grid [kg/m^3] 
      exner             ! Exner function.         [-]

    real, intent(in) ::  & 
      wpthlp_sfc,   & ! w' theta_l' at surface.   [(m K)/s]
      wprtp_sfc,    & ! w' r_t' at surface.       [(kg m)/( kg s)]
      upwp_sfc,     & ! u'w' at surface.          [m^2/s^2]
      vpwp_sfc        ! v'w' at surface.          [m^2/s^2]

    real, intent(in),  dimension(sclr_dim) ::  & 
      wpsclrp_sfc      ! Scalar flux at surface           [units m/s]

    real, intent(in),  dimension(edsclr_dim) ::  & 
      wpedsclrp_sfc    ! Eddy-Scalar flux at surface      [units m/s]

    ! Input/Output
    ! These are prognostic or are planned to be in the future
    real, intent(inout), dimension(gr%nnzp) ::  & 
      um,         & ! u wind.                       [m/s]
      upwp,       & ! u'w'.                         [m^2/s^2]
      vm,         & ! v wind.                       [m/s]
      vpwp,       & ! u'w'.                         [m^2/s^2]
      up2,        & ! u'^2                          [m^2/s^2]
      vp2,        & ! v'^2                          [m^2/s^2]
      rtm,        & ! r_t Total water mixing ratio. [kg/kg]
      wprtp,      & ! w' r_t'.                      [(m kg)/(s kg)]
      thlm,       & ! th_l Liquid potential temp.   [K]
      wpthlp,     & ! w' th_l'.                     [(m K)/s]
      wpthvp,     & ! w' th_v'.                     [(m K)/s]
      Kh_zt,      & ! Eddy-diffusivity              [m^2/s]
      wp2,        & ! w'^2.                         [m^2/s^2]
      wp3,        & ! w'^3.                         [m^3/s^3]
      sigma_sqd_w,& ! sigma_sqd_w on moment. grid.           [-]
      rtp2,       & ! r_t'^2.                       [(kg/kg)^2]
      thlp2,      & ! th_l'^2.                      [K^2]
      rtpthlp,    & ! r_t' th_l'.                   [(kg K)/kg]
      tau_zm,     & ! Tau on moment. grid.          [s]
      rcm           ! Liquid water mixing ratio.    [kg/kg]

    ! Needed for output for host models
    real, intent(inout), dimension(gr%nnzp) ::  & 
      cf ! Cloud fraction.     [%]

    ! Diagnostic, for if some calculation goes amiss.
    integer, intent(inout) :: err_code

    ! Passive scalar variables
    real, intent(inout), dimension(gr%nnzp,sclr_dim) :: & 
      sclrm,         & ! Passive scalar mean.           [units vary]
      sclrp2,        & ! Passive scalar variance.       [{units vary}^2]
      sclrprtp,      & ! sclr'rt'                       [{units vary}^2]
      sclrpthlp,     & ! sclr'thl'                      [{units vary}^2]
      wpsclrp          ! w'sclr'                        [units vary m/s]

    real, intent(in), dimension(gr%nnzp,sclr_dim) :: & 
      sclrm_forcing    ! Passive scalar forcing.        [{units vary}/s]

    real, intent(inout), dimension(gr%nnzp,edsclr_dim) :: & 
      edsclrm          ! Eddy passive scalar mean.      [units vary]

    real, intent(in), dimension(gr%nnzp,edsclr_dim) :: & 
      edsclrm_forcing    ! Eddy passive scalar forcing. [{units vary}/s]

    type(pdf_parameter), intent(inout) :: & 
      pdf_params ! PDF parameters [units vary]

    ! Local Variables
    integer :: i, k

    real, dimension(gr%nnzp) :: &
      tmp1, gamma_Skw_fnc

    real, dimension(gr%nnzp,sclr_dim) :: & 
      sclr_tmp1, sclr_tmp2, sclr_tmp3, sclr_tmp4 ! for PDF closure

    integer :: &
      wprtp_cl_num,   & ! Instance of w'r_t' clipping (1st or 3rd).
      wpthlp_cl_num,  & ! Instance of w'th_l' clipping (1st or 3rd).
      wpsclrp_cl_num, & ! Instance of w'sclr' clipping (1st or 3rd).
      upwp_cl_num,    & ! Instance of u'w' clipping (1st or 2nd).
      vpwp_cl_num       ! Instance of v'w' clipping (1st or 2nd).

    ! ldgrant June 2009
    logical, parameter :: &
      l_trapezoid_rule = .false., &     ! Logical flag used to turn the trapezoidal 
                                        ! rule on or off.  

      l_call_pdf_closure_twice = .true. ! If the trapezoidal rule is used, this logical 
                                        ! flag determines if the variables originally
                                        ! on the thermodynamic grid will be brought
                                        ! to the momentum grid for use with the trapezoidal
                                        ! rule by calling the subroutine pdf_closure a 
                                        ! second time (.true.) or by interpolation 
                                        ! (.false.).  Calling pdf_closure twice is more 
                                        ! expensive but produces better results.

    ! ldgrant July 2009
    real, dimension(gr%nnzp) :: &
      cloud_cover,  & ! Cloud cover                               [%]
      rcm_in_layer    ! liquid water mixing ratio in cloud layer  [kg/kg]

    !----- Begin Code -----

    !----------------------------------------------------------------
    ! Test input variables
    !----------------------------------------------------------------
    if ( clubb_at_least_debug_level( 2 ) ) then
      call parameterization_check & 
           ( thlm_forcing, rtm_forcing, wm_zm, wm_zt, p_in_Pa, rho_zm, & ! intent(in)
           rho, exner, wpthlp_sfc, wprtp_sfc,                & ! intent(in)
           upwp_sfc, vpwp_sfc, um, upwp, vm, vpwp,           & ! intent(in)
           up2, vp2, rtm, wprtp, thlm,                       & ! intent(in)
           wpthlp, wp2, wp3, sigma_sqd_w, rtp2, thlp2,       & ! intent(in)
           rtpthlp, tau_zm, rcm, cf, "beginning of ",        & ! intent(in)
           wpsclrp_sfc, wpedsclrp_sfc,                       & ! intent(in)
           sclrm, sclrm_forcing, edsclrm, edsclrm_forcing )    ! intent(in)
    end if
    !-----------------------------------------------------------------------

    ! SET SURFACE VALUES OF FLUXES (BROUGHT IN)
    wpthlp(1) = wpthlp_sfc
    wprtp(1)  = wprtp_sfc
    upwp(1)   = upwp_sfc
    vpwp(1)   = vpwp_sfc

    ! Set fluxes for passive scalars (if enabled)
    if ( sclr_dim > 0 ) then
      wpsclrp(1,1:sclr_dim)   = wpsclrp_sfc(1:sclr_dim)
    end if

    if ( edsclr_dim > 0 ) then
      wpedsclrp(1,1:edsclr_dim) = wpedsclrp_sfc(1:edsclr_dim)
    end if

    !----------------------------------------------------------------
    ! Set Surface variances
    !----------------------------------------------------------------

!      Surface variances should be set here, before the call to advance_xp2_xpyp.
!      The reasons that surface variances can be set here are because the
!      only variables that are the input into surface variances are the
!      surface values of wpthlp, wprtp, upwp, and vpwp.  The surface values
!      of all those variables are set in the surface forcings section of the
!      GCSS cases listed in the main timestep above.  Even if they weren't
!      set there, the updates to wpthlp, wprtp, upwp, and vpwp are at the
!      end of the closure loop, right before the code loops back around to
!      this point at the top of the closure loop.
!      Surface variances need to be set here for two reasons.  One reason is
!      that the values of rtp2, thlp2, and rtpthlp at the surface will be
!      used to find the diffusional term and the mean advection term in each
!      predictive equation for those respective terms.  The other reason is
!      that if the correct surface variances are not set here and advance_xp2_xpyp
!      outputs it's own value for them, it will results in a faulty value for
!      sigma_sqd_w at the surface.  Brian Griffin.  December 18, 2005.

!      Surface effects should not be included with any case where the lowest
!      level is not the ground level.  Brian Griffin.  December 22, 2005.
    IF ( gr%zm(1) == 0.0 ) THEN
      call sfc_var( upwp(1), vpwp(1), wpthlp(1), wprtp(1),   & ! intent(in)
                    wpsclrp(1,1:sclr_dim),                   & ! intent(in)
                    wp2(1), up2(1), vp2(1),                  & ! intent(out)
                    thlp2(1), rtp2(1), rtpthlp(1), err_code, & ! intent(out)
                    sclrp2(1,1:sclr_dim),                    & ! intent(out)
                    sclrprtp(1,1:sclr_dim),                  & ! intent(out) 
                    sclrpthlp(1,1:sclr_dim) )                  ! intent(out)
      ! Subroutine may produce NaN values, and if so, exit
      ! gracefully.
      ! Joshua Fasching March 2008
      if( err_code == clubb_var_equals_NaN ) return

    ELSE
      ! Variances for cases where the lowest level is not at the surface.
      ! Eliminate surface effects on lowest level variances.
      wp2(1)     = wtol_sqd
      up2(1)     = wtol_sqd
      vp2(1)     = wtol_sqd
      thlp2(1)   = 0.0
      rtp2(1)    = 0.0
      rtpthlp(1) = 0.0

      DO i = 1, sclr_dim, 1
        sclrp2(1,i)    = 0.0
        sclrprtp(1,i)  = 0.0
        sclrpthlp(1,i) = 0.0
      END DO
    END IF


    !----------------------------------------------------------------
    ! Interpolate wp2 & wp3, and then compute Skw for m & t grid
    !----------------------------------------------------------------

    wp2_zt = max( zm2zt( wp2 ), wtol_sqd ) ! Positive definite quantity

    Skw_zt(1:gr%nnzp) = Skw_func( wp2_zt(1:gr%nnzp), wp3(1:gr%nnzp) )
    Skw_zm(1:gr%nnzp) = zt2zm( Skw_zt(1:gr%nnzp) ) 

    !----------------------------------------------------------------
    ! Diagnose variances
    !----------------------------------------------------------------

    ! We also found that certain cases require a time tendency to run
    ! at shorter timesteps.
    ! This requires us to store in memory sigma_sqd_w and tau_zm between timesteps.

    ! We found that if we call advance_xp2_xpyp first, we can use a longer timestep.
    call advance_xp2_xpyp( tau_zm, wm_zm, rtm, wprtp,     & ! intent(in)
                           thlm, wpthlp, wpthvp, um, vm,  & ! intent(in)
                           wp2, wp2_zt, wp3, upwp, vpwp,  & ! intent(in)
                           sigma_sqd_w, Skw_zm, Kh_zt,    & ! intent(in)
 ! Vince Larson used prognostic timestepping of variances 
 !    in order to increase numerical stability.  17 Jul 2007
 !                          .false., dt,                   & ! intent(in)
                           .true., dt,                    & ! intent(in)
                           sclrm, wpsclrp,                & ! intent(in) 
                           rtp2, thlp2, rtpthlp,          & ! intent(inout)
                           up2, vp2,                      & ! intent(inout)
                           err_code,                      & ! intent(out)
                           sclrp2, sclrprtp, sclrpthlp  )   ! intent(out)

    ! Iterpolate variances to the zt grid (statistics and closure)
    thlp2_zt   = max( zm2zt( thlp2 ), thltol**2 )  ! Positive def. quantity
    rtp2_zt    = max( zm2zt( rtp2 ), rttol**2 )   ! Positive def. quantity
    rtpthlp_zt = zm2zt( rtpthlp )

    ! Check stability
    ! Changed from a logical flag to an integer indicating nature of
    ! error.
    ! Joshua Fasching March 2008
    if ( lapack_error( err_code ) ) return


    !----------------------------------------------------------------
    ! Covariance clipping for wprtp, wpthlp, wpsclrp, upwp, and vpwp
    ! after subroutine advance_xp2_xpyp updated xp2.
    !----------------------------------------------------------------

    wprtp_cl_num   = 1 ! First instance of w'r_t' clipping.
    wpthlp_cl_num  = 1 ! First instance of w'th_l' clipping.
    wpsclrp_cl_num = 1 ! First instance of w'sclr' clipping.
    upwp_cl_num    = 1 ! First instance of u'w' clipping.
    vpwp_cl_num    = 1 ! First instance of v'w' clipping.

    call clip_covariances_denom( dt, rtp2, thlp2, up2, vp2, wp2, &
                                 sclrp2, wprtp_cl_num, wpthlp_cl_num, &
                                 wpsclrp_cl_num, upwp_cl_num, vpwp_cl_num, &
                                 wprtp, wpthlp, upwp, vpwp, wpsclrp )


    ! The right hand side of this conjunction is only for reducing cpu time,
    ! since the more complicated formula is mathematically equivalent
    if ( l_gamma_Skw .and. ( gamma_coef /= gamma_coefb ) ) then
      !----------------------------------------------------------------
      ! Compute gamma as a function of Skw  - 14 April 06 dschanen
      !----------------------------------------------------------------

      gamma_Skw_fnc = gamma_coefb + (gamma_coef-gamma_coefb) & 
            *exp( -(1.0/2.0) * (Skw_zm/gamma_coefc)**2 )

    else
      gamma_Skw_fnc = gamma_coef

    end if

    !----------------------------------------------------------------
    ! Compute sigma_sqd_w with new formula from Vince
    !----------------------------------------------------------------

    sigma_sqd_w = gamma_Skw_fnc * &
      ( 1.0 - min( & 
                  max( ( wpthlp / ( sqrt( wp2 * thlp2 )  & 
                      + 0.01 * wtol * thltol ) )**2, & 
                       ( wprtp / ( sqrt( wp2 * rtp2 )  & 
                      + 0.01 * wtol * rttol ) )**2 &
                     ), & ! max  
             1.0 ) & ! min 
       )

    sigma_sqd_w_zt = max( zm2zt( sigma_sqd_w ), zero_threshold )  ! Pos. def. quantity

    !----------------------------------------------------------------
    ! Call closure scheme
    !----------------------------------------------------------------

    ! Put passive scalar input on the t grid for the PDF
    do i = 1, sclr_dim, 1
      sclr_tmp1(:,i) = zm2zt( wpsclrp(:,i) )
      sclr_tmp2(:,i) = zm2zt( sclrprtp(:,i) )
      sclr_tmp3(:,i) = max( zm2zt( sclrp2(:,i) ), zero_threshold ) ! Pos. def. quantity
      sclr_tmp4(:,i) = zm2zt( sclrpthlp(:,i) )
    end do ! i = 1, sclr_dim

    do k = 1, gr%nnzp, 1
      call pdf_closure & 
         ( p_in_Pa(k), exner(k), wm_zt(k), wp2_zt(k), wp3(k), sigma_sqd_w_zt(k), & ! intent(in)
           Skw_zt(k), rtm(k), rtp2_zt(k), zm2zt( wprtp, k ),       & ! intent(in)
           thlm(k), thlp2_zt(k), zm2zt( wpthlp, k ),               & ! intent(in)
           rtpthlp_zt(k), sclrm(k,:), sclr_tmp1(k,:),              & ! intent(in)
           sclr_tmp3(k,:),sclr_tmp2(k,:), sclr_tmp4(k,:), k,       & ! intent(in)
           wp4(k), wprtp2(k), wp2rtp(k),                           & ! intent(out)
           wpthlp2(k), wp2thlp(k), wprtpthlp(k),                   & ! intent(out)
           cf(k), rcm(k), wpthvp(k), wp2thvp(k), rtpthvp(k),       & ! intent(out)
           thlpthvp(k), wprcp(k), wp2rcp(k), rtprcp(k), thlprcp(k),& ! intent(out)
           rcp2(k), pdf_params,                                    & ! intent(out)
           err_code,                                               & ! intent(out)
           wpsclrprtp(k,:), wpsclrp2(k,:), sclrpthvp(k,:),         & ! intent(out)
           wpsclrpthlp(k,:), sclrprcp(k,:), wp2sclrp(k,:) )          ! intent(out)

      ! Subroutine may produce NaN values, and if so, exit
      ! gracefully.
      ! Joshua Fasching March 2008

      if ( err_code == clubb_var_equals_NaN ) then
        write(fstderr,*) "At grid level = ",k
        return
      end if

    end do ! k = 2, gr%nnzp-1

    ! Interpolate momentum variables back to momentum grid.
    ! Since top momentum level is higher than top thermo level,
    ! Set variables at top momentum level to 0.

    ! Only do this for wp4 and rcp2 if we're saving stats, since they are not
    ! used elsewhere in the parameterization 
    if ( iwp4 > 0 ) then
      wp4 = max( zt2zm( wp4 ), zero_threshold )  ! Pos. def. quantity
      wp4(gr%nnzp)  = 0.0
    end if

    if ( ircp2 > 0 ) then
      rcp2 = max( zt2zm( rcp2 ), zero_threshold )  ! Pos. def. quantity
      rcp2(gr%nnzp) = 0.0
    end if

    wpthvp            = zt2zm( wpthvp )
    wpthvp(gr%nnzp)   = 0.0
    thlpthvp          = zt2zm( thlpthvp )
    thlpthvp(gr%nnzp) = 0.0
    rtpthvp           = zt2zm( rtpthvp )
    rtpthvp(gr%nnzp)  = 0.0
    wprcp             = zt2zm( wprcp )
    wprcp(gr%nnzp)    = 0.0
    rtprcp            = zt2zm( rtprcp )
    rtprcp(gr%nnzp)   = 0.0
    thlprcp           = zt2zm( thlprcp )
    thlprcp(gr%nnzp)  = 0.0

    ! Interpolate passive scalars back onto the m grid
    do i = 1, sclr_dim
      sclrpthvp(:,i)       = zt2zm( sclrpthvp(:,i) )
      sclrpthvp(gr%nnzp,i) = 0.0
      sclrprcp(:,i)        = zt2zm( sclrprcp(:,i) )
      sclrprcp(gr%nnzp,i)  = 0.0
    end do ! i=1, sclr_dim

    ! If logical flag is true, use the trapezoidal rule by calling the subroutine
    ! ldgrant June 2009
    if ( l_trapezoid_rule ) then
      call trapezoidal_rule &
         ( l_call_pdf_closure_twice,                    & ! intent(in)
           p_in_Pa, exner, wm_zm, wp2, wp3,             & ! intent(in)
           sigma_sqd_w, Skw_zm, rtm, rtp2,              & ! intent(in)
           wprtp, thlm, thlp2, wpthlp, rtpthlp,         & ! intent(in)
           sclrm, wpsclrp, sclrp2, sclrprtp, sclrpthlp, & ! intent(in)
           wprtp2, wp2rtp, wpthlp2, wp2thlp,            & ! intent(inout)
           wprtpthlp, cf, rcm, wp2thvp, wp2rcp,         & ! intent(inout)
           wpsclrprtp, wpsclrp2, wpsclrpthlp,           & ! intent(inout)
           wp2sclrp, pdf_params, err_code )               ! intent(inout)
    end if

    ! Compute variables cloud_cover and rcm_in_layer
    ! ldgrant July 2009
    call compute_cloud_cover &
       ( pdf_params, cf, rcm,        & ! intent(in)
         cloud_cover, rcm_in_layer )   ! intent(out)

    !----------------------------------------------------------------
    ! Compute thvm
    !----------------------------------------------------------------

    thvm = thlm + ep1 * T0 * rtm + ( Lv/(Cp*exner) - ep2 * T0 ) * rcm

    !----------------------------------------------------------------
    ! Compute tke (turbulent kinetic energy)
    !----------------------------------------------------------------

    if ( .not. l_tke_aniso ) then
      ! tke is assumed to be 3/2 of wp2
      em = 1.5 * wp2
    else
      em = 0.5 * ( wp2 + vp2 + up2 )
    end if

    !----------------------------------------------------------------
    ! Compute mixing length
    !----------------------------------------------------------------

    call compute_length( thvm, thlm, rtm, rcm, em, p_in_Pa, exner, &    ! intent(in)
                         err_code, &                                    ! intent(inout)
                         Lscale )                                       ! intent(out)

    ! Subroutine may produce NaN values, and if so, exit
    ! gracefully.
    ! Joshua Fasching March 2008
    if( err_code == clubb_var_equals_NaN ) return


    !----------------------------------------------------------------
    ! Dissipation time
    !----------------------------------------------------------------
! Vince Larson replaced the cutoff of emin by wtol**2.  7 Jul 2007
!     This is to prevent tau from being too large (producing little damping)
!     in stably stratified layers with little turbulence.
!       tmp1 = SQRT( MAX( emin, zm2zt( em ) ) )
!       tau_zt = MIN( Lscale / tmp1, taumax )
!       tau_zm &
!       = MIN( ( zt2zm( Lscale ) / SQRT( MAX( emin, em ) ) ), taumax )
!   Addition by Brian:  Model constant emin is now set to (3/2)*wtol_sqd.
!                       Thus, emin can replace wtol_sqd here.
    tmp1   = SQRT( MAX( emin, zm2zt( em ) ) )
    tau_zt = MIN( Lscale / tmp1, taumax )
    tau_zm = MIN( ( MAX( zt2zm( Lscale ), zero_threshold )  & 
                   / SQRT( MAX( emin, em ) ) ), taumax )
! End Vince Larson's replacement.

    ! Modification to damp noise in stable region
! Vince Larson commented out because it may prevent turbulence from
!    initiating in unstable regions.  7 Jul 2007
!       do k = 1, gr%nnzp
!         if ( wp2(k) <= 0.005 ) then
!           tau_zt(k) = taumin
!           tau_zm(k) = taumin
!         end if
!       end do
! End Vince Larson's commenting.

    !----------------------------------------------------------------
    ! Eddy diffusivity coefficient
    !----------------------------------------------------------------
    ! c_K is 0.548 usually (Duynkerke and Driedonks 1987)

    Kh_zt = c_K * Lscale * tmp1
    Kh_zm = c_K * max( zt2zm( Lscale ), zero_threshold )  & 
                * sqrt( max( em, emin ) )

    !#######################################################################
    !############## ADVANCE PROGNOSTIC VARIABLES ONE TIMESTEP ##############
    !#######################################################################


    ! Store the saturation mixing ratio for output purposes.  Brian
    if ( irsat > 0 ) then
      rsat = sat_mixrat_liq( p_in_Pa, thlm2T_in_K( thlm, exner, rcm ) )
    end if

    !----------------------------------------------------------------
    ! Advance rtm/wprtp and thlm/wpthlp one time step
    !----------------------------------------------------------------
    call advance_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, wp3, & ! intent(in)
                          Kh_zt, tau_zm, Skw_zm, rtpthvp,          & ! intent(in)
                          rtm_forcing, thlpthvp, rtm_ref, thlm_ref, & ! intent(in)
                          thlm_forcing, rtp2, thlp2, wp2_zt,       & ! intent(in)
                          pdf_params, l_implemented,               & ! intent(in)
                          sclrpthvp, sclrm_forcing, sclrp2,        & ! intent(in)
                          rtm, wprtp, thlm, wpthlp,                & ! intent(inout)
                          err_code,                                & ! intent(inout)
                          sclrm, wpsclrp )                           ! intent(inout)

    ! Wrapped LAPACK procedures may report errors, and if so, exit
    ! gracefully.
    ! Joshua Fasching March 2008

    if ( lapack_error( err_code ) ) return


    ! Vince Larson clipped rcm in order to prevent rvm < 0.  5 Apr 2008.
    ! This code won't work unless rtm >= 0 !!!
    do k = 1, gr%nnzp
      if ( rtm(k) < rcm(k) ) then

        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,*) 'rtm < rcm in advance_xm_wpxp at k=', k, '.', & 
            '  Clipping rcm.'

        end if ! clubb_at_least_debug_level(1)

        rcm(k) = max( zero_threshold, rtm(k) - eps )

      end if ! rtm(k) < rcm(k)

    end do ! k=1..gr%nnzp


    !----------------------------------------------------------------
    ! Advance wp2/wp3 one timestep
    !----------------------------------------------------------------

    call advance_wp2_wp3 &
         ( dt, sigma_sqd_w, wm_zm, wm_zt, wpthvp, wp2thvp,  & ! intent(in)
           um, vm, upwp, vpwp, up2, vp2, Kh_zm, Kh_zt,      & ! intent(in)
           tau_zm, tau_zt, Skw_zm, Skw_zt, pdf_params%a,    & ! intent(in)
           wp2_zt, wp2, wp3, err_code )                       ! intent(inout)

    ! Wrapped LAPACK procedures may report errors, and if so, exit
    ! gracefully.
    ! Joshua Fasching March 2008
    if ( lapack_error( err_code ) ) return


    !----------------------------------------------------------------
    ! Covariance clipping for wprtp, wpthlp, wpsclrp, upwp, and vpwp
    ! after subroutine advance_wp2_wp3 updated wp2.
    !----------------------------------------------------------------

    wprtp_cl_num   = 3 ! Third instance of w'r_t' clipping.
    wpthlp_cl_num  = 3 ! Third instance of w'th_l' clipping.
    wpsclrp_cl_num = 3 ! Third instance of w'sclr' clipping.
    upwp_cl_num    = 2 ! Second instance of u'w' clipping.
    vpwp_cl_num    = 2 ! Second instance of v'w' clipping.

    call clip_covariances_denom( dt, rtp2, thlp2, up2, vp2, wp2, &
                                 sclrp2, wprtp_cl_num, wpthlp_cl_num, &
                                 wpsclrp_cl_num, upwp_cl_num, vpwp_cl_num, &
                                 wprtp, wpthlp, upwp, vpwp, wpsclrp )


    !----------------------------------------------------------------
    ! Advance um, vm, and edsclrm one time step
    !----------------------------------------------------------------

    call advance_windm_edsclrm( dt, wm_zt, Kh_zm, ug, vg, um_ref, vm_ref,  & ! In
                                wp2, up2, vp2, um_forcing, vm_forcing, edsclrm_forcing, & ! In
                                upwp_sfc, vpwp_sfc, wpedsclrp_sfc, fcor,  &  !  In
                                l_implemented, um, vm, edsclrm, &
                                upwp, vpwp, wpedsclrp, err_code )

    ! Wrapped LAPACK procedures may report errors, and if so, exit
    ! gracefully.
    ! Joshua Fasching March 2008
    if ( lapack_error( err_code ) ) return


    ! Compute Shear Production  -Brian
    ! This is a non-interative diagnostic, for statistical purposes
    if ( ishear > 1  ) then
      do k = 1, gr%nnzp-1, 1
        shear(k) = -upwp(k) * ( um(k+1) - um(k) ) * gr%dzm(k) & 
                   -vpwp(k) * ( vm(k+1) - vm(k) ) * gr%dzm(k)
      end do
      shear(gr%nnzp) = 0.0
    end if

    !#######################################################################
    !#############            ACCUMULATE STATISTICS            #############
    !#######################################################################

    if ( iwprtp_zt > 0 ) then
      wpthlp_zt  = zm2zt( wpthlp )
    end if

    if ( iwpthlp_zt > 0 ) then
      wprtp_zt   = zm2zt( wprtp )
    end if

    call stats_accumulate & 
         ( um, vm, upwp, vpwp, up2, vp2, thlm,                 & ! intent(in)
           rtm, wprtp, wpthlp, wpthvp,                         & ! intent(in) 
           wp2, wp3, rtp2, thlp2, rtpthlp,                     & ! intent(in)
           p_in_Pa, exner, rho, rho_zm, Kh_zt,                 & ! intent(in)
           wm_zt, sigma_sqd_w, tau_zm, rcm, cf,                & ! intent(in)
           cloud_cover, rcm_in_layer,                          & ! intent(in)
           pdf_params,                                         & ! intent(in)
           sclrm, sclrp2, sclrprtp, sclrpthlp, sclrm_forcing,  & ! intent(in)
           wpsclrp, edsclrm, edsclrm_forcing )                   ! intent(in)


    if ( clubb_at_least_debug_level( 2 ) ) then
      call parameterization_check & 
           ( thlm_forcing, rtm_forcing, wm_zm, wm_zt, p_in_Pa, rho_zm, & ! intent(in)
             rho, exner, wpthlp_sfc, wprtp_sfc,                  & ! intent(in)
             upwp_sfc, vpwp_sfc, um, upwp, vm, vpwp,             & ! intent(in)
             up2, vp2, rtm, wprtp, thlm,                         & ! intent(in)
             wpthlp, wp2, wp3, sigma_sqd_w, rtp2, thlp2,         & ! intent(in)
             rtpthlp, tau_zm, rcm, cf, "end of ",                & ! intent(in)
             wpsclrp_sfc, wpedsclrp_sfc,                         & ! intent(in)
             sclrm, sclrm_forcing, edsclrm, edsclrm_forcing )      ! intent(in)
    end if

    return
  end subroutine advance_clubb_core

  !-----------------------------------------------------------------------
  subroutine setup_clubb_core & 
             ( nzmax, T0_in, ts_nudge_in, & ! In
               hydromet_dim_in, sclr_dim_in, & ! In
               sclrtol_in, edsclr_dim_in, params,  &  ! In
               l_soil_veg, & ! In
               l_uv_nudge, l_tke_aniso, saturation_formula, &  ! In
               l_implemented, grid_type, deltaz, zm_init, zm_top, &  ! In
               momentum_heights, thermodynamic_heights,  &  ! In
               host_dx, host_dy, & ! In
               err_code ) ! Out
    !
    !   Description:
    !     Subroutine to set up the model for execution.
    !
    !-------------------------------------------------------------------------
    use grid_class, only: & 
      setup_grid, & ! Procedure
      gr ! Variable(s)

    use parameter_indices, only:  & 
      nparams ! Variable(s)

    use parameters_tunable, only: & 
      setup_parameters ! Procedure

    use parameters_model, only: & 
      setup_parameters_model ! Procedure

    use variables_diagnostic_module, only: & 
      setup_diagnostic_variables ! Procedure

    use variables_prognostic_module, only: & 
      setup_prognostic_variables ! Procedure

    use constants, only:  & 
      fstderr  ! Variable(s)

    use error_code, only:  & 
      clubb_var_out_of_bounds ! Variable(s)

    use model_flags, only: & 
      setup_model_flags ! Subroutine

    implicit none

    ! Input Variables

    ! Grid definition
    integer, intent(in) :: nzmax  ! Vertical grid levels            [#]
    !                      Only true when used in a host model
    !                      CLUBB determines what nzmax should be
    !                      given zm_init and zm_top when
    !                      running in standalone mode.

    ! Flag to see if CLUBB is running on it's own,
    ! or if it's implemented as part of a host model.
    logical, intent(in) :: l_implemented   ! (T/F)

    ! If CLUBB is running on it's own, this option determines
    ! if it is using:
    ! 1) an evenly-spaced grid,
    ! 2) a stretched (unevenly-spaced) grid entered on the
    !    thermodynamic grid levels (with momentum levels set
    !    halfway between thermodynamic levels), or
    ! 3) a stretched (unevenly-spaced) grid entered on the
    !    momentum grid levels (with thermodynamic levels set
    !    halfway between momentum levels).
    integer, intent(in) :: grid_type

    ! If the CLUBB model is running by itself, and is using an
    ! evenly-spaced grid (grid_type = 1), it needs the vertical
    ! grid spacing, momentum-level starting altitude, and maximum 
    ! altitude as input.
    real, intent(in) :: & 
      deltaz,   & ! Change in altitude per level           [m]
      zm_init,  & ! Initial grid altitude (momentum level) [m]
      zm_top      ! Maximum grid altitude (momentum level) [m]

    ! If the CLUBB parameterization is implemented in a host model,
    ! it needs to use the host model's momentum level altitudes
    ! and thermodynamic level altitudes.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on thermodynamic levels (grid_type = 2),
    ! it needs to use the thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on momentum levels (grid_type = 3),
    ! it needs to use the momentum level altitudes as input.
    real, intent(in), dimension(nzmax) :: & 
      momentum_heights,      & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights    ! Thermodynamic level altitudes (input) [m]

    ! Host model horizontal grid spacing, if part of host model.
    real, intent(in) :: & 
      host_dx,  & ! East-West horizontal grid spacing     [m]
      host_dy     ! North-South horizontal grid spacing   [m]

    ! Model parameters
    real, intent(in) ::  & 
      T0_in, ts_nudge_in

    integer, intent(in) :: & 
      hydromet_dim_in,  & ! Number of hydrometeor species
      sclr_dim_in,      & ! Number of passive scalars
      edsclr_dim_in       ! Number of eddy-diff. passive scalars

    real, intent(in), dimension(sclr_dim_in) :: & 
      sclrtol_in    ! Thresholds for passive scalars

    real, intent(in), dimension(nparams) :: & 
      params  ! Including C1, nu1, nu2, etc.

    ! Flags 
    logical, intent(in) ::  & 
      l_soil_veg,     & ! Simple surface scheme
      l_uv_nudge,     & ! Wind nudging
      l_tke_aniso       ! For anisotropic turbulent kinetic energy,
                        !   i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)

    character(len=*), intent(in) :: &
      saturation_formula ! "bolton" approx. or "flatau" approx.

    ! Output variables
    integer, intent(out) :: & 
      err_code   ! Diagnostic for a problem with the setup

    ! Local variables
    real :: Lscale_max
    integer :: begin_height, end_height

    !----- Begin Code -----

    ! Sanity check for the saturation formula
    select case ( trim( saturation_formula ) )
    case ( "bolton", "Bolton" )
      ! Using the Bolton 1980 approximations for SVP over vapor/ice

    case ( "flatau", "Flatau" )
      ! Using the Flatau, et al. polynomial approximation for SVP over vapor/ice

      ! Add new cases after this
    case default
      write(fstderr,*) "Error in setup_clubb_core."
      write(fstderr,*) "Unknown approx. of saturation vapor pressure: "// &
        trim( saturation_formula )
      stop
    end select

    ! Setup grid
    call setup_grid( nzmax, l_implemented, grid_type,           & ! intent(in)
                     deltaz, zm_init, zm_top, momentum_heights, & ! intent(in)
                     thermodynamic_heights,                     & ! intent(in)
                     begin_height, end_height )                   ! intent(in)

    ! Setup flags

    call setup_model_flags & 
         ( l_soil_veg, l_uv_nudge,  & ! intent(in)
           l_tke_aniso, saturation_formula )     ! intent(in)

    ! Determine the maximum allowable value for Lscale (in meters).
    if ( l_implemented ) then
      Lscale_max = 0.25 * min( host_dx, host_dy )
    else
      Lscale_max = 1.0e5
    end if

    ! Define model constant parameters
    call setup_parameters_model( T0_in, ts_nudge_in, &      ! In
                                 hydromet_dim_in, &  ! in
                                 sclr_dim_in, sclrtol_in, edsclr_dim_in, &! In
                                 Lscale_max )   ! In

    ! Define tunable constant parameters
    call setup_parameters & 
         ( deltaz, params, gr%nnzp, l_implemented,               & ! intent(in)
           grid_type, momentum_heights(begin_height:end_height), & ! intent(in)
           thermodynamic_heights(begin_height:end_height),       & ! intent(in)
           err_code )                                              ! intent(out)

    ! Error Report
    ! Joshua Fasching February 2008
    if ( err_code == clubb_var_out_of_bounds ) then

      write(fstderr,*) "Error in setup_clubb_core"

      write(fstderr,*) "Intent(in)"

      write(fstderr,*) "deltaz = ", deltaz
      write(fstderr,*) "zm_init = ", zm_init
      write(fstderr,*) "zm_top = ", zm_top
      write(fstderr,*) "momentum_heights = ", momentum_heights
      write(fstderr,*) "thermodynamic_heights = ",  & 
        thermodynamic_heights
      write(fstderr,*) "T0_in = ", T0_in
      write(fstderr,*) "ts_nudge_in = ", ts_nudge_in
      write(fstderr,*) "params = ", params

      return

    end if

    if ( .not. l_implemented ) then
      call setup_prognostic_variables( gr%nnzp )        ! intent(in)
    end if

    ! Both prognostic variables and diagnostic variables need to be
    ! declared, allocated, initialized, and deallocated whether CLUBB
    ! is part of a larger model or not.
    call setup_diagnostic_variables( gr%nnzp )


    return
  end subroutine setup_clubb_core

  !-----------------------------------------------------------------------
  subroutine cleanup_clubb_core( l_implemented )
    !
    !  Description:
    !    Frees memory used by the model itself.
    !
    !--------------------------------------------------------------------
    use parameters_model, only: sclrtol ! Variable

    use variables_diagnostic_module, only: & 
      cleanup_diagnostic_variables ! Procedure
    use variables_prognostic_module, only: & 
      cleanup_prognostic_variables ! Procedure

    implicit none

    ! Flag to see if CLUBB is running on it's own,
    ! or if it's implemented as part of a host model.
    logical, intent(in) :: l_implemented   ! (T/F)

    !----- Begin Code -----

    if ( .not. l_implemented ) then
      call cleanup_prognostic_variables( )
    end if

    ! Both prognostic variables and diagnostic variables need to be
    ! declared, allocated, initialized, and deallocated whether CLUBB
    ! is part of a larger model or not.
    call cleanup_diagnostic_variables( )

    ! De-allocate the array for the passive scalar tolerances
    deallocate( sclrtol )

    return
  end subroutine cleanup_clubb_core

  !-----------------------------------------------------------------------
  subroutine trapezoidal_rule &
             ( l_call_pdf_closure_twice,                    & ! intent(in)
               p_in_Pa, exner, wm_zm, wp2, wp3,             & ! intent(in)
               sigma_sqd_w, Skw_zm, rtm, rtp2,              & ! intent(in)
               wprtp, thlm, thlp2, wpthlp, rtpthlp,         & ! intent(in)
               sclrm, wpsclrp, sclrp2, sclrprtp, sclrpthlp, & ! intent(in)
               wprtp2, wp2rtp, wpthlp2, wp2thlp,            & ! intent(inout)
               wprtpthlp, cf, rcm, wp2thvp, wp2rcp,         & ! intent(inout)
               wpsclrprtp, wpsclrp2, wpsclrpthlp,           & ! intent(inout)
               wp2sclrp, pdf_params, err_code )               ! intent(inout)
    !
    ! Description:  This subroutine takes the output variables on the thermo.
    ! grid and either interpolates them to the momentum grid or calls the 
    ! subroutine pdf_closure a second time, determined by the variable 
    ! l_call_pdf_closure_twice.  It then calls the function trapezoid to 
    ! recompute the variables on the thermo. grid.
    ! ldgrant June 2009
    !-----------------------------------------------------------------------
    use grid_class, only: &
      gr, & ! Variable
      zt2zm ! Procedure

    use parameters_model, only: &
      sclr_dim ! Number of passive scalar variables

    use variables_prognostic_module, only: &
      pdf_parameter ! Derived data type

    use pdf_closure_module, only:  & 
      ! Procedure 
      pdf_closure     ! Prob. density function

    use constants, only:  & 
      fstderr  ! Variable(s)

    use error_code, only :  & 
      clubb_var_equals_NaN ! Variable(s)

    implicit none

    ! Input variables
    logical, intent(in) :: l_call_pdf_closure_twice

    real, dimension(gr%nnzp), intent(in) :: &
      p_in_Pa,     & ! Pressure                       [Pa]
      exner,       & ! Exner function.                [-]
      wm_zm,       & ! mean w on momentum levels      [m/s]
      wp2,         & ! w'^2                           [m^2/s^2]
      wp3,         & ! w'^3                           [m^3/s^3]
      sigma_sqd_w, & ! Width of individual w plumes   [-]
      Skw_zm,      & ! Skewness of w on momentum grid [-]
      rtm,         & ! Total water mixing ratio       [kg/kg]
      rtp2,        & ! rt'^2                          [kg^2/kg^2]
      wprtp,       & ! w' r_t'                        [(kg m)(kg s)]
      thlm,        & ! Mean th_l                      [K]
      thlp2,       & ! th_l'^2                        [K^2]
      wpthlp,      & ! w' th_l'                       [(m K)/s]
      rtpthlp        ! r_t' th_l'                     [(K kg)/kg]

    real, dimension(gr%nnzp, sclr_dim), intent(in) :: &
      sclrm,       & ! Mean passive scalar        [units vary]
      wpsclrp,     & ! w' sclr'                   [units vary]
      sclrp2,      & ! sclr'^2                    [units vary]
      sclrprtp,    & ! sclr' r_t'                 [units vary]
      sclrpthlp      ! sclr' th_l'                [units vary]

    ! Input/Output variables
    integer, intent(inout) :: err_code ! Are the outputs usable numbers?

    real, dimension(gr%nnzp), intent(inout) :: &
      wprtp2,      & ! w'rt'^2                   [m kg^2/kg^2]
      wp2rtp,      & ! w'^2 rt'                  [m^2 kg/kg]
      wpthlp2,     & ! w'thl'^2                  [m K^2/s]
      wp2thlp,     & ! w'^2 thl'                 [m^2 K/s^2]
      wprtpthlp,   & ! w'rt'thl'                 [m kg K/kg s]
      cf,          & ! Cloud Fraction            [%]
      rcm,         & ! Liquid water mixing ratio [kg/kg]
      wp2thvp,     & ! w'^2 th_v'                [m^2 K/s^2]
      wp2rcp         ! w'^2 rc'                  [m^2 kg/kg s^2]    

    real, dimension(gr%nnzp,sclr_dim), intent(inout) :: & 
      wpsclrprtp,  & ! w'sclr'rt' 
      wpsclrp2,    & ! w'sclr'^2
      wpsclrpthlp, & ! w'sclr'thl'
      wp2sclrp       ! w'^2 sclr'

    type (pdf_parameter), intent(inout) :: &
      pdf_params ! PDF parameters [units vary]

    ! Local variables
    integer :: i, k

    real, dimension(gr%nnzp) :: &
      wprtp2_zm,    & ! w'rt'^2 on momentum grid                   [m kg^2/kg^2]
      wp2rtp_zm,    & ! w'^2 rt' on momentum grid                  [m^2 kg/kg]
      wpthlp2_zm,   & ! w'thl'^2 on momentum grid                  [m K^2/s]
      wp2thlp_zm,   & ! w'^2 thl' on momentum grid                 [m^2 K/s^2]
      wprtpthlp_zm, & ! w'rt'thl' on momentum grid                 [m kg K/kg s]
      cf_zm,        & ! Cloud Fraction on momentum grid            [%]
      rcm_zm,       & ! Liquid water mixing ratio on momentum grid [kg/kg]
      wp2thvp_zm,   & ! w'^2 th_v' on momentum grid                [m^2 K/s^2]
      wp2rcp_zm,    & ! w'^2 rc' on momentum grid                  [m^2 kg/kg s^2]

      wp4_zm,       & ! w'^4 on momentum grid          [m^4/s^4]
      wpthvp_zm,    & ! Buoyancy flux on momentum grid [(K m)/s]
      rtpthvp_zm,   & ! r_t' th_v' on momentum grid    [(kg K)/kg]
      thlpthvp_zm,  & ! th_l' th_v' on momentum grid   [K^2]
      wprcp_zm,     & ! w' r_c' on momentum grid       [(m kg)/(s kg)]
      rtprcp_zm,    & ! r_t' r_c' on momentum grid     [(kg^2)/(kg^2)]
      thlprcp_zm,   & ! th_l' r_c' on momentum grid    [(K kg)/kg]
      rcp2_zm         ! r_c'^2 on momentum grid        [(kg^2)/(kg^2)]
      

    real, dimension(gr%nnzp,sclr_dim) :: & 
      wpsclrprtp_zm,  & ! w'sclr'rt' on momentum grid 
      wpsclrp2_zm,    & ! w'sclr'^2 on momentum grid 
      wpsclrpthlp_zm, & ! w'sclr'thl' on momentum grid 
      wp2sclrp_zm,    & ! w'^2 sclr' on momentum grid

      sclrm_zm,       & ! Passive scalar mean on momentum grid
      sclrpthvp_zm,   & ! sclr'th_v' on momentum grid
      sclrprcp_zm       ! sclr'rc' on momentum grid

    ! Components of PDF_parameters on the momentum grid (_zm) and on the thermo. grid (_zt)
    real, dimension(gr%nnzp) :: &
      w1_zt,        & ! Mean of w for 1st normal distribution                 [m/s]
      w1_zm,        & ! Mean of w for 1st normal distribution                 [m/s]
      w2_zm,        & ! Mean of w for 2nd normal distribution                 [m/s]
      w2_zt,        & ! Mean of w for 2nd normal distribution                 [m/s]
      sw1_zm,       & ! Variance of w for 1st normal distribution         [m^2/s^2]
      sw1_zt,       & ! Variance of w for 1st normal distribution         [m^2/s^2]
      sw2_zm,       & ! Variance of w for 2nd normal distribution         [m^2/s^2]
      sw2_zt,       & ! Variance of w for 2nd normal distribution         [m^2/s^2]
      rt1_zm,       & ! Mean of r_t for 1st normal distribution             [kg/kg]
      rt1_zt,       & ! Mean of r_t for 1st normal distribution             [kg/kg]
      rt2_zm,       & ! Mean of r_t for 2nd normal distribution             [kg/kg]
      rt2_zt,       & ! Mean of r_t for 2nd normal distribution             [kg/kg]
      srt1_zm,      & ! Variance of r_t for 1st normal distribution     [kg^2/kg^2]
      srt1_zt,      & ! Variance of r_t for 1st normal distribution     [kg^2/kg^2]
      srt2_zm,      & ! Variance of r_t for 2nd normal distribution     [kg^2/kg^2]
      srt2_zt,      & ! Variance of r_t for 2nd normal distribution     [kg^2/kg^2]
      crt1_zm,      & ! Coefficient for s'                                      [-]
      crt1_zt,      & ! Coefficient for s'                                      [-]
      crt2_zm,      & ! Coefficient for s'                                      [-]
      crt2_zt,      & ! Coefficient for s'                                      [-]
      cthl1_zm,     & ! Coefficient for s'                                    [1/K]
      cthl1_zt,     & ! Coefficient for s'                                    [1/K]
      cthl2_zm,     & ! Coefficient for s'                                    [1/K]
      cthl2_zt,     & ! Coefficient for s'                                    [1/K]
      thl1_zm,      & ! Mean of th_l for 1st normal distribution                [K]
      thl1_zt,      & ! Mean of th_l for 1st normal distribution                [K]
      thl2_zm,      & ! Mean of th_l for 2nd normal distribution                [K]
      thl2_zt,      & ! Mean of th_l for 2nd normal distribution                [K]
      sthl1_zm,     & ! Variance of th_l for 1st normal distribution          [K^2]
      sthl1_zt,     & ! Variance of th_l for 1st normal distribution          [K^2]
      sthl2_zm,     & ! Variance of th_l for 2nd normal distribution          [K^2]
      sthl2_zt        ! Variance of th_l for 2nd normal distribution          [K^2]

    ! Continuation of PDF_parameters above (split into two sections of variable
    ! declarations so there are not more than 39 continuation lines of code).
    real, dimension(gr%nnzp) :: &
      a_zm,           & ! Weight of 1st normal distribution (Sk_w dependent)      [-]
      a_zt,           & ! Weight of 1st normal distribution (Sk_w dependent)      [-]
      rc1_zm,         & ! Mean of r_c for 1st normal distribution             [kg/kg]
      rc1_zt,         & ! Mean of r_c for 1st normal distribution             [kg/kg]
      rc2_zm,         & ! Mean of r_c for 2nd normal distribution             [kg/kg]
      rc2_zt,         & ! Mean of r_c for 2nd normal distribution             [kg/kg]
      rsl1_zm,        & ! Mean of r_sl for 1st normal distribution            [kg/kg]
      rsl1_zt,        & ! Mean of r_sl for 1st normal distribution            [kg/kg]
      rsl2_zm,        & ! Mean of r_sl for 2nd normal distribution            [kg/kg]
      rsl2_zt,        & ! Mean of r_sl for 2nd normal distribution            [kg/kg]
      cloud_frac1_zm, & ! Cloud fraction for 1st normal distribution              [-]
      cloud_frac1_zt, & ! Cloud fraction for 1st normal distribution              [-]
      cloud_frac2_zm, & ! Cloud fraction for 2nd normal distribution              [-]
      cloud_frac2_zt, & ! Cloud fraction for 2nd normal distribution              [-]
      s1_zm,          & ! Mean of s for 1st normal distribution               [kg/kg]
      s1_zt,          & ! Mean of s for 1st normal distribution               [kg/kg]
      s2_zm,          & ! Mean of s for 2nd normal distribution               [kg/kg]
      s2_zt,          & ! Mean of s for 2nd normal distribution               [kg/kg]
      ss1_zm,         & ! Standard deviation of s for 1st normal distribution [kg/kg]
      ss1_zt,         & ! Standard deviation of s for 1st normal distribution [kg/kg]
      ss2_zm,         & ! Standard deviation of s for 2nd normal distribution [kg/kg]
      ss2_zt,         & ! Standard deviation of s for 2nd normal distribution [kg/kg]
      rrtthl_zm,      & ! Within-a-normal correlation of r_t and th_l             [-]
      rrtthl_zt,      & ! Within-a-normal correlation of r_t and th_l             [-]
      alpha_thl_zm,   & ! Factor relating to normalized variance for th_l         [-]
      alpha_thl_zt,   & ! Factor relating to normalized variance for th_l         [-]
      alpha_rt_zm,    & ! Factor relating to normalized variance for r_t          [-]
      alpha_rt_zt       ! Factor relating to normalized variance for r_t          [-]

    !----------------------- Begin Code -----------------------------

    ! Store components of pdf_params in the locally declared variables, since
    ! pdf_params may be overwritten if l_call_pdf_closure_twice is true.
    w1_zt          = pdf_params%w1
    w2_zt          = pdf_params%w2
    sw1_zt         = pdf_params%sw1
    sw2_zt         = pdf_params%sw2
    rt1_zt         = pdf_params%rt1
    rt2_zt         = pdf_params%rt2
    srt1_zt        = pdf_params%srt1
    srt2_zt        = pdf_params%srt2
    crt1_zt        = pdf_params%crt1
    crt2_zt        = pdf_params%crt2
    cthl1_zt       = pdf_params%cthl1
    cthl2_zt       = pdf_params%cthl2
    thl1_zt        = pdf_params%thl1
    thl2_zt        = pdf_params%thl2
    sthl1_zt       = pdf_params%sthl1
    sthl2_zt       = pdf_params%sthl2
    a_zt           = pdf_params%a
    rc1_zt         = pdf_params%rc1
    rc2_zt         = pdf_params%rc2
    rsl1_zt        = pdf_params%rsl1
    rsl2_zt        = pdf_params%rsl2
    cloud_frac1_zt = pdf_params%cloud_frac1
    cloud_frac2_zt = pdf_params%cloud_frac2
    s1_zt          = pdf_params%s1
    s2_zt          = pdf_params%s2
    ss1_zt         = pdf_params%ss1
    ss2_zt         = pdf_params%ss2
    rrtthl_zt      = pdf_params%rrtthl
    alpha_thl_zt   = pdf_params%alpha_thl
    alpha_rt_zt    = pdf_params%alpha_rt

    ! If l_call_pdf_closure_twice is true, use the subroutine pdf_closure to compute
    ! the variables on the momentum grid
    if (l_call_pdf_closure_twice) then

      ! Interpolate sclrm to the momentum level for use in pdf_closure
      do i = 1, sclr_dim
        sclrm_zm(:,i) = zt2zm( sclrm(:,i) )
      end do ! i = 1, sclr_dim

      ! Call pdf_closure to compute the variables on the momentum grid for use
      ! in the trapezoid rule
      do k = 1, gr%nnzp, 1
        call pdf_closure & 
           ( zt2zm( p_in_Pa, k ), zt2zm( exner, k ), wm_zm(k),                      & ! intent(in)
               wp2(k), zt2zm( wp3, k ), sigma_sqd_w(k),                             & ! intent(in)
             Skw_zm(k), zt2zm( rtm, k ), rtp2(k), wprtp(k),                         & ! intent(in)
             zt2zm( thlm, k ), thlp2(k), wpthlp(k),                                 & ! intent(in)
             rtpthlp(k), sclrm_zm(k,:), wpsclrp(k,:),                               & ! intent(in)
             sclrp2(k,:), sclrprtp(k,:), sclrpthlp(k,:), k,                         & ! intent(in)
             wp4_zm(k), wprtp2_zm(k), wp2rtp_zm(k),                                 & ! intent(out)
             wpthlp2_zm(k), wp2thlp_zm(k), wprtpthlp_zm(k),                         & ! intent(out)
             cf_zm(k), rcm_zm(k), wpthvp_zm(k), wp2thvp_zm(k), rtpthvp_zm(k),       & ! intent(out)
             thlpthvp_zm(k), wprcp_zm(k), wp2rcp_zm(k), rtprcp_zm(k), thlprcp_zm(k),& ! intent(out)
             rcp2_zm(k), pdf_params,                                                & ! intent(out)
             err_code,                                                              & ! intent(out)
             wpsclrprtp_zm(k,:), wpsclrp2_zm(k,:), sclrpthvp_zm(k,:),               & ! intent(out)
             wpsclrpthlp_zm(k,:), sclrprcp_zm(k,:), wp2sclrp_zm(k,:) )                ! intent(out)

        ! Subroutine may produce NaN values, and if so, exit
        ! gracefully.
        ! Joshua Fasching March 2008

        if ( err_code == clubb_var_equals_NaN ) then
          write(fstderr,*) "At grid level = ",k
          return
        end if
      end do ! k = 2, gr%nnzp-1

      ! Store, in locally declared variables, the pdf_params output 
      ! from the second call to pdf_closure
      w1_zm          = pdf_params%w1
      w2_zm          = pdf_params%w2
      sw1_zm         = pdf_params%sw1
      sw2_zm         = pdf_params%sw2
      rt1_zm         = pdf_params%rt1
      rt2_zm         = pdf_params%rt2
      srt1_zm        = pdf_params%srt1
      srt2_zm        = pdf_params%srt2
      crt1_zm        = pdf_params%crt1
      crt2_zm        = pdf_params%crt2
      cthl1_zm       = pdf_params%cthl1
      cthl2_zm       = pdf_params%cthl2
      thl1_zm        = pdf_params%thl1
      thl2_zm        = pdf_params%thl2
      sthl1_zm       = pdf_params%sthl1
      sthl2_zm       = pdf_params%sthl2
      a_zm           = pdf_params%a
      rc1_zm         = pdf_params%rc1
      rc2_zm         = pdf_params%rc2
      rsl1_zm        = pdf_params%rsl1
      rsl2_zm        = pdf_params%rsl2
      cloud_frac1_zm = pdf_params%cloud_frac1
      cloud_frac2_zm = pdf_params%cloud_frac2
      s1_zm          = pdf_params%s1
      s2_zm          = pdf_params%s2
      ss1_zm         = pdf_params%ss1
      ss2_zm         = pdf_params%ss2
      rrtthl_zm      = pdf_params%rrtthl
      alpha_thl_zm   = pdf_params%alpha_thl
      alpha_rt_zm    = pdf_params%alpha_rt

    else      

      ! Interpolate thermodynamic variables to the momentum grid.
      ! Since top momentum level is higher than top thermo. level,
      ! set variables at top momentum level to 0.
      wprtp2_zm             = zt2zm( wprtp2 )
      wprtp2_zm(gr%nnzp)    = 0.0
      wp2rtp_zm             = zt2zm( wp2rtp )
      wp2rtp_zm(gr%nnzp)    = 0.0
      wpthlp2_zm            = zt2zm( wpthlp2 )
      wpthlp2_zm(gr%nnzp)   = 0.0
      wp2thlp_zm            = zt2zm( wp2thlp )
      wp2thlp_zm(gr%nnzp)   = 0.0
      wprtpthlp_zm          = zt2zm( wprtpthlp )
      wprtpthlp_zm(gr%nnzp) = 0.0
      cf_zm                 = zt2zm( cf )
      cf_zm(gr%nnzp)        = 0.0
      rcm_zm                = zt2zm( rcm )
      rcm_zm(gr%nnzp)       = 0.0
      wp2thvp_zm            = zt2zm( wp2thvp )
      wp2thvp_zm(gr%nnzp)   = 0.0
      wp2rcp_zm             = zt2zm( wp2rcp )
      wp2rcp_zm(gr%nnzp)    = 0.0

      do i = 1, sclr_dim
        wpsclrprtp_zm(:,i)        = zt2zm( wpsclrprtp(:,i) )
        wpsclrprtp_zm(gr%nnzp,i)  = 0.0
        wpsclrp2_zm(:,i)          = zt2zm( wpsclrp2(:,i) )
        wpsclrp2_zm(gr%nnzp,i)    = 0.0
        wpsclrpthlp_zm(:,i)       = zt2zm( wpsclrpthlp(:,i) )
        wpsclrpthlp_zm(gr%nnzp,i) = 0.0
        wp2sclrp_zm(:,i)          = zt2zm( wp2sclrp(:,i) )
        wp2sclrp_zm(gr%nnzp,i)    = 0.0
      end do ! i = 1, sclr_dim

      w1_zm                   = zt2zm( pdf_params%w1 )
      w1_zm(gr%nnzp)          = 0.0
      w2_zm                   = zt2zm( pdf_params%w2 )
      w2_zm(gr%nnzp)          = 0.0
      sw1_zm                  = zt2zm( pdf_params%sw1 )
      sw1_zm(gr%nnzp)         = 0.0   
      sw2_zm                  = zt2zm( pdf_params%sw2 )
      sw2_zm(gr%nnzp)         = 0.0
      rt1_zm                  = zt2zm( pdf_params%rt1 )
      rt1_zm(gr%nnzp)         = 0.0
      rt2_zm                  = zt2zm( pdf_params%rt2 )
      rt2_zm(gr%nnzp)         = 0.0
      srt1_zm                 = zt2zm( pdf_params%srt1 )
      srt1_zm(gr%nnzp)        = 0.0
      srt2_zm                 = zt2zm( pdf_params%srt2 )
      srt2_zm(gr%nnzp)        = 0.0
      crt1_zm                 = zt2zm( pdf_params%crt1 )
      crt1_zm(gr%nnzp)        = 0.0
      crt2_zm                 = zt2zm( pdf_params%crt2 )
      crt2_zm(gr%nnzp)        = 0.0
      cthl1_zm                = zt2zm( pdf_params%cthl1 )
      cthl1_zm(gr%nnzp)       = 0.0
      cthl2_zm                = zt2zm( pdf_params%cthl2 )
      cthl2_zm(gr%nnzp)       = 0.0
      thl1_zm                 = zt2zm( pdf_params%thl1 )
      thl1_zm(gr%nnzp)        = 0.0
      thl2_zm                 = zt2zm( pdf_params%thl2 )
      thl2_zm(gr%nnzp)        = 0.0
      sthl1_zm                = zt2zm( pdf_params%sthl1 )
      sthl1_zm(gr%nnzp)       = 0.0
      sthl2_zm                = zt2zm( pdf_params%sthl2 )
      sthl2_zm(gr%nnzp)       = 0.0
      a_zm                    = zt2zm( pdf_params%a )
      a_zm(gr%nnzp)           = 0.0
      rc1_zm                  = zt2zm( pdf_params%rc1 )
      rc1_zm(gr%nnzp)         = 0.0
      rc2_zm                  = zt2zm( pdf_params%rc2 )
      rc2_zm(gr%nnzp)         = 0.0
      rsl1_zm                 = zt2zm( pdf_params%rsl1 )
      rsl1_zm(gr%nnzp)        = 0.0
      rsl2_zm                 = zt2zm( pdf_params%rsl2 )
      rsl2_zm(gr%nnzp)        = 0.0
      cloud_frac1_zm          = zt2zm( pdf_params%cloud_frac1 )
      cloud_frac1_zm(gr%nnzp) = 0.0
      cloud_frac2_zm          = zt2zm( pdf_params%cloud_frac2 )
      cloud_frac2_zm(gr%nnzp) = 0.0
      s1_zm                   = zt2zm( pdf_params%s1 )
      s1_zm(gr%nnzp)          = 0.0
      s2_zm                   = zt2zm( pdf_params%s2 )
      s2_zm(gr%nnzp)          = 0.0
      ss1_zm                  = zt2zm( pdf_params%ss1 )
      ss1_zm(gr%nnzp)         = 0.0
      ss2_zm                  = zt2zm( pdf_params%ss2 )
      ss2_zm(gr%nnzp)         = 0.0
      rrtthl_zm               = zt2zm( pdf_params%rrtthl )
      rrtthl_zm(gr%nnzp)      = 0.0
      alpha_thl_zm            = zt2zm( pdf_params%alpha_thl )
      alpha_thl_zm(gr%nnzp)   = 0.0
      alpha_rt_zm             = zt2zm( pdf_params%alpha_rt )
      alpha_rt_zm(gr%nnzp)    = 0.0
    end if ! l_call_pdf_closure_twice

    ! Use the trapezoidal rule to recompute the variables on the zt level
    wprtp2    = trapezoid( wprtp2, wprtp2_zm )
    wp2rtp    = trapezoid( wp2rtp, wp2rtp_zm )
    wpthlp2   = trapezoid( wpthlp2, wpthlp2_zm )
    wp2thlp   = trapezoid( wp2thlp, wp2thlp_zm )
    wprtpthlp = trapezoid( wprtpthlp, wprtpthlp_zm )
    cf        = trapezoid( cf, cf_zm )
    rcm       = trapezoid( rcm, rcm_zm )
    wp2thvp   = trapezoid( wp2thvp, wp2thvp_zm )
    wp2rcp    = trapezoid( wp2rcp, wp2rcp_zm ) 

    do i = 1, sclr_dim 
      wpsclrprtp(:,i)  = trapezoid( wpsclrprtp(:,i), wpsclrprtp_zm(:,i) )
      wpsclrp2(:,i)    = trapezoid( wpsclrp2(:,i), wpsclrp2_zm(:,i) )
      wpsclrpthlp(:,i) = trapezoid( wpsclrpthlp(:,i), wpsclrpthlp_zm(:,i) )
      wp2sclrp(:,i)    = trapezoid( wp2sclrp(:,i), wp2sclrp_zm(:,i) )
    end do ! i = 1, sclr_dim

    pdf_params%w1          = trapezoid( w1_zt, w1_zm )
    pdf_params%w2          = trapezoid( w2_zt, w2_zm )
    pdf_params%sw1         = trapezoid( sw1_zt, sw1_zm )
    pdf_params%sw2         = trapezoid( sw2_zt, sw2_zm )
    pdf_params%rt1         = trapezoid( rt1_zt, rt1_zm )
    pdf_params%rt2         = trapezoid( rt2_zt, rt2_zm )
    pdf_params%srt1        = trapezoid( srt1_zt, srt1_zm )
    pdf_params%srt2        = trapezoid( srt2_zt, srt2_zm )
    pdf_params%crt1        = trapezoid( crt1_zt, crt1_zm )
    pdf_params%crt2        = trapezoid( crt2_zt, crt2_zm )
    pdf_params%cthl1       = trapezoid( cthl1_zt, cthl1_zm )
    pdf_params%cthl2       = trapezoid( cthl2_zt, cthl2_zm )
    pdf_params%thl1        = trapezoid( thl1_zt, thl1_zm )
    pdf_params%thl2        = trapezoid( thl2_zt, thl2_zm )
    pdf_params%sthl1       = trapezoid( sthl1_zt, sthl1_zm )
    pdf_params%sthl2       = trapezoid( sthl2_zt, sthl2_zm )
    pdf_params%a           = trapezoid( a_zt, a_zm )
    pdf_params%rc1         = trapezoid( rc1_zt, rc1_zm )
    pdf_params%rc2         = trapezoid( rc2_zt, rc2_zm )
    pdf_params%rsl1        = trapezoid( rsl1_zt, rsl1_zm )
    pdf_params%rsl2        = trapezoid( rsl2_zt, rsl2_zm )
    pdf_params%cloud_frac1 = trapezoid( cloud_frac1_zt, cloud_frac1_zm )
    pdf_params%cloud_frac2 = trapezoid( cloud_frac2_zt, cloud_frac2_zm )
    pdf_params%s1          = trapezoid( s1_zt, s1_zm )
    pdf_params%s2          = trapezoid( s2_zt, s2_zm )
    pdf_params%ss1         = trapezoid( ss1_zt, ss1_zm )
    pdf_params%ss2         = trapezoid( ss2_zt, ss2_zm )
    pdf_params%rrtthl      = trapezoid( rrtthl_zt, rrtthl_zm )
    pdf_params%alpha_thl   = trapezoid( alpha_thl_zt, alpha_thl_zm )
    pdf_params%alpha_rt    = trapezoid( alpha_rt_zt, alpha_rt_zm )
    ! End of trapezoidal rule

    return
  end subroutine trapezoidal_rule

  !-----------------------------------------------------------------------
  pure function trapezoid( variable_zt, variable_zm )
    !
    ! Description: Function which uses the trapezoidal rule from calculus
    ! to recompute the values for the variables on the thermo. grid which
    ! are output from the first call to pdf_closure in module clubb_core.
    ! ldgrant June 2009
    !--------------------------------------------------------------------

    use grid_class, only: gr ! Variable

    implicit none

    ! Input Variables
    real, dimension(gr%nnzp), intent(in) :: &
      variable_zt, &
      variable_zm

    ! Result
    real, dimension(gr%nnzp) :: trapezoid

    ! Local Variable
    integer :: k
 
    !------------ Begin Code --------------

    trapezoid(1) = variable_zt(1)

    do k = 2, gr%nnzp
      trapezoid(k) =  0.5 * ( variable_zm(k) + variable_zt(k) ) &
                          * ( gr%zm(k) - gr%zt(k) ) * gr%dzt(k) &
                    + 0.5 * ( variable_zt(k) + variable_zm(k-1) ) &
                          * ( gr%zt(k) - gr%zm(k-1) ) * gr%dzt(k)
    end do ! k = 2, gr%nnzp

    return 
  end function trapezoid

  !-----------------------------------------------------------------------
  subroutine compute_cloud_cover &
           ( pdf_params, cf, rcm,        & ! intent(in)
             cloud_cover, rcm_in_layer )   ! intent(out)
    !
    ! Description:  Subroutine to compute cloud cover (the amount of sky
    ! covered by cloud) and rcm in layer (liquid water mixing ratio in
    ! the portion of the grid box filled by cloud).
    ! ldgrant July 2009 
    !---------------------------------------------------------------------

    use constants, only: rc_tol ! Variable

    use grid_class, only: gr ! Variable
  
    use variables_prognostic_module, only: &
      pdf_parameter ! Derived data type

    implicit none

    ! External functions
    intrinsic :: abs, min, max

    ! Input variables
    real, dimension(gr%nnzp), intent(in) :: &
      cf,  & ! Cloud fraction             [%]
      rcm    ! Liquid water mixing ratio  [kg/kg]

    type (pdf_parameter), intent(in) :: &
      pdf_params ! PDF Parameters  [units vary]

    ! Output variables
    real, dimension(gr%nnzp), intent(out) :: &
      cloud_cover,  & ! Cloud cover                               [%]
      rcm_in_layer    ! Liquid water mixing ratio in cloud layer  [kg/kg]

    ! Local variables
    real, dimension(gr%nnzp) :: &
      s_mean,                & ! Mean of the two Gaussian distributions
      vert_cloud_frac_upper, & ! Fraction of cloud in top half of grid box
      vert_cloud_frac_lower, & ! Fraction of cloud in bottom half of grid box
      vert_cloud_frac          ! Fraction of cloud filling the grid box in the vertical

    integer :: k

    ! ------------ Begin code ---------------

    do k = 1, gr%nnzp

      s_mean(k) = pdf_params%a(k) * pdf_params%s1(k) + &
                  (1.0-pdf_params%a(k)) * pdf_params%s2(k)

    end do

    do k = 2, gr%nnzp-1, 1

      if ( rcm(k) < rc_tol ) then ! No cloud at this level

        cloud_cover(k)  = cf(k)
        rcm_in_layer(k) = rcm(k)
  
      else if ( ( rcm(k+1) > rc_tol ) .and. ( rcm(k-1) > rc_tol ) ) then
        ! There is cloud above and below, 
        !   so assume cloud fills grid box from top to bottom
  
        cloud_cover(k) = cf(k)
        rcm_in_layer(k) = rcm(k)

      else if ( ( rcm(k+1) < rc_tol ) .or. ( rcm(k-1) < rc_tol) ) then 
        ! Cloud may fail to reach gridbox top or base or both

        ! First let the cloud fill the entire grid box, then overwrite
        ! vert_cloud_frac_upper(k) and/or vert_cloud_frac_lower(k)
        ! for a cloud top, cloud base, or one-point cloud.
        vert_cloud_frac_upper(k) = 0.5
        vert_cloud_frac_lower(k) = 0.5

        if ( rcm(k+1) < rc_tol ) then ! Cloud top

          vert_cloud_frac_upper(k) = &
                   ( ( 0.5 / gr%dzm(k) ) / ( gr%zm(k) - gr%zt(k) ) ) * &
                   ( rcm(k) / ( rcm(k) + abs( s_mean(k+1) ) ) ) 

          vert_cloud_frac_upper(k) = min( 0.5, vert_cloud_frac_upper(k) ) 

        else

          vert_cloud_frac_upper(k) = 0.5

        end if

        if ( rcm(k-1) < rc_tol ) then ! Cloud base

          vert_cloud_frac_lower(k) = &
                   ( ( 0.5 / gr%dzm(k-1) ) / ( gr%zt(k) - gr%zm(k-1) ) ) * &
                   ( rcm(k) / ( rcm(k) + abs( s_mean(k-1) ) ) )

          vert_cloud_frac_lower(k) = min( 0.5, vert_cloud_frac_lower(k) )

        else

          vert_cloud_frac_lower(k) = 0.5

        end if

        vert_cloud_frac(k) = &
          vert_cloud_frac_upper(k) + vert_cloud_frac_lower(k)

        vert_cloud_frac(k) = &
          max( cf(k), min( 1.0, vert_cloud_frac(k) ) )

        cloud_cover(k)  = cf(k) / vert_cloud_frac(k)
        rcm_in_layer(k) = rcm(k) / vert_cloud_frac(k)

      else

        print*, "Error: Should not arrive here in computation of cloud_cover"
        stop

      end if

    end do ! k = 2, gr%nnzp-1, 1

    cloud_cover(1)       = cf(1)
    cloud_cover(gr%nnzp) = cf(gr%nnzp)

    rcm_in_layer(1)       = rcm(1)
    rcm_in_layer(gr%nnzp) = rcm(gr%nnzp)

    return
  end subroutine compute_cloud_cover
  !-----------------------------------------------------------------------

end module clubb_core
