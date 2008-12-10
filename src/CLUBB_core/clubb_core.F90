!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------
module clubb_core

! Description:
!   The module containing the `core' of the HOC model.

! References:
!   None
!-----------------------------------------------------------------------

  implicit none

  public ::  & 
    setup_clubb_core, & 
    advance_clubb_core, & 
    cleanup_clubb_core

  private :: latin_hypercube_sampling

  private ! Default Scope

  contains

  !-----------------------------------------------------------------------
  subroutine advance_clubb_core & 
             ( iter, l_implemented, dt, fcor, & 
               thlm_forcing, rtm_forcing, um_forcing, vm_forcing, wm_zm, wm_zt, & 
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, & 
               p_in_Pa, rho_zm, rho, exner, & 
               wpsclrp_sfc, wpedsclrp_sfc,    & ! Optional
               um, vm, upwp, vpwp, up2, vp2, & 
               thlm, rtm, wprtp, wpthlp, wp2, wp3, & 
               rtp2, thlp2, rtpthlp, & 
               sigma_sqd_w, tau_zm, rcm, cf, & 
               err_code,  & 
               sclrm, sclrm_forcing, edsclrm,  & ! Optional
               wpsclrp )

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
      sclr_dim

    use model_flags, only: & 
      l_LH_on,  & ! Variable(s)
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
      wpthvp, & 
      thlpthvp, & 
      rtpthvp, & 
      wprcp, & 
      rtprcp, & 
      thlprcp, & 
      pdf_parms, & 
      rcp2, & 
      rsat, & 
      shear, & 
      Kh_zt, & 
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
      wpedsclrp, & 
      sclrpthvp,   & ! sclr'th_v'
      sclrprtp,    & ! sclr'rt'
      sclrp2,      & ! sclr'^2
      sclrpthlp,   & ! sclr'th_l'
      sclrprcp,    & ! sclr'rc'
      wp2sclrp,    & ! w'^2 sclr'
      wpsclrp2,    & ! w'sclr'^2
      wpsclrprtp,  & ! w'sclr'rt'
      wpsclrpthlp    ! w'sclr'thl'

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
      lapack_error,  & ! Procedure(s)
      clubb_at_least_debug_level

    use Skw_module, only:  & 
      Skw_func ! Procedure

    use clip_explicit, only: & 
      clip_covariances_denom ! Procedure(s)

    use permute_height_time_mod, only:  & 
      permute_height_time ! Procedure

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

    intrinsic :: sqrt, min, max, exp, mod

    ! Input
    integer, intent(in) ::  & 
      iter      ! Closure iteration number

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

    ! Optional Input Variables
    real, intent(in),  dimension(sclr_dim) ::  & 
      wpsclrp_sfc,   & ! Scalar flux at surface           [units m/s]
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

    ! Optional Input/Output Variables
    real, intent(inout), dimension(gr%nnzp,sclr_dim) :: & 
      sclrm,         & ! Passive scalar mean.           [units vary]
      sclrm_forcing, & ! Passive scalar forcing.        [units vary/s]
      edsclrm,       & ! Eddy passive scalar mean.      [units vary]
      wpsclrp          ! w'sclr'                        [units vary m/s]


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

    !------- Local variables for Latin Hypercube sampling ------------------

    integer i_rmd

    ! Number of variables to sample
    integer, parameter :: d_variables = 5

    ! n = number of calls to microphysics per timestep (normally=2)
    integer, parameter :: n_micro_call = 12

    ! sequence_length = nt/n = number of timesteps before sequence repeats.
    integer, parameter :: sequence_length = 1

    ! nt = number of random samples before sequence of repeats (normally=10)
    integer, parameter :: nt_repeat = n_micro_call * sequence_length

    ! A true/false flag that determines whether
    !     the PDF allows us to construct a sample
    !       logical sample_flag

    integer, dimension(gr%nnzp, nt_repeat, d_variables+1)  & 
    :: p_height_time ! matrix of rand ints

    ! coeffs of s from pdf_closure
    real :: crt1, crt2, cthl1, cthl2
    !-------- End Latin hypercube section ----------------------------------

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
           wpthlp, wp2, wp3, sigma_sqd_w, rtp2, thlp2,               & ! intent(in)
           rtpthlp, tau_zm, rcm, cf, "beginning of ",        & ! intent(in)
           wpsclrp_sfc, wpedsclrp_sfc,                       & ! intent(in)
           sclrm, sclrm_forcing, edsclrm )                     ! intent(in)
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
      wpedsclrp(1,1:sclr_dim) = wpedsclrp_sfc(1:sclr_dim)
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

    ! Latin hypercube sample generation
    ! Generate p_height_time, an nnzp x nt_repeat x d_variables array of random integers
    if ( l_LH_on ) then
      i_rmd = mod( iter-1, sequence_length )
      if ( i_rmd == 0) then
        call permute_height_time( gr%nnzp, nt_repeat, d_variables+1, & ! intent(in)
                                  p_height_time )                   ! intent(out)
      end if
    end if
    ! End Latin hypercube sample generation

    ! print*, 'hoc.F: i_rmd=', i_rmd

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

    do k = 2, gr%nnzp, 1
      call pdf_closure & 
         ( p_in_Pa(k), exner(k), wm_zt(k), wp2_zt(k), wp3(k), sigma_sqd_w_zt(k), & ! intent(in)
           Skw_zt(k), rtm(k), rtp2_zt(k), zm2zt( wprtp, k ),       & ! intent(in)
           thlm(k), thlp2_zt(k), zm2zt( wpthlp, k ),               & ! intent(in)
           rtpthlp_zt(k), sclrm(k,:), sclr_tmp1(k,:),              & ! intent(in)
           sclr_tmp3(k,:),sclr_tmp2(k,:), sclr_tmp4(k,:),          & ! intent(in)
           wp4(k), wprtp2(k), wp2rtp(k),                           & ! intent(out)
           wpthlp2(k), wp2thlp(k), wprtpthlp(k),                   & ! intent(out)
           cf(k), rcm(k), wpthvp(k), wp2thvp(k), rtpthvp(k),       & ! intent(out)
           thlpthvp(k), wprcp(k), wp2rcp(k), rtprcp(k), thlprcp(k),& ! intent(out)
           rcp2(k), pdf_parms(k, :), crt1,                         & ! intent(out)
           crt2, cthl1, cthl2, err_code,                           & ! intent(out)
           wpsclrprtp(k,:), wpsclrp2(k,:), sclrpthvp(k,:),         & ! intent(out)
           wpsclrpthlp(k,:), sclrprcp(k,:), wp2sclrp(k,:) )          ! intent(out)

      ! Subroutine may produce NaN values, and if so, exit
      ! gracefully.
      ! Joshua Fasching March 2008

      if ( err_code == clubb_var_equals_NaN ) then
        write(0,*) "At grid level = ",k
        return
      end if

      !--------------------------------------------------------------
      ! Latin hypercube sampling
      !--------------------------------------------------------------
      if ( l_LH_on ) then
        ! call latin_hypercube_sampling
        !           ( k, n_micro_call, d_variables, nt_repeat, i_rmd, &
        !             crt1, crt2, cthl1, cthl2, hydromet(:,1), &
        !             cf, gr%nnzp, sample_flag, p_height_time ) &
      end if

    end do ! k = 2, nz-1

    !            print*, 'hoc.F: AKm=', AKm
    !            print*, 'hoc.F: AKm_est=', AKm_est

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
    enddo ! i=1, sclr_dim

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
    call advance_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, wp3,   & ! intent(in)
                          Kh_zt, tau_zm, Skw_zm, rtpthvp,    & ! intent(in)
                          rtm_forcing, thlpthvp,             & ! intent(in)
                          thlm_forcing, rtp2, thlp2, wp2_zt, & ! intent(in)
                          l_implemented,                     & ! intent(in)
                          sclrpthvp, sclrm_forcing, sclrp2,  & ! intent(in)
                          rtm, wprtp, thlm, wpthlp,          & ! intent(inout)
                          err_code,                          & ! intent(inout)
                          sclrm, wpsclrp )                     ! intent(inout)

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
           tau_zm, tau_zt, Skw_zm, Skw_zt, pdf_parms(:,13), & ! intent(in)
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

    call advance_windm_edsclrm( dt, wm_zt, Kh_zm, ug, vg, um_ref, vm_ref,  &
                                wp2, up2, vp2, um_forcing, vm_forcing, &
                                upwp_sfc, vpwp_sfc, wpedsclrp_sfc, fcor,  &
                                l_implemented, um, vm, edsclrm,  &
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
           rtm, wprtp, wpthlp, wp2, wp3, rtp2, thlp2, rtpthlp, & ! intent(in)
           p_in_Pa, exner, rho, rho_zm,                        & ! intent(in)
           wm_zt, sigma_sqd_w, tau_zm, rcm, cf,                        & ! intent(in)
           sclrm, edsclrm, sclrm_forcing, wpsclrp )              ! intent(in)


    if ( clubb_at_least_debug_level( 2 ) ) then
      call parameterization_check & 
           ( thlm_forcing, rtm_forcing, wm_zm, wm_zt, p_in_Pa, rho_zm, & ! intent(in)
             rho, exner, wpthlp_sfc, wprtp_sfc,                  & ! intent(in)
             upwp_sfc, vpwp_sfc, um, upwp, vm, vpwp,             & ! intent(in)
             up2, vp2, rtm, wprtp, thlm,                         & ! intent(in)
             wpthlp, wp2, wp3, sigma_sqd_w, rtp2, thlp2,         & ! intent(in)
             rtpthlp, tau_zm, rcm, cf, "end of ",                & ! intent(in)
             wpsclrp_sfc, wpedsclrp_sfc,                         & ! intent(in)
             sclrm, sclrm_forcing, edsclrm )                       ! intent(in)
    end if

!-----------------------------------------------------------------------

    return
  end subroutine advance_clubb_core


  !-----------------------------------------------------------------------
  subroutine latin_hypercube_sampling & 
             ( k, n, dvar, nt, i_rmd, & 
               crt1, crt2, cthl1, cthl2, & 
               rrainm, cf, grid, l_sflag, p_height_time )
    !       Description:
    !       Estimate using Latin Hypercubes.  This is usually disabled by default.
    !       The actual generation of a random matrix is done in a call from the
    !       subroutine hoc_initialize to permute_height_time()
    !       References:
    !-----------------------------------------------------------------------

    use variables_diagnostic_module, only:  & 
        pdf_parms,  & ! Variable(s) 
        AKm_est,  & 
        AKm, & 
        AKstd, & 
        AKstd_cld, & 
        AKm_rcm, & 
        AKm_rcc, & 
        rcm_est

    use lh_sampler_mod, only: & 
        lh_sampler ! Procedure

    use micro_calcs_mod, only: & 
        micro_calcs ! Procedure

    implicit none


    ! Input Variables
    integer, intent(in) :: k  ! index
    integer, intent(in) :: n, dvar, i_rmd, nt, grid
    logical, intent(out) :: l_sflag

    ! coeffs of s from pdf_closure
    real, intent(in) :: crt1, crt2, cthl1, cthl2

    real, dimension(grid), intent(in) ::  & 
    rrainm,  & ! Rain water mixing ratio  [kg/kg]
    cf   ! Cloud fraction           [%]

    integer, dimension(1:grid, 1:nt, 1:(dvar+1) ), intent(in) :: & 
    p_height_time ! matrix of rand ints

    ! Local Variables

    integer :: p_matrix(n, dvar+1)
    ! Sample drawn from uniform distribution
    double precision, dimension(1:n,1:(dvar+1)) :: X_u

    ! Sample that is transformed ultimately to normal-lognormal
    double precision, dimension(1:n,1:dvar) :: X_nl

    ! Choose which rows of LH sample to feed into closure.
    p_matrix(1:n,1:(dvar+1)) = & 
    p_height_time( k,n*i_rmd+1:n*i_rmd+n, 1:(dvar+1) )

    ! print*, 'hoc.F: got past p_matrix'

    ! Generate LH sample, represented by X_u and X_nl, for level k
    call lh_sampler( n, nt, dvar, p_matrix,       & ! intent(in)
                     cf(k), pdf_parms(k, :),      & ! intent(in)
                     crt1, crt2, cthl1, cthl2,    & ! intent(in)
                     rrainm(k),                   & ! intent(in)
                     X_u, X_nl, l_sflag )           ! intent(out)

    !print *, 'hoc.F: got past lh_sampler'

    ! Perform LH and analytic microphysical calculations
    call micro_calcs( n, dvar, X_u, X_nl, l_sflag,                & ! intent(in)
                      pdf_parms(k,:),                             & ! intent(in)
                      AKm_est(k), AKm(k), AKstd(k), AKstd_cld(k), & ! intent(out)
                      AKm_rcm(k), AKm_rcc(k), rcm_est(k) )          ! intent(out)

    !print*, 'k, AKm_est=', k, AKm_est(k)
    !print*, 'k, AKm=', k, AKm(k)

    return
  end subroutine latin_hypercube_sampling

  !-----------------------------------------------------------------------
  subroutine setup_clubb_core & 
             ( nzmax, T0_in, ts_nudge_in, hydromet_dim_in,  & 
               sclr_dim_in, sclrtol_in, params,  & 
               l_bugsrad, l_surface_scheme, l_kk_rain, l_icedfs, l_coamps_micro, & 
               l_cloud_sed, l_uv_nudge, l_tke_aniso,  & 
               l_implemented, grid_type, deltaz, zm_init, & 
               momentum_heights, thermodynamic_heights,  & 
               host_dx, host_dy, err_code )

    use grid_class, only: & 
      gridsetup ! Procedure

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
    ! grid spacing and momentum-level starting altitude as input.
    real, intent(in) :: & 
      deltaz,   & ! Change in altitude per level           [m]
      zm_init     ! Initial grid altitude (momentum level) [m]

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
      sclr_dim_in         ! Number of passive scalars

    real, intent(in), dimension(sclr_dim_in) :: & 
      sclrtol_in    ! Thresholds for passive scalars

    real, intent(in), dimension(nparams) :: & 
      params  ! Including C1, nu1, nu2, etc.

    ! Flags 
    logical, intent(in) ::  & 
      l_bugsrad,      & ! BUGSrad interactive radiation scheme
      l_surface_scheme,& ! Simple surface scheme
      l_kk_rain,      & ! K & K rain microphysics
      l_icedfs,       & ! Simplified ice scheme
      l_coamps_micro, & ! COAMPS microphysics scheme
      l_cloud_sed,    & ! Cloud water droplet sedimentation
      l_uv_nudge,     & ! Wind nudging
      l_tke_aniso       ! For anisotropic turbulent kinetic energy,
                        !   i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)

    ! Output variables
    integer, intent(out) :: & 
      err_code   ! Diagnostic for a problem with the setup

    ! Local variables
    real :: Lscale_max

    !----- Begin Code -----

    ! Setup flags

    call setup_model_flags & 
         ( l_bugsrad, l_surface_scheme, l_kk_rain, l_cloud_sed,     & ! intent(in)
           l_icedfs, l_coamps_micro, l_uv_nudge,  & ! intent(in)
           l_tke_aniso )                            ! intent(in)


    ! Determine the maximum allowable value for Lscale (in meters).
    if ( l_implemented ) then
      Lscale_max = 0.25 * min( host_dx, host_dy )
    else
      Lscale_max = 1.0e5
    end if

    ! Define model constant parameters
    call setup_parameters_model( T0_in, ts_nudge_in, hydromet_dim_in, & ! in
                                 sclr_dim_in, sclrtol_in, Lscale_max )

    ! Define tunable constant parameters
    call setup_parameters & 
         ( deltaz, params, nzmax, l_implemented, & ! intent(in)
           grid_type, momentum_heights, thermodynamic_heights, & ! intent(in)
           err_code )                                            ! intent(out)

    ! Error Report
    ! Joshua Fasching February 2008
    if ( err_code == clubb_var_out_of_bounds ) then

      write(fstderr,*) "Error in setup_clubb_core"

      write(fstderr,*) "Intent(in)"

      write(fstderr,*) "deltaz = ", deltaz
      write(fstderr,*) "zm_init = ", zm_init
      write(fstderr,*) "momentum_heights = ", momentum_heights
      write(fstderr,*) "thermodynamic_heights = ",  & 
        thermodynamic_heights
      write(fstderr,*) "T0_in = ", T0_in
      write(fstderr,*) "ts_nudge_in = ", ts_nudge_in
      write(fstderr,*) "params = ", params

      return

    end if

    ! if ( .not. l_implemented ) then
    !   call setup_diagnostic_variables( nzmax )
    ! end if

    ! Both prognostic variables and diagnostic variables need to be
    ! declared, allocated, initialized, and deallocated whether HOC
    ! is part of a larger model or not.
    call setup_prognostic_variables( nzmax )        ! intent(in)
    call setup_diagnostic_variables( nzmax )        ! intent(in)

    ! Setup grid
    call gridsetup( nzmax, l_implemented, grid_type,    & ! intent(in)
                    deltaz, zm_init, momentum_heights,  & ! intent(in)
                    thermodynamic_heights )               ! intent(in)

    return
  end subroutine setup_clubb_core

  !-----------------------------------------------------------------------
  subroutine cleanup_clubb_core( )

    use parameters_model, only: sclrtol ! Variable

    use variables_diagnostic_module, only: & 
      cleanup_diagnostic_variables ! Procedure
    use variables_prognostic_module, only: & 
      cleanup_prognostic_variables ! Procedure

    implicit none

    !----- Begin Code -----

    !if ( .not. l_implemented ) then
    !  call cleanup_diagnostic_variables( )
    !end if

    ! Both prognostic variables and diagnostic variables need to be
    ! declared, allocated, initialized, and deallocated whether HOC
    ! is part of a larger model or not.
    call cleanup_prognostic_variables( )
    call cleanup_diagnostic_variables( )

    ! De-allocate the array for the passive scalar tolerances
    deallocate( sclrtol )

    return
  end subroutine cleanup_clubb_core


end module clubb_core
