!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module pdf_closure_module

  ! Options for the two component normal (double Gaussian) PDF type to use for
  ! the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate,
  ! two-component PDF.
  use model_flags, only: &
      iiPDF_ADG1,       & ! ADG1 PDF
      iiPDF_ADG2,       & ! ADG2 PDF
      iiPDF_3D_Luhar,   & ! 3D Luhar PDF
      iiPDF_new,        & ! new PDF
      iiPDF_TSDADG,     & ! new TSDADG PDF
      iiPDF_LY93,       & ! Lewellen and Yoh (1993)
      iiPDF_new_hybrid    ! new hybrid PDF

  implicit none

  public :: pdf_closure, &
            calc_wp4_pdf, &
            calc_wp2xp_pdf, &
            calc_wpxp2_pdf, &
            calc_wpxpyp_pdf, &
            calc_w_up_in_cloud, &
            pdf_closure_driver

  private ! Set Default Scope

  contains
!------------------------------------------------------------------------

  !#######################################################################
  !#######################################################################
  ! If you change the argument list of pdf_closure you also have to
  ! change the calls to this function in the host models CAM, WRF, SAM
  ! and GFDL.
  !#######################################################################
  !#######################################################################
  subroutine pdf_closure( nz, ngrdcol, sclr_dim, sclr_tol, gr,        &
                          hydromet_dim, p_in_Pa, exner, thv_ds,       &
                          wm, wp2, wp3,                               &
                          Skw, Skthl_in, Skrt_in, Sku_in, Skv_in,     &
                          rtm, rtp2, wprtp,                           &
                          thlm, thlp2, wpthlp,                        &
                          um, up2, upwp,                              &
                          vm, vp2, vpwp,                              &
                          rtpthlp,                                    &
                          sclrm, wpsclrp, sclrp2,                     &
                          sclrprtp, sclrpthlp, Sksclr_in,             &
                          gamma_Skw_fnc,                              &
#ifdef GFDL
                          RH_crit, do_liquid_only_in_clubb,           & ! h1g, 2010-06-15
#endif
                          wphydrometp, wp2hmp,                        &
                          rtphmp, thlphmp,                            &
                          clubb_params, mixt_frac_max_mag,            &
                          saturation_formula,                         &
                          stats,                                      &
                          iiPDF_type,                                 &
                          l_mix_rat_hm,                               &
                          sigma_sqd_w,                                &
                          pdf_params, pdf_implicit_coefs_terms,       &
                          err_info,                                   &
                          wpup2, wpvp2,                               &
                          wp2up2, wp2vp2, wp4,                        &
                          wprtp2, wp2rtp,                             &
                          wpthlp2, wp2thlp, wprtpthlp,                &
                          cloud_frac, ice_supersat_frac,              &
                          rcm, wpthvp, wp2thvp, wp2up, rtpthvp,       &
                          thlpthvp, wprcp, wp2rcp, rtprcp,            &
                          thlprcp, rcp2,                              &
                          uprcp, vprcp,                               &
                          w_up_in_cloud, w_down_in_cloud,             &
                          cloudy_updraft_frac, cloudy_downdraft_frac, &
                          wpsclrprtp, wpsclrp2, sclrpthvp,            &
                          wpsclrpthlp, sclrprcp, wp2sclrp,            &
                          rc_coef                                     )


    ! Description:
    ! Subroutine that computes pdf parameters analytically.
    !
    ! Based of the original formulation, but with some tweaks
    ! to remove some of the less realistic assumptions and
    ! improve transport terms.

    !   Corrected version that should remove inconsistency

    ! References:
    !   The shape of CLUBB's PDF is given by the expression in
    !   https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:clubb_pdf

    !   Eqn. 29, 30, 31, 32 & 33  on p. 3547 of
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    !----------------------------------------------------------------------

    use constants_clubb, only: &  ! Constants
        three,          & ! 3
        one,            & ! 1
        one_half,       & ! 1/2
        zero,           & ! 0
        Cp,             & ! Dry air specific heat at constant p [J/kg/K]
        Lv,             & ! Latent heat of vaporization         [J/kg]
        ep1,            & ! (1.0-ep)/ep; ep1 = 0.61             [-]
        ep2,            & ! 1.0/ep;      ep2 = 1.61             [-]
        rt_tol,         & ! Tolerance for r_t                   [kg/kg]
        thl_tol,        & ! Tolerance for th_l                  [K]
        fstderr,        &
        zero_threshold, &
        eps, &
        w_tol

    use grid_class, only: &
        grid

    use parameter_indices, only: &
        nparams,                       & ! Variable(s)
        ibeta,                         &
        iSkw_denom_coef,               &
        ipdf_component_stdev_factor_w, &
        islope_coef_spread_DG_means_w

    use pdf_parameter_module, only: &
        pdf_parameter,        & ! Variable Type
        implicit_coefs_terms

    use new_pdf_main, only: &
        new_pdf_driver    ! Procedure(s)

    use new_hybrid_pdf_main, only: &
        new_hybrid_pdf_driver    ! Procedure(s)

    use adg1_adg2_3d_luhar_pdf, only: &
        ADG1_pdf_driver,     & ! Procedure(s)
        ADG2_pdf_driver,     &
        Luhar_3D_pdf_driver

    use new_tsdadg_pdf, only: &
        tsdadg_pdf_driver    ! Procedure(s)

    use LY93_pdf, only: &
        LY93_driver    ! Procedure(s)

    use pdf_utilities, only: &
        calc_comp_corrs_binormal, & ! Procedure(s)
        calc_corr_chi_x,          &
        calc_corr_eta_x

    use model_flags, only: &
        l_explicit_turbulent_adv_xpyp ! Variable(s)

    use numerical_check, only: &
        pdf_closure_check ! Procedure(s)

    use saturation, only: &
        sat_mixrat_liq_api, & ! Procedure(s)
        sat_mixrat_ice

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api,  & ! Procedure
        clubb_fatal_error              ! Constant

    use stats_netcdf, only: &
      stats_type, &
      stats_update, &
      var_on_stats_list

    use err_info_type_module, only: &
      err_info_type     ! Type

    implicit none

    !----------------------------- Input Variables -----------------------------
    integer, intent(in) :: &
      nz,           & ! Number of vertical levels
      ngrdcol,      & ! Number of grid columns
      hydromet_dim, & ! Number of hydrometeor species
      sclr_dim        ! Number of passive scalars
      
    real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: &
      sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

    type( grid ), intent(in) :: &
      gr

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      p_in_Pa,     & ! Pressure                                   [Pa]
      exner,       & ! Exner function                             [-]
      thv_ds,      & ! Dry, base-state theta_v (ref. th_l here)   [K]
      wm,          & ! mean w-wind component (vertical velocity)  [m/s]
      wp2,         & ! w'^2                                       [m^2/s^2]
      wp3,         & ! w'^3                                       [m^3/s^3]
      Skw,         & ! Skewness of w                              [-]
      Skthl_in,    & ! Skewness of thl                            [-]
      Skrt_in,     & ! Skewness of rt                             [-]
      Sku_in,      & ! Skewness of u                              [-]
      Skv_in,      & ! Skewness of v                              [-]
      rtm,         & ! Mean total water mixing ratio              [kg/kg]
      rtp2,        & ! r_t'^2                                     [(kg/kg)^2]
      wprtp,       & ! w'r_t'                                     [(kg/kg)(m/s)]
      thlm,        & ! Mean liquid water potential temperature    [K]
      thlp2,       & ! th_l'^2                                    [K^2]
      wpthlp,      & ! w'th_l'                                    [K(m/s)]
      rtpthlp        ! r_t'th_l'                                  [K(kg/kg)]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      um,          & ! Grid-mean eastward wind     [m/s]
      up2,         & ! u'^2                        [(m/s)^2]
      upwp,        & ! u'w'                        [(m/s)^2]
      vm,          & ! Grid-mean northward wind    [m/s]
      vp2,         & ! v'^2                        [(m/s)^2]
      vpwp           ! v'w'                        [(m/s)^2]

    real( kind = core_rknd ), dimension(ngrdcol,nz, sclr_dim), intent(in) :: &
      sclrm,       & ! Mean passive scalar        [units vary]
      wpsclrp,     & ! w' sclr'                   [units vary]
      sclrp2,      & ! sclr'^2                    [units vary]
      sclrprtp,    & ! sclr' r_t'                 [units vary]
      sclrpthlp,   & ! sclr' th_l'                [units vary]
      Sksclr_in      ! Skewness of sclr           [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      gamma_Skw_fnc    ! Gamma as a function of skewness            [-]

#ifdef  GFDL
    ! critial relative humidity for nucleation
    real( kind = core_rknd ), dimension(ngrdcol, nz, min(1,sclr_dim), 2 ), & ! h1g, 2010-06-15
       intent(in) :: & ! h1g, 2010-06-15
       RH_crit     ! critical relative humidity for droplet and ice nucleation
! ---> h1g, 2012-06-14
    logical, intent(in) :: do_liquid_only_in_clubb
! <--- h1g, 2012-06-14
#endif

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(in) :: &
      wphydrometp, & ! Covariance of w and a hydrometeor    [(m/s) <hm units>]
      wp2hmp,      & ! Third-order moment:  < w'^2 hm' >    [(m/s)^2 <hm units>]
      rtphmp,      & ! Covariance of rt and a hydrometeor   [(kg/kg) <hm units>]
      thlphmp        ! Covariance of thl and a hydrometeor  [K <hm units>]

    real( kind = core_rknd ), intent(in) :: &
      mixt_frac_max_mag     ! Maximum allowable mag. of mixt_frac   [-]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    type(stats_type), intent(inout) :: &
      stats

    !----------------------------- InOut Variables -----------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: &
      ! If iiPDF_type == iiPDF_ADG2, this gets overwritten. Therefore,
      ! intent(inout). Otherwise it should be intent(in)
      sigma_sqd_w   ! Width of individual w plumes               [-]

    type(pdf_parameter), intent(inout) :: &
      pdf_params     ! pdf paramters         [units vary]

    type(implicit_coefs_terms), intent(inout) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    type(err_info_type), intent(inout) :: &
      err_info      ! err_info struct containing err_code and err_header

    !----------------------------- Output Variables -----------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      wpup2,                 & ! w'u'^2                     [m^3/s^3]
      wpvp2,                 & ! w'v'^2                     [m^3/s^3]
      wp2up2,                & ! w'^2u'^2                   [m^2/s^4]
      wp2vp2,                & ! w'^2v'^2                   [m^2/s^4]
      wp4,                   & ! w'^4                       [m^4/s^4]
      wprtp2,                & ! w' r_t'                    [(m kg)/(s kg)]
      wp2rtp,                & ! w'^2 r_t'                  [(m^2 kg)/(s^2 kg)]
      wpthlp2,               & ! w' th_l'^2                 [(m K^2)/s]
      wp2thlp,               & ! w'^2 th_l'                 [(m^2 K)/s^2]
      cloud_frac,            & ! Cloud fraction             [-]
      ice_supersat_frac,     & ! Ice cloud fracion          [-]
      rcm,                   & ! Mean liquid water          [kg/kg]
      wpthvp,                & ! Buoyancy flux              [(K m)/s] 
      wp2thvp,               & ! w'^2 th_v'                 [(m^2 K)/s^2]
      wp2up,                 & ! w'^2 u'                    [m^3/s^3]
      rtpthvp,               & ! r_t' th_v'                 [(kg K)/kg]
      thlpthvp,              & ! th_l' th_v'                [K^2]
      wprcp,                 & ! w' r_c'                    [(m kg)/(s kg)]
      wp2rcp,                & ! w'^2 r_c'                  [(m^2 kg)/(s^2 kg)]
      rtprcp,                & ! r_t' r_c'                  [(kg^2)/(kg^2)]
      thlprcp,               & ! th_l' r_c'                 [(K kg)/kg]
      rcp2,                  & ! r_c'^2                     [(kg^2)/(kg^2)]
      wprtpthlp,             & ! w' r_t' th_l'              [(m kg K)/(s kg)]
      w_up_in_cloud,         & ! cloudy updraft vel         [m/s]
      w_down_in_cloud,       & ! cloudy downdraft vel       [m/s]
      cloudy_updraft_frac,   & ! cloudy updraft fraction    [-]
      cloudy_downdraft_frac    ! cloudy downdraft fraction  [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      uprcp,              & ! u' r_c'               [(m kg)/(s kg)]
      vprcp                 ! v' r_c'               [(m kg)/(s kg)]

    ! Output (passive scalar variables)
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nz,sclr_dim) :: &
      sclrpthvp, &
      sclrprcp, &
      wpsclrp2, &
      wpsclrprtp, &
      wpsclrpthlp, &
      wp2sclrp

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      rc_coef    ! Coefficient on X'r_c' in X'th_v' equation    [K/(kg/kg)]

    !----------------------------- Local Variables -----------------------------

    ! Variables that are stored in derived data type pdf_params.
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      u_1,           & ! Mean of eastward wind (1st PDF component)         [m/s]
      u_2,           & ! Mean of eastward wind (2nd PDF component)         [m/s]
      varnce_u_1,    & ! Variance of u (1st PDF component)             [m^2/s^2]
      varnce_u_2,    & ! Variance of u (2nd PDF component)             [m^2/s^2]
      v_1,           & ! Mean of northward wind (1st PDF component)        [m/s]
      v_2,           & ! Mean of northward wind (2nd PDF component)        [m/s]
      varnce_v_1,    & ! Variance of v (1st PDF component)             [m^2/s^2]
      varnce_v_2,    & ! Variance of v (2nd PDF component)             [m^2/s^2]
      alpha_u,       & ! Factor relating to normalized variance for u        [-]
      alpha_v          ! Factor relating to normalized variance for v        [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      corr_u_w_1,      & ! Correlation of u and w   (1st PDF component)      [-]
      corr_u_w_2,      & ! Correlation of u and w   (2nd PDF component)      [-]
      corr_v_w_1,      & ! Correlation of v and w   (1st PDF component)      [-]
      corr_v_w_2         ! Correlation of v and w   (2nd PDF component)      [-]

    ! Note:  alpha coefficients = 0.5 * ( 1 - correlations^2 ).
    !        These are used to calculate the scalar widths
    !        varnce_thl_1, varnce_thl_2, varnce_rt_1, and varnce_rt_2 as in
    !        Eq. (34) of Larson and Golaz (2005)

    ! Passive scalar local variables

    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim) :: &
      sclr1, sclr2,  &
      varnce_sclr1, varnce_sclr2, &
      alpha_sclr,  &
      corr_sclr_thl_1, corr_sclr_thl_2, &
      corr_sclr_rt_1, corr_sclr_rt_2, &
      corr_w_sclr_1, corr_w_sclr_2

    logical :: &
      l_scalar_calc, &          ! True if sclr_dim > 0
      l_calc_ice_supersat_frac  ! True if we should calculate ice_supersat_frac

    ! Quantities needed to predict higher order moments
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      tl1, tl2

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      sqrt_wp2, & ! Square root of wp2          [m/s]
      Skthl,    & ! Skewness of thl             [-]
      Skrt,     & ! Skewness of rt              [-]
      Sku,      & ! Skewness of u               [-]
      Skv         ! Skewness of v               [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz,sclr_dim) :: &
      Sksclr      ! Skewness of rt              [-]

    ! Thermodynamic quantity

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      wprcp_contrib_comp_1,   & ! <w'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
      wprcp_contrib_comp_2,   & ! <w'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]
      wp2rcp_contrib_comp_1,  & ! <w'^2rc'> contrib. (1st comp) [m^2/s^2(kg/kg)]
      wp2rcp_contrib_comp_2,  & ! <w'^2rc'> contrib. (2nd comp) [m^2/s^2(kg/kg)]
      rtprcp_contrib_comp_1,  & ! <rt'rc'> contrib. (1st PDF comp.)  [kg^2/kg^2]
      rtprcp_contrib_comp_2,  & ! <rt'rc'> contrib. (2nd PDF comp.)  [kg^2/kg^2]
      thlprcp_contrib_comp_1, & ! <thl'rc'> contrib. (1st PDF comp.)  [K(kg/kg)]
      thlprcp_contrib_comp_2, & ! <thl'rc'> contrib. (2nd PDF comp.)  [K(kg/kg)]
      uprcp_contrib_comp_1,   & ! <u'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
      uprcp_contrib_comp_2,   & ! <u'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]
      vprcp_contrib_comp_1,   & ! <v'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
      vprcp_contrib_comp_2      ! <v'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]

    ! variables for computing ice cloud fraction
    real( kind = core_rknd), dimension(ngrdcol,nz) :: &
      rc_1_ice, rc_2_ice
    
    ! To test pdf parameters
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      wm_clubb_pdf,    &
      rtm_clubb_pdf,   &
      thlm_clubb_pdf,  &
      wp2_clubb_pdf,   &
      rtp2_clubb_pdf,  &
      thlp2_clubb_pdf, &
      wp3_clubb_pdf,   &
      rtp3_clubb_pdf,  &
      thlp3_clubb_pdf, &
      Skw_clubb_pdf,   &
      Skrt_clubb_pdf,  &
      Skthl_clubb_pdf, &
      rsatl_1, &
      rsatl_2

    real( kind = core_rknd ) :: &
      Skw_denom_coef     ! CLUBB tunable parameter beta
      
    logical, parameter :: &
      l_liq_ice_loading_test = .false. ! Temp. flag liq./ice water loading test

    integer :: k, i, sclr, hm_idx   ! Indices

#ifdef GFDL
    real ( kind = core_rknd ), parameter :: t1_combined = 273.16, &
                                            t2_combined = 268.16, &
                                            t3_combined = 238.16 
#endif

    !----------------------------- Begin Code -----------------------------

    !$acc enter data create( u_1, u_2, varnce_u_1, varnce_u_2, v_1, v_2, &
    !$acc                 varnce_v_1, varnce_v_2, alpha_u, alpha_v, &
    !$acc                 corr_u_w_1, corr_u_w_2, corr_v_w_1, corr_v_w_2, &
    !$acc                 tl1, tl2, sqrt_wp2, Skthl, &
    !$acc                 Skrt, Sku, Skv, wprcp_contrib_comp_1, wprcp_contrib_comp_2, &
    !$acc                 wp2rcp_contrib_comp_1, wp2rcp_contrib_comp_2, &
    !$acc                 rtprcp_contrib_comp_1, rtprcp_contrib_comp_2, &
    !$acc                 thlprcp_contrib_comp_1, thlprcp_contrib_comp_2, &
    !$acc                 uprcp_contrib_comp_1, uprcp_contrib_comp_2, &
    !$acc                 vprcp_contrib_comp_1, vprcp_contrib_comp_2, &
    !$acc                 rc_1_ice, rc_2_ice, rsatl_1, rsatl_2 )

    !$acc enter data if( sclr_dim > 0 ) &
    !$acc            create( sclr1, sclr2, varnce_sclr1, varnce_sclr2, &
    !$acc                    alpha_sclr, corr_sclr_thl_1, corr_sclr_thl_2, &
    !$acc                    corr_sclr_rt_1, corr_sclr_rt_2, corr_w_sclr_1, &
    !$acc                    corr_w_sclr_2, Sksclr )

    ! Check whether the passive scalars are present.
    if ( sclr_dim > 0 ) then
      l_scalar_calc = .true.
    else
      l_scalar_calc = .false.
    end if

    ! Initialize to default values to prevent a runtime error
    if ( ( iiPDF_type /= iiPDF_ADG1 ) .and. ( iiPDF_type /= iiPDF_ADG2 ) ) then

      do k = 1, nz
        do i = 1, ngrdcol
          pdf_params%alpha_thl(i,k) = one_half
          pdf_params%alpha_rt(i,k) = one_half
        end do
      end do
      
      ! This allows for skewness to be clipped locally without passing the updated
      ! value back out.
      do k = 1, nz
        do i = 1, ngrdcol
          Skrt(i,k) = Skrt_in(i,k)
          Skthl(i,k) = Skthl_in(i,k)
          Sku(i,k) = Sku_in(i,k)
          Skv(i,k) = Skv_in(i,k)
        end do
      end do

      do sclr = 1, sclr_dim
        do k = 1, nz
          do i = 1, ngrdcol
            
            Sksclr(i,k,sclr) = Sksclr_in(i,k,sclr)
            
            if ( l_scalar_calc ) then
                alpha_sclr(i,k,sclr) = one_half
            end if
            
          end do
        end do
      end do

    end if

    ! To avoid recomputing
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        sqrt_wp2(i,k) = sqrt( wp2(i,k) )
      end do
    end do
    !$acc end parallel loop

    ! Select the PDF closure method for the two-component PDF used by CLUBB for
    ! w, rt, theta-l, and passive scalar variables.
    ! Calculate the mixture fraction for the multivariate PDF, as well as both
    ! PDF component means and both PDF component variances for each of w, rt,
    ! theta-l, and passive scalar variables.
    if ( iiPDF_type == iiPDF_ADG1 ) then ! use ADG1

      call ADG1_pdf_driver( nz, ngrdcol, sclr_dim, sclr_tol,                        & ! In
                            wm, rtm, thlm, um, vm,                                  & ! In
                            wp2, rtp2, thlp2, up2, vp2,                             & ! In
                            Skw, wprtp, wpthlp, upwp, vpwp, sqrt_wp2,               & ! In
                            sigma_sqd_w, clubb_params(:,ibeta), mixt_frac_max_mag,  & ! In
                            sclrm, sclrp2, wpsclrp, l_scalar_calc,                  & ! In
                            pdf_params%w_1, pdf_params%w_2,                         & ! Out
                            pdf_params%rt_1, pdf_params%rt_2,                       & ! Out
                            pdf_params%thl_1, pdf_params%thl_2,                     & ! Out
                            u_1, u_2, v_1, v_2,                                     & ! Out
                            pdf_params%varnce_w_1, pdf_params%varnce_w_2,           & ! Out
                            pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,         & ! Out
                            pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,       & ! Out
                            varnce_u_1, varnce_u_2,                                 & ! Out
                            varnce_v_1, varnce_v_2,                                 & ! Out
                            pdf_params%mixt_frac,                                   & ! Out
                            pdf_params%alpha_rt, pdf_params%alpha_thl,              & ! Out
                            alpha_u, alpha_v,                                       & ! Out
                            sclr1, sclr2, varnce_sclr1,                             & ! Out
                            varnce_sclr2, alpha_sclr )                                ! Out

    elseif ( iiPDF_type == iiPDF_ADG2 ) then ! use ADG2

      call ADG2_pdf_driver( nz, ngrdcol, sclr_dim, sclr_tol,                      & ! In
                            wm, rtm, thlm, wp2, rtp2, thlp2,                      & ! In
                            Skw, wprtp, wpthlp, sqrt_wp2, clubb_params(:,ibeta),  & ! In
                            sclrm, sclrp2, wpsclrp, l_scalar_calc,                & ! In
                            pdf_params%w_1, pdf_params%w_2,                       & ! Out
                            pdf_params%rt_1, pdf_params%rt_2,                     & ! Out
                            pdf_params%thl_1, pdf_params%thl_2,                   & ! Out
                            pdf_params%varnce_w_1, pdf_params%varnce_w_2,         & ! Out
                            pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,       & ! Out
                            pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,     & ! Out
                            pdf_params%mixt_frac,                                 & ! Out
                            pdf_params%alpha_rt, pdf_params%alpha_thl,            & ! Out
                            sigma_sqd_w, sclr1, sclr2,                            & ! Out
                            varnce_sclr1, varnce_sclr2, alpha_sclr )                ! Out

    elseif ( iiPDF_type == iiPDF_3D_Luhar ) then ! use 3D Luhar
      do i = 1, ngrdcol
        call Luhar_3D_pdf_driver( nz, &
                           wm(i,:), rtm(i,:), thlm(i,:), wp2(i,:), rtp2(i,:), thlp2(i,:), & ! In
                           Skw(i,:), Skrt(i,:), Skthl(i,:), wprtp(i,:), wpthlp(i,:),      & ! In
                           pdf_params%w_1(i,:), pdf_params%w_2(i,:),                      & ! Out
                           pdf_params%rt_1(i,:), pdf_params%rt_2(i,:),                    & ! Out
                           pdf_params%thl_1(i,:), pdf_params%thl_2(i,:),                  & ! Out
                           pdf_params%varnce_w_1(i,:), pdf_params%varnce_w_2(i,:),        & ! Out
                           pdf_params%varnce_rt_1(i,:), pdf_params%varnce_rt_2(i,:),      & ! Out
                           pdf_params%varnce_thl_1(i,:), pdf_params%varnce_thl_2(i,:),    & ! Out
                           pdf_params%mixt_frac(i,:) )                                      ! Out
      end do
    elseif ( iiPDF_type == iiPDF_new ) then ! use new PDF
      call new_pdf_driver( nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2,      & ! In
                           Skw,                                               & ! In
                           wprtp, wpthlp, rtpthlp,                            & ! In
                           clubb_params,                                      & ! In
                           Skrt, Skthl,                                       & ! In/Out
                           pdf_params%w_1, pdf_params%w_2,                    & ! Out
                           pdf_params%rt_1, pdf_params%rt_2,                  & ! Out
                           pdf_params%thl_1, pdf_params%thl_2,                & ! Out
                           pdf_params%varnce_w_1, pdf_params%varnce_w_2,      & ! Out
                           pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,    & ! Out
                           pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,  & ! Out
                           pdf_params%mixt_frac,                              & ! Out
                           pdf_implicit_coefs_terms )                           ! Out
    elseif ( iiPDF_type == iiPDF_TSDADG ) then
      do i = 1, ngrdcol
        call tsdadg_pdf_driver( nz, &
                          wm(i,:), rtm(i,:), thlm(i,:), wp2(i,:), rtp2(i,:), thlp2(i,:),  & ! In
                          Skw(i,:), Skrt(i,:), Skthl(i,:), wprtp(i,:), wpthlp(i,:),       & ! In
                          pdf_params%w_1(i,:), pdf_params%w_2(i,:),                       & ! Out
                          pdf_params%rt_1(i,:), pdf_params%rt_2(i,:),                     & ! Out
                          pdf_params%thl_1(i,:), pdf_params%thl_2(i,:),                   & ! Out
                          pdf_params%varnce_w_1(i,:), pdf_params%varnce_w_2(i,:),         & ! Out
                          pdf_params%varnce_rt_1(i,:), pdf_params%varnce_rt_2(i,:),       & ! Out
                          pdf_params%varnce_thl_1(i,:), pdf_params%varnce_thl_2(i,:),     & ! Out
                          pdf_params%mixt_frac(i,:) )                                       ! Out
      end do
    elseif ( iiPDF_type == iiPDF_LY93 ) then ! use LY93
      do i = 1, ngrdcol
        call LY93_driver( nz, wm(i,:), rtm(i,:), thlm(i,:), wp2(i,:), rtp2(i,:),       & ! In
                          thlp2(i,:), Skw(i,:), Skrt(i,:), Skthl(i,:),                 & ! In
                          pdf_params%w_1(i,:), pdf_params%w_2(i,:),                    & ! Out
                          pdf_params%rt_1(i,:), pdf_params%rt_2(i,:),                  & ! Out
                          pdf_params%thl_1(i,:), pdf_params%thl_2(i,:),                & ! Out
                          pdf_params%varnce_w_1(i,:), pdf_params%varnce_w_2(i,:),      & ! Out
                          pdf_params%varnce_rt_1(i,:), pdf_params%varnce_rt_2(i,:),    & ! Out
                          pdf_params%varnce_thl_1(i,:), pdf_params%varnce_thl_2(i,:),  & ! Out
                          pdf_params%mixt_frac(i,:) )                                    ! Out
      end do
    elseif ( iiPDF_type == iiPDF_new_hybrid ) then ! use new hybrid PDF
      call new_hybrid_pdf_driver( nz, ngrdcol, sclr_dim,                          & ! In
                                  wm, rtm, thlm, um, vm,                          & ! In
                                  wp2, rtp2, thlp2, up2, vp2,                     & ! In
                                  Skw, wprtp, wpthlp, upwp, vpwp,                 & ! In
                                  sclrm, sclrp2, wpsclrp,                         & ! In
                                  gamma_Skw_fnc,                                  & ! In
                                  clubb_params(:,islope_coef_spread_DG_means_w),  & ! In
                                  clubb_params(:,ipdf_component_stdev_factor_w),  & ! In
                                  Skrt, Skthl, Sku, Skv, Sksclr,                  & ! I/O
                                  pdf_params%w_1, pdf_params%w_2,                 & ! Out
                                  pdf_params%rt_1, pdf_params%rt_2,               & ! Out
                                  pdf_params%thl_1, pdf_params%thl_2,             & ! Out
                                  u_1, u_2, v_1, v_2,                             & ! Out
                                  pdf_params%varnce_w_1,                          & ! Out
                                  pdf_params%varnce_w_2,                          & ! Out
                                  pdf_params%varnce_rt_1,                         & ! Out
                                  pdf_params%varnce_rt_2,                         & ! Out
                                  pdf_params%varnce_thl_1,                        & ! Out
                                  pdf_params%varnce_thl_2,                        & ! Out
                                  varnce_u_1, varnce_u_2,                         & ! Out
                                  varnce_v_1, varnce_v_2,                         & ! Out
                                  sclr1, sclr2,                                   & ! Out
                                  varnce_sclr1, varnce_sclr2,                     & ! Out
                                  pdf_params%mixt_frac, sigma_sqd_w,              & ! Out
                                  pdf_implicit_coefs_terms )                        ! Out

    end if ! iiPDF_type

    ! Calculate the PDF component correlations of rt and thl.
    call calc_comp_corrs_binormal( nz, ngrdcol,                                           & ! In
                                   rtpthlp, rtm, thlm,                                    & ! In
                                   pdf_params%rt_1, pdf_params%rt_2,                      & ! In
                                   pdf_params%thl_1, pdf_params%thl_2,                    & ! In
                                   pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,        & ! In
                                   pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,      & ! In
                                   pdf_params%mixt_frac,                                  & ! In
                                   pdf_params%corr_rt_thl_1, pdf_params%corr_rt_thl_2 )     ! Out

    if ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
         .or. iiPDF_type == iiPDF_new_hybrid ) then

      ! These PDF types define corr_w_rt_1, corr_w_rt_2, corr_w_thl_1, and
      ! corr_w_thl_2 to all have a value of 0, so skip the calculation.
      ! The values of corr_u_w_1, corr_u_w_2, corr_v_w_1, and corr_v_w_2 are
      ! all defined to be 0, as well.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          pdf_params%corr_w_rt_1(i,k)  = zero
          pdf_params%corr_w_rt_2(i,k)  = zero
          pdf_params%corr_w_thl_1(i,k) = zero
          pdf_params%corr_w_thl_2(i,k) = zero
          corr_u_w_1(i,k)   = zero
          corr_u_w_2(i,k)   = zero
          corr_v_w_1(i,k)   = zero
          corr_v_w_2(i,k)   = zero
        end do
      end do
      !$acc end parallel loop

    else

      ! Calculate the PDF component correlations of w and rt.
      call calc_comp_corrs_binormal( nz, ngrdcol,                                      & ! In
                                     wprtp, wm, rtm,                                   & ! In
                                     pdf_params%w_1, pdf_params%w_2,                   & ! In
                                     pdf_params%rt_1, pdf_params%rt_2,                 & ! In
                                     pdf_params%varnce_w_1, pdf_params%varnce_w_2,     & ! In
                                     pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,   & ! In
                                     pdf_params%mixt_frac,                             & ! In
                                     pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2 )    ! Out

      ! Calculate the PDF component correlations of w and thl.
      call calc_comp_corrs_binormal( nz, ngrdcol,                                          & ! In
                                     wpthlp, wm, thlm,                                     & ! In
                                     pdf_params%w_1, pdf_params%w_2,                       & ! In
                                     pdf_params%thl_1, pdf_params%thl_2,                   & ! In
                                     pdf_params%varnce_w_1, pdf_params%varnce_w_2,         & ! In
                                     pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,     & ! In
                                     pdf_params%mixt_frac,                                 & ! In
                                     pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2 )      ! Out
    end if

    if ( l_scalar_calc ) then

      ! Calculate the PDF component correlations of a passive scalar and thl.
      do sclr = 1, sclr_dim
        call calc_comp_corrs_binormal( nz, ngrdcol,                                         & ! In
                                       sclrpthlp(:,:,sclr), sclrm(:,:,sclr), thlm,          & ! In
                                       sclr1(:,:,sclr), sclr2(:,:,sclr),                    & ! In
                                       pdf_params%thl_1, pdf_params%thl_2,                  & ! In
                                       varnce_sclr1(:,:,sclr), varnce_sclr2(:,:,sclr),      & ! In
                                       pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,    & ! In
                                       pdf_params%mixt_frac,                                & ! In
                                       corr_sclr_thl_1(:,:,sclr), corr_sclr_thl_2(:,:,sclr) ) ! Out
      end do

      ! Calculate the PDF component correlations of a passive scalar and rt.
      do sclr = 1, sclr_dim
        call calc_comp_corrs_binormal( nz, ngrdcol,                                       & ! In
                                       sclrprtp(:,:,sclr), sclrm(:,:,sclr), rtm,          & ! In
                                       sclr1(:,:,sclr), sclr2(:,:,sclr),                  & ! In
                                       pdf_params%rt_1, pdf_params%rt_2,                  & ! In
                                       varnce_sclr1(:,:,sclr), varnce_sclr2(:,:,sclr),    & ! In
                                       pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,    & ! In
                                       pdf_params%mixt_frac,                              & ! In
                                       corr_sclr_rt_1(:,:,sclr), corr_sclr_rt_2(:,:,sclr) ) ! Out
      end do

      if ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
           .or. iiPDF_type == iiPDF_new_hybrid ) then

        ! These PDF types define all PDF component correlations involving w
        ! to have a value of 0, so skip the calculation.
        !$acc parallel loop gang vector collapse(2) default(present)
        do sclr = 1, sclr_dim
          do k = 1, nz
            do i = 1, ngrdcol
              corr_w_sclr_1(i,k,sclr) = zero
              corr_w_sclr_2(i,k,sclr) = zero
            end do
          end do
        end do
        !$acc end parallel loop

      else

        ! Calculate the PDF component correlations of w and a passive scalar.
        do sclr = 1, sclr_dim
          call calc_comp_corrs_binormal( nz, ngrdcol,                                      & ! In
                                         wpsclrp(:,:,sclr), wm, sclrm(:,:,sclr),           & ! In
                                         pdf_params%w_1, pdf_params%w_2,                   & ! In
                                         sclr1(:,:,sclr), sclr2(:,:,sclr),                 & ! In
                                         pdf_params%varnce_w_1, pdf_params%varnce_w_2,     & ! In
                                         varnce_sclr1(:,:,sclr), varnce_sclr2(:,:,sclr),   & ! In
                                         pdf_params%mixt_frac,                             & ! In
                                         corr_w_sclr_1(:,:,sclr), corr_w_sclr_2(:,:,sclr) )  ! Out

        end do
        
      end if

    end if


    ! Compute higher order moments (these are interactive)
    call calc_wp2xp_pdf( nz, ngrdcol,                                       &
                         wm, rtm, pdf_params%w_1, pdf_params%w_2,           &
                         pdf_params%rt_1, pdf_params%rt_2,                  &
                         pdf_params%varnce_w_1, pdf_params%varnce_w_2,      &
                         pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,    &
                         pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2,    &
                         pdf_params%mixt_frac,                              &
                         wp2rtp )

    call calc_wp2xp_pdf( nz, ngrdcol,                                          &
                         wm, thlm, pdf_params%w_1, pdf_params%w_2,             &
                         pdf_params%thl_1, pdf_params%thl_2,                   &
                         pdf_params%varnce_w_1, pdf_params%varnce_w_2,         &
                         pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,     &
                         pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2,     &
                         pdf_params%mixt_frac,                                 &
                         wp2thlp )

    call calc_wp2xp_pdf( nz, ngrdcol,                                          &
                         wm, um, pdf_params%w_1, pdf_params%w_2,               &
                         u_1, u_2,                                             &
                         pdf_params%varnce_w_1, pdf_params%varnce_w_2,         &
                         varnce_u_1, varnce_u_2,                               &
                         corr_u_w_1, corr_u_w_2,                               &
                         pdf_params%mixt_frac,                                 &
                         wp2up )

    ! Compute higher order moments (these may be interactive)
    call calc_wpxp2_pdf( nz, ngrdcol, &
                         wm, um, pdf_params%w_1, pdf_params%w_2, &
                         u_1, u_2, &
                         pdf_params%varnce_w_1, pdf_params%varnce_w_2, &
                         varnce_u_1, varnce_u_2, &
                         corr_u_w_1, corr_u_w_2, &
                         pdf_params%mixt_frac, &
                         wpup2 )

    call calc_wpxp2_pdf( nz, ngrdcol, &
                         wm, vm, pdf_params%w_1, pdf_params%w_2, &
                         v_1, v_2, &
                         pdf_params%varnce_w_1, pdf_params%varnce_w_2, &
                         varnce_v_1, varnce_v_2, &
                         corr_v_w_1, corr_v_w_2, &
                         pdf_params%mixt_frac, &
                         wpvp2 )

    call calc_wp2xp2_pdf( nz, ngrdcol, &
                          wm, um, pdf_params%w_1, pdf_params%w_2, &
                          u_1, u_2, &
                          pdf_params%varnce_w_1, pdf_params%varnce_w_2, &
                          varnce_u_1, varnce_u_2, &
                          corr_u_w_1, corr_u_w_2, &
                          pdf_params%mixt_frac, &
                          wp2up2 )

    call calc_wp2xp2_pdf( nz, ngrdcol, &
                          wm, vm, pdf_params%w_1, pdf_params%w_2, &
                          v_1, v_2, &
                          pdf_params%varnce_w_1, pdf_params%varnce_w_2, &
                          varnce_v_1, varnce_v_2, &
                          corr_v_w_1, corr_v_w_2, &
                          pdf_params%mixt_frac, &
                          wp2vp2 )

    call calc_wp4_pdf( nz, ngrdcol, &
                       wm, pdf_params%w_1, pdf_params%w_2, &
                       pdf_params%varnce_w_1, pdf_params%varnce_w_2, &
                       pdf_params%mixt_frac, &
                       wp4 )

    if ( l_explicit_turbulent_adv_xpyp .or. &
         ( var_on_stats_list( stats, "wprtp2" )) ) then
      call calc_wpxp2_pdf( nz, ngrdcol, &
                           wm, rtm, pdf_params%w_1, pdf_params%w_2,        &
                           pdf_params%rt_1, pdf_params%rt_2,               &
                           pdf_params%varnce_w_1, pdf_params%varnce_w_2,   &
                           pdf_params%varnce_rt_1, pdf_params%varnce_rt_2, &
                           pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2, &
                           pdf_params%mixt_frac, &
                           wprtp2 )
    end if

    if ( l_explicit_turbulent_adv_xpyp .or. &
         ( var_on_stats_list( stats, "wpthlp2" )) ) then
      call calc_wpxp2_pdf( nz, ngrdcol, &
                           wm, thlm, pdf_params%w_1, pdf_params%w_2,          &
                           pdf_params%thl_1, pdf_params%thl_2,                &
                           pdf_params%varnce_w_1, pdf_params%varnce_w_2,      &
                           pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,  &
                           pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2,  &
                           pdf_params%mixt_frac, &
                           wpthlp2 )
    end if

    if ( l_explicit_turbulent_adv_xpyp .or. &
         ( var_on_stats_list( stats, "wprtpthlp" )) ) then
      
      call calc_wpxpyp_pdf( nz, ngrdcol, &
                            wm, rtm, thlm, pdf_params%w_1, pdf_params%w_2,      &
                            pdf_params%rt_1, pdf_params%rt_2,                   &
                            pdf_params%thl_1, pdf_params%thl_2,                 &
                            pdf_params%varnce_w_1, pdf_params%varnce_w_2,       &
                            pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,     &
                            pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,   &
                            pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2,     &
                            pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2,   &
                            pdf_params%corr_rt_thl_1, pdf_params%corr_rt_thl_2, &
                            pdf_params%mixt_frac, &
                            wprtpthlp )
    end if


    ! Scalar Addition to higher order moments
    if ( l_scalar_calc ) then

      do sclr = 1, sclr_dim
        call calc_wp2xp_pdf( nz, ngrdcol,                                          &
                             wm, sclrm(:,:,sclr), pdf_params%w_1, pdf_params%w_2,  &
                             sclr1(:,:,sclr), sclr2(:,:,sclr),                     &
                             pdf_params%varnce_w_1, pdf_params%varnce_w_2,         &
                             varnce_sclr1(:,:,sclr), varnce_sclr2(:,:,sclr),       &
                             corr_w_sclr_1(:,:,sclr), corr_w_sclr_2(:,:,sclr),     &
                             pdf_params%mixt_frac,                                 &
                             wp2sclrp(:,:,sclr) )
      end do

      do sclr = 1, sclr_dim
        call calc_wpxp2_pdf( nz, ngrdcol, &
                             wm, sclrm(:,:,sclr), pdf_params%w_1, pdf_params%w_2,  &
                             sclr1(:,:,sclr), sclr2(:,:,sclr),                     &
                             pdf_params%varnce_w_1, pdf_params%varnce_w_2,         &
                             varnce_sclr1(:,:,sclr), varnce_sclr2(:,:,sclr),       &
                             corr_w_sclr_1(:,:,sclr), corr_w_sclr_2(:,:,sclr),     &
                             pdf_params%mixt_frac, &
                             wpsclrp2(:,:,sclr) )
      end do

      do sclr = 1, sclr_dim
        call calc_wpxpyp_pdf( nz, ngrdcol, &
                              wm, sclrm(:,:,sclr), rtm, pdf_params%w_1, pdf_params%w_2,  &
                              sclr1(:,:,sclr), sclr2(:,:,sclr),                          &
                              pdf_params%rt_1, pdf_params%rt_2,                          &
                              pdf_params%varnce_w_1, pdf_params%varnce_w_2,              &
                              varnce_sclr1(:,:,sclr), varnce_sclr2(:,:,sclr),            &
                              pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,            &
                              corr_w_sclr_1(:,:,sclr), corr_w_sclr_2(:,:,sclr),          &
                              pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2,            &
                              corr_sclr_rt_1(:,:,sclr), corr_sclr_rt_2(:,:,sclr),        &
                              pdf_params%mixt_frac, &
                              wpsclrprtp(:,:,sclr) )
      end do

      do sclr = 1, sclr_dim
        call calc_wpxpyp_pdf( nz, ngrdcol, &
                              wm, sclrm(:,:,sclr), thlm, pdf_params%w_1, pdf_params%w_2,   &
                              sclr1(:,:,sclr), sclr2(:,:,sclr),                            &
                              pdf_params%thl_1, pdf_params%thl_2,                          &
                              pdf_params%varnce_w_1, pdf_params%varnce_w_2,                &
                              varnce_sclr1(:,:,sclr), varnce_sclr2(:,:,sclr),              &
                              pdf_params%varnce_thl_1, pdf_params%varnce_thl_2,            &
                              corr_w_sclr_1(:,:,sclr), corr_w_sclr_2(:,:,sclr),            &
                              pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2,            &
                              corr_sclr_thl_1(:,:,sclr), corr_sclr_thl_2(:,:,sclr),        &
                              pdf_params%mixt_frac, &
                              wpsclrpthlp(:,:,sclr) )
      end do

    end if

    ! Compute higher order moments that include theta_v.

    ! First compute some preliminary quantities.
    ! "1" denotes first Gaussian; "2" denotes 2nd Gaussian
    ! liq water temp (Sommeria & Deardorff 1977 (SD), eqn. 3)
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        tl1(i,k)  = pdf_params%thl_1(i,k)*exner(i,k)
        tl2(i,k)  = pdf_params%thl_2(i,k)*exner(i,k)
      end do
    end do
    !$acc end parallel loop

#ifdef GFDL
    if ( sclr_dim > 0  .and.  (.not. do_liquid_only_in_clubb) ) then ! h1g, 2010-06-16 begin mod

      do i = 1, ngrdcol
        where ( tl1(i,:) > t1_combined )
          pdf_params%rsatl_1(i,:) = sat_mixrat_liq_api( nz, p_in_Pa(i,:), tl1(i,:) )
        elsewhere ( tl1(i,:) > t2_combined )
          pdf_params%rsatl_1(i,:) = sat_mixrat_liq_api( nz, p_in_Pa(i,:), tl1(i,:) ) &
                    * (tl1(i,:) - t2_combined)/(t1_combined - t2_combined) &
                    + sat_mixrat_ice( nz, p_in_Pa(i,:), tl1(i,:) ) &
                      * (t1_combined - tl1(i,:))/(t1_combined - t2_combined)
        elsewhere ( tl1(i,:) > t3_combined )
          pdf_params%rsatl_1(i,:) = sat_mixrat_ice( nz, p_in_Pa(i,:), tl1(i,:) ) &
                    + sat_mixrat_ice( nz, p_in_Pa(i,:), tl1(i,:) ) * (RH_crit(i, :, 1, 1) -one ) &
                      * ( t2_combined -tl1(i,:))/(t2_combined - t3_combined)
        elsewhere
          pdf_params%rsatl_1(i,:) = sat_mixrat_ice( nz, p_in_Pa(i,:), tl1(i,:) ) &
                                    * RH_crit(i, :, 1, 1)
        endwhere

        where ( tl2(i,:) > t1_combined )
          pdf_params%rsatl_2(i,:) = sat_mixrat_liq_api( nz, p_in_Pa(i,:), tl2(i,:) )
        elsewhere ( tl2(i,:) > t2_combined )
          pdf_params%rsatl_2(i,:) = sat_mixrat_liq_api( nz, p_in_Pa(i,:), tl2(i,:) ) &
                    * (tl2(i,:) - t2_combined)/(t1_combined - t2_combined) &
                    + sat_mixrat_ice( nz, p_in_Pa(i,:), tl2(i,:) ) &
                      * (t1_combined - tl2(i,:))/(t1_combined - t2_combined)
        elsewhere ( tl2(i,:) > t3_combined )
          pdf_params%rsatl_2(i,:) = sat_mixrat_ice( nz, p_in_Pa(i,:), tl2(i,:) ) &
                    + sat_mixrat_ice( nz, p_in_Pa(i,:), tl2(i,:) )* (RH_crit(i, :, 1, 2) -one) &
                      * ( t2_combined -tl2(i,:))/(t2_combined - t3_combined)
        elsewhere
          pdf_params%rsatl_2(i,:) = sat_mixrat_ice( nz, p_in_Pa(i,:), tl2(i,:) ) &
                                    * RH_crit(i, :, 1, 2)
        endwhere

      end do

    else ! sclr_dim <= 0  or  do_liquid_only_in_clubb = .T.

      pdf_params%rsatl_1 = sat_mixrat_liq_api( nz, ngrdcol, p_in_Pa, tl1, saturation_formula )
      pdf_params%rsatl_2 = sat_mixrat_liq_api( nz, ngrdcol, p_in_Pa, tl2, saturation_formula )

    end if !sclr_dim > 0

    ! Determine whether to compute ice_supersat_frac. We do not compute
    ! ice_supersat_frac for GFDL (unless do_liquid_only_in_clubb is true),
    ! because liquid and ice are both fed into rtm, ruining the calculation.
    if (do_liquid_only_in_clubb) then
      l_calc_ice_supersat_frac = .true.
    else
      l_calc_ice_supersat_frac = .false.
    end if

#else
    rsatl_1 = sat_mixrat_liq_api( nz, ngrdcol, gr, p_in_Pa, tl1, saturation_formula  )
    rsatl_2 = sat_mixrat_liq_api( nz, ngrdcol, gr, p_in_Pa, tl2, &
                                  saturation_formula  ) ! h1g, 2010-06-16 end mod

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        pdf_params%rsatl_1(i,k) = rsatl_1(i,k)
        pdf_params%rsatl_2(i,k) = rsatl_2(i,k)
      end do
    end do
    !$acc end parallel loop

    l_calc_ice_supersat_frac = .true.
#endif

    call transform_pdf_chi_eta_component( nz, ngrdcol, &
                                          tl1, pdf_params%rsatl_1, pdf_params%rt_1, exner,  & ! In
                                          pdf_params%varnce_thl_1, pdf_params%varnce_rt_1,  & ! In
                                          pdf_params%corr_rt_thl_1, pdf_params%chi_1,       & ! In
                                          pdf_params%crt_1, pdf_params%cthl_1,              & ! Out
                                          pdf_params%stdev_chi_1, pdf_params%stdev_eta_1,   & ! Out
                                          pdf_params%covar_chi_eta_1,                       & ! Out
                                          pdf_params%corr_chi_eta_1 )                         ! Out


    ! Calculate cloud fraction component for pdf 1
    call calc_liquid_cloud_frac_component( nz, ngrdcol, &
                                           pdf_params%chi_1, pdf_params%stdev_chi_1,    & ! In
                                           pdf_params%cloud_frac_1, pdf_params%rc_1 )     ! Out

    ! Calc ice_supersat_frac
    if ( l_calc_ice_supersat_frac ) then

      call calc_ice_cloud_frac_component( nz, ngrdcol, &
                                          pdf_params%chi_1, pdf_params%stdev_chi_1, &
                                          pdf_params%rc_1, pdf_params%cloud_frac_1, &
                                          p_in_Pa, tl1, &
                                          pdf_params%rsatl_1, pdf_params%crt_1, &
                                          saturation_formula, &
                                          pdf_params%ice_supersat_frac_1, rc_1_ice )
    end if

    call transform_pdf_chi_eta_component( nz, ngrdcol, &
                                          tl2, pdf_params%rsatl_2, pdf_params%rt_2, exner,  & ! In
                                          pdf_params%varnce_thl_2, pdf_params%varnce_rt_2,  & ! In
                                          pdf_params%corr_rt_thl_2, pdf_params%chi_2,       & ! In
                                          pdf_params%crt_2, pdf_params%cthl_2,              & ! Out
                                          pdf_params%stdev_chi_2, pdf_params%stdev_eta_2,   & ! Out
                                          pdf_params%covar_chi_eta_2,                       & ! Out
                                          pdf_params%corr_chi_eta_2 )                         ! Out


    ! Calculate cloud fraction component for pdf 2
    call calc_liquid_cloud_frac_component( nz, ngrdcol, &
                                           pdf_params%chi_2, pdf_params%stdev_chi_2,    & ! In
                                           pdf_params%cloud_frac_2, pdf_params%rc_2 )     ! Out

    ! Calc ice_supersat_frac
    if ( l_calc_ice_supersat_frac ) then

      call calc_ice_cloud_frac_component( nz, ngrdcol, &
                                          pdf_params%chi_2, pdf_params%stdev_chi_2, &
                                          pdf_params%rc_2, pdf_params%cloud_frac_2, &
                                          p_in_Pa, tl2, &
                                          pdf_params%rsatl_2, pdf_params%crt_2, &
                                          saturation_formula, &
                                          pdf_params%ice_supersat_frac_2, rc_2_ice )

      ! Compute ice cloud fraction, ice_supersat_frac
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          ice_supersat_frac(i,k) = pdf_params%mixt_frac(i,k) &
                                   * pdf_params%ice_supersat_frac_1(i,k) &
                                   + ( one - pdf_params%mixt_frac(i,k) ) &
                                     * pdf_params%ice_supersat_frac_2(i,k)
        end do
      end do
      !$acc end parallel loop

    else 

      ! ice_supersat_frac will be garbage if computed as above
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          ice_supersat_frac(i,k) = 0.0_core_rknd
        end do
      end do
      !$acc end parallel loop

      if (clubb_at_least_debug_level_api( 1 )) then
          write(fstderr,*) "Warning: ice_supersat_frac has garbage values if &
                          & do_liquid_only_in_clubb = .false."
      end if

    end if ! l_calc_ice_supersat_frac


    ! Compute cloud fraction and mean cloud water mixing ratio.
    ! Reference:
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:anl_int_cloud_terms
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        cloud_frac(i,k) = pdf_params%mixt_frac(i,k) * pdf_params%cloud_frac_1(i,k) &
                     + ( one - pdf_params%mixt_frac(i,k) ) * pdf_params%cloud_frac_2(i,k)
        rcm(i,k) = pdf_params%mixt_frac(i,k) * pdf_params%rc_1(i,k) &
                   + ( one - pdf_params%mixt_frac(i,k) ) * pdf_params%rc_2(i,k)
        rcm(i,k) = max( zero_threshold, rcm(i,k) )
      end do
    end do
    !$acc end parallel loop

    if ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
         .or. iiPDF_type == iiPDF_new_hybrid ) then

      ! corr_w_rt and corr_w_thl are zero for these pdf types so
      ! corr_w_chi and corr_w_eta are zero as well
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nz
        do i = 1, ngrdcol
          pdf_params%corr_w_chi_1(i,k) = zero
          pdf_params%corr_w_chi_2(i,k) = zero
          pdf_params%corr_w_eta_1(i,k) = zero
          pdf_params%corr_w_eta_2(i,k) = zero
        end do
      end do
      !$acc end parallel loop

    else
      ! Correlation of w and chi for each component.
      pdf_params%corr_w_chi_1 &
      = calc_corr_chi_x( pdf_params%crt_1, pdf_params%cthl_1, &
                         sqrt(pdf_params%varnce_rt_1), sqrt(pdf_params%varnce_thl_1), &
                         pdf_params%stdev_chi_1, &
                         pdf_params%corr_w_rt_1, pdf_params%corr_w_thl_1 )

      pdf_params%corr_w_chi_2 &
      = calc_corr_chi_x( pdf_params%crt_2, pdf_params%cthl_2, &
                         sqrt(pdf_params%varnce_rt_2), sqrt(pdf_params%varnce_thl_2), &
                         pdf_params%stdev_chi_2, pdf_params%corr_w_rt_2, &
                         pdf_params%corr_w_thl_2 )

      ! Correlation of w and eta for each component.
      pdf_params%corr_w_eta_1 &
      = calc_corr_eta_x( pdf_params%crt_1, pdf_params%cthl_1, &
                         sqrt(pdf_params%varnce_rt_1), sqrt(pdf_params%varnce_thl_1), &
                         pdf_params%stdev_eta_1, pdf_params%corr_w_rt_1, &
                         pdf_params%corr_w_thl_1 )

      pdf_params%corr_w_eta_2 &
      = calc_corr_eta_x( pdf_params%crt_2, pdf_params%cthl_2, &
                         sqrt(pdf_params%varnce_rt_2), sqrt(pdf_params%varnce_thl_2), &
                         pdf_params%stdev_eta_2, pdf_params%corr_w_rt_2, &
                         pdf_params%corr_w_thl_2 )

    end if

    
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
    ! 
    ! Reference:
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:anl_int_buoy_terms
    !
    ! Calculate the contributions to <w'rc'>, <w'^2 rc'>, <rt'rc'>, and
    ! <thl'rc'> from the 1st PDF component.
    call calc_xprcp_component( nz, ngrdcol,                                          & ! In
                               wm, rtm, thlm, um, vm, rcm,                           & ! In
                               pdf_params%w_1, pdf_params%rt_1,                      & ! In
                               pdf_params%thl_1, u_1, v_1,                           & ! In
                               pdf_params%varnce_w_1, pdf_params%chi_1,              & ! In
                               pdf_params%stdev_chi_1, pdf_params%stdev_eta_1,       & ! In
                               pdf_params%corr_w_chi_1, pdf_params%corr_chi_eta_1,   & ! In
                               pdf_params%crt_1, pdf_params%cthl_1,                  & ! In
                               pdf_params%rc_1, pdf_params%cloud_frac_1, iiPDF_type, & ! In
                               wprcp_contrib_comp_1, wp2rcp_contrib_comp_1,          & ! Out
                               rtprcp_contrib_comp_1, thlprcp_contrib_comp_1,        & ! Out
                               uprcp_contrib_comp_1, vprcp_contrib_comp_1 )            ! Out

    call calc_xprcp_component( nz, ngrdcol,                                          & ! In
                               wm, rtm, thlm, um, vm, rcm,                           & ! In
                               pdf_params%w_2, pdf_params%rt_2,                      & ! In
                               pdf_params%thl_2, u_2, v_2,                           & ! In
                               pdf_params%varnce_w_2, pdf_params%chi_2,              & ! In
                               pdf_params%stdev_chi_2, pdf_params%stdev_eta_2,       & ! In
                               pdf_params%corr_w_chi_2, pdf_params%corr_chi_eta_2,   & ! In
                               pdf_params%crt_2, pdf_params%cthl_2,                  & ! In
                               pdf_params%rc_2, pdf_params%cloud_frac_2, iiPDF_type, & ! In
                               wprcp_contrib_comp_2, wp2rcp_contrib_comp_2,          & ! Out
                               rtprcp_contrib_comp_2, thlprcp_contrib_comp_2,        & ! Out
                               uprcp_contrib_comp_2, vprcp_contrib_comp_2 )            ! Out


    ! Calculate rc_coef, which is the coefficient on <x'rc'> in the <x'thv'> equation.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        rc_coef(i,k) = Lv / ( exner(i,k) * Cp ) - ep2 * thv_ds(i,k)

        ! Calculate <w'rc'>, <w'^2 rc'>, <rt'rc'>, and <thl'rc'>.
        wprcp(i,k) = pdf_params%mixt_frac(i,k) * wprcp_contrib_comp_1(i,k) &
                     + ( one - pdf_params%mixt_frac(i,k) ) * wprcp_contrib_comp_2(i,k)

        wp2rcp(i,k) = pdf_params%mixt_frac(i,k) * wp2rcp_contrib_comp_1(i,k) &
                      + ( one - pdf_params%mixt_frac(i,k) ) * wp2rcp_contrib_comp_2(i,k)

        rtprcp(i,k) = pdf_params%mixt_frac(i,k) * rtprcp_contrib_comp_1(i,k) &
                      + ( one - pdf_params%mixt_frac(i,k) ) * rtprcp_contrib_comp_2(i,k)

        thlprcp(i,k) = pdf_params%mixt_frac(i,k) * thlprcp_contrib_comp_1(i,k) &
                       + ( one - pdf_params%mixt_frac(i,k) ) * thlprcp_contrib_comp_2(i,k)

        uprcp(i,k) = pdf_params%mixt_frac(i,k) * uprcp_contrib_comp_1(i,k) &
                     + ( one - pdf_params%mixt_frac(i,k) ) * uprcp_contrib_comp_2(i,k)

        vprcp(i,k) = pdf_params%mixt_frac(i,k) * vprcp_contrib_comp_1(i,k) &
                     + ( one - pdf_params%mixt_frac(i,k) ) * vprcp_contrib_comp_2(i,k)
      end do
    end do
    !$acc end parallel loop

    ! Calculate <w'thv'>, <w'^2 thv'>, <rt'thv'>, and <thl'thv'>.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        wpthvp(i,k)  = wpthlp(i,k)  + ep1 * thv_ds(i,k) * wprtp(i,k)   + rc_coef(i,k) * wprcp(i,k)
        wp2thvp(i,k) = wp2thlp(i,k) + ep1 * thv_ds(i,k) * wp2rtp(i,k)  + rc_coef(i,k) * wp2rcp(i,k)
        rtpthvp(i,k) = rtpthlp(i,k) + ep1 * thv_ds(i,k) * rtp2(i,k)    + rc_coef(i,k) * rtprcp(i,k)
        thlpthvp(i,k)= thlp2(i,k)   + ep1 * thv_ds(i,k) * rtpthlp(i,k) + rc_coef(i,k) * thlprcp(i,k)
      end do
    end do
    !$acc end parallel loop

    ! Add the precipitation loading term in the <x'thv'> equation.
    if ( l_liq_ice_loading_test ) then

       do hm_idx = 1, hydromet_dim, 1

          if ( l_mix_rat_hm(hm_idx) ) then
            !$acc parallel loop gang vector collapse(2) default(present)
            do k = 1, nz
              do i = 1, ngrdcol
                wp2thvp(i,k)  = wp2thvp(i,k)  - thv_ds(i,k) * wp2hmp(i,k,hm_idx)
                wpthvp(i,k)   = wpthvp(i,k)   - thv_ds(i,k) * wphydrometp(i,k,hm_idx)
                thlpthvp(i,k) = thlpthvp(i,k) - thv_ds(i,k) * thlphmp(i,k,hm_idx)
                rtpthvp(i,k)  = rtpthvp(i,k)  - thv_ds(i,k) * rtphmp(i,k,hm_idx)
              end do
            end do
            !$acc end parallel loop
          end if

       end do

    end if

    ! Account for subplume correlation of scalar, theta_v.
    ! See Eqs. A13, A8 from Larson et al. (2002) ``Small-scale...''
    !  where the ``scalar'' in this paper is w.
    if ( l_scalar_calc ) then

      !$acc parallel loop gang vector collapse(3) default(present)
      do sclr = 1, sclr_dim
        do k = 1, nz
          do i = 1, ngrdcol

            sclrprcp(i,k,sclr) &
            = pdf_params%mixt_frac(i,k) * ( ( sclr1(i,k,sclr) - sclrm(i,k,sclr) ) &
                                            * pdf_params%rc_1(i,k) ) &
              + ( one - pdf_params%mixt_frac(i,k) ) * ( ( sclr2(i,k,sclr) - sclrm(i,k,sclr) ) &
                                                        * pdf_params%rc_2(i,k) ) &
              + pdf_params%mixt_frac(i,k) * corr_sclr_rt_1(i,k,sclr) * pdf_params%crt_1(i,k) &
                * sqrt( varnce_sclr1(i,k,sclr) * pdf_params%varnce_rt_1(i,k) ) &
                * pdf_params%cloud_frac_1(i,k) &
              + ( one - pdf_params%mixt_frac(i,k) ) * corr_sclr_rt_2(i,k,sclr) &
                * pdf_params%crt_2(i,k) &
                * sqrt( varnce_sclr2(i,k,sclr) * pdf_params%varnce_rt_2(i,k) ) &
                * pdf_params%cloud_frac_2(i,k) &
              - pdf_params%mixt_frac(i,k) * corr_sclr_thl_1(i,k,sclr) * pdf_params%cthl_1(i,k) &
                * sqrt( varnce_sclr1(i,k,sclr) * pdf_params%varnce_thl_1(i,k) ) &
                * pdf_params%cloud_frac_1(i,k) &
              - ( one - pdf_params%mixt_frac(i,k) ) * corr_sclr_thl_2(i,k,sclr) &
                * pdf_params%cthl_2(i,k) &
                * sqrt( varnce_sclr2(i,k,sclr) * pdf_params%varnce_thl_2(i,k) ) &
                * pdf_params%cloud_frac_2(i,k)

            sclrpthvp(i,k,sclr) = sclrpthlp(i,k,sclr) + ep1*thv_ds(i,k)*sclrprtp(i,k,sclr) &
                             + rc_coef(i,k)*sclrprcp(i,k,sclr)

          end do
        end do
      end do ! sclr=1, sclr_dim
      !$acc end parallel loop

    end if ! l_scalar_calc


    ! Compute variance of liquid water mixing ratio.
    ! This is not needed for closure.  Statistical Analysis only.

#ifndef CLUBB_CAM
      ! if CLUBB is used in CAM we want this variable computed no matter what
      if ( var_on_stats_list( stats, "rcp2" ) ) then
#endif
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1,nz
      do i = 1, ngrdcol
        rcp2(i,k) = pdf_params%mixt_frac(i,k) &
                    * ( pdf_params%chi_1(i,k)*pdf_params%rc_1(i,k) &
                        + pdf_params%cloud_frac_1(i,k)*pdf_params%stdev_chi_1(i,k)**2 ) &
                    + ( one-pdf_params%mixt_frac(i,k) ) &
                      * ( pdf_params%chi_2(i,k)*pdf_params%rc_2(i,k) &
                          + pdf_params%cloud_frac_2(i,k)*pdf_params%stdev_chi_2(i,k)**2 ) &
                    - rcm(i,k)**2
        rcp2(i,k) = max( zero_threshold, rcp2(i,k) )

      end do
    end do
    !$acc end parallel loop
#ifndef CLUBB_CAM
      !  if CLUBB is used in CAM we want this variable computed no matter what
      end if
#endif

    if ( ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
           .or. iiPDF_type == iiPDF_new_hybrid ) &
         .and. ( var_on_stats_list( stats, "w_up_in_cloud" )       &
               .or. var_on_stats_list( stats, "w_down_in_cloud" )) ) then
                 
      call calc_w_up_in_cloud( nz, ngrdcol, &                                      ! In
                               pdf_params%mixt_frac, &                             ! In
                               pdf_params%cloud_frac_1, pdf_params%cloud_frac_2, & ! In
                               pdf_params%w_1, pdf_params%w_2, &                   ! In
                               pdf_params%varnce_w_1, pdf_params%varnce_w_2, &     ! In
                               w_up_in_cloud, w_down_in_cloud, &                   ! Out
                               cloudy_updraft_frac, cloudy_downdraft_frac )        ! Out

    else
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1,nz
        do i = 1, ngrdcol
          w_up_in_cloud(i,k) = zero
          w_down_in_cloud(i,k) = zero
          cloudy_updraft_frac(i,k) = zero
          cloudy_downdraft_frac(i,k) = zero
        end do
      end do
    end if

#ifdef TUNER

    !$acc update host( pdf_params, pdf_params%thl_1, pdf_params%thl_2 )

    ! Check the first levels (and first gridcolumn) for reasonable temperatures
    ! greater than 190K and less than 1000K
    ! This is necessary because for certain parameter sets we can get floating point errors
    do i=1, min( 10, size(pdf_params%thl_1(1,:)) )
        if ( pdf_params%thl_1(1,i) < 190. ) then
            write(fstderr, *) err_info%err_header(1)
            write(fstderr,*) "Fatal error: pdf_params%thl_1 =", pdf_params%thl_1(1,i), &
                             " < 190K at first grid column and grid level i = ", i
          ! Error in grid column 1 -> set 1st entry to clubb_fatal_error
          err_info%err_code(1) = clubb_fatal_error
        end if
        if ( pdf_params%thl_2(1,i) < 190. ) then
            write(fstderr, *) err_info%err_header(1)
            write(fstderr,*) "Fatal error: pdf_params%thl_2 =", pdf_params%thl_2(1,i), &
                             " < 190K at first grid column and grid level i = ", i
          ! Error in grid column 1 -> set 1st entry to clubb_fatal_error
          err_info%err_code(1) = clubb_fatal_error
        end if
        if ( pdf_params%thl_1(1,i) > 1000. ) then
            write(fstderr, *) err_info%err_header(1)
            write(fstderr,*) "Fatal error: pdf_params%thl_1 =", pdf_params%thl_1(1,i), &
                             " > 1000K at first grid column and grid level i = ", i
          ! Error in grid column 1 -> set 1st entry to clubb_fatal_error
          err_info%err_code(1) = clubb_fatal_error
        end if
        if ( pdf_params%thl_2(1,i) > 1000. ) then
            write(fstderr, *) err_info%err_header(1)
            write(fstderr,*) "Fatal error: pdf_params%thl_2 =", pdf_params%thl_2(1,i), &
                             " > 1000K at first grid column and grid level i = ", i
          ! Error in grid column 1 -> set 1st entry to clubb_fatal_error
          err_info%err_code(1) = clubb_fatal_error
        end if
    end do
    
    if ( any(err_info%err_code == clubb_fatal_error) ) return

#endif /*TUNER*/

    if ( clubb_at_least_debug_level_api( 2 ) ) then

      !$acc update host( wp4, wprtp2, wp2rtp, wpthlp2, wp2thlp, cloud_frac, &
      !$acc              rcm, wpthvp, wp2thvp, rtpthvp, thlpthvp, wprcp, wp2rcp, &
      !$acc              rtprcp, thlprcp, rcp2, wprtpthlp, wp2up, &
      !$acc              pdf_params%w_1, pdf_params%w_2, &
      !$acc              pdf_params%varnce_w_1, pdf_params%varnce_w_2, &
      !$acc              pdf_params%rt_1, pdf_params%rt_2, &
      !$acc              pdf_params%varnce_rt_1, pdf_params%varnce_rt_2,  &
      !$acc              pdf_params%thl_1, pdf_params%thl_2, &
      !$acc              pdf_params%varnce_thl_1, pdf_params%varnce_thl_2, &
      !$acc              pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2,  &
      !$acc              pdf_params%corr_w_thl_1, pdf_params%corr_w_thl_2, &
      !$acc              pdf_params%corr_rt_thl_1, pdf_params%corr_rt_thl_2,&
      !$acc              pdf_params%alpha_thl, pdf_params%alpha_rt, &
      !$acc              pdf_params%crt_1, pdf_params%crt_2, pdf_params%cthl_1, &
      !$acc              pdf_params%cthl_2, pdf_params%chi_1, &
      !$acc              pdf_params%chi_2, pdf_params%stdev_chi_1, &
      !$acc              pdf_params%stdev_chi_2, pdf_params%stdev_eta_1, &
      !$acc              pdf_params%stdev_eta_2, pdf_params%covar_chi_eta_1, &
      !$acc              pdf_params%covar_chi_eta_2, pdf_params%corr_w_chi_1, &
      !$acc              pdf_params%corr_w_chi_2, pdf_params%corr_w_eta_1, &
      !$acc              pdf_params%corr_w_eta_2, pdf_params%corr_chi_eta_1, &
      !$acc              pdf_params%corr_chi_eta_2, pdf_params%rsatl_1, &
      !$acc              pdf_params%rsatl_2, pdf_params%rc_1, pdf_params%rc_2, &
      !$acc              pdf_params%cloud_frac_1, pdf_params%cloud_frac_2,  &
      !$acc              pdf_params%mixt_frac, pdf_params%ice_supersat_frac_1, &
      !$acc              pdf_params%ice_supersat_frac_2 )

      !$acc update host( sclrpthvp, sclrprcp, wpsclrp2, wpsclrprtp, wpsclrpthlp, wp2sclrp ) &
      !$acc if ( sclr_dim > 0 )

      do i = 1, ngrdcol

        call pdf_closure_check( &
               nz, sclr_dim, &
               wp4(i,:), wprtp2(i,:), wp2rtp(i,:), wpthlp2(i,:), & ! intent(in)
               wp2thlp(i,:), cloud_frac(i,:), rcm(i,:), wpthvp(i,:), wp2thvp(i,:), & ! intent(in)
               wp2up(i,:), & ! intent(in)
               rtpthvp(i,:), thlpthvp(i,:), wprcp(i,:), wp2rcp(i,:), & ! intent(in)
               rtprcp(i,:), thlprcp(i,:), rcp2(i,:), wprtpthlp(i,:), & ! intent(in)
               pdf_params%crt_1(i,:), pdf_params%crt_2(i,:), & ! intent(in)
               pdf_params%cthl_1(i,:), pdf_params%cthl_2(i,:), & ! intent(in)
               pdf_params, & ! intent(in)
               sclrpthvp(i,:,:), sclrprcp(i,:,:), wpsclrp2(i,:,:), & ! intent(in)
               wpsclrprtp(i,:,:), wpsclrpthlp(i,:,:), wp2sclrp(i,:,:), & ! intent(in)
               stats,         & ! intent(in)
               err_info ) ! intent(inout)

        if ( err_info%err_code(i) == clubb_fatal_error ) then
          write(fstderr, *) err_info%err_header(i)
        endif
      end do
    end if

    ! Error Reporting
    ! Joshua Fasching February 2008
    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then

        !$acc update host( p_in_Pa, exner, thv_ds, wm, wp2, wp3, sigma_sqd_w, &
        !$acc              rtm, rtp2, wprtp, thlm, thlp2, wpthlp, rtpthlp, ice_supersat_frac )

        !$acc update host( sclrm, wpsclrp, sclrp2, sclrprtp, sclrpthlp ) if ( sclr_dim > 0 )

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

        write(fstderr,*) "Intent(out)"

        write(fstderr,*) "wp4 = ", wp4
        if ( l_explicit_turbulent_adv_xpyp .or. var_on_stats_list( stats, "wprtp2" ) ) then
          write(fstderr,*) "wprtp2 = ", wprtp2
        end if
        write(fstderr,*) "wp2rtp = ", wp2rtp
        if ( l_explicit_turbulent_adv_xpyp .or. var_on_stats_list( stats, "wpthlp2" ) ) then
          write(fstderr,*) "wpthlp2 = ", wpthlp2
        end if
        write(fstderr,*) "cloud_frac = ", cloud_frac
        write(fstderr,*) "ice_supersat_frac = ", ice_supersat_frac
        write(fstderr,*) "rcm = ", rcm
        write(fstderr,*) "wpthvp = ", wpthvp
        write(fstderr,*) "wp2thvp = ", wp2thvp
        write(fstderr,*) "wp2up = ", wp2up
        write(fstderr,*) "rtpthvp = ", rtpthvp
        write(fstderr,*) "thlpthvp = ", thlpthvp
        write(fstderr,*) "wprcp = ", wprcp
        write(fstderr,*) "wp2rcp = ", wp2rcp
        write(fstderr,*) "rtprcp = ", rtprcp
        write(fstderr,*) "thlprcp = ", thlprcp
#ifndef CLUBB_CAM
        ! if CLUBB is used in CAM we want this variable computed no matter what
        if ( var_on_stats_list( stats, "rcp2" ) ) then
#endif
          write(fstderr,*) "rcp2 = ", rcp2
#ifndef CLUBB_CAM
        end if
#endif
        if ( l_explicit_turbulent_adv_xpyp .or. var_on_stats_list( stats, "wprtpthlp" ) ) then
          write(fstderr,*) "wprtpthlp = ", wprtpthlp
        end if
        write(fstderr,*) "rcp2 = ", rcp2
        write(fstderr,*) "wprtpthlp = ", wprtpthlp
        write(fstderr,*) "pdf_params%w_1 = ", pdf_params%w_1
        write(fstderr,*) "pdf_params%w_2 = ", pdf_params%w_2
        write(fstderr,*) "pdf_params%varnce_w_1 = ", pdf_params%varnce_w_1
        write(fstderr,*) "pdf_params%varnce_w_2 = ", pdf_params%varnce_w_2
        write(fstderr,*) "pdf_params%rt_1 = ", pdf_params%rt_1
        write(fstderr,*) "pdf_params%rt_2 = ", pdf_params%rt_2
        write(fstderr,*) "pdf_params%varnce_rt_1 = ", pdf_params%varnce_rt_1
        write(fstderr,*) "pdf_params%varnce_rt_2 = ", pdf_params%varnce_rt_2
        write(fstderr,*) "pdf_params%thl_1 = ", pdf_params%thl_1
        write(fstderr,*) "pdf_params%thl_2 = ", pdf_params%thl_2
        write(fstderr,*) "pdf_params%varnce_thl_1 = ", pdf_params%varnce_thl_1
        write(fstderr,*) "pdf_params%varnce_thl_2 = ", pdf_params%varnce_thl_2
        write(fstderr,*) "pdf_params%corr_w_rt_1 = ", pdf_params%corr_w_rt_1
        write(fstderr,*) "pdf_params%corr_w_rt_2 = ", pdf_params%corr_w_rt_2
        write(fstderr,*) "pdf_params%corr_w_thl_1 = ", pdf_params%corr_w_thl_1
        write(fstderr,*) "pdf_params%corr_w_thl_2 = ", pdf_params%corr_w_thl_2
        write(fstderr,*) "pdf_params%corr_rt_thl_1 = ", pdf_params%corr_rt_thl_1
        write(fstderr,*) "pdf_params%corr_rt_thl_2 = ", pdf_params%corr_rt_thl_2
        write(fstderr,*) "pdf_params%alpha_thl = ", pdf_params%alpha_thl
        write(fstderr,*) "pdf_params%alpha_rt = ", pdf_params%alpha_rt
        write(fstderr,*) "pdf_params%crt_1 = ", pdf_params%crt_1
        write(fstderr,*) "pdf_params%crt_2 = ", pdf_params%crt_2
        write(fstderr,*) "pdf_params%cthl_1 = ", pdf_params%cthl_1
        write(fstderr,*) "pdf_params%cthl_2 = ", pdf_params%cthl_2
        write(fstderr,*) "pdf_params%chi_1 = ", pdf_params%chi_1
        write(fstderr,*) "pdf_params%chi_2 = ", pdf_params%chi_2
        write(fstderr,*) "pdf_params%stdev_chi_1 = ", pdf_params%stdev_chi_1
        write(fstderr,*) "pdf_params%stdev_chi_2 = ", pdf_params%stdev_chi_2
        write(fstderr,*) "pdf_params%stdev_eta_1 = ", pdf_params%stdev_eta_1
        write(fstderr,*) "pdf_params%stdev_eta_2 = ", pdf_params%stdev_eta_2
        write(fstderr,*) "pdf_params%covar_chi_eta_1 = ", &
                         pdf_params%covar_chi_eta_1
        write(fstderr,*) "pdf_params%covar_chi_eta_2 = ", &
                         pdf_params%covar_chi_eta_2
        write(fstderr,*) "pdf_params%corr_w_chi_1 = ", pdf_params%corr_w_chi_1
        write(fstderr,*) "pdf_params%corr_w_chi_2 = ", pdf_params%corr_w_chi_2
        write(fstderr,*) "pdf_params%corr_w_eta_1 = ", pdf_params%corr_w_eta_1
        write(fstderr,*) "pdf_params%corr_w_eta_2 = ", pdf_params%corr_w_eta_2
        write(fstderr,*) "pdf_params%corr_chi_eta_1 = ", &
                         pdf_params%corr_chi_eta_1
        write(fstderr,*) "pdf_params%corr_chi_eta_2 = ", &
                         pdf_params%corr_chi_eta_2
        write(fstderr,*) "pdf_params%rsatl_1 = ", pdf_params%rsatl_1
        write(fstderr,*) "pdf_params%rsatl_2 = ", pdf_params%rsatl_2
        write(fstderr,*) "pdf_params%rc_1 = ", pdf_params%rc_1
        write(fstderr,*) "pdf_params%rc_2 = ", pdf_params%rc_2
        write(fstderr,*) "pdf_params%cloud_frac_1 = ", pdf_params%cloud_frac_1
        write(fstderr,*) "pdf_params%cloud_frac_2 = ", pdf_params%cloud_frac_2
        write(fstderr,*) "pdf_params%mixt_frac = ", pdf_params%mixt_frac
        write(fstderr,*) "pdf_params%ice_supersat_frac_1 = ", &
                         pdf_params%ice_supersat_frac_1
        write(fstderr,*) "pdf_params%ice_supersat_frac_2 = ", &
                         pdf_params%ice_supersat_frac_2

        if ( sclr_dim > 0 )then
          write(fstderr,*) "sclrpthvp = ", sclrpthvp
          write(fstderr,*) "sclrprcp = ", sclrprcp
          write(fstderr,*) "wpsclrp2 = ", wpsclrp2
          write(fstderr,*) "wpsclrprtp = ", wpsclrprtp
          write(fstderr,*) "wpsclrpthlp = ", wpsclrpthlp
          write(fstderr,*) "wp2sclrp = ", wp2sclrp
        end if

        return

      end if ! Fatal error
          
      do i = 1, ngrdcol

        ! Error check pdf parameters and moments to ensure consistency
        if ( iiPDF_type == iiPDF_3D_Luhar ) then

          ! Means
          wm_clubb_pdf(i,:) = pdf_params%mixt_frac(i,:) * pdf_params%w_1(i,:) &
                         + ( one - pdf_params%mixt_frac(i,:) ) * pdf_params%w_2(i,:)

          do k = 1, nz, 1
             if ( abs( ( wm_clubb_pdf(i,k) - wm(i,k) ) &
                       / max( wm(i,k), eps ) ) > .05_core_rknd ) then
                write(fstderr,*) "wm error at thlm = ", thlm(i,k), &
                                 ( ( wm_clubb_pdf(i,k) - wm(i,k) ) &
                                   / max( wm(i,k), eps ) )
             end if
          end do ! k = 1, nz, 1

          rtm_clubb_pdf(i,:) = pdf_params%mixt_frac(i,:) * pdf_params%rt_1(i,:) &
                          + ( one - pdf_params%mixt_frac(i,:) ) * pdf_params%rt_2(i,:)

          do k = 1, nz, 1
             if ( abs( ( rtm_clubb_pdf(i,k) - rtm(i,k) ) &
                       / max( rtm(i,k), eps ) ) > .05_core_rknd ) then
                write(fstderr,*) "rtm error at thlm = ", thlm(i,k), &
                                 ( ( rtm_clubb_pdf(i,k) - rtm(i,k) ) &
                                   / max( rtm(i,k), eps ) )
             end if
          end do ! k = 1, nz, 1

          thlm_clubb_pdf(i,:) = pdf_params%mixt_frac(i,:) * pdf_params%thl_1(i,:) &
                           + ( one - pdf_params%mixt_frac(i,:) ) * pdf_params%thl_2(i,:)

          do k = 1, nz, 1
             if ( abs( ( thlm_clubb_pdf(i,k) - thlm(i,k) ) / thlm(i,k) ) &
                  > .05_core_rknd ) then
                write(fstderr,*) "thlm error at thlm = ", thlm(i,k), &
                                 ( ( thlm_clubb_pdf(i,k) - thlm(i,k) ) / thlm(i,k) )
             end if
          end do ! k = 1, nz, 1

          ! Variances
          wp2_clubb_pdf(i,:) = pdf_params%mixt_frac(i,:) &
                          * ( ( pdf_params%w_1(i,:) - wm(i,:) )**2 + pdf_params%varnce_w_1(i,:) ) &
                          + ( one - pdf_params%mixt_frac(i,:) ) &
                            * ( ( pdf_params%w_2(i,:) - wm(i,:) )**2 + pdf_params%varnce_w_2(i,:) )

          do k = 1, nz, 1
             if ( wp2(i,k) > w_tol**2 ) then
                if ( abs( ( wp2_clubb_pdf(i,k) - wp2(i,k) ) / wp2(i,k) ) &
                     > .05_core_rknd ) then
                   write(fstderr,*) "wp2 error at thlm = ", thlm(i,k), &
                                    ( ( wp2_clubb_pdf(i,k) - wp2(i,k) ) / wp2(i,k) )
                end if
             end if
          end do ! k = 1, nz, 1

          rtp2_clubb_pdf(i,:) &
          = pdf_params%mixt_frac(i,:) &
            * ( ( pdf_params%rt_1(i,:) - rtm(i,:) )**2 + pdf_params%varnce_rt_1(i,:) ) &
            + ( one - pdf_params%mixt_frac(i,:) ) &
              * ( ( pdf_params%rt_2(i,:) - rtm(i,:) )**2 + pdf_params%varnce_rt_2(i,:) )

          do k = 1, nz, 1
             if ( rtp2(i,k) > rt_tol**2 ) then
                if ( abs( ( rtp2_clubb_pdf(i,k) - rtp2(i,k) ) / rtp2(i,k) ) &
                     > .05_core_rknd ) then
                   write(fstderr,*) "rtp2 error at thlm = ", thlm(i,k), &
                   "Error = ", ( ( rtp2_clubb_pdf(i,k) - rtp2(i,k) ) / rtp2(i,k) )
                end if
             end if
          end do ! k = 1, nz, 1

          thlp2_clubb_pdf(i,:) &
          = pdf_params%mixt_frac(i,:) &
            * ( ( pdf_params%thl_1(i,:) - thlm(i,:) )**2 + pdf_params%varnce_thl_1(i,:) ) &
            + ( one - pdf_params%mixt_frac(i,:) ) &
              * ( ( pdf_params%thl_2(i,:) - thlm(i,:) )**2 + pdf_params%varnce_thl_2(i,:) )

          do k = 1, nz, 1
             if( thlp2(i,k) > thl_tol**2 ) then
                if ( abs( ( thlp2_clubb_pdf(i,k) - thlp2(i,k) ) / thlp2(i,k) ) &
                     > .05_core_rknd ) then
                   write(fstderr,*) "thlp2 error at thlm = ", thlm(i,k), &
                   "Error = ", ( ( thlp2_clubb_pdf(i,k) - thlp2(i,k) ) / thlp2(i,k) )
                end if
             end if
          end do ! k = 1, nz, 1

          ! Third order moments
          wp3_clubb_pdf(i,:) &
          = pdf_params%mixt_frac(i,:) * ( pdf_params%w_1(i,:) - wm(i,:) ) &
            * ( ( pdf_params%w_1(i,:) - wm(i,:) )**2 + three * pdf_params%varnce_w_1(i,:) ) &
            + ( one - pdf_params%mixt_frac(i,:) ) * ( pdf_params%w_2(i,:) - wm(i,:) ) &
              * ( ( pdf_params%w_2(i,:) - wm(i,:) )**2 + three * pdf_params%varnce_w_2(i,:) )

          rtp3_clubb_pdf(i,:) &
          = pdf_params%mixt_frac(i,:) * ( pdf_params%rt_1(i,:) - rtm(i,:) ) &
            * ( ( pdf_params%rt_1(i,:) - rtm(i,:) )**2 + three * pdf_params%varnce_rt_1(i,:) ) &
            + ( one - pdf_params%mixt_frac(i,:) ) * ( pdf_params%rt_2(i,:) - rtm(i,:) ) &
              * ( ( pdf_params%rt_2(i,:) - rtm(i,:) )**2 + three * pdf_params%varnce_rt_2(i,:) )

          thlp3_clubb_pdf(i,:) &
          = pdf_params%mixt_frac(i,:) * ( pdf_params%thl_1(i,:) - thlm(i,:) ) &
            * ( ( pdf_params%thl_1(i,:) - thlm(i,:) )**2 + three * pdf_params%varnce_thl_1(i,:) ) &
            + ( one - pdf_params%mixt_frac(i,:) ) * ( pdf_params%thl_2(i,:) - thlm(i,:) ) &
              * ( ( pdf_params%thl_2(i,:) - thlm(i,:) )**2 + three * pdf_params%varnce_thl_2(i,:) )

          ! Skewness
          Skw_denom_coef = clubb_params(i,iSkw_denom_coef)

          Skw_clubb_pdf(i,:) &
          = wp3_clubb_pdf(i,:) &
            / ( wp2_clubb_pdf(i,:) + Skw_denom_coef * w_tol**2 )**1.5_core_rknd

          do k = 1, nz, 1
             if ( Skw(i,k) > .05_core_rknd ) then
                if( abs( ( Skw_clubb_pdf(i,k) - Skw(i,k) ) / Skw(i,k) ) &
                    > .25_core_rknd ) then
                   write(fstderr,*) "Skw error at thlm = ", thlm(i,k), &
                   "Error = ", ( ( Skw_clubb_pdf(i,k) - Skw(i,k) ) / Skw(i,k) ), &
                   Skw_clubb_pdf(i,k), Skw(i,k)
                end if
             end if
          end do ! k = 1, nz, 1

          Skrt_clubb_pdf(i,:) &
          = rtp3_clubb_pdf(i,:) &
            / ( rtp2_clubb_pdf(i,:) + Skw_denom_coef * rt_tol**2 )**1.5_core_rknd

          do k = 1, nz, 1
             if ( Skrt(i,k) > .05_core_rknd ) then
                if( abs( ( Skrt_clubb_pdf(i,k) - Skrt(i,k) ) / Skrt(i,k) ) &
                    > .25_core_rknd ) then
                   write(fstderr,*) "Skrt error at thlm = ", thlm(i,k), &
                   "Error = ", ( ( Skrt_clubb_pdf(i,k) - Skrt(i,k) ) / Skrt(i,k) ), &
                   Skrt_clubb_pdf(i,k), Skrt(i,k)
                end if
             end if
          end do ! k = 1, nz, 1

          Skthl_clubb_pdf(i,:) &
          = thlp3_clubb_pdf(i,:) &
            / ( thlp2_clubb_pdf(i,:) + Skw_denom_coef * thl_tol**2 )**1.5_core_rknd

          do k = 1, nz, 1
             if ( Skthl(i,k) > .05_core_rknd ) then
                if ( abs( ( Skthl_clubb_pdf(i,k) - Skthl(i,k) ) / Skthl(i,k) ) &
                     > .25_core_rknd ) then
                   write(fstderr,*) "Skthl error at thlm = ", thlm(i,k), &
                   "Error = ", ( ( Skthl_clubb_pdf(i,k) - Skthl(i,k) ) / Skthl(i,k) ), &
                   Skthl_clubb_pdf(i,k), Skthl(i,k)
                end if
             end if
          end do ! k = 1, nz, 1

        end if ! iiPDF_type == iiPDF_3D_Luhar
        
      end do

    end if ! clubb_at_least_debug_level_api

    !$acc exit data delete( u_1, u_2, varnce_u_1, varnce_u_2, v_1, v_2, &
    !$acc                   varnce_v_1, varnce_v_2, alpha_u, alpha_v, &
    !$acc                   corr_u_w_1, corr_u_w_2, corr_v_w_1, corr_v_w_2, &
    !$acc                   tl1, tl2, sqrt_wp2, Skthl, &
    !$acc                   Skrt, Sku, Skv, wprcp_contrib_comp_1, wprcp_contrib_comp_2, &
    !$acc                   wp2rcp_contrib_comp_1, wp2rcp_contrib_comp_2, &
    !$acc                   rtprcp_contrib_comp_1, rtprcp_contrib_comp_2, &
    !$acc                   thlprcp_contrib_comp_1, thlprcp_contrib_comp_2, &
    !$acc                   uprcp_contrib_comp_1, uprcp_contrib_comp_2, &
    !$acc                   vprcp_contrib_comp_1, vprcp_contrib_comp_2, &
    !$acc                   rc_1_ice, rc_2_ice, rsatl_1, rsatl_2 )

    !$acc exit data if( sclr_dim > 0 ) &
    !$acc           delete( sclr1, sclr2, varnce_sclr1, varnce_sclr2, &
    !$acc                   alpha_sclr, corr_sclr_thl_1, corr_sclr_thl_2, &
    !$acc                   corr_sclr_rt_1, corr_sclr_rt_2, corr_w_sclr_1, &
    !$acc                   corr_w_sclr_2, Sksclr )

    return

  end subroutine pdf_closure

  !===============================================================================================
  subroutine transform_pdf_chi_eta_component( nz, ngrdcol, &
                                              tl, rsatl, rt, exner,     & ! In
                                              varnce_thl, varnce_rt,    & ! In
                                              corr_rt_thl, chi,         & ! In
                                              crt, cthl,                & ! Out
                                              stdev_chi, stdev_eta,     & ! Out
                                              covar_chi_eta,            & ! Out
                                              corr_chi_eta )              ! Out

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        one, two, &
        ep, Lv, Rd, Cp, &
        chi_tol, &
        eta_tol

    use pdf_utilities, only: &
        smooth_corr_quotient

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! ----------- Input Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      tl, &
      rsatl, &
      rt, &
      varnce_thl, &
      varnce_rt, &
      corr_rt_thl, &
      exner

    ! ----------- Output Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      chi, &            ! s from Lewellen and Yoh 1993 (LY) eqn. 1
      crt, &            ! Coefficients for s'
      cthl, &           ! Coefficients for s'
      stdev_chi, &      ! Standard deviation of chi for each component.
      stdev_eta, &      ! Standard deviation of eta for each component.
      covar_chi_eta, &  ! Covariance of chi and eta for each component.
      corr_chi_eta      ! Correlation of chi and eta for each component.

    ! ----------- Local Variables -----------
    real( kind = core_rknd ) :: &
      varnce_rt_term, &
      corr_rt_thl_term, &
      varnce_thl_term, &
      varnce_chi, &
      varnce_eta, &
      beta, &
      invrs_beta_rsatl_p1

    real( kind = core_rknd ), parameter :: &
      denom_thresh = chi_tol*eta_tol, &
      Cp_on_Lv = Cp / Lv

    real( kind = core_rknd ), dimension(ngrdcol, nz) :: &
      denominator

    ! Loop variable
    integer :: k, i

    ! ----------- Begin Code -----------

    !$acc enter data create( denominator )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        ! SD's beta (eqn. 8)
        beta = ep * Lv**2 / ( Rd * Cp * tl(i,k)**2 )

        invrs_beta_rsatl_p1 = one / ( one + beta * rsatl(i,k) )

        ! s from Lewellen and Yoh 1993 (LY) eqn. 1
        chi(i,k) = ( rt(i,k) - rsatl(i,k) ) * invrs_beta_rsatl_p1

        ! For each normal distribution in the sum of two normal distributions,
        ! s' = crt * rt'  +  cthl * thl';
        ! therefore, x's' = crt * x'rt'  +  cthl * x'thl'.
        ! Larson et al. May, 2001.
        crt(i,k)  = invrs_beta_rsatl_p1
        cthl(i,k) = ( one + beta * rt(i,k) ) * invrs_beta_rsatl_p1**2 &
                    * Cp_on_Lv * beta * rsatl(i,k) * exner(i,k)

      end do
    end do
    !$acc end parallel loop

    ! Calculate covariance, correlation, and standard deviation of
    ! chi and eta for each component
    ! Include subplume correlation of qt, thl
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        varnce_rt_term = crt(i,k)**2 * varnce_rt(i,k)
        varnce_thl_term = cthl(i,k)**2 * varnce_thl(i,k)

        covar_chi_eta(i,k) = varnce_rt_term - varnce_thl_term

        corr_rt_thl_term = two * corr_rt_thl(i,k) * crt(i,k) * cthl(i,k) &
                           * sqrt( varnce_rt(i,k) * varnce_thl(i,k) )

        varnce_chi = varnce_rt_term - corr_rt_thl_term + varnce_thl_term
        varnce_eta = varnce_rt_term + corr_rt_thl_term + varnce_thl_term

        stdev_chi(i,k) = sqrt( varnce_chi )
        stdev_eta(i,k) = sqrt( varnce_eta )

        denominator(i,k) = stdev_chi(i,k) * stdev_eta(i,k)

      end do
    end do
    !$acc end parallel loop

    call smooth_corr_quotient( ngrdcol, nz, covar_chi_eta, denominator, denom_thresh, corr_chi_eta )

    !$acc exit data delete( denominator )

    return

  end subroutine transform_pdf_chi_eta_component

  !=============================================================================
  subroutine calc_wp4_pdf( nz, ngrdcol, &
                           wm, w_1, w_2, &
                           varnce_w_1, varnce_w_2, &
                           mixt_frac, &
                           wp4 )

    ! Description:
    ! Calculates <w'^4> by integrating over the PDF of w.  The integral is:
    !
    ! <w'^4> = INT(-inf:inf) ( w - <w> )^4 P(w) dw;
    !
    ! where <w> is the overall mean of w and P(w) is a two-component normal
    ! distribution of w.  The integrated equation is:
    !
    ! <w'^4> = mixt_frac * ( 3 * sigma_w_1^4
    !                        + 6 * ( mu_w_1 - <w> )^2 * sigma_w_1^2
    !                        + ( mu_w_1 - <w> )^4 )
    !          + ( 1 - mixt_frac ) * ( 3 * sigma_w_2^4
    !                                  + 6 * ( mu_w_2 - <w> )^2 * sigma_w_2^2
    !                                  + ( mu_w_2 - <w> )^4 );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, sigma_w_1 is the standard deviation of w in
    ! the 1st PDF component, sigma_w_2 is the standard deviation of w in the 2nd
    ! PDF component, and mixt_frac is the mixture fraction, which is the weight
    ! of the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        six,   & ! Variable(s)
        three, &
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wm,         & ! Mean of w (overall)                           [m/s]
      w_1,        & ! Mean of w (1st PDF component)                 [m/s]
      w_2,        & ! Mean of w (2nd PDF component)                 [m/s]
      varnce_w_1, & ! Variance of w (1st PDF component)             [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)             [m^2/s^2]
      mixt_frac     ! Mixture fraction                              [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      wp4    ! <w'^4>                   [m^4/s^4]

    ! Local Variables
    integer :: i, k

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        ! Calculate <w'^4> by integrating over the PDF.
        wp4(i,k) = mixt_frac(i,k) * ( three * varnce_w_1(i,k)**2 &
                            + six * ( ( w_1(i,k) - wm(i,k) )**2 ) * varnce_w_1(i,k) &
                            + ( w_1(i,k) - wm(i,k) )**4 ) &
                   + ( one - mixt_frac(i,k) ) * ( three * varnce_w_2(i,k)**2 &
                                        + six * ( (w_2(i,k) - wm(i,k) )**2 )*varnce_w_2(i,k) &
                                        + ( w_2(i,k) - wm(i,k) )**4 )
      end do
    end do
    !$acc end parallel loop

    return

  end subroutine calc_wp4_pdf

  !=============================================================================
  subroutine calc_wp2xp2_pdf( nz, ngrdcol,             &
                              wm, xm, w_1,             &
                              w_2, x_1, x_2,           &
                              varnce_w_1, varnce_w_2,  &
                              varnce_x_1, varnce_x_2,  &
                              corr_w_x_1, corr_w_x_2,  &
                              mixt_frac, &
                              wp2xp2 )

    ! Description:
    ! Calculates <w'^2x'^2> by integrating over the PDF of w and x.  The
    ! integral
    ! is:
    !
    ! <w'^2x'^2>
    ! = INT(-inf:inf) INT(-inf:inf) ( w - <w> )^2 ( x - <x> )^2 P(w,x) dx dw;
    !
    ! where <w> is the overall mean of w, <x> is the overall mean of x, and
    ! P(w,x) is a two-component bivariate normal distribution of w and x.  The
    ! integrated equation is:
    !
    ! <w'^2x'^2>
    !   = mixt_frac
    !      * ( ( mu_w_1 - <w> )**2 * ( ( mu_x_1 - <x> )**2 + sigma_x_1^2 )
    !      + four * corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_x_1 - <x> ) * (
    !      mu_w_1 - <w> )
    !      + ( ( mu_x_1 - <x> )**2 + ( 1 + 2*corr_w_x_1**2 ) * sigma_x_1^2 ) *
    !      sigma_w_1^2 )
    !    + ( one - mixt_frac )
    !      * ( ( mu_w_2 - <w> )**2 * ( ( mu_x_2 - <x> )**2 + sigma_x_2^2 )
    !      + four * corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_x_2 - <x> ) * (
    !      mu_w_2 - <w> )
    !      + ( ( mu_x_2 - <x> )**2 + ( 1 + 2*corr_w_x_2**2 ) * sigma_x_2^2 ) *
    !      sigma_w_2^2 )
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    ! component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
    ! the standard deviation of w in the 1st PDF component, sigma_w_2 is the
    ! standard deviation of w in the 2nd PDF component, sigma_x_1 is the
    ! standard deviation of x in the 1st PDF component, sigma_x_2 is the
    ! standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
    ! correlation of w and x in the 1st PDF component, corr_w_x_2 is the
    ! correlation of w and x in the 2nd PDF component, and mixt_frac is the
    ! mixture fraction, which is the weight of the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,   & ! Variable(s)
        four

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wm,         & ! Mean of w (overall)                       [m/s]
      xm,         & ! Mean of x (overall)                       [units vary]
      w_1,        & ! Mean of w (1st PDF component)             [m/s]
      w_2,        & ! Mean of w (2nd PDF component)             [m/s]
      x_1,        & ! Mean of x (1st PDF component)             [units vary]
      x_2,        & ! Mean of x (2nd PDF component)             [units vary]
      varnce_w_1, & ! Variance of w (1st PDF component)         [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)         [m^2/s^2]
      varnce_x_1, & ! Variance of x (1st PDF component)         [(units vary)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)         [(units vary)^2]
      corr_w_x_1, & ! Correlation of w and x (1st PDF comp.)    [-]
      corr_w_x_2, & ! Correlation of w and x (2nd PDF comp.)    [-]
      mixt_frac     ! Mixture fraction                          [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      wp2xp2        ! <w'^2x'^2>                   [m^2/s^2 (units vary)^2]

    ! Local Variable
    integer :: i, k

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        ! Calculate <w'x'^2> by integrating over the PDF.
        wp2xp2(i,k) = mixt_frac(i,k) &
                * ( ( w_1(i,k) - wm(i,k) )**2 * ( ( x_1(i,k) - xm(i,k) )**2 + varnce_x_1(i,k) ) &
                + four * corr_w_x_1(i,k) * sqrt( varnce_w_1(i,k) * varnce_x_1(i,k) ) &
                                         * ( x_1(i,k) - xm(i,k) ) * ( w_1(i,k) - wm(i,k) ) &
                + ( ( x_1(i,k) - xm(i,k) )**2 &
                + ( 1 + 2*corr_w_x_1(i,k)**2 ) * varnce_x_1(i,k) ) * varnce_w_1(i,k) ) &
                + ( one - mixt_frac(i,k) ) &
                * ( ( w_2(i,k) - wm(i,k) )**2 * ( ( x_2(i,k) - xm(i,k) )**2 + varnce_x_2(i,k) ) &
                + four * corr_w_x_2(i,k) * sqrt( varnce_w_2(i,k) * varnce_x_2(i,k) ) &
                                         * ( x_2(i,k) - xm(i,k) ) * ( w_2(i,k) - wm(i,k) ) &
                + ( ( x_2(i,k) - xm(i,k) )**2 &
                + ( 1 + 2*corr_w_x_2(i,k)**2 ) * varnce_x_2(i,k) ) * varnce_w_2(i,k) )
      end do
    end do
    !$acc end parallel loop

    return

  end subroutine calc_wp2xp2_pdf

  !=============================================================================
  subroutine calc_wp2xp_pdf( nz, ngrdcol,             &
                             wm, xm, w_1, w_2,        &
                             x_1, x_2,                &
                             varnce_w_1, varnce_w_2,  &
                             varnce_x_1, varnce_x_2,  &
                             corr_w_x_1, corr_w_x_2,  &
                             mixt_frac,               &
                             wp2xp )

    ! Description:
    ! Calculates <w'^2 x'> by integrating over the PDF of w and x.  The integral
    ! is:
    !
    ! <w'^2 x'>
    ! = INT(-inf:inf) INT(-inf:inf) ( w - <w> )^2 ( x - <x> ) P(w,x) dx dw;
    !
    ! where <w> is the overall mean of w, <x> is the overall mean of x, and
    ! P(w,x) is a two-component bivariate normal distribution of w and x.  The
    ! integrated equation is:
    !
    ! <w'^2 x'>
    ! = mixt_frac * ( ( mu_x_1 - <x> ) * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
    !                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
    !                   * ( mu_w_1 - <w> ) )
    !   + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )
    !                           * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 )
    !                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
    !                             * ( mu_w_2 - <w> ) );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    ! component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
    ! the standard deviation of w in the 1st PDF component, sigma_w_2 is the
    ! standard deviation of w in the 2nd PDF component, sigma_x_1 is the
    ! standard deviation of x in the 1st PDF component, sigma_x_2 is the
    ! standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
    ! correlation of w and x in the 1st PDF component, corr_w_x_2 is the
    ! correlation of w and x in the 2nd PDF component, and mixt_frac is the
    ! mixture fraction, which is the weight of the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,   & ! Variable(s)
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wm,         & ! Mean of w (overall)                       [m/s]
      xm,         & ! Mean of x (overall)                       [units vary]
      w_1,        & ! Mean of w (1st PDF component)             [m/s]
      w_2,        & ! Mean of w (2nd PDF component)             [m/s]
      x_1,        & ! Mean of x (1st PDF component)             [units vary]
      x_2,        & ! Mean of x (2nd PDF component)             [units vary]
      varnce_w_1, & ! Variance of w (1st PDF component)         [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)         [m^2/s^2]
      varnce_x_1, & ! Variance of x (1st PDF component)         [(units vary)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)         [(units vary)^2]
      corr_w_x_1, & ! Correlation of w and x (1st PDF comp.)    [-]
      corr_w_x_2, & ! Correlation of w and x (2nd PDF comp.)    [-]
      mixt_frac     ! Mixture fraction                          [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      wp2xp    ! <w'^2 x'>                   [m^2/s^2 (units vary)]

    ! Local Variables
    integer :: i, k


    ! Calculate <w'^2 x'> by integrating over the PDF.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        
        wp2xp(i,k)  = mixt_frac(i,k) &
                   * ( ( ( w_1(i,k) - wm(i,k) )**2 + varnce_w_1(i,k) ) * ( x_1(i,k) - xm(i,k) ) &
                       + two * corr_w_x_1(i,k) * sqrt( varnce_w_1(i,k) * varnce_x_1(i,k) ) &
                         * ( w_1(i,k) - wm(i,k) ) ) &
                   + ( one - mixt_frac(i,k) ) &
                     * ( ( ( w_2(i,k) - wm(i,k) )**2 + varnce_w_2(i,k) ) * ( x_2(i,k) - xm(i,k) ) &
                         + two * corr_w_x_2(i,k) * sqrt( varnce_w_2(i,k) * varnce_x_2(i,k) ) &
                           * ( w_2(i,k) - wm(i,k) ) )
      end do
    end do
    !$acc end parallel loop

    return

  end subroutine calc_wp2xp_pdf

  !=============================================================================
  subroutine calc_wpxp2_pdf( nz, ngrdcol,             &
                             wm, xm, w_1,             &
                             w_2, x_1, x_2,           &
                             varnce_w_1, varnce_w_2,  &
                             varnce_x_1, varnce_x_2,  &
                             corr_w_x_1, corr_w_x_2,  &
                             mixt_frac,               &
                             wpxp2 )

    ! Description:
    ! Calculates <w'x'^2> by integrating over the PDF of w and x.  The integral
    ! is:
    !
    ! <w'x'^2>
    ! = INT(-inf:inf) INT(-inf:inf) ( w - <w> ) ( x - <x> )^2 P(w,x) dx dw;
    !
    ! where <w> is the overall mean of w, <x> is the overall mean of x, and
    ! P(w,x) is a two-component bivariate normal distribution of w and x.  The
    ! integrated equation is:
    !
    ! <w'x'^2>
    ! = mixt_frac * ( ( mu_w_1 - <w> ) * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
    !                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
    !                   * ( mu_x_1 - <x> ) )
    !   + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )
    !                           * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 )
    !                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
    !                             * ( mu_x_2 - <x> ) );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    ! component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
    ! the standard deviation of w in the 1st PDF component, sigma_w_2 is the
    ! standard deviation of w in the 2nd PDF component, sigma_x_1 is the
    ! standard deviation of x in the 1st PDF component, sigma_x_2 is the
    ! standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
    ! correlation of w and x in the 1st PDF component, corr_w_x_2 is the
    ! correlation of w and x in the 2nd PDF component, and mixt_frac is the
    ! mixture fraction, which is the weight of the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        two,   & ! Variable(s)
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wm,         & ! Mean of w (overall)                       [m/s]
      xm,         & ! Mean of x (overall)                       [units vary]
      w_1,        & ! Mean of w (1st PDF component)             [m/s]
      w_2,        & ! Mean of w (2nd PDF component)             [m/s]
      x_1,        & ! Mean of x (1st PDF component)             [units vary]
      x_2,        & ! Mean of x (2nd PDF component)             [units vary]
      varnce_w_1, & ! Variance of w (1st PDF component)         [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)         [m^2/s^2]
      varnce_x_1, & ! Variance of x (1st PDF component)         [(units vary)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)         [(units vary)^2]
      corr_w_x_1, & ! Correlation of w and x (1st PDF comp.)    [-]
      corr_w_x_2, & ! Correlation of w and x (2nd PDF comp.)    [-]
      mixt_frac     ! Mixture fraction                          [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      wpxp2    ! <w'x'^2>                   [m/s (units vary)^2]
      
    ! Local Variables
    integer :: i, k

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        ! Calculate <w'x'^2> by integrating over the PDF.
        wpxp2(i,k) = mixt_frac(i,k) &
                * ( ( w_1(i,k) - wm(i,k) ) * ( ( x_1(i,k) - xm(i,k) )**2 + varnce_x_1(i,k) ) &
                    + two * corr_w_x_1(i,k) * sqrt( varnce_w_1(i,k) * varnce_x_1(i,k) ) &
                      * ( x_1(i,k) - xm(i,k) ) ) &
                + ( one - mixt_frac(i,k) ) &
                  * ( ( w_2(i,k) - wm(i,k) ) * ( ( x_2(i,k) - xm(i,k) )**2 + varnce_x_2(i,k) ) &
                      + two * corr_w_x_2(i,k) * sqrt( varnce_w_2(i,k) * varnce_x_2(i,k) ) &
                        * ( x_2(i,k) - xm(i,k) ) )
      end do
    end do
    !$acc end parallel loop

    return

  end subroutine calc_wpxp2_pdf

  !=============================================================================
  subroutine calc_wpxpyp_pdf( nz, ngrdcol, &
                              wm, xm, ym, w_1, w_2,   &
                              x_1, x_2,               &
                              y_1, y_2,               &
                              varnce_w_1, varnce_w_2, &
                              varnce_x_1, varnce_x_2, &
                              varnce_y_1, varnce_y_2, &
                              corr_w_x_1, corr_w_x_2, &
                              corr_w_y_1, corr_w_y_2, &
                              corr_x_y_1, corr_x_y_2, &
                              mixt_frac, &
                              wpxpyp )

    ! Description:
    ! Calculates <w'x'y'> by integrating over the PDF of w, x, and y.  The
    ! integral is:
    !
    ! <w'x'y'>
    ! = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> ) ( x - <x> ) ( y - <y> ) P(w,x,y) dy dx dw;
    !
    ! where <w> is the overall mean of w, <x> is the overall mean of x, <y> is
    ! the overall mean of y, and P(w,x,y) is a two-component trivariate normal
    ! distribution of w, x, and y.  The integrated equation is:
    !
    ! <w'x'y'>
    ! = mixt_frac 
    !   * ( ( mu_w_1 - <w> ) * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
    !       + corr_x_y_1 * sigma_x_1 * sigma_y_1 * ( mu_w_1 - <w> )
    !       + corr_w_y_1 * sigma_w_1 * sigma_y_1 * ( mu_x_1 - <x> )
    !       + corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_y_1 - <y> ) )
    !   + ( 1 - mixt_frac )
    !     * ( ( mu_w_2 - <w> ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
    !         + corr_x_y_2 * sigma_x_2 * sigma_y_2 * ( mu_w_2 - <w> )
    !         + corr_w_y_2 * sigma_w_2 * sigma_y_2 * ( mu_x_2 - <x> )
    !         + corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_y_2 - <y> ) );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    ! component, mu_x_2 is the mean of x in the 2nd PDF component, mu_y_1 is the
    ! mean of y in the 1st PDF component, mu_y_2 is the mean of y in the 2nd PDF
    ! component, sigma_w_1 is the standard deviation of w in the 1st PDF
    ! component, sigma_w_2 is the standard deviation of w in the 2nd PDF
    ! component, sigma_x_1 is the standard deviation of x in the 1st PDF
    ! component, sigma_x_2 is the standard deviation of x in the 2nd PDF
    ! component, sigma_y_1 is the standard deviation of y in the 1st PDF
    ! component, sigma_y_2 is the standard deviation of y in the 2nd PDF
    ! component, corr_w_x_1 is the correlation of w and x in the 1st PDF
    ! component, corr_w_x_2 is the correlation of w and x in the 2nd PDF
    ! component, corr_w_y_1 is the correlation of w and y in the 1st PDF
    ! component, corr_w_y_2 is the correlation of w and y in the 2nd PDF
    ! component, corr_x_y_1 is the correlation of x and y in the 1st PDF
    ! component, corr_x_y_2 is the correlation of x and y in the 2nd PDF
    ! component, and mixt_frac is the mixture fraction, which is the weight of
    ! the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one    ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wm,         & ! Mean of w (overall)                          [m/s]
      xm,         & ! Mean of x (overall)                          [x units]
      ym,         & ! Mean of y (overall)                          [y units]
      w_1,        & ! Mean of w (1st PDF component)                [m/s]
      w_2,        & ! Mean of w (2nd PDF component)                [m/s]
      x_1,        & ! Mean of x (1st PDF component)                [x units]
      x_2,        & ! Mean of x (2nd PDF component)                [x units]
      y_1,        & ! Mean of y (1st PDF component)                [y units]
      y_2,        & ! Mean of y (2nd PDF component)                [y units]
      varnce_w_1, & ! Variance of w (1st PDF component)            [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)            [m^2/s^2]
      varnce_x_1, & ! Variance of x (1st PDF component)            [(x units)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)            [(x units)^2]
      varnce_y_1, & ! Variance of y (1st PDF component)            [(y units)^2]
      varnce_y_2, & ! Variance of y (2nd PDF component)            [(y units)^2]
      corr_w_x_1, & ! Correlation of w and x (1st PDF component)   [-]
      corr_w_x_2, & ! Correlation of w and x (2nd PDF component)   [-]
      corr_w_y_1, & ! Correlation of w and y (1st PDF component)   [-]
      corr_w_y_2, & ! Correlation of w and y (2nd PDF component)   [-]
      corr_x_y_1, & ! Correlation of x and y (1st PDF component)   [-]
      corr_x_y_2, & ! Correlation of x and y (2nd PDF component)   [-]
      mixt_frac     ! Mixture fraction                             [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      wpxpyp    ! <w'x'y'>                   [m/s (units vary)]

    ! Local Variables
    integer :: i, k


    ! Calculate <w'x'y'> by integrating over the PDF.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        wpxpyp(i,k) &
        = mixt_frac(i,k) &
          * ( ( w_1(i,k) - wm(i,k) ) * ( x_1(i,k) - xm(i,k) ) * ( y_1(i,k) - ym(i,k) ) &
              + corr_x_y_1(i,k)*sqrt( varnce_x_1(i,k)*varnce_y_1(i,k) )*( w_1(i,k)-wm(i,k) ) &
              + corr_w_y_1(i,k)*sqrt( varnce_w_1(i,k)*varnce_y_1(i,k) )*( x_1(i,k)-xm(i,k) ) &
              + corr_w_x_1(i,k)*sqrt( varnce_w_1(i,k)*varnce_x_1(i,k) )*( y_1(i,k)-ym(i,k) ) ) &
          + ( one - mixt_frac(i,k) ) &
            * ( ( w_2(i,k) - wm(i,k) )*( x_2(i,k) - xm(i,k) ) * ( y_2(i,k) - ym(i,k) ) &
                + corr_x_y_2(i,k)*sqrt( varnce_x_2(i,k)*varnce_y_2(i,k) )*( w_2(i,k)-wm(i,k) ) &
                + corr_w_y_2(i,k)*sqrt( varnce_w_2(i,k)*varnce_y_2(i,k) )*( x_2(i,k)-xm(i,k) ) &
                + corr_w_x_2(i,k)*sqrt( varnce_w_2(i,k)*varnce_x_2(i,k) )*( y_2(i,k)-ym(i,k) ) )
      end do
    end do
    !$acc end parallel loop

    return

  end subroutine calc_wpxpyp_pdf

  !=============================================================================
  subroutine calc_liquid_cloud_frac_component( nz, ngrdcol, &
                                               mean_chi, stdev_chi, &
                                               cloud_frac, rc )
    ! Description:
    ! Calculates the PDF component cloud water mixing ratio, rc_i, and cloud
    ! fraction, cloud_frac_i, for the ith PDF component.
    !
    ! The equation for cloud water mixing ratio, rc, at any point is:
    !
    ! rc = chi * H(chi);
    !
    ! and the equation for cloud fraction at a point, fc, is:
    !
    ! fc = H(chi);
    !
    ! where where extended liquid water mixing ratio, chi, is equal to cloud
    ! water mixing ratio, rc, when positive.  When the atmosphere is saturated
    ! at this point, cloud water is found, and rc = chi, while fc = 1.
    ! Otherwise, clear air is found at this point, and rc = fc = 0.
    !
    ! The mean of rc and fc is calculated by integrating over the PDF, such
    ! that:
    !
    ! <rc> = INT(-inf:inf) chi * H(chi) * P(chi) dchi; and
    !
    ! cloud_frac = <fc> = INT(-inf:inf) H(chi) * P(chi) dchi.
    !
    ! This can be rewritten as:
    !
    ! <rc> = INT(0:inf) chi * P(chi) dchi; and
    !
    ! cloud_frac = <fc> = INT(0:inf) P(chi) dchi;
    !
    ! and further rewritten as:
    !
    ! <rc> = SUM(i=1,N) mixt_frac_i INT(0:inf) chi * P_i(chi) dchi; and
    !
    ! cloud_frac = SUM(i=1,N) mixt_frac_i INT(0:inf) P_i(chi) dchi;
    !
    ! where N is the number of PDF components.  The equation for mean rc in the
    ! ith PDF component is:
    !
    ! rc_i = INT(0:inf) chi * P_i(chi) dchi;
    !
    ! and the equation for cloud fraction in the ith PDF component is:
    ! 
    ! cloud_frac_i = INT(0:inf) P_i(chi) dchi.
    !
    ! The component values are related to the overall values by:
    !
    ! <rc> = SUM(i=1,N) mixt_frac_i * rc_i; and
    !
    ! cloud_frac = SUM(i=1,N) mixt_frac_i * cloud_frac_i.

    ! References:
    !----------------------------------------------------------------------

    use constants_clubb, only: &
        chi_tol,        & ! Tolerance for pdf parameter chi       [kg/kg]
        sqrt_2pi,       & ! sqrt(2*pi)
        sqrt_2,         & ! sqrt(2)
        one,            & ! 1
        one_half,       & ! 1/2
        zero,           & ! 0
        max_num_stdevs, &
        eps

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    !----------- Input Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      mean_chi,  & ! Mean of chi (old s) (ith PDF component)           [kg/kg]
      stdev_chi    ! Standard deviation of chi (ith PDF component)     [kg/kg]

    !----------- Output Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      cloud_frac, & ! Cloud fraction (ith PDF component)               [-]
      rc            ! Mean cloud water mixing ratio (ith PDF comp.)    [kg/kg]

    !----------- Local Variables -----------
    real( kind = core_rknd), parameter :: &
      invrs_sqrt_2 = one / sqrt_2, &
      invrs_sqrt_2pi = one / sqrt_2pi

    real( kind = core_rknd ) :: &
      zeta

    integer :: k, i    ! Vertical loop index

    !----------- Begin Code -----------
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        if ( ( abs( mean_chi(i,k) ) <= eps .and. stdev_chi(i,k) <= chi_tol ) &
               .or. ( mean_chi(i,k) < - max_num_stdevs * stdev_chi(i,k) ) ) then

            ! The mean of chi is at saturation and does not vary in the ith PDF component
            cloud_frac(i,k) = zero
            rc(i,k)         = zero

        elseif ( mean_chi(i,k) > max_num_stdevs * stdev_chi(i,k) ) then

            ! The mean of chi is multiple standard deviations above the saturation point.
            ! Thus, all cloud in the ith PDF component.
            cloud_frac(i,k) = one
            rc(i,k)         = mean_chi(i,k)

        else

            ! The mean of chi is within max_num_stdevs of the saturation point.
            ! Thus, layer is partly cloudy, requires calculation.

            zeta = mean_chi(i,k) / stdev_chi(i,k)

            cloud_frac(i,k) = one_half * ( one + erf( zeta * invrs_sqrt_2 )  )

            rc(i,k) = mean_chi(i,k) * cloud_frac(i,k) &
                      + stdev_chi(i,k) * exp( - one_half * zeta**2 ) * invrs_sqrt_2pi

        end if

      end do
    end do
    !$acc end parallel loop

    return

  end subroutine calc_liquid_cloud_frac_component

  !=============================================================================
  subroutine calc_ice_cloud_frac_component( nz, ngrdcol, &
                                            mean_chi, stdev_chi, &
                                            rc_in, cloud_frac, &
                                            p_in_Pa, tl, &
                                            rsatl, crt, &
                                            saturation_formula, &
                                            ice_supersat_frac, rc )
  ! Description:
  !   A version of the cloud fraction calculation modified to work
  !   for layers that are potentially below freezing. If there are
  !   no below freezing levels, the ice_supersat_frac calculation is 
  !   the same as cloud_frac. 
  !
  !   For the below freezing levels, the saturation point will be
  !   non-zero, thus we need to calculate chi_at_ice_sat.
  !
  !   The description of the equations are located in the description
  !   of calc_liquid_cloud_frac_component.
  !----------------------------------------------------------------------

    use constants_clubb, only: &
        chi_tol,        & ! Tolerance for pdf parameter chi       [kg/kg]
        T_freeze_K,     & ! Freezing point of water             [K]
        sqrt_2pi,       & ! sqrt(2*pi)
        sqrt_2,         & ! sqrt(2)
        one,            & ! 1
        one_half,       & ! 1/2
        zero,           & ! 0
        max_num_stdevs, &
        eps

    use clubb_precision, only: &
        core_rknd     ! Precision

    use saturation, only: &
        sat_mixrat_ice

    implicit none

    ! ---------------------- Input Variables ----------------------
    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      mean_chi,   & ! Mean of chi (old s) (ith PDF component)           [kg/kg]
      stdev_chi,  & ! Standard deviation of chi (ith PDF component)     [kg/kg]
      rc_in,      & ! Mean cloud water mixing ratio (ith PDF comp.)     [kg/kg]
      cloud_frac, & ! Cloud fraction                                    [-]
      p_in_Pa,    & ! Pressure                                          [Pa]
      rsatl,      & ! Saturation mixing ratio of liquid                 [kg/kg]
      crt,        & ! r_t coef. in chi/eta eqns.                        [-]
      tl            ! Quantities needed to predict higher order moments
                    ! tl = thl*exner

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    ! ---------------------- Output Variables ----------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      ice_supersat_frac,  & ! Ice supersaturation fraction                [-]
      rc                    ! Mean cloud ice mixing ratio (ith PDF comp.) [kg/kg]

    ! ---------------------- Local Variables----------------------
    real( kind = core_rknd), parameter :: &
      invrs_sqrt_2 = one / sqrt_2, &
      invrs_sqrt_2pi = one / sqrt_2pi

    real( kind = core_rknd ) :: &
      zeta, &
      chi_at_ice_sat

    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      rsat_ice

    integer :: k, i    ! Loop indices

    logical :: &
      l_any_below_freezing

    ! ---------------------- Begin Code ----------------------

    l_any_below_freezing = .false.

    ! If a grid boxes is above freezing, then the calculation is the 
    ! same as the cloud_frac calculation
    !$acc parallel loop gang vector collapse(2) default(present) &
    !$acc          reduction(.or.:l_any_below_freezing)
    do k = 1, nz
      do i = 1, ngrdcol 
        if ( tl(i,k) > T_freeze_K ) then
          ice_supersat_frac(i,k) = cloud_frac(i,k)
          rc(i,k)                = rc_in(i,k)
        else
          l_any_below_freezing = .true.
        end if
      end do
    end do
    !$acc end parallel loop

    ! If all grid boxes are above freezing, then the calculation is the
    ! same as the cloud_frac calculation
    if ( .not. l_any_below_freezing ) then
      return
    end if

    !$acc data create( rsat_ice )

    ! Calculate the saturation mixing ratio of ice
    rsat_ice = sat_mixrat_ice( nz, ngrdcol, p_in_Pa, tl, saturation_formula )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        if ( tl(i,k) <= T_freeze_K ) then

          ! Temperature is freezing, we must compute chi_at_ice_sat and
          ! calculate the new cloud_frac_component
          chi_at_ice_sat = crt(i,k) * ( rsat_ice(i,k) - rsatl(i,k) )

          if ( ( abs( mean_chi(i,k)-chi_at_ice_sat ) <= eps .and. stdev_chi(i,k) <= chi_tol ) &
           .or. ( mean_chi(i,k)-chi_at_ice_sat < - max_num_stdevs * stdev_chi(i,k) ) ) then

            ! The mean of chi is at saturation and does not vary in the ith PDF component
            ice_supersat_frac(i,k) = zero
            rc(i,k)         = zero

          elseif ( mean_chi(i,k)-chi_at_ice_sat > max_num_stdevs * stdev_chi(i,k) ) then

            ! The mean of chi is multiple standard deviations above the saturation point.
            ! Thus, all cloud in the ith PDF component.
            ice_supersat_frac(i,k) = one
            rc(i,k)         = mean_chi(i,k)-chi_at_ice_sat

          else

            ! The mean of chi is within max_num_stdevs of the saturation point.
            ! Thus, layer is partly cloudy, requires calculation.

            zeta = (mean_chi(i,k)-chi_at_ice_sat) / stdev_chi(i,k)

            ice_supersat_frac(i,k) = one_half * ( one + erf( zeta * invrs_sqrt_2 ) )

            rc(i,k) = (mean_chi(i,k)-chi_at_ice_sat) * ice_supersat_frac(i,k) &
                      + stdev_chi(i,k) * exp( - one_half * zeta**2 ) * invrs_sqrt_2pi

          end if

        end if

      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end subroutine calc_ice_cloud_frac_component

  !=============================================================================
  subroutine calc_xprcp_component( nz, ngrdcol,                                     & ! In
                                   wm, rtm, thlm, um, vm, rcm,                      & ! In
                                   w_i, rt_i,                                       & ! In
                                   thl_i, u_i, v_i,                                 & ! In
                                   varnce_w_i, chi_i,                               & ! In
                                   stdev_chi_i, stdev_eta_i,                        & ! In
                                   corr_w_chi_i, corr_chi_eta_i,                    & ! In
!                                  corr_u_w_i, corr_v_w_i,                          & ! In
                                   crt_i, cthl_i,                                   & ! In
                                   rc_i, cloud_frac_i, iiPDF_type,                  & ! In
                                   wprcp_contrib_comp_i, wp2rcp_contrib_comp_i,     & ! Out
                                   rtprcp_contrib_comp_i, thlprcp_contrib_comp_i,   & ! Out
                                   uprcp_contrib_comp_i, vprcp_contrib_comp_i )       ! Out

    ! Description:
    ! Calculates the contribution to <w'rc'>, <w'^2 rc'>, <rt'rc'>, and
    ! <thl'rc'> from the ith PDF component.
    !
    !
    ! <w'rc'>
    ! -------
    !
    ! The value of <w'rc'> is calculated by integrating over the PDF:
    !
    ! <w'rc'>
    ! = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> ) ( rc - <rc> ) P(w,rt,thl) dthl drt dw;
    !
    ! where <w> is the overall mean of w, <rc> is the overall mean of rc, and
    ! P(w,rt,thl) is a two-component trivariate normal distribution of w, rt,
    ! and thl.  This equation is rewritten as:
    !
    ! <w'rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !     ( w - <w> ) ( rc - <rc> ) P_1(w,rt,thl) dthl drt dw
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !       ( w - <w> ) ( rc - <rc> ) P_2(w,rt,thl) dthl drt dw;
    !
    ! where mixt_frac is the mixture fraction, which is the weight of the 1st
    ! PDF component, and where P_1(w,rt,thl) and P_2(w,rt,thl) are the equations
    ! for the trivariate normal PDF of w, rt, and thl in the 1st and 2nd PDF
    ! components, respectively.  The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> ) ( rc - <rc> ) P_i(w,rt,thl) dthl drt dw;
    !
    ! where P_i(w,rt,thl) is the trivariate normal PDF of w, rt, and thl in the
    ! ith PDF component.  The PDF undergoes a PDF transformation in each PDF
    ! component, which is a change of variables and a translation, stretching,
    ! and rotation of the axes.  The PDF becomes a trivariate normal PDF that is
    ! written in terms of w, chi, and eta coordinates.  Cloud water mixing
    ! ratio, rc, is written in terms of extended liquid water mixing ratio, chi,
    ! such that:
    !
    ! rc = chi H(chi);
    !
    ! where H(chi) is the Heaviside step function.  The contribution from the
    ! ith PDF component to <w'rc'> can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> ) ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw;
    !
    ! where P_i(w,chi) is the bivariate normal PDF of w and chi in the ith PDF
    ! component.  The solved equation for the <w'rc'> contribution from the ith
    ! PDF component (wprcp_contrib_comp_i) is:
    !
    ! wprcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> ) ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw
    ! = ( mu_w_i - <w> )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + corr_w_chi_i * sigma_w_i * sigma_chi_i
    !     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
    !
    ! where mu_w_i is the mean of w in the ith PDF component, mu_chi_i is the
    ! mean of chi in the ith PDF component, sigma_w_i is the standard deviation
    ! of w in the ith PDF component, sigma_chi_i is the standard deviation of
    ! chi in the ith PDF component, and corr_w_chi_i is the correlation of w and
    ! chi in the ith PDF component.
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! wprcp_contrib_comp_i
    ! = | ( mu_w_i - <w> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
    !   | ( mu_w_i - <w> ) * ( -<rc> ); when mu_chi_i <= 0.
    !
    !
    ! <w'^2 rc'>
    ! ----------
    !
    ! The value of <w'^2 rc'> is calculated by integrating over the PDF:
    !
    ! <w'^2 rc'>
    ! = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> )^2 ( rc - <rc> ) P(w,rt,thl) dthl drt dw.
    !
    ! This equation is rewritten as:
    !
    ! <w'^2 rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !     ( w - <w> )^2 ( rc - <rc> ) P_1(w,rt,thl) dthl drt dw
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !       ( w - <w> )^2 ( rc - <rc> ) P_2(w,rt,thl) dthl drt dw.
    !
    ! The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> )^2 ( rc - <rc> ) P_i(w,rt,thl) dthl drt dw.
    !
    ! The PDF undergoes a PDF transformation in each PDF component, and becomes
    ! a trivariate normal PDF that is written in terms of w, chi, and eta
    ! coordinates.  The contribution from the ith PDF component to <w'^2 rc'>
    ! can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> )^2 ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw.
    !
    ! The solved equation for the <w'^2 rc'> contribution from the ith PDF
    ! component (wp2rcp_contrib_comp_i) is:
    !
    ! wp2rcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> )^2 ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw
    ! = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + ( mu_w_i - <w> ) * corr_w_chi_i * sigma_w_i * sigma_chi_i
    !     * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !   + 1/sqrt(2*pi) * corr_w_chi_i^2 * sigma_w_i^2 * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }.
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! wp2rcp_contrib_comp_i
    ! = | ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( mu_chi_i - <rc> );
    !   |     when mu_chi_i > 0;
    !   | ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( -<rc> );
    !   |     when mu_chi_i <= 0.
    !
    !
    ! <rt'rc'>
    ! --------
    !
    ! The value of <rt'rc'> is calculated by integrating over the PDF:
    !
    ! <rt'rc'>
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( rt - <rt> ) ( rc - <rc> ) P(rt,thl) dthl drt;
    !
    ! where <rt> is the overall mean of rt, and where P(rt,thl) is a
    ! two-component bivariate normal distribution of rt and thl.  This equation
    ! is rewritten as:
    !
    ! <rt'rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf)
    !     ( rt - <rt> ) ( rc - <rc> ) P_1(rt,thl) dthl drt
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf)
    !       ( rt - <rt> ) ( rc - <rc> ) P_2(rt,thl) dthl drt;
    !
    ! where P_1(rt,thl) and P_2(rt,thl) are the equations for the bivariate
    ! normal PDF of rt and thl in the 1st and 2nd PDF components, respectively.
    ! The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( rt - <rt> ) ( rc - <rc> ) P_i(rt,thl) dthl drt;
    !
    ! where P_i(rt,thl) is the bivariate normal PDF of rt and thl in the ith PDF
    ! component.  The PDF undergoes a PDF transformation in each PDF component,
    ! and becomes a bivariate normal PDF that is written in terms of chi and
    ! eta coordinates.  Total water mixing ratio, rt, is rewritten in terms of
    ! chi and eta by:
    !
    ! rt = mu_rt_i
    !      + ( ( eta - mu_eta_i ) + ( chi - mu_chi_i ) ) / ( 2 * crt_i );
    !
    ! where mu_rt_i is the mean of rt in the ith PDF component, mu_eta_i is the
    ! mean of eta in the ith PDF component, and crt_i is a coefficient on rt in
    ! the chi/eta transformation equations.  The contribution from the ith PDF
    ! component to <rt'rc'> can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( mu_rt_i - <rt> + ( eta - mu_eta_i ) / ( 2 * crt_i )
    !   + ( chi - mu_chi_i ) / ( 2 * crt_i ) )
    ! * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi;
    !
    ! where P_i(chi,eta) is the bivariate normal PDF of chi and eta in the ith
    ! PDF component.  The solved equation for the <rt'rc'> contribution from the
    ! ith PDF component (rtprcp_contrib_comp_i) is:
    !
    ! rtprcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( mu_rt_i - <rt> + ( eta - mu_eta_i ) / ( 2 * crt_i )
    !     + ( chi - mu_chi_i ) / ( 2 * crt_i ) )
    !   * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi
    ! = ( mu_rt_i - <rt> )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i + sigma_chi_i ) / ( 2 * crt_i )
    !     * sigma_chi_i
    !     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
    !
    ! where sigma_eta_i is the standard deviation of eta in the ith PDF
    ! component and corr_chi_eta_i is the correlation of chi and eta in the ith
    ! PDF component.
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! rtprcp_contrib_comp_i
    ! = | ( mu_rt_i - <rt> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
    !   | ( mu_rt_i - <rt> ) * ( -<rc> ); when mu_chi_i <= 0.
    !
    !
    ! <thl'rc'>
    ! ---------
    !
    ! The value of <thl'rc'> is calculated by integrating over the PDF:
    !
    ! <thl'rc'>
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( thl - <thl> ) ( rc - <rc> ) P(rt,thl) dthl drt;
    !
    ! where <thl> is the overall mean of thl.  This equation is rewritten as:
    !
    ! <thl'rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf)
    !     ( thl - <thl> ) ( rc - <rc> ) P_1(rt,thl) dthl drt
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf)
    !       ( thl - <thl> ) ( rc - <rc> ) P_2(rt,thl) dthl drt.
    !
    ! The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( thl - <thl> ) ( rc - <rc> ) P_i(rt,thl) dthl drt.
    !
    ! The PDF undergoes a PDF transformation in each PDF component, and becomes
    ! a bivariate normal PDF that is written in terms of chi and eta
    ! coordinates.  Liquid water potential temperature, thl, is rewritten in
    ! terms of chi and eta by:
    !
    ! thl = mu_thl_i
    !       + ( ( eta - mu_eta_i ) - ( chi - mu_chi_i ) ) / ( 2 * cthl_i );
    !
    ! where mu_thl_i is the mean of thl in the ith PDF component and cthl_i is a
    ! coefficient on thl in the chi/eta transformation equations.  The
    ! contribution from the ith PDF component to <thl'rc'> can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( mu_thl_i - <thl> + ( eta - mu_eta_i ) / ( 2 * cthl_i )
    !   - ( chi - mu_chi_i ) / ( 2 * cthl_i ) )
    ! * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi.
    !
    ! The solved equation for the <thl'rc'> contribution from the ith PDF
    ! component (thlprcp_contrib_comp_i) is:
    !
    ! thlprcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( mu_thl_i - <thl> + ( eta - mu_eta_i ) / ( 2 * cthl_i )
    !     - ( chi - mu_chi_i ) / ( 2 * cthl_i ) )
    !   * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi
    ! = ( mu_thl_i - <thl> )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i - sigma_chi_i ) / ( 2 * cthl_i )
    !     * sigma_chi_i
    !     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! thlprcp_contrib_comp_i
    ! = | ( mu_thl_i - <thl> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
    !   | ( mu_thl_i - <thl> ) * ( -<rc> ); when mu_chi_i <= 0.
    !
    !
    ! Use equations for PDF component cloud fraction cloud water mixing ratio
    ! -----------------------------------------------------------------------
    !
    ! The equation for cloud fraction in the ith PDF component, fc_i, is:
    !
    ! fc_i = 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
    !
    ! In the special case that sigma_chi_i = 0, the equation becomes:
    !
    ! fc_i = | 1; when mu_chi_i > 0;
    !        | 0; when mu_chi_i <= 0.
    !
    ! The equation for mean cloud water mixing ratio in the ith PDF component,
    ! rc_i, is:
    !
    ! rc_i
    ! = mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !   + 1/sqrt(2*pi) * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }
    ! = mu_chi_i * fc_i
    !   + 1/sqrt(2*pi) * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }.
    !
    ! In the special case that sigma_chi_i = 0, the equation becomes:
    !
    ! rc_i = | mu_chi_i; when mu_chi_i > 0;
    !        | 0; when mu_chi_i <= 0.
    !
    ! The above equations can be substituted into the equations for
    ! wprcp_contrib_comp_i, wp2rcp_contrib_comp_i, rtprcp_contrib_comp_i, and
    ! thlprcp_contrib_comp_i.  The new equations are:
    !
    ! wprcp_contrib_comp_i
    ! = ( mu_w_i - <w> ) * ( rc_i - <rc> )
    !   + corr_w_chi_i * sigma_w_i * sigma_chi_i * fc_i;
    !
    ! wp2rcp_contrib_comp_i
    ! = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( rc_i - <rc> )
    !   + 2 * ( mu_w_i - <w> ) * corr_w_chi_i * sigma_w_i * sigma_chi_i * fc_i
    !   + 1/sqrt(2*pi) * corr_w_chi_i^2 * sigma_w_i^2 * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) };
    !
    ! rtprcp_contrib_comp_i
    ! = ( mu_rt_i - <rt> ) * ( rc_i - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i + sigma_chi_i ) / ( 2 * crt_i )
    !     * sigma_chi_i * fc_i; and
    !
    ! thlprcp_contrib_comp_i
    ! = ( mu_thl_i - <thl> ) * ( rc_i - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i - sigma_chi_i ) / ( 2 * cthl_i )
    !     * sigma_chi_i * fc_i.
    !
    ! While the above equations reduce to their listed versions in the special
    ! case that sigma_chi_i = 0, those versions are faster to calculate.  When
    ! mu_chi_i > 0, they are:
    !
    ! wprcp_contrib_comp_i = ( mu_w_i - <w> ) * ( mu_chi_i - <rc> );
    ! wp2rcp_contrib_comp_i
    ! = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( mu_chi_i - <rc> );
    ! rtprcp_contrib_comp_i = ( mu_rt_i - <rt> ) * ( mu_chi_i - <rc> ); and
    ! thlprcp_contrib_comp_i = ( mu_thl_i - <thl> ) * ( mu_chi_i - <rc> );
    !
    ! and when mu_chi_i <= 0, they are:
    !
    ! wprcp_contrib_comp_i = - ( mu_w_i - <w> ) * <rc>;
    ! wp2rcp_contrib_comp_i = - ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * <rc>;
    ! rtprcp_contrib_comp_i = - ( mu_rt_i - <rt> ) * <rc>; and
    ! thlprcp_contrib_comp_i = - ( mu_thl_i - <thl> ) * <rc>.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        sqrt_2pi,       & ! Variable(s)
        two,            &
        zero,           &
        chi_tol

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      wm,             & ! Mean of w (overall)                          [m/s]
      rtm,            & ! Mean of rt (overall)                         [kg/kg]
      thlm,           & ! Mean of thl (overall)                        [K]
      um,             & ! Mean of eastward wind (overall)              [m/s]
      vm,             & ! Mean of northward wind (overall)             [m/s]
      rcm,            & ! Mean of rc (overall)                         [kg/kg]
      w_i,            & ! Mean of w (ith PDF component)                [m/s]
      rt_i,           & ! Mean of rt (ith PDF component)               [kg/kg]
      thl_i,          & ! Mean of thl (ith PDF component)              [K]
      u_i,            & ! Mean of eastward wind (ith PDF component)    [m/s]
      v_i,            & ! Mean of northward wind (ith PDF component)   [m/s]
      varnce_w_i,     & ! Variance of w (ith PDF component)            [m^2/s^2]
      chi_i,          & ! Mean of chi (ith PDF component)              [kg/kg]
      stdev_chi_i,    & ! Standard deviation of chi (ith PDF comp.)    [kg/kg]
      stdev_eta_i,    & ! Standard deviation of eta (ith PDF comp.)    [kg/kg]
      corr_w_chi_i,   & ! Correlation of w and chi (ith PDF component) [-]
      corr_chi_eta_i, & ! Correlation of chi and eta (ith PDF comp.)   [-]
!     corr_u_w_i,     & ! Correlation of u and w (ith PDF component)   [-]
!     corr_v_w_i,     & ! Correlation of v and w (ith PDF component)   [-]
      crt_i,          & ! Coef. on rt in chi/eta eqns. (ith PDF comp.) [-]
      cthl_i,         & ! Coef. on thl: chi/eta eqns. (ith PDF comp.)  [kg/kg/K]
      rc_i,           & ! Mean of rc (ith PDF component)               [kg/kg]
      cloud_frac_i      ! Cloud fraction (ith PDF component)           [-]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      wprcp_contrib_comp_i,   & ! <w'rc'> contrib. (ith PDF comp.)  [m/s(kg/kg)]
      wp2rcp_contrib_comp_i,  & ! <w'^2rc'> contrib. (ith comp) [m^2/s^2(kg/kg)]
      rtprcp_contrib_comp_i,  & ! <rt'rc'> contrib. (ith PDF comp.)  [kg^2/kg^2]
      thlprcp_contrib_comp_i, & ! <thl'rc'> contrib. (ith PDF comp.)  [K(kg/kg)]
      uprcp_contrib_comp_i,   & ! <u'rc'> contrib. (ith PDF comp.)  [m/s(kg/kg)]
      vprcp_contrib_comp_i      ! <v'rc'> contrib. (ith PDF comp.)  [m/s(kg/kg)]

    ! Local Variables
    integer :: i, k

    ! ---------------------- Begin Code ------------------
    
    ! Changing these conditionals may result in inconsistencies with the conditional
    ! statements located in calc_cloud_frac_component
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        wprcp_contrib_comp_i(i,k) = ( w_i(i,k) - wm(i,k) ) * ( rc_i(i,k) - rcm(i,k) )

        wp2rcp_contrib_comp_i(i,k) = ( ( w_i(i,k) - wm(i,k) )**2 + varnce_w_i(i,k) ) &
                                     * ( rc_i(i,k) - rcm(i,k) )

        rtprcp_contrib_comp_i(i,k) = ( rt_i(i,k) - rtm(i,k) ) * ( rc_i(i,k) - rcm(i,k) ) &
                                + ( corr_chi_eta_i(i,k) * stdev_eta_i(i,k) + stdev_chi_i(i,k) ) &
                                  / ( two * crt_i(i,k) ) * stdev_chi_i(i,k) * cloud_frac_i(i,k)

        thlprcp_contrib_comp_i(i,k) = ( thl_i(i,k) - thlm(i,k) ) * ( rc_i(i,k) - rcm(i,k) ) &
                                 + ( corr_chi_eta_i(i,k) * stdev_eta_i(i,k) - stdev_chi_i(i,k) ) &
                                   / ( two * cthl_i(i,k) ) * stdev_chi_i(i,k) * cloud_frac_i(i,k)

        uprcp_contrib_comp_i(i,k) = ( u_i(i,k) - um(i,k) ) * ( rc_i(i,k) - rcm(i,k) )

        vprcp_contrib_comp_i(i,k) = ( v_i(i,k) - vm(i,k) ) * ( rc_i(i,k) - rcm(i,k) )

      end do
    end do
    !$acc end parallel loop

    ! If iiPDF_type isn't iiPDF_ADG1, iiPDF_ADG2, or iiPDF_new_hybrid, so
    ! corr_w_chi_i /= 0 (and perhaps corr_u_w_i /= 0).
    if ( .not. ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 &
                 .or. iiPDF_type == iiPDF_new_hybrid ) ) then

        ! Chi varies significantly in the ith PDF component (stdev_chi > chi_tol)
        ! and there is some cloud (0 < cloud_frac <= 1)
        do k = 1, nz
          do i = 1, ngrdcol
            if ( stdev_chi_i(i,k) > chi_tol .and. cloud_frac_i(i,k) > zero ) then

              wprcp_contrib_comp_i(i,k) = wprcp_contrib_comp_i(i,k) &
                                          + corr_w_chi_i(i,k) * sqrt( varnce_w_i(i,k) ) &
                                            * stdev_chi_i(i,k) * cloud_frac_i(i,k)

              wp2rcp_contrib_comp_i(i,k) = wp2rcp_contrib_comp_i(i,k) &
                                           + two * ( w_i(i,k) - wm(i,k) ) * corr_w_chi_i(i,k) &
                                             * sqrt( varnce_w_i(i,k) ) * stdev_chi_i(i,k) &
                                             * cloud_frac_i(i,k) &
                                           + corr_w_chi_i(i,k)**2 * varnce_w_i(i,k) &
                                             * stdev_chi_i(i,k) &
                                             * exp( -chi_i(i,k)**2 / ( two*stdev_chi_i(i,k)**2 ) ) &
                                               / sqrt_2pi

            ! In principle, uprcp_contrib_comp_i might depend on corr_u_w_i here.
          end if
        end do
      end do
    end if 

    return

  end subroutine calc_xprcp_component

  !=============================================================================
  subroutine calc_w_up_in_cloud( nz, ngrdcol, &
                                 mixt_frac, cloud_frac_1, cloud_frac_2, &
                                 w_1, w_2, varnce_w_1, varnce_w_2, &
                                 w_up_in_cloud, w_down_in_cloud, &
                                 cloudy_updraft_frac, cloudy_downdraft_frac )

    ! Description:
    ! Subroutine that computes the mean cloudy updraft (and also calculates
    ! the mean cloudy downdraft).
    !
    ! In order to activate aerosol, we'd like to feed the activation scheme
    ! a vertical velocity that's representative of cloudy updrafts. For skewed
    ! layers, like cumulus layers, this might be an improvement over the square
    ! root of wp2 that's currently used. At the same time, it would be simpler
    ! and less expensive than feeding SILHS samples into the aerosol code
    ! (see larson-group/e3sm#19 and larson-group/e3sm#26).
    !
    ! The formulas are only valid for certain PDFs in CLUBB (ADG1, ADG2,
    ! new hybrid), hence we omit calculation if another PDF type is used.
    !
    ! References: https://www.overleaf.com/project/614a136d47846639af22ae34
    !----------------------------------------------------------------------

    use constants_clubb, only: &
        sqrt_2pi,       & ! sqrt(2*pi)
        sqrt_2,         & ! sqrt(2)
        one,            & ! 1
        one_half,       & ! 1/2
        zero,           & ! 0
        max_num_stdevs, &
        eps

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nz          ! Number of vertical level

    !----------- Input Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      mixt_frac, &      ! mixture fraction                             [-]
      cloud_frac_1, &   ! cloud fraction (1st PDF component)           [-]
      cloud_frac_2, &   ! cloud fraction (2nd PDF component)           [-]
      w_1, &            ! upward velocity (1st PDF component)          [m/s]
      w_2, &            ! upward velocity (2nd PDF component)          [m/s]
      varnce_w_1, &     ! standard deviation of w (1st PDF component)  [m^2/s^2]
      varnce_w_2        ! standard deviation of w (2nd PDF component)  [m^2/s^2]

    !----------- Output Variables -----------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      w_up_in_cloud,         & ! mean cloudy updraft speed             [m/s]
      w_down_in_cloud,       & ! mean cloudy downdraft speed           [m/s]
      cloudy_updraft_frac,   & ! cloudy updraft fraction               [-]
      cloudy_downdraft_frac    ! cloudy downdraft fraction             [-]

    !----------- Local Variables -----------
    real( kind = core_rknd ) :: &
      w_up_1, w_up_2, &        ! integral of w and Heaviside fnc, where w > 0
      w_down_1, w_down_2, &    ! integral of w and Heaviside fnc, where w < 0
      stdev_w_1, stdev_w_2, &  ! Standard deviation of w
      ratio_w_1, &             ! mu_w_1 / ( sqrt(2) * sigma_w_1 )
      ratio_w_2, &             ! mu_w_2 / ( sqrt(2) * sigma_w_2 )
      erf_ratio_w_1, &         ! erf( ratio_w_1 )
      erf_ratio_w_2, &         ! erf( ratio_w_2 )
      exp_neg_ratio_w_1_sqd, & ! exp( -ratio_w_1^2 )
      exp_neg_ratio_w_2_sqd, & ! exp( -ratio_w_2^2 )
      updraft_frac_1, &        ! Fraction of 1st PDF comp. where w > 0
      updraft_frac_2, &        ! Fraction of 2nd PDF comp. where w > 0
      downdraft_frac_1, &      ! Fraction of 1st PDF comp. where w < 0
      downdraft_frac_2         ! Fraction of 2nd PDF comp. where w < 0

    integer :: i, k

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol

        stdev_w_1 = sqrt(varnce_w_1(i,k))
        stdev_w_2 = sqrt(varnce_w_2(i,k))

        ! Calculate quantities in the 1st PDF component.
        if ( w_1(i,k) > max_num_stdevs * stdev_w_1 ) then

           ! The mean of w in the 1st PDF component is more than
           ! max_num_stdevs standard deviations above 0.
           ! The entire 1st PDF component is found in an updraft (w > 0).
           w_up_1 = w_1(i,k)
           updraft_frac_1 = one
           w_down_1 = zero
           downdraft_frac_1 = zero

        elseif ( w_1(i,k) < - max_num_stdevs * stdev_w_1 ) then

           ! The mean of w in the 1st PDF component is more than
           ! max_num_stdevs standard deviations below 0.
           ! The entire 1st PDF component is found in a downdraft (w < 0).
           w_up_1 = zero
           updraft_frac_1 = zero
           w_down_1 = w_1(i,k)
           downdraft_frac_1 = one

        else

           ! The 1st PDF component contains both updraft and downdraft.
           ratio_w_1 = w_1(i,k) / ( sqrt_2 * max(eps, stdev_w_1) )
           erf_ratio_w_1 = erf( ratio_w_1 )
           exp_neg_ratio_w_1_sqd = exp( -ratio_w_1**2 )

           w_up_1 &
           = one_half * w_1(i,k) * ( one + erf_ratio_w_1 ) &
             + ( stdev_w_1 / sqrt_2pi ) * exp_neg_ratio_w_1_sqd

           updraft_frac_1 = one_half * ( one + erf_ratio_w_1 )

           w_down_1 &
           = one_half * w_1(i,k) * ( one - erf_ratio_w_1 ) &
             - ( stdev_w_1 / sqrt_2pi ) * exp_neg_ratio_w_1_sqd

           !downdraft_frac_1 = one_half * ( one - erf_ratio_w_1 )
           downdraft_frac_1 = one - updraft_frac_1

        endif

        ! Calculate quantities in the 2nd PDF component.
        if ( w_2(i,k) > max_num_stdevs * stdev_w_2 ) then

           ! The mean of w in the 2nd PDF component is more than
           ! max_num_stdevs standard deviations above 0.
           ! The entire 2nd PDF component is found in an updraft (w > 0).
           w_up_2 = w_2(i,k)
           updraft_frac_2 = one
           w_down_2 = zero
           downdraft_frac_2 = zero

        elseif ( w_2(i,k) < - max_num_stdevs * stdev_w_2 ) then

           ! The mean of w in the 2nd PDF component is more than
           ! max_num_stdevs standard deviations below 0.
           ! The entire 2nd PDF component is found in a downdraft (w < 0).
           w_up_2 = zero
           updraft_frac_2 = zero
           w_down_2 = w_2(i,k)
           downdraft_frac_2 = one

        else

           ! The 2nd PDF component contains both updraft and downdraft.
           ratio_w_2 = w_2(i,k) / ( sqrt_2 * max(eps, stdev_w_2) )
           erf_ratio_w_2 = erf( ratio_w_2 )
           exp_neg_ratio_w_2_sqd = exp( -ratio_w_2**2 )

           w_up_2 &
           = one_half * w_2(i,k) * ( one + erf_ratio_w_2 ) &
             + ( stdev_w_2 / sqrt_2pi ) * exp_neg_ratio_w_2_sqd

           updraft_frac_2 = one_half * ( one + erf_ratio_w_2 )

           w_down_2 &
           = one_half * w_2(i,k) * ( one - erf_ratio_w_2 ) &
             - ( stdev_w_2 / sqrt_2pi ) * exp_neg_ratio_w_2_sqd

           !downdraft_frac_2 = one_half * ( one - erf_ratio_w_2 )
           downdraft_frac_2 = one - updraft_frac_2

        endif

        ! Calculate the total cloudy updraft fraction.
        cloudy_updraft_frac(i,k) &
        = mixt_frac(i,k) * cloud_frac_1(i,k) * updraft_frac_1 &
          + ( one - mixt_frac(i,k) ) * cloud_frac_2(i,k) * updraft_frac_2

        ! Calculate the total cloudy downdraft fraction.
        cloudy_downdraft_frac(i,k) &
        = mixt_frac(i,k) * cloud_frac_1(i,k) * downdraft_frac_1 &
          + ( one - mixt_frac(i,k) ) * cloud_frac_2(i,k) * downdraft_frac_2

        ! Calculate the mean vertical velocity found in a cloudy updraft.
        w_up_in_cloud(i,k) &
        = ( mixt_frac(i,k) * cloud_frac_1(i,k) * w_up_1 &
            + ( one - mixt_frac(i,k) ) * cloud_frac_2(i,k) * w_up_2 ) &
          / max( eps, cloudy_updraft_frac(i,k) )

        ! Calculate the mean vertical velocity found in a cloudy downdraft.
        w_down_in_cloud(i,k) &
        = ( mixt_frac(i,k) * cloud_frac_1(i,k) * w_down_1 &
            + ( one - mixt_frac(i,k) ) * cloud_frac_2(i,k) * w_down_2 ) &
          / max( eps, cloudy_downdraft_frac(i,k) )

      end do
    end do
    !$acc end parallel loop

    return

  end subroutine calc_w_up_in_cloud

  !=============================================================================
  function interp_var_array( n_points, nz, k, z_vals, var )

  ! Description:
  !   Interpolates a variable to an array of values about a given level

  ! References
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd           ! Constant

    implicit none

  ! Input Variables
    integer, intent(in) :: &
      n_points, & ! Number of points to interpolate to (must be odd and >= 3)
      nz,       & ! Total number of vertical levels
      k           ! Center of interpolation array

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      z_vals, &         ! Height at each vertical level           [m]
      var               ! Variable values on grid                 [units vary]

  ! Output Variables
    real( kind = core_rknd ), dimension(n_points) :: &
      interp_var_array  ! Interpolated values of variable         [units vary]

  ! Local Variables
    real( kind = core_rknd ) :: &
      dz    ! Distance between vertical levels

    real( kind = core_rknd ) :: &
      z_val             ! Height at some sub-grid level

    integer :: &
      i, &                      ! Loop iterator

      subgrid_lev_count         ! Number of refined grid points located between
                              ! two defined grid levels

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Place a point at each of k-1, k, and k+1.
    interp_var_array(1) = var_value_integer_height( nz, k-1, z_vals, var )
    interp_var_array((n_points+1)/2) = var_value_integer_height( nz, k, z_vals, var )
    interp_var_array(n_points) = var_value_integer_height( nz, k+1, z_vals, var )

    subgrid_lev_count = (n_points - 3) / 2

    ! Lower half
    if ( k == 1 ) then
      dz = (z_vals(2) - z_vals(1)) / real( subgrid_lev_count+1, kind=core_rknd )
    else
      dz = (z_vals(k) - z_vals(k-1)) / real( subgrid_lev_count+1, kind=core_rknd )
    end if
    do i=1, subgrid_lev_count
      z_val = z_vals(k) - real( i, kind=core_rknd ) * dz
      interp_var_array(1+i) &
      = var_subgrid_interp( nz, k, z_vals, var, z_val, l_below=.true. )
    end do

    ! Upper half
    if ( k == nz ) then
      dz = ( z_vals(nz) - z_vals(nz-1) ) / real( subgrid_lev_count+1, kind=core_rknd )
    else
      dz = ( z_vals(k+1) - z_vals(k) ) / real( subgrid_lev_count+1, kind=core_rknd )
    end if
    do i=1, (n_points-3)/2
      z_val = z_vals(k) + real( i, kind=core_rknd ) * dz
      interp_var_array((n_points+1)/2+i) &
      = var_subgrid_interp( nz, k, z_vals, var, z_val, l_below=.false. )
    end do

    return
  end function interp_var_array

  !=============================================================================
  function var_value_integer_height( nz, k, z_vals, var_grid_value ) result( var_value )

  ! Description
  !   Returns the value of a variable at an integer height between 0 and
  !   nz+1 inclusive, using extrapolation when k==0 or k==nz+1

  ! References
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd       ! Constant

    use interpolation, only: &
        mono_cubic_interp  ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,    & ! Total number of vertical levels
      k        ! Level to resolve variable value

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      z_vals,            & ! Height at each vertical level                  [m]
      var_grid_value       ! Value of variable at each grid level           [units vary]

    ! Output Variables
    real( kind = core_rknd ) :: &
      var_value            ! Value of variable at height level              [units vary]

    ! Local Variables
    integer :: km1, k00, kp1, kp2
  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( k >= 1 .and. k <= nz ) then
      ! This is the simple case. No extrapolation necessary.
      var_value = var_grid_value(k)
    else if ( k == 0 ) then
      ! Extrapolate below the lower boundary
      km1 = nz
      k00 = 1
      kp1 = 2
      kp2 = 3
      var_value = mono_cubic_interp( z_vals(1)-(z_vals(2)-z_vals(1)), &
                                     km1, k00, kp1, kp2, &
                                     z_vals(km1), z_vals(k00), z_vals(kp1), z_vals(kp2), &
                                     var_grid_value(km1), var_grid_value(k00), &
                                     var_grid_value(kp1), var_grid_value(kp2) )
    else if ( k == nz+1 ) then
      ! Extrapolate above the upper boundary
      km1 = nz
      k00 = nz-1
      kp1 = nz
      kp2 = nz
      var_value = mono_cubic_interp( z_vals(nz)+(z_vals(nz)-z_vals(nz-1)), &
                                     km1, k00, kp1, kp2, &
                                     z_vals(km1), z_vals(k00), z_vals(kp1), z_vals(kp2), &
                                     var_grid_value(km1), var_grid_value(k00), &
                                     var_grid_value(kp1), var_grid_value(kp2) )
    else
      ! Invalid height requested
      var_value = -999._core_rknd
    end if ! k > 1 .and. k < nz
    return
  end function var_value_integer_height

  !=============================================================================
  function var_subgrid_interp( nz, k, z_vals, var, z_interp, l_below ) result( var_value )

  ! Description
  !   Interpolates (or extrapolates) a variable to a value between grid
  !   levels

  ! References
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd       ! Constant

    use interpolation, only: &
        mono_cubic_interp   ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &         ! Number of vertical levels
      k             ! Grid level near interpolation target

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      z_vals, &     ! Height at each grid level          [m]
      var           ! Variable values at grid levels     [units vary]

    real( kind = core_rknd ), intent(in) :: &
      z_interp      ! Interpolation target height        [m]

    logical, intent(in) :: &
      l_below       ! True if z_interp < z_vals(k), false otherwise

    ! Output Variable
    real( kind = core_rknd ) :: &
      var_value     ! Interpolated value of variable     [units vary]

    ! Local Variables
    integer :: km1, k00, kp1, kp2 ! Parameters for call to mono_cubic_interp
  !----------------------------------------------------------------------

    !----- Begin Code -----
    if ( l_below ) then

      if ( k == 1 ) then ! Extrapolation
        km1 = nz
        k00 = 1
        kp1 = 2
        kp2 = 3
      else if ( k == 2 ) then
        km1 = 1
        k00 = 1
        kp1 = 2
        kp2 = 3
      else if ( k == nz ) then
        km1 = nz-2
        k00 = nz-1
        kp1 = nz
        kp2 = nz
      else
        km1 = k-2
        k00 = k-1
        kp1 = k
        kp2 = k+1
      end if ! k == 1

    else ! .not. l_below

      if ( k == 1 ) then
        km1 = 1
        k00 = 1
        kp1 = 2
        kp2 = 3
      else if ( k == nz-1 ) then
        km1 = nz-2
        k00 = nz-1
        kp1 = nz
        kp2 = nz
      else if ( k == nz ) then ! Extrapolation
        km1 = nz
        k00 = nz-1
        kp1 = nz
        kp2 = nz
      else
        km1 = k-1
        k00 = k
        kp1 = k+1
        kp2 = k+2
      end if ! k == 1

    end if ! l_below

    ! Now perform the interpolation
    var_value = mono_cubic_interp( z_interp, km1, k00, kp1, kp2, &
                                   z_vals(km1), z_vals(k00), z_vals(kp1), z_vals(kp2), &
                                   var(km1), var(k00), var(kp1), var(kp2) )

    return

  end function var_subgrid_interp

  !=============================================================================

  subroutine pdf_closure_driver( gr, nzm, nzt, ngrdcol,                   & ! Intent(in)
                                 dt, hydromet_dim, sclr_dim, sclr_tol,    & ! Intent(in)
                                 wprtp, thlm, wpthlp, rtp2, rtp3,         & ! Intent(in)
                                 thlp2, thlp3, rtpthlp, wp2,              & ! Intent(in)
                                 wp3, wm_zm, wm_zt,                       & ! Intent(in)
                                 um, up2, upwp, up3,                      & ! Intent(in)
                                 vm, vp2, vpwp, vp3,                      & ! Intent(in)
                                 p_in_Pa, exner,                          & ! Intent(in)
                                 thv_ds_zm, thv_ds_zt, rtm_ref,           & ! Intent(in)
                                 wphydrometp,                             & ! Intent(in)
                                 wp2hmp, rtphmp_zt, thlphmp_zt,           & ! Intent(in)
                                 sclrm, wpsclrp, sclrp2,                  & ! Intent(in)
                                 sclrprtp, sclrpthlp, sclrp3,             & ! Intent(in)
                                 p_sfc, l_samp_stats_in_pdf_call,         & ! Intent(in)
                                 mixt_frac_max_mag, ts_nudge,             & ! Intent(in)
                                 rtm_min, rtm_nudge_max_altitude,         & ! Intent(in)
                                 clubb_params,                            & ! Intent(in)
                                 iiPDF_type,                              & ! Intent(in)
                                 saturation_formula,                      & ! Intent(in)
                                 l_predict_upwp_vpwp,                     & ! Intent(in)
                                 l_rtm_nudge,                             & ! Intent(in)
                                 l_trapezoidal_rule_zt,                   & ! Intent(in)
                                 l_trapezoidal_rule_zm,                   & ! Intent(in)
                                 l_call_pdf_closure_twice,                & ! Intent(in)
                                 l_use_cloud_cover,                       & ! Intent(in)
                                 l_rcm_supersat_adj,                      & ! Intent(in
                                 l_mix_rat_hm,                            & ! Intent(in)
                                 stats,                                   & ! Intent(inout)
                                 rtm,                                     & ! Intent(inout)
                                 pdf_implicit_coefs_terms,                & ! Intent(inout)
                                 pdf_params, pdf_params_zm, err_info,     & ! Intent(inout)
#ifdef GFDL
                                 RH_crit(k, : , :),                       & ! Intent(inout)
                                 do_liquid_only_in_clubb,                 & ! Intent(in)
#endif
                                 rcm, cloud_frac,                         & ! Intent(out)
                                 ice_supersat_frac, wprcp,                & ! Intent(out)
                                 sigma_sqd_w, wpthvp, wp2thvp, wp2up,     & ! Intent(out)
                                 rtpthvp, thlpthvp, rc_coef,              & ! Intent(out)
                                 rcm_in_layer, cloud_cover,               & ! Intent(out)
                                 rcp2_zt, thlprcp,                        & ! Intent(out)
                                 rc_coef_zm, sclrpthvp,                   & ! Intent(out)
                                 wpup2, wpvp2,                            & ! Intent(out)
                                 wp2up2, wp2vp2, wp4,                     & ! Intent(out)
                                 wp2rtp, wprtp2, wp2thlp,                 & ! Intent(out)
                                 wpthlp2, wprtpthlp, wp2rcp,              & ! Intent(out)
                                 rtprcp, rcp2,                            & ! Intent(out)
                                 uprcp, vprcp,                            & ! Intent(out)
                                 w_up_in_cloud, w_down_in_cloud,          & ! Intent(out)
                                 cloudy_updraft_frac,                     & ! Intent(out)
                                 cloudy_downdraft_frac,                   & ! intent(out)
                                 Skw_velocity,                            & ! Intent(out)
                                 cloud_frac_zm,                           & ! Intent(out)
                                 ice_supersat_frac_zm,                    & ! Intent(out)
                                 rtm_zm, thlm_zm, rcm_zm,                 & ! Intent(out)
                                 rcm_supersat_adj,                        & ! Intent(out)
                                 wp2sclrp, wpsclrp2, sclrprcp,            & ! Intent(out)
                                 wpsclrprtp, wpsclrpthlp )                  ! Intent(out)

    use grid_class, only: &
        grid,       & ! Type
        zt2zm_api,  & ! Procedure(s)
        zm2zt_api,  &
        zm2zt2zm

    use constants_clubb, only: &
        one_half,       & ! Variable(s)
        w_tol,          & 
        w_tol_sqd,      &
        rt_tol,         &
        thl_tol,        &
        p0,             &
        kappa,          &
        fstderr,        &
        zero,           &
        zero_threshold, &
        eps

    use pdf_parameter_module, only: &
        pdf_parameter,        & ! Variable Type
        implicit_coefs_terms, &  ! Variable Type
        init_pdf_implicit_coefs_terms_api ! Procedure

    use parameter_indices, only: &
        nparams,         & ! Variable(s)
        igamma_coef,     &
        igamma_coefb,    &
        igamma_coefc

    use Skx_module, only: &
        Skx_func    ! Procedure(s)

    use sigma_sqd_w_module, only: &
        compute_sigma_sqd_w    ! Procedure(s)

    use pdf_utilities, only: &
        compute_mean_binormal    ! Procedure(s)

    use T_in_K_module, only: &
        thlm2T_in_K_api    ! Procedure(s)

    use saturation, only:  &
        sat_mixrat_liq_api    ! Procedure(s)

    use model_flags, only: &
        l_gamma_Skw,      & ! Variable(s)
        iiPDF_new,        & ! new PDF
        iiPDF_new_hybrid    ! new hybrid PDF

    use error_code, only: &
        clubb_at_least_debug_level_api,  & ! Procedure
        clubb_fatal_error              ! Constant

    use stats_netcdf, only: &
        stats_type, &
        stats_update, &
        var_on_stats_list

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use err_info_type_module, only: &
      err_info_type        ! Type

    use clip_explicit, only: &
      clip_rcm ! Procedure(s)

    implicit none

    !------------------------------- Input Variables -------------------------------
    type (grid), intent(in) :: &
      gr

    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol

    real( kind = core_rknd ), intent(in) ::  &
      dt  ! Current timestep duration    [s]

    integer, intent(in) :: &
      hydromet_dim,   & ! Total number of hydrometeor species       [#]
      sclr_dim          ! Number of passive scalars                 [#]

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) :: & 
      sclr_tol          ! Threshold(s) on the passive scalars  [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
      !rtm,       & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      thlm,      & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      rtp3,      & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp3,     & ! th_l'^3 (thermodynamic levels)                 [K^3]
      wp3,       & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
      wm_zt,     & ! w mean wind component on thermo. levels        [m/s]
      p_in_Pa,   & ! Air pressure (thermodynamic levels)            [Pa]
      exner,     & ! Exner function (thermodynamic levels)          [-]
      thv_ds_zt, & ! Dry, base-state theta_v on thermo. levs.       [K]
      rtm_ref      ! Initial total water mixing ratio               [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) ::  &
      wprtp,     & ! w' r_t' (momentum levels)                      [(kg/kg)m/s]
      wpthlp,    & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2,      & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,     & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp,   & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2,       & ! w'^2 (momentum levels)                         [m^2/s^2]
      wm_zm,     & ! w mean wind component on momentum levels       [m/s]
      thv_ds_zm    ! Dry, base-state theta_v on momentum levs.      [K]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
      um,          & ! Grid-mean eastward wind     [m/s]
      up3,         & ! u'^3                        [(m/s)^3]
      vm,          & ! Grid-mean northward wind    [m/s]
      vp3            ! v'^3                        [(m/s)^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) ::  &
      up2,         & ! u'^2                        [(m/s)^2]
      upwp,        & ! u'w'                        [(m/s)^2]
      vp2,         & ! v'^2                        [(m/s)^2]
      vpwp           ! v'w'                        [(m/s)^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      wp2hmp,      & ! Third-order moment:  < w'^2 hm' >    [(m/s)^2 <hm units>]
      rtphmp_zt,   & ! Covariance of rt and hm (on t-levs.) [(kg/kg) <hm units>]
      thlphmp_zt     ! Covariance of thl and hm (on t-levs.)      [K <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance of w and a hydrometeor      [(m/s) <hm units>]

    ! Passive scalar variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim), intent(in) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      sclrp3       ! sclr'^3 (thermodynamic levels)       [{units vary}^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(in) :: &
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      p_sfc        ! Pressure at surface                  [Pa]

    logical, intent(in) :: &
      l_samp_stats_in_pdf_call    ! Sample stats in this call to this subroutine

    real( kind = core_rknd ), intent(in) :: &
      mixt_frac_max_mag, &    ! Maximum allowable mag. of mixt_frac   [-]
      ts_nudge, &             ! Timescale of u/v nudging             [s]
      rtm_min, &              ! Value below which rtm will be nudged [kg/kg]
      rtm_nudge_max_altitude  ! Highest altitude at which to nudge rtm [m]

    real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      iiPDF_type,      & ! Selected option for the two-component normal (double
                         ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                         ! w, chi, and eta) portion of CLUBB's multivariate,
                         ! two-component PDF.
      saturation_formula ! Integer that stores the saturation formula to be used

    logical, dimension(hydromet_dim), intent(in) :: &
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    logical, intent(in) :: &
      l_predict_upwp_vpwp,      & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                  ! alongside the advancement of <rt>, <w'rt'>, <thl>, <wpthlp>,
                                  ! <sclr>, and <w'sclr'> in subroutine advance_xm_wpxp.
                                  ! Otherwise, <u'w'> and <v'w'> are still approximated by eddy
                                  ! diffusivity when <u> and <v> are advanced in subroutine
                                  ! advance_windm_edsclrm.
      l_rtm_nudge,              & ! For rtm nudging
      l_trapezoidal_rule_zt,    & ! If true, the trapezoidal rule is called for the
                                  ! thermodynamic-level variables output from pdf_closure.
      l_trapezoidal_rule_zm,    & ! If true, the trapezoidal rule is called for three
                                  ! momentum-level variables – wpthvp, thlpthvp, and rtpthvp -
                                  ! output from pdf_closure.
      l_call_pdf_closure_twice, & ! This logical flag determines whether or not to call subroutine
                                  ! pdf_closure twice.  If true, pdf_closure is called first on
                                  ! thermodynamic levels and then on momentum levels so that each
                                  ! variable is computed on its native level.  If false,
                                  ! pdf_closure is only called on thermodynamic levels, and
                                  ! variables which belong on momentum levels are interpolated.
      l_use_cloud_cover,        & ! Use cloud_cover and rcm_in_layer to help boost cloud_frac and
                                  ! rcm to help increase cloudiness at coarser grid resolutions.
      l_rcm_supersat_adj          ! Add excess supersaturated vapor to cloud water

    type(stats_type), intent(inout) :: &
      stats

    !------------------------------- InOut Variables -------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) ::  &
      rtm    ! total water mixing ratio, r_t (thermo. levels) [kg/kg]

    type(implicit_coefs_terms), intent(inout) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Variable being passed back to and out of advance_clubb_core.
    type(pdf_parameter), intent(inout) :: &
      pdf_params,    & ! PDF parameters                           [units vary]
      pdf_params_zm    ! PDF parameters                           [units vary]

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), dimension(ngrdcol,nz, min(1,sclr_dim) , 2), intent(inout) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
! ---> h1g, 2012-06-14
    logical, intent(in)                 ::  do_liquid_only_in_clubb
! <--- h1g, 2012-06-14
#endif

    !------------------------------- Output Variables -------------------------------
    ! Variables being passed back to and out of advance_clubb_core.
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      rcm,               & ! mean r_c (thermodynamic levels)        [kg/kg]
      cloud_frac,        & ! cloud fraction (thermodynamic levels)  [-]
      ice_supersat_frac, & ! ice supersat. frac. (thermo. levels)   [-]
      wp2thvp,           & ! < w'^2 th_v' > (thermodynamic levels)  [m^2/s^2 K]
      wp2up,             & ! < w'^2 u' > (thermodynamic levels)     [m^3/s^3]
      rc_coef,           & ! Coefficient of X'r_c' (thermo. levs.)  [K/(kg/kg)]
      rcm_in_layer,      & ! rcm in cloud layer                     [kg/kg]
      cloud_cover,       & ! cloud cover                            [-]
      rcp2_zt              ! r_c'^2 (on thermo. grid)               [kg^2/kg^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) ::  &
      wprcp,             & ! < w'r_c' > (momentum levels)           [m/s kg/kg]
      sigma_sqd_w,       & ! PDF width parameter (momentum levels)  [-]
      wpthvp,            & ! < w' th_v' > (momentum levels)         [kg/kg K]
      rtpthvp,           & ! < r_t' th_v' > (momentum levels)       [kg/kg K]
      thlpthvp,          & ! < th_l' th_v' > (momentum levels)      [K^2]
      thlprcp,           & ! < th_l' r_c' > (momentum levels)       [K kg/kg]
      rc_coef_zm           ! Coefficient of X'r_c' on m-levs.       [K/(kg/kg)]

    ! Variable being passed back to and out of advance_clubb_core.
    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(out) :: &
      sclrpthvp    ! < sclr' th_v' > (momentum levels)   [units vary]

    ! Variables being passed back to only advance_clubb_core (for statistics).
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      wpup2,     & ! < w'u'^2 > (thermodynamic levels)        [m^3/s^3]
      wpvp2,     & ! < w'v'^2 > (thermodynamic levels)        [m^3/s^3]
      wp2rtp,    & ! < w'^2 r_t' > (thermodynamic levels)     [m^2/s^2 kg/kg]
      wprtp2,    & ! < w' r_t'^2 > (thermodynamic levels)     [m/s kg^2/kg^2]
      wp2thlp,   & ! < w'^2 th_l' > (thermodynamic levels)    [m^2/s^2 K]
      wpthlp2,   & ! < w' th_l'^2 > (thermodynamic levels)    [m/s K^2]
      wprtpthlp, & ! < w' r_t' th_l' > (thermodynamic levels) [m/s kg/kg K]
      wp2rcp       ! < w'^2 r_c' > (thermodynamic levels)     [m^2/s^2 kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) ::  &
      wp2up2,    & ! < w'^2u'^2 > (momentum levels)           [m^4/s^4]
      wp2vp2,    & ! < w'^2v'^2 > (momentum levels)           [m^4/s^4]
      wp4,       & ! < w'^4 > (momentum levels)               [m^4/s^4]
      rtprcp,    & ! < r_t' r_c' > (momentum levels)          [kg^2/kg^2]
      rcp2         ! Variance of r_c (momentum levels)        [kg^2/kg^2]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) ::  &
      uprcp,                 & ! < u' r_c' >                [(m kg)/(s kg)]
      vprcp                    ! < v' r_c' >                [(m kg)/(s kg)]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      w_up_in_cloud,         & ! mean cloudy updraft vel    [m/s]
      w_down_in_cloud,       & ! mean cloudy downdraft vel  [m/s]
      cloudy_updraft_frac,   & ! cloudy updraft fraction    [-]
      cloudy_downdraft_frac    ! cloudy downdraft fraction  [-]

    ! Variables being passed back to only advance_clubb_core (for statistics).
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) ::  &
      Skw_velocity,         & ! Skewness velocity                        [m/s]
      cloud_frac_zm,        & ! Cloud Fraction on momentum levels        [-]
      ice_supersat_frac_zm, & ! Ice supersat. frac. on momentum levels   [-]
      rtm_zm,               & ! Total water mixing ratio at mom. levs.   [kg/kg]
      thlm_zm,              & ! Liquid water pot. temp. at mom. levs.    [K]
      rcm_zm                  ! rcm at momentum levels                   [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
      rcm_supersat_adj        ! Adjust. to rcm due to spurious supersat. [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim), intent(out) :: &
      wp2sclrp,    & ! < w'^2 sclr' > (thermodynamic levels)      [units vary]
      wpsclrp2,    & ! < w' sclr'^2 > (thermodynamic levels)      [units vary]
      wpsclrprtp,  & ! < w' sclr' r_t' > (thermodynamic levels)   [units vary]
      wpsclrpthlp    ! < w' sclr' th_l' > (thermodynamic levels)  [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(out) :: &
      sclrprcp       ! < sclr' r_c' > (momentum levels)           [units vary]

    !------------------------------- Local Variables -------------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      wp2_zt,           & ! wp2 interpolated to thermodynamic levels   [m^2/s^2]
      rtp2_zt,          & ! rtp2 interpolated to thermodynamic levels  [kg^2/kg^2]
      thlp2_zt,         & ! thlp2 interpolated to thermodynamic levels [K^2]
      wprtp_zt,         & ! wprtp interpolated to thermodynamic levels [m/s kg/kg]
      wpthlp_zt,        & ! wpthlp interpolated to thermodynamic levs. [m/s K]
      rtpthlp_zt,       & ! rtpthlp interp. to thermodynamic levels    [kg/kg K]
      up2_zt,           & ! up2 interpolated to thermodynamic levels   [m^2/s^2]
      vp2_zt,           & ! vp2 interpolated to thermodynamic levels   [m^2/s^2]
      upwp_zt,          & ! upwp interpolated to thermodynamic levels  [m^2/s^2]
      vpwp_zt,          & ! vpwp interpolated to thermodynamic levels  [m^2/s^2]
      gamma_Skw_fnc_zt, & ! Gamma as a function of skewness (t-levs.)  [-]
      sigma_sqd_w_zt,   & ! PDF width parameter (thermodynamic levels) [-]
      Skw_zt,           & ! Skewness of w on thermodynamic levels      [-]
      Skrt_zt,          & ! Skewness of rt on thermodynamic levels     [-]
      Skthl_zt,         & ! Skewness of thl on thermodynamic levels    [-]
      Sku_zt,           & ! Skewness of u on thermodynamic levels      [-]
      Skv_zt              ! Skewness of v on thermodynamic levels      [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      wp3_zm,           & ! wp3 interpolated to momentum levels        [m^3/s^3]
      rtp3_zm,          & ! rtp3 interpolated to momentum levels       [kg^3/kg^3]
      thlp3_zm,         & ! thlp3 interpolated to momentum levels      [K^3]
      up3_zm,           & ! up3 interpolated to momentum levels        [m^3/s^3]
      vp3_zm,           & ! vp3 interpolated to momentum levels        [m^3/s^3]
      gamma_Skw_fnc,    & ! Gamma as a function of skewness            [-]
      Skw_zm,           & ! Skewness of w on momentum levels           [-]
      Skrt_zm,          & ! Skewness of rt on momentum levels          [-]
      Skthl_zm,         & ! Skewness of thl on momentum levels         [-]
      Sku_zm,           & ! Skewness of u on momentum levels           [-]
      Skv_zm              ! Skewness of v on momentum levels           [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      w_up_in_cloud_zm,         & ! Avg. cloudy updraft velocity; m-levs   [m/s]
      w_down_in_cloud_zm,       & ! Avg. cloudy downdraft velocity; m-levs [m/s]
      cloudy_updraft_frac_zm,   & ! cloudy updraft fraction; m-levs        [-]
      cloudy_downdraft_frac_zm    ! cloudy downdraft fraction; m-levs      [-]

    ! Interpolated values for optional second call to PDF closure.
    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      p_in_Pa_zm, & ! Pressure interpolated to momentum levels  [Pa]
      exner_zm      ! Exner interpolated to momentum levels     [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim) :: &
      wphydrometp_zt    ! Covariance of w and hm (on t-levs.) [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,hydromet_dim) :: &
      wp2hmp_zm,      & ! Moment <w'^2 hm'> (on m-levs.)    [(m/s)^2 <hm units>]
      rtphmp,         & ! Covariance of rt and hm           [(kg/kg) <hm units>]
      thlphmp           ! Covariance of thl and hm                [K <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim) :: &
      wpsclrp_zt,   & ! w' sclr' interpolated to thermo. levels
      sclrp2_zt,    & ! sclr'^2 interpolated to thermo. levels
      sclrprtp_zt,  & ! sclr' r_t' interpolated to thermo. levels
      sclrpthlp_zt, & ! sclr' th_l' interpolated thermo. levels
      Sksclr_zt       ! Skewness of sclr on thermodynamic levels      [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim) :: &
      sclrp3_zm,    & ! sclr'^3 interpolated to momentum levels
      Sksclr_zm       ! Skewness of sclr on momentum levels           [-]

    ! These local variables are declared because they originally belong on the
    ! momentum grid levels, but pdf_closure outputs them on the thermodynamic
    ! grid levels.
    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      wpup2_zm,    & ! w'u'^2 (on momentum grid)        [m^3/s^3]
      wpvp2_zm       ! w'v'^2 (on momentum grid)        [m^3/s^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      wp2up2_zt,   & ! w'^2u'^2 (on thermo. grid)       [m^4/s^4]
      wp2vp2_zt,   & ! w'^2v'^2 (on thermo. grid)       [m^4/s^4]
      wp4_zt,      & ! w'^4 (on thermo. grid)           [m^4/s^4]
      wpthvp_zt,   & ! Buoyancy flux (on thermo. grid)  [(K m)/s]
      rtpthvp_zt,  & ! r_t' th_v' (on thermo. grid)     [(kg K)/kg]
      thlpthvp_zt, & ! th_l' th_v' (on thermo. grid)    [K^2]
      wprcp_zt,    & ! w' r_c' (on thermo. grid)        [(m kg)/(s kg)]
      rtprcp_zt,   & ! r_t' r_c' (on thermo. grid)      [(kg^2)/(kg^2)]
      thlprcp_zt,  & ! th_l' r_c' (on thermo. grid)     [(K kg)/kg]
      uprcp_zt,    & ! u' r_c' (on thermo. grid)        [(m kg)/(s kg)]
      vprcp_zt       ! v' r_c' (on thermo. grid)        [(m kg)/(s kg)]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim) :: &
      sclrpthvp_zt, & ! sclr'th_v' (on thermo. grid)
      sclrprcp_zt     ! sclr'rc' (on thermo. grid)

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      wprtp2_zm,    & ! < w' r_t'^2 > on momentum levels      [m/s kg^2/kg^2]
      wp2rtp_zm,    & ! < w'^2 r_t' > on momentum levels      [m^2/s^2 kg/kg]
      wpthlp2_zm,   & ! < w' th_l'^2 > on momentum levels     [m/s K^2]
      wp2thlp_zm,   & ! < w'^2 th_l' > on momentum levels     [m^2/s^2 K]
      wprtpthlp_zm, & ! < w' r_t' th_l' > on momentum levels  [m/s kg/kg K]
      wp2thvp_zm,   & ! < w'^2 th_v' > on momentum levels     [m^2/s^2 K]
      wp2up_zm,     & ! < w'^2 u' > on momentum levels        [m^3/s^3]
      wp2rcp_zm       ! < w'^2 r_c' > on momentum levles      [m^2/s^2 kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim) :: &
      wpsclrprtp_zm,  & ! w'sclr'rt' on momentum grid
      wpsclrp2_zm,    & ! w'sclr'^2 on momentum grid
      wpsclrpthlp_zm, & ! w'sclr'thl' on momentum grid
      wp2sclrp_zm,    & ! w'^2 sclr' on momentum grid
      sclrm_zm          ! Passive scalar mean on momentum grid

    type(implicit_coefs_terms) :: &
      pdf_implicit_coefs_terms_zm

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rsat,             & ! Saturation mixing ratio from mean rt and thl.
      rel_humidity        ! Relative humidity after PDF closure [-]

    real( kind = core_rknd ) :: &
      gamma_coef,     & ! CLUBB tunable parameter gamma_coef
      gamma_coefb,    & ! CLUBB tunable parameter gamma_coefb
      gamma_coefc       ! CLUBB tunable parameter gamma_coefc
      
    real( kind = core_rknd ), dimension(ngrdcol, nzm) :: &
      um_zm, &
      vm_zm, &
      sigma_sqd_w_tmp

    real( kind = core_rknd ), dimension(ngrdcol, nzt) :: &
      T_in_K

    logical :: l_spur_supersat   ! Spurious supersaturation?

    integer :: i, k, j, sclr

    !-------------------------------------- Begin Code --------------------------------------

    !$acc enter data create( wp2_zt,wp3_zm, rtp2_zt,rtp3_zm, thlp2_zt,  thlp3_zm, &
    !$acc                    wprtp_zt, wpthlp_zt, rtpthlp_zt, up2_zt, up3_zm, &
    !$acc                    vp2_zt, vp3_zm, upwp_zt, vpwp_zt, gamma_Skw_fnc, &
    !$acc                    gamma_Skw_fnc_zt,sigma_sqd_w_zt,  Skw_zt, Skw_zm, &
    !$acc                    Skrt_zt, Skrt_zm, Skthl_zt, Skthl_zm, Sku_zt, &
    !$acc                    Sku_zm, Skv_zt, Skv_zm, wp2up2_zt, &
    !$acc                    wp2vp2_zt, wp4_zt, wpthvp_zt, rtpthvp_zt, thlpthvp_zt, &
    !$acc                    wprcp_zt, rtprcp_zt, thlprcp_zt, uprcp_zt, vprcp_zt, &
    !$acc                    rsat, rel_humidity, um_zm, vm_zm, T_in_K, sigma_sqd_w_tmp )

    !$acc enter data if( l_call_pdf_closure_twice ) &
    !$acc            create( w_up_in_cloud_zm, wpup2_zm, wpvp2_zm, &
    !$acc                    w_down_in_cloud_zm, cloudy_updraft_frac_zm,  &
    !$acc                    cloudy_downdraft_frac_zm, p_in_Pa_zm, exner_zm, &
    !$acc                    wprtp2_zm, wp2rtp_zm, wpthlp2_zm, &
    !$acc                    wp2thlp_zm, wprtpthlp_zm, wp2thvp_zm, wp2rcp_zm, wp2up_zm )

    !$acc enter data if( sclr_dim > 0 ) &
    !$acc            create( wpsclrp_zt, sclrp2_zt, sclrp3_zm, sclrprtp_zt, sclrpthlp_zt, &
    !$acc                    Sksclr_zt, Sksclr_zm, sclrpthvp_zt, sclrprcp_zt, wpsclrprtp_zm, &
    !$acc                    wpsclrp2_zm, wpsclrpthlp_zm, wp2sclrp_zm, sclrm_zm )

    !$acc enter data if( hydromet_dim > 0 ) create( wphydrometp_zt, wp2hmp_zm, rtphmp, thlphmp )

    !---------------------------------------------------------------------------
    ! Interpolate wp3, rtp3, thlp3, up3, vp3, and sclrp3 to momentum levels, and
    ! wp2, rtp2, thlp2, up2, vp2, and sclrp2 to thermodynamic levels, and then
    ! compute Skw, Skrt, Skthl, Sku, Skv, and Sksclr for both the momentum and
    ! thermodynamic grid levels.
    !---------------------------------------------------------------------------

    ! Positive definite quantity
    wp2_zt(:,:)   = zm2zt_api( nzm, nzt, ngrdcol, gr, wp2(:,:), w_tol_sqd )
    wp3_zm(:,:)   = zt2zm_api( nzm, nzt, ngrdcol, gr, wp3(:,:) )
    
    ! Positive definite quantity
    thlp2_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, thlp2(:,:), thl_tol**2 )
    thlp3_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, thlp3(:,:) )

    ! Positive definite quantity
    rtp2_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, rtp2(:,:), rt_tol**2 )
    rtp3_zm(:,:)  = zt2zm_api( nzm, nzt, ngrdcol, gr, rtp3(:,:) )

    ! Positive definite quantity
    up2_zt(:,:)   = zm2zt_api( nzm, nzt, ngrdcol, gr, up2(:,:), w_tol_sqd )
    up3_zm(:,:)   = zt2zm_api( nzm, nzt, ngrdcol, gr, up3(:,:) )

    ! Positive definite quantity
    vp2_zt(:,:)   = zm2zt_api( nzm, nzt, ngrdcol, gr, vp2(:,:), w_tol_sqd )
    vp3_zm(:,:)   = zt2zm_api( nzm, nzt, ngrdcol, gr, vp3(:,:) )

    do sclr = 1, sclr_dim, 1
      sclrp2_zt(:,:,sclr) = zm2zt_api( nzm, nzt, ngrdcol, gr, sclrp2(:,:,sclr), &
                                       sclr_tol(sclr)**2 ) ! Pos. def. quantity
      sclrp3_zm(:,:,sclr) = zt2zm_api( nzm, nzt, ngrdcol, gr, sclrp3(:,:,sclr) )
    end do ! sclr = 1, sclr_dim, 1

    call Skx_func( nzt, ngrdcol, wp2_zt, wp3, &
                   w_tol, clubb_params, &
                   Skw_zt )
                   
    call Skx_func( nzm, ngrdcol, wp2, wp3_zm, &
                   w_tol, clubb_params, &
                   Skw_zm )    
                   
    call Skx_func( nzt, ngrdcol, thlp2_zt, thlp3, &
                   thl_tol, clubb_params, &
                   Skthl_zt )  
                   
    call Skx_func( nzm, ngrdcol, thlp2, thlp3_zm, &
                   thl_tol, clubb_params, &
                   Skthl_zm )  
                   
    call Skx_func( nzt, ngrdcol, rtp2_zt, rtp3, &
                   rt_tol, clubb_params, &
                   Skrt_zt )   
                   
    call Skx_func( nzm, ngrdcol, rtp2, rtp3_zm, &
                   rt_tol, clubb_params, &
                   Skrt_zm )   
                   
    call Skx_func( nzt, ngrdcol, up2_zt, up3, &
                   w_tol, clubb_params, &
                   Sku_zt )   
                                      
    call Skx_func( nzm, ngrdcol, up2, up3_zm, &
                   w_tol, clubb_params, &
                   Sku_zm )   
                   
    call Skx_func( nzt, ngrdcol, vp2_zt, vp3, &
                   w_tol, clubb_params, &
                   Skv_zt )   
                   
    call Skx_func( nzm, ngrdcol, vp2, vp3_zm, &
                   w_tol, clubb_params, &
                   Skv_zm )      

    do sclr = 1, sclr_dim
      
      call Skx_func( nzt, ngrdcol, sclrp2_zt(:,:,sclr), sclrp3(:,:,sclr), &
                     sclr_tol(sclr), clubb_params, &
                     Sksclr_zt(:,:,sclr) )   
                     
      call Skx_func( nzm, ngrdcol, sclrp2(:,:,sclr), sclrp3_zm(:,:,sclr), &
                     sclr_tol(sclr), clubb_params, &
                     Sksclr_zm(:,:,sclr) )  
                      
    end do ! sclr = 1, sclr_dim, 1

    if ( stats%l_sample .and. l_samp_stats_in_pdf_call ) then
      !$acc update host( Skw_zt, Skw_zm, Skthl_zt, Skrt_zt, Skrt_zm, Skthl_zm )
      call stats_update( "Skw_zt", Skw_zt, stats )
      call stats_update( "Skw_zm", Skw_zm, stats )
      call stats_update( "Skthl_zt", Skthl_zt, stats )
      call stats_update( "Skthl_zm", Skthl_zm, stats )
      call stats_update( "Skrt_zt", Skrt_zt, stats )
      call stats_update( "Skrt_zm", Skrt_zm, stats )
    end if

    ! The right hand side of this conjunction is only for reducing cpu time,
    ! since the more complicated formula is mathematically equivalent
    if ( l_gamma_Skw ) then

      !----------------------------------------------------------------
      ! Compute gamma as a function of Skw  - 14 April 06 dschanen
      !----------------------------------------------------------------
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol

          gamma_coef  = clubb_params(i,igamma_coef)
          gamma_coefb = clubb_params(i,igamma_coefb)
          gamma_coefc = clubb_params(i,igamma_coefc)

          if ( abs( gamma_coef - gamma_coefb ) > abs( gamma_coef + gamma_coefb ) * eps/2 ) then

            gamma_Skw_fnc(i,k) = gamma_coefb &
                                 + ( gamma_coef - gamma_coefb ) &
                                   * exp( -one_half * ( Skw_zm(i,k) / gamma_coefc )**2 )

          else

            gamma_Skw_fnc(i,k) = gamma_coef

          end if

        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol

          gamma_coef  = clubb_params(i,igamma_coef)
          gamma_coefb = clubb_params(i,igamma_coefb)
          gamma_coefc = clubb_params(i,igamma_coefc)

          if ( abs( gamma_coef - gamma_coefb ) > abs( gamma_coef + gamma_coefb ) * eps/2 ) then

            gamma_Skw_fnc_zt(i,k) = gamma_coefb &
                                    + ( gamma_coef - gamma_coefb ) &
                                      * exp( -one_half * ( Skw_zt(i,k) / gamma_coefc )**2 )

          else

            gamma_Skw_fnc_zt(i,k) = gamma_coef

          end if

        end do
      end do
      !$acc end parallel loop

    else
      
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          gamma_Skw_fnc(i,k) = clubb_params(i,igamma_coef)
        end do
      end do
      !$acc end parallel loop

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          gamma_Skw_fnc_zt(i,k) = clubb_params(i,igamma_coef)
        end do
      end do
      !$acc end parallel loop

    end if

    if ( stats%l_sample .and. l_samp_stats_in_pdf_call ) then
      !$acc update host(gamma_Skw_fnc)
      call stats_update( "gamma_Skw_fnc", gamma_Skw_fnc, stats )
    end if

    ! Compute sigma_sqd_w (dimensionless PDF width parameter)
    call compute_sigma_sqd_w( nzm, ngrdcol, &
                              gamma_Skw_fnc, wp2, thlp2, rtp2, &
                              up2, vp2, wpthlp, wprtp, upwp, vpwp, &
                              l_predict_upwp_vpwp, &
                              sigma_sqd_w_tmp )

    ! Smooth in the vertical using interpolation
    sigma_sqd_w(:,:) = zm2zt2zm( nzm, nzt, ngrdcol, gr, sigma_sqd_w_tmp(:,:), &
                                 zero_threshold ) ! Pos. def. quantity


    ! Interpolate the the stats_zt grid
    ! Pos. def. quantity
    sigma_sqd_w_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, sigma_sqd_w(:,:), zero_threshold )

    !---------------------------------------------------------------------------
    ! Interpolate thlp2, rtp2, and rtpthlp to thermodynamic levels,
    !---------------------------------------------------------------------------

    ! Interpolate variances to the stats_zt grid (statistics and closure)
    ! Positive def. quantity
    rtp2_zt(:,:)    = zm2zt_api( nzm, nzt, ngrdcol, gr, rtp2(:,:), rt_tol**2 )
    ! Positive def. quantity
    thlp2_zt(:,:)   = zm2zt_api( nzm, nzt, ngrdcol, gr, thlp2(:,:), thl_tol**2 )
    ! Positive def. quantity
    up2_zt(:,:)     = zm2zt_api( nzm, nzt, ngrdcol, gr, up2(:,:), w_tol_sqd )
    ! Positive def. quantity
    vp2_zt(:,:)     = zm2zt_api( nzm, nzt, ngrdcol, gr, vp2(:,:), w_tol_sqd )
    wprtp_zt(:,:)   = zm2zt_api( nzm, nzt, ngrdcol, gr, wprtp(:,:) )
    wpthlp_zt(:,:)  = zm2zt_api( nzm, nzt, ngrdcol, gr, wpthlp(:,:) )
    rtpthlp_zt(:,:) = zm2zt_api( nzm, nzt, ngrdcol, gr, rtpthlp(:,:) )
    upwp_zt(:,:)    = zm2zt_api( nzm, nzt, ngrdcol, gr, upwp(:,:) )
    vpwp_zt(:,:)    = zm2zt_api( nzm, nzt, ngrdcol, gr, vpwp(:,:) )

    ! Compute skewness velocity for stats output purposes
    if ( var_on_stats_list( stats, "Skw_velocity" ) ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          Skw_velocity(i,k) = ( 1.0_core_rknd / ( 1.0_core_rknd - sigma_sqd_w(i,k) ) ) &
                       * ( wp3_zm(i,k) / max( wp2(i,k), w_tol_sqd ) )
        end do
      end do
      !$acc end parallel loop
    end if

    !----------------------------------------------------------------
    ! Call closure scheme
    !----------------------------------------------------------------

    ! Put passive scalar input on the t grid for the PDF
    do sclr = 1, sclr_dim
      wpsclrp_zt(:,:,sclr)   = zm2zt_api( nzm, nzt, ngrdcol, gr, wpsclrp(:,:,sclr) )
      sclrp2_zt(:,:,sclr)    = zm2zt_api( nzm, nzt, ngrdcol, gr, sclrp2(:,:,sclr), &
                                          sclr_tol(sclr)**2 ) ! Pos. def. quantity
      sclrprtp_zt(:,:,sclr)  = zm2zt_api( nzm, nzt, ngrdcol, gr, sclrprtp(:,:,sclr) )
      sclrpthlp_zt(:,:,sclr) = zm2zt_api( nzm, nzt, ngrdcol, gr, sclrpthlp(:,:,sclr) )
    end do ! sclr = 1, sclr_dim, 1

    ! Interpolate hydrometeor mixed moments to momentum levels.
    do j = 1, hydromet_dim
      wphydrometp_zt(:,:,j) = zm2zt_api( nzm, nzt, ngrdcol, gr, wphydrometp(:,:,j) )
    end do ! i = 1, hydromet_dim, 1

    call pdf_closure( nzt, ngrdcol, sclr_dim, sclr_tol, gr, & ! intent(in)
           hydromet_dim, p_in_Pa, exner, thv_ds_zt,         & ! intent(in)
           wm_zt, wp2_zt, wp3,                              & ! intent(in)
           Skw_zt, Skthl_zt, Skrt_zt, Sku_zt, Skv_zt,       & ! intent(in)
           rtm, rtp2_zt, wprtp_zt,                          & ! intent(in)
           thlm, thlp2_zt, wpthlp_zt,                       & ! intent(in)
           um, up2_zt, upwp_zt,                             & ! intent(in)
           vm, vp2_zt, vpwp_zt,                             & ! intent(in)
           rtpthlp_zt,                                      & ! intent(in)
           sclrm, wpsclrp_zt, sclrp2_zt,                    & ! intent(in)
           sclrprtp_zt, sclrpthlp_zt, Sksclr_zt,            & ! intent(in)
           gamma_Skw_fnc_zt,                                & ! intent(in)
#ifdef GFDL
           RH_crit,                                         & ! intent(inout)
           do_liquid_only_in_clubb,                         & ! intent(in)
#endif

           wphydrometp_zt, wp2hmp,                          & ! intent(in)
           rtphmp_zt, thlphmp_zt,                           & ! intent(in)
           clubb_params, mixt_frac_max_mag,                 & ! intent(in)
           saturation_formula,                              & ! intent(in)
           stats,                                           & ! intent(inout)
           iiPDF_type,                                      & ! intent(in)
           l_mix_rat_hm,                                    & ! intent(in)
           sigma_sqd_w_zt,                                  & ! intent(inout)
           pdf_params, pdf_implicit_coefs_terms,            & ! intent(inout)
           err_info,                                        & ! intent(inout)
           wpup2, wpvp2,                                    & ! intent(out)
           wp2up2_zt, wp2vp2_zt, wp4_zt,                    & ! intent(out)
           wprtp2, wp2rtp,                                  & ! intent(out)
           wpthlp2, wp2thlp, wprtpthlp,                     & ! intent(out)
           cloud_frac, ice_supersat_frac,                   & ! intent(out)
           rcm, wpthvp_zt, wp2thvp, wp2up, rtpthvp_zt,      & ! intent(out)
           thlpthvp_zt, wprcp_zt, wp2rcp, rtprcp_zt,        & ! intent(out)
           thlprcp_zt, rcp2_zt,                             & ! intent(out)
           uprcp_zt, vprcp_zt,                              & ! intent(out)
           w_up_in_cloud, w_down_in_cloud,                  & ! intent(out)
           cloudy_updraft_frac, cloudy_downdraft_frac,      & ! intent(out)
           wpsclrprtp, wpsclrp2, sclrpthvp_zt,              & ! intent(out)
           wpsclrpthlp, sclrprcp_zt, wp2sclrp,              & ! intent(out)
           rc_coef                                          ) ! intent(out)

    ! Subroutine may produce NaN values, and if so, return
    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "After first call to pdf_closure in pdf_closure_driver"
        return
      endif
    endif

    if( l_rtm_nudge ) then
      ! Nudge rtm to prevent excessive drying
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          if ( rtm(i,k) < rtm_min .and. gr%zt(i,k) < rtm_nudge_max_altitude ) then
            rtm(i,k) = rtm(i,k) + (rtm_ref(i,k) - rtm(i,k)) * ( dt / ts_nudge )
          end if
        end do
      end do
      !$acc end parallel loop
    end if

    if ( l_call_pdf_closure_twice ) then

      ! Call pdf_closure a second time on momentum levels, to
      ! output (rather than interpolate) the variables which
      ! belong on the momentum levels.

      ! Interpolate sclrm to the momentum level for use in
      ! the second call to pdf_closure
      do sclr = 1, sclr_dim
        ! Clip if extrap. causes sclrm_zm to be less than sclr_tol
        sclrm_zm(:,:,sclr) = zt2zm_api( nzm, nzt, ngrdcol, gr, sclrm(:,:,sclr), sclr_tol(sclr) )
      end do ! sclr = 1, sclr_dim

      ! Interpolate pressure, p_in_Pa, to momentum levels.
      ! Since the surface (or model lower boundary) is located at momentum level
      ! k = 1, the pressure there is p_sfc.
      p_in_Pa_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, p_in_Pa(:,:) )

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        p_in_Pa_zm(i,gr%k_lb_zm) = p_sfc(i)

        ! Clip pressure if the extrapolation leads to a negative value of pressure
        p_in_Pa_zm(i,gr%k_ub_zm) &
        = max( p_in_Pa_zm(i,gr%k_ub_zm), 0.5_core_rknd * p_in_Pa(i,gr%k_ub_zt) )
      end do
      !$acc end parallel loop

      ! Set exner at momentum levels, exner_zm, based on p_in_Pa_zm.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          exner_zm(i,k) = (p_in_Pa_zm(i,k)/p0)**kappa
        end do
      end do
      !$acc end parallel loop

      ! Clip if extrapolation at the top level causes rtm_zm to be < rt_tol
      rtm_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, rtm(:,:), rt_tol )

      ! Clip if extrapolation at the top level causes thlm_zm to be < thl_tol
      thlm_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, thlm(:,:), thl_tol )

      ! Interpolate hydrometeor mixed moments to momentum levels.
      do j = 1, hydromet_dim
        rtphmp(:,:,j)    = zt2zm_api( nzm, nzt, ngrdcol, gr, rtphmp_zt(:,:,j) )
        thlphmp(:,:,j)   = zt2zm_api( nzm, nzt, ngrdcol, gr, thlphmp_zt(:,:,j) )
        wp2hmp_zm(:,:,j) = zt2zm_api( nzm, nzt, ngrdcol, gr, wp2hmp(:,:,j) )
      end do ! i = 1, hydromet_dim, 1

      um_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, um(:,:) )
      vm_zm(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, vm(:,:) )
      
      ! pdf_implicit_coefs_terms is only used in the iiPDF_new and iiPDF_new_hybrid closures.
      ! So we only need to initialize our local _zm version if we're working with one of those.
      if ( iiPDF_type == iiPDF_new .or. iiPDF_type == iiPDF_new_hybrid ) then
        call init_pdf_implicit_coefs_terms_api( nzm, ngrdcol, sclr_dim, &     ! Intent(in)
                                                pdf_implicit_coefs_terms_zm ) ! Intent(out)
      end if

      ! Call pdf_closure to output the variables which belong on the momentum grid.
      call pdf_closure( nzm, ngrdcol, sclr_dim, sclr_tol, gr,      & ! intent(in)
             hydromet_dim, p_in_Pa_zm, exner_zm, thv_ds_zm,        & ! intent(in)
             wm_zm, wp2, wp3_zm,                                   & ! intent(in)
             Skw_zm, Skthl_zm, Skrt_zm, Sku_zm, Skv_zm,            & ! intent(in)
             rtm_zm, rtp2, wprtp,                                  & ! intent(in)
             thlm_zm, thlp2, wpthlp,                               & ! intent(in)
             um_zm, up2, upwp,                                     & ! intent(in)
             vm_zm, vp2, vpwp,                                     & ! intent(in)
             rtpthlp,                                              & ! intent(in)
             sclrm_zm, wpsclrp, sclrp2,                            & ! intent(in)
             sclrprtp, sclrpthlp, Sksclr_zm,                       & ! intent(in)
             gamma_Skw_fnc,                                        & ! intent(in)
#ifdef GFDL
             RH_crit,                                              & ! intent(inout)
             do_liquid_only_in_clubb,                              & ! intent(in)
#endif
             wphydrometp, wp2hmp_zm,                               & ! intent(in)
             rtphmp, thlphmp,                                      & ! intent(in)
             clubb_params, mixt_frac_max_mag,                      & ! intent(in)
             saturation_formula,                                   & ! intent(in)
             stats,                                                & ! intent(inout)
             iiPDF_type,                                           & ! intent(in)
             l_mix_rat_hm,                                         & ! intent(in)
             sigma_sqd_w,                                          & ! intent(inout)
             pdf_params_zm, pdf_implicit_coefs_terms_zm,           & ! intent(inout)
             err_info,                                             & ! intent(inout)
             wpup2_zm, wpvp2_zm,                                   & ! intent(out)
             wp2up2, wp2vp2, wp4,                                  & ! intent(out)
             wprtp2_zm, wp2rtp_zm,                                 & ! intent(out)
             wpthlp2_zm, wp2thlp_zm, wprtpthlp_zm,                 & ! intent(out)
             cloud_frac_zm, ice_supersat_frac_zm,                  & ! intent(out)
             rcm_zm, wpthvp, wp2thvp_zm, wp2up_zm, rtpthvp,        & ! intent(out)
             thlpthvp, wprcp, wp2rcp_zm, rtprcp,                   & ! intent(out)
             thlprcp, rcp2,                                        & ! intent(out)
             uprcp, vprcp,                                         & ! intent(out)
             w_up_in_cloud_zm, w_down_in_cloud_zm,                 & ! intent(out)
             cloudy_updraft_frac_zm, cloudy_downdraft_frac_zm,     & ! intent(out)
             wpsclrprtp_zm, wpsclrp2_zm, sclrpthvp,                & ! intent(out)
             wpsclrpthlp_zm, sclrprcp, wp2sclrp_zm,                & ! intent(out)
             rc_coef_zm                                            ) ! intent(out)

      ! Subroutine may produce NaN values, and if so, return
      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr,*) "After second call to pdf_closure in pdf_closure_driver"
          return
        endif
      endif

    else ! l_call_pdf_closure_twice is false
      
      ! Interpolate momentum variables output from the first call to
      ! pdf_closure back to momentum grid.
      ! Pos. def. quantity
      wp4(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, wp4_zt(:,:), zero_threshold )

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        ! Since top momentum level is higher than top thermo level,
        ! set variables at top momentum level to 0.
        wp4(i,gr%k_ub_zm) = zero
        ! Set wp4 to 0 at the lowest momentum level (momentum level 1).
        wp4(i,gr%k_lb_zm) = zero
      end do
      !$acc end parallel loop

#ifndef CLUBB_CAM
      ! CAM-CLUBB needs cloud water variance thus always compute this
      if ( var_on_stats_list( stats, "rcp2" ) ) then
#endif
        ! Pos. def. quantity
        rcp2(:,:) = zt2zm_api( nzm, nzt, ngrdcol, gr, rcp2_zt(:,:), zero_threshold )
#ifndef CLUBB_CAM
        !$acc parallel loop gang vector default(present) 
        do i = 1, ngrdcol
          rcp2(i,gr%k_ub_zm) = zero
        end do
        !$acc end parallel loop
      endif
#endif

      wpthvp(:,:)      = zt2zm_api( nzm, nzt, ngrdcol, gr, wpthvp_zt(:,:) )
      thlpthvp(:,:)    = zt2zm_api( nzm, nzt, ngrdcol, gr, thlpthvp_zt(:,:) )
      rtpthvp(:,:)     = zt2zm_api( nzm, nzt, ngrdcol, gr, rtpthvp_zt(:,:) )
      wprcp(:,:)       = zt2zm_api( nzm, nzt, ngrdcol, gr, wprcp_zt(:,:) )
      rc_coef_zm(:,:)  = zt2zm_api( nzm, nzt, ngrdcol, gr, rc_coef(:,:) )
      rtprcp(:,:)      = zt2zm_api( nzm, nzt, ngrdcol, gr, rtprcp_zt(:,:) )
      thlprcp(:,:)     = zt2zm_api( nzm, nzt, ngrdcol, gr, thlprcp_zt(:,:) )
      uprcp(:,:)       = zt2zm_api( nzm, nzt, ngrdcol, gr, uprcp_zt(:,:) )
      vprcp(:,:)       = zt2zm_api( nzm, nzt, ngrdcol, gr, vprcp_zt(:,:) )
      wp2up2(:,:)      = zt2zm_api( nzm, nzt, ngrdcol, gr, wp2up2_zt(:,:) )
      wp2vp2(:,:)      = zt2zm_api( nzm, nzt, ngrdcol, gr, wp2vp2_zt(:,:) )

      !$acc parallel loop gang vector default(present) 
      do i = 1, ngrdcol 
        wpthvp(i,gr%k_ub_zm)     = 0.0_core_rknd
        thlpthvp(i,gr%k_ub_zm)   = 0.0_core_rknd
        rtpthvp(i,gr%k_ub_zm)    = 0.0_core_rknd
        wprcp(i,gr%k_ub_zm)      = 0.0_core_rknd
        rc_coef_zm(i,gr%k_ub_zm) = 0.0_core_rknd
        rtprcp(i,gr%k_ub_zm)     = 0.0_core_rknd
        thlprcp(i,gr%k_ub_zm)    = 0.0_core_rknd
        uprcp(i,gr%k_ub_zm)      = 0.0_core_rknd
        vprcp(i,gr%k_ub_zm)      = 0.0_core_rknd
        wp2up2(i,gr%k_ub_zm)     = 0.0_core_rknd
        wp2vp2(i,gr%k_ub_zm)     = 0.0_core_rknd
      end do
      !$acc end parallel loop

      ! Initialize variables to avoid uninitialized variables.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzm
        do i = 1, ngrdcol
          cloud_frac_zm(i,k)        = 0.0_core_rknd
          ice_supersat_frac_zm(i,k) = 0.0_core_rknd
          rcm_zm(i,k)               = 0.0_core_rknd
          rtm_zm(i,k)               = 0.0_core_rknd
          thlm_zm(i,k)              = 0.0_core_rknd
        end do
      end do
      !$acc end parallel loop

      ! Interpolate passive scalars back onto the m grid
      do sclr = 1, sclr_dim
        sclrpthvp(:,:,sclr)       = zt2zm_api( nzm, nzt, ngrdcol, gr, sclrpthvp_zt(:,:,sclr) )
        sclrprcp(:,:,sclr)        = zt2zm_api( nzm, nzt, ngrdcol, gr, sclrprcp_zt(:,:,sclr) )

        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          sclrpthvp(i,gr%k_ub_zm,sclr) = 0.0_core_rknd
          sclrprcp(i,gr%k_ub_zm,sclr)  = 0.0_core_rknd
        end do
        !$acc end parallel loop

      end do ! sclr=1, sclr_dim

    end if ! l_call_pdf_closure_twice
    
    if ( stats%l_sample .and. l_samp_stats_in_pdf_call ) then
      !$acc update host( uprcp, vprcp )
      call stats_update( "uprcp", uprcp, stats )
      call stats_update( "vprcp", vprcp, stats )
    end if
    
    ! If l_trapezoidal_rule_zt is true, call trapezoidal_rule_zt for
    ! thermodynamic-level variables output from pdf_closure.
    ! ldgrant June 2009
    if ( l_trapezoidal_rule_zt ) then
      call trapezoidal_rule_zt( nzm, nzt, ngrdcol, sclr_dim, gr,             & ! intent(in)
                                l_call_pdf_closure_twice,                    & ! intent(in)
                                stats,                                       & ! intent(in)
                                wprtp2, wpthlp2,                             & ! intent(inout)
                                wprtpthlp, cloud_frac, ice_supersat_frac,    & ! intent(inout)
                                rcm, wp2thvp, wp2up, wpsclrprtp, wpsclrp2,   & ! intent(inout)
                                wpsclrpthlp,                                 & ! intent(inout)
                                wprtp2_zm, wpthlp2_zm,                       & ! intent(inout)
                                wprtpthlp_zm, cloud_frac_zm,                 & ! intent(inout)
                                ice_supersat_frac_zm, rcm_zm, wp2thvp_zm,    & ! intent(inout)
                                wp2up_zm,                                    & ! intent(inout)
                                wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm )   ! intent(inout)
    else ! l_trapezoidal_rule_zt
      cloud_frac_zm = zt2zm_api( nzm, nzt, ngrdcol, gr, cloud_frac )
      ! Since top momentum level is higher than top thermo. level,
      ! set variables at top momentum level to 0.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        cloud_frac_zm(i,gr%k_ub_zm)  = 0.0_core_rknd
      end do
      !$acc end parallel loop
    end if ! l_trapezoidal_rule_zt

    ! If l_trapezoidal_rule_zm is true, call trapezoidal_rule_zm for
    ! the important momentum-level variabes output from pdf_closure.
    ! ldgrant Feb. 2010
    if ( l_trapezoidal_rule_zm ) then
      call trapezoidal_rule_zm( nzm, nzt, ngrdcol, gr,              & ! intent(in)
                                wpthvp_zt, thlpthvp_zt, rtpthvp_zt, & ! intent(in)
                                wpthvp, thlpthvp, rtpthvp )           ! intent(inout)
    end if ! l_trapezoidal_rule_zm


    ! Vince Larson clipped rcm in order to prevent rvm < 0.  5 Apr 2008.
    ! This code won't work unless rtm >= 0 !!!
    ! We do not clip rcm_in_layer because rcm_in_layer only influences
    ! radiation, and we do not want to bother recomputing it.
    ! Code is duplicated from below to ensure that relative humidity
    ! is calculated properly.  3 Sep 2009
    call clip_rcm( nzt, ngrdcol, rtm,             & ! intent(in)
                   'rtm < rcm after pdf_closure', & ! intent(in)
                   rcm )                            ! intent(inout)

    ! Compute variables cloud_cover and rcm_in_layer.
    ! Added July 2009
    call compute_cloud_cover( gr, nzt, ngrdcol,            & ! intent(in)
                              pdf_params, cloud_frac, rcm, & ! intent(in)
                              err_info,                    & ! intent(inout)
                              cloud_cover, rcm_in_layer )    ! intent(out)

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "calling compute_cloud_cover in pdf_closure_driver"
        return
      endif
    endif

    if ( l_use_cloud_cover ) then
      ! Use cloud_cover and rcm_in_layer to help boost cloud_frac and rcm to help
      ! increase cloudiness at coarser grid resolutions.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          cloud_frac(i,k) = cloud_cover(i,k)
          rcm(i,k) = rcm_in_layer(i,k)
        end do
      end do
    !$acc end parallel loop
    end if

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        ! Clip cloud fraction here if it still exceeds 1.0 due to round off
        cloud_frac(i,k) = min( 1.0_core_rknd, cloud_frac(i,k) )
        ! Ditto with ice cloud fraction
        ice_supersat_frac(i,k) = min( 1.0_core_rknd, ice_supersat_frac(i,k) )
      end do
    end do
    !$acc end parallel loop

    T_in_K = thlm2T_in_K_api( nzt, ngrdcol, thlm, exner, rcm )
    rsat = sat_mixrat_liq_api( nzt, ngrdcol, gr, p_in_Pa, T_in_K, saturation_formula )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        rel_humidity(i,k) = (rtm(i,k) - rcm(i,k)) / rsat(i,k)
        rcm_supersat_adj(i,k) = zero
      end do
    end do
    !$acc end parallel loop
      
    if ( l_rcm_supersat_adj ) then
      ! +PAB mods, take remaining supersaturation that may exist
      !   after CLUBB PDF call and add it to rcm.  Supersaturation
      !   may exist after PDF call due to issues with calling PDF on the
      !   thermo grid and momentum grid and the interpolation between the two
      l_spur_supersat = .false.
      

      !$acc parallel loop gang vector collapse(2) default(present) reduction(.or.:l_spur_supersat)
      do k = 1, nzt
        do i = 1, ngrdcol
          if (rel_humidity(i,k) > 1.0_core_rknd) then
            rcm_supersat_adj(i,k) = (rtm(i,k) - rcm(i,k)) - rsat(i,k)
            rcm(i,k) = rcm(i,k) + rcm_supersat_adj(i,k)
            l_spur_supersat = .true.
          end if
        end do
      end do
      !$acc end parallel loop

      if ( clubb_at_least_debug_level_api( 1 ) .and. l_spur_supersat ) then
        write(fstderr,*) 'Warning: spurious supersaturation was removed after pdf_closure!'
      end if

    end if ! l_rcm_supersat_adj

    !$acc exit data delete( wp2_zt,wp3_zm, rtp2_zt,rtp3_zm, thlp2_zt,  thlp3_zm, &
    !$acc                   wprtp_zt, wpthlp_zt, rtpthlp_zt, up2_zt, up3_zm, &
    !$acc                   vp2_zt, vp3_zm, upwp_zt, vpwp_zt, gamma_Skw_fnc, &
    !$acc                   gamma_Skw_fnc_zt,sigma_sqd_w_zt,  Skw_zt, Skw_zm, &
    !$acc                   Skrt_zt, Skrt_zm, Skthl_zt, Skthl_zm, Sku_zt, &
    !$acc                   Sku_zm, Skv_zt, Skv_zm, wp2up2_zt, &
    !$acc                   wp2vp2_zt, wp4_zt, wpthvp_zt, rtpthvp_zt, thlpthvp_zt, &
    !$acc                   wprcp_zt, rtprcp_zt, thlprcp_zt, uprcp_zt, vprcp_zt, &
    !$acc                   rsat, rel_humidity, um_zm, vm_zm, T_in_K, sigma_sqd_w_tmp )

    !$acc exit data if( l_call_pdf_closure_twice ) &
    !$acc           delete( w_up_in_cloud_zm, wpup2_zm, wpvp2_zm, &
    !$acc                   w_down_in_cloud_zm, cloudy_updraft_frac_zm,  &
    !$acc                   cloudy_downdraft_frac_zm, p_in_Pa_zm, exner_zm, &
    !$acc                   wprtp2_zm, wp2rtp_zm, wpthlp2_zm, &
    !$acc                   wp2thlp_zm, wprtpthlp_zm, wp2thvp_zm, wp2rcp_zm, wp2up_zm )

    !$acc exit data if( sclr_dim > 0 ) &
    !$acc           delete( wpsclrp_zt, sclrp2_zt, sclrp3_zm, sclrprtp_zt, sclrpthlp_zt, &
    !$acc                   Sksclr_zt, Sksclr_zm, sclrpthvp_zt, sclrprcp_zt, wpsclrprtp_zm, &
    !$acc                   wpsclrp2_zm, wpsclrpthlp_zm, wp2sclrp_zm, sclrm_zm )

    !$acc exit data if( hydromet_dim > 0 ) delete( wphydrometp_zt, wp2hmp_zm, rtphmp, thlphmp )

    return

  end subroutine pdf_closure_driver

  !-----------------------------------------------------------------------
  subroutine trapezoidal_rule_zt( nzm, nzt, ngrdcol, sclr_dim, gr,             & ! intent(in)
                                  l_call_pdf_closure_twice,                    & ! intent(in)
                                  stats,                                       & ! intent(in)
                                  wprtp2, wpthlp2,                             & ! intent(inout)
                                  wprtpthlp, cloud_frac, ice_supersat_frac,    & ! intent(inout)
                                  rcm, wp2thvp, wp2up, wpsclrprtp, wpsclrp2,   & ! intent(inout)
                                  wpsclrpthlp,                                 & ! intent(inout)
                                  wprtp2_zm, wpthlp2_zm,                       & ! intent(inout)
                                  wprtpthlp_zm, cloud_frac_zm,                 & ! intent(inout)
                                  ice_supersat_frac_zm, rcm_zm, wp2thvp_zm,    & ! intent(inout)
                                  wp2up_zm,                                    & ! intent(inout)
                                  wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm )   ! intent(inout)
               
    !
    ! Description:
    !   This subroutine takes the output variables on the thermo.
    !   grid and either: interpolates them to the momentum grid, or uses the
    !   values output from the second call to pdf_closure on momentum levels if
    !   l_call_pdf_closure_twice is true.  It then calls the function
    !   trapezoid_zt to recompute the variables on the thermo. grid.
    !
    !   ldgrant June 2009
    !
    ! Note:
    !   The argument variables in the last 5 lines of the subroutine
    !   (wprtp2_zm through pdf_params_zm) are declared intent(inout) because
    !   if l_call_pdf_closure_twice is true, these variables will already have
    !   values from pdf_closure on momentum levels and will not be altered in
    !   this subroutine.  However, if l_call_pdf_closure_twice is false, these
    !   variables will not have values yet and will be interpolated to
    !   momentum levels in this subroutine.
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid,     & ! Type
        zt2zm_api   ! Procedure

    use pdf_parameter_module, only: &
        pdf_parameter ! Derived data type

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_netcdf, only: &
      stats_type, &
      var_on_stats_list

    implicit none

    !------------------------ Input variables ------------------------
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol, &
      sclr_dim

    type (grid), intent(in) :: &
      gr

    logical, intent(in) :: &
      l_call_pdf_closure_twice

    type(stats_type), intent(in) :: &
      stats

    !------------------------ Input/Output variables ------------------------
    ! Thermodynamic level variables output from the first call to pdf_closure
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
      wprtp2,             & ! w'rt'^2                   [m kg^2/kg^2]
      wpthlp2,            & ! w'thl'^2                  [m K^2/s]
      wprtpthlp,          & ! w'rt'thl'                 [m kg K/kg s]
      cloud_frac,         & ! Cloud Fraction            [-]
      ice_supersat_frac,  & ! Ice Cloud Fraction        [-]
      rcm,                & ! Liquid water mixing ratio [kg/kg]
      wp2thvp,            & ! w'^2 th_v'                [m^2 K/s^2]
      wp2up                 ! w'^2 u'                   [m^3/s^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,sclr_dim), intent(inout) :: &
      wpsclrprtp,  & ! w'sclr'rt'
      wpsclrp2,    & ! w'sclr'^2
      wpsclrpthlp    ! w'sclr'thl'

    ! Thermo. level variables brought to momentum levels either by
    ! interpolation (in subroutine trapezoidal_rule_zt) or by
    ! the second call to pdf_closure (in subroutine advance_clubb_core)
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
      wprtp2_zm,            & ! w'rt'^2 on momentum grid                   [m kg^2/kg^2]
      wpthlp2_zm,           & ! w'thl'^2 on momentum grid                  [m K^2/s]
      wprtpthlp_zm,         & ! w'rt'thl' on momentum grid                 [m kg K/kg s]
      cloud_frac_zm,        & ! Cloud Fraction on momentum grid            [-]
      ice_supersat_frac_zm, & ! Ice Cloud Fraction on momentum grid        [-]
      rcm_zm,               & ! Liquid water mixing ratio on momentum grid [kg/kg]
      wp2thvp_zm,           & ! w'^2 th_v' on momentum grid                [m^2 K/s^2]
      wp2up_zm                ! w'^2 u' on momentum grid                   [m^3/s^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzm,sclr_dim), intent(inout) :: &
      wpsclrprtp_zm,  & ! w'sclr'rt' on momentum grid
      wpsclrp2_zm,    & ! w'sclr'^2 on momentum grid
      wpsclrpthlp_zm    ! w'sclr'thl' on momentum grid

    !------------------------ Local variables ------------------------

    integer :: i, sclr

    !----------------------- Begin Code -----------------------------

    ! Store components of pdf_params in the locally declared variables
    ! We only apply the trapezoidal rule to these when
    ! l_apply_rule_to_pdf_params is true.  This is because when we apply the
    ! rule to the final result of pdf_closure rather than the intermediate
    ! results it can lead to an inconsistency in how we determine which
    ! PDF component a point is in and whether the point is in or out of cloud,
    ! which is turn will break the latin hypercube code that samples
    ! preferentially in cloud. -dschanen 13 Feb 2012


    ! If l_call_pdf_closure_twice is true, the _zm variables already have
    ! values from the second call to pdf_closure in advance_clubb_core.
    ! If it is false, the variables are interpolated to the _zm levels.
    if ( .not. l_call_pdf_closure_twice ) then

      ! Interpolate thermodynamic variables to the momentum grid.
      wprtp2_zm                   = zt2zm_api( nzm, nzt, ngrdcol, gr, wprtp2 )
      wpthlp2_zm                  = zt2zm_api( nzm, nzt, ngrdcol, gr, wpthlp2 )
      wprtpthlp_zm                = zt2zm_api( nzm, nzt, ngrdcol, gr, wprtpthlp )
      cloud_frac_zm               = zt2zm_api( nzm, nzt, ngrdcol, gr, cloud_frac )
      ice_supersat_frac_zm        = zt2zm_api( nzm, nzt, ngrdcol, gr, ice_supersat_frac )
      rcm_zm                      = zt2zm_api( nzm, nzt, ngrdcol, gr, rcm )
      wp2thvp_zm                  = zt2zm_api( nzm, nzt, ngrdcol, gr, wp2thvp )
      wp2up_zm                    = zt2zm_api( nzm, nzt, ngrdcol, gr, wp2up )

      ! Since top momentum level is higher than top thermo. level,
      ! set variables at top momentum level to 0.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        wprtp2_zm(i,gr%k_ub_zm)            = 0.0_core_rknd
        wpthlp2_zm(i,gr%k_ub_zm)           = 0.0_core_rknd
        wprtpthlp_zm(i,gr%k_ub_zm)         = 0.0_core_rknd
        cloud_frac_zm(i,gr%k_ub_zm)        = 0.0_core_rknd
        ice_supersat_frac_zm(i,gr%k_ub_zm) = 0.0_core_rknd
        rcm_zm(i,gr%k_ub_zm)               = 0.0_core_rknd
        wp2thvp_zm(i,gr%k_ub_zm)           = 0.0_core_rknd
        wp2up_zm(i,gr%k_ub_zm)             = 0.0_core_rknd
      end do
      !$acc end parallel loop

      do sclr = 1, sclr_dim
        wpsclrprtp_zm(:,:,sclr)   = zt2zm_api( nzm, nzt, ngrdcol, gr, wpsclrprtp(:,:,sclr) )
        wpsclrp2_zm(:,:,sclr)     = zt2zm_api( nzm, nzt, ngrdcol, gr, wpsclrp2(:,:,sclr) )
        wpsclrpthlp_zm(:,:,sclr)  = zt2zm_api( nzm, nzt, ngrdcol, gr, wpsclrpthlp(:,:,sclr) )

        !$acc parallel loop gang vector default(present)
        do i = 1, ngrdcol
          wpsclrprtp_zm(i,gr%k_ub_zm,sclr)  = 0.0_core_rknd
          wpsclrp2_zm(i,gr%k_ub_zm,sclr)    = 0.0_core_rknd
          wpsclrpthlp_zm(i,gr%k_ub_zm,sclr) = 0.0_core_rknd
        end do
        !$acc end parallel loop
      end do ! sclr = 1, sclr_dim

    end if ! .not. l_call_pdf_closure_twice

    if ( stats%enabled ) then

      ! Use the trapezoidal rule to recompute the variables on the stats_zt level
      if ( var_on_stats_list( stats, "wprtp2" ) ) then
        call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                wprtp2_zm, &
                                wprtp2 )
      end if

      if ( var_on_stats_list( stats, "wpthlp2" ) ) then
        call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                wpthlp2_zm, &
                                wpthlp2 )
      end if

      if ( var_on_stats_list( stats, "wprtpthlp" ) ) then
        call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                wprtpthlp_zm, &
                                wprtpthlp )
      end if

      do sclr = 1, sclr_dim
        call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                wpsclrprtp_zm(:,:,sclr), &
                                wpsclrprtp(:,:,sclr) )
        call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                wpsclrpthlp_zm(:,:,sclr), &
                                wpsclrpthlp(:,:,sclr) )
        call calc_trapezoid_zt( nzm, nzt, ngrdcol,  gr, &
                                wpsclrp2_zm(:,:,sclr), &
                                wpsclrp2(:,:,sclr) )
      end do ! sclr = 1, sclr_dim
      
    end if

    call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                            cloud_frac_zm, &
                            cloud_frac )
                            
    call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                            ice_supersat_frac_zm, &
                            ice_supersat_frac )

    call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                            rcm_zm, &
                            rcm )

    call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                            wp2thvp_zm, &
                            wp2thvp )

     call calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                              wp2up_zm, &
                              wp2up )

    ! End of trapezoidal rule

    return
  end subroutine trapezoidal_rule_zt
  
  !-----------------------------------------------------------------------
  subroutine trapezoidal_rule_zm( nzm, nzt, ngrdcol, gr,              & ! intent(in)
                                  wpthvp_zt, thlpthvp_zt, rtpthvp_zt, & ! intent(in)
                                  wpthvp, thlpthvp, rtpthvp )           ! intent(inout)
    !
    ! Description:
    !   This subroutine recomputes three variables on the
    !   momentum grid from pdf_closure -- wpthvp, thlpthvp, and
    !   rtpthvp -- by calling the function trapezoid_zm.  Only these three
    !   variables are used in this subroutine because they are the only
    !   pdf_closure momentum variables used elsewhere in CLUBB.
    !
    !   The _zt variables are output from the first call to pdf_closure.
    !   The _zm variables are output from the second call to pdf_closure
    !   on the momentum levels.
    !   This is done before the call to this subroutine.
    !
    !   ldgrant Feb. 2010
    !
    !  References:
    !    None
    !-----------------------------------------------------------------------

    use grid_class, only: grid

    use clubb_precision, only: &
      core_rknd ! variable(s)

    implicit none

    ! ----------------------- Input variables -----------------------
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol

    type (grid), intent(in) :: gr
  
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      wpthvp_zt,   & ! Buoyancy flux (on thermo. grid)  [(K m)/s]
      thlpthvp_zt, & ! th_l' th_v' (on thermo. grid)    [K^2]
      rtpthvp_zt     ! r_t' th_v' (on thermo. grid)     [(kg K)/kg]

    ! ----------------------- Input/Output variables -----------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
      wpthvp,   & ! Buoyancy flux   [(K m)/s]
      thlpthvp, & ! th_l' th_v'     [K^2]
      rtpthvp     ! r_t' th_v'      [(kg K)/kg]

    ! ----------------------- Begin Code -----------------------

    ! Use the trapezoidal rule to recompute the variables on the zm level
    call calc_trapezoid_zm( nzm, nzt, ngrdcol, gr, wpthvp_zt, &         ! Intent(in) 
                            wpthvp )                                    ! Intent(inout)
                       
    call calc_trapezoid_zm( nzm, nzt, ngrdcol, gr, thlpthvp_zt, &       ! Intent(in)
                            thlpthvp )                                  ! Intent(inout)
                       
    call calc_trapezoid_zm( nzm, nzt, ngrdcol, gr, rtpthvp_zt, &        ! Intent(in)
                            rtpthvp )                                   ! Intent(inout)

    return
  end subroutine trapezoidal_rule_zm

  !-----------------------------------------------------------------------
  subroutine calc_trapezoid_zt( nzm, nzt, ngrdcol, gr, &
                                variable_zm, &
                                variable_zt )
    !
    ! Description:
    !   Function which uses the trapezoidal rule from calculus
    !   to recompute the values for the variables on the thermo. grid which
    !   are output from the first call to pdf_closure in module clubb_core.
    !
    !   ldgrant June 2009
    !--------------------------------------------------------------------

    use constants_clubb, only: &
        one_half

    use grid_class, only: &
        grid

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! ---------------- Input Variables ----------------
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol

    type (grid), intent(in) :: gr

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      variable_zm    ! Variable on the zm grid

    ! ---------------- Input/Output Variable ----------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
      variable_zt    ! Variable on the zt grid

    ! ---------------- Local Variables ----------------
    integer :: i, k, k_zm, k_zmp1 ! Loop index

    ! ---------------- Begin Code ----------------

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = gr%k_lb_zt, gr%k_ub_zt, gr%grid_dir_indx
      do i = 1, ngrdcol

        if ( gr%grid_dir_indx > 0 ) then
          ! Ascending grid
          k_zmp1 = k+1
          k_zm = k
        else ! gr%grid_dir_indx < 0
          ! Descending grid
          k_zmp1 = k
          k_zm = k+1
        endif

        ! Trapezoidal rule from calculus
        variable_zt(i,k) &
        = one_half &
          * ( variable_zm(i,k_zmp1) + variable_zt(i,k) ) &
          * ( gr%zm(i,k_zmp1) - gr%zt(i,k) ) &
          * gr%grid_dir * gr%invrs_dzt(i,k) &
          + one_half &
            * ( variable_zt(i,k) + variable_zm(i,k_zm) ) &
            * ( gr%zt(i,k) - gr%zm(i,k_zm) ) &
            * gr%grid_dir * gr%invrs_dzt(i,k)

      end do
    end do ! k = gr%k_lb_zt, gr%k_ub_zt, gr%grid_dir_indx
    !$acc end parallel loop

    return
  end subroutine calc_trapezoid_zt

  !-----------------------------------------------------------------------
  subroutine calc_trapezoid_zm( nzm, nzt, ngrdcol, gr, variable_zt, &
                                variable_zm )
    !
    ! Description:
    !   Function which uses the trapezoidal rule from calculus
    !   to recompute the values for the important variables on the momentum
    !   grid which are output from pdf_closure in module clubb_core.
    !   These momentum variables only include wpthvp, thlpthvp, and rtpthvp.
    !
    !   ldgrant Feb. 2010
    !--------------------------------------------------------------------

    use constants_clubb, only: &
        one_half

    use grid_class, only: &
        grid

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! -------------------- Input Variables --------------------
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol

    type (grid), intent(in) :: gr

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      variable_zt    ! Variable on the zt grid

    ! -------------------- Input/Output Variable --------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
      variable_zm    ! Variable on the zm grid

    ! -------------------- Local Variables --------------------
    integer :: i, k, k_zt, k_ztm1 ! Loop index

    ! -------------------- Begin Code --------------------

    ! Boundary conditions: trapezoidal rule not valid at top zm level, nzmax.
    ! Trapezoidal rule also not used at zm level 1.

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = gr%k_lb_zm+gr%grid_dir_indx, gr%k_ub_zm-gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol

        if ( gr%grid_dir_indx > 0 ) then
          ! Ascending grid
          k_zt = k
          k_ztm1 = k-1
        else ! gr%grid_dir_indx < 0
          ! Descending grid
          k_zt = k-1
          k_ztm1 = k
        endif

        ! Trapezoidal rule from calculus
        variable_zm(i,k) &
        = one_half &
          * ( variable_zt(i,k_zt) + variable_zm(i,k) ) &
          * ( gr%zt(i,k_zt) - gr%zm(i,k) ) &
          * gr%grid_dir * gr%invrs_dzm(i,k) &
          + one_half &
            * ( variable_zm(i,k) + variable_zt(i,k_ztm1) ) &
            * ( gr%zm(i,k) - gr%zt(i,k_ztm1) ) &
            * gr%grid_dir * gr%invrs_dzm(i,k)

      end do
    end do 
    !$acc end parallel loop

    return
  end subroutine calc_trapezoid_zm

  !-----------------------------------------------------------------------
  subroutine compute_cloud_cover( gr, nzt, ngrdcol, &
                                  pdf_params, cloud_frac, rcm, & ! intent(in)
                                  err_info,                    & ! intent(inout)
                                  cloud_cover, rcm_in_layer )    ! intent(out)
    !
    ! Description:
    !   Subroutine to compute cloud cover (the amount of sky
    !   covered by cloud) and rcm in layer (liquid water mixing ratio in
    !   the portion of the grid box filled by cloud).
    !
    ! References:
    !   Definition of 's' comes from:
    !   ``The Gaussian Cloud Model Relations'' G. L. Mellor (1977)
    !   JAS, Vol. 34, pp. 356--358.
    !
    ! Notes:
    !   Added July 2009
    !---------------------------------------------------------------------

    use constants_clubb, only: &
        rc_tol, & ! Variable(s)
        fstderr, &
        unused_var

    use grid_class, only: grid

    use pdf_parameter_module, only: &
        pdf_parameter ! Derived data type

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
      clubb_at_least_debug_level_api,  & ! Procedure
      clubb_fatal_error              ! Constant

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    !------------------------ Input variables ------------------------
    integer, intent(in) :: &
      ngrdcol,  & ! Number of grid columns
      nzt         ! Number of thermodynamic vertical levels

    type (grid), intent(in) :: gr

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      cloud_frac, & ! Cloud fraction             [-]
      rcm           ! Liquid water mixing ratio  [kg/kg]

    type (pdf_parameter), intent(in) :: &
      pdf_params    ! PDF Parameters  [units vary]

    !------------------------ Input/Output variables ------------------------
  type(err_info_type), intent(inout) :: &
    err_info        ! err_info struct containing err_code and err_header

    !------------------------ Output variables ------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      cloud_cover,  & ! Cloud cover                               [-]
      rcm_in_layer    ! Liquid water mixing ratio in cloud layer  [kg/kg]

    !------------------------ Local variables ------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      chi_mean,              & ! Mean extended cloud water mixing ratio of the
                               ! two Gaussian distributions
      vert_cloud_frac_upper, & ! Fraction of cloud in top half of grid box
      vert_cloud_frac_lower, & ! Fraction of cloud in bottom half of grid box
      vert_cloud_frac          ! Fraction of cloud filling the grid box in the vertical

    integer :: i, k, kp1, km1, k_zmp1, k_zm

    !------------------------ Begin code ------------------------

    !$acc enter data create( chi_mean, vert_cloud_frac_upper, &
    !$acc                    vert_cloud_frac_lower, vert_cloud_frac )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol

        chi_mean(i,k) &
        = pdf_params%mixt_frac(i,k) * pdf_params%chi_1(i,k) &
          + ( 1.0_core_rknd - pdf_params%mixt_frac(i,k) ) * pdf_params%chi_2(i,k)

      end do
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = gr%k_lb_zt, gr%k_ub_zt-gr%grid_dir_indx, gr%grid_dir_indx
      do i = 1, ngrdcol

        if ( gr%grid_dir_indx > 0 ) then
           km1 = max( k-1, 1 )
           kp1 = min( k+1, nzt )
        else ! gr%grid_dir_indx < 0
           km1 = min( k+1, nzt )
           kp1 = max( k-1, 1 )
        endif

        if ( rcm(i,k) < rc_tol ) then ! No cloud at this level

          cloud_cover(i,k)  = cloud_frac(i,k)
          rcm_in_layer(i,k) = rcm(i,k)

        else if ( ( rcm(i,kp1) >= rc_tol ) .and. ( rcm(i,km1) >= rc_tol ) ) then
          ! There is cloud above and below,
          !   so assume cloud fills grid box from top to bottom

          cloud_cover(i,k) = cloud_frac(i,k)
          rcm_in_layer(i,k) = rcm(i,k)

        else if ( ( rcm(i,kp1) < rc_tol ) .or. ( rcm(i,km1) < rc_tol) ) then
          ! Cloud may fail to reach gridbox top or base or both

          ! First let the cloud fill the entire grid box, then overwrite
          ! vert_cloud_frac_upper(k) and/or vert_cloud_frac_lower(k)
          ! for a cloud top, cloud base, or one-point cloud.
          vert_cloud_frac_upper(i,k) = 0.5_core_rknd
          vert_cloud_frac_lower(i,k) = 0.5_core_rknd

          if ( gr%grid_dir_indx > 0 ) then
            ! Ascending grid
            k_zmp1 = k+1
            k_zm = k
          else ! gr%grid_dir_indx < 0
            ! Descending grid
            k_zmp1 = k
            k_zm = k+1
          endif

          if ( rcm(i,kp1) < rc_tol ) then ! Cloud top

            vert_cloud_frac_upper(i,k) &
            = ( ( 0.5_core_rknd / ( gr%grid_dir * gr%invrs_dzm(i,k_zmp1) ) ) &
                / ( gr%zm(i,k_zmp1) - gr%zt(i,k) ) ) &
              * ( rcm(i,k) / ( rcm(i,k) + abs( chi_mean(i,kp1) ) ) )

            vert_cloud_frac_upper(i,k) = min( 0.5_core_rknd, vert_cloud_frac_upper(i,k) )

            ! Make the transition in cloudiness more gradual than using
            ! the above min statement alone.
            vert_cloud_frac_upper(i,k) = vert_cloud_frac_upper(i,k) + &
              ( ( rcm(i,kp1)/rc_tol )*( 0.5_core_rknd -vert_cloud_frac_upper(i,k) ) )

          else

            vert_cloud_frac_upper(i,k) = 0.5_core_rknd

          end if

          if ( rcm(i,km1) < rc_tol ) then ! Cloud base

            vert_cloud_frac_lower(i,k) &
            = ( ( 0.5_core_rknd / ( gr%grid_dir * gr%invrs_dzm(i,k_zm) ) ) &
                / ( gr%zt(i,k) - gr%zm(i,k_zm) ) ) &
              * ( rcm(i,k) / ( rcm(i,k) + abs( chi_mean(i,km1) ) ) )

            vert_cloud_frac_lower(i,k) = min( 0.5_core_rknd, vert_cloud_frac_lower(i,k) )

            ! Make the transition in cloudiness more gradual than using
            ! the above min statement alone.
            vert_cloud_frac_lower(i,k) = vert_cloud_frac_lower(i,k) + &
              ( ( rcm(i,km1)/rc_tol )*( 0.5_core_rknd -vert_cloud_frac_lower(i,k) ) )

          else

            vert_cloud_frac_lower(i,k) = 0.5_core_rknd

          end if

          vert_cloud_frac(i,k) = &
            vert_cloud_frac_upper(i,k) + vert_cloud_frac_lower(i,k)

          vert_cloud_frac(i,k) = &
            max( cloud_frac(i,k), min( 1.0_core_rknd, vert_cloud_frac(i,k) ) )

          cloud_cover(i,k)  = cloud_frac(i,k) / vert_cloud_frac(i,k)
          rcm_in_layer(i,k) = rcm(i,k) / vert_cloud_frac(i,k)

        else

          ! This case should not be entered
          ! This case should be literally unreachable
          ! since all possible options are covered in the above cases
          cloud_cover(i,k) = unused_var
          rcm_in_layer(i,k) = unused_var
          ! Error in column i -> set ith entry to clubb_fatal_error
          err_info%err_code(i) = clubb_fatal_error
          write(fstderr, *) err_info%err_header(i)
          write(fstderr, *) "in compute_cloud_cover"
          write(fstderr, *) "invalid rcm values"

        end if ! rcm(k) < rc_tol

      end do ! i = 1, ngrdcol
    end do ! k = 1, gr%nzt-1, 1
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      cloud_cover(i,gr%k_ub_zt) = cloud_frac(i,gr%k_ub_zt)
      rcm_in_layer(i,gr%k_ub_zt) = rcm(i,gr%k_ub_zt)
    end do
    !$acc end parallel loop

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      !$acc update host( err_info%err_code )
      if ( any(err_info%err_code == clubb_fatal_error) ) then

        !$acc update host( pdf_params%mixt_frac, pdf_params%chi_1, pdf_params%chi_2, &
        !$acc              cloud_frac, rcm )

        write(fstderr,*)  &
           "ERROR: compute_cloud_cover entered a conditional case it should not have "

        write(fstderr,*) "pdf_params%mixt_frac = ", pdf_params%mixt_frac
        write(fstderr,*) "pdf_params%chi_1 = ", pdf_params%chi_1
        write(fstderr,*) "pdf_params%chi_2 = ", pdf_params%chi_2
        write(fstderr,*) "cloud_frac = ", cloud_frac
        write(fstderr,*) "rcm = ", rcm
      end if
    end if 

    !$acc exit data delete( chi_mean, vert_cloud_frac_upper, &
    !$acc                   vert_cloud_frac_lower, vert_cloud_frac )

    return

  end subroutine compute_cloud_cover

end module pdf_closure_module
