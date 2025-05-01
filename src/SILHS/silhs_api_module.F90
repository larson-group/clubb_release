!--------------------------------------------------------------------------------------------------
! $Id$
!==================================================================================================
!
!       ######## ########### ###        ###    ###  ########            ###     ######### #########
!     ###    ###    ###     ###        ###    ### ###    ###         ### ###   ###    ###   ###
!    ###           ###     ###        ###    ### ###               ###   ###  ###    ###   ###
!   ##########    ###     ###        ########## ##########       ########### #########    ###
!         ###    ###     ###        ###    ###        ###       ###     ### ###          ###
! ###    ###    ###     ###        ###    ### ###    ###       ###     ### ###          ###
! ######## ########### ########## ###    ###  ########        ###     ### ###       #########
!
! The SILHS API serves as the doorway through which external models can interact with SILHS.
!
!               PLEASE REMEMBER, IF ANY CODE IS CHANGED IN THIS DOCUMENT,
!                   THE CHANGES MUST BE PROPOGATED TO ALL HOST MODELS.
!
!
! Cloud Layers Unified By Binormals (CLUBB) user license 
! agreement.
!
! Thank you for your interest in CLUBB. We work hard to create a
! code that implements the best software engineering practices,
! is supported to the extent allowed by our limited resources,
! and is available without cost to non-commercial users. You may
! use CLUBB if, in return, you abide by these conditions:
!
! 1. Please cite CLUBB in presentations and publications that
!  contain results obtained using CLUBB.
!
! 2. You may not use any part of CLUBB to create or modify
!  another single-column (1D) model that is not called CLUBB.
!  However, you may modify or augment CLUBB or parts of CLUBB if
!  you include "CLUBB" in the name of the resulting single-column
!  model. For example, a user at MIT might modify CLUBB and call
!  the modified version "CLUBB-MIT." Or, for example, a user of
!  the CLM land-surface model might interface CLM to CLUBB and
!  call it "CLM-CLUBB." This naming convention recognizes the
!  contributions of both sets of developers.
!
! 3. You may implement CLUBB as a parameterization in a large-
!  scale host model that has 2 or 3 spatial dimensions without 
!  including "CLUBB" in the combined model name, but please 
!  acknowledge in presentations and publications that CLUBB has 
!  been included as a parameterization.
!
! 4. You may not provide all or part of CLUBB to anyone without 
!  prior permission from Vincent Larson (vlarson@uwm.edu). If 
!  you wish to share CLUBB with your collaborators without 
!  seeking permission, please ask your collaborators to register 
!  as CLUBB users at https://carson.math.uwm.edu/larson-group/clubb_site/ and to 
!  download CLUBB from there.
!
! 5. You may not use CLUBB for commercial purposes unless you 
!  receive permission from Vincent Larson.
!
! 6. You may not re-license all or any part of CLUBB.
!
! 7. CLUBB is provided "as is" and without warranty.
!
! We hope that CLUBB will develop into a community resource. We 
! encourage users to contribute their CLUBB modifications or 
! extensions to the CLUBB development group. We will then 
! consider them for inclusion in CLUBB. Such contributions will 
! benefit all CLUBB users. We would be pleased to acknowledge 
! contributors and list their CLUBB-related papers on our "About 
! CLUBB" webpage (https://carson.math.uwm.edu/larson-group/clubb_site/about.html) for 
! those contributors who so desire.
!
! Thanks so much and best wishes for your research!
!
! The CLUBB Development Group
! (Present and past contributors to the source code include 
! Vincent Larson, Chris Golaz, David Schanen, Brian Griffin, 
! Joshua Fasching, Adam Smith, and Michael Falk).
!------------------------------------------------------------------

module silhs_api_module

#ifdef SILHS

  use parameters_silhs, only: &
    silhs_config_flags_type,   & ! Type(s)
    eight_cluster_presc_probs, & ! Variable(s)
    vert_decorr_coef, &
    importance_prob_thresh

  use latin_hypercube_driver_module, only : &
    stats_accumulate_lh_api

  ! latin_hypercube_2D_output_api - Creates and opens the SILHS 2D output files.
  use latin_hypercube_driver_module, only: &
    latin_hypercube_2D_output_api, &
    latin_hypercube_2D_close_api

  use est_kessler_microphys_module, only: &
    est_kessler_microphys_api

  use lh_microphys_var_covar_module, only: &
    lh_microphys_var_covar_driver_api

  use parameters_silhs, only: &
    set_default_silhs_config_flags_api, &
    initialize_silhs_config_flags_type_api, &
    print_silhs_config_flags_api

#endif

  implicit none

  private

#ifdef SILHS

  public  &
    generate_silhs_sample_api, & ! Procedure(s)
    stats_accumulate_lh_api, &
    est_kessler_microphys_api, &
    clip_transform_silhs_output_api, &
    lh_microphys_var_covar_driver_api, &
    set_default_silhs_config_flags_api, &
    initialize_silhs_config_flags_type_api, &
    print_silhs_config_flags_api, &
    silhs_config_flags_type, & ! Type(s)
    eight_cluster_presc_probs, & ! Variable(s)
    vert_decorr_coef, &
    importance_prob_thresh

  public  & ! to print 2D lh samples
    latin_hypercube_2D_output_api, &
    latin_hypercube_2D_close_api

  private &
    generate_silhs_sample_api_single_col, &
    generate_silhs_sample_api_multi_col, &
    clip_transform_silhs_output_api_single_col, &
    clip_transform_silhs_output_api_multi_col

  interface generate_silhs_sample_api
    module procedure generate_silhs_sample_api_single_col
    module procedure generate_silhs_sample_api_multi_col
  end interface

  interface clip_transform_silhs_output_api
    module procedure clip_transform_silhs_output_api_single_col
    module procedure clip_transform_silhs_output_api_multi_col
  end interface

contains

  !================================================================================================
  ! generate_silhs_sample - Generates sample points of moisture, temperature, et cetera.
  !================================================================================================

  subroutine generate_silhs_sample_api_single_col( &
    iter, pdf_dim, num_samples, sequence_length, nzt, & ! In
    l_calc_weights_all_levs_itime, &
    gr, pdf_params, delta_zm, Lscale, & ! In
    lh_seed, hm_metadata, & ! In
    !rho_ds_zt, & ! In
    mu1, mu2, sigma1, sigma2, & ! In
    corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! In
    precip_fracs, silhs_config_flags, & ! In
    vert_decorr_coef, & ! In
    stats_metadata, & ! In
    stats_lh_zt, stats_lh_sfc, err_info, & ! intent(inout)
    X_nl_all_levs, X_mixt_comp_all_levs, & ! Out
    lh_sample_point_weights ) ! Out

    use grid_class, only: &
        grid   ! Type(s)

    use latin_hypercube_driver_module, only : &
      generate_silhs_sample

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use hydromet_pdf_parameter_module, only: &
      precipitation_fractions      ! Type

    use parameters_silhs, only: &
      silhs_config_flags_type ! Type

    use clubb_precision, only: &
      core_rknd
      
    use mt95, only: &
      genrand_intg

    use stats_type, only: &
      stats ! Type

    use stats_variables, only: &
      stats_metadata_type

    use corr_varnce_module, only: &
      hm_metadata_type

    use err_info_type_module, only: &
      err_info_type         ! Type

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      iter,            & ! Model iteration number
      pdf_dim,         & ! Number of variables to sample
      num_samples,     & ! Number of samples per variable
      sequence_length, & ! nt_repeat/num_samples; number of timesteps before sequence repeats.
      nzt                ! Number of vertical model levels

    type(grid), intent(in) :: &
      gr    ! Grid variable type

    type(pdf_parameter), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
      delta_zm     ! Difference in moment. altitudes    [m]

    !real( kind = core_rknd ), dimension(nzt), intent(in) :: &
    !  rho_ds_zt    ! Dry, static density on thermo. levels    [kg/m^3]

    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
      Lscale       ! Turbulent Mixing Length  [m]
      
    integer( kind = genrand_intg ), intent(in) :: &
      lh_seed      ! Random number generator seed

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(num_samples,nzt,pdf_dim) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(out), dimension(num_samples,nzt) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(out), dimension(num_samples,nzt) :: &
      lh_sample_point_weights

    ! More Input Variables!
    real( kind = core_rknd ), dimension(nzt,pdf_dim,pdf_dim), intent(in) :: &
      corr_cholesky_mtx_1, & ! Correlations Cholesky matrix (1st comp.)  [-]
      corr_cholesky_mtx_2    ! Correlations Cholesky matrix (2nd comp.)  [-]

    real( kind = core_rknd ), dimension(nzt,pdf_dim), intent(in) :: &
      mu1,    & ! Means of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>)  [units vary]
      mu2,    & ! Means of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>)  [units vary]
      sigma1, & ! Stdevs of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>) [units vary]
      sigma2    ! Stdevs of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>) [units vary]

    logical, intent(in) :: &
      l_calc_weights_all_levs_itime ! determines if vertically correlated sample points are needed
      
    type(precipitation_fractions), intent(in) :: &
      precip_fracs           ! Precipitation fractions      [-]

    type(silhs_config_flags_type), intent(in) :: &
      silhs_config_flags

    real( kind = core_rknd ), intent(in) :: &
      vert_decorr_coef    ! Empirically defined de-correlation constant [-]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    ! Input/Output Variables
    type(stats), intent(inout) :: &
      stats_lh_zt, &
      stats_lh_sfc

    type(err_info_type), intent(inout) :: &
      err_info          ! err_info struct containing err_code and err_header

    ! -------------- Local Variables --------------
      
    real( kind = core_rknd ), dimension(1,nzt) :: &
      delta_zm_col    ! Difference in moment. altitudes, with column dimension 1    [m]
      
    !real( kind = core_rknd ), dimension(1,nzt) :: &
    !  rho_ds_zt_col    ! Dry, static density on thermo. levels, with column dimension 1    [kg/m^3]
      
    real( kind = core_rknd ), dimension(1,nzt) :: &
      Lscale_col       ! Turbulent Mixing Length, with column dimension 1  [m]

    ! More Input Variables!
    real( kind = core_rknd ), dimension(1,nzt,pdf_dim,pdf_dim) :: &
      corr_cholesky_mtx_1_col, & ! Correlations Cholesky matrix, with column dimension 1
      corr_cholesky_mtx_2_col    ! Correlations Cholesky matrix, with column dimension 1

    real( kind = core_rknd ), dimension(1,nzt,pdf_dim) :: &
      mu1_col,    & ! Means of the hydrometeors, with column dimension 1
      mu2_col,    & ! Means of the hydrometeors, with column dimension 1
      sigma1_col, & ! Stdevs of the hydrometeors, with column dimension 1
      sigma2_col    ! Stdevs of the hydrometeors, with column dimension 1

    ! Output Variables
    real( kind = core_rknd ), dimension(1,num_samples,nzt,pdf_dim) :: &
      X_nl_all_levs_col ! Sample that is transformed ultimately to normal-lognormal

    integer, dimension(1,num_samples,nzt) :: &
      X_mixt_comp_all_levs_col ! Which mixture component we're in

    real( kind = core_rknd ), dimension(1,num_samples,nzt) :: &
      lh_sample_point_weights_col

    type(stats), dimension(1) :: &
      stats_lh_zt_col, &
      stats_lh_sfc_col

    ! -------------- Begin Code --------------

    delta_zm_col(1,:)                 = delta_zm
    !rho_ds_zt_col(1,:)                = rho_ds_zt
    Lscale_col(1,:)                   = Lscale
    corr_cholesky_mtx_1_col(1,:,:,:)  = corr_cholesky_mtx_1
    corr_cholesky_mtx_2_col(1,:,:,:)  = corr_cholesky_mtx_2
    mu1_col(1,:,:)                    = mu1
    mu2_col(1,:,:)                    = mu2
    sigma1_col(1,:,:)                 = sigma1
    sigma2_col(1,:,:)                 = sigma2
    stats_lh_zt_col(1)                = stats_lh_zt
    stats_lh_sfc_col(1)               = stats_lh_sfc

    !$acc data copyin( precip_fracs, precip_fracs%precip_frac_1, precip_fracs%precip_frac_2, &
    !$acc              delta_zm_col, Lscale_col, corr_cholesky_mtx_1_col, &
    !$acc              corr_cholesky_mtx_2_col, mu1_col, mu2_col, sigma1_col, sigma2_col ) &
    !$acc     copyout( X_nl_all_levs, X_mixt_comp_all_levs, lh_sample_point_weights )

    call generate_silhs_sample( &
      iter, pdf_dim, num_samples, sequence_length, nzt, 1, & ! In
      l_calc_weights_all_levs_itime, & ! In
      gr, pdf_params, delta_zm_col, Lscale_col, & ! In
      lh_seed, hm_metadata, & ! In
      !rho_ds_zt_col, & ! Unused
      mu1_col, mu2_col, sigma1_col, sigma2_col, & ! In
      corr_cholesky_mtx_1_col, corr_cholesky_mtx_2_col, & ! In
      precip_fracs, silhs_config_flags, & ! In
      vert_decorr_coef, & ! In
      stats_metadata, & ! In
      stats_lh_zt_col, stats_lh_sfc_col, err_info, & ! intent(inout)
      X_nl_all_levs_col, X_mixt_comp_all_levs_col, & ! Out
      lh_sample_point_weights_col ) ! Out

    !$acc end data

    X_nl_all_levs = X_nl_all_levs_col(1,:,:,:)
    X_mixt_comp_all_levs = X_mixt_comp_all_levs_col(1,:,:)
    lh_sample_point_weights = lh_sample_point_weights_col(1,:,:)

  end subroutine generate_silhs_sample_api_single_col
!=================================================================================================!
  subroutine generate_silhs_sample_api_multi_col( &
    iter, pdf_dim, num_samples, sequence_length, nzt, ngrdcol, & ! In
    l_calc_weights_all_levs_itime, &
    gr, pdf_params, delta_zm, Lscale, & ! In
    lh_seed, hm_metadata, & ! In
    !rho_ds_zt, & ! In
    mu1, mu2, sigma1, sigma2, & ! In
    corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! In
    precip_fracs, silhs_config_flags, & ! In
    vert_decorr_coef, & ! In
    stats_metadata, & ! In
    stats_lh_zt, stats_lh_sfc, err_info, & ! intent(inout)
    X_nl_all_levs, X_mixt_comp_all_levs, & ! Out
    lh_sample_point_weights ) ! Out

    use grid_class, only: &
        grid    ! Type(s)

    use latin_hypercube_driver_module, only : generate_silhs_sample

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use hydromet_pdf_parameter_module, only: &
      precipitation_fractions      ! Type

    use parameters_silhs, only: &
      silhs_config_flags_type ! Type

    use clubb_precision, only: &
      core_rknd
      
    use mt95, only: &
      genrand_intg

    use stats_type, only: &
      stats ! Type

    use stats_variables, only: &
      stats_metadata_type

    use corr_varnce_module, only: &
      hm_metadata_type

    use err_info_type_module, only: &
      err_info_type         ! Type

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      iter,            & ! Model iteration number
      pdf_dim,         & ! Number of variables to sample
      num_samples,     & ! Number of samples per variable
      sequence_length, & ! nt_repeat/num_samples; number of timesteps before sequence repeats.
      nzt,             & ! Number of vertical model levels
      ngrdcol            ! Number of grid columns

    type(grid), intent(in) :: &
      gr    ! Grid variable type

    type(pdf_parameter), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      delta_zm     ! Difference in moment. altitudes    [m]

    !real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
    !  rho_ds_zt    ! Dry, static density on thermo. levels    [kg/m^3]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      Lscale       ! Turbulent Mixing Length  [m]
      
    integer( kind = genrand_intg ), intent(in) :: &
      lh_seed      ! Random number generator seed

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    ! InOut Variables
    type(stats), dimension(ngrdcol), intent(inout) :: &
      stats_lh_zt, &
      stats_lh_sfc

    type(err_info_type), intent(inout) :: &
      err_info          ! err_info struct containing err_code and err_header

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,num_samples,nzt,pdf_dim) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(out), dimension(ngrdcol,num_samples,nzt) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,num_samples,nzt) :: &
      lh_sample_point_weights

    ! More Input Variables!
    real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), intent(in) :: &
      corr_cholesky_mtx_1, & ! Correlations Cholesky matrix (1st comp.)  [-]
      corr_cholesky_mtx_2    ! Correlations Cholesky matrix (2nd comp.)  [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(in) :: &
      mu1,    & ! Means of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>)  [units vary]
      mu2,    & ! Means of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>)  [units vary]
      sigma1, & ! Stdevs of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>) [units vary]
      sigma2    ! Stdevs of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>) [units vary]

    logical, intent(in) :: &
      l_calc_weights_all_levs_itime ! determines if vertically correlated sample points are needed
      
    type(precipitation_fractions), intent(in) :: &
      precip_fracs           ! Precipitation fractions      [-]

    type(silhs_config_flags_type), intent(in) :: &
      silhs_config_flags

    real( kind = core_rknd ), intent(in) :: &
      vert_decorr_coef    ! Empirically defined de-correlation constant [-]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !$acc data copyin( precip_fracs, precip_fracs%precip_frac_1, precip_fracs%precip_frac_2, &
    !$acc              delta_zm, Lscale, corr_cholesky_mtx_1, &
    !$acc              corr_cholesky_mtx_2, mu1, mu2, sigma1, sigma2 ) &
    !$acc     copyout( X_nl_all_levs, X_mixt_comp_all_levs, lh_sample_point_weights )

    call generate_silhs_sample( &
      iter, pdf_dim, num_samples, sequence_length, nzt, ngrdcol, & ! In
      l_calc_weights_all_levs_itime, & ! In
      gr, pdf_params, delta_zm, Lscale, & ! In
      lh_seed, hm_metadata, & ! In
      !rho_ds_zt, &
      mu1, mu2, sigma1, sigma2, & ! In
      corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! In
      precip_fracs, silhs_config_flags, & ! In
      vert_decorr_coef, & ! In
      stats_metadata, & ! In
      stats_lh_zt, stats_lh_sfc, err_info, & ! intent(inout)
      X_nl_all_levs, X_mixt_comp_all_levs, & ! Out
      lh_sample_point_weights ) ! Out

    !$acc end data

  end subroutine generate_silhs_sample_api_multi_col

  !================================================================================================
  ! clip_transform_silhs_output - Computes extra SILHS sample variables, such as rt and thl.
  !================================================================================================

  subroutine clip_transform_silhs_output_api_single_col( &
                                              nzt, num_samples,                   & ! In
                                              pdf_dim, hydromet_dim, hm_metadata, & ! In
                                              X_mixt_comp_all_levs,               & ! In
                                              X_nl_all_levs,                      & ! Inout
                                              pdf_params, l_use_Ncn_to_Nc,        & ! In
                                              lh_rt_clipped, lh_thl_clipped,      & ! Out
                                              lh_rc_clipped, lh_rv_clipped,       & ! Out
                                              lh_Nc_clipped                       ) ! Out

    use grid_class, only: grid ! Type

    use latin_hypercube_driver_module, only : clip_transform_silhs_output

    use clubb_precision, only: &
      core_rknd       ! Our awesome generalized precision (constant)

    use pdf_parameter_module, only: &
      pdf_parameter

    use corr_varnce_module, only: &
      hm_metadata_type

    implicit none

    ! Input Variables
    logical, intent(in) :: l_use_Ncn_to_Nc

    integer, intent(in) :: &
      nzt,           & ! Number of vertical levels
      num_samples,  & ! Number of SILHS sample points
      pdf_dim,      & ! Number of variates in X_nl_one_lev
      hydromet_dim    ! Number of hydrometeor species

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    integer, dimension(num_samples,nzt), intent(in) :: &
      X_mixt_comp_all_levs   ! Which component this sample is in (1 or 2)

    real( kind = core_rknd ), dimension(num_samples,nzt,pdf_dim), intent(inout) :: &
      X_nl_all_levs         ! SILHS sample points    [units vary]

    type(pdf_parameter), intent(in) :: &
      pdf_params             ! **The** PDF parameters!

    ! Output Variables
    real( kind = core_rknd ), dimension(num_samples,nzt), intent(out) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points
      
    ! -------------- Local Variables --------------
    
    integer, dimension(1,num_samples,nzt) :: &
      X_mixt_comp_all_levs_col   ! Which component this sample is in (1 or 2)

    real( kind = core_rknd ), dimension(1,num_samples,nzt,pdf_dim) :: &
      X_nl_all_levs_col         ! SILHS sample points    [units vary]

    ! Output Variables
    real( kind = core_rknd ), dimension(1,num_samples,nzt) :: &
      lh_rt_clipped_col,  & ! rt generated from silhs sample points
      lh_thl_clipped_col, & ! thl generated from silhs sample points
      lh_rc_clipped_col,  & ! rc generated from silhs sample points
      lh_rv_clipped_col,  & ! rv generated from silhs sample points
      lh_Nc_clipped_col     ! Nc generated from silhs sample points
    
    ! -------------- Begin Code --------------
    
    X_mixt_comp_all_levs_col(1,:,:) = X_mixt_comp_all_levs
    X_nl_all_levs_col(1,:,:,:)      = X_nl_all_levs

    !$acc data copyin( pdf_params, pdf_params%rt_1, pdf_params%thl_1, & 
    !$acc              pdf_params%rt_2, pdf_params%thl_2, pdf_params%crt_1, pdf_params%cthl_1, & 
    !$acc              pdf_params%crt_2, pdf_params%cthl_2, pdf_params%chi_1, pdf_params%chi_2, &
    !$acc              X_mixt_comp_all_levs_col, X_nl_all_levs_col, hm_metadata ) &
    !$acc     copyout( lh_rt_clipped_col, lh_thl_clipped_col,  lh_rc_clipped_col, &
    !$acc              lh_rv_clipped_col, lh_Nc_clipped_col )

    call clip_transform_silhs_output( nzt, 1, num_samples,                 & ! In
                                      pdf_dim, hydromet_dim, hm_metadata,  & ! In
                                      X_mixt_comp_all_levs_col,               & ! In
                                      X_nl_all_levs_col,                      & ! In
                                      pdf_params, l_use_Ncn_to_Nc,            & ! In
                                      lh_rt_clipped_col, lh_thl_clipped_col,  & ! Out
                                      lh_rc_clipped_col, lh_rv_clipped_col,   & ! Out
                                      lh_Nc_clipped_col                       ) ! Out
    !$acc end data
                                      
    lh_rt_clipped     = lh_rt_clipped_col(1,:,:)
    lh_thl_clipped    = lh_thl_clipped_col(1,:,:)
    lh_rc_clipped     = lh_rc_clipped_col(1,:,:)
    lh_rv_clipped     = lh_rv_clipped_col(1,:,:)
    lh_Nc_clipped     = lh_Nc_clipped_col(1,:,:)
                                      
  end subroutine clip_transform_silhs_output_api_single_col
!=======================================================================================!
  subroutine clip_transform_silhs_output_api_multi_col( &
                                              nzt, ngrdcol, num_samples,              & ! In
                                              pdf_dim, hydromet_dim, hm_metadata,     & ! In
                                              X_mixt_comp_all_levs,                   & ! In
                                              X_nl_all_levs,                          & ! Inout
                                              pdf_params, l_use_Ncn_to_Nc,            & ! In
                                              lh_rt_clipped, lh_thl_clipped,          & ! Out
                                              lh_rc_clipped, lh_rv_clipped,           & ! Out
                                              lh_Nc_clipped                           ) ! Out

    use grid_class, only: grid ! Type

    use latin_hypercube_driver_module, only : clip_transform_silhs_output

    use clubb_precision, only: &
      core_rknd       ! Our awesome generalized precision (constant)

    use pdf_parameter_module, only: &
      pdf_parameter

    use corr_varnce_module, only: &
      hm_metadata_type

    implicit none

    ! Input Variables
    logical, intent(in) :: l_use_Ncn_to_Nc

    integer, intent(in) :: &
      nzt,           & ! Number of vertical levels
      ngrdcol,      & ! Number of grid columns
      num_samples,  & ! Number of SILHS sample points
      pdf_dim,      & ! Number of variates in X_nl_one_lev
      hydromet_dim    ! Number of hydrometeor species

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    integer, dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      X_mixt_comp_all_levs   ! Which component this sample is in (1 or 2)

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt,pdf_dim), intent(inout) :: &
      X_nl_all_levs         ! SILHS sample points    [units vary]

    type(pdf_parameter), intent(in) :: &
      pdf_params             ! **The** PDF parameters!

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt), intent(out) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    !$acc data copyin( pdf_params, pdf_params%rt_1, pdf_params%thl_1, & 
    !$acc              pdf_params%rt_2, pdf_params%thl_2, pdf_params%crt_1, pdf_params%cthl_1, & 
    !$acc              pdf_params%crt_2, pdf_params%cthl_2, pdf_params%chi_1, pdf_params%chi_2, &
    !$acc              X_mixt_comp_all_levs, X_nl_all_levs, hm_metadata ) &
    !$acc     copyout( lh_rt_clipped, lh_thl_clipped,  lh_rc_clipped, &
    !$acc              lh_rv_clipped, lh_Nc_clipped )

    call clip_transform_silhs_output( nzt, ngrdcol, num_samples,           & ! In
                                      pdf_dim, hydromet_dim, hm_metadata,  & ! In
                                      X_mixt_comp_all_levs,                   & ! In
                                      X_nl_all_levs,                          & ! In
                                      pdf_params, l_use_Ncn_to_Nc,            & ! In
                                      lh_rt_clipped, lh_thl_clipped,          & ! Out
                                      lh_rc_clipped, lh_rv_clipped,           & ! Out
                                      lh_Nc_clipped                           ) ! Out
   !$acc end data

  end subroutine clip_transform_silhs_output_api_multi_col

#endif

end module silhs_api_module
