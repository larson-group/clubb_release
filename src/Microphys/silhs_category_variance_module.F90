!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module silhs_category_variance_module

  implicit none

  public :: silhs_category_variance_driver

  private ! Set Default Scope

  contains

  !-----------------------------------------------------------------------
  subroutine silhs_category_variance_driver( &
               nzt, num_samples, pdf_dim, hydromet_dim, hm_metadata, &
               X_nl_all_levs, &
               X_mixt_comp_all_levs, microphys_stats_vars_all, &
               lh_hydromet_mc_all, lh_sample_point_weights, pdf_params, &
               precip_fracs, &
               stats, icol )

  ! Description:
  !   Computes the variance of a microphysics variable in each importance
  !   category!

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd

    use microphys_stats_vars_module, only: &
      microphys_stats_vars_type       ! Type

    use corr_varnce_module, only: &
      hm_metadata_type

    use pdf_parameter_module, only: &
      pdf_parameter     ! Type

    use hydromet_pdf_parameter_module, only: &
      precipitation_fractions

    use stats_netcdf, only: &
      stats_type

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nzt,             &      ! Number of height levels
      num_samples,     &      ! Number of SILHS sample points
      pdf_dim,         &      ! Number of variates in X_nl
      hydromet_dim            ! Number of elements of hydromet array

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), dimension(num_samples,nzt,pdf_dim), intent(in) :: &
      X_nl_all_levs           ! SILHS samples at all height levels

    integer, dimension(nzt,num_samples), intent(in) :: &
      X_mixt_comp_all_levs    ! Mixture component (1 or 2) of each sample point

    type(microphys_stats_vars_type), dimension(num_samples), intent(in) :: &
      microphys_stats_vars_all! The statistics objects to sample from, for each sample point

    real( kind = core_rknd ), dimension(num_samples,nzt,hydromet_dim), intent(in) :: &
      lh_hydromet_mc_all      ! Tendencies of hydometeors at all sample points

    real( kind = core_rknd ), dimension(num_samples,nzt), intent(in) :: &
      lh_sample_point_weights ! Weight of SILHS sample points

    type(pdf_parameter), intent(in) :: &
      pdf_params              ! The PDF parameters!

    type(precipitation_fractions), intent(in) :: &
      precip_fracs           ! Precipitation fractions      [-]

    !--------------------------- InOut Variables ---------------------------
    type(stats_type), intent(inout) :: &
      stats

    integer, intent(in) :: &
      icol

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ), dimension(num_samples,nzt) :: &
      samples_all

    !--------------------------- Begin Code ---------------------------

    ! Sample rrm_mc from the hydrometeor tendencies generated for each SILHS point.
    samples_all = lh_hydromet_mc_all(:,:,hm_metadata%iirr)

    call silhs_sample_category_variance( &
           nzt, num_samples, pdf_dim, X_nl_all_levs, X_mixt_comp_all_levs, &
           samples_all, lh_sample_point_weights, pdf_params, precip_fracs, &
           hm_metadata, stats, icol )

    return
  end subroutine silhs_category_variance_driver
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine silhs_sample_category_variance( &
               nzt, num_samples, pdf_dim, X_nl_all_levs, X_mixt_comp_all_levs, &
               samples_all, lh_sample_point_weights, pdf_params, precip_fracs, &
               hm_metadata, stats, icol )

  ! Description:
  !   Computes the variance of a microphysics variable in each importance
  !   category!

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      zero         ! Constant

    use silhs_importance_sample_module, only: &
      num_importance_categories, &       ! Constant
      define_importance_categories, &
      compute_category_real_probs, &
      importance_category_type, &
      determine_sample_categories

    use pdf_parameter_module, only: &
      pdf_parameter       ! Type

    use hydromet_pdf_parameter_module, only: &
      precipitation_fractions

    use stats_netcdf, only: &
      stats_type, &
      stats_update

    use corr_varnce_module, only: &
      hm_metadata_type

    implicit none

    !---------------------- Input Variables ----------------------
    integer, intent(in) :: &
      nzt,             & ! Number of height levels
      num_samples,     & ! Number of SILHS sample points
      pdf_dim            ! Number of variates in X_nl

    real( kind = core_rknd ), dimension(num_samples,nzt,pdf_dim), intent(in) :: &
      X_nl_all_levs           ! SILHS samples at all height levels

    integer, dimension(num_samples,nzt), intent(in) :: &
      X_mixt_comp_all_levs    ! Mixture component (1 or 2) of each sample point

    real( kind = core_rknd ), dimension(num_samples,nzt), intent(in) :: &
      samples_all             ! Sample points of variable to compute variance of

    real( kind = core_rknd ), dimension(num_samples,nzt), intent(in) :: &
      lh_sample_point_weights ! Weight of SILHS sample points

    type(pdf_parameter), intent(in) :: &
      pdf_params              ! The PDF parameters!

    type(precipitation_fractions), intent(in) :: &
      precip_fracs           ! Precipitation fractions      [-]

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    type(stats_type), intent(inout) :: &
      stats

    integer, intent(in) :: &
      icol

    !---------------------- Local Variables ----------------------
    type(importance_category_type), dimension(num_importance_categories) :: &
      importance_categories

    integer, dimension(num_samples) :: &
      int_sample_category     ! Category of each sample point

    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_real_probs     ! PDF probability of each category

    real( kind = core_rknd ), dimension(num_importance_categories,nzt) :: &
      root_weight_mean_sq_cat

    integer :: isample, icat, k
    character(len=32) :: cat_name
    character(len=2) :: cat_num

    !---------------------- Begin Code----------------------

    root_weight_mean_sq_cat = zero

    ! We want to make sure that the output from the determine_sample_categories
    ! function is consistent with the categories from the
    ! define_importance_categories function.
    importance_categories   = define_importance_categories( )

    do k = 1, nzt

      int_sample_category = determine_sample_categories( &
                              num_samples, pdf_dim, hm_metadata, &
                              X_nl_all_levs(:,k,:), &
                              X_mixt_comp_all_levs(:,k), importance_categories )

      category_real_probs = &
        compute_category_real_probs( importance_categories, &
                                     pdf_params%cloud_frac_1(1,k), pdf_params%cloud_frac_2(1,k), &
                                     pdf_params%mixt_frac(1,k), &
                                     precip_fracs%precip_frac_1(1,k), &
                                     precip_fracs%precip_frac_2(1,k) )

      do isample=1, num_samples

        icat = int_sample_category(isample)

        root_weight_mean_sq_cat(icat,k) = root_weight_mean_sq_cat(icat,k) + &
          lh_sample_point_weights(isample,k) * ( samples_all(isample,k) ** 2 )

      end do ! isample=1, num_samples

      root_weight_mean_sq_cat(:,k) = root_weight_mean_sq_cat(:,k) / &
                                     real( num_samples, kind=core_rknd )

      where ( category_real_probs > zero )
        root_weight_mean_sq_cat(:,k) = sqrt( root_weight_mean_sq_cat(:,k) / &
                                             category_real_probs(:) )
      else where
        root_weight_mean_sq_cat(:,k) = -999._core_rknd
      end where

    end do ! k = 1, nzt

    if ( stats%l_sample ) then
      do icat=1, num_importance_categories
        write(cat_num,'(I1)') icat
        cat_name = "silhs_var_cat_" // trim(cat_num)
        call stats_update( cat_name, root_weight_mean_sq_cat(icat,:), stats, icol )
      end do
    end if

    return
  end subroutine silhs_sample_category_variance
  !-----------------------------------------------------------------------

end module silhs_category_variance_module
