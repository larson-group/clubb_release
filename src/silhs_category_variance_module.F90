!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module silhs_category_variance_module

  implicit none

  public :: silhs_category_variance_driver

  private ! Set Default Scope

  contains

  !-----------------------------------------------------------------------
  subroutine silhs_category_variance_driver &
             ( nz, num_samples, d_variables, hydromet_dim, X_nl_all_levs, &
               X_mixt_comp_all_levs, microphys_stats_vars_all, &
               lh_hydromet_mc_all, lh_sample_point_weights, pdf_params, &
               hydromet_pdf_params )

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
      microphys_stats_vars_type, &    ! Type
      microphys_get_index             ! Procedure

    use stats_variables, only: &
      irrm_auto

    use array_index, only: &
      iirrm

    use pdf_parameter_module, only: &
      pdf_parameter     ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter

    use corr_varnce_module, only: &
      iiPDF_chi

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,              &      ! Number of height levels
      num_samples,     &      ! Number of SILHS sample points
      d_variables,     &      ! Number of variates in X_nl
      hydromet_dim            ! Number of elements of hydromet array

    real( kind = core_rknd ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs           ! SILHS samples at all height levels

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs    ! Mixture component (1 or 2) of each sample point

    type(microphys_stats_vars_type), dimension(num_samples), intent(in) :: &
      microphys_stats_vars_all! The statistics objects to sample from, for each sample point

    real( kind = core_rknd ), dimension(nz,num_samples,hydromet_dim), intent(in) :: &
      lh_hydromet_mc_all      ! Tendencies of hydometeors at all sample points

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! Weight of SILHS sample points

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params              ! The PDF parameters!

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params

    ! Local Variables
    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      samples_all

    integer :: structure_index, isample

    integer :: istat_var  ! Statistics variable in microphys_stats_vars being sampled

  !-----------------------------------------------------------------------

    !----- Begin Code -----


    if ( .false. ) then

      ! Sample rrm_auto?
      istat_var = irrm_auto

      ! It is assumed here that the structure_index will be the same for all the
      ! microphys_stats_vars objects.
      structure_index = microphys_get_index( istat_var, microphys_stats_vars_all(1) )

      do isample=1, num_samples
        samples_all(:,isample) = microphys_stats_vars_all(isample)%output_values &
                                   (:,structure_index)
      end do

    else if ( .false. ) then

      samples_all = X_nl_all_levs(:,:,iiPDF_chi)

    else ! .true.

      samples_all = lh_hydromet_mc_all(:,:,iirrm) ! Sample rrm_mc

    end if ! .false.

    call silhs_sample_category_variance &
         ( nz, num_samples, d_variables, X_nl_all_levs, X_mixt_comp_all_levs, &
           samples_all, lh_sample_point_weights, pdf_params, hydromet_pdf_params )

    return
  end subroutine silhs_category_variance_driver
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine silhs_sample_category_variance &
             ( nz, num_samples, d_variables, X_nl_all_levs, X_mixt_comp_all_levs, &
               samples_all, lh_sample_point_weights, pdf_params, hydromet_pdf_params )

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
      importance_category_type

    use stats_variables, only: &
      stats_lh_zt, &                  ! Variable(s)
      isilhs_variance_category, &
      l_stats_samp

    use stats_type_utilities, only: &
      stat_update_var                 ! Procedure

    use pdf_parameter_module, only: &
      pdf_parameter       ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,              &      ! Number of height levels
      num_samples,     &      ! Number of SILHS sample points
      d_variables             ! Number of variates in X_nl

    real( kind = core_rknd ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs           ! SILHS samples at all height levels

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs    ! Mixture component (1 or 2) of each sample point

    real( kind = core_rknd ), dimension(nz,num_samples), intent(in) :: &
      samples_all             ! Sample points of variable to compute variance of

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! Weight of SILHS sample points

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params              ! The PDF parameters!

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params

    ! Local Variables
    type(importance_category_type), dimension(num_importance_categories) :: &
      importance_categories

    integer, dimension(num_samples) :: &
      int_sample_category     ! Category of each sample point

    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_real_probs     ! PDF probability of each category

    real( kind = core_rknd ), dimension(nz,num_importance_categories) :: &
      root_weight_mean_sq_cat

    integer :: isample, icat, k

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    root_weight_mean_sq_cat = zero

    ! We want to make sure that the output from the determine_sample_categories
    ! function is consistent with the categories from the
    ! define_importance_categories function.
    importance_categories   = define_importance_categories( )


    do k=2, nz

      int_sample_category = determine_sample_categories &
                            ( num_samples, d_variables, X_nl_all_levs(k,:,:), &
                              X_mixt_comp_all_levs(k,:), importance_categories )

      category_real_probs = &
        compute_category_real_probs( importance_categories, pdf_params(k), &
                                     hydromet_pdf_params(k) )

      do isample=1, num_samples

        icat = int_sample_category(isample)

        root_weight_mean_sq_cat(k,icat) = root_weight_mean_sq_cat(k,icat) + &
          lh_sample_point_weights(isample) * ( samples_all(k,isample) ** 2 )

      end do ! isample=1, num_samples

      root_weight_mean_sq_cat(k,:) = root_weight_mean_sq_cat(k,:) / &
                                     real( num_samples, kind=core_rknd )

      root_weight_mean_sq_cat(k,:) = sqrt( root_weight_mean_sq_cat(k,:) * &
                                           category_real_probs(:) )

    end do ! k=2, nz

    ! Microphysics is not run on the lowest thermodynamic grid level.
    root_weight_mean_sq_cat(1,:) = zero

    if ( l_stats_samp ) then
      do icat=1, num_importance_categories
        call stat_update_var( isilhs_variance_category(icat), &
                              root_weight_mean_sq_cat(:,icat), stats_lh_zt )
      end do
    end if

    return
  end subroutine silhs_sample_category_variance
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  function determine_sample_categories( num_samples, d_variables, X_nl_one_lev, &
                                        X_mixt_comp_one_lev, importance_categories ) &
  result( int_sample_category )

  ! Description:
  !   Determines the importance category of each sample.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd       ! Constant

    use constants_clubb, only: &
      zero     ! Constant(s)

    use silhs_importance_sample_module, only: &
      num_importance_categories, &     ! Constant
      importance_category_type         ! Type

    use corr_varnce_module, only: &
      iiPDF_chi, &
      iiPDF_rr

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples, &        ! Number of SILHS sample points
      d_variables           ! Number of variates in X_nl

    real( kind = core_rknd ), dimension(num_samples,d_variables), intent(in) :: &
      X_nl_one_lev          ! SILHS sample vector at one height level

    integer, dimension(num_samples), intent(in) :: &
      X_mixt_comp_one_lev   ! Category of each sample point

    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories

    ! Output Variable
    integer, dimension(num_samples) :: &
      int_sample_category   ! Category of each sample

    type(importance_category_type) :: sample_category

    integer :: isample, icategory, found_category_index

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    do isample=1, num_samples

      if ( X_nl_one_lev(isample,iiPDF_chi) < zero ) then
        sample_category%l_in_cloud = .false.
      else
        sample_category%l_in_cloud = .true.
      end if

      if ( iiPDF_rr == -1 ) then
        stop "iiPDF_rr must be greater than zero for the category sampler to work."
      end if

      if ( X_nl_one_lev(isample,iiPDF_rr) > zero ) then
        sample_category%l_in_precip = .true.
      else
        sample_category%l_in_precip = .false.
      end if

      if ( X_mixt_comp_one_lev(isample) == 1 ) then
        sample_category%l_in_component_1 = .true.
      else
        sample_category%l_in_component_1 = .false.
      end if

      found_category_index = -1

      do icategory=1, num_importance_categories

        if ( (importance_categories(icategory)%l_in_cloud .eqv. sample_category%l_in_cloud) .and. &
             (importance_categories(icategory)%l_in_precip .eqv. sample_category%l_in_precip) .and.&
             (importance_categories(icategory)%l_in_component_1 .eqv. &
              sample_category%l_in_component_1) ) then

          found_category_index = icategory
          exit

        end if

      end do ! icategory=1, num_importance_categories

      if ( found_category_index == -1 ) then
        stop "Fatal error determining category in determine_sample_categories"
      end if

      int_sample_category(isample) = found_category_index

    end do ! isample=1, num_samples

    return
  end function determine_sample_categories
  !-----------------------------------------------------------------------

end module silhs_category_variance_module
