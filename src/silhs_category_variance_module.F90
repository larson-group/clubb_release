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
             ( nz, num_samples, d_variables, X_nl_all_levs, X_mixt_comp_all_levs, &
               istat_var, microphys_stats_vars_all, microphys_stats_vars_mean, &
               lh_sample_point_weights )

  ! Description:
  !   Computes the variance of a microphysics variable in each importance
  !   category!

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      dp,        & ! Constant(s)
      core_rknd

    use constants_clubb, only: &
      zero         ! Constant

    use microphys_stats_vars_module, only: &
      microphys_stats_vars_type, &    ! Type
      microphys_get_index             ! Procedure

    use latin_hypercube_driver_module, only: &
      num_importance_categories       ! Constant

    use stats_variables, only: &
      stats_lh_zt, &                  ! Variable(s)
      isilhs_variance_category, &
      l_stats_samp

    use stats_type_utilities, only: &
      stat_update_var                 ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,              &      ! Number of height levels
      num_samples,     &      ! Number of SILHS sample points
      d_variables             ! Number of variates in X_nl

    real( kind = dp ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs           ! SILHS samples at all height levels

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs    ! Mixture component (1 or 2) of each sample point

    integer, intent(in) :: &
      istat_var               ! Statistics variable in microphys_stats_vars
                              ! being sampled

    type(microphys_stats_vars_type), dimension(num_samples), intent(in) :: &
      microphys_stats_vars_all    ! The statistics objects to sample from, for each sample point

    type(microphys_stats_vars_type), intent(in) :: &
      microphys_stats_vars_mean   ! Overall mean of each microphysics statistics

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! Weight of SILHS sample points

    ! Local Variables
    integer, dimension(num_samples) :: &
      int_sample_category     ! Category of each sample point

    real( kind = core_rknd ) :: &
      var_value,      &
      overall_mean

    real( kind = core_rknd ), dimension(nz,num_importance_categories) :: &
      category_variance

    integer :: structure_index, isample, icat, k

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! It is assumed here that the structure_index will be the same for all the
    ! microphys_stats_vars objects.
    structure_index = microphys_get_index( istat_var, microphys_stats_vars_mean )

    category_variance    =    zero

    do k=2, nz

      overall_mean = microphys_stats_vars_mean%output_values(k,structure_index)

      int_sample_category = determine_sample_categories &
                            ( num_samples, d_variables, X_nl_all_levs(k,:,:), &
                              X_mixt_comp_all_levs(k,:) )

      do isample=1, num_samples

        icat = int_sample_category(isample)
        var_value = microphys_stats_vars_all(isample)%output_values(k,structure_index)
        category_variance(k,icat) = category_variance(k,icat) + &
          lh_sample_point_weights(isample) * ( ( var_value - overall_mean ) ** 2 )

      end do ! isample=1, num_samples

      category_variance(k,:) = category_variance(k,:) / real( num_samples, kind=core_rknd )

    end do ! k=2, nz

    ! Microphysics is not run on the lowest thermodynamic grid level.
    category_variance(1,:) = zero

    if ( l_stats_samp ) then
      do icat=1, num_importance_categories
        call stat_update_var( isilhs_variance_category(icat), &
                              category_variance(:,icat), stats_lh_zt )
      end do
    end if

    return
  end subroutine silhs_category_variance_driver
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  function determine_sample_categories( num_samples, d_variables, X_nl_one_lev, &
                                        X_mixt_comp_one_lev ) &
    result( int_sample_category )

  ! Description:
  !   Determines the importance category of each sample.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      dp       ! Constant

    use constants_clubb, only: &
      zero, &  ! Constant(s)
      zero_dp

    use latin_hypercube_driver_module, only: &
      num_importance_categories, &     ! Constant
      importance_category_type,  &     ! Type
      define_importance_categories     ! Procedure

    use corr_varnce_module, only: &
      iiPDF_chi, &
      iiPDF_rr

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples, &        ! Number of SILHS sample points
      d_variables           ! Number of variates in X_nl

    real( kind = dp ), dimension(num_samples,d_variables), intent(in) :: &
      X_nl_one_lev          ! SILHS sample vector at one height level

    integer, dimension(num_samples), intent(in) :: &
      X_mixt_comp_one_lev   ! Category of each sample point

    ! Output Variable
    integer, dimension(num_samples) :: &
      int_sample_category   ! Category of each sample

    ! Local Variables
    type(importance_category_type), dimension(num_importance_categories) :: &
      importance_categories

    type(importance_category_type) :: sample_category

    integer :: isample, icategory, found_category_index

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! We want to make sure that the output from the determine_sample_categories
    ! function is consistent with the categories from the define_importance_categories
    ! function.
    importance_categories = define_importance_categories( )

    do isample=1, num_samples

      if ( X_nl_one_lev(isample,iiPDF_chi) < zero_dp ) then
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

        if ( importance_categories(icategory)%l_in_cloud .eqv. sample_category%l_in_cloud .and. &
             importance_categories(icategory)%l_in_cloud .eqv. sample_category%l_in_cloud .and. &
             importance_categories(icategory)%l_in_cloud .eqv. sample_category%l_in_cloud ) then

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
