!-------------------------------------------------------------------------------
! $Id$
!===============================================================================

module silhs_importance_sample_module

  implicit none

  integer, parameter, public :: &
    num_importance_categories = 8 ! Number of importance sampling categories
                                  ! ( e.g., (cloud,precip,comp1) )

  private ! Default scope

  public :: importance_category_type, importance_sampling_driver, define_importance_categories, &
            compute_category_real_probs, generate_strat_uniform_variate, pick_sample_categories, &
            cloud_weighted_sampling_driver

  type importance_category_type

    logical :: &
      l_in_cloud,  &
      l_in_precip, &
      l_in_component_1

  end type importance_category_type

  contains

!-----------------------------------------------------------------------
  subroutine importance_sampling_driver &
             ( num_samples, pdf_params, hydromet_pdf_params,        &
               X_u_chi_one_lev, X_u_dp1_one_lev, X_u_dp2_one_lev,   &
               lh_sample_point_weights )

  ! Description:
  !   Applies importance sampling to a single vertical level !

  ! References:
  !   clubb:ticket:736 !
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd, &      ! Constant(s)
      dp

    use constants_clubb, only: &
      fstderr           ! Constant

    use pdf_parameter_module, only: &
      pdf_parameter     ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter ! Type

    use error_code, only: &
      clubb_at_least_debug_level  ! Procedure

    implicit none

    ! Local Parameters
    logical, parameter :: &
      l_use_prescribed_probs   = .true., &   ! Use prescribed probability importance sampling
      l_use_clustered_sampling = .false.     ! Use clustered category importance sampling

    ! Cluster allocation strategies!!!
    integer, parameter :: &
      ! All eight categories, effectively no clustering
      eight_cluster_allocation_opt = 1, &
      ! Four clusters for the combinations of cloud/no cloud and component 1/2.
      ! Precipitation fraction is ignored.
      four_cluster_allocation_opt  = 2

    integer, parameter :: &
      cluster_allocation_strategy = 1

    ! Input Variables
    integer, intent(in) :: &
      num_samples       ! Number of SILHS sample points

    type(pdf_parameter), intent(in) :: &
      pdf_params

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params

    ! Input/Output Variables
    real( kind = dp ), dimension(num_samples), intent(inout) :: &
      X_u_chi_one_lev,  & ! These uniform variates are scaled by the importance
      X_u_dp1_one_lev,  & ! sampling process.
      X_u_dp2_one_lev

    ! Output Variables
    real( kind = core_rknd ), dimension(num_samples), intent(out) :: &
      lh_sample_point_weights

    ! Local Variables
    type(importance_category_type), dimension(num_importance_categories) :: &
      importance_categories ! A vector containing the different importance categories

    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_real_probs,       & ! The real PDF probabilities for each category
      category_prescribed_probs, & ! Prescribed probability for each category
      category_sample_weights      ! Sample weight for each category

    real( kind = core_rknd ), dimension(num_samples) :: &
      rand_vect

    integer, dimension(num_samples) :: &
      int_sample_category  ! An integer for each sample corresponding to the
                           ! category picked for the sample

    integer :: sample
    logical :: l_error

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    l_error = .false.

    importance_categories = define_importance_categories( )

    category_real_probs = compute_category_real_probs &
                          ( importance_categories, pdf_params, hydromet_pdf_params )

    if ( l_use_prescribed_probs ) then

      ! Prescribe the probabilities
      category_prescribed_probs = prescribe_importance_probs &
                                  ( importance_categories, category_real_probs )

    else if ( l_use_clustered_sampling ) then

      ! A cluster allocation strategy is selected based on the value of the
      ! integer parameter.
      select case ( cluster_allocation_strategy )
      case ( eight_cluster_allocation_opt )
        category_prescribed_probs = eight_cluster_allocation &
                                    ( importance_categories, category_real_probs )
      case ( four_cluster_allocation_opt )
        stop "Not implemented"
      case default
        write(fstderr,*) "Unsupported allocation strategy:", cluster_allocation_strategy
        stop "Fatal error in importance_sampling_driver"
      end select

    else

      ! Importance sampling strategy that places half of all sample points in cloud!
      category_prescribed_probs = cloud_importance_sampling &
                                  ( importance_categories, category_real_probs, &
                                    pdf_params )

    end if ! l_use_prescribed_probs

    ! Compute weight of each sample category
    category_sample_weights = compute_category_sample_weights &
                              ( category_real_probs, category_prescribed_probs )

    ! Generate a stratified sample to be used to pick categories for the samples
    rand_vect = generate_strat_uniform_variate( num_samples )

    ! Pick the sample points!
    int_sample_category = pick_sample_categories( num_samples, category_prescribed_probs, &
                                                  rand_vect )

    ! Loop over sample points
    do sample=1, num_samples

      ! Scale and translate point to reside in its respective cateogry
      call scale_sample_to_category &
           ( importance_categories(int_sample_category(sample)), & ! In
             pdf_params, hydromet_pdf_params, & ! In
             X_u_chi_one_lev(sample), X_u_dp1_one_lev(sample), X_u_dp2_one_lev(sample) ) ! In/Out

      ! Pick a weight for the sample point
      lh_sample_point_weights(sample) = &
        category_sample_weights(int_sample_category(sample))

    end do ! sample=1, num_samples

    if ( clubb_at_least_debug_level( 2 ) ) then

      call importance_sampling_assertions &
           ( num_samples, importance_categories, category_real_probs, & ! In
             category_prescribed_probs, category_sample_weights, X_u_chi_one_lev, & ! In
             X_u_dp1_one_lev, X_u_dp2_one_lev, int_sample_category, & ! In
             pdf_params, hydromet_pdf_params, & ! In
             l_error ) ! Out

      if ( l_error ) then
        stop "Fatal error in importance_sampling_driver"
      end if ! l_error

    end if

    return
  end subroutine importance_sampling_driver
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function define_importance_categories( ) &

  result( importance_categories )

  ! Description:
  !   Creates a vector of size num_importance_categories that defines the eight
  !   importance sampling categories.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    implicit none

    ! Output Variable
    type(importance_category_type), dimension(num_importance_categories) :: &
      importance_categories ! A vector containing the different importance categories

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    importance_categories(1)%l_in_cloud       = .true.
    importance_categories(1)%l_in_precip      = .true.
    importance_categories(1)%l_in_component_1 = .true.

    importance_categories(2)%l_in_cloud       = .true.
    importance_categories(2)%l_in_precip      = .true.
    importance_categories(2)%l_in_component_1 = .false.

    importance_categories(3)%l_in_cloud       = .false.
    importance_categories(3)%l_in_precip      = .true.
    importance_categories(3)%l_in_component_1 = .true.

    importance_categories(4)%l_in_cloud       = .false.
    importance_categories(4)%l_in_precip      = .true.
    importance_categories(4)%l_in_component_1 = .false.

    importance_categories(5)%l_in_cloud       = .true.
    importance_categories(5)%l_in_precip      = .false.
    importance_categories(5)%l_in_component_1 = .true.

    importance_categories(6)%l_in_cloud       = .true.
    importance_categories(6)%l_in_precip      = .false.
    importance_categories(6)%l_in_component_1 = .false.

    importance_categories(7)%l_in_cloud       = .false.
    importance_categories(7)%l_in_precip      = .false.
    importance_categories(7)%l_in_component_1 = .true.

    importance_categories(8)%l_in_cloud       = .false.
    importance_categories(8)%l_in_precip      = .false.
    importance_categories(8)%l_in_component_1 = .false.

    return
  end function define_importance_categories
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function compute_category_real_probs( importance_categories, pdf_params, &
                                        hydromet_pdf_params ) &

  result( category_real_probs )

  ! Description:
  !   Computes the real PDF probability associated with each importance sampling
  !   category.
  !
  !   For example, if a category is in cloud, out of precipitation, and in mixture
  !   component two, then without importance sampling, the probability that a
  !   point will appear in that category is:
  !   P(cloud,noprecip,comp2) = cloud_frac_2 * (1-precip_frac_2) * (1-mixt_frac)

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      one    ! Constant

    use pdf_parameter_module, only: &
      pdf_parameter           ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter  ! Type

    implicit none

    ! Input Variables
    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories  ! A list of importance categories

    type(pdf_parameter), intent(in) :: &
      pdf_params             ! The PDF parameters!

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params    ! The hydrometeor PDF parameters!

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_real_probs ! The real PDF probabilities for each category

    ! Local Variables
    real( kind = core_rknd ) :: &
      cloud_frac_i,        &
      precip_frac_i,       &
      cloud_factor,        &
      precip_factor,       &
      component_factor

    integer :: icategory

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    do icategory = 1, num_importance_categories

      ! Determine component of category
      if ( importance_categories(icategory)%l_in_component_1 ) then
        cloud_frac_i     = pdf_params%cloud_frac_1
        precip_frac_i    = hydromet_pdf_params%precip_frac_1
        component_factor = pdf_params%mixt_frac
      else
        cloud_frac_i     = pdf_params%cloud_frac_2
        precip_frac_i    = hydromet_pdf_params%precip_frac_2
        component_factor = (one-pdf_params%mixt_frac)
      end if

      ! Determine cloud factor
      if ( importance_categories(icategory)%l_in_cloud ) then
        cloud_factor = cloud_frac_i
      else
        cloud_factor = (one-cloud_frac_i)
      end if

      ! Determine precip factor
      if ( importance_categories(icategory)%l_in_precip ) then
        precip_factor = precip_frac_i
      else
        precip_factor = (one-precip_frac_i)
      end if

      ! Compute the category probability
      category_real_probs(icategory) = component_factor * cloud_factor * precip_factor

    end do ! icategory = 1, num_importance_categories

    return
  end function compute_category_real_probs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function compute_category_sample_weights( category_real_probs, category_prescribed_probs ) &

  result( category_sample_weights )

  ! Description:
  !   Compute the sample point weights for a sample point in each category based
  !   on the PDF probability and the modified probability from importance
  !   sampling

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd     ! Constant

    use constants_clubb, only: &
      zero, &       ! Constant
      unused_var

    implicit none

    ! Input Variables

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs,  &   ! The actual PDF probability of each category
      category_prescribed_probs ! The modified probability of each category due to
                                ! importance sampling

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_sample_weights   ! Sample weight for each category

    ! Local Variable
    integer :: icategory

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    do icategory=1, num_importance_categories

      if ( category_prescribed_probs(icategory) == zero ) then
        ! If a category has no probability of being sampled, then its weight is irrevelant.
        category_sample_weights(icategory) = unused_var
      else
        category_sample_weights(icategory) = &
          category_real_probs(icategory) / category_prescribed_probs(icategory)
      end if

    end do

    return
  end function compute_category_sample_weights
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function pick_sample_categories( num_samples, category_prescribed_probs, rand_vect ) &

  result( int_sample_category )

  ! Description:
  !   Picks a category for each sample point, based on the given probabilities,
  !   such that the distribution of categories of the sample points
  !   approximates the probabilities that are given.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      dp,           &    ! Konstant(s)
      core_rknd

    use constants_clubb, only: &
      fstderr, &
      one, &
      zero

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples          ! Number of sample points to be picked

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_prescribed_probs ! Prescribed probability for each category

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      rand_vect            ! A sample of num_samples values from the uniform distribution
                           ! in the range (0,1). This will be used to pick the
                           ! categories.

    ! Output Variable
    integer, dimension(num_samples) :: &
      int_sample_category  ! An integer for each sample corresponding to the
                           ! category picked for the sample

    ! Local Variables
    integer :: sample,category      ! Looping variable(s)

    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_cumulative_probs

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    !--------------------------------------------------------------------------
    ! In order to facilitate picking categories for the sample points, a new
    ! array, category_cumulative_probs, is created.
    !
    ! Each element of category_cumulative_probs is simply the sum of the
    ! probabilities of all the previous categories. For example,
    !
    ! category_cumulative_probs(1) = 0
    ! category_cumulative_probs(2) = category_prescribed_probs(1)
    ! category_cumulative_probs(3) = category_prescribed_probs(1) + category_prescribed_probs(2)
    ! ...
    ! category_cumulative_probs(num_importance_categories) =
    !             category_prescribed_probs(1) + category_prescribed_probs(2) + ... +
    !             category_prescribed_probs(num_importance_categories-1)
    !--------------------------------------------------------------------------
    category_cumulative_probs(1) = zero
    do category=2, num_importance_categories
      category_cumulative_probs(category) = category_cumulative_probs(category-1) + &
                                            category_prescribed_probs(category-1)
    end do

    !--------------------------------------------------------------------------
    ! Pick categories based on the values of rand_vect.
    !--------------------------------------------------------------------------
    do sample=1, num_samples
      ! Initialize int_sample_category(sample) for error checking purposes.
      int_sample_category(sample) = 0
      do category=1, num_importance_categories

        if ( category < num_importance_categories ) then

          if ( rand_vect(sample) >= category_cumulative_probs(category) .and. &
               rand_vect(sample) <  category_cumulative_probs(category+1) ) then

            int_sample_category(sample) = category
            exit   ! Break out of the loop over categories, since we have found the category

          end if

        else if ( category == num_importance_categories ) then

          ! If the random number is greater than category_cumulative_probs(num_imp_categories-1)
          ! and less than 1, then it belongs in the last category
          if ( rand_vect(sample) >= category_cumulative_probs(category) .and. &
               rand_vect(sample) < one ) then
            int_sample_category(sample) = category
          end if

        end if ! category < num_importance_categories

      end do ! category=1, num_importance_categories-1

      ! We should have picked a category by now.
      if ( int_sample_category(sample) == 0 ) then
        write(fstderr,*) "Invalid rand_vect number in pick_sample_categories"
        write(fstderr,*) "rand_vect(sample) = ", rand_vect(sample)
        stop "Fatal error"
      end if

    end do ! sample=1, num_samples

    return
  end function pick_sample_categories
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine scale_sample_to_category( category, pdf_params, hydromet_pdf_params, &
                                       X_u_chi, X_u_dp1, X_u_dp2 )

  ! Description:
  !   Scale and transpose a sample point to reside in the specified category

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      dp,           &    ! Konstant(s)
      core_rknd

    use constants_clubb, only: &
      one_dp             ! Constant

    use pdf_parameter_module, only: &
      pdf_parameter

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter

    implicit none

    ! Input Variables
    type(importance_category_type), intent(in) :: &
      category             ! Scale the sample point to reside in this category

    type(pdf_parameter), intent(in) :: &
      pdf_params

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params

    ! Input/Output Variable
    ! These uniform samples, upon input, are uniformly distributed in the
    ! range (0,1).
    real( kind = dp ), intent(inout) :: &
      X_u_chi, & ! Uniform sample of extend cloud water mixing ratio
      X_u_dp1, & ! Uniform sample of the d+1 variate
      X_u_dp2    ! Uniform sample of the d+2 variate

    ! Local Variables
    real( kind = dp ) :: &
      cloud_frac_i, & 
      precip_frac_i, &
      mixt_frac

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    mixt_frac = real( pdf_params%mixt_frac, kind=dp )

    !--------------------------------------------------------
    ! Scale dp1 variate to be in component 1 or 2
    !--------------------------------------------------------
    if ( category%l_in_component_1 ) then

      ! Samples in component 1 have a dp1 variate that satisfies
      ! 0 < X_u_dp1 < mixt_frac.
      ! Scale X_u_dp1 to lie in (0,mixt_frac)
      X_u_dp1 = X_u_dp1 * mixt_frac

      ! Choose appropriate component cloud and precipitation fractions
      cloud_frac_i  = real( pdf_params%cloud_frac_1, kind=dp )
      precip_frac_i = real( hydromet_pdf_params%precip_frac_1, kind=dp )

    else  ! in component 2

      ! Scale and translate X_u_dp1 to lie in (mixt_frac, 1)
      X_u_dp1 = X_u_dp1 * (one_dp - mixt_frac) + mixt_frac

      ! Choose appropriate component cloud and precipitation fractions
      cloud_frac_i  = real( pdf_params%cloud_frac_2, kind=dp )
      precip_frac_i = real( hydromet_pdf_params%precip_frac_2, kind=dp )

    end if ! category%l_in_component_1

    !--------------------------------------------------------
    ! Scale dp2 variate to be in or out of precipitation
    !--------------------------------------------------------
    if ( category%l_in_precip ) then
      ! Samples in precipitation have a dp2 variate that satisfies
      ! 0 < X_u_dp2 < precip_frac_i
      ! Scale X_u_dp2 to lie in (0,precip_frac_i)
      X_u_dp2 = X_u_dp2 * precip_frac_i
    else
      ! Scale and translate X_u_dp2 to lie in (precip_frac_i,1)
      X_u_dp2 = X_u_dp2 * (one_dp - precip_frac_i) + precip_frac_i
    end if

    !--------------------------------------------------------
    ! Scale chi variate to be in or out of cloud
    !--------------------------------------------------------
    if ( category%l_in_cloud ) then
      ! Samples in cloud have a chi variate that satisfies
      ! (1.0 - cloud_frac_i) < X_u_chi < 1
      ! Scale and translate X_u_chi to lie in ( (1 - cloud_frac_i) , 1 )
      X_u_chi = X_u_chi * cloud_frac_i + (one_dp - cloud_frac_i)
    else
      ! Scale X_u_chi to lie in ( 0, (1 - cloud_frac_i) )
      X_u_chi = X_u_chi * (one_dp - cloud_frac_i)
    end if

    return
  end subroutine scale_sample_to_category
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function prescribe_importance_probs( importance_categories, category_real_probs ) &

  result( category_prescribed_probs )

  ! Description:
  !   Outputs a prescribed value for the probability of each importance category

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd     ! Constant

    use constants_clubb, only: &
      zero          ! Constant

    implicit none

    ! Local Constants
    real( kind = core_rknd ), parameter :: &
      prob_thresh = 5.0e-3_core_rknd

    ! Input Variables
    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories   ! A list of importance categories

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs     ! The real probability for each category

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs ! The prescribed probability for each category

    ! Local Variables
    integer :: icategory, jcategory
    real( kind = core_rknd ) :: nonzero_real_prob_sum, presc_prob_difference
    logical :: l_in_cloud, l_in_component_1, l_in_precip

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    do icategory=1, num_importance_categories

      l_in_cloud       = importance_categories(icategory)%l_in_cloud
      l_in_component_1 = importance_categories(icategory)%l_in_component_1
      l_in_precip      = importance_categories(icategory)%l_in_precip

      if ( l_in_cloud .and. l_in_precip .and. l_in_component_1 ) then
        category_prescribed_probs(icategory) = 0.15_core_rknd

      else if ( l_in_cloud .and. l_in_precip .and. (.not. l_in_component_1) ) then
        category_prescribed_probs(icategory) = 0.15_core_rknd

      else if ( l_in_cloud .and. (.not. l_in_precip) .and. l_in_component_1 ) then
        category_prescribed_probs(icategory) = 0.15_core_rknd

      else if ( l_in_cloud .and. (.not. l_in_precip) .and. (.not. l_in_component_1) ) then
        category_prescribed_probs(icategory) = 0.15_core_rknd

      else if ( (.not. l_in_cloud) .and. l_in_precip .and. l_in_component_1 ) then
        category_prescribed_probs(icategory) = 0.15_core_rknd

      else if ( (.not. l_in_cloud) .and. l_in_precip .and. (.not. l_in_component_1) ) then
        category_prescribed_probs(icategory) = 0.15_core_rknd

      else if ( (.not. l_in_cloud) .and. (.not. l_in_precip) .and. l_in_component_1 ) then
        category_prescribed_probs(icategory) = 0.05_core_rknd

      else if ( (.not. l_in_cloud) .and. (.not. l_in_precip) .and. (.not. l_in_component_1) ) then
        category_prescribed_probs(icategory) = 0.05_core_rknd

      else
        stop "Invalid category in prescribe_importance_probs"
      end if

    end do ! icategory=1, num_importance_categories

    !-----------------------------------------------------------------------------
    ! The following code ensures that we do not sample from categories that have
    ! no PDF weight
    !-----------------------------------------------------------------------------

    ! This value needs to be computed only once, so if it is needed below, it is computed
    ! on demand. For now, it is just set to zero.
    nonzero_real_prob_sum = zero

    do icategory=1, num_importance_categories

      if ( category_real_probs(icategory) < prob_thresh .and. &
           category_prescribed_probs(icategory) > category_real_probs(icategory) ) then

        ! Transfer all prescribed mass of this category (minus its PDF probability) to
        ! other categories.
        presc_prob_difference = category_prescribed_probs(icategory) - &
                                category_real_probs(icategory)

        ! Compute the sum of all categories with non-zero (greater than prob_thresh) PDF
        ! probability, iff this value is not already computed.
        if ( nonzero_real_prob_sum == zero ) then
          do jcategory=1, num_importance_categories
            if ( category_real_probs(jcategory) >= prob_thresh ) then
              nonzero_real_prob_sum = nonzero_real_prob_sum + category_real_probs(jcategory)
            end if
          end do
        end if

        ! An amount of prescribed probability equal to presc_prob_difference must be
        ! transferred to the other categories. We transfer to the other categories
        ! proportionally to the PDF probabilities of the other categories.
        do jcategory=1, num_importance_categories

          if ( category_real_probs(jcategory) >= prob_thresh ) then
            category_prescribed_probs(jcategory) = category_prescribed_probs(jcategory) + &
              ( presc_prob_difference * category_real_probs(jcategory) / nonzero_real_prob_sum )
          end if

        end do

        ! Do not perform importance_sampling on this category
        category_prescribed_probs(icategory) = category_real_probs(icategory)

      end if ! category_real_probs(icategory) < prob_thresh .and. ...

    end do ! icategory=1, num_importance_categories

    return
  end function prescribe_importance_probs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function eight_cluster_allocation( importance_categories, category_real_probs ) &

  result( category_prescribed_probs )

  ! Description:
  !   Clusters importance categories such that each of the eight importance
  !   categories has its own cluster. Effectively, there are no clusters.

  ! References:
  !   clubb:ticket:752
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd     ! Constant

    implicit none

    ! Local Constants
    integer, parameter :: &
      num_clusters = 8, &
      num_categories_per_cluster = 1

    !!! Prescribed probability definitions
    real( kind = core_rknd ), parameter :: &
      cloud_precip_comp1      = 0.15_core_rknd, &
      cloud_precip_comp2      = 0.15_core_rknd, &
      nocloud_precip_comp1    = 0.15_core_rknd, &
      nocloud_precip_comp2    = 0.15_core_rknd, &
      cloud_noprecip_comp1    = 0.15_core_rknd, &
      cloud_noprecip_comp2    = 0.15_core_rknd, &
      nocloud_noprecip_comp1  = 0.05_core_rknd, &
      nocloud_noprecip_comp2  = 0.05_core_rknd

    ! Input Variables
    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories   ! A list of importance categories

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs     ! The real probability for each category

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs ! The prescribed probability for each category

    ! Local Variables
    integer, dimension(num_clusters,num_categories_per_cluster) :: &
      cluster_categories

    real( kind = core_rknd ), dimension(num_clusters) :: &
      cluster_prescribed_probs

    logical :: l_in_cloud, l_in_component_1, l_in_precip

    integer :: icategory

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    do icategory=1, num_importance_categories

      cluster_categories(icategory,1) = icategory

      l_in_cloud       = importance_categories(icategory)%l_in_cloud
      l_in_component_1 = importance_categories(icategory)%l_in_component_1
      l_in_precip      = importance_categories(icategory)%l_in_precip

      if ( l_in_cloud .and. l_in_precip .and. l_in_component_1 ) then
        cluster_prescribed_probs(icategory) = cloud_precip_comp1

      else if ( l_in_cloud .and. l_in_precip .and. (.not. l_in_component_1) ) then
        cluster_prescribed_probs(icategory) = cloud_precip_comp2

      else if ( (.not. l_in_cloud) .and. l_in_precip .and. l_in_component_1 ) then
        cluster_prescribed_probs(icategory) = nocloud_precip_comp1

      else if ( (.not. l_in_cloud) .and. l_in_precip .and. (.not. l_in_component_1) ) then
        cluster_prescribed_probs(icategory) = nocloud_precip_comp2

      else if ( l_in_cloud .and. (.not. l_in_precip) .and. l_in_component_1 ) then
        cluster_prescribed_probs(icategory) = cloud_noprecip_comp1

      else if ( l_in_cloud .and. (.not. l_in_precip) .and. (.not. l_in_component_1) ) then
        cluster_prescribed_probs(icategory) = cloud_noprecip_comp2

      else if ( (.not. l_in_cloud) .and. (.not. l_in_precip) .and. l_in_component_1 ) then
        cluster_prescribed_probs(icategory) = nocloud_noprecip_comp1

      else if ( (.not. l_in_cloud) .and. (.not. l_in_precip) .and. (.not. l_in_component_1) ) then
        cluster_prescribed_probs(icategory) = nocloud_noprecip_comp2

      else
        stop "Invalid category in eight_cluster_allocation"
      end if

    end do ! icategory=1, num_importance_categories

    category_prescribed_probs = compute_clust_category_probs &
                                ( category_real_probs, num_clusters, num_categories_per_cluster, &
                                  cluster_categories, cluster_prescribed_probs )

    return
  end function eight_cluster_allocation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function compute_clust_category_probs &
           ( category_real_probs, num_clusters, num_categories_per_cluster, &
             cluster_categories, cluster_prescribed_probs ) &

  result( category_prescribed_probs )

  ! Description:
  !   This is a generalized algorithm that takes as input a set of "clusters"
  !   of the importance categories and a prescribed probability for each
  !   cluster, and computes the prescribed probabilities for each category
  !   such that the sum of the prescribed probabilities of every category
  !   within a cluster is equal to the prescribed probability of the cluster.

  ! References:
  !   clubb:ticket:752
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd       ! Constant

    use constants_clubb, only: &
      zero            ! Constant

    implicit none

    ! Local Constants
    real( kind = core_rknd ), parameter :: &
      prob_thresh = 5.0e-3_core_rknd

    ! Input Variables
    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs           ! The real probability for each category

    integer, intent(in) :: &
      num_clusters, &               ! The number of clusters to sample from
      num_categories_per_cluster    ! The number of categories in each cluster

    integer, dimension(num_clusters,num_categories_per_cluster), intent(in) :: &
      cluster_categories            ! An integer matrix containing indices corresponding
                                    ! to the members of the clusters

    real( kind = core_rknd ), dimension(num_clusters), intent(in) :: &
      cluster_prescribed_probs      ! Prescribed probability sum for each cluster

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs     ! Resulting prescribed probability for each individual category

    ! Local Variables
    real( kind = core_rknd ), dimension(num_clusters) :: &
      cluster_real_probs            ! Total PDF probability for each cluster

    real( kind = core_rknd ), dimension(num_clusters) :: &
      cluster_prescribed_probs_mod  ! Prescribed probability sum for each cluster, modified to
                                    ! take real probability thresholding into account

    logical, dimension(num_clusters) :: &
      l_cluster_presc_prob_modified ! Whether a given cluster prescribed probability was modified
                                    ! due to thresholding

    real( kind = core_rknd ) :: &
      nonzero_real_clust_sum, &     ! Sum of PDF probabilities for each non-modified cluster
      presc_prob_difference         ! Extra prescribed mass for given cluster to be distributed

    integer :: icluster, jcluster, icategory, cat_idx

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    ! Compute the total PDF probability for each cluster.
    cluster_real_probs(:) = zero
    do icluster=1, num_clusters
      do icategory=1, num_categories_per_cluster
        cluster_real_probs(icluster) = cluster_real_probs(icluster) + &
            category_real_probs(cluster_categories(icluster,icategory))
      end do
    end do

    ! Apply thresholding to ensure that clusters with extremely small PDF
    ! probability are not importance sampled.
    do icluster=1, num_clusters
      if ( cluster_real_probs(icluster) < prob_thresh .and. &
           cluster_prescribed_probs(icluster) > cluster_real_probs(icluster) ) then
        ! Thresholding is necessary for this cluster. The prescribed probability for
        ! this cluster will be set equal to the PDF probability of the cluster (that
        ! is, no importance sampling).
        cluster_prescribed_probs_mod(icluster) = cluster_real_probs(icluster)
        l_cluster_presc_prob_modified(icluster) = .true.
      else
        ! Thresholding is not necessary
        cluster_prescribed_probs_mod(icluster) = cluster_prescribed_probs(icluster)
        l_cluster_presc_prob_modified(icluster) = .false.
      end if ! cluster_real_probs(icluster) < prob_thresh .and. ...
    end do ! icluster=1, num_clusters

    ! Distribute any "extra" prescribed probability weight from thresholding to
    ! other clusters.
    if ( any( l_cluster_presc_prob_modified ) ) then

      ! Compute the sum of prescribed probabilities for all non-modified clusters
      nonzero_real_clust_sum = zero
      do icluster=1, num_clusters
        if ( .not. l_cluster_presc_prob_modified(icluster) ) then
          nonzero_real_clust_sum = nonzero_real_clust_sum + cluster_prescribed_probs_mod(icluster)
        end if
      end do ! icluster=1, num_clusters

      ! Transfer extra prescribed probability mass to other clusters.
      do icluster=1, num_clusters
        if ( l_cluster_presc_prob_modified(icluster) ) then

          presc_prob_difference = cluster_prescribed_probs_mod(icluster) - &
                                  cluster_prescribed_probs(icluster)

          do jcluster=1, num_clusters
            if ( .not. l_cluster_presc_prob_modified(jcluster) ) then
              cluster_prescribed_probs_mod(jcluster) = cluster_prescribed_probs_mod(jcluster) + &
                ( presc_prob_difference * cluster_real_probs(jcluster) / nonzero_real_clust_sum )
            end if
          end do

        end if ! l_cluster_presc_prob_modified(icluster)
      end do ! icluster=1, num_clusters

    end if ! any( l_cluster_presc_prob_modified )

    ! Finally, compute the prescribed probabilities for each category based on the cluster
    ! probabilities.
    do icluster=1, num_clusters
      do icategory=1, num_categories_per_cluster
        cat_idx = cluster_categories(icluster,icategory)
        if ( l_cluster_presc_prob_modified(icluster) ) then
          ! No scaling needs to be done, since the cluster's prescribed probability
          ! equals its PDF probability.
          category_prescribed_probs(cat_idx) = category_real_probs(cat_idx)
        else

          ! Scale category probability based on the cluster probability
          category_prescribed_probs(cat_idx) = category_real_probs(cat_idx) * &
            ( cluster_prescribed_probs_mod(icluster) / cluster_real_probs(icluster) )

        end if ! l_cluster_presc_prob_modified(icluster)
      end do ! icategory=1, num_categories_per_cluster
    end do ! icluster=1, num_clusters

    return
  end function compute_clust_category_probs
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  function cloud_importance_sampling( importance_categories, category_real_probs, &
                                      pdf_params ) &

  result( category_prescribed_probs )

  ! Description:
  !   Applies cloud weighted sampling such that approximately half of all
  !   sample points land in cloud and half land out of cloud!

  ! References:
  !   None :(
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd  ! Constant

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use constants_clubb, only: &
      one,  &    ! Constant(s)
      two,  &
      zero, &
      fstderr

    use pdf_utilities, only: &
      compute_mean_binormal ! Procedure

    implicit none

    ! Local Constants
    real( kind = core_rknd ), parameter :: &
      cloud_frac_min_samp = 0.001_core_rknd, &  ! Minimum cloud fraction for sampling
                                                ! preferentially within cloud
      cloud_frac_max_samp = 0.5_core_rknd       ! Maximum cloud fraction (exclusive) for
                                                ! sampling preferentially within cloud

    ! Input Variables
    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories   ! A list of importance categories

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs     ! The actual PDF probability for each category

    type(pdf_parameter), intent(in) :: &
      pdf_params

    ! Output Variable
    real( kind = core_rknd ), dimension(num_importance_categories) :: &
      category_prescribed_probs ! Probability of each category, scaled such that approximately half
                                ! of all sample points will appear in cloud

    ! Local Variables
    real( kind = core_rknd ) :: &
      cloud_frac
    integer :: icategory

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    cloud_frac = compute_mean_binormal( pdf_params%cloud_frac_1, pdf_params%cloud_frac_2, &
                                        pdf_params%mixt_frac )

    if ( cloud_frac >= cloud_frac_min_samp .and. cloud_frac < cloud_frac_max_samp ) then

      ! In-cloud categories ought to be divided by 2*cloud_frac
      ! Out-of-cloud categories should be divided by 2*(1-cloud_frac)

      do icategory=1, num_importance_categories
        if ( importance_categories(icategory)%l_in_cloud ) then
          category_prescribed_probs(icategory) = &
            category_real_probs(icategory) / (two*cloud_frac)
        else
          category_prescribed_probs(icategory) = &
            category_real_probs(icategory) / (two*(one-cloud_frac))
        end if ! importance_categories(icategory)%l_in_cloud
      end do

    else ! cloud_frac < cloud_frac_min_samp .or. cloud_frac >= cloud_frac_max_samp

      ! Do not perform cloud weighted sampling. Let the probabilities remain unmodified.
      category_prescribed_probs = category_real_probs

    end if ! cloud_frac < cloud_frac_min_samp .or. cloud_frac >= cloud_frac_max_samp

    return
  end function cloud_importance_sampling
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine importance_sampling_assertions &
             ( num_samples, importance_categories, category_real_probs, &
               category_prescribed_probs, category_sample_weights, X_u_chi_one_lev, &
               X_u_dp1_one_lev, X_u_dp2_one_lev, int_sample_category, &
               pdf_params, hydromet_pdf_params, &
               l_error )

  ! Description:
  !   Various assertion checks for importance sampling are performed here.

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd, &          ! Precision(s)
      dp

    use constants_clubb, only: &
      zero, &               ! Constant(s)
      one, &
      fstderr

    use pdf_parameter_module, only: &
      pdf_parameter

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter

    implicit none

    ! Input Variables
    integer :: &
      num_samples                         ! Number of SILHS sample points

    type(importance_category_type), dimension(num_importance_categories), intent(in) :: &
      importance_categories               ! The defined importance categories

    real( kind = core_rknd ), dimension(num_importance_categories), intent(in) :: &
      category_real_probs,        &       ! The real PDF probabilities for each category
      category_prescribed_probs,  &       ! Prescribed probability for each category
      category_sample_weights             ! Sample weight for each category

    real( kind = dp ), dimension(num_samples), intent(in) :: &
      X_u_chi_one_lev, &                  ! Samples of chi in uniform space
      X_u_dp1_one_lev, &                  ! Samples of the dp1 variate
      X_u_dp2_one_lev                     ! Samples of the dp2 variate

    integer, dimension(num_samples), intent(in) :: &
      int_sample_category                 ! An integer for each sample corresponding to the
                                          ! category picked for the sample

    type(pdf_parameter), intent(in) :: &
      pdf_params

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params

    ! Output Variables
    logical, intent(out) :: &
      l_error                             ! True if the assertion check fails.

    ! Local Variables
    real( kind = core_rknd ) :: &
      category_sum, tolerance

    real( kind = dp ) :: &
      cloud_frac_i, precip_frac_i

    integer :: isample

    type(importance_category_type) :: category

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    l_error = .false.

    !----------------------------------------------------------
    ! Assert that each of the probability groups sum to 1.0
    !----------------------------------------------------------
    tolerance = real( num_importance_categories, kind=core_rknd ) * epsilon( category_sum )

    category_sum = sum( category_real_probs )
    if ( abs( category_sum - one ) > tolerance ) then
      write(fstderr,*) "The real category probabilities do not sum to one."
      write(fstderr,*) "sum( category_real_probs ) = ", category_sum
      l_error = .true.
    end if

    category_sum = sum( category_prescribed_probs )
    if ( abs( category_sum - one ) > tolerance ) then
      write(fstderr,*) "The prescribed category probabilities do not sum to one."
      write(fstderr,*) "sum( category_prescribed_probs ) = ", category_sum
      l_error = .true.
    end if

    category_sum = sum( category_sample_weights * category_prescribed_probs )
    if ( abs( category_sum - one ) > tolerance ) then
      write(fstderr,*) "The weighted prescribed category probabilities do not sum to one."
      write(fstderr,*) "sum( category_sample_weights * category_prescribed_probs ) = ", category_sum
      l_error = .true.
    end if

    !---------------------------------------------------------------------
    ! Verify that samples have been correctly scaled to reside in the
    ! appropriate category
    !---------------------------------------------------------------------
    do isample = 1, num_samples

      category = importance_categories(int_sample_category(isample))

      ! Verification of component
      if ( category%l_in_component_1 ) then
        if ( X_u_dp1_one_lev(isample) > real( pdf_params%mixt_frac, kind=dp ) ) then
          write(fstderr,*) "The component of a sample is incorrect."
          l_error = .true.
        end if
        cloud_frac_i = real( pdf_params%cloud_frac_1, kind=dp )
        precip_frac_i = real( hydromet_pdf_params%precip_frac_1, kind=dp )
      else ! .not. category%l_in_component_1
        if ( X_u_dp1_one_lev(isample) < real( pdf_params%mixt_frac, kind=dp ) ) then
          write(fstderr,*) "The component of a sample is incorrect."
          l_error = .true.
        end if
        cloud_frac_i = real( pdf_params%cloud_frac_2, kind=dp )
        precip_frac_i = real( hydromet_pdf_params%precip_frac_2, kind=dp )
      end if ! category%l_in_component_1

      ! Verification of cloud
      if ( category%l_in_cloud ) then
        if ( X_u_chi_one_lev(isample) < (one - cloud_frac_i) ) then
          write(fstderr,*) "The chi element of a sample is incorrect."
          l_error = .true.
        end if
      else
        if ( X_u_chi_one_lev(isample) > (one - cloud_frac_i) ) then
          write(fstderr,*) "The chi element of a sample is incorrect."
          l_error = .true.
        end if
      end if

      ! Verification of precipitation
      if ( category%l_in_precip ) then
        if ( X_u_dp2_one_lev(isample) > precip_frac_i ) then
          write(fstderr,*) "The in-precipitation status of a sample is incorrect."
          l_error = .true.
        end if
      else
        if ( X_u_dp2_one_lev(isample) < precip_frac_i ) then
          write(fstderr,*) "The in-precipitation status of a sample is incorrect."
          l_error = .true.
        end if
      end if

    end do ! isample = 1, num_samples

    return
  end subroutine importance_sampling_assertions
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine cloud_weighted_sampling_driver &
             ( num_samples, p_matrix_chi, p_matrix_dp1, &
               cloud_frac_1, cloud_frac_2, mixt_frac, &
               X_u_chi, X_u_dp1, &
               lh_sample_point_weights, l_half_in_cloud )

  ! Description:
  !   Performs importance sampling such that half of sample points are in cloud

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd, &          ! Precision(s)
      dp

    use generate_lh_sample_module, only: &
      generate_uniform_sample ! Procedure

    use pdf_utilities, only: &
      compute_mean_binormal

    use constants_clubb, only: &
      one

    implicit none

    ! Parameter Constants
    real( kind = core_rknd ), parameter :: &
      cloud_frac_min_samp = 0.001_core_rknd, &  ! Minimum cloud fraction for sampling
                                                ! preferentially within cloud
      cloud_frac_max_samp = 0.5_core_rknd       ! Maximum cloud fraction (exclusive) for
                                                ! sampling preferentiallywithin cloud

    ! Input Variables
    integer :: &
      num_samples                         ! Number of SILHS sample points

    integer, dimension(num_samples), intent(in) :: &
      p_matrix_chi, &                     ! Permutation of the integers 0..num_samples
                                          ! from p_matrix for chi
      p_matrix_dp1                        ! Elements from p_matrix for dp1 element

    real( kind = core_rknd ), intent(in) :: &
      cloud_frac_1, &                     ! Cloud fraction in PDF component 1
      cloud_frac_2, &                     ! Cloud fraction in PDF component 2
      mixt_frac                           ! Weight of first gaussian component

    ! Input/Output Variables
    real( kind = dp ), dimension(num_samples), intent(inout) :: &
      X_u_chi, &                          ! Samples of chi in uniform space
      X_u_dp1                             ! Samples of the dp1 variate for determining mixture
                                          ! component

    ! Output Variables
    real( kind = core_rknd ), dimension(num_samples), intent(out) :: &
      lh_sample_point_weights             ! Weight of SILHS sample points (these must be applied
                                          ! when averaging results from, e.g., calling microphysics

    logical, intent(out) :: &
      l_half_in_cloud                     ! True if half of sample points are in cloud. Used
                                          ! for assertion checks later in the code

    ! Local Variables
    real( kind = core_rknd ) :: &
      cloud_frac                          ! Cloud fraction at k_lh_start

    real( kind = dp ), dimension(num_samples/2,1) :: &
      mixt_rand_cloud, &                  ! Stratified random variable for determining mixture
                                          ! component in cloud
      mixt_rand_clear                     ! Stratified random variable for determining mixture
                                          ! component in clear air

    real( kind = dp ) :: &
      mixt_rand_element

    real( kind = core_rknd ) :: &
      lh_sample_cloud_weight,     &       ! Weight of a cloudy sample point
      lh_sample_clear_weight              ! Weight of a clear  sample point

    integer, dimension(num_samples/2,1) :: &
      mixt_permuted_cloud, &              ! Permuted random numbers for mixt_rand_cloud
      mixt_permuted_clear                 ! Permuted random numbers for mixt_rand_clear

    integer :: sample, n_cloudy_samples, n_clear_samples
    logical :: l_cloudy_sample

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    cloud_frac = compute_mean_binormal( cloud_frac_1, cloud_frac_2, mixt_frac )

    if ( cloud_frac >= cloud_frac_min_samp .and. cloud_frac < cloud_frac_max_samp ) then

      n_cloudy_samples = 0
      n_clear_samples  = 0
      ! Pick two stratified random numbers for determining a mixture fraction.
      do sample=1, num_samples
        ! We arbitrarily choose that p_matrix elements less than num_samples/2
        ! will be used for mixt_rand_cloud
        if ( p_matrix_dp1(sample) < ( num_samples / 2 ) ) then
          n_cloudy_samples = n_cloudy_samples + 1
          mixt_permuted_cloud(n_cloudy_samples,1) = p_matrix_dp1(sample)
        else if ( p_matrix_dp1(sample) >= ( num_samples / 2 ) ) then
          n_clear_samples = n_clear_samples + 1
          mixt_permuted_clear(n_clear_samples,1) = p_matrix_dp1(sample) - (num_samples / 2)
        end if ! p_matrix_dp1(sample) < ( num_samples / 2 )
      end do
      ! Generate the stratified uniform numbers!
      call generate_uniform_sample( num_samples/2, num_samples/2, 1, mixt_permuted_cloud, & !In
                                    mixt_rand_cloud ) !Out
      call generate_uniform_sample( num_samples/2, num_samples/2, 1, mixt_permuted_clear, & !In
                                    mixt_rand_clear ) !Out

      n_cloudy_samples = 0
      n_clear_samples  = 0

      ! Compute weights for each type of point
      lh_sample_cloud_weight = 2._core_rknd * cloud_frac
      lh_sample_clear_weight = 2._core_rknd - lh_sample_cloud_weight

      do sample=1, num_samples

        ! Detect which half of the sample points are in clear air and which half are in
        ! the cloudy air
        if ( p_matrix_chi(sample) < ( num_samples / 2 ) ) then

          l_cloudy_sample = .false.
          lh_sample_point_weights(sample) = lh_sample_clear_weight
          n_clear_samples = n_clear_samples + 1
          mixt_rand_element = mixt_rand_clear(n_clear_samples,1)
        else

          l_cloudy_sample = .true.
          lh_sample_point_weights(sample) = lh_sample_cloud_weight
          n_cloudy_samples = n_cloudy_samples + 1
          mixt_rand_element = mixt_rand_cloud(n_cloudy_samples,1)
        end if

        ! Transpose and scale the points to be in or out of cloud
        call choose_X_u_scaled &
             ( l_cloudy_sample, & ! In
               p_matrix_chi(sample), num_samples, & ! In
               cloud_frac_1, cloud_frac_2, & ! In
               mixt_frac, mixt_rand_element, & !In
               X_u_dp1(sample), X_u_chi(sample) ) ! Out

      end do ! sample=1, num_samples

      l_half_in_cloud = .true.

    else ! cloud_frac < cloud_frac_min_samp .or. cloud_frac >= cloud_frac_max_samp

      ! Do not perform cloud weighted sampling.
      l_half_in_cloud = .false.
      lh_sample_point_weights(:) = one

    end if ! cloud_frac >= cloud_frac_min_samp .and. cloud_frac < cloud_frac_max_samp

    return
  end subroutine cloud_weighted_sampling_driver
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  function generate_strat_uniform_variate( num_samples ) &

  result( stratified_variate )

  ! Description:
  !   Generates a stratified uniform sample for a single variable

  ! References:
  !   None
  !-----------------------------------

    use clubb_precision, only: &
      core_rknd    ! Constant

    use permute_height_time_module, only: &
      rand_permute ! Procedure

    use generate_lh_sample_module, only: &
      choose_permuted_random

    implicit none

    ! Input Variable
    integer, intent(in) :: &
      num_samples ! Number of SILHS sample points

    ! Output Variable
    real( kind = core_rknd ), dimension(num_samples) :: &
      stratified_variate      ! Uniform samples stratified in (0,1)

    ! Local Variables

    ! Vector of the integers 0,1,2,...,num_samples in random order
    integer, dimension(num_samples) :: pvect

    integer :: sample

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    !-------------------------------------------------------------------
    ! Draw the integers 0,1,2,...,num_samples in random order
    !-------------------------------------------------------------------
    call rand_permute( num_samples, pvect )

    !----------------------------------------------------------------------
    ! For each permuted integer (each box), determine a random real number
    !----------------------------------------------------------------------
    do sample=1, num_samples
      stratified_variate(sample) = &
        real( choose_permuted_random( num_samples, pvect(sample) ), kind = core_rknd )
    end do

    return
  end function generate_strat_uniform_variate
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine choose_X_u_scaled &
             ( l_cloudy_sample, &
               p_matrix_element, num_samples, &
               cloud_frac_1, cloud_frac_2, &
               mixt_frac, mixt_rand_element, &
               X_u_dp1_element, X_u_chi_element )

! Description:
!   Find a clear or cloudy point for sampling.
!
! References:
!   None
!-------------------------------------------------------------------------------

    use generate_lh_sample_module, only: &
      choose_permuted_random    ! Procedure

    use clubb_precision, only: &
      core_rknd, & ! Variable(s)
      dp

    implicit none

    ! Input Variables
    logical, intent(in) :: &
      l_cloudy_sample ! Whether his is a cloudy or clear air sample point

    integer, intent(in) :: &
      p_matrix_element, & ! Integer from 0..num_samples for this sample
      num_samples       ! Total number of calls to the microphysics

    real( kind = core_rknd ), intent(in) :: &
      cloud_frac_1, &    ! Cloud fraction associated with mixture component 1     [-]
      cloud_frac_2, &    ! Cloud fraction associated with mixture component 2     [-]
      mixt_frac         ! Mixture fraction                                       [-]

    real( kind = dp ), intent(in) :: &
      mixt_rand_element ! Random number (0,1) for determining mixture component

    ! Output Variables
    real(kind=dp), intent(out) :: &
      X_u_dp1_element, X_u_chi_element ! Elements from X_u (uniform dist.)

    ! Local Variables
    real(kind=dp) :: cloud_frac_i, conditional_mixt_frac

    real( kind = dp ) :: &
      cld_comp1_frac,    &          ! Fraction of points in component 1 and cloud
      cld_comp2_frac,    &          ! Fraction of points in component 2 and cloud
      nocld_comp1_frac,  &          ! Fraction of points in component 1 and clear air
      nocld_comp2_frac              ! Fraction of points in component 2 and clear air

    integer :: X_mixt_comp_one_lev, p_matrix_element_ranged

    real( kind = dp ) :: mixt_rand_element_scaled, mixt_frac_dp, chi_rand_element

    ! ---- Begin code ----

    mixt_frac_dp = real( mixt_frac, kind=dp )

    cld_comp1_frac = real( mixt_frac*cloud_frac_1, kind=dp )
    cld_comp2_frac = real( (1._core_rknd-mixt_frac)*cloud_frac_2, kind=dp )

    !---------------------------------------------------------------------
    ! Determine the conditional mixture fraction, given whether we are in
    ! cloud
    !---------------------------------------------------------------------
    if ( l_cloudy_sample ) then

      conditional_mixt_frac = cld_comp1_frac / ( cld_comp1_frac + cld_comp2_frac )

    else ! .not. l_cloudy_sample

      nocld_comp1_frac = mixt_frac_dp - cld_comp1_frac
      nocld_comp2_frac = (1._dp - mixt_frac_dp) - cld_comp2_frac

      conditional_mixt_frac = nocld_comp1_frac / (nocld_comp1_frac + nocld_comp2_frac)

    end if ! l_cloudy_sample

    !---------------------------------------------------------------------
    ! Determine mixture component given the conditional mixture fraction
    !---------------------------------------------------------------------
!    if ( in_mixt_comp_1( mixt_rand_element, conditional_mixt_frac ) ) then
    if ( mixt_rand_element < conditional_mixt_frac ) then
      X_mixt_comp_one_lev = 1
    else
      X_mixt_comp_one_lev = 2
    end if

    !---------------------------------------------------------------------
    ! Determine dp1 element given mixture component
    !---------------------------------------------------------------------
    if ( X_mixt_comp_one_lev == 1 ) then
      ! mixt_rand_element is scaled to give real number stratified in (0,1)
      mixt_rand_element_scaled = mixt_rand_element / conditional_mixt_frac
      X_u_dp1_element = mixt_rand_element_scaled * mixt_frac_dp

    else if ( X_mixt_comp_one_lev == 2 ) then
      ! mixt_rand_element is scaled to give real number stratified in (0,1)
      mixt_rand_element_scaled = (mixt_rand_element - conditional_mixt_frac) / &
                                   (1._dp - conditional_mixt_frac) 
      X_u_dp1_element = mixt_rand_element_scaled * (1._dp - mixt_frac_dp) + mixt_frac_dp

    else
      stop "Should not be here"
    end if ! X_mixt_comp_one_lev == 1

    !---------------------------------------------------------------------
    ! Determine chi element given mixture component and l_cloudy_sample
    !---------------------------------------------------------------------
    ! Get p_matrix_element in the proper range
    if ( p_matrix_element >= ( num_samples / 2 ) ) then
      p_matrix_element_ranged = p_matrix_element - ( num_samples / 2 )
    else
      p_matrix_element_ranged = p_matrix_element
    end if
    
    ! Get a stratified number in (0,1) using the element of p_matrix!
    chi_rand_element = real( choose_permuted_random( num_samples/2, p_matrix_element_ranged ), &
                              kind=dp )

    ! Determine cloud fraction
    if ( X_mixt_comp_one_lev == 1 ) then
      cloud_frac_i = real( cloud_frac_1, kind=dp )
    else
      cloud_frac_i = real( cloud_frac_2, kind=dp )
    end if

    if ( l_cloudy_sample ) then
      ! Scale and translate sample point to reside in cloud
      X_u_chi_element = cloud_frac_i * chi_rand_element + &
                          (1._dp - cloud_frac_i )
    else
      ! Scale and translate sample point to reside in clear air (no cloud)
      X_u_chi_element = (1._dp-cloud_frac_i) * chi_rand_element
    end if ! l_cloudy_sample

    return
  end subroutine choose_X_u_scaled
!----------------------------------------------------------------------

end module silhs_importance_sample_module
