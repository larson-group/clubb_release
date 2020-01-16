!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module parameters_silhs

! Description:
!   Parameters for SILHS!

! References:
!   None
!-------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd     ! Constant

  implicit none

  ! Cluster allocation strategies!!!
  integer, parameter, public :: &
    ! All eight categories, effectively no clustering
    eight_cluster_allocation_opt = 1, &
    ! Four clusters for the combinations of cloud/no cloud and component 1/2.
    ! Precipitation fraction is ignored.
    four_cluster_allocation_opt  = 2, &
    ! Two clusters, one containing all categories with either cloud or precip,
    ! and the other containing categories with neither
    two_cluster_cp_nocp_opt      = 3

  ! The following type defines parameters that control the sample point
  ! allocation for the clustered sampling scheme
  ! (l_lh_clustered_sampling = .true.).
  type eight_cluster_presc_probs_type

    real( kind = core_rknd ) :: &
      cloud_precip_comp1      = 0.15_core_rknd, &
      cloud_precip_comp2      = 0.15_core_rknd, &
      nocloud_precip_comp1    = 0.15_core_rknd, &
      nocloud_precip_comp2    = 0.15_core_rknd, &
      cloud_noprecip_comp1    = 0.15_core_rknd, &
      cloud_noprecip_comp2    = 0.15_core_rknd, &
      nocloud_noprecip_comp1  = 0.05_core_rknd, &
      nocloud_noprecip_comp2  = 0.05_core_rknd

  end type eight_cluster_presc_probs_type

  ! Flags for the SILHS sampling code
  type silhs_config_flags_type

    integer :: &
      cluster_allocation_strategy
    logical :: &
        l_lh_importance_sampling,   &  ! Do importance sampling
        l_Lscale_vert_avg,          &  ! Vertically average Lscale
        l_lh_straight_mc,           &  ! Do not apply LH or importance sampling at all
        l_lh_clustered_sampling,    &  ! Use prescribed probability sampling with clusters
        l_rcm_in_cloud_k_lh_start,  &  ! Determine k_lh_start based on maximum within-cloud rcm
        l_random_k_lh_start,        &  ! k_lh_start found randomly between max rcm and rcm_in_cloud
        l_max_overlap_in_cloud,     &  ! Use maximum vertical overlap in cloud
        l_lh_instant_var_covar_src, &  ! Produce instantaneous var/covar tendencies
        l_lh_limit_weights,         &  ! Ensure weights stay under a given value
        l_lh_var_frac,              &  ! Prescribe variance fractions
        l_lh_normalize_weights         ! Normalize weights to sum to num_samples

  end type silhs_config_flags_type

  type(eight_cluster_presc_probs_type), public, save :: &
    eight_cluster_presc_probs                 ! Prescribed probabilities for
                                              ! l_lh_clustered_sampling = .true.

  !$omp threadprivate( eight_cluster_presc_probs )

  real( kind = core_rknd ), public :: &
    importance_prob_thresh = 1.0e-8_core_rknd, & ! Minimum PDF probability of category for
                                                 ! importance sampling
    vert_decorr_coef       = 0.03_core_rknd      ! Empirically defined de-correlation constant [-]

  !$omp threadprivate( importance_prob_thresh, vert_decorr_coef )
  
  
  real( kind = core_rknd ), public, parameter :: &
    single_prec_thresh   = 3.e-8_core_rknd       ! Uniform samples are expected to be in the range
                                                 ! [3.e-8_core_rknd,1-3.e-8_core_rknd] since the 
                                                 ! algorithm used to calculate the inverse cdf is 
                                                 ! only accurate for single precision values

  private ! Default Scope

  public :: eight_cluster_presc_probs_type, silhs_config_flags_type, &
            set_default_silhs_config_flags, initialize_silhs_config_flags_type, &
            print_silhs_config_flags

  contains

!-------------------------------------------------------------------------------
  subroutine set_default_silhs_config_flags( cluster_allocation_strategy, &
                                             l_lh_importance_sampling, &
                                             l_Lscale_vert_avg, &
                                             l_lh_straight_mc, &
                                             l_lh_clustered_sampling, &
                                             l_rcm_in_cloud_k_lh_start, &
                                             l_random_k_lh_start, &
                                             l_max_overlap_in_cloud, &
                                             l_lh_instant_var_covar_src, &
                                             l_lh_limit_weights, &
                                             l_lh_var_frac, &
                                             l_lh_normalize_weights )

    ! Description:
    !   Sets all SILHS flags to a default setting.

    ! References:
    !   None
    !---------------------------------------------------------------------------

    implicit none

    ! Output variables
    integer, intent(out) :: &
      cluster_allocation_strategy   ! Two clusters, one containing all categories with either
                                    ! cloud or precip, and the other containing categories with
                                    ! neither

    logical, intent(out) :: &
      l_lh_importance_sampling, &   ! Limit noise by performing importance sampling
      l_Lscale_vert_avg, &          ! Calculate Lscale_vert_avg in generate_silhs_sample
      l_lh_straight_mc, &           ! Use true Monte Carlo sampling with no Latin
                                    !  hypercube sampling and no importance sampling
      l_lh_clustered_sampling, &    ! Use the "new" SILHS importance sampling
                                    !  scheme with prescribed probabilities
      l_rcm_in_cloud_k_lh_start, &  ! Determine k_lh_start based on maximum within-cloud rcm
      l_random_k_lh_start, &        ! Place k_lh_start at a random grid level between
                                    !  maximum rcm and maximum rcm_in_cloud
      l_max_overlap_in_cloud, &     ! Assume maximum vertical overlap when grid-box rcm
                                    !  exceeds cloud threshold
      l_lh_instant_var_covar_src, & ! Produces "instantaneous" variance-covariance
                                    !  microphysical source terms, ignoring
                                    !  discretization effects
      l_lh_limit_weights, &         ! Limit SILHS sample point weights for stability
      l_lh_var_frac, &              ! Prescribe variance fractions
      l_lh_normalize_weights        ! Scale sample point weights to sum to num_samples
                                    ! (the "ratio estimate")

!-----------------------------------------------------------------------
    ! Begin code

    cluster_allocation_strategy = two_cluster_cp_nocp_opt
    l_lh_importance_sampling = .true.   ! Limit noise by performing importance sampling
    l_Lscale_vert_avg = .true.          ! Calculate Lscale_vert_avg in generate_silhs_sample
    l_lh_straight_mc = .false.          ! Use true Monte Carlo sampling with no Latin
                                        !  hypercube sampling and no importance sampling
    l_lh_clustered_sampling = .true.    ! Use the "new" SILHS importance sampling
                                        !  scheme with prescribed probabilities
    l_rcm_in_cloud_k_lh_start = .true.  ! Determine k_lh_start based on maximum within-cloud rcm
    l_random_k_lh_start = .false.       ! Place k_lh_start at a random grid level between
                                        !  maximum rcm and maximum rcm_in_cloud
    l_max_overlap_in_cloud = .true.     ! Assume maximum vertical overlap when grid-box rcm
                                        !  exceeds cloud threshold
    l_lh_instant_var_covar_src = .true. ! Produces "instantaneous" variance-covariance
                                        !  microphysical source terms, ignoring
                                        !  discretization effects
    l_lh_limit_weights = .true.         ! Limit SILHS sample point weights for stability
    l_lh_var_frac = .false.             ! Prescribe variance fractions
    l_lh_normalize_weights = .true.     ! Scale sample point weights to sum to num_samples
                                        ! (the "ratio estimate")

    return
  end subroutine set_default_silhs_config_flags
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine initialize_silhs_config_flags_type( cluster_allocation_strategy, &
                                                 l_lh_importance_sampling, &
                                                 l_Lscale_vert_avg, &
                                                 l_lh_straight_mc, &
                                                 l_lh_clustered_sampling, &
                                                 l_rcm_in_cloud_k_lh_start, &
                                                 l_random_k_lh_start, &
                                                 l_max_overlap_in_cloud, &
                                                 l_lh_instant_var_covar_src, &
                                                 l_lh_limit_weights, &
                                                 l_lh_var_frac, &
                                                 l_lh_normalize_weights, &
                                                 silhs_config_flags )

    ! Description:
    !   Initialize the silhs_config_flags_type.

    ! References:
    !   None
    !---------------------------------------------------------------------------

    implicit none

    ! Input variables
    integer, intent(in) :: &
      cluster_allocation_strategy   ! Two clusters, one containing all categories with either
                                    ! cloud or precip, and the other containing categories with
                                    ! neither

    logical, intent(in) :: &
      l_lh_importance_sampling, &   ! Limit noise by performing importance sampling
      l_Lscale_vert_avg, &          ! Calculate Lscale_vert_avg in generate_silhs_sample
      l_lh_straight_mc, &           ! Use true Monte Carlo sampling with no Latin
                                    !  hypercube sampling and no importance sampling
      l_lh_clustered_sampling, &    ! Use the "new" SILHS importance sampling
                                    !  scheme with prescribed probabilities
      l_rcm_in_cloud_k_lh_start, &  ! Determine k_lh_start based on maximum within-cloud rcm
      l_random_k_lh_start, &        ! Place k_lh_start at a random grid level between
                                    !  maximum rcm and maximum rcm_in_cloud
      l_max_overlap_in_cloud, &     ! Assume maximum vertical overlap when grid-box rcm
                                    !  exceeds cloud threshold
      l_lh_instant_var_covar_src, & ! Produces "instantaneous" variance-covariance
                                    !  microphysical source terms, ignoring
                                    !  discretization effects
      l_lh_limit_weights, &         ! Limit SILHS sample point weights for stability
      l_lh_var_frac, &              ! Prescribe variance fractions
      l_lh_normalize_weights        ! Scale sample point weights to sum to num_samples
                                    ! (the "ratio estimate")

    ! Output variables
    type(silhs_config_flags_type), intent(out) :: &
      silhs_config_flags            ! Derived type holding all configurable SILHS flags

!-----------------------------------------------------------------------
    ! Begin code

    silhs_config_flags%cluster_allocation_strategy = cluster_allocation_strategy
    silhs_config_flags%l_lh_importance_sampling    = l_lh_importance_sampling
    silhs_config_flags%l_Lscale_vert_avg           = l_Lscale_vert_avg
    silhs_config_flags%l_lh_straight_mc            = l_lh_straight_mc
    silhs_config_flags%l_lh_clustered_sampling     = l_lh_clustered_sampling
    silhs_config_flags%l_rcm_in_cloud_k_lh_start   = l_rcm_in_cloud_k_lh_start
    silhs_config_flags%l_random_k_lh_start         = l_random_k_lh_start
    silhs_config_flags%l_max_overlap_in_cloud      = l_max_overlap_in_cloud
    silhs_config_flags%l_lh_instant_var_covar_src  = l_lh_instant_var_covar_src
    silhs_config_flags%l_lh_limit_weights          = l_lh_limit_weights
    silhs_config_flags%l_lh_var_frac               = l_lh_var_frac
    silhs_config_flags%l_lh_normalize_weights      = l_lh_normalize_weights

    return
  end subroutine initialize_silhs_config_flags_type
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine print_silhs_config_flags( iunit, silhs_config_flags )

    ! Description:
    !   Prints the silhs_config_flags.

    ! References:
    !   None
    !---------------------------------------------------------------------------

    implicit none

    ! Input variables
    integer, intent(in) :: &
      iunit ! The file to write to

    type(silhs_config_flags_type), intent(in) :: &
      silhs_config_flags ! Derived type holding all configurable SILHS flags

!-----------------------------------------------------------------------
    ! Begin code

    write(iunit,*) "cluster_allocation_strategy = ", silhs_config_flags%cluster_allocation_strategy
    write(iunit,*) "l_lh_importance_sampling = ", silhs_config_flags%l_lh_importance_sampling
    write(iunit,*) "l_Lscale_vert_avg = ", silhs_config_flags%l_Lscale_vert_avg
    write(iunit,*) "l_lh_straight_mc = ", silhs_config_flags%l_lh_straight_mc
    write(iunit,*) "l_lh_clustered_sampling = ", silhs_config_flags%l_lh_clustered_sampling
    write(iunit,*) "l_rcm_in_cloud_k_lh_start = ", silhs_config_flags%l_rcm_in_cloud_k_lh_start
    write(iunit,*) "l_random_k_lh_start = ", silhs_config_flags%l_random_k_lh_start
    write(iunit,*) "l_max_overlap_in_cloud = ", silhs_config_flags%l_max_overlap_in_cloud
    write(iunit,*) "l_lh_instant_var_covar_src = ", silhs_config_flags%l_lh_instant_var_covar_src
    write(iunit,*) "l_lh_limit_weights = ", silhs_config_flags%l_lh_limit_weights
    write(iunit,*) "l_lh_var_frac = ", silhs_config_flags%l_lh_var_frac
    write(iunit,*) "l_lh_normalize_weights = ", silhs_config_flags%l_lh_normalize_weights

    return
  end subroutine print_silhs_config_flags
!-----------------------------------------------------------------------

end module parameters_silhs
