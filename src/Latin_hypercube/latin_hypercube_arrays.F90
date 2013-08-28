! $Id$
module latin_hypercube_arrays

  use clubb_precision, only: &
    dp, & ! double precision
    core_rknd

  implicit none

  public :: cleanup_latin_hypercube_arrays

  private

  logical, public :: &
    l_fixed_corr_initialized = .false.

!$omp threadprivate(l_fixed_corr_initialized)

  real( kind = dp ), allocatable, dimension(:,:), target, public :: &
    corr_stw_cloud_Cholesky, & ! Cholesky factorization of the correlation matrix
    corr_stw_below_Cholesky    ! Cholesky factorization of the correlation matrix

!$omp threadprivate(corr_stw_cloud_Cholesky, corr_stw_below_Cholesky)

  real( kind = dp ), allocatable, dimension(:), public :: &
    corr_stw_cloud_scaling, & ! Scaling factors for the correlation matrix [-]
    corr_stw_below_scaling    ! Scaling factors for the correlation matrix [-]

!$omp threadprivate(corr_stw_cloud_scaling, corr_stw_below_scaling)

  logical, public :: &
    l_corr_stw_cloud_scaling, & ! Whether we're scaling the correlation matrix
    l_corr_stw_below_scaling

!$omp threadprivate(l_corr_stw_cloud_scaling, l_corr_stw_below_scaling)

  integer, allocatable, dimension(:,:,:), public :: & 
    height_time_matrix ! matrix of rand ints

!$omp threadprivate(height_time_matrix)

  contains

  !-----------------------------------------------------------------------------
  subroutine cleanup_latin_hypercube_arrays( )

    ! Description:
    !   De-allocate latin hypercube arrays
    ! References:
    !   None
    !---------------------------------------------------------------------------
    implicit none

    ! External
    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( corr_stw_cloud_Cholesky ) ) then
      deallocate( corr_stw_cloud_Cholesky )
    end if

    if ( allocated( corr_stw_below_Cholesky ) ) then
      deallocate( corr_stw_below_Cholesky )
    end if

    if ( allocated( corr_stw_cloud_scaling ) ) then
      deallocate( corr_stw_cloud_scaling )
    end if

    if ( allocated( corr_stw_below_scaling ) ) then
      deallocate( corr_stw_below_scaling )
    end if

    if ( allocated( height_time_matrix ) ) then
      deallocate( height_time_matrix )
    end if

    return
  end subroutine cleanup_latin_hypercube_arrays

end module latin_hypercube_arrays
