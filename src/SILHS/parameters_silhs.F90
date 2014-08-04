!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module parameters_silhs

! Description:
!   Parameters for SILHS!

! References:
!   None
!-------------------------------------------------------------------------

  implicit none

  ! Flags for the Latin Hypercube sampling code 
  logical, public :: &
    l_lh_cloud_weighted_sampling  = .true., & ! Limit noise by sampling in-cloud
    l_lh_vert_overlap = .true.                ! Assume maximum overlap for s_mellor (chi)

  !$omp threadprivate( l_lh_cloud_weighted_sampling, l_lh_vert_overlap )

  private ! Default Scope

end module parameters_silhs
