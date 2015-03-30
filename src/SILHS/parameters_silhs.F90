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
    l_Lscale_vert_avg = .true.               ! Calculate Lscale_vert_avg in lh_subcolumn_generator

  !$omp threadprivate( l_lh_cloud_weighted_sampling, l_Lscale_vert_avg )

  private ! Default Scope

end module parameters_silhs
