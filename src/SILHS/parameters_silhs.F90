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

  ! Flags for the SILHS sampling code 
  logical, public :: &
    l_lh_importance_sampling  = .true., &     ! Limit noise by performing importance sampling
    l_Lscale_vert_avg         = .true., &     ! Calculate Lscale_vert_avg in lh_subcolumn_generator
    l_lh_straight_mc          = .false.       ! Use true Monte Carlo sampling with no Latin
                                              ! hypercube sampling and no importance sampling.

  !$omp threadprivate( l_lh_importance_sampling, l_Lscale_vert_avg, l_lh_straight_mc )

  private ! Default Scope

end module parameters_silhs
