!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module parameters_silhs

! Description:
!   Parameters for SILHS!

! References:
!   None
!-------------------------------------------------------------------------

  use mt95, only: genrand_intg

  implicit none

  ! Constant Parameters
  integer, parameter, public :: &
    lh_microphys_interactive     = 1, & ! Feed the samples into the microphysics and allow feedback
    lh_microphys_non_interactive = 2, & ! Feed the samples into the microphysics with no feedback
    lh_microphys_disabled        = 3    ! Disable Latin hypercube entirely

  ! Flags for the Latin Hypercube sampling code 
  logical, public :: &
    l_lh_cloud_weighted_sampling  = .true., & ! Limit noise by sampling in-cloud
    l_lh_vert_overlap = .true.                ! Assume maximum overlap for s_mellor (chi)

  !$omp threadprivate( l_lh_cloud_weighted_sampling, l_lh_vert_overlap )

  integer, public :: &
    lh_microphys_calls = 2, & ! Number of latin hypercube samples to call the microphysics with
    lh_sequence_length = 1    ! Number of timesteps before the latin hypercube seq. repeats

  integer(kind=genrand_intg), public :: &
    lh_seed = 5489_genrand_intg ! Seed for the Mersenne Twister

!$omp threadprivate( lh_microphys_calls, lh_sequence_length, lh_seed )

  ! Determines how the latin hypercube samples should be used (microphysics/radiation)
  integer, public :: &
    lh_microphys_type = 3 ! Default to non-interactive

!$omp threadprivate( lh_microphys_type )

  private ! Default Scope

end module parameters_silhs
