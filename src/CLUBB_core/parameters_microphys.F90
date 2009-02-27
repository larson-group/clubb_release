!-------------------------------------------------------------------------------
! $Id$
module parameters_microphys

! Description:
!   Parameters for microphysical schemes

! References:
! None
!-------------------------------------------------------------------------------
  implicit none

  logical, public :: & 
    l_cloud_sed,     & ! Cloud water sedimentation
    l_kk_rain,       & ! K&K microphysics
    l_coamps_micro,  & ! COAMPS microphysics
    l_icedfs           ! Simplified ice

  character(len=10), dimension(:), allocatable, public :: & 
    hydromet_list

  private ! Default Scope

!$omp threadprivate(l_cloud_sed, l_kk_rain, l_coamps_micro, l_icedfs, &
!$omp   hydromet_list)

end module parameters_microphys
