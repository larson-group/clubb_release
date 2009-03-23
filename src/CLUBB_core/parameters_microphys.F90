!-------------------------------------------------------------------------------
! $Id$
module parameters_microphys

! Description:
!   Parameters for microphysical schemes

! References:
! None
!-------------------------------------------------------------------------------
  use stats_precision, only: &
    time_precision

  implicit none

  logical, public :: & 
    l_cloud_sed,       & ! Cloud water sedimentation
    l_ice_micro,       & ! Compute ice
    l_graupel,         & ! Compute graupel
    l_hail,            & ! 
    l_seifert_beheng,  & ! Use Seifert and Behneng warm drizzle
    l_predictnc,       & ! Predict cloud droplet number conc
    l_specify_aerosol, & ! Specify aerosol
    l_subgrid_w,       & ! Use subgrid w 
    l_arctic_nucl,     & ! Use MPACE observations
    l_cloud_edge_activation, &! Activate on cloud edges
    l_fix_pgam  ! Fix pgam

!$omp threadprivate(l_cloud_sed, l_ice_micro, l_graupel, l_hail, &
!$omp   l_seifert_beheng, l_predictnc, l_specify_aerosol, l_subgrid_w, &
!$omp   l_arctic_nucl, l_cloud_edge_activation, l_fix_pgam)


  character(len=50), public :: &
    micro_scheme ! khairoutdinv_kogan, simplified_ice, coamps, etc.

!$omp threadprivate(micro_scheme)

  character(len=10), dimension(:), allocatable, public :: & 
    hydromet_list

!$omp threadprivate(hydromet_list)

  real(kind=time_precision), public :: &
    microphys_start_time  ! When to start the microphysics      [s]

!$omp threadprivate(microphys_start_time)

  real, public :: &
    Ncm_initial ! Initial cloud droplet number concentration

!$omp threadprivate(Ncm_initial)

  private ! Default Scope


end module parameters_microphys
