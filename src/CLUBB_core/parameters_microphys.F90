!-------------------------------------------------------------------------------
! $Id$
module parameters_microphys

! Description:
!   Parameters for microphysical schemes

! References:
!   None
!-------------------------------------------------------------------------------
  use stats_precision, only: &
    time_precision

  implicit none

  logical, public :: & 
    l_cloud_sed,       & ! Cloud water sedimentation (K&K/No microphysics)
    l_ice_micro,       & ! Compute ice (COAMPS/Morrison)
    l_graupel,         & ! Compute graupel (COAMPS/Morrison)
    l_hail,            & ! Assumption about graupel/hail? (Morrison)
    l_seifert_beheng,  & ! Use Seifert and Behneng warm drizzle (Morrison)
    l_predictnc,       & ! Predict cloud droplet number conc (Morrison)
    l_specify_aerosol, & ! Specify aerosol (Morrison)
    l_subgrid_w,       & ! Use subgrid w (Morrison)
    l_arctic_nucl,     & ! Use MPACE observations (Morrison)
    l_fix_pgam           ! Fix pgam (Morrison)

!$omp threadprivate(l_cloud_sed, l_ice_micro, l_graupel, l_hail, &
!$omp   l_seifert_beheng, l_predictnc, l_specify_aerosol, l_subgrid_w, &
!$omp   l_arctic_nucl, l_fix_pgam)

  logical, public :: & 
    l_cloud_edge_activation,  & ! Activate on cloud edges (Morrison)
    l_latin_hypercube_sampling  ! Latin Hypercube Sampling (K&K)
!$omp threadprivate(l_cloud_edge_activation, l_latin_hypercube_sampling)


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

  logical, parameter, public :: &
    l_local_kk = .false.  ! Local drizzle for Khairoutdinov & Kogan microphysics


  private ! Default Scope


end module parameters_microphys
