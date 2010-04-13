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
    l_fix_pgam,        & ! Fix pgam (Morrison)
    l_in_cloud_Nc_diff   ! Use in cloud values of Nc for diffusion

!$omp threadprivate(l_cloud_sed, l_ice_micro, l_graupel, l_hail, &
!$omp   l_seifert_beheng, l_predictnc, l_specify_aerosol, l_subgrid_w, &
!$omp   l_arctic_nucl, l_fix_pgam, l_in_cloud_Nc_diff)

  logical, public :: & 
    l_cloud_edge_activation,    & ! Activate on cloud edges (Morrison)
    l_local_kk                    ! Local drizzle for Khairoutdinov & Kogan microphysics

!$omp threadprivate(l_cloud_edge_activation, l_local_kk)

  ! Flags for the Latin Hypercube sampling code (currently configured to run
  ! in diagnostic mode)
  logical, public :: & 
    l_latin_hypercube_sampling, &   ! Latin Hypercube Sampling (K&K)
    l_lh_cloud_weighted_sampling, & ! Limit noise by sampling in-cloud
    l_lh_vert_overlap               ! Assume maximum overlap for s_mellor

!$omp threadprivate(l_latin_hypercube_sampling, l_lh_cloud_weighted_sampling, &
!$omp   l_lh_vert_overlap)

  real, allocatable, dimension(:), public :: &
    LH_sample_point_weights

!$omp threadprivate(LH_sample_point_weights)

  integer, public :: &
    LH_microphys_calls, & ! Number of latin hypercube samples to call the microphysics with
    LH_sequence_length    ! Number of timesteps before the latin hypercube seq. repeats
!$omp threadprivate(LH_microphys_calls,LH_sequence_length)

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
    Ncm_initial ! Initial cloud droplet number concentration [#/cc]

!$omp threadprivate(Ncm_initial)

  ! Statistical rain parameters        .

  ! Parameters for in-cloud (from SAM RF02 DO).
  real, public :: &       ! RF02 value
    rrp2_on_rrainm2_cloud, & ! 0.766
    Nrp2_on_Nrm2_cloud,    & ! 0.429
    Ncp2_on_Ncm2_cloud,    & ! 0.003
    corr_rrNr_LL_cloud, & ! 0.786
    corr_srr_NL_cloud,  & ! 0.242
    corr_sNr_NL_cloud,  & ! 0.285
    corr_sNc_NL_cloud     ! 0.433

  ! Parameters for below-cloud (from SAM RF02 DO).
  real, public :: &       ! RF02 value
    rrp2_on_rrainm2_below, & ! 8.97
    Nrp2_on_Nrm2_below,    & ! 12.03
    Ncp2_on_Ncm2_below,    & ! 0.00 ! Not applicable below cloud.
    corr_rrNr_LL_below, & ! 0.886
    corr_srr_NL_below,  & ! 0.056
    corr_sNr_NL_below,  & ! 0.015
    corr_sNc_NL_below     ! 0.00 ! Not applicable below cloud.

  ! Other needed parameters
  real, public :: C_evap ! 0.86    ! Khairoutdinov and Kogan (2000) ratio of
  ! drizzle drop mean geometric radius to
  ! drizzle drop mean volume radius.
  ! Khairoutdinov and Kogan (2000); p. 233.
  !real, public :: C_evap = 0.86*0.2 ! COAMPS value of KK C_evap
  !real, public :: C_evap = 0.55     ! KK 2000, Marshall-Palmer (1948) value.

  real, public :: r_0 ! 25.0e-6   ! Assumed radius of all new drops; m.
  ! Value specified in KK (2000); p. 235.
  ! Vince Larson set r_0=28mum to agree with COAMPS-LES formula. 15 April 2005
  !REAL, PARAMETER:: r_0 = 28.0e-6   ! Assumed radius of all new drops; m.
  !                                  ! Value that COAMPS LES has in it.
  !REAL, PARAMETER:: r_0 = 30.0e-6   ! Assumed radius of all new drops; m.
  !                                  ! Khairoutdinov said it was okay!
  ! End Vince Larson's change.

!$omp threadprivate( rrp2_on_rrainm2_cloud, Nrp2_on_Nrm2_cloud, Ncp2_on_Ncm2_cloud, &
!$omp   corr_rrNr_LL_cloud, corr_srr_NL_cloud,  corr_sNr_NL_cloud,  corr_sNc_NL_cloud, &
!$omp   rrp2_on_rrainm2_below, Nrp2_on_Nrm2_below, Ncp2_on_Ncm2_below, &
!$omp   corr_rrNr_LL_below, corr_srr_NL_below, corr_sNr_NL_below, corr_sNc_NL_below, &
!$omp   C_evap, r_0 )

  private ! Default Scope


end module parameters_microphys
