!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module parameters_microphys

! Description:
!   Parameters for microphysical schemes

! References:
!   None
!-------------------------------------------------------------------------

  use clubb_precision, only: &
    time_precision, &
    core_rknd

  use mt95, only: &
    genrand_intg

  use constants_clubb, only: &
    one, &
    one_third, &
    two_thirds

  implicit none

  ! Constant Parameters
  integer, parameter, public :: &
    lh_microphys_interactive     = 1, & ! Feed the samples into the microphysics and allow feedback
    lh_microphys_non_interactive = 2, & ! Feed the samples into the microphysics with no feedback
    lh_microphys_disabled        = 3    ! Disable Latin hypercube entirely

  ! Morrison aerosol parameters
  integer, parameter, public :: &
    morrison_no_aerosol = 0, &
    morrison_power_law  = 1, &
    morrison_lognormal  = 2

  ! Local Variables
  logical, public :: & 
    l_cloud_sed = .false.,         & ! Cloud water sedimentation (K&K/No microphysics)
    l_ice_microphys  = .false.,    & ! Compute ice (COAMPS/Morrison)
    l_upwind_diff_sed = .false.,   & ! Use upwind differencing approx for sedimentation (K&K/COAMPS)
    l_graupel  = .false.,          & ! Compute graupel (COAMPS/Morrison)
    l_hail  = .false.,             & ! Assumption about graupel/hail? (Morrison)
    l_seifert_beheng  = .false.,   & ! Use Seifert and Behneng warm drizzle (Morrison)
    l_predict_Nc = .false.,        & ! Predict cloud droplet conconcentration (Morrison)
    l_subgrid_w = .true.,          & ! Use subgrid w (Morrison)
    l_arctic_nucl = .false.,       & ! Use MPACE observations (Morrison)
    l_fix_pgam = .false.,          & ! Fix pgam (Morrison)
    l_in_cloud_Nc_diff = .true.,   & ! Use in cloud values of Nc for diffusion
    l_var_covar_src = .false.        ! Flag for using upscaled microphysics source terms
                                     ! for predictive variances and covariances (KK)

  ! KK microphysics has a feature that clips rrm_source down to the maximum allowable value if it
  ! grows so large during a timestep that it over-depletes from cloud water (making rcm negative).
  ! In KK upscaled microphysics, in which the microphysics equations are analytically integrated
  ! over the entire grid-box, this clipping is only performed for the mean value of rrm_source.
  ! However, in Latin hypercube sampling, this adjustment is performed for every sample point. Thus,
  ! in order to test convergence of SILHS to the KK analytic solution, we must turn this adjustment
  ! off for individual sample points, and only do it for the mean. This flag does that when it is
  ! set to true. (ticket:558)
  logical, public :: &
    l_silhs_KK_convergence_adj_mean = .false.

!$omp threadprivate( l_cloud_sed, l_ice_microphys, l_graupel, l_hail, &
!$omp                l_upwind_diff_sed, l_seifert_beheng, l_predict_Nc, &
!$omp                l_subgrid_w, l_arctic_nucl, &
!$omp                l_fix_pgam, l_in_cloud_Nc_diff, l_var_covar_src, &
!$omp                l_silhs_KK_convergence_adj_mean )

  logical, public :: & 
    l_cloud_edge_activation = .false., & ! Activate on cloud edges (Morrison)
    l_local_kk              = .false.    ! Local drizzle for Khairoutdinov & Kogan microphysics

!$omp threadprivate(l_cloud_edge_activation, l_local_kk)

  character(len=30), public :: &
    specify_aerosol = "morrison_lognormal"  ! Specify aerosol (Morrison)

  integer, public :: &
    lh_num_samples = 2, &     ! Number of latin hypercube samples to call the microphysics (or
                              ! radiation) scheme with
    lh_sequence_length = 1    ! Number of timesteps before the latin hypercube seq. repeats

  integer(kind=genrand_intg), public :: &
    lh_seed = 5489_genrand_intg ! Seed for the Mersenne

!$omp threadprivate( lh_num_samples, lh_sequence_length, lh_seed )

  ! Determines how the latin hypercube samples should be used with the microphysics
  integer, public :: &
    lh_microphys_type = 3

!$omp threadprivate( lh_microphys_type )

  character(len=50), public :: &
    microphys_scheme = "none" ! khairoutdinv_kogan, simplified_ice, coamps, etc.

!$omp threadprivate( microphys_scheme )

  logical, dimension(:), allocatable, public :: &
    l_hydromet_sed    ! Flag to sediment mean hydrometeor fields

!$omp threadprivate( l_hydromet_sed )

  logical, public :: l_gfdl_activation    ! Flag for GFDL activation code

!$omp threadprivate( l_gfdl_activation )

  real(kind=time_precision), public :: &
    microphys_start_time = 0._time_precision  ! When to start the microphysics      [s]

!$omp threadprivate( microphys_start_time )

  real( kind = core_rknd ), public :: &
    Nc0_in_cloud = 100.e6_core_rknd    ! Initial cloud droplet concentration  [num/m^3]

!$omp threadprivate( Nc0_in_cloud )

  real( kind = core_rknd ), public :: &
    sigma_g  = 1.5_core_rknd ! Geometric std. dev. of cloud droplets falling in a stokes regime.

!$omp threadprivate( sigma_g )

  private ! Default Scope

end module parameters_microphys
