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
    LH_microphys_interactive     = 1, & ! Feed the samples into the microphysics and allow feedback
    LH_microphys_non_interactive = 2, & ! Feed the samples into the microphysics with no feedback
    LH_microphys_disabled        = 3    ! Disable Latin hypercube entirely

  ! Morrison aerosol parameters
  integer, parameter, public :: &
    morrison_no_aerosol = 0, &
    morrison_power_law  = 1, &
    morrison_lognormal  = 2

  ! Local Variables
  logical, public :: & 
    l_cloud_sed = .false.,         & ! Cloud water sedimentation (K&K/No microphysics)
    l_ice_micro  = .false.,        & ! Compute ice (COAMPS/Morrison)
    l_upwind_diff_sed = .false.,   & ! Use upwind differencing approx for sedimentation (K&K/COAMPS)
    l_graupel  = .false.,          & ! Compute graupel (COAMPS/Morrison)
    l_hail  = .false.,             & ! Assumption about graupel/hail? (Morrison)
    l_seifert_beheng  = .false.,   & ! Use Seifert and Behneng warm drizzle (Morrison)
    l_predictnc  = .true.,         & ! Predict cloud droplet conconcentration (Morrison)
    l_const_Nc_in_cloud = .false., & ! Use a constant cloud droplet conc. within cloud (K&K)
    l_subgrid_w = .true.,          & ! Use subgrid w (Morrison)
    l_arctic_nucl = .false.,       & ! Use MPACE observations (Morrison)
    l_fix_pgam = .false.,          & ! Fix pgam (Morrison)
    l_in_cloud_Nc_diff = .true.,   & ! Use in cloud values of Nc for diffusion
    l_var_covar_src = .false.        ! Flag for using upscaled microphysics source terms
                                     ! for predictive variances and covariances (KK micro)

  ! KK microphysics has a feature that clips rrainm_source down to the maximum allowable value if it
  ! grows so large during a timestep that it over-depletes from cloud water (making rcm negative).
  ! In KK upscaled microphysics, in which the microphysics equations are analytically integrated
  ! over the entire grid-box, this clipping is only performed for the mean value of rrainm_source.
  ! However, in Latin hypercube sampling, this adjustment is performed for every sample point. Thus,
  ! in order to test convergence of SILHS to the KK analytic solution, we must turn this adjustment
  ! off for individual sample points, and only do it for the mean. This flag does that when it is
  ! set to true. (ticket:558)
  logical, public :: &
    l_silhs_KK_convergence_adj_mean = .false.

!$omp threadprivate( l_cloud_sed, l_ice_micro, l_graupel, l_hail, &
!$omp                l_upwind_diff_sed, l_seifert_beheng, l_predictnc, &
!$omp                l_const_Nc_in_cloud, l_subgrid_w, l_arctic_nucl, &
!$omp                l_fix_pgam, l_in_cloud_Nc_diff, l_var_covar_src, &
!$omp                l_silhs_KK_convergence_adj_mean )

  logical, public :: & 
    l_cloud_edge_activation = .false., & ! Activate on cloud edges (Morrison)
    l_local_kk              = .false.    ! Local drizzle for Khairoutdinov & Kogan microphysics

!$omp threadprivate(l_cloud_edge_activation, l_local_kk)

  character(len=30), public :: &
    specify_aerosol = "morrison_lognormal"  ! Specify aerosol (Morrison)

  ! Flags for the Latin Hypercube sampling code 
  logical, public :: &
    l_fix_s_t_correlations = .true., &        ! Use a fixed correlation for s and t Mellor
    l_lh_cloud_weighted_sampling  = .true., & ! Limit noise by sampling in-cloud
    l_lh_vert_overlap = .true.                ! Assume maximum overlap for s_mellor

!$omp threadprivate( l_fix_s_t_correlations, l_lh_cloud_weighted_sampling, &
!$omp                l_lh_vert_overlap )

  integer, public :: &
    LH_microphys_calls = 2, & ! Number of latin hypercube samples to call the microphysics with
    LH_sequence_length = 1    ! Number of timesteps before the latin hypercube seq. repeats

  integer(kind=genrand_intg), public :: &
    LH_seed = 5489_genrand_intg ! Seed for the Mersenne

!$omp threadprivate( LH_microphys_calls, LH_sequence_length, LH_seed )

  ! Determines how the latin hypercube samples should be used with the microphysics
  integer, public :: &
    LH_microphys_type = 0

!$omp threadprivate( LH_microphys_type )

  character(len=50), public :: &
    micro_scheme = "none" ! khairoutdinv_kogan, simplified_ice, coamps, etc.

!$omp threadprivate( micro_scheme )

  character(len=10), dimension(:), allocatable, public :: & 
    hydromet_list

!$omp threadprivate( hydromet_list )

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

  ! Statistical rain parameters        .

  ! Prescribed parameters for hydrometeor in-precip values of
  ! sigma_hm_i^2 / mu_hm_i^2 at grid levels that have some cloud.
  ! They can be set based on values for a given case in the CASE_model.in file.
  real( kind = core_rknd ), public :: &
    rr_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_rr_i^2/mu_rr_i^2  [-]
    Nr_sigma2_on_mu2_ip_cloud = 1.0_core_rknd    ! sigma_Nr_i^2/mu_Nr_i^2  [-]

!$omp threadprivate( rr_sigma2_on_mu2_ip_cloud, Nr_sigma2_on_mu2_ip_cloud )

  ! Parameters added for ice microphysics and latin hypercube sampling
  real( kind = core_rknd ), public :: &
    rs_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_rs_i^2/mu_rs_i^2  [-]
    Ns_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_Ns_i^2/mu_Ns_i^2  [-]
    ri_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_ri_i^2/mu_ri_i^2  [-]
    Ni_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_Ni_i^2/mu_Ni_i^2  [-]
    rg_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_rg_i^2/mu_rg_i^2  [-]
    Ng_sigma2_on_mu2_ip_cloud = 1.0_core_rknd    ! sigma_Ng_i^2/mu_Ng_i^2  [-]

!$omp threadprivate( rs_sigma2_on_mu2_ip_cloud, Ns_sigma2_on_mu2_ip_cloud, & 
!$omp                ri_sigma2_on_mu2_ip_cloud, Ni_sigma2_on_mu2_ip_cloud, &
!$omp                rg_sigma2_on_mu2_ip_cloud, Ng_sigma2_on_mu2_ip_cloud )

  ! Prescribed parameters for hydrometeor in-precip values of
  ! sigma_hm_i^2 / mu_hm_i^2 at grid levels that are entirely clear.
  ! They can be set based on values for a given case in the CASE_model.in file.
  real( kind = core_rknd ), public :: &
    rr_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_rr_i^2/mu_rr_i^2  [-]
    Nr_sigma2_on_mu2_ip_below = 1.0_core_rknd    ! sigma_Nr_i^2/mu_Nr_i^2  [-]

!$omp threadprivate( rr_sigma2_on_mu2_ip_below, Nr_sigma2_on_mu2_ip_below )

  ! Parameters added for ice microphysics and latin hypercube sampling
  real( kind = core_rknd ), public :: &
    rs_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_rs_i^2/mu_rs_i^2  [-]
    Ns_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_Ns_i^2/mu_Ns_i^2  [-]
    ri_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_ri_i^2/mu_ri_i^2  [-]
    Ni_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_Ni_i^2/mu_Ni_i^2  [-]
    rg_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_rg_i^2/mu_rg_i^2  [-]
    Ng_sigma2_on_mu2_ip_below = 1.0_core_rknd    ! sigma_Ng_i^2/mu_Ng_i^2  [-]

!$omp threadprivate( rs_sigma2_on_mu2_ip_below, Ns_sigma2_on_mu2_ip_below, & 
!$omp                ri_sigma2_on_mu2_ip_below, Ni_sigma2_on_mu2_ip_below, &
!$omp                rg_sigma2_on_mu2_ip_below, Ng_sigma2_on_mu2_ip_below )

  ! Prescribed parameter for <N_cn'^2> / <N_cn>^2 at any grid level.
  ! The parameter can be set based on the value for a given case in the
  ! CASE_model.in file.
  real( kind = core_rknd ), public :: &
    Ncnp2_on_Ncnm2 = 1.0_core_rknd   ! Prescribed ratio <N_cn'^2>/<N_cn>^2 [-]

!$omp threadprivate( Ncnp2_on_Ncnm2 )

  ! Other needed parameters
  ! Khairoutdinov and Kogan (2000) ratio of drizzle drop mean geometric
  ! radius to drizzle drop mean volume radius.
  ! Khairoutdinov and Kogan (2000); p. 233
  real( kind = core_rknd ), public :: C_evap = 0.86_core_rknd ! Khairoutdinov and Kogan (2000) 

  !real, public :: C_evap = 0.86*0.2 ! COAMPS value of KK C_evap
  !real, public :: C_evap = 0.55     ! KK 2000, Marshall-Palmer (1948) value.

  real( kind = core_rknd ), public :: r_0 = 25.0e-6_core_rknd ! Assumed radius of all new drops; m.
  ! Value specified in KK (2000); p. 235.
  ! Vince Larson set r_0=28mum to agree with COAMPS-LES formula. 15 April 2005
  !REAL, PARAMETER:: r_0 = 28.0e-6   ! Assumed radius of all new drops; m.
  !                                  ! Value that COAMPS LES has in it.
  !REAL, PARAMETER:: r_0 = 30.0e-6   ! Assumed radius of all new drops; m.
  !                                  ! Khairoutdinov said it was okay!
  ! End Vince Larson's change.

!$omp threadprivate( C_evap, r_0 )

  ! Values of exponents in KK microphysics
  real( kind = core_rknd ), public :: &
    KK_evap_Supersat_exp = one,  & ! Exponent on Supersaturation (S) in KK evap. eq.; 1
    KK_evap_rr_exp = one_third,  & ! Exponent on r_r in KK evaporation eq.; 1/3
    KK_evap_Nr_exp = two_thirds, & ! Exponent on N_r in KK evaporation eq.; 2/3
    KK_auto_rc_exp = 2.47_core_rknd, & ! Exponent on r_c in KK autoconversion eq.; 2.47
    KK_auto_Nc_exp = -1.79_core_rknd, & ! Exponent on N_c in KK autoconversion eq.; -1.79
    KK_accr_rc_exp = 1.15_core_rknd, & ! Exponent on r_c in KK accretion eq.; 1.15
    KK_accr_rr_exp = 1.15_core_rknd, & ! Exponent on r_r in KK accretion eq.; 1.15
    KK_mvr_rr_exp  = one_third, & ! Exponent on r_r in KK mean volume radius
    KK_mvr_Nr_exp  = -one_third, & ! Exponent on N_r in KK mean volume radius
    KK_Nrm_evap_nu = one    ! Exponent (parameter) in <N_r> evaporation eq.; 1

!$omp threadprivate( KK_evap_Supersat_exp, KK_evap_rr_exp, KK_evap_Nr_exp, &
!$omp                KK_auto_rc_exp, KK_auto_Nc_exp, KK_accr_rc_exp, &
!$omp                KK_accr_rr_exp, KK_mvr_rr_exp, KK_mvr_Nr_exp, &
!$omp                KK_Nrm_evap_nu )

  real( kind = core_rknd ), dimension(:), allocatable, public :: &
    hydromet_tol    ! Tolerance values for all hydrometeors    [units vary]

!$omp threadprivate( hydromet_tol )   

  private ! Default Scope


end module parameters_microphys
