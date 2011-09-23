!------------------------------------------------------------------------------
! $Id$
module microphys_driver

  ! Description:
  ! Call a microphysical scheme to compute hydrometeor species, and advect,
  ! sediment, & diffuse using a tridiagonal system.

  ! References:
  !  H. Morrison, J. A. Curry, and V. I. Khvorostyanov, 2005: A new double-
  !  moment microphysics scheme for application in cloud and
  !  climate models. Part 1: Description. J. Atmos. Sci., 62, 1665-1677.

  !  Khairoutdinov, M. and Kogan, Y.: A new cloud physics parameterization in a
  !  large-eddy simulation model of marine stratocumulus, Mon. Wea. Rev., 128,
  !  229-243, 2000.
  !-----------------------------------------------------------------------------

  use parameters_microphys, only: &
    l_in_cloud_Nc_diff,           & ! Use in cloud values of Nc for diffusion
    l_cloud_sed,                  & ! Cloud water sedimentation (K&K or no microphysics)
    l_ice_micro,                  & ! Compute ice (COAMPS / Morrison)
    l_graupel,                    & ! Compute graupel (Morrison)
    l_hail,                       & ! See module_mp_graupel for a description
    l_seifert_beheng,             & ! Use Seifert and Beheng (2001) warm drizzle (Morrison)
    l_predictnc,                  & ! Predict cloud droplet number conc (Morrison)
    specify_aerosol,              & ! Specify aerosol (Morrison)
    l_subgrid_w,                  & ! Use subgrid w  (Morrison)
    l_arctic_nucl,                & ! Use MPACE observations (Morrison)
    l_cloud_edge_activation,      & ! Activate on cloud edges (Morrison)
    l_fix_pgam,                   & ! Fix pgam (Morrison)
    l_fix_s_t_correlations,       & ! Use a fixed correlation for s and t Mellor (Latin Hypercube)
    l_lh_vert_overlap,            & ! Assume maximum overlap for s_mellor (Latin Hypercube)
    l_lh_cloud_weighted_sampling, & ! Sample preferentially within cloud (Latin Hypercube)
    LH_microphys_calls,           & ! # of latin hypercube samples to call the microphysics with 
    LH_sequence_length,           & ! Number of timesteps before the latin hypercube seq. repeats
    l_local_kk,                   & ! Use local formula for K&K
    l_upwind_diff_sed,            & ! Use the upwind differencing approx. for sediementation
    micro_scheme,                 & ! The microphysical scheme in use
    hydromet_list,                & ! Names of the hydrometeor species
    microphys_start_time,         & ! When to start the microphysics [s]
    Ncm_initial                     ! Initial value for Ncm (K&K, l_cloud_sed, Morrison)

  use parameters_microphys, only: &
    LH_microphys_interactive,     & ! Feed the subcolumns into the microphysics and allow feedback
    LH_microphys_non_interactive, & ! Feed the subcolumns into the microphysics with no feedback
    LH_microphys_disabled           ! Disable latin hypercube entirely

  use constants_clubb, only: &
    cloud_frac_min

  use phys_buffer, only: & ! Used for placing wp2_zt in morrison-gettelman microphysics
    pbuf_init,           &
    pbuf_add,            &
    pbuf_allocate,       &
    pbuf_deallocate

  implicit none

  ! Subroutines
  public :: advance_microphys, init_microphys, cleanup_microphys

  private :: microphys_lhs, microphys_solve
  private :: adj_microphys_tndcy

  ! Functions
  private :: sed_centered_diff_lhs, sed_upwind_diff_lhs

  ! Variables
  logical, private, allocatable, dimension(:) :: l_hydromet_sed ! Whether to sediment variables
  logical :: l_gfdl_activation ! Flag for GFDL activation code

  private ! Default Scope

  contains

  !-------------------------------------------------------------------------------
  subroutine init_microphys( iunit, runtype, namelist_file, case_info_file, &
                             hydromet_dim )

    ! Description:
    ! Set indices to the various hydrometeor species and define hydromet_dim for
    ! the purposes of allocating memory.

    ! References:
    ! None
    !---------------------------------------------------------------------------

    use grid_class, only: & 
      gr

    use array_index, only: & 
      iirrainm, iiNrm, iirsnowm, iiricem, iirgraupelm, & ! Variables
      iiNcm, iiNsnowm, iiNim, iiNgraupelm

    use parameters_microphys, only: &
      morrison_no_aerosol, &  ! Constants
      morrison_power_law,  &
      morrison_lognormal

    use parameters_microphys, only: &
      rrp2_on_rrainm2_cloud, & ! Variable(s)
      Nrp2_on_Nrm2_cloud,    &
      Ncp2_on_Ncm2_cloud,    &
      corr_rrNr_LL_cloud,    &
      corr_rtrr_NL_cloud,    &
      corr_rtNr_NL_cloud,    &
      corr_rtNc_NL_cloud,    &
      corr_thlrr_NL_cloud,   &
      corr_thlNr_NL_cloud,   &
      corr_thlNc_NL_cloud,   &
      corr_srr_NL_cloud,     &
      corr_sNr_NL_cloud,     &
      corr_sNc_NL_cloud,     &
      rrp2_on_rrainm2_below, &
      Nrp2_on_Nrm2_below,    &
      Ncp2_on_Ncm2_below,    &
      corr_rrNr_LL_below,    &
      corr_rtrr_NL_below,    &
      corr_rtNr_NL_below,    &
      corr_rtNc_NL_below,    &
      corr_thlrr_NL_below,   &
      corr_thlNr_NL_below,   &
      corr_thlNc_NL_below,   &
      corr_srr_NL_below,     &
      corr_sNr_NL_below,     &
      corr_sNc_NL_below,     &
      C_evap,                &
      r_0

    use parameters_microphys, only: &
      rsnowp2_on_rsnowm2_cloud, & ! Variables
      Nsnowp2_on_Nsnowm2_cloud, & 
      ricep2_on_ricem2_cloud, & 
      Nicep2_on_Nicem2_cloud, &
      rsnowp2_on_rsnowm2_below, & 
      Nsnowp2_on_Nsnowm2_below, & 
      ricep2_on_ricem2_below, & 
      Nicep2_on_Nicem2_below

    use parameters_microphys, only: &
      LH_microphys_type_int => LH_microphys_type ! Determines how the LH samples are used

    ! The version of the Morrison 2005 microphysics that is in SAM.
    use module_mp_GRAUPEL, only: &
      Nc0, ccnconst, ccnexpnt, & ! Variables
      aer_rm1, aer_rm2, aer_n1, aer_n2, &
      aer_sig1, aer_sig2, pgam_fixed, &
      doicemicro, &         ! use ice species (snow/cloud ice/graupel)
      dograupel, &          ! use graupel
      dohail, &             ! make graupel species have properties of hail
      dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization
      dopredictNc, &        ! prediction of cloud droplet number
      aerosol_mode, &       ! specify two modes of (sulfate) aerosol
      dosubgridw, &         ! input estimate of subgrid w to microphysics
      doarcticicenucl, &    ! use arctic parameter values for ice nucleation
      docloudedgeactivation,& ! activate cloud droplets throughout the cloud
      dofix_pgam            ! option to fix value of pgam (exponent in cloud water gamma distn)

    use module_mp_Graupel, only: &
      GRAUPEL_INIT ! Subroutine
      
    use cldwat2m_micro, only: &
      ini_micro ! Subroutine
      
    use microp_aero, only: &
      ini_microp_aero ! Subroutine

    use constants_clubb, only: &
      fstderr,   & ! Constant
      cm3_per_m3, &
      micron_per_m

    use text_writer, only: &
      write_text   ! Used to write microphysics settings to setup.txt file

    use error_code, only: clubb_at_least_debug_level ! Function

    use gfdl_activation, only: nooc, sul_concen, & ! Variables
      low_concen, high_concen, &
      lowup, highup, lowup2, highup2, &
      lowmass2, highmass2, lowmass3, highmass3, &
      lowmass4, highmass4, lowmass5, highmass5, &
      lowT2, highT2, aeromass_value

    use aer_ccn_act_k_mod, only: aer_ccn_act_k_init ! Procedure

    use gfdl_activation, only: Loading ! Procedure

    use stats_precision, only:  & 
        time_precision ! Variable(s)

#ifdef LATIN_HYPERCUBE
    use latin_hypercube_arrays, only: &
      setup_corr_varnce_array   ! Procedure(s)
#endif /* LATIN_HYPERCUBE */

    implicit none

    ! Constant Parameters
    logical, parameter :: &
      l_write_to_file = .true. ! If true, will write case information to a file.

#ifdef LATIN_HYPERCUBE
    character(len=21), parameter :: &
      LH_input_path = "../input/case_setups/"
#endif

    ! External
    intrinsic :: trim

    ! Input variables
    integer, intent(in) :: &
      iunit ! File unit

    character(len=*), intent(in) :: &
      runtype

    character(len=*), intent(in) :: &
      namelist_file ! File name

    character(len=*), intent(in) :: &
      case_info_file ! Name of simulation info file (plain text)

    ! Output variables
    integer, intent(out) :: & 
      hydromet_dim ! Number of hydrometeor fields.

    ! Local variables 
    character(len=30) :: LH_microphys_type
    integer, parameter :: res = 20   ! Used for lookup tables with GFDL activation
    integer, parameter :: res2 = 20  ! Used for lookup tables with GFDL activation
    real, dimension( res, res, res, res, res ) :: &
      droplets, droplets2            ! Used for lookup tables with GFDL activation

    namelist /microphysics_setting/ &
      micro_scheme, l_cloud_sed, &
      l_ice_micro, l_graupel, l_hail, l_upwind_diff_sed, &
      l_seifert_beheng, l_predictnc, specify_aerosol, l_subgrid_w, &
      l_arctic_nucl, l_cloud_edge_activation, l_fix_pgam, l_in_cloud_Nc_diff, &
      LH_microphys_type, l_local_kk, LH_microphys_calls, LH_sequence_length, &
      l_lh_cloud_weighted_sampling, l_fix_s_t_correlations, l_lh_vert_overlap, &
      rrp2_on_rrainm2_cloud, Nrp2_on_Nrm2_cloud, Ncp2_on_Ncm2_cloud, &
      corr_rrNr_LL_cloud, corr_rtrr_NL_cloud, corr_rtNr_NL_cloud, &
      corr_rtNc_NL_cloud, corr_thlrr_NL_cloud, corr_thlNr_NL_cloud, &
      corr_thlNc_NL_cloud, corr_srr_NL_cloud, corr_sNr_NL_cloud, &
      corr_sNc_NL_cloud, rrp2_on_rrainm2_below, &
      Nrp2_on_Nrm2_below, Ncp2_on_Ncm2_below, &
      corr_rrNr_LL_below, corr_rtrr_NL_below, corr_rtNr_NL_below, &
      corr_rtNc_NL_below, corr_thlrr_NL_below, corr_thlNr_NL_below, &
      corr_thlNc_NL_below, corr_srr_NL_below, &
      corr_sNr_NL_below, corr_sNc_NL_below, &
      C_evap, r_0, microphys_start_time, &
      Ncm_initial, ccnconst, ccnexpnt, aer_rm1, aer_rm2, &
      aer_n1, aer_n2, aer_sig1, aer_sig2, pgam_fixed

    namelist /gfdl_activation_setting/ &
      nooc, sul_concen, low_concen, high_concen, &
      lowup, highup, lowup2, highup2, lowmass2, &
      highmass2, lowmass3, highmass3,  &
      lowmass4, highmass4, lowmass5, highmass5, &
      lowT2, highT2, aeromass_value, l_gfdl_activation

    ! ---- Begin Code ----

    ! Set default values, then read in the namelist

    micro_scheme = "none"

    l_gfdl_activation = .false.

    ! Cloud water sedimentation from the RF02 case
    ! This has no effect on Morrison's cloud water sedimentation
    l_cloud_sed = .false.

    !---------------------------------------------------------------------------
    ! Parameters for Khairoutdinov and Kogan microphysics
    !---------------------------------------------------------------------------
    l_local_kk = .false. ! Use the local parameterization for K&K

    C_evap = 0.86    ! Khairoutdinov and Kogan (2000) ratio of
    ! drizzle drop mean geometric radius to
    ! drizzle drop mean volume radius.
    ! Khairoutdinov and Kogan (2000); p. 233.
    !C_evap = 0.86*0.2 ! COAMPS value of KK C_evap
    !C_evap = 0.55     ! KK 2000, Marshall-Palmer (1948) value.

    r_0 = 25.0e-6   ! Assumed radius of all new drops; m.

    !---------------------------------------------------------------------------
    ! Parameters for Khairoutdinov and Kogan microphysics analytic solution
    ! (local_kk = .false.), or latin hypercube sampling (using either Khairoutdinov
    !  Kogan or Morrison microphysics).
    !---------------------------------------------------------------------------
    ! Parameters for in-cloud (from SAM RF02 DO).
    rrp2_on_rrainm2_cloud = 0.766
    Nrp2_on_Nrm2_cloud    = 0.429
    Ncp2_on_Ncm2_cloud    = 0.003
    corr_rrNr_LL_cloud    = 0.786
    corr_rtrr_NL_cloud    = 0.268
    corr_rtNr_NL_cloud    = 0.295
    corr_rtNc_NL_cloud    = 0.413
    corr_thlrr_NL_cloud   = -0.185
    corr_thlNr_NL_cloud   = -0.258
    corr_thlNc_NL_cloud   = -0.420
!    corr_srr_NL_cloud     = 0.242
!    corr_sNr_NL_cloud     = 0.285
!    corr_sNc_NL_cloud     = 0.433
    ! Parameters for below-cloud (from SAM RF02 DO).
    rrp2_on_rrainm2_below = 8.97
    Nrp2_on_Nrm2_below    = 12.03
    Ncp2_on_Ncm2_below    = 0.00  ! Not applicable below cloud.
    corr_rrNr_LL_below    = 0.886
    corr_rtrr_NL_below    = -0.005
    corr_rtNr_NL_below    = -0.041
    corr_rtNc_NL_below    = 0.00  ! Not applicable below cloud.
    corr_thlrr_NL_below   = -0.241
    corr_thlNr_NL_below   = -0.217
    corr_thlNc_NL_below   = 0.00  ! Not applicable below cloud.
!    corr_srr_NL_below     = 0.056
!    corr_sNr_NL_below     = 0.015
!    corr_sNc_NL_below     = 0.00  ! Not applicable below cloud.
    ! Other needed parameters

    ! Made up values for the variance of ice/snow, since we currently lack data
    ! for this.
    rsnowp2_on_rsnowm2_cloud = 0.766
    Nsnowp2_on_Nsnowm2_cloud = 0.429
    ricep2_on_ricem2_cloud = 1.0
    Nicep2_on_Nicem2_cloud = 1.0

    rsnowp2_on_rsnowm2_below = 0.766
    Nsnowp2_on_Nsnowm2_below = 0.429
    ricep2_on_ricem2_below = 1.0
    Nicep2_on_Nicem2_below = 1.0

    ! MPACE-A values for the correlation of ice/snow
!   corr_srsnow_NL_cloud     = 0.24
!   corr_srsnow_NL_below     = corr_srsnow_NL_cloud
!   corr_sNsnow_NL_cloud     = 0.29
!   corr_sNsnow_NL_below     = corr_sNsnow_NL_cloud
!   corr_rsnowNsnow_LL_cloud = 0.88
!   corr_rsnowNsnow_LL_below = corr_rsnowNsnow_LL_cloud
!   corr_srice_NL_cloud      = 0.31
!   corr_srice_NL_below      = corr_srice_NL_cloud
!   corr_sNi_NL_cloud        = 0.36
!   corr_sNi_NL_below        = corr_sNi_NL_cloud
!   corr_riceNi_LL_cloud     = 0.01
!   corr_riceNi_LL_below     = corr_riceNi_LL_cloud
!   corr_sNc_NL_cloud        = 0.67
!   corr_sNc_NL_below        = 0.67

!   corr_sw_NN_cloud         = 0.09

!   corr_wrice_NL_cloud      = 0.64
!   corr_wrice_NL_below      = corr_wrice_NL_cloud
!   corr_wNi_NL_cloud        = -0.35
!   corr_wNi_NL_below        = corr_wNi_NL_cloud
!   corr_wrsnow_NL_cloud     = 0.45
!   corr_wrsnow_NL_below     = corr_wrsnow_NL_cloud
!   corr_wNsnow_NL_cloud     = 0.58
!   corr_wNsnow_NL_below     = corr_wNsnow_NL_cloud
!   corr_wNc_NL_cloud        = 0.38
!   corr_wNc_NL_below        = 0.38


    ! MPACE-B values for the correlation of ice/snow
!   corr_srsnow_NL_cloud     = 0.43
!   corr_srsnow_NL_below     = corr_srsnow_NL_cloud
!   corr_sNsnow_NL_cloud     = 0.55
!   corr_sNsnow_NL_below     = corr_sNsnow_NL_cloud
!   corr_rsnowNsnow_LL_cloud = 0.91
!   corr_rsnowNsnow_LL_below = corr_rsnowNsnow_LL_cloud
!   corr_srice_NL_cloud      = 0.61
!   corr_srice_NL_below      = corr_srice_NL_cloud
!   corr_sNi_NL_cloud        = 0.90
!   corr_sNi_NL_below        = corr_sNi_NL_cloud
!   corr_riceNi_LL_cloud     = 0.63
!   corr_riceNi_LL_below     = corr_riceNi_LL_cloud
!   corr_sNc_NL_cloud        = 0.94
!   corr_sNc_NL_below        = 0.94

!   corr_sw_NN_cloud         = 0.20

!   corr_wrice_NL_cloud      = 0.01
!   corr_wrice_NL_below      = corr_wrice_NL_cloud
!   corr_wNi_NL_cloud        = 0.25
!   corr_wNi_NL_below        = corr_wNi_NL_cloud
!   corr_wrsnow_NL_cloud     = 0.01
!   corr_wrsnow_NL_below     = corr_wrsnow_NL_cloud
!   corr_wNsnow_NL_cloud     = 0.02
!   corr_wNsnow_NL_below     = corr_wNsnow_NL_cloud
!   corr_wNc_NL_cloud        = 0.24
!   corr_wNc_NL_below        = corr_wNc_NL_cloud

! NOTE: These values are commented out because they have been
!       moved to *_corr_array_cloud.in and *_corr_array_below.in
!       files. - Kenneth Connor, August 22, 2011

    ! ISDAC values for the correlation of ice/snow
!    corr_srsnow_NL_cloud     = 0.06
!    corr_srsnow_NL_below     = corr_srsnow_NL_cloud
!    corr_sNsnow_NL_cloud     = 0.04
!    corr_sNsnow_NL_below     = corr_sNsnow_NL_cloud
!    corr_rsnowNsnow_LL_cloud = 0.95
!    corr_rsnowNsnow_LL_below = corr_rsnowNsnow_LL_cloud
!    corr_srice_NL_cloud      = -0.08
!    corr_srice_NL_below      = corr_srice_NL_cloud
!    corr_sNi_NL_cloud        = 0.28
!    corr_sNi_NL_below        = corr_sNi_NL_cloud
!    corr_riceNi_LL_cloud     = 0.77
!    corr_riceNi_LL_below     = corr_riceNi_LL_cloud
!   corr_sNc_NL_cloud        = 0.09
!   corr_sNc_NL_below        = corr_sNc_NL_cloud

!    corr_sw_NN_cloud         = 0.09

!    corr_wrice_NL_cloud      = 0.44
!    corr_wrice_NL_below      = corr_wrice_NL_cloud
!    corr_wNi_NL_cloud        = 0.55
!    corr_wNi_NL_below        = corr_wNi_NL_cloud
!    corr_wrsnow_NL_cloud     = 0.65
!    corr_wrsnow_NL_below     = corr_wrsnow_NL_cloud
!    corr_wNsnow_NL_cloud     = 0.73
!    corr_wNsnow_NL_below     = corr_wNsnow_NL_cloud
!    corr_wNc_NL_cloud        = 0.34
!    corr_wNc_NL_below        = corr_wNc_NL_cloud

    !---------------------------------------------------------------------------
    ! Parameters for Morrison and COAMPS microphysics
    !---------------------------------------------------------------------------
    l_ice_micro = .true. ! Enable non-sedimenting ice and snow
    l_graupel = .true.   ! Enable graupel formation

    !---------------------------------------------------------------------------
    ! Parameters for Khairoutdinov & Kogan and COAMPS microphysics
    !---------------------------------------------------------------------------
    ! Enable to use an upwind differencing approximation for sedimentation 
    ! rather than a 3 point difference approximation.
    l_upwind_diff_sed = .false.

    !---------------------------------------------------------------------------
    ! Parameters for Morrison microphysics only
    !---------------------------------------------------------------------------
    l_hail = .false. ! Graupel will have properties of hail when true

    ! Enable to Use Seifert and Beheng (2001) warm rain parameterization 
    ! rather than Khairoutdinov Kogan (2000)
    l_seifert_beheng = .false. 
    l_predictnc = .true. ! Predict cloud droplet number concentration
    specify_aerosol = "morrison_lognormal"
    l_subgrid_w = .true.
    l_arctic_nucl = .false.
    l_fix_pgam  = .false.
    l_cloud_edge_activation = .true.

    ! Aerosol for RF02 from Mikhail Ovtchinnikov
    aer_rm1  = 0.011 ! Mean geometric radius                 [μ]
    aer_rm2  = 0.06

    aer_sig1 = 1.2   ! Std dev of aerosol size distribution  [-]
    aer_sig2 = 1.7

    aer_n1   = 125.  ! Aerosol contentration                 [#/cm3]
    aer_n2   = 65.

    ! Other parameters, set as in SAM
    ccnconst = 120. ! Parameter for powerlaw CCN [#/cm3]
    ccnexpnt = 0.4

    pgam_fixed = 5.

    !---------------------------------------------------------------------------
    ! Parameters for Morrison microphysics and Khairoutdinov & Kogan microphysics
    !---------------------------------------------------------------------------
    Ncm_initial = 100. ! #/cm^3

    ! We fix the correlations between s and t to avoid factorizing the
    ! correlation matrix more than once per simulation. 
    ! -dschanen 24 Jan 2010
    l_fix_s_t_correlations = .true. 

    l_lh_cloud_weighted_sampling = .false.
    l_lh_vert_overlap = .false.
    LH_microphys_calls = 2
    LH_sequence_length = 1

    LH_microphys_type = "disabled"
    !---------------------------------------------------------------------------
    ! Parameters for all microphysics schemes
    !---------------------------------------------------------------------------
    microphys_start_time = 0.0_time_precision ! [s]

    l_in_cloud_Nc_diff = .false. ! Don't use in cloud values of Nc for diffusion

    ! The next three lines open the cases model.in file and replace values of
    ! the parameters if they exist in the file.
    open(unit=iunit, file=namelist_file, status='old', action='read')
    read(iunit, nml=microphysics_setting)
    close(unit=iunit)

    ! Printing Microphysics inputs
    if ( clubb_at_least_debug_level( 1 ) ) then

      ! This will open the cases setup.txt file and append it to include the
      ! parameters in the microphysics_setting namelist. This file was created
      ! and written to from clubb_driver previously.
      if ( l_write_to_file ) open(unit=iunit, file=case_info_file, &
          status='old', action='write', position='append')

      ! Write to file all parameters from the namelist microphysics_seting.
      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )
      call write_text( "&microphysics_setting", l_write_to_file, iunit )
      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )

      call write_text ( "micro_scheme = " //  micro_scheme, l_write_to_file, &
        iunit )
      call write_text ( "l_cloud_sed = ", l_cloud_sed, l_write_to_file, iunit )
      call write_text ( "l_graupel = ", l_graupel, l_write_to_file, iunit )
      call write_text ( "l_hail = ", l_hail, l_write_to_file, iunit )
      call write_text ( "l_seifert_beheng = ", l_seifert_beheng, &
        l_write_to_file, iunit )
      call write_text ( "l_predictnc = ", l_predictnc, l_write_to_file, iunit )
      call write_text ( "specify_aerosol = "// specify_aerosol, &
        l_write_to_file, iunit )
      call write_text ( "l_subgrid_w = ", l_subgrid_w, l_write_to_file, iunit )
      call write_text ( "l_arctic_nucl = ", l_arctic_nucl, l_write_to_file, &
        iunit )
      call write_text ( "l_cloud_edge_activation = ", l_cloud_edge_activation, &
        l_write_to_file, iunit )
      call write_text ( "l_fix_pgam = ", l_fix_pgam, l_write_to_file, iunit )
      call write_text ( "l_in_cloud_Nc_diff = ", l_in_cloud_Nc_diff, &
        l_write_to_file, iunit )
      call write_text ( "l_upwind_diff_sed = ", l_upwind_diff_sed, &
        l_write_to_file, iunit )
      call write_text ( "LH_microphys_type = " // &
        trim( LH_microphys_type ), l_write_to_file, iunit )
      call write_text ( "LH_microphys_calls = ", LH_microphys_calls, &
        l_write_to_file, iunit )
      call write_text ( "LH_sequence_length = ", LH_sequence_length, &
        l_write_to_file, iunit )
      call write_text ( "l_lh_cloud_weighted_sampling = ", &
        l_lh_cloud_weighted_sampling, l_write_to_file, iunit )
      call write_text ( "l_fix_s_t_correlations = ", l_fix_s_t_correlations, &
        l_write_to_file, iunit )
      call write_text ( "l_lh_vert_overlap = ", l_lh_vert_overlap, &
        l_write_to_file, iunit )
      call write_text ( "rrp2_on_rrainm2_cloud = ", rrp2_on_rrainm2_cloud, &
        l_write_to_file, iunit )
      call write_text ( "Nrp2_on_Nrm2_cloud = ", Nrp2_on_Nrm2_cloud, &
        l_write_to_file, iunit )
      call write_text ( "Ncp2_on_Ncm2_cloud = ", Ncp2_on_Ncm2_cloud, &
        l_write_to_file, iunit )
      call write_text ( "corr_rrNr_LL_cloud = ", corr_rrNr_LL_cloud, &
        l_write_to_file, iunit )
      call write_text ( "corr_rtrr_NL_cloud = ", corr_rtrr_NL_cloud, &
        l_write_to_file, iunit )
      call write_text ( "corr_rtNr_NL_cloud = ", corr_rtNr_NL_cloud, &
        l_write_to_file, iunit )
      call write_text ( "corr_rtNc_NL_cloud = ", corr_rtNc_NL_cloud, &
        l_write_to_file, iunit )
      call write_text ( "corr_thlrr_NL_cloud = ", corr_thlrr_NL_cloud, &
        l_write_to_file, iunit )
      call write_text ( "corr_thlNr_NL_cloud = ", corr_thlNr_NL_cloud, &
        l_write_to_file, iunit )
      call write_text ( "corr_thlNc_NL_cloud = ", corr_thlNc_NL_cloud, &
        l_write_to_file, iunit )
      call write_text ( "corr_srr_NL_cloud = ", corr_srr_NL_cloud, &
        l_write_to_file, iunit )
      call write_text ( "corr_sNr_NL_cloud = ", corr_sNr_NL_cloud, &
        l_write_to_file, iunit )
      call write_text ( "corr_sNc_NL_cloud = ", corr_sNc_NL_cloud, &
        l_write_to_file, iunit )
      call write_text ( "rrp2_on_rrainm2_below = ", rrp2_on_rrainm2_below, &
        l_write_to_file, iunit )
      call write_text ( "Nrp2_on_Nrm2_below = ", Nrp2_on_Nrm2_below, &
        l_write_to_file, iunit )
      call write_text ( "Ncp2_on_Ncm2_below = ", Ncp2_on_Ncm2_below, &
        l_write_to_file, iunit )
      call write_text ( "corr_rrNr_LL_below = ", corr_rrNr_LL_below, &
        l_write_to_file, iunit )
      call write_text ( "corr_rtrr_NL_below = ", corr_rtrr_NL_below, &
        l_write_to_file, iunit )
      call write_text ( "corr_rtNr_NL_below = ", corr_rtNr_NL_below, &
        l_write_to_file, iunit )
      call write_text ( "corr_rtNc_NL_below = ", corr_rtNc_NL_below, &
        l_write_to_file, iunit )
      call write_text ( "corr_thlrr_NL_below = ", corr_thlrr_NL_below, &
        l_write_to_file, iunit )
      call write_text ( "corr_thlNr_NL_below = ", corr_thlNr_NL_below, &
        l_write_to_file, iunit )
      call write_text ( "corr_thlNc_NL_below = ", corr_thlNc_NL_below, &
        l_write_to_file, iunit )
      call write_text ( "corr_srr_NL_below = ", corr_srr_NL_below, &
        l_write_to_file, iunit )
      call write_text ( "corr_sNr_NL_below = ", corr_sNr_NL_below, &
        l_write_to_file, iunit )
      call write_text ( "corr_sNc_NL_below = ", corr_sNc_NL_below, &
        l_write_to_file, iunit )
      call write_text ( "C_evap = ", C_evap, l_write_to_file, iunit )
      call write_text ( "r_0 = ", r_0, l_write_to_file, iunit )
      call write_text ( "microphys_start_time = ", real( microphys_start_time ), &
        l_write_to_file, iunit )
      call write_text ( "Ncm_initial = ", Ncm_initial, l_write_to_file, iunit )
      call write_text ( "ccnconst = ", ccnconst, l_write_to_file, iunit )
      call write_text ( "ccnexpnt = ", ccnexpnt, l_write_to_file, iunit )
      call write_text ( "aer_rm1 = ", aer_rm1, l_write_to_file, iunit )
      call write_text ( "aer_rm2 = ", aer_rm2, l_write_to_file, iunit )
      call write_text ( "aer_n1 = ", aer_n1, l_write_to_file, iunit )
      call write_text ( "aer_n2 = ", aer_n2, l_write_to_file, iunit )
      call write_text ( "aer_sig1 = ", aer_sig1, l_write_to_file, iunit )
      call write_text ( "aer_sig2 = ", aer_sig2, l_write_to_file, iunit )
      call write_text ( "pgam_fixed = ", pgam_fixed, l_write_to_file, iunit )

      if ( trim( micro_scheme ) == "khairoutdinov_kogan" .and. l_predictnc ) then
        call write_text("WARNING: l_predictnc has been set to .false. for KK microphysics.", &
          l_write_to_file, iunit )
      end if

      if ( l_write_to_file ) close(unit=iunit);

    end if ! clubb_at_least_debug_level(1)

    ! Change l_predictnc if we are using KK microphysics
    ! This flag is typically used only in Morrison microphysics, but has been adapted
    ! for use with KK microphysics.  It is used to turn off the budget terms for Ncm
    ! since Ncm is essentially clipped every timestep when using KK microphysics.
    if ( trim( micro_scheme ) == "khairoutdinov_kogan" ) then
      l_predictnc = .false.
    end if

    ! Read in the name list for initialization, if it exists
    open(unit=iunit, file=namelist_file, status='old', action='read')
    read(iunit, nml=gfdl_activation_setting)
    close(unit=iunit)

    ! Initialize the GFDL activation code, if necessary
    if( l_gfdl_activation ) then
      ! Ensure a microphysics that has Ncm is being used
      if( trim( micro_scheme ) == "coamps" .or. trim( micro_scheme ) == "morrison" & 
            .or. trim( micro_scheme ) == "morrison-gettelman") then

          ! Read in the lookup tables
          call Loading( droplets, droplets2 )

          ! Initialize the activation variables
          call aer_ccn_act_k_init                            &
                 ( droplets, droplets2, res, res2, nooc,     &
                   sul_concen, low_concen, high_concen,      &
                   lowup, highup, lowup2, highup2, lowmass2, &
                   highmass2, lowmass3, highmass3,           &
                   lowmass4, highmass4, lowmass5, highmass5, &
                   lowT2, highT2 )
      end if ! coamps .or. morrison .or. khairoutdinov_kogan
    end if ! l_gfdl_activation

    ! The location of the fields in the hydromet array are arbitrary,
    ! and don't need to be set consistently among schemes so long as
    ! the 'i' indices point to the correct parts of the array.

    select case ( trim( micro_scheme ) )
    case ( "morrison" )
      iirrainm    = 1
      iirsnowm    = 2
      iiricem     = 3
      iirgraupelm = 4

      iiNrm       = 5
      iiNsnowm    = 6
      iiNim       = 7
      iiNgraupelm = 8
      iiNcm       = 9

      hydromet_dim = 9

      allocate( hydromet_list(hydromet_dim) )

      hydromet_list(iirrainm)    = "rrainm"
      hydromet_list(iirsnowm)    = "rsnowm"
      hydromet_list(iiricem)     = "ricem"
      hydromet_list(iirgraupelm) = "rgraupelm"

      hydromet_list(iiNrm)       = "Nrm"
      hydromet_list(iiNsnowm)    = "Nsnowm"
      hydromet_list(iiNim)       = "Nim"
      hydromet_list(iiNgraupelm) = "Ngraupelm"
      hydromet_list(iiNcm)       = "Ncm"

      ! Set Nc0 in the Morrison code (module_MP_graupel) based on Ncm_initial
      Nc0 = Ncm_initial

      ! Set flags from the Morrison scheme as in GRAUPEL_INIT
      if ( l_predictnc ) then
        dopredictNc = .true.
      else
        dopredictNc = .false.
      end if

      ! Set the mode of aerosol to be used
      if ( trim( specify_aerosol ) == "morrison_no_aerosol" ) then
        aerosol_mode = morrison_no_aerosol
      else if ( trim( specify_aerosol ) == "morrison_power_law" ) then
        aerosol_mode = morrison_power_law
      else if ( trim( specify_aerosol ) == "morrison_lognormal" ) then
        aerosol_mode = morrison_lognormal
      else
        stop "Unknown Morrison aerosol mode." 
      end if

      if ( l_cloud_edge_activation ) then
        docloudedgeactivation = .true.
      else
        docloudedgeactivation = .false.
      end if

      if ( l_ice_micro ) then
        doicemicro = .true.
      else
        doicemicro = .false.
      end if

      if ( l_arctic_nucl ) then
        doarcticicenucl = .true.
      else
        doarcticicenucl = .false.
      end if

      if ( l_graupel ) then
        dograupel = .true.
      else
        dograupel = .false.
      end if

      if ( l_hail ) then
        dohail = .true.
      else
        dohail = .false.
      end if

      if ( l_seifert_beheng ) then
        dosb_warm_rain = .true.
      else
        dosb_warm_rain = .false.
      end if

      if ( l_fix_pgam ) then
        dofix_pgam = .true.
      else
        dofix_pgam = .false.
      end if

      if ( l_subgrid_w ) then
        dosubgridw = .true.
      else
        dosubgridw = .false.
      end if

      if ( l_cloud_sed ) then
        write(fstderr,*) "Morrison microphysics has seperate code for cloud water sedimentation,"
        write(fstderr,*) " therefore l_cloud_sed should be set to .false."
        stop "Fatal error."
      end if

      if ( .not. l_fix_s_t_correlations .and. l_ice_micro &
           .and. trim( LH_microphys_type ) /= "disabled" ) then
        write(fstderr,*) "The variable l_fix_s_t_correlations must be true in order to "// &
          "enable latin hypercube sampling and ice microphysics."
        stop "Fatal error."
      end if
      allocate( l_hydromet_sed(hydromet_dim) )
      ! Sedimentation is handled within the Morrison microphysics
      l_hydromet_sed(iiNrm) = .false.
      l_hydromet_sed(iiNim) = .false.
      l_hydromet_sed(iiNcm) = .false.
      l_hydromet_sed(iiNgraupelm) = .false.
      l_hydromet_sed(iiNsnowm) = .false.

      l_hydromet_sed(iirrainm)    = .false.
      l_hydromet_sed(iirsnowm)    = .false.
      l_hydromet_sed(iiricem)     = .false.
      l_hydromet_sed(iirgraupelm) = .false.

      ! Convert from μ to m as in SAM
      aer_rm1 = aer_rm1 / micron_per_m
      aer_rm2 = aer_rm2 / micron_per_m
      ! Convert from #/cm3 to #/m3
      aer_n1 = cm3_per_m3 * aer_n1
      aer_n2 = cm3_per_m3 * aer_n2

      ! Setup the Morrison scheme
      call GRAUPEL_INIT()
      
    case ( "morrison-gettelman" )
      iirrainm    = -1
      iirsnowm    = -1
      iiricem     = 1
      iirgraupelm = -1

      iiNrm       = -1
      iiNsnowm    = -1
      iiNim       = 2
      iiNgraupelm = -1
      iiNcm       = 3

      hydromet_dim = 3
      
      allocate( hydromet_list(hydromet_dim) )
      
      hydromet_list(iiricem)     = "ricem"
      hydromet_list(iiNim)       = "Nim"
      hydromet_list(iiNcm)       = "Ncm"
      
      allocate( l_hydromet_sed(hydromet_dim) )
      ! Sedimentation is handled within the MG microphysics
      l_hydromet_sed(iiricem)  = .false.
      l_hydromet_sed(iiNim)    = .false.
      l_hydromet_sed(iiNcm)    = .false.
      
      ! Initialize constants for aerosols
      call ini_microp_aero()

      ! Setup the MG scheme
      call ini_micro()
      call pbuf_init()
      call pbuf_add('WP2', 1, gr%nzmax, 1)
      call pbuf_allocate()
      
    case ( "coamps" )
      iirrainm    = 1
      iirsnowm    = 2
      iiricem     = 3
      iirgraupelm = 4

      iiNrm       = 5
      ! Nsnowm is computed diagnostically in the subroutine coamps_micro_driver
      iiNsnowm    = -1
      iiNim       = 6
      iiNgraupelm = -1
      iiNcm       = 7

      hydromet_dim = 7

      allocate( l_hydromet_sed(hydromet_dim) )

      allocate( hydromet_list(hydromet_dim) )

      hydromet_list(iirrainm)    = "rrainm"
      hydromet_list(iirsnowm)    = "rsnowm"
      hydromet_list(iiricem)     = "ricem"
      hydromet_list(iirgraupelm) = "rgraupelm"

      hydromet_list(iiNrm)       = "Nrm"
      hydromet_list(iiNcm)       = "Ncm"
      hydromet_list(iiNim)       = "Nim"

      l_hydromet_sed(iiNrm) = .true.
      l_hydromet_sed(iiNim) = .false.
      l_hydromet_sed(iiNcm) = .false.

      l_hydromet_sed(iirrainm)    = .true.
      l_hydromet_sed(iirsnowm)    = .true.
      l_hydromet_sed(iiricem)     = .true.
      l_hydromet_sed(iirgraupelm) = .true.

    case ( "khairoutdinov_kogan" )
      iirrainm    = 1
      iirsnowm    = -1
      iiricem     = -1
      iirgraupelm = -1

      iiNrm       = 2
      iiNsnowm    = -1
      iiNim       = -1
      iiNgraupelm = -1
      iiNcm       = 3

      hydromet_dim = 3

      allocate( l_hydromet_sed(hydromet_dim) )

      allocate( hydromet_list(hydromet_dim) )

      hydromet_list(iirrainm) = "rrainm"
      hydromet_list(iiNrm) = "Nrm"
      hydromet_list(iiNcm) = "Ncm"

      l_hydromet_sed(iirrainm) = .true.
      l_hydromet_sed(iiNrm)    = .true.
      l_hydromet_sed(iiNcm)    = .false.

    case ( "simplified_ice" )
      iirrainm    = -1
      iirsnowm    = -1
      iiricem     = -1
      iirgraupelm = -1

      iiNrm       = -1
      iiNsnowm    = -1
      iiNim       = -1
      iiNgraupelm = -1
      iiNcm       = -1

      hydromet_dim = 0

    case ( "none" )
      if ( l_cloud_sed ) then
        iiNcm       = 1
        hydromet_dim = 1

        allocate( l_hydromet_sed(hydromet_dim) )
        allocate( hydromet_list(hydromet_dim) )
        hydromet_list(iiNcm) = "Ncm"
      else
        iiNcm = -1
        hydromet_dim = 0
      end if

      iirrainm    = -1
      iirsnowm    = -1
      iiricem     = -1
      iirgraupelm = -1

      iiNrm       = -1
      iiNsnowm    = -1
      iiNim       = -1
      iiNgraupelm = -1

    case default
      write(fstderr,*) "Unknown micro_scheme"// trim( micro_scheme )
      stop

    end select

    ! Sanity check
    if ( l_lh_cloud_weighted_sampling .and. .not.  l_lh_vert_overlap ) then
      write(fstderr,*) "Error in init_microphys: "// &
        "l_lh_cloud_weighted_sampling requires l_lh_vert_overlap."
      stop
    end if

    select case ( trim( LH_microphys_type ) )
    case ( "interactive" )
      LH_microphys_type_int = LH_microphys_interactive

    case ( "non-interactive" )
      LH_microphys_type_int = LH_microphys_non_interactive

    case ( "disabled" )
      LH_microphys_type_int = LH_microphys_disabled

    case default
      stop "Error determining LH_microphys_type"

    end select

    ! Setup index variables for latin hypercube sampling
    if ( LH_microphys_type_int /= LH_microphys_disabled ) then

#ifdef LATIN_HYPERCUBE
      ! Allocate and set the arrays containing the correlations
      ! and the X'^2 / X'^2 terms
      call setup_corr_varnce_array( iiNcm, iirrainm, iiNrm, iiricem, iiNim, iirsnowm, iiNsnowm, &
                                    LH_input_path, iunit, runtype )

#endif

    end if

    return
  end subroutine init_microphys

!-------------------------------------------------------------------------------
  subroutine advance_microphys & 
             ( iter, runtype, dt, time_current,  & 
               thlm, p_in_Pa, exner, rho, rho_zm, rtm, rcm, cloud_frac, & 
               wm_zt, wm_zm, Kh_zm, pdf_params, & 
               wp2_zt, rho_ds_zt, rho_ds_zm, LH_sample_point_weights, &
               X_nl_all_levs, X_mixt_comp_all_levs, LH_rt, LH_thl, &
               Ncnm, hydromet, & 
               rvm_mc, rcm_mc, thlm_mc, &
               err_code )

    ! Description:
    ! Compute pristine ice, snow, graupel, & rain hydrometeor fields.
    ! Uses implicit discretization.

    ! References:
    !   None
    !---------------------------------------------------------------------------

    use grid_class, only: & 
        gr,  & ! Variable(s)
        zm2zt,  & ! Procedure(s)
        zt2zm

    use KK_microphys_module, only: & 
        KK_micro_driver  ! Procedure(s)

    use morrison_micro_driver_module, only: &
        morrison_micro_driver
      
    use mg_micro_driver_module, only: &
      mg_microphys_driver

#ifdef LATIN_HYPERCUBE
    use latin_hypercube_driver_module, only: &
      LH_microphys_driver ! Procedure

!   use latin_hypercube_arrays, only: &
!     X_nl_all_levs, & ! Variable(s)
!     X_mixt_comp_all_levs, &
!     LH_rt, LH_thl

#endif /*LATIN_HYPERCUBE*/

    use ice_dfsn_module, only: & 
        ice_dfsn ! Procedure(s)

    use T_in_K_module, only: &
        thlm2T_in_K ! Procedure

    use gfdl_activation, only: &
        aer_act_clubb_quadrature_Gauss, & ! Procedure
        aeromass_value                    ! Variable

    use parameters_tunable, only: & 
        c_Krrainm,  & ! Variable(s) 
        nu_r_vert_res_dep

    use parameters_model, only: & 
        hydromet_dim   ! Integer

    use constants_clubb, only: & 
        Lv, & ! Constant(s)
        Ls, &
        Cp, & 
        rho_lw, & 
        rc_tol, &
        fstderr, & 
        zero_threshold, &
        sec_per_day, &
        cm3_per_m3, &
        mm_per_m

    use model_flags, only: &
      l_hole_fill ! Variable(s)

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use error_code, only:  & 
        fatal_error, & ! Procedure
        clubb_at_least_debug_level, &
        clubb_no_error  ! Constant

#ifdef COAMPS_MICRO
    use coamps_micro_driver_module, only:  & 
        coamps_micro_driver ! Procedure
#endif

    use pdf_parameter_module, only:  &
        pdf_parameter  ! Type

    use array_index, only:  & 
        iirrainm, iirsnowm, iiricem, iirgraupelm, & ! Variable(s)
        iiNrm, iiNim, iiNcm

    use stats_variables, only: & 
      iVrr,  & ! Variable(s)
      iVNr, & 
      iVsnow, & 
      iVice, & 
      iVgraupel, & 
      irain_rate_zt, & 
      iFprec, & 
      irrainm_bt, & 
      irrainm_mc, & 
      irrainm_cond_adj, & 
      irrainm_cl, & 
      iNrm_bt, & 
      iNrm_mc, & 
      iNrm_cond_adj, & 
      iNrm_cl, &
      iNcm_in_cloud, &
      iNc_activated, &
      iNcnm, & 
      iricem_bt, & 
      iricem_mc, & 
      iricem_cl, & 
      irgraupelm_bt, & 
      irgraupelm_mc, & 
      irgraupelm_cl, & 
      irsnowm_bt, & 
      irsnowm_mc, & 
      irsnowm_cl, & 
      irain_rate_sfc, & 
      irain_flux_sfc, & 
      irrainm_sfc

    use stats_variables, only: & 
      iNim_bt, &
      iNim_cl, &
      iNim_mc, &
      iNsnowm_bt, &
      iNsnowm_cl, &
      iNsnowm_mc, &
      iNgraupelm_bt, &
      iNgraupelm_cl, &
      iNgraupelm_mc, &
      iNcm_bt, & 
      iNcm_cl, &
      iNcm_mc, &
      iNcm_act

    use stats_subs, only: & 
     stats_accumulate_LH_tend ! Procedure(s)

    use stats_variables, only: & 
      iLH_Vrr, &
      iLH_VNr

    use stats_variables, only: & 
        zt, &  ! Variables
        zm, & 
        sfc, & 
        l_stats_samp

    use stats_type, only: & 
        stat_update_var, stat_update_var_pt, & ! Procedure(s)
        stat_begin_update, stat_end_update

    use stats_subs, only: &
        stats_accumulate_hydromet

    use fill_holes, only: &
      vertical_avg, & ! Procedure(s)
      fill_holes_driver

    use parameters_microphys, only: &
      LH_microphys_type, & ! Determines how the LH samples are used
      LH_microphys_interactive,     & ! Feed the subcolumns into the microphysics and allow feedback
      LH_microphys_non_interactive, & ! Feed the subcolumns into the microphysics with no feedback
      LH_microphys_disabled           ! Disable latin hypercube entirely

#ifdef LATIN_HYPERCUBE
    use latin_hypercube_arrays, only: &
      d_variables ! Variable(s)
#else
#define d_variables 0
#endif /* LATIN_HYPERCUBE */

    implicit none

    ! Input Variables

    integer, intent(in) :: iter ! Model iteration number

    character(len=*), intent(in) :: & 
      runtype ! Name of the run, for case specific effects.

    real(kind=time_precision), intent(in) ::  & 
      dt           ! Timestep         [s]

    real(kind=time_precision), intent(in) ::  & 
      time_current ! Current time     [s]

    real, dimension(gr%nzmax), intent(in) :: & 
      thlm,       & ! Liquid potential temp.                 [K]
      p_in_Pa,    & ! Pressure                               [Pa]
      exner,      & ! Exner function                         [-]
      rho,        & ! Density on thermo. grid                [kg/m^3]
      rho_zm,     & ! Density on moment. grid                [kg/m^3]
      rtm,        & ! Total water mixing ratio               [kg/kg]
      rcm,        & ! Liquid water mixing ratio              [kg/kg]
      cloud_frac, & ! Cloud fraction                         [-]
      wm_zt,      & ! w wind on moment. grid                 [m/s]
      wm_zm,      & ! w wind on thermo. grid                 [m/s]
      Kh_zm         ! Kh Eddy diffusivity on momentum grid   [m^2/s]

    type(pdf_parameter), dimension(gr%nzmax), intent(in) :: & 
      pdf_params     ! PDF parameters

    real, dimension(gr%nzmax), intent(in) :: & 
      wp2_zt,    & ! w'^2 on the thermo. grid               [m^2/s^2]
      rho_ds_zm, & ! Dry, static density on moment. levels  [kg/m^3]
      rho_ds_zt    ! Dry, static density on thermo. levels  [kg/m^3]

    double precision, dimension(gr%nzmax,LH_microphys_calls,d_variables), intent(in) :: &
      X_nl_all_levs ! Lognormally distributed hydrometeors

    integer, dimension(gr%nzmax,LH_microphys_calls), intent(in) :: &
      X_mixt_comp_all_levs ! Which mixture component the sample is in

    real, dimension(gr%nzmax,LH_microphys_calls), intent(in) :: &
      LH_rt, LH_thl ! Samples of rt, thl	[kg/kg,K]

    real, dimension(LH_microphys_calls), intent(in) :: &
      LH_sample_point_weights ! Weights for cloud weighted sampling

    ! Note:
    ! K & K only uses Ncm, while for COAMPS Ncnm is initialized
    ! and Nim & Ncm are computed within subroutine adjtg.
    real, dimension(gr%nzmax), intent(inout) :: & 
      Ncnm       ! Cloud nuclei number concentration     [count/m^3]

    real, dimension(gr%nzmax,hydromet_dim), intent(inout) :: & 
      hydromet      ! Array of rain, prist. ice, graupel, etc. [units vary]

    real, dimension(gr%nzmax), intent(inout) :: & 
      rcm_mc,  & ! Microphysics contributions to liquid water           [kg/kg/s]
      rvm_mc,  & ! Microphysics contributions to vapor water            [kg/kg/s]
      thlm_mc    ! Microphysics contributions to liquid potential temp. [K/s]

    integer, intent(out) :: err_code ! Exit code returned from subroutine

    ! Local Variables
    real, dimension(3,gr%nzmax) :: & 
      lhs ! Left hand side of tridiagonal matrix

    real, dimension(gr%nzmax,hydromet_dim) :: &
      hydromet_vel,    & ! Contains vel. of the hydrometeors        [m/s]
      hydromet_vel_zt

    real, dimension(gr%nzmax,hydromet_dim) :: & 
      hydromet_mc  ! Change in hydrometeors due to microphysics  [units/s]

    real, dimension(gr%nzmax) :: &
      delta_zt  ! Difference in thermo. height levels     [m]

   real, dimension(gr%nzmax) :: &
     T_in_K  ! Temperature          [K]

    real, dimension(1,1,gr%nzmax) :: & 
      cond ! COAMPS stat for condesation/evap of rcm

    ! Eddy diffusivity for rain and rain drop concentration.
    ! It is also used for the other hydrometeor variables.
    ! Kr = Constant * Kh_zm; Constant is named c_Krrainm.
    real, dimension(gr%nzmax) :: Kr   ! [m^2/s]

    ! Variable needed to handle correction to rtm and thlm microphysics
    ! tendency arrays, as well rrainm_cond and Nrm_cond statistical
    ! tendency arrays, due to a negative result being produced by
    ! over-evaporation of rain water over the course of a timestep.
    ! Brian Griffin.  April 14, 2007.
    real :: overevap_rate ! Absolute value of negative evap. rate.

    real, dimension(gr%nzmax) :: &
      wtmp,    & ! Standard dev. of w                   [m/s]
      s_mellor   ! The variable 's' in Mellor (1977)    [kg/kg]

    real :: &
      max_velocity, & ! Maximum sedimentation velocity [m/s]
      temp            ! Temporary variables     [units vary]

    integer :: i, k ! Loop iterators / Array indices

    integer :: ixrm_cl, ixrm_bt, ixrm_mc

    real, dimension(gr%nzmax) :: Ndrop_max  ! GFDL droplet activation concentration [#/kg]

    !Input aerosol mass concentration: the unit is 10^12 ug/m3.
    !For example, aeromass=2.25e-12 means that the aerosol mass concentration is 2.25 ug/m3.
    !This value of aeromass was recommended by Huan Guo
    !See http://carson.math.uwm.edu/trac/climate_process_team/ticket/46#comment:12
    real, dimension(gr%nzmax, 4) :: aeromass ! ug/m^3

    real, dimension(gr%nzmax) :: Ncm_in_cloud ! cloud average droplet concentration [#/kg]

    character(len=10) :: hydromet_name

    logical :: l_local_kk_input, l_latin_hypercube_input, l_sed

!-------------------------------------------------------------------------------

    ! ---- Begin code ----

    Ndrop_max = 0.
    ! Set the value of the aerosol mass array to a constant
    aeromass = aeromass_value

    err_code = clubb_no_error  ! Initialize to the value for no errors

    ! Make some compiler warnings go away for external users
    if ( runtype == "" .or. iter == -1 ) then
      stop "Runtype is null or iter is -1, neither of which should happen."
    end if

    ! Return if there is delay between the model start time and start of the
    ! microphysics
    if ( time_current < microphys_start_time ) return

    ! Solve for the value of Kr, the hydrometeor eddy diffusivity.
    do k = 1, gr%nzmax, 1
      Kr(k) = c_Krrainm * Kh_zm(k)
    end do

    ! Determine 's' from Mellor (1977)
    s_mellor(:) = pdf_params(:)%mixt_frac * pdf_params(:)%s1 &
                + (1.0-pdf_params(:)%mixt_frac) * pdf_params(:)%s2

    ! Compute standard deviation of vertical velocity in the grid column
    wtmp(:) = sqrt( wp2_zt(:) )

    ! Compute difference in thermodynamic height levels
    delta_zt(1:gr%nzmax) = 1./gr%invrs_dzm(1:gr%nzmax)

    ! Calculate T_in_K
    T_in_K = thlm2T_in_K( thlm, exner, rcm )

    ! Start Ncm budget to capture Ncm_act term
    if( l_stats_samp ) then
      call stat_begin_update( iNcm_bt, hydromet(:, iiNcm) / real( dt ), zt )
    end if
    
    ! Call GFDL activation code
    if( l_gfdl_activation ) then
      ! Ensure a microphysics that has Ncm is being used
      ! Note: KK micro is not used because it resets Ncm every timestep 
      if( trim( micro_scheme ) == "coamps" .or. trim( micro_scheme ) == "morrison" & 
            .or. trim( micro_scheme ) == "morrison-gettelman" ) then

        ! Save the initial Ncm value for the Ncm_act term
        if( l_stats_samp ) then
          call stat_begin_update( iNcm_act, hydromet(:, iiNcm) / real( dt ), zt )
        end if

        call aer_act_clubb_quadrature_Gauss( aeromass, T_in_K, Ndrop_max )

        ! Convert to #/kg
        Ndrop_max = Ndrop_max * cm3_per_m3 / rho

        if( l_stats_samp ) then
          call stat_update_var( iNc_activated, Ndrop_max, zt)
        end if

        ! Clip Ncm values that are outside of cloud by CLUBB standards
        do k=1, gr%nzmax

! ---> h1g, 2011-04-20,   no liquid drop nucleation if T < -40 C
          if( T_in_K(k) <= 233.15 )  Ndrop_max(k) = 0.0  ! if T<-40C, no liquid drop nucleation
! <--- h1g, 2011-04-20

          ! Clip Ncm to only be where there is cloud to avoid divide by zero errors
          if( cloud_frac(k) > cloud_frac_min ) then
                hydromet(k, iiNcm) = max( Ndrop_max(k), hydromet(k, iiNcm) )
          else
                hydromet(k, iiNcm) = 0.0
          end if
        end do

        ! Update the Ncm_act term
        if( l_stats_samp ) then
          call stat_end_update( iNcm_act, hydromet(:,iiNcm) / real( dt ), zt )
        end if

      else

        stop "Unsupported microphysics scheme for GFDL activation."

      end if  ! coamps .or. morrison .or. khairoutdinov_kogan
    end if ! l_gfdl_activation

    ! Begin by calling Brian Griffin's implementation of the
    ! Khairoutdinov and Kogan microphysics (analytic or local formulas),
    ! the COAMPS implementation of Rutlege and Hobbes, or the Morrison
    ! microphysics.
    ! Note: COAMPS appears to have some K&K elements to it as well.

    select case ( trim( micro_scheme ) )

    case ( "coamps" )

      ! Initialize tendencies to zero
      hydromet_mc(:,:) = 0.0
      hydromet_vel(:,:) = 0.0
      hydromet_vel_zt(:,:) = 0.0
#ifdef COAMPS_MICRO
      call coamps_micro_driver & 
           ( runtype, time_current, dt, & 
             rtm, wm_zm, p_in_Pa, exner, rho, & 
             thlm, hydromet(:,iiricem), hydromet(:,iirrainm),  & 
             hydromet(:,iirgraupelm), hydromet(:,iirsnowm), & 
             rcm, hydromet(:,iiNcm), hydromet(:,iiNrm), Ncnm, hydromet(:,iiNim), cond, & 
             hydromet_vel_zt(:,iirsnowm), hydromet_vel_zt(:,iiricem), & 
             hydromet_vel_zt(:,iirrainm), hydromet_vel_zt(:,iiNrm),  & 
             hydromet_vel_zt(:,iirgraupelm), & 
             hydromet_mc(:,iiricem), hydromet_mc(:,iirrainm), & 
             hydromet_mc(:,iirgraupelm), hydromet_mc(:,iirsnowm), & 
             hydromet_mc(:,iiNrm), & 
             rvm_mc, rcm_mc, thlm_mc )
#else
      stop "Not compiled with COAMPS microphysics"
      cond = -999.
      if ( cond(1,1,1) /= cond(1,1,1) ) stop
#endif

      if ( l_stats_samp ) then

        ! Sedimentation velocity for rrainm
        call stat_update_var(iVrr, zt2zm( hydromet_vel_zt(:,iirrainm) ), zm)

        ! Sedimentation velocity for Nrm
        call stat_update_var(iVNr, zt2zm( hydromet_vel(:,iiNrm) ), zm )

        ! Sedimentation velocity for snow
        call stat_update_var(iVsnow, zt2zm( hydromet_vel(:,iirsnowm) ), zm )

        ! Sedimentation velocity for pristine ice
        call stat_update_var( iVice, zt2zm( hydromet_vel(:,iiricem) ), zm )

        ! Sedimentation velocity for graupel
        call stat_update_var( iVgraupel, zt2zm( hydromet_vel(:,iirgraupelm) ), zm )
      end if ! l_stats_samp

    case ( "morrison" )

      ! Initialize tendencies to zero
      hydromet_mc(:,:) = 0.0
      rcm_mc(:) = 0.0
      rvm_mc(:) = 0.0
      thlm_mc(:) = 0.0

      if ( LH_microphys_type /= LH_microphys_disabled ) then
#ifdef LATIN_HYPERCUBE
        call LH_microphys_driver &
             ( real( dt ), gr%nzmax, LH_microphys_calls, d_variables, & ! In
               X_nl_all_levs, LH_rt, LH_thl, LH_sample_point_weights, & ! In
               pdf_params, p_in_Pa, exner, rho, & ! In
               rcm, wtmp, delta_zt, cloud_frac, & ! In
               hydromet, X_mixt_comp_all_levs, & !In 
               hydromet_mc, hydromet_vel_zt, & ! In/Out
               rcm_mc, rvm_mc, thlm_mc,  & ! Out
               morrison_micro_driver )  ! Procedure
#else
        stop "Latin hypercube was not enabled at compile time"
        ! Get rid of compiler warnings
        if ( .false. .and. size( X_nl_all_levs ) < 1 ) then
           rcm_mc(1) = + LH_rt(1,1) + LH_thl(1,1) &
             + LH_sample_point_weights(1) + real( X_mixt_comp_all_levs(1,1) )
        end if
#endif
        call stats_accumulate_LH_tend( hydromet_mc, thlm_mc, rvm_mc, rcm_mc )

      end if ! LH isn't disabled

      ! Based on YSU PBL interface to the Morrison scheme WRF driver, the standard dev. of w
      ! will be clipped to be between 0.1 m/s and 4.0 m/s in WRF.  -dschanen 23 Mar 2009
!     wtmp(:) = max( 0.1, wtmp ) ! Disabled for now
!     wtmp(:) = min( 4., wtmp )

!     wtmp = 0.5 ! %% debug
      ! Call the microphysics if we don't want to have feedback effects from the
      ! latin hypercube result (above)
      if ( LH_microphys_type /= LH_microphys_interactive ) then
        l_local_kk_input = .false.
        l_latin_hypercube_input = .false.
        call morrison_micro_driver & 
             ( real( dt ), gr%nzmax, l_stats_samp, l_local_kk_input, l_latin_hypercube_input, &
               thlm, p_in_Pa, exner, rho, pdf_params, &
               wm_zt, wtmp, delta_zt, rcm, s_mellor, rtm-rcm, hydromet, hydromet_mc, &
               hydromet_vel_zt, rcm_mc, rvm_mc, thlm_mc )
      end if
      

    case ( "morrison-gettelman" )

      ! Initialize tendencies to zero
      hydromet_mc(:,:) = 0.0
      rcm_mc(:) = 0.0
      rvm_mc(:) = 0.0
      thlm_mc(:) = 0.0

      if ( LH_microphys_type /= LH_microphys_disabled ) then
#ifdef LATIN_HYPERCUBE
!       call LH_microphys_driver &
!            ( real( dt ), gr%nzmax, LH_microphys_calls, d_variables, & ! In
!              X_nl_all_levs, LH_rt, LH_thl, & ! In
!              pdf_params, p_in_Pa, exner, rho, & ! In
!              rcm, wtmp, delta_zt, cloud_frac, & ! In
!              hydromet, X_mixt_comp_all_levs, & !In 
!              hydromet_mc, hydromet_vel_zt, & ! In/Out
!              rcm_mc, rvm_mc, thlm_mc,  & ! Out
!              mg_microphys_driver )  ! Procedure
        stop "Latin hypercube is not setup for MG yet"
#else
        stop "Latin hypercube was not enabled at compile time"
#endif
        call stats_accumulate_LH_tend( hydromet_mc, thlm_mc, rvm_mc, rcm_mc )

      end if ! LH isn't disabled

      ! Call the microphysics if we don't want to have feedback effects from the
      ! latin hypercube result (above)
      if ( LH_microphys_type /= LH_microphys_interactive ) then
        l_local_kk_input = .false.
        l_latin_hypercube_input = .false.
        call mg_microphys_driver &
          ( real( dt ), gr%nzmax, l_stats_samp, l_local_kk_input, l_latin_hypercube_input, &
              thlm, p_in_Pa, exner, rho, pdf_params, &
              rcm, rtm-rcm, Ncnm, hydromet, hydromet_mc, &
              hydromet_vel_zt, rcm_mc, rvm_mc, thlm_mc, wp2_zt)
      end if
          
    case ( "khairoutdinov_kogan" )

      ! Initialize tendencies to zero
      hydromet_mc(:,:) = 0.0
      rcm_mc(:) = 0.0
      rvm_mc(:) = 0.0
      thlm_mc(:) = 0.0

      if ( LH_microphys_type /= LH_microphys_disabled ) then

#ifdef LATIN_HYPERCUBE
        call LH_microphys_driver &
             ( real( dt ), gr%nzmax, LH_microphys_calls, d_variables, & ! In
               X_nl_all_levs, LH_rt, LH_thl, LH_sample_point_weights, & ! In
               pdf_params, p_in_Pa, exner, rho, & ! In
               rcm, wtmp, delta_zt, cloud_frac, & ! In
               hydromet, X_mixt_comp_all_levs, & !In 
               hydromet_mc, hydromet_vel_zt, & ! In/Out
               rcm_mc, rvm_mc, thlm_mc,  & ! Out
               KK_micro_driver ) ! Procedure
#else
        stop "Latin hypercube was not enabled at compile time"
#endif

        call stats_accumulate_LH_tend( hydromet_mc, thlm_mc, rvm_mc, rcm_mc )

        if ( l_stats_samp ) then
          ! Latin hypercube estimate for sedimentation velocities
          call stat_update_var( iLH_Vrr, zt2zm( hydromet_vel_zt(:,iirrainm) ), zt )

          call stat_update_var( iLH_VNr, zt2zm( hydromet_vel_zt(:,iiNrm) ), zt )

        end if

      end if ! LH isn't disabled

      ! Call the microphysics if we don't want to have feedback effects from the
      ! latin hypercube result (above)
      if ( LH_microphys_type /= LH_microphys_interactive ) then

!         call KK_micro_driver( real( dt ), l_local_kk, thlm, rho, p_in_Pa, &
!                               exner, s_mellor, rcm, hydromet, &
!                               pdf_params, hydromet_mc, hydromet_vel_zt, &
!                               rcm_mc, rvm_mc, thlm_mc )
         l_latin_hypercube_input = .false.
         call KK_micro_driver &
                 ( real( dt ), gr%nzmax, l_stats_samp, l_local_kk, &
                   l_latin_hypercube_input, thlm, p_in_Pa, exner, rho, &
                   pdf_params, wm_zt, wtmp, delta_zt, rcm, s_mellor, &
                   rtm-rcm, hydromet, hydromet_mc, hydromet_vel_zt, &
                   rcm_mc, rvm_mc, thlm_mc )

      end if

      if ( l_stats_samp ) then
        ! Sedimentation velocity for rrainm
        call stat_update_var( iVrr, zt2zm( hydromet_vel_zt(:,iirrainm) ), zm )

        ! Sedimentation velocity for Nrm
        call stat_update_var( iVNr, zt2zm( hydromet_vel_zt(:,iiNrm) ), zm )
      end if

    case default

    end select ! micro_scheme

    !-----------------------------------------------------------------------
    !       Loop over all hydrometeor species and apply sedimentation,
    !       advection and diffusion.
    !-----------------------------------------------------------------------

    if ( hydromet_dim > 0 ) then

      do i = 1, hydromet_dim

        ! Initializing max_velocity in order to avoid a compiler warning.
        ! Regardless of the case, it will be reset in the 'select case'
        ! statement immediately below.
        max_velocity = 0.0

        select case ( trim( hydromet_list(i) ) )
        case ( "rrainm" )
          ixrm_bt = irrainm_bt
          ixrm_cl = irrainm_cl
          ixrm_mc = irrainm_mc

          max_velocity = -9.1 ! m/s

        case ( "ricem" )
          ixrm_bt = iricem_bt
          ixrm_cl = iricem_cl
          ixrm_mc = iricem_mc

          max_velocity = -1.2 ! m/s

        case ( "rsnowm" )
          ixrm_bt = irsnowm_bt
          ixrm_cl = irsnowm_cl
          ixrm_mc = irsnowm_mc

          ! Morrison limit
!         max_velocity = -1.2 ! m/s
          ! Made up limit.  The literature suggests that it is quite possible
          ! that snow flake might achieve a terminal velocity of 2 m/s, and this
          ! happens in the COAMPS microphysics -dschanen 29 Sept 2009
          max_velocity = -2.0 ! m/s

        case ( "rgraupelm" )
          ixrm_bt = irgraupelm_bt
          ixrm_cl = irgraupelm_cl
          ixrm_mc = irgraupelm_mc

          max_velocity = -20. ! m/s

        case ( "Nrm" )
          ixrm_bt = iNrm_bt
          ixrm_cl = iNrm_cl
          ixrm_mc = iNrm_mc

          max_velocity = -9.1 ! m/s

        case ( "Nim" )
          ixrm_bt = iNim_bt
          ixrm_cl = iNim_cl
          ixrm_mc = iNim_mc

          max_velocity = -1.2 ! m/s

        case ( "Nsnowm" )
          ixrm_bt = iNsnowm_bt
          ixrm_cl = iNsnowm_cl
          ixrm_mc = iNsnowm_mc

          ! Morrison limit
!         max_velocity = -1.2 ! m/s
          ! Made up limit.  The literature suggests that it is quite possible
          ! that snow flake might achieve a terminal velocity of 2 m/s, and this
          ! happens in the COAMPS microphysics -dschanen 29 Sept 2009
          max_velocity = -2.0 ! m/s

        case ( "Ngraupelm" )
          ixrm_bt = iNgraupelm_bt
          ixrm_cl = iNgraupelm_cl
          ixrm_mc = iNgraupelm_mc

          max_velocity = -20. ! m/s

        case ( "Ncm" )
          ixrm_bt = iNcm_bt
          ixrm_cl = iNcm_cl
          ixrm_mc = iNcm_mc

          ! Use the rain water limit, since Morrison has no explicit limit on
          ! cloud water.  Presumably these numbers are never large.
          ! -dschanen 28 Sept 2009
          max_velocity = -9.1 ! m/s

        case default
          ixrm_bt = 0
          ixrm_cl = 0
          ixrm_mc = 0

          max_velocity = -9.1 ! m/s

        end select

        if ( l_stats_samp ) then

          ! Update explicit contributions to the hydrometeor species
          call stat_update_var( ixrm_mc, hydromet_mc(:,i), zt )

          ! Save prior value of the hydrometeors for determining total time
          ! tendency
          ! This kludge is to allow Ncm_bt to include the Ncm_act budget term
          ! that is calculated before microphysics is called.  The stat_begin_update
          ! subroutine is called prior to the GFDL activation code, so we can't call
          ! it here again. -meyern
          if ( ixrm_bt /= iNcm_bt ) then
            call stat_begin_update( ixrm_bt, hydromet(:,i) / real( dt ), zt )
          end if

        end if

        ! Set realistic limits on sedimentation velocities, following the
        ! numbers in the Morrison microphysics.

        do k = 1, gr%nzmax
          if ( clubb_at_least_debug_level( 1 ) ) then
            ! Print a warning if the velocity has a large magnitude or the
            ! velocity is in the wrong direction.
            if ( hydromet_vel_zt(k,i) < max_velocity .or. &   
                 hydromet_vel_zt(k,i) > zero_threshold ) then  

              write(fstderr,*) trim( hydromet_list(i) )// &
                " velocity at k = ", k, " = ", hydromet_vel_zt(k,i), "m/s"
            end if
          end if
          hydromet_vel_zt(k,i) = min( max( hydromet_vel_zt(k,i), max_velocity ), zero_threshold )
        end do ! k = 1..gr%nzmax

        ! Interpolate velocity to the momentum grid for a centered difference
        ! approximation of the sedimenation term.
        hydromet_vel(:,i) = zt2zm( hydromet_vel_zt(:,i) )
        hydromet_vel(gr%nzmax,i) = 0.0 ! Upper boundary condition

        ! Don't calculate if we aren't trying to predict Ncm in a meaningful way
        if ( trim( hydromet_list(i) ) /= "Ncm" .or. l_predictnc ) then

          ! Add implicit terms to the LHS matrix
          call microphys_lhs & 
               ( trim( hydromet_list(i) ), l_hydromet_sed(i), & ! In
                 dt, Kr, cloud_frac, nu_r_vert_res_dep, wm_zt, &  ! In
                 hydromet_vel(:,i), hydromet_vel_zt(:,i), & ! In
                 lhs ) ! Out

          if ( trim( hydromet_list(i) ) == "Ncm" ) then
            Ncm_in_cloud = hydromet(:,iiNcm) / max( cloud_frac, cloud_frac_min )
            call microphys_solve &
                 ( trim( hydromet_list(i) ), l_hydromet_sed(i), dt, cloud_frac, lhs, &
                   hydromet_mc(:,i)/max( cloud_frac, cloud_frac_min ),  & 
                   Ncm_in_cloud, err_code )

            hydromet(:,iiNcm) = Ncm_in_cloud * max( cloud_frac, cloud_frac_min )

          else
            call microphys_solve & 
                 ( trim( hydromet_list(i) ), l_hydromet_sed(i), dt, cloud_frac, &
                   lhs, hydromet_mc(:,i), hydromet(:,i), err_code )
          end if

        end if ! hydromet /= Ncm .or. l_predictnc

        if ( trim( micro_scheme ) == "khairoutdinov_kogan" ) then
          if ( i == iirrainm ) then
            ! Handle over-evaporation of rrainm and adjust rt and theta-l
            ! hydrometeor tendency arrays accordingly.
            do k = 1, gr%nzmax, 1
              if ( hydromet(k,i) < 0.0 ) then
                l_sed = .true.
                call adj_microphys_tndcy & 
                   ( hydromet_mc(:,i), wm_zt, hydromet_vel(:,i), hydromet_vel_zt(:,i), & 
                     Kr, nu_r_vert_res_dep, dt, k, l_sed, & 
                     hydromet(:,i), overevap_rate )

                ! overevap_rate is defined as positive.
                ! It is a correction factor.
                rcm_mc(k)  = rcm_mc(k) - overevap_rate

                thlm_mc(k) = thlm_mc(k) + ( Lv / ( Cp*exner(k) ) ) * overevap_rate

                ! Moved from adj_microphys_tndcy
                if ( l_stats_samp ) then

                  call stat_update_var_pt( irrainm_cond_adj, k,  & 
                                           overevap_rate, zt )
                end if

              else


                if ( l_stats_samp ) then

                  call stat_update_var_pt( irrainm_cond_adj, k,  & 
                                      0.0, zt )
                end if
                ! Joshua Faschinj December 2007
              end if

            end do ! k=1..gr%nzmax

          else if ( i == iiNrm ) then
            ! Handle over-evaporation similar to rrainm.  However, in the case
            ! of Nrm there is no effect on rtm or on thlm.
            ! Brian Griffin.  April 14, 2007.
            do k = 1, gr%nzmax, 1
              if ( hydromet(k,i) < 0.0 ) then
                l_sed = .true.
                call adj_microphys_tndcy & 
                     ( hydromet_mc(:,i), wm_zt, hydromet_vel(:,i), hydromet_vel_zt(:,i), & 
                       Kr, nu_r_vert_res_dep, dt, k, l_sed, & 
                       hydromet(:,i), overevap_rate )

                ! Moved from adj_microphys_tndcy
                if ( l_stats_samp ) then
                  call stat_update_var_pt( iNrm_cond_adj, k,  & 
                                           overevap_rate, zt )
                end if
              else
                if ( l_stats_samp ) then
                  call stat_update_var_pt( iNrm_cond_adj,k, 0.0, zt )
                end if

              end if ! Nrm(k) < 0
              ! Joshua Fasching December 2007
            end do ! k = 1..gr%nzmax

          end if ! i == rrainm else if i == Nrm
        end if ! trim( micro_scheme  ) == khairoutdinov_kogan

        if ( l_stats_samp ) then

          call stat_begin_update & 
             ( ixrm_cl, hydromet(:,i) / real( dt ), zt )

        end if

        ! Don't clip Ncm if we aren't trying to predict Ncm in a meaningful way
        if ( trim( hydromet_list(i) ) /= "Ncm" .or. l_predictnc ) then

          ! Clip all hydrometeor species to be >= zero
          if ( any( hydromet(:,i) < zero_threshold ) ) then
            hydromet_name = hydromet_list(i)
            if ( clubb_at_least_debug_level( 1 ) ) then
              do k = 1, gr%nzmax
                if ( hydromet(k,i) < zero_threshold ) then
                  write(fstderr,*) trim( hydromet_name ) //" < ", zero_threshold, &
                    " in advance_microphys at k= ", k
                end if
              end do
            end if

            ! If we're dealing with a mixing ratio and hole filling is enabled,
            ! then we apply the hole filling algorithm
            if ( hydromet_name(1:1) == "r" .and. l_hole_fill ) then
              ! Apply the hole filling algorithm
              call fill_holes_driver( 2, zero_threshold, "zt", &
                                     rho_ds_zt, rho_ds_zm, &
                                     hydromet(:,i) )

              ! If the hole filling algorithm failed, then we attempt to fill
              ! the missing mass with water vapor mixing ratio.
              ! We noticed this is needed for ASEX A209, particularly if Latin
              ! hypercube sampling is enabled.  -dschanen 11 Nov 2010
              do k = 2, gr%nzmax
                if ( hydromet(k,i) < zero_threshold ) then

                  ! Set temp to the time tendency applied to vapor and removed
                  ! from the hydrometeor.
                  temp = hydromet(k,i) / real( dt )

                  ! Adjust the tendency rvm_mc accordingly
                  rvm_mc(k) = rvm_mc(k) + temp

                  ! Adjust the tendency of thlm_mc according to whether the effect
                  ! is an evaporation or sublimation tendency.
                  select case ( trim( hydromet_name ) )
                  case( "rrainm" )
                    thlm_mc(k) = thlm_mc(k) - temp * ( Lv / ( Cp*exner(k) ) )
                  case( "ricem", "rsnowm", "rgraupelm" )
                    thlm_mc(k) = thlm_mc(k) - temp * ( Ls / ( Cp*exner(k) ) )
                  case default
                    stop "Fatal error in microphys_driver"
                  end select

                  ! Set the mixing ratio to 0
                  hydromet(k,i) = zero_threshold

                end if ! hydromet(k,i) < 0
              end do ! k = 2..gr%nzmax

              ! Boundary condition
              ! Rain, snow and graupel which is at the ghost point has presumably
              ! sedimented out of the model domain, and is not conserved.
              if ( hydromet(1,i) < zero_threshold ) then
                hydromet(1,i) = zero_threshold
              end if

            else
            ! This includes the case where the variable is a number
            ! concentration and is therefore not conserved.
            where ( hydromet(:,i) < zero_threshold ) hydromet(:,i) = zero_threshold

            end if ! Variable is a mixing ratio and l_hole_fill is true

          end if ! hydromet(:,i) < 0

        end if ! hydromet_list(i) /= Ncm .or. l_predictnc

        if ( l_stats_samp ) then

          ! Effects of clipping
          call stat_end_update( ixrm_cl, hydromet(:,i) / real( dt ), zt )

          ! Total time tendency
          call stat_end_update( ixrm_bt, hydromet(:,i) / real( dt ), zt )

        end if ! l_stats_samp

      end do ! i=1..hydromet_dim

    end if ! hydromet_dim > 0


    ! Call the ice diffusion scheme
    if ( trim( micro_scheme ) == "simplified_ice" ) then
      call ice_dfsn( dt, thlm, rcm, exner, p_in_Pa, rho, rcm_mc, thlm_mc )
    end if

    if ( l_stats_samp ) then
      if ( iiNcm > 0 ) then
        call stat_update_var( iNcm_in_cloud, &
               hydromet(:, iiNcm) / max( cloud_frac, cloud_frac_min ) , zt )
      end if

      call stat_update_var( iNcnm, Ncnm, zt )


    end if ! l_stats_samp

    if ( l_stats_samp .and. iirrainm > 0 ) then
      ! Rainfall rate (mm/day) should be defined on thermodynamic
      ! levels.  -Brian
      ! The absolute value of Vrr is taken because rainfall rate
      ! is a scalar quantity, and is therefore positive.
      call stat_update_var( irain_rate_zt,  & 
         ( hydromet(:,iirrainm)  & 
           * zm2zt( abs( hydromet_vel(:,iirrainm) ) ) ) & 
           * ( rho / rho_lw ) & 
           * real( sec_per_day ) * mm_per_m, zt )

      ! Precipitation Flux (W/m^2) should be defined on
      ! momentum levels.  -Brian
      ! Normally, a flux is a vector quantity.  Since rain obviously
      ! falls downward, the sign of the flux would normally be negative.
      ! However, it is generally a convention in meteorology to show
      ! Precipitation Flux as a positive downward quantity.  Thus, the
      ! absolute value of vrr is taken.
      call stat_update_var( iFprec,  & 
           ( zt2zm( hydromet(:,iirrainm) )  & 
             * abs( hydromet_vel(:,iirrainm) ) ) & 
             * ( rho_zm / rho_lw ) * rho_lw * Lv, zm )

      ! Store values of surface fluxes for statistics
      ! See notes above.
      call stat_update_var_pt( irain_rate_sfc, 1,  & 
           ( hydromet(2,iirrainm) & 
             * abs( zm2zt( hydromet_vel(:,iirrainm), 2 ) ) ) & 
             * ( rho(2) / rho_lw ) & 
             * real( sec_per_day ) * 1000.0, sfc )

      call stat_update_var_pt( irain_flux_sfc, 1, & 
           ( zt2zm( hydromet(:,iirrainm), 1 )  & 
             * abs( hydromet_vel(1,iirrainm) ) ) * ( rho_zm(1) / rho_lw )  & 
             * rho_lw * Lv, sfc )

      ! Also store the value of surface rain water mixing ratio.
      call stat_update_var_pt( irrainm_sfc, 1,  & 
           ( zt2zm( hydromet(:,iirrainm), 1 ) ), sfc )

    end if ! l_stats_samp

    call stats_accumulate_hydromet( hydromet )

!       Error Report
!       Joshua Fasching Feb 2008

    if ( fatal_error( err_code ) .and.  &
         clubb_at_least_debug_level( 1 ) ) then

      write(fstderr,*) "Error in advance_microphys"

      write(fstderr,*) "Intent(in)"

      write(fstderr,*) "thlm = ", thlm
      write(fstderr,*) "p_in_Pa = ", p_in_Pa
      write(fstderr,*) "exner = ", exner
      write(fstderr,*) "rho = ", rho
      write(fstderr,*) "rho_zm = ", rho_zm
      write(fstderr,*) "rtm = ", rtm
      write(fstderr,*) "rcm = ", rcm
      write(fstderr,*) "wm_zt = ", wm_zt
      write(fstderr,*) "wm_zm = ", wm_zm
      write(fstderr,*) "Kh_zm = ", Kh_zm
      write(fstderr,*) "pdf_params%thl1 = ", pdf_params%thl1
      write(fstderr,*) "pdf_params%thl2 = ", pdf_params%thl2
      write(fstderr,*) "pdf_params%mixt_frac = ", pdf_params%mixt_frac
      write(fstderr,*) "pdf_params%rc1 = ", pdf_params%rc1
      write(fstderr,*) "pdf_params%rc2 = ", pdf_params%rc2
      write(fstderr,*) "pdf_params%s1 = ", pdf_params%s1
      write(fstderr,*) "pdf_params%s2 = ", pdf_params%s2
      write(fstderr,*) "pdf_params%stdev_s1 = ", pdf_params%stdev_s1
      write(fstderr,*) "pdf_params%stdev_s2 = ", pdf_params%stdev_s2

      write(fstderr,*) "Intent(inout)"

      write(fstderr,*) "Ncnm = ", Ncnm
      write(fstderr,*) "hydromet = ", hydromet
      write(fstderr,*) "Intent(out)"
      write(fstderr,*) "rcm_mc = ", rcm_mc
      write(fstderr,*) "rvm_mc = ", rvm_mc
      write(fstderr,*) "thlm_mc = ", thlm_mc

    end if

    return
  end subroutine advance_microphys

!===============================================================================
  subroutine microphys_solve( solve_type, l_sed, dt, cloud_frac, lhs, & 
                                xrm_tndcy, xrm, err_code )

    ! Description:
    ! Solve the tridiagonal system for hydrometeor variable.

    ! References:
    !  None
    !---------------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use lapack_wrap, only:  & 
        tridag_solve !,& ! Procedure(s)
!       band_solve

    use error_code, only: &
        clubb_no_error ! Constant

    use stats_variables, only: & 
        zt,  & ! Variable(s)
        irrainm_ma, & 
        irrainm_sd, & 
        irrainm_dff, & 
        iricem_ma, & 
        iricem_sd, & 
        iricem_dff, & 
        irsnowm_ma, & 
        irsnowm_sd, & 
        irsnowm_dff, & 
        irgraupelm_ma, & 
        irgraupelm_sd, & 
        irgraupelm_dff, & 
        l_stats_samp, & 
        ztscr01, & 
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
        ztscr06, & 
        ztscr07, & 
        ztscr08, & 
        ztscr09

    use stats_variables, only: & 
        iNrm_ma, & 
        iNrm_sd, & 
        iNrm_dff, & 
        iNim_ma, & 
        iNim_sd, & 
        iNim_dff, & 
        iNsnowm_ma, & 
        iNsnowm_sd, & 
        iNsnowm_dff, & 
        iNgraupelm_ma, & 
        iNgraupelm_sd, & 
        iNgraupelm_dff

    use stats_variables, only: & 
        iNcm_ma, & 
        iNcm_dff


    use stats_type, only: stat_update_var_pt ! Procedure(s)

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: solve_type

    real(kind=time_precision), intent(in) :: dt ! Timestep     [s]

    logical, intent(in) :: &
      l_sed ! Whether sedimentation is included in lhs (T/F)

    ! Explicit contrbution to the hydrometeor, e.g. evaporation
    ! from Brian Griffin's K & K microphysics implementation
    real, intent(in), dimension(gr%nzmax) :: & 
      xrm_tndcy, &  !                                 [units/s]
      cloud_frac    ! Cloud fraction                  [-]

    ! Input/Output Variables
    real, intent(inout), dimension(3,gr%nzmax) :: & 
      lhs ! Left hand side

    real, intent(inout), dimension(gr%nzmax) :: & 
      xrm ! Hydrometeor being solved for              [units vary]

    ! Output Variables
    integer, intent(out) :: err_code

    ! Local Variables
    real, dimension(gr%nzmax) :: & 
      rhs ! Right hand side

    integer :: k, kp1, km1 ! Array indices

    integer :: & 
      ixrm_ma,  & ! Mean advection budget stats toggle
      ixrm_sd,  & ! Sedimentation budget stats toggle
      ixrm_dff    ! Diffusion budget stats toggle

    err_code = clubb_no_error  ! Initialize to the value for no errors

    ! Initializing ixrm_ma, ixrm_sd, and ixrm_dff in order to avoid compiler
    ! warnings.
    ixrm_ma  = 0
    ixrm_sd  = 0
    ixrm_dff = 0

    select case( solve_type )
    case( "rrainm" )
      ixrm_ma  = irrainm_ma
      ixrm_sd  = irrainm_sd
      ixrm_dff = irrainm_dff
    case( "ricem" )
      ixrm_ma  = iricem_ma
      ixrm_sd  = iricem_sd
      ixrm_dff = iricem_dff
    case( "rsnowm" )
      ixrm_ma  = irsnowm_ma
      ixrm_sd  = irsnowm_sd
      ixrm_dff = irsnowm_dff
    case( "rgraupelm" )
      ixrm_ma  = irgraupelm_ma
      ixrm_sd  = irgraupelm_sd
      ixrm_dff = irgraupelm_dff
    case( "Ncm" )
      ixrm_ma  = iNcm_ma
      ixrm_sd  = 0
      ixrm_dff = iNcm_dff
    case( "Nrm" )
      ixrm_ma  = iNrm_ma
      ixrm_sd  = iNrm_sd
      ixrm_dff = iNrm_dff
    case( "Nim" )
      ixrm_ma  = iNim_ma
      ixrm_sd  = iNim_sd
      ixrm_dff = iNim_dff
    case( "Nsnowm" )
      ixrm_ma  = iNsnowm_ma
      ixrm_sd  = iNsnowm_sd
      ixrm_dff = iNsnowm_dff
    case( "Ngraupelm" )
      ixrm_ma  = iNgraupelm_ma
      ixrm_sd  = iNgraupelm_sd
      ixrm_dff = iNgraupelm_dff
    case default
      ixrm_ma  = 0
      ixrm_sd  = 0
      ixrm_dff = 0
    end select


    ! RHS of equation, following Brian's method from the rain subroutine
    rhs(2:gr%nzmax-1)  & 
      = (xrm(2:gr%nzmax-1) / real( dt )  & ! Time tendency
      + xrm_tndcy(2:gr%nzmax-1))


    ! Boundary condition on the RHS
    rhs(1) = xrm(1) / real( dt )
    rhs(gr%nzmax) =  & 
       ( xrm(gr%nzmax) / real( dt ) + xrm_tndcy(gr%nzmax) )


    ! Solve system using tridag_solve. This uses LAPACK sgtsv,
    ! which relies on Gaussian elimination to decompose the matrix.
    call tridag_solve & 
         ( solve_type, gr%nzmax, 1, lhs(1,:), lhs(2,:), lhs(3,:), & 
           rhs, xrm, err_code )

    if ( l_stats_samp ) then

      do k = 1, gr%nzmax, 1

        km1 = max( k-1, 1 )
        kp1 = min( k+1, gr%nzmax )

        ! Finalize implicit contributions

        ! xrm term ma is completely implicit; call stat_update_var_pt.
        if ( solve_type == "Ncm" ) then
          ! For Ncm, we divide by cloud_frac when entering the subroutine, but
          ! do not multiply until we return from the subroutine, so we must
          ! account for this here for the budget to balance
          call stat_update_var_pt( ixrm_ma, k, & 
                ztscr01(k) * xrm(km1) * max( cloud_frac(k), cloud_frac_min ) & 
                + ztscr02(k) * xrm(k) * max( cloud_frac(k), cloud_frac_min ) & 
                + ztscr03(k) * xrm(kp1) * max( cloud_frac(k), cloud_frac_min ), zt)
        else
          call stat_update_var_pt( ixrm_ma, k, & 
                ztscr01(k) * xrm(km1) & 
                + ztscr02(k) * xrm(k) & 
                + ztscr03(k) * xrm(kp1), zt)
        end if

        ! xrm term sd is completely implicit; call stat_update_var_pt.
        if ( l_sed ) then
          call stat_update_var_pt( ixrm_sd, k, & 
                ztscr04(k) * xrm(km1) & 
                + ztscr05(k) * xrm(k) & 
                + ztscr06(k) * xrm(kp1), zt )
        end if

        ! xrm term dff is completely implicit; call stat_update_var_pt.
        if ( solve_type == "Ncm" ) then
          ! For Ncm, we divide by cloud_frac when entering the subroutine, but
          ! do not multiply until we return from the subroutine, so we must
          ! account for this here for the budget to balance
          call stat_update_var_pt( ixrm_dff, k, & 
                ztscr07(k) * xrm(km1) * max( cloud_frac(k), cloud_frac_min ) & 
                + ztscr08(k) * xrm(k) * max( cloud_frac(k), cloud_frac_min ) & 
                + ztscr09(k) * xrm(kp1) * max( cloud_frac(k), cloud_frac_min ), zt )
        else
          call stat_update_var_pt( ixrm_dff, k, & 
                ztscr07(k) * xrm(km1) & 
                + ztscr08(k) * xrm(k) & 
                + ztscr09(k) * xrm(kp1), zt )
        end if

      enddo ! 1..gr%nzmax

    end if ! l_stats_samp


    ! Boundary conditions on results
    !xrm(1) = xrm(2)
    ! Michael Falk, 7 Sep 2007, made this change to eliminate problems
    ! with anomalous rain formation at the top boundary.
    !        xrm(gr%nzmax) = 0
    !xrm(gr%nzmax) = xrm(gr%nzmax-1)
    ! eMFc

    return
  end subroutine microphys_solve

!===============================================================================
  subroutine microphys_lhs & 
             ( solve_type, l_sed, dt, Kr, cloud_frac, nu, wm_zt, &
               V_hm, V_hmt, &
               lhs )

  ! Description:
  !   Setup the matrix of implicit contributions to a term.
  !   Can include the effects of sedimentation, diffusion, and advection.
  !   The Morrison microphysics has an explicit sedimentation code, which is
  !   handled elsewhere.
  !
  ! Notes:
  !   Setup for tridiagonal system and boundary conditions should be the same as
  !   the original rain subroutine code.
  !-------------------------------------------------------------------------------

    use grid_class, only:  & 
        gr,    & ! Variable(s)
        zm2zt, & ! Procedure(s)
        zt2zm    ! Procedure(s)

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use diffusion, only:  & 
        diffusion_zt_lhs, & ! Procedure(s)
        diffusion_cloud_frac_zt_lhs

    use mean_adv, only:  & 
        term_ma_zt_lhs ! Procedure(s)

    use constants_clubb, only: &
      sec_per_day ! Variable(s)

    use stats_variables, only: & 
        irrainm_ma,   & ! Variable(s)
        irrainm_sd, & 
        irrainm_dff, & 
        iNrm_ma, & 
        iNrm_sd, & 
        iNrm_dff, & 
        iNcm_ma, &
        iNcm_dff, &
        iricem_ma, & 
        iricem_sd, & 
        iricem_dff, & 
        irsnowm_ma, & 
        irsnowm_sd, & 
        irsnowm_dff, & 
        irgraupelm_ma, & 
        irgraupelm_sd, & 
        irgraupelm_dff, & 
        ztscr01, & 
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
        ztscr06, & 
        ztscr07, & 
        ztscr08, & 
        ztscr09, & 
        l_stats_samp

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

    real, parameter :: &
      cloud_frac_thresh = 1.e-3 ! Minimum threshold on cloud fraction

    ! Input Variables
    character(len=*), intent(in) :: &
      solve_type  ! Description of which hydrometeor is being solved for.

    logical, intent(in) ::  & 
      l_sed    ! Whether to add a hydrometeor sedimentation term.

    real(kind=time_precision), intent(in) ::  & 
      dt       ! Model timestep                                           [s]

    real, dimension(gr%nzmax), intent(in) ::  & 
      nu       ! Background diffusion coefficient                         [m^2/s]

    real, intent(in), dimension(gr%nzmax) ::  & 
      cloud_frac, & ! Cloud fraction                                          [-]
      wm_zt,      & ! w wind component on thermodynamic levels                [m/s]
      V_hm,       & ! Sedimentation velocity of hydrometeor (momentum levels) [m/s]
      V_hmt,      & ! Sedimentation velocity of hydrometeor (thermo. levels)  [m/s]
      Kr            ! Eddy diffusivity for hydrometeor on momentum levels     [m^2/s]

    real, intent(out), dimension(3,gr%nzmax) :: & 
      lhs      ! Left hand side of tridiagonal matrix.

    ! Local Variables
    real, dimension(3) :: tmp

    real, dimension(gr%nzmax) :: & 
      cloud_frac_zt, & ! Cloud fraction on thermodynamic levels  [-]
      cloud_frac_zm    ! Cloud fraction on momentum levels       [-]

    ! Array indices
    integer :: k, km1, kp1

    !integer kp1

    integer :: & 
      ixrm_ma,  & ! Mean advection budget stats toggle
      ixrm_sd,  & ! Sedimentation budget stats toggle
      ixrm_dff    ! Diffusion budget stats toggle

    ! ----- Begin Code -----

    ! Initializing ixrm_ma, ixrm_sd, and ixrm_dff in order to avoid compiler
    ! warnings.
    ixrm_ma  = 0
    ixrm_sd  = 0
    ixrm_dff = 0

    select case( solve_type )
    case ( "rrainm" )
      ixrm_ma  = irrainm_ma
      ixrm_sd  = irrainm_sd
      ixrm_dff = irrainm_dff
    case ( "Nrm" )
      ixrm_ma  = iNrm_ma
      ixrm_sd  = iNrm_sd
      ixrm_dff = iNrm_dff
    case ( "ricem" )
      ixrm_ma  = iricem_ma
      ixrm_sd  = iricem_sd
      ixrm_dff = iricem_dff
    case ( "rsnowm" )
      ixrm_ma  = irsnowm_ma
      ixrm_sd  = irsnowm_sd
      ixrm_dff = irsnowm_dff
    case ( "rgraupelm" )
      ixrm_ma  = irgraupelm_ma
      ixrm_sd  = irgraupelm_sd
      ixrm_dff = irgraupelm_dff
    case ( "Ncm" )
      ixrm_ma  = iNcm_ma
      ixrm_sd  = 0
      ixrm_dff = iNcm_dff
    case default
      ixrm_ma  = 0
      ixrm_sd  = 0
      ixrm_dff = 0
    end select

    ! Determine cloud fraction for diffusion of Ncm
    if ( solve_type == "Ncm".and. l_in_cloud_Nc_diff ) then
      ! Impose a threshold on cloud fract to avoid a divide by 0.
      cloud_frac_zt = max( cloud_frac, cloud_frac_thresh )
!     cloud_frac_zt = 1.0 ! %% Debug
      ! Don't impose a threshold on cloud_frac_zm in the numerator.
      cloud_frac_zm = zt2zm( cloud_frac )
    end if

    ! Reset LHS Matrix for current timestep.
    lhs = 0.0

    ! Setup LHS Matrix
    do k = 2, gr%nzmax-1, 1

      km1 = max( k-1, 1 )
      kp1 = min( k+1, gr%nzmax )

      ! Main diagonal

      ! LHS time tendency.
      lhs(k_tdiag,k) = lhs(k_tdiag,k) + ( 1.0 / real( dt ) )


      ! All diagonals

      ! LHS eddy-diffusion term.
      if ( solve_type == "Ncm" .and. l_in_cloud_Nc_diff ) then
        lhs(kp1_tdiag:km1_tdiag,k) & 
          = lhs(kp1_tdiag:km1_tdiag,k) & 
          + diffusion_cloud_frac_zt_lhs &
            ( Kr(k), Kr(km1), cloud_frac_zt(k), cloud_frac_zt(k-1), &
              cloud_frac_zt(k+1), cloud_frac_zm(k), &
              cloud_frac_zm(k-1), &
              nu, gr%invrs_dzm(km1), gr%invrs_dzm(k), gr%invrs_dzt(k), k )
      else ! All other cases
        lhs(kp1_tdiag:km1_tdiag,k) & 
          = lhs(kp1_tdiag:km1_tdiag,k) & 
          + diffusion_zt_lhs( Kr(k), Kr(km1), nu,  & 
                              gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                              gr%invrs_dzt(k), k )
      end if

      ! LHS mean advection term.
      lhs(kp1_tdiag:km1_tdiag,k) & 
        = lhs(kp1_tdiag:km1_tdiag,k) & 
        + term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), gr%invrs_dzm(km1) )

      ! LHS hydrometeor sedimentation term.
      ! Note: the Morrison microphysics has its own sedimentation code, which
      ! is applied through the _mc terms to each hydrometeor species.
      ! Therefore, l_sed will always be false when the Morrison microphysics 
      ! is enabled.  -dschanen 24 Jan 2011
      if ( l_sed ) then
        if ( .not. l_upwind_diff_sed ) then
          lhs(kp1_tdiag:km1_tdiag,k) & 
            = lhs(kp1_tdiag:km1_tdiag,k) & 
            + sed_centered_diff_lhs( V_hm(k), V_hm(km1), gr%invrs_dzt(k), k )
        else
          lhs(kp1_tdiag:km1_tdiag,k) & 
            = lhs(kp1_tdiag:km1_tdiag,k) & 
            + sed_upwind_diff_lhs( V_hmt(k), V_hmt(kp1), gr%invrs_dzm(k), k )
        end if
      end if

      if ( l_stats_samp ) then

        ! Statistics:  implicit contributions to hydrometeor xrm.

        if ( solve_type == "Ncm" .and. ixrm_ma > 0 ) then
          tmp(1:3) = term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), &
            gr%invrs_dzm(km1) )

          ztscr01(k) = -tmp(3)
          ztscr02(k) = -tmp(2)
          ztscr03(k) = -tmp(1)

        else if ( ixrm_ma > 0 ) then
          tmp(1:3) = term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), &
            gr%invrs_dzm(km1) )

          ztscr01(k) = -tmp(3)
          ztscr02(k) = -tmp(2)
          ztscr03(k) = -tmp(1)
        end if

        if ( ixrm_sd > 0 .and. l_sed ) then
          if ( .not. l_upwind_diff_sed ) then
            tmp(1:3) = sed_centered_diff_lhs( V_hm(k), V_hm(km1), gr%invrs_dzt(k), k )
          else
            tmp(1:3) = sed_upwind_diff_lhs( V_hmt(k), V_hmt(kp1), gr%invrs_dzm(k), k )
          end if

          ztscr04(k) = -tmp(3)
          ztscr05(k) = -tmp(2)
          ztscr06(k) = -tmp(1)
        end if

        if ( solve_type == "Ncm" .and. l_in_cloud_Nc_diff &
             .and. ixrm_dff > 0 ) then
          tmp(1:3) &
            = diffusion_cloud_frac_zt_lhs &
              ( Kr(k), Kr(km1), cloud_frac_zt(k), cloud_frac_zt(k-1), &
                cloud_frac_zt(k+1), cloud_frac_zm(k), &
                cloud_frac_zm(k-1), &
                nu, gr%invrs_dzm(km1), gr%invrs_dzm(k), gr%invrs_dzt(k), k )
          ztscr07(k) = -tmp(3)
          ztscr08(k) = -tmp(2)
          ztscr09(k) = -tmp(1)

        else if ( solve_type == "Ncm" .and. ixrm_dff > 0 ) then
          tmp(1:3) & 
            = diffusion_zt_lhs( Kr(k), Kr(km1), nu,  & 
                                gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                                gr%invrs_dzt(k), k )
          ztscr07(k) = -tmp(3)
          ztscr08(k) = -tmp(2)
          ztscr09(k) = -tmp(1)

        else if ( ixrm_dff > 0 ) then
          tmp(1:3) & 
            = diffusion_zt_lhs( Kr(k), Kr(km1), nu,  & 
                                gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                                gr%invrs_dzt(k), k )
          ztscr07(k) = -tmp(3)
          ztscr08(k) = -tmp(2)
          ztscr09(k) = -tmp(1)
        end if

      end if ! l_stats_samp

    end do ! 2..gr%nzmax-1


    ! Boundary Conditions

    ! The hydrometeor eddy-diffusion term has zero-flux boundary conditions, meaning
    ! that amounts of a hydrometeor are not allowed to escape the model boundaries
    ! through the process of eddy-diffusion.  It should be noted that amounts of a
    ! hydrometeor are allowed to leave the model at the lower boundary through the
    ! process of hydrometeor sedimentation.  However, only the eddy-diffusion term
    ! contributes to the LHS matrix at the k=1 and k=gr%nzmax levels.  Thus, function
    ! diffusion_zt_lhs needs to be called at both the upper boundary level and the
    ! lower boundary level.


    ! Lower Boundary
    k   = 1
    km1 = max( k-1, 1 )
    kp1 = k+1
    ! Note:  In function diffusion_zt_lhs, at the k=1 (lower boundary) level,
    !        variables referenced at the km1 level don't factor into the equation.

    ! LHS time tendency at the lower boundary.
    lhs(k_tdiag,k) = lhs(k_tdiag,k) + ( 1.0 / real( dt ) )

    ! LHS eddy-diffusion term at the lower boundary.
    lhs(kp1_tdiag:km1_tdiag,k) &
      = lhs(kp1_tdiag:km1_tdiag,k) &
      + diffusion_zt_lhs( Kr(k), Kr(km1), nu,  &
                          gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                          gr%invrs_dzt(k), k )

    ! Here we apply the upwind differencing at the lower boundary.
    if ( l_sed .and. l_upwind_diff_sed ) then
      lhs(kp1_tdiag:km1_tdiag,k) & 
        = lhs(kp1_tdiag:km1_tdiag,k) & 
        + sed_upwind_diff_lhs( V_hmt(k), V_hmt(kp1), gr%invrs_dzm(k), k )
    end if

    if ( l_stats_samp ) then

      ! Statistics:  implicit contributions to hydrometeor xrm.

      if ( ixrm_dff > 0 ) then
        tmp(1:3) & 
          = diffusion_zt_lhs( Kr(k), Kr(km1), nu,  & 
                              gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                              gr%invrs_dzt(k), k )
        ztscr07(k) = -tmp(3)
        ztscr08(k) = -tmp(2)
        ztscr09(k) = -tmp(1)
      end if

      if ( ixrm_sd > 0 .and. l_sed .and. l_upwind_diff_sed ) then
        tmp(1:3) = sed_upwind_diff_lhs( V_hmt(k), V_hmt(kp1), gr%invrs_dzm(k), k )

        ztscr04(k) = -tmp(3)
        ztscr05(k) = -tmp(2)
        ztscr06(k) = -tmp(1)
      end if

    end if  ! l_stats_samp


    ! Upper Boundary
    k   = gr%nzmax
    km1 = max( k-1, 1 )

    ! LHS time tendency at the upper boundary.
    lhs(k_tdiag,k) = lhs(k_tdiag,k) + ( 1.0 / real( dt ) )

    ! LHS eddy-diffusion term at the upper boundary.
    lhs(kp1_tdiag:km1_tdiag,k) &
      = lhs(kp1_tdiag:km1_tdiag,k) &
      + diffusion_zt_lhs( Kr(k), Kr(km1), nu,  &
                          gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                          gr%invrs_dzt(k), k )

    if ( l_stats_samp ) then

      ! Statistics:  implicit contributions to hydrometeor xrm.

      if ( ixrm_dff > 0 ) then
        tmp(1:3) & 
          = diffusion_zt_lhs( Kr(k), Kr(km1), nu,  & 
                              gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                              gr%invrs_dzt(k), k )
        ztscr07(k) = -tmp(3)
        ztscr08(k) = -tmp(2)
        ztscr09(k) = -tmp(1)
      end if

    end if  ! l_stats_samp

    return
  end subroutine microphys_lhs

!===============================================================================
  pure function sed_centered_diff_lhs( V_hm, V_hmm1, invrs_dzt, level ) & 
    result( lhs )

    ! Description:
    ! Sedimentation of a hydrometeor:  implicit portion of the code.
    !
    ! The variable "hm" stands for one of the five hydrometeor variables
    ! currently in the code:  mean rain mixing ratio (rrainm), mean rain drop
    ! concentration (Nrm), mean ice mixing ratio (ricem), mean snow mixing
    ! ratio (rsnowm), or mean graupel mixing ratio (rgraupelm).  The variable
    ! "V_hm" stands for the sedimentation velocity of the appropriate
    ! hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - d(V_hm*hm)/dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - d( V_hm(t) * hm(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of the term is changed to
    !        a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d(hm)/dt equation.
    !
    ! This term is discretized as follows when using the centered-difference
    ! approximation:
    !
    ! The values of hm are found on the thermodynamic levels, while the values
    ! of V_hm are found on the momentum levels.  The variable hm is
    ! interpolated to the intermediate momentum levels.  At the intermediate
    ! momentum levels, the interpolated values of hm are multiplied by the
    ! values of V_hm.  Then, the derivative of (hm*V_hm) is taken over the
    ! central thermodynamic level.
    !
    ! -----hmp1------------------------------------------------ t(k+1)
    !
    ! =============hm(interp)=====V_hm========================= m(k)
    !
    ! -----hm--------------------------------d(V_hm*hm)/dz----- t(k)
    !
    ! =============hm(interp)=====V_hmm1======================= m(k-1)
    !
    ! -----hmm1------------------------------------------------ t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1),
    ! respectively.  The letter "t" is used for thermodynamic levels and the
    ! letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )
    !
    !
    ! Conservation Properties:
    !
    ! When a hydrometeor is sedimented to the ground (or out the lower
    ! boundary of the model), it is removed from the atmosphere (or from the
    ! model domain).  Thus, the quantity of the hydrometeor over the entire
    ! vertical domain should not be conserved due to the process of
    ! sedimentation.  Thus, not all of the column totals in the left-hand side
    ! matrix should be equal to 0.  Instead, the sum of all the column totals
    ! should equal the flux of hm out the bottom (zm(1) level) of the domain,
    ! -V_hm(1) * ( D(2)*hm(1) + C(2)*hm(2) ), where the factor in parentheses
    ! is the interpolated value of hm at the zm(1) level.  Furthermore, most
    ! of the individual column totals should sum to 0, but the 1st and 2nd
    ! (from the left) columns should combine to sum to the flux out the bottom
    ! of the domain.
    !
    ! To see that this modified conservation law is satisfied, compute the
    ! sedimentation of hm and integrate vertically.  In discretized matrix
    ! notation (where "i" stands for the matrix column and "j" stands for the
    ! matrix row):
    !
    !  -V_hm(1) * ( D(2)*hm(1) + C(2)*hm(2) )
    !     =
    !  Sum_j Sum_i ( 1/invrs_dzt )_i ( d (V_hm * weights_hm) / dz )_ij hm_j.
    !
    ! The left-hand side matrix, ( d (V_hm * weights_hm) / dz )_ij, is
    ! partially written below.  The sum over i in the above equation removes
    ! invrs_dzt everywhere from the matrix below.  The sum over j leaves the
    ! column totals and the flux at zm(1) that are desired.
    !
    ! Left-hand side matrix contributions from the sedimentation term (only);
    ! first four vertical levels:
    !
    !     ------------------------------------------------------------------->
    !k=1 |           0                     0                       0
    !    |
    !k=2 |   -invrs_dzt(k)       +invrs_dzt(k)           +invrs_dzt(k)
    !    |     *V_hm(k-1)*D(k)     *[ V_hm(k)*B(k)         *V_hm(k)*A(k)
    !    |                           -V_hm(k-1)*C(k) ]
    !    |
    !k=3 |           0           -invrs_dzt(k)           +invrs_dzt(k)
    !    |                         *V_hm(k-1)*D(k)         *[ V_hm(k)*B(k)
    !    |                                                   -V_hm(k-1)*C(k) ]
    !    |
    !k=4 |           0                     0             -invrs_dzt(k)
    !    |                                                 *V_hm(k-1)*D(k)
    !    |
    !   \ /
    !
    ! The variables A(k), B(k), C(k), and D(k) are weights of interpolation
    ! around the central thermodynamic level (k), such that:
    !
    ! A(k) = ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ),
    ! B(k) = 1 - [ ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
    !      = 1 - A(k);
    ! C(k) = ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) ), and
    ! D(k) = 1 - [ ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
    !      = 1 - C(k).
    !
    ! Furthermore, for all intermediate thermodynamic grid levels (as long as
    ! k /= gr%nzmax and k /= 1), the four weighting factors have the following
    ! relationships:  A(k) = C(k+1) and B(k) = D(k+1).
    !
    ! Note:  The superdiagonal term from level 3 and both the main diagonal
    !        and superdiagonal terms from level 4 are not shown on this
    !        diagram.

    ! References:
    ! None

    ! Notes:  
    !   Both COAMPS Microphysics and Brian Griffin's implementation use
    !   Khairoutdinov and Kogan (2000) for the calculation of rain
    !   mixing ratio and rain droplet number concentration sedimentation
    !   velocities, but COAMPS has only the local parameterization.
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

    integer, parameter :: & 
      t_above = 1,    & ! Index for upper thermodynamic level grid weight.
      t_below = 2       ! Index for lower thermodynamic level grid weight.

    ! Input Variables
    real, intent(in) :: & 
      V_hm,      & ! Sedimentation velocity of hydrometeor (k)     [m/s]
      V_hmm1,    & ! Sedimentation velocity of hydrometeor (k-1)   [m/s]
      invrs_dzt    ! Inverse of grid spacing (k)                   [m]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real, dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      mk,    & ! Momentum level directly above central thermodynamic level.
      mkm1     ! Momentum level directly below central thermodynamic level.

    ! ---- Begin Code ----

    ! Momentum level (k) is between thermodynamic level (k+1)
    ! and thermodynamic level (k).
    mk   = level
    ! Momentum level (k-1) is between thermodynamic level (k)
    ! and thermodynamic level (k-1).
    mkm1 = level - 1

    ! Note:  The code is now written so that V_hm has been pulled inside of the
    !        derivative.  The sedimentation term is now of the form -d(V_hm*hm)/dz,
    !        rather than of the form -V_hm d(hm)/dz.  The term has been
    !        re-discretized in a conservative manner and the results are listed
    !        below, with the old code commented out.

    if ( level == 1 ) then

      ! k = 1 (bottom level); lower boundary level; no effects.

      ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
      lhs(kp1_tdiag) = 0.0

      ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
      lhs(k_tdiag)   = 0.0

      ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
      lhs(km1_tdiag) = 0.0


    else if ( level > 1 .and. level < gr%nzmax ) then

      ! Most of the interior model; normal conditions.

      ! Vince Larson pulled V_hm inside derivative to make conservative.
      ! 13 Dec 2007
      !
      ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
!     lhs(kp1_tdiag)  &
!       = + V_hmzt * invrs_dzt * gr%weights_zt2zm(t_above,mk)

!     ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
!     lhs(k_tdiag)  &
!       = + V_hmzt * invrs_dzt * (   gr%weights_zt2zm(t_below,mk)  &
!                                  - gr%weights_zt2zm(t_above,mkm1)  )
!
!     ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
!     lhs(km1_tdiag)  &
!       = - V_hmzt * invrs_dzt * gr%weights_zt2zm(t_below,mkm1)

      ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
      lhs(kp1_tdiag)  & 
        = + invrs_dzt * V_hm * gr%weights_zt2zm(t_above,mk)

      ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
      lhs(k_tdiag)  & 
        = + invrs_dzt * (   V_hm * gr%weights_zt2zm(t_below,mk) & 
                          - V_hmm1 * gr%weights_zt2zm(t_above,mkm1)  )

      ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
      lhs(km1_tdiag)  & 
        = - invrs_dzt * V_hmm1 * gr%weights_zt2zm(t_below,mkm1)

      !  End Vince Larson change


    else if ( level == gr%nzmax ) then

      ! k = gr%nzmax (top level); upper boundary level; no flux.

      ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
      lhs(kp1_tdiag) = 0.0

      ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
      lhs(k_tdiag)   = 0.0

      ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
      lhs(km1_tdiag) = 0.0


    end if

    return
  end function sed_centered_diff_lhs

!---------------------------------------------------------------------------
  pure function sed_upwind_diff_lhs( V_hmt, V_hmtp1, invrs_dzm, level ) & 
    result( lhs )

  ! Description:
  !   Setup the LHS matrix (implicit component) for the upwind difference
  !   approximation for the sedimentation of a hydrometeor.
  ! 
  ! References:
  !   None
  !---------------------------------------------------------------------------
    use grid_class, only:  & 
        gr ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

    ! Input Variables
    real, intent(in) :: & 
      V_hmt,   & ! Sedimentation velocity of hydrometeor (thermo. levels) (k)   [m/s]
      V_hmtp1, & ! Sedimentation velocity of hydrometeor (thermo. levels) (k+1) [m/s]
      invrs_dzm  ! Inverse of grid spacing (k)                   [m]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real, dimension(3) :: lhs

    ! ---- Begin Code ----

    ! Sedimention is always a downward process, so we omit the upward case
    ! (i.e. the V_hmt variable will always be negative).
    if ( level == gr%nzmax ) then

      ! k = gr%nzmax (top level); upper boundary level; no flux.

      ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
      lhs(kp1_tdiag) = 0.0

      ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
      lhs(k_tdiag)   = 0.0

      ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
      lhs(km1_tdiag) = 0.0


    else  

      ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
      lhs(kp1_tdiag)  & 
        = + invrs_dzm * V_hmtp1

      ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
      lhs(k_tdiag)  & 
        = - invrs_dzm * V_hmt

      ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
      lhs(km1_tdiag) = 0.0

    end if

    return
  end function sed_upwind_diff_lhs
!===============================================================================
  subroutine adj_microphys_tndcy( xrm_tndcy, wm_zt, V_hm, V_hmt, Kr, nu, & 
                                  dt, level, l_sed, & 
                                  xrm, overevap_rate )

    ! DESCRIPTION:  Correction for the over-evaporation of a hydrometeor.
    !
    ! If a small amount of a hydrometeor (such as rain water) gets diffused
    ! into an area that is very dry (such as right above the cloud top), the
    ! hydrometeor (rain water) will have a very high rate of evaporation and
    ! will evaporate entirely in a short amount of time.  However, the
    ! evaporation rate is computed instantaneously at a given moment in time.
    !  This rate is then projected over the entire length of the given
    ! timestep.  Therefore, a high-enough rate of evaporation combined with a
    ! small-enough amount of the hydrometeor (rain water) and a long-enough
    ! timestep will cause the hydrometeor value (rain water mixing ratio) to
    ! be negative by the end of the timestep.  Therefore, a correction factor
    ! needs to be imposed on the evaporation rate so that the amount of the
    ! hydrometeor (rain water mixing ratio) does not fall below 0.
    !
    ! Besides over-evaporation of a hydrometeor, other factors may come into
    ! play that cause the value of a hydrometeor to fall below 0.  These
    ! factors are due to the nature of implicit discretization and numerical
    ! errors.  In a nutshell, the eddy diffusion parameter used currently in
    ! this model smooths out the entire hydrometeor profile as a whole at
    ! every timestep.  This smoothing may cause negative values at certain
    ! levels.  Also, mean advection and hydrometeor sedimentation can cause
    ! negative values to occur in the hydrometeor.  This can happen in places
    ! where the profile abruptly goes from a large positive value to 0 (such
    ! as at cloud top).  The nature of the discretization of taking a
    ! derivative at these levels may cause negative values of a hydrometeor.
    !
    ! This subroutine is called only if a hydrometeor at a certain level
    ! contains a negative value.  First, this subroutine uses the same methods
    ! that the model statistical code uses in computing budget terms in order
    ! to determine what factors effected the value of the given hydrometeor
    ! during the timestep that was just solved for.  The mean advection,
    ! sedimentation, and diffusion budget terms are all computed.  These three
    ! terms are then added together to make up the total transport and
    ! sedimentation tendency.  This tendency is then added to the total
    ! microphysical tendency to find the overall hydrometeor tendency.  The
    ! overall hydrometeor tendency is then multiplied by the timestep length
    ! to find the net change in the hydrometeor over the last timestep.  This
    ! net change is then added to the current value of the hydrometeor in
    ! order to find the value of the hydrometeor at the previous timestep.
    ! This method has been well tested and produces accurate results.
    !
    ! Once the value of the hydrometeor at the previous timestep has been
    ! found, the net change in the hydrometeor due to ONLY mean advection,
    ! diffusion, and sedimentation is calculated.  This net change is added to
    ! the value of the hydrometeor at the previous timestep.  If the new value
    ! is below zero, then the negative value of the hydrometeor was caused by
    ! the mean advection, diffusion, and sedimentation terms.  The
    ! microphysical terms (evaporation) did not cause the negative value.
    ! There was no over-evaporation and the evaporation rate can be set to 0.
    ! However, if the new value of the hydrometeor is greater than or equal to
    ! 0, then the microphysical tendencies (evaporation) did cause the
    ! hydrometeor array to have negative values.  The amount of hydrometeor
    ! evaporated is set equal to the amount that was left-over after the
    ! transport and sedimentation effects were added in.  The evaporation rate
    ! is that amount divided by the timestep.  This can be viewed as the
    ! timestep-average evaporation rate, whereas the rate previously
    ! calculated can be viewed as the instantaneous evaporation rate.  The
    ! amount of the hydrometeor that was over-evaporated is the amount of the
    ! hydrometeor that is negative.  The over-evaporation rate is that amount
    ! divided by the timestep.
    !
    ! It should be noted that this is important because the rain water mixing
    ! ratio time tendency (drr/dt) due to microphysics at every level is
    ! incorporated into the total water mixing ratio (rtm) and liquid water
    ! potential temperature (thlm) equations.  Any artificial excess in
    ! evaporation will artificially increase water vapor, and thus rtm, and
    ! artificially decrease thlm (due to evaporative cooling).  This may
    ! result in an artificial increase in cloud water.
    !
    ! rrainm_mc_tndcy = rrainm_cond + rrainm_auto + rrainm_accr
    ! rtm_mc  = - rrainm_mc_tndcy
    ! thlm_mc = ( Lv / (Cp*exner) ) * rrainm_mc_tndcy
    !
    ! Anyplace where rrainm drops below zero due to microphysics, there is too
    ! much evaporation rate for the timestep, so rrainm_cond is too negative.
    ! We must add in the over-evaporated amount of rrainm/dt to make the rate
    ! accurate.  The over-evaporated amount is being defined as a positive
    ! scalar, so that:  overevap_rrainm = -rrainm (where rrainm < 0) -- this
    ! makes overevap_rrainm positive.
    !
    ! New cond/evap rate = rrainm_cond + overevap_rrainm/dt
    ! (overevap_rate = overevap_rrainm/dt)
    ! -- since rrainm_cond can only be negative (we don't allow rain droplets
    ! to grow by condensation) and overevap_rrainm/dt can only be positive (we
    ! define it that way), the new cond/evap rate will be less negative, which
    ! is what we want.
    !
    ! To update the effects of microphysics on rtm and thl:
    !
    ! rtm_mc = rtm_mc - overevap_rate
    ! thlm_mc = thlm_mc + ( Lv / (Cp*exner) ) * overevap_rate
    !
    ! This is done in the subroutine which calls this one.
    !
    ! If the hydrometeor is negative due to reasons besides over-evaporation,
    ! the value is clipped.  This is statistically stored in the clipping
    ! array.  This is also done in the subroutine which calls this one.
    !
    ! Brian Griffin.

    use grid_class, only:  & 
        gr,  & ! Variable(s) 
        zm2zt ! Procedure(s)

    use stats_precision, only: & 
        time_precision ! Variable(s)

    use diffusion, only:  & 
        diffusion_zt_lhs ! Procedure(s)

    use mean_adv, only:  & 
        term_ma_zt_lhs ! Procedure(s)

    use parameters_microphys, only: &
        l_upwind_diff_sed ! Variable(s)

    implicit none

    ! Input variables.

    real, dimension(gr%nzmax), intent(in) :: &
      xrm_tndcy, & ! Hydrometeor microphysical tendency.                      [hm_units/s]
      wm_zt,     & ! Vertical velocity (thermo. levels).                      [m/s]
      V_hm,      & ! Sedimentation velocity (interpolated to moment. levels). [m/s]
      V_hmt,     & ! Sedimentation velocity (thermo. levels).                 [m/s]
      Kr           ! Eddy diffusivity for hydrometeors (m-lev).               [m^2/s]

    real, dimension(gr%nzmax), intent(in) :: nu  ! Diffusion coefficient      [m^2/s]

    real(kind=time_precision), intent(in) :: dt  ! Timestep   [s]

    integer, intent(in) :: level  ! Vertical grid index

    logical, intent(in) :: l_sed   ! Whether to add a sedimentation term


    ! Input/output variable.

    real, dimension(gr%nzmax), intent(inout) :: &
      xrm  ! Hydrometeor.  [hm_units]

    ! Output variable.

    ! Excess evaporation rate.
    real, intent(out) :: overevap_rate                 ! [hm_units/s]

    ! Local variables.
    real :: ma_subdiag   ! Term to be multiplied by xrm(k-1) in m.a. eq.
    real :: ma_maindiag  ! Term to be multiplied by xrm(k) in m.a. eq.
    real :: ma_supdiag   ! Term to be multiplied by xrm(k+1) in m.a. eq.
    real :: sd_subdiag   ! Term to be multiplied by xrm(k-1) in sed. eq.
    real :: sd_maindiag  ! Term to be multiplied by xrm(k) in sed. eq.
    real :: sd_supdiag   ! Term to be multiplied by xrm(k+1) in sed. eq.
    real :: df_subdiag   ! Term to be multiplied by xrm(k-1) in diff. eq.
    real :: df_maindiag  ! Term to be multiplied by xrm(k) in diff. eq.
    real :: df_supdiag   ! Term to be multiplied by xrm(k+1) in diff. eq.

    real :: ma_tndcy     ! Mean advection tendency  [hm_units/s]
    real :: sd_tndcy     ! Sedimentation tendency   [hm_units/s]
    real :: df_tndcy     ! Diffusion tendency       [hm_units/s]

    real :: trnsprt_sed_tndcy ! Total transport and sedimentation tendency.
    real :: tot_tndcy         ! Overall hydrometeor total tendency.
    real :: xrm_chge          ! Total change in hydrometeor over last t.s.
    real :: xrm_old           ! Value of hydrometeor at previous timestep.
    real :: xrm_chge_trsed    ! Net change in hm. due to only transport/sed.
    real :: xrm_trsed_only    ! New hm. val. due only to transport/sed.

!   real :: evap_amt          ! The actual evaporation amount over the t.s.
!   real :: evap_rate         ! The time-averaged rate.
    real :: overevap_amt      ! The amount of h.m. that was over-evap.

    real, dimension(1:3) :: tmp

    integer :: k, km1, kp1

    !
    !integer ::  &
    !ixrm_cond_adj  ! Adjustment to xrm evaporation rate due to over-evap.

    !select case( solve_type )
    !case( "rrainm" )
    !  ixrm_cond_adj  = irrainm_cond_adj
    !case( "Nrm" )
    !  ixrm_cond_adj  = iNrm_cond_adj
    !end select
    !
    ! Joshua Fasching 2007

    k = level
    km1 = max( k-1, 1 )
    kp1 = min( k+1, gr%nzmax )


    ! Mean advection tendency component

    ! The implicit (LHS) value of the mean advection component of the equation
    ! used during the timestep that was just solved for.
    tmp(1:3) = term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), gr%invrs_dzm(km1) )

    ma_subdiag  = -tmp(3) ! subdiagonal
    ma_maindiag = -tmp(2) ! main diagonal
    ma_supdiag  = -tmp(1) ! superdiagonal

    ma_tndcy = & 
      + ma_subdiag  * xrm(km1) & 
      + ma_maindiag * xrm(k) & 
      + ma_supdiag  * xrm(kp1)


    ! Sedimentation tendency component

    if ( l_sed ) then

      if ( .not. l_upwind_diff_sed ) then
        ! The implicit (LHS) value of the sedimentation component of the equation
        ! used during the timestep that was just solved for.
        tmp(1:3) = sed_centered_diff_lhs( V_hm(k), V_hm(km1), gr%invrs_dzt(k), k )

        sd_subdiag  = -tmp(3) ! subdiagonal
        sd_maindiag = -tmp(2) ! main diagonal
        sd_supdiag  = -tmp(1) ! superdiagonal

        sd_tndcy =  & 
          + sd_subdiag  * xrm(km1) & 
          + sd_maindiag * xrm(k) & 
          + sd_supdiag  * xrm(kp1)

      else ! Upwind differencing approximation

        ! The implicit (LHS) value of the sedimentation component of the equation
        ! used during the timestep that was just solved for.
        tmp(1:3) = sed_upwind_diff_lhs( V_hmt(k), V_hmt(kp1), gr%invrs_dzm(k), k )

        sd_maindiag = -tmp(2) ! main diagonal
        sd_supdiag  = -tmp(1) ! superdiagonal

        sd_tndcy =  & 
          + sd_supdiag  * xrm(kp1) & 
          + sd_maindiag * xrm(k)

      end if ! ~l_upwind_diff_sed 

    else

      sd_tndcy = 0.0

    end if


    ! Diffusion tendency component

    ! The implicit (LHS) value of the diffusion component of the equation used
    ! during the timestep that was just solved for.
    tmp(1:3) & 
      = diffusion_zt_lhs( Kr(k), Kr(km1), nu,  & 
                          gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                          gr%invrs_dzt(k), k )

    df_subdiag  = -tmp(3) ! subdiagonal
    df_maindiag = -tmp(2) ! main diagonal
    df_supdiag  = -tmp(1) ! superdiagonal

    df_tndcy =  & 
      + df_subdiag  * xrm(km1) & 
      + df_maindiag * xrm(k) & 
      + df_supdiag  * xrm(kp1)


    ! Total transport and sedimentation tendency
    trnsprt_sed_tndcy = ma_tndcy + df_tndcy + sd_tndcy

    ! Overall hydrometeor tendency
    tot_tndcy = trnsprt_sed_tndcy + xrm_tndcy(k)

    ! The net amount of change in the hydrometeor over the last timestep.
    xrm_chge = tot_tndcy * real( dt )

    ! The value of xrm at the previous timestep.
    xrm_old = xrm(k) - xrm_chge

    ! The net amount of change in the hydrometeor due to only the transport
    ! (mean advection and diffusion) and sedimentation terms.
    xrm_chge_trsed = trnsprt_sed_tndcy * real( dt )

    ! The new value of the hydrometeor at this timestep due to only
    ! the transport and sedimentation terms.
    xrm_trsed_only = xrm_old + xrm_chge_trsed

    if ( xrm_trsed_only >= 0.0 ) then
      ! The negative value of hydrometeor (xrm) is due ONLY to microphysical
      ! tendencies, namely the over-evaporation of xrm.
      ! Find the actual amount of the hydrometeor that evaporated during the
      ! timestep to make the value of xrm go to 0.
      !evap_amt = -xrm_trsed_only
      ! Divide by the timestep to find the actual evaporation rate.
      !evap_rate = evap_amt / dt
      ! The amount of the hydrometeor that was artificially excessively
      ! evaporated.  Define as positive.
      overevap_amt = -xrm(k)
      ! Divide by the timestep to find the over-evaporation rate.  Define as
      ! positive.  This rate should also be the difference between the
      ! computed evaporation rate (xrm_tndcy) and the actual evaporation rate
      ! (evap_rate).
      overevap_rate = overevap_amt / real( dt )
      ! Reset the value of the hydrometeor (xrm) to 0.
      xrm(k) = 0.0
    else
      ! The negative value of hydrometeor (xrm) is due to transport (mean
      ! advection and diffusion) and sedimentation.  Find the actual amount of
      ! the hydrometeor that evaporated during the timestep to make the value
      ! of xrm go to 0.  Even though the microphysical tendency portion of the
      ! code may have computed an evaporation rate, we figure that the
      ! transport and sedimentation terms made the value of the hydrometeor
      ! negative, so we say that the evaporation amount and rate is 0.
      !evap_amt = 0.0
      !evap_rate = 0.0
      ! The amount of the hydrometeor that was artificially excessively
      ! evaporated.  Define as positive.  In this case, any evaporation that
      ! was computed is considered to be over-evaporation.
      overevap_amt = -xrm_tndcy(k) * real( dt )
      overevap_rate = -xrm_tndcy(k)
      ! Currently reset xrm to xrm_trsed_only.  This is done to make the
      ! statistical budget for xrm balance correctly.  The value of xrm(k)
      ! will still be negative at this this point.  However, it will be less
      ! negative because it has been adjusted for over-evaporation.  The
      ! remaining negative value of hydrometeor xrm, which is due to transport
      ! and sedimentation, will be zeroed out in clipping in the subroutine
      ! that calls this one.
      xrm(k) = xrm_trsed_only
    end if

    return
  end subroutine adj_microphys_tndcy

!===============================================================================

  subroutine cleanup_microphys( )

    ! Description:
    !   De-allocate arrays used by the microphysics
    ! References:
    !   None
    !-------------------------------------------------------------------------

    implicit none

    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( hydromet_list ) ) then
      deallocate( hydromet_list )
    end if

    if ( allocated( l_hydromet_sed ) ) then
      deallocate( l_hydromet_sed )
    end if

    if ( trim( micro_scheme ) == "morrison-gettelman" ) then
      call pbuf_deallocate()
    end if

    return
  end subroutine cleanup_microphys
!===============================================================================

end module microphys_driver
