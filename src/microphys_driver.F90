! $Id$
!===============================================================================
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
  !-------------------------------------------------------------------------

  use parameters_microphys, only: &
    l_in_cloud_Nc_diff,           & ! Use in cloud values of Nc for diffusion
    l_cloud_sed,                  & ! Cloud water sedimentation (K&K or no microphysics)
    l_ice_micro,                  & ! Compute ice (COAMPS / Morrison)
    l_graupel,                    & ! Compute graupel (Morrison)
    l_hail,                       & ! See module_mp_graupel for a description
    l_seifert_beheng,             & ! Use Seifert and Beheng (2001) warm drizzle (Morrison)
    l_predictnc,                  & ! Predict cloud droplet number conc (Morrison)
    l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
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
    LH_seed,                      & ! Integer seed for the Mersenne twister
    l_local_kk,                   & ! Use local formula for K&K
    l_upwind_diff_sed,            & ! Use the upwind differencing approx. for sediementation
    l_var_covar_src,              & ! Flag for using variance and covariance src terms
    micro_scheme,                 & ! The microphysical scheme in use
    hydromet_list,                & ! Names of the hydrometeor species
    microphys_start_time,         & ! When to start the microphysics [s]
    sigma_g,                      & ! Parameter used in the cloud droplet sedimentation code
    Ncnm_initial,                 & ! Initial value for Ncnm (K&K)
    Nc0_in_cloud                    ! Initial value for Ncm (K&K, l_cloud_sed, Morrison)

  use parameters_microphys, only: &
    LH_microphys_interactive,     & ! Feed the subcolumns into the microphysics and allow feedback
    LH_microphys_non_interactive, & ! Feed the subcolumns into the microphysics with no feedback
    LH_microphys_disabled           ! Disable latin hypercube entirely

  use parameters_microphys, only: &
    l_silhs_KK_convergence_adj_mean ! Clip source adjustment terms on mean instead of individual
                                    ! sample points to test convergence with KK analytic

  use constants_clubb, only: &
    cloud_frac_min

  use phys_buffer, only: & ! Used for placing wp2_zt in morrison_gettelman microphysics
    pbuf_init,           &
    pbuf_deallocate

  implicit none

  ! Subroutines
  public :: advance_microphys, &
            init_microphys, &
            cleanup_microphys

  private :: microphys_solve, &
             microphys_lhs, &
             microphys_rhs

  ! Functions
  private :: sed_centered_diff_lhs, &
             sed_upwind_diff_lhs, &
             term_turb_sed_lhs, &
             term_turb_sed_rhs

  ! Variables
  logical, private, allocatable, dimension(:) :: &
    l_hydromet_sed ! Whether to sediment variables

!$omp threadprivate( l_hydromet_sed )

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

! Adding coefficient variable for clex9_oct14 case to reduce NNUCCD and NNUCCC
    use module_mp_graupel, only: &
        NNUCCD_REDUCE_COEF, &
        NNUCCC_REDUCE_COEF
! Change by Marc Pilon on 11/16/11

    use array_index, only: & 
        iirrainm,    & ! Variables
        iiNrm,       &
        iirsnowm,    &
        iiricem,     &
        iirgraupelm, &
        iiNcnm,      &
        iiNsnowm,    & 
        iiNim,       &
        iiNgraupelm, &
        iiNcm          ! Note: Ncm is not part of CLUBB's PDF.

    use parameters_microphys, only: &
        morrison_no_aerosol, &  ! Constants
        morrison_power_law,  &
        morrison_lognormal

    use parameters_microphys, only: &
        KK_evap_Supersat_exp, & ! Variable(s)
        KK_evap_rr_exp,       &
        KK_evap_Nr_exp,       &
        KK_auto_rc_exp,       &
        KK_auto_Nc_exp,       &
        KK_accr_rc_exp,       &
        KK_accr_rr_exp,       &
        KK_mvr_rr_exp,        &
        KK_mvr_Nr_exp,        &
        KK_Nrm_evap_nu,       &
        rrp2_on_rrm2_cloud,   &
        Nrp2_on_Nrm2_cloud,   &
        Ncnp2_on_Ncnm2_cloud, &
        rrp2_on_rrm2_below,   &
        Nrp2_on_Nrm2_below,   &
        Ncnp2_on_Ncnm2_below, &
        C_evap,               &
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

    use constants_clubb, only: &
        one,        & ! Constant(s)
        two_thirds, &
        one_third,  &
        zero, &
        cm3_per_m3

    use KK_fixed_correlations, only: &
        setup_KK_corr ! Procedure(s)

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
        ini_micro,             & ! Subroutine
        l_use_CLUBB_pdf_in_mg    ! Flag

    use microp_aero, only: &
        ini_microp_aero ! Subroutine

    use constants_clubb, only: &
        fstderr    ! Constant

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

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

#ifdef LATIN_HYPERCUBE
    use latin_hypercube_arrays, only: &
        setup_corr_varnce_array   ! Procedure(s)

    use mt95, only: &
        genrand_intg
#endif /* LATIN_HYPERCUBE */

    implicit none

    ! Constant Parameters
    logical, parameter :: &
      l_write_to_file = .true. ! If true, will write case information to a file.

    character(len=*), parameter :: &
      corr_input_path = "../input/case_setups/", & ! Path to correlation files
      cloud_file_ext  = "_corr_array_cloud.in", & ! File extensions for correlation files
      below_file_ext  = "_corr_array_below.in"

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
    real( kind = core_rknd ), dimension( res, res, res, res, res ) :: &
      droplets, droplets2            ! Used for lookup tables with GFDL activation

    character(len=128) :: corr_file_path_cloud, corr_file_path_below

    namelist /microphysics_setting/ &
      micro_scheme, l_cloud_sed, sigma_g, &
      l_ice_micro, l_graupel, l_hail, l_var_covar_src, l_upwind_diff_sed, &
      l_seifert_beheng, l_predictnc, l_const_Nc_in_cloud, specify_aerosol, &
      l_subgrid_w, l_arctic_nucl, l_cloud_edge_activation, l_fix_pgam, &
      l_in_cloud_Nc_diff, LH_microphys_type, l_local_kk, LH_microphys_calls, &
      LH_sequence_length, LH_seed, l_lh_cloud_weighted_sampling, &
      l_fix_s_t_correlations, l_lh_vert_overlap, l_silhs_KK_convergence_adj_mean, &
      rrp2_on_rrm2_cloud, Nrp2_on_Nrm2_cloud, Ncnp2_on_Ncnm2_cloud, &
      rrp2_on_rrm2_below, Nrp2_on_Nrm2_below, &
      Ncnp2_on_Ncnm2_below, C_evap, r_0, microphys_start_time, &
      Ncnm_initial, Nc0_in_cloud, ccnconst, ccnexpnt, aer_rm1, aer_rm2, &
      aer_n1, aer_n2, aer_sig1, aer_sig2, pgam_fixed

    namelist /gfdl_activation_setting/ &
      nooc, sul_concen, low_concen, high_concen, &
      lowup, highup, lowup2, highup2, lowmass2, &
      highmass2, lowmass3, highmass3,  &
      lowmass4, highmass4, lowmass5, highmass5, &
      lowT2, highT2, aeromass_value, l_gfdl_activation

    ! ---- Begin Code ----

    ! Set default values, then read in the namelist.
    ! Note: many parameters are set in the microphys_parameters module.

    l_gfdl_activation = .false.

    ! Cloud water sedimentation from the RF02 case
    ! This has no effect on Morrison's cloud water sedimentation
    l_cloud_sed = .false.
    sigma_g = 1.5_core_rknd ! Parameter for cloud droplet sedimentation code (RF02 value)

    !--------------------------------------------------------------------------
    ! Parameters for NNUCCD & NNUCCC coefficients on clex9_oct14 case
    !--------------------------------------------------------------------------

    select case (trim( runtype ))
    case ( "clex9_oct14")
      NNUCCD_REDUCE_COEF = .01 ! Reduce NNUCCD by factor of 100 for clex9_oct14
      NNUCCC_REDUCE_COEF = .01 ! Reduce NNUCCC by factor of 100 for clex9_oct14
    end select
    ! end change by Marc Pilon 11/16/11

    ! Aerosol for RF02 from Mikhail Ovtchinnikov
    aer_rm1  = 0.011e-6 ! Mean geometric radius  [m]
    aer_rm2  = 0.06e-6

    aer_sig1 = 1.2   ! Std dev of aerosol size distribution  [-]
    aer_sig2 = 1.7

    aer_n1   = 125.e6  ! Aerosol contentration                 [#/m3]
    aer_n2   = 65.e6

    ! Other parameters, set as in SAM
    ccnconst = 120. ! Parameter for powerlaw CCN [#/cm3]
    ccnexpnt = 0.4

    pgam_fixed = 5.

    LH_microphys_type = "disabled"

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
      call write_text ( "sigma_g = ", sigma_g, l_write_to_file, iunit )
      call write_text ( "l_graupel = ", l_graupel, l_write_to_file, iunit )
      call write_text ( "l_hail = ", l_hail, l_write_to_file, iunit )
      call write_text ( "l_seifert_beheng = ", l_seifert_beheng, &
        l_write_to_file, iunit )
      call write_text ( "l_predictnc = ", l_predictnc, l_write_to_file, iunit )
      call write_text ( "l_const_Nc_in_cloud = ", l_const_Nc_in_cloud, &
        l_write_to_file, iunit )
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
      call write_text ( "l_var_covar_src = ", l_var_covar_src, &
        l_write_to_file, iunit )
      call write_text ( "l_upwind_diff_sed = ", l_upwind_diff_sed, &
        l_write_to_file, iunit )
      call write_text ( "LH_microphys_type = " // &
        trim( LH_microphys_type ), l_write_to_file, iunit )
      call write_text ( "LH_microphys_calls = ", LH_microphys_calls, &
        l_write_to_file, iunit )
      call write_text ( "LH_sequence_length = ", LH_sequence_length, &
        l_write_to_file, iunit )
      call write_text ( "LH_seed = ", LH_seed, l_write_to_file, iunit )
      call write_text ( "l_lh_cloud_weighted_sampling = ", &
        l_lh_cloud_weighted_sampling, l_write_to_file, iunit )
      call write_text ( "l_fix_s_t_correlations = ", l_fix_s_t_correlations, &
        l_write_to_file, iunit )
      call write_text ( "l_silhs_KK_convergence_adj_mean = ", l_silhs_KK_convergence_adj_mean, &
        l_write_to_file, iunit)
      call write_text ( "l_lh_vert_overlap = ", l_lh_vert_overlap, &
        l_write_to_file, iunit )
      call write_text ( "rrp2_on_rrm2_cloud = ", rrp2_on_rrm2_cloud, &
        l_write_to_file, iunit )
      call write_text ( "Nrp2_on_Nrm2_cloud = ", Nrp2_on_Nrm2_cloud, &
        l_write_to_file, iunit )
      call write_text ( "Ncnp2_on_Ncnm2_cloud = ", Ncnp2_on_Ncnm2_cloud, &
        l_write_to_file, iunit )
      call write_text ( "rrp2_on_rrm2_below = ", rrp2_on_rrm2_below, &
        l_write_to_file, iunit )
      call write_text ( "Nrp2_on_Nrm2_below = ", Nrp2_on_Nrm2_below, &
        l_write_to_file, iunit )
      call write_text ( "Ncnp2_on_Ncnm2_below = ", Ncnp2_on_Ncnm2_below, &
        l_write_to_file, iunit )
      call write_text ( "C_evap = ", C_evap, l_write_to_file, iunit )
      call write_text ( "r_0 = ", r_0, l_write_to_file, iunit )
      call write_text ( "microphys_start_time = ", real( microphys_start_time, kind = core_rknd ), &
        l_write_to_file, iunit )
      call write_text ( "Ncnm_initial = ", Ncnm_initial, l_write_to_file, iunit )
      call write_text ( "Nc0_in_cloud = ", Nc0_in_cloud, l_write_to_file, iunit )
      call write_text ( "ccnconst = ", real(ccnconst, kind = core_rknd), l_write_to_file, iunit )
      call write_text ( "ccnexpnt = ", real(ccnexpnt, kind = core_rknd), l_write_to_file, iunit )
      call write_text ( "aer_rm1 = ", real(aer_rm1, kind = core_rknd), l_write_to_file, iunit )
      call write_text ( "aer_rm2 = ", real(aer_rm2, kind = core_rknd), l_write_to_file, iunit )
      call write_text ( "aer_n1 = ", real(aer_n1, kind = core_rknd), l_write_to_file, iunit )
      call write_text ( "aer_n2 = ", real(aer_n2, kind = core_rknd), l_write_to_file, iunit )
      call write_text ( "aer_sig1 = ", real(aer_sig1, kind = core_rknd), l_write_to_file, iunit )
      call write_text ( "aer_sig2 = ", real(aer_sig2, kind = core_rknd), l_write_to_file, iunit )
      call write_text ( "pgam_fixed = ", real(pgam_fixed, kind = core_rknd),l_write_to_file,iunit)

      if ( l_write_to_file ) close(unit=iunit)

    end if ! clubb_at_least_debug_level(1)

    ! Read in the name list for initialization, if it exists
    open(unit=iunit, file=namelist_file, status='old', action='read')
    read(iunit, nml=gfdl_activation_setting)
    close(unit=iunit)

    ! Initialize the GFDL activation code, if necessary
    if( l_gfdl_activation ) then
      ! Ensure a microphysics that has Ncm is being used
      if( trim( micro_scheme ) == "coamps" .or. trim( micro_scheme ) == "morrison" & 
            .or. trim( micro_scheme ) == "morrison_gettelman") then

        ! Read in the lookup tables
        call Loading( droplets, droplets2 )

        ! Initialize the activation variables
        call aer_ccn_act_k_init                            &
               ( real( droplets ), real( droplets2 ), res, res2, nooc,     &
                 real( sul_concen ), real( low_concen ), real( high_concen ),      &
                 real( lowup ), real( highup ), real( lowup2 ), real( highup2 ), real( lowmass2 ), &
                 real( highmass2 ), real( lowmass3 ), real( highmass3 ),           &
                 real( lowmass4 ), real( highmass4 ), real( lowmass5 ), real( highmass5 ), &
                 real( lowT2 ),real( highT2 ) )
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
      if ( l_predictnc ) then
        iiNcm       = 9
        hydromet_dim = 9
      else
        iiNcm        = -1
        hydromet_dim = 8
      end if


      allocate( hydromet_list(hydromet_dim) )

      hydromet_list(iirrainm)    = "rrainm"
      hydromet_list(iirsnowm)    = "rsnowm"
      hydromet_list(iiricem)     = "ricem"
      hydromet_list(iirgraupelm) = "rgraupelm"

      hydromet_list(iiNrm)       = "Nrm"
      hydromet_list(iiNsnowm)    = "Nsnowm"
      hydromet_list(iiNim)       = "Nim"
      hydromet_list(iiNgraupelm) = "Ngraupelm"
      if ( l_predictnc ) then
        hydromet_list(iiNcm)     = "Ncm"
      endif

      ! Set Nc0 in the Morrison code (module_MP_graupel) based on Nc0_in_cloud
      Nc0 = real( Nc0_in_cloud / cm3_per_m3 ) ! Units on Nc0 are per cm^3

      ! Set flags from the Morrison scheme as in GRAUPEL_INIT
      if ( l_predictnc ) then
        dopredictNc = .true.
      else
        dopredictNc = .false.
      end if

      ! Set the mode of aerosol to be used
      select case ( trim( specify_aerosol ) )
      case ( "morrison_no_aerosol" )
        aerosol_mode = morrison_no_aerosol

      case ( "morrison_power_law" )
        aerosol_mode = morrison_power_law

      case ( "morrison_lognormal" )
        aerosol_mode = morrison_lognormal

      case default
        stop "Unknown Morrison aerosol mode."

      end select

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
        write(fstderr,*) "Therefore l_ice_micro must be set to false to use this option."
        stop "Fatal error."
      end if
      allocate( l_hydromet_sed(hydromet_dim) )
      ! Sedimentation is handled within the Morrison microphysics
      l_hydromet_sed(iiNrm) = .false.
      l_hydromet_sed(iiNim) = .false.
      if ( l_predictnc ) then
        l_hydromet_sed(iiNcm) = .false.
      end if
      l_hydromet_sed(iiNgraupelm) = .false.
      l_hydromet_sed(iiNsnowm) = .false.

      l_hydromet_sed(iirrainm)    = .false.
      l_hydromet_sed(iirsnowm)    = .false.
      l_hydromet_sed(iiricem)     = .false.
      l_hydromet_sed(iirgraupelm) = .false.

      call GRAUPEL_INIT()

    case ( "morrison_gettelman" )
      iirrainm    = -1
      iirsnowm    = -1
      iiricem     = 1
      iirgraupelm = -1

      iiNrm       = -1
      iiNsnowm    = -1
      iiNim       = 2
      iiNgraupelm = -1
      if ( l_predictnc ) then
        write(fstderr,*) "Morrison-Gettelman microphysics is not currently configured for "// &
          "l_predictnc = T"
        stop "Fatal error."
        iiNcm       = 3
        hydromet_dim = 3
      else
        iiNcm       = -1
        hydromet_dim = 2
      end if

      if ( l_cloud_sed ) then
        write(fstderr,*) "Morrison-Gettelman microphysics has seperate code for cloud water"
        write(fstderr,*) "sedimentation, therefore l_cloud_sed should be set to .false."
        stop "Fatal error."
      endif

      allocate( hydromet_list(hydromet_dim) )

      hydromet_list(iiricem) = "ricem"
      hydromet_list(iiNim)   = "Nim"
      if ( l_predictnc ) then
        hydromet_list(iiNcm) = "Ncm"
      end if

      allocate( l_hydromet_sed(hydromet_dim) )
      ! Sedimentation is handled within the MG microphysics
      l_hydromet_sed(iiricem) = .false.
      l_hydromet_sed(iiNim)   = .false.
      if ( l_predictnc ) then
        l_hydromet_sed(iiNcm) = .false.
      endif

      ! Initialize constants for aerosols
      call ini_microp_aero()

      ! Setup the MG scheme
      call ini_micro()
      call pbuf_init()

      if ( l_use_CLUBB_pdf_in_mg ) then

         corr_file_path_cloud = corr_input_path//trim( runtype )//cloud_file_ext
         corr_file_path_below = corr_input_path//trim( runtype )//below_file_ext

         call setup_KK_corr( iunit, corr_file_path_cloud, corr_file_path_below, & ! In
                             l_write_to_file, case_info_file ) ! In

      endif

    case ( "coamps" )
      if ( .not. l_predictnc ) then
        write(fstderr,*) "COAMPS microphysics does not support l_predictnc = F"
        stop "Fatal Error"
      end if
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
      if ( l_predictnc ) then
        write(fstderr,*) "Khairoutdinov-Kogan microphysics does not support l_predictnc = T"
        stop "Fatal Error"
      end if
      iirrainm    = 1
      iirsnowm    = -1
      iiricem     = -1
      iirgraupelm = -1

      iiNrm       = 2
      iiNsnowm    = -1
      iiNim       = -1
      iiNgraupelm = -1
      iiNcm       = -1

      hydromet_dim = 2

      allocate( l_hydromet_sed(hydromet_dim) )

      allocate( hydromet_list(hydromet_dim) )

      hydromet_list(iirrainm) = "rrainm"
      hydromet_list(iiNrm) = "Nrm"

      l_hydromet_sed(iirrainm) = .true.
      l_hydromet_sed(iiNrm)    = .true.

      corr_file_path_cloud = corr_input_path//trim( runtype )//cloud_file_ext
      corr_file_path_below = corr_input_path//trim( runtype )//below_file_ext

      call setup_KK_corr( iunit, corr_file_path_cloud, corr_file_path_below, & ! In
                          l_write_to_file, case_info_file ) ! In

    case ( "simplified_ice", "none" )
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

    case default
      write(fstderr,*) "Unknown micro_scheme: "// trim( micro_scheme )
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
    
    ! Make sure the user didn't select LH sampling using
    ! coamps, morrison-gettelman, or simplified_ice microphysics
    ! (ticket:539)
    if ( (.not. (LH_microphys_type_int == LH_microphys_disabled)) .and. ( &
      trim( micro_scheme ) == "coamps" .or. &
      trim( micro_scheme ) == "morrison_gettelman" .or. &
      trim( micro_scheme ) == "simplified_ice" &
    ) ) then
      stop "LH sampling can not be enabled when using coamps, morrison_gettelman, or &
           &simplified_ice microphysics types"
    end if

    ! Make sure user hasn't selected l_silhs_KK_convergence_adj_mean when SILHS is disabled
    ! (ticket:558)
    if ( l_silhs_KK_convergence_adj_mean .and. &
      LH_microphys_type_int == LH_microphys_disabled) then
      stop "l_silhs_KK_convergence_adj_mean requires LH microphysics to be enabled."
    end if

    ! Make sure user hasn't selected l_silhs_KK_convergence_adj_mean when using
    ! a microphysics scheme other than khairoutdinov_kogan (KK)
    if ( l_silhs_KK_convergence_adj_mean .and. &
          trim( micro_scheme ) /= "khairoutdinov_kogan" ) then
      stop "l_silhs_KK_convergence_adj_mean requires khairoutdinov_kogan microphysics"
    end if

    ! Stop the run if user has selected l_silhs_KK_convergence_adj_mean and not
    ! included the preprocessor flag SILHS_KK_CONVERGENCE_TEST
#ifndef SILHS_KK_CONVERGENCE_TEST
    if (l_silhs_KK_convergence_adj_mean) then
      stop "Use of the l_silhs_KK_convergence_adj_mean flag requires the &
           &preprocessor flag SILHS_KK_CONVERGENCE_TEST"
    end if
#endif

    ! Setup index variables for latin hypercube sampling
    if ( LH_microphys_type_int /= LH_microphys_disabled ) then

#ifdef LATIN_HYPERCUBE
      corr_file_path_cloud = corr_input_path//trim( runtype )//cloud_file_ext
      corr_file_path_below = corr_input_path//trim( runtype )//below_file_ext
      ! Allocate and set the arrays containing the correlations
      ! and the X'^2 / X'^2 terms
      call setup_corr_varnce_array( iirrainm, iiNrm, iiricem, iiNim, iirsnowm, iiNsnowm, &
                                    l_ice_micro, &
                                    corr_file_path_cloud, corr_file_path_below, iunit )
#endif

    end if

    return
  end subroutine init_microphys

!-------------------------------------------------------------------------------
  subroutine advance_microphys & 
             ( iter, runtype, dt, time_current,                         & ! Intent(in)
               thlm, p_in_Pa, exner, rho, rho_zm, rtm, rcm, cloud_frac, & ! Intent(in)
               wm_zt, wm_zm, Kh_zm, pdf_params,                         & ! Intent(in)
               wp2_zt, rho_ds_zt, rho_ds_zm, invrs_rho_ds_zt,           & ! Intent(in)
               LH_sample_point_weights,                                 & ! Intent(in)
               X_nl_all_levs, X_mixt_comp_all_levs, LH_rt, LH_thl,      & ! Intent(in)
               n_variables, corr_array_1, corr_array_2,                 & ! Intent(in)
               mu_x_1, mu_x_2, sigma_x_1, sigma_x_2,                    & ! Intent(in)
               hydromet_pdf_params,                                     & ! Intent(in)
               Ncnm, hydromet, wphydrometp,                             & ! Intent(inout)
               rvm_mc, rcm_mc, thlm_mc,                                 & ! Intent(inout)
               wprtp_mc_tndcy, wpthlp_mc_tndcy,                         & ! Intent(inout)
               rtp2_mc_tndcy, thlp2_mc_tndcy, rtpthlp_mc_tndcy,         & ! Intent(inout)
               err_code )                                                 ! Intent(out)

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
        KK_local_micro_driver, &  ! Procedure(s)
        KK_upscaled_micro_driver

    use cloud_sed_module, only: cloud_drop_sed ! Procedure(s)

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

    use advance_windm_edsclrm_module, only : &
        xpwp_fnc  ! Procedure(s)

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
        one, &
        one_half, &
        zero, &
        zero_threshold, &
        sec_per_day, &
        cm3_per_m3, &
        mm_per_m

    use model_flags, only: &
        l_hole_fill, & ! Variable(s)
        l_evaporate_cold_rcm, &
        l_morr_xp2_mc_tndcy

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        dp, &
        core_rknd

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

    use hydromet_pdf_parameter_module, only:  &
        hydromet_pdf_parameter  ! Type

    use array_index, only:  & 
        iirrainm, iirsnowm, iiricem, iirgraupelm, & ! Variable(s)
        iiNrm, iiNim, iiNcm

    use stats_variables, only: & 
        iVrr,  & ! Variable(s)
        iVNr, & 
        iVrsnow, & 
        iVrice, & 
        iVrgraupel, &
        iwprrp, &
        iwprip, &
        iwprsp, &
        iwprgp, &
        iwpNrp, &
        iwpNip, &
        iwpNsp, &
        iwpNgp, &
        iwpNcp, &
        iVrrprrp, &
        iVNrpNrp, & 
        irain_rate_zt, & 
        iFprec

    use stats_variables, only: & 
        irrainm_bt,       & ! Variable(s)
        irrainm_mc,       & 
        irrainm_hf,       & 
        irrainm_wvhf,     &
        irrainm_cl,       &
        iricem_bt,        & 
        iricem_mc,        &
        iricem_hf,        &
        iricem_wvhf,      &
        iricem_cl,        &
        irgraupelm_bt,    &
        irgraupelm_mc,    &
        irgraupelm_hf,    &
        irgraupelm_wvhf,  &
        irgraupelm_cl,    &
        irsnowm_bt,       & 
        irsnowm_mc,       &
        irsnowm_hf,       &
        irsnowm_wvhf,     &
        irsnowm_cl,       &
        irain_rate_sfc,   & 
        irain_flux_sfc,   & 
        irrainm_sfc

    use stats_variables, only: & 
        iNrm_bt,       & ! Variable(s)
        iNrm_mc,       &
        iNrm_cl,       &
        iNim_bt,       &
        iNim_cl,       &
        iNim_mc,       &
        iNsnowm_bt,    &
        iNsnowm_cl,    &
        iNsnowm_mc,    &
        iNgraupelm_bt, &
        iNgraupelm_cl, &
        iNgraupelm_mc, &
        iNc_in_cloud, &
        iNc_activated, &
        iNcnm,         &
        iNcm_bt,       &
        iNcm_cl,       &
        iNcm_mc,       &
        iNcm_act,      &
        iNcm

    use stats_subs, only: & 
        stats_accumulate_LH_tend ! Procedure(s)

    use stats_variables, only: & 
        LH_zt, & ! Variable(s)
        iLH_Vrr, &
        iLH_VNr

    use stats_variables, only: & 
        zt, &  ! Variables
        zm, & 
        sfc, & 
        l_stats_samp

    use stats_type, only: & 
        stat_update_var,      & ! Procedure(s)
        stat_update_var_pt,   &
        stat_begin_update,    &
        stat_begin_update_pt, &
        stat_end_update,      &
        stat_end_update_pt

    use stats_subs, only: &
        stats_accumulate_hydromet

    use fill_holes, only: &
        vertical_avg, & ! Procedure(s)
        fill_holes_driver

    use phys_buffer, only: & ! Used for placing wp2_zt in morrison_gettelman microphysics
        pbuf_add,            &
        pbuf_allocate,       &
        pbuf_setval

    use shr_kind_mod, only: r8 => shr_kind_r8

    use parameters_microphys, only: &
        LH_microphys_type, & ! Determines how the LH samples are used
        LH_microphys_interactive,     & ! Feed the subcols into the microphys and allow feedback
        LH_microphys_non_interactive, & ! Feed the subcols into the microphys with no feedback
        LH_microphys_disabled           ! Disable latin hypercube entirely

    use clubb_precision, only: &
        core_rknd ! Variable(s)

#ifdef LATIN_HYPERCUBE
    use latin_hypercube_arrays, only: &
        d_variables ! Variable(s)
#else
#define d_variables 0
#endif /* LATIN_HYPERCUBE */

    implicit none

    ! Input Variables

    integer, intent(in) :: &
    iter,       & ! Model iteration number
    n_variables   ! Number of variables in the correlation arrays

    character(len=*), intent(in) :: & 
      runtype ! Name of the run, for case specific effects.

    real(kind=time_precision), intent(in) ::  & 
      dt           ! Timestep         [s]

    real(kind=time_precision), intent(in) ::  & 
      time_current ! Current time     [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      thlm,       & ! Liquid potential temp.                    [K]
      p_in_Pa,    & ! Pressure                                  [Pa]
      exner,      & ! Exner function                            [-]
      rho,        & ! Density on thermodynamic levels           [kg/m^3]
      rho_zm,     & ! Density on momentum levels                [kg/m^3]
      rtm,        & ! Total water mixing ratio                  [kg/kg]
      rcm,        & ! Liquid water mixing ratio                 [kg/kg]
      cloud_frac, & ! Cloud fraction                            [-]
      wm_zt,      & ! w wind component on thermodynamic levels  [m/s]
      wm_zm,      & ! w wind component on momentum levels       [m/s]
      Kh_zm         ! Kh Eddy diffusivity on momentum grid      [m^2/s]

    type(pdf_parameter), dimension(gr%nz), intent(in) :: & 
      pdf_params     ! PDF parameters

    type(hydromet_pdf_parameter), dimension(gr%nz), intent(in) :: &
      hydromet_pdf_params     ! PDF parameters

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      wp2_zt,          & ! w'^2 on the thermo. grid                 [m^2/s^2]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = dp ), dimension(gr%nz,LH_microphys_calls,d_variables), intent(in) :: &
      X_nl_all_levs ! Lognormally distributed hydrometeors

    integer, dimension(gr%nz,LH_microphys_calls), intent(in) :: &
      X_mixt_comp_all_levs ! Which mixture component the sample is in

    real( kind = core_rknd ), dimension(gr%nz,LH_microphys_calls), intent(in) :: &
      LH_rt, LH_thl ! Samples of rt, thl        [kg/kg,K]

    real( kind = core_rknd ), dimension(LH_microphys_calls), intent(in) :: &
      LH_sample_point_weights ! Weights for cloud weighted sampling

    real( kind = core_rknd ), dimension(n_variables, n_variables, gr%nz), intent(in) :: &
      corr_array_1, & ! Correlation matrix for the first pdf component    [-]
      corr_array_2    ! Correlation matrix for the second pdf component   [-]

    real( kind = core_rknd ), dimension(n_variables, gr%nz), intent(in) :: &
      mu_x_1,    & ! Mean array for the 1st PDF component                 [units vary]
      mu_x_2,    & ! Mean array for the 2nd PDF component                 [units vary]
      sigma_x_1, & ! Standard deviation array for the 1st PDF component   [units vary]
      sigma_x_2    ! Standard deviation array for the 2nd PDF component   [units vary]

    ! Input/Output Variables

    ! Note:
    ! For COAMPS Ncnm is initialized and Nim & Ncm are computed within
    ! subroutine adjtg.
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: & 
      Ncnm    ! Cloud nuclei concentration     [num/m^3]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(inout) :: &
      hydromet,    & ! Hydrometeor mean, < h_m > (thermodynamic levels)  [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)      [(m/s)units]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: & 
      rcm_mc,  & ! Microphysics contributions to liquid water           [kg/kg/s]
      rvm_mc,  & ! Microphysics contributions to vapor water            [kg/kg/s]
      thlm_mc    ! Microphysics contributions to liquid potential temp. [K/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      wprtp_mc_tndcy,   & ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc_tndcy,  & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtp2_mc_tndcy,    & ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc_tndcy,   & ! Microphysics tendency for <thl'^2>  [K^2/s]
      rtpthlp_mc_tndcy    ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]

    ! Output Variables

    integer, intent(out) :: err_code ! Exit code returned from subroutine

    ! Local Variables

    real( kind = core_rknd ), dimension(3,gr%nz) :: & 
      lhs    ! Left hand side of tridiagonal matrix

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      rhs    ! Right hand side vector

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim) :: &
      hydromet_vel,    & ! Contains vel. of the hydrometeors        [m/s]
      hydromet_vel_zt

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim) :: & 
      hydromet_mc  ! Change in hydrometeors due to microphysics  [units/s]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim) :: &
      hydromet_vel_covar,    & ! Covariance of V_xx & x_x (m-levs)  [units(m/s)]
      hydromet_vel_covar_zt    ! Covariance of V_xx & x_x (t-levs)  [units(m/s)]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim) :: &
      hydromet_vel_covar_zt_impc, & ! Imp. comp. <V_xx'x_x'> t-levs [m/s]
      hydromet_vel_covar_zt_expc    ! Exp. comp. <V_xx'x_x'> t-levs [units(m/s)]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      delta_zt  ! Difference in thermo. height levels     [m]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      T_in_K, & ! Temperature              [K]
      Ncm,    & ! Mean cloud droplet number concentration [#/kg]
      rvm,    & ! Vapor water mixing ratio [kg/kg]
      thlm_morr, & ! Thlm fed into morrison microphysics [K]
      rcm_morr, &  ! rcm fed into morrison microphysics [kg/kg]
      cloud_frac_morr  ! Cloud fraction fed into morrision microphysics []

    real( kind = core_rknd ), dimension(gr%nz) :: & 
      rrainm_auto, & ! Autoconversion rate for rrainm      [kg/kg/s]
      rrainm_accr, & ! Accretion rate for rrainm           [kg/kg/s]
      rrainm_evap, & ! Evaporation rate for rrainm         [kg/kg/s]
      Nrm_auto,    & ! Change in Nrm due to autoconversion [num/kg/s]
      Nrm_evap       ! Change in Nrm due to evaporation    [num/kg/s]

    real( kind = core_rknd ), dimension(1,1,gr%nz) :: & 
      cond ! COAMPS stat for condesation/evap of rcm

    ! Eddy diffusivity for rain and rain drop concentration.
    ! It is also used for the other hydrometeor variables.
    ! Kr = Constant * Kh_zm; Constant is named c_Krrainm.
    real( kind = core_rknd ), dimension(gr%nz) :: Kr   ! [m^2/s]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      wtmp,    & ! Standard dev. of w                   [m/s]
      s_mellor   ! The variable 's' in Mellor (1977)    [kg/kg]

    real( kind = core_rknd ) :: &
      max_velocity, & ! Maximum sedimentation velocity [m/s]
      temp            ! Temporary variables     [units vary]

    integer :: i, k ! Loop iterators / Array indices

    integer :: ixrm_hf, ixrm_wvhf, ixrm_cl, &
               ixrm_bt, ixrm_mc, iwpxrp

    real( kind = core_rknd ), dimension(gr%nz) :: &
      Ndrop_max  ! GFDL droplet activation concentration [#/kg]

    !Input aerosol mass concentration: the unit is 10^12 ug/m3.
    !For example, aeromass=2.25e-12 means that the aerosol mass
    !concentration is 2.25 ug/m3.
    !This value of aeromass was recommended by Huan Guo
    !See http://carson.math.uwm.edu/trac/climate_process_team/ticket/46#comment:12
    real( kind = core_rknd ), dimension(gr%nz, 4) :: &
      aeromass ! ug/m^3

    real( kind = core_rknd ), dimension(gr%nz) :: &
      Nc_in_cloud ! cloud average droplet concentration [#/kg]

    character(len=10) :: hydromet_name

    logical :: l_latin_hypercube_input

!-------------------------------------------------------------------------------

    ! ---- Begin code ----

    Ndrop_max = zero
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

    ! Initialize the values of the implicit and explicit components used to 
    ! calculate the covariances of hydrometeor sedimentation velocities and
    ! their associated hydrometeors (for example, <V_rr'r_r'> and <V_Nr'N_r'>)
    ! to 0.
    if ( hydromet_dim > 0 ) then
       hydromet_vel_covar_zt_impc = zero
       hydromet_vel_covar_zt_expc = zero
    endif

    ! Solve for the value of Kr, the hydrometeor eddy diffusivity.
    do k = 1, gr%nz, 1
      Kr(k) = c_Krrainm * Kh_zm(k)
    end do

    ! Determine 's' from Mellor (1977)
    s_mellor(:) = pdf_params(:)%mixt_frac * pdf_params(:)%s1 &
                  + ( one - pdf_params(:)%mixt_frac ) * pdf_params(:)%s2

    ! Compute standard deviation of vertical velocity in the grid column
    wtmp(:) = sqrt( wp2_zt(:) )

    ! Compute difference in thermodynamic height levels
    delta_zt(1:gr%nz) = one / gr%invrs_dzm(1:gr%nz)

    ! Calculate T_in_K
    T_in_K = thlm2T_in_K( thlm, exner, rcm )

    ! Start Ncm budget to capture Ncm_act term
    if ( l_stats_samp .and. iiNcm > 0 ) then
       if ( l_hydromet_sed(iiNcm) .and. l_upwind_diff_sed ) then
          call stat_begin_update( iNcm_bt, &
                                   hydromet(:,iiNcm) &
                                   / real( dt, kind = core_rknd ), zt )
       else
          ! The mean hydrometeor field at thermodynamic level k = 1 is simply
          ! set to the value of the mean hydrometeor field at k = 2.  Don't
          ! include budget stats for level k = 1.
          do k = 2, gr%nz, 1
             call stat_begin_update_pt( iNcm_bt, k, &
                                        hydromet(k,iiNcm) &
                                        / real( dt, kind = core_rknd ), zt )
          enddo
       endif
    endif ! l_stats_samp and iiNcm > 0

    ! Call GFDL activation code
    if( l_gfdl_activation ) then
      ! Ensure a microphysics that has Ncm is being used
      ! Note: KK micro is not used because it resets Ncm every timestep
      if ( l_predictnc ) then

        ! Save the initial Ncm value for the Ncm_act term
        if ( l_stats_samp ) then
          call stat_begin_update( iNcm_act, hydromet(:, iiNcm) / real( dt, kind = core_rknd ), zt )
        end if

        call aer_act_clubb_quadrature_Gauss( aeromass, T_in_K, Ndrop_max )

        ! Convert to #/kg
        Ndrop_max = Ndrop_max * cm3_per_m3 / rho

        if( l_stats_samp ) then
          call stat_update_var( iNc_activated, Ndrop_max, zt)
        end if

        ! Clip Ncm values that are outside of cloud by CLUBB standards
        do k=1, gr%nz

! ---> h1g, 2011-04-20,   no liquid drop nucleation if T < -40 C
          if( T_in_K(k) <= 233.15_core_rknd )  Ndrop_max(k) = &
             zero  ! if T<-40C, no liquid drop nucleation
! <--- h1g, 2011-04-20

          ! Clip Ncm to only be where there is cloud to avoid divide by zero errors
          if( cloud_frac(k) > cloud_frac_min ) then
            hydromet(k, iiNcm) = max( Ndrop_max(k), hydromet(k, iiNcm) )
          else
            hydromet(k, iiNcm) = zero
          end if
        end do

        ! Update the Ncm_act term
        if( l_stats_samp ) then
          call stat_end_update( iNcm_act, hydromet(:,iiNcm) / real( dt, kind = core_rknd ), zt )
        end if

      else

        stop "Unsupported microphysics scheme for GFDL activation."

      end if  ! l_predictnc
    end if ! l_gfdl_activation

    ! Determine how Ncm and Nc_in_cloud will be computed
    select case ( trim( micro_scheme ) )

    case ( "coamps", "morrison", "morrison_gettelman", "khairoutdinov_kogan" )

       if ( l_predictnc ) then

          Ncm = hydromet(:,iiNcm)

       else ! Compute the fixed value by multiplying by cloud fraction

          if ( .not. l_const_Nc_in_cloud ) then
             where ( rcm >= rc_tol )
                Ncm = cloud_frac * ( Nc0_in_cloud / rho ) ! Convert to #/kg air
             else where
                Ncm = zero
             end where

             Nc_in_cloud = Ncm / max( cloud_frac, cloud_frac_min )

          else  ! Constant, specified value of Nc within cloud

             Nc_in_cloud = Nc0_in_cloud / rho ! Convert to #/kg air.

             Ncm = cloud_frac * ( Nc0_in_cloud / rho ) ! Convert to #/kg air

          endif

      endif

    case default

      if ( l_cloud_sed ) then
        where ( rcm >= rc_tol )
          Ncm = cloud_frac * ( Nc0_in_cloud / rho ) ! Convert to #/kg air
        else where
          Ncm = zero
        end where
        Nc_in_cloud = Ncm / max( cloud_frac, cloud_frac_min )

      else
        ! These quantites are undefined
        Ncm = -999._core_rknd
        Nc_in_cloud = -999._core_rknd
      end if

    end select ! micro_scheme

    ! Begin by calling Brian Griffin's implementation of the
    ! Khairoutdinov and Kogan microphysics (analytic or local formulas),
    ! the COAMPS implementation of Rutlege and Hobbes, or the Morrison
    ! microphysics.
    ! Note: COAMPS appears to have some K&K elements to it as well.

    select case ( trim( micro_scheme ) )

    case ( "coamps" )

      ! Initialize tendencies to zero
      hydromet_mc(:,:) = zero
      hydromet_vel(:,:) = zero
      hydromet_vel_zt(:,:) = zero
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
      cond = -999._core_rknd
      if ( cond(1,1,1) /= cond(1,1,1) ) stop
#endif

      if ( l_stats_samp ) then

        ! Sedimentation velocity for rrainm
        call stat_update_var(iVrr, zt2zm( hydromet_vel_zt(:,iirrainm) ), zm)

        ! Sedimentation velocity for Nrm
        call stat_update_var(iVNr, zt2zm( hydromet_vel(:,iiNrm) ), zm )

        ! Sedimentation velocity for snow
        call stat_update_var(iVrsnow, zt2zm( hydromet_vel(:,iirsnowm) ), zm )

        ! Sedimentation velocity for pristine ice
        call stat_update_var( iVrice, zt2zm( hydromet_vel(:,iiricem) ), zm )

        ! Sedimentation velocity for graupel
        call stat_update_var( iVrgraupel, zt2zm( hydromet_vel(:,iirgraupelm) ), zm )
      end if ! l_stats_samp

    case ( "morrison" )

      ! Initialize tendencies to zero
      hydromet_mc(:,:) = zero
      rcm_mc(:) = zero
      rvm_mc(:) = zero
      thlm_mc(:) = zero

      rcm_morr(:) = rcm(:)
      cloud_frac_morr(:) = cloud_frac(:)

      if ( l_evaporate_cold_rcm  ) then
        ! Convert liquid to vapor at temperatures colder than -37C
        where ( T_in_K(:) < 236.15_core_rknd )
          rcm_morr(:) = 0.0_core_rknd
          cloud_frac_morr(:) = 0.0_core_rknd
          hydromet(:,iiNcm) = 0.0_core_rknd
        end where
      end if
 
      if ( LH_microphys_type /= LH_microphys_disabled ) then
#ifdef LATIN_HYPERCUBE
        call LH_microphys_driver &
             ( dt, gr%nz, LH_microphys_calls, d_variables, & ! In
               X_nl_all_levs, LH_rt, LH_thl, LH_sample_point_weights, & ! In
               pdf_params, p_in_Pa, exner, rho, & ! In
               rcm_morr, wtmp, delta_zt, cloud_frac_morr, & ! In
               hydromet, X_mixt_comp_all_levs, Nc_in_cloud, & !In 
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
#endif /* LATIN_HYPERCUBE */
        call stats_accumulate_LH_tend( hydromet_mc, thlm_mc, rvm_mc, rcm_mc )

      end if ! LH isn't disabled

      ! Call the microphysics if we don't want to have feedback effects from the
      ! latin hypercube result (above)
      if ( LH_microphys_type /= LH_microphys_interactive ) then
        l_latin_hypercube_input = .false.
        
        if ( l_morr_xp2_mc_tndcy ) then
          !Use the moister rt1/rt2 rather than rtm in morrison micro
          !Also use the colder of thl1/thl2
          where ( pdf_params%rt1 > pdf_params%rt2 )
            rvm = pdf_params%rt1 - pdf_params%rc1
          else where
            rvm = pdf_params%rt2 - pdf_params%rc2
          end where

          where ( pdf_params%thl1 < pdf_params%thl2 )
            thlm_morr = pdf_params%thl1
          else where
            thlm_morr = pdf_params%thl2
          end where
        else
          rvm = rtm - rcm
          thlm_morr = thlm
        end if

        call morrison_micro_driver & 
             ( dt, gr%nz, l_stats_samp, &
               l_latin_hypercube_input, thlm_morr, wm_zt, p_in_Pa, &
               exner, rho, cloud_frac_morr, pdf_params, wtmp, &
               delta_zt, rcm_morr, Ncm, s_mellor, rvm, hydromet, &
               hydromet_mc, hydromet_vel_zt, &
               rcm_mc, rvm_mc, thlm_mc, &
               rtp2_mc_tndcy, thlp2_mc_tndcy, &
               wprtp_mc_tndcy, wpthlp_mc_tndcy, rtpthlp_mc_tndcy, &
               rrainm_auto, rrainm_accr, rrainm_evap, &
               Nrm_auto, Nrm_evap )

        ! Output rain sedimentation velocity
        if ( l_stats_samp ) then
          call stat_update_var(iVrr, zt2zm( hydromet_vel_zt(:,iirrainm) ), zm)
        end if

      end if ! LH_microphys_type /= interactive

    case ( "morrison_gettelman" )

      ! Initialize tendencies to zero
      hydromet_mc(:,:) = zero
      rcm_mc(:) = zero
      rvm_mc(:) = zero
      thlm_mc(:) = zero

      hydromet_vel = zero
      hydromet_vel_zt = zero

      ! Place wp2 into the dummy phys_buffer module to import it into microp_aero_ts.
      ! Placed here because parameters cannot be changed on mg_microphys_driver with
      ! the way LH is currently set up.
      call pbuf_add( 'WP2', 1, gr%nz, 1 )
      call pbuf_allocate()
      call pbuf_setval( 'WP2', real( wp2_zt, kind=r8 ) )

      rvm = rtm - rcm
      call mg_microphys_driver &
          ( dt, gr%nz, l_stats_samp, gr%invrs_dzt, thlm, p_in_Pa, exner, &
            rho, cloud_frac, rcm, Ncm, rvm, Ncnm, pdf_params, hydromet, &
            hydromet_mc, hydromet_vel_zt, rcm_mc, rvm_mc, thlm_mc )

    case ( "khairoutdinov_kogan" )

      ! Initialize tendencies to zero
      hydromet_mc(:,:) = zero
      rcm_mc(:) = zero
      rvm_mc(:) = zero
      thlm_mc(:) = zero

      if ( LH_microphys_type /= LH_microphys_disabled ) then

#ifdef LATIN_HYPERCUBE
        call LH_microphys_driver &
             ( dt, gr%nz, LH_microphys_calls, d_variables, & ! In
               X_nl_all_levs, LH_rt, LH_thl, LH_sample_point_weights, & ! In
               pdf_params, p_in_Pa, exner, rho, & ! In
               rcm, wtmp, delta_zt, cloud_frac, & ! In
               hydromet, X_mixt_comp_all_levs, Nc_in_cloud, & !In 
               hydromet_mc, hydromet_vel_zt, & ! In/Out
               rcm_mc, rvm_mc, thlm_mc,  & ! Out
               KK_local_micro_driver ) ! Procedure
#else
        stop "Latin hypercube was not enabled at compile time"
#endif /* LATIN_HYPERCUBE */

        call stats_accumulate_LH_tend( hydromet_mc, thlm_mc, rvm_mc, rcm_mc )

        if ( l_stats_samp ) then
          ! Latin hypercube estimate for sedimentation velocities
          call stat_update_var( iLH_Vrr, hydromet_vel_zt(:,iirrainm), LH_zt )

          call stat_update_var( iLH_VNr, hydromet_vel_zt(:,iiNrm), LH_zt )

        end if

      end if ! LH isn't disabled

      ! Call the microphysics if we don't want to have feedback effects from the
      ! latin hypercube result (above)
      if ( LH_microphys_type /= LH_microphys_interactive ) then

        l_latin_hypercube_input = .false.
        rvm = rtm - rcm

        if ( l_local_kk ) then

          call KK_local_micro_driver( dt, gr%nz, l_stats_samp, &
                                      l_latin_hypercube_input, thlm, wm_zt, &
                                      p_in_Pa, exner, rho, cloud_frac, &
                                      pdf_params, wtmp, delta_zt, rcm, &
                                      Ncm, s_mellor, rvm, &
                                      hydromet, &
                                      hydromet_mc, hydromet_vel_zt, &
                                      rcm_mc, rvm_mc, thlm_mc, &
                                      rtp2_mc_tndcy, thlp2_mc_tndcy, &
                                      wprtp_mc_tndcy, wpthlp_mc_tndcy, &
                                      rtpthlp_mc_tndcy,  &
                                      rrainm_auto, rrainm_accr, rrainm_evap, &
                                      Nrm_auto, Nrm_evap )

        else

          call KK_upscaled_micro_driver( dt, gr%nz, l_stats_samp, thlm, wm_zt,    & ! Intent(in)
                                         p_in_Pa, exner, rho, cloud_frac,         & ! Intent(in)
                                         pdf_params, wtmp, rcm, Ncnm,             & ! Intent(in)
                                         s_mellor, Nc_in_cloud,                   & ! Intent(in)
                                         hydromet, wphydrometp,                   & ! Intent(in)
                                         n_variables, corr_array_1, corr_array_2, & ! Intent(in)
                                         mu_x_1, mu_x_2, sigma_x_1, sigma_x_2,    & ! Intent(in)
                                         hydromet_pdf_params,                     & ! Intent(in)
                                         hydromet_mc, hydromet_vel_zt,            & ! Intent(inout)
                                         rcm_mc, rvm_mc, thlm_mc,                 & ! Intent(out)
                                         hydromet_vel_covar_zt_impc,              & ! Intent(out)
                                         hydromet_vel_covar_zt_expc,              & ! Intent(out)
                                         wprtp_mc_tndcy, wpthlp_mc_tndcy,         & ! Intent(out)
                                         rtp2_mc_tndcy, thlp2_mc_tndcy,           & ! Intent(out)
                                         rtpthlp_mc_tndcy )                         ! Intent(out)

        endif

      end if ! LH_microphys_type /= interactive

      if ( l_stats_samp ) then

        ! Sedimentation velocity for rrainm
        call stat_update_var( iVrr, zt2zm( hydromet_vel_zt(:,iirrainm) ), zm )

        ! Sedimentation velocity for Nrm
        call stat_update_var( iVNr, zt2zm( hydromet_vel_zt(:,iiNrm) ), zm )

      end if ! l_stats_samp

    case ( "simplified_ice" )

      ! Call the simplified ice diffusion scheme
      call ice_dfsn( dt, thlm, rcm, exner, p_in_Pa, rho, rcm_mc, thlm_mc )

    case default
      ! Do nothing

    end select ! micro_scheme


    ! Re-compute Ncm and Nc_in_cloud if needed
    if ( iiNcm > 0 ) then
      Ncm = hydromet(:,iiNcm)
      Nc_in_cloud = Ncm / max( cloud_frac, cloud_frac_min )
    end if


    !-----------------------------------------------------------------------
    !       Loop over all hydrometeor species and apply sedimentation,
    !       advection and diffusion.
    !-----------------------------------------------------------------------

    do i = 1, hydromet_dim

      ! Initializing max_velocity in order to avoid a compiler warning.
      ! Regardless of the case, it will be reset in the 'select case'
      ! statement immediately below.
      max_velocity = zero

      select case ( trim( hydromet_list(i) ) )
      case ( "rrainm" )
        ixrm_bt   = irrainm_bt
        ixrm_hf   = irrainm_hf
        ixrm_wvhf = irrainm_wvhf
        ixrm_cl   = irrainm_cl
        ixrm_mc   = irrainm_mc

        iwpxrp    = iwprrp

        max_velocity = -9.1_core_rknd ! m/s

      case ( "ricem" )
        ixrm_bt   = iricem_bt
        ixrm_hf   = iricem_hf
        ixrm_wvhf = iricem_wvhf
        ixrm_cl   = iricem_cl
        ixrm_mc   = iricem_mc

        iwpxrp    = iwprip

        max_velocity = -1.2_core_rknd ! m/s

      case ( "rsnowm" )
        ixrm_bt   = irsnowm_bt
        ixrm_hf   = irsnowm_hf
        ixrm_wvhf = irsnowm_wvhf
        ixrm_cl   = irsnowm_cl
        ixrm_mc   = irsnowm_mc

        iwpxrp    = iwprsp

        ! Morrison limit
!         max_velocity = -1.2_core_rknd ! m/s
        ! Made up limit.  The literature suggests that it is quite possible
        ! that snow flake might achieve a terminal velocity of 2 m/s, and this
        ! happens in the COAMPS microphysics -dschanen 29 Sept 2009
        max_velocity = -2.0_core_rknd ! m/s

      case ( "rgraupelm" )
        ixrm_bt   = irgraupelm_bt
        ixrm_hf   = irgraupelm_hf
        ixrm_wvhf = irgraupelm_wvhf
        ixrm_cl   = irgraupelm_cl
        ixrm_mc   = irgraupelm_mc

        iwpxrp    = iwprgp

        max_velocity = -20._core_rknd ! m/s

      case ( "Nrm" )
        ixrm_bt   = iNrm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = iNrm_cl
        ixrm_mc   = iNrm_mc

        iwpxrp    = iwpNrp

        max_velocity = -9.1_core_rknd ! m/s

      case ( "Nim" )
        ixrm_bt   = iNim_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = iNim_cl
        ixrm_mc   = iNim_mc

        iwpxrp    = iwpNip

        max_velocity = -1.2_core_rknd ! m/s

      case ( "Nsnowm" )
        ixrm_bt   = iNsnowm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = iNsnowm_cl
        ixrm_mc   = iNsnowm_mc

        iwpxrp    = iwpNsp

        ! Morrison limit
!         max_velocity = -1.2_core_rknd ! m/s
        ! Made up limit.  The literature suggests that it is quite possible
        ! that snow flake might achieve a terminal velocity of 2 m/s, and this
        ! happens in the COAMPS microphysics -dschanen 29 Sept 2009
        max_velocity = -2.0_core_rknd ! m/s

      case ( "Ngraupelm" )
        ixrm_bt   = iNgraupelm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = iNgraupelm_cl
        ixrm_mc   = iNgraupelm_mc

        iwpxrp    = iwpNgp

        max_velocity = -20._core_rknd ! m/s

      case ( "Ncm" )
        ixrm_bt   = iNcm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = iNcm_cl
        ixrm_mc   = iNcm_mc

        iwpxrp    = iwpNcp

        ! Use the rain water limit, since Morrison has no explicit limit on
        ! cloud water.  Presumably these numbers are never large.
        ! -dschanen 28 Sept 2009
        max_velocity = -9.1_core_rknd ! m/s

      case default
        ixrm_bt   = 0
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = 0
        ixrm_mc   = 0

        iwpxrp    = 0

        max_velocity = -9.1_core_rknd ! m/s

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
           if ( l_hydromet_sed(i) .and. l_upwind_diff_sed ) then
              call stat_begin_update( ixrm_bt, &
                                      hydromet(:,i) &
                                      / real( dt, kind = core_rknd ), zt )
           else
              ! The mean hydrometeor field at thermodynamic level k = 1 is
              ! simply set to the value of the mean hydrometeor field at k = 2.
              ! Don't include budget stats for level k = 1.
              do k = 2, gr%nz, 1
                 call stat_begin_update_pt( ixrm_bt, k, &
                                            hydromet(k,i) &
                                            / real( dt, kind = core_rknd ), zt )
              enddo
           endif
        endif

      endif ! l_stats_samp

      ! Set realistic limits on sedimentation velocities, following the
      ! numbers in the Morrison microphysics.

      do k = 1, gr%nz
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
      end do ! k = 1..gr%nz

      ! Interpolate velocity to the momentum grid for a centered difference
      ! approximation of the sedimenation term.
      hydromet_vel(:,i) = zt2zm( hydromet_vel_zt(:,i) )
      hydromet_vel(gr%nz,i) = zero ! Upper boundary condition

      ! Solve for < w'h_m' > at all intermediate (momentum) grid levels, using
      ! a down-gradient approximation:  < w'h_m' > = - K * d< h_m >/dz.
      ! A Crank-Nicholson time-stepping scheme is used for this variable.
      ! This is the portion of the calculation using < h_m > from timestep t. 
      wphydrometp(1:gr%nz-1,i) &
      = - one_half * xpwp_fnc( Kr(1:gr%nz-1)+nu_r_vert_res_dep(1:gr%nz-1), &
                               hydromet(1:gr%nz-1,i), hydromet(2:gr%nz,i), &
                               gr%invrs_dzm(1:gr%nz-1) )

      ! A zero-flux boundary condition at the top of the model is used for
      ! hydrometeors.
      wphydrometp(gr%nz,i) = zero

      ! Add implicit terms to the LHS matrix
      call microphys_lhs & 
           ( trim( hydromet_list(i) ), l_hydromet_sed(i), & ! In
             dt, Kr, nu_r_vert_res_dep, wm_zt, &  ! In
             hydromet_vel(:,i), hydromet_vel_zt(:,i), & ! In
             hydromet_vel_covar_zt_impc(:,i), & ! In
             rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, & ! In
             lhs ) ! Out

      ! Set up explicit term in the RHS vector
      if ( trim( hydromet_list(i) ) == "Ncm" .and. l_in_cloud_Nc_diff ) then

         call microphys_rhs( trim( hydromet_list(i) ), dt, l_hydromet_sed(i), &
                             Nc_in_cloud, &
                             hydromet_mc(:,i)/max(cloud_frac,cloud_frac_min), &
                             Kr, nu_r_vert_res_dep, cloud_frac, &
                             hydromet_vel_covar_zt_expc(:,i), &
                             rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                             rhs )

      else

         call microphys_rhs( trim( hydromet_list(i) ), dt, l_hydromet_sed(i), &
                             hydromet(:,i), hydromet_mc(:,i), &
                             Kr, nu_r_vert_res_dep, cloud_frac, &
                             hydromet_vel_covar_zt_expc(:,i), &
                             rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                             rhs )

      endif


      !!!!! Advance hydrometeor one time step.
      if ( trim( hydromet_list(i) ) == "Ncm" .and. l_in_cloud_Nc_diff ) then

         ! Ncm in cloud is computed above
         call microphys_solve( trim( hydromet_list(i) ), l_hydromet_sed(i), &
                               cloud_frac, &
                               lhs, rhs, Nc_in_cloud, err_code )

         hydromet(:,iiNcm) = Nc_in_cloud * max( cloud_frac, cloud_frac_min )

      else

         call microphys_solve( trim( hydromet_list(i) ), l_hydromet_sed(i), &
                               cloud_frac, &
                               lhs, rhs, hydromet(:,i), err_code )

      endif


      ! Print warning message if any hydrometeor species has a value < 0.
      if ( any( hydromet(:,i) < zero_threshold ) ) then

         hydromet_name = hydromet_list(i)

         if ( clubb_at_least_debug_level( 1 ) ) then
            do k = 1, gr%nz
               if ( hydromet(k,i) < zero_threshold ) then
                  write(fstderr,*) trim( hydromet_name ) //" < ", &
                                   zero_threshold, &
                                   " in advance_microphys at k= ", k
               endif
            enddo
         endif

      endif ! hydromet(:,i) < 0


      ! Store the previous value of the hydrometeor for the effect of the
      ! hole-filling scheme.
      if ( l_stats_samp ) then
         call stat_begin_update( ixrm_hf, hydromet(:,i) &
                                          / real( dt, kind = core_rknd ), zt )
      endif

      ! If we're dealing with a mixing ratio and hole filling is enabled,
      ! then we apply the hole filling algorithm
      if ( any( hydromet(:,i) < zero_threshold ) ) then

         if ( hydromet_name(1:1) == "r" .and. l_hole_fill ) then

            ! Apply the hole filling algorithm
            call fill_holes_driver( 2, zero_threshold, "zt", &
                                    rho_ds_zt, rho_ds_zm, &
                                    hydromet(:,i) )

         endif ! Variable is a mixing ratio and l_hole_fill is true

      endif ! hydromet(:,i) < 0

      ! Enter the new value of the hydrometeor for the effect of the
      ! hole-filling scheme.
      if ( l_stats_samp ) then
         call stat_end_update( ixrm_hf, hydromet(:,i) &
                                        / real( dt, kind = core_rknd ), zt )
      endif

      ! Store the previous value of the hydrometeor for the effect of the water
      ! vapor hole-filling scheme.
      if ( l_stats_samp ) then
         call stat_begin_update( ixrm_wvhf, hydromet(:,i) &
                                            / real( dt, kind = core_rknd ), zt )
      endif

      if ( any( hydromet(:,i) < zero_threshold ) ) then

         if ( hydromet_name(1:1) == "r" .and. l_hole_fill ) then

            ! If the hole filling algorithm failed, then we attempt to fill
            ! the missing mass with water vapor mixing ratio.
            ! We noticed this is needed for ASEX A209, particularly if Latin
            ! hypercube sampling is enabled.  -dschanen 11 Nov 2010
            do k = 2, gr%nz, 1

               if ( hydromet(k,i) < zero_threshold ) then

                  ! Set temp to the time tendency applied to vapor and removed
                  ! from the hydrometeor.
                  temp = hydromet(k,i) / real( dt, kind = core_rknd )

                  ! Adjust the tendency rvm_mc accordingly
                  rvm_mc(k) = rvm_mc(k) + temp

                  ! Adjust the tendency of thlm_mc according to whether the
                  ! effect is an evaporation or sublimation tendency.
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

               endif ! hydromet(k,i) < 0

            enddo ! k = 2..gr%nz

         endif ! Variable is a mixing ratio and l_hole_fill is true

      endif ! hydromet(:,i) < 0

      ! Enter the new value of the hydrometeor for the effect of the water vapor
      ! hole-filling scheme.
      if ( l_stats_samp ) then
         call stat_end_update( ixrm_wvhf, hydromet(:,i) &
                                          / real( dt, kind = core_rknd ), zt )
      endif

      ! Store the previous value of the hydrometeor for the effect of clipping.
      if ( l_stats_samp ) then
         call stat_begin_update( ixrm_cl, hydromet(:,i) & 
                                          / real( dt, kind = core_rknd ), zt )
      endif

      if ( any( hydromet(:,i) < zero_threshold ) ) then

           ! Clip any remaining negative values of hydrometeors to 0.
           ! This includes the case where the variable is a number
           ! concentration and is therefore not conserved.
           where ( hydromet(:,i) < zero_threshold )
              hydromet(:,i) = zero_threshold
           end where

      endif ! hydromet(:,i) < 0

      ! Enter the new value of the hydrometeor for the effect of clipping.
      if ( l_stats_samp ) then
         call stat_end_update( ixrm_cl, hydromet(:,i) &
                                        / real( dt, kind = core_rknd ), zt )
      endif


      ! Lower boundary condition
      ! Hydrometeors that are below the model lower boundary level have
      ! sedimented out of the model domain, and is not conserved.
      if ( hydromet(1,i) < zero_threshold ) then
         hydromet(1,i) = zero_threshold
      endif


      ! Solve for < w'h_m' > at all intermediate (momentum) grid levels, using
      ! a down-gradient approximation:  < w'h_m' > = - K * d< h_m >/dz.
      ! A Crank-Nicholson time-stepping scheme is used for this variable.
      ! This is the portion of the calculation using < h_m > from timestep t+1. 
      wphydrometp(1:gr%nz-1,i) &
      = - one_half * xpwp_fnc( Kr(1:gr%nz-1)+nu_r_vert_res_dep(1:gr%nz-1), &
                               hydromet(1:gr%nz-1,i), hydromet(2:gr%nz,i), &
                               gr%invrs_dzm(1:gr%nz-1) )

      ! A zero-flux boundary condition at the top of the model is used for
      ! hydrometeors.
      wphydrometp(gr%nz,i) = zero

      !!! Calculate the covariance of hydrometeor sedimentation velocity and the
      !!! hydrometeor, which is solved semi-implicitly on thermodynamic levels.
      hydromet_vel_covar_zt(:,i) &
      = hydromet_vel_covar_zt_impc(:,i) * hydromet(:,i) &
        + hydromet_vel_covar_zt_expc(:,i)

      ! Boundary conditions for < V_hm'hm' >|_zt.
      hydromet_vel_covar_zt(1,i) = hydromet_vel_covar_zt(2,i)

      !!! Calculate the covariance of hydrometeor sedimentation velocity and the
      !!! hydrometeor, < V_hm'h_m' >, by interpolating the thermodynamic level
      !!! results to momentum levels.
      hydromet_vel_covar(:,i) = zt2zm( hydromet_vel_covar_zt(:,i) )

      ! Boundary conditions for < V_hm'hm' >.
      hydromet_vel_covar(1,i)     = hydromet_vel_covar_zt(2,i)
      hydromet_vel_covar(gr%nz,i) = zero

      ! Statistics for all covariances involving hydrometeors:  < w'h_m' >,
      ! <V_rr'r_r'>, and <V_Nr'N_r'>.
      if ( l_stats_samp ) then

         if ( iwpxrp > 0 ) then

            ! Covariance of vertical velocity and the hydrometeor.
            call stat_update_var( iwpxrp, wphydrometp(:,i), zm )

         endif

         if ( trim( hydromet_list(i) ) == "rrainm" .and. iVrrprrp > 0 ) then

            ! Covariance of sedimentation velocity of r_r and r_r.
            call stat_update_var( iVrrprrp, hydromet_vel_covar(:,iirrainm), zm )

         elseif ( trim( hydromet_list(i) ) == "Nrm" .and. iVNrpNrp > 0 ) then

            ! Covariance of sedimentation velocity of N_r and N_r.
            call stat_update_var( iVNrpNrp, hydromet_vel_covar(:,iiNrm), zm )

         endif

      endif ! l_stats_samp


      if ( l_stats_samp ) then

         ! Total time tendency
         if ( l_hydromet_sed(i) .and. l_upwind_diff_sed ) then
            call stat_end_update( ixrm_bt, &
                                  hydromet(:,i) &
                                  / real( dt, kind = core_rknd ), zt )
         else
            ! The mean hydrometeor field at thermodynamic level k = 1 is simply
            ! set to the value of the mean hydrometeor field at k = 2.  Don't
            ! include budget stats for level k = 1.
            do k = 2, gr%nz, 1
               call stat_end_update_pt( ixrm_bt, k, &
                                        hydromet(k,i) &
                                        / real( dt, kind = core_rknd ), zt )
            enddo
         endif

      endif ! l_stats_samp


    enddo ! i=1..hydromet_dim


    if ( l_cloud_sed ) then

      ! For the GCSS DYCOMS II RF02 case the following lines of code are 
      ! currently setup to use a specified cloud droplet number
      ! concentration (Ncm).  The cloud droplet concentration has
      ! been moved here instead of being stated in KK_microphys
      ! for the following reasons:
      !    a) The effects of cloud droplet sedimentation can be computed
      !       without having to call the precipitation scheme.
      !    b) Ncm tends to be a case-specific parameter.  Therefore, it
      !       is appropriate to declare in the same place as other
      !       case-specific parameters.
      !
      ! Someday, we could move the setting of Ncm to pdf_closure
      ! for the following reasons:
      !    a) The cloud water mixing ratio (rcm) is computed using the
      !       PDF scheme.  Perhaps someday Ncm can also be computed by
      !       the same scheme.
      !    b) It seems more appropriate to declare Ncm in the same place
      !       where rcm is computed.
      !
      ! Since cloud base (z_cloud_base) is determined by the mixing ratio rc_tol,
      ! so will cloud droplet number concentration (Ncm).

      call cloud_drop_sed( rcm, Ncm, rho_zm, rho, & ! Intent(in)
                           exner, sigma_g,  &       ! Intent(in)
                           rcm_mc, thlm_mc )        ! Intent(inout)
    end if ! l_cloud_sed

    if ( l_stats_samp ) then
      call stat_update_var( iNcnm, Ncnm, zt )

      ! In the case where Ncm is neither a fixed number nor predicted these will
      ! be set to -999.
      call stat_update_var( iNcm, Ncm, zt )
      call stat_update_var( iNc_in_cloud, Nc_in_cloud, zt )
    end if ! l_stats_samp

    if ( l_stats_samp .and. iirrainm > 0 ) then

      ! Rainfall rate (mm/day) is defined on thermodynamic levels.  The rainfall
      ! rate is given by the equation:
      ! Rainfall rate (mm/day) = - V_rr * r_r * ( rho_a / rho_lw )
      !                            * ( 86400 s/day ) * ( 1000 mm/m ).
      ! The level average rainfall rate is given by:
      ! < Rainfall rate (mm/day) > = - < V_rr * r_r > * ( rho_a / rho_lw )
      !                                * ( 86400 s/day ) * ( 1000 mm/m );
      ! which can also be written as:
      ! < Rainfall rate (mm/day) > = - ( < V_rr > * < r_r > + < V_rr'r_r' > )
      !                                * ( rho_a / rho_lw )
      !                                * ( 86400 s/day ) * ( 1000 mm/m ).
      ! Rainfall rate is defined as positive.  Since V_rr is always negative,
      ! the minus (-) sign is necessary.
      call stat_update_var( irain_rate_zt,  & 
                            max( - ( hydromet(:,iirrainm) &
                                     * hydromet_vel_zt(:,iirrainm) &
                                     + hydromet_vel_covar_zt(:,iirrainm) ), &
                                 zero ) &
                            * ( rho / rho_lw ) & 
                            * real( sec_per_day, kind = core_rknd ) &
                            * mm_per_m, zt )

      ! Precipitation Flux (W/m^2) is defined on momentum levels.  The
      ! precipitation flux is given by the equation:
      ! Precip. flux (W/m^2) = - V_rr * r_r * rho_a * L_v.
      ! The level average precipitation flux is given by:
      ! < Precip. flux (W/m^2) > = - < V_rr * r_r > * rho_a * L_v;
      ! which can also be written as:
      ! < Precip. flux (W/m^2) > = - ( < V_rr > * < r_r > + < V_rr'r_r' > )
      !                              * rho_a * L_v.
      ! It is generally a convention in meteorology to show Precipitation Flux
      ! as a positive downward quantity, so the minus (-) sign is necessary.
      call stat_update_var( iFprec,  & 
                            max( - ( zt2zm( hydromet(:,iirrainm) )  & 
                                     * hydromet_vel(:,iirrainm) &
                                     + hydromet_vel_covar(:,iirrainm) ), &
                                 zero ) &
                            * rho_zm * Lv, zm )

      ! Store values of surface fluxes for statistics
      ! See notes above.
      call stat_update_var_pt( irain_rate_sfc, 1,  & 
                               max( - ( hydromet(2,iirrainm) &
                                        * hydromet_vel_zt(2,iirrainm) &
                                        + hydromet_vel_covar_zt(2,iirrainm) ), &
                                    zero ) &
                                * ( rho(2) / rho_lw ) & 
                                * real( sec_per_day, kind = core_rknd ) &
                                * mm_per_m, sfc )

      call stat_update_var_pt( irain_flux_sfc, 1, & 
                               max( - ( zt2zm( hydromet(:,iirrainm), 1 )  & 
                                        * hydromet_vel(1,iirrainm) &
                                        + hydromet_vel_covar(1,iirrainm) ), &
                                    zero ) &
                               * rho_zm(1) * Lv, sfc )

      ! Also store the value of surface rain water mixing ratio.
      call stat_update_var_pt( irrainm_sfc, 1,  & 
                               ( zt2zm( hydromet(:,iirrainm), 1 ) ), sfc )

    endif ! l_stats_samp

    call stats_accumulate_hydromet( hydromet, rho_ds_zt )

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

  !=============================================================================
  subroutine microphys_solve( solve_type, l_sed, &
                              cloud_frac, &
                              lhs, rhs, hmm, err_code )

    ! Description:
    ! Solve the tridiagonal system for hydrometeor variable.

    ! References:
    !  None
    !---------------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    use clubb_precision, only:  & 
        core_rknd    ! Variable(s)

    use lapack_wrap, only:  & 
        tridag_solve    ! Procedure(s)

    use error_code, only: &
        clubb_no_error ! Constant

    use stats_variables, only: & 
        zt,  & ! Variable(s)
        irrainm_ma, & 
        irrainm_sd, & 
        irrainm_ts, & 
        irrainm_ta, & 
        iricem_ma, & 
        iricem_sd, & 
        iricem_ta, & 
        irsnowm_ma, & 
        irsnowm_sd, & 
        irsnowm_ta, & 
        irgraupelm_ma, & 
        irgraupelm_sd, & 
        irgraupelm_ta, & 
        l_stats_samp, & 
        ztscr01, & 
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
        ztscr06, & 
        ztscr07, & 
        ztscr08, & 
        ztscr09, &
        ztscr10, &
        ztscr11, &
        ztscr12

    use stats_variables, only: & 
        iNrm_ma, & 
        iNrm_sd, & 
        iNrm_ts, & 
        iNrm_ta, & 
        iNim_ma, & 
        iNim_sd, & 
        iNim_ta, & 
        iNsnowm_ma, & 
        iNsnowm_sd, & 
        iNsnowm_ta, & 
        iNgraupelm_ma, & 
        iNgraupelm_sd, & 
        iNgraupelm_ta

    use stats_variables, only: & 
        iNcm_ma, & 
        iNcm_ta

    use stats_type, only: &
        stat_update_var_pt, & ! Procedure(s)
        stat_end_update_pt

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      solve_type  ! Description of which hydrometeor is being solved for.

    logical, intent(in) ::  & 
      l_sed    ! Whether to add a hydrometeor sedimentation term.

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      cloud_frac    ! Cloud fraction (thermodynamic levels)        [-]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(3,gr%nz) :: & 
      lhs    ! Left hand side

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: & 
      rhs    ! Right hand side vector

    real( kind = core_rknd ), intent(inout), dimension(gr%nz) :: & 
      hmm    ! Mean value of hydrometeor (thermodynamic levels)    [units vary]

    ! Output Variables
    integer, intent(out) :: err_code

    ! Local Variables
    integer :: k, km1, kp1  ! Array indices

    integer :: & 
      ihmm_ma,  & ! Mean advection budget stats toggle
      ihmm_ta,  & ! Turbulent advection budget stats toggle
      ihmm_sd,  & ! Mean sedimentation budget stats toggle
      ihmm_ts     ! Turbulent sedimentation budget stats toggle


    err_code = clubb_no_error  ! Initialize to the value for no errors

    ! Initializing ihmm_ma, ihmm_ta, ihmm_sd, ihmm_ts, and in order to avoid
    ! compiler warnings.
    ihmm_ma = 0
    ihmm_ta = 0
    ihmm_sd = 0
    ihmm_ts = 0

    select case( solve_type )
    case( "rrainm" )
      ihmm_ma = irrainm_ma
      ihmm_ta = irrainm_ta
      ihmm_sd = irrainm_sd
      ihmm_ts = irrainm_ts
    case( "ricem" )
      ihmm_ma = iricem_ma
      ihmm_ta = iricem_ta
      ihmm_sd = iricem_sd
      ihmm_ts = 0
    case( "rsnowm" )
      ihmm_ma = irsnowm_ma
      ihmm_ta = irsnowm_ta
      ihmm_sd = irsnowm_sd
      ihmm_ts = 0
    case( "rgraupelm" )
      ihmm_ma = irgraupelm_ma
      ihmm_ta = irgraupelm_ta
      ihmm_sd = irgraupelm_sd
      ihmm_ts = 0
    case( "Ncm" )
      ihmm_ma = iNcm_ma
      ihmm_ta = iNcm_ta
      ihmm_sd = 0
      ihmm_ts = 0
    case( "Nrm" )
      ihmm_ma = iNrm_ma
      ihmm_ta = iNrm_ta
      ihmm_sd = iNrm_sd
      ihmm_ts = iNrm_ts
    case( "Nim" )
      ihmm_ma = iNim_ma
      ihmm_ta = iNim_ta
      ihmm_sd = iNim_sd
      ihmm_ts = 0
    case( "Nsnowm" )
      ihmm_ma = iNsnowm_ma
      ihmm_ta = iNsnowm_ta
      ihmm_sd = iNsnowm_sd
      ihmm_ts = 0
    case( "Ngraupelm" )
      ihmm_ma = iNgraupelm_ma
      ihmm_ta = iNgraupelm_ta
      ihmm_sd = iNgraupelm_sd
      ihmm_ts = 0
    case default
      ihmm_ma = 0
      ihmm_ta = 0
      ihmm_sd = 0
      ihmm_ts = 0
    end select

    ! Solve system using tridag_solve. This uses LAPACK sgtsv,
    ! which relies on Gaussian elimination to decompose the matrix.
    call tridag_solve( solve_type, gr%nz, 1, lhs(1,:), lhs(2,:), lhs(3,:), & 
                       rhs, hmm, err_code )


    ! Statistics
    if ( l_stats_samp ) then

       do k = 1, gr%nz, 1

          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nz )

          ! Finalize implicit contributions

          ! hmm term ma is completely implicit; call stat_update_var_pt.
          if ( solve_type == "Ncm" .and. l_in_cloud_Nc_diff ) then

             ! For Ncm, we divide by cloud_frac when entering the subroutine,
             ! but do not multiply until we return from the subroutine, so we
             ! must account for this here for the budget to balance.
             call stat_update_var_pt( ihmm_ma, k, & 
               ztscr01(k) * hmm(km1) * max( cloud_frac(k), cloud_frac_min ) & 
               + ztscr02(k) * hmm(k) * max( cloud_frac(k), cloud_frac_min ) & 
               + ztscr03(k) * hmm(kp1) * max( cloud_frac(k), cloud_frac_min ), &
                                      zt )

          else

             call stat_update_var_pt( ihmm_ma, k, & 
                                      ztscr01(k) * hmm(km1) & 
                                      + ztscr02(k) * hmm(k) & 
                                      + ztscr03(k) * hmm(kp1), zt)

          endif

          ! hmm term sd is completely implicit; call stat_update_var_pt.
          if ( l_sed ) then
             call stat_update_var_pt( ihmm_sd, k, & 
                                      ztscr04(k) * hmm(km1) & 
                                      + ztscr05(k) * hmm(k) & 
                                      + ztscr06(k) * hmm(kp1), zt )
          endif

          ! hmm term ts has both implicit and explicit components; call
          ! stat_end_update_pt.
          if ( l_sed .and. k > 1 ) then
             call stat_end_update_pt( ihmm_ts, k, & 
                                      ztscr07(k) * hmm(km1) & 
                                      + ztscr08(k) * hmm(k) & 
                                      + ztscr09(k) * hmm(kp1), zt )
          endif

          ! hmm term ta has both implicit and explicit components; call
          ! stat_end_update_pt.
          if ( k > 1 ) then

             if ( solve_type == "Ncm" .and. l_in_cloud_Nc_diff ) then

                ! For Ncm, we divide by cloud_frac when entering the subroutine,
                ! but do not multiply until we return from the subroutine, so we
                ! must account for this here for the budget to balance.
                call stat_end_update_pt( ihmm_ta, k, & 
               ztscr10(k) * hmm(km1) * max( cloud_frac(k), cloud_frac_min ) & 
               + ztscr11(k) * hmm(k) * max( cloud_frac(k), cloud_frac_min ) & 
               + ztscr12(k) * hmm(kp1) * max( cloud_frac(k), cloud_frac_min ), &
                                         zt )

             else

                call stat_end_update_pt( ihmm_ta, k, & 
                                         ztscr10(k) * hmm(km1) & 
                                         + ztscr11(k) * hmm(k) & 
                                         + ztscr12(k) * hmm(kp1), zt )

             endif

          endif ! k > 1

       enddo ! 1..gr%nz

    endif ! l_stats_samp


    ! Boundary conditions on results
    !hmm(1) = hmm(2)
    ! Michael Falk, 7 Sep 2007, made this change to eliminate problems
    ! with anomalous rain formation at the top boundary.
    !        hmm(gr%nz) = 0
    !hmm(gr%nz) = hmm(gr%nz-1)
    ! eMFc

    return

  end subroutine microphys_solve

  !=============================================================================
  subroutine microphys_lhs & 
             ( solve_type, l_sed, dt, Kr, nu, wm_zt, &
               V_hm, V_hmt, &
               Vhmphmp_zt_impc, &
               rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
               lhs )

    ! Description:
    ! Setup the matrix of implicit contributions to a term.
    ! Can include the effects of sedimentation, diffusion, and advection.
    ! The Morrison microphysics has an explicit sedimentation code, which is
    ! handled elsewhere.
    !
    ! Notes:
    ! Setup for tridiagonal system and boundary conditions should be the same as
    ! the original rain subroutine code.
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr,    & ! Variable(s)
        zm2zt, & ! Procedure(s)
        zt2zm    ! Procedure(s)

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    use diffusion, only:  & 
        diffusion_zt_lhs ! Procedure(s)

    use mean_adv, only:  & 
        term_ma_zt_lhs ! Procedure(s)

    use constants_clubb, only: &
        one,         & ! Constant(s)
        one_half,    &
        zero,        & 
        sec_per_day

    use stats_variables, only: & 
        irrainm_ma, & ! Variable(s)
        irrainm_sd, & 
        irrainm_ts, & 
        irrainm_ta, & 
        iNrm_ma, & 
        iNrm_sd, & 
        iNrm_ts, & 
        iNrm_ta, & 
        iNcm_ma, &
        iNcm_ta, &
        iricem_ma, & 
        iricem_sd, & 
        iricem_ta, & 
        irsnowm_ma, & 
        irsnowm_sd, & 
        irsnowm_ta, & 
        irgraupelm_ma, & 
        irgraupelm_sd, & 
        irgraupelm_ta, & 
        iNim_ma, & 
        iNim_sd, & 
        iNim_ta, & 
        iNsnowm_ma, & 
        iNsnowm_sd, & 
        iNsnowm_ta, & 
        iNgraupelm_ma, & 
        iNgraupelm_sd, & 
        iNgraupelm_ta

    use stats_variables, only: & 
        ztscr01, & 
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
        ztscr06, & 
        ztscr07, & 
        ztscr08, & 
        ztscr09, & 
        ztscr10, & 
        ztscr11, & 
        ztscr12, & 
        l_stats_samp

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    ! Input Variables
    character(len=*), intent(in) :: &
      solve_type  ! Description of which hydrometeor is being solved for.

    logical, intent(in) ::  & 
      l_sed    ! Whether to add a hydrometeor sedimentation term.

    real(kind=time_precision), intent(in) ::  & 
      dt       ! Model timestep                                          [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      nu       ! Background diffusion coefficient                        [m^2/s]

    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  & 
      wm_zt, & ! w wind component on thermodynamic levels                [m/s]
      V_hm,  & ! Sedimentation velocity of hydrometeor (momentum levels) [m/s]
      V_hmt, & ! Sedimentation velocity of hydrometeor (thermo. levels)  [m/s]
      Kr       ! Eddy diffusivity for hydrometeor on momentum levels     [m^2/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: & 
      Vhmphmp_zt_impc, & ! Implicit comp. of <V_hm'h_m'> on t-levs  [units(m/s)]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    real( kind = core_rknd ), intent(out), dimension(3,gr%nz) :: & 
      lhs      ! Left hand side of tridiagonal matrix.

    ! Local Variables
    real( kind = core_rknd ), dimension(3) :: tmp

    real( kind = core_rknd ), dimension(gr%nz) :: &
      Vhmphmp_impc    ! Implicit comp. <V_hm'h_m'>: interp. m-levs  [units(m/s)]

    integer :: k, km1, kp1  ! Array indices
    integer :: diff_k_in

    integer :: & 
      ihmm_ma, & ! Mean advection budget stats toggle
      ihmm_ta, & ! Turbulent advection budget stats toggle
      ihmm_sd, & ! Mean sedimentation budget stats toggle
      ihmm_ts    ! Turbulent sedimentation budget stats toggle


    ! Initializing ihmm_ma, ihmm_sd, ihmm_ts, and ihmm_ta in order to avoid
    ! compiler warnings.
    ihmm_ma = 0
    ihmm_ta = 0
    ihmm_sd = 0
    ihmm_ts = 0

    select case( solve_type )
    case ( "rrainm" )
      ihmm_ma = irrainm_ma
      ihmm_ta = irrainm_ta
      ihmm_sd = irrainm_sd
      ihmm_ts = irrainm_ts
    case ( "Nrm" )
      ihmm_ma = iNrm_ma
      ihmm_ta = iNrm_ta
      ihmm_sd = iNrm_sd
      ihmm_ts = iNrm_ts
    case ( "ricem" )
      ihmm_ma = iricem_ma
      ihmm_ta = iricem_ta
      ihmm_sd = iricem_sd
      ihmm_ts = 0
    case ( "rsnowm" )
      ihmm_ma = irsnowm_ma
      ihmm_ta = irsnowm_ta
      ihmm_sd = irsnowm_sd
      ihmm_ts = 0
    case ( "rgraupelm" )
      ihmm_ma = irgraupelm_ma
      ihmm_ta = irgraupelm_ta
      ihmm_sd = irgraupelm_sd
      ihmm_ts = 0
    case ( "Ncm" )
      ihmm_ma = iNcm_ma
      ihmm_ta = iNcm_ta
      ihmm_sd = 0
      ihmm_ts = 0
    case( "Nim" )
      ihmm_ma = iNim_ma
      ihmm_ta = iNim_ta
      ihmm_sd = iNim_sd
      ihmm_ts = 0
    case( "Nsnowm" )
      ihmm_ma = iNsnowm_ma
      ihmm_ta = iNsnowm_ta
      ihmm_sd = iNsnowm_sd
      ihmm_ts = 0
    case( "Ngraupelm" )
      ihmm_ma = iNgraupelm_ma
      ihmm_ta = iNgraupelm_ta
      ihmm_sd = iNgraupelm_sd
      ihmm_ts = 0
    case default
      ihmm_ma = 0
      ihmm_ta = 0
      ihmm_sd = 0
      ihmm_ts = 0
    end select


    ! Interpolate the implicit component of < V_hm'h_m' >, a momentum-level
    ! variable that is calculated on thermodynamic levels, from thermodynamic
    ! levels to momentum levels.
    Vhmphmp_impc = zt2zm( Vhmphmp_zt_impc )

    ! Reset LHS Matrix for current timestep.
    lhs = zero

    ! Setup LHS Matrix
    do k = 2, gr%nz-1, 1

       km1 = max( k-1, 1 )
       kp1 = min( k+1, gr%nz )

       ! LHS time tendency.
       lhs(k_tdiag,k) = lhs(k_tdiag,k) + ( one / real( dt, kind = core_rknd ) )

       ! LHS turbulent advection term.
       ! - (1/rho_ds) * d( rho_ds * <w'h_m'> ) / dz.
       ! Note:  a down gradient closure approximation is made for < w'h_m' >, so
       !        the turbulent advection term is solved as an eddy-diffusion
       !        term:  + (1/rho_ds) * d( rho_ds * K_hm * (dh_m/dz) ) / dz.
       ! A Crank-Nicholson time-stepping scheme is used for this term.
       if ( k == 2 ) then
          ! The lower boundary condition needs to be applied here at level 2.
          ! The lower boundary condition is a zero-flux boundary condition.
          ! A hydrometeor is not allowed to be fluxed through the model lower
          ! boundary by the processes of mean or turbulent advection.  Only
          ! mean or turbulent sedimentation can flux a hydrometeor through the
          ! lower boundary.  Subroutine diffusion_zt_lhs is set-up to apply a
          ! zero-flux boundary condition at thermodynamic level 1.  In order to
          ! apply the same boundary condition code here at level 2, an adjuster
          ! needs to be used to tell diffusion_zt_lhs to use the code at level 2
          ! that it normally uses at level 1.
          diff_k_in = 1
       else
          diff_k_in = k
       endif
       lhs(kp1_tdiag:km1_tdiag,k) & 
       = lhs(kp1_tdiag:km1_tdiag,k) &
         + one_half &
           * invrs_rho_ds_zt(k) &
           * diffusion_zt_lhs( rho_ds_zm(k) * Kr(k), &
                               rho_ds_zm(km1) * Kr(km1), nu, & 
                               gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                               gr%invrs_dzt(k), diff_k_in )

       ! LHS mean advection term.
       lhs(kp1_tdiag:km1_tdiag,k) & 
       = lhs(kp1_tdiag:km1_tdiag,k) & 
         + term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), &
                           gr%invrs_dzm(km1) )

       if ( l_sed ) then

          ! LHS mean sedimentation term.
          ! Note: the Morrison microphysics has its own sedimentation code,
          ! which is applied through the _mc terms to each hydrometeor species.
          ! Therefore, l_sed will always be false when the Morrison microphysics
          ! is enabled.  -dschanen 24 Jan 2011
          if ( .not. l_upwind_diff_sed ) then

             ! Sedimentation (both mean and turbulent) uses centered
             ! differencing.  This is the default method.
             lhs(kp1_tdiag:km1_tdiag,k) & 
             = lhs(kp1_tdiag:km1_tdiag,k) & 
               + sed_centered_diff_lhs( V_hm(k), V_hm(km1), rho_ds_zm(k), &
                                        rho_ds_zm(km1), invrs_rho_ds_zt(k), &
                                        gr%invrs_dzt(k), k )

          else

             ! Sedimentation (both mean and turbulent) uses "upwind"
             ! differencing.
             lhs(kp1_tdiag:km1_tdiag,k) & 
             = lhs(kp1_tdiag:km1_tdiag,k) & 
               + sed_upwind_diff_lhs( V_hmt(k), V_hmt(kp1), rho_ds_zt(k), &
                                      rho_ds_zt(kp1), invrs_rho_ds_zt(k), &
                                      gr%invrs_dzm(k), k )

          endif

          ! LHS turbulent sedimentation term.
          lhs(kp1_tdiag:km1_tdiag,k) & 
          = lhs(kp1_tdiag:km1_tdiag,k) & 
             + term_turb_sed_lhs( Vhmphmp_impc(k), Vhmphmp_impc(km1), &
                                  Vhmphmp_zt_impc(kp1), Vhmphmp_zt_impc(k), &
                                  rho_ds_zm(k), rho_ds_zm(km1), &
                                  rho_ds_zt(kp1), rho_ds_zt(k), &
                                  gr%invrs_dzt(k), gr%invrs_dzm(k), &
                                  invrs_rho_ds_zt(k), k )

       endif ! l_sed

       if ( l_stats_samp ) then

          ! Statistics:  implicit contributions to hydrometeor hmm.

          if ( ihmm_ma > 0 ) then
             tmp(1:3) &
             = term_ma_zt_lhs( wm_zt(k), gr%invrs_dzt(k), k, gr%invrs_dzm(k), &
                               gr%invrs_dzm(km1) )

             ztscr01(k) = -tmp(3)
             ztscr02(k) = -tmp(2)
             ztscr03(k) = -tmp(1)
          endif

          if ( ihmm_sd > 0 .and. l_sed ) then
             if ( .not. l_upwind_diff_sed ) then
                tmp(1:3) &
                = sed_centered_diff_lhs( V_hm(k), V_hm(km1), rho_ds_zm(k), &
                                         rho_ds_zm(km1), invrs_rho_ds_zt(k), &
                                         gr%invrs_dzt(k), k )
             else
                tmp(1:3) &
                = sed_upwind_diff_lhs( V_hmt(k), V_hmt(kp1), rho_ds_zt(k), &
                                       rho_ds_zt(kp1), invrs_rho_ds_zt(k), &
                                       gr%invrs_dzm(k), k )
             endif

             ztscr04(k) = -tmp(3)
             ztscr05(k) = -tmp(2)
             ztscr06(k) = -tmp(1)

          endif

          if ( ihmm_ts > 0 .and. l_sed ) then
             tmp(1:3) &
             = term_turb_sed_lhs( Vhmphmp_impc(k), Vhmphmp_impc(km1), &
                                  Vhmphmp_zt_impc(kp1), Vhmphmp_zt_impc(k), &
                                  rho_ds_zm(k), rho_ds_zm(km1), &
                                  rho_ds_zt(kp1), rho_ds_zt(k), &
                                  gr%invrs_dzt(k), gr%invrs_dzm(k), &
                                  invrs_rho_ds_zt(k), k )
             ztscr07(k) = -tmp(3)
             ztscr08(k) = -tmp(2)
             ztscr09(k) = -tmp(1)
          endif

          if ( ihmm_ta > 0 ) then
             tmp(1:3) &
             = one_half &
               * invrs_rho_ds_zt(k) & 
               * diffusion_zt_lhs( rho_ds_zm(k) * Kr(k), &
                                   rho_ds_zm(km1) * Kr(km1), nu,  & 
                                   gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                                   gr%invrs_dzt(k), diff_k_in )
             ztscr10(k) = -tmp(3)
             ztscr11(k) = -tmp(2)
             ztscr12(k) = -tmp(1)
          endif

       endif ! l_stats_samp

    enddo ! 2..gr%nz-1


    ! Boundary Conditions

    ! Lower Boundary
    k   = 1
    km1 = max( k-1, 1 )
    kp1 = k+1

    if ( l_sed .and. l_upwind_diff_sed ) then

       ! LHS time tendency at the lower boundary.
       lhs(k_tdiag,k) = lhs(k_tdiag,k) + ( one / real( dt, kind = core_rknd ) )

       ! Here we apply the upwind differencing at the lower boundary.
       lhs(kp1_tdiag:km1_tdiag,k) & 
       = lhs(kp1_tdiag:km1_tdiag,k) & 
       + sed_upwind_diff_lhs( V_hmt(k), V_hmt(kp1), rho_ds_zt(k), &
                              rho_ds_zt(kp1), invrs_rho_ds_zt(k), &
                              gr%invrs_dzm(k), k )

    else

       ! This is set so that < h_m > at thermodynamic level k = 1, which is
       ! below the model lower boundary, is equal to < h_m > at k = 2.
       lhs(k_tdiag,k)   = one
       lhs(kp1_tdiag,k) = -one

    endif  ! l_sed and l_upwind_diff_sed

    if ( l_stats_samp ) then

       ! Statistics:  implicit contributions to hydrometeor hmm.

       if ( ihmm_sd > 0 .and. l_sed .and. l_upwind_diff_sed ) then
          tmp(1:3) &
          = sed_upwind_diff_lhs( V_hmt(k), V_hmt(kp1), rho_ds_zt(k), &
                                 rho_ds_zt(kp1), invrs_rho_ds_zt(k), &
                                 gr%invrs_dzm(k), k )

          ztscr04(k) = -tmp(3)
          ztscr05(k) = -tmp(2)
          ztscr06(k) = -tmp(1)
      endif

    endif  ! l_stats_samp


    ! Upper Boundary
    k   = gr%nz
    km1 = max( k-1, 1 )

    ! LHS time tendency at the upper boundary.
    lhs(k_tdiag,k) = lhs(k_tdiag,k) + ( one / real( dt, kind = core_rknd ) )

    ! LHS eddy-diffusion term at the upper boundary.
    lhs(kp1_tdiag:km1_tdiag,k) &
    = lhs(kp1_tdiag:km1_tdiag,k) &
      + one_half &
        * invrs_rho_ds_zt(k) &
        * diffusion_zt_lhs( rho_ds_zm(k) * Kr(k), &
                            rho_ds_zm(km1) * Kr(km1), nu, & 
                            gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                            gr%invrs_dzt(k), k )

    if ( l_stats_samp ) then

       ! Statistics:  implicit contributions to hydrometeor hmm.

       if ( ihmm_ta > 0 ) then
          tmp(1:3) & 
          = one_half &
            * invrs_rho_ds_zt(k) &
             * diffusion_zt_lhs( rho_ds_zm(k) * Kr(k), &
                                 rho_ds_zm(km1) * Kr(km1), nu, & 
                                 gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                                 gr%invrs_dzt(k), k )
          ztscr10(k) = -tmp(3)
          ztscr11(k) = -tmp(2)
          ztscr12(k) = -tmp(1)
       endif

    endif  ! l_stats_samp


    return

  end subroutine microphys_lhs

  !=============================================================================
  subroutine microphys_rhs( solve_type, dt, l_sed, &
                            hmm, hmm_tndcy, &
                            Kr, nu, cloud_frac, &
                            Vhmphmp_zt_expc, &
                            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt, &
                            rhs )

    ! Description:
    ! Compute RHS vector for a given hydrometeor.
    ! This subroutine computes the explicit portion of the predictive equation
    ! for a given hydrometeor.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr,    & ! Variable(s)
        zt2zm    ! Procedure(s)

    use constants_clubb, only: &
        one_half, & ! Constant(s)
        zero

    use diffusion, only:  & 
        diffusion_zt_lhs ! Procedure(s)

    use clubb_precision, only:  & 
        time_precision, & ! Variable(s)
        core_rknd

    use stats_variables, only: & 
        irrainm_ta, & ! Variable(s)
        irrainm_ts, &
        iNrm_ta, &
        iNrm_ts, &
        zt, &
        l_stats_samp

    use stats_variables, only: &
        iNim_ta, &
        iNsnowm_ta, &
        iNgraupelm_ta, &
        iNcm_ta, &
        iricem_ta, &
        irsnowm_ta, &
        irgraupelm_ta

    use stats_type, only: &
        stat_begin_update_pt ! Procedure(s)

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      solve_type  ! Description of which hydrometeor is being solved for.

    real( kind = time_precision ), intent(in) :: &
      dt    ! Duration of model timestep     [s]

    logical, intent(in) :: &
      l_sed    ! Flag for hydrometeor sedimentation

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      hmm,             & ! Mean value of hydrometeor (t-levs.)      [units]
      hmm_tndcy,       & ! Microphysics tendency (thermo. levels)   [units/s]
      Kr,              & ! Eddy diffusivity for hydromet on m-levs  [m^2/s]
      nu,              & ! Background diffusion coefficient         [m^2/s]
      cloud_frac,      & ! Cloud fraction                           [-]
      Vhmphmp_zt_expc, & ! Explicit comp. of <V_hm'h_m'> on t-levs  [units(m/s)]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zt    ! Inv. dry, static density @ thermo. levs. [m^3/kg]

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: & 
      rhs   ! Right hand side

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      Vhmphmp_expc    ! Explicit comp. <V_hm'h_m'>: interp. m-levs  [units(m/s)]

    real( kind = core_rknd ), dimension(3) :: &
      rhs_diff    ! For use in Crank-Nicholson eddy diffusion

    integer :: k, kp1, km1  ! Array indices
    integer :: diff_k_in

    integer :: ihmm_ta, & ! Turbulent advection budget toggle.
               ihmm_ts    ! Turbulent sedimentation budget toggle.
 

    ! Initializing ihmm_ta and ihmm_ts in order to avoid compiler warnings.
    ihmm_ta = 0
    ihmm_ts = 0

    select case( solve_type )
    case ( "rrainm" )
      ihmm_ta = irrainm_ta
      ihmm_ts = irrainm_ts
    case ( "Nrm" )
      ihmm_ta = iNrm_ta
      ihmm_ts = iNrm_ts
    case( "ricem" )
      ihmm_ta = iricem_ta
      ihmm_ts = 0
    case( "rsnowm" )
      ihmm_ta = irsnowm_ta
      ihmm_ts = 0
    case( "rgraupelm" )
      ihmm_ta = irgraupelm_ta
      ihmm_ts = 0
    case( "Ncm" )
      ihmm_ta = iNcm_ta
      ihmm_ts = 0
    case( "Nim" )
      ihmm_ta = iNim_ta
      ihmm_ts = 0
    case( "Nsnowm" )
      ihmm_ta = iNsnowm_ta
      ihmm_ts = 0
    case( "Ngraupelm" )
      ihmm_ta = iNgraupelm_ta
      ihmm_ts = 0
    case default
      ihmm_ta = 0
      ihmm_ts = 0
    end select


    ! Interpolate the explicit component of < V_hm'h_m' >, a momentum-level
    ! variable that is calculated on thermodynamic levels, from thermodynamic
    ! levels to momentum levels.
    Vhmphmp_expc = zt2zm( Vhmphmp_zt_expc )

    ! Initialize right-hand side vector to 0.
    rhs = zero

    ! Hydrometeor right-hand side (explicit portion of the code).
    do k = 2, gr%nz, 1

       km1 = max( k-1, 1 )
       kp1 = min( k+1, gr%nz )

       ! RHS time tendency.
       rhs(k) = hmm(k) / real( dt, kind = core_rknd )

       ! RHS microphysics tendency term (autoconversion, accretion, evaporation,
       ! etc.).
       rhs(k) = rhs(k) + hmm_tndcy(k)

       ! RHS turbulent advection term.
       ! - (1/rho_ds) * d( rho_ds * <w'h_m'> ) / dz.
       ! Note:  a down gradient closure approximation is made for < w'h_m' >, so
       !        the turbulent advection term is solved as an eddy-diffusion
       !        term:  + (1/rho_ds) * d( rho_ds * K_hm * (dh_m/dz) ) / dz.
       ! A Crank-Nicholson time-stepping scheme is used for this term.
       if ( k == 2 ) then
          ! The lower boundary condition needs to be applied here at level 2.
          ! The lower boundary condition is a zero-flux boundary condition.
          ! A hydrometeor is not allowed to be fluxed through the model lower
          ! boundary by the processes of mean or turbulent advection.  Only
          ! mean or turbulent sedimentation can flux a hydrometeor through the
          ! lower boundary.  Subroutine diffusion_zt_lhs is set-up to apply a
          ! zero-flux boundary condition at thermodynamic level 1.  In order to
          ! apply the same boundary condition code here at level 2, an adjuster
          ! needs to be used to tell diffusion_zt_lhs to use the code at level 2
          ! that it normally uses at level 1.
          diff_k_in = 1
       else
          diff_k_in = k
       endif
       rhs_diff(1:3) &
       = one_half &
         * invrs_rho_ds_zt(k) &
         * diffusion_zt_lhs( rho_ds_zm(k) * Kr(k), &
                             rho_ds_zm(km1) * Kr(km1), nu, & 
                             gr%invrs_dzm(km1), gr%invrs_dzm(k), &
                             gr%invrs_dzt(k), diff_k_in )

       rhs(k) &
       = rhs(k) &
         - rhs_diff(3) * hmm(km1) &
         - rhs_diff(2) * hmm(k) &
         - rhs_diff(1) * hmm(kp1)

       ! RHS turbulent sedimentation term.
       if ( l_sed ) then
          rhs(k) &
          = rhs(k) &
            + term_turb_sed_rhs( Vhmphmp_expc(k), Vhmphmp_expc(km1), &
                                 Vhmphmp_zt_expc(kp1), Vhmphmp_zt_expc(k), &
                                 rho_ds_zm(k), rho_ds_zm(km1), &
                                 rho_ds_zt(kp1), rho_ds_zt(k), &
                                 gr%invrs_dzt(k), gr%invrs_dzm(k), &
                                 invrs_rho_ds_zt(k), k )
       endif

       
       if ( l_stats_samp ) then

          ! Statistics: explicit contributions for the hydrometeor.

          ! hmm term ta has both implicit and explicit components; call
          ! stat_begin_update_pt.  Since stat_begin_update_pt automatically
          ! subtracts the value sent in, reverse the sign on the right-hand side
          ! turbulent advection component.
          if ( ihmm_ta > 0 ) then

             if ( solve_type == "Ncm" .and. l_in_cloud_Nc_diff ) then

                ! For Ncm, we divide by cloud_frac when entering the subroutine,
                ! but do not multiply until we return from the subroutine, so we
                ! must account for this here for the budget to balance.
                call stat_begin_update_pt( ihmm_ta, k, & 
              rhs_diff(3) * hmm(km1) * max( cloud_frac(k), cloud_frac_min ) & 
              + rhs_diff(2) * hmm(k) * max( cloud_frac(k), cloud_frac_min ) & 
              + rhs_diff(1) * hmm(kp1) * max( cloud_frac(k), cloud_frac_min ), &
                                           zt )

             else

                call stat_begin_update_pt( ihmm_ta, k, & 
                                           rhs_diff(3) * hmm(km1) &
                                           + rhs_diff(2) * hmm(k)   &
                                           + rhs_diff(1) * hmm(kp1), zt )

             endif

          endif

          ! hmm term ts has both implicit and explicit components; call
          ! stat_update_var_pt.  Since stat_begin_update_pt automatically
          ! subtracts the value sent in, reverse the sign on term_turb_sed_rhs.
          if ( ihmm_ts > 0 .and. l_sed ) then
             call stat_begin_update_pt( ihmm_ts, k, &
                 -term_turb_sed_rhs( Vhmphmp_expc(k), Vhmphmp_expc(km1), &
                                     Vhmphmp_zt_expc(kp1), Vhmphmp_zt_expc(k), &
                                     rho_ds_zm(k), rho_ds_zm(km1), &
                                     rho_ds_zt(kp1), rho_ds_zt(k), &
                                     gr%invrs_dzt(k), gr%invrs_dzm(k), &
                                     invrs_rho_ds_zt(k), k ), &
                                        zt )
          endif ! ihmm_ts > 0 and l_sed

       endif ! l_stats_samp

    enddo ! k = 2, gr%nz, 1
    

    ! Lower boundary conditions on the RHS

    if ( l_sed .and. l_upwind_diff_sed ) then

       ! RHS time tendency at the lower boundary.
       rhs(1) = hmm(1) / real( dt, kind = core_rknd )

    else

       ! This is set so that < h_m > at thermodynamic level k = 1, which is
       ! below the model lower boundary, is equal to < h_m > at k = 2.
       rhs(1) = zero

    endif  ! l_sed and l_upwind_diff_sed


    return

  end subroutine microphys_rhs

  !=============================================================================
  pure function sed_centered_diff_lhs( V_hm, V_hmm1, rho_ds_zm, &
                                       rho_ds_zmm1, invrs_rho_ds_zt, &
                                       invrs_dzt, level ) &
    result( lhs )

    ! Description:
    ! Mean sedimentation of a hydrometeor:  implicit portion of the code, using
    ! the centered difference approximation to the vertical derivative.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The mean sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm >|_(t) * < hm >|_(t+1) ) / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of the term is changed to
    !        a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d<hm>/dt equation.
    !
    ! This term is discretized as follows when using the centered-difference
    ! approximation:
    !
    ! The values of <hm> are found on the thermodynamic levels, while the values
    ! of <V_hm> are found on the momentum levels.  Additionally, the values of
    ! rho_ds_zm are found on the momentum levels, and the values of
    ! invrs_rho_ds_zt are found on the thermodynamic levels.  The variable <hm>
    ! is interpolated to the intermediate momentum levels.  At the intermediate
    ! momentum levels, the interpolated values of <hm> are multiplied by the
    ! values of <V_hm> and the values of rho_ds_zm.  Then, the derivative of
    ! (rho_ds*<V_hm>*<hm>) is taken over the central thermodynamic level, where
    ! it is multiplied by invrs_rho_ds_zt.
    !
    ! -----hmp1------------------------------------------------ t(k+1)
    !
    ! =============hm(interp)=====V_hm=====rho_ds_zm=========== m(k)
    !
    ! -----hm--------invrs_rho_ds_zt----d(rho_ds*V_hm*hm)/dz--- t(k)
    !
    ! =============hm(interp)=====V_hmm1===rho_ds_zmm1========= m(k-1)
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
    ! When a hydrometeor is sedimented to the ground (or out the lower boundary
    ! of the model), it is removed from the atmosphere (or from the model
    ! domain).  Thus, the quantity of the hydrometeor over the entire vertical
    ! domain should not be conserved due to the process of sedimentation.  Thus,
    ! not all of the column totals in the left-hand side matrix should be equal
    ! to 0. Instead, the sum of all the column totals should equal the flux of
    ! <hm> out the bottom (zm(1) level) of the domain,
    ! -rho_ds_zm(1) * V_hm(1) * ( D(2)*hm(1) + C(2)*hm(2) ), where the factor in
    ! parentheses is the interpolated value of hm at the zm(1) level.
    ! Furthermore, most of the individual column totals should sum to 0, but the
    ! 1st and 2nd (from the left) columns should combine to sum to the flux out
    ! the bottom of the domain.
    !
    ! To see that this modified conservation law is satisfied, compute the
    ! sedimentation of hm and integrate vertically.  In discretized matrix
    ! notation (where "i" stands for the matrix column and "j" stands for the
    ! matrix row):
    !
    ! - rho_ds_zm(1) * V_hm(1) * ( D(2)*hm(1) + C(2)*hm(2) )
    ! = Sum_j Sum_i
    !   ( 1 / invrs_rho_ds_zt )_i * ( 1 / invrs_dzt )_i
    !   * ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij * hm_j.
    !
    ! The left-hand side matrix,
    ! ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij, is
    ! partially written below.  The sum over i in the above equation removes
    ! invrs_rho_ds_zt and invrs_dzt everywhere from the matrix below.  The sum
    ! over j leaves the column totals and the flux at zm(1) that are desired.
    !
    ! Left-hand side matrix contributions from the sedimentation term (only);
    ! first four vertical levels:
    !
    !     -------------------------------------------------------------------->
    !k=1 |           0                     0                       0
    !    |
    !k=2 |   -invrs_rho_ds_zt(k)  +invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    !    |    *invrs_dzt(k)        *invrs_dzt(k)          *invrs_dzt(k)
    !    |    *rho_ds_zm(k-1)      *[ rho_ds_zm(k)        *rho_ds_zm(k)
    !    |    *V_hm(k-1)*D(k)         *V_hm(k)*B(k)       *V_hm(k)*A(k)
    !    |                           -rho_ds_zm(k-1)
    !    |                            *V_hm(k-1)*C(k) ]
    !    |
    !k=3 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    !    |                         *invrs_dzt(k)          *invrs_dzt(k)
    !    |                         *rho_ds_zm(k-1)        *[ rho_ds_zm(k)
    !    |                         *V_hm(k-1)*D(k)           *V_hm(k)*B(k)
    !    |                                                  -rho_ds_zm(k-1)
    !    |                                                   *V_hm(k-1)*C(k) ]
    !    |
    !k=4 |           0                     0             -invrs_rho_ds_zt(k)
    !    |                                                *invrs_dzt(k)
    !    |                                                *rho_ds_zm(k-1)
    !    |                                                *V_hm(k-1)*D(k)
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
    ! k /= gr%nz and k /= 1), the four weighting factors have the following
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

    use constants_clubb, only: &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    integer, parameter :: & 
      t_above = 1, & ! Index for upper thermodynamic level grid weight.
      t_below = 2    ! Index for lower thermodynamic level grid weight.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      V_hm,            & ! Sedimentation velocity of hydrometeor (k)    [m/s]
      V_hmm1,          & ! Sedimentation velocity of hydrometeor (k-1)  [m/s]
      rho_ds_zm,       & ! Dry, static density at momentum level (k)    [kg/m^3]
      rho_ds_zmm1,     & ! Dry, static density at momentum level (k-1)  [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. level (k) [m^3/kg]
      invrs_dzt          ! Inverse of grid spacing (k)                  [m]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      mk,   & ! Momentum level directly above central thermodynamic level.
      mkm1    ! Momentum level directly below central thermodynamic level.

    ! ---- Begin Code ----

    ! Momentum level (k) is between thermodynamic level (k+1)
    ! and thermodynamic level (k).
    mk   = level
    ! Momentum level (k-1) is between thermodynamic level (k)
    ! and thermodynamic level (k-1).
    mkm1 = level - 1

    if ( level == 1 ) then

       ! k = 1 (bottom level); lower boundary level; no effects.

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag) = zero

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)   = zero

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag) = zero


    elseif ( level > 1 .and. level < gr%nz ) then

       ! Most of the interior model; normal conditions.

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag)  & 
       = + invrs_rho_ds_zt * invrs_dzt &
           * rho_ds_zm * V_hm * gr%weights_zt2zm(t_above,mk)

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)  & 
       = + invrs_rho_ds_zt &
           * invrs_dzt &
           * ( rho_ds_zm * V_hm * gr%weights_zt2zm(t_below,mk) & 
               - rho_ds_zmm1 * V_hmm1 * gr%weights_zt2zm(t_above,mkm1)  )

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag)  & 
       = - invrs_rho_ds_zt * invrs_dzt &
           * rho_ds_zmm1 * V_hmm1 * gr%weights_zt2zm(t_below,mkm1)


    elseif ( level == gr%nz ) then

       ! k = gr%nz (top level); upper boundary level; no flux.

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag) = zero

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)   = zero

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag) = zero


    endif


    return

  end function sed_centered_diff_lhs

  !=============================================================================
  pure function sed_upwind_diff_lhs( V_hmt, V_hmtp1, rho_ds_zt, &
                                     rho_ds_ztp1, invrs_rho_ds_zt, &
                                     invrs_dzm, level ) &
    result( lhs )

    ! Description:
    ! Mean sedimentation of a hydrometeor:  implicit portion of the code, using
    ! the "upwind" difference approximation to the vertical derivative.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The mean sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm >|_(t) * < hm >|_(t+1) ) / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of the term is changed to
    !        a "+".
    !
    ! Timestep index (t) stands for the index of the current timestep, while
    ! timestep index (t+1) stands for the index of the next timestep, which is
    ! being advanced to in solving the d<hm>/dt equation.
    !
    ! This term is discretized as follows when using the upwind-difference
    ! approximation:
    !
    ! The values of <hm> and the values of V_hmt are found on the thermodynamic
    ! levels.  Additionally, the values of rho_ds_zt and the values of
    ! invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    ! thermodynamic levels, the values of <hm> are multiplied by the values of
    ! V_hmt and the values of rho_ds_zt.  Then, the derivative of
    ! (rho_ds*<V_hm>*<hm>) is taken between the thermodynamic level above the
    ! central thermodynamic level and the central thermodynamic level.  The
    ! derivative is multiplied by invrs_rho_ds_zt.
    !
    ! --hmp1--V_hmtp1--rho_ds_ztp1--------------------------------------- t(k+1)
    !
    ! =================================================================== m(k)
    !
    ! --hm----V_hmt----rho_ds_zt--invrs_rho_ds_zt--d(rho_ds*V_hm*hm)/dz-- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )
    !
    !
    ! Conservation Properties:
    !
    ! When a hydrometeor is sedimented to the ground (or out the lower boundary
    ! of the model), it is removed from the atmosphere (or from the model
    ! domain).  Thus, the quantity of the hydrometeor over the entire vertical
    ! domain should not be conserved due to the process of sedimentation.  Thus,
    ! not all of the column totals in the left-hand side matrix should be equal
    ! to 0. Instead, the sum of all the column totals should equal the flux of
    ! <hm> out the bottom (zm(1) level, approximated by the zt(1) level for the
    ! "upwind" sedimentation option) of the domain,
    ! -rho_ds_zt(1) * V_hmt(1) * hm(1).  Furthermore, most of the individual
    ! column totals should sum to 0, but the 2nd (from the left) column should
    ! be equal to the flux out the bottom of the domain.
    !
    ! To see that this modified conservation law is satisfied, compute the
    ! sedimentation of hm and integrate vertically.  In discretized matrix
    ! notation (where "i" stands for the matrix column and "j" stands for the
    ! matrix row):
    !
    ! - rho_ds_zt(1) * V_hmt(1) * hm(1)
    ! = Sum_j Sum_i
    !   ( 1 / invrs_rho_ds_zt )_i * ( 1 / invrs_dzm )_i
    !   * ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij * hm_j.
    !
    ! The left-hand side matrix,
    ! ( invrs_rho_ds_zt * d(rho_ds_zm * V_hm * weights_hm) / dz )_ij, is
    ! partially written below.  The sum over i in the above equation removes
    ! invrs_rho_ds_zt and invrs_dzm everywhere from the matrix below.  The sum
    ! over j leaves the column totals and the flux at zt(1) that are desired.
    !
    ! Left-hand side matrix contributions from the sedimentation term (only);
    ! first three vertical levels:
    !
    !     -------------------------------------------------------------------->
    !k=1 | -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)              0
    !    |  *invrs_dzm(k)          *invrs_dzm(k)
    !    |  *rho_ds_zt(k)          *rho_ds_zt(k+1)
    !    |  *V_hmt(k)              *V_hmt(k+1)
    !    |
    !k=2 |           0            -invrs_rho_ds_zt(k)    +invrs_rho_ds_zt(k)
    !    |                         *invrs_dzm(k)          *invrs_dzm(k)
    !    |                         *rho_ds_zt(k)          *rho_ds_zt(k+1)
    !    |                         *V_hmt(k)              *V_hmt(k+1)
    !    |
    !k=3 |           0                     0             -invrs_rho_ds_zt(k)
    !    |                                                *invrs_dzm(k)
    !    |                                                *rho_ds_zt(k)
    !    |                                                *V_hmt(k)
    !   \ /
    !
    ! Note:  The superdiagonal term from level 3 is not shown on this diagram.

    ! References:
    ! None

    ! Notes:
    ! Both COAMPS Microphysics and Brian Griffin's implementation use
    ! Khairoutdinov and Kogan (2000) for the calculation of rain
    ! mixing ratio and rain droplet number concentration sedimentation
    ! velocities, but COAMPS has only the local parameterization.
    !
    ! Please note that "upwind" sedimentation is only 1st-order accurate and
    ! highly diffusive. 
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable(s)

    use constants_clubb, only: &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      V_hmtp1,         & ! Sed. velocity of hydrometeor at t-lev (k+1)  [m/s]
      V_hmt,           & ! Sed. velocity of hydrometeor at t-lev (k)    [m/s]
      rho_ds_ztp1,     & ! Dry, static density at thermo. level (k+1)   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density at thermo. level (k)     [kg/m^3]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. level (k) [m^3/kg]
      invrs_dzm          ! Inverse of grid spacing over m-lev. (k)      [m]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! ---- Begin Code ----

    ! Sedimention is always a downward process, so we omit the upward case
    ! (i.e. the V_hmt variable will always be negative).
    if ( level == gr%nz ) then

       ! k = gr%nz (top level); upper boundary level; no flux.

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag) = zero

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)   = zero

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag) = zero


    else

       ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
       lhs(kp1_tdiag) = + invrs_rho_ds_zt * invrs_dzm * rho_ds_ztp1 * V_hmtp1

       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
       lhs(k_tdiag)   = - invrs_rho_ds_zt * invrs_dzm * rho_ds_zt * V_hmt

       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
       lhs(km1_tdiag) = zero


    endif


    return

  end function sed_upwind_diff_lhs

  !=============================================================================
  pure function term_turb_sed_lhs( Vhmphmp_impcm, Vhmphmp_impcm1, &
                                   Vhmphmp_zt_impcp1, Vhmphmp_zt_impc, &
                                   rho_ds_zm, rho_ds_zmm1, &
                                   rho_ds_ztp1, rho_ds_zt, &
                                   invrs_dzt, invrs_dzm, &
                                   invrs_rho_ds_zt, level ) &
    result( lhs )

    ! Description:
    ! Turbulent sedimentation of a hydrometeor:  implicit portion of the code.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The turbulent sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! This term is solved for semi-implicitly by rewritting < V_hm'hm' > based
    ! on < hm > in the manner:
    !
    ! < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
    !
    ! This term can also be solved for completely explicitly (it's original
    ! form) by setting Vhmphmp_inc to 0 and setting Vhmphmp_expc to
    ! < V_hm'hm' >.  The equation becomes:
    !
    ! - (1/rho_ds) 
    !   * d( rho_ds * ( Vhmphmp_impc * < hm >(t+1) + Vhmphmp_expc ) ) / dz;
    !
    ! where the timestep index (t+1) means that the value of < hm > being used
    ! is from the next timestep, which is being advanced to in solving the
    ! d<hm>/dt equation.  Implicit and explicit portions of this term are
    ! produced.  The implicit portion of this term is:
    !
    ! - (1/rho_ds) * d( rho_ds * Vhmphmp_impc * < hm >(t+1) ) / dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign is
    !        reversed and the leading "-" in front of the d[ ] / dz term is
    !        changed to a "+".
    !
    ! This term can be discretized using the centered-difference approximation
    ! (which is preferred), or else using the "upwind"-difference approximation.
    !
    ! The implicit portion of this term is discretized as follows when using
    ! the centered-difference approximation:
    ! 
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_impc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zm are found on
    ! the momentum levels, and the values of invrs_rho_ds_zt are found on the
    ! thermodynamic levels.  The variables < hm > and Vhmphmp_zt_impc are both
    ! interpolated to the intermediate momentum levels.  At the momentum levels,
    ! the values of interpolated < hm > and interpolated Vhmphmp_zt_impc are
    ! multiplied together, and their products are multiplied by the values of
    ! rho_ds_zm.  The mathematical expression F is the product of these three
    ! variables at momentum levels.  Then, the derivative dF/dz is taken over
    ! the central thermodynamic level, where it is multiplied by
    ! invrs_rho_ds_zt.  In this function, the value of F is as follows:
    !
    ! F = rho_ds_zm * Vhmphmp_impc(interp) * hmm(interp).
    !
    !
    ! ----hmmp1--------Vhmphmp_zt_impcp1--------------------------------- t(k+1)
    !
    ! =====hmm(interp)=====Vhmphmp_impc(interp)=======rho_ds_zm========== m(k)
    !
    ! ----hmm----------Vhmphmp_zt_impc-----invrs_rho_ds_zt-----dF/dz----- t(k)
    !
    ! =====hmm(interp)=====Vhmphmp_impcm1(interp)=====rho_ds_zmm1======== m(k-1)
    !
    ! ----hmmm1--------Vhmphmp_zt_impcm1--------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is
    ! used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) ).
    !
    ! The implicit portion of this term is discretized as follows when using
    ! the upwind-difference approximation:
    !
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_impc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zt and the
    ! values of invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    ! thermodynamic levels, the values of < hm > and Vhmphmp_zt_impc are
    ! multiplied together, and their products are multiplied by the values of
    ! rho_ds_zt.  The mathematical expression F is the product of these three
    ! variables at thermodynamic levels.  Then, the derivative dF/dz is taken
    ! between the thermodynamic level above the central thermodynamic level and
    ! the central thermodynamic level.  The derivative is multiplied
    ! by invrs_rho_ds_zt.  In this function, the value of F is as follows:
    !
    ! F = rho_ds_zt * Vhmphmp_zt_impc * hmm.
    !
    !
    ! --hmmp1---Vhmphmp_zt_impcp1---rho_ds_ztp1-------------------------- t(k+1)
    !
    ! =================================================================== m(k)
    !
    ! --hmm-----Vhmphmp_zt_impc-----rho_ds_zt---invrs_rho_ds_zt---dF/dz-- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) ).

    ! References:
    !  None
    !
    ! Notes:
    ! Please note that "upwind" sedimentation is only 1st-order accurate and
    ! highly diffusive. 
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable(s)

    use constants_clubb, only: &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    integer, parameter :: & 
      t_above = 1, & ! Index for upper thermodynamic level grid weight.
      t_below = 2    ! Index for lower thermodynamic level grid weight.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Vhmphmp_impcm,     & ! Imp. comp. <V_hm'h_m'> interp. m-lev (k)   [vary]
      Vhmphmp_impcm1,    & ! Imp. comp. <V_hm'h_m'> interp. m-lev (k-1) [vary]
      Vhmphmp_zt_impcp1, & ! Imp. comp. <V_hm'h_m'>|_zt; t-lev (k+1)    [vary]
      Vhmphmp_zt_impc,   & ! Imp. comp. <V_hm'h_m'>|_zt; t-lev (k)      [vary]
      rho_ds_zm,         & ! Dry, static density at moment. lev (k)     [kg/m^3]
      rho_ds_zmm1,       & ! Dry, static density at moment. lev (k-1)   [kg/m^3]
      rho_ds_ztp1,       & ! Dry, static density at thermo. level (k+1) [kg/m^3]
      rho_ds_zt,         & ! Dry, static density at thermo. level (k)   [kg/m^3]
      invrs_dzt,         & ! Inverse of grid spacing over t-levs. (k)   [1/m]
      invrs_dzm,         & ! Inverse of grid spacing over m-levs. (k)   [1/m]
      invrs_rho_ds_zt      ! Inv dry, static density @ thermo lev (k)   [m^3/kg]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ), dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      mk,   & ! Momentum level directly above central thermodynamic level.
      mkm1    ! Momentum level directly below central thermodynamic level.


    ! Momentum level (k) is between thermodynamic level (k+1)
    ! and thermodynamic level (k).
    mk   = level
    ! Momentum level (k-1) is between thermodynamic level (k)
    ! and thermodynamic level (k-1).
    mkm1 = level - 1


    ! LHS (implicit component) of turbulent sedimentation term,
    ! - (1/rho_ds) * d( rho_ds * < V_hm'h_m' > ) / dz.
    ! = - (1/rho_ds)
    !     * d( rho_ds * ( Vhmphmp_impc * < h_m > + Vhmphmp_expc ) ) / dz.
    ! Implicit component:
    ! - (1/rho_ds) * d( rho_ds * Vhmphmp_impc * < h_m > ) / dz.
    if ( .not. l_upwind_diff_sed ) then

       ! Sedimentation (both mean and turbulent) uses centered differencing.
       if ( level == 1 ) then

          ! k = 1 (bottom level); lower boundary level; no effects.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) = zero

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)   = zero

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       elseif ( level > 1 .and. level < gr%nz ) then

          ! Most of the interior model; normal conditions.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag)  & 
          = + invrs_rho_ds_zt &
              * invrs_dzt &
              * rho_ds_zm * Vhmphmp_impcm * gr%weights_zt2zm(t_above,mk)

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)  & 
          = + invrs_rho_ds_zt &
              * invrs_dzt &
              * ( rho_ds_zm * Vhmphmp_impcm &
                            * gr%weights_zt2zm(t_below,mk) & 
                  - rho_ds_zmm1 * Vhmphmp_impcm1 &
                                * gr%weights_zt2zm(t_above,mkm1)  )

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag)  & 
          = - invrs_rho_ds_zt &
              * invrs_dzt &
              * rho_ds_zmm1 * Vhmphmp_impcm1 * gr%weights_zt2zm(t_below,mkm1)


       elseif ( level == gr%nz ) then

          ! k = gr%nz (top level); upper boundary level; no flux.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) = zero

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)   = zero

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       endif  ! level


    else ! l_upwind_diff_sed

       ! Sedimentation (both mean and turbulent) uses "upwind" differencing.
       if ( level == 1 ) then

          ! k = 1 (bottom level); lower boundary level; no effects.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) = zero

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)   = zero

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       elseif ( level > 1 .and. level < gr%nz ) then

          ! Most of the interior model; normal conditions.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) &
          = + invrs_rho_ds_zt &
              * invrs_dzm * rho_ds_ztp1 * Vhmphmp_zt_impcp1

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag) &
          = - invrs_rho_ds_zt &
              * invrs_dzm * rho_ds_zt * Vhmphmp_zt_impc

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       elseif ( level == gr%nz ) then

          ! k = gr%nz (top level); upper boundary level; no flux.

          ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
          lhs(kp1_tdiag) = zero

          ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
          lhs(k_tdiag)   = zero

          ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
          lhs(km1_tdiag) = zero


       endif  ! level


    endif


    return

  end function term_turb_sed_lhs

  !=============================================================================
  pure function term_turb_sed_rhs( Vhmphmp_expcm, Vhmphmp_expcm1, &
                                   Vhmphmp_zt_expcp1, Vhmphmp_zt_expc, &
                                   rho_ds_zm, rho_ds_zmm1, &
                                   rho_ds_ztp1, rho_ds_zt, &
                                   invrs_dzt, invrs_dzm, &
                                   invrs_rho_ds_zt, level ) &
    result( rhs )

    ! Description:
    ! Turbulent sedimentation of a hydrometeor:  explicit portion of the code.
    !
    ! The variable "hm" stands for a hydrometeor variable.  The variable "V_hm"
    ! stands for the sedimentation velocity of the aforementioned hydrometeor.
    !
    ! The d(hm)/dt equation contains a sedimentation term:
    !
    ! - (1/rho_ds) * d( rho_ds * V_hm * hm ) / dz.
    !
    ! The variables hm and V_hm in the sedimentation term are divided into mean
    ! and turbulent components, and the term is averaged, resulting in:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm > * < hm > ) / dz
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! The turbulent sedimentation term in the d<hm>/dt equation is:
    !
    ! - (1/rho_ds) * d( rho_ds * < V_hm'hm' > ) / dz.
    !
    ! This term is solved for semi-implicitly by rewritting < V_hm'hm' > based
    ! on < hm > in the manner:
    !
    ! < V_hm'hm' > = Vhmphmp_impc * < hm > + Vhmphmp_expc.
    !
    ! This term can also be solved for completely explicitly (it's original
    ! form) by setting Vhmphmp_inc to 0 and setting Vhmphmp_expc to
    ! < V_hm'hm' >.  The equation becomes:
    !
    ! - (1/rho_ds) 
    !   * d( rho_ds * ( Vhmphmp_impc * < hm >(t+1) + Vhmphmp_expc ) ) / dz;
    !
    ! where the timestep index (t+1) means that the value of < hm > being used
    ! is from the next timestep, which is being advanced to in solving the
    ! d<hm>/dt equation.  Implicit and explicit portions of this term are
    ! produced.  The explicit portion of this term is:
    !
    ! - (1/rho_ds) * d( rho_ds * Vhmphmp_expc ) / dz.
    !
    ! This term can be discretized using the centered-difference approximation
    ! (which is preferred), or else using the "upwind"-difference approximation.
    !
    ! The explicit portion of this term is discretized as follows when using
    ! the centered-difference approximation:
    !
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_expc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zm are found on
    ! the momentum levels, and the values of invrs_rho_ds_zt are found on the
    ! thermodynamic levels.  The variable Vhmphmp_zt_expc is interpolated to the
    ! intermediate momentum levels.  At the momentum levels, the values of
    ! interpolated Vhmphmp_zt_expc are multiplied by the values of rho_ds_zm.
    ! Then, the derivative d(rho_ds*Vhmphmp_zt_expc)/dz is taken over the
    ! central thermodynamic level, where it is multiplied by invrs_rho_ds_zt.
    !
    ! ---Vhmphmp_zt_expcp1----------------------------------------------- t(k+1)
    !
    ! ======Vhmphmp_expc(interp)=========rho_ds_zm======================= m(k)
    !
    ! ---Vhmphmp_zt_expc--invrs_rho_ds_zt--d(rho_ds*Vhmphmp_zt_expc)/dz-- t(k)
    !
    ! ======Vhmphmp_expcm1(interp)=======rho_ds_zmm1===================== m(k-1)
    !
    ! ---Vhmphmp_zt_expcm1----------------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is
    ! used for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) ).
    !
    ! The explicit portion of this term is discretized as follows when using
    ! the upwind-difference approximation:
    !
    ! The values of < hm > and the values of <V_hm'hm'>|_zt are found on the
    ! thermodynamic levels.  The values of Vhmphmp_zt_expc are also found on the
    ! thermodynamic levels.  Additionally, the values of rho_ds_zt and the
    ! values of invrs_rho_ds_zt are found on the thermodynamic levels.  At the
    ! thermodynamic levels, the values of Vhmphmp_zt_expc are multiplied by the
    ! values of rho_ds_zt.  The mathematical expression F is the product of
    ! these variables at thermodynamic levels.  Then, the derivative dF/dz is
    ! taken between the thermodynamic level above the central thermodynamic
    ! level and the central thermodynamic level.  The derivative is multiplied
    ! by invrs_rho_ds_zt.  In this function, the value of F is as follows:
    !
    ! F = rho_ds_zt * Vhmphmp_zt_expc.
    !
    !
    ! -----Vhmphmp_zt_expcp1---rho_ds_ztp1------------------------------- t(k+1)
    !
    ! =================================================================== m(k)
    !
    ! -----Vhmphmp_zt_expc-----rho_ds_zt----invrs_rho_ds_zt----dF/dz----- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) ).
    
    ! References:
    !  None
    !
    ! Notes:
    ! Please note that "upwind" sedimentation is only 1st-order accurate and
    ! highly diffusive. 
    !-----------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable(s)

    use constants_clubb, only: &
        zero  ! Constant(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Vhmphmp_expcm,     & ! Exp. comp. <V_hm'h_m'> interp. m-lev (k)   [vary]
      Vhmphmp_expcm1,    & ! Exp. comp. <V_hm'h_m'> interp. m-lev (k-1) [vary]
      Vhmphmp_zt_expcp1, & ! Exp. comp. <V_hm'h_m'>|_zt; t-lev (k+1)    [vary]
      Vhmphmp_zt_expc,   & ! Exp. comp. <V_hm'h_m'>|_zt; t-lev (k)      [vary]
      rho_ds_zm,         & ! Dry, static density at moment. lev (k)     [kg/m^3]
      rho_ds_zmm1,       & ! Dry, static density at moment. lev (k-1)   [kg/m^3]
      rho_ds_ztp1,       & ! Dry, static density at thermo. level (k+1) [kg/m^3]
      rho_ds_zt,         & ! Dry, static density at thermo. level (k)   [kg/m^3]
      invrs_dzt,         & ! Inverse of grid spacing over t-levs. (k)   [1/m]
      invrs_dzm,         & ! Inverse of grid spacing over m-levs. (k)   [1/m]
      invrs_rho_ds_zt      ! Inv dry, static density @ thermo lev (k)   [m^3/kg]

    integer, intent(in) ::  & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real( kind = core_rknd ) :: rhs


    ! RHS (explicit component) of turbulent sedimentation term,
    ! - (1/rho_ds) * d( rho_ds * < V_hm'h_m' > ) / dz.
    ! = - (1/rho_ds)
    !     * d( rho_ds * ( Vhmphmp_impc * < h_m > + Vhmphmp_expc ) ) / dz.
    ! Explicit component:  - (1/rho_ds) * d( rho_ds * Vhmphmp_expc ) / dz.
    if ( .not. l_upwind_diff_sed ) then

       ! Sedimentation (both mean and turbulent) uses centered differencing.
       if ( level == 1 ) then

          ! k = 1 (bottom level); lower boundary level; no effects.
          rhs = zero


       elseif ( level > 1 .and. level < gr%nz ) then

          ! Most of the interior model; normal conditions.
          rhs &
          = - invrs_rho_ds_zt &
              * invrs_dzt * ( rho_ds_zm * Vhmphmp_expcm &
                              - rho_ds_zmm1 * Vhmphmp_expcm1 )


       elseif ( level == gr%nz ) then

          ! k = gr%nz (top level); upper boundary level; no flux.
          rhs = zero


       endif


    else ! l_upwind_diff_sed

       ! Sedimentation (both mean and turbulent) uses "upwind" differencing.
       if ( level == 1 ) then

          ! k = 1 (bottom level); lower boundary level; no effects.
          rhs = zero


       elseif ( level > 1 .and. level < gr%nz ) then

          ! Most of the interior model; normal conditions.
          rhs &
          = - invrs_rho_ds_zt &
              * invrs_dzm * ( rho_ds_ztp1 * Vhmphmp_zt_expcp1 &
                              - rho_ds_zt * Vhmphmp_zt_expc )


       elseif ( level == gr%nz ) then

          ! k = gr%nz (top level); upper boundary level; no flux.
          rhs = zero


       endif


    endif


    return

  end function term_turb_sed_rhs

  !=============================================================================
  subroutine cleanup_microphys( )

    ! Description:
    ! De-allocate arrays used by the microphysics

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( hydromet_list ) ) then
      deallocate( hydromet_list )
    end if

    if ( allocated( l_hydromet_sed ) ) then
      deallocate( l_hydromet_sed )
    end if

    if ( trim( micro_scheme ) == "morrison_gettelman" ) then
      call pbuf_deallocate()
    end if

    return

  end subroutine cleanup_microphys
!===============================================================================

end module microphys_driver
