! $Id$
!===============================================================================
module microphys_init_cleanup

  ! Description:
  ! Initialize and cleanup code pertaining to microphysics schemes.

  ! References:
  !  H. Morrison, J. A. Curry, and V. I. Khvorostyanov, 2005: A new double-
  !  moment microphysics scheme for application in cloud and
  !  climate models. Part 1: Description. J. Atmos. Sci., 62, 1665-1677.

  !  Khairoutdinov, M. and Kogan, Y.: A new cloud physics parameterization in a
  !  large-eddy simulation model of marine stratocumulus, Mon. Wea. Rev., 128,
  !  229-243, 2000.
  !-------------------------------------------------------------------------

  implicit none

  public :: init_microphys,    &
            cleanup_microphys

  private ! Default Scope

  contains

  !=============================================================================
  subroutine init_microphys( iunit, runtype, namelist_file, case_info_file, &
                             hydromet_dim )

    ! Description:
    ! Set indices to the various hydrometeor species and define hydromet_dim for
    ! the purposes of allocating memory.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use parameters_microphys, only: &
        l_in_cloud_Nc_diff,           & ! Use in cloud values of Nc for diffusion
        l_cloud_sed,                  & ! Cloud water sedimentation (K&K or no microphysics)
        l_ice_microphys,              & ! Compute ice (COAMPS / Morrison)
        l_graupel,                    & ! Compute graupel (Morrison)
        l_hail,                       & ! See module_mp_graupel for a description
        l_seifert_beheng,             & ! Use Seifert and Beheng (2001) warm drizzle (Morrison)
        l_predict_Nc,                 & ! Predict cloud droplet number conc
        specify_aerosol,              & ! Specify aerosol (Morrison)
        l_subgrid_w,                  & ! Use subgrid w  (Morrison)
        l_arctic_nucl,                & ! Use MPACE observations (Morrison)
        l_cloud_edge_activation,      & ! Activate on cloud edges (Morrison)
        l_fix_pgam,                   & ! Fix pgam (Morrison)
        lh_sequence_length,           & ! Number of timesteps before the SILHS seq. repeats
        lh_seed,                      & ! Integer seed for the Mersenne twister
        l_local_kk,                   & ! Use local formula for K&K
        l_upwind_diff_sed,            & ! Use the upwind differencing approx. for sedimentation
        l_var_covar_src,              & ! Flag for using variance and covariance src terms
        microphys_scheme,             & ! The microphysical scheme in use
        l_hydromet_sed,               & ! Flag to sediment a hydrometeor
        l_gfdl_activation,            & ! Flag to use GFDL activation scheme
        microphys_start_time,         & ! When to start the microphysics [s]
        sigma_g,                      & ! Parameter used in the cloud droplet sedimentation code
        Nc0_in_cloud                    ! Initial value for Nc (K&K, l_cloud_sed, Morrison)

    use parameters_silhs, only: &
        l_lh_importance_sampling, &     ! Do importance sampling (SILHS)
        l_lh_straight_mc,         &     ! Do not apply LH or importance sampling at all (SILHS)
        l_lh_clustered_sampling,  &     ! Use prescribed probability sampling with clusters (SILHS)
        eight_cluster_presc_probs,&     ! Sampling probabilities for prescribed mode (SILHS)
        l_rcm_in_cloud_k_lh_start,&     ! Determine k_lh_start based on maximum within-cloud rcm
        l_random_k_lh_start,      &     ! k_lh_start found randomly between max rcm and rcm_in_cloud
        importance_prob_thresh,   &     ! Minimum PDF probability for importance sampling
        l_lh_limit_weights,       &     ! Ensure weights stay under a given value
        cluster_allocation_strategy, &  ! Strategy for distributing sample points
        l_lh_var_frac                   ! Prescribe variance fractions

    use parameters_microphys, only: &
        lh_num_samples                  ! SILHS sample points

    use parameters_microphys, only: &
        lh_microphys_interactive,     & ! Feed the subcolumns into microphysics and allow feedback
        lh_microphys_non_interactive, & ! Feed the subcolumns into microphysics with no feedback
        lh_microphys_disabled           ! Disable latin hypercube entirely

    use parameters_microphys, only: &
        l_silhs_KK_convergence_adj_mean ! Clip source adjustment terms on mean instead of individual
                                        ! sample points to test convergence with KK analytic
    use parameters_microphys, only: &
        morrison_no_aerosol, &  ! Constants
        morrison_power_law,  &
        morrison_lognormal

    use parameters_KK, only: &
        C_evap,              &
        r_0

    use parameters_microphys, only: &
        lh_microphys_type_int => lh_microphys_type ! Determines how the LH samples are used

    use array_index, only: &
        hydromet_list, & ! Names of the hydrometeor species
        hydromet_tol     ! Variable(s)

    use phys_buffer, only: & ! Used for placing wp2_zt in morrison_gettelman microphysics
        pbuf_init

    ! Adding coefficient variable for clex9_oct14 case to reduce NNUCCD
    ! and NNUCCC
    use module_mp_graupel, only: &
        NNUCCD_REDUCE_COEF, &
        NNUCCC_REDUCE_COEF
    ! Change by Marc Pilon on 11/16/11

    use array_index, only: & 
        l_frozen_hm, & ! Variables
        l_mix_rat_hm,&
        iirrm, &
        iiNrm, &
        iirsm, &
        iirim, &
        iirgm, &
        iiNsm, & 
        iiNim, &
        iiNgm

    use constants_clubb, only: &
        cm3_per_m3, &
        rr_tol,     &
        ri_tol,     &
        rs_tol,     &
        rg_tol,     &
        Nr_tol,     &
        Ni_tol,     &
        Ns_tol,     &
        Ng_tol,     &
        zero

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
        ini_micro                ! Subroutine

    use microp_aero, only: &
        ini_microp_aero ! Subroutine

    use constants_clubb, only: &
        fstdout, & ! Constant(s)
        fstderr

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
        core_rknd

    use precipitation_fraction, only: &
        precip_frac_calc_type  ! Variable(s)

    use corr_varnce_module, only: &
        hmp2_ip_on_hmm2_ip_ratios_type, & ! Type(s)
        hmp2_ip_on_hmm2_ip,             & ! Variable(s)
        Ncnp2_on_Ncnm2,                 &
        d_variables,                    &
        iiPDF_Ncn,                      &
        corr_array_n_cloud,             &
        corr_array_n_below,             &
        setup_pdf_indices,              & ! Procedure(s)
        setup_corr_varnce_array

    use pdf_utilities, only: &
        stdev_L2N    ! Procedure(s)

    use setup_clubb_pdf_params, only: &
        denorm_transform_corr    ! Procedure(s)

    use parameters_tunable, only: &
        omicron,        & ! Procedure(s)
        zeta_vrnce_rat

    use index_mapping, only: &
        pdf2hydromet_idx    ! Procedure(s)

    use matrix_operations, only: &
        mirror_lower_triangular_matrix    ! Procedure(s)

    use model_flags, only: &
        l_diagnose_correlations, &
        l_evaporate_cold_rcm, &
        l_morr_xp2_mc, &
        l_const_Nc_in_cloud, &  ! Use a constant cloud droplet conc. within cloud (K&K)
        l_fix_chi_eta_correlations  ! Use a fixed correlation for chi/eta(s/t Mellor) (SILHS)

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
    character(len=30) :: lh_microphys_type
    integer, parameter :: res = 20   ! Used for lookup tables with GFDL activation
    integer, parameter :: res2 = 20  ! Used for lookup tables with GFDL activation
    real( kind = core_rknd ), dimension( :, :, :, :, : ), allocatable :: &
      droplets, droplets2            ! Used for lookup tables with GFDL activation

    character(len=128) :: &
     corr_file_path_cloud, &
     corr_file_path_below

    type(hmp2_ip_on_hmm2_ip_ratios_type) :: hmp2_ip_on_hmm2_ip_ratios

    real( kind = core_rknd ), dimension(:), allocatable :: &
      sigma_x_n_cloud, & ! Std. devs. (normal space): PDF vars (in cloud) [u.v.]
      sigma_x_n_below    ! Std. devs. (normal space): PDF vars (below cl) [u.v.]

    real( kind = core_rknd ), dimension(:), allocatable :: &
      sigma2_on_mu2_ip_cloud, & ! Ratio array sigma_x_cloud^2/mu_x_cloud^2   [-]
      sigma2_on_mu2_ip_below    ! Ratio array sigma_x_below^2/mu_x_below^2   [-]

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      corr_array_cloud, & ! Correlation array of PDF vars. (in cloud)        [-]
      corr_array_below    ! Correlation array of PDF vars. (below cloud)     [-]

    integer :: ivar  ! Loop index


    namelist /microphysics_setting/ &
      microphys_scheme, l_cloud_sed, sigma_g, &
      l_ice_microphys, l_graupel, l_hail, l_var_covar_src, l_upwind_diff_sed, &
      l_seifert_beheng, l_predict_Nc, l_const_Nc_in_cloud, specify_aerosol, &
      l_subgrid_w, l_arctic_nucl, l_cloud_edge_activation, l_fix_pgam, &
      l_in_cloud_Nc_diff, lh_microphys_type, l_local_kk, lh_num_samples, &
      lh_sequence_length, lh_seed, l_lh_importance_sampling, &
      l_fix_chi_eta_correlations, l_silhs_KK_convergence_adj_mean, l_lh_straight_mc, &
      l_lh_clustered_sampling, eight_cluster_presc_probs, &
      l_rcm_in_cloud_k_lh_start, l_random_k_lh_start, importance_prob_thresh, &
      l_lh_limit_weights, cluster_allocation_strategy, l_lh_var_frac, &
      hmp2_ip_on_hmm2_ip_ratios, Ncnp2_on_Ncnm2, &
      C_evap, r_0, microphys_start_time, &
      Nc0_in_cloud, ccnconst, ccnexpnt, aer_rm1, aer_rm2, &
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

    lh_microphys_type = "disabled"

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

       call write_text ( "microphys_scheme = " //  microphys_scheme, l_write_to_file, &
                         iunit )
       call write_text ( "l_cloud_sed = ", l_cloud_sed, l_write_to_file, iunit )
       call write_text ( "sigma_g = ", sigma_g, l_write_to_file, iunit )
       call write_text ( "l_graupel = ", l_graupel, l_write_to_file, iunit )
       call write_text ( "l_hail = ", l_hail, l_write_to_file, iunit )
       call write_text ( "l_seifert_beheng = ", l_seifert_beheng, &
                         l_write_to_file, iunit )
       call write_text ( "l_predict_Nc = ", l_predict_Nc, l_write_to_file, iunit )
       call write_text ( "l_const_Nc_in_cloud = ", l_const_Nc_in_cloud, &
                         l_write_to_file, iunit )
       call write_text ( "specify_aerosol = "// specify_aerosol, &
                         l_write_to_file, iunit )
       call write_text ( "l_subgrid_w = ", l_subgrid_w, l_write_to_file, iunit )
       call write_text ( "l_arctic_nucl = ", l_arctic_nucl, l_write_to_file, &
                          iunit )
       call write_text ( "l_cloud_edge_activation = ", &
                         l_cloud_edge_activation, l_write_to_file, iunit )
       call write_text ( "l_fix_pgam = ", l_fix_pgam, l_write_to_file, iunit )
       call write_text ( "l_in_cloud_Nc_diff = ", l_in_cloud_Nc_diff, &
                         l_write_to_file, iunit )
       call write_text ( "l_var_covar_src = ", l_var_covar_src, &
                         l_write_to_file, iunit )
       call write_text ( "l_upwind_diff_sed = ", l_upwind_diff_sed, &
                         l_write_to_file, iunit )
       call write_text ( "lh_microphys_type = " // &
                         trim( lh_microphys_type ), l_write_to_file, iunit )
       call write_text ( "lh_num_samples = ", lh_num_samples, &
                         l_write_to_file, iunit )
       call write_text ( "lh_sequence_length = ", lh_sequence_length, &
                         l_write_to_file, iunit )
       call write_text ( "lh_seed = ", lh_seed, l_write_to_file, iunit )
       call write_text ( "l_lh_importance_sampling = ", &
                         l_lh_importance_sampling, l_write_to_file, iunit )
       call write_text ( "l_fix_chi_eta_correlations = ", l_fix_chi_eta_correlations, &
                         l_write_to_file, iunit )
       call write_text ( "l_silhs_KK_convergence_adj_mean = ", &
                         l_silhs_KK_convergence_adj_mean, &
                         l_write_to_file, iunit )
       call write_text ( "l_lh_straight_mc = ", l_lh_straight_mc, l_write_to_file, iunit )
       call write_text ( "l_lh_clustered_sampling = ", l_lh_clustered_sampling, l_write_to_file, &
                         iunit )
       call write_text ( "l_rcm_in_cloud_k_lh_start = ", l_rcm_in_cloud_k_lh_start, &
                         l_write_to_file, iunit )
       call write_text ( "l_random_k_lh_start = ", l_random_k_lh_start, l_write_to_file, iunit )
       call write_text ( "importance_prob_thresh = ", importance_prob_thresh, l_write_to_file, &
                         iunit )
       call write_text ( "l_lh_limit_weights = ", l_lh_limit_weights, l_write_to_file, iunit )
       call write_text ( "cluster_allocation_strategy = ", cluster_allocation_strategy, &
                         l_write_to_file, iunit )
       call write_text ( "l_lh_var_frac = ", l_lh_var_frac, l_write_to_file, iunit )
       call write_text ( "rrp2_ip_on_rrm2_ip = ", &
                         hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip, &
                         l_write_to_file, iunit )
       call write_text ( "Nrp2_ip_on_Nrm2_ip = ", &
                         hmp2_ip_on_hmm2_ip_ratios%Nrp2_ip_on_Nrm2_ip, &
                         l_write_to_file, iunit )
       call write_text ( "rip2_ip_on_rim2_ip = ", &
                         hmp2_ip_on_hmm2_ip_ratios%rip2_ip_on_rim2_ip, &
                         l_write_to_file, iunit )
       call write_text ( "Nip2_ip_on_Nim2_ip = ", &
                         hmp2_ip_on_hmm2_ip_ratios%Nip2_ip_on_Nim2_ip, &
                         l_write_to_file, iunit )
       call write_text ( "rsp2_ip_on_rsm2_ip = ", &
                         hmp2_ip_on_hmm2_ip_ratios%rsp2_ip_on_rsm2_ip, &
                         l_write_to_file, iunit )
       call write_text ( "Nsp2_ip_on_Nsm2_ip = ", &
                         hmp2_ip_on_hmm2_ip_ratios%Nsp2_ip_on_Nsm2_ip, &
                         l_write_to_file, iunit )
       call write_text ( "rgp2_ip_on_rgm2_ip = ", &
                         hmp2_ip_on_hmm2_ip_ratios%rgp2_ip_on_rgm2_ip, &
                         l_write_to_file, iunit )
       call write_text ( "Ngp2_ip_on_Ngm2_ip = ", &
                         hmp2_ip_on_hmm2_ip_ratios%Ngp2_ip_on_Ngm2_ip, &
                         l_write_to_file, iunit )
       call write_text ( "Ncnp2_on_Ncnm2 = ", Ncnp2_on_Ncnm2, &
                         l_write_to_file, iunit )
       call write_text ( "C_evap = ", C_evap, l_write_to_file, iunit )
       call write_text ( "r_0 = ", r_0, l_write_to_file, iunit )
       call write_text ( "microphys_start_time = ", &
                         real( microphys_start_time, kind = core_rknd ), &
                         l_write_to_file, iunit )
       call write_text ( "Nc0_in_cloud = ", Nc0_in_cloud, &
                         l_write_to_file, iunit )
       call write_text ( "ccnconst = ", real(ccnconst, kind = core_rknd), &
                         l_write_to_file, iunit )
       call write_text ( "ccnexpnt = ", real(ccnexpnt, kind = core_rknd), &
                         l_write_to_file, iunit )
       call write_text ( "aer_rm1 = ", real(aer_rm1, kind = core_rknd), &
                         l_write_to_file, iunit )
       call write_text ( "aer_rm2 = ", real(aer_rm2, kind = core_rknd), &
                         l_write_to_file, iunit )
       call write_text ( "aer_n1 = ", real(aer_n1, kind = core_rknd), &
                         l_write_to_file, iunit )
       call write_text ( "aer_n2 = ", real(aer_n2, kind = core_rknd), &
                         l_write_to_file, iunit )
       call write_text ( "aer_sig1 = ", real(aer_sig1, kind = core_rknd), &
                         l_write_to_file, iunit )
       call write_text ( "aer_sig2 = ", real(aer_sig2, kind = core_rknd), &
                         l_write_to_file, iunit )
       call write_text ( "pgam_fixed = ", real(pgam_fixed, kind = core_rknd), &
                         l_write_to_file, iunit )
       call write_text ( "precip_frac_calc_type = ", precip_frac_calc_type, &
                         l_write_to_file, iunit )

       if ( l_write_to_file ) close(unit=iunit)

    endif ! clubb_at_least_debug_level(1)

    ! Read in the name list for initialization, if it exists
    open(unit=iunit, file=namelist_file, status='old', action='read')
    read(iunit, nml=gfdl_activation_setting)
    close(unit=iunit)

    ! Initialize the GFDL activation code, if necessary
    if( l_gfdl_activation ) then
       ! Ensure a microphysics that has Ncm is being used
       if( trim( microphys_scheme ) == "coamps" .or. &
           trim( microphys_scheme ) == "morrison" & 
           .or. trim( microphys_scheme ) == "morrison_gettelman") then

          ! Read in the lookup tables
          call Loading( droplets, droplets2 )
          allocate( droplets(res,res,res,res,res), droplets2(res,res,res,res,res) )
          ! Initialize the activation variables
          call aer_ccn_act_k_init                            &
               ( real( droplets ), real( droplets2 ), res, res2, nooc,     &
                 real( sul_concen ), real( low_concen ), real( high_concen ),      &
                 real( lowup ), real( highup ), real( lowup2 ), real( highup2 ), real( lowmass2 ), &
                 real( highmass2 ), real( lowmass3 ), real( highmass3 ),           &
                 real( lowmass4 ), real( highmass4 ), real( lowmass5 ), real( highmass5 ), &
                 real( lowT2 ),real( highT2 ) )
          deallocate( droplets, droplets2 )
       endif ! coamps .or. morrison .or. khairoutdinov_kogan
    endif ! l_gfdl_activation

    ! The location of the fields in the hydromet array are arbitrary,
    ! and don't need to be set consistently among schemes so long as
    ! the 'i' indices point to the correct parts of the array.

    select case ( trim( microphys_scheme ) )

    case ( "morrison" )

       iirrm = 1
       iiNrm = 2

       if ( l_ice_microphys ) then

          iirim = 3
          iiNim = 4
          iirsm = 5
          iiNsm = 6

          doicemicro = .true.

          if ( l_graupel ) then

             iirgm = 7
             iiNgm = 8

             hydromet_dim = 8

             dograupel = .true.

          else ! l_graupel disabled

             iirgm = -1
             iiNgm = -1

             hydromet_dim = 6

             dograupel = .false.

          endif ! l_graupel

       else ! l_ice_microphys disabled

          iirsm = -1
          iirim = -1
          iiNsm = -1
          iiNim = -1
          iirgm = -1
          iiNgm = -1

          hydromet_dim = 2

          doicemicro = .false.
          dograupel = .false.

       endif ! l_ice_microphys

       ! Set Nc0 in the Morrison code (module_MP_graupel) based on Nc0_in_cloud
       Nc0 = real( Nc0_in_cloud / cm3_per_m3 ) ! Units on Nc0 are per cm^3

       ! Set flags from the Morrison scheme as in GRAUPEL_INIT
       if ( l_predict_Nc ) then
          dopredictNc = .true.
       else
          dopredictNc = .false.
       endif

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
       endif

       if ( l_arctic_nucl ) then
          doarcticicenucl = .true.
       else
          doarcticicenucl = .false.
       endif

       if ( l_hail ) then
          dohail = .true.
       else
          dohail = .false.
       endif

       if ( l_seifert_beheng ) then
          dosb_warm_rain = .true.
       else
          dosb_warm_rain = .false.
       endif

       if ( l_fix_pgam ) then
          dofix_pgam = .true.
       else
          dofix_pgam = .false.
       endif

       if ( l_subgrid_w ) then
          dosubgridw = .true.
       else
          dosubgridw = .false.
       endif

       if ( l_cloud_sed ) then
          write(fstderr,*) "Morrison microphysics has seperate code for" &
                           //" cloud water sedimentation, therefore " &
                           //"l_cloud_sed should be set to .false."
          stop "Fatal error."
       endif

       if ( .not. l_fix_chi_eta_correlations .and. l_ice_microphys &
            .and. trim( lh_microphys_type ) /= "disabled" ) then
          write(fstderr,*) "The flag l_fix_chi_eta_correlations must be true" &
                           // " in order to enable latin hypercube sampling" &
                           // " and ice microphysics."
          write(fstderr,*) "The flag l_ice_microphys must be set" &
                           // " to false to use this option."
          stop "Fatal error."
       endif

       allocate( l_hydromet_sed(hydromet_dim) )

       ! Sedimentation is handled within the Morrison microphysics
       l_hydromet_sed(:) = .false.

       call GRAUPEL_INIT()

    case ( "morrison_gettelman" )

       iirrm = -1
       iirsm = -1
       iirim = 1
       iirgm = -1

       iiNrm = -1
       iiNsm = -1
       iiNim = 2
       iiNgm = -1

       hydromet_dim = 2

       if ( l_predict_Nc ) then
          write(fstderr,*) "Morrison-Gettelman microphysics is not currently" &
                           // " configured for l_predict_Nc = T"
          stop "Fatal error."
       endif

       if ( l_cloud_sed ) then
          write(fstderr,*) "Morrison-Gettelman microphysics has seperate code" &
                           // " for cloud water sedimentation, therefore" &
                           // " l_cloud_sed should be set to .false."
          stop "Fatal error."
       endif

       allocate( l_hydromet_sed(hydromet_dim) )

       ! Sedimentation is handled within the MG microphysics
       l_hydromet_sed(iirim) = .false.
       l_hydromet_sed(iiNim)   = .false.

       ! Initialize constants for aerosols
       call ini_microp_aero()

       ! Setup the MG scheme
       call ini_micro()
       call pbuf_init()

    case ( "coamps" )

       if ( .not. l_predict_Nc ) then
          write(fstderr,*) "COAMPS microphysics" &
                           // " does not support l_predict_Nc = F"
          stop "Fatal Error"
       endif

       iirrm = 1
       iirsm = 2
       iirim = 3
       iirgm = 4

       iiNrm = 5
       ! Nsm is computed diagnostically in the subroutine coamps_microphys_driver
       iiNsm = -1
       iiNim = 6
       iiNgm = -1

       hydromet_dim = 6

       allocate( l_hydromet_sed(hydromet_dim) )

       l_hydromet_sed(iiNrm) = .true.
       l_hydromet_sed(iiNim) = .false.

       l_hydromet_sed(iirrm) = .true.
       l_hydromet_sed(iirsm) = .true.
       l_hydromet_sed(iirim) = .true.
       l_hydromet_sed(iirgm) = .true.

    case ( "khairoutdinov_kogan" )

       if ( l_predict_Nc ) then
          write(fstderr,*) "Khairoutdinov-Kogan microphysics" &
                           // " does not support l_predict_Nc = T"
          stop "Fatal Error"
       endif

       iirrm = 1
       iirsm = -1
       iirim = -1
       iirgm = -1

       iiNrm = 2
       iiNsm = -1
       iiNim = -1
       iiNgm = -1

       hydromet_dim = 2

       allocate( l_hydromet_sed(hydromet_dim) )

       l_hydromet_sed(iirrm) = .true.
       l_hydromet_sed(iiNrm) = .true.

    case ( "simplified_ice", "none" )

       iirrm = -1
       iirsm = -1
       iirim = -1
       iirgm = -1

       iiNrm = -1
       iiNsm = -1
       iiNim = -1
       iiNgm = -1

       hydromet_dim = 0

       l_predict_Nc = .false.

    case default

       write(fstderr,*) "Unknown microphys_scheme: "// trim( microphys_scheme )
       stop

    end select

    ! Set up predictive precipitating hydrometeor arrays.
    allocate( hydromet_list(hydromet_dim) )
    allocate( hydromet_tol(hydromet_dim) )
    allocate( l_mix_rat_hm(hydromet_dim) )
    allocate( l_frozen_hm(hydromet_dim) )
    allocate( hmp2_ip_on_hmm2_ip(hydromet_dim) )
    if ( iirrm > 0 ) then
       ! The microphysics scheme predicts rain water mixing ratio, rr.
       hydromet_list(iirrm)      = "rrm"
       l_mix_rat_hm(iirrm)       = .true.
       l_frozen_hm(iirrm)        = .false.
       hydromet_tol(iirrm)       = rr_tol
       hmp2_ip_on_hmm2_ip(iirrm) = hmp2_ip_on_hmm2_ip_ratios%rrp2_ip_on_rrm2_ip
    endif
    if ( iirim > 0 ) then
       ! The microphysics scheme predicts ice mixing ratio, ri.
       hydromet_list(iirim)      = "rim"
       l_mix_rat_hm(iirim)       = .true.
       l_frozen_hm(iirim)        = .true.
       hydromet_tol(iirim)       = ri_tol
       hmp2_ip_on_hmm2_ip(iirim) = hmp2_ip_on_hmm2_ip_ratios%rip2_ip_on_rim2_ip
    endif
    if ( iirsm > 0 ) then
       ! The microphysics scheme predicts snow mixing ratio, rs.
       hydromet_list(iirsm)      = "rsm"
       l_mix_rat_hm(iirsm)       = .true.
       l_frozen_hm(iirsm)        = .true.
       hydromet_tol(iirsm)       = rs_tol
       hmp2_ip_on_hmm2_ip(iirsm) = hmp2_ip_on_hmm2_ip_ratios%rsp2_ip_on_rsm2_ip
    endif
    if ( iirgm > 0 ) then
       ! The microphysics scheme predicts graupel mixing ratio, rg.
       hydromet_list(iirgm)      = "rgm"
       l_mix_rat_hm(iirgm)       = .true.
       l_frozen_hm(iirgm)        = .true.
       hydromet_tol(iirgm)       = rg_tol
       hmp2_ip_on_hmm2_ip(iirgm) = hmp2_ip_on_hmm2_ip_ratios%rgp2_ip_on_rgm2_ip
    endif
    if ( iiNrm > 0 ) then
       ! The microphysics scheme predicts rain drop concentration, Nr.
       hydromet_list(iiNrm)      = "Nrm"
       l_frozen_hm(iiNrm)        = .false.
       l_mix_rat_hm(iiNrm)       = .false.
       hydromet_tol(iiNrm)       = Nr_tol
       hmp2_ip_on_hmm2_ip(iiNrm) = hmp2_ip_on_hmm2_ip_ratios%Nrp2_ip_on_Nrm2_ip
    endif
    if ( iiNim > 0 ) then
       ! The microphysics scheme predicts ice concentration, Ni.
       hydromet_list(iiNim)      = "Nim"
       l_mix_rat_hm(iiNim)       = .false.
       l_frozen_hm(iiNim)        = .true.
       hydromet_tol(iiNim)       = Ni_tol
       hmp2_ip_on_hmm2_ip(iiNim) = hmp2_ip_on_hmm2_ip_ratios%Nip2_ip_on_Nim2_ip
    endif
    if ( iiNsm > 0 ) then
       ! The microphysics scheme predicts snowflake concentration, Ns.
       hydromet_list(iiNsm)      = "Nsm"
       l_mix_rat_hm(iiNsm)       = .false.
       l_frozen_hm(iiNsm)        = .true.
       hydromet_tol(iiNsm)       = Ns_tol
       hmp2_ip_on_hmm2_ip(iiNsm) = hmp2_ip_on_hmm2_ip_ratios%Nsp2_ip_on_Nsm2_ip
    endif
    if ( iiNgm > 0 ) then
       ! The microphysics scheme predicts graupel concentration, Ng.
       hydromet_list(iiNgm)      = "Ngm"
       l_mix_rat_hm(iiNgm)       = .false.
       l_frozen_hm(iiNgm)        = .true.
       hydromet_tol(iiNgm)       = Ng_tol
       hmp2_ip_on_hmm2_ip(iiNgm) = hmp2_ip_on_hmm2_ip_ratios%Ngp2_ip_on_Ngm2_ip
    endif

    select case ( trim( lh_microphys_type ) )
    case ( "interactive" )
       lh_microphys_type_int = lh_microphys_interactive

    case ( "non-interactive" )
       lh_microphys_type_int = lh_microphys_non_interactive

    case ( "disabled" )
       lh_microphys_type_int = lh_microphys_disabled

    case default
       stop "Error determining lh_microphys_type"

    end select

    ! Make sure the user didn't select LH sampling using
    ! coamps, morrison-gettelman, or simplified_ice microphysics
    if ( ( .not. ( lh_microphys_type_int == lh_microphys_disabled ) ) &
           .and. ( trim( microphys_scheme ) == "coamps" .or. &
                   trim( microphys_scheme ) == "morrison_gettelman" .or. &
                   trim( microphys_scheme ) == "simplified_ice" ) ) then
       stop "LH sampling can not be enabled when using coamps," &
            // " morrison_gettelman, or simplified_ice microphysics types"
    endif

    ! Make sure user hasn't selected l_silhs_KK_convergence_adj_mean when using
    ! a microphysics scheme other than khairoutdinov_kogan (KK)
    if ( l_silhs_KK_convergence_adj_mean .and. &
          trim( microphys_scheme ) /= "khairoutdinov_kogan" ) then
       stop "l_silhs_KK_convergence_adj_mean requires khairoutdinov_kogan microphysics"
    endif

    !The algorithm for diagnosing the correlations only works with the KK
    !microphysics by now. 
    !<Changes by janhft 02/19/13>
    if ( l_diagnose_correlations &
         .and. ( ( trim( microphys_scheme ) /= "khairoutdinov_kogan" ) &
         .and. ( lh_microphys_type_int == lh_microphys_disabled ) ) ) then
       write(fstderr,*) "Error: The diagnose_corr algorithm only works " &
                        // "for KK microphysics by now."
       stop
    endif

    if ( ( .not. l_local_kk) .and. &
         ( trim( microphys_scheme ) == "khairoutdinov_kogan" ) .and. &
         ( lh_microphys_type_int == lh_microphys_interactive ) ) then
       write(fstderr,*) "Error:  KK upscaled microphysics " &
                        // "(l_local_kk = .false.) and interactive Latin " &
                        // "Hypercube (lh_microphys_type = interactive) " &
                        // "are incompatible."
       stop
    endif

    if ( l_morr_xp2_mc .and. &
         ( lh_microphys_type_int /= lh_microphys_disabled ) ) then
       write(fstderr,*) "Error:  The code to include the effects of rain " &
                        // "evaporation on rtp2 and thlp2 in Morrison " &
                        // "microphysics (l_morr_xp2_mc = .true.) and " &
                        // "Latin Hypercube are incompatible."
       stop
    endif

    if ( l_morr_xp2_mc .and. l_var_covar_src ) then
       write(fstderr,*) "Error: The code l_morr_xp2_mc and " &
                        // "l_var_covar_src are incompatible, since " &
                        // "they both are used to determine the effect " &
                        // "of microphysics on variances."
       stop
    endif

    if ( l_morr_xp2_mc .and. l_evaporate_cold_rcm ) then
       write(fstderr,*) "Error: l_morr_xp2_mc and l_evaporate_cold_rcm " &
                        //  "are currently incompatible."
       stop
    endif

    call setup_pdf_indices( hydromet_dim, iirrm, iiNrm, &
                            iirim, iiNim, iirsm, iiNsm, &
                            iirgm, iiNgm )

    corr_file_path_cloud = corr_input_path//trim( runtype )//cloud_file_ext
    corr_file_path_below = corr_input_path//trim( runtype )//below_file_ext

    ! Allocate and set the arrays containing the correlations
    call setup_corr_varnce_array( corr_file_path_cloud, corr_file_path_below, &
                                  iunit ) ! Intent(in)

    ! Print the in-cloud and below-cloud actual (real-space) correlation arrays.
    ! This should only be done when zeta_vrnce_rat = 0.  Even when this is true,
    ! it is still possible to have a correlation that is other than these
    ! values.  The code that sets component in-precip means and variances of
    ! hydrometeors has an emergency situation (where one component mean is
    ! being pushed negative) that requires a different hm_sigma2_on_mu than
    ! hmp2_ip_on_hmm2_ip * omicron.  These printed arrays should be used as a
    ! GUIDE.  I still recommend using the GrADS or netCDF output file.
    if ( clubb_at_least_debug_level( 1 ) &
         .and. zeta_vrnce_rat == zero &
         .and. trim( microphys_scheme ) /= "none" ) then

       ! Allocate variables.
       allocate( sigma2_on_mu2_ip_cloud(d_variables) )
       allocate( sigma2_on_mu2_ip_below(d_variables) )
       allocate( sigma_x_n_cloud(d_variables) )
       allocate( sigma_x_n_below(d_variables) )
       allocate( corr_array_cloud(d_variables,d_variables) )
       allocate( corr_array_below(d_variables,d_variables) )

       ! Initialize variables.
       sigma2_on_mu2_ip_cloud(d_variables) = zero
       sigma2_on_mu2_ip_below(d_variables) = zero
       sigma_x_n_cloud = zero
       sigma_x_n_below = zero
       corr_array_cloud = zero
       corr_array_below = zero

       ! Ncn:  sigma_Ncn_i^2/mu_Ncn_i^2
       if ( .not. l_const_Nc_in_cloud ) then
          sigma2_on_mu2_ip_cloud(iiPDF_Ncn) = Ncnp2_on_Ncnm2
          sigma2_on_mu2_ip_below(iiPDF_Ncn) = Ncnp2_on_Ncnm2
       else
          sigma2_on_mu2_ip_cloud(iiPDF_Ncn) = zero
          sigma2_on_mu2_ip_below(iiPDF_Ncn) = zero
       endif
       ! Ncn:  sigma_Ncn_i_n
       sigma_x_n_cloud(iiPDF_Ncn) &
       = stdev_L2N( sigma2_on_mu2_ip_cloud(iiPDF_Ncn) )
       sigma_x_n_below(iiPDF_Ncn) = sigma_x_n_cloud(iiPDF_Ncn)

       ! Loop over all hydrometeors.
       do ivar = iiPDF_Ncn+1, d_variables, 1
          ! Hydrometeor sigma_hm_i^2/mu_hm_i^2
          sigma2_on_mu2_ip_cloud(ivar) &
          = omicron * hmp2_ip_on_hmm2_ip(pdf2hydromet_idx(ivar))
          sigma2_on_mu2_ip_below(ivar) = sigma2_on_mu2_ip_cloud(ivar)
          ! Hydrometeor sigma_hm_i_n
          sigma_x_n_cloud(ivar) = stdev_L2N( sigma2_on_mu2_ip_cloud(ivar) )
          sigma_x_n_below(ivar) = sigma_x_n_cloud(ivar)
       enddo ! i = 1, hydromet_dim, 1

       ! Calculate the correlations given the normal space correlations.
       call denorm_transform_corr( d_variables, &
                                   sigma_x_n_cloud, sigma_x_n_below, &
                                   sigma2_on_mu2_ip_cloud, &
                                   sigma2_on_mu2_ip_below, &
                                   corr_array_n_cloud, &
                                   corr_array_n_below, &
                                   corr_array_cloud, corr_array_below )

       call mirror_lower_triangular_matrix( d_variables, corr_array_cloud(:,:) )
       call mirror_lower_triangular_matrix( d_variables, corr_array_below(:,:) )

       ! Print the correlation arrays to the screen.
       write(fstdout,'(1x,A)') "Correlation array (approximate); in cloud:"
       do ivar = 1, d_variables, 1
          write(fstdout,'(12F7.3)') corr_array_cloud(ivar,:)
       enddo ! ivar = 1, d_variables, 1

       write(fstdout,'(1x,A)') "Correlation array (approximate); below cloud:"
       do ivar = 1, d_variables, 1
          write(fstdout,'(12F7.3)') corr_array_below(ivar,:)
       enddo ! ivar = 1, d_variables, 1

       ! This will open the cases setup.txt file and append it to include the
       ! parameters in the microphysics_setting namelist. This file was created
       ! and written to from clubb_driver previously.
       if ( l_write_to_file ) then

          open( unit=iunit, file=case_info_file, status='old', action='write', &
                position='append' )

          write(iunit,'(1x,A)') "Correlation array (approximate); in cloud:"
          do ivar = 1, d_variables, 1
             write(iunit,'(12F7.3)') corr_array_cloud(ivar,:)
          enddo ! ivar = 1, d_variables, 1

          write(iunit,'(1x,A)') "Correlation array (approximate); below cloud:"
          do ivar = 1, d_variables, 1
             write(iunit,'(12F7.3)') corr_array_below(ivar,:)
          enddo ! ivar = 1, d_variables, 1

          close( unit=iunit )

       endif ! l_write_to_file

       ! Deallocate variables.
       deallocate( sigma2_on_mu2_ip_cloud )
       deallocate( sigma2_on_mu2_ip_below )
       deallocate( sigma_x_n_cloud )
       deallocate( sigma_x_n_below )
       deallocate( corr_array_cloud )
       deallocate( corr_array_below )

    endif ! clubb_at_least_debug_level( 1 )
          ! and zeta_vrnce_rat = 0
          ! and microphys_scheme /= "none"


    return

  end subroutine init_microphys

  !=============================================================================
  subroutine cleanup_microphys( )

    ! Description:
    ! De-allocate arrays used by the microphysics

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use parameters_microphys, only: &
        l_hydromet_sed, &
        microphys_scheme

    use array_index, only: &
        hydromet_list,  & ! Variable(s)
        hydromet_tol,   &
        l_mix_rat_hm, & ! Variable(s)
        l_frozen_hm

    use corr_varnce_module, only: &
        hmp2_ip_on_hmm2_ip

    use phys_buffer, only: & ! Used for placing wp2_zt in MG microphys.
        pbuf_deallocate

    implicit none

    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( hydromet_list ) ) then
       deallocate( hydromet_list )
    endif

    if ( allocated( l_hydromet_sed ) ) then
       deallocate( l_hydromet_sed )
    endif

    if ( allocated( l_mix_rat_hm ) ) then
       deallocate( l_mix_rat_hm )
    endif

    if ( allocated( l_frozen_hm ) ) then
       deallocate( l_frozen_hm )
    endif

    if ( allocated( hydromet_tol ) ) then
       deallocate( hydromet_tol )
    endif

    if ( allocated( hmp2_ip_on_hmm2_ip ) ) then
       deallocate( hmp2_ip_on_hmm2_ip )
    endif

    if ( trim( microphys_scheme ) == "morrison_gettelman" ) then
       call pbuf_deallocate()
    endif


    return

  end subroutine cleanup_microphys

  !=============================================================================

end module microphys_init_cleanup
