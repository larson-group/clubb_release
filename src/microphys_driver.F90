!-------------------------------------------------------------------------------
! $Id$
module microphys_driver

! Description:
!   Call a microphysical scheme to compute hydrometeor species, and advect,
!   sediment, & diffuse using a tridiagonal system.

! References:
!  H. Morrison, J. A. Curry, and V. I. Khvorostyanov, 2005: A new double-
!  moment microphysics scheme for application in cloud and
!  climate models. Part 1: Description. J. Atmos. Sci., 62, 1665-1677.

! Khairoutdinov, M. and Kogan, Y.: A new cloud physics parameteri-
! zation in a large-eddy simulation model of marine stratocumulus,
! Mon. Wea. Rev., 128, 229-243, 2000. 
!-------------------------------------------------------------------------------
  use parameters_microphys, only: &
    l_cloud_sed,                & ! Cloud water sedimentation (K&K or no microphysics)
    l_ice_micro,                & ! Compute ice (COAMPS / Morrison)
    l_graupel,                  & ! Compute graupel (Morrison)
    l_hail,                     & ! See module_mp_graupel for a description
    l_seifert_beheng,           & ! Use Seifert and Beheng (2001) warm drizzle (Morrison)
    l_predictnc,                & ! Predict cloud droplet number conc (Morrison)
    l_specify_aerosol,          & ! Specify aerosol (Morrison)
    l_subgrid_w,                & ! Use subgrid w  (Morrison)
    l_arctic_nucl,              & ! Use MPACE observations (Morrison)
    l_cloud_edge_activation,    & ! Activate on cloud edges (Morrison)
    l_fix_pgam,                 & ! Fix pgam (Morrison)
    l_latin_hypercube_sampling, & ! Use Latin Hypercube Sampling (K&K only)
    LH_microphys_calls,         & ! Number of latin hypercube samples to call the microphysics with 
    LH_sequence_length,         & ! Number of timesteps before the latin hypercube seq. repeats
    l_local_kk,                 & ! Use local formula for K&K
    micro_scheme,               & ! The microphysical scheme in use
    hydromet_list,              & ! Names of the hydrometeor species
    microphys_start_time,       & ! When to start the microphysics [s]
    Ncm_initial                   ! Initial value for Ncm (K&K, l_cloud_sed, Morrison)

  implicit none

  ! Subroutines
  public :: advance_microphys, init_microphys

  private :: microphys_lhs, microphys_solve
  private :: adj_microphys_tndcy

  ! Functions
  private :: sedimentation

  ! Variables
  logical, private, allocatable, dimension(:) :: hydromet_sed ! Whether to sediment variables

  private ! Default Scope

  contains

!-------------------------------------------------------------------------------
  subroutine init_microphys( iunit, namelist_file, hydromet_dim )

! Description:
!   Set indices to the various hydrometeor species and define hydromet_dim for 
!   the purposes of allocating memory.

! References:
!   None
!-------------------------------------------------------------------------------

    use array_index, only: & 
      iirrainm, iiNrm, iirsnowm, iiricem, iirgraupelm, & ! Variables
      iiNcm, iiNsnowm, iiNim, iiNgraupelm

    use array_index, only: &
      iiLH_rt, & ! Variables
      iiLH_thl, &
      iiLH_w

    use array_index, only: &
      iiLH_rrain, & ! Variables
      iiLH_rsnow, &
      iiLH_rice, &
      iiLH_rgraupel, &
      iiLH_Nr, &
      iiLH_Nsnow, &
      iiLH_Ni, &
      iiLH_Ngraupel, &
      iiLH_Nc

    use KK_microphys_module, only: &
      rrp2_rrainm2_cloud, Nrp2_Nrm2_cloud, Ncp2_Ncm2_cloud, & ! Variables
      corr_rrNr_LL_cloud, corr_srr_NL_cloud, corr_sNr_NL_cloud, &
      corr_sNc_NL_cloud, rrp2_rrainm2_below, &
      Nrp2_Nrm2_below, Ncp2_Ncm2_below, &
      corr_rrNr_LL_below, corr_srr_NL_below, &
      corr_sNr_NL_below, corr_sNc_NL_below, &
      C_evap, r_0

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
      dospecifyaerosol, &   ! specify two modes of (sulfate) aerosol
      dosubgridw, &         ! input estimate of subgrid w to microphysics
      doarcticicenucl, &    ! use arctic parameter values for ice nucleation
      docloudedgeactivation,& ! activate cloud droplets throughout the cloud
      dofix_pgam            ! option to fix value of pgam (exponent in cloud water gamma distn)

    use module_mp_Graupel, only: &
      GRAUPEL_INIT ! Subroutine

    use constants, only: &
      fstderr,   & ! Constant
      cm3_per_m3

    implicit none

    ! External
    intrinsic :: trim

    ! Input variables
    integer, intent(in) :: &
      iunit ! File unit

    character(len=*), intent(in) :: &
      namelist_file ! File name

    ! Output variables
    integer, intent(out) :: & 
      hydromet_dim ! Number of hydrometeor fields.

    namelist /microphysics_setting/ &
      micro_scheme, l_cloud_sed, &
      l_cloud_sed, l_ice_micro, l_graupel, l_hail, &
      l_seifert_beheng, l_predictnc, l_specify_aerosol, l_subgrid_w, &
      l_arctic_nucl, l_cloud_edge_activation, l_fix_pgam, &
      l_latin_hypercube_sampling, LH_microphys_calls, LH_sequence_length, &
      rrp2_rrainm2_cloud, Nrp2_Nrm2_cloud, Ncp2_Ncm2_cloud, &
      corr_rrNr_LL_cloud, corr_srr_NL_cloud, corr_sNr_NL_cloud, &
      corr_sNc_NL_cloud, rrp2_rrainm2_below, &
      Nrp2_Nrm2_below, Ncp2_Ncm2_below, &
      corr_rrNr_LL_below, corr_srr_NL_below, &
      corr_sNr_NL_below, corr_sNc_NL_below, &
      C_evap, r_0, microphys_start_time, &
      Ncm_initial, ccnconst, ccnexpnt, aer_rm1, aer_rm2, &
      aer_n1, aer_n2, aer_sig1, aer_sig2, pgam_fixed

    ! Local variables
    integer :: i

    ! ---- Begin Code ----

    ! Set default values, then read in the namelist
    micro_scheme = "none"

    ! Cloud water sedimentation from the RF02 case
    ! This has no effect on Morrison's cloud water sedimentation
    l_cloud_sed = .false.

    !---------------------------------------------------------------------------
    ! Parameters for Khairoutdinov and Kogan microphysics 
    !---------------------------------------------------------------------------

    ! Parameters for in-cloud (from SAM RF02 DO).
    rrp2_rrainm2_cloud = 0.766
    Nrp2_Nrm2_cloud    = 0.429
    Ncp2_Ncm2_cloud    = 0.003
    corr_rrNr_LL_cloud = 0.786
    corr_srr_NL_cloud  = 0.242
    corr_sNr_NL_cloud  = 0.285
    corr_sNc_NL_cloud  = 0.433
    ! Parameters for below-cloud (from SAM RF02 DO).
    rrp2_rrainm2_below = 8.97
    Nrp2_Nrm2_below    = 12.03
    Ncp2_Ncm2_below    = 0.00 ! Not applicable below cloud.
    corr_rrNr_LL_below = 0.886
    corr_srr_NL_below  = 0.056
    corr_sNr_NL_below  = 0.015
    corr_sNc_NL_below  = 0.00 ! Not applicable below cloud.
    ! Other needed parameters
    C_evap = 0.86    ! Khairoutdinov and Kogan (2000) ratio of
    ! drizzle drop mean geometric radius to
    ! drizzle drop mean volume radius.
    ! Khairoutdinov and Kogan (2000); p. 233.
    !C_evap = 0.86*0.2 ! COAMPS value of KK C_evap
    !C_evap = 0.55     ! KK 2000, Marshall-Palmer (1948) value.

    r_0 = 25.0e-6   ! Assumed radius of all new drops; m.

    !---------------------------------------------------------------------------
    ! Parameters for Morrison and COAMPS microphysics 
    !---------------------------------------------------------------------------
    l_ice_micro = .true.
    l_graupel = .true.

    !---------------------------------------------------------------------------
    ! Parameters for Morrison microphysics only
    !---------------------------------------------------------------------------
    l_hail = .false.
    l_seifert_beheng = .false.
    l_predictnc = .true.
    l_specify_aerosol = .true.
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

    l_latin_hypercube_sampling = .false.
    LH_microphys_calls = 2
    LH_sequence_length = 1

    !---------------------------------------------------------------------------
    ! Parameters for all microphysics schemes
    !---------------------------------------------------------------------------
    microphys_start_time = 0.0 ! [s]

    open(unit=iunit, file=namelist_file, status='old',action='read')
    read(iunit, nml=microphysics_setting)
    close(unit=iunit)

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

      if ( l_specify_aerosol ) then
        dospecifyaerosol = .true.
      else
        dospecifyaerosol = .false.
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

      allocate( hydromet_sed(hydromet_dim) )
      ! Sedimentation is handled within the Morrison microphysics
      hydromet_sed(iiNrm) = .false.
      hydromet_sed(iiNim) = .false.
      hydromet_sed(iiNcm) = .false.
      hydromet_sed(iiNgraupelm) = .false.
      hydromet_sed(iiNsnowm) = .false.

      hydromet_sed(iirrainm)    = .false.
      hydromet_sed(iirsnowm)    = .false.
      hydromet_sed(iiricem)     = .false.
      hydromet_sed(iirgraupelm) = .false.

      ! Convert from μ to m as in SAM
      aer_rm1 = 1.e-6*aer_rm1 
      aer_rm2 = 1.e-6*aer_rm2 
      ! Convert from #/cm3 to #/m3
      aer_n1 = cm3_per_m3 * aer_n1 
      aer_n2 = cm3_per_m3 * aer_n2

      ! Setup the Morrison scheme
      call GRAUPEL_INIT()

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

      allocate( hydromet_sed(hydromet_dim) )

      allocate( hydromet_list(hydromet_dim) )

      hydromet_list(iirrainm)    = "rrainm"
      hydromet_list(iirsnowm)    = "rsnowm"
      hydromet_list(iiricem)     = "ricem"
      hydromet_list(iirgraupelm) = "rgraupelm"

      hydromet_list(iiNrm)       = "Nrm"
      hydromet_list(iiNcm)       = "Ncm"
      hydromet_list(iiNim)       = "Nim"

      hydromet_sed(iiNrm) = .true.
      hydromet_sed(iiNim) = .false.
      hydromet_sed(iiNcm) = .false.

      hydromet_sed(iirrainm)    = .true.
      hydromet_sed(iirsnowm)    = .true.
      hydromet_sed(iiricem)     = .true.
      hydromet_sed(iirgraupelm) = .true.

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

      allocate( hydromet_sed(hydromet_dim) )

      allocate( hydromet_list(hydromet_dim) )

      hydromet_list(iirrainm) = "rrainm"
      hydromet_list(iiNrm) = "Nrm"
      hydromet_list(iiNcm) = "Ncm"

      hydromet_sed(iirrainm) = .true.
      hydromet_sed(iiNrm)    = .true.
      hydromet_sed(iiNcm)    = .false.

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

        allocate( hydromet_sed(hydromet_dim) )
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

    ! Setup index variables for latin hypercube sampling
    if ( l_latin_hypercube_sampling ) then

      iiLH_rt  = 1
      iiLH_thl = 2
      iiLH_w   = 3

      i = iiLH_w

      call return_LH_index( iiNcm, i, iiLH_Nc )
      call return_LH_index( iiNrm, i, iiLH_Nr )
      call return_LH_index( iiNsnowm, i, iiLH_Nsnow )
      call return_LH_index( iiNim, i, iiLH_Ni )
      call return_LH_index( iiNgraupelm, i, iiLH_Ngraupel )

      call return_LH_index( iirrainm, i, iiLH_rrain )
      call return_LH_index( iirsnowm, i, iiLH_rsnow )
      call return_LH_index( iiricem, i, iiLH_rice )
      call return_LH_index( iirgraupelm, i, iiLH_rgraupel )

    end if

    return
  end subroutine init_microphys

!-------------------------------------------------------------------------------
  subroutine advance_microphys & 
             ( iter, runtype, dt, time_current,  & 
               thlm, p_in_Pa, exner, rho, rho_zm, rtm, rcm, cf, & 
               wm_zt, wm_zm, Kh_zm, pdf_params, & 
               wp2_zt, &
               Ncnm, hydromet, & 
               rtm_forcing, thlm_forcing, &
               err_code )

! Description:
!   Compute pristine ice, snow, graupel, & rain hydrometeor fields.
!   Uses implicit discretization.

! References:
!   None
!-------------------------------------------------------------------------------

    use grid_class, only: & 
        gr,  & ! Variable(s)
        zm2zt,  & ! Procedure(s)
        zt2zm

    use KK_microphys_module, only: & 
        KK_microphys ! Procedure(s)

    use KK_microphys_module, only: & 
      rrp2_rrainm2_cloud, & ! Variable(s)
      Nrp2_Nrm2_cloud, & 
      Ncp2_Ncm2_cloud  

    use morrison_micro_driver_mod, only: &
      morrison_micro_driver

    use latin_hypercube_mod, only: &
      latin_hypercube_driver ! Procedure

    use ice_dfsn_mod, only: & 
        ice_dfsn ! Procedure(s)

    use parameters_tunable, only: & 
        c_Krrainm,  & ! Variable(s) 
        nu_r

    use parameters_model, only: & 
        hydromet_dim   ! Integer

    use constants, only:  & 
        Lv,   & ! Constant(s)
        Cp,  & 
        rho_lw,  & 
        fstderr, & 
        zero_threshold, &
        sec_per_day

    use stats_precision, only:  & 
        time_precision ! Variable(s)

    use error_code, only:  & 
        lapack_error, & ! Procedure
        clubb_at_least_debug_level

    use coamps_micro_driver_mod, only:  & 
        coamps_micro_driver ! Procedure

    use variables_prognostic_module, only:  &
        pdf_parameter  ! type

    use array_index, only:  & 
        iirrainm, iirsnowm, iiricem, iirgraupelm, &
        iiNrm, iiNsnowm, iiNim, iiNgraupelm, iiNcm

    use stats_variables, only: & 
      iVrr,  & ! Variable(s)
      iVNr, & 
      iVsnow, & 
      iVice, & 
      iVgraupel, & 
      ithlm_mc, & 
      irtm_mc, & 
      irain_rate, & 
      iFprec, & 
      irrainm_bt, & 
      irrainm_mc, & 
      irrainm_cond_adj, & 
      irrainm_cl, & 
      iNrm_bt, & 
      iNrm_mc, & 
      iNrm_cond_adj, & 
      iNrm_cl, & 
      iNcnm, & 
      irrainm, & 
      irsnowm, & 
      iricem, & 
      irgraupelm, & 
      iricem_bt, & 
      iricem_mc, & 
      iricem_cl, & 
      irgraupelm_bt, & 
      irgraupelm_mc, & 
      irgraupelm_cl, & 
      irsnowm_bt, & 
      irsnowm_mc, & 
      irsnowm_cl, & 
      irain, & 
      irain_flux, & 
      irrainm_sfc

    use stats_variables, only: & 
      iNcm_bt, & 
      iNcm_cl

    use stats_variables, only: & 
      iNcm, & 
      iNim, & 
      iNrm, & 
      iNsnowm, &
      iNgraupelm

!   use stats_variables, only: & 
!     iNim_mc, & 
!     iNcm_mc, & 
!     iNrm_mc, & 
!     iNsnowm_mc, &
!     iNgraupelm_mc

    use stats_variables, only: & 
      iLH_rcm_mc, &
      iLH_rvm_mc, &
      iLH_thlm_mc, &
      iLH_rrainm_mc, &
      iLH_Nrm_mc

    use stats_variables, only: & 
      iLH_Vrr, &
      iLH_VNr

    use stats_variables, only: & 
        zt, &  ! Variables
        zm, & 
        sfc, & 
        l_stats_samp

    use stats_type, only:  & 
        stat_update_var, stat_update_var_pt,  & ! Procedure(s)
        stat_begin_update, stat_end_update

    implicit none

    ! Input Variables

    integer, intent(in) :: iter ! Model iteration number

    character(len=*), intent(in) :: & 
      runtype ! Name of the run, for case specific effects.

    real(kind=time_precision), intent(in) ::  & 
      dt           ! Timestep         [s]

    real(kind=time_precision), intent(in) ::  & 
      time_current ! Current time     [s]

    real, dimension(gr%nnzp), intent(in) :: & 
      thlm,    & ! Liquid potential temp.                 [K]
      p_in_Pa, & ! Pressure                               [Pa]
      exner,   & ! Exner function                         [-]
      rho,     & ! Density on thermo. grid                [kg/m^3]
      rho_zm,  & ! Density on moment. grid                [kg/m^3]
      rtm,     & ! Total water mixing ratio               [kg/kg]
      rcm,     & ! Liquid water mixing ratio              [kg/kg]
      cf,      & ! Cloud fraction                         [%]
      wm_zt,   & ! w wind on moment. grid                 [m/s]
      wm_zm,   & ! w wind on thermo. grid                 [m/s]
      Kh_zm      ! Kh Eddy diffusivity on momentum grid   [m^2/s]

    type(pdf_parameter), intent(in) :: & 
      pdf_params     ! PDF parameters

    real, dimension(gr%nnzp), intent(in) :: & 
      wp2_zt ! w'^2 on the thermo. grid         [m^2/s^2]

    ! Note:
    ! K & K only uses Ncm, while for COAMPS Ncnm is initialized
    ! and Nim & Ncm are computed within subroutine adjtg.
    real, dimension(gr%nnzp), intent(inout) :: & 
      Ncnm       ! Cloud nuclei number concentration     [count/m^3]

    real, dimension(gr%nnzp,hydromet_dim), intent(inout) :: & 
      hydromet      ! Array of rain, prist. ice, graupel, etc. [units vary]

    real, dimension(gr%nnzp), intent(inout) :: & 
      rtm_forcing,  & ! Imposed contributions to total water        [kg/kg/s]
      thlm_forcing    ! Imposed contributions to liquid potential temp. [K/s]

    integer, intent(out) :: err_code ! Exit code returned from subroutine

    ! Local Variables
    real, dimension(3,gr%nnzp) :: & 
      lhs ! Left hand side of tridiagonal matrix

    real, dimension(gr%nnzp,hydromet_dim) :: &
      hydromet_vel ! Contains vel. of the hydrometeors        [m/s]
    !  Can contain:
    !  Vrr      ! Rain mixing ratio sedimentation velocity     [m/s]
    !  VNr      ! Rain number conc. sedimentation velocity     [m/s]
    !  Vice     ! Ice mixing ratio sedimentation velocity      [m/s]
    !  Vsnow    ! Snow mixing ratio sedimentation velocity     [m/s]
    !  Vgraupel ! Graupel mixing ratio sedimentation velocity  [m/s]

    real, dimension(gr%nnzp) :: & 
      rtm_mc,   & ! Change in total water due to microphysics   [kg/kg/s]
      rvm_mc,   & ! Change in vapor water due to microphysics   [kg/kg/s]
      rcm_mc,   & ! Change in liquid water due to microphysics  [kg/kg/s]
      thlm_mc     ! Change in liquid potential temperature due to microphysics [K/s]

    real, dimension(gr%nnzp,hydromet_dim) :: & 
      hydromet_mc  ! Change in hydrometeors due to microphysics  [units/s]
      
    real, dimension(hydromet_dim,hydromet_dim) :: & 
      hydromet_corr  ! Array for correlations   [-]

    real, dimension(gr%nnzp) :: &
      dzq         ! Difference in height levels                   [m]

!   real, dimension(gr%nnzp) :: & 
!     T_in_K  ! Temperature          [K]
      
    real, dimension(1,1,gr%nnzp) :: & 
      cond ! COAMPS stat for condesation/evap of rcm

    ! Eddy diffusivity for rain and rain drop concentration.
    ! It is also used for the other hydrometeor variables.
    ! Kr = Constant * Kh_zm; Constant is named c_Krrainm.
    real, dimension(gr%nnzp) :: Kr   ! [m^2/s]

    ! Variable needed to handle correction to rtm and thlm microphysics
    ! tendency arrays, as well rrainm_cond and Nrm_cond statistical
    ! tendency arrays, due to a negative result being produced by
    ! over-evaporation of rain water over the course of a timestep.
    ! Brian Griffin.  April 14, 2007.
    real :: overevap_rate ! Absolute value of negative evap. rate.

    real, dimension(gr%nnzp) :: &
      wtmp   ! Standard dev. of w       [m/s]

    integer :: i, k ! Array index

    integer :: d_variables

    integer :: ixrm_cl, ixrm_bt

!-------------------------------------------------------------------------------

    ! ---- Begin code ----

    ! Return if there is delay between the model start time and start of the
    ! microphysics
    if ( time_current < microphys_start_time ) return

    ! Solve for the value of Kr, the hydrometeor eddy diffusivity.
    do k = 1, gr%nnzp, 1
      Kr(k) = c_Krrainm * Kh_zm(k)
    end do

    ! Determine temperature in K for the microphysics

    ! Compute standard deviation of vertical velocity in the grid column
    wtmp(:) = sqrt( wp2_zt(:) )

    ! Compute difference in height levels
    dzq(1:gr%nnzp) = 1./gr%dzm(1:gr%nnzp)

   ! Begin by calling Brian Griffin's implementation of the
   ! Khairoutdinov and Kogan microphysical scheme, 
   ! the Rutlege and Hobbes scheme from COAMPS, or the Morrison scheme.
   ! Note: COAMPS appears to have some K&K elements to it as well.

    select case ( trim( micro_scheme ) )

    case ( "coamps" ) 

      ! Initialize tendencies to zero
      hydromet_mc(:,:) = 0.0

      call coamps_micro_driver & 
           ( runtype, time_current, dt, & 
             rtm, wm_zm, p_in_Pa, exner, rho, & 
             thlm, hydromet(:,iiricem), hydromet(:,iirrainm),  & 
             hydromet(:,iirgraupelm), hydromet(:,iirsnowm), & 
             rcm, hydromet(:,iiNcm), hydromet(:,iiNrm), Ncnm, hydromet(:,iiNim), cond, & 
             hydromet_vel(:,iirsnowm), hydromet_vel(:,iiricem), & 
             hydromet_vel(:,iirrainm), hydromet_vel(:,iiNrm),  & 
             hydromet_vel(:,iirgraupelm), & 
             hydromet_mc(:,iiricem), hydromet_mc(:,iirrainm), & 
             hydromet_mc(:,iirgraupelm), hydromet_mc(:,iirsnowm), & 
             hydromet_mc(:,iiNrm), & 
             rtm_mc, thlm_mc )

      if ( l_stats_samp ) then

        ! Sedimentation velocity for rrainm
        call stat_update_var(iVrr, hydromet_vel(:,iirrainm), zm)

        ! Sedimentation velocity for Nrm
        call stat_update_var(iVNr, hydromet_vel(:,iiNrm), zm )

        ! Sedimentation velocity for snow
        call stat_update_var(iVsnow, hydromet_vel(:,iirsnowm), zm )

        ! Sedimentation velocity for pristine ice
        call stat_update_var( iVice, hydromet_vel(:,iiricem), zm )

        ! Sedimentation velocity for graupel
        call stat_update_var( iVgraupel,  & 
                              hydromet_vel(:,iirgraupelm), zm )

        ! Sum total of rrainm microphysics
        call stat_update_var( irrainm_mc, hydromet_mc(:,iirrainm), zt )

        ! Sum total of Nrm microphysics
        call stat_update_var( iNrm_mc, hydromet_mc(:,iiNrm), zt )

        ! Sum total of pristine ice microphysics
        call stat_update_var( iricem_mc, hydromet_mc(:,iiricem), zt )

        ! Sum total of graupel microphysics
        call stat_update_var( irgraupelm_mc, & 
                              hydromet_mc(:,iirgraupelm), zt )

        ! Sum total of snow microphysical processeses
        call stat_update_var( irsnowm_mc,  & 
                              hydromet_mc(:,iirsnowm), zt )

      end if ! l_stats_samp

    case ( "morrison" )

      ! Initialize tendencies to zero
      hydromet_mc(:,:) = 0.0
      rcm_mc(:) = 0.0
      rvm_mc(:) = 0.0
      thlm_mc(:) = 0.0

      if ( l_latin_hypercube_sampling ) then

        d_variables = 3 + hydromet_dim
!       d_variables = 5

        ! For latin hypercube sampling
        hydromet_corr(:,:) = 0.0 ! Initialize to 0
        hydromet_corr(iiNcm,iiNcm)       = Ncp2_Ncm2_cloud
        hydromet_corr(iiNrm,iiNrm)       = Nrp2_Nrm2_cloud
        hydromet_corr(iirrainm,iirrainm) = rrp2_rrainm2_cloud

        call latin_hypercube_driver &
             ( real( dt ), iter, d_variables, LH_microphys_calls, &
               LH_sequence_length, gr%nnzp, &
               cf, thlm, p_in_Pa, exner, &
               rho, pdf_params, wm_zt, wtmp, dzq, rcm, rtm-rcm, &
               hydromet, hydromet_corr, hydromet_mc, hydromet_vel, rcm_mc, &
               rvm_mc, thlm_mc, morrison_micro_driver )

        if ( l_stats_samp ) then

          ! Latin hypercube estimate for cloud water mixing ratio microphysical tendency
          call stat_update_var( iLH_rcm_mc, rcm_mc, zt )

          ! Latin hypercube estimate for vapor water mixing ratio microphysical tendency
          call stat_update_var( iLH_rvm_mc, rvm_mc, zt )

          ! Latin hypercube estimate for rain water mixing ratio microphysical tendency
          call stat_update_var( iLH_rrainm_mc, hydromet_mc(:,iirrainm), zt )

          ! Latin hypercube estimate for rain water number concentration microphysical tendency
          call stat_update_var( iLH_Nrm_mc, hydromet_mc(:,iiNrm), zt )

          ! Latin hypercube estimate for liquid potential temperature
          call stat_update_var( iLH_thlm_mc, thlm_mc, zt )

        end if
      end if ! l_latin_hypercube_sampling

      ! Based on YSU PBL interface to the Morrison scheme WRF driver, the standard dev. of w
      ! will be clipped to be between 0.1 m/s and 4.0 m/s in WRF.  -dschanen 23 Mar 2009
!     wtmp(:) = max( 0.1, wtmp ) ! Disabled for now
!     wtmp(:) = min( 4., wtmp )  

!     wtmp = 0.5 ! %% debug
      call morrison_micro_driver & 
           ( real( dt ), gr%nnzp, l_stats_samp, .false., thlm, p_in_Pa, exner, rho, pdf_params, &
             wm_zt, wtmp, dzq, rcm, rtm-rcm, hydromet, hydromet_mc, &
             hydromet_vel, rcm_mc, rvm_mc, thlm_mc )
 
      ! Update total water tendency
      rtm_mc = rcm_mc + rvm_mc
 
    case ( "khairoutdinov_kogan" )

      ! Initialize tendencies to zero
      hydromet_mc(:,:) = 0.0
      rcm_mc(:) = 0.0
      rvm_mc(:) = 0.0
      thlm_mc(:) = 0.0

      if ( l_latin_hypercube_sampling ) then

        d_variables = 3 + hydromet_dim
!       d_variables = 5

        ! For latin hypercube sampling
        hydromet_corr(:,:) = 0.0 ! Initialize to 0
        hydromet_corr(iiNcm,iiNcm)       = Ncp2_Ncm2_cloud
        hydromet_corr(iiNrm,iiNrm)       = Nrp2_Nrm2_cloud
        hydromet_corr(iirrainm,iirrainm) = rrp2_rrainm2_cloud

        call latin_hypercube_driver &
             ( real( dt ), iter, d_variables, LH_microphys_calls, &
               LH_sequence_length, gr%nnzp, &
               cf, thlm, p_in_Pa, exner, &
               rho, pdf_params, wm_zt, wtmp, dzq, rcm, rtm-rcm, &
               hydromet, hydromet_corr, hydromet_mc, hydromet_vel, rcm_mc, &
               rvm_mc, thlm_mc, KK_microphys )

        if ( l_stats_samp ) then

          ! Latin hypercube estimate for cloud water mixing ratio microphysical tendency
          call stat_update_var( iLH_rcm_mc, rcm_mc, zt )

          ! Latin hypercube estimate for vapor water mixing ratio microphysical tendency
          call stat_update_var( iLH_rvm_mc, rvm_mc, zt )

          ! Latin hypercube estimate for rain water mixing ratio microphysical tendency
          call stat_update_var( iLH_rrainm_mc, hydromet_mc(:,iirrainm), zt )

          ! Latin hypercube estimate for rain water number concentration microphysical tendency
          call stat_update_var( iLH_Nrm_mc, hydromet_mc(:,iiNrm), zt )

          ! Latin hypercube estimate for liquid potential temperature
          call stat_update_var( iLH_thlm_mc, thlm_mc, zt )

          ! Latin hypercube estimate for sedimentation velocities
          call stat_update_var( iLH_Vrr, zt2zm( hydromet_vel(:,iirrainm) ), zt )

          call stat_update_var( iLH_VNr, zt2zm( hydromet_vel(:,iiNrm) ), zt )

        end if

      end if ! l_latin_hypercube_sampling

      call KK_microphys & 
           ( real( dt ), gr%nnzp, l_stats_samp, l_local_kk, thlm, p_in_Pa, exner, rho, pdf_params, &
             wm_zt, wtmp, dzq, rcm, rtm-rcm, hydromet, hydromet_mc, &
             hydromet_vel, rcm_mc, rvm_mc, thlm_mc )

      ! Interpolate velocity to the momentum grid
      ! Note that the interpolation functions are written such that we can
      ! safely have hydromet_vel reference itself here.
      do i = 1, hydromet_dim
        hydromet_vel(:,i) = zt2zm( hydromet_vel(:,i) )
      end do
      hydromet_vel(gr%nnzp,1:hydromet_dim) = 0.0

      rtm_mc = rcm_mc + rvm_mc

      if ( l_stats_samp ) then
        ! Sedimentation velocity for rrainm
        call stat_update_var( iVrr, hydromet_vel(:,iirrainm), zm )

        ! Sedimentation velocity for Nrm
        call stat_update_var( iVNr, hydromet_vel(:,iiNrm), zm )
      end if

    case default

    end select ! micro_scheme

!-----------------------------------------------------------------------
!       Loop over all hydrometeor species and apply sedimentation,
!       advection and diffusion.
!-----------------------------------------------------------------------
    if ( hydromet_dim > 0 ) then

      do i = 1, hydromet_dim

        select case( trim( hydromet_list(i) ) )
        case( "rrainm" )
          ixrm_bt = irrainm_bt
          ixrm_cl = irrainm_cl
        case( "Nrm" )
          ixrm_bt = iNrm_bt
          ixrm_cl = iNrm_cl
        case( "Ncm" )
          ixrm_bt = iNcm_bt
          ixrm_cl = iNcm_cl
        case( "ricem" )
          ixrm_bt = iricem_bt
          ixrm_cl = iricem_cl
        case( "rsnowm" )
          ixrm_bt = irsnowm_bt
          ixrm_cl = irsnowm_cl
        case( "rgraupelm" )
          ixrm_bt = irgraupelm_bt
          ixrm_cl = irgraupelm_cl
        case default
          ixrm_bt = 0
          ixrm_cl = 0
        end select

        if ( l_stats_samp ) then
          call stat_begin_update & 
             ( ixrm_bt, real(hydromet(:,i) / dt), zt )
        end if



        call microphys_lhs & 
             ( trim( hydromet_list(i) ), hydromet_sed(i), dt, Kr, nu_r, wm_zt, & 
               hydromet_vel(:,i), lhs )

        call microphys_solve & 
             ( trim( hydromet_list(i) ), hydromet_sed(i), dt, lhs, hydromet_mc(:,i),  & 
               hydromet(:,i), err_code )

        if ( i == iirrainm ) then
          ! Handle over-evaporation of rrainm and adjust rt and theta-l
          ! hydrometeor tendency arrays accordingly.
          do k = 1, gr%nnzp, 1
            if ( hydromet(k,i) < 0.0 ) then

              call adj_microphys_tndcy & 
                 ( hydromet_mc(:,i), wm_zt, hydromet_vel(:,i),  & 
                   Kr, nu_r, dt, k, .true., & 
                   hydromet(:,i), overevap_rate )

              ! overevap_rate is defined as positive.
              ! It is a correction factor.
              rtm_mc(k)  = rtm_mc(k) - overevap_rate

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

          end do ! k=1..gr%nnzp

        else if ( i == iiNrm ) then
          ! Handle over-evaporation similar to rrainm.  However, in the case
          ! of Nrm there is no effect on rtm or on thlm.
          ! Brian Griffin.  April 14, 2007.
          do k = 1, gr%nnzp, 1
            if ( hydromet(k,i) < 0.0 ) then

              call adj_microphys_tndcy & 
                   ( hydromet_mc(:,i), wm_zt, hydromet_vel(:,i),  & 
                     Kr, nu_r, dt, k, .true., & 
                     hydromet(:,i), overevap_rate )

              ! Moved from adj_microphys_tndcy
              if( l_stats_samp ) then 
                 call stat_update_var_pt( iNrm_cond_adj, k,  & 
                                       overevap_rate, zt )
              endif
            else
              if( l_stats_samp ) then
                call stat_update_var_pt( iNrm_cond_adj,k, 0.0, zt )
              end if

            end if ! ! Nrm(k) < 0
            ! Joshua Fasching December 2007
          end do

        end if

        if ( l_stats_samp ) then

          call stat_begin_update & 
             ( ixrm_cl, real( hydromet(:,i) / dt ), zt )

        end if

        ! Clip all hydrometeor species to be > zero
        where ( hydromet(:,i) < zero_threshold ) hydromet(:,i) = zero_threshold

          if ( l_stats_samp ) then

            ! Effects of clipping
            call stat_end_update( ixrm_cl, real( hydromet(:,i) / dt ), zt )

            ! Total time tendency
            call stat_end_update( ixrm_bt, real(hydromet(:,i) / dt), zt )

          end if ! l_stats_samp

        end do ! i=1..hydromet_dim

      end if ! hydromet_dim > 0


      ! Call the ice diffusion scheme
      if ( trim( micro_scheme ) == "simplified_ice" ) then
        call ice_dfsn( dt, thlm, rcm, exner, p_in_Pa, rho, rtm_mc, thlm_mc )
      end if

      if ( l_stats_samp ) then

        if ( iiNcm > 0 ) then
          call stat_update_var( iNcm, hydromet(:,iiNcm), zt )
        end if

        call stat_update_var( iNcnm, Ncnm, zt )

        if ( iirrainm > 0 ) then
          call stat_update_var( irrainm, hydromet(:,iirrainm), zt )
        end if

        if ( iirsnowm > 0 ) then
          call stat_update_var( irsnowm, hydromet(:,iirsnowm), zt )
        end if

        if ( iiricem > 0 ) then
          call stat_update_var( iricem, hydromet(:,iiricem), zt )
        end if

        if ( iirgraupelm > 0 ) then
          call stat_update_var( irgraupelm,  & 
                                hydromet(:,iirgraupelm), zt )
        end if


        if ( iiNim > 0 ) then
          call stat_update_var( iNim, hydromet(:,iiNim), zt )
        end if

        if ( iiNrm > 0 ) then
          call stat_update_var( iNrm, hydromet(:,iiNrm), zt )
        end if

        if ( iiNsnowm > 0 ) then
          call stat_update_var( iNsnowm, hydromet(:,iiNsnowm), zt )
        end if

        if ( iiNgraupelm > 0 ) then
          call stat_update_var( iNgraupelm, hydromet(:,iiNgraupelm), zt )
        end if

        call stat_update_var( ithlm_mc, thlm_mc(:), zt )

        call stat_update_var( irtm_mc, rtm_mc(:), zt )

      end if

      ! Add microphysical tendencies to large-scale and radiative tendencies
      rtm_forcing  = rtm_forcing + rtm_mc
      thlm_forcing = thlm_forcing + thlm_mc


      if ( l_stats_samp .and. iirrainm > 0 ) then
        ! Rainfall rate (mm/day) should be defined on thermodynamic
        ! levels.  -Brian
        ! The absolute value of Vrr is taken because rainfall rate
        ! is a scalar quantity, and is therefore positive.
        call stat_update_var( irain_rate,  & 
           ( hydromet(:,iirrainm)  & 
             * zm2zt( abs( hydromet_vel(:,iirrainm) ) ) ) & 
              * ( rho / rho_lw ) & 
            * real( sec_per_day * 1000.0 ), zt )

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
        call stat_update_var_pt( irain, 1,  & 
            ( hydromet(2,iirrainm) & 
              * abs( zm2zt( hydromet_vel(:,iirrainm), 2 ) ) ) & 
                * ( rho(2) / rho_lw ) & 
            * real( sec_per_day * 1000.0 ), sfc )

        call stat_update_var_pt( irain_flux, 1, & 
           ( zt2zm( hydromet(:,iirrainm), 1 )  & 
             * abs( hydromet_vel(1,iirrainm) ) ) * ( rho_zm(1) / rho_lw )  & 
             * rho_lw * Lv, sfc )

        ! Also store the value of surface rain water mixing ratio.
        call stat_update_var_pt( irrainm_sfc, 1,  & 
               ( zt2zm( hydromet(:,iirrainm), 1 ) ), sfc )

      end if ! l_stats_samp


!       Error Report
!       Joshua Fasching Feb 2008

      if ( lapack_error( err_code ) .and.  &
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
        write(fstderr,*) "pdf_params%a = ", pdf_params%a
        write(fstderr,*) "pdf_params%rc1 = ", pdf_params%rc1
        write(fstderr,*) "pdf_params%rc2 = ", pdf_params%rc2
        write(fstderr,*) "pdf_params%s1 = ", pdf_params%s1
        write(fstderr,*) "pdf_params%s2 = ", pdf_params%s2
        write(fstderr,*) "pdf_params%ss1 = ", pdf_params%ss1
        write(fstderr,*) "pdf_params%ss2 = ", pdf_params%ss2

        write(fstderr,*) "Intent(inout)"

        write(fstderr,*) "Ncnm = ", Ncnm
        write(fstderr,*) "hydromet = ", hydromet
        write(fstderr,*) "Intent(out)"
        write(fstderr,*) "rtm_mc = ", rtm_mc
        write(fstderr,*) "thlm_mc = ", thlm_mc

      end if

      return
    end subroutine advance_microphys

!===============================================================================
    subroutine microphys_solve( solve_type, l_sed, dt, lhs, & 
                                xrm_tndcy, xrm, err_code )

! Description:

! References:
!-----------------------------------------------------------------------

      use grid_class, only: & 
          gr ! Variable(s)

      use stats_precision, only:  & 
          time_precision ! Variable(s)

      use lapack_wrap, only:  & 
          tridag_solve !,& ! Procedure(s)
!         band_solve


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
      real, intent(in), dimension(gr%nnzp) :: & 
        xrm_tndcy !                                     [units/s]

      ! Input/Output Variables
      real, intent(inout), dimension(3,gr%nnzp) :: & 
        lhs ! Left hand side

      real, intent(inout), dimension(gr%nnzp) :: & 
        xrm ! Hydrometeor being solved for              [units vary]

      ! Output Variables
      integer, intent(out) :: err_code

      ! Local Variables
      real, dimension(gr%nnzp) :: & 
        rhs ! Right hand side

      integer :: k, kp1, km1 ! Array indices


      integer :: & 
        ixrm_ma,  & ! Mean advection budget stats toggle
        ixrm_sd,  & ! Sedimentation budget stats toggle
        ixrm_dff    ! Diffusion budget stats toggle


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
      rhs(2:gr%nnzp-1)  & 
      = real((xrm(2:gr%nnzp-1) / dt )  & ! Time tendency
      + xrm_tndcy(2:gr%nnzp-1))


      ! Boundary condition on the RHS
      rhs(1) = real( xrm(1) / dt )
      rhs(gr%nnzp) =  & 
         real( ( xrm(gr%nnzp) / dt ) + xrm_tndcy(gr%nnzp-1) )


      ! Solve system using tridag_solve. This uses LAPACK sgtsv,
      ! which relies on Gaussian elimination to decompose the matrix.
      call tridag_solve & 
           ( solve_type, gr%nnzp, 1, lhs(1,:), lhs(2,:), lhs(3,:), & 
             rhs, xrm, err_code )

      if ( l_stats_samp ) then

        do k = 1, gr%nnzp, 1

          km1 = max( k-1, 1 )
          kp1 = min( k+1, gr%nnzp )

          ! Finalize implicit contributions

          ! xrm term ma is completely implicit; call stat_update_var_pt.
          call stat_update_var_pt( ixrm_ma, k, & 
                ztscr01(k) * xrm(km1) & 
              + ztscr02(k) * xrm(k) & 
              + ztscr03(k) * xrm(kp1), zt)

          ! xrm term sd is completely implicit; call stat_update_var_pt.
          if ( l_sed ) then
            call stat_update_var_pt( ixrm_sd, k, & 
                  ztscr04(k) * xrm(km1) & 
                + ztscr05(k) * xrm(k) & 
                + ztscr06(k) * xrm(kp1), zt )
          end if

          ! xrm term dff is completely implicit; call stat_update_var_pt.
          call stat_update_var_pt( ixrm_dff, k, & 
                ztscr07(k) * xrm(km1) & 
              + ztscr08(k) * xrm(k) & 
              + ztscr09(k) * xrm(kp1), zt )

        enddo ! 1..gr%nnzp

      endif ! l_stats_samp


! Boundary conditions on results
!xrm(1) = xrm(2)
! Michael Falk, 7 Sep 2007, made this change to eliminate problems
! with anomalous rain formation at the top boundary.
!        xrm(gr%nnzp) = 0
!xrm(gr%nnzp) = xrm(gr%nnzp-1)
! eMFc

      return
    end subroutine microphys_solve

!===============================================================================
    subroutine microphys_lhs & 
               ( solve_type, l_sed, dt, Kr, nu, wm_zt, V_hm, lhs )

! Description:
! Setup the matrix of implicit contributions to a term.
! Includes the effects of sedimentation, diffusion, and advection.
!
! Notes:
! Setup for tridiagonal system and boundary conditions should be the same as
! the original rain subroutine code.
!-----------------------------------------------------------------------

      use grid_class, only:  & 
          gr,  & ! Variable(s)
          zm2zt ! Procedure(s)

      use stats_precision, only:  & 
          time_precision ! Variable(s)

      use diffusion, only:  & 
          diffusion_zt_lhs ! Procedure(s)

      use mean_adv, only:  & 
          term_ma_zt_lhs ! Procedure(s)

      use constants, only: sec_per_day ! Variable(s)

      use stats_variables, only: & 
          irrainm_ma,   & ! Variable(s)
          irrainm_sd, & 
          irrainm_dff, & 
          iNrm_ma, & 
          iNrm_sd, & 
          iNrm_dff, & 
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

      ! Input Variables
      character(len=*), intent(in) :: &
        solve_type  ! Description of which hydrometeor is being solved for.

      logical, intent(in) ::  & 
        l_sed    ! Whether to add a hydrometeor sedimentation term.

      real(kind=time_precision), intent(in) ::  & 
        dt       ! Model timestep                                           [s]

      real, intent(in) ::  & 
        nu       ! Background diffusion coefficient                         [m^2/s]

      real, intent(in), dimension(gr%nnzp) ::  & 
        wm_zt, & ! w wind component on thermodynamic levels                 [m/s]
        V_hm,  & ! Sedimentation velocity of hydrometeor (momentum levels)  [m/s]
        Kr       ! Eddy diffusivity for hydrometeor on momentum levels      [m^2/s]

      real, intent(out), dimension(3,gr%nnzp) :: & 
        lhs      ! Left hand side of tridiagonal matrix.

      ! Local Variables
      real, dimension(3) :: tmp

      ! Array indices
      integer :: k, km1

      !integer kp1


      integer :: & 
        ixrm_ma,  & ! Mean advection budget stats toggle
        ixrm_sd,  & ! Sedimentation budget stats toggle
        ixrm_dff    ! Diffusion budget stats toggle

      ! ----- Begin Code -----

      select case( solve_type )
      case( "rrainm" )
        ixrm_ma  = irrainm_ma
        ixrm_sd  = irrainm_sd
        ixrm_dff = irrainm_dff
      case( "Nrm" )
        ixrm_ma  = iNrm_ma
        ixrm_sd  = iNrm_sd
        ixrm_dff = iNrm_dff
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
      case default
        ixrm_ma  = 0
        ixrm_sd  = 0
        ixrm_dff = 0
      end select


      ! Reset LHS Matrix for current timestep.
      lhs = 0.0

      ! Setup LHS Matrix
      do k = 2, gr%nnzp-1, 1

        km1 = max( k-1, 1 )
        !kp1 = min( k+1, gr%nnzp )

        ! Main diagonal

        ! LHS time tendency.
        lhs(k_tdiag,k) = real( lhs(k_tdiag,k) + ( 1.0 / dt ) )


        ! All diagonals

        ! LHS eddy-diffusion term.
        lhs(kp1_tdiag:km1_tdiag,k) & 
        = lhs(kp1_tdiag:km1_tdiag,k) & 
        + diffusion_zt_lhs( Kr(k), Kr(km1), nu,  & 
                            gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )

        ! LHS mean advection term.
        lhs(kp1_tdiag:km1_tdiag,k) & 
        = lhs(kp1_tdiag:km1_tdiag,k) & 
        + term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )

        ! LHS hydrometeor sedimentation term.
        ! Note: originally pristine ice did not sediment, so it was
        ! setup to be disabled as needed, but now pristine ice does
        ! sediment in the COAMPS case. -dschanen 12 Feb 2007
        if ( l_sed ) then
          lhs(kp1_tdiag:km1_tdiag,k) & 
          = lhs(kp1_tdiag:km1_tdiag,k) & 
          + sedimentation( V_hm(k), V_hm(km1), gr%dzt(k), k )
        endif

        if ( l_stats_samp ) then

          ! Statistics:  implicit contributions to hydrometeor xrm.

          if ( ixrm_ma > 0 ) then
            tmp(1:3) = term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )
            ztscr01(k) = -tmp(3)
            ztscr02(k) = -tmp(2)
            ztscr03(k) = -tmp(1)
          endif

          if ( ixrm_sd > 0 ) then
            tmp(1:3) = sedimentation( V_hm(k), V_hm(km1), gr%dzt(k), k )
            ztscr04(k) = -tmp(3)
            ztscr05(k) = -tmp(2)
            ztscr06(k) = -tmp(1)
          endif

          if ( ixrm_dff > 0 .and. l_sed ) then
            tmp(1:3) & 
            = diffusion_zt_lhs( Kr(k), Kr(km1), nu,  & 
                                gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
            ztscr07(k) = -tmp(3)
            ztscr08(k) = -tmp(2)
            ztscr09(k) = -tmp(1)
          endif

        endif ! l_stats_samp

      enddo ! 2..gr%nnzp-1


      ! Boundary Conditions

      ! The hydrometeor eddy-diffusion term has zero-flux boundary conditions, meaning
      ! that amounts of a hydrometeor are not allowed to escape the model boundaries
      ! through the process of eddy-diffusion.  It should be noted that amounts of a
      ! hydrometeor are allowed to leave the model at the lower boundary through the
      ! process of hydrometeor sedimentation.  However, only the eddy-diffusion term
      ! contributes to the LHS matrix at the k=1 and k=gr%nnzp levels.  Thus, function
      ! diffusion_zt_lhs needs to be called at both the upper boundary level and the
      ! lower boundary level.


      ! Lower Boundary
      k   = 1
      km1 = max( k-1, 1 )
      ! Note:  In function diffusion_zt_lhs, at the k=1 (lower boundary) level,
      !        variables referenced at the km1 level don't factor into the equation.

      ! LHS time tendency at the lower boundary.
      lhs(k_tdiag,k) = real( lhs(k_tdiag,k) + ( 1.0 / dt ) )

      ! LHS eddy-diffusion term at the lower boundary.
      lhs(kp1_tdiag:km1_tdiag,k) &
      = lhs(kp1_tdiag:km1_tdiag,k) &
      + diffusion_zt_lhs( Kr(k), Kr(km1), nu,  &
                          gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )

      if ( l_stats_samp ) then

        ! Statistics:  implicit contributions to hydrometeor xrm.

        if ( ixrm_dff > 0 ) then
          tmp(1:3) & 
          = diffusion_zt_lhs( Kr(k), Kr(km1), nu,  & 
                              gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
          ztscr07(k) = -tmp(3)
          ztscr08(k) = -tmp(2)
          ztscr09(k) = -tmp(1)
        endif

      endif  ! l_stats_samp


      ! Upper Boundary
      k   = gr%nnzp
      km1 = max( k-1, 1 )

      ! LHS time tendency at the upper boundary.
      lhs(k_tdiag,k) = real( lhs(k_tdiag,k) + ( 1.0 / dt ) )

      ! LHS eddy-diffusion term at the upper boundary.
      lhs(kp1_tdiag:km1_tdiag,k) &
      = lhs(kp1_tdiag:km1_tdiag,k) &
      + diffusion_zt_lhs( Kr(k), Kr(km1), nu,  &
                          gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )

      if ( l_stats_samp ) then

        ! Statistics:  implicit contributions to hydrometeor xrm.

        if ( ixrm_dff > 0 ) then
          tmp(1:3) & 
          = diffusion_zt_lhs( Kr(k), Kr(km1), nu,  & 
                              gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )
          ztscr07(k) = -tmp(3)
          ztscr08(k) = -tmp(2)
          ztscr09(k) = -tmp(1)
        endif

      endif  ! l_stats_samp


      return
    end subroutine microphys_lhs

!===============================================================================
    pure function sedimentation( V_hm, V_hmm1, dzt, level ) & 
      result( lhs )

! Description:
! Sedimentation of a hydrometeor:  implicit portion of the code.
!
! The variable "hm" stands for one of the five hydrometeor variables currently
! in the code:  mean rain mixing ratio (rrainm), mean rain drop
! concentration (Nrm), mean ice mixing ratio (ricem), mean snow mixing
! ratio (rsnowm), or mean graupel mixing ratio (rgraupelm).  The variable
! "V_hm" stands for the sedimentation velocity of the appropriate hydrometeor.
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
! reversed and the leading "-" in front of the term is changed to a "+".
!
! Timestep index (t) stands for the index of the current timestep, while
! timestep index (t+1) stands for the index of the next timestep, which is being
! advanced to in solving the d(hm)/dt equation.
!
! This term is discretized as follows:
!
! The values of hm are found on the thermodynamic levels, while the values of
! V_hm are found on the momentum levels.  The variable hm is interpolated to the
! intermediate momentum levels.  At the intermediate momentum levels, the
! interpolated values of hm are multiplied by the values of V_hm.  Then, the
! derivative of (hm*V_hm) is taken over the central thermodynamic level.
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
! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond with
! altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.  The
! letter "t" is used for thermodynamic levels and the letter "m" is used for
! momentum levels.
!
! dzt(k) = 1 / ( zm(k) - zm(k-1) )
!
!
! Conservation Properties:
!
! When a hydrometeor is sedimented to the ground (or out the lower boundary of
! the model), it is removed from the atmosphere (or from the model domain).
! Thus, the quantity of the hydrometeor over the entire vertical domain should
! not be conserved due to the process of sedimentation.  Thus, not all of the
! column totals in the left-hand side matrix should be equal to 0.  Instead, the
! sum of all the column totals should equal the flux of hm out the bottom
! (zm(1) level) of the domain, -V_hm(1) * ( D(2)*hm(1) + C(2)*hm(2) ), where the
! factor in parentheses is the interpolated value of hm at the zm(1) level.
! Furthermore, most of the individual column totals should sum to 0, but the 1st
! and 2nd (from the left) columns should combine to sum to the flux out the
! bottom of the domain.
!
! To see that this modified conservation law is satisfied, compute the
! sedimentation of hm and integrate vertically.  In discretized matrix notation
! (where "i" stands for the matrix column and "j" stands for the matrix row):
!
!  -V_hm(1) * ( D(2)*hm(1) + C(2)*hm(2) )
!     =
!  Sum_j Sum_i ( 1/dzt )_i ( d (V_hm * weights_hm) / dz )_ij hm_j.
!
! The left-hand side matrix, ( d (V_hm * weights_hm) / dz )_ij, is partially
! written below.  The sum over i in the above equation removes dzt everywhere
! from the matrix below.  The sum over j leaves the column totals and the flux
! at zm(1) that are desired.
!
! Left-hand side matrix contributions from the sedimentation term (only); first
! four vertical levels:
!
!     -------------------------------------------------------------------------------->
!k=1 |            0                          0                            0
!    |
!k=2 | -dzt(k)*V_hm(k-1)*D(k)   +dzt(k)*[ V_hm(k)*B(k)       +dzt(k)*V_hm(k)*A(k)
!    |                                   -V_hm(k-1)*C(k) ]
!    |
!k=3 |            0             -dzt(k)*V_hm(k-1)*D(k)       +dzt(k)*[ V_hm(k)*B(k)
!    |                                                                -V_hm(k-1)*C(k) ]
!    |
!k=4 |            0                          0               -dzt(k)*V_hm(k-1)*D(k)
!    |
!   \ /
!
! The variables A(k), B(k), C(k), and D(k) are weights of interpolation around
! the central thermodynamic level (k), such that:
! A(k) = ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ),
! B(k) = 1 - [ ( zm(k) - zt(k) ) / ( zt(k+1) - zt(k) ) ]
!      = 1 - A(k);
! C(k) = ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) ), and
! D(k) = 1 - [ ( zm(k-1) - zt(k-1) ) / ( zt(k) - zt(k-1) ) ]
!      = 1 - C(k).
!
! Note:  The superdiagonal term from level 3 and both the main diagonal and
!        superdiagonal terms from level 4 are not shown on this diagram.

! References:
! None

! Notes:
! Both COAMPS Microphysics and Brian Griffin's implementation use Khairoutdinov
! and Kogan (2000) for the calculation of rain mixing ratio and rain droplet
! number concentration sedimentation velocities.
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
        V_hm,    & ! Sedimentation velocity of hydrometeor (k)                [m/s]
        V_hmm1,  & ! Sedimentation velocity of hydrometeor (k-1)              [m/s]
        dzt        ! Inverse of grid spacing (k)                              [m]

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


      elseif ( level > 1 .and. level < gr%nnzp ) then

        ! Most of the interior model; normal conditions.
        
        ! Vince Larson pulled V_hm inside derivative to make conservative.
        ! 13 Dec 2007
        !
        ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
!       lhs(kp1_tdiag)  &
!       = + V_hmzt * dzt * gr%weights_zt2zm(t_above,mk)
      
!       ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
!       lhs(k_tdiag)  &
!       = + V_hmzt * dzt * (   gr%weights_zt2zm(t_below,mk)  &
!                          - gr%weights_zt2zm(t_above,mkm1)  )
!
!       ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
!       lhs(km1_tdiag)  &
!       = - V_hmzt * dzt * gr%weights_zt2zm(t_below,mkm1)

        ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
        lhs(kp1_tdiag)  & 
        = + dzt * V_hm * gr%weights_zt2zm(t_above,mk)

        ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
        lhs(k_tdiag)  & 
        = + dzt * (   V_hm * gr%weights_zt2zm(t_below,mk) & 
                    - V_hmm1 * gr%weights_zt2zm(t_above,mkm1)  )

        ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
        lhs(km1_tdiag)  & 
        = - dzt * V_hmm1 * gr%weights_zt2zm(t_below,mkm1)

        !  End Vince Larson change


      elseif ( level == gr%nnzp ) then

        ! k = gr%nnzp (top level); upper boundary level; no flux.

        ! Thermodynamic superdiagonal: [ x hm(k+1,<t+1>) ]
        lhs(kp1_tdiag) = 0.0

        ! Thermodynamic main diagonal: [ x hm(k,<t+1>) ]
        lhs(k_tdiag)   = 0.0

        ! Thermodynamic subdiagonal: [ x hm(k-1,<t+1>) ]
        lhs(km1_tdiag) = 0.0


      endif


      return
    end function sedimentation

!===============================================================================
    subroutine adj_microphys_tndcy( xrm_tndcy, wm_zt, V_hm, Kr, nu, & 
                                    dt, level, l_sed, & 
                                    xrm, overevap_rate )

! DESCRIPTION:  Correction for the over-evaporation of a hydrometeor.
!
! If a small amount of a hydrometeor (such as rain water) gets diffused into an
! area that is very dry (such as right above the cloud top), the hydrometeor
! (rain water) will have a very high rate of evaporation and will evaporate
! entirely in a short amount of time.  However, the evaporation rate is computed
! instantaneously at a given moment in time.  This rate is then projected over
! the entire length of the given timestep.  Therefore, a high-enough rate of
! evaporation combined with a small-enough amount of the hydrometeor (rain
! water) and a long-enough timestep will cause the hydrometeor value (rain water
! mixing ratio) to be negative by the end of the timestep.  Therefore, a
! correction factor needs to be imposed on the evaporation rate so that the
! amount of the hydrometeor (rain water mixing ratio) does not fall below 0.
!
! Besides over-evaporation of a hydrometeor, other factors may come into play
! that cause the value of a hydrometeor to fall below 0.  These factors are due
! to the nature of implicit discretization and numerical errors.  In a nutshell,
! the eddy diffusion parameter used currently in this model smooths out the
! entire hydrometeor profile as a whole at every timestep.  This smoothing may
! cause negative values at certain levels.  Also, mean advection and hydrometeor
! sedimentation can cause negative values to occur in the hydrometeor.  This can
! happen in places where the profile abruptly goes from a large positive value
! to 0 (such as at cloud top).  The nature of the discretization of taking a
! derivative at these levels may cause negative values of a hydrometeor.
!
! This subroutine is called only if a hydrometeor at a certain level contains
! a negative value.  First, this subroutine uses the same methods that the model
! statistical code uses in computing budget terms in order to determine what
! factors effected the value of the given hydrometeor during the timestep that
! was just solved for.  The mean advection, sedimentation, and diffusion budget
! terms are all computed.  These three terms are then added together to make up
! the total transport and sedimentation tendency.  This tendency is then added
! to the total microphysical tendency to find the overall hydrometeor tendency.
! The overall hydrometeor tendency is then multiplied by the timestep length to
! find the net change in the hydrometeor over the last timestep.  This net
! change is then added to the current value of the hydrometeor in order to find
! the value of the hydrometeor at the previous timestep.  This method has been
! well tested and produces accurate results.
!
! Once the value of the hydrometeor at the previous timestep has been found, the
! net change in the hydrometeor due to ONLY mean advection, diffusion, and
! sedimentation is calculated.  This net change is added to the value of the
! hydrometeor at the previous timestep.  If the new value is below zero, then
! the negative value of the hydrometeor was caused by the mean advection,
! diffusion, and sedimentation terms.  The microphysical terms (evaporation)
! did not cause the negative value.  There was no over-evaporation and the
! evaporation rate can be set to 0.  However, if the new value of the
! hydrometeor is greater than or equal to 0, then the microphysical tendencies
! (evaporation) did cause the hydrometeor array to have negative values.  The
! amount of hydrometeor evaporated is set equal to the amount that was left-over
! after the transport and sedimentation effects were added in.  The evaporation
! rate is that amount divided by the timestep.  This can be viewed as the
! timestep-average evaporation rate, whereas the rate previously calculated can
! be viewed as the instantaneous evaporation rate.  The amount of the
! hydrometeor that was over-evaporated is the amount of the hydrometeor that is
! negative.  The over-evaporation rate is that amount divided by the timestep.
!
! It should be noted that this is important because the rain water mixing ratio
! time tendency (drr/dt) due to microphysics at every level is incorporated into
! the total water mixing ratio (rtm) and liquid water potential temperature
! (thlm) equations.  Any artificial excess in evaporation will artificially
! increase water vapor, and thus rtm, and artificially decrease thlm (due to
! evaporative cooling).  This may result in an artificial increase in cloud
! water.
!
! rrainm_mc_tndcy = rrainm_cond + rrainm_auto + rrainm_accr
! rtm_mc  = - rrainm_mc_tndcy
! thlm_mc = ( Lv / (Cp*exner) ) * rrainm_mc_tndcy
!
! Anyplace where rrainm drops below zero due to microphysics, there is too much
! evaporation rate for the timestep, so rrainm_cond is too negative.  We must
! add in the over-evaporated amount of rrainm/dt to make the rate accurate.  The
! over-evaporated amount is being defined as a positive scalar, so that:
! overevap_rrainm = -rrainm (where rrainm < 0) -- this makes overevap_rrainm
! positive.
!
! New cond/evap rate = rrainm_cond + overevap_rrainm/dt
! (overevap_rate = overevap_rrainm/dt)
! -- since rrainm_cond can only be negative (we don't allow rain droplets to
!    grow by condensation) and overevap_rrainm/dt can only be positive (we
!    define it that way), the new cond/evap rate will be less negative, which
!    is what we want.
!
! To update the effects of microphysics on rtm and thl:
!
! rtm_mc = rtm_mc - overevap_rate
! thlm_mc = thlm_mc + ( Lv / (Cp*exner) ) * overevap_rate
!
! This is done in the subroutine which calls this one.
!
! If the hydrometeor is negative due to reasons besides over-evaporation, the
! value is clipped.  This is statistically stored in the clipping array.  This
! is also done in the subroutine which calls this one.
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

      implicit none

      ! Input variables.

      real, dimension(gr%nnzp), intent(in) :: &
        xrm_tndcy, & ! Hydrometeor microphysical tendency.                      [hm_units/s]
        wm_zt,     & ! Vertical velocity (thermo. levels).                      [m/s]
        V_hm,      & ! Sedimentation velocity (interpolated to thermo. levels). [m/s]
        Kr           ! Eddy diffusivity for hydrometeors (m-lev).               [m^2/s]

      real, intent(in) :: nu  ! Diffusion coefficient      [m^2/s]

      real(kind=time_precision), intent(in) :: dt  ! Timestep   [s]

      integer, intent(in) :: level  ! Vertical grid index

      logical, intent(in) :: l_sed   ! Whether to add a sedimentation term


      ! Input/output variable.

      real, dimension(gr%nnzp), intent(inout) :: &
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

!     real :: evap_amt          ! The actual evaporation amount over the t.s.
!     real :: evap_rate         ! The time-averaged rate.
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
      kp1 = min( k+1, gr%nnzp )


! Mean advection tendency component

! The implicit (LHS) value of the mean advection component of the equation used
! during the timestep that was just solved for.
      tmp(1:3) = term_ma_zt_lhs( wm_zt(k), gr%dzt(k), k )

      ma_subdiag  = -tmp(3) ! subdiagonal
      ma_maindiag = -tmp(2) ! main diagonal
      ma_supdiag  = -tmp(1) ! superdiagonal

      ma_tndcy = & 
         + ma_subdiag  * xrm(km1) & 
         + ma_maindiag * xrm(k) & 
         + ma_supdiag  * xrm(kp1)


! Sedimentation tendency component

      if ( l_sed ) then

        ! The implicit (LHS) value of the sedimentation component of the equation
        ! used during the timestep that was just solved for.
        tmp(1:3) = sedimentation( V_hm(k), V_hm(km1),  & 
                                  gr%dzt(k), k )

        sd_subdiag  = -tmp(3) ! subdiagonal
        sd_maindiag = -tmp(2) ! main diagonal
        sd_supdiag  = -tmp(1) ! superdiagonal

        sd_tndcy =  & 
           + sd_subdiag  * xrm(km1) & 
           + sd_maindiag * xrm(k) & 
           + sd_supdiag  * xrm(kp1)

      else

        sd_tndcy = 0.0

      endif


! Diffusion tendency component

! The implicit (LHS) value of the diffusion component of the equation used
! during the timestep that was just solved for.
      tmp(1:3) & 
         = diffusion_zt_lhs( Kr(k), Kr(km1), nu,  & 
                             gr%dzm(km1), gr%dzm(k), gr%dzt(k), k )

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
      xrm_chge = real( tot_tndcy * dt )

! The value of xrm at the previous timestep.
      xrm_old = xrm(k) - xrm_chge

! The net amount of change in the hydrometeor due to only the transport (mean
! advection and diffusion) and sedimentation terms.
      xrm_chge_trsed = real( trnsprt_sed_tndcy * dt )

! The new value of the hydrometeor at this timestep due to only
! the transport and sedimentation terms.
      xrm_trsed_only = xrm_old + xrm_chge_trsed

      if ( xrm_trsed_only >= 0.0 ) then
        ! The negative value of hydrometeor (xrm) is due ONLY to microphysical
        ! tendencies, namely the over-evaporation of xrm.
        ! Find the actual amount of the hydrometeor that evaporated during the
        ! timestep to make the value of xrm go to 0.
!    evap_amt = -xrm_trsed_only
        ! Divide by the timestep to find the actual evaporation rate.
!    evap_rate = evap_amt / dt
        ! The amount of the hydrometeor that was artificially excessively evaporated.
        ! Define as positive.
        overevap_amt = -xrm(k)
        ! Divide by the timestep to find the over-evaporation rate.  Define as
        ! positive.  This rate should also be the difference between the computed
        ! evaporation rate (xrm_tndcy) and the actual evaporation rate (evap_rate).
        overevap_rate = real( overevap_amt / dt )
        ! Reset the value of the hydrometeor (xrm) to 0.
        xrm(k) = 0.0
      else
        ! The negative value of hydrometeor (xrm) is due to transport (mean advection
        ! and diffusion) and sedimentation.  Find the actual amount of the
        ! hydrometeor that evaporated during the timestep to make the value of xrm go
        ! to 0.  Even though the microphysical tendency portion of the code may have
        ! computed an evaporation rate, we figure that the transport and
        ! sedimentation terms made the value of the hydrometeor negative, so we say
        ! that the evaporation amount and rate is 0.
!       evap_amt = 0.0
!       evap_rate = 0.0
        ! The amount of the hydrometeor that was artificially excessively evaporated.
        ! Define as positive.  In this case, any evaporation that was computed is
        ! considered to be over-evaporation.  Define as positive.
        overevap_amt = real( -xrm_tndcy(k) * dt )
        overevap_rate = -xrm_tndcy(k)
        ! Currently reset xrm to xrm_trsed_only.  This is done to make the
        ! statistical budget for xrm balance correctly.  The value of xrm(k) will
        ! still be negative at this this point.  However, it will be less negative
        ! because it has been adjusted for over-evaporation.  The remaining negative
        ! value of hydrometeor xrm, which is due to transport and sedimentation, will
        ! be zeroed out in clipping in the subroutine that calls this one.
        xrm(k) = xrm_trsed_only
      endif

      return
    end subroutine adj_microphys_tndcy

!-------------------------------------------------------------------------------
    subroutine return_LH_index( hydromet_index, LH_count, LH_index )
! Description:
!   Set the Latin hypercube variable index if the hydrometeor exists
! References:
!   None
!-------------------------------------------------------------------------------

      implicit none

      ! Input Variables
      integer, intent(in) :: &
        hydromet_index

      ! Input/Output Variables
      integer, intent(inout) :: &
        LH_count

      ! Output Variables
      integer, intent(out) :: &
        LH_index

      ! ---- Begin Code ----

      if ( hydromet_index > 0 ) then
        LH_count = LH_count + 1
        LH_index = LH_count
      else
        LH_index = -1
      end if

      return
    end subroutine return_LH_index
!===============================================================================

  end module microphys_driver
