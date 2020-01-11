! $Id$
module morrison_microphys_module

  implicit none

  public :: morrison_microphys_driver

  private

  contains
!-------------------------------------------------------------------------------
  subroutine morrison_microphys_driver &
             ( dt, nz, &
               l_latin_hypercube, thlm, wm_zt, p_in_Pa, &
               exner, rho, cloud_frac, w_std_dev, &
               dzq, rcm, Ncm, chi, rvm, hydromet, &
               hydromet_mc, hydromet_vel_zt, Ncm_mc, &
               rcm_mc, rvm_mc, thlm_mc, &
               microphys_stats_zt, microphys_stats_sfc )

    ! Description:
    ! Wrapper for the Morrison microphysics
    !
    ! References:
    ! None
    !-----------------------------------------------------------------------

    use parameters_model, only: hydromet_dim

    ! The version of the Morrison 2005 microphysics that is in SAM.
    use module_MP_graupel, only: &
        M2005MICRO_GRAUPEL  ! Procedure

    use grid_class, only: &
        gr

    use constants_clubb, only: &
        Lv,   & ! Constants
        Ls,   &
        Cp,   &
        grav, &
        zero

    use stats_variables, only: & 
        irsm_sd_morr, & ! Variables
        irim_sd_mg_morr, & 
        irrm_sd_morr, & 
        irsm_sd_morr, &
        irgm_sd_morr, &
        ircm_sd_mg_morr, &
        irrm_auto, &
        irrm_accr, &
        irrm_evap, &
        irsm_sd_morr_int, &
        ihl_on_Cp_residual, &
        iqto_residual

    use stats_variables, only: & 
        ieff_rad_cloud, & ! Variables
        ieff_rad_ice, &
        ieff_rad_snow, &
        ieff_rad_rain, &
        ieff_rad_graupel

    use stats_variables, only: & 
        iprecip_rate_sfc, & !Variables
        imorr_snow_rate, &
        iNrm_auto, &
        iNrm_evap

    use stats_variables, only: &
        iPSMLT, & ! Variable(s)
        iEVPMS, &
        iPRACS, &
        iEVPMG, &
        iPRACG, &
        iPGMLT, &
        iMNUCCC, &
        iPSACWS, &
        iPSACWI, &
        iQMULTS, &
        iQMULTG, &
        iPSACWG, &
        iPGSACW, &
        iPRD, &
        iPRCI, &
        iPRAI, &
        iQMULTR, &
        iQMULTRG, &
        iMNUCCD, &
        iPRACI, &
        iPRACIS, &
        iEPRD, &
        iMNUCCR, &
        iPIACR, &
        iPIACRS, &
        iPGRACS, &
        iPRDS, &
        iEPRDS, &
        iPSACR, &
        iPRDG, &
        iEPRDG

    use stats_variables, only: &
        iNGSTEN, & ! Lots of variables
        iNRSTEN, &
        iNISTEN, &
        iNSSTEN, &
        iNCSTEN, &
        iNPRC1,  &
        iNRAGG,  &
        iNPRACG, &
        iNSUBR,  &
        iNSMLTR, &
        iNGMLTR, &
        iNPRACS, &
        iNNUCCR, &
        iNIACR,  &
        iNIACRS, &
        iNGRACS, &
        iNSMLTS, &
        iNSAGG, &
        iNPRCI, &
        iNSCNG, &
        iNSUBS, &
        iPRA, &
        iPRC, &
        iPRE

    use stats_variables, only: &
        iPCC, & ! Even more variables
        iNNUCCC, & 
        iNPSACWS, &
        iNPRA, &
        iNPRC, &
        iNPSACWI, &
        iNPSACWG, &
        iNPRAI, &
        iNMULTS, & 
        iNMULTG, &
        iNMULTR, &
        iNMULTRG, &
        iNNUCCD, &
        iNSUBI, &
        iNGMLTG, &
        iNSUBG, &
        iNACT, &
        iSIZEFIX_NR, &
        iSIZEFIX_NC, &
        iSIZEFIX_NI, &
        iSIZEFIX_NS, &
        iSIZEFIX_NG, &
        iNEGFIX_NR, &
        iNEGFIX_NC, &
        iNEGFIX_NI, &
        iNEGFIX_NS, &
        iNEGFIX_NG, &
        iNIM_MORR_CL, &
        iQC_INST, & 
        iQR_INST, &
        iQI_INST, &
        iQS_INST, &
        iQG_INST, &
        iNC_INST, &
        iNR_INST, &
        iNI_INST, & 
        iNS_INST, &
        iNG_INST, &
        iT_in_K_mc

    use T_in_K_module, only: &
        T_in_K2thlm, & ! Procedure(s)
        thlm2T_in_K

    use model_flags, only: &
        l_evaporate_cold_rcm  ! Flag(s)

    use parameters_microphys, only: &
        l_ice_microphys, & ! Flag(s)
        l_graupel

    use array_index, only:  & 
        iirr, iirs, iiri, iirg, &
        iiNr, iiNs, iiNi, iiNg

    use constants_clubb, only: &
        sec_per_day

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use variables_prognostic_module, only: &
        rho_ds_zt

    use error_code, only: &
        clubb_at_least_debug_level   ! Procedure

    use fill_holes, only: &
        vertical_integral

    use microphys_stats_vars_module, only: &
        microphys_stats_vars_type, & ! Type
        microphys_stats_alloc, & ! Procedure
        microphys_put_var

    implicit none

    ! External
    intrinsic :: max, real, maxval

    ! Constant parameters
    integer, parameter :: &
      num_output_zt = 200, & ! Overestimate of the number of variables to be output to zt grid
      num_output_sfc = 10    ! Same, for sfc grid


    real( kind = core_rknd ), parameter :: &
      w_thresh = 0.1_core_rknd ! Minimum value w for latin hypercube [m/s]

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: dt ! Model timestep        [s]

    integer, intent(in) :: nz ! Points in the Vertical        [-]

    logical, intent(in) :: &
      l_latin_hypercube   ! Whether we're using latin hypercube sampling

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      thlm,       & ! Liquid potential temperature       [K]
      p_in_Pa,    & ! Pressure                           [Pa]
      exner,      & ! Exner function                     [-]
      rho,        & ! Density on thermo. grid            [kg/m^3]
      cloud_frac    ! Cloud fraction                     [-]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      wm_zt, &     ! Mean vertical velocity on the thermo grid     [m/s]
      w_std_dev, & ! Standard deviation of vertical vel. [m/s]
      dzq          ! Change in altitude                  [m]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm,          & ! Cloud water mixing ratio                  [kg/kg]
      Ncm,          & ! Grid mean value for cloud droplet conc.    [#/kg]
      chi,     & ! The variable 's' from Mellor              [kg/kg]
      rvm             ! Vapor water mixing ratio                  [kg/kg]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel_zt ! Hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      Ncm_mc    ! Cloud droplet concentration time tendency    [num/kg/s]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    type(microphys_stats_vars_type), intent(out) :: &
      microphys_stats_zt, & ! Variables output for statistical sampling (zt grid)
      microphys_stats_sfc   ! Variables output for statistical sampling (sfc grid)

    ! Local Variables

    real( kind = core_rknd ), dimension(nz) :: &
      rrm_auto,     &  ! Autoconversion rate                 [kg/kg/s]
      rrm_accr,     &  ! Accretion rate                      [kg/kg/s]
      rrm_evap         ! Rain evaporation rate               [kg/kg/s]

    real, dimension(nz) :: & 
      effc, effi, effg, effs, effr ! Effective droplet radii [μ]

    real, dimension(nz) :: & 
      T_in_K,           & ! Temperature                                      [K]
      T_in_K_mc,        & ! Temperature tendency                           [K/s]
      rcm_r4,           & ! Temporary array for cloud water mixing ratio [kg/kg]
      Ncm_r4,           & ! Temporary array for cloud number conc.        [#/kg]
      Ncm_mc_r4,        & ! Temporary array for cloud number conc.      [#/kg/s]
      rvm_r4,           & ! Temporary array for vapor water mixing ratio [kg/kg]
      rcm_sten,         & ! Mean rc sedimentation tendency             [kg/kg/s]
      rrm_sten,      & ! Mean rr sedimentation tendency             [kg/kg/s]
      rim_sten,       & ! Mean ri sedimentation tendency             [kg/kg/s]
      rsm_sten,      & ! Mean rs sedimentation tendency             [kg/kg/s]
      rgm_sten,   & ! Mean rg sedimentation tendency             [kg/kg/s]
      morr_rain_vel_r4, & ! Rain fall velocity from Morrison microphysics  [m/s]
      cloud_frac_in       ! Cloud frac used as input for the Morrison scheme [-]

    ! In the comments below, by "adds to" we mean that if the quantity is
    ! positive, it adds positively to the prognostic variable, but if the
    ! quantity is negative, it subtracts from the prognostic variable.
    real, dimension(nz) :: &
      PSMLT,  & ! Freezing of rain to form snow.
                !    Adds to rsm, subtracts from rrm [kg/kg/s]
      EVPMS,  & ! Evaporation of melted snow.
                !    Adds to rsm, subtracts from rvm [kg/kg/s]
      PRACS,  & ! Collection of rain by snow.
                !    Adds to rsm, subtracts from rrm [kg/kg/s]
      EVPMG,  & ! Evaporation of melted graupel.
                !    Adds to rgm, subtracts from rvm [kg/kg/s]
      PRACG,  & ! Negative of collection of rain by graupel.
                !    Adds to rrm, subtracts from rgm [kg/kg/s]
      PGMLT,  & ! Negative of melting of graupel.
                !    Adds to rgm, subtracts from rrm [kg/kg/s]
      MNUCCC, & ! Contact freezing of cloud droplets.
                !    Adds to rim, subtracts from rcm [kg/kg/s]
      PSACWS, & ! Collection of cloud water by snow.
                !    Adds to rsm, subtracts from rcm [kg/kg/s]
      PSACWI, & ! Collection of cloud water by cloud ice.
                !    Adds to rim, subtracts from rcm [kg/kg/s]
      QMULTS, & ! Splintering from cloud droplets accreted onto snow.
                !    Adds to rim, subtracts from rcm [kg/kg/s]
      QMULTG, & ! Splintering from droplets accreted onto graupel.
                !    Adds to rim, subtracts from rcm [kg/kg/s]
      PSACWG, & ! Collection of cloud water by graupel.
                !    Adds to rgm, subtracts from rcm [kg/kg/s]
      PGSACW, & ! Reclassification of rimed snow as graupel.
                !    Adds to rgm, subtracts from rcm [kg/kg/s]
      PRD,    & ! Depositional growth of cloud ice.
                !    Adds to rim, subtracts from rcm [kg/kg/s]
      PRCI,   & ! Autoconversion of cloud ice to snow.
                !    Adds to rsm, subtracts from rim [kg/kg/s]
      PRAI,   & ! Collection of cloud ice by snow.
                !    Adds to rsm, subtracts from rim [kg/kg/s]
      QMULTR, & ! Splintering from rain droplets accreted onto snow.
                !    Adds to rim, subtracts from rrm [kg/kg/s]
      QMULTRG,& ! Splintering from rain droplets accreted onto graupel.
                !    Adds to rim, subtracts from rrm [kg/kg/s]
      MNUCCD, & ! Freezing of aerosol.
                !    Adds to rim, subtracts from rvm [kg/kg/s]
      PRACI,  & ! Collection of cloud ice by rain.
                !    Adds to rgm, subtracts from rim [kg/kg/s]
      PRACIS, & ! Collection of cloud ice by rain.
                !    Adds to rsm, subtracts from rim [kg/kg/s]
      EPRD,   & ! Negative of sublimation of cloud ice.
                !    Adds to rim, subtracts from rvm [kg/kg/s]
      MNUCCR, & ! Contact freezing of rain droplets.
                !    Adds to rgm, subtracts from rrm [kg/kg/s]
      PIACR,  & ! Collection of cloud ice by rain.
                !    Adds to rgm, subtracts from rrm [kg/kg/s]
      PIACRS, & ! Collection of cloud ice by rain.
                !    Adds to rsm, subtracts from rrm [kg/kg/s]
      PGRACS, & ! Collection of rain by snow.
                !    Adds to rgm, subtracts from rrm [kg/kg/s]
      PRDS,   & ! Depositional growth of snow.
                !    Adds to rsm, subtracts from rvm [kg/kg/s]
      EPRDS,  & ! Negative of sublimation of snow.
                !    Adds to rsm, subtracts from rvm [kg/kg/s]
      PSACR,  & ! Collection of snow by rain.
                !    Adds to rgm, subtracts from rsm [kg/kg/s]
      PRDG,   & ! Depositional growth of graupel.
                !    Adds to rgm, subtracts from rvm [kg/kg/s]
      EPRDG     ! Negative of sublimation of graupel.
                !    Adds to rgm, subtracts from rvm [kg/kg/s]

    real, dimension(nz) :: &
      NGSTEN, & ! Graupel sedimentation tendency [#/kg/s]
      NRSTEN, & ! Rain sedimentation tendency [#/kg/s]
      NISTEN, & ! Cloud ice sedimentation tendency [#/kg/s]
      NSSTEN, & ! Snow sedimentation tendency [#/kg/s]
      NCSTEN, & ! Cloud water sedimentation tendency [#/kg/s]
      NPRC1,  & ! Change in Nrm due to autoconversion of droplets. Adds to Nrm [#/kg/s]
      NRAGG,  & ! Change in Nrm due to self-collection of raindrops. Adds to Nrm [#/kg/s]
      NPRACG, & ! Collection of rainwater by graupel. Subtracts from Nrm [#/kg/s]
      NSUBR,  & ! Loss of Nrm by evaporation. Adds to Nrm [#/kg/s]
      NSMLTR, & ! Melting of snow to form rain. Subtracts from Nrm [#/kg/s]
      NGMLTR, & ! Melting of graupel to form rain. Subtracts from Nrm [#/kg/s]
      NPRACS, & ! Collection of rainwater by snow. Subtracts from Nrm [#/kg/s]
      NNUCCR, & ! Contact freezing of rain. Adds to Ngm, subtracts from Nrm [#/kg/s]
      NIACR,  & ! Collection of cloud ice by rain.
                !    Adds to Ngm, subtracts from Nrm and Nim [#/kg/s] 
      NIACRS, & ! Collection of cloud ice by rain.
                !    Adds to Nsm, subtracts from Nrm and Nim [#/kg/s]
      NGRACS, & ! Collection of rain by snow.
                !    Adds to Ngm, subtracts from Nrm and Nsm [#/kg/s]
      NSMLTS, & ! Melting of snow
                !    Adds to Nsm [#/kg/s]
      NSAGG, &  ! Self collection of snow
                !    Adds to Nsm [#/kg/s]
      NPRCI, &  ! Autoconversion of cloud ice to snow
                !    Adds to Nsm, subtracts from Nim [#/kg/s]
      NSCNG, &  ! Conversion of snow to graupel
                !    Adds to Ngm, subtracts from Nsm [#/kg/s]
      NSUBS, &  ! Loss of Nsm due to sublimation
                !    Adds to Nsm [#/kg/s]
      PRA,   &  ! Accretion. Adds to rrm, subtracts from rcm [kg/kg/s]
      PRC,   &  ! Autoconversion. Adds to rrm, subtracts from rcm [kg/kg/s]
      PRE       ! Rain evaporation. Subtracts from rrm [kg/kg/s]              
          
    real, dimension(nz) :: &
      PCC, &    ! Saturation adjustment 
                !    Adds to rcm, substracts from rvm [kg/kg/s]
      NNUCCC, & ! Contact freezing of cloud drops
                !    Adds to Nim, subtracts from Ncm [#/kg/s]
      NPSACWS,& ! Droplet accretion by snow. Subtracts from Ncm [#/kg/s] 
      NPRA, &   ! Droplet accretion by rain. Subtracts from Ncm [#/kg/s]
      NPRC, &   ! Autoconversion of cloud drops. Subtracts from Ncm [#/kg/s]
      NPSACWI,& ! Droplet accretion by cloud ice. Subtracts from Ncm [#/kg/s]
      NPSACWG,& ! Collection of cloud drops by graupel. Subtracts from Ncm [#/kg/s]
      NPRAI, &  ! Accretion of cloud ice by snow. Subtracts from Nim [#/kg/s]
      NMULTS, & ! Ice mult. due to riming of cloud drops by snow. Adds to Nim [#/kg/s]
      NMULTG, & ! Ice mult. due to accretion of cloud drops by graupel. Adds to Nim [#/kg/s]
      NMULTR, & ! Ice mult. due to riming of rain by snow. Adds to Nim [#/kg/s]
      NMULTRG,& ! Ice mult. due to accretion of rain by graupel. Adds to Nim [#/kg/s]
      NNUCCD, & ! Primary ice nucleation, freezing of aerosol. Adds to Nim [#/kg/s]
      NSUBI, &  ! Loss of ice due to sublimation. Subtracts from Nim [#/kg/s]
      NGMLTG, & ! Loss of graupel due to melting. Subtracts from Ngm [#/kg/s]
      NSUBG, &  ! Loss of graupel due to sublimation. Subtracts from Ngm [#/kg/s]
      NACT,  &  ! Cloud droplet formation by aerosol activation. Adds to Ncm [#/kg/s]
      SIZEFIX_NR, &  ! Adjustment to rain drop number concentration for large/small drops 
      SIZEFIX_NC, &  ! Adjustment to cloud drop number concentration for large/small drops
      SIZEFIX_NI, &  ! Adjustment to ice number concentration for large/small drops
      SIZEFIX_NS, &  ! Adjustment to snow number concentration for large/small drops
      SIZEFIX_NG, &  ! Adjustment to graupel number concentration for large/small drops
      NEGFIX_NI, &    ! Removal of negative ice number concentration 
      NEGFIX_NS, &    ! Removal of negative snow number concentration 
      NEGFIX_NC, &    ! Removal of negative cloud drop number concentration 
      NEGFIX_NR, &    ! Removal of negative rain drop number concentration 
      NEGFIX_NG, &    ! Removal of negative graupel number concentration 
      NIM_MORR_CL, &   ! Clipping of large ice number concentrations
      QC_INST, & ! Change in cloud mixing ratio due to instantaneous processes
      QR_INST, & ! Change in rain mixing ratio due to instantaneous processes
      QI_INST, & ! Change in ice mixing ratio due to instantaneous processes
      QS_INST, & ! Change in snow mixing ratio due to instantaneous processes
      QG_INST, & ! Change in graupel mixing ratio due to instantaneous processes
      NC_INST, & ! Change in cloud number concentration due to instantaneous processes
      NR_INST, & ! Change in rain number concentration due to instantaneous processes
      NI_INST, & ! Change in ice number concentration due to instantaneous processes
      NS_INST, & ! Change in snow number concentration due to instantaneous processes
      NG_INST    ! Change in graupel number concentration due to instantaneous processes

    real( kind = core_rknd ), dimension(nz) :: & ! Temporary variables 
      rrm,    & ! Mean rain water mixing ratio            [kg/kg]
      rim,     & ! Mean ice mixing ratio                   [kg/kg]
      rsm,    & ! Mean snow mixing ratio                  [kg/kg]
      rgm    ! Mean graupel mixing ratio               [kg/kg]

    real :: Morr_snow_rate, Morr_precip_rate

    real, dimension(nz,hydromet_dim) :: &
      hydromet_r4,    & ! Temporary variable
      hydromet_mc_r4

    real, dimension(nz) :: & ! Temporary variables
      rrm_r4,       & ! Mean rain water mixing ratio            [kg/kg]
      Nrm_r4,          & ! Mean rain drop concentration            [num/kg]
      rim_r4,        & ! Mean ice mixing ratio                   [kg/kg]
      Nim_r4,          & ! Mean ice crystal concentration          [num/kg]
      rsm_r4,       & ! Mean snow mixing ratio                  [kg/kg]
      Nsm_r4,       & ! Mean snow flake concentration           [num/kg]
      rgm_r4,    & ! Mean graupel mixing ratio               [kg/kg]
      Ngm_r4,    & ! Mean graupel concentration              [num/kg]
      rrm_mc_r4,    & ! Mean rain water mixing ratio tendency   [kg/kg/s]
      Nrm_mc_r4,       & ! Mean rain drop concentration tendency   [num/kg/s]
      rim_mc_r4,     & ! Mean ice mixing ratio tendency          [kg/kg/s]
      Nim_mc_r4,       & ! Mean ice crystal concentration tendency [num/kg/s]
      rsm_mc_r4,    & ! Mean snow mixing ratio tendency         [kg/kg/s]
      Nsm_mc_r4,    & ! Mean snow flake concentration tendency  [num/kg/s]
      rgm_mc_r4, & ! Mean graupel mixing ratio tendency      [kg/kg/s]
      Ngm_mc_r4    ! Mean graupel concentration tendency     [num/kg/s]

    real, dimension(nz) :: &
      rcm_mc_r4,    &
      rvm_mc_r4,    &
      P_in_pa_r4,   &
      rho_r4,       &
      dzq_r4,       &
      wm_zt_r4,     & ! Mean vertical velocity on the thermo grid  [m/s]
      w_std_dev_r4    ! Standard deviation of w                    [m/s]

    integer :: i, k

    real( kind = core_rknd ), dimension(nz) :: &
      Nrm_auto, & ! Change in Nrm due to autoconversion               [num/kg/s]
      Nrm_evap    ! Change in Nrm due to evaporation                  [num/kg/s]

    ! Local Variables
    real( kind = core_rknd ) :: rsm_sd_morr_int

    real( kind = core_rknd ), dimension(nz) :: &
      hl_before, &
      qto_before, &
      hl_after, &
      hl_on_Cp_residual, &
      qto_after, &
      qto_residual


    ! ---- Begin Code ----

    ! Get rid of compiler warnings
    if ( .false. ) then
       print *, "chi = ", chi
    endif

    call microphys_stats_alloc( nz, num_output_zt, microphys_stats_zt )
    call microphys_stats_alloc( 1, num_output_sfc, microphys_stats_sfc )


    ! Determine temperature
    T_in_K = real( thlm2T_in_K( thlm, exner, rcm ) )

    if ( l_latin_hypercube ) then
      ! Don't use sgs cloud fraction to weight the tendencies
      cloud_frac_in(1:nz) = 0.0

      wm_zt_r4 = real( max( wm_zt, w_thresh ) ) ! Impose a minimum value on w
      w_std_dev_r4 = 0. ! Don't add in a standard deviation for aerosol activation

    else 
      ! Use sgs cloud fraction to weight tendencies
      cloud_frac_in(1:nz) = real( cloud_frac(1:nz) )

      wm_zt_r4 = real( wm_zt ) ! Use the mean value without a threshold
      w_std_dev_r4 = real( w_std_dev ) ! Add in a standard deviation
    end if

    rcm_r4 = real( rcm )
    rvm_r4 = real( rvm )

    ! Note: The Ncm_r4 variable is only used if INUM = 0 in the Morrison code;
    ! otherwise the NDCNST variable is used as the fixed value.
    Ncm_r4 = real( Ncm )
    
    ! Eric Raut added to remove compiler warning (initial value not used)
    i = 1

    forall ( i = 1:hydromet_dim )
      hydromet_r4(1:nz,i) = real( hydromet(1:nz,i) )
    end forall

    if ( l_evaporate_cold_rcm ) then
       ! Convert liquid to vapor at temperatures colder than -37C
       where ( T_in_K < 236.15 )
          rcm_r4 = 0.0
          cloud_frac_in = 0.0
          Ncm_r4 = 0.0
       end where
    end if
    
    ! Initialize tendencies to zero
    T_in_K_mc(1:nz) = 0.0
    rcm_mc(1:nz) = 0.0_core_rknd
    rvm_mc(1:nz) = 0.0_core_rknd
    hydromet_mc(1:nz,:) = 0.0_core_rknd
    Ncm_mc(1:nz) = 0.0_core_rknd
    rcm_sten = 0.0
    rrm_sten(1:nz) = 0.0
    rim_sten(1:nz) = 0.0
    rsm_sten(1:nz) = 0.0
    rgm_sten(1:nz) = 0.0

    ! Initialize effective radius to zero
    effc = 0.0
    effi = 0.0
    effg = 0.0
    effs = 0.0
    effr = 0.0

    ! Initialize Morrison budgets to zero
    PSMLT = 0.0
    EVPMS = 0.0 
    PRACS = 0.0
    EVPMG = 0.0
    PRACG = 0.0
    PGMLT = 0.0
    MNUCCC = 0.0
    PSACWS = 0.0
    PSACWI = 0.0
    QMULTS = 0.0
    QMULTG = 0.0
    PSACWG = 0.0 
    PGSACW = 0.0
    PRD = 0.0
    PRCI = 0.0
    PRAI = 0.0
    QMULTR = 0.0
    QMULTRG = 0.0
    MNUCCD = 0.0
    PRACI = 0.0
    PRACIS = 0.0
    EPRD = 0.0
    MNUCCR = 0.0
    PIACR = 0.0
    PIACRS = 0.0
    PGRACS = 0.0
    PRDS = 0.0
    EPRDS = 0.0
    PSACR = 0.0
    PRDG = 0.0
    EPRDG = 0.0

    ! Initialize these Morrison budgets to zero too
    NGSTEN = 0.0
    NRSTEN = 0.0
    NISTEN = 0.0
    NSSTEN = 0.0
    NCSTEN = 0.0
    NPRC1 = 0.0
    NRAGG = 0.0
    NPRACG = 0.0
    NSUBR = 0.0
    NSMLTR = 0.0
    NGMLTR = 0.0
    NPRACS = 0.0
    NNUCCR = 0.0
    NIACR = 0.0
    NIACRS = 0.0
    NGRACS = 0.0
    NSMLTS = 0.0
    NSAGG = 0.0
    NPRCI = 0.0
    NSCNG = 0.0
    NSUBS = 0.0
    PRA = 0.0
    PRC = 0.0
    PRE = 0.0
    PCC = 0.0 
    NNUCCC = 0.0
    NPSACWS = 0.0
    NPRA = 0.0
    NPRC = 0.0
    NPSACWI = 0.0
    NPSACWG = 0.0
    NPRAI = 0.0
    NMULTS = 0.0
    NMULTG = 0.0
    NMULTR = 0.0
    NMULTRG = 0.0
    NNUCCD = 0.0
    NSUBI = 0.0
    NGMLTG = 0.0
    NSUBG = 0.0
    NACT = 0.0
    SIZEFIX_NR = 0.0
    SIZEFIX_NC = 0.0 
    SIZEFIX_NI = 0.0 
    SIZEFIX_NS = 0.0 
    SIZEFIX_NG = 0.0 
    NEGFIX_NR = 0.0 
    NEGFIX_NC = 0.0 
    NEGFIX_NI = 0.0 
    NEGFIX_NS = 0.0 
    NEGFIX_NG = 0.0 
    NIM_MORR_CL = 0.0
    QC_INST = 0.0
    QR_INST = 0.0
    QI_INST = 0.0
    QS_INST = 0.0
    QG_INST = 0.0
    NC_INST = 0.0
    NR_INST = 0.0
    NI_INST = 0.0
    NS_INST = 0.0
    NG_INST = 0.0


    hydromet_mc_r4 = real( hydromet_mc )
    Ncm_mc_r4 = real( Ncm_mc )
    rcm_mc_r4 = real( rcm_mc )
    rvm_mc_r4 = real( rvm_mc )
    P_in_pa_r4 = real( P_in_Pa )
    rho_r4 = real( rho )
    dzq_r4 = real( dzq )


    ! Unpack hydrometeor arrays.
    rrm = hydromet(:,iirr)

    rrm_r4 = hydromet_r4(:,iirr)
    Nrm_r4    = hydromet_r4(:,iiNr)

    rrm_mc_r4 = hydromet_mc_r4(:,iirr)
    Nrm_mc_r4    = hydromet_mc_r4(:,iiNr)

    if ( l_ice_microphys ) then

       rim  = hydromet(:,iiri)
       rsm = hydromet(:,iirs)

       rim_r4  = hydromet_r4(:,iiri)
       Nim_r4    = hydromet_r4(:,iiNi)
       rsm_r4 = hydromet_r4(:,iirs)
       Nsm_r4 = hydromet_r4(:,iiNs)

       rim_mc_r4  = hydromet_mc_r4(:,iiri)
       Nim_mc_r4    = hydromet_mc_r4(:,iiNi)
       rsm_mc_r4 = hydromet_mc_r4(:,iirs)
       Nsm_mc_r4 = hydromet_mc_r4(:,iiNs)

       if ( l_graupel ) then

          rgm = hydromet(:,iirg)

          rgm_r4 = hydromet_r4(:,iirg)
          Ngm_r4 = hydromet_r4(:,iiNg)

          rgm_mc_r4 = hydromet_mc_r4(:,iirg)
          Ngm_mc_r4 = hydromet_mc_r4(:,iiNg)

       else ! l_graupel disabled

          rgm = zero

          rgm_r4 = 0.0
          Ngm_r4 = 0.0

          rgm_mc_r4 = 0.0
          Ngm_mc_r4 = 0.0

       endif ! l_graupel

    else ! l_ice_microphys disabled

       rim     = zero
       rsm    = zero
       rgm = zero

       rim_r4     = 0.0
       Nim_r4       = 0.0
       rsm_r4    = 0.0
       Nsm_r4    = 0.0
       rgm_r4 = 0.0
       Ngm_r4 = 0.0

       rim_mc_r4     = 0.0
       Nim_mc_r4       = 0.0
       rsm_mc_r4    = 0.0
       Nsm_mc_r4    = 0.0
       rgm_mc_r4 = 0.0
       Ngm_mc_r4 = 0.0

    endif ! l_ice_microphys

    hl_before = Cp * real( T_in_K, kind = core_rknd ) + grav * gr%zt &
                - Lv * ( rcm + rrm ) &
                - Ls * ( rim + rsm + rgm )

    qto_before = rvm + rcm + rrm + rim + rsm + rgm

    ! Call the Morrison microphysics
    call M2005MICRO_GRAUPEL &
         ( rcm_mc_r4, rim_mc_r4, rsm_mc_r4, &
           rrm_mc_r4, Ncm_mc_r4, &
           Nim_mc_r4, Nsm_mc_r4, &
           Nrm_mc_r4, rcm_r4, rim_r4, &
           rsm_r4, rrm_r4, Ncm_r4, &
           Nim_r4, Nsm_r4, Nrm_r4, &
           T_in_K_mc, rvm_mc_r4, T_in_K, rvm_r4, P_in_pa_r4, rho_r4, dzq_r4, &
           wm_zt_r4, w_std_dev_r4, morr_rain_vel_r4, &
           Morr_precip_rate, Morr_snow_rate, effc, effi, effs, effr, real( dt ), &
           1,1, 1,1, 1,nz, 1,1, 1,1, 2,nz, &
           rgm_mc_r4, Ngm_mc_r4, &
           rgm_r4, Ngm_r4, effg, &
           rgm_sten, rrm_sten, &
           rim_sten, rsm_sten, &
           rcm_sten, &
           NGSTEN, NRSTEN, NISTEN, NSSTEN, NCSTEN, &
           cloud_frac_in, &
           PRC, PRA, &
           PSMLT, EVPMS, PRACS, EVPMG, &
           PRACG, PRE,PGMLT,MNUCCC, &
           PSACWS, PSACWI, QMULTS, QMULTG, &
           PSACWG, PGSACW, PRD, PRCI, &
           PRAI, QMULTR, QMULTRG, MNUCCD, &
           PRACI, PRACIS, EPRD, MNUCCR, &
           PIACR, PIACRS, PGRACS, PRDS, &
           EPRDS, PSACR, PRDG, EPRDG, &
           NPRC1, NRAGG, NPRACG, NSUBR, NSMLTR, NGMLTR, NPRACS, NNUCCR, NIACR, &
           NIACRS, NGRACS, NSMLTS, NSAGG, NPRCI, NSCNG, NSUBS, &
           PCC, NNUCCC, NPSACWS, NPRA, NPRC, NPSACWI, NPSACWG, NPRAI, &
           NMULTS, NMULTG, NMULTR, NMULTRG, NNUCCD, NSUBI, NGMLTG, NSUBG, NACT, &
           SIZEFIX_NR, SIZEFIX_NC, SIZEFIX_NI, SIZEFIX_NS, SIZEFIX_NG, &
           NEGFIX_NR, NEGFIX_NC, NEGFIX_NI, NEGFIX_NS, NEGFIX_NG, &
           NIM_MORR_CL, QC_INST, QR_INST, QI_INST, QS_INST, QG_INST, &
           NC_INST, NR_INST, NI_INST, NS_INST, NG_INST )

    hl_after = Cp * real( T_in_K, kind = core_rknd ) + grav * gr%zt &
               - Lv * ( real( rcm_r4, kind = core_rknd) &
                        + real( rrm_r4, kind = core_rknd ) ) &
               - Ls * ( real( rim_r4, kind = core_rknd ) &
                        + real( rsm_r4, kind = core_rknd ) &
                        + real( rgm_r4, kind = core_rknd ) )

    hl_on_Cp_residual &
    = ( hl_after - hl_before &
        - dt * Lv * ( real( rcm_sten, kind = core_rknd ) &
                      + real( rrm_sten, kind = core_rknd ) ) &
        - dt * Ls * ( real( rim_sten, kind = core_rknd ) &
                      + real( rsm_sten, kind = core_rknd ) &
                      + real( rgm_sten, kind = core_rknd ) ) ) / Cp

    qto_after = real( rvm_r4, kind = core_rknd ) &
                + real( rcm_r4, kind = core_rknd ) &
                + real( rrm_r4, kind = core_rknd ) &
                + real( rim_r4, kind = core_rknd ) &
                + real( rsm_r4, kind = core_rknd ) &
                + real( rgm_r4, kind = core_rknd )

    qto_residual = qto_after - qto_before &
                   - dt * ( real( rcm_sten, kind = core_rknd ) &
                            + real( rrm_sten, kind = core_rknd ) &
                            + real( rim_sten, kind = core_rknd ) &
                            + real( rsm_sten, kind = core_rknd ) &
                            + real( rgm_sten, kind = core_rknd) )

    ! Pack hydrometeor arrays.
    hydromet_r4(:,iirr) = rrm_r4
    hydromet_r4(:,iiNr)    = Nrm_r4

    hydromet_mc_r4(:,iirr) = rrm_mc_r4
    hydromet_mc_r4(:,iiNr)    = Nrm_mc_r4

    if ( l_ice_microphys ) then

       hydromet_r4(:,iiri)  = rim_r4
       hydromet_r4(:,iiNi)    = Nim_r4
       hydromet_r4(:,iirs) = rsm_r4
       hydromet_r4(:,iiNs) = Nsm_r4

       hydromet_mc_r4(:,iiri)  = rim_mc_r4
       hydromet_mc_r4(:,iiNi)    = Nim_mc_r4
       hydromet_mc_r4(:,iirs) = rsm_mc_r4
       hydromet_mc_r4(:,iiNs) = Nsm_mc_r4

       if ( l_graupel ) then

          hydromet_r4(:,iirg) = rgm_r4
          hydromet_r4(:,iiNg) = Ngm_r4

          hydromet_mc_r4(:,iirg) = rgm_mc_r4
          hydromet_mc_r4(:,iiNg) = Ngm_mc_r4

       endif

    endif

    !hydromet_mc = real( hydromet_mc_r4, kind = core_rknd )
    rcm_mc = real( rcm_mc_r4, kind = core_rknd )
    rvm_mc = real( rvm_mc_r4, kind = core_rknd )

    rrm_auto = real( PRC, kind = core_rknd )
    rrm_accr = real( PRA, kind = core_rknd )
    rrm_evap = real( PRE, kind = core_rknd )

    Nrm_auto = real( NPRC1, kind = core_rknd )
    Nrm_evap = real( NSUBR, kind = core_rknd )
    ! Update hydrometeor tendencies
    ! This done because the hydromet_mc arrays that are produced by
    ! M2005MICRO_GRAUPEL don't include the clipping term.
    do i = 1, hydromet_dim, 1
       hydromet_mc(2:,i) = ( real( hydromet_r4(2:,i), kind = core_rknd ) - &
                             hydromet(2:,i) ) / dt
    end do

    hydromet_mc(1,:) = 0.0_core_rknd ! Boundary condition

    Ncm_mc(2:) = ( real( Ncm_r4(2:), kind = core_rknd ) - Ncm(2:) ) &
                 / dt

    Ncm_mc(1) = 0.0_core_rknd ! Boundary condition

    ! Update thetal based on absolute temperature
    thlm_mc = ( T_in_K2thlm( real( T_in_K, kind = core_rknd ), exner, &
                real( rcm_r4, kind = core_rknd ) ) - thlm ) / dt

    ! Sedimentation is handled within the Morrison microphysics
    hydromet_vel_zt(:,:) = 0.0_core_rknd

    ! Output rain sedimentation velocity
    ! Multiply by -1 so that negative is associated with falling precip
    morr_rain_vel_r4(:) = morr_rain_vel_r4(:) * (-1.0)
    do k = 1, nz, 1
      hydromet_vel_zt(k,iirr) = real( morr_rain_vel_r4(k), kind = core_rknd )
    end do
    

    rsm_sd_morr_int = vertical_integral( (nz - 2 + 1), rho_ds_zt(2:nz), &
                            real( rsm_sten(2:nz), kind=core_rknd ), &
                            gr%dzt(2:nz) )

    call microphys_put_var( irsm_sd_morr_int, (/rsm_sd_morr_int/), microphys_stats_sfc )

    if ( clubb_at_least_debug_level( 1 ) ) then
        if ( rsm_sd_morr_int > maxval( real( rsm_sten(2:nz), &
                                               kind=core_rknd ) ) ) then
            print *, "Warning: rsm_sd_morr was not conservative!" // &
                   " rsm_sd_morr_verical_integr = ", rsm_sd_morr_int
        endif
    endif
    
    call microphys_put_var( irrm_auto, rrm_auto, microphys_stats_zt )
    call microphys_put_var( irrm_accr, rrm_accr, microphys_stats_zt )
    call microphys_put_var( irrm_evap, rrm_evap, microphys_stats_zt )

    call microphys_put_var( iNrm_auto,    Nrm_auto,    microphys_stats_zt )
    call microphys_put_var( iNrm_evap,    Nrm_evap,    microphys_stats_zt )

    ! Update Morrison budgets
    call microphys_put_var( ihl_on_Cp_residual, hl_on_Cp_residual, microphys_stats_zt )
    call microphys_put_var( iqto_residual, qto_residual, microphys_stats_zt )
    call microphys_put_var( irgm_sd_morr, &
              real( rgm_sten, kind = core_rknd ), microphys_stats_zt )
    call microphys_put_var( irrm_sd_morr, &
              real( rrm_sten, kind = core_rknd ), microphys_stats_zt )
    call microphys_put_var( irsm_sd_morr, &
              real( rsm_sten, kind = core_rknd ), microphys_stats_zt )
    call microphys_put_var( irim_sd_mg_morr, &
              real( rim_sten, kind = core_rknd ), microphys_stats_zt )
    call microphys_put_var( ircm_sd_mg_morr, &
              real( rcm_sten, kind = core_rknd), microphys_stats_zt )
    call microphys_put_var( iPRC,real(PRC,kind=core_rknd),microphys_stats_zt )
    call microphys_put_var( iPRA,real(PRA,kind=core_rknd),microphys_stats_zt )
    call microphys_put_var( iPRE,real(PRE,kind=core_rknd),microphys_stats_zt )
    call microphys_put_var( iPSMLT, real( PSMLT, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iEVPMS, real( EVPMS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPRACS, real( PRACS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iEVPMG, real( EVPMG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPRACG, real( PRACG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPGMLT, real( PGMLT, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iMNUCCC, real( MNUCCC, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPSACWS, real( PSACWS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPSACWI, real( PSACWI, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iQMULTS, real( QMULTS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iQMULTG, real( QMULTG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPSACWG, real( PSACWG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPGSACW, real( PGSACW, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPRD, real( PRD, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPRCI, real( PRCI, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPRAI, real( PRAI, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iQMULTR, real( QMULTR, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iQMULTRG, real( QMULTRG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iMNUCCD, real( MNUCCD, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPRACI, real( PRACI, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPRACIS, real( PRACIS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iEPRD, real( EPRD, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iMNUCCR, real( MNUCCR, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPIACR, real( PIACR, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPIACRS, real( PIACRS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPGRACS, real( PGRACS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPRDS, real( PRDS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iEPRDS, real( EPRDS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPSACR, real( PSACR, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iPRDG, real( PRDG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iEPRDG, real( EPRDG, kind=core_rknd ), microphys_stats_zt )
      
    ! Update more Morrison budgets
    call microphys_put_var( iNGSTEN, real( NGSTEN, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNRSTEN, real( NRSTEN, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNISTEN, real( NISTEN, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNSSTEN, real( NSSTEN, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNCSTEN, real( NCSTEN, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNPRC1, real( NPRC1, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNRAGG, real( NRAGG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNPRACG, real( NPRACG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNSUBR, real( NSUBR, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNSMLTR, real( NSMLTR, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNGMLTR, real( NGMLTR, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNPRACS, real( NPRACS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNNUCCR, real( NNUCCR, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNIACR, real( NIACR, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNIACRS, real( NIACRS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNGRACS, real( NGRACS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNSMLTS, real( NSMLTS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNSAGG, real( NSAGG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNPRCI, real( NPRCI, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNSCNG, real( NSCNG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNSUBS, real( NSUBS, kind=core_rknd ), microphys_stats_zt )

    call microphys_put_var( iPCC, real( PCC, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNNUCCC, real( NNUCCC, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNPSACWS, real( NPSACWS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNPRA, real( NPRA, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNPRC, real( NPRC, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNPSACWI, real( NPSACWI, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNPSACWG, real( NPSACWG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNPRAI, real( NPRAI, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNMULTS, real( NMULTS, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNMULTG, real( NMULTG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNMULTR, real( NMULTR, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNMULTRG, real( NMULTRG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNNUCCD, real( NNUCCD, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNSUBI, real( NSUBI, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNGMLTG, real( NGMLTG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNSUBG, real( NSUBG, kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iNACT, real( NACT,kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iSIZEFIX_NR, real( SIZEFIX_NR,kind=core_rknd),microphys_stats_zt )
    call microphys_put_var( iSIZEFIX_NC, real( SIZEFIX_NC,kind=core_rknd),microphys_stats_zt )
    call microphys_put_var( iSIZEFIX_NI, real( SIZEFIX_NI,kind=core_rknd),microphys_stats_zt )
    call microphys_put_var( iSIZEFIX_NS, real( SIZEFIX_NS,kind=core_rknd),microphys_stats_zt )
    call microphys_put_var( iSIZEFIX_NG, real( SIZEFIX_NG,kind=core_rknd),microphys_stats_zt )
    call microphys_put_var( iNEGFIX_NR, real( NEGFIX_NR,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iNEGFIX_NC, real( NEGFIX_NC,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iNEGFIX_NI, real( NEGFIX_NI,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iNEGFIX_NS, real( NEGFIX_NS,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iNEGFIX_NG, real( NEGFIX_NG,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iNIM_MORR_CL, real( NIM_MORR_CL,kind=core_rknd ), microphys_stats_zt )
    call microphys_put_var( iQC_INST, real( QC_INST,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iQR_INST, real( QR_INST,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iQI_INST, real( QI_INST,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iQS_INST, real( QS_INST,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iQG_INST, real( QG_INST,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iNC_INST, real( NC_INST,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iNR_INST, real( NR_INST,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iNI_INST, real( NI_INST,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iNS_INST, real( NS_INST,kind=core_rknd ),microphys_stats_zt )
    call microphys_put_var( iNG_INST, real( NG_INST,kind=core_rknd ),microphys_stats_zt )

    call microphys_put_var( iT_in_K_mc, real( T_in_K_mc, kind=core_rknd ),microphys_stats_zt )

    ! --- Number concentrations ---
    ! No budgets for sedimentation are output
    ! Effective radii of hydrometeor species
    call microphys_put_var( ieff_rad_cloud, real( effc(:), kind = core_rknd ), microphys_stats_zt )
    call microphys_put_var( ieff_rad_ice, real( effi(:), kind = core_rknd ), microphys_stats_zt )
    call microphys_put_var( ieff_rad_snow, real( effs(:), kind = core_rknd ), microphys_stats_zt )
    call microphys_put_var( ieff_rad_rain, real( effr(:), kind = core_rknd ), microphys_stats_zt )
    call microphys_put_var( ieff_rad_graupel, real( effg(:), kind=core_rknd ), microphys_stats_zt )

      ! Snow and Rain rates at the bottom of the domain, in mm/day
    call microphys_put_var( iprecip_rate_sfc, &
      (/ real(Morr_precip_rate, kind = core_rknd) * &
      sec_per_day / dt /), microphys_stats_sfc )

    call microphys_put_var( imorr_snow_rate, &
      (/ real(Morr_snow_rate, kind = core_rknd) * &
      sec_per_day / dt /), microphys_stats_sfc )

    return
  end subroutine morrison_microphys_driver

end module morrison_microphys_module
