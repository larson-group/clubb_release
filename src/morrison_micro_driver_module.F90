! $Id$
module morrison_micro_driver_module

  implicit none

  public :: morrison_micro_driver

  private

  contains
!-------------------------------------------------------------------------------
  subroutine morrison_micro_driver &
             ( dt, nz, l_stats_samp, &
               l_latin_hypercube, thlm, wm_zt, p_in_Pa, &
               exner, rho, cloud_frac, w_std_dev, &
               dzq, rcm, Ncm, s_mellor, rvm, hydromet, lh_stat_sample_weight, &
               hydromet_mc, hydromet_vel_zt, &
               rcm_mc, rvm_mc, thlm_mc, &
               rrainm_auto, rrainm_accr, rrainm_evap, &
               Nrm_auto, Nrm_evap )

! Description:
!   Wrapper for the Morrison microphysics
!
! References:
!   None
!-------------------------------------------------------------------------------


    use parameters_model, only: hydromet_dim

    ! The version of the Morrison 2005 microphysics that is in SAM.
    use module_MP_graupel, only: &
      M2005MICRO_GRAUPEL  ! Procedure

    use module_MP_graupel, only: &
      cloud_frac_thresh ! Constant

    use stats_variables, only: &
      zt,  & ! Variables
      LH_sfc, &
      sfc

    use stats_variables, only: & 
      irsnowm_sd_morr, & ! Variables
      iricem_sd_mg_morr, & 
      irrainm_sd_morr, & 
      irsnowm_sd_morr, &
      irgraupelm_sd_morr, &
      ircm_sd_mg_morr, &
      ircm_in_cloud, &
      irrainm_auto, &
      irrainm_accr, &
      irrainm_cond

    use stats_variables, only: & 
      ieff_rad_cloud, & ! Variables
      ieff_rad_ice, &
      ieff_rad_snow, &
      ieff_rad_rain, &
      ieff_rad_graupel

    use stats_variables, only: & 
      iprecip_rate_sfc, & !Variables
      imorr_snow_rate, &
      iLH_morr_snow_rate

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
      iNGSTEN, & ! Lots of variable(s)
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
      iPCC, &
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


    use stats_type, only:  & 
        stat_update_var, stat_update_var_pt  ! Procedure(s)

    use T_in_K_module, only: &
      T_in_K2thlm, & ! Procedure(s)
      thlm2T_in_K

    use array_index, only:  & 
      iirrainm, iirsnowm, iiricem, iirgraupelm, &
      iiNrm, iiNsnowm, iiNim, iiNgraupelm, iiNcm

    use constants_clubb, only: &
      sec_per_day, &
      zero

    use clubb_precision, only: &
      core_rknd, & ! Variable(s)
      time_precision

    implicit none

    ! External
    intrinsic :: max, real

    ! Constant parameters
    real( kind = core_rknd ), parameter :: &
      w_thresh = 0.1_core_rknd ! Minimum value w for latin hypercube [m/s]

    ! Input Variables
    real( kind = time_precision ), intent(in) :: dt ! Model timestep        [s]

    integer, intent(in) :: nz ! Points in the Vertical        [-]

    logical, intent(in) :: &
      l_stats_samp,     & ! Whether to accumulate statistics [T/F]
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
      s_mellor,     & ! The variable 's' from Mellor              [kg/kg]
      rvm             ! Vapor water mixing ratio                  [kg/kg]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), intent(in) :: &
      lh_stat_sample_weight

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), target, intent(out) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel_zt ! Hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rrainm_auto,     & ! Autoconversion rate                 [kg/kg/s]
      rrainm_accr,     & ! Accretion rate                      [kg/kg/s]
      rrainm_evap        ! Rain evaporation rate               [kg/kg/s]

    ! Local Variables
    real, dimension(nz) :: & 
      effc, effi, effg, effs, effr ! Effective droplet radii [μ]

    real, dimension(nz) :: & 
      T_in_K,        & ! Temperature                        [K]
      T_in_K_mc,     & ! Temperature tendency               [K/s]
      rcm_r4,        & ! Temporary array for cloud water mixing ratio  [kg/kg]
      Ncm_r4,        & ! Temporary array for cloud number conc.        [#/kg]
      Ncm_mc_r4,     & ! Temporary array for cloud number conc.        [#/kg/s]
      rvm_r4,        & ! Temporary array for vapor water mixing ratio  [kg/kg]
      rcm_sten,      & ! Cloud dropet sedimentation tendency           [kg/kg/s]
      morr_rain_vel_r4, & ! Rain fall velocity from Morrison microphysics [m/s] 
      cloud_frac_in    ! Cloud fraction used as input for the Morrison scheme [-]

    ! In the comments below, by "adds to" we mean that if the quantity is
    ! positive, it adds positively to the prognostic variable, but if the
    ! quantity is negative, it subtracts from the prognostic variable.
    real, dimension(nz) :: &
      PSMLT,  & ! Freezing of rain to form snow.
                !    Adds to rsnowm, subtracts from rrainm [kg/kg/s]
      EVPMS,  & ! Evaporation of melted snow.
                !    Adds to rsnowm, subtracts from rvm [kg/kg/s]
      PRACS,  & ! Collection of rain by snow.
                !    Adds to rsnowm, subtracts from rrainm [kg/kg/s]
      EVPMG,  & ! Evaporation of melted graupel.
                !    Adds to rgraupelm, subtracts from rvm [kg/kg/s]
      PRACG,  & ! Negative of collection of rain by graupel.
                !    Adds to rrainm, subtracts from rgraupelm [kg/kg/s]
      PGMLT,  & ! Negative of melting of graupel.
                !    Adds to rgraupelm, subtracts from rrainm [kg/kg/s]
      MNUCCC, & ! Contact freezing of cloud droplets.
                !    Adds to ricem, subtracts from rcm [kg/kg/s]
      PSACWS, & ! Collection of cloud water by snow.
                !    Adds to rsnowm, subtracts from rcm [kg/kg/s]
      PSACWI, & ! Collection of cloud water by cloud ice.
                !    Adds to ricem, subtracts from rcm [kg/kg/s]
      QMULTS, & ! Splintering from cloud droplets accreted onto snow.
                !    Adds to ricem, subtracts from rcm [kg/kg/s]
      QMULTG, & ! Splintering from droplets accreted onto graupel.
                !    Adds to ricem, subtracts from rcm [kg/kg/s]
      PSACWG, & ! Collection of cloud water by graupel.
                !    Adds to rgraupelm, subtracts from rcm [kg/kg/s]
      PGSACW, & ! Reclassification of rimed snow as graupel.
                !    Adds to rgraupelm, subtracts from rcm [kg/kg/s]
      PRD,    & ! Depositional growth of cloud ice.
                !    Adds to ricem, subtracts from rcm [kg/kg/s]
      PRCI,   & ! Autoconversion of cloud ice to snow.
                !    Adds to rsnowm, subtracts from ricem [kg/kg/s]
      PRAI,   & ! Collection of cloud ice by snow.
                !    Adds to rsnowm, subtracts from ricem [kg/kg/s]
      QMULTR, & ! Splintering from rain droplets accreted onto snow.
                !    Adds to ricem, subtracts from rrainm [kg/kg/s]
      QMULTRG,& ! Splintering from rain droplets accreted onto graupel.
                !    Adds to ricem, subtracts from rrainm [kg/kg/s]
      MNUCCD, & ! Freezing of aerosol.
                !    Adds to ricem, subtracts from rvm [kg/kg/s]
      PRACI,  & ! Collection of cloud ice by rain.
                !    Adds to rgraupelm, subtracts from ricem [kg/kg/s]
      PRACIS, & ! Collection of cloud ice by rain.
                !    Adds to rsnowm, subtracts from ricem [kg/kg/s]
      EPRD,   & ! Negative of sublimation of cloud ice.
                !    Adds to ricem, subtracts from rvm [kg/kg/s]
      MNUCCR, & ! Contact freezing of rain droplets.
                !    Adds to rgraupelm, subtracts from rrainm [kg/kg/s]
      PIACR,  & ! Collection of cloud ice by rain.
                !    Adds to rgraupelm, subtracts from rrainm [kg/kg/s]
      PIACRS, & ! Collection of cloud ice by rain.
                !    Adds to rsnowm, subtracts from rrainm [kg/kg/s]
      PGRACS, & ! Collection of rain by snow.
                !    Adds to rgraupelm, subtracts from rrainm [kg/kg/s]
      PRDS,   & ! Depositional growth of snow.
                !    Adds to rsnowm, subtracts from rvm [kg/kg/s]
      EPRDS,  & ! Negative of sublimation of snow.
                !    Adds to rsnowm, subtracts from rvm [kg/kg/s]
      PSACR,  & ! Collection of snow by rain.
                !    Adds to rgraupelm, subtracts from rsnowm [kg/kg/s]
      PRDG,   & ! Depositional growth of graupel.
                !    Adds to rgraupelm, subtracts from rvm [kg/kg/s]
      EPRDG     ! Negative of sublimation of graupel.
                !    Adds to rgraupelm, subtracts from rvm [kg/kg/s]

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
      NNUCCR, & ! Contact freezing of rain. Adds to Ngraupelm, subtracts from Nrm [#/kg/s]
      NIACR,  & ! Collection of cloud ice by rain.
                !    Adds to Ngraupelm, subtracts from Nrm and Nim [#/kg/s] 
      NIACRS, & ! Collection of cloud ice by rain.
                !    Adds to Nsnowm, subtracts from Nrm and Nim [#/kg/s]
      NGRACS, & ! Collection of rain by snow.
                !    Adds to Ngraupelm, subtracts from Nrm and Nsnowm [#/kg/s]
      NSMLTS, & ! Melting of snow
                !    Adds to Nsnowm [#/kg/s]
      NSAGG, &  ! Self collection of snow
                !    Adds to Nsnowm [#/kg/s]
      NPRCI, &  ! Autoconversion of cloud ice to snow
                !    Adds to Nsnowm, subtracts from Nim [#/kg/s]
      NSCNG, &  ! Conversion of snow to graupel
                !    Adds to Ngraupelm, subtracts from Nsnowm [#/kg/s]
      NSUBS, &  ! Loss of Nsnowm due to sublimation
                !    Adds to Nsnowm [#/kg/s]
      PRA,   &  ! Accretion. Adds to rrainm, subtracts from rcm [kg/kg/s]
      PRC,   &  ! Autoconversion. Adds to rrainm, subtracts from rcm [kg/kg/s]
      PRE       ! Rain evaporation. Subtracts from rrainm [kg/kg/s]              
          
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
      NGMLTG, & ! Loss of graupel due to melting. Subtracts from Ngraupelm [#/kg/s]
      NSUBG, &  ! Loss of graupel due to sublimation. Subtracts from Ngraupelm [#/kg/s]
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


    real( kind = core_rknd ), dimension(nz) :: & 
      rcm_in_cloud     ! Liquid water in cloud           [kg/kg]

    real :: Morr_snow_rate, Morr_precip_rate

    real, dimension(nz,hydromet_dim) :: &
      hydromet_r4, &   ! Temporary variable
      hydromet_sten, & ! Hydrometeor sedimentation tendency [(units vary)/s]
      hydromet_mc_r4

    real, dimension(nz) :: &
      rcm_mc_r4, &
      rvm_mc_r4, &
      P_in_pa_r4, &
      rho_r4, &
      dzq_r4, &
      wm_zt_r4, &     ! Mean vertical velocity on the thermo grid   [m/s]
      w_std_dev_r4 ! Standard deviation of w  [m/s]

    integer :: i, k

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      Nrm_auto, & ! Change in Nrm due to autoconversion               [num/kg/s]
      Nrm_evap    ! Change in Nrm due to evaporation                  [num/kg/s]


    ! ---- Begin Code ----

    ! Get rid of compiler warnings
    if ( .false. ) then
       print *, "s_mellor = ", s_mellor
    endif

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
    
    ! Initialize tendencies to zero
    T_in_K_mc(1:nz) = 0.0
    rcm_mc(1:nz) = 0.0_core_rknd
    rvm_mc(1:nz) = 0.0_core_rknd
    hydromet_mc(1:nz,:) = 0.0_core_rknd
    hydromet_sten(1:nz,:) = 0.0
    rcm_sten = 0.0
    Ncm_mc_r4 = 0.0

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
    rcm_mc_r4 = real( rcm_mc )
    rvm_mc_r4 = real( rvm_mc )
    P_in_pa_r4 = real( P_in_Pa )
    rho_r4 = real( rho )
    dzq_r4 = real( dzq )

    ! Call the Morrison microphysics
    call M2005MICRO_GRAUPEL &
         ( rcm_mc_r4, hydromet_mc_r4(:,iiricem), hydromet_mc_r4(:,iirsnowm), &
           hydromet_mc_r4(:,iirrainm), Ncm_mc_r4(:), &
           hydromet_mc_r4(:,iiNim), hydromet_mc_r4(:,iiNsnowm), &
           hydromet_mc_r4(:,iiNrm), rcm_r4, hydromet_r4(:,iiricem), &
           hydromet_r4(:,iirsnowm), hydromet_r4(:,iirrainm), Ncm_r4(:), &
           hydromet_r4(:,iiNim), hydromet_r4(:,iiNsnowm), hydromet_r4(:,iiNrm), &
           T_in_K_mc, rvm_mc_r4, T_in_K, rvm_r4, P_in_pa_r4, rho_r4, dzq_r4, &
           wm_zt_r4, w_std_dev_r4, morr_rain_vel_r4, &
           Morr_precip_rate, Morr_snow_rate, effc, effi, effs, effr, real( dt ), &
           1,1, 1,1, 1,nz, 1,1, 1,1, 2,nz, &
           hydromet_mc_r4(:,iirgraupelm), hydromet_mc_r4(:,iiNgraupelm), &
           hydromet_r4(:,iirgraupelm), hydromet_r4(:,iiNgraupelm), effg, &
           hydromet_sten(:,iirgraupelm), hydromet_sten(:,iirrainm), &
           hydromet_sten(:,iiricem), hydromet_sten(:,iirsnowm), &
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

    !hydromet_mc = real( hydromet_mc_r4, kind = core_rknd )
    rcm_mc = real( rcm_mc_r4, kind = core_rknd )
    rvm_mc = real( rvm_mc_r4, kind = core_rknd )

    rrainm_auto = real( PRC, kind = core_rknd )
    rrainm_accr = real( PRA, kind = core_rknd )
    rrainm_evap = real( PRE, kind = core_rknd )

    Nrm_auto = real( NPRC1, kind = core_rknd )
    Nrm_evap = real( NSUBR, kind = core_rknd )
    ! Update hydrometeor tendencies
    ! This done because the hydromet_mc arrays that are produced by
    ! M2005MICRO_GRAUPEL don't include the clipping term.
    do i = 1, hydromet_dim, 1
      if ( i == iiNcm ) then
        hydromet_mc(2:,iiNcm) = ( real( Ncm_r4(2:), kind = core_rknd ) - &
                                hydromet(2:,iiNcm) ) / real( dt, kind=core_rknd )
      else
        hydromet_mc(2:,i) = ( real( hydromet_r4(2:,i), kind = core_rknd ) - &
                              hydromet(2:,i)  ) / real( dt, kind = core_rknd )
      end if
    end do


    hydromet_mc(1,:) = 0.0_core_rknd ! Boundary condition

    ! Update thetal based on absolute temperature
    thlm_mc = ( T_in_K2thlm( real( T_in_K, kind = core_rknd ), exner, &
                real( rcm_r4, kind = core_rknd ) ) - thlm ) / real( dt, kind = core_rknd )

    ! Sedimentation is handled within the Morrison microphysics
    hydromet_vel_zt(:,:) = 0.0_core_rknd

    ! Output rain sedimentation velocity
    ! Multiply by -1 so that negative is associated with falling precip
    morr_rain_vel_r4(:) = morr_rain_vel_r4(:) * (-1.0)
    do k = 1, nz, 1
      hydromet_vel_zt(k,iirrainm) = real( morr_rain_vel_r4(k), kind = core_rknd )
    end do

    if ( .not. l_latin_hypercube .and. l_stats_samp ) then

      where ( cloud_frac(:) > real( cloud_frac_thresh, kind = core_rknd ) ) 
        rcm_in_cloud(:) = rcm / cloud_frac
      else where
        rcm_in_cloud(:) = rcm
      end where

      call stat_update_var( ircm_in_cloud, rcm_in_cloud, zt )

      call stat_update_var( irrainm_auto, rrainm_auto, zt )
      call stat_update_var( irrainm_accr, rrainm_accr, zt )
      call stat_update_var( irrainm_cond, rrainm_evap, zt )
    end if ! ( .not. l_latin_hypercube .and. l_stats_samp )

    if ( l_stats_samp ) then
      call stat_update_var( irgraupelm_sd_morr, lh_stat_sample_weight  &
                * real( hydromet_sten(:,iirgraupelm), kind = core_rknd ), zt )
      call stat_update_var( irrainm_sd_morr,lh_stat_sample_weight &
                * real( hydromet_sten(:,iirrainm), kind = core_rknd ), zt )
      call stat_update_var( irsnowm_sd_morr,lh_stat_sample_weight &
                * real( hydromet_sten(:,iirsnowm), kind = core_rknd ), zt )
      call stat_update_var( iricem_sd_mg_morr,lh_stat_sample_weight &
                * real( hydromet_sten(:,iiricem), kind = core_rknd ), zt )
      call stat_update_var( ircm_sd_mg_morr,lh_stat_sample_weight &
                * real( rcm_sten, kind = core_rknd), zt )
      call stat_update_var( iPRC,lh_stat_sample_weight*real(PRC,kind=core_rknd),zt)
      call stat_update_var( iPRA,lh_stat_sample_weight*real(PRA,kind=core_rknd),zt)
      call stat_update_var( iPRE,lh_stat_sample_weight*real(PRE,kind=core_rknd),zt)
      call stat_update_var( iPSMLT, lh_stat_sample_weight*real( PSMLT, kind=core_rknd ), zt )
      call stat_update_var( iEVPMS, lh_stat_sample_weight*real( EVPMS, kind=core_rknd ), zt )
      call stat_update_var( iPRACS, lh_stat_sample_weight*real( PRACS, kind=core_rknd ), zt )
      call stat_update_var( iEVPMG, lh_stat_sample_weight*real( EVPMG, kind=core_rknd ), zt )
      call stat_update_var( iPRACG, lh_stat_sample_weight*real( PRACG, kind=core_rknd ), zt )
      call stat_update_var( iPGMLT, lh_stat_sample_weight*real( PGMLT, kind=core_rknd ), zt )
      call stat_update_var( iMNUCCC, lh_stat_sample_weight*real( MNUCCC, kind=core_rknd ), zt )
      call stat_update_var( iPSACWS, lh_stat_sample_weight*real( PSACWS, kind=core_rknd ), zt )
      call stat_update_var( iPSACWI, lh_stat_sample_weight*real( PSACWI, kind=core_rknd ), zt )
      call stat_update_var( iQMULTS, lh_stat_sample_weight*real( QMULTS, kind=core_rknd ), zt )
      call stat_update_var( iQMULTG, lh_stat_sample_weight*real( QMULTG, kind=core_rknd ), zt )
      call stat_update_var( iPSACWG, lh_stat_sample_weight*real( PSACWG, kind=core_rknd ), zt )
      call stat_update_var( iPGSACW, lh_stat_sample_weight*real( PGSACW, kind=core_rknd ), zt )
      call stat_update_var( iPRD, lh_stat_sample_weight*real( PRD, kind=core_rknd ), zt )
      call stat_update_var( iPRCI, lh_stat_sample_weight*real( PRCI, kind=core_rknd ), zt )
      call stat_update_var( iPRAI, lh_stat_sample_weight*real( PRAI, kind=core_rknd ), zt )
      call stat_update_var( iQMULTR, lh_stat_sample_weight*real( QMULTR, kind=core_rknd ), zt )
      call stat_update_var( iQMULTRG, lh_stat_sample_weight*real( QMULTRG, kind=core_rknd ), zt )
      call stat_update_var( iMNUCCD, lh_stat_sample_weight*real( MNUCCD, kind=core_rknd ), zt )
      call stat_update_var( iPRACI, lh_stat_sample_weight*real( PRACI, kind=core_rknd ), zt )
      call stat_update_var( iPRACIS, lh_stat_sample_weight*real( PRACIS, kind=core_rknd ), zt )
      call stat_update_var( iEPRD, lh_stat_sample_weight*real( EPRD, kind=core_rknd ), zt )
      call stat_update_var( iMNUCCR, lh_stat_sample_weight*real( MNUCCR, kind=core_rknd ), zt )
      call stat_update_var( iPIACR, lh_stat_sample_weight*real( PIACR, kind=core_rknd ), zt )
      call stat_update_var( iPIACRS, lh_stat_sample_weight*real( PIACRS, kind=core_rknd ), zt )
      call stat_update_var( iPGRACS, lh_stat_sample_weight*real( PGRACS, kind=core_rknd ), zt )
      call stat_update_var( iPRDS, lh_stat_sample_weight*real( PRDS, kind=core_rknd ), zt )
      call stat_update_var( iEPRDS, lh_stat_sample_weight*real( EPRDS, kind=core_rknd ), zt )
      call stat_update_var( iPSACR, lh_stat_sample_weight*real( PSACR, kind=core_rknd ), zt )
      call stat_update_var( iPRDG, lh_stat_sample_weight*real( PRDG, kind=core_rknd ), zt )
      call stat_update_var( iEPRDG, lh_stat_sample_weight*real( EPRDG, kind=core_rknd ), zt )
        
      ! Update more Morrison budgets
      call stat_update_var( iNGSTEN, lh_stat_sample_weight*real( NGSTEN, kind=core_rknd ), zt )
      call stat_update_var( iNRSTEN, lh_stat_sample_weight*real( NRSTEN, kind=core_rknd ), zt )
      call stat_update_var( iNISTEN, lh_stat_sample_weight*real( NISTEN, kind=core_rknd ), zt )
      call stat_update_var( iNSSTEN, lh_stat_sample_weight*real( NSSTEN, kind=core_rknd ), zt )
      call stat_update_var( iNCSTEN, lh_stat_sample_weight*real( NCSTEN, kind=core_rknd ), zt )
      call stat_update_var( iNPRC1, lh_stat_sample_weight*real( NPRC1, kind=core_rknd ), zt )
      call stat_update_var( iNRAGG, lh_stat_sample_weight*real( NRAGG, kind=core_rknd ), zt )
      call stat_update_var( iNPRACG, lh_stat_sample_weight*real( NPRACG, kind=core_rknd ), zt )
      call stat_update_var( iNSUBR, lh_stat_sample_weight*real( NSUBR, kind=core_rknd ), zt )
      call stat_update_var( iNSMLTR, lh_stat_sample_weight*real( NSMLTR, kind=core_rknd ), zt )
      call stat_update_var( iNGMLTR, lh_stat_sample_weight*real( NGMLTR, kind=core_rknd ), zt )
      call stat_update_var( iNPRACS, lh_stat_sample_weight*real( NPRACS, kind=core_rknd ), zt )
      call stat_update_var( iNNUCCR, lh_stat_sample_weight*real( NNUCCR, kind=core_rknd ), zt )
      call stat_update_var( iNIACR, lh_stat_sample_weight*real( NIACR, kind=core_rknd ), zt )
      call stat_update_var( iNIACRS, lh_stat_sample_weight*real( NIACRS, kind=core_rknd ), zt )
      call stat_update_var( iNGRACS, lh_stat_sample_weight*real( NGRACS, kind=core_rknd ), zt )
      call stat_update_var( iNSMLTS, lh_stat_sample_weight*real( NSMLTS, kind=core_rknd ), zt )
      call stat_update_var( iNSAGG, lh_stat_sample_weight*real( NSAGG, kind=core_rknd ), zt )
      call stat_update_var( iNPRCI, lh_stat_sample_weight*real( NPRCI, kind=core_rknd ), zt )
      call stat_update_var( iNSCNG, lh_stat_sample_weight*real( NSCNG, kind=core_rknd ), zt )
      call stat_update_var( iNSUBS, lh_stat_sample_weight*real( NSUBS, kind=core_rknd ), zt )

      call stat_update_var( iPCC, lh_stat_sample_weight*real( PCC, kind=core_rknd ), zt )
      call stat_update_var( iNNUCCC, lh_stat_sample_weight*real( NNUCCC, kind=core_rknd ), zt )
      call stat_update_var( iNPSACWS, lh_stat_sample_weight*real( NPSACWS, kind=core_rknd ), zt )
      call stat_update_var( iNPRA, lh_stat_sample_weight*real( NPRA, kind=core_rknd ), zt )
      call stat_update_var( iNPRC, lh_stat_sample_weight*real( NPRC, kind=core_rknd ), zt )
      call stat_update_var( iNPSACWI, lh_stat_sample_weight*real( NPSACWI, kind=core_rknd ), zt )
      call stat_update_var( iNPSACWG, lh_stat_sample_weight*real( NPSACWG, kind=core_rknd ), zt )
      call stat_update_var( iNPRAI, lh_stat_sample_weight*real( NPRAI, kind=core_rknd ), zt )
      call stat_update_var( iNMULTS, lh_stat_sample_weight*real( NMULTS, kind=core_rknd ), zt )
      call stat_update_var( iNMULTG, lh_stat_sample_weight*real( NMULTG, kind=core_rknd ), zt )
      call stat_update_var( iNMULTR, lh_stat_sample_weight*real( NMULTR, kind=core_rknd ), zt )
      call stat_update_var( iNMULTRG, lh_stat_sample_weight*real( NMULTRG, kind=core_rknd ), zt )
      call stat_update_var( iNNUCCD, lh_stat_sample_weight*real( NNUCCD, kind=core_rknd ), zt )
      call stat_update_var( iNSUBI, lh_stat_sample_weight*real( NSUBI, kind=core_rknd ), zt )
      call stat_update_var( iNGMLTG, lh_stat_sample_weight*real( NGMLTG, kind=core_rknd ), zt )
      call stat_update_var( iNSUBG, lh_stat_sample_weight*real( NSUBG, kind=core_rknd ), zt )
      call stat_update_var( iNACT, lh_stat_sample_weight*real( NACT,kind=core_rknd ), zt )
      call stat_update_var( iSIZEFIX_NR, lh_stat_sample_weight*real( SIZEFIX_NR,kind=core_rknd),zt)
      call stat_update_var( iSIZEFIX_NC, lh_stat_sample_weight*real( SIZEFIX_NC,kind=core_rknd),zt)
      call stat_update_var( iSIZEFIX_NI, lh_stat_sample_weight*real( SIZEFIX_NI,kind=core_rknd),zt)
      call stat_update_var( iSIZEFIX_NS, lh_stat_sample_weight*real( SIZEFIX_NS,kind=core_rknd),zt)
      call stat_update_var( iSIZEFIX_NG, lh_stat_sample_weight*real( SIZEFIX_NG,kind=core_rknd),zt)
      call stat_update_var( iNEGFIX_NR, lh_stat_sample_weight*real( NEGFIX_NR,kind=core_rknd ),zt)
      call stat_update_var( iNEGFIX_NC, lh_stat_sample_weight*real( NEGFIX_NC,kind=core_rknd ),zt)
      call stat_update_var( iNEGFIX_NI, lh_stat_sample_weight*real( NEGFIX_NI,kind=core_rknd ),zt)
      call stat_update_var( iNEGFIX_NS, lh_stat_sample_weight*real( NEGFIX_NS,kind=core_rknd ),zt)
      call stat_update_var( iNEGFIX_NG, lh_stat_sample_weight*real( NEGFIX_NG,kind=core_rknd ),zt)
      call stat_update_var( iNIM_MORR_CL, lh_stat_sample_weight &
                *real( NIM_MORR_CL,kind=core_rknd ), zt )
      call stat_update_var( iQC_INST, lh_stat_sample_weight*real( QC_INST,kind=core_rknd ),zt)
      call stat_update_var( iQR_INST, lh_stat_sample_weight*real( QR_INST,kind=core_rknd ),zt)
      call stat_update_var( iQI_INST, lh_stat_sample_weight*real( QI_INST,kind=core_rknd ),zt)
      call stat_update_var( iQS_INST, lh_stat_sample_weight*real( QS_INST,kind=core_rknd ),zt)
      call stat_update_var( iQG_INST, lh_stat_sample_weight*real( QG_INST,kind=core_rknd ),zt)
      call stat_update_var( iNC_INST, lh_stat_sample_weight*real( NC_INST,kind=core_rknd ),zt)
      call stat_update_var( iNR_INST, lh_stat_sample_weight*real( NR_INST,kind=core_rknd ),zt)
      call stat_update_var( iNI_INST, lh_stat_sample_weight*real( NI_INST,kind=core_rknd ),zt)
      call stat_update_var( iNS_INST, lh_stat_sample_weight*real( NS_INST,kind=core_rknd ),zt)
      call stat_update_var( iNG_INST, lh_stat_sample_weight*real( NG_INST,kind=core_rknd ),zt)

      call stat_update_var( iT_in_K_mc, lh_stat_sample_weight*real( T_in_K_mc, kind=core_rknd ),zt)

    end if ! l_stats_samp

      ! --- Number concentrations ---
      ! No budgets for sedimentation are output
    if ( .not. l_latin_hypercube .and. l_stats_samp ) then
      ! Effective radii of hydrometeor species
      call stat_update_var( ieff_rad_cloud, real( effc(:), kind = core_rknd ), zt )
      call stat_update_var( ieff_rad_ice, real( effi(:), kind = core_rknd ), zt )
      call stat_update_var( ieff_rad_snow, real( effs(:), kind = core_rknd ), zt )
      call stat_update_var( ieff_rad_rain, real( effr(:), kind = core_rknd ), zt )
      call stat_update_var( ieff_rad_graupel, real( effg(:), kind = core_rknd ), zt )

      ! Snow and Rain rates at the bottom of the domain, in mm/day
      call stat_update_var_pt( iprecip_rate_sfc, 1, &
        real(Morr_precip_rate, kind = core_rknd) * &
        real( sec_per_day, kind = core_rknd) / real( dt, kind = core_rknd ), sfc )

      call stat_update_var_pt( imorr_snow_rate, 1, &
        real(Morr_snow_rate, kind = core_rknd) * &
        real( sec_per_day, kind = core_rknd) / real( dt, kind = core_rknd ), sfc )

    end if ! .not. l_latin_hypercube .and l_stats_samp

    if ( l_latin_hypercube .and. l_stats_samp ) then
      ! Snow and Rain rates at the bottom of the domain, in mm/day
      call stat_update_var_pt( iprecip_rate_sfc, 1, &
        lh_stat_sample_weight*real( Morr_precip_rate, kind = core_rknd ) * &
        real( sec_per_day, kind = core_rknd) / real( dt, kind = core_rknd ), sfc )

      call stat_update_var_pt( iLH_morr_snow_rate, 1, &
        lh_stat_sample_weight*real( Morr_snow_rate, kind = core_rknd ) * &
        real( sec_per_day, kind = core_rknd) / real( dt, kind = core_rknd ), LH_sfc )

    end if ! l_latin_hypercube .and. l_stats_samp

    return
  end subroutine morrison_micro_driver

end module morrison_micro_driver_module
