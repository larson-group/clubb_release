! $Id$
module morrison_microphys_module

  implicit none

  public :: morrison_microphys_driver

  private

  contains
!-------------------------------------------------------------------------------
  subroutine morrison_microphys_driver( &
               gr, ngrdcol, dt, nzt, &
               hydromet_dim, hm_metadata, &
               l_latin_hypercube, thlm, wm_zt, p_in_Pa, &
               exner, rho, cloud_frac, w_std_dev, &
               dzq, rcm, Ncm, chi, rvm, hydromet, &
               saturation_formula, &
               sample_weight, stats, &
               hydromet_mc, hydromet_vel_zt, Ncm_mc, &
               rcm_mc, rvm_mc, thlm_mc, &
               rrm_auto_diag, rrm_accr_diag, rrm_evap_diag, &
               Nrm_auto_diag, Nrm_evap_diag )

    ! Description:
    ! Wrapper for the Morrison microphysics, applied to all columns.
    !
    ! References:
    ! None
    !-----------------------------------------------------------------------

    ! The version of the Morrison 2005 microphysics that is in SAM.
    use module_MP_graupel, only: &
        M2005MICRO_GRAUPEL  ! Procedure

    use constants_clubb, only: &
        Lv,   & ! Constants
        Ls,   &
        Cp,   &
        grav, &
        zero

    use T_in_K_module, only: &
        T_in_K2thlm_api, & ! Procedure(s)
        thlm2T_in_K_api

    use model_flags, only: &
        l_evaporate_cold_rcm  ! Flag(s)

    use parameters_microphys, only: &
        l_ice_microphys, & ! Flag(s)
        l_graupel, &
        lh_microphys_type, &
        lh_microphys_non_interactive

    use corr_varnce_module, only: &
        hm_metadata_type

    use constants_clubb, only: &
        sec_per_day

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api   ! Procedure

    use stats_netcdf, only: &
        stats_type, &
        stats_update

    use grid_class, only: &
        grid ! Type

    implicit none

    ! External
    intrinsic :: max, real, maxval

    real( kind = core_rknd ), parameter :: &
      w_thresh = 0.1_core_rknd ! Minimum value w for latin hypercube [m/s]

    ! Input Variables
    type(grid), intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in) :: dt ! Model timestep        [s]

    integer, intent(in) :: &
      ngrdcol,    & ! Number of model columns       [-]
      nzt,        & ! Points in the Vertical      [-]
      hydromet_dim  ! Number of hydrometeor species [-]

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    logical, intent(in) :: &
      l_latin_hypercube   ! Whether we're using latin hypercube sampling

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      thlm,       & ! Liquid potential temperature       [K]
      p_in_Pa,    & ! Pressure                           [Pa]
      exner,      & ! Exner function                     [-]
      rho,        & ! Density on thermo. grid            [kg/m^3]
      cloud_frac    ! Cloud fraction                     [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      wm_zt, &     ! Mean vertical velocity on the thermo grid     [m/s]
      w_std_dev, & ! Standard deviation of vertical vel. [m/s]
      dzq          ! Change in altitude                  [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      rcm,          & ! Cloud water mixing ratio                  [kg/kg]
      Ncm,          & ! Grid mean value for cloud droplet conc.    [#/kg]
      chi,          & ! The variable 's' from Mellor              [kg/kg]
      rvm             ! Vapor water mixing ratio                  [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      sample_weight ! SILHS sample weight; unity for an ordinary call

    type(stats_type), intent(inout) :: &
      stats

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(out) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel_zt ! Hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      Ncm_mc    ! Cloud droplet concentration time tendency    [num/kg/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      rrm_auto_diag, &
      rrm_accr_diag, &
      rrm_evap_diag, &
      Nrm_auto_diag, &
      Nrm_evap_diag

    ! Local Variables

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rrm_auto,     &  ! Autoconversion rate                 [kg/kg/s]
      rrm_accr,     &  ! Accretion rate                      [kg/kg/s]
      rrm_evap,     &  ! Rain evaporation rate               [kg/kg/s]
      T_in_K_core      ! Temperature at CLUBB precision      [K]

    real, dimension(ngrdcol,nzt) :: &
      effc, effi, effg, effs, effr ! Effective droplet radii [μ]

    real, dimension(ngrdcol,nzt) :: &
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
    real, dimension(ngrdcol,nzt) :: &
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

    real, dimension(ngrdcol,nzt) :: &
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
          
    real, dimension(ngrdcol,nzt) :: &
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

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: & ! Temporary variables
      rrm,    & ! Mean rain water mixing ratio            [kg/kg]
      rim,     & ! Mean ice mixing ratio                   [kg/kg]
      rsm,    & ! Mean snow mixing ratio                  [kg/kg]
      rgm    ! Mean graupel mixing ratio               [kg/kg]

    real, dimension(ngrdcol) :: Morr_snow_rate, Morr_precip_rate

    real, dimension(ngrdcol,nzt,hydromet_dim) :: &
      hydromet_r4,    & ! Temporary variable
      hydromet_mc_r4

    real, dimension(ngrdcol,nzt) :: & ! Temporary variables
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

    real, dimension(ngrdcol,nzt) :: &
      rcm_mc_r4,    &
      rvm_mc_r4,    &
      P_in_pa_r4,   &
      rho_r4,       &
      dzq_r4,       &
      wm_zt_r4,     & ! Mean vertical velocity on the thermo grid  [m/s]
      w_std_dev_r4    ! Standard deviation of w                    [m/s]

    integer :: i, k

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      Nrm_auto, & ! Change in Nrm due to autoconversion               [num/kg/s]
      Nrm_evap    ! Change in Nrm due to evaporation                  [num/kg/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol) :: rsm_sd_morr_int

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      hl_before, &
      qto_before, &
      hl_after, &
      hl_on_Cp_residual, &
      qto_after, &
      qto_residual

    ! Just to avoid typing hm_metadata%iixx everywhere
    integer ::  & 
      iirr, &
      iirs, &
      iiri, &
      iirg, &
      iiNr, &
      iiNs, &
      iiNi, &
      iiNg

    ! ---- Begin Code ----

    iirr = hm_metadata%iirr
    iirs = hm_metadata%iirs
    iiri = hm_metadata%iiri
    iirg = hm_metadata%iirg
    iiNr = hm_metadata%iiNr
    iiNs = hm_metadata%iiNs
    iiNi = hm_metadata%iiNi
    iiNg = hm_metadata%iiNg

    ! Determine temperature
    !$acc data copyin( thlm, exner, rcm ) copyout( T_in_K_core )
    T_in_K_core = thlm2T_in_K_api( nzt, ngrdcol, thlm, exner, rcm )
    !$acc end data
    T_in_K = real( T_in_K_core )

    if ( l_latin_hypercube ) then
      ! Don't use sgs cloud fraction to weight the tendencies
      cloud_frac_in = 0.0

      wm_zt_r4 = real( max( wm_zt, w_thresh ) ) ! Impose a minimum value on w
      w_std_dev_r4 = 0. ! Don't add in a standard deviation for aerosol activation

    else 
      ! Use sgs cloud fraction to weight tendencies
      cloud_frac_in = real( cloud_frac )

      wm_zt_r4 = real( wm_zt ) ! Use the mean value without a threshold
      w_std_dev_r4 = real( w_std_dev ) ! Add in a standard deviation
    end if

    rcm_r4 = real( rcm )
    rvm_r4 = real( rvm )

    ! Note: The Ncm_r4 variable is only used if INUM = 0 in the Morrison code;
    ! otherwise the NDCNST variable is used as the fixed value.
    Ncm_r4 = real( Ncm )
    
    hydromet_r4 = real( hydromet )

    if ( l_evaporate_cold_rcm ) then
       ! Convert liquid to vapor at temperatures colder than -37C
       do k = 1, nzt
          do i = 1, ngrdcol
             if ( T_in_K(i,k) < 236.15 ) then
                rcm_r4(i,k) = 0.0
                cloud_frac_in(i,k) = 0.0
                Ncm_r4(i,k) = 0.0
             endif
          enddo
       enddo
    end if
    
    ! Initialize tendencies to zero
    T_in_K_mc = 0.0
    rcm_mc = 0.0_core_rknd
    rvm_mc = 0.0_core_rknd
    hydromet_mc = 0.0_core_rknd
    Ncm_mc = 0.0_core_rknd

    ! Initialize process rates not guaranteed to be set on Morrison's early-exit path.
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
    NSAGG = 0.0
    NPRCI = 0.0
    NSCNG = 0.0
    NSUBS = 0.0
    NNUCCC = 0.0
    NPSACWS = 0.0
    NPSACWI = 0.0
    NPSACWG = 0.0
    NPRAI = 0.0
    NMULTS = 0.0
    NMULTG = 0.0
    NMULTR = 0.0
    NMULTRG = 0.0
    NNUCCD = 0.0
    NSUBI = 0.0
    NSUBG = 0.0

    hydromet_mc_r4 = real( hydromet_mc )
    Ncm_mc_r4 = real( Ncm_mc )
    rcm_mc_r4 = real( rcm_mc )
    rvm_mc_r4 = real( rvm_mc )
    P_in_pa_r4 = real( P_in_Pa )
    rho_r4 = real( rho )
    dzq_r4 = real( dzq )


    ! Unpack hydrometeor arrays.
    rrm = hydromet(:,:,iirr)

    rrm_r4 = hydromet_r4(:,:,iirr)
    Nrm_r4 = hydromet_r4(:,:,iiNr)

    rrm_mc_r4 = hydromet_mc_r4(:,:,iirr)
    Nrm_mc_r4 = hydromet_mc_r4(:,:,iiNr)

    if ( l_ice_microphys ) then

       rim = hydromet(:,:,iiri)
       rsm = hydromet(:,:,iirs)

       rim_r4 = hydromet_r4(:,:,iiri)
       Nim_r4 = hydromet_r4(:,:,iiNi)
       rsm_r4 = hydromet_r4(:,:,iirs)
       Nsm_r4 = hydromet_r4(:,:,iiNs)

       rim_mc_r4 = hydromet_mc_r4(:,:,iiri)
       Nim_mc_r4 = hydromet_mc_r4(:,:,iiNi)
       rsm_mc_r4 = hydromet_mc_r4(:,:,iirs)
       Nsm_mc_r4 = hydromet_mc_r4(:,:,iiNs)

       if ( l_graupel ) then

          rgm = hydromet(:,:,iirg)

          rgm_r4 = hydromet_r4(:,:,iirg)
          Ngm_r4 = hydromet_r4(:,:,iiNg)

          rgm_mc_r4 = hydromet_mc_r4(:,:,iirg)
          Ngm_mc_r4 = hydromet_mc_r4(:,:,iiNg)

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

    do k = 1, nzt
      do i = 1, ngrdcol
        hl_before(i,k) = Cp * real( T_in_K(i,k), kind = core_rknd ) &
                         + grav * gr%zt(i,k) &
                         - Lv * ( rcm(i,k) + rrm(i,k) ) &
                         - Ls * ( rim(i,k) + rsm(i,k) + rgm(i,k) )
      enddo
    enddo

    qto_before = rvm + rcm + rrm + rim + rsm + rgm

    ! Call the one-column Morrison microphysics core.
    do i = 1, ngrdcol
      ! Call the Morrison microphysics
      call M2005MICRO_GRAUPEL &
         ( rcm_mc_r4(i,:), rim_mc_r4(i,:), rsm_mc_r4(i,:), &
           rrm_mc_r4(i,:), Ncm_mc_r4(i,:), &
           Nim_mc_r4(i,:), Nsm_mc_r4(i,:), &
           Nrm_mc_r4(i,:), rcm_r4(i,:), rim_r4(i,:), &
           rsm_r4(i,:), rrm_r4(i,:), Ncm_r4(i,:), &
           Nim_r4(i,:), Nsm_r4(i,:), Nrm_r4(i,:), &
           T_in_K_mc(i,:), rvm_mc_r4(i,:), T_in_K(i,:), rvm_r4(i,:), &
           P_in_pa_r4(i,:), rho_r4(i,:), dzq_r4(i,:), &
           wm_zt_r4(i,:), w_std_dev_r4(i,:), morr_rain_vel_r4(i,:), &
           Morr_precip_rate(i), Morr_snow_rate(i), effc(i,:), effi(i,:), &
           effs(i,:), effr(i,:), real( dt ), &
           1,1, 1,1, 1,nzt, 1,1, 1,1, 1,nzt, &
           rgm_mc_r4(i,:), Ngm_mc_r4(i,:), &
           rgm_r4(i,:), Ngm_r4(i,:), effg(i,:), &
           rgm_sten(i,:), rrm_sten(i,:), &
           rim_sten(i,:), rsm_sten(i,:), &
           rcm_sten(i,:), &
           NGSTEN(i,:), NRSTEN(i,:), NISTEN(i,:), NSSTEN(i,:), NCSTEN(i,:), &
           cloud_frac_in(i,:), &
           PRC(i,:), PRA(i,:), &
           PSMLT(i,:), EVPMS(i,:), PRACS(i,:), EVPMG(i,:), &
           PRACG(i,:), PRE(i,:), PGMLT(i,:), MNUCCC(i,:), &
           PSACWS(i,:), PSACWI(i,:), QMULTS(i,:), QMULTG(i,:), &
           PSACWG(i,:), PGSACW(i,:), PRD(i,:), PRCI(i,:), &
           PRAI(i,:), QMULTR(i,:), QMULTRG(i,:), MNUCCD(i,:), &
           PRACI(i,:), PRACIS(i,:), EPRD(i,:), MNUCCR(i,:), &
           PIACR(i,:), PIACRS(i,:), PGRACS(i,:), PRDS(i,:), &
           EPRDS(i,:), PSACR(i,:), PRDG(i,:), EPRDG(i,:), &
           NPRC1(i,:), NRAGG(i,:), NPRACG(i,:), NSUBR(i,:), &
           NSMLTR(i,:), NGMLTR(i,:), NPRACS(i,:), NNUCCR(i,:), NIACR(i,:), &
           NIACRS(i,:), NGRACS(i,:), NSMLTS(i,:), NSAGG(i,:), NPRCI(i,:), &
           NSCNG(i,:), NSUBS(i,:), &
           PCC(i,:), NNUCCC(i,:), NPSACWS(i,:), NPRA(i,:), NPRC(i,:), &
           NPSACWI(i,:), NPSACWG(i,:), NPRAI(i,:), &
           NMULTS(i,:), NMULTG(i,:), NMULTR(i,:), NMULTRG(i,:), NNUCCD(i,:), &
           NSUBI(i,:), NGMLTG(i,:), NSUBG(i,:), NACT(i,:), &
           SIZEFIX_NR(i,:), SIZEFIX_NC(i,:), SIZEFIX_NI(i,:), &
           SIZEFIX_NS(i,:), SIZEFIX_NG(i,:), &
           NEGFIX_NR(i,:), NEGFIX_NC(i,:), NEGFIX_NI(i,:), &
           NEGFIX_NS(i,:), NEGFIX_NG(i,:), &
           NIM_MORR_CL(i,:), QC_INST(i,:), QR_INST(i,:), QI_INST(i,:), &
           QS_INST(i,:), QG_INST(i,:), NC_INST(i,:), NR_INST(i,:), &
           NI_INST(i,:), NS_INST(i,:), NG_INST(i,:) )

    enddo

    if ( clubb_at_least_debug_level_api( 2 ) ) then
      call print_morr_error_output()
    endif ! clubb_at_least_debug_level_api( 2 )

    do k = 1, nzt
      do i = 1, ngrdcol
        hl_after(i,k) = Cp * real( T_in_K(i,k), kind = core_rknd ) &
                        + grav * gr%zt(i,k) &
                        - Lv * ( real( rcm_r4(i,k), kind = core_rknd) &
                                 + real( rrm_r4(i,k), kind = core_rknd ) ) &
                        - Ls * ( real( rim_r4(i,k), kind = core_rknd ) &
                                 + real( rsm_r4(i,k), kind = core_rknd ) &
                                 + real( rgm_r4(i,k), kind = core_rknd ) )
      enddo
    enddo

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
    hydromet_r4(:,:,iirr) = rrm_r4
    hydromet_r4(:,:,iiNr) = Nrm_r4

    hydromet_mc_r4(:,:,iirr) = rrm_mc_r4
    hydromet_mc_r4(:,:,iiNr) = Nrm_mc_r4

    if ( l_ice_microphys ) then

       hydromet_r4(:,:,iiri) = rim_r4
       hydromet_r4(:,:,iiNi) = Nim_r4
       hydromet_r4(:,:,iirs) = rsm_r4
       hydromet_r4(:,:,iiNs) = Nsm_r4

       hydromet_mc_r4(:,:,iiri) = rim_mc_r4
       hydromet_mc_r4(:,:,iiNi) = Nim_mc_r4
       hydromet_mc_r4(:,:,iirs) = rsm_mc_r4
       hydromet_mc_r4(:,:,iiNs) = Nsm_mc_r4

       if ( l_graupel ) then

          hydromet_r4(:,:,iirg) = rgm_r4
          hydromet_r4(:,:,iiNg) = Ngm_r4

          hydromet_mc_r4(:,:,iirg) = rgm_mc_r4
          hydromet_mc_r4(:,:,iiNg) = Ngm_mc_r4

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
    hydromet_mc = ( real( hydromet_r4, kind = core_rknd ) - hydromet ) / dt

    Ncm_mc = ( real( Ncm_r4, kind = core_rknd ) - Ncm ) / dt

    ! Update thetal based on absolute temperature
    thlm_mc = ( T_in_K2thlm_api( real( T_in_K, kind = core_rknd ), exner, &
                real( rcm_r4, kind = core_rknd ) ) - thlm ) / dt

    ! Sedimentation is handled within the Morrison microphysics
    hydromet_vel_zt = 0.0_core_rknd

    ! Output rain sedimentation velocity
    ! Multiply by -1 so that negative is associated with falling precip
    morr_rain_vel_r4 = morr_rain_vel_r4 * (-1.0)
    hydromet_vel_zt(:,:,iirr) = real( morr_rain_vel_r4, kind = core_rknd )
    

    rsm_sd_morr_int = zero
    do k = 1, nzt
      do i = 1, ngrdcol
        rsm_sd_morr_int(i) = rsm_sd_morr_int(i) &
                             + rho(i,k) &
                               * real( rsm_sten(i,k), kind=core_rknd ) &
                               * gr%dzt(i,k)
      enddo
    enddo

    call update_microphys_stat_sfc( "rs_sd_morr_int", rsm_sd_morr_int )

    if ( clubb_at_least_debug_level_api( 1 ) .and. &
         any( rsm_sd_morr_int > &
              maxval( real( rsm_sten, kind=core_rknd ), dim=2 ) ) ) then
      print *, "Warning: rsm_sd_morr was not conservative!" // &
               " rsm_sd_morr_verical_integr = ", rsm_sd_morr_int
    endif
    
    call update_microphys_stat( "rrm_auto", rrm_auto )
    call update_microphys_stat( "rrm_accr", rrm_accr )
    call update_microphys_stat( "rrm_evap", rrm_evap )

    call update_microphys_stat( "Nrm_auto", Nrm_auto )
    call update_microphys_stat( "Nrm_evap", Nrm_evap )

    ! Update Morrison budgets
    call update_microphys_stat( "hl_on_Cp_residual", hl_on_Cp_residual )
    call update_microphys_stat( "qto_residual", qto_residual )
    call update_microphys_stat( "rgm_sd_morr", real( rgm_sten, kind = core_rknd ) )
    call update_microphys_stat( "rrm_sd_morr", real( rrm_sten, kind = core_rknd ) )
    call update_microphys_stat( "rsm_sd_morr", real( rsm_sten, kind = core_rknd ) )
    call update_microphys_stat( "rim_sd_mg_morr", real( rim_sten, kind = core_rknd ) )
    call update_microphys_stat( "rcm_sd_mg_morr", real( rcm_sten, kind = core_rknd) )
    call update_microphys_stat("PRC", real(PRC,kind=core_rknd) )
    call update_microphys_stat("PRA", real(PRA,kind=core_rknd) )
    call update_microphys_stat("PRE", real(PRE,kind=core_rknd) )
    call update_microphys_stat( "PSMLT", real( PSMLT, kind=core_rknd ) )
    call update_microphys_stat( "EVPMS", real( EVPMS, kind=core_rknd ) )
    call update_microphys_stat( "PRACS", real( PRACS, kind=core_rknd ) )
    call update_microphys_stat( "EVPMG", real( EVPMG, kind=core_rknd ) )
    call update_microphys_stat( "PRACG", real( PRACG, kind=core_rknd ) )
    call update_microphys_stat( "PGMLT", real( PGMLT, kind=core_rknd ) )
    call update_microphys_stat( "MNUCCC", real( MNUCCC, kind=core_rknd ) )
    call update_microphys_stat( "PSACWS", real( PSACWS, kind=core_rknd ) )
    call update_microphys_stat( "PSACWI", real( PSACWI, kind=core_rknd ) )
    call update_microphys_stat( "QMULTS", real( QMULTS, kind=core_rknd ) )
    call update_microphys_stat( "QMULTG", real( QMULTG, kind=core_rknd ) )
    call update_microphys_stat( "PSACWG", real( PSACWG, kind=core_rknd ) )
    call update_microphys_stat( "PGSACW", real( PGSACW, kind=core_rknd ) )
    call update_microphys_stat( "PRD", real( PRD, kind=core_rknd ) )
    call update_microphys_stat( "PRCI", real( PRCI, kind=core_rknd ) )
    call update_microphys_stat( "PRAI", real( PRAI, kind=core_rknd ) )
    call update_microphys_stat( "QMULTR", real( QMULTR, kind=core_rknd ) )
    call update_microphys_stat( "QMULTRG", real( QMULTRG, kind=core_rknd ) )
    call update_microphys_stat( "MNUCCD", real( MNUCCD, kind=core_rknd ) )
    call update_microphys_stat( "PRACI", real( PRACI, kind=core_rknd ) )
    call update_microphys_stat( "PRACIS", real( PRACIS, kind=core_rknd ) )
    call update_microphys_stat( "EPRD", real( EPRD, kind=core_rknd ) )
    call update_microphys_stat( "MNUCCR", real( MNUCCR, kind=core_rknd ) )
    call update_microphys_stat( "PIACR", real( PIACR, kind=core_rknd ) )
    call update_microphys_stat( "PIACRS", real( PIACRS, kind=core_rknd ) )
    call update_microphys_stat( "PGRACS", real( PGRACS, kind=core_rknd ) )
    call update_microphys_stat( "PRDS", real( PRDS, kind=core_rknd ) )
    call update_microphys_stat( "EPRDS", real( EPRDS, kind=core_rknd ) )
    call update_microphys_stat( "PSACR", real( PSACR, kind=core_rknd ) )
    call update_microphys_stat( "PRDG", real( PRDG, kind=core_rknd ) )
    call update_microphys_stat( "EPRDG", real( EPRDG, kind=core_rknd ) )
      
    ! Update more Morrison budgets
    call update_microphys_stat( "NGSTEN", real( NGSTEN, kind=core_rknd ) )
    call update_microphys_stat( "NRSTEN", real( NRSTEN, kind=core_rknd ) )
    call update_microphys_stat( "NISTEN", real( NISTEN, kind=core_rknd ) )
    call update_microphys_stat( "NSSTEN", real( NSSTEN, kind=core_rknd ) )
    call update_microphys_stat( "NCSTEN", real( NCSTEN, kind=core_rknd ) )
    call update_microphys_stat( "NPRC1", real( NPRC1, kind=core_rknd ) )
    call update_microphys_stat( "NRAGG", real( NRAGG, kind=core_rknd ) )
    call update_microphys_stat( "NPRACG", real( NPRACG, kind=core_rknd ) )
    call update_microphys_stat( "NSUBR", real( NSUBR, kind=core_rknd ) )
    call update_microphys_stat( "NSMLTR", real( NSMLTR, kind=core_rknd ) )
    call update_microphys_stat( "NGMLTR", real( NGMLTR, kind=core_rknd ) )
    call update_microphys_stat( "NPRACS", real( NPRACS, kind=core_rknd ) )
    call update_microphys_stat( "NNUCCR", real( NNUCCR, kind=core_rknd ) )
    call update_microphys_stat( "NIACR", real( NIACR, kind=core_rknd ) )
    call update_microphys_stat( "NIACRS", real( NIACRS, kind=core_rknd ) )
    call update_microphys_stat( "NGRACS", real( NGRACS, kind=core_rknd ) )
    call update_microphys_stat( "NSMLTS", real( NSMLTS, kind=core_rknd ) )
    call update_microphys_stat( "NSAGG", real( NSAGG, kind=core_rknd ) )
    call update_microphys_stat( "NPRCI", real( NPRCI, kind=core_rknd ) )
    call update_microphys_stat( "NSCNG", real( NSCNG, kind=core_rknd ) )
    call update_microphys_stat( "NSUBS", real( NSUBS, kind=core_rknd ) )

    call update_microphys_stat( "PCC", real( PCC, kind=core_rknd ) )
    call update_microphys_stat( "NNUCCC", real( NNUCCC, kind=core_rknd ) )
    call update_microphys_stat( "NPSACWS", real( NPSACWS, kind=core_rknd ) )
    call update_microphys_stat( "NPRA", real( NPRA, kind=core_rknd ) )
    call update_microphys_stat( "NPRC", real( NPRC, kind=core_rknd ) )
    call update_microphys_stat( "NPSACWI", real( NPSACWI, kind=core_rknd ) )
    call update_microphys_stat( "NPSACWG", real( NPSACWG, kind=core_rknd ) )
    call update_microphys_stat( "NPRAI", real( NPRAI, kind=core_rknd ) )
    call update_microphys_stat( "NMULTS", real( NMULTS, kind=core_rknd ) )
    call update_microphys_stat( "NMULTG", real( NMULTG, kind=core_rknd ) )
    call update_microphys_stat( "NMULTR", real( NMULTR, kind=core_rknd ) )
    call update_microphys_stat( "NMULTRG", real( NMULTRG, kind=core_rknd ) )
    call update_microphys_stat( "NNUCCD", real( NNUCCD, kind=core_rknd ) )
    call update_microphys_stat( "NSUBI", real( NSUBI, kind=core_rknd ) )
    call update_microphys_stat( "NGMLTG", real( NGMLTG, kind=core_rknd ) )
    call update_microphys_stat( "NSUBG", real( NSUBG, kind=core_rknd ) )
    call update_microphys_stat( "NACT", real( NACT,kind=core_rknd ) )
    call update_microphys_stat( "SIZEFIX_NR", real( SIZEFIX_NR,kind=core_rknd) )
    call update_microphys_stat( "SIZEFIX_NC", real( SIZEFIX_NC,kind=core_rknd) )
    call update_microphys_stat( "SIZEFIX_NI", real( SIZEFIX_NI,kind=core_rknd) )
    call update_microphys_stat( "SIZEFIX_NS", real( SIZEFIX_NS,kind=core_rknd) )
    call update_microphys_stat( "SIZEFIX_NG", real( SIZEFIX_NG,kind=core_rknd) )
    call update_microphys_stat( "NEGFIX_NR", real( NEGFIX_NR,kind=core_rknd ) )
    call update_microphys_stat( "NEGFIX_NC", real( NEGFIX_NC,kind=core_rknd ) )
    call update_microphys_stat( "NEGFIX_NI", real( NEGFIX_NI,kind=core_rknd ) )
    call update_microphys_stat( "NEGFIX_NS", real( NEGFIX_NS,kind=core_rknd ) )
    call update_microphys_stat( "NEGFIX_NG", real( NEGFIX_NG,kind=core_rknd ) )
    call update_microphys_stat( "NIM_MORR_CL", real( NIM_MORR_CL,kind=core_rknd ) )
    call update_microphys_stat( "QC_INST", real( QC_INST,kind=core_rknd ) )
    call update_microphys_stat( "QR_INST", real( QR_INST,kind=core_rknd ) )
    call update_microphys_stat( "QI_INST", real( QI_INST,kind=core_rknd ) )
    call update_microphys_stat( "QS_INST", real( QS_INST,kind=core_rknd ) )
    call update_microphys_stat( "QG_INST", real( QG_INST,kind=core_rknd ) )
    call update_microphys_stat( "NC_INST", real( NC_INST,kind=core_rknd ) )
    call update_microphys_stat( "NR_INST", real( NR_INST,kind=core_rknd ) )
    call update_microphys_stat( "NI_INST", real( NI_INST,kind=core_rknd ) )
    call update_microphys_stat( "NS_INST", real( NS_INST,kind=core_rknd ) )
    call update_microphys_stat( "NG_INST", real( NG_INST,kind=core_rknd ) )

    call update_microphys_stat( "T_in_K_mc", real( T_in_K_mc, kind=core_rknd ) )

    ! --- Number concentrations ---
    ! No budgets for sedimentation are output
    ! Effective radii of hydrometeor species
    call update_microphys_stat( "eff_rad_cloud", real( effc, kind = core_rknd ) )
    call update_microphys_stat( "eff_rad_ice", real( effi, kind = core_rknd ) )
    call update_microphys_stat( "eff_rad_snow", real( effs, kind = core_rknd ) )
    call update_microphys_stat( "eff_rad_rain", real( effr, kind = core_rknd ) )
    call update_microphys_stat( "eff_rad_graupel", real( effg, kind=core_rknd ) )

      ! Snow and Rain rates at the bottom of the domain, in mm/day
    call update_microphys_stat_sfc( "precip_rate_sfc", &
      real(Morr_precip_rate, kind = core_rknd) * sec_per_day / dt )

    call update_microphys_stat_sfc( "morr_snow_rate", &
      real(Morr_snow_rate, kind = core_rknd) * sec_per_day / dt )

    rrm_auto_diag = rrm_auto
    rrm_accr_diag = rrm_accr
    rrm_evap_diag = rrm_evap
    Nrm_auto_diag = real( Nrm_auto, kind=core_rknd )
    Nrm_evap_diag = real( Nrm_evap, kind=core_rknd )

    return

  contains

    subroutine update_microphys_stat( name, value )

      implicit none

      character(len=*), intent(in) :: name
      real(kind=core_rknd), dimension(:,:), intent(in) :: value

      character(len=64) :: output_name

      if ( .not. stats%l_sample ) then
        return
      endif

      output_name = name
      if ( l_latin_hypercube .and. &
           lh_microphys_type == lh_microphys_non_interactive ) then
        select case ( trim(name) )
        case ( "rrm_auto", "rrm_accr", "rrm_evap", "Nrm_auto", "Nrm_evap" )
          output_name = "lh_" // trim(name)
        case default
          return
        end select
      endif

      if ( l_latin_hypercube ) then
        call stats_update( trim(output_name), value * sample_weight, stats, &
                           sub_timestep_average=.true. )
      else
        call stats_update( trim(output_name), value, stats )
      endif

    end subroutine update_microphys_stat

    subroutine update_microphys_stat_sfc( name, value )

      implicit none

      character(len=*), intent(in) :: name
      real(kind=core_rknd), dimension(:), intent(in) :: value

      if ( .not. stats%l_sample .or. &
           ( l_latin_hypercube .and. &
             lh_microphys_type == lh_microphys_non_interactive ) ) then
        return
      endif

      if ( l_latin_hypercube ) then
        call stats_update( name, value * sample_weight(:,1), stats, &
                           sub_timestep_average=.true. )
      else
        call stats_update( name, value, stats )
      endif

    end subroutine update_microphys_stat_sfc

    !=============================================================================
    subroutine print_morr_error_output()

      ! Description:
      ! Prints all the relevant Morrison microphysics inputs and outputs when
      ! a non-finite value has been detected in an output tendency from Morrison.
      ! Read the wrapper's arrays through host association to retain the full
      ! diagnostic context without duplicating its argument list.
      !-----------------------------------------------------------------------

      use constants_clubb, only: &
        fstderr ! Standard error output unit [-]

      use, intrinsic :: ieee_arithmetic, only: &
        ieee_is_finite

      implicit none

      !--------------------------- Local Variables ---------------------------
      logical :: nonfinite_at_lev ! Whether this column/level has a non-finite tendency [-]
      integer :: i, k             ! Column and vertical level indices [-]

      !--------------------------- Begin Code ---------------------------

      ! Check for a non-finite value in any of the major microphysics tendency variables
      ! that are output from Morrison microphysics for use in the code.
      if ( any( .not. ieee_is_finite( real( rcm_mc_r4, kind=core_rknd ) ) ) &
           .or. any( .not. ieee_is_finite( real( rim_mc_r4, kind=core_rknd ) ) ) &
           .or. any( .not. ieee_is_finite( real( rsm_mc_r4, kind=core_rknd ) ) ) &
           .or. any( .not. ieee_is_finite( real( rrm_mc_r4, kind=core_rknd ) ) ) &
           .or. any( .not. ieee_is_finite( real( rgm_mc_r4, kind=core_rknd ) ) ) &
           .or. any( .not. ieee_is_finite( real( Ncm_mc_r4, kind=core_rknd ) ) ) &
           .or. any( .not. ieee_is_finite( real( Nim_mc_r4, kind=core_rknd ) ) ) &
           .or. any( .not. ieee_is_finite( real( Nsm_mc_r4, kind=core_rknd ) ) ) &
           .or. any( .not. ieee_is_finite( real( Nrm_mc_r4, kind=core_rknd ) ) ) &
           .or. any( .not. ieee_is_finite( real( Ngm_mc_r4, kind=core_rknd ) ) ) &
           .or. any( .not. ieee_is_finite( real( rvm_mc_r4, kind=core_rknd ) ) ) &
           .or. any( .not. ieee_is_finite( real( T_in_K_mc, kind=core_rknd ) ) ) ) then
        write(fstderr,*) "non-finite detected in a Morrison microphysics tendency"
        do k = 1, nzt, 1
          do i = 1, ngrdcol
            nonfinite_at_lev = .false.
            if ( .not. ieee_is_finite( real( rcm_mc_r4(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in rcm_mc_r4 at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( .not. ieee_is_finite( real( rim_mc_r4(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in rim_mc_r4 at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( .not. ieee_is_finite( real( rsm_mc_r4(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in rsm_mc_r4 at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( .not. ieee_is_finite( real( rrm_mc_r4(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in rrm_mc_r4 at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( .not. ieee_is_finite( real( rgm_mc_r4(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in rgm_mc_r4 at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( .not. ieee_is_finite( real( Ncm_mc_r4(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in Ncm_mc_r4 at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( .not. ieee_is_finite( real( Nim_mc_r4(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in Nim_mc_r4 at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( .not. ieee_is_finite( real( Nsm_mc_r4(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in Nsm_mc_r4 at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( .not. ieee_is_finite( real( Nrm_mc_r4(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in Nrm_mc_r4 at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( .not. ieee_is_finite( real( Ngm_mc_r4(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in Ngm_mc_r4 at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( .not. ieee_is_finite( real( rvm_mc_r4(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in rvm_mc_r4 at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( .not. ieee_is_finite( real( T_in_K_mc(i,k), kind=core_rknd ) ) ) then
               write(fstderr,*) "non-finite detected in T_in_K_mc at column = ", i, ", k = ", k, &
                                "altitude (m) = ", gr%zt(i,k)
               nonfinite_at_lev = .true.
            endif
            if ( nonfinite_at_lev ) then
               write(fstderr,*) "At column = ", i, ", k = ", k, "Altitude (m) = ", gr%zt(i,k)
               write(fstderr,*) "thlm (in) = ", thlm(i,k)
               write(fstderr,*) "rvm (in) = ", rvm(i,k)
               write(fstderr,*) "rcm (in) = ", rcm(i,k)
               write(fstderr,*) "rrm (in) = ", hydromet(i,k,hm_metadata%iirr)
               write(fstderr,*) "rim (in) = ", hydromet(i,k,hm_metadata%iiri)
               write(fstderr,*) "rsm (in) = ", hydromet(i,k,hm_metadata%iirs)
               write(fstderr,*) "rgm (in) = ", hydromet(i,k,hm_metadata%iirg)
               write(fstderr,*) "Ncm (in) = ", Ncm(i,k)
               write(fstderr,*) "Nrm (in) = ", hydromet(i,k,hm_metadata%iiNr)
               write(fstderr,*) "Nim (in) = ", hydromet(i,k,hm_metadata%iiNi)
               write(fstderr,*) "Nsm (in) = ", hydromet(i,k,hm_metadata%iiNs)
               write(fstderr,*) "Ngm (in) = ", hydromet(i,k,hm_metadata%iiNg)
               write(fstderr,*) "rgm_sten = ", rgm_sten(i,k)
               write(fstderr,*) "rrm_sten = ", rrm_sten(i,k)
               write(fstderr,*) "rim_sten = ", rim_sten(i,k)
               write(fstderr,*) "rsm_sten = ", rsm_sten(i,k)
               write(fstderr,*) "rcm_sten = ", rcm_sten(i,k)
               write(fstderr,*) "NGSTEN = ", NGSTEN(i,k)
               write(fstderr,*) "NRSTEN = ", NRSTEN(i,k)
               write(fstderr,*) "NISTEN = ", NISTEN(i,k)
               write(fstderr,*) "NSSTEN = ", NSSTEN(i,k)
               write(fstderr,*) "NCSTEN = ", NCSTEN(i,k)
               write(fstderr,*) "cloud_frac_in = ", cloud_frac_in(i,k)
               write(fstderr,*) "PRC = ", PRC(i,k)
               write(fstderr,*) "PRA = ", PRA(i,k)
               write(fstderr,*) "PSMLT = ", PSMLT(i,k)
               write(fstderr,*) "EVPMS = ", EVPMS(i,k)
               write(fstderr,*) "PRACS = ", PRACS(i,k)
               write(fstderr,*) "EVPMG = ", EVPMG(i,k)
               write(fstderr,*) "PRACG = ", PRACG(i,k)
               write(fstderr,*) "PRE = ", PRE(i,k)
               write(fstderr,*) "PGMLT = ", PGMLT(i,k)
               write(fstderr,*) "MNUCCC = ", MNUCCC(i,k)
               write(fstderr,*) "PSACWS = ", PSACWS(i,k)
               write(fstderr,*) "PSACWI = ", PSACWI(i,k)
               write(fstderr,*) "QMULTS = ", QMULTS(i,k)
               write(fstderr,*) "QMULTG = ", QMULTG(i,k)
               write(fstderr,*) "PSACWG = ", PSACWG(i,k)
               write(fstderr,*) "PGSACW = ", PGSACW(i,k)
               write(fstderr,*) "PRD = ", PRD(i,k)
               write(fstderr,*) "PRCI = ", PRCI(i,k)
               write(fstderr,*) "PRAI = ", PRAI(i,k)
               write(fstderr,*) "QMULTR = ", QMULTR(i,k)
               write(fstderr,*) "QMULTRG = ", QMULTRG(i,k)
               write(fstderr,*) "MNUCCD = ", MNUCCD(i,k)
               write(fstderr,*) "PRACI = ", PRACI(i,k)
               write(fstderr,*) "PRACIS = ", PRACIS(i,k)
               write(fstderr,*) "EPRD = ", EPRD(i,k)
               write(fstderr,*) "MNUCCR = ", MNUCCR(i,k)
               write(fstderr,*) "PIACR = ", PIACR(i,k)
               write(fstderr,*) "PIACRS = ", PIACRS(i,k)
               write(fstderr,*) "PGRACS = ", PGRACS(i,k)
               write(fstderr,*) "PRDS = ", PRDS(i,k)
               write(fstderr,*) "EPRDS = ", EPRDS(i,k)
               write(fstderr,*) "PSACR = ", PSACR(i,k)
               write(fstderr,*) "PRDG = ", PRDG(i,k)
               write(fstderr,*) "EPRDG = ", EPRDG(i,k)
               write(fstderr,*) "NPRC1 = ", NPRC1(i,k)
               write(fstderr,*) "NRAGG = ", NRAGG(i,k)
               write(fstderr,*) "NPRACG = ", NPRACG(i,k)
               write(fstderr,*) "NSUBR = ", NSUBR(i,k)
               write(fstderr,*) "NSMLTR = ", NSMLTR(i,k)
               write(fstderr,*) "NGMLTR = ", NGMLTR(i,k)
               write(fstderr,*) "NPRACS = ", NPRACS(i,k)
               write(fstderr,*) "NNUCCR = ", NNUCCR(i,k)
               write(fstderr,*) "NIACR = ", NIACR(i,k)
               write(fstderr,*) "NIACRS = ", NIACRS(i,k)
               write(fstderr,*) "NGRACS = ", NGRACS(i,k)
               write(fstderr,*) "NSMLTS = ", NSMLTS(i,k)
               write(fstderr,*) "NSAGG = ", NSAGG(i,k)
               write(fstderr,*) "NPRCI = ", NPRCI(i,k)
               write(fstderr,*) "NSCNG = ", NSCNG(i,k)
               write(fstderr,*) "NSUBS = ", NSUBS(i,k)
               write(fstderr,*) "PCC = ", PCC(i,k)
               write(fstderr,*) "NNUCCC = ", NNUCCC(i,k)
               write(fstderr,*) "NPSACWS = ", NPSACWS(i,k)
               write(fstderr,*) "NPRA = ", NPRA(i,k)
               write(fstderr,*) "NPRC = ", NPRC(i,k)
               write(fstderr,*) "NPSACWI = ", NPSACWI(i,k)
               write(fstderr,*) "NPSACWG = ", NPSACWG(i,k)
               write(fstderr,*) "NPRAI = ", NPRAI(i,k)
               write(fstderr,*) "NMULTS = ", NMULTS(i,k)
               write(fstderr,*) "NMULTG = ", NMULTG(i,k)
               write(fstderr,*) "NMULTR = ", NMULTR(i,k)
               write(fstderr,*) "NMULTRG = ", NMULTRG(i,k)
               write(fstderr,*) "NNUCCD = ", NNUCCD(i,k)
               write(fstderr,*) "NSUBI = ", NSUBI(i,k)
               write(fstderr,*) "NGMLTG = ", NGMLTG(i,k)
               write(fstderr,*) "NSUBG = ", NSUBG(i,k)
               write(fstderr,*) "NACT = ", NACT(i,k)
               write(fstderr,*) "SIZEFIX_NR = ", SIZEFIX_NR(i,k)
               write(fstderr,*) "SIZEFIX_NC = ", SIZEFIX_NC(i,k)
               write(fstderr,*) "SIZEFIX_NI = ", SIZEFIX_NI(i,k)
               write(fstderr,*) "SIZEFIX_NS = ", SIZEFIX_NS(i,k)
               write(fstderr,*) "SIZEFIX_NG = ", SIZEFIX_NG(i,k)
               write(fstderr,*) "NEGFIX_NR = ", NEGFIX_NR(i,k)
               write(fstderr,*) "NEGFIX_NC = ", NEGFIX_NC(i,k)
               write(fstderr,*) "NEGFIX_NI = ", NEGFIX_NI(i,k)
               write(fstderr,*) "NEGFIX_NS = ", NEGFIX_NS(i,k)
               write(fstderr,*) "NEGFIX_NG = ", NEGFIX_NG(i,k)
               write(fstderr,*) "NIM_MORR_CL = ", NIM_MORR_CL(i,k)
               write(fstderr,*) "QC_INST = ", QC_INST(i,k)
               write(fstderr,*) "QR_INST = ", QR_INST(i,k)
               write(fstderr,*) "QI_INST = ", QI_INST(i,k)
               write(fstderr,*) "QS_INST = ", QS_INST(i,k)
               write(fstderr,*) "QG_INST = ", QG_INST(i,k)
               write(fstderr,*) "NC_INST = ", NC_INST(i,k)
               write(fstderr,*) "NR_INST = ", NR_INST(i,k)
               write(fstderr,*) "NI_INST = ", NI_INST(i,k)
               write(fstderr,*) "NS_INST = ", NS_INST(i,k)
               write(fstderr,*) "NG_INST = ", NG_INST(i,k)
               write(fstderr,*) "---------------------------------------------"
            endif ! nonfinite_at_lev
          enddo ! i = 1, ngrdcol
        enddo ! k = 1, nzt, 1
      endif


      return

    end subroutine print_morr_error_output

  end subroutine morrison_microphys_driver

  !=============================================================================

end module morrison_microphys_module
