! $Id$
module morrison_micro_driver_module

  implicit none

  public :: morrison_micro_driver

  private

  contains
!-------------------------------------------------------------------------------
  subroutine morrison_micro_driver &
             ( dt, nz, l_stats_samp, &
               l_latin_hypercube, thlm, wm, p_in_Pa, &
               exner, rho, cloud_frac, pdf_params, w_std_dev, &
               dzq, rcm, Ncm, s_mellor, rvm, Ncm_in_cloud, hydromet, &
               hydromet_mc, hydromet_vel_zt, &
               rcm_mc, rvm_mc, thlm_mc, &
               wprtp_mc_tndcy, wpthlp_mc_tndcy, &
               rtp2_mc_tndcy, thlp2_mc_tndcy, rtpthlp_mc_tndcy, &
               rrainm_auto, rrainm_accr )

! Description:
!   Wrapper for the Morrison microphysics
!
! References:
!   None
!-------------------------------------------------------------------------------

    use pdf_parameter_module, only: pdf_parameter

    use parameters_model, only: hydromet_dim

    ! The version of the Morrison 2005 microphysics that is in SAM.
    use module_MP_graupel, only: &
      M2005MICRO_GRAUPEL  ! Procedure

    use module_MP_graupel, only: &
      cloud_frac_thresh ! Constant

    use advance_xp2_xpyp_module, only: &
      update_xp2_mc_tndcy

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
      imorr_rain_rate, & !Variables
      imorr_snow_rate, &
      iLH_morr_rain_rate, &
      iLH_morr_snow_rate

    use stats_variables, only: &
      irrainm_PSMLT, & ! Variable(s)
      iNicem_EVPMS, &
      irrainm_PRACS, &
      irgraupelm_EVPMG, &
      irrainm_PRACG, &
      irrainm_PGMLT, &
      ircm_MNUCCC, &
      ircm_PSACWS, &
      ircm_PSACWI, &
      ircm_QMULTS, &
      ircm_QMULTG, &
      ircm_PSACWG, &
      ircm_PGSACW, &
      iricem_PRD, &
      iricem_PRCI, &
      iricem_PRAI, &
      irrainm_QMULTR, &
      irrainm_QMULTRG, &
      iricem_MNUCCD, &
      iricem_PRACI, &
      iricem_PRACIS, &
      iricem_EPRD, &
      irrainm_MNUCCR, &
      irrainm_PIACR, &
      irrainm_PIACRS, &
      irrainm_PGRACS, &
      iNicem_PRDS, &
      iNicem_EPRDS, &
      iNicem_PSACR, &
      irgraupelm_PRDG, &
      irgraupelm_EPRDG

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
      zero, &
      zero_threshold

    use clubb_precision, only: &
      core_rknd, & ! Variable(s)
      time_precision

    use model_flags, only: &
      l_morr_xp2_mc_tndcy

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

    type(pdf_parameter), target, dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      wm, &        ! Mean vertial velocity               [m/s]
      w_std_dev, & ! Standard deviation of vertical vel. [m/s]
      dzq          ! Change in altitude                  [m]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm,          & ! Cloud water mixing ratio                  [kg/kg]
      Ncm,          & ! In cloud value for cloud droplet conc.    [#/kg]
      s_mellor,     & ! The variable 's' from Mellor              [kg/kg]
      rvm,          & ! Vapor water mixing ratio                  [kg/kg]
      Ncm_in_cloud    ! Constant cloud droplet conc. within cloud [#/kg] 

    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), &
    target, intent(inout) :: &
      hydromet_mc,  & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel_zt ! Hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      wprtp_mc_tndcy,   & ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc_tndcy,  & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtp2_mc_tndcy,    & ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc_tndcy,   & ! Microphysics tendency for <thl'^2>  [K^2/s]
      rtpthlp_mc_tndcy, & ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]
      rrainm_auto,      & ! Autoconversion rate     [kg/kg/s]
      rrainm_accr         ! Accretion rate         [kg/kg/s]

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
      cloud_frac_in, & ! Cloud fraction used as input for the Morrison scheme [-]
      rrainm_auto_r4,& ! Autoconversion rate     [kg/kg/s]
      rrainm_accr_r4   ! Accretion rate         [kg/kg/s]

    real, dimension(nz) :: &
      rrainm_PSMLT_r4, &       ! Change Q melting snow to rain [kg/kg/s]
      Nicem_EVPMS_r4, &        ! Change Q melting snow evaporating [kg/kg/s]
      rrainm_PRACS_r4, &       ! Change Q rain-snow collection [kg/kg/s]
      rgraupelm_EVPMG_r4, &    ! Change Q melting of graupel and evap [kg/kg/s]
      rrainm_PRACG_r4, &       ! Change Q collection rain by graupel [kg/kg/s]
      rrainm_PGMLT_r4, &       ! Change Q melting of graupel [kg/kg/s]
      rcm_MNUCCC_r4, &         ! Change Q contact freeze droplets [kg/kg/s]
      rcm_PSACWS_r4, &         ! Change Q droplet accretion by snow [kg/kg/s]
      rcm_PSACWI_r4, &         ! Change Q droplet accretion by cloud ice [kg/kg/s]
      rcm_QMULTS_r4, &         ! Change Q due to ice mult droplets/snow [kg/kg/s]
      rcm_QMULTG_r4, &         ! Change Q due to ice mult droplets/graupel [kg/kg/s]
      rcm_PSACWG_r4, &         ! Change Q collection droplets by graupel [kg/kg/s]
      rcm_PGSACW_r4, & ! Conversion Q to graupel due to collection droplets by snow [kg/kg/s]
      ricem_PRD_r4, &          ! Dep cloud ice [kg/kg/s]
      ricem_PRCI_r4, &         ! Change Q autoconversion cloud ice by snow [kg/kg/s]
      ricem_PRAI_r4, &         ! Change Q accretion cloud ice [kg/kg/s]
      rrainm_QMULTR_r4, &      ! Change Q due to ice rain/snow [kg/kg/s]
      rrainm_QMULTRG_r4, &     ! Change Q due to ice mult rain/graupel [kg/kg/s]
      ricem_MNUCCD_r4, &       ! Change Q freezing areosol [kg/kg/s]
      ricem_PRACI_r4, &        ! Change QI ice-rain collection [kg/kg/s]
      ricem_PRACIS_r4, &       ! Change QI ice-rain collision, added to snow [kg/kg/s]
      ricem_EPRD_r4, &         ! Sublimation cloud ice [kg/kg/s]
      rrainm_MNUCCR_r4, &      ! Change Q due to contact freeze rain [kg/kg/s]
      rrainm_PIACR_r4, &       ! Change QR ice-rain collection [kg/kg/s]
      rrainm_PIACRS_r4, &      ! Change QR ice-rain collision, added to snow [kg/kg/s]
      rrainm_PGRACS_r4, &      ! Conversion Q to graupel due to collection rain by snow [kg/kg/s]
      Nicem_PRDS_r4, &         ! Dep of snow [kg/kg/s]
      Nicem_EPRDS_r4, &        ! Sublimation of snow [kg/kg/s]
      Nicem_PSACR_r4, &        ! Conversion due to coll of snow by rain [kg/kg/s]
      rgraupelm_PRDG_r4, &     ! Dep of graupel [kg/kg/s]
      rgraupelm_EPRDG_r4       ! Sublimation of graupel [kg/kg/s]

    real( kind = core_rknd ), dimension(nz) :: &
      rrainm_PSMLT, &       ! See above for units 
      Nicem_EVPMS, &         
      rrainm_PRACS, &        
      rgraupelm_EVPMG, &     
      rrainm_PRACG, &        
      rrainm_PGMLT, &        
      rcm_MNUCCC, &          
      rcm_PSACWS, &          
      rcm_PSACWI, &          
      rcm_QMULTS, &          
      rcm_QMULTG, &          
      rcm_PSACWG, &           
      rcm_PGSACW, &          
      ricem_PRD, &           
      ricem_PRCI, &          
      ricem_PRAI, &          
      rrainm_QMULTR, &       
      rrainm_QMULTRG, &      
      ricem_MNUCCD, &        
      ricem_PRACI, &         
      ricem_PRACIS, &        
      ricem_EPRD, &          
      rrainm_MNUCCR, &       
      rrainm_PIACR, &        
      rrainm_PIACRS, &       
      rrainm_PGRACS, &       
      Nicem_PRDS, &          
      Nicem_EPRDS, &         
      Nicem_PSACR, &         
      rgraupelm_PRDG, &      
      rgraupelm_EPRDG       

    real( kind = core_rknd ), dimension(nz) :: & 
      rcm_in_cloud     ! Liquid water in cloud           [kg/kg]

    real( kind = core_rknd ), pointer, dimension(:,:) :: &
      dummy

    real( kind = core_rknd ), pointer, dimension(:) :: &
      dummy_1D

    real :: Morr_snow_rate, Morr_rain_rate

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
      wm_r4, &     ! Mean vertical velocity   [m/s]
      w_std_dev_r4 ! Standard deviation of w  [m/s]

    integer :: i, k

    !variables needed to computer rtp2_mc_tndcy when l_morr_xp2_mc_tndcy = .true.
    real( kind = core_rknd ), dimension(nz) :: &
      rrainm_evap         !Evaporation of rain   [kg/kg/s]

    real, dimension(nz) :: &
      rrainm_evap_r4

    ! ---- Begin Code ----

    ! Some dummy assignments to make compiler warnings go away...
    if ( .false. ) then
      dummy => hydromet
      dummy => hydromet_vel_zt
      dummy => hydromet_mc
      dummy_1D => pdf_params(:)%cloud_frac1
      rcm_in_cloud = dummy(:,1)
      rcm_in_cloud = s_mellor
      rcm_in_cloud = dummy_1D
      rcm_in_cloud = Ncm_in_cloud
    end if

    ! Determine temperature
    T_in_K = real( thlm2T_in_K( thlm, exner, rcm ) )

    if ( l_latin_hypercube ) then
      ! Don't use sgs cloud fraction to weight the tendencies
      cloud_frac_in(1:nz) = 0.0

      wm_r4 = real( max( wm, w_thresh ) ) ! Impose a minimum value on w
      w_std_dev_r4 = 0. ! Don't add in a standard deviation for aerosol activation

    else 
      ! Use sgs cloud fraction to weight tendencies
      cloud_frac_in(1:nz) = real( cloud_frac(1:nz) )

      wm_r4 = real( wm ) ! Use the mean value without a threshold
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

    ! Initialize autoconversion/accretion to zero
    rrainm_auto_r4 = 0.0
    rrainm_accr_r4 = 0.0
    rrainm_evap_r4 = 0.0

    ! Initialize Morrison budgets to zero
    rrainm_PSMLT_r4 = 0.0
    Nicem_EVPMS_r4 = 0.0 
    rrainm_PRACS_r4 = 0.0
    rgraupelm_EVPMG_r4 = 0.0
    rrainm_PRACG_r4 = 0.0
    rrainm_PGMLT_r4 = 0.0
    rcm_MNUCCC_r4 = 0.0
    rcm_PSACWS_r4 = 0.0
    rcm_PSACWI_r4 = 0.0
    rcm_QMULTS_r4 = 0.0
    rcm_QMULTG_r4 = 0.0
    rcm_PSACWG_r4 = 0.0 
    rcm_PGSACW_r4 = 0.0
    ricem_PRD_r4 = 0.0
    ricem_PRCI_r4 = 0.0
    ricem_PRAI_r4 = 0.0
    rrainm_QMULTR_r4 = 0.0
    rrainm_QMULTRG_r4 = 0.0
    ricem_MNUCCD_r4 = 0.0
    ricem_PRACI_r4 = 0.0
    ricem_PRACIS_r4 = 0.0
    ricem_EPRD_r4 = 0.0
    rrainm_MNUCCR_r4 = 0.0
    rrainm_PIACR_r4 = 0.0
    rrainm_PIACRS_r4 = 0.0
    rrainm_PGRACS_r4 = 0.0
    Nicem_PRDS_r4 = 0.0
    Nicem_EPRDS_r4 = 0.0
    Nicem_PSACR_r4 = 0.0
    rgraupelm_PRDG_r4 = 0.0
    rgraupelm_EPRDG_r4 = 0.0

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
           wm_r4, w_std_dev_r4, morr_rain_vel_r4, &
           Morr_rain_rate, Morr_snow_rate, effc, effi, effs, effr, real( dt ), &
           1,1, 1,1, 1,nz, 1,1, 1,1, 2,nz, &
           hydromet_mc_r4(:,iirgraupelm), hydromet_mc_r4(:,iiNgraupelm), &
           hydromet_r4(:,iirgraupelm), hydromet_r4(:,iiNgraupelm), effg, &
           hydromet_sten(:,iirgraupelm), hydromet_sten(:,iirrainm), &
           hydromet_sten(:,iiricem), hydromet_sten(:,iirsnowm), &
           rcm_sten, cloud_frac_in, &
           rrainm_auto_r4, rrainm_accr_r4, &
           rrainm_PSMLT_r4, Nicem_EVPMS_r4, rrainm_PRACS_r4, rgraupelm_EVPMG_r4, &
           rrainm_PRACG_r4, rrainm_evap_r4,rrainm_PGMLT_r4,rcm_MNUCCC_r4, &
           rcm_PSACWS_r4, rcm_PSACWI_r4, rcm_QMULTS_r4, rcm_QMULTG_r4, &
           rcm_PSACWG_r4, rcm_PGSACW_r4, ricem_PRD_r4, ricem_PRCI_r4, &
           ricem_PRAI_r4, rrainm_QMULTR_r4, rrainm_QMULTRG_r4, ricem_MNUCCD_r4, &
           ricem_PRACI_r4, ricem_PRACIS_r4, ricem_EPRD_r4, rrainm_MNUCCR_r4, &
           rrainm_PIACR_r4, rrainm_PIACRS_r4, rrainm_PGRACS_r4, Nicem_PRDS_r4, &
           Nicem_EPRDS_r4, Nicem_PSACR_r4, rgraupelm_PRDG_r4, rgraupelm_EPRDG_r4 )
           

    !hydromet_mc = real( hydromet_mc_r4, kind = core_rknd )
    rcm_mc = real( rcm_mc_r4, kind = core_rknd )
    rvm_mc = real( rvm_mc_r4, kind = core_rknd )

    rrainm_auto = real( rrainm_auto_r4, kind = core_rknd )
    rrainm_accr = real( rrainm_accr_r4, kind = core_rknd ) 
    rrainm_evap = real( rrainm_evap_r4, kind = core_rknd )

    rrainm_PSMLT = real( rrainm_PSMLT_r4, kind = core_rknd )
    Nicem_EVPMS = real( Nicem_EVPMS_r4, kind = core_rknd ) 
    rrainm_PRACS = real( rrainm_PRACS_r4, kind = core_rknd )
    rgraupelm_EVPMG = real( rgraupelm_EVPMG_r4, kind = core_rknd )
    rrainm_PRACG = real( rrainm_PRACG_r4, kind = core_rknd )
    rrainm_PGMLT = real( rrainm_PGMLT_r4, kind = core_rknd )
    rcm_MNUCCC = real( rcm_MNUCCC_r4, kind = core_rknd )
    rcm_PSACWS = real( rcm_PSACWS_r4, kind = core_rknd )
    rcm_PSACWI = real( rcm_PSACWI_r4, kind = core_rknd )
    rcm_QMULTS = real( rcm_QMULTS_r4, kind = core_rknd )
    rcm_QMULTG = real( rcm_QMULTG_r4, kind = core_rknd )
    rcm_PSACWG = real( rcm_PSACWG_r4, kind = core_rknd ) 
    rcm_PGSACW = real( rcm_PGSACW_r4, kind = core_rknd )
    ricem_PRD = real( ricem_PRD_r4, kind = core_rknd )
    ricem_PRCI = real( ricem_PRCI_r4, kind = core_rknd )
    ricem_PRAI = real( ricem_PRAI_r4, kind = core_rknd )
    rrainm_QMULTR = real( rrainm_QMULTR_r4, kind = core_rknd )
    rrainm_QMULTRG = real( rrainm_QMULTRG_r4, kind = core_rknd )
    ricem_MNUCCD = real( ricem_MNUCCD_r4, kind = core_rknd )
    ricem_PRACI = real( ricem_PRACI_r4, kind = core_rknd )
    ricem_PRACIS = real( ricem_PRACIS_r4, kind = core_rknd )
    ricem_EPRD = real( ricem_EPRD_r4, kind = core_rknd )
    rrainm_MNUCCR = real( rrainm_MNUCCR_r4, kind = core_rknd )
    rrainm_PIACR = real( rrainm_PIACR_r4, kind = core_rknd )
    rrainm_PIACRS = real( rrainm_PIACRS_r4, kind = core_rknd )
    rrainm_PGRACS = real( rrainm_PGRACS_r4, kind = core_rknd )
    Nicem_PRDS = real( Nicem_PRDS_r4, kind = core_rknd )
    Nicem_EPRDS = real( Nicem_EPRDS_r4, kind = core_rknd )
    Nicem_PSACR = real( Nicem_PSACR_r4, kind = core_rknd )
    rgraupelm_PRDG = real( rgraupelm_PRDG_r4, kind = core_rknd )
    rgraupelm_EPRDG = real( rgraupelm_EPRDG_r4, kind = core_rknd )


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

    ! Set microphysics tendencies for model variances and covariances to 0.
    wprtp_mc_tndcy   = zero
    wpthlp_mc_tndcy  = zero
    rtp2_mc_tndcy    = zero
    thlp2_mc_tndcy   = zero
    rtpthlp_mc_tndcy = zero

    if ( l_morr_xp2_mc_tndcy ) then
      call update_xp2_mc_tndcy( nz, dt, cloud_frac, rcm, rvm, thlm, &
                                exner, rrainm_evap, pdf_params,     &
                                rtp2_mc_tndcy, thlp2_mc_tndcy       )

    end if

    if ( .not. l_latin_hypercube .and. l_stats_samp ) then

      ! -------- Sedimentation tendency from Morrison microphysics --------

      ! --- Mixing ratios ---

      call stat_update_var( irgraupelm_sd_morr, &
                 real( hydromet_sten(:,iirgraupelm), kind = core_rknd ), zt )

      call stat_update_var( irrainm_sd_morr, &
                 real( hydromet_sten(:,iirrainm), kind = core_rknd ), zt )

      call stat_update_var( irsnowm_sd_morr, &
                 real( hydromet_sten(:,iirsnowm), kind = core_rknd ), zt )

      call stat_update_var( iricem_sd_mg_morr, &
                 real( hydromet_sten(:,iiricem), kind = core_rknd ), zt )

      call stat_update_var( ircm_sd_mg_morr, &
                 real( rcm_sten, kind = core_rknd), zt )

      where ( cloud_frac(:) > real( cloud_frac_thresh, kind = core_rknd ) ) 
        rcm_in_cloud(:) = rcm / cloud_frac
      else where
        rcm_in_cloud(:) = rcm
      end where

      call stat_update_var( ircm_in_cloud, rcm_in_cloud, zt )

      call stat_update_var( irrainm_auto, real( rrainm_auto, kind=core_rknd ), zt )
      call stat_update_var( irrainm_accr, real( rrainm_accr, kind=core_rknd ), zt )
      call stat_update_var( irrainm_cond, real( rrainm_evap, kind=core_rknd ), zt )

      call stat_update_var( irrainm_PSMLT, real( rrainm_PSMLT, kind=core_rknd ), zt )
      call stat_update_var( iNicem_EVPMS, real( Nicem_EVPMS, kind=core_rknd ), zt )
      call stat_update_var( irrainm_PRACS, real( rrainm_PRACS, kind=core_rknd ), zt )
      call stat_update_var( irgraupelm_EVPMG, real( rgraupelm_EVPMG, kind=core_rknd ), zt )
      call stat_update_var( irrainm_PRACG, real( rrainm_PRACG, kind=core_rknd ), zt )
      call stat_update_var( irrainm_PGMLT, real( rrainm_PGMLT, kind=core_rknd ), zt )
      call stat_update_var( ircm_MNUCCC, real( rcm_MNUCCC, kind=core_rknd ), zt )
      call stat_update_var( ircm_PSACWS, real( rcm_PSACWS, kind=core_rknd ), zt )
      call stat_update_var( ircm_PSACWI, real( rcm_PSACWI, kind=core_rknd ), zt )
      call stat_update_var( ircm_QMULTS, real( rcm_QMULTS, kind=core_rknd ), zt )
      call stat_update_var( ircm_QMULTG, real( rcm_QMULTG, kind=core_rknd ), zt )
      call stat_update_var( ircm_PSACWG, real( rcm_PSACWG, kind=core_rknd ), zt )
      call stat_update_var( ircm_PGSACW, real( rcm_PGSACW, kind=core_rknd ), zt )
      call stat_update_var( iricem_PRD, real( ricem_PRD, kind=core_rknd ), zt )
      call stat_update_var( iricem_PRCI, real( ricem_PRCI, kind=core_rknd ), zt )
      call stat_update_var( iricem_PRAI, real( ricem_PRAI, kind=core_rknd ), zt )
      call stat_update_var( irrainm_QMULTR, real( rrainm_QMULTR, kind=core_rknd ), zt )
      call stat_update_var( irrainm_QMULTRG, real( rrainm_QMULTRG, kind=core_rknd ), zt )
      call stat_update_var( iricem_MNUCCD, real( ricem_MNUCCD, kind=core_rknd ), zt )
      call stat_update_var( iricem_PRACI, real( ricem_PRACI, kind=core_rknd ), zt )
      call stat_update_var( iricem_PRACIS, real( ricem_PRACIS, kind=core_rknd ), zt )
      call stat_update_var( iricem_EPRD, real( ricem_EPRD, kind=core_rknd ), zt )
      call stat_update_var( irrainm_MNUCCR, real( rrainm_MNUCCR, kind=core_rknd ), zt )
      call stat_update_var( irrainm_PIACR, real( rrainm_PIACR, kind=core_rknd ), zt )
      call stat_update_var( irrainm_PIACRS, real( rrainm_PIACRS, kind=core_rknd ), zt )
      call stat_update_var( irrainm_PGRACS, real( rrainm_PGRACS, kind=core_rknd ), zt )
      call stat_update_var( iNicem_PRDS, real( Nicem_PRDS, kind=core_rknd ), zt )
      call stat_update_var( iNicem_EPRDS, real( Nicem_EPRDS, kind=core_rknd ), zt )
      call stat_update_var( iNicem_PSACR, real( Nicem_PSACR, kind=core_rknd ), zt )
      call stat_update_var( irgraupelm_PRDG, real( rgraupelm_PRDG, kind=core_rknd ), zt )
      call stat_update_var( irgraupelm_EPRDG, real( rgraupelm_EPRDG, kind=core_rknd ), zt )

      ! --- Number concentrations ---
      ! No budgets for sedimentation are output

      ! Effective radii of hydrometeor species
      call stat_update_var( ieff_rad_cloud, real( effc(:), kind = core_rknd ), zt )
      call stat_update_var( ieff_rad_ice, real( effi(:), kind = core_rknd ), zt )
      call stat_update_var( ieff_rad_snow, real( effs(:), kind = core_rknd ), zt )
      call stat_update_var( ieff_rad_rain, real( effr(:), kind = core_rknd ), zt )
      call stat_update_var( ieff_rad_graupel, real( effg(:), kind = core_rknd ), zt )

      ! Snow and Rain rates at the bottom of the domain, in mm/day
      call stat_update_var_pt( imorr_rain_rate, 1, &
        real(Morr_rain_rate, kind = core_rknd) * &
        real( sec_per_day, kind = core_rknd) / real( dt, kind = core_rknd ), sfc )

      call stat_update_var_pt( imorr_snow_rate, 1, &
        real(Morr_snow_rate, kind = core_rknd) * &
        real( sec_per_day, kind = core_rknd) / real( dt, kind = core_rknd ), sfc )

    end if ! .not. l_latin_hypercube .and l_stats_samp

    if ( l_stats_samp ) then
      ! Snow and Rain rates at the bottom of the domain, in mm/day
      call stat_update_var_pt( iLH_morr_rain_rate, 1, &
        real( Morr_rain_rate, kind = core_rknd ) * &
        real( sec_per_day, kind = core_rknd) / real( dt, kind = core_rknd ), LH_sfc )

      call stat_update_var_pt( iLH_morr_snow_rate, 1, &
        real( Morr_snow_rate, kind = core_rknd ) * &
        real( sec_per_day, kind = core_rknd) / real( dt, kind = core_rknd ), LH_sfc )

    end if ! l_stats_samp

    return

    return
  end subroutine morrison_micro_driver

end module morrison_micro_driver_module
