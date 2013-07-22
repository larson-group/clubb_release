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
               dzq, rcm, Ncm, s_mellor, rvm, hydromet, &
               hydromet_mc, hydromet_vel_zt, &
               rcm_mc, rvm_mc, thlm_mc, &
               rtp2_mc_tndcy, thlp2_mc_tndcy, &
               rrainm_auto, rrainm_accr, rrainm_evap )

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
      iNGRACS

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
      rvm             ! Vapor water mixing ratio                  [kg/kg]

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
      rtp2_mc_tndcy,  & ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc_tndcy, & ! Microphysics tendency for <thl'^2>  [K^2/s]
      rrainm_auto,    & ! Autoconversion rate                 [kg/kg/s]
      rrainm_accr       ! Accretion rate                      [kg/kg/s]

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

    ! In the comments below, by "adds to" we mean that if the quantity is
    ! positive, it adds positively to the prognostic variable, but if the
    ! quantity is negative, it subtracts from the prognostic variable.
    real, dimension(nz) :: &
      PSMLT,  & ! Freezing of rain to form snow.
                !    Adds to rsnowm, subtracts from rrainm [kg/kg/s]
      EVPMS,  & ! Depositional growth of snow.
                !    Adds to rsnowm, subtracts from rvm [kg/kg/s]
      PRACS,  & ! Collection of rain by snow.
                !    Adds to rsnowm, subtracts from rrainm [kg/kg/s]
      EVPMG,  & ! Depositional growth of graupel.
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
      NGRACS    ! Collection of rain by snow.
                !    Adds to Ngraupelm, subtracts from Nrm and Nsnowm [#/kg/s]

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
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
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
           rcm_sten, &
           NGSTEN, NRSTEN, NISTEN, NSSTEN, NCSTEN, &
           cloud_frac_in, &
           rrainm_auto_r4, rrainm_accr_r4, &
           PSMLT, EVPMS, PRACS, EVPMG, &
           PRACG, rrainm_evap_r4,PGMLT,MNUCCC, &
           PSACWS, PSACWI, QMULTS, QMULTG, &
           PSACWG, PGSACW, PRD, PRCI, &
           PRAI, QMULTR, QMULTRG, MNUCCD, &
           PRACI, PRACIS, EPRD, MNUCCR, &
           PIACR, PIACRS, PGRACS, PRDS, &
           EPRDS, PSACR, PRDG, EPRDG, &
           NPRC1, NRAGG, NPRACG, NSUBR, NSMLTR, NGMLTR, NPRACS, NNUCCR, NIACR, &
           NIACRS, NGRACS )
           

    !hydromet_mc = real( hydromet_mc_r4, kind = core_rknd )
    rcm_mc = real( rcm_mc_r4, kind = core_rknd )
    rvm_mc = real( rvm_mc_r4, kind = core_rknd )

    rrainm_auto = real( rrainm_auto_r4, kind = core_rknd )
    rrainm_accr = real( rrainm_accr_r4, kind = core_rknd )
    rrainm_evap = real( rrainm_evap_r4, kind = core_rknd )
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

    if ( l_morr_xp2_mc_tndcy ) then

       call update_xp2_mc_tndcy( nz, dt, cloud_frac, rcm, rvm, thlm, &
                                 exner, rrainm_evap, pdf_params,     &
                                 rtp2_mc_tndcy, thlp2_mc_tndcy       )

    else

       ! Set microphysics tendencies for model variances to 0.
       rtp2_mc_tndcy  = zero
       thlp2_mc_tndcy = zero

    endif

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
      call stat_update_var( irrainm_cond, rrainm_evap, zt )

      call stat_update_var( iPSMLT, real( PSMLT, kind=core_rknd ), zt )
      call stat_update_var( iEVPMS, real( EVPMS, kind=core_rknd ), zt )
      call stat_update_var( iPRACS, real( PRACS, kind=core_rknd ), zt )
      call stat_update_var( iEVPMG, real( EVPMG, kind=core_rknd ), zt )
      call stat_update_var( iPRACG, real( PRACG, kind=core_rknd ), zt )
      call stat_update_var( iPGMLT, real( PGMLT, kind=core_rknd ), zt )
      call stat_update_var( iMNUCCC, real( MNUCCC, kind=core_rknd ), zt )
      call stat_update_var( iPSACWS, real( PSACWS, kind=core_rknd ), zt )
      call stat_update_var( iPSACWI, real( PSACWI, kind=core_rknd ), zt )
      call stat_update_var( iQMULTS, real( QMULTS, kind=core_rknd ), zt )
      call stat_update_var( iQMULTG, real( QMULTG, kind=core_rknd ), zt )
      call stat_update_var( iPSACWG, real( PSACWG, kind=core_rknd ), zt )
      call stat_update_var( iPGSACW, real( PGSACW, kind=core_rknd ), zt )
      call stat_update_var( iPRD, real( PRD, kind=core_rknd ), zt )
      call stat_update_var( iPRCI, real( PRCI, kind=core_rknd ), zt )
      call stat_update_var( iPRAI, real( PRAI, kind=core_rknd ), zt )
      call stat_update_var( iQMULTR, real( QMULTR, kind=core_rknd ), zt )
      call stat_update_var( iQMULTRG, real( QMULTRG, kind=core_rknd ), zt )
      call stat_update_var( iMNUCCD, real( MNUCCD, kind=core_rknd ), zt )
      call stat_update_var( iPRACI, real( PRACI, kind=core_rknd ), zt )
      call stat_update_var( iPRACIS, real( PRACIS, kind=core_rknd ), zt )
      call stat_update_var( iEPRD, real( EPRD, kind=core_rknd ), zt )
      call stat_update_var( iMNUCCR, real( MNUCCR, kind=core_rknd ), zt )
      call stat_update_var( iPIACR, real( PIACR, kind=core_rknd ), zt )
      call stat_update_var( iPIACRS, real( PIACRS, kind=core_rknd ), zt )
      call stat_update_var( iPGRACS, real( PGRACS, kind=core_rknd ), zt )
      call stat_update_var( iPRDS, real( PRDS, kind=core_rknd ), zt )
      call stat_update_var( iEPRDS, real( EPRDS, kind=core_rknd ), zt )
      call stat_update_var( iPSACR, real( PSACR, kind=core_rknd ), zt )
      call stat_update_var( iPRDG, real( PRDG, kind=core_rknd ), zt )
      call stat_update_var( iEPRDG, real( EPRDG, kind=core_rknd ), zt )
      
      ! Update more Morrison budgets
      call stat_update_var( iNGSTEN, real( NGSTEN, kind=core_rknd ), zt )
      call stat_update_var( iNRSTEN, real( NRSTEN, kind=core_rknd ), zt )
      call stat_update_var( iNISTEN, real( NISTEN, kind=core_rknd ), zt )
      call stat_update_var( iNSSTEN, real( NSSTEN, kind=core_rknd ), zt )
      call stat_update_var( iNCSTEN, real( NCSTEN, kind=core_rknd ), zt )
      call stat_update_var( iNPRC1, real( NPRC1, kind=core_rknd ), zt )
      call stat_update_var( iNRAGG, real( NRAGG, kind=core_rknd ), zt )
      call stat_update_var( iNPRACG, real( NPRACG, kind=core_rknd ), zt )
      call stat_update_var( iNSUBR, real( NSUBR, kind=core_rknd ), zt )
      call stat_update_var( iNSMLTR, real( NSMLTR, kind=core_rknd ), zt )
      call stat_update_var( iNGMLTR, real( NGMLTR, kind=core_rknd ), zt )
      call stat_update_var( iNPRACS, real( NPRACS, kind=core_rknd ), zt )
      call stat_update_var( iNNUCCR, real( NNUCCR, kind=core_rknd ), zt )
      call stat_update_var( iNIACR, real( NIACR, kind=core_rknd ), zt )
      call stat_update_var( iNIACRS, real( NIACRS, kind=core_rknd ), zt )
      call stat_update_var( iNGRACS, real( NGRACS, kind=core_rknd ), zt )

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
