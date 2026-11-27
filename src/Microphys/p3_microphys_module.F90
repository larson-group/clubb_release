!-------------------------------------------------------------------------------
! p3_microphys_module: connect the P3 microphysics scheme to CLUBB
!
! P3 (Predicted Particle Properties) is a bulk microphysics scheme with
! two-moment cloud droplets and rain and free ice categories whose particle
! properties evolve with time. The P3 source lives in
! P3_microphys/module_mp_p3.F90; this module converts CLUBB's fields to P3's
! inputs, calls P3 for all columns, and returns P3's changes to CLUBB as time
! tendencies.
!
! Usage:
!   - p3_microphys_init(): read the P3 lookup table and set the warm-rain
!     scheme. Call once during microphysics initialization.
!   - p3_microphys_driver(): call P3 for one timestep and return microphysics
!     tendencies. Call each timestep in place of another microphysics scheme,
!     or pass it to lh_microphys_driver, which calls it for each SILHS sample.
!
! Understanding the coupling:
!   - Configuration: one ice category with two-moment ice. Rain mass and
!     number (rrm, Nrm) are always CLUBB hydrometeors. With l_ice_microphys,
!     P3's ice mass, number, rime mass and rime volume are the hydrometeors
!     rim, Nim, qirim and birim; rime mass is part of rim. Without it, ice is
!     passed to P3 as zero and p3_microphys_driver stops if P3 forms ice. P3
!     only forms ice below 269 K (freezing) or 258 K (deposition nucleation).
!   - Cloud fraction: P3 runs without its subgrid cloud and precipitation
!     fraction option (SCPF), as in CM1. It treats each grid box as fully
!     cloudy, so cloud processes such as autoconversion use the grid means
!     as in-cloud values, and condensation and rain evaporation use the
!     grid-mean supersaturation. In partly cloudy layers this underestimates
!     autoconversion: for given grid means, the Seifert and Beheng (2001)
!     rate grows about as 1 / cloud fraction. CLUBB's cloud fraction is not
!     used. init_microphys prints a warning unless SILHS is interactive.
!     SCPF was tested with CLUBB and not adopted: its diagnosed cloud
!     fraction barely changed autoconversion, and rain then fell through
!     subsaturated air below cloud with little evaporation.
!   - SILHS: with lh_microphys_type = "interactive", SILHS draws sample
!     subcolumns from CLUBB's subgrid PDF and P3 runs on each. A sample is
!     fully cloudy or clear at each level, as a grid box in CM1, so P3 gets
!     in-cloud values and separate cloudy and clear vapor without a cloud
!     fraction. The tendencies are averaged with the sample weights. P3
!     sediments within each sample subcolumn. With ice enabled, sample rime
!     mass and volume follow sampled total ice at the grid-mean rime fraction
!     and density. Ice number is sampled separately. Cloud droplet number
!     must be prescribed (l_predict_Nc = .false.); the predicted-droplet-number
!     budget with SILHS has not been validated.
!   - Cloud water: CLUBB diagnoses rcm from its subgrid PDF every timestep.
!     P3 receives the grid-mean rcm, and its changes to cloud water and vapor
!     return as rcm_mc and rvm_mc. clubb_driver adds both to rtm_forcing at
!     the start of the next timestep, before advancing rtm. Condensation and
!     evaporation of cloud water therefore cancel in CLUBB's variables; only
!     conversion to rain, rain evaporation and droplet sedimentation change
!     rtm and thlm.
!   - Droplet number: P3 predicts droplet number and receives the grid-mean
!     Ncm, Nc_in_cloud times cloud_frac. Ncm_mc is used only when CLUBB predicts Ncm
!     (l_predict_Nc). Otherwise P3's change to droplet number is discarded,
!     and each timestep starts again from the prescribed Nc_in_cloud.
!   - Sedimentation: P3 sediments cloud droplets, rain and ice itself, so
!     hydromet_vel_zt is returned as zero and CLUBB adds no sedimentation.
!   - Precision: P3 uses default (single precision) reals. Inputs are
!     converted to that kind and results back to core_rknd.
!
! References:
!   Morrison, H. and J. A. Milbrandt, 2015: Parameterization of cloud
!   microphysics based on the prediction of bulk ice particle properties.
!   Part I: Scheme description and idealized tests. J. Atmos. Sci., 72,
!   287-311.
!-------------------------------------------------------------------------------
module p3_microphys_module

  use clubb_precision, only: &
      core_rknd ! Variable(s)

  implicit none

  public :: p3_microphys_init, &
            p3_microphys_driver, &
            p3_ice_tol

  private ! Default Scope

  ! Smallest ice mass mixing ratio P3 keeps (qsmall in module_mp_p3.F90).
  ! Newly nucleated ice is far below CLUBB's ri_tol, so P3's ice tolerances use
  ! this value instead [kg/kg]
  real( kind = core_rknd ), parameter :: &
    p3_ice_tol = 1.0e-14_core_rknd

  ! Warm-rain scheme P3 was initialized with; P3 initializes only once per process
  integer, save :: &
    p3_iparam_initialized = -1

!$omp threadprivate( p3_iparam_initialized )

  contains

  !=============================================================================
  subroutine p3_microphys_init( p3_iparam, p3_lookup_table_dir )

    ! Description:
    ! Set P3's warm-rain scheme and read its ice lookup table. P3 initializes
    ! only once per process, so a later call must use the same warm-rain scheme.
    !
    ! References:
    ! None
    !-----------------------------------------------------------------------

    use module_mp_p3, only: &
        p3_init,    & ! Procedure(s)
        iparam_host   ! Variable(s)

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use error_code, only: &
        clubb_at_least_debug_level_api ! Procedure(s)

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      p3_iparam ! Warm-rain scheme: 1 = SB2001, 2 = Beheng 1994, 3 = KK2000 [-]

    character(len=*), intent(in) :: &
      p3_lookup_table_dir ! Directory holding the P3 ice lookup table

    !--------------------------- Local Variables ---------------------------
    integer :: &
      p3_status ! Status returned by p3_init; 0 is success [-]

    !--------------------------- Begin Code ---------------------------

    if ( p3_iparam_initialized > 0 ) then

      if ( p3_iparam /= p3_iparam_initialized ) then
        write(fstderr,*) "P3 was already initialized with p3_iparam = ", &
                         p3_iparam_initialized, " and cannot change to ", p3_iparam
        error stop "Fatal error in p3_microphys_init."
      end if

      return

    end if

    iparam_host = p3_iparam

    call p3_init( lookup_file_dir = trim( p3_lookup_table_dir ), & ! In
                  nCat = 1, trplMomI = .false., model = 'WRF',   & ! In
                  abort_on_err = .true.,                         & ! In
                  dowr = clubb_at_least_debug_level_api( 1 ),    & ! In
                  stat = p3_status )                               ! Out

    if ( p3_status /= 0 ) then
      write(fstderr,*) "p3_init failed with status ", p3_status, &
                       " reading the lookup table from ", trim( p3_lookup_table_dir )
      error stop "Fatal error in p3_microphys_init."
    end if

    p3_iparam_initialized = p3_iparam

  end subroutine p3_microphys_init

  !=============================================================================
  subroutine p3_microphys_driver( gr, ngrdcol, dt, nzt,                 & ! In
                                  hydromet_dim, hm_metadata,            & ! In
                                  l_latin_hypercube,                    & ! In
                                  thlm, wm_zt, p_in_Pa,                 & ! In
                                  exner, rho, cloud_frac, w_std_dev,    & ! In
                                  dzq_in, rcm, Ncm, chi, rvm, hydromet, & ! In
                                  saturation_formula,                   & ! In
                                  sample_weight, stats,                 & ! In, InOut
                                  hydromet_mc, hydromet_vel_zt, Ncm_mc, & ! Out
                                  rcm_mc, rvm_mc, thlm_mc,              & ! Out
                                  rrm_auto_diag, rrm_accr_diag,         & ! Out
                                  rrm_evap_diag, Nrm_auto_diag,         & ! Out
                                  Nrm_evap_diag )                         ! Out

    ! Description:
    ! Call P3 for all columns for one timestep and return its changes as
    ! tendencies of CLUBB's fields. P3 updates the fields by process rates and
    ! sedimentation over dt; each tendency is (P3 result - input) / dt.
    !
    ! The argument list follows microphys_interface.inc, so SILHS can call
    ! this routine once per sample subcolumn (l_latin_hypercube = .true.).
    ! A sample is either fully cloudy or clear at each level, which is what
    ! P3 assumes without its cloud fraction option. SILHS averages the
    ! tendencies over the samples; statistics are then accumulated as
    ! weighted sub-timestep averages.
    !
    ! The thlm tendency is computed with CLUBB's constants from the change in
    ! vapor plus cloud water and from the water that freezes,
    !   thlm_mc = - Lv / ( Cp * exner ) * ( rcm_mc + rvm_mc )
    !             + Lf / ( Cp * exner ) * ri_frz.
    ! Liquid phase changes alter theta_l only when water enters or leaves
    ! vapor plus cloud water (rain formation, rain evaporation and droplet
    ! sedimentation). This matches the KK form, thlm_mc = Lv / ( Cp * exner )
    ! * rrm_mc, and avoids a spurious theta_l change from cloud condensation
    ! that would follow from P3's slightly different latent heat. The freezing
    ! rate ri_frz, the net rate at which vapor and liquid turn into ice, cannot
    ! be taken from the change in ice mass, which includes ice falling between
    ! levels. It follows from P3's heat budget, Cp_P3 * dT = -Lv_P3 * dqv +
    ! Lf_P3 * dri_frz, using P3's temperature and vapor changes.
    !
    ! References:
    ! Morrison and Milbrandt (2015), see the module description.
    !-----------------------------------------------------------------------

    use module_mp_p3, only: &
        p3_main,    & ! Procedure(s)
        n_p3diag_3d   ! Constant(s)

    use grid_class, only: &
        grid ! Type

    use corr_varnce_module, only: &
        hm_metadata_type ! Type

    use constants_clubb, only: &
        Lv,           & ! Constant(s)
        Lf,           &
        Cp,           &
        zero,         &
        sec_per_day,  &
        mm_per_m,     &
        micron_per_m, &
        fstderr

    use T_in_K_module, only: &
        thlm2T_in_K_api ! Procedure(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use stats_netcdf, only: &
        stats_type,  & ! Type
        stats_update   ! Procedure(s)

    use parameters_microphys, only: &
        lh_microphys_type,           & ! Variable(s)
        lh_microphys_non_interactive   ! Constant(s)

    implicit none

    !--------------------------- Local Constants ---------------------------
    integer, parameter :: &
      nCat      = 1, & ! Number of P3 ice categories [-]
      n_diag_2d = 1, & ! Number of P3 user-defined 2D diagnostics [-]
      it_p3     = 2    ! P3 timestep counter; any value > 1 skips P3's initial
                       ! saturation adjustment, which CLUBB's diagnosed rcm makes
                       ! unnecessary [-]

    ! P3's heat capacity and latent heats (module_mp_p3.F90: cp in p3_init,
    ! xxlv and xlf in p3_main), used to recover the freezing rate
    real( kind = core_rknd ), parameter :: &
      Cp_P3 = 1005.0_core_rknd,                                         & ! [J/(kg K)]
      Lv_P3 = 3.1484e6_core_rknd - 2370.0_core_rknd * 273.15_core_rknd, & ! [J/kg]
      Lf_P3 = 0.3337e6_core_rknd                                          ! [J/kg]

    !--------------------------- Input Variables ---------------------------
    type(grid), intent(in) :: &
      gr ! CLUBB grid

    integer, intent(in) :: &
      ngrdcol,     & ! Number of grid columns [-]
      nzt,         & ! Number of thermodynamic levels [-]
      hydromet_dim   ! Number of hydrometeor species [-]

    real( kind = core_rknd ), intent(in) :: &
      dt ! Model timestep [s]

    type(hm_metadata_type), intent(in) :: &
      hm_metadata ! Hydrometeor indices and properties

    logical, intent(in) :: &
      l_latin_hypercube ! Called by SILHS for one sample subcolumn

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      thlm,       & ! Liquid water potential temperature        [K]
      wm_zt,      & ! Vertical velocity                         [m/s]
      p_in_Pa,    & ! Pressure                                  [Pa]
      exner,      & ! Exner function                            [-]
      rho,        & ! Air density; unused                       [kg/m^3]
      cloud_frac, & ! Cloud fraction; unused (P3 runs without SCPF) [-]
      w_std_dev,  & ! Standard deviation of w; unused           [m/s]
      dzq_in,     & ! Distance between thermodynamic levels; unused, since
                    ! P3 needs the layer thickness gr%dzt         [m]
      rcm,        & ! Cloud water mixing ratio (grid mean or sample) [kg/kg]
      Ncm,        & ! Cloud droplet number (grid mean or sample)     [num/kg]
      chi,        & ! Extended liquid water mixing ratio; unused     [kg/kg]
      rvm           ! Vapor water mixing ratio                  [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeors (rain, and ice with l_ice_microphys) [units vary]

    integer, intent(in) :: &
      saturation_formula ! Saturation formula; unused (P3 uses its own) [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      sample_weight ! SILHS sample weight; unity for an ordinary call [-]

    !--------------------------- Input/Output Variables ---------------------------
    type(stats_type), intent(inout) :: &
      stats ! Statistics output

    !--------------------------- Output Variables ---------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(out) :: &
      hydromet_mc,   & ! Hydrometeor tendencies, including sedimentation [(units vary)/s]
      hydromet_vel_zt  ! Hydrometeor sedimentation velocity; zero        [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      Ncm_mc,  & ! Cloud droplet number tendency    [num/kg/s]
      rcm_mc,  & ! Cloud water tendency             [kg/kg/s]
      rvm_mc,  & ! Vapor water tendency             [kg/kg/s]
      thlm_mc    ! Liquid potential temp. tendency  [K/s]

    ! P3 does not return rain evaporation or rain number source rates, so
    ! those three are zero
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      rrm_auto_diag, & ! Autoconversion rate             [kg/kg/s]
      rrm_accr_diag, & ! Accretion rate                  [kg/kg/s]
      rrm_evap_diag, & ! Rain evaporation rate; zero     [kg/kg/s]
      Nrm_auto_diag, & ! Rain number from autoconversion; zero [num/kg/s]
      Nrm_evap_diag    ! Rain number evaporation; zero   [num/kg/s]

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      T_in_K, & ! Temperature [K]
      ri_frz    ! Net rate of freezing and deposition into ice; ice only [kg/kg/s]

    ! P3 inputs and outputs in P3's default real kind
    real, dimension(ngrdcol,nzt) :: &
      qc,     & ! Cloud water mixing ratio                     [kg/kg]
      nc,     & ! Cloud droplet number                         [num/kg]
      qr,     & ! Rain water mixing ratio                      [kg/kg]
      nr,     & ! Rain drop number                             [num/kg]
      th_old, & ! Potential temperature at the start of dt     [K]
      th,     & ! Potential temperature                        [K]
      th_in,  & ! Potential temperature passed to P3           [K]
      qv_old, & ! Vapor mixing ratio at the start of dt        [kg/kg]
      qv,     & ! Vapor mixing ratio                           [kg/kg]
      qv_in,  & ! Vapor mixing ratio passed to P3              [kg/kg]
      ssat,   & ! Supersaturation, qv - qvs; diagnosed by P3   [kg/kg]
      uzpl,   & ! Vertical velocity                            [m/s]
      pres,   & ! Pressure                                     [Pa]
      dzq,    & ! Thermodynamic layer thickness                [m]
      diag_ze,   & ! Equivalent reflectivity                   [dBZ]
      diag_effc, & ! Cloud droplet effective radius            [m]
      SCF_out      ! Cloud fraction from P3's SCPF option; unused [-]

    real, dimension(ngrdcol,nzt,nCat) :: &
      qitot,     & ! Total ice mass mixing ratio  [kg/kg]
      qirim,     & ! Rime ice mass mixing ratio   [kg/kg]
      nitot,     & ! Total ice number             [num/kg]
      birim,     & ! Rime ice volume              [m^3/kg]
      diag_effi, & ! Ice effective radius         [m]
      diag_vmi,  & ! Ice mass-weighted fall speed [m/s]
      diag_di,   & ! Mean ice diameter            [m]
      diag_rhoi    ! Bulk ice density             [kg/m^3]

    real, dimension(ngrdcol) :: &
      prt_liq, & ! Surface liquid precipitation rate [m/s]
      prt_sol    ! Surface solid precipitation rate  [m/s]

    real, dimension(ngrdcol,n_diag_2d) :: &
      diag_2d ! P3 user-defined 2D diagnostics; unused [units vary]

    ! Slots 3 and 4 hold the autoconversion and accretion rates [kg/kg/s]
    real, dimension(ngrdcol,nzt,n_p3diag_3d) :: &
      diag_3d ! P3 3D diagnostics [units vary]

    integer :: &
      iirr,    & ! Hydrometeor index of rain water mixing ratio [-]
      iiNr,    & ! Hydrometeor index of rain drop number        [-]
      iiri,    & ! Hydrometeor index of ice mixing ratio        [-]
      iiNi,    & ! Hydrometeor index of ice number              [-]
      iiqirim, & ! Hydrometeor index of rime ice mixing ratio   [-]
      iibirim    ! Hydrometeor index of rime ice volume         [-]

    logical :: &
      l_ice ! CLUBB carries P3's ice

    !--------------------------- Begin Code ---------------------------

    iirr    = hm_metadata%iirr
    iiNr    = hm_metadata%iiNr
    iiri    = hm_metadata%iiri
    iiNi    = hm_metadata%iiNi
    iiqirim = hm_metadata%iiqirim
    iibirim = hm_metadata%iibirim

    l_ice = iiri > 0

    ! Convert CLUBB's fields to P3's inputs. P3 updates th and qv; th_old and
    ! qv_old give the start-of-step values.
    T_in_K = thlm2T_in_K_api( nzt, ngrdcol, thlm, exner, rcm )

    th     = real( T_in_K / exner )
    th_old = th
    qv     = real( rvm )
    qv_old = qv
    qc     = real( rcm )
    nc     = real( Ncm )
    qr     = real( hydromet(:,:,iirr) )
    nr     = real( hydromet(:,:,iiNr) )
    uzpl   = real( wm_zt )
    pres   = real( p_in_Pa )
    dzq    = real( gr%dzt )
    th_in  = th
    qv_in  = qv

    if ( l_ice ) then
      qitot(:,:,1) = real( hydromet(:,:,iiri) )
      nitot(:,:,1) = real( hydromet(:,:,iiNi) )
      qirim(:,:,1) = real( hydromet(:,:,iiqirim) )
      birim(:,:,1) = real( hydromet(:,:,iibirim) )
    else
      ! Without l_ice_microphys, P3 starts each timestep without ice
      qitot = 0.0
      qirim = 0.0
      nitot = 0.0
      birim = 0.0
    end if

    ssat = 0.0

    call p3_main( qc = qc, nc = nc, qr = qr, nr = nr,                        & ! InOut
                  th_old = th_old, th = th, qv_old = qv_old, qv = qv,        & ! InOut
                  dt = real( dt ),                                           & ! In
                  qitot = qitot, qirim = qirim, nitot = nitot, birim = birim, & ! InOut
                  ssat = ssat,                                               & ! InOut
                  uzpl = uzpl, pres = pres, dzq = dzq, it = it_p3,           & ! In
                  prt_liq = prt_liq, prt_sol = prt_sol,                      & ! Out
                  its = 1, ite = ngrdcol, kts = 1, kte = nzt, nCat = nCat,   & ! In
                  diag_ze = diag_ze, diag_effc = diag_effc,                  & ! Out
                  diag_effi = diag_effi, diag_vmi = diag_vmi,                & ! Out
                  diag_di = diag_di, diag_rhoi = diag_rhoi,                  & ! Out
                  n_diag_2d = n_diag_2d,                                     & ! In
                  diag_2d = diag_2d,                                         & ! Out
                  n_diag_3d = n_p3diag_3d,                                   & ! In
                  diag_3d = diag_3d,                                         & ! Out
                  log_predictNc = .true., typeDiags_ON = .false.,            & ! In
                  model = 'WRF', clbfact_dep = 1.0, clbfact_sub = 1.0,       & ! In
                  debug_on = .false., scpf_on = .false.,                     & ! In
                  scpf_pfrac = 0.0, scpf_resfact = 0.0,                      & ! In
                  SCF_out = SCF_out )                                          ! Out

    if ( .not. l_ice .and. ( any( qitot > 0.0 ) .or. any( nitot > 0.0 ) ) ) then
      write(fstderr,*) "P3 formed ice; set l_ice_microphys to carry P3's ice in CLUBB."
      write(fstderr,*) "Maximum ice mass mixing ratio [kg/kg]: ", maxval( qitot )
      write(fstderr,*) "Maximum ice number [num/kg]: ", maxval( nitot )
      error stop "Fatal error in p3_microphys_driver."
    end if

    ! Return P3's changes as tendencies. Sedimentation is included in
    ! hydromet_mc, rcm_mc and Ncm_mc, so no sedimentation velocity is returned.
    ! P3 overwrites th_old and qv_old, so th_in and qv_in hold the inputs.
    hydromet_mc = zero
    hydromet_vel_zt = zero

    hydromet_mc(:,:,iirr) = ( real( qr, kind = core_rknd ) - hydromet(:,:,iirr) ) / dt
    hydromet_mc(:,:,iiNr) = ( real( nr, kind = core_rknd ) - hydromet(:,:,iiNr) ) / dt

    if ( l_ice ) then
      hydromet_mc(:,:,iiri) &
      = ( real( qitot(:,:,1), kind = core_rknd ) - hydromet(:,:,iiri) ) / dt
      hydromet_mc(:,:,iiNi) &
      = ( real( nitot(:,:,1), kind = core_rknd ) - hydromet(:,:,iiNi) ) / dt
      hydromet_mc(:,:,iiqirim) &
      = ( real( qirim(:,:,1), kind = core_rknd ) - hydromet(:,:,iiqirim) ) / dt
      hydromet_mc(:,:,iibirim) &
      = ( real( birim(:,:,1), kind = core_rknd ) - hydromet(:,:,iibirim) ) / dt
    end if

    Ncm_mc = ( real( nc, kind = core_rknd ) - Ncm ) / dt
    rcm_mc = ( real( qc, kind = core_rknd ) - rcm ) / dt
    rvm_mc = ( real( qv, kind = core_rknd ) - rvm ) / dt

    thlm_mc = - Lv / ( Cp * exner ) * ( rcm_mc + rvm_mc )

    if ( l_ice ) then
      ri_frz = ( Cp_P3 * real( th - th_in, kind = core_rknd ) * exner &
                 + Lv_P3 * real( qv - qv_in, kind = core_rknd ) ) / ( Lf_P3 * dt )
      thlm_mc = thlm_mc + Lf / ( Cp * exner ) * ri_frz
    end if

    rrm_auto_diag = real( diag_3d(:,:,3), kind = core_rknd )
    rrm_accr_diag = real( diag_3d(:,:,4), kind = core_rknd )
    rrm_evap_diag = zero
    Nrm_auto_diag = zero
    Nrm_evap_diag = zero

    call update_p3_stat( "rrm_auto", rrm_auto_diag )
    call update_p3_stat( "rrm_accr", rrm_accr_diag )

    call update_p3_stat( "eff_rad_cloud", real( diag_effc, kind = core_rknd ) * micron_per_m )

    ! P3 gives the solid rate as liquid-equivalent depth
    call update_p3_stat_sfc( "precip_rate_sfc", &
                             real( prt_liq + prt_sol, kind = core_rknd ) &
                             * mm_per_m * sec_per_day )

    if ( l_ice ) then
      call update_p3_stat( "eff_rad_ice", &
                           real( diag_effi(:,:,1), kind = core_rknd ) * micron_per_m )
    end if

  contains

    !===========================================================================
    subroutine update_p3_stat( name, value )

      ! Description:
      ! Output a profile statistic. For a SILHS sample, add the weighted value
      ! to a sub-timestep average over the samples. For SILHS without feedback
      ! (non-interactive), only the rates are output, under their lh_ names,
      ! since the ordinary call outputs the rest.
      !-------------------------------------------------------------------------

      implicit none

      character(len=*), intent(in) :: &
        name ! Statistics variable name

      real( kind = core_rknd ), dimension(:,:), intent(in) :: &
        value ! Profile to output [units vary]

      character(len=64) :: &
        output_name ! Name the value is output under

      if ( .not. stats%l_sample ) return

      output_name = name

      if ( l_latin_hypercube .and. &
           lh_microphys_type == lh_microphys_non_interactive ) then
        select case ( trim( name ) )
        case ( "rrm_auto", "rrm_accr" )
          output_name = "lh_" // trim( name )
        case default
          return
        end select
      end if

      if ( l_latin_hypercube ) then
        call stats_update( trim( output_name ), value * sample_weight, stats, &
                           sub_timestep_average = .true. )
      else
        call stats_update( trim( output_name ), value, stats )
      end if

    end subroutine update_p3_stat

    !===========================================================================
    subroutine update_p3_stat_sfc( name, value )

      ! Description:
      ! Output a surface statistic, averaged over SILHS samples as in
      ! update_p3_stat with the weight of the lowest level. Not output by
      ! SILHS without feedback.
      !-------------------------------------------------------------------------

      implicit none

      character(len=*), intent(in) :: &
        name ! Statistics variable name

      real( kind = core_rknd ), dimension(:), intent(in) :: &
        value ! Surface value to output [units vary]

      if ( .not. stats%l_sample .or. &
           ( l_latin_hypercube .and. &
             lh_microphys_type == lh_microphys_non_interactive ) ) return

      if ( l_latin_hypercube ) then
        call stats_update( name, value * sample_weight(:,1), stats, &
                           sub_timestep_average = .true. )
      else
        call stats_update( name, value, stats )
      end if

    end subroutine update_p3_stat_sfc

  end subroutine p3_microphys_driver

end module p3_microphys_module
