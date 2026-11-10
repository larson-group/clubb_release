!----------------------------------------------------------------------
! $Id$

module inputfields
!  Description:
!    This file contains information and subroutines needed for restart
!    simulations and simulations where the variables are reset at every
!    timestep, called a "clubb_inputfields" simulation.
!===============================================================================

  use clubb_precision, only: &
    core_rknd ! variable(s)

  use input_netcdf, only: &
    clubb_day, &
    clubb_month, &
    clubb_year

  implicit none

  ! Run information
  character(len=120), dimension(:), allocatable, public :: & 
    stat_files ! a list of all files used for input

  character(len=10), public :: input_type

  logical, public :: &
    l_input_um = .false., l_input_vm = .false., l_input_rtm = .false., l_input_thlm = .false., & 
    l_input_wp2 = .false., l_input_wprtp = .false., l_input_wpthlp = .false.,  & 
    l_input_wp3 = .false., l_input_rtp2 = .false., l_input_rtp3 = .false., &
    l_input_thlp2 = .false., l_input_thlp3 = .false., l_input_rtpthlp = .false., &
    l_input_upwp = .false., l_input_vpwp = .false., &
    l_input_ug = .false., l_input_vg = .false., l_input_rcm = .false.,  & 
    l_input_wm_zt = .false., l_input_exner = .false., l_input_em = .false., & 
    l_input_p = .false., l_input_rho = .false., l_input_rho_zm = .false., & 
    l_input_rho_ds_zm = .false., l_input_rho_ds_zt = .false., &
    l_input_thv_ds_zm = .false., l_input_thv_ds_zt = .false., &
    l_input_Lscale = .false., l_input_Lscale_up = .false., l_input_Lscale_down = .false., & 
    l_input_Kh_zt = .false., l_input_Kh_zm = .false., &
    l_input_tau_zm = .false., l_input_tau_zt = .false., & 
    l_input_wpthvp = .false., l_input_wp2thvp = .false., l_input_wp2up = .false., &
    l_input_rtpthvp = .false., l_input_thlpthvp = .false., &
    l_input_wp2rtp = .false., l_input_wp2thlp = .false., &
    l_input_uprcp = .false., l_input_vprcp = .false., &
    l_input_rc_coef_zm = .false., l_input_wp4 = .false., &
    l_input_wpup2 = .false., l_input_wpvp2 = .false., &
    l_input_wp2up2 = .false., l_input_wp2vp2 = .false., l_input_iss_frac = .false., &
    l_input_radht = .false., &
    l_input_w_1 = .false., l_input_w_2 = .false., &
    l_input_varnce_w_1 = .false., l_input_varnce_w_2 = .false., &
    l_input_rt_1 = .false., l_input_rt_2 = .false., &
    l_input_varnce_rt_1 = .false., l_input_varnce_rt_2 = .false., &
    l_input_thl_1 = .false., l_input_thl_2 = .false., &
    l_input_varnce_thl_1 = .false., l_input_varnce_thl_2 = .false., &
    l_input_mixt_frac = .false., &
    l_input_chi_1 = .false., l_input_chi_2 = .false., &
    l_input_stdev_chi_1 = .false., l_input_stdev_chi_2 = .false., &
    l_input_rc_1 = .false., l_input_rc_2 = .false., &
    l_input_w_1_zm = .false., l_input_w_2_zm = .false., &
    l_input_varnce_w_1_zm = .false., l_input_varnce_w_2_zm = .false., &
    l_input_mixt_frac_zm = .false., &
    l_input_thvm = .false., l_input_rrm = .false., &
    l_input_Nrm = .false.,  l_input_Ncm = .false.,  & 
    l_input_rsm = .false., l_input_rim = .false., &
    l_input_Nsm = .false., l_input_Ngm = .false., &
    l_input_rgm = .false., l_input_Nccnm = .false., l_input_Nim = .false., &
    l_input_rrp2 = .false., l_input_Nrp2 = .false., &
    l_input_wprrp = .false., l_input_wpNrp = .false., &
    l_input_thlm_forcing = .false., l_input_rtm_forcing = .false., & 
    l_input_up2 = .false., l_input_vp2 = .false., l_input_sigma_sqd_w = .false., & 
    l_input_cloud_frac = .false., l_input_sigma_sqd_w_zt = .false., &
    l_input_veg_T_in_K = .false., l_input_deep_soil_T_in_K = .false., &
    l_input_sfc_soil_T_in_K = .false., l_input_wprtp_forcing = .false., &
    l_input_wpthlp_forcing = .false., l_input_rtp2_forcing = .false., &
    l_input_thlp2_forcing = .false., l_input_rtpthlp_forcing = .false., &
    l_input_thlprcp = .false., l_input_rcm_mc = .false., l_input_rvm_mc = .false., &
    l_input_thlm_mc = .false., l_input_wprtp_mc = .false., l_input_wpthlp_mc = .false., &
    l_input_rtp2_mc = .false., l_input_thlp2_mc = .false., l_input_rtpthlp_mc = .false.

  integer, parameter, private :: &
    coamps_input_type = 1, &
    clubb_input_type =  2, &
    sam_input_type = 3, &
    rams_input_type = 4

  ! These represent the index of each file in the stat_files dimension for each input type
  integer, parameter, private :: &
    clubb_zt = 1, &
    clubb_zm = 2, &
    clubb_sfc = 3, &
    coamps_sm = 1, &
    coamps_sw = 2, &
    sam_file = 1, &
    rams_file =1

  integer, parameter, private :: &
    num_sam_inputfields = 70, & ! The number of input fields for SAM
    num_coamps = 61, & ! The number of input fields for coamps
    num_rams_inputfields =  4 ! The number of input fields for the RAMS LES case

  integer, private :: &
    stats_input_type = -999

  ! Procedures
  public  :: stat_fields_reader, &
             compute_timestep, &
             inputfields_init, &
             set_filenames, &
             cleanup_input_fields, &
             get_clubb_variable_interpolated

  private ! Default Scope

  type input_field
 
    logical :: l_input_var ! If .true., this variable should be read

    character(len = 50) :: &
      input_name, & ! The name of the variable in SAM 
      clubb_name  ! The name of the variable in CLUBB

    character(len = 3) :: &
      clubb_grid_type  ! The type of the grid for clubb, either "zm" or "zt"

    integer :: &
      input_file_index ! The index of the file this variable is located in

    real( kind = core_rknd ), dimension(:,:), pointer :: &
      clubb_var  ! The clubb variable to store the result in

    real( kind = core_rknd ) :: &
      adjustment ! The SAM variable will be multiplied by this amount
                 ! to convert the variable to the correct units.

  end type input_field

  contains

!===============================================================================


!-----------------------------------------------------------------------
  subroutine set_filenames( file_prefix )

! Description: Set the names of the netCDF files to be used.
!   Used by clubb_inputfields and clubb_restart.
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only: fstderr ! Constants

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variables
    character(len=*), intent(in) :: &
      file_prefix

    ! Local Variables
    logical :: l_stats_netcdf

    ! ---- Begin Code ----

    if ( allocated(stat_files) ) then
      deallocate(stat_files)
    end if

    stats_input_type = -1

    select case ( trim( input_type ) )

    case ( "coamps_les", "COAMPS_LES", "coamps" )
      allocate(stat_files(1:2))
      stat_files(coamps_sm) = trim( file_prefix )//"_coamps_sm.nc"
      stat_files(coamps_sw) = trim( file_prefix )//"_coamps_sw.nc"

      stats_input_type = coamps_input_type

    case ( "clubb_scm", "CLUBB_SCM", "clubb" )
      allocate(stat_files(1:3))
      stat_files(clubb_zt) = trim( file_prefix )//"_zt.nc"
      stat_files(clubb_zm) = trim( file_prefix )//"_zm.nc"
      stat_files(clubb_sfc) = trim( file_prefix )//"_sfc.nc"

      stats_input_type = clubb_input_type

      inquire(file=stat_files(clubb_zt), exist=l_stats_netcdf)
      
      if ( .not. l_stats_netcdf ) then
        write(fstderr,*) "inputfields: Expected netCDF stats file not found: ", &
                         trim(file_prefix)//"_stats.nc"
    
        stat_files(clubb_zt) = trim( file_prefix )// "_stats.nc"
        stat_files(clubb_zm) = trim( file_prefix )// "_stats.nc"
        stat_files(clubb_sfc) = trim( file_prefix )// "_stats.nc"
      end if
      
    case ( "SAM", "sam" )
      allocate(stat_files(1:1))
      ! SAM uses one netCDF file for all of its output, so the entire filename
      ! should be provided in file_prefix
      stat_files(sam_file) = trim( file_prefix )
      stats_input_type = sam_input_type

    case ( "rams_les", "RAMS_LES", "rams" )
      allocate(stat_files(1:1))
      ! RAMS uses one file for all of it's output.
      stat_files(rams_file) = trim( file_prefix )//"_rams.nc"

      stats_input_type = rams_input_type

    case default
      write(fstderr,*) "Don't know how to handle input_type = "// & 
        input_type
      error stop

    end select

    return
  end subroutine set_filenames
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine stat_fields_reader( gr, ngrdcol, timestep, hydromet_dim, & ! In
                                 hm_metadata, & ! In
                                 microphys_scheme, l_predict_Nc, & ! In
                                 um, upwp, vm, vpwp, & ! InOut
                                 up2, vp2, rtm, & ! InOut
                                 wprtp, thlm, wpthlp, & ! InOut
                                 rtp2, rtp3, & ! InOut
                                 thlp2, thlp3, rtpthlp, & ! InOut
                                 wp2, wp3, & ! InOut
                                 p_in_Pa, exner, rcm, cloud_frac, & ! InOut
                                 wpthvp, wp2thvp, wp2up, rtpthvp, thlpthvp, & ! InOut
                                 wp2rtp, wp2thlp, uprcp, vprcp, & ! InOut
                                 rc_coef_zm, wp4, wpup2, & ! InOut
                                 wpvp2, wp2up2, & ! InOut
                                 wp2vp2, ice_supersat_frac, & ! InOut
                                 wm_zt, rho, rho_zm, rho_ds_zm, & ! InOut
                                 rho_ds_zt, thv_ds_zm, thv_ds_zt, & ! InOut
                                 thlm_forcing, rtm_forcing, wprtp_forcing, & ! InOut
                                 wpthlp_forcing, rtp2_forcing, & ! InOut
                                 thlp2_forcing, rtpthlp_forcing, & ! InOut
                                 hydromet, hydrometp2, wphydrometp, & ! InOut
                                 Ncm, Nccnm, thvm, em, & ! InOut
                                 tau_zm, tau_zt, & ! InOut
                                 Kh_zt, Kh_zm, ug, vg, & ! InOut
                                 thlprcp, & ! InOut
                                 sigma_sqd_w, sigma_sqd_w_zt, radht, & ! InOut
                                 deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K, & ! InOut
                                 pdf_params, pdf_params_zm ) ! InOut

! Description:
!   Read shared statistical profiles and interpolate them onto every column.
!   Preserve fields outside each column's LES domain and apply its density conversions.

! References:
!   None
!-----------------------------------------------------------------------

    use grid_class, only: & 
        zt2zm_api ! Procedure(s)

    use constants_clubb, only:  &
        rt_tol,    & ! Variable(s)
        thl_tol,   &
        w_tol_sqd, &
        em_min,     &
        fstderr, &
        pascal_per_mb, &
        g_per_kg, &
        sec_per_day

    use pdf_parameter_module, only: &
        pdf_parameter    ! Type(s)

    use corr_varnce_module, only: &
        hm_metadata_type

    use stat_file_utils, only: & 
        LES_grid_to_CLUBB_grid, & ! Procedure(s)
        CLUBB_levels_within_LES_domain

    use extrapolation, only: &
        lin_ext_zt_bottom, &
        lin_ext_zm_bottom

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use grid_class, only: grid

    implicit none


    ! External
    intrinsic :: max, trim, any

    !--------------------------- Input Variables ---------------------------

    type (grid), target, intent(in) :: &
      gr

    integer, intent(in) :: &
      ngrdcol, & ! Number of columns [-]
      timestep

    integer, intent(in) :: &
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    character(len=*), intent(in) :: &
      microphys_scheme

    logical, intent(in) :: &
      l_predict_Nc

    !--------------------------- Input/Output Variables ---------------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), target, intent(inout) :: &
      um,                & ! eastward grid-mean wind component (thermo. levs.)  [m/s]
      vm,                & ! northward grid-mean wind component (thermo. levs.) [m/s]
      rtm,               & ! total water mixing ratio, r_t (thermo. levels)     [kg/kg]
      thlm,              & ! liq. water pot. temp., th_l (thermo. levels)       [K]
      rtp3,              & ! r_t'^3 (thermodynamic levels)                      [(kg/kg)^3]
      thlp3,             & ! th_l'^3 (thermodynamic levels)                     [K^3]
      wp3,               & ! w'^3 (thermodynamic levels)                        [m^3/s^3]
      p_in_Pa,           & ! Air pressure (thermodynamic levels)                [Pa]
      exner,             & ! Exner function (thermodynamic levels)              [-]
      rcm,               & ! cloud water mixing ratio, r_c (thermo. levels)     [kg/kg]
      cloud_frac,        & ! cloud fraction (thermodynamic levels)              [-]
      wp2thvp,           & ! < w'^2 th_v' > (thermodynamic levels)              [m^2/s^2 K]
      wp2up,             & ! < w'^2 u' > (thermodynamic levels)                 [m^3/s^3]
      wp2rtp,            & ! w'^2 rt' (thermodynamic levels)                    [m^2/s^2 kg/kg]
      wp2thlp,           & ! w'^2 thl' (thermodynamic levels)                   [m^2/s^2 K]
      wpup2,             & ! w'u'^2 (thermodynamic levels)                      [m^3/s^3]
      wpvp2,             & ! w'v'^2 (thermodynamic levels)                      [m^3/s^3]
      ice_supersat_frac    ! ice cloud fraction (thermo. levels)                [-]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), target, intent(inout) :: &
      upwp,       & ! u'w' (momentum levels)                         [m^2/s^2]
      vpwp,       & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,        & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,        & ! v'^2 (momentum levels)                         [m^2/s^2]
      wprtp,      & ! w' r_t' (momentum levels)                      [kg/kg m/s]
      wpthlp,     & ! w'th_l' (momentum levels)                      [(m/s) K]
      rtp2,       & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      thlp2,      & ! th_l'^2 (momentum levels)                      [K^2]
      rtpthlp,    & ! r_t'th_l' (momentum levels)                    [(kg/kg) K]
      wp2,        & ! w'^2 (momentum levels)                         [m^2/s^2]
      wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp,   & ! < th_l' th_v' > (momentum levels)              [K^2]
      uprcp,      & ! < u' r_c' > (momentum levels)                  [(m/s)(kg/kg)]
      vprcp,      & ! < v' r_c' > (momentum levels)                  [(m/s)(kg/kg)]
      rc_coef_zm, & ! Coef of X'r_c' in Eq. (34) (m-levs.)           [K/(kg/kg)]
      wp4,        & ! w'^4 (momentum levels)                         [m^4/s^4]
      wp2up2,     & ! w'^2 u'^2 (momentum levels)                    [m^4/s^4]
      wp2vp2        ! w'^2 v'^2 (momentum levels)                    [m^4/s^4]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), target, intent(inout) :: &
      wm_zt,     & ! vertical mean wind component on thermo. levels  [m/s]
      rho,       & ! Air density on thermodynamic levels             [kg/m^3]
      rho_ds_zt, & ! Dry, static density on thermo. levels           [kg/m^3]
      thv_ds_zt    ! Dry, base-state theta_v on thermo levels        [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), target, intent(inout) :: &
      rho_zm,    & ! Air density on momentum levels                  [kg/m^3]
      rho_ds_zm, & ! Dry, static density on momentum levels          [kg/m^3]
      thv_ds_zm    ! Dry, base-state theta_v on momentum levels      [K]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), target, intent(inout) :: &
      thlm_forcing,    & ! liquid potential temp. forcing (thermo. levels) [K/s]
      rtm_forcing        ! total water forcing (thermo. levels)      [(kg/kg)/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), target, intent(inout) :: &
      wprtp_forcing,   & ! total water turbulent flux forcing (m-levs) [m*K/s^2]
      wpthlp_forcing,  & ! liq pot temp turb flux forcing (m-levs)[m(kg/kg)/s^2]
      rtp2_forcing,    & ! total water variance forcing (m-levs)   [(kg/kg)^2/s]
      thlp2_forcing,   & ! liq pot temp variance forcing (m-levs)  [K^2/s]
      rtpthlp_forcing    ! <r_t'th_l'> covariance forcing (m-levs) [K*(kg/kg)/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,hydromet_dim), intent(inout) :: &
      hydromet       ! Array of hydrometeors                [hm units]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm,hydromet_dim), intent(inout) :: &
      hydrometp2,  & ! Variance of a hydrometeor (m-levs.)  [<hm units>^2]
      wphydrometp    ! Covariance of w and a hydrometeor    [(m/s) <hm units>]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      Ncm       ! Mean cloud droplet concentration, <N_c> (t-levs.)    [num/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), target, intent(inout) :: &
      Nccnm,          & ! Cloud condensation nuclei concentration (COAMPS/MG)  [num/kg]
      thvm,           & ! Virtual potential temperature                        [K]
      tau_zt,         & ! Eddy dissipation time scale on thermodynamic levels  [s]
      Kh_zt,          & ! Eddy diffusivity coefficient on thermodynamic levels [m^2/s]
      ug,             & ! u geostrophic wind                                   [m/s]
      vg,             & ! v geostrophic wind                                   [m/s]
      sigma_sqd_w_zt, & ! PDF width parameter interpolated to t-levs.          [-]
      radht             ! SW + LW heating rate                                 [K/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), target, intent(inout) :: &
      thlprcp,     & ! thl'rc'                                              [K kg/kg]
      sigma_sqd_w, & ! PDF width parameter (momentum levels)                [-]
      em,          & ! Turbulent Kinetic Energy (TKE)                       [m^2/s^2]
      tau_zm,      & ! Eddy dissipation time scale on momentum levels       [s]
      Kh_zm          ! Eddy diffusivity coefficient on momentum levels      [m^2/s]
    
    real( kind = core_rknd ), dimension(ngrdcol), intent(inout) :: &
      deep_soil_T_in_K, & ! Deep soil temperature [K]
      sfc_soil_T_in_K, & ! Surface soil temperature [K]
      veg_T_in_K ! Vegetation temperature [K]

    type(pdf_parameter), intent(inout) :: &
      pdf_params,    & ! PDF parameters (thermodynamic levels)    [units vary]
      pdf_params_zm    ! PDF parameters on momentum levels        [units vary]

    !--------------------------- Local Variables ---------------------------
    logical :: l_read_error, l_fatal_error

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) :: tmp1 ! Read/interpolation scratch [units vary]

    integer, dimension(ngrdcol,size(stat_files)) ::  &
      k_lowest_zt, &  ! The lowest CLUBB thermodynamic level that's within the LES domain.
      k_highest_zt, & ! The highest CLUBB thermodynamic level that's within the LES domain.
      k_lowest_zm, &  ! The lowest CLUBB momentum level that's within the LES domain.
      k_highest_zm    ! The highest CLUBB momentum level that's within the LES domain.

    real( kind = core_rknd ), dimension(ngrdcol,1) :: surface_heights ! Surface input height [m]

    integer :: i, k ! Column and vertical/field indices
 
    real( kind = core_rknd), dimension(ngrdcol,gr%nzt), target :: &
      temp_Nrm, temp_Ncm, temp_rgm, temp_rim, temp_Ngm, temp_Nsm, &
      temp_rrm, temp_rsm, temp_Nim, up2_zt, vp2_zt ! temp. variables

    real( kind = core_rknd), dimension(ngrdcol,gr%nzm), target :: &
      temp_tke, temp_wpup_sgs, temp_wpvp_sgs, &
      temp_rrp2, temp_Nrp2, temp_wprrp, temp_wpNrp ! temp. variables

    type (input_field), dimension(:), allocatable :: &
      SAM_variables, & ! A list of SAM variables to read in.
      coamps_variables, & ! A list of coamps variables to read in.
      RAMS_variables ! A list of RAMS LES variables to read in.


    integer :: & 
      iirr, iiNr, iirs, iiri, iirg, iiNi, iiNg, iiNs

    !--------------------------- Begin Code ---------------------------

    surface_heights = 0.0_core_rknd

    iirr = hm_metadata%iirr
    iiNr = hm_metadata%iiNr
    iirs = hm_metadata%iirs
    iiri = hm_metadata%iiri
    iirg = hm_metadata%iirg
    iiNi = hm_metadata%iiNi
    iiNg = hm_metadata%iiNg
    iiNs = hm_metadata%iiNs

    select case ( stats_input_type )


    !-------------------------------------
    ! CLUBB stats data
    !-------------------------------------
    case ( clubb_input_type )

      !  Thermo grid - zt file

      ! Initialize l_fatal_error for case clubb_input_type
      l_fatal_error = .false.

      call get_clubb_variable_interpolated &
           ( l_input_um, stat_files(clubb_zt), "um", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             um, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vm, stat_files(clubb_zt), "vm", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             vm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtm, stat_files(clubb_zt), "rtm", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             rtm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_thlm, stat_files(clubb_zt), "thlm", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             thlm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wp3, stat_files(clubb_zt), "wp3", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             wp3, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_tau_zt, stat_files(clubb_zt), "tau_zt", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             tau_zt, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rrm, stat_files(clubb_zt), "rrm", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             tmp1(:,1:gr%nzt), & ! InOut
             l_read_error ) ! Out
      if ( l_input_rrm ) then
        hydromet(:,1:gr%nzt,iirr) = tmp1(:,1:gr%nzt)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rsm, stat_files(clubb_zt), "rsm", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             tmp1(:,1:gr%nzt), & ! InOut
             l_read_error ) ! Out
      if ( l_input_rsm ) then
        hydromet(:,1:gr%nzt,iirs) = tmp1(:,1:gr%nzt)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rim, stat_files(clubb_zt), "rim", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             tmp1(:,1:gr%nzt), & ! InOut
             l_read_error ) ! Out
      if ( l_input_rim ) then
        hydromet(:,1:gr%nzt,iiri) = tmp1(:,1:gr%nzt)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rgm, stat_files(clubb_zt), "rgm", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             tmp1(:,1:gr%nzt), & ! InOut
             l_read_error ) ! Out
      if ( l_input_rgm ) then
        hydromet(:,1:gr%nzt,iirg) = tmp1(:,1:gr%nzt)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

!--------------------------------------------------------
! Added variables for clubb_restart
      call get_clubb_variable_interpolated &
           ( l_input_p, stat_files(clubb_zt), "p_in_Pa", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             p_in_Pa, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_exner, stat_files(clubb_zt), "exner", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             exner, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_ug, stat_files(clubb_zt), "ug", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             ug, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vg, stat_files(clubb_zt), "vg", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             vg, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rcm, stat_files(clubb_zt), "rcm", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             rcm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wm_zt, stat_files(clubb_zt), "wm_zt", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             wm_zt, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_rho, stat_files(clubb_zt), "rho", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             rho, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rho_ds_zt, stat_files(clubb_zt), "rho_ds_zt", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             rho_ds_zt, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thv_ds_zt, stat_files(clubb_zt), "thv_ds_zt", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             thv_ds_zt, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Kh_zt, stat_files(clubb_zt), "Kh_zt", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             Kh_zt, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thvm, stat_files(clubb_zt), "thvm", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             thvm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thlm_forcing, stat_files(clubb_zt), "thlm_forcing", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             thlm_forcing, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtm_forcing, stat_files(clubb_zt), "rtm_forcing", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             rtm_forcing, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Ncm, stat_files(clubb_zt), "Ncm", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             tmp1(:,1:gr%nzt), & ! InOut
             l_read_error ) ! Out
      if ( l_input_Ncm ) then
         Ncm(:,1:gr%nzt) = tmp1(:,1:gr%nzt)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Nccnm, stat_files(clubb_zt), "Nccnm", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             Nccnm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Nim, stat_files(clubb_zt), "Nim", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             tmp1(:,1:gr%nzt), & ! InOut
             l_read_error ) ! Out
      if ( l_input_Nim ) then
        hydromet(:,1:gr%nzt, iiNi) = tmp1(:,1:gr%nzt)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_cloud_frac, stat_files(clubb_zt), "cloud_frac", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             cloud_frac, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Nrm, stat_files(clubb_zt), "Nrm", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             tmp1(:,1:gr%nzt), & ! InOut
             l_read_error ) ! Out
      if ( l_input_Nrm ) then
        hydromet(:,1:gr%nzt, iiNr) = tmp1(:,1:gr%nzt)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_sigma_sqd_w_zt, stat_files(clubb_zt), "sigma_sqd_w_zt", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             sigma_sqd_w_zt, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wp2thvp, stat_files(clubb_zt), "wp2thvp", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             wp2thvp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wp2up, stat_files(clubb_zt), "wp2up", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             wp2up, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wp2rtp, stat_files(clubb_zt), "wp2rtp", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             wp2rtp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wp2thlp, stat_files(clubb_zt), "wp2thlp", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             wp2thlp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wpup2, stat_files(clubb_zt), "wpup2", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             wpup2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wpvp2, stat_files(clubb_zt), "wpvp2", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             wpvp2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_iss_frac, stat_files(clubb_zt), "ice_supersat_frac", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             ice_supersat_frac, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_radht, stat_files(clubb_zt), "radht", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             radht, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      ! PDF Parameters (needed for CLUBB restarts)
      call get_clubb_variable_interpolated &
           ( l_input_w_1, stat_files(clubb_zt), "w_1", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%w_1, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_w_2, stat_files(clubb_zt), "w_2", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%w_2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_varnce_w_1, stat_files(clubb_zt), "varnce_w_1", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%varnce_w_1, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_varnce_w_2, stat_files(clubb_zt), "varnce_w_2", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%varnce_w_2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rt_1, stat_files(clubb_zt), "rt_1", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%rt_1, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rt_2, stat_files(clubb_zt), "rt_2", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%rt_2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_varnce_rt_1, stat_files(clubb_zt), "varnce_rt_1", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%varnce_rt_1, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_varnce_rt_2, stat_files(clubb_zt), "varnce_rt_2", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%varnce_rt_2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thl_1, stat_files(clubb_zt), "thl_1", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%thl_1, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thl_2, stat_files(clubb_zt), "thl_2", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%thl_2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_varnce_thl_1, stat_files(clubb_zt), "varnce_thl_1", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%varnce_thl_1, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_varnce_thl_2, stat_files(clubb_zt), "varnce_thl_2", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%varnce_thl_2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_mixt_frac, stat_files(clubb_zt), "mixt_frac", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%mixt_frac, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_chi_1, stat_files(clubb_zt), "chi_1", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%chi_1, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_chi_2, stat_files(clubb_zt), "chi_2", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%chi_2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_stdev_chi_1, stat_files(clubb_zt), "stdev_chi_1", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%stdev_chi_1, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_stdev_chi_2, stat_files(clubb_zt), "stdev_chi_2", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%stdev_chi_2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rc_1, stat_files(clubb_zt), "rc_1", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%rc_1, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rc_2, stat_files(clubb_zt), "rc_2", & ! In
             gr%nzt, ngrdcol, timestep, gr%zt, & ! In
             pdf_params%rc_2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error




!--------------------------------------------------------

      ! Read in the zm file

      call get_clubb_variable_interpolated &
           ( l_input_wp2, stat_files(clubb_zm), "wp2", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             wp2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wprtp, stat_files(clubb_zm), "wprtp", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             wprtp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wpthlp, stat_files(clubb_zm), "wpthlp", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             wpthlp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wpthvp, stat_files(clubb_zm), "wpthvp", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             wpthvp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtpthvp, stat_files(clubb_zm), "rtpthvp", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             rtpthvp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thlpthvp, stat_files(clubb_zm), "thlpthvp", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             thlpthvp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtp2, stat_files(clubb_zm), "rtp2", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             rtp2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thlp2, stat_files(clubb_zm), "thlp2", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             thlp2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtpthlp, stat_files(clubb_zm), "rtpthlp", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             rtpthlp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_upwp, stat_files(clubb_zm), "upwp", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             upwp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vpwp, stat_files(clubb_zm), "vpwp", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             vpwp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

!-----------------------------------------------------------
      call get_clubb_variable_interpolated &
           ( l_input_em, stat_files(clubb_zm), "em", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             em, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rho_zm, stat_files(clubb_zm), "rho_zm", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             rho_zm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rho_ds_zm, stat_files(clubb_zm), "rho_ds_zm", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             rho_ds_zm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thv_ds_zm, stat_files(clubb_zm), "thv_ds_zm", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             thv_ds_zm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Kh_zm, stat_files(clubb_zm), "Kh_zm", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             Kh_zm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_tau_zm, stat_files(clubb_zm), "tau_zm", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             tau_zm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_up2, stat_files(clubb_zm), "up2", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             up2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vp2, stat_files(clubb_zm), "vp2", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             vp2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wp4, stat_files(clubb_zm), "wp4", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             wp4, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_uprcp, stat_files(clubb_zm), "uprcp", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             uprcp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vprcp, stat_files(clubb_zm), "vprcp", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             vprcp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wp2up2, stat_files(clubb_zm), "wp2up2", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             wp2up2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wp2vp2, stat_files(clubb_zm), "wp2vp2", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             wp2vp2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_sigma_sqd_w, stat_files(clubb_zm), "sigma_sqd_w", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             sigma_sqd_w, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wprtp_forcing, stat_files(clubb_zm), "wprtp_forcing", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             wprtp_forcing, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_wpthlp_forcing, stat_files(clubb_zm), "wpthlp_forcing", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             wpthlp_forcing, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_rtp2_forcing, stat_files(clubb_zm), "rtp2_forcing", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             rtp2_forcing, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_thlp2_forcing, stat_files(clubb_zm), "thlp2_forcing", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             thlp2_forcing, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtpthlp_forcing, stat_files(clubb_zm), "rtpthlp_forcing", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             rtpthlp_forcing, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_thlprcp, stat_files(clubb_zm), "thlprcp", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             thlprcp, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rc_coef_zm, stat_files(clubb_zm), "rc_coef_zm", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             rc_coef_zm, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      ! PDF Parameters (needed for CLUBB restarts)
      call get_clubb_variable_interpolated &
           ( l_input_w_1_zm, stat_files(clubb_zm), "w_1_zm", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             pdf_params_zm%w_1, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_w_2_zm, stat_files(clubb_zm), "w_2_zm", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             pdf_params_zm%w_2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_varnce_w_1_zm, stat_files(clubb_zm), "varnce_w_1_zm", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             pdf_params_zm%varnce_w_1, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_varnce_w_2_zm, stat_files(clubb_zm), "varnce_w_2_zm", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             pdf_params_zm%varnce_w_2, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_mixt_frac_zm, stat_files(clubb_zm), "mixt_frac_zm", & ! In
             gr%nzm, ngrdcol, timestep, gr%zm, & ! In
             pdf_params_zm%mixt_frac, & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error


!-----------------------------------------------------------

      tmp1(:,1) = veg_T_in_K

      call get_clubb_variable_interpolated &
           ( l_input_veg_T_in_K, stat_files(clubb_sfc), "veg_T_in_K", & ! In
             1, ngrdcol, timestep, surface_heights, & ! In
             tmp1(:,1:1), & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error
      veg_T_in_K = tmp1(:,1)

      tmp1(:,1) = deep_soil_T_in_K

      call get_clubb_variable_interpolated &
           ( l_input_deep_soil_T_in_K, stat_files(clubb_sfc), "deep_soil_T_in_K", & ! In
             1, ngrdcol, timestep, surface_heights, & ! In
             tmp1(:,1:1), & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error
      deep_soil_T_in_K = tmp1(:,1)

      tmp1(:,1) = sfc_soil_T_in_K

      call get_clubb_variable_interpolated &
           ( l_input_sfc_soil_T_in_K, stat_files(clubb_sfc), "sfc_soil_T_in_K", & ! In
             1, ngrdcol, timestep, surface_heights, & ! In
             tmp1(:,1:1), & ! InOut
             l_read_error ) ! Out

      l_fatal_error = l_fatal_error .or. l_read_error
      sfc_soil_T_in_K = tmp1(:,1)



      if ( l_fatal_error ) error stop "oops, get_grads_var failed in stat_fields_reader"


    !---------------------------------------
    ! COAMPS LES stats data
    !--------------------------------------

    case ( coamps_input_type )    
      l_fatal_error = .false.

      allocate( coamps_variables( num_coamps ) )

      k = 1

      coamps_variables(k)%l_input_var = l_input_um
      coamps_variables(k)%input_name = "um"
      coamps_variables(k)%clubb_var => um
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_vm
      coamps_variables(k)%input_name = "vm"
      coamps_variables(k)%clubb_var => vm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rtm
      coamps_variables(k)%input_name = "qtm"
      coamps_variables(k)%clubb_var => rtm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_thlm
      coamps_variables(k)%input_name = "thlm"
      coamps_variables(k)%clubb_var => thlm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      ! We obtain wp2 from stats_sw

      coamps_variables(k)%l_input_var = l_input_wp3
      coamps_variables(k)%input_name = "wp3"
      coamps_variables(k)%clubb_var => wp3
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_wprtp
      coamps_variables(k)%input_name = "wpqtp"
      coamps_variables(k)%clubb_var => wprtp
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_wpthlp
      coamps_variables(k)%input_name = "wpthlp"
      coamps_variables(k)%clubb_var => wpthlp
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rtp2
      coamps_variables(k)%input_name = "qtp2"
      coamps_variables(k)%clubb_var => rtp2
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_thlp2
      coamps_variables(k)%input_name = "thlp2"
      coamps_variables(k)%clubb_var => thlp2
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rtpthlp
      coamps_variables(k)%input_name = "qtpthlp"
      coamps_variables(k)%clubb_var => rtpthlp
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_ug
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "ug"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_vg
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "vg"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rcm
      coamps_variables(k)%input_name = "qcm"
      coamps_variables(k)%clubb_var => rcm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_wm_zt
      coamps_variables(k)%input_name = "wlsm"
      coamps_variables(k)%clubb_var => wm_zt
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_exner
      coamps_variables(k)%input_name = "ex0"
      coamps_variables(k)%clubb_var => exner
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      ! Note that em (clubb) = em (coamps) + tke (coamps), so just store tke
      ! in a temporary variable for now and add them together later.
      coamps_variables(k)%l_input_var = l_input_em
      coamps_variables(k)%input_name = "em"
      coamps_variables(k)%clubb_var => em
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      temp_tke = 0.0_core_rknd! Initialize temp to 0.0

      coamps_variables(k)%l_input_var = l_input_em
      coamps_variables(k)%input_name = "tke"
      coamps_variables(k)%clubb_var => temp_tke
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_p
      coamps_variables(k)%input_name = "pm"
      coamps_variables(k)%clubb_var => p_in_Pa
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rho
      coamps_variables(k)%input_name = "dn0"
      coamps_variables(k)%clubb_var => rho
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1


      ! rho_zm is in stats_sw


      coamps_variables(k)%l_input_var = l_input_Kh_zt
      coamps_variables(k)%input_name = "kh"
      coamps_variables(k)%clubb_var => Kh_zt
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_Kh_zm
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "Kh_zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_tau_zt
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "tau_zt"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_tau_zm
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "tau_zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_wpthvp
      coamps_variables(k)%input_name = "wpthvp"
      coamps_variables(k)%clubb_var => wpthvp
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_thl_1
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "thl_1"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_thl_2
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "thl_2"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_mixt_frac
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "mixt_frac"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_chi_1
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "chi_1"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_chi_2
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "chi_2"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_stdev_chi_1
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "stdev_chi_1"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_stdev_chi_2
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "stdev_chi_2"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rc_1
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "rc_1"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rc_2
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "rc_2"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_thvm
      coamps_variables(k)%input_name = "thvm"
      coamps_variables(k)%clubb_var => thvm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      temp_rrm = 0.0_core_rknd! initialize to 0.0

      if ( l_input_rrm ) then
        if ( iirr < 1 ) then
            write(fstderr,*) "Rain water mixing ratio cannot be input with"// &
              " microphys_scheme = "//microphys_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_rrm
          coamps_variables(k)%input_name = "qrm"
          coamps_variables(k)%clubb_var => temp_rrm
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%clubb_grid_type = "zt"
          coamps_variables(k)%input_file_index = coamps_sm

          k = k + 1
        end if
      end if

      temp_Nrm = 0.0_core_rknd! Initialize to 0.0

      if ( l_input_Nrm ) then
        if ( iiNr < 1 ) then
            write(fstderr,*) "Rain droplet number conc. cannot be input with"// &
              " microphys_scheme = "//microphys_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_Nrm
          coamps_variables(k)%input_name = "nrm"
          coamps_variables(k)%clubb_var => temp_Nrm
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%clubb_grid_type = "zt"
          coamps_variables(k)%input_file_index = coamps_sm

          k = k + 1
        end if
      end if

      temp_Ncm = 0.0_core_rknd! initialize to 0.0

      if ( l_input_Ncm ) then
        if ( .not. l_predict_Nc ) then
            write(fstderr,*) "Cloud droplet number conc. cannot be input with"// &
              " microphys_scheme = "//microphys_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_Ncm
          coamps_variables(k)%input_name = "ncm"
          coamps_variables(k)%clubb_var => temp_Ncm
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%clubb_grid_type = "zt"
          coamps_variables(k)%input_file_index = coamps_sm

          k = k + 1
        end if
      end if

      temp_rsm = 0.0_core_rknd! initialize to 0.0

      if ( l_input_rsm ) then
        if ( iirs < 1 ) then
            write(fstderr,*) "Snow mixing ratio cannot be input with"// &
              " microphys_scheme = "//microphys_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_rsm
          coamps_variables(k)%input_name = "qsm"
          coamps_variables(k)%clubb_var => temp_rsm
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%clubb_grid_type = "zt"
          coamps_variables(k)%input_file_index = coamps_sm

          k = k + 1
        end if
      end if

      temp_rim = 0.0_core_rknd! initialize to 0.0

      if ( l_input_rim ) then
        if ( iiri < 1 ) then
            write(fstderr,*) "Ice mixing ratio cannot be input with"// &
              " microphys_scheme = "//microphys_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_rim
          coamps_variables(k)%input_name = "qim"
          coamps_variables(k)%clubb_var => temp_rim
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%clubb_grid_type = "zt"
          coamps_variables(k)%input_file_index = coamps_sm

          k = k + 1
        end if
      end if

      temp_rgm = 0.0_core_rknd! initialize to 0.0

      if ( l_input_rgm ) then
        if ( iirg < 1 ) then
            write(fstderr,*) "Graupel mixing ratio cannot be input with"// &
              " microphys_scheme = "//microphys_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_rgm
          coamps_variables(k)%input_name = "qgm"
          coamps_variables(k)%clubb_var => temp_rgm
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%clubb_grid_type = "zt"
          coamps_variables(k)%input_file_index = coamps_sm

          k = k + 1
        end if
      end if

      coamps_variables(k)%l_input_var = l_input_Nccnm
      coamps_variables(k)%input_name = "nccnm"
      coamps_variables(k)%clubb_var => Nccnm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zt"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      temp_Nim = 0.0_core_rknd! initialize to 0.0

      if ( l_input_Nim ) then
        if ( iiNi < 1 ) then
            write(fstderr,*) "Ice number conc. cannot be input with"// &
              " microphys_scheme = "//microphys_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_Nim
          coamps_variables(k)%input_name = "nim"
          coamps_variables(k)%clubb_var => temp_Nim
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%clubb_grid_type = "zt"
          coamps_variables(k)%input_file_index = coamps_sm

          k = k + 1
        end if
      end if

      coamps_variables(k)%l_input_var = l_input_thlm_forcing
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "thlm_forcing"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rtm_forcing
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "rtm_forcing"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_up2
      coamps_variables(k)%input_name = "up2"
      coamps_variables(k)%clubb_var => up2
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_vp2
      coamps_variables(k)%input_name = "vp2"
      coamps_variables(k)%clubb_var => vp2
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sm

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_sigma_sqd_w
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "sigma_sqd_w"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_sigma_sqd_w_zt
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "sigma_sqd_w_zt"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_cloud_frac
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "cloud_frac"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_veg_T_in_K
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "veg_T_in_K"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_deep_soil_T_in_K
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "deep_soil_T_in_K"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_sfc_soil_t_in_K
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "sfc_soil_T_in_K"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_radht
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "radht"

      k = k + 1

      ! Note:  wpup_sgs and wpvp_sgs must be added to make the u'w' and v'w' terms
      !        as they are in CLUBB.
      coamps_variables(k)%l_input_var = l_input_upwp
      coamps_variables(k)%input_name = "wpup"
      coamps_variables(k)%clubb_var => upwp
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sw

      k = k + 1

      temp_wpup_sgs = 0.0_core_rknd! initialize to 0.0

      coamps_variables(k)%l_input_var = l_input_upwp
      coamps_variables(k)%input_name = "wpup_sgs"
      coamps_variables(k)%clubb_var => temp_wpup_sgs
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sw

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_vpwp
      coamps_variables(k)%input_name = "wpvp"
      coamps_variables(k)%clubb_var => vpwp
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sw

      k = k + 1

      temp_wpvp_sgs = 0.0_core_rknd! initialize to 0.0

      coamps_variables(k)%l_input_var = l_input_vpwp
      coamps_variables(k)%input_name = "wpvp_sgs"
      coamps_variables(k)%clubb_var => temp_wpvp_sgs
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sw

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_wp2
      coamps_variables(k)%input_name = "wp2"
      coamps_variables(k)%clubb_var => wp2
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sw

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rho_zm
      coamps_variables(k)%input_name = "dn0"
      coamps_variables(k)%clubb_var => rho_zm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%clubb_grid_type = "zm"
      coamps_variables(k)%input_file_index = coamps_sw

      ! Initialize l_read_error for case ( "les" )
      l_read_error = .false.

      if ( l_fatal_error ) error stop "oops, bad input for inputfields"
      
      call get_input_variables_interp &
           ( gr, ngrdcol, k, coamps_variables, timestep, & ! In
             k_lowest_zt, k_highest_zt, & ! Out
             k_lowest_zm, k_highest_zm, l_fatal_error ) ! Out

      if ( l_input_rrm ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,coamps_sm) .or. k > k_highest_zt(i,coamps_sm) ) cycle
            hydromet(i,k,iirr) =  temp_rrm(i,k)
          end do
        end do
      end if

      if ( l_input_Nrm ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,coamps_sm) .or. k > k_highest_zt(i,coamps_sm) ) cycle
            hydromet(i,k,iiNr) =  temp_Nrm(i,k)
          end do
        end do
      end if

      if ( l_input_Ncm ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,coamps_sm) .or. k > k_highest_zt(i,coamps_sm) ) cycle
            Ncm(i,k) =  temp_Ncm(i,k)
          end do
        end do
      end if

      if ( l_input_rsm ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,coamps_sm) .or. k > k_highest_zt(i,coamps_sm) ) cycle
            hydromet(i,k,iirs) =  temp_rsm(i,k)
          end do
        end do
      end if

      if ( l_input_rim ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,coamps_sm) .or. k > k_highest_zt(i,coamps_sm) ) cycle
            hydromet(i,k,iiri) =  temp_rim(i,k)
          end do
        end do
      end if

      if ( l_input_rgm ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,coamps_sm) .or. k > k_highest_zt(i,coamps_sm) ) cycle
            hydromet(i,k,iirg) =  temp_rgm(i,k)
          end do
        end do
      end if

      if ( l_input_Nim ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,coamps_sm) .or. k > k_highest_zt(i,coamps_sm) ) cycle
            hydromet(i,k,iiNi) =  temp_Nim(i,k)
          end do
        end do
      end if

      if ( l_input_wprtp ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wprtp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like surface_varnce.
        do i = 1, ngrdcol
          if ( k_lowest_zm(i,coamps_sm) == 2 ) then
            wprtp(i,1) = lin_ext_zm_bottom( wprtp(i,3), wprtp(i,2), &
            gr%zm(i,3), gr%zm(i,2), gr%zm(i,1) )
          end if
        end do
      endif

      if ( l_input_wpthlp ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wpthlp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like surface_varnce.
        do i = 1, ngrdcol
          if ( k_lowest_zm(i,coamps_sm) == 2 ) then
            wpthlp(i,1) = lin_ext_zm_bottom( wpthlp(i,3), wpthlp(i,2), &
            gr%zm(i,3), gr%zm(i,2), gr%zm(i,1) )
          end if
        end do
      endif

      if ( l_input_rtp2 ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, set the value of rtp2
        ! at momentum level 1 to the value at momentum level 2.
        ! Using a linear extension here resulted in negatives.
        do i = 1, ngrdcol
          if ( k_lowest_zm(i,coamps_sm) == 2 ) then
            rtp2(i,1) =  rtp2(i,2)
          end if
        end do

        where ( rtp2 < rt_tol**2 )
          rtp2 = rt_tol**2
        end where
      endif

      if ( l_input_thlp2 ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, set the value of thlp2
        ! at momentum level 1 to the value at momentum level 2.
        ! Using a linear extension here resulted in negatives.
        do i = 1, ngrdcol
          if ( k_lowest_zm(i,coamps_sm) == 2 ) then
            thlp2(i,1) = thlp2(i,2)
          end if
        end do
        
        where ( thlp2 < thl_tol**2 )
          thlp2 = thl_tol**2
        end where
      endif

      if ( l_input_rtpthlp) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! rtpthlp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like surface_varnce.
        do i = 1, ngrdcol
          if ( k_lowest_zm(i,coamps_sm) == 2 ) then
            rtpthlp(i,1) = lin_ext_zm_bottom( rtpthlp(i,3), rtpthlp(i,2), &
            gr%zm(i,3), gr%zm(i,2), gr%zm(i,1) )
          end if
        end do
      endif

      if( l_input_em ) then
        do k = 1, gr%nzm
          do i = 1, ngrdcol
            if ( k < k_lowest_zm(i,coamps_sm) .or. k > k_highest_zm(i,coamps_sm) ) cycle
            em(i,k) = em(i,k) + temp_tke(i,k)
          end do
        end do

        where ( em < em_min ) em = em_min
      end if

      if ( l_input_up2 ) then
        ! Clip up2 to be no smaller than w_tol_sqd
        where ( up2 < w_tol_sqd ) 
          up2 = w_tol_sqd
        end where
      endif

      if ( l_input_vp2 ) then
         where ( vp2 < w_tol_sqd ) 
           vp2 = w_tol_sqd
         end where
      endif

       if ( l_input_upwp) then
        do k = 1, gr%nzm
          do i = 1, ngrdcol
            if ( k < k_lowest_zm(i,coamps_sw) .or. k > k_highest_zm(i,coamps_sw) ) cycle
            upwp(i,k) = temp_wpup_sgs(i,k)
          end do
        end do
      endif

      if ( l_input_vpwp) then
        do k = 1, gr%nzm
          do i = 1, ngrdcol
            if ( k < k_lowest_zm(i,coamps_sw) .or. k > k_highest_zm(i,coamps_sw) ) cycle
            vpwp(i,k) = temp_wpvp_sgs(i,k)
          end do
        end do
      endif

      if ( l_input_wp2 ) then
        where ( wp2 < w_tol_sqd )
          wp2 = w_tol_sqd
        end where
      end if

      deallocate(coamps_variables)
    


    !-------------------------------------------
    ! SAM stats data
    !-------------------------------------------
    case( sam_input_type )

      allocate( SAM_variables(num_sam_inputfields) )

      l_fatal_error = .false.

      k = 1

      SAM_variables(k)%l_input_var = l_input_um
      SAM_variables(k)%input_name = "U"
      SAM_variables(k)%clubb_var => um
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_vm
      SAM_variables(k)%input_name = "V"
      SAM_variables(k)%clubb_var => vm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rtm
      SAM_variables(k)%input_name = "RTM"
      SAM_variables(k)%clubb_var => rtm
      SAM_variables(k)%adjustment = 1.0_core_rknd !1.0e-3_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thlm
      SAM_variables(k)%input_name = "THLM"
      SAM_variables(k)%clubb_var => thlm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_wp2
      SAM_variables(k)%input_name = "WP2"
      SAM_variables(k)%clubb_var => wp2
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rho
      SAM_variables(k)%input_name = "RHO"
      SAM_variables(k)%clubb_var => rho
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      ! Note that this needs to be adjusted by 1/(RHO * Lv)
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_wprtp
      SAM_variables(k)%input_name = "WPRTP"
      SAM_variables(k)%clubb_var => wprtp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      ! Note that this needs to be adjusted by 1/(RHO * Cp)
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_wpthlp
      SAM_variables(k)%input_name = "WPTHLP"
      SAM_variables(k)%clubb_var => wpthlp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_wp3
      SAM_variables(k)%input_name = "W3"
      SAM_variables(k)%clubb_var => wp3
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rtp2
      SAM_variables(k)%input_name = "RTP2"
      SAM_variables(k)%clubb_var => rtp2
      SAM_variables(k)%adjustment = 1.0_core_rknd !1.0e-6_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rtp3
      SAM_variables(k)%input_name = "QTO3"
      SAM_variables(k)%clubb_var => rtp3
      SAM_variables(k)%adjustment = 1.0e-9_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thlp2
      SAM_variables(k)%input_name = "THLP2"
      SAM_variables(k)%clubb_var => thlp2
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thlp3
      SAM_variables(k)%input_name = "THEL3"
      SAM_variables(k)%clubb_var => thlp3
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_upwp
      SAM_variables(k)%input_name = "UW"
      SAM_variables(k)%clubb_var => upwp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_vpwp
      SAM_variables(k)%input_name = "VW"
      SAM_variables(k)%clubb_var => vpwp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_ug
      SAM_variables(k)%input_name = "UOBS"
      SAM_variables(k)%clubb_var => ug
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_vg
      SAM_variables(k)%input_name = "VOBS"
      SAM_variables(k)%clubb_var => vg
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rcm
      SAM_variables(k)%input_name = "QCL"
      SAM_variables(k)%clubb_var => rcm
      SAM_variables(k)%adjustment = 1.0e-3_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_wm_zt
      SAM_variables(k)%input_name = "WOBS"
      SAM_variables(k)%clubb_var => wm_zt
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_exner
      SAM_variables(k)%clubb_name = "exner"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_em
      SAM_variables(k)%input_name = "TKE"
      SAM_variables(k)%clubb_var => em
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_p
      SAM_variables(k)%input_name = "PRES"
      SAM_variables(k)%clubb_var => p_in_Pa
      SAM_variables(k)%adjustment = pascal_per_mb
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rho_zm
      SAM_variables(k)%clubb_name = "rho_zm"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rho_ds_zm
      SAM_variables(k)%clubb_name = "rho_ds_zm"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rho_ds_zt
      SAM_variables(k)%clubb_name = "rho_ds_zt"
      SAM_variables(k)%input_name = "none"

      k = k + 1
       
      SAM_variables(k)%l_input_var = l_input_thv_ds_zm
      SAM_variables(k)%clubb_name = "thv_ds_zm"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thv_ds_zt
      SAM_variables(k)%clubb_name = "thv_ds_zt"
      SAM_variables(k)%input_name = "none" 

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_Kh_zt
      SAM_variables(k)%input_name = "TKH"
      SAM_variables(k)%clubb_var => Kh_zt
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_Kh_zm
      SAM_variables(k)%clubb_name = "Kh_zm"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_tau_zm
      SAM_variables(k)%clubb_name = "tau_zm"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_tau_zt
      SAM_variables(k)%clubb_name = "tau_zt"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thl_1
      SAM_variables(k)%clubb_name = "thl_1"
      SAM_variables(k)%input_name = "none"

      k = k + 1
      SAM_variables(k)%l_input_var = l_input_thl_2
      SAM_variables(k)%clubb_name = "thl_2"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_mixt_frac
      SAM_variables(k)%clubb_name = "mixt_frac"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_chi_1
      SAM_variables(k)%clubb_name = "chi_1"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_chi_2
      SAM_variables(k)%clubb_name = "chi_2"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_stdev_chi_1
      SAM_variables(k)%clubb_name = "stdev_chi_1"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_stdev_chi_2
      SAM_variables(k)%clubb_name = "stdev_chi_2"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rc_1
      SAM_variables(k)%clubb_name = "rc_1"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rc_2
      SAM_variables(k)%clubb_name = "rc_2"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thvm
      SAM_variables(k)%input_name = "THETAV"
      SAM_variables(k)%clubb_var => thvm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_rrm = 0.0_core_rknd! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_rrm
      SAM_variables(k)%input_name = "QPL"
      SAM_variables(k)%clubb_var => temp_rrm
      SAM_variables(k)%adjustment = 1.0_core_rknd/g_per_kg
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_Nrm = 0.0_core_rknd! initialize to 0.0

      ! Note that this needs to be adjusted by 1e6/RHO
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_Nrm
      SAM_variables(k)%input_name = "NR" ! SAM Morrison microphysics
      !SAM_variables(k)%input_name = "CONP" ! SAM KK microphysics
      SAM_variables(k)%clubb_var => temp_Nrm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1
      temp_Ncm = 0.0_core_rknd! initialize to 0.0

      ! Note that this needs to be adjusted by 1e6/RHO
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_Ncm
      SAM_variables(k)%input_name = "NC"
      SAM_variables(k)%clubb_var => temp_Ncm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_rsm = 0.0_core_rknd! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_rsm
      SAM_variables(k)%input_name = "QS"
      SAM_variables(k)%clubb_var => temp_rsm
      SAM_variables(k)%adjustment = 1.0_core_rknd/g_per_kg
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_rim = 0.0_core_rknd! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_rim
      SAM_variables(k)%input_name = "QI"
      SAM_variables(k)%clubb_var => temp_rim
      SAM_variables(k)%adjustment = 1.0_core_rknd/g_per_kg
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_rgm = 0.0_core_rknd ! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_rgm
      SAM_variables(k)%input_name = "QG"
      SAM_variables(k)%clubb_var => temp_rgm
      SAM_variables(k)%adjustment = 1.0_core_rknd/g_per_kg
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_Nccnm
      SAM_variables(k)%clubb_name = "Nccnm"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      temp_Nim = 0.0_core_rknd ! initialize to 0.0

      ! Note that this needs to be adjusted by 1e6/RHO
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_Nim
      SAM_variables(k)%input_name = "NI"
      SAM_variables(k)%clubb_var => temp_Nim
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_Nsm = 0.0_core_rknd ! initialize to 0.0

      ! Note that this needs to be adjusted by 1e6/RHO
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_Nsm
      SAM_variables(k)%input_name = "NS"
      SAM_variables(k)%clubb_var => temp_Nsm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_Ngm = 0.0_core_rknd ! initialize to 0.0

      ! Note that this needs to be adjusted by 1e6/RHO
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_Ngm
      SAM_variables(k)%input_name = "NG"
      SAM_variables(k)%clubb_var => temp_Ngm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_rrp2 = 0.0_core_rknd ! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_rrp2
      SAM_variables(k)%input_name = "rrp2" ! SAM Morrison microphysics
      !SAM_variables(k)%input_name = "RRP2" ! SAM KK microphysics
      SAM_variables(k)%clubb_var => temp_rrp2
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_Nrp2 = 0.0_core_rknd ! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_Nrp2
      SAM_variables(k)%input_name = "Nrp2" ! SAM Morrison microphysics
      !SAM_variables(k)%input_name = "NRP2" ! SAM KK microphysics
      SAM_variables(k)%clubb_var => temp_Nrp2
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_wprrp = 0.0_core_rknd ! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_wprrp
      SAM_variables(k)%input_name = "WPRRP"
      SAM_variables(k)%clubb_var => temp_wprrp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_wpNrp = 0.0_core_rknd ! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_wpNrp
      SAM_variables(k)%input_name = "WPNRP"
      SAM_variables(k)%clubb_var => temp_wpNrp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thlm_forcing
      SAM_variables(k)%clubb_name = "thlm_forcing"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rtm_forcing
      SAM_variables(k)%clubb_name = "rtm_forcing"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      ! Note that U2 is on the thermodynamic grid in SAM while up2 is
      ! on the momentum grid in CLUBB. This will have to be interpolated.
      SAM_variables(k)%l_input_var = l_input_up2
      SAM_variables(k)%input_name = "U2"
      SAM_variables(k)%clubb_var => up2_zt
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      ! Note that V2 is on the thermodynamic grid in SAM while vp2 is
      ! on the momentum grid in CLUBB. This will have to be interpolated.
      SAM_variables(k)%l_input_var = l_input_vp2
      SAM_variables(k)%input_name = "V2"
      SAM_variables(k)%clubb_var => vp2_zt
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_sigma_sqd_w
      SAM_variables(k)%clubb_name = "sigma_sqd_w"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_cloud_frac
      SAM_variables(k)%input_name = "CLD"
      SAM_variables(k)%clubb_var => cloud_frac
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_sigma_sqd_w_zt
      SAM_variables(k)%clubb_name = "sigma_sqd_w_zt"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_veg_T_in_K
      SAM_variables(k)%clubb_name = "veg_T_in_K"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_deep_soil_T_in_K
      SAM_variables(k)%clubb_name = "deep_soil_T_in_K"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_sfc_soil_T_in_K
      SAM_variables(k)%clubb_name = "sfc_soil_T_in_K"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_radht
      SAM_variables(k)%clubb_name = "radht"
      SAM_variables(k)%input_name = "RADQR"
      SAM_variables(k)%clubb_var => radht
      SAM_variables(k)%adjustment = 1.0_core_rknd/sec_per_day
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rtpthlp
      SAM_variables(k)%clubb_name = "rtpthlp"
      SAM_variables(k)%input_name = "RTPTHLP"
      SAM_variables(k)%clubb_var => rtpthlp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file


      up2_zt = 0.0_core_rknd
      vp2_zt = 0.0_core_rknd

      call get_input_variables_interp &
           ( gr, ngrdcol, k, SAM_variables, timestep, & ! In
             k_lowest_zt, k_highest_zt, & ! Out
             k_lowest_zm, k_highest_zm, l_fatal_error ) ! Out

      if( l_fatal_error ) then
        STOP "Error reading SAM input fields"
      end if

      ! Convert number concentrations using each column's density.
      if ( l_input_Nrm ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            temp_Nrm(i,k) = temp_Nrm(i,k) * (1.0e6_core_rknd / rho(i,k))
          end do
        end do
      end if

      if ( l_input_Ncm ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            temp_Ncm(i,k) = temp_Ncm(i,k) * (1.0e6_core_rknd / rho(i,k))
          end do
        end do
      end if

      if ( l_input_Nim ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            temp_Nim(i,k) = temp_Nim(i,k) * (1.0e7_core_rknd / rho(i,k))
          end do
        end do
      end if

      if ( l_input_Nsm ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            temp_Nsm(i,k) = temp_Nsm(i,k) * (1.0e7_core_rknd / rho(i,k))
          end do
        end do
      end if

      if ( l_input_Ngm ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            temp_Ngm(i,k) = temp_Ngm(i,k) * (1.0e7_core_rknd / rho(i,k))
          end do
        end do
      end if

      ! Interpolate all columns, preserving levels outside each LES domain.
      if ( l_input_up2 ) then
        tmp1 = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, up2_zt ) ! In
        do k = 1, gr%nzm
          do i = 1, ngrdcol
            if ( k >= k_lowest_zm(i,sam_file) .and. k <= k_highest_zm(i,sam_file) ) &
              up2(i,k) = tmp1(i,k)
          end do
        end do
      end if

      if ( l_input_vp2 ) then
        tmp1 = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, vp2_zt ) ! In
        do k = 1, gr%nzm
          do i = 1, ngrdcol
            if ( k >= k_lowest_zm(i,sam_file) .and. k <= k_highest_zm(i,sam_file) ) &
              vp2(i,k) = tmp1(i,k)
          end do
        end do
      end if

      ! Add hydrometeor variables.
      if( l_input_Nrm .and. iiNr > 0) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            hydromet(i,k,iiNr) =  temp_Nrm(i,k)
          end do
        end do
      end if

      if ( l_input_Ncm .and. l_predict_Nc ) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            Ncm(i,k) =  temp_Ncm(i,k)
          end do
        end do
      end if

      if( l_input_Nim .and. iiNi > 0) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            hydromet(i,k,iiNi) =  temp_Nim(i,k)
          end do
        end do
      end if

      if( l_input_rrm .and. iirr > 0) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            hydromet(i,k,iirr) =  temp_rrm(i,k)
          end do
        end do
      end if

      if( l_input_Ngm .and. iiNg > 0) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            hydromet(i,k,iiNg) =  temp_Ngm(i,k)
          end do
        end do
      end if

      if( l_input_rgm .and. iirg > 0) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            hydromet(i,k,iirg) =  temp_rgm(i,k)
          end do
        end do
      end if

      if( l_input_rim .and. iiri > 0) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            hydromet(i,k,iiri) =  temp_rim(i,k)
          end do
        end do
      end if

      if( l_input_Nsm .and. iiNs > 0) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            hydromet(i,k,iiNs) =  temp_Nsm(i,k)
          end do
        end do
      end if

      if( l_input_rsm .and. iirs > 0) then
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
            hydromet(i,k,iirs) =  temp_rsm(i,k)
          end do
        end do
      end if

      if ( l_input_rrp2 .and. iirr > 0 ) then
         do k = 1, gr%nzt
           do i = 1, ngrdcol
             if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
             hydrometp2(i,k,iirr)  = temp_rrp2(i,k)
           end do
         end do
      endif

      if ( l_input_Nrp2 .and. iiNr > 0 ) then
         do k = 1, gr%nzt
           do i = 1, ngrdcol
             if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
             hydrometp2(i,k,iiNr)  = temp_Nrp2(i,k)
           end do
         end do
      endif

      if ( l_input_wprrp .and. iirr > 0 ) then
         do k = 1, gr%nzt
           do i = 1, ngrdcol
             if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
             wphydrometp(i,k,iirr)  = temp_wprrp(i,k)
           end do
         end do
      endif

      if ( l_input_wpNrp .and. iiNr > 0 ) then
         do k = 1, gr%nzt
           do i = 1, ngrdcol
             if ( k < k_lowest_zt(i,sam_file) .or. k > k_highest_zt(i,sam_file) ) cycle
             wphydrometp(i,k,iiNr)  = temp_wpNrp(i,k)
           end do
         end do
      endif

      if ( l_fatal_error ) error stop "Failed to read inputfields for SAM."

      deallocate( SAM_variables )



    !---------------------------------------------------
    ! RAMS LES stats data
    !---------------------------------------------------
    case ( rams_input_type )
      allocate( RAMS_variables(num_rams_inputfields) )

      l_fatal_error = .false.

      k = 1

      RAMS_variables(k)%l_input_var = l_input_wp3
      RAMS_variables(k)%input_name = "wp3"
      RAMS_variables(k)%clubb_var => wp3
      RAMS_variables(k)%adjustment = 1.0_core_rknd
      RAMS_variables(k)%clubb_grid_type = "zt"
      RAMS_variables(k)%input_file_index = rams_file

      k = k + 1

      RAMS_variables(k)%l_input_var = l_input_wp2
      RAMS_variables(k)%input_name = "wp2"
      RAMS_variables(k)%clubb_var => wp2
      RAMS_variables(k)%adjustment = 1.0_core_rknd
      RAMS_variables(k)%clubb_grid_type = "zm"
      RAMS_variables(k)%input_file_index = rams_file

      k = k + 1

      RAMS_variables(k)%l_input_var = l_input_up2
      RAMS_variables(k)%input_name = "up2"

      RAMS_variables(k)%clubb_var => up2
      RAMS_variables(k)%adjustment = 1.0_core_rknd
      RAMS_variables(k)%clubb_grid_type = "zm"
      RAMS_variables(k)%input_file_index = rams_file

      k = k + 1

      RAMS_variables(k)%l_input_var = l_input_vp2
      RAMS_variables(k)%input_name = "vp2"
      RAMS_variables(k)%clubb_var => vp2
      RAMS_variables(k)%adjustment = 1.0_core_rknd
      RAMS_variables(k)%clubb_grid_type = "zm"
      RAMS_variables(k)%input_file_index = rams_file

      call get_input_variables_interp &
           ( gr, ngrdcol, k, RAMS_variables, timestep, & ! In
             k_lowest_zt, k_highest_zt, & ! Out
             k_lowest_zm, k_highest_zm, l_fatal_error ) ! Out

      if( l_fatal_error ) then
        STOP "Error reading RAMS LES input fields"
      end if

      deallocate( RAMS_variables )

    end select

    ! Clipping on the variance of u, v and w 
    if ( l_input_wp2 ) then 
      where ( wp2 < w_tol_sqd ) 
        wp2 = w_tol_sqd 
      end where 
    end if 

    if ( l_input_up2 ) then 
      where ( up2 < w_tol_sqd )  
        up2 = w_tol_sqd 
      end where 
    end if 

    if ( l_input_vp2 ) then 
      where ( vp2 < w_tol_sqd ) 
        vp2 = w_tol_sqd 
      end where 
    end if 

    return
  end subroutine stat_fields_reader

!===============================================================================                     
  subroutine compute_timestep( filename, l_restart, &                                        
                               time, nearest_timestep )                                            
                                                                                                     
    ! Description:                                                                                   
    ! Given a time 'time', determines the closest output time in a NetCDF file.                      
                                                                                                     
    !-----------------------------------------------------------------------                         
                                                                                                     
    use stat_file_module, only: &                                                                    
        stat_file     ! Type                                                                         
    use stat_file_utils, only: &
        l_netcdf_file

    use constants_clubb, only:  &                                                                    
        sec_per_min, & ! Variable(s)                                                                 
        fstderr ! Constant(s)                                                                        
                                                                                                     
    use clubb_precision, only:  &                                                                    
        time_precision                                                                               
                                                                                                     
    use input_netcdf, only: open_netcdf_read, close_netcdf_read ! Procedure(s)                       
                                                                                                     
    implicit none                                                                                    
                                                                                                     
    ! Input Variable(s)                                                                              
    character(len=*), intent(in) ::filename ! Path to the file and its name                          
                                                                                                     
    logical, intent(in) :: l_restart      ! Whether this is a restart run                            
                                                                                                     
    real( kind = time_precision ), intent(in) ::  &                                                  
      time ! Time near which we want to find netCDF output,                                          
    ! e.g. time_restart      [s]                                                                     
                                                                                                     
                                                                                                     
    ! Output Variable(s)                                                                             
    integer, intent(out) ::  &                                                                       
      nearest_timestep ! Nearest integer output time to time [units vary]                            
                                                                                                     
    ! Local Variables                                                                                
    type (stat_file) :: fread_var                                                                    
                                                                                                     
    real( kind = core_rknd ) :: delta_time   ! In seconds                                            
                                                                                                     
    logical :: l_error                                                                               
                                                                                                     
    ! ---- Begin Code ----                                                                           
    
    ! Initialize error state
    l_error = .false.

    ! Explicitly reject non-netCDF files
    if ( .not. l_netcdf_file( filename ) ) then
      write(fstderr,*) "compute_timestep: Expected netCDF file, got: "//trim(filename)
      error stop "Fatal error: Non-netCDF input."
    end if
                                                                                                     
    if( stats_input_type == sam_input_type ) then                                                  
      call open_netcdf_read( 'U', trim( filename ), & ! In                                         
                             fread_var, & ! In/Out                                                 
                             l_error ) ! Out                                                       
    else                                                                                           
      call open_netcdf_read( 'thlm', trim( filename ), & ! In                                      
                             fread_var, & ! In/Out                                                 
                             l_error ) ! Out                                                       
    end if                                                                                         
                                                                                                     
    if ( l_error ) then                                                                              
      write(fstderr,*) "Error reading file " // trim( filename )                                     
      error stop "Fatal error"                                                                       
    end if                                                                                           
                                                                                                     
    !  Eric   Raut      June  2013                                                                   
    delta_time = real(time - fread_var%time, kind=core_rknd)                                         
    !   Joshua Fasching March 2008                                                                   
!     .        time - fread_var%time                                                                 
                                                                                                     
    ! Reporting                                                                                      
    if ( l_restart ) then                                                                            
      print *, "Initial time of netCDF reference file ", &                                           
               "[seconds since midnight]: ",  &                                                      
               fread_var%time                                                                        
      print *, "Model restart time [s]: ", time                                                      
      print *, "Elapsed time between ", &                                                            
               "initial time of ref file and restart time [s]: ",  &                                 
               delta_time                                                                            
      print *, "netCDF file output time interval [s]: ",  &                                          
               fread_var%dtwrite                                                                     
                                                                                                     
      if ( ( mod( delta_time , fread_var%dtwrite )  > 1e-8_core_rknd) .or.  &                        
           ( mod( delta_time, fread_var%dtwrite ) < -1e-8_core_rknd) ) then                          
        write(fstderr,*) "Error: Elapsed time is not a multiple ", &                                 
                "of the reference netCDF output time interval."                                      
        write(fstderr,*) "Elapsed time [s] = ", delta_time                                           
        write(fstderr,*) "netCDF output time interval = ", fread_var%dtwrite                         
        error stop                                                                                   
      end if                                                                                         
                                                                                                     
      if ( mod( delta_time , sec_per_min ) > 1e-8_core_rknd&                                         
            .or. mod( delta_time, sec_per_min ) < -1e-8_core_rknd) then                              
        write(fstderr,*) "Error: Elapsed time is not a multiple ", &                                 
                "of one minute."                                                                     
        write(fstderr,*) "Elapsed time [s] = ", delta_time                                           
        error stop                                                                                   
      end if                                                                                         
                                                                                                     
    end if ! l_restart                                                                               
                                                                                                     
    ! Determines the closest recorded timestep to the restart                                        
    ! time.                                                                                          
    nearest_timestep = nint( delta_time / fread_var%dtwrite )                                        
                                                                                                     
    if ( l_restart ) then                                                                            
      print *, "Elapsed time between ", &                                                            
               "initial time of ref file and restart time ", &                                       
               "rounded to nearest minute: ",  &                                                     
               nearest_timestep                                                                      
                                                                                                     
      ! Print the actual record being recalled.                                                      
      ! Joshua Fasching March 2008                                                                   
      print *, "Nearest netCDF output time iteration [ ]: ", &                                       
               nint( real( nearest_timestep, kind=core_rknd) /  &                                    
                     fread_var%dtwrite/sec_per_min ) - 1                                             
    end if ! l_restart                                                                               
                                                                                                     
    call close_netcdf_read( fread_var )                                                              
                                                                                                     
  end subroutine compute_timestep

!===============================================================================
  subroutine get_clubb_variable_interpolated &
           ( l_input_var, filename, varname, & ! In
             vardim, ngrdcol, timestep, clubb_heights, & ! In
             variable_interpolated, & ! InOut
             l_read_error ) ! Out
! Description:
!   Read one shared file profile and interpolate it onto all column grids.
! References:
!   None
!--------------------------------------------------------------------------------
    use stat_file_utils, only: stat_file_average ! Procedure(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: &
      npower = 1 ! Power to raise the profile to (not applied)

    ! Input Variables
    logical, intent(in) :: &
      l_input_var !  Whether to read the variable in

    character(len=*), intent(in) :: &
      filename, & ! Filename
      varname     ! Variable name

    integer, intent(in) :: &
      vardim, & ! Number of vertical levels [-]
      ngrdcol, & ! Number of columns [-]
      timestep  ! Timestep to read in

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,vardim) :: &
      clubb_heights ! Altitudes         [m]

    !--------------------------- Input/Output Variables ---------------------------
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,vardim) :: &
      variable_interpolated     ! The variable after interpolation

    !--------------------------- Output Variables ---------------------------
    logical, intent(out) :: &
      l_read_error ! Whether there was a read error

    !--------------------------- Local Variables ---------------------------
    logical, dimension(ngrdcol) :: l_spec_bound_cond ! Lower ghost policy per column [-]

    !--------------------------- Begin Code ---------------------------

    if ( l_input_var ) then
      l_spec_bound_cond = clubb_heights(:,1) < 0.0_core_rknd

      variable_interpolated = stat_file_average( filename, vardim, ngrdcol, & ! In
                                                 timestep, timestep, clubb_heights, & ! In
                                                 varname, & ! In
                                                 npower, l_spec_bound_cond, & ! In
                                                 l_read_error ) ! Out

    else
      l_read_error = .false.

    end if

    return
  end subroutine get_clubb_variable_interpolated

  !--------------------------------------------------------------------------
  subroutine get_input_variables_interp( gr, ngrdcol, nvars, variables , timestep, &
                                         k_lowest_zt, k_highest_zt,   &
                                         k_lowest_zm, k_highest_zm, l_error )

  ! Description:
  !   Obtains profiles for COAMPS or SAM variables from a NetCDF file and
  !   interpolates if needed.
  ! References:
  !   None
  !--------------------------------------------------------------------------------

    use stat_file_module, only: stat_file

    use input_netcdf, only: &
      open_netcdf_read, &
      get_netcdf_var, &
      close_netcdf_read

    use stat_file_utils, only: & 
       LES_grid_to_CLUBB_grid, & ! Procedure(s)
       CLUBB_levels_within_LES_domain

    use constants_clubb, only: &
      fstderr  ! Variable(s)

    use grid_class, only: grid

    implicit none

    type (grid), target, intent(in) :: gr

    integer, intent(in) :: &
      ngrdcol, & ! Number of columns [-]
      nvars, & ! The size of the variables array
      timestep ! The timestep to read the variables at.

    type (input_field), dimension(nvars), intent(in) :: &
      variables ! The list of external variables to be read in.
  
    ! These variables are output by this subroutine in case certain variables
    ! need special adjustments.
    integer, dimension(ngrdcol,size(stat_files)), intent(out) ::  &
      k_lowest_zt, &  ! The lowest CLUBB thermodynamic level that's within the LES domain.
      k_highest_zt, & ! The highest CLUBB thermodynamic level that's within the LES domain.
      k_lowest_zm, &  ! The lowest CLUBB momentum level that's within the LES domain.
      k_highest_zm    ! The highest CLUBB momentum level that's within the LES domain.

    logical, intent(out) :: &
      l_error ! Error flag

    !--------------------------- Local Variables ---------------------------

    real( kind = core_rknd), dimension(:,:), allocatable :: &
      LES_tmp ! Temporary variable

    integer :: &
      i, &          ! Field index
      file_index, & ! index of the files
      min_ia, &     ! minimum ia value from all files
      max_iz        ! maximum iz value from all files

    logical :: l_internal_error

    type (input_field) :: &
      current_var ! The current variable

    type (stat_file), dimension(1:size(stat_files)) :: &
      fread_vars

    logical :: &
      l_convert_to_MKS ! convert inputs to MKS units

    !--------------------------- BEGIN CODE ---------------------------------

    l_error = .false.

      do file_index=1, size(stat_files)
        call open_netcdf_read( "U", stat_files(file_index),  &
                              fread_vars(file_index), l_internal_error )
        l_error = l_error .or. l_internal_error
      end do ! file_index=1, size(stat_files)

    if (l_error) then
      write(fstderr,*) "A fatal error occured while reading the input file."
      error stop "Fatal error in input_fields.get_input_variables_interp"
    end if

    ! initialize to boundary values to determine min and max later
    min_ia = 9999
    max_iz = -1

    do file_index=1, size(stat_files)
      ! Each column can have different levels within the shared file domain.
      call CLUBB_levels_within_LES_domain( fread_vars(file_index), gr%nzt, ngrdcol, gr%zt, & ! In
                                           k_lowest_zt(:,file_index), & ! Out
                                           k_highest_zt(:,file_index) ) ! Out

      call CLUBB_levels_within_LES_domain( fread_vars(file_index), gr%nzm, ngrdcol, gr%zm, & ! In
                                           k_lowest_zm(:,file_index), & ! Out
                                           k_highest_zm(:,file_index) ) ! Out

      ! update the min and max variables to allocate LES_tmp
      min_ia = min(min_ia, fread_vars(file_index)%ia)
      max_iz = max(max_iz, fread_vars(file_index)%iz)

    end do ! file_index=1, size(stat_files)

    ! Allocate the LES_tmp variables
    allocate( LES_tmp(1:size(stat_files), min_ia:max_iz) )

    do i=1, nvars

      current_var = variables(i)

      if( current_var%input_name == "none" .and. &
          current_var%l_input_var ) then
        write(fstderr,*) "The variable " //trim( current_var%clubb_name ) &
           //" is not setup for input_type " //trim( input_type )
        l_error = .true.
        cycle
      end if

      if( current_var%l_input_var ) then
        file_index = current_var%input_file_index

          l_convert_to_MKS = .false.
          call get_netcdf_var( fread_vars(file_index), current_var%input_name, timestep, &
                      l_convert_to_MKS,&
                      LES_tmp(file_index, fread_vars(file_index)%ia:fread_vars(file_index)%iz), &
                      l_internal_error )

        l_error = l_error .or. l_internal_error

        if( current_var%clubb_grid_type == "zt" ) then
          call inputfields_interp_and_adjust( gr, ngrdcol, current_var%clubb_grid_type, & ! In
                                              fread_vars(file_index), current_var%adjustment, & ! In
                                              LES_tmp(file_index, :), & ! In
                                              k_lowest_zt(:,file_index), & ! In
                                              k_highest_zt(:,file_index), & ! In
                                              current_var%clubb_var ) ! InOut
        else if( current_var%clubb_grid_type == "zm" ) then
          call inputfields_interp_and_adjust( gr, ngrdcol, current_var%clubb_grid_type, & ! In
                                              fread_vars(file_index), current_var%adjustment, & ! In
                                              LES_tmp(file_index, :), & ! In
                                              k_lowest_zm(:,file_index), & ! In
                                              k_highest_zm(:,file_index), & ! In
                                              current_var%clubb_var ) ! InOut
        end if ! current_var%clubb_grid_type == "zt"

      end if ! current_var%l_input_var
    end do ! i=1, nvars

      do file_index=1, size(stat_files)
        call close_netcdf_read( fread_vars(file_index) )
      end do ! file_index=1, size(stat_files)

    deallocate(LES_tmp)

  end subroutine get_input_variables_interp
        
  !----------------------------------------------------------------------------------------
  subroutine inputfields_interp_and_adjust( gr, ngrdcol, clubb_grid_type, & ! In
                                            fread_var, adjustment, & ! In
                                            LES_tmp, & ! In
                                            k_lowest_input, & ! In
                                            k_highest_input, & ! In
                                            clubb_var ) ! InOut

  ! Description:
  !   Interpolates an input field as needed and adjusts it for unit conversion.
  !
  !  References:
  !    None
  !----------------------------------------------------------------------------------------
    
    use interpolation, only: &
      lin_interpolate_two_points ! Procedure(s)

    use stat_file_utils, only: & 
      LES_grid_to_CLUBB_grid, & ! Procedure(s)
      CLUBB_levels_within_LES_domain

    use stat_file_module, only: &
      stat_file ! Type

    use grid_class, only: grid ! Type

    implicit none

    type (grid), target, intent(in) :: gr

    integer, intent(in) :: ngrdcol ! Number of columns [-]

    !--------------------------- Input Variables ---------------------------

    character(len = 3), intent(in) :: &
      clubb_grid_type ! The type of the grid, either "zm" or "zt"

    type (stat_file), intent(in) :: &
      fread_var

    real( kind = core_rknd ), intent(in) :: &
      adjustment ! The value used to adjust the variable for unit conversions

    real( kind = core_rknd ), dimension(fread_var%ia:fread_var%iz), intent(in) :: &
      LES_tmp ! The values that were read in from the file.

    integer, dimension(ngrdcol), intent(in) :: &
      k_lowest_input, &  ! The lowest CLUBB level that's within the LES domain.
      k_highest_input ! The highest CLUBB level that's within the LES domain.

    !--------------------------- Input/Output Variables ---------------------------
    
    real( kind = core_rknd ), dimension(:,:), intent(inout) :: &
      clubb_var ! The clubb variable to store the result in

    ! Local Variable(s)

    integer :: i, k ! Column and vertical indices
    integer :: exact_lev_idx, lower_lev_idx, upper_lev_idx ! File level indices [-]
    logical :: l_lin_int ! Whether interpolation is required [-]

    real( kind = core_rknd ), pointer, dimension(:,:) :: &
      grid_var ! All column heights on the requested grid [m]

    !--------------------------- Begin Code ---------------------------

    if ( clubb_grid_type == "zt" ) then
      grid_var => gr%zt
    else if ( clubb_grid_type == "zm" ) then
      grid_var => gr%zm
    end if

    do k = minval(k_lowest_input), maxval(k_highest_input)
      do i = 1, ngrdcol
        if ( k < k_lowest_input(i) .or. k > k_highest_input(i) ) cycle

        call LES_grid_to_CLUBB_grid( fread_var, grid_var(i,k), & ! In
                                     exact_lev_idx, lower_lev_idx, & ! Out
                                     upper_lev_idx, l_lin_int ) ! Out

        if ( l_lin_int ) then
          clubb_var(i,k) = lin_interpolate_two_points( grid_var(i,k), &
                            fread_var%z(upper_lev_idx), fread_var%z(lower_lev_idx), &
                            LES_tmp(upper_lev_idx), LES_tmp(lower_lev_idx) )
        else
          clubb_var(i,k) = LES_tmp(exact_lev_idx)
        end if

        clubb_var(i,k) = clubb_var(i,k) * adjustment
      end do
    end do

  end subroutine inputfields_interp_and_adjust

!-------------------------------------------------------------------------------
  subroutine inputfields_init( iunit, namelist_filename, day, month, year )

!-------------------------------------------------------------------------------
    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      namelist_filename

    integer, intent(in) :: iunit

    integer, intent(in) ::  & 
      day, month, year ! Start time the of simulation

    ! Local variables
    character(len=120) :: datafile

    ! Namelist definitions
    namelist /setfields/ datafile, input_type, & 
      l_input_um, l_input_vm, l_input_rtm, l_input_thlm, & 
      l_input_wp2, l_input_wprtp, l_input_wpthlp,  & 
      l_input_wp3, l_input_rtp2, l_input_rtp3, l_input_thlp2, l_input_thlp3,  &
      l_input_rtpthlp, l_input_upwp, l_input_vpwp, & 
      l_input_ug, l_input_vg, l_input_rcm,  & 
      l_input_wm_zt, l_input_exner, l_input_em, & 
      l_input_p, l_input_rho, l_input_rho_zm, & 
      l_input_Lscale, l_input_Lscale_up, l_input_Lscale_down, & 
      l_input_Kh_zt, l_input_Kh_zm, l_input_tau_zm, l_input_tau_zt, & 
      l_input_wpthvp, l_input_radht, &
      l_input_thl_1, l_input_thl_2, l_input_mixt_frac, l_input_chi_1, l_input_chi_2, &
      l_input_stdev_chi_1, l_input_stdev_chi_2, l_input_rc_1, l_input_rc_2, &
      l_input_thvm, l_input_rrm,l_input_Nrm,  & 
      l_input_rsm, l_input_rim, l_input_rgm,  &
      l_input_rrp2, l_input_Nrp2, l_input_wprrp, l_input_wpNrp, & 
      l_input_thlm_forcing, l_input_rtm_forcing, & 
      l_input_up2, l_input_vp2, l_input_sigma_sqd_w, l_input_Ncm,  & 
      l_input_Nccnm, l_input_Nim, l_input_Ngm, l_input_Nsm, &
      l_input_cloud_frac, l_input_sigma_sqd_w_zt, &
      l_input_veg_T_in_K, l_input_deep_soil_T_in_K, &
      l_input_sfc_soil_T_in_K, l_input_wprtp_forcing, l_input_wpthlp_forcing, &
      l_input_rtp2_forcing, l_input_thlp2_forcing, l_input_rtpthlp_forcing, l_input_thlprcp , &
      l_input_rcm_mc, l_input_rvm_mc, l_input_thlm_mc, l_input_wprtp_mc, &
      l_input_wpthlp_mc, l_input_rtp2_mc, l_input_thlp2_mc, l_input_rtpthlp_mc

    ! --- Begin Code ---

    clubb_day = day
    clubb_month = month
    clubb_year = year

    ! Pick some initial values
    datafile = ''
    input_type = 'clubb'

    l_input_um = .false.
    l_input_vm = .false.
    l_input_rtm = .false.
    l_input_thlm  = .false.
    l_input_wp2 = .false.
    l_input_wprtp = .false.
    l_input_wpthlp = .false.
    l_input_wp3 = .false.
    l_input_rtp2 = .false.
    l_input_rtp3 = .false.
    l_input_thlp2 = .false.
    l_input_thlp3 = .false.
    l_input_rtpthlp = .false.
    l_input_upwp = .false.
    l_input_vpwp = .false.
    l_input_ug = .false.
    l_input_vg = .false.
    l_input_rcm = .false.
    l_input_wm_zt = .false.
    l_input_exner = .false.
    l_input_em = .false.
    l_input_p = .false.
    l_input_rho = .false.
    l_input_Lscale = .false.
    l_input_Lscale_up = .false.
    l_input_Lscale_down = .false.
    l_input_Kh_zt = .false.
    l_input_Kh_zm = .false.
    l_input_tau_zm = .false.
    l_input_tau_zt = .false.
    l_input_wpthvp = .false.
    l_input_thl_1 = .false.
    l_input_thl_2 = .false.
    l_input_mixt_frac = .false.
    l_input_chi_1 = .false.
    l_input_chi_2 = .false.
    l_input_stdev_chi_1 = .false.
    l_input_stdev_chi_2 = .false.
    l_input_rc_1 = .false.
    l_input_rc_2 = .false.
    l_input_thvm = .false.
    l_input_rrm = .false.
    l_input_Nrm = .false.
    l_input_rsm = .false.
    l_input_rim = .false.
    l_input_rgm = .false.
    l_input_rrp2 = .false.
    l_input_Nrp2 = .false.
    l_input_wprrp = .false.
    l_input_wpNrp = .false.
    l_input_thlm_forcing = .false.
    l_input_rtm_forcing = .false.
    l_input_up2 = .false.
    l_input_vp2 = .false.
    l_input_sigma_sqd_w = .false.
    l_input_Ncm = .false.
    l_input_Nccnm = .false.
    l_input_Nim = .false.
    l_input_Nsm = .false.
    l_input_Ngm = .false.
    l_input_cloud_frac = .false.
    l_input_sigma_sqd_w_zt = .false.
    l_input_veg_T_in_K = .false.
    l_input_deep_soil_T_in_K = .false.
    l_input_sfc_soil_T_in_K  = .false.
    l_input_wprtp_forcing = .false.
    l_input_wpthlp_forcing = .false.
    l_input_rtp2_forcing = .false.
    l_input_thlp2_forcing = .false.
    l_input_rtpthlp_forcing = .false.
    l_input_thlprcp = .false.
    l_input_rcm_mc = .false.
    l_input_rvm_mc = .false.
    l_input_thlm_mc = .false.
    l_input_wprtp_mc = .false.
    l_input_wpthlp_mc = .false.
    l_input_rtp2_mc = .false.
    l_input_thlp2_mc = .false.
    l_input_rtpthlp_mc = .false.

    print *, namelist_filename
    ! Read in our namelist
    open(unit=iunit, file=namelist_filename, status='old', action='read')

    read(unit=iunit, nml=setfields)

    close(unit=iunit)

    ! Setup the netCDF file reader
    call set_filenames( datafile )    
      

    return
  end subroutine inputfields_init

  !-------------------------------------------------------------------------------
  subroutine cleanup_input_fields()

  ! Description:
  !   Cleans up any inputfields variables that were allocated during setup and now
  !   need to be deallocated. This must be called before CLUBB ends or there will
  !   be memory leaks.
  !
  !  References:
  !    None
  !-------------------------------------------------------------------------------

    implicit none

    deallocate(stat_files)

  end subroutine cleanup_input_fields
  !-----------------------------------------------------------------------

!===============================================================================

end module inputfields
