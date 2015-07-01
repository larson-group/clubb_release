!-----------------------------------------------------------------------
! $Id$

module inputfields
!  Description:
!    This file contains information and subroutines needed for restart
!    simulations and simulations where the variables are reset at every
!    timestep, called a "clubb_inputfields" simulation.
!===============================================================================

  use clubb_precision, only: &
    core_rknd ! variable(s)

  implicit none

  ! Run information
  character(len=120), dimension(:), allocatable, public :: & 
    stat_files ! a list of all files used for input

  character(len=10), public :: input_type

  logical, public :: &
    l_input_um = .false., l_input_vm = .false., l_input_rtm = .false., l_input_thlm = .false., & 
    l_input_wp2 = .false., l_input_wprtp = .false., l_input_wpthlp = .false.,  & 
    l_input_wp3 = .false., l_input_rtp2 = .false., l_input_thlp2 = .false.,  & 
    l_input_rtpthlp = .false., l_input_upwp = .false., l_input_vpwp = .false., & 
    l_input_ug = .false., l_input_vg = .false., l_input_rcm = .false.,  & 
    l_input_wm_zt = .false., l_input_exner = .false., l_input_em = .false., & 
    l_input_p = .false., l_input_rho = .false., l_input_rho_zm = .false., & 
    l_input_rho_ds_zm = .false., l_input_rho_ds_zt = .false., &
    l_input_thv_ds_zm = .false., l_input_thv_ds_zt = .false., &
    l_input_Lscale = .false., l_input_Lscale_up = .false., l_input_Lscale_down = .false., & 
    l_input_Kh_zt = .false., l_input_Kh_zm = .false., &
    l_input_tau_zm = .false., l_input_tau_zt = .false., & 
    l_input_wpthvp = .false., l_input_radht = .false., &
    l_input_thl_1 = .false., l_input_thl_2 = .false., l_input_mixt_frac = .false., &
    l_input_chi_1 = .false., l_input_chi_2 = .false., &
    l_input_stdev_chi_1 = .false., l_input_stdev_chi_2 = .false., &
    l_input_rc_1 = .false., l_input_rc_2 = .false., &
    l_input_thvm = .false., l_input_rrm = .false., &
    l_input_Nrm = .false.,  l_input_Ncm = .false.,  & 
    l_input_rsm = .false., l_input_rim = .false., &
    l_input_Nsm = .false., l_input_Ngm = .false., &
    l_input_rgm = .false., l_input_Nccnm = .false., l_input_Nim = .false., &
    l_input_rrp2 = .false., l_input_Nrp2 = .false., &
    l_input_wprrp = .false., l_input_wpNrp = .false., &
    l_input_precip_frac = .false., &
    l_input_rtprrp = .false., l_input_thlprrp = .false., &
    l_input_rtpNrp = .false., l_input_thlpNrp = .false., &
    l_input_rrpNrp = .false., &
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
    num_sam_inputfields = 74, & ! The number of input fields for SAM
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

    real( kind = core_rknd ), dimension(:), pointer :: &
      clubb_var  ! The clubb variable to store the result in

    real( kind = core_rknd ) :: &
      adjustment ! The SAM variable will be multiplied by this amount
                 ! to convert the variable to the correct units.

  end type input_field

  contains

!===============================================================================


!-----------------------------------------------------------------------
  subroutine set_filenames( file_prefix )

! Description: Set the names of the GrADS files to be used.
!   Used by clubb_inputfields and clubb_restart.
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only: fstderr, fstdout ! Constants

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variables
    character(len=*), intent(in) :: &
      file_prefix

    ! Local Variables
    logical :: l_grads

    ! ---- Begin Code ----

    stats_input_type = -1

    select case ( trim( input_type ) )

    case ( "coamps_les", "COAMPS_LES", "coamps" )
      allocate(stat_files(1:2))
      stat_files(coamps_sm) = trim( file_prefix )//"_coamps_sm.ctl"
      stat_files(coamps_sw) = trim( file_prefix )//"_coamps_sw.ctl"

      stats_input_type = coamps_input_type

    case ( "clubb_scm", "CLUBB_SCM", "clubb" )
      allocate(stat_files(1:3))
      stat_files(clubb_zt) = trim( file_prefix )//"_zt.ctl"
      stat_files(clubb_zm) = trim( file_prefix )//"_zm.ctl"
      stat_files(clubb_sfc) = trim( file_prefix )//"_sfc.ctl"

      stats_input_type = clubb_input_type

      inquire(file=stat_files(clubb_zt),exist=l_grads)

      if ( .not. l_grads ) then
        write(fstdout,*) "inputfields: Cannot find GrADS ctl file, assuming netCDF input"
        stat_files(clubb_zt) = trim( file_prefix )// "_zt.nc"
        stat_files(clubb_zm) = trim( file_prefix )// "_zm.nc"
        stat_files(clubb_sfc) = trim( file_prefix )// "_sfc.nc"
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
      stat_files(rams_file) = trim( file_prefix )//"_rams.ctl"

      stats_input_type = rams_input_type

    case default
      write(fstderr,*) "Don't know how to handle input_type = "// & 
        input_type
      stop

    end select

    return
  end subroutine set_filenames
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine stat_fields_reader( timestep, hydromet_pdf_params, &
                                 rtphmp_zt, thlphmp_zt, hmxphmyp_zt )
! Description:
!   Reads in variables for the model from statistical data

! References:
!   None
!-----------------------------------------------------------------------

    use variables_prognostic_module, only: & 
        um,  & ! Variable(s)
        vm, & 
        rtm, & 
        thlm, & 
        wp2, & 
        wp3, & 
        wprtp, & 
        wpthlp, & 
        rtp2, & 
        thlp2, & 
        rtpthlp, & 
        upwp, & 
        vpwp, & 
        p_in_Pa, & 
        exner, & 
        rcm, & 
        wm_zt, & 
        rho, & 
        rho_zm, & 
        rho_ds_zm, &
        rho_ds_zt, &
        thv_ds_zm, &
        thv_ds_zt, &
        thlm_forcing, & 
        rtm_forcing, & 
        cloud_frac, & 
        up2, & 
        vp2, & 
        sigma_sqd_w, &
        wprtp_forcing, &
        wpthlp_forcing, &
        rtp2_forcing, &
        thlp2_forcing, &
        rtpthlp_forcing

    use variables_diagnostic_module, only: & 
        hydromet, & ! Variable(s)
        hydrometp2, &
        wphydrometp, &
        Ncm, &
        tau_zm, &
        tau_zt, & 
        ug, & 
        vg, & 
        Lscale, & 
        Lscale_up, & 
        Lscale_down, & 
        Kh_zm, & 
        Kh_zt, &
        thvm, &
        wpthvp, &
        Nccnm, & 
        sigma_sqd_w_zt, & 
        em, &
        radht, &
        !Frad, &
        thlprcp

    use variables_prognostic_module, only: & 
        pdf_params ! Variable(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter

    use grid_class, only: & 
        gr,  & ! Variable(s)
        zt2zm ! Procedure(s)

    use constants_clubb, only:  &
        rt_tol,    & ! Variable(s)
        thl_tol,   &
        w_tol_sqd, &
        em_min,     &
        fstderr, &
        Cp, &
        Lv, &
        pascal_per_mb, &
        g_per_kg, &
        sec_per_day

    use array_index, only:  & 
        iirrm, iiNrm, iirsm, iirim, iirgm, iiNim, iiNgm, iiNsm

    use stat_file_utils, only: & 
        LES_grid_to_CLUBB_grid, & ! Procedure(s)
        CLUBB_levels_within_LES_domain

    use extrapolation, only: &
        lin_ext_zt_bottom, &
        lin_ext_zm_bottom

    use parameters_microphys, only: &
        microphys_scheme, & ! Variable(s)
        l_predict_Nc

    use parameters_model, only: &
        hydromet_dim

    use soil_vegetation, only: deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K

    use clubb_precision, only: &
        core_rknd ! Variable(s)


    implicit none

    ! External
    intrinsic :: max, trim, any

    ! Arguments
    integer, intent(in) :: timestep

    type(hydromet_pdf_parameter), dimension(gr%nz), intent(inout) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters      [units vary]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(inout) :: &
      rtphmp_zt,  &
      thlphmp_zt

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim,hydromet_dim), &
    intent(inout) :: &
      hmxphmyp_zt

    ! Local Variables
    logical :: l_read_error, l_fatal_error

    real( kind = core_rknd ), dimension(gr%nz+1) :: tmp1

    integer, dimension(1:size(stat_files)) ::  &
      k_lowest_zt, &  ! The lowest CLUBB thermodynamic level that's within the LES domain.
      k_highest_zt, & ! The highest CLUBB thermodynamic level that's within the LES domain.
      k_lowest_zm, &  ! The lowest CLUBB momentum level that's within the LES domain.
      k_highest_zm    ! The highest CLUBB momentum level that's within the LES domain.

    integer :: k  ! Array index
 
    real( kind = core_rknd), dimension(gr%nz), target :: &
      temp_Nrm, temp_Ncm, temp_rgm, temp_rim, temp_Ngm, temp_Nsm, &
      temp_rrm, temp_rsm, temp_tke, temp_wpup_sgs, temp_wpvp_sgs, &
      temp_Nim, temp_rrp2, temp_Nrp2, temp_wprrp, temp_wpNrp, &
      temp_precip_frac, temp_rtprrp, temp_rtpNrp, temp_thlprrp, &
      temp_thlpNrp, temp_rrpNrp ! temp. variables

    type (input_field), dimension(:), allocatable :: &
      SAM_variables, & ! A list of SAM variables to read in.
      coamps_variables, & ! A list of coamps variables to read in.
      RAMS_variables ! A list of RAMS LES variables to read in.

    ! ---- Begin Code ----

    select case ( stats_input_type )


    !-------------------------------------
    ! CLUBB stats data
    !-------------------------------------
    case ( clubb_input_type )

      !  Thermo grid - zt file

      ! Initialize l_fatal_error for case clubb_input_type
      l_fatal_error = .false.

      call get_clubb_variable_interpolated &
           ( l_input_um, stat_files(clubb_zt), "um", gr%nz, timestep, &
             gr%zt, um, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vm, stat_files(clubb_zt), "vm", gr%nz, timestep, &
             gr%zt, vm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtm, stat_files(clubb_zt), "rtm", gr%nz, timestep, &
             gr%zt, rtm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_thlm, stat_files(clubb_zt), "thlm", gr%nz, timestep, &
             gr%zt, thlm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wp3, stat_files(clubb_zt), "wp3", gr%nz, timestep, &
             gr%zt, wp3, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_tau_zt, stat_files(clubb_zt), "tau_zt", gr%nz, timestep, &
             gr%zt, tau_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rrm, stat_files(clubb_zt), "rrm", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_rrm ) then
        hydromet(1:gr%nz,iirrm) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rsm, stat_files(clubb_zt), "rsm", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_rsm ) then
        hydromet(1:gr%nz,iirsm) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rim, stat_files(clubb_zt), "rim", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_rim ) then
        hydromet(1:gr%nz,iirim) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rgm, stat_files(clubb_zt), "rgm", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_rgm ) then
        hydromet(1:gr%nz,iirgm) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

!--------------------------------------------------------
! Added variables for clubb_restart
      call get_clubb_variable_interpolated &
           ( l_input_p, stat_files(clubb_zt), "p_in_Pa", gr%nz, timestep, &
             gr%zt, p_in_Pa, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_exner, stat_files(clubb_zt), "exner", gr%nz, timestep, &
             gr%zt, exner, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_ug, stat_files(clubb_zt), "ug", gr%nz, timestep, &
             gr%zt, ug, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vg, stat_files(clubb_zt), "vg", gr%nz, timestep, &
             gr%zt, vg, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rcm, stat_files(clubb_zt), "rcm", gr%nz, timestep, &
             gr%zt, rcm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wm_zt, stat_files(clubb_zt), "wm", gr%nz, timestep, &
             gr%zt, wm_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_rho, stat_files(clubb_zt), "rho", gr%nz, timestep, &
             gr%zt, rho, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rho_ds_zt, stat_files(clubb_zt), "rho_ds_zt", gr%nz, timestep, &
             gr%zt, rho_ds_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thv_ds_zt, stat_files(clubb_zt), "thv_ds_zt", gr%nz, timestep, &
             gr%zt, thv_ds_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Lscale, stat_files(clubb_zt), "Lscale", gr%nz, timestep, &
             gr%zt, Lscale, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Lscale_up, stat_files(clubb_zt), "Lscale_up", gr%nz, timestep, &
             gr%zt, Lscale_up, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Lscale_down, stat_files(clubb_zt), "Lscale_down", gr%nz, timestep, &
             gr%zt, Lscale_down, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Kh_zt, stat_files(clubb_zt), "Kh_zt", gr%nz, timestep, &
             gr%zt, Kh_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thvm, stat_files(clubb_zt), "thvm", gr%nz, timestep, &
             gr%zt, thvm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thlm_forcing, stat_files(clubb_zt), "thlm_forcing", gr%nz, timestep, &
             gr%zt, thlm_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtm_forcing, stat_files(clubb_zt), "rtm_forcing", gr%nz, timestep, &
             gr%zt, rtm_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Ncm, stat_files(clubb_zt), "Ncm", gr%nz, timestep, &
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_Ncm ) then
         Ncm(1:gr%nz) = tmp1(1:gr%nz)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Nccnm, stat_files(clubb_zt), "Nccnm", gr%nz, timestep, &
             gr%zt, Nccnm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Nim, stat_files(clubb_zt), "Nim", gr%nz, timestep, &
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_Nim ) then
        hydromet(1:gr%nz, iiNim) = tmp1(1:gr%nz)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_cloud_frac, stat_files(clubb_zt), "cloud_frac", gr%nz, timestep, &
             gr%zt, cloud_frac, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Nrm, stat_files(clubb_zt), "Nrm", gr%nz, timestep, &
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_Nrm ) then
        hydromet(1:gr%nz, iiNrm) = tmp1(1:gr%nz)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_sigma_sqd_w_zt, stat_files(clubb_zt), "sigma_sqd_w_zt", gr%nz, timestep, &
             gr%zt, sigma_sqd_w_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_radht, stat_files(clubb_zt), "radht", gr%nz, timestep, &
             gr%zt, radht, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      ! PDF Parameters (needed for K&K microphysics)
      call get_clubb_variable_interpolated &
           ( l_input_thl_1, stat_files(clubb_zt), "thl_1", gr%nz, timestep, &
             gr%zt, pdf_params%thl_1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thl_2, stat_files(clubb_zt), "thl_2", gr%nz, timestep, &
             gr%zt, pdf_params%thl_2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_mixt_frac, stat_files(clubb_zt), "mixt_frac", gr%nz, timestep, &
             gr%zt, pdf_params%mixt_frac, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_chi_1, stat_files(clubb_zt), "chi_1", gr%nz, timestep, &
             gr%zt, pdf_params%chi_1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_chi_2, stat_files(clubb_zt), "chi_2", gr%nz, timestep, &
             gr%zt, pdf_params%chi_2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_stdev_chi_1, stat_files(clubb_zt), "stdev_chi_1", gr%nz, timestep, &
             gr%zt, pdf_params%stdev_chi_1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_stdev_chi_2, stat_files(clubb_zt), "stdev_chi_2", gr%nz, timestep, &
             gr%zt, pdf_params%stdev_chi_2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rc_1, stat_files(clubb_zt), "rc_1", gr%nz, timestep, &
             gr%zt, pdf_params%rc_1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rc_2, stat_files(clubb_zt), "rc_2", gr%nz, timestep, &
             gr%zt, pdf_params%rc_2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error




!--------------------------------------------------------

      ! Read in the zm file

      call get_clubb_variable_interpolated &
           ( l_input_wp2, stat_files(clubb_zm), "wp2", gr%nz, timestep, &
             gr%zm, wp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wprtp, stat_files(clubb_zm), "wprtp", gr%nz, timestep, &
             gr%zm, wprtp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wpthlp, stat_files(clubb_zm), "wpthlp", gr%nz, timestep, &
             gr%zm, wpthlp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wpthvp, stat_files(clubb_zm), "wpthvp", gr%nz, timestep, &
             gr%zm, wpthvp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtp2, stat_files(clubb_zm), "rtp2", gr%nz, timestep, &
             gr%zm, rtp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thlp2, stat_files(clubb_zm), "thlp2", gr%nz, timestep, &
             gr%zm, thlp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtpthlp, stat_files(clubb_zm), "rtpthlp", gr%nz, timestep, &
             gr%zm, rtpthlp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_upwp, stat_files(clubb_zm), "upwp", gr%nz, timestep, &
             gr%zm, upwp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vpwp, stat_files(clubb_zm), "vpwp", gr%nz, timestep, &
             gr%zm, vpwp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

!-----------------------------------------------------------
      call get_clubb_variable_interpolated &
           ( l_input_em, stat_files(clubb_zm), "em", gr%nz, timestep, &
             gr%zm, em, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rho_zm, stat_files(clubb_zm), "rho_zm", gr%nz, timestep, &
             gr%zm, rho_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rho_ds_zm, stat_files(clubb_zm), "rho_ds_zm", gr%nz, timestep, &
             gr%zm, rho_ds_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thv_ds_zm, stat_files(clubb_zm), "thv_ds_zm", gr%nz, timestep, &
             gr%zm, thv_ds_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Kh_zm, stat_files(clubb_zm), "Kh_zm", gr%nz, timestep, &
             gr%zm, Kh_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_tau_zm, stat_files(clubb_zm), "tau_zm", gr%nz, timestep, &
             gr%zm, tau_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_up2, stat_files(clubb_zm), "up2", gr%nz, timestep, &
             gr%zm, up2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vp2, stat_files(clubb_zm), "vp2", gr%nz, timestep, &
             gr%zm, vp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_sigma_sqd_w, stat_files(clubb_zm), "sigma_sqd_w", gr%nz, timestep, &
             gr%zm, sigma_sqd_w, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wprtp_forcing, stat_files(clubb_zm), "wprtp_forcing", gr%nz, timestep, &
             gr%zm, wprtp_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_wpthlp_forcing, stat_files(clubb_zm), "wpthlp_forcing", gr%nz, timestep, &
             gr%zm, wpthlp_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_rtp2_forcing, stat_files(clubb_zm), "rtp2_forcing", gr%nz, timestep, &
             gr%zm, rtp2_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_thlp2_forcing, stat_files(clubb_zm), "thlp2_forcing", gr%nz, timestep, &
             gr%zm, thlp2_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtpthlp_forcing, stat_files(clubb_zm), "rtpthlp_forcing", gr%nz, timestep, &
             gr%zm, rtpthlp_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_thlprcp, stat_files(clubb_zm), "thlprcp", gr%nz, timestep, &
             gr%zm, thlprcp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


!-----------------------------------------------------------

      call get_clubb_variable_interpolated &
           ( l_input_veg_T_in_K, stat_files(clubb_sfc), "veg_T_in_K", 1, timestep, &
             (/0._core_rknd/), tmp1(1), l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error
      veg_T_in_K = tmp1(1)

      call get_clubb_variable_interpolated &
           ( l_input_deep_soil_T_in_K, stat_files(clubb_sfc), "deep_soil_T_in_K", 1, timestep, &
             (/0._core_rknd/), tmp1(1), l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error
      deep_soil_T_in_K = tmp1(1)

      call get_clubb_variable_interpolated &
           ( l_input_sfc_soil_T_in_K, stat_files(clubb_sfc), "sfc_soil_T_in_K", 1, timestep, &
             (/0._core_rknd/), tmp1(1), l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error
      sfc_soil_T_in_K = tmp1(1)



      if ( l_fatal_error ) stop "oops, get_grads_var failed in stat_fields_reader"


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


      coamps_variables(k)%l_input_var = l_input_Lscale
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "Lscale"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_Lscale_up
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "Lscale_up"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_Lscale_down
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "Lscale_down"

      k = k + 1

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
        if ( iirrm < 1 ) then
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
        if ( iiNrm < 1 ) then
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
        if ( iirsm < 1 ) then
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
        if ( iirim < 1 ) then
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
        if ( iirgm < 1 ) then
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
        if ( iiNim < 1 ) then
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

      if ( l_fatal_error ) stop "oops, bad input for inputfields"
      
      call get_input_variables_interp &
             (k, coamps_variables, timestep, &
              k_lowest_zt, k_highest_zt, &
              k_lowest_zm, k_highest_zm, l_fatal_error )

      if ( l_input_um ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, use the
        ! values of um at thermodynamic levels 3 and 2 to find the value at
        ! thermodynamic level 1 through the use of a linear extension.
        if ( k_lowest_zt(coamps_sm) == 2 ) then
          um(1)  & 
          = lin_ext_zt_bottom( um(3), um(2), & 
                               gr%zt(3), gr%zt(2), gr%zt(1) )

        endif
      endif

      if ( l_input_vm ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, use the
        ! values of vm at thermodynamic levels 3 and 2 to find the value at
        ! thermodynamic level 1 through the use of a linear extension.
        if ( k_lowest_zt(coamps_sm) == 2 ) then
          vm(1) = lin_ext_zt_bottom( vm(3), vm(2), & 
                                     gr%zt(3), gr%zt(2), gr%zt(1) )
        endif
      endif

      if ( l_input_rtm ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of rtm at thermodynamic level 1 to the value at thermodynamic
        ! level 2, as it is done in advance_xm_wpxp.
        if ( k_lowest_zt(coamps_sm) == 2 ) then
          rtm(1) = rtm(2)
        endif
      endif

      if ( l_input_thlm ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of thlm at thermodynamic level 1 to the value at thermodynamic
        ! level 2, as it is done in advance_xm_wpxp.
        if ( k_lowest_zt(coamps_sm) == 2 ) then
          thlm(1) = thlm(2)
        endif
      endif

      if ( l_input_wp3 ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of wp3 at thermodynamic level 1 to 0, as it is done in
        ! advance_wp2_wp3.
        if ( k_lowest_zt(coamps_sm) == 2 ) then
          wp3(1) = 0.0_core_rknd
        endif
      endif

      if ( l_input_rrm ) then
        hydromet(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm),iirrm) = &
                    temp_rrm(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm))
      end if

      if ( l_input_Nrm ) then
        hydromet(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm),iiNrm) = &
                    temp_Nrm(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm))
      end if

      if ( l_input_Ncm ) then
        Ncm(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm)) = &
                    temp_Ncm(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm))
      end if

      if ( l_input_rsm ) then
        hydromet(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm),iirsm) = &
                    temp_rsm(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm))
      end if

      if ( l_input_rim ) then
        hydromet(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm),iirim) = &
                    temp_rim(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm))
      end if

      if ( l_input_rgm ) then
        hydromet(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm),iirgm) = &
                    temp_rgm(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm))
      end if

      if ( l_input_Nim ) then
        hydromet(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm),iiNim) = &
                    temp_Nim(k_lowest_zt(coamps_sm):k_highest_zt(coamps_sm))
      end if

      if ( l_input_wprtp ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wprtp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like surface_varnce.
        if ( k_lowest_zm(coamps_sm) == 2 ) then
          wprtp(1) = lin_ext_zm_bottom( wprtp(3), wprtp(2), & 
                                        gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if ( l_input_wpthlp ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wpthlp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like surface_varnce.
        if ( k_lowest_zm(coamps_sm) == 2 ) then
          wpthlp(1) = lin_ext_zm_bottom( wpthlp(3), wpthlp(2), & 
                                         gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if ( l_input_rtp2 ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, set the value of rtp2
        ! at momentum level 1 to the value at momentum level 2.
        ! Using a linear extension here resulted in negatives.
        if ( k_lowest_zm(coamps_sm) == 2 ) then
          rtp2(1) =  rtp2(2)
        endif

        where ( rtp2 < rt_tol**2 )
          rtp2 = rt_tol**2
        end where
      endif

      if ( l_input_thlp2 ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, set the value of thlp2
        ! at momentum level 1 to the value at momentum level 2.
        ! Using a linear extension here resulted in negatives.
        if ( k_lowest_zm(coamps_sm) == 2 ) then
          thlp2(1) = thlp2(2)
        endif
        
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
        if ( k_lowest_zm(coamps_sm) == 2 ) then
          rtpthlp(1) = lin_ext_zm_bottom( rtpthlp(3), rtpthlp(2), & 
                                          gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if( l_input_em ) then
        do k=k_lowest_zm(coamps_sm), k_highest_zm(coamps_sm), 1
          em(k) = em(k) + temp_tke(k)
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
        do k=k_lowest_zm(coamps_sw), k_highest_zm(coamps_sw), 1
          upwp(k) = temp_wpup_sgs(k)
        end do
      endif

      if ( l_input_vpwp) then
        do k=k_lowest_zm(coamps_sw), k_highest_zm(coamps_sw), 1
          vpwp(k) = temp_wpvp_sgs(k)
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
      SAM_variables(k)%input_name = "QT"
      SAM_variables(k)%clubb_var => rtm
      SAM_variables(k)%adjustment = 1.0e-3_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thlm
      SAM_variables(k)%input_name = "THETAL"
      SAM_variables(k)%clubb_var => thlm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_wp2
      SAM_variables(k)%input_name = "W2"
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
      SAM_variables(k)%input_name = "QTFLUX"
      SAM_variables(k)%clubb_var => wprtp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      ! Note that this needs to be adjusted by 1/(RHO * Cp)
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_wpthlp
      SAM_variables(k)%input_name = "TLFLUX"
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
      SAM_variables(k)%input_name = "QT2"
      SAM_variables(k)%clubb_var => rtp2
      SAM_variables(k)%adjustment = 1.0e-6_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thlp2
      SAM_variables(k)%input_name = "THEL2"
      SAM_variables(k)%clubb_var => thlp2
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
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

      SAM_variables(k)%l_input_var = l_input_Lscale
      SAM_variables(k)%clubb_name = "Lscale"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_Lscale_up
      SAM_variables(k)%clubb_name = "Lscale_up"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_Lscale_down
      SAM_variables(k)%clubb_name = "Lscale_down"
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
      !SAM_variables(k)%input_name = "NR" ! SAM Morrison microphysics
      SAM_variables(k)%input_name = "CONP" ! SAM KK microphysics
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
      SAM_variables(k)%input_name = "RRP2"
      SAM_variables(k)%clubb_var => temp_rrp2
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_Nrp2 = 0.0_core_rknd ! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_Nrp2
      SAM_variables(k)%input_name = "NRP2"
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

      temp_precip_frac = 0.0_core_rknd

      SAM_variables(k)%l_input_var = l_input_precip_frac
      SAM_variables(k)%input_name = "PREC_FRAC"
      SAM_variables(k)%clubb_var => temp_precip_frac
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_rtprrp = 0.0_core_rknd

      SAM_variables(k)%l_input_var = l_input_rtprrp
      SAM_variables(k)%input_name = "RTPRRP"
      SAM_variables(k)%clubb_var => temp_rtprrp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_rtpNrp = 0.0_core_rknd

      SAM_variables(k)%l_input_var = l_input_rtpNrp
      SAM_variables(k)%input_name = "RTPNRP"
      SAM_variables(k)%clubb_var => temp_rtpNrp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_thlprrp = 0.0_core_rknd

      SAM_variables(k)%l_input_var = l_input_thlprrp
      SAM_variables(k)%input_name = "THLPRRP"
      SAM_variables(k)%clubb_var => temp_thlprrp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_thlpNrp = 0.0_core_rknd

      SAM_variables(k)%l_input_var = l_input_thlpNrp
      SAM_variables(k)%input_name = "THLPNRP"
      SAM_variables(k)%clubb_var => temp_thlpNrp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      temp_rrpNrp = 0.0_core_rknd

      SAM_variables(k)%l_input_var = l_input_rrpNrp
      SAM_variables(k)%input_name = "RRPNRP"
      SAM_variables(k)%clubb_var => temp_rrpNrp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zt"
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
      SAM_variables(k)%clubb_var => up2
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file

      k = k + 1

      ! Note that V2 is on the thermodynamic grid in SAM while vp2 is
      ! on the momentum grid in CLUBB. This will have to be interpolated.
      SAM_variables(k)%l_input_var = l_input_vp2
      SAM_variables(k)%input_name = "V2"
      SAM_variables(k)%clubb_var => vp2
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
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
      SAM_variables(k)%input_name = "TQ"
      SAM_variables(k)%clubb_var => rtpthlp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%clubb_grid_type = "zm"
      SAM_variables(k)%input_file_index = sam_file


      call get_input_variables_interp &
             (k, SAM_variables, timestep, &
              k_lowest_zt, k_highest_zt, &
              k_lowest_zm, k_highest_zm, l_fatal_error )

      if( l_fatal_error ) then
        STOP "Error reading SAM input fields"
      end if

      ! Perform non-constant adjustments
      do k=k_lowest_zt(sam_file), k_highest_zt(sam_file)
        if(l_input_Nrm) then
          ! Nrm = NR * (1.0e6 / RHO)
          !temp_Nrm(k) = temp_Nrm(k) * (1.0e6_core_rknd / rho(k))
          temp_Nrm(k) = temp_Nrm(k) * 1.0e6_core_rknd
        end if

        if(l_input_Ncm) then
          ! Ncm = NC * (1.0e6 / RHO)
          temp_Ncm(k) = temp_Ncm(k) * (1.0e6_core_rknd / rho(k))
        end if

        if(l_input_Nim) then
          ! Nim = NI * (1.0e6 / RHO)
          temp_Nim(k) = temp_Nim(k) * (1.0e7_core_rknd / rho(k))
        end if

        if(l_input_Nsm) then
          ! Nsm = NI * (1.0e6 / RHO)
          temp_Nsm(k) = temp_Nsm(k) * (1.0e7_core_rknd / rho(k))
        end if

        if(l_input_Ngm) then
          ! Ngm = NG * (1.0e6 / RHO)
          temp_Ngm(k) = temp_Ngm(k) * (1.0e7_core_rknd / rho(k))
        end if
 
      end do

      do k=k_lowest_zm(sam_file), k_highest_zm(sam_file)
        if(l_input_wprtp) then
          ! wprtp = QTFLUX / (RHO*Lv)
          wprtp(k) = wprtp(k) / (rho(k)*Lv)
        end if

        if(l_input_wpthlp) then
          ! wpthlp = TLFLUX / (RHO*Cp)
          wpthlp(k) = wpthlp(k) / (rho(k)*Cp)
        end if
        
        ! Use linear interpolation to convert up2 and vp2 from the thermodynamic
        ! grid to the momentum grid. See CLUBB Equations 170 and 171.
        if( k < gr%nz) then
          if(l_input_up2) then
            up2(k) = ( (up2(k+1) - up2(k))/(gr%zt(k+1) - gr%zt(k)) ) &
                     * (gr%zm(k) - gr%zt(k)) + up2(k)
          end if

          if(l_input_vp2) then
            vp2(k) = ( (vp2(k+1) - vp2(k))/(gr%zt(k+1) - gr%zt(k)) ) &
                     * (gr%zm(k) - gr%zt(k)) + vp2(k)
          end if
        else ! k = gr%nz
          if(l_input_up2) then
            up2(k) = ( (up2(k) - up2(k-1))/(gr%zt(k) - gr%zt(k-1)) ) &
                     * (gr%zm(k) - gr%zt(k)) + up2(k)
          end if

          if(l_input_vp2) then
            vp2(k) = ( (vp2(k) - vp2(k-1))/(gr%zt(k) - gr%zt(k-1)) ) &
                     * (gr%zm(k) - gr%zt(k)) + vp2(k)
          end if
        end if
      end do
       
      ! Add hydrometeor variables.
      if( l_input_Nrm .and. iiNrm > 0) then
        hydromet(k_lowest_zt(sam_file):k_highest_zt(sam_file),iiNrm) = &
                   temp_Nrm(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      end if

      if ( l_input_Ncm .and. l_predict_Nc ) then
        Ncm(k_lowest_zt(sam_file):k_highest_zt(sam_file)) = &
                   temp_Ncm(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      end if

      if( l_input_Nim .and. iiNim > 0) then
        hydromet(k_lowest_zt(sam_file):k_highest_zt(sam_file),iiNim) = &
                 temp_Nim(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      end if

      if( l_input_rrm .and. iirrm > 0) then
        hydromet(k_lowest_zt(sam_file):k_highest_zt(sam_file),iirrm) = &
                   temp_rrm(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      end if

      if( l_input_Ngm .and. iiNgm > 0) then
        hydromet(k_lowest_zt(sam_file):k_highest_zt(sam_file),iiNgm) = &
                 temp_Ngm(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      end if

      if( l_input_rgm .and. iirgm > 0) then
        hydromet(k_lowest_zt(sam_file):k_highest_zt(sam_file),iirgm) = &
                   temp_rgm(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      end if

      if( l_input_rim .and. iirim > 0) then
        hydromet(k_lowest_zt(sam_file):k_highest_zt(sam_file),iirim) = &
                   temp_rim(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      end if

      if( l_input_Nsm .and. iiNsm > 0) then
        hydromet(k_lowest_zt(sam_file):k_highest_zt(sam_file),iiNsm) = &
                 temp_Nsm(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      end if

      if( l_input_rsm .and. iirsm > 0) then
        hydromet(k_lowest_zt(sam_file):k_highest_zt(sam_file),iirsm) = &
                   temp_rsm(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      end if

      if ( l_input_rrp2 .and. iirrm > 0 ) then
         hydrometp2(k_lowest_zt(sam_file):k_highest_zt(sam_file),iirrm) &
         = temp_rrp2(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      endif

      if ( l_input_Nrp2 .and. iiNrm > 0 ) then
         hydrometp2(k_lowest_zt(sam_file):k_highest_zt(sam_file),iiNrm) &
         = temp_Nrp2(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      endif

      if ( l_input_wprrp .and. iirrm > 0 ) then
         wphydrometp(k_lowest_zt(sam_file):k_highest_zt(sam_file),iirrm) &
         = temp_wprrp(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      endif

      if ( l_input_wpNrp .and. iiNrm > 0 ) then
         wphydrometp(k_lowest_zt(sam_file):k_highest_zt(sam_file),iiNrm) &
         = temp_wpNrp(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      endif

      if ( l_input_rtprrp .and. iirrm > 0 ) then
         rtphmp_zt(k_lowest_zt(sam_file):k_highest_zt(sam_file),iirrm) &
         = temp_rtprrp(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      endif

      if ( l_input_rtpNrp .and. iiNrm > 0 ) then
         rtphmp_zt(k_lowest_zt(sam_file):k_highest_zt(sam_file),iiNrm) &
         = temp_rtpNrp(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      endif

      if ( l_input_thlprrp .and. iirrm > 0 ) then
         thlphmp_zt(k_lowest_zt(sam_file):k_highest_zt(sam_file),iirrm) &
         = temp_thlprrp(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      endif

      if ( l_input_thlpNrp .and. iiNrm > 0 ) then
         thlphmp_zt(k_lowest_zt(sam_file):k_highest_zt(sam_file),iiNrm) &
         = temp_thlpNrp(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      endif

      if ( l_input_rrpNrp .and. iirrm > 0 .and. iiNrm > 0 ) then
         hmxphmyp_zt(k_lowest_zt(sam_file):k_highest_zt(sam_file),iirrm,iiNrm) &
         = temp_rrpNrp(k_lowest_zt(sam_file):k_highest_zt(sam_file))
         hmxphmyp_zt(k_lowest_zt(sam_file):k_highest_zt(sam_file),iiNrm,iirrm) &
         = temp_rrpNrp(k_lowest_zt(sam_file):k_highest_zt(sam_file))
      endif

      hydromet_pdf_params(k_lowest_zt(sam_file):k_highest_zt(sam_file))%precip_frac &
      = temp_precip_frac(k_lowest_zt(sam_file):k_highest_zt(sam_file))

      if ( l_fatal_error ) stop "Failed to read inputfields for SAM."

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
             (k, RAMS_variables, timestep, &
              k_lowest_zt, k_highest_zt, &
              k_lowest_zm, k_highest_zm, l_fatal_error )

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
  subroutine compute_timestep( iunit, filename, l_restart, & 
                               time, nearest_timestep )

    ! Description:
    ! Given a time 'time', determines the closest output time in a GrADS (or
    ! NetCDF) file.

    !-----------------------------------------------------------------------

    use stat_file_module, only: & 
        stat_file     ! Type

    use stat_file_utils, only: &
        l_netcdf_file

    use input_grads, only: &
        open_grads_read,  & ! Procedure(s)
        close_grads_read

    use constants_clubb, only:  & 
        sec_per_min, & ! Variable(s)
        fstderr ! Constant(s)

    use clubb_precision, only:  & 
        time_precision

#ifdef NETCDF
    use input_netcdf, only: open_netcdf_read, close_netcdf_read ! Procedure(s)
#endif

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O unit

    character(len=*), intent(in) ::filename ! Path to the file and its name

    logical, intent(in) :: l_restart     ! Whether this is a restart run

    real( kind = time_precision ), intent(in) ::  & 
      time ! Time near which we want to find GrADS output,
    ! e.g. time_restart     [s]


    ! Output Variable(s)
    integer, intent(out) ::  & 
      nearest_timestep ! Nearest integer output time to time [units vary]

    ! Local Variables
    type (stat_file) :: fread_var

    real( kind = core_rknd ) :: delta_time   ! In seconds

    logical :: l_grads_file, l_error

    ! ---- Begin Code ----

    l_grads_file = .not. l_netcdf_file( filename )

    if ( l_grads_file ) then
      ! Read in the control file
      call open_grads_read( iunit, trim( filename ), & ! In
                            fread_var,  & ! In/Out
                            l_error ) ! Out
    
    else
#ifdef NETCDF
      if( stats_input_type == sam_input_type ) then
        call open_netcdf_read( 'U', trim( filename ), & ! In
                               fread_var, & ! In/Out
                               l_error ) ! Out
      else
        call open_netcdf_read( 'thlm', trim( filename ), & ! In
                               fread_var, & ! In/Out
                               l_error ) ! Out
      end if
#else
      write(fstderr,*) "This version of CLUBB was not compiled with netCDF support"
#endif
    end if

    if ( l_error ) then
      write(fstderr,*) "Error reading file " // trim( filename )
      stop "Fatal error"
    end if

    if (l_grads_file) then
      ! (restart time) - (initial time)
      delta_time = real(time - (fread_var%time - &
                   real(fread_var%dtwrite, kind=time_precision)), kind=core_rknd)
    else
      !  Eric   Raut     June  2013
      delta_time = real(time - fread_var%time, kind=core_rknd)
    !    Joshua Fasching March 2008
!     .        time - fread_var%time
    end if

    ! Reporting
    if ( l_restart ) then
      print *, "Initial time of GrADS reference file ", & 
               "[seconds since midnight]: ",  & 
               fread_var%time
      print *, "Model restart time [s]: ", time
      print *, "Elapsed time between ", & 
               "initial time of ref file and restart time [s]: ",  & 
               delta_time
      print *, "GrADS file output time interval [s]: ",  & 
               fread_var%dtwrite

      if ( ( mod( delta_time , fread_var%dtwrite )  > 1e-8_core_rknd) .or.  & 
           ( mod( delta_time, fread_var%dtwrite ) < -1e-8_core_rknd) ) then
        write(fstderr,*) "Error: Elapsed time is not a multiple ", & 
                "of the reference GrADS output time interval."
        write(fstderr,*) "Elapsed time [s] = ", delta_time
        write(fstderr,*) "GrADS output time interval = ", fread_var%dtwrite
        stop
      end if

      if ( mod( delta_time , sec_per_min ) > 1e-8_core_rknd& 
            .or. mod( delta_time, sec_per_min ) < -1e-8_core_rknd) then
        write(fstderr,*) "Error: Elapsed time is not a multiple ", & 
                "of one minute."
        write(fstderr,*) "Elapsed time [s] = ", delta_time
        stop
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
      print *, "Nearest GrADS output time iteration [ ]: ", & 
               nint( real( nearest_timestep, kind=core_rknd) /  & 
                     fread_var%dtwrite/sec_per_min ) - 1
    end if ! l_restart

    if ( l_grads_file ) then
      call close_grads_read( fread_var )
    else
#ifdef NETCDF
      call close_netcdf_read( fread_var )
#endif
    end if

  end subroutine compute_timestep
!===============================================================================
  subroutine get_clubb_variable_interpolated &
           ( l_input_var, filename, varname, vardim, timestep, &
             clubb_heights, variable_interpolated, l_read_error )
! Description:
!   Obtain a profile of a CLUBB variable from a GrADS file and interpolate if
!   needed.
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
      vardim, & ! Dimension of the profile in CLUBB
      timestep  ! Timestep to read in

    real( kind = core_rknd ), intent(in), dimension(vardim) :: &
      clubb_heights ! Altitudes         [m]

    real( kind = core_rknd ), intent(inout), dimension(vardim) :: &
      variable_interpolated     ! The variable after interpolation

    logical, intent(out) :: &
      l_read_error ! Whether there was a read error

    logical :: l_spec_bound_cond
    ! Local Variables

    ! ---- Begin Code ----

    if ( l_input_var ) then
      if ( clubb_heights(1) < 0.0_core_rknd ) then
        l_spec_bound_cond = .true.
      else
        l_spec_bound_cond = .false.
      end if

      variable_interpolated = stat_file_average( filename, vardim, timestep, &
        timestep, clubb_heights, varname, npower, l_spec_bound_cond, l_read_error )

    else
      l_read_error = .false.

    end if

    return
  end subroutine get_clubb_variable_interpolated

  !--------------------------------------------------------------------------
  subroutine get_input_variables_interp( nvars, variables , timestep, &
                                         k_lowest_zt, k_highest_zt,   &
                                         k_lowest_zm, k_highest_zm, l_error )

  ! Description:
  !   Obtains profiles for COAMPS or SAM variables from a GrADS or NetCDF file and
  !   interpolates if needed.
  ! References:
  !   None
  !--------------------------------------------------------------------------------

    use stat_file_module, only: stat_file

    use grid_class, only: gr

    use input_grads, only: &
      open_grads_read, &
      get_grads_var, &
      close_grads_read

#ifdef NETCDF
    use input_netcdf, only: &
      open_netcdf_read, &
      get_netcdf_var, &
      close_netcdf_read
#endif

    use stat_file_utils, only: & 
       LES_grid_to_CLUBB_grid, & ! Procedure(s)
       CLUBB_levels_within_LES_domain, &
       l_netcdf_file

    use constants_clubb, only: &
      fstderr  ! Variable(s)

    implicit none

    integer, intent(in) :: &
      nvars, & ! The size of the variables array
      timestep ! The timestep to read the variables at.

    type (input_field), dimension(nvars), intent(in) :: &
      variables ! The list of external variables to be read in.
  
    ! These variables are output by this subroutine in case certain variables
    ! need special adjustments.
    integer, dimension(1:size(stat_files)), intent(out) ::  &
      k_lowest_zt, &  ! The lowest CLUBB thermodynamic level that's within the LES domain.
      k_highest_zt, & ! The highest CLUBB thermodynamic level that's within the LES domain.
      k_lowest_zm, &  ! The lowest CLUBB momentum level that's within the LES domain.
      k_highest_zm    ! The highest CLUBB momentum level that's within the LES domain.

    logical, intent(out) :: &
      l_error ! Error flag

    ! Local variables

    real( kind = core_rknd), dimension(:,:), allocatable :: &
      LES_tmp ! Temporary variable

    integer :: &
      i, &          ! index variable
      file_index, & ! index of the files
      min_ia, &     ! minimum ia value from all files
      max_iz        ! maximum iz value from all files

    logical :: l_grads_file, & ! use grads or netcdf
      l_internal_error

    type (input_field) :: &
      current_var ! The current variable

    type (stat_file), dimension(1:size(stat_files)) :: &
      fread_vars

    integer, parameter :: &
      base_unit_number = 15

    logical :: &
      l_convert_to_MKS ! convert inputs to MKS units

    !--------------------------- BEGIN CODE ---------------------------------

    l_error = .false.

    l_grads_file = .not. l_netcdf_file(stat_files(1))

    if( l_grads_file ) then
      do file_index=1, size(stat_files)
        call open_grads_read( base_unit_number + (file_index-1), stat_files(file_index), &
                              fread_vars(file_index), l_internal_error )
        l_error = l_error .or. l_internal_error
      end do ! file_index=1, size(stat_files)
    else

#ifdef NETCDF
      do file_index=1, size(stat_files)
        call open_netcdf_read( "U", stat_files(file_index),  &
                              fread_vars(file_index), l_internal_error )
        l_error = l_error .or. l_internal_error
      end do ! file_index=1, size(stat_files)
#else
      write(fstderr,*) "This version of CLUBB was not compiled with netCDF support"
      stop "Fatal error"

#endif
    end if

    if (l_error) then
      write(fstderr,*) "A fatal error occured while reading the input file."
      stop "Fatal error in input_fields.get_input_variables_interp"
    end if

    ! initialize to boundary values to determine min and max later
    min_ia = 9999
    max_iz = -1

    do file_index=1, size(stat_files)
      ! Find the lowest and highest indices of CLUBB thermodynamic levels that
      ! fall within the domain of the LES output.
      call CLUBB_levels_within_LES_domain( fread_vars(file_index), gr%zt,  &
                      k_lowest_zt(file_index), k_highest_zt(file_index) )

      ! Find the lowest and highest indices of CLUBB momentum levels that
      ! fall within the domain of the LES output.
      call CLUBB_levels_within_LES_domain( fread_vars(file_index), gr%zm,  &
                      k_lowest_zm(file_index), k_highest_zm(file_index) )

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

        if( l_grads_file ) then
          call get_grads_var( fread_vars(file_index), current_var%input_name, timestep, &
                      LES_tmp(file_index, fread_vars(file_index)%ia:fread_vars(file_index)%iz), &
                      l_internal_error )
        else
#ifdef NETCDF
          l_convert_to_MKS = .false.
          call get_netcdf_var( fread_vars(file_index), current_var%input_name, timestep, &
                      l_convert_to_MKS,&
                      LES_tmp(file_index, fread_vars(file_index)%ia:fread_vars(file_index)%iz), &
                      l_internal_error )
#else
          write(fstderr,*) "This version of CLUBB was not compiled with netCDF support"
          l_internal_error = .true.
#endif
        end if ! l_grads_file

        l_error = l_error .or. l_internal_error

        if( current_var%clubb_grid_type == "zt" ) then
          call inputfields_interp_and_adjust( current_var%clubb_grid_type, & ! In
                    fread_vars(file_index), current_var%adjustment, & ! In
                    LES_tmp(file_index, :), k_lowest_zt(file_index), k_highest_zt(file_index), &!In
                      current_var%clubb_var ) ! In/out
        else if( current_var%clubb_grid_type == "zm" ) then
          call inputfields_interp_and_adjust( current_var%clubb_grid_type, & ! In
                    fread_vars(file_index), current_var%adjustment, & ! In
                    LES_tmp(file_index, :), k_lowest_zm(file_index), k_highest_zm(file_index), &!In
                      current_var%clubb_var ) ! In/out
        end if ! current_var%clubb_grid_type == "zt"

      end if ! current_var%l_input_var
    end do ! i=1, nvars

    if( l_grads_file ) then
      do file_index=1, size(stat_files)
        call close_grads_read( fread_vars(file_index) )
      end do ! file_index=1, size(stat_files)
    else

#ifdef NETCDF
      do file_index=1, size(stat_files)
        call close_netcdf_read( fread_vars(file_index) )
      end do ! file_index=1, size(stat_files)
#else
      write(fstderr,*) "This version of CLUBB was not compiled with netCDF support"
      stop "Fatal error"

#endif
    end if ! l_grads_file

    deallocate(LES_tmp)

  end subroutine get_input_variables_interp
        
  !----------------------------------------------------------------------------------------
  subroutine inputfields_interp_and_adjust( clubb_grid_type, fread_var, adjustment, &
                                            LES_tmp, k_lowest_input, k_highest_input, &
                                            clubb_var )

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

    use grid_class, only: &
      gr ! Type

    implicit none

    ! Input Variable(s)

    character(len = 3), intent(in) :: &
      clubb_grid_type ! The type of the grid, either "zm" or "zt"

    type (stat_file), intent(in) :: &
      fread_var

    real( kind = core_rknd ), intent(in) :: &
      adjustment ! The value used to adjust the variable for unit conversions

    real( kind = core_rknd ), dimension(fread_var%ia:fread_var%iz), intent(in) :: &
      LES_tmp ! The values that were read in from the file.

    integer, intent(in) :: &
      k_lowest_input, &  ! The lowest CLUBB level that's within the LES domain.
      k_highest_input ! The highest CLUBB level that's within the LES domain.

    ! Output Variable(s)
    
    real( kind = core_rknd ), dimension(:), pointer :: &
      clubb_var ! The clubb variable to store the result in

    ! Local Variable(s)

    ! Variables used to reconcile CLUBB levels with LES vertical levels.
    integer, dimension(k_lowest_input:k_highest_input) ::  &
      exact_lev_idx, & ! In case of an exact match, index of LES level that is
                          ! exactly even with CLUBB thermodynamic level k.
      lower_lev_idx, & ! In case linear interpolation is needed, index of LES
                          ! level that is immediately below CLUBB thermo. level k.
      upper_lev_idx    ! In case linear interpolation is needed, index of LES
                          ! level that is immediately above CLUBB thermo. level k.

    logical, dimension(k_lowest_input:k_highest_input) ::  &
      l_lin_int  ! Flag that is turned on if linear interpolation is needed.

    integer :: &
      k ! loop counter
    
    real( kind = core_rknd ), pointer, dimension(:) :: &
      grid ! the momentum or thermo grid

    !----------------BEGIN CODE----------------

    if ( clubb_grid_type == "zt" ) then
      grid => gr%zt
    else if ( clubb_grid_type == "zm" ) then
      grid => gr%zm
    end if

    ! For all CLUBB levels, k, that are within the LES domain,
    ! find either the index of the LES level that exactly matches the altitude
    ! of the CLUBB level, or find the two indices of the LES levels that are on
    ! either side of the CLUBB level.
    do k = k_lowest_input, k_highest_input, 1

      ! CLUBB vertical level k is found at an altitude that is within the
      ! domain of the LES output.
      call LES_grid_to_CLUBB_grid( fread_var, grid, k,  &
                                   exact_lev_idx(k), lower_lev_idx(k),  &
                                   upper_lev_idx(k), l_lin_int(k) )
    enddo

    !LES_tmp is the value of the variable from the LES GrADS file.
    do k = k_lowest_input, k_highest_input, 1
      if( l_lin_int(k) ) then
        ! CLUBB level k is found at an altitude that is between two
        ! LES levels.  Linear interpolation is required.
         clubb_var(k) = lin_interpolate_two_points( grid(k), &
                          fread_var%z(upper_lev_idx(k)), &
                          fread_var%z(lower_lev_idx(k)), &
                          LES_tmp(upper_lev_idx(k)), &
                          LES_tmp(lower_lev_idx(k)) )

      else
        ! CLUBB level k is found at an altitude that is an exact
        ! match with an LES level altitude.
        clubb_var(k) = LES_tmp(exact_lev_idx(k))
      end if
      if( adjustment /= 1.0_core_rknd ) then
        clubb_var(k) = clubb_var(k) &
                         * adjustment
      end if
    end do

    return
  end subroutine inputfields_interp_and_adjust

!-------------------------------------------------------------------------------
  subroutine inputfields_init( iunit, namelist_filename )

!-------------------------------------------------------------------------------
    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      namelist_filename

    integer, intent(in) :: iunit

    ! Local variables
    character(len=120) :: datafile

    ! Namelist definitions
    namelist /setfields/ datafile, input_type, & 
      l_input_um, l_input_vm, l_input_rtm, l_input_thlm, & 
      l_input_wp2, l_input_wprtp, l_input_wpthlp,  & 
      l_input_wp3, l_input_rtp2, l_input_thlp2,  & 
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
      l_input_precip_frac, l_input_rtprrp, l_input_thlprrp, & 
      l_input_rtpNrp, l_input_thlpNrp, l_input_rrpNrp, & 
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
    l_input_thlp2 = .false.
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

    ! Setup the GrADS file reader
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
