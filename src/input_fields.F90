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
  character(len=100), public :: & 
    stat_file_zt, stat_file_zm, stat_file_sfc

  character(len=10), public :: input_type

  logical, public :: &
    l_input_um, l_input_vm, l_input_rtm, l_input_thlm, & 
    l_input_wp2, l_input_wprtp, l_input_wpthlp,  & 
    l_input_wp3, l_input_rtp2, l_input_thlp2,  & 
    l_input_rtpthlp, l_input_upwp, l_input_vpwp, & 
    l_input_ug, l_input_vg, l_input_rcm,  & 
    l_input_wm_zt, l_input_exner, l_input_em, & 
    l_input_p, l_input_rho, l_input_rho_zm, & 
    l_input_rho_ds_zm, l_input_rho_ds_zt, &
    l_input_thv_ds_zm, l_input_thv_ds_zt, &
    l_input_Lscale, l_input_Lscale_up, l_input_Lscale_down, & 
    l_input_Kh_zt, l_input_Kh_zm, l_input_tau_zm, l_input_tau_zt, & 
    l_input_wpthvp, l_input_radht, &
    l_input_thl1, l_input_thl2, l_input_mixt_frac, l_input_s1, l_input_s2, &
    l_input_stdev_s1, l_input_stdev_s2, l_input_rc1, l_input_rc2, &
    l_input_thvm, l_input_rrainm, l_input_Nrm,  l_input_Ncm,  & 
    l_input_rsnowm, l_input_ricem, l_input_rgraupelm, l_input_Ncnm, l_input_Nim, & 
    l_input_thlm_forcing, l_input_rtm_forcing, & 
    l_input_up2, l_input_vp2, l_input_sigma_sqd_w, & 
    l_input_cloud_frac, l_input_sigma_sqd_w_zt, &
    l_input_veg_T_in_K, l_input_deep_soil_T_in_K, &
    l_input_sfc_soil_T_in_K

  integer, parameter, private :: &
    coamps_input_type = 1, &
    clubb_input_type =  2, &
    sam_input_type = 3

  integer, parameter, private :: &
    num_sam_inputfields = 62, & ! The number of input fields for SAM
    num_coamps = 61 ! The number of input fields for coamps

  integer, private :: &
    stats_input_type

  ! Procedures
  public  :: stat_fields_reader, &
             compute_timestep, &
             inputfields_init, &
             set_filenames

  private ! Default Scope

  type input_field
 
    logical :: l_input_var ! If .true., this variable should be read

    character(len = 50) :: &
      input_name, & ! The name of the variable in SAM 
      clubb_name  ! The name of the variable in CLUBB

    character(len = 3) :: &
      grid_type ! The type of the grid, either "zm" or "zt"

    real( kind = core_rknd ), dimension(:), pointer :: &
      clubb_var ! The clubb variable to store the result in

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
      stat_file_zt = trim( file_prefix )//"_coamps_sm.ctl"
      stat_file_zm = trim( file_prefix )//"_coamps_sw.ctl"
      stat_file_sfc = trim( file_prefix )//"_coamps_sfc.ctl"

      stats_input_type = coamps_input_type

    case ( "clubb_scm", "CLUBB_SCM", "clubb" )
      stat_file_zt = trim( file_prefix )//"_zt.ctl"
      stat_file_zm = trim( file_prefix )//"_zm.ctl"
      stat_file_sfc = trim( file_prefix )//"_sfc.ctl"

      stats_input_type = clubb_input_type

      inquire(file=stat_file_zt,exist=l_grads)

      if ( .not. l_grads ) then
        write(fstdout,*) "inputfields: Cannot find GrADS ctl file, assuming netCDF input"
        stat_file_zt = trim( file_prefix )// "_zt.nc"
        stat_file_zm = trim( file_prefix )// "_zm.nc"
        stat_file_sfc = trim( file_prefix )// "_sfc.nc"
      end if

    case ( "SAM", "sam" )
      ! SAM uses one netCDF file for all of its output, so the entire filename
      ! should be provided in file_prefix
      stat_file_zt = trim( file_prefix )
      stat_file_zm = trim( file_prefix )
      stat_file_sfc = trim( file_prefix )

      stats_input_type = sam_input_type

    case default
      write(fstderr,*) "Don't know how to handle input_type = "// & 
        input_type
      stop

    end select

    return
  end subroutine set_filenames
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine stat_fields_reader( timestep )
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
        sigma_sqd_w

    use variables_diagnostic_module, only: & 
        hydromet,  & ! Variable(s)
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
        Ncnm, & 
        sigma_sqd_w_zt, & 
        em, &
        radht

    use variables_prognostic_module, only: & 
        pdf_params ! Variable(s)

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
        g_per_kg

    use array_index, only:  & 
        iirrainm, iiNrm, iirsnowm, iiricem, iirgraupelm, iiNim, iiNcm

    use stat_file_utils, only: & 
      LES_grid_to_CLUBB_grid, & ! Procedure(s)
      CLUBB_levels_within_LES_domain

    use extrapolation, only: &
      lin_ext_zt_bottom, &
      lin_ext_zm_bottom

    use parameters_microphys, only: &
      micro_scheme ! Variable(s)

    use soil_vegetation, only: deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: max, trim, any

    ! Arguments
    integer, intent(in) :: timestep

    ! Local Variables
    logical :: l_read_error, l_fatal_error

    real( kind = core_rknd ), dimension(gr%nz+1) :: tmp1

    integer ::  &
      k_lowest_zt_input, &  ! The lowest CLUBB thermodynamic level that's within the LES domain.
      k_highest_zt_input, & ! The highest CLUBB thermodynamic level that's within the LES domain.
      k_lowest_zm_input, &  ! The lowest CLUBB momentum level that's within the LES domain.
      k_highest_zm_input    ! The highest CLUBB momentum level that's within the LES domain.

    integer :: k  ! Array index
 
    real( kind = core_rknd), dimension(gr%nz), target :: &
      temp_Nrm, temp_Ncm, temp_rgraupelm, temp_ricem, &
      temp_rrainm, temp_rsnowm, temp_tke, temp_wpup_sgs, temp_wpvp_sgs, &
      temp_Nim ! temporary variable

    type (input_field), dimension(:), allocatable :: &
      SAM_variables, & ! A list of SAM variables to read in.
      coamps_variables ! A list of coamps variables to read in.


    ! ---- Begin Code ----

    select case ( stats_input_type )

    case ( clubb_input_type )

      !  Thermo grid - zt file

      ! Initialize l_fatal_error for case clubb_input_type
      l_fatal_error = .false.

      call get_clubb_variable_interpolated &
           ( l_input_um, stat_file_zt, "um", gr%nz, timestep, &
             gr%zt, um, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vm, stat_file_zt, "vm", gr%nz, timestep, &
             gr%zt, vm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtm, stat_file_zt, "rtm", gr%nz, timestep, &
             gr%zt, rtm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_thlm, stat_file_zt, "thlm", gr%nz, timestep, &
             gr%zt, thlm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wp3, stat_file_zt, "wp3", gr%nz, timestep, &
             gr%zt, wp3, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_tau_zt, stat_file_zt, "tau_zt", gr%nz, timestep, &
             gr%zt, tau_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rrainm, stat_file_zt, "rrainm", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_rrainm ) then
        hydromet(1:gr%nz,iirrainm) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rsnowm, stat_file_zt, "rsnowm", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_rsnowm ) then
        hydromet(1:gr%nz,iirsnowm) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_ricem, stat_file_zt, "ricem", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_ricem ) then
        hydromet(1:gr%nz,iiricem) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rgraupelm, stat_file_zt, "rgraupelm", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_rgraupelm ) then
        hydromet(1:gr%nz,iirgraupelm) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

!--------------------------------------------------------
! Added variables for clubb_restart
      call get_clubb_variable_interpolated &
           ( l_input_p, stat_file_zt, "p_in_Pa", gr%nz, timestep, &
             gr%zt, p_in_Pa, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_exner, stat_file_zt, "exner", gr%nz, timestep, &
             gr%zt, exner, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_ug, stat_file_zt, "ug", gr%nz, timestep, &
             gr%zt, ug, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vg, stat_file_zt, "vg", gr%nz, timestep, &
             gr%zt, vg, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rcm, stat_file_zt, "rcm", gr%nz, timestep, &
             gr%zt, rcm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wm_zt, stat_file_zt, "wm", gr%nz, timestep, &
             gr%zt, wm_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( l_input_rho, stat_file_zt, "rho", gr%nz, timestep, &
             gr%zt, rho, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rho_ds_zt, stat_file_zt, "rho_ds_zt", gr%nz, timestep, &
             gr%zt, rho_ds_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thv_ds_zt, stat_file_zt, "thv_ds_zt", gr%nz, timestep, &
             gr%zt, thv_ds_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Lscale, stat_file_zt, "Lscale", gr%nz, timestep, &
             gr%zt, Lscale, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Lscale_up, stat_file_zt, "Lscale_up", gr%nz, timestep, &
             gr%zt, Lscale_up, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Lscale_down, stat_file_zt, "Lscale_down", gr%nz, timestep, &
             gr%zt, Lscale_down, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Kh_zt, stat_file_zt, "Kh_zt", gr%nz, timestep, &
             gr%zt, Kh_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thvm, stat_file_zt, "thvm", gr%nz, timestep, &
             gr%zt, thvm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thlm_forcing, stat_file_zt, "thlm_forcing", gr%nz, timestep, &
             gr%zt, thlm_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtm_forcing, stat_file_zt, "rtm_forcing", gr%nz, timestep, &
             gr%zt, rtm_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Ncm, stat_file_zt, "Ncm", gr%nz, timestep, &
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_Ncm ) then
        hydromet(1:gr%nz,iiNcm) = tmp1(1:gr%nz)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Ncnm, stat_file_zt, "Ncnm", gr%nz, timestep, &
             gr%zt, Ncnm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Nim, stat_file_zt, "Nim", gr%nz, timestep, &
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_Nim ) then
        hydromet(1:gr%nz, iiNcm) = tmp1(1:gr%nz)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_cloud_frac, stat_file_zt, "cloud_frac", gr%nz, timestep, &
             gr%zt, cloud_frac, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Nrm, stat_file_zt, "Nrm", gr%nz, timestep, &
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( l_input_Nrm ) then
        hydromet(1:gr%nz, iiNrm) = tmp1(1:gr%nz)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_sigma_sqd_w_zt, stat_file_zt, "sigma_sqd_w_zt", gr%nz, timestep, &
             gr%zt, sigma_sqd_w_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_radht, stat_file_zt, "radht", gr%nz, timestep, &
             gr%zt, radht, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      ! PDF Parameters (needed for K&K microphysics)
      call get_clubb_variable_interpolated &
           ( l_input_thl1, stat_file_zt, "thl1", gr%nz, timestep, &
             gr%zt, pdf_params%thl1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thl2, stat_file_zt, "thl2", gr%nz, timestep, &
             gr%zt, pdf_params%thl2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_mixt_frac, stat_file_zt, "mixt_frac", gr%nz, timestep, &
             gr%zt, pdf_params%mixt_frac, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_s1, stat_file_zt, "s1", gr%nz, timestep, &
             gr%zt, pdf_params%s1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_s2, stat_file_zt, "s2", gr%nz, timestep, &
             gr%zt, pdf_params%s2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_stdev_s1, stat_file_zt, "stdev_s1", gr%nz, timestep, &
             gr%zt, pdf_params%stdev_s1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_stdev_s2, stat_file_zt, "stdev_s2", gr%nz, timestep, &
             gr%zt, pdf_params%stdev_s2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rc1, stat_file_zt, "rc1", gr%nz, timestep, &
             gr%zt, pdf_params%rc1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rc2, stat_file_zt, "rc2", gr%nz, timestep, &
             gr%zt, pdf_params%rc2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

!--------------------------------------------------------

      ! Read in the zm file

      call get_clubb_variable_interpolated &
           ( l_input_wp2, stat_file_zm, "wp2", gr%nz, timestep, &
             gr%zm, wp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wprtp, stat_file_zm, "wprtp", gr%nz, timestep, &
             gr%zm, wprtp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wpthlp, stat_file_zm, "wpthlp", gr%nz, timestep, &
             gr%zm, wpthlp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_wpthvp, stat_file_zm, "wpthvp", gr%nz, timestep, &
             gr%zm, wpthvp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtp2, stat_file_zm, "rtp2", gr%nz, timestep, &
             gr%zm, rtp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thlp2, stat_file_zm, "thlp2", gr%nz, timestep, &
             gr%zm, thlp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rtpthlp, stat_file_zm, "rtpthlp", gr%nz, timestep, &
             gr%zm, rtpthlp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_upwp, stat_file_zm, "upwp", gr%nz, timestep, &
             gr%zm, upwp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vpwp, stat_file_zm, "vpwp", gr%nz, timestep, &
             gr%zm, vpwp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

!-----------------------------------------------------------
      call get_clubb_variable_interpolated &
           ( l_input_em, stat_file_zm, "em", gr%nz, timestep, &
             gr%zm, em, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rho_zm, stat_file_zm, "rho_zm", gr%nz, timestep, &
             gr%zm, rho_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_rho_ds_zm, stat_file_zm, "rho_ds_zm", gr%nz, timestep, &
             gr%zm, rho_ds_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_thv_ds_zm, stat_file_zm, "thv_ds_zm", gr%nz, timestep, &
             gr%zm, thv_ds_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_Kh_zm, stat_file_zm, "Kh_zm", gr%nz, timestep, &
             gr%zm, Kh_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_tau_zm, stat_file_zm, "tau_zm", gr%nz, timestep, &
             gr%zm, tau_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_up2, stat_file_zm, "up2", gr%nz, timestep, &
             gr%zm, up2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_vp2, stat_file_zm, "vp2", gr%nz, timestep, &
             gr%zm, vp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( l_input_sigma_sqd_w, stat_file_zm, "sigma_sqd_w", gr%nz, timestep, &
             gr%zm, sigma_sqd_w, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


!-----------------------------------------------------------

      call get_clubb_variable_interpolated &
           ( l_input_veg_T_in_K, stat_file_sfc, "veg_T_in_K", 1, timestep, &
             (/0._core_rknd/), tmp1(1), l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error
      veg_T_in_K = tmp1(1)

      call get_clubb_variable_interpolated &
           ( l_input_deep_soil_T_in_K, stat_file_sfc, "deep_soil_T_in_K", 1, timestep, &
             (/0._core_rknd/), tmp1(1), l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error
      deep_soil_T_in_K = tmp1(1)

      call get_clubb_variable_interpolated &
           ( l_input_sfc_soil_T_in_K, stat_file_sfc, "sfc_soil_T_in_K", 1, timestep, &
             (/0._core_rknd/), tmp1(1), l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error
      sfc_soil_T_in_K = tmp1(1)



      if ( l_fatal_error ) stop "oops, get_grads_var failed in stat_fields_reader"

    case ( coamps_input_type ) ! COAMPS LES stats data     
      l_fatal_error = .false.

      allocate( coamps_variables( num_coamps ) )

      k = 1

      coamps_variables(k)%l_input_var = l_input_um
      coamps_variables(k)%input_name = "um"
      coamps_variables(k)%clubb_var => um
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zt"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_vm
      coamps_variables(k)%input_name = "vm"
      coamps_variables(k)%clubb_var => vm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zt"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rtm
      coamps_variables(k)%input_name = "qtm"
      coamps_variables(k)%clubb_var => rtm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zt"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_thlm
      coamps_variables(k)%input_name = "thlm"
      coamps_variables(k)%clubb_var => thlm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zt"

      k = k + 1

      ! We obtain wp2 from stats_sw

      coamps_variables(k)%l_input_var = l_input_wp3
      coamps_variables(k)%input_name = "wp3"
      coamps_variables(k)%clubb_var => wp3
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zt"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_wprtp
      coamps_variables(k)%input_name = "wpqtp"
      coamps_variables(k)%clubb_var => wprtp
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_wpthlp
      coamps_variables(k)%input_name = "wpthlp"
      coamps_variables(k)%clubb_var => wpthlp
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rtp2
      coamps_variables(k)%input_name = "qtp2"
      coamps_variables(k)%clubb_var => rtp2
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_thlp2
      coamps_variables(k)%input_name = "thlp2"
      coamps_variables(k)%clubb_var => thlp2
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rtpthlp
      coamps_variables(k)%input_name = "qtpthlp"
      coamps_variables(k)%clubb_var => rtpthlp
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

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
      coamps_variables(k)%grid_type = "zt"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_wm_zt
      coamps_variables(k)%input_name = "wlsm"
      coamps_variables(k)%clubb_var => wm_zt
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zt"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_exner
      coamps_variables(k)%input_name = "ex0"
      coamps_variables(k)%clubb_var => exner
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zt"

      k = k + 1

      ! Note that em (clubb) = em (coamps) + tke (coamps), so just store tke
      ! in a temporary variable for now and add them together later.
      coamps_variables(k)%l_input_var = l_input_em
      coamps_variables(k)%input_name = "em"
      coamps_variables(k)%clubb_var => em
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      temp_tke = 0.0_core_rknd! Initialize temp to 0.0

      coamps_variables(k)%l_input_var = l_input_em
      coamps_variables(k)%input_name = "tke"
      coamps_variables(k)%clubb_var => temp_tke
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_p
      coamps_variables(k)%input_name = "pm"
      coamps_variables(k)%clubb_var => p_in_Pa
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zt"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rho
      coamps_variables(k)%input_name = "dn0"
      coamps_variables(k)%clubb_var => rho
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zt"

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
      coamps_variables(k)%grid_type = "zt"

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
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_thl1
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "thl1"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_thl2
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "thl2"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_mixt_frac
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "mixt_frac"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_s1
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "s1"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_s2
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "s2"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_stdev_s1
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "stdev_s1"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_stdev_s2
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "stdev_s2"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rc1
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "rc1"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rc2
      coamps_variables(k)%input_name = "none"
      coamps_variables(k)%clubb_name = "rc2"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_thvm
      coamps_variables(k)%input_name = "thvm"
      coamps_variables(k)%clubb_var => thvm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zt"

      k = k + 1

      temp_rrainm = 0.0_core_rknd! initialize to 0.0

      if ( l_input_rrainm ) then
        if ( iirrainm < 1 ) then
            write(fstderr,*) "Rain water mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_rrainm
          coamps_variables(k)%input_name = "qrm"
          coamps_variables(k)%clubb_var => temp_rrainm
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%grid_type = "zt"

          k = k + 1
        end if
      end if

      temp_Nrm = 0.0_core_rknd! Initialize to 0.0

      if ( l_input_Nrm ) then
        if ( iiNrm < 1 ) then
            write(fstderr,*) "Rain droplet number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_Nrm
          coamps_variables(k)%input_name = "nrm"
          coamps_variables(k)%clubb_var => temp_Nrm
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%grid_type = "zt"

          k = k + 1
        end if
      end if

      temp_Ncm = 0.0_core_rknd! initialize to 0.0

      if ( l_input_Ncm ) then
        if ( iiNcm < 1 ) then
            write(fstderr,*) "Cloud droplet number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_Ncm
          coamps_variables(k)%input_name = "ncm"
          coamps_variables(k)%clubb_var => temp_Ncm
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%grid_type = "zt"

          k = k + 1
        end if
      end if

      temp_rsnowm = 0.0_core_rknd! initialize to 0.0

      if ( l_input_rsnowm ) then
        if ( iirsnowm < 1 ) then
            write(fstderr,*) "Snow mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_rsnowm
          coamps_variables(k)%input_name = "qsm"
          coamps_variables(k)%clubb_var => temp_rsnowm
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%grid_type = "zt"

          k = k + 1
        end if
      end if

      temp_ricem = 0.0_core_rknd! initialize to 0.0

      if ( l_input_ricem ) then
        if ( iiricem < 1 ) then
            write(fstderr,*) "Ice mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_ricem
          coamps_variables(k)%input_name = "qim"
          coamps_variables(k)%clubb_var => temp_ricem
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%grid_type = "zt"

          k = k + 1
        end if
      end if

      temp_rgraupelm = 0.0_core_rknd! initialize to 0.0

      if ( l_input_rgraupelm ) then
        if ( iirgraupelm < 1 ) then
            write(fstderr,*) "Graupel mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_rgraupelm
          coamps_variables(k)%input_name = "qgm"
          coamps_variables(k)%clubb_var => temp_rgraupelm
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%grid_type = "zt"

          k = k + 1
        end if
      end if

      coamps_variables(k)%l_input_var = l_input_Ncnm
      coamps_variables(k)%input_name = "ncnm"
      coamps_variables(k)%clubb_var => Ncnm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zt"

      k = k + 1

      temp_Nim = 0.0_core_rknd! initialize to 0.0

      if ( l_input_Nim ) then
        if ( iiNim < 1 ) then
            write(fstderr,*) "Ice number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          coamps_variables(k)%l_input_var = l_input_Nim
          coamps_variables(k)%input_name = "nim"
          coamps_variables(k)%clubb_var => temp_Nim
          coamps_variables(k)%adjustment = 1.0_core_rknd
          coamps_variables(k)%grid_type = "zt"

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
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_vp2
      coamps_variables(k)%input_name = "vp2"
      coamps_variables(k)%clubb_var => vp2
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

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

      ! Note:  wpup_sgs and wpvp_sgs must be added to make the u'w' and v'w' terms
      !        as they are in CLUBB.
      coamps_variables(k)%l_input_var = l_input_upwp
      coamps_variables(k)%input_name = "wpup"
      coamps_variables(k)%clubb_var => upwp
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      temp_wpup_sgs = 0.0_core_rknd! initialize to 0.0

      coamps_variables(k)%l_input_var = l_input_upwp
      coamps_variables(k)%input_name = "wpup_sgs"
      coamps_variables(k)%clubb_var => temp_wpup_sgs
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_vpwp
      coamps_variables(k)%input_name = "wpvp"
      coamps_variables(k)%clubb_var => vpwp
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      temp_wpvp_sgs = 0.0_core_rknd! initialize to 0.0

      coamps_variables(k)%l_input_var = l_input_vpwp
      coamps_variables(k)%input_name = "wpvp_sgs"
      coamps_variables(k)%clubb_var => temp_wpvp_sgs
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_wp2
      coamps_variables(k)%input_name = "wp2"
      coamps_variables(k)%clubb_var => wp2
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      k = k + 1

      coamps_variables(k)%l_input_var = l_input_rho_zm
      coamps_variables(k)%input_name = "dn0"
      coamps_variables(k)%clubb_var => rho_zm
      coamps_variables(k)%adjustment = 1.0_core_rknd
      coamps_variables(k)%grid_type = "zm"

      ! Initialize l_read_error for case ( "les" )
      l_read_error = .false.

      if ( l_fatal_error ) stop "oops, bad input for inputfields"
      
      call get_input_variables_interp &
             (k, coamps_variables, stat_file_zt, timestep, &
              k_lowest_zt_input, k_highest_zt_input, &
              k_lowest_zm_input, k_highest_zm_input, l_fatal_error )

      if ( l_input_um ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, use the
        ! values of um at thermodynamic levels 3 and 2 to find the value at
        ! thermodynamic level 1 through the use of a linear extension.
        if ( k_lowest_zt_input == 2 ) then
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
        if ( k_lowest_zt_input == 2 ) then
          vm(1) = lin_ext_zt_bottom( vm(3), vm(2), & 
                                     gr%zt(3), gr%zt(2), gr%zt(1) )
        endif
      endif

      if ( l_input_rtm ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of rtm at thermodynamic level 1 to the value at thermodynamic
        ! level 2, as it is done in advance_xm_wpxp.
        if ( k_lowest_zt_input == 2 ) then
          rtm(1) = rtm(2)
        endif
      endif

      if ( l_input_thlm ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of thlm at thermodynamic level 1 to the value at thermodynamic
        ! level 2, as it is done in advance_xm_wpxp.
        if ( k_lowest_zt_input == 2 ) then
          thlm(1) = thlm(2)
        endif
      endif

      if ( l_input_wp3 ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of wp3 at thermodynamic level 1 to 0, as it is done in
        ! advance_wp2_wp3.
        if ( k_lowest_zt_input == 2 ) then
          wp3(1) = 0.0_core_rknd
        endif
      endif

      if ( l_input_rrainm ) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iirrainm) = &
                    temp_rrainm(k_lowest_zt_input:k_highest_zt_input)
      end if

      if ( l_input_Nrm ) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iiNrm) = &
                    temp_Nrm(k_lowest_zt_input:k_highest_zt_input)
      end if

      if ( l_input_Ncm ) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iiNcm) = &
                    temp_Ncm(k_lowest_zt_input:k_highest_zt_input)
      end if

      if ( l_input_rsnowm ) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iirsnowm) = &
                    temp_rsnowm(k_lowest_zt_input:k_highest_zt_input)
      end if

      if ( l_input_ricem ) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iiricem) = &
                    temp_ricem(k_lowest_zt_input:k_highest_zt_input)
      end if

      if ( l_input_rgraupelm ) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iirgraupelm) = &
                    temp_rgraupelm(k_lowest_zt_input:k_highest_zt_input)
      end if

      if ( l_input_Nim ) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iiNim) = &
                    temp_Nim(k_lowest_zt_input:k_highest_zt_input)
      end if

      if ( l_input_wprtp ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wprtp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like surface_varnce.
        if ( k_lowest_zm_input == 2 ) then
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
        if ( k_lowest_zm_input == 2 ) then
          wpthlp(1) = lin_ext_zm_bottom( wpthlp(3), wpthlp(2), & 
                                         gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if ( l_input_rtp2 ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, set the value of rtp2
        ! at momentum level 1 to the value at momentum level 2.
        ! Using a linear extension here resulted in negatives.
        if ( k_lowest_zm_input == 2 ) then
          rtp2(1) =  rtp2(2)
        endif
        if ( any( rtp2(1:gr%nz) < rt_tol**2 ) ) then
          do k=1, gr%nz
            rtp2(k) = max(rtp2(k), rt_tol**2)
          enddo
        endif
      endif

      if ( l_input_thlp2 ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, set the value of thlp2
        ! at momentum level 1 to the value at momentum level 2.
        ! Using a linear extension here resulted in negatives.
        if ( k_lowest_zm_input == 2 ) then
          thlp2(1) = thlp2(2)
        endif
        if ( any( thlp2(1:gr%nz) < thl_tol**2 ) ) then
          do k=1, gr%nz
            thlp2(k) = max(thlp2(k), thl_tol**2)
          enddo
        endif
      endif

      if ( l_input_rtpthlp) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! rtpthlp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like surface_varnce.
        if ( k_lowest_zm_input == 2 ) then
          rtpthlp(1) = lin_ext_zm_bottom( rtpthlp(3), rtpthlp(2), & 
                                          gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if( l_input_em ) then
        do k=k_lowest_zm_input, k_highest_zm_input, 1
          em(k) = em(k) + temp_tke(k)
        end do

        where ( em < em_min ) em = em_min
      end if

      if ( l_input_up2 ) then
        ! Clip up2 to be no smaller than w_tol_sqd
        where ( up2 < w_tol_sqd ) up2 = w_tol_sqd
      endif

      if ( l_input_vp2 ) then
         where ( vp2 < w_tol_sqd ) vp2 = w_tol_sqd
      endif

       if ( l_input_upwp) then
        do k=k_lowest_zm_input, k_highest_zm_input, 1
          upwp(k) = upwp(k) + temp_wpup_sgs(k)
        end do
      endif

      if ( l_input_vpwp) then
        do k=k_lowest_zm_input, k_highest_zm_input, 1
          vpwp(k) = temp_wpvp_sgs(k)
        end do
      endif

      if ( l_input_wp2 ) then
         if ( any( wp2(1:gr%nz) < w_tol_sqd ) ) then
          do k=1, gr%nz
            wp2(k) = max( wp2(k), w_tol_sqd )
          enddo
        endif
      endif

      deallocate(coamps_variables)

    case( sam_input_type )

      allocate( SAM_variables(num_sam_inputfields) )

      l_fatal_error = .false.

      k = 1

      SAM_variables(k)%l_input_var = l_input_um
      SAM_variables(k)%input_name = "U"
      SAM_variables(k)%clubb_var => um
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_vm
      SAM_variables(k)%input_name = "V"
      SAM_variables(k)%clubb_var => vm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rtm
      SAM_variables(k)%input_name = "QT"
      SAM_variables(k)%clubb_var => rtm
      SAM_variables(k)%adjustment = 1.0e-3_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thlm
      SAM_variables(k)%input_name = "THETAL"
      SAM_variables(k)%clubb_var => thlm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_wp2
      SAM_variables(k)%input_name = "W2"
      SAM_variables(k)%clubb_var => wp2
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zm"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rho
      SAM_variables(k)%input_name = "RHO"
      SAM_variables(k)%clubb_var => rho
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      ! Note that this needs to be adjusted by 1/(RHO * Lv)
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_wprtp
      SAM_variables(k)%input_name = "QTFLUX"
      SAM_variables(k)%clubb_var => wprtp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zm"

      k = k + 1

      ! Note that this needs to be adjusted by 1/(RHO * Cp)
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_wpthlp
      SAM_variables(k)%input_name = "TLFLUX"
      SAM_variables(k)%clubb_var => wpthlp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zm"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_wp3
      SAM_variables(k)%input_name = "W3"
      SAM_variables(k)%clubb_var => wp3
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rtp2
      SAM_variables(k)%input_name = "QT2"
      SAM_variables(k)%clubb_var => rtp2
      SAM_variables(k)%adjustment = 1.0e-6_core_rknd
      SAM_variables(k)%grid_type = "zm"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thlp2
      SAM_variables(k)%input_name = "TL2"
      SAM_variables(k)%clubb_var => thlp2
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zm"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rtpthlp
      SAM_variables(k)%clubb_name = "rtpthlp"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_upwp
      SAM_variables(k)%input_name = "UW"
      SAM_variables(k)%clubb_var => upwp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zm"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_vpwp
      SAM_variables(k)%input_name = "VW"
      SAM_variables(k)%clubb_var => vpwp
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zm"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_ug
      SAM_variables(k)%input_name = "UOBS"
      SAM_variables(k)%clubb_var => ug
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_vg
      SAM_variables(k)%input_name = "VOBS"
      SAM_variables(k)%clubb_var => vg
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rcm
      SAM_variables(k)%input_name = "QCL"
      SAM_variables(k)%clubb_var => rcm
      SAM_variables(k)%adjustment = 1.0e-3_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_wm_zt
      SAM_variables(k)%input_name = "WOBS"
      SAM_variables(k)%clubb_var => wm_zt
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_exner
      SAM_variables(k)%clubb_name = "exner"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_em
      SAM_variables(k)%input_name = "TKE"
      SAM_variables(k)%clubb_var => em
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zm"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_p
      SAM_variables(k)%input_name = "PRES"
      SAM_variables(k)%clubb_var => p_in_Pa
      SAM_variables(k)%adjustment = pascal_per_mb
      SAM_variables(k)%grid_type = "zt"

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
      SAM_variables(k)%grid_type = "zt"

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

      SAM_variables(k)%l_input_var = l_input_thl1
      SAM_variables(k)%clubb_name = "thl1"
      SAM_variables(k)%input_name = "none"

      k = k + 1
      SAM_variables(k)%l_input_var = l_input_thl2
      SAM_variables(k)%clubb_name = "thl2"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_mixt_frac
      SAM_variables(k)%clubb_name = "mixt_frac"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_s1
      SAM_variables(k)%clubb_name = "s1"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_s2
      SAM_variables(k)%clubb_name = "s2"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_stdev_s1
      SAM_variables(k)%clubb_name = "stdev_s1"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_stdev_s2
      SAM_variables(k)%clubb_name = "stdev_s2"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rc1
      SAM_variables(k)%clubb_name = "rc1"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_rc2
      SAM_variables(k)%clubb_name = "rc2"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_thvm
      SAM_variables(k)%input_name = "THETAV"
      SAM_variables(k)%clubb_var => thvm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      temp_rrainm = 0.0_core_rknd! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_rrainm
      SAM_variables(k)%input_name = "QPL"
      SAM_variables(k)%clubb_var => temp_rrainm
      SAM_variables(k)%adjustment = 1.0_core_rknd/g_per_kg
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      temp_Nrm = 0.0_core_rknd! initialize to 0.0

      ! Note that this needs to be adjusted by 1e6/RHO
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_Nrm
      SAM_variables(k)%input_name = "NR"
      SAM_variables(k)%clubb_var => temp_Nrm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1
      temp_Ncm = 0.0_core_rknd! initialize to 0.0

      ! Note that this needs to be adjusted by 1e6/RHO
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_Ncm
      SAM_variables(k)%input_name = "NC"
      SAM_variables(k)%clubb_var => temp_Ncm
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      temp_rsnowm = 0.0_core_rknd! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_rsnowm
      SAM_variables(k)%input_name = "QS"
      SAM_variables(k)%clubb_var => temp_rsnowm
      SAM_variables(k)%adjustment = 1.0_core_rknd/g_per_kg
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      temp_ricem = 0.0_core_rknd! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_ricem
      SAM_variables(k)%input_name = "QI"
      SAM_variables(k)%clubb_var => temp_ricem
      SAM_variables(k)%adjustment = 1.0_core_rknd/g_per_kg
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      temp_rgraupelm = 0.0_core_rknd ! initialize to 0.0

      SAM_variables(k)%l_input_var = l_input_rgraupelm
      SAM_variables(k)%input_name = "QG"
      SAM_variables(k)%clubb_var => temp_rgraupelm
      SAM_variables(k)%adjustment = 1.0_core_rknd/g_per_kg
      SAM_variables(k)%grid_type = "zt"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_Ncnm
      SAM_variables(k)%clubb_name = "Ncnm"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      temp_Nim = 0.0_core_rknd ! initialize to 0.0

      ! Note that this needs to be adjusted by 1e6/RHO
      ! This will need to be adjusted outside of get_sam_variable_interpolated
      SAM_variables(k)%l_input_var = l_input_Nim
      SAM_variables(k)%input_name = "NI"
      SAM_variables(k)%clubb_var => temp_Nim
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

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
      SAM_variables(k)%grid_type = "zm"

      k = k + 1

      ! Note that V2 is on the thermodynamic grid in SAM while vp2 is
      ! on the momentum grid in CLUBB. This will have to be interpolated.
      SAM_variables(k)%l_input_var = l_input_vp2
      SAM_variables(k)%input_name = "V2"
      SAM_variables(k)%clubb_var => vp2
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zm"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_sigma_sqd_w
      SAM_variables(k)%clubb_name = "sigma_sqd_w"
      SAM_variables(k)%input_name = "none"

      k = k + 1

      SAM_variables(k)%l_input_var = l_input_cloud_frac
      SAM_variables(k)%input_name = "CLD"
      SAM_variables(k)%clubb_var => cloud_frac
      SAM_variables(k)%adjustment = 1.0_core_rknd
      SAM_variables(k)%grid_type = "zt"

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

      call get_input_variables_interp &
             (k, SAM_variables, stat_file_zt, timestep, &
              k_lowest_zt_input, k_highest_zt_input, &
              k_lowest_zm_input, k_highest_zm_input, l_fatal_error )

      if( l_fatal_error ) then
        STOP "Error reading SAM input fields"
      end if

      ! Perform non-constant adjustments
      do k=k_lowest_zt_input, k_highest_zt_input
        if(l_input_Nrm) then
          ! Nrm = NR * (1.0e6 / RHO)
          temp_Nrm(k) = temp_Nrm(k) * (1.0e6_core_rknd / rho(k))
        end if

        if(l_input_Ncm) then
          ! Ncm = NC * (1.0e6 / RHO)
          temp_Ncm(k) = temp_Ncm(k) * (1.0e6_core_rknd / rho(k))
        end if

        if(l_input_Nim) then
          ! Nim = NI * (1.0e6 / RHO)
          temp_Nim(k) = temp_Nim(k) * (1.0e7_core_rknd / rho(k))
        end if
 
      end do

      do k=k_lowest_zm_input, k_highest_zm_input
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
       
      ! Add hydromet variables.
      if( l_input_Nrm .and. iiNrm > 0) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iiNrm) = &
                   temp_Nrm(k_lowest_zt_input:k_highest_zt_input)
      end if

      if( l_input_Ncm .and. iiNcm > 0) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iiNcm) = &
                   temp_Ncm(k_lowest_zt_input:k_highest_zt_input)
      end if

      if( l_input_Nim .and. iiNim > 0) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iiNim) = &
                 temp_Nim(k_lowest_zt_input:k_highest_zt_input)
      end if

      if( l_input_rrainm .and. iirrainm > 0) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iirrainm) = &
                   temp_rrainm(k_lowest_zt_input:k_highest_zt_input)
      end if

      if( l_input_rgraupelm .and. iirgraupelm > 0) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iirgraupelm) = &
                   temp_rgraupelm(k_lowest_zt_input:k_highest_zt_input)
      end if

      if( l_input_ricem .and. iiricem > 0) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iiricem) = &
                   temp_ricem(k_lowest_zt_input:k_highest_zt_input)
      end if

      if( l_input_rsnowm .and. iirsnowm > 0) then
        hydromet(k_lowest_zt_input:k_highest_zt_input,iirsnowm) = &
                   temp_rsnowm(k_lowest_zt_input:k_highest_zt_input)
      end if

      if ( l_fatal_error ) stop "Failed to read inputfields for SAM."

      deallocate( SAM_variables )

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
    ! Given a time 'time', determines the closest output time in a GrADS file.

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

    real(kind=time_precision), intent(in) ::  & 
      time ! Time near which we want to find GrADS output,
    ! e.g. time_restart     [s]


    ! Output Variable(s)
    integer, intent(out) ::  & 
      nearest_timestep ! Nearest GrADS output time to time [min]

    ! Local Variables
    type (stat_file) :: fread_var

    real(kind=time_precision) :: delta_time   ! In seconds

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

    ! (restart time) - (initial time)
    delta_time = time - (fread_var%time - fread_var%dtwrite)

    !    Joshua Fasching March 2008
!     .        time - fread_var%time

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

      if ( ( mod( delta_time , fread_var%dtwrite )  > 1e-8_time_precision ) .or.  & 
           ( mod( delta_time, fread_var%dtwrite ) < -1e-8_time_precision ) ) then
        write(fstderr,*) "Error: Elapsed time is not a multiple ", & 
                "of the reference GrADS output time interval."
        write(fstderr,*) "Elapsed time [s] = ", delta_time
        write(fstderr,*) "GrADS output time interval = ", fread_var%dtwrite
        stop
      end if

      if ( mod( delta_time , sec_per_min ) > 1e-8_time_precision & 
            .or. mod( delta_time, sec_per_min ) < -1e-8_time_precision ) then
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
               nint( real( nearest_timestep, kind=time_precision ) /  & 
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
  subroutine get_input_variables_interp( nvars, variables, filename, timestep, &
                                   k_lowest_zt_input, k_highest_zt_input, &
                                   k_lowest_zm_input, k_highest_zm_input, l_error )

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

    character(len=*), intent(in) :: &
      filename ! The name of the file that the variables are located in.
  
    ! These variables are output by this subroutine in case certain variables
    ! need special adjustments.
    integer, intent(out) ::  &
      k_lowest_zt_input, &  ! The lowest CLUBB thermodynamic level that's within the LES domain.
      k_highest_zt_input, & ! The highest CLUBB thermodynamic level that's within the LES domain.
      k_lowest_zm_input, &  ! The lowest CLUBB momentum level that's within the LES domain.
      k_highest_zm_input    ! The highest CLUBB momentum level that's within the LES domain.

    logical, intent(out) :: &
      l_error ! Error flag

    ! Local variables

    real( kind = core_rknd), dimension(:), allocatable :: &
      LES_tmp ! Temporary variable

    integer :: &
      i ! index variable


    logical l_grads_file, & ! use grads or netcdf
            l_internal_error, &
            l_convert_to_MKS ! convert inputs to MKS units

    type (input_field) :: &
      current_var ! The current variable

    type (stat_file) :: &
      fread_var

    integer, parameter :: &
      unit_number = 15


    !--------------------------- BEGIN CODE ---------------------------------

    l_error = .false.

    l_grads_file = .not. l_netcdf_file(filename)

    if( l_grads_file ) then
      call open_grads_read( unit_number, filename, &
                           fread_var, l_internal_error )

    else

#ifdef NETCDF
      call open_netcdf_read( "THETAL", filename,  &
                            fread_var, l_internal_error )
#else
      write(fstderr,*) "This version of CLUBB was not compiled with netCDF support"
      write(fstderr,*) "Error reading file "// trim( filename )
      stop "Fatal error"

#endif
    end if

    l_error = l_error .or. l_internal_error

    ! Find the lowest and highest indices of CLUBB thermodynamic levels that
    ! fall within the domain of the LES output.
    call CLUBB_levels_within_LES_domain( fread_var, gr%zt,  &
                                       k_lowest_zt_input, k_highest_zt_input )

    ! Find the lowest and highest indices of CLUBB momentum levels that
    ! fall within the domain of the LES output.
    call CLUBB_levels_within_LES_domain( fread_var, gr%zm,  &
                                       k_lowest_zm_input, k_highest_zm_input )

    allocate( LES_tmp(fread_var%ia:fread_var%iz) )

    do i=1, nvars

      current_var = variables(i);

      if( current_var%input_name == "none" .and. &
          current_var%l_input_var ) then
        write(fstderr,*) "The variable " //trim( current_var%clubb_name ) &
           //" is not setup for input_type " //trim( input_type )
        l_error = .true.
        cycle
      end if

      if( current_var%l_input_var ) then

        if( l_grads_file ) then
          call get_grads_var( fread_var, current_var%input_name, timestep, &
                        LES_tmp(fread_var%ia:fread_var%iz), l_error )
        else
#ifdef NETCDF
          l_convert_to_MKS = .false.
          call get_netcdf_var( fread_var, current_var%input_name, timestep, l_convert_to_MKS, &
                        LES_tmp(fread_var%ia:fread_var%iz), l_error )
#else
          write(fstderr,*) "This version of CLUBB was not compiled with netCDF support"
          l_error = .true.
#endif
        end if
 
        l_error = l_error .or. l_internal_error

        call inputfields_interp_and_adjust( current_var%clubb_var, current_var%grid_type, &
                                            fread_var, current_var%adjustment, &
                                            LES_tmp, k_lowest_zt_input, k_highest_zt_input, &
                                            k_lowest_zm_input, k_highest_zm_input )

      end if
    end do

    if( l_grads_file ) then
      call close_grads_read( fread_var )

    else

#ifdef NETCDF
      call close_netcdf_read( fread_var )
#else
      write(fstderr,*) "This version of CLUBB was not compiled with netCDF support"
      write(fstderr,*) "Error reading file "// trim( filename )
      stop "Fatal error"

#endif
    end if

  end subroutine get_input_variables_interp
        
  !----------------------------------------------------------------------------------------
  subroutine inputfields_interp_and_adjust( clubb_var, grid_type, fread_var, adjustment, &
                                       LES_tmp, k_lowest_zt_input, k_highest_zt_input,& 
                                       k_lowest_zm_input, k_highest_zm_input )
  !
  !  Interpolates an input field as needed and adjusts it for unit conversion.
  !
  !  References:
  !    none
  !----------------------------------------------------------------------------------------
    
    use interpolation, only: &
      lin_int ! Procedure(s)

    use stat_file_utils, only: & 
      LES_grid_to_CLUBB_grid, & ! Procedure(s)
      CLUBB_levels_within_LES_domain

    use stat_file_module, only: &
      stat_file ! Type

    use grid_class, only: &
      gr ! Type

    use constants_clubb, only: &
      fstderr  ! Variable(s)

    implicit none

    real( kind = core_rknd ), dimension(:), pointer, intent(in) :: &
      clubb_var ! The clubb variable to store the result in

    character(len = 3), intent(in) :: &
      grid_type ! The type of the grid, either "zm" or "zt"

    type (stat_file), intent(in) :: &
      fread_var

    real( kind = core_rknd ), intent(in) :: &
      adjustment ! The value used to adjust the variable for unit conversions

    real( kind = core_rknd ), dimension(fread_var%ia:fread_var%iz), intent(in) :: &
      LES_tmp ! The values that were read in from the file.

    integer, intent(in) :: &
      k_lowest_zt_input, &  ! The lowest CLUBB thermodynamic level that's within the LES domain.
      k_highest_zt_input, & ! The highest CLUBB thermodynamic level that's within the LES domain.
      k_lowest_zm_input, &  ! The lowest CLUBB momentum level that's within the LES domain.
      k_highest_zm_input    ! The highest CLUBB momentum level that's within the LES domain.

    ! Variables used to reconcile CLUBB thermodynamic levels with LES vertical levels.
    integer, dimension(k_lowest_zt_input:k_highest_zt_input) ::  &
      exact_lev_idx_zt, & ! In case of an exact match, index of LES level that is
                          ! exactly even with CLUBB thermodynamic level k.
      lower_lev_idx_zt, & ! In case linear interpolation is needed, index of LES
                          ! level that is immediately below CLUBB thermo. level k.
      upper_lev_idx_zt    ! In case linear interpolation is needed, index of LES
                          ! level that is immediately above CLUBB thermo. level k.

    logical, dimension(k_lowest_zt_input:k_highest_zt_input) ::  &
      l_lin_int_zt  ! Flag that is turned on if linear interpolation is needed.

    ! Variables used to reconcile CLUBB momentum levels with LES vertical levels.
    integer, dimension(k_lowest_zm_input:k_highest_zm_input) ::  &
      exact_lev_idx_zm, & ! In case of an exact match, index of LES level that is
                          ! exactly even with CLUBB momentum level k.
      lower_lev_idx_zm, & ! In case linear interpolation is needed, index of LES
                          ! level that is immediately below CLUBB momentum level k
      upper_lev_idx_zm    ! In case linear interpolation is needed, index of LES
                          ! level that is immediately above CLUBB momentum level k

    logical, dimension(k_lowest_zm_input:k_highest_zm_input) ::  &
      l_lin_int_zm  ! Flag that is turned on if linear interpolation is needed.

    integer :: &
      k ! loop counter

    ! For all CLUBB thermodynamic levels, k, that are within the LES domain,
    ! find either the index of the LES level that exactly matches the altitude
    ! of the CLUBB level, or find the two indices of the LES levels that are on
    ! either side of the CLUBB level.
    do k = k_lowest_zt_input, k_highest_zt_input, 1

      ! CLUBB vertical level k is found at an altitude that is within the
      ! domain of the LES output.
      call LES_grid_to_CLUBB_grid( fread_var, gr%zt, k,  &
                                 exact_lev_idx_zt(k), lower_lev_idx_zt(k),  &
                                 upper_lev_idx_zt(k), l_lin_int_zt(k) )
    enddo

    ! For all CLUBB momentum levels, k, that are within the LES domain,
    ! find either the index of the LES level that exactly matches the altitude
    ! of the CLUBB level, or find the two indices of the LES levels that are on
    ! either side of the CLUBB level.
    do k = k_lowest_zm_input, k_highest_zm_input, 1
      ! CLUBB vertical level k is found at an altitude that is within the
      ! domain of the LES output.
      call LES_grid_to_CLUBB_grid( fread_var, gr%zm, k,  &
                                 exact_lev_idx_zm(k), lower_lev_idx_zm(k),  &
                                 upper_lev_idx_zm(k), l_lin_int_zm(k) )
    enddo

    !LES_tmp is the value of the variable from the LES GrADS file.
    if( grid_type == "zt" ) then
      do k = k_lowest_zt_input, k_highest_zt_input, 1
        if( l_lin_int_zt(k) ) then
          ! CLUBB level k is found at an altitude that is between two
          ! LES levels.  Linear interpolation is required.
           clubb_var(k) = lin_int( gr%zt(k), &
                            fread_var%z(upper_lev_idx_zt(k)), &
                            fread_var%z(lower_lev_idx_zt(k)), &
                            LES_tmp(upper_lev_idx_zt(k)), &
                            LES_tmp(lower_lev_idx_zt(k)) )
        else
          ! CLUBB level k is found at an altitude that is an exact
          ! match with an LES level altitude.
          clubb_var(k) = LES_tmp(exact_lev_idx_zt(k))
        end if
        if( adjustment /= 1.0_core_rknd ) then
          clubb_var(k) = clubb_var(k) &
                           * adjustment
        end if
      end do
    else if( grid_type == "zm" ) then
      do k = k_lowest_zm_input, k_highest_zm_input, 1
        if( l_lin_int_zm(k) ) then
          ! CLUBB level k is found at an altitude that is between two
          ! LES levels.  Linear interpolation is required.
          clubb_var(k) = lin_int( gr%zm(k), &
                           fread_var%z(upper_lev_idx_zm(k)), &
                           fread_var%z(lower_lev_idx_zm(k)), &
                           LES_tmp(upper_lev_idx_zm(k)), &
                           LES_tmp(lower_lev_idx_zm(k)) )
        else
          ! CLUBB level k is found at an altitude that is an exact
          ! match with an LES level altitude.
          clubb_var(k) = LES_tmp(exact_lev_idx_zm(k))
        end if
        if( adjustment /= 1.0_core_rknd ) then
          clubb_var(k) = clubb_var(k) &
                           * adjustment
        end if
      end do
    else
      write(fstderr,*) "Invalid grid_type specified in input field: "// grid_type
      STOP "Grid type does not exist."
    end if

  end subroutine inputfields_interp_and_adjust

!-------------------------------------------------------------------------------
  subroutine inputfields_init( iunit, namelist_filename )

! Description:
!   This subroutine reads in a namelist to determine which variables are to be 
!   read in, etc.

! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      namelist_filename

    integer, intent(in) :: iunit

    ! Local variables
    character(len=98) :: datafile

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
      l_input_thl1, l_input_thl2, l_input_mixt_frac, l_input_s1, l_input_s2, &
      l_input_stdev_s1, l_input_stdev_s2, l_input_rc1, l_input_rc2, &
      l_input_thvm, l_input_rrainm,l_input_Nrm,  & 
      l_input_rsnowm, l_input_ricem, l_input_rgraupelm,  & 
      l_input_thlm_forcing, l_input_rtm_forcing, & 
      l_input_up2, l_input_vp2, l_input_sigma_sqd_w, l_input_Ncm,  & 
      l_input_Ncnm, l_input_Nim, l_input_cloud_frac, l_input_sigma_sqd_w_zt, &
      l_input_veg_T_in_K, l_input_deep_soil_T_in_K, &
      l_input_sfc_soil_T_in_K

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
    l_input_thl1 = .false.
    l_input_thl2 = .false.
    l_input_mixt_frac = .false.
    l_input_s1 = .false.
    l_input_s2 = .false.
    l_input_stdev_s1 = .false.
    l_input_stdev_s2 = .false.
    l_input_rc1 = .false.
    l_input_rc2 = .false.
    l_input_thvm = .false.
    l_input_rrainm = .false.
    l_input_Nrm = .false.
    l_input_rsnowm = .false.
    l_input_ricem = .false.
    l_input_rgraupelm = .false.
    l_input_thlm_forcing = .false.
    l_input_rtm_forcing = .false.
    l_input_up2 = .false.
    l_input_vp2 = .false.
    l_input_sigma_sqd_w = .false.
    l_input_Ncm = .false.
    l_input_Ncnm = .false.
    l_input_Nim = .false.
    l_input_cloud_frac = .false.
    l_input_sigma_sqd_w_zt = .false.
    l_input_veg_T_in_K = .false.
    l_input_deep_soil_T_in_K = .false.
    l_input_sfc_soil_T_in_K  = .false.

    print *, namelist_filename
    ! Read in our namelist
    open(unit=iunit, file=namelist_filename, status='old', action='read')

    read(unit=iunit, nml=setfields)

    close(unit=iunit)

    ! Setup the GrADS file reader
    call set_filenames( datafile )    
      

    return
  end subroutine inputfields_init
!-----------------------------------------------------------------------

!===============================================================================

end module inputfields
