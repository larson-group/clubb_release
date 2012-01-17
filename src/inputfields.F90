!-----------------------------------------------------------------------
! $Id$

module inputfields
!  Description:
!    This file contains information and subroutines needed for restart
!    simulations and simulations where the variables are reset at every
!    timestep, called a "clubb_inputfields" simulation.
!===============================================================================

  implicit none

  ! Run information
  character(len=100), public :: & 
    stat_file_zt, stat_file_zm, stat_file_sfc

  character(len=10), public :: input_type

  logical, public :: &
    input_um, input_vm, input_rtm, input_thlm, & 
    input_wp2, input_wprtp, input_wpthlp,  & 
    input_wp3, input_rtp2, input_thlp2,  & 
    input_rtpthlp, input_upwp, input_vpwp, & 
    input_ug, input_vg, input_rcm,  & 
    input_wm_zt, input_exner, input_em, & 
    input_p, input_rho, input_rho_zm, & 
    input_rho_ds_zm, input_rho_ds_zt, &
    input_thv_ds_zm, input_thv_ds_zt, &
    input_Lscale, input_Lscale_up, input_Lscale_down, & 
    input_Kh_zt, input_Kh_zm, input_tau_zm, input_tau_zt, & 
    input_wpthvp, input_radht, &
    input_thl1, input_thl2, input_mixt_frac, input_s1, input_s2, &
    input_stdev_s1, input_stdev_s2, input_rc1, input_rc2, &
    input_thvm, input_rrainm, input_Nrm,  input_Ncm,  & 
    input_rsnowm, input_ricem, input_rgraupelm, input_Ncnm, input_Nim, & 
    input_thlm_forcing, input_rtm_forcing, & 
    input_up2, input_vp2, input_sigma_sqd_w, & 
    input_cloud_frac, input_sigma_sqd_w_zt, &
    input_veg_T_in_K, input_deep_soil_T_in_K, &
    input_sfc_soil_T_in_K

  integer, parameter, private :: &
    coamps_input_type = 1, &
    clubb_input_type =  2

  integer, private :: &
    stats_input_type

  ! Procedures
  public  :: stat_fields_reader, &
             compute_timestep, &
             inputfields_init, &
             set_filenames

  private ! Default Scope

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
        fstderr

    use array_index, only:  & 
        iirrainm, iiNrm, iirsnowm, iiricem, iirgraupelm, iiNim, iiNcm

    use stat_file_module, only: & 
        stat_file     ! Type

    use stat_file_utils, only: & 
       LES_grid_to_CLUBB_grid, & ! Procedure(s)
       CLUBB_levels_within_LES_domain

    use input_grads, only: & 
        open_grads_read, & ! Procedure(s)
        close_grads_read

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

    type (stat_file) :: fread_var

    real( kind = core_rknd ), dimension(:), allocatable :: LES_tmp1

    real( kind = core_rknd ), dimension(gr%nz+1) :: tmp1

    integer ::  &
      k_lowest_zt_input, &  ! The lowest CLUBB thermodynamic level that's within the LES domain.
      k_highest_zt_input, & ! The highest CLUBB thermodynamic level that's within the LES domain.
      k_lowest_zm_input, &  ! The lowest CLUBB momentum level that's within the LES domain.
      k_highest_zm_input    ! The highest CLUBB momentum level that's within the LES domain.

    ! Variables used to reconcile CLUBB thermodynamic levels with LES vertical levels.
    integer, dimension(:), allocatable ::  &
      exact_lev_idx_zt, & ! In case of an exact match, index of LES level that is
                        ! exactly even with CLUBB thermodynamic level k.
      lower_lev_idx_zt, & ! In case linear interpolation is needed, index of LES
                        ! level that is immediately below CLUBB thermo. level k.
      upper_lev_idx_zt  ! In case linear interpolation is needed, index of LES
                        ! level that is immediately above CLUBB thermo. level k.

    logical, dimension(:), allocatable ::  &
      l_lin_int_zt  ! Flag that is turned on if linear interpolation is needed.

    ! Variables used to reconcile CLUBB momentum levels with LES vertical levels.
    integer, dimension(:), allocatable ::  &
      exact_lev_idx_zm, & ! In case of an exact match, index of LES level that is
                        ! exactly even with CLUBB momentum level k.
      lower_lev_idx_zm, & ! In case linear interpolation is needed, index of LES
                        ! level that is immediately below CLUBB momentum level k
      upper_lev_idx_zm    ! In case linear interpolation is needed, index of LES
    ! level that is immediately above CLUBB momentum level k

    logical, dimension(:), allocatable ::  &
      l_lin_int_zm  ! Flag that is turned on if linear interpolation is needed.

    integer :: k, &  ! Array index
               unit_number ! file unit number

    real( kind = core_rknd ), dimension(gr%nz) :: &
      temp ! temporary variable

    ! ---- Begin Code ----

    select case ( stats_input_type )

    case ( clubb_input_type )

      !  Thermo grid - zt file

      ! Initialize l_fatal_error for case clubb_input_type
      l_fatal_error = .false.

      call get_clubb_variable_interpolated &
           ( input_um, stat_file_zt, "um", gr%nz, timestep, &
             gr%zt, um, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_vm, stat_file_zt, "vm", gr%nz, timestep, &
             gr%zt, vm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rtm, stat_file_zt, "rtm", gr%nz, timestep, &
             gr%zt, rtm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( input_thlm, stat_file_zt, "thlm", gr%nz, timestep, &
             gr%zt, thlm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_wp3, stat_file_zt, "wp3", gr%nz, timestep, &
             gr%zt, wp3, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_tau_zt, stat_file_zt, "tau_zt", gr%nz, timestep, &
             gr%zt, tau_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rrainm, stat_file_zt, "rrainm", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( input_rrainm ) then
        hydromet(1:gr%nz,iirrainm) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rsnowm, stat_file_zt, "rsnowm", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( input_rsnowm ) then
        hydromet(1:gr%nz,iirsnowm) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_ricem, stat_file_zt, "ricem", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( input_ricem ) then
        hydromet(1:gr%nz,iiricem) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rgraupelm, stat_file_zt, "rgraupelm", gr%nz, timestep, & 
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( input_rgraupelm ) then
        hydromet(1:gr%nz,iirgraupelm) = tmp1(1:gr%nz)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

!--------------------------------------------------------
! Added variables for clubb_restart
      call get_clubb_variable_interpolated &
           ( input_p, stat_file_zt, "p_in_Pa", gr%nz, timestep, &
             gr%zt, p_in_Pa, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_exner, stat_file_zt, "exner", gr%nz, timestep, &
             gr%zt, exner, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_ug, stat_file_zt, "ug", gr%nz, timestep, &
             gr%zt, ug, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_vg, stat_file_zt, "vg", gr%nz, timestep, &
             gr%zt, vg, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rcm, stat_file_zt, "rcm", gr%nz, timestep, &
             gr%zt, rcm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_wm_zt, stat_file_zt, "wm", gr%nz, timestep, &
             gr%zt, wm_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( input_rho, stat_file_zt, "rho", gr%nz, timestep, &
             gr%zt, rho, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rho_ds_zt, stat_file_zt, "rho_ds_zt", gr%nz, timestep, &
             gr%zt, rho_ds_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_thv_ds_zt, stat_file_zt, "thv_ds_zt", gr%nz, timestep, &
             gr%zt, thv_ds_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Lscale, stat_file_zt, "Lscale", gr%nz, timestep, &
             gr%zt, Lscale, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Lscale_up, stat_file_zt, "Lscale_up", gr%nz, timestep, &
             gr%zt, Lscale_up, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Lscale_down, stat_file_zt, "Lscale_down", gr%nz, timestep, &
             gr%zt, Lscale_down, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Kh_zt, stat_file_zt, "Kh_zt", gr%nz, timestep, &
             gr%zt, Kh_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_thvm, stat_file_zt, "thvm", gr%nz, timestep, &
             gr%zt, thvm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_thlm_forcing, stat_file_zt, "thlm_forcing", gr%nz, timestep, &
             gr%zt, thlm_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rtm_forcing, stat_file_zt, "rtm_forcing", gr%nz, timestep, &
             gr%zt, rtm_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Ncm, stat_file_zt, "Ncm", gr%nz, timestep, &
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( input_Ncm ) then
        hydromet(1:gr%nz,iiNcm) = tmp1(1:gr%nz)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Ncnm, stat_file_zt, "Ncnm", gr%nz, timestep, &
             gr%zt, Ncnm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Nim, stat_file_zt, "Nim", gr%nz, timestep, &
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( input_Nim ) then
        hydromet(1:gr%nz, iiNcm) = tmp1(1:gr%nz)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_cloud_frac, stat_file_zt, "cloud_frac", gr%nz, timestep, &
             gr%zt, cloud_frac, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Nrm, stat_file_zt, "Nrm", gr%nz, timestep, &
             gr%zt, tmp1(1:gr%nz), l_read_error )
      if ( input_Nrm ) then
        hydromet(1:gr%nz, iiNrm) = tmp1(1:gr%nz)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_sigma_sqd_w_zt, stat_file_zt, "sigma_sqd_w_zt", gr%nz, timestep, &
             gr%zt, sigma_sqd_w_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_radht, stat_file_zt, "radht", gr%nz, timestep, &
             gr%zt, radht, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      ! PDF Parameters (needed for K&K microphysics)
      call get_clubb_variable_interpolated &
           ( input_thl1, stat_file_zt, "thl1", gr%nz, timestep, &
             gr%zt, pdf_params%thl1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_thl2, stat_file_zt, "thl2", gr%nz, timestep, &
             gr%zt, pdf_params%thl2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_mixt_frac, stat_file_zt, "mixt_frac", gr%nz, timestep, &
             gr%zt, pdf_params%mixt_frac, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_s1, stat_file_zt, "s1", gr%nz, timestep, &
             gr%zt, pdf_params%s1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_s2, stat_file_zt, "s2", gr%nz, timestep, &
             gr%zt, pdf_params%s2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_stdev_s1, stat_file_zt, "stdev_s1", gr%nz, timestep, &
             gr%zt, pdf_params%stdev_s1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_stdev_s2, stat_file_zt, "stdev_s2", gr%nz, timestep, &
             gr%zt, pdf_params%stdev_s2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rc1, stat_file_zt, "rc1", gr%nz, timestep, &
             gr%zt, pdf_params%rc1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rc2, stat_file_zt, "rc2", gr%nz, timestep, &
             gr%zt, pdf_params%rc2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

!--------------------------------------------------------

      ! Read in the zm file

      call get_clubb_variable_interpolated &
           ( input_wp2, stat_file_zm, "wp2", gr%nz, timestep, &
             gr%zm, wp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_wprtp, stat_file_zm, "wprtp", gr%nz, timestep, &
             gr%zm, wprtp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_wpthlp, stat_file_zm, "wpthlp", gr%nz, timestep, &
             gr%zm, wpthlp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_wpthvp, stat_file_zm, "wpthvp", gr%nz, timestep, &
             gr%zm, wpthvp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rtp2, stat_file_zm, "rtp2", gr%nz, timestep, &
             gr%zm, rtp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_thlp2, stat_file_zm, "thlp2", gr%nz, timestep, &
             gr%zm, thlp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rtpthlp, stat_file_zm, "rtpthlp", gr%nz, timestep, &
             gr%zm, rtpthlp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_upwp, stat_file_zm, "upwp", gr%nz, timestep, &
             gr%zm, upwp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_vpwp, stat_file_zm, "vpwp", gr%nz, timestep, &
             gr%zm, vpwp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

!-----------------------------------------------------------
      call get_clubb_variable_interpolated &
           ( input_em, stat_file_zm, "em", gr%nz, timestep, &
             gr%zm, em, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rho_zm, stat_file_zm, "rho_zm", gr%nz, timestep, &
             gr%zm, rho_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rho_ds_zm, stat_file_zm, "rho_ds_zm", gr%nz, timestep, &
             gr%zm, rho_ds_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_thv_ds_zm, stat_file_zm, "thv_ds_zm", gr%nz, timestep, &
             gr%zm, thv_ds_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Kh_zm, stat_file_zm, "Kh_zm", gr%nz, timestep, &
             gr%zm, Kh_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_tau_zm, stat_file_zm, "tau_zm", gr%nz, timestep, &
             gr%zm, tau_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_up2, stat_file_zm, "up2", gr%nz, timestep, &
             gr%zm, up2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_vp2, stat_file_zm, "vp2", gr%nz, timestep, &
             gr%zm, vp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_sigma_sqd_w, stat_file_zm, "sigma_sqd_w", gr%nz, timestep, &
             gr%zm, sigma_sqd_w, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


!-----------------------------------------------------------

      call get_clubb_variable_interpolated &
           ( input_veg_T_in_K, stat_file_sfc, "veg_T_in_K", 1, timestep, &
             (/0._core_rknd/), tmp1(1), l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error
      veg_T_in_K = tmp1(1)

!     if ( input_veg_T_in_K ) then
!       call get_grads_var( fread_var, "veg_T_in_K", & 
!                     timestep, & 
!                     tmp1, l_read_error )
!       l_fatal_error = l_fatal_error .or. l_read_error
!       veg_T_in_K = tmp1(1)
!       print *, "Veg T = ", veg_T_in_K
!     endif
      call get_clubb_variable_interpolated &
           ( input_deep_soil_T_in_K, stat_file_sfc, "deep_soil_T_in_K", 1, timestep, &
             (/0._core_rknd/), tmp1(1), l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error
      deep_soil_T_in_K = tmp1(1)

!     if ( input_deep_soil_T_in_K ) then
!       call get_grads_var( fread_var, "deep_soil_T_in_", & 
!                     timestep, & 
!                     tmp1, l_read_error )
!       l_fatal_error = l_fatal_error .or. l_read_error
!       deep_soil_T_in_K = tmp1(1)
!       print *,"Deep soil = ",deep_soil_T_in_K
!     endif
      call get_clubb_variable_interpolated &
           ( input_sfc_soil_T_in_K, stat_file_sfc, "sfc_soil_T_in_K", 1, timestep, &
             (/0._core_rknd/), tmp1(1), l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error
      sfc_soil_T_in_K = tmp1(1)

!     if ( input_sfc_soil_T_in_K ) then
!       call get_grads_var( fread_var, "sfc_soil_T_in_K", & 
!                     timestep, & 
!                     tmp1, l_read_error )
!       l_fatal_error = l_fatal_error .or. l_read_error
!       sfc_soil_T_in_K = tmp1(1)
!       print *,"surface_soil = ", sfc_soil_T_in_K
!     endif


      if ( l_fatal_error ) stop "oops, get_grads_var failed in stat_fields_reader"

    case ( coamps_input_type ) ! COAMPS LES stats data     
      unit_number = 15

      ! stats_sm
      call open_grads_read( unit_number, stat_file_zt,  & 
                            fread_var, l_read_error )

      if ( l_read_error ) then
        write(fstderr,*) "Error reading file "// trim( stat_file_zt )
        stop "Fatal error"
      end if

      l_fatal_error = .false.

      ! Temporarily store LES output in variable array LES_tmp1.
      ! Allocate LES_tmp1 based on lowest and highest vertical indices of LES
      ! output.
      allocate( LES_tmp1(fread_var%ia:fread_var%iz) )

      ! Find the lowest and highest indices of CLUBB thermodynamic levels that
      ! fall within the domain of the LES output.
      call CLUBB_levels_within_LES_domain( fread_var, gr%zt,  &
                                           k_lowest_zt_input, k_highest_zt_input )

      allocate( exact_lev_idx_zt(k_lowest_zt_input:k_highest_zt_input) )
      allocate( lower_lev_idx_zt(k_lowest_zt_input:k_highest_zt_input) )
      allocate( upper_lev_idx_zt(k_lowest_zt_input:k_highest_zt_input) )
      allocate( l_lin_int_zt(k_lowest_zt_input:k_highest_zt_input) )

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

      ! Find the lowest and highest indices of CLUBB momentum levels that fall
      ! within the domain of the LES output.
      call CLUBB_levels_within_LES_domain( fread_var, gr%zm,  &
                                           k_lowest_zm_input, k_highest_zm_input )

      allocate( exact_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( lower_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( upper_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( l_lin_int_zm(k_lowest_zm_input:k_highest_zm_input) )

      ! For all CLUBB momentum levels, k, that are within the LES domain, find
      ! either the index of the LES level that exactly matches the altitude of the
      ! CLUBB level, or find the two indices of the LES levels that are on either
      ! side of the CLUBB level.
      do k = k_lowest_zm_input, k_highest_zm_input, 1
        ! CLUBB vertical level k is found at an altitude that is within the
        ! domain of the LES output.
        call LES_grid_to_CLUBB_grid( fread_var, gr%zm, k,  &
                                     exact_lev_idx_zm(k), lower_lev_idx_zm(k),  &
                                     upper_lev_idx_zm(k), l_lin_int_zm(k) )
      enddo


      ! Initialize l_read_error for case ( "les" )
      l_read_error = .false.

      call get_coamps_variable_interp( &
              input_um, fread_var, "um", timestep, gr%nz, &               ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              um, l_read_error )                           ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_um ) then
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

      call get_coamps_variable_interp( &
              input_vm, fread_var, "vm", timestep, gr%nz, &               ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              vm, l_read_error )                           ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_vm ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, use the
        ! values of vm at thermodynamic levels 3 and 2 to find the value at
        ! thermodynamic level 1 through the use of a linear extension.
        if ( k_lowest_zt_input == 2 ) then
          vm(1)  & 
          = lin_ext_zt_bottom( vm(3), vm(2), & 
                               gr%zt(3), gr%zt(2), gr%zt(1) )
        endif
      endif


      call get_coamps_variable_interp( &
              input_rtm, fread_var, "qtm", timestep, gr%nz, &             ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              rtm, l_read_error )                          ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_rtm ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of rtm at thermodynamic level 1 to the value at thermodynamic
        ! level 2, as it is done in advance_xm_wpxp.
        if ( k_lowest_zt_input == 2 ) then
          rtm(1) = rtm(2)
        endif
      endif


      call get_coamps_variable_interp( &
              input_thlm, fread_var, "thlm", timestep, gr%nz, &           ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              thlm, l_read_error )                         ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_thlm ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of thlm at thermodynamic level 1 to the value at thermodynamic
        ! level 2, as it is done in advance_xm_wpxp.
        if ( k_lowest_zt_input == 2 ) then
          thlm(1) = thlm(2)
        endif
      endif

      ! We obtain wp2 from stats_sw

      call get_coamps_variable_interp( &
              input_wp3, fread_var, "wp3", timestep, gr%nz, &             ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              wp3, l_read_error )                          ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_wp3 ) then
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of wp3 at thermodynamic level 1 to 0, as it is done in
        ! advance_wp2_wp3.
        if ( k_lowest_zt_input == 2 ) then
          wp3(1) = 0.0_core_rknd
        endif
      endif

      call get_coamps_variable_interp( &
              input_wprtp, fread_var, "wpqtp", timestep, gr%nz, &         ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              wprtp, l_read_error )                        ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_wprtp ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wprtp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like surface_varnce.
        if ( k_lowest_zm_input == 2 ) then
          wprtp(1)  & 
          = lin_ext_zm_bottom( wprtp(3), wprtp(2), & 
                               gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      call get_coamps_variable_interp( &
              input_wpthlp, fread_var, "wpthltp", timestep, gr%nz, &      ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              wpthlp, l_read_error )                       ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_wpthlp ) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wpthlp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like surface_varnce.
        if ( k_lowest_zm_input == 2 ) then
          wpthlp(1)  & 
          = lin_ext_zm_bottom( wpthlp(3), wpthlp(2), & 
                               gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif


      call get_coamps_variable_interp( &
              input_rtp2, fread_var, "qtp2", timestep, gr%nz, &           ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              rtp2, l_read_error )                         ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_rtp2 ) then
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

      call get_coamps_variable_interp( &
              input_thlp2, fread_var, "thlp2", timestep, gr%nz, &         ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              thlp2, l_read_error )                        ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_thlp2 ) then
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


      call get_coamps_variable_interp( &
              input_rtpthlp, fread_var, "qtpthlp", timestep, gr%nz, &     ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              rtpthlp, l_read_error )                      ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_rtpthlp) then
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! rtpthlp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like surface_varnce.
        if ( k_lowest_zm_input == 2 ) then
          rtpthlp(1)  & 
          = lin_ext_zm_bottom( rtpthlp(3), rtpthlp(2), & 
                               gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      ! upwp/vpwp in stats_sw
      if ( input_ug ) then
        write(fstderr,*) "The variable ug is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_vg ) then
        write(fstderr,*) "The variable vg is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if


      call get_coamps_variable_interp( &
              input_rcm, fread_var, "qcm", timestep, gr%nz, &             ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              rcm, l_read_error )                          ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_coamps_variable_interp( &
              input_wm_zt, fread_var, "wlsm", timestep, gr%nz, &          ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              wm_zt, l_read_error )                        ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_coamps_variable_interp( &
              input_exner, fread_var, "ex0", timestep, gr%nz, &           ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              exner, l_read_error )                        ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_coamps_variable_interp( &
              input_em, fread_var, "em", timestep, gr%nz, &               ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              em, l_read_error )                           ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      temp = 0.0_core_rknd ! Initialize temp to 0.0

      call get_coamps_variable_interp( &
              input_em, fread_var, "tke", timestep, gr%nz, &              ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              temp, l_read_error )                         ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if( input_em ) then
        do k=k_lowest_zm_input, k_highest_zm_input, 1
          em(k) = em(k) + temp(k)
        end do

        where ( em < em_min ) em = em_min
      end if



      call get_coamps_variable_interp( &
              input_p, fread_var, "pm", timestep, gr%nz, &                ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              p_in_Pa, l_read_error )                      ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error



      call get_coamps_variable_interp( &
              input_rho, fread_var, "dn0", timestep, gr%nz, &             ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              rho, l_read_error )                          ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      ! rho_zm is in stats_sw

      if ( input_Lscale ) then
        write(fstderr,*) "The variable Lscale is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      endif

      if ( input_Lscale_up ) then
        write(fstderr,*) "The variable Lscale_up is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_Lscale_down ) then
        write(fstderr,*) "The variable Lscale_down is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if



      call get_coamps_variable_interp( &
              input_Kh_zt, fread_var, "kh", timestep, gr%nz, &            ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              Kh_zt, l_read_error )                        ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error


      if ( input_Kh_zm ) then
        write(fstderr,*) "The variable Kh_zm is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_tau_zt ) then
        write(fstderr,*) "The variable tau_zt is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_tau_zm ) then
        write(fstderr,*) "The variable tau_zm is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if


      call get_coamps_variable_interp( &
              input_wpthvp, fread_var, "wpthvp", timestep, gr%nz, &       ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              wpthvp, l_read_error )                       ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error


      if ( input_thl1 ) then
        write(fstderr,*) "The variable thl1 is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_thl2 ) then
        write(fstderr,*) "The variable thl2 is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_mixt_frac ) then
        write(fstderr,*) "The variable mixt_frac is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_s1 ) then
        write(fstderr,*) "The variable s1 is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_s2 ) then
        write(fstderr,*) "The variable s2 is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_stdev_s1 ) then
        write(fstderr,*) "The variable stdev_s1 is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_stdev_s2 ) then
        write(fstderr,*) "The variable stdev_s2 is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_rc1 ) then
        write(fstderr,*) "The variable rc1 is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_rc2 ) then
        write(fstderr,*) "The variable rc2 is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if


      call get_coamps_variable_interp( &
              input_thvm, fread_var, "thvm", timestep, gr%nz, &           ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              thvm, l_read_error )                         ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error


      if ( input_rrainm ) then
        if ( iirrainm < 1 ) then
            write(fstderr,*) "Rain water mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          call get_coamps_variable_interp( &
                  input_rrainm, fread_var, "qrm", timestep, gr%nz, &          ! Intent(in)
                  gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
                  upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
                  hydromet(:,iirrainm), l_read_error )         ! Intent(in/out), Intent(out)

          l_fatal_error = l_fatal_error .or. l_read_error
        end if
      end if

      if ( input_Nrm ) then
        if ( iiNrm < 1 ) then
            write(fstderr,*) "Rain droplet number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          call get_coamps_variable_interp( &
                  input_Nrm, fread_var, "nrm", timestep, gr%nz, &          ! Intent(in)
                  gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
                  upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
                  hydromet(:,iiNrm), l_read_error )         ! Intent(in/out), Intent(out)

          l_fatal_error = l_fatal_error .or. l_read_error
        end if
      end if


      if ( input_Ncm ) then
        if ( iiNcm < 1 ) then
            write(fstderr,*) "Cloud droplet number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          call get_coamps_variable_interp( &
                  input_Ncm, fread_var, "ncm", timestep, gr%nz, &          ! Intent(in)
                  gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
                  upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
                  hydromet(:,iiNcm), l_read_error )         ! Intent(in/out), Intent(out)

          l_fatal_error = l_fatal_error .or. l_read_error
        end if
      end if


      if ( input_rsnowm ) then
        if ( iirsnowm < 1 ) then
            write(fstderr,*) "Snow mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          call get_coamps_variable_interp( &
                  input_rsnowm, fread_var, "qsm", timestep, gr%nz, &          ! Intent(in)
                  gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
                  upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
                  hydromet(:,iirsnowm), l_read_error )         ! Intent(in/out), Intent(out)

          l_fatal_error = l_fatal_error .or. l_read_error
        end if
      end if


      if ( input_ricem ) then
        if ( iiricem < 1 ) then
            write(fstderr,*) "Ice mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          call get_coamps_variable_interp( &
                  input_ricem, fread_var, "qim", timestep, gr%nz, &          ! Intent(in)
                  gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
                  upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
                  hydromet(:,iiricem), l_read_error )         ! Intent(in/out), Intent(out)

          l_fatal_error = l_fatal_error .or. l_read_error
        end if
      end if


      if ( input_rgraupelm ) then
        if ( iirgraupelm < 1 ) then
            write(fstderr,*) "Graupel mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          call get_coamps_variable_interp( &
                  input_rgraupelm, fread_var, "qgm", timestep, gr%nz, &          ! Intent(in)
                  gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
                  upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
                  hydromet(:,iirgraupelm), l_read_error )         ! Intent(in/out), Intent(out)

          l_fatal_error = l_fatal_error .or. l_read_error
        end if
      end if


      call get_coamps_variable_interp( &
              input_Ncnm, fread_var, "ncnm", timestep, gr%nz, &           ! Intent(in)
              gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
              upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
              Ncnm, l_read_error )                         ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error


      if ( input_Nim ) then
        if ( iiNim < 1 ) then
            write(fstderr,*) "Ice number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
        else
          call get_coamps_variable_interp( &
                  input_Nim, fread_var, "nim", timestep, gr%nz, &          ! Intent(in)
                  gr%zt, k_lowest_zt_input, k_highest_zt_input, l_lin_int_zt, & ! Intent(in)
                  upper_lev_idx_zt, lower_lev_idx_zt, exact_lev_idx_zt, &       ! Intent(in)
                  hydromet(:,iiNim), l_read_error )         ! Intent(in/out), Intent(out)

          l_fatal_error = l_fatal_error .or. l_read_error
        end if
      end if


      if ( input_thlm_forcing ) then
        write(fstderr,*) "The variable thlm_forcing is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_rtm_forcing ) then
        write(fstderr,*) "The variable rtm_forcing is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if


      call get_coamps_variable_interp( &
              input_up2, fread_var, "up2", timestep, gr%nz, &             ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              up2, l_read_error )                          ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_up2 ) then
        ! Clip up2 to be no smaller than w_tol_sqd
        where ( up2 < w_tol_sqd ) up2 = w_tol_sqd
      endif


      call get_coamps_variable_interp( &
              input_vp2, fread_var, "vp2", timestep, gr%nz, &             ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              vp2, l_read_error )                          ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_vp2 ) then
         where ( vp2 < w_tol_sqd ) vp2 = w_tol_sqd
      endif

      if ( input_sigma_sqd_w ) then
        write(fstderr,*) "The variable sigma_sqd_w is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_sigma_sqd_w_zt ) then
        write(fstderr,*) "The variable sigma_sqd_w_zt is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_cloud_frac ) then
        write(fstderr,*) "The variable cloud_frac is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_veg_T_in_K ) then
        write(fstderr,*) "The variable veg_T_in_K is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_deep_soil_T_in_K ) then
        write(fstderr,*) "The variable deep_soil_T_in_K is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_sfc_soil_T_in_K ) then
        write(fstderr,*) "The variable sfc_soil_T_in_K is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( input_radht ) then
        write(fstderr,*) "The variable radht is not setup for input_type" &
          //trim( input_type )
        l_fatal_error = .true.
      end if

      if ( l_fatal_error ) stop "oops, get_grads_var failed in stat_fields_reader"

      ! Deallocate temporary storage variable LES_tmp1.
      deallocate( LES_tmp1 )

      deallocate( exact_lev_idx_zt )
      deallocate( lower_lev_idx_zt )
      deallocate( upper_lev_idx_zt )
      deallocate( l_lin_int_zt )

      deallocate( exact_lev_idx_zm )
      deallocate( lower_lev_idx_zm )
      deallocate( upper_lev_idx_zm )
      deallocate( l_lin_int_zm )

      call close_grads_read( fread_var )

    end select

    select case ( stats_input_type )

    case ( coamps_input_type )
      unit_number = 15

      ! stats_sw
      call open_grads_read( unit_number, stat_file_zm, & ! In
                            fread_var, & ! In/Out
                            l_read_error ) ! Out

      if ( l_read_error ) then
        write(fstderr,*) "Error reading file "// trim( stat_file_zm )
        stop "Fatal error"
      end if

      ! Temporarily store LES output in variable array LES_tmp1.
      ! Allocate LES_tmp1 based on lowest and highest vertical indices of LES
      ! output.
      allocate( LES_tmp1(fread_var%ia:fread_var%iz) )

      ! Find the lowest and highest indices of CLUBB momentum levels that fall
      ! within the domain of the LES output.
      call CLUBB_levels_within_LES_domain( fread_var, gr%zm,  &
                                           k_lowest_zm_input, k_highest_zm_input )

      allocate( exact_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( lower_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( upper_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( l_lin_int_zm(k_lowest_zm_input:k_highest_zm_input) )

      ! For all CLUBB momentum levels, k, that are within the LES domain, find
      ! either the index of the LES level that exactly matches the altitude of the
      ! CLUBB level, or find the two indices of the LES levels that are on either
      ! side of the CLUBB level.
      do k = k_lowest_zm_input, k_highest_zm_input, 1
        ! CLUBB vertical level k is found at an altitude that is within the
        ! domain of the LES output.
        call LES_grid_to_CLUBB_grid( fread_var, gr%zm, k,  &
                                     exact_lev_idx_zm(k), lower_lev_idx_zm(k),  &
                                     upper_lev_idx_zm(k), l_lin_int_zm(k) )
      enddo


      ! Note:  l_fatal_error has already been initialized

      ! Note:  wpup_sgs and wpvp_sgs must be added to make the u'w' and v'w' terms
      !        as they are in CLUBB.

      call get_coamps_variable_interp( &
              input_upwp, fread_var, "wpup", timestep, gr%nz, &           ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              upwp, l_read_error )                         ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      temp = 0.0_core_rknd  ! clear temp
      call get_coamps_variable_interp( &
              input_upwp, fread_var, "wpup_sgs", timestep, gr%nz, &       ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              temp, l_read_error )                         ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_upwp) then
        do k=k_lowest_zm_input, k_highest_zm_input, 1
          upwp(k) = temp(k)
        end do
      endif


      call get_coamps_variable_interp( &
              input_vpwp, fread_var, "wpvp", timestep, gr%nz, &           ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              vpwp, l_read_error )                         ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      temp = 0.0_core_rknd  ! clear temp
      call get_coamps_variable_interp( &
              input_vpwp, fread_var, "wpvp_sgs", timestep, gr%nz, &       ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              temp, l_read_error )                         ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_vpwp) then
        do k=k_lowest_zm_input, k_highest_zm_input, 1
          vpwp(k) = temp(k)
        end do
      endif


      call get_coamps_variable_interp( &
              input_wp2, fread_var, "wp2", timestep, gr%nz, &            ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              wp2, l_read_error )                          ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error

      if ( input_wp2 ) then
         if ( any( wp2(1:gr%nz) < w_tol_sqd ) ) then
          do k=1, gr%nz
            wp2(k) = max( wp2(k), w_tol_sqd )
          enddo
        endif
      endif


      call get_coamps_variable_interp( &
              input_rho_zm, fread_var, "dn0", timestep, gr%nz, &          ! Intent(in)
              gr%zm, k_lowest_zm_input, k_highest_zm_input, l_lin_int_zm, & ! Intent(in)
              upper_lev_idx_zm, lower_lev_idx_zm, exact_lev_idx_zm, &       ! Intent(in)
              rho_zm, l_read_error )                       ! Intent(in/out), Intent(out)

      l_fatal_error = l_fatal_error .or. l_read_error


      if ( l_fatal_error ) stop "get_grads_var failed for stats_sw in stat_fields_reader"

      ! Deallocate temporary storage variable LES_tmp1.
      deallocate( LES_tmp1 )

      deallocate( exact_lev_idx_zm )
      deallocate( lower_lev_idx_zm )
      deallocate( upper_lev_idx_zm )
      deallocate( l_lin_int_zm )

      call close_grads_read( fread_var )

    end select

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
      call open_netcdf_read( 'thlm', trim( filename ), & ! In
                             fread_var, & ! In/Out
                             l_error ) ! Out
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
    nearest_timestep = nint( delta_time / sec_per_min )

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
  subroutine get_coamps_variable_interp( &
                 l_input_var, fread_var, var_name, timestep, vardim, &
                 clubb_heights, k_lowest, k_highest, l_lin_int, &
                 upper_lev_idx, lower_lev_idx, exact_lev_idx, &
                 variable_interpolated, l_error )

  ! Description:
  !   Obtain a profile of a COAMPS variable from a GrADS file and interpolate if
  !   needed.
  ! References:
  !   None
  !--------------------------------------------------------------------------------


    use stat_file_module, only: stat_file

    use input_grads, only: get_grads_var

    use interpolation, only: lin_int

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    logical, intent(in) :: &
      l_input_var ! If .true., this variable needs to be read in.

    type (stat_file), intent(in) :: &
      fread_var ! The grads data.

    character(len=*), intent(in) :: &
      var_name ! The name of the variable to read and interpolate

    integer, intent(in) :: &
      timestep, & ! The current timestep
      vardim, &   ! The dimension of the variable
      k_lowest, & ! The lowest level from the LES output
      k_highest   ! The highest level from the LES output

    real( kind = core_rknd ), dimension(vardim), intent(in) :: &
      clubb_heights ! The heights for each level

    logical, dimension(k_lowest:k_highest), intent(in) :: &
      l_lin_int ! Flag that is turned on if linear interpolation is needed.

    integer, dimension(k_lowest:k_highest), intent(in) :: &
      upper_lev_idx, &
      lower_lev_idx, &
      exact_lev_idx

    real( kind = core_rknd ), dimension(vardim), intent(inout) :: &
      variable_interpolated ! The resulting interpolated variable
  
    logical, intent(out) :: &
      l_error ! Error flag

    ! Local variables

    real( kind = core_rknd ), dimension(fread_var%ia:fread_var%iz) :: &
      LES_tmp ! Temporary variable

    integer :: &
      k ! index variable

    !--------------------------- BEGIN CODE ---------------------------------

    l_error = .false. ! Initialize to false

    if ( l_input_var ) then
      call get_grads_var( fread_var, var_name, timestep, &
                    LES_tmp(fread_var%ia:fread_var%iz), l_error )
      !LES_tmp is the value of the variable from the LES GrADS file.
      do k = k_lowest, k_highest, 1
        if( l_lin_int(k) ) then
          ! CLUBB level k is found at an altitude that is between two
          ! LES levels.  Linear interpolation is required.
          variable_interpolated(k) = lin_int( clubb_heights(k), &
                                        fread_var%z(upper_lev_idx(k)), &
                                        fread_var%z(lower_lev_idx(k)), &
                                        LES_tmp(upper_lev_idx(k)), &
                                        LES_tmp(lower_lev_idx(k)) )
        else
          ! CLUBB level k is found at an altitude that is an exact
          ! match with an LES level altitude.
          variable_interpolated(k) = LES_tmp(exact_lev_idx(k))
        end if
      end do

    end if

  end subroutine get_coamps_variable_interp
        


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
      input_um, input_vm, input_rtm, input_thlm, & 
      input_wp2, input_wprtp, input_wpthlp,  & 
      input_wp3, input_rtp2, input_thlp2,  & 
      input_rtpthlp, input_upwp, input_vpwp, & 
      input_ug, input_vg, input_rcm,  & 
      input_wm_zt, input_exner, input_em, & 
      input_p, input_rho, input_rho_zm, & 
      input_Lscale, input_Lscale_up, input_Lscale_down, & 
      input_Kh_zt, input_Kh_zm, input_tau_zm, input_tau_zt, & 
      input_wpthvp, input_radht, &
      input_thl1, input_thl2, input_mixt_frac, input_s1, input_s2, &
      input_stdev_s1, input_stdev_s2, input_rc1, input_rc2, &
      input_thvm, input_rrainm,input_Nrm,  & 
      input_rsnowm, input_ricem, input_rgraupelm,  & 
      input_thlm_forcing, input_rtm_forcing, & 
      input_up2, input_vp2, input_sigma_sqd_w, input_Ncm,  & 
      input_Ncnm, input_Nim, input_cloud_frac, input_sigma_sqd_w_zt, &
      input_veg_T_in_K, input_deep_soil_T_in_K, &
      input_sfc_soil_T_in_K

    ! --- Begin Code ---

    ! Pick some initial values
    datafile = ''
    input_type = 'clubb'

    input_um = .false.
    input_vm = .false.
    input_rtm = .false.
    input_thlm  = .false.
    input_wp2 = .false.
    input_wprtp = .false.
    input_wpthlp = .false.
    input_wp3 = .false.
    input_rtp2 = .false.
    input_thlp2 = .false.
    input_rtpthlp = .false.
    input_upwp = .false.
    input_vpwp = .false.
    input_ug = .false.
    input_vg = .false.
    input_rcm = .false.
    input_wm_zt = .false.
    input_exner = .false.
    input_em = .false.
    input_p = .false.
    input_rho = .false.
    input_Lscale = .false.
    input_Lscale_up = .false.
    input_Lscale_down = .false.
    input_Kh_zt = .false.
    input_Kh_zm = .false.
    input_tau_zm = .false.
    input_tau_zt = .false.
    input_wpthvp = .false.
    input_thl1 = .false.
    input_thl2 = .false.
    input_mixt_frac = .false.
    input_s1 = .false.
    input_s2 = .false.
    input_stdev_s1 = .false.
    input_stdev_s2 = .false.
    input_rc1 = .false.
    input_rc2 = .false.
    input_thvm = .false.
    input_rrainm = .false.
    input_Nrm = .false.
    input_rsnowm = .false.
    input_ricem = .false.
    input_rgraupelm = .false.
    input_thlm_forcing = .false.
    input_rtm_forcing = .false.
    input_up2 = .false.
    input_vp2 = .false.
    input_sigma_sqd_w = .false.
    input_Ncm = .false.
    input_Ncnm = .false.
    input_Nim = .false.
    input_cloud_frac = .false.
    input_sigma_sqd_w_zt = .false.
    input_veg_T_in_K = .false.
    input_deep_soil_T_in_K = .false.
    input_sfc_soil_T_in_K  = .false.

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
