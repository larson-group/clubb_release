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

  logical, public :: input_um, input_vm, input_rtm, input_thlm, & 
                     input_wp2, input_wprtp, input_wpthlp,  & 
                     input_wp3, input_rtp2, input_thlp2,  & 
                     input_rtpthlp, input_upwp, input_vpwp, & 
                     input_ug, input_vg, input_rcm,  & 
                     input_wm_zt, input_exner, input_em, & 
                     input_p, input_rho, input_rho_zm, & 
                     input_Lscale, input_Lscale_up, input_Lscale_down, & 
                     input_Kh_zt, input_Kh_zm, input_tau_zm, input_tau_zt, & 
                     input_wpthvp, &
                     input_thl1, input_thl2, input_a, input_s1, input_s2, &
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

    use constants, only: fstderr, fstdout ! Constants

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
        wpthvp, & 
        rtp2, & 
        thlp2, & 
        rtpthlp, & 
        upwp, & 
        vpwp, & 
        Kh_zt, & 
        p_in_Pa, & 
        exner, & 
        rcm, & 
        wm_zt, & 
        rho, & 
        rho_zm, & 
        thlm_forcing, & 
        rtm_forcing, & 
        cloud_frac, & 
        tau_zm, & 
        up2, & 
        vp2, & 
        sigma_sqd_w

    use variables_diagnostic_module, only: & 
        hydromet,  & ! Variable(s)
        tau_zt, & 
        ug, & 
        vg, & 
        Lscale, & 
        Lscale_up, & 
        Lscale_down, & 
        Kh_zm, & 
        thvm, & 
        Ncnm, & 
        sigma_sqd_w_zt, & 
        em

    use variables_prognostic_module, only: & 
        pdf_params ! Variable(s)

    use grid_class, only: & 
        gr,  & ! Variable(s)
        zt2zm ! Procedure(s)

    use constants, only:  &
        rttol,    & ! Variable(s)
        thltol,   &
        wtol_sqd, &
        emin,     &
        fstderr

    use array_index, only:  & 
        iirrainm, iiNrm, iirsnowm, iiricem, iirgraupelm, iiNim, iiNcm

    use stat_file_module, only: & 
        stat_file     ! Type

    use stat_file_utils, only: & 
       LES_grid_to_CLUBB_grid, & ! Procedure(s)
       CLUBB_levels_within_LES_domain

    use inputfile_class, only: & 
        get_grads_var,  & ! Procedure(s)
        open_grads_read, & 
        close_grads_read

    use interpolation, only: &
        lin_int ! Procedure(s)

    use extrapolation, only: &
      lin_ext_zt_bottom, &
      lin_ext_zm_bottom

    use parameters_microphys, only: &
      micro_scheme ! Variable(s)

    use soil_vegetation, only: deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K

    implicit none

    ! External
    intrinsic :: max, trim, any

    ! Arguments
    integer, intent(in) :: timestep

    ! Local Variables
    logical :: l_read_error, l_fatal_error

    type (stat_file) :: fread_var

    real, dimension(:), allocatable :: LES_tmp1

    real, dimension(gr%nnzp+1) :: tmp1

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

    integer :: k  ! Array index

    ! ---- Begin Code ----

    select case ( stats_input_type )

    case ( clubb_input_type )

      !  Thermo grid - zt file

      ! Initialize l_fatal_error for case clubb_input_type
      l_fatal_error = .false.

      call get_clubb_variable_interpolated &
           ( input_um, stat_file_zt, "um", gr%nnzp, timestep, &
             gr%zt, um, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_vm, stat_file_zt, "vm", gr%nnzp, timestep, &
             gr%zt, vm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rtm, stat_file_zt, "rtm", gr%nnzp, timestep, &
             gr%zt, rtm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( input_thlm, stat_file_zt, "thlm", gr%nnzp, timestep, &
             gr%zt, thlm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_wp3, stat_file_zt, "wp3", gr%nnzp, timestep, &
             gr%zt, wp3, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_tau_zt, stat_file_zt, "tau_zt", gr%nnzp, timestep, &
             gr%zt, tau_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rrainm, stat_file_zt, "rrainm", gr%nnzp, timestep, & 
             gr%zt, tmp1(1:gr%nnzp), l_read_error )
      if ( input_rrainm ) then
        hydromet(1:gr%nnzp,iirrainm) = tmp1(1:gr%nnzp)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rsnowm, stat_file_zt, "rsnowm", gr%nnzp, timestep, & 
             gr%zt, tmp1(1:gr%nnzp), l_read_error )
      if ( input_rsnowm ) then
        hydromet(1:gr%nnzp,iirsnowm) = tmp1(1:gr%nnzp)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_ricem, stat_file_zt, "ricem", gr%nnzp, timestep, & 
             gr%zt, tmp1(1:gr%nnzp), l_read_error )
      if ( input_ricem ) then
        hydromet(1:gr%nnzp,iiricem) = tmp1(1:gr%nnzp)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rgraupelm, stat_file_zt, "rgraupelm", gr%nnzp, timestep, & 
             gr%zt, tmp1(1:gr%nnzp), l_read_error )
      if ( input_rgraupelm ) then
        hydromet(1:gr%nnzp,iirgraupelm) = tmp1(1:gr%nnzp)
      end if
      l_fatal_error = l_fatal_error .or. l_read_error

!--------------------------------------------------------
! Added variables for clubb_restart
      call get_clubb_variable_interpolated &
           ( input_p, stat_file_zt, "p_in_Pa", gr%nnzp, timestep, &
             gr%zt, p_in_Pa, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_exner, stat_file_zt, "exner", gr%nnzp, timestep, &
             gr%zt, exner, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_ug, stat_file_zt, "ug", gr%nnzp, timestep, &
             gr%zt, ug, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_vg, stat_file_zt, "vg", gr%nnzp, timestep, &
             gr%zt, vg, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rcm, stat_file_zt, "rcm", gr%nnzp, timestep, &
             gr%zt, rcm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_wm_zt, stat_file_zt, "wm", gr%nnzp, timestep, &
             gr%zt, wm_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


      call get_clubb_variable_interpolated &
           ( input_rho, stat_file_zt, "rho", gr%nnzp, timestep, &
             gr%zt, rho, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Lscale, stat_file_zt, "Lscale", gr%nnzp, timestep, &
             gr%zt, Lscale, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Lscale_up, stat_file_zt, "Lscale_up", gr%nnzp, timestep, &
             gr%zt, Lscale_up, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Lscale_down, stat_file_zt, "Lscale_down", gr%nnzp, timestep, &
             gr%zt, Lscale_down, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Kh_zt, stat_file_zt, "Kh_zt", gr%nnzp, timestep, &
             gr%zt, Kh_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_thvm, stat_file_zt, "thvm", gr%nnzp, timestep, &
             gr%zt, thvm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_thlm_forcing, stat_file_zt, "thlm_f", gr%nnzp, timestep, &
             gr%zt, thlm_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rtm_forcing, stat_file_zt, "rtm_f", gr%nnzp, timestep, &
             gr%zt, rtm_forcing, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Ncm, stat_file_zt, "Ncm", gr%nnzp, timestep, &
             gr%zt, tmp1(1:gr%nnzp), l_read_error )
      if ( input_Ncm ) then
        hydromet(1:gr%nnzp,iiNcm) = tmp1(1:gr%nnzp)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Ncnm, stat_file_zt, "Ncnm", gr%nnzp, timestep, &
             gr%zt, Ncnm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Nim, stat_file_zt, "Nim", gr%nnzp, timestep, &
             gr%zt, tmp1(1:gr%nnzp), l_read_error )
      if ( input_Nim ) then
        hydromet(1:gr%nnzp, iiNcm) = tmp1(1:gr%nnzp)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_cloud_frac, stat_file_zt, "cloud_frac", gr%nnzp, timestep, &
             gr%zt, cloud_frac, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Nrm, stat_file_zt, "Nrm", gr%nnzp, timestep, &
             gr%zt, tmp1(1:gr%nnzp), l_read_error )
      if ( input_Nrm ) then
        hydromet(1:gr%nnzp, iiNrm) = tmp1(1:gr%nnzp)
      end if

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_sigma_sqd_w_zt, stat_file_zt, "sigma_sqd_w_zt", gr%nnzp, timestep, &
             gr%zt, sigma_sqd_w_zt, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      ! PDF Parameters (needed for K&K microphysics)
      call get_clubb_variable_interpolated &
           ( input_thl1, stat_file_zt, "thl1", gr%nnzp, timestep, &
             gr%zt, pdf_params%thl1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_thl2, stat_file_zt, "thl2", gr%nnzp, timestep, &
             gr%zt, pdf_params%thl2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_a, stat_file_zt, "a", gr%nnzp, timestep, &
             gr%zt, pdf_params%a, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_s1, stat_file_zt, "s1", gr%nnzp, timestep, &
             gr%zt, pdf_params%s1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_s2, stat_file_zt, "s2", gr%nnzp, timestep, &
             gr%zt, pdf_params%s2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_stdev_s1, stat_file_zt, "stdev_s1", gr%nnzp, timestep, &
             gr%zt, pdf_params%stdev_s1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_stdev_s2, stat_file_zt, "stdev_s2", gr%nnzp, timestep, &
             gr%zt, pdf_params%stdev_s2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rc1, stat_file_zt, "rc1", gr%nnzp, timestep, &
             gr%zt, pdf_params%rc1, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rc2, stat_file_zt, "rc2", gr%nnzp, timestep, &
             gr%zt, pdf_params%rc2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

!--------------------------------------------------------

      ! Read in the zm file

      call get_clubb_variable_interpolated &
           ( input_wp2, stat_file_zm, "wp2", gr%nnzp, timestep, &
             gr%zm, wp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_wprtp, stat_file_zm, "wprtp", gr%nnzp, timestep, &
             gr%zm, wprtp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_wpthlp, stat_file_zm, "wpthlp", gr%nnzp, timestep, &
             gr%zm, wpthlp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_wpthvp, stat_file_zm, "wpthvp", gr%nnzp, timestep, &
             gr%zm, wpthvp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rtp2, stat_file_zm, "rtp2", gr%nnzp, timestep, &
             gr%zm, rtp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_thlp2, stat_file_zm, "thlp2", gr%nnzp, timestep, &
             gr%zm, thlp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rtpthlp, stat_file_zm, "rtpthlp", gr%nnzp, timestep, &
             gr%zm, rtpthlp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_upwp, stat_file_zm, "upwp", gr%nnzp, timestep, &
             gr%zm, upwp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_vpwp, stat_file_zm, "vpwp", gr%nnzp, timestep, &
             gr%zm, vpwp, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

!-----------------------------------------------------------
      call get_clubb_variable_interpolated &
           ( input_em, stat_file_zm, "em", gr%nnzp, timestep, &
             gr%zm, em, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_rho_zm, stat_file_zm, "rho_zm", gr%nnzp, timestep, &
             gr%zm, rho_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_Kh_zm, stat_file_zm, "Kh_zm", gr%nnzp, timestep, &
             gr%zm, Kh_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_tau_zm, stat_file_zm, "tau_zm", gr%nnzp, timestep, &
             gr%zm, tau_zm, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_up2, stat_file_zm, "up2", gr%nnzp, timestep, &
             gr%zm, up2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_vp2, stat_file_zm, "vp2", gr%nnzp, timestep, &
             gr%zm, vp2, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error

      call get_clubb_variable_interpolated &
           ( input_sigma_sqd_w, stat_file_zm, "sigma_sqd_w", gr%nnzp, timestep, &
             gr%zm, sigma_sqd_w, l_read_error )

      l_fatal_error = l_fatal_error .or. l_read_error


!-----------------------------------------------------------

      call get_clubb_variable_interpolated &
           ( input_veg_T_in_K, stat_file_sfc, "veg_T_in_K", 1, timestep, &
             (/0./), tmp1(1), l_read_error )

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
             (/0./), tmp1(1), l_read_error )

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
             (/0./), tmp1(1), l_read_error )

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

      ! stats_sm
      call open_grads_read( 15, stat_file_zt,  & 
                            fread_var )
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

      if ( input_um ) then
        call get_grads_var( fread_var, "um", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of um from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            um(k) = lin_int( gr%zt(k),  &
                             fread_var%z(upper_lev_idx_zt(k)),  &
                             fread_var%z(lower_lev_idx_zt(k)),  &
                             LES_tmp1(upper_lev_idx_zt(k)),  &
                             LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            um(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
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

      if ( input_vm ) then
        call get_grads_var( fread_var, "vm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of vm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            vm(k) = lin_int( gr%zt(k),  &
                             fread_var%z(upper_lev_idx_zt(k)),  &
                             fread_var%z(lower_lev_idx_zt(k)),  &
                             LES_tmp1(upper_lev_idx_zt(k)),  &
                             LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            vm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
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

      if ( input_rtm ) then
        call get_grads_var( fread_var, "qtm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rtm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            rtm(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            rtm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of rtm at thermodynamic level 1 to the value at thermodynamic
        ! level 2, as it is done in advance_xm_wpxp.
        if ( k_lowest_zt_input == 2 ) then
          rtm(1) = rtm(2)
        endif
      endif

      if ( input_thlm ) then
        call get_grads_var( fread_var, "thlm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of thlm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            thlm(k) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            thlm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of thlm at thermodynamic level 1 to the value at thermodynamic
        ! level 2, as it is done in advance_xm_wpxp.
        if ( k_lowest_zt_input == 2 ) then
          thlm(1) = thlm(2)
        endif
      endif

      ! We obtain wp2 from stats_sw

      if ( input_wp3 ) then
        call get_grads_var( fread_var, "wp3", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wp3 from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            wp3(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            wp3(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of wp3 at thermodynamic level 1 to 0, as it is done in
        ! advance_wp2_wp3.
        if ( k_lowest_zt_input == 2 ) then
          wp3(1) = 0.0
        endif
      endif

      if ( input_wprtp ) then
        call get_grads_var( fread_var, "wpqtp", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wprtp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            wprtp(k) = lin_int( gr%zm(k),  &
                                fread_var%z(upper_lev_idx_zm(k)),  &
                                fread_var%z(lower_lev_idx_zm(k)),  &
                                LES_tmp1(upper_lev_idx_zm(k)),  &
                                LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            wprtp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wprtp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like sfc_var.
        if ( k_lowest_zm_input == 2 ) then
          wprtp(1)  & 
          = lin_ext_zm_bottom( wprtp(3), wprtp(2), & 
                               gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if ( input_wpthlp ) then
        call get_grads_var( fread_var, "wpthlp", timestep,  &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wpthlp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            wpthlp(k) = lin_int( gr%zm(k),  &
                                 fread_var%z(upper_lev_idx_zm(k)),  &
                                 fread_var%z(lower_lev_idx_zm(k)),  &
                                 LES_tmp1(upper_lev_idx_zm(k)),  &
                                 LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            wpthlp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wpthlp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like sfc_var.
        if ( k_lowest_zm_input == 2 ) then
          wpthlp(1)  & 
          = lin_ext_zm_bottom( wpthlp(3), wpthlp(2), & 
                               gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if ( input_rtp2 ) then
        call get_grads_var( fread_var, "qtp2", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rtp2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            rtp2(k) = lin_int( gr%zm(k),  &
                               fread_var%z(upper_lev_idx_zm(k)),  &
                               fread_var%z(lower_lev_idx_zm(k)),  &
                               LES_tmp1(upper_lev_idx_zm(k)),  &
                               LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            rtp2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, set the value of rtp2
        ! at momentum level 1 to the value at momentum level 2.
        ! Using a linear extension here resulted in negatives.
        if ( k_lowest_zm_input == 2 ) then
          rtp2(1) =  rtp2(2)
        endif
        if ( any( rtp2(1:gr%nnzp) < rttol**2 ) ) then
          do k=1, gr%nnzp
            rtp2(k) = max(rtp2(k), rttol**2)
          enddo
        endif
      endif

      if ( input_thlp2 ) then
        call get_grads_var( fread_var, "thlp2", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of thlp2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            thlp2(k) = lin_int( gr%zm(k),  &
                                fread_var%z(upper_lev_idx_zm(k)),  &
                                fread_var%z(lower_lev_idx_zm(k)),  &
                                LES_tmp1(upper_lev_idx_zm(k)),  &
                                LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            thlp2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, set the value of thlp2
        ! at momentum level 1 to the value at momentum level 2.
        ! Using a linear extension here resulted in negatives.
        if ( k_lowest_zm_input == 2 ) then
          thlp2(1) = thlp2(2)
        endif
        if ( any( thlp2(1:gr%nnzp) < thltol**2 ) ) then
          do k=1, gr%nnzp
            thlp2(k) = max(thlp2(k), thltol**2)
          enddo
        endif
      endif

      if ( input_rtpthlp) then
        call get_grads_var( fread_var, "qtpthlp", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rtpthlp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            rtpthlp(k) = lin_int( gr%zm(k),  &
                                  fread_var%z(upper_lev_idx_zm(k)),  &
                                  fread_var%z(lower_lev_idx_zm(k)),  &
                                  LES_tmp1(upper_lev_idx_zm(k)),  &
                                  LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            rtpthlp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! rtpthlp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like sfc_var.
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

      if ( input_rcm ) then
        call get_grads_var( fread_var, "qcm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rcm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            rcm(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            rcm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_wm_zt ) then
        call get_grads_var( fread_var, "wlsm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wm_zt from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            wm_zt(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            wm_zt(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_exner ) then
        call get_grads_var( fread_var, "ex0", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of exner from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            exner(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            exner(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_em ) then

        ! Read in SGS TKE
        call get_grads_var( fread_var, "em", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of em from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            em(k) = lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            em(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        ! Read in Resolved TKE
        call get_grads_var( fread_var, "tke", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of em from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            em(k) = em(k) + lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            em(k) = em(k) + LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        where ( em < emin ) em = emin

      endif ! input_em

      if ( input_p ) then
        call get_grads_var( fread_var, "pm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of pressure from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            p_in_Pa(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            p_in_Pa(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_rho ) then
        call get_grads_var( fread_var, "dn0", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rho from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            rho(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            rho(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

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

      if ( input_Kh_zt ) then
        call get_grads_var( fread_var, "kh", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )
        ! LES_tmp1 is the value of mixing length from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            Kh_zt(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            Kh_zt(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

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

      if ( input_wpthvp ) then
        call get_grads_var( fread_var, "wpthvp", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_fatal_error )
        ! LES_tmp1 is the value of wpthvp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            wpthvp(k) = lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            wpthvp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
      endif

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

      if ( input_a ) then
        write(fstderr,*) "The variable a is not setup for input_type" &
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

      if ( input_thvm ) then
        call get_grads_var( fread_var, "thvm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of thvm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            thvm(k) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            thvm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_rrainm ) then
        call get_grads_var( fread_var, "qrm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1

          if ( iirrainm < 1 ) then
            write(fstderr,*) "Rain water mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if

          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iirrainm) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iirrainm) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_Nrm ) then
        call get_grads_var( fread_var, "nrm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1

          if ( iiNrm < 1 ) then
            write(fstderr,*) "Rain droplet number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if

          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iiNrm) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iiNrm) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_Ncm ) then
        call get_grads_var( fread_var, "ncm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1

          if ( iiNcm < 1 ) then
            write(fstderr,*) "Cloud droplet number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if

          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iiNcm) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iiNcm) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif


      if ( input_rsnowm ) then
        call get_grads_var( fread_var, "qsm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( iirsnowm < 1 ) then
            write(fstderr,*) "Snow mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iirsnowm) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iirsnowm) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_ricem ) then
        call get_grads_var( fread_var, "qim", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( iiricem < 1 ) then
            write(fstderr,*) "Ice mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iiricem) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iiricem) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_rgraupelm ) then
        call get_grads_var( fread_var, "qgm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( iirgraupelm < 1 ) then
            write(fstderr,*) "Graupel mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iirgraupelm) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iirgraupelm) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_Ncnm ) then
        call get_grads_var( fread_var, "ncnm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            Ncnm(k) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            Ncnm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_Nim ) then
        call get_grads_var( fread_var, "nim", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( iiNim < 1 ) then
            write(fstderr,*) "Ice number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iiNim) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iiNim) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

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

      if ( input_up2 ) then
        call get_grads_var( fread_var, "up2", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of up2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            up2(k) = lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            up2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        ! Clip up2 to be no smaller than wtol_sqd
        where ( up2 < wtol_sqd ) up2 = wtol_sqd

      endif

      if ( input_vp2 ) then
        call get_grads_var( fread_var, "vp2", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of vp2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            vp2(k) = lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            vp2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        where ( vp2 < wtol_sqd ) vp2 = wtol_sqd

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

      ! stats_sw
      call open_grads_read( 15, stat_file_zm,  fread_var )

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

      if ( input_upwp) then

        call get_grads_var( fread_var, "wpup", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of upwp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            upwp(k) = lin_int( gr%zm(k),  &
                               fread_var%z(upper_lev_idx_zm(k)),  &
                               fread_var%z(lower_lev_idx_zm(k)),  &
                               LES_tmp1(upper_lev_idx_zm(k)),  &
                               LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            upwp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        call get_grads_var( fread_var, "wpup_sgs", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of upwp_sgs from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            upwp(k) = lin_int( gr%zm(k),  &
                               fread_var%z(upper_lev_idx_zm(k)),  &
                               fread_var%z(lower_lev_idx_zm(k)),  &
                               LES_tmp1(upper_lev_idx_zm(k)),  &
                               LES_tmp1(lower_lev_idx_zm(k)) )  &
                      + upwp(k)
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            upwp(k) = LES_tmp1(exact_lev_idx_zm(k)) + upwp(k)
          endif
        enddo

      endif

      if ( input_vpwp) then

        call get_grads_var( fread_var, "wpvp", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of vpwp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            vpwp(k) = lin_int( gr%zm(k),  &
                               fread_var%z(upper_lev_idx_zm(k)),  &
                               fread_var%z(lower_lev_idx_zm(k)),  &
                               LES_tmp1(upper_lev_idx_zm(k)),  &
                               LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            vpwp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        call get_grads_var( fread_var, "wpvp_sgs", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of vpwp_sgs from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            vpwp(k) = lin_int( gr%zm(k),  &
                               fread_var%z(upper_lev_idx_zm(k)),  &
                               fread_var%z(lower_lev_idx_zm(k)),  &
                               LES_tmp1(upper_lev_idx_zm(k)),  &
                               LES_tmp1(lower_lev_idx_zm(k)) )  &
                      + vpwp(k)
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            vpwp(k) = LES_tmp1(exact_lev_idx_zm(k)) + vpwp(k)
          endif
        enddo

      endif

      if ( input_wp2 ) then
        call get_grads_var( fread_var, "wp2", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wp2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            wp2(k) = lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            wp2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        if ( any( wp2(1:gr%nnzp) < wtol_sqd ) ) then
          do k=1, gr%nnzp
            wp2(k) = max( wp2(k), wtol_sqd )
          enddo
        endif
      endif

      if ( input_rho_zm ) then
        call get_grads_var( fread_var, "dn0", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rho_zm from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            rho_zm(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            rho_zm(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
      endif

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

    use inputfile_class, only: &
        open_grads_read,  & ! Procedure(s)
        close_grads_read

    use constants, only:  & 
        sec_per_min ! Variable(s)

    use stats_precision, only:  & 
        time_precision

#ifdef NETCDF
    use input_netcdf, only: open_netcdf_read, close_netcdf_read ! Procedure(s)
#endif

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O unit

    character(len=*), intent(in) ::filename ! Path to the file and its name

    logical, intent(in) :: l_restart ! Whether this is a restart run

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
      call open_grads_read( iunit, trim( filename ), fread_var )
    
    else
#ifdef NETCDF
      call open_netcdf_read( 'thlm', trim( filename ), fread_var, l_error )
#else
      write(fstderr,*) "This version of CLUBB was not compiled with netCDF support"
#endif
    endif

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

      if ( ( mod( delta_time , fread_var%dtwrite )  > 1e-8 ) .or.  & 
           ( mod( delta_time, fread_var%dtwrite ) < -1e-8 ) ) then
        print*, "Error: Elapsed time is not a multiple ", & 
                "of the reference GrADS output time interval."
        print*, "Elapsed time [s] = ", delta_time
        print*, "GrADS output time interval = ", fread_var%dtwrite
        stop
      end if

      if ( mod( delta_time , sec_per_min ) > 1e-8 & 
            .or. mod( delta_time, sec_per_min ) < -1e-8 ) then
        print*, "Error: Elapsed time is not a multiple ", & 
                "of one minute."
        print*, "Elapsed time [s] = ", delta_time
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
               nint( nearest_timestep /  & 
                     (fread_var%dtwrite/sec_per_min) ) - 1
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

    real, intent(in), dimension(vardim) :: &
      clubb_heights ! Altitudes         [m]

    real, intent(inout), dimension(vardim) :: &
      variable_interpolated     ! The variable after interpolation

    logical, intent(out) :: &
      l_read_error ! Whether there was a read error

    logical :: l_spec_bound_cond
    ! Local Variables

    ! ---- Begin Code ----

    if ( l_input_var ) then
      if ( clubb_heights(1) < 0.0 ) then
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
      input_wpthvp, &
      input_thl1, input_thl2, input_a, input_s1, input_s2, &
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
    input_a = .false.
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
