!-----------------------------------------------------------------------
! $Id$

program clubb_inputfields

! Description:
!   This is a minimalist frontend for the run_clubb subroutine.
!   This version is modified to allow the input of LES, CLUBB or
!   some other pre-calculated input data.
!-----------------------------------------------------------------------
  use clubb_driver, only: run_clubb

  use inputfields, only: input_type, & ! Variable(s)
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
    input_ss1, input_ss2, input_rc1, input_rc2, &
    input_thvm, input_rrainm,input_Nrm,  & 
    input_rsnowm, input_ricem, input_rgraupelm,  & 
    input_thlm_forcing, input_rtm_forcing, & 
    input_up2, input_vp2, input_sigma_sqd_w, input_Ncm,  & 
    input_Ncnm, input_Nim, input_cf, input_sigma_sqd_w_zt, &
    input_veg_T_in_K, input_deep_soil_T_in_K, &
    input_sfc_soil_T_in_K 


  use inputfields, only: set_filenames ! Procedure(s)

  use parameter_indices, only: nparams ! Variable(s)

  use parameters_tunable, only: read_parameters ! Procedure(s)

  use error_code, only: clubb_no_error ! Variable(s)

  use error_code, only: fatal_error ! Procedure(s)

  implicit none

  ! Constant parameters
  integer, parameter :: iunit = 10

  character(len=14), parameter :: finputfields = "inputfields.in"

  ! Run information
  real, dimension(nparams) :: & 
    params  ! Array of the model constants

  character(len=50) :: run_file ! Text file with the namelists

  character(len=80) :: datafile

  logical :: stdout    ! Whether to print iteration number, etc.

  integer :: err_code   ! Numerical diagnostic


  ! Namelist definitions
  namelist /model/ run_file, stdout

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
    input_ss1, input_ss2, input_rc1, input_rc2, &
    input_thvm, input_rrainm,input_Nrm,  & 
    input_rsnowm, input_ricem, input_rgraupelm,  & 
    input_thlm_forcing, input_rtm_forcing, & 
    input_up2, input_vp2, input_sigma_sqd_w, input_Ncm,  & 
    input_Ncnm, input_Nim, input_cf, input_sigma_sqd_w_zt, &
    input_veg_T_in_K, input_deep_soil_T_in_K, &
    input_sfc_soil_T_in_K

!-----------------------------------------------------------------------

  ! --- Begin Code ---

  ! Pick some initial values
  datafile = ''
  input_type = 'hoc'

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
  input_ss1 = .false.
  input_ss2 = .false.
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
  input_cf = .false.
  input_sigma_sqd_w_zt = .false.
  input_veg_T_in_K = .false.
  input_deep_soil_T_in_K = .false.
  input_sfc_soil_T_in_K  = .false.

  ! Read in model constant values
  call read_parameters( iunit, finputfields, params )

  ! Read in our namelists
  open(unit=iunit, file=finputfields, status='old', action='read')

  read(unit=iunit, nml=model)
  read(unit=iunit, nml=setfields)

  close(unit=iunit)

  ! Initialize the status of the run
  err_code = clubb_no_error

  ! Setup the GrADS file reader
  call set_filenames( datafile )

  ! Run the model
  call run_clubb( params, trim( run_file ), err_code, stdout, .true. )

  if ( fatal_error( err_code ) ) then
    stop "Model wasn't valid, check your constants and field input"

  else
    stop "Program exited normally"

  end if

end program clubb_inputfields
!-----------------------------------------------------------------------
