! parameters_tunable_wrapper.F90 — wrappers extracted from util_wrappers.F90 for module parameters_tunable

subroutine f2py_check_parameters(ngrdcol, clubb_params, lmin)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: stored_err_info
  use parameters_tunable, only: check_parameters_api

  implicit none

  integer, intent(in) :: ngrdcol
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  real(core_rknd), intent(in) :: lmin

  call check_parameters_api(ngrdcol, clubb_params, lmin, stored_err_info)

end subroutine f2py_check_parameters

subroutine f2py_get_nparams(nparams_out)

  use parameter_indices, only: nparams

  implicit none

  integer, intent(out) :: nparams_out

  nparams_out = nparams

end subroutine f2py_get_nparams

subroutine f2py_init_clubb_params(ngrdcol, iunit, filename, nparams_in, clubb_params)

  use clubb_precision, only: core_rknd
  use parameters_tunable, only: init_clubb_params_api
  use parameter_indices, only: nparams

  implicit none

  integer, intent(in) :: ngrdcol, iunit, nparams_in
  character(len=*), intent(in) :: filename
  real(core_rknd), dimension(ngrdcol, nparams_in), intent(out) :: clubb_params

  if ( nparams_in /= nparams ) then
    error stop "f2py_init_clubb_params: nparams mismatch"
  end if

  call init_clubb_params_api(ngrdcol, iunit, trim(filename), clubb_params)

end subroutine f2py_init_clubb_params

subroutine f2py_get_param_names(nparams_in, param_names)

  use parameters_tunable, only: params_list
  use parameter_indices, only: nparams

  implicit none

  integer, intent(in) :: nparams_in
  character(kind=1), dimension(nparams_in, 28), intent(out) :: param_names
  integer :: i, j

  if ( nparams_in /= nparams ) then
    error stop "f2py_get_param_names: nparams mismatch"
  end if

  param_names = ' '

  do i = 1, nparams_in
    do j = 1, 28
      param_names(i, j) = params_list(i)(j:j)
    end do
  end do

end subroutine f2py_get_param_names

subroutine f2py_init_clubb_params_file(ngrdcol, filename, flen, nparams_in, clubb_params) &
  bind(C, name="f2py_init_clubb_params_file_")

  use clubb_precision, only: core_rknd
  use parameters_tunable, only: init_clubb_params_api
  use parameter_indices, only: nparams

  implicit none

  integer, intent(in) :: ngrdcol
  integer, intent(in) :: flen
  integer, intent(in) :: nparams_in
  character(kind=1), dimension(flen), intent(in) :: filename
  real(core_rknd), dimension(ngrdcol, nparams_in), intent(out) :: clubb_params

  character(len=flen) :: fname_str
  integer :: i

  if ( nparams_in /= nparams ) then
    error stop "f2py_init_clubb_params_file: nparams mismatch"
  end if

  ! Convert C character array to Fortran string
  do i = 1, flen
    fname_str(i:i) = filename(i)
  end do

  call init_clubb_params_api(ngrdcol, 10, fname_str, clubb_params)

end subroutine f2py_init_clubb_params_file

subroutine f2py_calc_derrived_params(ngrdcol, grid_type, deltaz, &
    clubb_params, l_prescribed_avg_deltaz, lmin, mixt_frac_max_mag) &
  bind(C, name="f2py_calc_derrived_params_")

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_grid, stored_nu_vert_res_dep
  use parameters_tunable, only: calc_derrived_params
  use parameter_indices, only: nparams

  implicit none

  integer, intent(in) :: ngrdcol, grid_type
  real(core_rknd), dimension(ngrdcol), intent(in) :: deltaz
  real(core_rknd), dimension(ngrdcol, nparams), intent(in) :: clubb_params
  logical, intent(in) :: l_prescribed_avg_deltaz
  real(core_rknd), intent(out) :: lmin, mixt_frac_max_mag

  lmin = 0.0_core_rknd
  mixt_frac_max_mag = 0.0_core_rknd

  call calc_derrived_params(stored_grid, ngrdcol, grid_type, deltaz, &
                            clubb_params, l_prescribed_avg_deltaz, &
                            stored_nu_vert_res_dep, lmin, mixt_frac_max_mag)

end subroutine f2py_calc_derrived_params
