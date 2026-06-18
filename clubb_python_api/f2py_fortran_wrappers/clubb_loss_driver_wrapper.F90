! clubb_loss_driver_wrapper.F90 -- F2PY wrappers for reusable CLUBB loss routines

subroutine f2py_init_clubb_loss(runfile, num_variables, total_param_sets, num_time_windows, clubb_var_names)

  use clubb_loss_driver, only: init_clubb_loss, get_loss_time_window_count, loss_name_len, max_loss_variables
  use clubb_driver, only: get_runtime_batch_config
  use err_info_type_module, only: err_info_type

  implicit none

  character(len=*), intent(in) :: runfile
  integer, intent(out) :: num_variables, total_param_sets, num_time_windows
  character(kind=1), dimension(max_loss_variables, loss_name_len), intent(out) :: clubb_var_names
  character(len=loss_name_len), allocatable, dimension(:) :: local_names
  type(err_info_type) :: err_info
  integer :: batch_size
  integer :: i, j

  call init_clubb_loss(trim(runfile), local_names, err_info)
  call get_runtime_batch_config(total_param_sets, batch_size)
  call get_loss_time_window_count(num_time_windows)

  num_variables = size(local_names)
  clubb_var_names = ' '
  do i = 1, num_variables
    do j = 1, loss_name_len
      clubb_var_names(i, j) = local_names(i)(j:j)
    end do
  end do

end subroutine f2py_init_clubb_loss

subroutine f2py_get_clubb_params_all(total_param_sets, nparams_in, clubb_params)

  use clubb_precision, only: core_rknd
  use clubb_driver, only: get_clubb_params_all
  use parameter_indices, only: nparams

  implicit none

  integer, intent(in) :: total_param_sets, nparams_in
  real(core_rknd), dimension(total_param_sets, nparams_in), intent(out) :: clubb_params

  real(core_rknd), allocatable, dimension(:,:) :: local_params

  if ( nparams_in /= nparams ) then
    error stop "f2py_get_clubb_params_all: nparams mismatch"
  end if

  call get_clubb_params_all(local_params)

  if ( size(local_params, 1) /= total_param_sets ) then
    error stop "f2py_get_clubb_params_all: total_param_sets mismatch"
  end if

  clubb_params = local_params

end subroutine f2py_get_clubb_params_all

subroutine f2py_clubb_get_loss_for_params(total_param_sets, nparams_in, num_variables, &
                                          num_time_windows, clubb_params_all, scaled_rmse, correlation, &
                                          std_ratio, centered_rmse_norm, bias_norm)

  use clubb_precision, only: core_rknd
  use clubb_loss_driver, only: clubb_get_loss_for_params
  use err_info_type_module, only: err_info_type
  use parameter_indices, only: nparams

  implicit none

  integer, intent(in) :: total_param_sets, nparams_in, num_variables, num_time_windows
  real(core_rknd), dimension(total_param_sets, nparams_in), intent(in) :: clubb_params_all
  real(core_rknd), dimension(num_time_windows, num_variables, total_param_sets), intent(out) :: scaled_rmse
  real(core_rknd), dimension(num_time_windows, num_variables, total_param_sets), intent(out) :: correlation
  real(core_rknd), dimension(num_time_windows, num_variables, total_param_sets), intent(out) :: std_ratio
  real(core_rknd), dimension(num_time_windows, num_variables, total_param_sets), intent(out) :: centered_rmse_norm
  real(core_rknd), dimension(num_time_windows, num_variables, total_param_sets), intent(out) :: bias_norm

  real(core_rknd), allocatable, dimension(:,:) :: local_params
  real(core_rknd), allocatable, dimension(:,:,:) :: local_scaled_rmse
  real(core_rknd), allocatable, dimension(:,:,:) :: local_correlation
  real(core_rknd), allocatable, dimension(:,:,:) :: local_std_ratio
  real(core_rknd), allocatable, dimension(:,:,:) :: local_centered_rmse_norm
  real(core_rknd), allocatable, dimension(:,:,:) :: local_bias_norm
  type(err_info_type) :: err_info

  if ( nparams_in /= nparams ) then
    error stop "f2py_clubb_get_loss_for_params: nparams mismatch"
  end if

  allocate( local_params(total_param_sets, nparams_in) )
  local_params = clubb_params_all

  call clubb_get_loss_for_params(local_params, local_scaled_rmse, local_correlation, local_std_ratio, &
                                 local_centered_rmse_norm, local_bias_norm, err_info)

  if ( size(local_scaled_rmse, 1) /= num_time_windows .or. &
       size(local_scaled_rmse, 2) /= num_variables .or. &
       size(local_scaled_rmse, 3) /= total_param_sets ) then
    error stop "f2py_clubb_get_loss_for_params: scaled_rmse shape mismatch"
  end if

  if ( size(local_correlation, 1) /= num_time_windows .or. &
       size(local_correlation, 2) /= num_variables .or. &
       size(local_correlation, 3) /= total_param_sets ) then
    error stop "f2py_clubb_get_loss_for_params: correlation shape mismatch"
  end if

  if ( size(local_std_ratio, 1) /= num_time_windows .or. &
       size(local_std_ratio, 2) /= num_variables .or. &
       size(local_std_ratio, 3) /= total_param_sets ) then
    error stop "f2py_clubb_get_loss_for_params: std_ratio shape mismatch"
  end if

  if ( size(local_centered_rmse_norm, 1) /= num_time_windows .or. &
       size(local_centered_rmse_norm, 2) /= num_variables .or. &
       size(local_centered_rmse_norm, 3) /= total_param_sets ) then
    error stop "f2py_clubb_get_loss_for_params: centered_rmse_norm shape mismatch"
  end if

  if ( size(local_bias_norm, 1) /= num_time_windows .or. &
       size(local_bias_norm, 2) /= num_variables .or. &
       size(local_bias_norm, 3) /= total_param_sets ) then
    error stop "f2py_clubb_get_loss_for_params: bias_norm shape mismatch"
  end if

  scaled_rmse = local_scaled_rmse
  correlation = local_correlation
  std_ratio = local_std_ratio
  centered_rmse_norm = local_centered_rmse_norm
  bias_norm = local_bias_norm

end subroutine f2py_clubb_get_loss_for_params

subroutine f2py_finalize_clubb_loss()

  use clubb_loss_driver, only: finalize_clubb_loss
  use err_info_type_module, only: err_info_type

  implicit none
  type(err_info_type) :: err_info

  call finalize_clubb_loss(err_info)

end subroutine f2py_finalize_clubb_loss
