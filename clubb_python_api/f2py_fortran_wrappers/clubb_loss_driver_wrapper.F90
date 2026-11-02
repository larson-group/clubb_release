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

subroutine f2py_get_clubb_params_all(total_param_sets, nparams_in, clubb_params_all_out)

  use clubb_precision, only: core_rknd
  use clubb_driver, only: get_clubb_params_all
  use parameter_indices, only: nparams

  implicit none

  integer, intent(in) :: total_param_sets, nparams_in
  real(core_rknd), dimension(total_param_sets, nparams_in), intent(out) :: clubb_params_all_out

  real(core_rknd), allocatable, dimension(:,:) :: local_params

  if ( nparams_in /= nparams ) then
    error stop "f2py_get_clubb_params_all: nparams mismatch"
  end if

  call get_clubb_params_all(local_params)

  if ( size(local_params, 1) /= total_param_sets ) then
    error stop "f2py_get_clubb_params_all: total_param_sets mismatch"
  end if

  clubb_params_all_out = local_params

end subroutine f2py_get_clubb_params_all

subroutine f2py_clubb_get_loss_for_params(total_param_sets, nparams_in, num_variables, &
                                        num_time_windows, clubb_params_all, &
                                        rmse, correlation, std_ratio, bias_norm, case_loss)

  ! Description:
  !   Expose four profile metrics and one per-candidate case loss with explicit
  !   F2PY output shapes. Request weights are fixed by init_clubb_loss.

  use clubb_precision, only: core_rknd
  use clubb_loss_driver, only: clubb_get_loss_for_params
  use err_info_type_module, only: err_info_type
  use parameter_indices, only: nparams

  implicit none

  !--------------------------- Input Variables ---------------------------
  integer, intent(in) :: total_param_sets, nparams_in ! Parameter-matrix dimensions.
  integer, intent(in) :: num_variables, num_time_windows ! Profile-metric dimensions.
  real(core_rknd), dimension(total_param_sets, nparams_in), intent(in) :: clubb_params_all ! Candidate parameters.

  !--------------------------- Output Variables ---------------------------
  real(core_rknd), dimension(num_time_windows, num_variables, total_param_sets), intent(out) :: rmse ! Profile error.
  real(core_rknd), dimension(num_time_windows, num_variables, total_param_sets), intent(out) :: correlation ! Shape.
  real(core_rknd), dimension(num_time_windows, num_variables, total_param_sets), intent(out) :: std_ratio ! Variability.
  real(core_rknd), dimension(num_time_windows, num_variables, total_param_sets), intent(out) :: bias_norm ! Mean bias.
  real(core_rknd), dimension(total_param_sets), intent(out) :: case_loss ! RMS over fields and windows.

  !--------------------------- Local Variables ---------------------------
  real(core_rknd), allocatable, dimension(:,:) :: local_params ! Source parameter matrix.
  real(core_rknd), allocatable, dimension(:,:,:) :: local_rmse ! Source profile RMSE.
  real(core_rknd), allocatable, dimension(:,:,:) :: local_correlation ! Source correlations.
  real(core_rknd), allocatable, dimension(:,:,:) :: local_std_ratio ! Source std ratios.
  real(core_rknd), allocatable, dimension(:,:,:) :: local_bias_norm ! Source normalized biases.
  real(core_rknd), allocatable, dimension(:) :: local_case_loss ! Source per-candidate losses.
  type(err_info_type) :: err_info ! Column-local runtime error state.

  !--------------------------- Begin Code ---------------------------
  if ( nparams_in /= nparams ) then
    error stop "f2py_clubb_get_loss_for_params: nparams mismatch"
  end if

  allocate( local_params(total_param_sets, nparams_in) )
  local_params = clubb_params_all

  call clubb_get_loss_for_params(local_params, & ! In
                                 local_rmse, local_correlation, local_std_ratio, local_bias_norm, & ! Out
                                 local_case_loss, err_info) ! Out

  if ( any(shape(local_rmse) /= [num_time_windows, num_variables, total_param_sets]) .or. &
       any(shape(local_correlation) /= shape(local_rmse)) .or. &
       any(shape(local_std_ratio) /= shape(local_rmse)) .or. &
       any(shape(local_bias_norm) /= shape(local_rmse)) .or. &
       size(local_case_loss) /= total_param_sets ) then
    error stop "f2py_clubb_get_loss_for_params: profile/case output shape mismatch"
  end if

  rmse = local_rmse
  correlation = local_correlation
  std_ratio = local_std_ratio
  bias_norm = local_bias_norm
  case_loss = local_case_loss

end subroutine f2py_clubb_get_loss_for_params

subroutine f2py_finalize_clubb_loss()

  use clubb_loss_driver, only: finalize_clubb_loss
  use err_info_type_module, only: err_info_type

  implicit none
  type(err_info_type) :: err_info

  call finalize_clubb_loss(err_info)

end subroutine f2py_finalize_clubb_loss
