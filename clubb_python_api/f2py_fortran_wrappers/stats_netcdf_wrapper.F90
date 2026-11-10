! stats_netcdf_wrapper.F90 — F2PY wrappers for CLUBB stats system
!
! Standalone subroutines (NOT in a module) so they get standard Fortran
! name mangling (f2py_stats_init_) which F2PY expects. This also avoids
! the bind(C) restriction on character(len>1) arguments.
!
! All subroutines operate on stored_stats from derived_type_storage.
! Logicals are passed as integers (0/1). Optional args use sentinel values.

!---------------------------------------------------------------------------
! 1. Initialize stats with registry
!---------------------------------------------------------------------------
subroutine f2py_stats_init(registry_path, output_path, ncol, &
    stats_tsamp, stats_tout, dt_main, &
    day_in, month_in, year_in, time_initial, &
    stats_tstart, stats_tend, nzt, zt, nzm, zm, sclr_dim, edsclr_dim)

  use clubb_precision, only: core_rknd, time_precision
  use derived_type_storage, only: stored_stats, stored_err_info
  use stats_netcdf, only: stats_init_api

  implicit none

  character(len=256), intent(in) :: registry_path, output_path
  integer, intent(in) :: ncol, day_in, month_in, year_in
  real(core_rknd), intent(in) :: stats_tsamp, stats_tout, dt_main
  real(time_precision), intent(in) :: time_initial
  real(time_precision), intent(in) :: stats_tstart, stats_tend
  integer, intent(in) :: nzt, nzm, sclr_dim, edsclr_dim
  real(core_rknd), dimension(nzt), intent(in) :: zt
  real(core_rknd), dimension(nzm), intent(in) :: zm
  real(time_precision), parameter :: missing_time = -1.0e30_time_precision
  logical :: has_stats_tstart, has_stats_tend

  has_stats_tstart = stats_tstart > 0.5_time_precision * missing_time
  has_stats_tend = stats_tend > 0.5_time_precision * missing_time

  if ( has_stats_tstart .and. has_stats_tend ) then
    call stats_init_api( &
      trim(registry_path), trim(output_path), ncol, &
      stats_tsamp, stats_tout, dt_main, &
      day_in, month_in, year_in, time_initial, zt, zm, stored_stats, stored_err_info, &
      sclr_dim=sclr_dim, edsclr_dim=edsclr_dim, &
      stats_tstart=stats_tstart, stats_tend=stats_tend)
  else if ( has_stats_tstart ) then
    call stats_init_api( &
      trim(registry_path), trim(output_path), ncol, &
      stats_tsamp, stats_tout, dt_main, &
      day_in, month_in, year_in, time_initial, zt, zm, stored_stats, stored_err_info, &
      sclr_dim=sclr_dim, edsclr_dim=edsclr_dim, &
      stats_tstart=stats_tstart)
  else if ( has_stats_tend ) then
    call stats_init_api( &
      trim(registry_path), trim(output_path), ncol, &
      stats_tsamp, stats_tout, dt_main, &
      day_in, month_in, year_in, time_initial, zt, zm, stored_stats, stored_err_info, &
      sclr_dim=sclr_dim, edsclr_dim=edsclr_dim, &
      stats_tend=stats_tend)
  else
    call stats_init_api( &
      trim(registry_path), trim(output_path), ncol, &
      stats_tsamp, stats_tout, dt_main, &
      day_in, month_in, year_in, time_initial, zt, zm, stored_stats, stored_err_info, &
      sclr_dim=sclr_dim, edsclr_dim=edsclr_dim)
  end if

end subroutine f2py_stats_init

subroutine f2py_stats_init_with_params(registry_path, output_path, ncol, &
    stats_tsamp, stats_tout, dt_main, &
    day_in, month_in, year_in, time_initial, &
    stats_tstart, stats_tend, nzt, zt, nzm, zm, sclr_dim, edsclr_dim, &
    clubb_params, param_names)

  use clubb_precision, only: core_rknd, time_precision
  use derived_type_storage, only: stored_stats, stored_err_info
  use stats_netcdf, only: stats_init_api
  use parameter_indices, only: nparams

  implicit none

  character(len=256), intent(in) :: registry_path, output_path
  integer, intent(in) :: ncol, day_in, month_in, year_in
  real(core_rknd), intent(in) :: stats_tsamp, stats_tout, dt_main
  real(time_precision), intent(in) :: time_initial
  real(time_precision), intent(in) :: stats_tstart, stats_tend
  integer, intent(in) :: nzt, nzm, sclr_dim, edsclr_dim
  real(core_rknd), dimension(nzt), intent(in) :: zt
  real(core_rknd), dimension(nzm), intent(in) :: zm
  real(core_rknd), dimension(ncol, nparams), intent(in) :: clubb_params
  character(kind=1), dimension(nparams, 32), intent(in) :: param_names
  character(len=32), dimension(nparams) :: param_names_fortran
  real(time_precision), parameter :: missing_time = -1.0e30_time_precision
  integer :: i, j
  logical :: has_stats_tstart, has_stats_tend

  do i = 1, nparams
    param_names_fortran(i) = ' '
    do j = 1, 32
      param_names_fortran(i)(j:j) = param_names(i, j)
    end do
  end do

  has_stats_tstart = stats_tstart > 0.5_time_precision * missing_time
  has_stats_tend = stats_tend > 0.5_time_precision * missing_time

  if ( has_stats_tstart .and. has_stats_tend ) then
    call stats_init_api( &
      trim(registry_path), trim(output_path), ncol, &
      stats_tsamp, stats_tout, dt_main, &
      day_in, month_in, year_in, time_initial, zt, zm, stored_stats, stored_err_info, &
      clubb_params=clubb_params, param_names=param_names_fortran, &
      sclr_dim=sclr_dim, edsclr_dim=edsclr_dim, &
      stats_tstart=stats_tstart, stats_tend=stats_tend)
  else if ( has_stats_tstart ) then
    call stats_init_api( &
      trim(registry_path), trim(output_path), ncol, &
      stats_tsamp, stats_tout, dt_main, &
      day_in, month_in, year_in, time_initial, zt, zm, stored_stats, stored_err_info, &
      clubb_params=clubb_params, param_names=param_names_fortran, &
      sclr_dim=sclr_dim, edsclr_dim=edsclr_dim, &
      stats_tstart=stats_tstart)
  else if ( has_stats_tend ) then
    call stats_init_api( &
      trim(registry_path), trim(output_path), ncol, &
      stats_tsamp, stats_tout, dt_main, &
      day_in, month_in, year_in, time_initial, zt, zm, stored_stats, stored_err_info, &
      clubb_params=clubb_params, param_names=param_names_fortran, &
      sclr_dim=sclr_dim, edsclr_dim=edsclr_dim, &
      stats_tend=stats_tend)
  else
    call stats_init_api( &
      trim(registry_path), trim(output_path), ncol, &
      stats_tsamp, stats_tout, dt_main, &
      day_in, month_in, year_in, time_initial, zt, zm, stored_stats, stored_err_info, &
      clubb_params=clubb_params, param_names=param_names_fortran, &
      sclr_dim=sclr_dim, edsclr_dim=edsclr_dim)
  end if

end subroutine f2py_stats_init_with_params

!---------------------------------------------------------------------------
! 2. Finalize stats
!---------------------------------------------------------------------------
subroutine f2py_stats_finalize()

  use derived_type_storage, only: stored_stats, stored_err_info
  use stats_netcdf, only: stats_finalize_api

  implicit none

  call stats_finalize_api(stored_stats, stored_err_info)

end subroutine f2py_stats_finalize

!---------------------------------------------------------------------------
! 3. Get stats configuration scalars
!---------------------------------------------------------------------------
subroutine f2py_get_stats_config(l_enabled, ncol, nvars, &
    stats_nsamp, stats_nout, samples_per_write, &
    time_index, l_sample, l_last_sample)

  use derived_type_storage, only: stored_stats

  implicit none

  integer, intent(out) :: ncol, nvars
  logical, intent(out) :: l_enabled
  integer, intent(out) :: stats_nsamp, stats_nout, samples_per_write
  integer, intent(out) :: time_index
  logical, intent(out) :: l_sample, l_last_sample

  l_enabled       = stored_stats%enabled
  ncol             = stored_stats%ncol_batch
  nvars            = stored_stats%nvars
  stats_nsamp      = stored_stats%stats_nsamp
  stats_nout       = stored_stats%stats_nout
  samples_per_write = stored_stats%samples_per_write
  time_index       = stored_stats%time_index
  l_sample        = stored_stats%l_sample
  l_last_sample   = stored_stats%l_last_sample

end subroutine f2py_get_stats_config

!---------------------------------------------------------------------------
! 4. Get metadata for one stats variable (1-based index)
!---------------------------------------------------------------------------
subroutine f2py_get_stats_var_meta(ivar, name, grid_name, units, &
    long_name, grid_id, nz)

  use derived_type_storage, only: stored_stats

  implicit none

  integer, intent(in) :: ivar
  character(len=64), intent(out) :: name
  character(len=16), intent(out) :: grid_name
  character(len=32), intent(out) :: units
  character(len=128), intent(out) :: long_name
  integer, intent(out) :: grid_id, nz

  name      = stored_stats%vars(ivar)%name
  grid_name = stored_stats%vars(ivar)%grid
  units     = stored_stats%vars(ivar)%units
  long_name = stored_stats%vars(ivar)%long_name
  grid_id   = stored_stats%vars(ivar)%grid_id
  nz        = stored_stats%vars(ivar)%nz

end subroutine f2py_get_stats_var_meta

!---------------------------------------------------------------------------
! 5. Get buffer data for one stats variable
!---------------------------------------------------------------------------
subroutine f2py_get_stats_var_data(ivar, ncol, nz, buffer, nsamples)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats

  implicit none

  integer, intent(in) :: ivar, ncol, nz
  real(core_rknd), dimension(ncol, nz), intent(out) :: buffer
  integer, dimension(ncol, nz), intent(out) :: nsamples

  buffer   = real(stored_stats%vars(ivar)%buffer, core_rknd)
  nsamples = stored_stats%vars(ivar)%nsamples

end subroutine f2py_get_stats_var_data

!---------------------------------------------------------------------------
! 6. Set buffer data for one stats variable
!---------------------------------------------------------------------------
subroutine f2py_set_stats_var_data(ivar, ncol, nz, buffer, nsamples)

  use clubb_precision, only: core_rknd, stat_rknd
  use derived_type_storage, only: stored_stats

  implicit none

  integer, intent(in) :: ivar, ncol, nz
  real(core_rknd), dimension(ncol, nz), intent(in) :: buffer
  integer, dimension(ncol, nz), intent(in) :: nsamples

  stored_stats%vars(ivar)%buffer   = real(buffer, stat_rknd)
  stored_stats%vars(ivar)%nsamples = nsamples

end subroutine f2py_set_stats_var_data

!---------------------------------------------------------------------------
! 7. Begin timestep
!---------------------------------------------------------------------------
subroutine f2py_stats_begin_timestep(itime)

  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_begin_timestep_api

  implicit none

  integer, intent(in) :: itime

  call stats_begin_timestep_api(itime, stored_stats)

end subroutine f2py_stats_begin_timestep

!---------------------------------------------------------------------------
! 8. End timestep
!---------------------------------------------------------------------------
subroutine f2py_stats_end_timestep(time_value)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats, stored_err_info
  use stats_netcdf, only: stats_end_timestep_api

  implicit none

  real(core_rknd), intent(in) :: time_value

  call stats_end_timestep_api(time_value, stored_stats, stored_err_info)

end subroutine f2py_stats_end_timestep

subroutine f2py_stats_update_sfc( name, ncol, values ) ! In

  ! Description:
  !   Accumulate a complete surface column vector.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_update

  implicit none

  !--------------------------- Input Variables ---------------------------
  character(len=64), intent(in) :: name ! Registered variable name to update [-]
  integer, intent(in) :: ncol ! Number of grid columns [-]
  ! Surface values indexed by column [var-dependent]
  real(core_rknd), dimension(ncol), intent(in) :: values

  !--------------------------- Begin Code ---------------------------

  call stats_update(trim(name), values, & ! In
    stored_stats) ! InOut

end subroutine f2py_stats_update_sfc

subroutine f2py_stats_update_var( name, ncol, nz, values ) ! In

  ! Description:
  !   Accumulate a complete profile column/level batch.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_update

  implicit none

  !--------------------------- Input Variables ---------------------------
  character(len=64), intent(in) :: name ! Registered variable name to update [-]
  integer, intent(in) :: ncol ! Number of grid columns [-]
  integer, intent(in) :: nz ! Number of vertical levels [-]
  ! 2D payload (col,level) in variable units [var-dependent]
  real(core_rknd), dimension(ncol, nz), intent(in) :: values

  !--------------------------- Begin Code ---------------------------

  call stats_update(trim(name), values, & ! In
    stored_stats) ! InOut

end subroutine f2py_stats_update_var

subroutine f2py_stats_begin_budget( name, ncol, nz, values ) ! In

  ! Description:
  !   Call stats_begin_budget with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_begin_budget

  implicit none

  !--------------------------- Input Variables ---------------------------
  character(len=64), intent(in) :: name ! Registered budget variable name [-]
  integer, intent(in) :: ncol ! Number of grid columns [-]
  integer, intent(in) :: nz ! Number of vertical levels [-]
  ! 2D state payload (col,level) in variable units [var-dependent]
  real(core_rknd), dimension(ncol, nz), intent(in) :: values

  !--------------------------- Begin Code ---------------------------

  call stats_begin_budget(trim(name), values, stored_stats)

end subroutine f2py_stats_begin_budget

subroutine f2py_stats_update_budget( name, ncol, nz, values ) ! In

  ! Description:
  !   Accumulate a complete profile column/level batch.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_update_budget

  implicit none

  !--------------------------- Input Variables ---------------------------
  character(len=64), intent(in) :: name ! Registered budget variable name [-]
  integer, intent(in) :: ncol ! Number of grid columns [-]
  integer, intent(in) :: nz ! Number of vertical levels [-]
  ! 2D incremental budget contribution [var-dependent]
  real(core_rknd), dimension(ncol, nz), intent(in) :: values

  !--------------------------- Begin Code ---------------------------

  call stats_update_budget(trim(name), values, & ! In
    stored_stats) ! InOut

end subroutine f2py_stats_update_budget

subroutine f2py_stats_finalize_budget( name, ncol, nz, values, l_count ) ! In

  ! Description:
  !   Call stats_finalize_budget with complete column batches through the F2PY bridge.

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_finalize_budget

  implicit none

  !--------------------------- Input Variables ---------------------------
  character(len=64), intent(in) :: name ! Registered budget variable name [-]
  integer, intent(in) :: ncol ! Number of grid columns [-]
  integer, intent(in) :: nz ! Number of vertical levels [-]
  ! 2D closing state payload in variable units [var-dependent]
  real(core_rknd), dimension(ncol, nz), intent(in) :: values
  logical, intent(in) :: l_count ! l count [var-dependent]

  !--------------------------- Begin Code ---------------------------

  call stats_finalize_budget( trim(name), values, & ! In
                              stored_stats, & ! InOut
                              l_count_sample=l_count ) ! Optional In

end subroutine f2py_stats_finalize_budget

!---------------------------------------------------------------------------
! 21. var_on_stats_list — query if variable is registered
!---------------------------------------------------------------------------
subroutine f2py_var_on_stats_list(name, l_result)

  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: var_on_stats_list

  implicit none

  character(len=64), intent(in) :: name
  logical, intent(out) :: l_result

  l_result = var_on_stats_list(stored_stats, trim(name))

end subroutine f2py_var_on_stats_list

subroutine f2py_stats_update_grid(ncol, nzt, nzm, nrho, &
    zt_src, zm_src, rho_vals, rho_levels, p_sfc)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_update_grid

  implicit none

  integer, intent(in) :: ncol, nzt, nzm, nrho
  real(core_rknd), dimension(ncol, nzt), intent(in) :: zt_src
  real(core_rknd), dimension(ncol, nzm), intent(in) :: zm_src
  real(core_rknd), dimension(ncol, nrho), intent(in) :: rho_vals, rho_levels
  real(core_rknd), dimension(ncol), intent(in) :: p_sfc

  call stats_update_grid(zt_src, zm_src, rho_vals, rho_levels, p_sfc, &
                         stored_stats)

end subroutine f2py_stats_update_grid

!---------------------------------------------------------------------------
! 23. stats_lh_samples_init — define SILHS sample output variables
!---------------------------------------------------------------------------
subroutine f2py_stats_lh_samples_init(num_samples, nzt, &
    n_nl, n_u, nl_names, u_names, zt_vals)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats, stored_err_info
  use stats_netcdf, only: stats_lh_samples_init

  implicit none

  integer, intent(in) :: num_samples, nzt, n_nl, n_u
  character(kind=1), dimension(n_nl, 64), intent(in) :: nl_names
  character(kind=1), dimension(n_u, 64), intent(in) :: u_names
  real(core_rknd), dimension(nzt), intent(in) :: zt_vals
  character(len=64), dimension(n_nl) :: nl_names_fortran
  character(len=64), dimension(n_u) :: u_names_fortran
  integer :: i, j

  do i = 1, n_nl
    nl_names_fortran(i) = ' '
    do j = 1, 64
      nl_names_fortran(i)(j:j) = nl_names(i, j)
    end do
  end do

  do i = 1, n_u
    u_names_fortran(i) = ' '
    do j = 1, 64
      u_names_fortran(i)(j:j) = u_names(i, j)
    end do
  end do

  call stats_lh_samples_init(num_samples, nzt, nl_names_fortran, u_names_fortran, &
                             zt_vals, stored_stats, stored_err_info)

end subroutine f2py_stats_lh_samples_init

!---------------------------------------------------------------------------
! 24. stats_lh_samples_write_lognormal
!---------------------------------------------------------------------------
subroutine f2py_stats_lh_samples_write_ln(ncol, nsamp, nzt, nvars, samples)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats, stored_err_info
  use stats_netcdf, only: stats_lh_samples_write_lognormal

  implicit none

  integer, intent(in) :: ncol, nsamp, nzt, nvars
  real(core_rknd), dimension(ncol, nsamp, nzt, nvars), intent(in) :: samples

  call stats_lh_samples_write_lognormal(samples, stored_stats, stored_err_info)

end subroutine f2py_stats_lh_samples_write_ln

!---------------------------------------------------------------------------
! 25. stats_lh_samples_write_uniform
!---------------------------------------------------------------------------
subroutine f2py_stats_lh_samples_write_u(ncol, nsamp, nzt, nvars, &
    uniform_vals, mixture_comp, sample_weights)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats, stored_err_info
  use stats_netcdf, only: stats_lh_samples_write_uniform

  implicit none

  integer, intent(in) :: ncol, nsamp, nzt, nvars
  real(core_rknd), dimension(ncol, nsamp, nzt, nvars), intent(in) :: uniform_vals
  integer, dimension(ncol, nsamp, nzt), intent(in) :: mixture_comp
  real(core_rknd), dimension(ncol, nsamp, nzt), intent(in) :: sample_weights

  call stats_lh_samples_write_uniform(uniform_vals, mixture_comp, &
    sample_weights, stored_stats, stored_err_info)

end subroutine f2py_stats_lh_samples_write_u
