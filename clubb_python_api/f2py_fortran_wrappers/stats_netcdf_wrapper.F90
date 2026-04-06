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
    nzt, zt, nzm, zm, sclr_dim, edsclr_dim)

  use clubb_precision, only: core_rknd, time_precision
  use derived_type_storage, only: stored_stats, stored_err_info
  use stats_netcdf, only: stats_init_api

  implicit none

  character(len=256), intent(in) :: registry_path, output_path
  integer, intent(in) :: ncol, day_in, month_in, year_in
  real(core_rknd), intent(in) :: stats_tsamp, stats_tout, dt_main
  real(time_precision), intent(in) :: time_initial
  integer, intent(in) :: nzt, nzm, sclr_dim, edsclr_dim
  real(core_rknd), dimension(nzt), intent(in) :: zt
  real(core_rknd), dimension(nzm), intent(in) :: zm

  call stats_init_api( &
    trim(registry_path), trim(output_path), ncol, &
    stats_tsamp, stats_tout, dt_main, &
    day_in, month_in, year_in, time_initial, &
    zt, zm, stored_stats, stored_err_info, &
    sclr_dim=sclr_dim, edsclr_dim=edsclr_dim)

end subroutine f2py_stats_init

subroutine f2py_stats_init_with_params(registry_path, output_path, ncol, &
    stats_tsamp, stats_tout, dt_main, &
    day_in, month_in, year_in, time_initial, &
    nzt, zt, nzm, zm, sclr_dim, edsclr_dim, clubb_params, param_names)

  use clubb_precision, only: core_rknd, time_precision
  use derived_type_storage, only: stored_stats, stored_err_info
  use stats_netcdf, only: stats_init_api
  use parameter_indices, only: nparams

  implicit none

  character(len=256), intent(in) :: registry_path, output_path
  integer, intent(in) :: ncol, day_in, month_in, year_in
  real(core_rknd), intent(in) :: stats_tsamp, stats_tout, dt_main
  real(time_precision), intent(in) :: time_initial
  integer, intent(in) :: nzt, nzm, sclr_dim, edsclr_dim
  real(core_rknd), dimension(nzt), intent(in) :: zt
  real(core_rknd), dimension(nzm), intent(in) :: zm
  real(core_rknd), dimension(ncol, nparams), intent(in) :: clubb_params
  character(len=28), dimension(nparams), intent(in) :: param_names

  call stats_init_api( &
    trim(registry_path), trim(output_path), ncol, &
    stats_tsamp, stats_tout, dt_main, &
    day_in, month_in, year_in, time_initial, &
    zt, zm, stored_stats, stored_err_info, &
    clubb_params=clubb_params, param_names=param_names, &
    sclr_dim=sclr_dim, edsclr_dim=edsclr_dim)

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
  ncol             = stored_stats%ncol
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
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_end_timestep_api

  implicit none

  real(core_rknd), intent(in) :: time_value

  call stats_end_timestep_api(time_value, stored_stats)

end subroutine f2py_stats_end_timestep

!---------------------------------------------------------------------------
! 9. stats_update — scalar variant
!    icol=0 and level=0 mean "not present"
!---------------------------------------------------------------------------
subroutine f2py_stats_update_scalar(name, values, icol, level)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_update

  implicit none

  character(len=64), intent(in) :: name
  real(core_rknd), intent(in) :: values
  integer, intent(in) :: icol, level

  if (level > 0 .and. icol > 0) then
    call stats_update(trim(name), values, stored_stats, icol=icol, level=level)
  else if (icol > 0) then
    call stats_update(trim(name), values, stored_stats, icol=icol)
  end if

end subroutine f2py_stats_update_scalar

!---------------------------------------------------------------------------
! 10. stats_update — 1d variant
!     icol=0 means "not present" (updates all columns for surface)
!---------------------------------------------------------------------------
subroutine f2py_stats_update_1d(name, n, values, icol)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_update

  implicit none

  character(len=64), intent(in) :: name
  integer, intent(in) :: n, icol
  real(core_rknd), dimension(n), intent(in) :: values

  if (icol > 0) then
    call stats_update(trim(name), values, stored_stats, icol=icol)
  else
    call stats_update(trim(name), values, stored_stats)
  end if

end subroutine f2py_stats_update_1d

!---------------------------------------------------------------------------
! 11. stats_update — 2d variant
!---------------------------------------------------------------------------
subroutine f2py_stats_update_2d(name, ncol, nz, values)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_update

  implicit none

  character(len=64), intent(in) :: name
  integer, intent(in) :: ncol, nz
  real(core_rknd), dimension(ncol, nz), intent(in) :: values

  call stats_update(trim(name), values, stored_stats)

end subroutine f2py_stats_update_2d

!---------------------------------------------------------------------------
! 12. stats_begin_budget — scalar
!---------------------------------------------------------------------------
subroutine f2py_stats_begin_budget_scalar(name, values, icol)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_begin_budget

  implicit none

  character(len=64), intent(in) :: name
  real(core_rknd), intent(in) :: values
  integer, intent(in) :: icol

  if (icol > 0) then
    call stats_begin_budget(trim(name), values, stored_stats, icol=icol)
  else
    call stats_begin_budget(trim(name), values, stored_stats)
  end if

end subroutine f2py_stats_begin_budget_scalar

!---------------------------------------------------------------------------
! 13. stats_begin_budget — 1d
!---------------------------------------------------------------------------
subroutine f2py_stats_begin_budget_1d(name, n, values, icol)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_begin_budget

  implicit none

  character(len=64), intent(in) :: name
  integer, intent(in) :: n, icol
  real(core_rknd), dimension(n), intent(in) :: values

  if (icol > 0) then
    call stats_begin_budget(trim(name), values, stored_stats, icol=icol)
  else
    call stats_begin_budget(trim(name), values, stored_stats)
  end if

end subroutine f2py_stats_begin_budget_1d

!---------------------------------------------------------------------------
! 14. stats_begin_budget — 2d
!---------------------------------------------------------------------------
subroutine f2py_stats_begin_budget_2d(name, ncol, nz, values)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_begin_budget

  implicit none

  character(len=64), intent(in) :: name
  integer, intent(in) :: ncol, nz
  real(core_rknd), dimension(ncol, nz), intent(in) :: values

  call stats_begin_budget(trim(name), values, stored_stats)

end subroutine f2py_stats_begin_budget_2d

!---------------------------------------------------------------------------
! 15. stats_update_budget — scalar
!---------------------------------------------------------------------------
subroutine f2py_stats_update_budget_scalar(name, values, icol, level)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_update_budget

  implicit none

  character(len=64), intent(in) :: name
  real(core_rknd), intent(in) :: values
  integer, intent(in) :: icol, level

  if (level > 0 .and. icol > 0) then
    call stats_update_budget(trim(name), values, stored_stats, icol=icol, level=level)
  else if (icol > 0) then
    call stats_update_budget(trim(name), values, stored_stats, icol=icol)
  end if

end subroutine f2py_stats_update_budget_scalar

!---------------------------------------------------------------------------
! 16. stats_update_budget — 1d
!---------------------------------------------------------------------------
subroutine f2py_stats_update_budget_1d(name, n, values, icol)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_update_budget

  implicit none

  character(len=64), intent(in) :: name
  integer, intent(in) :: n, icol
  real(core_rknd), dimension(n), intent(in) :: values

  if (icol > 0) then
    call stats_update_budget(trim(name), values, stored_stats, icol=icol)
  else
    call stats_update_budget(trim(name), values, stored_stats)
  end if

end subroutine f2py_stats_update_budget_1d

!---------------------------------------------------------------------------
! 17. stats_update_budget — 2d
!---------------------------------------------------------------------------
subroutine f2py_stats_update_budget_2d(name, ncol, nz, values)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_update_budget

  implicit none

  character(len=64), intent(in) :: name
  integer, intent(in) :: ncol, nz
  real(core_rknd), dimension(ncol, nz), intent(in) :: values

  call stats_update_budget(trim(name), values, stored_stats)

end subroutine f2py_stats_update_budget_2d

!---------------------------------------------------------------------------
! 18. stats_finalize_budget — scalar
!---------------------------------------------------------------------------
subroutine f2py_stats_finalize_budget_scalar(name, values, icol, l_count_sample)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_finalize_budget

  implicit none

  character(len=64), intent(in) :: name
  real(core_rknd), intent(in) :: values
  integer, intent(in) :: icol
  logical, intent(in) :: l_count_sample

  if (icol > 0) then
    call stats_finalize_budget(trim(name), values, stored_stats, &
      icol=icol, l_count_sample=l_count_sample)
  else
    call stats_finalize_budget(trim(name), values, stored_stats, &
      l_count_sample=l_count_sample)
  end if

end subroutine f2py_stats_finalize_budget_scalar

!---------------------------------------------------------------------------
! 19. stats_finalize_budget — 1d
!---------------------------------------------------------------------------
subroutine f2py_stats_finalize_budget_1d(name, n, values, icol, l_count)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_finalize_budget

  implicit none

  character(len=64), intent(in) :: name
  integer, intent(in) :: n, icol
  logical, intent(in) :: l_count
  real(core_rknd), dimension(n), intent(in) :: values

  if (icol > 0) then
    call stats_finalize_budget(trim(name), values, stored_stats, &
      icol=icol, l_count_sample=l_count)
  else
    call stats_finalize_budget(trim(name), values, stored_stats, &
      l_count_sample=l_count)
  end if

end subroutine f2py_stats_finalize_budget_1d

!---------------------------------------------------------------------------
! 20. stats_finalize_budget — 2d
!---------------------------------------------------------------------------
subroutine f2py_stats_finalize_budget_2d(name, ncol, nz, values, l_count)

  use clubb_precision, only: core_rknd
  use derived_type_storage, only: stored_stats
  use stats_netcdf, only: stats_finalize_budget

  implicit none

  character(len=64), intent(in) :: name
  integer, intent(in) :: ncol, nz
  logical, intent(in) :: l_count
  real(core_rknd), dimension(ncol, nz), intent(in) :: values

  call stats_finalize_budget(trim(name), values, stored_stats, &
    l_count_sample=l_count)

end subroutine f2py_stats_finalize_budget_2d

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

!---------------------------------------------------------------------------
! 22. stats_update_grid — update adaptive grid remapping inputs
!---------------------------------------------------------------------------
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
  character(len=64), dimension(n_nl), intent(in) :: nl_names
  character(len=64), dimension(n_u), intent(in) :: u_names
  real(core_rknd), dimension(nzt), intent(in) :: zt_vals

  call stats_lh_samples_init(num_samples, nzt, nl_names, u_names, &
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
