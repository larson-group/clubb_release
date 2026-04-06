! precipitation_fraction_wrapper.F90 - wrappers for module precipitation_fraction

subroutine f2py_precip_fraction(nzt, ngrdcol, hydromet_dim, hydromet_dim_transport, &
    hydromet, cloud_frac, cloud_frac_1, l_mix_rat_hm, l_frozen_hm, hydromet_tol, &
    cloud_frac_2, ice_supersat_frac, ice_supersat_frac_1, ice_supersat_frac_2, &
    mixt_frac, clubb_params, precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol)

  use clubb_precision, only: core_rknd
  use parameter_indices, only: nparams
  use derived_type_storage, only: stored_grid, stored_stats, stored_err_info
  use precipitation_fraction, only: precip_fraction

  implicit none

  integer, intent(in) :: nzt, ngrdcol, hydromet_dim, hydromet_dim_transport
  real(core_rknd), dimension(ngrdcol, nzt, hydromet_dim_transport), intent(in) :: hydromet
  logical, dimension(hydromet_dim_transport), intent(in) :: l_mix_rat_hm, l_frozen_hm
  real(core_rknd), dimension(hydromet_dim_transport), intent(in) :: hydromet_tol
  real(core_rknd), dimension(ngrdcol, nzt), intent(in) :: &
    cloud_frac, cloud_frac_1, cloud_frac_2, ice_supersat_frac, &
    ice_supersat_frac_1, ice_supersat_frac_2, mixt_frac
  real(core_rknd), dimension(nparams), intent(in) :: clubb_params
  real(core_rknd), dimension(ngrdcol, nzt), intent(out) :: precip_frac, precip_frac_1, precip_frac_2
  real(core_rknd), dimension(ngrdcol), intent(out) :: precip_frac_tol

  call precip_fraction( &
    stored_grid, nzt, ngrdcol, hydromet_dim, &
    hydromet(:, :, 1:hydromet_dim), cloud_frac, cloud_frac_1, &
    l_mix_rat_hm(1:hydromet_dim), l_frozen_hm(1:hydromet_dim), hydromet_tol(1:hydromet_dim), &
    cloud_frac_2, ice_supersat_frac, ice_supersat_frac_1, ice_supersat_frac_2, &
    mixt_frac, clubb_params, stored_err_info, &
    precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol, stored_stats)

end subroutine f2py_precip_fraction
