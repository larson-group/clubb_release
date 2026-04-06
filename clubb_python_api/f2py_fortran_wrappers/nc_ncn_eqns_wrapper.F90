subroutine f2py_nc_in_cloud_to_ncnm(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, mixt_frac, &
    nc_in_cloud, cloud_frac_1, cloud_frac_2, const_ncnp2_on_ncnm2, const_corr_chi_ncn, ncnm)

  use clubb_precision, only: core_rknd
  use nc_ncn_eqns, only: nc_in_cloud_to_ncnm

  implicit none

  real(core_rknd), intent(in) :: &
    mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, mixt_frac, nc_in_cloud, &
    cloud_frac_1, cloud_frac_2, const_ncnp2_on_ncnm2, const_corr_chi_ncn
  real(core_rknd), intent(out) :: ncnm

  ncnm = nc_in_cloud_to_ncnm(mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, mixt_frac, &
    nc_in_cloud, cloud_frac_1, cloud_frac_2, const_ncnp2_on_ncnm2, const_corr_chi_ncn)

end subroutine f2py_nc_in_cloud_to_ncnm
