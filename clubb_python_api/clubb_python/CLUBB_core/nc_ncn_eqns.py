"""User-facing wrappers for routines from CLUBB_core/Nc_Ncn_eqns.F90."""

import clubb_f2py


def nc_in_cloud_to_ncnm(
    mu_chi_1: float, mu_chi_2: float, sigma_chi_1: float, sigma_chi_2: float, mixt_frac: float,
    nc_in_cloud: float, cloud_frac_1: float, cloud_frac_2: float,
    const_ncnp2_on_ncnm2: float, const_corr_chi_ncn: float,
) -> float:
    """Convert in-cloud droplet concentration to mean simplified cloud nuclei concentration."""
    return float(
        clubb_f2py.f2py_nc_in_cloud_to_ncnm(
            float(mu_chi_1), float(mu_chi_2), float(sigma_chi_1), float(sigma_chi_2), float(mixt_frac),
            float(nc_in_cloud), float(cloud_frac_1), float(cloud_frac_2),
            float(const_ncnp2_on_ncnm2), float(const_corr_chi_ncn),
        )
    )
