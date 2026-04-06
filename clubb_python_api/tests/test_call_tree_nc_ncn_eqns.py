"""Test wrappers for routines from CLUBB_core/Nc_Ncn_eqns.F90."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_nc_in_cloud_to_ncnm_returns_nc_in_cloud_for_clear_or_constant_case():
    """The fallback branch should return Nc_in_cloud unchanged."""
    out = clubb_api.nc_in_cloud_to_ncnm(
        mu_chi_1=-1.0,
        mu_chi_2=-1.0,
        sigma_chi_1=0.1,
        sigma_chi_2=0.1,
        mixt_frac=0.5,
        nc_in_cloud=42.0,
        cloud_frac_1=0.0,
        cloud_frac_2=0.0,
        const_ncnp2_on_ncnm2=0.0,
        const_corr_chi_ncn=0.0,
    )
    np.testing.assert_allclose(out, 42.0)
