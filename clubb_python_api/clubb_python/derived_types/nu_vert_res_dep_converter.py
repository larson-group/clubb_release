"""Push/pull converters for NuVertResDep <-> Fortran module storage."""

import numpy as np

import clubb_f2py

from clubb_python.derived_types.nu_vert_res_dep import NuVertResDep


def set_fortran_nu_vert_res_dep(nu: NuVertResDep):
    """Push nu coefficients into Fortran module storage."""
    clubb_f2py.set_nu_vert_res_dep(
        np.asfortranarray(nu.nu1),
        np.asfortranarray(nu.nu2),
        np.asfortranarray(nu.nu6),
        np.asfortranarray(nu.nu8),
        np.asfortranarray(nu.nu9),
        np.asfortranarray(nu.nu10),
        np.asfortranarray(nu.nu_hm),
    )


def get_fortran_nu_vert_res_dep(nzm: int) -> NuVertResDep:
    """Pull nu coefficients from Fortran module storage."""
    nu1, nu2, nu6, nu8, nu9, nu10, nu_hm = clubb_f2py.get_nu_vert_res_dep(nzm)
    return NuVertResDep(
        nzm=nzm,
        nu1=nu1,
        nu2=nu2,
        nu6=nu6,
        nu8=nu8,
        nu9=nu9,
        nu10=nu10,
        nu_hm=nu_hm,
    )
