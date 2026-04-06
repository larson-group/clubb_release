"""User-facing wrappers for routines from CLUBB_core/fill_holes.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py


def fill_holes_vertical(
    nz: int, ngrdcol: int, threshold: float,
    lower_hf_level: int, upper_hf_level: int,
    dz, rho_ds, grid_dir_indx: int, fill_holes_type: int, field,
):
    """Fill negative holes in a vertical field using Python level indices."""
    return clubb_f2py.f2py_fill_holes_vertical(
        threshold, int(lower_hf_level) + 1, int(upper_hf_level) + 1,
        f_arr(dz), f_arr(rho_ds), grid_dir_indx, fill_holes_type,
        f_arr(field), nz=int(nz), ngrdcol=int(ngrdcol))


def fill_holes_wp2_from_horz_tke(
    nz: int, ngrdcol: int, threshold: float,
    lower_hf_level: int, upper_hf_level: int,
    wp2, up2, vp2,
):
    """Fill wp2 holes by borrowing TKE from up2 and vp2 using Python level indices."""
    return clubb_f2py.f2py_fill_holes_wp2_from_horz_tke(
        threshold, int(lower_hf_level) + 1, int(upper_hf_level) + 1,
        f_arr(wp2), f_arr(up2), f_arr(vp2), nz=int(nz), ngrdcol=int(ngrdcol))
