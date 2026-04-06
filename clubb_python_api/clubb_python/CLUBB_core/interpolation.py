"""User-facing wrappers for routines from CLUBB_core/interpolation.F90."""

from numpy import asfortranarray as f_arr

import clubb_f2py


def lin_interpolate_two_points(height_int: float, height_high: float,
                               height_low: float, var_high: float,
                               var_low: float) -> float:
    """Linearly interpolate between two points."""
    return float(clubb_f2py.f2py_lin_interpolate_two_points(
        float(height_int), float(height_high), float(height_low),
        float(var_high), float(var_low)))


def mono_cubic_interp(
    z_in: float, km1: int, k00: int, kp1: int, kp2: int,
    zm1: float, z00: float, zp1: float, zp2: float,
    fm1: float, f00: float, fp1: float, fp2: float,
) -> float:
    """Monotone cubic interpolation with Python zero-based level indices."""
    return float(clubb_f2py.f2py_mono_cubic_interp(
        float(z_in), int(km1) + 1, int(k00) + 1, int(kp1) + 1, int(kp2) + 1,
        float(zm1), float(z00), float(zp1), float(zp2),
        float(fm1), float(f00), float(fp1), float(fp2),
    ))


def zlinterp_fnc(dim_out: int, dim_src: int, grid_out, grid_src, var_src):
    """Linearly interpolate a 1D profile onto a target grid."""
    return clubb_f2py.f2py_zlinterp_fnc(
        f_arr(grid_out), f_arr(grid_src), f_arr(var_src), dim_out=int(dim_out), dim_src=int(dim_src)
    )
