"""User-facing wrappers for simplified radiation helper routines."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.grid_class import Grid
from clubb_python.derived_types.grid_class_converter import set_fortran_grid


def _fa_2d(arr):
    """Convert an array to Fortran-contiguous 2D float64."""
    out = f_arr(arr, dtype=np.float64)
    if out.ndim == 1:
        out = f_arr(out.reshape(1, -1))
    if out.ndim != 2:
        raise ValueError(f"Expected a 1D or 2D array, got rank-{out.ndim}.")
    return out


def _padded_1d(values, fill: float, size: int = 20):
    """Convert values to a fixed-size 1D float64 array with fill."""
    arr = np.asarray(values, dtype=np.float64).reshape(-1)
    out = np.full((size,), fill, dtype=np.float64)
    ncopy = int(min(arr.size, size))
    if ncopy > 0:
        out[:ncopy] = arr[:ncopy]
    return out, ncopy


def set_simplified_radiation_params(
    f0: float,
    f1: float,
    kappa: float,
    eff_drop_radius: float,
    alvdr: float,
    gc: float,
    omega: float,
    l_rad_above_cloud: bool,
    l_sw_radiation: bool,
    l_fix_cos_solar_zen: bool,
    fs_values,
    cos_solar_zen_values,
    cos_solar_zen_times,
):
    """Initialize Fortran simplified-radiation module parameters."""
    fs_values_use, n_fs = _padded_1d(fs_values, fill=0.0)
    cos_values_use, n_cos = _padded_1d(cos_solar_zen_values, fill=-999.0)
    cos_times_use, n_times = _padded_1d(cos_solar_zen_times, fill=-999.0)
    nparam = max(1, n_fs, n_cos, n_times)

    clubb_f2py.f2py_set_simplified_radiation_params(
        float(f0), float(f1), float(kappa),
        float(eff_drop_radius), float(alvdr), float(gc), float(omega),
        bool(l_rad_above_cloud), bool(l_sw_radiation), bool(l_fix_cos_solar_zen),
        int(nparam),
        f_arr(fs_values_use), f_arr(cos_values_use), f_arr(cos_times_use),
    )


def cos_solar_zen(day: int, month: int, year: int, current_time: float,
                  lat_vals: float, lon_vals: float) -> float:
    """Compute cosine of solar zenith angle from date/time and lat/lon."""
    return float(clubb_f2py.f2py_cos_solar_zen(
        int(day), int(month), int(year),
        float(current_time), float(lat_vals), float(lon_vals),
    ))


def sunray_sw(gr: Grid, nzt: int, fs0: float, amu0: float,
              rho, rcm, ngrdcol: int | None = None):
    """Compute shortwave radiative flux Frad_SW on momentum levels.

    Returns Frad_SW with shape (ngrdcol, nzt+1). The caller is responsible
    for computing radht_SW from the flux, matching the pattern in
    radiation_module.F90.
    """
    set_fortran_grid(gr)
    rho_use = _fa_2d(rho)
    rcm_use = _fa_2d(rcm)
    if rho_use.shape != rcm_use.shape:
        raise ValueError(f"rho/rcm shape mismatch: {rho_use.shape} vs {rcm_use.shape}")
    ncol = int(rho_use.shape[0] if ngrdcol is None else ngrdcol)
    return clubb_f2py.f2py_sunray_sw(
        fs0=float(fs0), amu0=float(amu0),
        rho=rho_use, rcm=rcm_use,
        ngrdcol=ncol, nzt=int(nzt),
    )
