"""User-facing wrappers for routines from CLUBB_core/stats_netcdf.F90."""

from __future__ import annotations

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info
from clubb_python.string_conversion import python_strings_to_fortran_char_matrix


def _as_stats_array(values) -> np.ndarray:
    return np.asarray(values, dtype=np.float64)


def _as_level_axis(values) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim == 1:
        return arr
    if arr.ndim >= 2:
        return np.asarray(arr[0, :], dtype=np.float64)
    raise ValueError("Expected a 1D or 2D vertical coordinate array.")


def stats_update(name: str, values):
    """Update a complete surface column vector or a complete profile batch."""
    arr = _as_stats_array(values)
    if arr.ndim == 1:
        clubb_f2py.f2py_stats_update_sfc(name, f_arr(arr), ncol=arr.shape[0])
    elif arr.ndim == 2:
        clubb_f2py.f2py_stats_update_var(name, f_arr(arr), ncol=arr.shape[0], nz=arr.shape[1])
    else:
        raise ValueError("stats_update requires a surface vector or a (col,level) batch.")


def stats_begin_budget(name: str, values):
    """Begin accumulation for a complete (col,level) budget field."""
    arr = _as_stats_array(values)
    if arr.ndim != 2:
        raise ValueError("stats_begin_budget requires a complete (col,level) batch.")
    clubb_f2py.f2py_stats_begin_budget(name, f_arr(arr), ncol=arr.shape[0], nz=arr.shape[1])


def stats_update_budget(name: str, values):
    """Update a complete (col,level) budget field."""
    arr = _as_stats_array(values)
    if arr.ndim != 2:
        raise ValueError("stats_update_budget requires a complete (col,level) batch.")
    clubb_f2py.f2py_stats_update_budget(name, f_arr(arr), ncol=arr.shape[0], nz=arr.shape[1])


def stats_finalize_budget(name: str, values, l_count_sample: bool = True):
    """Finalize a complete (col,level) budget field and optionally count a sample."""
    arr = _as_stats_array(values)
    if arr.ndim != 2:
        raise ValueError("stats_finalize_budget requires a complete (col,level) batch.")
    clubb_f2py.f2py_stats_finalize_budget(
        name, f_arr(arr), l_count_sample, ncol=arr.shape[0], nz=arr.shape[1])


def init_stats(registry_path: str, output_path: str, ncol_batch: int,
               stats_tsamp: float, stats_tout: float, dt_main: float,
               day_in: int, month_in: int, year_in: int,
               time_initial: float,
               zt: np.ndarray, zm: np.ndarray,
               *,
               err_info: ErrInfo,
               clubb_params: np.ndarray | None = None,
               param_names: list[str] | np.ndarray | None = None,
               rad_zt=None, rad_zm=None, hydromet_list=None,
               sclr_dim: int = 0, edsclr_dim: int = 0,
               output_zt=None, output_zm=None, grid_remap_method=None,
               stats_tstart: float | None = None,
               stats_tend: float | None = None):
    """Initialize the stats system with a registry file."""
    if err_info is None or clubb_params is None or param_names is None:
        raise ValueError("init_stats requires err_info, clubb_params, and param_names.")
    zt_levels = _as_level_axis(zt)
    zm_levels = _as_level_axis(zm)
    set_fortran_err_info(err_info)

    missing_time = -1.0e30
    stats_tstart_value = missing_time if stats_tstart is None else float(stats_tstart)
    stats_tend_value = missing_time if stats_tend is None else float(stats_tend)

    if len(param_names) != clubb_params.shape[1]:
        raise ValueError("param_names length must match clubb_params.shape[1].")
    encoded_param_names = python_strings_to_fortran_char_matrix(param_names, width=32)
    # The compiled f2py wrapper currently exposes the legacy initializer
    # surface and does not accept the newer optional remap/radiation args.
    clubb_f2py.f2py_stats_init_with_params(
        registry_path, output_path,
        stats_tsamp, stats_tout, dt_main,
        day_in, month_in, year_in, time_initial,
        stats_tstart=stats_tstart_value,
        stats_tend=stats_tend_value,
        zt=zt_levels,
        zm=zm_levels,
        sclr_dim=int(sclr_dim), edsclr_dim=int(edsclr_dim),
        clubb_params=f_arr(clubb_params), param_names=encoded_param_names,
        ncol=ncol_batch, nzt=zt_levels.shape[0], nzm=zm_levels.shape[0],
    )
    return get_fortran_err_info()


def finalize_stats(err_info: ErrInfo):
    """Finalize the stats system and release resources."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_stats_finalize()
    return get_fortran_err_info()


def get_stats_config():
    """Get stats config scalars."""
    config = list(clubb_f2py.f2py_get_stats_config())
    for logical_idx in (0, 7, 8):
        config[logical_idx] = 1 if config[logical_idx] else 0
    return tuple(config)


def stats_l_sample():
    """Return the current stats sampling flag."""
    return bool(get_stats_config()[7])


def get_stats_var_meta(ivar: int):
    """Get metadata for a stats variable using Python indexing."""
    return clubb_f2py.f2py_get_stats_var_meta(int(ivar) + 1)


def get_stats_var_data(ivar: int, ncol: int, nz: int):
    """Get buffer data for a stats variable using Python indexing."""
    return clubb_f2py.f2py_get_stats_var_data(int(ivar) + 1, ncol, nz)


def set_stats_var_data(ivar: int, ncol: int, nz: int,
                       buffer: np.ndarray, nsamples: np.ndarray):
    """Set buffer data for a stats variable using Python indexing."""
    clubb_f2py.f2py_set_stats_var_data(
        int(ivar) + 1, f_arr(buffer), f_arr(nsamples), ncol=ncol, nz=nz,
    )


def stats_begin_timestep(itime: int):
    """Begin a stats sampling timestep using Python indexing."""
    clubb_f2py.f2py_stats_begin_timestep(int(itime) + 1)


def stats_end_timestep(time_value: float, err_info: ErrInfo):
    """End a stats sampling timestep."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_stats_end_timestep(time_value)
    return get_fortran_err_info()


def var_on_stats_list(name: str) -> bool:
    """Check if a variable name is registered in the stats system."""
    return bool(clubb_f2py.f2py_var_on_stats_list(name))


def stats_update_grid(zt_src: np.ndarray, zm_src: np.ndarray,
                      rho_vals: np.ndarray, rho_levels: np.ndarray,
                      p_sfc: np.ndarray):
    """Update adaptive grid remapping inputs."""
    clubb_f2py.f2py_stats_update_grid(
        f_arr(zt_src), f_arr(zm_src),
        f_arr(rho_vals), f_arr(rho_levels),
        f_arr(p_sfc), ncol=zt_src.shape[0], nzt=zt_src.shape[1], nzm=zm_src.shape[1], nrho=rho_vals.shape[1],
    )


def stats_lh_samples_init(num_samples: int, nzt: int,
                          nl_var_names: list, u_var_names: list,
                          zt_vals: np.ndarray,
                          err_info: ErrInfo):
    """Initialize SILHS sample output variables."""
    if err_info is None:
        raise ValueError("stats_lh_samples_init requires err_info.")
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_stats_lh_samples_init(
        num_samples,
        python_strings_to_fortran_char_matrix(nl_var_names, width=64),
        python_strings_to_fortran_char_matrix(u_var_names, width=64),
        zt_vals,
        nzt=nzt,
        n_nl=len(nl_var_names),
        n_u=len(u_var_names),
    )


def stats_lh_samples_write_lognormal(samples: np.ndarray, err_info: ErrInfo):
    """Write lognormal SILHS sample data."""
    if err_info is None:
        raise ValueError("stats_lh_samples_write_lognormal requires err_info.")
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_stats_lh_samples_write_ln(
        f_arr(samples), ncol=samples.shape[0], nsamp=samples.shape[1], nzt=samples.shape[2], nvars=samples.shape[3])


def stats_lh_samples_write_uniform(uniform_vals: np.ndarray,
                                   mixture_comp: np.ndarray,
                                   sample_weights: np.ndarray,
                                   err_info: ErrInfo):
    """Write uniform SILHS sample data."""
    if err_info is None:
        raise ValueError("stats_lh_samples_write_uniform requires err_info.")
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_stats_lh_samples_write_u(
        f_arr(uniform_vals),
        f_arr(mixture_comp),
        f_arr(sample_weights), ncol=uniform_vals.shape[0], nsamp=uniform_vals.shape[1],
        nzt=uniform_vals.shape[2], nvars=uniform_vals.shape[3],
    )
