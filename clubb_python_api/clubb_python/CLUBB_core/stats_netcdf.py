"""User-facing wrappers for routines from CLUBB_core/stats_netcdf.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.err_info import ErrInfo
from clubb_python.derived_types.err_info_converter import get_fortran_err_info, set_fortran_err_info


def _as_stats_array(values) -> np.ndarray:
    return np.asarray(values, dtype=np.float64)


def stats_update(name: str, values, icol: int | None = None,
                 level: int | None = None):
    """Update a stats variable using scalar, 1D, or 2D payloads."""
    arr = _as_stats_array(values)
    icol_idx = 0 if icol is None else int(icol) + 1
    level_idx = 0 if level is None else int(level) + 1

    if arr.ndim == 0:
        if icol_idx <= 0:
            raise ValueError("Scalar stats_update requires icol > 0.")
        clubb_f2py.f2py_stats_update_scalar(name, float(arr), icol_idx, level_idx)
        return

    if arr.ndim == 1:
        if icol_idx > 0 and level_idx > 0:
            raise ValueError("rank-1 stats_update accepts icol or level, not both.")
        if level_idx > 0:
            raise ValueError("rank-1 stats_update does not support level-specific dispatch.")
        clubb_f2py.f2py_stats_update_1d(name, arr, icol_idx, n=arr.shape[0])
        return

    if arr.ndim == 2:
        if icol_idx > 0 or level_idx > 0:
            raise ValueError("rank-2 stats_update does not support icol or level.")
        clubb_f2py.f2py_stats_update_2d(name, f_arr(arr), ncol=arr.shape[0], nz=arr.shape[1])
        return

    raise ValueError(f"stats_update only supports rank-0, rank-1, or rank-2 values; got rank-{arr.ndim}.")


def stats_begin_budget(name: str, values, icol: int | None = None):
    """Begin budget accumulation using scalar, 1D, or 2D payloads."""
    arr = _as_stats_array(values)
    icol_idx = 0 if icol is None else int(icol) + 1

    if arr.ndim == 0:
        if icol_idx <= 0:
            raise ValueError("Scalar stats_begin_budget requires icol > 0.")
        clubb_f2py.f2py_stats_begin_budget_scalar(name, float(arr), icol_idx)
        return

    if arr.ndim == 1:
        clubb_f2py.f2py_stats_begin_budget_1d(name, arr, icol_idx, n=arr.shape[0])
        return

    if arr.ndim == 2:
        if icol_idx > 0:
            raise ValueError("rank-2 stats_begin_budget does not support icol.")
        clubb_f2py.f2py_stats_begin_budget_2d(name, f_arr(arr), ncol=arr.shape[0], nz=arr.shape[1])
        return

    raise ValueError(
        f"stats_begin_budget only supports rank-0, rank-1, or rank-2 values; got rank-{arr.ndim}."
    )


def stats_update_budget(name: str, values, icol: int | None = None,
                        level: int | None = None):
    """Update an active budget using scalar, 1D, or 2D payloads."""
    arr = _as_stats_array(values)
    icol_idx = 0 if icol is None else int(icol) + 1
    level_idx = 0 if level is None else int(level) + 1

    if arr.ndim == 0:
        if icol_idx <= 0:
            raise ValueError("Scalar stats_update_budget requires icol > 0.")
        clubb_f2py.f2py_stats_update_budget_scalar(name, float(arr), icol_idx, level_idx)
        return

    if arr.ndim == 1:
        if icol_idx > 0 and level_idx > 0:
            raise ValueError("rank-1 stats_update_budget accepts icol or level, not both.")
        if level_idx > 0:
            raise ValueError("rank-1 stats_update_budget does not support level-specific dispatch.")
        clubb_f2py.f2py_stats_update_budget_1d(name, arr, icol_idx, n=arr.shape[0])
        return

    if arr.ndim == 2:
        if icol_idx > 0 or level_idx > 0:
            raise ValueError("rank-2 stats_update_budget does not support icol or level.")
        clubb_f2py.f2py_stats_update_budget_2d(name, f_arr(arr), ncol=arr.shape[0], nz=arr.shape[1])
        return

    raise ValueError(
        f"stats_update_budget only supports rank-0, rank-1, or rank-2 values; got rank-{arr.ndim}."
    )


def stats_finalize_budget(name: str, values, icol: int | None = None,
                          l_count_sample: bool = True):
    """Finalize an active budget using scalar, 1D, or 2D payloads."""
    arr = _as_stats_array(values)
    icol_idx = 0 if icol is None else int(icol) + 1

    if arr.ndim == 0:
        if icol_idx <= 0:
            raise ValueError("Scalar stats_finalize_budget requires icol > 0.")
        clubb_f2py.f2py_stats_finalize_budget_scalar(name, float(arr), icol_idx, l_count_sample)
        return

    if arr.ndim == 1:
        clubb_f2py.f2py_stats_finalize_budget_1d(name, arr, icol_idx, l_count_sample, n=arr.shape[0])
        return

    if arr.ndim == 2:
        if icol_idx > 0:
            raise ValueError("rank-2 stats_finalize_budget does not support icol.")
        clubb_f2py.f2py_stats_finalize_budget_2d(
            name, f_arr(arr), l_count_sample, ncol=arr.shape[0], nz=arr.shape[1]
        )
        return

    raise ValueError(
        f"stats_finalize_budget only supports rank-0, rank-1, or rank-2 values; got rank-{arr.ndim}."
    )


def init_stats(registry_path: str, output_path: str, ncol: int,
               stats_tsamp: float, stats_tout: float, dt_main: float,
               day_in: int, month_in: int, year_in: int,
               time_initial: float,
               nzt: int, zt: np.ndarray, nzm: int, zm: np.ndarray,
               sclr_dim: int, edsclr_dim: int,
               err_info: ErrInfo,
               clubb_params: np.ndarray,
               param_names: list[str] | np.ndarray):
    """Initialize the stats system with a registry file."""
    set_fortran_err_info(err_info)

    encoded_param_names = np.asarray(
        [str(name).encode("ascii", errors="replace") for name in param_names],
        dtype="S28",
    )
    if encoded_param_names.shape != (clubb_params.shape[1],):
        raise ValueError("param_names length must match clubb_params.shape[1].")
    clubb_f2py.f2py_stats_init_with_params(
        registry_path, output_path,
        stats_tsamp, stats_tout, dt_main,
        day_in, month_in, year_in, time_initial,
        zt, zm, sclr_dim, edsclr_dim, f_arr(clubb_params), encoded_param_names,
        ncol=ncol, nzt=nzt, nzm=nzm,
    )
    return get_fortran_err_info()


def finalize_stats(err_info: ErrInfo):
    """Finalize the stats system and release resources."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_stats_finalize()
    return get_fortran_err_info()


def get_stats_config():
    """Get stats config scalars."""
    return clubb_f2py.f2py_get_stats_config()


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


def stats_update_scalar(name: str, values: float, icol: int | None = None,
                        level: int | None = None):
    """Update a scalar stats variable using Python indexing."""
    stats_update(name, values, icol=icol, level=level)


def stats_update_1d(name: str, n: int, values: np.ndarray, icol: int | None = None):
    """Update a 1D profile or surface stats variable using Python indexing."""
    stats_update(name, np.asarray(values)[:n], icol=icol)


def stats_update_2d(name: str, ncol: int, nz: int, values: np.ndarray):
    """Update a 2D field stats variable."""
    stats_update(name, np.asarray(values)[:ncol, :nz])


def stats_begin_budget_scalar(name: str, values: float, icol: int | None = None):
    """Begin budget tracking for a scalar variable using Python indexing."""
    stats_begin_budget(name, values, icol=icol)


def stats_begin_budget_1d(name: str, n: int, values: np.ndarray, icol: int | None = None):
    """Begin budget tracking for a 1D variable using Python indexing."""
    stats_begin_budget(name, np.asarray(values)[:n], icol=icol)


def stats_begin_budget_2d(name: str, ncol: int, nz: int, values: np.ndarray):
    """Begin budget tracking for a 2D variable."""
    stats_begin_budget(name, np.asarray(values)[:ncol, :nz])


def stats_update_budget_scalar(name: str, values: float, icol: int | None = None,
                               level: int | None = None):
    """Update budget accumulation for a scalar variable using Python indexing."""
    stats_update_budget(name, values, icol=icol, level=level)


def stats_update_budget_1d(name: str, n: int, values: np.ndarray, icol: int | None = None):
    """Update budget accumulation for a 1D variable using Python indexing."""
    stats_update_budget(name, np.asarray(values)[:n], icol=icol)


def stats_update_budget_2d(name: str, ncol: int, nz: int, values: np.ndarray):
    """Update budget accumulation for a 2D variable."""
    stats_update_budget(name, np.asarray(values)[:ncol, :nz])


def stats_finalize_budget_scalar(name: str, values: float, icol: int | None,
                                 l_count_sample: bool):
    """Finalize budget for a scalar variable using Python indexing."""
    stats_finalize_budget(name, values, icol=icol, l_count_sample=l_count_sample)


def stats_finalize_budget_1d(name: str, n: int, values: np.ndarray, icol: int | None,
                             l_count: bool):
    """Finalize budget for a 1D variable using Python indexing."""
    stats_finalize_budget(name, np.asarray(values)[:n], icol=icol, l_count_sample=l_count)


def stats_finalize_budget_2d(name: str, ncol: int, nz: int, values: np.ndarray,
                             l_count: bool):
    """Finalize budget for a 2D variable."""
    stats_finalize_budget(name, np.asarray(values)[:ncol, :nz], l_count_sample=l_count)


def var_on_stats_list(name: str) -> bool:
    """Check if a variable name is registered in the stats system."""
    return bool(clubb_f2py.f2py_var_on_stats_list(name))


def stats_update_grid(ncol: int, nzt: int, nzm: int, nrho: int,
                      zt_src: np.ndarray, zm_src: np.ndarray,
                      rho_vals: np.ndarray, rho_levels: np.ndarray,
                      p_sfc: np.ndarray):
    """Update adaptive grid remapping inputs."""
    clubb_f2py.f2py_stats_update_grid(
        f_arr(zt_src), f_arr(zm_src),
        f_arr(rho_vals), f_arr(rho_levels),
        f_arr(p_sfc), ncol=ncol, nzt=nzt, nzm=nzm, nrho=nrho,
    )


def stats_lh_samples_init(num_samples: int, nzt: int,
                          n_nl: int, n_u: int,
                          nl_names: list, u_names: list,
                          zt_vals: np.ndarray,
                          err_info: ErrInfo):
    """Initialize SILHS sample output variables."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_stats_lh_samples_init(
        num_samples, nl_names, u_names, zt_vals, nzt=nzt, n_nl=n_nl, n_u=n_u,
    )


def stats_lh_samples_write_lognormal(
    ncol: int, nsamp: int, nzt: int, nvars: int,
    samples: np.ndarray,
    err_info: ErrInfo,
):
    """Write lognormal SILHS sample data."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_stats_lh_samples_write_ln(
        f_arr(samples), ncol=ncol, nsamp=nsamp, nzt=nzt, nvars=nvars)


def stats_lh_samples_write_uniform(ncol: int, nsamp: int, nzt: int, nvars: int,
                                   uniform_vals: np.ndarray,
                                   mixture_comp: np.ndarray,
                                   sample_weights: np.ndarray,
                                   err_info: ErrInfo):
    """Write uniform SILHS sample data."""
    set_fortran_err_info(err_info)
    clubb_f2py.f2py_stats_lh_samples_write_u(
        f_arr(uniform_vals),
        f_arr(mixture_comp),
        f_arr(sample_weights), ncol=ncol, nsamp=nsamp, nzt=nzt, nvars=nvars,
    )
