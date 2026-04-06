"""Python radiation driver for simplified CLUBB radiation schemes.

Supported schemes:
- none
- simplified
- simplified_bomex

The BUGSrad and LBA branches from Fortran radiation_module are intentionally
excluded; those schemes are handled by the Fortran driver.
"""

import numpy as np

from clubb_python import clubb_api


_CP = 1004.67
_EPS = np.finfo(np.float64).eps
_LS_DIV = 3.75e-6


def simple_rad_bomex(zt: np.ndarray):
    """Compute BOMEX simplified radiation heating rate."""
    return np.where(
        zt < 1500.0,
        -2.315e-5,
        np.where(
            zt < 2500.0,
            -2.315e-5 + 2.315e-5 * (zt - 1500.0) / 1000.0,
            0.0,
        ),
    )


def advance_radiation(state: dict, time_current: float, l_sample: bool = False):
    """Advance radiation tendencies for the current timestep."""
    scheme = str(state['rad_scheme']).strip().lower()

    if scheme == "none":
        state['radht'].fill(0.0)
    elif scheme == "simplified_bomex":
        state['radht'][:] = simple_rad_bomex(state['gr'].zt)
    elif scheme == "simplified":
        _advance_simplified_radiation(state, time_current, l_sample=l_sample)
    else:
        raise ValueError(
            f"Unsupported rad_scheme={scheme!r} in Python radiation driver. "
            "Supported: none, simplified, simplified_bomex."
        )

    if l_sample:
        clubb_api.stats_update("radht", state['radht'])


def _advance_simplified_radiation(state: dict, time_current: float, l_sample: bool = False):
    """Python port of the simplified branch from radiation_module."""
    cfg = state['cfg']
    gr = state['gr']
    ngrdcol = state['ngrdcol']
    nzm = state['nzm']
    nzt = state['nzt']

    l_sw_radiation = bool(cfg.get('l_sw_radiation', False))

    frad_sw = np.zeros((ngrdcol, nzm), dtype=np.float64)
    radht_sw = np.zeros((ngrdcol, nzt), dtype=np.float64)

    amu0 = _compute_amu0(state, time_current)
    if l_sw_radiation and amu0 > 0.0:
        fs0 = _compute_fs0(cfg, amu0)

        # Call sunray_sw to get Frad_SW on momentum levels (ngrdcol, nzm).
        frad_sw = clubb_api.sunray_sw(
            gr=gr,
            ngrdcol=ngrdcol,
            nzt=nzt,
            fs0=fs0,
            amu0=amu0,
            rho=state['rho'],
            rcm=state['rcm'],
        )

        # Compute radht_SW from Frad_SW: ddzm(Frad_SW) / (rho * Cp)
        # This matches the pattern in radiation_module.F90.
        radht_sw = (
            -(frad_sw[:, 1:] - frad_sw[:, :-1]) * gr.invrs_dzt
            / (state['rho'] * _CP)
        )

    frad_lw, radht_lw = _simple_rad_lw(
        gr=gr,
        ngrdcol=ngrdcol,
        rho=state['rho'],
        rho_zm=state['rho_zm'],
        rtm=state['rtm'],
        rcm=state['rcm'],
        exner=state['exner'],
        f0=float(cfg.get('f0', 0.0)),
        f1=float(cfg.get('f1', 0.0)),
        kappa=float(cfg.get('kappa', 0.0)),
        l_rad_above_cloud=bool(cfg.get('l_rad_above_cloud', False)),
        l_sample=l_sample,
    )

    frad_total = frad_sw + frad_lw
    state['radht'][:] = radht_sw + radht_lw
    state['Frad'] = frad_total
    state['Frad_SW'] = frad_sw
    state['Frad_LW'] = frad_lw
    state['radht_SW'] = radht_sw
    state['radht_LW'] = radht_lw

    if l_sample:
        clubb_api.stats_update("Frad", frad_total)
        clubb_api.stats_update("Frad_SW", frad_sw)
        clubb_api.stats_update("Frad_LW", frad_lw)
        clubb_api.stats_update("radht_SW", radht_sw)
        clubb_api.stats_update("radht_LW", radht_lw)


def _simple_rad_lw(
    gr,
    ngrdcol: int,
    rho: np.ndarray,
    rho_zm: np.ndarray,
    rtm: np.ndarray,
    rcm: np.ndarray,
    exner: np.ndarray,
    f0: float,
    f1: float,
    kappa: float,
    l_rad_above_cloud: bool,
    l_sample: bool,
):
    """Port of simple_rad from simple_rad_module.F90.

    Operates on full 2D (ngrdcol, nz) arrays matching the Fortran version.
    """
    nzm = gr.zm.shape[1]
    nzt = gr.zt.shape[1]

    # Liquid water path on momentum levels — accumulates top-down.
    lwp = _liq_water_path(ngrdcol, nzm, nzt, rho, rcm, gr.invrs_dzt)

    # Longwave radiative flux.
    if f1 > _EPS:
        frad_lw = f0 * np.exp(-kappa * lwp) + f1 * np.exp(-kappa * (lwp[:, 0:1] - lwp))
    else:
        frad_lw = f0 * np.exp(-kappa * lwp)

    if l_rad_above_cloud:
        # Find z_i per column where rtm crosses 8 g/kg.
        z_i = np.zeros(ngrdcol, dtype=np.float64)
        for i in range(ngrdcol):
            k_iso = 0
            while k_iso < nzt and rtm[i, k_iso] > 8.0e-3:
                k_iso += 1
            if k_iso == 0 or k_iso > nzt:
                z_i[i] = 0.0
                continue
            if k_iso == nzt:
                k_iso = nzt - 1
            z_i[i] = _linear_interp(
                8.0e-3, rtm[i, k_iso], rtm[i, k_iso - 1],
                gr.zt[i, k_iso], gr.zt[i, k_iso - 1],
            )

        # Heaviside step function and above-cloud contribution.
        # z_i shape (ngrdcol,) -> broadcast against zm (ngrdcol, nzm).
        dz = gr.zm - z_i[:, np.newaxis]
        heaviside = np.where(dz < -_EPS, 0.0, np.where(dz > _EPS, 1.0, 0.5))

        pos = heaviside > 0.0
        if np.any(pos):
            dz_pos = np.maximum(dz[pos], 0.0)
            z_i_broad = np.broadcast_to(z_i[:, np.newaxis], (ngrdcol, nzm))[pos]
            frad_lw[pos] += (
                rho_zm[pos] * _CP * _LS_DIV * heaviside[pos]
                * (0.25 * (dz_pos ** (4.0 / 3.0)) + z_i_broad * (dz_pos ** (1.0 / 3.0)))
            )

        if l_sample:
            clubb_api.stats_update("z_inversion", z_i)

    # Radiative heating rate on thermodynamic levels.
    radht_lw = (
        (1.0 / exner) * (-1.0 / (_CP * rho))
        * (frad_lw[:, 1:] - frad_lw[:, :-1]) * gr.invrs_dzt
    )

    return frad_lw, radht_lw


def _liq_water_path(ngrdcol: int, nzm: int, nzt: int,
                    rho: np.ndarray, rcm: np.ndarray,
                    invrs_dzt: np.ndarray) -> np.ndarray:
    """Compute liquid water path on momentum levels, matching liq_water_path in Fortran."""
    lwp = np.zeros((ngrdcol, nzm), dtype=np.float64)
    for k in range(nzm - 2, -1, -1):
        lwp[:, k] = lwp[:, k + 1] + rcm[:, k] * rho[:, k] / invrs_dzt[:, k]
    return lwp


def _linear_interp(x: float, x_high: float, x_low: float,
                   y_high: float, y_low: float) -> float:
    """Linear interpolation in the same argument order used by Fortran."""
    denom = x_high - x_low
    if abs(denom) < _EPS:
        return 0.5 * (y_high + y_low)
    return y_low + (x - x_low) * (y_high - y_low) / denom


def _compute_amu0(state: dict, time_current: float) -> float:
    """Compute cosine of solar zenith angle for simplified radiation."""
    cfg = state['cfg']
    l_fix = bool(cfg.get('l_fix_cos_solar_zen', False))

    cos_vals = _as_1d_float(cfg.get('cos_solar_zen_values', [0.0]))
    if l_fix:
        if cos_vals.size == 1:
            return float(cos_vals[0])
        times = _as_1d_float(cfg.get('cos_solar_zen_times', [time_current]))
        if times.size != cos_vals.size:
            return float(cos_vals[0])
        idx = int(np.searchsorted(times, time_current, side='left'))
        if idx >= times.size:
            raise ValueError("time_current exceeds provided cos_solar_zen_times range.")
        return float(cos_vals[idx])

    return clubb_api.cos_solar_zen(
        int(cfg.get('day', 1)),
        int(cfg.get('month', 1)),
        int(cfg.get('year', 2000)),
        float(time_current),
        float(cfg.get('lat_vals', 0.0)),
        float(cfg.get('lon_vals', 0.0)),
    )


def _compute_fs0(cfg: dict, amu0: float) -> float:
    """Interpolate solar flux on cosine solar zenith lookup table."""
    fs_values = _as_1d_float(cfg.get('fs_values', [0.0]))
    cos_values = _as_1d_float(cfg.get('cos_solar_zen_values', [0.0]))

    if fs_values.size == 0:
        return 0.0
    if fs_values.size == 1 or cos_values.size != fs_values.size:
        return float(fs_values[0])

    return float(np.interp(amu0, cos_values, fs_values))


def _as_1d_float(value):
    arr = np.asarray(value, dtype=np.float64)
    if arr.ndim == 0:
        arr = arr.reshape(1)
    return arr
