"""Read CLUBB sounding files (*_sounding.in).

Format: whitespace-separated columns. First line is a header describing
columns. Missing values are -999.9. Typical columns:

    z[m]  thm[K]  'rt[kg/kg]'  'u[m/s]'  'v[m/s]'  'w[m/s]'  'ug[m/s]'  'vg[m/s]'

Header tokens determine what type of temperature and altitude are given:
  - Altitude: "z[m]" or "Press[Pa]"
  - Temperature: "thm[K]" (potential), "thlm[K]" (liquid potential), "T[K]" (absolute)
  - Subsidence: "w[m/s]" or "omega[Pa/s]"

This module reads the raw columns and returns them as numpy arrays,
along with the header metadata needed to interpret them.
"""
import numpy as np
from pathlib import Path

MISSING = -999.9


def read_sounding(path: str) -> dict:
    """Read a sounding file and return a dict of column arrays + metadata.

    Returns:
        dict with keys:
            'alt_type': str ('z[m]' or 'Press[Pa]')
            'theta_type': str ('thm[K]', 'thlm[K]', or 'T[K]')
            'subs_type': str ('w[m/s]' or 'omega[Pa/s]')
            'z': np.ndarray of altitudes/pressures
            'theta': np.ndarray of temperature values
            'rt': np.ndarray of total water mixing ratio
            'u': np.ndarray of eastward wind
            'v': np.ndarray of northward wind
            'w': np.ndarray of vertical velocity / omega
            'ug': np.ndarray of geostrophic u wind
            'vg': np.ndarray of geostrophic v wind
            'columns': list of column names from header
            'data': np.ndarray of shape (nlevels, ncols) raw data
    """
    lines = Path(path).read_text().splitlines()

    # Skip comment lines (starting with !)
    data_lines = [l for l in lines if l.strip() and not l.strip().startswith('!')]

    # First non-comment line is the header
    header_line = data_lines[0]
    # Parse column names — may be quoted with single quotes
    cols = header_line.split()
    col_names = [c.strip("'\"") for c in cols]

    # Read data
    rows = []
    for line in data_lines[1:]:
        vals = line.split()
        rows.append([float(v) for v in vals])

    data = np.array(rows)

    def _find_idx(options: tuple[str, ...]) -> int | None:
        for i, name in enumerate(col_names):
            for opt in options:
                if name.lower() == opt.lower():
                    return i
        return None

    # Identify required column types from names (order in file can vary by case).
    alt_idx = _find_idx(("z[m]", "press[pa]"))
    if alt_idx is None:
        raise ValueError(f"Could not find altitude column in sounding header: {col_names}")
    alt_type = col_names[alt_idx]

    theta_idx = _find_idx(("thm[k]", "thlm[k]", "t[k]"))
    if theta_idx is None:
        raise ValueError(f"Could not find temperature column in sounding header: {col_names}")
    theta_type = col_names[theta_idx]

    rt_idx = _find_idx(("rt[kg/kg]",))
    u_idx = _find_idx(("u[m/s]",))
    v_idx = _find_idx(("v[m/s]",))
    ug_idx = _find_idx(("ug[m/s]",))
    vg_idx = _find_idx(("vg[m/s]",))

    # Find subsidence column
    subs_idx = _find_idx(("w[m/s]", "omega[pa/s]"))
    subs_type = col_names[subs_idx] if subs_idx is not None else "w[m/s]"

    result = {
        'alt_type': alt_type,
        'theta_type': theta_type,
        'subs_type': subs_type,
        'columns': col_names,
        'data': data,
        'z': data[:, alt_idx],
        'theta': data[:, theta_idx],
        'rt': data[:, rt_idx] if rt_idx is not None else np.full(len(data), MISSING),
        'u': data[:, u_idx] if u_idx is not None else np.full(len(data), MISSING),
        'v': data[:, v_idx] if v_idx is not None else np.full(len(data), MISSING),
        'w': data[:, subs_idx] if subs_idx is not None else np.full(len(data), MISSING),
        'ug': data[:, ug_idx] if ug_idx is not None else np.full(len(data), MISSING),
        'vg': data[:, vg_idx] if vg_idx is not None else np.full(len(data), MISSING),
    }

    return result


def _steffen_interp_1d(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    """Steffen monotone cubic interpolation (matches CLUBB's mono_cubic_interp path)."""
    n = len(x)
    if n < 3:
        return np.interp(x_new, x, y)

    out = np.empty_like(x_new, dtype=np.float64)
    for j, xj in enumerate(x_new):
        # Find first level above (or equal to) the target.
        k_above = int(np.searchsorted(x, xj, side='left'))
        if k_above <= 0:
            # Linear extrapolation below first point.
            out[j] = y[0] + (xj - x[0]) * (y[1] - y[0]) / (x[1] - x[0])
            continue
        if k_above >= n:
            # Linear extrapolation above last point.
            out[j] = y[-2] + (xj - x[-2]) * (y[-1] - y[-2]) / (x[-1] - x[-2])
            continue

        # Convert Fortran-style stencil selection to 0-based indexing.
        if k_above == 1:
            km1, k00, kp1, kp2 = 0, 0, 1, 2
        else:
            km1 = k_above - 2
            k00 = k_above - 1
            kp1 = k_above
            kp2 = min(k_above + 1, n - 1)
            if x[k_above] >= x[-1]:
                km1, k00, kp1, kp2 = n - 3, n - 2, n - 1, n - 1

        zm1, z00, zp1, zp2 = x[km1], x[k00], x[kp1], x[kp2]
        fm1, f00, fp1, fp2 = y[km1], y[k00], y[kp1], y[kp2]

        h00 = zp1 - z00
        if h00 == 0.0:
            out[j] = f00
            continue

        if km1 == k00:
            s00 = (fp1 - f00) / (zp1 - z00)
            sp1 = (fp2 - fp1) / (zp2 - zp1)
            dfdx00 = s00
            pp1 = (s00 * (zp2 - zp1) + sp1 * h00) / (h00 + (zp2 - zp1))
            dfdxp1 = (np.sign(s00) + np.sign(sp1)) * min(
                abs(s00), abs(sp1), 0.5 * abs(pp1)
            )
        elif kp1 == kp2:
            hm1 = z00 - zm1
            s00 = (fp1 - f00) / (zp1 - z00)
            sm1 = (f00 - fm1) / (z00 - zm1)
            p00 = (sm1 * h00 + s00 * hm1) / (hm1 + h00)
            dfdx00 = (np.sign(sm1) + np.sign(s00)) * min(
                abs(sm1), abs(s00), 0.5 * abs(p00)
            )
            dfdxp1 = s00
        else:
            hm1 = z00 - zm1
            hp1 = zp2 - zp1
            sm1 = (f00 - fm1) / (z00 - zm1)
            s00 = (fp1 - f00) / (zp1 - z00)
            sp1 = (fp2 - fp1) / (zp2 - zp1)
            p00 = (sm1 * h00 + s00 * hm1) / (hm1 + h00)
            pp1 = (s00 * hp1 + sp1 * h00) / (h00 + hp1)
            dfdx00 = (np.sign(sm1) + np.sign(s00)) * min(
                abs(sm1), abs(s00), 0.5 * abs(p00)
            )
            dfdxp1 = (np.sign(s00) + np.sign(sp1)) * min(
                abs(s00), abs(sp1), 0.5 * abs(pp1)
            )

        c1 = (dfdx00 + dfdxp1 - 2.0 * s00) / (h00**2)
        c2 = (3.0 * s00 - 2.0 * dfdx00 - dfdxp1) / h00
        zprime = xj - z00
        out[j] = f00 + zprime * (dfdx00 + zprime * (c2 + zprime * c1))

    return out


def _fill_blanks(z_snd: np.ndarray, vals: np.ndarray) -> np.ndarray:
    """Fill missing (-999.9) values by linear interpolation at original z levels.

    Mirrors Fortran fill_blanks_one_dim_vars / zlinterp_fnc behaviour:
    - Within the range of known values: linear interpolation
    - Below or above the known range: set to 0.0
    """
    valid = vals > MISSING + 1.0
    n_valid = np.sum(valid)
    if n_valid == 0:
        return np.zeros_like(vals)
    if n_valid == len(vals):
        return vals.copy()

    x_known = z_snd[valid]
    y_known = vals[valid]

    filled = np.zeros_like(vals)
    for i, z in enumerate(z_snd):
        if z < x_known[0] or z > x_known[-1]:
            filled[i] = 0.0
        else:
            filled[i] = np.interp(z, x_known, y_known)

    return filled


def interpolate_sounding(snd: dict, grid_z: np.ndarray, use_cubic: bool = False) -> dict:
    """Interpolate sounding profiles onto a model grid.

    Mirrors the Fortran pipeline: first fill_blanks (linear fill at original
    sounding z-levels), then interpolate (linear or cubic) from the full set
    of sounding levels to the model grid.

    Args:
        snd: dict from read_sounding()
        grid_z: 1D array of target altitudes [m]
        use_cubic: use Steffen monotone cubic interpolation

    Returns:
        dict with same keys as snd but interpolated to grid_z.
    """
    result = {'z': grid_z, 'alt_type': snd['alt_type'],
              'theta_type': snd['theta_type'], 'subs_type': snd['subs_type']}

    z_snd = snd['z']

    for key in ('theta', 'rt', 'u', 'v', 'w', 'ug', 'vg'):
        vals = snd[key]
        # Mask out missing values
        valid = vals > MISSING + 1.0  # not exactly -999.9
        if np.sum(valid) < 2:
            result[key] = np.full_like(grid_z, 0.0)
            continue

        # Step 1: fill blanks at original sounding z-levels (matches Fortran)
        filled = _fill_blanks(z_snd, vals)

        # Step 2: interpolate from filled sounding to model grid
        if use_cubic:
            result[key] = _steffen_interp_1d(z_snd, filled, grid_z)
        else:
            result[key] = np.interp(grid_z, z_snd, filled)

    return result


def read_scalar_sounding(path: str, nvars: int) -> dict:
    """Read a scalar sounding file (*_sclr_sounding.in or *_edsclr_sounding.in).

    Format: header names on the first non-comment line, followed by data rows.
    Unlike *_sounding.in, these files do not include a height column; profiles
    are interpreted on the same source levels as the main sounding file.

    Returns:
        dict with keys:
            'names': list[str] of length nvars
            'data': np.ndarray of shape (nlevels, nvars)
    """
    lines = Path(path).read_text().splitlines()
    data_lines = [l for l in lines if l.strip() and not l.strip().startswith('!')]
    if not data_lines:
        raise ValueError(f"Scalar sounding file is empty: {path}")

    header = [c.strip("'\"") for c in data_lines[0].split()]
    if len(header) < nvars:
        raise ValueError(
            f"Scalar sounding header has {len(header)} columns, expected at least {nvars}: {path}"
        )
    names = header[:nvars]

    rows = []
    for line in data_lines[1:]:
        vals = [float(v) for v in line.split()]
        if len(vals) < nvars:
            raise ValueError(
                f"Scalar sounding row has {len(vals)} values, expected at least {nvars}: {path}"
            )
        rows.append(vals[:nvars])

    data = np.asarray(rows, dtype=np.float64)
    return {'names': names, 'data': data}


def interpolate_scalar_sounding(
    z_snd: np.ndarray, scalar_data: np.ndarray, grid_z: np.ndarray, use_cubic: bool = False
) -> np.ndarray:
    """Interpolate scalar-sounding columns from sounding levels to model zt levels.

    The scalar source profile is padded with zeros if it has fewer rows than
    `z_snd`, matching Fortran's zero-initialized storage before assignment.
    """
    n_src = int(len(z_snd))
    ncols = int(scalar_data.shape[1])
    prof = np.zeros((n_src, ncols), dtype=np.float64)
    n_copy = min(n_src, int(scalar_data.shape[0]))
    if n_copy > 0:
        prof[:n_copy, :] = scalar_data[:n_copy, :]

    out = np.zeros((int(len(grid_z)), ncols), dtype=np.float64)
    for j in range(ncols):
        if use_cubic:
            out[:, j] = _steffen_interp_1d(z_snd, prof[:, j], grid_z)
        else:
            out[:, j] = np.interp(grid_z, z_snd, prof[:, j])

    return out
