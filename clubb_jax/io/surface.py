"""Read CLUBB surface files (*_sfc.in).

Format: whitespace-separated columns with a header line.
Typical columns:
    Time[s]  'wpqtp_sfc[(kg/kg)m/s]'  'wpthlp_sfc[mK/s]'  T_sfc[K]

Values are linearly interpolated in time during the simulation.
"""
import numpy as np
from pathlib import Path


def read_surface(path: str) -> dict:
    """Read a surface file and return time series of surface fluxes.

    Returns:
        dict with keys:
            'time': np.ndarray of times [s]
            'wprtp_sfc': np.ndarray of moisture flux [(kg/kg)(m/s)]
            'wpthlp_sfc': np.ndarray of heat flux [m K/s]
            't_sfc': np.ndarray of surface temperature [K]
            'columns': list of column names
    """
    lines = Path(path).read_text().splitlines()
    data_lines = [l for l in lines if l.strip() and not l.strip().startswith('!')]

    header_line = data_lines[0]
    cols = header_line.split()
    col_names = [c.strip("'\"") for c in cols]

    rows = []
    for line in data_lines[1:]:
        vals = line.split()
        rows.append([float(v) for v in vals])

    data = np.array(rows)

    result = {
        'columns': col_names,
        'time': data[:, 0],
        'wpqtp_sfc': data[:, 1] if data.shape[1] > 1 else np.zeros(len(data)),
        'wpthlp_sfc': data[:, 2] if data.shape[1] > 2 else np.zeros(len(data)),
        't_sfc': data[:, 3] if data.shape[1] > 3 else np.full(len(data), 300.0),
    }

    return result


def interp_surface(sfc: dict, time: float) -> dict:
    """Interpolate surface values to a given simulation time.

    Returns:
        dict with scalar values for each surface variable at time t.
    """
    return {
        'wpqtp_sfc': float(np.interp(time, sfc['time'], sfc['wpqtp_sfc'])),
        'wpthlp_sfc': float(np.interp(time, sfc['time'], sfc['wpthlp_sfc'])),
        't_sfc': float(np.interp(time, sfc['time'], sfc['t_sfc'])),
    }
