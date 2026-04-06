"""Read CLUBB grid files (*.grd).

Grid files are simple single-column files containing altitude values [m],
one per line. Used for stretched grids (grid_type 2 or 3).

For grid_type = 1 (evenly-spaced), no grid file is needed — the grid is
generated from deltaz, zm_init, and zm_top.
"""
import numpy as np
from pathlib import Path


def read_grid_file(path: str) -> np.ndarray:
    """Read a grid file and return altitude array.

    Args:
        path: Path to the .grd file

    Returns:
        1D numpy array of altitude values [m]
    """
    lines = Path(path).read_text().splitlines()
    values = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith('!') or line.startswith('#'):
            continue
        values.append(float(line))
    return np.array(values)
