"""JAX port of extended_atmosphere_module.F90 — BUGSrad extended-atmosphere grid bounds.

Mirrors clubb_release/src/Radiation/extended_atmosphere_module.F90: `determine_extended_atmos_bounds` finds the
bottom/top levels of the standard-atmosphere extension above the model top that the radiation grid must span,
plus the linear-interpolation buffer count. (The Fortran `finalize_extended_atm` deallocation has no JAX
analog.) Used by the BUGSrad driver; bugsrad_driver.py imports it back, mirroring the Fortran
`use extended_atmosphere_module`.

Pure-numpy → bit-identical. Validated by `tests/test_bugsrad.py`.
"""

from __future__ import annotations

import numpy as np

# extended_atmosphere_module.F90:117-119 `use constants_clubb, only: pascal_per_mb`
from clubb_jax.src.CLUBB_core.constants_clubb import pascal_per_mb as PASCAL_PER_MB


def determine_extended_atmos_bounds(zm_grid, zm_grid_spacing, p_in_Pa_zm, radiation_top, ext):
    """Port of extended_atmosphere_module.F90:determine_extended_atmos_bounds. `zm_grid`=(nzm,) momentum
    altitudes [m]; `zm_grid_spacing`=(nzm-1,) Δz; `p_in_Pa_zm`=(nzm,) momentum-level pressure [Pa];
    `radiation_top` [m]; `ext`=load_extended_std_atm() dict. Returns (bottom_level, top_level, range_size,
    lin_int_buffer) — bottom/top are 1-based Fortran indices into the std-atm arrays."""
    extended_alt = np.asarray(ext["alt"]); extended_p_in_mb = np.asarray(ext["p_in_mb"])
    ext_dim = extended_alt.shape[0]
    grid_size = zm_grid.shape[0]
    zm_top = float(zm_grid[grid_size - 1])
    p_top = float(p_in_Pa_zm[grid_size - 1])

    if radiation_top < zm_top:
        raise ValueError("top of the radiation grid is below the top of the computational grid")

    j = 1  # 1-based
    while extended_alt[j - 1] < zm_top and j < ext_dim:
        j += 1
    while p_top < extended_p_in_mb[j - 1] * PASCAL_PER_MB:
        j += 1
    if extended_alt[j - 1] < zm_top:
        raise ValueError("Extended atmosphere is below the top of the computational grid")
    if extended_alt[ext_dim - 1] < radiation_top:
        raise ValueError("extension data does not reach radiation_top")
    if p_top < extended_p_in_mb[j - 1] * PASCAL_PER_MB:
        raise ValueError("pressure at top of computational grid less than pressure at base of radiative grid")

    k = 1
    if j <= ext_dim:
        while extended_alt[k - 1] < radiation_top and k < ext_dim:
            k += 1
        if extended_alt[k - 1] > radiation_top:
            k -= 1
    else:
        k = j

    bottom_level, top_level = j, k
    range_size = k - j + 1
    if range_size < 1:
        raise ValueError("radiation top below computational grid")

    extended_bottom = float(extended_alt[bottom_level - 1])
    dz10 = (extended_bottom - zm_top) / 10.0
    dz_model = float(zm_grid_spacing[grid_size - 2])      # zm_grid_spacing(grid_size-1), 1-based
    dz = max(dz10, dz_model)
    lin_int_buffer = int((extended_bottom - zm_top) / dz)  # Fortran int() truncates toward zero
    return bottom_level, top_level, range_size, lin_int_buffer
