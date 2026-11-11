"""Shared Fortran-reference setup for the JAX F2PY unit tests.

The F2PY unit tests call compiled CLUBB through the Python API (`clubb_python`)
or, where no API wrapper exists, the raw `clubb_f2py` binding. Several Fortran
routines read grid, error-info or configuration state from module storage
rather than from their argument lists, so that state must be set before the
call. Use the helpers here instead of calling raw bindings on uninitialized
storage, which can corrupt the heap and crash a later, unrelated check.
"""

import numpy as np

from clubb_python.CLUBB_core.err_info_type_module import init_err_info
from clubb_python.CLUBB_core.grid_class import setup_grid as fortran_setup_grid
from clubb_python.CLUBB_core.model_flags import get_default_config_flags, init_config_flags
from clubb_python.derived_types.err_info import ErrInfo

from clubb_jax.src.CLUBB_core.grid_class import setup_grid as jax_setup_grid
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats


def init_fortran_state(ngrdcol):
    """Allocate Fortran error info and load the default configuration flags."""
    init_err_info(ngrdcol)
    init_config_flags(get_default_config_flags())


def setup_matching_grids(ngrdcol, deltaz, zm_top, zm_init=0.0, grid_type=1,
                         momentum_heights=None, thermodynamic_heights=None):
    """Build a JAX grid and the identical Fortran grid in module storage.

    The Fortran grid is built from the JAX heights (grid_type 2 for an evenly
    spaced JAX grid, otherwise the given type and heights), so both
    implementations see the same levels. Returns the JAX grid and the Fortran
    `Grid` returned by `setup_grid`.
    """
    jgr = jax_setup_grid(ngrdcol=ngrdcol, deltaz=deltaz, zm_init=zm_init, zm_top=zm_top,
                         grid_type=grid_type, momentum_heights=momentum_heights,
                         thermodynamic_heights=thermodynamic_heights)
    zm = np.asfortranarray(np.asarray(jgr.zm))
    zt = np.asfortranarray(np.asarray(jgr.zt))
    nzm = zm.shape[1]
    init_fortran_state(ngrdcol)
    fgr, err = fortran_setup_grid(
        nzmax=nzm, ngrdcol=ngrdcol, sfc_elevation=np.zeros(ngrdcol), l_implemented=False,
        l_ascending_grid=True, grid_type=2 if grid_type == 1 else grid_type,
        deltaz=np.full(ngrdcol, float(deltaz)), zm_init=np.full(ngrdcol, float(zm_init)),
        zm_top=zm[:, -1].copy(), momentum_heights=zm, thermodynamic_heights=zt,
        err_info=ErrInfo(ngrdcol=ngrdcol))
    assert not err.is_fatal(), "Fortran setup_grid reported a fatal error"
    return jgr, fgr


def no_stats(ngrdcol):
    """Return a non-sampling JAX stats object; the Fortran reference runs without stats output."""
    return JaxStats.empty(l_sample=False, names=(), ncol=ngrdcol, max_nlev=1)


def _checked_pair(name, jax_out, fortran_out):
    """Return both outputs as float64 arrays after checking that shapes agree and values are finite."""
    j, f = np.asarray(jax_out, dtype=np.float64), np.asarray(fortran_out, dtype=np.float64)
    assert j.shape == f.shape, f"{name}: JAX shape {j.shape} != Fortran shape {f.shape}"
    n_j, n_f = int(np.count_nonzero(~np.isfinite(j))), int(np.count_nonzero(~np.isfinite(f)))
    assert n_j == 0 and n_f == 0, f"{name}: {n_j} nonfinite JAX values, {n_f} nonfinite Fortran values"
    return j, f


def max_abs_diff(name, jax_out, fortran_out):
    """Return the largest |JAX - Fortran| difference after checking shapes and finiteness.

    Accepts arrays or matching tuples/lists of arrays. The shape check keeps an
    accidental scalar or single-column result from broadcasting into a pass.
    Nonfinite values fail here because Python's max(worst, nan) keeps `worst`,
    so a NaN difference accumulated in a test loop would otherwise pass.
    """
    if isinstance(fortran_out, (tuple, list)):
        assert isinstance(jax_out, (tuple, list)) and len(jax_out) == len(fortran_out), (
            f"{name}: JAX returned {type(jax_out).__name__} of length "
            f"{len(jax_out) if isinstance(jax_out, (tuple, list)) else '-'}, Fortran {len(fortran_out)} outputs")
        return max((max_abs_diff(f"{name}[{i}]", j, f) for i, (j, f) in enumerate(zip(jax_out, fortran_out))),
                   default=0.0)
    j, f = _checked_pair(name, jax_out, fortran_out)
    if j.size == 0:
        return 0.0
    return float(np.max(np.abs(j - f)))


def max_rel_diff(name, jax_out, fortran_out, floor):
    """Return max |JAX - Fortran| / max(|Fortran|, floor) after checking shapes and finiteness."""
    j, f = _checked_pair(name, jax_out, fortran_out)
    if j.size == 0:
        return 0.0
    return float(np.max(np.abs(j - f) / np.maximum(np.abs(f), floor)))
