"""Sampling storage from latin_hypercube_arrays.F90.

JAX adaptation: threadprivate permutation storage and prior_iter are explicit
pytree leaves carried by the caller; compiled routines never mutate globals.
"""

from typing import NamedTuple


class LatinHypercubeArrays(NamedTuple):
    """Case-owned replacement for the source threadprivate sampling storage."""

    one_height_time_matrix: object  # Matrix of random integers: (nt_repeat, n_vars).
    prior_iter: object              # Prior model iteration, for diagnostics.


# -----------------------------------------------------------------------------
def cleanup_latin_hypercube_arrays():
    """Deallocate Latin-hypercube arrays in the source interface.

    JAX adaptation: immutable storage is released with its owning case state;
    there is no module-global allocation to deallocate here.
    """
    return None
