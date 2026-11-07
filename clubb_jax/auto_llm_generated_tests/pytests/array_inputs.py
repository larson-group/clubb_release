"""Supply JAX array inputs at numerical-test call sites.

The numerical ports require callers to create their arrays. Independent NumPy
references still build host arrays, so call_jax converts those fixtures before
calling a port. Existing JAX arrays/tracers pass through unchanged. Python real
scalars become scalar arrays; static integer/boolean controls and model objects'
static pytree metadata keep their types. No axes or precision are changed here.
"""
import jax
import jax.numpy as jnp
import numpy as np


def call_jax(routine, *args, **kwargs):
    """Call a numerical port with JAX leaves while retaining independent host fixtures."""
    def array_input(value):
        if isinstance(value, (np.ndarray, np.floating, float, complex)):
            return jnp.asarray(value)
        return value
    args, kwargs = jax.tree.map(array_input, (args, kwargs))
    return routine(*args, **kwargs)
