"""Tracer-aware NumPy compatibility for legacy radiation and microphysics.

Porting deviations:
- This is JAX-only infrastructure with no corresponding Fortran source file.
  TODO(port-mirror): Remove it after the remaining radiation and microphysics
  callers use JAX arrays and functional updates directly.
"""
import numpy as np
import jax
import jax.numpy as jnp

_np_asarray = np.asarray


def _is_tracer_arg(a):
    if isinstance(a, jax.core.Tracer):
        return True
    if isinstance(a, (list, tuple)):
        return any(isinstance(e, jax.core.Tracer) for e in a)
    return False


def _asarray(x, dtype=None):
    """Tracer-transparent ``np.asarray``: jnp under a JAX trace, else exactly numpy (bit-identical)."""
    if isinstance(x, jax.core.Tracer):
        return jnp.asarray(x, dtype)
    return _np_asarray(x, dtype)


def _safe_sqrt(x):
    """``sqrt(max(x,0))`` with a finite (0) gradient at x<=0 (double-where) — the bare
    ``jnp.sqrt(jnp.maximum(0,x))`` has an inf reverse-mode gradient at the clip. Forward-identical."""
    xp = jnp.maximum(x, 0.0)
    safe = jnp.where(xp > 0.0, xp, 1.0)
    return jnp.where(xp > 0.0, jnp.sqrt(safe), 0.0)


class _TracerXP:
    """numpy shim: ``_xp.<fn>(...)`` dispatches to jnp when any arg (or list/tuple element) is a tracer,
    else to numpy (bit-identical). Non-callables / dtypes / constants fall through to numpy."""
    def __getattr__(self, name):
        np_fn = getattr(np, name)
        jnp_fn = getattr(jnp, name, None)
        if jnp_fn is None or not callable(np_fn):
            return np_fn

        def dispatch(*args, **kw):
            return (jnp_fn if any(_is_tracer_arg(a) for a in args) else np_fn)(*args, **kw)
        return dispatch


_xp = _TracerXP()
