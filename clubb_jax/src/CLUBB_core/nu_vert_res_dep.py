"""JAX-side representation of `nu_vert_res_dep`."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp


_NU_STATIC_FIELDS = ("nzm",)

_NU_ARRAY_FIELDS = (
    "nu1",
    "nu2",
    "nu6",
    "nu8",
    "nu9",
    "nu10",
    "nu_hm",
)


@jax.tree_util.register_pytree_node_class
class NuVertResDep(NamedTuple):
    """Vertical-resolution-dependent nu coefficient arrays."""

    nzm: int
    nu1: jnp.ndarray
    nu2: jnp.ndarray
    nu6: jnp.ndarray
    nu8: jnp.ndarray
    nu9: jnp.ndarray
    nu10: jnp.ndarray
    nu_hm: jnp.ndarray

    def replace(self, **kwargs):
        return self._replace(**kwargs)

    def tree_flatten(self):
        children = tuple(getattr(self, name) for name in _NU_ARRAY_FIELDS)
        aux_data = tuple(getattr(self, name) for name in _NU_STATIC_FIELDS)
        return children, aux_data

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        data = dict(zip(_NU_STATIC_FIELDS, aux_data))
        data.update(dict(zip(_NU_ARRAY_FIELDS, children)))
        return cls(**data)
