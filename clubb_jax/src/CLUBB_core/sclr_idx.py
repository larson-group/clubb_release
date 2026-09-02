"""JAX-side representation of `sclr_idx_type`."""

from __future__ import annotations

from typing import NamedTuple

import jax


_SCLR_IDX_STATIC_FIELDS = (
    "iisclr_rt",
    "iisclr_thl",
    "iisclr_CO2",
    "iiedsclr_rt",
    "iiedsclr_thl",
    "iiedsclr_CO2",
)


@jax.tree_util.register_pytree_node_class
class SclrIdx(NamedTuple):
    """Scalar and extra diagnostic scalar index mapping."""

    iisclr_rt: int
    iisclr_thl: int
    iisclr_CO2: int
    iiedsclr_rt: int
    iiedsclr_thl: int
    iiedsclr_CO2: int

    def replace(self, **kwargs):
        return self._replace(**kwargs)

    def tree_flatten(self):
        aux_data = tuple(getattr(self, name) for name in _SCLR_IDX_STATIC_FIELDS)
        return (), aux_data

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        del children
        return cls(**dict(zip(_SCLR_IDX_STATIC_FIELDS, aux_data)))
