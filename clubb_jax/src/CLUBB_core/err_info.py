"""JAX-side representation of `err_info_type`."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.error_code import CLUBB_FATAL_ERROR, CLUBB_NO_ERROR
from clubb_jax.src.CLUBB_core.err_info_codes import ERR_NONE, messages_for


_ERR_INFO_STATIC_FIELDS = (
    "ngrdcol",
    "chunk_idx",
    "mpi_rank",
)

_ERR_INFO_ARRAY_FIELDS = (
    "lat",
    "lon",
    "err_code",
    "reason_code",
)


@jax.tree_util.register_pytree_node_class
class ErrInfo(NamedTuple):
    """Bridgeable subset of CLUBB `err_info_type` with functional updates."""

    ngrdcol: int
    chunk_idx: int = 1
    mpi_rank: int = 0
    lat: jnp.ndarray | None = None
    lon: jnp.ndarray | None = None
    err_code: jnp.ndarray | None = None
    reason_code: jnp.ndarray | None = None

    @classmethod
    def initialized(
        cls,
        ngrdcol: int,
        chunk_idx: int = 1,
        mpi_rank: int = 0,
        lat=None,
        lon=None,
    ):
        return cls(
            ngrdcol=ngrdcol,
            chunk_idx=chunk_idx,
            mpi_rank=mpi_rank,
            lat=lat,
            lon=lon,
            err_code=jnp.full((int(ngrdcol),), CLUBB_NO_ERROR, dtype=jnp.int32),
            reason_code=jnp.full((int(ngrdcol),), ERR_NONE, dtype=jnp.int32),
        )

    def replace(self, **kwargs):
        return self._replace(**kwargs)

    def tree_flatten(self):
        children = tuple(getattr(self, name) for name in _ERR_INFO_ARRAY_FIELDS)
        aux_data = tuple(getattr(self, name) for name in _ERR_INFO_STATIC_FIELDS)
        return children, aux_data

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        data = dict(zip(_ERR_INFO_STATIC_FIELDS, aux_data))
        data.update(dict(zip(_ERR_INFO_ARRAY_FIELDS, children)))
        return cls(**data)

    def err_code_or_default(self):
        if self.err_code is None:
            return jnp.full((int(self.ngrdcol),), CLUBB_NO_ERROR, dtype=jnp.int32)
        return jnp.asarray(self.err_code, dtype=jnp.int32)

    def reason_code_or_default(self):
        if self.reason_code is None:
            return jnp.full((int(self.ngrdcol),), ERR_NONE, dtype=jnp.int32)
        return jnp.asarray(self.reason_code, dtype=jnp.int32)

    def fatal_mask(self):
        return self.err_code_or_default() == CLUBB_FATAL_ERROR

    def any_fatal(self):
        return jnp.any(self.fatal_mask())

    def has_fatal(self):
        return self.any_fatal()

    def is_fatal(self) -> bool:
        return self.is_fatal_host()

    def is_fatal_host(self) -> bool:
        return bool(self.any_fatal())

    def has_fatal_host(self) -> bool:
        return self.is_fatal_host()

    def set_fatal(self, mask=None):
        err_code = self.err_code_or_default()

        if mask is None:
            err_code = jnp.full_like(err_code, CLUBB_FATAL_ERROR)
        else:
            mask = jnp.asarray(mask, dtype=bool)
            if mask.ndim == 0:
                mask = jnp.broadcast_to(mask, err_code.shape)
            else:
                mask = jnp.reshape(mask, err_code.shape)
            err_code = jnp.where(mask, CLUBB_FATAL_ERROR, err_code)

        return self._replace(err_code=err_code)

    def mark_fatal(self, mask=None):
        return self.set_fatal(mask=mask)

    def set_reason(self, reason_code, mask=None):
        current_reason = self.reason_code_or_default()
        new_reason = jnp.asarray(reason_code, dtype=jnp.int32)

        if mask is None:
            mask = current_reason == ERR_NONE
        else:
            mask = jnp.asarray(mask, dtype=bool)
            if mask.ndim == 0:
                mask = jnp.broadcast_to(mask, current_reason.shape)
            else:
                mask = jnp.reshape(mask, current_reason.shape)
            mask = mask & (current_reason == ERR_NONE)

        reason = jnp.where(mask, new_reason, current_reason)
        return self._replace(reason_code=reason)

    def reason_messages_host(self) -> list[str]:
        return messages_for(self.reason_code_or_default())

    def reset_code(self):
        return self._replace(
            err_code=jnp.full((int(self.ngrdcol),), CLUBB_NO_ERROR, dtype=jnp.int32),
            reason_code=jnp.full((int(self.ngrdcol),), ERR_NONE, dtype=jnp.int32),
        )
