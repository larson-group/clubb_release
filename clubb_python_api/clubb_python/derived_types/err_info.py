"""Python representation of err_info_type derived-type data."""

from typing import NamedTuple, Optional

from clubb_python.derived_types.common import Array


class ErrInfo(NamedTuple):
    """Subset of err_info_type that is currently bridgeable via F2PY wrappers."""

    ngrdcol: int
    chunk_idx: int = 1
    mpi_rank: int = 0
    lat: Optional[Array] = None
    lon: Optional[Array] = None
    err_code: Optional[Array] = None
