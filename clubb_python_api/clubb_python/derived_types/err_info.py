"""Python representation of err_info_type derived-type data."""

from typing import NamedTuple, Optional

import numpy as np

from clubb_python.derived_types.common import Array


class ErrInfo(NamedTuple):
    """Subset of err_info_type that is currently bridgeable via F2PY wrappers."""

    ngrdcol: int
    chunk_idx: int = 1
    mpi_rank: int = 0
    lat: Optional[Array] = None
    lon: Optional[Array] = None
    err_code: Optional[Array] = None

    def is_fatal(self) -> bool:
        """Return whether any column has a fatal CLUBB error code."""
        if self.err_code is None:
            return False
        from clubb_python.CLUBB_core.error_code import CLUBB_FATAL_ERROR

        return bool(np.any(np.asarray(self.err_code) == CLUBB_FATAL_ERROR))

    def has_fatal(self) -> bool:
        """Alias for is_fatal."""
        return self.is_fatal()

    def set_fatal(self, mask=None):
        """Return a copy with fatal error codes set.

        If mask is omitted, all columns are marked fatal. If mask is provided,
        only true mask entries are marked fatal.
        """
        from clubb_python.CLUBB_core.error_code import (
            CLUBB_FATAL_ERROR,
            CLUBB_NO_ERROR,
        )

        err_code = (
            np.full(int(self.ngrdcol), CLUBB_NO_ERROR, dtype=np.int32)
            if self.err_code is None
            else np.asarray(self.err_code, dtype=np.int32).copy()
        )
        if mask is None:
            err_code[:] = CLUBB_FATAL_ERROR
        else:
            mask_arr = np.asarray(mask, dtype=bool).reshape(-1)
            if mask_arr.size != int(self.ngrdcol):
                raise ValueError(
                    f"fatal mask length mismatch: expected {self.ngrdcol}, "
                    f"got {mask_arr.size}"
                )
            err_code[mask_arr] = CLUBB_FATAL_ERROR
        return self._replace(err_code=err_code)

    def mark_fatal(self, mask=None):
        """Alias for set_fatal."""
        return self.set_fatal(mask=mask)

    def reset_code(self):
        """Return a copy with all error codes reset to no-error."""
        from clubb_python.CLUBB_core.error_code import CLUBB_NO_ERROR

        return self._replace(
            err_code=np.full(int(self.ngrdcol), CLUBB_NO_ERROR, dtype=np.int32)
        )
