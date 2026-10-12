"""Python representation of nu_vert_res_dep type (parameters_tunable.F90)."""

from typing import NamedTuple

import numpy as np

from clubb_python.derived_types.common import Array


class NuVertResDep(NamedTuple):
    """Vertical-resolution-dependent nu coefficient arrays."""

    nzm: int
    nu1: Array
    nu2: Array
    nu6: Array
    nu8: Array
    nu9: Array
    nu10: Array
    nu_hm: Array

    @classmethod
    def zeros(cls, ngrdcol: int) -> "NuVertResDep":
        """Zero coefficients for ngrdcol columns, e.g. to pass to calc_derived_params."""
        ngrdcol = int(ngrdcol)
        return cls(ngrdcol, *(np.zeros(ngrdcol, dtype=np.float64) for _ in range(7)))
