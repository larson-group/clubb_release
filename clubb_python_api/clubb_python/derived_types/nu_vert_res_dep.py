"""Python representation of nu_vert_res_dep type (parameters_tunable.F90)."""

from typing import NamedTuple

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
