"""Python representation of scalar indices from array_index.F90."""

from typing import NamedTuple


class SclrIdx(NamedTuple):
    """Scalar and extra diagnostic scalar index mapping."""

    iisclr_rt: int
    iisclr_thl: int
    iisclr_CO2: int
    iiedsclr_rt: int
    iiedsclr_thl: int
    iiedsclr_CO2: int
