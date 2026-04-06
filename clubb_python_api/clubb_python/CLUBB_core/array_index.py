"""User-facing wrappers for scalar-index data using CLUBB_core/array_index.F90 conventions."""

from clubb_python.derived_types.sclr_idx import SclrIdx
from clubb_python.derived_types.sclr_idx_converter import get_fortran_sclr_idx, set_fortran_sclr_idx


def set_sclr_idx(
    iisclr_rt: int, iisclr_thl: int, iisclr_co2: int = 0,
    iiedsclr_rt: int = 0, iiedsclr_thl: int = 0, iiedsclr_co2: int = 0,
):
    """Set scalar index mapping in Fortran module storage."""
    set_fortran_sclr_idx(
        SclrIdx(
            iisclr_rt=int(iisclr_rt),
            iisclr_thl=int(iisclr_thl),
            iisclr_CO2=int(iisclr_co2),
            iiedsclr_rt=int(iiedsclr_rt),
            iiedsclr_thl=int(iiedsclr_thl),
            iiedsclr_CO2=int(iiedsclr_co2),
        )
    )


def get_sclr_idx():
    """Get scalar index mapping from Fortran module storage."""
    idx = get_fortran_sclr_idx()
    return tuple(int(v) for v in idx)
