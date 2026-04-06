"""Push/pull converters for SclrIdx <-> Fortran module storage."""

import clubb_f2py

from clubb_python.derived_types.sclr_idx import SclrIdx


def set_fortran_sclr_idx(idx: SclrIdx):
    """Push scalar indices into Fortran module storage."""
    clubb_f2py.set_sclr_idx(*idx)


def get_fortran_sclr_idx() -> SclrIdx:
    """Pull scalar indices from Fortran module storage."""
    return SclrIdx(*clubb_f2py.get_sclr_idx())
