"""User-facing wrappers for routines from CLUBB_core/file_functions.F90."""

import clubb_f2py


def file_read_1d(file_unit: int, path_and_filename: str, num_datapts: int, entries_per_line: int):
    """Read a 1D data file with the CLUBB file_functions layout."""
    return clubb_f2py.f2py_file_read_1d(
        int(file_unit), str(path_and_filename), int(num_datapts), int(entries_per_line)
    )
