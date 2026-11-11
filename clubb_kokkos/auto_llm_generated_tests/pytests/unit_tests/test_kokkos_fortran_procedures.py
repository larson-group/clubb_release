"""Check Native declarations read by the permanent namelist generators.

The real declaration reader must accept Fortran's bare subroutine form and
still reject a malformed header rather than treating END SUBROUTINE as one.
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "clubb_kokkos/src/f2k"))
from fortran_procedures import read_procedure


def test_kokkos_bare_zero_argument_metadata_procedure():
    """Read the actual error-header initializer without adding fake arguments."""
    procedure = read_procedure(
        ROOT / "src/CLUBB_core/error_code.F90", "initialize_error_headers"
    )
    assert procedure.arguments == []
    assert procedure.result is None


def test_kokkos_end_subroutine_is_not_a_header(tmp_path):
    """A malformed start cannot silently become an empty contract at its END."""
    source = tmp_path / "invalid.F90"
    source.write_text(
        "module invalid\ncontains\nsubroutine broken(argument\n"
        "end subroutine broken\nend module invalid\n"
    )
    with pytest.raises(ValueError, match="missing procedure"):
        read_procedure(source, "broken")
