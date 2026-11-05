"""Check that port_fortran_declarations_to_jax_docstrings.py copies Fortran
declarations correctly.

Usage:
    human_reviewed_tests/run_pytests.sh -unit -k test_port_fortran_declarations_to_jax_docstrings

Understanding:
    The test writes a small Fortran module and JAX module to a temporary
    directory, so it does not depend on the repository sources. Whether the
    repository's JAX docstrings match the Fortran is checked by the
    "Fortran declaration match" stage of the clubb_jax Jenkins job.

    It covers the parsing cases that broke during development: "type( grid )"
    with a space, a comment line inside a continued declaration, a function
    result without intent, local variables that must be left out (including
    an array whose bounds name an argument), and a generic interface whose
    specific is chosen by argument names. It checks the written section,
    that -check then passes, that a second -write_declaration changes
    nothing, and that only the docstring changed.
"""

import ast

from utilities.port_fortran_declarations_to_jax_docstrings import \
    port_fortran_declarations_to_jax_docstrings as sync


FORTRAN = """\
module demo

  implicit none

  interface scale
    module procedure scale_1D, scale_2D
  end interface

  contains

  function scale_1D( nz, x ) result( y )

    ! Description:
    ! Scales one column.
    implicit none

    integer, intent(in) :: nz
    real, dimension(nz), intent(in) :: x
    real, dimension(nz) :: y

    y = 2.0 * x
  end function scale_1D

  function scale_2D( nz, ngrdcol, gr, x ) result( y )

    ! Description:
    ! Scales every column.

    ! References:
    ! None
    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nz, &      ! Number of levels
      ! Columns are the leading dimension.
      ngrdcol    ! Number of columns

    type( grid ), intent(in) :: &
      gr         ! Grid

    real, dimension(ngrdcol,nz), intent(in) :: &
      x          ! Input field [-]

    !--------------------------- Output Variables ---------------------------
    real, dimension(ngrdcol,nz) :: &
      y          ! Scaled field [-]

    !--------------------------- Local Variables ---------------------------
    integer :: k
    real :: work(ngrdcol,nz)

    y = 2.0 * x
  end function scale_2D

end module demo
"""

JAX = '''\
def scale(nz, ngrdcol, gr, x):
    """Scale x."""
    # Doubles the field.
    return 2.0 * x


def _helper(x):
    return x
'''

EXPECTED_SECTION = [
    "Fortran declarations (scale_2D in {path}):",
    "",
    "    !--------------------------- Input Variables ---------------------------",
    "    integer, intent(in) :: &",
    "      nz, &      ! Number of levels",
    "      ! Columns are the leading dimension.",
    "      ngrdcol    ! Number of columns",
    "",
    "    type( grid ), intent(in) :: &",
    "      gr         ! Grid",
    "",
    "    real, dimension(ngrdcol,nz), intent(in) :: &",
    "      x          ! Input field [-]",
    "",
    "    !--------------------------- Output Variables ---------------------------",
    "    real, dimension(ngrdcol,nz) :: &",
    "      y          ! Scaled field [-]",
]


def test_write_copies_declarations_of_matching_specific(tmp_path):
    fortran_dir, jax_dir = tmp_path / "fortran", tmp_path / "jax"
    fortran_dir.mkdir()
    jax_dir.mkdir()
    (fortran_dir / "demo.F90").write_text(FORTRAN)
    jax_file = jax_dir / "demo.py"
    jax_file.write_text(JAX)
    options = ["-fortran_dir", str(fortran_dir), "-jax_dir", str(jax_dir)]

    assert sync.main(options) == 1
    assert sync.main(["-write_declaration"] + options) == 0
    written = jax_file.read_text()

    function = ast.parse(written).body[0]
    docstring = ast.get_docstring(function)
    section = sync.current_section(docstring)
    assert section == [line.format(path=fortran_dir / "demo.F90") for line in EXPECTED_SECTION]
    assert docstring.startswith("Scale x.\n\nFortran declarations (")

    assert sync.main(options) == 0
    assert sync.main(["-write_declaration"] + options) == 0
    assert jax_file.read_text() == written
    assert sync.strip_docstrings(ast.parse(JAX)) == sync.strip_docstrings(ast.parse(written))
