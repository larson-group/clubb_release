"""Require the compiled Fortran reference for the JAX F2PY unit tests.

Run these tests with `bash human_reviewed_tests/run_pytests.sh -jax_f2py_unit_tests`,
which puts a compatible F2PY build on the import path.
A missing or unloadable build is a setup error rather than a skip: silently
skipped Fortran comparisons are how the previous generation of these checks
stopped comparing anything.
"""

import sys

import jax

jax.config.update("jax_enable_x64", True)

try:
    # The package import preloads a compatible NetCDF stack before clubb_f2py.
    import clubb_python  # noqa: F401
    import clubb_f2py  # noqa: F401
except ImportError as error:
    raise ImportError(
        "The JAX F2PY unit tests need the compiled CLUBB Python API. Build it with "
        "`./compile.py -debug -python` and run `bash human_reviewed_tests/run_pytests.sh -jax_f2py_unit_tests` "
        "(set CLUBB_F2PY_DIR to use a build other than install/latest/python). "
        "The build must use the same Python minor version as the JAX environment "
        f"(Python {sys.version_info.major}.{sys.version_info.minor} here): "
        "a clubb_f2py.cpython-3XX-*.so built by another Python cannot be imported. "
        "Compare .venv-python/bin/python -V with .venv-jax/bin/python -V, or set "
        "PYTHON to one interpreter when creating both environments."
    ) from error
