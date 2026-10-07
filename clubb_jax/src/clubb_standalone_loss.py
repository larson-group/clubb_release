"""Command-line JAX loss evaluation using the shared Fortran-style table.

Adapted from src/clubb_standalone_loss.F90. The existing Python CLI requires
one namelist argument and accepts help instead of using the native clubb.in
default. Runtime setup belongs to clubb_jax/run_jax.py.
"""

import sys

from clubb_jax.src.clubb_loss_driver import clubb_get_loss
from utilities.loss_metrics import print_loss_matrix


def main():
    """Print the shared loss table and retain the native completion status."""
    if len(sys.argv) != 2 or sys.argv[1] in ("-h", "-help", "--help"):
        print("Usage: python -m clubb_jax.src.clubb_standalone_loss <namelist_path>")
        return 0 if any(arg in ("-h", "-help", "--help") for arg in sys.argv[1:]) else 1

    print_loss_matrix(*clubb_get_loss(sys.argv[1]))
    # run_scm_loss/run_case uses 6 for a completed loss evaluation.
    return 6


if __name__ == "__main__":
    sys.exit(main())
