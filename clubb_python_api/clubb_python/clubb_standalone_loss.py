"""Command-line front end to the Fortran loss API.

Adapted from src/clubb_standalone_loss.F90. The existing Python CLI requires
one namelist argument and accepts help instead of using the native clubb.in
default. Shared table formatting preserves the existing output contract.
"""

import sys

from clubb_python import clubb_api
from utilities.loss_metrics import print_loss_matrix


def main():
    """Evaluate one namelist and return the native completion status."""
    if len(sys.argv) != 2 or sys.argv[1] in ("-h", "-help", "--help"):
        print("Usage: python -m clubb_python.clubb_standalone_loss <namelist_path>")
        return 0 if any(arg in ("-h", "-help", "--help") for arg in sys.argv[1:]) else 1

    print_loss_matrix(*clubb_api.clubb_get_loss(sys.argv[1]))
    # run_scm_loss/run_case uses 6 for a completed loss evaluation.
    return 6


if __name__ == "__main__":
    sys.exit(main())
