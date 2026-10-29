"""Test wrappers for routines from CLUBB_core/file_functions.F90."""
from utilities.output_paths import REPO_ROOT as _REPO_ROOT
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(_REPO_ROOT / "clubb_python_api"))

from clubb_python import clubb_api


def test_file_read_1d_reads_flat_rows(tmp_path):
    """file_read_1d should read values row-wise into a flat vector."""
    path = tmp_path / "sample_1d.dat"
    path.write_text("1 2 3 4\n5 6\n")

    out = clubb_api.file_read_1d(
        file_unit=99,
        path_and_filename=str(path),
        num_datapts=6,
        entries_per_line=4,
    )

    np.testing.assert_allclose(out, np.array([1, 2, 3, 4, 5, 6], dtype=np.float64))
