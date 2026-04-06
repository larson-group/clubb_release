"""Test wrappers for routines from CLUBB_core/corr_varnce_module.F90."""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def test_assert_corr_symmetric_accepts_identity_and_flags_bad_matrix(run_quiet_python):
    """assert_corr_symmetric should pass identity and flag nonsymmetric inputs."""
    ngrdcol = 1
    clubb_api.init_err_info(ngrdcol)
    err_info = ErrInfo(ngrdcol=ngrdcol)

    clubb_api.reset_err_code()
    clubb_api.assert_corr_symmetric(3, np.eye(3, dtype=np.float64), err_info)
    assert int(clubb_api.get_err_code(ngrdcol)[0]) == 0

    clubb_api.reset_err_code()
    bad = np.array(
        [[1.0, 0.2, 0.0],
         [0.1, 1.0, 0.0],
         [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    proc = run_quiet_python(
        """
import json
import numpy as np
from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo

ngrdcol = 1
clubb_api.init_err_info(ngrdcol)
err_info = ErrInfo(ngrdcol=ngrdcol)
clubb_api.reset_err_code()
bad = np.array(
    [[1.0, 0.2, 0.0],
     [0.1, 1.0, 0.0],
     [0.0, 0.0, 1.0]],
    dtype=np.float64,
)
clubb_api.assert_corr_symmetric(3, bad, err_info)
print(json.dumps({"err_code": int(clubb_api.get_err_code(ngrdcol)[0])}))
"""
    )
    assert proc.returncode == 0, proc.stderr
    assert json.loads(proc.stdout)["err_code"] != 0
