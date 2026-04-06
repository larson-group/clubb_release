"""Test wrappers for routines from CLUBB_core/interpolation.F90."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_mono_cubic_interp_linear_profile():
    """A linear profile should interpolate to the linear midpoint value."""
    out = clubb_api.mono_cubic_interp(
        z_in=1.5,
        km1=0,
        k00=1,
        kp1=2,
        kp2=3,
        zm1=0.0,
        z00=1.0,
        zp1=2.0,
        zp2=3.0,
        fm1=0.0,
        f00=1.0,
        fp1=2.0,
        fp2=3.0,
    )
    np.testing.assert_allclose(out, 1.5)


def test_zlinterp_fnc_zero_outside_domain():
    """zlinterp_fnc should interpolate inside the domain and zero outside it."""
    grid_src = np.array([0.0, 1.0, 2.0], dtype=np.float64)
    var_src = np.array([0.0, 10.0, 20.0], dtype=np.float64)
    grid_out = np.array([-1.0, 0.5, 2.5], dtype=np.float64)

    out = clubb_api.zlinterp_fnc(
        dim_out=grid_out.size,
        dim_src=grid_src.size,
        grid_out=grid_out,
        grid_src=grid_src,
        var_src=var_src,
    )

    np.testing.assert_allclose(out, np.array([0.0, 5.0, 0.0], dtype=np.float64))
