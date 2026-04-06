"""Test wrappers for saturation call-tree routines."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_sat_mixrat_ice_returns_finite_positive_values():
    """sat_mixrat_ice should produce finite positive values below freezing."""
    p_in_pa = np.array([[90000.0, 80000.0, 70000.0]], dtype=np.float64)
    t_in_k = np.array([[260.0, 250.0, 240.0]], dtype=np.float64)

    rsat = clubb_api.sat_mixrat_ice(
        nz=t_in_k.shape[1],
        ngrdcol=t_in_k.shape[0],
        p_in_Pa=p_in_pa,
        T_in_K=t_in_k,
        saturation_formula=3,
    )

    assert rsat.shape == t_in_k.shape
    assert np.all(np.isfinite(rsat))
    assert np.all(rsat > 0.0)
