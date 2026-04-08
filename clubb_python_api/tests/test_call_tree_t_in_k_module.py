"""Test wrappers for T_in_K_module call-tree routines."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_thlm2t_in_k_returns_finite_physical_values():
    """thlm2t_in_k should return a finite positive temperature profile."""
    thlm = np.linspace(290.0, 310.0, 20, dtype=np.float64)
    exner = np.linspace(1.0, 0.85, 20, dtype=np.float64)
    rcm = np.linspace(0.002, 0.0, 20, dtype=np.float64)

    result = clubb_api.thlm2t_in_k(
        nz=thlm.shape[0],
        thlm=thlm,
        exner=exner,
        rcm=rcm,
    )

    assert result.shape == thlm.shape
    assert np.all(np.isfinite(result))
    assert np.all(result > 0.0)


def test_thlm2t_in_k_increases_with_thlm():
    """Holding exner and condensate fixed, warmer thlm should yield warmer T."""
    exner = np.full(8, 0.95, dtype=np.float64)
    rcm = np.full(8, 0.001, dtype=np.float64)
    cooler_thlm = np.full(8, 295.0, dtype=np.float64)
    warmer_thlm = np.full(8, 300.0, dtype=np.float64)

    cooler = clubb_api.thlm2t_in_k(
        nz=cooler_thlm.shape[0],
        thlm=cooler_thlm,
        exner=exner,
        rcm=rcm,
    )
    warmer = clubb_api.thlm2t_in_k(
        nz=warmer_thlm.shape[0],
        thlm=warmer_thlm,
        exner=exner,
        rcm=rcm,
    )

    assert np.all(warmer > cooler)


def test_thlm2t_in_k_increases_with_rcm():
    """Holding thlm and exner fixed, more condensate should raise T."""
    thlm = np.full(8, 300.0, dtype=np.float64)
    exner = np.full(8, 0.95, dtype=np.float64)
    drier_rcm = np.zeros(8, dtype=np.float64)
    wetter_rcm = np.full(8, 0.002, dtype=np.float64)

    drier = clubb_api.thlm2t_in_k(
        nz=thlm.shape[0],
        thlm=thlm,
        exner=exner,
        rcm=drier_rcm,
    )
    wetter = clubb_api.thlm2t_in_k(
        nz=thlm.shape[0],
        thlm=thlm,
        exner=exner,
        rcm=wetter_rcm,
    )

    assert np.all(wetter > drier)
