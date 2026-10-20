"""Check the JAX Luhar_3D_pdf_driver + backsolve_Luhar_params (adg1_adg2_3d_luhar_pdf.F90)
against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.adg1_adg2_3d_luhar_pdf as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import Luhar_3D_pdf_driver

NG, NZ = 2, 8


def _inputs(rng):
    sh = (NG, NZ)
    return dict(wm=rng.uniform(-1, 1, sh), rtm=rng.uniform(0, 0.02, sh), thlm=rng.uniform(285, 305, sh),
                wp2=rng.uniform(0.05, 2, sh), rtp2=rng.uniform(1e-7, 1e-6, sh), thlp2=rng.uniform(0.05, 1, sh),
                Skw=rng.uniform(-2, 2, sh), Skrt=rng.uniform(-2, 2, sh), Skthl=rng.uniform(-2, 2, sh),
                wprtp=rng.uniform(-1e-3, 1e-3, sh), wpthlp=rng.uniform(-0.5, 0.5, sh))


def test_f2py_oracle():
    rng = np.random.default_rng(11)
    order = ('wm', 'rtm', 'thlm', 'wp2', 'rtp2', 'thlp2', 'Skw', 'Skrt', 'Skthl', 'wprtp', 'wpthlp')
    worst = 0.0
    for _ in range(4):
        a = _inputs(rng)
        f = fortran.luhar_3d_pdf_driver(NZ, *[a[k] for k in order])
        g = Luhar_3D_pdf_driver(*[a[k] for k in order])
        worst = max(worst, max_abs_diff("Luhar_3D_pdf_driver", tuple(g), tuple(f)))
    assert worst < 1e-9, f"luhar_3d_pdf_driver Fortran mismatch {worst:.2e}"
