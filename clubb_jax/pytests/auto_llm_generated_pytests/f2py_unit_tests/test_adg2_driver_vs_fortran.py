"""Check the JAX ADG2_pdf_driver port (adg1_adg2_3d_luhar_pdf.F90)
against Fortran."""

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.adg1_adg2_3d_luhar_pdf as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import ADG2_pdf_driver as _ADG2_pdf_driver

NG, NZ = 2, 8
_BETA = 1.75


def ADG2_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, wprtp, wpthlp, sqrt_wp2, beta):
    ng, nz = np.asarray(wm).shape
    sclr_dim = 1
    sclr_tol = jnp.full((sclr_dim,), 1e-8)
    sclrm = jnp.zeros((ng, nz, sclr_dim), dtype=jnp.float64)
    sclrp2 = jnp.ones((ng, nz, sclr_dim), dtype=jnp.float64) * 1e-6
    wpsclrp = jnp.zeros((ng, nz, sclr_dim), dtype=jnp.float64)
    return _ADG2_pdf_driver(
        nz,
        ng,
        sclr_dim,
        sclr_tol,
        wm,
        rtm,
        thlm,
        wp2,
        rtp2,
        thlp2,
        Skw,
        wprtp,
        wpthlp,
        sqrt_wp2,
        beta,
        sclrm,
        sclrp2,
        wpsclrp,
        False,
    )


def _inputs(rng):
    sh = (NG, NZ)
    wp2 = rng.uniform(0.05, 2, sh)
    return dict(wm=rng.uniform(-1, 1, sh), rtm=rng.uniform(0, 0.02, sh), thlm=rng.uniform(285, 305, sh),
                wp2=wp2, rtp2=rng.uniform(1e-7, 1e-6, sh), thlp2=rng.uniform(0.05, 1, sh),
                Skw=rng.uniform(-2, 2, sh), wprtp=rng.uniform(-1e-3, 1e-3, sh), wpthlp=rng.uniform(-0.5, 0.5, sh),
                sqrt_wp2=np.sqrt(wp2))


def test_f2py_oracle():
    rng = np.random.default_rng(11)
    worst = 0.0
    for _ in range(3):
        a = _inputs(rng)
        beta_col = np.full(NG, _BETA)
        # f2py rejects 0-sized dims, so use sclr_dim=1 with dummy scalars and l_scalar_calc=False
        # (the scalar branch is then skipped and only the 16 w/rt/thl outputs are compared).
        sd = 1
        sclr_tol = np.full(sd, 1e-8)
        sm = np.zeros((NG, NZ, sd)); sp = np.ones((NG, NZ, sd)) * 1e-6; ws = np.zeros((NG, NZ, sd))
        f = fortran.adg2_pdf_driver(
            NZ, NG, sd, sclr_tol, a['wm'], a['rtm'], a['thlm'], a['wp2'], a['rtp2'], a['thlp2'],
            a['Skw'], a['wprtp'], a['wpthlp'], a['sqrt_wp2'], beta_col, sm, sp, ws, False, None)
        g = ADG2_pdf_driver(a['wm'], a['rtm'], a['thlm'], a['wp2'], a['rtp2'], a['thlp2'],
                            a['Skw'], a['wprtp'], a['wpthlp'], a['sqrt_wp2'], beta_col[:, None])
        for k in range(16):   # the 16 non-scalar PDF-param outputs
            worst = max(worst, max_abs_diff("adg2_driver", np.asarray(g[k]), np.asarray(f[k])))
    assert worst < 1e-10, f"adg2_pdf_driver f2py mismatch {worst:.2e}"
