"""Check the JAX close_Luhar_pdf port (adg1_adg2_3d_luhar_pdf.F90)
against Fortran."""

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.adg1_adg2_3d_luhar_pdf as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import close_Luhar_pdf, ADG1_w_closure

NG, NZ = 2, 6
_X_TOL_SQD = 1.0e-8


def _inputs(seed):
    rng = np.random.default_rng(seed)
    xm = rng.uniform(-2, 2, (NG, NZ))
    xp2 = rng.uniform(0.05, 2.0, (NG, NZ))      # all > x_tol_sqd
    mixt_frac = rng.uniform(0.2, 0.8, (NG, NZ))
    small_m = rng.uniform(0.05, 1.0, (NG, NZ))
    wpxp = rng.uniform(-1, 1, (NG, NZ))
    return xm, xp2, mixt_frac, small_m, wpxp


def test_f2py_oracle():
    worst = 0.0
    for seed in (11, 22, 33):
        xm, xp2, mf, m, wpxp = _inputs(seed)
        f = fortran.close_luhar_pdf(NZ, xm, xp2, mf, m, wpxp, _X_TOL_SQD)
        g = close_Luhar_pdf(xm, xp2, mf, m, wpxp, _X_TOL_SQD)
        assert len(f) == len(g)
        for fi, gi in zip(f, g):
            worst = max(worst, max_abs_diff("close_luhar_pdf", gi, fi))
    assert worst < 1e-11, f"close_luhar_pdf f2py mismatch {worst:.2e}"


def test_ADG1_w_closure_f2py():
    """ADG1_w_closure (adg1_adg2_3d_luhar_pdf.F90) — the ADG1 two-component w-PDF closure (w_1/w_2, normal-space
    w_1_n/w_2_n, varnce_w_1/2, mixt_frac from wm/wp2/Skw/sigma_sqd_w). A core every-step routine that had only a
    grad test; bit-shadow its 7 outputs vs the f2py oracle."""
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(40):
        ng, nz = 2, 8
        wp2 = rng.uniform(1e-3, 2.0, (ng, nz)); wm = rng.uniform(-1.0, 1.0, (ng, nz))
        Skw = rng.uniform(-3.0, 3.0, (ng, nz)); ssw = rng.uniform(0.0, 0.6, (ng, nz)); mfmm = 0.4
        jr = ADG1_w_closure(jnp.asarray(wm), jnp.asarray(wp2), jnp.asarray(Skw), jnp.asarray(ssw),
                            jnp.asarray(np.sqrt(wp2)), mfmm)
        fr = fortran.adg1_w_closure(nz, ng, wm, wp2, Skw, ssw, np.sqrt(wp2), mfmm, None)
        for i in range(7):
            worst = max(worst, max_abs_diff(f"ADG1_w_closure[{i}]", jr[i], fr[i]))
    assert worst < 1e-12, f"ADG1_w_closure f2py mismatch {worst:.2e}"


def test_ADG1_pdf_driver_f2py():
    """ADG1_pdf_driver (adg1_adg2_3d_luhar_pdf.F90) — THE core every-step ADG1 PDF-parameter driver: from the mean
    fields + 2nd moments it produces the two-component (w, rt, thl, u, v) means/variances + mixt_frac + alpha
    factors (25 non-scalar outputs). Case-validated end-to-end but had no unit test. Bit-shadow all 25 outputs vs
    the f2py oracle (sclr_dim=0, matching the JAX driver). The JAX returns a dict; f2py outputs 0..24 are these names in order."""
    from clubb_jax.src.CLUBB_core.adg1_adg2_3d_luhar_pdf import ADG1_pdf_driver
    fn = ['w_1', 'w_2', 'rt_1', 'rt_2', 'thl_1', 'thl_2', 'u_1', 'u_2', 'v_1', 'v_2',
          'varnce_w_1', 'varnce_w_2', 'varnce_rt_1', 'varnce_rt_2', 'varnce_thl_1', 'varnce_thl_2',
          'varnce_u_1', 'varnce_u_2', 'varnce_v_1', 'varnce_v_2', 'mixt_frac', 'alpha_rt', 'alpha_thl',
          'alpha_u', 'alpha_v']
    rng = np.random.default_rng(1); ng, nz = 2, 8; worst = 0.0
    for _ in range(20):
        wm = rng.uniform(-1, 1, (ng, nz)); rtm = rng.uniform(1e-3, 1.5e-2, (ng, nz)); thlm = rng.uniform(290, 320, (ng, nz))
        um = rng.uniform(-5, 5, (ng, nz)); vm = rng.uniform(-5, 5, (ng, nz)); wp2 = rng.uniform(1e-2, 2, (ng, nz))
        rtp2 = rng.uniform(1e-7, 1e-5, (ng, nz)); thlp2 = rng.uniform(1e-3, 1, (ng, nz))
        up2 = rng.uniform(1e-2, 1, (ng, nz)); vp2 = rng.uniform(1e-2, 1, (ng, nz))
        Skw = rng.uniform(-3, 3, (ng, nz)); wprtp = rng.uniform(-1e-3, 1e-3, (ng, nz)); wpthlp = rng.uniform(-1, 1, (ng, nz))
        upwp = rng.uniform(-0.1, 0.1, (ng, nz)); vpwp = rng.uniform(-0.1, 0.1, (ng, nz))
        sqw = np.sqrt(wp2); ssw = rng.uniform(0, 0.5, (ng, nz)); beta = np.full(ng, 1.0); mfmm = 0.4
        jr = ADG1_pdf_driver(jnp.asarray(wm), jnp.asarray(rtm), jnp.asarray(thlm), jnp.asarray(um), jnp.asarray(vm),
            jnp.asarray(wp2), jnp.asarray(rtp2), jnp.asarray(thlp2), jnp.asarray(up2), jnp.asarray(vp2),
            jnp.asarray(Skw), jnp.asarray(wprtp), jnp.asarray(wpthlp), jnp.asarray(upwp), jnp.asarray(vpwp),
            jnp.asarray(sqw), jnp.asarray(ssw), jnp.asarray(beta)[:, None], mfmm)
        # CLUBB sizes scalar arrays with max(sclr_dim, 1); the JAX driver has no scalars.
        sm = np.zeros((ng, nz, 1))
        fr = fortran.adg1_pdf_driver(nz, ng, 0, np.zeros(1), wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2,
                                     Skw, wprtp, wpthlp, upwp, vpwp, sqw, ssw, beta, mfmm, sm, sm, sm, False, None)
        for i, nm in enumerate(fn):
            worst = max(worst, max_abs_diff(f"ADG1_pdf_driver {nm}", jr[nm], fr[i]))
    assert worst < 1e-12, f"ADG1_pdf_driver f2py mismatch {worst:.2e}"
