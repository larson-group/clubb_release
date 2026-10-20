"""Check the JAX diffusion_zt_lhs / diffusion_zm_lhs ports (diffusion.F90) against Fortran."""

import jax.numpy as jnp
import numpy as np

import clubb_python.CLUBB_core.diffusion as fortran

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.diffusion import diffusion_zt_lhs, diffusion_zm_lhs


def call_diffusion_zt_lhs(K_zm, nu, invrs_rho_ds_zt, rho_ds_zm, gr):
    K_zt = jnp.zeros((gr.ngrdcol, gr.nzt), dtype=jnp.float64)
    return diffusion_zt_lhs(
        gr.nzm, gr.nzt, gr.ngrdcol, gr, K_zm, K_zt, nu, invrs_rho_ds_zt, rho_ds_zm
    )


def call_diffusion_zm_lhs(K_zt, nu, invrs_rho_ds_zm, rho_ds_zt, gr):
    K_zm = jnp.zeros((gr.ngrdcol, gr.nzm), dtype=jnp.float64)
    return diffusion_zm_lhs(
        gr.nzm, gr.nzt, gr.ngrdcol, gr, K_zt, K_zm, nu, invrs_rho_ds_zm, rho_ds_zt
    )


def test_f2py_oracle():
    """Validate diffusion_zt_lhs / diffusion_zm_lhs against the f2py Fortran oracle on a STRETCHED grid with
    varying K_zm/K_zt (beyond the uniform-grid constant-coefficient analytic check above). The Fortran takes both
    k_zm and k_zt; the JAX uses only K_zm (the k_zt term is the l_upwind_Kh_dp_term branch, off by default) — passing
    random K_zt and matching confirms the default path ignores it.
    """
    NG, DZ, ZTOP = 2, 40.0, 1200.0
    jgr, fgr = setup_matching_grids(NG, DZ, ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    rng = np.random.default_rng(7)
    K_zm = rng.uniform(0.1, 5.0, (ng, nzm)); K_zt = rng.uniform(0.1, 5.0, (ng, nzt))
    nu = rng.uniform(0.5, 2.0, (ng,))
    rho_zm = rng.uniform(0.8, 1.2, (ng, nzm)); rho_zt = rng.uniform(0.8, 1.2, (ng, nzt))
    irho_zt = 1.0 / rho_zt; irho_zm = 1.0 / rho_zm

    ref_zt = np.asarray(fortran.diffusion_zt_lhs(nzm, nzt, ng, fgr, K_zm, K_zt, nu, irho_zt, rho_zm))
    got_zt = np.asarray(call_diffusion_zt_lhs(jnp.asarray(K_zm), jnp.asarray(nu), jnp.asarray(irho_zt), jnp.asarray(rho_zm), jgr))
    worst_zt = float(max_abs_diff("diffusion", ref_zt, got_zt))
    assert worst_zt < 1e-12, f"diffusion_zt_lhs f2py mismatch {worst_zt:.2e}"

    ref_zm = np.asarray(fortran.diffusion_zm_lhs(nzm, nzt, ng, fgr, K_zt, K_zm, nu, irho_zm, rho_zt))
    got_zm = np.asarray(call_diffusion_zm_lhs(jnp.asarray(K_zt), jnp.asarray(nu), jnp.asarray(irho_zm), jnp.asarray(rho_zt), jgr))
    worst_zm = float(max_abs_diff("diffusion", ref_zm, got_zm))
    assert worst_zm < 1e-12, f"diffusion_zm_lhs f2py mismatch {worst_zm:.2e}"
