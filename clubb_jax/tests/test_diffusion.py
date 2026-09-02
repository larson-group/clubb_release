"""Tests for JAX diffusion LHS operators vs. Fortran oracle expectations.

Verifies diffusion_zt_lhs and diffusion_zm_lhs against:
  1. Output shape
  2. Boundary conditions (super_top=0, sub_bot=0)
  3. Conservation: weighted column sums = 0 (zero-flux BC property)
  4. Exact numerical values for uniform grid with constant coefficients
  5. Main-diagonal = -(super + sub) conservation identity at every level

Run:
    python clubb_jax/tests/test_diffusion.py
"""
from __future__ import annotations

import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

try:
    import jax
    import jax.numpy as jnp
    jax.config.update("jax_enable_x64", True)
    HAS_JAX = True
except ImportError:
    HAS_JAX = False
    print("JAX not available; skipping tests.")

from clubb_jax.src.CLUBB_core.diffusion import (
    diffusion_zt_lhs, diffusion_zm_lhs,
)
from clubb_jax.src.CLUBB_core.grid_class import setup_grid
# term_dp1_lhs / xp2_xpyp_lhs live in their Fortran home advance_xp2_xpyp_module.py
# (advance_xp2_xpyp_module.F90), not diffusion.F90 — mirror-refactor iter 228 removed the
# dead duplicates from diffusion.py.
from clubb_jax.src.CLUBB_core.advance_xp2_xpyp_module import (
    term_dp1_lhs as _term_dp1_lhs,
    xp2_xpyp_lhs as _xp2_xpyp_lhs,
)


# ──────────────────────────────────────────────────────────────────────────────
# Grid helper (reused from test_grid_ops pattern)
# ──────────────────────────────────────────────────────────────────────────────

def make_even_grid(nzm=5, deltaz=100.0, zm_init=50.0, ngrdcol=1):
    """Evenly spaced ascending grid matching Fortran ascending convention."""
    return setup_grid(
        ngrdcol=ngrdcol,
        deltaz=deltaz,
        zm_init=zm_init,
        zm_top=zm_init + (nzm - 1) * deltaz,
        grid_type=1,
    )


def _to_jax(gr):
    return gr


def term_dp1_lhs(Cn, invrs_tau):
    """Test helper for current explicit term_dp1_lhs."""
    ngrdcol, nzm = Cn.shape
    gr = make_even_grid(nzm=nzm, ngrdcol=ngrdcol)
    return _term_dp1_lhs(nzm, ngrdcol, gr, Cn, invrs_tau)


def xp2_xpyp_lhs(ta, ma, diff, dp1, dt=60.0, gamma=1.5):
    """Test helper for current explicit xp2_xpyp_lhs."""
    if gamma != 1.5:
        raise ValueError("current source uses gamma_over_implicit_ts=1.5")
    ngrdcol, nzm = dp1.shape
    gr = make_even_grid(nzm=nzm, ngrdcol=ngrdcol)
    return _xp2_xpyp_lhs(nzm, ngrdcol, dt, gr, ta, ma, diff, dp1)


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


# ──────────────────────────────────────────────────────────────────────────────
# diffusion_zt_lhs tests
# ──────────────────────────────────────────────────────────────────────────────

def test_diffusion_zt_lhs_shape():
    """Output shape is (3, ngrdcol, nzt)."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, nzt, ngrdcol = 5, 4, 2
    gr = _to_jax(make_even_grid(nzm=nzm, ngrdcol=ngrdcol))
    K_zm = jnp.ones((ngrdcol, nzm))
    nu = jnp.zeros(ngrdcol)
    invrs_rho_ds_zt = jnp.ones((ngrdcol, nzt))
    rho_ds_zm = jnp.ones((ngrdcol, nzm))
    lhs = call_diffusion_zt_lhs(K_zm, nu, invrs_rho_ds_zt, rho_ds_zm, gr)
    assert lhs.shape == (3, ngrdcol, nzt), f"Wrong shape: {lhs.shape}"
    print(f"  shape = {lhs.shape}  PASS")


def test_diffusion_zt_lhs_boundary_zeros():
    """super_top = 0 and sub_bot = 0 (zero-flux BC)."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, nzt, ngrdcol = 6, 5, 1
    gr = _to_jax(make_even_grid(nzm=nzm, ngrdcol=ngrdcol))
    K_zm = 5.0 * jnp.ones((ngrdcol, nzm))
    nu = 1.0 * jnp.ones(ngrdcol)
    invrs_rho_ds_zt = jnp.ones((ngrdcol, nzt))
    rho_ds_zm = 1.2 * jnp.ones((ngrdcol, nzm))
    lhs = np.asarray(call_diffusion_zt_lhs(K_zm, nu, invrs_rho_ds_zt, rho_ds_zm, gr))
    # super at top level (k=nzt-1) must be 0
    assert abs(lhs[0, 0, -1]) < 1e-15, f"super_top = {lhs[0,0,-1]}, expected 0"
    # sub at bottom level (k=0) must be 0
    assert abs(lhs[2, 0, 0]) < 1e-15, f"sub_bot = {lhs[2,0,0]}, expected 0"
    print(f"  super_top={lhs[0,0,-1]:.2e}  sub_bot={lhs[2,0,0]:.2e}  PASS")


def test_diffusion_zt_lhs_conservation():
    """Zero-flux BCs: weighted column sums of LHS = 0.

    Column sum for column j:
      dzt[j-1]*super[j-1] + dzt[j]*main[j] + dzt[j+1]*sub[j+1] = 0
    (using Python 0-based indexing; boundary terms absorb missing neighbors)
    """
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, nzt, ngrdcol = 7, 6, 3
    gr = make_even_grid(nzm=nzm, ngrdcol=ngrdcol)
    gr_j = _to_jax(gr)
    K_zm = 3.0 * jnp.ones((ngrdcol, nzm))
    nu = 0.5 * jnp.ones(ngrdcol)
    invrs_rho_ds_zt = 0.8 * jnp.ones((ngrdcol, nzt))
    rho_ds_zm = 1.25 * jnp.ones((ngrdcol, nzm))
    lhs = np.asarray(call_diffusion_zt_lhs(K_zm, nu, invrs_rho_ds_zt, rho_ds_zm, gr_j))
    super_d, main_d, sub_d = lhs[0], lhs[1], lhs[2]    # (ngrdcol, nzt)
    dzt = np.asarray(1.0 / gr.invrs_dzt)                 # (ngrdcol, nzt)

    # Weighted column sum: col_sum[j] = sum_k{ dzt[k] * lhs_kj }
    col_sum = dzt * main_d
    col_sum[:, 1:] += dzt[:, :-1] * super_d[:, :-1]   # super at k → col k+1
    col_sum[:, :-1] += dzt[:, 1:] * sub_d[:, 1:]      # sub at k → col k-1

    max_err = np.max(np.abs(col_sum))
    print(f"  max column sum = {max_err:.3e}  PASS" if max_err < 1e-12
          else f"  FAIL max_err={max_err}")
    assert max_err < 1e-12, f"Conservation violated: max col sum = {max_err}"


def test_diffusion_zt_lhs_values():
    """Exact numerical values for nzm=5, uniform grid, K=10, nu=0, rho=1.

    With deltaz=100m: invrs_dzt=invrs_dzm=0.01, K_zm_nu=10.
    C = 0.01 * 1 * 10 * 1 * 0.01 = 0.001
    Expected lhs (nzt=4):
      super = [-C, -C, -C,  0]
      main  = [+C, 2C, 2C, +C]
      sub   = [ 0, -C, -C, -C]
    """
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, nzt, ngrdcol = 5, 4, 1
    gr = _to_jax(make_even_grid(nzm=nzm, deltaz=100.0, ngrdcol=ngrdcol))
    K_zm = 10.0 * jnp.ones((ngrdcol, nzm))
    nu = jnp.zeros(ngrdcol)
    invrs_rho_ds_zt = jnp.ones((ngrdcol, nzt))
    rho_ds_zm = jnp.ones((ngrdcol, nzm))
    lhs = np.asarray(call_diffusion_zt_lhs(K_zm, nu, invrs_rho_ds_zt, rho_ds_zm, gr))
    C = 0.001
    expected_super = np.array([-C, -C, -C,  0.0])
    expected_main  = np.array([ C, 2*C, 2*C, C])
    expected_sub   = np.array([ 0.0, -C, -C, -C])
    err = max(
        np.max(np.abs(lhs[0, 0] - expected_super)),
        np.max(np.abs(lhs[1, 0] - expected_main)),
        np.max(np.abs(lhs[2, 0] - expected_sub)),
    )
    print(f"  super={lhs[0,0]}, main={lhs[1,0]}, sub={lhs[2,0]}")
    print(f"  max_err = {err:.3e}  PASS" if err < 1e-15 else f"  FAIL err={err}")
    assert err < 1e-15


def test_diffusion_zt_lhs_main_diag_identity():
    """main[k] = -(super[k] + sub[k]) at every level (conservation identity)."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, nzt, ngrdcol = 8, 7, 2
    gr = _to_jax(make_even_grid(nzm=nzm, ngrdcol=ngrdcol))
    rng = np.random.default_rng(7)
    K_zm = jnp.array(rng.uniform(0.1, 5.0, (ngrdcol, nzm)))
    nu = jnp.array(rng.uniform(0.0, 0.5, ngrdcol))
    invrs_rho_ds_zt = jnp.array(rng.uniform(0.5, 2.0, (ngrdcol, nzt)))
    rho_ds_zm = jnp.array(rng.uniform(0.8, 1.5, (ngrdcol, nzm)))
    lhs = np.asarray(call_diffusion_zt_lhs(K_zm, nu, invrs_rho_ds_zt, rho_ds_zm, gr))
    err = np.max(np.abs(lhs[1] + lhs[0] + lhs[2]))
    print(f"  max |main + super + sub| = {err:.3e}  PASS" if err < 1e-13
          else f"  FAIL err={err}")
    assert err < 1e-13


# ──────────────────────────────────────────────────────────────────────────────
# diffusion_zm_lhs tests
# ──────────────────────────────────────────────────────────────────────────────

def test_diffusion_zm_lhs_shape():
    """Output shape is (3, ngrdcol, nzm)."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, nzt, ngrdcol = 5, 4, 2
    gr = _to_jax(make_even_grid(nzm=nzm, ngrdcol=ngrdcol))
    K_zt = jnp.ones((ngrdcol, nzt))
    nu = jnp.zeros(ngrdcol)
    invrs_rho_ds_zm = jnp.ones((ngrdcol, nzm))
    rho_ds_zt = jnp.ones((ngrdcol, nzt))
    lhs = call_diffusion_zm_lhs(K_zt, nu, invrs_rho_ds_zm, rho_ds_zt, gr)
    assert lhs.shape == (3, ngrdcol, nzm), f"Wrong shape: {lhs.shape}"
    print(f"  shape = {lhs.shape}  PASS")


def test_diffusion_zm_lhs_boundary_zeros():
    """super_top = 0 and sub_bot = 0 (zero-flux BC)."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, nzt, ngrdcol = 6, 5, 1
    gr = _to_jax(make_even_grid(nzm=nzm, ngrdcol=ngrdcol))
    K_zt = 4.0 * jnp.ones((ngrdcol, nzt))
    nu = 0.2 * jnp.ones(ngrdcol)
    invrs_rho_ds_zm = jnp.ones((ngrdcol, nzm))
    rho_ds_zt = 1.1 * jnp.ones((ngrdcol, nzt))
    lhs = np.asarray(call_diffusion_zm_lhs(K_zt, nu, invrs_rho_ds_zm, rho_ds_zt, gr))
    assert abs(lhs[0, 0, -1]) < 1e-15, f"super_top = {lhs[0,0,-1]}, expected 0"
    assert abs(lhs[2, 0, 0]) < 1e-15, f"sub_bot = {lhs[2,0,0]}, expected 0"
    print(f"  super_top={lhs[0,0,-1]:.2e}  sub_bot={lhs[2,0,0]:.2e}  PASS")


def test_diffusion_zm_lhs_conservation():
    """Zero-flux BCs: weighted column sums = 0 for interior+upper levels.

    The lower boundary (k=0) is intentionally excluded because it is not
    fed into the solver (per Fortran comment).

    Column sum for column j (j=1..nzm-1):
      dzm[j-1]*super[j-1] + dzm[j]*main[j] + dzm[j+1]*sub[j+1] = 0
    """
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, nzt, ngrdcol = 7, 6, 3
    gr = make_even_grid(nzm=nzm, ngrdcol=ngrdcol)
    gr_j = _to_jax(gr)
    K_zt = 2.5 * jnp.ones((ngrdcol, nzt))
    nu = 0.3 * jnp.ones(ngrdcol)
    invrs_rho_ds_zm = (1.0 / 1.2) * jnp.ones((ngrdcol, nzm))
    rho_ds_zt = 1.2 * jnp.ones((ngrdcol, nzt))
    lhs = np.asarray(call_diffusion_zm_lhs(K_zt, nu, invrs_rho_ds_zm, rho_ds_zt, gr_j))
    super_d, main_d, sub_d = lhs[0], lhs[1], lhs[2]   # (ngrdcol, nzm)
    dzm = np.asarray(1.0 / gr.invrs_dzm)                # (ngrdcol, nzm)

    # Check conservation for columns j=1..nzm-1 (lower boundary excluded)
    col_sum = dzm * main_d
    col_sum[:, 1:] += dzm[:, :-1] * super_d[:, :-1]
    col_sum[:, :-1] += dzm[:, 1:] * sub_d[:, 1:]

    # Only check j=1..nzm-1 (exclude j=0 which uses unused lower boundary row)
    max_err = np.max(np.abs(col_sum[:, 1:]))
    print(f"  max column sum (j>=1) = {max_err:.3e}  PASS" if max_err < 1e-12
          else f"  FAIL max_err={max_err}")
    assert max_err < 1e-12, f"Conservation violated: max col sum = {max_err}"


def test_diffusion_zm_lhs_values():
    """Exact numerical values for nzm=5, uniform grid, K_zt=1, nu=0, rho=1.

    With deltaz=100m: invrs_dzm=invrs_dzt=0.01, K_zt_nu=1.
    C = 0.01 * 1 * 1 * 1 * 0.01 = 0.0001
    Expected lhs (nzm=5):
      super = [-C, -C, -C, -C,  0]
      main  = [+C, 2C, 2C, 2C, +C]
      sub   = [ 0, -C, -C, -C, -C]
    """
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, nzt, ngrdcol = 5, 4, 1
    gr = _to_jax(make_even_grid(nzm=nzm, deltaz=100.0, ngrdcol=ngrdcol))
    K_zt = 1.0 * jnp.ones((ngrdcol, nzt))
    nu = jnp.zeros(ngrdcol)
    invrs_rho_ds_zm = jnp.ones((ngrdcol, nzm))
    rho_ds_zt = jnp.ones((ngrdcol, nzt))
    lhs = np.asarray(call_diffusion_zm_lhs(K_zt, nu, invrs_rho_ds_zm, rho_ds_zt, gr))
    C = 0.0001
    expected_super = np.array([-C, -C, -C, -C,  0.0])
    expected_main  = np.array([ C, 2*C, 2*C, 2*C, C])
    expected_sub   = np.array([ 0.0, -C, -C, -C, -C])
    err = max(
        np.max(np.abs(lhs[0, 0] - expected_super)),
        np.max(np.abs(lhs[1, 0] - expected_main)),
        np.max(np.abs(lhs[2, 0] - expected_sub)),
    )
    print(f"  super={lhs[0,0]}")
    print(f"  main ={lhs[1,0]}")
    print(f"  sub  ={lhs[2,0]}")
    print(f"  max_err = {err:.3e}  PASS" if err < 1e-15 else f"  FAIL err={err}")
    assert err < 1e-15


def test_diffusion_zm_lhs_main_diag_identity():
    """main[k] = -(super[k] + sub[k]) at every level."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, nzt, ngrdcol = 8, 7, 2
    gr = _to_jax(make_even_grid(nzm=nzm, ngrdcol=ngrdcol))
    rng = np.random.default_rng(13)
    K_zt = jnp.array(rng.uniform(0.1, 5.0, (ngrdcol, nzt)))
    nu = jnp.array(rng.uniform(0.0, 0.5, ngrdcol))
    invrs_rho_ds_zm = jnp.array(rng.uniform(0.5, 2.0, (ngrdcol, nzm)))
    rho_ds_zt = jnp.array(rng.uniform(0.8, 1.5, (ngrdcol, nzt)))
    lhs = np.asarray(call_diffusion_zm_lhs(K_zt, nu, invrs_rho_ds_zm, rho_ds_zt, gr))
    err = np.max(np.abs(lhs[1] + lhs[0] + lhs[2]))
    print(f"  max |main + super + sub| = {err:.3e}  PASS" if err < 1e-13
          else f"  FAIL err={err}")
    assert err < 1e-13


# ──────────────────────────────────────────────────────────────────────────────
# term_dp1_lhs tests (Iter 7)
# ──────────────────────────────────────────────────────────────────────────────

def test_term_dp1_lhs_shape():
    """Output shape is (ngrdcol, nzm)."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, ngrdcol = 6, 2
    Cn = jnp.ones((ngrdcol, nzm))
    invrs_tau = jnp.ones((ngrdcol, nzm))
    lhs = term_dp1_lhs(Cn, invrs_tau)
    assert lhs.shape == (ngrdcol, nzm), f"Expected ({ngrdcol},{nzm}), got {lhs.shape}"
    print(f"  shape = {lhs.shape}  PASS")


def test_term_dp1_lhs_boundary_zeros():
    """Boundaries k=0 and k=nzm-1 are zero."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, ngrdcol = 7, 2
    rng = np.random.default_rng(42)
    Cn = jnp.array(rng.uniform(0.1, 2.0, (ngrdcol, nzm)))
    invrs_tau = jnp.array(rng.uniform(0.01, 1.0, (ngrdcol, nzm)))
    lhs = np.asarray(term_dp1_lhs(Cn, invrs_tau))
    assert np.all(lhs[:, 0] == 0.0), f"Lower boundary not zero: {lhs[:, 0]}"
    assert np.all(lhs[:, -1] == 0.0), f"Upper boundary not zero: {lhs[:, -1]}"
    print(f"  lower={lhs[:, 0]}, upper={lhs[:, -1]}  PASS")


def test_term_dp1_lhs_interior_values():
    """Interior values equal Cn * invrs_tau_zm, bit-exact."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, ngrdcol = 5, 1
    Cn = jnp.array([[0.5, 0.7, 0.9, 1.1, 1.3]])
    invrs_tau = jnp.array([[0.1, 0.2, 0.3, 0.4, 0.5]])
    lhs = np.asarray(term_dp1_lhs(Cn, invrs_tau))
    # Boundaries must be zero
    assert lhs[0, 0] == 0.0 and lhs[0, -1] == 0.0
    # Interior: lhs[k] = Cn[k] * invrs_tau[k] for k=1..nzm-2
    expected_interior = np.asarray(Cn)[0, 1:-1] * np.asarray(invrs_tau)[0, 1:-1]
    err = np.max(np.abs(lhs[0, 1:-1] - expected_interior))
    print(f"  interior={lhs[0, 1:-1]}, expected={expected_interior}")
    print(f"  max_err = {err:.3e}  PASS" if err == 0.0 else f"  FAIL err={err}")
    assert err == 0.0, f"Interior mismatch: {err}"


def test_term_dp1_lhs_multi_col():
    """Correct for ngrdcol=3 with spatially-varying Cn and invrs_tau."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, ngrdcol = 8, 3
    rng = np.random.default_rng(99)
    Cn_np = rng.uniform(0.1, 3.0, (ngrdcol, nzm))
    tau_np = rng.uniform(0.01, 2.0, (ngrdcol, nzm))
    lhs = np.asarray(term_dp1_lhs(jnp.array(Cn_np), jnp.array(tau_np)))
    ref = Cn_np * tau_np
    ref[:, 0] = 0.0;  ref[:, -1] = 0.0
    err = np.max(np.abs(lhs - ref))
    print(f"  max_err = {err:.3e}  PASS" if err == 0.0 else f"  FAIL err={err}")
    assert err == 0.0


# ──────────────────────────────────────────────────────────────────────────────
# xp2_xpyp_lhs tests (Iter 8)
# ──────────────────────────────────────────────────────────────────────────────

def _make_random_3band(ngrdcol, nzm, seed):
    """Random (3, ngrdcol, nzm) tridiagonal with zero boundaries."""
    rng = np.random.default_rng(seed)
    arr = rng.uniform(-1.0, 1.0, (3, ngrdcol, nzm))
    arr[:, :, 0] = 0.0;  arr[:, :, -1] = 0.0
    return arr


def test_xp2_xpyp_lhs_shape():
    """Output shape is (3, ngrdcol, nzm)."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, ngrdcol = 7, 2
    ta = jnp.array(_make_random_3band(ngrdcol, nzm, 1))
    ma = jnp.array(_make_random_3band(ngrdcol, nzm, 2))
    diff = jnp.array(_make_random_3band(ngrdcol, nzm, 3))
    dp1 = jnp.zeros((ngrdcol, nzm))
    lhs = xp2_xpyp_lhs(ta, ma, diff, dp1, dt=60.0)
    assert lhs.shape == (3, ngrdcol, nzm), f"Expected (3,{ngrdcol},{nzm}), got {lhs.shape}"
    print(f"  shape = {lhs.shape}  PASS")


def test_xp2_xpyp_lhs_boundary_conditions():
    """Boundaries k=0 and k=nzm-1 are [0, 1, 0] for all columns."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, ngrdcol = 6, 2
    rng = np.random.default_rng(7)
    ta = jnp.array(rng.uniform(-1.0, 1.0, (3, ngrdcol, nzm)))
    ma = jnp.array(rng.uniform(-1.0, 1.0, (3, ngrdcol, nzm)))
    diff = jnp.array(rng.uniform(-1.0, 1.0, (3, ngrdcol, nzm)))
    dp1 = jnp.array(rng.uniform(0.0, 1.0, (ngrdcol, nzm)))
    lhs = np.asarray(xp2_xpyp_lhs(ta, ma, diff, dp1, dt=30.0))
    for bnd in [0, -1]:
        assert np.all(lhs[0, :, bnd] == 0.0), f"super at bnd={bnd}: {lhs[0,:,bnd]}"
        assert np.all(lhs[1, :, bnd] == 1.0), f"main at bnd={bnd}: {lhs[1,:,bnd]}"
        assert np.all(lhs[2, :, bnd] == 0.0), f"sub at bnd={bnd}: {lhs[2,:,bnd]}"
    print(f"  lower=[{lhs[0,0,0]},{lhs[1,0,0]},{lhs[2,0,0]}]  upper=[{lhs[0,0,-1]},{lhs[1,0,-1]},{lhs[2,0,-1]}]  PASS")


def test_xp2_xpyp_lhs_interior_formula():
    """Interior matches numpy reference assembly exactly."""
    if not HAS_JAX:
        print("  SKIP (JAX not available)")
        return
    nzm, ngrdcol = 8, 2
    gamma = 1.5
    dt = 60.0
    rng = np.random.default_rng(55)
    ta_np   = _make_random_3band(ngrdcol, nzm, 11)
    ma_np   = _make_random_3band(ngrdcol, nzm, 12)
    diff_np = _make_random_3band(ngrdcol, nzm, 13)
    dp1_np  = np.abs(rng.uniform(0.0, 0.5, (ngrdcol, nzm)))
    dp1_np[:, 0] = 0.0;  dp1_np[:, -1] = 0.0

    lhs_jax = np.asarray(xp2_xpyp_lhs(
        jnp.array(ta_np), jnp.array(ma_np), jnp.array(diff_np),
        jnp.array(dp1_np), dt=dt, gamma=gamma,
    ))

    # numpy reference for interior k=1..nzm-2
    super_ref = diff_np[0,:,1:-1] + ma_np[0,:,1:-1] + ta_np[0,:,1:-1] * gamma
    main_ref  = (diff_np[1,:,1:-1] + ma_np[1,:,1:-1] + ta_np[1,:,1:-1] * gamma
                 + dp1_np[:,1:-1] + 1.0/dt)
    sub_ref   = diff_np[2,:,1:-1] + ma_np[2,:,1:-1] + ta_np[2,:,1:-1] * gamma

    err = max(
        np.max(np.abs(lhs_jax[0, :, 1:-1] - super_ref)),
        np.max(np.abs(lhs_jax[1, :, 1:-1] - main_ref)),
        np.max(np.abs(lhs_jax[2, :, 1:-1] - sub_ref)),
    )
    print(f"  interior max_err = {err:.3e}  PASS" if err < 1e-14 else f"  FAIL err={err}")
    assert err < 1e-14


# ──────────────────────────────────────────────────────────────────────────────
# f2py Fortran-oracle bit-shadow (general stretched grid + varying K)
# ──────────────────────────────────────────────────────────────────────────────

def test_f2py_oracle():
    """Validate diffusion_zt_lhs / diffusion_zm_lhs against the f2py Fortran oracle on a STRETCHED grid with
    varying K_zm/K_zt (beyond the uniform-grid constant-coefficient analytic check above). The Fortran takes both
    k_zm and k_zt; the JAX uses only K_zm (the k_zt term is the l_upwind_Kh_dp_term branch, off by default) — passing
    random K_zt and matching confirms the default path ignores it. SKIPs if clubb_f2py / clubb_python is unbuilt.
    """
    for p in (os.path.join(os.path.dirname(__file__), '..', '..'),
              os.path.join(os.path.dirname(__file__), '..', '..', 'clubb_python_api')):
        ap = os.path.abspath(p)
        if ap not in sys.path:
            sys.path.append(ap)
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py diffusion oracle: SKIP ({type(e).__name__})")
        return
    NG, DZ, ZTOP = 2, 40.0, 1200.0
    jgr = setup_grid(ngrdcol=NG, deltaz=DZ, zm_init=0.0, zm_top=ZTOP, grid_type=1)  # stretched grid
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    clubb_api.init_err_info(ng)
    cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng),
                         l_implemented=False, l_ascending_grid=True, grid_type=2,
                         deltaz=np.full(ng, DZ), zm_init=np.zeros(ng), zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)),
                         err_info=ErrInfo(ngrdcol=ng))
    rng = np.random.default_rng(7)
    K_zm = rng.uniform(0.1, 5.0, (ng, nzm)); K_zt = rng.uniform(0.1, 5.0, (ng, nzt))
    nu = rng.uniform(0.5, 2.0, (ng,))
    rho_zm = rng.uniform(0.8, 1.2, (ng, nzm)); rho_zt = rng.uniform(0.8, 1.2, (ng, nzt))
    irho_zt = 1.0 / rho_zt; irho_zm = 1.0 / rho_zm

    ref_zt = np.asarray(clubb_f2py.f2py_diffusion_zt_lhs(
        np.asfortranarray(K_zm), np.asfortranarray(K_zt), nu, np.asfortranarray(irho_zt), np.asfortranarray(rho_zm)))
    got_zt = np.asarray(call_diffusion_zt_lhs(jnp.asarray(K_zm), jnp.asarray(nu), jnp.asarray(irho_zt), jnp.asarray(rho_zm), jgr))
    worst_zt = float(np.max(np.abs(ref_zt - got_zt)))
    assert worst_zt < 1e-12, f"diffusion_zt_lhs f2py mismatch {worst_zt:.2e}"

    ref_zm = np.asarray(clubb_f2py.f2py_diffusion_zm_lhs(
        np.asfortranarray(K_zt), np.asfortranarray(K_zm), nu, np.asfortranarray(irho_zm), np.asfortranarray(rho_zt)))
    got_zm = np.asarray(call_diffusion_zm_lhs(jnp.asarray(K_zt), jnp.asarray(nu), jnp.asarray(irho_zm), jnp.asarray(rho_zt), jgr))
    worst_zm = float(np.max(np.abs(ref_zm - got_zm)))
    assert worst_zm < 1e-12, f"diffusion_zm_lhs f2py mismatch {worst_zm:.2e}"
    print(f"  f2py diffusion zt/zm LHS: bit-match on stretched grid + varying K, worst {max(worst_zt, worst_zm):.2e}  PASS")


# ──────────────────────────────────────────────────────────────────────────────
# Main runner
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("=" * 60)
    print("CLUBB JAX diffusion LHS tests")
    print("=" * 60)

    if not HAS_JAX:
        print("ERROR: JAX not available. Install with: pip install 'jax[cpu]'")
        sys.exit(1)

    tests = [
        test_diffusion_zt_lhs_shape,
        test_diffusion_zt_lhs_boundary_zeros,
        test_diffusion_zt_lhs_conservation,
        test_diffusion_zt_lhs_values,
        test_diffusion_zt_lhs_main_diag_identity,
        test_diffusion_zm_lhs_shape,
        test_diffusion_zm_lhs_boundary_zeros,
        test_diffusion_zm_lhs_conservation,
        test_diffusion_zm_lhs_values,
        test_diffusion_zm_lhs_main_diag_identity,
        test_term_dp1_lhs_shape,
        test_term_dp1_lhs_boundary_zeros,
        test_term_dp1_lhs_interior_values,
        test_term_dp1_lhs_multi_col,
        test_xp2_xpyp_lhs_shape,
        test_xp2_xpyp_lhs_boundary_conditions,
        test_xp2_xpyp_lhs_interior_formula,
        test_f2py_oracle,
    ]

    passed = 0
    failed = 0
    for t in tests:
        print(f"\n{t.__name__}:")
        try:
            t()
            passed += 1
        except AssertionError as e:
            print(f"  FAIL: {e}")
            failed += 1
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  ERROR: {type(e).__name__}: {e}")
            failed += 1

    print("\n" + "=" * 60)
    print(f"Results: {passed}/{passed+failed} passed, {failed} failed")
    print("=" * 60)
    sys.exit(0 if failed == 0 else 1)
