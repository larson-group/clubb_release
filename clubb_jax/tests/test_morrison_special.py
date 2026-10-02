"""Checks of DERF1 retained by the active Morrison core."""
import numpy as np
import pytest
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp
from clubb_jax.src.Microphys.Morrison_microphys.module_mp_graupel import DERF1


def test_derf1_vs_scipy():
    """DERF1 (Ooura approximation) == scipy.special.erf to ~double precision over the full range,
    including across both branch boundaries (w=2.2, 6.9) and the saturated tail. The Ooura table is
    a near-double-precision fit, so the match should be tight (≪ the gate), confirming the 130
    coefficients + branch logic were transcribed correctly."""
    _scipy_erf = pytest.importorskip("scipy.special").erf
    x = np.concatenate([np.linspace(-7.0, 7.0, 600),
                        np.array([0.0, 2.2, 2.2 - 1e-9, 6.9, 6.9 - 1e-9, -2.2, -6.9, 1.0, -1.0])])
    j = np.array(DERF1(jnp.asarray(x)))
    s = _scipy_erf(x)
    err = np.abs(j - s)
    assert err.max() < 1e-12, f"DERF1 vs scipy max abs err {err.max():.2e}"
    print(f"  DERF1 == scipy.special.erf (max abs err {err.max():.1e}, incl. branch boundaries)  PASS")


def test_derf1_identities():
    """erf(0)=0, erf is odd, erf(±large)=±1, monotonic increasing."""
    assert abs(float(DERF1(jnp.array(0.0)))) < 1e-15, "erf(0)≠0"
    xs = jnp.array([0.3, 1.7, 3.5, 5.0])
    odd = np.abs(np.array(DERF1(xs)) + np.array(DERF1(-xs)))
    assert odd.max() < 1e-14, "erf not odd"
    assert abs(float(DERF1(jnp.array(8.0))) - 1.0) < 1e-15, "erf(8)≠1"
    assert abs(float(DERF1(jnp.array(-8.0))) + 1.0) < 1e-15, "erf(-8)≠-1"
    g = np.array(DERF1(jnp.linspace(-6.0, 6.0, 200)))
    assert np.all(np.diff(g) >= -1e-15), "erf not monotonic"
    print("  DERF1 identities (zero/odd/saturation/monotonic)  PASS")
