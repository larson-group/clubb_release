"""Differentiability / composability tests for the JAX CLUBB building blocks.

The project goal (DESIGN.md) is a *differentiable, composable* JAX CLUBB for ML
and autodiff workflows. The bit-faithful forward pass is verified by
`compare_runs.py`; this file verifies the other half — that the pure-JAX physics
modules support `jax.grad` (and that gradients are *correct*, via finite
differences), so they can be composed into differentiable pipelines.

Scope: the pure-JAX building blocks (saturation, the tridiagonal solver, the PDF
cloud-fraction core, Brunt-Vaisala). Whole-core differentiation is covered by
`clubb_jax/tests/run_jax_timestep_gradient_test.py`.

Run from the repo root: bash tests/run_pytests.sh -jax -include_generated -k test_differentiability
"""
import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq

_SQRT2 = jnp.sqrt(2.0)


def test_composability():
    """Gradient flows through a COMPOSITION of modules (saturation -> chi -> cloud frac),
    demonstrating the modules compose into a differentiable pipeline."""
    rt = 0.012                           # fixed total water
    def chi_of(T):
        rsat = sat_mixrat_liq(jnp.full_like(T, 9.0e4), T, 3)
        return (rt - rsat) / 1.0e-3      # crude chi (scaled excess)
    def pipeline(T):
        return jnp.sum(0.5 * (1.0 + jax.scipy.special.erf(chi_of(T) / _SQRT2)))
    T = jnp.linspace(280.0, 300.0, 40)
    g = jax.grad(pipeline)(T)
    assert bool(jnp.all(jnp.isfinite(g))) and float(jnp.sum(jnp.abs(g))) > 0
    # finite-difference check at the cloud edge (chi ~ 0), where the gradient is
    # well away from the erf-saturated (vanishing-gradient) tails.
    k_edge = int(jnp.argmin(jnp.abs(chi_of(T))))
    eps = 1e-6
    fd = (pipeline(T.at[k_edge].add(eps)) - pipeline(T.at[k_edge].add(-eps))) / (2 * eps)
    rel = abs(float(g[k_edge]) - float(fd)) / (abs(float(fd)) + 1e-30)
    assert rel < 1e-4, f"composed grad wrong: ad={float(g[k_edge]):.4e} fd={float(fd):.4e}"
    print(f"  composability (sat->chi->cloud_frac): grad flows + correct at cloud edge "
          f"(rel {rel:.1e})  PASS")
