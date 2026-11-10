#!/usr/bin/env python3
"""validate the JAX mass-conserving remapping port (remapping_module.F90)."""

from clubb_jax.auto_llm_generated_tests.pytests.array_inputs import call_jax


import numpy as np
import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.remapping_module import calc_mass_over_grid_intervals, remapping_matrix, matrix_vector_mult, remap_vals_to_target, remap_vals_to_target_helper, _pressure_levels

_NG, _DZ, _ZTOP = 2, 50.0, 1500.0


def _rho_spline(ng, ztop):
    # Density spline covering [0, ztop+pad], decreasing with height; ascending levels.
    levels = np.linspace(0.0, ztop + 100.0, 24)
    vals = 1.2 * np.exp(-levels / 8000.0)
    return np.tile(vals, (ng, 1)), np.tile(levels, (ng, 1))


def test_remapping_matrix_identity_and_conservation():
    rng = np.random.default_rng(2)
    # Strictly decreasing pressure levels (as built from hydrostatic integration).
    p = np.cumsum(np.concatenate([[1e5], -rng.uniform(500, 1500, 9)]))[None, :]
    R = np.asarray(remapping_matrix((p).shape[0], (p).shape[1], (p).shape[1], p, p))
    assert np.allclose(R[0], np.eye(9)), "same-level remapping matrix is not identity"
    # Refined target (split each source cell in two): each target row sums to 1 (consistency).
    pt = np.sort(np.concatenate([p[0], 0.5 * (p[0, :-1] + p[0, 1:])]))[::-1][None, :]
    R2 = np.asarray(remapping_matrix((p).shape[0], (p).shape[1], (pt).shape[1], p, pt))
    assert np.allclose(R2.sum(axis=2), 1.0), "remapping matrix rows do not sum to 1"
    print("  remapping_matrix: identity (same grid) + row-sum consistency  PASS")


def test_calc_mass_analytic():
    # Linear density rho(z) = a + b z over [0, 1000]; analytic mass over [z0,z1] = a(z1-z0) + b/2(z1^2-z0^2).
    a, b = 1.2, -3.0e-4
    levels = np.linspace(0, 1000, 11)
    vals = a + b * levels
    grid = np.array([[0.0, 250.0, 600.0, 1000.0]])
    mass = np.asarray(calc_mass_over_grid_intervals((vals[None, :]).shape[0], (vals[None, :]).shape[1], vals[None, :], levels[None, :], (grid).shape[1], grid))[0]
    edges = grid[0]
    exact = a * (edges[1:] - edges[:-1]) + 0.5 * b * (edges[1:] ** 2 - edges[:-1] ** 2)
    assert np.allclose(mass, exact, atol=1e-9), f"mass integral mismatch: {mass} vs {exact}"
    print(f"  calc_mass_over_grid_intervals: matches analytic linear-density integral (max |Δ| "
          f"{np.max(np.abs(mass - exact)):.2e})  PASS")


def test_matvec_and_grad():
    rng = np.random.default_rng(5)
    x = rng.uniform(-1, 1, (2, 4)); A = rng.uniform(-1, 1, (2, 3, 4))
    y = np.asarray(matrix_vector_mult((x).shape[0], (x).shape[1], (A).shape[1], x, A))
    assert np.allclose(y, np.einsum('ikj,ij->ik', A, x)), "matrix_vector_mult != einsum"
    a, b = 1.2, -3.0e-4
    levels = jnp.asarray(np.linspace(0, 1000, 11))
    grid = jnp.asarray([[0.0, 400.0, 1000.0]])
    g = np.asarray(jax.grad(lambda v: jnp.sum(
        calc_mass_over_grid_intervals((v[None, :]).shape[0], (v[None, :]).shape[1], v[None, :], levels[None, :], (grid).shape[1], grid) ** 2))(a + b * levels))
    assert np.isfinite(g).all(), "non-finite grad through calc_mass_over_grid_intervals"
    print(f"  matrix_vector_mult==einsum + finite jax.grad through calc_mass ({g.size} entries)  PASS")


def test_remap_vals_to_target_two_grids():
    """`remap_vals_to_target` (remapping_module.F90:remap_vals_to_target_helper, Ullrich-linear) is the GENERAL
    two-grid conservative remap — distinct source/target grids. The existing tests only cover its building blocks
    and same-grid identity checks. Validate the two defining
    properties: (1) mass conservation — when both grids span the same domain, the pressure-thickness-weighted
    integral is preserved exactly: sum(target·dp_tgt) == sum(source·dp_src); (2) identity — a target grid equal
    to the source grid returns source_values unchanged (the remap matrix is the identity). (iter 526)"""
    ng = _NG
    H = 1000.0
    # Source (10 cells) and target (7 cells) grids spanning the SAME [0, H] domain (required for exact conservation).
    src = np.tile(np.linspace(0.0, H, 11), (ng, 1))
    tgt = np.tile(np.array([0.0, 120.0, 300.0, 470.0, 640.0, 800.0, 930.0, H]), (ng, 1))
    rho_vals, rho_levels = _rho_spline(ng, H)
    p_sfc = np.full(ng, 1.0e5)
    rng = np.random.default_rng(526)
    source_values = rng.uniform(1.0, 5.0, (ng, src.shape[1] - 1))      # arbitrary field on source cells

    target_values = np.asarray(call_jax(remap_vals_to_target_helper, (src).shape[0], (src).shape[1], (tgt).shape[1], src, tgt, (rho_vals).shape[1], rho_vals, rho_levels, source_values, 1, p_sfc, 1))
    # Reconstruct the same internal pressure pipeline to form dp on each grid.
    p_src = np.asarray(_pressure_levels(calc_mass_over_grid_intervals((rho_vals).shape[0], (rho_vals).shape[1], rho_vals, rho_levels, (src).shape[1], src), p_sfc))
    p_tgt = np.asarray(_pressure_levels(calc_mass_over_grid_intervals((rho_vals).shape[0], (rho_vals).shape[1], rho_vals, rho_levels, (tgt).shape[1], tgt), p_sfc))
    dp_src = p_src[:, :-1] - p_src[:, 1:]
    dp_tgt = p_tgt[:, :-1] - p_tgt[:, 1:]
    mass_src = np.sum(source_values * dp_src, axis=1)
    mass_tgt = np.sum(target_values * dp_tgt, axis=1)
    rel = np.max(np.abs(mass_tgt - mass_src) / np.abs(mass_src))
    assert rel < 1e-10, f"two-grid remap not mass-conserving: rel {rel:.2e}"

    # Identity: target grid == source grid -> output == input.
    ident = np.asarray(call_jax(remap_vals_to_target_helper, (src).shape[0], (src).shape[1], (src).shape[1], src, src, (rho_vals).shape[1], rho_vals, rho_levels, source_values, 1, p_sfc, 1))
    iworst = np.max(np.abs(ident - source_values))
    assert iworst < 1e-9, f"same-grid remap_vals_to_target not identity: {iworst:.2e}"

    # Finite grad through the full two-grid pipeline.
    g = np.asarray(jax.grad(lambda v: jnp.sum(
        call_jax(remap_vals_to_target_helper, (src).shape[0], (src).shape[1], (tgt).shape[1], src, tgt, (rho_vals).shape[1], rho_vals, rho_levels, v, 1, p_sfc, 1) ** 2))(jnp.asarray(source_values)))
    assert np.isfinite(g).all(), "non-finite grad through remap_vals_to_target"
    print(f"  remap_vals_to_target (two grids): mass-conserving (rel {rel:.1e}) + identity + finite grad  PASS")


def test_remap_mask_keeps_inactive_bits_and_remaps_active_column():
    """A selected column follows conservative cell overlaps; its neighbor retains signed zero."""
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid

    source_heights = np.tile(np.array([0., 1., 2., 3.]), (2, 1))
    target_heights = source_heights.copy()
    target_heights[1] = [0., .5, 1.5, 3.]
    source_grid = setup_grid(2, 0., 0., 3., grid_type=3, momentum_heights=source_heights)
    target_grid = setup_grid(2, 0., 0., 3., grid_type=3, momentum_heights=target_heights)
    values = jnp.array([[-0., 2., 3.], [1., 3., 5.]])
    rho = jnp.ones((2, 2))
    z = jnp.array([[0., 3.], [0., 3.]])
    mask = jnp.array([False, True])
    l_zt_variable = True
    remap = jax.jit(lambda v: remap_vals_to_target(
        2, source_grid, target_grid, 3, v, 3, 2, rho, z,
        1, jnp.full(2, 100000.), 1, l_zt_variable, mask,
    ))
    result = np.asarray(remap(values))
    assert result.shape == (2, 3)
    assert result[0].tobytes() == np.asarray(values[0]).tobytes()
    np.testing.assert_allclose(result[1], [1., 2., 13./3.], rtol=1e-10, atol=1e-12)
