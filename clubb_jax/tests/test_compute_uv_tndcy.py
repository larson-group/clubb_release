#!/usr/bin/env python3
"""test_compute_uv_tndcy.py — pin the um/vm Coriolis + geostrophic + forcing tendency.

`compute_uv_tndcy` (advance_windm_edsclrm_module.py ↔ advance_windm_edsclrm_module.F90:compute_uv_tndcy, F90:2051)
builds the explicit momentum tendency that the wind advance seeds the implicit eddy-diffusion solve with. It is on the
live path (every standalone case advances u/v under rotation) but is validated only end-to-end (compare_cases) — its
Coriolis SIGN convention was never pinned in isolation, so a sign flip would only surface as a full-case divergence.

The Fortran computes one component per call via `solve_type`:
  windm_edsclrm_um:  xm_gf = −fcor·perp_wind_g(=vg);  xm_cf =  fcor·perp_wind_m(=vm)   → um_tndcy = fcor·(vm−vg) + um_f
  windm_edsclrm_vm:  xm_gf =  fcor·perp_wind_g(=ug);  xm_cf = −fcor·perp_wind_m(=um)   → vm_tndcy = fcor·(ug−um) + vm_f
The JAX consolidates both components into one call returning (um_tndcy, vm_tndcy). This pins both against an
INDEPENDENT transcription of the F90 sign convention (so the +fcor/−fcor placement is checked directly), plus a
geostrophic-balance property (zero tendency when wind == geostrophic and no forcing) and a finite jax.grad.
Oracle-independent (closed-form); never SKIPs. (iter 543)
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_windm_edsclrm_module import compute_uv_tndcy

_NG, _NZT = 2, 8


def _inputs(rng):
    fcor = rng.uniform(0.5e-4, 1.5e-4, (_NG,))           # mid-latitude Coriolis
    ug = rng.uniform(-12.0, 12.0, (_NG, _NZT))
    vg = rng.uniform(-12.0, 12.0, (_NG, _NZT))
    um = rng.uniform(-12.0, 12.0, (_NG, _NZT))
    vm = rng.uniform(-12.0, 12.0, (_NG, _NZT))
    um_f = rng.uniform(-1e-3, 1e-3, (_NG, _NZT))
    vm_f = rng.uniform(-1e-3, 1e-3, (_NG, _NZT))
    return fcor, ug, vg, um, vm, um_f, vm_f


def test_matches_f90_sign_convention():
    rng = np.random.default_rng(543)
    fcor, ug, vg, um, vm, um_f, vm_f = _inputs(rng)
    um_t, vm_t = compute_uv_tndcy(jnp.asarray(fcor), jnp.asarray(ug), jnp.asarray(vg),
                                  jnp.asarray(um), jnp.asarray(vm), jnp.asarray(um_f), jnp.asarray(vm_f))
    f = fcor[:, None]
    # Independent transcription of the F90 per-component formula (windm_edsclrm_um / _vm):
    ref_um = -f * vg + f * vm + um_f       # xm_gf(=−fcor·vg) + xm_cf(=+fcor·vm) + forcing
    ref_vm = f * ug - f * um + vm_f        # xm_gf(=+fcor·ug) + xm_cf(=−fcor·um) + forcing
    eu = float(np.max(np.abs(np.asarray(um_t) - ref_um)))
    ev = float(np.max(np.abs(np.asarray(vm_t) - ref_vm)))
    assert eu < 1e-15 and ev < 1e-15, f"um/vm tndcy mismatch vs F90 formula: {eu:.2e}/{ev:.2e}"
    print(f"  um/vm tndcy match the F90 Coriolis sign convention (max |Δ| {max(eu, ev):.1e})  PASS")


def test_geostrophic_balance():
    """With wind == geostrophic wind and zero forcing, the Coriolis+geostrophic tendency must vanish exactly."""
    rng = np.random.default_rng(7)
    fcor = rng.uniform(0.5e-4, 1.5e-4, (_NG,))
    ug = rng.uniform(-10.0, 10.0, (_NG, _NZT)); vg = rng.uniform(-10.0, 10.0, (_NG, _NZT))
    zero = np.zeros((_NG, _NZT))
    um_t, vm_t = compute_uv_tndcy(jnp.asarray(fcor), jnp.asarray(ug), jnp.asarray(vg),
                                  jnp.asarray(ug), jnp.asarray(vg), jnp.asarray(zero), jnp.asarray(zero))
    assert np.max(np.abs(np.asarray(um_t))) < 1e-18 and np.max(np.abs(np.asarray(vm_t))) < 1e-18, \
        "tendency nonzero at geostrophic balance (um=ug, vm=vg, no forcing)"
    print("  geostrophic balance: zero tendency when (um,vm)==(ug,vg) and no forcing  PASS")


def test_forcing_additive_and_grad():
    rng = np.random.default_rng(99)
    fcor, ug, vg, um, vm, um_f, vm_f = _inputs(rng)
    # Forcing enters additively: tndcy(forcing) − tndcy(0) == forcing.
    a_um, a_vm = compute_uv_tndcy(jnp.asarray(fcor), jnp.asarray(ug), jnp.asarray(vg),
                                  jnp.asarray(um), jnp.asarray(vm), jnp.asarray(um_f), jnp.asarray(vm_f))
    z = np.zeros((_NG, _NZT))
    b_um, b_vm = compute_uv_tndcy(jnp.asarray(fcor), jnp.asarray(ug), jnp.asarray(vg),
                                  jnp.asarray(um), jnp.asarray(vm), jnp.asarray(z), jnp.asarray(z))
    assert np.allclose(np.asarray(a_um) - np.asarray(b_um), um_f, atol=1e-15) and \
           np.allclose(np.asarray(a_vm) - np.asarray(b_vm), vm_f, atol=1e-15), "forcing not additive"
    g = jax.grad(lambda u: jnp.sum(compute_uv_tndcy(jnp.asarray(fcor), jnp.asarray(ug), jnp.asarray(vg),
                u, jnp.asarray(vm), jnp.asarray(um_f), jnp.asarray(vm_f))[0] ** 2))(jnp.asarray(um))
    assert np.all(np.isfinite(np.asarray(g))), "non-finite grad wrt um"
    print("  forcing additive + finite jax.grad wrt um  PASS")


def main():
    print("test_compute_uv_tndcy:")
    test_matches_f90_sign_convention()
    test_geostrophic_balance()
    test_forcing_additive_and_grad()
    print("All compute_uv_tndcy checks PASSED")


if __name__ == "__main__":
    main()
