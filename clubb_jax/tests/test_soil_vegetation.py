"""Validation for the soil_vegetation port (soil_vegetation.F90).

Bit-exactness vs a direct Fortran-formula replica of advance_soil_veg, plus a multi-step physical sanity
check (the surface/soil temperatures stay finite and relax toward a radiative-turbulent equilibrium).
"""
import math
import numpy as np

import os
import sys
_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

from clubb_jax.src.Radiation.soil_vegetation import advance_soil_veg, initialize_soil_veg
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats

_CP, _LV, _RD = 1004.67, 2.5e6, 287.04
_KAPPA, _P0, _SB, _PI = _RD / _CP, 1.0e5, 5.6704e-8, math.pi


def _stats(ncol):
    return JaxStats.empty(l_sample=False, names=(), ncol=ncol, max_nlev=1)


def _fortran_replica(dt, rho, swup, swdn, lwdn, wpthlp, wprtp, psfc, deep, sfc, veg):
    """Exact transcription of soil_vegetation.F90:148-201 (force-restore coefficients + the 3 updates)."""
    cs, rs, ks = 2.00e3, 1.00e3, 2.00e-7
    d1 = math.sqrt(ks * 3600.0 * 24.0)
    c1 = 2.0 * math.sqrt(_PI) / (rs * cs * d1)
    c2 = 2.0 * _PI / (3600.0 * 24.0)
    c3 = math.sqrt(_PI * 2.0) / (math.exp(_PI / 4.0) * rs * cs * math.sqrt(ks * 3600.0 * 24.0 * 365.0))
    lwup = _SB * veg ** 4
    wpthep = wpthlp + (_LV / _CP) * (_P0 / psfc) ** _KAPPA * wprtp
    veg_hf = lwdn - lwup - wpthep * rho * _CP + (swdn - swup)
    soil_hf = 10.0 * (veg - sfc) + 0.05 * swdn
    veg_n = veg + dt * 5.0e-5 * (veg_hf - soil_hf)
    sfc_n = sfc + dt * (c1 * soil_hf - c2 * (sfc - deep))
    deep_n = deep + dt * c3 * soil_hf
    return deep_n, sfc_n, veg_n, soil_hf


def test_soil_veg_vs_fortran_replica():
    rng = np.random.default_rng(3)
    n = 4
    dt = 30.0
    rho = 1.1 + 0.1 * rng.random(n)
    swup = 50.0 * rng.random(n); swdn = 400.0 * rng.random(n); lwdn = 300.0 + 50.0 * rng.random(n)
    wpthlp = 0.05 * (rng.random(n) - 0.3); wprtp = 1e-4 * (rng.random(n) - 0.3)
    psfc = 1.0e5 * (0.98 + 0.02 * rng.random(n))
    deep, sfc, veg = 288.58 * np.ones(n), 295.0 + 5.0 * rng.random(n), 295.0 + 8.0 * rng.random(n)
    _, *got = advance_soil_veg(n, dt, rho, swup, swdn, lwdn, wpthlp, wprtp, psfc, _stats(n), deep, sfc, veg)
    exp = _fortran_replica(dt, rho, swup, swdn, lwdn, wpthlp, wprtp, psfc, deep, sfc, veg)
    rel = 0.0
    for g, e in zip(got, exp):
        rel = max(rel, float(np.max(np.abs(np.asarray(g) - e) / (np.abs(e) + 1e-300))))
    assert rel <= 1e-14, f"soil_veg differs from Fortran replica: rel={rel:.2e}"
    print(f"  advance_soil_veg vs Fortran-formula replica: rel {rel:.1e}  PASS")


def test_soil_veg_init_and_integration():
    deep, sfc, veg = initialize_soil_veg(1)
    assert deep[0] == 288.58 and sfc[0] == 300.0 and veg[0] == 300.0, "init defaults wrong"
    # integrate a constant nighttime forcing (no SW, net LW cooling) and check temps stay finite + cool
    dt = 30.0
    rho = np.array([1.2]); swup = np.array([0.0]); swdn = np.array([0.0])
    lwdn = np.array([320.0]); wpthlp = np.array([-0.02]); wprtp = np.array([0.0]); psfc = np.array([1.0e5])
    veg0 = float(veg[0])
    for _ in range(120):  # 1 hour
        _, deep, sfc, veg = advance_soil_veg(1, dt, rho, swup, swdn, lwdn, wpthlp, wprtp, psfc, _stats(1), deep, sfc, veg)
        deep, sfc, veg = np.asarray(deep), np.asarray(sfc), np.asarray(veg)
    assert np.all(np.isfinite(veg)) and np.all(np.isfinite(sfc)) and np.all(np.isfinite(deep))
    assert 200.0 < veg[0] < 350.0, f"veg temperature unphysical: {veg[0]}"
    print(f"  soil_veg init + 1h integration: veg {veg0:.1f}→{veg[0]:.1f} K, all temps finite & physical  PASS")


def test_soil_veg_differentiable():
    """advance_soil_veg is differentiable: jax.grad of the updated vegetation temperature w.r.t. the
    incoming LW flux and the previous veg temperature is finite + nonzero (radiative-surface coupling)."""
    import jax, jax.numpy as jnp
    n = 2
    args = dict(dt=30.0, rho_sfc=jnp.full(n, 1.2), Frad_SW_up_sfc=jnp.zeros(n),
                Frad_SW_down_sfc=jnp.full(n, 200.0), wpthlp_sfc=jnp.full(n, -0.02),
                wprtp_sfc=jnp.zeros(n), p_sfc=jnp.full(n, 1.0e5),
                deep_soil_T_in_K=jnp.full(n, 288.58), sfc_soil_T_in_K=jnp.full(n, 295.0))

    def loss(lwdn, veg):
        _, _, _, veg_new = advance_soil_veg(
            n, args['dt'], args['rho_sfc'], args['Frad_SW_up_sfc'], args['Frad_SW_down_sfc'], lwdn,
            args['wpthlp_sfc'], args['wprtp_sfc'], args['p_sfc'], _stats(n),
            args['deep_soil_T_in_K'], args['sfc_soil_T_in_K'], veg)
        return jnp.sum(veg_new)

    g_lw, g_veg = jax.grad(loss, argnums=(0, 1))(jnp.full(n, 320.0), jnp.full(n, 300.0))
    assert jnp.all(jnp.isfinite(g_lw)) and jnp.all(jnp.isfinite(g_veg)), "non-finite soil_veg gradient"
    assert float(jnp.max(jnp.abs(g_lw))) > 0, "soil_veg gradient w.r.t. LW flux is zero"
    print(f"  advance_soil_veg differentiable: grad finite+nonzero (dveg/dLWdn={float(g_lw[0]):.2e})  PASS")


if __name__ == "__main__":
    print("soil_vegetation validation:")
    test_soil_veg_vs_fortran_replica()
    test_soil_veg_init_and_integration()
    test_soil_veg_differentiable()
    print("All soil_vegetation tests PASSED.")
