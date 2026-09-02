#!/usr/bin/env python3
"""test_set_sfc_value_of_flux_profiles.py — validate advance_clubb_core_module.py:set_sfc_value_of_flux_profiles.

`set_sfc_value_of_flux_profiles` (advance_clubb_core_module.F90:1586) sets — or clears — the surface (`k_lb_zm`)
value of every turbulent flux profile. It is pure index-assignment logic (no numerics → no f2py oracle needed),
so the oracle here is a direct, independent transcription of the F90 source. The validated standalone cases all run
`l_host_applies_sfc_fluxes=False` (CLUBB applies the surface flux), so the `=True` WRF-host branch (zero the surface)
and the `l_linearize_pbl_winds` perturbation path are otherwise unexercised — this pins both. Checks each branch
combination AND that interior levels (k != k_lb) are left untouched. Pure-Python (no JAX physics run, no Fortran
oracle), so it never SKIPs. (iter 489)
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

from clubb_jax.src.CLUBB_core.advance_clubb_core_module import set_sfc_value_of_flux_profiles
from clubb_jax.src.CLUBB_core.grid_class import setup_grid


def _ref(k_lb, ngrdcol, nzm, sclr_dim, edsclr_dim, l_host, l_lin,
         sfc, prof):
    """Independent transcription of the F90 source. `sfc` / `prof` are dicts of the surface-value / profile arrays."""
    out = {k: np.array(v, dtype=np.float64) for k, v in prof.items()}
    if not l_host:
        for f in ("wpthlp", "wprtp", "upwp", "vpwp"):
            out[f][:, k_lb] = sfc[f]
        if l_lin:
            out["upwp_pert"][:, k_lb] = sfc["upwp_pert"]
            out["vpwp_pert"][:, k_lb] = sfc["vpwp_pert"]
        for s in range(sclr_dim):
            out["wpsclrp"][:, k_lb, s] = sfc["wpsclrp"][:, s]
        out["wpedsclrp"] = np.zeros((ngrdcol, nzm, max(edsclr_dim, 1)))
        for e in range(edsclr_dim):
            out["wpedsclrp"][:, k_lb, e] = sfc["wpedsclrp"][:, e]
    else:
        for f in ("wpthlp", "wprtp", "upwp", "vpwp"):
            out[f][:, k_lb] = 0.0
        for s in range(sclr_dim):
            out["wpsclrp"][:, k_lb, s] = 0.0
        out["wpedsclrp"] = np.zeros((ngrdcol, nzm, max(edsclr_dim, 1)))
    return out


def test_against_transcription():
    gr = setup_grid(ngrdcol=2, deltaz=40.0, zm_init=0.0, zm_top=400.0, grid_type=1)
    ng, nzm = gr.zm.shape
    k_lb = int(gr.k_lb_zm)
    rng = np.random.default_rng(7)
    combos = [
        dict(l_host=False, l_lin=False, sclr_dim=0, edsclr_dim=0),
        dict(l_host=False, l_lin=True,  sclr_dim=2, edsclr_dim=3),
        dict(l_host=True,  l_lin=False, sclr_dim=2, edsclr_dim=3),
        dict(l_host=True,  l_lin=True,  sclr_dim=0, edsclr_dim=0),
        dict(l_host=False, l_lin=False, sclr_dim=1, edsclr_dim=0),
    ]
    worst = 0.0
    for c in combos:
        sclr_dim, edsclr_dim = c["sclr_dim"], c["edsclr_dim"]
        # random nonzero interior profiles (so interior preservation is actually exercised)
        prof = dict(
            wpthlp=rng.standard_normal((ng, nzm)), wprtp=rng.standard_normal((ng, nzm)),
            upwp=rng.standard_normal((ng, nzm)),   vpwp=rng.standard_normal((ng, nzm)),
            upwp_pert=rng.standard_normal((ng, nzm)), vpwp_pert=rng.standard_normal((ng, nzm)),
            wpsclrp=rng.standard_normal((ng, nzm, sclr_dim)),
            wpedsclrp=rng.standard_normal((ng, nzm, max(edsclr_dim, 1))),
        )
        sfc = dict(
            wpthlp=rng.standard_normal(ng), wprtp=rng.standard_normal(ng),
            upwp=rng.standard_normal(ng), vpwp=rng.standard_normal(ng),
            upwp_pert=rng.standard_normal(ng), vpwp_pert=rng.standard_normal(ng),
            wpsclrp=rng.standard_normal((ng, sclr_dim)),
            wpedsclrp=rng.standard_normal((ng, edsclr_dim)),
        )
        ref = _ref(k_lb, ng, nzm, sclr_dim, edsclr_dim, c["l_host"], c["l_lin"], sfc, prof)
        jout = set_sfc_value_of_flux_profiles(
            sclr_dim, edsclr_dim, gr, c["l_host"], c["l_lin"],
            sfc["wpthlp"], sfc["wprtp"], sfc["upwp"], sfc["vpwp"], sfc["upwp_pert"], sfc["vpwp_pert"],
            sfc["wpsclrp"], sfc["wpedsclrp"],
            prof["wpthlp"], prof["wprtp"], prof["upwp"], prof["vpwp"],
            prof["upwp_pert"], prof["vpwp_pert"], prof["wpsclrp"], prof["wpedsclrp"],
            ng, nzm)
        names = ("wpthlp", "wprtp", "upwp", "vpwp", "upwp_pert", "vpwp_pert", "wpsclrp", "wpedsclrp")
        for nm, jv in zip(names, jout):
            rv = ref[nm]
            jv = np.asarray(jv)
            assert jv.shape == rv.shape, f"{nm} shape {jv.shape} != ref {rv.shape} (combo {c})"
            worst = max(worst, float(np.max(np.abs(jv - rv))) if jv.size else 0.0)
        # explicit interior-preservation guard for the scalar flux profiles
        interior = [k for k in range(nzm) if k != k_lb]
        assert np.array_equal(np.asarray(jout[0])[:, interior], prof["wpthlp"][:, interior]), \
            "interior levels of wpthlp were modified"
    assert worst < 1e-15, f"set_sfc_value_of_flux_profiles vs transcription mismatch {worst:.2e}"
    print(f"  set_sfc_value_of_flux_profiles vs F90 transcription (5 branch combos, interior preserved): "
          f"exact, worst {worst:.2e}  PASS")


def main():
    print("test_set_sfc_value_of_flux_profiles:")
    test_against_transcription()
    print("All set_sfc_value_of_flux_profiles checks PASSED")


if __name__ == "__main__":
    main()
