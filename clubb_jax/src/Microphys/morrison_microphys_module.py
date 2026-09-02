"""JAX port of morrison_microphys_module.F90 — the CLUBB↔Morrison (M2005) microphysics interface driver.

Mirrors clubb_release/src/Microphys/morrison_microphys_module.F90: `morrison_microphys_driver` is the CLUBB-side
wrapper that runs the full single-column M2005 step (the upstream WRF scheme in
Microphys/Morrison_microphys/module_mp_graupel.py — `m2005_driver`), folds rain/ice/snow/cloud sedimentation
(`morrison_sedimentation`) into the hydrometeor tendencies, and returns the CLUBB-form `*_mc` tendencies.
morrison_microphys_step.py imports it, mirroring the Fortran `use morrison_microphys_module`.

Pure-jnp → differentiable. Validated by `tests/test_morrison_rates.py`.
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import Cp as _CP, Lv as _LV
from clubb_jax.src.Microphys.Morrison_microphys.module_mp_graupel import (
    m2005_driver, morrison_sedimentation, _sizefix_exp_number,
    _M_CONS_RAIN, _M_LAMMINR, _M_LAMMAXR,
)


def morrison_microphys_driver(rcm, Ncm, rrm, Nrm, rim, Nim, rsm, Nsm, rgm, Ngm,
                              thlm, rvm, T_in_K, exner, pres, rho, cf3d, dzq, dt, igraup=0):
    """JAX CLUBB-Morrison interface. Returns a dict of the CLUBB-form tendencies:
    {rcm_mc, rvm_mc, thlm_mc, Ncm_mc, rrm_mc, Nrm_mc, rim_mc, Nim_mc, rsm_mc, Nsm_mc, rgm_mc, Ngm_mc}.
    Inputs are the grid-mean CLUBB fields (rcm/Ncm cloud from the PDF; rrm…Ngm hydrometeors; thlm/rvm
    prognostics; T_in_K=absolute temperature). Sedimentation (rain+ice+snow) is folded into the
    hydrometeor *_mc via (field_final−field_initial)/dt; cloud-droplet sedimentation is the separate
    CLUBB cloud_drop_sed (not here)."""
    a = lambda x: jnp.asarray(x, dtype=jnp.float64)
    rcm, Ncm, rvm, thlm = a(rcm), a(Ncm), a(rvm), a(thlm)
    rrm, Nrm, rim, Nim = a(rrm), a(Nrm), a(rim), a(Nim)
    rsm, Nsm, rgm, Ngm = a(rsm), a(Nsm), a(rgm), a(Ngm)
    T_in_K, exner, pres, rho, dzq = a(T_in_K), a(exner), a(pres), a(rho), a(dzq)
    # process tendencies (grid-mean) from the full single-column driver
    ten = m2005_driver(rcm, Ncm, rrm, Nrm, rim, Nim, rsm, Nsm, rgm, Ngm,
                       rvm, T_in_K, pres, rho, cf3d, dt, igraup=igraup)
    # post-process fields (DUM = field + tendency·dt), then sedimentation on these. Cloud water is
    # sedimented too (its tendency folds into rcm_mc/Ncm_mc — QC3DTEN += QCSTEN, :4885).
    qr_p = rrm + ten['qr'] * dt; nr_p = Nrm + ten['nr'] * dt
    qi_p = rim + ten['qi'] * dt; ni_p = Nim + ten['ni'] * dt
    qs_p = rsm + ten['qni'] * dt; ns_p = Nsm + ten['ns'] * dt
    qc_p = rcm + ten['qc'] * dt; nc_p = Ncm + ten['nc'] * dt
    sed = morrison_sedimentation(qr_p, nr_p, qi_p, ni_p, qs_p, ns_p, rho, dzq, dt,
                                 qc=qc_p, nc=nc_p, T=T_in_K)
    clip = lambda x: jnp.maximum(x, 0.0)
    # final fields = post-process + sedimentation, clipped ≥ 0
    rrm_f = clip(qr_p + sed['qr'] * dt); Nrm_f = clip(nr_p + sed['nr'] * dt)
    # The slope clamps are applied PRE-RATE inside m2005_driver (faithful to F90:1881-2002 — affects the
    # rate inputs + all 5 species). The rain post-sed clamp is ALSO kept: the stored stats reflect the
    # Fortran's NEXT-step pre-rate clamp on the (unclamped) post-sed output, so a per-step driver must
    # clamp its OUTPUT to match (removing it regresses dycoms Nrm 1.0→1.15; Iter235/250).
    Nrm_f = _sizefix_exp_number(rrm_f, Nrm_f, _M_CONS_RAIN, _M_LAMMINR, _M_LAMMAXR)
    rim_f = clip(qi_p + sed['qi'] * dt); Nim_f = clip(ni_p + sed['ni'] * dt)
    rsm_f = clip(qs_p + sed['qs'] * dt); Nsm_f = clip(ns_p + sed['ns'] * dt)
    rgm_f = clip(rgm + ten['qg'] * dt); Ngm_f = clip(Ngm + ten['ng'] * dt)
    Ncm_f = clip(nc_p + sed['nc'] * dt)
    rcm_mc = ten['qc'] + sed['qc']   # rcm tendency = process (QC3DTEN, incl. PCC) + cloud sedimentation
    # M2005 integrates the fields at the end (QC3D+=QC3DTEN·dt, T3D+=T3DTEN·dt, :4911-4929), so rcm_r4 =
    # rcm + rcm_mc·dt and T_in_K = T + ten['T']·dt → thlm_mc = (ten['T'] − Lv/Cp·rcm_mc)/exner. The PCC
    # parts of ten['T'] and rcm_mc cancel (thlm conserved under condensation); the cloud-sed rcm change
    # remains — that is the strong WBF-like heating at cloud-top mixed-phase points (the 184-pt signal).
    # thlm_mc: the float64-exact form. The Fortran computes this through a single-precision
    # thlm<->T_in_K round-trip (morrison_microphys_module.F90:399/416/793, `real(...)`=REAL(4)), which
    # algebraically reduces to (ten['T'] − Lv/Cp·rcm_mc)/exner once the T_in_K_init and rcm terms cancel
    # (thlm is conserved under condensation; only the cloud-sed rcm change survives as the cloud-top heating).
    # The REFACTOR drops the deliberate `real*4` round-trip (its sole effect was a ~1e-7 single-precision
    # residual reproduced for bit-faithfulness, REFACTOR.md §3.1 A2) — float64 is simpler and more accurate;
    # the clear-air case (mpace_a) now correctly gives thlm_mc≈0 instead of the ~2.8e-7 artifact.
    thlm_mc = (ten['T'] - _LV / _CP * rcm_mc) / exner
    return {
        'rcm_mc': rcm_mc, 'rvm_mc': ten['qv'],
        'thlm_mc': thlm_mc,
        'Ncm_mc': (Ncm_f - Ncm) / dt,
        'rrm_mc': (rrm_f - rrm) / dt, 'Nrm_mc': (Nrm_f - Nrm) / dt,
        'rim_mc': (rim_f - rim) / dt, 'Nim_mc': (Nim_f - Nim) / dt,
        'rsm_mc': (rsm_f - rsm) / dt, 'Nsm_mc': (Nsm_f - Nsm) / dt,
        'rgm_mc': (rgm_f - rgm) / dt, 'Ngm_mc': (Ngm_f - Ngm) / dt,
    }
