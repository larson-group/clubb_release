"""Per-step Morrison microphysics wiring for the JAX CLUBB loop.

Mirrors the Fortran call sequence (clubb_driver.F90 → morrison_microphys_driver): each step computes
the CLUBB-form microphysics tendencies (rcm_mc/rvm_mc/thlm_mc + the hydrometeor *_mc) from the
post-advance state via the ported `morrison_microphys_driver`, stores rcm_mc/thlm_mc for the next
step's rtm/thlm forcings (applied in advance_clubb_to_end.py), and advances the hydrometeor fields.

The hydrometeor *_mc already include sedimentation (folded in as (field_final-field_initial)/dt, the
Iter204 form), so the Morrison sedimentation velocity for the CLUBB transport is zero.
"""
import numpy as np

from clubb_jax.src.CLUBB_core.tracer_numpy import _is_tracer_arg  # REFACTOR B5: detect a jax.grad trace
from clubb_jax.src.Microphys.morrison_microphys_module import morrison_microphys_driver


def advance_morrison_microphysics(state: dict):
    """Compute the Morrison microphysics tendencies from the post-advance state and store them.

    Stores `_morr_rcm_mc`, `_morr_thlm_mc` (→ next step's forcings), and updates the hydrometeor
    array + Ncm by the microphysics tendencies (which include sedimentation). First-pass field update
    is an explicit Euler integration by the *_mc; the full CLUBB hydrometeor transport (advection +
    eddy diffusion, with zero sedimentation velocity) is layered on subsequently.
    """
    # Detach-under-trace (REFACTOR B5): Morrison runs AFTER the core, storing _morr_*_mc for the NEXT step's
    # forcings — dead for a single-step gradient. Skip under a jax.grad trace (exact for single-step; detached
    # forcing for multi-step), same as KK microphysics / BUGSrad radiation.
    if _is_tracer_arg([state['thlm'], state['rcm'], state['cloud_frac']]):
        return

    hmm = state['hm_metadata']
    g = lambda k: np.asarray(state[k], np.float64)
    hydromet = g('hydromet')                      # (ngrdcol, nzt, 8)
    iirr, iiNr, iiri, iiNi = int(hmm.iirr), int(hmm.iiNr), int(hmm.iiri), int(hmm.iiNi)
    iirs, iiNs, iirg, iiNg = int(hmm.iirs), int(hmm.iiNs), int(hmm.iirg), int(hmm.iiNg)

    rcm = g('rcm'); thlm = g('thlm')
    # grid-mean cloud droplet number is DIAGNOSED as Nc_in_cloud · cloud_frac (the oracle's Ncm
    # matches this exactly, ratio 1.0); the JAX state['Ncm'] is stale (0) at this point in the loop.
    Ncm = g('Nc_in_cloud') * g('cloud_frac')
    rvm = g('rtm') - rcm                            # vapor = total water − cloud water
    exner = g('exner'); rho = g('rho'); pres = g('p_in_Pa'); cf = g('cloud_frac')
    T_in_K = thlm * exner + (_LV() / _CP()) * rcm   # thlm2T_in_K inverse of (T-Lv/Cp·rcm)/exner
    dt = float(state['dt_main'])
    gr = state['gr']
    # Morrison's dzq = delta_zt(k) = zt(k+1)−zt(k) = dzm(k+1) (the MOMENTUM-level spacing, NOT dzt;
    # microphys_driver.F90:413-419). On a stretched grid (nov11 is grid_type=2) dzm≠dzt → this matters.
    dzq = np.asarray(gr.dzm, np.float64)[:, 1:]   # (ngrdcol, nzt)

    pick = lambda i: hydromet[..., i]
    out = morrison_microphys_driver(
        rcm, Ncm, pick(iirr), pick(iiNr), pick(iiri), pick(iiNi), pick(iirs), pick(iiNs),
        pick(iirg), pick(iiNg), thlm, rvm, T_in_K, exner, pres, rho, cf, dzq, dt)
    out = {k: np.asarray(v, np.float64) for k, v in out.items()}

    # store the mean-field tendencies for the next step's forcings (rtm += rcm_mc, thlm += thlm_mc)
    state['_morr_rcm_mc'] = out['rcm_mc']
    state['_morr_thlm_mc'] = out['thlm_mc']
    state['_morr_rvm_mc'] = out['rvm_mc']

    # Advance each hydrometeor by the CLUBB transport — the implicit solve (mean advection wm_zt +
    # eddy diffusion K_hm+nu_hm), seeded by the microphysics+sedimentation source *_mc (Fortran
    # advance_microphys; the same `advance_one_hydrometeor` the KK path uses). The sedimentation is
    # already folded into *_mc so the solve runs with l_sed=False (zero sed velocity). The eddy
    # diffusivity is K_hm + nu_hm, with the VARIABLE K_hm = calculate_K_hm(.., hydrometp2). The overall
    # hydrometeor variance is hydrometp2 = ((ratio+1)/precip_frac − 1)·hm² (setup_clubb_pdf_params:449);
    # the bulk scheme has precip_frac≤0 → pf=1 → hydrometp2 = ratio·hm² (ratio=1.25 for rr/Nr,
    # corr_varnce slope=0/intrcpt=1.25), so K_hm is LARGE (NOT 0 — the Iter239 K_hm=0 was wrong; the
    # stored rrp2=0 ≠ this internal hydrometp2). This is what sustains the cloud-top rain.
    import jax.numpy as jnp
    from clubb_jax.src.Microphys.advance_microphys_module import advance_one_hydrometeor, calculate_K_hm
    from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_vertical
    from clubb_jax.src.CLUBB_core.grid_class import zt2zm
    from clubb_jax.src.CLUBB_core.Skx_module import Skx_func
    wm_zt = jnp.asarray(g('wm_zt'))
    rho_ds_zm = jnp.asarray(g('rho_ds_zm'))
    invrs_rho_ds_zt = jnp.asarray(g('invrs_rho_ds_zt'))
    w_above = jnp.asarray(np.asarray(gr.weights_zt2zm, np.float64)[:, :, 0])
    nu_hm = float(np.asarray(state['nu_vert_res_dep'].nu_hm).reshape(-1)[0])
    ngrdcol, nzm = rho_ds_zm.shape
    nzt = hydromet.shape[1]
    zero_zm = jnp.zeros((ngrdcol, nzm))               # zero sed velocities (sed is in *_mc)
    rho_ds_zt = np.asarray(state['rho_ds_zt'], np.float64)
    dzt = np.asarray(gr.dzt, np.float64)
    fh_type = int(state['flags'].fill_holes_type)
    # K_hm inputs (momentum levels): wp2, Kh_zm, Skw_zm.
    clubb_params = jnp.asarray(g('clubb_params'))
    w_tol = float(state.get('w_tol', 2.0e-3))
    wp2 = jnp.asarray(g('wp2'))
    Skw_zm = Skx_func(wp2, zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, jnp.asarray(g('wp3'))), w_tol, clubb_params)
    Kh_zm = jnp.asarray(g('Kh_zm'))
    hm_tol = np.asarray(hmm.hydromet_tol).reshape(-1)
    # ★ The overall hydrometeor variance is hydrometp2 = ((ratio_ip+1)/precip_frac − 1)·hm²
    # (setup_clubb_pdf_params:449). precip_frac is NOT 1 — it is the precipitation fraction, set by the
    # max-in-precip-mean limiter to ~hm/(mixt_frac·MAX_HM_IP_COMP_MEAN) ≈ 0.3 for dycoms — so hydrometp2
    # ≈ 5·hm² (not 1.25·hm²) and K_hm is ~2× larger than precip_frac=1 gives. Compute the real
    # precip_frac (the same precip_fraction the KK path uses) from the PDF component fields.
    from clubb_jax.src.CLUBB_core.precipitation_fraction import precip_fraction
    from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import hydrometp2_zt
    from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
    from clubb_jax.src.CLUBB_core.err_info import ErrInfo
    pdf = state['pdf_params']
    ga = lambda a: np.asarray(getattr(pdf, a), np.float64)
    stats = JaxStats.empty(
        l_sample=False,
        names=(),
        ncol=gr.ngrdcol,
        max_nlev=max(gr.nzm, gr.nzt),
    )
    _err, precip_frac, _pf1, _pf2, _pftol, _stats = precip_fraction(
        gr,
        gr.nzt,
        gr.ngrdcol,
        hydromet.shape[-1],
        hydromet,
        g('cloud_frac'),
        ga('cloud_frac_1'),
        hmm.l_mix_rat_hm,
        hmm.l_frozen_hm,
        hmm.hydromet_tol,
        ga('cloud_frac_2'),
        g('ice_supersat_frac'),
        ga('ice_supersat_frac_1'),
        ga('ice_supersat_frac_2'),
        ga('mixt_frac'),
        jnp.asarray(g('clubb_params')),
        ErrInfo.initialized(gr.ngrdcol),
        stats,
    )
    precip_frac = jnp.asarray(np.asarray(precip_frac, np.float64))
    ratio_ip = np.asarray(hmm.hmp2_ip_on_hmm2_ip, np.float64).reshape(-1)

    new_hm = hydromet.copy()
    for idx, key in ((iirr, 'rrm_mc'), (iiNr, 'Nrm_mc'), (iiri, 'rim_mc'), (iiNi, 'Nim_mc'),
                     (iirs, 'rsm_mc'), (iiNs, 'Nsm_mc'), (iirg, 'rgm_mc'), (iiNg, 'Ngm_mc')):
        hm = jnp.asarray(hydromet[..., idx])
        # hydrometp2 on zt = ((ratio_ip+1)/precip_frac − 1)·hm² where hm >= tol, else 0 (the per-species
        # guard, setup_clubb_pdf_params.F90:446-455), then interpolated to zm. ratio_ip = hmp2_ip_on_hmm2_ip.
        hmp2_zt = jnp.where(hm >= float(hm_tol[idx]),
                            hydrometp2_zt(hm, precip_frac, float(ratio_ip[idx])), 0.0)
        hmp2_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hmp2_zt)
        K_hm = calculate_K_hm(wp2, Kh_zm, Skw_zm, hm, hmp2_zm, float(hm_tol[idx]), gr)
        adv = advance_one_hydrometeor(
            dt, hm, jnp.asarray(out[key]), K_hm, nu_hm, wm_zt,
            zero_zm, zero_zm, zero_zm, rho_ds_zm, invrs_rho_ds_zt, gr, w_above, l_sed=False)
        adv = fill_holes_vertical(adv, rho_ds_zt, dzt, 0.0, 0, nzt - 1, fh_type)
        new_hm[..., idx] = np.maximum(np.asarray(adv), 0.0)
    state['hydromet'] = new_hm
    state['Ncm'] = np.maximum(Ncm + out['Ncm_mc'] * dt, 0.0)


def _LV():
    from clubb_jax.src.CLUBB_core.constants_clubb import Lv
    return float(Lv)


def _CP():
    from clubb_jax.src.CLUBB_core.constants_clubb import Cp
    return float(Cp)
