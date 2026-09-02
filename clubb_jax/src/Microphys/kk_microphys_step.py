"""Per-step KK (Khairoutdinov-Kogan) microphysics orchestration (Iter155).

Composes the oracle-validated KK pieces into the per-step call that the Fortran driver makes
(clubb_driver.F90: pdf_hydromet_microphys_prep -> calc_microphys_scheme_tendcies -> advance_microphys).
Mirrors `_cloud_drop_sed` in advance_clubb_to_end.py: extract from `state`, run the physics, store the
tendencies for the next step's forcings.

TWO-STAGE STRUCTURE (incremental-replacement / shadow-comparison discipline, see DESIGN.md):
  - The TENDENCY computation (precip_fraction -> compute_kk_microphysics) runs on LIVE state and STORES the
    result, validated against the rf02_do/rico oracles. The hydrometeor TRANSPORT (advance_one_hydrometeor +
    fill_holes, which needs the sed velocities from `prereqs`) and the feedback APPLICATION are gated behind
    `state['l_kk_micro_apply']`. That flag's model default is False, but `clubb_driver` sets it TRUE for every
    KK case (`l_kk_micro_apply=(microphys_scheme=='khairoutdinov_kogan')`, clubb_driver.py:1069) — so the
    transport+application stage IS live for the KK cases (rico), which is consequently FP-limited at precip
    onset and classified BLOCKED in compare_cases (DESIGN.md "rico"). Non-KK cases never reach this code
    (hydromet_dim==0 / microphys_scheme!='khairoutdinov_kogan'), so they are byte-for-byte unaffected.
"""
import numpy as np
import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.tracer_numpy import _is_tracer_arg  # REFACTOR B5: detect a jax.grad trace
from clubb_jax.src.CLUBB_core.precipitation_fraction import precip_fraction
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.Microphys.KK_microphys.kk_microphys_driver import compute_kk_microphysics

_UPSILON_IDX = 64   # clubb_params index of upsilon_precip_frac_rat (numerical_check._PNAME_IDX)


def _covar_driver_jit():
    """jit-compiled KK_upscaled_covar_driver (the eager 16-branch dispatch + parabolic-cylinder Dv is
    ~100s/step on a 160-level grid; jit fuses it). Cached on first use."""
    global _COVAR_DRIVER_JIT
    try:
        return _COVAR_DRIVER_JIT
    except NameError:
        from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_covariances import KK_upscaled_covar_driver
        _COVAR_DRIVER_JIT = jax.jit(KK_upscaled_covar_driver)
        return _COVAR_DRIVER_JIT


def _compute_kk_covar_mc(state, pdf, prereqs, precip_frac_1, precip_frac_2):
    """Compute the 5 second-moment microphysics tendencies (wprtp_mc/wpthlp_mc/rtp2_mc/thlp2_mc/
    rtpthlp_mc) via the validated `KK_upscaled_covar_driver`, on zt levels. Sources the ~70 component
    moments from `pdf` (the post-advance ADG1 PDF) + `prereqs` (rr/Nr moments, process tendencies +
    coefficients, Ncnm) + the prescribed normal-space correlations. For ADG1 corr_w_chi/corr_w_eta = 0
    (pdf_closure:1037); rico/dycoms use constant Nc so sigma_Ncn = 0. Returns 5 zt-level arrays."""
    from clubb_jax.src.CLUBB_core.corr_varnce_module import kk_prescribed_correlations
    KK_upscaled_covar_driver = _covar_driver_jit()
    J = jnp.asarray
    pc = kk_prescribed_correlations()
    one = J(np.ones_like(np.asarray(pdf.chi_1, np.float64)))
    z = 0.0 * one

    def f(a):           # pdf field -> jnp zt array
        return J(np.asarray(a, np.float64))

    Ncnm = prereqs['Ncnm']
    Ncnm_n = jnp.log(jnp.maximum(Ncnm, 1e-30))
    sw1 = jnp.sqrt(jnp.maximum(f(pdf.varnce_w_1), 0.0))
    sw2 = jnp.sqrt(jnp.maximum(f(pdf.varnce_w_2), 0.0))
    return KK_upscaled_covar_driver(
        J(np.asarray(state['wm_zt'], np.float64)), J(np.asarray(state['rtm'], np.float64)),
        J(np.asarray(state['thlm'], np.float64)), J(np.asarray(state['exner'], np.float64)),
        f(pdf.w_1), f(pdf.w_2), f(pdf.chi_1), f(pdf.chi_2), z, z,                 # mu_w, mu_chi, mu_eta(=0)
        prereqs['mu_rr_1'], prereqs['mu_rr_2'], prereqs['mu_Nr_1'], prereqs['mu_Nr_2'], Ncnm, Ncnm,
        prereqs['mu_rr_1_n'], prereqs['mu_rr_2_n'], prereqs['mu_Nr_1_n'], prereqs['mu_Nr_2_n'], Ncnm_n, Ncnm_n,
        sw1, sw2, f(pdf.stdev_chi_1), f(pdf.stdev_chi_2), f(pdf.stdev_eta_1), f(pdf.stdev_eta_2),
        prereqs['sigma_rr_1'], prereqs['sigma_rr_2'], prereqs['sigma_Nr_1'], prereqs['sigma_Nr_2'], z, z,
        prereqs['sigma_rr_1_n'], prereqs['sigma_rr_2_n'], prereqs['sigma_Nr_1_n'], prereqs['sigma_Nr_2_n'], z, z,
        z, z, pc['corr_w_rr'] * one, pc['corr_w_rr'] * one, pc['corr_w_Nr'] * one, pc['corr_w_Nr'] * one,
        pc['corr_w_Ncn'] * one, pc['corr_w_Ncn'] * one,                          # corr_w_chi=0, corr_w_*_n
        f(pdf.corr_chi_eta_1), f(pdf.corr_chi_eta_2),
        pc['corr_chi_rr'] * one, pc['corr_chi_rr'] * one, pc['corr_chi_Nr'] * one, pc['corr_chi_Nr'] * one,
        pc['corr_chi_Ncn'] * one, pc['corr_chi_Ncn'] * one,
        pc['corr_eta_rr'] * one, pc['corr_eta_rr'] * one, pc['corr_eta_Nr'] * one, pc['corr_eta_Nr'] * one,
        pc['corr_eta_Ncn'] * one, pc['corr_eta_Ncn'] * one, pc['corr_rr_Nr'] * one, pc['corr_rr_Nr'] * one,
        f(pdf.mixt_frac), J(np.asarray(precip_frac_1, np.float64)), J(np.asarray(precip_frac_2, np.float64)),
        prereqs['coef_evap'], prereqs['coef_auto'], prereqs['coef_accr'],
        prereqs['evap'], prereqs['auto'], prereqs['accr'],
        f(pdf.rt_1), f(pdf.rt_2), f(pdf.thl_1), f(pdf.thl_2),
        f(pdf.crt_1), f(pdf.crt_2), f(pdf.cthl_1), f(pdf.cthl_2))


def advance_kk_microphysics(state: dict):
    """Compute the KK microphysics tendencies (rcm_mc/thlm_mc/rrm_mc/Nrm_mc) from the post-advance
    state and store them. Transport + feedback application are gated (see module docstring)."""
    # Detach-under-trace (REFACTOR B5): KK microphysics runs AFTER the core, so the *_mc tendencies it stores
    # feed only the NEXT step's forcings — dead for a single-step gradient. Skip under a jax.grad trace
    # (exact for single-step; a detached forcing for multi-step rollouts), same rationale as BUGSrad radiation.
    if _is_tracer_arg([state['thlm'], state['rcm'], state['cloud_frac']]):
        return

    hmm = state['hm_metadata']
    pdf = state['pdf_params']
    iirr, iiNr = int(hmm.iirr), int(hmm.iiNr)
    hydromet = np.asarray(state['hydromet'], np.float64)          # (ngrdcol, nzt, hydromet_dim)
    rrm = hydromet[..., iirr]
    Nrm = hydromet[..., iiNr]

    cloud_frac = np.asarray(state['cloud_frac'], np.float64)
    isf = np.asarray(state['ice_supersat_frac'], np.float64)
    cf1, cf2 = np.asarray(pdf.cloud_frac_1, np.float64), np.asarray(pdf.cloud_frac_2, np.float64)
    isf1 = np.asarray(pdf.ice_supersat_frac_1, np.float64)
    isf2 = np.asarray(pdf.ice_supersat_frac_2, np.float64)
    mixt_frac = np.asarray(pdf.mixt_frac, np.float64)
    l_mix = np.asarray(state['l_mix_rat_hm'], bool)
    l_frozen = np.asarray(hmm.l_frozen_hm, bool)
    hm_tol = np.asarray(hmm.hydromet_tol, np.float64)
    upsilon = float(np.asarray(state['clubb_params']).reshape(-1)[_UPSILON_IDX]) \
        if np.asarray(state['clubb_params']).ndim == 1 else \
        float(np.asarray(state['clubb_params'])[0, _UPSILON_IDX])

    gr = state['gr']
    clubb_params = jnp.asarray(np.asarray(state['clubb_params'], np.float64))
    if clubb_params.ndim == 1:
        clubb_params = clubb_params[None, :]
    stats = JaxStats.empty(
        l_sample=False,
        names=(),
        ncol=gr.ngrdcol,
        max_nlev=max(gr.nzm, gr.nzt),
    )
    _err, precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol, _stats = precip_fraction(
        gr,
        gr.nzt,
        gr.ngrdcol,
        hydromet.shape[-1],
        hydromet,
        cloud_frac,
        cf1,
        l_mix,
        l_frozen,
        hm_tol,
        cf2,
        isf,
        isf1,
        isf2,
        mixt_frac,
        clubb_params,
        ErrInfo.initialized(gr.ngrdcol),
        stats,
    )
    precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol = (
        np.asarray(x) for x in (precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol))

    rho = np.asarray(state['rho'], np.float64)
    exner = np.asarray(state['exner'], np.float64)
    T_liq = np.asarray(state['thlm'], np.float64) * exner
    p_in_Pa = np.asarray(state['p_in_Pa'], np.float64)
    rcm = np.asarray(state['rcm'], np.float64)
    Nc_in_cloud = np.asarray(state['Nc_in_cloud'], np.float64)
    dt = float(state['dt_main'])
    rr_ratio = float(getattr(hmm, 'rr_ratio', 1.25))
    Nr_ratio = float(getattr(hmm, 'Nr_ratio', 1.25))

    (rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc), prereqs = compute_kk_microphysics(
        rrm, Nrm, np.asarray(pdf.chi_1), np.asarray(pdf.chi_2),
        np.asarray(pdf.stdev_chi_1), np.asarray(pdf.stdev_chi_2), mixt_frac,
        np.asarray(pdf.thl_1), np.asarray(pdf.thl_2), Nc_in_cloud, cf1, cf2,
        precip_frac, precip_frac_1, precip_frac_2, precip_frac_tol,
        rho, T_liq, p_in_Pa, exner, rcm, dt,
        rr_ratio=rr_ratio, Nr_ratio=Nr_ratio, l_return_vel_prereqs=True)

    state['_kk_rcm_mc'] = np.asarray(rcm_mc, np.float64)
    state['_kk_thlm_mc'] = np.asarray(thlm_mc, np.float64)
    state['_kk_rrm_mc'] = np.asarray(rrm_mc, np.float64)
    state['_kk_Nrm_mc'] = np.asarray(Nrm_mc, np.float64)
    state['_kk_precip_frac'] = precip_frac
    state['_kk_prereqs'] = prereqs

    # Second-moment microphysics source (KK_upscaled_covar_driver, Iter171-validated). Computed on zt,
    # then zt2zm with the boundary momentum levels zeroed (KK_microphys_module.F90:901-916), and stored
    # for the next step's second-moment forcings (l_var_covar_src; applied in advance_clubb_to_end.py).
    if state.get('l_var_covar_src', True):
        mc_zt = _compute_kk_covar_mc(state, pdf, prereqs, precip_frac_1, precip_frac_2)
        for nm, v_zt in zip(('wprtp_mc', 'wpthlp_mc', 'rtp2_mc', 'thlp2_mc', 'rtpthlp_mc'), mc_zt):
            v_zm = np.array(zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, v_zt), np.float64)
            v_zm[:, 0] = 0.0
            v_zm[:, -1] = 0.0
            state['_kk_' + nm] = v_zm

    if not state.get('l_kk_micro_apply', False):
        return

    # ---- Transport + application stage -----------------------------------------------
    # Per the Fortran advance_microphys (advance_microphys_module.F90): each hydrometeor is advanced by
    # the implicit transport solve (eddy diffusion K_hm + nu_hm, mean advection, mean+turbulent
    # sedimentation), seeded by the microphysics source tendency; then fill_holes; rcm_mc/thlm_mc feed
    # the next step's forcings (already wired in advance_clubb_to_end.py:67-68).
    from clubb_jax.src.Microphys.KK_microphys_module import KK_sedimentation
    from clubb_jax.src.CLUBB_core.setup_clubb_pdf_params import hydrometp2_zt
    from clubb_jax.src.Microphys.KK_microphys.KK_upscaled_turbulent_sed import KK_sed_vel_covars
    from clubb_jax.src.Microphys.advance_microphys_module import calculate_K_hm, advance_one_hydrometeor
    from clubb_jax.src.CLUBB_core.Skx_module import Skx_func
    from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_vertical, fill_holes_hydromet_clip_jax

    clubb_params = jnp.asarray(np.asarray(state['clubb_params'], np.float64))
    w_tol = float(state.get('w_tol', 2.0e-3))
    wp2 = jnp.asarray(np.asarray(state['wp2'], np.float64))            # (ng, nzm)
    wp3_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, jnp.asarray(np.asarray(state['wp3'], np.float64)))
    Skw_zm = Skx_func(wp2, wp3_zm, w_tol, clubb_params)            # (ng, nzm)
    Kh_zm = jnp.asarray(np.asarray(state['Kh_zm'], np.float64))
    wm_zt = jnp.asarray(np.asarray(state['wm_zt'], np.float64))
    rho_ds_zm = jnp.asarray(np.asarray(state['rho_ds_zm'], np.float64))
    invrs_rho_ds_zt = jnp.asarray(np.asarray(state['invrs_rho_ds_zt'], np.float64))
    w_above = jnp.asarray(np.asarray(gr.weights_zt2zm, np.float64)[:, :, 0])   # T_ABOVE
    nu_hm = float(np.asarray(state['nu_vert_res_dep'].nu_hm).reshape(-1)[0])
    hm_tol = float(np.asarray(hmm.hydromet_tol).reshape(-1)[iirr])

    # Mean sedimentation velocities (zt) from the rain-drop mean volume radius; to momentum (zt2zm).
    Vrr_zt, VNr_zt = KK_sedimentation(prereqs['mvr'], l_clip_positive_sed=True)
    sed = KK_sed_vel_covars(
        precip_frac_1 * prereqs['mu_rr_1'], precip_frac_2 * prereqs['mu_rr_2'],
        precip_frac_1 * prereqs['mu_Nr_1'], precip_frac_2 * prereqs['mu_Nr_2'], prereqs['mvr'],
        prereqs['mu_rr_1'], prereqs['mu_rr_2'], prereqs['mu_Nr_1'], prereqs['mu_Nr_2'],
        prereqs['mu_rr_1_n'], prereqs['mu_rr_2_n'], prereqs['mu_Nr_1_n'], prereqs['mu_Nr_2_n'],
        prereqs['sigma_rr_1'], prereqs['sigma_rr_2'], prereqs['sigma_Nr_1'], prereqs['sigma_Nr_2'],
        prereqs['sigma_rr_1_n'], prereqs['sigma_rr_2_n'], prereqs['sigma_Nr_1_n'], prereqs['sigma_Nr_2_n'],
        prereqs['corr_rr_Nr_n'], prereqs['corr_rr_Nr_n'], jnp.asarray(mixt_frac))

    def _advance_hm(hm, hm_tndcy, ratio, V_zt, Vi_zt, Ve_zt):
        hmp2_zm = zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, hydrometp2_zt(hm, jnp.asarray(precip_frac), ratio))
        K_hm = calculate_K_hm(wp2, Kh_zm, Skw_zm, hm, hmp2_zm, hm_tol, gr)
        out = advance_one_hydrometeor(
            float(dt), hm, hm_tndcy, K_hm, nu_hm, wm_zt,
            zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, V_zt),
            zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Vi_zt),
            zt2zm(gr.nzm, gr.nzt, gr.ngrdcol, gr, Ve_zt),
            rho_ds_zm, invrs_rho_ds_zt, gr, w_above)
        # fill_holes_vertical: faithful to fill_holes_driver_api (fill_holes.F90) — threshold is
        # zero_threshold (0), NOT hydromet_tol (the JAX previously over-filled tiny-positive edge values
        # below tol, diverging rrm at the cloud edge); full zt range begin=1..nzt (0-based 0..nzt-1). With
        # threshold 0 the call is a no-op when there are no negatives (matching the Fortran's `any(<0)` guard).
        dz = np.asarray(gr.dzt, np.float64)
        out = fill_holes_vertical(out, np.asarray(state['rho_ds_zt'], np.float64), dz,
                                      0.0, 0, np.asarray(hm).shape[-1] - 1,
                                      int(state['flags'].fill_holes_type))
        return np.asarray(out, np.float64)

    rrm_new = _advance_hm(jnp.asarray(rrm), jnp.asarray(rrm_mc), rr_ratio,
                          Vrr_zt, sed['Vrrprrp_impc'], sed['Vrrprrp_expc'])
    Nrm_new = _advance_hm(jnp.asarray(Nrm), jnp.asarray(Nrm_mc), Nr_ratio,
                          VNr_zt, sed['VNrpNrp_impc'], sed['VNrpNrp_expc'])

    # fill_holes_driver_api ALSO clips hydromet <= hydromet_tol to 0 (fill_holes.F90:2444-2476), returning the
    # removed mass to vapor + a latent thlm adjustment (number zeroed by clip_hydromet_conc_mvr). The JAX did only
    # the zero_threshold hole-fill, so sub-tol rrm ACCUMULATED → premature precip onset (the rico seed, Iter306).
    # Apply the tol clip on the post-advance rain — now via the fill_holes.py mirror (mirror-refactor iter 214).
    _rr_tol = float(np.asarray(hmm.hydromet_tol).reshape(-1)[iirr])
    rrm_new, Nrm_new, rvm_mc, thlm_mc = fill_holes_hydromet_clip_jax(
        rrm_new, Nrm_new, _rr_tol, rvm_mc, thlm_mc, state['exner'], float(dt))

    hydromet = hydromet.copy()
    hydromet[..., iirr] = rrm_new
    hydromet[..., iiNr] = Nrm_new
    state['hydromet'] = hydromet
    # the rt microphysics tendency is rcm_mc + rvm_mc (cloud→rain + rain→vapor), per
    # clubb_driver.F90:3337 `rtm_forcing += rcm_mc + rvm_mc` (the Morrison JAX path already does this; the KK
    # path was DROPPING rvm_mc — the rain-evaporation→vapor term, =0 until rain forms ~step 5, which is exactly
    # where rico started diverging). advance_clubb_to_end adds state['rcm_mc'] to rtm_forcing, so fold rvm_mc in.
    state['rcm_mc'] = np.asarray(rcm_mc, np.float64) + np.asarray(rvm_mc, np.float64)
    state['thlm_mc'] = np.asarray(thlm_mc, np.float64)
