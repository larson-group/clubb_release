"""JAX port of ``src/CLUBB_core/stats_clubb_utilities.F90``.

Porting deviations:
- Only ``stats_accumulate`` is ported here.  The Fortran
  ``stats_accumulate_hydromet_api`` and ``stats_accumulate_lh_tend`` entry
  points are not present because the current JAX core does not call them.
- Fortran ``stats_type`` mutation through ``stats_update`` is represented by
  functional updates to the JAX stats state.
- OpenACC host/device synchronization comments are omitted because JAX handles
  device placement through array tracing rather than explicit ``!$acc update``
  directives.
- Small per-column stats reductions remain Python loops over static grid
  columns; field math is expressed with JAX arrays.
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.T_in_K_module import thlm2T_in_K
from clubb_jax.src.CLUBB_core.Skx_module import Skx_func, compute_gamma_Skw
from clubb_jax.src.CLUBB_core.advance_helper_module import (
    calc_wp3_on_wp2,
    vertical_avg,
    vertical_integral,
)
from clubb_jax.src.CLUBB_core.constants_clubb import (
    cloud_frac_min,
    eps,
    rc_tol,
    w_tol,
)
from clubb_jax.src.CLUBB_core.model_flags import l_gamma_Skw
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_ice


def _calculate_spurious_source(
    integral_after,
    integral_before,
    flux_top,
    flux_sfc,
    integral_forcing,
    dt,
):
    """JAX equivalent of ``calculate_spurious_source`` from numerical_check."""
    return (integral_after - integral_before) / dt + flux_top - flux_sfc - integral_forcing


def stats_accumulate(
    nzm,
    nzt,
    ngrdcol,
    sclr_dim,
    edsclr_dim,
    gr,
    dt,
    l_implemented,
    l_host_applies_sfc_fluxes,
    l_stability_correct_tau_zm,
    clubb_params,
    um,
    vm,
    upwp,
    vpwp,
    up2,
    vp2,
    thlm,
    rtm,
    thlm_before,
    rtm_before,
    thlm_forcing,
    rtm_forcing,
    wpthlp_sfc,
    wprtp_sfc,
    wprtp,
    wpthlp,
    wp2,
    wp3,
    rtp2,
    rtp3,
    thlp2,
    thlp3,
    rtpthlp,
    p_in_Pa,
    exner,
    rho,
    rho_zm,
    rho_ds_zm,
    rho_ds_zt,
    thv_ds_zm,
    thv_ds_zt,
    wm_zt,
    wm_zm,
    rcm,
    cloud_frac,
    thvm,
    ug,
    vg,
    ddzt_umvm_sqd,
    stability_correction,
    Kh_zt,
    rsat,
    Kh_zm,
    em,
    sclrm,
    sclrp2,
    sclrprtp,
    sclrpthlp,
    sclrm_forcing,
    wpsclrp,
    wpedsclrp,
    edsclrm,
    edsclrm_forcing,
    saturation_formula,
    stats,
):
    """Accumulate per-timestep CLUBB statistics into ``JaxStats`` banks."""
    # ---- Begin Code ----

    # Sample fields
    if not stats.l_sample:
        return stats

    dzm = gr.grid_dir * gr.dzm
    dzt = gr.grid_dir * gr.dzt

    wp3_on_wp2, wp3_on_wp2_zt = calc_wp3_on_wp2(
        nzm, nzt, ngrdcol, gr, wp2, wp3,
    )

    if stats.var_on_stats_list("T_in_K") or stats.var_on_stats_list("rsati"):
        T_in_K = thlm2T_in_K(thlm, exner, rcm)
        stats = stats.update("T_in_K", T_in_K)
    if stats.var_on_stats_list("rsati"):
        rsati = sat_mixrat_ice(p_in_Pa, T_in_K)
        stats = stats.update("rsati", rsati)

    if stats.var_on_stats_list("rcm_in_cloud"):
        cloud_frac_safe = jnp.where(cloud_frac > cloud_frac_min, cloud_frac, 1.0)
        rcm_in_cloud = jnp.where(
            cloud_frac > cloud_frac_min,
            rcm / cloud_frac_safe,
            rcm,
        )
        stats = stats.update("rcm_in_cloud", rcm_in_cloud)

    if stats.var_on_stats_list("shear"):
        shear = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
        shear_int = (
            -upwp[:, 1:nzm - 1]
            * (um[:, 1:nzm - 1] - um[:, :nzm - 2])
            * gr.invrs_dzm[:, 1:nzm - 1]
            - vpwp[:, 1:nzm - 1]
            * (vm[:, 1:nzm - 1] - vm[:, :nzm - 2])
            * gr.invrs_dzm[:, 1:nzm - 1]
        )
        shear = shear.at[:, 1:nzm - 1].set(shear_int)
        stats = stats.update("shear", shear)

    for name, value in (
        # stats_zt variables
        ("thlm", thlm),
        ("thvm", thvm),
        ("rtm", rtm),
        ("rcm", rcm),
        ("um", um),
        ("vm", vm),
        ("wm_zt", wm_zt),
        ("ug", ug),
        ("vg", vg),
        ("p_in_Pa", p_in_Pa),
        ("exner", exner),
        ("rho_ds_zt", rho_ds_zt),
        ("thv_ds_zt", thv_ds_zt),
        ("wp3", wp3),
        ("Kh_zt", Kh_zt),
        ("rho", rho),
        ("rsat", rsat),
        ("thlp3", thlp3),
        ("rtp3", rtp3),
        ("wp3_on_wp2_zt", wp3_on_wp2_zt),
    ):
        stats = stats.update(name, value)

    if sclr_dim > 0:
        for sclr in range(sclr_dim):
            sclr_idx = sclr + 1
            stats = stats.update(f"sclr{sclr_idx}m", sclrm[:, :, sclr])
            stats = stats.update(f"sclr{sclr_idx}m_f", sclrm_forcing[:, :, sclr])

    if edsclr_dim > 0:
        for edsclr in range(edsclr_dim):
            edsclr_idx = edsclr + 1
            stats = stats.update(f"edsclr{edsclr_idx}m", edsclrm[:, :, edsclr])
            stats = stats.update(
                f"edsclr{edsclr_idx}m_f", edsclrm_forcing[:, :, edsclr],
            )

    for name, value in (
        # stats_zm variables
        ("wm_zm", wm_zm),
        ("ddzt_umvm_sqd", ddzt_umvm_sqd),
        ("wp2", wp2),
        ("rtp2", rtp2),
        ("thlp2", thlp2),
        ("rtpthlp", rtpthlp),
        ("wprtp", wprtp),
        ("wpthlp", wpthlp),
    ):
        stats = stats.update(name, value)
    if l_stability_correct_tau_zm:
        stats = stats.update("stability_correction", stability_correction)
    for name, value in (
        ("Kh_zm", Kh_zm),
        ("upwp", upwp),
        ("vpwp", vpwp),
        ("vp2", vp2),
        ("up2", up2),
        ("rho_zm", rho_zm),
        ("rho_ds_zm", rho_ds_zm),
        ("thv_ds_zm", thv_ds_zm),
        ("em", em),
        ("wp3_on_wp2", wp3_on_wp2),
        ("wp3_on_wp2_cfl_num", wp3_on_wp2 * dt / dzm),
    ):
        stats = stats.update(name, value)

    if stats.var_on_stats_list("gamma_Skw_fnc"):
        wp3_zm = zt2zm(nzm, nzt, ngrdcol, gr, wp3)
        Skw_zm = Skx_func(nzm, ngrdcol, wp2, wp3_zm, w_tol, clubb_params)
        gamma_Skw_fnc = compute_gamma_Skw(
            nzm, ngrdcol, Skw_zm, clubb_params, l_gamma_Skw,
        )
        stats = stats.update("gamma_Skw_fnc", gamma_Skw_fnc)

    if sclr_dim > 0:
        for sclr in range(sclr_dim):
            sclr_idx = sclr + 1
            stats = stats.update(f"sclr{sclr_idx}p2", sclrp2[:, :, sclr])
            stats = stats.update(f"sclr{sclr_idx}prtp", sclrprtp[:, :, sclr])
            stats = stats.update(f"sclr{sclr_idx}pthlp", sclrpthlp[:, :, sclr])
            stats = stats.update(f"wpsclr{sclr_idx}p", wpsclrp[:, :, sclr])

    if edsclr_dim > 0:
        for edsclr in range(edsclr_dim):
            edsclr_idx = edsclr + 1
            stats = stats.update(f"wpedsclr{edsclr_idx}p", wpedsclrp[:, :, edsclr])

    # stats_sfc variables
    for i in range(ngrdcol):
        stats = stats.update("cc", jnp.max(cloud_frac[i, :]), icol=i)

    if stats.var_on_stats_list("z_cloud_base"):
        for i in range(ngrdcol):
            cloudy = rcm[i, :] >= rc_tol
            k = jnp.argmax(cloudy)
            k_prev = jnp.maximum(k - 1, 0)
            rcm_k = rcm[i, k]
            rcm_km1 = rcm[i, k_prev]
            zt_k = gr.zt[i, k]
            zt_km1 = gr.zt[i, k_prev]
            denom = jnp.where(jnp.abs(rcm_k - rcm_km1) > eps, rcm_k - rcm_km1, 1.0)
            z_cloud_base = zt_km1 + (rc_tol - rcm_km1) * (zt_k - zt_km1) / denom
            z_cloud_base = jnp.where(k == 0, gr.zt[i, 0], z_cloud_base)
            z_cloud_base = jnp.where((jnp.any(cloudy)) & (k < nzt - 1), z_cloud_base, -10.0)
            z_cloud_base = jnp.where((jnp.any(cloudy)) & (k == 0), gr.zt[i, 0], z_cloud_base)
            stats = stats.update("z_cloud_base", z_cloud_base, icol=i)

    if stats.var_on_stats_list("lwp"):
        for i in range(ngrdcol):
            stats = stats.update(
                "lwp",
                vertical_integral(nzt, rho_ds_zt[i, :], rcm[i, :], dzt[i, :]),
                icol=i,
            )

    if stats.var_on_stats_list("vwp"):
        for i in range(ngrdcol):
            stats = stats.update(
                "vwp",
                vertical_integral(nzt, rho_ds_zt[i, :], rtm[i, :] - rcm[i, :], dzt[i, :]),
                icol=i,
            )

    for i in range(ngrdcol):
        stats = stats.update(
            "thlm_vert_avg",
            vertical_avg(nzt, rho_ds_zt[i, :], thlm[i, :], dzt[i, :]),
            icol=i,
        )
        stats = stats.update(
            "rtm_vert_avg",
            vertical_avg(nzt, rho_ds_zt[i, :], rtm[i, :], dzt[i, :]),
            icol=i,
        )
        stats = stats.update(
            "um_vert_avg",
            vertical_avg(nzt, rho_ds_zt[i, :], um[i, :], dzt[i, :]),
            icol=i,
        )
        stats = stats.update(
            "vm_vert_avg",
            vertical_avg(nzt, rho_ds_zt[i, :], vm[i, :], dzt[i, :]),
            icol=i,
        )
        stats = stats.update(
            "wp2_vert_avg",
            vertical_avg(nzm, rho_ds_zm[i, :], wp2[i, :], dzm[i, :]),
            icol=i,
        )
        stats = stats.update(
            "up2_vert_avg",
            vertical_avg(nzm, rho_ds_zm[i, :], up2[i, :], dzm[i, :]),
            icol=i,
        )
        stats = stats.update(
            "vp2_vert_avg",
            vertical_avg(nzm, rho_ds_zm[i, :], vp2[i, :], dzm[i, :]),
            icol=i,
        )
        stats = stats.update(
            "rtp2_vert_avg",
            vertical_avg(nzm, rho_ds_zm[i, :], rtp2[i, :], dzm[i, :]),
            icol=i,
        )
        stats = stats.update(
            "thlp2_vert_avg",
            vertical_avg(nzm, rho_ds_zm[i, :], thlp2[i, :], dzm[i, :]),
            icol=i,
        )

    if stats.var_on_stats_list("tot_vartn_normlzd_rtm"):
        for i in range(ngrdcol):
            span = jnp.abs(rtm[i, nzt - 1] - rtm[i, 0])
            val = jnp.where(
                span < eps,
                -999.0,
                jnp.sum(jnp.abs(rtm[i, 1:nzt] - rtm[i, :nzt - 1])) / span,
            )
            stats = stats.update("tot_vartn_normlzd_rtm", val, icol=i)

    if stats.var_on_stats_list("tot_vartn_normlzd_thlm"):
        for i in range(ngrdcol):
            span = jnp.abs(thlm[i, nzt - 1] - thlm[i, 0])
            val = jnp.where(
                span < eps,
                -999.0,
                jnp.sum(jnp.abs(thlm[i, 1:nzt] - thlm[i, :nzt - 1])) / span,
            )
            stats = stats.update("tot_vartn_normlzd_thlm", val, icol=i)

    if stats.var_on_stats_list("tot_vartn_normlzd_wprtp"):
        for i in range(ngrdcol):
            span = jnp.abs(wprtp[i, nzm - 1] - wprtp[i, 0])
            val = jnp.where(
                span < eps,
                -999.0,
                jnp.sum(jnp.abs(wprtp[i, 1:nzm] - wprtp[i, :nzm - 1])) / span,
            )
            stats = stats.update("tot_vartn_normlzd_wprtp", val, icol=i)

    for i in range(ngrdcol):
        rtm_flux_top = rho_ds_zm[i, gr.k_ub_zm] * wprtp[i, gr.k_ub_zm]
        rtm_flux_sfc = jnp.where(
            l_host_applies_sfc_fluxes,
            0.0,
            rho_ds_zm[i, gr.k_lb_zm] * wprtp_sfc[i],
        )
        # Get the vertical integral of rtm before this function begins
        # so that spurious source can be calculated.
        rtm_integral_before = vertical_integral(
            nzt, rho_ds_zt[i, :], rtm_before[i, :], dzt[i, :],
        )
        rtm_integral_after = vertical_integral(
            nzt, rho_ds_zt[i, :], rtm[i, :], dzt[i, :],
        )
        rtm_integral_forcing = vertical_integral(
            nzt, rho_ds_zt[i, :], rtm_forcing[i, :], dzt[i, :],
        )
        # Calculate the spurious source for rtm.
        rtm_spur_src = _calculate_spurious_source(
            rtm_integral_after,
            rtm_integral_before,
            rtm_flux_top,
            rtm_flux_sfc,
            rtm_integral_forcing,
            dt,
        )

        thlm_flux_top = rho_ds_zm[i, gr.k_ub_zm] * wpthlp[i, gr.k_ub_zm]
        thlm_flux_sfc = jnp.where(
            l_host_applies_sfc_fluxes,
            0.0,
            rho_ds_zm[i, gr.k_lb_zm] * wpthlp_sfc[i],
        )
        # Get the vertical integral of thlm before this function begins
        # so that spurious source can be calculated.
        thlm_integral_before = vertical_integral(
            nzt, rho_ds_zt[i, :], thlm_before[i, :], dzt[i, :],
        )
        thlm_integral_after = vertical_integral(
            nzt, rho_ds_zt[i, :], thlm[i, :], dzt[i, :],
        )
        thlm_integral_forcing = vertical_integral(
            nzt, rho_ds_zt[i, :], thlm_forcing[i, :], dzt[i, :],
        )
        # Calculate the spurious source for thlm.
        thlm_spur_src = _calculate_spurious_source(
            thlm_integral_after,
            thlm_integral_before,
            thlm_flux_top,
            thlm_flux_sfc,
            thlm_integral_forcing,
            dt,
        )

        l_compute_spurious_source = jnp.asarray(l_implemented) | (
            jnp.all(jnp.abs(wm_zt[i, :]) < eps) & jnp.all(jnp.abs(wm_zm[i, :]) < eps)
        )
        # Spurious source will only be calculated if rtm_ma and thlm_ma are zero.
        # Therefore, wm must be zero or l_implemented must be true.
        stats = stats.update(
            "rtm_spur_src",
            jnp.where(l_compute_spurious_source, rtm_spur_src, -9999.0),
            icol=i,
        )
        stats = stats.update(
            "thlm_spur_src",
            jnp.where(l_compute_spurious_source, thlm_spur_src, -9999.0),
            icol=i,
        )

    return stats


__all__ = ["stats_accumulate"]
