"""Apply case-specific large-scale and surface forcings.

TODO(port-scope): This standalone JAX driver intentionally omits the Fortran
CLEX cases and dycore-grid forcing/remapping branches. CLEX input readers have
not been ported, and this driver owns only the CLUBB physics grid rather than a
host model's dycore grid. Consequently the supported-path interface below is
narrower than the Fortran routine instead of retaining unusable arguments and
branches. Add those source paths back only with their data owners and end-to-end
correctness tests.

References:
    None.
"""

from __future__ import annotations

from functools import partial

import jax
import jax.numpy as jnp

from clubb_jax.src.Benchmark_cases.arm import arm_sfclyr
from clubb_jax.src.Benchmark_cases.arm_0003 import arm_0003_sfclyr
from clubb_jax.src.Benchmark_cases.arm_3year import arm_3year_sfclyr
from clubb_jax.src.Benchmark_cases.arm_97 import arm_97_sfclyr
from clubb_jax.src.Benchmark_cases import time_dependent_input
from clubb_jax.src.Benchmark_cases.astex_a209 import astex_a209_sfclyr
from clubb_jax.src.Benchmark_cases.atex import atex_sfclyr, atex_tndcy
from clubb_jax.src.Benchmark_cases.atex_long import atex_long_sfclyr, atex_long_tndcy
from clubb_jax.src.Benchmark_cases.bomex import bomex_sfclyr, bomex_tndcy
from clubb_jax.src.Benchmark_cases.cloud_feedback import cloud_feedback_sfclyr
from clubb_jax.src.Benchmark_cases.cobra import cobra_sfclyr
from clubb_jax.src.Benchmark_cases.dycoms2_rf01 import (
    dycoms2_rf01_sfclyr,
    dycoms2_rf01_tndcy,
)
from clubb_jax.src.Benchmark_cases.dycoms2_rf02 import (
    dycoms2_rf02_sfclyr,
    dycoms2_rf02_tndcy,
)
from clubb_jax.src.Benchmark_cases.ekman import ekman_sfclyr
from clubb_jax.src.Benchmark_cases.fire import fire_sfclyr
from clubb_jax.src.Benchmark_cases.gabls2 import gabls2_sfclyr, gabls2_tndcy
from clubb_jax.src.Benchmark_cases.gabls3 import gabls3_sfclyr
from clubb_jax.src.Benchmark_cases.gabls3_night import gabls3_night_sfclyr
from clubb_jax.src.Benchmark_cases.jun25 import jun25_altocu_read_t_dependent
from clubb_jax.src.Benchmark_cases.lba import lba_sfclyr, lba_tndcy
from clubb_jax.src.Benchmark_cases.mpace_a import mpace_a_sfclyr, mpace_a_tndcy
from clubb_jax.src.Benchmark_cases.mpace_b import mpace_b_sfclyr, mpace_b_tndcy
from clubb_jax.src.Benchmark_cases.neutral_case import neutral_case_sfclyr
from clubb_jax.src.Benchmark_cases.nov11 import (
    nov11_altocu_read_t_dependent,
    nov11_altocu_rtm_adjust,
)
from clubb_jax.src.Benchmark_cases.rico import rico_sfclyr, rico_tndcy
from clubb_jax.src.Benchmark_cases.sfc_flux import (
    compute_momentum_flux,
    compute_ubar,
    set_sclr_sfc_rtm_thlm,
)
from clubb_jax.src.Benchmark_cases.twp_ice import twp_ice_sfclyr
from clubb_jax.src.Benchmark_cases.wangara import wangara_sfclyr, wangara_tndcy
from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv, kappa, p0
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.grid_class import zt2zm
from clubb_jax.src.CLUBB_core.interpolation import mono_cubic_interp


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "sclr_dim",
        "edsclr_dim",
        "runtype",
        "sfctype",
        "l_modify_bc_for_cnvg_test",
        "saturation_formula",
    ),
)
def prescribe_forcings(
    gr, nzm, nzt, ngrdcol,
    sclr_dim, edsclr_dim, sclr_idx,
    runtype, sfctype,
    time_current, time_initial, dt,
    um, vm, thlm,
    p_in_Pa, exner, rho, rho_zm, thvm,
    veg_T_in_K,
    l_modify_bc_for_cnvg_test,
    saturation_formula,
    stats,
    rtm, wm_zm, wm_zt, ug, vg, um_ref, vm_ref,
    thlm_forcing, rtm_forcing, um_forcing,
    vm_forcing, wprtp_forcing, wpthlp_forcing,
    rtp2_forcing, thlp2_forcing, rtpthlp_forcing,
    wpsclrp, sclrm_forcing, edsclrm_forcing,
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
    T_sfc, p_sfc, sens_ht, latent_ht,
    wpsclrp_sfc, wpedsclrp_sfc, err_info,
):
    """Calculate tendency and surface variables.

    References:
        None.

    The dycore-only arguments are intentionally omitted as documented in the
    module TODO. In/out values are returned in their source argument order.
    """
    l_t_dependent = time_dependent_input.l_t_dependent
    l_ignore_forcings = time_dependent_input.l_ignore_forcings

    # -----------------------------------------------------------------------
    #                    FIND ALL DIAGNOSTIC VARIABLES
    # -----------------------------------------------------------------------
    l_compute_momentum_flux = False
    l_set_sclr_sfc_rtm_thlm = False
    l_fixed_flux = False

    # ----------------------------------------------------------------
    # Set vertical velocity, w, and compute large-scale forcings.
    # ----------------------------------------------------------------
    # These lines were added to reset the forcing arrays to 0 each iteration.
    # This was previously done in the <case>_tndcy subroutine.
    rtm_forcing = jnp.zeros_like(rtm_forcing)
    thlm_forcing = jnp.zeros_like(thlm_forcing)
    wprtp_forcing = jnp.zeros_like(wprtp_forcing)
    wpthlp_forcing = jnp.zeros_like(wpthlp_forcing)
    rtp2_forcing = jnp.zeros_like(rtp2_forcing)
    thlp2_forcing = jnp.zeros_like(thlp2_forcing)
    rtpthlp_forcing = jnp.zeros_like(rtpthlp_forcing)

    if l_t_dependent and not l_ignore_forcings:
        (
            thlm_forcing, rtm_forcing, um_ref, vm_ref, um_forcing, vm_forcing,
            wm_zt, wm_zm, ug, vg,
            sclrm_forcing, edsclrm_forcing,
        ) = time_dependent_input.apply_time_dependent_forcings(
            ngrdcol, gr.nzm, gr.nzt,
            sclr_dim, edsclr_dim, sclr_idx,
            gr, time_current, rtm, rho, exner,
            thlm_forcing, rtm_forcing, um_ref, vm_ref, um_forcing, vm_forcing,
            wm_zt, wm_zm, ug, vg,
            sclrm_forcing, edsclrm_forcing,
        )

        # Vince Larson set forcing to zero at the top point so that we don't need
        # so much sponge damping, which is associated with sawtooth noise
        # in the cloud_feedback cases. I don't know how it will affect
        # the other cases.
        rtm_forcing = rtm_forcing.at[:, nzt - 1].set(0.0)
        thlm_forcing = thlm_forcing.at[:, nzt - 1].set(0.0)
    else:
        # Legacy method of setting the forcings.
        if runtype == "atex":
            (
                err_info,
                wm_zt, wm_zm,
                thlm_forcing, rtm_forcing,
                sclrm_forcing, edsclrm_forcing,
            ) = atex_tndcy(
                ngrdcol, sclr_dim, edsclr_dim, sclr_idx,
                gr, time_current, time_initial,
                rtm,
                err_info,
            )
            # TODO(port-mirror): A traced fatal status cannot drive the source's
            # early return. Remove this guard when compiled diagnostics are carried.
            if (
                clubb_at_least_debug_level(0)
                and not isinstance(err_info.err_code, jax.core.Tracer)
                and err_info.is_fatal_host()
            ):
                return (
                    stats,
                    rtm, wm_zm, wm_zt, ug, vg, um_ref, vm_ref,
                    thlm_forcing, rtm_forcing, um_forcing,
                    vm_forcing, wprtp_forcing, wpthlp_forcing,
                    rtp2_forcing, thlp2_forcing, rtpthlp_forcing,
                    wpsclrp, sclrm_forcing, edsclrm_forcing,
                    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
                    T_sfc, p_sfc, sens_ht, latent_ht,
                    wpsclrp_sfc, wpedsclrp_sfc, err_info,
                )
        elif runtype == "atex_long":
            (
                wm_zt, wm_zm,
                thlm_forcing, rtm_forcing,
                sclrm_forcing, edsclrm_forcing,
            ) = atex_long_tndcy(
                ngrdcol, sclr_dim, edsclr_dim, sclr_idx,
                gr, time_current,
            )
        elif runtype == "bomex":
            (
                thlm_forcing, rtm_forcing,
                sclrm_forcing, edsclrm_forcing,
            ) = bomex_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr, rtm)
        elif runtype == "dycoms2_rf01":
            (
                thlm_forcing, rtm_forcing,
                sclrm_forcing, edsclrm_forcing,
            ) = dycoms2_rf01_tndcy(ngrdcol, gr, sclr_dim, edsclr_dim, sclr_idx)
        elif runtype == "dycoms2_rf02":
            (
                wm_zt, wm_zm,
                thlm_forcing, rtm_forcing,
                sclrm_forcing, edsclrm_forcing,
            ) = dycoms2_rf02_tndcy(
                ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr,
                wm_zt, wm_zm,
            )
        elif runtype in ("fire", "generic", "coriolis_test", "ekman"):
            # Analytic radiation is computed elsewhere.
            pass
        elif runtype == "gabls2":
            (
                wm_zt, wm_zm, thlm_forcing,
                rtm_forcing,
                sclrm_forcing, edsclrm_forcing,
            ) = gabls2_tndcy(
                ngrdcol, sclr_dim, edsclr_dim, sclr_idx,
                gr, time_current, time_initial,
            )
        elif runtype == "lba":
            (
                thlm_forcing, rtm_forcing,
                sclrm_forcing, edsclrm_forcing,
            ) = lba_tndcy(ngrdcol, sclr_dim, edsclr_dim, sclr_idx, gr)
        elif runtype == "mpace_a":
            (
                wm_zt, wm_zm, thlm_forcing, rtm_forcing,
                um_ref, vm_ref,
                sclrm_forcing, edsclrm_forcing,
            ) = mpace_a_tndcy(
                ngrdcol, sclr_dim, edsclr_dim, sclr_idx,
                gr, time_current, p_in_Pa,
            )
        elif runtype == "mpace_b":
            (
                wm_zt, wm_zm, thlm_forcing, rtm_forcing,
                sclrm_forcing, edsclrm_forcing,
            ) = mpace_b_tndcy(
                ngrdcol, sclr_dim, edsclr_dim, sclr_idx,
                gr, p_in_Pa, thvm,
            )
        elif runtype == "rico":
            (
                thlm_forcing, rtm_forcing,
                sclrm_forcing, edsclrm_forcing,
            ) = rico_tndcy(
                ngrdcol, sclr_dim, edsclr_dim, sclr_idx,
                gr, rtm, exner,
            )
        elif runtype == "neutral":
            sclrm_forcing = jnp.zeros_like(sclrm_forcing)
            edsclrm_forcing = jnp.zeros_like(edsclrm_forcing)
        elif runtype == "wangara":
            (
                wm_zt, wm_zm,
                thlm_forcing, rtm_forcing,
                sclrm_forcing, edsclrm_forcing,
            ) = wangara_tndcy(ngrdcol, gr, sclr_dim, edsclr_dim, sclr_idx)
        else:
            raise NotImplementedError(
                f"prescribe_forcings: unknown legacy forcing for runtype {runtype!r}"
            )

    # ----------------------------------------------------------------
    # Compute surface fluxes.
    # ----------------------------------------------------------------
    # Derive physical quantities at the bottom model level for the surface
    # boundary conditions.
    (
        z_bot, um_bot, vm_bot, rtm_bot,
        thlm_bot, rho_bot, exner_bot,
    ) = read_surface_var_for_bc(
        gr, ngrdcol,
        um, vm, rtm, thlm, rho_zm, exner,
        p_sfc, l_modify_bc_for_cnvg_test,
    )

    # Boundary conditions for the second-order moments.
    ubar = compute_ubar(ngrdcol, um_bot, vm_bot)

    upwp_sfc = jnp.asarray(upwp_sfc)
    vpwp_sfc = jnp.asarray(vpwp_sfc)
    wpthlp_sfc = jnp.asarray(wpthlp_sfc)
    wprtp_sfc = jnp.asarray(wprtp_sfc)
    T_sfc = jnp.asarray(T_sfc)
    ustar = jnp.zeros(ngrdcol)

    if runtype == "rico":
        l_set_sclr_sfc_rtm_thlm = True
        (
            upwp_sfc, vpwp_sfc, wpthlp_sfc,
            wprtp_sfc, ustar, T_sfc,
        ) = rico_sfclyr(
            ngrdcol, time_current, um_bot, vm_bot, thlm_bot, rtm_bot,
            z_bot, p_sfc, exner_bot,
            saturation_formula,
        )
    elif runtype == "gabls3":
        l_compute_momentum_flux = True
        wpthlp_sfc, wprtp_sfc, ustar = gabls3_sfclyr(
            ngrdcol, ubar, veg_T_in_K,
            thlm_bot, rtm_bot, z_bot, exner_bot,
            wpthlp_sfc, wprtp_sfc,
        )
    elif runtype == "gabls3_night":
        (
            upwp_sfc, vpwp_sfc,
            wpthlp_sfc, wprtp_sfc, ustar,
        ) = gabls3_night_sfclyr(
            ngrdcol, time_current, um_bot, vm_bot,
            thlm_bot, rtm_bot, z_bot,
        )
    elif runtype == "jun25_altocu":
        # There are no surface momentum or heat fluxes
        # for the Jun. 25 Altocumulus case.

        # Ensure ustar(i) is set
        ustar = jnp.zeros(ngrdcol)

        # Read in time dependent inputs
        sens_ht, latent_ht = jun25_altocu_read_t_dependent(time_current)
    elif runtype == "cobra":
        l_compute_momentum_flux = True
        (
            wpthlp_sfc, wprtp_sfc, ustar,
            wpsclrp_sfc, wpedsclrp_sfc, T_sfc,
        ) = cobra_sfclyr(
            ngrdcol, sclr_dim, edsclr_dim, sclr_idx,
            time_current, z_bot, rho_bot,
            thlm_bot, ubar,
        )
    elif runtype == "astex_a209":
        l_compute_momentum_flux = True
        wpthlp_sfc, wprtp_sfc, ustar, T_sfc = astex_a209_sfclyr(
            ngrdcol, time_current, ubar, rtm_bot,
            thlm_bot, z_bot, exner_bot, p_sfc,
            saturation_formula,
        )
    elif runtype == "nov11_altocu":
        # There are no surface momentum or heat fluxes. However, this case has
        # a one-time rtm adjustment one hour after initialization.
        rtm = nov11_altocu_rtm_adjust(
            ngrdcol, gr,
            time_current, time_initial, dt,
            rtm,
        )
        sens_ht, latent_ht = nov11_altocu_read_t_dependent(time_current)
    elif runtype in ("fire", "generic"):
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        l_fixed_flux = True
        wpthlp_sfc, wprtp_sfc, ustar, T_sfc = fire_sfclyr(
            ngrdcol, time_current, ubar, p_sfc,
            thlm_bot, rtm_bot, exner_bot,
            saturation_formula,
        )
    elif runtype in (
        "cloud_feedback_s6",
        "cloud_feedback_s6_p2k",
        "cloud_feedback_s11",
        "cloud_feedback_s11_p2k",
        "cloud_feedback_s12",
        "cloud_feedback_s12_p2k",
        "cgils_s6",
        "cgils_s6_p2k",
        "cgils_s11",
        "cgils_s11_p2k",
        "cgils_s12",
        "cgils_s12_p2k",
    ):
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        l_fixed_flux = True
        wpthlp_sfc, wprtp_sfc, ustar, T_sfc = cloud_feedback_sfclyr(
            ngrdcol, time_current, sfctype,
            thlm_bot, rtm_bot, z_bot,
            ubar, p_sfc,
            saturation_formula,
        )
    elif runtype == "arm":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar = arm_sfclyr(
            ngrdcol, time_current, z_bot,
            thlm_bot, ubar,
        )
    elif runtype == "arm_0003":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar = arm_0003_sfclyr(
            ngrdcol, time_current, z_bot,
            rho_bot, thlm_bot, ubar,
        )
    elif runtype == "arm_3year":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar = arm_3year_sfclyr(
            ngrdcol, time_current, z_bot, rho_bot,
            thlm_bot, ubar,
        )
    elif runtype in ("arm_97", "mc3e"):
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar = arm_97_sfclyr(
            ngrdcol, time_current, z_bot, rho_bot,
            thlm_bot, ubar,
        )
    elif runtype == "atex":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar, T_sfc = atex_sfclyr(
            ngrdcol, time_current, ubar,
            thlm_bot, rtm_bot, exner_bot,
        )
    elif runtype == "atex_long":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar, T_sfc = atex_long_sfclyr(
            ngrdcol, time_current, ubar,
            thlm_bot, rtm_bot, exner_bot, rho_bot,
        )
    elif runtype == "bomex":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar = bomex_sfclyr(ngrdcol, time_current, rtm_bot)
    elif runtype == "dycoms2_rf01":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar, T_sfc = dycoms2_rf01_sfclyr(
            ngrdcol, time_current, sfctype, p_sfc,
            exner_bot, ubar,
            thlm_bot, rtm_bot, rho_bot,
            saturation_formula,
        )
    elif runtype == "dycoms2_rf02":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar = dycoms2_rf02_sfclyr(ngrdcol, time_current)
    elif runtype == "gabls2":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar, T_sfc = gabls2_sfclyr(
            ngrdcol, time_current, time_initial,
            z_bot, p_sfc,
            ubar, thlm_bot, rtm_bot, exner_bot,
            saturation_formula,
        )
    elif runtype == "lba":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar = lba_sfclyr(
            ngrdcol, time_current, time_initial,
            z_bot, rho_bot, thlm_bot, ubar,
        )
    elif runtype == "mpace_a":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar = mpace_a_sfclyr(ngrdcol, time_current, rho_bot)
    elif runtype == "mpace_b":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar = mpace_b_sfclyr(ngrdcol, time_current, rho_bot)
    elif runtype == "neutral":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        (
            upwp_sfc, vpwp_sfc,
            wpthlp_sfc, wprtp_sfc, ustar,
        ) = neutral_case_sfclyr(
            ngrdcol, time_current,
            um_bot, vm_bot, ubar,
        )
    elif runtype == "ekman":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        (
            upwp_sfc, vpwp_sfc,
            wpthlp_sfc, wprtp_sfc, ustar,
        ) = ekman_sfclyr(
            ngrdcol,
            um_bot, vm_bot, ubar,
        )
    elif runtype == "twp_ice":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar, T_sfc = twp_ice_sfclyr(
            ngrdcol, time_current, z_bot, exner_bot,
            thlm_bot, ubar, rtm_bot, p_sfc,
            saturation_formula,
        )
    elif runtype == "wangara":
        l_compute_momentum_flux = True
        l_set_sclr_sfc_rtm_thlm = True
        wpthlp_sfc, wprtp_sfc, ustar = wangara_sfclyr(ngrdcol, time_current)
    elif runtype == "coriolis_test":
        l_compute_momentum_flux = True
        l_fixed_flux = True
        ustar = jnp.zeros(ngrdcol)
    else:
        raise NotImplementedError(f"invalid runtype {runtype!r}")

    # These have been placed here to help avoid repetition in the cases
    if l_compute_momentum_flux:
        upwp_sfc, vpwp_sfc = compute_momentum_flux(
            ngrdcol, um_bot, vm_bot, ubar, ustar,
        )

    if l_set_sclr_sfc_rtm_thlm:
        wpsclrp_sfc, wpedsclrp_sfc = set_sclr_sfc_rtm_thlm(
            ngrdcol, sclr_dim, edsclr_dim, sclr_idx,
            wpthlp_sfc, wprtp_sfc,
        )

    # If the surface type is 0, use fixed fluxes
    if sfctype == 0 and l_fixed_flux:
        wpthlp_sfc = jnp.full(ngrdcol, sens_ht)
        wprtp_sfc = jnp.full(ngrdcol, latent_ht)
        if sclr_idx.iisclr_thl > 0:
            wpsclrp = wpsclrp.at[:, :, sclr_idx.iisclr_thl - 1].set(sens_ht)
        if sclr_idx.iisclr_rt > 0:
            wpsclrp = wpsclrp.at[:, :, sclr_idx.iisclr_rt - 1].set(latent_ht)

    wpthlp_sfc = jnp.asarray(wpthlp_sfc)
    wprtp_sfc = jnp.asarray(wprtp_sfc)
    upwp_sfc = jnp.asarray(upwp_sfc)
    vpwp_sfc = jnp.asarray(vpwp_sfc)
    ustar = jnp.asarray(ustar)
    T_sfc = jnp.asarray(T_sfc)

    # Store values of surface fluxes for statistics
    if stats.l_sample:
        rho_zm_sfc = rho_zm[:, 0]
        stats = stats.update("sh", wpthlp_sfc * rho_zm_sfc * Cp)
        stats = stats.update("lh", wprtp_sfc * rho_zm_sfc * Lv)
        stats = stats.update("wpthlp_sfc", wpthlp_sfc)
        stats = stats.update("wprtp_sfc", wprtp_sfc)
        stats = stats.update("upwp_sfc", upwp_sfc)
        stats = stats.update("vpwp_sfc", vpwp_sfc)
        stats = stats.update("ustar", ustar)
        stats = stats.update("T_sfc", T_sfc)

    return (
        stats,
        rtm, wm_zm, wm_zt, ug, vg, um_ref, vm_ref,
        thlm_forcing, rtm_forcing, um_forcing,
        vm_forcing, wprtp_forcing, wpthlp_forcing,
        rtp2_forcing, thlp2_forcing, rtpthlp_forcing,
        wpsclrp, sclrm_forcing, edsclrm_forcing,
        wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
        T_sfc, p_sfc, sens_ht, latent_ht,
        wpsclrp_sfc, wpedsclrp_sfc, err_info,
    )


# -------------------------------------------------------------------------------
def read_surface_var_for_bc(
    gr,
    ngrdcol,
    um,
    vm,
    rtm,
    thlm,
    rho_zm,
    exner,
    p_sfc,
    l_modify_bc_for_cnvg_test,
):
    """Derives the physical quantities at the bottom model level for calculating
    surface fluxes (boundary conditions). The default option is to use the
    quantities at first/second model level. When l_modify_bc_for_cnvg_test is
    true, the quantities at a fixed model height (25m) are obtained via vertical
    interpolation and used for calculating the surface fluxes. The purpose is
    to eleminate the space-dependence of quantities in the default option when
    the model is refined vertically, which results in a space-dependence of
    surface fluxes. The modified option is found to be the correct treatment for
    evaluating space-time convergence in CLUBB-SCM.

    Author: Shixuan Zhang (Shixuan.Zhang@pnnl.gov).
    """
    constant_height_option = 2

    if not l_modify_bc_for_cnvg_test:
        # Default model setup in CLUBB-SCM
        z_bot = gr.zt[:, 0]
        um_bot = um[:, 0]
        vm_bot = vm[:, 0]
        rtm_bot = rtm[:, 0]
        thlm_bot = thlm[:, 0]
        rho_bot = rho_zm[:, 0]
        exner_bot = (p_sfc / p0) ** kappa
        return z_bot, um_bot, vm_bot, rtm_bot, thlm_bot, rho_bot, exner_bot

    # Modified option which find the values of physical quantities
    # at a fixed model height (25m)
    z_bot = jnp.full(ngrdcol, 25.0)
    k_min = jnp.argmin(jnp.abs(gr.zt - 25.0), axis=1)

    if constant_height_option == 1:
        column = jnp.arange(ngrdcol)
        return (
            z_bot,
            um[column, k_min],
            vm[column, k_min],
            rtm[column, k_min],
            thlm[column, k_min],
            rho_zm[column, k_min],
            exner[column, k_min],
        )

    # option 2 (interpolation)
    um_zm = zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, um)
    vm_zm = zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, vm)
    thlm_zm = zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, thlm)
    rtm_zm = zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, rtm)
    exner_zm = zt2zm(gr.nzm, gr.nzt, ngrdcol, gr, exner)
    exner_zm = exner_zm.at[:, 0].set((p_sfc / p0) ** kappa)

    def interpolate_column(
        z_target,
        k00,
        zm,
        rho,
        um_values,
        vm_values,
        rtm_values,
        thlm_values,
        exner_values,
    ):
        """Apply the source column body to one column; vmap supplies columns."""
        km1 = jnp.maximum(k00 - 1, 0)
        kp1 = k00 + 1
        kp2 = k00 + 2

        def interpolate(values):
            return mono_cubic_interp(
                z_target,
                km1,
                k00,
                kp1,
                kp2,
                zm[km1],
                zm[k00],
                zm[kp1],
                zm[kp2],
                values[km1],
                values[k00],
                values[kp1],
                values[kp2],
            )

        return (
            interpolate(um_values),
            interpolate(vm_values),
            interpolate(rtm_values),
            interpolate(thlm_values),
            rho[k00],
            interpolate(exner_values),
        )

    um_bot, vm_bot, rtm_bot, thlm_bot, rho_bot, exner_bot = jax.vmap(
        interpolate_column
    )(z_bot, k_min, gr.zm, rho_zm, um_zm, vm_zm, rtm_zm, thlm_zm, exner_zm)

    return (
        z_bot,
        um_bot,
        vm_bot,
        rtm_bot,
        thlm_bot,
        rho_bot,
        exner_bot,
    )


__all__ = ["prescribe_forcings", "read_surface_var_for_bc"]
