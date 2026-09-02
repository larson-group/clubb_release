"""JAX port of KK_microphys_module.F90 — Khairoutdinov-Kogan (KK) warm-rain microphysics driver routines.

Mirrors clubb_release/src/Microphys/KK_microphys_module.F90: `KK_microphys_adjust` (assemble the second-moment
state tendencies from the process rates, with the source/evaporation clip adjustments) and `KK_sedimentation`
(KK00 Eq. 37 mean rain-mass/number sedimentation velocities from the mean volume radius). The per-process
upscaled mean rates (autoconversion/accretion/evaporation) and the full `compute_kk_microphysics` orchestration
live in `Microphys/KK_microphys/kk_microphys_driver.py`, which imports these two back.

Pure-jnp → differentiable. Validated against the rico oracle (`tests/test_kk_rico_oracle.py`).
"""

from __future__ import annotations

import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq, SATURATION_FLATAU
from clubb_jax.src.CLUBB_core.constants_clubb import Lv, Rv, Rd, Cp, rho_lw, cm3_per_m3, rr_tol, mvr_rain_max
from clubb_jax.src.Microphys.KK_microphys.KK_utilities import G_T_p
from clubb_jax.src.Microphys.KK_microphys.parameters_KK import KK_auto_Nc_exp as KK_AUTO_NC_EXP  # `use parameters_KK`


def kk_evap_coef(T_liq_in_K, p_in_Pa, C_evap, saturation_formula=SATURATION_FLATAU):
    """KK evaporation coefficient. KK_microphys_module.F90:1177.

    = 3 C_evap G_T_p ((4/3)pi rho_lw)^(2/3) (1 + Beta_Tl r_sl)/r_sl,
    Beta_Tl = (Rd/Rv)(Lv/(Rd T))(Lv/(Cp T)), r_sl = sat_mixrat_liq(p,T). Computed inline in the
    Fortran KK_microphys_module.F90 (not a standalone subroutine); factored here into its Fortran-home
    module. G_T_p (the drop-radius growth coefficient) is its real KK_utilities.F90 sibling."""
    T = jnp.asarray(T_liq_in_K, dtype=jnp.float64)
    p = jnp.asarray(p_in_Pa, dtype=jnp.float64)
    r_sl = sat_mixrat_liq(p, T, saturation_formula)
    Beta_Tl = (Rd / Rv) * (Lv / (Rd * T)) * (Lv / (Cp * T))
    return (3.0 * C_evap * G_T_p(T, p, saturation_formula)
            * ((4.0 / 3.0) * jnp.pi * rho_lw) ** (2.0 / 3.0)
            * ((1.0 + Beta_Tl * r_sl) / r_sl))


def kk_auto_coef(rho):
    """KK autoconversion coefficient = 1350 (rho/cm3_per_m3)^KK_auto_Nc_exp. KK_microphys_module.F90:1182.

    Computed inline in the Fortran KK_microphys_module (then passed as the KK_auto_coef argument into the
    KK_upscaled_means routines); factored here into its Fortran-home module. The KK_auto_Nc_exp exponent
    lives with the other KK exponents in parameters_KK.py."""
    return 1350.0 * (rho / cm3_per_m3) ** KK_AUTO_NC_EXP


# KK accretion coefficient (constant); KK_tendency_coefs sets KK_accr_coef = 67.0 (KK_microphys_module.F90:1185),
# then passes it to KK_accr_local_mean / the upscaled accretion mean.
KK_ACCR_COEF = 67.0

# KK mean-volume-radius coefficient; KK_tendency_coefs sets KK_mvr_coef = ((4/3)·π·ρ_lw)^(-1/3)
# (KK_microphys_module.F90:1188). (full-precision jnp.pi, as the KK path is Tier-C — see DESIGN.md.)
KK_MVR_COEF = ((4.0 / 3.0) * jnp.pi * rho_lw) ** (-1.0 / 3.0)


def KK_microphys_adjust(dt, exner, rcm, rrm, Nrm,
                        KK_evap_tndcy, KK_auto_tndcy, KK_accr_tndcy,
                        KK_Nrm_evap_tndcy, KK_Nrm_auto_tndcy,
                        l_src_adj_enabled=True, l_evap_adj_enabled=True):
    """Assemble the KK microphysics state tendencies from the process rates.
    KK_microphys_module.F90:1196 (the upscaled path enables both adjustments).

    Source adjustment: limit auto+accr so they don't draw more cloud water than available
    (rate <= rcm/dt). Evaporation adjustment: limit so rain can't go negative (>= -rrm/dt,
    -Nrm/dt). Returns (rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc)."""
    from clubb_jax.src.Microphys.KK_microphys.KK_Nrm_tendencies import (
        KK_Nrm_auto_mean, KK_Nrm_evap_local_mean)
    from clubb_jax.src.CLUBB_core.constants_clubb import Lv, Cp
    Nr_tol = rr_tol / ((4.0 / 3.0) * jnp.pi * rho_lw * mvr_rain_max ** 3)   # constants_clubb.F90:306
    eps = jnp.finfo(jnp.float64).eps

    rrm_source = KK_auto_tndcy + KK_accr_tndcy
    Nrm_source = KK_Nrm_auto_tndcy

    if l_src_adj_enabled:
        # Over a long step auto+accr may over-deplete rcm; cap the total source at rcm/dt.
        over = (rrm_source * dt) > rcm
        rrm_src_max = rcm / dt
        src_safe = jnp.where(rrm_source != 0.0, rrm_source, 1.0)
        rrm_auto_ratio = KK_auto_tndcy / src_safe
        rrm_src_adj = rrm_src_max - rrm_source
        Nrm_src_adj = KK_Nrm_auto_mean(rrm_auto_ratio * rrm_src_adj)
        rrm_source = jnp.where(over, rrm_src_max, rrm_source)
        Nrm_source = jnp.where(over, Nrm_source + Nrm_src_adj, Nrm_source)

    if l_evap_adj_enabled:
        rrm_evap_net = jnp.maximum(KK_evap_tndcy, -rrm / dt)
        # recompute Nrm evap from the net rrm evap when the rrm evap was limited
        limited = (jnp.abs(KK_evap_tndcy - rrm_evap_net)
                   > jnp.abs(KK_evap_tndcy + rrm_evap_net) * eps / 2.0) \
                  & (rrm > rr_tol) & (Nrm > Nr_tol)
        Nrm_evap_recomp = KK_Nrm_evap_local_mean(rrm_evap_net, Nrm, rrm, dt)
        Nrm_evap_net = jnp.where(limited, Nrm_evap_recomp, KK_Nrm_evap_tndcy)
        Nrm_evap_net = jnp.maximum(Nrm_evap_net, -Nrm / dt)
    else:
        rrm_evap_net = KK_evap_tndcy
        Nrm_evap_net = KK_Nrm_evap_tndcy

    rrm_mc = rrm_evap_net + rrm_source
    Nrm_mc = Nrm_evap_net + Nrm_source
    rvm_mc = -rrm_evap_net
    rcm_mc = -rrm_source
    thlm_mc = (Lv / (Cp * exner)) * rrm_mc
    return rrm_mc, Nrm_mc, rvm_mc, rcm_mc, thlm_mc


def KK_sedimentation(mvr, cloud_top_level=None, l_clip_positive_sed=True):
    """Mean KK sedimentation velocities Vrr (rain mass) and VNr (rain number) from the rain-drop
    mean volume radius. Oracle KK_microphys_module.F90:1542 (KK00 Eq. 37).

      Vrr = -(0.012 * mvr_micron - 0.2),  VNr = -(0.007 * mvr_micron - 0.1)   [m/s, downward<=0]

    `mvr` is the mean volume radius in METERS (ngrdcol, nzt) or (nzt,). With l_clip_positive_sed
    (true for rico / non-SILHS) the velocities are clipped to <= 0 (no upward sedimentation), and
    the top model level is set to 0 (zero flux through the model top). If `cloud_top_level` (the
    0-based thermodynamic index of cloud top, per column) is given, velocities ABOVE cloud top are
    zeroed (faithful to the Fortran cloud_top_level+1:nzt-1 slice); for rico this is largely a
    no-op since mvr~0 above cloud already clips to 0. Returns (Vrr, VNr), differentiable in mvr."""
    from clubb_jax.src.CLUBB_core.constants_clubb import micron_per_m
    mvr = jnp.asarray(mvr, dtype=jnp.float64)
    mvr_micron = micron_per_m * mvr
    Vrr = -(0.012 * mvr_micron - 0.2)
    VNr = -(0.007 * mvr_micron - 0.1)
    if l_clip_positive_sed:
        Vrr = jnp.minimum(Vrr, 0.0)
        VNr = jnp.minimum(VNr, 0.0)
        if cloud_top_level is not None:
            nzt = mvr.shape[-1]
            k = jnp.arange(nzt)
            ctl = jnp.asarray(cloud_top_level)[..., None]   # 0-based cloud-top index per column
            above = (k[None, :] if mvr.ndim > 1 else k) > ctl
            above = above & (ctl > 0)
            Vrr = jnp.where(above, 0.0, Vrr)
            VNr = jnp.where(above, 0.0, VNr)
    # Zero-flux boundary condition at the model top (highest level).
    Vrr = Vrr.at[..., -1].set(0.0)
    VNr = VNr.at[..., -1].set(0.0)
    return Vrr, VNr
