"""JAX port of ``src/CLUBB_core/sfc_varnce_module.F90``.

Porting deviations:
- The Fortran source sets ``l_andre_1978 = .false.`` as a local parameter, so
  this port mirrors the active "Previous code" branch and omits the dead Andre
  et al. 1978 branch.
- Stats are threaded explicitly with JaxStats because the Fortran routine uses
  global stats side effects that are not JAX-compatible state.
- Fortran inout arguments are returned as updated arrays, which keeps the
  routine functional and JIT-friendly.
- The Fortran debug path calls ``sfc_varnce_check`` and prints diagnostics to
  stderr.  The JAX path preserves the fatal error state with finite checks but
  omits host-side debug printing inside the jitted routine.

References:
Andre, J. C., G. De Moor, P. Lacarrere, G. Therry, and R. Du Vachat, 1978.
  Modeling the 24-Hour Evolution of the Mean and Turbulent Structures of
  the Planetary Boundary Layer.  J. Atmos. Sci., 35, 1861 -- 1883.
"""

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import (
    eps,
    grav,
    max_mag_correlation_flux,
    one,
    one_half,
    one_third,
    rt_tol,
    thl_tol,
    w_tol_sqd,
    wp2_max,
    zero,
)
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core import ErrInfo, Grid, SclrIdx


# Vince Larson increased ufmin to stabilize arm_97.  24 Jul 2007
ufmin = 0.01  # Minimum allowable value of u* [m/s]


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "sclr_dim",
        "l_vary_convect_depth",
    ),
)
def calc_sfc_varnce(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    sclr_dim: int,
    sclr_idx: SclrIdx,
    gr: Grid,
    dt,
    sfc_elevation,
    upwp_sfc,
    vpwp_sfc,
    wpthlp,
    wprtp_sfc,
    um,
    vm,
    Lscale_up,
    wpsclrp_sfc,
    lhs_splat_wp2,
    tau_zm,
    l_vary_convect_depth: bool,
    T0,
    up2_sfc_coef,
    a_const,
    stats: JaxStats,
    wp2,
    up2,
    vp2,
    thlp2,
    rtp2,
    rtpthlp,
    sclrp2,
    sclrprtp,
    sclrpthlp,
    err_info: ErrInfo,
):
    """This subroutine computes estimate of the surface thermodynamic and wind
    component second order moments."""
    del nzt, um, vm, Lscale_up

    # Constant Parameters

    # Logical for Andre et al., 1978 parameterization.
    l_andre_1978 = False

    # Defined height of 1 meter                [m]
    z_const = one
    # This value is made up! - Vince Larson 12 Jul 2005
    sclr_var_coef = 0.4

    #-------------------------- Begin Code --------------------------

    # Reflect surface varnce changes in budget
    if stats.l_sample:
        stats = stats.begin_budget("thlp2_sf", thlp2 / dt)
        stats = stats.begin_budget("rtp2_sf", rtp2 / dt)
        stats = stats.begin_budget("rtpthlp_sf", rtpthlp / dt)
        stats = stats.begin_budget("up2_sf", up2 / dt)
        stats = stats.begin_budget("vp2_sf", vp2 / dt)
        stats = stats.begin_budget("wp2_sf", wp2 / dt)

    level_idx = jnp.asarray(
        tuple(
            range(
                int(gr.k_lb_zm),
                int(gr.k_ub_zm) + int(gr.grid_dir_indx),
                int(gr.grid_dir_indx),
            )
        )
    )
    wpthlp_path = wpthlp[:, level_idx]
    height_path = gr.zm[:, level_idx] - sfc_elevation[:, None]
    not_positive_layer = jnp.logical_not(
        (wpthlp_path > zero) & (height_path < 1000.0)
    )
    has_stop_level = jnp.any(not_positive_layer, axis=1)
    stop_level_pos = jnp.argmax(not_positive_layer, axis=1)
    stop_level_pos = jnp.where(
        has_stop_level,
        stop_level_pos,
        level_idx.shape[0] - 1,
    )

    # Find thickness of layer near surface with positive heat flux.
    # This is used when l_vary_convect_depth=.true. in order to determine wp2.
    # When sfc heat flux is negative, set depth to 1 m.
    # When sfc heat flux is positive, march up sounding until wpthlp 1st becomes negative.
    depth_pos_wpthlp = jnp.where(
        wpthlp[:, gr.k_lb_zm] <= zero,
        one,
        jnp.maximum(
            one,
            height_path[jnp.arange(ngrdcol), stop_level_pos],
        ),
    )

    # a_const used to be set here; now it is a tunable parameter
    #if ( .not. l_vary_convect_depth ) then
    #   a_const = 1.8_core_rknd
    #else
    #   a_const = 0.6_core_rknd
    #end if

    if l_andre_1978:
        raise NotImplementedError(
            "l_andre_1978 is a false Fortran parameter and is not active."
        )
    else:

        # Compute ustar^2
        ustar2 = jnp.sqrt(upwp_sfc * upwp_sfc + vpwp_sfc * vpwp_sfc)

        # Compute wstar following Andre et al., 1976
        if not l_vary_convect_depth:
            wstar = (
                one / T0
                * grav
                * jnp.maximum(wpthlp[:, gr.k_lb_zm], zero)
                * z_const
            ) ** one_third
        else:
            wstar = (
                one / T0
                * grav
                * jnp.maximum(wpthlp[:, gr.k_lb_zm], zero)
                * 0.2
                * depth_pos_wpthlp
            ) ** one_third
        wstar = jnp.where(wpthlp[:, gr.k_lb_zm] > zero, wstar, zero)

        # Surface friction velocity following Andre et al. 1978
        if not l_vary_convect_depth:
            uf = jnp.sqrt(ustar2 + 0.3 * wstar * wstar)
        else:
            uf = jnp.sqrt(ustar2 + wstar * wstar)

        uf = jnp.maximum(ufmin, uf)

        # Compute estimate for surface second order moments
        wp2_sfc = a_const * uf ** 2
        up2_sfc = up2_sfc_coef * a_const * uf ** 2  # From Andre, et al. 1978
        vp2_sfc = up2_sfc_coef * a_const * uf ** 2  # "  "

        # Notes:  1) With "a" having a value of 1.8, the surface correlations of
        #            both w & rt and w & thl have a value of about 0.878.
        #         2) The surface correlation of rt & thl is 0.5.
        # Brian Griffin; February 2, 2008.

        if not l_vary_convect_depth:
            thlp2_sfc = 0.4 * a_const * (wpthlp[:, gr.k_lb_zm] / uf) ** 2
            rtp2_sfc = 0.4 * a_const * (wprtp_sfc / uf) ** 2
            rtpthlp_sfc = (
                0.2
                * a_const
                * (wpthlp[:, gr.k_lb_zm] / uf)
                * (wprtp_sfc / uf)
            )
        else:
            thlp2_sfc = (
                (wpthlp[:, gr.k_lb_zm] / uf) ** 2
                / (max_mag_correlation_flux ** 2 * a_const)
            )
            rtp2_sfc = (
                (wprtp_sfc / uf) ** 2
                / (max_mag_correlation_flux ** 2 * a_const)
            )
            rtpthlp_sfc = (
                max_mag_correlation_flux
                * jnp.sqrt(thlp2_sfc * rtp2_sfc)
            )

        thlp2_sfc = jnp.maximum(thl_tol ** 2, thlp2_sfc)
        rtp2_sfc = jnp.maximum(rt_tol ** 2, rtp2_sfc)

        # Add effect of vertical compression of eddies on horizontal gustiness.
        # First, ensure that wp2 does not make the correlation
        #   of (w,thl) or (w,rt) outside (-1,1).
        # Perhaps in the future we should also ensure that the correlations
        #   of (w,u) and (w,v) are not outside (-1,1).
        min_wp2_sfc_val = jnp.maximum(
            w_tol_sqd,
            jnp.maximum(
                wprtp_sfc ** 2
                / (rtp2_sfc * max_mag_correlation_flux ** 2),
                wpthlp[:, gr.k_lb_zm] ** 2
                / (thlp2_sfc * max_mag_correlation_flux ** 2),
            ),
        )

        splat_wp2 = (
            tau_zm[:, gr.k_lb_zm]
            * lhs_splat_wp2[:, gr.k_lb_zm]
            * wp2_sfc
        )
        splatting_drives_wp2_small = (
            wp2_sfc - splat_wp2 < min_wp2_sfc_val
        )
        wp2_splat_sfc_correction = jnp.where(
            splatting_drives_wp2_small,
            -wp2_sfc + min_wp2_sfc_val,
            splat_wp2,
        )
        wp2_sfc = jnp.where(
            splatting_drives_wp2_small,
            min_wp2_sfc_val,
            wp2_sfc + wp2_splat_sfc_correction,
        )

        up2_sfc = up2_sfc - one_half * wp2_splat_sfc_correction
        vp2_sfc = vp2_sfc - one_half * wp2_splat_sfc_correction

        # Passive scalars
        for sclr in range(sclr_dim):
            # Vince Larson changed coeffs to make correlations between [-1,1].
            # 31 Jan 2008
            # sclrprtp(i) =
            # a * (wprtp_sfc / uf) * (wpsclrp(i) / uf)
            # sclrpthlp(i) =
            # a * (wpthlp / uf) * (wpsclrp(i) / uf)
            # sclrp2(i) =
            # sclr_var_coef * a * ( wpsclrp(i) / uf )**2
            # Notes:  1) With "a" having a value of 1.8 and "sclr_var_coef"
            #            having a value of 0.4, the surface correlation of
            #            w & sclr has a value of about 0.878.
            #         2) With "sclr_var_coef" having a value of 0.4, the
            #            surface correlations of both rt & sclr and
            #            thl & sclr are 0.5.
            # Brian Griffin; February 2, 2008.

            # We use the following if..then's to make sclr_rt and sclr_thl
            # close to the actual thlp2/rtp2 at the surface.
            # -dschanen 25 Sep 08
            if int(sclr_idx.iisclr_rt) == sclr + 1:
                # If we are trying to emulate rt with the scalar, then we
                # use the variance coefficient from above
                sclrprtp_sfc = (
                    0.4
                    * a_const
                    * (wprtp_sfc / uf)
                    * (wpsclrp_sfc[:, sclr] / uf)
                )
            else:
                sclrprtp_sfc = (
                    0.2
                    * a_const
                    * (wprtp_sfc / uf)
                    * (wpsclrp_sfc[:, sclr] / uf)
                )

            if int(sclr_idx.iisclr_thl) == sclr + 1:
                # As above, but for thetal
                sclrpthlp_sfc = (
                    0.4
                    * a_const
                    * (wpthlp[:, gr.k_lb_zm] / uf)
                    * (wpsclrp_sfc[:, sclr] / uf)
                )
            else:
                sclrpthlp_sfc = (
                    0.2
                    * a_const
                    * (wpthlp[:, gr.k_lb_zm] / uf)
                    * (wpsclrp_sfc[:, sclr] / uf)
                )

            sclrp2_sfc = (
                sclr_var_coef
                * a_const
                * (wpsclrp_sfc[:, sclr] / uf) ** 2
            )

            sclrprtp = sclrprtp.at[:, gr.k_lb_zm, sclr].set(sclrprtp_sfc)
            sclrpthlp = sclrpthlp.at[:, gr.k_lb_zm, sclr].set(sclrpthlp_sfc)
            sclrp2 = sclrp2.at[:, gr.k_lb_zm, sclr].set(sclrp2_sfc)

            # End Vince Larson's change.

    # Clip wp2 at wp2_max, same as in advance_wp2_wp3
    wp2_sfc = jnp.minimum(wp2_sfc, wp2_max)

    not_at_surface = (
        jnp.abs(gr.zm[:, gr.k_lb_zm] - sfc_elevation)
        > jnp.abs(gr.zm[:, gr.k_lb_zm] + sfc_elevation) * eps / 2
    )

    # Variances for cases where the lowest level is not at the surface.
    # Eliminate surface effects on lowest level variances.
    wp2_sfc = jnp.where(not_at_surface, w_tol_sqd, wp2_sfc)
    up2_sfc = jnp.where(not_at_surface, w_tol_sqd, up2_sfc)
    vp2_sfc = jnp.where(not_at_surface, w_tol_sqd, vp2_sfc)
    thlp2_sfc = jnp.where(not_at_surface, thl_tol ** 2, thlp2_sfc)
    rtp2_sfc = jnp.where(not_at_surface, rt_tol ** 2, rtp2_sfc)
    rtpthlp_sfc = jnp.where(not_at_surface, zero, rtpthlp_sfc)

    wp2 = wp2.at[:, gr.k_lb_zm].set(wp2_sfc)
    up2 = up2.at[:, gr.k_lb_zm].set(up2_sfc)
    vp2 = vp2.at[:, gr.k_lb_zm].set(vp2_sfc)
    thlp2 = thlp2.at[:, gr.k_lb_zm].set(thlp2_sfc)
    rtp2 = rtp2.at[:, gr.k_lb_zm].set(rtp2_sfc)
    rtpthlp = rtpthlp.at[:, gr.k_lb_zm].set(rtpthlp_sfc)

    for sclr in range(sclr_dim):
        sclrp2 = sclrp2.at[:, gr.k_lb_zm, sclr].set(
            jnp.where(not_at_surface, zero, sclrp2[:, gr.k_lb_zm, sclr])
        )
        sclrprtp = sclrprtp.at[:, gr.k_lb_zm, sclr].set(
            jnp.where(not_at_surface, zero, sclrprtp[:, gr.k_lb_zm, sclr])
        )
        sclrpthlp = sclrpthlp.at[:, gr.k_lb_zm, sclr].set(
            jnp.where(not_at_surface, zero, sclrpthlp[:, gr.k_lb_zm, sclr])
        )

    if stats.l_sample:
        stats = stats.finalize_budget("thlp2_sf", thlp2 / dt)
        stats = stats.finalize_budget("rtp2_sf", rtp2 / dt)
        stats = stats.finalize_budget("rtpthlp_sf", rtpthlp / dt)
        stats = stats.finalize_budget("up2_sf", up2 / dt)
        stats = stats.finalize_budget("vp2_sf", vp2 / dt)
        stats = stats.finalize_budget("wp2_sf", wp2 / dt)

    if clubb_at_least_debug_level(2):
        nan_mask = (
            jnp.isnan(wp2[:, gr.k_lb_zm])
            | jnp.isnan(up2[:, gr.k_lb_zm])
            | jnp.isnan(vp2[:, gr.k_lb_zm])
            | jnp.isnan(thlp2[:, gr.k_lb_zm])
            | jnp.isnan(rtp2[:, gr.k_lb_zm])
            | jnp.isnan(rtpthlp[:, gr.k_lb_zm])
        )

        for sclr in range(sclr_dim):
            nan_mask = (
                nan_mask
                | jnp.isnan(sclrp2[:, gr.k_lb_zm, sclr])
                | jnp.isnan(sclrprtp[:, gr.k_lb_zm, sclr])
                | jnp.isnan(sclrpthlp[:, gr.k_lb_zm, sclr])
            )

        err_info = err_info.set_fatal(nan_mask)

    return (
        wp2, up2, vp2,
        thlp2, rtp2, rtpthlp,
        sclrp2, sclrprtp, sclrpthlp,
        err_info, stats,
    )
