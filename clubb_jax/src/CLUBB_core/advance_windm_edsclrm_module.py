"""JAX-side entry point for ``src/CLUBB_core/advance_windm_edsclrm_module.F90``.

Description:
  Solves for both mean horizontal wind components, um and vm, and for the
  eddy-scalars (passive scalars that don't use the high-order closure).

Uses the LAPACK tridiagonal solver subroutine with 2 + # of scalar(s)
back substitutions (since the left hand side matrix is the same for all
input variables).

References:
  Eqn. 8 & 9 on p. 3545 of
  ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    Method and Model Description'' Golaz, et al. (2002)
    JAS, Vol. 59, pp. 3540--3551.

Porting deviations:
- Sponge damping blocks are unsupported here because clubb_case_initalization
  rejects sponge-enabled Python/JAX driver cases before this routine is called.
- The detailed Fortran error-print diagnostics and diagnostic-only early
  returns are reduced until full JAX diagnostic state is available; fatal
  conditions still mark err_info.
- Stats are threaded explicitly with JaxStats because the Fortran routine uses
  global stats side effects that are not JAX-compatible state.
- Fortran intent(out)/intent(inout) arguments are returned explicitly.
"""

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_helper_module import calc_xpwp, pvertinterp
from clubb_jax.src.CLUBB_core.clip_explicit import (
    clip_covar,
    clip_upwp,
    clip_vpwp,
    upwp_cl_max,
    vpwp_cl_max,
)
from clubb_jax.src.CLUBB_core.constants_clubb import (
    eps,
    one,
    one_half,
    zero,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.model_flags import l_force_descending_solves
from clubb_jax.src.CLUBB_core.parameter_indices import ic_K10, ic_K10h
from clubb_jax.src.CLUBB_core.diffusion import diffusion_zt_lhs
from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_vertical
from clubb_jax.src.CLUBB_core.grid_class import zm2zt
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.matrix_solver_wrapper import tridiag_solve
from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs
from clubb_jax.src.CLUBB_core import ErrInfo, Grid, NuVertResDep


windm_edsclrm_um = 1
"""Named constant to handle um solves."""
windm_edsclrm_vm = 2
"""Named constant to handle vm solves."""
windm_edsclrm_scalar = 3
"""Named constant to handle scalar solves."""

ndiags3 = 3


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "edsclr_dim",
        "tridiag_solve_method",
        "l_predict_upwp_vpwp",
        "l_upwind_xm_ma",
        "l_uv_nudge",
        "l_lmm_stepping",
        "l_linearize_pbl_winds",
        "l_do_expldiff_rtm_thlm",
        "fill_holes_type",
        "upwp_cl_num",
        "vpwp_cl_num",
        "l_implemented",
    ),
)
def advance_windm_edsclrm(
    nzm: int, nzt: int, ngrdcol: int, edsclr_dim: int, gr: Grid, dt,
    wm_zt, Kh_zm, clubb_params,
    ug, vg, um_ref, vm_ref,
    wp2, up2, vp2, um_forcing, vm_forcing,
    edsclrm_forcing, p_in_Pa,
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zt,
    fcor, l_implemented: bool,
    nu_vert_res_dep: NuVertResDep, ts_nudge,
    tridiag_solve_method: int,
    l_predict_upwp_vpwp: bool,
    l_upwind_xm_ma: bool,
    l_uv_nudge: bool,
    l_lmm_stepping: bool,
    l_linearize_pbl_winds: bool,
    l_do_expldiff_rtm_thlm: bool,
    fill_holes_type: int,
    upwp_cl_num: int,
    vpwp_cl_num: int,
    stats: JaxStats,
    um, vm, thlm, rtm, edsclrm,
    upwp, vpwp, wpedsclrp,
    um_pert, vm_pert, upwp_pert, vpwp_pert,
    err_info: ErrInfo,
):
    """Solves for both mean horizontal wind components, um and vm, and for the
    eddy-scalars (passive scalars that don't use the high-order closure).

    Uses the LAPACK tridiagonal solver subroutine with 2 + # of scalar(s)
    back substitutions (since the left hand side matrix is the same for all
    input variables).
    """
    nu_zero = jnp.zeros((ngrdcol,), dtype=jnp.float64)

    # Coefficient for momentum
    Km_zm = Kh_zm * clubb_params[:, ic_K10][:, None]
    # Coefficient for thermo
    Kmh_zm = Kh_zm * clubb_params[:, ic_K10h][:, None]
    Km_zm_p_nu10 = Km_zm + nu_vert_res_dep.nu10[:, None]

    if edsclr_dim > 1 and l_do_expldiff_rtm_thlm:
        edsclrm = edsclrm.at[:, :, edsclr_dim - 2].set(thlm)
        edsclrm = edsclrm.at[:, :, edsclr_dim - 1].set(rtm)

    l_perturbed_wind = (not l_predict_upwp_vpwp) and l_linearize_pbl_winds

    if not l_implemented:
        lhs_ma_zt = term_ma_zt_lhs(
            nzm, nzt, ngrdcol, wm_zt, gr.weights_zt2zm,
            gr.invrs_dzt, gr.invrs_dzm, l_upwind_xm_ma, gr.grid_dir,
        )
    else:
        lhs_ma_zt = jnp.zeros((ndiags3, ngrdcol, nzt), dtype=jnp.float64)

    lhs = jnp.zeros((ndiags3, ngrdcol, nzt), dtype=jnp.float64)
    rhs = jnp.zeros((ngrdcol, nzt, max(2, edsclr_dim)), dtype=jnp.float64)
    solution = jnp.zeros_like(rhs)
    upwp_chnge = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    vpwp_chnge = jnp.zeros((ngrdcol, nzm), dtype=jnp.float64)
    wind_speed = jnp.ones((ngrdcol, nzt), dtype=jnp.float64)
    u_star_sqd = jnp.zeros((ngrdcol,), dtype=jnp.float64)

    if not l_predict_upwp_vpwp:
        Km_zt = zm2zt(nzm, nzt, ngrdcol, gr, Km_zm, zero)

        # Calculate diffusion terms
        lhs_diff = diffusion_zt_lhs(
            nzm, nzt, ngrdcol, gr, Km_zm, Km_zt, nu_vert_res_dep.nu10,
            invrs_rho_ds_zt, rho_ds_zm,
        )

        if l_lmm_stepping:
            um_old = um
            vm_old = vm
        else:
            um_old = jnp.zeros_like(um)
            vm_old = jnp.zeros_like(vm)

        # ----------------------------------------------------------------
        # Prepare tridiagonal system for horizontal winds, um and vm
        # ----------------------------------------------------------------

        # Compute Coriolis, geostrophic, and other prescribed wind forcings for um.
        um_tndcy, stats = compute_uv_tndcy(
            nzt, ngrdcol, windm_edsclrm_um,
            fcor, vm, vg,
            um_forcing, l_implemented,
            stats,
        )

        # Compute Coriolis, geostrophic, and other prescribed wind forcings for vm.
        vm_tndcy, stats = compute_uv_tndcy(
            nzt, ngrdcol, windm_edsclrm_vm,
            fcor, um, ug,
            vm_forcing, l_implemented,
            stats,
        )

        # Momentum surface fluxes, u'w'|_sfc and v'w'|_sfc, are applied through
        # an implicit method, such that:
        #    x'w'|_sfc = - ( u_star(t)^2 / wind_speed(t) ) * xm(t+1).
        l_imp_sfc_momentum_flux = True

        # Compute wind speed (use threshold "eps" to prevent divide-by-zero
        # error).
        wind_speed = jnp.maximum(jnp.sqrt(um ** 2 + vm ** 2), eps)
        # Compute u_star_sqd according to the definition of u_star.
        u_star_sqd = jnp.sqrt(
            upwp[:, gr.k_lb_zm] ** 2 + vpwp[:, gr.k_lb_zm] ** 2
        )

        # Compute the explicit portion of the um equation.
        # Build the right-hand side vector.
        rhs_um, stats = windm_edsclrm_rhs(
            nzm, nzt, ngrdcol, gr, windm_edsclrm_um, dt,
            lhs_diff, um, um_tndcy,
            rho_ds_zm, invrs_rho_ds_zt,
            l_imp_sfc_momentum_flux, upwp[:, gr.k_lb_zm],
            stats,
        )
        rhs = rhs.at[:, :, windm_edsclrm_um - 1].set(rhs_um)

        # Compute the explicit portion of the vm equation.
        # Build the right-hand side vector.
        rhs_vm, stats = windm_edsclrm_rhs(
            nzm, nzt, ngrdcol, gr, windm_edsclrm_vm, dt,
            lhs_diff, vm, vm_tndcy,
            rho_ds_zm, invrs_rho_ds_zt,
            l_imp_sfc_momentum_flux, vpwp[:, gr.k_lb_zm],
            stats,
        )
        rhs = rhs.at[:, :, windm_edsclrm_vm - 1].set(rhs_vm)

        # Store momentum flux (explicit component)
        xpwp = calc_xpwp(nzm, nzt, ngrdcol, gr, Km_zm_p_nu10, um)
        # Solve for x'w' at all intermediate model levels.
        # A Crank-Nicholson timestep is used.
        upwp = upwp.at[:, 1:nzm - 1].set(-one_half * xpwp[:, 1:nzm - 1])

        xpwp = calc_xpwp(nzm, nzt, ngrdcol, gr, Km_zm_p_nu10, vm)
        vpwp = vpwp.at[:, 1:nzm - 1].set(-one_half * xpwp[:, 1:nzm - 1])

        # A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
        # means that x'w' at the top model level is 0,
        # since x'w' = - K_zm * d(xm)/dz.
        upwp = upwp.at[:, gr.k_ub_zm].set(zero)
        vpwp = vpwp.at[:, gr.k_ub_zm].set(zero)

        # Compute the implicit portion of the um and vm equations.
        # Build the left-hand side matrix.
        lhs = windm_edsclrm_lhs(
            nzm, nzt, ngrdcol, gr, dt,
            lhs_ma_zt, lhs_diff,
            wind_speed, u_star_sqd,
            rho_ds_zm, invrs_rho_ds_zt,
            l_implemented, l_imp_sfc_momentum_flux,
        )

        # Decompose and back substitute for um and vm
        nrhs = 2
        l_need_rcond = bool(
            stats.l_sample and stats.var_on_stats_list("windm_matrix_condt_num")
        )
        solution_wind, err_info, stats = windm_edsclrm_solve(
            nzt, ngrdcol, gr, nrhs,
            tridiag_solve_method,
            l_implemented,
            stats, l_need_rcond,
            lhs, rhs[:, :, :nrhs], err_info,
        )
        solution = solution.at[:, :, :nrhs].set(solution_wind)

        # ----------------------------------------------------------------
        # Update zonal (west-to-east) component of mean wind, um
        # ----------------------------------------------------------------
        um = solution[:, :, windm_edsclrm_um - 1]
        # ----------------------------------------------------------------
        # Update meridional (south-to-north) component of mean wind, vm
        # ----------------------------------------------------------------
        vm = solution[:, :, windm_edsclrm_vm - 1]

        if stats.l_sample:
            stats = windm_edsclrm_implicit_stats(
                nzm, nzt, ngrdcol, windm_edsclrm_um, gr,
                um, gr.invrs_dzt,
                lhs_diff, lhs_ma_zt,
                invrs_rho_ds_zt, u_star_sqd,
                rho_ds_zm, wind_speed,
                l_imp_sfc_momentum_flux,
                stats,
            )

            stats = windm_edsclrm_implicit_stats(
                nzm, nzt, ngrdcol, windm_edsclrm_vm, gr,
                vm, gr.invrs_dzt,
                lhs_diff, lhs_ma_zt,
                invrs_rho_ds_zt, u_star_sqd,
                rho_ds_zm, wind_speed,
                l_imp_sfc_momentum_flux,
                stats,
            )

        if l_lmm_stepping:
            um = one_half * (um_old + um)
            vm = one_half * (vm_old + vm)

        # Second part of momentum (implicit component)

        # Solve for x'w' at all intermediate model levels.
        # A Crank-Nicholson timestep is used.
        xpwp = calc_xpwp(nzm, nzt, ngrdcol, gr, Km_zm_p_nu10, um)
        upwp = upwp.at[:, 1:nzm - 1].add(-one_half * xpwp[:, 1:nzm - 1])

        xpwp = calc_xpwp(nzm, nzt, ngrdcol, gr, Km_zm_p_nu10, vm)
        vpwp = vpwp.at[:, 1:nzm - 1].add(-one_half * xpwp[:, 1:nzm - 1])

        # Adjust um and vm if nudging is turned on.
        if l_uv_nudge:
            # Reflect nudging in budget
            if stats.l_sample:
                stats = stats.begin_budget("um_ndg", um / dt)
                stats = stats.begin_budget("vm_ndg", vm / dt)

            um = um - ((um - um_ref) * (dt / ts_nudge))
            vm = vm - ((vm - vm_ref) * (dt / ts_nudge))

            if stats.l_sample:
                stats = stats.finalize_budget("um_ndg", um / dt)
                stats = stats.finalize_budget("vm_ndg", vm / dt)

        if stats.l_sample:
            stats = stats.update("um_ref", um_ref)
            stats = stats.update("vm_ref", vm_ref)

        # Clipping for u'w'
        #
        # Clipping u'w' at each vertical level, based on the
        # correlation of u and w at each vertical level, such that:
        # corr_(u,w) = u'w' / [ sqrt(u'^2) * sqrt(w'^2) ];
        # -1 <= corr_(u,w) <= 1.
        #
        # Since u'^2, w'^2, and u'w' are each advanced in different
        # subroutines from each other in advance_clubb_core, clipping for u'w'
        # has to be done three times during each timestep (once after each
        # variable has been updated).
        # This clip can be the first, middle, or last budget contribution,
        # depending on the relative advancement order in this timestep.
        if stats.l_sample and l_predict_upwp_vpwp:
            if upwp_cl_num == 0:
                stats = stats.begin_budget("upwp_cl", upwp / dt)
            else:
                stats = stats.update_budget("upwp_cl", -upwp / dt)
        upwp_cl_num = upwp_cl_num + 1
        upwp, upwp_chnge = clip_covar(
            nzm, ngrdcol, clip_upwp, wp2, up2, upwp,
        )
        if stats.l_sample and l_predict_upwp_vpwp:
            if upwp_cl_num == upwp_cl_max:
                stats = stats.finalize_budget("upwp_cl", upwp / dt)
            else:
                stats = stats.update_budget("upwp_cl", upwp / dt)

        # Clipping for v'w'
        #
        # Clipping v'w' at each vertical level, based on the
        # correlation of v and w at each vertical level, such that:
        # corr_(v,w) = v'w' / [ sqrt(v'^2) * sqrt(w'^2) ];
        # -1 <= corr_(v,w) <= 1.
        #
        # Since v'^2, w'^2, and v'w' are each advanced in different
        # subroutines from each other in advance_clubb_core, clipping for v'w'
        # has to be done three times during each timestep (once after each
        # variable has been updated).
        # This clip can be the first, middle, or last budget contribution,
        # depending on the relative advancement order in this timestep.
        if stats.l_sample and l_predict_upwp_vpwp:
            if vpwp_cl_num == 0:
                stats = stats.begin_budget("vpwp_cl", vpwp / dt)
            else:
                stats = stats.update_budget("vpwp_cl", -vpwp / dt)
        vpwp_cl_num = vpwp_cl_num + 1
        vpwp, vpwp_chnge = clip_covar(
            nzm, ngrdcol, clip_vpwp, wp2, vp2, vpwp,
        )
        if stats.l_sample and l_predict_upwp_vpwp:
            if vpwp_cl_num == vpwp_cl_max:
                stats = stats.finalize_budget("vpwp_cl", vpwp / dt)
            else:
                stats = stats.update_budget("vpwp_cl", vpwp / dt)

    if l_perturbed_wind:
        l_imp_sfc_momentum_flux = True

        # Compute wind speed (use threshold "eps" to prevent divide-by-zero
        # error).
        wind_speed_pert = jnp.maximum(jnp.sqrt(um_pert ** 2 + vm_pert ** 2), eps)
        # Compute u_star_sqd according to the definition of u_star.
        u_star_sqd_pert = jnp.sqrt(
            upwp_pert[:, gr.k_lb_zm] ** 2 + vpwp_pert[:, gr.k_lb_zm] ** 2
        )

        # Compute the explicit portion of the um equation.
        # Build the right-hand side vector.
        rhs_um, stats = windm_edsclrm_rhs(
            nzm, nzt, ngrdcol, gr, windm_edsclrm_um, dt,
            lhs_diff, um_pert, um_tndcy,
            rho_ds_zm, invrs_rho_ds_zt,
            l_imp_sfc_momentum_flux, upwp_pert[:, gr.k_lb_zm],
            stats,
        )
        rhs = rhs.at[:, :, windm_edsclrm_um - 1].set(rhs_um)

        # Compute the explicit portion of the vm equation.
        # Build the right-hand side vector.
        rhs_vm, stats = windm_edsclrm_rhs(
            nzm, nzt, ngrdcol, gr, windm_edsclrm_vm, dt,
            lhs_diff, vm_pert, vm_tndcy,
            rho_ds_zm, invrs_rho_ds_zt,
            l_imp_sfc_momentum_flux, vpwp_pert[:, gr.k_lb_zm],
            stats,
        )
        rhs = rhs.at[:, :, windm_edsclrm_vm - 1].set(rhs_vm)

        # Store momentum flux (explicit component)
        xpwp = calc_xpwp(nzm, nzt, ngrdcol, gr, Km_zm_p_nu10, um_pert)
        # Solve for x'w' at all intermediate model levels.
        # A Crank-Nicholson timestep is used.
        upwp_pert = upwp_pert.at[:, 1:nzm - 1].set(-one_half * xpwp[:, 1:nzm - 1])

        xpwp = calc_xpwp(nzm, nzt, ngrdcol, gr, Km_zm_p_nu10, vm_pert)
        vpwp_pert = vpwp_pert.at[:, 1:nzm - 1].set(-one_half * xpwp[:, 1:nzm - 1])

        # A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
        # means that x'w' at the top model level is 0,
        # since x'w' = - K_zm * d(xm)/dz.
        upwp_pert = upwp_pert.at[:, gr.k_ub_zm].set(zero)
        vpwp_pert = vpwp_pert.at[:, gr.k_ub_zm].set(zero)

        # Compute the implicit portion of the um and vm equations.
        # Build the left-hand side matrix.
        lhs = windm_edsclrm_lhs(
            nzm, nzt, ngrdcol, gr, dt,
            lhs_ma_zt, lhs_diff,
            wind_speed_pert, u_star_sqd_pert,
            rho_ds_zm, invrs_rho_ds_zt,
            l_implemented, l_imp_sfc_momentum_flux,
        )

        # Decompose and back substitute for um and vm
        nrhs = 2
        l_need_rcond = bool(
            stats.l_sample and stats.var_on_stats_list("windm_matrix_condt_num")
        )
        solution_wind, err_info, stats = windm_edsclrm_solve(
            nzt, ngrdcol, gr, nrhs,
            tridiag_solve_method,
            l_implemented,
            stats, l_need_rcond,
            lhs, rhs[:, :, :nrhs], err_info,
        )
        solution = solution.at[:, :, :nrhs].set(solution_wind)

        um_pert = solution[:, :, windm_edsclrm_um - 1]
        vm_pert = solution[:, :, windm_edsclrm_vm - 1]

        # Second part of momentum (implicit component)

        # Solve for x'w' at all intermediate model levels.
        # A Crank-Nicholson timestep is used.
        xpwp = calc_xpwp(nzm, nzt, ngrdcol, gr, Km_zm_p_nu10, um_pert)
        upwp_pert = upwp_pert.at[:, 1:nzm - 1].add(-one_half * xpwp[:, 1:nzm - 1])

        xpwp = calc_xpwp(nzm, nzt, ngrdcol, gr, Km_zm_p_nu10, vm_pert)
        vpwp_pert = vpwp_pert.at[:, 1:nzm - 1].add(-one_half * xpwp[:, 1:nzm - 1])

        # Clipping for u'w' and v'w'.
        upwp_pert, upwp_chnge = clip_covar(
            nzm, ngrdcol, clip_upwp, wp2, up2, upwp_pert,
        )
        vpwp_pert, vpwp_chnge = clip_covar(
            nzm, ngrdcol, clip_vpwp, wp2, vp2, vpwp_pert,
        )

    if edsclr_dim > 0:
        Kmh_zt = zm2zt(nzm, nzt, ngrdcol, gr, Kmh_zm, zero)

        # Calculate diffusion terms
        lhs_diff = diffusion_zt_lhs(
            nzm, nzt, ngrdcol, gr, Kmh_zm, Kmh_zt, nu_zero,
            invrs_rho_ds_zt, rho_ds_zm,
        )

        if l_lmm_stepping:
            edsclrm_old = edsclrm
        else:
            edsclrm_old = jnp.zeros_like(edsclrm)

        l_imp_sfc_momentum_flux = False

        for edsclr in range(edsclr_dim):
            # Compute the explicit portion of the scalar equation.
            # Build the right-hand side vector.
            rhs_scalar, stats = windm_edsclrm_rhs(
                nzm, nzt, ngrdcol, gr,
                windm_edsclrm_scalar, dt,
                lhs_diff, edsclrm[:, :, edsclr],
                edsclrm_forcing[:, :, edsclr],
                rho_ds_zm, invrs_rho_ds_zt,
                l_imp_sfc_momentum_flux,
                wpedsclrp[:, gr.k_lb_zm, edsclr],
                stats,
            )
            rhs = rhs.at[:, :, edsclr].set(rhs_scalar)

        for edsclr in range(edsclr_dim):
            # Store eddy scalar flux (explicit component)
            xpwp = calc_xpwp(
                nzm, nzt, ngrdcol, gr,
                Km_zm_p_nu10, edsclrm[:, :, edsclr],
            )
            wpedsclrp = wpedsclrp.at[:, 1:nzm - 1, edsclr].set(
                -one_half * xpwp[:, 1:nzm - 1]
            )

        # A zero-flux boundary condition at the top of the model, d(xm)/dz = 0,
        # means that x'w' at the top model level is 0,
        # since x'w' = - K_zm * d(xm)/dz.
        wpedsclrp = wpedsclrp.at[:, gr.k_ub_zm, :edsclr_dim].set(zero)

        # Compute the implicit portion of the scalar equations.
        # Build the left-hand side matrix.
        lhs = windm_edsclrm_lhs(
            nzm, nzt, ngrdcol, gr, dt,
            lhs_ma_zt, lhs_diff,
            wind_speed, u_star_sqd,
            rho_ds_zm, invrs_rho_ds_zt,
            l_implemented, l_imp_sfc_momentum_flux,
        )

        l_need_rcond = False
        solution_scalar, err_info, stats = windm_edsclrm_solve(
            nzt, ngrdcol, gr, edsclr_dim,
            tridiag_solve_method,
            l_implemented,
            stats, l_need_rcond,
            lhs, rhs[:, :, :edsclr_dim], err_info,
        )
        solution = solution.at[:, :, :edsclr_dim].set(solution_scalar)

        # Update eddy scalar quantities.
        edsclrm = edsclrm.at[:, :, :edsclr_dim].set(solution[:, :, :edsclr_dim])

        if l_lmm_stepping:
            edsclrm = edsclrm.at[:, :, :edsclr_dim].set(
                one_half * (edsclrm_old[:, :, :edsclr_dim] + edsclrm[:, :, :edsclr_dim])
            )

        for edsclr in range(edsclr_dim):
            # Second part of eddy scalar flux (implicit component)
            xpwp = calc_xpwp(
                nzm, nzt, ngrdcol, gr,
                Kmh_zm, edsclrm[:, :, edsclr],
            )
            wpedsclrp = wpedsclrp.at[:, 1:nzm - 1, edsclr].set(
                -one_half * xpwp[:, 1:nzm - 1]
            )

    if edsclr_dim > 1 and l_do_expldiff_rtm_thlm:
        thlm700 = pvertinterp(nzt, ngrdcol, gr, p_in_Pa, 70000.0, thlm)
        thlm1000 = pvertinterp(nzt, ngrdcol, gr, p_in_Pa, 100000.0, thlm)
        apply_explicit_diffusion = (thlm700 - thlm1000) < 20.0

        if stats.l_sample:
            thlm_ed = jnp.where(
                apply_explicit_diffusion[:, None],
                (edsclrm[:, :, edsclr_dim - 2] - thlm) / dt,
                zero,
            )
            rtm_ed = jnp.where(
                apply_explicit_diffusion[:, None],
                (edsclrm[:, :, edsclr_dim - 1] - rtm) / dt,
                zero,
            )
        else:
            thlm_ed = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
            rtm_ed = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)

        thlm = jnp.where(
            apply_explicit_diffusion[:, None],
            edsclrm[:, :, edsclr_dim - 2],
            thlm,
        )
        rtm = jnp.where(
            apply_explicit_diffusion[:, None],
            edsclrm[:, :, edsclr_dim - 1],
            rtm,
        )
    else:
        thlm_ed = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        rtm_ed = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)

    if stats.l_sample:
        stats = stats.update("thlm_ed", thlm_ed)
        stats = stats.update("rtm_ed", rtm_ed)

    if edsclr_dim > 0 and fill_holes_type != 0:
        for edsclr in range(edsclr_dim):
            edsclrm_filled = fill_holes_vertical(
                nzt, ngrdcol, zero_threshold,
                gr.k_lb_zt, gr.k_ub_zt,
                gr.dzt, rho_ds_zt, gr.grid_dir_indx,
                fill_holes_type,
                edsclrm[:, :, edsclr],
            )
            edsclrm = edsclrm.at[:, :, edsclr].set(edsclrm_filled)

    return (
        upwp_cl_num, vpwp_cl_num,
        um, vm, thlm, rtm, edsclrm, upwp, vpwp, wpedsclrp,
        um_pert, vm_pert, upwp_pert, vpwp_pert, err_info, stats,
    )


@partial(
    jax.jit,
    static_argnames=(
        "nzt",
        "ngrdcol",
        "nrhs",
        "tridiag_solve_method",
        "l_implemented",
        "l_need_rcond",
    ),
)
def windm_edsclrm_solve(
    nzt, ngrdcol, gr, nrhs,
    tridiag_solve_method,
    l_implemented,
    stats, l_need_rcond,
    lhs, rhs, err_info,
):
    """Solves the horizontal wind or eddy-scalar time-tendency equation, and
    diagnoses the turbulent flux.  A Crank-Nicholson time-stepping algorithm
    is used in solving the turbulent advection term and in diagnosing the
    turbulent flux.

    The rate of change of an eddy-scalar quantity, xm, is:

    d(xm)/dt = - w * d(xm)/dz - (1/rho_ds) * d( rho_ds * x'w' )/dz
               + xm_forcings.

    References:
    Eqn. 8 & 9 on p. 3545 of
    ``A PDF-Based Model for Boundary Layer Clouds. Part I:
      Method and Model Description'' Golaz, et al. (2002)
    JAS, Vol. 59, pp. 3540--3551.
    """
    lhs_solve = lhs
    rhs_solve = rhs

    # Matrix solves are bit-different between ascending and descending.
    # This ensures matrices are solved in the same (descending) order,
    # which is useful for ensuring BFBness between grid modes
    if l_force_descending_solves and gr.grid_dir_indx > 0:
        # We need to flip in both the vertical dimensions and the bands for the lhs
        lhs_solve = lhs_solve[::-1, :, ::-1]
        rhs_solve = rhs_solve[:, ::-1, :]

    err_info, solution, rcond = tridiag_solve(
        "windm_edsclrm", tridiag_solve_method, ngrdcol, nzt,
        lhs_solve, rhs_solve, err_info,
        use_rcond=l_need_rcond,
        l_implemented=l_implemented,
    )

    if l_need_rcond and stats.l_sample:
        stats = stats.update("windm_matrix_condt_num", one / rcond)

    # Flip the back to the ascending direction if we forced the solve
    # to be in descending mode
    if l_force_descending_solves and gr.grid_dir_indx > 0:
        solution = solution[:, ::-1, :]

    return solution, err_info, stats


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "solve_type",
        "l_imp_sfc_momentum_flux",
    ),
)
def windm_edsclrm_implicit_stats(
    nzm, nzt, ngrdcol, solve_type, gr,
    xm, invrs_dzt,
    lhs_diff, lhs_ma_zt,
    invrs_rho_ds_zt, u_star_sqd,
    rho_ds_zm, wind_speed,
    l_imp_sfc_momentum_flux,
    stats,
):
    """Compute implicit contributions to um and vm.

    References:
    None
    """
    del nzm
    if not stats.l_sample:
        return stats

    if solve_type == windm_edsclrm_um:
        name_ma = "um_ma"
        name_ta = "um_ta"
    elif solve_type == windm_edsclrm_vm:
        name_ma = "vm_ma"
        name_ta = "vm_ta"
    else:
        return stats

    imp_sfc_flux = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)

    if l_imp_sfc_momentum_flux and stats.var_on_stats_list(name_ta):
        # Statistics:  implicit contributions for um or vm.
        # xm term ta is modified at level 2 to include the effects of the
        # surface flux.  In this case, this affects the implicit portion of
        # the term, which handles the main diagonal for the turbulent advection term.
        imp_sfc_flux = imp_sfc_flux.at[:, gr.k_lb_zt].set(
            -gr.grid_dir
            * invrs_rho_ds_zt[:, gr.k_lb_zt]
            * invrs_dzt[:, gr.k_lb_zt]
            * rho_ds_zm[:, gr.k_lb_zm]
            * (u_star_sqd / wind_speed[:, gr.k_lb_zt])
        )

    stats_tmp = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    if gr.grid_dir_indx > 0:
        # Lower boundary conditions
        # xm mean advection
        # xm term ma is completely implicit; record with stats_update.
        stats_tmp = stats_tmp.at[:, 0].set(
            -lhs_ma_zt[1, :, 0] * xm[:, 0]
            -lhs_ma_zt[0, :, 0] * xm[:, 1]
        )
        # Interior domain
        # xm mean advection
        # xm term ma is completely implicit; record with stats_update.
        stats_tmp = stats_tmp.at[:, 1:-1].set(
            -lhs_ma_zt[2, :, 1:-1] * xm[:, :-2]
            -lhs_ma_zt[1, :, 1:-1] * xm[:, 1:-1]
            -lhs_ma_zt[0, :, 1:-1] * xm[:, 2:]
        )
        # Upper boundary conditions
        # xm mean advection
        # xm term ma is completely implicit; record with stats_update.
        stats_tmp = stats_tmp.at[:, -1].set(
            -lhs_ma_zt[2, :, -1] * xm[:, -2]
            -lhs_ma_zt[1, :, -1] * xm[:, -1]
        )
    else:
        stats_tmp = stats_tmp.at[:, -1].set(
            -lhs_ma_zt[1, :, -1] * xm[:, -1]
            -lhs_ma_zt[2, :, -1] * xm[:, -2]
        )
        stats_tmp = stats_tmp.at[:, 1:-1].set(
            -lhs_ma_zt[0, :, 1:-1] * xm[:, 2:]
            -lhs_ma_zt[1, :, 1:-1] * xm[:, 1:-1]
            -lhs_ma_zt[2, :, 1:-1] * xm[:, :-2]
        )
        stats_tmp = stats_tmp.at[:, 0].set(
            -lhs_ma_zt[0, :, 0] * xm[:, 1]
            -lhs_ma_zt[1, :, 0] * xm[:, 0]
        )
    stats = stats.update(name_ma, stats_tmp)

    stats_tmp = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    if gr.grid_dir_indx > 0:
        # Lower boundary conditions
        # xm turbulent transport (implicit component)
        # xm term ta has both implicit and explicit components;
        # finalize the implicit part with stats_finalize_budget.
        stats_tmp = stats_tmp.at[:, 0].set(
            (-one_half * lhs_diff[1, :, 0] + imp_sfc_flux[:, 0]) * xm[:, 0]
            -one_half * lhs_diff[0, :, 0] * xm[:, 1]
        )
        # Interior domain
        # xm turbulent transport (implicit component)
        # xm term ta has both implicit and explicit components;
        # finalize the implicit part with stats_finalize_budget.
        stats_tmp = stats_tmp.at[:, 1:-1].set(
            -one_half * lhs_diff[2, :, 1:-1] * xm[:, :-2]
            + (-one_half * lhs_diff[1, :, 1:-1] + imp_sfc_flux[:, 1:-1]) * xm[:, 1:-1]
            -one_half * lhs_diff[0, :, 1:-1] * xm[:, 2:]
        )
        # Upper boundary conditions
        # xm turbulent transport (implicit component)
        # xm term ta has both implicit and explicit components;
        # finalize the implicit part with stats_finalize_budget.
        stats_tmp = stats_tmp.at[:, -1].set(
            -one_half * lhs_diff[2, :, -1] * xm[:, -2]
            + (-one_half * lhs_diff[1, :, -1] + imp_sfc_flux[:, -1]) * xm[:, -1]
        )
    else:
        stats_tmp = stats_tmp.at[:, -1].set(
            (-one_half * lhs_diff[1, :, -1] + imp_sfc_flux[:, -1]) * xm[:, -1]
            -one_half * lhs_diff[2, :, -1] * xm[:, -2]
        )
        stats_tmp = stats_tmp.at[:, 1:-1].set(
            -one_half * lhs_diff[0, :, 1:-1] * xm[:, 2:]
            + (-one_half * lhs_diff[1, :, 1:-1] + imp_sfc_flux[:, 1:-1]) * xm[:, 1:-1]
            -one_half * lhs_diff[2, :, 1:-1] * xm[:, :-2]
        )
        stats_tmp = stats_tmp.at[:, 0].set(
            -one_half * lhs_diff[0, :, 0] * xm[:, 1]
            + (-one_half * lhs_diff[1, :, 0] + imp_sfc_flux[:, 0]) * xm[:, 0]
        )
    # Finalize implicit contributions for xm
    return stats.finalize_budget(name_ta, stats_tmp)


@partial(
    jax.jit,
    static_argnames=("nzt", "ngrdcol", "solve_type", "l_implemented"),
)
def compute_uv_tndcy(
    nzt, ngrdcol, solve_type,
    fcor, perp_wind_m, perp_wind_g,
    xm_forcing, l_implemented,
    stats,
):
    """Computes the explicit tendency for the um and vm wind components.

    The only explicit tendency that is involved in the d(um)/dt or d(vm)/dt
    equations is the Coriolis tendency.

    The d(um)/dt equation contains the term:

    - f * ( v_g - vm );

    where f is the Coriolis parameter and v_g is the v component of the
    geostrophic wind.

    Likewise, the d(vm)/dt equation contains the term:

    + f * ( u_g - um );

    where u_g is the u component of the geostrophic wind.

    This term is treated completely explicitly.  The values of um, vm, u_g,
    and v_g are all found on the thermodynamic levels.

    Wind forcing from the GCSS cases is also added here.
    """
    if not l_implemented:
        # Only compute the Coriolis term if the model is running on it's own,
        # and is not part of a larger, host model.
        if solve_type == windm_edsclrm_um:
            name_gf = "um_gf"
            name_cf = "um_cf"
            name_f = "um_f"
            xm_gf = -fcor[:, None] * perp_wind_g
            xm_cf = fcor[:, None] * perp_wind_m
        elif solve_type == windm_edsclrm_vm:
            name_gf = "vm_gf"
            name_cf = "vm_cf"
            name_f = "vm_f"
            xm_gf = fcor[:, None] * perp_wind_g
            xm_cf = -fcor[:, None] * perp_wind_m
        else:
            name_gf = ""
            name_cf = ""
            name_f = ""
            xm_gf = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
            xm_cf = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)

        xm_tndcy = xm_gf + xm_cf + xm_forcing

        if stats.l_sample:
            # Statistics:  explicit contributions for um or vm.
            # xm term gf is completely explicit; record with stats_update.
            stats = stats.update(name_gf, xm_gf)
            # xm term cf is completely explicit; record with stats_update.
            stats = stats.update(name_cf, xm_cf)
            # xm term F
            stats = stats.update(name_f, xm_forcing)
    else:
        xm_tndcy = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)

    return xm_tndcy, stats


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "l_implemented",
        "l_imp_sfc_momentum_flux",
    ),
)
def windm_edsclrm_lhs(
    nzm, nzt, ngrdcol, gr, dt,
    lhs_ma_zt, lhs_diff,
    wind_speed, u_star_sqd,
    rho_ds_zm, invrs_rho_ds_zt,
    l_implemented, l_imp_sfc_momentum_flux,
):
    """Calculate the implicit portion of the horizontal wind or eddy-scalar
    time-tendency equation.  See the description in subroutine
    windm_edsclrm_solve for more details.

      --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
      Simple changes to this procedure may adversely affect computational speed
          - Gunther Huebler, Aug. 2018, clubb:ticket:834
    """
    del nzm, ngrdcol
    # Calculate coefs of eddy diffusivity and inverse of dt
    invrs_dt = one / dt

    # Add terms to lhs
    lhs = jnp.zeros_like(lhs_diff)
    lhs = lhs.at[0].set(one_half * lhs_diff[0])
    lhs = lhs.at[1].set(one_half * lhs_diff[1] + invrs_dt)
    lhs = lhs.at[2].set(one_half * lhs_diff[2])

    # LHS mean advection term.
    if not l_implemented:
        if gr.grid_dir_indx > 0:
            lhs = lhs.at[:, :, :nzt - 1].add(lhs_ma_zt[:, :, :nzt - 1])
        else:
            lhs = lhs.at[:, :, 1:nzt].add(lhs_ma_zt[:, :, 1:nzt])

    if l_imp_sfc_momentum_flux:
        # LHS momentum surface flux.
        lhs = lhs.at[1, :, gr.k_lb_zt].add(
            gr.grid_dir
            * invrs_rho_ds_zt[:, gr.k_lb_zt]
            * gr.invrs_dzt[:, gr.k_lb_zt]
            * rho_ds_zm[:, gr.k_lb_zm]
            * (u_star_sqd / wind_speed[:, gr.k_lb_zt])
        )

    return lhs


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "solve_type",
        "l_imp_sfc_momentum_flux",
    ),
)
def windm_edsclrm_rhs(
    nzm, nzt, ngrdcol, gr, solve_type, dt,
    lhs_diff, xm, xm_tndcy,
    rho_ds_zm, invrs_rho_ds_zt,
    l_imp_sfc_momentum_flux, xpwp_sfc,
    stats,
):
    """Calculate the explicit portion of the horizontal wind or eddy-scalar
    time-tendency equation.  See the description in subroutine
    windm_edsclrm_solve for more details.

    References:
      None

      --- THIS SUBROUTINE HAS BEEN OPTIMIZED ---
      Simple changes to this procedure may adversely affect computational speed
          - Gunther Huebler, Aug. 2018, clubb:ticket:834
    """
    del nzm
    # Precalculate 1.0/dt to avoid redoing the divide
    invrs_dt = one / dt

    if solve_type == windm_edsclrm_um:
        name_ta = "um_ta"
    elif solve_type == windm_edsclrm_vm:
        name_ta = "vm_ta"
    else:
        name_ta = ""

    rhs = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    # Lower (ascending grid) or upper (descending grid) boundary calculation
    rhs = rhs.at[:, 0].set(
        one_half
        * (-lhs_diff[1, :, 0] * xm[:, 0] - lhs_diff[0, :, 0] * xm[:, 1])
        + xm_tndcy[:, 0]
        + invrs_dt * xm[:, 0]
    )

    # Non-boundary rhs calculation, this is a highly vectorized loop
    # Note: Some gymnastics are done here to produce a bit-for-bit match
    #       in the ascending vs. descending grid test at -O0 optimization.
    #       This code forces the same calculations to be done in the same
    #       order, regardless of grid direction.
    # For ascending:
    # ( - lhs_diff(3,i,k) * xm(i,k-1)
    #   - lhs_diff(2,i,k) * xm(i,k)
    #   - lhs_diff(1,i,k) * xm(i,k+1) )
    # For descending:
    # ( - lhs_diff(1,i,k) * xm(i,k+1)
    #   - lhs_diff(2,i,k) * xm(i,k)
    #   - lhs_diff(3,i,k) * xm(i,k-1) )
    if gr.grid_dir_indx > 0:
        rhs = rhs.at[:, 1:-1].set(
            one_half
            * (
                -lhs_diff[2, :, 1:-1] * xm[:, :-2]
                -lhs_diff[1, :, 1:-1] * xm[:, 1:-1]
                -lhs_diff[0, :, 1:-1] * xm[:, 2:]
            )
            + xm_tndcy[:, 1:-1]
            + invrs_dt * xm[:, 1:-1]
        )
    else:
        rhs = rhs.at[:, 1:-1].set(
            one_half
            * (
                -lhs_diff[0, :, 1:-1] * xm[:, 2:]
                -lhs_diff[1, :, 1:-1] * xm[:, 1:-1]
                -lhs_diff[2, :, 1:-1] * xm[:, :-2]
            )
            + xm_tndcy[:, 1:-1]
            + invrs_dt * xm[:, 1:-1]
        )

    # Upper (ascending grid) or lower (descending grid) boundary calculation
    rhs = rhs.at[:, -1].set(
        one_half
        * (-lhs_diff[2, :, -1] * xm[:, -2] - lhs_diff[1, :, -1] * xm[:, -1])
        + xm_tndcy[:, -1]
        + invrs_dt * xm[:, -1]
    )

    if stats.l_sample and stats.var_on_stats_list(name_ta):
        stats_tmp = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
        # Lower boundary
        stats_tmp = stats_tmp.at[:, 0].set(
            one_half
            * (lhs_diff[1, :, 0] * xm[:, 0] + lhs_diff[0, :, 0] * xm[:, 1])
        )

        if gr.grid_dir_indx > 0:
            stats_tmp = stats_tmp.at[:, 1:-1].set(
                one_half
                * (
                    lhs_diff[2, :, 1:-1] * xm[:, :-2]
                    + lhs_diff[1, :, 1:-1] * xm[:, 1:-1]
                    + lhs_diff[0, :, 1:-1] * xm[:, 2:]
                )
            )
        else:
            stats_tmp = stats_tmp.at[:, 1:-1].set(
                one_half
                * (
                    lhs_diff[0, :, 1:-1] * xm[:, 2:]
                    + lhs_diff[1, :, 1:-1] * xm[:, 1:-1]
                    + lhs_diff[2, :, 1:-1] * xm[:, :-2]
                )
            )

        # Upper boundary
        stats_tmp = stats_tmp.at[:, -1].set(
            one_half
            * (lhs_diff[2, :, -1] * xm[:, -2] + lhs_diff[1, :, -1] * xm[:, -1])
        )
        # Statistics:  explicit contributions for um or vm.
        # xm term ta has both implicit and explicit components; call
        # stats_begin_budget for the explicit part of turbulent advection.
        stats = stats.begin_budget(name_ta, stats_tmp)

    if not l_imp_sfc_momentum_flux:
        # RHS generalized surface flux.
        sfc_term = (
            gr.grid_dir
            * invrs_rho_ds_zt[:, gr.k_lb_zt]
            * gr.invrs_dzt[:, gr.k_lb_zt]
            * rho_ds_zm[:, gr.k_lb_zm]
            * xpwp_sfc
        )
        rhs = rhs.at[:, gr.k_lb_zt].add(sfc_term)

        if stats.l_sample and stats.var_on_stats_list(name_ta):
            stats_tmp = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
            stats_tmp = stats_tmp.at[:, gr.k_lb_zt].set(sfc_term)
            stats = stats.update_budget(name_ta, stats_tmp)

    return rhs, stats


__all__ = [
    "advance_windm_edsclrm",
    "windm_edsclrm_solve",
    "windm_edsclrm_implicit_stats",
    "compute_uv_tndcy",
    "windm_edsclrm_lhs",
    "windm_edsclrm_rhs",
]
