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
- Diagnostic-only early returns remain reduced; fatal conditions mark err_info
  and the host driver owns stopping after the compiled call.
- Non-zero edsclrm_forcing or surface flux for the thlm/rtm copies marks the
  affected columns fatal but does not return early as the Fortran does.
- Stats are threaded explicitly with JaxStats because the Fortran routine uses
  global stats side effects that are not JAX-compatible state.
- Fortran intent(out)/intent(inout) arguments are returned explicitly.
"""

from functools import partial

import jax
import sys

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
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.model_flags import l_force_descending_solves
from clubb_jax.src.CLUBB_core.parameter_indices import ic_K10, ic_K10h
from clubb_jax.src.CLUBB_core.diffusion import diffusion_zt_lhs
from clubb_jax.src.CLUBB_core.fill_holes import fill_holes_vertical
from clubb_jax.src.CLUBB_core.grid_class import zm2zt
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.matrix_solver_wrapper import tridiag_solve
from clubb_jax.src.CLUBB_core.mean_adv import term_ma_zt_lhs
from clubb_jax.src.CLUBB_core import ErrInfo, Grid, NuVertResDep
from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx


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
    nzm: int, nzt: int, ngrdcol: int, edsclr_dim: int, gr: Grid, sclr_idx: SclrIdx, dt,
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
    With l_do_expldiff_rtm_thlm, sclr_idx must identify distinct, valid
    thlm/rtm slots, which check_clubb_settings verifies. Those slots hold
    temporary copies for diffusion, so their edsclrm_forcing and surface
    flux (wpedsclrp at gr.k_lb_zm, i.e. wpedsclrp_sfc) must be zero.

    References:
    Eqn. 8 & 9 on p. 3545 of
    ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    Method and Model Description'' Golaz, et al. (2002)
    JAS, Vol. 59, pp. 3540--3551.

    Returns:
        (upwp_cl_num, vpwp_cl_num, um, vm, thlm, rtm, edsclrm, upwp, vpwp,
          wpedsclrp, um_pert, vm_pert, upwp_pert, vpwp_pert, err_info, stats)

    Fortran declarations (advance_windm_edsclrm in src/CLUBB_core/advance_windm_edsclrm_module.F90):

        ! ------------------------ Input Variables ------------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol, &
          edsclr_dim

        type (grid), intent(in) :: &
          gr

        type (sclr_idx_type), intent(in) :: &
          sclr_idx  ! Indices identifying the eddy scalars for thlm and rtm [-]

        real( kind = core_rknd ), intent(in) ::  &
          dt                 ! Model timestep                             [s]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
          wm_zt,           & ! w wind component on thermodynamic levels    [m/s]
          ug,              & ! u (west-to-east) geostrophic wind comp.     [m/s]
          vg,              & ! v (south-to-north) geostrophic wind comp.   [m/s]
          um_ref,          & ! Reference u wind component for nudging      [m/s]
          vm_ref,          & ! Reference v wind component for nudging      [m/s]
          um_forcing,      & ! u forcing                                   [m/s/s]
          vm_forcing,      & ! v forcing                                   [m/s/s]
          p_in_Pa,         & ! Air pressure (thermodynamic levels)          [Pa]
          rho_ds_zt,       & ! Dry, static density on thermodynamic levels  [kg/m^3]
          invrs_rho_ds_zt    ! Inv. dry, static density at thermo. levels  [m^3/kg]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) ::  &
          Kh_zm,           & ! Eddy diffusivity for scalars on momentum levs. [m^2/s]
          wp2,             & ! w'^2 (momentum levels)                      [m^2/s^2]
          up2,             & ! u'^2 (momentum levels)                      [m^2/s^2]
          vp2,             & ! v'^2 (momentum levels)                      [m^2/s^2]
          rho_ds_zm          ! Dry, static density on momentum levels      [kg/m^3]

        real( kind = core_rknd ), dimension(ngrdcol,nparams), intent(in) ::  &
          clubb_params     ! Array of CLUBB's tunable parameters            [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,edsclr_dim), intent(in) ::  &
          edsclrm_forcing  ! Eddy scalar large-scale forcing        [{units vary}/s]

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
          fcor           ! Coriolis parameter                            [s^-1]

        logical, intent(in) ::  &
          l_implemented  ! True if CLUBB is being implemented and run in a host model

        type(nu_vertical_res_dep), intent(in) :: &
          nu_vert_res_dep    ! Vertical resolution dependent nu values

        real( kind = core_rknd ), intent(in) ::  &
          ts_nudge               ! Timescale of u/v nudging             [s]

        integer, intent(in) :: &
          tridiag_solve_method  ! Specifier for method to solve tridiagonal systems

        logical, intent(in) :: &
          l_predict_upwp_vpwp,   & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v> alongside
                                   ! the advancement of <rt>, <w'rt'>, <thl>, <wpthlp>, <sclr>, and
                                   ! <w'sclr'> in subroutine advance_xm_wpxp.  Otherwise, <u'w'> and
                                   ! <v'w'> are still approximated by eddy diffusivity when <u> and <v>
                                   ! are advanced in subroutine advance_windm_edsclrm.
          l_upwind_xm_ma,        & ! This flag determines whether we want to use an upwind differencing
                                   ! approximation rather than a centered differencing for turbulent or
                                   ! mean advection terms. It affects rtm, thlm, sclrm, um and vm.
          l_uv_nudge,            & ! For wind speed nudging
          l_lmm_stepping,        & ! Apply Linear Multistep Method (LMM) Stepping
          l_linearize_pbl_winds,  & ! Flag (used by E3SM) to linearize PBL winds
          l_do_expldiff_rtm_thlm   ! Explicitly diffuse rtm and thlm with eddy scalars

        integer, intent(in) :: &
          fill_holes_type

        integer, intent(inout) :: &
          upwp_cl_num, &
          vpwp_cl_num

        ! ------------------------ Input/Output Variables ------------------------
        type(stats_type), intent(inout) :: &
          stats

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) ::  &
          um,   & ! Mean u (west-to-east) wind component         [m/s]
          vm,   & ! Mean v (south-to-north) wind component       [m/s]
          thlm, & ! Liquid potential temperature                 [K]
          rtm     ! Total water mixing ratio                     [kg/kg]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) ::  &
          upwp, & ! <u'w'> (momentum levels)                     [m^2/s^2]
          vpwp    ! <v'w'> (momentum levels)                     [m^2/s^2]

        ! Input/Output Variable for eddy-scalars
        real( kind = core_rknd ), dimension(ngrdcol,nzt,edsclr_dim), intent(inout) ::  &
          edsclrm        ! Mean eddy scalar quantity             [units vary]

        ! Output Variable for eddy-scalars
        real( kind = core_rknd ), dimension(ngrdcol,nzm,edsclr_dim), intent(inout) ::  &
          wpedsclrp      ! w'edsclr' (momentum levels)           [m/s {units vary}]

        ! Variables used to track perturbed version of winds.
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(inout) :: &
          um_pert,   & ! perturbed <u>    [m/s]
          vm_pert      ! perturbed <v>    [m/s]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(inout) :: &
          upwp_pert, & ! perturbed <u'w'> [m^2/s^2]
          vpwp_pert    ! perturbed <v'w'> [m^2/s^2]

        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header
    """
    # Retain immutable pre-solve values solely for the Fortran failure dump.
    um_before, vm_before = um, vm
    if edsclr_dim > 1 and l_do_expldiff_rtm_thlm:

        # The thlm and rtm eddy scalars hold copies of thlm and rtm, which
        # advance_xm_wpxp already advances with thlm_forcing and rtm_forcing.
        # Forcing the copies here would apply those forcings a second time.
        forced_copy = jnp.any(
            (edsclrm_forcing[:, :, sclr_idx.iiedsclr_thl - 1] != zero)
            | (edsclrm_forcing[:, :, sclr_idx.iiedsclr_rt - 1] != zero),
            axis=1,
        )
        err_info = err_info.set_fatal(forced_copy)

        def print_forced_copy(info):
            """Write the source forcing error and global header on the host."""
            print(info.global_error_header_host(), file=sys.stderr)
            print(
                "advance_windm_edsclrm: l_do_expldiff_rtm_thlm requires zero"
                " edsclrm_forcing for iiedsclr_thl and iiedsclr_rt, because"
                " thlm and rtm are already forced in advance_xm_wpxp",
                file=sys.stderr,
            )

        jax.lax.cond(
            jnp.any(forced_copy),
            lambda _: jax.debug.callback(print_forced_copy, err_info, ordered=True),
            lambda _: None,
            operand=None,
        )

        # advance_xm_wpxp also already applies wpthlp_sfc and wprtp_sfc to thlm
        # and rtm. A surface flux on the copies would apply it a second time.
        sfc_flux_on_copy = (
            (wpedsclrp[:, gr.k_lb_zm, sclr_idx.iiedsclr_thl - 1] != zero)
            | (wpedsclrp[:, gr.k_lb_zm, sclr_idx.iiedsclr_rt - 1] != zero)
        )
        err_info = err_info.set_fatal(sfc_flux_on_copy)

        def print_sfc_flux_on_copy(info):
            """Write the source surface-flux error and global header on the host."""
            print(info.global_error_header_host(), file=sys.stderr)
            print(
                "advance_windm_edsclrm: l_do_expldiff_rtm_thlm requires zero"
                " wpedsclrp_sfc for iiedsclr_thl and iiedsclr_rt, because"
                " wpthlp_sfc and wprtp_sfc are already applied in advance_xm_wpxp",
                file=sys.stderr,
            )

        jax.lax.cond(
            jnp.any(sfc_flux_on_copy),
            lambda _: jax.debug.callback(print_sfc_flux_on_copy, err_info, ordered=True),
            lambda _: None,
            operand=None,
        )
    nu_zero = jnp.zeros((ngrdcol,))

    # Coefficient for momentum
    Km_zm = Kh_zm * clubb_params[:, ic_K10][:, None]
    # Coefficient for thermo
    Kmh_zm = Kh_zm * clubb_params[:, ic_K10h][:, None]
    Km_zm_p_nu10 = Km_zm + nu_vert_res_dep.nu10[:, None]

    if edsclr_dim > 1 and l_do_expldiff_rtm_thlm:
        edsclrm = edsclrm.at[:, :, sclr_idx.iiedsclr_thl - 1].set(thlm)
        edsclrm = edsclrm.at[:, :, sclr_idx.iiedsclr_rt - 1].set(rtm)

    edsclrm_before = edsclrm
    l_perturbed_wind = (not l_predict_upwp_vpwp) and l_linearize_pbl_winds

    if not l_implemented:
        lhs_ma_zt = term_ma_zt_lhs(
            nzm, nzt, ngrdcol, wm_zt, gr.weights_zt2zm,
            gr.invrs_dzt, gr.invrs_dzm, l_upwind_xm_ma, gr.grid_dir,
        )
    else:
        lhs_ma_zt = jnp.zeros((ndiags3, ngrdcol, nzt))

    lhs = jnp.zeros((ndiags3, ngrdcol, nzt))
    rhs = jnp.zeros((ngrdcol, nzt, max(2, edsclr_dim)))
    solution = jnp.zeros_like(rhs)
    upwp_chnge = jnp.zeros((ngrdcol, nzm))
    vpwp_chnge = jnp.zeros((ngrdcol, nzm))
    wind_speed = jnp.ones((ngrdcol, nzt))
    u_star_sqd = jnp.zeros((ngrdcol,))

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
            stats.var_on_stats_list("windm_matrix_condt_num")
        )
        solution_wind, err_info, stats = windm_edsclrm_solve(
            nzt, ngrdcol, gr, nrhs,
            tridiag_solve_method,
            l_implemented,
            stats, l_need_rcond,
            lhs, rhs[:, :, :nrhs], err_info,
        )
        if clubb_at_least_debug_level(0):
            jax.lax.cond(
                err_info.any_fatal(),
                lambda _: jax.debug.callback(
                    lambda info: print(
                        info.global_error_header_host(),
                        'Fatal error solving for um/vm',
                        sep="\n",
                        file=sys.stderr,
                    ),
                    err_info, ordered=True,
                ),
                lambda _: None, operand=None,
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

        def accumulate_stats(stats):
            """Collect this diagnostic block only on sampled timesteps."""
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

            return stats

        stats = stats.sampled(accumulate_stats)

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
            def accumulate_stats(stats):
                """Collect this diagnostic block only on sampled timesteps."""
                stats = stats.begin_budget("um_ndg", um / dt)
                stats = stats.begin_budget("vm_ndg", vm / dt)

                return stats

            stats = stats.sampled(accumulate_stats)

            um = um - ((um - um_ref) * (dt / ts_nudge))
            vm = vm - ((vm - vm_ref) * (dt / ts_nudge))

            def accumulate_stats(stats):
                """Collect this diagnostic block only on sampled timesteps."""
                stats = stats.finalize_budget("um_ndg", um / dt)
                stats = stats.finalize_budget("vm_ndg", vm / dt)

                return stats

            stats = stats.sampled(accumulate_stats)

        def accumulate_stats(stats):
            """Collect this diagnostic block only on sampled timesteps."""
            stats = stats.update("um_ref", um_ref)
            stats = stats.update("vm_ref", vm_ref)

            return stats

        stats = stats.sampled(accumulate_stats)

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
        if l_predict_upwp_vpwp:
            def accumulate_stats(stats):
                """Collect this diagnostic block only on sampled timesteps."""
                if upwp_cl_num == 0:
                    stats = stats.begin_budget("upwp_cl", upwp / dt)
                else:
                    stats = stats.update_budget("upwp_cl", -upwp / dt)

                return stats

            stats = stats.sampled(accumulate_stats)
        upwp_cl_num = upwp_cl_num + 1
        upwp, upwp_chnge = clip_covar(
            nzm, ngrdcol, clip_upwp, wp2, up2, upwp,
        )
        if l_predict_upwp_vpwp:
            def accumulate_stats(stats):
                """Collect this diagnostic block only on sampled timesteps."""
                if upwp_cl_num == upwp_cl_max:
                    stats = stats.finalize_budget("upwp_cl", upwp / dt)
                else:
                    stats = stats.update_budget("upwp_cl", upwp / dt)

                return stats

            stats = stats.sampled(accumulate_stats)

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
        if l_predict_upwp_vpwp:
            def accumulate_stats(stats):
                """Collect this diagnostic block only on sampled timesteps."""
                if vpwp_cl_num == 0:
                    stats = stats.begin_budget("vpwp_cl", vpwp / dt)
                else:
                    stats = stats.update_budget("vpwp_cl", -vpwp / dt)

                return stats

            stats = stats.sampled(accumulate_stats)
        vpwp_cl_num = vpwp_cl_num + 1
        vpwp, vpwp_chnge = clip_covar(
            nzm, ngrdcol, clip_vpwp, wp2, vp2, vpwp,
        )
        if l_predict_upwp_vpwp:
            def accumulate_stats(stats):
                """Collect this diagnostic block only on sampled timesteps."""
                if vpwp_cl_num == vpwp_cl_max:
                    stats = stats.finalize_budget("vpwp_cl", vpwp / dt)
                else:
                    stats = stats.update_budget("vpwp_cl", vpwp / dt)

                return stats

            stats = stats.sampled(accumulate_stats)

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
            stats.var_on_stats_list("windm_matrix_condt_num")
        )
        solution_wind, err_info, stats = windm_edsclrm_solve(
            nzt, ngrdcol, gr, nrhs,
            tridiag_solve_method,
            l_implemented,
            stats, l_need_rcond,
            lhs, rhs[:, :, :nrhs], err_info,
        )
        if clubb_at_least_debug_level(0):
            jax.lax.cond(
                err_info.any_fatal(),
                lambda _: jax.debug.callback(
                    lambda info: print(
                        info.global_error_header_host(),
                        'Fatal error solving for um_pert/vm_pert',
                        sep="\n",
                        file=sys.stderr,
                    ),
                    err_info, ordered=True,
                ),
                lambda _: None, operand=None,
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
        if clubb_at_least_debug_level(0):
            jax.lax.cond(
                err_info.any_fatal(),
                lambda _: jax.debug.callback(
                    lambda info: print(
                        info.global_error_header_host(),
                        'Fatal error solving for edsclrm',
                        sep="\n",
                        file=sys.stderr,
                    ),
                    err_info, ordered=True,
                ),
                lambda _: None, operand=None,
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

    if clubb_at_least_debug_level(0):
        # Host I/O adaptation: source field order, complete arrays, one callback.
        def print_wind_error(inputs, outputs, old, scalar_forcing):
            import numpy as np
            with np.printoptions(threshold=np.inf):
                print("Error in advance_windm_edsclrm", file=sys.stderr)
                print("intent(in)", file=sys.stderr)
                print("dt = ", inputs["dt"], file=sys.stderr)
                print("wm_zt = ", inputs["wm_zt"], file=sys.stderr)
                print("Km_zm = ", inputs["Km_zm"], file=sys.stderr)
                print("ug = ", inputs["ug"], file=sys.stderr)
                print("vg = ", inputs["vg"], file=sys.stderr)
                print("um_ref = ", inputs["um_ref"], file=sys.stderr)
                print("vm_ref = ", inputs["vm_ref"], file=sys.stderr)
                print("wp2 = ", inputs["wp2"], file=sys.stderr)
                print("up2 = ", inputs["up2"], file=sys.stderr)
                print("vp2 = ", inputs["vp2"], file=sys.stderr)
                print("um_forcing = ", inputs["um_forcing"], file=sys.stderr)
                print("vm_forcing = ", inputs["vm_forcing"], file=sys.stderr)
                for edsclr in range(edsclr_dim):
                    print("edsclrm_forcing # = ", edsclr + 1, scalar_forcing, file=sys.stderr)
                print("fcor = ", inputs["fcor"], file=sys.stderr)
                print("l_implemented = ", l_implemented, file=sys.stderr)
                print("intent(inout)", file=sys.stderr)
                for name in ("um", "vm"):
                    if l_lmm_stepping:
                        print(name + " (pre-solve) = ", old[name], file=sys.stderr)
                    print(name + " = ", outputs[name], file=sys.stderr)
                for edsclr in range(edsclr_dim):
                    if l_lmm_stepping:
                        print(
                            "edsclrm (pre-solve) # ",
                            edsclr + 1,
                            "=",
                            old["edsclrm"][:, :, edsclr],
                            file=sys.stderr,
                        )
                    print("edsclrm # ", edsclr + 1, "=", outputs["edsclrm"][:, :, edsclr], file=sys.stderr)
                for name in ("upwp", "vpwp", "wpedsclrp"):
                    print(name + " = ", outputs[name], file=sys.stderr)

        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(
                print_wind_error,
                {"dt": dt, "wm_zt": wm_zt, "Km_zm": Km_zm, "ug": ug, "vg": vg, "um_ref": um_ref, "vm_ref": vm_ref, "wp2": wp2, "up2": up2, "vp2": vp2, "um_forcing": um_forcing, "vm_forcing": vm_forcing, "fcor": fcor},
                {"um": um, "vm": vm, "upwp": upwp, "vpwp": vpwp, "wpedsclrp": wpedsclrp, "edsclrm": edsclrm},
                {"um": um_before, "vm": vm_before, "edsclrm": edsclrm_before},
                edsclrm_forcing, ordered=True,
            ),
            lambda _: None, operand=None,
        )

    if edsclr_dim > 1 and l_do_expldiff_rtm_thlm:
        thlm700 = pvertinterp(nzt, ngrdcol, gr, p_in_Pa, 70000.0, thlm)
        thlm1000 = pvertinterp(nzt, ngrdcol, gr, p_in_Pa, 100000.0, thlm)
        apply_explicit_diffusion = (thlm700 - thlm1000) < 20.0

        def explicit_diffusion_stats(previous):
            """Compute the explicit-diffusion tendencies on sampled timesteps."""
            thlm_ed = jnp.where(
                apply_explicit_diffusion[:, None],
                (edsclrm[:, :, sclr_idx.iiedsclr_thl - 1] - thlm) / dt,
                zero,
            )
            rtm_ed = jnp.where(
                apply_explicit_diffusion[:, None],
                (edsclrm[:, :, sclr_idx.iiedsclr_rt - 1] - rtm) / dt,
                zero,
            )

            return thlm_ed, rtm_ed

        thlm_ed, rtm_ed = stats.sampled_value(
            explicit_diffusion_stats,
            (jnp.zeros((ngrdcol, nzt)), jnp.zeros((ngrdcol, nzt))),
        )

        thlm = jnp.where(
            apply_explicit_diffusion[:, None],
            edsclrm[:, :, sclr_idx.iiedsclr_thl - 1],
            thlm,
        )
        rtm = jnp.where(
            apply_explicit_diffusion[:, None],
            edsclrm[:, :, sclr_idx.iiedsclr_rt - 1],
            rtm,
        )
    else:
        thlm_ed = jnp.zeros((ngrdcol, nzt))
        rtm_ed = jnp.zeros((ngrdcol, nzt))

    def accumulate_stats(stats):
        """Collect this diagnostic block only on sampled timesteps."""
        stats = stats.update("thlm_ed", thlm_ed)
        stats = stats.update("rtm_ed", rtm_ed)

        return stats

    stats = stats.sampled(accumulate_stats)

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
    r"""Solves the horizontal wind or eddy-scalar time-tendency equation, and
    diagnoses the turbulent flux.  A Crank-Nicholson time-stepping algorithm
    is used in solving the turbulent advection term and in diagnosing the
    turbulent flux.

    The rate of change of an eddy-scalar quantity, xm, is:

    d(xm)/dt = - w * d(xm)/dz - (1/rho_ds) * d( rho_ds * x'w' )/dz
               + xm_forcings.

    The Turbulent Advection Term
    ----------------------------

    The above equation contains a turbulent advection term:

    - (1/rho_ds) * d( rho_ds * x'w' )/dz;

    where the momentum flux, x'w', is closed using a down gradient approach:

    x'w' = - K_zm * d(xm)/dz.

    The turbulent advection term becomes:

    + (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz;

    which is the same as a standard eddy-diffusion term (if "rho_ds * K_zm" in
    the term above is substituted for "K_zm" in a standard eddy-diffusion
    term, and if the standard eddy-diffusion term is multiplied by
    "1/rho_ds").  Thus, the turbulent advection term is treated and solved in
    the same way that a standard eddy-diffusion term would be solved.  The
    term is discretized as follows:

    The values of xm are found on the thermodynamic levels, while the values
    of K_zm are found on the momentum levels.  Additionally, the values of
    rho_ds_zm are found on the momentum levels, and the values of
    invrs_rho_ds_zt are found on the thermodynamic levels.  The
    derivatives (d/dz) of xm are taken over the intermediate momentum levels.
    At the intermediate momentum levels, d(xm)/dz is multiplied by K_zm and by
    rho_ds_zm.  Then, the derivative of the whole mathematical expression is
    taken over the central thermodynamic level, where it is multiplied by
    invrs_rho_ds_zt, which yields the desired result.

    ---xm(kp1)----------------------------------------------------- t(k+1)

    ===========d(xm)/dz===K_zm(k+1)===rho_ds_zm(k+1)=============== m(k+1)

    ---xm(k)---invrs_rho_ds_zt---d[rho_ds_zm*K_zm*d(xm)/dz]/dz----- t(k)

    ===========d(xm)/dz===K_zm(k)=====rho_ds_zm(k)================= m(k)

    ---xm(km1)----------------------------------------------------- t(k-1)

    The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
    with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
    The letter "t" is used for thermodynamic levels and the letter "m" is used
    for momentum levels.

    invrs_dzt(k)   = 1 / ( zm(k+1) - zm(k) )
    invrs_dzm(k+1) = 1 / ( zt(k+1) - zt(k) )
    invrs_dzm(k)   = 1 / ( zt(k) - zt(k-1) )

    The vertically discretized form of the turbulent advection term (treated
    as an eddy diffusion term) is written out as:

    + invrs_rho_ds_zt(k)
      * invrs_dzt(k)
        * [   rho_ds_zm(k+1) * K_zm(k+1)
              * invrs_dzm(k+1) * ( xm(k+1) - xm(k) )
            - rho_ds_zm(k) * K_zm(k)
              * invrs_dzm(k) * ( xm(k) - xm(k-1) ) ].

    For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme is
    used to solve the (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz
    eddy-diffusion term.  The discretized implicit portion of the term is
    written out as:

    + (1/2) * invrs_rho_ds_zt(k)
      * invrs_dzt(k)
        * [   rho_ds_zm(k+1) * K_zm(k+1)
              * invrs_dzm(k+1) * ( xm(k+1,<t+1>) - xm(k,<t+1>) )
            - rho_ds_zm(k) * K_zm(k)
              * invrs_dzm(k) * ( xm(k,<t+1>) - xm(k-1,<t+1>) ) ].

    Note:  When the implicit term is brought over to the left-hand side,
           the sign is reversed and the leading "+" in front of the term
           is changed to a "-".

    The discretized explicit portion of the term is written out as:

    + (1/2) * invrs_rho_ds_zt(k)
      * invrs_dzt(k)
        * [   rho_ds_zm(k+1) * K_zm(k+1)
              * invrs_dzm(k+1) * ( xm(k+1,<t>) - xm(k,<t>) )
            - rho_ds_zm(k) * K_zm(k)
              * invrs_dzm(k) * ( xm(k,<t>) - xm(k-1,<t>) ) ].

    Timestep index (t) stands for the index of the current timestep, while
    timestep index (t+1) stands for the index of the next timestep, which is
    being advanced to in solving the d(xm)/dt equation.

    Boundary Conditions:

    An eddy-scalar quantity is not allowed to flux out the upper boundary.
    Thus, a zero-flux boundary condition is used for the upper boundary in the
    eddy-diffusion equation.

    The lower boundary condition is much more complicated.  It is neither a
    zero-flux nor a fixed-point boundary condition.  Rather, it is a
    fixed-flux boundary condition.  This term is a turbulent advection term,
    but with the eddy-scalars, the only value of x'w' relevant in solving the
    d(xm)/dt equation is the value of x'w' at the surface (the first momentum
    level), which is written as x'w'|_sfc.

    1) x'w' surface flux; generalized explicit form

       The x'w' surface flux is applied to the d(xm)/dt equation through the
       turbulent advection term, which is:

       - (1/rho_ds) * d( rho_ds * x'w' )/dz.

       At most vertical levels, a substitution can be made for x'w', such
       that:

       x'w' = - K_zm * d(xm)/dz.

       However, the same substitution cannot be made at the surface (momentum
       level 1), as x'w'|_sfc is a surface flux that is explicitly computed
       elsewhere in the model code.

       The lower boundary condition, which in this case needs to be applied to
       the d(xm)/dt equation at level 1, is discretized as follows:

       --xm(2)------------------------------------------------------- t(2)

       ========[x'w'(2) = -K_zm(2)*d(xm)/dz]===rho_ds_zm(2)========== m(2)

       --xm(1)---invrs_rho_ds_zt(1)---d[rho_ds_zm*K_zm*d(xm)/dz]/dz-- t(1)

       ========[x'w'|_sfc]=====================rho_ds_zm(1)========== m(1) sfc

       where "sfc" is the level of the model surface or lower boundary.

       The vertically discretized form of the turbulent advection term
       (treated as an eddy diffusion term), with the explicit surface flux,
       x'w'|_sfc, in place, is written out as:

       - invrs_rho_ds_zt(1)
         * invrs_dzt(1) * [ rho_ds_zm(2) * x'w'(2) - rho_ds_zm(1) * x'w'|_sfc ];

       and after substituting x'w'(2) = -K_zm(2) * invrs_dzm(2) * ( xm(2) - xm(1) ),
       this expression can be re-written as:

       + invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * [   rho_ds_zm(2) * K_zm(2) * invrs_dzm(2) * ( xm(2) - xm(1) )
               + rho_ds_zm(1) * x'w'|_sfc ];

       which can be re-written again as:

       + invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * rho_ds_zm(2) * K_zm(2) * invrs_dzm(2) * ( xm(2) - xm(1) )
       + invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * rho_ds_zm(1) * x'w'|_sfc.

       For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme
       is used to solve the (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz
       eddy-diffusion term.  The discretized implicit portion of the term is
       written out as:

       + (1/2) * invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * [ rho_ds_zm(2) * K_zm(2)
               * invrs_dzm(2) * ( xm(2,<t+1>) - xm(1,<t+1>) ) ].

       Note:  When the implicit term is brought over to the left-hand side,
              the sign is reversed and the leading "+" in front of the term
              is changed to a "-".

       The discretized explicit portion of the term is written out as:

       + (1/2) * invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * [ rho_ds_zm(2) * K_zm(2)
               * invrs_dzm(2) * ( xm(2,<t>) - xm(1,<t>) ) ]
       + invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * rho_ds_zm(1) * x'w'|_sfc.

       Note:  The x'w'|_sfc portion of the term written above has been pulled
              away from the rest of the explicit form written above because
              the (1/2) factor due to Crank-Nicholson time_stepping does not
              apply to it, as there isn't an implicit portion for x'w'|_sfc.

       Timestep index (t) stands for the index of the current timestep, while
       timestep index (t+1) stands for the index of the next timestep, which
       is being advanced to in solving the d(xm)/dt equation.

    2) x'w' surface flux; implicit form for momentum fluxes u'w' and v'w'

       The x'w' surface flux is applied to the d(xm)/dt equation through the
       turbulent advection term, which is:

       - (1/rho_ds) * d( rho_ds * x'w' )/dz.

       At most vertical levels, a substitution can be made for x'w', such
       that:

       x'w' = - K_zm * d(xm)/dz.

       However, the same substitution cannot be made at the surface (momentum
       level 1), as x'w'|_sfc is a surface momentum flux that is found by the
       following equation:

       x'w'|_sfc = - [ u_star^2 / sqrt( um^2 + vm^2 ) ] * xm;

       where x'w'|_sfc and xm are either u'w'|_sfc and um, respectively, or
       v'w'|_sfc and vm, respectively (um and vm are located at the first
       thermodynamic level above the surface, which is thermodynamic level 1),
       sqrt( um^2 + vm^2 ) is the wind speed (also at thermodynamic level 1),
       and u_star is defined as:

       u_star = ( u'w'|_sfc^2 + v'w'|_sfc^2 )^(1/4);

       and thus u_star^2 is defined as:

       u_star^2 = sqrt( u'w'|_sfc^2 + v'w'|_sfc^2 ).

       The value of u_star is either set to a constant value or computed
       (through function diag_ustar) based on the surface wind speed, the
       height above surface of the surface wind speed (as compared to the
       roughness height), and the buoyancy flux at the surface.  Either way,
       u_star is computed elsewhere in the model, and the values of u'w'|_sfc
       and v'w'|_sfc are based on it and computed along with it.  The values
       of u'w'|_sfc and v'w'|_sfc are then passed into advance_clubb_core,
       and are eventually passed into advance_windm_edsclrm.  In subroutine
       advance_windm_edsclrm, the value of u_star_sqd is then recomputed
       based on u'w'|_sfc and v'w'|_sfc.  The value of sqrt( u_star_sqd ) is
       consistent with the value of the original computation of u_star.

       The equation listed above is substituted for x'w'|_sfc.  The lower
       boundary condition, which in this case needs to be applied to the
       d(xm)/dt equation at level 1, is discretized as follows:

       --xm(2)------------------------------------------------------- t(2)

       ===[x'w'(2) = -K_zm(2)*d(xm)/dz]=================rho_ds_zm(2)= m(2)

       --xm(1)---invrs_rho_ds_zt(1)---d[rho_ds_zm*K_zm*d(xm)/dz]/dz-- t(1)

       ===[x'w'|_sfc = -[u_star^2/sqrt(um^2+vm^2)]*xm]==rho_ds_zm(1)= m(1) sfc

       where "sfc" is the level of the model surface or lower boundary.

       The vertically discretized form of the turbulent advection term
       (treated as an eddy diffusion term), with the implicit surface momentum
       flux in place, is written out as:

       - invrs_rho_ds_zt(1)
         * invrs_dzt(1) * [ rho_ds_zm(2) * x'w'(2) - rho_ds_zm(1) * x'w'|_sfc ];

       which can be re-written as:

       - invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * [   rho_ds_zm(2)
                 * { - K_zm(2) * invrs_dzm(2) * ( xm(2) - xm(1) ) }
               - rho_ds_zm(1)
                 * { - [ u_star^2 / sqrt( um(1)^2 + vm(1)^2 ) ] * xm(1) } ];

       which can be re-written as:

       + invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * rho_ds_zm(2) * K_zm(2) * invrs_dzm(2) * ( xm(2) - xm(1) )
       - invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * rho_ds_zm(1) * [ u_star^2 / sqrt( um(1)^2 + vm(1)^2 ) ] * xm(1).

       For this equation, a Crank-Nicholson (semi-implicit) diffusion scheme
       is used to solve the (1/rho_ds) * d [ rho_ds * K_zm * d(xm)/dz ] / dz
       eddy-diffusion term.  The discretized implicit portion of the term is
       written out as:

       + (1/2) * invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * [ rho_ds_zm(2) * K_zm(2)
               * invrs_dzm(2) * ( xm(2,<t+1>) - xm(1,<t+1>) ) ]
       - invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * rho_ds_zm(1)
           * [u_star^2/sqrt( um(1,<t>)^2 + vm(1,<t>)^2 )] * xm(1,<t+1>).

       Note:  When the implicit term is brought over to the left-hand side,
              the signs are reversed and the leading "+" in front of the first
              part of the term is changed to a "-", while the leading "-" in
              front of the second part of the term is changed to a "+".

       Note:  The x'w'|_sfc portion of the term written above has been pulled
              away from the rest of the implicit form written above because
              the (1/2) factor due to Crank-Nicholson time_stepping does not
              apply to it.  The x'w'|_sfc portion of the term is treated
              completely implicitly in order to enhance numerical stability.

       The discretized explicit portion of the term is written out as:

       + (1/2) * invrs_rho_ds_zt(1)
         * invrs_dzt(1)
           * [ rho_ds_zm(2) * K_zm(2)
               * invrs_dzm(2) * ( xm(2,<t>) - xm(1,<t>) ) ].

       Timestep index (t) stands for the index of the current timestep, while
       timestep index (t+1) stands for the index of the next timestep, which
       is being advanced to in solving the d(xm)/dt equation.

    The lower boundary condition for the implicit and explicit portions of the
    turbulent advection term, without the x'w'|_sfc portion of the term, can
    easily be invoked by using the zero-flux boundary conditions found in the
    generalized diffusion function (function diffusion_zt_lhs), which is used
    for many other equations in this model.  Either the generalized explicit
    surface flux needs to be added onto the explicit term after the diffusion
    function has been called from subroutine windm_edsclrm_rhs, or the
    implicit momentum surface flux needs to be added onto the implicit term
    after the diffusion function has been called from subroutine
    windm_edsclrm_lhs.

    Conservation Properties:

    When a fixed-flux lower boundary condition is used (combined with a
    zero-flux upper boundary condition), this technique of discretizing the
    turbulent advection term (treated as an eddy-diffusion term) leads to
    conservative differencing.  When the implicit momentum surface flux is
    either zero or not used, the column totals for each column in the
    left-hand side matrix (for the turbulent advection term) should be equal
    to 0.  Otherwise, the column total for the second column will be equal to
    rho_ds_zm(1) * x'w'|_sfc<t+1>.   When the generalized explicit surface
    flux is either zero or not used, the column total for the right-hand side
    vector (for the turbulent advection term) should be equal to 0.
    Otherwise, the column total for the right-hand side vector (for the
    turbulent advection term) will be equal to rho_ds_zm(1) * x'w'|_sfc<t>.
    This ensures that the total amount of quantity xm over the entire vertical
    domain is only changed by the surface flux (neglecting any forcing terms).
    The total amount of change is equal to rho_ds_zm(1) * x'w'|_sfc.

    To see that this conservation law is satisfied by the left-hand side
    matrix, compute the turbulent advection (treated as eddy diffusion) of xm,
    neglecting any implicit momentum surface flux, multiply by rho_ds_zt, and
    integrate vertically.  In discretized matrix notation (where "i" stands
    for the matrix column and "j" stands for the matrix row):

     0 = Sum_j Sum_i
          (rho_ds_zt)_i ( dzt )_i
          ( 0.5 * (1/rho_ds_zt) * (1/dzt) * (rho_ds_zm*K_zm*(1/dzm)) )_ij (xm<t+1>)_j.

    The left-hand side matrix,
    ( 0.5 * (1/rho_ds_zt) * (1/dzt) * (rho_ds_zm*K_zm*(1/dzm)) )_ij, is
    partially written below.  The sum over i in the above equation removes
    (1/rho_ds_zt) and (1/dzt) everywhere from the matrix below.  The sum
    over j leaves the column totals that are desired, which are 0.

    Left-hand side matrix contributions from the turbulent advection term
    (treated as an eddy-diffusion term using a Crank-Nicholson timestep);
    first four vertical levels:

         -------------------------------------------------------------------------------->
    k=1 |  +0.5*                    -0.5*                                0
        |    (1/rho_ds_zt(k))*        (1/rho_ds_zt(k))*
        |    (1/dzt(k))*              (1/dzt(k))*
        |    rho_ds_zm(k+1)*          rho_ds_zm(k+1)*
        |    K_zm(k+1)*(1/dzm(k+1))   K_zm(k+1)*(1/dzm(k+1))
        |
    k=2 |  -0.5*                    +0.5*                      -0.5*
        |    (1/rho_ds_zt(k))*        (1/rho_ds_zt(k))*          (1/rho_ds_zt(k))*
        |    (1/dzt(k))*              (1/dzt(k))*                (1/dzt(k))*
        |    rho_ds_zm(k)*            [ rho_ds_zm(k+1)*          rho_ds_zm(k+1)*
        |    K_zm(k)*(1/dzm(k))         K_zm(k+1)*(1/dzm(k+1))   K_zm(k+1)*(1/dzm(k+1))
        |                              +rho_ds_zm(k)*
        |                               K_zm(k)*(1/dzm(k)) ]
        |
    k=3 |            0              -0.5*                      +0.5*
        |                             (1/rho_ds_zt(k))*          (1/rho_ds_zt(k))*
        |                             (1/dzt(k))*                (1/dzt(k))*
        |                             rho_ds_zm(k)*              [ rho_ds_zm(k+1)*
        |                             K_zm(k)*(1/dzm(k))           K_zm(k+1)*(1/dzm(k+1))
        |                                                         +rho_ds_zm(k)*
        |                                                          K_zm(k)*(1/dzm(k)) ]
        |
    k=4 |            0                        0                -0.5*
        |                                                        (1/rho_ds_zt(k))*
        |                                                        (1/dzt(k))*
        |                                                        rho_ds_zm(k)*
        |                                                        K_zm(k)*(1/dzm(k))
       \ /

    Note:  The superdiagonal term from level 3 and both the main diagonal and
           superdiagonal terms from level 4 are not shown on this diagram.

    Note:  If an implicit momentum surface flux is used, an additional term,
           + (1/rho_ds_zt(1)) * (1/dzt(1)) * rho_ds_zm(1)
             * [ u_star^2 / sqrt( um(1,<t>)^2 + vm(1,<t>)^2 ) ], is added to
           row 1 (k=1), column 1.

    To see that the above conservation law is satisfied by the right-hand side
    vector, compute the turbulent advection (treated as eddy diffusion) of xm,
    neglecting any generalized explicit surface flux, multiply by rho_ds_zt,
    and integrate vertically.  In discretized matrix notation (where "i"
    stands for the matrix column and "j" stands for the matrix row):

     0 = Sum_j Sum_i (rho_ds_zt)_i ( dzt )_i ( rhs_vector )_j.

    The right-hand side vector, ( rhs_vector )_j, is partially written below.
    The sum over i in the above equation removes (1/rho_ds_zt) and (1/dzt)
    everywhere from the vector below.  The sum over j leaves the column total
    that is desired, which is 0.

    Right-hand side vector contributions from the turbulent advection term
    (treated as an eddy-diffusion term using a Crank-Nicholson timestep);
    first four vertical levels:

         --------------------------------------------------
    k=1 | +0.5*(1/rho_ds_zt(k))*                           |
        |      (1/dzt(k))*                                 |
        |       [ rho_ds_zm(k+1)*K_zm(k+1)*                |
        |         (1/dzm(k+1))*(xm(k+1,<t>)-xm(k,<t>)) ]   |
        |                                                  |
    k=2 | +0.5*(1/rho_ds_zt(k))*                           |
        |      (1/dzt(k))*                                 |
        |       [ rho_ds_zm(k+1)*K_zm(k+1)*                |
        |         (1/dzm(k+1))*(xm(k+1,<t>)-xm(k,<t>))     |
        |        -rho_ds_zm(k)*K_zm(k)*                    |
        |         (1/dzm(k))*(xm(k,<t>)-xm(k-1,<t>)) ]     |
        |                                                  |
    k=3 | +0.5*(1/rho_ds_zt(k))*                           |
        |      (1/dzt(k))*                                 |
        |       [ rho_ds_zm(k+1)*K_zm(k+1)*                |
        |         (1/dzm(k+1))*(xm(k+1,<t>)-xm(k,<t>))     |
        |        -rho_ds_zm(k)*K_zm(k)*                    |
        |         (1/dzm(k))*(xm(k,<t>)-xm(k-1,<t>)) ]     |
        |                                                  |
    k=4 | +0.5*(1/rho_ds_zt(k))*                           |
        |      (1/dzt(k))*                                 |
        |       [ rho_ds_zm(k+1)*K_zm(k+1)*                |
        |         (1/dzm(k+1))*(xm(k+1,<t>)-xm(k,<t>))     |
        |        -rho_ds_zm(k)*K_zm(k)*                    |
        |         (1/dzm(k))*(xm(k,<t>)-xm(k-1,<t>)) ]     |
       \ /                                                \ /

    Note:  If a generalized explicit surface flux is used, an additional term,
           + (1/rho_ds_zt(1)) * (1/dzt(1)) * rho_ds_zm(1) * x'w'|_sfc, is
           added to row 1 (k=1).

    Note:  Only the contributions by the turbulent advection term are shown
           for both the left-hand side matrix and the right-hand side vector.
           There are more terms in the equation, and thus more factors to be
           added to both the left-hand side matrix (such as time tendency and
           mean advection) and the right-hand side vector (such as xm
           forcings).  The left-hand side matrix is set-up so that a singular
           matrix is not encountered.

    References:
    Eqn. 8 & 9 on p. 3545 of
    ``A PDF-Based Model for Boundary Layer Clouds. Part I:
      Method and Model Description'' Golaz, et al. (2002)
    JAS, Vol. 59, pp. 3540--3551.

    Returns:
        (solution, err_info, stats)

    Fortran declarations (windm_edsclrm_solve in src/CLUBB_core/advance_windm_edsclrm_module.F90):

        ! ------------------------ Input Variables ------------------------
        integer, intent(in) :: &
          nzt, &
          ngrdcol

        type (grid), intent(in) :: &
          gr

        integer, intent(in) :: &
          nrhs ! Number of right-hand side (explicit) vectors & Number of solution vectors.

        integer, intent(in) :: &
          tridiag_solve_method  ! Specifier for method to solve tridiagonal systems

        logical, intent(in) :: &
          l_implemented  ! True if CLUBB is being implemented and run in a host model

        type(stats_type), intent(inout) :: &
          stats

        logical, intent(in) :: &
          l_need_rcond

        ! ------------------------ Inout variables ------------------------
        real( kind = core_rknd ), dimension(3,ngrdcol,nzt), intent(inout) :: &
          lhs    ! Implicit contributions to um, vm, and eddy scalars  [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,nzt,nrhs), intent(inout) :: &
          rhs    ! Right-hand side (explicit) contributions.

        type(err_info_type), intent(inout) :: &
          err_info      ! err_info struct containing err_code and err_header

        ! ------------------------ Output variables ------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt,nrhs), intent(out) :: &
          solution ! Solution to the system of equations    [units vary]
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
        l_use_rcond=stats.l_sample if l_need_rcond else False,
        l_implemented=l_implemented,
    )

    if clubb_at_least_debug_level(0):
        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(
                lambda info: print(
                    info.global_error_header_host(),
                    'calling tridiag_solve in windm_edsclrm_solve',
                    sep="\n",
                    file=sys.stderr,
                ),
                err_info, ordered=True,
            ),
            lambda _: None, operand=None,
        )

    if l_need_rcond:
        def accumulate_stats(stats):
            """Collect this diagnostic block only on sampled timesteps."""
            stats = stats.update("windm_matrix_condt_num", one / rcond)

            return stats

        stats = stats.sampled(accumulate_stats)

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
    """Compute implicit contributions to um and vm

    Fortran declarations (windm_edsclrm_implicit_stats in src/CLUBB_core/advance_windm_edsclrm_module.F90):

        !---------------------- Input Variables ----------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        integer, intent(in) :: &
          solve_type     ! Desc. of what is being solved for

        type(grid), intent(in) :: &
          gr

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          xm,         & !  Computed value um or vm at <t+1>    [m/s]
          invrs_dzt     ! The inverse spacing between momentum grid levels;
                        ! centered over thermodynamic grid levels.

        real( kind = core_rknd ), dimension(3,ngrdcol,nzt), intent(in) :: &
          lhs_diff, & ! LHS diffustion terms
          lhs_ma_zt   ! LHS mean advection terms

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          wind_speed,      & ! wind speed; sqrt(u^2 + v^2)              [m/s]
          invrs_rho_ds_zt    ! Inv. dry, static density at thermo. levels  [m^3/kg]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          rho_ds_zm          ! Dry, static density on momentum levels      [kg/m^3]

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          u_star_sqd   ! Surface friction velocity, u_star, squared      [m/s]

        logical, intent(in) :: &
          l_imp_sfc_momentum_flux  ! Flag for implicit momentum surface fluxes.

        type(stats_type), intent(inout) :: &
          stats
    """
    del nzm
    def accumulate_stats(stats):
        """Collect diagnostics only on sampled timesteps."""

        if solve_type == windm_edsclrm_um:
            name_ma = "um_ma"
            name_ta = "um_ta"
        elif solve_type == windm_edsclrm_vm:
            name_ma = "vm_ma"
            name_ta = "vm_ta"
        else:
            return stats

        imp_sfc_flux = jnp.zeros((ngrdcol, nzt))

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

        stats_tmp = jnp.zeros((ngrdcol, nzt))
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

        stats_tmp = jnp.zeros((ngrdcol, nzt))
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

    return stats.sampled(accumulate_stats)


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

    Returns:
        (xm_tndcy, stats)

    Fortran declarations (compute_uv_tndcy in src/CLUBB_core/advance_windm_edsclrm_module.F90):

        ! -------------------------- Input Variables --------------------------
        integer, intent(in) :: &
          nzt, &
          ngrdcol

        integer, intent(in) ::  &
          solve_type      ! Description of what is being solved for

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
          fcor            ! Coriolis parameter     [s^-1]

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          perp_wind_m,  & ! Perpendicular component of the mean wind (e.g. v, for the u-eqn) [m/s]
          perp_wind_g,  & ! Perpendicular component of the geostropic wind (e.g. vg)         [m/s]
          xm_forcing      ! Prescribed wind forcing                                          [m/s/s]

        logical, intent(in) :: &
          l_implemented   ! True if CLUBB is being implemented and run in a host model

        type(stats_type), intent(inout) :: &
          stats

        ! -------------------------- Output Variables --------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) ::  &
          xm_tndcy        ! xm tendency            [m/s^2]
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
            xm_gf = jnp.zeros((ngrdcol, nzt))
            xm_cf = jnp.zeros((ngrdcol, nzt))

        xm_tndcy = xm_gf + xm_cf + xm_forcing

        def accumulate_stats(stats):
            """Collect this diagnostic block only on sampled timesteps."""
            # Statistics:  explicit contributions for um or vm.
            # xm term gf is completely explicit; record with stats_update.
            stats = stats.update(name_gf, xm_gf)
            # xm term cf is completely explicit; record with stats_update.
            stats = stats.update(name_cf, xm_cf)
            # xm term F
            stats = stats.update(name_f, xm_forcing)

            return stats

        stats = stats.sampled(accumulate_stats)
    else:
        xm_tndcy = jnp.zeros((ngrdcol, nzt))

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
    time-tendency equation.  See the description in function
    windm_edsclrm_solve for more details.

    Returns:
        lhs

    Fortran declarations (windm_edsclrm_lhs in src/CLUBB_core/advance_windm_edsclrm_module.F90):

        ! ----------------------- Input Variables -----------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        type (grid), intent(in) :: gr

        real( kind = core_rknd ), intent(in) :: &
          dt                 ! Model timestep                             [s]

        real( kind = core_rknd ), intent(in), dimension(3,ngrdcol,nzt) :: &
          lhs_diff, & ! LHS diffustion terms
          lhs_ma_zt   ! LHS mean advection terms

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          wind_speed,      & ! wind speed; sqrt( u^2 + v^2 )              [m/s]
          invrs_rho_ds_zt    ! Inv. dry, static density at thermo. levels [m^3/kg]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          rho_ds_zm          ! Dry, static density on momentum levels     [kg/m^3]

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          u_star_sqd    ! Surface friction velocity, u_*, squared  [m/s]

        logical, intent(in) ::  &
          l_implemented, & ! True if CLUBB is being implemented and run in a host model
          l_imp_sfc_momentum_flux  ! Flag for implicit momentum surface fluxes.

        ! ----------------------- Output Variable -----------------------
        real( kind = core_rknd ), dimension(3,ngrdcol,nzt), intent(out) :: &
          lhs           ! Implicit contributions to xm (tridiagonal matrix)
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
    time-tendency equation.  See the description in function
    windm_edsclrm_solve for more details.

    Returns:
        (rhs, stats)

    Fortran declarations (windm_edsclrm_rhs in src/CLUBB_core/advance_windm_edsclrm_module.F90):

        !------------------- Input Variables -------------------
        integer, intent(in) :: &
          nzm, &
          nzt, &
          ngrdcol

        type (grid), intent(in) :: &
          gr

        integer, intent(in) :: &
          solve_type ! Description of what is being solved for

        real( kind = core_rknd ), intent(in) :: &
          dt                 ! Model timestep                             [s]

        real( kind = core_rknd ), dimension(3,ngrdcol,nzt), intent(in) :: &
          lhs_diff   ! LHS diffusion terms

        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
          xm,              & ! Eddy-scalar variable, xm (thermo. levels)  [units vary]
          xm_tndcy,        & ! The explicit time-tendency acting on xm    [units vary]
          invrs_rho_ds_zt    ! Inv. dry, static density at thermo. levels [m^3/kg]

        real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
          rho_ds_zm          ! Dry, static density on momentum levels     [kg/m^3]

        real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
          xpwp_sfc     ! x'w' at the surface                              [units vary]

        logical, intent(in) :: &
          l_imp_sfc_momentum_flux  ! Flag for implicit momentum surface fluxes.

        !------------------- Inout Variable -------------------
        type(stats_type), intent(inout) :: &
          stats

        !------------------- Output Variable -------------------
        real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
          rhs          ! Right-hand side (explicit) contributions.
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

    rhs = jnp.zeros((ngrdcol, nzt))
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

    if stats.var_on_stats_list(name_ta):
        def accumulate_stats(stats):
            """Collect this diagnostic block only on sampled timesteps."""
            stats_tmp = jnp.zeros((ngrdcol, nzt))
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

            return stats

        stats = stats.sampled(accumulate_stats)

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

        if stats.var_on_stats_list(name_ta):
            def accumulate_stats(stats):
                """Collect this diagnostic block only on sampled timesteps."""
                stats_tmp = jnp.zeros((ngrdcol, nzt))
                stats_tmp = stats_tmp.at[:, gr.k_lb_zt].set(sfc_term)
                stats = stats.update_budget(name_ta, stats_tmp)

                return stats

            stats = stats.sampled(accumulate_stats)

    return rhs, stats
