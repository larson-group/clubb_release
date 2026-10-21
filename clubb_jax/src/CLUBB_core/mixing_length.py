"""Compute CLUBB mixing distances and turbulence dissipation time scales.

The mixing length estimates how far turbulent air motions exchange heat,
water, and momentum. This port follows ``src/CLUBB_core/mixing_length.F90``:
``Lscale_method`` selects the legacy buoyant-parcel method (1) or length
from a dissipation time scale (2).

Usage:
  - ``calc_Lscale`` computes the complete timestep diagnosis and statistics.
  - ``set_Lscale_max`` returns the length cap from host horizontal grid spacing,
    or 100000 m for standalone runs.
  - ``calc_Lscale_directly`` computes parcel method 1, including the
    source-controlled optional perturbed-launch calculation.
  - ``diagnose_Lscale_from_tau`` computes method 2 and its individual inverse
    damping-time contributions.
  - ``compute_mixing_length`` implements the legacy parcel kernel and is
    private in the Fortran module.

Understanding:
  Fields have shape ``(ngrdcol, nzt)`` on thermodynamic levels or
  ``(ngrdcol, nzm)`` on momentum levels. Lengths are in m, time scales in s,
  inverse times in 1/s, temperatures in K, and water mixing ratios in kg/kg.
  Column parameters have no height axis. Callers provide initialized grid,
  PDF (probability density function), statistics, and error structures.
  Updated statistics and per-column errors are returned as immutable state;
  callers must retain them. Parcel histories are local to each call.

JAX adaptations:
  The parcel kernel uses fixed-length scans and masks for source early exits.
  A scan offset counts level moves from launch, not meters. Physical-height
  ordering supports either grid orientation; interpolation retains native
  grid weights. The standalone constructor currently allows ascending grids.
  Launch grouping evaluates several launch levels together. Local callbacks
  implement JAX scan/conditional boundaries; public signatures and routine
  order follow Fortran, with outputs and inout state returned functionally.
  Fatal masks remain per column. Ordered host callbacks reproduce source
  error messages and bad-column array dumps without leaving JAX physics.
  The enclosing driver supplies the calc_Lscale timing region.

References:
  Golaz et al. (2002), JAS 59, 3540--3551 gives parcel-method
  background; Guo et al. (2021, JAMES) describes length diagnosed from tau.
"""

from __future__ import annotations

import jax
import sys

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_helper_module import (
    sqrt_clipped,
    calc_Ri_zm,
    smooth_heaviside_peskin,
    smooth_max,
    smooth_min,
)
from clubb_jax.src.CLUBB_core.constants_clubb import (
    Cp,
    Lv,
    Rd,
    em_min,
    ep,
    ep1,
    ep2,
    eps,
    grav,
    min_max_smth_mag,
    one,
    one_fourth,
    one_half,
    one_third,
    rt_tol,
    thl_tol,
    two,
    unused_var,
    vonk,
    zero,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.parameter_indices import (
    iC_invrs_tau_N2,
    iC_invrs_tau_N2_clear_wp3,
    iC_invrs_tau_N2_wp2,
    iC_invrs_tau_N2_wpxp,
    iC_invrs_tau_N2_xp2,
    iC_invrs_tau_bkgnd,
    iC_invrs_tau_sfc,
    iC_invrs_tau_shear,
    iC_invrs_tau_wpxp_N2_thresh,
    iC_invrs_tau_wpxp_Ri,
    ialtitude_threshold,
    iLscale_mu_entrain,
    iLscale_mu_coef,
    iLscale_pert_coef,
    itaumax,
    iwpxp_Ri_exp,
    iz_displace,
)
from clubb_jax.src.CLUBB_core.model_flags import Lscale_direct, Lscale_tau
from clubb_jax.src.CLUBB_core.error_code import CLUBB_FATAL_ERROR, clubb_at_least_debug_level
from clubb_jax.src.CLUBB_core.grid_class import zm2zt, zm2zt2zm, zt2zm
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq


#===============================================================================

def calc_Lscale(
    nzm, nzt, ngrdcol, gr, l_implemented, host_dx, host_dy,
    p_in_Pa, exner, rtm, thlm, thvm, thlp2, rtp2, rtpthlp,
    pdf_params, em, thv_ds_zt, lmin,
    upwp_sfc, vpwp_sfc, ddzt_umvm_sqd, ice_supersat_frac,
    ufmin, tau_const, sfc_elevation, clubb_params,
    saturation_formula, l_Lscale_plume_centered,
    Lscale_method,
    l_e3sm_config,
    l_smooth_Heaviside_tau_wpxp,
    l_modify_limiters_for_cnvg_test,
    l_use_invrs_tau_N2_iso,
    brunt_vaisala_freq_sqd_smth,
    stats, err_info,
):
    """Compute CLUBB's mixing length and dissipation time scales, including
    Richardson-number-dependent damping terms and Lscale-related statistics.

    Read thermodynamic fields on ``nzt`` levels, moments on ``nzm`` levels, and
    per-column parameters. Grid and PDF state must be initialized by the caller.
    ``Lscale_method`` selects methods 1 and 2. Lengths are
    in m, times in s, and inverse times in 1/s; ``_zt``/``_zm`` name the grids.

    Return ``(err_info, invrs_tau_zt, invrs_tau_zm, invrs_tau_xp2_zm,
    invrs_tau_wp3_zt, invrs_tau_C1_zm, invrs_tau_C4_zm, invrs_tau_C6_zm,
    invrs_tau_C14_zm, tau_max_zm, tau_max_zt, tau_zm, Lscale, Lscale_zm,
    Lscale_up, Lscale_down, stats)``. Retain the returned error/statistics state.

    References:
      Guo et al. (2021, JAMES) for the diagnosed Lscale-from-tau option.
    """
    l_smooth_min_max = False  # Fortran compile-time parameter; retain both branches.

    # Determine the maximum allowable value for Lscale (in meters).
    Lscale_max = set_Lscale_max(ngrdcol, l_implemented, host_dx, host_dy)

    # Interpolate rtpthlp here so the diagnostic stats value is generated at
    # the same point as the Lscale calculation that consumes it.
    rtpthlp_zt = zm2zt(nzm, nzt, ngrdcol, gr, rtpthlp)
    if stats.l_sample:
        stats = stats.update("rtpthlp_zt", rtpthlp_zt)

    # Calculate Richardson number Ri_zm.
    if l_modify_limiters_for_cnvg_test:
        Ri_zm = calc_Ri_zm(
            nzm, ngrdcol, brunt_vaisala_freq_sqd_smth, ddzt_umvm_sqd,
            0.0, 1.0e-12,
        )
        Ri_zm = zm2zt2zm(nzm, nzt, ngrdcol, gr, Ri_zm)
    else:
        if l_smooth_min_max:
            brunt_vaisala_freq_clipped = smooth_max(
                nzm, ngrdcol, 1.0e-7,
                brunt_vaisala_freq_sqd_smth,
                1.0e-4 * min_max_smth_mag,
            )
            ddzt_umvm_sqd_clipped = smooth_max(
                nzm, ngrdcol, ddzt_umvm_sqd,
                1.0e-7,
                1.0e-6 * min_max_smth_mag,
            )
            Ri_zm = calc_Ri_zm(
                nzm, ngrdcol, brunt_vaisala_freq_clipped,
                ddzt_umvm_sqd_clipped, 0.0, 0.0,
            )
        else:
            Ri_zm = calc_Ri_zm(
                nzm, ngrdcol, brunt_vaisala_freq_sqd_smth, ddzt_umvm_sqd,
                1.0e-7, 1.0e-7,
            )

    if stats.l_sample:
        stats = stats.update("Ri_zm", Ri_zm)

    if Lscale_method != Lscale_tau:
        # Compute Lscale first using the buoyant parcel calculation.

        # Interpolate thlp2 and rtp2 to thermodynamic levels.
        thlp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, thlp2,
                       thl_tol**2)  # Positive definite quantity
        rtp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, rtp2,
                      rt_tol**2)   # Positive definite quantity

        err_info, Lscale, Lscale_up, Lscale_down, stats = calc_Lscale_directly(
            ngrdcol, nzm, nzt, gr,  # In
            l_implemented, p_in_Pa,  # In
            exner, rtm, thlm, thvm,  # In
            rtp2_zt, thlp2_zt, rtpthlp_zt, pdf_params,  # In
            em, thv_ds_zt, Lscale_max, lmin,  # In
            clubb_params,  # In
            saturation_formula,  # In
            l_Lscale_plume_centered,  # In
            Lscale_method,  # In
            stats, err_info,  # InOut
        )
        # JAX adaptation: propagate per-column fatal status to the driver
        # instead of returning the whole batch at a host debug/print boundary.

        # Calculate CLUBB's turbulent eddy-turnover time scale as CLUBB's
        # length scale divided by a velocity scale.
        em_zt = zm2zt(nzm, nzt, ngrdcol, gr, em, em_min)
        tau_zt = jnp.minimum(
            Lscale / jnp.sqrt(em_zt),
            clubb_params[:, itaumax][:, None],
        )

        Lscale_zm = zt2zm(nzm, nzt, ngrdcol, gr, Lscale, zero_threshold)
        tau_zm = jnp.minimum(
            Lscale_zm / jnp.sqrt(jnp.maximum(em_min, em)),
            clubb_params[:, itaumax][:, None],
        )

        invrs_tau_zm = one / tau_zm
        invrs_tau_wp2_zm = invrs_tau_zm
        invrs_tau_xp2_zm = invrs_tau_zm
        invrs_tau_wpxp_zm = invrs_tau_zm
        invrs_tau_wp3_zm = invrs_tau_zm
        tau_max_zm = jnp.broadcast_to(clubb_params[:, itaumax][:, None], (ngrdcol, nzm))

        invrs_tau_zt = one / tau_zt
        invrs_tau_wp3_zt = invrs_tau_zt
        tau_max_zt = jnp.broadcast_to(clubb_params[:, itaumax][:, None], (ngrdcol, nzt))

        invrs_tau_N2_iso = jnp.zeros((ngrdcol, nzm))
    else:
        # Diagnose simple tau and Lscale.
        (
            err_info, invrs_tau_zt, invrs_tau_zm,
            invrs_tau_sfc, invrs_tau_no_N2_zm, invrs_tau_bkgnd,
            invrs_tau_shear, invrs_tau_N2_iso,
            invrs_tau_wp2_zm, invrs_tau_xp2_zm,
            invrs_tau_wp3_zm, invrs_tau_wp3_zt, invrs_tau_wpxp_zm,
            tau_max_zm, tau_max_zt, tau_zm, tau_zt,
            Lscale, Lscale_up, Lscale_down, stats,
        ) = diagnose_Lscale_from_tau(
            nzm, nzt, ngrdcol, gr,  # In
            upwp_sfc, vpwp_sfc, ddzt_umvm_sqd,  # In
            ice_supersat_frac,  # In
            em,  # In
            ufmin, tau_const,  # In
            sfc_elevation, Lscale_max,  # In
            clubb_params,  # In
            stats,  # InOut
            l_e3sm_config,  # In
            l_smooth_Heaviside_tau_wpxp,  # In
            brunt_vaisala_freq_sqd_smth, Ri_zm,  # In
            err_info,  # InOut
        )
        # The same per-column error propagation applies to the tau method.

    Lscale_zm = zt2zm(nzm, nzt, ngrdcol, gr, Lscale)

    invrs_tau_C6_zm = invrs_tau_wpxp_zm
    invrs_tau_C1_zm = invrs_tau_wp2_zm
    invrs_tau_C14_zm = invrs_tau_wp2_zm

    if (Lscale_method != Lscale_tau) and l_use_invrs_tau_N2_iso:
        jax.debug.callback(
            lambda: print("calc_Lscale: l_use_invrs_tau_N2_iso requires Lscale_method=2.",
                          file=sys.stderr), ordered=True,
        )
        err_info = err_info.set_fatal()

    if not l_use_invrs_tau_N2_iso:
        invrs_tau_C4_zm = invrs_tau_wp2_zm
    else:
        invrs_tau_C4_zm = invrs_tau_N2_iso

    if stats.l_sample:
        stats = stats.update("Lscale", Lscale)
        stats = stats.update("Lscale_up", Lscale_up)
        stats = stats.update("Lscale_down", Lscale_down)
        stats = stats.update("tau_zm", tau_zm)
        stats = stats.update("invrs_tau_zm", invrs_tau_zm)
        stats = stats.update("invrs_tau_xp2_zm", invrs_tau_xp2_zm)
        stats = stats.update("invrs_tau_wp2_zm", invrs_tau_wp2_zm)
        stats = stats.update("invrs_tau_wpxp_zm", invrs_tau_wpxp_zm)
        stats = stats.update("invrs_tau_wp3_zm", invrs_tau_wp3_zm)
        stats = stats.update("tau_zt", tau_zt)

        if Lscale_method == Lscale_tau:
            stats = stats.update("invrs_tau_no_N2_zm", invrs_tau_no_N2_zm)
            stats = stats.update("invrs_tau_bkgnd", invrs_tau_bkgnd)
            stats = stats.update("invrs_tau_sfc", invrs_tau_sfc)
            stats = stats.update("invrs_tau_shear", invrs_tau_shear)

    return (
        err_info,
        invrs_tau_zt,
        invrs_tau_zm,
        invrs_tau_xp2_zm,
        invrs_tau_wp3_zt,
        invrs_tau_C1_zm,
        invrs_tau_C4_zm,
        invrs_tau_C6_zm,
        invrs_tau_C14_zm,
        tau_max_zm,
        tau_max_zt,
        tau_zm,
        Lscale,
        Lscale_zm,
        Lscale_up,
        Lscale_down,
        stats,
    )


#===============================================================================

def set_Lscale_max(ngrdcol, l_implemented, host_dx, host_dy):
    """Return one maximum mixing length per column, in meters.

    For a host model, the cap is one quarter of the smaller horizontal spacing
    ``host_dx`` or ``host_dy`` [m]. Smaller host cells therefore limit parcel
    lengths and shorten damping times, increasing damping of the turbulent
    moments. The source rationale is that smaller host cells usually contain
    less variation in meteorological quantities than larger cells.
    Standalone calls return 100000 m for each of ``ngrdcol`` columns,
    independently of the supplied spacings, so this cap ordinarily does not bind.
    """
    # Determine the maximum allowable value for Lscale (in meters).
    if l_implemented:
        return 0.25 * jnp.minimum(host_dx, host_dy)
    return jnp.full((ngrdcol,), 1.0e5)


#===============================================================================

def compute_mixing_length(
    nzm, nzt, ngrdcol, gr, thvm, thlm,
    rtm, em, Lscale_max, p_in_Pa,
    exner, thv_ds, Lscale_mu_entrain, lmin,
    saturation_formula,
    l_implemented,
    err_info,
):
    """Compute the legacy moist, nonlocal parcel mixing length.

    Read environmental profiles on thermodynamic levels, momentum-level TKE
    [m^2/s^2], per-column entrainment rates [1/m], and length floors/caps [m].
    Return ``(err_info, Lscale, Lscale_up, Lscale_down)`` with lengths shaped
    ``(ngrdcol, nzt)`` in m. ``saturation_formula`` selects shared saturation.

    The layer-linear environment permits analytic entrainment updates. Initial
    TKE is exhausted by trapezoidal buoyancy work; a linear/quadratic root gives
    the distance within the first impassable layer. Previous launch endpoints
    extend the local lengths before the surface floor and geometric mean.

    Notes:
      The equation for the rate of change of theta_l and r_t of the parcel with
      respect to height, due to entrainment, is:

              d(thl_par)/dz = - Lscale_mu_entrain * ( thl_parcel - thl_environment );

              d(rt_par)/dz = - Lscale_mu_entrain * ( rt_parcel - rt_environment );

      where Lscale_mu_entrain is the entrainment rate,
      such that:

              Lscale_mu_entrain = (1/m)*(dm/dz);

      where m is the mass of the parcel.  The value of Lscale_mu_entrain is set to be a
      constant.

      The differential equations are solved for given the boundary condition
      and given the fact that the value of thl_environment and rt_environment
      are treated as changing linearly for a parcel of air from one grid level
      to the next.

      For the special case where entrainment rate, Lscale_mu_entrain, is set to 0,
      thl_parcel and rt_parcel remain constant


      The equation for Lscale_up is:

          INT(z_i:z_i+Lscale_up) g * ( thv_par - thvm ) / thvm dz = -em(z_i);

      and for Lscale_down

          INT(z_i-Lscale_down:z_i) g * ( thv_par - thvm ) / thvm dz = em(z_i);

      where thv_par is theta_v of the parcel, thvm is the mean
      environmental value of theta_v, z_i is the altitude that the parcel
      started from, and em is the mean value of TKE at
      altitude z_i (which gives the parcel its initial boost)

      The increment of CAPE (convective air potential energy) for any two
      successive vertical levels is:

          Upwards:
              CAPE_incr = INT(z_0:z_1) g * ( thv_par - thvm ) / thvm dz

          Downwards:
              CAPE_incr = INT(z_(-1):z_0) g * ( thv_par - thvm ) / thvm dz

      Thus, the derivative of CAPE with respect to height is:

              dCAPE/dz = g * ( thv_par - thvm ) / thvm.

      A purely trapezoidal rule is used between levels, and is considered
      to vary linearly at all altitudes.  Thus, dCAPE/dz is considered to be
      of the form:  A * (z-zo) + dCAPE/dz|_(z_0),
      where A = ( dCAPE/dz|_(z_1) - dCAPE/dz|_(z_0) ) / ( z_1 - z_0 )

      The integral is evaluated to find the CAPE increment between two
      successive vertical levels.  The result either adds to or depletes
      from the total amount of energy that keeps the parcel ascending/descending.


      IMPORTANT NOTE:
      This subroutine has been optimized by adding precalculations, rearranging
      equations to avoid divides, and modifying the algorithm entirely.
          -Gunther Huebler

      The algorithm previously used looped over every grid level, following a
      a parcel up from its initial grid level to its max. The very nature of this
      algorithm is an N^2

      The analytic zero-entrainment limit above is not implemented: this
      routine reports fatal status for abs(Lscale_mu_entrain) < eps.

    JAX adaptation:
      Fixed-length scans advance independent column/launch pairs by a common
      offset in grid levels. Masks retain the first energy-exhaustion result;
      scalar conditionals skip inactive work without changing scan length.
      Each carry holds parcel properties, remaining TKE, previous buoyancy,
      live/exit masks, and the level/energy/buoyancy data for the stopping root.
      Prefix extrema replace serial endpoint coupling. Forward working state is
      O(columns*levels); no launch-by-target trajectory matrix is returned.
      Checkpoints recompute layer intermediates during reverse differentiation.
      Blocks of eight offsets skip inactive tails. CPU batches above four columns
      use groups of four launches; accelerators evaluate all launches together.

    References:
      Larson's fifth moist length scale; Section 3b, Golaz et al. (2002),
      JAS 59, 3540--3551.
    """
    zlmin = 0.1
    Lscale_sfclyr_depth = 500.0
    invrs_Lscale_sfclyr_depth = one / Lscale_sfclyr_depth

    # Error in grid column i -> set ith entry to clubb_fatal_error.
    mu_bad = jnp.abs(Lscale_mu_entrain) < eps
    err_info = err_info.set_fatal(mask=mu_bad)
    def print_mu_error(values, info):
        """Print source column headers and invalid entrainment values at the host boundary."""
        import numpy as np
        for i in np.flatnonzero(np.abs(values) < eps):
            print(info.column_error_header_host(int(i)), file=sys.stderr)
            print("Lscale_mu_entrain = ", values[i], file=sys.stderr)
        print("Entrainment rate Lscale_mu_entrain cannot be 0", file=sys.stderr)
        print("Fatal error in subroutine compute_mixing_length", file=sys.stderr)
    jax.lax.cond(
        jnp.any(mu_bad),
        lambda _: jax.debug.callback(print_mu_error, Lscale_mu_entrain, err_info, ordered=True),
        lambda _: None, operand=None,
    )
    # JAX adaptation: preserve healthy columns rather than host-returning the
    # batch at debug >= 0. Fatal columns still need finite traced arithmetic.
    Lscale_mu_entrain = jnp.where(mu_bad, one, Lscale_mu_entrain)

    # Calculate initial turbulent kinetic energy for each grid level.
    tke_i = zm2zt(nzm, nzt, ngrdcol, gr, em)

    # Index adaptation: traverse physical height bottom-to-top on either storage
    # orientation. Interpolation above uses the original grid's native weights.
    step = gr.grid_dir_indx
    tke_i = tke_i[:, ::step]
    thvm = thvm[:, ::step]
    thlm = thlm[:, ::step]
    rtm = rtm[:, ::step]
    p_in_Pa = p_in_Pa[:, ::step]
    exner = exner[:, ::step]
    thv_ds = thv_ds[:, ::step]
    zt = gr.zt[:, ::step]
    dzm = gr.grid_dir * gr.dzm[:, ::step]
    invrs_dzm = gr.grid_dir * gr.invrs_dzm[:, ::step]

    # Initialize arrays and precalculate values for computational efficiency.
    Lscale_up = jnp.full((ngrdcol, nzt), zlmin)
    Lscale_down = jnp.full((ngrdcol, nzt), zlmin)
    grav_on_thvm = grav / thvm
    Lv_coef = Lv / (exner * Cp) - ep2 * thv_ds
    exp_mu_dzm = jnp.exp(-Lscale_mu_entrain[:, None] * dzm)
    invrs_dzm_on_mu = invrs_dzm / Lscale_mu_entrain[:, None]
    entrain_coef = (one - exp_mu_dzm) * invrs_dzm_on_mu
    Lv2_coef = ep * Lv**2 / (Rd * Cp)

    # ---------------- Upwards Length Scale Calculation ----------------

    # Precalculate values for upward Lscale, useful only beyond the first layer.
    thl_par_j_precalc = thlm[:, 1:-1] - thlm[:, :-2] * exp_mu_dzm[:, 1:nzt-1] \
        - (thlm[:, 1:-1] - thlm[:, :-2]) * entrain_coef[:, 1:nzt-1]
    rt_par_j_precalc = rtm[:, 1:-1] - rtm[:, :-2] * exp_mu_dzm[:, 1:nzt-1] \
        - (rtm[:, 1:-1] - rtm[:, :-2]) * entrain_coef[:, 1:nzt-1]

    # Calculate the initial change in TKE for each level in parallel.
    thl_par_1 = thlm[:, 1:-1] \
        - (thlm[:, 1:-1] - thlm[:, :-2]) * entrain_coef[:, 1:nzt-1]
    tl_par_1 = thl_par_1 * exner[:, 1:-1]
    rt_par_1 = rtm[:, 1:-1] \
        - (rtm[:, 1:-1] - rtm[:, :-2]) * entrain_coef[:, 1:nzt-1]
    # The source start_index selects this same subarray for saturation. Shared
    # saturation supports selectors 1--3; selector 4's lookup table remains
    # explicitly gated in saturation.py, outside this module's ownership.
    rsatl_par_1 = sat_mixrat_liq(p_in_Pa[:, 1:-1], tl_par_1, saturation_formula)
    tl_par_j_sqd = tl_par_1**2
    # Lewellen and Yoh (1993): multiply by tl^2 to avoid a division.
    s_par_1 = (rt_par_1 - rsatl_par_1) * tl_par_j_sqd \
        / (tl_par_j_sqd + Lv2_coef * rsatl_par_1)
    rc_par_1 = jnp.maximum(s_par_1, zero_threshold)
    thv_par_1 = thl_par_1 + ep1 * thv_ds[:, 1:-1] * rt_par_1 \
        + Lv_coef[:, 1:-1] * rc_par_1
    dCAPE_dz_1 = grav_on_thvm[:, 1:-1] * (thv_par_1 - thvm[:, 1:-1])
    # Trapezoidal estimate with dCAPE at the launch point equal to zero.
    CAPE_incr_1 = one_half * dCAPE_dz_1 * dzm[:, 1:nzt-1]

    if nzt > 2:
        k = jnp.arange(nzt - 2)[None, :]
        tke = tke_i[:, :-2] + CAPE_incr_1
        alive = tke > zero
        j_last = jnp.broadcast_to(k + 1, (ngrdcol, nzt - 2))
        initial = (tke, thl_par_1, rt_par_1, dCAPE_dz_1, alive, j_last,
                   jnp.zeros_like(alive), j_last, tke, dCAPE_dz_1, dCAPE_dz_1)

        def upward_launch_group(group):
            """Track one batch of upward launches; return their final first-passage data."""
            k, initial = group
            initial = (*initial, jnp.any(initial[4] & (k + 2 < nzt - 1)))
            def upward_layer(state, offset):
                """Advance the legacy carry by one upward level offset; return carry and
                no scan output.
                """
                tke, thl_par_j, rt_par_j, dCAPE_dz_j_minus_1, alive, j_last, \
                    exhausted, exit_j, tke_exit, dcape_exit_prev, dcape_exit_j = state[:-1]
                j = jnp.broadcast_to(jnp.minimum(k + offset, nzt - 2), state[4].shape)
                active = alive & (k + offset < nzt - 1)

                def advance(_):
                    """Return masked upward parcel updates and the first energy-exhaustion data."""
                    # thl, rt of parcel are conserved except for entrainment in a
                    # layer-linear environment: d(thl_par)/dz = -mu*(thl_par-thl_env).
                    thl_next = jnp.take_along_axis(thl_par_j_precalc, j-1, axis=1)
                    thl_next = thl_next + thl_par_j * jnp.take_along_axis(exp_mu_dzm, j, axis=1)
                    rt_next = jnp.take_along_axis(rt_par_j_precalc, j-1, axis=1) \
                        + rt_par_j * jnp.take_along_axis(exp_mu_dzm, j, axis=1)
                    # Include latent heating using Bougeault (1981), then Lewellen/Yoh.
                    tl_par_j = thl_next * jnp.take_along_axis(exner, j, axis=1)
                    rsatl_par_j = sat_mixrat_liq(
                        jnp.take_along_axis(p_in_Pa, j, axis=1), tl_par_j, saturation_formula,
                    )
                    tl_par_j_sqd = tl_par_j**2
                    s_par_j = (rt_next - rsatl_par_j) * tl_par_j_sqd \
                        / (tl_par_j_sqd + Lv2_coef * rsatl_par_j)
                    rc_par_j = jnp.maximum(s_par_j, zero_threshold)
                    thv_par_j = thl_next + ep1 * jnp.take_along_axis(thv_ds, j, axis=1) * rt_next \
                        + jnp.take_along_axis(Lv_coef, j, axis=1) * rc_par_j
                    dCAPE_dz_j = jnp.take_along_axis(grav_on_thvm, j, axis=1) \
                        * (thv_par_j - jnp.take_along_axis(thvm, j, axis=1))
                    CAPE_incr = one_half * (dCAPE_dz_j + dCAPE_dz_j_minus_1) * jnp.take_along_axis(dzm, j, axis=1)
                    # Preserve the first layer that exhausts TKE; later trips are masked.
                    stops = active & (tke + CAPE_incr <= zero)
                    passes = active & ~stops
                    return (
                        jnp.where(passes, tke + CAPE_incr, tke),
                        jnp.where(passes, thl_next, thl_par_j),
                        jnp.where(passes, rt_next, rt_par_j),
                        jnp.where(passes, dCAPE_dz_j, dCAPE_dz_j_minus_1),
                        alive & ~stops,
                        jnp.where(passes, j, j_last), exhausted | stops,
                        jnp.where(stops, j, exit_j), jnp.where(stops, tke, tke_exit),
                        jnp.where(stops, dCAPE_dz_j_minus_1, dcape_exit_prev),
                        jnp.where(stops, dCAPE_dz_j, dcape_exit_j),
                        jnp.any(passes & (k + offset + 1 < nzt - 1)),
                    ), None

                return jax.lax.cond(state[-1], advance, lambda _: (state, None), operand=None)

            def offset_block(state, start):
                """Skip an inactive eight-offset block, or return the carry after its
                upward steps.
                """
                def advance(s):
                    """Scan eight upward offsets from start and return the final parcel carry."""
                    return jax.lax.scan(jax.checkpoint(upward_layer, prevent_cse=False),
                                        s, start + jnp.arange(8))[0]
                return jax.lax.cond(state[-1], advance, lambda s: s, state), None

            final, _ = jax.lax.scan(offset_block, initial, jnp.arange(2, nzt - 1, 8))
            return final[:-1]

        # CPU launch groups avoid evaluating all stopped parcels whenever one
        # launch still has energy. Keep accelerator launches fully batched.
        if ngrdcol > 4 and jax.default_backend() == "cpu":
            launch_block = 4
            groups = (nzt - 2 + launch_block - 1) // launch_block
            padding = groups * launch_block - (nzt - 2)

            def pack_launches(a):
                """Pad column/launch data and return it as (group, column, launch-within-group)."""
                return jnp.pad(a, ((0, 0), (0, padding))).reshape(
                    a.shape[0], groups, launch_block).transpose(1, 0, 2)

            grouped = jax.tree.map(pack_launches, (k, initial))
            final = jax.lax.map(jax.checkpoint(upward_launch_group, prevent_cse=False), grouped)
            final = jax.tree.map(lambda a: a.transpose(1, 0, 2).reshape(
                ngrdcol, groups * launch_block)[:, :nzt-2], final)
        else:
            final = upward_launch_group((k, initial))
        _, _, _, _, _, j_last, exhausted, j, tke_exit, dCAPE_prev, dCAPE_j = final
        dCAPE_diff = dCAPE_j - dCAPE_prev
        linear = jnp.abs(dCAPE_diff) * two <= jnp.abs(dCAPE_j + dCAPE_prev) * eps
        # Guard operands before selecting: inactive sqrt/division derivatives
        # must not contaminate the first-passage gradient.
        invrs_dCAPE_diff = one / jnp.where(exhausted & ~linear, dCAPE_diff, one)
        frac_linear = -tke_exit / jnp.where(exhausted & linear, dCAPE_j, one)
        discriminant = dCAPE_prev**2 - two * tke_exit * jnp.take_along_axis(invrs_dzm, j, axis=1) * dCAPE_diff
        frac_quadratic = (-dCAPE_prev - sqrt_clipped(discriminant)) \
            * invrs_dCAPE_diff * jnp.take_along_axis(dzm, j, axis=1)
        frac_exhausted = jnp.where(exhausted, jnp.where(linear, frac_linear, frac_quadratic), zero)
        # Full passed distance plus the fractional impassable layer.
        full_length = zlmin + jnp.take_along_axis(zt, j_last, axis=1) - zt[:, :-2] + frac_exhausted
        frac_first = -sqrt_clipped(-two * tke_i[:, :-2] * dzm[:, 1:nzt-1] * dCAPE_dz_1) \
            / jnp.where(dCAPE_dz_1 != zero, dCAPE_dz_1, one)
        values = jnp.where(tke_i[:, :-2] + CAPE_incr_1 > zero, full_length, zlmin + frac_first)

        # If a lower-starting parcel can rise higher, extend this one's endpoint.
        parcel_top = zt[:, :-2] + values
        previous_top = jnp.maximum.accumulate(
            jnp.concatenate((jnp.zeros((ngrdcol, 1)), parcel_top), axis=1), axis=1,
        )[:, :-1]
        values = jnp.where(parcel_top < previous_top, previous_top - zt[:, :-2], values)
        Lscale_up = jnp.concatenate((values, Lscale_up[:, -2:]), axis=1)

    # ---------------- Downwards Length Scale Calculation ----------------

    # Precalculate values for downward Lscale beyond the first layer.
    thl_par_j_precalc = thlm[:, :-1] - thlm[:, 1:] * exp_mu_dzm[:, 1:nzt] \
        - (thlm[:, :-1] - thlm[:, 1:]) * entrain_coef[:, 1:nzt]
    rt_par_j_precalc = rtm[:, :-1] - rtm[:, 1:] * exp_mu_dzm[:, 1:nzt] \
        - (rtm[:, :-1] - rtm[:, 1:]) * entrain_coef[:, 1:nzt]
    thl_par_1 = thlm[:, :-1] - (thlm[:, :-1] - thlm[:, 1:]) * entrain_coef[:, 1:nzt]
    tl_par_1 = thl_par_1 * exner[:, :-1]
    rt_par_1 = rtm[:, :-1] - (rtm[:, :-1] - rtm[:, 1:]) * entrain_coef[:, 1:nzt]
    rsatl_par_1 = sat_mixrat_liq(p_in_Pa[:, :-1], tl_par_1, saturation_formula)
    tl_par_j_sqd = tl_par_1**2
    s_par_1 = (rt_par_1 - rsatl_par_1) * tl_par_j_sqd \
        / (tl_par_j_sqd + Lv2_coef * rsatl_par_1)
    rc_par_1 = jnp.maximum(s_par_1, zero_threshold)
    thv_par_1 = thl_par_1 + ep1 * thv_ds[:, :-1] * rt_par_1 + Lv_coef[:, :-1] * rc_par_1
    dCAPE_dz_1 = grav_on_thvm[:, :-1] * (thv_par_1 - thvm[:, :-1])
    CAPE_incr_1 = one_half * dCAPE_dz_1 * dzm[:, 1:nzt]

    k = jnp.arange(1, nzt)[None, :]
    tke = tke_i[:, 1:] - CAPE_incr_1
    alive = tke > zero
    j_last = jnp.broadcast_to(k - 1, (ngrdcol, nzt - 1))
    initial = (tke, thl_par_1, rt_par_1, dCAPE_dz_1, alive, j_last,
               jnp.zeros_like(alive), j_last, tke, dCAPE_dz_1, dCAPE_dz_1)

    def downward_launch_group(group):
        """Track one batch of downward launches; return their final first-passage data."""
        k, initial = group
        initial = (*initial, jnp.any(initial[4] & (k - 2 >= 0)))
        def downward_layer(state, offset):
            """Advance the legacy carry by one downward level offset; return carry and
            no scan output.
            """
            tke, thl_par_j, rt_par_j, dCAPE_dz_j_plus_1, alive, j_last, \
                exhausted, exit_j, tke_exit, dcape_exit_next, dcape_exit_j = state[:-1]
            j = jnp.broadcast_to(jnp.maximum(k - offset, 0), state[4].shape)
            active = alive & (k - offset >= 0)

            def advance(_):
                """Return masked downward parcel updates and the first energy-exhaustion data."""
                # d(thl_par)/dz and d(rt_par)/dz conserve properties except entrainment.
                thl_next = jnp.take_along_axis(thl_par_j_precalc, j, axis=1) \
                    + thl_par_j * jnp.take_along_axis(exp_mu_dzm, j+1, axis=1)
                rt_next = jnp.take_along_axis(rt_par_j_precalc, j, axis=1) \
                    + rt_par_j * jnp.take_along_axis(exp_mu_dzm, j+1, axis=1)
                tl_par_j = thl_next * jnp.take_along_axis(exner, j, axis=1)
                rsatl_par_j = sat_mixrat_liq(
                    jnp.take_along_axis(p_in_Pa, j, axis=1), tl_par_j, saturation_formula,
                )
                tl_par_j_sqd = tl_par_j**2
                s_par_j = (rt_next - rsatl_par_j) * tl_par_j_sqd \
                    / (tl_par_j_sqd + Lv2_coef * rsatl_par_j)
                rc_par_j = jnp.maximum(s_par_j, zero_threshold)
                thv_par_j = thl_next + ep1 * jnp.take_along_axis(thv_ds, j, axis=1) * rt_next + jnp.take_along_axis(Lv_coef, j, axis=1) * rc_par_j
                dCAPE_dz_j = jnp.take_along_axis(grav_on_thvm, j, axis=1) \
                    * (thv_par_j - jnp.take_along_axis(thvm, j, axis=1))
                CAPE_incr = one_half * (dCAPE_dz_j + dCAPE_dz_j_plus_1) * jnp.take_along_axis(dzm, j+1, axis=1)
                stops = active & (tke - CAPE_incr <= zero)
                passes = active & ~stops
                return (
                    jnp.where(passes, tke - CAPE_incr, tke),
                    jnp.where(passes, thl_next, thl_par_j),
                    jnp.where(passes, rt_next, rt_par_j),
                    jnp.where(passes, dCAPE_dz_j, dCAPE_dz_j_plus_1), alive & ~stops,
                    jnp.where(passes, j, j_last), exhausted | stops,
                    jnp.where(stops, j, exit_j), jnp.where(stops, tke, tke_exit),
                    jnp.where(stops, dCAPE_dz_j_plus_1, dcape_exit_next),
                    jnp.where(stops, dCAPE_dz_j, dcape_exit_j),
                    jnp.any(passes & (k - offset - 1 >= 0)),
                ), None

            return jax.lax.cond(state[-1], advance, lambda _: (state, None), operand=None)

        def offset_block(state, start):
            """Skip an inactive eight-offset block, or return the carry after its downward steps."""
            def advance(s):
                """Scan eight downward offsets from start and return the final parcel carry."""
                return jax.lax.scan(jax.checkpoint(downward_layer, prevent_cse=False),
                                    s, start + jnp.arange(8))[0]
            return jax.lax.cond(state[-1], advance, lambda s: s, state), None

        final, _ = jax.lax.scan(offset_block, initial, jnp.arange(2, nzt, 8))
        return final[:-1]

    if ngrdcol > 4 and jax.default_backend() == "cpu":
        launch_block = 4
        groups = (nzt - 1 + launch_block - 1) // launch_block
        padding = groups * launch_block - (nzt - 1)

        def pack_launches(a):
            """Pad column/launch data and return it as (group, column, launch-within-group)."""
            return jnp.pad(a, ((0, 0), (0, padding))).reshape(
                a.shape[0], groups, launch_block).transpose(1, 0, 2)

        grouped = jax.tree.map(pack_launches, (k, initial))
        final = jax.lax.map(jax.checkpoint(downward_launch_group, prevent_cse=False), grouped)
        final = jax.tree.map(lambda a: a.transpose(1, 0, 2).reshape(
            ngrdcol, groups * launch_block)[:, :nzt-1], final)
    else:
        final = downward_launch_group((k, initial))
    _, _, _, _, _, j_last, exhausted, j, tke_exit, dCAPE_next, dCAPE_j = final
    dCAPE_diff = dCAPE_j - dCAPE_next
    linear = jnp.abs(dCAPE_diff) * two <= jnp.abs(dCAPE_j + dCAPE_next) * eps
    invrs_dCAPE_diff = one / jnp.where(exhausted & ~linear, dCAPE_diff, one)
    frac_linear = tke_exit / jnp.where(exhausted & linear, dCAPE_j, one)
    discriminant = dCAPE_next**2 + two * tke_exit * jnp.take_along_axis(invrs_dzm, j+1, axis=1) * dCAPE_diff
    frac_quadratic = (-dCAPE_next + sqrt_clipped(discriminant)) \
        * invrs_dCAPE_diff * jnp.take_along_axis(dzm, j+1, axis=1)
    frac_exhausted = jnp.where(exhausted, jnp.where(linear, frac_linear, frac_quadratic), zero)
    full_length = zlmin + zt[:, 1:] - jnp.take_along_axis(zt, j_last, axis=1) + frac_exhausted
    frac_first = sqrt_clipped(two * tke_i[:, 1:] * dzm[:, 1:nzt] * dCAPE_dz_1) \
        / jnp.where(dCAPE_dz_1 != zero, dCAPE_dz_1, one)
    values = jnp.where(tke_i[:, 1:] - CAPE_incr_1 > zero, full_length, zlmin + frac_first)
    # Extend downward endpoints from higher launches, starting at the top lid.
    parcel_bottom = zt[:, 1:] - values
    previous_bottom = jnp.minimum.accumulate(
        jnp.concatenate((zt[:, -1:], parcel_bottom[:, ::-1]), axis=1), axis=1,
    )[:, :-1][:, ::-1]
    values = jnp.where(parcel_bottom > previous_bottom, zt[:, 1:] - previous_bottom, values)
    Lscale_down = jnp.concatenate((Lscale_down[:, :1], values), axis=1)

    # ---------------- Final Lscale Calculation ----------------

    # Make lminh linear from lmin at the bottom to zero at 500 m altitude.
    if l_implemented:
        # Within a host model, measure height above ground.
        lminh = jnp.maximum(zero_threshold, Lscale_sfclyr_depth - (zt - gr.zm[:, gr.k_lb_zm, None])) \
            * lmin * invrs_Lscale_sfclyr_depth
    else:
        # In standalone mode, measure height above mean sea level.
        lminh = jnp.maximum(zero_threshold, Lscale_sfclyr_depth - zt) * lmin * invrs_Lscale_sfclyr_depth
    Lscale_up = jnp.maximum(lminh, Lscale_up)
    Lscale_down = jnp.maximum(lminh, Lscale_down)
    # A geometric mean becomes small if either directional length is small.
    Lscale = jnp.sqrt(Lscale_up * Lscale_down)
    Lscale = Lscale.at[:, -1].set(Lscale[:, -2])
    # Allow the host model to take over deep convection (Larson, 13 Feb 2008).
    Lscale = jnp.minimum(Lscale, Lscale_max[:, None])
    Lscale, Lscale_up, Lscale_down = Lscale[:, ::step], Lscale_up[:, ::step], Lscale_down[:, ::step]

    # Ensure that no Lscale values are NaN. Keep the JAX finite predicate
    # per column; the source array dump runs only at the host I/O boundary.
    if clubb_at_least_debug_level(1):
        invalid_length = (~jnp.all(jnp.isfinite(Lscale), axis=1)
                          | ~jnp.all(jnp.isfinite(Lscale_up), axis=1)
                          | ~jnp.all(jnp.isfinite(Lscale_down), axis=1))
        err_info = err_info.set_fatal(mask=invalid_length)
        def print_length_errors(info, inputs, outputs):
            """Print source nonfinite-field messages and full bad-column arrays on the host."""
            import numpy as np
            from clubb_jax.src.CLUBB_core.numerical_check import length_check
            with np.printoptions(threshold=np.inf):
                for i in range(ngrdcol):
                    length_check(outputs["Lscale"][i], outputs["Lscale_up"][i], outputs["Lscale_down"][i])
                    if info.err_code[i] != CLUBB_FATAL_ERROR:
                        continue
                    print(info.column_error_header_host(i), file=sys.stderr)
                    print("Errors in compute_mixing_length subroutine in grid column ", i+1, file=sys.stderr)
                    print("Intent(in)", file=sys.stderr)
                    for name in ("thvm", "thlm", "rtm", "em", "exner", "p_in_Pa", "thv_ds"):
                        print(name + " = ", inputs[name][i, :], file=sys.stderr)
                    print("Intent(out)", file=sys.stderr)
                    for name in ("Lscale", "Lscale_up", "Lscale_down"):
                        print(name + " = ", outputs[name][i, :], file=sys.stderr)
        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(
                print_length_errors, err_info,
                {"thvm": thvm[:, ::step], "thlm": thlm[:, ::step],
                 "rtm": rtm[:, ::step], "em": em, "exner": exner[:, ::step],
                 "p_in_Pa": p_in_Pa[:, ::step], "thv_ds": thv_ds[:, ::step]},
                {"Lscale": Lscale, "Lscale_up": Lscale_up, "Lscale_down": Lscale_down},
                ordered=True,
            ),
            lambda _: None, operand=None,
        )
    return err_info, Lscale, Lscale_up, Lscale_down


#===============================================================================

def calc_Lscale_directly(
    ngrdcol, nzm, nzt, gr,
    l_implemented, p_in_Pa, exner, rtm,
    thlm, thvm, rtp2_zt, thlp2_zt, rtpthlp_zt,
    pdf_params, em, thv_ds_zt, Lscale_max, lmin,
    clubb_params,
    saturation_formula,
    l_Lscale_plume_centered,
    Lscale_method,
    stats, err_info,
):
    """Compute the legacy parcel mixing length (method 1).

    Environmental profiles, momentum-level TKE, PDF properties, and per-column
    parameters supply the launches. Return ``(err_info, Lscale, Lscale_up,
    Lscale_down, stats)`` with lengths in m. Retain both returned state objects.
    The source constant ``l_avg_Lscale`` enables optional perturbed launches;
    those calls precede the mean call so directional outputs describe the mean.
    """
    # This mirrors the Fortran compile-time parameter, not a runtime option.
    l_avg_Lscale = False
    if l_avg_Lscale and Lscale_method != Lscale_direct:
        jax.debug.callback(
            lambda: print("calc_Lscale_directly: l_avg_Lscale requires Lscale_method=1.",
                          file=sys.stderr), ordered=True,
        )
        err_info = err_info.set_fatal()
        return (err_info, jnp.zeros((ngrdcol, nzt)), jnp.zeros((ngrdcol, nzt)),
                jnp.zeros((ngrdcol, nzt)), stats)
    if clubb_at_least_debug_level(0):
        if l_Lscale_plume_centered and not l_avg_Lscale:
            # l_Lscale_plume_centered requires l_avg_Lscale.
            jax.debug.callback(
                lambda info: print(info.global_error_header_host(),
                                   "l_Lscale_plume_centered requires l_avg_Lscale",
                                   "Fatal error in calc_Lscale_directly",
                                   sep="\n", file=sys.stderr),
                err_info, ordered=True,
            )
            err_info = err_info.set_fatal()
            return (err_info, jnp.zeros((ngrdcol, nzt)), jnp.zeros((ngrdcol, nzt)),
                    jnp.zeros((ngrdcol, nzt)), stats)

    if l_avg_Lscale and not l_Lscale_plume_centered:
        # Call compute length twice with perturbed rtm/thlm before the mean call.
        sign_rtpthlp_zt = jnp.where(rtpthlp_zt >= zero, one, -one)
        rtm_pert_1 = rtm + clubb_params[:, iLscale_pert_coef, None] * jnp.sqrt(jnp.maximum(rtp2_zt, rt_tol**2))
        thlm_pert_1 = thlm + sign_rtpthlp_zt * clubb_params[:, iLscale_pert_coef, None] * jnp.sqrt(jnp.maximum(thlp2_zt, thl_tol**2))
        mu_pert_1 = clubb_params[:, iLscale_mu_entrain] / clubb_params[:, iLscale_mu_coef]
        err_info, Lscale_pert_1, Lscale_up, Lscale_down = compute_mixing_length(
            nzm, nzt, ngrdcol, gr, thvm, thlm_pert_1,  # In
            rtm_pert_1, em, Lscale_max, p_in_Pa,  # In
            exner, thv_ds_zt, mu_pert_1, lmin,  # In
            saturation_formula, l_implemented,  # In
            err_info,  # InOut
        )
        rtm_pert_2 = rtm - clubb_params[:, iLscale_pert_coef, None] * jnp.sqrt(jnp.maximum(rtp2_zt, rt_tol**2))
        thlm_pert_2 = thlm - sign_rtpthlp_zt * clubb_params[:, iLscale_pert_coef, None] * jnp.sqrt(jnp.maximum(thlp2_zt, thl_tol**2))
        mu_pert_2 = clubb_params[:, iLscale_mu_entrain] * clubb_params[:, iLscale_mu_coef]
        err_info, Lscale_pert_2, Lscale_up, Lscale_down = compute_mixing_length(
            nzm, nzt, ngrdcol, gr, thvm, thlm_pert_2,  # In
            rtm_pert_2, em, Lscale_max, p_in_Pa,  # In
            exner, thv_ds_zt, mu_pert_2, lmin,  # In
            saturation_formula, l_implemented,  # In
            err_info,  # InOut
        )
    elif l_avg_Lscale and l_Lscale_plume_centered:
        # Take values of thl/rt from the wetter/drier PDF plume respectively.
        sign_rtpthlp_zt = jnp.where(rtpthlp_zt >= zero, one, -one)
        Lscale_pert_coef = clubb_params[:, iLscale_pert_coef, None]
        first_wetter = pdf_params.rt_1 > pdf_params.rt_2
        rtm_pert_pos_rt = jnp.where(first_wetter,
            pdf_params.rt_1 + Lscale_pert_coef * jnp.sqrt(jnp.maximum(pdf_params.varnce_rt_1, rt_tol**2)),
            pdf_params.rt_2 + Lscale_pert_coef * jnp.sqrt(jnp.maximum(pdf_params.varnce_rt_2, rt_tol**2)))
        thlm_pert_pos_rt = jnp.where(first_wetter,
            pdf_params.thl_1 + sign_rtpthlp_zt * Lscale_pert_coef * jnp.sqrt(jnp.maximum(pdf_params.varnce_thl_1, thl_tol**2)),
            pdf_params.thl_2 + sign_rtpthlp_zt * Lscale_pert_coef * jnp.sqrt(jnp.maximum(pdf_params.varnce_thl_2, thl_tol**2)))
        thlm_pert_neg_rt = jnp.where(first_wetter,
            pdf_params.thl_2 - sign_rtpthlp_zt * Lscale_pert_coef * jnp.sqrt(jnp.maximum(pdf_params.varnce_thl_2, thl_tol**2)),
            pdf_params.thl_1 - sign_rtpthlp_zt * Lscale_pert_coef * jnp.sqrt(jnp.maximum(pdf_params.varnce_thl_1, thl_tol**2)))
        rtm_pert_neg_rt = jnp.where(first_wetter,
            pdf_params.rt_2 - Lscale_pert_coef * jnp.sqrt(jnp.maximum(pdf_params.varnce_rt_2, rt_tol**2)),
            pdf_params.rt_1 - Lscale_pert_coef * jnp.sqrt(jnp.maximum(pdf_params.varnce_rt_1, rt_tol**2)))
        mu_pert_pos_rt = clubb_params[:, iLscale_mu_entrain] / clubb_params[:, iLscale_mu_coef]
        mu_pert_neg_rt = clubb_params[:, iLscale_mu_entrain] * clubb_params[:, iLscale_mu_coef]
        err_info, Lscale_pert_1, Lscale_up, Lscale_down = compute_mixing_length(
            nzm, nzt, ngrdcol, gr, thvm, thlm_pert_pos_rt,  # In
            rtm_pert_pos_rt, em, Lscale_max, p_in_Pa,  # In
            exner, thv_ds_zt, mu_pert_pos_rt, lmin,  # In
            saturation_formula, l_implemented,  # In
            err_info,  # InOut
        )
        err_info, Lscale_pert_2, Lscale_up, Lscale_down = compute_mixing_length(
            nzm, nzt, ngrdcol, gr, thvm, thlm_pert_neg_rt,  # In
            rtm_pert_neg_rt, em, Lscale_max, p_in_Pa,  # In
            exner, thv_ds_zt, mu_pert_neg_rt, lmin,  # In
            saturation_formula, l_implemented,  # In
            err_info,  # InOut
        )
    else:
        Lscale_pert_1 = jnp.full((ngrdcol, nzt), unused_var)  # Undefined
        Lscale_pert_2 = jnp.full((ngrdcol, nzt), unused_var)  # Undefined

    if stats.l_sample:
        stats = stats.update("Lscale_pert_1", Lscale_pert_1)
        stats = stats.update("Lscale_pert_2", Lscale_pert_2)

    # ********** NOTE: **********
    # The mean call must be last, so directional stats describe the mean length.
    if Lscale_method == Lscale_direct:
        err_info, Lscale, Lscale_up, Lscale_down = compute_mixing_length(
            nzm, nzt, ngrdcol, gr, thvm, thlm,  # In
            rtm, em, Lscale_max, p_in_Pa,  # In
            exner, thv_ds_zt, clubb_params[:, iLscale_mu_entrain], lmin,  # In
            saturation_formula, l_implemented,  # In
            err_info,  # InOut
        )
    else:
        jax.debug.callback(
            lambda: print("calc_Lscale_directly cannot diagnose the selected Lscale_method.",
                          file=sys.stderr), ordered=True,
        )
        err_info = err_info.set_fatal()
        Lscale = jnp.zeros((ngrdcol, nzt))
        Lscale_up = jnp.zeros((ngrdcol, nzt))
        Lscale_down = jnp.zeros((ngrdcol, nzt))

    if l_avg_Lscale:
        if l_Lscale_plume_centered:
            Lscale = one_half * (Lscale_pert_1 + Lscale_pert_2)
        else:
            Lscale = one_third * (Lscale + Lscale_pert_1 + Lscale_pert_2)
    return err_info, Lscale, Lscale_up, Lscale_down, stats


#===============================================================================

def diagnose_Lscale_from_tau(
    nzm, nzt, ngrdcol, gr,
    upwp_sfc, vpwp_sfc, ddzt_umvm_sqd,
    ice_supersat_frac,
    em,
    ufmin, tau_const,
    sfc_elevation, Lscale_max,
    clubb_params,
    stats,
    l_e3sm_config,
    l_smooth_Heaviside_tau_wpxp,
    brunt_vaisala_freq_sqd_smth, Ri_zm,
    err_info,
):
    """Diagnose damping times and mixing length from stress, shear, and stability.

    Method 2 reads surface stresses, momentum-level TKE and gradient profiles,
    thermodynamic-level ice fractions, and per-column tunable parameters.
    Return ``err_info``, the separate inverse-time contributions [1/s], bounded
    time scales [s], combined/directional lengths [m], and ``stats`` in source
    output order. Directional lengths are zero placeholders for this method.
    Suffixes ``_zt`` and ``_zm`` identify thermodynamic and momentum grids.

    References:
      Guo et al. (2021, JAMES).
    """
    l_smooth_min_max = False  # Fortran compile-time parameter; retain both branches.
    heaviside_smth_range = 1.0

    em_zt = zm2zt(nzm, nzt, ngrdcol, gr, em, em_min)

    # Error in grid column i -> set ith entry to clubb_fatal_error
    lowest_zm_below_ground = (
        gr.zm[:, gr.k_lb_zm]
        - sfc_elevation
        + clubb_params[:, iz_displace]
        < eps
    )
    err_info = err_info.set_fatal(mask=lowest_zm_below_ground)
    def print_below_ground(mask, info):
        """Print source error headers for columns whose lowest momentum level is below ground."""
        import numpy as np
        for i in np.flatnonzero(mask):
            print(info.column_error_header_host(int(i)), file=sys.stderr)
    jax.lax.cond(
        jnp.any(lowest_zm_below_ground),
        lambda _: jax.debug.callback(print_below_ground, lowest_zm_below_ground, err_info, ordered=True),
        lambda _: None, operand=None,
    )
    if clubb_at_least_debug_level(0):
        jax.lax.cond(
            err_info.any_fatal(),
            lambda _: jax.debug.callback(lambda: print(
                "Lowest zm grid level is below ground in CLUBB.",
                file=sys.stderr,
            ), ordered=True),
            lambda _: None, operand=None,
        )

    surface_stress = upwp_sfc**2 + vpwp_sfc**2
    # Guard the inactive fourth-root operand so zero stress has a finite
    # reverse derivative before applying the friction-velocity floor.
    tmp_calc_ngrdcol = jnp.where(surface_stress > zero,
        jnp.where(surface_stress > zero, surface_stress, one)**one_fourth, zero)

    if l_smooth_min_max:
        # tmp_calc_ngrdcol used here because openacc has an issue with
        #  a variable being both input and output of a function
        # JAX smooth_max shares one metadata signature across array ranks;
        # the source rank-one overload supplies ngrdcol in both size slots.
        ustar = smooth_max(
            ngrdcol, ngrdcol, tmp_calc_ngrdcol, ufmin, min_max_smth_mag,
        )
    else:
        ustar = jnp.maximum(tmp_calc_ngrdcol, ufmin)

    invrs_tau_bkgnd = (
        clubb_params[:, iC_invrs_tau_bkgnd][:, None] / tau_const
        + jnp.zeros((ngrdcol, nzm))
    )

    norm_ddzt_umvm = sqrt_clipped(ddzt_umvm_sqd)
    smooth_norm_ddzt_umvm = zm2zt2zm(nzm, nzt, ngrdcol, gr, norm_ddzt_umvm)

    invrs_tau_shear_smooth = (
        clubb_params[:, iC_invrs_tau_shear][:, None] * smooth_norm_ddzt_umvm
    )

    # Enforce that invrs_tau_shear is positive.
    invrs_tau_shear = smooth_max(
        nzm, ngrdcol, invrs_tau_shear_smooth,
        zero_threshold, min_max_smth_mag,
    )

    z_eff = (
        gr.zm
        - sfc_elevation[:, None]
        + clubb_params[:, iz_displace][:, None]
    )
    # Fatal ground-level columns still pass through traced arithmetic; guard
    # the denominator shared by the surface and E3SM inverse-time terms.
    z_eff = jnp.where(jnp.abs(z_eff) < eps, eps, z_eff)
    invrs_tau_sfc = (
        clubb_params[:, iC_invrs_tau_sfc][:, None]
        * (ustar[:, None] / vonk)
        / z_eff
    )
    # C_invrs_tau_sfc * ( wp2 / vonk /ustar ) / ( gr%zm(1,:) -sfc_elevation + z_displace )

    invrs_tau_no_N2_zm = invrs_tau_bkgnd + invrs_tau_sfc + invrs_tau_shear

    if l_smooth_min_max:
        brunt_vaisala_freq_clipped = smooth_max(
            nzm, ngrdcol, zero_threshold,
            brunt_vaisala_freq_sqd_smth,
            1.0e-4 * min_max_smth_mag,
        )
        brunt_freq_pos = jnp.sqrt(brunt_vaisala_freq_clipped)
    else:
        brunt_freq_pos = sqrt_clipped(jnp.maximum(zero_threshold, brunt_vaisala_freq_sqd_smth))

    ice_supersat_frac_zm = zt2zm(
        nzm, nzt, ngrdcol, gr, ice_supersat_frac, zero_threshold,
    )

    if l_smooth_min_max:
        tmp_calc = one - ice_supersat_frac_zm / 0.001
        tmp_calc_max = smooth_max(nzm, ngrdcol, zero_threshold, tmp_calc,
                                  min_max_smth_mag)
        tmp_calc_min_max = smooth_min(nzm, ngrdcol, one,
                                      tmp_calc_max, min_max_smth_mag)
        brunt_freq_out_cloud = brunt_freq_pos * tmp_calc_min_max
    else:
        brunt_freq_out_cloud = brunt_freq_pos * jnp.minimum(
            one,
            jnp.maximum(zero_threshold, one - ice_supersat_frac_zm / 0.001),
        )

    brunt_freq_out_cloud = jnp.where(
        gr.zm < clubb_params[:, ialtitude_threshold][:, None],
        zero,
        brunt_freq_out_cloud,
    )

    # Write both bv extra terms to invrs_taus to disk
    if stats.l_sample:
        stats = stats.update("bv_freq_pos", brunt_freq_pos)
        stats = stats.update("bv_freq_out_cloud", brunt_freq_out_cloud)

    # This time scale is used optionally for the return-to-isotropy term. It
    # omits invrs_tau_sfc based on the rationale that the isotropization
    # rate shouldn't be enhanced near the ground.
    invrs_tau_N2_iso = (
        invrs_tau_bkgnd
        + invrs_tau_shear
        + clubb_params[:, iC_invrs_tau_N2_wp2][:, None] * brunt_freq_pos
    )

    invrs_tau_wp2_zm = (
        invrs_tau_no_N2_zm
        + clubb_params[:, iC_invrs_tau_N2][:, None] * brunt_freq_pos
        + clubb_params[:, iC_invrs_tau_N2_wp2][:, None] * brunt_freq_out_cloud
    )

    invrs_tau_zm = (
        invrs_tau_no_N2_zm
        + clubb_params[:, iC_invrs_tau_N2][:, None] * brunt_freq_pos
    )

    if l_e3sm_config:
        invrs_tau_zm = one_half * invrs_tau_zm
        invrs_tau_xp2_zm = (
            invrs_tau_no_N2_zm
            + clubb_params[:, iC_invrs_tau_N2_xp2][:, None] * brunt_freq_pos
            + clubb_params[:, iC_invrs_tau_sfc][:, None]
            * two
            * sqrt_clipped(em)
            / z_eff
        )
        if l_smooth_min_max:
            brunt_vaisala_freq_clipped = smooth_max(
                nzm, ngrdcol, 1.0e-7,
                brunt_vaisala_freq_sqd_smth,
                1.0e-4 * min_max_smth_mag,
            )
            tmp_calc = sqrt_clipped(ddzt_umvm_sqd / brunt_vaisala_freq_clipped)
            tmp_calc_max = smooth_max(nzm, ngrdcol, tmp_calc,
                                      0.3, 0.3 * min_max_smth_mag)
            tmp_calc_min_max = smooth_min(nzm, ngrdcol, tmp_calc_max,
                                          one, min_max_smth_mag)
            invrs_tau_xp2_zm = tmp_calc_min_max * invrs_tau_xp2_zm
        else:
            invrs_tau_xp2_zm = (
                jnp.minimum(
                    jnp.maximum(
                        sqrt_clipped(
                            ddzt_umvm_sqd
                            / jnp.maximum(1.0e-7, brunt_vaisala_freq_sqd_smth)
                        ),
                        0.3,
                    ),
                    one,
                )
                * invrs_tau_xp2_zm
            )
        invrs_tau_wpxp_zm = (
            two * invrs_tau_zm
            + clubb_params[:, iC_invrs_tau_N2_wpxp][:, None]
            * brunt_freq_out_cloud
        )
    else:
        # l_e3sm_config = false
        invrs_tau_xp2_zm = (
            invrs_tau_no_N2_zm
            + clubb_params[:, iC_invrs_tau_N2][:, None] * brunt_freq_pos
            + clubb_params[:, iC_invrs_tau_N2_xp2][:, None] * brunt_freq_out_cloud
        )
        ice_supersat_frac_zm = zt2zm(
            nzm, nzt, ngrdcol, gr, ice_supersat_frac, zero_threshold,
        )
        invrs_tau_wpxp_zm = (
            invrs_tau_no_N2_zm
            + clubb_params[:, iC_invrs_tau_N2][:, None] * brunt_freq_pos
            + clubb_params[:, iC_invrs_tau_N2_wpxp][:, None]
            * brunt_freq_out_cloud
        )

    if l_smooth_Heaviside_tau_wpxp:
        bvf_thresh = (
            brunt_vaisala_freq_sqd_smth
            / clubb_params[:, iC_invrs_tau_wpxp_N2_thresh][:, None]
            - one
        )
        H_invrs_tau_wpxp_N2 = smooth_heaviside_peskin(
            nzm, ngrdcol, bvf_thresh,
            heaviside_smth_range,
        )
    else:
        # l_smooth_Heaviside_tau_wpxp = .false.
        H_invrs_tau_wpxp_N2 = jnp.where(
            brunt_vaisala_freq_sqd_smth
            > clubb_params[:, iC_invrs_tau_wpxp_N2_thresh][:, None],
            one,
            zero,
        )

    if l_smooth_min_max:
        tmp_calc = clubb_params[:, iC_invrs_tau_wpxp_Ri, None] \
            * Ri_zm**clubb_params[:, iwpxp_Ri_exp, None]
        Ri_zm_smooth = smooth_min(nzm, ngrdcol, tmp_calc,
                                 12.0, 12.0 * min_max_smth_mag)
        invrs_tau_wpxp_zm = jnp.where(
            gr.zm > clubb_params[:, ialtitude_threshold, None],
            invrs_tau_wpxp_zm * (one + H_invrs_tau_wpxp_N2 * Ri_zm_smooth), invrs_tau_wpxp_zm)
    else:
        # Choose a zero boundary derivative for a clipped nonpositive Ri.
        # Guard the base before the power, including differentiation of its exponent.
        Ri_positive = Ri_zm > zero
        Ri_power = jnp.where(Ri_positive,
            jnp.where(Ri_positive, Ri_zm, one)**clubb_params[:, iwpxp_Ri_exp, None],
            jnp.where(clubb_params[:, iwpxp_Ri_exp, None] == zero, one, zero))
        invrs_tau_wpxp_zm = jnp.where(
            gr.zm > clubb_params[:, ialtitude_threshold][:, None],
            invrs_tau_wpxp_zm
            * (
                one
                + H_invrs_tau_wpxp_N2
                * jnp.minimum(
                    clubb_params[:, iC_invrs_tau_wpxp_Ri][:, None]
                    * Ri_power,
                    12.0,
                )
            ),
            invrs_tau_wpxp_zm,
        )

    invrs_tau_wp3_zm = (
        invrs_tau_wp2_zm
        + clubb_params[:, iC_invrs_tau_N2_clear_wp3][:, None]
        * brunt_freq_out_cloud
    )

    # Calculate the maximum allowable value of time-scale tau,
    # which depends of the value of Lscale_max.
    if l_smooth_min_max:
        em_clipped = smooth_max(nzm, ngrdcol, em, em_min, min_max_smth_mag)
        tau_max_zt = Lscale_max[:, None] / jnp.sqrt(em_zt)
        tau_max_zm = Lscale_max[:, None] / jnp.sqrt(em_clipped)
    else:
        tau_max_zt = Lscale_max[:, None] / jnp.sqrt(em_zt)
        tau_max_zm = Lscale_max[:, None] / jnp.sqrt(jnp.maximum(em, em_min))

    if l_smooth_min_max:
        tau_zm_unclipped = one / invrs_tau_zm
        tau_zm = smooth_min(nzm, ngrdcol, tau_zm_unclipped,
                           tau_max_zm, 1.0e3 * min_max_smth_mag)
        tau_zt_unclipped = zm2zt(nzm, nzt, ngrdcol, gr, tau_zm)
        tau_zt = smooth_min(nzt, ngrdcol, tau_zt_unclipped,
                           tau_max_zt, 1.0e3 * min_max_smth_mag)
    else:
        tau_zm = jnp.minimum(one / invrs_tau_zm, tau_max_zm)
        tau_zt = zm2zt(nzm, nzt, ngrdcol, gr, tau_zm)
        tau_zt = jnp.minimum(tau_zt, tau_max_zt)

    invrs_tau_zt = zm2zt(nzm, nzt, ngrdcol, gr, invrs_tau_zm)
    invrs_tau_wp3_zt = zm2zt(nzm, nzt, ngrdcol, gr, invrs_tau_wp3_zm)

    Lscale = tau_zt * jnp.sqrt(em_zt)

    # Lscale_up and Lscale_down aren't calculated with this option. They are
    # set to 0 for stats output.
    Lscale_up = jnp.zeros((ngrdcol, nzt))
    Lscale_down = jnp.zeros((ngrdcol, nzt))

    return (
        err_info,
        invrs_tau_zt,
        invrs_tau_zm,
        invrs_tau_sfc,
        invrs_tau_no_N2_zm,
        invrs_tau_bkgnd,
        invrs_tau_shear,
        invrs_tau_N2_iso,
        invrs_tau_wp2_zm,
        invrs_tau_xp2_zm,
        invrs_tau_wp3_zm,
        invrs_tau_wp3_zt,
        invrs_tau_wpxp_zm,
        tau_max_zm,
        tau_max_zt,
        tau_zm,
        tau_zt,
        Lscale,
        Lscale_up,
        Lscale_down,
        stats,
    )
