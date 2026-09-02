"""JAX port of `src/CLUBB_core/advance_helper_module.F90`.

Description:
  This module contains helper methods for the advance_* modules.

Porting deviations:
- The public routine names mirror the Fortran module where practical, while
  returning arrays instead of mutating Fortran out-arguments.
- The Fortran generic interfaces for ``smooth_min``, ``smooth_max``, and
  ``calc_xpwp`` are represented by JAX functions that rely on broadcasting or
  input rank instead of separate scalar/array overloads.
- Optional Fortran stats side effects are returned as diagnostic arrays from
  JAX helpers so callers can update stats outside JIT-compiled code.

These functions smooth the output of the min function
by introducing a varyingly steep path between the two input variables.
The degree to which smoothing is applied depends on the value of 'smth_coef'.
If 'smth_coef' goes toward 0, the output of the min function will be
       0.5 * ((a+b) - abs(a-b))
If a > b, then this comes out to be b. Likewise, if a < b, abs(a-b)=b-a so we get a.
Increasing the smoothing coefficient will lead to a greater degree of smoothing
in the smooth min and max functions. Generally, the coefficient should roughly scale
with the magnitude of data in the data structure that is to be smoothed, in order to
obtain a sensible degree of smoothing (not too much, not too little).

These functions smooth the output of the max functions
by introducing a varyingly steep path between the two input variables.
The degree to which smoothing is applied depends on the value of 'smth_coef'.
If 'smth_coef' goes toward 0, the output of the max function will be
       0.5 * ((a+b) + abs(a-b))
If a > b, then this comes out to be a. Likewise, if a < b, abs(a-b)=b-a so we get b.
Increasing the smoothing coefficient will lead to a greater degree of smoothing
in the smooth min and max functions. Generally, the coefficient should roughly scale
with the magnitude of data in the data structure that is to be smoothed, in order to
obtain a sensible degree of smoothing (not too much, not too little).
"""

from __future__ import annotations

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.T_in_K_module import thlm2T_in_K
from clubb_jax.src.CLUBB_core.constants_clubb import (
    Cp,
    Lv,
    Rd,
    ep,
    grav,
    min_max_smth_mag,
    one_half,
    rc_tol,
    three,
    w_tol_sqd,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.parameter_indices import (
    iCx_max,
    iCx_min,
    iRichardson_num_max,
    iRichardson_num_min,
    ithlp2_rad_coef,
)
from clubb_jax.src.CLUBB_core.grid_class import (
    ddzt,
    zm2zt,
    zm2zt2zm,
    zt2zm,
)
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq


def set_boundary_conditions_lhs(
    diag_index,
    low_bound,
    high_bound,
    lhs,
    diag_index2=None,
    low_bound2=None,
    high_bound2=None,
):
    """Set boundary conditions for a left-hand side LAPACK matrix.

    References:
      none

    Fortran indices are 1-based; this JAX/Python port converts them to
    0-based indices before updating `lhs`.
    """
    diag = int(diag_index) - 1
    low = int(low_bound) - 1
    high = int(high_bound) - 1

    # Set the lower boundaries for the first variable.
    lhs = lhs.at[:, low].set(0.0).at[diag, low].set(1.0)

    # Set the upper boundaries for the first variable.
    lhs = lhs.at[:, high].set(0.0).at[diag, high].set(1.0)

    if low_bound2 is not None:
        if diag_index2 is None:
            raise ValueError("Boundary index provided without diag_index.")
        low2 = int(low_bound2) - 1
        diag2 = int(diag_index2) - 1
        # Set the lower boundaries for the second variable, if it is provided.
        lhs = lhs.at[:, low2].set(0.0).at[diag2, low2].set(1.0)

    if high_bound2 is not None:
        if diag_index2 is None:
            raise ValueError("Boundary index provided without diag_index.")
        high2 = int(high_bound2) - 1
        diag2 = int(diag_index2) - 1
        # Set the upper boundaries for the second variable, if it is provided.
        lhs = lhs.at[:, high2].set(0.0).at[diag2, high2].set(1.0)

    return lhs


def set_boundary_conditions_rhs(
    low_value,
    low_bound,
    high_value,
    high_bound,
    rhs,
    low_value2=None,
    low_bound2=None,
    high_value2=None,
    high_bound2=None,
):
    """Set boundary conditions for a right-hand side LAPACK vector.

    References:
      None
    """
    # Set the lower and upper bounds for the first variable.
    rhs = rhs.at[int(low_bound) - 1].set(low_value)
    rhs = rhs.at[int(high_bound) - 1].set(high_value)

    if low_bound2 is not None:
        if low_value2 is None:
            raise ValueError("Boundary condition provided without value.")
        # If a lower bound was given for the second variable, set it.
        rhs = rhs.at[int(low_bound2) - 1].set(low_value2)

    if high_bound2 is not None:
        if high_value2 is None:
            raise ValueError("Boundary condition provided without value.")
        # If an upper bound was given for the second variable, set it.
        rhs = rhs.at[int(high_bound2) - 1].set(high_value2)

    return rhs


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def calc_ddzt_umvm_sqd(nzm: int, nzt: int, ngrdcol: int, gr, um, vm):
    """Computes the squared vertical norm of the derivative of the mean
    horizontal wind.

    `um` and `vm` are mean horizontal winds on thermodynamic levels [m/s].
    The result is squared vertical shear of horizontal wind [s^-2].

    References:
      None
    """
    ddzt_um = ddzt(nzm, nzt, ngrdcol, gr, um)
    ddzt_vm = ddzt(nzm, nzt, ngrdcol, gr, vm)
    return ddzt_um ** 2 + ddzt_vm ** 2


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def calc_wp3_on_wp2(nzm: int, nzt: int, ngrdcol: int, gr, wp2, wp3):
    """Computes a smoothed ratio of w'^3 to w'^2 on thermodynamic and momentum
    levels.

    `wp2` is the variance of vertical velocity on momentum levels [m^2/s^2].
    `wp3` is the third moment of vertical velocity on thermodynamic levels
    [m^3/s^3]. The returned ratios are on momentum levels and thermodynamic
    levels [m/s].

    References:
      None
    """
    wp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, wp2, w_tol_sqd)

    # Compute wp3 / wp2 on zt levels. Always use the interpolated value in the
    # denominator since it is less likely to create spikes.
    wp3_on_wp2_zt = wp3 / jnp.maximum(wp2_zt, w_tol_sqd)
    wp3_on_wp2_zt = jnp.where(
        wp3_on_wp2_zt < 0.0,
        jnp.maximum(-1000.0, wp3_on_wp2_zt),
        jnp.minimum(1000.0, wp3_on_wp2_zt),
    )

    wp3_on_wp2 = zt2zm(nzm, nzt, ngrdcol, gr, wp3_on_wp2_zt)
    wp3_on_wp2_zt = zm2zt(nzm, nzt, ngrdcol, gr, wp3_on_wp2)
    return wp3_on_wp2, wp3_on_wp2_zt


@partial(jax.jit, static_argnames=("nzm", "ngrdcol"))
def calc_stability_correction(
    nzm: int,
    ngrdcol: int,
    brunt_vaisala_freq_sqd,
    Lscale_zm,
    em,
    lambda0_stability_coef,
):
    """Stability factor.

    `brunt_vaisala_freq_sqd` is the Brunt-Vaisala frequency squared [1/s^2],
    `Lscale_zm` is the turbulent mixing length interpolated to momentum levels
    [m], and `em` is TKE [m^2/s^2].
    """
    del nzm, ngrdcol
    lambda0_stability = jnp.where(
        brunt_vaisala_freq_sqd > 0.0,
        lambda0_stability_coef[:, None],
        0.0,
    )
    return 1.0 + jnp.minimum(
        lambda0_stability * brunt_vaisala_freq_sqd * Lscale_zm ** 2 / em,
        three,
    )


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "saturation_formula",
        "l_brunt_vaisala_freq_moist",
        "l_use_thvm_in_bv_freq",
        "l_modify_limiters_for_cnvg_test",
    ),
)
def calc_brunt_vaisala_freq_sqd(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr,
    thlm,
    exner,
    rtm,
    rcm,
    p_in_Pa,
    thvm,
    ice_supersat_frac,
    saturation_formula: int,
    l_brunt_vaisala_freq_moist: bool,
    l_use_thvm_in_bv_freq: bool,
    l_modify_limiters_for_cnvg_test: bool,
    bv_efold,
    T0: float,
):
    """Calculate the Brunt-Vaisala frequency squared, N^2.

    `l_brunt_vaisala_freq_moist` uses a different formula for saturated
    atmospheres from Durran and Klemp (1982). `l_use_thvm_in_bv_freq` uses
    virtual potential temperature in the Brunt-Vaisala calculation.

    References:
      ?

    Flag to activate modifications on limiters for convergence test
    (smoothed max and min for Cx_fnc_Richardson in advance_helper_module.F90)
    (remove the clippings on brunt_vaisala_freq_sqd_smth in mixing_length.F90)
    (reduce threshold on limiters for Ri_zm in mixing_length.F90)

    `bv_efold` controls the inverse e-folding of cloud fraction in the mixed
    Brunt-Vaisala frequency. `T0` is the reference temperature, usually 300 K.
    """
    ddzt_thlm = ddzt(nzm, nzt, ngrdcol, gr, thlm)
    thvm_zm = zt2zm(nzm, nzt, ngrdcol, gr, thvm, zero_threshold)
    ddzt_thvm = ddzt(nzm, nzt, ngrdcol, gr, thvm)

    # Dry Brunt-Vaisala frequency
    if l_use_thvm_in_bv_freq:
        brunt_vaisala_freq_sqd = (grav / thvm_zm) * ddzt_thvm
    else:
        brunt_vaisala_freq_sqd = (grav / T0) * ddzt_thlm

    T_in_K = thlm2T_in_K(thlm, exner, rcm)
    T_in_K_zm = zt2zm(nzm, nzt, ngrdcol, gr, T_in_K, zero_threshold)
    rsat = sat_mixrat_liq(p_in_Pa, T_in_K, saturation_formula)
    rsat_zm = zt2zm(nzm, nzt, ngrdcol, gr, rsat, zero_threshold)
    ddzt_rsat = ddzt(nzm, nzt, ngrdcol, gr, rsat)

    thm = thlm + Lv / (Cp * exner) * rcm
    thm_zm = zt2zm(nzm, nzt, ngrdcol, gr, thm, zero_threshold)
    ddzt_thm = ddzt(nzm, nzt, ngrdcol, gr, thm)
    ddzt_rtm = ddzt(nzm, nzt, ngrdcol, gr, rtm)

    brunt_vaisala_freq_sqd_dry = (grav / thm_zm) * ddzt_thm

    # In-cloud Brunt-Vaisala frequency. This is Eq. (36) of Durran and
    # Klemp (1982)
    brunt_vaisala_freq_sqd_moist = grav * (
        (
            (1.0 + Lv * rsat_zm / (Rd * T_in_K_zm))
            / (1.0 + ep * (Lv ** 2) * rsat_zm / (Cp * Rd * T_in_K_zm ** 2))
        )
        * ((1.0 / thm_zm * ddzt_thm) + (Lv / (Cp * T_in_K_zm)) * ddzt_rsat)
        - ddzt_rtm
    )

    ice_supersat_frac_zm = zt2zm(nzm, nzt, ngrdcol, gr, ice_supersat_frac, zero_threshold)
    brunt_vaisala_freq_sqd_mixed = (
        brunt_vaisala_freq_sqd_moist
        + jnp.exp(-bv_efold[:, None] * ice_supersat_frac_zm)
        * (brunt_vaisala_freq_sqd_dry - brunt_vaisala_freq_sqd_moist)
    )

    # The min function below smooths the slope discontinuity in brunt freq
    #  and thereby allows tau to remain large in Sc layers in which thlm may
    #  be slightly stably stratified.
    if l_modify_limiters_for_cnvg_test:
        # Remove the limiters to improve the solution convergence
        brunt_vaisala_freq_sqd_smth = zm2zt2zm(
            nzm, nzt, ngrdcol, gr, brunt_vaisala_freq_sqd_mixed
        )
    else:
        brunt_vaisala_freq_clipped = jnp.minimum(
            brunt_vaisala_freq_sqd_mixed,
            1.0e8 * jnp.abs(brunt_vaisala_freq_sqd_mixed) ** 3,
        )
        brunt_vaisala_freq_sqd_smth = zm2zt2zm(
            nzm, nzt, ngrdcol, gr, brunt_vaisala_freq_clipped
        )

    if l_brunt_vaisala_freq_moist:
        brunt_vaisala_freq_sqd = brunt_vaisala_freq_sqd_moist

    # The Fortran routine updates statistics here when `stats` is present and
    # sampling. This JAX routine returns the diagnostic arrays to the caller so
    # stats side effects stay outside the compiled helper.
    return (
        brunt_vaisala_freq_sqd,
        brunt_vaisala_freq_sqd_mixed,
        brunt_vaisala_freq_sqd_smth,
        brunt_vaisala_freq_sqd_dry,
        brunt_vaisala_freq_sqd_moist,
    )


@partial(jax.jit, static_argnames=("nzm", "ngrdcol"))
def calc_Ri_zm(nzm: int, ngrdcol: int, bv_freq_sqd, shear, lim_bv: float, lim_shear: float):
    """Calculate Richardson number from clipped Brunt-Vaisala and shear.

    Calculate the Richardson number as a quotient of a clipped Brunt-Vaisala frequency
    and a clipped shear

    References:
      ?

    `shear` is usually norm squared of horizontal wind speeds.
    """
    del nzm, ngrdcol
    return jnp.maximum(bv_freq_sqd, lim_bv) / jnp.maximum(shear, lim_shear)


@partial(
    jax.jit,
    static_argnames=(
        "nzm",
        "nzt",
        "ngrdcol",
        "l_use_shear_Richardson",
        "l_modify_limiters_for_cnvg_test",
    ),
)
def compute_Cx_fnc_Richardson(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr,
    Lscale_zm,
    ddzt_umvm_sqd,
    rho_ds_zm,
    brunt_vaisala_freq_sqd,
    brunt_vaisala_freq_sqd_mixed,
    clubb_params,
    l_use_shear_Richardson: bool,
    l_modify_limiters_for_cnvg_test: bool,
):
    """Compute Cx as a function of the Richardson number.

    References:
      cam:ticket:59

    `brunt_vaisala_freq_sqd` is Brunt-Vaisala frequency squared, N^2 [1/s^2].
    `brunt_vaisala_freq_sqd_mixed` is a mixture of dry and moist N^2 [s^-2].
    `ddzt_umvm_sqd` is squared vertical norm of derivative of mean horizontal
    wind speed [s^-2]. `Lscale_zm` and `rho_ds_zm` are accepted for parity with
    the Fortran interface.

    Flag to activate modifications on limiters for convergence test
    (smoothed max and min for Cx_fnc_Richardson in advance_helper_module.F90)
    (remove the clippings on brunt_vaisala_freq_sqd_smth in mixing_length.F90)
    (reduce threshold on limiters for sqrt_Ri_zm in mixing_length.F90)

    The Fortran routine also contains a disabled path that vertically averages
    `Cx_fnc_Richardson` over a distance of `Lscale`. That flag is a compile-time
    `.false.` parameter there, so the JAX port keeps only the active path.
    """
    del Lscale_zm, rho_ds_zm
    Richardson_num_divisor_threshold = 1.0e-6
    invrs_num_div_thresh = 1.0 / Richardson_num_divisor_threshold

    if l_use_shear_Richardson:
        Ri_zm_Cx = calc_Ri_zm(
            nzm,
            ngrdcol,
            brunt_vaisala_freq_sqd_mixed,
            ddzt_umvm_sqd,
            1.0e-7,
            1.0e-7,
        )
    else:
        # Note1: We kind of want this calculation to be done in calc_Ri_zm, as well.
        # But multiplying by the inverse and dividing are not numerically equal.
        # So to preserve BFBness, we decided to keep this as is, for now.
        # Note2: Keep in mind, that this Brunt Vaisala variable can contain negative values.
        Ri_zm_Cx = brunt_vaisala_freq_sqd * invrs_num_div_thresh

    Richardson_num_max = clubb_params[:, iRichardson_num_max][:, None]
    Richardson_num_min = clubb_params[:, iRichardson_num_min][:, None]
    Cx_max = clubb_params[:, iCx_max][:, None]
    Cx_min = clubb_params[:, iCx_min][:, None]
    invrs_min_max_diff = 1.0 / (Richardson_num_max - Richardson_num_min)

    # Cx_fnc_Richardson is interpolated based on the value of Richardson_num
    # The min function ensures that Cx does not exceed Cx_max, regardless of the
    #     value of Richardson_num_max.
    if l_modify_limiters_for_cnvg_test:
        fnc_Richardson = (Ri_zm_Cx - Richardson_num_min) * invrs_min_max_diff
        fnc_Richardson_clipped = smooth_min(
            nzm, ngrdcol, 1.0, fnc_Richardson, min_max_smth_mag
        )
        fnc_Richardson_smooth = smooth_max(
            nzm, ngrdcol, 0.0, fnc_Richardson_clipped, min_max_smth_mag
        )
        # use smoothed max amd min to achive smoothed profile and avoid discontinuities
        Cx_fnc_interp = fnc_Richardson_smooth * (Cx_max - Cx_min) + Cx_min
        Cx_fnc_Richardson = zm2zt2zm(nzm, nzt, ngrdcol, gr, Cx_fnc_interp)
    else:
        Cx_fnc_Richardson = (
            jnp.maximum(
                jnp.minimum(Richardson_num_max, Ri_zm_Cx),
                Richardson_num_min,
            )
            - Richardson_num_min
        ) * invrs_min_max_diff * (Cx_max - Cx_min) + Cx_min

    # On some compilers, roundoff error can result in Cx_fnc_Richardson being
    # slightly outside the range [0,1]. Thus, it is clipped here.
    return jnp.maximum(0.0, jnp.minimum(1.0, Cx_fnc_Richardson))


@partial(jax.jit, static_argnames=("nzm", "ngrdcol", "smth_type"))
def Lscale_width_vert_avg(
    nzm: int,
    ngrdcol: int,
    gr,
    smth_type: int,
    var_profile,
    Lscale_zm,
    rho_ds_zm,
    var_below_ground_value: float,
):
    """Average a profile with a running mean of width Lscale_zm.

    References:
      cam:ticket:59

    The result is a vertically averaged profile on momentum levels. Values below
    ground use `var_below_ground_value`. `var_profile`, `Lscale_zm`, and
    `rho_ds_zm` are on momentum levels. The Fortran comment labels `rho_ds_zm`
    as dry static energy, but its use here is dry static density.
    """
    del nzm, ngrdcol
    if smth_type == 1:
        one_half_avg_width = jnp.maximum(Lscale_zm, 500.0)
    elif smth_type == 2:
        one_half_avg_width = jnp.full_like(var_profile, 60.0)
    else:
        raise ValueError(f"Unsupported Lscale_width_vert_avg smth_type={smth_type}")

    # Pre calculate numerator and denominator terms
    numer_terms = rho_ds_zm * gr.grid_dir * gr.dzm * var_profile
    denom_terms = rho_ds_zm * gr.grid_dir * gr.dzm

    # For every grid level
    # -----------------------------------------------------------------------
    # Hunt down all vertical levels with one_half_avg_width(k) of gr%zm(k).
    #
    # Note: Outdated explanation of version that improves CPU performance
    #       below. Reworked due to it requiring a k dependency. Now we
    #       begin looking for k_avg_upper and k_avg_lower starting at
    #       the kth level.
    #
    # Outdated but potentially useful note:
    #     k_avg_upper and k_avg_lower can be saved each loop iteration, this
    #     reduces iterations beacuse the kth values are likely to be within
    #     one or two grid levels of the k-1th values. Less searching is required
    #     by starting the search at the previous values and incrementing or
    #     decrement as needed.
    #
    # The vectorized JAX form builds the full window mask instead.
    # -----------------------------------------------------------------------
    zm_k = gr.zm[:, :, None]
    zm_j = gr.zm[:, None, :]
    in_window = jnp.abs(zm_k - zm_j) <= one_half_avg_width[:, :, None]

    numer_integral = jnp.sum(numer_terms[:, None, :] * in_window, axis=2)
    denom_integral = jnp.sum(denom_terms[:, None, :] * in_window, axis=2)

    # Compute the number of levels below ground to include.
    #
    # The number of below-ground levels included is equal to the distance
    # below the lowest level spanned by one_half_avg_width(k)
    # divided by the distance between vertical levels below ground; the
    # latter is assumed to be the same as the distance between the first and
    # second vertical levels.
    dz_lb = gr.zm[:, 1:2] - gr.zm[:, 0:1]
    below_dist = one_half_avg_width - (gr.zm - gr.zm[:, 0:1])
    n_below_ground_levels = jnp.where(
        below_dist > 0.0,
        jnp.floor(below_dist / dz_lb),
        0.0,
    )
    numer_integral = (
        numer_integral
        + n_below_ground_levels * denom_terms[:, 0:1] * var_below_ground_value
    )
    denom_integral = denom_integral + n_below_ground_levels * denom_terms[:, 0:1]

    # Add numerator and denominator terms for all above-ground levels
    denom_safe = jnp.where(denom_integral != 0.0, denom_integral, 1.0)
    return numer_integral / denom_safe


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def wp23_term_splat_lhs(
    nzm: int,
    nzt: int,
    ngrdcol: int,
    gr,
    C_wp2_splat,
    brunt_vaisala_freq_sqd_mixed,
    Lscale_zm,
    rho_ds_zm,
):
    """Calculate LHS coefficients for the wp2 and wp3 splatting terms.

    `brunt_vaisala_freq_sqd_mixed` is the mixture of dry and moist
    Brunt-Vaisala frequency squared [s^-2], `Lscale_zm` is the turbulent mixing
    length on momentum levels [m], and `rho_ds_zm` is dry static density on
    momentum levels [kg/m^3].
    """
    below_grnd_val = 0.01
    smth_type = 2

    # Calculate Brunt-Vaisala frequency used for splatting.
    brunt_vaisala_freq_sqd_splat = Lscale_width_vert_avg(
        nzm,
        ngrdcol,
        gr,
        smth_type,
        brunt_vaisala_freq_sqd_mixed,
        Lscale_zm,
        rho_ds_zm,
        below_grnd_val,
    )
    brunt_vaisala_freq_splat_clipped = jnp.sqrt(
        jnp.maximum(0.0, brunt_vaisala_freq_sqd_splat)
    )
    brunt_vaisala_freq_splat_smooth = zm2zt2zm(
        nzm, nzt, ngrdcol, gr, brunt_vaisala_freq_splat_clipped
    )
    brunt_vaisala_freq_splat_clipped_zt = zm2zt(
        nzm, nzt, ngrdcol, gr, brunt_vaisala_freq_splat_clipped
    )

    # Vertical compression of eddies causes gustiness (increase in up2 and vp2).
    lhs_splat_wp2 = C_wp2_splat[:, None] * brunt_vaisala_freq_splat_smooth

    # Vertical compression of eddies also diminishes w'3.
    lhs_splat_wp3 = (
        one_half * three * C_wp2_splat[:, None] * brunt_vaisala_freq_splat_clipped_zt
    )
    # The Fortran routine updates `bv_freq_sqd_splat` stats here when sampling.
    # Return it so the caller can handle that side effect outside this JIT.
    return lhs_splat_wp2, lhs_splat_wp3, brunt_vaisala_freq_sqd_splat


@partial(jax.jit, static_argnames=("nz", "ngrdcol"))
def smooth_min(nz: int, ngrdcol: int, input_var1, input_var2, smth_coef: float):
    """Compute a smoothed version of the min function.

    This single JAX implementation covers the Fortran scalar/array overloads.
    For more details, see the module-level interface notes.

    `input_var1` and `input_var2` have units that vary by caller, and the output
    has the same unit as those inputs. `smth_coef` is the smoothing intensity and
    should be of a similar magnitude to the input data structures.

    References:
      See clubb:ticket:894, updated version: 965
    """
    del nz, ngrdcol
    return one_half * (
        (input_var1 + input_var2)
        - jnp.sqrt((input_var1 - input_var2) ** 2 + smth_coef ** 2)
    )


@partial(jax.jit, static_argnames=("nz", "ngrdcol"))
def smooth_max(nz: int, ngrdcol: int, input_var1, input_var2, smth_coef: float):
    """Compute a smoothed version of the max function.

    This single JAX implementation covers the Fortran scalar/array overloads.
    For more details, see the module-level interface notes.

    `input_var1` and `input_var2` have units that vary by caller, and the output
    has the same unit as those inputs. `smth_coef` is the smoothing intensity and
    should be of a similar magnitude to the input data structures.

    References:
      See clubb:ticket:894, updated version: 965
    """
    del nz, ngrdcol
    return one_half * (
        (input_var1 + input_var2)
        + jnp.sqrt((input_var1 - input_var2) ** 2 + smth_coef ** 2)
    )


@partial(jax.jit, static_argnames=("nz", "ngrdcol"))
def smooth_heaviside_peskin(nz: int, ngrdcol: int, input, smth_range: float):
    """Compute a smoothed Heaviside function.

    Computes a smoothed heaviside function as in
        [Lin, Lee et al., 2005, A level set characteristic Galerkin
        finite element method for free surface flows], equation (2)

    `smth_range` defines the smooth Heaviside interval
    [-smth_range, smth_range]. The output has the same units as `input`.

    References:
      See clubb:ticket:965
    """
    del nz, ngrdcol
    input_over_smth_range = input / smth_range
    interior = one_half * (
        1.0
        + input_over_smth_range
        + (1.0 / jnp.pi) * jnp.sin(jnp.pi * input_over_smth_range)
    )
    # Note that this case will only ever be reached if smth_range != 0,
    # so this division is fine and should not cause any issues
    #
    # JAX evaluates both branches, so callers should use the same nonzero
    # smoothing range assumed by the active formula.
    return jnp.where(input < -smth_range, 0.0, jnp.where(input > smth_range, 1.0, interior))


@partial(jax.jit, static_argnames=("nzm", "nzt", "ngrdcol"))
def calc_xpwp(nzm: int, nzt: int, ngrdcol: int, gr, Km_zm, xm):
    """Compute x'w' from x<k>, x<k-1>, Km_zm, and invrs_dzm.

    References:
      None

    This combines the Fortran `calc_xpwp_1D` and `calc_xpwp_2D` overloads.
    """
    del nzt, ngrdcol
    if Km_zm.ndim == 1:
        xpwp = jnp.zeros((nzm,), dtype=jnp.float64)
        # Solve for x'w' at all intermediate model levels.
        interior = Km_zm[1:nzm - 1] * gr.invrs_dzm[0, 1:nzm - 1] * (xm[1:] - xm[:-1])
        return xpwp.at[1:nzm - 1].set(interior)

    xpwp = jnp.zeros((Km_zm.shape[0], nzm), dtype=jnp.float64)
    # Solve for x'w' at all intermediate model levels.
    interior = Km_zm[:, 1:nzm - 1] * gr.invrs_dzm[:, 1:nzm - 1] * (xm[:, 1:] - xm[:, :-1])
    return xpwp.at[:, 1:nzm - 1].set(interior)


def vertical_avg(total_idx: int, rho_ds, field, dz):
    """Compute the density-weighted vertical average of a field.

    The average value of a function, f, over a set domain, [a,b], is calculated
    by the equation:

      f_avg = ( INT(a:b) f*g ) / ( INT(a:b) g )

    as long as f is continous and g is nonnegative and integrable.  Therefore,
    the density-weighted (by dry, static, base-static density) vertical
    average value of any model field, x, is calculated by the equation:

      x_avg|_z = ( INT(z_bot:z_top) x rho_ds dz )
                 / ( INT(z_bot:z_top) rho_ds dz )

    For numerical purposes, both thermodynamic-level and momentum-level
    computations use sums of x(k) * rho_ds(k) * delta_z(k) for the numerator
    and rho_ds(k) * delta_z(k) for the denominator over the selected vertical
    levels. The caller chooses the included levels by passing arrays arranged
    from lowest to highest altitude and setting `total_idx`.

    Thermodynamic-level computation:
    For numerical purposes, INT(z_bot:z_top) x rho_ds dz, which is the numerator
    integral, is calculated as:

      SUM(k_bot:k_top) x(k) rho_ds(k) delta_z(k)

    where k is the index of the given thermodynamic level, x and rho_ds are both
    thermodynamic-level variables, and delta_z(k) = zm(k) - zm(k-1). The indices
    k_bot and k_top are the indices of the respective lower and upper
    thermodynamic levels involved in the integration.

    Likewise, INT(z_bot:z_top) rho_ds dz, which is the denominator integral, is
    calculated as:

      SUM(k_bot:k_top) rho_ds(k) delta_z(k)

    The first (k=1) thermodynamic level is below ground, or below the official
    lower boundary at the first momentum level, so it should not count in a
    vertical average for hole-filling or statistics. Begin no lower than level
    k=2, the first thermodynamic level above ground. For full-domain
    hole-filling or statistics, lower thermodynamic-level index k=2 and upper
    thermodynamic-level index k=gr%nz means the overall vertical domain is
    gr%zm(1,gr%nz) - gr%zm(1,1).

    Momentum-level computation:
    For numerical purposes, INT(z_bot:z_top) x rho_ds dz is calculated as:

      SUM(k_bot:k_top) x(k) rho_ds(k) delta_z(k)

    where k is the index of the given momentum level, x and rho_ds are both
    momentum-level variables, and delta_z(k) = zt(k+1) - zt(k). The indices
    k_bot and k_top are the respective lower and upper momentum levels involved
    in the integration. The denominator integral is:

      SUM(k_bot:k_top) rho_ds(k) delta_z(k)

    The first (k=1) momentum level is at ground level, or at the official lower
    boundary. Momentum-level variables that call the hole-filling scheme have
    set values at the surface, so those values should not be changed. For
    hole-filling, begin no lower than level k=2, the second momentum level above
    ground. The value at the model upper boundary (k=gr%nz) is also set for
    momentum-level variables and should not be changed.

    This function is also used for statistical vertical averages of certain
    variables. In that case, the vertical average needs to be taken over the
    entire vertical domain, level 1 to level gr%nz.

    In both thermodynamic-level and momentum-level computations, the numerator
    integral is divided by the denominator integral to find the average value of
    x over the vertical domain.

    References:
      None
    """
    # Initialize variable
    # Compute the numerator and denominator integral.
    # Multiply rho_ds at level k by the level thickness
    # at level k.  Then, sum over all vertical levels.
    vertical_avg_result = jnp.sum(rho_ds[:total_idx] * dz[:total_idx] * field[:total_idx]) / jnp.sum(
        rho_ds[:total_idx] * dz[:total_idx]
    )
    # Find the vertical average of 'field'.
    return vertical_avg_result


def vertical_integral(total_idx: int, rho_ds, field, dz):
    """Compute the vertical integral.

    `rho_ds`, `field`, and `dz` must all be of size `total_idx` and should all
    start at the same index. The values need to be arranged from lowest to
    highest altitude, with `field[0]` and `rho_ds[0]` corresponding to the
    original field and density at the caller's starting level.

    References:
      None
    """
    # Assertion checks: that k_start <= gr%nz - 1
    #                   that k_end   >= 2
    #                   that k_start <= k_end

    # Initializing vertical_integral to avoid a compiler warning.
    # Compute the integral.
    # Multiply the field at level k by rho_ds at level k and by
    # the level thickness at level k.  Then, sum over all vertical levels.
    # Note:  The values of the field and rho_ds are passed into this function
    #        so that field(1) and rho_ds(1) are actually the field and rho_ds
    #        at level k_start.
    return jnp.sum(field[:total_idx] * rho_ds[:total_idx] * dz[:total_idx])


def pvertinterp(nzt: int, ngrdcol: int, gr, p_mid, p_out: float, input_var):
    """Interpolate a variable from model pressure levels to target pressure.

    `p_mid` contains input pressure levels, `p_out` is the output pressure
    level, and `input_var` is the input array. The result is the interpolated
    output array.
    """
    del nzt, ngrdcol, gr
    # Initialize index array and logical flags
    # If we've fallen through the k=1,nz-1 loop, we cannot interpolate and
    # must extrapolate from the bottom or top data level for at least some
    # of the longitude points.
    #
    # Store level indices for interpolation.
    # If all indices for this level have been found,
    # do the interpolation
    return jax.vmap(lambda p_col, v_col: jnp.interp(p_out, p_col[::-1], v_col[::-1]))(
        p_mid,
        input_var,
    )


@partial(jax.jit, static_argnames=("ngrdcol", "nzm", "nzt"))
def calculate_thlp2_rad(
    ngrdcol: int,
    nzm: int,
    nzt: int,
    gr,
    rcm,
    thlprcp,
    radht,
    clubb_params,
    thlp2_forcing,
):
    """Computes the contribution of radiative cooling to thlp2.

    References:
      See clubb:ticket:632

    `rcm` is cloud water mixing ratio, `radht` is SW + LW heating rate [K/s],
    `thlprcp` is thl'rc' [K kg/kg], and `thlp2_forcing` is the
    <th_l'^2> forcing on momentum levels [K^2/s].
    """
    rcm_zm = zt2zm(nzm, nzt, ngrdcol, gr, rcm)
    radht_zm = zt2zm(nzm, nzt, ngrdcol, gr, radht)
    increment = (
        clubb_params[:, ithlp2_rad_coef][:, None]
        * 2.0
        * radht_zm
        / jnp.where(rcm_zm > rc_tol, rcm_zm, 1.0)
        * thlprcp
    )
    return jnp.where(rcm_zm > rc_tol, thlp2_forcing + increment, thlp2_forcing)
