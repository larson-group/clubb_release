"""JAX port of `src/CLUBB_core/clip_explicit.F90`.

Porting deviations:
- Fortran inout arguments are returned as values so clipping routines remain
  functional and JIT-friendly.
- Fortran loops over grid columns, levels, and scalar dimensions are expressed
  as JAX array operations or small Python loops over static dimensions.
- Fortran stats side effects are represented by functional updates to the JAX
  stats object.
- Fortran debug printing and hard stops are omitted from jitted routines; the
  active callers enforce thresholds before calling or carry error state
  separately.
"""

from __future__ import annotations

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.advance_helper_module import smooth_heaviside_peskin
from clubb_jax.src.CLUBB_core.constants_clubb import (
    max_mag_correlation,
    max_mag_correlation_flux,
    zero_threshold,
)


# Named constants to avoid string comparisons.
clip_rtp2 = 1         # Named constant for rtp2 clipping
clip_thlp2 = 2        # Named constant for thlp2 clipping
clip_rtpthlp = 3      # Named constant for rtpthlp clipping
clip_up2 = 5          # Named constant for up2 clipping
clip_vp2 = 6          # Named constant for vp2 clipping
clip_wprtp = 8        # Named constant for wprtp clipping
clip_wpthlp = 9       # Named constant for wpthlp clipping
clip_upwp = 10        # Named constant for upwp clipping
clip_vpwp = 11        # Named constant for vpwp clipping
clip_wp2 = 12         # Named constant for wp2 clipping
clip_wpsclrp = 13     # Named constant for wp scalar clipping
clip_sclrp2 = 14      # Named constant for sclrp2 clipping
clip_sclrprtp = 15    # Named constant for sclrprtp clipping
clip_sclrpthlp = 16   # Named constant for sclrpthlp clipping
clip_wphydrometp = 17 # Named constant for wphydrometp clipping

wprtp_cl_max = 3   # Number of wprtp clipping events per timestep
wpthlp_cl_max = 3  # Number of wpthlp clipping events per timestep
upwp_cl_max = 3    # Number of upwp clipping events per timestep
vpwp_cl_max = 3    # Number of vpwp clipping events per timestep

_F64_EPS = jnp.finfo(jnp.float64).eps


@partial(jax.jit, static_argnames=("nzm", "ngrdcol", "solve_type"))
def clip_covar(nzm: int, ngrdcol: int, solve_type: int, xp2, yp2, xpyp):
    """Clip covariance x'y' based on the correlation between x and y.

    The correlation between variables x and y is:

      corr_(x,y) = x'y' / [ sqrt(x'^2) * sqrt(y'^2) ]

    The correlation must stay between -1 and 1. Therefore, x'y' has upper and
    lower limits:

      x'y' <=  [ sqrt(x'^2) * sqrt(y'^2) ]
      x'y' >= -[ sqrt(x'^2) * sqrt(y'^2) ]

    The values of x'y', x'^2, and y'^2 are all found on momentum levels. The
    value of x'y' may need to be clipped whenever x'y', x'^2, or y'^2 is
    updated.

    References:
      None
    """
    del ngrdcol

    # When clipping for wprtp or wpthlp, use the special value for
    # max_mag_correlation_flux. For all other correlations, use
    # max_mag_correlation.
    if solve_type in (clip_wprtp, clip_wpthlp):
        max_mag_corr = max_mag_correlation_flux
    else:
        max_mag_corr = max_mag_correlation

    xpyp_bound = max_mag_corr * jnp.sqrt(xp2 * yp2)
    xpyp_clipped = jnp.clip(xpyp, -xpyp_bound, xpyp_bound)

    # The value of x'y' at the surface (or lower boundary) is set elsewhere,
    # and the value at the upper boundary is also set. Do not clip either
    # boundary. If clipping were applied at the lower boundary, momentum would
    # not be conserved, so it should never be added.
    interior_mask = (jnp.arange(nzm)[None, :] > 0) & (jnp.arange(nzm)[None, :] < nzm - 1)
    xpyp_out = jnp.where(interior_mask, xpyp_clipped, xpyp)
    xpyp_chnge = jnp.where(interior_mask, xpyp_out - xpyp, 0.0)
    return xpyp_out, xpyp_chnge


def clip_variance(
    nzm: int,
    ngrdcol: int,
    gr,
    solve_type: int,
    dt: float,
    threshold_lo,
    stats,
    xp2,
    threshold_hi=None,
):
    """Clip variance x'^2 based on minimum and optional maximum thresholds.

    The threshold value must be greater than or equal to 0. The values of x'^2
    are found on momentum levels.

    The following variances are found in the code:

    r_t'^2, th_l'^2, u'^2, v'^2, sclr'^2, computed in advance_xp2_xpyp; and
    w'^2, computed in advance_wp2_wp3.

    References:
      None
    """
    del ngrdcol, gr

    name_xp2_cl = ""
    if solve_type == clip_wp2:
        name_xp2_cl = "wp2_cl"
    elif solve_type == clip_rtp2:
        name_xp2_cl = "rtp2_cl"
    elif solve_type == clip_thlp2:
        name_xp2_cl = "thlp2_cl"
    elif solve_type == clip_up2:
        name_xp2_cl = "up2_cl"
    elif solve_type == clip_vp2:
        name_xp2_cl = "vp2_cl"
    l_sample = stats.l_sample
    if l_sample and name_xp2_cl.strip():
        stats = stats.begin_budget(name_xp2_cl, xp2 / dt)

    # Limit the value of x'^2 at threshold. The original code did not need to
    # invoke clipping at the surface, because the surface value is set
    # elsewhere. charlass changed this on 09/11/2013 so the surface level is
    # also clipped after slightly negative thlp2(1) and rtp2(1) values appeared
    # in WRF-CLUBB quarter_ss.
    active_mask = jnp.arange(nzm)[None, :] < nzm - 1
    xp2_out = jnp.where(active_mask, jnp.maximum(threshold_lo, xp2), xp2)

    # Optional clipping of large values. This was added so the <var>_cl budget
    # term also contains the upper clipping contributions.
    if threshold_hi is not None:
        xp2_out = jnp.where(active_mask, jnp.minimum(threshold_hi, xp2_out), xp2_out)

    if l_sample and name_xp2_cl.strip():
        stats = stats.finalize_budget(name_xp2_cl, xp2_out / dt)

    return xp2_out, stats


def clip_covars_denom(
    nzm: int,
    ngrdcol: int,
    sclr_dim: int,
    dt: float,
    rtp2,
    thlp2,
    up2,
    vp2,
    wp2,
    sclrp2,
    l_linearize_pbl_winds: bool,
    l_predict_upwp_vpwp: bool,
    stats,
    wprtp_cl_num: int,
    wpthlp_cl_num: int,
    upwp_cl_num: int,
    vpwp_cl_num: int,
    wprtp,
    wpthlp,
    upwp,
    vpwp,
    wpsclrp,
    upwp_pert,
    vpwp_pert,
):
    """Clip covariances after denominator terms have been altered.

    Some of the covariances found in the CLUBB model code need to be clipped
    multiple times during each timestep to ensure that the correlation between
    the two relevant variables stays between -1 and 1 at all times during the
    model run.  The covariances that need to be clipped multiple times are
    w'r_t', w'th_l', w'sclr', u'w', and v'w'.  One of the times that each one
    of these covariances is clipped is immediately after each one is set.
    However, each covariance still needs to be clipped two more times during
    each timestep (once after advance_xp2_xpyp is called and once after
    advance_wp2_wp3 is called).  This subroutine handles the times that the
    covariances are clipped away from the time that they are set.  In other
    words, this subroutine clips the covariances after the denominator terms
    in the relevant correlation equation have been altered, ensuring that
    all correlations will remain between -1 and 1 at all times.

    References:
      None
    """
    l_sample = stats.l_sample

    if l_sample:
        if wprtp_cl_num == 0:
            stats = stats.begin_budget("wprtp_cl", wprtp / dt)
        else:
            stats = stats.update_budget("wprtp_cl", -wprtp / dt)

        if wpthlp_cl_num == 0:
            stats = stats.begin_budget("wpthlp_cl", wpthlp / dt)
        else:
            stats = stats.update_budget("wpthlp_cl", -wpthlp / dt)

        if l_predict_upwp_vpwp:
            if upwp_cl_num == 0:
                stats = stats.begin_budget("upwp_cl", upwp / dt)
            else:
                stats = stats.update_budget("upwp_cl", -upwp / dt)

            if vpwp_cl_num == 0:
                stats = stats.begin_budget("vpwp_cl", vpwp / dt)
            else:
                stats = stats.update_budget("vpwp_cl", -vpwp / dt)

    wprtp_cl_num += 1
    wpthlp_cl_num += 1
    upwp_cl_num += 1
    vpwp_cl_num += 1

    # Clip w'r_t'. This routine handles the first and third instances of
    # w'r_t' clipping: after r_t'^2 is updated in advance_xp2_xpyp and after
    # w'^2 is updated in advance_wp2_wp3.
    wprtp, _ = clip_covar(nzm, ngrdcol, clip_wprtp, wp2, rtp2, wprtp)

    # Clip w'th_l'.
    wpthlp, _ = clip_covar(nzm, ngrdcol, clip_wpthlp, wp2, thlp2, wpthlp)

    # Clip w'sclr'. This handles the first and third instances of w'sclr'
    # clipping: after sclr'^2 is updated in advance_xp2_xpyp and after w'^2 is
    # updated in advance_wp2_wp3.
    for sclr in range(sclr_dim):
        wpsclrp_s, _ = clip_covar(
            nzm,
            ngrdcol,
            clip_wpsclrp,
            wp2,
            sclrp2[:, :, sclr],
            wpsclrp[:, :, sclr],
        )
        wpsclrp = wpsclrp.at[:, :, sclr].set(wpsclrp_s)

    # Clip u'w'. This handles the first and second instances of u'w' clipping:
    # after u'^2 is updated in advance_xp2_xpyp and after w'^2 is updated in
    # advance_wp2_wp3.
    upwp, _ = clip_covar(nzm, ngrdcol, clip_upwp, wp2, up2, upwp)
    if l_linearize_pbl_winds:
        upwp_pert, _ = clip_covar(nzm, ngrdcol, clip_upwp, wp2, up2, upwp_pert)

    # Clip v'w'. This handles the first and second instances of v'w' clipping:
    # after v'^2 is updated in advance_xp2_xpyp and after w'^2 is updated in
    # advance_wp2_wp3.
    vpwp, _ = clip_covar(nzm, ngrdcol, clip_vpwp, wp2, vp2, vpwp)
    if l_linearize_pbl_winds:
        vpwp_pert, _ = clip_covar(nzm, ngrdcol, clip_vpwp, wp2, vp2, vpwp_pert)

    if l_sample:
        if wprtp_cl_num == wprtp_cl_max:
            stats = stats.finalize_budget("wprtp_cl", wprtp / dt)
        else:
            stats = stats.update_budget("wprtp_cl", wprtp / dt)

        if wpthlp_cl_num == wpthlp_cl_max:
            stats = stats.finalize_budget("wpthlp_cl", wpthlp / dt)
        else:
            stats = stats.update_budget("wpthlp_cl", wpthlp / dt)

        if l_predict_upwp_vpwp:
            if upwp_cl_num == upwp_cl_max:
                stats = stats.finalize_budget("upwp_cl", upwp / dt)
            else:
                stats = stats.update_budget("upwp_cl", upwp / dt)

            if vpwp_cl_num == vpwp_cl_max:
                stats = stats.finalize_budget("vpwp_cl", vpwp / dt)
            else:
                stats = stats.update_budget("vpwp_cl", vpwp / dt)

    return (
        wprtp_cl_num,
        wpthlp_cl_num,
        upwp_cl_num,
        vpwp_cl_num,
        wprtp,
        wpthlp,
        upwp,
        vpwp,
        wpsclrp,
        upwp_pert,
        vpwp_pert,
        stats,
    )


def clip_skewness(
    nzt: int,
    ngrdcol: int,
    gr,
    dt: float,
    sfc_elevation,
    Skw_max_mag,
    wp2_zt,
    l_use_wp3_lim_with_smth_Heaviside: bool,
    stats,
    wp3,
):
    """Clipping the value of w'^3 based on the skewness of w, Sk_w.

    Aditionally, to prevent possible crashes due to wp3 growing too large,
    abs(wp3) will be clipped to 100.

    The skewness of w is:

    Sk_w = w'^3 / (w'^2)^(3/2).

    The value of Sk_w is limited to a range between an upper limit and a lower
    limit.  The values of the limits depend on whether the level altitude is
    within 100 meters of the surface.

    For altitudes less than or equal to 100 meters above ground level (AGL):

    -0.2_core_rknd*sqrt(2) <= Sk_w <= 0.2_core_rknd*sqrt(2);

    while for all altitudes greater than 100 meters AGL:

    -4.5_core_rknd <= Sk_w <= 4.5_core_rknd.

    Therefore, there is an upper limit on w'^3, such that:

    w'^3  <=  threshold_magnitude * (w'^2)^(3/2);

    and a lower limit on w'^3, such that:

    w'^3  >= -threshold_magnitude * (w'^2)^(3/2).

    The values of w'^3 are found on the thermodynamic levels, while the values
    of w'^2 are found on the momentum levels.  Therefore, the values of w'^2
    are interpolated to the thermodynamic levels before being used to
    calculate the upper and lower limits for w'^3.

    References:
      None
    """
    l_sample = stats.l_sample
    if l_sample:
        stats = stats.begin_budget("wp3_cl", wp3 / dt)

    wp3 = clip_skewness_core(
        nzt,
        ngrdcol,
        gr,
        sfc_elevation,
        Skw_max_mag,
        wp2_zt,
        l_use_wp3_lim_with_smth_Heaviside,
        wp3,
    )

    if l_sample:
        stats = stats.finalize_budget("wp3_cl", wp3 / dt)

    return wp3, stats


@partial(
    jax.jit,
    static_argnames=("nzt", "ngrdcol", "l_use_wp3_lim_with_smth_Heaviside"),
)
def clip_skewness_core(
    nzt: int,
    ngrdcol: int,
    gr,
    sfc_elevation,
    Skw_max_mag,
    wp2_zt,
    l_use_wp3_lim_with_smth_Heaviside: bool,
    wp3,
):
    """Core skewness clipping math for w'^3.

    The skewness of w is:

      Sk_w = w'^3 / (w'^2)^(3/2)

    The value of Sk_w is limited to a range whose limits depend on whether the
    level altitude is within 100 meters of the surface. For altitudes less than
    or equal to 100 meters AGL:

      -0.2*sqrt(2) <= Sk_w <= 0.2*sqrt(2)

    For all altitudes greater than 100 meters AGL:

      -4.5 <= Sk_w <= 4.5
    """
    wp3_max = 100.0

    # Compute the upper and lower limits of w'^3 at every level,
    # based on the skewness of w, Sk_w, such that:
    # Sk_w = w'^3 / (w'^2)^(3/2);
    # -4.5 <= Sk_w <= 4.5;
    # or, if the level altitude is within 100 meters of the surface,
    # -0.2*sqrt(2) <= Sk_w <= 0.2*sqrt(2).

    # The normal magnitude limit of skewness of w in the CLUBB code is 4.5.
    # However, according to Andre et al. (1976b & 1978), wp3 should not exceed
    # [2*(wp2^3)]^(1/2) at any level.  However, this term should be multiplied
    # by 0.2 close to the surface to include surface effects.  There already is
    # a wp3 clipping term in place for all other altitudes, but this term will
    # be included for the surface layer only.  Therefore, the lowest level wp3
    # should not exceed 0.2 * sqrt(2) * wp2^(3/2).  Brian Griffin.  12/18/05.

    # To lower compute time, square both sides of the equation and compute
    # wp2^3 only once.
    wp2_zt_cubed = wp2_zt ** 3

    if l_use_wp3_lim_with_smth_Heaviside:
        #implement a smoothed Heaviside function to avoid discontinuities
        zagl_thresh = (gr.zt - sfc_elevation[:, None]) / 100.0
        zagl_thresh = zagl_thresh - 1.0
        H_zagl = smooth_heaviside_peskin(nzt, ngrdcol, zagl_thresh, 0.6)
        wp3_lim_sqd = wp2_zt_cubed * (
            H_zagl * Skw_max_mag[:, None] ** 2
            + (1.0 - H_zagl) * 0.0021 * Skw_max_mag[:, None] ** 2
        )
    else:
        # Clip for 100 m. AGL.
        # Clip skewness consistently with a.
        wp3_lim_sqd = jnp.where(
            gr.zt - sfc_elevation[:, None] <= 100.0,
            0.0021 * Skw_max_mag[:, None] ** 2 * wp2_zt_cubed,
            Skw_max_mag[:, None] ** 2 * wp2_zt_cubed,
        )

    # Clipping for w'^3 at an upper and lower limit corresponding with
    # the appropriate value of Sk_w.
    # Set the magnitude to the wp3 limit and apply the sign of the current wp3
    wp3 = jnp.where(wp3 ** 2 > wp3_lim_sqd, jnp.sign(wp3) * jnp.sqrt(wp3_lim_sqd), wp3)

    # Clipping abs(wp3) to 100. This keeps wp3 from growing too large in some
    # convective cases, which helps prevent these cases from blowing up.
    return jnp.clip(wp3, -wp3_max, wp3_max)


@partial(jax.jit, static_argnames=("nzt", "ngrdcol", "message"))
def clip_rcm(nzt: int, ngrdcol: int, rtm, message: str, rcm):
    """Reduce cloud water rcm whenever it exceeds total water rtm.

    This avoids negative values of rvm, the water vapor mixing ratio. However,
    it will not ensure that rcm <= rtm if rtm <= 0.

    References:
      None
    """
    del nzt, ngrdcol, message
    # Vince Larson clipped rcm in order to prevent rvm < 0. 5 Apr 2008.
    # This code will not work unless rtm >= 0.
    # We do not clip rcm_in_layer because rcm_in_layer only influences
    # radiation, and we do not want to bother recomputing it. 6 Aug 2009.
    return jnp.where(
        rtm < rcm,
        jnp.maximum(zero_threshold, rtm - _F64_EPS),
        rcm,
    )
