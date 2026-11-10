"""JAX port of new_tsdadg_pdf.F90 - the new TSDADG PDF closure helpers.

Description:
The new trivariate, skewness-dependent, analytic double Gaussian (TSDADG)
PDF.

Porting deviations:
- The Fortran nz loop in `tsdadg_pdf_driver` is expressed as JAX vectorized
  three-way selection.
- Fortran intent(out) subroutines return tuples.
- The misspelled Fortran name `calc_respnder_parameters` is preserved for API
  parity.

!=============================================================================
!
! DESCRIPTION OF THE METHOD FOR THE VARIABLE THAT SETS THE MIXTURE FRACTION
! =========================================================================
!
! There are five PDF parameters that need to be calculated for the setting
! variable, which are mu_x_1 (the mean of x in the 1st PDF component), mu_x_2
! (the mean of x in the 2nd PDF component), sigma_x_1 (the standard deviation
! of x in the 1st PDF component), sigma_x_2 (the standard deviation of x in
! the 2nd PDF component), and mixt_frac (the mixture fraction).
!
! This method does NOT work for all values of L_x_1 and L_x_2 (where
! 0 <= L_x_1 <= 1 and 0 <= L_x_2 <= 1).  Only a subregion of this parameter
! space produces valid results.
!
!=============================================================================
!
! DESCRIPTION OF THE METHOD FOR EACH RESPONDING VARIABLE
! ======================================================
!
! In order to find equations for the four PDF parameters for each responding
! variable, which are mu_x_1, mu_x_2, sigma_x_1, and sigma_x_2 (where x stands
! for a responding variable here), four equations are needed.
"""
import jax
import sys
import jax.numpy as jnp


def calc_L_x_Skx_fnc(Skx, sgn_wpxp, small_l_x_1, small_l_x_2):
    """Calculates the values of big_L_x_1 and big_L_x_2 as functions of Skx.

    Returns:
        (big_L_x_1, big_L_x_2)

    Fortran declarations (calc_L_x_Skx_fnc in src/CLUBB_core/new_tsdadg_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          Skx,         & ! Skewness of x (overall)                               [-]
          sgn_wpxp,    & ! Sign of the covariance of w and x (overall)           [-]
          small_l_x_1, & ! Param. for the spread of the 1st PDF comp. mean of x  [-]
          small_l_x_2    ! Param. for the spread of the 2nd PDF comp. mean of x  [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), intent(out) :: &
          big_L_x_1, & ! Parameter for the spread of the 1st PDF comp. mean of x [-]
          big_L_x_2    ! Parameter for the spread of the 2nd PDF comp. mean of x [-]
    """

    factor = jnp.abs(Skx) / jnp.sqrt(4.0 + Skx ** 2)
    cond = Skx * sgn_wpxp >= 0.0
    big_L_x_1 = jnp.where(cond, small_l_x_1, small_l_x_2) * factor
    big_L_x_2 = jnp.where(cond, small_l_x_2, small_l_x_1) * factor
    return big_L_x_1, big_L_x_2


from clubb_jax.src.CLUBB_core.constants_clubb import eps as _EPS  # noqa: E402


# grad-safe sqrt(max(x,0)) — the canonical tracer-toolkit helper.
from clubb_jax.src.CLUBB_core.advance_helper_module import sqrt_clipped


def calc_setter_parameters(xm, xp2, Skx, sgn_wpxp, big_L_x_1, big_L_x_2):
    """Calculates the PDF component means, the PDF component standard deviations,
    and the mixture fraction for the variable that sets the PDF.

    Returns:
        (mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, mixt_frac, coef1, coef2)

    Fortran declarations (calc_setter_parameters in src/CLUBB_core/new_tsdadg_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          xm,        & ! Mean of x (overall)                            [units vary]
          xp2,       & ! Variance of x (overall)                    [(units vary)^2]
          Skx,       & ! Skewness of x                                           [-]
          sgn_wpxp,  & ! Sign of the covariance of w and x (overall)             [-]
          big_L_x_1, & ! Parameter for the spread of the 1st PDF comp. mean of x [-]
          big_L_x_2    ! Parameter for the spread of the 2nd PDF comp. mean of x [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), intent(out) :: &
          mu_x_1,        & ! Mean of x (1st PDF component)              [units vary]
          mu_x_2,        & ! Mean of x (2nd PDF component)              [units vary]
          sigma_x_1_sqd, & ! Variance of x (1st PDF component)      [(units vary)^2]
          sigma_x_2_sqd, & ! Variance of x (2nd PDF component)      [(units vary)^2]
          mixt_frac        ! Mixture fraction                b                   [-]

        real( kind = core_rknd ), intent(out) :: &
          coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
          coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]
    """

    # Calculate the factors in the PDF component mean equations.
    t = Skx * sgn_wpxp / jnp.sqrt(4.0 + Skx ** 2)
    factor_plus = 1.0 + t
    factor_minus = 1.0 - t
    # Calculate the normalized mean of x in the 1st PDF component.
    mu1n = big_L_x_1 * jnp.sqrt(factor_plus / factor_minus) * sgn_wpxp
    # Calculate the normalized mean of x in the 2nd PDF component.
    mu2n = -big_L_x_2 * jnp.sqrt(factor_minus / factor_plus) * sgn_wpxp
    sqrt_xp2 = sqrt_clipped(xp2)
    # Calculate the mean of x in the 1st PDF component.
    mu_x_1 = xm + mu1n * sqrt_xp2
    # Calculate the mean of x in the 2nd PDF component.
    mu_x_2 = xm + mu2n * sqrt_xp2

    # Calculate the mixture fraction.
    a1 = jnp.abs(mu1n) >= _EPS
    a2 = jnp.abs(mu2n) >= _EPS
    mu2n_safe = jnp.where(a2, mu2n, 1.0)
    mf_both = 1.0 / (1.0 + jnp.abs(mu1n / mu2n_safe))
    mf_1 = 1.0 / (1.0 + jnp.abs(mu1n / _EPS))
    mf_2 = 1.0 / (1.0 + jnp.abs(_EPS / mu2n_safe))
    mixt_frac = jnp.where(a1 & a2, mf_both,
                          jnp.where(a1 & (~a2), mf_1,
                                    jnp.where((~a1) & a2, mf_2, 0.5)))

    # Use a minimum magnitude value of mu_x_1_nrmlized in the denominator of a
    # term in the PDF component variance equations in order to prevent a
    # divide-by-zero error.
    mu1n_thresh = jnp.where(mu1n >= 0.0, jnp.maximum(mu1n, _EPS), jnp.minimum(mu1n, -_EPS))
    common = Skx / (3.0 * mixt_frac * mu1n_thresh) - mu1n ** 2 / 3.0 + mu2n ** 2 / 3.0
    base = 1.0 - mixt_frac * mu1n ** 2 - (1.0 - mixt_frac) * mu2n ** 2
    # Calculate the variance of x in the 1st PDF component.
    coef1 = base + (1.0 - mixt_frac) * common
    # Calculate the variance of x in the 2nd PDF component.
    coef2 = base - mixt_frac * common
    sigma_x_1_sqd = coef1 * xp2
    sigma_x_2_sqd = coef2 * xp2
    return mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, mixt_frac, coef1, coef2


def calc_respnder_parameters(xm, xp2, Skx, sgn_wpxp, mixt_frac, big_L_x_1):
    """Calculates the PDF component means, the PDF component standard deviations,
    and the mixture fraction for the variable that sets the PDF.

    The misspelled name matches the Fortran routine.

    Returns:
        (mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, coef1, coef2)

    Fortran declarations (calc_respnder_parameters in src/CLUBB_core/new_tsdadg_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          xm,        & ! Mean of x (overall)                            [units vary]
          xp2,       & ! Variance of x (overall)                    [(units vary)^2]
          Skx,       & ! Skewness of x                                           [-]
          sgn_wpxp,  & ! Sign of the covariance of w and x (overall)             [-]
          mixt_frac, & ! Mixture fraction                                        [-]
          big_L_x_1    ! Parameter for the spread of the 1st PDF comp. mean of x [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), intent(out) :: &
          mu_x_1,        & ! Mean of x (1st PDF component)              [units vary]
          mu_x_2,        & ! Mean of x (2nd PDF component)              [units vary]
          sigma_x_1_sqd, & ! Variance of x (1st PDF component)      [(units vary)^2]
          sigma_x_2_sqd    ! Variance of x (2nd PDF component)      [(units vary)^2]

        real( kind = core_rknd ), intent(out) :: &
          coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
          coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]
    """

    # Calculate the factors in the PDF component mean equations.
    t = Skx * sgn_wpxp / jnp.sqrt(4.0 + Skx ** 2)
    # Calculate the normalized mean of x in the 1st PDF component.
    mu1n = big_L_x_1 * jnp.sqrt((1.0 + t) / (1.0 - t)) * sgn_wpxp
    # Calculate the normalized mean of x in the 2nd PDF component.
    mu2n = -(mixt_frac / (1.0 - mixt_frac)) * mu1n
    sqrt_xp2 = sqrt_clipped(xp2)
    # Calculate the mean of x in the 1st PDF component.
    mu_x_1 = xm + mu1n * sqrt_xp2
    # Calculate the mean of x in the 2nd PDF component.
    mu_x_2 = xm + mu2n * sqrt_xp2

    # Use a minimum magnitude value of mu_x_1_nrmlized in the denominator of a
    # term in the PDF component variance equations in order to prevent a
    # divide-by-zero error.
    mu1n_thresh = jnp.where(mu1n >= 0.0, jnp.maximum(mu1n, _EPS), jnp.minimum(mu1n, -_EPS))
    common = Skx / (3.0 * mixt_frac * mu1n_thresh) - mu1n ** 2 / 3.0 + mu2n ** 2 / 3.0
    base = 1.0 - mixt_frac * mu1n ** 2 - (1.0 - mixt_frac) * mu2n ** 2
    # Calculate the variance of x in the 1st PDF component.
    coef1 = base + (1.0 - mixt_frac) * common
    # Calculate the variance of x in the 2nd PDF component.
    coef2 = base - mixt_frac * common
    sigma_x_1_sqd = coef1 * xp2
    sigma_x_2_sqd = coef2 * xp2
    return mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, coef1, coef2


def tsdadg_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl, wprtp, wpthlp):
    """Selects which variable is used to set the mixture fraction for the PDF
    ("the setter") and which variables are handled after that mixture fraction
    has been set ("the responders").  Traditionally, w has been used to set
    the PDF.  However, here, the variable with the greatest magnitude of
    skewness is used to set the PDF.

    Fortran declarations (tsdadg_pdf_driver in src/CLUBB_core/new_tsdadg_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz, & ! Number of vertical levels
          ngrdcol ! Number of grid columns

        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
          wm,     & ! Mean of w (overall)                [m/s]
          rtm,    & ! Mean of rt (overall)               [kg/kg]
          thlm,   & ! Mean of thl (overall)              [K]
          wp2,    & ! Variance of w (overall)            [m^2/s^2]
          rtp2,   & ! Variance of rt (overall)           [kg^2/kg^2]
          thlp2,  & ! Variance of thl (overall)          [K^2]
          Skw,    & ! Skewness of w (overall)            [-]
          Skrt,   & ! Skewness of rt (overall)           [-]
          Skthl,  & ! Skewness of thl (overall)          [-]
          wprtp,  & ! Covariance of w and rt (overall)   [(m/s)kg/kg]
          wpthlp    ! Covariance of w and thl (overall)  [(m/s)K]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
          mu_w_1,          & ! Mean of w (1st PDF component)        [m/s]
          mu_w_2,          & ! Mean of w (2nd PDF component)        [m/s]
          mu_rt_1,         & ! Mean of rt (1st PDF component)       [kg/kg]
          mu_rt_2,         & ! Mean of rt (2nd PDF component)       [kg/kg]
          mu_thl_1,        & ! Mean of thl (1st PDF component)      [K]
          mu_thl_2,        & ! Mean of thl (2nd PDF component)      [K]
          sigma_w_1_sqd,   & ! Variance of w (1st PDF component)    [m^2/s^2]
          sigma_w_2_sqd,   & ! Variance of w (2nd PDF component)    [m^2/s^2]
          sigma_rt_1_sqd,  & ! Variance of rt (1st PDF component)   [kg^2/kg^2]
          sigma_rt_2_sqd,  & ! Variance of rt (2nd PDF component)   [kg^2/kg^2]
          sigma_thl_1_sqd, & ! Variance of thl (1st PDF component)  [K^2]
          sigma_thl_2_sqd, & ! Variance of thl (2nd PDF component)  [K^2]
          mixt_frac          ! Mixture fraction                     [-]
    """

    # Calculate sgn( <w'rt'> ).
    sgn_wprtp = jnp.where(wprtp >= 0.0, 1.0, -1.0)
    # Calculate sgn( <w'thl'> ).
    sgn_wpthlp = jnp.where(wpthlp >= 0.0, 1.0, -1.0)
    # The sign of the variance of w is always positive.
    sgn_wp2 = jnp.ones_like(Skw)
    l1, l2 = 0.75, 0.5

    bLw1, bLw2 = calc_L_x_Skx_fnc(Skw, sgn_wp2, l1, l2)
    bLrt1, bLrt2 = calc_L_x_Skx_fnc(Skrt, sgn_wprtp, l1, l2)
    bLthl1, bLthl2 = calc_L_x_Skx_fnc(Skthl, sgn_wpthlp, l1, l2)

    aw, art, athl = jnp.abs(Skw), jnp.abs(Skrt), jnp.abs(Skthl)
    # The variable with the greatest magnitude of skewness will be the setter
    # variable and the other variables will be responder variables.
    w_set = (aw >= art) & (aw >= athl)
    rt_set = (art > aw) & (art >= athl)
    thl_set = (~w_set) & (~rt_set)

    # The variable w has the greatest magnitude of skewness.
    sw = calc_setter_parameters(wm, wp2, Skw, sgn_wp2, bLw1, bLw2)
    # The variable rt has the greatest magnitude of skewness.
    srt = calc_setter_parameters(rtm, rtp2, Skrt, sgn_wprtp, bLrt1, bLrt2)
    # The variable thl has the greatest magnitude of skewness.
    sthl = calc_setter_parameters(thlm, thlp2, Skthl, sgn_wpthlp, bLthl1, bLthl2)
    mixt_frac = jnp.where(w_set, sw[4], jnp.where(rt_set, srt[4], sthl[4]))

    rw = calc_respnder_parameters(wm, wp2, Skw, sgn_wp2, mixt_frac, bLw1)
    rrt = calc_respnder_parameters(rtm, rtp2, Skrt, sgn_wprtp, mixt_frac, bLrt1)
    rthl = calc_respnder_parameters(thlm, thlp2, Skthl, sgn_wpthlp, mixt_frac, bLthl1)

    def _sel(mask, setter, responder, i):
        return jnp.where(mask, setter[i], responder[i])

    mu_w_1 = _sel(w_set, sw, rw, 0); mu_w_2 = _sel(w_set, sw, rw, 1)
    sig_w_1 = _sel(w_set, sw, rw, 2); sig_w_2 = _sel(w_set, sw, rw, 3)
    mu_rt_1 = _sel(rt_set, srt, rrt, 0); mu_rt_2 = _sel(rt_set, srt, rrt, 1)
    sig_rt_1 = _sel(rt_set, srt, rrt, 2); sig_rt_2 = _sel(rt_set, srt, rrt, 3)
    mu_thl_1 = _sel(thl_set, sthl, rthl, 0); mu_thl_2 = _sel(thl_set, sthl, rthl, 1)
    sig_thl_1 = _sel(thl_set, sthl, rthl, 2); sig_thl_2 = _sel(thl_set, sthl, rthl, 3)

    # Host adaptation of source per-level warnings, before variance clipping.
    def print_negative_components(w1, w2, rt1, rt2, thl1, thl2):
        import numpy as np
        values = (w1, w2, rt1, rt2, thl1, thl2)
        for index in np.ndindex(w1.shape):
            for name, component, value in zip(("w", "w", "rt", "rt", "thl", "thl"),
                                             ("1st", "2nd", "1st", "2nd", "1st", "2nd"), values):
                if value[index] < 0.0:
                    print(
                        "WARNING:  New TSDADG PDF.  The variance of " + name + " in the " + component + " PDF component is negative and is being clipped to 0.",
                        file=sys.stderr,
                    )
                    number = "1" if component == "1st" else "2"
                    print(
                        "sigma_" + name + "_" + number + "^2 (before clipping) = ",
                        value[index],
                        file=sys.stderr,
                    )
    jax.lax.cond(
        jnp.any((sig_w_1 < 0) | (sig_w_2 < 0) | (sig_rt_1 < 0) | (sig_rt_2 < 0) | (sig_thl_1 < 0) | (sig_thl_2 < 0)),
        lambda _: jax.debug.callback(print_negative_components, sig_w_1, sig_w_2, sig_rt_1, sig_rt_2, sig_thl_1, sig_thl_2, ordered=True),
        lambda _: None, operand=None,
    )

    # Negative variances are clipped to 0.  The Fortran also writes warnings to
    # stderr before clipping.
    clip = lambda s: jnp.where(s < 0.0, 0.0, s)
    return (mu_w_1, mu_w_2, mu_rt_1, mu_rt_2, mu_thl_1, mu_thl_2,
            clip(sig_w_1), clip(sig_w_2), clip(sig_rt_1), clip(sig_rt_2),
            clip(sig_thl_1), clip(sig_thl_2), mixt_frac)
