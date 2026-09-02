"""JAX implementation of ADG1 PDF driver (adg1_adg2_3d_luhar_pdf.F90).

Contains code to calculate the PDF parameters using the Analytic Double
Gaussian 1 (ADG1) PDF closure, the Analytic Double Gaussian 2 (ADG2) PDF
closure, or the 3D Luhar PDF closure.  The ADG1 PDF is the original PDF used
by CLUBB in interactive runs.

Porting deviations:
- Fortran output and inout arguments are returned as tuples or dictionaries.
- The Fortran driver routines are split into smaller pure functions so each
  closure component is directly testable and differentiable.
- ``err_info`` debug/fatal mutation is not threaded through these low-level
  pure functions.  Callers are expected to enforce tolerance guards before
  entering the jitted closure path.
- ``calc_comp_corrs_binormal`` lives in its Fortran home
  ``pdf_utilities.F90`` and is imported/re-exported here because the Fortran
  ADG driver obtains it via ``use pdf_utilities``.
- Comments tied to commented-out or inactive Fortran branches are omitted when
  the corresponding JAX code is absent.

Implements:
  ADG1_w_closure             - mixture fraction and w PDF component params
  ADG1_ADG2_responder_params - PDF component params for rt/thl/u/v/sclr
  ADG1_pdf_driver            - top-level ADG1 PDF parameter driver
  calc_comp_corrs_binormal   - PDF component correlations (pdf_utilities.F90;
                                   _safe_sqrt-hardened grad variant of pdf_utilities.calc_comp_corrs_binormal)
  calc_Luhar_params / close_Luhar_pdf / max_cubic_root / backsolve_Luhar_params
  Luhar_3D_pdf_driver / ADG2_pdf_driver

The PDF moment-integral closures (calc_{wp2xp,wpxp2,wp2xp2,wp4,wpxpyp}_pdf and calc_w_up_in_cloud)
live in pdf_closure_module.py, mirroring their Fortran home pdf_closure_module.F90.
"""

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import rt_tol, thl_tol, w_tol_sqd, zero_threshold, eps


def _safe_sqrt(value):
    return jnp.sqrt(jnp.maximum(value, 0.0))


def ADG1_w_closure(wm, wp2, Skw, sigma_sqd_w, sqrt_wp2, mixt_frac_max_mag):
    """Calculates the mixture fraction, the PDF component means of w, and the PDF
    component variances of w for the Analytic Double Gaussian 1 (ADG1)
    closure.  It sets the widths of both w Gaussians to be the same, and
    furthermore, relates them to the overall variance of w (<w'^2>) by a
    parameter, sigma_sqd_w.  The equation is:

    sigma_w_1^2 = sigma_w_2^2 = sigma_sqd_w * <w'^2>;

    where sigma_w_1^2 is the variance of w in the 1st PDF component,
    sigma_w_2^2 is the variance of w in the 2nd PDF component, and parameter
    sigma_sqd_w must have a value within the range 0 <= sigma_sqd_w <= 1.

    References:
    Golaz, J-C., V. E. Larson, and W. R. Cotton, 2002a: A PDF-based model for
    boundary layer clouds. Part I: Method and model description. J. Atmos.
    Sci., 59, 3540-3551.

    Vincent E. Larson and Jean-Christophe Golaz, 2005: Using Probability
    Density Functions to Derive Consistent Closure Relationships among
    Higher-Order Moments. Mon. Wea. Rev., 133, 1023-1042.

    Args:
        wm:               (ngrdcol, nzt) — mean w [m/s]
        wp2:              (ngrdcol, nzt) — w variance [m^2/s^2]
        Skw:              (ngrdcol, nzt) — w skewness [-]
        sigma_sqd_w:      (ngrdcol, nzt) — PDF width parameter [-]
        sqrt_wp2:         (ngrdcol, nzt) — sqrt(wp2) [m/s]
        mixt_frac_max_mag: scalar — max |mixt_frac - 0.5| + 0.5

    Returns:
        w_1, w_2:           (ngrdcol, nzt) — w means of each PDF component
        w_1_n, w_2_n:       (ngrdcol, nzt) — normalized w means
        varnce_w_1, varnce_w_2: (ngrdcol, nzt) — w variances per component
        mixt_frac:          (ngrdcol, nzt) — mixture fraction (weight of 1st)
    """
    _SKW_TOL = 1.0e-5

    # Width (standard deviation) parameters are non-zero

    # The variable "mixt_frac" is the weight of the 1st PDF component.  The
    # weight of the 2nd PDF component is "1-mixt_frac".  If there isn't any
    # skewness of w (Sk_w = 0 because w'^3 = 0), mixt_frac = 0.5, and both
    # PDF components are equally weighted.  If there is positive skewness of
    # w (Sk_w > 0 because w'^3 > 0), 0 < mixt_frac < 0.5, and the 2nd PDF
    # component has greater weight than does the 1st PDF component.  If there
    # is negative skewness of w (Sk_w < 0 because w'^3 < 0),
    # 0.5 < mixt_frac < 1, and the 1st PDF component has greater weight than
    # does the 2nd PDF component.
    denom_sq = 4.0 * (1.0 - sigma_sqd_w) ** 3 + Skw ** 2
    mf_formula = 0.5 * (1.0 - Skw / jnp.sqrt(denom_sq))
    mixt_frac = jnp.where(jnp.abs(Skw) <= _SKW_TOL, 0.5, mf_formula)

    # Clip mixt_frac, and 1 - mixt_frac, to avoid dividing by a small number.
    # Formula for mixt_frac_max_mag =
    # 1 - ( 1/2 * ( 1 - Skw_max
    #                   / sqrt( 4*( 1 - sigma_sqd_w )^3 + Skw_max^2 ) ) ),
    # where sigma_sqd_w is fixed at 0.4 for this calculation.
    mixt_frac = jnp.clip(mixt_frac, 1.0 - mixt_frac_max_mag, mixt_frac_max_mag)

    # The normalized mean of w for Gaussian "plume" 1 is w_1_n.  It's value
    # will always be greater than 0.  As an example, a value of 1.0 would
    # indicate that the actual mean of w for Gaussian "plume" 1 is found 1.0
    # standard deviation above the overall mean for w.
    one_minus_mf = 1.0 - mixt_frac
    sigma_factor = 1.0 - sigma_sqd_w
    # _safe_sqrt (REFACTOR B5): sigma_factor = 1-sigma_sqd_w -> 0 in well-mixed/surface layers, so the bare
    # sqrt has an inf reverse-mode gradient there (the bomex thlm surface nan). Forward-identical (arg >= 0).
    w_1_n = _safe_sqrt(one_minus_mf / mixt_frac * sigma_factor)
    # The normalized mean of w for Gaussian "plume" 2 is w_2_n.  It's value
    # will always be less than 0.  As an example, a value of -0.5 would
    # indicate that the actual mean of w for Gaussian "plume" 2 is found 0.5
    # standard deviations below the overall mean for w.
    w_2_n = -_safe_sqrt(mixt_frac / one_minus_mf * sigma_factor)

    # The mean of w for Gaussian "plume" 1 is w_1.
    w_1 = wm + sqrt_wp2 * w_1_n
    # The mean of w for Gaussian "plume" 2 is w_2.
    w_2 = wm + sqrt_wp2 * w_2_n

    # The variance of w for Gaussian "plume" 1 for varnce_w_1.
    varnce_w_1 = sigma_sqd_w * wp2
    # The variance of w for Gaussian "plume" 2 for varnce_w_2.
    # The variance in both Gaussian "plumes" is defined to be the same.
    varnce_w_2 = sigma_sqd_w * wp2

    return w_1, w_2, w_1_n, w_2_n, varnce_w_1, varnce_w_2, mixt_frac


def ADG1_ADG2_responder_params(xm, xp2, wp2, sqrt_wp2, wpxp,
                                   w_1_n, w_2_n, mixt_frac,
                                   sigma_sqd_w, beta):
    """Calculates the PDF component means and PDF component variances for thl,
    rt, or sclr when either the ADG1 PDF or ADG2 PDF are used.

    The normalized variance for thl, rt, and sclr for "plume" 1 is:

    { 1 - [1/(1-sigma_sqd_w)]*[ <w'x'>^2 / (<w'^2> * <x'^2>) ] / mixt_frac }
    * { (1/3)*beta + mixt_frac*( 1 - (2/3)*beta ) };

    where "x" stands for thl, rt, or sclr; "mixt_frac" is the weight of
    Gaussian "plume" 1, and 0 <= beta <= 3.

    The factor { (1/3)*beta + mixt_frac*( 1 - (2/3)*beta ) } does not depend
    on which varable "x" stands for.  The factor is multiplied by 2 and
    defined as width_factor_1.

    The factor
    { 1 - [1/(1-sigma_sqd_w)]*[ (w'x')^2 / (w'^2 * x'^2) ] / mixt_frac }
    depends on which variable "x" stands for.  It is multiplied by 1/2 and
    defined as alpha_x, where "x" stands for thl, rt, or sclr.

    Vince Larson added a dimensionless factor so that the
    width of plumes in theta_l, rt can vary.
    beta is a constant defined in CLUBB's tunable parameters.
      Set 0<beta<3.
    beta=1.5_core_rknd recovers Chris Golaz' simplified formula.
    3 Nov 2003

    Args:
        xm:         (ngrdcol, nzt) — mean of x
        xp2:        (ngrdcol, nzt) — variance of x
        wp2:        (ngrdcol, nzt) — variance of w
        sqrt_wp2:   (ngrdcol, nzt) — sqrt(wp2)
        wpxp:       (ngrdcol, nzt) — covariance w'x'
        w_1_n, w_2_n: (ngrdcol, nzt) — normalized w component means
        mixt_frac:  (ngrdcol, nzt) — mixture fraction
        sigma_sqd_w: (ngrdcol, nzt) — PDF width parameter
        beta:       (ngrdcol,) or scalar — CLUBB tunable parameter beta

    Returns:
        x_1, x_2:           (ngrdcol, nzt) — component means
        varnce_x_1, varnce_x_2: (ngrdcol, nzt) — component variances
        alpha_x:            (ngrdcol, nzt) — normalized variance factor
    """
    x_1 = xm - wpxp / (sqrt_wp2 * w_2_n)
    x_2 = xm - wpxp / (sqrt_wp2 * w_1_n)

    alpha_x = 0.5 * (1.0 - wpxp ** 2 / ((1.0 - sigma_sqd_w) * wp2 * xp2))
    alpha_x = jnp.clip(alpha_x, zero_threshold, 1.0)

    # variables for a generalization of Chris Golaz' closure
    # varies width of plumes in theta_l, rt
    if hasattr(beta, 'shape') and len(beta.shape) == 1:
        beta_bc = beta[:, None]  # (ngrdcol, 1)
    else:
        beta_bc = beta
    two_thirds_beta = (2.0 / 3.0) * beta_bc
    width_factor_1 = two_thirds_beta + 2.0 * mixt_frac * (1.0 - two_thirds_beta)

    # Vince Larson multiplied original expressions by width_factor_1,2
    #   to generalize scalar skewnesses.  05 Nov 03
    varnce_x_1 = width_factor_1 * xp2 * alpha_x / mixt_frac
    varnce_x_2 = (2.0 - width_factor_1) * xp2 * alpha_x / (1.0 - mixt_frac)

    return x_1, x_2, varnce_x_1, varnce_x_2, alpha_x


def ADG1_pdf_driver(wm, rtm, thlm, um, vm,
                        wp2, rtp2, thlp2, up2, vp2,
                        Skw, wprtp, wpthlp, upwp, vpwp,
                        sqrt_wp2, sigma_sqd_w,
                        beta, mixt_frac_max_mag):
    """adg1_adg2_3d_luhar_pdf.F90:ADG1_pdf_driver — top-level ADG1 driver.

    All inputs on zt grid (ngrdcol, nzt).

    Args:
        wm:     (ngrdcol, nzt) — mean w
        rtm:    (ngrdcol, nzt) — mean total water
        thlm:   (ngrdcol, nzt) — mean theta-l
        um, vm: (ngrdcol, nzt) — mean winds
        wp2:    (ngrdcol, nzt) — w variance
        rtp2:   (ngrdcol, nzt) — rt variance
        thlp2:  (ngrdcol, nzt) — thl variance
        up2, vp2: (ngrdcol, nzt) — u,v variances
        Skw:    (ngrdcol, nzt) — w skewness
        wprtp, wpthlp, upwp, vpwp: (ngrdcol, nzt) — covariances
        sqrt_wp2: (ngrdcol, nzt) — sqrt(wp2)
        sigma_sqd_w: (ngrdcol, nzt) — PDF width parameter
        beta:   (ngrdcol,) — tunable parameter
        mixt_frac_max_mag: scalar

    Returns:
        dict with keys: w_1, w_2, varnce_w_1, varnce_w_2, mixt_frac,
                        rt_1, rt_2, varnce_rt_1, varnce_rt_2, alpha_rt,
                        thl_1, thl_2, varnce_thl_1, varnce_thl_2, alpha_thl,
                        u_1, u_2, varnce_u_1, varnce_u_2, alpha_u,
                        v_1, v_2, varnce_v_1, varnce_v_2, alpha_v
    """
    # w closure
    (w_1, w_2, w_1_n, w_2_n,
     varnce_w_1, varnce_w_2, mixt_frac) = ADG1_w_closure(
        wm, wp2, Skw, sigma_sqd_w, sqrt_wp2, mixt_frac_max_mag)

    # Responder parameters for rt
    (rt_1, rt_2, varnce_rt_1, varnce_rt_2, alpha_rt) = ADG1_ADG2_responder_params(
        rtm, rtp2, wp2, sqrt_wp2, wprtp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)

    # Responder parameters for thl
    (thl_1, thl_2, varnce_thl_1, varnce_thl_2, alpha_thl) = ADG1_ADG2_responder_params(
        thlm, thlp2, wp2, sqrt_wp2, wpthlp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)

    # Responder parameters for u
    (u_1, u_2, varnce_u_1, varnce_u_2, alpha_u) = ADG1_ADG2_responder_params(
        um, up2, wp2, sqrt_wp2, upwp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)

    # Responder parameters for v
    (v_1, v_2, varnce_v_1, varnce_v_2, alpha_v) = ADG1_ADG2_responder_params(
        vm, vp2, wp2, sqrt_wp2, vpwp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)

    return {
        'w_1': w_1, 'w_2': w_2,
        'w_1_n': w_1_n, 'w_2_n': w_2_n,
        'varnce_w_1': varnce_w_1, 'varnce_w_2': varnce_w_2,
        'mixt_frac': mixt_frac,
        'rt_1': rt_1, 'rt_2': rt_2,
        'varnce_rt_1': varnce_rt_1, 'varnce_rt_2': varnce_rt_2, 'alpha_rt': alpha_rt,
        'thl_1': thl_1, 'thl_2': thl_2,
        'varnce_thl_1': varnce_thl_1, 'varnce_thl_2': varnce_thl_2, 'alpha_thl': alpha_thl,
        'u_1': u_1, 'u_2': u_2,
        'varnce_u_1': varnce_u_1, 'varnce_u_2': varnce_u_2, 'alpha_u': alpha_u,
        'v_1': v_1, 'v_2': v_2,
        'varnce_v_1': varnce_v_1, 'varnce_v_2': varnce_v_2, 'alpha_v': alpha_v,
    }


# calc_comp_corrs_binormal lives in its Fortran home pdf_utilities.F90; imported here
# (the ADG1 PDF driver `use pdf_utilities` for it) and re-exported under its bare Fortran name.
from clubb_jax.src.CLUBB_core.pdf_utilities import (
    calc_comp_corrs_binormal,
)


def calc_Luhar_params(Skx, wpxp, xp2, x_tol_sqd):
    """For the Luhar closure, this subroutine takes Skx (and Skw) as input and
    outputs the mixture fraction, big_m, and small_m. This code was written
    using the equations and nomenclature of Larson et al. (2002) Appendix
    section e.

    The relationship between skewness of x (Skx), mixture fraction (a), and
    Luhar's small m (m) is given by:

    Skx^2 = ( m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3 )
            * ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).

    Luhar's large M (M) is used to more easily express the factor involving
    the m's:

    M = ( m^2 + 1 )^3 / ( m^2 * ( m^2 + 3 )^2 ).

    The equation involving skewness of x becomes:

    Skx^2 = ( 1 / M ) * ( 1 - 2*a )^2 / ( a * ( 1 - a ) );

    or:

    M * Skx^2 = ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).

    This equation can be rewritten as:

    ( a * ( 1 - a ) ) * M * Skx^2 = ( 1 - 2*a )^2;

    as well as:

    ( a - a^2 ) * M * Skx^2 = 1 - 4*a + 4*a^2;

    and eventually as:

    ( 4 + M * Skx^2 ) * a^2 - ( 4 + M * Skx^2 ) * a + 1 = 0.

    Solving the quadratic equation for a:

    a = (1/2) * ( 1 +- Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).

    Since by definition, mu_w_1 >= mu_w_2, a < 0.5 when Skw > 0, the equation
    for mixture fraction is:

    a = (1/2) * ( 1 - Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).

    For 3-D Luhar, the variable (w, rt, or theta-l) with the greatest
    magnitude of skewness is used to calculate mixture fraction.  Since it is
    desirable to still have a < 0.5 when Skw > 0 and a > 0.5 when Skw < 0, the
    sign function is used.  The value of Skx is replaced by:

    Skx|_adj = sgn( <w'x'> ) * Skx;

    where

    sgn( <w'x'> ) = | 1 when <w'x'> >= 0
                    | -1 when <w'x'> < 0.

    Since Skx|_adj^2 = ( sgn( <w'x'> ) * Skx )^2 = ( sgn( <w'x'> ) )^2 * Skx^2
    = Skx^2, the equation for mixture fraction is:

    a = (1/2) * ( 1 - sgn( <w'x'> ) * Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).

    When using the ADG2 closure or when using the 3-D Luhar closure when the
    variable with the greatest magnitude of skewness is w, Skw = Skx and
    sgn( <w'^2> ) is always equal to 1, reducing the equation to its previous
    form.

    References:
    Vincent E. Larson, Jean-Christophe Golaz, and William R. Cotton, 2002:
    Small-Scale and Mesoscale Variability in Cloudy Boundary Layers: Joint
    Probability Density Functions. J. Atmos. Sci., 59, 3519-3539.
    """
    Skx = jnp.asarray(Skx, dtype=jnp.float64)
    wpxp = jnp.asarray(wpxp, dtype=jnp.float64)
    xp2 = jnp.asarray(xp2, dtype=jnp.float64)

    # Width (standard deviation) parameters are non-zero

    # Calculate Luhar's m (small m).
    # If Skx is very small, then small_m will tend to zero which risks
    # divide-by-zero.  To ameliorate this problem, we enforce abs( x_1_n )
    # and abs( x_2_n ) > 0.05.
    # Note:  Luhar's small_m (m) is the only tunable parameter in the Luhar
    #        closure, so this equation can be changed.  However, the value
    #        of m should go toward 0 as Skx goes toward 0 so that the double
    #        Gaussian reduces to a single Gaussian when the distribution is
    #        unskewed.
    # Guarded cube root (cbrt'(0)=inf otherwise; masked by the max floor but would poison grad as 0*inf).
    absS = jnp.abs(Skx)
    cbrt = jnp.where(absS > 0.0, jnp.cbrt(jnp.where(absS > 0.0, absS, 1.0)), 0.0)
    small_m = jnp.maximum((2.0 / 3.0) * cbrt, 0.05)
    # Calculate m^2.
    m_sqd = small_m ** 2
    # Calculate Luhar's M (big M).
    big_m = (1.0 + m_sqd) ** 3 / ((3.0 + m_sqd) ** 2 * m_sqd)

    # Calculate sgn( <w'x'> ).
    sgn = jnp.where(wpxp >= 0.0, 1.0, -1.0)
    # Calculate mixture fraction.
    mixt_frac_vary = 0.5 * (1.0 - sgn * Skx * jnp.sqrt(1.0 / ((4.0 / big_m) + Skx ** 2)))

    vary = xp2 > x_tol_sqd
    mixt_frac = jnp.where(vary, mixt_frac_vary, 0.5)
    big_m = jnp.where(vary, big_m, 0.0)
    small_m = jnp.where(vary, small_m, 0.0)
    return mixt_frac, big_m, small_m


def close_Luhar_pdf(xm, xp2, mixt_frac, small_m, wpxp, x_tol_sqd):
    """For the Luhar closure, this subroutine takes Skx, xm, xp2, and mixt_frac,
    big_m, and small_m (calculated in calc_Luhar_params) as input and outputs
    the PDF component means and variances of a variable x in the joint-PDF
    according to Luhar et al. (1996).  This code was written using the
    equations and nomenclature of Larson et al. (2002) Appendix section e.

    References:
       Vincent E. Larson, Jean-Christophe Golaz, and William R. Cotton, 2002:
       Small-Scale and Mesoscale Variability in Cloudy Boundary Layers: Joint
       Probability Density Functions. J. Atmos. Sci., 59, 3519-3539.
    """
    xm = jnp.asarray(xm, dtype=jnp.float64); xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64); m = jnp.asarray(small_m, dtype=jnp.float64)
    wpxp = jnp.asarray(wpxp, dtype=jnp.float64)
    omf = 1.0 - mf
    msq = 1.0 + m ** 2
    sgn = jnp.where(wpxp >= 0.0, 1.0, -1.0)

    # Width (standard deviation) parameters are non-zero

    # Calculate sgn( <w'x'> ).

    # Calculate the square root of the overall variance of x.
    vary = xp2 > x_tol_sqd
    # Normalized width parameter of x in the 1st PDF component.
    sigma_sqd_x_1 = jnp.where(vary, omf / (mf * msq), 1.0 / msq)
    # Normalized width parameter of x in the 2nd PDF component.
    sigma_sqd_x_2 = jnp.where(vary, mf / (omf * msq), 1.0 / msq)
    # The variance of x in the 1st PDF component.
    varnce_x_1 = jnp.where(vary, sigma_sqd_x_1 * xp2, 0.0)
    # The variance of x in the 2nd PDF component.
    varnce_x_2 = jnp.where(vary, sigma_sqd_x_2 * xp2, 0.0)
    # Normalized mean of x in the 1st PDF component.
    x_1_n = sgn * m * _safe_sqrt(sigma_sqd_x_1)
    # Normalized mean of x in the 2nd PDF component.
    x_2_n = -sgn * m * _safe_sqrt(sigma_sqd_x_2)
    sqrt_xp2 = _safe_sqrt(xp2)
    # The mean of x in the 1st PDF component.
    x_1 = jnp.where(vary, xm + sqrt_xp2 * x_1_n, xm)
    # The mean of x in the 2nd PDF component.
    x_2 = jnp.where(vary, xm + sqrt_xp2 * x_2_n, xm)
    return sigma_sqd_x_1, sigma_sqd_x_2, varnce_x_1, varnce_x_2, x_1_n, x_2_n, x_1, x_2


def max_cubic_root(a_coef, b_coef, c_coef, d_coef):
    """Largest real root of a·x³ + b·x² + c·x + d = 0 (adg1_adg2_3d_luhar_pdf.F90:max_cubic_root). Falls back to
    the quadratic solver when |a| is negligible vs the other coefficients, and to the linear solution when |a|
    and |b| are both negligible. For a true cubic: if the 2nd/3rd roots are real (negligible imaginary part) the
    max over all three real parts, else the always-real first root. Pure-jnp (uses the iter-71-corrected
    cubic_solve) → differentiable away from the root-multiplicity branch points. Inputs broadcast elementwise."""
    from clubb_jax.src.CLUBB_core.calc_roots import cubic_solve, quadratic_solve
    from clubb_jax.src.CLUBB_core.constants_clubb import eps
    a = jnp.asarray(a_coef, dtype=jnp.float64); b = jnp.asarray(b_coef, dtype=jnp.float64)
    c = jnp.asarray(c_coef, dtype=jnp.float64); d = jnp.asarray(d_coef, dtype=jnp.float64)

    a_thresh = 0.001 * jnp.maximum(jnp.maximum(jnp.abs(b), jnp.abs(c)), jnp.abs(d))
    b_thresh = 0.001 * jnp.maximum(jnp.abs(c), jnp.abs(d))

    cr = cubic_solve(a, b, c, d)                      # (..., 3) complex
    r0, r1, r2 = jnp.real(cr[..., 0]), jnp.real(cr[..., 1]), jnp.real(cr[..., 2])
    is_3real = (jnp.abs(jnp.imag(cr[..., 1])) < eps) & (jnp.abs(jnp.imag(cr[..., 2])) < eps)
    cubic_max = jnp.where(is_3real, jnp.maximum(jnp.maximum(r0, r1), r2), r0)

    qr = quadratic_solve(b, c, d)                     # (..., 2) complex
    quad_max = jnp.maximum(jnp.real(qr[..., 0]), jnp.real(qr[..., 1]))

    c_safe = jnp.where(c != 0.0, c, 1.0)
    linear = -d / c_safe

    return jnp.where(jnp.abs(a) > a_thresh, cubic_max,
                     jnp.where(jnp.abs(b) > b_thresh, quad_max, linear))


def backsolve_Luhar_params(Sk_max, Skx, big_m_max, mixt_frac):
    """Backsolve Luhar's m (and M) for a responding variable given the setter's mixt_frac and the variable's
    skewness (adg1_adg2_3d_luhar_pdf.F90:backsolve_Luhar_params, cubic branch — l_use_cubic_backsolve is
    hardcoded true; Sk_max/big_m_max are then unused). When mixt_frac≈0.5 → m=0.05; when |Skx|<eps → m=0;
    otherwise alpha = mf(1−mf)/(1−2mf)²·Skx², and m² is the largest real root of
    (alpha−1)(m²)³ + (3alpha−6)(m²)² + (3alpha−9)(m²) + alpha = 0 (floored at 0.05²). Pure-jnp → differentiable.
    Returns (big_m_x, small_m_x)."""
    from clubb_jax.src.CLUBB_core.constants_clubb import eps
    Skx = jnp.asarray(Skx, dtype=jnp.float64)
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64)

    near_half = jnp.abs(mf - 0.5) < 0.001
    skx_zero = jnp.abs(Skx) < eps

    m_nh = 0.05
    bigm_nh = (1.0 + m_nh ** 2) ** 3 / ((3.0 + m_nh ** 2) ** 2 * m_nh ** 2)

    denom = jnp.where(near_half, 1.0, (1.0 - 2.0 * mf) ** 2)        # guard the 1/(1-2mf)^2 at mf=0.5
    alpha = mf * (1.0 - mf) / denom * Skx ** 2
    alpha_safe = jnp.where(alpha > 0.0, alpha, 1.0)
    bigm_else = 1.0 / alpha_safe
    m2 = max_cubic_root(alpha - 1.0, 3.0 * alpha - 6.0, 3.0 * alpha - 9.0, alpha)
    m_else = jnp.sqrt(jnp.maximum(m2, 0.05 ** 2))

    small_m_x = jnp.where(near_half, m_nh, jnp.where(skx_zero, 0.0, m_else))
    big_m_x = jnp.where(near_half, bigm_nh, jnp.where(skx_zero, jnp.inf, bigm_else))
    return big_m_x, small_m_x


def Luhar_3D_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, Skrt, Skthl, wprtp, wpthlp):
    """3-D Luhar PDF driver (adg1_adg2_3d_luhar_pdf.F90:Luhar_3D_pdf_driver). The variable with the greatest
    |skewness| sets the PDF (calc_Luhar_params → mixt_frac); the other two backsolve their m from the setter's
    mixt_frac (backsolve_Luhar_params); each variable's PDF component means/variances follow from close_Luhar_pdf.
    Pure-jnp (vectorized 3-way setter select) → differentiable. Returns (w_1, w_2, rt_1, rt_2, thl_1, thl_2,
    varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2, varnce_thl_1, varnce_thl_2, mixt_frac)."""
    Skw = jnp.asarray(Skw, dtype=jnp.float64); Skrt = jnp.asarray(Skrt, dtype=jnp.float64)
    Skthl = jnp.asarray(Skthl, dtype=jnp.float64)
    wp2 = jnp.asarray(wp2, dtype=jnp.float64); rtp2 = jnp.asarray(rtp2, dtype=jnp.float64); thlp2 = jnp.asarray(thlp2, dtype=jnp.float64)
    wprtp = jnp.asarray(wprtp, dtype=jnp.float64); wpthlp = jnp.asarray(wpthlp, dtype=jnp.float64)
    w_tsq, rt_tsq, thl_tsq = w_tol_sqd, rt_tol ** 2, thl_tol ** 2

    # calc_Luhar_params for each variable as a potential setter (w uses wp2 as the covariance -> sgn = +1).
    mf_w, bigm_w, m_w = calc_Luhar_params(Skw, wp2, wp2, w_tsq)
    mf_rt, bigm_rt, m_rt = calc_Luhar_params(Skrt, wprtp, rtp2, rt_tsq)
    mf_thl, bigm_thl, m_thl = calc_Luhar_params(Skthl, wpthlp, thlp2, thl_tsq)

    aw, art, athl = jnp.abs(Skw), jnp.abs(Skrt), jnp.abs(Skthl)
    w_set = (aw >= athl) & (aw >= art)
    thl_set = (~w_set) & (athl > aw) & (athl >= art)
    rt_set = (~w_set) & (~thl_set)

    mixt_frac = jnp.where(w_set, mf_w, jnp.where(thl_set, mf_thl, mf_rt))
    setter_Sk = jnp.where(w_set, Skw, jnp.where(thl_set, Skthl, Skrt))
    setter_bigm = jnp.where(w_set, bigm_w, jnp.where(thl_set, bigm_thl, bigm_rt))

    # Each variable: its own m if it is the setter, else backsolve from the setter.
    _, bm_w_resp = backsolve_Luhar_params(setter_Sk, Skw, setter_bigm, mixt_frac)
    _, bm_rt_resp = backsolve_Luhar_params(setter_Sk, Skrt, setter_bigm, mixt_frac)
    _, bm_thl_resp = backsolve_Luhar_params(setter_Sk, Skthl, setter_bigm, mixt_frac)
    m_w_f = jnp.where(w_set, m_w, bm_w_resp)
    m_rt_f = jnp.where(rt_set, m_rt, bm_rt_resp)
    m_thl_f = jnp.where(thl_set, m_thl, bm_thl_resp)

    _, _, vw1, vw2, _, _, w1, w2 = close_Luhar_pdf(wm, wp2, mixt_frac, m_w_f, wp2, w_tsq)
    _, _, vrt1, vrt2, _, _, rt1, rt2 = close_Luhar_pdf(rtm, rtp2, mixt_frac, m_rt_f, wprtp, rt_tsq)
    _, _, vthl1, vthl2, _, _, thl1, thl2 = close_Luhar_pdf(thlm, thlp2, mixt_frac, m_thl_f, wpthlp, thl_tsq)

    return (w1, w2, rt1, rt2, thl1, thl2, vw1, vw2, vrt1, vrt2, vthl1, vthl2, mixt_frac)


def ADG2_pdf_driver(nz, ngrdcol, sclr_dim, sclr_tol,
                    wm, rtm, thlm, wp2, rtp2, thlp2,
                    Skw, wprtp, wpthlp, sqrt_wp2, beta,
                    sclrm, sclrp2, wpsclrp, l_scalar_calc):
    """ADG2 PDF driver (adg1_adg2_3d_luhar_pdf.F90:ADG2_pdf_driver), PDF-parameter outputs for w/rt/thl
    and passive scalars. w's PDF comes from the Luhar closure (calc_Luhar_params + close_Luhar_pdf, with
    sigma_sqd_w overwritten to min(1/(1+m²), 0.99) for ADG1 consistency); rt, thl, and scalar variables are
    responders via ADG1_ADG2_responder_params. Pure-jnp → differentiable.
    Returns (w_1, w_2, rt_1, rt_2, thl_1, thl_2, varnce_w_1, varnce_w_2, varnce_rt_1, varnce_rt_2,
    varnce_thl_1, varnce_thl_2, mixt_frac, alpha_rt, alpha_thl, sigma_sqd_w, sclr_1, sclr_2, varnce_sclr_1,
    varnce_sclr_2, alpha_sclr)."""
    del nz, ngrdcol
    wm = jnp.asarray(wm, dtype=jnp.float64); wp2 = jnp.asarray(wp2, dtype=jnp.float64)
    Skw = jnp.asarray(Skw, dtype=jnp.float64); sqrt_wp2 = jnp.asarray(sqrt_wp2, dtype=jnp.float64)
    w_tsq = w_tol_sqd

    mixt_frac, _, small_m_w = calc_Luhar_params(Skw, wp2, wp2, w_tsq)
    _, _, vw1, vw2, w_1_n, w_2_n, w1, w2 = close_Luhar_pdf(wm, wp2, mixt_frac, small_m_w, wp2, w_tsq)
    sigma_sqd_w = jnp.minimum(1.0 / (1.0 + small_m_w ** 2), 0.99)

    rt1, rt2, vrt1, vrt2, alpha_rt = ADG1_ADG2_responder_params(
        rtm, rtp2, wp2, sqrt_wp2, wprtp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)
    thl1, thl2, vthl1, vthl2, alpha_thl = ADG1_ADG2_responder_params(
        thlm, thlp2, wp2, sqrt_wp2, wpthlp, w_1_n, w_2_n, mixt_frac, sigma_sqd_w, beta)

    if l_scalar_calc:
        sclr1, sclr2, varnce_sclr1, varnce_sclr2, alpha_sclr = jax.vmap(
            lambda sclrm_s, sclrp2_s, wpsclrp_s: ADG1_ADG2_responder_params(
                sclrm_s,
                sclrp2_s,
                wp2,
                sqrt_wp2,
                wpsclrp_s,
                w_1_n,
                w_2_n,
                mixt_frac,
                sigma_sqd_w,
                beta,
            ),
            in_axes=(2, 2, 2),
            out_axes=-1,
        )(sclrm, sclrp2, wpsclrp)
    else:
        sclr1 = jnp.zeros_like(sclrm)
        sclr2 = jnp.zeros_like(sclrm)
        varnce_sclr1 = jnp.zeros_like(sclrm)
        varnce_sclr2 = jnp.zeros_like(sclrm)
        alpha_sclr = jnp.zeros_like(sclrm)

    return (w1, w2, rt1, rt2, thl1, thl2, vw1, vw2, vrt1, vrt2, vthl1, vthl2,
            mixt_frac, alpha_rt, alpha_thl, sigma_sqd_w, sclr1, sclr2,
            varnce_sclr1, varnce_sclr2, alpha_sclr)
