"""JAX port of new_hybrid_pdf.F90 - new-hybrid PDF leaf routines.

Description:
The portion of CLUBB's multivariate, two-component PDF that is the
multivariate, two-component normal PDF of vertical velocity (w), total water
mixing ratio (rt), liquid water potential temperature (thl), and optionally,
the west-east horizontal wind component (u), the south-north horizontal wind
component (v), and passive scalars (sclr).

References:
Griffin and Larson (2020)

Porting deviations:
- The Fortran elemental routines are implemented as array-broadcasting JAX
  functions.
- `calculate_coef_wp4_implicit` and `calculate_mixture_fraction` are thin
  wrappers around the equivalent `new_pdf.py` formulas. The Fortran keeps
  separate implementations; JAX keeps one tested implementation and documents
  the alias.
- Guarded denominators and `_ssqrt` avoid tracing invalid inactive branches.
- Fortran subroutine outputs are returned as tuples.

!=============================================================================
!
! DESCRIPTION OF THE METHOD FOR THE VARIABLE THAT SETS THE MIXTURE FRACTION
! =========================================================================
!
! The variable that sets the mixture fraction for the PDF is w.  There are
! five PDF parameters that need to be calculated, which are mu_w_1 (the mean
! of w is the 1st PDF component), mu_w_2 (the mean of w in the 2nd PDF
! component), sigma_w_1 (the standard deviation of w in the 1st PDF
! component), sigma_w_2 (the standard deviation of w in the 2nd PDF
! component), and mixt_frac (the mixture fraction, which is the weight of the
! 1st PDF component).  In order to solve for these five parameters, five
! equations are needed.  These five equations are the equations for <w>,
! <w'^2>, and <w'^3> as found by integrating over the PDF.  Additionally, two
! more equations, which involve tunable parameters F_w and zeta_w, and which
! are used to help control the spread of the PDF component means and the size
! of the PDF component standard deviations compared to each other,
! respectively, are used in this equation set.
!
! Brian Griffin; September 2017.
!
!=============================================================================
!
! DESCRIPTION OF THE METHOD FOR EACH RESPONDING VARIABLE
! ======================================================
!
! In order to find equations for the four PDF parameters for each responding
! variable, which are mu_x_1, mu_x_2, sigma_x_1, and sigma_x_2 (where x stands
! for a responding variable here), four equations are needed.  These four
! equations are the equations for <x>, <x'^2>, <x'^3>, and <w'x'> as found by
! integrating over the PDF.
!
! where the correlations that are normally found in the <w'x'> equation,
! corr_w_x_1 and corr_w_x_2, have both been set to 0.
!
! Limits on F_w:
!
! The only limits placed on the value of F_w from the w equation set itself
! are 0 <= F_w <= 1.  However, use of the above equation set for responder
! variable x forces an additional limit to be placed on the value of F_w.
! That additional limit restricts the range of F_w to:
!
! <w'x'>^2 / ( <w'^2> * <x'^2> ) <= F_w <= 1.
!
! Brian Griffin; September 2019.
"""

import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.new_pdf import _ssqrt, calc_coef_wp4_implicit, calc_mixture_fraction


def calc_coef_wp2xp_implicit(wp2, mixt_frac, F_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd):
    """Calculate coef_wp2xp_implicit.

    ! Description:
    ! The predictive equation for <w'x'> contains a turbulent advection term of
    ! the form:
    !
    ! - ( 1 / rho_ds ) * d ( rho_ds * <w'^2 x'> ) / dz;
    !
    ! where z is height, rho_ds is the dry, base-state density, and <w'^2 x'> is
    ! calculated by integrating over the PDF.
    !
    ! This equation is of the form:
    !
    ! <w'^2 x'> = coef_wp2xp_implicit * <w'x'>;
    !
    ! In the special case that F_w = 0, <w'x'> must have a value of 0, in which
    ! case mu_x_1 - <x> = mu_x_2 - <x> = 0, and <w'^2 x'> = 0.
    """
    wp2 = jnp.asarray(wp2, dtype=jnp.float64)
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64)
    F = jnp.asarray(F_w, dtype=jnp.float64)
    c1 = jnp.asarray(coef_sigma_w_1_sqd, dtype=jnp.float64)
    c2 = jnp.asarray(coef_sigma_w_2_sqd, dtype=jnp.float64)
    omf = 1.0 - mf
    F_safe = jnp.where(F > 0.0, F, 1.0)
    # Calculate coef_wp2xp_implicit.
    coef = (_ssqrt(mf * omf)
            * (F * (omf / mf - mf / omf) + c1 - c2)
            * _ssqrt(wp2 / F_safe))
    return jnp.where(F > 0.0, coef, 0.0)


# --- new_hybrid_pdf.F90 aliases (Griffin & Larson 2020) ---
# These are byte-identical / sgn=+1 specializations of the new_pdf.F90 forms; CLUBB defines them
# separately in module new_hybrid_pdf with their own (capital-C) names. Kept as thin wrappers so callers of
# either module map to one tested implementation.

def calculate_coef_wp4_implicit(mixt_frac, F_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd):
    """Byte-identical to calc_coef_wp4_implicit.

    Porting deviation: Fortran repeats the full `calculate_coef_wp4_implicit`
    body in this module. JAX delegates to the shared `new_pdf.py` formula to
    avoid two copies of the same algebra.
    """
    return calc_coef_wp4_implicit(mixt_frac, F_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd)


def calculate_mixture_fraction(Skw, F_w, zeta_w):
    """Calculates mixture fraction.

    ! References:
    ! Griffin and Larson (2020)

    Porting deviation: Fortran keeps this as a separate elemental function;
    JAX calls the shared new_pdf mixture-fraction helper with sgn_wpxp = +1
    because w is the setting variable.
    """
    return calc_mixture_fraction(Skw, F_w, zeta_w, 1.0)


def calculate_w_params(wm, wp2, Skw, F_w, zeta_w):
    """Calculates the PDF component means, the PDF component standard deviations,
    and the mixture fraction for the variable that sets the PDF.

    ! References:
    ! Griffin and Larson (2020)
    """
    wm = jnp.asarray(wm, dtype=jnp.float64); wp2 = jnp.asarray(wp2, dtype=jnp.float64)
    Skw = jnp.asarray(Skw, dtype=jnp.float64); F_w = jnp.asarray(F_w, dtype=jnp.float64)
    zeta = jnp.asarray(zeta_w, dtype=jnp.float64)

    # Calculate the mixture fraction.
    mixt_frac = calculate_mixture_fraction(Skw, F_w, zeta)
    valid = (mixt_frac > 0.0) & (mixt_frac < 1.0)
    mf = jnp.where(valid, mixt_frac, 0.5)          # avoid div-by-zero in the invalid branch
    omf = 1.0 - mf

    # Calculate the mean of w in the 1st PDF component.
    mu_w_1_v = wm + _ssqrt(F_w * (omf / mf) * wp2)
    # Calculate the mean of w in the 2nd PDF component.
    mu_w_2_v = wm - (mf / omf) * (mu_w_1_v - wm)
    # Calculate the standard deviation of w in the 1st PDF component.
    # sigma_w_1 = sqrt( ( ( zeta_w + 1 ) * ( 1 - F_w ) )
    #                   / ( ( zeta_w + 2 ) * mixt_frac ) * <w'^2> )
    c1_v = ((zeta + 1.0) * (1.0 - F_w)) / ((zeta + 2.0) * mf)
    # Calculate the standard deviation of w in the 2nd PDF component.
    # sigma_w_2 = sqrt( ( mixt_frac * sigma_w_1^2 )
    #                   / ( ( 1 - mixt_frac ) * ( 1 + zeta_w ) ) )
    #           = sqrt( ( 1 - F_w )
    #                   / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ) * <w'^2> )
    c2_v = (1.0 - F_w) / ((zeta + 2.0) * omf)
    sig1_v = _ssqrt(c1_v * wp2); sig2_v = _ssqrt(c2_v * wp2)

    # The mixture fraction produced is invalid.  This should only happen in
    # the scenario where F_w = 0 and | Skw | > 0, where the value of
    # mixt_frac has been set to -1.  Set all output variables to 0 in this
    # scenario.  Since F_w is a function of skewness, the mixture fraction
    # and the PDF should always be valid, and this section of code shouldn't
    # be entered.
    z = jnp.zeros_like(mu_w_1_v)
    mu_w_1 = jnp.where(valid, mu_w_1_v, z)
    mu_w_2 = jnp.where(valid, mu_w_2_v, z)
    sigma_w_1 = jnp.where(valid, sig1_v, z)
    sigma_w_2 = jnp.where(valid, sig2_v, z)
    coef_sigma_w_1_sqd = jnp.where(valid, c1_v, z)
    coef_sigma_w_2_sqd = jnp.where(valid, c2_v, z)
    return (mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, mixt_frac,
            coef_sigma_w_1_sqd, coef_sigma_w_2_sqd)


def calculate_responder_params(xm, xp2, Skx, wpxp, wp2, F_w, mixt_frac):
    """Calculates the PDF component means and the PDF component standard
    deviations for a responding variable (a variable that is not used to set
    the mixture fraction).

    ! References:
    ! Griffin and Larson (2020)
    """
    xm = jnp.asarray(xm, dtype=jnp.float64); xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    Skx = jnp.asarray(Skx, dtype=jnp.float64); wpxp = jnp.asarray(wpxp, dtype=jnp.float64)
    wp2 = jnp.asarray(wp2, dtype=jnp.float64); F_w = jnp.asarray(F_w, dtype=jnp.float64)
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64)
    omf = 1.0 - mf

    gate = jnp.abs(wpxp) > 0.0
    # Note:  when |<w'x'>| > 0, F_w, <w'^2>, and <x'^2> must all have values
    #        greater than 0.
    wpxp_safe = jnp.where(gate, wpxp, 1.0)              # guard the /(3·wpxp) and /√(F_w·wp2) terms
    fw_wp2_safe = jnp.where(gate, F_w * wp2, 1.0)       # >0 whenever gate is true (Fortran note)
    den_safe = jnp.where(gate, 3.0 * F_w * wp2 * xp2, 1.0)

    # Calculate the mean of x in the 1st PDF component.
    mu_x_1_v = xm + _ssqrt(omf / mf) * wpxp / _ssqrt(fw_wp2_safe)
    # Calculate the mean of x in the 2nd PDF component.
    mu_x_2_v = xm - _ssqrt(mf / omf) * wpxp / _ssqrt(fw_wp2_safe)
    # Calculate the variance of x in the 1st PDF component.
    # sigma_x_1^2
    # = ( 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
    #         * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
    #     - ( ( 1 + mixt_frac ) / mixt_frac )
    #       * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ) )
    #   * <x'^2>
    c1_v = (1.0 + _ssqrt(omf / mf) * Skx * _ssqrt(fw_wp2_safe * xp2) / (3.0 * wpxp_safe)
            - ((1.0 + mf) / mf) * wpxp ** 2 / den_safe)
    # Mathematically, the value of coef_sigma_x_1_sqd cannot be less than 0.
    # Numerically, this can happen when numerical round off error causes an
    # epsilon-sized negative value.  When this happens, reset the value of
    # coef_sigma_x_1_sqd to 0.
    c1_v = jnp.maximum(c1_v, 0.0)
    # Calculate the variance of x in the 2nd PDF component.
    # sigma_x_2^2
    # = ( 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
    #         * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
    #     + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
    #       * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ) )
    #   * <x'^2>
    c2_v = (1.0 - _ssqrt(mf / omf) * Skx * _ssqrt(fw_wp2_safe * xp2) / (3.0 * wpxp_safe)
            + ((mf - 2.0) / omf) * wpxp ** 2 / den_safe)
    # Mathematically, the value of coef_sigma_x_2_sqd cannot be less than 0.
    # Numerically, this can happen when numerical round off error causes an
    # epsilon-sized negative value.  When this happens, reset the value of
    # coef_sigma_x_2_sqd to 0.
    c2_v = jnp.maximum(c2_v, 0.0)

    # When <w'x'> has a value of 0, the PDF becomes a single Gaussian.  This
    # only works when Skx = 0.  However, when Skx /= 0, the value of min_F_x
    # is greater than 0, preventing a problem where F_x = 0 but | Skx | > 0.
    mu_x_1 = jnp.where(gate, mu_x_1_v, xm)
    mu_x_2 = jnp.where(gate, mu_x_2_v, xm)
    coef_sigma_x_1_sqd = jnp.where(gate, c1_v, 1.0)
    coef_sigma_x_2_sqd = jnp.where(gate, c2_v, 1.0)
    sigma_x_1_sqd = coef_sigma_x_1_sqd * xp2
    sigma_x_2_sqd = coef_sigma_x_2_sqd * xp2
    return (mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd,
            coef_sigma_x_1_sqd, coef_sigma_x_2_sqd)


def calc_coefs_wpxp2_semiimpl(wp2, wpxp, mixt_frac, F_w, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd):
    """Calculate coef_wpxp2_implicit and term_wpxp2_explicit.

    ! Description:
    ! The predictive equation for <x'^2> contains a turbulent advection term of
    ! the form:
    !
    ! - ( 1 / rho_ds ) * d ( rho_ds * <w'x'^2> ) / dz;
    !
    ! where z is height, rho_ds is the dry, base-state density, and <w'x'^2> is
    ! calculated by integrating over the PDF.
    !
    ! <w'x'^2> = coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit;
    """
    wp2 = jnp.asarray(wp2, dtype=jnp.float64); wpxp = jnp.asarray(wpxp, dtype=jnp.float64)
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64); F_w = jnp.asarray(F_w, dtype=jnp.float64)
    cx1 = jnp.asarray(coef_sigma_x_1_sqd, dtype=jnp.float64); cx2 = jnp.asarray(coef_sigma_x_2_sqd, dtype=jnp.float64)
    omf = 1.0 - mf

    # Calculate coef_wpxp2_implicit and term_wpxp2_explicit.
    vary = (F_w > 0.0) & (wp2 > 0.0)
    sFwwp2 = _ssqrt(F_w * wp2)
    sFwwp2_safe = jnp.where(sFwwp2 > 0.0, sFwwp2, 1.0)
    base = _ssqrt(mf * omf)
    coef = base * sFwwp2 * (cx1 - cx2)
    term = base * wpxp ** 2 / sFwwp2_safe * (omf / mf - mf / omf)
    return jnp.where(vary, coef, 0.0), jnp.where(vary, term, 0.0)


def calc_coefs_wpxpyp_semiimpl(wp2, wpxp, wpyp, mixt_frac, F_w,
                               coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
                               coef_sigma_y_1_sqd, coef_sigma_y_2_sqd):
    """Calculate coef_wpxpyp_implicit and term_wpxpyp_explicit.

    ! Description:
    ! The predictive equation for <w'x'y'> contains a turbulent advection term
    ! of the form:
    !
    ! - ( 1 / rho_ds ) * d ( rho_ds * <w'x'y'> ) / dz;
    !
    ! where z is height, rho_ds is the dry, base-state density, and <w'x'y'> is
    ! calculated by integrating over the PDF.
    !
    ! <w'x'y'> = coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit;
    !
    ! There are also special cases for the above equations.
    """
    wp2 = jnp.asarray(wp2, dtype=jnp.float64); wpxp = jnp.asarray(wpxp, dtype=jnp.float64)
    wpyp = jnp.asarray(wpyp, dtype=jnp.float64)
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64); F_w = jnp.asarray(F_w, dtype=jnp.float64)
    cx1 = jnp.asarray(coef_sigma_x_1_sqd, dtype=jnp.float64); cx2 = jnp.asarray(coef_sigma_x_2_sqd, dtype=jnp.float64)
    cy1 = jnp.asarray(coef_sigma_y_1_sqd, dtype=jnp.float64); cy2 = jnp.asarray(coef_sigma_y_2_sqd, dtype=jnp.float64)
    omf = 1.0 - mf

    # Calculate coef_wpxpyp_implicit and term_wpxpyp_explicit.
    s1 = _ssqrt(cx1 * cy1); s2 = _ssqrt(cx2 * cy2)
    xy_vary = ((cx1 * cy1 > 0.0) | (cx2 * cy2 > 0.0)) & (F_w > 0.0) & (wp2 > 0.0)
    denom = mf * s1 + omf * s2
    denom_safe = jnp.where(denom > 0.0, denom, 1.0)
    # coefs_factor_xy
    # = ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    #     - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    #   / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    #       + ( 1 - mixt_frac )
    #         * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    f_xy = (s1 - s2) / denom_safe

    sFwwp2 = _ssqrt(F_w * wp2)
    sFwwp2_safe = jnp.where(sFwwp2 > 0.0, sFwwp2, 1.0)
    base = _ssqrt(mf * omf)

    coef_vary = base * sFwwp2 * f_xy
    term_vary = base * wpxp * wpyp / sFwwp2_safe * (omf / mf - mf / omf - f_xy)
    coef_else = jnp.where(F_w > 0.0, base * sFwwp2 * (omf / mf - mf / omf), 0.0)

    coef = jnp.where(xy_vary, coef_vary, coef_else)
    term = jnp.where(xy_vary, term_vary, 0.0)
    return coef, term


__all__ = [
    "calc_coef_wp2xp_implicit",
    "calculate_coef_wp4_implicit",
    "calculate_mixture_fraction",
    "calculate_w_params",
    "calculate_responder_params",
    "calc_coefs_wpxp2_semiimpl",
    "calc_coefs_wpxpyp_semiimpl",
]
