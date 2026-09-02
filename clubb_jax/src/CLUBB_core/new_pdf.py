"""JAX port of new_pdf.F90 - the "new hybrid" PDF closure helpers.

Description:
The portion of CLUBB's multivariate, two-component PDF that is the
trivariate, two-component normal PDF of vertical velocity (w), total water
mixing ratio (rt), and liquid water potential temperature (thl).

References:
Griffin and Larson (2018)

Porting deviations:
- The Fortran nz shape argument is omitted; JAX broadcasts over the shapes
  supplied by callers.
- Fortran subroutines with intent(out) arrays return tuples instead.
- Fortran WHERE blocks are represented with jnp.where and guarded
  denominators/square roots so the routines remain JIT friendly.
- calc_mixture_fraction returns the symmetric limit for F_x = 0; the Fortran
  writes to stderr and error-stops when F_x = 0 and | Skx | > 0. That
  impossible PDF state is guarded before these routines are used by the core.
- sort_roots uses jnp.sort instead of the explicit Fortran three-root sorting
  network; it returns the same values.
- _ssqrt is imported from pdf_utilities so sibling JAX PDF modules can share a
  grad-safe sqrt helper.

!=============================================================================
!
! DESCRIPTION OF THE METHOD FOR THE VARIABLE THAT SETS THE MIXTURE FRACTION
! =========================================================================
!
! Many times, w has been used as the variable that sets the mixture fraction
! for the PDF.  There are five PDF parameters that need to be calculated,
! which are mu_w_1 (the mean of w is the 1st PDF component), mu_w_2 (the mean
! of w in the 2nd PDF component), sigma_w_1 (the standard deviation of w in
! the 1st PDF component), sigma_w_2 (the standard deviation of w in the 2nd
! PDF component), and mixt_frac (the mixture fraction, which is the weight of
! the 1st PDF component).  In order to solve for these five parameters, five
! equations are needed.  These five equations are the equations for <w>,
! <w'^2>, and <w'^3> as found by integrating over the PDF.  Additionally, two
! more equations, which involve tunable parameters F_w and zeta_w, and which
! are used to help control the spread of the PDF component means and the size
! of the PDF component standard deviations compared to each other,
! respectively, are used in this equation set.  The five equations are:
!
! <w> = mixt_frac * mu_w_1 + ( 1 - mixt_frac ) * mu_w_2;
!
! <w'^2> = mixt_frac * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
!          + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 );
!
! <w'^3> = mixt_frac * ( mu_w_1 - <w> )
!                    * ( ( mu_w_1 - <w> )^2 + 3 * sigma_w_1^2 )
!          + ( 1 - mixt_frac ) * ( mu_w_2 - <w> )
!                              * ( ( mu_w_2 - <w> )^2 + 3 * sigma_w_2^2 );
!
! mu_w_1 - <w> = sqrt(F_w) * ( sqrt( 1 - mixt_frac ) / sqrt( mixt_frac ) )
!                * sqrt( <w'^2> );
!
! where 0 <= F_w <= 1; and
!
! 1 + zeta_w = ( mixt_frac * sigma_w_1^2 )
!              / ( ( 1 - mixt_frac ) * sigma_w_2^2 );
!
! where zeta_w > -1.
!
! Following convention for w, mu_w_1 is defined to be greater than or equal to
! mu_w_2 (and is also greater than or equal to <w>, while mu_w_2 is less than
! or equal to <w>).  This is relationship is found in the mu_w_1 - <w>
! equation above.
!
! The resulting equations for the five PDF parameters are:
!
! mixt_frac
! = ( 4 * F_w^3
!     + 18 * F_w * ( zeta_w + 1 ) * ( 1 - F_w ) / ( zeta_w + 2 )
!     + 6 * F_w^2 * ( 1 - F_w ) / ( zeta_w + 2 )
!     + Skw^2
!     - Skw * sqrt( 4 * F_w^3
!                   + 12 * F_w^2 * ( 1 - F_w )
!                   + 36 * F_w * ( zeta_w + 1 ) * ( 1 - F_w )^2
!                     / ( zeta_w + 2 )^2
!                   + Skw^2 ) )
!   / ( 2 * F_w * ( F_w - 3 )^2 + 2 * Skw^2 );
!
! mu_w_1 = <w> + sqrt( F_w * ( ( 1 - mixt_frac ) / mixt_frac ) * <w'^2> );
!
! mu_w_2 = <w> - ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_w_1 - <w> );
!
! sigma_w_1
! = sqrt( ( ( zeta_w + 1 ) * ( 1 - F_w ) )
!         / ( ( zeta_w + 2 ) * mixt_frac ) * <w'^2> ); and
!
! sigma_w_2
! = sqrt( ( mixt_frac * sigma_w_1^2 )
!         / ( ( 1 - mixt_frac ) * ( 1 + zeta_w ) ) );
!
! where Skw is the skewness of w, and Skw = <w'^3> / <w'^2>^(3/2).
!
! This method works for all values of F_w (where 0 <= F_w <= 1) and zeta_w
! (where zeta_w > -1).
!
!
! Generalized equations for any variable, x, that sets the mixture fraction:
!
! A slight alteration is made to the above equations in order to have any
! variable, x, set the mixture fraction.  The same five PDF parameters need to
! be calculated for the setting variable, which are mu_x_1 (the mean of x in
! the 1st PDF component), mu_x_2 (the mean of x in the 2nd PDF component),
! sigma_x_1 (the standard deviation of x in the 1st PDF component), sigma_x_2
! (the standard deviation of x in the 2nd PDF component), and mixt_frac (the
! mixture fraction).  Again, five equations are needed, and they are the
! equations for <x>, <x'^2>, and <x'^3> as found by integrating over the PDF,
! as well as the equations that involve tunable parameters F_x and zeta_x.
! However, the equation for F_x is multiplied by a new variable,
! sgn( <w'x'> ), where <w'x'> is the covariance of w and x, and sgn( <w'x'> )
! is given by:
!
! sgn( <w'x'> ) = |  1, when <w'x'> >= 0;
!                 | -1, when <w'x'> < 0.
!
! The five equations are:
!
! <x> = mixt_frac * mu_x_1 + ( 1 - mixt_frac ) * mu_x_2;
!
! <x'^2> = mixt_frac * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
!          + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 );
!
! <x'^3> = mixt_frac * ( mu_x_1 - <x> )
!                    * ( ( mu_x_1 - <x> )^2 + 3 * sigma_x_1^2 )
!          + ( 1 - mixt_frac ) * ( mu_x_2 - <x> )
!                              * ( ( mu_x_2 - <x> )^2 + 3 * sigma_x_2^2 );
!
! mu_x_1 - <x> = sqrt(F_x) * ( sqrt( 1 - mixt_frac ) / sqrt( mixt_frac ) )
!                * sqrt( <x'^2> ) * sgn( <w'x'> );
!
! where 0 <= F_x <= 1; and
!
! 1 + zeta_x = ( mixt_frac * sigma_x_1^2 )
!              / ( ( 1 - mixt_frac ) * sigma_x_2^2 );
!
! where zeta_x > -1.
!
! The only equations that are altered are the equation for mu_x_1 and the
! equation for mixt_frac, which now both contain a sgn( <w'x'> ).  The mu_x_2
! equation is not altered, but the sign of mu_x_2 - <x> will be the opposite
! of the sign of mu_x_1 - <x>.  The resulting equations for the five PDF
! parameters are:
!
! mixt_frac
! = ( 4 * F_x^3
!     + 18 * F_x * ( zeta_x + 1 ) * ( 1 - F_x ) / ( zeta_x + 2 )
!     + 6 * F_x^2 * ( 1 - F_x ) / ( zeta_x + 2 )
!     + Skx^2
!     - Skx * sgn( <w'x'> )
!           * sqrt( 4 * F_x^3
!                   + 12 * F_x^2 * ( 1 - F_x )
!                   + 36 * F_x * ( zeta_x + 1 ) * ( 1 - F_x )^2
!                     / ( zeta_x + 2 )^2
!                   + Skx^2 ) )
!   / ( 2 * F_x * ( F_x - 3 )^2 + 2 * Skx^2 );
!
! mu_x_1 = <x> + sqrt( F_x * ( ( 1 - mixt_frac ) / mixt_frac ) * <x'^2> )
!                * sgn( <w'x'> );
!
! mu_x_2 = <x> - ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_x_1 - <x> );
!
! sigma_x_1
! = sqrt( ( ( zeta_x + 1 ) * ( 1 - F_x ) )
!         / ( ( zeta_x + 2 ) * mixt_frac ) * <x'^2> ); and
!
! sigma_x_2
! = sqrt( ( mixt_frac * sigma_x_1^2 )
!         / ( ( 1 - mixt_frac ) * ( 1 + zeta_x ) ) );
!
! where Skx is the skewness of x, and Skx = <x'^3> / <x'^2>^(3/2).
!
! This method works for all values of F_x (where 0 <= F_x <= 1) and zeta_x
! (where zeta_x > -1).
!
! When the generalized form is solved for w (x = w), sgn( <w'^2> ) = 1, and
! the equations are unaltered from the equations listed above for w.
!
!
! Special case:
!
! When Skx = 0 and F_x = 0, the equation for mixt_frac is undefined.  The
! equation for mixture fraction in this scenario can be derived by using the
! above equation for mixture fraction and then setting Skx = 0.  The resulting
! equation becomes:
!
! mixt_frac
! = ( 4 * F_x^3
!     + 18 * F_x * ( zeta_x + 1 ) * ( 1 - F_x ) / ( zeta_x + 2 )
!     + 6 * F_x^2 * ( 1 - F_x ) / ( zeta_x + 2 ) )
!   / ( 2 * F_x * ( F_x - 3 )^2 ).
!
! All of the terms in the numerator and denominator contain a F_x, so this
! equation can be rewritten as:
!
! mixt_frac
! = ( 4 * F_x^2
!     + 18 * ( zeta_x + 1 ) * ( 1 - F_x ) / ( zeta_x + 2 )
!     + 6 * F_x * ( 1 - F_x ) / ( zeta_x + 2 ) )
!   / ( 2 * ( F_x - 3 )^2 ).
!
! Now setting F_x = 0, the equation becomes:
!
! mixt_frac = ( 18 * ( zeta_x + 1 ) / ( zeta_x + 2 ) ) / 18;
!
! which can be rewritten as:
!
! mixt_frac = ( zeta_x + 1 ) / ( zeta_x + 2 ).
!
! When F_x = 0, Skx must have a value of 0 in order for the PDF to function
! correctly.  When F_x = 0, mu_x_1 = mu_x_2.  When the two PDF component means
! are equal to each other (and to the overall mean, <x>), the only value of
! Skx that can be represented is a value of 0.  In the equation for mixture
! fraction, when F_x is set to 0, but | Skx | > 0, mixt_frac will either have
! a value of 0 or 1, depending on whether Skx is positive or negative,
! respectively.
!
! The value of F_x should be set as a function of Skx.  The value F_x should
! go toward 0 as | Skx | (or Skx^2) goes toward 0.  The value of F_x should
! go toward 1 as | Skx | (or Skx^2) goes to infinity.
!
!
! Tunable parameters:
!
! 1) F_x:  This parameter controls the spread of the PDF component means.  The
!          range of this parameter is 0 <= F_x <= 1.  When F_x = 0, the two
!          PDF component means (mu_x_1 and mu_x_2) are equal to each other
!          (and Skx must equal 0).  All of the variance of x is accounted for
!          by the PDF component standard deviations (sigma_x_1 and sigma_x_2).
!          When F_x = 1, mu_x_1 and mu_x_2 are spread as far apart as they can
!          be.  Both PDF component standard deviations (sigma_x_1 and
!          sigma_x_2) are equal to 0, and all of the variance of x is
!          accounted for by the spread of the PDF component means.
!
! 2) zeta_x:  This parameter controls the size of the PDF component standard
!             deviations compared to each other.
!
! Notes:
!
! When F_x = 0 (which can only happen when Skx = 0), mu_x_1 = mu_x_2, and
! mixt_frac = ( zeta_x + 1 ) / ( zeta_x + 2 ).  When these equations are
! substituted into the equations for sigma_x_1 and sigma_x_2, the result is
! sigma_x_1 = sigma_x_2 = sqrt( <x'^2> ).  This means that the distribution
! becomes a single Gaussian when F_x = 0 (and Skx = 0).  This happens
! regardless of the value of zeta_x.
!
! Brian Griffin; September 2017.
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()


def calc_coef_wp4_implicit(mixt_frac, F_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd):
    """Calculate coef_wp4_implicit.

    ! Description:
    ! The predictive equation for <w'^3> contains a turbulent advection term of
    ! the form:
    !
    ! - ( 1 / rho_ds ) * d ( rho_ds * <w'^4> ) / dz;
    !
    ! where z is height, rho_ds is the dry, base-state density, and <w'^4> is
    ! calculated by integrating over the PDF.  The equation for <w'^4> is:
    !
    ! <w'^4> = mixt_frac * ( 3 * sigma_w_1^4
    !                        + 6 * ( mu_w_1 - <w> )^2 * sigma_w_1^2
    !                        + ( mu_w_1 - <w> )^4 )
    !          + ( 1 - mixt_frac ) * ( 3 * sigma_w_2^4
    !                                  + 6 * ( mu_w_2 - <w> )^2 * sigma_w_2^2
    !                                  + ( mu_w_2 - <w> )^4 ).
    !
    ! The following substitutions are made into the above equation:
    !
    ! mu_w_1 - <w> = sqrt(F_w) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
    !                * sqrt( <w'^2> );
    !
    ! mu_w_2 - <w> = - sqrt(F_w) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
    !                  * sqrt( <w'^2> );
    !
    ! sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> ); and
    !
    ! sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> ).
    !
    ! The equation for <w'4> becomes:
    !
    ! <w'^4> = ( 3 * mixt_frac * coef_sigma_w_1_sqd^2
    !            + 6 * F_w * ( 1 - mixt_frac ) * coef_sigma_w_1_sqd
    !            + F_w^2 * ( 1 - mixt_frac )^2 / mixt_frac
    !            + 3 * ( 1 - mixt_frac ) * coef_sigma_w_2_sqd^2
    !            + 6 * F_w * mixt_frac * coef_sigma_w_2_sqd
    !            + F_w^2 * mixt_frac^2 / ( 1 - mixt_frac ) ) * <w'^2>^2.
    !
    ! This equation is of the form:
    !
    ! <w'^4> = coef_wp4_implicit * <w'^2>^2;
    !
    ! where:
    !
    ! coef_wp4_implicit = 3 * mixt_frac * coef_sigma_w_1_sqd^2
    !                     + 6 * F_w * ( 1 - mixt_frac ) * coef_sigma_w_1_sqd
    !                     + F_w^2 * ( 1 - mixt_frac )^2 / mixt_frac
    !                     + 3 * ( 1 - mixt_frac ) * coef_sigma_w_2_sqd^2
    !                     + 6 * F_w * mixt_frac * coef_sigma_w_2_sqd
    !                     + F_w^2 * mixt_frac^2 / ( 1 - mixt_frac ).
    !
    ! While the <w'^4> term is found in the <w'^3> predictive equation and not
    ! the <w'^2> predictive equation, the <w'^3> and <w'^2> predictive equations
    ! are solved together.  This allows the term containing <w'^4> to be solved
    ! implicitly.
    """
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64)
    F = jnp.asarray(F_w, dtype=jnp.float64)
    c1 = jnp.asarray(coef_sigma_w_1_sqd, dtype=jnp.float64)
    c2 = jnp.asarray(coef_sigma_w_2_sqd, dtype=jnp.float64)
    omf = 1.0 - mf
    # Calculate coef_wp4_implicit.
    return (3.0 * mf * c1 ** 2
            + 6.0 * F * omf * c1
            + F ** 2 * omf ** 2 / mf
            + 3.0 * omf * c2 ** 2
            + 6.0 * F * mf * c2
            + F ** 2 * mf ** 2 / omf)


def calc_mixture_fraction(Skx, F_x, zeta_x, sgn_wpxp):
    """Calculates mixture fraction.

    ! References:
    ! Griffin and Larson (2018)
    """
    Skx = jnp.asarray(Skx, dtype=jnp.float64)
    F = jnp.asarray(F_x, dtype=jnp.float64)
    zeta = jnp.asarray(zeta_x, dtype=jnp.float64)
    sgn = jnp.asarray(sgn_wpxp, dtype=jnp.float64)
    zp2 = zeta + 2.0

    sqrt_arg = (4.0 * F ** 3 + 12.0 * F ** 2 * (1.0 - F)
                + 36.0 * F * (zeta + 1.0) * (1.0 - F) ** 2 / zp2 ** 2 + Skx ** 2)
    num = (4.0 * F ** 3
           + 18.0 * F * (zeta + 1.0) * (1.0 - F) / zp2
           + 6.0 * F ** 2 * (1.0 - F) / zp2
           + Skx ** 2
           - Skx * sgn * jnp.sqrt(jnp.maximum(sqrt_arg, 0.0)))
    den = 2.0 * F * (F - 3.0) ** 2 + 2.0 * Skx ** 2
    # Calculate mixture fraction, which is the weight of the 1st PDF component.
    # The 2nd PDF component has a weight of 1 - mixt_frac.
    mf_pos = num / jnp.where(den > 0.0, den, 1.0)

    # F_x = 0 and Skx = 0
    mf_symmetric = (zeta + 1.0) / zp2
    return jnp.where(F > 0.0, mf_pos, mf_symmetric)


# grad-safe sqrt(max(x,0)) — the canonical tracer-toolkit helper (re-exported under _ssqrt
# so new_hybrid_pdf can keep `from new_pdf import _ssqrt`, mirroring the Fortran shared helper).
from clubb_jax.src.CLUBB_core.pdf_utilities import _safe_sqrt as _ssqrt


def calc_coef_wpxp2_implicit(wp2, xp2, wpxp, sgn_wpxp, mixt_frac, F_w, F_x,
                             coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
                             coef_sigma_x_1_sqd, coef_sigma_x_2_sqd):
    """Calculate coef_wpxp2_implicit.

    ! Description:
    ! The predictive equation for <x'^2> contains a turbulent advection term of
    ! the form:
    !
    ! - ( 1 / rho_ds ) * d ( rho_ds * <w'x'^2> ) / dz;
    !
    ! where z is height, rho_ds is the dry, base-state density, and <w'x'^2> is
    ! calculated by integrating over the PDF.
    !
    ! This equation is of the form:
    !
    ! <w'x'^2> = coef_wpxp2_implicit * <x'^2>;
    !
    ! In the special case that coef_sigma_w_1_sqd * coef_sigma_x_1_sqd = 0 and
    ! coef_sigma_w_2_sqd * coef_sigma_x_2_sqd = 0, the above equation is
    ! undefined.  However, the equation for this special case can be derived by
    ! taking the original equation for <w'x'^2> and setting both
    ! sigma_w_1 * sigma_x_1 = 0 and sigma_w_2 * sigma_x_2 = 0.
    """
    wp2 = jnp.asarray(wp2, dtype=jnp.float64); xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    wpxp = jnp.asarray(wpxp, dtype=jnp.float64); sgn = jnp.asarray(sgn_wpxp, dtype=jnp.float64)
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64)
    F_w = jnp.asarray(F_w, dtype=jnp.float64); F_x = jnp.asarray(F_x, dtype=jnp.float64)
    cw1 = jnp.asarray(coef_sigma_w_1_sqd, dtype=jnp.float64); cw2 = jnp.asarray(coef_sigma_w_2_sqd, dtype=jnp.float64)
    cx1 = jnp.asarray(coef_sigma_x_1_sqd, dtype=jnp.float64); cx2 = jnp.asarray(coef_sigma_x_2_sqd, dtype=jnp.float64)
    omf = 1.0 - mf

    F_spread = F_x * (omf / mf - mf / omf)
    base = _ssqrt(mf * omf) * _ssqrt(wp2)

    # ( coef_sigma_w_1_sqd * coef_sigma_x_1_sqd = 0
    #   and coef_sigma_w_2_sqd * coef_sigma_x_2_sqd = 0 )
    # or wp2 * xp2 = 0
    branchB = base * _ssqrt(F_w) * (F_spread + (cx1 - cx2))

    cwcx1 = cw1 * cx1
    cwcx2 = cw2 * cx2
    denom = mf * _ssqrt(cwcx1) + omf * _ssqrt(cwcx2)
    # Factor involving coef_sigma_... coefficients
    coefs_factor = (_ssqrt(cwcx1) - _ssqrt(cwcx2)) / jnp.where(denom > 0.0, denom, 1.0)
    denom2 = _ssqrt(wp2) * _ssqrt(xp2)
    wpxp_ratio = wpxp / jnp.where(denom2 > 0.0, denom2, 1.0)
    branchA = base * (_ssqrt(F_w) * F_x * (omf / mf - mf / omf)
                      + _ssqrt(F_w) * (cx1 - cx2)
                      + 2.0 * _ssqrt(F_x) * coefs_factor * sgn * wpxp_ratio
                      - 2.0 * _ssqrt(F_w) * F_x * coefs_factor)

    # Calculate coef_wpxp2_implicit.
    cond = ((cwcx1 > 0.0) | (cwcx2 > 0.0)) & (wp2 * xp2 > 0.0)
    return jnp.where(cond, branchA, branchB)


def calc_setter_var_params(xm, xp2, Skx, sgn_wpxp, F_x, zeta_x):
    """Calculates the PDF component means, the PDF component standard deviations,
    and the mixture fraction for the variable that sets the PDF.

    ! References:
    ! Griffin and Larson (2018)
    """
    xm = jnp.asarray(xm, dtype=jnp.float64); xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    Skx = jnp.asarray(Skx, dtype=jnp.float64); sgn = jnp.asarray(sgn_wpxp, dtype=jnp.float64)
    F_x = jnp.asarray(F_x, dtype=jnp.float64); zeta = jnp.asarray(zeta_x, dtype=jnp.float64)

    # Calculate the mixture fraction.
    mixt_frac = calc_mixture_fraction(Skx, F_x, zeta, sgn)
    omf = 1.0 - mixt_frac
    # Calculate the mean of x in the 1st PDF component.
    mu_x_1 = xm + sgn * _ssqrt(F_x * (omf / mixt_frac) * xp2)
    # Calculate the mean of x in the 2nd PDF component.
    mu_x_2 = xm - (mixt_frac / omf) * (mu_x_1 - xm)
    # Calculate the standard deviation of x in the 1st PDF component.
    # sigma_x_1 = sqrt( ( ( zeta_x + 1 ) * ( 1 - F_x ) )
    #                   / ( ( zeta_x + 2 ) * mixt_frac ) * <x'^2> )
    coef_sigma_x_1_sqd = ((zeta + 1.0) * (1.0 - F_x)) / ((zeta + 2.0) * mixt_frac)
    sigma_x_1 = _ssqrt(coef_sigma_x_1_sqd * xp2)
    # Calculate the standard deviation of x in the 2nd PDF component.
    # sigma_x_2 = sqrt( ( mixt_frac * sigma_x_1^2 )
    #                   / ( ( 1 - mixt_frac ) * ( 1 + zeta_x ) ) )
    #           = sqrt( ( 1 - F_x )
    #                   / ( ( zeta_x + 2 ) * ( 1 - mixt_frac ) ) * <x'^2> )
    coef_sigma_x_2_sqd = (1.0 - F_x) / ((zeta + 2.0) * omf)
    sigma_x_2 = _ssqrt(coef_sigma_x_2_sqd * xp2)
    return mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, mixt_frac, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd


def calc_responder_params(xm, xp2, Skx, sgn_wpxp, F_x, mixt_frac):
    """Calculates the PDF component means and the PDF component standard
    deviations for a responding variable (a variable that is not used to set
    the mixture fraction).

    ! DESCRIPTION OF THE METHOD FOR EACH RESPONDING VARIABLE
    ! ======================================================
    !
    ! In order to find equations for the four PDF parameters for each responding
    ! variable, which are mu_x_1, mu_x_2, sigma_x_1, and sigma_x_2 (where x stands
    ! for a responding variable here), four equations are needed.  These four
    ! equations are the equations for <x>, <x'^2>, and <x'^3> as found by
    ! integrating over the PDF.  Additionally, one more equation, which involves
    ! a tunable parameter F_x, and which is used to help control the spread of the
    ! PDF component means, is used in this equation set.
    !
    ! The resulting equations for the four PDF parameters are:
    !
    ! mu_x_1 = <x> + sqrt( F_x * ( ( 1 - mixt_frac ) / mixt_frac ) * <x'^2> )
    !                * sgn( <w'x'> );
    !
    ! mu_x_2 = <x> - ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_x_1 - <x> );
    !
    ! sigma_x_1^2
    ! = ( ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
    !       - ( 1 + mixt_frac ) * F_x^1.5 + 3 * mixt_frac * sqrt( F_x ) )
    !     / ( 3 * mixt_frac * sqrt( F_x ) ) )
    !   * <x'^2>; and
    !
    ! sigma_x_2^2 = ( ( 1 - F_x ) / ( 1 - mixt_frac ) ) * <x'^2>
    !               - ( mixt_frac / ( 1 - mixt_frac ) ) * sigma_x_1^2;
    !
    ! where Skx is the skewness of x, and Skx = <x'^3> / <x'^2>^(3/2).
    """
    xm = jnp.asarray(xm, dtype=jnp.float64); xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    Skx = jnp.asarray(Skx, dtype=jnp.float64); sgn = jnp.asarray(sgn_wpxp, dtype=jnp.float64)
    F_x = jnp.asarray(F_x, dtype=jnp.float64); mf = jnp.asarray(mixt_frac, dtype=jnp.float64)
    omf = 1.0 - mf
    F_safe = jnp.where(F_x > 0.0, F_x, 1.0)

    # Calculate the mean of x in the 1st PDF component.
    mu_x_1_v = xm + sgn * _ssqrt(F_x * (omf / mf) * xp2)
    # Calculate the mean of x in the 2nd PDF component.
    mu_x_2_v = xm - (mf / omf) * (mu_x_1_v - xm)
    # Calculate the variance of x in the 1st PDF component.
    # sigma_x_1^2
    # = ( ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
    #       - ( 1 + mixt_frac ) * F_x^1.5 + 3 * mixt_frac * sqrt( F_x ) )
    #     / ( 3 * mixt_frac * sqrt( F_x ) ) ) * <x'^2>
    # = ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
    #     / ( 3 * mixt_frac * sqrt( F_x ) )
    #     - ( 1 + mixt_frac ) * F_x / ( 3 * mixt_frac )
    #     + 1 ) * <x'^2>
    c1_v = (_ssqrt(mf * omf) * Skx * sgn / (3.0 * mf * _ssqrt(F_safe))
            - (1.0 + mf) * F_x / (3.0 * mf) + 1.0)
    # Calculate the variance of x in the 2nd PDF component.
    # sigma_x_2^2
    # = ( ( 1 - F_x ) / ( 1 - mixt_frac )
    #     - mixt_frac / ( 1 - mixt_frac )
    #       * ( sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> )
    #           / ( 3 * mixt_frac * sqrt( F_x ) )
    #           - ( 1 + mixt_frac ) * F_x / ( 3 * mixt_frac )
    #           + 1 ) ) * <x'^2>
    # = ( ( ( 1 - F_x ) - mixt_frac * coef_sigma_x_1_sqd )
    #     / ( 1 - mixt_frac ) ) * <x'^2>
    c2_v = ((1.0 - F_x) - mf * c1_v) / omf

    vary = F_x > 0.0
    # When F_x has a value of 0, the PDF becomes a single Gaussian.  This
    # only works when Skx = 0.  However, when Skx /= 0, the value of min_F_x
    # is greater than 0, preventing a problem where F_x = 0 but | Skx | > 0.
    mu_x_1 = jnp.where(vary, mu_x_1_v, xm)
    mu_x_2 = jnp.where(vary, mu_x_2_v, xm)
    coef_sigma_x_1_sqd = jnp.where(vary, c1_v, 1.0)
    coef_sigma_x_2_sqd = jnp.where(vary, c2_v, 1.0)
    sigma_x_1_sqd = coef_sigma_x_1_sqd * xp2
    sigma_x_2_sqd = coef_sigma_x_2_sqd * xp2
    return mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, coef_sigma_x_1_sqd, coef_sigma_x_2_sqd


def sort_roots(roots):
    """Sorts roots from smallest to largest.

    The Fortran implementation uses an explicit three-root sorting network with
    comments for each branch.  jnp.sort is the same operation for the returned
    values.
    """
    return jnp.sort(jnp.asarray(roots, dtype=jnp.float64), axis=-1)


def calc_limits_F_x_responder(mixt_frac, Skx, sgn_wpxp,
                              max_Skx2_pos_Skx_sgn_wpxp, max_Skx2_neg_Skx_sgn_wpxp):
    """Calculates the minimum and maximum allowable values for F_x for a
    responding variable.

    ! Limits on F_x:
    !
    ! Since the PDF parameters for this variable need to work with the mixture
    ! fraction that has been provided by the setting variable, the method does
    ! not work for all values of F_x and Skx.  However, the limits of Skx and F_x
    ! can always be calculated.  The limits are based on keeping the values of
    ! sigma_x_1 and sigma_x_2 greater than or equal to 0.  The equation for
    ! keeping the value of sigma_x_1 greater than or equal to 0 is:
    !
    ! - ( 1 + mixt_frac ) * sqrt( F_x )^3 + 3 * mixt_frac * sqrt( F_x )
    ! + sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ) >= 0.
    !
    ! The roots of sqrt( F_x ) can be solved by an equation of the form:
    !
    ! A * sqrt( F_x )^3 + B * sqrt( F_x )^2 + C * sqrt( F_x ) + D = 0;
    !
    ! where:
    !
    ! A = - ( 1 + mixt_frac );
    ! B = 0;
    ! C = 3 * mixt_frac; and
    ! D = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ).
    !
    ! The equation for keeping the value of sigma_x_2 greater than or equal to 0
    ! is:
    !
    ! - ( 2 - mixt_frac ) * sqrt( F_x )^3 + 3 * ( 1 - mixt_frac ) * sqrt( F_x )
    ! - sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ) >= 0.
    !
    ! The roots of sqrt( F_x ) can be solved by an equation of the form:
    !
    ! A * sqrt( F_x )^3 + B * sqrt( F_x )^2 + C * sqrt( F_x ) + D = 0;
    !
    ! where:
    !
    ! A = - ( 2 - mixt_frac );
    ! B = 0;
    ! C = 3 * ( 1 - mixt_frac ); and
    ! D = - sqrt( mixt_frac * ( 1 - mixt_frac ) ) * Skx * sgn( <w'x'> ).
    !
    ! The value of sqrt( F_x ) is also limited with a minimum of 0 and a maximum
    ! of 1.  The minimum and maximum allowable values of F_x are found by taking
    ! the square of the minimum and maximum allowable values of sqrt( F_x ),
    ! respectively.
    """
    from clubb_jax.src.CLUBB_core.calc_roots import cubic_solve
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64)
    Skx = jnp.asarray(Skx, dtype=jnp.float64)
    sgn = jnp.asarray(sgn_wpxp, dtype=jnp.float64)
    max_pos = jnp.asarray(max_Skx2_pos_Skx_sgn_wpxp, dtype=jnp.float64)
    max_neg = jnp.asarray(max_Skx2_neg_Skx_sgn_wpxp, dtype=jnp.float64)
    omf = 1.0 - mf

    Dterm = _ssqrt(mf * omf) * Skx * sgn

    # Set up the coefficients in the equation for the limit of sqrt(F_x) based
    # on the 1st PDF component standard deviation (sigma_x_1) being greater than
    # or equal to 0.  This equation has the form:
    # A * sqrt(F_x)^3 + B * sqrt(F_x)^2 + C * sqrt(F_x) + D = 0.
    r1 = cubic_solve(-(1.0 + mf), jnp.zeros_like(mf), 3.0 * mf, Dterm)
    # Solve for the roots (values of sqrt(F_x)) that satisfy the above equation.
    r2 = cubic_solve(-(2.0 - mf), jnp.zeros_like(mf), 3.0 * omf, -Dterm)
    r1_real = jnp.real(r1)
    r2_real = jnp.real(r2)
    # Sort the values of the roots (values of sqrt(F_x)) from smallest to
    # largest.  Ignore any complex component of the roots.  The code below that
    # uses sqrt_F_x_roots_1_sorted already factors the appropriate roots to use
    # into account.
    r1s = sort_roots(r1_real)
    # Set up the coefficients in the equation for the limit of sqrt(F_x) based
    # on the 2nd PDF component standard deviation (sigma_x_2) being greater than
    # or equal to 0.  This equation has the form:
    # A * sqrt(F_x)^3 + B * sqrt(F_x)^2 + C * sqrt(F_x) + D = 0.
    # Solve for the roots (values of sqrt(F_x)) that satisfy the above equation.
    # Sort the values of the roots (values of sqrt(F_x)) from smallest to
    # largest.  Ignore any complex component of the roots.  The code below that
    # uses sqrt_F_x_roots_2_sorted already factors the appropriate roots to use
    # into account.
    r2s = sort_roots(r2_real)

    Skx2 = Skx ** 2
    pos = Skx * sgn >= 0.0

    # Find the minimum and maximum allowable values of sqrt(F_x) based on Skx
    # and sgn( <w'x'> ).
    min_sqrt = jnp.where(pos, r2s[..., 1], r1s[..., 1])
    max_pos_branch = jnp.where(Skx2 > max_neg,
                               jnp.minimum(r1_real[..., 0], r2s[..., 2]),
                               jnp.minimum(r1s[..., 2], r2s[..., 2]))
    max_neg_branch = jnp.where(Skx2 > max_pos,
                               jnp.minimum(r2_real[..., 0], r1s[..., 2]),
                               jnp.minimum(r1s[..., 2], r2s[..., 2]))
    max_sqrt = jnp.where(pos, max_pos_branch, max_neg_branch)

    # The minimum and maximum are also limited by 0 and 1, respectively.
    min_sqrt = jnp.maximum(min_sqrt, 0.0)
    max_sqrt = jnp.minimum(max_sqrt, 1.0)
    # The minimum and maximum allowable values for F_x are the squares of the
    # minimum and maximum allowable values for sqrt(F_x).
    return min_sqrt ** 2, max_sqrt ** 2


def calc_coefs_wp2xp_semiimpl(wp2, xp2, sgn_wpxp, mixt_frac, F_w, F_x,
                              coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
                              coef_sigma_x_1_sqd, coef_sigma_x_2_sqd):
    """Calculate coef_wp2xp_implicit and term_wp2xp_explicit.

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
    ! <w'^2 x'> = coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit;
    !
    ! In the special case that coef_sigma_w_1_sqd * coef_sigma_x_1_sqd = 0 and
    ! coef_sigma_w_2_sqd * coef_sigma_x_2_sqd = 0, the above equation is
    ! undefined.  However, the equation for this special case can be derived by
    ! taking the original equation for <w'^2 x'> and setting both
    ! sigma_w_1 * sigma_x_1 = 0 and sigma_w_2 * sigma_x_2 = 0.
    """
    wp2 = jnp.asarray(wp2, dtype=jnp.float64); xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    sgn = jnp.asarray(sgn_wpxp, dtype=jnp.float64); mf = jnp.asarray(mixt_frac, dtype=jnp.float64)
    F_w = jnp.asarray(F_w, dtype=jnp.float64); F_x = jnp.asarray(F_x, dtype=jnp.float64)
    cw1 = jnp.asarray(coef_sigma_w_1_sqd, dtype=jnp.float64); cw2 = jnp.asarray(coef_sigma_w_2_sqd, dtype=jnp.float64)
    cx1 = jnp.asarray(coef_sigma_x_1_sqd, dtype=jnp.float64); cx2 = jnp.asarray(coef_sigma_x_2_sqd, dtype=jnp.float64)
    omf = 1.0 - mf

    base = _ssqrt(mf * omf)
    swp2 = _ssqrt(wp2)
    F_spread = omf / mf - mf / omf
    explicit_prefac = base * _ssqrt(F_x) * _ssqrt(xp2) * wp2 * sgn

    # Calculate coef_wp2xp_implicit and term_wp2xp_explicit.
    cwcx1 = cw1 * cx1
    cwcx2 = cw2 * cx2
    denom = mf * _ssqrt(cwcx1) + omf * _ssqrt(cwcx2)
    # Factor involving coef_sigma_... coefficients
    coefs_factor = (_ssqrt(cwcx1) - _ssqrt(cwcx2)) / jnp.where(denom > 0.0, denom, 1.0)

    coef_full = base * 2.0 * _ssqrt(F_w) * swp2 * coefs_factor
    term_full = explicit_prefac * (F_w * F_spread + (cw1 - cw2) - 2.0 * F_w * coefs_factor)

    coef_red = base * _ssqrt(F_w) * swp2 * F_spread
    term_red = explicit_prefac * (cw1 - cw2)

    # coef_sigma_w_1_sqd * coef_sigma_x_1_sqd = 0
    # and coef_sigma_w_2_sqd * coef_sigma_x_2_sqd = 0
    cond = (cwcx1 > 0.0) | (cwcx2 > 0.0)
    coef_wp2xp_implicit = jnp.where(cond, coef_full, coef_red)
    term_wp2xp_explicit = jnp.where(cond, term_full, term_red)
    return coef_wp2xp_implicit, term_wp2xp_explicit


def calc_coefs_wpxpyp_semiimpl(wp2, xp2, yp2, wpxp, wpyp, sgn_wpxp, sgn_wpyp, mixt_frac,
                               F_w, F_x, F_y, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd,
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
    ! This equation is of the form:
    !
    ! <w'^2 x'> = coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit;
    !
    ! There are also special cases for the above equations.
    """
    wp2 = jnp.asarray(wp2, dtype=jnp.float64); xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    yp2 = jnp.asarray(yp2, dtype=jnp.float64)
    wpxp = jnp.asarray(wpxp, dtype=jnp.float64); wpyp = jnp.asarray(wpyp, dtype=jnp.float64)
    sx = jnp.asarray(sgn_wpxp, dtype=jnp.float64); sy = jnp.asarray(sgn_wpyp, dtype=jnp.float64)
    mf = jnp.asarray(mixt_frac, dtype=jnp.float64)
    F_w = jnp.asarray(F_w, dtype=jnp.float64); F_x = jnp.asarray(F_x, dtype=jnp.float64); F_y = jnp.asarray(F_y, dtype=jnp.float64)
    cw1 = jnp.asarray(coef_sigma_w_1_sqd, dtype=jnp.float64); cw2 = jnp.asarray(coef_sigma_w_2_sqd, dtype=jnp.float64)
    cx1 = jnp.asarray(coef_sigma_x_1_sqd, dtype=jnp.float64); cx2 = jnp.asarray(coef_sigma_x_2_sqd, dtype=jnp.float64)
    cy1 = jnp.asarray(coef_sigma_y_1_sqd, dtype=jnp.float64); cy2 = jnp.asarray(coef_sigma_y_2_sqd, dtype=jnp.float64)
    omf = 1.0 - mf

    def _cfactor(p1, p2):
        s1, s2 = _ssqrt(p1), _ssqrt(p2)
        denom = mf * s1 + omf * s2
        # When coef_sigma_*_1_sqd * coef_sigma_*_1_sqd = 0 and
        # coef_sigma_*_2_sqd * coef_sigma_*_2_sqd = 0, the value of
        # coefs_factor_* is undefined.  However, setting coefs_factor_* to a
        # value of 0 in this scenario allows for the use of general form
        # equations below for coef_wpxpyp_implicit and term_wpxpyp_explicit.
        return jnp.where((p1 > 0.0) | (p2 > 0.0), (s1 - s2) / jnp.where(denom > 0.0, denom, 1.0), 0.0)

    # Calculate coefs_factor_wx.
    factor_wx = _cfactor(cw1 * cx1, cw2 * cx2)
    # Calculate coefs_factor_wy.
    factor_wy = _cfactor(cw1 * cy1, cw2 * cy2)
    # Calculate coefs_factor_xy.
    factor_xy = _cfactor(cx1 * cy1, cx2 * cy2)

    base = _ssqrt(mf * omf)
    sFw_wp2 = _ssqrt(F_w) * _ssqrt(wp2)
    Fx_xp2_sx = _ssqrt(F_x) * _ssqrt(xp2) * sx
    Fy_yp2_sy = _ssqrt(F_y) * _ssqrt(yp2) * sy
    spread = omf / mf - mf / omf

    coef_A = base * sFw_wp2 * factor_xy
    term_A = (base * sFw_wp2 * Fx_xp2_sx * Fy_yp2_sy * (spread - factor_xy - factor_wy - factor_wx)
              + base * Fx_xp2_sx * factor_wy * wpyp
              + base * Fy_yp2_sy * factor_wx * wpxp)
    coef_B = base * sFw_wp2 * (spread - factor_wy - factor_wx)
    term_B = base * Fx_xp2_sx * factor_wy * wpyp + base * Fy_yp2_sy * factor_wx * wpxp

    # Calculate coef_wpxpyp_implicit and term_wpxpyp_explicit.
    cond_A = (cx1 * cy1 > 0.0) | (cx2 * cy2 > 0.0)
    return jnp.where(cond_A, coef_A, coef_B), jnp.where(cond_A, term_A, term_B)
