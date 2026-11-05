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
- Guarded denominators and `sqrt_clipped` avoid tracing invalid inactive branches.
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

from clubb_jax.src.CLUBB_core.advance_helper_module import sqrt_clipped
from clubb_jax.src.CLUBB_core.new_pdf import calc_coef_wp4_implicit, calc_mixture_fraction


def calc_coef_wp2xp_implicit(wp2, mixt_frac, F_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd):
    """The predictive equation for <w'x'> contains a turbulent advection term of
    the form:

    - ( 1 / rho_ds ) * d ( rho_ds * <w'^2 x'> ) / dz;

    where z is height, rho_ds is the dry, base-state density, and <w'^2 x'> is
    calculated by integrating over the PDF.  The equation for <w'^2 x'> is:

    <w'^2 x'>
    = mixt_frac * ( ( mu_x_1 - <x> ) * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
                    + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
                      * ( mu_w_1 - <w> ) )
      + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )
                              * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 )
                              + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
                                * ( mu_w_2 - <w> ) ).

    The following substitutions are made into the above equation:

    mu_w_1 - <w> = sqrt(F_w) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
                   * sqrt( <w'^2> );

    mu_w_2 - <w> = - sqrt(F_w) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
                     * sqrt( <w'^2> );

    mu_x_1 - <x> = sqrt( ( 1 - mixt_frac ) / mixt_frac )
                   * <w'x'> / sqrt( F_w * <w'^2> );

    mu_x_2 - <x> = - sqrt( mixt_frac / ( 1 - mixt_frac ) )
                     * <w'x'> / sqrt( F_w * <w'^2> );

    sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> );

    sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> );

    sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and

    sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ).

    The equations for coef_sigma_w_1_sqd and coef_sigma_w_2_sqd are:

    coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
                         / ( ( zeta_w + 2 ) * mixt_frac ); and

    coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).

    The equations for coef_sigma_x_1_sqd and coef_sigma_x_2_sqd are:

    coef_sigma_x_1_sqd
    = 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
          * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
      - ( ( 1 + mixt_frac ) / mixt_frac )
        * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ); and

    coef_sigma_x_2_sqd
    = 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
          * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
      + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
        * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ).

    Additionally, corr_w_x_1 = corr_w_x_2 = 0.

    The equation for <w'^2 x'> becomes:

    <w'^2 x'> = sqrt( mixt_frac * ( 1 - mixt_frac ) )
                * ( F_w * ( ( 1 - mixt_frac ) / mixt_frac
                            - mixt_frac / ( 1 - mixt_frac ) )
                    + coef_sigma_w_1_sqd - coef_sigma_w_2_sqd )
                * sqrt( <w'^2> / F_w ) * <w'x'>.

    This equation is of the form:

    <w'^2 x'> = coef_wp2xp_implicit * <w'x'>;

    where:

    coef_wp2xp_implicit = sqrt( mixt_frac * ( 1 - mixt_frac ) )
                          * ( F_w * ( ( 1 - mixt_frac ) / mixt_frac
                                      - mixt_frac / ( 1 - mixt_frac ) )
                              + coef_sigma_w_1_sqd - coef_sigma_w_2_sqd )
                          * sqrt( <w'^2> / F_w ).

    In the special case that F_w = 0, <w'x'> must have a value of 0, in which
    case mu_x_1 - <x> = mu_x_2 - <x> = 0, and <w'^2 x'> = 0.  The value of
    coef_wp2xp_implicit when F_w = 0 and <w'x'> = 0 can be calculated since
    Skw must also have a value of 0 when F_w = 0.  When F_w = 0 and Skw = 0,
    mixt_frac = ( zeta_w + 1 ) / ( zeta_w + 2 ).  When this happens,
    coef_sigma_w_1_sqd - coef_sigma_w_2_sqd = 0.  The equation for
    coef_wp2xp_implicit becomes:

    coef_wp2xp_implicit = sqrt( mixt_frac * ( 1 - mixt_frac ) )
                          * ( ( 1 - mixt_frac ) / mixt_frac
                              - mixt_frac / ( 1 - mixt_frac ) )
                          * sqrt( F_w * <w'^2> );

    and since F_w = 0, coef_wp2xp_implicit = 0.

    Returns:
        coef_wp2xp_implicit

    Fortran declarations (calc_coef_wp2xp_implicit in src/CLUBB_core/new_hybrid_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real ( kind = core_rknd ), intent(in) :: &
          wp2,                & ! Variance of w (overall)                  [m^2/s^2]
          mixt_frac,          & ! Mixture fraction                               [-]
          F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
          coef_sigma_w_1_sqd, & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
          coef_sigma_w_2_sqd    ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]

        !--------------------------- Output Variables ---------------------------
        ! Coefficient: <w'^2 x'> = coef_wp2xp_implicit * <w'x'>
        real ( kind = core_rknd ) :: &
          coef_wp2xp_implicit    ! Coefficient that is multiplied by <w'x'>    [m/s]
    """
    omf = 1.0 - mixt_frac
    F_safe = jnp.where(F_w > 0.0, F_w, 1.0)
    # Calculate coef_wp2xp_implicit.
    coef = (sqrt_clipped(mixt_frac * omf)
            * (F_w * (omf / mixt_frac - mixt_frac / omf) + coef_sigma_w_1_sqd - coef_sigma_w_2_sqd)
            * sqrt_clipped(wp2 / F_safe))
    coef_wp2xp_implicit = jnp.where(F_w > 0.0, coef, 0.0)
    return coef_wp2xp_implicit


# --- new_hybrid_pdf.F90 aliases (Griffin & Larson 2020) ---
# These are byte-identical / sgn=+1 specializations of the new_pdf.F90 forms; CLUBB defines them
# separately in module new_hybrid_pdf with their own (capital-C) names. Kept as thin wrappers so callers of
# either module map to one tested implementation.

def calculate_coef_wp4_implicit(mixt_frac, F_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd):
    """The predictive equation for <w'^3> contains a turbulent advection term of
    the form:

    - ( 1 / rho_ds ) * d ( rho_ds * <w'^4> ) / dz;

    where z is height, rho_ds is the dry, base-state density, and <w'^4> is
    calculated by integrating over the PDF.  The equation for <w'^4> is:

    <w'^4> = mixt_frac * ( 3 * sigma_w_1^4
                           + 6 * ( mu_w_1 - <w> )^2 * sigma_w_1^2
                           + ( mu_w_1 - <w> )^4 )
             + ( 1 - mixt_frac ) * ( 3 * sigma_w_2^4
                                     + 6 * ( mu_w_2 - <w> )^2 * sigma_w_2^2
                                     + ( mu_w_2 - <w> )^4 ).

    The following substitutions are made into the above equation:

    mu_w_1 - <w> = sqrt(F_w) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
                   * sqrt( <w'^2> );

    mu_w_2 - <w> = - sqrt(F_w) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
                     * sqrt( <w'^2> );

    sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> ); and

    sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> ).

    The equations for coef_sigma_w_1_sqd and coef_sigma_w_2_sqd are:

    coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
                         / ( ( zeta_w + 2 ) * mixt_frac ); and

    coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).

    The equation for <w'4> becomes:

    <w'^4> = ( 3 * mixt_frac * coef_sigma_w_1_sqd^2
               + 6 * F_w * ( 1 - mixt_frac ) * coef_sigma_w_1_sqd
               + F_w^2 * ( 1 - mixt_frac )^2 / mixt_frac
               + 3 * ( 1 - mixt_frac ) * coef_sigma_w_2_sqd^2
               + 6 * F_w * mixt_frac * coef_sigma_w_2_sqd
               + F_w^2 * mixt_frac^2 / ( 1 - mixt_frac ) ) * <w'^2>^2.

    This equation is of the form:

    <w'^4> = coef_wp4_implicit * <w'^2>^2;

    where:

    coef_wp4_implicit = 3 * mixt_frac * coef_sigma_w_1_sqd^2
                        + 6 * F_w * ( 1 - mixt_frac ) * coef_sigma_w_1_sqd
                        + F_w^2 * ( 1 - mixt_frac )^2 / mixt_frac
                        + 3 * ( 1 - mixt_frac ) * coef_sigma_w_2_sqd^2
                        + 6 * F_w * mixt_frac * coef_sigma_w_2_sqd
                        + F_w^2 * mixt_frac^2 / ( 1 - mixt_frac ).

    While the <w'^4> term is found in the <w'^3> predictive equation and not
    the <w'^2> predictive equation, the <w'^3> and <w'^2> predictive equations
    are solved together.  This allows the term containing <w'^4> to be solved
    implicitly.

    Fortran repeats the full calculate_coef_wp4_implicit body in this module;
    the JAX function delegates to the shared formula in new_pdf.py instead of
    keeping a second copy of the same algebra.

    Returns:
        coef_wp4_implicit

    Fortran declarations (calculate_coef_wp4_implicit in src/CLUBB_core/new_hybrid_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real ( kind = core_rknd ), intent(in) :: &
          mixt_frac,          & ! Mixture fraction                               [-]
          F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
          coef_sigma_w_1_sqd, & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>      [-]
          coef_sigma_w_2_sqd    ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>      [-]

        !--------------------------- Output Variables ---------------------------
        real ( kind = core_rknd ) :: &
          coef_wp4_implicit    ! Coef.: <w'^4> = coef_wp4_implicit * <w'^2>^2    [-]
    """
    return calc_coef_wp4_implicit(mixt_frac, F_w, coef_sigma_w_1_sqd, coef_sigma_w_2_sqd)


def calculate_mixture_fraction(Skw, F_w, zeta_w):
    """Calculates mixture fraction.

    References:
    Griffin and Larson (2020)

    Fortran keeps this as a separate elemental function; the JAX function
    calls the shared new_pdf mixture-fraction function with sgn_wpxp = +1,
    because w is the setting variable.

    Returns:
        mixt_frac

    Fortran declarations (calculate_mixture_fraction in src/CLUBB_core/new_hybrid_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          Skw,    & ! Skewness of w                                            [-]
          F_w,    & ! Parameter for the spread of the PDF component means of w [-]
          zeta_w    ! Parameter for the PDF component variances of w           [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: &
          mixt_frac    ! Mixture fraction    [-]
    """
    return calc_mixture_fraction(Skw, F_w, zeta_w, 1.0)


def calculate_w_params(wm, wp2, Skw, F_w, zeta_w):
    """Calculates the PDF component means, the PDF component standard deviations,
    and the mixture fraction for the variable that sets the PDF.

    References:
    Griffin and Larson (2020)

    Returns:
        (mu_w_1, mu_w_2, sigma_w_1, sigma_w_2, mixt_frac, coef_sigma_w_1_sqd,
          coef_sigma_w_2_sqd)

    Fortran declarations (calculate_w_params in src/CLUBB_core/new_hybrid_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          wm,     & ! Mean of w (overall)                                     [m/s]
          wp2,    & ! Variance of w (overall)                             [m^2/s^2]
          Skw,    & ! Skewness of w                                             [-]
          F_w,    & ! Parameter for the spread of the PDF component means of w  [-]
          zeta_w    ! Parameter for the PDF component variances of w            [-]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), intent(out) :: &
          mu_w_1,    & ! Mean of w (1st PDF component)                  [m/s]
          mu_w_2,    & ! Mean of w (2nd PDF component)                  [m/s]
          sigma_w_1, & ! Standard deviation of w (1st PDF component)    [m/s]
          sigma_w_2, & ! Standard deviation of w (2nd PDF component)    [m/s]
          mixt_frac    ! Mixture fraction                               [-]

        real( kind = core_rknd ), intent(out) :: &
          coef_sigma_w_1_sqd, & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>     [-]
          coef_sigma_w_2_sqd    ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>     [-]
    """

    # Calculate the mixture fraction.
    mixt_frac = calculate_mixture_fraction(Skw, F_w, zeta_w)
    valid = (mixt_frac > 0.0) & (mixt_frac < 1.0)
    mf = jnp.where(valid, mixt_frac, 0.5)          # avoid div-by-zero in the invalid branch
    omf = 1.0 - mf

    # Calculate the mean of w in the 1st PDF component.
    mu_w_1_v = wm + sqrt_clipped(F_w * (omf / mf) * wp2)
    # Calculate the mean of w in the 2nd PDF component.
    mu_w_2_v = wm - (mf / omf) * (mu_w_1_v - wm)
    # Calculate the standard deviation of w in the 1st PDF component.
    # sigma_w_1 = sqrt( ( ( zeta_w + 1 ) * ( 1 - F_w ) )
    #                   / ( ( zeta_w + 2 ) * mixt_frac ) * <w'^2> )
    c1_v = ((zeta_w + 1.0) * (1.0 - F_w)) / ((zeta_w + 2.0) * mf)
    # Calculate the standard deviation of w in the 2nd PDF component.
    # sigma_w_2 = sqrt( ( mixt_frac * sigma_w_1^2 )
    #                   / ( ( 1 - mixt_frac ) * ( 1 + zeta_w ) ) )
    #           = sqrt( ( 1 - F_w )
    #                   / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ) * <w'^2> )
    c2_v = (1.0 - F_w) / ((zeta_w + 2.0) * omf)
    sig1_v = sqrt_clipped(c1_v * wp2); sig2_v = sqrt_clipped(c2_v * wp2)

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

    References:
    Griffin and Larson (2020)

    Returns:
        (mu_x_1, mu_x_2, sigma_x_1_sqd, sigma_x_2_sqd, coef_sigma_x_1_sqd,
          coef_sigma_x_2_sqd)

    Fortran declarations (calculate_responder_params in src/CLUBB_core/new_hybrid_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          xm,       & ! Mean of x (overall)                             [units vary]
          xp2,      & ! Variance of x (overall)                     [(units vary)^2]
          Skx,      & ! Skewness of x                                            [-]
          wpxp,     & ! Covariance of w and x (overall)            [m/s(units vary)]
          wp2,      & ! Variance of w (overall)                            [m^2/s^2]
          F_w,      & ! Parameter for the spread of the PDF component means of w [-]
          mixt_frac   ! Mixture fraction                                         [-]

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
    omf = 1.0 - mixt_frac

    gate = jnp.abs(wpxp) > 0.0
    # Note:  when |<w'x'>| > 0, F_w, <w'^2>, and <x'^2> must all have values
    #        greater than 0.
    wpxp_safe = jnp.where(gate, wpxp, 1.0)              # guard the /(3·wpxp) and /√(F_w·wp2) terms
    fw_wp2_safe = jnp.where(gate, F_w * wp2, 1.0)       # >0 whenever gate is true (Fortran note)
    den_safe = jnp.where(gate, 3.0 * F_w * wp2 * xp2, 1.0)

    # Calculate the mean of x in the 1st PDF component.
    mu_x_1_v = xm + sqrt_clipped(omf / mixt_frac) * wpxp / sqrt_clipped(fw_wp2_safe)
    # Calculate the mean of x in the 2nd PDF component.
    mu_x_2_v = xm - sqrt_clipped(mixt_frac / omf) * wpxp / sqrt_clipped(fw_wp2_safe)
    # Calculate the variance of x in the 1st PDF component.
    # sigma_x_1^2
    # = ( 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
    #         * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
    #     - ( ( 1 + mixt_frac ) / mixt_frac )
    #       * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ) )
    #   * <x'^2>
    c1_v = (1.0 + sqrt_clipped(omf / mixt_frac) * Skx * sqrt_clipped(fw_wp2_safe * xp2) / (3.0 * wpxp_safe)
            - ((1.0 + mixt_frac) / mixt_frac) * wpxp ** 2 / den_safe)
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
    c2_v = (1.0 - sqrt_clipped(mixt_frac / omf) * Skx * sqrt_clipped(fw_wp2_safe * xp2) / (3.0 * wpxp_safe)
            + ((mixt_frac - 2.0) / omf) * wpxp ** 2 / den_safe)
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
    """The predictive equation for <x'^2> contains a turbulent advection term of
    the form:

    - ( 1 / rho_ds ) * d ( rho_ds * <w'x'^2> ) / dz;

    where z is height, rho_ds is the dry, base-state density, and <w'x'^2> is
    calculated by integrating over the PDF.  The equation for <w'x'^2> is:

    <w'x'^2>
    = mixt_frac * ( ( mu_w_1 - <w> ) * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
                    + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
                      * ( mu_x_1 - <x> ) )
      + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )
                              * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 )
                              + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
                                * ( mu_x_2 - <x> ) ).

    The following substitutions are made into the above equation:

    mu_w_1 - <w> = sqrt(F_w) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
                   * sqrt( <w'^2> );

    mu_w_2 - <w> = - sqrt(F_w) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
                     * sqrt( <w'^2> );

    mu_x_1 - <x> = sqrt( ( 1 - mixt_frac ) / mixt_frac )
                   * <w'x'> / sqrt( F_w * <w'^2> );

    mu_x_2 - <x> = - sqrt( mixt_frac / ( 1 - mixt_frac ) )
                     * <w'x'> / sqrt( F_w * <w'^2> );

    sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> );

    sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> );

    sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and

    sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ).

    The equations for coef_sigma_w_1_sqd and coef_sigma_w_2_sqd are:

    coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
                         / ( ( zeta_w + 2 ) * mixt_frac ); and

    coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).

    The equations for coef_sigma_x_1_sqd and coef_sigma_x_2_sqd are:

    coef_sigma_x_1_sqd
    = 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
          * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
      - ( ( 1 + mixt_frac ) / mixt_frac )
        * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ); and

    coef_sigma_x_2_sqd
    = 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
          * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
      + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
        * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ).

    Additionally, corr_w_x_1 = corr_w_x_2 = 0.

    The equation for <w'x'^2> becomes:

    <w'x'^2>
    = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w * <w'^2> )
      * ( coef_sigma_x_1_sqd - coef_sigma_x_2_sqd ) * <x'^2>
      + sqrt( mixt_frac * ( 1 - mixt_frac ) )
        * <w'x'>^2 / sqrt( F_w * <w'^2> )
        * ( ( 1 - mixt_frac ) / mixt_frac  - mixt_frac / ( 1 - mixt_frac ) ).

    This equation is of the form:

    <w'x'^2> = coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit;

    where:

    coef_wpxp2_implicit
    = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w * <w'^2> )
      * ( coef_sigma_x_1_sqd - coef_sigma_x_2_sqd ); and

    term_wpxp2_explicit
    = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * <w'x'>^2 / sqrt( F_w * <w'^2> )
      * ( ( 1 - mixt_frac ) / mixt_frac  - mixt_frac / ( 1 - mixt_frac ) ).

    In the special case that F_w = 0, mu_w_1 - <w> = mu_w_2 - <w> = 0, and
    <w'x'^2> = 0.  The value of coef_wp2xp_implicit = 0 when F_w = 0.
    Likewise, the value of <w'x'> = 0 when F_w = 0, which makes
    term_wp2xp_explicit undefined in this scenario.  However, since both
    <w'x'^2> = 0 and coef_wp2xp_implicit = 0, term_wp2xp_explicit = 0 in this
    situation.

    Fortran declarations (calc_coefs_wpxp2_semiimpl in src/CLUBB_core/new_hybrid_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real ( kind = core_rknd ), intent(in) :: &
          wp2,                & ! Variance of w (overall)                  [m^2/s^2]
          wpxp,               & ! Covariance of w and x           [m/s (units vary)]
          mixt_frac,          & ! Mixture fraction                               [-]
          F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
          coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
          coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]

        !--------------------------- Output Variables ---------------------------
        ! Coefs.: <w'x'^2> = coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit
        real ( kind = core_rknd ), intent(out) :: &
          coef_wpxp2_implicit, & ! Coefficient that is multiplied by <x'^2>  [m/s]
          term_wpxp2_explicit    ! Term that is on the RHS    [m/s (units vary)^2]
    """
    omf = 1.0 - mixt_frac

    # Calculate coef_wpxp2_implicit and term_wpxp2_explicit.
    vary = (F_w > 0.0) & (wp2 > 0.0)
    sFwwp2 = sqrt_clipped(F_w * wp2)
    sFwwp2_safe = jnp.where(sFwwp2 > 0.0, sFwwp2, 1.0)
    base = sqrt_clipped(mixt_frac * omf)
    coef = base * sFwwp2 * (coef_sigma_x_1_sqd - coef_sigma_x_2_sqd)
    term = base * wpxp ** 2 / sFwwp2_safe * (omf / mixt_frac - mixt_frac / omf)
    return jnp.where(vary, coef, 0.0), jnp.where(vary, term, 0.0)


def calc_coefs_wpxpyp_semiimpl(wp2, wpxp, wpyp, mixt_frac, F_w,
                               coef_sigma_x_1_sqd, coef_sigma_x_2_sqd,
                               coef_sigma_y_1_sqd, coef_sigma_y_2_sqd):
    """The predictive equation for <w'x'y'> contains a turbulent advection term
    of the form:

    - ( 1 / rho_ds ) * d ( rho_ds * <w'x'y'> ) / dz;

    where z is height, rho_ds is the dry, base-state density, and <w'x'y'> is
    calculated by integrating over the PDF.  The equation for <w'x'y'> is:

    <w'x'y'>
    = mixt_frac
      * ( ( mu_w_1 - <w> ) * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
          + corr_x_y_1 * sigma_x_1 * sigma_y_1 * ( mu_w_1 - <w> )
          + corr_w_y_1 * sigma_w_1 * sigma_y_1 * ( mu_x_1 - <x> )
          + corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_y_1 - <y> ) )
      + ( 1 - mixt_frac )
        * ( ( mu_w_2 - <w> ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
            + corr_x_y_2 * sigma_x_2 * sigma_y_2 * ( mu_w_2 - <w> )
            + corr_w_y_2 * sigma_w_2 * sigma_y_2 * ( mu_x_2 - <x> )
            + corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_y_2 - <y> ) ).

    The following substitutions are made into the above equation:

    mu_w_1 - <w> = sqrt(F_w) * sqrt( ( 1 - mixt_frac ) / mixt_frac )
                   * sqrt( <w'^2> );

    mu_w_2 - <w> = - sqrt(F_w) * sqrt( mixt_frac / ( 1 - mixt_frac ) )
                     * sqrt( <w'^2> );

    mu_x_1 - <x> = sqrt( ( 1 - mixt_frac ) / mixt_frac )
                   * <w'x'> / sqrt( F_w * <w'^2> );

    mu_x_2 - <x> = - sqrt( mixt_frac / ( 1 - mixt_frac ) )
                     * <w'x'> / sqrt( F_w * <w'^2> );

    mu_y_1 - <y> = sqrt( ( 1 - mixt_frac ) / mixt_frac )
                   * <w'y'> / sqrt( F_w * <w'^2> );

    mu_y_2 - <y> = - sqrt( mixt_frac / ( 1 - mixt_frac ) )
                     * <w'y'> / sqrt( F_w * <w'^2> );

    sigma_w_1 = sqrt( coef_sigma_w_1_sqd * <w'^2> );

    sigma_w_2 = sqrt( coef_sigma_w_2_sqd * <w'^2> );

    sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> );

    sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> );

    sigma_y_1 = sqrt( coef_sigma_y_1_sqd * <y'^2> ); and

    sigma_y_2 = sqrt( coef_sigma_y_2_sqd * <y'^2> ).

    The equations for coef_sigma_w_1_sqd and coef_sigma_w_2_sqd are:

    coef_sigma_w_1_sqd = ( ( zeta_w + 1 ) * ( 1 - F_w ) )
                         / ( ( zeta_w + 2 ) * mixt_frac ); and

    coef_sigma_w_2_sqd = ( 1 - F_w ) / ( ( zeta_w + 2 ) * ( 1 - mixt_frac ) ).

    The equations for coef_sigma_x_1_sqd and coef_sigma_x_2_sqd are:

    coef_sigma_x_1_sqd
    = 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
          * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
      - ( ( 1 + mixt_frac ) / mixt_frac )
        * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ); and

    coef_sigma_x_2_sqd
    = 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
          * Skx * sqrt( F_w * <w'^2> * <x'^2> ) / ( 3 * <w'x'> )
      + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
        * <w'x'>^2 / ( 3 * F_w * <w'^2> * <x'^2> ).

    The equations for coef_sigma_y_1_sqd and coef_sigma_y_2_sqd are:

    coef_sigma_y_1_sqd
    = 1 + sqrt( ( 1 - mixt_frac ) / mixt_frac )
          * Sky * sqrt( F_w * <w'^2> * <y'^2> ) / ( 3 * <w'y'> )
      - ( ( 1 + mixt_frac ) / mixt_frac )
        * <w'y'>^2 / ( 3 * F_w * <w'^2> * <y'^2> ); and

    coef_sigma_y_2_sqd
    = 1 - sqrt( mixt_frac / ( 1 - mixt_frac ) )
          * Sky * sqrt( F_w * <w'^2> * <y'^2> ) / ( 3 * <w'y'> )
      + ( ( mixt_frac - 2 ) / ( 1 - mixt_frac ) )
        * <w'y'>^2 / ( 3 * F_w * <w'^2> * <y'^2> ).

    Additionally, corr_w_x_1 = corr_w_x_2 = corr_w_y_1 = corr_w_y_2 = 0; and:

    corr_x_y_1 = corr_x_y_2
    = ( <x'y'> - mixt_frac * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
               - ( 1 - mixt_frac ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> ) )
      / ( mixt_frac * sigma_x_1 * sigma_y_1
          + ( 1 - mixt_frac ) * sigma_x_2 * sigma_y_2 );

    where -1 <= corr_x_y_1 = corr_x_y_2 <= 1.  This equation can be rewritten
    as:

    corr_x_y_1 = corr_x_y_2
    = ( <x'y'>
        - sqrt( F_x ) * sqrt( F_y ) * sgn( <w'x'> ) * sgn( <w'y'> )
          * sqrt( <x'^2> ) * sqrt( <y'^2 > ) )
      / ( ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
            + ( 1 - mixt_frac )
              * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
          * sqrt( <x'^2> ) * sqrt( <y'^2> ) ).

    The equation for <w'x'y'> becomes:

    <w'x'y'>
    = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w * <w'^2> )
      * ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
                - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
        / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
            + ( 1 - mixt_frac )
              * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
      * <x'y'>
      + sqrt( mixt_frac * ( 1 - mixt_frac ) )
        * <w'x'> * <w'y'> / sqrt( F_w * <w'^2> )
        * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac )
            - ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
                - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
              / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
                  + ( 1 - mixt_frac )
                    * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) ) )

    This equation is of the form:

    <w'x'y'> = coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit;

    where:

    coef_wpxpyp_implicit
    = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w * <w'^2> )
      * ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
                - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
        / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
            + ( 1 - mixt_frac )
              * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) ); and

    term_wpxpyp_explicit
    = sqrt( mixt_frac * ( 1 - mixt_frac ) )
      * <w'x'> * <w'y'> / sqrt( F_w * <w'^2> )
      * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac )
          - ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
              - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
            / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
                + ( 1 - mixt_frac )
                  * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) ) ).

    There are also special cases for the above equations.  In the scenario
    that sigma_x_1 * sigma_y_1 = sigma_x_2 * sigma_y_2 = 0, and equations for
    coef_wpxpyp_implicit and term_wpxpyp_explicit become:

    coef_wpxpyp_implicit
    = sqrt( mixt_frac * ( 1 - mixt_frac ) ) * sqrt( F_w * <w'^2> )
      * ( ( 1 - mixt_frac ) / mixt_frac - mixt_frac / ( 1 - mixt_frac ) ); and

    term_wpxpyp_explicit = 0.

    In the scenario where F_w = 0, mu_w_1 - <w> = mu_w_2 - <w> = 0, and
    <w'x'> = <w'y'> = 0, which means mu_x_1 - <x> = mu_x_2 - <x> = 0, as well
    as mu_y_1 - <y> = mu_y_2 - <y> = 0, and as a result, <w'x'y'> = 0.  When
    F_w = 0, coef_wpxpyp_implicit = 0.  Since <w'x'y'> = 0 and
    coef_wpxpyp_implicit = 0 when F_w = 0, term_wpxpyp_explicit = 0.

    Returns:
        (coef, term)

    Fortran declarations (calc_coefs_wpxpyp_semiimpl in src/CLUBB_core/new_hybrid_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        real ( kind = core_rknd ), intent(in) :: &
          wp2,                & ! Variance of w (overall)                  [m^2/s^2]
          wpxp,               & ! Covariance of w and x (overall)     [m/s(x units)]
          wpyp,               & ! Covariance of w and y (overall)     [m/s(y units)]
          mixt_frac,          & ! Mixture fraction                               [-]
          F_w,                & ! Parameter: spread of the PDF comp. means of w  [-]
          coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
          coef_sigma_x_2_sqd, & ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]
          coef_sigma_y_1_sqd, & ! sigma_y_1^2 = coef_sigma_y_1_sqd * <y'^2>      [-]
          coef_sigma_y_2_sqd    ! sigma_y_2^2 = coef_sigma_y_2_sqd * <y'^2>      [-]

        !--------------------------- Output Variables ---------------------------
        ! Coefs.: <w'x'y'> = coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit
        real ( kind = core_rknd ), intent(out) :: &
          coef_wpxpyp_implicit, & ! Coefficient that is multiplied by <x'y'>   [m/s]
          term_wpxpyp_explicit    ! Term that is on the RHS  [m/s(x units)(y units)]
    """
    omf = 1.0 - mixt_frac

    # Calculate coef_wpxpyp_implicit and term_wpxpyp_explicit.
    s1 = sqrt_clipped(coef_sigma_x_1_sqd * coef_sigma_y_1_sqd)
    s2 = sqrt_clipped(coef_sigma_x_2_sqd * coef_sigma_y_2_sqd)
    xy_vary = (((coef_sigma_x_1_sqd * coef_sigma_y_1_sqd > 0.0)
                | (coef_sigma_x_2_sqd * coef_sigma_y_2_sqd > 0.0))
               & (F_w > 0.0) & (wp2 > 0.0))
    denom = mixt_frac * s1 + omf * s2
    denom_safe = jnp.where(denom > 0.0, denom, 1.0)
    # coefs_factor_xy
    # = ( sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    #     - sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    #   / ( mixt_frac * sqrt( coef_sigma_x_1_sqd * coef_sigma_y_1_sqd )
    #       + ( 1 - mixt_frac )
    #         * sqrt( coef_sigma_x_2_sqd * coef_sigma_y_2_sqd ) )
    f_xy = (s1 - s2) / denom_safe

    sFwwp2 = sqrt_clipped(F_w * wp2)
    sFwwp2_safe = jnp.where(sFwwp2 > 0.0, sFwwp2, 1.0)
    base = sqrt_clipped(mixt_frac * omf)

    coef_vary = base * sFwwp2 * f_xy
    term_vary = base * wpxp * wpyp / sFwwp2_safe * (omf / mixt_frac - mixt_frac / omf - f_xy)
    coef_else = jnp.where(F_w > 0.0, base * sFwwp2 * (omf / mixt_frac - mixt_frac / omf), 0.0)

    coef = jnp.where(xy_vary, coef_vary, coef_else)
    term = jnp.where(xy_vary, term_vary, 0.0)
    return coef, term
