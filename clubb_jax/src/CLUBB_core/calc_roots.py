"""JAX port of `src/CLUBB_core/calc_roots.F90`.

Porting deviations:
- The Fortran routines take an explicit ``nz`` and loop over 1D arrays.  The
  JAX routines broadcast over any leading array shape and stack roots on the
  last axis.
- ``cubic_solve`` follows the Fortran Cardano formula, but `_cardano_cbrt`
  takes the real, sign-preserving cube root when the Cardano argument is real.
  This is required for mathematically correct roots when ``D > 0`` and
  ``R - sqrt(D)`` is a negative real number; the plain principal complex branch
  can return non-roots for that case.  The deviation is covered by residual and
  set-match tests against `numpy.roots`.
"""
from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()   # cubic_solve uses complex128 (principal-branch ** and sqrt)
import jax.numpy as jnp

_ONE_THIRD = 1.0 / 3.0
_C = jnp.complex128


def _cardano_cbrt(z):
    """Cube root matching gfortran's complex ``z**(1/3)`` as used by the Fortran cubic_solve.

    Subtlety: when the determinant D >= 0, the Cardano arguments ``R ± sqrt(D)`` are EXACTLY real (imag == 0),
    and gfortran's complex power of a real argument returns the **real, sign-preserving** cube root
    (``x<0 -> -|x|^(1/3)``) — NOT the principal complex branch (``|x|^(1/3) e^{iπ/3}``) that a naive
    ``z**(1/3)`` produces. For genuinely complex arguments (D < 0) the principal branch is correct (the
    conjugate pair recombines to real roots). This reproduces the Fortran/f2py exactly. Finite-gradient at 0."""
    re = jnp.real(z)
    im = jnp.imag(z)
    re_abs = jnp.abs(re)
    real_cbrt = jnp.where(re_abs > 0.0,
                          jnp.sign(re) * jnp.where(re_abs > 0.0, re_abs, 1.0) ** _ONE_THIRD,
                          0.0).astype(_C)
    principal = z ** _C(_ONE_THIRD)
    return jnp.where(im == 0.0, real_cbrt, principal)


def cubic_solve(a_coef, b_coef, c_coef, d_coef):
    """Solve for the roots of x in a cubic equation.

    Fortran comments:
      Description:
        Solve for the roots of x in a cubic equation.

        The cubic equation has the form:

        f(x) = a*x^3 + b*x^2 + c*x + d;

        where a /= 0.  When f(x) = 0, the cubic formula is used to solve:

        a*x^3 + b*x^2 + c*x + d = 0.

        The cubic formula is also called Cardano's Formula.

        The three solutions for x are:

        x(1) = -(1/3)*(b/a) + ( S + T );
        x(2) = -(1/3)*(b/a) - (1/2) * ( S + T ) + (1/2)i * sqrt(3) * ( S - T );
        x(3) = -(1/3)*(b/a) - (1/2) * ( S + T ) - (1/2)i * sqrt(3) * ( S - T );

        where:

        S = ( R + sqrt( D ) )^(1/3); and
        T = ( R - sqrt( D ) )^(1/3).

        The determinant, D, is given by:

        D = R^2 + Q^3.

        The values of R and Q relate back to the a, b, c, and d coefficients:

        Q = ( 3*(c/a) - (b/a)^2 ) / 9; and
        R = ( 9*(b/a)*(c/a) - 27*(d/a) - 2*(b/a)^3 ) / 54.

        When D < 0, there are three unique, real-valued roots.  When D = 0, there
        are three real-valued roots, but one root is a double root or a triple
        root.  When D > 0, there is one real-valued root and there are two roots
        that are complex conjugates.

      References:
        http://mathworld.wolfram.com/CubicFormula.html
    """
    a = jnp.asarray(a_coef); b = jnp.asarray(b_coef)
    c = jnp.asarray(c_coef); d = jnp.asarray(d_coef)

    ba = b / a
    ca = c / a
    da = d / a

    # Find the value of the coefficient Q; where
    # Q = ( 3*(c/a) - (b/a)^2 ) / 9.
    cap_Q = (3.0 * ca - ba ** 2) / 9.0

    # Find the value of the coefficient R; where
    # R = ( 9*(b/a)*(c/a) - 27*(d/a) - 2*(b/a)^3 ) / 54.
    cap_R = (9.0 * ba * ca - 27.0 * da - 2.0 * ba ** 3) / 54.0

    # Find the value of the determinant D; where
    # D = R^2 + Q^3.
    determinant = cap_Q ** 3 + cap_R ** 2

    # Calculate the square root of the determinant.  This will be a complex
    # number.
    sqrt_det = jnp.sqrt(determinant.astype(_C))
    R_c = cap_R.astype(_C)
    one_third_c = _C(_ONE_THIRD)

    # Find the value of the coefficient S; where
    # S = ( R + sqrt( D ) )^(1/3).
    cap_S = _cardano_cbrt(R_c + sqrt_det)

    # Find the value of the coefficient T; where
    # T = ( R - sqrt( D ) )^(1/3).
    cap_T = _cardano_cbrt(R_c - sqrt_det)

    sqrt_3 = jnp.sqrt(jnp.asarray(3.0, dtype=_C))
    i_c = _C(1j)
    ba_c = ba.astype(_C)
    base = -one_third_c * ba_c                                          # -(1/3)(b/a)
    SpT = cap_S + cap_T
    SmT = cap_S - cap_T
    # Find the values of the roots.
    # This root is always real-valued.
    # x(1) = -(1/3)*(b/a) + ( S + T ).
    root1 = base + SpT

    # This root is real-valued when D < 0 (even though the square root of the
    # determinant is a complex number), as well as when D = 0 (when it is part
    # of a double or triple root).  When D > 0, this root is a complex number.
    # It is the complex conjugate of roots(3).
    # x(2) = -(1/3)*(b/a) - (1/2) * ( S + T ) + (1/2)i * sqrt(3) * ( S - T ).
    root2 = base - _C(0.5) * SpT + _C(0.5) * i_c * sqrt_3 * SmT

    # This root is real-valued when D < 0 (even though the square root of the
    # determinant is a complex number), as well as when D = 0 (when it is part
    # of a double or triple root).  When D > 0, this root is a complex number.
    # It is the complex conjugate of roots(2).
    # x(3) = -(1/3)*(b/a) - (1/2) * ( S + T ) - (1/2)i * sqrt(3) * ( S - T ).
    root3 = base - _C(0.5) * SpT - _C(0.5) * i_c * sqrt_3 * SmT
    return jnp.stack([root1, root2, root3], axis=-1)


def quadratic_solve(a_coef, b_coef, c_coef):
    """Solve for the roots of x in a quadratic equation.

    Fortran comments:
      Description:
        Solve for the roots of x in a quadratic equation.

        The equation has the form:

        f(x) = a*x^2 + b*x + c;

        where a /= 0.  When f(x) = 0, the quadratic formula is used to solve:

        a*x^2 + b*x + c = 0.

        The two solutions for x are:

        x(1) = ( -b + sqrt( b^2 - 4*a*c ) ) / (2*a); and
        x(2) = ( -b - sqrt( b^2 - 4*a*c ) ) / (2*a).

        The determinant, D, is given by:

        D = b^2 - 4*a*c.

        When D > 0, there are two unique, real-valued roots.  When D = 0, there
        are two real-valued roots, but they are a double root.  When D < 0, there
        there are two roots that are complex conjugates.

      References:
    """
    a = jnp.asarray(a_coef); b = jnp.asarray(b_coef); c = jnp.asarray(c_coef)

    # Find the value of the determinant D; where
    # D = b^2 - 4*a*c.
    determinant = b ** 2 - 4.0 * a * c

    # Calculate the square root of the determinant.  This will be a complex
    # number.
    sqrt_det = jnp.sqrt(determinant.astype(_C))
    b_c = b.astype(_C)
    two_a = (2.0 * a).astype(_C)

    # Find the values of the roots.
    # This root is real-valued when D > 0, as well as when D = 0 (when it is
    # part of a double root).  When D < 0, this root is a complex number.  It is
    # the complex conjugate of roots(2).
    # x(1) = ( -b + sqrt( b^2 - 4*a*c ) ) / (2*a); and
    root1 = (-b_c + sqrt_det) / two_a

    # This root is real-valued when D > 0, as well as when D = 0 (when it is
    # part of a double root).  When D < 0, this root is a complex number.  It is
    # the complex conjugate of roots(1).
    # x(2) = ( -b - sqrt( b^2 - 4*a*c ) ) / (2*a).
    root2 = (-b_c - sqrt_det) / two_a
    return jnp.stack([root1, root2], axis=-1)


def cube_root(x):
    """Calculates the cube root of x.

    Fortran comments:
      Description:
        Calculates the cube root of x.

        When x >= 0, this code simply calculates x^(1/3).  When x < 0, this code
        uses x^(1/3) = -|x|^(1/3).  This eliminates numerical errors when the
        exponent of 1/3 is not treated as exactly 1/3, which would sometimes
        result in values of NaN.

      References:
    """
    x = jnp.asarray(x)
    abs_cbrt = jnp.abs(x) ** _ONE_THIRD
    return jnp.where(x >= 0.0, abs_cbrt, -abs_cbrt)
