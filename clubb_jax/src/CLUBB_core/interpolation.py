"""JAX port of ``interpolation.F90``.

Porting deviations:
  * Fortran exposes ``binary_search`` and ``lin_interpolate_on_grid_api``.
    JAX does not keep them as public routines because callers use vectorized
    ``jnp.interp``-based routines instead.
  * Fortran debug builds stop on invalid linear-interpolation inputs and
    unsorted grids.  JAX routines assume validated monotone grids so they can
    remain jit-friendly.
  * Fortran ``l_quintic_poly_interp`` is a module flag; JAX passes
    ``l_quintic_poly_interp`` as an explicit argument.
  * ``zlinterp_fnc`` and ``lin_interp_between_grids`` replace explicit Fortran
    searches with vectorized JAX interpolation while preserving zero-fill or
    clamp semantics.
"""
import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

_EPS = 1.0e-10   # constants_clubb.F90 eps


def lin_interpolate_two_points(height_int, height_high, height_low, var_high, var_low):
    """This function computes a linear interpolation of the value of variable.

    Given two known values of a variable at two height values, the value
    of that variable at a height between those two height levels (rather
    than a height outside of those two height levels) is computed.

    Here is a diagram:

     ################################ Height high, know variable value



     -------------------------------- Height to be interpolated to; linear interpolation





     ################################ Height low, know variable value


    FORMULA:

    variable(@ Height interpolation) =

    [ (variable(@ Height high) - variable(@ Height low)) / (Height high - Height low) ]
    * (Height interpolation - Height low)  +  variable(@ Height low)

    Comments from WRF-HOC, Brian Griffin.

    References:
    None
    """
    # Compute linear interpolation
    return ((height_int - height_low) / (height_high - height_low)) * (var_high - var_low) + var_low


def mono_cubic_interp(z_in, km1, k00, kp1, kp2, zm1, z00, zp1, zp2, fm1, f00, fp1, fp2,
                      l_quintic_poly_interp=False):
    """Steffen's monotone cubic interpolation method.

    Returns monotone cubic interpolated value between x00 and xp1

    Original Author:
      Takanobu Yamaguchi
      tak.yamaguchi@noaa.gov

      This version has been modified slightly for CLUBB's coding standards and
      adds the 3/2 from eqn 21. -dschanen 26 Oct 2011
      We have also added a quintic polynomial option.

    References:
      M. Steffen, Astron. Astrophys. 239, 443-450 (1990)
    """
    coef1 = 1.0
    coef2 = 1.0

    def _lim(sa, sb, p):
        # Steffen's limited derivative: (sign(sa)+sign(sb)) * min(|sa|, |sb|, |p|/2 * coef2).
        return (coef1 * (jnp.sign(sa) + jnp.sign(sb))
                * jnp.minimum(jnp.minimum(jnp.abs(sa), jnp.abs(sb)), coef2 * 0.5 * jnp.abs(p)))

    # ---- Begin Code ----

    # TODO(port-mirror): lax.cond requires callables for the source's
    # data-dependent branches; remove them if JAX supports traced Python ifs.
    def interpolate(_):
        hm1 = z00 - zm1
        h00 = zp1 - z00
        hp1 = zp2 - zp1

        def lower_endpoint(_):
            s00 = (fp1 - f00) / (zp1 - z00)
            sp1 = (fp2 - fp1) / (zp2 - zp1)
            dfdx00 = s00
            pp1 = (s00 * hp1 + sp1 * h00) / (h00 + hp1)
            dfdxp1 = _lim(s00, sp1, pp1)
            return s00, dfdx00, dfdxp1

        def upper_endpoint(_):
            sm1 = (f00 - fm1) / (z00 - zm1)
            s00 = (fp1 - f00) / (zp1 - z00)
            p00 = (sm1 * h00 + s00 * hm1) / (hm1 + h00)
            dfdx00 = _lim(sm1, s00, p00)
            dfdxp1 = s00
            return s00, dfdx00, dfdxp1

        def interior(_):
            sm1 = (f00 - fm1) / (z00 - zm1)
            s00 = (fp1 - f00) / (zp1 - z00)
            sp1 = (fp2 - fp1) / (zp2 - zp1)
            p00 = (sm1 * h00 + s00 * hm1) / (hm1 + h00)
            pp1 = (s00 * hp1 + sp1 * h00) / (h00 + hp1)
            dfdx00 = _lim(sm1, s00, p00)
            dfdxp1 = _lim(s00, sp1, pp1)
            return s00, dfdx00, dfdxp1

        s00, dfdx00, dfdxp1 = jax.lax.cond(
            km1 == k00,
            lower_endpoint,
            lambda operand: jax.lax.cond(
                kp1 == kp2,
                upper_endpoint,
                interior,
                operand,
            ),
            operand=None,
        )

        if not l_quintic_poly_interp:
            # Old formula
            # f_out = c1 * ( (z_in - z00)**3 ) + c2 * ( (z_in - z00)**2 ) + c3 * (z_in - z00) + c4
            c1 = (dfdx00 + dfdxp1 - 2.0 * s00) / (h00 ** 2)
            c2 = (3.0 * s00 - 2.0 * dfdx00 - dfdxp1) / h00
            c3 = dfdx00
            c4 = f00
            # Faster nested multiplication
            zprime = z_in - z00
            return c4 + zprime * (c3 + zprime * (c2 + zprime * c1))
        else:
            # Use a quintic polynomial interpolation instead instead of the Steffen formula.
            # Unlike the formula above, this formula does not guarantee monotonicity.
            beta = 120.0 * ((fp1 - f00) - 0.5 * h00 * (dfdx00 + dfdxp1))

            # Prevent an underflow by using a linear interpolation
            alpha = (6.0 / jnp.where(jnp.abs(beta) < _EPS, 1.0, beta)) * h00 * (dfdxp1 - dfdx00) + 0.5
            zn = (z_in - z00) / h00
            quintic = (((beta / 20.0) * zn - (beta * (1.0 + alpha) / 12.0)) * zn
                       + (beta * alpha / 6.0)) * zn ** 2 * zn + dfdx00 * h00 * zn + f00
            linfall = lin_interpolate_two_points(z00, zp1, zm1, fp1, fm1)
            return jnp.where(jnp.abs(beta) < _EPS, linfall, quintic)

    def extrapolate(_):
        # Linear extrapolation.
        wp1 = (z_in - z00) / (zp1 - z00)
        w00 = 1.0 - wp1
        return wp1 * fp1 + w00 * f00

    return jax.lax.cond(km1 <= k00, interpolate, extrapolate, operand=None)


def linear_interp_factor(factor, var_high, var_low):
    """Determines the coefficient for a linear interpolation.

    References:
      None
    """
    return factor * (var_high - var_low) + var_low


def lin_interp_between_grids(interpolate_altitudes, current_altitudes, current_values):
    """Interpolates from the values (gr_current_values) on the current grid (gr_current) to the gr_interpolate grid."""
    interpolate_altitudes = jnp.asarray(interpolate_altitudes, dtype=jnp.float64)
    current_altitudes = jnp.asarray(current_altitudes, dtype=jnp.float64)
    current_values = jnp.asarray(current_values, dtype=jnp.float64)

    # --------------------- Begin Code ---------------------
    tol = 1.0e-6

    interp = jnp.interp(interpolate_altitudes, current_altitudes, current_values)

    # find nearest zt level of gr_current for zt(i) of gr_interpolate grid
    # if there is no gr_current zt level above or below the gr_interpolate grid level,
    # then the gr_current value of the closest level is taken
    abs_delta = jnp.abs(interpolate_altitudes[:, None] - current_altitudes[None, :])
    nearest = jnp.argmin(abs_delta, axis=1)
    exact = jnp.min(abs_delta, axis=1) < tol
    return jnp.where(exact, current_values[nearest], interp)


def zlinterp_fnc(grid_out, grid_src, var_src):
    """Do a linear interpolation in the vertical.

    Assumes values that
    are less than lowest source point are zero and above the highest
    source point are zero. Also assumes altitude increases linearly.

    References:
      function LIN_INT from WRF-HOC
    """
    grid_out = jnp.asarray(grid_out, dtype=jnp.float64)
    grid_src = jnp.asarray(grid_src, dtype=jnp.float64)
    var_src = jnp.asarray(var_src, dtype=jnp.float64)
    # ---- Begin Code ----
    return jnp.interp(grid_out, grid_src, var_src, left=0.0, right=0.0)


def plinterp_fnc(grid_out, grid_src, var_src):
    """Do a linear interpolation in the vertical with pressures.

    Assumes
    values that are less than lowest source point are zero and above the
    highest source point are zero. Also assumes altitude increases linearly.
    This function just calls zlinterp_fnc, but negates grid_out and grid_src.

    References:
      function LIN_INT from WRF-HOC
    """
    grid_out = jnp.asarray(grid_out, dtype=jnp.float64)
    grid_src = jnp.asarray(grid_src, dtype=jnp.float64)
    # ---- Begin Code ----
    return zlinterp_fnc(-grid_out, -grid_src, var_src)
