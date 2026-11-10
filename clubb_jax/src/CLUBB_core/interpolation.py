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
  * ``zlinterp`` and ``lin_interp_between_grids`` replace explicit Fortran
    searches with vectorized JAX interpolation while preserving zero-fill or
    clamp semantics.
"""
import jax
import sys
from clubb_jax.src.CLUBB_core.error_code import clubb_at_least_debug_level
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import eps


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

    Returns:
        lin_interpolate_two_points

    Fortran declarations (lin_interpolate_two_points in src/CLUBB_core/interpolation.F90):

        !--------------------------- Input Variables ---------------------------

        real( kind = core_rknd ), intent(in) :: &
          height_int,  & ! Height to be interpolated to     [m]
          height_high, & ! Height above the interpolation   [m]
          height_low,  & ! Height below the interpolation   [m]
          var_high,    & ! Variable above the interpolation [units vary]
          var_low        ! Variable below the interpolation [units vary]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: lin_interpolate_two_points
    """
    # Diagnostic-only restoration; retain the existing numerical/STOP boundary.
    jax.lax.cond(
        jnp.any(jnp.abs(height_low - height_high) < 1.0e-12),
        lambda _: jax.debug.callback(
            lambda: print(
                "lin_interpolate_two_points: height_high and height_low cannot be equal.",
                file=sys.stderr,
            ), ordered=True,
        ),
        lambda _: None, operand=None,
    )
    # Compute linear interpolation
    lin_interpolate_two_points = ((height_int - height_low) / (height_high - height_low)) * (var_high - var_low) + var_low
    return lin_interpolate_two_points


def mono_cubic_interp(z_in, km1, k00, kp1, kp2, zm1, z00, zp1, zp2, fm1, f00, fp1, fp2,
                      l_quintic_poly_interp=False):
    """Steffen's monotone cubic interpolation method
    Returns monotone cubic interpolated value between x00 and xp1

    References:
    M. Steffen, Astron. Astrophys. 239, 443-450 (1990)

    Returns:
        f_out

    Fortran declarations (mono_cubic_interp in src/CLUBB_core/interpolation.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          z_in ! The altitude to be interpolated to [m]

        ! k-levels;  their meaning depends on whether we're extrapolating or interpolating
        integer, intent(in) :: &
          km1, k00, kp1, kp2

        real( kind = core_rknd ), intent(in) :: &
          zm1, z00, zp1, zp2, & ! The altitudes for km1, k00, kp1, kp2      [m]
          fm1, f00, fp1, fp2    ! The field at km1, k00, kp1, and kp2       [units vary]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ) :: f_out ! The interpolated field
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
            alpha = (6.0 / jnp.where(jnp.abs(beta) < eps, 1.0, beta)) * h00 * (dfdxp1 - dfdx00) + 0.5
            zn = (z_in - z00) / h00
            quintic = (((beta / 20.0) * zn - (beta * (1.0 + alpha) / 12.0)) * zn
                       + (beta * alpha / 6.0)) * zn ** 2 * zn + dfdx00 * h00 * zn + f00
            linfall = lin_interpolate_two_points(z00, zp1, zm1, fp1, fm1)
            return jnp.where(jnp.abs(beta) < eps, linfall, quintic)

    def extrapolate(_):
        # Linear extrapolation.
        wp1 = (z_in - z00) / (zp1 - z00)
        w00 = 1.0 - wp1
        return wp1 * fp1 + w00 * f00

    return jax.lax.cond(km1 <= k00, interpolate, extrapolate, operand=None)


def linear_interp_factor(factor, var_high, var_low):
    """Determines the coefficient for a linear interpolation

    Returns:
        linear_interp_factor

    Fortran declarations (linear_interp_factor in src/CLUBB_core/interpolation.F90):

        !--------------------------- Input Variables ---------------------------
        real( kind = core_rknd ), intent(in) :: &
          factor,   & ! Factor                           [units vary]
          var_high, & ! Variable above the interpolation [units vary]
          var_low     ! Variable below the interpolation [units vary]
    """
    linear_interp_factor = factor * (var_high - var_low) + var_low
    return linear_interp_factor


def lin_interp_between_grids(ngrdcol, size_interpolate, size_current,
                             interpolate_altitudes, current_altitudes, current_values):
    """Interpolate all columns to the target grids, using the nearest source
    value outside each source grid and preserving the near-level tolerance.

    Returns:
        lin_interp_between_grids

    Fortran declarations (lin_interp_between_grids in src/CLUBB_core/interpolation.F90):

        !--------------------- Input Variables ---------------------
        integer, intent(in) :: &
            ngrdcol, & ! Number of grid columns [-]
            size_interpolate, & ! Number of target levels [-]
            size_current ! Number of source levels [-]

        real( kind = core_rknd ), dimension(ngrdcol,size_interpolate), intent(in) :: &
          interpolate_altitudes ! Target heights [m]

        real( kind = core_rknd ), dimension(ngrdcol,size_current), intent(in) :: &
          current_altitudes ! Source heights [m]

        real( kind = core_rknd ), dimension(ngrdcol,size_current), intent(in) :: &
          current_values ! given values on the current_altitudes altitudes that should be interpolated
                         ! to the interpolate_altitudes altitudes

        !--------------------- Return Variable ---------------------
        real( kind = core_rknd ), dimension(ngrdcol,size_interpolate) :: &
          lin_interp_between_grids ! Field on target heights [units vary]
    """

    # --------------------- Begin Code ---------------------
    tol = 1.0e-6

    if clubb_at_least_debug_level(0):
        jax.lax.cond(
            jnp.any(jnp.diff(current_altitudes) <= 0),
            lambda _: jax.debug.callback(
                lambda: print(
                    "current_altitudes must be sorted for lin_interp_between_grids.",
                    file=sys.stderr,
                ), ordered=True,
            ),
            lambda _: None, operand=None,
        )
        jax.lax.cond(
            jnp.any(jnp.min(interpolate_altitudes, axis=1) < current_altitudes[:, 0])
            | jnp.any(jnp.max(interpolate_altitudes, axis=1) > current_altitudes[:, -1]),
            lambda _: jax.debug.callback(
                lambda: print(
                    "lin_interp_between_grids: Value is out of range. Extrapolating...",
                    file=sys.stderr,
                ), ordered=True,
            ),
            lambda _: None, operand=None,
        )
    interp = jax.vmap(jnp.interp)(interpolate_altitudes, current_altitudes, current_values)

    # find nearest zt level of gr_current for zt(i) of gr_interpolate grid
    # if there is no gr_current zt level above or below the gr_interpolate grid level,
    # then the gr_current value of the closest level is taken
    abs_delta = jnp.abs(interpolate_altitudes[:, :, None] - current_altitudes[:, None, :])
    nearby = abs_delta < tol
    nearest = jnp.argmax(nearby, axis=2)
    exact = jnp.any(nearby, axis=2)
    near_values = jnp.take_along_axis(current_values, nearest, axis=1)
    lin_interp_between_grids = jnp.where(exact, near_values, interp)
    return lin_interp_between_grids


def zlinterp(ngrdcol, dim_out, dim_src, grid_out, grid_src, var_src, nout=None):
    """Interpolate every column independently in the vertical. Optional nout
    limits the active output levels when merged grids have different sizes.
    Assumes values that
    are less than lowest source point are zero and above the highest
    source point are zero. Also assumes altitude increases linearly.

    References:
    function LIN_INT from WRF-HOC

    Returns:
        var_out

    Fortran declarations (zlinterp in src/CLUBB_core/interpolation.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: ngrdcol ! Number of grid columns
        integer, intent(in) :: dim_out, dim_src

        real( kind = core_rknd ), dimension(ngrdcol,dim_src), intent(in) ::  &
          grid_src,  & ! [m]
          var_src      ! [units vary]

        real( kind = core_rknd ), dimension(ngrdcol,dim_out), intent(in) :: &
          grid_out ! [m]

        !--------------------------- Output Variables ---------------------------
        real( kind = core_rknd ), dimension(ngrdcol,dim_out) :: &
          var_out ! [units vary]

        !--------------------------- Optional Input Variables ---------------------------
        integer, intent(in), optional :: nout(ngrdcol) ! Active output levels in each column
    """
    if dim_src < 2:
        return jnp.zeros((ngrdcol, dim_out), dtype=var_src.dtype)

    # searchsorted requires one sorted vector; batching belongs inside this routine.
    upper = jax.vmap(lambda source, target: jnp.searchsorted(source, target, side="left"))(
        grid_src, grid_out,
    )
    upper = jnp.clip(upper, 1, dim_src - 1)
    lower = upper - 1
    z_hi = jnp.take_along_axis(grid_src, upper, axis=1)
    z_lo = jnp.take_along_axis(grid_src, lower, axis=1)
    value_hi = jnp.take_along_axis(var_src, upper, axis=1)
    value_lo = jnp.take_along_axis(var_src, lower, axis=1)
    last_out = jnp.full((ngrdcol,), dim_out) if nout is None else nout
    active = jnp.arange(dim_out)[None, :] < last_out[:, None]
    above_source = active & (grid_out > grid_src[:, -1:])
    active &= jnp.cumsum(above_source, axis=1) == 0
    active &= grid_out >= grid_src[:, :1]
    delta_z = jnp.where(active, grid_out - z_lo, 0.0)
    interval = jnp.where(active, z_hi - z_lo, 1.0)
    value_hi = jnp.where(active, value_hi, 0.0)
    value_lo = jnp.where(active, value_lo, 0.0)
    var_out = (delta_z / interval) * (value_hi - value_lo) + value_lo
    var_out = jnp.where(active, var_out, 0.0)
    return var_out
