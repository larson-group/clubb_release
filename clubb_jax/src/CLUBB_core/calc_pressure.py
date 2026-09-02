"""JAX implementations of selected routines from ``calc_pressure.F90``.

Porting deviations:
The repeated Fortran Exner-update branches are factored into ``_exner_step``
so the thermodynamic-level scan and momentum-level vector operation use the
same formula.  Fortran subroutines mutate output arrays and loop over explicit
dimensions; these JAX routines return arrays.  ``init_pressure`` takes the
grid object directly and uses ``lax.scan`` for the vertical recurrence.
"""

from functools import partial

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

from clubb_jax.src.CLUBB_core.constants_clubb import (
    Cp,
    Lv,
    ep1,
    ep2,
    grav,
    kappa,
    p0,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.grid_class import zt2zm


def _exner_step(exner_1, z_1, z_2, thvm_1, thvm_2):
    """One hydrostatic Exner integration step from the formulas in init_pressure."""
    g_ov_cp = grav / Cp
    diff = thvm_2 - thvm_1
    tol = jnp.finfo(jnp.asarray(thvm_2).dtype).eps * thvm_2
    use_log = jnp.abs(diff) > tol
    safe_diff = jnp.where(use_log, diff, 1.0)
    safe_ratio = jnp.where(use_log, thvm_2 / thvm_1, 1.0)

    # exner2
    # = exner1
    #     | ( grav / Cp )
    #     | * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    #   - | where thvm2 /= thvm1;
    log_step = (
        exner_1
        - g_ov_cp * (z_2 - z_1) / safe_diff * jnp.log(safe_ratio)
    )

    #     |
    #     | ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).
    constant_step = exner_1 - g_ov_cp * (z_2 - z_1) / thvm_2
    return jnp.where(use_log, log_step, constant_step)


def init_pressure(thvm, p_sfc, gr):
    """Calculate hydrostatic pressure and Exner on thermodynamic and momentum levels.

    Description:
    Calculates the initial pressure according to the hydrostatic
    approximation.  Combining the moist equation of state and the hydrostatic
    approximation, the change of pressure with respect to height can be
    calculated based on theta_v, such that:

    dp/dz = - p * grav / ( Rd * theta_v * exner );

    where exner = ( p / p0 )^(Rd/Cp);

    and where p0 is a reference pressure of 100000 Pa.

    The integral equation is set up to integrate over p on the left-hand side
    and integrate over z on the right-hand side.  The equation is:

    INT(p1:p2) p^(Rd/Cp-1) dp
    = - p0^(Rd/Cp) * ( grav / Rd ) * INT(z1:z2) (1/thvm) dz.

    The value of mean theta_v (thvm) is calculated at each thermodynamic grid
    level, and linear interpolation is used in the integral equation for all
    altitude in-between successive thermodynamic levels, such that:

    thvm(z) = ( ( thvm2 - thvm1 ) / ( z2 - z1 ) ) * ( z - z1 ) + thvm1.

    The integrals are solved, and the results for pressure can be rewritten
    in terms of exner, such that:

    exner2 - exner1
      | - ( grav / Cp )
      |   * ( ( z2 - z1 ) / ( thvm2 - thvm1 ) ) * ln( thvm2 / thvm1 );
    = | where thvm2 /= thvm1;
      |
      | - ( grav / Cp ) * ( z2 - z1 ) / thvm; where thvm2 = thvm1 (= thvm).

    The value of pressure (exner) can be calculated using the above equation
    at all vertical levels once the value of pressure is known at one level.
    Since the surface pressure is known at the initial time, that allows
    pressure to be calculated for the rest of the vertical profile.

    References:
    """
    thvm = jnp.asarray(thvm, dtype=jnp.float64)
    p_sfc = jnp.asarray(p_sfc, dtype=jnp.float64)

    # Interpolate theta_v to momentum levels.
    thvm_zm = zt2zm(
        gr.nzm,
        gr.nzt,
        gr.ngrdcol,
        gr,
        thvm,
        zero_threshold,
    )

    # The pressure (and exner) at the lowest momentum level is the surface
    # pressure (and exner based on the surface pressure).
    exner_sfc = (p_sfc / p0) ** kappa

    # Calculate exner at thermodynamic level 1 (first thermodynamic grid level
    # that is above the lower boundary).
    exner_first = _exner_step(
        exner_sfc,
        gr.zm[:, 0],
        gr.zt[:, 0],
        thvm_zm[:, 0],
        thvm[:, 0],
    )

    # Loop over all other thermodynamic vertical grid levels.
    def zt_step(exner_prev, k):
        # Calculate exner on thermodynamic levels.
        exner_next = _exner_step(
            exner_prev,
            gr.zt[:, k - 1],
            gr.zt[:, k],
            thvm[:, k - 1],
            thvm[:, k],
        )
        return exner_next, exner_next

    _, exner_rest = jax.lax.scan(
        zt_step,
        exner_first,
        jnp.arange(1, gr.nzt),
    )
    exner = jnp.concatenate([exner_first[None, :,], exner_rest], axis=0).T

    # Calculate pressure on thermodynamic levels.
    p_in_Pa = p0 * exner ** (1.0 / kappa)

    exner_zm_first = exner_sfc[:, None]

    # Loop over all momentum grid levels.
    # Calculate exner on momentum levels.
    exner_zm_rest = _exner_step(
        exner,
        gr.zt,
        gr.zm[:, 1:],
        thvm,
        thvm_zm[:, 1:],
    )
    exner_zm = jnp.concatenate([exner_zm_first, exner_zm_rest], axis=1)

    # Calculate pressure on momentum levels.
    p_in_Pa_zm = p0 * exner_zm ** (1.0 / kappa)

    return p_in_Pa, exner, p_in_Pa_zm, exner_zm


@partial(jax.jit, static_argnames=("nzt", "ngrdcol"))
def calculate_thvm(
    nzt: int,
    ngrdcol: int,
    thlm,
    rtm,
    rcm,
    exner,
    thv_ds_zt,
):
    """Calculate mean theta_v.

    Description:
    Calculates mean theta_v based on a linearized approximation to the theta_v
    equation.  This equation also includes liquid water loading.

    References:
    """
    del nzt, ngrdcol

    # ---------------------------- Begin Code ----------------------------

    # Calculate mean theta_v
    thvm = (
        thlm
        + ep1 * thv_ds_zt * rtm
        + (Lv / (Cp * exner) - ep2 * thv_ds_zt) * rcm
    )

    return thvm
