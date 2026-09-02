"""DISCLAIMER: this code appears to be correct but has not been
very thouroughly tested. If you do notice any anomalous behaviour then please
contact Andy and/or Bjorn.

Function diag_ustar returns the value of ustar using the below similarity
functions and a specified buoyancy flux (bflx) given in kinematic units.

phi_m (zeta > 0) =  (1 + am * zeta)
phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)

where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar).

Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface Layer, in
Workshop on Micormeteorology, pages 67-100.

Code writen March, 1999 by Bjorn Stevens.
"""

from __future__ import annotations

from functools import partial

import jax
import jax.numpy as jnp
from jax import lax

from clubb_jax.src.CLUBB_core.constants_clubb import pi, vonk


@partial(jax.jit)
def diag_ustar(z, bflx, wnd, z0):
    """Return the value of ustar using the module similarity functions."""
    am = 4.8
    bm = 19.3

    z, bflx, wnd, z0 = jnp.broadcast_arrays(z, bflx, wnd, z0)
    shape = z.shape

    # TODO(port-mirror): The source routine is scalar. This local vmap is the
    # JAX array-language expression of the caller's Fortran column loop; remove
    # it if all callers move their column vmap outside this routine.
    def column(z_col, bflx_col, wnd_col, z0_col):
        lnz = jnp.log(z_col / z0_col)
        klnz = vonk / lnz
        c1 = pi / 2.0 - 3.0 * jnp.log(2.0)

        ustar = wnd_col * klnz

        def apply_similarity(ustar_initial):
            def iterate(_, carry):
                ustar_iter, finished = carry

                def update(_):
                    lmo = -ustar_iter**3 / (vonk * bflx_col)
                    zeta = z_col / lmo

                    def stable(_):
                        return lax.cond(
                            zeta > 1.0e10,
                            lambda __: (jnp.asarray(1.0e-10), jnp.asarray(True)),
                            lambda __: (vonk * wnd_col / (lnz + am * zeta), jnp.asarray(False)),
                            operand=None,
                        )

                    def unstable(_):
                        x = jnp.sqrt(jnp.sqrt(1.0 - bm * zeta))
                        psi1 = (
                            2.0 * jnp.log(1.0 + x)
                            + jnp.log(1.0 + x * x)
                            - 2.0 * jnp.arctan(x)
                            + c1
                        )
                        return wnd_col * vonk / (lnz - psi1), jnp.asarray(False)

                    return lax.cond(zeta > 0.0, stable, unstable, operand=None)

                return lax.cond(
                    finished,
                    lambda _: (ustar_iter, finished),
                    update,
                    operand=None,
                )

            return lax.fori_loop(0, 4, iterate, (ustar_initial, jnp.asarray(False)))[0]

        return lax.cond(
            jnp.abs(bflx_col) > 1.0e-6,
            apply_similarity,
            lambda value: value,
            ustar,
        )

    return jax.vmap(column)(z.ravel(), bflx.ravel(), wnd.ravel(), z0.ravel()).reshape(shape)


__all__ = ["diag_ustar"]
