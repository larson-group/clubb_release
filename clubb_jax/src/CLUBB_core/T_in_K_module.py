"""JAX port of `src/CLUBB_core/T_in_K_module.F90`.

Porting deviations:
- The Fortran module exposes scalar, 1D, and 2D ``thlm2T_in_K_api``
  procedures plus an elemental ``T_in_K2thlm_api``.  The JAX functions use one
  array-broadcasting implementation for each formula because JAX arrays cover
  the scalar and array cases without separate overloads.
"""

from __future__ import annotations

from clubb_jax.src.CLUBB_core.constants_clubb import Cp, Lv


def thlm2T_in_K(thlm, exner, rcm):
    """Calculates absolute temperature from liquid water potential temperature.

    Fortran comments:
      Description:
        Calculates absolute temperature from liquid water potential
        temperature.  (Does not include ice.)

      References:
        Cotton and Anthes (1989), "Storm and Cloud Dynamics", Eqn. (2.51).

    JAX adaptation: the Fortran scalar, 1D, and 2D overloads share this
    broadcasting implementation.
    """
    # ---------------------- Begin Code ----------------------
    return thlm * exner + Lv * rcm / Cp


def T_in_K2thlm(T_in_K, exner, rcm):
    """Calculates liquid water potential temperature from absolute temperature.

    Fortran comments:
      Description:
        Calculates liquid water potential temperature from absolute temperature

      References:
        None
    """
    # ---- Begin Code ----
    return (T_in_K - Lv / Cp * rcm) / exner
