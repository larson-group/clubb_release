"""JAX port of `src/CLUBB_core/T_in_K_module.F90`.

Porting deviations:
- The Fortran forward conversion accepts full column/level batches; the
  inverse is elemental. JAX expresses both formulas with array operations,
  taking their dimensions from the supplied fields.
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

    JAX adaptation: array operations apply the formula to every column and level.
    """
    # ---------------------- Begin Code ----------------------
    T_in_K = thlm * exner + Lv * rcm / Cp
    return T_in_K


def T_in_K2thlm(T_in_K, exner, rcm):
    """Calculates liquid water potential temperature from absolute temperature.

    Fortran comments:
      Description:
        Calculates liquid water potential temperature from absolute temperature

      References:
        None
    """
    # ---- Begin Code ----
    thlm = (T_in_K - Lv / Cp * rcm) / exner
    return thlm
