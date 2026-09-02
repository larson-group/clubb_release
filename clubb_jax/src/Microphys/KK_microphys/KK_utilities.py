"""JAX port of KK_utilities.F90 — thermodynamic helpers for the KK microphysics.

So far: G_T_p, the drop-radius-growth coefficient G(T,p) [m^2/s] used in the KK
evaporation coefficient (KK_evap_coef = 3 C_evap G_T_p). From Rogers & Yau (1989)
Eqs. 7.17-7.18 / Khairoutdinov & Kogan (2000) Eq. 22. Oracle: KK_utilities.F90:169.

(The parabolic cylinder function Dv_fnc, also in KK_utilities.F90, lives in
parabolic_cylinder.py.)
"""
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()

from clubb_jax.src.CLUBB_core.saturation import (
    sat_vapor_press_liq, SATURATION_FLATAU,
)
from clubb_jax.src.CLUBB_core.constants_clubb import (
    Lv, Rv, T_freeze_K, rho_lw,
)


def G_T_p(T_in_K, p_in_Pa, saturation_formula=SATURATION_FLATAU):
    """G(T,p) [m^2/s], the drop-radius growth coefficient. KK_utilities.F90:169.

    Ka = thermal conductivity of air; Dv = vapor diffusion coefficient;
    Fk (heat-conduction term) and Fd (vapor-diffusion term) per Rogers & Yau 7.17-7.18;
    G = 1 / (Fk + Fd)."""
    T_in_K = jnp.asarray(T_in_K, dtype=jnp.float64)
    p_in_Pa = jnp.asarray(p_in_Pa, dtype=jnp.float64)
    Celsius = T_in_K - T_freeze_K

    # Thermal conductivity of air: cal/(cm s C) -> J/(m s K).
    Ka = (5.69 + 0.017 * Celsius) * 0.00001
    Ka = 4.1868 * 100.0 * Ka

    # Vapor diffusion coefficient: cm^2/s -> m^2/s.
    Dv = 0.221 * (T_in_K / 273.16) ** 1.94 * (101325.0 / p_in_Pa)
    Dv = Dv / 10000.0

    esatv = sat_vapor_press_liq(T_in_K, saturation_formula)

    Fk = (Lv / (Rv * T_in_K) - 1.0) * (Lv * rho_lw) / (Ka * T_in_K)
    Fd = (rho_lw * Rv * T_in_K) / (Dv * esatv)
    return 1.0 / (Fk + Fd)
