"""Test adapters for independently supplied thermodynamic/PDF inputs."""
from clubb_jax.src.Microphys.KK_microphys_module import KK_tendency_coefs
from clubb_jax.src.Microphys.KK_microphys import parameters_KK


def kk_auto_coef(rho):
    return KK_tendency_coefs(290., 1., 1.e5, rho, 1)[1]


def kk_evap_coef(T_liq, p_in_Pa, C_evap, saturation_formula=3):
    return KK_tendency_coefs(T_liq, 1., p_in_Pa, 1., saturation_formula)[0] * C_evap / parameters_KK.C_evap
