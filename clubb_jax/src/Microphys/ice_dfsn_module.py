"""JAX port of Microphys/ice_dfsn_module.F90 — depletion of cloud water by diffusional growth of ice.

Mirrors `ice_dfsn` (Larson et al. 2006; Rogers & Yau 1989, Eq. 9.4): in mixed-phase cloud below freezing,
ice grows by vapor diffusion at the expense of liquid, so `rcm` is depleted and `thlm` warmed. The single ice
crystal's mass is integrated DOWNWARD as it falls (`mass(k-1) = mass(k) + dmass(k)`), a strictly sequential
vertical recurrence — ported as a top-to-bottom `lax.scan`. The thermodynamic factors (saturation, S_i, the
diffusion denominator) do not depend on the carried mass, so they are precomputed vectorized; only the mass
integration and the mass-dependent rates live in the scan.

Single-column (1D, length nzt) to mirror the Fortran exactly (which hardcodes `gr%invrs_dzm(1,...)`); a
multi-column driver vmaps over columns. Pure jnp → differentiable. Validated in `tests/test_ice_dfsn.py`
against a literal NumPy transcription (rel ~1e-14), conservation of the rcm/thlm tendency coupling, the
in-cloud/below-freezing branch, the over-depletion cap, and a finite `jax.grad`.
"""
import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import (
    Cp, Lv, ep, Rv, Lf, Ls, T_freeze_K, cm_per_m)
from clubb_jax.src.CLUBB_core.T_in_K_module import thlm2T_in_K
from clubb_jax.src.CLUBB_core.saturation import sat_mixrat_liq

# Constant parameters (ice_dfsn_module.F90)
_N_I = 2000.0            # Number of ice crystals per unit volume of air [m^-3]
_MASS_INIT = 1.0e-11     # Initial ice particle mass at model top         [kg]
_RCM_THRESHOLD = 1.0e-5  # Min rcm to be considered "in cloud"            [kg/kg]
# Mitchell (1996) mass-diameter: mass = a_coef * (diam/1m)^b_expn
_A_COEF, _B_EXPN = 2.05e-3, 1.8
# Mitchell (1996) fallspeed-diameter: u_T = k_u * rho^-q * (diam/1m)^n
_K_U_COEF, _Q_EXPN, _N_EXPN = 55.0, 0.17, 0.70


def Diff_denom(T_in_K, p_in_Pa, e_i):
    """Denominator of the diffusional-growth equation (ice_dfsn_module.F90:Diff_denom; R&Y Eq. 9.4) [m s/kg]."""
    Celsius = T_in_K - T_freeze_K
    Ka = (5.69 + 0.017 * Celsius) * 0.00001          # cal/(cm s C)
    Ka = 4.1868 * 100.0 * Ka                          # J/(m s K)
    Dv = 0.221 * (T_in_K / T_freeze_K) ** 1.94 * (101325.0 / p_in_Pa)  # cm^2/s
    Dv = Dv / 10000.0                                 # m^2/s
    Fk = (Ls / (Rv * T_in_K) - 1.0) * Ls / (Ka * T_in_K)
    Fd = (Rv * T_in_K) / (Dv * e_i)
    return Fk + Fd


def ice_dfsn(gr, dt, thlm, rcm, exner, p_in_Pa, rho, saturation_formula):
    """Time tendencies of rcm and thlm from ice diffusional growth.

    Args:
        gr:                 JAX grid (uses gr.invrs_dzm, single column, shape (1, nzm)).
        dt:                 Model timestep [s].
        thlm, rcm, exner, p_in_Pa, rho: 1D arrays, length nzt (thermodynamic grid).
        saturation_formula: SATURATION_* integer for sat_mixrat_liq.

    Returns:
        (rcm_icedfsn, thlm_icedfsn): 1D tendencies, each length nzt [kg/kg/s], [K/s].
    """
    thlm = jnp.asarray(thlm)
    rcm = jnp.asarray(rcm)
    exner = jnp.asarray(exner)
    p_in_Pa = jnp.asarray(p_in_Pa)
    rho = jnp.asarray(rho)
    nzt = thlm.shape[0]

    # --- Vectorized thermodynamic factors (mass-independent) ---
    T_in_K = thlm2T_in_K(thlm, exner, rcm)
    in_cloud = (rcm >= _RCM_THRESHOLD) & (T_in_K < T_freeze_K)
    r_s = sat_mixrat_liq(p_in_Pa, T_in_K, saturation_formula)
    e_s = (r_s * p_in_Pa) / (ep + r_s)
    e_i = e_s / jnp.exp((Lf / (Rv * T_freeze_K)) * (T_freeze_K / T_in_K - 1.0))
    S_i = e_s / e_i
    Denom = Diff_denom(T_in_K, p_in_Pa, e_i)
    factor = 4.0 * (S_i - 1.0) / Denom   # common 4*(S_i-1)/Denom term

    # Lagged momentum-grid spacing used by dmass: Fortran gr%invrs_dzm(1,k-1) -> 0-based [k-2] = [j-1].
    inv_dzm = gr.invrs_dzm[0]
    lag_idx = jnp.clip(jnp.arange(nzt) - 1, 0, inv_dzm.shape[0] - 1)
    inv_dzm_lag = inv_dzm[lag_idx]

    # --- Sequential downward mass integration (top j=nzt-1 -> bottom j=0) ---
    def step(mass, inp):
        cloud, fac, rho_j, rcm_j, inv_dzm_j = inp
        base = mass / _A_COEF                                   # > 0 (mass grows monotonically)
        rate = -(_N_I / rho_j) * fac * base ** (1.0 / _B_EXPN)
        rcm_ice = jnp.where(cloud, rate, 0.0)
        # Ensure liquid is not over-depleted.
        rcm_ice = jnp.where(cloud & (rcm_j + rcm_ice * dt < 0.0), -rcm_j / dt, rcm_ice)
        dmass = (fac * (1.0 / _K_U_COEF) * rho_j ** _Q_EXPN
                 * base ** ((1.0 - _N_EXPN) / _B_EXPN) * (1.0 / inv_dzm_j))
        diam = jnp.where(cloud, base ** (1.0 / _B_EXPN), 0.0)
        u_T_cm = jnp.where(cloud, cm_per_m * _K_U_COEF * base ** (_N_EXPN / _B_EXPN)
                           * rho_j ** (-_Q_EXPN), 0.0)
        next_mass = jnp.where(cloud, mass + dmass, mass)
        return next_mass, (rcm_ice, mass, diam, u_T_cm)

    # Reverse so scan runs top -> bottom, then un-reverse the outputs.
    rev = lambda a: a[::-1]
    xs = (rev(in_cloud), rev(factor), rev(rho), rev(rcm), rev(inv_dzm_lag))
    _, (rcm_ice_r, mass_r, diam_r, u_T_cm_r) = jax.lax.scan(step, _MASS_INIT, xs)
    rcm_icedfsn = rev(rcm_ice_r)

    # thlm tendency (ice_dfsn_module.F90:305)
    thlm_icedfsn = -(Lv / (Cp * exner)) * rcm_icedfsn
    return rcm_icedfsn, thlm_icedfsn
