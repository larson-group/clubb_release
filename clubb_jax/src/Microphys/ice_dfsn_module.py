"""JAX port of Microphys/ice_dfsn_module.F90 — depletion of cloud water by diffusional growth of ice.

Mirrors `ice_dfsn` (Larson et al. 2006; Rogers & Yau 1989, Eq. 9.4): in mixed-phase cloud below freezing,
ice grows by vapor diffusion at the expense of liquid, so `rcm` is depleted and `thlm` warmed. The single ice
crystal's mass is integrated DOWNWARD as it falls (`mass(k-1) = mass(k) + dmass(k)`), a strictly sequential
vertical recurrence — ported as a top-to-bottom `lax.scan`. The thermodynamic factors (saturation, S_i, the
diffusion denominator) do not depend on the carried mass, so they are precomputed vectorized; only the mass
integration and the mass-dependent rates live in the scan.

The last axis is vertical; leading column axes are batched. Grid spacing is taken independently for each column. Stats and output tendencies are returned in order. Pure jnp → differentiable. Validated in `tests/test_ice_dfsn.py`
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


def ice_dfsn(gr, ngrdcol, dt, thlm, rcm,
    exner, p_in_Pa, rho, saturation_formula, stats):
    """Time tendencies of rcm and thlm from ice diffusional growth.

    Args:
        gr:                 JAX grid with gr.invrs_dzm shaped (ngrdcol, nzm).
        dt:                 Model timestep [s].
        thlm, rcm, exner, p_in_Pa, rho: Arrays shaped (ngrdcol, nzt) on the thermodynamic grid.
        saturation_formula: SATURATION_* integer for sat_mixrat_liq.

    Returns:
        (stats, rcm_icedfsn, thlm_icedfsn): updated statistics and tendencies
        shaped (ngrdcol, nzt), in [kg/kg/s] and [K/s].
    """
    # Description:
    #   This subroutine is based on a COAMPS subroutine (nov11_icedfs)
    #   written by Adam Smith and Vince Larson to calculate the
    #   depletion of cloud water by the diffusional growth of ice.
    #
    #---------------Brian's comment--------------------------------------!
    # This code does not use actual microphysics.  Diffusional growth of !
    # ice is supposed to be the growth of ice due to diffusion of water  !
    # vapor.  Liquid water is not involved in diffusional growth.        !
    # However, in mixed phase clouds (both ice and liquid water), most   !
    # of the water vapor condenses onto the liquid droplets due to the   !
    # fact that they have so much more available surface area.  This     !
    # brings the amount of water vapor in the atmosphere to the          !
    # saturation level with respect to liquid water.  However, since the !
    # saturation vapor pressure with respect to ice is less than the     !
    # saturation vapor pressure with respect to liquid water, a          !
    # saturated atmosphere with respect to liquid water is still         !
    # supersaturated with respect to ice.  As a result, ice still grows  !
    # due to diffusion.  When this happens, the environmental vapor      !
    # pressure drops to the point of saturation with respect to ice.     !
    # This leaves the atmosphere subsaturated with respect to liquid     !
    # water.  As a result, some of the liquid water evaporates until     !
    # the atmosphere becomes saturated with respect to liquid water      !
    # again.  The process then repeats itself.  As a result, the ice     !
    # essentially grows at the expense of the liquid water.  This is     !
    # why the diffusional growth of ice is being deducted from liquid    !
    # water in this subroutine.
    #-------------------------------------------------------------------------------
    # References:
    #   Section 4.2 of Larson et al. (2006), "What determines altocumulus
    #     dissipation time?", J. Geophys. Res., Vol. 111, D19207.
    #
    #   Mitchell, D. L. (1996), "Use of mass- and area- ...", J. Atmos. Sci.
    #     Vol. 53, 1710--1723.
    #
    #   Rogers and Yau (1989), "A Short Course in Cloud Physics", 3rd. Ed.
    #
    #   Fleishauer et al. (2002), "Observed microphysical structure of
    #     midlevel, mixed-phase clouds", J. Atmos. Sci., Vol. 59,
    #     pp. 1779--1804.
    #-------------------------------------------------------------------------------
    thlm = jnp.asarray(thlm)
    rcm = jnp.asarray(rcm)
    exner = jnp.asarray(exner)
    p_in_Pa = jnp.asarray(p_in_Pa)
    rho = jnp.asarray(rho)
    nzt = thlm.shape[-1]

    # --- Vectorized thermodynamic factors (mass-independent) ---
    T_in_K = thlm2T_in_K(thlm, exner, rcm)
    in_cloud = (rcm >= _RCM_THRESHOLD) & (T_in_K < T_freeze_K)
    r_s = sat_mixrat_liq(p_in_Pa, T_in_K, saturation_formula)
    e_s = (r_s * p_in_Pa) / (ep + r_s)
    e_i = e_s / jnp.exp((Lf / (Rv * T_freeze_K)) * (T_freeze_K / T_in_K - 1.0))
    S_i = e_s / e_i
    Denom = Diff_denom(T_in_K, p_in_Pa, e_i)
    factor = 4.0 * (S_i - 1.0) / Denom   # common 4*(S_i-1)/Denom term

    # Lagged momentum-grid spacing used by dmass: Fortran gr%invrs_dzm(icol,k-1) -> 0-based [k-2] = [j-1].
    inv_dzm = gr.invrs_dzm
    lag_idx = jnp.clip(jnp.arange(nzt) - 1, 0, inv_dzm.shape[-1] - 1)
    inv_dzm_lag = jnp.broadcast_to(inv_dzm[:, lag_idx], thlm.shape)

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

    # JAX scan-boundary adapters: the recurrent vertical axis must be first and
    # run top -> bottom. These two local transforms only change scan layout;
    # they do not wrap physics or replace source routines.
    rev = lambda a: jnp.moveaxis(a[..., ::-1], -1, 0)
    restore = lambda a: jnp.moveaxis(a, 0, -1)[..., ::-1]
    xs = (rev(in_cloud), rev(factor), rev(rho), rev(rcm), rev(inv_dzm_lag))
    _, (rcm_ice_r, mass_r, diam_r, u_T_cm_r) = jax.lax.scan(step, jnp.full(thlm.shape[:-1], _MASS_INIT), xs)
    rcm_icedfsn = restore(rcm_ice_r)
    if stats.l_sample:
        stats = stats.update('rcm_icedfs', rcm_icedfsn)
        stats = stats.update('diam', restore(diam_r))
        stats = stats.update('mass_ice_cryst', restore(mass_r))
        stats = stats.update('u_T_cm', restore(u_T_cm_r))

    # thlm tendency (ice_dfsn_module.F90:305)
    thlm_icedfsn = -(Lv / (Cp * exner)) * rcm_icedfsn
    return stats, rcm_icedfsn, thlm_icedfsn


def Diff_denom(T_in_K, p_in_Pa, e_i):
    """Denominator of the diffusional-growth equation (ice_dfsn_module.F90:Diff_denom; R&Y Eq. 9.4) [m s/kg]."""
    # Description:
    #   Compute denominator of diffusional growth equation
    #
    # References:
    #   Eqn. 9.4 of Rogers and Yau (1989), "A Short Course on Cloud Physics"
    #
    #-----------------------------------------------------------------------------
    # Reference:  Eqn. 9.4 of Rogers and Yau (1989), "A Short Course on Cloud Physics"
    # Constant Parameters
    #   real, parameter :: Ls = 2.834e6
    Celsius = T_in_K - T_freeze_K
    Ka = (5.69 + 0.017 * Celsius) * 0.00001          # cal/(cm s C)
    Ka = 4.1868 * 100.0 * Ka                          # J/(m s K)
    Dv = 0.221 * (T_in_K / T_freeze_K) ** 1.94 * (101325.0 / p_in_Pa)  # cm^2/s
    Dv = Dv / 10000.0                                 # m^2/s
    Fk = (Ls / (Rv * T_in_K) - 1.0) * Ls / (Ka * T_in_K)
    Fd = (Rv * T_in_K) / (Dv * e_i)
    return Fk + Fd
