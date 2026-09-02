"""JAX implementations of selected routines from ``Skx_module.F90``.

Porting deviations:
Fortran subroutines mutate output arrays and loop over explicit ``nz`` and
``ngrdcol`` dimensions; these JAX functions return arrays and broadcast over
the provided profile shapes.  ``Skx_func`` omits the disabled
``l_clipping_kluge = .false.`` block because the Fortran comments state that
clipping is no longer needed and the block is inactive in the source.
"""

from __future__ import annotations

from functools import partial

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import eps, one, one_half, w_tol_sqd
from clubb_jax.src.CLUBB_core.parameter_indices import (
    iSkw_denom_coef,
    ibeta,
    igamma_coef,
    igamma_coefb,
    igamma_coefc,
)


@partial(jax.jit, static_argnames=("nz", "ngrdcol"))
def Skx_func(nz: int, ngrdcol: int, xp2, xp3, x_tol: float, clubb_params):
    """Calculate the skewness of x.

    Description:
    Calculate the skewness of x

    References:
    None
    """
    del nz, ngrdcol
    # ---- Begin Code ----
    Skx_denom_tol = jnp.asarray(clubb_params)[:, iSkw_denom_coef, None] * x_tol ** 2

    # Calculation of skewness to help reduce the sensitivity of this value to
    # small values of xp2.
    return jnp.asarray(xp3) * jnp.sqrt(jnp.asarray(xp2) + Skx_denom_tol) ** (-3)


@partial(jax.jit, static_argnames=("nz", "ngrdcol", "l_gamma_Skw"))
def compute_gamma_Skw(nz: int, ngrdcol: int, Skw, clubb_params, l_gamma_Skw: bool):
    """Compute gamma as a function of Skw on a single vertical grid.

    Description:
      Compute gamma as a function of Skw on a single vertical grid.

    References:
      None
    """
    del nz, ngrdcol
    Skw = jnp.asarray(Skw, dtype=jnp.float64)
    clubb_params = jnp.asarray(clubb_params, dtype=jnp.float64)

    # First gamma(Skw) coefficient [-]
    gamma_coef = clubb_params[:, igamma_coef, None]
    # Second gamma(Skw) coefficient [-]
    gamma_coefb = clubb_params[:, igamma_coefb, None]
    # Third gamma(Skw) coefficient [-]
    gamma_coefc = clubb_params[:, igamma_coefc, None]

    #----------------------------- Begin Code ------------------------------
    if l_gamma_Skw:
        varying = gamma_coefb + (gamma_coef - gamma_coefb) * jnp.exp(
            -one_half * (Skw / gamma_coefc) ** 2
        )
        l_varying = (
            jnp.abs(gamma_coef - gamma_coefb)
            > jnp.abs(gamma_coef + gamma_coefb) * eps / 2.0
        )
        gamma_Skw_fnc = jnp.where(l_varying, varying, gamma_coef)
    else:
        gamma_Skw_fnc = gamma_coef + jnp.zeros_like(Skw)

    return gamma_Skw_fnc


@partial(jax.jit, static_argnames=("nz", "ngrdcol"))
def LG_2005_ansatz(
    nz: int,
    ngrdcol: int,
    Skw,
    wpxp,
    wp2,
    xp2,
    beta,
    sigma_sqd_w,
    x_tol: float,
):
    """Calculate the skewness of x using the diagnostic ansatz of Larson and Golaz (2005).

    Description:
    Calculate the skewness of x using the diagnostic ansatz of Larson and
    Golaz (2005).

    References:
    Vincent E. Larson and Jean-Christophe Golaz, 2005:  Using Probability
    Density Functions to Derive Consistent Closure Relationships among
    Higher-Order Moments.  Mon. Wea. Rev., 133, 1023–1042.
    """
    del nz, ngrdcol
    Skw = jnp.asarray(Skw, dtype=jnp.float64)
    wpxp = jnp.asarray(wpxp, dtype=jnp.float64)
    wp2 = jnp.asarray(wp2, dtype=jnp.float64)
    xp2 = jnp.asarray(xp2, dtype=jnp.float64)
    beta = jnp.asarray(beta, dtype=jnp.float64)
    sigma_sqd_w = jnp.asarray(sigma_sqd_w, dtype=jnp.float64)

    #--------------------------Begin Code --------------------------

    # weberjk, 8-July 2015. Commented this out for now. cgils was failing during some tests.

    # Larson and Golaz (2005) eq. 16
    one_minus_sigma_sqd_w = one - sigma_sqd_w
    nrmlzd_corr_wx = wpxp / jnp.sqrt(
        jnp.maximum(wp2, w_tol_sqd)
        * jnp.maximum(xp2, x_tol ** 2)
        * one_minus_sigma_sqd_w
    )

    # Larson and Golaz (2005) eq. 11
    nrmlzd_Skw = Skw / (one_minus_sigma_sqd_w * jnp.sqrt(one_minus_sigma_sqd_w))

    # Larson and Golaz (2005) eq. 33
    return nrmlzd_Skw * nrmlzd_corr_wx * (
        beta[:, None] + (one - beta[:, None]) * nrmlzd_corr_wx ** 2
    )


@partial(jax.jit, static_argnames=("nzt", "ngrdcol"))
def xp3_LG_2005_ansatz(
    nzt: int,
    ngrdcol: int,
    Skw_zt,
    wpxp_zt,
    wp2_zt,
    xp2_zt,
    sigma_sqd_w_zt,
    clubb_params,
    x_tol: float,
):
    """Calculate <x'^3> after calculating the skewness of x using the ansatz of Larson and Golaz (2005).

    Description:
    Calculate <x'^3> after calculating the skewness of x using the ansatz of
    Larson and Golaz (2005).

    References:
    """
    del ngrdcol
    clubb_params = jnp.asarray(clubb_params, dtype=jnp.float64)
    Skx_denom_tol = clubb_params[:, iSkw_denom_coef, None] * x_tol ** 2
    beta = clubb_params[:, ibeta]

    #-------------------------- Begin Code --------------------------

    # Calculate skewness of x using the ansatz of LG05.
    Skx_zt = LG_2005_ansatz(
        nzt,
        clubb_params.shape[0],
        Skw_zt,
        wpxp_zt,
        wp2_zt,
        xp2_zt,
        beta,
        sigma_sqd_w_zt,
        x_tol,
    )

    # Calculate <x'^3> using the reverse of the special sensitivity reduction
    # formula in function Skx_func above.
    return Skx_zt * (jnp.asarray(xp2_zt) + Skx_denom_tol) * jnp.sqrt(
        jnp.asarray(xp2_zt) + Skx_denom_tol
    )
