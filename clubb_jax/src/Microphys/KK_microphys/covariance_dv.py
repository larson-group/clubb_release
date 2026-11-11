"""Shared parabolic-cylinder values for the KK covariance driver.

Fortran evaluates Dv inside each integral. Here we batch the requests before
assembling the covariances so XLA does not inline the special-function algorithm
at every use. All values are local to this invocation and remain differentiable.
The standalone mean-tendency path does not use this batch.
"""
from typing import NamedTuple

import jax
import jax.numpy as jnp

from clubb_jax.src.CLUBB_core.constants_clubb import chi_tol
from clubb_jax.src.Microphys.KK_microphys.parameters_KK import (
    KK_auto_rc_exp, KK_auto_Nc_exp, KK_accr_rc_exp, KK_accr_rr_exp,
    KK_evap_Supersat_exp, KK_evap_rr_exp, KK_evap_Nr_exp,
)
from clubb_jax.src.Microphys.KK_microphys.parabolic_cylinder import _dvc


class CovarianceDv(NamedTuple):
    # Each covariance variant contains (D_{-(alpha+1)}, D_{-(alpha+2)}).
    covar: dict
    # Mean integrals at powers alpha and alpha+1 use different sigma guards.
    mean: tuple


def batch_covariance_dv(mu_chi, sigma_chi, sigma_Ncn_n, corr_chi_Ncn_n,
                        sigma_rr_n, corr_chi_rr_n, sigma_Nr_n, corr_chi_Nr_n):
    """Return (auto, accr, evap), each containing values for PDF components 1/2.

Inputs are pairs of component arrays (or scalars). The 64 requests cover both
orders, the reduced-distribution variants, and the distinct covariance/mean
denominator guards. They are shared across w, rt, and thl without identifying
requests by their runtime values or changing the integral arithmetic.
"""
    orders, arguments = [], []

    def request(order, argument):
        index = len(orders)
        orders.append(order)
        arguments.append(argument)
        return index

    def component(i, alpha, beta, sigma_y_n, corr_chi_y_n,
                  gamma=None, sigma_z_n=None, corr_chi_z_n=None):
        variants = []
        for guard, denominator in enumerate((
                jnp.where(sigma_chi[i] > chi_tol, sigma_chi[i], 1.0),
                jnp.maximum(sigma_chi[i], chi_tol))):
            r = mu_chi[i] / denominator
            if gamma is None:
                args = dict(full=-(r + corr_chi_y_n[i] * sigma_y_n[i] * beta),
                            const_y=-r)
            else:
                args = dict(
                    full=(r + corr_chi_y_n[i] * sigma_y_n[i] * beta
                          + corr_chi_z_n[i] * sigma_z_n[i] * gamma),
                    const_y=r + corr_chi_z_n[i] * sigma_z_n[i] * gamma,
                    const_z=r + corr_chi_y_n[i] * sigma_y_n[i] * beta,
                    const_yz=r)
            # Preserve the mean's (alpha+1)+1 evaluation order as well.
            upper_order = alpha + 2.0 if guard == 0 else (alpha + 1.0) + 1.0
            variants.append({name: (request(-(alpha + 1.0), arg),
                                    request(-upper_order, arg))
                             for name, arg in args.items()})
        covar, mean = variants
        return CovarianceDv(covar, tuple({name: pair[k] for name, pair in mean.items()}
                                        for k in range(2)))

    auto = tuple(component(i, KK_auto_rc_exp, KK_auto_Nc_exp,
                           sigma_Ncn_n, corr_chi_Ncn_n) for i in range(2))
    accr = tuple(component(i, KK_accr_rc_exp, KK_accr_rr_exp,
                           sigma_rr_n, corr_chi_rr_n) for i in range(2))
    evap = tuple(component(i, KK_evap_Supersat_exp, KK_evap_rr_exp,
                           sigma_rr_n, corr_chi_rr_n,
                           KK_evap_Nr_exp, sigma_Nr_n, corr_chi_Nr_n) for i in range(2))
    arguments = jnp.stack(jnp.broadcast_arrays(*arguments))
    orders = jnp.asarray(orders)
    # Compile one Dv body, but evaluate each request with its original array
    # shape. A single vectorized Dv call changes CPU arithmetic in 
    # coupled runs - accumulated roundoff made DYCOMS RF02 DS fail at 1e-7.
    values = jax.lax.map(lambda request: _dvc(*request), (orders, arguments))
    return jax.tree.map(lambda index: values[index], (auto, accr, evap))
