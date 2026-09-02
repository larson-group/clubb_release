"""JAX port of ``CLUBB_core/hydromet_pdf_parameter_module.F90``.

Description:
This module defines the derived type hydromet_pdf_parameter.

References:
  None

Porting deviations:
  * Fortran derived types are represented as frozen dataclasses that are JAX
    pytrees.
  * Fortran allocates and mutates ``precipitation_fractions``.  JAX returns new
    dataclass instances with fixed-shape arrays.
  * Fortran ``zero_precip_fracs_api`` mutates its inout argument.  JAX
    ``zero_precip_fracs`` returns a zeroed copy preserving dimensions.
"""
from dataclasses import dataclass

import jax

from clubb_jax.src.CLUBB_core.clubb_precision import configure_jax_precision
configure_jax_precision()
import jax.numpy as jnp

MAX_HYDROMET_DIM = 8


_HYDROMET_PDF_PARAMETER_FIELDS = (
    "hm_1",
    "hm_2",
    "mu_hm_1",
    "mu_hm_2",
    "sigma_hm_1",
    "sigma_hm_2",
    "corr_w_hm_1",
    "corr_w_hm_2",
    "corr_chi_hm_1",
    "corr_chi_hm_2",
    "corr_eta_hm_1",
    "corr_eta_hm_2",
    "corr_hmx_hmy_1",
    "corr_hmx_hmy_2",
    "mu_Ncn_1",
    "mu_Ncn_2",
    "sigma_Ncn_1",
    "sigma_Ncn_2",
)


@jax.tree_util.register_pytree_node_class
@dataclass(frozen=True)
class HydrometPdfParameter:
    """JAX representation of Fortran ``type hydromet_pdf_parameter``."""
    hm_1: jnp.ndarray          # Mean of hydrometeor, hm (1st PDF component)   [un vary]
    hm_2: jnp.ndarray          # Mean of hydrometeor, hm (2nd PDF component)   [un vary]
    mu_hm_1: jnp.ndarray       # Mean of hm (1st PDF component) in-precip (ip) [un vary]
    mu_hm_2: jnp.ndarray       # Mean of hm (2nd PDF component) ip             [un vary]
    sigma_hm_1: jnp.ndarray    # Standard deviation of hm (1st PDF comp.) ip   [un vary]
    sigma_hm_2: jnp.ndarray    # Standard deviation of hm (2nd PDF comp.) ip   [un vary]
    corr_w_hm_1: jnp.ndarray   # Correlation of w and hm (1st PDF component) ip      [-]
    corr_w_hm_2: jnp.ndarray   # Correlation of w and hm (2nd PDF component) ip      [-]
    corr_chi_hm_1: jnp.ndarray # Correlation of chi and hm (1st PDF component) ip    [-]
    corr_chi_hm_2: jnp.ndarray # Correlation of chi and hm (2nd PDF component) ip    [-]
    corr_eta_hm_1: jnp.ndarray # Correlation of eta and hm (1st PDF component) ip    [-]
    corr_eta_hm_2: jnp.ndarray # Correlation of eta and hm (2nd PDF component) ip    [-]
    corr_hmx_hmy_1: jnp.ndarray # Correlation of hmx and hmy (1st PDF component) ip  [-]
    corr_hmx_hmy_2: jnp.ndarray # Correlation of hmx and hmy (2nd PDF component) ip  [-]
    mu_Ncn_1: float    # Mean of Ncn (1st PDF component)                  [num/kg]
    mu_Ncn_2: float    # Mean of Ncn (2nd PDF component)                  [num/kg]
    sigma_Ncn_1: float # Standard deviation of Ncn (1st PDF component)    [num/kg]
    sigma_Ncn_2: float # Standard deviation of Ncn (2nd PDF component)    [num/kg]

    def tree_flatten(self):
        return tuple(getattr(self, name) for name in _HYDROMET_PDF_PARAMETER_FIELDS), None

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        del aux_data
        return cls(**dict(zip(_HYDROMET_PDF_PARAMETER_FIELDS, children)))


@jax.tree_util.register_pytree_node_class
@dataclass(frozen=True)
class PrecipitationFractions:
    """JAX representation of Fortran ``type precipitation_fractions``."""
    ngrdcol: int # Dimensions of variables
    nzt: int
    precip_frac: jnp.ndarray   # Precipitation fraction (overall)           [-]
    precip_frac_1: jnp.ndarray # Precipitation fraction (1st PDF component) [-]
    precip_frac_2: jnp.ndarray # Precipitation fraction (2nd PDF component) [-]

    def tree_flatten(self):
        children = (self.precip_frac, self.precip_frac_1, self.precip_frac_2)
        aux_data = (self.ngrdcol, self.nzt)
        return children, aux_data

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        ngrdcol, nzt = aux_data
        precip_frac, precip_frac_1, precip_frac_2 = children
        return cls(
            ngrdcol=ngrdcol,
            nzt=nzt,
            precip_frac=precip_frac,
            precip_frac_1=precip_frac_1,
            precip_frac_2=precip_frac_2,
        )


def init_hydromet_pdf_params():
    """Initialize the elements of hydromet_pdf_params."""
    vec = lambda: jnp.zeros(MAX_HYDROMET_DIM, dtype=jnp.float64)
    mat = lambda: jnp.zeros((MAX_HYDROMET_DIM, MAX_HYDROMET_DIM), dtype=jnp.float64)
    # Initialize hydromet_pdf_params.
    return HydrometPdfParameter(
        hm_1=vec(), hm_2=vec(), mu_hm_1=vec(), mu_hm_2=vec(),
        sigma_hm_1=vec(), sigma_hm_2=vec(), corr_w_hm_1=vec(), corr_w_hm_2=vec(),
        corr_chi_hm_1=vec(), corr_chi_hm_2=vec(), corr_eta_hm_1=vec(), corr_eta_hm_2=vec(),
        corr_hmx_hmy_1=mat(), corr_hmx_hmy_2=mat(),
        mu_Ncn_1=0.0, mu_Ncn_2=0.0, sigma_Ncn_1=0.0, sigma_Ncn_2=0.0)


def init_precip_fracs(nzt, ngrdcol):
    """Initialize the elements of precip_fracs."""
    # Allocate precip frac arrays
    z = jnp.zeros((ngrdcol, nzt), dtype=jnp.float64)
    # Set metadata values
    # Zero precip_fracs
    return PrecipitationFractions(ngrdcol=ngrdcol, nzt=nzt,
                                  precip_frac=z, precip_frac_1=z, precip_frac_2=z)


def zero_precip_fracs(precip_fracs):
    """Initialize the elements of precip_fracs."""
    # Initialize precip_fracs.
    z = jnp.zeros((precip_fracs.ngrdcol, precip_fracs.nzt), dtype=jnp.float64)
    return PrecipitationFractions(ngrdcol=precip_fracs.ngrdcol, nzt=precip_fracs.nzt,
                                  precip_frac=z, precip_frac_1=z, precip_frac_2=z)
