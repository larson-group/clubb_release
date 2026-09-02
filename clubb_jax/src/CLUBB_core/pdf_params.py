"""JAX-side representations from `pdf_parameter_module.F90`."""

from __future__ import annotations

from typing import NamedTuple

import jax
import jax.numpy as jnp


_IMPLICIT_COEFS_STATIC_FIELDS = (
    "ngrdcol",
    "nz",
    "sclr_dim",
)

_IMPLICIT_COEFS_ARRAY_FIELDS = (
    "coef_wp4_implicit",
    "coef_wp2rtp_implicit",
    "term_wp2rtp_explicit",
    "coef_wp2thlp_implicit",
    "term_wp2thlp_explicit",
    "coef_wp2up_implicit",
    "term_wp2up_explicit",
    "coef_wp2vp_implicit",
    "term_wp2vp_explicit",
    "coef_wprtp2_implicit",
    "term_wprtp2_explicit",
    "coef_wpthlp2_implicit",
    "term_wpthlp2_explicit",
    "coef_wprtpthlp_implicit",
    "term_wprtpthlp_explicit",
    "coef_wpup2_implicit",
    "term_wpup2_explicit",
    "coef_wpvp2_implicit",
    "term_wpvp2_explicit",
    "coef_wp2sclrp_implicit",
    "term_wp2sclrp_explicit",
    "coef_wpsclrp2_implicit",
    "term_wpsclrp2_explicit",
    "coef_wprtpsclrp_implicit",
    "term_wprtpsclrp_explicit",
    "coef_wpthlpsclrp_implicit",
    "term_wpthlpsclrp_explicit",
)


@jax.tree_util.register_pytree_node_class
class implicit_coefs_terms(NamedTuple):
    """Mirror of Fortran `type(implicit_coefs_terms)`."""

    ngrdcol: int
    nz: int
    sclr_dim: int

    coef_wp4_implicit: jnp.ndarray
    coef_wp2rtp_implicit: jnp.ndarray
    term_wp2rtp_explicit: jnp.ndarray
    coef_wp2thlp_implicit: jnp.ndarray
    term_wp2thlp_explicit: jnp.ndarray
    coef_wp2up_implicit: jnp.ndarray
    term_wp2up_explicit: jnp.ndarray
    coef_wp2vp_implicit: jnp.ndarray
    term_wp2vp_explicit: jnp.ndarray
    coef_wprtp2_implicit: jnp.ndarray
    term_wprtp2_explicit: jnp.ndarray
    coef_wpthlp2_implicit: jnp.ndarray
    term_wpthlp2_explicit: jnp.ndarray
    coef_wprtpthlp_implicit: jnp.ndarray
    term_wprtpthlp_explicit: jnp.ndarray
    coef_wpup2_implicit: jnp.ndarray
    term_wpup2_explicit: jnp.ndarray
    coef_wpvp2_implicit: jnp.ndarray
    term_wpvp2_explicit: jnp.ndarray

    coef_wp2sclrp_implicit: jnp.ndarray | None = None
    term_wp2sclrp_explicit: jnp.ndarray | None = None
    coef_wpsclrp2_implicit: jnp.ndarray | None = None
    term_wpsclrp2_explicit: jnp.ndarray | None = None
    coef_wprtpsclrp_implicit: jnp.ndarray | None = None
    term_wprtpsclrp_explicit: jnp.ndarray | None = None
    coef_wpthlpsclrp_implicit: jnp.ndarray | None = None
    term_wpthlpsclrp_explicit: jnp.ndarray | None = None

    def replace(self, **kwargs):
        return self._replace(**kwargs)

    def tree_flatten(self):
        children = tuple(
            getattr(self, name) for name in _IMPLICIT_COEFS_ARRAY_FIELDS
        )
        aux_data = tuple(
            getattr(self, name) for name in _IMPLICIT_COEFS_STATIC_FIELDS
        )
        return children, aux_data

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        data = dict(zip(_IMPLICIT_COEFS_STATIC_FIELDS, aux_data))
        data.update(dict(zip(_IMPLICIT_COEFS_ARRAY_FIELDS, children)))
        return cls(**data)


class pdf_parameter(NamedTuple):
    """Mirror of Fortran `type(pdf_parameter)`."""

    ngrdcol: int
    nz: int

    w_1: jnp.ndarray
    w_2: jnp.ndarray
    varnce_w_1: jnp.ndarray
    varnce_w_2: jnp.ndarray
    rt_1: jnp.ndarray
    rt_2: jnp.ndarray
    varnce_rt_1: jnp.ndarray
    varnce_rt_2: jnp.ndarray
    thl_1: jnp.ndarray
    thl_2: jnp.ndarray
    varnce_thl_1: jnp.ndarray
    varnce_thl_2: jnp.ndarray
    corr_w_rt_1: jnp.ndarray
    corr_w_rt_2: jnp.ndarray
    corr_w_thl_1: jnp.ndarray
    corr_w_thl_2: jnp.ndarray
    corr_rt_thl_1: jnp.ndarray
    corr_rt_thl_2: jnp.ndarray
    alpha_thl: jnp.ndarray
    alpha_rt: jnp.ndarray
    crt_1: jnp.ndarray
    crt_2: jnp.ndarray
    cthl_1: jnp.ndarray
    cthl_2: jnp.ndarray
    chi_1: jnp.ndarray
    chi_2: jnp.ndarray
    stdev_chi_1: jnp.ndarray
    stdev_chi_2: jnp.ndarray
    stdev_eta_1: jnp.ndarray
    stdev_eta_2: jnp.ndarray
    covar_chi_eta_1: jnp.ndarray
    covar_chi_eta_2: jnp.ndarray
    corr_w_chi_1: jnp.ndarray
    corr_w_chi_2: jnp.ndarray
    corr_w_eta_1: jnp.ndarray
    corr_w_eta_2: jnp.ndarray
    corr_chi_eta_1: jnp.ndarray
    corr_chi_eta_2: jnp.ndarray
    rsatl_1: jnp.ndarray
    rsatl_2: jnp.ndarray
    rc_1: jnp.ndarray
    rc_2: jnp.ndarray
    cloud_frac_1: jnp.ndarray
    cloud_frac_2: jnp.ndarray
    mixt_frac: jnp.ndarray
    ice_supersat_frac_1: jnp.ndarray
    ice_supersat_frac_2: jnp.ndarray

    def replace(self, **kwargs):
        return self._replace(**kwargs)


def init_pdf_implicit_coefs_terms_api(nz: int, ngrdcol: int, sclr_dim: int) -> implicit_coefs_terms:
    """JAX equivalent of pdf_parameter_module.F90:init_pdf_implicit_coefs_terms_api."""
    zero_2d = jnp.zeros((ngrdcol, nz), dtype=jnp.float64)

    if sclr_dim > 0:
        zero_3d = jnp.zeros((ngrdcol, nz, sclr_dim), dtype=jnp.float64)
        coef_wp2sclrp_implicit = zero_3d
        term_wp2sclrp_explicit = zero_3d
        coef_wpsclrp2_implicit = zero_3d
        term_wpsclrp2_explicit = zero_3d
        coef_wprtpsclrp_implicit = zero_3d
        term_wprtpsclrp_explicit = zero_3d
        coef_wpthlpsclrp_implicit = zero_3d
        term_wpthlpsclrp_explicit = zero_3d
    else:
        coef_wp2sclrp_implicit = None
        term_wp2sclrp_explicit = None
        coef_wpsclrp2_implicit = None
        term_wpsclrp2_explicit = None
        coef_wprtpsclrp_implicit = None
        term_wprtpsclrp_explicit = None
        coef_wpthlpsclrp_implicit = None
        term_wpthlpsclrp_explicit = None

    return implicit_coefs_terms(
        ngrdcol=ngrdcol,
        nz=nz,
        sclr_dim=sclr_dim,
        coef_wp4_implicit=zero_2d,
        coef_wp2rtp_implicit=zero_2d,
        term_wp2rtp_explicit=zero_2d,
        coef_wp2thlp_implicit=zero_2d,
        term_wp2thlp_explicit=zero_2d,
        coef_wp2up_implicit=zero_2d,
        term_wp2up_explicit=zero_2d,
        coef_wp2vp_implicit=zero_2d,
        term_wp2vp_explicit=zero_2d,
        coef_wprtp2_implicit=zero_2d,
        term_wprtp2_explicit=zero_2d,
        coef_wpthlp2_implicit=zero_2d,
        term_wpthlp2_explicit=zero_2d,
        coef_wprtpthlp_implicit=zero_2d,
        term_wprtpthlp_explicit=zero_2d,
        coef_wpup2_implicit=zero_2d,
        term_wpup2_explicit=zero_2d,
        coef_wpvp2_implicit=zero_2d,
        term_wpvp2_explicit=zero_2d,
        coef_wp2sclrp_implicit=coef_wp2sclrp_implicit,
        term_wp2sclrp_explicit=term_wp2sclrp_explicit,
        coef_wpsclrp2_implicit=coef_wpsclrp2_implicit,
        term_wpsclrp2_explicit=term_wpsclrp2_explicit,
        coef_wprtpsclrp_implicit=coef_wprtpsclrp_implicit,
        term_wprtpsclrp_explicit=term_wprtpsclrp_explicit,
        coef_wpthlpsclrp_implicit=coef_wpthlpsclrp_implicit,
        term_wpthlpsclrp_explicit=term_wpthlpsclrp_explicit,
    )


def zero_pdf_implicit_coefs_terms_api(coefs: implicit_coefs_terms) -> implicit_coefs_terms:
    """Return a zeroed implicit-coefficient container with the same dimensions."""
    return init_pdf_implicit_coefs_terms_api(coefs.nz, coefs.ngrdcol, coefs.sclr_dim)


def init_pdf_params(nz: int, ngrdcol: int) -> pdf_parameter:
    """JAX equivalent of pdf_parameter_module.F90:init_pdf_params."""
    zero_2d = jnp.zeros((ngrdcol, nz), dtype=jnp.float64)
    values = {
        field: zero_2d
        for field in pdf_parameter._fields
        if field not in ("ngrdcol", "nz")
    }
    return pdf_parameter(ngrdcol=ngrdcol, nz=nz, **values)


def zero_pdf_params_api(params: pdf_parameter) -> pdf_parameter:
    """Return a zeroed PDF-parameter container with the same dimensions."""
    return init_pdf_params(params.nz, params.ngrdcol)


def pack_pdf_params_api(params: pdf_parameter, nz: int):
    """Pack the first column of a pdf_parameter into a ``(nz, 47)`` array."""
    slots = []
    for field in pdf_parameter._fields[2:]:
        value = jnp.asarray(getattr(params, field), dtype=jnp.float64)
        if value.ndim == 1:
            slots.append(value[:nz])
        else:
            slots.append(value[0, :nz])
    return jnp.stack(slots, axis=1)


def unpack_pdf_params_api(r_param_array, nz: int, template: pdf_parameter) -> pdf_parameter:
    """Unpack a ``(nz, 47)`` array into the first column of a pdf_parameter."""
    r_param_array = jnp.asarray(r_param_array, dtype=jnp.float64)
    values = {}
    for idx, field in enumerate(pdf_parameter._fields[2:]):
        arr = jnp.zeros((template.ngrdcol, nz), dtype=jnp.float64)
        arr = arr.at[0, :].set(r_param_array[:nz, idx])
        values[field] = arr
    return pdf_parameter(ngrdcol=template.ngrdcol, nz=nz, **values)


__all__ = [
    "implicit_coefs_terms",
    "pdf_parameter",
    "init_pdf_implicit_coefs_terms_api",
    "zero_pdf_implicit_coefs_terms_api",
    "init_pdf_params",
    "zero_pdf_params_api",
    "pack_pdf_params_api",
    "unpack_pdf_params_api",
]
