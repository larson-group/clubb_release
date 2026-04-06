"""Python type mirrors and helper routines for `pdf_parameter_module.F90`."""

from __future__ import annotations

import sys
from typing import NamedTuple, Optional, TextIO

import numpy as np

from clubb_python.derived_types.common import Array


class pdf_parameter(NamedTuple):
    """Mirror of Fortran `type(pdf_parameter)`."""

    ngrdcol: int
    nz: int

    w_1: Array
    w_2: Array
    varnce_w_1: Array
    varnce_w_2: Array
    rt_1: Array
    rt_2: Array
    varnce_rt_1: Array
    varnce_rt_2: Array
    thl_1: Array
    thl_2: Array
    varnce_thl_1: Array
    varnce_thl_2: Array
    corr_w_rt_1: Array
    corr_w_rt_2: Array
    corr_w_thl_1: Array
    corr_w_thl_2: Array
    corr_rt_thl_1: Array
    corr_rt_thl_2: Array
    alpha_thl: Array
    alpha_rt: Array
    crt_1: Array
    crt_2: Array
    cthl_1: Array
    cthl_2: Array
    chi_1: Array
    chi_2: Array
    stdev_chi_1: Array
    stdev_chi_2: Array
    stdev_eta_1: Array
    stdev_eta_2: Array
    covar_chi_eta_1: Array
    covar_chi_eta_2: Array
    corr_w_chi_1: Array
    corr_w_chi_2: Array
    corr_w_eta_1: Array
    corr_w_eta_2: Array
    corr_chi_eta_1: Array
    corr_chi_eta_2: Array
    rsatl_1: Array
    rsatl_2: Array
    rc_1: Array
    rc_2: Array
    cloud_frac_1: Array
    cloud_frac_2: Array
    mixt_frac: Array
    ice_supersat_frac_1: Array
    ice_supersat_frac_2: Array


class implicit_coefs_terms(NamedTuple):
    """Mirror of Fortran `type(implicit_coefs_terms)`."""

    ngrdcol: int
    nz: int
    sclr_dim: int

    coef_wp4_implicit: Array
    coef_wp2rtp_implicit: Array
    term_wp2rtp_explicit: Array
    coef_wp2thlp_implicit: Array
    term_wp2thlp_explicit: Array
    coef_wp2up_implicit: Array
    term_wp2up_explicit: Array
    coef_wp2vp_implicit: Array
    term_wp2vp_explicit: Array
    coef_wprtp2_implicit: Array
    term_wprtp2_explicit: Array
    coef_wpthlp2_implicit: Array
    term_wpthlp2_explicit: Array
    coef_wprtpthlp_implicit: Array
    term_wprtpthlp_explicit: Array
    coef_wpup2_implicit: Array
    term_wpup2_explicit: Array
    coef_wpvp2_implicit: Array
    term_wpvp2_explicit: Array

    coef_wp2sclrp_implicit: Optional[Array] = None
    term_wp2sclrp_explicit: Optional[Array] = None
    coef_wpsclrp2_implicit: Optional[Array] = None
    term_wpsclrp2_explicit: Optional[Array] = None
    coef_wprtpsclrp_implicit: Optional[Array] = None
    term_wprtpsclrp_explicit: Optional[Array] = None
    coef_wpthlpsclrp_implicit: Optional[Array] = None
    term_wpthlpsclrp_explicit: Optional[Array] = None


def _zeros_2d(ngrdcol: int, nz: int) -> Array:
    return np.zeros((ngrdcol, nz), dtype=np.float64, order='F')


def _zeros_3d(ngrdcol: int, nz: int, sclr_dim: int) -> Array:
    return np.zeros((ngrdcol, nz, sclr_dim), dtype=np.float64, order='F')


def _as_2d(name: str, arr, ngrdcol: int, nz: int) -> Array:
    out = np.asfortranarray(arr, dtype=np.float64)
    if out.shape != (ngrdcol, nz):
        raise ValueError(f"{name} shape mismatch: expected {(ngrdcol, nz)}, got {out.shape}")
    return out


def _as_3d(name: str, arr, ngrdcol: int, nz: int, sclr_dim: int) -> Array:
    out = np.asfortranarray(arr, dtype=np.float64)
    if out.shape != (ngrdcol, nz, sclr_dim):
        raise ValueError(
            f"{name} shape mismatch: expected {(ngrdcol, nz, sclr_dim)}, got {out.shape}"
        )
    return out


def init_pdf_params(nz: int, ngrdcol: int) -> pdf_parameter:
    """Python equivalent of Fortran `init_pdf_params` (alloc + zero)."""
    return pdf_parameter(
        ngrdcol=ngrdcol,
        nz=nz,
        w_1=_zeros_2d(ngrdcol, nz),
        w_2=_zeros_2d(ngrdcol, nz),
        varnce_w_1=_zeros_2d(ngrdcol, nz),
        varnce_w_2=_zeros_2d(ngrdcol, nz),
        rt_1=_zeros_2d(ngrdcol, nz),
        rt_2=_zeros_2d(ngrdcol, nz),
        varnce_rt_1=_zeros_2d(ngrdcol, nz),
        varnce_rt_2=_zeros_2d(ngrdcol, nz),
        thl_1=_zeros_2d(ngrdcol, nz),
        thl_2=_zeros_2d(ngrdcol, nz),
        varnce_thl_1=_zeros_2d(ngrdcol, nz),
        varnce_thl_2=_zeros_2d(ngrdcol, nz),
        corr_w_rt_1=_zeros_2d(ngrdcol, nz),
        corr_w_rt_2=_zeros_2d(ngrdcol, nz),
        corr_w_thl_1=_zeros_2d(ngrdcol, nz),
        corr_w_thl_2=_zeros_2d(ngrdcol, nz),
        corr_rt_thl_1=_zeros_2d(ngrdcol, nz),
        corr_rt_thl_2=_zeros_2d(ngrdcol, nz),
        alpha_thl=_zeros_2d(ngrdcol, nz),
        alpha_rt=_zeros_2d(ngrdcol, nz),
        crt_1=_zeros_2d(ngrdcol, nz),
        crt_2=_zeros_2d(ngrdcol, nz),
        cthl_1=_zeros_2d(ngrdcol, nz),
        cthl_2=_zeros_2d(ngrdcol, nz),
        chi_1=_zeros_2d(ngrdcol, nz),
        chi_2=_zeros_2d(ngrdcol, nz),
        stdev_chi_1=_zeros_2d(ngrdcol, nz),
        stdev_chi_2=_zeros_2d(ngrdcol, nz),
        stdev_eta_1=_zeros_2d(ngrdcol, nz),
        stdev_eta_2=_zeros_2d(ngrdcol, nz),
        covar_chi_eta_1=_zeros_2d(ngrdcol, nz),
        covar_chi_eta_2=_zeros_2d(ngrdcol, nz),
        corr_w_chi_1=_zeros_2d(ngrdcol, nz),
        corr_w_chi_2=_zeros_2d(ngrdcol, nz),
        corr_w_eta_1=_zeros_2d(ngrdcol, nz),
        corr_w_eta_2=_zeros_2d(ngrdcol, nz),
        corr_chi_eta_1=_zeros_2d(ngrdcol, nz),
        corr_chi_eta_2=_zeros_2d(ngrdcol, nz),
        rsatl_1=_zeros_2d(ngrdcol, nz),
        rsatl_2=_zeros_2d(ngrdcol, nz),
        rc_1=_zeros_2d(ngrdcol, nz),
        rc_2=_zeros_2d(ngrdcol, nz),
        cloud_frac_1=_zeros_2d(ngrdcol, nz),
        cloud_frac_2=_zeros_2d(ngrdcol, nz),
        mixt_frac=_zeros_2d(ngrdcol, nz),
        ice_supersat_frac_1=_zeros_2d(ngrdcol, nz),
        ice_supersat_frac_2=_zeros_2d(ngrdcol, nz),
    )


def zero_pdf_params_api(pdf_params: pdf_parameter) -> pdf_parameter:
    """Python equivalent of Fortran `zero_pdf_params_api`."""
    return init_pdf_params(pdf_params.nz, pdf_params.ngrdcol)


def init_pdf_implicit_coefs_terms_api(
    nz: int,
    ngrdcol: int,
    sclr_dim: int,
) -> implicit_coefs_terms:
    """Python equivalent of Fortran `init_pdf_implicit_coefs_terms_api`."""
    if sclr_dim > 0:
        coef_wp2sclrp_implicit = _zeros_3d(ngrdcol, nz, sclr_dim)
        term_wp2sclrp_explicit = _zeros_3d(ngrdcol, nz, sclr_dim)
        coef_wpsclrp2_implicit = _zeros_3d(ngrdcol, nz, sclr_dim)
        term_wpsclrp2_explicit = _zeros_3d(ngrdcol, nz, sclr_dim)
        coef_wprtpsclrp_implicit = _zeros_3d(ngrdcol, nz, sclr_dim)
        term_wprtpsclrp_explicit = _zeros_3d(ngrdcol, nz, sclr_dim)
        coef_wpthlpsclrp_implicit = _zeros_3d(ngrdcol, nz, sclr_dim)
        term_wpthlpsclrp_explicit = _zeros_3d(ngrdcol, nz, sclr_dim)
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
        coef_wp4_implicit=_zeros_2d(ngrdcol, nz),
        coef_wp2rtp_implicit=_zeros_2d(ngrdcol, nz),
        term_wp2rtp_explicit=_zeros_2d(ngrdcol, nz),
        coef_wp2thlp_implicit=_zeros_2d(ngrdcol, nz),
        term_wp2thlp_explicit=_zeros_2d(ngrdcol, nz),
        coef_wp2up_implicit=_zeros_2d(ngrdcol, nz),
        term_wp2up_explicit=_zeros_2d(ngrdcol, nz),
        coef_wp2vp_implicit=_zeros_2d(ngrdcol, nz),
        term_wp2vp_explicit=_zeros_2d(ngrdcol, nz),
        coef_wprtp2_implicit=_zeros_2d(ngrdcol, nz),
        term_wprtp2_explicit=_zeros_2d(ngrdcol, nz),
        coef_wpthlp2_implicit=_zeros_2d(ngrdcol, nz),
        term_wpthlp2_explicit=_zeros_2d(ngrdcol, nz),
        coef_wprtpthlp_implicit=_zeros_2d(ngrdcol, nz),
        term_wprtpthlp_explicit=_zeros_2d(ngrdcol, nz),
        coef_wpup2_implicit=_zeros_2d(ngrdcol, nz),
        term_wpup2_explicit=_zeros_2d(ngrdcol, nz),
        coef_wpvp2_implicit=_zeros_2d(ngrdcol, nz),
        term_wpvp2_explicit=_zeros_2d(ngrdcol, nz),
        coef_wp2sclrp_implicit=coef_wp2sclrp_implicit,
        term_wp2sclrp_explicit=term_wp2sclrp_explicit,
        coef_wpsclrp2_implicit=coef_wpsclrp2_implicit,
        term_wpsclrp2_explicit=term_wpsclrp2_explicit,
        coef_wprtpsclrp_implicit=coef_wprtpsclrp_implicit,
        term_wprtpsclrp_explicit=term_wprtpsclrp_explicit,
        coef_wpthlpsclrp_implicit=coef_wpthlpsclrp_implicit,
        term_wpthlpsclrp_explicit=term_wpthlpsclrp_explicit,
    )


def zero_pdf_implicit_coefs_terms_api(
    pdf_implicit_coefs_terms: implicit_coefs_terms,
) -> implicit_coefs_terms:
    """Python equivalent of Fortran `zero_pdf_implicit_coefs_terms_api`."""
    return init_pdf_implicit_coefs_terms_api(
        pdf_implicit_coefs_terms.nz,
        pdf_implicit_coefs_terms.ngrdcol,
        pdf_implicit_coefs_terms.sclr_dim,
    )


def print_pdf_params(pdf_params: pdf_parameter, ngrdcol: int, stream: TextIO | None = None) -> None:
    """Python equivalent of Fortran `print_pdf_params`."""
    out = sys.stderr if stream is None else stream
    ncol = min(int(ngrdcol), int(pdf_params.ngrdcol))
    for i in range(ncol):
        print(f"PDF parameters for column {i + 1}", file=out)
        print("---------------------------------------------------", file=out)
        print("w_1(i,:) =", pdf_params.w_1[i, :], file=out)
        print("w_2(i,:) =", pdf_params.w_2[i, :], file=out)
        print("varnce_w_1(i,:) =", pdf_params.varnce_w_1[i, :], file=out)
        print("varnce_w_2(i,:) =", pdf_params.varnce_w_2[i, :], file=out)
        print("rt_1(i,:) =", pdf_params.rt_1[i, :], file=out)
        print("rt_2(i,:) =", pdf_params.rt_2[i, :], file=out)
        print("varnce_rt_1(i,:) =", pdf_params.varnce_rt_1[i, :], file=out)
        print("varnce_rt_2(i,:) =", pdf_params.varnce_rt_2[i, :], file=out)
        print("thl_1(i,:) =", pdf_params.thl_1[i, :], file=out)
        print("thl_2(i,:) =", pdf_params.thl_2[i, :], file=out)
        print("varnce_thl_1(i,:) =", pdf_params.varnce_thl_1[i, :], file=out)
        print("varnce_thl_2(i,:) =", pdf_params.varnce_thl_2[i, :], file=out)
        print("corr_w_rt_1(i,:) =", pdf_params.corr_w_rt_1[i, :], file=out)
        print("corr_w_rt_2(i,:) =", pdf_params.corr_w_rt_2[i, :], file=out)
        print("corr_w_thl_1(i,:) =", pdf_params.corr_w_thl_1[i, :], file=out)
        print("corr_w_thl_2(i,:) =", pdf_params.corr_w_thl_2[i, :], file=out)
        print("corr_rt_thl_1(i,:) =", pdf_params.corr_rt_thl_1[i, :], file=out)
        print("corr_rt_thl_2(i,:) =", pdf_params.corr_rt_thl_2[i, :], file=out)
        print("alpha_thl(i,:) =", pdf_params.alpha_thl[i, :], file=out)
        print("alpha_rt(i,:) =", pdf_params.alpha_rt[i, :], file=out)
        print("crt_1(i,:) =", pdf_params.crt_1[i, :], file=out)
        print("crt_2(i,:) =", pdf_params.crt_2[i, :], file=out)
        print("cthl_1(i,:) =", pdf_params.cthl_1[i, :], file=out)
        print("cthl_2(i,:) =", pdf_params.cthl_2[i, :], file=out)
        print("chi_1(i,:) =", pdf_params.chi_1[i, :], file=out)
        print("chi_2(i,:) =", pdf_params.chi_2[i, :], file=out)
        print("stdev_chi_1(i,:) =", pdf_params.stdev_chi_1[i, :], file=out)
        print("stdev_chi_2(i,:) =", pdf_params.stdev_chi_2[i, :], file=out)
        print("stdev_eta_1(i,:) =", pdf_params.stdev_eta_1[i, :], file=out)
        print("stdev_eta_2(i,:) =", pdf_params.stdev_eta_2[i, :], file=out)
        print("covar_chi_eta_1(i,:) =", pdf_params.covar_chi_eta_1[i, :], file=out)
        print("covar_chi_eta_2(i,:) =", pdf_params.covar_chi_eta_2[i, :], file=out)
        print("corr_w_chi_1(i,:) =", pdf_params.corr_w_chi_1[i, :], file=out)
        print("corr_w_chi_2(i,:) =", pdf_params.corr_w_chi_2[i, :], file=out)
        print("corr_w_eta_1(i,:) =", pdf_params.corr_w_eta_1[i, :], file=out)
        print("corr_w_eta_2(i,:) =", pdf_params.corr_w_eta_2[i, :], file=out)
        print("corr_chi_eta_1(i,:) =", pdf_params.corr_chi_eta_1[i, :], file=out)
        print("corr_chi_eta_2(i,:) =", pdf_params.corr_chi_eta_2[i, :], file=out)
        print("rsatl_1(i,:) =", pdf_params.rsatl_1[i, :], file=out)
        print("rsatl_2(i,:) =", pdf_params.rsatl_2[i, :], file=out)
        print("rc_1(i,:) =", pdf_params.rc_1[i, :], file=out)
        print("rc_2(i,:) =", pdf_params.rc_2[i, :], file=out)
        print("cloud_frac_1(i,:) =", pdf_params.cloud_frac_1[i, :], file=out)
        print("cloud_frac_2(i,:) =", pdf_params.cloud_frac_2[i, :], file=out)
        print("mixt_frac(i,:) =", pdf_params.mixt_frac[i, :], file=out)
        print("ice_supersat_frac_1(i,:) =", pdf_params.ice_supersat_frac_1[i, :], file=out)
        print("ice_supersat_frac_2(i,:) =", pdf_params.ice_supersat_frac_2[i, :], file=out)


def pack_pdf_params(pdf_params: pdf_parameter) -> Array:
    """Pack `pdf_parameter` to shape `(ngrdcol, nz, 47)` for Fortran wrapper APIs."""
    ncol = int(pdf_params.ngrdcol)
    nz = int(pdf_params.nz)
    packed = np.empty((ncol, nz, 47), dtype=np.float64, order='F')

    packed[:, :, 0] = _as_2d('w_1', pdf_params.w_1, ncol, nz)
    packed[:, :, 1] = _as_2d('w_2', pdf_params.w_2, ncol, nz)
    packed[:, :, 2] = _as_2d('varnce_w_1', pdf_params.varnce_w_1, ncol, nz)
    packed[:, :, 3] = _as_2d('varnce_w_2', pdf_params.varnce_w_2, ncol, nz)
    packed[:, :, 4] = _as_2d('rt_1', pdf_params.rt_1, ncol, nz)
    packed[:, :, 5] = _as_2d('rt_2', pdf_params.rt_2, ncol, nz)
    packed[:, :, 6] = _as_2d('varnce_rt_1', pdf_params.varnce_rt_1, ncol, nz)
    packed[:, :, 7] = _as_2d('varnce_rt_2', pdf_params.varnce_rt_2, ncol, nz)
    packed[:, :, 8] = _as_2d('thl_1', pdf_params.thl_1, ncol, nz)
    packed[:, :, 9] = _as_2d('thl_2', pdf_params.thl_2, ncol, nz)
    packed[:, :, 10] = _as_2d('varnce_thl_1', pdf_params.varnce_thl_1, ncol, nz)
    packed[:, :, 11] = _as_2d('varnce_thl_2', pdf_params.varnce_thl_2, ncol, nz)
    packed[:, :, 12] = _as_2d('corr_w_rt_1', pdf_params.corr_w_rt_1, ncol, nz)
    packed[:, :, 13] = _as_2d('corr_w_rt_2', pdf_params.corr_w_rt_2, ncol, nz)
    packed[:, :, 14] = _as_2d('corr_w_thl_1', pdf_params.corr_w_thl_1, ncol, nz)
    packed[:, :, 15] = _as_2d('corr_w_thl_2', pdf_params.corr_w_thl_2, ncol, nz)
    packed[:, :, 16] = _as_2d('corr_rt_thl_1', pdf_params.corr_rt_thl_1, ncol, nz)
    packed[:, :, 17] = _as_2d('corr_rt_thl_2', pdf_params.corr_rt_thl_2, ncol, nz)
    packed[:, :, 18] = _as_2d('alpha_thl', pdf_params.alpha_thl, ncol, nz)
    packed[:, :, 19] = _as_2d('alpha_rt', pdf_params.alpha_rt, ncol, nz)
    packed[:, :, 20] = _as_2d('crt_1', pdf_params.crt_1, ncol, nz)
    packed[:, :, 21] = _as_2d('crt_2', pdf_params.crt_2, ncol, nz)
    packed[:, :, 22] = _as_2d('cthl_1', pdf_params.cthl_1, ncol, nz)
    packed[:, :, 23] = _as_2d('cthl_2', pdf_params.cthl_2, ncol, nz)
    packed[:, :, 24] = _as_2d('chi_1', pdf_params.chi_1, ncol, nz)
    packed[:, :, 25] = _as_2d('chi_2', pdf_params.chi_2, ncol, nz)
    packed[:, :, 26] = _as_2d('stdev_chi_1', pdf_params.stdev_chi_1, ncol, nz)
    packed[:, :, 27] = _as_2d('stdev_chi_2', pdf_params.stdev_chi_2, ncol, nz)
    packed[:, :, 28] = _as_2d('stdev_eta_1', pdf_params.stdev_eta_1, ncol, nz)
    packed[:, :, 29] = _as_2d('stdev_eta_2', pdf_params.stdev_eta_2, ncol, nz)
    packed[:, :, 30] = _as_2d('covar_chi_eta_1', pdf_params.covar_chi_eta_1, ncol, nz)
    packed[:, :, 31] = _as_2d('covar_chi_eta_2', pdf_params.covar_chi_eta_2, ncol, nz)
    packed[:, :, 32] = _as_2d('corr_w_chi_1', pdf_params.corr_w_chi_1, ncol, nz)
    packed[:, :, 33] = _as_2d('corr_w_chi_2', pdf_params.corr_w_chi_2, ncol, nz)
    packed[:, :, 34] = _as_2d('corr_w_eta_1', pdf_params.corr_w_eta_1, ncol, nz)
    packed[:, :, 35] = _as_2d('corr_w_eta_2', pdf_params.corr_w_eta_2, ncol, nz)
    packed[:, :, 36] = _as_2d('corr_chi_eta_1', pdf_params.corr_chi_eta_1, ncol, nz)
    packed[:, :, 37] = _as_2d('corr_chi_eta_2', pdf_params.corr_chi_eta_2, ncol, nz)
    packed[:, :, 38] = _as_2d('rsatl_1', pdf_params.rsatl_1, ncol, nz)
    packed[:, :, 39] = _as_2d('rsatl_2', pdf_params.rsatl_2, ncol, nz)
    packed[:, :, 40] = _as_2d('rc_1', pdf_params.rc_1, ncol, nz)
    packed[:, :, 41] = _as_2d('rc_2', pdf_params.rc_2, ncol, nz)
    packed[:, :, 42] = _as_2d('cloud_frac_1', pdf_params.cloud_frac_1, ncol, nz)
    packed[:, :, 43] = _as_2d('cloud_frac_2', pdf_params.cloud_frac_2, ncol, nz)
    packed[:, :, 44] = _as_2d('mixt_frac', pdf_params.mixt_frac, ncol, nz)
    packed[:, :, 45] = _as_2d('ice_supersat_frac_1', pdf_params.ice_supersat_frac_1, ncol, nz)
    packed[:, :, 46] = _as_2d('ice_supersat_frac_2', pdf_params.ice_supersat_frac_2, ncol, nz)

    return packed


def unpack_pdf_params(packed) -> pdf_parameter:
    """Unpack shape `(ngrdcol, nz, 47)` to `pdf_parameter`."""
    p = np.asfortranarray(packed, dtype=np.float64)
    if p.ndim != 3 or p.shape[2] != 47:
        raise ValueError(f"Expected packed pdf shape (ngrdcol, nz, 47), got {p.shape}")
    ngrdcol, nz, _ = p.shape

    return pdf_parameter(
        ngrdcol=ngrdcol,
        nz=nz,
        w_1=p[:, :, 0],
        w_2=p[:, :, 1],
        varnce_w_1=p[:, :, 2],
        varnce_w_2=p[:, :, 3],
        rt_1=p[:, :, 4],
        rt_2=p[:, :, 5],
        varnce_rt_1=p[:, :, 6],
        varnce_rt_2=p[:, :, 7],
        thl_1=p[:, :, 8],
        thl_2=p[:, :, 9],
        varnce_thl_1=p[:, :, 10],
        varnce_thl_2=p[:, :, 11],
        corr_w_rt_1=p[:, :, 12],
        corr_w_rt_2=p[:, :, 13],
        corr_w_thl_1=p[:, :, 14],
        corr_w_thl_2=p[:, :, 15],
        corr_rt_thl_1=p[:, :, 16],
        corr_rt_thl_2=p[:, :, 17],
        alpha_thl=p[:, :, 18],
        alpha_rt=p[:, :, 19],
        crt_1=p[:, :, 20],
        crt_2=p[:, :, 21],
        cthl_1=p[:, :, 22],
        cthl_2=p[:, :, 23],
        chi_1=p[:, :, 24],
        chi_2=p[:, :, 25],
        stdev_chi_1=p[:, :, 26],
        stdev_chi_2=p[:, :, 27],
        stdev_eta_1=p[:, :, 28],
        stdev_eta_2=p[:, :, 29],
        covar_chi_eta_1=p[:, :, 30],
        covar_chi_eta_2=p[:, :, 31],
        corr_w_chi_1=p[:, :, 32],
        corr_w_chi_2=p[:, :, 33],
        corr_w_eta_1=p[:, :, 34],
        corr_w_eta_2=p[:, :, 35],
        corr_chi_eta_1=p[:, :, 36],
        corr_chi_eta_2=p[:, :, 37],
        rsatl_1=p[:, :, 38],
        rsatl_2=p[:, :, 39],
        rc_1=p[:, :, 40],
        rc_2=p[:, :, 41],
        cloud_frac_1=p[:, :, 42],
        cloud_frac_2=p[:, :, 43],
        mixt_frac=p[:, :, 44],
        ice_supersat_frac_1=p[:, :, 45],
        ice_supersat_frac_2=p[:, :, 46],
    )


def pack_implicit_coefs_2d(pdf_implicit_coefs_terms: implicit_coefs_terms) -> Array:
    """Pack implicit 2D fields to shape `(ngrdcol, nz, 19)`."""
    ncol = int(pdf_implicit_coefs_terms.ngrdcol)
    nz = int(pdf_implicit_coefs_terms.nz)
    packed = np.empty((ncol, nz, 19), dtype=np.float64, order='F')

    packed[:, :, 0] = _as_2d('coef_wp4_implicit', pdf_implicit_coefs_terms.coef_wp4_implicit, ncol, nz)
    packed[:, :, 1] = _as_2d('coef_wp2rtp_implicit', pdf_implicit_coefs_terms.coef_wp2rtp_implicit, ncol, nz)
    packed[:, :, 2] = _as_2d('term_wp2rtp_explicit', pdf_implicit_coefs_terms.term_wp2rtp_explicit, ncol, nz)
    packed[:, :, 3] = _as_2d('coef_wp2thlp_implicit', pdf_implicit_coefs_terms.coef_wp2thlp_implicit, ncol, nz)
    packed[:, :, 4] = _as_2d('term_wp2thlp_explicit', pdf_implicit_coefs_terms.term_wp2thlp_explicit, ncol, nz)
    packed[:, :, 5] = _as_2d('coef_wp2up_implicit', pdf_implicit_coefs_terms.coef_wp2up_implicit, ncol, nz)
    packed[:, :, 6] = _as_2d('term_wp2up_explicit', pdf_implicit_coefs_terms.term_wp2up_explicit, ncol, nz)
    packed[:, :, 7] = _as_2d('coef_wp2vp_implicit', pdf_implicit_coefs_terms.coef_wp2vp_implicit, ncol, nz)
    packed[:, :, 8] = _as_2d('term_wp2vp_explicit', pdf_implicit_coefs_terms.term_wp2vp_explicit, ncol, nz)
    packed[:, :, 9] = _as_2d('coef_wprtp2_implicit', pdf_implicit_coefs_terms.coef_wprtp2_implicit, ncol, nz)
    packed[:, :, 10] = _as_2d('term_wprtp2_explicit', pdf_implicit_coefs_terms.term_wprtp2_explicit, ncol, nz)
    packed[:, :, 11] = _as_2d('coef_wpthlp2_implicit', pdf_implicit_coefs_terms.coef_wpthlp2_implicit, ncol, nz)
    packed[:, :, 12] = _as_2d('term_wpthlp2_explicit', pdf_implicit_coefs_terms.term_wpthlp2_explicit, ncol, nz)
    packed[:, :, 13] = _as_2d('coef_wprtpthlp_implicit', pdf_implicit_coefs_terms.coef_wprtpthlp_implicit, ncol, nz)
    packed[:, :, 14] = _as_2d('term_wprtpthlp_explicit', pdf_implicit_coefs_terms.term_wprtpthlp_explicit, ncol, nz)
    packed[:, :, 15] = _as_2d('coef_wpup2_implicit', pdf_implicit_coefs_terms.coef_wpup2_implicit, ncol, nz)
    packed[:, :, 16] = _as_2d('term_wpup2_explicit', pdf_implicit_coefs_terms.term_wpup2_explicit, ncol, nz)
    packed[:, :, 17] = _as_2d('coef_wpvp2_implicit', pdf_implicit_coefs_terms.coef_wpvp2_implicit, ncol, nz)
    packed[:, :, 18] = _as_2d('term_wpvp2_explicit', pdf_implicit_coefs_terms.term_wpvp2_explicit, ncol, nz)

    return packed


def pack_implicit_coefs_3d(pdf_implicit_coefs_terms: implicit_coefs_terms) -> Optional[Array]:
    """Pack implicit 3D fields to shape `(ngrdcol, nz, sclr_dim, 8)` or `None`."""
    if int(pdf_implicit_coefs_terms.sclr_dim) <= 0:
        return None

    ncol = int(pdf_implicit_coefs_terms.ngrdcol)
    nz = int(pdf_implicit_coefs_terms.nz)
    sclr_dim = int(pdf_implicit_coefs_terms.sclr_dim)
    packed = np.empty((ncol, nz, sclr_dim, 8), dtype=np.float64, order='F')

    if pdf_implicit_coefs_terms.coef_wp2sclrp_implicit is None or pdf_implicit_coefs_terms.term_wp2sclrp_explicit is None:
        raise ValueError("3D implicit fields are required when sclr_dim > 0")
    if pdf_implicit_coefs_terms.coef_wpsclrp2_implicit is None or pdf_implicit_coefs_terms.term_wpsclrp2_explicit is None:
        raise ValueError("3D implicit fields are required when sclr_dim > 0")
    if pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit is None or pdf_implicit_coefs_terms.term_wprtpsclrp_explicit is None:
        raise ValueError("3D implicit fields are required when sclr_dim > 0")
    if pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit is None or pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit is None:
        raise ValueError("3D implicit fields are required when sclr_dim > 0")

    packed[:, :, :, 0] = _as_3d('coef_wp2sclrp_implicit', pdf_implicit_coefs_terms.coef_wp2sclrp_implicit, ncol, nz, sclr_dim)
    packed[:, :, :, 1] = _as_3d('term_wp2sclrp_explicit', pdf_implicit_coefs_terms.term_wp2sclrp_explicit, ncol, nz, sclr_dim)
    packed[:, :, :, 2] = _as_3d('coef_wpsclrp2_implicit', pdf_implicit_coefs_terms.coef_wpsclrp2_implicit, ncol, nz, sclr_dim)
    packed[:, :, :, 3] = _as_3d('term_wpsclrp2_explicit', pdf_implicit_coefs_terms.term_wpsclrp2_explicit, ncol, nz, sclr_dim)
    packed[:, :, :, 4] = _as_3d('coef_wprtpsclrp_implicit', pdf_implicit_coefs_terms.coef_wprtpsclrp_implicit, ncol, nz, sclr_dim)
    packed[:, :, :, 5] = _as_3d('term_wprtpsclrp_explicit', pdf_implicit_coefs_terms.term_wprtpsclrp_explicit, ncol, nz, sclr_dim)
    packed[:, :, :, 6] = _as_3d('coef_wpthlpsclrp_implicit', pdf_implicit_coefs_terms.coef_wpthlpsclrp_implicit, ncol, nz, sclr_dim)
    packed[:, :, :, 7] = _as_3d('term_wpthlpsclrp_explicit', pdf_implicit_coefs_terms.term_wpthlpsclrp_explicit, ncol, nz, sclr_dim)

    return packed


def unpack_implicit_coefs(data_2d, data_3d=None) -> implicit_coefs_terms:
    """Unpack packed arrays to `implicit_coefs_terms`."""
    d2 = np.asfortranarray(data_2d, dtype=np.float64)
    if d2.ndim != 3 or d2.shape[2] != 19:
        raise ValueError(f"Expected implicit 2D packed shape (ngrdcol, nz, 19), got {d2.shape}")
    ngrdcol, nz, _ = d2.shape

    if data_3d is None:
        sclr_dim = 0
        d3 = None
    else:
        d3 = np.asfortranarray(data_3d, dtype=np.float64)
        if d3.ndim != 4 or d3.shape[0] != ngrdcol or d3.shape[1] != nz or d3.shape[3] != 8:
            raise ValueError(
                f"Expected implicit 3D packed shape (ngrdcol, nz, sclr_dim, 8), got {d3.shape}"
            )
        sclr_dim = d3.shape[2]

    return implicit_coefs_terms(
        ngrdcol=ngrdcol,
        nz=nz,
        sclr_dim=sclr_dim,
        coef_wp4_implicit=d2[:, :, 0],
        coef_wp2rtp_implicit=d2[:, :, 1],
        term_wp2rtp_explicit=d2[:, :, 2],
        coef_wp2thlp_implicit=d2[:, :, 3],
        term_wp2thlp_explicit=d2[:, :, 4],
        coef_wp2up_implicit=d2[:, :, 5],
        term_wp2up_explicit=d2[:, :, 6],
        coef_wp2vp_implicit=d2[:, :, 7],
        term_wp2vp_explicit=d2[:, :, 8],
        coef_wprtp2_implicit=d2[:, :, 9],
        term_wprtp2_explicit=d2[:, :, 10],
        coef_wpthlp2_implicit=d2[:, :, 11],
        term_wpthlp2_explicit=d2[:, :, 12],
        coef_wprtpthlp_implicit=d2[:, :, 13],
        term_wprtpthlp_explicit=d2[:, :, 14],
        coef_wpup2_implicit=d2[:, :, 15],
        term_wpup2_explicit=d2[:, :, 16],
        coef_wpvp2_implicit=d2[:, :, 17],
        term_wpvp2_explicit=d2[:, :, 18],
        coef_wp2sclrp_implicit=None if d3 is None else d3[:, :, :, 0],
        term_wp2sclrp_explicit=None if d3 is None else d3[:, :, :, 1],
        coef_wpsclrp2_implicit=None if d3 is None else d3[:, :, :, 2],
        term_wpsclrp2_explicit=None if d3 is None else d3[:, :, :, 3],
        coef_wprtpsclrp_implicit=None if d3 is None else d3[:, :, :, 4],
        term_wprtpsclrp_explicit=None if d3 is None else d3[:, :, :, 5],
        coef_wpthlpsclrp_implicit=None if d3 is None else d3[:, :, :, 6],
        term_wpthlpsclrp_explicit=None if d3 is None else d3[:, :, :, 7],
    )


def pack_pdf_params_api(
    pdf_params: pdf_parameter,
    nz: int,
    k_start_in: Optional[int] = None,
    k_end_in: Optional[int] = None,
) -> Array:
    """Python equivalent of Fortran `pack_pdf_params_api` for column 1."""
    if int(nz) != int(pdf_params.nz):
        raise ValueError(f"nz mismatch: expected {pdf_params.nz}, got {nz}")

    k_start = 0 if k_start_in is None else int(k_start_in) - 1
    k_end = int(nz) if k_end_in is None else int(k_end_in)
    if k_start < 0 or k_end > int(nz) or k_start >= k_end:
        raise ValueError(f"Invalid k range: k_start={k_start + 1}, k_end={k_end}, nz={nz}")

    r_param_array = np.empty((k_end - k_start, 47), dtype=np.float64, order='F')
    r_param_array[:, 0] = pdf_params.w_1[0, k_start:k_end]
    r_param_array[:, 1] = pdf_params.w_2[0, k_start:k_end]
    r_param_array[:, 2] = pdf_params.varnce_w_1[0, k_start:k_end]
    r_param_array[:, 3] = pdf_params.varnce_w_2[0, k_start:k_end]
    r_param_array[:, 4] = pdf_params.rt_1[0, k_start:k_end]
    r_param_array[:, 5] = pdf_params.rt_2[0, k_start:k_end]
    r_param_array[:, 6] = pdf_params.varnce_rt_1[0, k_start:k_end]
    r_param_array[:, 7] = pdf_params.varnce_rt_2[0, k_start:k_end]
    r_param_array[:, 8] = pdf_params.thl_1[0, k_start:k_end]
    r_param_array[:, 9] = pdf_params.thl_2[0, k_start:k_end]
    r_param_array[:, 10] = pdf_params.varnce_thl_1[0, k_start:k_end]
    r_param_array[:, 11] = pdf_params.varnce_thl_2[0, k_start:k_end]
    r_param_array[:, 12] = pdf_params.corr_w_rt_1[0, k_start:k_end]
    r_param_array[:, 13] = pdf_params.corr_w_rt_2[0, k_start:k_end]
    r_param_array[:, 14] = pdf_params.corr_w_thl_1[0, k_start:k_end]
    r_param_array[:, 15] = pdf_params.corr_w_thl_2[0, k_start:k_end]
    r_param_array[:, 16] = pdf_params.corr_rt_thl_1[0, k_start:k_end]
    r_param_array[:, 17] = pdf_params.corr_rt_thl_2[0, k_start:k_end]
    r_param_array[:, 18] = pdf_params.alpha_thl[0, k_start:k_end]
    r_param_array[:, 19] = pdf_params.alpha_rt[0, k_start:k_end]
    r_param_array[:, 20] = pdf_params.crt_1[0, k_start:k_end]
    r_param_array[:, 21] = pdf_params.crt_2[0, k_start:k_end]
    r_param_array[:, 22] = pdf_params.cthl_1[0, k_start:k_end]
    r_param_array[:, 23] = pdf_params.cthl_2[0, k_start:k_end]
    r_param_array[:, 24] = pdf_params.chi_1[0, k_start:k_end]
    r_param_array[:, 25] = pdf_params.chi_2[0, k_start:k_end]
    r_param_array[:, 26] = pdf_params.stdev_chi_1[0, k_start:k_end]
    r_param_array[:, 27] = pdf_params.stdev_chi_2[0, k_start:k_end]
    r_param_array[:, 28] = pdf_params.stdev_eta_1[0, k_start:k_end]
    r_param_array[:, 29] = pdf_params.stdev_eta_2[0, k_start:k_end]
    r_param_array[:, 30] = pdf_params.covar_chi_eta_1[0, k_start:k_end]
    r_param_array[:, 31] = pdf_params.covar_chi_eta_2[0, k_start:k_end]
    r_param_array[:, 32] = pdf_params.corr_w_chi_1[0, k_start:k_end]
    r_param_array[:, 33] = pdf_params.corr_w_chi_2[0, k_start:k_end]
    r_param_array[:, 34] = pdf_params.corr_w_eta_1[0, k_start:k_end]
    r_param_array[:, 35] = pdf_params.corr_w_eta_2[0, k_start:k_end]
    r_param_array[:, 36] = pdf_params.corr_chi_eta_1[0, k_start:k_end]
    r_param_array[:, 37] = pdf_params.corr_chi_eta_2[0, k_start:k_end]
    r_param_array[:, 38] = pdf_params.rsatl_1[0, k_start:k_end]
    r_param_array[:, 39] = pdf_params.rsatl_2[0, k_start:k_end]
    r_param_array[:, 40] = pdf_params.rc_1[0, k_start:k_end]
    r_param_array[:, 41] = pdf_params.rc_2[0, k_start:k_end]
    r_param_array[:, 42] = pdf_params.cloud_frac_1[0, k_start:k_end]
    r_param_array[:, 43] = pdf_params.cloud_frac_2[0, k_start:k_end]
    r_param_array[:, 44] = pdf_params.mixt_frac[0, k_start:k_end]
    r_param_array[:, 45] = pdf_params.ice_supersat_frac_1[0, k_start:k_end]
    r_param_array[:, 46] = pdf_params.ice_supersat_frac_2[0, k_start:k_end]
    return r_param_array


def unpack_pdf_params_api(
    r_param_array,
    nz: int,
    pdf_params: pdf_parameter,
    k_start_in: Optional[int] = None,
    k_end_in: Optional[int] = None,
) -> pdf_parameter:
    """Python equivalent of Fortran `unpack_pdf_params_api` for column 1."""
    data = np.asfortranarray(r_param_array, dtype=np.float64)

    if int(nz) != int(pdf_params.nz):
        raise ValueError(f"nz mismatch: expected {pdf_params.nz}, got {nz}")

    k_start = 0 if k_start_in is None else int(k_start_in) - 1
    k_end = int(nz) if k_end_in is None else int(k_end_in)
    if k_start < 0 or k_end > int(nz) or k_start >= k_end:
        raise ValueError(f"Invalid k range: k_start={k_start + 1}, k_end={k_end}, nz={nz}")

    expected_rows = k_end - k_start
    if data.shape != (expected_rows, 47):
        raise ValueError(
            f"Expected r_param_array shape {(expected_rows, 47)}, got {data.shape}"
        )

    w_1 = np.asfortranarray(pdf_params.w_1)
    w_2 = np.asfortranarray(pdf_params.w_2)
    varnce_w_1 = np.asfortranarray(pdf_params.varnce_w_1)
    varnce_w_2 = np.asfortranarray(pdf_params.varnce_w_2)
    rt_1 = np.asfortranarray(pdf_params.rt_1)
    rt_2 = np.asfortranarray(pdf_params.rt_2)
    varnce_rt_1 = np.asfortranarray(pdf_params.varnce_rt_1)
    varnce_rt_2 = np.asfortranarray(pdf_params.varnce_rt_2)
    thl_1 = np.asfortranarray(pdf_params.thl_1)
    thl_2 = np.asfortranarray(pdf_params.thl_2)
    varnce_thl_1 = np.asfortranarray(pdf_params.varnce_thl_1)
    varnce_thl_2 = np.asfortranarray(pdf_params.varnce_thl_2)
    corr_w_rt_1 = np.asfortranarray(pdf_params.corr_w_rt_1)
    corr_w_rt_2 = np.asfortranarray(pdf_params.corr_w_rt_2)
    corr_w_thl_1 = np.asfortranarray(pdf_params.corr_w_thl_1)
    corr_w_thl_2 = np.asfortranarray(pdf_params.corr_w_thl_2)
    corr_rt_thl_1 = np.asfortranarray(pdf_params.corr_rt_thl_1)
    corr_rt_thl_2 = np.asfortranarray(pdf_params.corr_rt_thl_2)
    alpha_thl = np.asfortranarray(pdf_params.alpha_thl)
    alpha_rt = np.asfortranarray(pdf_params.alpha_rt)
    crt_1 = np.asfortranarray(pdf_params.crt_1)
    crt_2 = np.asfortranarray(pdf_params.crt_2)
    cthl_1 = np.asfortranarray(pdf_params.cthl_1)
    cthl_2 = np.asfortranarray(pdf_params.cthl_2)
    chi_1 = np.asfortranarray(pdf_params.chi_1)
    chi_2 = np.asfortranarray(pdf_params.chi_2)
    stdev_chi_1 = np.asfortranarray(pdf_params.stdev_chi_1)
    stdev_chi_2 = np.asfortranarray(pdf_params.stdev_chi_2)
    stdev_eta_1 = np.asfortranarray(pdf_params.stdev_eta_1)
    stdev_eta_2 = np.asfortranarray(pdf_params.stdev_eta_2)
    covar_chi_eta_1 = np.asfortranarray(pdf_params.covar_chi_eta_1)
    covar_chi_eta_2 = np.asfortranarray(pdf_params.covar_chi_eta_2)
    corr_w_chi_1 = np.asfortranarray(pdf_params.corr_w_chi_1)
    corr_w_chi_2 = np.asfortranarray(pdf_params.corr_w_chi_2)
    corr_w_eta_1 = np.asfortranarray(pdf_params.corr_w_eta_1)
    corr_w_eta_2 = np.asfortranarray(pdf_params.corr_w_eta_2)
    corr_chi_eta_1 = np.asfortranarray(pdf_params.corr_chi_eta_1)
    corr_chi_eta_2 = np.asfortranarray(pdf_params.corr_chi_eta_2)
    rsatl_1 = np.asfortranarray(pdf_params.rsatl_1)
    rsatl_2 = np.asfortranarray(pdf_params.rsatl_2)
    rc_1 = np.asfortranarray(pdf_params.rc_1)
    rc_2 = np.asfortranarray(pdf_params.rc_2)
    cloud_frac_1 = np.asfortranarray(pdf_params.cloud_frac_1)
    cloud_frac_2 = np.asfortranarray(pdf_params.cloud_frac_2)
    mixt_frac = np.asfortranarray(pdf_params.mixt_frac)
    ice_supersat_frac_1 = np.asfortranarray(pdf_params.ice_supersat_frac_1)
    ice_supersat_frac_2 = np.asfortranarray(pdf_params.ice_supersat_frac_2)

    w_1[0, k_start:k_end] = data[:, 0]
    w_2[0, k_start:k_end] = data[:, 1]
    varnce_w_1[0, k_start:k_end] = data[:, 2]
    varnce_w_2[0, k_start:k_end] = data[:, 3]
    rt_1[0, k_start:k_end] = data[:, 4]
    rt_2[0, k_start:k_end] = data[:, 5]
    varnce_rt_1[0, k_start:k_end] = data[:, 6]
    varnce_rt_2[0, k_start:k_end] = data[:, 7]
    thl_1[0, k_start:k_end] = data[:, 8]
    thl_2[0, k_start:k_end] = data[:, 9]
    varnce_thl_1[0, k_start:k_end] = data[:, 10]
    varnce_thl_2[0, k_start:k_end] = data[:, 11]
    corr_w_rt_1[0, k_start:k_end] = data[:, 12]
    corr_w_rt_2[0, k_start:k_end] = data[:, 13]
    corr_w_thl_1[0, k_start:k_end] = data[:, 14]
    corr_w_thl_2[0, k_start:k_end] = data[:, 15]
    corr_rt_thl_1[0, k_start:k_end] = data[:, 16]
    corr_rt_thl_2[0, k_start:k_end] = data[:, 17]
    alpha_thl[0, k_start:k_end] = data[:, 18]
    alpha_rt[0, k_start:k_end] = data[:, 19]
    crt_1[0, k_start:k_end] = data[:, 20]
    crt_2[0, k_start:k_end] = data[:, 21]
    cthl_1[0, k_start:k_end] = data[:, 22]
    cthl_2[0, k_start:k_end] = data[:, 23]
    chi_1[0, k_start:k_end] = data[:, 24]
    chi_2[0, k_start:k_end] = data[:, 25]
    stdev_chi_1[0, k_start:k_end] = data[:, 26]
    stdev_chi_2[0, k_start:k_end] = data[:, 27]
    stdev_eta_1[0, k_start:k_end] = data[:, 28]
    stdev_eta_2[0, k_start:k_end] = data[:, 29]
    covar_chi_eta_1[0, k_start:k_end] = data[:, 30]
    covar_chi_eta_2[0, k_start:k_end] = data[:, 31]
    corr_w_chi_1[0, k_start:k_end] = data[:, 32]
    corr_w_chi_2[0, k_start:k_end] = data[:, 33]
    corr_w_eta_1[0, k_start:k_end] = data[:, 34]
    corr_w_eta_2[0, k_start:k_end] = data[:, 35]
    corr_chi_eta_1[0, k_start:k_end] = data[:, 36]
    corr_chi_eta_2[0, k_start:k_end] = data[:, 37]
    rsatl_1[0, k_start:k_end] = data[:, 38]
    rsatl_2[0, k_start:k_end] = data[:, 39]
    rc_1[0, k_start:k_end] = data[:, 40]
    rc_2[0, k_start:k_end] = data[:, 41]
    cloud_frac_1[0, k_start:k_end] = data[:, 42]
    cloud_frac_2[0, k_start:k_end] = data[:, 43]
    mixt_frac[0, k_start:k_end] = data[:, 44]
    ice_supersat_frac_1[0, k_start:k_end] = data[:, 45]
    ice_supersat_frac_2[0, k_start:k_end] = data[:, 46]

    return pdf_parameter(
        ngrdcol=pdf_params.ngrdcol,
        nz=pdf_params.nz,
        w_1=w_1,
        w_2=w_2,
        varnce_w_1=varnce_w_1,
        varnce_w_2=varnce_w_2,
        rt_1=rt_1,
        rt_2=rt_2,
        varnce_rt_1=varnce_rt_1,
        varnce_rt_2=varnce_rt_2,
        thl_1=thl_1,
        thl_2=thl_2,
        varnce_thl_1=varnce_thl_1,
        varnce_thl_2=varnce_thl_2,
        corr_w_rt_1=corr_w_rt_1,
        corr_w_rt_2=corr_w_rt_2,
        corr_w_thl_1=corr_w_thl_1,
        corr_w_thl_2=corr_w_thl_2,
        corr_rt_thl_1=corr_rt_thl_1,
        corr_rt_thl_2=corr_rt_thl_2,
        alpha_thl=alpha_thl,
        alpha_rt=alpha_rt,
        crt_1=crt_1,
        crt_2=crt_2,
        cthl_1=cthl_1,
        cthl_2=cthl_2,
        chi_1=chi_1,
        chi_2=chi_2,
        stdev_chi_1=stdev_chi_1,
        stdev_chi_2=stdev_chi_2,
        stdev_eta_1=stdev_eta_1,
        stdev_eta_2=stdev_eta_2,
        covar_chi_eta_1=covar_chi_eta_1,
        covar_chi_eta_2=covar_chi_eta_2,
        corr_w_chi_1=corr_w_chi_1,
        corr_w_chi_2=corr_w_chi_2,
        corr_w_eta_1=corr_w_eta_1,
        corr_w_eta_2=corr_w_eta_2,
        corr_chi_eta_1=corr_chi_eta_1,
        corr_chi_eta_2=corr_chi_eta_2,
        rsatl_1=rsatl_1,
        rsatl_2=rsatl_2,
        rc_1=rc_1,
        rc_2=rc_2,
        cloud_frac_1=cloud_frac_1,
        cloud_frac_2=cloud_frac_2,
        mixt_frac=mixt_frac,
        ice_supersat_frac_1=ice_supersat_frac_1,
        ice_supersat_frac_2=ice_supersat_frac_2,
    )
