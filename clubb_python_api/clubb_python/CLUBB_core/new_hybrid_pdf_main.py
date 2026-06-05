"""User-facing wrappers for routines from CLUBB_core/new_hybrid_pdf_main.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.pdf_params import implicit_coefs_terms
from clubb_python.derived_types.pdf_params_converter import set_fortran_implicit_coefs


def new_hybrid_pdf_driver(
    nz: int | np.ndarray, ngrdcol: int | np.ndarray, sclr_dim: int | np.ndarray,
    wm=None, rtm=None, thlm=None, um=None, vm=None, wp2=None, rtp2=None, thlp2=None, up2=None, vp2=None,
    skw=None, wprtp=None, wpthlp=None, upwp=None, vpwp=None, sclrm=None, sclrp2=None, wpsclrp=None, clubb_params=None,
    slope_coef_spread_dg_means_w=None, pdf_component_stdev_factor_w=None, skrt=None, skthl=None, sku=None, skv=None,
    sksclr=None, mixt_frac=None, *,
    pdf_implicit_coefs_terms: implicit_coefs_terms, **compat_kwargs,
):
    """Compute full hybrid-PDF moments/mixture fraction for w/rt/thl/u/v/(sclr)."""
    if not np.isscalar(nz):
        gamma_skw_fnc_legacy = sclrm
        (
            wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2, skw, wprtp, wpthlp,
            upwp, vpwp, sclrm, sclrp2, wpsclrp, slope_coef_spread_dg_means_w,
            pdf_component_stdev_factor_w, skrt, skthl, sku, skv, sksclr,
        ) = (
            nz, ngrdcol, sclr_dim, wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2,
            skw, wprtp, wpthlp, upwp, vpwp, sclrp2, wpsclrp, clubb_params,
            slope_coef_spread_dg_means_w, pdf_component_stdev_factor_w, skrt, skthl,
        )
        compat_kwargs.setdefault("gamma_skw_fnc", gamma_skw_fnc_legacy)
        ngrdcol = int(np.asarray(wm).shape[0])
        nz = int(np.asarray(wm).shape[1])
        sclr_dim = int(np.asarray(sclrm).shape[2]) if np.asarray(sclrm).ndim == 3 else 0
    if mixt_frac is None:
        mixt_frac = np.zeros_like(wm)
    set_fortran_implicit_coefs(pdf_implicit_coefs_terms)
    gamma_skw_fnc = compat_kwargs.pop("gamma_skw_fnc", np.zeros_like(wm))
    result = clubb_f2py.f2py_new_hybrid_pdf_driver(
        sclr_dim, f_arr(wm), f_arr(rtm), f_arr(thlm), f_arr(um), f_arr(vm), f_arr(wp2), f_arr(rtp2), f_arr(thlp2),
        f_arr(up2), f_arr(vp2), f_arr(skw), f_arr(wprtp), f_arr(wpthlp), f_arr(upwp), f_arr(vpwp),
        f_arr(sclrm), f_arr(sclrp2), f_arr(wpsclrp), f_arr(gamma_skw_fnc), f_arr(slope_coef_spread_dg_means_w),
        f_arr(pdf_component_stdev_factor_w), f_arr(skrt), f_arr(skthl), f_arr(sku), f_arr(skv), f_arr(sksclr),
    )
    return result
