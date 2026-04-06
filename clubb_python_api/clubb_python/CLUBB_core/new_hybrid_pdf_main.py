"""User-facing wrappers for routines from CLUBB_core/new_hybrid_pdf_main.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.pdf_params import implicit_coefs_terms
from clubb_python.derived_types.pdf_params_converter import set_fortran_implicit_coefs


def new_hybrid_pdf_driver(
    wm, rtm, thlm, um, vm, wp2, rtp2, thlp2, up2, vp2, skw, wprtp, wpthlp, upwp, vpwp,
    sclrm, sclrp2, wpsclrp, gamma_skw_fnc, slope_coef_spread_dg_means_w,
    pdf_component_stdev_factor_w, skrt, skthl, sku, skv, sksclr,
    pdf_implicit_coefs_terms: implicit_coefs_terms,
):
    """Compute full hybrid-PDF moments/mixture fraction for w/rt/thl/u/v/(sclr)."""
    set_fortran_implicit_coefs(pdf_implicit_coefs_terms)
    sclr_dim = int(sclrm.shape[2]) if np.ndim(sclrm) == 3 else 0
    result = clubb_f2py.f2py_new_hybrid_pdf_driver(
        sclr_dim, f_arr(wm), f_arr(rtm), f_arr(thlm), f_arr(um), f_arr(vm), f_arr(wp2), f_arr(rtp2), f_arr(thlp2),
        f_arr(up2), f_arr(vp2), f_arr(skw), f_arr(wprtp), f_arr(wpthlp), f_arr(upwp), f_arr(vpwp),
        f_arr(sclrm), f_arr(sclrp2), f_arr(wpsclrp), f_arr(gamma_skw_fnc), f_arr(slope_coef_spread_dg_means_w),
        f_arr(pdf_component_stdev_factor_w), f_arr(skrt), f_arr(skthl), f_arr(sku), f_arr(skv), f_arr(sksclr),
    )
    return result
