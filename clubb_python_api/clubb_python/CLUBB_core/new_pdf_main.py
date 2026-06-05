"""User-facing wrappers for routines from CLUBB_core/new_pdf_main.F90."""

import numpy as np
from numpy import asfortranarray as f_arr

import clubb_f2py

from clubb_python.derived_types.pdf_params import implicit_coefs_terms
from clubb_python.derived_types.pdf_params_converter import set_fortran_implicit_coefs


def new_pdf_driver(
    nz: int | np.ndarray, ngrdcol: int | np.ndarray,
    wm=None, rtm=None, thlm=None, wp2=None, rtp2=None, thlp2=None, skw=None, wprtp=None, wpthlp=None,
    rtpthlp=None, clubb_params=None, skrt=None, skthl=None, *,
    pdf_implicit_coefs_terms: implicit_coefs_terms,
):
    """Compute full new-PDF moments/mixture fraction for w/rt/thl."""
    if not np.isscalar(nz):
        wm, rtm, thlm, wp2, rtp2, thlp2, skw, wprtp, wpthlp, rtpthlp, clubb_params, skrt, skthl = (
            nz, ngrdcol, wm, rtm, thlm, wp2, rtp2, thlp2, skw, wprtp, wpthlp, rtpthlp, clubb_params
        )
        ngrdcol = int(np.asarray(wm).shape[0])
        nz = int(np.asarray(wm).shape[1])
    set_fortran_implicit_coefs(pdf_implicit_coefs_terms)
    result = clubb_f2py.f2py_new_pdf_driver(
        f_arr(wm), f_arr(rtm), f_arr(thlm), f_arr(wp2), f_arr(rtp2), f_arr(thlp2), f_arr(skw),
        f_arr(wprtp), f_arr(wpthlp), f_arr(rtpthlp), f_arr(clubb_params), f_arr(skrt), f_arr(skthl),
    )
    return result
