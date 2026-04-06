"""Fortran module-storage push/pull converters for pdf_parameter/implicit_coefs_terms."""

import numpy as np

import clubb_f2py

from clubb_python.derived_types.pdf_params import implicit_coefs_terms, pdf_parameter


def _as_2d(name: str, arr, ngrdcol: int, nz: int):
    out = np.asfortranarray(arr, dtype=np.float64)
    if out.shape != (ngrdcol, nz):
        raise ValueError(f"{name} shape mismatch: expected {(ngrdcol, nz)}, got {out.shape}")
    return out


def _as_3d(name: str, arr, ngrdcol: int, nz: int, sclr_dim: int):
    out = np.asfortranarray(arr, dtype=np.float64)
    if out.shape != (ngrdcol, nz, sclr_dim):
        raise ValueError(
            f"{name} shape mismatch: expected {(ngrdcol, nz, sclr_dim)}, got {out.shape}"
        )
    return out


def get_fortran_pdf_params() -> pdf_parameter:
    """Pull `stored_pdf_params` (zt) from Fortran module storage."""
    ngrdcol, nz, _, _ = clubb_f2py.get_pdf_params_dims()
    (
        w_1,
        w_2,
        varnce_w_1,
        varnce_w_2,
        rt_1,
        rt_2,
        varnce_rt_1,
        varnce_rt_2,
        thl_1,
        thl_2,
        varnce_thl_1,
        varnce_thl_2,
        corr_w_rt_1,
        corr_w_rt_2,
        corr_w_thl_1,
        corr_w_thl_2,
        corr_rt_thl_1,
        corr_rt_thl_2,
        alpha_thl,
        alpha_rt,
        crt_1,
        crt_2,
        cthl_1,
        cthl_2,
        chi_1,
        chi_2,
        stdev_chi_1,
        stdev_chi_2,
        stdev_eta_1,
        stdev_eta_2,
        covar_chi_eta_1,
        covar_chi_eta_2,
        corr_w_chi_1,
        corr_w_chi_2,
        corr_w_eta_1,
        corr_w_eta_2,
        corr_chi_eta_1,
        corr_chi_eta_2,
        rsatl_1,
        rsatl_2,
        rc_1,
        rc_2,
        cloud_frac_1,
        cloud_frac_2,
        mixt_frac,
        ice_supersat_frac_1,
        ice_supersat_frac_2,
    ) = clubb_f2py.get_pdf_params_fields(ngrdcol, nz, 1)

    return pdf_parameter(
        ngrdcol=int(ngrdcol),
        nz=int(nz),
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


def set_fortran_pdf_params(p: pdf_parameter):
    """Push `stored_pdf_params` (zt) into Fortran module storage."""
    ncol = int(p.ngrdcol)
    nz = int(p.nz)
    clubb_f2py.set_pdf_params_fields(
        1,
        _as_2d('w_1', p.w_1, ncol, nz),
        _as_2d('w_2', p.w_2, ncol, nz),
        _as_2d('varnce_w_1', p.varnce_w_1, ncol, nz),
        _as_2d('varnce_w_2', p.varnce_w_2, ncol, nz),
        _as_2d('rt_1', p.rt_1, ncol, nz),
        _as_2d('rt_2', p.rt_2, ncol, nz),
        _as_2d('varnce_rt_1', p.varnce_rt_1, ncol, nz),
        _as_2d('varnce_rt_2', p.varnce_rt_2, ncol, nz),
        _as_2d('thl_1', p.thl_1, ncol, nz),
        _as_2d('thl_2', p.thl_2, ncol, nz),
        _as_2d('varnce_thl_1', p.varnce_thl_1, ncol, nz),
        _as_2d('varnce_thl_2', p.varnce_thl_2, ncol, nz),
        _as_2d('corr_w_rt_1', p.corr_w_rt_1, ncol, nz),
        _as_2d('corr_w_rt_2', p.corr_w_rt_2, ncol, nz),
        _as_2d('corr_w_thl_1', p.corr_w_thl_1, ncol, nz),
        _as_2d('corr_w_thl_2', p.corr_w_thl_2, ncol, nz),
        _as_2d('corr_rt_thl_1', p.corr_rt_thl_1, ncol, nz),
        _as_2d('corr_rt_thl_2', p.corr_rt_thl_2, ncol, nz),
        _as_2d('alpha_thl', p.alpha_thl, ncol, nz),
        _as_2d('alpha_rt', p.alpha_rt, ncol, nz),
        _as_2d('crt_1', p.crt_1, ncol, nz),
        _as_2d('crt_2', p.crt_2, ncol, nz),
        _as_2d('cthl_1', p.cthl_1, ncol, nz),
        _as_2d('cthl_2', p.cthl_2, ncol, nz),
        _as_2d('chi_1', p.chi_1, ncol, nz),
        _as_2d('chi_2', p.chi_2, ncol, nz),
        _as_2d('stdev_chi_1', p.stdev_chi_1, ncol, nz),
        _as_2d('stdev_chi_2', p.stdev_chi_2, ncol, nz),
        _as_2d('stdev_eta_1', p.stdev_eta_1, ncol, nz),
        _as_2d('stdev_eta_2', p.stdev_eta_2, ncol, nz),
        _as_2d('covar_chi_eta_1', p.covar_chi_eta_1, ncol, nz),
        _as_2d('covar_chi_eta_2', p.covar_chi_eta_2, ncol, nz),
        _as_2d('corr_w_chi_1', p.corr_w_chi_1, ncol, nz),
        _as_2d('corr_w_chi_2', p.corr_w_chi_2, ncol, nz),
        _as_2d('corr_w_eta_1', p.corr_w_eta_1, ncol, nz),
        _as_2d('corr_w_eta_2', p.corr_w_eta_2, ncol, nz),
        _as_2d('corr_chi_eta_1', p.corr_chi_eta_1, ncol, nz),
        _as_2d('corr_chi_eta_2', p.corr_chi_eta_2, ncol, nz),
        _as_2d('rsatl_1', p.rsatl_1, ncol, nz),
        _as_2d('rsatl_2', p.rsatl_2, ncol, nz),
        _as_2d('rc_1', p.rc_1, ncol, nz),
        _as_2d('rc_2', p.rc_2, ncol, nz),
        _as_2d('cloud_frac_1', p.cloud_frac_1, ncol, nz),
        _as_2d('cloud_frac_2', p.cloud_frac_2, ncol, nz),
        _as_2d('mixt_frac', p.mixt_frac, ncol, nz),
        _as_2d('ice_supersat_frac_1', p.ice_supersat_frac_1, ncol, nz),
        _as_2d('ice_supersat_frac_2', p.ice_supersat_frac_2, ncol, nz),
    )


def get_fortran_pdf_params_zm() -> pdf_parameter:
    """Pull `stored_pdf_params_zm` from Fortran module storage."""
    _, _, ngrdcol, nz = clubb_f2py.get_pdf_params_dims()
    (
        w_1,
        w_2,
        varnce_w_1,
        varnce_w_2,
        rt_1,
        rt_2,
        varnce_rt_1,
        varnce_rt_2,
        thl_1,
        thl_2,
        varnce_thl_1,
        varnce_thl_2,
        corr_w_rt_1,
        corr_w_rt_2,
        corr_w_thl_1,
        corr_w_thl_2,
        corr_rt_thl_1,
        corr_rt_thl_2,
        alpha_thl,
        alpha_rt,
        crt_1,
        crt_2,
        cthl_1,
        cthl_2,
        chi_1,
        chi_2,
        stdev_chi_1,
        stdev_chi_2,
        stdev_eta_1,
        stdev_eta_2,
        covar_chi_eta_1,
        covar_chi_eta_2,
        corr_w_chi_1,
        corr_w_chi_2,
        corr_w_eta_1,
        corr_w_eta_2,
        corr_chi_eta_1,
        corr_chi_eta_2,
        rsatl_1,
        rsatl_2,
        rc_1,
        rc_2,
        cloud_frac_1,
        cloud_frac_2,
        mixt_frac,
        ice_supersat_frac_1,
        ice_supersat_frac_2,
    ) = clubb_f2py.get_pdf_params_fields(ngrdcol, nz, 2)

    return pdf_parameter(
        ngrdcol=int(ngrdcol),
        nz=int(nz),
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


def set_fortran_pdf_params_zm(p: pdf_parameter):
    """Push `stored_pdf_params_zm` into Fortran module storage."""
    ncol = int(p.ngrdcol)
    nz = int(p.nz)
    clubb_f2py.set_pdf_params_fields(
        2,
        _as_2d('w_1', p.w_1, ncol, nz),
        _as_2d('w_2', p.w_2, ncol, nz),
        _as_2d('varnce_w_1', p.varnce_w_1, ncol, nz),
        _as_2d('varnce_w_2', p.varnce_w_2, ncol, nz),
        _as_2d('rt_1', p.rt_1, ncol, nz),
        _as_2d('rt_2', p.rt_2, ncol, nz),
        _as_2d('varnce_rt_1', p.varnce_rt_1, ncol, nz),
        _as_2d('varnce_rt_2', p.varnce_rt_2, ncol, nz),
        _as_2d('thl_1', p.thl_1, ncol, nz),
        _as_2d('thl_2', p.thl_2, ncol, nz),
        _as_2d('varnce_thl_1', p.varnce_thl_1, ncol, nz),
        _as_2d('varnce_thl_2', p.varnce_thl_2, ncol, nz),
        _as_2d('corr_w_rt_1', p.corr_w_rt_1, ncol, nz),
        _as_2d('corr_w_rt_2', p.corr_w_rt_2, ncol, nz),
        _as_2d('corr_w_thl_1', p.corr_w_thl_1, ncol, nz),
        _as_2d('corr_w_thl_2', p.corr_w_thl_2, ncol, nz),
        _as_2d('corr_rt_thl_1', p.corr_rt_thl_1, ncol, nz),
        _as_2d('corr_rt_thl_2', p.corr_rt_thl_2, ncol, nz),
        _as_2d('alpha_thl', p.alpha_thl, ncol, nz),
        _as_2d('alpha_rt', p.alpha_rt, ncol, nz),
        _as_2d('crt_1', p.crt_1, ncol, nz),
        _as_2d('crt_2', p.crt_2, ncol, nz),
        _as_2d('cthl_1', p.cthl_1, ncol, nz),
        _as_2d('cthl_2', p.cthl_2, ncol, nz),
        _as_2d('chi_1', p.chi_1, ncol, nz),
        _as_2d('chi_2', p.chi_2, ncol, nz),
        _as_2d('stdev_chi_1', p.stdev_chi_1, ncol, nz),
        _as_2d('stdev_chi_2', p.stdev_chi_2, ncol, nz),
        _as_2d('stdev_eta_1', p.stdev_eta_1, ncol, nz),
        _as_2d('stdev_eta_2', p.stdev_eta_2, ncol, nz),
        _as_2d('covar_chi_eta_1', p.covar_chi_eta_1, ncol, nz),
        _as_2d('covar_chi_eta_2', p.covar_chi_eta_2, ncol, nz),
        _as_2d('corr_w_chi_1', p.corr_w_chi_1, ncol, nz),
        _as_2d('corr_w_chi_2', p.corr_w_chi_2, ncol, nz),
        _as_2d('corr_w_eta_1', p.corr_w_eta_1, ncol, nz),
        _as_2d('corr_w_eta_2', p.corr_w_eta_2, ncol, nz),
        _as_2d('corr_chi_eta_1', p.corr_chi_eta_1, ncol, nz),
        _as_2d('corr_chi_eta_2', p.corr_chi_eta_2, ncol, nz),
        _as_2d('rsatl_1', p.rsatl_1, ncol, nz),
        _as_2d('rsatl_2', p.rsatl_2, ncol, nz),
        _as_2d('rc_1', p.rc_1, ncol, nz),
        _as_2d('rc_2', p.rc_2, ncol, nz),
        _as_2d('cloud_frac_1', p.cloud_frac_1, ncol, nz),
        _as_2d('cloud_frac_2', p.cloud_frac_2, ncol, nz),
        _as_2d('mixt_frac', p.mixt_frac, ncol, nz),
        _as_2d('ice_supersat_frac_1', p.ice_supersat_frac_1, ncol, nz),
        _as_2d('ice_supersat_frac_2', p.ice_supersat_frac_2, ncol, nz),
    )


def get_fortran_implicit_coefs() -> implicit_coefs_terms:
    """Pull `stored_pdf_implicit_coefs_terms` from Fortran module storage."""
    ngrdcol, nz, sclr_dim = clubb_f2py.get_implicit_coefs_dims()
    sclr_dim_transport = max(int(sclr_dim), 1)
    (
        coef_wp4_implicit,
        coef_wp2rtp_implicit,
        term_wp2rtp_explicit,
        coef_wp2thlp_implicit,
        term_wp2thlp_explicit,
        coef_wp2up_implicit,
        term_wp2up_explicit,
        coef_wp2vp_implicit,
        term_wp2vp_explicit,
        coef_wprtp2_implicit,
        term_wprtp2_explicit,
        coef_wpthlp2_implicit,
        term_wpthlp2_explicit,
        coef_wprtpthlp_implicit,
        term_wprtpthlp_explicit,
        coef_wpup2_implicit,
        term_wpup2_explicit,
        coef_wpvp2_implicit,
        term_wpvp2_explicit,
        coef_wp2sclrp_implicit,
        term_wp2sclrp_explicit,
        coef_wpsclrp2_implicit,
        term_wpsclrp2_explicit,
        coef_wprtpsclrp_implicit,
        term_wprtpsclrp_explicit,
        coef_wpthlpsclrp_implicit,
        term_wpthlpsclrp_explicit,
    ) = clubb_f2py.get_implicit_coefs_fields(ngrdcol, nz, sclr_dim, sclr_dim_transport=sclr_dim_transport)

    if sclr_dim <= 0:
        coef_wp2sclrp_implicit = None
        term_wp2sclrp_explicit = None
        coef_wpsclrp2_implicit = None
        term_wpsclrp2_explicit = None
        coef_wprtpsclrp_implicit = None
        term_wprtpsclrp_explicit = None
        coef_wpthlpsclrp_implicit = None
        term_wpthlpsclrp_explicit = None

    return implicit_coefs_terms(
        ngrdcol=int(ngrdcol),
        nz=int(nz),
        sclr_dim=int(sclr_dim),
        coef_wp4_implicit=coef_wp4_implicit,
        coef_wp2rtp_implicit=coef_wp2rtp_implicit,
        term_wp2rtp_explicit=term_wp2rtp_explicit,
        coef_wp2thlp_implicit=coef_wp2thlp_implicit,
        term_wp2thlp_explicit=term_wp2thlp_explicit,
        coef_wp2up_implicit=coef_wp2up_implicit,
        term_wp2up_explicit=term_wp2up_explicit,
        coef_wp2vp_implicit=coef_wp2vp_implicit,
        term_wp2vp_explicit=term_wp2vp_explicit,
        coef_wprtp2_implicit=coef_wprtp2_implicit,
        term_wprtp2_explicit=term_wprtp2_explicit,
        coef_wpthlp2_implicit=coef_wpthlp2_implicit,
        term_wpthlp2_explicit=term_wpthlp2_explicit,
        coef_wprtpthlp_implicit=coef_wprtpthlp_implicit,
        term_wprtpthlp_explicit=term_wprtpthlp_explicit,
        coef_wpup2_implicit=coef_wpup2_implicit,
        term_wpup2_explicit=term_wpup2_explicit,
        coef_wpvp2_implicit=coef_wpvp2_implicit,
        term_wpvp2_explicit=term_wpvp2_explicit,
        coef_wp2sclrp_implicit=coef_wp2sclrp_implicit,
        term_wp2sclrp_explicit=term_wp2sclrp_explicit,
        coef_wpsclrp2_implicit=coef_wpsclrp2_implicit,
        term_wpsclrp2_explicit=term_wpsclrp2_explicit,
        coef_wprtpsclrp_implicit=coef_wprtpsclrp_implicit,
        term_wprtpsclrp_explicit=term_wprtpsclrp_explicit,
        coef_wpthlpsclrp_implicit=coef_wpthlpsclrp_implicit,
        term_wpthlpsclrp_explicit=term_wpthlpsclrp_explicit,
    )


def set_fortran_implicit_coefs(ic: implicit_coefs_terms):
    """Push `stored_pdf_implicit_coefs_terms` into Fortran module storage."""
    ncol = int(ic.ngrdcol)
    nz = int(ic.nz)
    sclr_dim = int(ic.sclr_dim)
    sclr_dim_transport = max(sclr_dim, 1)
    coef_wp2sclrp_implicit = np.zeros((ncol, nz, sclr_dim_transport), dtype=np.float64, order='F')
    term_wp2sclrp_explicit = np.zeros((ncol, nz, sclr_dim_transport), dtype=np.float64, order='F')
    coef_wpsclrp2_implicit = np.zeros((ncol, nz, sclr_dim_transport), dtype=np.float64, order='F')
    term_wpsclrp2_explicit = np.zeros((ncol, nz, sclr_dim_transport), dtype=np.float64, order='F')
    coef_wprtpsclrp_implicit = np.zeros((ncol, nz, sclr_dim_transport), dtype=np.float64, order='F')
    term_wprtpsclrp_explicit = np.zeros((ncol, nz, sclr_dim_transport), dtype=np.float64, order='F')
    coef_wpthlpsclrp_implicit = np.zeros((ncol, nz, sclr_dim_transport), dtype=np.float64, order='F')
    term_wpthlpsclrp_explicit = np.zeros((ncol, nz, sclr_dim_transport), dtype=np.float64, order='F')

    if sclr_dim > 0:
        if ic.coef_wp2sclrp_implicit is None or ic.term_wp2sclrp_explicit is None:
            raise ValueError('3D implicit fields are required when sclr_dim > 0')
        if ic.coef_wpsclrp2_implicit is None or ic.term_wpsclrp2_explicit is None:
            raise ValueError('3D implicit fields are required when sclr_dim > 0')
        if ic.coef_wprtpsclrp_implicit is None or ic.term_wprtpsclrp_explicit is None:
            raise ValueError('3D implicit fields are required when sclr_dim > 0')
        if ic.coef_wpthlpsclrp_implicit is None or ic.term_wpthlpsclrp_explicit is None:
            raise ValueError('3D implicit fields are required when sclr_dim > 0')

        coef_wp2sclrp_implicit = _as_3d('coef_wp2sclrp_implicit', ic.coef_wp2sclrp_implicit, ncol, nz, sclr_dim)
        term_wp2sclrp_explicit = _as_3d('term_wp2sclrp_explicit', ic.term_wp2sclrp_explicit, ncol, nz, sclr_dim)
        coef_wpsclrp2_implicit = _as_3d('coef_wpsclrp2_implicit', ic.coef_wpsclrp2_implicit, ncol, nz, sclr_dim)
        term_wpsclrp2_explicit = _as_3d('term_wpsclrp2_explicit', ic.term_wpsclrp2_explicit, ncol, nz, sclr_dim)
        coef_wprtpsclrp_implicit = _as_3d('coef_wprtpsclrp_implicit', ic.coef_wprtpsclrp_implicit, ncol, nz, sclr_dim)
        term_wprtpsclrp_explicit = _as_3d('term_wprtpsclrp_explicit', ic.term_wprtpsclrp_explicit, ncol, nz, sclr_dim)
        coef_wpthlpsclrp_implicit = _as_3d('coef_wpthlpsclrp_implicit', ic.coef_wpthlpsclrp_implicit, ncol, nz, sclr_dim)
        term_wpthlpsclrp_explicit = _as_3d('term_wpthlpsclrp_explicit', ic.term_wpthlpsclrp_explicit, ncol, nz, sclr_dim)

    clubb_f2py.set_implicit_coefs_fields(
        sclr_dim,
        _as_2d('coef_wp4_implicit', ic.coef_wp4_implicit, ncol, nz),
        _as_2d('coef_wp2rtp_implicit', ic.coef_wp2rtp_implicit, ncol, nz),
        _as_2d('term_wp2rtp_explicit', ic.term_wp2rtp_explicit, ncol, nz),
        _as_2d('coef_wp2thlp_implicit', ic.coef_wp2thlp_implicit, ncol, nz),
        _as_2d('term_wp2thlp_explicit', ic.term_wp2thlp_explicit, ncol, nz),
        _as_2d('coef_wp2up_implicit', ic.coef_wp2up_implicit, ncol, nz),
        _as_2d('term_wp2up_explicit', ic.term_wp2up_explicit, ncol, nz),
        _as_2d('coef_wp2vp_implicit', ic.coef_wp2vp_implicit, ncol, nz),
        _as_2d('term_wp2vp_explicit', ic.term_wp2vp_explicit, ncol, nz),
        _as_2d('coef_wprtp2_implicit', ic.coef_wprtp2_implicit, ncol, nz),
        _as_2d('term_wprtp2_explicit', ic.term_wprtp2_explicit, ncol, nz),
        _as_2d('coef_wpthlp2_implicit', ic.coef_wpthlp2_implicit, ncol, nz),
        _as_2d('term_wpthlp2_explicit', ic.term_wpthlp2_explicit, ncol, nz),
        _as_2d('coef_wprtpthlp_implicit', ic.coef_wprtpthlp_implicit, ncol, nz),
        _as_2d('term_wprtpthlp_explicit', ic.term_wprtpthlp_explicit, ncol, nz),
        _as_2d('coef_wpup2_implicit', ic.coef_wpup2_implicit, ncol, nz),
        _as_2d('term_wpup2_explicit', ic.term_wpup2_explicit, ncol, nz),
        _as_2d('coef_wpvp2_implicit', ic.coef_wpvp2_implicit, ncol, nz),
        _as_2d('term_wpvp2_explicit', ic.term_wpvp2_explicit, ncol, nz),
        coef_wp2sclrp_implicit,
        term_wp2sclrp_explicit,
        coef_wpsclrp2_implicit,
        term_wpsclrp2_explicit,
        coef_wprtpsclrp_implicit,
        term_wprtpsclrp_explicit,
        coef_wpthlpsclrp_implicit,
        term_wpthlpsclrp_explicit,
    )
