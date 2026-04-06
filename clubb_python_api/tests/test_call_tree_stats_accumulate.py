"""Test wrapper for stats_clubb_utilities::stats_accumulate."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def test_stats_accumulate_smoke():
    """stats_accumulate wrapper should run without setting err_code."""
    ngrdcol = 1
    nzm = 4
    nzt = 3
    sclr_dim = 1
    edsclr_dim = 1

    clubb_api.init_err_info(ngrdcol)
    clubb_api.reset_err_code()
    pdf_params = clubb_api.init_pdf_params(nzt, ngrdcol)
    pdf_params_zm = clubb_api.init_pdf_params_zm(nzm, ngrdcol)
    err_info = ErrInfo(ngrdcol=ngrdcol)
    registry = str(Path(__file__).resolve().parent / "test_stats_registry.in")
    clubb_params = clubb_api.init_clubb_params(ngrdcol, iunit=10, filename="")

    zt_levels = np.linspace(50.0, 250.0, nzt, dtype=np.float64)
    zm_levels = np.linspace(0.0, 300.0, nzm, dtype=np.float64)
    clubb_api.init_stats(
        registry_path=registry,
        output_path="",
        ncol=ngrdcol,
        stats_tsamp=60.0,
        stats_tout=120.0,
        dt_main=60.0,
        day_in=1,
        month_in=1,
        year_in=2000,
        time_initial=0.0,
        nzt=nzt,
        zt=zt_levels,
        nzm=nzm,
        zm=zm_levels,
        sclr_dim=sclr_dim,
        edsclr_dim=edsclr_dim,
        err_info=err_info,
        clubb_params=clubb_params,
        param_names=clubb_api.get_param_names(),
    )

    def zt(val):
        return np.full((ngrdcol, nzt), val, dtype=np.float64, order="F")

    def zm(val):
        return np.full((ngrdcol, nzm), val, dtype=np.float64, order="F")

    def zts(val):
        return np.full((ngrdcol, nzt, sclr_dim), val, dtype=np.float64, order="F")

    def zms(val):
        return np.full((ngrdcol, nzm, sclr_dim), val, dtype=np.float64, order="F")

    def zte(val):
        return np.full((ngrdcol, nzt, edsclr_dim), val, dtype=np.float64, order="F")

    def zme(val):
        return np.full((ngrdcol, nzm, edsclr_dim), val, dtype=np.float64, order="F")

    try:
        clubb_api.stats_accumulate(
            nzm=nzm, nzt=nzt, ngrdcol=ngrdcol, sclr_dim=sclr_dim, edsclr_dim=edsclr_dim,
            invrs_dzm=zm(1.0), zt=zt(100.0), dzm=zm(100.0), dzt=zt(100.0), dt=60.0,
            um=zt(1.0), vm=zt(0.0), upwp=zm(0.0), vpwp=zm(0.0), up2=zm(0.01), vp2=zm(0.01),
            thlm=zt(300.0), rtm=zt(0.01), wprtp=zm(0.0), wpthlp=zm(0.0),
            wp2=zm(0.01), wp3=zt(0.0), rtp2=zm(1.0e-6), rtp3=zt(0.0), thlp2=zm(0.01), thlp3=zt(0.0),
            rtpthlp=zm(0.0),
            wpthvp=zm(0.0), wp2thvp=zt(0.0), wp2up=zt(0.0), rtpthvp=zm(0.0), thlpthvp=zm(0.0),
            p_in_pa=zt(90000.0), exner=zt(1.0), rho=zt(1.0), rho_zm=zm(1.0),
            rho_ds_zm=zm(1.0), rho_ds_zt=zt(1.0), thv_ds_zm=zm(300.0), thv_ds_zt=zt(300.0),
            wm_zt=zt(0.0), wm_zm=zm(0.0), rcm=zt(0.0), wprcp=zm(0.0), rc_coef=zt(0.0),
            rc_coef_zm=zm(0.0),
            rcm_zm=zm(0.0), rtm_zm=zm(0.01), thlm_zm=zm(300.0), cloud_frac=zt(0.0),
            ice_supersat_frac=zt(0.0),
            cloud_frac_zm=zm(0.0), ice_supersat_frac_zm=zm(0.0), rcm_in_layer=zt(0.0),
            cloud_cover=zt(0.0), rcm_supersat_adj=zt(0.0), sigma_sqd_w=zm(0.01),
            thvm=zt(300.0), ug=zt(1.0), vg=zt(0.0), lscale=zt(100.0), wpthlp2=zt(0.0), wp2thlp=zt(0.0),
            wprtp2=zt(0.0), wp2rtp=zt(0.0),
            lscale_up=zt(100.0), lscale_down=zt(100.0), tau_zt=zt(60.0), kh_zt=zt(1.0), wp2rcp=zt(0.0),
            wprtpthlp=zt(0.0), sigma_sqd_w_zt=zt(0.01), rsat=zt(0.01), wp2_zt=zt(0.01),
            thlp2_zt=zt(0.01),
            wpthlp_zt=zt(0.0), wprtp_zt=zt(0.0), rtp2_zt=zt(1.0e-6), rtpthlp_zt=zt(0.0),
            up2_zt=zt(0.01),
            vp2_zt=zt(0.01), upwp_zt=zt(0.0), vpwp_zt=zt(0.0), wpup2=zt(0.0), wpvp2=zt(0.0),
            wp2up2=zm(0.0), wp2vp2=zm(0.0), wp4=zm(0.0),
            tau_zm=zm(60.0), kh_zm=zm(1.0), thlprcp=zm(0.0),
            rtprcp=zm(0.0), rcp2=zm(0.0), em=zm(0.0), a3_coef=zm(0.0), a3_coef_zt=zt(0.0),
            wp3_zm=zm(0.0), wp3_on_wp2=zm(0.0), wp3_on_wp2_zt=zt(0.0), skw_velocity=zm(0.0),
            w_up_in_cloud=zt(0.0), w_down_in_cloud=zt(0.0),
            cloudy_updraft_frac=zt(0.0), cloudy_downdraft_frac=zt(0.0),
            pdf_params=pdf_params,
            pdf_params_zm=pdf_params_zm,
            sclrm=zts(0.0), sclrp2=zms(0.0),
            sclrprtp=zms(0.0), sclrpthlp=zms(0.0), sclrm_forcing=zts(0.0), sclrpthvp=zms(0.0),
            wpsclrp=zms(0.0), sclrprcp=zms(0.0), wp2sclrp=zts(0.0), wpsclrp2=zts(0.0),
            wpsclrprtp=zts(0.0),
            wpsclrpthlp=zts(0.0), wpedsclrp=zme(0.0), edsclrm=zte(0.0),
            edsclrm_forcing=zte(0.0),
            saturation_formula=1,
            l_call_pdf_closure_twice=False,
        )
    finally:
        clubb_api.finalize_stats(err_info=err_info)

    assert int(clubb_api.get_err_code(ngrdcol)[0]) == 0
