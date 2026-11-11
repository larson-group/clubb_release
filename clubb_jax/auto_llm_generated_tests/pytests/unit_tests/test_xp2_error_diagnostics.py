"""The unsupported multiple-LHS path reports its failure without stored messages.

Synthetic inputs exercise the JAX guard at the single/multiple-LHS dispatch
in advance_xp2_xpyp_module.F90; the source supports both solve branches.
"""

import jax
import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core import ErrInfo, NuVertResDep, SclrIdx
from clubb_jax.src.CLUBB_core.advance_xp2_xpyp_module import advance_xp2_xpyp
from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.model_flags import iiPDF_ADG1
from clubb_jax.src.CLUBB_core.parameter_indices import iC2thl, nparams
from clubb_jax.src.CLUBB_core.pdf_params import init_pdf_implicit_coefs_terms_api


def test_multiple_lhs_diagnostic_uses_runtime_parameters(capsys):
    """Cached calls print on failure, retain incoming fields, and flag all columns."""
    gr = setup_grid(ngrdcol=2, deltaz=100., zm_init=0., zm_top=500., grid_type=1)
    nc, nt, nm = gr.ngrdcol, gr.nzt, gr.nzm
    zt = lambda v: jnp.full((nc, nt), v, dtype=jnp.float64)
    zm = lambda v: jnp.full((nc, nm), v, dtype=jnp.float64)
    scalar = jnp.zeros((nc, nm, 0))
    args = dict(
        nzm=nm, nzt=nt, ngrdcol=nc, sclr_dim=0, sclr_tol=jnp.zeros(0), gr=gr,
        sclr_idx=SclrIdx(-1, -1, -1, -1, -1, -1),
        invrs_tau_xp2_zm=zm(.01), invrs_tau_C4_zm=zm(.01), invrs_tau_C14_zm=zm(.01),
        wm_zm=zm(0), rtm=zt(.01), wprtp=zm(0), thlm=zt(300),
        wpthlp=zm(0), wpthvp=zm(0), um=zt(1), vm=zt(0),
        wp2=zm(.1), wp3=zt(0), upwp=zm(0), vpwp=zm(0), sigma_sqd_w=zm(.1),
        wprtp2=zt(0), wpthlp2=zt(0), wprtpthlp=zt(0), Kh_zt=zt(0),
        rtp2_forcing=zm(0), thlp2_forcing=zm(0), rtpthlp_forcing=zm(0),
        rho_ds_zm=zm(1), rho_ds_zt=zt(1), invrs_rho_ds_zm=zm(1),
        thv_ds_zm=zm(300), cloud_frac=zt(0),
        pdf_implicit_coefs_terms=init_pdf_implicit_coefs_terms_api(nt, nc, 0),
        dt=60., fcor_y=jnp.zeros(nc), sclrm=jnp.zeros((nc, nt, 0)),
        wpsclrp=scalar, wpsclrp2=jnp.zeros((nc, nt, 0)),
        wpsclrprtp=jnp.zeros((nc, nt, 0)), wpsclrpthlp=jnp.zeros((nc, nt, 0)),
        lhs_splat_wp2=zm(0),
        nu_vert_res_dep=NuVertResDep(nm, *(jnp.zeros(nc) for _ in range(7))),
        iiPDF_type=iiPDF_ADG1, tridiag_solve_method=2, fill_holes_type=0,
        l_ho_nontrad_coriolis=False, l_min_xp2_from_corr_wx=False,
        l_C2_cloud_frac=False, l_upwind_xpyp_ta=True, l_godunov_upwind_xpyp_ta=False,
        l_lmm_stepping=False, l_implemented=True,
        stats=JaxStats.empty(l_sample=False, names=(), ncol=nc, max_nlev=nm),
        rtp2=zm(.001), thlp2=zm(.01), rtpthlp=zm(0), up2=zm(.1), vp2=zm(.1),
        sclrp2=scalar, sclrprtp=scalar, sclrpthlp=scalar, err_info=ErrInfo.initialized(nc),
    )
    params = jnp.ones((nc, nparams))
    for failed in (False, True, True, False):
        args["clubb_params"] = params.at[1, iC2thl].set(2.) if failed else params
        result = advance_xp2_xpyp(**args)
        jax.block_until_ready(result)
        jax.effects_barrier()
        text = capsys.readouterr().err
        np.testing.assert_array_equal(result[-2].fatal_mask(), [failed, failed])
        if failed:
            assert "requires the multiple-LHS solve" in text
            assert "(l_explicit_turbulent_adv_xpyp or iiPDF_type == iiPDF_ADG1)" in text
            assert "in every column, C2thl and C2rtthl must each match C2rt" in text
            assert "within the relative epsilon tolerance" in text
            for actual, name in zip(result[:8], (
                "rtp2", "thlp2", "rtpthlp", "up2", "vp2", "sclrp2", "sclrprtp", "sclrpthlp",
            )):
                np.testing.assert_array_equal(actual, args[name])
        else:
            assert text == ""
