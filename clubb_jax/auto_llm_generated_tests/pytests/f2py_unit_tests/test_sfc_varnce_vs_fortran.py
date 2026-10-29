"""Check the sfc_varnce_module.F90 port (calc_sfc_varnce) against Fortran."""

import numpy as np

import clubb_python.CLUBB_core.sfc_varnce_module as fortran

from clubb_jax.auto_llm_generated_tests.pytests.f2py_unit_test_helpers import max_abs_diff, no_stats, setup_matching_grids
from clubb_jax.src.CLUBB_core.sfc_varnce_module import calc_sfc_varnce

_NG, _DZ, _ZTOP = 2, 40.0, 1200.0


_NAMES = ("wp2", "up2", "vp2", "thlp2", "rtp2", "rtpthlp", "sclrp2", "sclrprtp", "sclrpthlp")


def test_f2py_oracle():
    """calc_sfc_varnce surface second moments (incl. one scalar tracking rt) match Fortran."""
    from clubb_python.derived_types.sclr_idx import SclrIdx as FortranSclrIdx
    from clubb_python.derived_types.err_info import ErrInfo as FortranErrInfo
    from clubb_jax.src.CLUBB_core.err_info import ErrInfo
    from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx
    jgr, fgr = setup_matching_grids(_NG, _DZ, _ZTOP)
    ng, nzm = jgr.zm.shape
    nzt = nzm - 1
    sd = 1
    idx = dict(iisclr_rt=1, iisclr_thl=-1, iisclr_CO2=-1, iiedsclr_rt=-1, iiedsclr_thl=-1, iiedsclr_CO2=-1)
    sfce = np.zeros(ng)
    rng = np.random.default_rng(2)
    worst = 0.0
    for it in range(20):
        upwp = rng.uniform(-0.1, 0.1, ng); vpwp = rng.uniform(-0.1, 0.1, ng)
        wpthlp = rng.uniform(-0.1, 0.1, (ng, nzm)); wprtp = rng.uniform(-1e-4, 1e-4, ng)
        um = rng.uniform(-5, 5, (ng, nzt)); vm = rng.uniform(-5, 5, (ng, nzt)); lsu = rng.uniform(10, 500, (ng, nzt))
        wpsclrp = rng.uniform(-1e-4, 1e-4, (ng, sd))
        splat = rng.uniform(0.0, 1e-3, (ng, nzm)); tau = rng.uniform(100.0, 1000.0, (ng, nzm))
        T0 = 300.0; up2c = np.full(ng, 2.0); ac = np.full(ng, 1.8)
        wp2 = rng.uniform(1e-2, 1, (ng, nzm)); up2 = rng.uniform(1e-2, 1, (ng, nzm)); vp2 = rng.uniform(1e-2, 1, (ng, nzm))
        thlp2 = rng.uniform(1e-3, 1, (ng, nzm)); rtp2 = rng.uniform(1e-7, 1e-5, (ng, nzm))
        rtpthlp = rng.uniform(-1e-3, 1e-3, (ng, nzm))
        sclr = [rng.uniform(1e-7, 1e-5, (ng, nzm, sd)) for _ in range(3)]
        l_vary = bool(it % 2)
        common = (sfce, upwp, vpwp, wpthlp, wprtp, um, vm, lsu, wpsclrp, splat, tau, l_vary, T0, up2c, ac)
        moments = (wp2, up2, vp2, thlp2, rtp2, rtpthlp, *sclr)
        fr = fortran.calc_sfc_varnce(nzm, nzt, ng, sd, FortranSclrIdx(**idx), fgr, 60.0, *common,
                                     *(m.copy() for m in moments), FortranErrInfo(ngrdcol=ng))
        jr = calc_sfc_varnce(nzm, nzt, ng, sd, SclrIdx(**idx), jgr, 60.0, *common, no_stats(ng),
                             *moments, ErrInfo(ngrdcol=ng))
        assert not fr[-1].is_fatal()
        worst = max(worst, max_abs_diff(f"calc_sfc_varnce(l_vary_convect_depth={l_vary})",
                                        tuple(jr[:len(_NAMES)]), tuple(fr[:len(_NAMES)])))
    assert worst < 1e-12, f"calc_sfc_varnce Fortran mismatch {worst:.2e}"
