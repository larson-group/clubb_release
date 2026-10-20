"""Check the new-hybrid PDF driver port (new_hybrid_pdf_main.F90)
against Fortran."""

import os

import numpy as np

from utilities.output_paths import REPO_ROOT

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.new_hybrid_pdf_main import new_hybrid_pdf_driver

NG, NZ, SD = 2, 6, 1


def _inputs(rng):
    sh = (NG, NZ)
    return dict(
        wm=rng.uniform(-1, 1, sh), rtm=rng.uniform(0.005, 0.02, sh), thlm=rng.uniform(285, 305, sh),
        um=rng.uniform(-5, 5, sh), vm=rng.uniform(-5, 5, sh),
        wp2=rng.uniform(0.1, 2, sh), rtp2=rng.uniform(1e-7, 1e-6, sh), thlp2=rng.uniform(0.05, 1, sh),
        up2=rng.uniform(0.1, 2, sh), vp2=rng.uniform(0.1, 2, sh),
        Skw=rng.uniform(-2, 2, sh),
        wprtp=rng.uniform(-2e-4, 2e-4, sh), wpthlp=rng.uniform(-0.2, 0.2, sh),
        upwp=rng.uniform(-0.3, 0.3, sh), vpwp=rng.uniform(-0.3, 0.3, sh),
        gamma_Skw_fnc=rng.uniform(0.1, 0.9, sh),
        slope_coef_spread_DG_means_w=rng.uniform(1.0, 3.0, NG),
        pdf_component_stdev_factor_w=rng.uniform(0.8, 1.5, NG),
        Skrt=rng.uniform(-3, 3, sh), Skthl=rng.uniform(-3, 3, sh),
        Sku=rng.uniform(-3, 3, sh), Skv=rng.uniform(-3, 3, sh))


def _sclr(rng):
    return dict(sclrm=rng.uniform(0, 1, (NG, NZ, SD)), sclrp2=rng.uniform(1e-3, 1e-2, (NG, NZ, SD)),
                wpsclrp=rng.uniform(-1e-2, 1e-2, (NG, NZ, SD)), Sksclr=rng.uniform(-3, 3, (NG, NZ, SD)))
_TUNABLE_PARAMS = os.path.join(REPO_ROOT, "input", "parameter_and_flag_configs", "default", "tunable_parameters.in")
_ORDER = ['Skrt', 'Skthl', 'Sku', 'Skv', 'Sksclr', 'mu_w_1', 'mu_w_2', 'mu_rt_1', 'mu_rt_2',
          'mu_thl_1', 'mu_thl_2', 'mu_u_1', 'mu_u_2', 'mu_v_1', 'mu_v_2', 'sigma_w_1_sqd', 'sigma_w_2_sqd',
          'sigma_rt_1_sqd', 'sigma_rt_2_sqd', 'sigma_thl_1_sqd', 'sigma_thl_2_sqd', 'sigma_u_1_sqd',
          'sigma_u_2_sqd', 'sigma_v_1_sqd', 'sigma_v_2_sqd', 'mu_sclr_1', 'mu_sclr_2', 'sigma_sclr_1_sqd',
          'sigma_sclr_2_sqd', 'mixt_frac', 'sigma_sqd_w']


def test_f2py_oracle():
    """new_hybrid_pdf_driver PDF parameters and implicit/explicit turbulent-advection terms match Fortran.

    Fortran derives gamma_Skw_fnc from clubb_params (new_hybrid_pdf_main.F90 calls compute_gamma_Skw); the
    JAX port takes it as an argument, so it is computed from the same default parameters here.
    """
    from clubb_python.CLUBB_core.new_hybrid_pdf_main import new_hybrid_pdf_driver as fortran_driver
    from clubb_python.CLUBB_core.pdf_parameter_module import init_pdf_implicit
    from clubb_python.derived_types.pdf_params_converter import get_fortran_implicit_coefs
    from clubb_jax.src.CLUBB_core.pdf_params import init_pdf_implicit_coefs_terms_api
    from clubb_jax.src.CLUBB_core.Skx_module import compute_gamma_Skw
    from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params
    rng = np.random.default_rng(17)
    cp = np.asarray(init_clubb_params(NG, filename=_TUNABLE_PARAMS))
    worst = {}
    for _ in range(3):
        a = _inputs(rng); s = _sclr(rng)
        f = fortran_driver(
            NZ, NG, SD, a['wm'], a['rtm'], a['thlm'], a['um'], a['vm'], a['wp2'], a['rtp2'], a['thlp2'],
            a['up2'], a['vp2'], a['Skw'], a['wprtp'], a['wpthlp'], a['upwp'], a['vpwp'],
            s['sclrm'], s['sclrp2'], s['wpsclrp'], cp,
            a['slope_coef_spread_DG_means_w'], a['pdf_component_stdev_factor_w'],
            a['Skrt'], a['Skthl'], a['Sku'], a['Skv'], s['Sksclr'],
            pdf_implicit_coefs_terms=init_pdf_implicit(NZ, NG, SD))
        f_coefs = get_fortran_implicit_coefs()
        g = new_hybrid_pdf_driver(
            a['wm'], a['rtm'], a['thlm'], a['um'], a['vm'], a['wp2'], a['rtp2'], a['thlp2'],
            a['up2'], a['vp2'], a['Skw'], a['wprtp'], a['wpthlp'], a['upwp'], a['vpwp'],
            compute_gamma_Skw(NZ, NG, a['Skw'], cp, True),
            a['slope_coef_spread_DG_means_w'], a['pdf_component_stdev_factor_w'],
            a['Skrt'], a['Skthl'], a['Sku'], a['Skv'], init_pdf_implicit_coefs_terms_api(NZ, NG, SD),
            s['sclrm'], s['sclrp2'], s['wpsclrp'], s['Sksclr'])
        pairs = list(zip(_ORDER, f, (g[k] for k in _ORDER)))
        pairs += [(n, getattr(f_coefs, n), getattr(g["pdf_implicit_coefs_terms"], n))
                  for n in f_coefs._fields if getattr(f_coefs, n) is not None and np.ndim(getattr(f_coefs, n)) >= 2]
        for name, fi, gi in pairs:
            worst[name] = max(worst.get(name, 0.0), max_abs_diff(name, gi, fi))
    bad = {k: v for k, v in worst.items() if not v < 1e-11}
    assert not bad, f"new_hybrid_pdf_driver Fortran mismatch: {bad}"
    assert len(worst) > len(_ORDER), "no implicit coefficient fields were compared"
