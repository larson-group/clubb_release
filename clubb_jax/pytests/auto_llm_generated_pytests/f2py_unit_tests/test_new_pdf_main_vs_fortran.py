"""Check the JAX new_pdf_driver port (new_pdf_main.F90) against Fortran."""

import numpy as np

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff
from clubb_jax.src.CLUBB_core.new_pdf_main import new_pdf_driver


_DRIVER_OUTPUTS = ("Skrt", "Skthl", "mu_w_1", "mu_w_2", "mu_rt_1", "mu_rt_2", "mu_thl_1", "mu_thl_2",
                   "sigma_w_1_sqd", "sigma_w_2_sqd", "sigma_rt_1_sqd", "sigma_rt_2_sqd",
                   "sigma_thl_1_sqd", "sigma_thl_2_sqd", "mixt_frac")


def test_driver_f2py():
    """new_pdf_driver PDF parameters and implicit/explicit turbulent-advection terms match Fortran."""
    from clubb_python.CLUBB_core.new_pdf_main import new_pdf_driver as fortran_new_pdf_driver
    from clubb_python.CLUBB_core.pdf_parameter_module import init_pdf_implicit
    from clubb_python.derived_types.pdf_params_converter import get_fortran_implicit_coefs
    from clubb_jax.src.CLUBB_core.pdf_params import init_pdf_implicit_coefs_terms_api
    nparams, ng, nz = 102, 2, 8
    rng = np.random.default_rng(11)
    # Plausible tunable params: slope>0, stdev_factor>0, coef_spread in (0,1].
    params = np.zeros((ng, nparams))
    params[:, 52] = 1.4    # islope_coef_spread_DG_means_w
    params[:, 53] = 1.2    # ipdf_component_stdev_factor_w
    params[:, 54] = 0.7    # icoef_spread_DG_means_rt
    params[:, 55] = 0.7    # icoef_spread_DG_means_thl
    worst = {}
    for _ in range(3):
        wm = rng.uniform(-1, 1, (ng, nz)); rtm = rng.uniform(0, 0.02, (ng, nz)); thlm = rng.uniform(285, 305, (ng, nz))
        wp2 = rng.uniform(0.05, 2, (ng, nz)); rtp2 = rng.uniform(1e-7, 1e-6, (ng, nz)); thlp2 = rng.uniform(0.05, 1, (ng, nz))
        Skw = rng.uniform(-2, 2, (ng, nz)); Skrt = rng.uniform(-2, 2, (ng, nz)); Skthl = rng.uniform(-2, 2, (ng, nz))
        wprtp = rng.uniform(-1e-3, 1e-3, (ng, nz)); wpthlp = rng.uniform(-0.5, 0.5, (ng, nz))
        rtpthlp = rng.uniform(-1e-4, 1e-4, (ng, nz))
        f = fortran_new_pdf_driver(nz, ng, wm, rtm, thlm, wp2, rtp2, thlp2, Skw, wprtp, wpthlp, rtpthlp,
                                   params, Skrt.copy(), Skthl.copy(),
                                   pdf_implicit_coefs_terms=init_pdf_implicit(nz, ng, 0))
        f_coefs = get_fortran_implicit_coefs()
        g = new_pdf_driver(wm, rtm, thlm, wp2, rtp2, thlp2, Skw, wprtp, wpthlp, rtpthlp, params, Skrt, Skthl,
                           init_pdf_implicit_coefs_terms_api(nz, ng, 0))
        pairs = list(zip(_DRIVER_OUTPUTS, f, (g[k] for k in _DRIVER_OUTPUTS)))
        for name in f_coefs._fields:
            ref = getattr(f_coefs, name)
            if ref is not None and np.asarray(ref).ndim == 2:
                pairs.append((name, ref, getattr(g["pdf_implicit_coefs_terms"], name)))
        for name, fi, gi in pairs:
            worst[name] = max(worst.get(name, 0.0), max_abs_diff(name, gi, fi))
    bad = {k: v for k, v in worst.items() if not v < 1e-11}
    assert not bad, f"new_pdf_driver Fortran mismatch: {bad}"
    assert len(worst) > len(_DRIVER_OUTPUTS), "no implicit coefficient fields were compared"
