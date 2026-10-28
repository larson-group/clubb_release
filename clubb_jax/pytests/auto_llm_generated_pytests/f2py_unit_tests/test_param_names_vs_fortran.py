"""Check JAX tunable-parameter names, bounds, defaults and derived parameters against Fortran."""

from input.scripts.create_case_namelist import create_case_namelist_file

import os

import numpy as np

from utilities.output_paths import REPO_ROOT as _REPO_ROOT

from clubb_jax.pytests.f2py_unit_test_helpers import max_abs_diff, setup_matching_grids
from clubb_jax.src.CLUBB_core.parameters_tunable import (
    get_param_names,
    get_parameter_hard_bounds,
    init_clubb_params,
)

_ROOT = str(_REPO_ROOT)


def test_param_names_match_f2py():
    """get_param_names() equals the Fortran clubb_params name order, name-for-name."""
    import clubb_f2py
    jax_names = list(get_param_names())
    n = int(clubb_f2py.f2py_get_nparams())
    raw = np.asarray(clubb_f2py.f2py_get_param_names(n))   # (nparams, namelen) of single chars
    f_names = [b"".join(raw[i]).decode("ascii", "ignore").strip() for i in range(raw.shape[0])]
    assert len(jax_names) == n, f"nparams mismatch: JAX {len(jax_names)} vs Fortran {n}"
    mism = [(i, jax_names[i], f_names[i]) for i in range(n) if jax_names[i] != f_names[i]]
    assert not mism, f"param-name order diverged from Fortran in {len(mism)} slot(s): {mism[:8]}"


def test_parameter_hard_bounds_match_f2py():
    """The full JAX hard-bound catalog matches Fortran slot-for-slot."""
    import clubb_f2py
    jax_bounds = get_parameter_hard_bounds(len(get_param_names()))
    fortran_bounds = np.asarray(
        clubb_f2py.f2py_get_parameter_hard_bounds(len(get_param_names()))
    )
    assert np.array_equal(jax_bounds, fortran_bounds), "JAX hard-bound table diverged from Fortran"


def test_param_values_match_fortran():
    """The JAX `init_clubb_params` VALUES (read from tunable_parameters.in) equal the Fortran's, slot-for-slot — so not
    only the NAMES/order (above) but the actual default tunable values (C1=2.4, beta, …) mirror the oracle. Uses the
    `clubb_python clubb_api.init_clubb_params` binding (the raw `f2py_init_clubb_params_file` mishandles the filename
    char-arg). ngrdcol=1 to avoid the binding's col>0 broadcast quirk (it doesn't replicate every param across columns;
    the JAX correctly does — col0 matches exactly either way)."""
    from clubb_python import clubb_api
    fn = str(create_case_namelist_file("bomex"))
    if not os.path.exists(fn):
        raise AssertionError("  param-values vs Fortran: required input missing (tunable_parameters.in absent)")
    from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params
    j = np.asarray(init_clubb_params(1, filename=fn))
    f = np.asarray(clubb_api.init_clubb_params(1, 50, fn))
    assert j.shape == f.shape, f"params shape mismatch: JAX {j.shape} vs Fortran {f.shape}"
    worst = float(max_abs_diff("param_names", j, f))
    assert worst == 0.0, f"init_clubb_params values diverge from Fortran (max abs diff {worst:.2e})"


def test_calc_derived_params_f2py():
    """`calc_derived_params` (parameters_tunable.F90:calc_derived_params_api) computes the derived tunable params
    from the base clubb_params + grid: `lmin` (the mixing-length floor) and `mixt_frac_max_mag` (the PDF mixture-fraction
    cap) — both load-bearing (lmin clips Lscale; mixt_frac_max_mag bounds the ADG1 closure). Validate the two scalars vs
    `f2py_calc_derived_params`."""
    from clubb_python.CLUBB_core.parameters_tunable import calc_derived_params as fortran_calc_derived_params
    from clubb_python.derived_types.nu_vert_res_dep import NuVertResDep
    fn = str(create_case_namelist_file("bomex"))
    if not os.path.exists(fn):
        raise AssertionError("  f2py calc_derived_params: required input missing (tunable_parameters.in absent)")
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    from clubb_jax.src.CLUBB_core.parameters_tunable import calc_derived_params, init_clubb_params
    NG, DZ, ZTOP = 2, 40.0, 1200.0
    jgr, fgr = setup_matching_grids(NG, DZ, ZTOP)
    ng, nzm = jgr.zm.shape
    cp = np.asarray(init_clubb_params(ng, filename=fn))
    deltaz = np.full(ng, DZ)
    _, j_lmin, j_mfmm = calc_derived_params(jgr, ng, 1, deltaz, cp, False)
    _, f_lmin, f_mfmm = fortran_calc_derived_params(fgr, ng, 1, deltaz, cp, False, NuVertResDep.zeros(ng))
    d_lmin = float(max_abs_diff("param_names", np.asarray(j_lmin), np.asarray(f_lmin)))
    d_mfmm = float(max_abs_diff("param_names", float(j_mfmm), f_mfmm))
    assert d_lmin < 1e-12 and d_mfmm < 1e-12, f"calc_derived_params mismatch: lmin {d_lmin:.2e}, mixt_frac_max_mag {d_mfmm:.2e}"
