#!/usr/bin/env python3
"""test_param_names.py — lock in the tunable-parameter name list + index mirror (parameters_tunable.py).

The order of `PARAM_NAMES` (returned by `get_param_names`) is load-bearing: every `i<name>` index constant in
parameter_indices.py (e.g. `ibeta=60`, `imu=59`, `iSkw_denom_coef=72`) is a zero-based offset into the `clubb_params`
array, so if the JAX name order ever drifted from the Fortran, every parameter access would silently read the wrong
value. The bit-faithful cases imply it's correct, but nothing pinned it directly. This test pins it two ways:

  1. **f2py order match** — `get_param_names()` equals `f2py_get_param_names()` name-for-name, in order (the Fortran
     `clubb_params` ordering is the authority). SKIPs if clubb_f2py is unbuilt.
  2. **index-constant self-consistency** — each `i<name>` constant indexes `PARAM_NAMES[i]` to exactly `<name>`
     (needs no oracle; catches a drifted index constant). (iter 465)
"""
import os
import sys

_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
for _p in (_ROOT, _ROOT + "/clubb_python_api"):
    if _p not in sys.path:
        sys.path.append(_p)

import numpy as np

from clubb_jax.src.CLUBB_core.parameters_tunable import (
    PARAMETER_HARD_BOUNDS,
    PNAME_IDX,
    check_parameters,
    get_param_names,
    get_parameter_hard_bounds,
    init_clubb_params,
)
from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core import parameter_indices

# Critical zero-based index constants → the PARAM_NAMES entry they must point at.
_INDEX_CONSTANTS = {
    "ibeta": "beta", "imu": "mu", "iSkw_denom_coef": "Skw_denom_coef",
    "ilambda0_stability_coef": "lambda0_stability_coef", "iC1": "C1", "iC2rt": "C2rt",
    "iC8": "C8", "iC11": "C11", "iC14": "C14", "igamma_coef": "gamma_coef",
}


def test_index_constant_self_consistency():
    """Each i<name> index constant maps PARAM_NAMES[i] == <name> (no oracle needed)."""
    names = list(get_param_names())
    bad = []
    for const, expected in _INDEX_CONSTANTS.items():
        idx = getattr(parameter_indices, const, None)
        if idx is None:
            bad.append(f"{const}: not defined in parameter_indices"); continue
        got = names[idx] if 0 <= idx < len(names) else f"<out of range {idx}>"
        if got != expected:
            bad.append(f"{const}={idx} -> PARAM_NAMES[{idx}]={got!r}, expected {expected!r}")
    assert not bad, "param-index constant(s) drifted from PARAM_NAMES:\n  " + "\n  ".join(bad)
    print(f"  {len(_INDEX_CONSTANTS)} index constants (ibeta/imu/iSkw_denom_coef/…) all map correctly  PASS")


def test_param_names_match_f2py():
    """get_param_names() equals the Fortran clubb_params name order, name-for-name."""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py param-names oracle: SKIP ({type(e).__name__})")
        return
    jax_names = list(get_param_names())
    n = int(clubb_f2py.f2py_get_nparams())
    raw = np.asarray(clubb_f2py.f2py_get_param_names(n))   # (nparams, namelen) of single chars
    f_names = [b"".join(raw[i]).decode("ascii", "ignore").strip() for i in range(raw.shape[0])]
    assert len(jax_names) == n, f"nparams mismatch: JAX {len(jax_names)} vs Fortran {n}"
    mism = [(i, jax_names[i], f_names[i]) for i in range(n) if jax_names[i] != f_names[i]]
    assert not mism, f"param-name order diverged from Fortran in {len(mism)} slot(s): {mism[:8]}"
    print(f"  get_param_names() matches f2py_get_param_names() exactly ({n} params, in order)  PASS")


def test_parameter_hard_bounds_match_f2py():
    """The full JAX hard-bound catalog matches Fortran slot-for-slot."""
    try:
        import clubb_f2py
    except Exception as e:
        print(f"  f2py parameter-hard-bounds oracle: SKIP ({type(e).__name__})")
        return
    jax_bounds = get_parameter_hard_bounds(len(get_param_names()))
    fortran_bounds = np.asarray(
        clubb_f2py.f2py_get_parameter_hard_bounds(len(get_param_names()))
    )
    assert np.array_equal(jax_bounds, fortran_bounds), "JAX hard-bound table diverged from Fortran"


def test_parameter_validation_rejects_hard_bound_and_equality():
    params = init_clubb_params(
        1,
        filename=os.path.join(
            _ROOT, "input", "parameter_and_flag_configs", "default", "tunable_parameters.in"
        ),
    )
    valid = check_parameters(1, params, ErrInfo.initialized(1))
    assert not valid.is_fatal_host()

    bad_bound = params.copy()
    bad_bound[0, PNAME_IDX["C_uu_shr"]] = PARAMETER_HARD_BOUNDS[1, PNAME_IDX["C_uu_shr"]] + 0.01
    assert check_parameters(1, bad_bound, ErrInfo.initialized(1)).is_fatal_host()

    bad_pair = params.copy()
    bad_pair[0, PNAME_IDX["C6thl"]] += 0.1
    assert check_parameters(1, bad_pair, ErrInfo.initialized(1)).is_fatal_host()


def test_param_values_match_fortran():
    """The JAX `init_clubb_params` VALUES (read from tunable_parameters.in) equal the Fortran's, slot-for-slot — so not
    only the NAMES/order (above) but the actual default tunable values (C1=2.4, beta, …) mirror the oracle. Uses the
    `clubb_python clubb_api.init_clubb_params` binding (the raw `f2py_init_clubb_params_file` mishandles the filename
    char-arg). ngrdcol=1 to avoid the binding's col>0 broadcast quirk (it doesn't replicate every param across columns;
    the JAX correctly does — col0 matches exactly either way). SKIPs if clubb_python is unbuilt. (iter 467)"""
    try:
        from clubb_python import clubb_api
    except Exception as e:
        print(f"  param-values vs Fortran: SKIP ({type(e).__name__})")
        return
    fn = os.path.join(_ROOT, "input", "tunable_parameters", "tunable_parameters.in")
    if not os.path.exists(fn):
        print("  param-values vs Fortran: SKIP (tunable_parameters.in absent)")
        return
    from clubb_jax.src.CLUBB_core.parameters_tunable import init_clubb_params
    j = np.asarray(init_clubb_params(1, filename=fn))
    f = np.asarray(clubb_api.init_clubb_params(1, 50, fn))
    assert j.shape == f.shape, f"params shape mismatch: JAX {j.shape} vs Fortran {f.shape}"
    worst = float(np.max(np.abs(j - f)))
    assert worst == 0.0, f"init_clubb_params values diverge from Fortran (max abs diff {worst:.2e})"
    print(f"  init_clubb_params values match Fortran exactly ({j.shape[1]} params, max diff 0.0)  PASS")


def test_calc_derrived_params_f2py():
    """`calc_derrived_params` (parameters_tunable.F90:calc_Derrived_Params_api) computes the DERIVED tunable params
    from the base clubb_params + grid: `lmin` (the mixing-length floor) and `mixt_frac_max_mag` (the PDF mixture-fraction
    cap) — both load-bearing (lmin clips Lscale; mixt_frac_max_mag bounds the ADG1 closure). Validate the two scalars vs
    `f2py_calc_derrived_params`. SKIPs if clubb_f2py/clubb_python are unbuilt. (iter 468)"""
    try:
        import clubb_f2py
        from clubb_python import clubb_api
        from clubb_python.derived_types.err_info import ErrInfo
    except Exception as e:
        print(f"  f2py calc_derrived_params oracle: SKIP ({type(e).__name__})")
        return
    fn = os.path.join(_ROOT, "input", "tunable_parameters", "tunable_parameters.in")
    if not os.path.exists(fn):
        print("  f2py calc_derrived_params: SKIP (tunable_parameters.in absent)")
        return
    from clubb_jax.src.CLUBB_core.grid_class import setup_grid
    from clubb_jax.src.CLUBB_core.parameters_tunable import calc_derrived_params, init_clubb_params
    NG, DZ, ZTOP = 2, 40.0, 1200.0
    jgr = setup_grid(ngrdcol=NG, deltaz=DZ, zm_init=0.0, zm_top=ZTOP, grid_type=1)
    ng, nzm = jgr.zm.shape
    clubb_api.init_err_info(ng); cf = clubb_api.get_default_config_flags(); clubb_api.init_config_flags(cf)
    clubb_api.setup_grid(nzmax=nzm, ngrdcol=ng, sfc_elevation=np.zeros(ng), l_implemented=False,
                         l_ascending_grid=True, grid_type=2, deltaz=np.full(ng, DZ), zm_init=np.zeros(ng),
                         zm_top=np.full(ng, float(jgr.zm[0, -1])),
                         momentum_heights=np.asfortranarray(np.asarray(jgr.zm)),
                         thermodynamic_heights=np.asfortranarray(np.asarray(jgr.zt)), err_info=ErrInfo(ngrdcol=ng))
    cp = np.asarray(init_clubb_params(ng, filename=fn))
    deltaz = np.full(ng, DZ)
    _, j_lmin, j_mfmm = calc_derrived_params(jgr, ng, 1, deltaz, cp, False)
    f_lmin, f_mfmm = clubb_f2py.f2py_calc_derrived_params(1, deltaz, cp, 0)
    d_lmin = float(np.max(np.abs(np.asarray(j_lmin) - np.asarray(f_lmin))))
    d_mfmm = abs(float(j_mfmm) - float(np.asarray(f_mfmm).ravel()[0]))
    assert d_lmin < 1e-12 and d_mfmm < 1e-12, f"calc_derrived_params mismatch: lmin {d_lmin:.2e}, mixt_frac_max_mag {d_mfmm:.2e}"
    print(f"  f2py calc_derrived_params: lmin + mixt_frac_max_mag match exactly (diff {max(d_lmin, d_mfmm):.1e})  PASS")


def main():
    print("test_param_names:")
    test_index_constant_self_consistency()
    test_param_names_match_f2py()
    test_parameter_hard_bounds_match_f2py()
    test_parameter_validation_rejects_hard_bound_and_equality()
    test_param_values_match_fortran()
    test_calc_derrived_params_f2py()
    print("All param-name checks PASSED")


if __name__ == "__main__":
    main()
