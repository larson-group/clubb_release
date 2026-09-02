"""f2py advance_xm_wpxp oracle — direct-`.so`-call smoke test (Iter147).

The `clubb_python` Python wrapper (`clubb_python/CLUBB_core/advance_xm_wpxp.py`) is
OUT OF SYNC with its compiled `clubb_f2py.*.so`: the wrapper passes `wp3, kh_zt`
where the recompiled `.so` now expects `wp3_on_wp2, wp3_on_wp2_zt, kh_zt` plus an
extra `skw_zm`, so `clubb_python`'s own `test_call_tree_advance_xm_wpxp` raises a
TypeError. This made the f2py oracle look "blocked" for several iterations.

It is NOT blocked: `clubb_f2py.f2py_advance_xm_wpxp.__doc__` gives the exact `.so`
signature, and the routine is callable DIRECTLY (after pushing the UDTs via the
`set_fortran_*` converters). This test proves that, and is the template for the
definitive rico stretched-grid **input-matched comparison**: capture rico's matched
step-1 `advance_xm_wpxp` inputs (eager) and feed them here, then diff the Fortran
vs JAX outputs to localise the grid_type=2 divergence (see DESIGN "rico diagnosis").

Requires the clubb_python_api build on PYTHONPATH:
  export PYTHONPATH=$PWD:$PWD/clubb_python_api
"""
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

# The .so signature (from clubb_f2py.f2py_advance_xm_wpxp.__doc__), in order — the
# 3 args the stale wrapper omitted are marked. Returned: the advanced state arrays.
_SO_NEW_ARGS = ("wp3_on_wp2", "wp3_on_wp2_zt", "skw_zm")  # absent from the old wrapper

_RET_NAMES = ("wprtp_cl_num wpthlp_cl_num upwp_cl_num vpwp_cl_num rtm wprtp thlm wpthlp "
              "sclrm wpsclrp um upwp vm vpwp um_pert vm_pert upwp_pert vpwp_pert").split()


def call_f2py_advance_xm_wpxp(gr, nu, pdf_ic, err, a):
    """Call the f2py advance_xm_wpxp `.so` DIRECTLY (bypassing the broken wrapper).

    `a` is the full argument dict (see clubb_python test `_make_args`) PLUS the three
    keys in `_SO_NEW_ARGS`. Returns a dict {name: array} of the advanced outputs.
    """
    import clubb_f2py
    from clubb_python.derived_types.grid_class_converter import set_fortran_grid
    from clubb_python.derived_types.nu_vert_res_dep_converter import set_fortran_nu_vert_res_dep
    from clubb_python.derived_types.pdf_params_converter import set_fortran_implicit_coefs
    from clubb_python.derived_types.err_info_converter import set_fortran_err_info

    set_fortran_grid(gr)
    set_fortran_nu_vert_res_dep(nu)
    set_fortran_implicit_coefs(pdf_ic)
    set_fortran_err_info(err)

    result = clubb_f2py.f2py_advance_xm_wpxp(
        a["sclr_dim"], a["sclr_tol"], a["dt"], a["sigma_sqd_w"], a["wm_zm"], a["wm_zt"],
        a["wp2"], a["lscale_zm"], a["wp3_on_wp2"], a["wp3_on_wp2_zt"], a["kh_zt"], a["kh_zm"],
        a["stability_correction"], a["invrs_tau_c6_zm"], a["tau_max_zm"], a["skw_zm"], a["wp2rtp"],
        a["rtpthvp"], a["rtm_forcing"], a["wprtp_forcing"], a["rtm_ref"], a["wp2thlp"],
        a["thlpthvp"], a["thlm_forcing"], a["wpthlp_forcing"], a["thlm_ref"], a["rho_ds_zm"],
        a["rho_ds_zt"], a["invrs_rho_ds_zm"], a["invrs_rho_ds_zt"], a["thv_ds_zm"], a["rtp2"],
        a["thlp2"], a["w_1_zm"], a["w_2_zm"], a["varnce_w_1_zm"], a["varnce_w_2_zm"],
        a["mixt_frac_zm"], a["l_implemented"], a["em"], a["wp2sclrp"], a["sclrpthvp"],
        a["sclrm_forcing"], a["sclrp2"], a["cx_fnc_richardson"], a["um_forcing"], a["vm_forcing"],
        a["ug"], a["vg"], a["wpthvp"], a["fcor"], a["fcor_y"], a["um_ref"], a["vm_ref"],
        a["up2"], a["vp2"], a["uprcp"], a["vprcp"], a["rc_coef_zm"], a["clubb_params"],
        a["ts_nudge"], a["iipdf_type"], a["penta_solve_method"], a["tridiag_solve_method"],
        a["fill_holes_type"], a["l_predict_upwp_vpwp"], a["l_ho_nontrad_coriolis"],
        a["l_ho_trad_coriolis"], a["l_diffuse_rtm_and_thlm"], a["l_stability_correct_kh_n2_zm"],
        a["l_godunov_upwind_wpxp_ta"], a["l_upwind_xm_ma"], a["l_uv_nudge"],
        a["l_diag_lscale_from_tau"], a["l_use_c7_richardson"], a["l_lmm_stepping"],
        a["l_enable_relaxed_clipping"], a["l_linearize_pbl_winds"], a["l_mono_flux_lim_thlm"],
        a["l_mono_flux_lim_rtm"], a["l_mono_flux_lim_um"], a["l_mono_flux_lim_vm"],
        a["l_mono_flux_lim_spikefix"], a["wprtp_cl_num"], a["wpthlp_cl_num"], a["upwp_cl_num"],
        a["vpwp_cl_num"], a["rtm"], a["wprtp"], a["thlm"], a["wpthlp"], a["sclrm"], a["wpsclrp"],
        a["um"], a["upwp"], a["vm"], a["vpwp"], a["um_pert"], a["vm_pert"], a["upwp_pert"],
        a["vpwp_pert"], nzm=a["nzm"], nzt=a["nzt"], ngrdcol=a["ngrdcol"],
    )
    return {nm: np.asarray(arr) for nm, arr in zip(_RET_NAMES, result)}


def test_f2py_so_directly_callable():
    """The f2py advance_xm_wpxp `.so` runs and returns finite advanced arrays."""
    api_root = Path(__file__).resolve().parents[2] / "clubb_python_api"
    if not (api_root / "clubb_f2py.pyf").exists():
        print("SKIP: clubb_python_api not present")
        return
    sys.path.insert(0, str(api_root))
    sys.path.insert(0, str(api_root / "clubb_python"))
    sys.path.insert(0, str(api_root / "tests"))
    import test_call_tree_advance_xm_wpxp as T

    tmp = Path(tempfile.mkdtemp())
    gr, flags, clubb_params, nu, pdf_ic, err = T._setup_env(tmp, ngrdcol=1)
    a = T._make_args(gr, flags, clubb_params, nu, pdf_ic, err)
    ng, nzt, nzm = gr.ngrdcol, gr.nzt, gr.nzm

    def full(shape, val):
        return np.full(shape, val, dtype=np.float64, order="F")

    # the 3 inputs the recompiled .so needs that the old wrapper omits
    a["wp3_on_wp2"] = full((ng, nzm), 0.0)
    a["wp3_on_wp2_zt"] = full((ng, nzt), 0.0)
    a["skw_zm"] = full((ng, nzm), 0.0)

    out = call_f2py_advance_xm_wpxp(gr, nu, pdf_ic, err, a)
    for nm in ("rtm", "wprtp", "thlm", "wpthlp", "um", "upwp", "vm", "vpwp"):
        assert np.all(np.isfinite(out[nm])), f"{nm} not finite"
    print("PASS f2py advance_xm_wpxp .so is directly callable (oracle UNBLOCKED)")


if __name__ == "__main__":
    test_f2py_so_directly_callable()
