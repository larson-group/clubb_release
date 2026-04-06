"""Test wrappers for numerical_check call-tree routines."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo


def test_calculate_spurious_source_matches_formula():
    """Spurious source should match the conservation residual formula."""
    integral_after = 10.0
    integral_before = 6.0
    flux_top = 0.5
    flux_sfc = 0.2
    integral_forcing = 0.1
    dt = 2.0

    got = clubb_api.calculate_spurious_source(
        integral_after, integral_before, flux_top, flux_sfc, integral_forcing, dt)
    expected = (integral_after - integral_before) / dt + flux_top - flux_sfc - integral_forcing
    np.testing.assert_allclose(got, expected)


def test_sfc_varnce_check_sets_error_for_nan(run_quiet_python):
    """sfc_varnce_check should leave finite inputs clean and flag NaN inputs."""
    ngrdcol = 1
    clubb_api.init_err_info(ngrdcol)
    err_info = ErrInfo(ngrdcol=ngrdcol)
    clubb_api.reset_err_code()

    clubb_api.sfc_varnce_check(
        sclr_dim=1,
        wp2_sfc=0.1, up2_sfc=0.2, vp2_sfc=0.3, thlp2_sfc=0.4,
        rtp2_sfc=0.5, rtpthlp_sfc=0.6,
        sclrp2_sfc=np.array([0.7], dtype=np.float64),
        sclrprtp_sfc=np.array([0.8], dtype=np.float64),
        sclrpthlp_sfc=np.array([0.9], dtype=np.float64),
        err_info=err_info,
    )
    assert int(clubb_api.get_err_code(ngrdcol)[0]) == 0

    proc = run_quiet_python(
        """
import json
import numpy as np
from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo

ngrdcol = 1
clubb_api.init_err_info(ngrdcol)
err_info = ErrInfo(ngrdcol=ngrdcol)
clubb_api.reset_err_code()
clubb_api.sfc_varnce_check(
    sclr_dim=1,
    wp2_sfc=np.nan, up2_sfc=0.2, vp2_sfc=0.3, thlp2_sfc=0.4,
    rtp2_sfc=0.5, rtpthlp_sfc=0.6,
    sclrp2_sfc=np.array([0.7], dtype=np.float64),
    sclrprtp_sfc=np.array([0.8], dtype=np.float64),
    sclrpthlp_sfc=np.array([0.9], dtype=np.float64),
    err_info=err_info,
)
print(json.dumps({"err_code": int(clubb_api.get_err_code(ngrdcol)[0])}))
"""
    )
    assert proc.returncode == 0, proc.stderr
    assert json.loads(proc.stdout)["err_code"] != 0


def test_parameterization_check_sets_error_for_nan(run_quiet_python):
    """parameterization_check should set err_code when given NaN state input."""
    ngrdcol = 1
    clubb_api.init_err_info(ngrdcol)
    err_info = ErrInfo(ngrdcol=ngrdcol)

    nzt = 6
    nzm = 7
    sclr_dim = 1
    edsclr_dim = 1

    args = dict(
        nzm=nzm, nzt=nzt, sclr_dim=sclr_dim, edsclr_dim=edsclr_dim,
        thlm_forcing=np.zeros(nzt, dtype=np.float64),
        rtm_forcing=np.zeros(nzt, dtype=np.float64),
        um_forcing=np.zeros(nzt, dtype=np.float64),
        vm_forcing=np.zeros(nzt, dtype=np.float64),
        wm_zm=np.zeros(nzm, dtype=np.float64),
        wm_zt=np.zeros(nzt, dtype=np.float64),
        p_in_pa=np.full(nzt, 1.0e5, dtype=np.float64),
        rho_zm=np.ones(nzm, dtype=np.float64),
        rho=np.ones(nzt, dtype=np.float64),
        exner=np.ones(nzt, dtype=np.float64),
        rho_ds_zm=np.ones(nzm, dtype=np.float64),
        rho_ds_zt=np.ones(nzt, dtype=np.float64),
        invrs_rho_ds_zm=np.ones(nzm, dtype=np.float64),
        invrs_rho_ds_zt=np.ones(nzt, dtype=np.float64),
        thv_ds_zm=np.full(nzm, 300.0, dtype=np.float64),
        thv_ds_zt=np.full(nzt, 300.0, dtype=np.float64),
        wpthlp_sfc=0.0,
        wprtp_sfc=0.0,
        upwp_sfc=0.0,
        vpwp_sfc=0.0,
        p_sfc=1.0e5,
        um=np.zeros(nzt, dtype=np.float64),
        upwp=np.zeros(nzm, dtype=np.float64),
        vm=np.zeros(nzt, dtype=np.float64),
        vpwp=np.zeros(nzm, dtype=np.float64),
        up2=np.full(nzm, 0.1, dtype=np.float64),
        vp2=np.full(nzm, 0.1, dtype=np.float64),
        rtm=np.full(nzt, 0.01, dtype=np.float64),
        wprtp=np.zeros(nzm, dtype=np.float64),
        thlm=np.full(nzt, 300.0, dtype=np.float64),
        wpthlp=np.zeros(nzm, dtype=np.float64),
        wp2=np.full(nzm, 0.1, dtype=np.float64),
        wp3=np.zeros(nzt, dtype=np.float64),
        rtp2=np.full(nzm, 1.0e-4, dtype=np.float64),
        thlp2=np.full(nzm, 0.2, dtype=np.float64),
        rtpthlp=np.zeros(nzm, dtype=np.float64),
        prefix="middle of ",
        wpsclrp_sfc=np.zeros(sclr_dim, dtype=np.float64),
        wpedsclrp_sfc=np.zeros(edsclr_dim, dtype=np.float64),
        sclrm=np.zeros((nzt, sclr_dim), dtype=np.float64),
        wpsclrp=np.zeros((nzm, sclr_dim), dtype=np.float64),
        sclrp2=np.zeros((nzm, sclr_dim), dtype=np.float64),
        sclrprtp=np.zeros((nzm, sclr_dim), dtype=np.float64),
        sclrpthlp=np.zeros((nzm, sclr_dim), dtype=np.float64),
        sclrm_forcing=np.zeros((nzt, sclr_dim), dtype=np.float64),
        edsclrm=np.zeros((nzt, edsclr_dim), dtype=np.float64),
        edsclrm_forcing=np.zeros((nzt, edsclr_dim), dtype=np.float64),
        err_info=err_info,
    )

    clubb_api.reset_err_code()
    clubb_api.parameterization_check(**args)
    assert int(clubb_api.get_err_code(ngrdcol)[0]) == 0

    proc = run_quiet_python(
        """
import json
import numpy as np
from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo

ngrdcol = 1
clubb_api.init_err_info(ngrdcol)
err_info = ErrInfo(ngrdcol=ngrdcol)

nzt = 6
nzm = 7
sclr_dim = 1
edsclr_dim = 1

args = dict(
    nzm=nzm, nzt=nzt, sclr_dim=sclr_dim, edsclr_dim=edsclr_dim,
    thlm_forcing=np.zeros(nzt, dtype=np.float64),
    rtm_forcing=np.zeros(nzt, dtype=np.float64),
    um_forcing=np.zeros(nzt, dtype=np.float64),
    vm_forcing=np.zeros(nzt, dtype=np.float64),
    wm_zm=np.zeros(nzm, dtype=np.float64),
    wm_zt=np.zeros(nzt, dtype=np.float64),
    p_in_pa=np.full(nzt, 1.0e5, dtype=np.float64),
    rho_zm=np.ones(nzm, dtype=np.float64),
    rho=np.ones(nzt, dtype=np.float64),
    exner=np.ones(nzt, dtype=np.float64),
    rho_ds_zm=np.ones(nzm, dtype=np.float64),
    rho_ds_zt=np.ones(nzt, dtype=np.float64),
    invrs_rho_ds_zm=np.ones(nzm, dtype=np.float64),
    invrs_rho_ds_zt=np.ones(nzt, dtype=np.float64),
    thv_ds_zm=np.full(nzm, 300.0, dtype=np.float64),
    thv_ds_zt=np.full(nzt, 300.0, dtype=np.float64),
    wpthlp_sfc=0.0,
    wprtp_sfc=0.0,
    upwp_sfc=0.0,
    vpwp_sfc=0.0,
    p_sfc=1.0e5,
    um=np.zeros(nzt, dtype=np.float64),
    upwp=np.zeros(nzm, dtype=np.float64),
    vm=np.zeros(nzt, dtype=np.float64),
    vpwp=np.zeros(nzm, dtype=np.float64),
    up2=np.full(nzm, 0.1, dtype=np.float64),
    vp2=np.full(nzm, 0.1, dtype=np.float64),
    rtm=np.full(nzt, 0.01, dtype=np.float64),
    wprtp=np.zeros(nzm, dtype=np.float64),
    thlm=np.full(nzt, 300.0, dtype=np.float64),
    wpthlp=np.zeros(nzm, dtype=np.float64),
    wp2=np.full(nzm, 0.1, dtype=np.float64),
    wp3=np.zeros(nzt, dtype=np.float64),
    rtp2=np.full(nzm, 1.0e-4, dtype=np.float64),
    thlp2=np.full(nzm, 0.2, dtype=np.float64),
    rtpthlp=np.zeros(nzm, dtype=np.float64),
    prefix="middle of ",
    wpsclrp_sfc=np.zeros(sclr_dim, dtype=np.float64),
    wpedsclrp_sfc=np.zeros(edsclr_dim, dtype=np.float64),
    sclrm=np.zeros((nzt, sclr_dim), dtype=np.float64),
    wpsclrp=np.zeros((nzm, sclr_dim), dtype=np.float64),
    sclrp2=np.zeros((nzm, sclr_dim), dtype=np.float64),
    sclrprtp=np.zeros((nzm, sclr_dim), dtype=np.float64),
    sclrpthlp=np.zeros((nzm, sclr_dim), dtype=np.float64),
    sclrm_forcing=np.zeros((nzt, sclr_dim), dtype=np.float64),
    edsclrm=np.zeros((nzt, edsclr_dim), dtype=np.float64),
    edsclrm_forcing=np.zeros((nzt, edsclr_dim), dtype=np.float64),
    err_info=err_info,
)

clubb_api.reset_err_code()
args["thlm"] = args["thlm"].copy()
args["thlm"][0] = np.nan
clubb_api.parameterization_check(**args)
print(json.dumps({"err_code": int(clubb_api.get_err_code(ngrdcol)[0])}))
"""
    )
    assert proc.returncode == 0, proc.stderr
    assert json.loads(proc.stdout)["err_code"] != 0


def test_pdf_closure_check_sets_error_for_nan(run_quiet_python):
    """pdf_closure_check should pass finite inputs and flag NaN in required terms."""
    ngrdcol = 1
    nz = 5
    sclr_dim = 1

    clubb_api.init_err_info(ngrdcol)
    err_info = ErrInfo(ngrdcol=ngrdcol)
    pdf_params = clubb_api.init_pdf_params(nz, ngrdcol)

    base = np.full(nz, 0.1, dtype=np.float64)
    sclr_base = np.full((nz, sclr_dim), 0.01, dtype=np.float64)

    clubb_api.reset_err_code()
    clubb_api.pdf_closure_check(
        nz, sclr_dim,
        base, base, base, base, base, base, base, base, base, base,
        base, base, base, base, base, base, base, base, base, base,
        base, base, sclr_base, sclr_base, sclr_base, sclr_base, sclr_base, sclr_base,
        pdf_params=pdf_params, err_info=err_info,
    )
    assert int(clubb_api.get_err_code(ngrdcol)[0]) == 0

    proc = run_quiet_python(
        """
import json
import numpy as np
from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo

ngrdcol = 1
nz = 5
sclr_dim = 1
clubb_api.init_err_info(ngrdcol)
err_info = ErrInfo(ngrdcol=ngrdcol)
pdf_params = clubb_api.init_pdf_params(nz, ngrdcol)

base = np.full(nz, 0.1, dtype=np.float64)
sclr_base = np.full((nz, sclr_dim), 0.01, dtype=np.float64)
bad_wp2rtp = base.copy()
bad_wp2rtp[0] = np.nan

clubb_api.reset_err_code()
clubb_api.pdf_closure_check(
    nz, sclr_dim,
    base, base, bad_wp2rtp, base, base, base, base, base, base, base,
    base, base, base, base, base, base, base, base, base, base,
    base, base, sclr_base, sclr_base, sclr_base, sclr_base, sclr_base, sclr_base,
    pdf_params=pdf_params, err_info=err_info,
)
print(json.dumps({"err_code": int(clubb_api.get_err_code(ngrdcol)[0])}))
"""
    )
    assert proc.returncode == 0, proc.stderr
    assert json.loads(proc.stdout)["err_code"] != 0


def test_length_check_sets_error_for_nan(run_quiet_python):
    """length_check should set err_code when a mixing-length array contains NaN."""
    ngrdcol = 1
    clubb_api.init_err_info(ngrdcol)
    err_info = ErrInfo(ngrdcol=ngrdcol)

    clubb_api.reset_err_code()
    clubb_api.length_check(
        nzt=4,
        lscale=np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64),
        lscale_up=np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64),
        lscale_down=np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64),
        err_info=err_info,
    )
    assert int(clubb_api.get_err_code(ngrdcol)[0]) == 0

    proc = run_quiet_python(
        """
import json
import numpy as np
from clubb_python import clubb_api
from clubb_python.derived_types.err_info import ErrInfo

ngrdcol = 1
clubb_api.init_err_info(ngrdcol)
err_info = ErrInfo(ngrdcol=ngrdcol)
clubb_api.reset_err_code()
clubb_api.length_check(
    nzt=4,
    lscale=np.array([1.0, np.nan, 3.0, 4.0], dtype=np.float64),
    lscale_up=np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64),
    lscale_down=np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float64),
    err_info=err_info,
)
print(json.dumps({"err_code": int(clubb_api.get_err_code(ngrdcol)[0])}))
"""
    )
    assert proc.returncode == 0, proc.stderr
    assert json.loads(proc.stdout)["err_code"] != 0
