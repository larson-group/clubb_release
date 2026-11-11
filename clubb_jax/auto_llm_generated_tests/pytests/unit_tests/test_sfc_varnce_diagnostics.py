"""Compiled surface diagnostics mirror sfc_varnce_module.F90's failure dump.

Checks runtime-only column selection, scalar fields, complete arrays, and the
source debug gate. Physics regression comparisons remain separate workflows.
"""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core.error_code import set_debug_level
from clubb_jax.src.CLUBB_core.grid_class import setup_grid
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.sclr_idx import SclrIdx
from clubb_jax.src.CLUBB_core.sfc_varnce_module import calc_sfc_varnce


@pytest.fixture(autouse=True)
def diagnostic_debug_level():
    from clubb_jax.src.CLUBB_core import error_code

    original_level = error_code._debug_level
    set_debug_level(2)
    calc_sfc_varnce.clear_cache()
    yield
    jax.effects_barrier()
    set_debug_level(original_level)
    calc_sfc_varnce.clear_cache()


def _inputs(sclr_dim=0, nlevels=5, coordinates=True):
    gr = setup_grid(2, 100.0, 0.0, (nlevels - 1) * 100.0)
    zm = jnp.ones((2, gr.nzm))
    zt = jnp.ones((2, gr.nzt))
    scalars = jnp.ones((2, gr.nzm, sclr_dim))
    err_info = ErrInfo.initialized(
        2, chunk_idx=3, mpi_rank=7,
        lat=jnp.array([10.0, 20.0]) if coordinates else None,
        lon=jnp.array([-30.0, -40.0]) if coordinates else None,
    )
    return dict(
        nzm=gr.nzm, nzt=gr.nzt, ngrdcol=2, sclr_dim=sclr_dim,
        sclr_idx=SclrIdx(1, 2, 0, 0, 0, 0), gr=gr, dt=1.0,
        sfc_elevation=jnp.zeros(2), upwp_sfc=jnp.full(2, 0.1),
        vpwp_sfc=jnp.full(2, 0.2), wpthlp=zm * 0.01,
        wprtp_sfc=jnp.full(2, 0.001), um=zt, vm=zt, Lscale_up=zt,
        wpsclrp_sfc=jnp.full((2, sclr_dim), 0.001),
        lhs_splat_wp2=zm * 0.0, tau_zm=zm,
        l_vary_convect_depth=False, T0=300.0,
        up2_sfc_coef=jnp.ones(2), a_const=jnp.full(2, 1.8),
        stats=JaxStats.empty(l_sample=False, names=(), ncol=2, max_nlev=gr.nzm),
        wp2=zm, up2=zm, vp2=zm, thlp2=zm, rtp2=zm, rtpthlp=zm,
        sclrp2=scalars, sclrprtp=scalars, sclrpthlp=scalars, err_info=err_info,
    )


def _run(inputs, capsys):
    result = calc_sfc_varnce(**inputs)
    jax.block_until_ready(result)
    jax.effects_barrier()
    captured = capsys.readouterr()
    assert captured.out == ""
    return result, captured.err


@pytest.mark.parametrize("sclr_dim", [0, 2])
def test_compiled_dump_selects_fatal_columns_and_preserves_source_fields(sclr_dim, capsys):
    inputs = _inputs(sclr_dim)
    # Healthy first call compiles the same signature without printing.
    _, text = _run(inputs, capsys)
    assert text == ""
    inputs["upwp_sfc"] = inputs["upwp_sfc"].at[1].set(jnp.nan)
    for _ in range(2):
        result, text = _run(inputs, capsys)
        np.testing.assert_array_equal(result[-2].err_code, [0, 99])
        assert "Grid column index i = 1," not in text
        assert "Grid column index i = 2, Latitude 20.00 / Longitude -40.00" in text
        assert "MPI Process 7 / Chunk index 3" in text
        labels = ["Error in calc_sfc_varnce", "Intent(in)", "upwp_sfc = ",
                  "vpwp_sfc = ", "wpthlp = ", "wprtp_sfc = "]
        if sclr_dim:
            labels += ["wpsclrp = "]
        labels += ["Intent(out)", "wp2 = ", "up2 = ", "vp2 = ",
                   "thlp2 = ", "rtp2 = ", "rtpthlp = "]
        if sclr_dim:
            labels += ["sclrp2 = ", "sclrprtp = ", "sclrpthlp = "]
        else:
            assert "wpsclrp = " not in text and "sclrp2 = " not in text
        positions = [text.index(label) for label in labels]
        assert positions == sorted(positions)
        with np.printoptions(threshold=np.inf):
            assert "wp2 =  " + str(np.asarray(result[0])) in text


@pytest.mark.parametrize("debug_level", [-1, 0, 1])
def test_lower_debug_levels_suppress_fatal_dump(debug_level, capsys):
    set_debug_level(debug_level)
    inputs = _inputs()
    inputs["err_info"] = inputs["err_info"].set_fatal(jnp.array([False, True]))
    result, text = _run(inputs, capsys)
    np.testing.assert_array_equal(result[-2].err_code, [0, 99])
    assert text == ""


def test_scalar_only_nan_and_existing_fatal_headers_use_column_order(capsys):
    inputs = _inputs(2, coordinates=False)
    inputs["wpsclrp_sfc"] = inputs["wpsclrp_sfc"].at[1, 1].set(jnp.nan)
    inputs["err_info"] = inputs["err_info"].set_fatal(jnp.array([True, False]))
    result, text = _run(inputs, capsys)
    np.testing.assert_array_equal(result[-2].err_code, [99, 99])
    assert text.index("Grid column index i = 1,") < text.index("Grid column index i = 2,")
    assert "Latitude -999.00 / Longitude -999.00" in text
    assert "sclrp2 = " in text and "nan" in text


def test_large_field_dump_does_not_abbreviate_arrays(capsys):
    inputs = _inputs(nlevels=601)
    inputs["err_info"] = inputs["err_info"].set_fatal(jnp.array([False, True]))
    inputs["wp2"] = inputs["wp2"].at[0, 300].set(12345.0)
    result, text = _run(inputs, capsys)
    assert "..." not in text
    with np.printoptions(threshold=np.inf):
        assert "wp2 =  " + str(np.asarray(result[0])) in text
