"""Runtime diagnostics restored from CLUBB_core Fortran source.

Synthetic failures verify host callbacks, source index order, full dumps and
unchanged failure status without needing a full atmospheric case.
"""

import io

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src.CLUBB_core.err_info import ErrInfo
from clubb_jax.src.CLUBB_core import error_code
from clubb_jax.src.CLUBB_core.numerical_check import check_nan, check_negative
from clubb_jax.src.CLUBB_core.matrix_solver_wrapper import (
    check_nan_soln_1d, check_nan_soln_2d, check_nan_soln_3d,
)
from clubb_jax.src.CLUBB_core.model_flags import get_default_config_flags, print_clubb_config_flags_api
from clubb_jax.src.CLUBB_core.pdf_params import init_pdf_params, print_pdf_params
from clubb_jax.src.CLUBB_core.clip_explicit import clip_rcm


@pytest.fixture(autouse=True)
def debug_configuration():
    original = error_code._debug_level
    error_code.set_debug_level(3)
    clip_rcm.clear_cache()
    yield
    jax.effects_barrier()
    error_code.set_debug_level(original)
    clip_rcm.clear_cache()


def _ready(result, capsys):
    jax.block_until_ready(result)
    jax.effects_barrier()
    return capsys.readouterr()


def test_negative_runtime_prints_preserve_level_then_column_order(capsys):
    info = ErrInfo.initialized(2, mpi_rank=5, chunk_idx=3)
    run = jax.jit(lambda values: check_negative(values, "rtm", "advance_clubb_core", info))
    healthy = run(jnp.ones((2, 3)))
    assert _ready(healthy, capsys).err == ""
    values = jnp.array([[1.0, -2.0, 1.0], [-3.0, 1.0, -4.0]])
    for _ in range(2):
        result = run(values)
        text = _ready(result, capsys).err
        np.testing.assert_array_equal(result.err_code, [99, 99])
        assert "MPI Process 5 / Chunk index 3" in text
        indices = [" at i =  2  and k =  1", " at i =  1  and k =  2", " at i =  2  and k =  3"]
        assert [text.index(value) for value in indices] == sorted(text.index(value) for value in indices)


@pytest.mark.parametrize("rank", [1, 2, 3])
def test_solver_nan_locations_and_existing_column_status(rank, capsys):
    info = ErrInfo.initialized(2)
    if rank == 1:
        values = jnp.array([1.0, jnp.nan, 2.0])
        run = jax.jit(lambda v: check_nan_soln_1d("test", 2, 3, v, info))
        location, codes = "soln( 2 )", [99, 99]
    elif rank == 2:
        values = jnp.ones((2, 3)).at[1, 2].set(jnp.nan)
        run = jax.jit(lambda v: check_nan_soln_2d("test", 2, 2, 3, v, info))
        location, codes = "soln( 2 , 3 )", [0, 99]
    else:
        values = jnp.ones((2, 3, 2)).at[1, 2, 1].set(jnp.nan)
        run = jax.jit(lambda v: check_nan_soln_3d("test", 2, 2, 3, 2, v, info))
        location, codes = "soln( 2 , 3 , 2 )", [0, 99]
    result = run(values)
    text = _ready(result, capsys).err
    assert location in text and "test NaNs in solution" in text
    np.testing.assert_array_equal(result.err_code, codes)


def test_nonfinite_check_keeps_source_message_and_fatal_status(capsys):
    info = ErrInfo.initialized(2)
    result = jax.jit(lambda v: check_nan(v, "rho", "input", info))(jnp.array([1.0, jnp.inf]))
    text = _ready(result, capsys).err
    assert "rho  is NaN in  input" in text
    np.testing.assert_array_equal(result.err_code, [99, 99])


def test_rcm_clipping_debug_gate_and_preclipping_values(capsys):
    total, cloud = jnp.array([[0.1, 0.2]]), jnp.array([[0.3, 0.1]])
    result = clip_rcm(2, 1, total, "rtm < rcm", cloud)
    text = _ready(result, capsys).err
    assert " at k= 1  at i= 1" in text and "rcm(k) =  0.3" in text
    expected = np.asarray(result)
    error_code.set_debug_level(2)
    clip_rcm.clear_cache()
    result = clip_rcm(2, 1, total, "rtm < rcm", cloud)
    assert _ready(result, capsys).err == ""
    np.testing.assert_array_equal(result, expected)


def test_source_config_and_pdf_printers_use_available_fields(capsys):
    stream = io.StringIO()
    print_clubb_config_flags_api(stream, get_default_config_flags())
    assert stream.getvalue().startswith("iiPDF_type = ")
    assert "l_add_dycore_grid = " in stream.getvalue()
    params = init_pdf_params(3, 2)
    jax.jit(lambda values: print_pdf_params(values, 2))(params)
    jax.effects_barrier()
    captured = capsys.readouterr()
    assert "PDF parameters for column  1" in captured.err
    assert "PDF parameters for column  2" in captured.err
    assert "ice_supersat_frac_2(i,:) = " in captured.err
    assert "..." not in captured.err
