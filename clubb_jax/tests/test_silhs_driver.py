"""SILHS-enabled driver checks; native randoms need not match Fortran's stream."""

from pathlib import Path
import jax
import jax.numpy as jnp
import numpy as np
import pytest
from netCDF4 import Dataset
from utilities.create_case_namelist import create_case_namelist_file
from clubb_jax.src.clubb_case_initalization import init_clubb_case, clean_up_clubb
from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats


@pytest.mark.parametrize("mode", ["disabled", "interactive"])
@pytest.mark.parametrize("stats", [
    "input/stats/multi_col_stats.in",
    "input/stats/standard_stats.in",
    "none",
])
def test_silhs_output_coordinates_initialized_without_sample_fields(mode, stats, tmp_path):
    # Reduced stats requests no SILHS fields. Fortran still defines the SILHS
    # coordinates during initialization, even with optional sample output off.
    path = create_case_namelist_file(
        "rico_silhs",
        tmp_path,
        stats=stats,
        debug="0",
        max_iters=1,
        override=f'lh_microphys_type="{mode}"',
    )
    state = init_clubb_case(str(path))
    output_path = state["stats_output_path"]
    zt = np.asarray(state["gr"].zt[0, :]).copy()
    num_samples = state["X_nl_all_levs"].shape[1]
    try:
        assert not state["err_info"].is_fatal()
        if stats == "none":
            assert state["stats_writer"] is None
    finally:
        clean_up_clubb(state)

    if stats == "none":
        assert not Path(output_path).exists()
        return
    with Dataset(output_path) as ds:
        assert ("lh_sample_number" in ds.dimensions) == (mode == "interactive")
        if mode == "interactive":
            assert len(ds.dimensions["lh_sample_number"]) == num_samples
            np.testing.assert_array_equal(ds.variables["lh_zt"][:], zt)
            assert ds.variables["lh_zt"].units == "m"
        elif stats == "input/stats/multi_col_stats.in":
            assert "lh_zt" not in ds.variables
        assert not any(name.startswith(("lh_nl_", "lh_u_")) for name in ds.variables)


@pytest.mark.parametrize("case", ["rico_silhs", "lba"])
@pytest.mark.parametrize("mode", ["interactive", "non-interactive"])
def test_real_silhs_driver_modes(case, mode, tmp_path):
    # Debug 2 exercises the source cloud/category consistency assertions and
    # Kessler diagnostics. Stats are disabled here to isolate real physics.
    path = create_case_namelist_file(
        case,
        tmp_path,
        stats="none",
        debug="2",
        max_iters=3,
        override=f'lh_microphys_type="{mode}"',
    )
    state = init_clubb_case(str(path))
    initial = state["sampling_state"]
    try:
        advance_clubb_to_end(state, l_stdout=False, max_steps=3)
        assert not state["err_info"].is_fatal()
        ns = state["X_nl_all_levs"].shape[1]
        assert ns > 0
        for name in (
            "hydromet",
            "rtm",
            "thlm",
            "lh_rt_clipped",
            "lh_rc_clipped",
            "lh_rv_clipped",
            "lh_Nc_clipped",
            "lh_sample_point_weights",
        ):
            assert np.isfinite(state[name]).all(), name
        np.testing.assert_allclose(
            jnp.sum(state["lh_sample_point_weights"], axis=1), ns, atol=1.0e-12
        )
        np.testing.assert_allclose(
            state["lh_rt_clipped"],
            state["lh_rc_clipped"] + state["lh_rv_clipped"],
            atol=1.0e-18,
        )
        assert np.all(np.asarray(state["lh_rv_clipped"]) >= 0.0)
        assert np.all(np.asarray(state["lh_rc_clipped"]) >= 0.0)
        np.testing.assert_array_equal(initial.one_height_time_matrix, 0.0)
        assert int(state["sampling_state"].prior_iter) == 1
    finally:
        clean_up_clubb(state)
        jax.clear_caches()


def test_subtimestep_statistics_preserve_previous_window_and_count_once():
    previous = JaxStats.empty(l_sample=True, names=("rrm_auto",), ncol=2, max_nlev=3)
    previous = previous.update("rrm_auto", jnp.full((2, 3), 10.0))

    def run(st):
        for value in (1.0, 2.0, 6.0):
            st = st.update("rrm_auto", jnp.full((2, 3), value))
        return st.average_subtimesteps(previous, 3)

    stats = jax.jit(run)(previous)
    np.testing.assert_array_equal(stats.buffers[0], 13.0)
    np.testing.assert_array_equal(stats.nsamples[0], 2)


@pytest.mark.parametrize("case", ["rico_silhs", "lba"])
def test_noninteractive_microphysics_has_no_feedback(case, tmp_path):
    # Two otherwise identical runs must follow the same grid-mean trajectory;
    # sample tendencies are diagnostics when feedback is disabled.
    outputs = []
    for mode in ("disabled", "non-interactive"):
        path = create_case_namelist_file(
            case,
            tmp_path / mode,
            stats="none",
            debug="2",
            max_iters=3,
            override=f'lh_microphys_type="{mode}"',
        )
        state = init_clubb_case(str(path))
        try:
            advance_clubb_to_end(state, l_stdout=False, max_steps=3)
            outputs.append(
                {
                    name: np.asarray(state[name]).copy()
                    for name in ("rtm", "thlm", "hydromet", "rtp2", "thlp2")
                }
            )
        finally:
            clean_up_clubb(state)
    for name in outputs[0]:
        np.testing.assert_array_equal(outputs[0][name], outputs[1][name])
    jax.clear_caches()


def test_multicol_sequence_and_all_stats(tmp_path):
    # No importance sampling permits the source multi-timestep permutation
    # reuse. A real multi-column run also checks lh_zt/lh_sfc output grids.
    path = create_case_namelist_file(
        "rico_silhs",
        tmp_path,
        stats="input/stats/all_stats.in",
        debug="2",
        max_iters=3,
        tout=300,
        override="l_lh_importance_sampling=.false.,lh_sequence_length=3",
    )
    with path.open("a") as handle:
        handle.write("\n&multicol_def\n ngrdcol = 2\n/\n")
    state = init_clubb_case(str(path))
    try:
        advance_clubb_to_end(state, l_stdout=False, max_steps=3)
        assert not state["err_info"].is_fatal()
        assert state["X_nl_all_levs"].shape[0] == 2
        np.testing.assert_array_equal(state["lh_sample_point_weights"], 1.0)
        assert int(state["sampling_state"].prior_iter) == 3
    finally:
        clean_up_clubb(state)
        jax.clear_caches()
