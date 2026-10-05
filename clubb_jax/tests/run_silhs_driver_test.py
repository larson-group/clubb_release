#!/usr/bin/env python3
"""Check SILHS initialization, feedback modes and sampling across real short case runs."""
from pathlib import Path
import argparse
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument('-help', '-h', action='help')
    parser.parse_args()
    from clubb_jax.run_jax import ensure_environment
    ensure_environment()

import jax
import jax.numpy as jnp
import numpy as np
from netCDF4 import Dataset
from utilities.create_case_namelist import create_case_namelist_file
from clubb_jax.src.clubb_case_initalization import init_clubb_case, clean_up_clubb
from clubb_jax.src.advance_clubb_to_end import advance_clubb_to_end


def check_silhs_output_coordinates_initialized_without_sample_fields(mode, stats, tmp_path):
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


def check_real_silhs_driver_modes(case, mode, tmp_path):
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


def check_noninteractive_microphysics_has_no_feedback(case, tmp_path):
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


def check_multicol_sequence_and_all_stats(tmp_path):
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


def main():
    checks = []
    for stats in ('input/stats/multi_col_stats.in', 'input/stats/standard_stats.in', 'none'):
        for mode in ('disabled', 'interactive'):
            checks.append((f'coordinates: {mode}, {stats}',
                           check_silhs_output_coordinates_initialized_without_sample_fields, (mode, stats)))
    for case in ('rico_silhs', 'lba'):
        for mode in ('interactive', 'non-interactive'):
            checks.append((f'driver: {case}, {mode}', check_real_silhs_driver_modes, (case, mode)))
        checks.append((f'feedback: {case}', check_noninteractive_microphysics_has_no_feedback, (case,)))
    checks.append(('two columns, sequence length 3, all stats', check_multicol_sequence_and_all_stats, ()))
    for label, check, arguments in checks:
        print(f'SILHS check: {label}', flush=True)
        try:
            # Each real run gets independent namelist/output storage and compilation caches.
            with tempfile.TemporaryDirectory(prefix='clubb-silhs-') as directory:
                check(*arguments, Path(directory))
        finally:
            jax.clear_caches()
        print(f'SILHS passed: {label}', flush=True)
    print(f'SILHS: {len(checks)} real-case checks passed', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
