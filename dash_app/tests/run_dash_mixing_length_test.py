#!/usr/bin/env python3
"""Check the Dash parcel-trajectory replica against an explicitly supplied Fortran ARM record."""
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument('-help', '-h', action='help')
    parser.add_argument('-stats_file', type=Path, required=True, help='ARM all_stats NetCDF output')
    parser.add_argument('-record', type=int, default=612, help='Zero-based saved record (default: 612)')
    args = parser.parse_args()
    if not args.stats_file.is_file():
        parser.error(f'Stats file not found: {args.stats_file}')
    from utilities.setup_python_venv import ensure_python_venv
    ensure_python_venv('dash')

import numpy as np
from dash_app.misc_tab.mixing_length_trajectories.analysis import (
    load_dataset_record, compute_record, profile_metrics,
)
from dash_app.misc_tab.mixing_length_trajectories.figures import parcel_state_figure

def check_replica_matches_fortran_arm_profiles_closely(stats_file, record_index):
    record = load_dataset_record(stats_file, record_index, 0)
    result = compute_record(record)

    for calculated, reference in (
        (result.lscale_up, record.fortran_up),
        (result.lscale_down, record.fortran_down),
        (result.lscale, record.fortran_lscale),
    ):
        metrics = profile_metrics(calculated, reference)
        assert metrics["rmse"] < 1.0e-6
        assert metrics["max_abs"] < 1.0e-5



def check_paths_begin_at_interpolated_tke_and_end_nonnegative(stats_file, record_index):
    record = load_dataset_record(stats_file, record_index, 0)
    result = compute_record(record)
    for path in (*result.upward_paths, *result.downward_paths):
        assert path.energy[0] == result.tke[path.launch_index]
        assert path.buoyancy[0] == 0.0
        assert path.buoyancy.shape == path.altitude.shape
        assert path.parcel_rt.shape == path.altitude.shape
        assert path.parcel_thl.shape == path.altitude.shape
        assert np.all(path.energy >= -1.0e-13)
        assert np.all(np.isfinite(path.energy))
        assert np.all(np.isfinite(path.buoyancy))
        assert np.all(np.isfinite(path.parcel_rt))
        assert np.all(np.isfinite(path.parcel_thl))



def check_parcel_thv_figure_is_consistent_with_buoyancy(stats_file, record_index):
    record = load_dataset_record(stats_file, record_index, 0)
    result = compute_record(record)
    figure = parcel_state_figure(record, result, "thv")
    path = result.upward_paths[0]
    environment = np.interp(path.altitude, result.z, record.thvm)

    np.testing.assert_allclose(
        np.asarray(figure.data[0].x),
        environment * (1.0 + path.buoyancy / 9.81),
    )
    np.testing.assert_allclose(np.asarray(figure.data[-1].x), record.thvm)


if __name__ == '__main__':
    check_replica_matches_fortran_arm_profiles_closely(args.stats_file, args.record)
    check_paths_begin_at_interpolated_tke_and_end_nonnegative(args.stats_file, args.record)
    check_parcel_thv_figure_is_consistent_with_buoyancy(args.stats_file, args.record)
    print(f'Dash ARM parcel-trajectory checks passed at record {args.record}')
