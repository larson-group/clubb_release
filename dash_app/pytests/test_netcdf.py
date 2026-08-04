from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset

from dash_app.plot_tab import benchmark_overlay
from dash_app.plot_tab.plot_types import shared as plot_shared

from dash_app.shared.netcdf import (
    coordinate_values,
    file_metadata_signature,
    file_signature,
    find_dimension,
    time_units_factor,
    time_values_to_seconds,
)


def test_file_signatures_cover_cache_and_serializable_forms(tmp_path):
    path = tmp_path / "sample.nc"
    path.write_bytes(b"netcdf-ish")

    cache_signature = file_signature(path)
    metadata_signature = file_metadata_signature(path)

    assert cache_signature == (path, path.stat().st_mtime_ns, path.stat().st_size)
    assert metadata_signature == {
        "path": path,
        "size_bytes": path.stat().st_size,
        "mtime_ns": path.stat().st_mtime_ns,
    }
    missing = tmp_path / "missing.nc"
    assert file_signature(missing) == (missing, None, None)


def test_find_dimension_is_case_insensitive_and_preserves_dataset_name():
    assert find_dimension(("column", "Time", "ZM"), {"t", "time"}) == "Time"
    assert find_dimension(("column", "Time", "ZM"), {"z", "zt"}) is None


def test_coordinate_values_reads_variable_or_uses_index_fallback(tmp_path):
    path = tmp_path / "coordinates.nc"
    with Dataset(path, "w") as dataset:
        dataset.createDimension("zm", 3)
        coordinate = dataset.createVariable("zm", "f8", ("zm",))
        coordinate[:] = [50.0, 100.0, 150.0]
        assert np.array_equal(
            coordinate_values(dataset, "zm", 3),
            np.array([50.0, 100.0, 150.0]),
        )
        assert np.array_equal(
            coordinate_values(dataset, "missing", 3), np.array([0.0, 1.0, 2.0])
        )


def test_time_conversions_preserve_plotting_conventions():
    assert time_units_factor("seconds") == 1.0
    assert time_units_factor("minutes since 2000-01-01") == 60.0
    assert time_units_factor("hours") == 3600.0
    assert time_units_factor("fortnights") == 1.0
    assert np.array_equal(
        time_values_to_seconds([0.0, 2.0], "minutes"),
        np.array([0.0, 120.0]),
    )
    assert np.array_equal(
        time_values_to_seconds(
            [0.0, 0.5], "hours since 2017-07-01 01:30:00"
        ),
        np.array([5400.0, 7200.0]),
    )


def test_existing_modules_delegate_to_shared_netcdf_helpers():
    assert benchmark_overlay._file_signature is file_signature
    assert benchmark_overlay._find_dim is find_dimension
    assert benchmark_overlay._time_values_seconds is time_values_to_seconds
    assert plot_shared._file_signature is file_signature
    assert plot_shared.find_dim is find_dimension
    assert plot_shared.time_units_factor is time_units_factor


def test_shared_netcdf_consumers_import_as_dashboard_packages():
    repo_root = Path(__file__).resolve().parents[2]
    command = (
        "from dash_app.plot_tab import benchmark_overlay; "
        "from dash_app.plot_tab.plot_types import shared; "
    )

    result = subprocess.run(
        [sys.executable, "-c", command],
        cwd=repo_root,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
