"""Focused parcel-trajectory discovery and presentation contracts."""

import numpy as np
import pytest
import netCDF4 as nc
from dash import Dash, dcc

from dash_app.misc_tab.mixing_length_trajectories.analysis import (
    inspect_dataset,
    profile_metrics,
)
from dash_app.misc_tab.mixing_length_trajectories.callbacks import (
    register_callbacks,
)
from dash_app.misc_tab.mixing_length_trajectories.layout import build_layout


def _empty_stats_file(path):
    from dash_app.misc_tab.mixing_length_trajectories.analysis import REQUIRED_VARIABLES

    with nc.Dataset(path, "w") as dataset:
        dataset.createDimension("time", None)
        for name in REQUIRED_VARIABLES:
            dimensions = ("time",) if name == "time" else ()
            dataset.createVariable(name, "f8", dimensions)
    return path


def _registered_callback(app, name):
    return next(
        entry["callback"].__wrapped__
        for entry in app.callback_map.values()
        if entry["callback"].__name__ == name
    )


def test_empty_dataset_is_reported_before_record_indexing(tmp_path):
    path = _empty_stats_file(tmp_path / "empty_stats.nc")

    with pytest.raises(ValueError, match="contains no time records yet"):
        inspect_dataset(path)


def test_dataset_discovery_can_skip_loading_the_time_axis(tmp_path):
    path = _empty_stats_file(tmp_path / "records_stats.nc")
    with nc.Dataset(path, "a") as dataset:
        dataset["time"][:] = [0.0, 60.0, 120.0]

    metadata = inspect_dataset(path, read_times=False)

    assert metadata["record_count"] == 3
    assert metadata["times"].size == 0


def test_empty_dataset_callbacks_do_not_raise_dash_errors(tmp_path):
    path = _empty_stats_file(tmp_path / "empty_stats.nc")
    app = Dash(__name__, suppress_callback_exceptions=True)
    register_callbacks(app)

    update_mu = _registered_callback(app, "update_mu_control")
    assert update_mu(str(path), 0) == (3.0e-3, {}, 1.0e-3)

    update_diagnostics = _registered_callback(app, "update_diagnostics")
    result = update_diagnostics(str(path), 0, 0, 1.0e-3)
    assert len(result) == 10
    assert result[8].className == "mlt-error"
    assert "contains no time records yet" in result[9]


def test_layout_keeps_figures_live_without_loading_replacement():
    layout = build_layout()
    descendants = []

    def collect(component):
        descendants.append(component)
        children = getattr(component, "children", None)
        if children is None:
            return
        for child in children if isinstance(children, (list, tuple)) else [children]:
            if hasattr(child, "children"):
                collect(child)

    collect(layout)
    assert not any(isinstance(component, dcc.Loading) for component in descendants)
    spaghetti = [
        component
        for component in descendants
        if getattr(component, "className", None) == "mlt-spaghetti-grid"
    ]
    assert len(spaghetti) == 2
    assert [child.id for child in spaghetti[0].children] == [
        "mlt-upward-figure",
        "mlt-downward-figure",
    ]
    assert [child.id for child in spaghetti[1].children] == [
        "mlt-upward-buoyancy-figure",
        "mlt-downward-buoyancy-figure",
    ]
    parcel_states = next(
        component
        for component in descendants
        if getattr(component, "className", None) == "mlt-parcel-state-grid"
    )
    assert [child.id for child in parcel_states.children] == [
        "mlt-parcel-thv-figure",
        "mlt-parcel-rt-figure",
        "mlt-parcel-thl-figure",
    ]


def test_profile_metrics_marks_constant_profile_correlation_undefined():
    with np.testing.assert_no_warnings():
        metrics = profile_metrics(np.ones(4), np.arange(4, dtype=float))

    assert np.isnan(metrics["correlation"])
