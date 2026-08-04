"""Tests for the exact parcel-trajectory replica."""

from pathlib import Path

import numpy as np
import pytest
from dash import Dash, dcc

from dash_app.misc_tab.mixing_length_trajectories import SUBTAB
from dash_app.misc_tab.mixing_length_trajectories.analysis import (
    compute_record,
    load_dataset_record,
    profile_metrics,
)
from dash_app.misc_tab.mixing_length_trajectories.callbacks import (
    register_callbacks,
)
from dash_app.misc_tab.mixing_length_trajectories.figures import (
    parcel_state_figure,
)
from dash_app.misc_tab.mixing_length_trajectories.layout import build_layout

REPO_ROOT = Path(__file__).resolve().parents[2]
VALIDATION_FILE = REPO_ROOT / "output" / "dash_default" / "arm_stats.nc"


@pytest.mark.skipif(
    not VALIDATION_FILE.exists(),
    reason="Local ARM validation output is not present.",
)
def test_replica_matches_fortran_arm_profiles_closely():
    record = load_dataset_record(VALIDATION_FILE, 612, 0)
    result = compute_record(record)

    for calculated, reference in (
        (result.lscale_up, record.fortran_up),
        (result.lscale_down, record.fortran_down),
        (result.lscale, record.fortran_lscale),
    ):
        metrics = profile_metrics(calculated, reference)
        assert metrics["rmse"] < 1.0e-6
        assert metrics["max_abs"] < 1.0e-5


@pytest.mark.skipif(
    not VALIDATION_FILE.exists(),
    reason="Local ARM validation output is not present.",
)
def test_paths_begin_at_interpolated_tke_and_end_nonnegative():
    record = load_dataset_record(VALIDATION_FILE, 612, 0)
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


@pytest.mark.skipif(
    not VALIDATION_FILE.exists(),
    reason="Local ARM validation output is not present.",
)
def test_parcel_thv_figure_is_consistent_with_buoyancy():
    record = load_dataset_record(VALIDATION_FILE, 612, 0)
    result = compute_record(record)
    figure = parcel_state_figure(record, result, "thv")
    path = result.upward_paths[0]
    environment = np.interp(path.altitude, result.z, record.thvm)

    np.testing.assert_allclose(
        np.asarray(figure.data[0].x),
        environment * (1.0 + path.buoyancy / 9.81),
    )
    np.testing.assert_allclose(np.asarray(figure.data[-1].x), record.thvm)


def test_tab_builds_and_registers_callbacks():
    app = Dash(__name__, suppress_callback_exceptions=True)
    register_callbacks(app)
    assert SUBTAB.title == "Mixing Length Trajectories"
    assert SUBTAB.page_value == "misc-mixing-length-trajectories"
    assert len(app.callback_map) == 5


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
