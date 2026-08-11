"""Focused contracts for the dedicated fill-holes diagnostics and viewer."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from dash import dcc
from netCDF4 import Dataset

from dash_app.misc_tab.fill_holes_viewer import SUBTAB
from dash_app.misc_tab.fill_holes_viewer.analysis import (
    below_threshold_counts,
    inspect_dataset,
    load_moment,
)
from dash_app.misc_tab.fill_holes_viewer.layout import build_layout
from dash_app.misc_tab.fill_holes_viewer.figures import profile_figure


def _walk(component):
    yield component
    children = getattr(component, "children", None)
    if children is None:
        return
    if not isinstance(children, (list, tuple)):
        children = [children]
    for child in children:
        if hasattr(child, "to_plotly_json"):
            yield from _walk(child)


def _write_stats(path: Path):
    with Dataset(path, "w") as dataset:
        dataset.createDimension("time", 2)
        dataset.createDimension("zm", 3)
        dataset.createDimension("col", 1)
        dataset.createVariable("time", "f8", ("time",))[:] = [0.0, 60.0]
        dataset.createVariable("zm", "f8", ("zm",))[:] = [10.0, 20.0, 30.0]

        before = dataset.createVariable("wp2_before_hf", "f8", ("time", "zm", "col"))
        before.units = "m^2/s^2"
        before[:] = np.array([[[2.0e-4], [8.0e-4], [1.0e-4]], [[8.0e-4], [1.2e-3], [1.6e-3]]])
        after = dataset.createVariable("wp2_after_hf", "f8", ("time", "zm", "col"))
        after[:] = np.array([[[5.0e-4], [8.0e-4], [2.0e-4]], [[8.0e-4], [1.2e-3], [1.6e-3]]])


def test_analysis_reads_direct_before_after_and_hardcoded_threshold(tmp_path):
    path = tmp_path / "arm_stats.nc"
    _write_stats(path)

    metadata = inspect_dataset(path)
    data = load_moment(path, "wp2")
    before_count, after_count = below_threshold_counts(data, 0)

    assert metadata["moments"] == ["wp2"]
    assert metadata["record_count"] == 2
    assert metadata["column_count"] == 1
    assert data["before"].shape == (2, 3, 1)
    assert data["threshold"].shape == (2, 1)
    assert data["threshold"].tolist() == [[0.0004], [0.0004]]
    assert before_count.tolist() == [2, 0]
    assert after_count.tolist() == [1, 0]


def test_viewer_has_separate_output_folder_and_case_pickers():
    components = list(_walk(build_layout()))
    by_id = {getattr(component, "id", None): component for component in components}

    assert SUBTAB.page_value == "misc-fill-holes-viewer"
    assert isinstance(by_id["hf-output-folder"], dcc.Dropdown)
    assert isinstance(by_id["hf-case-file"], dcc.Dropdown)
    assert isinstance(by_id["hf-field-scale"], dcc.RadioItems)
    assert by_id["hf-field-scale"].value == "log"
    assert "hf-profile-figure" in by_id
    assert "hf-map-figure" in by_id
    assert "hf-summary-figure" in by_id
    assert "hf-threshold-panel" in by_id
    playback_row = next(
        component
        for component in components
        if getattr(component, "className", None) == "hf-playback-row"
    )
    assert any(getattr(child, "id", None) == "hf-speed" for child in playback_row.children[0].children)
    assert getattr(playback_row.children[1], "id", None) == "hf-play"
    layout_text = str(build_layout().to_plotly_json())
    assert "w_tol_sqd" in layout_text
    assert "src/CLUBB_core/constants_clubb.F90" in layout_text
    assert any(
        getattr(component, "className", None) == "hf-right-pane"
        for component in components
    )


def test_profile_field_scale_defaults_to_log_and_can_be_linear(tmp_path):
    path = tmp_path / "arm_stats.nc"
    _write_stats(path)
    data = load_moment(path, "wp2")

    assert profile_figure(data, "wp2", 0, 0).layout.xaxis.type == "log"
    assert profile_figure(data, "wp2", 0, 0).layout.legend.y < 0
    assert profile_figure(data, "wp2", 0, 0).layout.margin.b >= 90
    assert (
        profile_figure(data, "wp2", 0, 0, "linear").layout.xaxis.type
        == "linear"
    )


def test_all_stats_registers_complete_diagnostics_for_clubb_moments():
    text = (Path(__file__).resolve().parents[2] / "input" / "stats" / "all_stats.in").read_text()
    for moment in ("wp2", "rtp2", "thlp2", "up2", "vp2", "rtm", "thlm"):
        for suffix in ("_before_hf", "_after_hf"):
            assert f'"{moment}{suffix}' in text
    assert '"wp2_hf_before' not in text
    assert '"wp2_hf_after' not in text
