import csv
from pathlib import Path

import pytest
from dash import Dash, html
from dash.development.base_component import Component

from dash_app.profile_tab.callbacks import (
    load_profile_plot_data,
    profile_rename_available,
    profile_selection_preferences,
    profile_control_rows,
    reconcile_profile_selection,
)
from dash_app.profile_tab.figures import (
    METRICS,
    decomposition_figure,
    profile_figure,
    relative_figure,
    timer_options,
    variability_figure,
)
from dash_app.profile_tab.library import (
    delete_profile_library_entry,
    enrich_summary_rows,
    profile_option,
)
from dash_app.profile_tab.layout import active_profile_items, available_profile_buttons, build_layout
from dash_app.profile_tab.runtime import (
    normalize_profile_settings,
    overwrite_confirmation,
    profile_command,
    profile_results_complete,
    profile_save_target,
    read_profile_data,
    read_profile_results,
)
from dash_app.profile_tab.tab import build_tab
from dash_app.shared import actions
from utilities.timing_profiles import BATCH_FIELDS, GROUP_WALL_TIMER, TIMING_FIELDS, write_profile_manifest


def settings(tmp_path: Path, executable: str = "") -> dict[str, object]:
    return {
        "case_name": "arm",
        "processes": "1,4",
        "columns": "2,8",
        "warmups": 1,
        "repetitions": 3,
        "output": str(tmp_path / "timing"),
        "name": "CPU baseline",
        "timeout": 30,
        "continue_on_error": True,
        "config": "default",
        "override": "C2=2.0",
        "executable": executable,
        "extra_args": "-max_iters 20 -dt_main 4",
    }


def _profile_record(run_id: str, label: str, case_name: str, started: str) -> dict[str, str]:
    return {
        "run_id": run_id,
        "label": label,
        "case_name": case_name,
        "started_utc": started,
        "git_revision": "1234567890",
        "status": "complete",
    }


def test_profile_picker_shows_three_recent_unselected_profiles_and_simple_labels():
    records = [
        _profile_record("newest", "GPU baseline", "arm", "2026-08-12T12:00:00Z"),
        _profile_record("middle", "CPU baseline", "bomex", "2026-08-11T12:00:00Z"),
        _profile_record("older", "Scaling", "rico", "2026-08-10T12:00:00Z"),
        _profile_record("oldest", "Paper run", "atex", "2026-08-09T12:00:00Z"),
    ]

    buttons = available_profile_buttons(records, [], expanded=False)

    assert profile_option(records[0]) == {"label": "GPU baseline · arm", "value": "newest"}
    assert [button.children for button in buttons] == [
        "GPU baseline · arm",
        "CPU baseline · bomex",
        "Scaling · rico",
    ]
    assert buttons[0].id == {"type": "profile-add-run", "run_id": "newest"}
    assert "2026-08-12 12:00" in buttons[0].title


def test_profile_picker_promotes_next_profile_and_renders_removable_active_pills():
    records = [
        _profile_record("newest", "GPU baseline", "arm", "2026-08-12T12:00:00Z"),
        _profile_record("older", "CPU baseline", "bomex", "2026-08-11T12:00:00Z"),
    ]

    buttons = available_profile_buttons(records, ["newest"], expanded=False)
    active = active_profile_items(records, ["newest"])[0]

    assert [button.id["run_id"] for button in buttons] == ["older"]
    assert active.children[0].children == "GPU baseline · arm"
    assert active.children[1].id == {"type": "profile-remove-run", "run_id": "newest"}
    assert active.children[1].children == "×"


def test_expanded_profile_picker_adds_guarded_delete_buttons():
    records = [
        _profile_record("newest", "GPU baseline", "arm", "2026-08-12T12:00:00Z"),
        _profile_record("older", "CPU baseline", "bomex", "2026-08-11T12:00:00Z"),
    ]

    available = available_profile_buttons(
        records,
        ["newest"],
        expanded=True,
        delete_confirmation={"run_id": "older"},
    )[0]
    active = active_profile_items(
        records,
        ["newest"],
        expanded=True,
    )[0]

    assert available.children[0].id == {"type": "profile-add-run", "run_id": "older"}
    assert available.children[1].id == {"type": "profile-delete-run", "run_id": "older"}
    assert available.children[1].children == "Confirm"
    assert active.children[1].id == {"type": "profile-delete-run", "run_id": "newest"}
    assert active.children[2].id == {"type": "profile-remove-run", "run_id": "newest"}


def test_profile_library_deletion_is_limited_to_recognized_children(tmp_path):
    library = tmp_path / "timing"
    profile = library / "cpu_baseline"
    profile.mkdir(parents=True)
    write_profile_manifest(
        profile,
        {"run_id": "cpu_baseline", "label": "CPU baseline", "case_name": "arm"},
    )

    assert delete_profile_library_entry(library, "cpu_baseline") == profile
    assert not profile.exists()

    write_profile_manifest(
        library,
        {"run_id": "timing", "label": "Library root", "case_name": "arm"},
    )
    with pytest.raises(ValueError, match="inside the selected library"):
        delete_profile_library_entry(library, "timing")
    assert library.exists()


def component_ids(component) -> set[str]:
    if not isinstance(component, Component):
        return set()
    found = {component.id} if isinstance(getattr(component, "id", None), str) else set()
    children = getattr(component, "children", None)
    if children is None:
        return found
    for child in children if isinstance(children, (list, tuple)) else [children]:
        found.update(component_ids(child))
    return found


def find_component(component, component_id: str):
    if not isinstance(component, Component):
        return None
    if getattr(component, "id", None) == component_id:
        return component
    children = getattr(component, "children", None)
    for child in children if isinstance(children, (list, tuple)) else [children]:
        found = find_component(child, component_id)
        if found is not None:
            return found
    return None


def test_profile_command_covers_timer_and_forwarded_run_settings(tmp_path):
    executable = tmp_path / "clubb_standalone"
    executable.write_text("fake", encoding="utf-8")

    command = profile_command(settings(tmp_path, str(executable)))

    assert command[2].endswith("utilities/time_clubb.py")
    assert command[3] == "arm"
    assert command[command.index("-processes") + 1] == "1,4"
    assert command[command.index("-batch_sizes") + 1] == "2,8"
    assert command[command.index("-config") + 1] == "default"
    assert command[command.index("-override") + 1] == "C2=2.0"
    assert command[command.index("-exe") + 1] == str(executable.resolve())
    assert command[-4:] == ["-max_iters", "20", "-dt_main", "4"]

    overwrite_settings = settings(tmp_path, str(executable))
    overwrite_settings["overwrite"] = True
    overwrite_command = profile_command(overwrite_settings)
    assert "-overwrite" in overwrite_command
    assert normalize_profile_settings(overwrite_settings)["run_id"] == "CPU_baseline"


@pytest.mark.parametrize(("implementation", "flag"), (("python", "-python"), ("jax", "-jax")))
def test_profile_command_selects_implementation_and_install(tmp_path, implementation, flag):
    install = tmp_path / "install"
    install.mkdir()
    request = settings(tmp_path)
    request.update({"implementation": implementation, "install_dir": str(install)})

    command = profile_command(request)

    assert flag in command
    assert command[command.index("-install_dir") + 1] == str(install.resolve())


def test_existing_profile_requires_overwrite_confirmation(tmp_path):
    request = settings(tmp_path)
    target = profile_save_target(request)
    profile = Path(target["path"])
    profile.mkdir(parents=True)
    write_profile_manifest(
        profile,
        {"run_id": target["run_id"], "label": request["name"], "case_name": "arm"},
    )

    refreshed = profile_save_target(request)
    assert refreshed["exists"] and refreshed["replaceable"]
    assert "Replace the saved profile" in overwrite_confirmation(request)

    conflict_request = {**request, "name": "ordinary folder"}
    Path(profile_save_target(conflict_request)["path"]).mkdir()
    with pytest.raises(ValueError, match="unrecognized directory"):
        overwrite_confirmation(conflict_request)


def test_profile_rename_requires_a_changed_nonempty_label():
    pending = {"name": "CPU baseline"}

    assert profile_rename_available("CPU baseline", pending) is False
    assert profile_rename_available("  CPU baseline  ", pending) is False
    assert profile_rename_available("", pending) is False
    assert profile_rename_available("CPU baseline 2", pending) is True


def test_profile_settings_reject_invalid_sweeps(tmp_path):
    request = settings(tmp_path)
    request["processes"] = "1,0"

    with pytest.raises(ValueError, match="processes"):
        normalize_profile_settings(request)

    request = settings(tmp_path)
    request["extra_args"] = "-batch_size 4"
    with pytest.raises(ValueError, match="managed by the Profile benchmark"):
        normalize_profile_settings(request)


def test_profile_results_select_only_rows_created_after_launch(tmp_path):
    output = tmp_path / "timing"
    output.mkdir()
    for index, run_id in enumerate(("old", "new")):
        profile = output / run_id
        profile.mkdir()
        write_profile_manifest(
            profile,
            {
                "run_id": run_id,
                "label": run_id,
                "case_name": "arm",
                "started_utc": f"2026-08-12T12:00:0{index}+00:00",
                "status": "complete",
                "benchmark": {"processes": [1], "batch_sizes": [1]},
                "results": {"timer_backends": ["cpu_time"], "time_bases": ["cpu_time"]},
            },
        )
        with (profile / "batches.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=BATCH_FIELDS)
            writer.writeheader()
            writer.writerow(
                {
                    "batch_id": "p0001_b000001",
                    "process_count": 1,
                    "batch_size": 1,
                    "total_batch_size": 1,
                    "status": "complete",
                    "warmups_completed": 0,
                    "measurements_completed": 1,
                    "failed_runs": 0,
                }
            )
        with (profile / "timings.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=TIMING_FIELDS)
            writer.writeheader()
            writer.writerow(
                {
                    "batch_id": "p0001_b000001",
                    "phase": "measured",
                    "repetition": 1,
                    "process": 0,
                    "thread": 0,
                    "timer": "advance_clubb_to_end",
                    "calls": 1,
                    "inclusive_s": 1,
                    "exclusive_s": 1,
                    "status": "success",
                    "return_code": 0,
                    "message": "",
                }
            )

    run_id, rows = read_profile_results(
        {
            "output": str(output),
            "existing_run_ids": ["old"],
        }
    )

    assert run_id == "new"
    assert {row["run_id"] for row in rows} == {"new"}

    overwritten_id, overwritten_rows = read_profile_results(
        {"output": str(output), "run_id": "old", "existing_run_ids": ["old", "new"]}
    )
    assert overwritten_id == "old"
    assert {row["run_id"] for row in overwritten_rows} == {"old"}

    data_id, summary_rows, process_rows = read_profile_data(
        {"output": str(output), "run_id": "new"}
    )
    assert data_id == "new"
    assert {row["run_id"] for row in summary_rows} == {"new"}
    assert {row["run_id"] for row in process_rows} == {"new"}


def test_live_plot_data_replaces_selected_disk_snapshot(monkeypatch):
    loaded = {}
    stored_summary = [
        {"profile_id": "baseline", "timer_max_seconds": 3.0},
        {"profile_id": "active", "timer_max_seconds": 99.0},
    ]
    stored_processes = [
        {"profile_id": "baseline", "inclusive_seconds": "3.0"},
        {"profile_id": "active", "inclusive_seconds": "99.0"},
    ]
    live_summary = [{"profile_id": "active", "timer_max_seconds": 1.0}]
    live_processes = [{"profile_id": "active", "inclusive_seconds": "1.0"}]
    monkeypatch.setattr(
        "dash_app.profile_tab.callbacks.load_profile_library_data",
        lambda _output, selected: (
            loaded.update(selected=list(selected)) or stored_summary,
            stored_processes,
        ),
    )
    monkeypatch.setattr(
        "dash_app.profile_tab.callbacks.read_profile_data",
        lambda _job: ("active", live_summary, live_processes),
    )

    summaries, processes = load_profile_plot_data(
        "/tmp/profiles",
        ["baseline", "active"],
        {"state": "running", "run_id": "active"},
    )

    assert [(row["profile_id"], row["timer_max_seconds"]) for row in summaries] == [
        ("baseline", 3.0),
        ("active", 1.0),
    ]
    assert [(row["profile_id"], row["inclusive_seconds"]) for row in processes] == [
        ("baseline", "3.0"),
        ("active", "1.0"),
    ]
    assert loaded["selected"] == ["baseline"]


def test_live_browser_store_keeps_only_timer_control_metadata():
    summary_rows = [
        {"profile_id": "active", "timer_name": "advance_clubb_to_end", "timer_max_seconds": 1.0},
        {"profile_id": "active", "timer_name": "advance_clubb_to_end", "timer_max_seconds": 1.1},
        {"profile_id": "active", "timer_name": "clean_up_clubb", "timer_max_seconds": 0.1},
    ]
    process_rows = [
        {
            "profile_id": "active",
            "process_count": 12,
            "timer_name": "advance_clubb_to_end",
            "inclusive_seconds": "1.0",
        },
        {
            "profile_id": "active",
            "process_count": 12,
            "timer_name": "clean_up_clubb",
            "inclusive_seconds": "0.1",
        },
    ]

    timers = profile_control_rows(summary_rows)

    assert timers == [
        {"timer_name": "advance_clubb_to_end"},
        {"timer_name": "clean_up_clubb"},
    ]


def test_profile_results_require_every_primary_measurement():
    job = {
        "settings": {
            "processes": [1, 2],
            "columns": [1, 4],
            "repetitions": 2,
        }
    }
    rows = [
        {"timer_name": "advance_clubb_to_end", "status": "success"}
        for _ in range(8)
    ]

    assert profile_results_complete(job, rows)
    assert not profile_results_complete(job, rows[:-1])
    assert not profile_results_complete(job, [{**rows[0], "status": "failed"}, *rows[1:]])


def test_profile_figure_averages_repetitions_by_process_count():
    rows = [
        {
            "timer_name": "advance_clubb_to_end",
            "status": "success",
            "process_count": "2",
            "columns_per_process": "4",
            "total_columns": "8",
            "timer_max_seconds": value,
            "throughput_columns_per_second": "10",
        }
        for value in ("1", "3")
    ]

    figure = profile_figure(rows, "advance_clubb_to_end", "timer_max_seconds")

    assert len(figure.data) == 1
    assert list(figure.data[0].x) == [8]
    assert list(figure.data[0].y) == [2.0]
    assert list(figure.data[0].error_y.array) == [1.0]
    assert timer_options(rows, "advance_clubb_to_end")[0]["value"] == "advance_clubb_to_end"


def test_profile_comparison_figures_keep_runs_distinct_and_legends_below():
    rows = []
    for run_id, label, seconds in (("base", "Baseline", 2.0), ("new", "New build", 2.5)):
        rows.append(
            {
                "profile_id": run_id,
                "profile_label": label,
                "timer_name": "advance_clubb_to_end",
                "status": "success",
                "process_count": "2",
                "columns_per_process": "4",
                "total_columns": "8",
                "timer_max_seconds": seconds,
            }
        )

    comparison = profile_figure(rows, "advance_clubb_to_end", "timer_max_seconds")
    relative = relative_figure(rows, "advance_clubb_to_end", "timer_max_seconds", "base")

    assert len(comparison.data) == 2
    assert relative.data[0].y[0] == pytest.approx(25.0)
    assert comparison.layout.legend.y < 0
    assert relative.layout.legend.y < 0
    for trace in (*comparison.data, *relative.data):
        assert trace.hovertemplate.startswith("%{fullData.name}<br>")
        assert trace.hovertemplate.endswith("<extra></extra>")


def test_profile_decomposition_uses_one_critical_process_and_exclusive_costs():
    rows = []
    for process_index, primary_time, core_cost in ((0, 1.0, 0.6), (1, 1.2, 0.8)):
        for timer, inclusive, exclusive in (
            ("advance_clubb_to_end", primary_time, 0.2),
            ("advance_clubb_core_api", core_cost, core_cost),
        ):
            rows.append(
                {
                    "profile_id": "run",
                    "profile_label": "Run",
                    "process_count": "2",
                    "columns_per_process": "4",
                    "total_columns": "8",
                    "repetition": "1",
                    "process_index": str(process_index),
                    "thread": "0",
                    "status": "success",
                    "timer_name": timer,
                    "inclusive_seconds": inclusive,
                    "exclusive_seconds": exclusive,
                }
            )

    figure = decomposition_figure(rows, "run", "advance_clubb_to_end", 2, top_n=8)
    spread = variability_figure(rows, "run", "advance_clubb_to_end", 2)

    values_by_timer = {trace.name: trace.y[0] for trace in figure.data}
    assert values_by_timer["advance_clubb_core_api"] == pytest.approx(0.8)
    assert values_by_timer["advance_clubb_to_end"] == pytest.approx(0.2)
    assert sorted(spread.data[0].y) == pytest.approx([1.0, 1.2])
    assert figure.layout.legend.y < 0
    assert spread.layout.legend.y < 0
    assert all(trace.hovertemplate.startswith("%{fullData.name}<br>") for trace in figure.data)


def test_profile_derived_metrics_and_selection():
    rows = [
        {
            "run_id": "run",
            "profile_id": "run",
            "process_count": "2",
            "columns_per_process": "4",
            "total_columns": "8",
            "repetition": "1",
            "timer_name": "advance_clubb_core_api",
            "calls": "5",
            "timer_max_seconds": "2",
            "timer_mean_seconds": "1.6",
        },
        {
            "run_id": "run",
            "profile_id": "run",
            "process_count": "2",
            "columns_per_process": "4",
            "total_columns": "8",
            "repetition": "1",
            "timer_name": "advance_clubb_to_end",
            "calls": "1",
            "timer_max_seconds": "4",
            "timer_mean_seconds": "3.2",
            "vertical_levels": "100",
        },
    ]
    enriched = enrich_summary_rows(rows)
    primary = enriched[1]
    assert primary["throughput_column_steps_per_second"] == pytest.approx(10.0)
    assert primary["throughput_grid_box_iterations_per_second"] == pytest.approx(1000.0)
    assert primary["process_imbalance_ratio"] == pytest.approx(1.25)

    catalog = [
        {"run_id": "old", "case_name": "arm", "git_revision": "a", "omp_num_threads": "1", "backends": "cpu_time", "time_bases": "cpu_time", "forwarded_args": "-max_iters 5"},
        {"run_id": "new", "case_name": "bomex", "git_revision": "b", "omp_num_threads": "1", "backends": "gptl", "time_bases": "wallclock", "forwarded_args": "-max_iters 10"},
    ]
    selected, baseline, detail = reconcile_profile_selection(catalog, ["missing"], None, None, ["new"])
    assert selected == ["new"]
    assert baseline == detail == "new"


def test_profile_metric_labels_and_defaults_are_the_streamlined_throughput_set(monkeypatch):
    monkeypatch.setattr(
        "dash_app.profile_tab.layout.build_selector_trigger",
        lambda component_id: html.Button(id=component_id),
    )
    layout = build_layout(
        {
            "cases": ["arm"],
            "default_case": "arm",
            "configs": [],
            "default_config": "default",
            "executables": [],
        }
    )

    assert METRICS == {
        "throughput_columns_per_second": "Throughput (runs / second)",
        "throughput_column_steps_per_second": "Throughput (iterations / second)",
        "throughput_grid_box_iterations_per_second": "Throughput (grid box iterations / second)",
        "timer_max_seconds": "Runtime (seconds)",
        "process_imbalance_ratio": "Process imbalance (max / mean)",
    }
    assert find_component(layout, "profile-result-metric").value == "throughput_columns_per_second"
    assert find_component(layout, "profile-detail-process") is None
    assert find_component(layout, "profile-process-results") is None
    assert timer_options(
        [{"timer_name": GROUP_WALL_TIMER}], GROUP_WALL_TIMER
    ) == [{"label": "Process-group wall time", "value": GROUP_WALL_TIMER}]


def test_replacement_profile_supersedes_selected_name_versions():
    catalog = [
        _profile_record("baseline", "CPU baseline", "arm", "2026-08-13T12:00:00Z"),
        _profile_record("baseline_old", "CPU baseline", "arm", "2026-08-12T12:00:00Z"),
        _profile_record("baseline_older", "CPU baseline", "arm", "2026-08-11T12:00:00Z"),
        _profile_record("gpu", "GPU baseline", "arm", "2026-08-10T12:00:00Z"),
    ]

    selected, baseline, detail = reconcile_profile_selection(
        catalog,
        ["baseline_old", "baseline_older", "gpu"],
        "baseline_old",
        "baseline_older",
        ["baseline"],
        "baseline",
    )

    assert selected == ["gpu", "baseline"]
    assert baseline == "baseline"
    assert detail == "baseline"

    legacy_selected, legacy_baseline, legacy_detail = reconcile_profile_selection(
        catalog,
        ["baseline_older", "baseline_old", "baseline", "gpu"],
        "baseline_older",
        "baseline_old",
    )
    assert legacy_selected == ["baseline", "gpu"]
    assert legacy_baseline == "baseline"
    assert legacy_detail == "baseline"


def test_empty_profile_selection_is_only_defaulted_for_a_new_library():
    catalog = [
        _profile_record("newest", "GPU", "arm", "2026-08-13T12:00:00Z"),
        _profile_record("older", "CPU", "arm", "2026-08-12T12:00:00Z"),
    ]

    selected, baseline, detail = reconcile_profile_selection(
        catalog, [], None, None, select_default=True
    )
    assert selected == ["newest"]
    assert baseline == detail == "newest"

    selected, baseline, detail = reconcile_profile_selection(catalog, [], None, None)
    assert selected == []
    assert baseline is None
    assert detail is None


def test_only_active_jobs_and_new_import_events_are_preferred():
    action = {"kind": "imported", "run_ids": ["imported"]}
    preferred, replacement = profile_selection_preferences(
        {"state": "finished", "run_id": "rtx3080"},
        action,
        library_action_triggered=False,
    )
    assert preferred == []
    assert replacement is None

    preferred, replacement = profile_selection_preferences(
        {"state": "running", "run_id": "rtx3080"},
        action,
        library_action_triggered=True,
    )
    assert preferred == ["imported", "rtx3080"]
    assert replacement == "rtx3080"


def test_profile_layout_contains_settings_panel_graph_and_lifecycle_controls():
    initial = {
        "cases": ["arm", "bomex"],
        "default_case": "arm",
        "configs": [{"label": "default", "value": "default", "path": "/tmp/default"}],
        "default_config": "default",
        "executables": [{"label": "Auto", "value": ""}],
    }
    ids = component_ids(build_layout(initial))

    assert {
        "profile-graph",
        "profile-relative-graph",
        "profile-decomposition-graph",
        "profile-variability-graph",
        "profile-selected-runs",
        "profile-picker-toggle",
        "profile-picker-expanded",
        "profile-delete-confirm",
        "profile-delete-expiry",
        "profile-available-list",
        "profile-active-list",
        "profile-library-import",
        "profile-library-export",
        "profile-case",
        "profile-processes",
        "profile-columns",
        "profile-warmups",
        "profile-repeats",
        "profile-output",
        "profile-config",
        "profile-override",
        "profile-executable",
        "profile-extra-args",
        "profile-selected-build-badge",
        "profile-name-control",
        "profile-name-status",
        "profile-selection-library",
        "profile-pending-run",
        "profile-overwrite-modal",
        "profile-overwrite-name",
        "profile-overwrite-message",
        "profile-overwrite-button",
        "profile-rename-button",
        "profile-overwrite-cancel-button",
        "profile-start",
        "profile-cancel",
    } <= ids
    assert "profile-heatmap-graph" not in ids
    assert "profile-visible-charts" not in ids
    name_input = find_component(build_layout(initial), "profile-name")
    assert name_input.debounce is False


def test_profile_layout_separates_benchmark_output_and_view_controls_with_help():
    initial = {
        "cases": ["arm"],
        "default_case": "arm",
        "configs": [{"label": "default", "value": "default"}],
        "default_config": "default",
        "executables": [{"label": "Auto", "value": ""}],
    }
    layout = build_layout(initial)
    serialized = str(layout.to_plotly_json())

    assert "profile-benchmark-panel" in serialized
    assert serialized.count("profile-benchmark-group") >= 4
    assert "profile-benchmark-identity" in serialized
    assert "profile-benchmark-workload" in serialized
    assert "profile-benchmark-sampling" in serialized
    assert "profile-benchmark-launch" in serialized
    assert "Choose compiler &amp; run" in serialized or "Choose compiler & run" in serialized
    assert serialized.index("profile-selected-build-badge") < serialized.index("profile-start")
    assert serialized.index("profile-start") < serialized.index("profile-cancel")
    assert serialized.index("profile-cancel") < serialized.index("profile-status")
    assert serialized.count("profile-sidebar-section") == 2
    assert "profile-output-section" in serialized
    assert "profile-view-section" in serialized
    assert "profile-chart-controls" in serialized
    assert serialized.index("profile-start") < serialized.index("profile-graph")
    assert serialized.index("profile-output") > serialized.index("profile-graph")
    assert serialized.index("profile-output") < serialized.index("profile-selected-runs")
    assert serialized.index("profile-selected-runs") < serialized.index("profile-result-timer")
    assert "profile-results-state" in serialized
    assert "Create a timing profile" not in serialized
    assert "Outputs &amp; view" not in serialized and "Outputs & view" not in serialized
    assert "Compare profiles" not in serialized
    assert "profile-chart-subtitle" not in serialized
    assert "profile-comparison-warnings" not in serialized
    assert serialized.count("profile-help-icon") >= 20
    assert "title" in serialized

    stylesheet = (Path(__file__).parents[1] / "assets" / "15_tab_profile_theme.css").read_text(encoding="utf-8")
    plot_grid = stylesheet.split(".profile-plot-grid {", 1)[1].split("}", 1)[0]
    assert "grid-template-columns: repeat(2, minmax(0, 1fr));" in plot_grid
    workspace = stylesheet.split(".profile-workspace {", 1)[1].split("}", 1)[0]
    results_panel = stylesheet.split(".profile-results-panel {", 1)[1].split("}", 1)[0]
    assert "overflow: visible;" in workspace
    assert "overflow: visible;" in results_panel
    assert ".profile-controls-panel {\n  padding: 13px;\n  overflow: visible;" in stylesheet


def test_profile_tab_registers_callbacks(monkeypatch):
    monkeypatch.setattr(
        "dash_app.profile_tab.tab.discover_profile_state",
        lambda: {
            "cases": ["arm"],
            "default_case": "arm",
            "configs": [{"label": "default", "value": "default"}],
            "default_config": "default",
            "executables": [{"label": "Auto", "value": ""}],
        },
    )
    app = Dash(__name__)

    tab = build_tab(app)

    assert tab.label == "Profile"
    assert tab.value == "profile"
    assert any("profile-action.data" in key for key in app.callback_map)
    assert any("profile-overwrite-modal.className" in key for key in app.callback_map)
    assert "profile-command-preview.children" in app.callback_map
    assert any("profile-name.value" in key for key in app.callback_map)
    assert not any("profile-name.className" in key for key in app.callback_map)
    assert "profile-graph.figure" in next(
        key for key in app.callback_map if "profile-graph.figure" in key
    )
    data_callback = next(
        entry
        for key, entry in app.callback_map.items()
        if key == "profile-results.data"
    )
    input_ids = {item["id"] for item in data_callback["inputs"]}
    assert "profile-interval" in input_ids
    assert "profile-active-results" not in input_ids
    graph_callback = next(
        (key, entry)
        for key, entry in app.callback_map.items()
        if "profile-graph.figure" in key
    )
    graph_key, graph_callback = graph_callback
    assert "profile-status.children" not in graph_key
    graph_input_ids = {item["id"] for item in graph_callback["inputs"]}
    assert "profile-job" in graph_input_ids
    assert "profile-selected-runs" in graph_input_ids
    assert "profile-interval" not in graph_input_ids
    assert "profile-results" not in graph_input_ids
    assert "profile-process-results" not in graph_input_ids
    lifecycle_callback = next(
        entry
        for key, entry in app.callback_map.items()
        if "profile-job.data" in key and "profile-status.children" in key
    )
    lifecycle_input_ids = {item["id"] for item in lifecycle_callback["inputs"]}
    assert "dashboard-broker-jobs" in lifecycle_input_ids
    assert "profile-selected-runs" not in lifecycle_input_ids


def test_profile_launch_is_owned_and_watched_by_broker(monkeypatch):
    captured = {}
    job = {
        "state": "running",
        "pid": 42,
        "output": "/tmp/profile",
        "command_display": "time_clubb.py arm",
    }
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"profile": None})
    monkeypatch.setattr(actions, "start_profile_process", lambda _settings: dict(job))
    monkeypatch.setattr(
        actions,
        "set_broker_job",
        lambda kind, record: captured.update(kind=kind, record=record),
    )
    monkeypatch.setattr(actions, "publish_event", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        actions,
        "_background",
        lambda target, record: captured.update(target=target, watched=record),
    )

    result = actions.launch_profile_request({"case_name": "arm"})

    assert result["status"] == "started"
    assert captured["kind"] == "profile"
    assert captured["record"]["broker_managed"] is True
    assert captured["target"] is actions._watch_profile
    assert captured["watched"]["pid"] == 42


def test_profile_stop_uses_graceful_timer_shutdown(monkeypatch):
    captured = {}
    job = {"state": "running", "pid": 42, "log": "/tmp/profile.log"}
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"profile": job})
    monkeypatch.setattr(actions, "stop_profile_process", lambda record: captured.update(job=record))
    monkeypatch.setattr(
        actions,
        "update_broker_job",
        lambda kind, **updates: captured.update(kind=kind, updates=updates),
    )
    monkeypatch.setattr(actions, "publish_event", lambda *_args, **_kwargs: None)

    result = actions.stop_profile()

    assert result == {"status": "stop_requested", "pid": 42}
    assert captured["job"] == job
    assert captured["kind"] == "profile"
    assert captured["updates"] == {"state": "stopping"}


def test_profile_watcher_publishes_terminal_results(monkeypatch):
    updates = []
    rows = [{"timer_name": "advance_clubb_to_end", "status": "success"}]
    job = {
        "pid": 42,
        "output": "/tmp/profile",
        "log": "/tmp/profile.log",
        "settings": {
            "processes": [1],
            "columns": [1],
            "repetitions": 1,
        },
    }
    monkeypatch.setattr(actions, "read_profile_results", lambda _job: ("run-1", rows))
    monkeypatch.setattr(actions, "profile_process_status", lambda _job: (False, 0))
    monkeypatch.setattr(actions, "read_profile_log_tail", lambda _path: "complete")
    monkeypatch.setattr(actions, "broker_jobs", lambda: {"profile": {"state": "running"}})
    monkeypatch.setattr(
        actions,
        "update_broker_job",
        lambda kind, **values: updates.append((kind, values)),
    )
    monkeypatch.setattr(actions, "release_profile_process", lambda _pid: None)
    monkeypatch.setattr(actions, "publish_event", lambda *_args, **_kwargs: None)

    actions._watch_profile(job)

    assert updates[-1][0] == "profile"
    assert updates[-1][1]["state"] == "finished"
    assert updates[-1][1]["returncode"] == 0
    assert updates[-1][1]["run_id"] == "run-1"


def test_manager_shutdown_includes_active_profile_job(monkeypatch):
    stopped = []

    class Store:
        def list_kind(self, _kind):
            return []

    monkeypatch.setattr(
        actions,
        "broker_jobs",
        lambda: {
            "compile": None,
            "profile": {"state": "running", "pid": 42},
            "runs": {},
            "tune": None,
            "loss_runs": {},
        },
    )
    monkeypatch.setattr(actions, "stop_profile", lambda: stopped.append("profile"))
    monkeypatch.setattr(actions, "publish_event", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(actions, "_JOB_STORE", Store())
    actions._BROKER_SHUTTING_DOWN.clear()
    try:
        result = actions.stop_all_broker_work(reason="test shutdown")
    finally:
        actions._BROKER_SHUTTING_DOWN.clear()

    assert stopped == ["profile"]
    assert result["requested"] == ["Profile"]
    assert result["errors"] == []
