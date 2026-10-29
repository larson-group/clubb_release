"""Focused contract tests for reusable CLUBB output discovery."""

import os

from dash_app.shared.output_discovery import (
    discover_stats_directories,
    discover_stats_files,
)


def _write_signature(path):
    path.write_bytes(b"CDF\x01\x00\x00\x00\x00")


def test_trajectory_discovery_reuses_only_unchanged_file_checks(tmp_path, monkeypatch):
    from dash_app.misc_tab.mixing_length_trajectories import analysis

    output = tmp_path / "output"
    output.mkdir()
    path = output / "arm_stats.nc"
    _write_signature(path)
    calls = []
    monkeypatch.setattr(analysis, "inspect_dataset", lambda path, **kwargs: calls.append(path))
    analysis._trajectory_compatible.cache_clear()
    assert len(analysis.discover_netcdf_file_records(tmp_path)) == 1
    assert len(analysis.discover_netcdf_file_records(tmp_path)) == 1
    assert len(calls) == 1
    path.write_bytes(path.read_bytes() + b"updated")
    assert len(analysis.discover_netcdf_file_records(tmp_path)) == 1
    assert len(calls) == 2
    replacement = output / "replacement.nc"
    replacement.write_bytes(path.read_bytes())
    replacement.replace(path)
    assert len(analysis.discover_netcdf_file_records(tmp_path)) == 1
    assert len(calls) == 3


def test_trajectory_discovery_retries_transient_io_failures(tmp_path, monkeypatch):
    from dash_app.misc_tab.mixing_length_trajectories import analysis

    output = tmp_path / "output"
    output.mkdir()
    _write_signature(output / "arm_stats.nc")
    attempts = []
    def inspect(path, **kwargs):
        attempts.append(path)
        if len(attempts) == 1:
            raise OSError("writer temporarily holds the file")
    monkeypatch.setattr(analysis, "inspect_dataset", inspect)
    analysis._trajectory_compatible.cache_clear()
    assert analysis.discover_netcdf_file_records(tmp_path) == []
    assert len(analysis.discover_netcdf_file_records(tmp_path)) == 1
    assert len(attempts) == 2


def test_directory_and_file_discovery_share_direct_stats_contract(tmp_path):
    run = tmp_path / "nested" / "run"
    run.mkdir(parents=True)
    _write_signature(run / "arm_stats.nc")
    _write_signature(run / "bomex_stats.nc")

    directories = discover_stats_directories([tmp_path])
    files = discover_stats_files([tmp_path])

    assert len(directories) == 1
    assert directories[0]["path"] == str(run)
    assert directories[0]["case_names"] == ["arm", "bomex"]
    assert directories[0]["case_count"] == 2
    assert directories[0]["cases"]["arm"] == str(run / "arm_stats.nc")
    assert directories[0]["available"] is True
    assert directories[0]["is_newest"] is True
    assert directories[0]["recency_index"] == 0
    assert directories[0]["case_fingerprints"]["arm"]["size_bytes"] == 8
    assert directories[0]["latest_stats_modified_ns"] > 0
    assert [item["path"] for item in files] == [
        str(run / "arm_stats.nc"),
        str(run / "bomex_stats.nc"),
    ]
    assert len({item["key"] for item in files}) == 2
    assert len({item["label"] for item in files}) == 2


def test_labels_widen_when_multiple_roots_have_the_same_relative_path(tmp_path):
    first = tmp_path / "one" / "output"
    second = tmp_path / "two" / "output"
    for root in (first, second):
        run = root / "run"
        run.mkdir(parents=True)
        _write_signature(run / "arm_stats.nc")

    files = discover_stats_files([first, second])

    assert len(files) == 2
    assert len({item["label"] for item in files}) == 2
    assert all(item["label"].startswith("/") for item in files)


def test_excluded_roots_are_not_descended_into(tmp_path):
    ordinary = tmp_path / "dash_default"
    artifact = tmp_path / "agent_artifacts" / "run"
    ordinary.mkdir()
    artifact.mkdir(parents=True)
    _write_signature(ordinary / "arm_stats.nc")
    _write_signature(artifact / "arm_stats.nc")

    directories = discover_stats_directories([tmp_path])

    assert [record["path"] for record in directories] == [str(ordinary)]


def test_directory_limit_is_applied_after_newest_stats_ordering(tmp_path):
    runs = []
    for index in range(3):
        run = tmp_path / f"run-{index}"
        run.mkdir()
        stats = run / "arm_stats.nc"
        _write_signature(stats)
        os.utime(stats, ns=(index + 1, index + 1))
        runs.append(run)

    directories = discover_stats_directories([tmp_path], max_directories=2)

    assert [record["path"] for record in directories] == [str(runs[2]), str(runs[1])]
    assert [record["recency_index"] for record in directories] == [0, 1]


def test_invalid_stats_files_are_not_counted(tmp_path):
    run = tmp_path / "mixed"
    run.mkdir()
    _write_signature(run / "arm_stats.nc")
    (run / "bomex_stats.nc").write_text("not netcdf")

    record = discover_stats_directories([tmp_path])[0]

    assert record["case_count"] == 1
    assert record["case_names"] == ["arm"]
