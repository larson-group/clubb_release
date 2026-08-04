"""Focused contract tests for reusable CLUBB output discovery."""

from dash_app.shared.output_discovery import (
    discover_stats_directories,
    discover_stats_files,
)


def _write_signature(path):
    path.write_bytes(b"CDF\x01\x00\x00\x00\x00")


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
    assert directories[0]["cases"]["arm"] == str(run / "arm_stats.nc")
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
