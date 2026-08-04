"""Focused tests for the output-directory discovery service."""

import os

from dash_app.services.profiles import discover_output_directories


def _write_netcdf_signature(path):
    path.write_bytes(b"CDF\x01\x00\x00\x00\x00")


def test_discover_output_directories_finds_direct_stats_files_and_sorts_newest(tmp_path):
    older = tmp_path / "older"
    newer = tmp_path / "nested" / "newer"
    ignored = tmp_path / "not-a-run"
    older.mkdir()
    newer.mkdir(parents=True)
    ignored.mkdir()
    _write_netcdf_signature(older / "arm_stats.nc")
    _write_netcdf_signature(newer / "bomex_stats.nc")
    (ignored / "notes.txt").write_text("not output")
    os.utime(older, (100.0, 100.0))
    os.utime(newer, (200.0, 200.0))

    records = discover_output_directories(tmp_path)

    assert [record["relative_path"] for record in records] == [
        "output/nested/newer",
        "output/older",
    ]
    assert records[0]["case_names"] == ["bomex"]
    assert records[1]["case_names"] == ["arm"]


def test_discover_output_directories_rejects_unreadable_stats_signature(tmp_path):
    invalid = tmp_path / "invalid"
    invalid.mkdir()
    (invalid / "arm_stats.nc").write_text("not netcdf")

    assert discover_output_directories(tmp_path) == []


def test_discover_output_directories_skips_agent_artifacts_even_when_they_contain_stats(tmp_path):
    ordinary = tmp_path / "dash_default"
    artifact = tmp_path / "agent_artifacts" / "investigation" / "baseline"
    ordinary.mkdir()
    artifact.mkdir(parents=True)
    _write_netcdf_signature(ordinary / "arm_stats.nc")
    _write_netcdf_signature(artifact / "arm_stats.nc")

    records = discover_output_directories(tmp_path)

    assert [record["relative_path"] for record in records] == ["output/dash_default"]
