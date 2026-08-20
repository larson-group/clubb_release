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
    os.utime(older / "arm_stats.nc", (100.0, 100.0))
    os.utime(newer / "bomex_stats.nc", (200.0, 200.0))

    records = discover_output_directories(tmp_path)

    assert [record["relative_path"] for record in records] == [
        "output/nested/newer",
        "output/older",
    ]
    assert records[0]["case_names"] == ["bomex"]
    assert records[1]["case_names"] == ["arm"]
    assert records[0]["label"] == "output/nested/newer"
    assert records[0]["case_count"] == 1


def test_discover_output_directories_rejects_unreadable_stats_signature(tmp_path):
    invalid = tmp_path / "invalid"
    invalid.mkdir()
    (invalid / "arm_stats.nc").write_text("not netcdf")

    assert discover_output_directories(tmp_path) == []


def test_discover_output_directories_skips_agent_and_mcp_artifacts(tmp_path):
    ordinary = tmp_path / "dash_default"
    artifact = tmp_path / "agent_artifacts" / "investigation" / "baseline"
    mcp_run = tmp_path / "mcp_runs" / "agent-run"
    ordinary.mkdir()
    artifact.mkdir(parents=True)
    mcp_run.mkdir(parents=True)
    _write_netcdf_signature(ordinary / "arm_stats.nc")
    _write_netcdf_signature(artifact / "arm_stats.nc")
    _write_netcdf_signature(mcp_run / "arm_stats.nc")

    records = discover_output_directories(tmp_path)

    assert [record["relative_path"] for record in records] == ["output/dash_default"]


def test_explicitly_selected_mcp_output_remains_available(tmp_path):
    mcp_run = tmp_path / "mcp_runs" / "agent-run"
    mcp_run.mkdir(parents=True)
    _write_netcdf_signature(mcp_run / "arm_stats.nc")

    records = discover_output_directories(tmp_path, selected_dirs=[str(mcp_run)])

    assert [record["path"] for record in records] == [str(mcp_run.resolve())]
    assert records[0]["available"] is True


def test_selected_external_directory_is_included_with_direct_case_metadata(tmp_path):
    output_root = tmp_path / "output"
    external = tmp_path / "external-run"
    output_root.mkdir()
    external.mkdir()
    _write_netcdf_signature(external / "arm_stats.nc")

    records = discover_output_directories(output_root, selected_dirs=[str(external)])

    assert len(records) == 1
    assert records[0]["path"] == str(external.resolve())
    assert records[0]["label"] == str(external.resolve())
    assert records[0]["available"] is True
    assert records[0]["case_count"] == 1


def test_selected_missing_directory_remains_as_unavailable_catalog_record(tmp_path):
    output_root = tmp_path / "output"
    output_root.mkdir()
    missing = tmp_path / "missing"

    records = discover_output_directories(output_root, selected_dirs=[str(missing)])

    assert records[0]["path"] == str(missing.resolve())
    assert records[0]["available"] is False
    assert records[0]["case_count"] == 0
