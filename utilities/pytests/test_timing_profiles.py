import csv
import io
import json
import zipfile
from pathlib import Path

import pytest

from utilities.timing_profiles import (
    BATCH_FIELDS,
    TIMING_FIELDS,
    ProfileFormatError,
    create_profile_manifest,
    discover_profiles,
    export_profiles,
    import_profiles,
    load_profiles,
    profile_directory,
    reserve_profile_directory,
    update_profile_manifest,
    write_batches,
    write_profile_manifest,
    write_profile_readme,
)


def write_csv(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_profile(root: Path, run_id="baseline"):
    profile = profile_directory(root, run_id)
    profile.mkdir(parents=True)
    manifest = create_profile_manifest(
        run_id=run_id,
        label="Baseline",
        case_name="arm",
        started_utc="2026-08-12T12:00:00+00:00",
        git_revision="abc123",
        git_dirty=False,
        forwarded_args="-max_iters 5",
        omp_num_threads="1",
        benchmark={
            "processes": [2],
            "batch_sizes": [4],
            "warmups": 1,
            "repetitions": 1,
        },
        environment={"hostname": "test-host", "cpu_model": "test-cpu"},
        build={"resolved_executable": "/tmp/clubb", "executable_sha256": "sha"},
    )
    manifest["results"].update(
        {
            "timing_rows": 5,
            "batches_completed": 1,
            "observed_model_steps": [3],
            "timer_backends": ["cpu_time"],
            "time_bases": ["cpu_time"],
        }
    )
    write_profile_manifest(profile, manifest)
    write_profile_readme(profile, manifest)
    write_batches(
        profile,
        [
            {
                "batch_id": "p0002_b000004",
                "process_count": 2,
                "batch_size": 4,
                "total_batch_size": 8,
                "status": "complete",
                "warmups_completed": 1,
                "measurements_completed": 1,
                "failed_runs": 0,
                "representative_phase": "measured",
                "representative_repetition": 1,
                "representative_process": 0,
            }
        ],
    )
    rows = []
    for phase, repetition, process, value in (
        ("warmup", 1, 0, 2.5),
        ("measured", 1, 0, 2.0),
        ("measured", 1, 1, 2.2),
    ):
        rows.append(
            {
                "batch_id": "p0002_b000004",
                "phase": phase,
                "repetition": repetition,
                "process": process,
                "thread": 0,
                "timer": "advance_clubb_to_end",
                "calls": 1,
                "inclusive_s": value,
                "exclusive_s": value / 2,
                "status": "success",
                "return_code": 0,
                "message": "",
            }
        )
    rows.append(
        {
            "batch_id": "p0002_b000004",
            "phase": "measured",
            "repetition": 1,
            "process": "",
            "thread": "",
            "timer": "__process_group_wall__",
            "calls": 1,
            "inclusive_s": 2.4,
            "exclusive_s": "",
            "status": "success",
            "return_code": "",
            "message": "",
        }
    )
    write_csv(profile / "timings.csv", TIMING_FIELDS, rows)
    logs = profile / "logs" / "p0002_b000004"
    logs.mkdir(parents=True)
    for suffix in ("in", "log", "timing"):
        (logs / f"arm.{suffix}").write_text(suffix, encoding="utf-8")
    (logs / "arm_setup.txt").write_text("setup", encoding="utf-8")
    update_profile_manifest(profile, manifest, status="complete", completed=True)
    return profile


def test_profile_directory_is_readable_and_collision_safe(tmp_path):
    first_id, first = reserve_profile_directory(tmp_path, "CPU baseline")
    second_id, second = reserve_profile_directory(tmp_path, "CPU baseline")

    assert first_id == "CPU_baseline"
    assert second_id.startswith("CPU_baseline_")
    assert first.is_dir() and second.is_dir()


def test_profile_directory_overwrite_reuses_exact_id_and_rejects_unknown_content(tmp_path):
    profile = make_profile(tmp_path, "CPU_baseline")
    marker = profile / "old-result.txt"
    marker.write_text("old", encoding="utf-8")

    run_id, replacement = reserve_profile_directory(tmp_path, "CPU baseline", overwrite=True)

    assert run_id == "CPU_baseline"
    assert replacement == profile
    assert replacement.is_dir()
    assert not marker.exists()
    assert not list(replacement.iterdir())

    conflict = tmp_path / "not_a_profile"
    conflict.mkdir()
    with pytest.raises(ProfileFormatError, match="unrecognized profile directory"):
        reserve_profile_directory(tmp_path, "not a profile", overwrite=True)


def test_compact_profile_derives_dash_summary_and_ignores_warmups(tmp_path):
    make_profile(tmp_path)

    summaries, processes = load_profiles(tmp_path, ["baseline"])

    assert len(summaries) == 1
    assert summaries[0]["columns_per_process"] == 4
    assert summaries[0]["total_columns"] == 8
    assert summaries[0]["timer_mean_seconds"] == pytest.approx(2.1)
    assert summaries[0]["timer_max_seconds"] == pytest.approx(2.2)
    assert summaries[0]["group_wall_seconds"] == pytest.approx(2.4)
    assert len(processes) == 2
    assert all(row["repetition"] == 1 for row in processes)
    assert processes[0]["timing_file"].endswith("arm.timing")


def test_profile_export_import_round_trip_preserves_compact_files(tmp_path):
    source = tmp_path / "source"
    make_profile(source)

    payload = export_profiles(source, ["baseline"])
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        names = archive.namelist()
        assert "clubb-profile-export.json" in names
        assert "profiles/baseline/profile.json" in names
        assert "profiles/baseline/batches.csv" in names
        assert "profiles/baseline/timings.csv" in names
        assert "profiles/baseline/logs/p0002_b000004/arm.timing" in names

    destination = tmp_path / "destination"
    imported = import_profiles(destination, "saved-profile.zip", payload)
    assert imported == ["baseline"]
    summaries, processes = load_profiles(destination, imported)
    assert summaries[0]["timer_max_seconds"] == pytest.approx(2.2)
    assert len(processes) == 2
    assert discover_profiles(destination)[0]["imported_from"] == "saved-profile.zip"

    duplicate = import_profiles(destination, "saved-profile.zip", payload)
    assert duplicate[0].startswith("baseline_imported_")
    assert len(discover_profiles(destination)) == 2


def test_import_rejects_csv_without_profile_metadata(tmp_path):
    with pytest.raises(ProfileFormatError, match="complete CLUBB timing profile ZIP"):
        import_profiles(tmp_path, "timings.csv", b"batch_id,timer\na,t\n")


def test_import_rejects_archive_path_traversal(tmp_path):
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w") as archive:
        archive.writestr(
            "profiles/run/profile.json",
            json.dumps({"format": "clubb-timing-profile", "format_version": 2, "run_id": "run"}),
        )
        archive.writestr("profiles/run/../outside", "unsafe")

    with pytest.raises(ProfileFormatError, match="unsafe archive path"):
        import_profiles(tmp_path, "unsafe.zip", stream.getvalue())


def test_discovery_reads_direct_profile_directory(tmp_path):
    profile = make_profile(tmp_path)

    records = discover_profiles(profile)

    assert records[0]["run_id"] == "baseline"
    assert records[0]["storage"] == "compact"
    assert records[0]["batch_sizes"] == "4"
    assert (profile / "README.md").read_text(encoding="utf-8").startswith("# Baseline")
