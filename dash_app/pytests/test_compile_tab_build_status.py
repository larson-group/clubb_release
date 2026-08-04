import json

from dash_app.compile_tab import callbacks
from dash_app.compile_tab import runtime


def _configured_build(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    src = repo / "src"
    src.mkdir(parents=True)
    source_file = src / "foo.F90"
    config_file = repo / "CMakeLists.txt"
    source_file.write_text("program foo\nend program foo\n")
    config_file.write_text("cmake_minimum_required(VERSION 3.20)\n")
    build = tmp_path / "build"
    build.mkdir()
    configure_marker = build / "build.ninja"
    configure_marker.write_text("# ninja\n")
    artifact = build / "clubb_standalone"
    artifact.write_text("#!/bin/sh\n")
    artifact.chmod(0o755)
    monkeypatch.setattr(runtime, "REPO_ROOT", str(repo))
    return build, source_file, config_file, configure_marker, artifact


def test_directory_scan_current_when_artifact_is_newer_than_sources(tmp_path, monkeypatch):
    build, source_file, config_file, configure_marker, artifact = _configured_build(tmp_path, monkeypatch)
    runtime.os.utime(config_file, ns=(1_000_000_000, 1_000_000_000))
    runtime.os.utime(source_file, ns=(2_000_000_000, 2_000_000_000))
    runtime.os.utime(configure_marker, ns=(3_000_000_000, 3_000_000_000))
    runtime.os.utime(artifact, ns=(4_000_000_000, 4_000_000_000))

    status = runtime._detect_target_status(str(build), {}, {}, ["clubb_standalone"])

    assert status["status"] == "current"


def test_directory_scan_source_newer_than_artifact_needs_rebuild(tmp_path, monkeypatch):
    build, source_file, config_file, configure_marker, artifact = _configured_build(tmp_path, monkeypatch)
    runtime.os.utime(config_file, ns=(1_000_000_000, 1_000_000_000))
    runtime.os.utime(configure_marker, ns=(2_000_000_000, 2_000_000_000))
    runtime.os.utime(artifact, ns=(3_000_000_000, 3_000_000_000))
    runtime.os.utime(source_file, ns=(4_000_000_000, 4_000_000_000))

    status = runtime._detect_target_status(str(build), {}, {}, ["clubb_standalone"])

    assert status["status"] == "needs_rebuild"


def test_install_artifact_allows_timestamp_granularity_and_rewritten_binary(tmp_path):
    build_artifact = tmp_path / "build_exe"
    install_artifact = tmp_path / "install_exe"
    build_artifact.write_bytes(b"build executable bytes")
    install_artifact.write_bytes(b"install rpath rewrite!")

    build_ns = 10_900_000_000
    install_ns = 10_000_000_000
    runtime.os.utime(build_artifact, ns=(build_ns, build_ns))
    runtime.os.utime(install_artifact, ns=(install_ns, install_ns))

    assert runtime._artifact_is_stale(build_artifact, install_artifact) is False


def test_install_artifact_older_than_tolerance_is_stale(tmp_path):
    build_artifact = tmp_path / "build_exe"
    install_artifact = tmp_path / "install_exe"
    build_artifact.write_bytes(b"build executable bytes")
    install_artifact.write_bytes(b"install rpath rewrite!")

    build_ns = 12_000_000_000
    install_ns = 10_000_000_000
    runtime.os.utime(build_artifact, ns=(build_ns, build_ns))
    runtime.os.utime(install_artifact, ns=(install_ns, install_ns))

    assert runtime._artifact_is_stale(build_artifact, install_artifact) is True


def test_rebuild_command_preserves_discovered_tuning_flag():
    options = runtime.compile_options_from_build(
        {
            "build_type": "Debug",
            "name": "gcc_DEBUG_PRECdouble_PYTHON_OPENMP_TUNING_GPTL",
            "precision": "double",
            "gpu": "none",
            "python": "ON",
            "openmp": "ON",
            "tuning": "ON",
            "gptl": "ON",
            "toolchain": "/tmp/toolchain.cmake",
            "install_prefix": "/tmp/gcc_DEBUG_PRECdouble_PYTHON_OPENMP_TUNING_GPTL",
        },
        {"toolchains": [{"path": "/tmp/toolchain.cmake"}]},
        module_stack=[],
    )

    argv = runtime.build_compile_argv(options)

    assert "-tuning" in argv
    assert "-gptl" in argv
    assert "-install" in argv


def test_rebuild_command_drops_mismatched_install_prefix():
    options = runtime.compile_options_from_build(
        {
            "build_type": "Debug",
            "name": "gcc_DEBUG_PRECdouble_PYTHON_OPENMP_GPTL",
            "precision": "double",
            "gpu": "none",
            "python": "ON",
            "openmp": "ON",
            "tuning": "OFF",
            "gptl": "ON",
            "toolchain": "/tmp/toolchain.cmake",
            "install_prefix": "/tmp/gcc_DEBUG_PRECdouble_PYTHON_OPENMP_TUNING_GPTL",
        },
        {"toolchains": [{"path": "/tmp/toolchain.cmake"}]},
        module_stack=[],
    )

    argv = runtime.build_compile_argv(options)

    assert "-install" not in argv
    assert "-tuning" not in argv
    assert "-gptl" in argv


def test_resolved_rebuild_failure_is_hidden_after_current_status():
    failures = {"/tmp/build": {"returncode": 1, "failed_at": 10.0}}
    statuses = {
        "checked_at": 20.0,
        "statuses": {"/tmp/build": {"status": "current"}},
    }

    assert callbacks.visible_failed_rebuild_paths(failures, statuses) == set()


def test_rebuild_failure_stays_visible_until_fresh_current_status():
    failures = {"/tmp/build": {"returncode": 1, "failed_at": 20.0}}
    statuses = {
        "checked_at": 10.0,
        "statuses": {"/tmp/build": {"status": "current"}},
    }

    assert callbacks.visible_failed_rebuild_paths(failures, statuses) == {"/tmp/build"}


def test_rebuild_returncode_recovers_progress_complete(tmp_path):
    progress = tmp_path / "progress.json"
    progress.write_text(
        json.dumps(
            {
                "current_path": None,
                "completed_paths": ["/tmp/build-a", "/tmp/build-b"],
                "queued_paths": [],
                "state": "complete",
            }
        )
    )
    job = {"kind": "rebuild", "status": "running", "progress": str(progress), "pid": 99999999}

    assert runtime.poll_compile_job(job) == 0
    assert runtime.job_process_is_live(job) is False


def test_rebuild_returncode_recovers_progress_failure(tmp_path):
    progress = tmp_path / "progress.json"
    progress.write_text(json.dumps({"state": "failed"}))
    job = {"kind": "rebuild", "status": "running", "progress": str(progress), "pid": 99999999}

    assert runtime.poll_compile_job(job) == 1
    assert runtime.job_process_is_live(job) is False


def test_toolchain_helpers_choose_the_matching_host_compiler():
    discovery = {
        "toolchains": [
            {"path": "/toolchains/gcc.cmake", "compiler": "gcc", "matches_host": True},
            {"path": "/toolchains/intel.cmake", "compiler": "intel", "matches_host": True},
        ]
    }

    assert runtime.toolchain_compiler(discovery, "/toolchains/gcc.cmake") == "gcc"
    assert runtime.preferred_toolchain_for_compiler(discovery, "intel") == "/toolchains/intel.cmake"


def test_toolchain_defaults_to_auto_and_preserves_an_explicit_override():
    options = [
        {"label": "Auto from environment", "value": "auto"},
        {"label": "gcc", "value": "/toolchains/gcc.cmake"},
        {"label": "intel", "value": "/toolchains/intel.cmake"},
    ]

    assert callbacks.retained_toolchain_value(options, None) == "auto"
    assert callbacks.retained_toolchain_value(options, "/toolchains/intel.cmake") == "/toolchains/intel.cmake"
    assert callbacks.toolchain_selection_mode("auto") == "auto"
    assert callbacks.toolchain_selection_mode("/toolchains/gcc.cmake") == "manual"


def test_ui_gptl_flag_uses_the_shared_compile_option_path():
    options = callbacks.collect_options("double", "none", "auto", ["gptl"], "")

    assert options["gptl"] is True
    assert "-gptl" in runtime.build_compile_argv(options)
