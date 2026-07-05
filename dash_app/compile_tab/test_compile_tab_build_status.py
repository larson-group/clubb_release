from types import SimpleNamespace

from dash_app.compile_tab import runtime


def _mock_run(monkeypatch, output, returncode=0):
    def fake_run(*_args, **_kwargs):
        return SimpleNamespace(stdout=output, returncode=returncode)

    monkeypatch.setattr(runtime.subprocess, "run", fake_run)


def _ninja_status():
    return runtime._detect_target_status(
        "/tmp/build",
        {},
        {"CMAKE_GENERATOR": "Ninja", "CMAKE_MAKE_PROGRAM": "/bin/echo"},
        ["clubb_standalone"],
    )


def test_ninja_no_work_takes_priority_over_pgcuda_diagnostics(monkeypatch):
    _mock_run(
        monkeypatch,
        "ninja: Entering directory `/tmp/build`\n"
        "pgcuda diagnostic from cached toolchain state\n"
        "ninja: no work to do.\n",
    )

    status = _ninja_status()

    assert status["status"] == "current"


def test_ninja_planned_pgcuda_work_is_needs_rebuild(monkeypatch):
    _mock_run(
        monkeypatch,
        "ninja: Entering directory `/tmp/build`\n"
        "[1/1] pgcuda -c generated_device_file.cuf\n",
    )

    status = _ninja_status()

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
