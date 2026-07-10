"""Compile subprocess launching, cancellation, and log helpers."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import time

from .discovery import (
    canonical_compiler_from_path,
    command_in_env,
    compiler_from_env,
    parse_cmake_cache,
    resolve_lmod_stack,
    install_prefix_matches_build_name,
)
from .state import BUILD_DIR, COMPILE_LOCK, COMPILE_PROC, MAX_UI_LOG_LINES, REPO_ROOT

INSTALL_MTIME_TOLERANCE_NS = 1_000_000_000
RUNTIME_BUILD_TARGETS = (
    "clubb_driver_lib",
    "clubb_standalone",
    "clubb_thread_test",
    "clubb_tuner",
    "G_unit_tests",
    "clubb_driver_test",
    "clubb_standalone_loss",
    "clubb_loss_driver_test",
)
PYTHON_BUILD_TARGETS = ("clubb_f2py",)
RUNTIME_INSTALL_ARTIFACTS = {
    "clubb_driver_lib": ("src/libclubb_driver_lib.a", "libclubb_driver_lib.a"),
    "clubb_standalone": ("src/clubb_standalone", "clubb_standalone"),
    "clubb_thread_test": ("src/clubb_thread_test", "clubb_thread_test"),
    "clubb_tuner": ("src/clubb_tuner", "clubb_tuner"),
    "G_unit_tests": ("src/G_unit_tests", "G_unit_tests"),
    "clubb_driver_test": ("src/clubb_driver_test", "clubb_driver_test"),
    "clubb_standalone_loss": ("src/clubb_standalone_loss", "clubb_standalone_loss"),
    "clubb_loss_driver_test": ("src/clubb_loss_driver_test", "clubb_loss_driver_test"),
}
CLUBB_STANDARDS_SCRIPT = Path(REPO_ROOT) / "utilities" / "CLUBBStandardsCheck.py"
MAX_SOURCE_CHECK_LOG_BYTES = 2_000_000
CLUBB_STANDARDS_CHECK_TARGETS = (
    Path(REPO_ROOT) / "src" / "clubb_driver.F90",
    Path(REPO_ROOT) / "src" / "CLUBB_core",
    Path(REPO_ROOT) / "src" / "SILHS",
    Path(REPO_ROOT) / "src" / "Benchmark_cases",
    Path(REPO_ROOT) / "src" / "Radiation",
    Path(REPO_ROOT) / "src" / "Microphys",
    Path(REPO_ROOT) / "src" / "Microphys" / "KK_microphys",
    Path(REPO_ROOT) / "src" / "G_unit_test_types",
)

_REBUILD_HELPER = r"""
import json
import os
import subprocess
import sys

compile_specs = json.loads(sys.argv[1])
progress_path = sys.argv[2]


def write_progress(current_path, completed_paths, state):
    payload = {
        "current_path": current_path,
        "completed_paths": completed_paths,
        "queued_paths": [
            spec["build_path"]
            for spec in compile_specs
            if spec["build_path"] != current_path and spec["build_path"] not in completed_paths
        ],
        "state": state,
    }
    try:
        with open(progress_path, "w", encoding="utf-8") as progress_file:
            json.dump(payload, progress_file)
    except OSError:
        pass


def display_command(cmd):
    return "+ " + " ".join(cmd)


def write_line(build_log, line):
    print(line, flush=True)
    build_log.write(line + "\n")
    build_log.flush()


def stream_command(cmd, build_log, env):
    write_line(build_log, display_command(cmd))
    proc = subprocess.Popen(
        cmd,
        cwd=os.getcwd(),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    for line in proc.stdout:
        print(line, end="", flush=True)
        build_log.write(line)
        build_log.flush()
    return proc.wait()


completed_paths = []
write_progress(None, completed_paths, "queued")

for index, spec in enumerate(compile_specs, start=1):
    cmd = spec["command"]
    env = spec.get("env") or os.environ.copy()
    build_path = spec["build_path"]
    build_log_path = spec["build_log_path"]
    write_progress(build_path, completed_paths, "running")
    header_lines = [
        f"=== [{index}/{len(compile_specs)}] Reconfiguring and rebuilding {build_path} ===",
        f"=== Writing build log: {build_log_path} ===",
    ]
    returncode = 0
    try:
        with open(build_log_path, "a", encoding="utf-8", errors="replace") as build_log:
            build_log.write("\n" + "\n".join(header_lines) + "\n")
            build_log.flush()
            for line in header_lines:
                print(line, flush=True)
            returncode = stream_command(cmd, build_log, env)
            if returncode == 0:
                write_line(build_log, f"=== Rebuild complete for {build_path} ===")
                completed_paths.append(build_path)
                write_progress(None, completed_paths, "queued")
    except OSError as exc:
        returncode = 1
        print(f"=== Rebuild failed for {build_path}: {exc} ===", flush=True)

    if returncode != 0:
        failure = f"=== Rebuild failed for {build_path} with exit {returncode} ==="
        print(failure, flush=True)
        try:
            with open(build_log_path, "a", encoding="utf-8", errors="replace") as build_log:
                build_log.write(failure + "\n")
        except OSError:
            pass
        sys.exit(returncode)

write_progress(None, completed_paths, "complete")
print("=== Rebuild complete ===", flush=True)
"""


def _bool_flag(options, key, flag):
    return [flag] if options.get(key) else []


def clubb_standards_failed_message(log_path):
    """Return the same source-check failure warning text that compile.py prints."""
    return "\n".join(
        [
            "=============================================================",
            "CLUBBStandardsCheck FAILED",
            "  THIS IS PRINTED IN ALL RED, CAPITAL LETTERS, AND USES",
            "  AN EXCLAMATION MARK TO ENSURE THE DEVELOPERS FEEL SHAME!",
            "  IF YOU ARE ONE OF THESE \"DEVELOPERS\" CHECK THE",
            f"  LOG FILE FOR DETAILS: {log_path}",
            "=============================================================",
        ]
    )


def run_clubb_standards_check(log_path=None):
    """Run CLUBBStandardsCheck.py over the same sources as compile.py."""
    if log_path is None:
        log_file = tempfile.NamedTemporaryFile(delete=False, prefix="clubb_source_check_", suffix=".log", dir="/tmp")
        log_path = log_file.name
        log_file.close()
    retcode = 0
    with open(log_path, "a", encoding="utf-8", errors="replace") as log:
        for target in CLUBB_STANDARDS_CHECK_TARGETS:
            if target.is_dir():
                files = [str(path) for path in target.glob("*.F90")]
            elif target.is_file():
                files = [str(target)]
            else:
                log.write(f"No matches for {target}\n")
                continue
            if not files:
                log.write(f"No matches for {target}\n")
                continue
            cmd = [sys.executable, str(CLUBB_STANDARDS_SCRIPT), *files]
            log.write("\n================= Running: " + " ".join(cmd) + " =================\n")
            log.flush()
            try:
                process = subprocess.Popen(
                    cmd,
                    cwd=REPO_ROOT,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                for line in process.stdout or ():
                    log.write(line)
                    log.flush()
                process.wait()
                retcode += process.returncode
            except OSError as exc:
                log.write(f"Failed to run CLUBBStandardsCheck.py: {exc}\n")
                retcode += 1
            log.write("\n")
            log.flush()
    return {
        "returncode": retcode,
        "log": str(log_path),
        "warning": clubb_standards_failed_message(log_path) if retcode != 0 else "",
        "checked_at": time.time(),
    }


def validated_source_check_log_path(log_path):
    """Return a source-check log path that is safe for the Dash app to open."""
    if not log_path:
        raise RuntimeError("No source-check log file is available.")
    path = Path(log_path).expanduser().resolve()
    # Source-check logs are written to /tmp; macOS resolves that directory as /private/tmp.
    tmp_dirs = {Path(tempfile.gettempdir()).resolve(), Path("/tmp").resolve()}
    if (
        path.parent not in tmp_dirs
        or not path.name.startswith("clubb_source_check_")
        or path.suffix != ".log"
    ):
        raise RuntimeError(f"Refusing to open non-source-check log path: {path}")
    if not path.is_file():
        raise RuntimeError(f"Source-check log file does not exist: {path}")
    return path


def read_source_check_log(log_path):
    """Read a source-check log file for display inside the Dash app."""
    path = validated_source_check_log_path(log_path)
    size = path.stat().st_size
    try:
        with path.open("rb") as log_file:
            if size > MAX_SOURCE_CHECK_LOG_BYTES:
                log_file.seek(size - MAX_SOURCE_CHECK_LOG_BYTES)
                raw = log_file.read()
                content = raw.decode("utf-8", errors="replace")
                content = (
                    f"[Showing final {MAX_SOURCE_CHECK_LOG_BYTES:,} bytes of "
                    f"{size:,} total bytes.]\n\n{content}"
                )
            else:
                content = log_file.read().decode("utf-8", errors="replace")
    except OSError as exc:
        raise RuntimeError(f"Could not read source-check log: {exc}") from exc
    return {"path": str(path), "content": content}


COMPILE_FEATURE_FLAG_ORDER = [
    "debug",
    "python",
    "openmp",
    "tuning",
    "gptl"
]


def cmake_enabled(value):
    """Return whether a CMake cache value represents ON."""
    return str(value or "").strip().upper() in {"1", "ON", "TRUE", "YES"}


def cmake_disabled(value):
    """Return whether a CMake cache value represents OFF."""
    return str(value or "").strip().upper() in {"0", "OFF", "FALSE", "NO"}


def _toolchain_compiler(toolchain_path):
    """Return the compiler suffix from a CLUBB toolchain path."""
    stem = Path(toolchain_path or "").stem
    if "_" not in stem:
        return ""
    return stem.rsplit("_", 1)[-1]


def _build_compiler_family(build):
    """Infer the canonical compiler family for one discovered build."""
    compiler = canonical_compiler_from_path((build or {}).get("fortran_compiler", ""))
    if compiler:
        return compiler
    return _toolchain_compiler((build or {}).get("toolchain", ""))


def lmod_compiler_from_build(build, discovery=None):
    """Return an Lmod compiler module that safely matches a discovered build."""
    modules = [
        item
        for item in (discovery or {}).get("lmod", {}).get("compiler_modules", [])
        if item.get("enabled") is not False and item.get("module")
    ]
    compiler_path = str((build or {}).get("fortran_compiler") or "")
    if compiler_path:
        exact_matches = [
            item
            for item in modules
            if item.get("diagnostics", {}).get("FC") == compiler_path
        ]
        if len(exact_matches) == 1:
            return exact_matches[0]["module"]

    family = _build_compiler_family(build)
    if not family:
        return None
    family_matches = [item for item in modules if item.get("canonical") == family]
    if len(family_matches) == 1:
        return family_matches[0]["module"]
    return None


def compile_settings_from_build(build, discovery=None, lmod_compiler=None):
    """Return compile.py option values recoverable from one discovered build."""
    valid_precisions = {"single", "double"}
    valid_gpus = {"none", "openacc", "openmp"}
    precision = str((build or {}).get("precision") or "double").strip().lower()
    gpu = str((build or {}).get("gpu") or "none").strip().lower()

    if precision not in valid_precisions:
        precision = "double"
    if gpu not in valid_gpus:
        gpu = "none"

    flags = []
    if str((build or {}).get("build_type") or "").strip().lower() == "debug":
        flags.append("debug")
    if cmake_enabled((build or {}).get("python")):
        flags.append("python")
    if cmake_enabled((build or {}).get("openmp")):
        flags.append("openmp")
    if cmake_enabled((build or {}).get("tuning")):
        flags.append("tuning")
    if cmake_enabled((build or {}).get("gptl")):
        flags.append("gptl")

    toolchain = (build or {}).get("toolchain") or "auto"
    valid_toolchains = {"auto"}
    valid_toolchains.update(
        item.get("path")
        for item in (discovery or {}).get("toolchains", [])
        if item.get("path")
    )
    if toolchain not in valid_toolchains:
        toolchain = "auto"

    return {
        "toolchain": toolchain,
        "precision": precision,
        "gpu": gpu,
        "flags": [flag for flag in COMPILE_FEATURE_FLAG_ORDER if flag in flags],
        "extra_args": "",
        "lmod_compiler": lmod_compiler or lmod_compiler_from_build(build, discovery),
    }


def compile_options_from_build(build, discovery=None, module_stack=None):
    """Return compile.py options for rerunning the build through compile.py."""
    settings = compile_settings_from_build(build, discovery, None)
    if module_stack is None:
        lmod_compiler = settings.get("lmod_compiler")
        module_stack = [lmod_compiler] if lmod_compiler else []
    options = {
        "precision": settings["precision"],
        "gpu": settings["gpu"],
        "toolchain": settings["toolchain"],
        "extra_args": "",
        "module_stack": module_stack,
        "debug": "debug" in settings["flags"],
        "run_tests": False,
        "python": "python" in settings["flags"],
        "fresh": False,
        "openmp": "openmp" in settings["flags"],
        "tuning": "tuning" in settings["flags"],
        "gptl": "gptl" in settings["flags"],
    }
    install_prefix = (build or {}).get("install_prefix") or ""
    build_name = (build or {}).get("name") or Path((build or {}).get("path") or "").name
    if install_prefix and install_prefix_matches_build_name(build_name, install_prefix):
        options["install"] = install_prefix
    return options


def selected_environment(discovery, env_id):
    """Return discovery metadata for a selected environment id."""
    for environment in (discovery or {}).get("environments", []):
        if environment.get("id") == env_id:
            return environment
    environments = (discovery or {}).get("environments", [])
    return environments[0] if environments else {"id": "current", "kind": "current", "canonical": "gcc"}


def build_compile_argv(options):
    """Build the compile.py argument vector after the Python executable."""
    argv = [sys.executable, "-u", "compile.py"]
    argv.append("-skip_source_checks")
    precision = options.get("precision") or "double"
    gpu = options.get("gpu") or "none"
    if precision != "double":
        argv.extend(["-precision", precision])
    if gpu != "none":
        argv.extend(["-gpu", gpu])
    toolchain = options.get("toolchain") or "auto"
    if toolchain != "auto":
        argv.extend(["-toolchain", toolchain])
    install = (options.get("install") or "").strip()
    if install:
        argv.extend(["-install", install])
    argv.extend(_bool_flag(options, "debug", "-debug"))
    argv.extend(_bool_flag(options, "run_tests", "-run_tests"))
    argv.extend(_bool_flag(options, "python", "-python"))
    argv.extend(_bool_flag(options, "fresh", "-fresh"))
    argv.extend(_bool_flag(options, "openmp", "-openmp"))
    argv.extend(_bool_flag(options, "tuning", "-tuning"))
    argv.extend(_bool_flag(options, "gptl", "-gptl"))
    extra_args = (options.get("extra_args") or "").strip()
    if extra_args:
        argv.extend(shlex.split(extra_args))
    return argv


def display_module_load(modules):
    """Return a shell-like display string for a module stack."""
    return "module load " + " ".join(shlex.quote(module) for module in modules)


def matching_toolchain_compilers(discovery):
    """Return host-matching toolchain compiler names."""
    return {
        item.get("compiler")
        for item in (discovery or {}).get("toolchains", [])
        if item.get("matches_host")
    }


def validated_build_path(build_path):
    """Return a safe absolute build path inside build/."""
    build_root = Path(BUILD_DIR).resolve()
    path = Path(build_path or "").resolve()
    if not path.is_relative_to(build_root):
        raise RuntimeError(f"Build path is outside {build_root}: {path}")
    if not path.is_dir():
        raise RuntimeError(f"Build directory not found: {path}")
    if not (path / "CMakeCache.txt").is_file():
        raise RuntimeError(f"Build directory has no CMakeCache.txt: {path}")
    return str(path)


def _status(status, label=None, detail="", output="", command=None, returncode=None):
    """Build one freshness status payload."""
    return {
        "status": status,
        "label": label or status.replace("_", " "),
        "detail": detail,
        "output": (output or "")[-4000:],
        "command": " ".join(shlex.quote(part) for part in command) if command else "",
        "returncode": returncode,
        "checked_at": time.time(),
    }


def _build_cache_for(path, build):
    """Return CMake cache metadata from discovery or CMakeCache.txt."""
    cache = {
        "CMAKE_GENERATOR": (build or {}).get("generator", ""),
        "CMAKE_MAKE_PROGRAM": (build or {}).get("make_program", ""),
        "CMAKE_INSTALL_PREFIX": (build or {}).get("install_prefix", ""),
        "ENABLE_F2PY": (build or {}).get("python", ""),
    }
    if not all(cache.values()):
        parsed = parse_cmake_cache(Path(path) / "CMakeCache.txt")
        for key, value in parsed.items():
            cache.setdefault(key, value)
            if not cache[key]:
                cache[key] = value
    return cache


def _runtime_targets(cache):
    """Return CLUBB targets whose freshness determines runtime usability."""
    targets = list(RUNTIME_BUILD_TARGETS)
    if cache.get("ENABLE_F2PY") == "ON":
        targets.extend(PYTHON_BUILD_TARGETS)
    return targets


def _configure_marker_mtime(path):
    """Return the newest CMake-generated marker for the last configure step."""
    build_path = Path(path)
    marker_paths = [
        build_path / "CMakeFiles" / "cmake.check_cache",
        build_path / "build.ninja",
        build_path / "Makefile",
        build_path / "CMakeFiles" / "Makefile.cmake",
        build_path / "cmake_install.cmake",
        build_path / "CMakeCache.txt",
    ]
    newest = (0.0, "")
    for marker_path in marker_paths:
        try:
            mtime = marker_path.stat().st_mtime
        except OSError:
            continue
        if mtime > newest[0]:
            newest = (mtime, str(marker_path))
    return newest


def _is_cmake_config_file(path):
    """Return whether a source-tree file should be treated as CMake config."""
    return path.name == "CMakeLists.txt" or path.suffix == ".cmake"


def _detect_target_status(path, build, cache, targets):
    """Ask whether CLUBB runtime targets need work using directory timestamp scanning."""
    # Get maximum modification time of the build artifacts (executables and libraries)
    build_path_obj = Path(path)
    max_build_time = 0.0
    try:
        for p in build_path_obj.rglob('*'):
            if p.is_file():
                # Skip intermediate files and folders generated by CMake
                if p.suffix in ['.o', '.mod', '.stamp', '.make', '.cmake', '.txt', '.log', '.ts', '.json']:
                    continue
                if 'CMakeFiles' in p.parts:
                    continue
                # If it's a library (.a, .so, .dylib) or an executable file
                if p.suffix in ['.a', '.so', '.dylib'] or (p.stat().st_mode & 0o111):
                    max_build_time = max(max_build_time, p.stat().st_mtime)
    except OSError as exc:
        return _status("unknown", detail=f"Error scanning build directory: {exc}")

    # Scan whitelisted directories and files to find max source and config mtimes
    whitelist_dirs = ["src", "clubb_python_api", "cmake"]
    whitelist_files = ["CMakeLists.txt"]
    
    max_src_time = 0.0
    max_config_time = 0.0
    repo_root = Path(REPO_ROOT)
    
    try:
        # 1. Check whitelisted files at the root
        for f_name in whitelist_files:
            p = repo_root / f_name
            if p.is_file():
                mtime = p.stat().st_mtime
                if _is_cmake_config_file(p):
                    max_config_time = max(max_config_time, mtime)
                else:
                    max_src_time = max(max_src_time, mtime)
                    
        # 2. Check whitelisted directories recursively
        for d_name in whitelist_dirs:
            d_path = repo_root / d_name
            if d_path.is_dir():
                for p in d_path.rglob("*"):
                    if p.is_file():
                        mtime = p.stat().st_mtime
                        # Only track relevant source, header, config, or scripting files
                        if p.suffix in [".f90", ".F90", ".c", ".h", ".cpp", ".F", ".txt", ".cmake", ".py"]:
                            if _is_cmake_config_file(p):
                                max_config_time = max(max_config_time, mtime)
                            else:
                                max_src_time = max(max_src_time, mtime)
    except OSError as exc:
        return _status("unknown", detail=f"Error scanning source directory: {exc}")

    configure_time, configure_marker = _configure_marker_mtime(path)

    if max_config_time > configure_time:
        return _status(
            "needs_configure",
            "needs configure",
            detail="CMake configuration files have changed.",
            output=(
                f"Max config modification time: {time.ctime(max_config_time)}\n"
                f"Last configure marker: {configure_marker or 'not found'}\n"
                f"Last configure marker modification time: {time.ctime(configure_time)}"
            ),
        )
    elif max_build_time == 0.0:
        return _status(
            "needs_rebuild",
            "needs rebuild",
            detail="Build artifacts not found.",
            output="No build artifacts (libraries or executables) found in the build directory.",
        )
    elif max_src_time > max_build_time:
        return _status(
            "needs_rebuild",
            "needs rebuild",
            detail="Source files have changed since the last build.",
            output=f"Max source modification time: {time.ctime(max_src_time)}\nMax build modification time: {time.ctime(max_build_time)}",
        )
    
    # Generate user-facing mock command explanation for compatibility
    command_str = f"Directory scan of source directories: {', '.join(whitelist_dirs)}"
    return _status(
        "current",
        detail="CLUBB runtime targets are current.",
        output=(
            "Directory scan results:\n"
            f"  Max source modification time: {time.ctime(max_src_time)}\n"
            f"  Max build modification time: {time.ctime(max_build_time)}\n"
            f"  Last configure marker: {configure_marker or 'not found'}\n"
            f"  Last configure marker modification time: {time.ctime(configure_time)}"
        ),
        command=[command_str],
        returncode=0,
    )


def _artifact_is_stale(build_artifact, install_artifact):
    """Return whether an installed artifact is missing or older than its build artifact."""
    if not build_artifact.is_file():
        return None
    if not install_artifact.is_file():
        return True
    try:
        build_stat = build_artifact.stat()
        install_stat = install_artifact.stat()
    except OSError:
        return True
    if install_stat.st_size != build_stat.st_size:
        return True
    if install_stat.st_mtime_ns + INSTALL_MTIME_TOLERANCE_NS >= build_stat.st_mtime_ns:
        return False
    return not _files_have_same_contents(build_artifact, install_artifact)


def _files_have_same_contents(left, right):
    """Return whether two same-sized files have identical bytes."""
    try:
        with left.open("rb") as left_file, right.open("rb") as right_file:
            while True:
                left_chunk = left_file.read(1024 * 1024)
                right_chunk = right_file.read(1024 * 1024)
                if left_chunk != right_chunk:
                    return False
                if not left_chunk:
                    return True
    except OSError:
        return False


def _python_install_artifacts(path, install_prefix):
    """Yield Python/F2PY build and install artifact pairs for comparison."""
    runtime_dir = Path(path) / "clubb_python_api" / "f2py_runtime"
    install_dir = Path(install_prefix) / "python"
    for pattern in ("clubb_f2py*", "libclubb_f2py_backend*"):
        for build_artifact in runtime_dir.glob(pattern):
            yield build_artifact, install_dir / build_artifact.name


def _detect_install_status(path, cache, targets):
    """Verify installed CLUBB runtime artifacts are at least as fresh as build artifacts."""
    install_prefix = cache.get("CMAKE_INSTALL_PREFIX", "")
    if not install_prefix:
        return _status("needs_install", "needs install", detail="CMAKE_INSTALL_PREFIX is not set.")
    install_root = Path(install_prefix)
    if not install_root.is_dir():
        return _status("needs_install", "needs install", detail=f"Install directory is missing: {install_prefix}.")

    stale = []
    unknown = []
    for target in targets:
        artifact_pair = RUNTIME_INSTALL_ARTIFACTS.get(target)
        if not artifact_pair:
            continue
        build_rel, install_rel = artifact_pair
        build_artifact = Path(path) / build_rel
        install_artifact = install_root / install_rel
        artifact_stale = _artifact_is_stale(build_artifact, install_artifact)
        if artifact_stale is None:
            unknown.append(str(build_artifact))
        elif artifact_stale:
            stale.append(str(install_artifact))

    if cache.get("ENABLE_F2PY") == "ON":
        python_pairs = list(_python_install_artifacts(path, install_prefix))
        if not python_pairs:
            unknown.append(str(Path(path) / "clubb_python_api" / "f2py_runtime"))
        for build_artifact, install_artifact in python_pairs:
            if _artifact_is_stale(build_artifact, install_artifact):
                stale.append(str(install_artifact))

    if stale:
        return _status(
            "needs_install",
            "needs install",
            detail=f"Installed artifact is stale or missing: {stale[0]}",
        )
    if unknown:
        return _status(
            "unknown",
            detail=f"Build artifact was not found for install freshness check: {unknown[0]}",
        )
    return _status("current", detail="Installed CLUBB runtime artifacts are current.")


def detect_build_status(build):
    """Return whether CLUBB runtime artifacts need rebuild or install work."""
    build_path = build.get("path") if isinstance(build, dict) else build
    try:
        path = validated_build_path(build_path)
    except RuntimeError as exc:
        return _status("unknown", detail=str(exc))

    cache = _build_cache_for(path, build if isinstance(build, dict) else {})
    targets = _runtime_targets(cache)
    install_prefix = cache.get("CMAKE_INSTALL_PREFIX", "")
    build_name = (build.get("name") if isinstance(build, dict) else "") or Path(path).name
    if install_prefix and not install_prefix_matches_build_name(build_name, install_prefix):
        return _status(
            "needs_configure",
            "needs configure",
            detail="CMake install prefix points at a different build directory.",
            output=(
                f"Build name: {build_name}\n"
                f"CMAKE_INSTALL_PREFIX: {install_prefix}"
            ),
        )
    target_status = _detect_target_status(path, build, cache, targets)
    if target_status.get("status") != "current":
        return target_status

    install_status = _detect_install_status(path, cache, targets)
    if install_status.get("status") != "current":
        return {
            **install_status,
            "output": target_status.get("output", ""),
            "command": target_status.get("command", ""),
            "returncode": target_status.get("returncode"),
        }

    return {
        **target_status,
        "detail": "CLUBB runtime targets and installed artifacts are current.",
    }



def detect_build_statuses(discovery):
    """Return dry-run freshness status for all discovered builds."""
    return {
        "checked_at": time.time(),
        "statuses": {
            build.get("path"): detect_build_status(build)
            for build in (discovery or {}).get("builds", [])
            if build.get("path")
        },
    }


def build_compile_spec(discovery, env_id, options):
    """Build subprocess parameters and display command for a compile request."""
    environment = selected_environment(discovery, env_id)
    argv = build_compile_argv(options)
    child_env = os.environ.copy()
    shell = False
    command = argv
    display_prefix = ""
    load_error = ""
    module_stack = [module for module in options.get("module_stack", []) if module]

    if module_stack:
        result = resolve_lmod_stack(module_stack)
        env = result.get("env", {})
        environment = {
            "id": "lmod:" + "|".join(module_stack),
            "kind": "lmod_stack",
            "label": "Lmod module stack",
            "detail": " ".join(module_stack),
            "modules": module_stack,
            "canonical": compiler_from_env(env),
            "compiler_path": env.get("FC", ""),
            "diagnostics": {
                "FC": env.get("FC", ""),
                "CC": env.get("CC", ""),
                "LMOD_FAMILY_COMPILER": env.get("LMOD_FAMILY_COMPILER", ""),
                "LOADEDMODULES": env.get("LOADEDMODULES", ""),
                "cmake": command_in_env("cmake", env) or "",
                "nf-config": command_in_env("nf-config", env) or "",
                "nc-config": command_in_env("nc-config", env) or "",
            },
            "load_ok": bool(result.get("ok")),
            "load_stderr": result.get("stderr", ""),
        }
        if result.get("ok"):
            child_env = env
        else:
            load_error = result.get("stderr", "Cannot load selected module stack.")
        display_prefix = f"{display_module_load(module_stack)} && "
    elif environment.get("kind") == "native" and environment.get("compiler_path"):
        child_env["FC"] = environment["compiler_path"]
        display_prefix = f"FC={shlex.quote(environment['compiler_path'])} "

    display = display_prefix + " ".join(shlex.quote(str(part)) for part in argv)
    return {
        "command": command,
        "display": display,
        "env": child_env,
        "shell": shell,
        "environment": environment,
        "load_error": load_error,
    }


def build_warnings(discovery, env_id, options):
    """Return user-visible compatibility warnings for the current selection."""
    warnings = []
    module_stack = [module for module in options.get("module_stack", []) if module]
    if (discovery or {}).get("lmod", {}).get("available") and not module_stack:
        warnings.append("Select a compiler module before compiling with Lmod.")

    spec = build_compile_spec(discovery, env_id, options)
    environment = spec.get("environment", {})
    if spec.get("load_error"):
        warnings.append(f"Cannot load selected module stack: {spec['load_error']}")
    canonical = environment.get("canonical") or ""
    gpu = options.get("gpu") or "none"
    if gpu != "none" and canonical in {"gcc", "intel"}:
        warnings.append(f"{canonical} toolchains currently reject GPU builds.")
    if environment.get("enabled") is False:
        warnings.append(environment.get("disabled_reason") or "Selected environment is not currently available.")
    if module_stack:
        diagnostics = environment.get("diagnostics", {})
        if not canonical:
            warnings.append("No compiler was detected from the selected module stack.")
        elif (options.get("toolchain") or "auto") == "auto" and canonical not in matching_toolchain_compilers(discovery):
            warnings.append(f"No matching CLUBB toolchain was found for compiler family {canonical}.")
        if not diagnostics.get("cmake"):
            warnings.append("cmake was not found on PATH for the selected module stack.")
        if not (diagnostics.get("nf-config") or diagnostics.get("nc-config")):
            warnings.append("NetCDF config tools were not found in the selected module stack.")
    elif not (discovery or {}).get("tools", {}).get("cmake"):
        warnings.append("cmake was not found on PATH.")
    return warnings


def read_log_increment(log_path, offset):
    """Read new log text from a file starting at byte offset."""
    if not log_path or not os.path.exists(log_path):
        return "", 0
    try:
        with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
            handle.seek(0, os.SEEK_END)
            end_pos = handle.tell()
            if offset < 0 or offset > end_pos:
                offset = 0
            if offset == end_pos:
                return "", end_pos
            handle.seek(offset)
            return handle.read(), handle.tell()
    except OSError:
        return "", offset


def append_log_tail(existing, chunk, max_lines=MAX_UI_LOG_LINES):
    """Append log text and keep only the final N lines for browser state."""
    text = f"{existing or ''}{chunk or ''}"
    if not text:
        return ""
    lines = text.splitlines(keepends=True)
    if len(lines) <= max_lines:
        return text
    return "".join(lines[-max_lines:])


def format_runtime(seconds):
    """Format elapsed runtime."""
    total = max(0, int(round(float(seconds or 0))))
    hours, rem = divmod(total, 3600)
    mins, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h {mins:02d}m {secs:02d}s"
    if mins:
        return f"{mins}m {secs:02d}s"
    return f"{secs}s"


def update_active_job(job, updates):
    """Update stored metadata for the active process and return a job copy."""
    with COMPILE_LOCK:
        active = COMPILE_PROC.get("job")
        if active and job and active.get("pid") == job.get("pid"):
            active = dict(active)
            active.update(updates)
            COMPILE_PROC["job"] = active
            return dict(active)
    updated = dict(job or {})
    updated.update(updates)
    return updated


def start_compile_job(discovery, env_id, options):
    """Start one compile.py subprocess and return job metadata."""
    spec = build_compile_spec(discovery, env_id, options)
    log_file = tempfile.NamedTemporaryFile(delete=False, prefix="clubb_compile_", suffix=".log", dir="/tmp")
    log_path = log_file.name
    with COMPILE_LOCK:
        proc = COMPILE_PROC.get("proc")
        if proc is not None and proc.poll() is None:
            raise RuntimeError("A compile or rebuild job is already running.")
        popen_kwargs = {
            "cwd": REPO_ROOT,
            "env": spec["env"],
            "stdout": log_file,
            "stderr": subprocess.STDOUT,
            "text": True,
            "preexec_fn": os.setsid,
        }
        if spec.get("load_error"):
            raise RuntimeError(f"Cannot load selected module stack: {spec['load_error']}")
        if spec["shell"]:
            popen_kwargs["executable"] = "/bin/bash"
            proc = subprocess.Popen(spec["command"], shell=True, **popen_kwargs)
        else:
            proc = subprocess.Popen(spec["command"], **popen_kwargs)
        log_file.close()
        start_time = time.time()
        job = {
            "pid": proc.pid,
            "log": log_path,
            "start_time": start_time,
            "last_output_time": start_time,
            "command": spec["display"],
            "status": "running",
            "returncode": None,
            "kind": "compile",
            "environment": spec["environment"],
            "options": dict(options),
        }
        COMPILE_PROC["proc"] = proc
        COMPILE_PROC["job"] = job
        return dict(job)


def _compile_specs_from_builds(builds, discovery):
    """Build serialized compile.py specs for reconfiguring existing builds."""
    specs = []
    for build in builds or []:
        build_path = validated_build_path(build.get("path"))
        options = compile_options_from_build(build, discovery)
        spec = build_compile_spec(discovery, "current", options)
        if spec.get("load_error"):
            raise RuntimeError(f"Cannot load selected module stack: {spec['load_error']}")
        specs.append(
            {
                "build_path": build_path,
                "build_log_path": str(Path(build_path) / "cmake_build_output.txt"),
                "command": spec["command"],
                "display": spec["display"],
                "env": spec["env"],
            }
        )
    return specs


def start_rebuild_job(builds, discovery, label="selected builds"):
    """Start a sequential compile.py job for existing build configurations."""
    compile_specs = _compile_specs_from_builds(builds, discovery)
    if not compile_specs:
        raise RuntimeError("No build directories were selected for rebuild.")
    progress_file = tempfile.NamedTemporaryFile(delete=False, prefix="clubb_rebuild_progress_", suffix=".json", dir="/tmp")
    progress_path = progress_file.name
    progress_file.close()
    command = [sys.executable, "-u", "-c", _REBUILD_HELPER, json.dumps(compile_specs), progress_path]
    build_paths = [spec["build_path"] for spec in compile_specs]
    build_count = len(compile_specs)
    display_commands = [spec["display"] for spec in compile_specs]
    display = "\n".join(display_commands)
    log_file = tempfile.NamedTemporaryFile(delete=False, prefix="clubb_rebuild_", suffix=".log", dir="/tmp")
    log_path = log_file.name
    with COMPILE_LOCK:
        proc = COMPILE_PROC.get("proc")
        if proc is not None and proc.poll() is None:
            raise RuntimeError("A compile or rebuild job is already running.")
        proc = subprocess.Popen(
            command,
            cwd=REPO_ROOT,
            env=os.environ.copy(),
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
            preexec_fn=os.setsid,
        )
        log_file.close()
        start_time = time.time()
        job = {
            "pid": proc.pid,
            "log": log_path,
            "start_time": start_time,
            "last_output_time": start_time,
            "command": display,
            "status": "running",
            "returncode": None,
            "kind": "rebuild",
            "build_paths": build_paths,
            "progress": progress_path,
            "label": label,
            "build_count": build_count,
        }
        COMPILE_PROC["proc"] = proc
        COMPILE_PROC["job"] = job
        return dict(job)


def active_job():
    """Return current active job metadata if a process is still running."""
    with COMPILE_LOCK:
        proc = COMPILE_PROC.get("proc")
        job = COMPILE_PROC.get("job")
        if proc is None or job is None or proc.poll() is not None:
            return None
        return dict(job)


def pid_is_live(pid):
    """Return whether pid appears to still be running."""
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return False
    if pid <= 0:
        return False
    proc_stat = Path("/proc") / str(pid) / "stat"
    try:
        fields = proc_stat.read_text(encoding="utf-8", errors="replace").split()
        if len(fields) > 2 and fields[2] == "Z":
            return False
    except OSError:
        pass
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def job_process_is_live(job):
    """Return whether a stored job still has a running process."""
    if not job or job.get("status") != "running":
        return False
    if job.get("kind") == "rebuild" and rebuild_returncode_from_log(job.get("log")) is not None:
        return False
    with COMPILE_LOCK:
        proc = COMPILE_PROC.get("proc")
        active = COMPILE_PROC.get("job")
        if proc and active and active.get("pid") == job.get("pid"):
            return proc.poll() is None
    return pid_is_live(job.get("pid"))


def rebuild_returncode_from_log(log_path):
    """Recover a finished rebuild helper return code from its terminal log line."""
    if not log_path:
        return None
    try:
        text = Path(log_path).read_text(encoding="utf-8", errors="replace")[-8000:]
    except OSError:
        return None
    if "=== Rebuild complete ===" in text:
        return 0
    for line in reversed(text.splitlines()):
        if not line.startswith("=== Rebuild failed"):
            continue
        if " with exit " not in line:
            return 1
        exit_text = line.rsplit(" with exit ", 1)[1].split()[0]
        try:
            return int(exit_text)
        except ValueError:
            return 1
    return None


def read_rebuild_progress(job):
    """Return queued/current/completed paths for a running rebuild job."""
    if not job or job.get("kind") != "rebuild":
        return {"current_path": None, "queued_paths": [], "completed_paths": [], "state": ""}
    progress_path = job.get("progress")
    if progress_path:
        try:
            payload = json.loads(Path(progress_path).read_text(encoding="utf-8"))
            return {
                "current_path": payload.get("current_path"),
                "queued_paths": list(payload.get("queued_paths") or []),
                "completed_paths": list(payload.get("completed_paths") or []),
                "state": payload.get("state", ""),
            }
        except (OSError, json.JSONDecodeError):
            pass
    return {
        "current_path": None,
        "queued_paths": list(job.get("build_paths") or []),
        "completed_paths": [],
        "state": "queued",
    }


def rebuild_failed_path_from_log(log_path):
    """Recover the build path named by a failed rebuild helper log line."""
    if not log_path:
        return None
    try:
        text = Path(log_path).read_text(encoding="utf-8", errors="replace")[-8000:]
    except OSError:
        return None
    prefix = "=== Rebuild failed for "
    for line in reversed(text.splitlines()):
        if not line.startswith(prefix):
            continue
        failure_text = line[len(prefix):]
        if " with exit " in failure_text:
            failure_text = failure_text.rsplit(" with exit ", 1)[0]
        return failure_text.strip()
    return None


def poll_compile_job(job):
    """Return the latest process return code for a stored job."""
    if job and job.get("kind") == "rebuild":
        recovered = rebuild_returncode_from_log(job.get("log"))
        if recovered is not None:
            return recovered
    with COMPILE_LOCK:
        proc = COMPILE_PROC.get("proc")
        active = COMPILE_PROC.get("job")
        if proc and active and job and active.get("pid") == job.get("pid"):
            return proc.poll()
    if job and job.get("status") == "running" and not pid_is_live(job.get("pid")):
        return "lost"
    return None


def finish_compile_job(job, returncode):
    """Mark the global compile job finished."""
    with COMPILE_LOCK:
        active = COMPILE_PROC.get("job")
        if active and job and active.get("pid") == job.get("pid"):
            active = dict(active)
            active["status"] = "completed" if returncode == 0 else "failed"
            active["returncode"] = returncode
            COMPILE_PROC["job"] = active
            COMPILE_PROC["proc"] = None
            return active
    updated = dict(job or {})
    updated["status"] = "completed" if returncode == 0 else "failed"
    updated["returncode"] = returncode
    return updated


def cancel_compile_job(job):
    """Terminate the active compile subprocess."""
    with COMPILE_LOCK:
        proc = COMPILE_PROC.get("proc")
        active = COMPILE_PROC.get("job")
        if not proc or not active or not job or active.get("pid") != job.get("pid"):
            return False
        try:
            if proc.poll() is None:
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                except Exception:
                    proc.terminate()
                try:
                    proc.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                    except Exception:
                        proc.kill()
        finally:
            COMPILE_PROC["proc"] = None
            active = dict(active)
            active["status"] = "cancelled"
            active["returncode"] = 1
            COMPILE_PROC["job"] = active
        return True
