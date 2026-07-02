"""Compile subprocess launching, cancellation, and log helpers."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shlex
import signal
import subprocess
import sys
import tempfile
import time

from .discovery import command_in_env, compiler_from_env, resolve_lmod_stack
from .state import BUILD_DIR, COMPILE_LOCK, COMPILE_PROC, MAX_UI_LOG_LINES, REPO_ROOT


BUILD_STATUS_TIMEOUT = 8

_REBUILD_HELPER = r"""
import json
import os
import subprocess
import sys

build_paths = json.loads(sys.argv[1])
parallel_jobs = sys.argv[2]


def display_command(cmd):
    return "+ " + " ".join(cmd)


def write_line(build_log, line):
    print(line, flush=True)
    build_log.write(line + "\n")
    build_log.flush()


def stream_command(cmd, build_log):
    write_line(build_log, display_command(cmd))
    proc = subprocess.Popen(
        cmd,
        cwd=os.getcwd(),
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


for index, build_path in enumerate(build_paths, start=1):
    cmd = ["cmake", "--build", build_path, "--target", "install", "--parallel", parallel_jobs]
    build_log_path = os.path.join(build_path, "cmake_build_output.txt")
    header_lines = [
        f"=== [{index}/{len(build_paths)}] Rebuilding {build_path} ===",
        f"=== Writing build log: {build_log_path} ===",
    ]
    returncode = 0
    try:
        with open(build_log_path, "a", encoding="utf-8", errors="replace") as build_log:
            build_log.write("\n" + "\n".join(header_lines) + "\n")
            build_log.flush()
            for line in header_lines:
                print(line, flush=True)
            returncode = stream_command(cmd, build_log)
            if returncode == 0:
                write_line(build_log, f"=== Rebuild complete for {build_path} ===")
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

print("=== Rebuild complete ===", flush=True)
"""


def _bool_flag(options, key, flag):
    return [flag] if options.get(key) else []


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
    precision = options.get("precision") or "double"
    gpu = options.get("gpu") or "none"
    if precision != "double":
        argv.extend(["-precision", precision])
    if gpu != "none":
        argv.extend(["-gpu", gpu])
    toolchain = options.get("toolchain") or "auto"
    if toolchain != "auto":
        argv.extend(["-toolchain", toolchain])
    argv.extend(_bool_flag(options, "debug", "-debug"))
    argv.extend(_bool_flag(options, "run_tests", "-run_tests"))
    argv.extend(_bool_flag(options, "python", "-python"))
    argv.extend(_bool_flag(options, "fresh", "-fresh"))
    argv.extend(_bool_flag(options, "disable_netcdf", "-disable_netcdf"))
    argv.extend(_bool_flag(options, "disable_silhs", "-disable_silhs"))
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


def detect_build_status(build):
    """Return whether the build system would do work for a build directory."""
    build_path = build.get("path") if isinstance(build, dict) else build
    try:
        path = validated_build_path(build_path)
    except RuntimeError as exc:
        return {
            "status": "unknown",
            "label": "unknown",
            "detail": str(exc),
            "output": "",
            "checked_at": time.time(),
        }

    command = ["cmake", "--build", path, "--", "-n"]
    try:
        proc = subprocess.run(
            command,
            cwd=REPO_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
            timeout=BUILD_STATUS_TIMEOUT,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {
            "status": "unknown",
            "label": "unknown",
            "detail": f"Dry-run timed out after {BUILD_STATUS_TIMEOUT}s.",
            "output": "",
            "checked_at": time.time(),
        }
    except OSError as exc:
        return {
            "status": "unknown",
            "label": "unknown",
            "detail": str(exc),
            "output": "",
            "checked_at": time.time(),
        }

    output = proc.stdout or ""
    lower_output = output.lower()
    detail_lines = [line.strip() for line in output.splitlines() if line.strip()]
    detail = detail_lines[-1] if detail_lines else ""
    if proc.returncode != 0:
        status = "unknown"
        label = "unknown"
        detail = detail or f"Dry-run failed with exit {proc.returncode}."
    elif "no work to do" in lower_output or "nothing to be done" in lower_output:
        status = "current"
        label = "current"
        detail = "No build work pending."
    elif "re-running cmake" in lower_output:
        status = "needs_configure"
        label = "needs configure"
    elif "pgcuda" in lower_output:
        status = "toolchain_dirty"
        label = "toolchain dirty"
    else:
        status = "needs_rebuild"
        label = "needs rebuild"

    return {
        "status": status,
        "label": label,
        "detail": detail,
        "output": output[-4000:],
        "command": " ".join(shlex.quote(part) for part in command),
        "returncode": proc.returncode,
        "checked_at": time.time(),
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
    if options.get("python") and options.get("disable_netcdf"):
        warnings.append("-python requires NetCDF; compile.py will reject this combination.")
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
        if not options.get("disable_netcdf") and not (diagnostics.get("nf-config") or diagnostics.get("nc-config")):
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


def start_rebuild_job(build_paths, label="selected builds"):
    """Start a sequential rebuild/install job for existing CMake build dirs."""
    safe_paths = [validated_build_path(path) for path in build_paths or []]
    if not safe_paths:
        raise RuntimeError("No build directories were selected for rebuild.")
    nproc = str(os.cpu_count() or 4)
    command = [sys.executable, "-u", "-c", _REBUILD_HELPER, json.dumps(safe_paths), nproc]
    build_count = len(safe_paths)
    display_commands = [
        " ".join(shlex.quote(part) for part in ["cmake", "--build", path, "--target", "install", "--parallel", nproc])
        for path in safe_paths
    ]
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
            "build_paths": safe_paths,
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
