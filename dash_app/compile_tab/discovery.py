"""Environment, toolchain, and existing-build discovery for the compile tab."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import shutil
from pathlib import Path

from .state import (
    BUILD_DIR,
    COMPILER_COMMANDS,
    COMPILER_NAME_MAP,
    INSTALL_DIR,
    REPO_ROOT,
    TOOLCHAIN_DIR,
)

COMPILER_HINTS = {
    "gcc": ("gcc", "gnu", "gfortran"),
    "intel": ("intel", "oneapi", "ifx", "ifort"),
    "nvhpc": ("nvhpc", "nvidia", "nvfortran", "pgi"),
    "cce": ("cce", "cray", "crayftn"),
}

COMPILER_MODULE_BASES = {
    "gcc",
    "gfortran",
    "intel",
    "intel-oneapi",
    "intel-oneapi-compilers",
    "nvhpc",
    "nvfortran",
    "oneapi",
    "cce",
    "cray",
    "crayftn",
}

NETCDF_HINTS = ("netcdf", "netcdf-fortran", "netcdf_fortran")
LMOD_HELPER_TIMEOUT = 30

COMPILER_SUFFIX_ALIASES = {
    "gcc": ("gcc", "gnu", "gfortran"),
    "intel": ("intel", "oneapi", "intel-oneapi"),
    "nvhpc": ("nvhpc", "nvidia", "pgi"),
    "cce": ("cce", "cray", "crayftn"),
}

COMPILER_SUFFIX_MAP = {
    alias: canonical
    for canonical, aliases in COMPILER_SUFFIX_ALIASES.items()
    for alias in aliases
}

_LMOD_LOAD_HELPER = r"""
import json
import os
import subprocess
import sys

modules = json.loads(sys.argv[1])
lmod_cmd = os.environ.get("LMOD_CMD")
if not lmod_cmd:
    print(json.dumps({"returncode": 1, "stderr": "LMOD_CMD is not set.", "env": dict(os.environ)}))
    raise SystemExit(0)

proc = subprocess.run(
    [lmod_cmd, "python", "load", *modules],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
)
exec_error = ""
if proc.returncode == 0 and proc.stdout:
    try:
        exec(proc.stdout, {"os": os})
    except Exception as exc:
        exec_error = str(exc)
        proc = subprocess.CompletedProcess(proc.args, 1, proc.stdout, proc.stderr)

print(json.dumps({
    "returncode": proc.returncode,
    "stderr": (proc.stderr or "") + (("\n" + exec_error) if exec_error else ""),
    "env": dict(os.environ),
}))
"""


def canonical_compiler(name):
    """Map compiler executable or module-family names to CLUBB toolchain names."""
    return COMPILER_NAME_MAP.get((name or "").lower(), (name or "").lower())


def canonical_compiler_from_path(path):
    """Map a compiler executable path to a CLUBB toolchain compiler name."""
    return canonical_compiler(os.path.basename(path or ""))


def compiler_from_env(env):
    """Infer the canonical compiler from an environment."""
    lmod_family = (env or {}).get("LMOD_FAMILY_COMPILER", "")
    if lmod_family:
        return canonical_compiler(lmod_family)
    fc = (env or {}).get("FC", "")
    if fc:
        return canonical_compiler_from_path(fc)
    return ""


def command_in_env(command, env):
    """Return the executable path for command under env's PATH."""
    return shutil.which(command, path=(env or {}).get("PATH"))


def current_platform_key():
    """Return the platform prefix used by default toolchain filenames."""
    uname = os.uname()
    return f"{uname.sysname.lower()}_{uname.machine}"


def parse_toolchain(path):
    """Return display metadata for one CMake toolchain file."""
    file_path = Path(path)
    stem = file_path.stem
    parts = stem.split("_")
    compiler = parts[-1] if parts else stem
    platform = stem[: -(len(compiler) + 1)] if len(stem) > len(compiler) else ""
    return {
        "name": stem,
        "path": str(file_path),
        "platform": platform,
        "compiler": compiler,
        "matches_host": platform == current_platform_key(),
    }


def discover_toolchains():
    """Discover checked-in toolchain files."""
    toolchain_dir = Path(TOOLCHAIN_DIR)
    if not toolchain_dir.exists():
        return []
    return [
        parse_toolchain(path)
        for path in sorted(toolchain_dir.glob("*.cmake"))
    ]


def discover_native_compilers(toolchains):
    """Discover compiler executables available directly on PATH."""
    host_key = current_platform_key()
    matching_compilers = {
        item["compiler"]
        for item in toolchains
        if item.get("platform") == host_key
    }
    environments = [
        {
            "id": "current",
            "kind": "current",
            "label": "Current environment",
            "detail": "Use FC/LMOD_FAMILY_COMPILER if set, otherwise compile.py falls back to gfortran.",
            "canonical": canonical_compiler(os.environ.get("LMOD_FAMILY_COMPILER") or os.environ.get("FC") or "gfortran"),
            "compiler_path": "",
            "enabled": True,
        }
    ]

    seen = set()
    for command in COMPILER_COMMANDS:
        path = shutil.which(command)
        if not path:
            continue
        canonical = canonical_compiler(command)
        key = (canonical, path)
        if key in seen:
            continue
        seen.add(key)
        has_toolchain = canonical in matching_compilers
        environments.append(
            {
                "id": f"native:{command}:{path}",
                "kind": "native",
                "label": f"{command} ({canonical})",
                "detail": path,
                "canonical": canonical,
                "compiler_path": path,
                "enabled": has_toolchain,
                "disabled_reason": "" if has_toolchain else f"No {host_key}_{canonical}.cmake toolchain found.",
            }
        )
    return environments


def find_lmod_command():
    """Return the Lmod command path if this process has Lmod initialized."""
    lmod_cmd = os.environ.get("LMOD_CMD")
    if lmod_cmd and os.path.exists(lmod_cmd):
        return lmod_cmd
    return shutil.which("lmod") or ""


def module_support_available():
    """Return whether Lmod is initialized enough for direct Python use."""
    return bool(find_lmod_command() and os.environ.get("MODULEPATH"))


def parse_lmod_module_entries(text):
    """Parse terse Lmod avail/spider output into module entries."""
    entries = []
    seen = set()
    section_path = ""
    for raw_line in (text or "").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("-"):
            continue
        if line.startswith("/") and line.endswith(":"):
            section_path = line[:-1]
            continue
        if line.startswith("/"):
            continue
        if "=" in line or line.startswith(("export ", "echo ", "unset ")):
            continue
        token = line.split()[0].strip("'\";")
        if not token or token.endswith(":") or token.endswith("/"):
            continue
        if token in seen:
            continue
        seen.add(token)
        entries.append({"module": token, "path": section_path})
    return entries


def list_lmod_available_module_entries(env=None):
    """Return module entries visible to Lmod in env."""
    base_env = dict(env or os.environ)
    lmod_cmd = base_env.get("LMOD_CMD") or find_lmod_command()
    if not lmod_cmd:
        return []
    try:
        proc = subprocess.run(
            [lmod_cmd, "--terse", "avail"],
            env=base_env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=LMOD_HELPER_TIMEOUT,
        )
    except (OSError, subprocess.TimeoutExpired):
        return []
    return parse_lmod_module_entries(f"{proc.stderr}\n{proc.stdout}")


def list_lmod_available_modules(env=None):
    """Return modules visible to Lmod in env."""
    return [
        entry["module"]
        for entry in list_lmod_available_module_entries(env)
    ]


def resolve_lmod_stack(modules, base_env=None):
    """Load modules in an isolated helper and return the resulting environment."""
    cleaned_modules = [module for module in modules or [] if module]
    env = dict(base_env or os.environ)
    lmod_cmd = env.get("LMOD_CMD") or find_lmod_command()
    if not (lmod_cmd and env.get("MODULEPATH")):
        return {
            "ok": False,
            "returncode": 1,
            "stderr": "Lmod is not initialized in this process.",
            "env": env,
            "modules": cleaned_modules,
        }
    env["LMOD_CMD"] = lmod_cmd
    try:
        proc = subprocess.run(
            [sys.executable, "-B", "-c", _LMOD_LOAD_HELPER, json.dumps(cleaned_modules)],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=LMOD_HELPER_TIMEOUT,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {
            "ok": False,
            "returncode": 1,
            "stderr": str(exc),
            "env": env,
            "modules": cleaned_modules,
        }
    try:
        payload = json.loads(proc.stdout)
    except json.JSONDecodeError:
        return {
            "ok": False,
            "returncode": proc.returncode or 1,
            "stderr": (proc.stderr or proc.stdout or "Lmod helper did not return JSON.").strip(),
            "env": env,
            "modules": cleaned_modules,
        }
    stderr = "\n".join(part for part in [payload.get("stderr", ""), proc.stderr] if part).strip()
    return {
        "ok": payload.get("returncode") == 0,
        "returncode": payload.get("returncode", proc.returncode),
        "stderr": stderr,
        "env": payload.get("env", env),
        "modules": cleaned_modules,
    }


def module_base_name(module_name):
    """Return a normalized package name for a module entry."""
    return (module_name or "").split("/", 1)[0].lower().replace("_", "-")


def module_matches_hints(module_name, hints):
    """Return whether module name looks relevant for one of the hint strings."""
    normalized = module_base_name(module_name)
    return any(hint in normalized for hint in hints)


def module_is_compiler_like(module_name):
    """Return whether a module name looks like a compiler-family module."""
    return module_base_name(module_name) in COMPILER_MODULE_BASES


def compiler_hint_rank(module_name):
    """Rank module names that are likely to provide a Fortran compiler."""
    normalized = module_base_name(module_name)
    for rank, hints in enumerate(COMPILER_HINTS.values()):
        if any(hint in normalized for hint in hints):
            return rank
    return len(COMPILER_HINTS)


def summarize_lmod_env(env):
    """Return compact diagnostics for a loaded module environment."""
    return {
        "LMOD_FAMILY_COMPILER": env.get("LMOD_FAMILY_COMPILER", ""),
        "LOADEDMODULES": env.get("LOADEDMODULES", ""),
        "FC": env.get("FC", ""),
        "CC": env.get("CC", ""),
        "cmake": command_in_env("cmake", env) or "",
        "nf-config": command_in_env("nf-config", env) or "",
        "nc-config": command_in_env("nc-config", env) or "",
    }


def module_version(module_name):
    """Return the version portion of a module name, without build suffixes."""
    if "/" not in (module_name or ""):
        return ""
    version = module_name.split("/", 1)[1]
    return version.split("-", 1)[0]


def active_compiler_version(modules, env):
    """Infer the selected compiler version from the loaded compiler module."""
    canonical = compiler_from_env(env)
    if not canonical:
        return ""
    aliases = COMPILER_SUFFIX_ALIASES.get(canonical, ())
    for module in modules or []:
        base = module_base_name(module)
        if base in aliases or any(alias in base for alias in aliases):
            return module_version(module)
    return ""


def compiler_suffix_for_module(module_name):
    """Return the explicit compiler suffix metadata embedded in a module name."""
    if "/" not in (module_name or ""):
        return None
    version_text = module_name.split("/", 1)[1]
    parts = version_text.split("-")
    for index, part in enumerate(parts[:-1]):
        canonical = COMPILER_SUFFIX_MAP.get(part)
        if not canonical:
            continue
        package_version = "-".join(parts[:index])
        suffix_version = "-".join(parts[index + 1:])
        if package_version and suffix_version:
            return {
                "compiler": canonical,
                "version": suffix_version,
                "package_version": package_version,
            }
    return None


def module_alias_key(module_name):
    """Return a de-duplication key that ignores matching compiler suffix aliases."""
    if "/" not in (module_name or ""):
        return (module_name or "").lower()
    base, version_text = module_name.split("/", 1)
    suffix = compiler_suffix_for_module(module_name)
    version = suffix.get("package_version") if suffix else version_text
    return f"{base.lower()}/{version.lower()}"


def dedupe_module_names(module_names):
    """Collapse compiler-suffixed aliases, preferring the cleaner module name."""
    best_by_key = {}
    order_by_key = {}
    for index, module in enumerate(module_names or []):
        key = module_alias_key(module)
        suffix = compiler_suffix_for_module(module)
        score = (1 if suffix else 0, index)
        if key not in best_by_key or score < best_by_key[key][0]:
            best_by_key[key] = (score, module)
            order_by_key.setdefault(key, index)
    return [
        best_by_key[key][1]
        for key in sorted(order_by_key, key=order_by_key.get)
    ]


def module_matches_active_compiler(module_name, modules, env):
    """Return whether a module is compatible with the loaded compiler stack."""
    canonical = compiler_from_env(env)
    if not canonical:
        return True
    suffix = compiler_suffix_for_module(module_name)
    if not suffix:
        return True
    if suffix["compiler"] != canonical:
        return False
    version = active_compiler_version(modules, env)
    return not version or suffix["version"] == version


def filter_modules_for_active_compiler(module_names, modules, env):
    """Drop compiler-suffixed modules that conflict with the active compiler."""
    return dedupe_module_names([
        module for module in module_names
        if module_matches_active_compiler(module, modules, env)
    ])


def module_path_entries(env):
    """Return the current env's module paths as a set."""
    return {
        path for path in (env or {}).get("MODULEPATH", "").split(os.pathsep)
        if path
    }


def filter_module_entries_for_active_compiler(entries, modules, env):
    """Keep modules from the active compiler tree or matching compiler suffix."""
    base_paths = module_path_entries(os.environ)
    active_paths = module_path_entries(env) - base_paths
    if not active_paths:
        return filter_modules_for_active_compiler(
            [entry["module"] for entry in entries],
            modules,
            env,
        )

    filtered = []
    for entry in entries:
        module = entry["module"]
        suffix = compiler_suffix_for_module(module)
        if entry.get("path") in active_paths:
            if module_matches_active_compiler(module, modules, env):
                filtered.append(module)
        elif suffix and module_matches_active_compiler(module, modules, env):
            filtered.append(module)
    return dedupe_module_names(filtered)


def discover_lmod_state(toolchains):
    """Discover Lmod module candidates for a hierarchical module-stack UI."""
    lmod_cmd = find_lmod_command()
    if not module_support_available():
        return {
            "available": False,
            "lmod_cmd": lmod_cmd,
            "message": "Lmod is not initialized. Native compiler detection is available.",
            "base_modules": [],
            "compiler_modules": [],
        }

    host_key = current_platform_key()
    matching_compilers = {
        item["compiler"]
        for item in toolchains
        if item.get("platform") == host_key
    }
    base_modules = dedupe_module_names(list_lmod_available_modules())
    likely_modules = [
        module for module in base_modules
        if module_is_compiler_like(module)
    ]
    compiler_modules = []
    seen = set()
    for module in likely_modules[:80]:
        result = resolve_lmod_stack([module])
        env = result.get("env", {})
        canonical = compiler_from_env(env)
        enabled = bool(result.get("ok") and canonical in matching_compilers)
        if module in seen:
            continue
        seen.add(module)
        compiler_modules.append(
            {
                "module": module,
                "label": module,
                "canonical": canonical,
                "enabled": enabled,
                "rank": compiler_hint_rank(module),
                "diagnostics": summarize_lmod_env(env),
                "disabled_reason": "" if enabled else (
                    result.get("stderr")
                    or ("No matching CLUBB toolchain was inferred." if canonical else "Module did not expose a Fortran compiler.")
                ),
            }
        )

    compiler_modules.sort(
        key=lambda item: (
            item.get("rank", 99),
            not item.get("enabled", False),
            item.get("module", ""),
        )
    )
    return {
        "available": True,
        "lmod_cmd": lmod_cmd,
        "message": "Lmod detected. Select a compiler module, then exposed library modules.",
        "base_modules": base_modules,
        "compiler_modules": compiler_modules,
    }


def available_modules_after_stack(modules):
    """Return currently visible modules after loading a partial module stack."""
    result = resolve_lmod_stack(modules)
    if not result.get("ok"):
        return []
    env = result.get("env", {})
    return filter_module_entries_for_active_compiler(
        list_lmod_available_module_entries(env),
        modules,
        env,
    )


def netcdf_modules_after_stack(modules):
    """Return NetCDF-like modules visible after loading a partial stack."""
    return [
        module for module in available_modules_after_stack(modules)
        if module_matches_hints(module, NETCDF_HINTS)
    ]


def parse_cmake_cache(cache_path):
    """Parse a CMakeCache.txt file into a small key/value mapping."""
    values = {}
    try:
        for raw_line in Path(cache_path).read_text(encoding="utf-8", errors="replace").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("//") or line.startswith("#") or "=" not in line:
                continue
            left, value = line.split("=", 1)
            key = left.split(":", 1)[0]
            values[key] = value
    except OSError:
        return values
    return values


def discover_existing_builds():
    """Discover existing CMake build directories and their important options."""
    build_root = Path(BUILD_DIR)
    install_root = Path(INSTALL_DIR)
    latest_target = ""
    selected_target = ""
    latest_path = install_root / "latest"
    selected_path = install_root / "selected"
    try:
        if latest_path.exists() or latest_path.is_symlink():
            latest_target = str(latest_path.resolve())
    except OSError:
        latest_target = ""
    try:
        if selected_path.exists() or selected_path.is_symlink():
            selected_target = str(selected_path.resolve())
    except OSError:
        selected_target = ""

    if not build_root.exists():
        return []

    builds = []
    for build_dir in sorted(path for path in build_root.iterdir() if path.is_dir()):
        cache_path = build_dir / "CMakeCache.txt"
        if not cache_path.exists():
            continue
        cache = parse_cmake_cache(cache_path)
        install_prefix = cache.get("CMAKE_INSTALL_PREFIX", "")
        install_exists = bool(install_prefix and Path(install_prefix).exists())
        build_log = build_dir / "cmake_build_output.txt"
        builds.append(
            {
                "name": build_dir.name,
                "path": str(build_dir),
                "install_prefix": install_prefix,
                "install_exists": install_exists,
                "is_latest": bool(install_prefix and latest_target and str(Path(install_prefix).resolve()) == latest_target),
                "is_selected": bool(install_prefix and selected_target and str(Path(install_prefix).resolve()) == selected_target),
                "has_log": build_log.exists(),
                "build_type": cache.get("CMAKE_BUILD_TYPE", ""),
                "toolchain": cache.get("CMAKE_TOOLCHAIN_FILE", ""),
                "precision": cache.get("PRECISION", ""),
                "gpu": cache.get("GPU", ""),
                "netcdf": cache.get("USE_NetCDF", ""),
                "openmp": cache.get("ENABLE_OMP", ""),
                "python": cache.get("ENABLE_F2PY", ""),
                "silhs": cache.get("SILHS", ""),
                "modified": int(build_dir.stat().st_mtime),
            }
        )
    builds.sort(key=lambda item: item.get("modified", 0), reverse=True)
    return builds


def discover_compile_state():
    """Collect all static-ish data needed by the compile tab."""
    toolchains = discover_toolchains()
    environments = discover_native_compilers(toolchains)
    lmod = discover_lmod_state(toolchains)
    return {
        "repo_root": REPO_ROOT,
        "host_platform": current_platform_key(),
        "toolchains": toolchains,
        "environments": environments,
        "lmod": lmod,
        "builds": discover_existing_builds(),
        "env": {
            "FC": os.environ.get("FC", ""),
            "CC": os.environ.get("CC", ""),
            "LMOD_FAMILY_COMPILER": os.environ.get("LMOD_FAMILY_COMPILER", ""),
            "LOADEDMODULES": os.environ.get("LOADEDMODULES", ""),
            "MODULEPATH": os.environ.get("MODULEPATH", ""),
        },
        "tools": {
            "cmake": shutil.which("cmake") or "",
            "ninja": shutil.which("ninja") or "",
            "python": shutil.which("python") or "",
        },
        "module_support": lmod.get("available", False),
    }
