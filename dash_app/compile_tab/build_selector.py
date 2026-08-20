"""Compile-owned selected-build triggers and the shared floating selector."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from dash import ALL, Input, Output, dcc, html


REPO_ROOT = Path(__file__).resolve().parents[2]
BUILD_SELECTOR_TRIGGER_IDS = (
    "run-selected-build-badge",
    "tune-selected-build-badge",
    "profile-selected-build-badge",
)
RUN_IMPLEMENTATIONS = ("fortran", "python", "jax")


def _cmake_cache(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return values
    for raw_line in lines:
        line = raw_line.strip()
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        left, value = line.split("=", 1)
        values[left.split(":", 1)[0]] = value
    return values


def selected_build_info(repo_root: Path = REPO_ROOT) -> dict[str, Any]:
    """Resolve the default install exactly as run_scm.py does."""
    root = Path(repo_root)
    install_root = root / "install"
    selected = install_root / "selected"
    latest = install_root / "latest"
    if os.path.lexists(selected):
        alias = selected
        source = "selected"
    elif os.path.lexists(latest):
        alias = latest
        source = "latest"
    else:
        return {
            "name": "No build selected",
            "status": "missing",
            "source": "none",
            "target": "",
            "title": "No default CLUBB install was found. Compile and select a build before running CLUBB.",
        }

    try:
        target = alias.resolve()
    except OSError:
        target = alias
    name = target.name or alias.name
    target_exists = target.is_dir()
    build_dir = root / "build" / name
    cache = _cmake_cache(build_dir / "CMakeCache.txt")
    compiler_path = cache.get("CMAKE_Fortran_COMPILER", "")
    compiler_name = Path(compiler_path).name if compiler_path else "unknown"
    status = "ready" if target_exists else "broken"
    if status == "ready" and source == "latest":
        status = "fallback"

    lines = [
        f"Default CLUBB build: {name}",
        f"Selection: install/{source} → {target}",
        f"CMake build directory: {build_dir}",
        f"Fortran compiler: {compiler_name}",
    ]
    if compiler_path:
        lines.append(f"Compiler path: {compiler_path}")
    for label, key in (
        ("Build type", "CMAKE_BUILD_TYPE"),
        ("Precision", "PRECISION"),
        ("GPU backend", "GPU"),
        ("OpenMP", "ENABLE_OMP"),
        ("GPTL", "USE_GPTL"),
    ):
        if cache.get(key):
            lines.append(f"{label}: {cache[key]}")
    if source == "latest":
        lines.append("install/selected is not set, so run_scm.py will use install/latest.")
    if not target_exists:
        lines.append("The selected install target is missing; runs using the default will fail.")
    lines.append("Explicit -exe or -install_dir runner options override this default.")
    return {
        "name": name,
        "status": status,
        "source": source,
        "target": str(target),
        "title": "\n".join(lines),
    }


def normalize_run_implementation(value: Any) -> str:
    implementation = str(value or "fortran").strip().lower()
    return implementation if implementation in RUN_IMPLEMENTATIONS else "fortran"


def build_implementation_capability(
    install_dir: str | Path,
    implementation: str,
    repo_root: Path = REPO_ROOT,
) -> tuple[bool, str]:
    """Return whether an install can launch one implementation today."""
    install = Path(install_dir)
    implementation = normalize_run_implementation(implementation)
    if implementation == "fortran":
        available = (install / "clubb_standalone").is_file()
        return available, "" if available else "Installed clubb_standalone is missing"

    runtime = install / "python"
    extension = any(runtime.glob("clubb_f2py*.so"))
    backend = any(
        (runtime / name).is_file()
        for name in (
            "libclubb_f2py_backend.so",
            "libclubb_f2py_backend.dylib",
            "clubb_f2py_backend.dll",
        )
    )
    python_package = (runtime / "clubb_python").is_dir()
    if not (runtime.is_dir() and extension and backend and python_package):
        return False, "Python runtime is incomplete; rebuild with Python API enabled"
    if implementation == "jax" and not (Path(repo_root) / "clubb_jax" / "clubb_standalone.py").is_file():
        return False, "JAX standalone driver is missing"
    return True, ""


def selected_launch_target(implementation: Any, repo_root: Path = REPO_ROOT) -> dict[str, str]:
    """Resolve and validate the implementation/install pair at submission time."""
    implementation = normalize_run_implementation(implementation)
    info = selected_build_info(repo_root)
    install_dir = str(info.get("target") or "")
    if not install_dir:
        raise ValueError("select a CLUBB build before running")
    available, reason = build_implementation_capability(install_dir, implementation, repo_root)
    if not available:
        raise ValueError(f"{implementation.title()} cannot use {info['name']}: {reason}")
    return {
        "implementation": implementation,
        "install_dir": install_dir,
        "build_name": str(info["name"]),
    }


def build_selector_trigger(component_id: str):
    """Render one lightweight request button for the shared selector."""
    info = selected_build_info()
    return html.Button(
        build_selector_trigger_children(info["name"], "fortran"),
        id=component_id,
        type="button",
        n_clicks=0,
        className=f"selected-build-badge selected-build-badge-{info['status']}",
        title="Choose or rebuild a CLUBB build",
        **{"aria-haspopup": "dialog"},
    )


def build_selector_trigger_children(name: str, implementation: Any = "fortran"):
    implementation = normalize_run_implementation(implementation)
    return [
        html.Span(implementation.upper(), className="selected-build-badge-label"),
        html.Code(name, className="selected-build-badge-name"),
        html.Span("▾", className="selected-build-badge-chevron", **{"aria-hidden": "true"}),
    ]


def build_selector_overlay():
    """Render the one shared popover used by every execution tab."""
    return html.Div(
        [
            dcc.Store(id="compile-build-selector-anchor", data=None),
            dcc.Store(
                id="compile-run-implementation",
                data="fortran",
                storage_type="local",
            ),
            html.Div(
                [
                    html.Button(
                        "",
                        id="compile-build-selector-dismiss",
                        type="button",
                        n_clicks=0,
                        className="compile-build-selector-backdrop",
                        **{"aria-label": "Close build selector"},
                    ),
                    html.Div(
                        id="compile-build-selector-menu",
                        className="compile-build-selector-menu",
                        role="dialog",
                        **{"aria-label": "Choose or rebuild a CLUBB build"},
                    ),
                ],
                id="compile-build-selector-popover",
                className="compile-build-selector-popover compile-build-selector-popover-closed",
            ),
        ],
        id="compile-build-selector-root",
    )


def selector_popover_style(anchor: dict[str, Any] | None) -> dict[str, Any]:
    """Place the fixed popover beside its requesting trigger."""
    if not anchor:
        return {}
    viewport_width = max(320.0, float(anchor.get("viewport_width") or 0))
    viewport_height = max(240.0, float(anchor.get("viewport_height") or 0))
    width = min(520.0, viewport_width - 16.0)
    left = min(max(8.0, float(anchor.get("left") or 0)), viewport_width - width - 8.0)
    top = float(anchor.get("top") or 0)
    bottom = float(anchor.get("bottom") or 0)
    space_below = viewport_height - bottom - 8.0
    space_above = top - 8.0
    style: dict[str, Any] = {"left": f"{left}px", "width": f"{width}px"}
    if space_below < 220.0 and space_above > space_below:
        style.update({"bottom": f"{viewport_height - top + 6.0}px", "maxHeight": f"{max(120.0, space_above - 6.0)}px"})
    else:
        style.update({"top": f"{bottom + 6.0}px", "maxHeight": f"{max(120.0, space_below - 6.0)}px"})
    return style


def register_build_selector_position_callback(app) -> None:
    """Capture trigger coordinates in the browser and close on outside/select."""
    app.clientside_callback(
        """
        function(runClicks, tuneClicks, profileClicks, dismissClicks, selectClicks) {
            const ctx = window.dash_clientside.callback_context || {};
            const triggered = ctx.triggered_id;
            const event = (ctx.triggered || [])[0] || {};
            if (!triggered || !(Number(event.value) > 0)) {
                return window.dash_clientside.no_update;
            }
            if (triggered === "compile-build-selector-dismiss" ||
                (typeof triggered === "object" && triggered.type === "compile-build-selector-select")) {
                return null;
            }
            if (typeof triggered !== "string") {
                return window.dash_clientside.no_update;
            }
            const button = document.getElementById(triggered);
            if (!button) {
                return window.dash_clientside.no_update;
            }
            const rect = button.getBoundingClientRect();
            return {
                trigger_id: triggered,
                left: rect.left,
                top: rect.top,
                bottom: rect.bottom,
                viewport_width: window.innerWidth,
                viewport_height: window.innerHeight
            };
        }
        """,
        Output("compile-build-selector-anchor", "data"),
        Input(BUILD_SELECTOR_TRIGGER_IDS[0], "n_clicks"),
        Input(BUILD_SELECTOR_TRIGGER_IDS[1], "n_clicks"),
        Input(BUILD_SELECTOR_TRIGGER_IDS[2], "n_clicks"),
        Input("compile-build-selector-dismiss", "n_clicks"),
        Input({"type": "compile-build-selector-select", "index": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
