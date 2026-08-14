"""Read-only selected CLUBB build badges shared by execution tabs."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from dash import Input, Output, html


REPO_ROOT = Path(__file__).resolve().parents[2]
BADGE_IDS = (
    "run-selected-build-badge",
    "tune-selected-build-badge",
    "profile-selected-build-badge",
)


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
            "title": (
                "No default CLUBB install was found.\n"
                f"Expected {selected} or fallback {latest}.\n"
                "Compile and select a build before running CLUBB."
            ),
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
        "title": "\n".join(lines),
    }


def selected_build_badge_properties(repo_root: Path = REPO_ROOT):
    info = selected_build_info(repo_root)
    children = [
        html.Span("BUILD", className="selected-build-badge-label"),
        html.Code(info["name"], className="selected-build-badge-name"),
    ]
    class_name = f"selected-build-badge selected-build-badge-{info['status']}"
    return children, info["title"], class_name


def selected_build_badge(component_id: str):
    children, title, class_name = selected_build_badge_properties()
    return html.Div(
        children,
        id=component_id,
        className=class_name,
        title=title,
        **{"aria-label": title},
    )


def register_selected_build_callback(app) -> None:
    """Keep every execution-tab badge synchronized with Compile selection."""
    outputs = []
    for component_id in BADGE_IDS:
        outputs.extend(
            (
                Output(component_id, "children"),
                Output(component_id, "title"),
                Output(component_id, "className"),
            )
        )

    @app.callback(
        *outputs,
        Input("selected-build-refresh", "n_intervals"),
        Input("compile-discovery", "data"),
    )
    def refresh_selected_build(_ticks, _compile_discovery):
        properties = selected_build_badge_properties()
        return tuple(value for _component_id in BADGE_IDS for value in properties)
