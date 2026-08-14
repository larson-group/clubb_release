from pathlib import Path

from dash import Dash, dcc, html
from dash.development.base_component import Component

from dash_app.run_tab.layout import build_run_action_section
from dash_app.shared.selected_build import (
    BADGE_IDS,
    register_selected_build_callback,
    selected_build_badge,
    selected_build_info,
)
from dash_app.tune_tab.layout import build_top_controls


def component_ids(component) -> set[str]:
    if not isinstance(component, Component):
        return set()
    found = {component.id} if isinstance(getattr(component, "id", None), str) else set()
    children = getattr(component, "children", None)
    if children is None:
        return found
    for child in children if isinstance(children, (list, tuple)) else [children]:
        found.update(component_ids(child))
    return found


def write_cache(path: Path) -> None:
    path.parent.mkdir(parents=True)
    path.write_text(
        "\n".join(
            (
                "CMAKE_BUILD_TYPE:STRING=Release",
                "CMAKE_Fortran_COMPILER:FILEPATH=/opt/nvidia/bin/nvfortran",
                "PRECISION:STRING=double",
                "GPU:STRING=openacc",
                "ENABLE_OMP:BOOL=OFF",
                "USE_GPTL:BOOL=ON",
            )
        ),
        encoding="utf-8",
    )


def test_selected_build_matches_run_scm_alias_order_and_reports_cmake_details(tmp_path):
    name = "nvhpc_GPUopenacc_PRECdouble"
    install = tmp_path / "install" / name
    install.mkdir(parents=True)
    write_cache(tmp_path / "build" / name / "CMakeCache.txt")
    (tmp_path / "install" / "selected").symlink_to(install, target_is_directory=True)

    info = selected_build_info(tmp_path)

    assert info["name"] == name
    assert info["source"] == "selected"
    assert info["status"] == "ready"
    assert "Fortran compiler: nvfortran" in info["title"]
    assert "GPU backend: openacc" in info["title"]
    assert "GPTL: ON" in info["title"]
    assert "Explicit -exe or -install_dir" in info["title"]

    (tmp_path / "install" / "selected").unlink()
    (tmp_path / "install" / "latest").symlink_to(install, target_is_directory=True)
    fallback = selected_build_info(tmp_path)
    assert fallback["source"] == "latest"
    assert fallback["status"] == "fallback"
    assert "run_scm.py will use install/latest" in fallback["title"]


def test_selected_build_reports_missing_and_broken_defaults(tmp_path):
    missing = selected_build_info(tmp_path)
    assert missing["status"] == "missing"
    assert missing["name"] == "No build selected"

    install = tmp_path / "install"
    install.mkdir()
    (install / "selected").symlink_to(install / "deleted-build", target_is_directory=True)
    broken = selected_build_info(tmp_path)
    assert broken["status"] == "broken"
    assert broken["name"] == "deleted-build"
    assert "selected install target is missing" in broken["title"]


def test_execution_controls_include_read_only_selected_build_badges():
    assert "run-selected-build-badge" in component_ids(build_run_action_section())
    assert "tune-selected-build-badge" in component_ids(
        build_top_controls(
            {
                "batch_size": 8,
                "max_workers": 4,
                "random_max_samples": 2000,
                "resolve_spacing": 0.1,
                "simann_max_iters": 200,
                "simann_initial_temp": 1.0,
                "simann_final_temp": 1.0e-12,
            }
        )
    )
    badge = selected_build_badge("profile-selected-build-badge")
    assert badge.id == "profile-selected-build-badge"
    assert badge.title
    assert not hasattr(badge, "n_clicks")


def test_shared_refresh_callback_updates_all_execution_badges():
    app = Dash(__name__, suppress_callback_exceptions=True)
    app.layout = html.Div(
        [
            dcc.Interval(id="selected-build-refresh"),
            dcc.Store(id="compile-discovery", data={}),
            *[selected_build_badge(component_id) for component_id in BADGE_IDS],
        ]
    )
    register_selected_build_callback(app)

    callback_key = next(iter(app.callback_map))
    assert all(f"{component_id}.children" in callback_key for component_id in BADGE_IDS)
