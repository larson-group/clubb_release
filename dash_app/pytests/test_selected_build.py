from pathlib import Path

from dash import Dash, dcc, html
from dash.development.base_component import Component

from dash_app.run_tab.layout import build_run_action_section
from dash_app.compile_tab.build_selector import (
    BUILD_SELECTOR_TRIGGER_IDS,
    build_implementation_capability,
    build_selector_overlay,
    build_selector_trigger,
    register_build_selector_position_callback,
    selected_build_info,
    selector_popover_style,
)
from dash_app.compile_tab import callbacks
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


def component_id_count(component, target: str) -> int:
    if not isinstance(component, Component):
        return 0
    count = int(getattr(component, "id", None) == target)
    children = getattr(component, "children", None)
    if children is None:
        return count
    return count + sum(
        component_id_count(child, target)
        for child in (children if isinstance(children, (list, tuple)) else [children])
    )


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


def test_execution_controls_include_shared_build_selector_triggers():
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
    badge = build_selector_trigger("profile-selected-build-badge")
    assert badge.id == "profile-selected-build-badge"
    assert badge.title
    assert badge.n_clicks == 0
    assert badge.type == "button"
    assert badge.children[0].children == "FORTRAN"


def test_shared_overlay_is_single_and_position_callback_owns_its_anchor():
    app = Dash(__name__, suppress_callback_exceptions=True)
    app.layout = html.Div(
        [
            build_selector_overlay(),
            *[build_selector_trigger(component_id) for component_id in BUILD_SELECTOR_TRIGGER_IDS],
        ]
    )
    register_build_selector_position_callback(app)

    callback_key = next(iter(app.callback_map))
    assert "compile-build-selector-anchor.data" in callback_key
    assert component_id_count(app.layout, "compile-build-selector-root") == 1


def test_selector_position_stays_in_the_viewport_and_can_open_above():
    below = selector_popover_style(
        {"left": 900, "top": 50, "bottom": 90, "viewport_width": 1000, "viewport_height": 800}
    )
    above = selector_popover_style(
        {"left": 20, "top": 700, "bottom": 740, "viewport_width": 1000, "viewport_height": 800}
    )

    assert below["left"] == "472.0px"
    assert below["top"] == "96.0px"
    assert above["bottom"] == "106.0px"


def test_compact_selector_shows_only_build_names_and_rebuild_controls(monkeypatch):
    monkeypatch.setattr(callbacks, "job_process_is_live", lambda _job: False)
    build = {
        "name": "gcc_PRECdouble_PYTHON",
        "path": "/build/gcc_PRECdouble_PYTHON",
        "install_prefix": "/install/gcc_PRECdouble_PYTHON",
        "install_exists": True,
        "install_prefix_mismatch": False,
        "is_selected": True,
    }
    menu = callbacks.render_compact_build_selector(
        {"builds": [build]},
        {"statuses": {build["path"]: {"status": "current", "label": "current"}}},
    )

    implementation_panel, row = menu
    assert [button.children for button in implementation_panel.children[1].children] == [
        "Fortran", "Python", "JAX"
    ]
    select, rebuild = row.children
    assert select.children == build["name"]
    assert select.id["type"] == "compile-build-selector-select"
    assert rebuild.id["type"] == "compile-build-selector-rebuild"
    assert "compile-build-card-current" in row.className
    assert "selected" not in str(select.children).lower()


def test_build_capabilities_follow_installed_runtime_contents(tmp_path):
    install = tmp_path / "install" / "test"
    install.mkdir(parents=True)
    assert build_implementation_capability(install, "fortran", tmp_path)[0] is False
    (install / "clubb_standalone").touch()
    assert build_implementation_capability(install, "fortran", tmp_path)[0] is True

    runtime = install / "python"
    (runtime / "clubb_python").mkdir(parents=True)
    (runtime / "clubb_f2py.test.so").touch()
    (runtime / "libclubb_f2py_backend.so").touch()
    assert build_implementation_capability(install, "python", tmp_path)[0] is True
    assert build_implementation_capability(install, "jax", tmp_path)[0] is False
    jax_driver = tmp_path / "clubb_jax" / "clubb_standalone.py"
    jax_driver.parent.mkdir()
    jax_driver.touch()
    assert build_implementation_capability(install, "jax", tmp_path)[0] is True
