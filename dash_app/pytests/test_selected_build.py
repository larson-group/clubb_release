import json
import subprocess
from pathlib import Path

import pytest

from dash import Dash, dcc, html
from dash.development.base_component import Component

from dash_app.run_tab.layout import build_run_action_section
from dash_app.compile_tab.build_selector import (
    BUILD_SELECTOR_TRIGGER_IDS,
    build_implementation_capability,
    build_selector_overlay,
    build_selector_trigger,
    inspect_jax_runtime_profiles,
    register_build_selector_position_callback,
    selected_build_info,
    selected_launch_target,
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


def component_text(component) -> str:
    if isinstance(component, str):
        return component
    if not isinstance(component, Component):
        return ""
    children = getattr(component, "children", None)
    values = children if isinstance(children, (list, tuple)) else [children]
    return " ".join(component_text(child) for child in values if child is not None)


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
    assert component_id_count(app.layout, "compile-build-selector-help") == 1


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
    jax_root = tmp_path / "clubb_jax"
    jax_driver = jax_root / "src" / "clubb_standalone.py"
    jax_driver.parent.mkdir(parents=True)
    jax_driver.touch()
    wrapper = jax_root / "run_jax.py"
    wrapper.touch(mode=0o755)
    (jax_root / "requirements.txt").touch()
    assert build_implementation_capability(install, "jax", tmp_path)[0] is True


def test_jax_launch_target_does_not_require_a_compiled_install(tmp_path):
    jax_root = tmp_path / "clubb_jax"
    (jax_root / "src").mkdir(parents=True)
    (jax_root / "src" / "clubb_standalone.py").touch()
    (jax_root / "requirements.txt").touch()
    (jax_root / "run_jax.py").touch(mode=0o755)

    target = selected_launch_target("jax", tmp_path)

    assert target == {
        "implementation": "jax",
        "jax_profile": "cpu",
        "install_dir": "",
        "build_name": "CPU environment",
    }


def test_jax_selector_replaces_compiled_build_rows_with_managed_runtime():
    menu = callbacks.render_compact_build_selector(
        {"builds": []},
        implementation="jax",
    )

    assert len(menu) == 2
    assert menu[0].children[2] is None
    toolbar = menu[0].children[0]
    assert toolbar.children[0].children == "Run with"
    assert toolbar.children[1].children == "?"
    assert toolbar.children[1].className == "plots-card-help"
    assert menu[1].children[0].children == "Compute"
    buttons = menu[1].children[1].children
    assert [button.id["index"] for button in buttons] == ["cpu", "gpu"]
    assert component_text(buttons[0]).startswith("CPU ")
    assert component_text(buttons[1]).startswith("GPU ")
    assert buttons[0].to_plotly_json()["props"]["aria-pressed"] == "true"
    assert buttons[1].to_plotly_json()["props"]["aria-pressed"] == "false"
    assert "compile-run-implementation-choice-selected" in buttons[0].className


def test_jax_selector_displays_wrapper_metadata_and_disables_unavailable_gpu():
    runtime_info = {
        "cpu": {
            "status": "ready",
            "selectable": True,
            "hardware": {"cpu": {"model": "Test CPU", "logical_cpus": 16}},
            "runtime": {
                "python": {"version": "3.12.4", "installed": True},
                "jax": {"required": "0.11.0", "installed": "0.11.0"},
            },
        },
        "gpu": {
            "status": "unavailable",
            "selectable": False,
            "reason": "CUDA 13 requires NVIDIA driver 580 or newer; detected 470.239",
            "hardware": {
                "cpu": {},
                "gpus": [
                    {
                        "name": "Test GPU",
                        "memory_mib": 24576,
                        "driver_version": "470.239",
                        "compute_capability": "8.0",
                    }
                ],
            },
            "runtime": {
                "cuda_major": 13,
                "python": {"version": "3.12.4", "installed": True},
                "jax": {"required": "0.11.0", "installed": "0.11.0"},
            },
        },
    }

    menu = callbacks.render_compact_build_selector(
        {"builds": []}, implementation="jax", jax_runtime_info=runtime_info
    )

    buttons = menu[1].children[1].children
    assert buttons[0].disabled is False
    assert buttons[1].disabled is True
    text = component_text(menu[1])
    assert "Test CPU | 16 logical CPUs" in text
    assert "Test GPU · 24 GiB" in text
    assert "Unavailable" in text
    assert "Python 3.12.4" not in text
    assert "requires NVIDIA driver 580" in buttons[1].title
    help_text = component_text(callbacks.render_build_selector_help(runtime_info))
    assert "Test GPU | 24.0 GiB | driver 470.239 | compute 8.0" in help_text
    assert "Python 3.12.4 | JAX 0.11.0 | CUDA 13 | Unavailable" in help_text
    assert "requires NVIDIA driver 580" in help_text
    assert "CUDA_VISIBLE_DEVICES" in help_text


def test_profile_card_describes_selected_gpu_not_first_inventory_entry():
    info = {
        "status": "ready",
        "hardware": {
            "gpus": [{"name": "Test GPU A"}, {"name": "Test GPU B"}],
            "selected_gpu": {"name": "Test GPU B", "memory_mib": 8192},
        },
    }
    card = callbacks.render_jax_profile_info("gpu", info, compact=True)
    assert component_text(card) == "GPU Ready Test GPU B · 8 GiB"
    del info["hardware"]["selected_gpu"]
    assert "2 NVIDIA GPUs available" in component_text(
        callbacks.render_jax_profile_info("gpu", info, compact=True)
    )


def test_runtime_help_ignores_mount_events_and_opens_and_closes(monkeypatch):
    app = Dash(__name__, suppress_callback_exceptions=True)
    callbacks.register_compile_callbacks(app)
    toggle = app.callback_map["compile-build-selector-help.children"]["callback"].__wrapped__
    monkeypatch.setattr(callbacks, "clicked_trigger_id", lambda: None)
    assert toggle([0], [], {}) is callbacks.no_update
    monkeypatch.setattr(callbacks, "clicked_trigger_id", lambda: {
        "type": "compile-selector-help-open", "index": "runtime",
    })
    help_card = toggle([1], [], {})
    assert help_card.className == "shared-notecard-overlay"
    assert "Choosing a runtime" in component_text(help_card)
    assert "Close" in component_text(help_card)
    monkeypatch.setattr(callbacks, "clicked_trigger_id", lambda: {
        "type": "compile-selector-help-close", "index": "runtime",
    })
    assert toggle([1], [1], {}) == ""


def test_jax_selector_checks_each_profile_source_independently(monkeypatch):
    monkeypatch.setattr(
        callbacks,
        "build_implementation_capability",
        lambda _install, _implementation, *, jax_profile: (
            (True, "")
            if jax_profile == "cpu"
            else (False, "GPU requirements are missing")
        ),
    )

    menu = callbacks.render_compact_build_selector(
        {"builds": []}, implementation="jax", jax_profile="cpu"
    )

    buttons = menu[1].children[1].children
    assert buttons[0].disabled is False
    assert buttons[1].disabled is True
    assert buttons[1].title == "GPU requirements are missing"


def test_dash_runtime_inspection_accepts_only_wrapper_schema(tmp_path, monkeypatch):
    wrapper = tmp_path / "clubb_jax" / "run_jax.py"
    wrapper.parent.mkdir(parents=True)
    wrapper.touch(mode=0o755)

    def fake_run(command, **_kwargs):
        profile = command[1].split("=", 1)[1]
        payload = {
            "schema_version": 1,
            "profile": profile,
            "selectable": profile == "cpu",
            "status": "ready" if profile == "cpu" else "unavailable",
        }
        return subprocess.CompletedProcess(command, 0, json.dumps(payload), "")

    monkeypatch.setattr("dash_app.compile_tab.build_selector.subprocess.run", fake_run)

    info = inspect_jax_runtime_profiles(tmp_path)

    assert info["cpu"]["selectable"] is True
    assert info["gpu"]["selectable"] is False


def test_jax_gpu_capability_and_launch_target_use_native_requirements(tmp_path):
    import platform

    jax_root = tmp_path / "clubb_jax"
    (jax_root / "src").mkdir(parents=True)
    (jax_root / "src" / "clubb_standalone.py").touch()
    (jax_root / "requirements.txt").touch()
    (jax_root / "run_jax.py").touch(mode=0o755)

    available, reason = build_implementation_capability(
        "", "jax", tmp_path, jax_profile="gpu"
    )
    assert available is False
    assert "requirements" in reason.lower()

    requirements_name = (
        "requirements-metal.txt"
        if platform.system() == "Darwin"
        else "requirements-cuda13.txt"
    )
    (jax_root / requirements_name).touch()
    target = selected_launch_target("jax", tmp_path, jax_profile="GPU")
    assert target["jax_profile"] == "gpu"
    assert target["build_name"] == "GPU environment"


@pytest.mark.parametrize("failure", ["exit", "timeout", "json", "schema"])
def test_dash_runtime_probe_failure_is_visible(tmp_path, monkeypatch, failure):
    wrapper = tmp_path / "clubb_jax" / "run_jax.py"
    wrapper.parent.mkdir()
    wrapper.touch(mode=0o755)

    def fake_run(command, **kwargs):
        if failure == "timeout":
            raise subprocess.TimeoutExpired(command, kwargs["timeout"])
        if failure == "exit":
            return subprocess.CompletedProcess(command, 1, "", "probe failed")
        return subprocess.CompletedProcess(command, 0, "not json" if failure == "json" else "{}", "")

    monkeypatch.setattr("dash_app.compile_tab.build_selector.subprocess.run", fake_run)
    for profile, info in inspect_jax_runtime_profiles(tmp_path).items():
        assert info["status"] == "unknown"
        assert info["selectable"] is False
        assert info["reason"]
        card = component_text(callbacks.render_jax_profile_info(profile, info, compact=True))
        assert "Checking" not in card
        assert "Unknown" in card


def test_metal_preallocation_is_disabled():
    info = {"gpu": {"status": "setup_required", "selectable": True,
                    "runtime": {"accelerator": "metal"}}}
    assert not callbacks.jax_preallocation_available(info, "gpu", "")
