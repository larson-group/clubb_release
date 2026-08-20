"""Focused tests for the small named-config save helper."""

from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest
from dash import Dash, dcc

from dash_app.run_tab import callbacks_settings
from dash_app.run_tab.namelist import read_namelist_entries
from dash_app.run_tab.callbacks_settings import register_settings_callbacks, settings_resolution_for_save
from dash_app.run_tab.config_state import build_tunable_config_state
from dash_app.run_tab.layout import build_config_save_dialog, build_run_config_buttons
from utilities.save_tunable_config import build_tuned_config_overrides, save_tunable_config


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_CONFIG_ROOT = REPO_ROOT / "input" / "parameter_and_flag_configs"


@pytest.fixture
def config_root(tmp_path):
    root = tmp_path / "configs"
    root.mkdir()
    shutil.copytree(SOURCE_CONFIG_ROOT / "default", root / "default")
    return root


def _value(path, setting):
    values = {entry["name"]: entry["value"] for entry in read_namelist_entries(path)}
    return values[setting]


def test_save_clones_config_applies_overrides_and_writes_readme(config_root):
    result = save_tunable_config(
        "default",
        {
            "flags": {"l_use_cloud_cover": ".true."},
            "tunable": {"C4": "2.5"},
        },
        "my_config",
        "ARM experiment",
        config_root=config_root,
    )

    saved = config_root / "my_config"
    assert result["name"] == "my_config"
    assert _value(saved / "configurable_model_flags.in", "l_use_cloud_cover") == ".true."
    assert _value(saved / "tunable_parameters.in", "C4") == "2.5"
    assert _value(config_root / "default" / "tunable_parameters.in", "C4") == "2.000000"
    readme = (saved / "README.md").read_text(encoding="utf-8")
    assert "Cloned from: `default`" in readme
    assert "ARM experiment" in readme
    assert "CLUBB commit:" in readme
    assert "`tunable_parameters.in:C4`" in readme


def test_tuned_config_overrides_apply_candidate_after_scm_override(config_root):
    overrides = build_tuned_config_overrides(
        "default",
        "l_use_cloud_cover=.true.,C4=2.5",
        {"C4": 3.25},
        config_root=config_root,
    )
    save_tunable_config(
        "default",
        overrides,
        "tuned",
        config_root=config_root,
    )

    saved = config_root / "tuned"
    assert _value(saved / "configurable_model_flags.in", "l_use_cloud_cover") == ".true."
    assert _value(saved / "tunable_parameters.in", "C4") == "3.25"


def test_tuned_config_rejects_overrides_not_stored_in_named_configs(config_root):
    with pytest.raises(ValueError, match="not stored in a tunable config"):
        build_tuned_config_overrides(
            "default",
            "dt_main=10",
            {"C4": 3.25},
            config_root=config_root,
        )


def test_overwrite_requires_force(config_root):
    save_tunable_config("default", {}, "saved", config_root=config_root)

    with pytest.raises(FileExistsError, match="already exists"):
        save_tunable_config("default", {}, "saved", config_root=config_root)

    result = save_tunable_config(
        "default",
        {"tunable": {"C4": "3.0"}},
        "saved",
        force=True,
        config_root=config_root,
    )
    assert result["overwritten"] is True
    assert _value(config_root / "saved" / "tunable_parameters.in", "C4") == "3.0"

    save_tunable_config(
        "saved",
        {"tunable": {"C4": "4.0"}},
        "saved",
        force=True,
        config_root=config_root,
    )
    assert _value(config_root / "saved" / "tunable_parameters.in", "C4") == "4.0"


@pytest.mark.parametrize("name", ["", "../escape", "two words", ".hidden", "name/child"])
def test_save_rejects_unsafe_names(config_root, name):
    with pytest.raises(ValueError, match="Config name"):
        save_tunable_config("default", {}, name, config_root=config_root)


def test_save_button_is_last_and_visually_distinct():
    choices = build_run_config_buttons(
        [{"label": "default", "value": "default"}],
        "default",
    )
    save_button = choices.children[-1]
    assert save_button.id == "run-config-save"
    assert save_button.children == "Save config"
    assert "dashed" in save_button.style["border"]


def _component_by_id(root, component_id):
    pending = [root]
    while pending:
        component = pending.pop()
        if getattr(component, "id", None) == component_id:
            return component
        children = getattr(component, "children", None)
        if isinstance(children, (list, tuple)):
            pending.extend(children)
        elif children is not None:
            pending.append(children)
    raise AssertionError(f"Missing component: {component_id}")


def test_save_dialog_has_name_note_and_footer_actions():
    dialog = build_config_save_dialog()

    assert "run-config-save-modal--hidden" in dialog.className
    assert isinstance(_component_by_id(dialog, "run-config-save-name"), dcc.Input)
    assert isinstance(_component_by_id(dialog, "run-config-save-note"), dcc.Textarea)
    assert _component_by_id(dialog, "run-config-save-submit").children == "Save config"
    assert _component_by_id(dialog, "run-config-save-cancel").children == "Cancel"


def test_existing_name_requires_second_overwrite_click(monkeypatch):
    app = Dash(__name__, suppress_callback_exceptions=True)
    register_settings_callbacks(app)
    callback = next(
        entry["callback"].__wrapped__
        for entry in app.callback_map.values()
        if entry["callback"].__name__ == "manage_config_save_dialog"
    )
    monkeypatch.setattr(
        callbacks_settings,
        "callback_context",
        SimpleNamespace(triggered_id="run-config-save-submit"),
    )
    configs = [{"label": "default", "value": "default"}]

    warning = callback(1, 0, 1, "default", "my note", configs, None)
    assert warning[4] is callbacks_settings.no_update
    assert warning[5] == "default"
    assert warning[6] == "Overwrite config"
    assert "danger" in warning[7]

    confirmed = callback(1, 0, 2, "default", "my note", configs, "default")
    assert "run-config-save-modal--hidden" in confirmed[0]
    assert confirmed[4]["name"] == "default"
    assert confirmed[4]["note"] == "my note"
    assert confirmed[4]["force"] is True


def test_save_resolution_uses_live_control_values():
    state = build_tunable_config_state("default")
    schema = state["settings_schema"]
    initial = state["settings_resolution"]
    flag_ids = [{"name": name} for name in schema["flag_defaults"]]
    flag_values = [["on"] if initial["normalized_flags"][item["name"]] else [] for item in flag_ids]
    param_ids = [
        {"file": item["file"], "name": item["name"]}
        for item in state["param_meta"]
    ]
    param_values = [
        initial["normalized_parameters"][item["file"]][item["name"]]
        for item in param_ids
    ]
    c4_index = param_ids.index({"file": "tunable", "name": "C4"})
    param_values[c4_index] = "2.5"

    resolution = settings_resolution_for_save(
        schema,
        flag_ids,
        flag_values,
        param_ids,
        param_values,
    )

    assert resolution["overrides"]["tunable"]["C4"] == "2.5"


def test_save_resolution_uses_visible_linked_control_values():
    state = build_tunable_config_state("default")
    schema = state["settings_schema"]
    initial = state["settings_resolution"]
    flag_ids = [{"name": name} for name in schema["flag_defaults"]]
    flag_values = [["on"] if initial["normalized_flags"][item["name"]] else [] for item in flag_ids]
    param_ids = [
        {"file": item["file"], "name": item["name"]}
        for item in state["param_meta"]
    ]
    param_values = [
        initial["normalized_parameters"][item["file"]][item["name"]]
        for item in param_ids
    ]

    resolution = settings_resolution_for_save(
        schema,
        flag_ids,
        flag_values,
        param_ids,
        param_values,
        [{"group": "C6rt=C6thl"}],
        ["3.25"],
    )

    assert resolution["overrides"]["tunable"]["C6rt"] == "3.25"
    assert resolution["overrides"]["tunable"]["C6thl"] == "3.25"


def test_save_completion_refreshes_choices_without_reselecting_config(monkeypatch):
    app = Dash(__name__, suppress_callback_exceptions=True)
    register_settings_callbacks(app)
    callback_entry = next(
        entry
        for entry in app.callback_map.values()
        if entry["callback"].__name__ == "save_named_config"
    )
    callback = callback_entry["callback"].__wrapped__
    configs = [
        {"label": "default", "value": "default"},
        {"label": "saved", "value": "saved"},
    ]
    monkeypatch.setattr(
        callbacks_settings,
        "settings_resolution_for_save",
        lambda *_args: {"issues": [], "overrides": {"flags": {"l_predict_upwp_vpwp": ".true."}}},
    )
    monkeypatch.setattr(
        callbacks_settings,
        "save_tunable_config",
        lambda *_args, **_kwargs: {"name": "saved", "overwritten": False, "edits": ["flag"]},
    )
    monkeypatch.setattr(callbacks_settings, "available_tunable_configs", lambda: configs)

    result = callback(
        {"name": "saved", "nonce": 1},
        "default",
        {},
        [],
        [],
        [],
        [],
        [],
        [],
    )

    output_ids = [output.component_id for output in callback_entry["output"]]
    assert "run-config-reset-signal" not in output_ids
    assert "run-config-save-selection" not in output_ids
    assert result[0] == result[1] == configs
    assert [button.id for button in result[2]][-2:] == [
        {"type": "run-config-button", "name": "saved"},
        "run-config-save",
    ]
