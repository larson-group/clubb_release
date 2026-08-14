"""Static filesystem discovery for the Profile tab."""

from __future__ import annotations

from pathlib import Path

from dash_app.run_tab.discovery import load_available_cases
from dash_app.shared.tunable_configs import (
    available_tunable_configs,
    default_tunable_config_name,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
INSTALL_ROOT = REPO_ROOT / "install"


def available_executables() -> list[dict[str, str]]:
    """Return installed standalone executables without duplicating symlink targets."""
    options = [{"label": "Auto (selected/latest install)", "value": ""}]
    seen: set[Path] = set()
    if not INSTALL_ROOT.is_dir():
        return options
    for executable in sorted(INSTALL_ROOT.glob("*/clubb_standalone")):
        if not executable.is_file():
            continue
        resolved = executable.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        options.append(
            {
                "label": executable.parent.name,
                "value": str(resolved),
            }
        )
    return options


def discover_profile_state() -> dict[str, object]:
    cases = load_available_cases()
    configs = available_tunable_configs()
    return {
        "cases": cases,
        "default_case": "arm" if "arm" in cases else (cases[0] if cases else None),
        "configs": configs,
        "default_config": default_tunable_config_name(configs),
        "executables": available_executables(),
    }
