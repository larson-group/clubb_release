"""Clone a named CLUBB config and apply Run-tab overrides."""

from datetime import datetime, timezone
from pathlib import Path
import re
import shutil
import subprocess

from dash_app.run_tab.namelist import apply_updates_to_lines, read_namelist_entries
from utilities.create_case_namelist import normalize_override_string, parse_override_pairs


REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_ROOT = REPO_ROOT / "input" / "parameter_and_flag_configs"
CONFIG_FILES = {
    "flags": "configurable_model_flags.in",
    "tunable": "tunable_parameters.in",
    "silhs": "silhs_parameters.in",
}
CONFIG_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


def normalize_config_name(value):
    """Return a safe config directory name or raise a user-facing error."""
    name = str(value or "").strip()
    if not CONFIG_NAME_RE.fullmatch(name):
        raise ValueError(
            "Config name must start with a letter or digit and contain only "
            "letters, digits, underscores, or hyphens."
        )
    return name


def _git_hash():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=5,
        ).strip()
    except (OSError, subprocess.SubprocessError):
        return "unavailable"


def build_tuned_config_overrides(
    source_config,
    scm_override,
    tuned_params,
    *,
    config_root=None,
):
    """Group a Tune run's SCM override and selected parameters by config file."""
    root = Path(config_root or CONFIG_ROOT)
    source_name = normalize_config_name(source_config)
    source = root / source_name
    locations = {}
    for group, filename in CONFIG_FILES.items():
        for entry in read_namelist_entries(source / filename):
            locations.setdefault(entry["name"].casefold(), []).append((group, entry["name"]))

    assignments = dict(parse_override_pairs(normalize_override_string(scm_override)))
    assignments.update(dict(tuned_params or {}))
    grouped = {group: {} for group in CONFIG_FILES}
    unknown = []
    ambiguous = []
    for raw_name, value in assignments.items():
        name = str(raw_name).strip()
        matches = locations.get(name.casefold(), [])
        if not matches:
            unknown.append(name)
            continue
        if len(matches) != 1:
            ambiguous.append(name)
            continue
        group, canonical_name = matches[0]
        grouped[group][canonical_name] = str(value).strip()
    if unknown:
        raise ValueError(
            "SCM override(s) are not stored in a tunable config: "
            + ", ".join(sorted(unknown))
        )
    if ambiguous:
        raise ValueError(
            "Config setting(s) occur in more than one namelist: "
            + ", ".join(sorted(ambiguous))
        )
    return grouped


def save_tunable_config(
    source_config,
    overrides,
    name,
    note="",
    force=False,
    *,
    config_root=None,
):
    """Copy ``source_config`` to ``name``, apply overrides, and add metadata."""
    root = Path(config_root or CONFIG_ROOT)
    source_name = normalize_config_name(source_config)
    destination_name = normalize_config_name(name)
    source = root / source_name
    destination = root / destination_name

    missing = [filename for filename in CONFIG_FILES.values() if not (source / filename).is_file()]
    if missing:
        raise ValueError(f"Source config {source_name} is incomplete: {', '.join(missing)}")
    if destination.is_symlink() or (destination.exists() and not destination.is_dir()):
        raise ValueError(f"Config destination is not a normal directory: {destination}")
    overwritten = destination.exists()
    if overwritten and not force:
        raise FileExistsError(f"Config {destination_name} already exists")

    requested = dict(overrides or {})
    unknown_groups = sorted(set(requested) - set(CONFIG_FILES))
    if unknown_groups:
        raise ValueError(f"Unknown config override group: {', '.join(unknown_groups)}")
    if source.resolve() != destination.resolve():
        shutil.copytree(source, destination, dirs_exist_ok=bool(force))

    edits = []
    for group, filename in CONFIG_FILES.items():
        updates = dict(requested.get(group) or {})
        if not updates:
            continue
        path = destination / filename
        old_values = {entry["name"]: entry["value"] for entry in read_namelist_entries(path)}
        with path.open("r", encoding="utf-8") as handle:
            lines = handle.readlines()
        with path.open("w", encoding="utf-8") as handle:
            handle.writelines(apply_updates_to_lines(lines, updates))
        edits.extend(
            (filename, setting, old_values.get(setting, "(missing)"), str(value))
            for setting, value in updates.items()
        )

    readme = [
        f"# CLUBB config: {destination_name}",
        "",
        f"- Cloned from: `{source_name}`",
        f"- Saved at: `{datetime.now(timezone.utc).replace(microsecond=0).isoformat()}`",
        f"- CLUBB commit: `{_git_hash()}`",
        f"- Overwrite: `{'yes' if overwritten else 'no'}`",
    ]
    if note:
        readme.extend(["", "## Note", "", str(note).strip()])
    readme.extend(["", "## Edits", ""])
    readme.extend(
        [f"- `{filename}:{setting}`: `{old}` → `{new}`" for filename, setting, old, new in edits]
        or ["No values differed from the cloned config."]
    )
    (destination / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    return {
        "name": destination_name,
        "path": str(destination),
        "overwritten": overwritten,
        "edits": edits,
    }
