"""Save parameter and flag edits as named JSON preset overrides.

normalize_config_name() checks UI-provided names. build_tuned_config_overrides()
assigns tuner edits to model parameters, SILHS parameters or flags.
save_tunable_config() resolves the source preset, applies typed changes and
stores sparse overrides relative to each catalog's defaults. Both catalogs
record source, note, revision, time and edits under preset_metadata. Replacing
an existing name requires force=True. No namelist fragments are created.
"""
from datetime import datetime, timezone
from pathlib import Path
import re
import subprocess

from input.scripts.namelist import updated_settings, atomic_text, settings_patch, updated_json_text
from input.scripts.parameters_json_to_namelist import parameter_settings
from input.scripts.flags_json_to_namelist import flag_settings
from input.scripts.create_case_namelist import CATALOG_FILES, DEFAULT_CONFIG, normalize_override_string, parse_override_pairs

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_ROOT = REPO_ROOT / "input"
CONFIG_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


def normalize_config_name(value):
    """Return a safe config name or raise a user-facing validation error."""
    name = str(value or "").strip()
    if not CONFIG_NAME_RE.fullmatch(name):
        raise ValueError("Config name must start with a letter or digit and contain only letters, digits, underscores, or hyphens.")
    return name


def _git_hash():
    """Return the current source revision for snapshot provenance when available."""
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT,
            text=True, stderr=subprocess.DEVNULL, timeout=5).strip()
    except (OSError, subprocess.SubprocessError):
        return "unavailable"


def _source_settings(source_config, root):
    """Resolve the baseline or matching presets, inheriting missing sections."""
    from input.scripts.create_case_namelist import load_catalog
    name = normalize_config_name(source_config)
    if name != DEFAULT_CONFIG and not any(name in load_catalog(kind, root / CATALOG_FILES[kind])["presets"]
               for kind in ("parameters", "flags")):
        raise ValueError(f"Unknown tunable config {name!r}")
    return (parameter_settings(config=name, json_file=root / "clubb_params.jsonc"),
            flag_settings(config=name, json_file=root / "clubb_config_flags.jsonc"))


def build_tuned_config_overrides(source_config, scm_override, tuned_params, *, config_root=None):
    """Group SCM edits followed by candidate values; reject unknown/ambiguous names."""
    params, flags = _source_settings(source_config, Path(config_root or CONFIG_ROOT))
    locations = {}
    for group, values in (("flags", {**flags["model"], **flags["silhs"]}),
                           ("tunable", params["model"]), ("silhs", params["silhs"])):
        for name in values:
            locations.setdefault(name.casefold(), []).append((group, name))
    assignments = dict(parse_override_pairs(normalize_override_string(scm_override)))
    assignments.update(dict(tuned_params or {}))
    grouped = {group: {} for group in ("flags", "tunable", "silhs")}
    for raw_name, value in assignments.items():
        matches = locations.get(str(raw_name).strip().casefold(), [])
        if len(matches) != 1:
            raise ValueError(f"SCM override(s) are not stored in a tunable config or are ambiguous: {raw_name}")
        group, name = matches[0]
        grouped[group][name] = str(value).strip()
    return grouped


def save_tunable_config(source_config, overrides, name, note="", force=False, *, config_root=None):
    """Save typed JSON snapshots and provenance after validating every requested edit."""
    root = Path(config_root or CONFIG_ROOT)
    source_name = normalize_config_name(source_config)
    destination_name = normalize_config_name(name)
    if destination_name == DEFAULT_CONFIG:
        raise ValueError(f"Config name {DEFAULT_CONFIG!r} is reserved for catalog defaults; choose another name.")
    from input.scripts.create_case_namelist import load_catalog
    catalogs = {kind: load_catalog(kind, root / CATALOG_FILES[kind]) for kind in ("parameters", "flags")}
    overwritten = any(destination_name in catalog["presets"] for catalog in catalogs.values())
    if overwritten and not force:
        raise FileExistsError(f"Config {destination_name} already exists")
    params, flags = _source_settings(source_name, root)
    requested = dict(overrides or {})
    if set(requested) - {"flags", "tunable", "silhs"}:
        raise ValueError("Unknown config override group")
    edits = []
    for group in ("flags", "tunable", "silhs"):
        updates = dict(requested.get(group) or {})
        settings = flags if group == "flags" else params
        categories = ("model", "silhs") if group == "flags" else ("model",) if group == "tunable" else ("silhs",)
        changed = updated_settings(settings, updates, categories)
        filename = "clubb_config_flags.jsonc" if group == "flags" else "clubb_params.jsonc"
        for setting, value in updates.items():
            category, key = next((c, k) for c in categories for k in settings[c] if k.casefold() == setting.casefold())
            edits.append((filename, key, settings[category][key], changed[category][key]))
        if group == "flags":
            flags = changed
        else:
            params = changed
    metadata = {"source": source_name, "note": str(note).strip(), "commit": _git_hash(),
                "saved_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
                "edits": edits}
    for kind, settings in (("parameters", params), ("flags", flags)):
        catalog = catalogs[kind]
        catalog["presets"][destination_name] = settings_patch(catalog["defaults"], settings)
        catalog.setdefault("preset_metadata", {})[destination_name] = metadata
        path = root / CATALOG_FILES[kind]
        atomic_text(path, updated_json_text(path, catalog))
    return {"name": destination_name, "path": str(root), "overwritten": overwritten, "edits": edits}
