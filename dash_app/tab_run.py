import ast
import json
import os
import re
import resource
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading
import time

from dash import dcc, html, Input, Output, State, ALL, callback_context, no_update

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CASE_DIR = os.path.join(REPO_ROOT, "input", "case_setups")
STATS_DIR = os.path.join(REPO_ROOT, "input", "stats")
DEFAULT_STATS_NAME = "standard_stats.in"
NO_STATS_NAME = "none"
FLAGS_FILE = os.path.join(REPO_ROOT, "input", "tunable_parameters", "configurable_model_flags.in")
TUNABLE_FILE = os.path.join(REPO_ROOT, "input", "tunable_parameters", "tunable_parameters.in")
SILHS_FILE = os.path.join(REPO_ROOT, "input", "tunable_parameters", "silhs_parameters.in")
RUN_PROCS = {}
RUN_ACTIVE_CASES = {}
RUN_STATUS = {}
RUN_FINALIZED = set()
RUN_LOCK = threading.Lock()
RUN_STREAM_LOCK = threading.Lock()
RUN_SCM_ALL = os.path.join(REPO_ROOT, "run_scripts", "run_scm_all.py")
ITERATION_RE = re.compile(
    r"iteration:\s*(\d+)\s*/\s*(\d+)\s*--\s*time\s*=\s*([-+0-9.eE]+)\s*/\s*([-+0-9.eE]+)",
    re.IGNORECASE,
)


def _physical_core_count():
    # Prefer physical cores where possible so the dashboard doesn't oversubscribe.
    if sys.platform.startswith("linux"):
        try:
            pairs = set()
            phys_id = None
            core_id = None
            with open("/proc/cpuinfo", "r", encoding="utf-8", errors="replace") as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if not line:
                        if phys_id is not None and core_id is not None:
                            pairs.add((phys_id, core_id))
                        phys_id = None
                        core_id = None
                        continue
                    if line.startswith("physical id"):
                        phys_id = line.split(":", 1)[1].strip()
                    elif line.startswith("core id"):
                        core_id = line.split(":", 1)[1].strip()
            if phys_id is not None and core_id is not None:
                pairs.add((phys_id, core_id))
            if pairs:
                return len(pairs)
        except Exception:
            pass

    if sys.platform == "darwin":
        try:
            out = subprocess.check_output(
                ["sysctl", "-n", "hw.physicalcpu"], text=True
            ).strip()
            cores = int(out)
            if cores > 0:
                return cores
        except Exception:
            pass

    return max(1, os.cpu_count() or 1)


MAX_RUN_PROCS = max(1, _physical_core_count() - 2)
MAX_UI_LOG_LINES = 2000
RUN_STACK_KB = int(os.environ.get("CLUBB_RUN_STACK_KB", "8192000"))
ENABLE_CUDA_MPS = os.environ.get("CLUBB_ENABLE_CUDA_MPS", "1").strip().lower() not in {"0", "false", "no"}
_uid = os.getuid() if hasattr(os, "getuid") else "user"
CUDA_MPS_PIPE_DIR = os.path.join("/tmp", f"clubb_mps_pipe_{_uid}")
CUDA_MPS_LOG_DIR = os.path.join("/tmp", f"clubb_mps_log_{_uid}")
MPS_LOCK = threading.Lock()
MPS_STATE = {
    "attempted": False,
    "enabled": False,
    "message": "",
    "env": {},
}


def _set_child_stack_limit():
    """Apply a larger stack size in the child before exec."""
    try:
        soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
        target = max(soft, RUN_STACK_KB * 1024)
        if hard != resource.RLIM_INFINITY:
            target = min(target, hard)
        resource.setrlimit(resource.RLIMIT_STACK, (target, hard))
    except Exception:
        # Keep launch robust; run with default limits if changing stack fails.
        pass


def _ensure_cuda_mps():
    """Try to start CUDA MPS once. Returns env overrides if enabled, else {}."""
    if not ENABLE_CUDA_MPS:
        return {}

    with MPS_LOCK:
        if MPS_STATE["attempted"]:
            return dict(MPS_STATE["env"])
        MPS_STATE["attempted"] = True

        mps_ctl = shutil.which("nvidia-cuda-mps-control")
        if not mps_ctl:
            MPS_STATE["message"] = "nvidia-cuda-mps-control not found; continuing without MPS."
            print(f"[run-tab] {MPS_STATE['message']}")
            return {}

        try:
            os.makedirs(CUDA_MPS_PIPE_DIR, exist_ok=True)
            os.makedirs(CUDA_MPS_LOG_DIR, exist_ok=True)
        except OSError as exc:
            MPS_STATE["message"] = f"failed to create MPS directories: {exc}; continuing without MPS."
            print(f"[run-tab] {MPS_STATE['message']}")
            return {}

        mps_env = os.environ.copy()
        mps_env["CUDA_MPS_PIPE_DIRECTORY"] = CUDA_MPS_PIPE_DIR
        mps_env["CUDA_MPS_LOG_DIRECTORY"] = CUDA_MPS_LOG_DIR

        try:
            start = subprocess.run(
                [mps_ctl, "-d"],
                env=mps_env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
        except Exception as exc:
            MPS_STATE["message"] = f"failed to launch MPS daemon: {exc}; continuing without MPS."
            print(f"[run-tab] {MPS_STATE['message']}")
            return {}

        if start.returncode != 0:
            err = (start.stderr or "").strip()
            suffix = f" ({err})" if err else ""
            MPS_STATE["message"] = f"MPS daemon start failed with code {start.returncode}{suffix}; continuing without MPS."
            print(f"[run-tab] {MPS_STATE['message']}")
            return {}

        MPS_STATE["enabled"] = True
        MPS_STATE["env"] = {
            "CUDA_MPS_PIPE_DIRECTORY": CUDA_MPS_PIPE_DIR,
            "CUDA_MPS_LOG_DIRECTORY": CUDA_MPS_LOG_DIR,
        }
        MPS_STATE["message"] = "CUDA MPS enabled for dashboard-launched runs."
        print(f"[run-tab] {MPS_STATE['message']}")
        return dict(MPS_STATE["env"])


def _run_child_env():
    child_env = os.environ.copy()
    child_env.update(_ensure_cuda_mps())
    return child_env


def _read_log_increment(log_path, offset):
    if not log_path or not os.path.exists(log_path):
        return "", 0

    try:
        with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
            handle.seek(0, os.SEEK_END)
            end_pos = handle.tell()
            if offset < 0 or offset > end_pos:
                offset = 0
            if offset == end_pos:
                return "", end_pos
            handle.seek(offset)
            return handle.read(), handle.tell()
    except Exception:
        return "", offset


def _cleanup_log_file(log_path):
    if not log_path:
        return
    try:
        if os.path.exists(log_path):
            os.remove(log_path)
    except Exception:
        pass


def _append_log_tail(existing, chunk, max_lines=MAX_UI_LOG_LINES):
    text = f"{existing or ''}{chunk or ''}"
    if not text:
        return ""
    lines = text.splitlines(keepends=True)
    if len(lines) <= max_lines:
        return text
    return "".join(lines[-max_lines:])


def _mark_case_started(case_name, proc, proc_data):
    with RUN_LOCK:
        RUN_PROCS[proc.pid] = proc
        RUN_ACTIVE_CASES[case_name] = dict(proc_data)
        RUN_STATUS.pop(case_name, None)
        RUN_FINALIZED.discard(case_name)


def _get_proc(pid):
    if not pid:
        return None
    with RUN_LOCK:
        return RUN_PROCS.get(pid)


def _record_case_finish(case_name, pid, status):
    with RUN_LOCK:
        already_finalized = case_name in RUN_FINALIZED
        if pid:
            RUN_PROCS.pop(pid, None)
        RUN_ACTIVE_CASES.pop(case_name, None)
        RUN_STATUS[case_name] = int(status)
        RUN_FINALIZED.add(case_name)
    return already_finalized


def _get_cached_status(case_name):
    with RUN_LOCK:
        return RUN_STATUS.get(case_name)


def _clear_case_status(case_name):
    with RUN_LOCK:
        RUN_STATUS.pop(case_name, None)
        RUN_FINALIZED.discard(case_name)


def _snapshot_status_lists():
    with RUN_LOCK:
        completed = sorted([name for name, status in RUN_STATUS.items() if status == 0])
        failed = sorted([name for name, status in RUN_STATUS.items() if status != 0])
    return completed, failed


def _is_case_active(case_name):
    with RUN_LOCK:
        return case_name in RUN_ACTIVE_CASES


def _snapshot_active_cases():
    with RUN_LOCK:
        return {name: dict(data) for name, data in RUN_ACTIVE_CASES.items()}


def _list_cases():
    cases = []
    if not os.path.isdir(CASE_DIR):
        return cases
    for entry in os.listdir(CASE_DIR):
        if entry.endswith("_model.in"):
            cases.append(entry[: -len("_model.in")])
    return sorted(cases)


def _list_stats_files():
    files = []
    if not os.path.isdir(STATS_DIR):
        return files
    for entry in os.listdir(STATS_DIR):
        if entry.endswith(".in"):
            files.append(entry)
    return sorted(files)


def _parse_line(line):
    if "=" not in line:
        return None
    if line.lstrip().startswith(("!", "&")):
        return None
    if line.strip() == "/":
        return None
    comment = ""
    idx = line.find("!")
    if idx != -1:
        comment = line[idx:].rstrip("\n")
        body = line[:idx]
    else:
        body = line
    match = re.match(r"^(\s*)([A-Za-z_]\w*)\s*=\s*(.*)$", body)
    if not match:
        return None
    indent, name, rest = match.groups()
    rest = rest.rstrip()
    comma = rest.endswith(",")
    if comma:
        rest = rest[:-1].rstrip()
    return {
        "indent": indent,
        "name": name,
        "value": rest.strip(),
        "comment": comment,
        "comma": comma,
    }


def _read_namelist_entries(path):
    entries = []
    if not os.path.exists(path):
        return entries
    with open(path, "r", encoding="utf-8") as handle:
        for idx, line in enumerate(handle.readlines()):
            parsed = _parse_line(line)
            if not parsed:
                continue
            parsed["line_index"] = idx
            entries.append(parsed)
    return entries


BOOL_TRUE = {".true.", "true", "t"}
BOOL_FALSE = {".false.", "false", "f"}


def _is_bool_value(value):
    return value.strip().lower() in BOOL_TRUE | BOOL_FALSE


def _is_true(value):
    return value.strip().lower() in BOOL_TRUE


def _normalize_numeric_display(value):
    text = str(value).strip()
    if not text:
        return text
    lower = text.lower()
    if lower in BOOL_TRUE or lower in BOOL_FALSE:
        return text

    # Keep exponent formatting untouched; trim trailing zeros only for plain decimals.
    if "e" in lower:
        return text
    if "." not in text:
        return text

    # Strip trailing zeros but keep one decimal place for integer-valued decimals (e.g., 1.00 -> 1.0).
    trimmed = text.rstrip("0")
    if trimmed.endswith("."):
        trimmed += "0"
    return trimmed


def _update_namelist_value(path, name, new_value):
    if not os.path.exists(path):
        return
    with open(path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    lines = _apply_updates_to_lines(lines, {name: new_value})
    with open(path, "w", encoding="utf-8") as handle:
        handle.writelines(lines)


def _apply_updates_to_lines(lines, updates):
    remaining = dict(updates or {})
    for idx, line in enumerate(lines):
        parsed = _parse_line(line)
        if not parsed:
            continue
        name = parsed["name"]
        if name not in remaining:
            continue
        new_value = remaining.pop(name)
        comma = "," if parsed["comma"] else ""
        comment = parsed["comment"]
        sep = " " if comment else ""
        updated = f"{parsed['indent']}{name} = {new_value}{comma}{sep}{comment}".rstrip()
        lines[idx] = updated + "\n"
        if not remaining:
            break
    for name, new_value in remaining.items():
        lines.append(f"{name} = {new_value}\n")
    return lines


def _write_temp_namelist(base_path, updates, prefix):
    if not updates:
        return None
    if not os.path.exists(base_path):
        return None
    with open(base_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    lines = _apply_updates_to_lines(lines, updates)
    temp_file = tempfile.NamedTemporaryFile(
        delete=False, prefix=prefix, suffix=".in", dir="/tmp", mode="w", encoding="utf-8"
    )
    temp_file.writelines(lines)
    temp_file.close()
    return temp_file.name


def _cleanup_temp_files(paths):
    if not paths:
        return
    for path in paths:
        if not path:
            continue
        try:
            if os.path.exists(path):
                os.remove(path)
        except Exception:
            pass


def _build_override_updates(flag_values, param_values, defaults_data, flag_names_data, param_meta):
    updates = {"flags": {}, "tunable": {}, "silhs": {}}
    if not defaults_data:
        return updates
    for name, values in zip(flag_names_data or [], flag_values or []):
        current = bool(values)
        default = bool(defaults_data["flags"].get(name))
        if current != default:
            updates["flags"][name] = ".true." if current else ".false."
    for meta, value in zip(param_meta or [], param_values or []):
        file_key = meta.get("file")
        name = meta.get("name")
        if not file_key or not name:
            continue
        default_value = defaults_data["params"][file_key].get(name)
        current_value = _normalize_numeric_display(value)
        if current_value != _normalize_numeric_display(default_value):
            updates[file_key][name] = current_value
    return updates


def _field_style(changed):
    base = {"padding": "4px", "display": "flex", "alignItems": "center"}
    if changed:
        base.update({"outline": "2px solid #f59e0b", "borderRadius": "4px"})
    return base


def _case_button_style(color, selected=False):
    border_color = "#f59e0b" if selected else "transparent"
    border_width = "3px" if selected else "1px"
    style = {
        "color": "#ffffff",
        "border": f"{border_width} solid {border_color}",
        "padding": "6px 10px",
        "margin": "4px",
        "borderRadius": "4px",
        "cursor": "pointer",
        "boxSizing": "border-box",
        "opacity": "1",
    }
    # Always set both keys so a previous gradient/plain value cannot linger.
    style["background"] = color
    if "gradient" in color:
        style["backgroundColor"] = "#dc2626"
    else:
        style["backgroundColor"] = color
    return style


def _run_action_button_style(color, disabled=False):
    style = {
        "backgroundColor": color,
        "color": "#ffffff",
        "border": "none",
        "padding": "10px 16px",
        "margin": "4px",
        "borderRadius": "6px",
        "cursor": "pointer",
        "fontSize": "14px",
        "fontWeight": "600",
    }
    if disabled:
        style.update(
            {
                "backgroundColor": "#9ca3af",
                "color": "#f3f4f6",
                "cursor": "not-allowed",
            }
        )
    return style


def _stats_button_style(selected=False):
    return {
        "backgroundColor": "#0ea5e9" if selected else "#e2e8f0",
        "color": "#0f172a" if selected else "#334155",
        "border": "1px solid #94a3b8",
        "padding": "6px 10px",
        "margin": "4px",
        "borderRadius": "4px",
        "cursor": "pointer",
        "fontWeight": "600" if selected else "500",
    }


def _case_dom_id(case_name):
    return re.sub(r"[^A-Za-z0-9_-]", "_", case_name)


def _format_runtime(seconds):
    if seconds is None:
        return ""
    total = max(0, int(round(float(seconds))))
    hours, rem = divmod(total, 3600)
    mins, secs = divmod(rem, 60)
    if hours > 0:
        return f"{hours}h {mins:02d}m {secs:02d}s"
    if mins > 0:
        return f"{mins}m {secs:02d}s"
    return f"{secs}s"


def _format_eta(seconds):
    if seconds is None:
        return ""
    total = max(0, int(round(float(seconds))))
    mins, secs = divmod(total, 60)
    return f"{mins}m {secs:02d}s"


def _extract_progress(log_text):
    if not log_text:
        return None
    last_match = None
    for match in ITERATION_RE.finditer(log_text):
        last_match = match
    if not last_match:
        return None
    try:
        current = int(last_match.group(1))
        total = int(last_match.group(2))
    except (TypeError, ValueError):
        return None
    if total <= 0:
        return None
    return max(0.0, min(1.0, current / total))


def _load_case_groups(available_cases):
    groups = {
        "all": [],
        "standard": [],
        "priority": [],
        "minimum": [],
        "short": [],
    }
    if not os.path.isfile(RUN_SCM_ALL):
        return groups
    try:
        with open(RUN_SCM_ALL, "r", encoding="utf-8") as handle:
            tree = ast.parse(handle.read(), filename=RUN_SCM_ALL)
    except Exception:
        return groups
    mapping = {
        "ALL_CASES": "all",
        "STANDARD_CASES": "standard",
        "PRIORITY_CASES": "priority",
        "MIN_CASES": "minimum",
        "SHORT_CASES": "short",
    }
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id in mapping:
                try:
                    value = ast.literal_eval(node.value)
                except Exception:
                    continue
                if isinstance(value, list):
                    groups[mapping[target.id]] = [
                        case for case in value if case in available_cases
                    ]
    return groups

def _clean_cli_option(value):
    if value is None:
        return ""
    return str(value).strip()


def _build_case_command(case_name, stats_name, cli_options=None):
    stats_value = str(stats_name).strip() if stats_name is not None else DEFAULT_STATS_NAME
    if stats_value.lower() == NO_STATS_NAME:
        stats_arg = NO_STATS_NAME
    else:
        stats_arg = os.path.join("input", "stats", stats_value)
    cmd = [sys.executable, "-u", "run_scripts/run_scm.py", "-stats", stats_arg]
    cli_options = cli_options or {}
    for flag, key in (("-ngrdcol", "ngrdcol"), ("-max_iters", "max_iters"), ("-debug", "debug")):
        value = _clean_cli_option(cli_options.get(key))
        if value:
            cmd.extend([flag, value])
    cmd.append(case_name)
    return " ".join(shlex.quote(str(part)) for part in cmd)


def _start_case_process(case_name, stats_name, overrides, cli_options=None):
    stats_value = str(stats_name).strip() if stats_name is not None else DEFAULT_STATS_NAME
    if stats_value.lower() == NO_STATS_NAME:
        stats_arg = NO_STATS_NAME
    else:
        stats_arg = os.path.join(STATS_DIR, stats_value)
    params_path = _write_temp_namelist(TUNABLE_FILE, overrides.get("tunable"), "clubb_params_")
    flags_path = _write_temp_namelist(FLAGS_FILE, overrides.get("flags"), "clubb_flags_")
    silhs_path = _write_temp_namelist(SILHS_FILE, overrides.get("silhs"), "clubb_silhs_")
    cmd = [sys.executable, "-u", "run_scripts/run_scm.py", "-stats", stats_arg]
    cli_options = cli_options or {}
    for flag, key in (("-ngrdcol", "ngrdcol"), ("-max_iters", "max_iters"), ("-debug", "debug")):
        value = _clean_cli_option(cli_options.get(key))
        if value:
            cmd.extend([flag, value])
    if params_path:
        cmd.extend(["-params", params_path])
    if flags_path:
        cmd.extend(["-flags", flags_path])
    if silhs_path:
        cmd.extend(["-silhs_params", silhs_path])
    cmd.append(case_name)
    log_file = tempfile.NamedTemporaryFile(
        delete=False, prefix="clubb_run_", suffix=".log", dir="/tmp"
    )
    log_path = log_file.name
    proc = subprocess.Popen(
        cmd,
        cwd=REPO_ROOT,
        env=_run_child_env(),
        stdout=log_file,
        stderr=subprocess.STDOUT,
        text=True,
        preexec_fn=_set_child_stack_limit,
    )
    log_file.close()
    temp_files = [path for path in (params_path, flags_path, silhs_path) if path]
    proc_data = {
        "pid": proc.pid,
        "log": log_path,
        "start_time": time.time(),
        "temp_files": temp_files,
    }
    _mark_case_started(case_name, proc, proc_data)
    return proc_data


def _launch_from_queue(running, queued, logs):
    queue = list(queued or [])
    launched = False
    while queue and len(running) < MAX_RUN_PROCS:
        item = queue.pop(0)
        case_name = item.get("case")
        stats_name = item.get("stats") or DEFAULT_STATS_NAME
        overrides = item.get("overrides") or {"flags": {}, "tunable": {}, "silhs": {}}
        cli_options = item.get("cli_options") or {}
        if (
            not case_name
            or case_name in running
            or _is_case_active(case_name)
            or _get_cached_status(case_name) is not None
        ):
            continue
        running[case_name] = _start_case_process(case_name, stats_name, overrides, cli_options)
        logs[case_name] = f"--- Running {case_name} ({stats_name}) ---\n"
        launched = True
    return queue, launched


def build_tab(app):
    cases = _list_cases()
    case_groups = _load_case_groups(cases)
    stats_files = _list_stats_files()
    if DEFAULT_STATS_NAME in stats_files:
        default_stats_name = DEFAULT_STATS_NAME
    elif stats_files:
        default_stats_name = stats_files[0]
    else:
        default_stats_name = DEFAULT_STATS_NAME

    flag_entries = _read_namelist_entries(FLAGS_FILE)
    flag_bools = [entry for entry in flag_entries if _is_bool_value(entry["value"])]
    flag_params = [entry for entry in flag_entries if entry not in flag_bools]

    tunable_entries = _read_namelist_entries(TUNABLE_FILE)
    silhs_entries = _read_namelist_entries(SILHS_FILE)

    flag_names = [entry["name"] for entry in flag_bools]
    all_config_names = (
        flag_names
        + [entry["name"] for entry in flag_params]
        + [entry["name"] for entry in tunable_entries]
        + [entry["name"] for entry in silhs_entries]
    )
    max_name_len = max((len(name) for name in all_config_names), default=16)
    label_width_px = max(170, min(560, int(5 + max_name_len * 7.0)))
    value_width_px = 65
    right_pane_width_px = max(340, min(960, label_width_px + value_width_px + 80))

    param_entries = (
        [{"file": "flags", **entry} for entry in flag_params]
        + [{"file": "tunable", **entry} for entry in tunable_entries]
        + [{"file": "silhs", **entry} for entry in silhs_entries]
    )

    defaults = {
        "flags": {entry["name"]: _is_true(entry["value"]) for entry in flag_bools},
        "params": {
            "flags": {entry["name"]: entry["value"] for entry in flag_params},
            "tunable": {entry["name"]: entry["value"] for entry in tunable_entries},
            "silhs": {entry["name"]: entry["value"] for entry in silhs_entries},
        },
    }
    file_map = {"flags": FLAGS_FILE, "tunable": TUNABLE_FILE, "silhs": SILHS_FILE}

    case_buttons = []
    for case_name in cases:
        case_buttons.append(
            html.Button(
                case_name,
                id={"type": "run-case-button", "name": case_name},
                n_clicks=0,
                style=_case_button_style("#2563eb", False),
            )
        )

    stats_buttons = []
    for stats_name in stats_files:
        stats_buttons.append(
            html.Button(
                stats_name,
                id={"type": "run-stats-button", "name": stats_name},
                n_clicks=0,
                style=_stats_button_style(stats_name == default_stats_name),
            )
        )
    if NO_STATS_NAME not in stats_files:
        stats_buttons.append(
            html.Button(
                "none",
                id={"type": "run-stats-button", "name": NO_STATS_NAME},
                n_clicks=0,
                style=_stats_button_style(NO_STATS_NAME == default_stats_name),
            )
        )

    flag_controls = []
    for entry in flag_bools:
        flag_controls.append(
            html.Div(
                [
                    dcc.Checklist(
                        id={"type": "run-flag", "name": entry["name"]},
                        className="run-flag-checklist",
                        options=[{"label": entry["name"], "value": "on"}],
                        value=["on"] if _is_true(entry["value"]) else [],
                        labelStyle={"display": "inline-flex", "alignItems": "center", "gap": "6px"},
                    )
                ],
                id={"type": "run-flag-container", "name": entry["name"]},
                className="run-param-container run-flag-container",
                style=_field_style(False),
            )
        )

    def _param_input(entry):
        display_value = _normalize_numeric_display(entry["value"])
        return html.Div(
            [
                html.Label(
                    entry["name"],
                    style={"marginRight": "4px", "whiteSpace": "nowrap", "minWidth": f"{label_width_px}px"},
                ),
                dcc.Input(
                    id={"type": "run-param", "file": entry["file"], "name": entry["name"]},
                    type="text",
                    value=display_value,
                    debounce=True,
                    style={"width": f"{value_width_px}px"},
                ),
            ],
            id={"type": "run-param-container", "file": entry["file"], "name": entry["name"]},
            style=_field_style(False),
            className="run-param-container",
        )

    param_sections = []
    if flag_params:
        param_sections.extend(
            [
                html.H4("Flag vals", className="run-settings-heading"),
                html.Div(
                    [
                        _param_input({"file": "flags", **entry})
                        for entry in flag_params
                    ],
                    className="run-param-list",
                ),
            ]
        )
    param_sections.extend(
        [
            html.H4("Flags", className="run-settings-heading"),
            html.Div(flag_controls, className="run-param-list"),
        ]
    )
    if tunable_entries:
        param_sections.extend(
            [
                html.H4("Tunables", className="run-settings-heading"),
                html.Div(
                    [
                        _param_input({"file": "tunable", **entry})
                        for entry in tunable_entries
                    ],
                    className="run-param-list",
                ),
            ]
        )
    if silhs_entries:
        param_sections.extend(
            [
                html.H4("SILHS", className="run-settings-heading"),
                html.Div(
                    [
                        _param_input({"file": "silhs", **entry})
                        for entry in silhs_entries
                    ],
                    className="run-param-list",
                ),
            ]
        )

    run_action_buttons = [
        html.Button(
            "Run selected",
            id="run-button",
            n_clicks=0,
            className="run-button-run-selected",
            style=_run_action_button_style("#111827"),
        ),
        html.Button(
            "Cancel runs",
            id="run-cancel",
            n_clicks=0,
            style=_run_action_button_style("#b91c1c"),
        ),
        html.Button(
            "Clear",
            id="run-clear",
            n_clicks=0,
            style=_run_action_button_style("#374151"),
        ),
    ]

    select_action_buttons = [
        html.Button(
            "Deselect",
            id="run-deselect",
            n_clicks=0,
            style=_run_action_button_style("#6b7280"),
        ),
        html.Button(
            "Select all",
            id={"type": "run-group-button", "name": "all"},
            n_clicks=0,
            disabled=not case_groups["all"],
            className="run-select-group-button",
            style=_run_action_button_style("#111827", disabled=not case_groups["all"]),
        ),
        html.Button(
            "Select standard",
            id={"type": "run-group-button", "name": "standard"},
            n_clicks=0,
            disabled=not case_groups["standard"],
            className="run-select-group-button",
            style=_run_action_button_style("#111827", disabled=not case_groups["standard"]),
        ),
        html.Button(
            "Select priority",
            id={"type": "run-group-button", "name": "priority"},
            n_clicks=0,
            disabled=not case_groups["priority"],
            className="run-select-group-button",
            style=_run_action_button_style("#111827", disabled=not case_groups["priority"]),
        ),
        html.Button(
            "Select minimum",
            id={"type": "run-group-button", "name": "minimum"},
            n_clicks=0,
            disabled=not case_groups["minimum"],
            className="run-select-group-button",
            style=_run_action_button_style("#111827", disabled=not case_groups["minimum"]),
        ),
        html.Button(
            "Select short",
            id={"type": "run-group-button", "name": "short"},
            n_clicks=0,
            disabled=not case_groups["short"],
            className="run-select-group-button",
            style=_run_action_button_style("#111827", disabled=not case_groups["short"]),
        ),
    ]

    layout = html.Div(
        [
            dcc.Store(id="run-defaults", data=defaults),
            dcc.Store(id="run-flag-names", data=flag_names),
            dcc.Store(id="run-param-meta", data=[{"file": entry["file"], "name": entry["name"]} for entry in param_entries]),
            dcc.Store(id="run-file-map", data=file_map),
            dcc.Store(id="run-dirty", data=True),
            dcc.Store(id="run-selected-cases", data=[]),
            dcc.Store(id="run-selected-stats-file", data=default_stats_name),
            dcc.Store(id="run-completed-cases", data=[]),
            dcc.Store(id="run-failed-cases", data=[]),
            dcc.Store(id="run-running-cases", data={}),
            dcc.Store(id="run-queued-cases", data=[]),
            dcc.Store(id="run-case-logs", data={}),
            dcc.Store(id="run-case-commands", data={}),
            dcc.Store(id="run-case-runtimes", data={}),
            dcc.Store(id="run-log-offsets", data={}),
            dcc.Store(id="run-case-order", data=[]),
            dcc.Interval(id="run-interval", interval=500, disabled=True),
            html.Div(
                [
                    html.Div(
                        [
                            html.Div(
                                select_action_buttons,
                                className="run-select-actions",
                                style={
                                    "display": "flex",
                                    "flexWrap": "wrap",
                                    "gap": "8px",
                                    "marginBottom": "6px",
                                },
                            ),
                            html.Div("Cases:", className="run-section-title"),
                            html.Div(case_buttons, className="run-case-buttons"),
                            html.Div(
                                [
                                    html.Div("Stats file:", className="run-section-title"),
                                    html.Div(stats_buttons, className="run-stats-buttons"),
                                ],
                                className="run-stats-section",
                                style={"marginTop": "6px"},
                            ),
                            html.Div(
                                [
                                    html.Div("Optional run args:", className="run-section-title"),
                                    html.Div(
                                        [
                                            dcc.Input(
                                                id="run-opt-ngrdcol",
                                                type="text",
                                                value="",
                                                placeholder="ngrdcol",
                                                style={"width": "130px"},
                                            ),
                                            dcc.Input(
                                                id="run-opt-max-iters",
                                                type="text",
                                                value="",
                                                placeholder="max_iters",
                                                style={"width": "130px"},
                                            ),
                                            dcc.Input(
                                                id="run-opt-debug",
                                                type="text",
                                                value="",
                                                placeholder="debug",
                                                style={"width": "130px"},
                                            ),
                                        ],
                                        style={
                                            "display": "flex",
                                            "flexWrap": "wrap",
                                            "gap": "8px",
                                            "marginTop": "4px",
                                        },
                                    ),
                                ],
                                style={"marginTop": "6px"},
                            ),
                            html.Div(
                                run_action_buttons,
                                className="run-action-buttons",
                                style={
                                    "display": "flex",
                                    "flexWrap": "wrap",
                                    "gap": "8px",
                                    "marginTop": "6px",
                                },
                            ),
                        ],
                        className="run-left-header",
                        style={"marginBottom": "10px"},
                    ),
                    html.Div(
                        id="run-console-container",
                        className="run-console-container",
                        style={"display": "flex", "flexDirection": "column", "gap": "10px"},
                    ),
                ],
                className="run-left-pane",
            ),
            html.Div(id="run-pane-divider", className="run-pane-divider"),
            html.Div(
                param_sections,
                id="run-right-pane",
                className="run-right-pane",
                style={
                    "paddingLeft": "16px",
                    "height": "calc(100vh - 96px)",
                    "minHeight": 0,
                    "overflowY": "auto",
                    "overflowX": "auto",
                },
            ),
        ],
        id="run-tab-layout",
        className="run-tab-layout",
        style={
            "display": "grid",
            "gridTemplateColumns": f"minmax(0,1fr) 8px {right_pane_width_px}px",
            "gap": "16px",
            "padding": "10px",
            "overflowX": "auto",
        },
    )

    @app.callback(
        Output({"type": "run-flag-container", "name": ALL}, "style"),
        Output({"type": "run-param-container", "file": ALL, "name": ALL}, "style"),
        Output("run-dirty", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Input({"type": "run-flag", "name": ALL}, "value"),
        Input({"type": "run-param", "file": ALL, "name": ALL}, "value"),
        State("run-defaults", "data"),
        State("run-flag-names", "data"),
        State("run-param-meta", "data"),
        prevent_initial_call=True,
    )
    def _sync_settings(
        flag_values,
        param_values,
        defaults_data,
        flag_names_data,
        param_meta,
    ):
        changed_flags = []
        for name, values in zip(flag_names_data, flag_values):
            current = bool(values)
            changed_flags.append(current != defaults_data["flags"].get(name))

        changed_params = []
        for meta, value in zip(param_meta, param_values):
            default_value = defaults_data["params"][meta["file"]].get(meta["name"])
            changed_params.append(
                _normalize_numeric_display(value) != _normalize_numeric_display(default_value)
            )

        dirty = any(changed_flags) or any(changed_params)
        flag_styles = [_field_style(changed) for changed in changed_flags]
        param_styles = [_field_style(changed) for changed in changed_params]
        if dirty:
            with RUN_LOCK:
                RUN_STATUS.clear()
                RUN_FINALIZED.clear()
        completed_update = [] if dirty else no_update
        failed_update = [] if dirty else no_update
        return flag_styles, param_styles, dirty, completed_update, failed_update

    @app.callback(
        Output({"type": "run-case-button", "name": ALL}, "style"),
        Input("run-selected-cases", "data"),
        Input("run-completed-cases", "data"),
        Input("run-failed-cases", "data"),
        Input("run-running-cases", "data"),
        Input("run-queued-cases", "data"),
        State({"type": "run-case-button", "name": ALL}, "id"),
    )
    def _update_case_button_styles(
        selected_cases, completed_cases, failed_cases, running_cases, queued_cases, ids
    ):
        if not ids:
            return []
        completed_global, failed_global = _snapshot_status_lists()
        selected = set(selected_cases or [])
        completed = set(completed_cases or []) | set(completed_global)
        failed = set(failed_cases or []) | set(failed_global)
        running = set((running_cases or {}).keys())
        queued = {item.get("case") for item in (queued_cases or []) if item.get("case")}
        styles = []
        for case_id in ids:
            name = case_id.get("name")
            if name in failed:
                color = "#dc2626"
            elif name in completed:
                color = "#16a34a"
            elif name in running or name in queued:
                color = "linear-gradient(135deg, #dc2626 0%, #16a34a 100%)"
            else:
                color = "#2563eb"
            is_selected = name in selected
            styles.append(_case_button_style(color, is_selected))
        return styles

    @app.callback(
        Output({"type": "run-stats-button", "name": ALL}, "style"),
        Input("run-selected-stats-file", "data"),
        State({"type": "run-stats-button", "name": ALL}, "id"),
    )
    def _update_stats_button_styles(selected_stats, ids):
        if not ids:
            return []
        return [_stats_button_style(btn_id.get("name") == selected_stats) for btn_id in ids]

    @app.callback(
        Output("run-selected-stats-file", "data", allow_duplicate=True),
        Input({"type": "run-stats-button", "name": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def _select_stats_file(_n_clicks):
        if not callback_context.triggered:
            return no_update
        trigger_id = callback_context.triggered_id
        if isinstance(trigger_id, dict):
            return trigger_id.get("name", no_update)
        return no_update

    @app.callback(
        Output("run-console-container", "children"),
        Output("run-case-order", "data"),
        Input("run-case-logs", "data"),
        Input("run-selected-cases", "data"),
        Input("run-completed-cases", "data"),
        Input("run-failed-cases", "data"),
        Input("run-running-cases", "data"),
        Input("run-case-runtimes", "data"),
        Input("run-case-commands", "data"),
        State("run-case-order", "data"),
    )
    def _render_consoles(
        case_logs,
        selected_cases,
        completed_cases,
        failed_cases,
        running_cases,
        case_runtimes,
        case_commands,
        case_order,
    ):
        logs = case_logs or {}
        runtimes = case_runtimes or {}
        commands = case_commands or {}
        running_data = running_cases or {}
        now = time.time()
        selected = set(selected_cases or [])
        completed = set(completed_cases or [])
        failed = set(failed_cases or [])
        running = set(running_data.keys())
        candidates = [
            case
            for case in cases
            if case in logs or case in selected or case in completed or case in failed or case in running
        ]
        if not candidates:
            return html.Div(
                "No runs yet.",
                className="run-empty-message",
                style={"fontStyle": "italic", "padding": "8px 0"},
            ), []

        candidate_set = set(candidates)
        prev_order = [name for name in (case_order or []) if name in candidate_set]
        for name in candidates:
            if name not in prev_order:
                prev_order.append(name)
        running_now = [name for name in prev_order if name in running]
        not_running_now = [name for name in prev_order if name not in running]
        display_cases = running_now + not_running_now

        panels = []
        for case_name in display_cases:
            if case_name in running:
                color = "linear-gradient(90deg, #dc2626, #16a34a)"
                status = "running"
            elif case_name in failed:
                color = "#dc2626"
                status = "failed"
            elif case_name in completed:
                color = "#16a34a"
                status = "completed"
            elif case_name in selected:
                color = "#f59e0b"
                status = "selected"
            else:
                color = "#2563eb"
                status = "pending"
            status_style = {"color": color, "fontSize": "12px"}
            if case_name in running:
                status_style = {
                    "backgroundImage": color,
                    "WebkitBackgroundClip": "text",
                    "WebkitTextFillColor": "transparent",
                    "fontSize": "14px",
                }
            else:
                status_style = {"color": color, "fontSize": "14px"}
            if case_name in running:
                start_time = running_data.get(case_name, {}).get("start_time")
                runtime_txt = _format_runtime(now - float(start_time)) if start_time is not None else ""
            else:
                runtime_txt = _format_runtime(runtimes.get(case_name))
            progress = _extract_progress(logs.get(case_name, ""))
            if progress is None:
                if case_name in completed:
                    progress = 1.0
                elif case_name in failed:
                    progress = 1.0
                else:
                    progress = 0.0
            eta_txt = ""
            if case_name in running and progress and progress > 0:
                start_time = running_data.get(case_name, {}).get("start_time")
                if start_time is not None:
                    elapsed = now - float(start_time)
                    remaining = max(0.0, elapsed * (1.0 / progress - 1.0))
                    eta_txt = _format_eta(remaining)
            summary_children = [
                html.Span(case_name, style={"fontWeight": "600", "marginRight": "8px"}),
                html.Span(status, style=status_style),
            ]
            if runtime_txt:
                summary_children.append(
                    html.Span(
                        f"{'elapsed' if case_name in running else 'runtime'}: {runtime_txt}",
                        className="run-muted-text",
                        style={"marginLeft": "10px", "fontSize": "14px"},
                    )
                )
            bar_color = "#16a34a"
            progress_bar = html.Div(
                [
                    html.Span(
                        f"eta {eta_txt}" if eta_txt else "eta --",
                        className="run-muted-text",
                        style={"fontSize": "14px"},
                    ),
                    html.Div(
                        [
                            html.Div(
                                className="run-progress-fill",
                                style={
                                    "width": f"{int(progress * 100)}%",
                                    "height": "16px",
                                    "background": bar_color,
                                    "borderRadius": "999px",
                                    "transition": "width 0.3s ease",
                                }
                            )
                        ],
                        className="run-progress-track",
                        style={
                            "width": "240px",
                            "height": "16px",
                            "borderRadius": "999px",
                            "overflow": "hidden",
                        },
                    ),
                ],
                style={
                    "display": "flex",
                    "alignItems": "center",
                    "gap": "8px",
                    "marginLeft": "auto",
                    "flexShrink": 0,
                },
            )
            summary = html.Summary(
                [
                    html.Div(
                        summary_children,
                        className="run-summary-left",
                        style={
                            "display": "flex",
                            "alignItems": "center",
                            "flexWrap": "wrap",
                            "gap": "6px",
                        },
                    ),
                    progress_bar,
                ],
                className="run-console-summary",
                style={
                    "cursor": "pointer",
                    "display": "flex",
                    "alignItems": "center",
                    "gap": "12px",
                    "fontSize": "14px",
                },
            )
            console_id = f"run-console-{_case_dom_id(case_name)}"
            pre_class = "run-console run-console-active" if case_name in running else "run-console"
            panel_children = [summary]
            command_text = commands.get(case_name)
            if command_text:
                panel_children.append(
                    html.Div(
                        [
                            html.Div(command_text, className="run-command-text"),
                            html.Span("Copy run command", className="run-copy-command-label"),
                            dcc.Clipboard(
                                content=command_text,
                                title="Copy run command",
                                className="run-copy-command-button",
                            ),
                        ],
                        className="run-command-row",
                    )
                )
            panel_children.append(
                html.Pre(
                    logs.get(case_name, "Not run yet."),
                    id=console_id,
                    className=pre_class,
                    style={
                        "padding": "12px",
                        "borderRadius": "6px",
                        "maxHeight": "300px",
                        "overflowY": "auto",
                        "marginTop": "6px",
                    },
                )
            )
            panels.append(
                html.Details(
                    panel_children,
                    open=case_name in running,
                    className="run-console-panel",
                    style={
                        "borderRadius": "6px",
                        "padding": "8px 10px",
                    },
                )
            )
        return panels, display_cases

    @app.callback(
        Output("run-case-logs", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Output("run-case-runtimes", "data", allow_duplicate=True),
        Output("run-log-offsets", "data", allow_duplicate=True),
        Output("run-case-order", "data", allow_duplicate=True),
        Output("run-case-commands", "data", allow_duplicate=True),
        Input("run-clear", "n_clicks"),
        State("run-running-cases", "data"),
        State("run-case-logs", "data"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        State("run-selected-cases", "data"),
        State("run-case-runtimes", "data"),
        State("run-log-offsets", "data"),
        State("run-case-order", "data"),
        State("run-case-commands", "data"),
        prevent_initial_call=True,
    )
    def _clear_non_running_tabs(
        _n_clicks,
        running_cases,
        case_logs,
        completed_cases,
        failed_cases,
        selected_cases,
        case_runtimes,
        log_offsets,
        case_order,
        case_commands,
    ):
        running = set((running_cases or {}).keys())
        if not running and not (case_logs or completed_cases or failed_cases or selected_cases):
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update

        logs_in = dict(case_logs or {})
        completed_in = list(completed_cases or [])
        failed_in = list(failed_cases or [])
        selected_in = list(selected_cases or [])
        runtimes_in = dict(case_runtimes or {})
        offsets_in = dict(log_offsets or {})
        order_in = list(case_order or [])
        commands_in = dict(case_commands or {})

        logs_out = {k: v for k, v in logs_in.items() if k in running}
        completed_out = [k for k in completed_in if k in running]
        failed_out = [k for k in failed_in if k in running]
        runtimes_out = {k: v for k, v in runtimes_in.items() if k in running}
        order_out = [k for k in order_in if k in running]
        commands_out = {k: v for k, v in commands_in.items() if k in running}

        running_logs = {
            v.get("log")
            for v in (running_cases or {}).values()
            if isinstance(v, dict) and v.get("log")
        }
        offsets_out = {}
        for key, value in offsets_in.items():
            if key in running or key in running_logs:
                offsets_out[key] = value

        # Also clear global run status cache for any case we removed from the UI,
        # otherwise poll callbacks repopulate completed/failed lists on next tick.
        known_cases = (
            set(logs_in.keys())
            | set(completed_in)
            | set(failed_in)
            | set(selected_in)
            | set(runtimes_in.keys())
            | set(order_in)
            | set(commands_in.keys())
        )
        removed_cases = [case_name for case_name in known_cases if case_name not in running]
        for case_name in removed_cases:
            _clear_case_status(case_name)

        return (
            logs_out if logs_out != logs_in else no_update,
            completed_out if completed_out != completed_in else no_update,
            failed_out if failed_out != failed_in else no_update,
            runtimes_out if runtimes_out != runtimes_in else no_update,
            offsets_out if offsets_out != offsets_in else no_update,
            order_out if order_out != order_in else no_update,
            commands_out if commands_out != commands_in else no_update,
        )

    @app.callback(
        Output("run-selected-cases", "data", allow_duplicate=True),
        Input({"type": "run-case-button", "name": ALL}, "n_clicks"),
        State("run-selected-cases", "data"),
        State("run-running-cases", "data"),
        prevent_initial_call=True,
    )
    def _select_case(_n_clicks, selected_cases, running_cases):
        if not callback_context.triggered:
            return no_update
        trigger_id = callback_context.triggered_id
        case_name = trigger_id.get("name") if isinstance(trigger_id, dict) else None
        if case_name is None:
            trigger = callback_context.triggered[0]["prop_id"].split(".")[0]
            try:
                trigger_id = json.loads(trigger)
            except Exception:
                trigger_id = None
            case_name = trigger_id.get("name") if trigger_id else None
        if not case_name:
            return no_update
        if case_name in (running_cases or {}):
            return no_update
        selected = list(selected_cases or [])
        if case_name in selected:
            selected.remove(case_name)
        else:
            selected.append(case_name)
        return selected

    @app.callback(
        Output("run-case-logs", "data", allow_duplicate=True),
        Output("run-running-cases", "data", allow_duplicate=True),
        Output("run-queued-cases", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Output("run-interval", "disabled", allow_duplicate=True),
        Output("run-interval", "n_intervals"),
        Output("run-case-runtimes", "data", allow_duplicate=True),
        Output("run-case-commands", "data", allow_duplicate=True),
        Input("run-button", "n_clicks"),
        State("run-selected-cases", "data"),
        State("run-running-cases", "data"),
        State("run-queued-cases", "data"),
        State("run-case-logs", "data"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        State("run-selected-stats-file", "data"),
        State("run-opt-ngrdcol", "value"),
        State("run-opt-max-iters", "value"),
        State("run-opt-debug", "value"),
        State("run-case-runtimes", "data"),
        State("run-case-commands", "data"),
        State({"type": "run-flag", "name": ALL}, "value"),
        State({"type": "run-param", "file": ALL, "name": ALL}, "value"),
        State("run-defaults", "data"),
        State("run-flag-names", "data"),
        State("run-param-meta", "data"),
        prevent_initial_call=True,
    )
    def _run_selected_cases(
        _selected_clicks,
        selected_cases,
        running_cases,
        queued_cases,
        case_logs,
        completed_cases,
        failed_cases,
        selected_stats,
        opt_ngrdcol,
        opt_max_iters,
        opt_debug,
        case_runtimes,
        case_commands,
        flag_values,
        param_values,
        defaults_data,
        flag_names_data,
        param_meta,
    ):
        trigger_id = callback_context.triggered_id
        if trigger_id == "run-button":
            cases_to_run = list(selected_cases or [])
        else:
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
            )
        if not cases_to_run:
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
            )

        running = dict(running_cases or {})
        queued = list(queued_cases or [])
        queued_names = {item.get("case") for item in queued}
        logs = dict(case_logs or {})
        completed = list(completed_cases or [])
        failed = list(failed_cases or [])
        runtimes = dict(case_runtimes or {})
        commands = dict(case_commands or {})
        stats_name = selected_stats or DEFAULT_STATS_NAME
        overrides = _build_override_updates(
            flag_values, param_values, defaults_data, flag_names_data, param_meta
        )
        cli_options = {}
        ngrdcol_value = _clean_cli_option(opt_ngrdcol)
        max_iters_value = _clean_cli_option(opt_max_iters)
        debug_value = _clean_cli_option(opt_debug)
        if ngrdcol_value:
            cli_options["ngrdcol"] = ngrdcol_value
        if max_iters_value:
            cli_options["max_iters"] = max_iters_value
        if debug_value:
            cli_options["debug"] = debug_value

        for case_name in cases_to_run:
            if case_name in running or case_name in queued_names or _is_case_active(case_name):
                continue
            _clear_case_status(case_name)
            if case_name in completed:
                completed.remove(case_name)
            if case_name in failed:
                failed.remove(case_name)
            runtimes.pop(case_name, None)
            overrides_copy = {key: dict(value) for key, value in overrides.items()}
            cli_options_copy = dict(cli_options)
            queued.append(
                {
                    "case": case_name,
                    "stats": stats_name,
                    "overrides": overrides_copy,
                    "cli_options": cli_options_copy,
                }
            )
            queued_names.add(case_name)
            logs[case_name] = f"--- Queued {case_name} ({stats_name}) ---\n"
            commands[case_name] = _build_case_command(case_name, stats_name, cli_options_copy)

        queued, started_any = _launch_from_queue(running, queued, logs)
        interval_disabled = not bool(running or queued)
        n_intervals = 0 if (started_any or queued) else no_update
        return logs, running, queued, completed, failed, interval_disabled, n_intervals, runtimes, commands

    @app.callback(
        Output("run-selected-cases", "data", allow_duplicate=True),
        Input({"type": "run-group-button", "name": ALL}, "n_clicks"),
        State("run-selected-cases", "data"),
        prevent_initial_call=True,
    )
    def _select_case_group(_n_clicks, selected_cases):
        trigger_id = callback_context.triggered_id
        if not isinstance(trigger_id, dict):
            return no_update
        group_name = trigger_id.get("name")
        if not group_name:
            return no_update
        group_cases = case_groups.get(group_name, [])
        if not group_cases:
            return no_update
        return list(group_cases)

    @app.callback(
        Output("run-selected-cases", "data", allow_duplicate=True),
        Input("run-deselect", "n_clicks"),
        prevent_initial_call=True,
    )
    def _deselect_all(_n_clicks):
        return []

    @app.callback(
        Output("run-running-cases", "data", allow_duplicate=True),
        Output("run-queued-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Output("run-case-logs", "data", allow_duplicate=True),
        Output("run-interval", "disabled", allow_duplicate=True),
        Output("run-case-runtimes", "data", allow_duplicate=True),
        Input("run-cancel", "n_clicks"),
        State("run-running-cases", "data"),
        State("run-queued-cases", "data"),
        State("run-failed-cases", "data"),
        State("run-case-logs", "data"),
        State("run-case-runtimes", "data"),
        prevent_initial_call=True,
    )
    def _cancel_runs(_n_clicks, running_cases, queued_cases, failed_cases, case_logs, case_runtimes):
        with RUN_STREAM_LOCK:
            running = dict(running_cases or {})
            for case_name, proc_data in _snapshot_active_cases().items():
                if case_name not in running:
                    running[case_name] = proc_data
            queued = list(queued_cases or [])
            if not running and not queued:
                return no_update, no_update, no_update, no_update, True, no_update
            failed = list(failed_cases or [])
            logs = dict(case_logs or {})
            runtimes = dict(case_runtimes or {})
            now = time.time()
            for case_name, proc_data in running.items():
                pid = proc_data.get("pid")
                log_path = proc_data.get("log")
                _cleanup_temp_files(proc_data.get("temp_files"))
                proc = _get_proc(pid)
                if proc is not None:
                    proc.terminate()
                    if proc.poll() is None:
                        proc.kill()
                _record_case_finish(case_name, pid, 1)
                _cleanup_log_file(log_path)
                if case_name not in failed:
                    failed.append(case_name)
                runtime_secs = now - float(proc_data.get("start_time", now))
                runtimes[case_name] = runtime_secs
                runtime_txt = _format_runtime(runtime_secs)
                existing = logs.get(case_name, "")
                if existing and not existing.endswith("\n"):
                    existing += "\n"
                logs[case_name] = _append_log_tail(
                    existing, f"--- Cancelled (runtime: {runtime_txt}) ---\n"
                )
            for item in queued:
                case_name = item.get("case")
                if not case_name:
                    continue
                _record_case_finish(case_name, None, 1)
                if case_name not in failed:
                    failed.append(case_name)
                existing = logs.get(case_name, "")
                if existing and not existing.endswith("\n"):
                    existing += "\n"
                logs[case_name] = _append_log_tail(
                    existing, "--- Cancelled (never started) ---\n"
                )
            _, failed_global = _snapshot_status_lists()
            return {}, [], failed_global, logs, True, runtimes

    @app.callback(
        Output("run-case-logs", "data", allow_duplicate=True),
        Output("run-interval", "disabled", allow_duplicate=True),
        Output("run-running-cases", "data", allow_duplicate=True),
        Output("run-queued-cases", "data", allow_duplicate=True),
        Output("run-completed-cases", "data", allow_duplicate=True),
        Output("run-failed-cases", "data", allow_duplicate=True),
        Output("run-selected-cases", "data", allow_duplicate=True),
        Output("run-case-runtimes", "data", allow_duplicate=True),
        Output("run-log-offsets", "data", allow_duplicate=True),
        Input("run-interval", "n_intervals"),
        State("run-running-cases", "data"),
        State("run-queued-cases", "data"),
        State("run-case-logs", "data"),
        State("run-completed-cases", "data"),
        State("run-failed-cases", "data"),
        State("run-selected-cases", "data"),
        State("run-case-runtimes", "data"),
        State("run-log-offsets", "data"),
        prevent_initial_call=True,
    )
    def _stream_output(
        _tick,
        running_cases,
        queued_cases,
        case_logs,
        completed_cases,
        failed_cases,
        selected_cases,
        case_runtimes,
        log_offsets,
    ):
        if not RUN_STREAM_LOCK.acquire(blocking=False):
            return no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update, no_update

        try:
            running = dict(running_cases or {})
            queued = list(queued_cases or [])
            if not running and not queued:
                return no_update, True, no_update, no_update, no_update, no_update, no_update, no_update, no_update
            active_cases = _snapshot_active_cases()
            active_names = set(active_cases.keys())
            for case_name, proc_data in active_cases.items():
                if case_name not in running:
                    running[case_name] = proc_data
            if queued:
                queued = [
                    item
                    for item in queued
                    if item.get("case")
                    and item.get("case") not in active_names
                    and _get_cached_status(item.get("case")) is None
                ]
            running_in = dict(running_cases or {})
            queued_in = list(queued_cases or [])
            completed_in = list(completed_cases or [])
            failed_in = list(failed_cases or [])
            selected_in = list(selected_cases or [])
            runtimes_in = dict(case_runtimes or {})
            logs = dict(case_logs or {})
            completed = list(completed_in)
            failed = list(failed_in)
            selected = list(selected_in)
            runtimes = dict(runtimes_in)
            offsets = dict(log_offsets or {})
            finished = []
            now = time.time()
            logs_changed = False
            for case_name, proc_data in running.items():
                pid = proc_data.get("pid")
                log_path = proc_data.get("log")
                offset_key = log_path or case_name
                chunk, new_offset = _read_log_increment(log_path, int(offsets.get(offset_key, 0)))
                offsets[offset_key] = new_offset
                if chunk:
                    logs[case_name] = _append_log_tail(logs.get(case_name, ""), chunk)
                    logs_changed = True
                proc = _get_proc(pid)
                if proc is None:
                    cached_status = _get_cached_status(case_name)
                    if cached_status is None:
                        misses = int(proc_data.get("missing_polls", 0)) + 1
                        proc_data["missing_polls"] = misses
                        running[case_name] = proc_data
                        if misses >= 5:
                            runtime_secs = now - float(proc_data.get("start_time", now))
                            finished.append((case_name, 1, runtime_secs))
                        continue
                    runtime_secs = now - float(proc_data.get("start_time", now))
                    finished.append((case_name, cached_status, runtime_secs))
                    continue
                status = proc.poll()
                if status is not None:
                    runtime_secs = now - float(proc_data.get("start_time", now))
                    finished.append((case_name, status, runtime_secs))
            for case_name, status, runtime_secs in finished:
                proc_data = running.pop(case_name, None) or {}
                runtimes[case_name] = runtime_secs
                runtime_txt = _format_runtime(runtime_secs)
                pid = proc_data.get("pid")
                log_path = proc_data.get("log")
                _cleanup_temp_files(proc_data.get("temp_files"))
                offset_key = log_path or case_name
                chunk, _ = _read_log_increment(log_path, int(offsets.get(offset_key, 0)))
                if chunk:
                    logs[case_name] = _append_log_tail(logs.get(case_name, ""), chunk)
                    logs_changed = True
                offsets.pop(offset_key, None)
                already_finalized = _record_case_finish(case_name, pid, status)
                _cleanup_log_file(log_path)
                if not already_finalized:
                    if status == 0:
                        if case_name not in completed:
                            completed.append(case_name)
                        if case_name in failed:
                            failed.remove(case_name)
                    else:
                        if case_name not in failed:
                            failed.append(case_name)
                        if case_name in completed:
                            completed.remove(case_name)
                    existing = logs.get(case_name, "")
                    if existing and not existing.endswith("\n"):
                        existing += "\n"
                    result_txt = "completed" if status == 0 else f"failed (exit {status})"
                    logs[case_name] = _append_log_tail(
                        existing, f"--- {result_txt}; runtime: {runtime_txt} ---\n"
                    )
                    logs_changed = True
            running_before = dict(running)
            queued_before = list(queued)
            queued, _ = _launch_from_queue(running, queued, logs)
            if running != running_before or queued != queued_before:
                logs_changed = True
            interval_disabled = not bool(running or queued)
            logs_out = logs if logs_changed else no_update
            running_out = running if running != running_in else no_update
            queued_out = queued if queued != queued_in else no_update
            completed_global, failed_global = _snapshot_status_lists()
            completed_out = completed_global if completed_global != completed_in else no_update
            failed_out = failed_global if failed_global != failed_in else no_update
            selected_out = selected if selected != selected_in else no_update
            runtimes_out = runtimes if runtimes != runtimes_in else no_update
            offsets_out = offsets if offsets != (log_offsets or {}) else no_update
            return (
                logs_out,
                interval_disabled,
                running_out,
                queued_out,
                completed_out,
                failed_out,
                selected_out,
                runtimes_out,
                offsets_out,
            )
        finally:
            RUN_STREAM_LOCK.release()

    return dcc.Tab(label="Run", value="run", children=layout)
