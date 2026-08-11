#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$script_dir"
cd "$repo_root"

usage() {
  cat <<'EOF'
Usage: ./launch_dashboard.sh [dash-app-args...]

Creates or reuses a local Dash virtualenv, installs dash_app/requirements.txt,
then launches the foreground dashboard manager. Arguments are passed through
to the Dash app. The manager owns the runtime broker and restarts Dash every
10 seconds for up to 5 minutes after a crash.

Environment:
  CLUBB_DASH_VENV  Virtualenv path. Default: .venv-dash
  PYTHON           Python executable used to create the venv.

The launcher keeps Dash stdout/stderr live in this terminal and also writes
them to the private rotating runtime log listed in connection.json.

Examples:
  ./launch_dashboard.sh
  ./launch_dashboard.sh --port 8060 -debug
  ./launch_dashboard.sh --restart-runtime
  CLUBB_DASH_VENV=.venv-dash-dev ./launch_dashboard.sh
EOF
}

if [[ "${1:-}" == "--launcher-help" ]]; then
  usage
  exit 0
fi

requirements="$repo_root/dash_app/requirements.txt"
venv_dir="${CLUBB_DASH_VENV:-$repo_root/.venv-dash}"

python_supports_dashboard() {
  "$1" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 10) else 1)' >/dev/null 2>&1
}

if [[ -n "${PYTHON:-}" ]]; then
  python_cmd="$PYTHON"
elif [[ -z "${python_cmd:-}" ]]; then
  python_cmd=""
  for candidate in python3.13 python3.12 python3.11 python3.10 python3 python; do
    if command -v "$candidate" >/dev/null 2>&1 && python_supports_dashboard "$candidate"; then
      python_cmd="$candidate"
      break
    fi
  done
fi

if [[ -z "${python_cmd:-}" ]] || ! command -v "$python_cmd" >/dev/null 2>&1; then
  echo "Could not find Python executable: $python_cmd" >&2
  echo "Set PYTHON=/path/to/python3.10-or-newer and try again." >&2
  exit 1
fi

if ! python_supports_dashboard "$python_cmd"; then
  echo "Dash requires Python 3.10 or newer; selected Python is: $("$python_cmd" --version 2>&1)" >&2
  echo "Set PYTHON=/path/to/python3.10-or-newer and try again." >&2
  exit 1
fi

if [[ ! -f "$requirements" ]]; then
  echo "Missing Dash requirements file: $requirements" >&2
  exit 1
fi

if [[ -x "$venv_dir/bin/python" ]] && ! python_supports_dashboard "$venv_dir/bin/python"; then
  echo "Existing Dash virtualenv uses unsupported Python: $("$venv_dir/bin/python" --version 2>&1)" >&2
  echo "Remove $venv_dir or set CLUBB_DASH_VENV to a new path, then relaunch." >&2
  exit 1
fi

if [[ ! -x "$venv_dir/bin/python" ]]; then
  echo "Creating Dash virtualenv: $venv_dir"
  "$python_cmd" -m venv "$venv_dir"
fi

venv_python="$venv_dir/bin/python"

echo "Installing Dash dependencies from dash_app/requirements.txt"
"$venv_python" -m pip install -r "$requirements"

echo "Starting CLUBB Dash manager"
dash_log="$($venv_python -m dash_app.shared.runtime_logging prepare --repo "$repo_root")"

# Keep both streams live in the launching terminal while forwarding them to
# the private rotating runtime log. Separate relays preserve stdout/stderr.
exec > >("$venv_python" -m dash_app.shared.runtime_logging relay --path "$dash_log" --stream stdout)
exec 2> >("$venv_python" -m dash_app.shared.runtime_logging relay --path "$dash_log" --stream stderr)
exec "$venv_python" -m dash_app.manager "$@"
