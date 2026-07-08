#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$script_dir"
cd "$repo_root"

usage() {
  cat <<'EOF'
Usage: ./launch_dashboard.sh [dash-app-args...]

Creates or reuses a local Dash virtualenv, installs dash_app/requirements.txt,
then launches dash_app/app.py. Arguments are passed through to the Dash app.

Environment:
  CLUBB_DASH_VENV  Virtualenv path. Default: .venv-dash
  PYTHON           Python executable used to create the venv.

Examples:
  ./launch_dashboard.sh
  ./launch_dashboard.sh --port 8060 -debug
  CLUBB_DASH_VENV=.venv-dash-dev ./launch_dashboard.sh
EOF
}

if [[ "${1:-}" == "--launcher-help" ]]; then
  usage
  exit 0
fi

requirements="$repo_root/dash_app/requirements.txt"
app="$repo_root/dash_app/app.py"
venv_dir="${CLUBB_DASH_VENV:-$repo_root/.venv-dash}"

if [[ -n "${PYTHON:-}" ]]; then
  python_cmd="$PYTHON"
elif command -v python3 >/dev/null 2>&1; then
  python_cmd="python3"
else
  python_cmd="python"
fi

if ! command -v "$python_cmd" >/dev/null 2>&1; then
  echo "Could not find Python executable: $python_cmd" >&2
  echo "Set PYTHON=/path/to/python and try again." >&2
  exit 1
fi

if [[ ! -f "$requirements" ]]; then
  echo "Missing Dash requirements file: $requirements" >&2
  exit 1
fi

if [[ ! -x "$venv_dir/bin/python" ]]; then
  echo "Creating Dash virtualenv: $venv_dir"
  "$python_cmd" -m venv "$venv_dir"
fi

venv_python="$venv_dir/bin/python"

echo "Installing Dash dependencies from dash_app/requirements.txt"
"$venv_python" -m pip install -r "$requirements"

echo "Starting CLUBB Dash app"
exec "$venv_python" "$app" "$@"
