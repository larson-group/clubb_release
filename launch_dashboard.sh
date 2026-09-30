#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$script_dir"
cd "$repo_root"

usage() {
  cat <<'EOF'
Usage: ./launch_dashboard.sh [dash-app-args...]

Creates or reuses the shared CLUBB virtualenv with uv, installs Dash requirements,
then launches the foreground dashboard manager. Arguments are passed through
to the Dash app. The manager owns the runtime broker and restarts Dash every
10 seconds for up to 5 minutes after a crash.

Environment:
  CLUBB_PYTHON_VENV  Virtualenv path. Default: .venv-python
  PYTHON           Python 3.12+ executable used to create the venv.

The launcher keeps Dash stdout/stderr live in this terminal and also writes
them to the private rotating runtime log listed in connection.json.

Examples:
  ./launch_dashboard.sh
  ./launch_dashboard.sh --port 8060 -debug
  ./launch_dashboard.sh --restart-runtime
  CLUBB_PYTHON_VENV=.venv-clubb-dev ./launch_dashboard.sh
EOF
}

if [[ "${1:-}" == "--launcher-help" ]]; then
  usage
  exit 0
fi

venv_dir="${CLUBB_PYTHON_VENV:-$repo_root/.venv-python}"
python3 "$repo_root/utilities/setup_python_venv.py" dash

venv_python="$venv_dir/bin/python"

echo "Starting CLUBB Dash manager"
dash_log="$($venv_python -m dash_app.shared.runtime_logging prepare --repo "$repo_root")"

# Keep both streams live in the launching terminal while forwarding them to
# the private rotating runtime log. Separate relays preserve stdout/stderr.
exec > >("$venv_python" -m dash_app.shared.runtime_logging relay --path "$dash_log" --stream stdout)
exec 2> >("$venv_python" -m dash_app.shared.runtime_logging relay --path "$dash_log" --stream stderr)
exec "$venv_python" -m dash_app.manager "$@"
