#!/usr/bin/env bash
#
# Run the repository's pytest suites without accidentally launching expensive
# SCM, GPU, or compiled-API tests.
#
# Usage:
#   tests/run_pytests.sh -unit [pytest options]
#   tests/run_pytests.sh -dash [pytest options]
#   tests/run_pytests.sh -api [pytest options]
#   tests/run_pytests.sh -all [pytest options]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODE="${1:--unit}"
shift || true
cd "$REPO_ROOT"

run_unit() {
    python3 -m pytest utilities/pytests tuner/pytests "$@"
}

run_dash() {
    local dash_python="${CLUBB_DASH_PYTHON:-$REPO_ROOT/.venv-dash/bin/python}"
    if [[ ! -x "$dash_python" ]]; then
        echo "Dash Python was not found: $dash_python" >&2
        echo "Create the Dash environment with ./launch_dashboard.sh or set CLUBB_DASH_PYTHON." >&2
        return 2
    fi
    if ! "$dash_python" -c 'import pytest' >/dev/null 2>&1; then
        echo "Pytest is not installed in the Dash environment: $dash_python" >&2
        echo "Install Dash dependencies with: $dash_python -m pip install -r dash_app/requirements.txt" >&2
        return 2
    fi
    "$dash_python" -m pytest dash_app "$@"
}

run_api() {
    bash clubb_python_api/run_pytests.sh "$@"
}

case "$MODE" in
    -unit|unit)
        run_unit "$@"
        ;;
    -dash|dash)
        run_dash "$@"
        ;;
    -api|api)
        run_api "$@"
        ;;
    -all|all)
        run_unit "$@"
        run_dash "$@"
        run_api "$@"
        ;;
    -h|--help|help)
        sed -n '3,10p' "$0"
        ;;
    *)
        echo "Unknown suite: $MODE" >&2
        echo "Use -unit, -dash, -api, or -all. Run $0 --help for details." >&2
        exit 2
        ;;
esac
