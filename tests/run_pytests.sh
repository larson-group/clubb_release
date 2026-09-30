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
    local python_bin="${CLUBB_PYTHON:-}"
    if [[ -z "$python_bin" ]]; then
        python3 "$REPO_ROOT/utilities/setup_python_venv.py" python
        python_bin="${CLUBB_PYTHON_VENV:-$REPO_ROOT/.venv-python}/bin/python"
    fi
    "$python_bin" -m pytest utilities/pytests tuner/pytests "$@"
}

run_dash() {
    python3 "$REPO_ROOT/utilities/setup_python_venv.py" dash
    local dash_python="${CLUBB_PYTHON_VENV:-$REPO_ROOT/.venv-python}/bin/python"
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
