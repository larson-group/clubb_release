#!/bin/bash
#
# Run the CLUBB Python API pytest suite.
#
# Usage:
#   bash clubb_python_api/run_pytests.sh
#     Run the full suite under clubb_python_api/tests/.
#
#   bash clubb_python_api/run_pytests.sh -q
#     Run with quieter pytest output.
#
#   bash clubb_python_api/run_pytests.sh -v
#     Run with more verbose pytest output.
#
#   bash clubb_python_api/run_pytests.sh tests/test_python_port_api_coverage.py -v
#     Run one specific test file with verbose output.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$SCRIPT_DIR"
export PYTHONPATH="$REPO_ROOT:$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}"

python3 -m pytest tests/ "$@"
