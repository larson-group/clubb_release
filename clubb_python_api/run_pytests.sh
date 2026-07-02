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

F2PY_DIR="${CLUBB_F2PY_DIR:-$REPO_ROOT/install/latest/python}"
if [[ ! -d "$F2PY_DIR" ]]; then
  echo "Python runtime directory not found: $F2PY_DIR" >&2
  echo "Rebuild with ./compile.py -python, or set CLUBB_F2PY_DIR." >&2
  exit 1
fi
if ! compgen -G "$F2PY_DIR/clubb_f2py*.so" > /dev/null; then
  echo "No clubb_f2py extension found in: $F2PY_DIR" >&2
  echo "Rebuild with ./compile.py -python, or set CLUBB_F2PY_DIR." >&2
  exit 1
fi
if [[ ! -f "$F2PY_DIR/libclubb_f2py_backend.so" ]]; then
  echo "libclubb_f2py_backend.so not found in: $F2PY_DIR" >&2
  echo "Rebuild with ./compile.py -python, or set CLUBB_F2PY_DIR." >&2
  exit 1
fi
if [[ ! -d "$F2PY_DIR/clubb_python" ]]; then
  echo "clubb_python package not found in: $F2PY_DIR" >&2
  echo "Rebuild with ./compile.py -python, or set CLUBB_F2PY_DIR." >&2
  exit 1
fi

cd "$SCRIPT_DIR"
export PYTHONPATH="$F2PY_DIR:$REPO_ROOT:$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}"

python3 -m pytest tests/ "$@"
