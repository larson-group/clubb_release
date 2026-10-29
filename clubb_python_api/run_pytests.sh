#!/bin/bash
#
# Run the CLUBB Python API pytest suite.
#
# Usage:
#   bash clubb_python_api/run_pytests.sh
#     Run the tests under clubb_python_api/auto_llm_generated_tests/pytests/.
#
#   bash clubb_python_api/run_pytests.sh -q
#     Run with quieter pytest output.
#
#   bash clubb_python_api/run_pytests.sh -v
#     Run with more verbose pytest output.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VENV_DIR="${CLUBB_PYTHON_VENV:-$REPO_ROOT/.venv-python}"
if [[ "$VENV_DIR" != /* ]]; then
  VENV_DIR="$REPO_ROOT/$VENV_DIR"
fi
PYTHON_BIN="${CLUBB_PYTHON:-}"
if [[ -z "$PYTHON_BIN" ]]; then
  python3 "$REPO_ROOT/utilities/setup_python_venv.py" python
  PYTHON_BIN="$VENV_DIR/bin/python"
fi
if [[ "$PYTHON_BIN" == */* && "$PYTHON_BIN" != /* ]]; then
  PYTHON_BIN="$REPO_ROOT/$PYTHON_BIN"
fi

F2PY_DIR="${CLUBB_F2PY_DIR:-$REPO_ROOT/install/latest/python}"
if [[ "$F2PY_DIR" != /* ]]; then
  F2PY_DIR="$REPO_ROOT/$F2PY_DIR"
fi
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
if ! compgen -G "$F2PY_DIR/libclubb_f2py_backend.*" > /dev/null; then
  echo "libclubb_f2py_backend library not found in: $F2PY_DIR" >&2
  echo "Rebuild with ./compile.py -python, or set CLUBB_F2PY_DIR." >&2
  exit 1
fi
if [[ ! -d "$F2PY_DIR/clubb_python" ]]; then
  echo "clubb_python package not found in: $F2PY_DIR" >&2
  echo "Rebuild with ./compile.py -python, or set CLUBB_F2PY_DIR." >&2
  exit 1
fi

cd "$REPO_ROOT"
export PYTHONPATH="$F2PY_DIR:$REPO_ROOT:$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}"
PYTEST_ARGS=()
has_test_target=false
for arg in "$@"; do
  target="${arg%%::*}"
  # Accept suite-relative targets while running file-reading checks from the repo root.
  if [[ "$arg" != -* && ! -e "$target" && -e "$SCRIPT_DIR/$target" ]]; then
    arg="$SCRIPT_DIR/$arg"
    target="$SCRIPT_DIR/$target"
  fi
  PYTEST_ARGS+=("$arg")
  if [[ "$arg" != -* && ( -f "$target" || -d "$target" ) ]]; then
    has_test_target=true
  fi
done
if [[ "$has_test_target" == false ]]; then
  PYTEST_ARGS=("$SCRIPT_DIR/auto_llm_generated_tests/pytests/" ${PYTEST_ARGS[@]+"${PYTEST_ARGS[@]}"})
fi

# Load the selected extension before pytest can add source paths to sys.path.
"$PYTHON_BIN" -c '
import sys
sys.path.insert(0, sys.argv[1])
import netCDF4
import clubb_f2py
import pytest
raise SystemExit(pytest.main(sys.argv[2:]))
' "$F2PY_DIR" ${PYTEST_ARGS[@]+"${PYTEST_ARGS[@]}"}
