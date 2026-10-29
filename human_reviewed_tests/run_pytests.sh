#!/usr/bin/env bash
# Run the auto-LLM-generated pytest suites in their owned environments.
# Usage: human_reviewed_tests/run_pytests.sh [-unit|-dash|-api|-jax|-jax_f2py_unit_tests|-all] [pytest options]
# Each mode runs one component's auto_llm_generated_tests/pytests folder.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODE=-unit
case "${1:-}" in
    -unit|-dash|-api|-jax|-jax_f2py_unit_tests|-all) MODE="$1"; shift ;;
    -h|-help) sed -n '2,4p' "$0"; exit 0 ;;
esac
PYTEST_ARGS=("$@")
cd "$REPO_ROOT"

run_pytest() {
    local python_bin="$1"
    shift
    "$python_bin" -m pytest "$@" ${PYTEST_ARGS[@]+"${PYTEST_ARGS[@]}"}
}

run_unit() {
    local python_bin="${CLUBB_PYTHON:-}"
    if [[ -z "$python_bin" ]]; then
        python3 utilities/setup_python_venv.py python
        local venv="${CLUBB_PYTHON_VENV:-$REPO_ROOT/.venv-python}"
        [[ "$venv" == /* ]] || venv="$REPO_ROOT/$venv"
        python_bin="$venv/bin/python"
    fi
    run_pytest "$python_bin" input/auto_llm_generated_tests/pytests utilities/auto_llm_generated_tests/pytests \
        tuner/auto_llm_generated_tests/pytests human_reviewed_tests/auto_llm_generated_tests/pytests
}

run_dash() {
    local python_bin="${CLUBB_PYTHON:-}"
    if [[ -z "$python_bin" ]]; then
        python3 utilities/setup_python_venv.py dash
        local venv="${CLUBB_PYTHON_VENV:-$REPO_ROOT/.venv-python}"
        [[ "$venv" == /* ]] || venv="$REPO_ROOT/$venv"
        python_bin="$venv/bin/python"
    fi
    run_pytest "$python_bin" dash_app/auto_llm_generated_tests/pytests
}

run_api() {
    bash clubb_python_api/run_pytests.sh ${PYTEST_ARGS[@]+"${PYTEST_ARGS[@]}"}
}

# Prepare the CPU JAX environment and set JAX_PYTHON to its interpreter.
prepare_jax() {
    python3 clubb_jax/run_jax.py -profile=cpu -init_env
    local venv="${CLUBB_JAX_VENV:-$REPO_ROOT/.venv-jax}"
    [[ "$venv" == /* ]] || venv="$REPO_ROOT/$venv"
    JAX_PYTHON="$venv/bin/python"
}

# JAX checks that need no compiled Fortran. The f2py_unit_tests folder needs the
# Python API build and runs with -jax_f2py_unit_tests instead.
run_jax() {
    prepare_jax
    local PYTEST_ARGS=(--ignore-glob='*/f2py_unit_tests' ${PYTEST_ARGS[@]+"${PYTEST_ARGS[@]}"})
    run_pytest "$JAX_PYTHON" clubb_jax/auto_llm_generated_tests/pytests
}

# JAX checks against the same Fortran routines called through the Python API build.
run_jax_f2py_unit_tests() {
    local f2py_dir="${CLUBB_F2PY_DIR:-$REPO_ROOT/install/latest/python}"
    [[ "$f2py_dir" == /* ]] || f2py_dir="$REPO_ROOT/$f2py_dir"
    if ! compgen -G "$f2py_dir/clubb_f2py*.so" > /dev/null || [[ ! -d "$f2py_dir/clubb_python" ]]; then
        echo "No clubb_f2py build with clubb_python in: $f2py_dir" >&2
        echo "Rebuild with ./compile.py -debug -python, or set CLUBB_F2PY_DIR." >&2
        return 1
    fi
    prepare_jax
    PYTHONPATH="$f2py_dir:$REPO_ROOT${PYTHONPATH:+:$PYTHONPATH}" run_pytest "$JAX_PYTHON" clubb_jax/auto_llm_generated_tests/pytests/f2py_unit_tests
}

case "$MODE" in
    -unit) run_unit ;;
    -dash) run_dash ;;
    -api) run_api ;;
    -jax) run_jax ;;
    -jax_f2py_unit_tests) run_jax_f2py_unit_tests ;;
    -all) run_unit; run_dash; run_api; run_jax; run_jax_f2py_unit_tests ;;
esac
