#!/usr/bin/env bash
# Run focused pytest suites in their owned environments.
# Usage: tests/run_pytests.sh [-unit|-dash|-api|-jax|-jax_f2py_unit_tests|-all] [-include_generated] [pytest options]
# Generated tests remain provisional even when explicitly included in a run.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODE=-unit
case "${1:-}" in
    -unit|-dash|-api|-jax|-jax_f2py_unit_tests|-all) MODE="$1"; shift ;;
    -h|-help) sed -n '2,4p' "$0"; exit 0 ;;
esac
PYTEST_ARGS=()
for arg in "$@"; do
    if [[ "$arg" == -include_generated ]]; then
        export CLUBB_PYTEST_INCLUDE_GENERATED=1
    else
        PYTEST_ARGS+=("$arg")
    fi
done
cd "$REPO_ROOT"

have_tests() {
    local path target
    for target in "$@"; do
        while IFS= read -r path; do
            if [[ "${CLUBB_PYTEST_INCLUDE_GENERATED:-0}" == 1 || "$path" != */auto_llm_generated_pytests/* ]]; then
                return 0
            fi
        done < <(find "$target" -type f \( -name 'test_*.py' -o -name '*_test.py' \))
    done
    return 1
}

run_pytest() {
    local python_bin="$1"
    shift
    if ! have_tests "$@"; then
        echo "No admitted pytests in $*; no tests executed. Use -include_generated for provisional coverage."
        return
    fi
    local collection_args=()
    # Admission belongs to the wrapper; pytest itself has no repository-wide policy.
    if [[ "${CLUBB_PYTEST_INCLUDE_GENERATED:-0}" != 1 ]]; then
        collection_args=(--ignore-glob='*/auto_llm_generated_pytests/*')
    fi
    "$python_bin" -m pytest ${collection_args[@]+"${collection_args[@]}"} \
        "$@" ${PYTEST_ARGS[@]+"${PYTEST_ARGS[@]}"}
}

run_unit() {
    local python_bin="${CLUBB_PYTHON:-}"
    if [[ -z "$python_bin" ]]; then
        python3 utilities/setup_python_venv.py python
        local venv="${CLUBB_PYTHON_VENV:-$REPO_ROOT/.venv-python}"
        [[ "$venv" == /* ]] || venv="$REPO_ROOT/$venv"
        python_bin="$venv/bin/python"
    fi
    run_pytest "$python_bin" input/pytests utilities/pytests tuner/pytests tests/pytests
}

run_dash() {
    local python_bin="${CLUBB_PYTHON:-}"
    if [[ -z "$python_bin" ]]; then
        python3 utilities/setup_python_venv.py dash
        local venv="${CLUBB_PYTHON_VENV:-$REPO_ROOT/.venv-python}"
        [[ "$venv" == /* ]] || venv="$REPO_ROOT/$venv"
        python_bin="$venv/bin/python"
    fi
    run_pytest "$python_bin" dash_app/pytests
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

# JAX checks that need no compiled Fortran. The f2py_unit_tests folders need the
# Python API build and run with -jax_f2py_unit_tests instead.
run_jax() {
    prepare_jax
    local PYTEST_ARGS=(--ignore-glob='*/f2py_unit_tests' ${PYTEST_ARGS[@]+"${PYTEST_ARGS[@]}"})
    run_pytest "$JAX_PYTHON" clubb_jax/pytests
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
    # Reviewed tests move from the generated folder to pytests/f2py_unit_tests.
    local dir dirs=()
    for dir in clubb_jax/pytests/f2py_unit_tests clubb_jax/pytests/auto_llm_generated_pytests/f2py_unit_tests; do
        if [[ -d "$dir" ]]; then
            dirs+=("$dir")
        fi
    done
    PYTHONPATH="$f2py_dir:$REPO_ROOT${PYTHONPATH:+:$PYTHONPATH}" run_pytest "$JAX_PYTHON" "${dirs[@]}"
}

case "$MODE" in
    -unit) run_unit ;;
    -dash) run_dash ;;
    -api) run_api ;;
    -jax) run_jax ;;
    -jax_f2py_unit_tests) run_jax_f2py_unit_tests ;;
    -all) run_unit; run_dash; run_api; run_jax; run_jax_f2py_unit_tests ;;
esac
