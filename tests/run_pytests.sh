#!/usr/bin/env bash
# Run focused pytest suites in their owned environments.
# Usage: tests/run_pytests.sh [-unit|-dash|-api|-jax|-all] [-include_generated] [pytest options]
# Generated tests remain provisional even when explicitly included in a run.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODE=-unit
case "${1:-}" in
    -unit|-dash|-api|-jax|-all) MODE="$1"; shift ;;
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
        collection_args=(--ignore-glob='*/auto_llm_generated_pytests')
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
    run_pytest "$python_bin" utilities/pytests tuner/pytests tests/pytests
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

run_jax() {
    for arg in ${PYTEST_ARGS[@]+"${PYTEST_ARGS[@]}"}; do
        case "$arg" in
            --junitxml|--junitxml=*|--junit-xml|--junit-xml=*)
                echo "JAX writes one XML report per module; select its directory with CLUBB_PYTEST_OUTPUT_DIR." >&2
                return 2
                ;;
        esac
    done
    if ! have_tests clubb_jax/pytests; then
        echo "No admitted JAX pytests; no tests executed. Use -include_generated for provisional coverage."
        return
    fi
    python3 clubb_jax/run_jax.py -profile=cpu -init_env
    local venv="${CLUBB_JAX_VENV:-$REPO_ROOT/.venv-jax}"
    [[ "$venv" == /* ]] || venv="$REPO_ROOT/$venv"
    local python_bin="$venv/bin/python"
    local reports="${CLUBB_PYTEST_OUTPUT_DIR:-$REPO_ROOT/output/tests/pytests/jax}"
    [[ "$reports" == /* ]] || reports="$REPO_ROOT/$reports"
    local path relative report result failures=0 selected=0 executed=0
    # Isolate module-level flags and JIT caches between independently selected modules.
    trap 'exit 130' INT
    trap 'exit 143' TERM
    while IFS= read -r path; do
        if [[ "${CLUBB_PYTEST_INCLUDE_GENERATED:-0}" != 1 && "$path" == */auto_llm_generated_pytests/* ]]; then
            continue
        fi
        relative="${path#$REPO_ROOT/clubb_jax/pytests/}"
        report="$reports/${relative%.py}.xml"
        mkdir -p "$(dirname "$report")"
        rm -f "$report"
        selected=$((selected + 1))
        echo "JAX pytest module: $relative"
        if "$python_bin" -m pytest "$path" \
            --capture=sys --junitxml="$report" \
            ${PYTEST_ARGS[@]+"${PYTEST_ARGS[@]}"}; then
            executed=$((executed + 1))
        else
            result=$?
            case "$result" in
                130|143) return "$result" ;;
                5) ;;  # A forwarded -k/-m filter can deselect every check in this module.
                *)
                    echo "JAX pytest module failed: $relative (exit $result)" >&2
                    failures=$((failures + 1))
                    # Record a process error if pytest exits before finalizing its report.
                    if [[ ! -s "$report" ]]; then
                        "$python_bin" - "$report" "$relative" "$result" <<'PY_REPORT'
from pathlib import Path
import sys
import xml.etree.ElementTree as ET
report, module, result = sys.argv[1:]
suite = ET.Element("testsuite", name=module, tests="1", errors="1", failures="0", skipped="0")
case = ET.SubElement(suite, "testcase", classname=module, name="pytest_module_process")
message = f"Module process exited {result} before writing its pytest report; see the console."
ET.SubElement(case, "error", message=message).text = message
ET.ElementTree(suite).write(Path(report), encoding="utf-8", xml_declaration=True)
PY_REPORT
                    fi
                    ;;
            esac
        fi
    done < <(find "$REPO_ROOT/clubb_jax/pytests" -type f \( -name 'test_*.py' -o -name '*_test.py' \) | sort)
    echo "JAX pytest commands: $selected modules selected, $failures module commands failed."
    if [[ "$failures" -gt 0 ]]; then
        return 1
    fi
    if [[ "$executed" -eq 0 ]]; then
        echo "No JAX checks matched the pytest selection; no tests executed." >&2
        return 5
    fi
}

case "$MODE" in
    -unit) run_unit ;;
    -dash) run_dash ;;
    -api) run_api ;;
    -jax) run_jax ;;
    -all) run_unit; run_dash; run_api; run_jax ;;
esac
