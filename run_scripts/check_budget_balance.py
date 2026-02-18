#!/usr/bin/env python3
"""Run CLUBB budget-balance checks across a fixed set of SCM cases."""

import argparse
import os
import subprocess
import sys
from pathlib import Path

RUN_SCRIPTS = Path(__file__).resolve().parent
DEFAULT_CLUBB_ROOT = RUN_SCRIPTS.parent


CASES = [
    "arm",
    "atex",
    "bomex",
    "dycoms2_rf01",
    "dycoms2_rf02_do",
    "dycoms2_rf02_ds",
    "dycoms2_rf02_nd",
    "dycoms2_rf02_so",
    "fire",
    "gabls2",
    "gabls3",
    "mpace_a",
    "mpace_b",
    "rico",
    "wangara",
]


def run_cmd(cmd, *, cwd=None, env=None):
    print(f"+ {' '.join(str(part) for part in cmd)}")
    return subprocess.run(cmd, cwd=cwd, env=env, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run SCM cases and then run postprocessing/check_budgets_balance/checkBudget.py."
    )
    parser.add_argument(
        "clubb_source",
        nargs="?",
        default=None,
        help="Path to CLUBB source root (default: derived from this script location).",
    )
    args = parser.parse_args()

    clubb_source = Path(args.clubb_source).resolve() if args.clubb_source else DEFAULT_CLUBB_ROOT.resolve()
    env = os.environ.copy()
    env["LAPACK"] = "/usr/lib64"

    run_scripts_dir = clubb_source / "run_scripts"
    run_scm_py = run_scripts_dir / "run_scm.py"
    budget_dir = clubb_source / "postprocessing" / "check_budgets_balance"
    if not run_scm_py.exists():
        raise SystemExit(f"Missing run_scm.py at: {run_scm_py}")
    if not budget_dir.exists():
        raise SystemExit(f"Missing budget check directory: {budget_dir}")

    for case in CASES:
        run_cmd(
            [
                sys.executable,
                str(run_scm_py),
                "-dt_main",
                "300",
                "-dt_rad",
                "300",
                "-tout",
                "300",
                # Keep l_lmm_stepping disabled and ensure sampling cadence matches output cadence.
                "-override",
                "l_lmm_stepping=.false.,stats_tsamp=300.0",
                case,
            ],
            cwd=str(run_scripts_dir),
            env=env,
        )

    result = subprocess.run(
        [sys.executable, "checkBudget.py", "all", "0"],
        cwd=str(budget_dir),
        env=env,
        check=False,
    )
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
