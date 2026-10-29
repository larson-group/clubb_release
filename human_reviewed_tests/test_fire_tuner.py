#!/usr/bin/env python3
"""Configure and run the FIRE tuning workflow."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Configure and run FIRE tuning.", add_help=False, allow_abbrev=False)
    parser.add_argument("-h", "-help", action="help", help="Show this help and exit.")
    parser.add_argument(
        "clubb_source",
        nargs="?",
        default=".",
        help="Path to CLUBB source root (default: current directory).",
    )
    args = parser.parse_args()

    clubb_source = Path(args.clubb_source).resolve()
    run_tuner_path = clubb_source / "run_scripts" / "run_tuner.py"
    if not run_tuner_path.is_file():
        print(f"{run_tuner_path} does not exist")
        return 1

    result = subprocess.run(
        [
            sys.executable,
            str(run_tuner_path),
            "-run_type",
            "single",
            "-run_case",
            "fire",
            "-stats_tune",
            "var:cloud_frac,var:rcm",
            "-override",
            'debug_level=0,microphys_scheme="khairoutdinov_kogan",dt_main=60.,dt_rad=60.,stats_tsamp=60.,stats_tout=60.',
        ],
        cwd=str(clubb_source),
        check=False,
    )

    if result.returncode == 0:
        print("Fire tuner successful!")
    else:
        print("Fire tuner failed!")

    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
