#!/usr/bin/env python3
"""Configure and run the FIRE tuning workflow."""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


FIRE_STATS_ENTRIES = [
    'entry(1) = "cloud_frac | zt | - | Cloud fraction (between 0 and 1)"',
    'entry(2) = "rcm | zt | kg/kg | Cloud water mixing ratio"',
]


def replace_or_fail(path: Path, pattern: str, replacement: str, *, count: int = 0) -> None:
    text = path.read_text(encoding="utf-8")
    new_text, num_subs = re.subn(pattern, replacement, text, count=count, flags=re.MULTILINE)
    if num_subs == 0:
        raise ValueError(f"Could not find pattern in {path}: {pattern}")
    path.write_text(new_text, encoding="utf-8")


def write_fire_stats(path: Path) -> None:
    lines = ["&clubb_stats_nl"]
    lines.extend(f"  {entry}" for entry in FIRE_STATS_ENTRIES)
    lines.extend(["/", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Configure and run FIRE tuning.")
    parser.add_argument(
        "clubb_source",
        nargs="?",
        default=".",
        help="Path to CLUBB source root (default: current directory).",
    )
    args = parser.parse_args()

    clubb_source = Path(args.clubb_source).resolve()
    run_tuner_path = clubb_source / "run_scripts" / "run_tuner.py"
    fire_model_path = clubb_source / "input" / "case_setups" / "fire_model.in"
    fire_stats_path = clubb_source / "input" / "stats" / "fire_stats.in"

    try:
        replace_or_fail(
            fire_model_path,
            r"^\s*debug_level\s*=\s*\d+\s*$",
            "debug_level   = 0",
            count=1,
        )
        replace_or_fail(
            fire_model_path,
            r'^\s*microphys_scheme\s*=.*$',
            'microphys_scheme = "khairoutdinov_kogan"',
            count=1,
        )
        replace_or_fail(
            fire_model_path,
            r"^\s*dt_main\s*=.*$",
            "dt_main    = 60.0",
            count=1,
        )
        replace_or_fail(
            fire_model_path,
            r"^\s*dt_rad\s*=.*$",
            "dt_rad = 60.0",
            count=1,
        )
        replace_or_fail(
            fire_model_path,
            r"^\s*stats_tsamp\s*=.*$",
            "stats_tsamp  = 60",
            count=1,
        )
        replace_or_fail(
            fire_model_path,
            r"^\s*stats_tout\s*=.*$",
            "stats_tout   = 60",
            count=1,
        )
        write_fire_stats(fire_stats_path)
    except (OSError, ValueError) as exc:
        print(str(exc))
        return 1

    if not run_tuner_path.is_file():
        print(f"{run_tuner_path} does not exist")
        return 1

    result = subprocess.run(
        [
            sys.executable,
            str(run_tuner_path),
            "--run-type",
            "single",
            "--run-case",
            "fire",
            "--stats-tune",
            str(fire_stats_path),
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
