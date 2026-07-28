#!/usr/bin/env python3
"""Run fresh PDF-9 SCM cases, then rebuild the Misc-lab cache.

This is the reproducible command behind the Misc ``Generate / refresh``
button.  It intentionally leaves the immutable raw SAM snapshots in place and
regenerates only CLUBB's parcel ledger and the compact dashboard cache.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
SAM_ROOT = Path(
    "/home/pub/les_and_clubb_benchmark_runs/sam_benchmark_runs/"
    "JULY_2017_3D_RECREATIONS"
)
RUNS_DIR = REPO_ROOT / "output_pdf_bakeoff/pdf9_parcel_runs"
STATS_FILE = REPO_ROOT / "input/stats/les_advance_stats.in"
CACHE_DIR = REPO_ROOT / "output_pdf_bakeoff/pdf9_parcel_ensemble_lab"


def _run(command: list[str]) -> None:
    print("+ " + " ".join(command), flush=True)
    subprocess.run(command, cwd=REPO_ROOT, check=True)


def refresh(*, run_cases: bool = True) -> None:
    """Create fresh parcel-ledger statistics and the dashboard cache."""

    if run_cases:
        for case_name in ("arm", "bomex"):
            output_dir = RUNS_DIR / case_name
            _run(
                [
                    sys.executable,
                    "-u",
                    "run_scripts/run_scm.py",
                    "-out_dir",
                    str(output_dir),
                    "-stats",
                    str(STATS_FILE),
                    "-override",
                    "iiPDF_type=9",
                    case_name,
                ]
            )

    cases = []
    for case_name in ("arm", "bomex"):
        stats_path = RUNS_DIR / case_name / f"{case_name}_stats.nc"
        sam_run = SAM_ROOT / f"{case_name}_3d"
        cases.extend(("--case", f"{case_name}::{stats_path}::{sam_run}"))
    _run(
        [
            sys.executable,
            "-u",
            "utilities/generate_pdf9_parcel_ensemble_diagnostics.py",
            "--output-dir",
            str(CACHE_DIR),
            *cases,
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache-only",
        action="store_true",
        help="Rebuild the dashboard cache from existing PDF-9 statistics.",
    )
    args = parser.parse_args()
    refresh(run_cases=not args.cache_only)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
