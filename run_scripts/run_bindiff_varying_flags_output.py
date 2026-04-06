#!/usr/bin/python3

import argparse
import subprocess
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare two run_clubb_w_varying_flags.py output directories by "
            "diffing each matching flag-set subdirectory."
        )
    )
    parser.add_argument(
        "-v",
        "--verbose",
        type=int,
        default=1,
        help="Verbosity level passed through to run_bindiff_all.py.",
    )
    parser.add_argument(
        "dirs",
        nargs=2,
        help="Two top-level output directories produced by run_clubb_w_varying_flags.py.",
    )
    return parser.parse_args()


def get_flag_dirs(root):
    return {path.name: path for path in root.iterdir() if path.is_dir()}


def get_nc_files(root):
    return {str(path.relative_to(root)) for path in root.rglob("*.nc")}


def print_name_diff(kind, names):
    print(f"{kind}:")
    for name in sorted(names):
        print(f"  {name}")


def main():
    args = parse_args()
    dir1 = Path(args.dirs[0]).resolve()
    dir2 = Path(args.dirs[1]).resolve()

    if not dir1.is_dir() or not dir2.is_dir():
        print("Both inputs must be existing directories.")
        return 2

    if dir1 == dir2:
        print("Input paths resolve to the same directory.")
        return 2

    flag_dirs_1 = get_flag_dirs(dir1)
    flag_dirs_2 = get_flag_dirs(dir2)

    if not flag_dirs_1 or not flag_dirs_2:
        print("Expected both inputs to contain per-flag output subdirectories.")
        return 2

    only_in_1 = set(flag_dirs_1) - set(flag_dirs_2)
    only_in_2 = set(flag_dirs_2) - set(flag_dirs_1)
    common = sorted(set(flag_dirs_1) & set(flag_dirs_2))

    had_difference = False

    if only_in_1:
        print_name_diff(f"Flag directories only in {dir1}", only_in_1)
        had_difference = True

    if only_in_2:
        print_name_diff(f"Flag directories only in {dir2}", only_in_2)
        had_difference = True

    bindiff_script = Path(__file__).resolve().with_name("run_bindiff_all.py")

    for flag_name in common:
        flag_dir_1 = flag_dirs_1[flag_name]
        flag_dir_2 = flag_dirs_2[flag_name]
        nc_files_1 = get_nc_files(flag_dir_1)
        nc_files_2 = get_nc_files(flag_dir_2)

        missing_in_2 = nc_files_1 - nc_files_2
        missing_in_1 = nc_files_2 - nc_files_1

        if missing_in_2:
            print_name_diff(
                f"NetCDF files in {flag_dir_1} but not {flag_dir_2}",
                missing_in_2,
            )
            had_difference = True

        if missing_in_1:
            print_name_diff(
                f"NetCDF files in {flag_dir_2} but not {flag_dir_1}",
                missing_in_1,
            )
            had_difference = True

        print(f"\nComparing flag set: {flag_name}")
        result = subprocess.run(
            [
                sys.executable,
                str(bindiff_script),
                "-v",
                str(args.verbose),
                str(flag_dir_1),
                str(flag_dir_2),
            ],
            check=False,
        )

        if result.returncode != 0:
            had_difference = True

    if had_difference:
        print("\nSUMMARY: differences were detected across varying-flags output.")
        return 1

    print("\nSUMMARY: no differences detected across varying-flags output.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
