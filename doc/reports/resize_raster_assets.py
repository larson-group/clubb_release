#!/usr/bin/env python3
"""Normalize published report PNG/JPEG assets to the shared 1200 px cap."""

from __future__ import annotations

from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from dash_app.reports_tab.publisher import MAX_RASTER_DIMENSION, cap_raster_image


def main() -> None:
    report_root = Path(__file__).resolve().parent
    changed = 0
    for path in sorted(report_root.glob("*/figures/*")):
        if path.suffix.lower() not in {".jpg", ".jpeg", ".png"}:
            continue
        if cap_raster_image(path, path, max_dimension=MAX_RASTER_DIMENSION):
            changed += 1
            print(f"resized {path.relative_to(REPO_ROOT)}")
    print(f"normalized {changed} report raster asset{'s' if changed != 1 else ''}")


if __name__ == "__main__":
    main()
