"""Python CLUBB standalone driver package.

This package lives at the repository root, while the Python API package lives
under ``clubb_python_api/``.  Add that sibling directory to ``sys.path`` so
``clubb_python`` imports work when the driver is launched via
``python -m clubb_python_driver.clubb_standalone`` from the repo root.
"""

from pathlib import Path
import sys

_REPO_ROOT = Path(__file__).resolve().parent.parent
_API_ROOT = _REPO_ROOT / "clubb_python_api"

if str(_API_ROOT) not in sys.path:
    sys.path.insert(0, str(_API_ROOT))
