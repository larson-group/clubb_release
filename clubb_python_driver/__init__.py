"""Python CLUBB standalone driver package.

This package lives at the repository root, while the Python API package lives
under ``clubb_python_api/``. Add that sibling directory as a fallback so
``clubb_python`` imports work when the driver is launched directly. An installed
runtime placed on PYTHONPATH by run_scm.py must take precedence, including its
matching F2PY extension.
"""

from pathlib import Path
import sys

_REPO_ROOT = Path(__file__).resolve().parent.parent
_API_ROOT = _REPO_ROOT / "clubb_python_api"

if str(_API_ROOT) not in sys.path:
    sys.path.append(str(_API_ROOT))
