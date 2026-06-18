"""CLUBB loss-driver and tuning utilities."""

from pathlib import Path
import sys

_REPO_ROOT = Path(__file__).resolve().parent.parent
_API_ROOT = _REPO_ROOT / "clubb_python_api"

if str(_API_ROOT) not in sys.path:
    sys.path.insert(0, str(_API_ROOT))
