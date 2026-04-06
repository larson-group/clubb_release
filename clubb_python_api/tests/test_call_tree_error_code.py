"""Test wrappers for routines from CLUBB_core/error_code.F90."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clubb_python import clubb_api


def test_initialize_error_headers_callable():
    """initialize_error_headers should be callable without error."""
    clubb_api.initialize_error_headers()
