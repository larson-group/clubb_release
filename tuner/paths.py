"""Shared filesystem paths for the CLUBB tuner package."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
RUN_SCRIPTS = REPO_ROOT / "run_scripts"
UTILITIES = REPO_ROOT / "utilities"
INPUT_DIR = REPO_ROOT / "input"
CREATE_CASE_NAMELIST = INPUT_DIR / "scripts" / "create_case_namelist.py"
CREATE_MULTI_COL_PARAMS = INPUT_DIR / "scripts" / "create_multi_col_params.py"
DEFAULT_TUNER_STATS = "all_tuning_stats"
DEFAULT_TUNABLE_CONFIG = "default"
DEFAULT_TUNABLE_PARAMS = INPUT_DIR / "clubb_params.jsonc"
