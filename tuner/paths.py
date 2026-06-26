"""Shared filesystem paths for the CLUBB tuner package."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
RUN_SCRIPTS = REPO_ROOT / "run_scripts"
UTILITIES = REPO_ROOT / "utilities"
INPUT_DIR = REPO_ROOT / "input"
CREATE_CASE_NAMELIST = UTILITIES / "create_case_namelist.py"
CREATE_MULTI_COL_PARAMS = UTILITIES / "create_multi_col_params.py"
DEFAULT_TUNER_STATS = INPUT_DIR / "stats" / "tuning_stats.in"
DEFAULT_TUNABLE_PARAMS = INPUT_DIR / "tunable_parameters" / "tunable_parameters.in"
