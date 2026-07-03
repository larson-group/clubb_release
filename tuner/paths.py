"""Shared filesystem paths for the CLUBB tuner package."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
RUN_SCRIPTS = REPO_ROOT / "run_scripts"
UTILITIES = REPO_ROOT / "utilities"
INPUT_DIR = REPO_ROOT / "input"
CREATE_CASE_NAMELIST = UTILITIES / "create_case_namelist.py"
CREATE_MULTI_COL_PARAMS = UTILITIES / "create_multi_col_params.py"
DEFAULT_TUNER_STATS = INPUT_DIR / "stats" / "all_tuning_stats.in"
DEFAULT_TUNABLE_CONFIG = INPUT_DIR / "parameter_and_flag_configs" / "default"
DEFAULT_TUNABLE_PARAMS = DEFAULT_TUNABLE_CONFIG / "tunable_parameters.in"
