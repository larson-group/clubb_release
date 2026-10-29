from pathlib import Path

import pytest

from dash_app.run_tab.settings import (
    write_temp_settings,
)
from dash_app.shared.runtime import read_file_tail, read_latest_run_progress


REPO_ROOT = Path(__file__).resolve().parents[4]
CONFIG_ROOT = REPO_ROOT / "input"


def test_config_update_places_model_flag_inside_clubb_flags():
    import json
    from input.scripts.flags_json_to_namelist import flag_settings
    base = flag_settings()
    path = write_temp_settings(base, {"Lscale_method": "2"}, "clubb_flags_test_")
    try:
        values = json.loads(Path(path).read_text())
    finally:
        Path(path).unlink(missing_ok=True)
    assert values["model"]["Lscale_method"] == 2
    assert "Lscale_method" not in values["silhs"]


def test_namelist_updater_rejects_unknown_settings():
    with pytest.raises(ValueError, match="not present"):
        write_temp_settings({"model": {}, "silhs": {}}, {"unknown_setting": "1"}, "clubb_flags_test_")


from input.scripts.create_case_namelist import available_configs


@pytest.mark.parametrize("config", available_configs())
def test_all_configs_accept_common_tunable_parameter(config):
    """Every selectable JSON config supports a typed parameter edit."""
    import json
    from input.scripts.parameters_json_to_namelist import parameter_settings
    base = parameter_settings(config=config)
    path = write_temp_settings(base, {"C8": "0.4"}, "clubb_params_test_")
    try:
        assert json.loads(Path(path).read_text())["model"]["C8"] == 0.4
    finally:
        Path(path).unlink(missing_ok=True)


def test_latest_progress_reads_only_the_last_iteration_in_a_log_tail(tmp_path):
    path = tmp_path / "run.log"
    path.write_text(
        "old\niteration: 1 / 40 -- time = 60s / 2400s\n"
        + "x" * 9000
        + "\niteration: 17 / 40 -- time = 1020s / 2400s\n",
        encoding="utf-8",
    )

    assert read_latest_run_progress(path) == (17, 40)


def test_file_tail_returns_only_the_requested_lines_and_eof_cursor(tmp_path):
    path = tmp_path / "run.log"
    path.write_text("".join(f"line {index}\n" for index in range(6000)), encoding="utf-8")

    chunk, cursor, eof = read_file_tail(path, 5000)

    assert chunk.decode().splitlines() == [f"line {index}" for index in range(1000, 6000)]
    assert cursor == path.stat().st_size
    assert eof is True
