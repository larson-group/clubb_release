from pathlib import Path

import pytest

from dash_app.run_tab.namelist import (
    apply_updates_to_lines,
    read_namelist_entries,
    write_temp_namelist,
)
from dash_app.shared.runtime import read_file_tail, read_latest_run_progress


REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG_ROOT = REPO_ROOT / "input/parameter_and_flag_configs"


def test_config_update_places_model_flag_inside_clubb_flags():
    base = REPO_ROOT / "input/parameter_and_flag_configs/default/configurable_model_flags.in"
    path = write_temp_namelist(
        base,
        {"l_diag_Lscale_from_tau": ".true."},
        "clubb_flags_test_",
    )
    try:
        values = {entry["name"]: entry["value"] for entry in read_namelist_entries(path)}
        rendered = Path(path).read_text(encoding="utf-8")
    finally:
        Path(path).unlink(missing_ok=True)

    clubb_flags = rendered.split("&configurable_clubb_flags_nl", 1)[1].split("/", 1)[0]
    assert values["l_diag_Lscale_from_tau"] == ".true."
    assert "l_diag_Lscale_from_tau = .true." in clubb_flags


def test_namelist_updater_rejects_unknown_settings():
    with pytest.raises(ValueError, match="not present"):
        apply_updates_to_lines(["&example\n", "/\n"], {"unknown_setting": "1"})


@pytest.mark.parametrize(
    "tunable_path",
    sorted(CONFIG_ROOT.glob("*/tunable_parameters.in")),
    ids=lambda path: path.parent.name,
)
def test_all_configs_accept_common_tunable_parameter(tunable_path):
    """Every selectable config supports a common tunable parameter."""
    path = write_temp_namelist(
        tunable_path,
        {"C8": "0.4"},
        "tunable_parameters_test_",
    )
    try:
        values = {entry["name"]: entry["value"] for entry in read_namelist_entries(path)}
    finally:
        Path(path).unlink(missing_ok=True)

    assert values["C8"] == "0.4"


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
