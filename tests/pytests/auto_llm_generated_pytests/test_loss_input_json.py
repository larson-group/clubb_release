"""Provisional checks for JSON parameter preparation in the shared loss workflow.

The Python Jenkins unit stage owns these bounded parser/generator checks. No
SCM case, native executable or benchmark data is used.
"""
import json
from pathlib import Path
import sys

import pytest

from tests import run_loss_output_consistency as loss
from tests import check_mirrored_multi_col_output as mirror


@pytest.mark.parametrize("selection", [None, "e3sm", "custom"])
def test_loss_workflow_resolves_presets_and_files_before_generating_columns(tmp_path, monkeypatch, selection):
    """Use the real public parser and column generator with shared-config inheritance."""
    params = selection
    if selection == "custom":
        path = tmp_path / "custom.json"
        path.write_text('{"model":{"C11":0.75}}')
        params = str(path)
    command = ["run_loss_output_consistency.py", "bomex", "-config", "e3sm",
               "-flags", "e3sm", "-output_root", str(tmp_path)]
    if params is not None:
        command.extend(["-params", params])
    monkeypatch.setattr(sys, "argv", command)
    args = loss.parse_args()
    assert args.flags == "e3sm"
    assert args.params == params
    output = tmp_path / "columns"
    output.mkdir()
    result = Path(loss.build_shared_multicol_params(args, str(output)))
    values = json.loads(result.read_text())
    assert values["multicol"] == {"ngrdcol": 4, "batch_size": 4}
    assert values["model"]["C11"] == [0.75 if selection == "custom" else 0.5] * 4
    assert len(values["model"]["C8"]) == 4
    assert (output / "base_parameters.json").is_file()


def test_mirror_wrapper_rejects_forwarded_parameter_alias():
    """The new spelling must not override the mirror test's generated column order."""
    with pytest.raises(SystemExit):
        mirror.parse_args(["-case", "bomex", "--", "-params", "custom.json"])
