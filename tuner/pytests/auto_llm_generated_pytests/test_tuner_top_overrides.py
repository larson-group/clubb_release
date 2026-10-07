"""Top-result reruns must retain candidate parameter values and physics overrides.

src/clubb_loss_driver.F90::clubb_get_loss_for_params replaces initialized
parameter values with candidate vectors. Use the actual parameter-file writer
and namelist override owner to check that both normal and window reruns
reproduce that same effective configuration for each supported override format.
"""

import json
from pathlib import Path
import re

import pytest

from run_scripts import run_tuner_job as cli
from utilities.create_case_namelist import override_value


@pytest.mark.parametrize("jax_backend", [False, True])
@pytest.mark.parametrize("override_format", ["assignments", "json", "json_file", "case_json"])
def test_top_reruns_preserve_tuned_columns_and_other_overrides(
    tmp_path, monkeypatch, jax_backend, override_format,
):
    cases = ["bomex"]
    expected_flags = {"bomex": ".true."}
    settings = {"clubb_params_nl.c8": 0.45, "C11": 0.4, "l_diag_Lscale_from_tau": True}
    if override_format == "assignments":
        override = "C8=0.45,C11=0.4,l_diag_Lscale_from_tau=.true."
    elif override_format == "json":
        override = json.dumps(settings)
    elif override_format == "json_file":
        path = tmp_path / "overrides.json"
        path.write_text(json.dumps(settings))
        override = str(path)
    else:
        cases.append("atex")
        expected_flags["atex"] = ".false."
        override = json.dumps({
            "all": {"clubb_params_nl.c8": 0.45, "C11": 0.4},
            "bomex": {"l_diag_Lscale_from_tau": True},
            "atex": {"l_diag_Lscale_from_tau": False},
        })
    args = cli.parse_args([
        "-cases", *cases, "-fields", "wp2", "-param_ranges", "C8:0.2:0.8",
        "-strategy", "random:2", "-top_n", "2",
        "-output_run_dir", str(tmp_path / "reruns"), "-override", override,
        *(["-jax"] if jax_backend else []),
    ])
    result_path = tmp_path / "results.json"
    result_path.write_text(json.dumps({
        "request": cli.build_request(args),
        "best_results": [
            {"selected_params": {"C8": 0.7}},
            {"selected_params": {"C8": 0.5}},
        ],
    }))
    effective_inputs = []

    def capture_rerun(command, **kwargs):
        params_index = command.index("-params_file") + 1
        text = Path(command[params_index]).read_text()
        text += "\n&configurable_model_flags\nl_diag_Lscale_from_tau=.false.\n/\n"
        overrides = [
            command[index + 1]
            for index, token in enumerate(command) if token == "-override"
        ]
        selected_cases = (
            command[command.index("-cases") + 1].split(",")
            if "-cases" in command else [command[params_index + 1]]
        )
        for case_name in selected_cases:
            effective_inputs.append((case_name, override_value(overrides, text, case_name)))
        return 0

    monkeypatch.setattr(cli, "run_command", capture_rerun)
    assert cli.run_top_results(args, tmp_path / "request.json", result_path, "both") == 0
    assert len(effective_inputs) == 2 * len(cases)
    for case_name, text in effective_inputs:
        match = re.search(r"(?im)^\s*C8\s*=\s*(.+)$", text)
        assert match is not None
        values = tuple(float(value.strip()) for value in match.group(1).rstrip(",").split(","))
        assert values == (0.7, 0.5)
        assert re.search(r"(?im)^\s*C11\s*=\s*0\.4,?\s*$", text)
        logical = re.search(r"(?im)^\s*l_diag_Lscale_from_tau\s*=\s*(\.\w+\.)", text)
        assert logical is not None
        assert logical.group(1) == expected_flags[case_name]
