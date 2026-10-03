import pytest

from utilities.create_case_namelist import override_value


def test_bare_override_replaces_existing_value():
    result = override_value("dt_main=2.0", "&model_setting\n  dt_main = 1.0,\n/\n")

    assert "dt_main = 2.0" in result


def test_grouped_override_adds_missing_value():
    result = override_value(
        "microphysics_setting.lh_num_samples=8",
        "&microphysics_setting\n  microphys_scheme = 'none',\n/\n",
    )

    assert "lh_num_samples = 8," in result


def test_grouped_override_only_updates_its_namelist():
    result = override_value(
        "second_setting.value=3",
        "&first_setting\n  value = 1,\n/\n&second_setting\n  value = 2,\n/\n",
    )

    assert "&first_setting\n  value = 1," in result
    assert "&second_setting\n  value = 3," in result


def test_missing_bare_override_explains_how_to_add_a_key():
    with pytest.raises(ValueError, match="could not be matched for find-and-replace") as exc_info:
        override_value("lh_num_samples=8", "&microphysics_setting\n/\n")

    assert "[namelist.]key=value" in str(exc_info.value)


@pytest.mark.parametrize("case_name", ["mc3e", "bomex"])
@pytest.mark.parametrize("as_file", [False, True])
def test_flat_json_applies_to_every_case(tmp_path, case_name, as_file):
    import json

    data = json.dumps({"dt_main": 2.0, "model_setting.enabled": False})
    if as_file:
        path = tmp_path / "shared.json"
        path.write_text(data)
        data = str(path)
    result = override_value(data, "&model_setting\n dt_main = 1.0,\n/\n", case_name)
    assert "dt_main = 2.0" in result
    assert "enabled = .false." in result


def test_case_json_selects_only_current_case_and_preserves_override_order():
    profile = '{"mc3e": {"dt_main": 3.0}}'
    namelist = "&model_setting\n dt_main = 1.0,\n/\n"
    assert "dt_main = 3.0" in override_value(["dt_main=2.0", profile], namelist, "mc3e")
    assert "dt_main = 2.0" in override_value(["dt_main=2.0", profile], namelist, "bomex")
    assert "dt_main = 4.0" in override_value([profile, "dt_main=4.0"], namelist, "mc3e")


@pytest.mark.parametrize("data", [[], {"typo": {}}, {"mc3e": []}, {"mc3e": {"dt_main": None}},
                                   {"mc3e": {}, "dt_main": 1}, {"dt_main": None}, {"dt_main": [1]}])
def test_invalid_json_override_is_rejected(tmp_path, data):
    import json

    path = tmp_path / "invalid.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        override_value(str(path), "&model_setting\n dt_main = 1.0,\n/\n", "mc3e")


def test_missing_json_override_reports_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        override_value(str(tmp_path / "missing.json"), "&model_setting\n/\n", "bomex")


def test_assignment_column_lists_still_work():
    result = override_value("C8=0.8,0.7,C11=1.0,1.1", "&params\n C8 = 0.0,\n C11 = 0.0,\n/\n")
    assert "C8 = 0.8, 0.7" in result
    assert "C11 = 1.0, 1.1" in result


@pytest.mark.parametrize("profile_name", ["core", "host"])
@pytest.mark.parametrize("case_name", ["cgils_s6", "cgils_s11", "cgils_s12", "mc3e", "bomex"])
def test_grid_json_generates_the_same_namelist_as_legacy_assignments(tmp_path, profile_name, case_name):
    import json
    from pathlib import Path
    from utilities.create_case_namelist import create_case_namelist_file

    repo = Path(__file__).resolve().parents[2]
    import re
    pipeline = (repo / "jenkins_tests/clubb_generalized_grid/Jenkinsfile").read_text()
    profiles = re.findall(r"-override '(\{.*?\})'", pipeline, re.S)
    profile = profiles[0 if profile_name == "core" else 1]
    settings = {}
    if case_name in {"cgils_s6", "cgils_s11", "cgils_s12"}:
        settings.update({"thlm_sponge_damp_settings%l_sponge_damping": False,
                         "rtm_sponge_damp_settings%l_sponge_damping": False})
    if profile_name == "core" and case_name == "mc3e":
        settings["time_final"] = 1944000.0
    expected_overrides = ",".join(f"{key}={'.false.' if value is False else value}" for key, value in settings.items())
    source = repo / "input/case_setups" / f"{case_name}_model.in"
    original = source.read_bytes()
    output = tmp_path / "output"
    expected = create_case_namelist_file(case_name, output, override=expected_overrides).read_text()
    actual = create_case_namelist_file(case_name, output, override=profile).read_text()
    assert actual == expected
    assert source.read_bytes() == original


def test_run_scm_forwards_mixed_overrides_to_the_real_generator(tmp_path):
    import json
    from pathlib import Path
    from types import SimpleNamespace
    from run_scripts import run_scm

    profile = tmp_path / "cases.json"
    profile.write_text(json.dumps({"bomex": {"dt_main": 3.0}}))
    fields = ("config params flags silhs_params stats multicol batch_size zt_grid zm_grid "
              "nzmax debug max_iters dt_main dt_rad tout stats_tstart stats_tend").split()
    args = SimpleNamespace(**dict.fromkeys(fields), case_name="bomex",
                           override=["dt_main=2.0", str(profile)])
    result = Path(run_scm.create_case_namelist(args, str(tmp_path)))
    assert "dt_main = 3.0" in result.read_text()


def test_varying_flags_retains_flag_set_with_a_forwarded_json_profile(tmp_path, monkeypatch):
    import json
    from types import SimpleNamespace
    from pathlib import Path
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[2] / "run_scripts"))
    from run_scripts import run_clubb_w_varying_flags as runner
    from utilities.create_case_namelist import create_case_namelist_file

    profile = tmp_path / "cases.json"
    profile.write_text(json.dumps({"mc3e": {"time_final": 1944000.0}}))
    args = SimpleNamespace(tout=None, max_iters=None, run_scm_extra_args=["-override", str(profile)])
    tasks = runner.build_tasks(str(tmp_path), {"flag1": {"penta_solve_method": 1}}, ["mc3e", "bomex"], args)
    for task in tasks:
        cmd = task["cmd"]
        overrides = [cmd[i + 1] for i, arg in enumerate(cmd) if arg == "-override"]
        assert overrides == ["penta_solve_method=1", str(profile)]
        text = create_case_namelist_file(task["case"], tmp_path / task["case"], override=overrides).read_text()
        assert "penta_solve_method = 1" in text
        if task["case"] == "mc3e":
            assert "time_final = 1944000.0" in text


def test_batch_runner_forwards_the_same_json_argument_to_every_case(tmp_path, monkeypatch):
    import sys
    from concurrent.futures import Future
    from run_scripts import run_scm_all as runner

    calls = []
    profile = str(tmp_path / "handled_by_generator.json")
    class Pool:
        def __init__(self, max_workers):
            assert max_workers == 2
        def __enter__(self):
            return self
        def __exit__(self, *_args):
            pass
        def submit(self, _fn, case, options, verbose):
            calls.append((case, list(options)))
            future = Future()
            future.set_result((case, 0, ""))
            return future
    monkeypatch.setattr(runner, "ProcessPoolExecutor", Pool)
    monkeypatch.setattr(sys, "argv", ["run_scm_all.py", "-cases", "mc3e,bomex", "-override", profile])
    assert runner.main() == 0
    assert calls == [("mc3e", ["-override", profile]), ("bomex", ["-override", profile])]


@pytest.mark.parametrize("case_args", [[], ["mc3e"]])
def test_flag_runner_does_not_mistake_json_filename_for_case(case_args, monkeypatch):
    import sys
    from pathlib import Path
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[2] / "run_scripts"))
    from run_scripts import run_clubb_w_varying_flags as runner

    monkeypatch.setattr(sys, "argv", ["run_clubb_w_varying_flags.py", "-override", "cases.json",
                                     "-override", "penta_solve_method=1", *case_args])
    args = runner.get_cli_args()
    assert args.case_name == ("mc3e" if case_args else None)
    assert args.run_scm_extra_args == ["-override", "cases.json", "-override", "penta_solve_method=1"]


@pytest.mark.parametrize("as_file", [False, True])
def test_all_defaults_case_settings_and_later_arguments_have_ordered_precedence(tmp_path, as_file):
    import json
    data = json.dumps({"all": {"dt_main": 2.0, "model_setting.enabled": False},
                       "mc3e": {"dt_main": 3.0, "model_setting.enabled": True}})
    if as_file:
        path = tmp_path / "defaults.json"
        path.write_text(data)
        data = str(path)
    namelist = "&model_setting\n dt_main = 1.0,\n/\n"
    selected = override_value(data, namelist, "mc3e")
    assert "dt_main = 3.0" in selected and "enabled = .true." in selected
    unlisted = override_value(data, namelist, "bomex")
    assert "dt_main = 2.0" in unlisted and "enabled = .false." in unlisted
    assert "dt_main = 4.0" in override_value([data, "dt_main=4.0"], namelist, "mc3e")


def test_all_only_json_works_without_a_case_name():
    assert "dt_main = 2.0" in override_value('{"all":{"dt_main":2.0}}',
                                            "&model_setting\n dt_main = 1.0,\n/\n")


def test_all_section_does_not_hide_unknown_cases():
    with pytest.raises(ValueError, match="Unknown case"):
        override_value('{"all":{"dt_main":2.0},"typo":{"dt_main":3.0}}',
                       "&model_setting\n dt_main = 1.0,\n/\n", "bomex")
