"""Input assembly regressions owned by the shared Python unit-test stage.

Check final file formatting and in-memory column preparation without running a
model. Benchmark conversion is stubbed only for the loss assembly check.
"""
import copy

import pytest

from input.scripts import create_case_namelist as builder
from input.scripts.create_multi_col_params import multi_col_settings


def test_final_writer_preserves_sections_indices_and_quotes(tmp_path):
    """Serialize shuffled groups with arrays and stats into the established format."""
    settings = {
        "model_setting": {"title": 'A "quoted" case', "x": [[1, 2], [3, 4]],
                          "_start_index": {"x": [1, 0]}},
        "clubb_stats_nl": {"entry": [' thlm | zt | K | A "quoted" label ']},
        "configurable_clubb_flags_nl": {"enabled": True},
        "clubb_params_nl": {"C8": 0.5},
    }
    before = copy.deepcopy(settings)
    text = builder.write_case_namelist(settings, tmp_path / "case.in").read_text()
    assert text == (
        '&clubb_params_nl\n  C8 = 0.5\n/\n\n&silhs_params_nl\n/\n'
        '&configurable_clubb_flags_nl\n  enabled = .true.\n/\n\n'
        '&configurable_silhs_flags_nl\n/\n'
        '&model_setting\n  title = "A ""quoted"" case"\n'
        '  x(1:2,0) = 1, 2\n  x(1:2,1) = 3, 4\n/\n'
        '&clubb_stats_nl\n  entry(1) = "thlm | zt | K | A ""quoted"" label"\n/\n'
    )
    assert settings == before


@pytest.mark.parametrize("mode,expected", [
    ("duplicate", [2.0, 2.0, 2.0]), ("dup_tweak", [1.0, 2.5, 4.0]),
])
def test_column_modes_preserve_source_and_silhs(mode, expected):
    """Count modes copy SILHS, preserve their input and clip only batch width."""
    source = {"model": {"C11": 2.0, "C8": 0.5}, "silhs": {"coefficient": 0.1}}
    before = copy.deepcopy(source)
    result = multi_col_settings(source, "3", batch_size=5, mode=mode)
    assert result == {"model": {"C11": expected, "C8": [0.5] * 3},
                      "silhs": before["silhs"], "multicol": {"ngrdcol": 3, "batch_size": 3}}
    result["silhs"]["coefficient"] = 99
    assert source == before


@pytest.mark.parametrize("spec", ["C8/NaN:1/2", "C8/0:Infinity/2"])
def test_generated_columns_reject_nonfinite_ranges(spec):
    """Keep numeric validation after generation even though no JSON file is reread."""
    with pytest.raises(ValueError, match="finite"):
        multi_col_settings({"model": {"C8": 0.5}, "silhs": {}}, spec)


def test_loss_batch_preparation_writes_only_final_namelist(tmp_path, monkeypatch):
    """Assemble an exact-duplicate tuning batch without intermediate parameter files."""
    monkeypatch.setattr(builder, "create_normalized_benchmark_file",
                        lambda *args: tmp_path / "benchmark.nc")
    path, fields, _ = builder.create_loss_case_namelist(
        "bomex", tmp_path, ["rtm", "thlm"], config="e3sm",
        case_defaults={"les_stats_file": "unused.nc", "time_average_range": [0, 600],
                       "altitude_comparison_range": [0, 2000]},
        batch_size=3, duplicate_params_for_batch=True, num_time_windows=2,
    )
    from clubb_python_driver.io.namelist import _read_namelist_groups
    groups = _read_namelist_groups(path)
    assert groups["clubb_params_nl"]["c11"] == [0.5] * 3
    assert groups["multicol_def"] == {"ngrdcol": 3, "batch_size": 3}
    assert groups["model_setting"]["time_final"] == 600
    assert fields == ["rtm", "thlm"]
    assert list(tmp_path.iterdir()) == [path]
