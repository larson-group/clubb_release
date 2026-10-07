"""Check runner selection and effective timing for stats consistency runs."""

import sys

import pytest

from tests import run_stats_output_consistency as harness


@pytest.mark.parametrize("options, case, jax_option", [
    (["bomex"], "bomex", None),
    (["-jax", "bomex"], "bomex", "-jax"),
    (["rico_silhs", "-jax=cpu"], "rico_silhs", "-jax=cpu"),
    (["-jax=gpu,xla_prealloc", "bomex"], "bomex", "-jax=gpu,xla_prealloc"),
])
def test_all_run_variants_forward_the_runner(monkeypatch, tmp_path, options, case, jax_option):
    monkeypatch.setattr(sys, "argv", [
        "run_stats_output_consistency.py", *options, "-config", "default",
        "-max_iters", "24", "-dt_main", "60", "-coarse_touts", "120,240",
        "-output_root", str(tmp_path),
    ])
    args = harness.parse_args()
    assert args.case_name == case
    window = harness.derive_stats_window(args)
    reference, batches, coarse, window_fine, window_coarse = harness.run_specs(args, window)
    for spec in [reference, *batches, *coarse, window_fine, window_coarse]:
        cmd = harness.build_run_command(args, spec)
        assert cmd[-1] == case
        if jax_option is None:
            assert not any(option.startswith("-jax") for option in cmd)
        else:
            assert cmd.count(jax_option) == 1
        assert cmd[cmd.index("-max_iters") + 1] == "24"
        assert cmd[cmd.index("-dt_main") + 1] == "60"
        assert cmd[cmd.index("-batch_size") + 1] == str(spec.batch_size)
        assert cmd[cmd.index("-config") + 1] == "default"
        if spec.stats_window is not None:
            assert cmd[cmd.index("-stats_tstart") + 1] == str(window[0])
            assert cmd[cmd.index("-stats_tend") + 1] == str(window[1])


def test_automatic_window_uses_shortened_effective_duration(monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "run_stats_output_consistency.py", "-jax", "bomex",
        "-max_iters", "24", "-dt_main", "60", "-coarse_touts", "120,240",
    ])
    args = harness.parse_args()
    monkeypatch.setattr(harness, "read_model_timing", lambda _case: {
        "time_initial": 1000, "time_final": 50000, "dt_main": 300,
    })
    assert harness.derive_stats_window(args) == (1360, 2080)
    args.window_start, args.window_end = 1000, 2800
    with pytest.raises(RuntimeError, match="Invalid stats window"):
        harness.derive_stats_window(args)


def test_payload_mismatch_produces_a_failed_test(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(sys, "argv", [
        "run_stats_output_consistency.py", "bomex", "-jax", "-max_iters", "24",
        "-coarse_touts", "120,240", "-output_root", str(tmp_path),
    ])
    monkeypatch.setattr(harness, "clear_dir", lambda _path: None)
    monkeypatch.setattr(harness, "run_subprocess", lambda _cmd, **_kwargs: None)
    monkeypatch.setattr(harness.os.path, "isfile", lambda _path: True)
    monkeypatch.setattr(harness, "compare_exact_stats", lambda *_args, **_kwargs: (
        ["thlm: payload mismatch"], "different model data",
    ))
    monkeypatch.setattr(harness, "compare_tout_average", lambda *_args, **_kwargs: ([], "match"))
    monkeypatch.setattr(harness, "compare_window_subset", lambda *_args, **_kwargs: ([], "match"))
    with pytest.raises(SystemExit) as exc:
        harness.main()
    assert exc.value.code == 1
    assert "thlm: payload mismatch" in capsys.readouterr().out


@pytest.mark.parametrize("options", [["-max_iters", "0"], ["-jax", "-jax=cpu"]])
def test_invalid_runner_or_duration_is_rejected(monkeypatch, options):
    monkeypatch.setattr(sys, "argv", ["run_stats_output_consistency.py", *options])
    with pytest.raises(SystemExit) as exc:
        harness.parse_args()
    assert exc.value.code == 2


def test_iteration_cap_does_not_extend_the_native_case(monkeypatch):
    """run_scm caps duration; a larger cap must not move the stats window."""
    monkeypatch.setattr(sys, "argv", [
        "run_stats_output_consistency.py", "bomex", "-jax=cpu",
        "-max_iters", "1000", "-coarse_touts", "120,240",
    ])
    args = harness.parse_args()
    monkeypatch.setattr(harness, "read_model_timing", lambda _case: {
        "time_initial": 1000, "time_final": 2440, "dt_main": 60,
    })
    assert harness.derive_stats_window(args) == (1360, 2080)
    args.window_start, args.window_end = 1000, 2800
    with pytest.raises(RuntimeError, match="Invalid stats window"):
        harness.derive_stats_window(args)


@pytest.mark.parametrize('options', [
    ['-jax_tolerance', '3e-9'],
    ['-jax', '-jax_tolerance', '0'],
    ['-jax', '-jax_tolerance', '-1'],
    ['-jax', '-jax_tolerance', 'nan'],
    ['-jax', '-jax_tolerance', 'inf'],
])
def test_tolerance_is_positive_finite_and_jax_only(monkeypatch, options):
    monkeypatch.setattr(sys, 'argv', ['stats', *options])
    with pytest.raises(SystemExit) as exc:
        harness.parse_args()
    assert exc.value.code == 2


def write_stats_fixture(path, values, *, coarse=False):
    import netCDF4
    with netCDF4.Dataset(path, 'w') as ds:
        ds.createDimension('time', 1 if coarse else 2)
        ds.createDimension('col', 1)
        ds.createDimension('bnds', 2)
        ds.createVariable('time', 'f8', ('time',))[:] = [120] if coarse else [60, 120]
        bounds = ds.createVariable('time_bnds', 'f8', ('time', 'bnds'))
        bounds.interval_semantics = '(start, end]'
        bounds[:] = [[0, 120]] if coarse else [[0, 60], [60, 120]]
        ds.createVariable('col', 'i4', ('col',))[:] = [1]
        field = ds.createVariable('thlm', 'f8', ('time', 'col'), fill_value=0)
        field.units = 'K'
        field[:] = [[value] for value in values]


def test_relaxed_batch_and_window_values_keep_native_comparisons_exact(tmp_path):
    reference, test = tmp_path/'reference.nc', tmp_path/'test.nc'
    write_stats_fixture(reference, [0, 1])
    write_stats_fixture(test, [1e-17, 1 + 1e-9])
    assert harness.compare_exact_stats(str(reference), str(test))[0]
    assert harness.compare_window_subset(str(reference), str(test), (0, 120))[0]
    assert not harness.compare_exact_stats(str(reference), str(test), tolerance=3e-9)[0]
    assert not harness.compare_window_subset(str(reference), str(test), (0, 120), tolerance=3e-9)[0]


@pytest.mark.parametrize('fault', ['payload', 'coordinate', 'integer', 'units', 'nan', 'inf', 'missing'])
def test_relaxed_comparisons_still_detect_bad_output(tmp_path, fault):
    import netCDF4
    reference, test = tmp_path/'reference.nc', tmp_path/'test.nc'
    write_stats_fixture(reference, [0, 1])
    write_stats_fixture(test, [1e-17, 1 + 1e-9])
    with netCDF4.Dataset(test, 'a') as ds:
        if fault == 'payload':
            ds['thlm'][1, 0] = 1.001
        elif fault == 'coordinate':
            ds['time'][1] = 120 + 1e-9
        elif fault == 'integer':
            ds['col'][0] = 2
        elif fault == 'units':
            ds['thlm'].units = 'degC'
        elif fault == 'missing':
            ds['thlm'].missing_value = -999
            ds['thlm'][1, 0] = -999
        else:
            ds['thlm'][1, 0] = float(fault)
    assert harness.compare_exact_stats(str(reference), str(test), tolerance=3e-9)[0]
    assert harness.compare_window_subset(str(reference), str(test), (0, 120), tolerance=3e-9)[0]


def test_relaxed_average_accepts_rounding_but_rejects_a_real_error(tmp_path):
    import netCDF4
    fine, coarse = tmp_path/'fine.nc', tmp_path/'coarse.nc'
    write_stats_fixture(fine, [0, 2])
    write_stats_fixture(coarse, [1 + 1e-9], coarse=True)
    assert harness.compare_tout_average(str(fine), str(coarse))[0]
    assert not harness.compare_tout_average(str(fine), str(coarse), tolerance=3e-9)[0]
    with netCDF4.Dataset(coarse, 'a') as ds:
        ds['thlm'][0, 0] = 1.001
    assert harness.compare_tout_average(str(fine), str(coarse), tolerance=3e-9)[0]
