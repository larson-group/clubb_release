"""Column reversal and failure detection for the standalone mirror test."""
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
from types import SimpleNamespace

import netCDF4
import numpy as np
import pytest

from tests import check_mirrored_multi_col_output as mirror


def read_params(path):
    text = Path(path).read_text().split('&clubb_params_nl\n')[1].split('/')[0]
    return {name.strip(): np.array([float(v) for v in values.split(',')])
            for name, values in (line.split('=', 1) for line in text.splitlines() if '=' in line)}


def write_stats(path, values, name='rtm'):
    values = np.ma.array(values, dtype=float)
    with netCDF4.Dataset(path, 'w') as ds:
        ds.createDimension('time', 1)
        ds.createDimension('col', len(values))
        var = ds.createVariable(name, 'f8', ('time', 'col'), fill_value=-999.)
        var[0, :] = values


def test_reverse_generated_hypergrid_exactly(tmp_path):
    base = tmp_path / 'base.in'
    base.write_text('&clubb_params_nl\nC8=0.3\nC11=1.0\nC1=2.0\n/\n')
    forward, reverse = mirror.generate_parameter_files(tmp_path, 'C8/0.2:0.8/4,C11/1:2/3', base)
    a, b = read_params(forward), read_params(reverse)
    assert len(a['C8']) == 12
    assert 'ngrdcol = 12' in reverse.read_text()
    assert set(a) == set(b)
    for key in a:
        assert a[key].tobytes() == b[key][::-1].tobytes()
    assert np.all(a['C1'] == 2.)
    # Independently check that all combinations were generated, not just one range.
    assert len(set(zip(a['C8'], a['C11']))) == 12


@pytest.mark.parametrize('defect', [False, True])
@pytest.mark.parametrize('separator', [[], ['--']])
def test_runner_detects_first_column_state_leak(tmp_path, monkeypatch, capsys, defect, separator):
    monkeypatch.setattr(mirror, 'build_test_executable', lambda log: tmp_path / 'clubb_standalone')
    actual_run = mirror.subprocess.run
    commands = []

    def run(command, **kwargs):
        if Path(command[1]).name == 'create_multi_col_params.py':
            return actual_run(command, **kwargs)
        commands.append(command)
        params = read_params(command[command.index('-params') + 1])
        values = params['C8'].copy()
        if defect:
            values += params['C8'][0]  # Accidental dependency on the first column.
        out = Path(command[command.index('-out_dir') + 1])
        out.mkdir()
        write_stats(out / 'rico_stats.nc', values)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(mirror.subprocess, 'run', run)
    assert mirror.main(['-case', 'rico', '-n', '4', '-out_dir', str(tmp_path),
                        '-override', 'microphys_scheme="morrison"', *separator, '-debug', '1']) == int(defect)
    assert len(commands) == 2
    for command in commands:
        assert Path(command[1]).name == 'run_scm.py'
        assert command[-1] == 'rico'
        assert command[command.index('-debug') + 1] == '1'
        assert command[command.index('-override') + 1] == 'microphys_scheme="morrison",l_lh_straight_mc=.true.'
    assert len(list(tmp_path.glob('run_*/forward.log'))) == 1
    output = capsys.readouterr().out
    stages = ['[1/4] Configuring and compiling', '[2/4] ABC run',
              '[3/4] CBA run', '[4/4] Comparison']
    assert [output.index(stage) for stage in stages] == sorted(output.index(stage) for stage in stages)
    assert output.count('Cases: rico') == 2
    assert 'Grid: -hr C8/0.2:0.8/4' in output
    assert 'C8 = [0.2, 0.4, 0.6, 0.8]' in output
    assert 'C8 = [0.8, 0.6, 0.4, 0.2]' in output
    assert f"{'FAIL' if defect else 'PASS'} rico" in output
    assert 'Summary for' not in output
    assert 'avg_abs_diff' not in output
    assert ('run_bindiff_all.py' in output) == defect
    if defect:
        command = shlex.split(output.splitlines()[-1])
        assert command[-2:] == ['-v', '2']
        assert Path(command[1]).name == 'run_bindiff_all.py'
        with netCDF4.Dataset(Path(command[3]) / 'rico_stats.nc') as aligned:
            np.testing.assert_allclose(aligned['rtm'][0], [.2+.8, .4+.8, .6+.8, .8+.8])
        assert 'avg_abs_diff' in next(tmp_path.glob('run_*/differences/comparison.log')).read_text()
        diagnostic = actual_run(command, capture_output=True, text=True)
        assert diagnostic.returncode == 1, diagnostic.stderr
        assert 'rtm' in diagnostic.stdout
        assert "Fail: 1 cases ['rico']" in diagnostic.stdout
    else:
        assert not list(tmp_path.glob('run_*/differences'))


def test_hypergrid_display_shows_actual_column_order(tmp_path, capsys):
    base = tmp_path / 'base.in'
    base.write_text('&clubb_params_nl\nC8=0.3\nC6rt=1.0\nC6thl=1.0\n/\n')
    spec = 'C8/0.2:0.8/3,C6rt=C6thl/1:2/2'
    forward, reverse = mirror.generate_parameter_files(tmp_path, spec, base)
    mirror.print_run_parameters(forward, spec)
    output = capsys.readouterr().out
    assert 'C8 = [0.2, 0.2, 0.5, 0.5, 0.8, 0.8]' in output
    for name in ('C6rt', 'C6thl'):
        assert f'{name} = [1, 2, 1, 2, 1, 2]' in output
    mirror.print_run_parameters(reverse, spec)
    output = capsys.readouterr().out
    assert 'C8 = [0.8, 0.8, 0.5, 0.5, 0.2, 0.2]' in output
    for name in ('C6rt', 'C6thl'):
        assert f'{name} = [2, 1, 2, 1, 2, 1]' in output


def test_diff_copies_align_raw_values_and_parameters_without_editing_originals(tmp_path, capsys):
    a, b = tmp_path / 'forward', tmp_path / 'reverse'
    a.mkdir()
    b.mkdir()
    write_stats(a / 'rico_stats.nc', [1, 2, 3])
    reverse_file = b / 'rico_stats.nc'
    write_stats(reverse_file, [3, 2, 1.5])
    with netCDF4.Dataset(reverse_file, 'r+') as ds:
        ds['rtm'].units = 'test units'
        ds.createVariable('col', 'i4', ('col',))[:] = [1, 2, 3]
        ds.createDimension('param', 1)
        ds.createVariable('clubb_params', 'f8', ('param', 'col'))[:] = [[.8, .5, .2]]
        ds.createVariable('packed', 'i2', ('col', 'time'), fill_value=-999)
        ds['packed'].scale_factor = .01
        ds['packed'][:] = np.ma.array([[3], [2], [1]], mask=[[False], [True], [False]])
    original = reverse_file.read_bytes()
    assert mirror.compare_directories(a, b, report_dir=tmp_path / 'report') == 1
    output = capsys.readouterr().out
    command = shlex.split(output.splitlines()[-1])
    with netCDF4.Dataset(Path(command[3]) / 'rico_stats.nc') as ds:
        np.testing.assert_array_equal(ds['rtm'][:], [[1.5, 2, 3]])
        np.testing.assert_array_equal(ds['col'][:], [1, 2, 3])
        np.testing.assert_array_equal(ds['clubb_params'][:], [[.2, .5, .8]])
        assert ds['rtm'].units == 'test units'
        assert ds['packed'][:].mask.tolist() == [[False], [True], [False]]
        ds.set_auto_maskandscale(False)
        np.testing.assert_array_equal(ds['packed'][:], [[100], [-999], [300]])
    assert reverse_file.read_bytes() == original


def test_missing_output_reports_case_failure_without_misleading_bindiff(tmp_path, capsys):
    a, b = tmp_path / 'forward', tmp_path / 'reverse'
    a.mkdir()
    b.mkdir()
    write_stats(a / 'rico_stats.nc', [1, 2, 3])
    assert mirror.compare_directories(a, b, report_dir=tmp_path / 'report') == 1
    output = capsys.readouterr().out
    assert 'FAIL rico' in output
    assert 'run_bindiff_all.py' not in output
    assert 'Missing CBA output' in (tmp_path / 'report' / 'comparison.log').read_text()


def test_failed_run_does_not_compare_stale_output(tmp_path, monkeypatch):
    monkeypatch.setattr(mirror, 'build_test_executable', lambda log: tmp_path / 'clubb_standalone')
    actual_run = mirror.subprocess.run
    commands = []

    def run(command, **kwargs):
        if Path(command[1]).name == 'create_multi_col_params.py':
            return actual_run(command, **kwargs)
        commands.append(command)
        return SimpleNamespace(returncode=7)

    monkeypatch.setattr(mirror.subprocess, 'run', run)
    monkeypatch.setattr(mirror, 'compare_directories', lambda *a: pytest.fail('compared after a failed run'))
    for _ in range(2):
        assert mirror.main(['-cases', 'rico,bomex', '-out_dir', str(tmp_path)]) == 1
    assert len(commands) == 2
    assert commands[0][commands[0].index('-cases') + 1] == 'rico,bomex'
    assert commands[0][commands[0].index('-nproc') + 1] == '2'
    assert len(list(tmp_path.glob('run_*'))) == 2


@pytest.mark.parametrize('values,name', [
    ([np.nan, 2, 1], 'rtm'),
    ([np.inf, 2, 1], 'rtm'),
    (np.ma.array([3, 2, 1], mask=[True, False, False]), 'rtm'),
    ([3, 2, 1], 'other_field'),
    ([3, 2], 'rtm'),
])
def test_invalid_output_fails(tmp_path, values, name):
    a, b = tmp_path / 'a.nc', tmp_path / 'b.nc'
    write_stats(a, [1, 2, 3])
    write_stats(b, values, name)
    assert mirror.check_file(a, b, 0., False)


def test_saved_output_comparison_and_tolerance(tmp_path):
    a, b = tmp_path / 'a', tmp_path / 'b'
    a.mkdir()
    b.mkdir()
    write_stats(a / 'rico_stats.nc', [1, 2, 3])
    write_stats(b / 'rico_stats.nc', [3, 2, 1])
    assert mirror.main([str(a), str(b)]) == 0
    write_stats(b / 'rico_stats.nc', [3, 2, 1.0001])
    assert mirror.main([str(a), str(b)]) == 1
    assert mirror.main(['-t', '0.001', str(a), str(b)]) == 0
    write_stats(b / 'extra_stats.nc', [3, 2, 1])
    assert mirror.main([str(a), str(b)]) == 1


@pytest.mark.parametrize('args', [
    ['-n', '1'], ['-hr', 'C8/1:1/3'], ['-hr', 'C8/0:1/1'],
    ['-hr', 'C8/0:nan/3'], ['-t', 'nan'], ['-nproc', '0'],
    ['-case', 'rico', '--', '-multicol', '5'],
    ['a', 'b', '-case', 'rico'],
])
def test_invalid_configuration_is_rejected(args):
    with pytest.raises(SystemExit) as error:
        mirror.parse_args(args)
    assert error.value.code == 2


@pytest.mark.parametrize('argv,expected', [
    (['-case', 'rico', '-debug', '1'], ['-debug', '1']),
    (['-debug', '1', '-case', 'rico', '-tout', '60', '-n', '4', '-nzmax', '128'],
     ['-debug', '1', '-tout', '60', '-nzmax', '128']),
    (['-exe', '/tmp/model with spaces', '-case', 'rico'], ['-exe', '/tmp/model with spaces']),
    (['-exe=/tmp/model', '-python'], ['-exe=/tmp/model', '-python']),
    (['-stats_tstart', '-1', '-stats_tend', '120', '-case', 'rico'],
     ['-stats_tstart', '-1', '-stats_tend', '120']),
    (['-case', 'rico', '-debug=1', '--', '-tout', '60'], ['-debug=1', '-tout', '60']),
    (['-n3', '-t1e-8', '-debug', '1'], ['-debug', '1']),
    (['-case', 'rico', '--', '-v'], ['-v']),
])
def test_run_options_are_forwarded_without_consuming_values_as_directories(argv, expected):
    args, forwarded = mirror.parse_args(argv)
    assert not args.directories
    assert forwarded == expected


@pytest.mark.parametrize('separator', [[], ['--']])
@pytest.mark.parametrize('option', ['-multicol=5', '-batch_size=2', '-all', '-min_cases',
                                   '-short_cases', '-priority_cases'])
def test_forwarding_cannot_override_mirror_setup(separator, option):
    with pytest.raises(SystemExit) as error:
        mirror.parse_args(['-case', 'rico', *separator, option])
    assert error.value.code == 2


@pytest.mark.parametrize('extra', [['-debug', '1'], ['--', '-exe', '/model']])
def test_saved_comparison_rejects_run_options(extra):
    with pytest.raises(SystemExit) as error:
        mirror.parse_args(['forward', 'reverse', *extra])
    assert error.value.code == 2


@pytest.mark.parametrize('forwarded', [
    ['-exe', '/custom/model'], ['-exe=/custom/model'],
    ['-install_dir', '/custom/install'], ['-driver_test'], ['-python'], ['-jax'],
])
def test_explicit_runtime_bypasses_build(tmp_path, monkeypatch, forwarded):
    monkeypatch.setattr(mirror, 'build_test_executable', lambda log: pytest.fail('unexpected build'))
    _, parsed = mirror.parse_args(['-case', 'rico', *forwarded])
    assert mirror.prepare_runtime(parsed, tmp_path / 'build.log') == forwarded


def test_saved_comparison_never_builds(tmp_path, monkeypatch):
    a, b = tmp_path / 'forward', tmp_path / 'reverse'
    a.mkdir()
    b.mkdir()
    write_stats(a / 'rico_stats.nc', [1, 2, 3])
    write_stats(b / 'rico_stats.nc', [3, 2, 1])
    monkeypatch.setattr(mirror, 'build_test_executable', lambda log: pytest.fail('unexpected build'))
    assert mirror.main([str(a), str(b)]) == 0


@pytest.mark.parametrize('failure_phase', [0, 1])
def test_build_failure_never_uses_stale_executable(tmp_path, monkeypatch, capsys, failure_phase):
    build = tmp_path / 'build'
    stale = build / 'src' / 'clubb_standalone'
    stale.parent.mkdir(parents=True)
    stale.write_text('stale executable')
    monkeypatch.setattr(mirror, 'BUILD_DIR', build)
    monkeypatch.setattr(mirror.shutil, 'which', lambda name: '/usr/bin/' + name)
    monkeypatch.setattr(mirror.platform, 'system', lambda: 'Linux')
    monkeypatch.setattr(mirror.platform, 'machine', lambda: 'x86_64')
    monkeypatch.setattr(mirror, 'generate_parameter_files', lambda *a: (tmp_path / 'a', tmp_path / 'b'))
    calls = []

    def run(command, **kwargs):
        calls.append(command)
        assert Path(command[0]).name == 'cmake', 'model ran despite failed build'
        kwargs['stdout'].write('diagnostic from cmake\n')
        return SimpleNamespace(returncode=9 if len(calls) == failure_phase + 1 else 0)

    monkeypatch.setattr(mirror.subprocess, 'run', run)
    assert mirror.main(['-case', 'rico', '-out_dir', str(tmp_path)]) == 1
    assert len(calls) == failure_phase + 1
    assert 'diagnostic from cmake' in capsys.readouterr().err


def test_missing_build_tool_has_actionable_error(tmp_path, monkeypatch):
    monkeypatch.setattr(mirror.shutil, 'which', lambda name: None)
    with pytest.raises(RuntimeError, match='CMake and gfortran on PATH'):
        mirror.build_test_executable(tmp_path / 'build.log')


def test_incremental_build_reuses_objects_and_picks_up_source_changes(tmp_path, monkeypatch):
    # Exercise real CMake with a tiny C target, avoiding a full CLUBB build.
    if not all(shutil.which(tool) for tool in ('cmake', 'gfortran', 'cc')):
        pytest.skip('CMake, gfortran, and a C compiler are required')
    source = tmp_path / 'source'
    (source / 'src').mkdir(parents=True)
    toolchains = source / 'cmake' / 'toolchains'
    toolchains.mkdir(parents=True)
    name = f'{platform.system().lower()}_{platform.machine().lower()}_gcc.cmake'
    (toolchains / name).write_text('# Minimal toolchain for the build-helper test.\n')
    (source / 'CMakeLists.txt').write_text(
        'cmake_minimum_required(VERSION 3.18)\nproject(probe C)\nadd_subdirectory(src)\n')
    (source / 'src' / 'CMakeLists.txt').write_text('add_executable(clubb_standalone main.c)\n')
    main = source / 'src' / 'main.c'
    contents = ('#include <stdio.h>\n#ifndef SILHS_MULTI_COL_RAND_DUPLICATE\n'
                '#error Missing SILHS test definition\n#endif\n'
                '#ifdef CONFIG_CHANGED\n#define MESSAGE "configured"\n'
                '#else\n#define MESSAGE "first"\n#endif\n'
                'int main(void) { puts(MESSAGE); return 0; }\n')
    main.write_text(contents)
    monkeypatch.setattr(mirror, 'REPO_ROOT', source)
    monkeypatch.setattr(mirror, 'BUILD_DIR', tmp_path / 'build')
    log = tmp_path / 'build.log'
    exe = mirror.build_test_executable(log)
    assert subprocess.check_output([str(exe)], text=True).strip() == 'first'
    wrapper = mirror.BUILD_DIR / 'toolchain.cmake'
    timestamps = (exe.stat().st_mtime_ns, wrapper.stat().st_mtime_ns)
    mirror.build_test_executable(log)
    assert (exe.stat().st_mtime_ns, wrapper.stat().st_mtime_ns) == timestamps
    main.write_text(contents.replace('first', 'second'))
    mirror.build_test_executable(log)
    assert subprocess.check_output([str(exe)], text=True).strip() == 'second'
    with (source / 'src' / 'CMakeLists.txt').open('a') as f:
        f.write('target_compile_definitions(clubb_standalone PRIVATE CONFIG_CHANGED)\n')
    mirror.build_test_executable(log)
    assert subprocess.check_output([str(exe)], text=True).strip() == 'configured'
