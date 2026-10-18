"""Check case-runner process reuse, asynchronous submission and cleanup.

Usage:
    tests/run_pytests.sh -unit -include_generated -k test_clubb_case_runner
    make_case() writes synthetic JSON inputs; make_runner() creates the matching
    worker. pytest runs the test_* checks of output, status and failure behavior.

Understanding:
    No scientific cases run. A small Python child mirrors the READY/path/DONE
    exchange in src/clubb_standalone.F90, deliberately writing large output,
    sleeping or exiting early to exercise pipe drainage, timeouts and closure.
    Separate checks cover one-shot exit status 6 and runtime-selection rejection.
"""

import json
import sys
import time
from pathlib import Path

import pytest

from utilities.clubb_case_runner import clubb_case_runner


WORKER = '''
import json, sys, time
from pathlib import Path
print('CLUBB_WORKER_READY', flush=True)
for line in sys.stdin:
    path = Path(line.strip())
    case = json.loads(path.read_text())
    if case.get('verbose'):
        for _ in range(1000):
            print('x'*1000, flush=True)
    time.sleep(case.get('delay', 0))
    if case.get('fail'):
        print('model diagnostic', flush=True)
        sys.exit(9)
    (path.parent/(path.stem+'.result')).write_text(str(case['value']))
    print('CLUBB_WORKER_DONE', flush=True)
'''


def make_case(tmp_path, name='case', **settings):
    """Write a synthetic case description and return its path."""
    path = tmp_path / (name + '.in')
    path.write_text(json.dumps({'value': 3, **settings}))
    return path


def make_runner(tmp_path, **kwargs):
    """Create an unstarted synthetic worker with captured output and a bounded timeout."""
    script = tmp_path / 'worker.py'
    script.write_text(WORKER)
    return clubb_case_runner(command=[sys.executable, str(script)], cwd=tmp_path,
                             persistent=True, log_file=tmp_path/'worker.log', timeout=3, **kwargs)


def test_lazy_runner_reuses_pid_for_different_complete_cases(tmp_path):
    """Submit two synthetic cases and require lazy startup, equal PIDs and distinct results."""
    with make_runner(tmp_path) as runner:
        assert runner.status()['state'] == 'new'
        assert runner.status()['pid'] is None
        first = runner.run_and_wait(make_case(tmp_path, 'first', value=11))
        second = runner.run_and_wait(make_case(tmp_path, 'second', value=27))
        assert first['state'] == second['state'] == 'completed'
        assert first['pid'] == second['pid']
        assert (tmp_path/'first.result').read_text() == '11'
        assert (tmp_path/'second.result').read_text() == '27'
    assert runner.status()['state'] == 'closed'
    with pytest.raises(RuntimeError, match='closed'):
        runner.run(make_case(tmp_path))


def test_run_is_nonblocking_and_output_drains_without_wait(tmp_path):
    """Produce more than a pipe buffer of output and require completion without calling wait."""
    with make_runner(tmp_path) as runner:
        runner.init()
        runner.run(make_case(tmp_path, verbose=True, delay=.2))
        assert runner.status()['state'] == 'running'
        deadline = time.monotonic() + 3
        while runner.status()['state'] == 'running' and time.monotonic() < deadline:
            time.sleep(.01)
        assert runner.status()['state'] == 'completed'
        assert (tmp_path/'worker.log').stat().st_size > 1_000_000
        assert runner.wait()['return_code'] == 0


def test_busy_runner_rejects_another_case_and_wait_timeout_is_explicit(tmp_path):
    """Keep one case active and require busy rejection and an explicit wait timeout."""
    with make_runner(tmp_path) as runner:
        path = make_case(tmp_path, delay=.2)
        runner.run(path)
        with pytest.raises(RuntimeError, match='active case'):
            runner.run(path)
        with pytest.raises(TimeoutError):
            runner.wait(.01)
        assert runner.wait()['state'] == 'completed'


def test_failure_is_reported_and_never_retried(tmp_path):
    """Fail a synthetic case and require its exit status, diagnostics and refusal of another run."""
    with make_runner(tmp_path) as runner:
        result = runner.run_and_wait(make_case(tmp_path, fail=True))
        assert result['state'] == 'failed'
        assert result['return_code'] == 9
        assert 'model diagnostic' in (tmp_path/'worker.log').read_text()
        with pytest.raises(RuntimeError, match='status 9'):
            runner.run(make_case(tmp_path, 'next'))
        assert not (tmp_path/'next.result').exists()


def test_one_shot_command_gets_the_namelist_without_worker_flags(tmp_path):
    """Run two one-shot commands and require only the namelist argument and native success code."""
    script = tmp_path/'single.py'
    script.write_text('import sys; from pathlib import Path; '
                      'Path(sys.argv[1]+".args").write_text(repr(sys.argv[1:])); sys.exit(6)')
    path = make_case(tmp_path)
    with clubb_case_runner(command=[sys.executable, str(script)],
                           persistent=False) as runner:
        assert runner.run_and_wait(path)['state'] == 'completed'
        assert runner.run_and_wait(path)['return_code'] == 6
    assert Path(str(path)+'.args').read_text() == repr([str(path)])


def test_close_terminates_active_process(tmp_path):
    """Close a sleeping worker and require the child to exit before close returns."""
    runner = make_runner(tmp_path)
    runner.run(make_case(tmp_path, delay=30))
    pid = runner.status()['pid']
    runner.close()
    assert not Path(f'/proc/{pid}').exists()
    assert runner.status()['return_code'] != 0
    runner.close()


@pytest.mark.parametrize('options', [
    {'python': True, 'jax': True},
    {'driver_test': True, 'exe': '/not/a/model'},
    {'gdb': True, 'jax': True},
])
def test_object_rejects_conflicting_runtime_options_before_launch(options):
    """Reject conflicting backend/debugger choices before creating a model process."""
    with clubb_case_runner(options) as runner:
        with pytest.raises(ValueError):
            runner.init()
        assert runner.status()['pid'] is None


def test_one_shot_runner_rejects_a_second_submission_while_busy(tmp_path):
    """Submit a slow one-shot case and require a second submission to fail immediately."""
    script = tmp_path / 'slow_case.py'
    script.write_text('import time; time.sleep(.2)')
    with clubb_case_runner(command=[sys.executable, str(script)], persistent=False) as runner:
        case = make_case(tmp_path)
        runner.run(case)
        with pytest.raises(RuntimeError, match='active case'):
            runner.run(case)
        assert runner.wait()['state'] == 'completed'


def test_run_requires_a_namelist_before_starting_a_runtime(tmp_path):
    """Reject a missing namelist without starting the worker."""
    with make_runner(tmp_path) as runner:
        with pytest.raises(ValueError, match='namelist'):
            runner.run(None)
        assert runner.status()['pid'] is None


def test_late_runtime_selection_preserves_one_shot_loss_defaults(tmp_path):
    """Supply the executable at init and require a one-shot loss run to accept native success."""
    script = tmp_path / 'loss_driver'
    script.write_text(f'#!{sys.executable}\nimport sys; sys.exit(6)\n')
    script.chmod(0o755)
    with clubb_case_runner(loss=True) as runner:
        runner.init({'exe': str(script)})
        assert runner.run_and_wait(make_case(tmp_path))['state'] == 'completed'
