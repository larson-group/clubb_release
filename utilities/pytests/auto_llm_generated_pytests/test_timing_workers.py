"""Check reused timing workers and complete warmup/repetition boundaries.

Usage:
    tests/run_pytests.sh -unit -include_generated -k test_timing_workers
    cohort() creates two synthetic peers; pytest runs the scheduling and failure
    checks against the real case-runner process interface.

Understanding:
    No scientific cases run. Peers implement the native standalone's
    READY/path/DONE protocol, emit simple timing tables and have unequal delays.
    Assertions require overlapping runs, stable PIDs and completion of both
    peers before the next execution starts. Failure/timeout injection must close
    the whole group rather than strand another peer.
"""

import json
import time
import subprocess
import sys

import pytest

from utilities.time_clubb import BenchmarkOptions, ProcessGroup
from utilities.clubb_case_runner import clubb_case_runner
from utilities.pytests.timing_test_inputs import prepare_fake_case


RUNNER = r'''
import json, os, sys, time
from pathlib import Path
index = int(sys.argv[1])
print('CLUBB_WORKER_READY', flush=True)
for line in sys.stdin:
    namelist = Path(line.rstrip('\r\n'))
    payload = json.loads(namelist.read_text())
    directory = namelist.parent
    directory.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    if 'FAIL_PEER=true' in payload['argv']:
        if index == 1:
            raise SystemExit(9)
        time.sleep(30)
    time.sleep(0.03 + 0.1*index)
    (directory/'execution.json').write_text(json.dumps(
        {'pid':os.getpid(), 'start':started, 'end':time.monotonic()}))
    (directory/'case.timing').write_text("""
========================= CLUBB TIMER SUMMARY =========================
Backend: cpu_time
Timer Calls Inclusive(s) Exclusive(s) Avg Inc(s) Excl Pct
init_clubb_case 1 0.010000 0.010000 0.010000 25.00
advance_clubb_to_end 1 0.020000 0.020000 0.020000 50.00
clean_up_clubb 1 0.010000 0.010000 0.010000 25.00
Total exclusive timed cpu_time seconds: 0.040000
""")
    print('CLUBB_WORKER_DONE', flush=True)
'''


def cohort(tmp_path, forwarded=(), timeout=5):
    """Create two synthetic peers with unequal runtimes and optional failure injection."""
    runner = tmp_path / 'runner.py'
    runner.write_text(RUNNER)
    options = BenchmarkOptions('case', (2,), (1,), 1, 2, tmp_path, 'check', timeout, False)
    def factory(args, **kwargs):
        """Build one synthetic case runner for the peer selected by its log directory."""
        index = int(kwargs['log_file'].parent.name.split('_')[-1])
        return clubb_case_runner(command=[sys.executable, str(runner), str(index)],
                                 cwd=tmp_path, **kwargs)

    def prepare(args, directory, **kwargs):
        """Generate synthetic inputs, or exercise the preparation timeout when requested."""
        if 'PREPARE_TIMEOUT=true' in (args.override or []):
            subprocess.run([sys.executable, '-c', 'import time; time.sleep(30)'],
                           timeout=kwargs['timeout'], check=True)
        return prepare_fake_case(args, directory, **kwargs)

    return ProcessGroup(
        options=options, forwarded=list(forwarded), process_count=2,
        batch_size=1, group_dir=tmp_path/'startup', runner_factory=factory, prepare_case=prepare,
    )


def test_warmup_and_repetitions_keep_pids_and_wait_for_all_peers(tmp_path):
    """Require overlapping peers, stable PIDs and complete phase barriers across three runs."""
    histories = []
    with cohort(tmp_path) as group:
        for phase in ('warmup', 'measured_1', 'measured_2'):
            results, elapsed = group.run(tmp_path/phase)
            assert all(result.status == 'success' for result in results)
            assert all(
                {row.timer_name for row in result.timers} ==
                {'init_clubb_case', 'advance_clubb_to_end', 'clean_up_clubb'}
                for result in results
            )
            records = [
                json.loads(path.read_text())
                for path in sorted((tmp_path/phase).glob('*/execution.json'))
            ]
            assert max(r['start'] for r in records) < min(r['end'] for r in records)
            assert elapsed >= max(r['end'] for r in records)-min(r['start'] for r in records)
            histories.append(records)
    for index in range(2):
        assert len({phase[index]['pid'] for phase in histories}) == 1
    for previous, following in zip(histories, histories[1:]):
        assert max(r['end'] for r in previous) <= min(r['start'] for r in following)


@pytest.mark.parametrize('forwarded,timeout,status', [
    (['-override', 'FAIL_PEER=true'], 5, 'failed'),
    (['-override', 'PREPARE_TIMEOUT=true'], 0.1, 'timeout'),
])
def test_failed_or_timed_out_peer_cannot_strand_cohort(tmp_path, forwarded, timeout, status):
    """Inject a child failure or preparation timeout and require the entire cohort to close."""
    started = time.monotonic()
    with cohort(tmp_path, forwarded, timeout) as group:
        results, _ = group.run(tmp_path/'failed')
        assert any(result.status == status for result in results)
        assert group.broken
        assert all(worker['runner'].status()['state'] == 'closed' for worker in group.workers)
    assert time.monotonic()-started < 5
