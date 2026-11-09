"""Configure and report the tuner's curated-function experiments.

The owned pytest wrapper loads these hooks for tuner tests. --tuner-trials
(default 10) controls repeated independent trials per function/method pair;
--tuner-seed (default 7) starts their consecutive deterministic seeds. These
are pytest plugin options, using pytest's supported double-dash syntax.

pytest_terminal_summary displays diagnostic global-minimum hit percentages,
improvement percentages, loss gaps and evaluation counts even with normal
pytest capture or -q. The tests record their own result data before assertions,
so failures to improve remain visible alongside global-hit measurements.
Incomplete trials retain their original exception and appear as ERROR here.
No algorithm settings or budgets are tuned by these hooks.
"""

import argparse
from statistics import fmean, median

import pytest


_RESULTS = pytest.StashKey[list]()


def _positive_trials(value):
    """Parse a positive trial count, rejecting invalid runs before test execution."""
    count = int(value)
    if count < 1:
        raise argparse.ArgumentTypeError('tuner trials must be >= 1')
    return count


def _nonnegative_seed(value):
    """Reject seeds unsupported by NumPy before any algorithm trial starts."""
    seed = int(value)
    if seed < 0:
        raise argparse.ArgumentTypeError('tuner seed must be >= 0')
    return seed


def pytest_addoption(parser):
    """Expose manual repetition/seed controls through pytest's plugin parser."""
    group = parser.getgroup('tuner curated-function benchmarks')
    group.addoption('--tuner-trials', type=_positive_trials, default=10,
                    help='Independent trials per curated function/algorithm (default: 10).')
    group.addoption('--tuner-seed', type=_nonnegative_seed, default=7,
                    help='First seed; trial i uses this value + i (default: 7).')


def pytest_configure(config):
    """Create one result list for the current pytest invocation."""
    config.stash[_RESULTS] = []


@pytest.fixture
def tuner_benchmark(request):
    """Supply shared reporting plus the configured trial count and first seed."""
    return {'results': request.config.stash[_RESULTS],
            'trials': request.config.getoption('--tuner-trials'),
            'seed': request.config.getoption('--tuner-seed')}


@pytest.hookimpl(trylast=True)
def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Print completed and failed benchmark metrics independently of stdout capture."""
    reports = config.stash[_RESULTS]
    if not reports:
        return
    terminalreporter.write_sep('=', 'Tuner curated functions: improvement and global-minimum hits')
    trials = config.getoption('--tuner-trials')
    seed = config.getoption('--tuner-seed')
    terminalreporter.write_line(f'Trials per pair: {trials}; seeds {seed}..{seed + trials - 1}. '
                               'Global-hit rules match Fortran G_unit; global percentages are diagnostic only.')
    terminalreporter.write_line('All methods require interface/numerical checks; optimizers additionally require at least one trial to lower the starting loss.')
    terminalreporter.write_line('Algorithms use registered defaults, batch size 8 and one worker; evaluation budgets differ.')
    headers = ['Function', 'Algorithm', 'Category', 'Global hits/runs', 'Global %', 'Improved %', 'Start loss', 'Best loss', 'Median gap', 'Avg evals', 'Result']
    rows = []
    for report in reports:
        completed = report['trials']
        hits = sum(trial['global_minimum_hit'] for trial in completed)
        percent = f'{100. * hits / len(completed):.1f}%' if completed else 'n/a'
        improved = f'{100. * sum(trial["improved"] for trial in completed) / len(completed):.1f}%' if completed else 'n/a'
        runs = f'{hits}/{len(completed)}'
        if len(completed) != report['requested_trials']:
            runs += f' of {report["requested_trials"]}'
        rows.append([
            report['objective'], report['algorithm'], report['category'], runs, percent, improved,
            f'{report["initial_loss"]:.6g}',
            f'{min(trial["best_loss"] for trial in completed):.6g}' if completed else 'n/a',
            f'{median(trial["gap"] for trial in completed):.6g}' if completed else 'n/a',
            f'{fmean(trial["evaluations"] for trial in completed):.0f}' if completed else 'n/a',
            report['status'],
        ])
    widths = [max(len(str(row[index])) for row in [headers, *rows]) for index in range(len(headers))]
    terminalreporter.write_line('  '.join(value.ljust(width) for value, width in zip(headers, widths)))
    terminalreporter.write_line('  '.join('-' * width for width in widths))
    for row in rows:
        terminalreporter.write_line('  '.join(value.ljust(width) for value, width in zip(row, widths)),
                                   red=row[-1] != 'PASS', green=row[-1] == 'PASS')
    terminalreporter.write_line('Global-minimum hit definitions (loss f; coordinates x):')
    definitions = {report['objective']: report['hit_condition'] for report in reports}
    for name, condition in definitions.items():
        terminalreporter.write_line(f'  {name}: {condition}')
    terminalreporter.write_line('Median gap = median(best loss - known global minimum). '
                               'ERROR means a trial did not finish; its original traceback is retained.')
