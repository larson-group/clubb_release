"""Check native driver-test ordering and cleanup on JAX advance failures.

src/clubb_driver_test.F90 defines reinitialization and reset/rerun ordering.
Mock resources expose reuse after cleanup without running a model case. The
real counterpart and its NetCDF comparison run in the clubb_driver Jenkins job.
"""
import sys
from types import SimpleNamespace

import pytest

from clubb_jax.src import clubb_driver_test


@pytest.mark.parametrize('arguments', [[], [''], ['generated.in']])
@pytest.mark.parametrize('fail_advance', [False, True])
def test_reinitialization_rerun_and_cleanup(monkeypatch, arguments, fail_advance):
    events = []
    states = []

    def initialize(path):
        state = dict(total_param_sets=4, ngrdcol=2, closed=False,
                     err_info=SimpleNamespace(is_fatal=lambda: False))
        states.append(state)
        events.append(('initialize', path))
        return state

    def cleanup(state):
        state['closed'] = True
        events.append(('cleanup',))

    def reset(state, batch_num=None):
        assert not state['closed']
        events.append(('reset', batch_num))

    def advance(state, l_stdout, l_suppress_stats=False):
        assert not state['closed']
        assert l_stdout
        events.append(('advance', l_suppress_stats))
        if fail_advance:
            raise RuntimeError('advance failed')

    monkeypatch.setattr(sys, 'argv', ['clubb_driver_test', *arguments])
    monkeypatch.setattr(clubb_driver_test, 'init_clubb_case', initialize)
    monkeypatch.setattr(clubb_driver_test, 'clean_up_clubb', cleanup)
    monkeypatch.setattr(clubb_driver_test, 'set_case_initial_conditions', reset)
    monkeypatch.setattr(clubb_driver_test, 'advance_clubb_to_end', advance)
    if fail_advance:
        with pytest.raises(RuntimeError, match='advance failed'):
            clubb_driver_test.main()
    else:
        assert clubb_driver_test.main() == 6
    path = arguments[0] if arguments and arguments[0] else 'clubb.in'
    expected = [('initialize', path), ('cleanup',), ('initialize', path), ('advance', True)]
    if not fail_advance:
        expected += [('reset', None), ('advance', False), ('reset', 2), ('advance', False)]
    assert events == [*expected, ('cleanup',)]
    assert len(states) == 2 and all(state['closed'] for state in states)
