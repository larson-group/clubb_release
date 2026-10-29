"""Explicit loss-window choices replace inherited tuning window defaults."""

from tuner.case_defaults import apply_case_overrides


def defaults():
    return dict(
        les_stats_file='truth.nc', altitude_comparison_range=[20., 2940.],
        time_average_range=[0, 600], average_time_seconds=300, num_time_windows=2,
    )


def test_explicit_window_count_replaces_inherited_average_interval():
    original = defaults()
    result = apply_case_overrides('bomex', original, {'num_time_windows': 1})
    assert result['num_time_windows'] == 1
    assert 'average_time_seconds' not in result
    assert result['time_average_range'] == [0, 600]
    assert original == defaults()


def test_explicit_average_interval_still_selects_window_count():
    result = apply_case_overrides('bomex', defaults(), {'average_time_seconds': 600})
    assert result['num_time_windows'] == 1
    assert result['average_time_seconds'] == 600


def test_duration_only_override_retains_inherited_average_interval():
    result = apply_case_overrides('bomex', defaults(), {'time_average_range': [0, 1200]})
    assert result['num_time_windows'] == 4
    assert result['average_time_seconds'] == 300
