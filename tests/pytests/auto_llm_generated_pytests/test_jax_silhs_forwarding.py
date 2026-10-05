"""The comparison harness applies each SILHS configuration to both model commands."""

import pytest
from tests import run_jax_vs_fortran_cases as harness
from tests.pytests.jax_comparison_inputs import run_comparison


@pytest.mark.parametrize("case", ["rico_silhs", "lba_kk_silhs", "lba_silhs"])
def test_silhs_overrides_and_thresholds_forwarded_to_both_runs(tmp_path, case):
    config = next(config for config in harness.DEFAULT_CASES if config.case == case)
    result, commands, task = run_comparison(tmp_path, config)
    jax, fortran, bindiff = commands
    assert result.status == "match"
    assert "-strict" in bindiff
    assert bindiff[bindiff.index("-percent_threshold") + 1] == str(config.percent_threshold)
    assert jax.count("-max_iters") == fortran.count("-max_iters") == 1
    assert jax.count("-dt_main") == fortran.count("-dt_main") == 1
    start = jax.index("-multicol") + 2
    assert task.run_scm_args == jax[start:start + len(task.run_scm_args)]
    assert task.run_scm_args == fortran[start:start + len(task.run_scm_args)]
    assert "-jax" in jax and "-jax" not in fortran
    assert jax[-1] == fortran[-1] == config.source_case
    override = jax[jax.index("-override") + 1]
    assert override == fortran[fortran.index("-override") + 1]
    assert "l_lh_deterministic_test=.true." in override
    assert "l_lh_importance_sampling=.false." in override
    assert "l_random_k_lh_start=.false." in override
