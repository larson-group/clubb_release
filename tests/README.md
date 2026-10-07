# tests

This directory contains repo-level test harnesses and test-specific validators.
General run helpers live in `run_scripts/`; scripts here should be used when the
goal is to verify a specific behavior or regression.

Examples below assume they are run from the repo root. Most tests require CLUBB
to have already been compiled with the executable or Python API variant they
exercise.

## Pytest suites

Pytests are fast, specific checks. Actual SCM cases and full CLI/application
workflows belong in their component's `tests/` folder; shared workflows live here.
The [pytest workflow](../LLM_prompts/pytest_workflow.md) defines test quality,
admission and review standards.

Each suite has a `pytests/` directory with admitted modules in its parent and
provisional agent-written modules in `auto_llm_generated_pytests/`. New agent
coverage goes in that subdirectory. A human promotes useful coverage by moving
or merging it into the parent and removing the provisional copy. Fixtures and
shared input helpers stay with the suite owner. All existing API and JAX
pytests are initially provisional; this classification does not discard them.

`tests/run_pytests.sh` is the intended pytest entry point. Run from the repository root:

```sh
bash tests/run_pytests.sh -unit
bash tests/run_pytests.sh -dash
./compile.py -debug -python
bash tests/run_pytests.sh -api -include_generated
bash tests/run_pytests.sh -jax -include_generated
```

`-unit` selects utilities, tuner and `tests/pytests` harness contracts. `-dash`
prepares the shared Python environment with Dash dependencies. Its browser
callback checks also use Node. Provide Node 22 on `PATH`; Jenkins workers use
an installation in the Jenkins account's `$HOME/.local/bin`, prepared once
outside the checkout. The pipeline does not download tools or require npm.
`-api` requires a compatible build in `install/latest/python` (or
`CLUBB_F2PY_DIR`). `-jax`
prepares the CPU JAX environment and runs independent NumPy/analytic and
JAX contract checks. It requires no compiled Fortran library or Python API. `-all`
runs unit, Dash, API and JAX in that order. Pytest options such as `-q`, `-k` and
`--durations=10` are forwarded unchanged. JAX runs each module in a fresh
process so module-level flags and JIT caches cannot contaminate other modules.
Its per-module XML reports are saved under
`output/tests/pytests/jax` (`CLUBB_PYTEST_OUTPUT_DIR` selects another location).
Filters apply per module; a filter that selects nothing returns pytest's usual
exit code 5.

The wrappers exclude generated subdirectories by default; `-include_generated`
opts in for a review run. They print an explicit "no tests executed" message
for an empty admitted API/JAX suite; that is not a passing-test count.
Suite selection and admission belong to these entry points. There is no root
pytest configuration: bare `python -m pytest` uses pytest's own discovery and
does not apply the wrappers' environment setup or generated-test exclusion.

| Pytest owner | Jenkins production job | Provisional coverage |
| --- | --- | --- |
| `utilities/pytests`, `tuner/pytests`, `tests/pytests` | `clubb_python` unit stage | Set `INCLUDE_GENERATED_PYTESTS` for a review build |
| `clubb_python_api/pytests` | `clubb_python` F2PY stage | Explicitly included during the initial migration |
| `clubb_jax/pytests` | `clubb_jax` first testing stage | Explicitly included during the initial migration |
| `dash_app/pytests` | `clubb_dash` | Set `INCLUDE_GENERATED_PYTESTS` for a review build |

The legacy SensMatrix tests and their local configuration remain unchanged in
`utilities/sens_matrix/`. They are outside this migration and are not selected
by these Jenkins stages.

The corresponding `clubb_branch_*` jobs use the same Jenkinsfiles and a `BRANCH`
parameter. API/JAX's explicit provisional inclusion preserves the old API
coverage and adds the previously unrun JAX suite while human review proceeds.
CI execution never promotes a module. New maintained suites need a Jenkins
owner; provisional additions need an explicit review route through that owner.

The entry points prepare their owned Python environments automatically, except
when `CLUBB_PYTHON` explicitly selects a prepared interpreter. Compilation is
separate. Required CI dependencies/builds fail setup when missing; expected
optional dependencies use real pytest skips rather than printing "SKIP" and
returning. Runtime reports distinguish executed checks from skips.

### Component-specific workflows

Real-case checks live with their component: [JAX manual case checks](../clubb_jax/README.md#manual-case-checks),
[Dash ARM comparison](../dash_app/README.md#development-and-tests), and the
[API argument-contract audit](../clubb_python_api/tests/argument_list_enforcer/argument_contract_audit.py).
Run these explicitly; they are separate from the focused pytest stages in Jenkins.
The root `tests/` directory holds shared entry points, comparisons across
implementations and general CLUBB regressions.

## Test Scripts

### `check_budget_balance.py`

Runs a fixed list of SCM cases and then checks CLUBB budget closure with
`postprocessing/check_budgets_balance/checkBudget.py`.

Examples:

- `python3 tests/check_budget_balance.py`
  Runs the budget-balance suite for the current checkout.

- `python3 tests/check_budget_balance.py /path/to/clubb`
  Runs the same suite against an explicit CLUBB source tree.

### `check_mirrored_multi_col_output.py`

Tests column independence: reordering columns should not change their
individual results. The script generates distinct parameter columns, runs
the same case in forward and reverse order, and compares matching columns.

```sh
python3 tests/check_mirrored_multi_col_output.py -case rico -multicol 3
```

`-multicol 3` varies `C8` evenly from 0.2 to 0.8:

| Run | First column | Second column | Third column |
| --- | --- | --- | --- |
| Forward (ABC) | A: `C8 = 0.2` | B: `C8 = 0.5` | C: `C8 = 0.8` |
| Reverse (CBA) | C: `C8 = 0.8` | B: `C8 = 0.5` | A: `C8 = 0.2` |

Column A must match A across runs, and likewise for B and C; the three columns should
differ from each other. This exposes accidental dependencies on column order,
such as reading column `1` instead of `i`, or carrying temporary values
between columns. Any bug that could cause information from one column to
infect other columns should cause this test to fail, because changing
the column order should change which column is the infectious one.

Use `-multicol PARAM/MIN:MAX/NPOINTS` for custom parameter ranges. Multiple ranges
form a grid; this example generates six columns, then reverses their order:

```sh
python3 tests/check_mirrored_multi_col_output.py -case rico -multicol 'C8/0.2:0.8/3,C11/0.2:0.8/2'
```

The runner automatically configures and builds a dedicated gfortran CPU
executable with the required SILHS sampling settings. Rerunning reuses the
build and rebuilds changed sources as needed. Have CMake, gfortran, and the
usual CLUBB build dependencies available; no manual toolchain edits are needed.

Omit `-case` to run the standard case set, or select several with
`-cases rico,rico_silhs,mc3e`. Defaults are three columns, 200 timesteps
(`-max_iters`), and half the available logical CPUs for concurrent cases (`-workers`). Use `-config` / `-params_file`
for alternate base parameters. Additional `run_scm.py` options are forwarded,
for example `-debug 1`; the `--` separator is optional. See `-help` for options and
the [script header](check_mirrored_multi_col_output.py) for build and diagnostic
details.

### `run_G_unit_tests.py`

Writes a temporary `G_unit_tests.in` namelist and runs the compiled
`install/latest/G_unit_tests` executable.

Examples:

- `python3 tests/run_G_unit_tests.py`
  Runs the built-in default G-unit test set.

- `python3 tests/run_G_unit_tests.py -all`
  Enables every G-unit test flag.

- `python3 tests/run_G_unit_tests.py -KK_unit_tests`
  Runs only the KK unit tests.

- `python3 tests/run_G_unit_tests.py -smooth_heaviside_test -smooth_min_max_test`
  Runs only the selected smooth-function tests.

### `run_benchmark_converter_test.py`

Smoke-tests the benchmark normalizer by creating a tiny SAM-like NetCDF file,
converting it, and checking aliases and formulas.

Example:

- `python3 tests/run_benchmark_converter_test.py`
  Runs the self-contained benchmark-converter smoke test.

### `run_bindiff_w_flags.py`

Clones two or more git refs, compiles each clone, runs
`run_scripts/run_clubb_w_varying_flags.py` in each clone, and compares the
resulting output trees with `run_scripts/run_bindiff_all.py -flag_sets`.

Examples:

- `python3 tests/run_bindiff_w_flags.py -branches master,my_branch -d /tmp/clubb_bindiff`
  Compares `master` and `my_branch` using the default core flag set.

- `python3 tests/run_bindiff_w_flags.py -overwrite_existing -branches master,my_branch -d /tmp/clubb_bindiff -max_iters 360`
  Recreates existing clone directories without prompting and forwards
  `-max_iters 360` to the case runs.

- `python3 tests/run_bindiff_w_flags.py -branches master,my_branch -flag_config_file input/flag_sets/run_bindiff_w_flags_config_example.json -priority_cases -workers 4`
  Uses an explicit flag config and forwards case-selection and worker-count
  options to `run_clubb_w_varying_flags.py`.

- `python3 tests/run_bindiff_w_flags.py -branches old_ref,new_ref -no_compile -skip_default_flags`
  Reuses existing builds and compares only alternate flag sets.

### `run_jax_comparison_mutation_test.py`

Checks that the real JAX-versus-Fortran comparison can detect model errors.
Runs BOMEX for four 60-second timesteps and four columns, using
`input/stats/multi_col_stats.in` and the comparison harness's default tolerances.
The numerical mutations use independent physical-signal checks with the
effective tolerances in the comparison report.
`CASE` and `TIMESTEPS_TO_RUN` are configured near the top of the script;
the current mutation sites and run length are validated for BOMEX.
An unmodified control must match, followed by four independent JAX mutations:

- `heating`: replace zero heating with `1e-6 K/s` (about `6e-5 K` per step).
- `parameter_handling`: read C8 `[0.2, 0.4, 0.6, 0.8]` as
  `[0.2, 0.3, 0.4, 0.5]`, preserving column 1 and all other parameters.
- `first_column`: use column 1's previous `wp3` in every column's timestep
  equation, simulating a `(1,k)` instead of `(i,k)` indexing error.
- `missing_stat`: omit `wprtp` from the JAX NetCDF output while still
  calculating it, so strict bindiff must detect a missing variable.

Each original/replacement pair is defined and briefly explained at the top of
the script. Each mutation must make the comparison exit 1 with both models
completing all four steps. The three numerical mutations must report their
expected physical field (`thlm` or `wp3`) above threshold. For the two column
mutations, column 1 must still match while each later column fails; a
parameter-metadata difference alone does not count. The missing-stat mutation
must report exactly `wprtp` absent from JAX output, with every other saved
variable unchanged. Saved C8 values are checked independently. A crash,
missing/empty/nonfinite required physical output, or a passing mutated
comparison fails this test. The final summary shows each mutation's result
and the output directory.

Run from the repository root:

    python3 tests/run_jax_comparison_mutation_test.py

Uses CPU/double precision and the existing selected/latest Fortran install;
build/install that executable first. The normal JAX launcher prepares its
environment. Five JAX compilations are needed even though the runs are short.
Source copies, commands, NetCDF files, logs and `mutation_test_summary.json`
are retained in a new `output/tests/jax_comparison_mutation_*` directory.
`-output_dir PATH` selects a new destination; `-timeout SECONDS` adjusts the
300-second limit per comparison. `-mutations parameter_handling first_column`
runs only those mutations, sharing one control. Working source and normal
comparison outputs are preserved. These are representative numerical and
output-schema errors through the entire runner, not exhaustive coverage of
every field or scheme.

### `run_clubb_conv_test.py`

Runs one case at several timesteps and checks convergence for a selected output
variable.

Examples:

- `python3 tests/run_clubb_conv_test.py`
  Runs the default BOMEX convergence check for `rcm`.

- `python3 tests/run_clubb_conv_test.py -case rico -var cloud_frac`
  Checks convergence of `cloud_frac` for `rico`.

- `python3 tests/run_clubb_conv_test.py -plot_result -case bomex -var rcm`
  Runs the check and writes a convergence plot plus final-profile comparisons
  for every timestep under `output/`.

- `python3 tests/run_clubb_conv_test.py -config Lscale`
  Forwards unrecognized options to every `run_scm.py` invocation, allowing the
  convergence test to use a parameter-and-flag configuration.

### `run_jax_vs_fortran_cases.py`

Runs selected cases with both the JAX driver and the Fortran standalone driver,
then compares outputs with `run_bindiff_all.py`.

Every case runs once per flag set. The unmodified "default" flag set is always
included for every run. Results are written as

    output/tests/jax_driver_test_results/
      jax_output/<flag set>/*.nc
      fortran_output/<flag set>/*.nc
      logs/<flag set>/<case>_{run_jax,run_fortran,bindiff}.log
      logs/<flag set>/<case>_bindiff.json
      case_compare_summary.json
      final_bindiff.log

so the two output roots can be compared directly with
`run_bindiff_all.py -flag_sets`, which is what the final combined diff does.

`DEFAULT_CASES` in the script defines the case list and any per-case overrides.
Pass `run_scm.py` options such as `-max_iters`, `-dt_main`, `-stats`, and
`-debug` directly; they apply to both models and override curated case settings.
Effective step and timestep settings and forwarded options are recorded in the
results JSON. The harness controls driver selection, output paths, columns, and
case/flag overrides.

Examples:

- `python3 tests/run_jax_vs_fortran_cases.py -cases bomex -workers 1`
  Runs a single serial JAX-vs-Fortran comparison for easier debugging.

- `python3 tests/run_jax_vs_fortran_cases.py -cases bomex atex -max_iters 3`
  Runs two cases with a short iteration limit.

- `python3 tests/run_jax_vs_fortran_cases.py -bindiff_verbose 2 -bindiff_threshold 1e-12`
  Runs the default case set and prints detailed strict bindiff output.

- `python3 tests/run_jax_vs_fortran_cases.py -flag_config_file input/flag_sets/run_bindiff_w_flags_config_example.json`
  Also runs every case under each flag set in the JSON file, using the same
  config format as `run_scripts/run_clubb_w_varying_flags.py`.

### `run_loss_output_consistency.py`

Runs normal CLUBB and the loss driver for one case, then checks that loss-driver
stats output and printed metrics are consistent. The `-jax` option selects JAX for both
runs. This check requests one full-window record even when tuning defaults
use several subwindows. Profiles and printed metrics retain existing criteria.
Jenkins adds a CPU JAX loss consistency stage.


Examples:

- `python3 tests/run_loss_output_consistency.py bomex`
  Runs the default consistency check for `bomex`.

- `python3 tests/run_loss_output_consistency.py arm -fields cloud_frac rcm`
  Checks selected fields for `arm`.

- `python3 tests/run_loss_output_consistency.py bomex -output_root output/loss_check`
  Writes all generated output under a custom root directory.

- `python3 tests/run_loss_output_consistency.py bomex -config default`
  Runs with an explicit named tunable config.

### `run_python_vs_fortran_cases.py`

Runs selected cases with both the Python standalone driver and the Fortran
standalone driver, then compares outputs with `run_bindiff_all.py`.

Examples:

- `python3 tests/run_python_vs_fortran_cases.py -cases bomex -workers 1`
  Runs one serial Python-vs-Fortran comparison for debugging.

- `python3 tests/run_python_vs_fortran_cases.py -cases bomex atex -max_iters 3`
  Runs two cases with a short iteration limit.

- `python3 tests/run_python_vs_fortran_cases.py -bindiff_verbose 2 -bindiff_threshold 1e-12`
  Runs the default case set with detailed strict bindiff output.

- `python3 tests/run_python_vs_fortran_cases.py -keep_existing`
  Reuses existing comparison output directories.

### `run_restart_test.py`

Runs a full case, moves its output into `restart/`, reruns from a restart time,
and compares the final timestep of a selected variable bit-for-bit.

Examples:

- `python3 tests/run_restart_test.py bomex`
  Runs the restart test for `bomex`, comparing `thlm`.

- `python3 tests/run_restart_test.py rico_silhs -var rcm`
  Runs the restart test for `rico_silhs`, comparing `rcm`.

- `python3 tests/run_restart_test.py bomex -keep_artifacts`
  Keeps generated `output/` and `restart/` files after the test.

Use `-jax` or `-jax=cpu` for both restart integrations. The JAX check compares
all saved columns; native first-column behavior is unchanged. The restart uses
a saved interior record nearest the effective run midpoint, including odd
shortened runs. Jenkins adds CPU normal and SILHS sequence restart stages.

### `run_silhs_test.py`

Checks SILHS convergence by comparing a small-sample run against a large-sample
run for one case.

Examples:

- `python3 tests/run_silhs_test.py`
  Runs the default SILHS convergence check.

- `python3 tests/run_silhs_test.py -case rico_silhs -n_small 8 -n_large 1000`
  Runs an explicit case and sample-count pair.

- `python3 tests/run_silhs_test.py -stats input/stats/all_stats.in -show_output`
  Uses an explicit stats file and prints the underlying `run_scm.py` output.

- `python3 tests/run_silhs_test.py -keep_outputs`
  Leaves the small and large output directories in place for inspection.

### `run_stats_output_consistency.py`

Runs one case several ways and verifies unified stats output is consistent
across batch sizes, output intervals, and stats windows.

Examples:

- `python3 tests/run_stats_output_consistency.py`
  Runs the default BOMEX stats-output consistency suite.

- `python3 tests/run_stats_output_consistency.py arm -stats input/stats/standard_stats.in`
  Runs the suite for `arm` with an explicit stats registry.

- `python3 tests/run_stats_output_consistency.py bomex -batch_sizes 4,2,1 -coarse_touts 300,600`
  Checks the default multicol setup with explicit batch sizes and coarse output
  intervals.

- `python3 tests/run_stats_output_consistency.py bomex -window_start 7200 -window_end 14400`
  Uses an explicit stats window for the windowing checks.

- `python3 tests/run_stats_output_consistency.py bomex -config default -debug 0`
  Forwards a named config and debug level to `run_scm.py`.

### `run_thread_test.py`

Runs the thread-safety regression test using the compiled standalone and thread
test executables.

Examples:

- `python3 tests/run_thread_test.py`
  Runs with the default OpenMP thread count.

- `python3 tests/run_thread_test.py -threads 4`
  Runs with four OpenMP threads.

### `run_timestep_tests.py`

Runs the standard SCM case list repeatedly with increasing `dt_main` and
`dt_rad` values while stats output is disabled.

Example:

- `python3 tests/run_timestep_tests.py`
  Runs the timestep sweep for the built-in case list.

### `test_fatal_error_handling.py`

Temporarily patches tagged lines in CLUBB source to inject NaNs, recompiles, and
checks that CLUBB reports the expected fatal errors. The script restores edited
files on normal exit and signal handling.

Example:

- `python3 tests/test_fatal_error_handling.py`
  Runs all configured fatal-error injection tests.

### `test_fire_tuner.py`

Configures the FIRE case for tuning, writes a focused FIRE stats file, and runs
`run_scripts/run_tuner.py`.

Examples:

- `python3 tests/test_fire_tuner.py`
  Runs the FIRE tuner test against the current checkout.

- `python3 tests/test_fire_tuner.py /path/to/clubb`
  Runs the FIRE tuner test against an explicit CLUBB source tree.

### `test_monoflux_limiter_GPU.py`

Runs the monotonic turbulent flux limiter GPU-vs-CPU PCAST test. It temporarily
patches the NVHPC toolchain and source marker, compiles with OpenACC PCAST
settings, runs a multi-column `mc3e` case, checks for PCAST differences, and
restores edited files.

Example:

- `python3 tests/test_monoflux_limiter_GPU.py`
  Runs the full GPU PCAST regression test.

Branch bindiff accepts the current option names and translates them for cloned
revisions that still expose the older interfaces. This preserves comparisons
against master and older refs without retaining old aliases in current scripts.
