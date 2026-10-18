## Utilities Directory Contents ## 

`case_json_to_namelist.py` resolves `input/case_setups/case_definitions.json`
and local saved-case snapshots into model namelists. For example,
`python utilities/case_json_to_namelist.py bomex -output_file /tmp/bomex_model.in`.
`create_case_namelist.py` looks for `NAME_model.in` first and uses JSON when no
model file exists. Saved custom cases are run by name with the normal CLI. See
[case setups](../input/case_setups/README) for inheritance and custom-case rules.
The converter serializes values without enforcing model rules or checking assets.

The [CLUBBStandardsCheck.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/CLUBBStandardsCheck.py) python script can be used to check Fortran 
source files to determine if they follow certain good software engineering 
practices which are meant to be enforced for CLUBB source files. 

[check_for_errors.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/check_for_errors.py) checks CLUBB source code files for common mistakes or
CLUBB standards violations. It checks for uninitialized output variables,
magic flags, and magic numbers. An output variable is only initialized if it
is set in every part of an if or select case statement or if it is set outside
of these statements. A magic number is any number that appears
in a subroutine or function call or appears in an equation on the right-hand
side of an equals sign. Integers from -6 to +6 are not considered magic numbers,
and nothing in an if statement is considered a magic number. If the argument -w
or --show-warnings is provided, warnings will be printed when a variable is set
by a subroutine or function call. If no files are provided in the argument list,
the script will check all .F90 files under the src directory (recursively).
USAGE: "python3 check_for_errors.py [-w or --show-warnings] [<filename>.F90 ...]"

[print_tunable_parameters_table.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/print_tunable_parameters_table.py) prints out a table of the tunable parameters
with their values in markdown format. The output can be copied and pasted directly
into a github comment and will be formatted as a table there. If multiple files
are listed as arguments, then the output will have a column for each file listing
the values of the parameters if they are present or an indicator if commented out.
USAGE: "python print_tunable_parameters_table.py [options]
               <list of tunable_parameters.in files>"

Other files in this directory:

[__init__.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/__init__.py) makes utilities/ a Python package, so other code can use
"from utilities.<module> import ...".

[benchmark_converter.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/benchmark_converter.py) converts SAM/COAMPS LES benchmark NetCDF files into a
small NetCDF that uses CLUBB stats names and CLUBB-style time. Used by the
tuner/loss driver (through [create_case_namelist.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/create_case_namelist.py)) and the Dash plot and
tune tabs.

[clubb_settings_validation.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/clubb_settings_validation.py) is a hand-maintained Python copy of CLUBB's
parameter bounds and flag/parameter compatibility rules, so Dash and the tuner
can check settings without a compiled CLUBB. [tests/run_clubb_settings_validation_test.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/tests/run_clubb_settings_validation_test.py)
checks it against the compiled Fortran.

[convert_acc_to_omp.bash](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/convert_acc_to_omp.bash) adds OpenMP directives next to the OpenACC directives
in [src/](https://github.com/larson-group/clubb/tree/bc8d70735e710ef36f81e42f58843908ae3ddcbc/src), using Intel's OpenACC-to-OpenMP migration tool (cloned on first use).
Used by the nvhpc GPU-vs-CPU Jenkins test.

[convert_to_async.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/convert_to_async.py) adds async(1)/wait clauses to selected OpenACC directives
in [src/](https://github.com/larson-group/clubb/tree/bc8d70735e710ef36f81e42f58843908ae3ddcbc/src). Used by the GPU async output-matching Jenkins test.

[create_case_namelist.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/create_case_namelist.py) combines a case's model file, tunable parameters,
flags, SILHS settings, and stats list into the single <case>.in namelist that
CLUBB reads, applying any overrides. Called by [run_scripts/run_scm.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/run_scripts/run_scm.py),
[run_scm_loss.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/run_scripts/run_scm_loss.py), and the tuner.

[create_multi_col_params.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/create_multi_col_params.py) writes a multicolumn tunable-parameter file
(duplicated, perturbed, mirrored, or hypergrid columns). Called by
[create_case_namelist.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/create_case_namelist.py), the tuner, and several GPU Jenkins tests.

[flag_sets.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/flag_sets.py) reads a flag-set JSON file (see [input/flag_sets/](https://github.com/larson-group/clubb/tree/bc8d70735e710ef36f81e42f58843908ae3ddcbc/input/flag_sets)) and turns each
set into a namelist override string. Shared by
[run_scripts/run_clubb_w_varying_flags.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/run_scripts/run_clubb_w_varying_flags.py) and [tests/run_jax_vs_fortran_cases.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/tests/run_jax_vs_fortran_cases.py).

[les_chi_moments.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/les_chi_moments.py) derives chi (extended liquid water) moments from LES
fields using CLUBB's thermodynamic constants and saturation formula. Used by
[benchmark_converter.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/benchmark_converter.py) and [sam_3d_reference.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/sam_3d_reference.py).

[loss_metrics.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/loss_metrics.py) is a Python copy of the profile loss and Taylor metrics
computed in [src/clubb_loss_driver.F90](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/src/clubb_loss_driver.F90). Used by the Dash plot tab and the
loss-driver tests. It also owns the shared loss-table formatter used by the
Fortran Python API and JAX standalone loss front ends.

[output_paths.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/output_paths.py) turns a user-given output name into a directory under output/
(and rejects ".."). Shared by [run_scm.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/run_scripts/run_scm.py), [run_scm_loss.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/run_scripts/run_scm_loss.py), and Dash.

[sam_3d_reference.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/sam_3d_reference.py) reads raw 3D SAM snapshots from the shared benchmark-run
directory for Dash diagnostics and experimental PDF analyses.

[save_tunable_config.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/save_tunable_config.py) copies a config directory under
[input/parameter_and_flag_configs/](https://github.com/larson-group/clubb/tree/bc8d70735e710ef36f81e42f58843908ae3ddcbc/input/parameter_and_flag_configs) and applies Run/Tune-tab overrides to the
copy. Used only by Dash (it imports [dash_app](https://github.com/larson-group/clubb/tree/bc8d70735e710ef36f81e42f58843908ae3ddcbc/dash_app)).

[split_stats_to_legacy.py](https://github.com/larson-group/clubb/blob/bc8d70735e710ef36f81e42f58843908ae3ddcbc/utilities/split_stats_to_legacy.py) splits one combined stats NetCDF file into
legacy-style per-grid files (e.g. <case>_zt.nc, <case>_zm.nc). Nothing in the
repository calls it.

### Reusable case runners

[`clubb_case_runner`](clubb_case_runner.py) lets a Python script select a native
Fortran, Python or JAX runtime and run complete cases. Keeping the standalone
process alive lets JAX reuse compiled programs in memory. Each request still
initializes and cleans up its model; model allocations are not retained.
[`run_scm.py`](../run_scripts/run_scm.py) uses this interface for one case, and
[`time_clubb.py`](time_clubb.py) uses it for repeated timing runs.

For example, run this from a repository Python script. It generates two short
BOMEX inputs in separate output directories and runs both in the same JAX
process:

```python
from pathlib import Path
from utilities.clubb_case_runner import clubb_case_runner, parse_case_arguments
from utilities.create_case_namelist import create_case_namelist_from_args

options = parse_case_arguments(["bomex", "-jax=cpu", "-max_iters", "20"])
with clubb_case_runner(options, echo=True) as runner:
    for index in range(2):
        output = Path(f"output/runner_example/run_{index}").resolve()
        namelist = create_case_namelist_from_args(options, output)
        result = runner.run_and_wait(namelist, log_file=output / "bomex_log")
        if result["state"] != "completed":
            raise RuntimeError(result["error"])
```

Omit `-jax=cpu` to select the native executable, or use `-python` for the Python
SCM. Native/Python runs require a current corresponding build; JAX's
[launcher](../clubb_jax/run_jax.py) prepares its own environment. By default,
native selection uses `install/selected` if present, otherwise `install/latest`.
`-install_dir` chooses another install. The shared
[namelist utility](create_case_namelist.py) owns model options and overrides;
`create_case_namelist_from_args` takes a resolved output directory and returns
the generated file path.

The runner also accepts a mapping of runtime options, such as
`{"jax": True, "jax_options": "cpu"}`. Its lifecycle is:

| Method | Behavior |
| --- | --- |
| `init(runtime_options=None)` | Resolve the runtime; start a persistent worker and wait for readiness. Optional before the first run. |
| `run(namelist, log_file=None)` | Initialize on first use, submit an existing aggregate namelist file and return the runner. |
| `wait(timeout=None)` | Wait for the submitted case and return its status dictionary. |
| `run_and_wait(namelist, log_file=None, timeout=None)` | Submit a case and return the same status as `wait()`. |
| `status()` | Return the latest status without waiting. |
| `close()` | Release the process and log. A context manager also calls this on exceptions. |

Construction and context entry do not start a process. A first `run()` may wait
for runtime startup; after submission it does not wait for the case. Submit one
case at a time. Runtime selection is fixed after initialization. `log_file`
appends console output, and `echo=True` also prints it.
A background reader drains non-debugger output even while the caller only polls
`status()`.

`status()` includes `state`, `pid`, `namelist`, `return_code`, `error`, `log_file`
and `finished_at`. States are `new`, `starting`, `ready`, `running`, `completed`,
`failed` and `closed`. A completed worker can still be alive. `finished_at` is a
monotonic timestamp in seconds, not a duration. Model failures return
`state="failed"`; startup/submission errors raise. The constructor's `timeout`
sets the default startup/wait limit in seconds; a wait timeout raises
`TimeoutError` and leaves
the case active for another wait or `close()`. Failed cases are not retried.

Ordinary standalones default to process reuse. `persistent=False` starts a new
process for every case. Explicit commands, `-exe`, `-gdb`, driver tests and loss
frontends use that one-shot mode by default; their `init()` resolves a command
and `run()` launches it. Closing an idle worker requests normal EOF shutdown;
closing an active case terminates its process group. Replace a closed or failed
runner before submitting another case.

`loss=True` selects a one-shot frontend for
[`run_scm_loss.py`](../run_scripts/run_scm_loss.py), which prints its loss table.
Parameter arrays and retained allocations remain in the tuner's existing F2PY
or JAX loss API (`init_clubb_loss`, `clubb_get_loss_for_params`,
`finalize_clubb_loss`). They do not pass through this ordinary case interface.

The [case runner lifecycle test](../tests/case_runner_lifecycle_test.py) checks
process reuse and exact numeric stats against fresh-process controls. It runs
Fortran, Python and CPU JAX in the first `clubb_driver` Jenkins stage. See the
[test instructions](../tests/README.md#case_runner_lifecycle_testpy) for the
case sequence, build requirements and standalone commands.

### Timing sweeps

[time_clubb.py](time_clubb.py) measures complete CLUBB runs at selected process
counts and column counts per process. It writes a profile that can be opened
in the Dash Profile tab. From the repository root, with a current native build:

```bash
./utilities/time_clubb.py bomex
```

The defaults are one process, one column, one full-case warmup and three measured
repetitions, with stats disabled and debug 0. To compare several workloads:

```bash
./utilities/time_clubb.py bomex -process_counts 1,2 -batch_sizes 1,4 \
    -warmups 1 -repeats 3 -max_iters 20 \
    -output_root output/timing -name bomex_short
```

A workload is one process-count/batch-size pair; the example runs all four
pairs. `-process_counts` selects concurrent model processes, while `-batch_sizes`
selects columns advanced together by each process. `-max_iters 20` deliberately
limits this example to 20 steps; omit it for the configured case duration.
Runtime/model options use the same parsers as `run_scm.py`. The timing script
assigns `-multicol`, `-batch_size` and per-run output directories, so those
options cannot be forwarded. Negative debug levels disable timers and are
rejected here. Timing requires ordinary standalones with the worker entry point;
an explicit `-exe` must implement that protocol too. Debugger and driver-test
modes cannot be used for this persistent timing sweep.

Each workload keeps its workers alive through all full-case warmups and
repetitions, with fresh initialization and cleanup each time. Every peer must
finish before the next run is submitted, including the transition from warmup
to measurements. Peers are submitted sequentially without waiting between
them; this coordinates phases rather than guaranteeing identical start times.

The script prints the profile directory under `output/timing` by default.
`profile.json` records settings and provenance, `batches.csv` records workloads,
and `timings.csv` holds raw `startup`, `warmup` and `measured` observations.
Representative inputs, console output and `.timing` files are saved under
`logs/`; other per-execution outputs are temporary. Named `.timing` regions keep
their backend's clock (CPU time or GPTL wall time). Process-group wall time
measures submission through all case completions, including model setup and
cleanup, but excluding worker startup, input preparation and inter-run waiting.
Startup is recorded separately. Older profiles included new-process startup in
each repetition; `benchmark.execution_model=persistent_workers` identifies this
methodology.

`-timeout SECONDS` bounds startup, input preparation and execution stages;
it is not a single deadline for the whole sweep. A failed or timed-out peer
stops the workload. `-continue_on_error` proceeds to the next workload without
retrying the failed one. Existing profile names get a version suffix unless
`-overwrite` is selected, which replaces that named profile.

JAX `.timing` output is added separately in
[#1393](https://github.com/larson-group/clubb/pull/1393); this branch alone cannot
produce JAX timing profiles. The Python SCM does not yet emit timing summaries.
[timing_profiles.py](timing_profiles.py) reads and writes this profile format
for the timing command and Dash.
