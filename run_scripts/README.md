# run_scripts

This directory contains scripts for running CLUBB, generating output, plotting
output, and comparing output. Examples below assume they are run from the repo
root unless the example explicitly changes directories.

Most model-running scripts expect CLUBB to have already been compiled with
`./compile.py`, unless the script has its own compile step.

## Top-Level Scripts

### `run_scm.py`

Runs individual CLUBB cases, inherited case groups, or named run lists with
native Fortran, the Python SCM or JAX. Selections come from
[`case_definitions.jsonc`](../input/case_definitions.jsonc); flag presets and
collections come from [`clubb_config_flags.jsonc`](../input/clubb_config_flags.jsonc).
It prepares each aggregate namelist in memory, then submits it to the shared
[`clubb_case_runner`](../utilities/clubb_case_runner.py). A single resolved run
uses no worker pool or preparation subprocess. Multiple runs use at most
`-workers` concurrent tasks (default: half the logical CPUs); `-workers 1`
executes them sequentially. The model itself runs in its normal process.

Examples from the repository root:

```sh
python3 run_scripts/run_scm.py bomex
python3 run_scripts/run_scm.py -cases priority_cases -workers 4
python3 run_scripts/run_scm.py -cases clex_group bomex -max_iters 5
python3 run_scripts/run_scm.py -cases arm,bomex,atex -workers 1 -show_output
python3 run_scripts/run_scm.py -jax=cpu -output_dir jax_bomex -max_iters 20 bomex
python3 run_scripts/run_scm.py -config compatible_r8029 arm
python3 run_scripts/run_scm.py -flags run_bindiff_w_flags_config_example_group bomex
```

`standard_cases`, `priority_cases`, `min_cases`, `short_cases` and `all_cases`
retain the former batch runner's exact membership and order. `all_cases` includes
unmaintained cases. Case groups use names ending in `_group`, such as
`clex_group`. Lists select cases without changing their settings. Repeated case
selections run once in first-use order. Unknown selectors and list cycles fail
before model execution. Supply an explicit case selection; no arguments do not
start the standard suite.

Generated inputs, console logs and configured stats go under `output/` or
`-output_dir`. Flag-collection sweeps use a separate subdirectory per qualified
member; they run exactly the selected members. Batch output is captured in each
case's log; `-show_output` also echoes it. Independent runs continue after a
failure, and any failed run makes the command return nonzero.

`-override` accepts assignments, inline JSON or a JSON filename. JSON may map
settings directly, or contain case-specific objects plus an optional `all`
object. Shared values apply first, then case-specific values. Repeated
`-override` arguments apply in order, with later values taking precedence:

```sh
python3 run_scripts/run_scm.py -cases arm,bomex -override '{"all":{"dt_main":60},"bomex":{"dt_main":30}}'
python3 run_scripts/run_scm.py bomex -override C8=0.4 -override ~/experiment.json
python3 run_scripts/run_scm.py arm -stats core+radiation -debug 0
python3 run_scripts/run_scm.py arm -zm_grid input/grid/31_level_zm_grid.grd -nzmax 31
python3 run_scripts/run_scm.py arm -gdb
```

JAX device/environment options remain owned by its launcher. Use `-python` for
the Python SCM with a Python-enabled build. Timing workflows use the same
runner through [`time_clubb.py`](../utilities/time_clubb.py).

### `run_clubb_w_varying_flags.py`

Runs a selected case list under each independent member of a flag collection.
`-flag_group` selects a group, list or preset from the flag catalog. An
unmodified `default` run is prepended unless `-skip_default_flags` is supplied.
Local names such as `flag1` remain unchanged within each collection, including
their output subdirectory names. When a selection combines collections with
colliding member names, qualified labels distinguish their outputs.

```sh
python3 run_scripts/run_clubb_w_varying_flags.py bomex
python3 run_scripts/run_clubb_w_varying_flags.py -cases priority_cases -workers 4
python3 run_scripts/run_clubb_w_varying_flags.py -flag_group run_bindiff_w_flags_config_host_flags_group -cases min_cases
python3 run_scripts/run_clubb_w_varying_flags.py -flag_group run_bindiff_w_flags_config_example_group -skip_default_flags -max_iters 10 bomex
```

The default collection is `run_bindiff_w_flags_config_core_flags_group`. Sweep
settings apply after the selected baseline, followed by user `-override`
arguments. Members never accumulate changes from previous runs. For comparing
checkouts, `-clubb_root PATH` uses the current catalog definitions while running
that checkout's SCM script and compiled model.

### `run_bindiff_all.py`

Compares matching CLUBB NetCDF files in two output directories. Variables
found on only one side and variable shape mismatches are warnings by default;
`-strict` makes them failures too. Comparisons fail in either mode if a file pair
has no comparable numeric stats fields (matching coordinates alone are not
enough), or if the files have different numbers of timesteps.

Examples:

- `python3 run_scripts/run_bindiff_all.py output_a output_b`
  Compares every matching case in two output directories with default reporting.

- `python3 run_scripts/run_bindiff_all.py -verbose 2 output_a output_b`
  Prints detailed per-variable numerical differences.

- `python3 run_scripts/run_bindiff_all.py -strict output_a output_b`
  Also fails unless every variable has a counterpart with the same shape in
  the paired file. Normally these are warnings to prevent metadata changes
  in the netcdf or added/removed variables from counting as differences. If
  every variable has a shape/mismatch issue, the test will fail with or without
  the -strict flag.

- `python3 run_scripts/run_bindiff_all.py -case bomex -threshold 1e-12 -percent_threshold 1e-12 output_a output_b`
  Compares only `bomex` with strict absolute and percent thresholds.

- `python3 run_scripts/run_bindiff_all.py -log_mode replace -verbose 2 output_a output_b`
  Writes per-case diff logs under `output/bindiffs`, replacing older logs.

- `python3 run_scripts/run_bindiff_all.py -result_json comparison.json output_a output_b`
  Writes a compact machine-readable report with input and log paths, case and
  file status, variable names grouped by exclusive comparison outcome, and
  the first failing saved record. This is mainly for usage in other scripts.

- `python3 run_scripts/run_bindiff_all.py -flag_sets old_output new_output`
  Treats each immediate child directory as a flag-set output directory and
  compares matching flag sets.

### `run_clubb_and_pyplotgen.py`

Compiles CLUBB, runs selected cases, and plots those outputs with pyplotgen.
Existing output directories for the selected run name are removed before the
run unless `-dry_run` is used.

Examples:

- `python3 run_scripts/run_clubb_and_pyplotgen.py`
  Compiles CLUBB, runs the default case list, writes to `output/new`, and writes
  plots to `output/pyplots_new`.

- `python3 run_scripts/run_clubb_and_pyplotgen.py -name smoke -cases bomex arm -max_iters 10`
  Runs and plots only `bomex` and `arm` for 10 iterations under
  `output/smoke` and `output/pyplots_smoke`.

- `python3 run_scripts/run_clubb_and_pyplotgen.py -skip_compile -name trial -compare baseline`
  Uses the existing executable, writes a new `trial` output, and overlays the
  existing `output/baseline` data in the plots.

- `python3 run_scripts/run_clubb_and_pyplotgen.py -dry_run -name check -cases bomex`
  Prints the compile, run, and plotting commands without creating output.

### `run_scm_loss.py`

Runs the in-memory CLUBB loss driver for one or more cases. It prepares the
comparison namelist and uses the shared case runner to select the native,
Python or JAX frontend. Each case runs once and prints its loss table against
the configured benchmark data. Benchmark files and default comparison windows
come from [`tuner/case_defaults.json`](../tuner/case_defaults.json). This command
finishes and releases each loss run; the tuner separately retains initialized
allocations through its loss API.

Examples:

- `python3 run_scripts/run_scm_loss.py bomex -fields cloud_frac,rcm`
  Runs the loss driver for `bomex` using the case-default comparison settings.

- `python3 run_scripts/run_scm_loss.py -cases bomex,arm -fields cloud_frac,rcm`
  Runs the loss driver for two cases with the same requested fields.

- `python3 run_scripts/run_scm_loss.py bomex -fields cloud_frac -num_time_windows 4`
  Splits the case-default time range into four loss windows.

- `python3 run_scripts/run_scm_loss.py bomex -fields cloud_frac -multicol C8/0.2:0.8/4 -batch_size 2`
  Builds a multi-column parameter sweep and evaluates two columns at a time.

- `python3 run_scripts/run_scm_loss.py -python bomex -fields cloud_frac`
  Runs the Python loss-driver front end instead of the compiled executable.

- `python3 run_scripts/run_scm_loss.py -jax=cpu bomex -fields cloud_frac`
  Runs the JAX loss frontend through its managed CPU launcher.

### `run_tuner_job.py`

Runs the newer Dash-style tuner workflow from the command line. It creates a
job directory, writes request/control/status files, starts worker processes, and
can optionally run the top parameter sets afterward.

Examples:

- `python3 run_scripts/run_tuner_job.py -cases bomex -fields cloud_frac -param_ranges C8:0.2:0.8`
  Runs a small random tuner job for one case, one field, and one parameter.

- `python3 run_scripts/run_tuner_job.py -cases arm:10800:21600:10800 bomex:7200:18000:2700 -fields cloud_frac rcm -param_ranges C8:0.2:0.8 C11:0.1:1.0 -strategy resolve:0.1`
  Tunes two cases with explicit time windows, two fields, two parameters, and a
  resolved grid strategy.

- `python3 run_scripts/run_tuner_job.py -cases bomex -fields cloud_frac -param_ranges C8:0.2:0.8 -dry_run`
  Writes the job files but does not launch the tuner.

- `python3 run_scripts/run_tuner_job.py -cases bomex -fields cloud_frac -param_ranges C8:0.2:0.8 -run_top complete -top_n 3`
  Runs the tuner, then runs standalone CLUBB output for the top three complete
  parameter sets.

- `python3 run_scripts/run_tuner_job.py -list_algorithms`
  Lists the [algorithm catalog](../tuner/algorithms/algorithms.json)'s typed
  options, defaults and backend requirements as JSON without launching a job.

- `python3 run_scripts/run_tuner_job.py -jax -cases bomex:0:600:300 -fields wp2 -param_ranges C8:0.2:0.8 -strategy adam_gradient:max_updates=2:learning_rate=0.01 -batch_size 1 -workers 1 -run_top never`
  Applies two bounded Adam updates using JAX gradients. This method enables
  gradients automatically and requires `-jax`. The existing `adam` method uses
  SPSA estimates. See the [tuner workflow](../tuner/README.md) for multiple cases,
  linked parameters and extending the algorithm registry.

## `convergence_run/` Scripts

The `convergence_run/` directory contains an older convergence-test workflow.
These scripts are more environment-specific than the top-level runners and may
need local module paths, output paths, or compiler settings updated before use.

### `convergence_run/convergence_config.py`

Builds and runs one convergence-test configuration. It modifies generated input
for a requested case, grid spacing/refinement, timestep, and convergence-study
option set.

Examples:

- `(cd run_scripts && python3 convergence_run/convergence_config.py bomex -dt 60 -dz 20 -skip_check)`
  Configures and runs a `bomex` convergence case with 60-second timesteps and a
  20 m uniform grid.

- `(cd run_scripts && python3 convergence_run/convergence_config.py rico -dt 30 -ref 2 -rad_off -micro_off -skip_check)`
  Runs `rico` with a refined grid and radiation/microphysics disabled.

- `(cd run_scripts && python3 convergence_run/convergence_config.py dycoms2_rf02_nd -dt 10 -ref 3 -new_ic -new_bc -skip_check)`
  Runs a DYCOMS-II RF02 ND convergence setup with modified initial and boundary
  conditions.

### `convergence_run/convergence_function.py`

Helper module imported by `convergence_config.py`, mainly for modifying initial
condition profiles. It is not intended to be called directly.

Example:

- `python3 run_scripts/convergence_run/convergence_config.py ...`
  Use `convergence_config.py`; it imports this helper as needed.

### `convergence_run/plot_l2_convergence.py`

Template-style plotting script for L2 convergence figures. The multi-case C
shell scripts fill in placeholders such as case name, experiment name, timestep
values, and refinement levels before plotting.

Example:

- `(cd run_scripts && python3 convergence_run/plot_l2_convergence.py)`
  Runs the plot script after its placeholders have been replaced by the
  convergence workflow.

### `convergence_run/run_cnvg_test_multi_cases_default.csh`

Runs the multi-case convergence workflow with the default convergence
configuration.

Examples:

- `(cd run_scripts/convergence_run && csh run_cnvg_test_multi_cases_default.csh)`
  Runs the default convergence workflow for all configured cases.

- `(cd run_scripts/convergence_run && csh run_cnvg_test_multi_cases_default.csh dycoms)`
  Runs only the DYCOMS-II subset in the default convergence workflow.

### `convergence_run/run_cnvg_test_multi_cases_baseline.csh`

Runs the multi-case convergence workflow with the baseline convergence
configuration.

Examples:

- `(cd run_scripts/convergence_run && csh run_cnvg_test_multi_cases_baseline.csh)`
  Runs the baseline convergence workflow for all configured cases.

- `(cd run_scripts/convergence_run && csh run_cnvg_test_multi_cases_baseline.csh dycoms)`
  Runs only the DYCOMS-II subset in the baseline convergence workflow.

### `convergence_run/run_cnvg_test_multi_cases_revall.csh`

Runs the multi-case convergence workflow with the full revised convergence
configuration.

Examples:

- `(cd run_scripts/convergence_run && csh run_cnvg_test_multi_cases_revall.csh)`
  Runs the revised-all convergence workflow for all configured cases.

- `(cd run_scripts/convergence_run && csh run_cnvg_test_multi_cases_revall.csh dycoms)`
  Runs only the DYCOMS-II subset in the revised-all convergence workflow.

### `convergence_run/run_clubb_default`

Small `at` wrapper for launching the default convergence workflow and writing
`clubb_output_default.txt`.

Example:

- `(cd run_scripts/convergence_run && at now -f run_clubb_default)`
  Submits the default convergence workflow through `at`.

### `convergence_run/run_clubb_baseline`

Small `at` wrapper for launching the baseline convergence workflow and writing
`clubb_output_baseline.txt`.

Example:

- `(cd run_scripts/convergence_run && at now -f run_clubb_baseline)`
  Submits the baseline convergence workflow through `at`.

### `convergence_run/run_clubb_revall`

Small `at` wrapper for launching the revised-all convergence workflow and
writing `clubb_output_revall.txt`.

Example:

- `(cd run_scripts/convergence_run && at now -f run_clubb_revall)`
  Submits the revised-all convergence workflow through `at`.

## Shared command-line conventions

Use single-dash options with underscores and `-h` or `-help` for help.
CPU worker pools use `-workers`, defaulting to half the logical CPUs available
through CPU affinity (at least one). GPU comparison runs remain serial.
`-threads` controls OpenMP threads; `-process_counts` selects a benchmark sweep.
Pyplotgen chooses its own pool size; `-serial` disables that pool.

Use `-params_file` for model parameters, `-silhs_params_file` for SILHS parameters,
and `-param_ranges` for tuner search ranges. `-batch_size` controls model column
batching, including parameter samples represented by columns in tuner workers.
`-multicol NUM` selects count-based generation; `-multicol C8/0.2:0.8/4` selects
an explicit grid. The standalone generator also preserves advanced `-mode`
settings; `-mode hypergrid` interprets NUM as points per parameter dimension.

Output locations use `-output_dir`, `-output_root`, or `-output_file`.
Tuner destinations are `-output_job_dir` and `-output_run_dir`; its internal
scheduler's `-job_dir` selects an existing job. Relative-path resolution follows
each script's documented behavior. Case lists accept commas and spaces.
Model options pass through wrappers to the SCM runner; namelist definitions and
JSON `-override` resolution are shared in [input/scripts/create_case_namelist.py](../input/scripts/create_case_namelist.py).
The [input guide](../input/README.md) describes defaults, shared `-config` settings
and explicit case/parameter/flag/stats presets or files.

JAX runtime settings travel in the existing opaque `-jax=VALUE` selection,
for example `-jax=gpu,device=DEVICE,prealloc_gpu_mem=false`.
This works with `run_scm.py`, `run_scm_loss.py`, and `run_tuner_job.py`.
Add `profile` to record a full execution profile, for example `-jax=profile`
or `-jax=gpu,profile`; see [JAX profiling](../clubb_jax/README.md#timing-and-experimental-execution-traces)
for capture and viewer commands. SCM scripts forward the value unchanged; the
JAX launcher owns its syntax and validation. `prealloc_gpu_mem` controls CUDA
GPU memory preallocation; both `true` and `false` override the inherited setting.
Use `clubb_jax/run_jax.py -profile=gpu -info=json` to inspect available devices
and runtime readiness. Omitting these settings preserves the inherited device
selection and launcher memory policy; the JAX backend owns their interpretation.
