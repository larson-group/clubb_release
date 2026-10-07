# CLUBB Tuner Package

This package contains the Python-side orchestration for the in-memory CLUBB loss
driver and Dash tuning workflow. It is separate from `clubb_python_driver/`,
which is the Python SCM driver. The tuner package does not advance CLUBB itself;
it prepares loss-driver namelists, launches reusable loss sessions, proposes
parameter sets, and records ranked results.

## JAX tuning

Use the same tuning CLI with `-jax` (CPU) or `-jax=gpu`:

```bash
python run_scripts/run_tuner_job.py -jax -cases bomex \
  -fields cloud_frac wp2 -param_ranges C8:0.2:0.8 \
  -strategy random:32 -batch_size 8 -workers 1 -run_top never
```

This prepares the managed JAX environment and runs the existing strategies,
multi-case scheduler, baselines, checkpointing, and results through the JAX
loss driver. A Fortran build is unnecessary. Job-directory callers set
`"backend": "jax"` and optionally `"jax_options": "cpu"` or `"gpu"` in
`request.json`; omitting `backend` retains the Fortran default. Dash's Tune
runtime chooser writes these same fields into each revision. Optional
`jax_gpu` selects a full GPU UUID, and `jax_xla_prealloc` explicitly controls
CUDA preallocation for that job. Loading or continuing a revision and replaying
its results retain the saved settings. CLI callers pass one selection value,
for example `-jax=gpu,device=DEVICE,prealloc_gpu_mem=false`.
The JAX launcher validates these selections and lets the backend configure the
child environment; the tuner does not set GPU environment variables.

For a standalone loss check:

```bash
python run_scripts/run_scm_loss.py -jax -cases bomex -fields cloud_frac,wp2 \
  -multicol C8/0.2:0.8/8 -batch_size 4 -disable_stats_storage
```

For direct reusable calls, use `clubb_jax.src.clubb_loss_driver`:

```python
from clubb_jax.src import clubb_loss_driver

names, defaults = clubb_loss_driver.init_clubb_loss(runfile, return_default_params=True)
try:
    metrics = clubb_loss_driver.clubb_get_loss_for_params(defaults)
    # Pass another (candidate_count, 102) matrix using the canonical parameter order.
    next_metrics = clubb_loss_driver.clubb_get_loss_for_params(next_params)
finally:
    clubb_loss_driver.finalize_clubb_loss()
```

Each of the five returned arrays has shape `(window, field, candidate)`.
Initialization prepares benchmark interpolation and an immutable initial-state
snapshot once. Reruns reset all state and refresh parameter-derived quantities.
Post-tuning runs retain the request's other physics overrides; winning parameter
values take precedence over initial overrides for those tuned parameters.
Parameter values are dynamic JAX arrays; changing them reuses compiled model
kernels. The runtime batch width stays fixed, with a short final batch padded
and trimmed. Different cases, grid sizes, flags, stats layouts, or runtime batch
widths can compile different kernels. Initial compilation also includes the
configured sampling/radiation paths. No caches are cleared between candidates.

Variable candidate counts require in-memory stats (`stats_output_filename = ""`),
as used by the scheduler. File-backed loss runs retain their configured total
column count. Tuning from a restart is currently rejected; standalone restart
tests remain supported separately.

The existing JAX standalone feature gates apply to tuning. Source-compatible
`lmin` and maximum mixture-fraction values still use the last runtime column,
as in Fortran; varying those parameters across columns retains that source
behavior. GPU precision follows the JAX launcher's normal settings. Adam uses
the existing SPSA strategy; this workflow does not require long-run gradients.
Port outline and adaptations are in [`TUNER_PORT.md`](../clubb_jax/TUNER_PORT.md).

## Entry Points

- `python run_scripts/run_tuner_job.py ...` is the friendly command-line wrapper.
  It builds a Dash-style tuning request from compact flags, launches
  `tuner.tune_clubb`, prints live status, handles Ctrl-C graceful stop, and can
  optionally run the top result afterward.  Example:
  `python run_scripts/run_tuner_job.py -cases bomex -fields cloud_frac -param_ranges C8:0.2:0.8 -strategy random:8`.
  Case specs may include Dash-style timing as `case:t_start:t_end:t_interval`.
- `python -m tuner.tune_clubb -job_dir <dir>` runs one tuning job from a job
  directory. This is the entry point used by the Dash tuning tab.

Ad-hoc loss checks from the Dash result table go through
`run_scripts/run_scm_loss.py`. The reusable Fortran loss API and its Python
front end belong to [clubb_python_api/](../clubb_python_api/README.md#loss-evaluation);
JAX loss implementation and its front end belong to
[clubb_jax/](../clubb_jax/README.md#standalone-loss-evaluation).

## High-Level Flow

The normal Dash tuning path is:

1. `dash_app/tune_tab/runtime.py` creates a unique job directory under
   `output/tuner/`.
2. Dash writes `request.json`, `control.json`, and an initial `status.json`.
3. Dash starts `python -m tuner.tune_clubb -job_dir <dir>` as a subprocess and
   logs stdout/stderr to `worker.log`.
4. `tuner.tune_clubb` validates the request with `tuner.request.load_request`.
5. `tuner.tuning_scheduler.run_scheduler` starts one worker process per case for
   initialization. Each initial worker evaluates the CLUBB-default baseline and,
   when an override exists, the config-plus-override default baseline before
   opening its reusable candidate session.
6. Each worker uses `utilities/create_case_namelist.py` to build the
   case-specific loss namelist and normalized LES benchmark file in its worker
   directory.
7. Workers initialize the selected reusable loss driver (Fortran through
   `clubb_python`, or the independent JAX port),
   evaluate parameter batches, and return loss matrices to the scheduler.
8. The scheduler converts the per-field Taylor diagnostics into a smart loss,
   ranks samples by that value, keeps the best result rows, and periodically
   updates `status.json` and `results.json`.

## Job Directory Contract

A tuning job communicates through files in one job directory:

- `request.json`: immutable input request from Dash or another launcher.
- `control.json`: mutable control file. Currently it contains
  `{"stop_requested": true|false}` for graceful stopping.
- `status.json`: lightweight live status for polling. It contains job state,
  sample counts, elapsed time, worker counts, and a short top-results summary
  ranked by smart loss.
- `results.json`: retained result data, including selected parameters, full
  parameter rows, smart losses, `scaled_rmse_sum`, field diagnostics, baseline
  diagnostics, and per-result improvement percentages.
- `worker.log`: subprocess log for the top-level tuning process.
- `workers/<case>_<id>/`: per-worker files, including generated aggregate
  namelists, duplicated multicol parameter files, and normalized benchmark
  NetCDF files.

All JSON writes use atomic replacement so Dash does not read partially-written
files.

Dash and `run_tuner_job.py` create jobs with a controller keepalive lease in
`control.json`. The controller renews that heartbeat while polling; if it stops
renewing for 300 seconds, the scheduler treats the expired lease like a graceful
stop request. Direct `python -m tuner.tune_clubb -job_dir ...` runs are not
leased unless their `control.json` explicitly enables keepalive.

## Request Shape

`request.json` is normalized by `tuner.request`. The important fields are:

- `case_configs`: list of per-case comparison configs. Each entry contains
  `case_name`, `altitude_comparison_range`, `time_average_range`, and
  `num_time_windows`. A value of `num_time_windows = 1` is the old
  single-average behavior; larger values split the case time range into equal
  windows. An explicit window count replaces an inherited
  `average_time_seconds`; supplying an averaging interval derives the count.
- `cases`: legacy list of case names. A single legacy `case_name` is also
  accepted and normalized into `case_configs`.
- `selected_fields`: CLUBB-facing field names to compare.
- `parameter_ranges`: logical sampling coordinates with `name`, `min`, `max`,
  and optional physical `targets`.  Omitting `targets` means `[name]`; linked
  targets receive the same sampled value and must be unique across the request.
- `preset`: optional provenance name from `tuner/presets.json`.  Saved requests
  retain it alongside their fully expanded cases, fields, ranges, and override.
- `batch_size`: number of parameter columns evaluated per loss-driver call.
- `max_workers`: maximum concurrent case evaluations.
- `strategy`: tuning algorithm config.
- `case_weights` and `field_weights`: optional non-negative loss weights.
- `case_overrides`: legacy optional per-case overrides for
  `altitude_comparison_range`, `time_average_range`, and `num_time_windows`.
- `seed`: optional random seed for random, SimAnn, and Adam initialization.

Case defaults are read from `tuner/case_defaults.json`. LES benchmark files are
owned by that file only and are not request-overridable. Fields are selected
separately by the request and must be normalized CLUBB-facing names supported by
the benchmark converter.

## Presets and linked command-line ranges

`run_scripts/run_tuner_job.py -list_presets` lists the checked-in experiment
presets.  A preset supplies its normal cases, fields, parameter coordinates, and required
override; explicitly supplied `-cases`, `-fields`, or `-param_ranges` replace that
piece.  For example:

```text
python run_scripts/run_tuner_job.py -preset wpxp -strategy random:2000
```

Use `PARAM:MIN:MAX` for an ordinary range and
`PARAM=PARAM:MIN:MAX` for an equality-constrained linked range, for example
`-param_ranges C6rt=C6thl:0:4`.  The sampler treats it as one coordinate while the
saved result and generated top-result namelist retain both physical names.
Dash exposes the same request shape as either an ordinary row or a visibly
bracketed locked group with one shared range.

## Benchmark Normalization

Workers do not pass raw SAM or COAMPS variable names to the loss driver. Before a
worker initializes its loss session, the shared loss-namelist builder calls
`utilities.benchmark_converter.convert_benchmark_file` to create a normalized
NetCDF file in the worker directory. The generated `&tuner_loss_nl` then uses the
selected CLUBB field names for both `clubb_var_names` and `benchmark_var_name`.

The loss driver can therefore read the converted benchmark file as if it used
CLUBB naming conventions. Time averaging and altitude comparisons remain in the
selected loss driver (Fortran or JAX).

## Scheduler And Workers

`tuner.tuning_scheduler.TuningScheduler` owns the master process state:

- Starts and monitors worker processes using Python `multiprocessing` with the
  `spawn` start method.
- Builds a strategy object from `tuner.tuning_strategy`.
- Maintains pending samples, packed multicol batches, queued case jobs, active
  worker assignments, completed samples, and ranked best results.
- Handles graceful stop requests by stopping workers and checkpointing strategy,
  random, pending-sample, and incomplete-batch state for continuation.

The ranking loss is selected by the request's versioned Python loss policy.  The
default policy uses `loss_mode = shape_first` and
`aggregation_mode = quantile_weighted`.  Each active time-window loss is sorted
best-to-worst, divided into four equally populated bins, and the bin means are
combined with normalized best-to-worst weights `0.1, 0.4, 0.4, 0.1`.
`time_window_aggregation_scope = overall` pools all active case/field/window
losses before this calculation; `by_case` applies it within each case and then
takes the case-weighted mean.  Requests retain both the supplied weights and
scope.  Legacy `mean_max` and `mean_worst_quantile` requests remain supported
for reproducibility.

Other selectable loss modes use the explicit diagnostics returned by the
Fortran loss driver: `scaled_rmse`, correlation, standard-deviation ratio,
centered RMSE, and bias. Results retain every per-field mode score, component
contribution, and sanitization flag. `loss` and `smart_loss` are aliases for the
selected Python loss mode. `scaled_rmse_sum` is kept for comparison and
debugging, but it only ranks tuner samples when `loss_mode = scaled_rmse`.

Baselines use temporary one-column sessions and pass through this same
aggregation path. `improvement_percent = 100 * (baseline_loss - candidate_loss)
/ baseline_loss`; positive values are improvements, while nonpositive or
nonfinite baselines leave the score unavailable. Baselines do not count as
candidate samples.

`tuner.tuning_worker.worker_main` owns one initialized loss session for one case.
It receives `evaluate_batch` messages containing a full parameter matrix, calls
`clubb_api.clubb_get_loss_for_params`, and sends explicit loss-metric arrays
back to the scheduler.

## Tuning Strategies

Strategies live in `tuner.tuning_strategy` and share a small interface:

- `fill(pending_samples, capacity)` proposes samples until the pending queue is
  full or the strategy is exhausted.
- `tell(completed_samples)` receives completed samples and advances adaptive
  strategies.
- `is_exhausted()` reports whether no more samples are available.
- `estimated_sample_count()` returns a finite count when known.

Current strategies:

- `random`: uniform random samples inside each selected parameter range, with
  optional `max_samples`.
- `resolve`: deterministic full-grid sampling using a requested spacing.
- `simann`: independent enhanced simulated-annealing chains.
- `adam`: projected Adam in normalized coordinate space using averaged SPSA
  gradient pairs. Dash presents learning and perturbation radii as percentages
  of each configured range. CLI syntax is
  `adam:MAX_UPDATES:LEARNING_RATE:PERTURBATION:SPSA_PAIRS`, using normalized
  fractions. Request normalization derives chain and concurrent-batch counts.

## Module Layout

- `tune_clubb.py`: thin CLI around request loading, scheduler launch, and
  top-level error handling.
- `job_runtime.py`: shared `TunerJob` wrapper for creating job directories,
  launching `tune_clubb`, reading status/results, and requesting graceful stop.
- `request.py`: request validation and parsing of case default files.
- `status.py`: atomic JSON I/O, status/results writers, and stop-control helpers.
- `tuning_scheduler.py`: master process scheduling and result aggregation.
- `tuning_worker.py`: child worker loop for one case-specific reusable loss
  session.
- `tuning_strategy.py`: parameter proposal algorithms.
- `adam_spsa_strategy.py`: multi-chain Adam with SPSA gradients and seeded
  Latin-hypercube starts.
- `paths.py`: shared repository paths used by tuner modules.
