# Tune CLUBB parameters

The tuner runs candidate parameter sets in selected cloud cases, compares model
profiles with LES reference profiles, and retains the candidates with the lowest
loss. Use the CLI below or [Dash’s Tune tab](../dash_app/tune_tab/tab.py) to select
cases, parameters, a search algorithm and the Fortran or JAX backend.
See [algorithm choices](algorithms/README.md) for their differences and requirements.

## Quick start

Run these commands from the repository root. This evaluates 16 Latin hypercube
points for C8 using CPU JAX, two model columns at a time:

```bash
python run_scripts/run_tuner_job.py -jax=cpu -cases bomex -fields wp2 \
  -param_ranges C8:0.2:0.8 -strategy latin_hypercube:num_points=16 \
  -batch_size 2 -workers 1 -seed 7 -run_top never
```

The [tuner CLI](../run_scripts/run_tuner_job.py) prepares its managed Python/JAX
environment. For the Fortran backend, first build its Python loss API, then omit
`-jax`:

```bash
python compile.py -debug -python
python run_scripts/run_tuner_job.py -cases bomex -fields wp2 \
  -param_ranges C8:0.2:0.8 -strategy random:max_samples=16 \
  -batch_size 2 -workers 1 -run_top never
```

Both commands print progress and paths to the saved request, results and logs.
`-run_top never` ends after tuning; `complete`, `window` or `both` rerun selected
best candidates for model output. `-top_n` controls how many to rerun, and
`-output_run_dir` chooses their output directory. The default `-run_top ask`
prompts after tuning. Do not interpret a short example as calibrated convergence
or a scientifically validated parameter choice.

## Choose cases, fields and parameter ranges

| Setting | Meaning and example |
| --- | --- |
| `-cases` | One or more `CASE`, `CASE:START:END` or `CASE:START:END:INTERVAL` specifications. Times are absolute model seconds. For example, `bomex:10800:21600:3600` scores three consecutive one-hour averaging windows after spinup. |
| `-fields` | Model field names, separated by spaces or commas; the same selection applies to every case. Each case can have its own height/time comparison ranges through its saved case config. |
| `-param_ranges` | `PARAM:MIN:MAX` in physical units. Use `C6rt=C6thl:0:4` to drive equal physical parameters with one search coordinate. Untuned parameters keep configured defaults. Fixed ranges use equal endpoints. |
| `-config` | Shared input preset name, default `default`. |
| `-override` | CLUBB parameter/flag overrides such as `C11=0.41`. Required parameter links and inactive parameters are validated against these settings. |
| `-preset` | Named experiment from [presets.json](presets.json). Explicit cases, fields, ranges and config replace the corresponding preset pieces. |
| `-strategy` | Algorithm name and colon-separated options, for example `adam_gradient:max_updates=30:learning_rate=0.01`. Compact positional spelling such as `random:200` is also supported. |
| `-seed` | Reproducible randomized starts/proposals. |

Discover presets, algorithms and the complete option reference without launching
model work:

```bash
python run_scripts/run_tuner_job.py -list_presets
python run_scripts/run_tuner_job.py -list_algorithms
python run_scripts/run_tuner_job.py -help
```

[case_defaults.json](case_defaults.json) supplies each case’s LES file, height
range and default time settings. A raw request uses a `cases` list or
`case_configs` objects with `case_name`, `altitude_comparison_range`,
`time_average_range`, and either `average_time_seconds` or `num_time_windows`.
Window boundaries must align with model sampling and radiation timesteps.
The model advances through spinup before scoring; it compares time-averaged
profiles, rather than averaging instantaneous RMSE values.

## Run gradient Adam

`adam_gradient` automatically requests loss gradients and requires JAX. For
example, run two independent chains over two cases and two fields:

```bash
python run_scripts/run_tuner_job.py -jax=cpu \
  -cases bomex:0:600:300 atex:0:600:300 -fields wp2 thlm \
  -param_ranges C8:0.2:0.8 C11:0.1:1.0 -override C11=0.41 \
  -strategy adam_gradient:max_updates=2:learning_rate=0.01 \
  -gradient_mode forward -batch_size 2 -workers 1 -run_top never
```

Here each case runs ten 60-second steps and supplies two averaging windows.
The first chain starts at clipped defaults, with C11 set to 0.41; additional
chains use Latin hypercube starts. Gradient Adam’s chain count is always
`workers × batch_size`. Both Adam methods default to 30 updates, but SPSA uses
a different probe layout and evaluation budget; see the algorithm guide.

C11=0.41 avoids the existing `C11 == C11b` shortcut at the default equality,
whose executed branch has zero derivative with respect to C11. Automatic
differentiation follows executed model branches. A finite or zero derivative
therefore needs interpretation before treating it as a useful search direction.

`-gradient_mode forward` differentiates selected parameter directions; `reverse`
uses a backward pass with timestep recomputation. `auto` chooses reverse when
more than one physical parameter targets the case loss, otherwise forward.
This heuristic does not predict runtime or memory. Reverse mode can use much
more memory, especially with multiple cases/chains. Start with small batches
and worker counts. `-jax=gpu` selects the managed GPU runtime; GPU-gradient
convergence/performance is not established by the CPU tests.

Value-only methods work with Fortran or JAX. An explicit `-grad` also records
JAX gradients for those methods without changing how they choose candidates.
Dash enables gradients from the algorithm requirement; Fortran with a gradient
method shows an incompatibility message and disables Start.

## Control the loss being minimized

![Candidate proposals, case evaluation, score aggregation and feedback.](diagrams/tuning_loop.svg)

Every candidate runs in every selected case. The [loss driver](../utilities/clubb_loss_runner.py)
compares each field/window’s mean vertical profile with the LES reference,
then returns one scalar loss per case. The tuner combines those case losses
into the final `total_loss` used for ranking and algorithm feedback.

`-loss_weights` changes the positive mean-offset, standard-deviation and
correlation error contributions, defaulting to one each. With normalized bias
b, model/reference standard-deviation ratio s and correlation r:

```math
\mathrm{rmse}^2 = w_{\mathrm{mean}}b^2 + w_{\mathrm{std}}(s-1)^2
                 + w_{\mathrm{corr}}2s(1-r).
```

These moments are across compared height levels after time averaging. Weights
do not alter the raw bias, ratio or correlation diagnostics. With unit weights,
the result is profile RMSE divided by reference standard deviation. Flat
profiles use the loss kernel’s defined moment-based objective; undefined
diagnostics remain NaN (JSON `null`) rather than invented correlation values.
The [shared kernel](../utilities/loss_metrics.py) defines those edge cases.

Fields combine by weighted RMS, with equal weight per window inside a case.
Raw requests can specify nonnegative `field_weights`; zero excludes a field.
Cases receive one scalar each regardless of their window count. Nonnegative
`case_weights` in a raw request control case contributions; zero excludes a case.
At least one field/case contribution must remain positive.

| `-aggregation_mode` | Across-case score |
| --- | --- |
| `quantile_weighted` (default) | Sort cases by loss and combine up to four case-weighted bin means; `-aggregation_weights` gives best-to-worst weights, default `0.1 0.4 0.4 0.1`. Empty bins are omitted and weights renormalized. |
| `mean_max` | 60% case-weighted mean and 40% worst case loss. |
| `mean_worst_quantile` | 70% case-weighted mean and 30% case-weighted mean of the worst quarter (at least one case). |

For example, change the Taylor weights and across-case policy with:

```bash
python run_scripts/run_tuner_job.py -jax=cpu \
  -cases bomex:10800:21600:3600 arm:68400:72000 -fields wp2 thlm \
  -param_ranges C8:0.2:0.8 -strategy random:8 \
  -loss_weights '{"mean":2,"std":0.5,"corr":1}' -aggregation_mode mean_max \
  -batch_size 2 -workers 1 -run_top never
```

JAX computes each case’s scalar-loss derivative through the model trajectory,
spinup, statistics and profile comparison. [Gradient helpers](gradients.py)
combine those case derivatives using the actual ranking policy, then sum linked
partials and scale parameter ranges for Adam. Sorting/ties can make the objective
nonsmooth; rank changes with unequal weights can be discontinuous. No smoothing
is applied. Nonfinite contributing derivatives stop gradient tuning.
The [case-loss diagram](diagrams/case_loss.svg) shows the differentiated path.

## Concurrency, outputs and continuation

`-batch_size` is the number of candidate model columns per evaluation (default 8).
`-workers` limits concurrent candidate-case evaluations (default half the available
logical CPUs). Cases have separate warm processes; initialization can keep more
resident processes than the active-evaluation limit. Worker count, batch width,
case count and JAX differentiation all affect memory consumption.

Jobs normally live under `output/tuner/`. `-output_root` changes that root;
`-output_job_dir` selects an exact new job directory. `-dry_run` writes the prepared
job without running it. Inspect the saved request before launching expensive work:

```bash
python run_scripts/run_tuner_job.py -jax=cpu -cases bomex -fields wp2 \
  -param_ranges C8:0.2:0.8 -strategy latin_hypercube:200 \
  -batch_size 2 -workers 1 -output_job_dir output/tuner/my_prepared_job -dry_run
.venv-python/bin/python -m tuner.tune_clubb -job_dir output/tuner/my_prepared_job
```

The [job-directory entry point](tune_clubb.py) also supports `-help`. It runs an
existing request using its recorded backend; the high-level CLI owns creation
and post-tuning reruns.

| File/directory | Contents |
| --- | --- |
| `request.json` | Experiment settings, algorithm options, cases and ranges. |
| `status.json` | Live progress, worker counts and a bounded best-results summary. |
| `results.json` | Up to 16 best candidates overall/per case, diagnostics, baselines and optional gradients. |
| `sample_history/` | All completed candidates in compressed NPZ chunks, including full parameters, case losses, window/field metrics and optional gradients. |
| `control.json` | Graceful-stop request and controller keepalive. |
| `resume_checkpoint.pkl` | Saved algorithm/pending-work state after a graceful stop, used by Continue. |
| `worker.log`, `workers/` | Controller and per-worker setup/log files. |

Histories label physical gradient columns with `gradient_parameters`. Per-case
loss gradients have candidate × physical-parameter axes; the combined gradient
has the same order without a case axis. Profile metrics retain case/window/field
labels. Public JSON uses `null` for numeric NaNs. The first completed history set
publishes immediately; pending rows then publish every 500 samples or 20 seconds,
whichever comes first. Stopping flushes remaining rows. Dash’s landscape updates
from these published chunks. Its per-field summary uses window-RMS error and
other diagnostics from the worst-RMSE window; full windows stay in history.

Dash workspace Continue restores a gracefully stopped search’s saved state.
This is distinct from JAX timestep checkpoints, which control differentiation
memory. Controllers renew a keepalive lease; a 300-second lapse requests a
graceful stop. Direct job-directory runs are not leased unless enabled in their
control file. Current requests/checkpoints/history paths are the supported format:
**old strategy import paths, scalar algorithm options, global window settings,
job-root history chunks and `output_tuner/` discovery are no longer supported.**
No backward-compatibility adapters are provided.

## Use the loss runner without tuning

The [standalone loss CLI](../run_scripts/run_scm_loss.py) evaluates the same case
loss without an optimizer:

```bash
python run_scripts/run_scm_loss.py -jax=cpu -cases bomex -fields wp2
python run_scripts/run_scm_loss.py -jax=cpu -cases bomex -fields wp2,thlm \
  -multicol C8/0.2:0.8/8 -batch_size 4 -disable_stats_storage
```

Programmatic callers initialize `clubb_loss_runner`, repeatedly call
`evaluate(params, grad=..., gradient_mode=..., gradient_target=...)`, then close
it. Each parameter row follows `runner.param_names`. The normal loss target
returns one loss per candidate and one gradient per candidate/selected physical
parameter. Explicit `gradient_target="metrics"` returns separate field/window
metric derivatives for diagnostics and costs more backward work.

## Adding an algorithm

The [JSON catalog](algorithms/algorithms.json) lists available algorithms, their
categories, options and requirements. The [registry loader](algorithm_registry.py)
uses it for request validation, strategy construction, CLI parsing and Dash
controls. Inspect the catalog without launching a job:

```bash
python run_scripts/run_tuner_job.py -list_algorithms
```

To add a method:

1. Add its Python file under [algorithms/](algorithms/), subclassing
   [BaseTuningStrategy](algorithms/base.py). Its constructor accepts
   `param_names`, `default_params_row`, `parameter_ranges`, `seed` and its
   declared options. Implement `fill()`, `tell()` and `is_exhausted()`; use
   `_sample_from_selected_params()` to expand linked coordinates and assign
   sample IDs. Feedback can arrive partially or out of order. Keep adaptive
   context by sample ID and make its state pickleable for continuation.
2. Add its entry to `algorithms.json`. This minimal example describes the
   existing gradient Adam class; the real entry also supplies labels, bounds
   and its additional optional controls:

   ```json
   {
     "adam_gradient": {
       "implementation": "tuner.algorithms.adam_gradient:AdamGradientStrategy",
       "category": "optimizer",
       "requires_gradients": true,
       "options": {
         "learning_rate": {"type": "float", "default": 0.01},
         "max_updates": {"type": "int", "default": 30}
       },
       "internal_options": ["chain_count"]
     }
   }
   ```

Every algorithm entry requires `category`: use `optimizer` when losses or
gradients guide subsequent proposals, or `sampler` for loss-independent point
selection. The curated-function benchmark runs every method and uses this
category to require loss improvement only from optimizers.

An option's type is `int`, `float`, `string` or `bool`. Optional fields include
`label`, `description`, `minimum`, `maximum`, `exclusive_minimum`,
`exclusive_maximum`, `choices`, `required` and `ui_default`. Labels default to
the option name with spaces. `default: null` leaves an option unspecified.
Values use the request's units; learning rate `0.01` is a 1% normalized step.

Execution logic stays in Python. Override `prepare_options(options, request)`
to derive scheduling values, and `sample_count(options, ranges)` to report a
budget before construction. Their inherited defaults need no extra hookup.
Declare any derived, non-public option names in the JSON `internal_options`
list, as [Adam with SPSA](algorithms/adam_spsa.py) does.

CLI and Dash discover catalog entries on their next start. Dash generates
numeric/text entries and boolean/choice menus from the JSON types; no extra
frontend registration or callbacks are needed. Entries sit beside their labels;
numeric entries use manual input without spinner arrows. The gradient, RMSE
weight and across-case aggregation sections have question-mark help notecards.
Saved requests use the same
`strategy.options` mapping, including the typed agent API's `strategy_options`.
Requests use a strategy object (`{"name": "random", "options": {"max_samples": 200}}`)
and a `cases` list or per-case `case_configs`. Window settings belong inside each
case config. Checkpoints and sample history use the current module paths and
`sample_history/` directory; older formats are not migrated.

Use `-strategy method_name:option=value:other=value` for any registered method.
`cli_order` optionally defines a compact positional spelling.
The scheduler returns derivatives of the executed numerical program; optimizers
should not assume that every active parameter has a nonzero derivative.
`tell()` receives `sample_id`, `selected_params`, `all_params`, `total_loss`,
`case_loss` and diagnostics. With gradients enabled it also receives
`gradient_parameters`, `gradient[case]["loss"]` and `total_gradient`.
Gradients are physical partial derivatives; use `range_gradient()` when an
optimizer works in linked or normalized coordinates. Algorithms may cap
in-flight batches through `pending_batch_limit(batch_size)`; workers and model
batch width remain scheduler concerns.

## Validation

The [algorithm guide](algorithms/README.md#assessing-the-methods) describes the
fast objective benchmarks and their limits. The Jenkins tuner job runs all tuner
unit tests first, before compilation, and publishes their JUnit results. It also
runs [the gradient integration check](auto_llm_generated_tests/integration_tests/run_adam_gradient_test.py):

```bash
python tuner/auto_llm_generated_tests/integration_tests/run_adam_gradient_test.py \
  -output_dir output/tests/adam_gradient -gradient_mode forward -batch_size 2
python tuner/auto_llm_generated_tests/integration_tests/run_adam_gradient_test.py \
  -output_dir output/tests/adam_reverse -gradient_mode reverse
```

Each output directory must be new. The bounded BOMEX/ATEX check verifies first
updates and transported gradients against independent finite differences, and
requires the retained best loss to be no worse than the first chain’s starting
loss. Jenkins runs CPU forward mode with two chains to keep memory bounded;
reverse remains an explicit manual option. These short runs establish mechanics,
not long-duration scientific convergence or GPU performance.

Diagram sources are [tuning_loop.dot](diagrams/tuning_loop.dot) and
[case_loss.dot](diagrams/case_loss.dot). Render them with Graphviz:

```bash
dot -Tsvg tuner/diagrams/tuning_loop.dot -o tuner/diagrams/tuning_loop.svg
dot -Tsvg tuner/diagrams/case_loss.dot -o tuner/diagrams/case_loss.svg
```
