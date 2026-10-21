# CLUBB JAX Driver

CLUBB JAX is a JAX port of the Fortran CLUBB single-column model (SCM). Its
Python modules follow the Fortran source layout and translate the supported
model initialization, timestep loop, CLUBB core, radiation, and statistics
code into JAX operations.

The supported standalone path is JAX-owned and runs without a compiled Fortran
library. The Fortran standalone is needed only when comparing JAX results
against the original implementation.

The main entry point is [`src/clubb_standalone.py`](./src/clubb_standalone.py),
and translated source lives under [`src/`](./src/):

- `src/CLUBB_core/` contains the supported CLUBB core and JAX statistics code.
- `src/Input_fields/` contains namelist, sounding, surface, and grid readers.
- `src/Radiation/` contains the supported radiation code.
- `src/Microphys/` contains KK and Morrison microphysics, including cloud
  sedimentation and feedback to CLUBB. The strict comparison suite covers KK
  and short Morrison runs; atmospheric ice and graupel validation remains
  incomplete.
- `src/SILHS/` contains native JAX Latin-hypercube sampling and its KK/Morrison
  microphysics hookup. Interactive and non-interactive modes are enabled;
  random-stream matching with Fortran is deferred. See
  [`SILHS_PORT_NOTES.md`](./SILHS_PORT_NOTES.md) for configuration and validation.
- `src/Radiation/BUGSrad/` remains disconnected from the supported standalone path.

The detailed support boundary and conversion workflow are documented in
[`JAX_CONVERSION_PLAN.md`](./JAX_CONVERSION_PLAN.md).

## Running

CLUBB JAX uses the normal SCM runner and case names. A JAX-only run does not
require a Fortran build or any manual dependency setup; the command prepares
what it needs automatically and writes results under `output/`.

### Basic Run

Run a case from the repository root by adding `-jax`:

```bash
./run_scripts/run_scm.py -jax arm
```

Bare `-jax` selects CPU unless the legacy `CLUBB_JAX_ACCELERATOR` environment
variable is set. The explicit forms `-jax=cpu` and `-jax=gpu` are also
accepted, case-insensitively.

### Run Options

The statistics registry controls which model fields are collected and written.
Statistics are particularly expensive in JAX, so selecting only the output a
run needs can substantially reduce its runtime. The default is
`input/stats/standard_stats.in`, which provides broad output coverage.

```bash
# Faster: retain a small set of important multi-column fields.
./run_scripts/run_scm.py -jax -stats input/stats/multi_col_stats.in arm

# Fastest: run without collecting or writing statistics.
./run_scripts/run_scm.py -jax -stats none arm

# Most comprehensive: collect every registered statistic.
./run_scripts/run_scm.py -jax -stats input/stats/all_stats.in arm
```

`multi_col_stats.in` retains 15 central state and turbulence fields and is
substantially faster than the default. `none` has the lowest overhead but
produces no statistics output, so it is useful for timing rather than validating
model results. `all_stats.in` gives maximal output coverage at the highest cost.

Debug checks also affect performance. Many cases default to `debug_level = 2`;
passing `-debug 0` typically saves about 10% of runtime and is appropriate when
the additional checks are not needed. `-max_iters` can limit a run to a smaller
number of timesteps as well, but reducing iterations makes testing fundamentally
more permissive, so use it only when a partial run is sufficient:

```bash
./run_scripts/run_scm.py -jax \
  -stats input/stats/multi_col_stats.in -debug 0 -max_iters 120 arm
```

### Reusable case processes

[`clubb_case_runner`](../utilities/README.md#reusable-case-runners) keeps the
standalone process alive across complete cases. Each request initializes,
advances and cleans up its model; compiled JAX programs stay in memory.
Device/environment setup continues to belong to `run_jax.py`.
Tuning retains model allocations through the existing JAX loss API.

### Timing and experimental execution traces

JAX writes `<fname_prefix>.timing` in the case's output directory at cleanup.
The basic timer runs automatically at debug levels `0` and above; negative
debug levels disable timing, matching the Fortran driver. For example:

```bash
./run_scripts/run_scm.py -jax=cpu -debug 0 -max_iters 20 bomex
```

The file uses Fortran's summary table format: call counts, inclusive and
exclusive seconds, average inclusive time and exclusive percentage. Its
`Backend: wall_clock` label identifies elapsed host time measured with Python's
`perf_counter()`, including compilation during timed calls. Basic timing adds
no device synchronization: a region can end before its submitted work completes,
so these are host-region timings rather than completed device execution times.
Fortran's default `cpu_time` summaries measure a different clock. Use untraced
process-group wall time for completed-run speed comparisons: in Dash, select
**Process-group wall time** from the Timer menu. The default
`advance_clubb_to_end` selection uses each implementation's own regional clock.
The existing [`time_clubb.py`](../utilities/time_clubb.py) utility and Dash's
Profile tab read both formats:

```bash
./utilities/time_clubb.py bomex -jax=cpu -max_iters 20 \
    -process_counts 1 -batch_sizes 1 -warmups 1 -repeats 2 \
    -output_root output/timing -name jax_bomex
```

The timing utility keeps each workload's workers alive for full-case warmups
and repetitions. Measured repetitions include fresh initialization and cleanup,
while reusing in-memory compiled programs. Workers wait for all peers before
starting each run. Process-group wall time measures completed model execution;
worker startup is recorded separately. Persistent caching also helps with the
initial warmup.

Add `profile` to `-jax` to record initialization, advancement and cleanup
with JAX's full profiler in addition to the basic summary:

```bash
./run_scripts/run_scm.py -jax=profile -debug 0 -max_iters 20 bomex
./run_scripts/run_scm.py -jax=gpu,profile -debug 0 -max_iters 20 bomex
```

`-jax=profile` uses the same default backend as bare `-jax`; an explicit CPU
selection is `-jax=cpu,profile`. Deep tracing is off by default and follows the
same negative-debug disabling rule. Tracing starts before the namelist is read;
if that namelist selects negative debug, the initial capture is stopped and
marked interrupted. Each profiled initialization saves a separate
capture beneath `output/jax_traces/`. At the end of the case, the driver prints
its directory and commands for opening that exact capture with either viewer.
It never starts a viewer itself.

Use [`view_profile.py`](view_profile.py) to open a saved capture. With no
arguments it opens the latest completed run in XProf. Select Perfetto or supply
a capture directory to revisit a particular run:

```bash
./clubb_jax/view_profile.py
./clubb_jax/view_profile.py -viewer perfetto
./clubb_jax/view_profile.py -viewer xprof output/jax_traces/trace_7h9mr__h
```

XProf is prepared automatically in an isolated tool environment, independently
of the model's JAX environment. Its server stays in the foreground until Ctrl+C.
Open **Trace Viewer** and expand the CPU process or use **Expand All** to reveal
its thread timelines. The Perfetto command opens a local page; click **Open in
Perfetto** to load the saved trace. Its temporary server exits once the browser
has fetched the file. Both commands print the browser URL for remote/headless
use; remote viewing requires forwarding the local server port. `-port` selects
another port; defaults are 8791 for XProf and 9001 for Perfetto. XProf also uses the next port for its worker service
(the preceding port when the UI uses 65535).

The same capture contains native XPlane data for XProf and
`perfetto_trace.json.gz` for [Perfetto](https://ui.perfetto.dev/), so either viewer
can be used without rerunning the model. The Perfetto file can also be opened
manually in its UI. `trace_info.json` records capture settings and separates
capture duration from export time. When deep tracing is enabled, the `.timing`
summary includes its recording overhead; leave deep tracing off for baseline
timing comparisons.

Each model timestep has a step marker for XProf's step analysis. These regions
preserve asynchronous execution: their host durations are not synchronized
device timings. Some XProf analysis tools have limited data for JAX CPU runs;
an empty accelerator chart does not mean the trace is empty. See the
[JAX profiling guide](https://docs.jax.dev/en/latest/profiling.html).
With Python function tracing enabled, Perfetto can report overlaps between
Python calls and named timer regions. It displays these on overflow tracks;
the tested viewer reports no data loss from those overlaps.
Large captures can also truncate JAX's JSON export: a tested cold CPU run with
Python tracing reached 5 million events and omitted later steps in Perfetto,
while XProf retained all steps in its native data. Turning off Python function
tracing reduced that capture enough to retain every timestep in both formats.

Normal profiling needs no environment variables. Advanced capture overrides
remain available: `CLUBB_JAX_TRACE_DIR` changes the parent directory,
`CLUBB_JAX_TRACE_HOST_LEVEL` selects host detail (0–3, default 2), and
`CLUBB_JAX_TRACE_PYTHON_LEVEL` controls Python function tracing (0 or 1, default
1). The earlier `CLUBB_JAX_TRACE=1` enable switch also remains usable.
Reducing host detail to level 1 still records many CPU operation events in the
tested JAX version. Turning off Python tracing does not guarantee a small trace.

Capture buffers and exports can be large, especially with Python tracing during
compilation. Recording can also substantially slow CPU execution; leave it off
for baseline speed measurements. Use an existing iteration limit for initial
exploration. Named regions mirror the timers in `clubb_driver.F90`, with an
additional region identifying each timestep; the corresponding JAX driver
is split between `clubb_case_initalization.py` and `advance_clubb_to_end.py`.
The timer module owns timing, capture and error cleanup. Deep trace finalization
waits for outstanding case work rather than synchronizing every region.
Instrumentation is outside compiled physics functions and does not change model
arguments.

### GPU Running

Use `-jax=gpu` to run on an NVIDIA or AMD GPU on Linux or an Apple Silicon GPU on macOS:

```bash
./run_scripts/run_scm.py -jax=gpu -stats none -debug 0 arm
```

The launcher checks GPU compatibility and prepares a separate environment on
first use. It reports an error if the requested GPU backend is unavailable.
Apple Metal runs use float32 because the plugin does not support float64.
See [Inspect Runtime Support](#inspect-runtime-support) to check your setup
without starting a run.

On NVIDIA systems, select a card using its index from `nvidia-smi`:

```bash
CUDA_VISIBLE_DEVICES=1 ./run_scripts/run_scm.py -jax=gpu arm
```

### AMD ROCm on Linux

Run JAX on CPU or an AMD GPU. GPU runs require a complete ROCm 7.2.x HIP SDK
(7.2.4 tested), Python 3.12–3.14, and access to `/dev/kfd` and `/dev/dri`.
Set `ROCM_PATH` for a nonstandard installation. The launcher installs the
matching JAX packages into `.venv-jax-rocm` on first use.

```bash
./run_scripts/run_scm.py -jax=cpu arm
CLUBB_JAX_ACCELERATOR=rocm ./run_scripts/run_scm.py -jax=gpu arm
./clubb_jax/run_jax.py -accelerator=rocm -info
```

Use `ROCR_VISIBLE_DEVICES` to select a GPU; `rocminfo` and `rocm-smi` inspect it.
GPU runs retain double precision and do not silently fall back to CPU.

Known issues:

- Radeon 8060S (`gfx1151`) needs crash workarounds: the launcher disables
  autotuning and command capture and uses native JAX Cholesky on the GPU.
  These may reduce performance; explicit settings override the defaults.
- Morrison's single-precision arithmetic can exceed JAX-versus-Fortran
  comparison tolerances on AMD GPUs. No Morrison math workarounds are applied.
- Large SILHS random permutations fail on `gfx1151` in the tested ROCm/JAX
  runtime. Use CPU for SILHS until this GPU issue is resolved.

CUDA memory preallocation is off by default unless enabled through the
`XLA_PYTHON_CLIENT_PREALLOCATE` environment variable. This lets JAX allocate
memory as needed when sharing a GPU. To enable up-front memory reservation:

```bash
./run_scripts/run_scm.py -jax=gpu,xla_prealloc arm
```

This option applies only to NVIDIA CUDA. See
[Advanced GPU Options](#advanced-gpu-options) for device selection details
and memory trade-offs.

## Mixing-length methods

The mixing length sets a distance for turbulent mixing and helps determine
turbulence damping times. [`mixing_length.py`](src/CLUBB_core/mixing_length.py)
implements all three values of the `Lscale_method` flag:

| Value | Method |
| --- | --- |
| `1` | Legacy buoyant-parcel calculation; the default |
| `2` | Diagnose length from the dissipation time scale, tau |

Select a method through the normal [SCM runner](../run_scripts/run_scm.py),
from the repository root:

```bash
./run_scripts/run_scm.py -jax bomex
./run_scripts/run_scm.py -jax -override 'Lscale_method=2' bomex
```

Method 1 uses `Lscale_mu_entrain` for parcel entrainment, formerly named `mu`.
The integer `Lscale_method` replaces `l_diag_Lscale_from_tau`: use 1 for the
former false/default setting and 2 for the former true setting.

### Measure kernel performance

[`benchmark_mixing_length.py`](tests/benchmark_mixing_length.py) measures
forward calls and parameter gradients on idealized stratified columns. It
reports compilation separately from synchronized, warmed execution. Run both
methods with the default sizes, or select method 1 and save JSON:

```bash
python3 clubb_jax/tests/benchmark_mixing_length.py
python3 clubb_jax/tests/benchmark_mixing_length.py -methods 1 \
  -columns 1 4 8 32 128 -levels 88 198 -samples 20 \
  -output_file mixing_length_timings.json
```

Rows printed to the terminal and saved in JSON identify the method, shape,
operation, compilation time, median execution time, and temporary memory.
The script prepares its managed JAX environment. It measures the full
`calc_Lscale` kernel with statistics disabled; whole-model runtime also
includes the other physics, statistics, and host work.

On a compatible NVIDIA host, select the card's physical index and CUDA profile:

```bash
CUDA_VISIBLE_DEVICES=0 CLUBB_JAX_ACCELERATOR=cuda13 \
  python3 clubb_jax/tests/benchmark_mixing_length.py -methods 1 -columns 1 4 32 128
```

## Differentiability

[`run_jax_timestep_gradient_test.py`](tests/run_jax_timestep_gradient_test.py)
checks forward and reverse derivatives of one normal driver timestep for
BOMEX and ATEX, using an artificial scalar loss and the configured default
mixing-length method. Run both cases or select one, from the repository root:

```bash
python3 clubb_jax/tests/run_jax_timestep_gradient_test.py
python3 clubb_jax/tests/run_jax_timestep_gradient_test.py -cases bomex
```

The test checks finite gradients, agreement between forward and reverse
derivatives, finite differences away from clipping boundaries, unchanged
inputs, and compiled/eager output agreement. Its file header gives the
settings and numerical limits. Passing this one-step check does not establish
long-run or all-configuration differentiability, or validate a tuning objective.

The legacy parcel method uses fixed-length scans with masks that preserve its
physical stopping rule. This supports reverse differentiation while retaining
the first energy-exhaustion point. Clipping, floors, and cutoff boundaries can still
change derivatives abruptly. Runtime parameter changes reuse compilation at
the same shapes and static settings. Standalone grid construction currently
requires ascending grids; kernel orientation checks do not remove that gate.

## Testing

The main regression test runs the same SCM cases through JAX and the original
Fortran standalone, then compares their NetCDF statistics with bindiff. This
requires a compiled Fortran standalone even though ordinary JAX runs do not.

### Basic Comparison

Compile Fortran, then compare one case:

```bash
./compile.py
./tests/run_jax_vs_fortran_cases.py -cases arm
```

The harness writes its logs, separate JAX and Fortran outputs, and final bindiff
report under `output/tests/jax_driver_test_results/`.

### Comparison Options

The harness defaults to `standard_stats.in`. Single-dash `run_scm.py` options
such as `-stats`, `-debug`, `-max_iters`, and `-dt_main` are forwarded to both
model runs, so JAX and Fortran use the same settings.
As with a normal run, reducing the statistics registry or debug level can make
comparisons faster. Do not use `-stats none` for numerical validation because
it leaves no statistics output to compare.

The case list and per-case settings are defined in
[`DEFAULT_CASES`](../tests/run_jax_vs_fortran_cases.py). Comparisons may override
native case settings; the effective overrides are logged and recorded in the
results JSON. Forwarded `-max_iters` and `-dt_main` override the configured
step limit and main timestep for both models.

Harness-specific options select cases and control parallelism or bindiff. For
example, this runs two shortened cases in parallel with faster model settings:

```bash
./tests/run_jax_vs_fortran_cases.py \
  -cases arm bomex -workers 2 \
  -stats input/stats/multi_col_stats.in -debug 0 -max_iters 120
```

`-bindiff_threshold` changes the numerical difference threshold and
`-bindiff_verbose` controls the final report detail. Each run starts with a
fresh results directory. GPU comparisons force one case worker at a time so
multiple processes do not contend for the same device.

### Focused Python Tests

Focused checks live in `pytests/auto_llm_generated_pytests/` pending human
review. The first testing stage of the CPU Jenkins job `clubb_jax` explicitly
runs this provisional coverage in the JAX environment, without Python API or
F2PY bindings. Independent analytic/NumPy references and focused runtime,
conservation, derivative and error-path checks complement the separate
Fortran executable comparisons.
From the repo root:

```bash
bash tests/run_pytests.sh -jax -include_generated --durations=15
```

Routine-level ground truth lives in
`pytests/auto_llm_generated_pytests/f2py_unit_tests/`: each check calls a JAX
routine and the same Fortran routine through the Python API (`clubb_python`)
with identical inputs, and asserts matching shapes and values. The exception is
`diagnose_upxp`, which is private in Fortran; its check evaluates the Fortran
formula on Fortran's `ddzt`. Build the Python API first; these tests fail setup
instead of skipping when the build is missing:

```bash
./compile.py -debug -python
bash tests/run_pytests.sh -jax_f2py_unit_tests -include_generated
```

Set `CLUBB_F2PY_DIR` to use a build other than `install/latest/python`.
The runner executes these tests as one pytest session, and `-jax` ignores them, and Fortran module storage
(grid, error info, PDF parameters) persists across calls in that process. Call
any Fortran routine that reads that storage through its `clubb_python` wrapper,
which sets the storage from its arguments on every call; raw `clubb_f2py` calls
are for routines that read no stored state. The wrappers also loop over
single-column routines. `pytests/f2py_unit_test_helpers.py` builds matching JAX
and Fortran grids and compares outputs with shape checks.

See [tests/README.md](../tests/README.md) for admission and Jenkins ownership.

### Manual component checks

JAX-specific case checks live in this directory's `tests/` folder. Run them
explicitly from the repository root; they prepare the JAX environment and are
separate from the focused Jenkins pytest stage.

| Check | Command |
| --- | --- |
| SILHS initialization, feedback modes and multicolumn sampling | `python3 clubb_jax/tests/run_silhs_driver_test.py` |
| BOMEX/ATEX driver gradients, JVPs and finite differences | `python3 clubb_jax/tests/run_jax_timestep_gradient_test.py` |
| Mixing-length forward and parameter-gradient kernel timings | `python3 clubb_jax/tests/benchmark_mixing_length.py` |
| Neutral/ARM initialization and timesteps with binding imports blocked | `python3 clubb_jax/tests/run_standalone_no_bindings_test.py` |
| Morrison eager/JIT species interfaces and compiled fatal diagnostics | `python3 clubb_jax/tests/run_morrison_interface_test.py` |
| Mixed-moment PDF integrals against large Monte Carlo samples | `python3 clubb_jax/tests/run_pdf_monte_carlo_test.py` |
| Rico microphysics and precipitation against Fortran statistics | `python3 clubb_jax/tests/run_rico_microphysics_oracle_test.py -stats_file output/tests/rico_oracle/rico_stats.nc` |

The Morrison and Monte Carlo checks can take several minutes. They retain the
expensive validation removed from focused pytest; they are manual checks and
are not added to the nightly Jenkins job.

The Rico check requires fresh `all_stats` output from a ten-step Fortran Rico
run; it does not generate the reference run. Missing input or a failed
numerical assertion fails the check. The gradient check's settings and limits
are described in [Differentiability](#differentiability).

## Reusable cases and runtime batches

Initialization captures the model's initial state. `set_case_initial_conditions`
resets fields/statistics and recomputes parameter-dependent derived values for
reruns or one-based runtime batches. Optional absolute iteration bounds in
`advance_clubb_to_end` support successive windows without resetting the model.
The standalone runs all configured parameter batches.

The native driver-test counterpart is `clubb_jax/src/clubb_driver_test.py`.
It accepts the standalone namelist, initializes/cleans up/reinitializes, then
advances, resets and advances again. Only the second advance writes statistics.
For example, run both frontends and compare their complete NetCDF output:

```sh
python3 run_scripts/run_scm.py bomex -jax -multicol 4 -max_iters 8 -tout 60 \
  -override 'stats_tsamp=60.' -output_dir driver_test/standalone_jax
python3 run_scripts/run_scm.py bomex -jax -driver_test -multicol 4 -max_iters 8 \
  -tout 60 -override 'stats_tsamp=60.' -output_dir driver_test/driver_jax
python3 run_scripts/run_bindiff_all.py output/driver_test/standalone_jax \
  output/driver_test/driver_jax -strict
```

The `clubb_driver` Jenkins job runs BOMEX and ATEX with distinct C8 columns and
two-column runtime batches alongside the native driver checks. Its additional
JAX-only window/batch assertions are runnable with
`python3 clubb_jax/tests/run_driver_extensions_test.py`. Restart input and
loss/tuning interfaces are separate functionality changes.

## Restarts

Use `python3 tests/run_restart_test.py bomex -jax -multicol 4` or
`python3 tests/run_restart_test.py rico_silhs -jax` to test the existing NetCDF
restart format through JAX. Initialization restores source-selected state and
absolute iteration/seeds, including retained SILHS permutations. Jenkins adds
separate normal and sequence-restart CPU stages. See
[RESTART_PORT_NOTES.md](./RESTART_PORT_NOTES.md) for reader/state limits.

## Standalone loss evaluation

`python3 run_scripts/run_scm_loss.py -jax -cases bomex -fields thlm,wp2`
runs `clubb_jax.src.clubb_standalone_loss` in the managed JAX runtime. This
small front end mirrors `src/clubb_standalone_loss.F90` and calls the loss
implementation directly; runtime setup remains in `run_jax.py`. For repeated
calls, import `clubb_jax.src.clubb_loss_driver` directly. Its five arrays have shape
(window, field, candidate). Runtime batches reuse initialized state and
compiled shapes. Shared tuner orchestration is a separate functionality PR.

Run `python3 tests/run_loss_output_consistency.py bomex -jax -fields thlm wp2`
for saved-profile and independently recomputed metric checks. Jenkins adds this
CPU check alongside the native stats/loss workflows. The explicit native-oracle
workflow requires `CLUBB_TEST_LOSS_EXE` or the installed loss executable:

```bash
.venv-jax/bin/python -m pytest -q clubb_jax/tests/test_loss_driver.py \
  clubb_jax/tests/test_silhs_loss_driver.py
```

See [TUNER_PORT.md](./TUNER_PORT.md) for loss-window and batch restrictions.

## Managed JAX tuning

Use `python3 run_scripts/run_tuner_job.py -jax` with the normal cases, fields,
parameter ranges and strategy options. Jobs retain shared baseline evaluation,
workers, durable errors, checkpoints and top-result reruns. `-jax=gpu` forwards
GPU selection to the existing launcher; GPU model execution is unvalidated.
Jenkins adds a bounded CPU JAX tuner stage to the existing tuner job.
See [tuner/README.md](../tuner/README.md) for examples and request fields.

## Requirements And Environments

No Fortran build is required for a JAX-only run. Runtime and test dependencies
are declared in [`requirements.txt`](./requirements.txt) for CPU and
[`requirements-cuda13.txt`](./requirements-cuda13.txt) for NVIDIA CUDA 13.
Apple Silicon uses the isolated [`requirements-metal.txt`](./requirements-metal.txt)
profile with Python 3.11 or 3.12 and JAX/JAXLIB 0.4.34 because Apple's experimental
plugin requires its own compatible JAX line. CPU and CUDA use JAX/JAXLIB 0.11.0
on Python 3.12 or newer; Python 3.11 is supported with JAX/JAXLIB 0.10.0.

### Automatic Setup With uv

The launcher handles the default environment automatically:

```bash
./clubb_jax/run_jax.py -init_env
./clubb_jax/run_jax.py -profile=gpu -init_env
```

It performs the following steps:

1. Reuses a supported Python already installed on the machine when possible.
2. Reuses `uv` from `PATH`, or downloads the pinned `uv` version into
   `.clubb-jax-tools/`.
3. Downloads Python 3.12 through `uv` only when no supported Python is present.
4. Creates `.venv-jax` for CPU, `.venv-jax-cuda13` for CUDA 13, or
   `.venv-jax-metal` for Apple Metal.
5. Installs and validates the matching requirements file.

Managed Python installations, the `uv` cache, and the virtualenv stay inside
the repository. The launcher does not modify the system Python or shell setup.
It hashes the requirements file and reuses a valid environment on later runs.

The managed locations and interpreter can be overridden:

```bash
PYTHON=python3.12 \
CLUBB_JAX_VENV=/path/to/clubb-jax-venv \
CLUBB_JAX_TOOLS_DIR=/path/to/clubb-jax-tools \
  ./clubb_jax/run_jax.py -init_env
```

### Inspect Runtime Support

The launcher can inspect the selected profile without creating an environment,
installing packages, or initializing JAX:

```bash
./clubb_jax/run_jax.py -profile=cpu -info
./clubb_jax/run_jax.py -profile=gpu -info
```

The report shows the detected devices, Python and JAX versions, environment
path, and whether setup is needed. Hardware discovery works before the JAX
environment is installed. The report does not initialize JAX or test available
GPU memory; the run itself reports the devices JAX actually uses.

For machine-readable output:

```bash
./clubb_jax/run_jax.py -profile=gpu -info=json
```

The CUDA 13 profile checks for an NVIDIA driver version of at least 580 and
compute capability of at least 7.5 on each exposed GPU. JAX's installed packages
supply the CUDA runtime libraries, so a local CUDA toolkit is not required.
Metal support remains experimental.

### Create A Virtualenv With uv

To create the environment yourself with `uv`:

```bash
uv venv --python 3.12 /path/to/clubb-jax-venv
uv pip install \
  --python /path/to/clubb-jax-venv/bin/python \
  -r clubb_jax/requirements.txt

CLUBB_JAX_VENV=/path/to/clubb-jax-venv \
  ./run_scripts/run_scm.py -jax arm
```

Use `requirements-cuda13.txt` on NVIDIA Linux or `requirements-metal.txt` on
Apple Silicon. In either case, `-jax=gpu` selects the native GPU profile.

### Use A Standard Virtualenv

A virtualenv created without `uv` can also be used:

```bash
python3.12 -m venv /path/to/clubb-jax-venv
/path/to/clubb-jax-venv/bin/python -m pip install \
  -r clubb_jax/requirements.txt

CLUBB_JAX_VENV=/path/to/clubb-jax-venv \
  ./run_scripts/run_scm.py -jax arm
```

The launcher never clears a custom `CLUBB_JAX_VENV`. It validates the selected
Python and installed packages, then uses `uv` to repair missing or incompatible
requirements if necessary.

Compiled JAX programs are cached under
`~/.cache/clubb-jax/compilation/<backend>`, outside the repository checkout.
Model subprocesses and later Jenkins builds can reuse matching executables.
Set `JAX_COMPILATION_CACHE_DIR` to choose another cache location; an existing
setting is preserved. JAX manages compatibility using the compiled computation,
compiler settings and device configuration. Initialization still clears the
in-memory cache when changing physics; the disk cache remains available.


### Advanced GPU Options

Select one CUDA device for a job using the identifier returned by
`./clubb_jax/run_jax.py -profile=gpu -info=json`:

```bash
./clubb_jax/run_jax.py \
  -options=gpu,device=GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx,prealloc_gpu_mem=false \
  input/case_setups/bomex_model.in
```

Replace the example identifier with a full GPU UUID from the report. The same
settings work for inspection with `-info=json`. With SCM wrappers, use
`-jax=gpu,device=UUID,prealloc_gpu_mem=true|false` as one selection value.
SCM scripts forward it unchanged; only the JAX launcher interprets modifiers.
An omitted device inherits visibility. Omitting `prealloc_gpu_mem` preserves
the launcher default/environment. Explicit single-device selection currently
supports CUDA whole GPUs. Other backends retain their existing selection rules.
The launcher and backend own validation and environment setup; UI and tuner
callers only forward the request and consume runtime information.

`CUDA_VISIBLE_DEVICES` accepts physical indices from `nvidia-smi`, full GPU
UUIDs, or unique UUID prefixes. Find UUIDs with `nvidia-smi -L`. The launcher
converts selections to full UUIDs so CUDA's enumeration order cannot change
which cards are selected. Lists retain their order: `1,0` exposes physical
GPU 1 as JAX device 0. Invalid, ambiguous, duplicate, or empty selections are
rejected before setup. Selection supports whole GPUs, not MIG instances;
all exposed GPUs must meet the CUDA requirements.

Without `xla_prealloc`, the launcher respects an existing
`XLA_PYTHON_CLIENT_PREALLOCATE` setting and otherwise defaults to `false`.
The modifier overrides that variable to `true`. Direct launcher users can pass
`-profile=gpu -xla_prealloc`.

Preallocation can reduce allocation overhead and fragmentation when a run
has the GPU to itself. Leaving it disabled lowers the initial memory footprint,
which helps on shared GPUs. JAX can still cache allocated memory; neither
setting limits total memory use or guarantees that a run will fit. The launch
report shows the effective setting.
