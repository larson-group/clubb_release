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
- `src/Microphys/` and `src/Radiation/BUGSrad/` contain work in progress that is
  not connected to the supported standalone path.

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

### GPU Running

Use `-jax=gpu` to run on an NVIDIA GPU on Linux or an Apple Silicon GPU on macOS:

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

CUDA memory preallocation is off by default unless enabled through the
`XLA_PYTHON_CLIENT_PREALLOCATE` environment variable. This lets JAX allocate
memory as needed when sharing a GPU. To enable up-front memory reservation:

```bash
./run_scripts/run_scm.py -jax=gpu,xla_prealloc arm
```

This option applies only to NVIDIA CUDA. See
[Advanced GPU Options](#advanced-gpu-options) for device selection details
and memory trade-offs.

## Differentiability: initial test

We have a quick check that JAX can differentiate one full driver timestep for
BOMEX and ATEX. It checks forward and reverse derivatives using an artificial
scalar loss. **This is not a meaningful tuning objective or a working tuning
workflow yet**, and it does not establish that long runs or every configuration
can be differentiated.

Run it from the repository root:

```bash
.venv-jax/bin/python -m pytest -q clubb_jax/tests/test_full_timestep_grad.py
```

The [test](tests/test_full_timestep_grad.py) sets up the normal driver with:

- `l_diag_Lscale_from_tau=.true.`: currently required for this reverse-mode
  test. The default parcel-based mixing-length loops are still unsupported.
  This selects different mixing-length physics; the default remains unchanged.
- `debug=-1` and statistics disabled: keeps host error checks and file output
  out of differentiation. The test checks returned errors afterward.
- `l_stdout=False`: turns off progress printing.

Clipping and branch changes still need care when interpreting gradients. The
comments in the test explain its setup and finite-difference checks.

## Testing

The main regression test runs the same SCM cases through JAX and the original
Fortran standalone, then compares their NetCDF statistics with bindiff. This
requires a compiled Fortran standalone even though ordinary JAX runs do not.

### Basic Comparison

Compile Fortran, then compare one case:

```bash
./compile.py
./tests/run_jax_vs_fortran_cases.py --cases arm
```

The harness writes its logs, separate JAX and Fortran outputs, and final bindiff
report under `output/tests/jax_driver_test_results/`.

### Comparison Options

The harness defaults to `standard_stats.in`. Its `-stats`, `-debug`, and
`--max-iters` options are forwarded to both internal `run_scm.py` calls, so the
JAX and Fortran runs use the same model settings. As with a normal run, reducing
the statistics registry or debug level can make comparisons faster. Do not use
`-stats none` for numerical validation because it leaves no statistics output
to compare.

Harness-specific options select cases and control parallelism or bindiff. For
example, this runs two shortened cases in parallel with faster model settings:

```bash
./tests/run_jax_vs_fortran_cases.py \
  --cases arm bomex -j 2 \
  -stats input/stats/multi_col_stats.in -debug 0 --max-iters 120
```

`--bindiff-threshold` changes the numerical difference threshold,
`--bindiff-verbose` controls the final report detail, and `--keep-existing`
preserves earlier output directories. GPU comparisons force one case worker at
a time so multiple processes do not contend for the same device.

### Focused Python Tests

After a JAX run has prepared the managed environment, focused tests can be run
directly with pytest. For example:

```bash
.venv-jax/bin/python -m pytest clubb_jax/tests/test_jit_compile.py
```

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
./clubb_jax/run_jax.py --init_env
./clubb_jax/run_jax.py --profile=gpu --init_env
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
  ./clubb_jax/run_jax.py --init_env
```

### Inspect Runtime Support

The launcher can inspect the selected profile without creating an environment,
installing packages, or initializing JAX:

```bash
./clubb_jax/run_jax.py --profile=cpu --info
./clubb_jax/run_jax.py --profile=gpu --info
```

The report shows the detected devices, Python and JAX versions, environment
path, and whether setup is needed. Hardware discovery works before the JAX
environment is installed. The report does not initialize JAX or test available
GPU memory; the run itself reports the devices JAX actually uses.

For machine-readable output:

```bash
./clubb_jax/run_jax.py --profile=gpu --info=json
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

### Advanced GPU Options

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
`--profile=gpu --xla-prealloc`.

Preallocation can reduce allocation overhead and fragmentation when a run
has the GPU to itself. Leaving it disabled lowers the initial memory footprint,
which helps on shared GPUs. JAX can still cache allocated memory; neither
setting limits total memory use or guarantees that a run will fit. The launch
report shows the effective setting.
