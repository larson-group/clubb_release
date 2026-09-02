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
number of timesteps as well, but reducing iterations makes testing fundamentally more permissive, so changing it should be considered only for rapid smoke tests:

```bash
./run_scripts/run_scm.py -jax \
  -stats input/stats/multi_col_stats.in -debug 0 -max_iters 120 arm
```

### GPU Running

To use an NVIDIA GPU with the CUDA 13 JAX packages, select the accelerator when
initializing the environment and running:

```bash
CLUBB_JAX_ACCELERATOR=cuda13 ./clubb_jax/run_jax_wrapper.sh --init_env

CLUBB_JAX_ACCELERATOR=cuda13 CUDA_VISIBLE_DEVICES=0 \
  ./run_scripts/run_scm.py -jax -stats none -debug 0 arm
```

The launcher verifies that JAX initialized the CUDA backend and fails rather
than silently falling back to CPU. `CUDA_VISIBLE_DEVICES` selects the GPU.

By default, JAX reserves 75% of the GPU's memory when the first JAX operation
runs. This reduces allocation overhead and memory fragmentation when JAX owns
the device, but it can leave too little memory for a desktop display or other
processes using a shared GPU and can cause an out-of-memory error at startup.
Disable preallocation when that reservation causes contention:

```bash
CLUBB_JAX_ACCELERATOR=cuda13 CUDA_VISIBLE_DEVICES=0 \
XLA_PYTHON_CLIENT_PREALLOCATE=false \
  ./run_scripts/run_scm.py -jax arm
```

With preallocation disabled, JAX allocates memory as the run needs it. This
usually lowers its initial footprint, but it is more vulnerable to memory
fragmentation; keep the default when the run owns the GPU and needs most of its
memory.

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
Python 3.12 or newer uses JAX/JAXLIB 0.11.0; Python 3.11 is supported with
JAX/JAXLIB 0.10.0.

### Automatic Setup With uv

The launcher handles the default environment automatically:

```bash
./clubb_jax/run_jax_wrapper.sh --init_env
```

It performs the following steps:

1. Reuses a supported Python already installed on the machine when possible.
2. Reuses `uv` from `PATH`, or downloads the pinned `uv` version into
   `.clubb-jax-tools/`.
3. Downloads Python 3.12 through `uv` only when no supported Python is present.
4. Creates `.venv-jax` for CPU or `.venv-jax-cuda13` for CUDA 13.
5. Installs and validates the matching requirements file.

Managed Python installations, the `uv` cache, and the virtualenv stay inside
the repository. The launcher does not modify the system Python or shell setup.
It hashes the requirements file and reuses a valid environment on later runs.

The managed locations and interpreter can be overridden:

```bash
PYTHON=python3.12 \
CLUBB_JAX_VENV=/path/to/clubb-jax-venv \
CLUBB_JAX_TOOLS_DIR=/path/to/clubb-jax-tools \
  ./clubb_jax/run_jax_wrapper.sh --init_env
```

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

Use `requirements-cuda13.txt` and set `CLUBB_JAX_ACCELERATOR=cuda13` for a GPU
environment.

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