# CLUBB Python Driver

This directory contains the Python SCM driver for CLUBB.

Its main purpose is to serve as a test bench for the Python API: it runs full SCM cases through Python and is intended to produce the same output as the normal Fortran standalone for supported configurations.

The driver uses the routine-level API in [`clubb_python_api/`](/home/guntherhuebler/Code/clubb_f2py/clubb_python_api), but it is a separate layer with a different job:

- `clubb_python_api/` exposes callable CLUBB routines to Python
- `clubb_python_driver/` reads case inputs, initializes state, runs the timestep loop, and drives a full case

Implementation-wise, the driver is a hybrid path:

- some parts of the Fortran standalone workflow have been ported directly to Python
- unported sections still call into the compiled Fortran through `clubb_python.clubb_api`

That makes the driver useful both as a correctness harness for the Python API and as a staging ground for further top-down Fortran-to-Python ports.

## Build

From the repo root:

```bash
./compile.py [-debug] -python
```

The `-python` flag is required because the driver depends on the compiled `clubb_f2py` extension provided by `clubb_python_api/`.

## Run

The normal entry point is still the existing SCM script:

```bash
./run_scripts/run_scm.py -python arm
```

With `-python`, [`run_scm.py`](/home/guntherhuebler/Code/clubb_f2py/run_scripts/run_scm.py) switches the execution path away from the compiled Fortran standalone binary and instead launches:

```bash
python -m clubb_python_driver.clubb_standalone
```

against the generated aggregate case namelist. In other words, `-python` keeps the normal SCM workflow but swaps the driver implementation from the Fortran executable to the Python standalone driver.

## Python vs Fortran Comparison

From the repo root:

```bash
./run_scripts/run_python_vs_fortran_cases.py
```

This runs selected SCM cases through both:

- the Python driver path
- the normal Fortran standalone path

and diffs the outputs with `run_bindiff_all.py`. Results are written under `python_driver_test_results/`.

## Unsupported Features

The Python driver is not a complete replacement for the Fortran standalone yet. The current feature gate in [`clubb_case_initalization.py`](/home/guntherhuebler/Code/clubb_f2py/clubb_python_driver/clubb_case_initalization.py) rejects these configurations:

- microphysics other than `microphys_scheme = 'none'`
- cloud water sedimentation via `l_cloud_sed = true`
- radiation schemes other than `rad_scheme = none`, `simplified`, or `simplified_bomex`
- sponge damping for any of: `thlm`, `rtm`, `uv`, `wp2`, `wp3`, `up2_vp2`
- SILHS / Latin hypercube sampling via `lh_microphys_type != disabled`
- `l_silhs_rad = true`
- `l_soil_veg = true`
- `l_restart = true`
- `l_input_fields = true`
- `l_test_grid_generalization = true`
- adaptive gridding via `grid_adapt_in_time_method > 0`

## Code Layout

```text
clubb_python_driver/
├── clubb_standalone.py
├── clubb_case_initalization.py
├── advance_clubb_to_end.py
├── advance_clubb_core.py
├── radiation.py
├── clubb_constants.py
└── io/
```

Key files:

- `clubb_standalone.py`: top-level driver entry point
- `clubb_case_initalization.py`: partially ported from `src/clubb_driver.F90`; reads namelists and case inputs, initializes model state
- `advance_clubb_to_end.py`: partially ported from `src/clubb_driver.F90`; runs the main timestep loop
- `advance_clubb_core.py`: Python port of `advance_clubb_core`
- `radiation.py`: partially ported Fortran radiation-driver logic used by the timestep loop
- `io/`: namelist, grid, sounding, and surface-file readers

Aside from the `io/` readers, the driver files here should generally be thought of as partially ported pieces of the original Fortran driver path. Only `advance_clubb_core.py` is a full port, and is not actually driver code, but has been included as a test/proof of concept - it's usage could be replaced by the api equivalent if desired.

## Relationship to `clubb_python_api`

The driver imports and uses `clubb_python.clubb_api` for lower-level routine calls. It does not duplicate the wrapper layer; it builds full-case execution on top of it and uses the API to fill gaps where the standalone workflow has not yet been ported fully into Python.

If you want the callable routine API itself, use the API package directly. If you want to run a full CLUBB SCM case through Python, use this driver layer.
