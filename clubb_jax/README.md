# CLUBB JAX Driver

This directory is the staging area for a JAX-based CLUBB SCM driver.

At the moment, `clubb_jax/` is mostly a package-local clone of
[clubb_python_driver/](../clubb_python_driver/). Its current purpose is not to
provide a distinct JAX implementation yet, but to create a separate place where
the driver can be rewritten incrementally without destabilizing the existing
Python-driver path.

The intended long-term direction is different from
[clubb_python_driver/](../clubb_python_driver/):

- replace the timestep-loop internals with JAX code
- move the main timestep path away from `clubb_python_api/`
- possibly remove the Python API dependency entirely if the JAX path becomes
  self-sufficient

Because that transition has not happened yet, the current `clubb_jax/` code
should be understood as a scaffold rather than a finished driver.

The recommended conversion workflow for the eventual Fortran-to-JAX port is
documented in [JAX_CONVERSION_PLAN.md](./JAX_CONVERSION_PLAN.md).

## Current Status

Right now, the package is still functionally close to
[clubb_python_driver/](../clubb_python_driver/):

- the package structure is the same
- the SCM entry point is present at [clubb_standalone.py](./clubb_standalone.py)
- much of the execution path still depends on the compiled Python API in
  [clubb_python_api/](../clubb_python_api/)

If you need the detailed explanation of the existing driver structure, feature
coverage, and current implementation style, start with
[clubb_python_driver/README.md](../clubb_python_driver/README.md). That
documentation remains the best reference for how the cloned code currently
behaves.

## Build

From the repo root:

```bash
./compile.py [-debug] -python
```

The `-python` flag is still required today because the current `clubb_jax/`
path continues to rely on the compiled `clubb_f2py` extension provided by
[clubb_python_api/](../clubb_python_api/).

## Run

The normal SCM entry point is the existing runner script, now with `-jax`:

```bash
./run_scripts/run_scm.py -jax arm
```

This keeps the standard namelist aggregation and SCM workflow, but launches:

```bash
python -m clubb_jax.clubb_standalone
```

instead of the Fortran standalone or the Python-driver standalone.

You can also run the module directly if you already have an aggregate case
namelist:

```bash
python -m clubb_jax.clubb_standalone output/arm.in
```

## Test and Compare

The main comparison harness for this directory is:

```bash
./run_scripts/run_jax_vs_fortran_cases.py
```

That script runs a curated set of SCM cases through both:

- the `clubb_jax` standalone path
- the normal Fortran standalone path

and compares the outputs with `run_bindiff_all.py`. Results are written under
`jax_driver_test_results/`.

You can also run individual SCM cases through the regular SCM entry point:

```bash
./run_scripts/run_scm.py -jax bomex
```

## Near-Term Plan

The practical plan for this directory is:

1. keep the top-level SCM workflow usable through `-jax`
2. replace cloned driver internals with JAX implementations incrementally
3. reduce reliance on `clubb_python_api/`, especially inside the timestep loop
4. eventually make `clubb_jax/` a genuinely separate execution path rather than
   a renamed clone of the Python driver

Until that work is done, this directory should be treated as an active
development branch for a future JAX driver, not as a separate mature backend.
