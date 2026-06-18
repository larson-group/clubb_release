# CLUBB Dash App

The Dash app is a browser interface for common CLUBB workflows. It provides a
run tab for launching CLUBB cases and a plots tab for inspecting CLUBB NetCDF
output. The app is intended to be run from an existing CLUBB checkout.

## Install

From the repository root:

```bash
python -m pip install -r dash_app/requirements.txt
```

Compile CLUBB before using the run tab:

```bash
./compile.py
```

Plotting existing NetCDF output does not require a fresh compile.

The tune tab uses CLUBB's Python/F2PY interface, so compile with `-python`
before running tuner workflows:

```bash
./compile.py -python
```

## Run

From the repository root:

```bash
python dash_app/app.py
```

By default the app opens in a browser at port `8050`, or the next available
port. Use `python dash_app/app.py --help` for host, port, debug, and threading
options.

## Basic Workflows

- **Run tab:** choose benchmark cases and settings, launch CLUBB, and watch the
  run output in the browser.
- **Plots tab:** load one or more CLUBB output directories and make profile,
  time-height, time-series, budget, and subcolumn plots from the NetCDF files.
- **Tune tab:** configure and monitor tuner runs when the branch and local build
  support the tuner workflow. This requires a CLUBB build compiled with
  `./compile.py -python`.

## LES Benchmark Overlays

The plots tab can overlay LES benchmark data for cases that define SAM or
COAMPS benchmark files in:

```text
postprocessing/pyplotgen/config/Case_definitions.py
```

Currently those case definitions contain absolute paths for the shared CLUBB
benchmark archive, typically under:

```text
/home/pub/les_and_clubb_benchmark_runs/
```

The app only shows LES overlay options when the configured files exist on the
local machine. If the archive is mounted somewhere else, update the relevant
`sam_benchmark_file` or `coamps_benchmark_file` entries in
`Case_definitions.py`.
